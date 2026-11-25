"""Relax the welded R4 crystal and coherent fields with the mouth held at checkpoint 10."""
import argparse
from pathlib import Path
from time import perf_counter

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_exchange_carrier import BIAS
from mouth_r4_aligned_crystal import AlignedCrystal
from mouth_r4_aligned_order_sign import order_work
from mouth_r4_aligned_rematch import FIELD_TOL, FIELDS, STRENGTH, front
from mouth_r4_depth_cell_work import SOURCE
from mouth_r4_eulerian_relax_labels import FORCE_TOL
from mouth_r4_eulerian_release import order_grid
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_newton_held import MAX_BURST
from mouth_r4_slaved_melt import FLOW
from mouth_r4_smooth_endpoint import BASELINE, MODES, certified_state, compare
from mouth_r4_smooth_fields import GaugeGeometry, coherent, read_fields
from mouth_r4_smooth_held import SmoothHeld
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json, log, save_fields


class CollarFaces:
    """Follow the active collar nodes and the endpoints of native radial quadrature rays."""
    def __init__(self, solid, labels, order):
        self.points = solid.points
        body = solid.body
        radial = body.buffers['x']
        self.weld = radial[2*(len(radial)-1)//3]
        self.planes = np.unique(np.r_[body.buffers['t'], body.buffers['reg_t']])[1:-1]
        self.z = solid.axes[1] * body.shape.length
        self.last = self.banks(labels, order)
        self.history = [np.zeros_like(b[0], dtype='int8') for b in self.last]
        self.repeat = 0

    def banks(self, labels, order):
        def index(t, x, chi):
            return np.searchsorted(self.planes, t, side='right'), (x > self.weld) & (chi > 0)
        result = [index(labels[...,1], labels[...,0], order)]
        for t in self.points:
            for a in self.points:
                def interp(v):
                    return (1-t)*((1-a)*v[:,:-1,:-1] + a*v[:,:-1,1:]) + t*((1-a)*v[:,1:,:-1] + a*v[:,1:,1:])
                result.append(index(interp(labels[...,1]), interp(labels[...,0]), interp(order)))
        return result

    def observe(self, labels, order, peak_node=None):
        next_banks = self.banks(labels, order)
        witnesses = []
        counts = dict(nodes=0, ray_endpoints=0, reversals=0, one_twelfth=0,
                      peak_neighborhood_crossings=0, peak_neighborhood_reversals=0)
        for k, ((before, was_active), (after, is_active)) in enumerate(zip(self.last, next_banks)):
            direction = np.sign(after-before).astype('int8')
            crossed = (direction != 0) & was_active & is_active
            reverse = crossed & (self.history[k] == -direction)
            counts['nodes' if k == 0 else 'ray_endpoints'] += int(crossed.sum())
            counts['reversals'] += int(reverse.sum())
            if peak_node is not None:
                near = tuple(slice(max(0,i-1), i+2) for i in peak_node)
                counts['peak_neighborhood_crossings'] += int(crossed[near].sum())
                counts['peak_neighborhood_reversals'] += int(reverse[near].sum())
            lo, hi = np.minimum(before, after), np.maximum(before, after)
            j = int(np.argmin(abs(self.planes-1/12)))
            counts['one_twelfth'] += int(np.sum(crossed & (lo <= j) & (j < hi)))
            for node in np.argwhere(reverse)[:2]:
                i = tuple(node)
                witnesses.append(dict(kind='node' if k == 0 else 'radial_ray_endpoint',
                    index=node.tolist(), ray_rule_index=k-1 if k else None, depth=float(self.z[node[1]]
                    if k == 0 else self.z[node[1]]+self.points[(k-1)//len(self.points)]*np.diff(self.z)[node[1]]),
                    reference_planes=self.planes[lo[i]:hi[i]].tolist(), sign=int(direction[i])))
            self.history[k][~is_active] = 0
            self.history[k][crossed] = direction[crossed]
        self.last = next_banks
        self.repeat += counts['reversals']
        return dict(**counts, cumulative_reversals=self.repeat, examples=witnesses[:12])


def source_state(study, case):
    if case == 'primary':
        old = study.parent/'r4_rematch39_01'
        for name in ('primary_field_guard','primary_sweep_guard'):
            verify_guard(old/name,0)
        source = old/'primary/rematched.npz'
        gate = read(old/'primary/gate1.json')
        if not gate['passed'] or not gate['rematched_sha256'] == digest(source) == SOURCE['primary']:
            raise ValueError('Primary must start at the certified, freshly matched step 108')
        return source, 108
    entry, _, _ = certified_state(study,'primary',STRENGTH['primary'],SOURCE['primary'])
    return study/'primary/stationary_state.npz', entry['outer']


def read_joint(driver, trial, solid, *, compatible=False):
    order = order_grid(trial,driver.fields)
    joint = AlignedCrystal(driver,trial,solid,driver.labels,order)
    # Axis labels must obey the new physical trace after an order update. Raw
    # coefficients may differ only on sites with zero elastic support.
    if compatible:
        old = driver.labels
        changed = joint.compatible(old)
        if not np.array_equal(changed,old):
            if not np.array_equal(changed[joint.opt.stencil],old[joint.opt.stencil]):
                raise ArithmeticError('A crystal registry update moved an active label')
            driver.labels = changed
    value, gradient, material = joint.evaluate(driver.labels)
    physical, load, field = read_fields(driver,trial,solid,driver.labels,driver.fields)
    if abs(field['gauge_energy']-material['gauge_energy']) > 1e-15 or not np.isfinite(physical):
        raise ArithmeticError('The field and crystal evaluations of the aligned gauge disagree')
    if not all(np.isfinite(field[k]) and field[k] <= FIELD_TOL for k in FIELDS):
        raise ArithmeticError(f'Fresh fixed-label field residuals exceed {FIELD_TOL}: {field}')
    return joint, value, gradient, material, physical, load, field


def sign_sweep(driver, trial, path, material):
    sign, data = order_work(trial,driver.labels,driver.fields,driver.epsilon)
    with path.with_suffix('.npz').open('xb') as stream:
        np.savez_compressed(stream, **data)
    sign['work_sha256'] = digest(path.with_suffix('.npz'))
    sign['gauge_energy_difference'] = driver.epsilon*sign['unit_gauge']-material['gauge_energy']
    sign['passed'] = bool(sign['passed'] and abs(sign['gauge_energy_difference'])<=1e-15)
    atomic_json(path,sign)
    if not sign['passed']:
        raise ArithmeticError(f'Freshly equilibrated boundary failed its all-node work sweep at {path}')
    return sign


def endpoint(driver, trial, solid, case, step, source_sha, material, physical, field, current, sign, smallest):
    order = order_grid(trial,driver.fields)
    if 'collar_coeff' in trial.b.buffers:
        from mouth_r4_aligned_derivatives import shape_bank
        # The smooth-reference audit includes the unsplit 2^-18 depth bank.
        shape = shape_bank(driver,trial,solid,order)
        shape['six_total'] = [row['analytic'] for row in shape['rows']]
    else:
        from mouth_r4_segmented_endpoint_forces import endpoint_shape_work
        shape = endpoint_shape_work(driver,trial,solid,driver.labels,driver.fields,driver.folder)
    state = driver.folder/'stationary_state.npz'
    save_fields(state, labels=driver.labels, field=driver.fields, x=driver.x,radius=4,winding=1,
        epsilon=driver.epsilon, outer=step, parent_sha256=CHECKPOINT_SHA, initial_sha256=source_sha)
    chi = driver.fields[2,trial.first]
    if (not np.array_equal(chi[trial.number],driver.fields[2]) or not np.all((chi>=0)&(chi<=1))
            or np.any(driver.fields[1,trial.model.feed]!=0) or np.any(driver.fields[1,trial.model.bulk]!=BIAS)):
        raise ValueError('Invalid order or phase registry at the held endpoint')
    gauge = material['gauge_energy']
    # The physical stopping tolerance is in actual-energy units. Preserve the
    # unscaled melt-gauge gradient as a separate numerical readout.
    force = material['max_free_label_force']
    melted = driver.epsilon*material['max_unit_melt_residual']
    residuals = dict(label=force, melted_label=melted, unit_melt=material['max_unit_melt_residual'],
                     **{k:field[k] for k in FIELDS})
    checks = dict(native=max(force,melted,*(field[k] for k in FIELDS))<=FORCE_TOL,
              orientation=material['chart']['chart_det_global_lower']>0 and material['min_Jinv_global']>0,
              signed_work=bool(sign['passed']),six_forces=bool(shape['passed']),
              current=(current['feed_current']>0 and abs(current['port_error'])<=FORCE_TOL
                           and abs(current['disk_minus_feed_current'])<=FORCE_TOL))
    orientation = dict(material['chart'])
    if 'collar_coeff' in trial.b.buffers:
        from mouth_r4_cubic_physical_orientation import physical_certificate
        certificate = physical_certificate(trial.b, driver.labels, orientation['chart_det_global_lower'])
        orientation['physical_certificate'] = certificate
        checks['physical_orientation'] = 0 < certificate['physical_global_lower'] <= material['min_Jinv_global']
    result = dict(case=case,outer=step,initial_sha256=source_sha,state_sha256=digest(state),
                  epsilon=driver.epsilon,shape=driver.x.tolist(),phase_drop=BIAS,depth=float(trial.shape.length),
                  exit_radius=float(trial.shape.exit),physical_energy=physical,total_energy=physical+gauge,
                  crystal_energy=material['elastic_energy'],gauge_energy=gauge,gauge_to_crystal=material['gauge_to_elastic'],
                  native_residuals=residuals,orientation=orientation,
                  min_sampled_physical_Jinv=material['min_Jinv_global'],min_sampled_ordered_Jinv=material['min_Jinv_active'],
                  smallest_accepted_active_fraction=smallest,order_boundary=sign,sampled_melted_front=front(trial,driver.fields),
                  order_max_by_depth=order.max(axis=(0,2)).tolist(),depth_knots=trial.model.z[1].tolist(),current=current,
                  shape_forces=shape['six_total'],shape_audit=shape,
                  cp10_comparison={m:dict(cp10=v,now=shape['six_total'][MODES.index(m)]) for m,v in BASELINE.items()},
                  checks={k:bool(v) for k,v in checks.items()},passed=all(checks.values()))
    atomic_json(driver.folder/'endpoint.json',result)
    if not result['passed']:
        raise ArithmeticError(f'Held endpoint failed verification: {checks}')
    atomic_json(driver.folder/'HELD.json',dict(saved=result['state_sha256'],endpoint=digest(driver.folder/'endpoint.json'),
                                            epsilon=driver.epsilon,passed=True))


def run(study, case, *, start=None, configure_reference=None, reviewed_match=None, rematch_start=False):
    path, step = start if start is not None else source_state(study,case)
    directory = study/case
    directory.mkdir()
    driver = SmoothHeld(study.parent,directory,STRENGTH[case])
    with np.load(path,allow_pickle=False) as state:
        if (int(state['radius'])!=4 or int(state['winding'])!=1 or str(state['parent_sha256'])!=CHECKPOINT_SHA
                or not np.array_equal(state['x'],driver.x) or int(state['outer'])!=step):
            raise ValueError('Wrong held shape or state')
        driver.labels, driver.fields = state['labels'], state['field']
    source_sha = digest(path)
    trial = driver.trial()
    driver.current_trial = trial
    trial.smooth_order_geometry = GaugeGeometry(trial,aligned=True)
    if configure_reference is not None:
        configure_reference(trial.b)
    solid = SpatialCrystal(trial.b)
    rematch = rematch_start or case == 'secondary'
    if rematch:
        driver.fields,_,_,initial_history = coherent(driver,trial,solid,driver.labels,driver.fields,field_target=FIELD_TOL,raw_target=FIELD_TOL)
    joint, value, g, m, physical, load, field = read_joint(driver,trial,solid,compatible=rematch)
    if rematch:
        driver.save_internal(f'matched_{step:05d}',trial,driver.labels,driver.fields,
                dict(outer=step,material=m,field=field,total_energy=physical+m['gauge_energy'],
                     field_work=initial_history['order_step_energy'],source_sha256=source_sha))
    faces = CollarFaces(solid,driver.labels,joint.order)
    smallest = None
    while True:
        total = physical+m['gauge_energy']
        native = max(m['max_free_label_force'],driver.epsilon*m['max_unit_melt_residual'])
        if reviewed_match is not None:
            frozen = reviewed_match['rematched']
            if (reviewed_match['state_sha256'] != source_sha or reviewed_match['epsilon'] != driver.epsilon
                    or total != frozen['total_energy'] or m['elastic_energy'] != frozen['crystal_energy']
                    or m['max_free_label_force'] != frozen['peak_active'] or m['chart'] != frozen['chart']
                    or any(field[k] != frozen['field'][k] for k in FIELDS)):
                raise ArithmeticError('The production state differs from the independently reviewed cubic rematch')
            sign = reviewed_match['sign']
        elif case=='primary' and step==108:
            # Gate 1's exact sweep is already on disk. Source digest and the
            # fresh raw residuals were verified above; never replay the sweep.
            sign = read(study.parent/'r4_rematch39_01/primary/gate1.json')
        else:
            sign = sign_sweep(driver,trial,directory/f'order_work_{step:05d}.json',m)
        current = {k:float(v) for k,v in trial.measure(driver.fields,load).items() if k in FLOW}
        if set(current)!=set(FLOW):
            raise ArithmeticError('Missing flow diagnostics')
        if reviewed_match is not None:
            if current != frozen['current'] or not sign['passed']:
                raise ArithmeticError('Changed initial current readings or one-sided order work')
            reviewed_match = None
        if case=='primary' and step==108:
            frozen = read(study.parent/'r4_rematch39_01/primary/rematched.json')
            original = read(study.parent/'r4_rematch39_01/primary/labels.json')['material']
            if (frozen['sha256']!=source_sha or frozen['phase_drop']!=BIAS or current!=frozen['current']
                    or total!=frozen['total'] or physical!=frozen['final']['physical_energy']
                    or any(field[k]!=frozen['final'][k] for k in FIELDS)
                    or any(m[k]!=original[k] for k in ('elastic_energy','gauge_energy','max_free_label_force'))):
                raise ArithmeticError('The held worker differs from the frozen Gate-1 work or fields')
        if native<=FORCE_TOL:
            endpoint(driver,trial,solid,case,step,source_sha,m,physical,field,current,sign,smallest)
            return
        log('held_fresh_fields',case=case,step=step,physical_energy=physical,total_energy=total,
            native_label=native,field={k:field[k] for k in FIELDS},current=current,front=front(trial,driver.fields))
        burst_start = native
        for _ in range(MAX_BURST):
            start = perf_counter()
            labels, value, g, m, search = joint.step(driver.labels,value,g,m)
            work, alpha = search['energy_change'], search['alpha']
            if (not (work < 0 and alpha>1e-4 and m['chart']['chart_det_global_lower']>0 and m['min_Jinv_global']>0)
                    or not np.isfinite(work)):
                raise ArithmeticError('A proposed accepted step failed work, fraction or orientation')
            driver.labels = labels
            smallest = min(smallest or 1.,alpha)
            step += 1
            total += work
            native = max(m['max_free_label_force'],driver.epsilon*m['max_unit_melt_residual'])
            chatter = faces.observe(labels,joint.order,m['worst_free_active']['node'])
            details = dict(outer=step,labels_fixed_order=True,total_energy=total,material=m,search=search,
                           reference_crossings=chatter,source_sha256=source_sha,seconds=perf_counter()-start)
            driver.save_internal(f'step_{step:05d}',trial,labels,driver.fields,details)
            log('held_crystal_step',outer=step,force=native,worst=m['worst_free_active'],alpha=alpha,
                work=work,cg_reason=search['cg_reason'],cg_iterations=search['cg_iterations'],chart=m['chart'],
                shape_locked=True,reference_planes=chatter,seconds=details['seconds'])
            if native <= max(FORCE_TOL,burst_start/5):
                break
        driver.fields,_,_,history = coherent(driver,trial,solid,labels,driver.fields,field_target=FIELD_TOL,raw_target=FIELD_TOL)
        joint,value,g,m,physical,load,field = read_joint(driver,trial,solid,compatible=True)
        faces.observe(driver.labels,joint.order)
        explicit=physical+m['gauge_energy']
        if abs(explicit-(total+history['order_step_energy']['delta']))>64*abs(np.spacing(total)):
            raise ArithmeticError('Cumulative signed crystal and field work differ from the fresh total energy')
        driver.save_internal(f'matched_{step:05d}',trial,driver.labels,driver.fields,
            dict(outer=step,material=m,field=field,total_energy=explicit,source_sha256=source_sha))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('case',choices=['primary','secondary','pair'])
    options=parser.parse_args()
    study=options.study.resolve()
    try:
        if options.case=='pair':
            compare(study,primary_initial=SOURCE['primary'],aligned=True)
        else:
            run(study,options.case)
    except Exception as error:
        root = study/options.case
        if root.is_dir():
            details = dict(type=type(error).__name__,reason=str(error))
            for key in ('diagnostics','details','witness'):
                if hasattr(error,key):
                    details[key]=getattr(error,key)
            atomic_json(root/'ISSUE.json',details)
        raise
