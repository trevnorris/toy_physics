"""Read the two saved R4 projected-stop field banks with their common labels at fixed geometry."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import mpmath as mp
import numpy as np


EXPECTED = {
    'last_active': 'c96fcc61880f6f1e3f2bd0e05194d85344b470ef2b05af153b19856ffc0005c2',
    'held_stop': 'ae9dd7452d5291c6af95cdcfa23ae1cb6d0705fe29fc7d32cfc596608854f008',
}
MANIFEST = '625c07dc0ed83f93df23959c4310c336c49e3ef61cd751fdbe6df0f9dbf37c49'


def digest(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def peak(x):
    return float(np.max(np.abs(x), initial=0.))


def location(trial, axes, idx):
    x, t, a = (axis[i] for axis, i in zip(axes, idx))
    z = t*trial.shape.length
    y = (t-.06)/.56
    bump = 16*y**2*(1-y)**2 if 0 < y < 1 else 0.
    width = trial.shape.radius(z)[0]
    radius = x*width + x**3*bump*trial.polar_shift*(3*np.cos(a)**2-1)/2
    return dict(node=list(map(int, idx)), radius=float(radius), depth=float(z), theta=float(a))


def read_forces(trial, opt, material, unit_gradient, order, epsilon):
    rows = []
    for c, cm in enumerate(opt._components()):
        elastic = -cm.sum(material[..., c])
        gauge = -epsilon*cm.sum(unit_gradient[..., c])
        groups = np.flatnonzero(cm.free)
        strongest = {}
        for name, force in [('physical', elastic), ('total', elastic+gauge)]:
            group = groups[np.argmax(np.abs(force[groups]))]
            ids = np.argwhere((cm.groups == group).reshape(order.shape))
            strongest[name] = dict(group=int(group), signed=float(force[group]),
                                   elastic=float(elastic[group]), gauge=float(gauge[group]),
                                   locations=[dict(**location(trial,opt.axes,tuple(i)), chi=float(order[tuple(i)]))
                                              for i in ids])
        rows.append(dict(component=('radial','axial','angular')[c], strongest=strongest,
                         active_free=len(groups), dead_free=int(cm.dead.sum()),
                         max_dead_elastic=peak(elastic[cm.dead]),
                         max_unit_melt_residual=peak(cm.sum(unit_gradient[..., c])[cm.dead])))
    return rows


def order_variations(trial, driver, solid, labels, fields, joined_grad, cell_groups, chosen, epsilon):
    from mouth_cp10_label_gauge import work
    from mouth_r4_barrier_shape import unit_work
    from mouth_r4_eulerian_release import block_axes, melted_cells, order_grid

    order = order_grid(trial, fields)
    mask = melted_cells(order)
    axes = block_axes(trial.b)
    result = {}
    for name, node in [('previous_blocking_node',(16,10,1)), ('gauge_witness',chosen)]:
        ix = int(trial.number[trial.b.ids[node[0],node[2],node[1]]])
        full = np.flatnonzero(trial.number == ix)
        start = fields[2, full]
        if not trial.free[ix] or np.any(start != start[0]) or start[0] >= 1:
            raise ValueError('Order probe requires a free joined node with room on the positive side')
        if name == 'gauge_witness' and start[0] != 0:
            raise ValueError('Witness is already ordered')
        old = float(start[0])
        lost = np.zeros_like(mask)
        lost[mask] = np.any(cell_groups == ix, axis=1)
        args = (trial.shape, trial.polar_shift, driver.reference_depth, axes, labels, lost, np.zeros(6))
        unit = float(unit_work(*args)) if lost.any() else 0.
        harmonic = float(work(*args)) if lost.any() else 0.
        jump = -epsilon*unit
        samples = []
        for h in (1e-9, 1e-12):
            # Integrate the smooth physical directional derivative before combining it
            # with the finite (step-independent) change in its adjacent gauge cells.
            probe = fields.copy()
            probe[2,full] = old+h
            h_actual = float(probe[2,full[0]]-old)
            new_order = order_grid(trial,probe)
            new_mask = melted_cells(new_order)
            if np.any(new_mask & ~mask) or not np.array_equal(mask & ~new_mask,lost):
                raise ArithmeticError('Predicted and actual incident mask exchange differ')
            if h == 1e-9:
                tested = solid.evaluate(labels,new_order,False)[3]
            slopes = [float(joined_grad[ix])]
            for fraction in (.5,1.):
                probe[2,full] = old+fraction*h_actual
                _,g = trial.model.evaluate(probe)
                slopes.append(float(np.sum(g[2,full])))
            smooth = h_actual*(slopes[0]+4*slopes[1]+slopes[2])/6
            samples.append(dict(h=h_actual,physical_slopes=slopes,smooth_physical_work=smooth,
                                difference_from_trapezoid=smooth-h_actual*(slopes[0]+slopes[2])/2,
                                full_work=smooth+jump))
        result[name] = dict(**location(trial,axes,node), order=old, group=ix, joined_carrier_nodes=len(full),
                            cells_removed=np.argwhere(lost).tolist(), physical_slope=slopes[0],
                            harmonic_jump=-epsilon*harmonic, determinant_jump=-epsilon*(unit-harmonic),
                            total_gauge_jump=jump, perturbed_chart_lower=tested['min_chart_Jinv_bound'],
                            min_positive_order_Jinv_after_probe=tested['min_Jinv_active'],samples=samples)
    return result


def inspect(trial, driver, crystal, labels, fields, epsilon, workdir, probe_mask):
    from mouth_r4_eulerian_corner_audit import point
    from mouth_r4_eulerian_relax_labels import LabelRelaxation, chart_orientation
    from mouth_r4_eulerian_release import melted_cells, order_grid
    from mouth_r4_fold_order import all_low_after
    from mouth_r4_melt_barrier import MeltBarrier
    from mouth_r4_refinement_review import _collar_reference
    from mouth_r4_slaved_melt import FLOW

    model = trial.model
    if (np.any(fields[1,model.feed] != 0.) or np.any(fields[1,model.bulk] != .4510084282857874)
            or not np.array_equal(fields[2,trial.first][trial.number], fields[2])):
        raise ValueError('Changed prepared bath or split order node')
    order = order_grid(trial,fields)
    mask = melted_cells(order)
    opt = LabelRelaxation(trial.b,labels,order,epsilon)
    barrier = MeltBarrier(trial.b,opt,workdir)
    unit, gg, parts, dets = barrier.sample(labels)
    gauge = float(epsilon*unit)
    elastic, ge, source, physical = crystal.evaluate(labels,order)
    inner = np.bincount(trial.b.ids.ravel(),weights=source.transpose(0,2,1).ravel(),minlength=model.size)
    model.material_load = trial.ext+inner
    energy,g = model.evaluate(fields)
    joined = trial.pull.T @ g[2]
    chi = fields[2,trial.first]
    feasible = np.where(chi==0,np.minimum(joined,0.),np.where(chi==1,np.maximum(joined,0.),joined))
    ringsite = (order.shape[0]-2,order.shape[1]-2)
    ir,iz = ringsite
    # An Eulerian label at a fixed spatial ring. Direct elastic work depends
    # only on incident cells; report their support and the complete raw gradient.
    ring = dict(radius_depth_indices=[ir,iz], depth=float(opt.axes[1][iz]*trial.shape.length),
                spatial_radius=float(location(trial,opt.axes,(ir,iz,0))['radius']),
                order_range=[float(order[ir,iz].min()),float(order[ir,iz].max())],
                max_order_in_incident_cells=float(order[ir-1:ir+2,iz-1:iz+2].max()),
                physical_label_force_by_angle=(-ge[ir,iz]).tolist(),
                physical_label_force_sum=(-ge[ir,iz].sum(axis=0)).tolist())
    flows = trial.measure(fields, model.material_load)
    _,_,n,current,velocity,_ = model.evaluate(fields,True)
    _,_,_,_,_,_,s,cut,_,area,_ = model.D.T
    face = cut == -1
    dense = n > .5
    subsonic = (np.sum(velocity**2,axis=1)+(model.winding/s)**2) < n**4
    positive = area*np.maximum(current[:,2],0.)
    shares = {}
    for name, sector in [('dense_subsonic', dense & subsonic), ('dense_sonic_or_supersonic', dense & ~subsonic),
                         ('low_density_subsonic',~dense & subsonic),('low_density_sonic_or_supersonic',~dense & ~subsonic)]:
        value = float(positive[face & sector].sum())
        shares[name] = dict(positive_current=value, fraction=value/float(positive[face].sum()))
    depths = opt.axes[1]*trial.shape.length
    section = np.max(order,axis=(0,2))
    rows = read_forces(trial,opt,ge,gg,order,epsilon)
    with mp.workdps(60):
        old_cell = {name:point(trial.b,opt.axes,labels,order,(15,9,1),vertex,_collar_reference)
                    for name,vertex in [('old_blocker',(1,1,0)),('opposite',(0,0,1))]}
    record = dict(physical_energy=float(energy), elastic_energy=elastic,
                  gauge_energy=gauge,total_energy=float(np.longdouble(energy)+epsilon*unit),
                  gauge_over_elastic=gauge/elastic,unit_gauge_parts=parts,
                  order_elastic_deposition_difference=float(elastic-inner@fields[2]),
                  free_label_forces=rows, max_active_label_force=max(abs(r['strongest']['total']['signed']) for r in rows),
                  max_active_elastic_force=max(abs(r['strongest']['physical']['signed']) for r in rows),
                  unscaled_melt_residual=max(r['max_unit_melt_residual'] for r in rows),
                  order_raw=peak(joined[trial.free]), order_kkt=peak(feasible[trial.free]),
                  order_projected=peak(trial.projected(chi,joined)),
                  phase_raw=peak(g[1,model.masks[1]]),amplitude_raw=peak(g[0,model.masks[0]]),
                  sampled_front=all_low_after(depths,order.transpose(0,2,1),0.), depth_knots=depths.tolist(),
                  max_order_on_sampled_section=section.tolist(),melt_cells=int(mask.sum()),
                  whole_cell_chart=chart_orientation(opt.axes,labels,opt.active),
                  sampled_min_active=physical['worst_active'],previous_blocking_cell=old_cell,
                  previous_cell_vertices=order[15:17,9:11,1:3].tolist(), terminal_interior_spatial_ring=ring,
                  currents={k:float(flows[k]) for k in (*FLOW,'minimum_total_mach')},positive_current_sectors=shares)
    if (abs(shares['dense_subsonic']['fraction']-flows['dense_total_subsonic_share'])>3e-13
            or abs(float(energy)-flows['excess_energy'])>1e-9 or not np.array_equal(mask,~opt.active)):
        raise ArithmeticError('Physical current, energy or mask readback mismatch')
    if probe_mask:
        # Select a free interior angular/radial/depth node off the axis. Adding
        # support here leaves the symmetry constraints at the axis unchanged.
        groups = trial.number[trial.b.ids.transpose(0,2,1)].ravel()[barrier.ids]
        ijk = np.array(np.unravel_index(barrier.ids,order.shape))
        allowed = trial.free[groups] & (ijk[0]>1) & (ijk[0]<order.shape[0]-1)
        allowed &= (ijk[1]>0) & (ijk[1]<order.shape[1]-1) & (ijk[2]>0) & (ijk[2]<order.shape[2]-1)
        if not allowed.any():
            raise ValueError('No free interior order node in fully melted cells')
        cell_cost = barrier.weight*np.sum(dets-1-np.log(dets),axis=1)
        position = np.unravel_index(np.argmax(np.where(allowed,cell_cost[:,None],-1.)),allowed.shape)
        selected = tuple(map(int,ijk[(slice(None),)+position]))
        record['one_sided_order_work'] = order_variations(trial,driver,crystal,labels,fields,joined,groups,selected,epsilon)
    return record, mask


def main(area, output):
    pin = area/'pin/software/stage1_solver/tools'
    if digest(area/'sources.json') != MANIFEST or any(digest(pin/name)!=sha for name,sha in read(area/'sources.json').items()):
        raise ValueError('Projected run source pin differs')
    sys.path.insert(0,str(pin))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_order import KAPPA,WELL
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_r4_eulerian_release import Driver
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_stage2_outcome import verify_guard

    verify_guard(area/'primary_guard',1)
    data = {}
    for name,sha in EXPECTED.items():
        manifest = read(area/'primary'/(name+'.json'))
        path = area/'primary'/manifest['file']
        if digest(path) != sha or sha != manifest['sha256']:
            raise ValueError('Wrong saved field bank')
        with np.load(path) as st:
            if (int(st['radius']),int(st['winding']),int(st['iteration']),float(st['epsilon']),
                    str(st['parent_sha256'])) != (4,1,10,1e-8,CHECKPOINT_SHA):
                raise ValueError('Wrong preparation')
            data[name] = tuple(st[k] for k in ('field','labels','x'))
    before, after = data.values()
    if not np.array_equal(before[1],after[1]) or not np.array_equal(before[2],after[2]):
        raise ValueError('Saved material or shape moved during final field rematch')
    output.parent.mkdir(parents=True,exist_ok=True)
    result = dict(measurement_source_sha256=digest(Path(__file__)), source_manifest_sha256=MANIFEST,
                  inputs={n:dict(path=str((area/'primary'/(n+'.npz')).relative_to(area.parents[2])),sha256=v) for n,v in EXPECTED.items()},
                  command=[str(Path(__file__).name),str(area),str(output)],
                  parent_checkpoint_sha256=CHECKPOINT_SHA,epsilon=1e-8,radius=4,winding=1,held_coordinates=after[2].tolist())
    masks = []
    with TemporaryDirectory(dir=area.parent) as build:
        driver = Driver(area.parent,Path(build),1e-8)
        if not np.array_equal(driver.x,after[2]):
            raise ValueError('Checkpoint-10 geometry differs')
        shape = driver.shape(driver.x)
        trial = ShapeTrial(driver.base,driver.prep,driver.x,8,24,6,3,2,precondition=False,
                           axes=axes_for_depth(shape.length,driver.reference_depth))
        trial.model.extended_energy_sum = True
        trial.diag = (sum(g.power(2).T @ (KAPPA*trial.model.hq) for g in trial.G)
                      + trial.V.power(2).T @ (2*WELL*trial.model.hq))
        solid = SpatialCrystal(trial.b)
        for name,(fields,labels,x) in data.items():
            directory = Path(build)/name
            directory.mkdir()
            result[name],mask = inspect(trial,driver,solid,labels,fields,1e-8,directory,name=='held_stop')
            masks.append(mask)
    first,last = result['last_active'],result['held_stop']
    result['rematch'] = dict(max_abs_field_change=[peak(a-b) for a,b in zip(before[0],after[0])],
                             changed_entries=[int(np.count_nonzero(a!=b)) for a,b in zip(before[0],after[0])],
                             max_abs_label_change=peak(before[1]-after[1]),
                             entered_melt=int((masks[1]&~masks[0]).sum()),left_melt=int((masks[0]&~masks[1]).sum()),
                             physical_energy_change=last['physical_energy']-first['physical_energy'],
                             gauge_energy_change=last['gauge_energy']-first['gauge_energy'])
    doc = read(area/'primary/held_stop.json')
    diffs = {k:last[k]-doc['field'][k] for k in ('phase_raw','amplitude_raw','order_raw','order_kkt','order_projected')}
    diffs.update({k:v-doc['current'][k] for k,v in last['currents'].items()})
    diffs.update(physical_energy=last['physical_energy']-doc['current']['excess_energy'],
                 elastic_energy=last['elastic_energy']-doc['elastic_energy'],front=last['sampled_front']-doc['front'],
                 last_active_force=first['max_active_label_force']-read(area/'primary/last_active.json')['search']['active_force'])
    result['stored_minus_recomputed_sign'] = 'recomputed minus stored'
    result['recomputed_minus_stored'] = diffs
    summary = {k:last[k] for k in ('max_active_label_force','unscaled_melt_residual','order_kkt','sampled_front')}
    summary.update(signed_mask_work={k:v['total_gauge_jump'] for k,v in last['one_sided_order_work'].items()},
                   max_checkpoint_discrepancy=max(abs(v) for v in diffs.values()))
    print(json.dumps(summary,allow_nan=False),flush=True)
    if any(abs(v)>2e-10 for v in diffs.values()):
        raise ArithmeticError('Endpoint field readback differs from saved data')
    with output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('area','output'):
        parser.add_argument(key,type=Path)
    args=parser.parse_args()
    main(args.area.resolve(),args.output.resolve())
