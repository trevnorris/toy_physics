"""Five-step nested shape pilot from the independently checked R4 transported objective."""
import argparse
import gc
import hashlib
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter
from types import SimpleNamespace

import numpy as np

from mouth_cp10_field_shape import Bath
from mouth_cp10_nested_forces import integral, kernel
from mouth_exact_checkpoint_control import digest, read
from mouth_r4_adaptive_cycles import AdaptiveCycles, FrozenOrder
from mouth_r4_aligned_crystal import AlignedGauge
from mouth_r4_aligned_rematch import FIELD_TOL, FIELDS
from mouth_r4_eulerian_release import block_axes, order_grid
from mouth_r4_eulerian_relax_labels import InfeasibleLabels
from mouth_r4_eulerian_solid import InvalidSpatialCrystal
from mouth_r4_newton_held import MAX_BURST
from mouth_r4_cubic_physical_orientation import positive_profile
from mouth_r4_release_pilot import (GeometryFailure, ROOT, ReleaseModel, current_functional_only,
                                    forces, start_state, verify_sources)
from mouth_r4_segmented_held import read_joint, sign_sweep
from mouth_r4_smooth_endpoint import LD, promoted_bath
from mouth_r4_smooth_fields import coherent
from mouth_whole_r4_relax import atomic_json, save_fields

TIGHT = 3e-5
SHAPE_TOL = 4e-4


def analytic_forces(d,t,solid):
    """Production shape-bank first variation without its signed audit banks."""
    order = order_grid(t,d.fields)
    energy,grad,_,_ = solid.evaluate(d.labels,order)
    gauge = AlignedGauge(t,d.reference_depth,order)
    bath = Bath(t,d.base,d.fields,d.reference_depth)
    L = LD(t.shape.length)
    t_axis = gauge.axes[1]
    lift = np.where(t_axis*L<L-.12,(t_axis*L/(L-.12))**25,1.)-t_axis
    with TemporaryDirectory(dir=d.folder) as temp:
        native = integral(kernel(Path(temp)),t.b,gauge.axes,np.asarray(lift,float),d.labels,order,
                          (solid.points,solid.weights),True)
    if abs(native['energy']-energy)>1e-11:
        raise ArithmeticError('Unmatched base for shape derivative')
    result = np.array(native['analytic'],dtype=LD)
    result[0] -= np.sum(grad[...,1]*lift[None,:,None],dtype=LD)
    for k,tangent in enumerate(np.eye(6,dtype=np.clongdouble)*LD(1e-24)*1j):
        bath_force = -np.sum(promoted_bath(bath,t,d.reference_depth,tangent)).imag/LD(1e-24)
        gauge_force = -np.sum(gauge.evaluate(d.labels,tangent)).imag/LD(1e-24)
        result[k] += bath_force + d.epsilon*gauge_force
    return result.astype(float)


def save_snapshot(path, snap, epsilon):
    save_fields(path,labels=snap.labels,field=snap.fields,x=snap.q,t_axis=snap.t_axis,
                epsilon=epsilon,radius=4,winding=1)
    return digest(path)


def physical_stop(row):
    flow=row['current']
    if row['min_radius']<1 or row['exit_radius']<1:
        return 'radius below 1'
    if row['sampled_melted_front']==0:
        return 'melt reached entrance'
    if flow['feed_current']<5:
        return 'feed current below 5'
    if flow['dense_total_subsonic_share'] is None or flow['dense_total_subsonic_share']<=0:
        return 'zero dense/subsonic share'


def observe_trial(d,t,solid,cert,folder):
    try:
        return d.observe(t,solid,cert)
    except TypeError:
        # The native solver defines the dense current fraction only at positive flux.
        flow = t.measure(d.fields,t.model.material_load)
        if flow['dense_total_subsonic_share'] is not None or flow['feed_current']>=5:
            raise
        record=dict(shape=d.x.tolist(),depth=t.shape.length,exit=t.shape.exit,
                    feed_current=flow['feed_current'],dense_share=None,certificate=cert)
        atomic_json(folder/'no_intake.json',record)
        raise RuntimeError(f'Flow stop: feed current {flow["feed_current"]}; draining share undefined')


def signed_check(d,t,solid,expected,folder):
    verified,curv,bank = forces(d,t,solid)
    atomic_json(folder/'force_audit.json',bank)
    if max(abs(verified-expected))>1e-12:
        raise ArithmeticError(f'Shape audit and production analytic forces differ: {verified-expected}')
    return curv[:,0]


def matched_snapshot(d,t,solid,row,folder,name,cycle):
    path=folder/(name+'.npz')
    row['state_sha256']=save_snapshot(path,d.snapshot(t),d.epsilon)
    row['matched_cycles']=cycle
    atomic_json(folder/(name+'.json'),row)


def settle(d,t,solid,folder,outer,*,tight=False):
    """Each fixed-label field rematch counts toward the 25-cycle cap, including the first."""
    folder.mkdir()
    events=[]
    work_before = None
    for cycle in range(1,26):
        d.fields,_,_,history = coherent(d,t,solid,d.labels,d.fields,field_target=FIELD_TOL,raw_target=FIELD_TOL)
        joint,value,g,m,physical,load,field = read_joint(d,t,solid,compatible=True)
        total = physical+m['gauge_energy']
        if work_before is not None and abs(work_before+history['order_step_energy']['delta']-total)>64*abs(np.spacing(total)):
            raise ArithmeticError('Internal signed work disagrees with matched total')
        cert=d.orientation(t)
        row=observe_trial(d,t,solid,cert,folder)
        if row['total_energy']!=total or any(row['native_residuals'][k]>FIELD_TOL for k in FIELDS):
            raise ArithmeticError('Independent native energy/field mismatch after rematch')
        matched_snapshot(d,t,solid,row,folder,f'matched_{cycle:02d}',cycle)
        flow=row['current']
        if max(abs(flow['port_error']),abs(flow['disk_minus_feed_current']))>TIGHT:
            raise ArithmeticError(f'Fresh matched flow fails feed/disk/bulk conservation: {flow}')
        stop=physical_stop(row)
        if stop:
            raise RuntimeError(f'Physical pause at matched trial {outer}, cycle {cycle}: {stop}')
        target=TIGHT if tight else 1.e-4
        force=None
        if row['max_internal']<=target:
            force=analytic_forces(d,t,solid)
            if max(abs(force))<=.02:
                target=TIGHT
        if force is not None and row['max_internal']<=target:
            if tight or max(abs(force))<=SHAPE_TOL:
                sign=sign_sweep(d,t,folder/'final_order_sign.json',m)
                if sign['total']['positive']!=row['free_zero']:
                    raise ArithmeticError('Zero-node sweep failed strict positivity')
                row['order_sign']=sign
            return row,force
        if cycle==25:
            if row['max_internal']<=1e-4 and force is None:
                force=analytic_forces(d,t,solid)
            sign=sign_sweep(d,t,folder/'pause_order_sign.json',m)
            atomic_json(folder/'pause.json',dict(matched=row,order_sign=sign,
                        shape_force=None if force is None else force.tolist(),target=target))
            raise RuntimeError(f'25 field rematches exhausted, max internal={row["max_internal"]}, target={target}')
        burst=AdaptiveCycles('primary')
        burst.frozen=FrozenOrder(t,d.fields)
        cached=burst.frozen.residual(t,d.labels,t.model.material_load,d.epsilon)
        if abs(cached-field['order_kkt'])>2e-12:
            raise ArithmeticError('Cached frozen-order force differs from the matched constrained residual')
        work_before=total
        for inner in range(1,MAX_BURST+1):
            labels,value,g,m,search=joint.step(d.labels,value,g,m)
            d.labels=labels
            event=burst.accepted(d,t,solid,inner,m,search)
            if (event['cg_residual_divergence'] or not search['energy_change']<0
                    or not search['alpha']>1e-4 or not m['chart']['chart_det_global_lower']>0
                    or not m['min_Jinv_global']>0 or joint.last_hvp.get('method')!='forward_tangent'):
                raise ArithmeticError(f'Internal Newton or orientation stop: {event}')
            work_before+=search['energy_change']
            events.append(dict(cycle=cycle,**event))
            save_snapshot(folder/'last_crystal.npz',d.snapshot(t),d.epsilon)
            atomic_json(folder/'crystal_steps.json',events)
            if event['trigger'] is not None:
                break
        del joint
    raise AssertionError('Unreachable cycle')


def try_open(d,saved,q):
    sh=d.shape(q)
    scaled=sh.reference.c.T*sh.unit_length**np.arange(3,-1,-1)
    params=[sh.exit/4-1.5,*sh.params[2:],sh.extra]
    lo,hi=positive_profile(SimpleNamespace(buffers=dict(spline=scaled,shape=params)),False)
    if not float(lo)>3*abs(q[5])+.05 or not float(hi)*np.exp(.125)<27.5:
        raise GeometryFailure(f'Wall/collar clearance failed before assembly: {float(lo)}, {float(hi)}')
    try:
        return d.open(saved,q)
    except (InfeasibleLabels,InvalidSpatialCrystal) as exc:
        raise GeometryFailure(str(exc)) from exc
    except ValueError as exc:
        known=('Folded radial chart','Nonpositive quadrature Jacobian','Exterior crystal left',
               'The side mode requires','New axes alter','Order band or native cut indices changed')
        if str(exc).startswith(known):
            raise GeometryFailure(str(exc)) from exc
        raise


def relaxed_half_check(d,saved,force,folder):
    values=[]
    h=2.**-19
    for k,direction in enumerate(np.eye(6)):
        energies=[]
        for sign in (-1,1):
            t,solid,cert=try_open(d,saved,saved.q+sign*h*direction)
            energies.append(d.observe(t,solid,cert)['total_energy'])
            d.current_trial=None
            del t,solid
            gc.collect()
        values.append((energies[0]-energies[1])/(2*h))
    passed=max(abs(np.asarray(values)-force))<=1e-4 and max(abs(np.asarray(values)))<=SHAPE_TOL
    atomic_json(folder/'half_probe.json',dict(h=h,transported_at_stationarity=values,passed=bool(passed)))
    if not passed:
        raise ArithmeticError('Equilibrium half-probe shape-force test failed')


def pilot(folder,audit_path):
    verify_sources()
    audit=read(audit_path)
    gate=read(audit_path.parent/'release_review.json')
    if not (audit.get('passed') and gate['passed'] and gate['audit_sha256']==digest(audit_path)
            and gate['outer_sha256']==digest(Path(__file__))
            and gate['wrapper_sha256']==audit['code_sha256']==digest(Path(__file__).with_name('mouth_r4_release_pilot.py'))):
        raise ValueError('Pilot requires reviewed phase-1 evidence and immutable driver sources')
    # A review receipt alone cannot admit production before the Phase 1 commit.
    for name, sha in (('mouth_r4_release_outer.py',gate['outer_sha256']),
                      ('mouth_r4_release_pilot.py',gate['wrapper_sha256'])):
        rel=Path(__file__).with_name(name).resolve().relative_to(ROOT)
        blob=subprocess.check_output(['git','show',f'{gate["phase1_commit"]}:{rel}'],cwd=ROOT)
        if hashlib.sha256(blob).hexdigest()!=sha:
            raise ValueError(f'Phase 1 commit differs from the reviewed source: {rel}')
    model=ReleaseModel(folder,1e-8)
    current,prior=start_state(model)
    force=np.asarray(audit['shape_force'])
    curvature=np.asarray(audit['curvatures_wide_narrow'])[:,0]
    diagonal=lambda c: np.diag(1/np.maximum(1.,abs(c)))
    metric=diagonal(curvature)
    base=audit['baseline']
    rows=[dict(**base,outer=0,shape_force=force.tolist(),max_shape_force=max(abs(force)),
               delta_q_inf=0.,alpha=None,delta_energy_trans=0.,delta_energy_eff=0.,
               ratio=None,matched_cycles=0,wall_seconds=0.,state_sha256=prior['state_sha256'])]
    atomic_json(folder/'ledger.json',rows)
    for outer in range(1,6):
        start=perf_counter()
        work=folder/f'outer_{outer:02d}'
        work.mkdir()
        direction=metric@force
        if not force@direction>0:
            metric=diagonal(curvature)
            direction=metric@force
        direction*=min(1.,.035/max(abs(direction)))
        slope=-float(force@direction)
        if not slope<0:
            raise ArithmeticError('No strict outer descent direction')
        tests=[]
        for power in range(21):
            alpha=2.**-power
            delta=alpha*direction
            try:
                t,solid,cert=try_open(model,current,current.q+delta)
                probe=observe_trial(model,t,solid,cert,work)
                trans=probe['total_energy']-base['total_energy']
                valid=trans<1e-4*alpha*slope and trans<0
                test=dict(alpha=alpha,delta_energy_trans=trans,certified_Jinv=probe['certified_Jinv'],
                          slope=slope,armijo=bool(valid))
            except (GeometryFailure,InfeasibleLabels,InvalidSpatialCrystal) as exc:
                valid=False
                test=dict(alpha=alpha,geometric_rejection=str(exc))
            tests.append(test)
            atomic_json(work/'line_search.json',tests)
            if valid:
                break
            model.current_trial=None
            t=solid=None
            gc.collect()
        else:
            raise RuntimeError('Shape line search exhausted 2^-20')
        stop=physical_stop(probe)
        if stop:
            atomic_json(work/'transport_pause.json',probe)
            raise RuntimeError(f'Transported trial triggered {stop}')
        matched,new_force=settle(model,t,solid,work/'internal',outer,tight=outer==5)
        if outer==1:
            curvature=signed_check(model,t,solid,new_force,work)
        effective=matched['total_energy']-base['total_energy']
        if not effective<0 or matched['total_energy']>probe['total_energy']+64*abs(np.spacing(probe['total_energy'])):
            raise ArithmeticError(f'Relaxed energy failed majorization: transported={trans}, relaxed={effective}')
        current=model.snapshot(t)
        sha=save_snapshot(work/'accepted.npz',current,model.epsilon)
        row=dict(**matched,outer=outer,shape_force=new_force.tolist(),max_shape_force=max(abs(new_force)),
                 delta_q_inf=max(abs(delta)),alpha=alpha,delta_energy_trans=trans,delta_energy_eff=effective,
                 ratio=effective/trans,wall_seconds=perf_counter()-start)
        row['state_sha256']=sha
        rows.append(row)
        atomic_json(folder/'ledger.json',rows)
        model.current_trial=None
        del t,solid
        gc.collect()
        if row['max_shape_force']<=SHAPE_TOL:
            if row['max_internal']>TIGHT:
                raise ArithmeticError('Shape equilibrium lacks tight internal settlement')
            relaxed_half_check(model,current,new_force,work)
            model.epsilon=1e-10
            p,s,c=try_open(model,current,current.q)
            other,f2=settle(model,p,s,work/'secondary',outer,tight=True)
            signed_check(model,p,s,f2,work/'secondary')
            secondary=model.snapshot(p)
            secondary_sha=save_snapshot(work/'secondary_endpoint.npz',secondary,model.epsilon)
            dchi=max(abs((current.fields[2]-secondary.fields[2]).ravel()))
            differences=np.asarray(f2)-new_force
            paired=dict(primary=row,secondary=other,secondary_shape_force=f2.tolist(),
                   secondary_sha256=secondary_sha,delta_chi_max=dchi,shape_force_delta=differences.tolist(),
                   support_changes=int(np.count_nonzero((secondary.fields[2]==0)!=(current.fields[2]==0))))
            physical_delta=(other['total_energy']-other['gauge_energy'])-(row['total_energy']-row['gauge_energy'])
            paired['energy_delta']={key:other[key]-row[key] for key in ('total_energy','crystal_energy','gauge_energy')}
            paired['energy_delta']['physical_energy']=physical_delta
            paired['current_delta']={key:other['current'][key]-row['current'][key] for key in row['current']}
            paired['checks']=dict(secondary_force=max(abs(f2))<=SHAPE_TOL,
                gauge_force=max(abs(differences))<=TIGHT,secondary_internal=other['max_internal']<=TIGHT,
                same_front=other['sampled_melted_front']==row['sampled_melted_front'],order_and_support=dchi<=TIGHT,
                sign=other['order_sign']['passed'],energy=max(abs(paired['energy_delta'][k]) for k in
                    ('total_energy','physical_energy','crystal_energy'))<=1e-6,
                current=max(abs(v) for k,v in paired['current_delta'].items() if k!='total_axial_stress')<=TIGHT)
            paired['checks']={k:bool(v) for k,v in paired['checks'].items()}
            paired['passed']=all(paired['checks'].values())
            atomic_json(folder/'paired.json',paired)
            if not paired['passed']:
                raise ArithmeticError('Secondary equilibrium differs beyond tolerance; see paired.json')
            break
        y=force-new_force
        sy=float(delta@y)
        if sy>0:
            update=np.eye(6)-np.outer(delta,y)/sy
            metric=update@metric@update.T+np.outer(delta,delta)/sy
        else:
            metric=diagonal(curvature)
        force,base=new_force,matched
    atomic_json(folder/'finished.json',dict(accepted=len(rows)-1,ledger_sha256=digest(folder/'ledger.json'),
                reason='paired equilibrium' if (folder/'paired.json').exists() else 'five-step bound'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    parser.add_argument('audit',type=Path)
    args=parser.parse_args()
    output=args.folder.resolve()
    output.mkdir(parents=True,exist_ok=False)
    try:
        with current_functional_only():
            pilot(output,args.audit.resolve())
    except Exception as exc:
        record=dict(type=type(exc).__name__,reason=str(exc))
        for key in ('witness','details','diagnostics'):
            if hasattr(exc,key):
                record[key]=getattr(exc,key)
        atomic_json(output/'issue.json',record)
        raise
