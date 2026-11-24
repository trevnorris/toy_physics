"""Check the one-edge material refinement on the saved feasible R4 state, before descent."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, PORTS, digest, read
from mouth_exchange_carrier import BIAS
from mouth_joined_front_audit import Front
from mouth_last_rings_fluid import GRIDS
from mouth_whole_r4_relax import atomic_json

STALL_SHA = 'c99706009efdaf3ce94f5c8a626a737927ff2ce4705cb030c843b5862d908cef'
TARGET = (5, 18, 0, 0)
EPS = np.finfo(float).eps


def inputs(root):
    old = root/'r4_feasible_descent01'
    source = root/'r4_melted_survival01'
    exact = root/'r4_exact_support01'
    for base, manifest in ((old, 'sources.json'), (exact, 'inputs.json')):
        if any(digest(base/p) != h for p, h in read(base/manifest).items()):
            raise ValueError('Prior production sources disagree with their manifest')
    if any(digest(source/'pin'/p) != h for p, h in read(source/'sources.json').items()):
        raise ValueError('Carrier, shape or material source changed')
    outcome = read(old/'guard/outcome.json')
    if not (outcome['limitsVerified'] and outcome['exitCode'] == outcome['childOutcome']['exitCode'] == 1
            and outcome['childOutcome']['guardReason'] is None
            and read(old/'coarse/issue.json')['error'] == 'MaterialStall'):
        raise ValueError('Previous stop does not match the recorded material stall')
    st = read(old/'coarse/stall.json')
    archive = old/'coarse'/st['stall_state']
    cp = read(source/'coarse/latest.json')
    if not (digest(archive) == st['sha256'] == STALL_SHA and
            digest(source/'coarse'/cp['file']) == cp['sha256'] == CHECKPOINT_SHA and cp['iteration'] == 10):
        raise ValueError('Incorrect restart state')
    previous = Path(read(source/'coarse/invocation.json')['command'][-3])/'coarse'
    cmd = read(previous/'invocation.json')['command']
    original = Path(cmd[-5])/'coarse'
    proof = read(Path(cmd[-3]))['results']['coarse']
    state = Path(read(original/'invocation.json')['command'][5])
    if digest(state) != proof['preparation_sha256'] or digest(original/'endpoint.json') != proof['endpoint_sha256']:
        raise ValueError('Wrong prepared crystal or base geometry')
    with np.load(archive, allow_pickle=False) as d:
        if int(d['radius']) != 4 or int(d['winding']) != 1 or str(d['state_sha']) != digest(state):
            raise ValueError('Wrong physical preparation in snapshot')
        x, q, u, offsets = (d[k].copy() for k in ('x','q','field','solid'))
    if BIAS != .4510084282857874:
        raise ValueError('Reservoir bias changed')
    return SimpleNamespace(base=Front(state, GRIDS['double'], 1), args=SimpleNamespace(state=state),
        discretization=proof['solid_discretization'], reference_depth=read(original/'endpoint.json')['depth'],
        x=x, q=q, u=u, offsets=offsets, cp=cp, stall=st,
        authority=read(root/'r4_feasible_stall_readback01/readout.json'),
        proposal=read(root/'r4_exact_fulltrial01/proposal.json'))


def discrepancy(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise ArithmeticError('Incompatible or non-finite conservation readouts')
    # Each component has its own absolute scale. A large bath entry cannot
    # hide a changed force or a changed small regional elastic energy.
    error = np.abs(a-b)
    tolerance = 64*EPS*np.maximum(1., np.maximum(abs(a), abs(b)))
    return dict(max_abs=float(np.max(error)),max_fraction_tolerance=float(np.max(error/tolerance)),
                pass_roundoff=bool(np.all(error <= tolerance)))


def spatial_neighborhood():
    # All cells sharing at least a vertex with the target, including the next
    # polar layer. Six change basis; each neighboring integral is checked.
    from mouth_r4_conforming_prisms import STAR
    region=set(STAR)
    for j in range(4,7):
        for i in range(17,20):
            for k in range(2):
                for triangle in range(2):
                    region.add((j,i,k,triangle))
    return sorted(region)


def total_fields(trial, body, q, fields):
    e, g, _, _ = body.evaluate(q, fields[2], load=False)
    _, _, deposition, _ = body.lower_bound_load(q, fields[2])
    load = deposition + trial.ext
    trial.model.material_load = load
    energy, gradient = trial.model.evaluate(fields)
    flow = trial.measure(fields, load)
    return dict(total=energy, elastic_energy=e, forces=-g, deposition=deposition, field_gradient=gradient,
                currents=np.array([flow[k] for k in PORTS])),flow


def run_audit(owner, trial, out):
    from mouth_r4_conforming_prisms import geometry_audit, refined_orientation
    from mouth_r4_refined_solid import RefinedSolid
    parent, q, u = trial.b, owner.q, owner.u
    if not np.array_equal(parent.unpack(q), owner.offsets):
        raise ValueError('Displacement or reference labels differ from the saved state')
    if not np.array_equal(u, trial.transported(u)):
        raise ValueError('Saved spatial fields changed under zero transport')
    child = RefinedSolid(parent)
    extq = np.r_[q, np.zeros(3)]
    old, flow = total_fields(trial, parent, q, u)
    new, _ = total_fields(trial, child, extq, u)
    archive = owner.authority
    if not (archive['stalled_snapshot_sha256'] == STALL_SHA and old['total'] == archive['energy'] and
            old['elastic_energy'] == archive['solid_energy'] and flow['feed_current'] == archive['flow']['feed'] and
            flow['port_error'] < 3e-5 and abs(flow['disk_minus_feed_current']) < 3e-5):
        raise ValueError('Parent readout fails independent archived state or port checks')
    checks = {k:discrepancy(v, new[k][:len(q)] if k == 'forces' else new[k]) for k,v in old.items()}
    regional = []
    for key in spatial_neighborhood():
        data=[child.evaluate_region(extq, u[2], key, refined=s, load=True, gradient=True)
              for s in (False,True)]
        a,b=data
        row=dict(cell=list(key),parent_energy=a[0], child_sum_energy=b[0],
                 energy=discrepancy(a[0],b[0]),old_coordinate_force=discrepancy(a[1][:len(q)],b[1][:len(q)]),
                 order_load=discrepancy(a[2],b[2]))
        regional.append(row)
    geo = geometry_audit(parent, owner.offsets, u[2])
    cert=refined_orientation(parent,owner.offsets,extq[-3:],u[2])
    control=geo['negative_control']
    geo_ok=(cert['status']=='certified' and cert['counts']['refined_children_checked']==12
            and geo['exact_max_reference_or_placement_delta']==geo['exact_max_reference_or_placement_derivative_delta']==0
            and geo['split_angle']==child.seam and geo['split_parent_fraction']==child.seam_fraction
            and control['child_status']==control['global_status']=='active_fold'
            and control['order']==.001 and control['chart_J']<=0)
    finer=[]
    for key in (TARGET,'star'):
        a,b=[child.evaluate_region(extq,u[2],key,refined=s,load=True,angular_subdivisions=2)
             for s in (False,True)]
        finer.append(dict(region=key,old=a[0],refined=b[0],energy=discrepancy(a[0],b[0]),
                           force=discrepancy(a[1][:len(q)],b[1][:len(q)]),load=discrepancy(a[2],b[2])))
    report = dict(orientation=cert, geometry_pass=bool(geo_ok), reintegrated=finer,parent_sha256=STALL_SHA,baseline_energy=old['total'],baseline_solid_energy=old['elastic_energy'],
        quadrature_seam_angle=float(child.seam), seam_fraction=float(child.seam_fraction),
        total_checks=checks, regional=regional, old_coordinate_count=len(q),
        new_signed_generalized_force=new['forces'][-3:].tolist(),
        currents=dict(zip(PORTS,map(float,old['currents']))), geometry=geo)
    ok=all(c['pass_roundoff'] for c in checks.values()) and all(r[k]['pass_roundoff']
                for r in regional for k in ('energy','old_coordinate_force','order_load'))
    ok=bool(ok and geo_ok and all(r[k]['pass_roundoff'] for r in finer for k in ('energy','force','load')))
    report['conservation_pass']=ok
    atomic_json(out/'preservation.json',report)
    return ok

if __name__=='__main__':
    import argparse
    from mouth_continuous_shape_relax import ContinuousRelease
    from mouth_exact_melt_trial import ExactMeltTrial
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('out',type=Path)
    args=parser.parse_args()
    args.out.mkdir(exist_ok=False,parents=True)
    owner=inputs(args.root.resolve())
    owner.trial_type=ExactMeltTrial
    run_audit(owner,ContinuousRelease.build(owner,owner.x,False),args.out)
