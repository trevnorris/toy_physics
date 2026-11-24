"""Test the genuine child displacement together with the existing R4 crystal coordinates."""
import argparse
from pathlib import Path

import numpy as np

from mouth_backtrack_r4_release import FeasibleOrderTrial
from mouth_continuous_shape_relax import ContinuousRelease
from mouth_feasible_crystal_step import feasible_crystal_block, TINY_ALPHA
from mouth_r4_conforming_prisms import STAR, child_witness, refined_orientation
from mouth_r4_refined_solid import RefinedSolid
from mouth_r4_refinement_audit import STALL_SHA, inputs
from mouth_r4_refinement_review import measure
from mouth_rim_pinned_comparison import digest, read
from mouth_whole_r4_relax import atomic_json, log, save_fields


class RefinedTrial(FeasibleOrderTrial):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.b=RefinedSolid(self.b)

    def certificate(self,q,order,collect):
        return refined_orientation(self.b.parent,self.b.unpack(q),q[-3:],order,collect_activation=collect)


def paired_point(trial,q,order,location):
    b=trial.b
    offsets=b.unpack(q)
    official=child_witness(b.parent,offsets,q[-3:],order,location)
    other=measure(b.parent,offsets,q[-3:],order,location,b.seam)
    tol=64*np.finfo(float).eps
    # The independent reader uses different radial interpolants in float64
    # for the physical shape. Compare relatively even for near-zero volumes.
    for k in ('J','chi','W','chi_W'):
        v=official[k]
        if v is None and other[k] is None:
            continue
        if v is None or other[k] is None or not np.isclose(v,other[k],rtol=tol,atol=0):
            raise ArithmeticError(f'Independent physical {k} disagrees: {v}, {other[k]}')
    return dict(primary=official, independent=other)


def read_points(trial,q,u,witnesses):
    rows=[]
    for w in witnesses:
        if w is None or (*w['cell_t_x_angle'],w['triangle']) not in STAR:
            continue
        rows.append(paired_point(trial,q,u[2],w))
    return rows


def hard_gate(rep, rows):
    last=rep['accepted'][-1] if rep['accepted'] else rep['last_attempt']
    if last is None:
        return False
    tiny=last['alpha'] < TINY_ALPHA and -rep['delta_total_energy'] <= rep['negligible_drop']
    return bool(tiny and any(0 < r['primary']['chi_W'] <= 1e-16 and
                             0 < r['primary']['J_interval'][1] <= 1e-13 for r in rows
                             if r['primary']['chi_W'] is not None))


def source_gate(area):
    audit=read(area/'audit/preservation.json')
    for folder in (area/'guard_audit',area.parent/'r4_refinement_native01/guard',
                   area.parent/'r4_conformal_witness_review01/guard'):
        g=read(folder/'outcome.json')
        if not(g['limitsVerified'] and g['exitCode']==g['childOutcome']['exitCode']==0
               and g['childOutcome']['guardReason'] is None):
            raise ValueError(f'Failed audited source in {folder}')
    if not audit['conservation_pass'] or audit['parent_sha256']!=STALL_SHA:
        raise ValueError('Pre-descent controls failed')
    return audit


def probe(area):
    audit=source_gate(area)
    owner=inputs(area.parent)
    owner.trial_type=RefinedTrial
    trial=ContinuousRelease.build(owner,owner.x,False)
    q=np.r_[owner.q,np.zeros(3)]
    u=owner.u
    original=owner.stall['material_blocks'][-1]['last_overlap']['witness']
    before=read_points(trial,q,u,[original])[0]
    history=[]
    initial_grad=trial.b.evaluate(q,u[2],load=False)[1]
    if not np.array_equal(-initial_grad[-3:],audit['new_signed_generalized_force']):
        raise ArithmeticError('New displacement forces disagree with pre-descent audit')
    direction=initial_grad/trial.b.metric
    full=-direction*min(1.,.025/np.max(abs(direction)))
    if not np.any(full[-3:]):
        raise ArithmeticError('Refined trial displacement is inactive')
    out=area/'first_step'
    out.mkdir(exist_ok=False)
    atomic_json(out/'direction.json',dict(before=before,full_child_increment=full[-3:].tolist(),
                                         maximum_joint_increment=float(max(abs(full)))))
    q1,rep=feasible_crystal_block(trial,q,u,history,steps=1,
        on_event=lambda row:log('refined_material_search',**row))
    cert=trial.certificate(q1,u[2],False)
    if cert['status']!='certified':
        raise ArithmeticError('Accepted refined state has uncertified positive-order support')
    locations=[original]
    if rep['last_overlap']:
        locations.append(rep['last_overlap']['witness'])
    if rep['first_overlap']:
        locations.append(rep['first_overlap']['witness'])
    rows=read_points(trial,q1,u,locations)
    blocked=hard_gate(rep,rows)
    solid=trial.b.unpack(q1)
    p=out/'material.npz'
    save_fields(p,x=owner.x,q=q1,surplus=q1[-3:],solid=solid,field=u,radius=4,winding=1,
                state_sha=digest(owner.args.state), energy=rep['energy_after'],
                history_step=np.asarray([r[0] for r in history]),
                history_gradient=np.asarray([r[1] for r in history]),
                history_inverse=np.asarray([r[2] for r in history]))
    load=trial.crystal_load_for_match(q1,u)+trial.ext
    flow=trial.measure(u,load)
    record=dict(source_sha256=STALL_SHA,state_sha256=digest(p),accepted_child_increment=q1[-3:].tolist(),
                state_kind='fixed spatial fields after first crystal search; preceding signed checkpoint 11',
                signed_material_search=rep,gate='refined_child_barrier' if blocked else 'continue',
                geometry=dict(exit=trial.shape.exit,depth=trial.shape.length),orientation=cert,
                physical_points_on_feasible_state=rows,flow=flow,
                fully_melted_depth=trial.shape.length*flow['shaft_order']['terminal_from_fraction']['disk']['0.0'])
    atomic_json(out/'result.json',record)
    log('refined_first_step',gate=record['gate'],alpha=rep['min_alpha'],delta=rep['delta_total_energy'],
        child_increment=q1[-3:].tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    probe(p.parse_args().area.resolve())
