"""Dense signed active-trace variation on both independently repaired R4 labels."""
import argparse
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest
from mouth_r4_barrier_gradient import evidence
from mouth_r4_eulerian_relax_labels import LabelRelaxation, chart_orientation
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_harmonic_gradient import differences
from mouth_r4_melt_barrier import MeltBarrier
from mouth_whole_r4_relax import atomic_json


WIDTHS=(2.e-6,1.e-6)
SEED=9256324


def mixed_direction(labels,barrier):
    """Address independent active registry groups in their native chart units."""
    rng=np.random.Generator(np.random.PCG64(SEED))
    move=np.zeros_like(labels)
    packed=[];stats=[]
    for c,cm in enumerate(barrier.components):
        axis=barrier.opt.axes[c];span=axis[-1]-axis[0]
        # Groups can include several physical nodes, in particular on the axis.
        # Bound both signs at every one of those nodes before constructing the
        # random probe. A free direction at an active box bound is unilateral.
        v=labels[...,c].ravel()
        slack=np.full(len(cm.free),float('inf'))
        np.minimum.at(slack,cm.groups,np.maximum(0.,np.minimum(v-axis[0],axis[-1]-v)))
        scale=np.minimum(.05*span,slack[cm.free]/(4*WIDTHS[0]))
        unit=rng.uniform(-1.,1.,int(cm.free.sum()))
        part=scale*unit
        groups=np.zeros(len(cm.free));groups[cm.free]=part
        move[...,c]=groups[cm.groups].reshape(labels.shape[:-1])
        touched=part!=0
        nnz=cm.cross.getnnz(axis=0)
        stats.append(dict(component=('radial','axial','angular')[c],free=int(cm.free.sum()),
             moved=int(touched.sum()),at_bound=int((~touched).sum()),
             cross_coupled_moved=int(np.count_nonzero(touched & (nnz>0))),l2=float(np.linalg.norm(part)),
             max_abs=float(np.max(np.abs(part)))))
        packed.append(part)
    if any(r['moved']<=r['free']/2 or not r['cross_coupled_moved'] for r in stats):
        raise ArithmeticError('Mixed direction lacks broad registry and interface coverage')
    direction=np.concatenate(packed)
    return move,direction,stats


def run(scratch,case,out):
    out.mkdir(exist_ok=False)
    state,repair,_=evidence(scratch,case)
    body,cp,axes=get_cp10(scratch)
    with np.load(state,allow_pickle=False) as f,np.load(repair,allow_pickle=False) as r:
        if (str(f['parent_sha256'])!=CHECKPOINT_SHA or int(f['radius'])!=4 or int(f['winding'])!=1 or
                int(f['iteration'])!=10 or not np.array_equal(f['x'],cp['x'])):
            raise ValueError('Wrong held mouth')
        eps=float(f['epsilon']); labels=r['labels'];original=f['labels']
        order=np.ascontiguousarray(f['field'][2,body.ids].transpose(0,2,1))
    opt=LabelRelaxation(body,labels,order,eps)
    if not np.array_equal(labels[opt.stencil],original[opt.stencil]):raise ValueError('Active repair trace changed')
    barrier=MeltBarrier(body,opt,out)
    solid=SpatialCrystal(body)
    elastic,physical_grad,_,_=solid.evaluate(labels,order)
    gauge,unit_grad,_,_=barrier.sample(labels)
    base_residual=max(abs(barrier.extract(unit_grad,'dead')))
    if base_residual>3e-8:raise ValueError('Baseline inner solve is unverified')
    move,direction,stats=mixed_direction(labels,barrier)
    np.save(out/'packed_direction.npy',direction,allow_pickle=False)
    zero=np.array([elastic,eps*gauge],dtype=np.longdouble)
    packed_gradient=barrier.extract(physical_grad+eps*unit_grad,'free')
    slope=float(np.dot(packed_gradient.astype(np.longdouble),direction))
    for c,s in enumerate(stats):
        s['analytic_physical']=float(np.sum(physical_grad[...,c]*move[...,c],dtype=np.longdouble))
        s['analytic_unit_gauge']=float(np.sum(unit_grad[...,c]*move[...,c],dtype=np.longdouble))
    direct=sum(s['analytic_physical']+eps*s['analytic_unit_gauge'] for s in stats)
    if abs(direct-slope)>1e-12:raise ArithmeticError('Group/nodal contraction differs')
    if abs(slope)<1e-6:raise ArithmeticError('Dense direction has unresolved net derivative')
    def bank(inc):
        start=labels+inc*move
        opt._pack(start)
        before=chart_orientation(axes,start,opt.active)
        key=('minus' if inc<0 else 'plus')+('_full' if abs(inc)==WIDTHS[0] else '_half')
        end,inner=barrier.solve(start,key,out)
        if not inner['stationary'] or not np.array_equal(start[opt.stencil],end[opt.stencil]):
            raise ArithmeticError('Transient trace or stationary melt check failed')
        value,_,_,physical=solid.evaluate(end,order,False)
        unit,grad,_,_=barrier.sample(end)
        residual=max(abs(barrier.extract(grad,'dead')))
        orient=chart_orientation(axes,end,opt.active)
        if residual>3e-8 or physical['min_Jinv_global']<=0 or physical['min_Jinv_active']<=0:
            raise ArithmeticError('Invalid bank orientation or melted residual')
        np.savez_compressed(out/(key+'.npz'),labels=end)
        return np.array([value,eps*unit],dtype=np.longdouble),dict(
            elastic=value,gauge_unit=float(unit),melted_label_change=float(np.max(np.abs(end-start))),
            predictor_chart=before['chart_det_global_lower'],min_chart=orient['chart_det_global_lower'],
            min_Jinv=physical['min_Jinv_global'],min_Jinv_active=physical['min_Jinv_active'],
            inner_residual=float(residual),iterations=inner['accepted_steps'],sha256=digest(out/(key+'.npz')))
    record=dict(case=case,epsilon=eps,source_sha256=digest(state),fixed_trace_sha256=digest(repair),
         source_test_sha256=digest(Path(__file__)),seed=SEED,packed_direction_sha256=digest(out/'packed_direction.npy'),
         component_work=stats,baseline_elastic=elastic,baseline_unit_gauge=float(gauge),baseline_inner=float(base_residual),
         analytic_from_nodes=direct,work=differences(bank,zero,slope,WIDTHS,1.))
    record['passed']=all(abs(s.get('central_minus_analytic',np.inf))<=1e-7
              and all(b['valid'] and b['iterations']>0 for b in s['banks'].values()) for s in record['work']['steps'])
    atomic_json(out/'readout.json',record)
    if not record['passed']:raise ArithmeticError('Dense three-component reduced gradient failed its signed controls')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('scratch',type=Path)
    p.add_argument('case',choices=('primary','secondary'));p.add_argument('out',type=Path)
    a=p.parse_args();run(a.scratch,a.case,a.out)
