"""Two-sided reduced-work and positive-orientation gate on fixed-trace barrier states."""
import argparse
from math import isfinite
from pathlib import Path
from tempfile import mkdtemp

import numpy as np

from mouth_continuous_depth_mesh import axes_for_depth
from mouth_continuous_shape_relax import ShapeTrial
from mouth_cp10_spatial_crystal import MODES
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_barrier_shape import BarrierShapeDriver, unit_work
from mouth_r4_eulerian_relax_labels import LabelRelaxation, chart_orientation
from mouth_r4_eulerian_release import block_axes,order_grid
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_harmonic_gradient import differences
from mouth_r4_melt_barrier import MeltBarrier
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json


def evidence(scratch, case):
    original = scratch/'r4_eulerian_stage2_01'/case
    repair = scratch/'r4_slaved_melt_barrier01'/case
    archive = scratch/'r4_spatial_stall_force01'/('secondary_fixed_field' if case=='secondary' else case)
    verify_guard(original.with_name(case+'_guard'),1)
    verify_guard(repair.with_name(case+'_guard'),0)
    verify_guard(archive.with_name(archive.name+'_guard'),0)
    endpoint,fix,analytic = (read(original/'label_stall.json'),read(repair/'readout.json'),read(archive/'readout.json'))
    state = original/endpoint['file']
    if not (digest(state) == endpoint['sha256'] == fix['sha256'] == analytic['state_sha256']
            and fix['inner']['stationary'] and analytic['complete'] and fix['static']['physical_force_delta']==0.):
        raise ValueError('Inconsistent source for reduced barrier work')
    return state, repair/'fixed_active.npz', analytic


def variation(scratch, case, mode, out):
    out.mkdir(parents=True,exist_ok=False)
    state,repair,previous = evidence(scratch,case)
    with np.load(state,allow_pickle=False) as f, np.load(repair,allow_pickle=False) as g:
        if str(f['parent_sha256'])!=CHECKPOINT_SHA or int(f['radius'])!=4 or int(f['iteration'])!=10:
            raise ValueError('Wrong held state')
        driver = BarrierShapeDriver(scratch,out,float(f['epsilon']))
        driver.fields,driver.x,driver.labels = f['field'].copy(),f['x'].copy(),g['labels'].copy()
        old = f['labels'].copy()
    sh = driver.shape(driver.x)
    trial = ShapeTrial(driver.base,driver.prep,driver.x,8,24,6,3,2,precondition=False,
                       axes=axes_for_depth(sh.length,driver.reference_depth))
    solid = SpatialCrystal(trial.b)
    chi = order_grid(trial,driver.fields)
    opt = LabelRelaxation(trial.b,driver.labels,chi,driver.epsilon)
    if not np.array_equal(old[opt.stencil],driver.labels[opt.stencil]):
        raise ValueError('Repaired starting trace has changed')
    record = dict(source_sha256=digest(state),fixed_trace_sha256=digest(repair),case=case,
                  epsilon=driver.epsilon,test=mode)
    if mode=='active':
        barrier = MeltBarrier(trial.b,opt,out)
        E,g,_,_=solid.evaluate(driver.labels,chi)
        unit,gg,_,_=barrier.sample(driver.labels)
        corner=read(scratch/'r4_eulerian_corner01'/(case+'_support.json'))
        node=tuple(corner['vertex'])
        # Vary the free radial group at exactly the limiting melted/ordered vertex.
        cm=barrier.components[0]
        group=int(cm.groups[np.ravel_multi_index(node,chi.shape)])
        if not cm.free[group] or not cm.cross.getcol(np.flatnonzero(np.flatnonzero(cm.free)==group)[0]).nnz:
            raise ValueError('Expected an active label which couples to the melted extension')
        move=np.zeros_like(driver.labels)
        move[...,0]=(cm.groups==group).reshape(chi.shape)
        slope=float(np.sum((g+driver.epsilon*gg)*move,dtype=np.longdouble))
        record.update(node=list(node),component='radial',analytic_physical=float(np.sum(g*move)))
        zero=np.array([E,driver.epsilon*unit],dtype=np.longdouble)
        def bank(inc):
            start=driver.labels+inc*move
            orient=chart_orientation(opt.axes,start,opt.active)
            folder=Path(mkdtemp(prefix='active_',dir=out))
            labels,inner=barrier.solve(start,'inner',folder)
            if not inner['stationary'] or not np.array_equal(start[opt.stencil],labels[opt.stencil]):
                raise ArithmeticError('Active-direction trial did not hold the trace and resolve the gauge')
            val=solid.evaluate(labels,chi,False)[0]
            guess=barrier.sample(labels)[0]
            return np.array([val,driver.epsilon*guess],dtype=np.longdouble),dict(
                min_chart=chart_orientation(opt.axes,labels,opt.active)['chart_det_global_lower'],
                predictor_chart=orient['chart_det_global_lower'],inner_residual=inner['residual'],
                melted_label_change=float(np.max(abs(labels-start))),iterations=inner['accepted_steps'])
        record['work']=differences(bank,zero,slope,(2.e-6,1.e-6),1.)
    else:
        i=MODES.index(mode)
        direction=np.eye(6)[i]
        axes=block_axes(trial.b)
        reduced=driver.parts_at(trial,solid)
        zero=reduced(np.zeros(6))
        qzero=np.zeros(6)
        gauge_base=unit_work(sh,trial.polar_shift,driver.reference_depth,axes,driver.labels,~opt.active,qzero)
        if abs(driver.epsilon*gauge_base-zero[-1])>driver.epsilon*1.e-9:
            raise ValueError('Analytic and assembled gauge energies disagree at the anchor')
        step=1.e-24
        analytic_unit=-unit_work(sh,trial.polar_shift,driver.reference_depth,axes,driver.labels,
                                 ~opt.active,1j*step*direction).imag/step
        physical=float(sum(previous['modes'][mode]['analytic'][:-1]))
        slope=physical+driver.epsilon*analytic_unit
        record.update(analytic_physical=physical,analytic_unit_gauge=float(analytic_unit),
                      gauge_anchor=float(zero[-1]),independent_unit_error=float(gauge_base-zero[-1]/driver.epsilon))
        def bank(inc):
            energy=reduced(inc*direction)
            info=reduced.last_info
            return energy,dict(min_chart=info['min_chart_Jinv_bound'],min_chart_active=info['min_chart_Jinv_active_bound'],
                   inner_residual=info['inner']['residual'],iterations=info['inner']['accepted_steps'],
                   melted_label_change=info['melted_label_change'],
                   one_sided_unit_gauge_force=float(-(energy[-1]-zero[-1])/(inc*driver.epsilon)))
        h=1.e-5 if i==0 else 1.e-4
        record['work']=differences(bank,zero,slope,(h,h/2),-1.)
    atomic_json(out/'readout.json',record)
    print(case,mode,record['work']['analytic'],
          [(s['h'],s.get('central_minus_analytic')) for s in record['work']['steps']],flush=True)


def collect(scratch,area):
    sources=Path(__file__).resolve().parent
    inputs={p:digest(sources/p) for p in ('mouth_r4_barrier_gradient.py','mouth_r4_barrier_shape.py',
              'mouth_r4_melt_barrier.py','mouth_r4_melt_barrier.cpp','mouth_r4_harmonic_gradient.py')}
    results={}
    failed=[]
    for case in ('primary','secondary'):
        state,repair,_=evidence(scratch,case)
        sha,repair_sha=digest(state),digest(repair)
        result=results[case]=dict(source_sha256=sha,fixed_trace_sha256=repair_sha,modes={})
        for name in ['active',*MODES]:
            run=area/(case+'_'+name)
            verify_guard(run.with_name(run.name+'_guard'),0)
            item=read(run/'readout.json')
            if item['source_sha256']!=sha or item['fixed_trace_sha256']!=repair_sha:
                raise ValueError('Gauge source differs across signed modes')
            if abs(item['work']['analytic'])==0 or len(item['work']['steps'])!=2:
                failed.append((case,name,'zero/absent derivative'))
            for entry in item['work']['steps']:
                if 'central_minus_analytic' not in entry or not abs(entry['central_minus_analytic'])<=1.e-7:
                    failed.append((case,name,entry['h'],'central check'))
                for sign,bank in entry['banks'].items():
                    if not (bank['valid'] and bank['min_chart']>0 and bank['inner_residual']<=3e-8):
                        failed.append((case,name,entry['h'],sign,'inner solve or orientation'))
            result['modes'][name]=item
    gate=dict(passed=not failed,failed=failed,sources=inputs,absolute_tolerance=1.e-7,cases=results)
    atomic_json(area/'GATE.json',gate)
    if failed:
        raise ArithmeticError(f'Barrier gradient clearance failed: {failed}')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('area',type=Path)
    sub=p.add_subparsers(dest='action',required=True)
    run=sub.add_parser('variation')
    run.add_argument('case',choices=('primary','secondary'))
    run.add_argument('mode',choices=['active',*MODES])
    sub.add_parser('collect')
    args=p.parse_args()
    if args.action=='collect':
        collect(args.scratch.resolve(),args.area.resolve())
    else:
        variation(args.scratch.resolve(),args.case,args.mode,args.area.resolve()/(args.case+'_'+args.mode))
