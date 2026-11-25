"""Held R4 alternation with a frozen-field order residual and a ten-cycle decision."""
import argparse
import math
from pathlib import Path
from time import perf_counter

import numpy as np

from mouth_exact_checkpoint_control import digest, read
from mouth_r4_aligned_rematch import FIELDS, STRENGTH, front
from mouth_r4_cubic_held import approved_reference
from mouth_r4_cubic_physical_orientation import physical_certificate
from mouth_r4_eulerian_relax_labels import FORCE_TOL, absmax
from mouth_r4_newton_held import MAX_BURST
from mouth_r4_segmented_held import run
from mouth_r4_smooth_endpoint import certified_state, compare
from mouth_r4_smooth_fields import material_load
from mouth_whole_r4_relax import atomic_json

STEP382_SHA = '346fc15d7dc678db51b420bbfde7e07ba8cbffac68aabc9d5cdbd4dd4b149f4c'


class FrozenOrder:
    """Only the deposited crystal load and the gauge change during a crystal burst."""
    def __init__(self,t,fields):
        self.chi=fields[2,t.first]
        load=t.model.material_load
        try:
            t.model.material_load=np.zeros_like(load)
            _,grad=t.model.evaluate(fields)
            self.background=t.pull.T@grad[2]
        finally:
            t.model.material_load=load
        geometry=t.smooth_order_geometry
        c=np.clip(geometry.V@self.chi,0.,1.)
        self.gauge_slope=-6*c*(1-c)

    def residual(self,t,labels,load,epsilon):
        geometry=t.smooth_order_geometry
        gauge=geometry.V.T@(self.gauge_slope*geometry.work_density(labels))
        g=self.background+t.pull.T@load+epsilon*gauge
        feasible=np.where(self.chi==0,np.minimum(g,0),np.where(self.chi==1,np.maximum(g,0),g))
        return absmax(feasible[t.free])


def cosine(x,y):
    magnitude=np.linalg.norm(x)*np.linalg.norm(y)
    return float(x@y/magnitude) if magnitude else None


def snapshot(joint,labels,chi):
    values=np.full(labels.size,np.nan)
    ordered=np.zeros(labels.size,bool)
    for c,component in enumerate(joint.components):
        reps=component.representative[~component.held]
        ids=3*reps+c
        # Store one displacement per independent label coordinate. The axis
        # registry ties the displacement, not the absolute reference position.
        values[ids]=(labels[...,c]-joint.opt.base[...,c]).ravel()[reps]
        ordered[ids]=component.free[~component.held]
    return values,ordered,chi.copy()


def movement(rows,states):
    """Compare on common physical DOFs. A persistent DOF was ordered at every endpoint."""
    registered=np.logical_and.reduce([np.isfinite(x) for x,_,_ in states])
    persistent=registered & np.logical_and.reduce([mask for _,mask,_ in states])
    previous=None
    for row,(x0,a,c0),(x1,b,c1) in zip(rows,states,states[1:]):
        shift=(x1[registered]-x0[registered])
        active=(a|b)[registered]
        full=np.where(active,shift,0.)
        permanent=shift[persistent[registered]]
        dc=c1-c0
        row['motion']=dict(persistent_ordered_dofs=int(persistent.sum()),common_free_dofs=int(registered.sum()),
            ordered_l2=float(np.linalg.norm(full)),persistent_l2=float(np.linalg.norm(permanent)),
            melt_l2=float(np.linalg.norm(shift[~active])),chi_l2=float(np.linalg.norm(dc)),
            ordered_cos=None,persistent_cos=None,chi_cos=None)
        if previous is not None:
            row['motion'].update(ordered_cos=cosine(full,previous[0]),persistent_cos=cosine(permanent,previous[1]),
                                 chi_cos=cosine(dc,previous[2]))
        previous=(full,permanent,dc)


class AdaptiveCycles:
    def __init__(self,case):
        self.case=case
        self.rows=[]
        self.samples=[]
        self.steps=[]
        self.last_toggle=None
        self.started=perf_counter()
        self.previous_chi=None
        self.cleared=case=='secondary'
        self.cycles=0

    def accepted(self,driver,t,solid,step,m,search):
        load,_=material_load(t,solid,driver.labels,driver.fields)
        order=self.frozen.residual(t,driver.labels,load,driver.epsilon)
        peak=m['max_free_label_force']
        limit=max(2*peak,FORCE_TOL)
        trigger=None
        if peak<=FORCE_TOL:
            trigger='crystal_target'
        elif order>limit:
            trigger='order_residual'
        elif len(self.steps)+1==MAX_BURST:
            trigger='twenty_steps'
        gap=abs(search['cg_relative_residual']-search['cg_recursive_relative_residual'])
        divergent=any(abs(r['explicit']-r['recursive'])>1e-10 for r in search['cg_residual_replacements'])
        if search['cg_reason']=='CG forcing attained':
            divergent |= gap>1e-10
        info=dict(step=step,force=peak,constrained_order=order,trigger_threshold=limit,trigger=trigger,
            delta_energy=search['energy_change'],alpha=search['alpha'],
            chart_lower=m['chart']['chart_det_global_lower'],min_sampled_physical_Jinv=m['min_Jinv_global'],
            cg_reason=search['cg_reason'],cg_iterations=search['cg_iterations'],cg_residual_divergence=bool(divergent))
        self.steps.append(info)
        return info

    def matched(self,driver,t,solid,joint,step,m,field,current,sign,total,converged):
        chi=driver.fields[2,t.first]
        self.frozen=FrozenOrder(t,driver.fields)
        observed=self.frozen.residual(t,driver.labels,t.model.material_load,driver.epsilon)
        if abs(observed-field['order_kkt'])>2e-12:
            raise ArithmeticError(f'Cached and native constrained order residuals differ: {observed}, {field["order_kkt"]}')
        if self.case=='primary' and not self.cleared:
            self.samples.append(snapshot(joint,driver.labels,chi[t.free]))
        geometry=physical_certificate(t.b,driver.labels,m['chart']['chart_det_global_lower'])
        common=dict(step=step,post_force=m['max_free_label_force'],post_node=m['worst_free_active'],
            post_field={k:field[k] for k in FIELDS},post_field_residual=max(field[k] for k in FIELDS),
            total_energy=total,crystal_energy=m['elastic_energy'],gauge_energy=m['gauge_energy'],
            physical_Jinv_certificate=geometry,current=current,front=front(t,driver.fields),free_zero=sign['free_zero'],
            positive_zero=sign['total']['positive'],max_gauge_over_positive=sign['max_abs_gauge_over_positive_physical'])
        paused=False
        if not self.steps:
            atomic_json(driver.folder/'baseline.json',dict(**common,source='fixed-label field rematch'))
            if converged and self.case=='primary':
                atomic_json(driver.folder/'benchmark_decision.json',{'cycles':0,'continue':False,
                    'reason':'Initial rematch already stationary; no measured second-half trend.'})
        else:
            self.cycles+=1
            previous=self.previous_chi
            toggled=(chi[t.free]>0).astype('int8')-(previous>0).astype('int8')
            burst=math.fsum(s['delta_energy'] for s in self.steps)
            record=read(driver.folder/f'matched_{step:05d}.json')
            rematch=record['field_work']['delta']
            row=dict(**common,cycle=self.cycles,crystal_steps=len(self.steps),start_force=self.start_force,
                pre_force=self.steps[-1]['force'],pre_order_kkt=self.steps[-1]['constrained_order'],
                rebound=common['post_force']-self.steps[-1]['force'],crystal_work=burst,rematch_work=rematch,
                cycle_work=burst+rematch,direct_total_change=total-self.start_energy,
                steps=list(self.steps),off_to_on=int(np.sum(toggled==1)),on_to_off=int(np.sum(toggled==-1)),
                previous_toggles_reversed=None,previous_toggle_reversal_fraction=None,
                previous_toggles_stayed=None,elapsed_seconds=perf_counter()-self.started,
                checkpoint_sha256=record['sha256'])
            if self.last_toggle is not None:
                both=(self.last_toggle!=0)&(toggled==-self.last_toggle)
                count=int(np.count_nonzero(self.last_toggle))
                row.update(previous_toggles_reversed=int(both.sum()),previous_toggles_stayed=count-int(both.sum()),
                           previous_toggle_reversal_fraction=float(both.sum()/count) if count else None)
            self.last_toggle=toggled
            self.rows.append(row)
            if self.case=='primary' and not self.cleared:
                movement(self.rows,self.samples)
                if self.cycles==10 or converged:
                    gate=decision(self.rows)
                    atomic_json(driver.folder/'benchmark_decision.json',gate)
                    self.cleared=gate['continue']
                    paused=not self.cleared
                    self.samples.clear()
            atomic_json(driver.folder/'cycles.json',self.rows)
        self.steps=[]
        self.previous_chi=chi[t.free].copy()
        self.start_force=common['post_force']
        self.start_energy=total
        self.started=perf_counter()
        return paused


def decision(rows):
    peaks=[r['post_force'] for r in rows]
    tail=peaks[len(peaks)//2:]
    descent=all(r['rematch_work']<=0 and r['cycle_work']<0 and r['direct_total_change']<0 and
                all(s['delta_energy']<0 for s in r['steps']) for r in rows)
    positive=all(r['physical_Jinv_certificate']['physical_global_lower']>0 and
                 all(s['chart_lower']>0 and s['min_sampled_physical_Jinv']>0 for s in r['steps']) for r in rows)
    clean=not any(s['cg_residual_divergence'] for r in rows for s in r['steps'])
    checks=dict(below_000150=peaks[-1]<.000150,second_half_strictly_down=len(tail)>=2 and
                all(a>b for a,b in zip(tail,tail[1:])),monotonic_energy=descent,
                positive_orientation=positive,no_cg_divergence_or_stop=clean)
    return {'cycles':len(rows),'post_rematch_peaks':peaks,'second_half':tail,'checks':checks,
            'continue':all(checks.values())}


def main(study,case):
    if case=='pair':
        compare(study,primary_initial=STEP382_SHA,aligned=True,configure_reference=approved_reference)
        return
    if case=='primary':
        path=study.parent/'r4_projector_repair01/primary/step_00382.npz'
        if digest(path)!=STEP382_SHA:
            raise ValueError('Unexpected step-382 source')
        start=(path,382)
    else:
        if not read(study/'primary/benchmark_decision.json')['continue']:
            raise ValueError('Secondary requires a positive benchmark decision')
        endpoint,_,_=certified_state(study,'primary',STRENGTH['primary'],STEP382_SHA)
        start=(study/'primary/stationary_state.npz',endpoint['outer'])
    try:
        run(study,case,start=start,configure_reference=approved_reference,rematch_start=True,
            cycles=AdaptiveCycles(case))
    except Exception as error:
        folder=study/case
        if folder.is_dir():
            data=dict(error=type(error).__name__,reason=str(error))
            for k in ('witness','details','diagnostics'):
                if hasattr(error,k):
                    data[k]=getattr(error,k)
            atomic_json(folder/'ISSUE.json',data)
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('case',choices=('primary','secondary','pair'))
    args=parser.parse_args()
    main(args.study.resolve(),args.case)
