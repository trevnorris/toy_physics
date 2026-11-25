"""Signed Hessian readback on the separately rematched, held R4 pair."""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_aligned_crystal import AlignedCrystal, cell_call
from mouth_r4_aligned_stencil import EDGES, MIDCOUNT, barrier_means, samples
from mouth_r4_edge_descent import CELL, NODES, STATE, edge_data
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, chart_orientation
from mouth_r4_eulerian_release import ShapeTrial, axes_for_depth, order_grid
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_r4_newton_audit import (FULL_ERROR_TOL, HALVING_TOL, SYMMETRY_TOL, WIDTHS,
                                  PromotedGradient, directions, norm, ratio, reduce_gradient)
from mouth_r4_smooth_held import SmoothHeld
from mouth_r4_smooth_static import BITS, Gauge, LD
from mouth_whole_r4_relax import atomic_json


def matched(study, case, work):
    eps = {'primary':1e-8,'secondary':1e-10}[case]
    path = study/case/'rematched.npz'
    state = read(path.with_suffix('.json'))
    if digest(path)!=state['sha256'] or state['source_sha256']!=STATE or not state['field_matched']:
        raise ValueError('Fresh rematched state provenance failed')
    driver = SmoothHeld(study.parent,work,eps)
    with np.load(path,allow_pickle=False) as f:
        if ((int(f['radius']),int(f['winding']),int(f['outer']),int(f['stencil_points']))!=(4,1,108,39)
                or str(f['parent_sha256'])!=CHECKPOINT_SHA or str(f['source_sha256'])!=STATE
                or float(f['epsilon'])!=eps or not np.array_equal(driver.x,f['x'])):
            raise ValueError('Held preparation mismatch')
        driver.labels,driver.fields = f['labels'],f['field']
    source = study.parent/'r4_newton_held01/primary/last_active.npz'
    if digest(source)!=STATE:
        raise ValueError('Step-108 source changed')
    with np.load(source,allow_pickle=False) as f:
        if not np.array_equal(f['labels'],driver.labels):
            raise ValueError('Label equality failed at the fixed-label rematch')
    trial = ShapeTrial(driver.base,driver.prep,driver.x,8,24,6,3,2,precondition=False,
                       axes=axes_for_depth(driver.shape(driver.x).length,driver.reference_depth))
    order = order_grid(trial,driver.fields)
    return state,driver,AlignedCrystal(driver,trial,SpatialCrystal(trial.b),driver.labels,order)


def provenance():
    here=Path(__file__).resolve().parent
    files=('mouth_r4_aligned_hessian.py','mouth_r4_aligned_crystal.py','mouth_r4_aligned_stencil.cpp',
           'mouth_r4_aligned_stencil.py','mouth_r4_newton_crystal.py','mouth_r4_newton_audit.py',
           'mouth_r4_smooth_static.py','mouth_r4_smooth_fields.py')
    return {name:digest(here/name) for name in files}


class ReferenceGauge(Gauge):
    def cell_barrier(self, offsets, widths, gradient):
        points=offsets+BITS*widths[:,None,:]
        value=barrier_means(points,widths)[1]
        if not gradient:
            return value,None
        # Independent differentiation of the reviewed scalar rule, including the
        # point location. No derivative of the compiled barrier enters the oracle.
        grad=np.zeros_like(points)
        h=LD('1.e-25')
        for vertex in range(8):
            for component in range(3):
                perturbed=points.astype(np.clongdouble)
                perturbed[:,vertex,component]+=1j*h
                grad[:,vertex,component]=barrier_means(perturbed,widths)[1].imag/h
        return value,grad


def edge_directions(joint,labels):
    widths=[Fraction(float(ax[i+1]))-Fraction(float(ax[i])) for ax,i in zip(joint.opt.axes,CELL)]
    t,jmin,normal,_=edge_data(labels,widths)
    out=[]
    mobility=[]
    for vertex in (0,1):
        node=NODES[vertex]
        cm=joint.components[0]
        group=cm.groups[np.ravel_multi_index(node,joint.order.shape)]
        if cm.held[group] or np.count_nonzero(cm.groups==group)!=1:
            raise ValueError('Blocking-edge radial coordinate is fixed or tied to another node')
        packed=np.zeros_like(joint.width)
        rank=np.flatnonzero(np.flatnonzero(~cm.held)==group).item()
        # Each direction opens the edge minimum by one initial minimum per unit
        # parameter, so signed probes at both epsilon values use the same scale.
        packed[rank]=float(jmin/normal[vertex][0])
        out.append((packed,joint.lift(packed)))
        mobility.append(dict(node=node,packed_index=rank,group=int(group),melt_only=not bool(joint.active[rank]),
                             direction=packed[rank],dJ_dxi=float(normal[vertex][0])))
    return dict(cell=list(CELL),edge_start=[0,0,0],edge_end=[0,0,1],fraction=float(t),
                chart_minimum=float(jmin),moving_nodes=mobility),out


def orientation(joint, labels, spec, dirs):
    width=np.array([[np.diff(ax)[i] for ax,i in zip(joint.opt.axes,CELL)]],LD)
    indices=tuple((np.asarray(CELL)+BITS).T)
    direction=dirs[0][1]
    eps=joint.epsilon
    packed=dirs[0][0]
    rows=[]
    for change in (0.,-.75,-.95,-4.):
        shifted=labels+change*direction
        points=shifted[indices][None].astype(LD)
        fixed,t,moving=samples(points,width)
        e=EDGES.index((0,1,2))
        row=dict(parameter=change,all_27_fixed_minimum=float(fixed.min()),eight_corner_minimum=float(
                 fixed[:,MIDCOUNT==0].min()),moving_fraction=float(t[0,e]),moving_J=float(moving[0,e]))
        if change==-4:
            # Exercise both the production pass over the full label array and
            # the cell evaluator on the actual target with 27 positive fixed nodes.
            try:
                joint.evaluate(shifted)
            except InfeasibleLabels as exc:
                row['whole_cell_rejection']=exc.witness
            try:
                cell_call('values',points,width,gradient=True)
            except InfeasibleLabels as exc:
                row['aligned_rejection']=exc.witness
            row['passed']=(row['moving_J']<0<row['all_27_fixed_minimum'] and
                           row.get('whole_cell_rejection',{}).get('chart_det',0)<0 and
                           row.get('aligned_rejection',{}).get('stencil_index')==27+e)
        else:
            _,cg=cell_call('values',points,width,gradient=True)
            unit_slope=float(np.sum(cg*direction[indices]))
            certified=chart_orientation(joint.opt.axes,shifted,joint.opt.active)
            row.update(local_unit_opening_gradient=unit_slope,
                       whole_chart_lower=certified['chart_det_global_lower'])
            try:
                _,g,_=joint.evaluate(shifted)
            except InvalidSpatialCrystal as exc:
                # Keep the separate production crystal's stricter enclosure
                # refusal visible even when the chart and moving node are positive.
                row.update(admitted=False,physical_rejection=dict(reason=str(exc),**exc.diagnostics))
            else:
                global_slope=float(g@packed)
                row.update(admitted=True,total_opening_gradient=global_slope,total_opening_force=-global_slope)
            row['passed']=(row['whole_chart_lower']>0 and row['moving_J']>0 and unit_slope<0 and
                           row['admitted'] and row['total_opening_gradient']<0)
            if change!=0:
                unit_signed_work=joint.weighted.work(labels,shifted)
                # This is an independent absolute-energy change on precisely the
                # touched cells; use promoted sums before subtraction.
                touched=np.any(direction.reshape(-1,3)[joint.weighted.ids],axis=(1,2))
                ref=ReferenceGauge(joint.trial,joint.gauge.reference,joint.order)
                zero=np.zeros(6,LD)
                left=ref.evaluate(labels,zero,enabled=touched)
                right=ref.evaluate(shifted,zero,enabled=touched)
                oracle=np.sum(right-left,dtype=LD)
                rel=float(abs((unit_signed_work-oracle)/oracle))
                row.update(signed_gauge_work=float(eps*unit_signed_work),
                           direct_gauge_work=float(eps*oracle),relative_work_error=rel)
                row['passed'] &= rel<=1e-8
        rows.append(row)
    positive=rows[:-1]
    admitted=[r for r in positive if r['admitted']]
    monotone=(all(b['local_unit_opening_gradient']<a['local_unit_opening_gradient']
                  for a,b in zip(positive,positive[1:])) and len(admitted)>=2 and
              all(b['total_opening_gradient']<a['total_opening_gradient'] for a,b in zip(admitted,admitted[1:])))
    return dict(cell=spec['cell'],perturbed_radial_node=spec['moving_nodes'][0],samples=rows,
                resolved_positive_chart_count=len(positive),positive_jointly_admitted_count=len(admitted),
                all_positive_chart_probes_admitted=all(r['admitted'] for r in positive),
                untouched_restoring_force=(positive[0]['admitted'] and positive[0]['total_opening_force']>0),
                stronger_restoring_gradient=monotone,passed=all(r['passed'] for r in rows) and monotone)



def audit_site(joint, driver, promoted, oracle, spec, dirs, widths):
    labels=driver.labels
    active=joint.opt.active.ravel()
    order_ld=np.ascontiguousarray(joint.order,LD)
    baseline=labels.astype(LD)
    got=[]; wanted=[]
    row=dict(location=spec,directions={},passed=True)
    for name,(packed,lifted) in zip(('u','v'),dirs):
        if not np.array_equal(joint.lift(packed),lifted):
            raise ValueError('Lifted and packed probes disagree')
        selected=np.any(lifted.reshape(-1,3)[joint.weighted.ids],axis=(1,2))
        reference=[]
        for h in widths:
            both=[]
            for sign in (-1,1):
                arm=baseline+LD(sign*h)*lifted.astype(LD)
                physical=promoted.gradient(arm,order_ld,selected&active)
                _,unit=oracle.evaluate(arm,np.zeros(6,LD),enabled=selected,gradient=True)
                both.append(reduce_gradient(joint.components,physical+LD(driver.epsilon)*unit))
            reference.append(np.asarray((both[1]-both[0])/(2*LD(h)),float))
        product=joint.hvp(labels,packed)
        exact=reference[-1]
        discrepancy=ratio(product-exact,exact)
        sequence=[ratio(a-b,b) for a,b in zip(reference,reference[1:])]
        row['directions'][name]=dict(steps=list(widths),relative_successive_differences=sequence,
             relative_Hv_error=discrepancy,reference_l2=norm(exact),production_l2=norm(product),
             error_linf=float(np.max(abs(product-exact))),incident_cells=int(np.sum(selected)),
             incident_ordered_cells=int(np.sum(selected&active)),production_secant=joint.last_hvp)
        row['passed'] &= (discrepancy<=FULL_ERROR_TOL and sequence[-1]<=HALVING_TOL)
        got.append(product); wanted.append(exact)
    u,v=(d[0] for d in dirs)
    for name,(hu,hv) in [('production',got),('promoted_reference',wanted)]:
        uv,vu=float(u@hv),float(v@hu)
        scale=max(abs(uv),abs(vu))
        if scale==0:raise ArithmeticError('Unresolved symmetric cross response')
        rel=abs(uv-vu)/scale
        row[name+'_symmetry']=dict(uHv=uv,vHu=vu,relative=rel)
        row['passed'] &= rel<=SYMMETRY_TOL
    return row


def run(study, case, operation):
    directory=study/case
    output=directory/(operation+'.json')
    if output.exists(): raise FileExistsError(output)
    with TemporaryDirectory(dir=directory) as tmp:
        state,driver,joint=matched(study,case,Path(tmp))
        result=dict(command=sys.orig_argv,source_sha256=provenance(),state_sha256=state['sha256'],
                    epsilon=driver.epsilon,held_shape=driver.x.tolist())
        if operation=='labels':
            _,_,m=joint.evaluate(driver.labels)
            delta=m['gauge_energy']-state['final']['gauge_energy']
            result.update(material=m,label_minus_field_gauge_energy=delta,passed=(
                abs(delta)<=2e-19 and m['chart']['chart_det_global_lower']>0 and m['min_Jinv_global']>0))
        else:
            result.update(limits=dict(relative_product=FULL_ERROR_TOL,relative_reference_halving=HALVING_TOL,
                                      relative_symmetry=SYMMETRY_TOL),mantissa=int(np.finfo(LD).nmant+1),sites={})
            oracle=ReferenceGauge(joint.trial,driver.reference_depth,joint.order)
            promoted=PromotedGradient(joint.trial.b,Path(tmp))
            edge,edge_dirs=edge_directions(joint,driver.labels)
            result['orientation']=orientation(joint,driver.labels,edge,edge_dirs)
            cases=[(k,s,d,WIDTHS) for k,s,d in directions(joint.trial,joint.opt,joint.order)]
            cases.append(('limiting_edge',edge,edge_dirs,(2.**-9,2.**-10,2.**-11)))
            for sector,spec,dirs,widths in cases:
                result['sites'][sector]=audit_site(joint,driver,promoted,oracle,spec,dirs,widths)
                atomic_json(output,result)
            result['passed']=result['orientation']['passed'] and all(r['passed'] for r in result['sites'].values())
        atomic_json(output,result)
        print(json.dumps(result,allow_nan=False),flush=True)
        if not result['passed']: raise ArithmeticError('Held label check failed, see '+str(output))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('study',type=Path)
    p.add_argument('case',choices=('primary','secondary'))
    p.add_argument('operation',choices=('labels','hessian'))
    a=p.parse_args()
    run(a.study.resolve(),a.case,a.operation)
