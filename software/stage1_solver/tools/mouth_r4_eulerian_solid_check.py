"""Independent finite-difference variations of the two-metric spatial solid on R4."""
import argparse
import json
from pathlib import Path
from time import monotonic

import numpy as np

from mouth_cp10_nested_forces import integral, kernel
from mouth_cp10_spatial_crystal import load_context
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal


def run(scratch, output):
    label_file=scratch/'r4_cp10_labels01/labels.npz'
    prep,_,body,_,_,_,labels,chi=load_context(scratch,label_file)
    solver=SpatialCrystal(body)
    x,t,a=np.meshgrid(*solver.axes,indexing='ij')
    output.mkdir(exist_ok=False,parents=True)
    ref=kernel(output)
    row=dict(checkpoint=CHECKPOINT_SHA,labels=digest(label_file),prep=digest(prep),degree=len(solver.points),cases={})
    def save():
        path=output/'record.tmp'
        path.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
        path.replace(output/'readout.json')
    slide=x*(1-x)*np.sin(np.pi*t)
    full=np.stack((slide*np.cos(a)**2,.3*slide*np.sin(a)**2,.25*slide*np.sin(2*a)),axis=-1)
    rim=full*np.exp(-((x-.86)/.08)**2-((t-.2)/.16)**2)[...,None]
    # Use exact zeros on constrained boundary coordinates.
    for d in (full,rim):
        d[[0,-1]]=0; d[:,[0,-1]]=0;d[:,:,[0,-1],2]=0
    budgets=[]
    for name,xi in [('checkpoint',labels),('C_infinity',np.stack((x+.12*x*(1.-x*x)*np.sin(np.pi*t),t,a),axis=-1))]:
        case=row['cases'][name]={}
        start=monotonic()
        E,g,loads,orientation=solver.evaluate(xi,chi)
        case.update(energy=E,gradient_seconds=monotonic()-start,order_identity_error=float(np.sum(chi*loads)-E),
                    legacy_energy=integral(ref,body,solver.axes,np.zeros_like(solver.axes[1]),xi,chi,(solver.points,solver.weights),False)['energy'],
                    min_Jinv=orientation['min_Jinv_global'],min_active_Jinv=orientation['min_Jinv_active'],
                    min_chart_bound=orientation['min_chart_Jinv_bound'],local_minimum=orientation['worst_global'],
                    local_active_minimum=orientation['worst_active'],signed={})
        case['legacy_error']=E-case['legacy_energy']
        budgets.extend([abs(case['legacy_error'])<1e-10,abs(case['order_identity_error'])<1e-13])
        save()
        for tag,d in [('whole',full),('weld_and_front',rim)]:
            projected=float(np.sum(g*d))
            result=case['signed'][tag]=dict(projected=projected,fd=[])
            for h in [2e-6,4e-7,8e-8]:
                start=monotonic(); plus=solver.evaluate(xi+h*d,chi,False)[0]; minus=solver.evaluate(xi-h*d,chi,False)[0]
                central=(plus-minus)/(2*h)
                result['fd'].append(dict(h=h,plus=plus,minus=minus,center=central,error=central-projected,
                                          one_sided_plus=(plus-E)/h,one_sided_minus=(E-minus)/h,seconds=monotonic()-start))
                save()
            budgets.append(abs(result['fd'][-1]['error'])<2e-7)
        delta=chi*(1-chi)*(1+np.cos(2*a))/2
        liquid=(chi==0)*(t>.5)*(x>.5)
        case['order_tests']={}
        for tag,d in [('interior',delta),('melted_positive_step',liquid)]:
            h=2.e-6; shifted=solver.evaluate(xi,chi+h*d,False)[0]
            work=float(np.sum(loads*d))
            case['order_tests'][tag]=dict(work=work,energy_change=shifted-E,error=shifted-E-h*work)
            budgets.append(abs(shifted-E-h*work)<2e-13)
        save()
    trial=labels.copy()
    # Reverse the radial label axis so the first active entrance cell has a real inversion.
    trial[...,0]=1-trial[...,0]
    try:
        solver.evaluate(trial,chi,False)
    except InvalidSpatialCrystal as exc:
        row['negative_control']=exc.diagnostics
    else:
        row['negative_control']={}
    d=row['negative_control']
    budgets.append(d.get('code')==2 and d['chi']>0 and d['chart_bound']<0)
    row['passed']=all(budgets)
    save()
    if not row['passed']:
        raise ValueError('Failed pre-descent spatial-solid variation control; see signed readout')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    run(args.scratch,args.output)
