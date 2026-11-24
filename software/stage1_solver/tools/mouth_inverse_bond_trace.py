"""Measure natural lengths on the two sides of one bonded, actively ordered R4 ring."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_free_shape_solid import Solid
from mouth_joined_angular_work import bump
from mouth_r4_conforming_prisms import _geometry, _Support
from mouth_r4_fold_order import axes_model
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_shape import preparation


def trace(root, path, expected, fixed_shape):
    prep=root/'supercore_joined_front01/adopt309/state.npz'
    if digest(path)!=expected:
        raise ValueError('Wrong input snapshot')
    with np.load(path,allow_pickle=False) as s, np.load(fixed_shape,allow_pickle=False) as ref_state:
        if not (np.array_equal(s['x'],ref_state['x']) and int(s['radius'])==4 and int(s['winding'])==1
                and str(s['state_sha'])==digest(prep)):
            raise ValueError('Wrong saved shape or physical preparation')
        sh=Meridian(4.,ref_state['shape_params'],float(ref_state['shape_extra']))
        eps=float(s['x'][5])
        d0=read(root/'r4_whole_mouth01/coarse/endpoint.json')['depth']
        body=Solid(prep,sh,eps,8,24,6,axes_model(sh,d0),degree=3,subdivisions=2)
        body.buffers['base'].fill(0.)
        order=_Support(body,s['field'][2])
        ref=preparation(4)[0]
        er,tr=body.buffers['reg_e'],body.buffers['reg_t']
        labels=body.buffers['labels'].reshape(len(tr),len(er),2)
        if (body.x[body.cut_node]!=body.cut
                or abs(-np.log(body.cut)/body.data.pressure-er[-1])>2e-15):
            raise ValueError('Probe is not the registered endpoint of the collar')
        rows=[]
        for j in range(9):
            # Tensor-product midpoint of the shared radial face (old and
            # refined chart have identical traces away from the angular split).
            mat,placed,_=_geometry(body,s['solid'],(j,body.cut_node,2,0))
            def midpoint(arr):
                return (arr[0,0]+arr[0,2]+arr[1,0]+arr[1,2])/4
            _,t,a=map(float,midpoint(mat))
            dest=midpoint(placed)
            x,tc,ac=map(float,dest)
            chi=float(order.at(dest))
            core=body.cut*body.data.core_ratio*float(ref.radius(t*ref.length)[0])
            collar=4.+np.interp(t,tr,labels[:,-1,0])
            z=tc*sh.length
            h=bump(np.array([z]),sh.length)[0][0]
            physical_r=x*float(sh.radius(z)[0])+eps*h*x**3*(3*np.cos(ac)**2-1)/2
            lever=physical_r*np.sin(ac)/np.sin(a)
            h_core,h_collar=lever/core,lever/collar
            if not min(core,collar,lever)>0:
                raise ArithmeticError('Degenerate hoop probe')
            rows.append(dict(material_t=t,material_theta=a,physical_r=physical_r,physical_w=z,
                       chi=chi,core_reference_radius=core,collar_reference_radius=collar,
                       core_hoop=h_core,collar_hoop=h_collar,
                       minimum_possible_max_relative_stretch_error=abs(h_core-h_collar)/(h_core+h_collar)))
        return dict(state_sha256=expected,bond=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    args=p.parse_args()
    root=args.root.resolve()
    cp=root/'r4_melted_survival01/coarse/state_00010.npz'
    area=root/'r4_local_refinement01/first_step'
    post=area/'material.npz'
    print(json.dumps(dict(cp10=trace(root,cp,CHECKPOINT_SHA,cp),
        post=trace(root,post,read(area/'result.json')['state_sha256'],
            root/'r4_feasible_descent01/coarse/stalled_00011.npz')),indent=2,allow_nan=False))
