"""Locate cells with exactly zero order in the receiving collar of a saved mouth."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_poly_shallow_release import CUT, Meridian, exterior, source
from mouth_supercore_front import check_guard


def sample(path, meshes):
    check_guard(path)
    with np.load(path, allow_pickle=False) as saved:
        if int(saved['radius'])!=4:
            raise ValueError('Expected entrance radius four')
        shape=Meridian(4.,saved['shape_params'],float(saved['shape_extra']))
        co=saved['registration']
    if np.any(co[1:]):
        raise ValueError('This sample uses a depth-constant registration')
    for mesh in meshes:
        result,(w,x,chi)=exterior(shape,mesh,CUT,return_side_field=True,
            side_order_work=lambda w,y:source(shape,co,w/shape.length,y,True))
        eta,t=np.meshgrid(-np.arange(1,101)/200,np.arange(1,200)/200)
        b=shape.radius(t*shape.length)[0]
        radial=b*np.expm1(-eta*np.exp(co[0]))/(CUT-b)
        collar=RegularGridInterpolator((w,x),chi)(np.stack([t*shape.length,radial],axis=-1))
        zero=np.all(np.stack([chi[1:,:-1],chi[1:,1:],chi[:-1,1:],chi[:-1,:-1]])==0,axis=0)
        wc,xc=np.meshgrid((w[1:]+w[:-1])/2,(x[1:]+x[:-1])/2,indexing='ij')
        bc=shape.radius(wc)[0]
        ec=-np.log1p(xc*(CUT-bc)/bc)/np.exp(co[0])
        inside=zero & (ec>=-.5)&(ec<0)
        print(json.dumps(dict(mesh=mesh,source=hashlib.sha256(path.read_bytes()).hexdigest(),ref_energy=result,
                  collar_min=collar.min(),collar_max=collar.max(),collar_fraction_zero=(collar==0).mean(),
                  cell_count=int(inside.sum()),cell_order_low=collar[80:120,-10:].max(),
                  w_over_depth_zero=(wc[inside]/shape.length).tolist()[::30],
                  edge_readouts=collar[:,-1][::20].tolist(),region_e=ec[inside].tolist()[::30])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state',type=Path)
    parser.add_argument('--mesh',nargs='+',type=int,default=[36,48])
    args=parser.parse_args()
    sample(args.state,args.mesh)
