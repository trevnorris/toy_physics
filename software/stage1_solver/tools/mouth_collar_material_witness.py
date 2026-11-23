"""Independent compactly supported work and material stretches of a saved collar field."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from numpy.polynomial.legendre import legval, legder

from mouth_collar_envelope_review import check_input, material_work
from mouth_internal_collar import CollarPair
from mouth_sealed_shape_relax import Meridian


def witness(path):
    check_input(path)
    with np.load(path) as z:
        if z['radius']!=4:
            raise ValueError('Expected radius 4')
        sh=Meridian(4.,z['shape_params'],float(z['shape_extra']))
        co=z['registration']; y=z['labels']
        na=(y.shape[1]-1)//2; nt=y.shape[0]-1
        raw=z['order']
        grid=(np.linspace(0,1,raw.shape[0]),np.linspace(0,1,raw.shape[1]),raw)
        ext=(z['outer_t'],z['outer_x'],z['outer_order'])
    p=CollarPair(sh,co,na,nt,quad=4)
    center=np.array(material_work(p,y,grid,ext))
    vals=y.reshape(-1,2)[p.id]
    inverse=np.stack([np.einsum('ki,kij->kj',g,vals) for g in (p.E,p.T)],axis=-1)
    S=sh.neck+(p.B*vals[:,:,0]).sum(axis=1)
    b,bw,_=sh.profile(p.z); slope=legval(2*p.z-1,co); pp=np.exp(slope)
    r=b*np.exp(-p.eta*pp); dr=-r*pp
    rt=r*(sh.length*bw/b - p.eta*pp*2*legval(2*p.z-1,legder(co)))
    current=np.stack((dr,rt,np.zeros_like(r),np.full_like(r,sh.length)),axis=-1).reshape(-1,2,2)
    inverse[:,1,:]*=-1
    F=current@np.linalg.inv(inverse)
    singular=np.linalg.svd(F,compute_uv=False)
    J=np.linalg.det(F)*(r/S)**2
    low=np.minimum(singular[:,1],r/S); high=np.maximum(singular[:,0],r/S)
    occupancy=np.zeros_like(r)
    for side,mesh in [(p.eta>0,grid),(p.eta<0,ext)]:
        interp=RegularGridInterpolator(mesh[:2],mesh[2])
        occupancy[side]=interp(np.column_stack([p.z[side],p.coordinate[side]]))
    bands={}
    for label,mask in [('ordered_chi_gt_half',occupancy>.5),('receiver_chi_gt_01',(p.eta<0)&(occupancy>.1)),
                       ('partial_chi_001_to_01',(occupancy<=.1)&(occupancy>.001)),
                       ('melted_chi_eq_zero',occupancy==0),('all',np.ones_like(r,dtype=bool))]:
        if np.any(mask):
            bands[label]=[float(low[mask].min()),float(high[mask].max()),float(J[mask].min()),float(J[mask].max())]
    print(json.dumps(dict(kind='bands',mesh=[na,nt],shape=sh.params.tolist(),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                          depth=sh.length,exit=sh.exit,central_energy=center.tolist(),stretch_min_max_j_min_max=bands)),flush=True)
    # Signed work using two smooth field variations that vanish on the ordered rim
    # and both ends. No solid–liquid extension or receiver label is changed.
    et,tt=np.meshgrid(p.a,p.t)
    shin=et>0
    h=np.where(shin,(np.sin(2*np.pi*et)*np.sin(np.pi*tt))**2,0)
    for mode in range(4):
        basis=h if mode<2 else h*np.cos(2*np.pi*tt)
        direction=np.zeros_like(y)
        direction[:,:,mode%2]=basis
        direction[~p.mask]=0
        for step in [.00002,.00001]:
            plus=np.array(material_work(p,y+step*direction,grid,ext))
            minus=np.array(material_work(p,y-step*direction,grid,ext))
            print(json.dumps(dict(kind='signed_compact_work',mode=mode,step=step,
                                  force=(minus-plus).sum()/(2*step),pos=(plus-center).sum(),
                                  neg=(minus-center).sum())),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state',type=Path)
    witness(parser.parse_args().state)
