"""Evaluate independent fixed-state virtual work on a saved live order profile."""
import json
import sys
import hashlib
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from closed_wave import Shape, arrays, mode
from shape_check import duct

def energy(d,chi):
    xy,ids,tri,grad,quads=d
    fields=chi[tri]
    g=np.einsum('ck,ckj->cj',fields,grad)
    density_squared=np.einsum('cd,cd->c',g,g)
    E=0.
    for N,r,weight in quads:
        c=fields@N
        E+=4*np.pi*np.dot(weight,.04**2*density_squared+(c*(1-c))**2)
    return E*75.

def measure(source,r,tag,n):
    inp=np.load(source); raw=inp['order']; target=float(inp['frequency'])
    rec=json.loads(Path('research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json').read_text())['training'][str(r)]
    state={'broad':rec['states'][-2], 'tight':rec['extension_states'][-1],'mid':rec['extension_states'][1]}[tag]
    X=np.asarray(state['coordinates']); grid=(np.linspace(0,1,len(raw)),np.linspace(0,1,raw.shape[1]))
    t,x=np.meshgrid(np.linspace(0,1,2*n+1),np.linspace(0,1,n+1),indexing='ij'); tx=np.c_[t.ravel(),x.ravel()]
    chi=RegularGridInterpolator(grid,raw)(tx)
    q=arrays(Shape(r,X),n)
    f,zfrac,_,wlist=mode(q,chi)
    nearest=np.argmin(abs(np.asarray(wlist)-target))
    print(json.dumps(dict(type='basis',radius=r,state=tag,res=n,source=str(source),data_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),primary_pole=target,independent_clamped=wlist,nearest=int(nearest))),flush=True)
    if nearest!=0: f,zfrac,_,wlist=mode(q,chi,index=nearest)
    h=.002
    for direction,i,scale in [('waist',2,1/r),('axial',0,1.)]:
        vals=[]
        for sign in (-1,1):
            at=X.copy();at[i]+=scale*sign*h
            shape=Shape(r,at); dat=arrays(shape,n)
            closed=mode(dat,chi,index=nearest)[0]; order=energy(dat,chi); flow=duct.solve(shape,n//2)
            vals.append([closed,order,flow])
        fd=((np.asarray(vals[1])-vals[0])/(2*h)).tolist()
        print(json.dumps(dict(type='work',R=r,point=tag,direction=direction,n=n,fixed_action_closed_wave=fd[0],resolved_core_order_at_transported_chi=fd[1],linear_sealed_carrier=fd[2],axial_wave_identity=-f*zfrac)),flush=True)
        if direction=='axial': assert abs(fd[0]+f*zfrac)<3e-4
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--radius', type=int, choices=(1, 2), required=True)
    parser.add_argument('--state', choices=('broad', 'mid', 'tight'), required=True)
    parser.add_argument('--mesh', type=int, required=True)
    args = parser.parse_args()
    measure(args.source, args.radius, args.state, args.mesh)
