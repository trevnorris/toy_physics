"""Prepare the draining carrier and its passive conversion layer at fixed exterior bias.

The reference R=1 bias is inherited from its archived intake-0.25 calculation.
Physical axial depth, slab and interface thickness, reservoir size, and material
parameters are common. Each new aperture selects its own inward current.
"""
import argparse
from pathlib import Path
import json

import numpy as np
from scipy.optimize import brentq

from mouth_carrier_depth_selection import DepthDomain
from mouth_carrier_radius_family import fixed_head
from mouth_free_order_diffuse import NU, Order, positive_factor, read
from mouth_free_order_residual_check import DifferenceOrder
from mouth_order_wave_pulse import ARCHIVE, ROOT, digest, emit


def domain(radius, mesh):
    return DepthDomain(radius,1.,mesh,9.,10.,1.)


def reference():
    path=ROOT/'mouth_flared_cubic_n34'/'field.npz'
    if digest(path)!=json.loads(ARCHIVE.read_text())['input_carriers_sha256'][str(path)]:
        raise ValueError('Unrecognized carrier reference')
    return json.loads((path.parent/'stdout').read_text().splitlines()[0])['bias']


def main(a):
    dom=domain(a.radius,a.mesh)
    bias=reference()
    carrier,trial=fixed_head(dom,bias)
    emit(test='carrier',radius=a.radius,mesh=a.mesh,bias=bias,**trial)
    delta=.04
    model=DifferenceOrder(dom,carrier,delta)
    order,stats=model.minimize(a.seed)
    original=Order(dom,carrier,delta)
    resid,jac=original.measure(order)
    fac=positive_factor(jac)
    err=float(np.max(abs(resid)))
    corr=float(np.max(abs(fac.solve(resid)))) if fac is not None else None
    # At the wall the native Q2 coordinates beyond the waist are (p=w,v=0).
    pg,vg,ids=dom.blocks[2]
    wall=order[ids.reshape(len(pg)-1,len(vg)-1,3,3)[:,0,:,0]]
    def rim(w):
        i=np.clip(np.searchsorted(pg,w)-1,0,len(pg)-2)
        s=2*(w-pg[i])/(pg[i+1]-pg[i])-1
        return np.dot([s*(s-1)/2,1-s*s,s*(s+1)/2],wall[i])-.5
    cuts=[brentq(rim,lo,hi) for lo,hi in zip(pg[:-1],pg[1:]) if rim(lo)*rim(hi)<0]
    emit(test='order',seed=a.seed,radius=a.radius,mesh=a.mesh,delta=delta,nu=NU,
         wall_half_order_depths=cuts,original_residual=err,newton_correction=corr,**stats)
    if corr is None or corr>1e-8 or err>1e-9 or trial['feed']<=0:
        raise ValueError('Carrier or stationary layer failed its independent check')
    dest=a.outdir/'state.npz'
    np.savez_compressed(dest, order=order, carrier=carrier, points=dom.points)
    emit(test='saved',path=str(dest),sha256=digest(dest))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--seed',type=float,default=1.55)
    p.add_argument('--outdir',type=Path,required=True)
    main(p.parse_args())
