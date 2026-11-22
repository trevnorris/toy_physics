"""Relax an internal receiver displacement at fixed physical shape, carrier and wave."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq, minimize_scalar

from mouth_polycrystal_material import (Meridian, RECORD, Side, archived, bump, exterior,
                                        fields, graph_values, stable_energy)


def diagnostic(radius, grade, directory):
    root = json.loads((directory/"outcome.json").read_text())
    raw = list(map(json.loads,(directory/"stdout").read_text().splitlines()))
    snapshot=directory/"coeffs.npz"
    fingerprint=hashlib.sha256(snapshot.read_bytes()).hexdigest()
    if (root['exitCode'] or root['childOutcome']['exitCode'] or not root['limitsVerified']
            or raw[-1]['sha256']!=fingerprint):
        raise ValueError("Unverified stationary collar")
    deep,chi,source = archived(radius,72,json.loads(RECORD.read_text()))
    if raw[0]['radius']!=radius or raw[0]['grade']!=grade:
        raise ValueError("Wrong material preparation")
    side=Side(radius,fields(deep,chi),-1,grade)
    with np.load(snapshot) as z:
        shape=Meridian(radius,z['coordinates'])
        c=z['receiver']
    r,rs,rz,h,j0,_=graph_values(side.a,side.t,shape,side.basis@c,side.deriv@c)
    L=shape.length
    z,zt=bump(side.t,.15,.85)
    b,ba=bump(-2*side.a,.04,.48)
    active=b*z>0
    g,gt,ga=(v[active] for v in (b*z,b*zt,-2*ba*z))
    t,r,rs,rz,h,j0=(v[active] for v in (side.t,r,rs,rz,h,j0))
    if min(side.chi[active])<=0: raise ArithmeticError("Positive order expected on support")
    volume=side.weight[active]/side.chi[active]
    dj=(-rs*gt-rz*ga)*h**2
    base_trace=rs**2+rz**2+2*h*h
    _, (wtab,xtab,chi)=exterior(shape,48 if grade==1 else 70,return_side_field=True)
    interp=RegularGridInterpolator((wtab,xtab),chi)
    def energy(e,live=False):
        ws,wz=e*ga,-L-e*gt
        jac=j0+e*dj
        if min(jac)<=0: raise ValueError('Inverted trial')
        density=stable_energy(base_trace+ws*ws+wz*wz,jac)
        if not live:
            return float(side.weight[active]@density)
        w=L*t+e*g
        width=shape.radius(w)[0]
        return float(volume@(interp(np.c_[w,(r-width)/(4.6-width)])*density))
    def derivative(e):
        jac=j0+e*dj
        return float(side.weight[active]@(e*ga*ga+(L+e*gt)*gt
                                          +(-1/jac+.5*(jac-1))*dj))
    limit=float(np.min(-j0[dj<0]/dj[dj<0]))
    bracket=(0., min(.77*limit,.12))
    root=brentq(derivative,*bracket)
    # Compute a genuine finite local lowering of the live-order energy, with
    # the receiver's Eulerian order held fixed along the path.
    live=minimize_scalar(lambda e:energy(e,True),bounds=bracket,method="bounded",
                         options={'xatol':1e-11})
    if not live.success: raise ArithmeticError(live.message)
    checks=[]
    for name,opt,live_flag in [('material_order',root,False),('eulerian_order',live.x,True)]:
        val=energy(opt,live_flag)
        check_step=min(2e-4,opt/8,(bracket[1]-opt)/8)
        checks.append(dict(order=name,amplitude=float(opt),energy_change=val-energy(0.,live_flag),
                           work_to_minus=energy(opt-check_step,live_flag)-val,
                           work_to_plus=energy(opt+check_step,live_flag)-val,
                           check_step=check_step,min_sampled_J_ratio=float(min(1+opt*dj/j0))))
    print(json.dumps(dict(radius=radius,grade=grade,source_sha256=fingerprint,
                          trial_interval=bracket,initial_frozen_derivative=derivative(0.),checks=checks),
                     allow_nan=False),flush=True)
    for row in checks:
        if (row['energy_change']>=0 or min(row['work_to_minus'],row['work_to_plus'])<=0 or
                row['check_step']<1.e-8 or row['min_sampled_J_ratio']<=0.):
            raise ArithmeticError('Finite lowering or interior directional minimum unresolved')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    p.add_argument('--grade',type=int,choices=(1,2),required=True)
    p.add_argument('--directory',type=Path,required=True)
    a=p.parse_args()
    diagnostic(a.radius,a.grade,a.directory)
