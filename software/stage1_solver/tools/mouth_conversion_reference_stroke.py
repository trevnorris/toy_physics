"""Test a conversion-dependent rest-reference stroke against passive contact renewal.

For g=s−p and y=1−H(chi), the independent material coefficient beta multiplies
−g D_t y in the reference-rate law. Chemical energy pays beta (t·g) D_t y.
For the same trial the contact contribution to the wall order variation is
k|g|²(1/2−beta chi H') per physical area. The stationary carrier is prescribed.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import linalg

from mouth_bend_renewal import ROOT, volume_cost, wall_data
from mouth_bend_static_shear import ANGULAR, archive_inputs
from mouth_curved_phase import quadratic
from mouth_supplied_contact_work import contact_readout,verified
from mouth_order_wave_pulse import digest,emit


def turnover(geo,chi,psi,edges,keep,quad):
    count=np.flatnonzero(keep)[-1]+1
    parent,vg,cells=geo.order_blocks[2]
    width=np.diff(edges)[:count,None]
    t,weight=leggauss(quad)
    t=(1+t)/2
    z=edges[:count,None]+width*t
    idx=np.searchsorted(parent,z,side='right')-1
    nodes=cells[idx*(len(vg)-1)][:,:,[0,3,6]]
    height=np.diff(parent)[idx]
    basis,derivative=quadratic(2*(z-parent[idx])/height-1)
    def sample(values):
        return (basis*values[nodes]).sum(axis=-1),(derivative*values[nodes]).sum(axis=-1)*2/height
    c,cz=sample(chi)
    f,fz=sample(psi)
    radius,rz=geo.radius(z)
    arc=np.hypot(1.,rz)
    speed_z=np.imag(fz/f)/arc**2
    ydot=-6*c*(1-c)*cz*speed_z
    line=width*weight/2
    basis=np.broadcast_to(np.stack([1-t,t],axis=-1),(*z.shape,2))
    grad=np.broadcast_to([-1.,1.],basis.shape)/width[:,:,None]
    ids=np.arange(count)[:,None]+np.arange(2)
    use=keep[ids].any(axis=1)
    if np.any(c[use]<=0) or np.any(speed_z<=0):
        raise ValueError('Invalid contact support or wall speed')
    measure=line/speed_z
    projection_weight=np.divide(measure*ydot,c,out=np.zeros_like(c),where=use[:,None])
    def gram(w,p,q):
        return np.einsum('eq,eqi,eqj->eij',w,p,q)
    mass=gram(measure,basis,basis)
    def project(w,mode):
        local=np.linalg.solve(mass,gram(w,basis,mode))
        matrix=np.zeros((2*count,len(edges)))
        row=(2*np.arange(count)[:,None]+np.arange(2))
        np.add.at(matrix,(row[:,:,None],ids[:,None,:]),local)
        return matrix[:,keep]
    Gc=np.vstack([project(projection_weight/radius**2,basis),project(projection_weight/(2*radius*arc),grad)])
    def work(m):
        field=np.zeros(len(edges))
        field[keep]=m
        tm=(basis@field[ids][:,:,None])[...,0]/radius**2
        ta=(grad@field[ids][:,:,None])[...,0]/(2*radius*arc)
        ext2=np.divide(tm*tm+2*ta*ta,c*c,out=np.zeros_like(c),where=use[:,None])
        area=ANGULAR*line*radius**2*arc
        return np.sum(area*c*ext2*ydot)
    return Gc,work


def main(a):
    records=verified(a.source)
    setup=next(r for r in records if r['test']=='calibration')
    mesh=setup['mesh']
    geo,_,carrier=archive_inputs(Path(setup['order_source']),ROOT/f'bend_static_shear_delcarrier{mesh}_01')
    file=Path(records[0]['source'])
    if digest(file)!=records[0]['sha256']:
        raise ValueError('Carrier changed')
    with np.load(file) as data:
        psi=data['psi']
    for record in records:
        if record.get('multiple') not in a.levels or record['test']!='supplied_order':
            continue
        file=Path(record['path'])
        if digest(file)!=record['sha256']:
            raise ValueError('Supplied order changed')
        with np.load(file) as data:
            chi=data['chi']
        edges,keep,bulk=volume_cost(geo,chi,a.sub,4,1e-9)
        contact,load,B,G,L,tau,_=wall_data(geo,chi,psi,edges,keep,4)
        Gc,work=turnover(geo,chi,psi,edges,keep,4)
        base=bulk+contact
        scale=1/np.sqrt(base.diagonal())
        inverseB=scale[:,None]*linalg.solve(scale[:,None]*base*scale,scale[:,None]*B,assume_a='pos')
        for beta in a.beta:
            rate=G/tau-beta*Gc
            history=linalg.solve(L,rate)
            matrix=base+B@history
            m=scale*linalg.solve(scale[:,None]*matrix*scale,scale*load[:,1])
            ref=history@m
            decay=np.min(linalg.eigvals(L+rate@inverseB).real)
            passive=contact_readout(geo,chi,psi,m,keep,edges,tau)
            act=beta*work(m)
            net=act+passive['order_work_on_contacts']
            # contact_readout infers slip work with passive reference evolution.
            # The active rate changes the true input through material slip by −act.
            slip=passive['slip_work_from_material_balance']-act
            balance=np.max(abs(L@ref-rate@m))
            emit(test='stroke',mesh=mesh,sub=a.sub,multiple=record['multiple'],head=record['head'],beta=beta,
                 active_reference_power=act,order_modulus_power=passive['order_work_on_contacts'],
                 net_conversion_power=net,renewal_power=passive['renewal'],material_slip_work=slip,
                 net_per_renewal=net/passive['renewal'],min_decay_rate=decay,dc_residual=balance,
                 contact_deep_holder=passive['generalized_contact_holder'],
                 two_term_deep_holder=passive['generalized_contact_holder']+carrier[1,1],
                 contact_force_on_material_segment=passive['segment_force_on_material_x'],
                 contact_force_on_wall_segment=passive['segment_force_on_wall_x'],
                 note='Per b² for power, per b for force. Trial beta is per converted constituent fraction. Straight-prepared entry.')
            if balance>2e-7 or any(not np.isfinite(v) for v in (decay,act,net)):
                raise ArithmeticError('Invalid stationary response')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--sub',type=int,default=3)
    p.add_argument('--levels',nargs='+',type=float,default=[.2,.5,.8])
    p.add_argument('--beta',nargs='+',type=float,default=[0.,.5,1.,2.,4.])
    main(p.parse_args())
