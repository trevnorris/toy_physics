"""Resolve the quadrature precision of the two extensive reservoir order energies at signed depth.

The carrier, chi and shape are the saved held state. The extension has no
arguments in either order term. Its already measured two signed banks are
retained when substituting the higher-precision difference of these terms.
"""
import argparse
from itertools import product
from pathlib import Path

import numpy as np

from mouth_continuous_depth_mesh import axes_for_depth
from mouth_continuous_shape_relax import ShapeTrial
from mouth_cp10_field_shape import Bath, BITS, GAUSS, corners
from mouth_exact_checkpoint_control import digest,read
from mouth_r4_barrier_gradient import evidence
from mouth_r4_eulerian_release import Driver
from mouth_whole_r4_relax import atomic_json


def extended(bath,held,length,ref):
    dtype=np.result_type(length,np.longdouble)
    def ar(x):
        return np.asarray(x,dtype=dtype)
    L=ar(length)
    nr=len(bath.model.rr[0])-1
    nz=len(bath.model.z[1])-nr-1
    s=ar(np.linspace(0.,1.,nr+1))
    v=ar(np.linspace(0.,1.,nz+1))
    width=ar(.12)
    orig=ar(ref)-width
    end=L-width
    lower=orig*s/(1+np.expm1(24*np.log(orig/end))*s**24)**(ar(1)/24)
    near=end+width*v
    upper=L+width*v
    far=L+width+(22-L-width)*(np.sinh(2*s)/np.sinh(ar(2)))**2
    zaxis=(ar(bath.model.z[0]),np.r_[lower[:-1],near],np.r_[upper[:-1],far])
    coef=ar(held.reference.c); knots=ar(held.reference.x)
    m,asym=ar(held.params[2:])
    e=ar(1.5)*np.expm1(ar(held.params[1]))
    extra=ar(held.extra)
    gx=ar(GAUSS); kappa=ar(bath.kappa);well=ar(bath.well)
    values=np.zeros(2,dtype=dtype)
    # Every arithmetic operand in both the differential volume and the order
    # strain is promoted upstream of multiplication, exponentials and division.
    for (slab,side), ids in bath.model.ids.items():
        zs=zaxis[slab]; xs=ar(bath.model.rr[side]); theta=ar(bath.model.a)
        dx=np.diff(xs)[:,None,None]; da=np.diff(theta)[None,:,None]; dz=np.diff(zs)[None,None,:]
        far_trace=1/(1+np.exp(zs/ar(bath.width)))
        if slab==2:far_trace[-1]=0
        chi=corners(np.where(bath.fixed[ids],far_trace,ar(bath.blocks[slab,side][2])))
        for pt in product(gx,repeat=3):
            p,a,w=pt
            x=xs[:-1,None,None]+p*dx
            angle=theta[None,:-1,None]+a*da
            t0=(zs[None,None,:-1]+w*dz)/L
            t=np.where(t0.real<0,0,np.where(t0.real>1,1,t0))
            u=t*ar(held.unit_length)
            idx=np.clip(np.searchsorted(knots.real,u.real,side='right')-1,0,len(knots)-2)
            pos=u-knots[idx]
            A,B,C,D=coef[:,idx]
            ref_radius=((A*pos+B)*pos+C)*pos+D
            ref_slope=(3*A*pos**2+2*B*pos+C)*ar(held.unit_length)
            v=t*(1-t)
            b=4*(ref_radius+e*t+4*m*v+4*asym*v*(2*t-1)+16*extra*v*v)
            bt=4*(ref_slope+e+4*m*(1-2*t)+4*asym*(-6*t*t+6*t-1)+32*extra*v*(1-2*t))
            bz=np.where((t0.real<0)|(t0.real>1),0,bt/L)
            y=(t0-ar(.06))/ar(.56)
            inside=(y.real>0)&(y.real<1)
            bump=np.where(inside,16*y*y*(1-y)**2,0)
            bump_z=np.where(inside,32*y*(1-y)*(1-2*y)/(ar(.56)*L),0)
            P=(3*np.cos(angle)**2-1)/2
            R=b+x*(28-b) if side else x*b
            rad_width=28-b if side else b
            shape_basis,shape_dx=((1-x,-1) if side else (x**3,3*x*x))
            shift=ar(bath.model.polar_shift)
            R=R+shift*bump*shape_basis*P
            Rx=rad_width+shift*bump*shape_dx*P
            Rtheta=-3*shift*bump*shape_basis*np.cos(angle)*np.sin(angle)
            Rz=bz*(1-x if side else x)+shift*bump_z*shape_basis*P
            weight=ar(.5*np.pi)*R*R*np.sin(angle)*Rx*dx*da*dz
            coords=ar(pt)
            fac=np.where(BITS,coords,1-coords)
            chi_q=np.einsum('i,ixat->xat',fac.prod(axis=1),chi)
            Dchi=[]
            for dim in range(3):
                basis=(2*np.array(BITS)[:,dim]-1)*fac[:,(dim+1)%3]*fac[:,(dim+2)%3]
                Dchi.append(np.einsum('i,ixat->xat',basis,chi))
            gr=Dchi[0]/dx/Rx
            gt=(Dchi[1]/da-Rtheta*gr)/R
            gw=Dchi[2]/dz-Rz*gr
            values[0]+=np.sum(weight*(kappa/2)*(gr*gr+gt*gt+gw*gw))
            values[1]+=np.sum(weight*well*(chi_q*(1-chi_q))**2)
    return values


def measure(sc,out):
    archive,_,previous=evidence(sc,'primary')
    other,_,secondary=evidence(sc,'secondary')
    with np.load(archive) as a,np.load(other) as b:
        if not np.array_equal(a['field'],b['field']) or not np.array_equal(a['x'],b['x']):
            raise ValueError('Depth controls require identical held fields and q')
        driver=Driver(sc,out.parent,1.e-8)
        driver.fields,driver.x=a['field'],a['x']
    shape=driver.shape(driver.x)
    trial=ShapeTrial(driver.base,driver.prep,driver.x,8,24,6,3,2,precondition=False,
                     axes=axes_for_depth(shape.length,driver.reference_depth))
    bath=Bath(trial,driver.base,driver.fields,driver.reference_depth)
    L=np.longdouble(shape.length)
    fixed=np.asarray(previous['modes']['log_depth']['analytic'][3:5])
    if not np.array_equal(fixed,secondary['modes']['log_depth']['analytic'][3:5]):
        raise ValueError('The analytic order-work components differ between gauges')
    zero_std=bath.energy(np.zeros(6))[3:5]
    zero_hp=extended(bath,shape,L,driver.reference_depth)
    imag=np.clongdouble(1j)*np.longdouble(1e-24)
    deriv=-extended(bath,shape,L*np.exp(imag),driver.reference_depth).imag/imag.imag
    record=dict(label='Common field, matched original signed depth grids; order-gradient and well only',
                schemas=['order_gradient','double_well'],longdouble_epsilon=float(np.finfo(np.longdouble).eps),
                analytic=fixed.tolist(),independent_80bit_analytic=deriv.astype(float).tolist(),
                energy_anchor_delta=(zero_hp-zero_std).astype(float).tolist(),banks={})
    if np.finfo(np.longdouble).eps>=np.finfo(float).eps or np.max(np.abs(deriv-fixed))>1.e-8 or np.max(abs(zero_hp-zero_std))>2.e-9:
        raise ArithmeticError(f'Extended-precision independent evaluator differs at anchor: {record}')
    for h in (4.e-5,2.e-5,1.e-5,5.e-6,2.5e-6):
        pairs={}
        for sign in (-1,1):
            shift=sign*h
            native=bath.energy(shift*np.eye(6)[0])[3:5]
            hp=extended(bath,shape,L*np.exp(np.longdouble(shift)),driver.reference_depth)
            pairs[str(sign)]=dict(native_work=((zero_std-native)/shift).astype(float).tolist(),
                                  extended_work=((zero_hp-hp)/shift).astype(float).tolist())
        fd=lambda name:sum(np.array(p[name]) for p in pairs.values())/2
        std=fd('native_work')
        high=fd('extended_work')
        record['banks'][str(h)]=dict(signed=pairs,native=std.tolist(),extended=high.tolist(),
            native_minus_analytic=(std-fixed).tolist(),extended_minus_analytic=(high-fixed).tolist(),
            correction=float((high-std).sum()))
    atomic_json(out,record)
    for h,item in record['banks'].items():
        print(h,item['native_minus_analytic'],item['extended_minus_analytic'],item['correction'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('out',type=Path)
    a=p.parse_args();measure(a.scratch.resolve(),a.out.resolve())
