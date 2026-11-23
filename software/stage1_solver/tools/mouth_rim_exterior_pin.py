"""Eulerian versus convected prescribed exterior order on the same axial bank strips."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from mouth_closed_rim_surface_review import quadrature
from mouth_disk_pillbox import dot
from mouth_last_ring_solid import interpolation
from mouth_poly_shallow_release import CUT
from mouth_rim_bank_work import Bank
from mouth_rim_disk_work_review import Disk
from mouth_sealed_shape_relax import DELTA
from mouth_supercore_joined_front import ORDER_INERTIA as I,check_guard,emit,file_sha


def gradient(model, values, x, z):
    ax=model.rr[1];zz=model.z[1]
    i=np.clip(ax.searchsorted(x,side='right')-1,0,len(ax)-2)
    j=np.clip(zz.searchsorted(z,side='right')-1,0,len(zz)-2)
    dx=np.diff(ax)[i];dz=np.diff(zz)[j]
    s=(x-ax[i])/dx;t=(z-zz[j])/dz
    # The input values have axes [case, physical depth node, outer radial node].
    lo=values[:,j,i];hi=values[:,j+1,i]
    da=values[:,j,i+1]-lo;db=values[:,j+1,i+1]-hi
    R,Rz=model.shape.radius(z)
    dr=((1-t)*da+t*db)/(dx*(28-R))
    return np.stack([dr,np.zeros_like(dr),(hi-lo+s*(db-da))/dz-dr*((1-x)*Rz)],axis=-1)


def pin_slopes(model,r,z,eta,side):
    # The receiver solves its own order table below the disk. At x<0 its
    # implementation explicitly clamps to the wall value of that table.
    ztest=z+side*1.e-9*eta
    b,slope=model.shape.radius(ztest)
    xx=(r-b)/(CUT-b)
    v,cx,cz=interpolation(*model.terminal_profile.side,ztest,np.clip(xx,0,1))
    inner=cz+np.where(xx<0,0.,cx)*(xx-1)*slope/(CUT-b)
    tail=expit(-ztest/DELTA)
    return eta*np.where(xx>1,-tail*(1-tail)/DELTA,inner)


def projection_energy(gp,gc,n):
    return -.5*I*n*dot(gp,gc)**2/(n+I*dot(gc,gc))


def audit(bank,n,baseline):
    m=bank.disk.model;sh=bank.disk.sh
    axes=m.rr[1];depth=m.z[1]
    xx,zz=np.meshgrid(axes,depth)
    b,_=sh.radius(zz);rr=b+xx*(28-b)
    eta=bank.motion(rr,zz)[0]
    profile=m.terminal_profile.outside
    c0=profile(rr,zz,False)
    # Original receiver data are independent of angle and equal in the two
    # histories on this bank. Make that a measured prerequisite for cancellation.
    for interp in bank.disk.reader.fields[1,1]:
        err=np.max(abs(interp.values[...,2]-c0.T[:,None,:]))
        if err>1e-12:raise ValueError(f'Exterior order differs across histories: {err}')
    signs=np.array([-1.,1.])
    slopes=np.array([pin_slopes(m,rr,zz,eta,s) for s in signs])
    eps=1.e-6
    steps=np.array([-eps,eps,-eps/2,eps/2])
    inc=np.array([(profile(rr,zz+h*eta,False)-c0)/h for h in steps])
    maxnodal=np.max(abs(inc-np.tile(slopes,(2,1,1))),axis=(1,2))
    lo,peak,hi=bank.strip
    zas,wz=quadrature(np.unique(np.r_[lo,peak,hi,depth[(depth>lo)&(depth<hi)]]),n)
    x,wx=quadrature(np.unique(np.r_[axes[axes<bank.edge],bank.edge]),n)
    an,wa=quadrature(m.a,n)
    x,a=np.meshgrid(x,an,indexing='ij')
    wx=wx[:,None]*wa
    cases=np.vstack([slopes,inc])
    weak=np.zeros(len(cases));signed=np.zeros(len(steps))
    for zc,dz in zip(zas,wz):
        b,_=sh.radius(zc);r=b+x*(28-b)
        _,ur,uz=bank.motion(r,zc)
        ur,uz=ur.ravel(),uz.ravel()
        weight=4*np.pi*dz*(28-b)*(wx*r**2*np.sin(a)).ravel()
        grad=gradient(m,cases,x.ravel(),zc)
        for sign,(state,*_) in zip((1.,-1.),bank.physical(x,a,zc)):
            gf,gp,gc,n0,swirl,chi=state
            pi=I*n0*dot(gp,gc)/(n0+I*dot(gc,gc))
            vel=gp-(I*dot(gp,gc)/(n0+I*dot(gc,gc)))[:,None]*gc
            weak+=sign*(np.einsum('kqi,qi->kq',grad,vel)@(weight*pi))
            for j,h in enumerate(steps):
                det=1+h*uz
                def pull(g):
                    out=g.copy()
                    out[:,2]/=det;out[:,0]-=h*ur*out[:,2]
                    return out
                # Direct effect of changing the boundary-order values with the
                # actual signed geometry while holding every other field fixed.
                p=pull(gp);ca=pull(gc);cb=pull(gc+h*grad[2+j])
                signed[j]-=sign*(weight@(det*(projection_energy(p,cb,n0)-projection_energy(p,ca,n0))))/h
    if np.max(abs(weak[4:]-weak[:2]))>5e-5:
        raise ArithmeticError('Independent nodal sampling and material derivative differ')
    if np.max(abs(signed[2:]-weak[4:]))>1e-3:
        raise ArithmeticError('Signed paired energy failed small-motion check')
    return dict(quadrature=n,signed_steps=steps.tolist(),profile_nodal_error=maxnodal.tolist(),
                independently_differentiated_nodal_work=weak[:2].tolist(),
                resampled_nodal_variation_work=weak[2:].tolist(),signed_full_energy_change=signed.tolist(),
                signed_external_order_cusp=float(weak[1]-weak[0]),
                prior_advected_residual=baseline,
                fixed_Eulerian_residual_one_sided=(baseline+weak[:2]).tolist(),
                centered_fixed_minus_traction=baseline+.5*(weak[0]+weak[1]))


def run(a):
    obj=Disk(a.state,a.wall,a.fields)
    check_guard(a.previous)
    rows=[json.loads(s) for s in a.previous.read_text().splitlines()]
    if rows[0]['field']!=obj.source_hash['field']:raise ValueError('Mismatched field')
    emit(kind='source',**obj.source_hash,bank_work=file_sha(a.previous))
    for row in rows:
        if row.get('bank')!='outside' or row['radial_cells']!=8:continue
        bank=Bank(obj,1,row['depth_knots'],row['normalized_radial_bound'])
        emit(kind='prescribed_profile',strip=row['strip'],**audit(bank,row['quadrature'],row['native_minus_wall']))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('state','wall','fields','previous'):parser.add_argument(key,type=Path)
    run(parser.parse_args())
