"""Static order work for a core stress-free at the shallow starting geometry, bonded to the collar."""
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legder, leggauss, legval

from mouth_sealed_shape_relax import Meridian


def _cut(shape, c, t):
    b, db_dw, _ = shape.profile(t)
    p = np.exp(legval(2*t-1, c))
    cut = b*np.exp(-.5*p)
    dlogcut_dt = shape.length*db_dw/b - p*legval(2*t-1, legder(c))
    return b, cut, dlogcut_dt


def core_stretches(shape, c, shape0, c0, t, r):
    """Return F_rr=F_hoop, F_rZ, F_wZ and det(F); r is the *current* radius."""
    _, cut, dt = _cut(shape, c, t)
    _, initial, dt0 = _cut(shape0, c0, t)
    q = cut/initial
    k = shape.length/shape0.length
    return q, r*(dt-dt0)/shape0.length, k, q**3*k


def _energy(q, shear, k):
    # Evaluate the stored energy without subtracting order-one terms near identity.
    alpha, eta = np.log(q), np.log(k)
    return (.5*(3*(np.expm1(2*alpha)-2*alpha) + np.expm1(2*eta)-2*eta + shear**2)
            + .25*np.expm1(3*alpha+eta)**2)


def core_source(shape, c, shape0, c0, t, x):
    """Return dE_el/dχ for the interior nodal order mesh on x=r/b(t), 0≤t≤1."""
    grid = np.zeros((len(t),len(x)))
    if t[0] != 0 or t[-1] != 1 or x[0] != 0 or x[-1] != 1:
        raise ValueError("Expected the full core order chart")
    if shape.neck != shape0.neck:
        raise ValueError("Different preparation and current entrance")
    dt, dx = np.diff(t)[:,None], np.diff(x)[None,:]
    if np.min(dt)<=0 or np.min(dx)<=0:
        raise ValueError("Expected increasing material chart")
    nodes, weights = leggauss(5)
    for z,wz in zip((1+nodes)/2,weights/2):
        ts=t[:-1,None]+dt*z
        b,cut,log_dt=_cut(shape,c,ts)
        _,initial,log_dt0=_cut(shape0,c0,ts)
        q=cut/initial
        k=shape.length/shape0.length
        jac=q**3*k
        left=x[None,:-1]
        width=np.maximum(0.,np.minimum(cut/b,x[None,1:])-left)
        for y,wy in zip((1+nodes)/2,weights/2):
            xx=left+width*y
            r=b*xx
            shear=r*(log_dt-log_dt0)/shape0.length
            density=_energy(q,shear,k)/jac
            measure=4*np.pi*r*r*b*width*shape.length*dt*wz*wy*density
            u=width*y/dx
            for iz,vz in ((0,1-z),(1,z)):
                for ix,vx in ((0,1-u),(1,u)):
                    grid[iz:iz+len(t)-1,ix:ix+len(x)-1]+=measure*vz*vx
    if not np.all(np.isfinite(grid)) or np.min(grid)<0:
        raise ArithmeticError("Invalid bonded-core order source")
    return grid


def _reference_integral(shape,c,shape0,c0):
    # Independent integration over initial radii and depth, with a bilinear χ
    # times the element shape density. This uses squared stretches directly.
    z,w = leggauss(172)
    ts=(1+z)[:,None]/2
    y,wy=leggauss(11)
    b,cs,dlogc=_cut(shape,c,ts)
    _,cinit,dlogc0=_cut(shape0,c0,ts)
    old_r=cinit*(1+y)/2
    new_r=old_r*cs/cinit
    q=cs/cinit
    k=shape.length/shape0.length
    beta=new_r*(dlogc-dlogc0)/shape0.length
    j=q**3*k
    trace=3*q*q+k*k+beta*beta
    energy=(trace-4)/2-np.log(j)+.25*(j-1)**2
    measure=4*np.pi*shape0.length*(w[:,None]/2)*old_r**2*cinit*wy/2
    xx=new_r/b
    chi=.31+.13*xx+.09*ts-.03*xx*ts
    return float(np.sum(measure*energy)),float(np.sum(measure*energy*chi))


def check():
    from mouth_strained_wave import core_factors

    for radius in (1,2):
        directory=Path(f'_scratch/s11c/polycrystal_material01/r{radius}g2')
        inp=directory/'coeffs.npz'
        out=json.loads((directory/'outcome.json').read_text())
        digest=hashlib.sha256(inp.read_bytes()).hexdigest()
        if (out['exitCode'] or not out['limitsVerified'] or digest not in (directory/'stdout').read_text()):
            raise ValueError('Unvalidated starting collar')
        with np.load(inp) as data:
            x0,c0=data['coordinates'],data['ordered']
        shape0=Meridian(radius,x0)
        t,x=np.linspace(0,1,145),np.linspace(0,1,73)
        chi=.31+.13*x[None,:]+.09*t[:,None]-.03*x[None,:]*t[:,None]
        for sign in (0,-1,1):
            coords=x0.copy()
            coords[[0,2]]+=np.array([.003,-.002])*sign
            c=c0.copy()
            c[[0,1]]+=np.array([.005,-.002])*sign
            shape=Meridian(radius,coords)
            coated=core_source(shape,c,shape0,c0,t,x)
            simple,weighted=_reference_integral(shape,c,shape0,c0)
            val,wt=float(np.sum(coated)),float(np.sum(chi*coated))
            points=np.array([.031,.19,.48,.69,.96])
            radii=np.array([.14,.26,.3,.12,.17])
            q,sz,k,j=core_stretches(shape,c,shape0,c0,points,radii)
            wave_q,wave_d,wave_k=core_factors(shape,c,(shape0,c0),points)
            tangent_difference=max(np.max(abs(q-wave_q)),np.max(abs(sz-radii*wave_d)),abs(k-wave_k))
            if tangent_difference>1.e-10:
                raise ArithmeticError("Order and wave use different core deformation")
            err=max(abs(simple-val),abs(weighted-wt))
            print(json.dumps(dict(radius=radius,sign=sign,source_sha256=digest,total=val,
                                  reference_total=simple,weighted=wt,reference_weighted=weighted,
                                  wave_factor_difference=tangent_difference,smallest_nodal_source=float(np.min(coated)))),flush=True)
            if (sign==0 and (val!=0 or wt!=0) or sign and min(val,wt)<=0 or
                    err>1.e-11+2.e-6*max(simple,weighted)):
                raise ArithmeticError('Bonded core reference/direct quadratures disagree')


if __name__=='__main__':
    check()
