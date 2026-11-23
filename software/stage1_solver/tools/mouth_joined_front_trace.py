"""Read saved bidirectional front fields on each face of the physical radius-four exit."""
import argparse
import contextlib
import io
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator

from mouth_joined_front_audit import Front
from mouth_supercore_joined_front import emit,check_guard,file_sha,WIDTH,KAPPA,WELL,ORDER_INERTIA
from mouth_sealed_shape_relax import Meridian
from mouth_wound_carrier import straight_vortex


def face(model,order,fields,sign,quadrature):
    slab=1 if sign<0 else 2
    j=-2 if sign<0 else 0
    xs,angles,zs=model.rr[0],model.a,model.z[slab]
    xn,an,zn=np.meshgrid(xs,angles,zs[j:j+2] if j>=0 else zs[-2:],indexing='ij')
    interp=RegularGridInterpolator((np.linspace(0,1,len(order)),np.linspace(0,1,order.shape[1])),order)
    ref=interp((zn/model.shape.length,xn)) if sign<0 else model.terminal_profile(xn,zn-model.shape.length)
    ids=model.ids[slab,0][:,:,-2:] if sign<0 else model.ids[slab,0][:,:,:2]
    v,p,c=fields[:,ids]
    c=c+ref
    x0,a0=np.meshgrid(xs[:-1],angles[:-1],indexing='ij')
    dx,da=np.meshgrid(np.diff(xs),np.diff(angles),indexing='ij')
    dz=(zs[-1]-zs[-2]) if sign<0 else zs[1]-zs[0]
    b,slope=model.shape.radius(model.shape.length)
    if sign>0:slope=0.
    bits=np.array([[i,k,h] for i in (0,1) for k in (0,1) for h in (0,1)])
    def nodal(q):
        return np.stack([q[i:i+len(xs)-1,k:k+len(angles)-1,h] for i,k,h in bits],-1)
    u=np.array(list(map(nodal,(v,p,c))))
    z=1 if sign<0 else 0
    g,wt=leggauss(quadrature);g=(1+g)/2;wt=wt/2
    rows=[];vortex=straight_vortex(model.winding,28.)
    for x,wx in zip(g,wt):
        for a,wa in zip(g,wt):
            factors=np.where(bits,np.array([x,a,z]),1-np.array([x,a,z]))
            vv=factors.prod(axis=1)
            dr,gt,dw=np.stack([(2*bits[:,i]-1)*factors[:,(i+1)%3]*factors[:,(i+2)%3] for i in range(3)])
            t=(x0+x*dx);r=b*t;th=a0+a*da;s=r*np.sin(th)
            weight=4*np.pi*r*r*np.sin(th)*b*dx*da*wx*wa
            V,phi,chi=u@vv
            F,Fp=vortex(s);n=(F*V)**2
            G=np.stack((u@dr/(dx*b),u@gt/(da*r),u@dw/dz-slope*t*(u@dr)/(dx*b)),-1)
            gv,gp,gc=G
            gf=F[...,None]*gv+V[...,None]*Fp[...,None]*np.stack((np.sin(th),np.cos(th),np.zeros_like(th)),-1)
            dot=np.sum(gp*gc,axis=-1);gg=np.sum(gc*gc,axis=-1)
            vel=gp-(ORDER_INERTIA*dot/(n+ORDER_INERTIA*gg))[...,None]*gc
            current=n[...,None]*vel
            sq=lambda arr:(arr*arr).sum(axis=-1)
            swirl=(model.winding/s)**2
            mach=np.sqrt(sq(vel)+swirl)/n**2
            energy=.5*(sq(gf)+n*swirl+KAPPA*gg+(gp*current).sum(axis=-1))+(n**5-5*n+4)/20+WELL*(chi*(1-chi))**2
            stress=n*vel[...,2]**2+gf[...,2]**2+KAPPA*gc[...,2]**2-energy
            rows.append(np.column_stack([f.ravel() for f in (weight,n,mach,current[...,2],stress)]))
    return np.vstack(rows).T


def main(a):
    check_guard(a.fields)
    base=Front(a.state,tuple(a.grid),a.winding)
    with np.load(a.fields) as values:
        u=values['fields'];rest=values['rest'];v=values['sharp_amplitude'];p=values['sharp_phase']
    m=base.joined(base.sh,False)
    emit(kind='sources',state_sha256=file_sha(a.state),field_sha256=file_sha(a.fields),grid=a.grid)
    for sign in (-1,1):
        for quad in (8,16):
            weights,n,mach,j,t=face(m,base.c,u,sign,quad)
            *_,t0=face(m,base.c,rest,sign,quad)
            dense=(n>.5)&(mach<1.)
            pos=np.maximum(j,0)
            emit(kind='physical_disk',side=sign,quadrature=quad,signed=weights@j,
                 positive=weights@pos,neg=weights@np.minimum(j,0),minimum_mach=mach.min(),
                 dense_subsonic=weights[dense]@pos[dense]/(weights@pos),
                 raw_axial=weights@t,drive_axial=weights@(t-t0),drive_dense=weights[dense]@(t-t0)[dense])
    if a.static_work:
        sharp=base.sharp(base.sh)
        with contextlib.redirect_stdout(io.StringIO()):
            v0,p0,_=sharp.solve(0.)
        for i in range(5):
            delta=.0004;rows=[]
            for sign in (-1,1):
                dq=sign*delta*np.eye(5)[i]*np.array([1.,.25,.25,.25,.25])
                sh=Meridian(4.,base.sh.params+dq[:4],extra=base.sh.extra+dq[4])
                live=base.joined(sh,False);old=base.sharp(sh)
                rows.append([live.evaluate(u)[0],live.evaluate(rest)[0],old.amplitude(v,p)[0],old.amplitude(v0,p0)[0]])
            f=-(np.array(rows[1])-rows[0])/(2*delta)
            emit(kind='static_dynamic_work',coordinate=i,static_correction=f[1]-f[3],
                 bias_induced_correction=(f[0]-f[1])-(f[2]-f[3]),total_correction=f[0]-f[2])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path);p.add_argument('fields',type=Path)
    p.add_argument('--grid',type=int,nargs=3,required=True)
    p.add_argument('--winding',type=int,choices=(0,1),required=True)
    p.add_argument('--static-work',action='store_true')
    main(p.parse_args())
