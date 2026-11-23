"""Signed polar wall work below the resolved terminal band with both baths rematched.

A smooth quadrupolar displacement has exactly zero angular integral. For
meridional solid and order fields the first-order static material work in this
test vanishes. Compute the remaining carrier work on the same six-domain state.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_joined_front_audit import Front
from mouth_supercore_front import check_guard
from mouth_wound_carrier import straight_vortex
from mouth_live_waist_fields import ORDER_INERTIA


def bump(z, L):
    t=(z/L-.06)/.56
    inside=(t>0)&(t<1)
    h=np.where(inside,16*(t*(1-t))**2,0.)
    dh=np.where(inside,32*t*(1-t)*(1-2*t)/(.56*L),0.)
    return h,dh


def deformed_work(m, state, eps):
    wt,r,theta,z,f0,_,s,_,live,_,_=m.D.T
    L=m.shape.length
    h,dz=bump(z,L)
    active=h>0
    if np.any(active&(live!=0)):
        raise ArithmeticError('The test overlaps the angularly mobile terminal order band')
    ref=m.ref[active]
    wt,r,theta,z,f0,s,h,dz=[v[active] for v in (wt,r,theta,z,f0,s,h,dz)]
    val=m.V[active]@state[0]
    gv,gp=(np.column_stack([g[active]@u for g in m.G]) for u in state[:2])
    gc=ref[:,1:]
    if np.max(np.abs(m.V[active]@state[2]))>1.e-12:
        raise ArithmeticError('Nonzero polar order correction below its finite band')
    b,bw=m.shape.radius(z)
    inner=r<b
    width=np.where(inner,b,28.-b)
    x=np.where(inner,r/b,(r-b)/(28.-b))
    k=np.where(inner,x**3,1-x)
    kx=np.where(inner,3*x*x,-1.)
    xw=-bw*np.where(inner,x,1-x)/width
    c=np.cos(theta)
    P=(3*c*c-1)/2
    Pa=-3*c*np.sin(theta)
    rr=r+eps*k*h*P
    rx=1+eps*kx/width*h*P
    ra=eps*k*h*Pa
    rw=eps*(kx*xw*h+k*dz)*P
    if np.min(rx)<=0 or np.min(rr)<=0:
        raise ArithmeticError('Inverted signed variation')
    vtx=straight_vortex(m.winding,28.)
    f1,_=vtx(rr*np.sin(theta))
    def grad(g):
        radial=g[:,0]/rx
        return np.column_stack((radial,(r*g[:,1]-ra*radial)/rr,g[:,2]-rw*radial))
    def density(f,v,p,cg):
        n=(f*val)**2
        cs=np.sum(cg*cg,axis=1)
        den=n+ORDER_INERTIA*cs
        dot=np.sum(p*cg,axis=1)
        phase=.5*n*(np.sum(p*p,axis=1)-ORDER_INERTIA*dot*dot/den)
        y=val*val
        pot=f**10*(y-1)**2*(y**3+2*y*y+3*y+4)/20
        return .5*f*f*np.sum(v*v,axis=1)+phase+pot
    delta=wt@(((rr/r)**2*rx)*density(f1,grad(gv),grad(gp),grad(gc))-density(f0,gv,gp,gc))
    # The straight-vortex factorization integrates by parts independently on
    # both sides of the sealed lateral face; a tilted normal has a polar part.
    zs=m.z[1]
    angles=m.a
    quad=(1+np.array([-1.,1.])/np.sqrt(3))/2
    face_delta=0.
    for side in (0,1):
        face=m.ids[1,side][0 if side else -1]
        uv=state[0,face]
        for a in quad:
            angle=angles[:-1,None]+a*np.diff(angles)[:,None]
            sn=np.sin(angle); co=np.cos(angle)
            for d in quad:
                w=zs[None,:-1]+d*np.diff(zs)[None,:]
                hh,_=bump(w,L)
                rad,_=m.shape.radius(w)
                r1=rad+eps*hh*(3*co*co-1)/2
                r1a=-eps*hh*3*co*sn
                ff,df=vtx(r1*sn)
                fo,do=vtx(rad*sn)
                v=((1-a)*((1-d)*uv[:-1,:-1]+d*uv[:-1,1:])+
                   a*((1-d)*uv[1:,:-1]+d*uv[1:,1:]))
                weight=.5*(1-2*side)*np.pi*np.diff(angles)[:,None]*np.diff(zs)[None,:]*(v*v-1)*sn
                face_delta+=np.sum(weight*(r1*r1*ff*df*(sn-r1a*co/r1)-rad*rad*fo*do*sn))
    return np.array([delta,face_delta])


def run(args):
    check_guard(args.fields)
    row=next(q for q in map(json.loads,(args.fields.parent/'stdout').read_text().splitlines())
             if q.get('kind')=='preparation')
    sha=hashlib.sha256(args.state.read_bytes()).hexdigest()
    if sha!=row['source_sha256']:
        raise ValueError('Mismatched material and reservoir fields')
    op=Front(args.state,tuple(row['grid']),row['winding'])
    m=op.joined(op.sh,False)
    with np.load(args.fields) as a:
        for key in ('fields','rest'):
            state=a[key]
            for step in (.0004,.0002):
                lo=deformed_work(m,state,-step)
                hi=deformed_work(m,state,step)
                force=(lo-hi)/(2*step)
                print(json.dumps(dict(kind='joined_angular_work',state_sha256=sha,
                    fields_sha256=hashlib.sha256(args.fields.read_bytes()).hexdigest(),
                    grid=row['grid'],winding=row['winding'],state=key,h=step,
                    depth_support=[.06,.62],volume_force=force[0],wall_force=force[1],
                    net_force=sum(force),signed_minus=lo.tolist(),signed_plus=hi.tolist()),
                    allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('fields',type=Path)
    run(p.parse_args())
