"""Independent W/J push-forward and paired order-work check on accepted R4 geometry."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_bonded_front_feedback import Relocated
from mouth_bonded_material_motion import Bonded
from mouth_joined_front_audit import Front
from mouth_supercore_front import check_guard


def interp_triangle(x,t,xs,ts,values):
    i=np.searchsorted(xs,x,side='right')-1;j=np.searchsorted(ts,t,side='right')-1
    i=np.clip(i,0,len(xs)-2);j=np.clip(j,0,len(ts)-2)
    a=(x-xs[i])/(xs[i+1]-xs[i]);b=(t-ts[j])/(ts[j+1]-ts[j])
    q0=values[j,i]*(1-a-b)[:,None]+values[j,i+1]*a[:,None]+values[j+1,i]*b[:,None]
    q1=values[j+1,i+1]*(a+b-1)[:,None]+values[j,i+1]*(1-b)[:,None]+values[j+1,i]*(1-a)[:,None]
    return np.where((a+b<1)[:,None],q0,q1)


def audit(args):
    for path in (args.state,args.motion,args.flow):check_guard(path)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    with np.load(args.motion) as a:
        assert a['joined_sha256'].item()==sha(args.flow) and a['source_sha256'].item()==sha(args.state)
        tabs={k:a[k].copy() for k in ('radial_fraction','depth_fraction','displacement')}
    props=next(r for r in map(json.loads,(args.flow.parent/'stdout').read_text().splitlines()) if r.get('kind')=='preparation')
    f=Front(args.state,tuple(props['grid']),1)
    sh=f.sh
    assert props['source_sha256']==sha(args.state)
    with np.load(args.flow) as a:solution=a['fields']
    L=sh.length;xs=tabs['radial_fraction'];ts=tabs['depth_fraction'];u=tabs['displacement']
    rel=Relocated(f,tabs)
    def chart(q):
        t=q[:,1]/L
        return q[:,0]/sh.profile(t)[0],t
    def forward(q,field):
        x,t=chart(q)
        return q+interp_triangle(x,t,xs,ts,field)
    def reference(q):
        x,t=chart(q)
        e=-np.log(x)/np.exp(f.co[0]);core=e>.5
        ref=np.zeros_like(q)
        ref[core,0]=x[core]*f.ref.profile(t[core])[0]
        ref[core,1]=t[core]*f.ref.length
        labels=interp_triangle(e[~core],t[~core],np.linspace(-.5,.5,f.labels.shape[1]),
                               np.linspace(0.,1.,f.labels.shape[0]),f.labels)
        ref[~core,0]=sh.neck+labels[:,0]
        ref[~core,1]=1.-labels[:,1]
        return ref
    def gradient(fn,q,eps):
        return np.stack([(fn(q+eps*v)-fn(q-eps*v))/(2*eps) for v in np.eye(2)],axis=-1)
    rng=np.random.default_rng(3092026)
    t=rng.uniform(.004,.996,3200);x=rng.uniform(.004,.997,len(t))
    q=np.column_stack((x*sh.profile(t)[0],L*t))
    ref=reference(q);p=forward(q,u);native=rel(p[:,0],p[:,1],0)
    for eps in (4e-7,2e-7):
        dr=gradient(reference,q,eps);du=gradient(lambda q:forward(q,u),q,eps)
        F=du@np.linalg.inv(dr)
        h=p[:,0]/ref[:,0];jac=np.abs(np.linalg.det(F))*h*h
        w=(.5*(np.sum(F*F,axis=(1,2))+2*h*h-4)-np.log(jac)+.25*(jac-1)**2)/jac
        err=np.abs(native-w)
        print(json.dumps(dict(kind='independent_pointwise_W_over_J',eps=eps,count=len(q),max_absolute=err.max(),
                              max_relative_to_one=np.max(err/np.maximum(1,abs(native))),
                              min_displacement_J=np.min(np.linalg.det(du)*(p[:,0]/q[:,0])**2),
                              native_range=[native.min(),native.max()])),flush=True)
        if np.max(err/np.maximum(1,abs(native)))>1e-5:
            raise ArithmeticError('Deformed material source does not agree at mapped points')
    owner=Bonded(args.state,(len(xs)-1)//3,len(ts)-1,args.flow)
    X=owner.free.copy();X[:]=u
    coeff=X[owner.mask]*owner.scale
    G=f.joined(sh,False)
    wt,r,a,wcoord,_,_,_,_,live,_,_=G.D.T
    dc=G.V@solution[2]
    mask=(wcoord<L)&(live>0)
    weights=wt[mask]*dc[mask]
    def read(alpha):
        rel.u=alpha*u
        on_grid=weights@rel(r[mask],wcoord[mask],a[mask])
        template=owner.front
        combined=owner.energy(alpha*coeff)[0]
        owner.front=None
        original=owner.energy(alpha*coeff)[0]
        owner.front=template
        # Actual angularly averaged native-motion energy change from adding delta chi.
        pullback=combined-original
        return np.array([combined, original, pullback, on_grid])
    for h in (4e-4,2e-4):
        rows=np.array([read(1+h),read(1-h)])
        derivative=(rows[0]-rows[1])/(2*h)
        print(json.dumps(dict(kind='reciprocal_work',step=h,
                 coordinates='scale all solved displacements, hold outer geometry and joined order fixed',
                 d_material_with_front=derivative[0],d_material_without_front=derivative[1],
                 d_delta_order_material_quadrature=derivative[2],d_delta_order_joined_quadrature=derivative[3],
                 source_work_defect=derivative[3]-derivative[2],native_energies=read(1).tolist())),flush=True)
    print(json.dumps(dict(kind='provenance',source_sha=sha(args.state),flow_sha=sha(args.flow),motion_sha=sha(args.motion))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('state','motion','flow'):p.add_argument(k,type=Path)
    audit(p.parse_args())
