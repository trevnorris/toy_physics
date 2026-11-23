"""Return the relaxed core/collar elastic load to the finite front and both reservoirs."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_internal_collar import CollarPair
from mouth_joined_front_audit import Front
from mouth_joint_collar_release import joined_collar_density
from mouth_poly_bonded_core import core_stretches
from mouth_poly_shallow_release import CUT
from mouth_polycrystal_material import stable_energy
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import Joined, emit
from mouth_supercore_terminal_order import terminal


class Relocated:
    def __init__(self, base, field):
        self.base=base
        self.x=field['radial_fraction'];self.t=field['depth_fraction']
        self.u=field['displacement']
        self.sh=base.sh
        self.co=base.co
        self.pair=CollarPair(self.sh,self.co,base.na,base.nt)
        self.original=joined_collar_density(self.sh,base.co,base.ref,self.pair,base.labels)
        if self.t[0]!=0 or self.t[-1]!=1 or len(self.co)!=1:
            raise ValueError('Expected whole-depth, constant registered initial spacing')
        if np.any(self.u[[0,-1]]) or np.any(self.u[:,-1,:]) or np.any(self.u[:,0,0]):
            raise ValueError('Movement touches an outer trace')

    def mapping(self,x,t):
        i=np.clip(np.searchsorted(self.x,x)-1,0,len(self.x)-2)
        j=np.clip(np.searchsorted(self.t,t)-1,0,len(self.t)-2)
        dx=self.x[i+1]-self.x[i];dt=self.t[j+1]-self.t[j]
        a=(x-self.x[i])/dx;v=(t-self.t[j])/dt
        upper=(a+v)>1
        corners=self.u[j+upper,i+upper]
        sign=np.where(upper,-1.,1.)
        du=(self.u[j+upper,i+(~upper)]-corners)*sign[:,None]/dx[:,None]
        dv=(self.u[j+(~upper),i+upper]-corners)*sign[:,None]/dt[:,None]
        motion=(corners+du*((a-upper)*dx)[:,None]+dv*((v-upper)*dt)[:,None])
        b,bw,_=self.sh.profile(t)
        r=b*x;w=self.sh.length*t
        B=np.zeros((len(x),2,2))
        B[:,:,0]=du/b[:,None]
        B[:,:,1]=dv/self.sh.length-B[:,:,0]*(x*bw)[:,None]
        B+=np.eye(2)
        return r,w,motion,B,b,bw

    def initial_strain(self,x,t,r,b,bw):
        press=np.exp(self.co[0])
        e=-np.log(x)/press
        col=e<.5
        F=np.zeros((len(x),2,2));h=np.zeros(len(x));J=h.copy()
        core=~col
        if np.any(core):
            q,s,k,j=core_stretches(self.sh,self.co,self.base.ref,self.co,t[core],r[core])
            F[core,0,0]=q;F[core,0,1]=s;F[core,1,1]=k;h[core]=q;J[core]=j
        if np.any(col):
            p=self.pair
            fx=(e[col]+.5)*2*p.na;ft=t[col]*p.nt
            i=np.clip(fx.astype(int),0,2*p.na-1);j=np.clip(ft.astype(int),0,p.nt-1)
            du,dv=fx-i,ft-j
            upper=du+dv>1
            vals=self.base.labels.reshape(-1,2)[p.triangles[2*(j*2*p.na+i)+upper]]
            du=np.where(upper,1-du,du);dv=np.where(upper,1-dv,dv)
            S=self.sh.neck+np.sum(np.column_stack((1-du-dv,du,dv))*vals[:,:,0],axis=1)
            sign=np.where(upper,-1.,1.)[:,None]
            inv=np.stack(((vals[:,1]-vals[:,0])*(2*p.na)*sign,(vals[:,2]-vals[:,0])*p.nt*sign),axis=2)
            spatial=np.zeros_like(inv)
            spatial[:,0,0]=-press*r[col]
            spatial[:,0,1]=r[col]*self.sh.length*bw[col]/b[col]
            spatial[:,1,1]=self.sh.length
            F[col]=spatial@np.linalg.inv(inv)
            h[col]=r[col]/S;J[col]=abs(np.linalg.det(F[col]))*h[col]**2
        return F,h,J

    def __call__(self,rnew,wnew,theta):
        rnew,wnew,theta=np.broadcast_arrays(rnew,wnew,theta)
        L=self.sh.length
        use=(wnew>0)&(wnew<L)&(rnew<=self.sh.radius(wnew)[0])&(rnew>0)
        ans=self.original(rnew,wnew,theta)
        if not np.any(use):return ans
        rr=rnew[use];ww=wnew[use]
        t=ww/L;x=rr/self.sh.radius(ww)[0]
        for _ in range(24):
            r,w,d,B,b,bw=self.mapping(x,t)
            error=np.column_stack((r-rr,w-ww))+d
            if np.max(abs(error))<2e-12:break
            # B maps current spatial corrections from the old chart to the new one.
            ds=np.linalg.solve(B,error[:,:,None])[:,:,0]
            t-=ds[:,1]/L;x-=(ds[:,0]-x*bw*ds[:,1])/b
            t=np.clip(t,0,1);x=np.clip(x,0,1)
        else:raise ArithmeticError('Inverse material movement unresolved')
        F,h,J=self.initial_strain(x,t,r,b,bw)
        F=B@F;h*=rr/r;J*=np.linalg.det(B)*(rr/r)**2
        if np.min(J)<=0:raise ValueError('Rearranged crystal folded')
        ans[use]=stable_energy(np.sum(F*F,axis=(1,2))+2*h*h,J)/J
        return ans


def run(a):
    from mouth_joined_front_stationarity import solve
    check_guard(a.flow);check_guard(a.motion)
    with np.load(a.motion) as fields:
        sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        if fields['source_sha256'].item()!=sha(a.state):
            raise ValueError('Mismatched material prepared for the displacement')
        flowrows=list(map(json.loads,(a.flow.parent/'stdout').read_text().splitlines()))
        setup=next(q for q in flowrows if q.get('kind')=='preparation')
        if setup['source_sha256']!=sha(a.state):raise ValueError('Mismatched flow geometry')
        base=Front(a.state,tuple(setup['grid']),setup['winding'])
        density=Relocated(base,fields)
        declared_order=fields['joined_sha256'].item()
    # Rebuilding the undeformed law must be independent of the map implementation.
    tt,xx=np.meshgrid(np.linspace(.013,.987,43),np.linspace(.005,.995,67))
    r=xx*base.sh.profile(tt)[0];w=base.sh.length*tt
    motion=density.u.copy();density.u.fill(0)
    mismatch=float(np.max(abs(density(r,w,0)-density.original(r,w,np.zeros_like(r)))))
    density.u=motion
    if mismatch>2e-10:raise ArithmeticError(f'Baseline density mismatch {mismatch}')
    pair=density.pair
    nr=base.order_mesh;t=np.linspace(0,1,2*nr+1);x=np.linspace(0,1,nr+1)
    out=pair.evaluate(base.labels,(t,x,base.c),((1-np.cos(np.pi*t))/2,.1*x+.9*x**2,base.outer))['outer_source']
    ref=terminal(base.sh,nr,CUT,lambda wa,xo:out)
    model=Joined(base.sh,base.c,base.N,base.grid,elastic_density=density,terminal_profile=ref)
    with np.load(a.flow) as f:
        saved=[f[k] for k in ('fields','rest')]
    emit(kind='preparation',grid=base.grid,winding=base.N,source_sha256=sha(a.state),
         flow_sha256=sha(a.flow),displacement_sha256=sha(a.motion),original_order_sha256=declared_order,
         zero_motion_error=mismatch)
    model.warm_path=a.flow
    solutions=[]
    for rest,initial in zip((False,True),saved):
        u,info=solve(model,rest)
        result=model.measure(u);solutions.append(u)
        emit(kind='rematched',rest=rest,**result,stationarity=info['stationarity'],
             max_additional_order_change=np.max(abs(u[2]-initial[2])))
        if result['port_error']>3e-5:raise ArithmeticError('Reservoir mismatch')
    np.savez_compressed(a.save,fields=solutions[0],rest=solutions[1])
    emit(kind='state',field_sha256=sha(a.save))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','motion','flow','save'):p.add_argument(name,type=Path)
    run(p.parse_args())
