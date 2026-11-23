"""Conforming displacement of bonded core and collar inside the ordered part of the radius-four mouth.

The two materials share physical nodes at their original join. A compact depth
window leaves the transmitting front and all carrier boundaries unchanged.
No interfacial spring or new material coefficient is used.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

from mouth_internal_collar_check import read
from mouth_poly_bonded_core import core_stretches
from mouth_supercore_shape import preparation


def emit(**row):
    print(json.dumps(row,allow_nan=False),flush=True)


class Bonded:
    def __init__(self, path, n, k, joined=None, *, shape=None, full_depth=False):
        self.chart,labels,inner,_=read(path)
        if shape is not None:
            self.chart=type(self.chart)(shape,self.chart.co,self.chart.na,self.chart.nt)
        p=self.chart
        if len(p.co)!=1 or p.shape.neck!=4:
            raise ValueError('This readout uses the fixed radius-4 reference spacing')
        pressure=np.exp(p.co[0])
        cut=np.exp(-.5*pressure)
        # The central axis has regular radial motion. The side and both cuts
        # are held for this internal test. The bond has no constrained DOFs.
        self.x=np.r_[np.linspace(0,cut,2*n+1)[:-1],np.linspace(cut,1,n+1)]
        self.t=np.linspace(0.,1.,k+1) if joined is not None or full_depth else np.linspace(.1,.65,k+1)
        x,t=np.meshgrid(self.x,self.t)
        self.mask=np.ones((*x.shape,2),bool)
        self.mask[[0,-1],:,:]=False
        self.mask[:,-1,:]=False
        self.mask[:,0,0]=False
        self.cut_node=2*n
        nx=len(self.x)
        nodes=(np.arange(k)[:,None]*nx+np.arange(nx-1)).ravel()
        corners=np.stack((nodes,nodes+1,nodes+nx,nodes+nx+1),axis=-1)
        triangles=np.stack((corners[:,[0,1,2]],corners[:,[3,2,1]]),axis=1).reshape(-1,3)
        q,w=leggauss(4)
        a,b=np.meshgrid((q+1)/2,(q+1)/2,indexing='ij')
        weights=(np.outer(w,w)/4*(1-a)).ravel()
        u=a.ravel();v=(b*(1-a)).ravel()
        self.basis=np.tile(np.stack((1-u-v,u,v),axis=-1),(len(triangles),1))
        self.id=np.repeat(triangles,len(u),axis=0)
        self.npts=t.size
        xv,tv=x.ravel()[self.id],t.ravel()[self.id]
        xx=np.sum(xv*self.basis,axis=1)
        tt=np.sum(tv*self.basis,axis=1)
        self.xx,self.tt=xx,tt
        dx=xv[:,1]-xv[:,0];dt=tv[:,2]-tv[:,0]
        gx=np.array([-1.,1.,0.])[None,:]/dx[:,None]
        gt=np.array([-1.,0.,1.])[None,:]/dt[:,None]
        sh=p.shape
        b,slope,_=sh.profile(tt)
        self.r=b*xx;self.w=sh.length*tt
        self.dr=gx/b[:,None]
        self.dw=gt/sh.length-(xx*slope)[:,None]*self.dr
        core=xx<cut
        ref,_,reg=preparation(4)
        q,beta,kk,j=core_stretches(sh,p.co,ref,reg,tt[core],self.r[core])
        F=np.zeros((len(tt),2,2))
        J=np.zeros_like(tt)
        hoop=J.copy()
        F[core,0,0]=q;F[core,0,1]=beta;F[core,1,1]=kk
        J[core]=j;hoop[core]=q
        incol=~core
        r=self.r[incol];tin=tt[incol];e=-np.log(xx[incol])/pressure
        fx=(e+.5)*2*p.na;ft=tin*p.nt
        i=np.clip(fx.astype(int),0,2*p.na-1);j=np.clip(ft.astype(int),0,p.nt-1)
        du,dv=fx-i,ft-j
        upper=du+dv>1
        vals=labels.reshape(-1,2)[p.triangles[2*(j*2*p.na+i)+upper]]
        du=np.where(upper,1-du,du);dv=np.where(upper,1-dv,dv)
        S=sh.neck+np.sum(np.column_stack((1-du-dv,du,dv))*vals[:,:,0],axis=1)
        sign=np.where(upper,-1.,1.)[:,None]
        inv=np.stack(((vals[:,1]-vals[:,0])*(2*p.na)*sign,(vals[:,2]-vals[:,0])*p.nt*sign),axis=2)
        spatial=np.zeros_like(inv)
        spatial[:,0,0]=-pressure*r
        spatial[:,0,1]=r*sh.length*slope[incol]/b[incol]
        spatial[:,1,1]=sh.length
        F[incol]=spatial@np.linalg.inv(inv)
        hoop[incol]=r/S
        J[incol]=abs(np.linalg.det(F[incol]))*hoop[incol]**2
        self.K=F@F.transpose((0,2,1))
        self.j=J; self.h=hoop
        self.measure=4*np.pi*self.r**2*b*sh.length*abs(dx*dt)*np.tile(weights,len(triangles))/J
        self.grid=inner[-1];self.L=sh.length
        self.front=None
        if joined is not None:
            from mouth_joined_front_audit import Front
            from mouth_supercore_front import check_guard
            check_guard(joined)
            rows=map(json.loads,(joined.parent/'stdout').read_text().splitlines())
            props=next(row for row in rows if row.get('kind')=='preparation')
            if props['source_sha256']!=hashlib.sha256(path.read_bytes()).hexdigest() or props['winding']!=1:
                raise ValueError('Material state and resolved front must match')
            geom=Front(path,tuple(props['grid']),1).joined(sh,False)
            with np.load(joined) as archive:
                dc=archive['fields'][2,geom.ids[1,0]]
            # Meridional strains multiply chi linearly. Integrate the actual
            # angle-dependent front correction against the 3D solid angle.
            a=geom.a[:-1]; btheta=geom.a[1:]
            right=(np.sin(btheta)-np.sin(a))/(btheta-a)-np.cos(btheta)
            left=np.cos(a)-np.cos(btheta)-right
            mean=np.einsum('raj,a->jr',dc[:,:-1,:],left)+np.einsum('raj,a->jr',dc[:,1:,:],right)
            self.front=(geom.z[1]/sh.length,geom.rr[0],mean)

        self.free=np.zeros(self.mask.shape)
        diag=np.bincount(self.id.ravel(),weights=(self.measure[:,None]*(self.dr**2+self.dw**2+
                                      self.basis**2)).ravel(),minlength=self.npts)
        self.scale=np.broadcast_to(np.sqrt(diag).reshape((*x.shape,1)),self.mask.shape)[self.mask]

    def energy(self,q,diagnose=False):
        u=self.free.copy();u[self.mask]=q/self.scale
        nodes=u.reshape(-1,2)[self.id]
        val=np.einsum('ij,ijk->ik',self.basis,nodes)
        radial=np.einsum('ij,ijk->ik',self.dr,nodes)
        axial=np.einsum('ij,ijk->ik',self.dw,nodes)
        B=np.stack((radial,axial),axis=2)+np.eye(2)
        det=np.linalg.det(B)
        r=self.r+val[:,0];z=self.w+val[:,1]
        hrat=r/self.r;j=self.j*det*hrat**2
        b,slope,_=self.chart.shape.profile(z/self.L)
        x=r/b
        if np.any((j<=0)|(hrat<=0)|(x<=0)|(x>=1)|(z<=0)|(z>=self.L)):
            raise ValueError('Nonpositive trial map or outside crystalline order grid')
        Knew=B@self.K@B.transpose(0,2,1)
        hnew=self.h*hrat
        W=.5*(np.trace(Knew,axis1=1,axis2=2)+2*hnew**2-4)-np.log(j)+.25*(j-1)**2
        nt,nx=np.subtract(self.grid.shape,1)
        t=z/self.L*nt;xx=x*nx
        k=t.astype(int);i=xx.astype(int)
        a,bb=xx-i,t-k
        c=self.grid
        low=c[k,i]; high=c[k+1,i]
        da=c[k,i+1]-low; db=c[k+1,i+1]-high
        chi=(1-bb)*(low+a*da)+bb*(high+a*db)
        cr=nx*((1-bb)*da+bb*db)/b
        cw=nt*(high-low+a*(db-da))/self.L-x*slope*cr
        if self.front is not None:
            tz,rx,dc=self.front
            zv=z/self.L
            kz=np.clip(np.searchsorted(tz,zv)-1,0,len(tz)-2)
            ix=np.clip(np.searchsorted(rx,x)-1,0,len(rx)-2)
            dt=tz[kz+1]-tz[kz];dx=rx[ix+1]-rx[ix]
            h=(zv-tz[kz])/dt;ux=(x-rx[ix])/dx
            lo=dc[kz,ix];hi=dc[kz+1,ix]
            da=dc[kz,ix+1]-lo;db=dc[kz+1,ix+1]-hi
            chi+=(1-h)*(lo+ux*da)+h*(hi+ux*db)
            dcdr=((1-h)*da+h*db)/(dx*b)
            cr+=dcdr
            cw+=(hi-lo+ux*(db-da))/(dt*self.L)-x*slope*dcdr
        value=(self.measure*chi)@W
        if diagnose:
            eig=np.linalg.eigvalsh(Knew)
            return dict(energy=value,min_J=j.min(),max_J=j.max(),min_stretch=min(hnew.min(),np.sqrt(eig[:,0]).min()),
                        max_stretch=max(hnew.max(),np.sqrt(eig[:,1]).max()),minimum_order=chi.min(),
                        peak_join_radial=abs(u[:,self.cut_node,0]).max(),peak_join_axial=abs(u[:,self.cut_node,1]).max(),
                        middle_join_motion=u[len(u)//2,self.cut_node].tolist())
        pressure=-1+.5*j*(j-1)
        P=B@self.K+pressure[:,None,None]*np.linalg.inv(B).transpose(0,2,1)
        force=np.column_stack((W*cr+2*chi*(hnew**2+pressure)/r,W*cw))
        g=np.zeros((self.npts,2))
        for d in (0,1):
            load=self.measure[:,None]*(force[:,d,None]*self.basis+chi[:,None]*(
                    P[:,d,0,None]*self.dr+P[:,d,1,None]*self.dw))
            g[:,d]=np.bincount(self.id.ravel(),weights=load.ravel(),minlength=self.npts)
        self.native=g.reshape(self.mask.shape)
        return value,self.native[self.mask]/self.scale


def run(a):
    b=Bonded(a.state,a.n,a.k,a.joined)
    y=np.zeros(b.mask.sum())
    def objective(y):
        try:return b.energy(y)
        except ValueError:return 1.e10+y@y,2*y
    e,g=objective(y)
    rng=np.random.default_rng(4072)
    d=rng.normal(size=len(y));d/=np.linalg.norm(d)
    read=b.energy(y,True)
    for h in (1.e-4,5.e-5):
        fd=(objective(y+h*d)[0]-objective(y-h*d)[0])/(2*h)
        emit(kind='directional_check',step=h,analytic=g@d,signed=fd)
        if abs(fd-g@d)>3e-6*max(1,abs(fd)):
            raise ValueError('Material gradient check failed')
    emit(kind='initial',sha256=hashlib.sha256(a.state.read_bytes()).hexdigest(),depth_span=[b.t[0],b.t[-1]],
         resolved_order_sha256=hashlib.sha256(a.joined.read_bytes()).hexdigest() if a.joined else None,
         radial_cells=3*a.n,axial_cells=a.k,**read)
    opt=minimize(objective,y,jac=True,method='L-BFGS-B',options=dict(maxiter=2500,gtol=1e-9,ftol=1e-15,maxls=60))
    en,gr=objective(opt.x)
    save_info={}
    if a.save is not None:
        displacement=np.zeros(b.mask.shape)
        displacement[b.mask]=opt.x/b.scale
        np.savez_compressed(a.save, displacement=displacement, radial_fraction=b.x, depth_fraction=b.t,
                            radius=4, source_sha256=hashlib.sha256(a.state.read_bytes()).hexdigest(),
                            joined_sha256=hashlib.sha256(a.joined.read_bytes()).hexdigest() if a.joined else '')
        save_info['displacement_sha256']=hashlib.sha256(a.save.read_bytes()).hexdigest()
    emit(kind='bonded_release',**save_info,**b.energy(opt.x,True),energy_change=en-e,native_free_residual=abs(b.native[b.mask]).max(),
         scaled_residual=abs(gr).max(),iterations=opt.nit,optimizer=str(opt.message),
         scope='Both bonded pieces mobile inside the indicated window, all its boundary coordinates held except radial symmetry axis.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('--n',type=int,required=True);p.add_argument('--k',type=int,required=True)
    p.add_argument('--joined',type=Path,help='Use whole depth with the matched angularly resolved crystalline front held at its measured location')
    p.add_argument('--save',type=Path,help='Save the physical displacement sampled at each shared and internal material node')
    run(p.parse_args())
