"""Matched crystal force and conversion-field load on the R4 terminal layer.

The material and order assembly share reference sites, volume, stored energy
and physical P1/Q1 interpolants. Azimuthal symmetry applies to this coupling
because neither angular displacement has support in the mobile terminal layer.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_bonded_front_adjoint import angular_weights
from mouth_internal_collar_check import read
from mouth_last_ring_solid import interpolation, source, tri
from mouth_poly_bonded_core import core_stretches
from mouth_supercore_joined_front import check_guard, file_sha
from mouth_supercore_shape import preparation


class TerminalQuadrature:
    def __init__(self, state, checkpoint, ax, z, subdivisions=32, halfwidth=.12):
        self.sh,self.d,_,_=source(state,checkpoint)
        self.chart,self.labels,self.inner,_=read(state)
        self.ax,self.z=np.asarray(ax),np.asarray(z)
        self.L=self.sh.length
        self.p=np.exp(self.chart.co[0])
        self.cut=np.exp(-self.p/2)
        self.ref,_,self.registration=preparation(4)
        nt,nx=self.d.shape[:2]
        n=(nx-1)//3
        self.xt=np.r_[np.linspace(0,self.cut,2*n+1)[:-1],np.linspace(self.cut,1,n+1)]
        self.tt=np.linspace(0,1,nt)
        self.subdivisions=subdivisions
        self.active_start=self.L-halfwidth
        if abs(self.z[-1]-self.L)>1e-12 or np.min(abs(self.z-self.active_start))>1e-12:
            raise ValueError('Terminal layer not present in duct grid')
        # Actual z is affine over each original material triangle. Thus all
        # material belonging to a mobile basis lies within these whole rows.
        zd=self.tt[:,None]*self.L+self.d[:,:,1]
        hi=np.maximum(zd[:-1].max(axis=1),zd[1:].max(axis=1))
        self.first=int(np.flatnonzero(hi>=self.active_start)[0])
        lo=zd[self.first:].min()
        if lo<.62*self.L or np.any(hi[:self.first]>=self.active_start):
            raise ValueError('Overlapping polar map or missed mobile material')
        self.tstart=self.tt[self.first]
        self.radial_knots=np.unique(np.r_[self.ax,self.inner[1]])
        self.depth_knots=np.unique(np.r_[self.z,self.inner[0]*self.L])
        self.g,self.wg=leggauss(5)
        self.g,self.wg=(self.g+1)/2,self.wg/2

    def _radial_lines(self,t):
        sh,xt=self.sh,self.xt
        j=np.clip(self.tt.searchsorted(t)-1,0,len(self.tt)-2)
        v=(t-self.tt[j])/(self.tt[j+1]-self.tt[j])
        rays=np.sort(np.r_[xt,xt[:-1]+np.diff(xt)*(1-v)])
        j0=np.clip(self.chart.t.searchsorted(t)-1,0,len(self.chart.t)-2)
        f=(t-self.chart.t[j0])/(self.chart.t[j0+1]-self.chart.t[j0])
        a=self.chart.a
        cracks=np.exp(-self.p*np.r_[a,a[:-1]+np.diff(a)*(1-f)])
        ids,bary,_,_=tri(xt,self.tt,len(xt),len(self.tt),rays,np.full(rays.size,t))
        disp=np.einsum('ij,ijk->ik',bary,self.d.reshape(-1,2)[ids])
        r=rays*sh.radius(t*self.L)[0]+disp[:,0]
        z=t*self.L+disp[:,1]
        xp=r/sh.radius(z)[0]
        # A positive derivative throughout each native line is necessary to
        # locate a radial Q1 knot uniquely. At this endpoint the map is monotone.
        dr,dz=np.diff(r),np.diff(z)
        for u in (0.,.25,.5,.75,1.):
            radius,slope=sh.radius(z[:-1]+u*dz)
            if np.any(dr*radius-(r[:-1]+u*dr)*slope*dz<=0.):
                raise ValueError('Radial lines fold under the saved displacement')
        knots=self.radial_knots[1:-1]
        k=np.clip(np.searchsorted(xp,knots)-1,0,len(xp)-2)
        u=(knots-xp[k])/(xp[k+1]-xp[k])
        for _ in range(4):
            radius,slope=sh.radius(z[k]+u*dz[k])
            u-=(r[k]+u*dr[k]-knots*radius)/(dr[k]-knots*slope*dz[k])
        if (np.any(u< -1e-10) or np.any(u>1+1e-10)
                or np.max(abs((r[k]+u*dr[k])/sh.radius(z[k]+u*dz[k])[0]-knots))>3e-13):
            raise ArithmeticError('Order radial-knot inverse failed')
        xcross=rays[k]+u*np.diff(rays)[k]
        hit=(self.depth_knots[None,:]>np.minimum(z[:-1],z[1:])[:,None])&(
              self.depth_knots[None,:]<np.maximum(z[:-1],z[1:])[:,None])
        k,key=np.nonzero(hit)
        zcross=rays[k]+np.diff(rays)[k]*(self.depth_knots[key]-z[k])/dz[k]
        return np.unique(np.r_[rays,cracks[(cracks>self.cut)&(cracks<1.)],xcross,zcross])

    def _material(self,x,t,wt):
        sh=self.sh
        b,slope=sh.radius(t*self.L)
        r0=x*b
        xt,tt=self.xt,self.tt
        ids,bary,dx,dt=tri(xt,tt,len(xt),len(tt),x,t)
        gx=dx[:,None]*np.array([[-1.,1.,0.]])/b[:,None]
        gt=dt[:,None]*np.array([[-1.,0.,1.]])/self.L-(x*slope)[:,None]*gx
        D=self.d.reshape(-1,2)[ids]
        disp=np.einsum('qi,qij->qj',bary,D)
        B=np.eye(2)+np.stack([np.einsum('qi,qij->qj',a,D) for a in (gx,gt)],axis=2)
        r,z=r0+disp[:,0],t*self.L+disp[:,1]
        F=np.zeros((len(t),2,2))
        refhoop=np.zeros_like(r)
        core=x<self.cut
        if core.any():
            q,beta,kk,_=core_stretches(sh,self.chart.co,self.ref,self.registration,t[core],r0[core])
            F[core,0,0]=refhoop[core]=q
            F[core,0,1]=beta
            F[core,1,1]=kk
        coll=~core
        if coll.any():
            e=-np.log(x[coll])/self.p
            ia,ba,g1,g2=tri(self.chart.a,self.chart.t,len(self.chart.a),len(self.chart.t),e,t[coll])
            ell=self.labels.reshape(-1,2)[ia]
            ref_r=sh.neck+np.einsum('qi,qi->q',ba,ell[:,:,0])
            g=np.stack([(ell[:,1]-ell[:,0])*g1[:,None],(ell[:,2]-ell[:,0])*g2[:,None]],axis=2)
            s=np.zeros_like(g)
            s[:,0,0]=-self.p*r0[coll]
            s[:,0,1]=x[coll]*self.L*slope[coll]
            s[:,1,1]=self.L
            F[coll]=s@np.linalg.inv(g)
            refhoop[coll]=r0[coll]/ref_r
        J0=abs(np.linalg.det(F))*refhoop**2
        BF=B@F
        K=BF@BF.transpose(0,2,1)
        hoop=refhoop*r/r0
        J=J0*np.linalg.det(B)*(r/r0)**2
        measure=4*np.pi*r0**2*b*self.L*wt/J0
        if np.min(J)<=0 or np.any(z>self.L) or np.any(r>sh.radius(z)[0]):
            raise ValueError('Invalid or out-of-domain material sites')
        W=.5*(K[:,0,0]+K[:,1,1]+2*hoop**2-4)-np.log(J)+.25*(J-1)**2
        return dict(r=r,z=z,x=r/sh.radius(z)[0],K=K,J=J,hoop=hoop,measure=measure,W=W,
                    B=B,ids=ids,bary=bary,gx=gx,gt=gt,partition=wt)

    def blocks(self):
        s=self.subdivisions
        breaks=np.linspace(self.tstart,1,(len(self.tt)-1-self.first)*s+1)
        breaks=np.unique(np.r_[breaks,self.chart.t[(self.chart.t>self.tstart)&(self.chart.t<1)]])
        xs=[];ts=[];ws=[];count=0
        for a,b in zip(breaks[:-1],breaks[1:]):
            for t,tw in zip(a+(b-a)*self.g,(b-a)*self.wg):
                edges=self._radial_lines(t)
                dx=np.diff(edges)
                x=(edges[:-1,None]+dx[:,None]*self.g).ravel()
                xs.append(x);ts.append(np.full(len(x),t));ws.append((tw*dx[:,None]*self.wg).ravel())
                count+=len(x)
                if count>=24000:
                    yield self._material(np.concatenate(xs),np.concatenate(ts),np.concatenate(ws))
                    xs,ts,ws,count=[],[],[],0
        if count:
            yield self._material(np.concatenate(xs),np.concatenate(ts),np.concatenate(ws))


def assemble(model, plane):
    if plane.shape!=(len(model.rr[0]),len(model.z[1])):
        raise ValueError('Wrong carrier grid')
    weight=plane[:,None,:]*angular_weights(model.a)[None,:,None]
    return np.bincount(model.ids[1,0].ravel(),weights=weight.ravel(),minlength=model.size)


def moving_energy(K,J,hoop,radial,axial,step):
    # The supplied gradient is in *current* spatial coordinates.
    jac=J*(1+step*axial)
    if np.min(jac)<=0:raise ValueError('Inverting virtual motion')
    zz=K[:,0,0]*(step*radial)**2+2*K[:,0,1]*(step*radial)*(1+step*axial)+K[:,1,1]*(1+step*axial)**2
    return .5*(K[:,0,0]+zz+2*hoop**2-4)-np.log(jac)+.25*(jac-1)**2


def integrate_motion(quad,delta_order,motion,h=1e-5):
    """Solid contribution with chi convected by motion(r,z)->(u,du/dr,du/dz)."""
    energies=np.zeros(3);analytic=0.
    for v in quad.blocks():
        chi=interpolation(quad.z,quad.ax,delta_order,v['z'],v['x'])[0]
        _,ur,uz=motion(v['r'],v['z'])
        measure=v['measure']*chi
        K,J=v['K'],v['J']
        analytic+=measure@(K[:,0,1]*ur+(K[:,1,1]-1+.5*J*(J-1))*uz)
        for i,s in enumerate((-h,0,h)):
            energies[i]+=measure@moving_energy(K,J,v['hoop'],ur,uz,s)
    return dict(analytic=float(analytic),signed=float((energies[2]-energies[0])/(2*h)),
                energies=energies.tolist(),step=h)


def deposit(plane,ax,z,xs,zs,weight):
    i=np.clip(ax.searchsorted(xs)-1,0,len(ax)-2)
    j=np.clip(z.searchsorted(zs)-1,0,len(z)-2)
    a=(xs-ax[i])/(ax[i+1]-ax[i]);b=(zs-z[j])/(z[j+1]-z[j])
    stride=plane.shape[1]
    for di,wi in ((0,1-a),(1,a)):
        for dj,wj in ((0,1-b),(1,b)):
            plane+=np.bincount((i+di)*stride+j+dj,weight*wi*wj,minlength=plane.size).reshape(plane.shape)


def probe(v,quad,order,direction,h):
    """Signed ring energy change at fixed spatial order; return independent derivative."""
    choice=v['ids']==31*len(quad.xt)+32
    sel=np.any(choice,axis=1)
    if not sel.any():return np.zeros(3)
    take=lambda key:v[key][sel]
    selhat=choice[sel].astype(float)
    hval=np.sum(take('bary')*selhat,axis=1)
    G=np.stack([np.sum(take(der)*selhat,axis=1) for der in ('gx','gt')],axis=1)
    G=np.einsum('qi,qij->qj',G,np.linalg.inv(take('B')))
    r,z=take('r'),take('z');K=take('K');J=take('J');hoop=take('hoop')
    b,slope=quad.sh.radius(z)
    x=r/b
    c,cx,cz=interpolation(quad.z,quad.ax,order,z,x)
    cprime=(cx/b,cz-x*slope*cx/b)[direction]
    hoopprime=hval/r if direction==0 else 0.
    dwork=np.sum(G*K[:,direction,:],axis=1)+G[:,direction]*(-1+.5*J*(J-1))
    dwork+=2*hoopprime*(hoop**2-1+.5*J*(J-1))
    V=take('measure')
    derivative=V@(c*dwork+take('W')*cprime*hval)
    changes=[]
    for signed in (-h,h):
        rr=r+signed*hval*(direction==0);zz=z+signed*hval*(direction==1)
        strain=np.broadcast_to(np.eye(2),K.shape).copy()
        strain[:,direction,:]+=signed*G
        Jnew=J*(1+signed*G[:,direction])*(rr/r)**2
        Knew=strain@K@strain.transpose(0,2,1)
        energy=.5*(Knew[:,0,0]+Knew[:,1,1]+2*(hoop*rr/r)**2-4)-np.log(Jnew)+.25*(Jnew-1)**2
        new=interpolation(quad.z,quad.ax,order,zz,rr/quad.sh.radius(zz)[0])[0]
        changes.append(V@(energy*new-take('W')*c))
    return np.array([derivative,*changes])


def measure(quad,angular_mean):
    plane=np.zeros((len(quad.ax),len(quad.z)))
    p=np.zeros((2,3)); unit=open_work=area=0.;samples=0
    h=1e-6
    ax,z=quad.ax,quad.z
    x,zeta=np.meshgrid(ax,z)
    test=np.where(zeta>quad.active_start,
                  np.sin(2*np.pi*(zeta-quad.active_start)/(quad.L-quad.active_start))*(1-x)*(1+.3*x),0)
    test[:, -1]=0.
    for v in quad.blocks():
        V=v['measure'];W=v['W']
        deposit(plane,ax,z,v['x'],v['z'],V*W)
        values=interpolation(z,ax,test,v['z'],v['x'])[0]
        open_work+=(V*W)@values
        unit+=np.sum(V*W)
        area+=np.sum(v['partition']);samples+=len(W)
        for direction in range(2):
            p[direction]+=probe(v,quad,angular_mean,direction,h)
    if abs(area+quad.tstart-1)>2e-12 or abs(plane.sum()-unit)>2e-10:
        raise ArithmeticError('Gap, overlap, or loss of partition of unity')
    adj=float(abs(np.sum(plane.T*test)-open_work))
    fd=(p[:,2]-p[:,1])/(2*h)
    if adj>2e-10 or np.max(abs(fd-p[:,0]))>2e-4:
        raise ArithmeticError(f'Failed weak-work check, {adj}, {p[:,0]}, {fd}')
    return plane,dict(points=samples,first_material_row=quad.first,logical_area=area,
                      analytic_ring_force_mobile=(-p[:,0]).tolist(),signed_ring_force_mobile=(-fd).tolist(),
                      signed_step=h,adjoint_residual=adj,field_probe_work=open_work,partition_work=unit)


def main(a):
    data=Path('_scratch/s11c/r4_last_rings_fluid01')/a.level/'solid_load.npz'
    check_guard(data)
    with np.load(data) as f:
        ax,z=f['inside_x'],f['duct_z']
        mean=np.einsum('raj,a->jr',f['duct_dc'],angular_weights(f['angles']))
        if str(f['checkpoint_sha'])!=file_sha(a.checkpoint):
            raise ValueError('Stored order from another wall')
    last=None
    for n in a.refinements:
        q=TerminalQuadrature(a.state,a.checkpoint,ax,z,n)
        plane,metrics=measure(q,mean)
        active=np.nonzero(z>q.active_start+1e-12)[0]
        if last is not None:
            delta=(plane-last)[:-1,active]
            metrics.update(free_column_max_abs_difference=float(np.max(abs(delta))),
                           free_l1_difference_relative=float(np.sum(abs(delta))/np.sum(abs(plane[:-1,active]))))
        print(json.dumps(dict(refinements=n,**metrics),allow_nan=False),flush=True)
        last=plane
    np.savez_compressed(a.output,plane=plane,ax=ax,z=z,checkpoint_sha=file_sha(a.checkpoint),
              original_load_field_sha=file_sha(data),source_sha=file_sha(a.state),
              subdivisions=a.refinements[-1],terminal_first_row=q.first)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ('state','checkpoint','output'):
        parser.add_argument(arg,type=Path)
    parser.add_argument('level',choices=['double','local'])
    parser.add_argument('--refinements',type=int,nargs='+',default=[16,32,64])
    main(parser.parse_args())
