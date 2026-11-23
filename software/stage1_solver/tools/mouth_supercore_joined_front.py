"""Global variational transmission of a linked vortex through a finite terminal layer.

Finite elements for both coherent fields extend from the feeding reservoir to
the receiving reservoir; there is one shared cap trace and no sharp contact term.
Only the order field uses finite-band matching to a supplied mechanical preparation.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import factorized
from scipy.special import expit

from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA, WELL, WIDTH
from mouth_exchange_carrier import BIAS, mapped_nodes
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_front import frozen_extra, check_guard, file_sha
from mouth_winding_resolved_patch import Patch
from mouth_wound_carrier import straight_vortex


def emit(**items):
    print(json.dumps(items, allow_nan=False), flush=True)


def state_from(path):
    check_guard(path)
    with np.load(path,allow_pickle=False) as f:
        if int(f['radius']) != 4:
            raise ValueError('Only the radius-4 preparation is released')
        return Meridian(4.,f['shape_params'],extra=frozen_extra(f)), f['order'].copy()


class Joined(Patch):
    """Resolved phase and density throughout; free crystalline front on finite band.

    elastic_density(r,w,theta) is the registered elastic energy per unit current
    volume for the inner material. reference_order(r/b,w/L) lives on the supplied
    meridional chart; the receiving cap continues the neutral outer order.
    """
    def __init__(self,shape,order,winding,grid=(9,7,7),halfwidth=.12,elastic_density=None,precondition=True,terminal_profile=None):
        self.shape,self.winding=shape,winding
        self.h=halfwidth
        self.precondition=precondition
        self.terminal_profile=terminal_profile
        nr,na,nz=grid
        if shape.neck!=4 or winding not in (0,1) or not 0<halfwidth<shape.length:
            raise ValueError('Invalid geometry or winding')
        self.rr,self.z=mapped_nodes(nr,shape.length,28.,22.)
        # Match both finite-band planes to grid nodes. Same topology on signed shapes.
        lower=np.linspace(0.,shape.length-halfwidth,max(3,nr)+1)
        middle=np.linspace(shape.length-halfwidth,shape.length,nz+1)
        above=np.linspace(shape.length,shape.length+halfwidth,nz+1)
        far=shape.length+halfwidth+(22.-shape.length-halfwidth)*(np.sinh(2*np.linspace(0,1,nr+1))/np.sinh(2))**2
        self.z=(self.z[0],np.r_[lower[:-1],middle],np.r_[above[:-1],far])
        self.a=np.pi/2*np.linspace(0.,1.,na+1)**1.35
        self.ids={}
        keys={}
        coords=[]
        zbase=0
        for slab,zs in enumerate(self.z):
            for side,xs in enumerate(self.rr):
                box=np.empty((len(xs),len(self.a),len(zs)),int)
                for i,k,j in np.ndindex(box.shape):
                    iz=zbase+j; ir=side*nr+i
                    # Below the terminal plane the exterior is sealed out of the duct.
                    # At the receiving rim all traces share the finite-element node.
                    feedbank=slab==0 or slab==1 and side==0
                    if slab==1 and j==len(zs)-1:
                        feedbank=False
                    key=(iz,ir,0 if ir==0 else k,feedbank)
                    if key not in keys:
                        keys[key]=len(coords);coords.append((slab,side,xs[i],self.a[k],zs[j]))
                    box[i,k,j]=keys[key]
                self.ids[slab,side]=box
            zbase+=len(zs)-1
        self.size=len(coords)
        self.feed=np.unique(self.ids[0,1][-1])
        self.bulk=np.unique(np.r_[self.ids[1,1][-1].ravel(),self.ids[2,1][-1].ravel(),
                                  self.ids[2,1][:,:,-1].ravel(),self.ids[2,0][:,:,-1].ravel()])
        self.seed=np.array([np.ones(self.size),np.zeros(self.size),np.zeros(self.size)])
        self.seed[1,self.bulk]=BIAS
        fmask=np.ones(self.size,bool);fmask[np.r_[self.feed,self.bulk]]=False
        omask=np.zeros(self.size,bool)
        omask[self.ids[1,0][:-1,:,-nz:-1]]=True
        omask[self.ids[2,0][:-1,:,:nz]]=True
        self.masks=np.array([fmask,fmask,omask]);self.slices=np.r_[0,np.cumsum(self.masks.sum(axis=1))]
        self.matrix(order,elastic_density)

    def matrix(self,order,elastic_density):
        V=[];G=[[],[],[]];out=[];refs=[];vortex=straight_vortex(self.winding,28.)
        interp=RegularGridInterpolator((np.linspace(0,1,len(order)),np.linspace(0,1,order.shape[1])),order)
        bits=np.array([(i,k,j) for i in (0,1) for k in (0,1) for j in (0,1)])
        quadr=(1+np.array([-1,1])/np.sqrt(3))/2
        self.blocks={}
        for (slab,side),ids in self.ids.items():
            xs=self.rr[side];zs=self.z[slab];a=self.a
            box=ids.shape; nr,na,nz=np.subtract(box,1)
            nodes=np.stack([ids[i:i+nr,k:k+na,j:j+nz].ravel() for i,k,j in bits],axis=1)
            rows=np.broadcast_to(np.arange(len(nodes))[:,None],nodes.shape).ravel()
            def mat(vals):
                return sparse.coo_matrix((np.broadcast_to(vals,nodes.shape).ravel(),(rows,nodes.ravel())),shape=(len(nodes),self.size)).tocsr()
            x0,a0,z0=np.meshgrid(xs[:-1],a[:-1],zs[:-1],indexing='ij')
            dx,da,dz=[v.ravel() for v in np.meshgrid(np.diff(xs),np.diff(a),np.diff(zs),indexing='ij')]
            xn,an,zn=np.meshgrid(xs,a,zs,indexing='ij')
            if (slab,side)==(1,0):
                ref=interp((zn/self.shape.length,xn))
            elif (slab,side)==(2,0):
                ref=(expit(-(zn-self.shape.length)/WIDTH) if self.terminal_profile is None
                     else self.terminal_profile(xn,zn-self.shape.length))
            elif side==1 and slab>0 and self.terminal_profile is not None:
                b,_=self.shape.radius(zn)
                ref=self.terminal_profile.outside(b+xn*(28.-b),zn,slab==2)
            else:
                ref=np.full(box,1. if slab==0 else 0.)
            cs=np.stack([ref[i:i+nr,k:k+na,j:j+nz].ravel() for i,k,j in bits],axis=1)
            self.blocks[slab,side]=(nodes,len(out))
            for xx in quadr:
                for aa in quadr:
                    for zz in quadr:
                        f=np.where(bits,[xx,aa,zz],1-np.array([xx,aa,zz]))
                        v=f.prod(axis=1); raw=np.stack([(2*bits[:,i]-1)*f[:,(i+1)%3]*f[:,(i+2)%3] for i in range(3)])
                        x=x0.ravel()+xx*dx;angle=a0.ravel()+aa*da;z=z0.ravel()+zz*dz
                        b,slope=self.shape.radius(z)
                        width=28.-b if side else b
                        r=b+x*width if side else x*width
                        ds=raw[0]/(dx*width)[:,None]
                        ga=raw[1]/(r*da)[:,None]
                        gz=raw[2]/dz[:,None]-(slope*(1-x if side else x))[:,None]*ds
                        d=[ds,ga,gz]
                        wt=4*np.pi*r*r*np.sin(angle)*width*dx*da*dz/8
                        s=r*np.sin(angle);F,Fp=vortex(s)
                        y=z-self.shape.length
                        cut=((slab==2 and side==0)&(z0.ravel()==self.shape.length)).astype(int)
                        cut-=((slab==1 and side==0)&(z0.ravel()+dz>self.shape.length-1.e-12)).astype(int)
                        live=(side==0 and slab in (1,2))&(abs(y)<self.h)
                        if elastic_density is not None and slab==1 and side==0:
                            Eel=elastic_density(r,z,angle)
                        else:Eel=np.zeros_like(r)
                        V.append(mat(v))
                        for i,gg in enumerate(d):G[i].append(mat(gg))
                        refs.append(np.column_stack((cs@v,*[np.sum(cs*g,axis=1) for g in d])))
                        out.append(np.column_stack((wt,r,angle,z,F,Fp,s,cut,live,wt/dz,Eel)))
        self.V=sparse.vstack(V).tocsr(); self.G=[sparse.vstack(g).tocsr() for g in G]
        self.D=np.vstack(out); self.ref=np.vstack(refs); self.hq=self.D[:,0]
        self.B=self.boundary(vortex,quadr)
        if not self.precondition:return
        weighted=self.hq*self.D[:,4]**2
        def metric(weights,mass):
            M=self.V.T@self.V.multiply(mass[:,None])
            for g in self.G:M+=g.T@g.multiply(weights[:,None])
            return M
        metric_v=metric(weighted,self.hq*self.D[:,4]**10*4.)
        metric_p=metric(weighted,0.*weighted)
        metric_c=metric(self.hq*KAPPA,2*WELL*self.hq)
        self.solvers=[factorized(M[m][:,m].tocsc()) for M,m in zip((metric_v,metric_p,metric_c),self.masks)]

    def boundary(self,vortex,quad):
        # Exact integration-by-parts boundary term of the straight-vortex factorization;
        # the cap traces now cancel pointwise, leaving only the two sealed lateral traces.
        B=sparse.csr_matrix((self.size,self.size))
        zs=self.z[1]; dzz=np.diff(zs)[None,:]
        angle=self.a[:-1,None];dt=np.diff(self.a)[:,None]
        for side in (0,1):
            face=self.ids[1,side][0 if side else -1]
            ids=np.stack([face[k:k+len(self.a)-1,j:j+len(zs)-1] for k in (0,1) for j in (0,1)],-1).reshape(-1,4)
            for c in quad:
                for b in quad:
                    r,_=self.shape.radius(zs[None,:-1]+dzz*b)
                    sine=np.sin(angle+c*dt);v,df=vortex(r*sine)
                    weight=((1-2*side)*4*np.pi*r*r*sine**2*dzz*dt*v*df/4).ravel()
                    vec=np.kron([1-c,c],[1-b,b])
                    arr=weight[:,None,None]*np.outer(vec,vec)
                    B+=sparse.coo_matrix((arr.ravel(),(np.broadcast_to(ids[:,:,None],arr.shape).ravel(),
                                          np.broadcast_to(ids[:,None,:],arr.shape).ravel())),shape=B.shape).tocsr()
        return B

    def evaluate(self,u,read=False):
        wt,r,a,z,F,Fp,s,cut,live,area,elastic=self.D.T
        v,p,dc=u; t=self.V@v; n=(F*t)**2
        chi=self.ref[:,0]+self.V@dc
        gv,gp,gc=(np.column_stack([G@q for G in self.G]) for q in u)
        gc+=self.ref[:,1:]
        square=lambda q:np.einsum('ij,ij->i',q,q)
        den=n+ORDER_INERTIA*square(gc); dot=np.einsum('ij,ij->i',gp,gc)
        pi=ORDER_INERTIA*n*dot/den
        vel=gp-(ORDER_INERTIA*dot/den)[:,None]*gc
        curr=vel*n[:,None]
        c0=self.ref[:,0]
        order=.5*KAPPA*(square(gc)-square(self.ref[:,1:]))+WELL*((chi*(1-chi))**2-(c0*(1-c0))**2)+elastic*(chi-c0)
        # Subtract the neutral vortex and the reference phase-field energy once.
        y=t*t
        interaction=F**10*(y-1)**2*(y**3+2*y*y+3*y+4)/20
        E=wt@(.5*F**2*square(gv)+interaction+.5*np.einsum('ij,ij->i',gp,curr)+order)+.5*((v-1)@(self.B@(v+1)))
        if read:
            gf=F[:,None]*gv+t[:,None]*Fp[:,None]*np.column_stack((np.sin(a),np.cos(a),np.zeros_like(a)))
            physical=.5*(square(gf)+n*(self.winding/s)**2+KAPPA*square(gc)+np.einsum('ij,ij->i',gp,curr))+(n**5-5*n+4)/20+WELL*(chi*(1-chi))**2
            T=n*vel[:,2]**2+gf[:,2]**2+KAPPA*gc[:,2]**2-physical
            return E,chi,n,curr,vel,T
        gr=np.stack((self.V.T@(wt*(.5*F**10*t*(t**8-1)+F**2*t*(square(gp)-ORDER_INERTIA**2*square(gc)*dot**2/den**2))),
                     np.zeros(self.size),
                     self.V.T@(wt*(2*WELL*chi*(1-chi)*(1-2*chi)+elastic))))
        gr[0]+=self.B@v
        for i,G in enumerate(self.G):
            gr[0]+=G.T@(wt*F**2*gv[:,i])
            gr[1]+=G.T@(wt*curr[:,i])
            gr[2]+=G.T@(wt*(KAPPA*gc[:,i]-pi*vel[:,i]))
        return E,gr

    def measure(self,u):
        E,chi,n,j,vel,T=self.evaluate(u,True)
        _,gr=self.evaluate(u); feed=-gr[1,self.feed].sum();bulk=gr[1,self.bulk].sum()
        wt,r,a,z,F,Fp,s,cut,live,area,el=self.D.T
        use=cut==-1; mtotal=np.sqrt(np.sum(vel**2,axis=1)+(self.winding/s)**2)/n**2
        dense=use&(n>.5)&(mtotal<1)
        jp=np.maximum(j[:,2],0.)
        flux=area[use]@j[use,2]
        rec=cut==1
        mask=np.zeros(self.size)
        mask[self.ids[1,0][:,:,-1]]=1.
        # Conservative integrated disk current, tested against a piecewise-linear
        # unit lift on the last core element. Includes the sloped radial weak flux.
        through=sum((wt[use]*(g@mask)[use])@j[use,i] for i,g in enumerate(self.G))
        delta_a=area[rec]@j[rec,2]
        return dict(excess_energy=E,feed_current=feed,bulk_current=bulk,port_error=abs(feed-bulk),
                    weak_disk_current=through,incident_slice_current=flux,receiving_slice_current=delta_a,positive_slice_current=area[use]@jp[use],
                    dense_total_subsonic_share=area[dense]@jp[dense]/(area[use]@jp[use]) if feed>1e-8 else None,
                    minimum_total_mach=float(mtotal[use].min()),total_axial_stress=area[use]@T[use],
                    crystalline_range_in_band=[float(chi[live>0].min()),float(chi[live>0].max())])

    def run_state(self,rest=False):
        u,info=self.solve(rest)
        reading=self.measure(u)
        emit(kind='joined',winding=self.winding,rest=rest,**reading,**info)
        if reading['port_error']>3e-5:raise ArithmeticError('Far reservoir current mismatch')
        return u,reading


def registered_density(shape,reference,co,co0):
    """W/J for the original registered collar and its bonded central material."""
    from numpy.polynomial.legendre import legval,legder
    from mouth_finite_bore_support import graph_values
    from mouth_poly_bonded_core import core_stretches,_energy
    from mouth_polycrystal_material import stable_energy
    def density(r,w,theta):
        t=w/shape.length
        b=shape.profile(t)[0]
        lnp=legval(2*t-1,co); pt=2*legval(2*t-1,legder(co))
        a=-np.log(r/b)/np.exp(lnp)
        _,rs,rz,h,j,_=graph_values(np.minimum(a,.5),t,shape,lnp,pt)
        W=stable_energy(rs*rs+rz*rz+shape.length**2+2*h*h,j)/j
        q,d,k,jc=core_stretches(shape,co,reference,co0,t,r)
        return np.where(a<=.5,W,_energy(q,d,k)/jc)
    return density


def main(args):
    shape,chi=state_from(args.state)
    if args.delta:
        params=shape.params.copy()
        if args.axis=='depth': params[0]+=np.log1p(args.delta/shape.length)
        elif args.axis=='exit': params[1]+=np.log1p(args.delta/shape.exit)
        else:raise ValueError('Signed change requires axis')
        shape=Meridian(4.,params,extra=shape.extra)
    from mouth_supercore_shape import preparation
    ref,_,co0=preparation(12)
    with np.load(args.state) as data:co=data['registration']
    model=Joined(shape,chi,args.winding,args.grid,args.width,registered_density(shape,ref,co,co0))
    emit(kind='setup',source=file_sha(args.state),depth=shape.length,exit=shape.exit,
         winding=args.winding,grid=args.grid,cells=len(model.hq)//8,half_layer=args.width)
    u,reading=model.run_state()
    if args.save:np.savez_compressed(args.save,fields=u)
    if args.rest:
        u0,rest=model.run_state(True)
        _,_,n,_,vel,t=model.evaluate(u,True)
        t-=model.evaluate(u0,True)[-1]
        wt,r,a,z,F,Fp,s,cut,live,area,elastic=model.D.T
        m=np.sqrt(np.sum(vel**2,axis=1)+(args.winding/s)**2)/n**2
        choose=cut==-1; dense=choose&(n>.5)&(m<1)
        for sign in (-1,1):
            choose=cut==sign;dense=choose&(n>.5)&(m<1)
            emit(kind='drive_minus_rest',sign=sign,energy=reading['excess_energy']-rest['excess_energy'],
                 axial=area[choose]@t[choose],dense=area[dense]@t[dense])


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('state',type=Path)
    ap.add_argument('--grid',type=int,nargs=3,default=(9,7,7))
    ap.add_argument('--winding',type=int,choices=(0,1),default=1)
    ap.add_argument('--width',type=float,default=.12)
    ap.add_argument('--rest',action='store_true')
    ap.add_argument('--save',type=Path)
    ap.add_argument('--axis',choices=('depth','exit'))
    ap.add_argument('--delta',type=float,default=0.)
    main(ap.parse_args())
