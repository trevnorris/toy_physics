"""Relax extra solid nodes in the six-triangle patch with its perimeter fixed.

Refine displacement DOFs independently of cubature while keeping the physical
patch, wall, material registration and saved carrier/order field fixed.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from scipy.optimize._numdiff import group_columns

from mouth_last_ring_solid import Ring, interpolation, order_tables, source, tri
from mouth_supercore_joined_front import file_sha


def rim_patch(state, sh, d, tables, subdivision, extend):
    """Extend two cells of the artificial patch perimeter into the same material."""
    root=Ring(state,sh,d,tables,(31,32),subdivision)
    if not extend:
        return root,[(0,0)]
    offsets=[(0,0),(-1,0),(0,-1),(-1,-1)]
    fields=('rbase','zbase','r','z','B','j0','refhoop','K','w')
    parts={name:[] for name in fields}
    dx,dt=root.vcell
    for k,(di,dj) in enumerate(offsets):
        r=root if k==0 else Ring(state,sh,d,tables,(31+dj,32+di),subdivision)
        t=r.zbase/sh.length
        x=r.rbase/sh.radius(r.zbase)[0]
        A=(x-(1-dx))/dx
        T=(t-31/32)/dt
        fresh=np.ones(A.size,dtype=bool)
        for oi,oj in offsets[:k]:
            # Complete base triangles are shared; quadrature lies strictly inside.
            fresh &= np.maximum.reduce([abs(A-oi),abs(T-oj),abs(A+T-oi-oj)])>1.
        for name in fields:
            parts[name].append(getattr(r,name)[fresh])
    parts={name:np.concatenate(values) for name,values in parts.items()}
    A=(parts['rbase']/sh.radius(parts['zbase'])[0]-(1-dx))/dx
    T=(parts['zbase']/sh.length-31/32)/dt
    phi=np.maximum(0.,1.-np.maximum.reduce([abs(A),abs(T),abs(A+T)]))
    return SimpleNamespace(**parts,sh=sh,order=tables,vcell=root.vcell,phi=phi),offsets


class Patch:
    def __init__(self, ring, divisions, offsets):
        self.ring=ring
        t=ring.zbase/ring.sh.length
        roof,slope=ring.sh.radius(ring.zbase)
        x=ring.rbase/roof
        dx,dt=ring.vcell
        cx=(1-dx)
        ct=31/32
        lo=np.min(offsets,axis=0)-1
        hi=np.max(offsets,axis=0)+1
        xa=np.arange(lo[0]*divisions,hi[0]*divisions+1)/divisions
        ta=np.arange(lo[1]*divisions,hi[1]*divisions+1)/divisions
        A,T=np.meshgrid(xa,ta)
        rho=np.maximum.reduce([abs(A),abs(T),abs(A+T)])
        # A node is free only if all six adjacent triangle centroids lie in the
        # union. This fixes the physical disk and wall in both perimeter tests.
        free=np.ones_like(A,dtype=bool)
        for da,dtc in [(1,1),(-1,-1),(2,-1),(-2,1),(1,-2),(-1,2)]:
            X=A+da/(3*divisions)
            Z=T+dtc/(3*divisions)
            interior=np.zeros_like(free)
            for oi,oj in offsets:
                interior |= np.maximum.reduce([abs(X-oi),abs(Z-oj),abs(X+Z-oi-oj)])<1.
            free &= interior
        free=free.ravel()
        self.hat=np.maximum(1-rho,0.).ravel()[free]
        xaxis=cx+xa*dx
        taxis=ct+ta*dt
        ids,N,gx,gt=tri(xaxis,taxis,len(xa),len(ta),x,t)
        rows=np.repeat(np.arange(len(x)),3)
        def operator(values):
            return sparse.coo_matrix((values.ravel(),(rows,ids.ravel())),
                                     shape=(len(x),free.size)).tocsr()[:,free]
        Nr=(gx/roof)[:,None]*np.array([[-1.,1.,0.]])
        Nz=(gt/ring.sh.length)[:,None]*np.array([[-1.,0.,1.]])-(x*slope)[:,None]*Nr
        self.V,self.Gr,self.Gz=operator(N),operator(Nr),operator(Nz)
        if np.max(abs(self.V@self.hat-ring.phi))>1e-10:
            raise ValueError('Refinement changed the physical probe')
        self.zero=np.zeros((self.hat.size,2))
        self.reference=None
        self.eval(self.zero)

    def eval(self,U):
        a=self.ring
        U=U.reshape(-1,2)
        offset=self.V@U
        r,z=a.r+offset[:,0],a.z+offset[:,1]
        R,slope=a.sh.radius(z)
        x=r/R
        if np.max(x)>1.000000001 or z.max()>a.sh.length:
            raise ValueError('Local patch left material domain')
        B=a.B+np.stack([self.Gr@U,self.Gz@U],axis=-1)
        J=a.j0*np.linalg.det(B)*(r/a.rbase)**2
        if min(J)<=0:
            raise ValueError('Inverted trial')
        hoop=a.refhoop*r/a.rbase
        K=B@a.K@B.transpose(0,2,1)
        W=.5*(K[:,0,0]+K[:,1,1]+2*hoop**2-4)-np.log(J)+.25*(J-1)**2
        chi,cx,ct=map(sum,zip(*(interpolation(*grid,z/a.sh.length,x) for grid in a.order)))
        cw=a.w*chi
        Wchi=W*chi
        if self.reference is None:self.reference=Wchi.copy()
        energy=a.w@(Wchi-self.reference)
        H=.5*J*(J-1)
        P=(B@a.K+(H-1)[:,None,None]*np.linalg.inv(B).transpose(0,2,1))
        direct=np.column_stack((cx/R,ct/a.sh.length-cx*x*slope/R))*(a.w*W)[:,None]
        direct[:,0]+=2*cw*(hoop**2-1+H)/r
        grad=self.V.T@direct+self.Gr.T@(cw[:,None]*P[:,:,0])+self.Gz.T@(cw[:,None]*P[:,:,1])
        self.last=(K,hoop,J)
        return float(energy),grad.ravel()

    def sample(self,U):
        energy,g=self.eval(U)
        K,hoop,J=self.last
        lam=np.linalg.eigvalsh(K)
        smallest=np.minimum(np.sqrt(lam[:,0]),hoop)
        site=int(np.argmin(smallest))
        mapped=self.V@U.reshape(-1,2)
        radius=self.ring.r[site]+mapped[site,0]
        depth=self.ring.z[site]+mapped[site,1]
        return dict(energy_change=energy,max_force=float(abs(g).max()),
                    original_probe_force=(-self.hat@g.reshape(-1,2)).tolist(),
                    min_stretch=float(min(np.sqrt(lam[:,0].min()),hoop.min())),
                    min_volume_ratio=float(J.min()),max_displacement=float(abs(U).max()),
                    min_stretch_site=[float(radius),float(depth/self.ring.sh.length)])


def descend(patch, metric):
    """Backtrack every material trial against domain, orientation and energy."""
    v=np.zeros(metric.shape[1])
    inverse=np.eye(v.size)
    e,gradient=patch.eval(metric@v)
    g=metric.T@gradient
    for iteration in range(350):
        if np.max(abs(gradient))<1e-6:
            return metric@v,iteration,'material force below 1e-6'
        direction=-inverse@g
        slope=g@direction
        if slope>=0:
            direction=-g
            slope=-g@g
        scale=1.
        while scale*np.max(abs(metric@direction))>3e-15:
            step=scale*direction
            try:
                trial,raw=patch.eval(metric@(v+step))
            except ValueError:
                scale*=.5
                continue
            if trial<=e+1e-4*scale*slope:
                break
            scale*=.5
        else:
            return metric@v,iteration,'no resolved feasible energy decrease'
        new=metric.T@raw
        y=new-g
        sy=step@y
        if sy>1e-12*np.linalg.norm(step)*np.linalg.norm(y):
            Hy=inverse@y
            inverse+=(1.+(y@Hy)/sy)*np.outer(step,step)/sy
            inverse-=(np.outer(step,Hy)+np.outer(Hy,step))/sy
        v+=step
        e,g,gradient=trial,new,raw
    return metric@v,iteration+1,'iteration limit'


def run(state,checkpoint,divisions,integration,extend):
    data=source(state,checkpoint)
    sh,d,_,_=data
    tables=order_tables(state,checkpoint,data,None)
    results=[]
    for q in integration:
        ring,offsets=rim_patch(state,sh,d,tables,q,extend)
        for n in divisions:
            patch=Patch(ring,n,offsets)
            start=patch.sample(patch.zero)
            U=patch.zero.ravel()
            # Hessian of the smooth positive elastic part gives a local coordinate
            # metric. L-BFGS line searches the full registered material energy.
            h=1e-7
            basis=patch.V.astype(bool)
            graph=sparse.kron(basis.T@basis,np.ones((2,2),dtype=bool),format='csc')
            colors=group_columns(graph)
            H=np.zeros((len(U),len(U)))
            for color in np.unique(colors):
                direction=(colors==color).astype(float)
                column=(patch.eval(U+h*direction)[1]-patch.eval(U-h*direction)[1])/(2*h)
                for k in np.flatnonzero(direction):
                    ids=graph.indices[graph.indptr[k]:graph.indptr[k+1]]
                    H[ids,k]=column[ids]
            eigen,V=np.linalg.eigh((H+H.T)/2)
            metric=V/np.sqrt(np.maximum(abs(eigen),1.))
            final,iterations,stop=descend(patch,metric)
            end=patch.sample(final)
            results.append(dict(subdivisions=n, integration=q,free_dofs=len(U),
                                tangent_min=float(eigen.min()),initial=start,final=end,
                                optimizer=stop,iterations=iterations))
            print(json.dumps(results[-1],allow_nan=False),flush=True)
    return dict(checkpoint=file_sha(checkpoint),material=file_sha(state),results=results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('checkpoint',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--divisions',type=int,nargs='+',default=[1,2,4])
    p.add_argument('--integration',type=int,nargs='+',default=[8,16])
    p.add_argument('--extend-perimeter',action='store_true')
    a=p.parse_args()
    out=run(a.state,a.checkpoint,a.divisions,a.integration,a.extend_perimeter)
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
