"""Separate the cell-interior and face-jump terms in the saved axial momentum balance.

On each bank two compact lifts have identical wall trace and different radial
support. Assemble the mapped axial Piola flux of the drive-minus-rest fields
from each cell's own Q1 slope, including angular flux in the vortex geometry.
Differentiate inside that cell; evaluate the jumps independently on its facets.
"""
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_disk_pillbox import lag
from mouth_rim_disk_work_review import Disk
from mouth_supercore_joined_front import KAPPA, file_sha


class Cells:
    def __init__(self,disk,bank,interval,edges):
        self.disk=disk
        self.bank=bank
        self.interval=np.array(interval)
        self.edges=np.array(edges)
        ax,angle,z=disk.model.rr[bank],disk.model.a,disk.model.z[1]
        xs=np.flatnonzero(ax[:-1]<max(edges)-1e-12) if bank else np.flatnonzero(ax[1:]>min(edges)+1e-12)
        zs=np.flatnonzero((z[:-1]>=interval[0]-1e-12)&(z[1:]<=interval[-1]+1e-12))
        self.axes=(ax,angle,z)
        self.ids=np.array(np.meshgrid(xs,np.arange(len(angle)-1),zs,indexing='ij')).reshape(3,-1).T
        self.widths=np.array([np.diff(arr)[indices] for arr,indices in zip(self.axes,self.ids.T)]).T
        self.lower=np.array([arr[indices] for arr,indices in zip(self.axes,self.ids.T)]).T
        self.bits=np.array(list(product((0,1),repeat=3)))
        indices=(self.ids[:,None,:]+self.bits)
        self.nodal=[interp.values[tuple(indices.transpose(2,0,1))] for interp in disk.reader.fields[1,bank]]
        self.wall_index=0 if bank else len(ax)-2
        self.theta_top=len(angle)-2

    def flux(self,u):
        u=np.asarray(u)
        X=self.lower+u*self.widths
        x,a,z=X.T
        radius,slope=self.disk.sh.radius(z)
        width=28-radius if self.bank else radius
        r=radius+x*width if self.bank else x*width
        rz=slope*(1-x if self.bank else x)
        sn,cs=np.sin(a),np.cos(a)
        # The axis has zero oriented angular area and is never a volume point.
        # Use a harmless radius only in axis swirl energy; overwrite its flux below.
        axis=(sn==0)
        swirl=1./np.where(axis,1.,r*sn)**2
        F,Fp=self.disk.reader.vortex(r*sn)
        w=np.where(self.bits,u,1-u)
        N=np.prod(w,axis=1)
        dN=np.array([(2*self.bits[:,d]-1)*w[:,(d+1)%3]*w[:,(d+2)%3] for d in range(3)])
        net=np.zeros((len(r),3,3))
        for sign,nodes in zip((1.,-1.),self.nodal):
            vals=np.einsum('ijc,j->ic',nodes,N)
            g=np.einsum('ijc,kj->ick',nodes,dN)/self.widths[:,None,:]
            g[:,:,0]/=width[:,None]
            g[:,:,2]-=g[:,:,0]*rz[:,None]
            g[:,:,1]/=r[:,None]
            amp,_,chi=vals.T
            gf=F[:,None]*g[:,0]+(amp*Fp)[:,None]*np.column_stack((sn,cs,0.*r))
            gp,gc=g[:,1],g[:,2]
            n=(amp*F)**2
            state=(gf,gp,gc,np.where(axis,1.,n),swirl,chi)
            energies,vel=lag(state)
            for i,(v,weight) in enumerate(((gf,1.),(vel,n),(gc,KAPPA))):
                pref=weight*v[:,2]
                Tz=pref*v[:,2]-energies[i]
                net[:,i,0]+=sign*(pref*v[:,0]-rz*Tz)
                net[:,i,1]+=sign*width*pref*v[:,1]/r
                net[:,i,2]+=sign*width*Tz
        return net*(4*np.pi*r*r*sn)[:,None,None]

    def test(self,u):
        x,_,z=(self.lower+u*self.widths).T
        lo,mid,hi=self.interval
        A=np.maximum(0,np.minimum((z-lo)/(mid-lo),(hi-z)/(hi-mid)))
        dA=np.where(z<mid,1/(mid-lo),-1/(hi-mid))
        test=[]
        for edge in self.edges:
            if self.bank:
                g=np.maximum(0,1-x/edge);dg=np.where(x<edge,-1/edge,0.)
            else:
                g=np.maximum(0,(x-edge)/(1-edge));dg=np.where(x>edge,1/(1-edge),0.)
            test.append((g*A,dg*A,g*dA))
        return np.array(test)

    def integrate(self,degree):
        s,ws=leggauss(degree);s=(s+1)/2;ws=ws/2
        # Accumulate each sector and each extension separately. Normalize using
        # reference coordinate lengths, because flux contains the physical area.
        bulk=np.zeros((len(self.ids),2,3));weak=bulk.copy()
        dV=self.widths.prod(axis=1)
        for i,j,k in product(range(degree),repeat=3):
            u=s[[i,j,k]]
            wt=dV*ws[i]*ws[j]*ws[k]
            eta,dx,dz=self.test(u).transpose(1,2,0)
            flux=self.flux(u)
            weak+=wt[:,None,None]*(dx[:,:,None]*flux[:,None,:,0]+dz[:,:,None]*flux[:,None,:,2])
            div=np.zeros((len(wt),3))
            for h in range(3):
                step=np.zeros(3);step[h]=1e-4
                diff=(self.flux(u+step)[:,:,h]-self.flux(u-step)[:,:,h])/(2e-4*self.widths[:,h,None])
                div+=diff
            bulk-=wt[:,None,None]*eta[:,:,None]*div[:,None,:]
        wall=np.zeros_like(weak);jump=np.zeros_like(weak);sym=np.zeros_like(weak)
        for d in range(3):
            h1,h2=(d+1)%3,(d+2)%3
            for i,j in product(range(degree),repeat=2):
                wt=ws[i]*ws[j]*self.widths[:,h1]*self.widths[:,h2]
                for face in (0,1):
                    u=np.zeros(3);u[[d,h1,h2]]=[face,s[i],s[j]]
                    eta=self.test(u)[:,0].T
                    f=self.flux(u)[:,:,d]
                    surface=(2*face-1)*wt[:,None,None]*f[:,None,:]*eta[:,:,None]
                    iswall=(d==0 and face==1-self.bank) & (self.ids[:,0]==self.wall_index)
                    issym=(d==1 and face==1) & (self.ids[:,1]==self.theta_top)
                    wall[iswall]+=surface[iswall]
                    sym[issym]+=surface[issym]
                    jump[~(iswall|issym)]+=surface[~(iswall|issym)]
        # Other sides of this measured patch have identically zero lift.
        packed={name:arr.sum(axis=0) for name,arr in dict(volume_work=weak,one_sided_wall=wall,
                interior_face_jumps=jump,minus_cell_divergence=bulk,equatorial_trace=sym).items()}
        packed['ibp_residual']=packed['volume_work']-packed['one_sided_wall']-packed['interior_face_jumps']-packed['minus_cell_divergence']-packed['equatorial_trace']
        return packed


def main(a):
    d=Disk(a.state,a.wall,a.fields)
    bankdoc=json.loads(a.prior.read_text())
    prior=bankdoc['data'][d.level]
    assert prior['guard_exit']==0 and prior['source']['field']==file_sha(a.fields)
    print(json.dumps(dict(kind='source',hashes=d.source_hash,paired_record=file_sha(a.prior))),flush=True)
    for name in ['near','middle','upstream']:
        for side,bank in enumerate(('duct','outside')):
            p=[r for r in prior['results'] if r['strip']==name and r['bank']==bank and r['quadrature']==4]
            p.sort(key=lambda row:row['radial_cells'])
            if [r['radial_cells'] for r in p]!=[4,8]:raise ValueError('Missing paired wall probes')
            edges=[r['normalized_radial_bound'] for r in p]
            cells=Cells(d,side,p[0]['depth_knots'],edges)
            for n in a.degree:
                terms=cells.integrate(n)
                diff={k:(v[1]-v[0]).tolist() for k,v in terms.items()}
                row=dict(kind='budget',strip=name,bank=bank,degree=n,radial_cells=[4,8],
                         field_sectors=['amplitude and vortex','phase and order inertia','static order'],
                         terms={k:v.tolist() for k,v in terms.items()},zero_trace_wide_minus_narrow=diff,
                         corrected_solid=[q['solid_reaction'] for q in p],
                         independent_signed_total=[q['total_native'] for q in p])
                if n==4:
                    err=np.max(abs(terms['volume_work']-np.array([q['canonical_sector_reactions'] for q in p])))
                    if err>1e-8:raise ArithmeticError(f'Independent volume work disagreement {err}')
                    if abs(terms['ibp_residual']).max()>2e-4:
                        raise ArithmeticError(f'Unresolved volume and face quadrature {terms["ibp_residual"]}')
                print(json.dumps(row,allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','wall','fields','prior'):p.add_argument(key,type=Path)
    p.add_argument('--degree',nargs='+',default=[2,4],type=int)
    main(p.parse_args())
