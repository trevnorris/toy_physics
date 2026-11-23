"""Independent sharp disk traction and weighted volume/face work on a clamped R4 disk.

Both sides share a disk displacement equal to one across the interior and
tapering to zero at its edge. A physical trace at the disk, a sharp cylinder,
and a finite-height virtual-work probe are different observables. Reconstruct
the Q1 fields in physical coordinates and separate their contributions.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_closed_rim_surface_review import Reader, quadrature, stress_energy
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_ring_surface_review import FieldChart
from mouth_last_rings_fluid import GRIDS, read_source
from mouth_supercore_joined_front import check_guard, emit, file_sha


class Disk:
    def __init__(self, state, wall, fields, original=False):
        x,d,_,_=read_source(state,wall)
        base=Front(state,GRIDS['native'],1)
        reference=Trial(base,state,x,(d.shape[1]-1)//3,d.shape[0]-1,False)
        self.sh=reference.shape
        check_guard(fields)
        rows=[json.loads(line) for line in (fields.parent/'stdout').read_text().splitlines()]
        digest=file_sha(fields)
        if not any(r.get('kind')=='fields' and r['sha256']==digest for r in rows):
            raise ValueError('Unverified matched fields')
        with np.load(fields) as doc:
            self.level=doc['level'].item()
            if original:
                if fields.parent.parent.name!='r4_last_rings_fluid01' or 'load_sha' in doc:
                    raise ValueError('Expected archived original load')
                load='archived native material deposition'
            else:
                load=doc['load_sha'].item()
            self.source_hash=dict(field=digest,wall=file_sha(wall),load=load)
            if tuple(doc['grid'])!=GRIDS[self.level] or doc['checkpoint_sha'].item()!=file_sha(wall):
                raise ValueError('Reconstructed topology or wall differs')
            m=FieldChart(self.sh,base.c,1,GRIDS[self.level],precondition=False,
                         terminal_profile=reference.model.terminal_profile)
            m.rr=tuple(doc[k] for k in ('inside_x','outside_x'))
            m.z=tuple(doc[k] for k in ('feed_z','duct_z','receiver_z'))
            m.a=doc['angles']
            self.reader=Reader(m,base.c,doc['wound'],doc['rest'])
        self.model=m
        orig=reference.model
        self.native_r=orig.rr[0]
        u=np.clip((self.native_r-.90)/.10,0.,1.)
        self.native_eta=1.-u*u*(3.-2*u)
        L=self.sh.length
        lower=np.flatnonzero(np.isclose(orig.z[1],L-orig.h,rtol=0,atol=1e-11)).item()
        upper=np.flatnonzero(np.isclose(orig.z[2],L+orig.h,rtol=0,atol=1e-11)).item()
        self.windows=[[orig.z[1][lower-k],orig.z[1][lower-k+1],
                       orig.z[2][upper+k-1],orig.z[2][upper+k]] for k in (2,3)]
        self.sharp_heights=[.16,.235]

    def weight(self,x):
        k=np.clip(np.searchsorted(self.native_r,x)-1,0,len(self.native_r)-2)
        return (np.interp(x,self.native_r,self.native_eta),
                (self.native_eta[k+1]-self.native_eta[k])/(self.native_r[k+1]-self.native_r[k]))

    def physical_state(self,slab,x,a,z):
        m=self.model
        x,a,z=np.broadcast_arrays(x,a,z)
        x,a,z=(t.ravel() for t in (x,a,z))
        b,slope=self.sh.radius(z)
        r=x*b
        F,Fp=self.reader.vortex(r*np.sin(a))
        # All probes are on the interior or its bulk continuation. The clamped
        # radial lift is identically zero on the exterior side of the wall.
        grids=(m.rr[0],m.a,m.z[slab])
        indices=[np.clip(np.searchsorted(g,q,side='right')-1,0,len(g)-2) for g,q in zip(grids,(x,a,z))]
        isx,ia,iz=indices
        gap=[g[k+1]-g[k] for g,k in zip(grids,indices)]
        lam=[(q-g[k])/h for g,k,h,q in zip(grids,indices,gap,(x,a,z))]
        diffs=[[-1/h,1/h] for h in gap]
        weights=[(1-p,p) for p in lam]
        for f in self.reader.fields[slab,0]:
            vals=np.zeros((len(x),3))
            grad=np.zeros((len(x),3,3))
            for i,j,k in np.ndindex(2,2,2):
                y=f.values[isx+i,ia+j,iz+k]
                a1,b1,c1=(arr[side] for arr,side in zip(weights,(i,j,k)))
                vals+=(a1*b1*c1)[:,None]*y
                grad[:,:,0]+=(diffs[0][i]*b1*c1/b)[:,None]*y
                grad[:,:,1]+=(a1*diffs[1][j]*c1/r)[:,None]*y
                grad[:,:,2]+=(a1*b1*diffs[2][k])[:,None]*y
            grad[:,:,2]-=(x*slope)[:,None]*grad[:,:,0]
            amp,phi,chi=vals.T
            grad[:,0]=F[:,None]*grad[:,0]+(amp*Fp)[:,None]*np.stack((np.sin(a),np.cos(a),0*a),axis=1)
            n=(amp*F)**2
            yield (grad[:,0],grad[:,1],grad[:,2],n,chi,1/(r*np.sin(a))**2)

    def stresses(self,slab,x,a,z):
        # Complex-step energy derivatives measure both physical fluxes. No
        # material stress or interpolated front-pin reaction is folded into them.
        increment=1j*1e-25
        axial=[0.,0.]
        for sign,state in zip((1.,-1.),self.physical_state(slab,x,a,z)):
            gs=list(state[:3])
            for axis in (0,2):
                jac=1+(increment if axis==2 else 0)
                mapped=[]
                for g in gs:
                    v=g.astype(complex)
                    v[:,axis]-=increment*g[:,2]/jac
                    mapped.append(v)
                axial[int(axis==2)]-=sign*np.imag(jac*stress_energy(*mapped,*state[3:]))/increment.imag
        return axial

    def verify_interpolation(self):
        errors=[]
        for slab,z in ((1,self.sh.length-.087),(2,self.sh.length+.07)):
            x=np.array([.325,.799,.945,.992])
            a=np.array([.412,.754,.566,1.004])
            radial,axial=self.stresses(slab,x,a,z)
            r=x*self.sh.radius(z)[0]
            for comp,axis in ((radial,0),(axial,2)):
                n=np.eye(3)[axis]
                ref=self.reader.sample(slab,0,r,a,z,n)[:,1]
                errors.append(np.max(abs(comp-ref)))
        if max(errors)>3e-10:
            raise ValueError(f'Independent physical interpolation differs {errors}')
        return float(max(errors))

    def cap(self,slab,z,order):
        x,wx=quadrature(self.model.rr[0],order)
        a,wa=quadrature(self.model.a,order)
        X,A=np.meshgrid(x,a,indexing='ij')
        R=self.sh.radius(z)[0]
        area=4*np.pi*(X*R)**2*R*np.sin(A)*wx[:,None]*wa
        _,tz=self.stresses(slab,X,A,z)
        q=self.weight(X)[0]
        return np.array([area.ravel()@tz,(area*q).ravel()@tz])

    def run(self,order):
        L=self.sh.length
        lower=self.cap(1,np.nextafter(L,-np.inf),order)
        upper=self.cap(2,np.nextafter(L,np.inf),order)
        # A true physical surface discontinuity is the value immediately above
        # minus the value immediately below, for an increasing-w orientation.
        emit(kind='terminal_trace',quadrature=order,
             aperture_and_clamped_lower=lower.tolist(),aperture_and_clamped_upper=upper.tolist(),
             aperture_and_clamped_jump=(upper-lower).tolist())
        x,wx=quadrature(self.model.rr[0],order)
        a,wa=quadrature(self.model.a,order)
        X,A=np.meshgrid(x,a,indexing='ij')
        q,dq=self.weight(X)
        weights=4*np.pi*X**2*np.sin(A)*wx[:,None]*wa
        windows=np.asarray([*self.windows,[L-self.model.h,L,L,L+self.model.h]])
        # [C, S] are work from closed axial cap tapers and from radial taper
        # inside the volume. Strong enclosing forces are also read independently.
        weak=np.zeros((len(windows),2))
        sharp_grad=np.zeros(2)
        for slab in (1,2):
            z0,z1=(windows[1,0],L) if slab==1 else (L,windows[1,3])
            zg=self.model.z[slab]
            cuts=np.r_[zg,windows.ravel(),L+np.array(self.sharp_heights),L-np.array(self.sharp_heights)]
            zs,wz=quadrature(np.r_[z0,cuts[(cuts>z0)&(cuts<z1)],z1],order)
            for z,dz in zip(zs,wz):
                r,slope=self.sh.radius(z)
                tr,tz=self.stresses(slab,X,A,z)
                ia=(weights*q).ravel()@tz*r**3*dz
                igrad=(weights*dq).ravel()@(tr-X.ravel()*slope*tz)*r**2*dz
                for j,(l0,l1,u0,u1) in enumerate(windows):
                    alpha=np.clip((z-l0)/(l1-l0),0,1)*np.clip((u1-z)/(u1-u0),0,1)
                    dalpha=(1/(l1-l0) if l0<z<l1 else -1/(u1-u0) if u0<z<u1 else 0.)
                    weak[j]-=[dalpha*ia,alpha*igrad]
                for j,h in enumerate(self.sharp_heights):
                    if abs(z-L)<h:
                        sharp_grad[j]-=igrad
        for j,h in enumerate(self.sharp_heights):
            lo=self.cap(1,L-h,order)[1]
            hi=self.cap(2,L+h,order)[1]
            _,whole=self.reader.integrate(4.55,h,order)
            slit=self.reader.sealed_increment(0.,h,order)[1]
            emit(kind='weighted_surface',quadrature=order, half_height=h, radius=4.55,
                 old_whole_physical=whole[1]+slit,sealed_slit=slit,
                 disk_low=-lo,disk_high=hi,clamped_closed=hi-lo,
                 rim_complement_closed=whole[1]+slit-hi+lo,
                 radial_work=sharp_grad[j],clamped_configurational_work=hi-lo+sharp_grad[j])
        for j,row in enumerate(weak):
            emit(kind='weak_sector',quadrature=order,window=('near','far','short')[j],
                 axial_cap_taper=row[0],interior_radial_gradient=row[1],convected_field_work=sum(row))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','wall','fields'):p.add_argument(name,type=Path)
    p.add_argument('--original-load',action='store_true')
    p.add_argument('--quadrature',type=int,nargs='+',default=[2,4])
    args=p.parse_args()
    obj=Disk(args.state,args.wall,args.fields,args.original_load)
    emit(kind='review_inputs',**obj.source_hash,level=obj.level,windows=obj.windows,
         independent_pointwise_difference=obj.verify_interpolation())
    for n in args.quadrature:obj.run(n)
