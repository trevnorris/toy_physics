"""Axial virtual work on both sides of the held bonded R4 disk.

A test displacement eta(r,w) e_w is supported a distance h either side of
the disk. The inner material is convected at its native quadrature. Drive
and rest share the same lattice, shape and physical test function. Surface
traces and volume weak work are evaluated separately.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator

from mouth_disk_mesh_rematch import digest, inputs
from mouth_bonded_front_adjoint import angular_weights, get_material_load, set_front
from mouth_supercore_joined_front import Joined, emit, KAPPA, ORDER_INERTIA, WELL
from mouth_wound_carrier import straight_vortex


def dot(a,b):
    return np.einsum('...i,...i->...',a,b)


def smooth(t):
    y = np.clip(t,0.,1.)
    return 1-y*y*(3-2*y), -6*y*(1-y)


def lift(shape,r,z,h,compact):
    y = (z-shape.length)/h
    az,daz = smooth(abs(y))
    daz *= np.sign(y)/h
    b,slope = shape.radius(z)
    if compact:
        radial,dr = smooth((r/b-.75)/.2)
        dr /= .2*b
    else:
        radial,dr = np.ones_like(r),np.zeros_like(r)
    return az*radial, az*dr, daz*radial-az*dr*(r/b)*slope


def physical(winding,r,angle,values,grads,F,Fp):
    t,_,chi = values
    gv,gp,gc = grads
    gf = F[:,None]*gv + (t*Fp)[:,None]*np.column_stack((np.sin(angle),np.cos(angle),0*r))
    n = (F*t)**2
    swirl = (winding/(r*np.sin(angle)))**2
    return gf,gp,gc,n,swirl,chi


def lag(q):
    gf,gp,gc,n,swirl,chi = q
    den = n+ORDER_INERTIA*dot(gc,gc)
    vel = gp - (ORDER_INERTIA*dot(gp,gc)/den)[:,None]*gc
    return (np.stack((.5*(dot(gf,gf)+n*swirl)+(n**5-5*n+4)/20,
                      .5*n*dot(gp,vel), .5*KAPPA*dot(gc,gc)+WELL*(chi*(1-chi))**2)),vel)


def stresses(q):
    gf,_,gc,n,_,_ = q
    energy,v = lag(q)
    # Each column is the axial row of the corresponding canonical momentum
    # flux. Eliminating order momentum turns j_i p_z-pi v_i chi_z into n v_i v_z.
    shear = np.stack((gf[:,0]*gf[:,2],n*v[:,0]*v[:,2],KAPPA*gc[:,0]*gc[:,2]))
    axial = np.stack((gf[:,2]**2,n*v[:,2]**2,KAPPA*gc[:,2]**2))-energy
    return shear,axial


def transformed(q,gr,gz,eps):
    # Pull back physical scalar fields by (r,a,z) -> (r,a,z+eps*eta).
    def grad(v):
        out=v.copy()
        out[:,2]=v[:,2]/(1+eps*gz)
        out[:,0]-=eps*gr*out[:,2]
        return out
    return lag(tuple(grad(v) for v in q[:3])+q[3:])[0].sum(axis=0)*(1+eps*gz)


class FieldTraces:
    def __init__(self,m,order,pair):
        self.m,self.pair=m,pair
        self.vortex=straight_vortex(m.winding,28.)
        self.nodal={}
        interp=RegularGridInterpolator(tuple(np.linspace(0,1,v) for v in order.shape),order)
        for slab in (1,2):
            x,_,z=np.meshgrid(m.rr[0],m.a,m.z[slab],indexing='ij')
            ref = interp((z/m.shape.length,x)) if slab==1 else m.terminal_profile(x,z-m.shape.length)
            self.nodal[slab]=[]
            for u in pair:
                nodal=u[:,m.ids[slab,0]].copy()
                nodal[2]+=ref
                self.nodal[slab].append(nodal)

    def sample(self,slab,x,a,z):
        m=self.m
        points=[v.ravel() for v in np.broadcast_arrays(x,a,z)]
        axis=(m.rr[0],m.a,m.z[slab])
        idx=[np.clip(np.searchsorted(ax,v,side='right')-1,0,len(ax)-2) for ax,v in zip(axis,points)]
        step=[np.diff(ax)[i] for ax,i in zip(axis,idx)]
        local=[(v-ax[i])/h for v,ax,i,h in zip(points,axis,idx,step)]
        saved=self.nodal[slab]
        values=np.zeros((2,3,len(points[0])))
        gradients=np.zeros((*values.shape,3))
        for bits in np.ndindex(2,2,2):
            f=[u if b else 1-u for b,u in zip(bits,local)]
            i,j,k=(index+bit for index,bit in zip(idx,bits))
            vals=np.stack([v[:,i,j,k] for v in saved])
            values+=vals*np.prod(f,axis=0)
            for j in range(3):
                gradients[...,j]+=vals*((2*bits[j]-1)*f[(j+1)%3]*f[(j+2)%3]/step[j])
        x,a,z=points
        b,s=m.shape.radius(z)
        if slab==2:s=0*s
        r=x*b
        gradients[...,0]/=b
        gradients[...,1]/=r
        gradients[...,2]-=x*s*gradients[...,0]
        F,Fp=self.vortex(r*np.sin(a))
        return [physical(m.winding,r,a,v,g,F,Fp) for v,g in zip(values,gradients)]


def gauss_cells(axis,n):
    g,w=leggauss(n)
    return ((axis[:-1,None]+np.diff(axis)[:,None]*(g+1)/2).ravel(),
            (np.diff(axis)[:,None]*w/2).ravel())


class Pillbox:
    def __init__(self,m,order,body,displacement,pair):
        self.m,self.b,self.d,self.pair=m,body,displacement,pair
        self.faces=FieldTraces(m,order,pair)
        wt,r,ang,z,F,Fp,*_=m.D.T
        b,_=m.shape.radius(z)
        self.mask=(r<b)&(abs(z-m.shape.length)<m.h+1e-12)
        mask=self.mask
        self.wt,self.r,self.z=wt[mask],r[mask],z[mask]
        self.field=[]
        for U in pair:
            vals=np.array([m.V[mask]@u for u in U])
            grads=np.array([np.column_stack([G[mask]@u for G in m.G]) for u in U])
            vals[2]+=m.ref[mask,0];grads[2]+=m.ref[mask,1:]
            self.field.append(physical(m.winding,r[mask],ang[mask],vals,grads,F[mask],Fp[mask]))
        self.T=tuple(a-b for a,b in zip(stresses(self.field[0]),stresses(self.field[1])))
        nodes=displacement.reshape(-1,2)[body.id]
        val=np.einsum('qi,qij->qj',body.basis,nodes)
        B=np.eye(2)+np.stack([np.einsum('qi,qij->qj',G,nodes) for G in (body.dr,body.dw)],axis=-1)
        self.sr,self.sz = body.r+val[:,0],body.w+val[:,1]
        self.K = B@body.K@B.transpose(0,2,1)
        self.hoop = body.h*self.sr/body.r
        self.J=body.j*np.linalg.det(B)*(self.sr/body.r)**2
        delta=(pair[0][2]-pair[1][2])[m.ids[1,0]]
        mean=np.einsum('raj,a->rj',delta,angular_weights(m.a))
        self.dc=RegularGridInterpolator((m.rr[0],m.z[1]),mean)((self.sr/m.shape.radius(self.sz)[0],self.sz))
        self.crystal_weight=body.measure*self.dc
        self.P=self.K[:,1,1]-1+.5*self.J*(self.J-1)
        self.ZR=self.K[:,1,0]
        self.material_grad=[]
        coeff=displacement[body.mask]*body.scale
        for U in pair:
            set_front(body,m,U)
            body.energy(coeff)
            self.material_grad.append(body.native.copy())
        # Order gradients for the full domain, with the same bonded material
        # deposition as the solved field and free and constrained nodes retained.
        load=get_material_load(body,coeff,m)
        point_source=m.V.T@(m.hq*m.D[:,-1])
        self.order_grad=[Joined.evaluate(m,U)[1][2]-point_source+load for U in pair]

    def solid_signed(self,gr,gz,eps):
        det=1+eps*gz
        if np.min(det)<=0:raise ValueError('Inverted virtual deformation')
        # For a fixed reference-volume crystal integral the Jacobian changes
        # in W, while its reference measure and the convected order stay fixed.
        k=self.K
        tr=k[:,0,0]*(1+(eps*gr)**2)+2*eps*gr*det*k[:,0,1]+det**2*k[:,1,1]+2*self.hoop**2
        j=self.J*det
        return self.crystal_weight@(.5*(tr-4)-np.log(j)+.25*(j-1)**2)

    def sides(self,slab,h,compact,n):
        m=self.m
        xs,wx=gauss_cells(m.rr[0],n)
        a,wa=gauss_cells(m.a,n)
        xx,aa=np.meshgrid(xs,a,indexing='ij')
        L=m.shape.length
        zz=L+np.zeros_like(xx)
        qs=self.faces.sample(slab,xx,aa,zz)
        A,C=stresses(qs[0]),stresses(qs[1])
        b=m.shape.exit
        weights=(4*np.pi*(b*xx)**2*np.sin(aa)*b*np.outer(wx,wa)).ravel()
        eta=lift(m.shape,b*xx.ravel(),zz.ravel(),h,compact)[0]
        axial=(A[1]-C[1])@(weights*eta)
        if compact:return axial,np.zeros(3)
        z0=L-h if slab==1 else L
        z1=L if slab==1 else L+h
        axis=m.z[slab]
        z,wz=gauss_cells(np.r_[z0,axis[(axis>z0+1e-12)&(axis<z1-1e-12)],z1],n)
        aa,zz=np.meshgrid(a,z,indexing='ij')
        b,s=m.shape.radius(zz)
        if slab==2:s=np.zeros_like(s)
        wall=self.faces.sample(slab,np.ones_like(aa),aa,zz)
        T,R=stresses(wall[0]),stresses(wall[1])
        eta=lift(m.shape,b.ravel(),zz.ravel(),h,False)[0]
        wt=(4*np.pi*b*b*np.sin(aa)*np.outer(wa,wz)).ravel()*eta
        return axial,((T[0]-R[0])-(T[1]-R[1])*s.ravel())@wt

    def audit(self,h,compact,n):
        m,b=self.m,self.b;L=m.shape.length
        if not 0<h<=m.h+1e-12:raise ValueError('Test extends outside the resolved finite band')
        eta,gr,gz=lift(m.shape,self.r,self.z,h,compact)
        _,sr,sz=lift(m.shape,self.sr,self.sz,h,compact)
        solid=-(self.crystal_weight@(self.ZR*sr+self.P*sz))
        stress=self.T[0]*gr+self.T[1]*gz
        out=[]
        for slab,sign in ((1,1),(2,-1)):
            side=(self.z<L) if slab==1 else (self.z>L)
            force=stress[:,side]@self.wt[side]
            raw,wall=self.sides(slab,h,compact,n)
            crystal=solid if slab==1 else 0.
            # The crystal ends at the disk. Its drive-minus-rest lateral
            # traction is zero because the order trace difference is zero.
            rim_delta=np.max(abs((self.pair[0][2]-self.pair[1][2])[m.ids[1,0][-1]]))
            if rim_delta>1e-13:raise ValueError('Need the drive-induced solid side traction on this boundary')
            out.append(dict(side=slab,raw_cap=raw.tolist(),wall=wall.tolist(),volume=force.tolist(),
                            crystal_volume=crystal,weak_cap=sign*(force.sum()+crystal-wall.sum())))
        signed={}
        for step in (2.e-5,1.e-5):
            def energy(eps):
                delta=transformed(self.field[0],gr,gz,eps)-transformed(self.field[1],gr,gz,eps)
                return self.wt@delta+self.solid_signed(sr,sz,eps)
            signed[str(step)]=float((energy(-step)-energy(step))/(2*step))
        analytic=self.wt@stress.sum(axis=0)+solid
        if max(abs(v-analytic) for v in signed.values())>1e-4*max(1,abs(analytic)):
            raise ArithmeticError(f'Axial stress and signed energy derivative disagree: {analytic}, {signed}')
        xx,tt=np.meshgrid(b.x,b.t)
        rz=b.chart.shape.radius(tt*L)[0]*xx+self.d[...,0]
        wz=tt*L+self.d[...,1]
        material_lift=lift(m.shape,rz,wz,h,compact)[0]
        reaction=(self.material_grad[0]-self.material_grad[1])[...,1]*material_lift
        free=b.mask[...,1]
        # Any supported prescribed order nodes are the side rim. Report the
        # available nodal reactions themselves; a normal slope has two traces.
        pinned=[]
        for slab in (1,2):
            ids=m.ids[slab,0]
            xx,_,zz=np.meshgrid(m.rr[0],m.a,m.z[slab],indexing='ij')
            lifts=lift(m.shape,xx*m.shape.radius(zz)[0],zz,h,compact)[0]
            pmask=~m.masks[2,ids] & (lifts>1e-13)
            pinned.extend(ids[pmask])
        unique=np.unique(np.asarray(pinned,dtype=int))
        order_delta=self.order_grad[0]-self.order_grad[1]
        return dict(kind='pillbox',h=h,radial='interior_to_.95' if compact else 'full',
                    quad=n,raw_and_weak=out,weak_gap=out[0]['weak_cap']-out[1]['weak_cap'],
                    volume_work=analytic,signed_volume_work=signed,
                    solid_dirichlet_axial_work=float(reaction[~free].sum()),
                    solid_free_axial_work=float(reaction[free].sum()),
                    supported_fixed_order_nodes=len(unique),fixed_order_delta_gradient_l1=abs(order_delta[unique]).sum())


def run(args):
    base,body,d,pair=inputs(args.state,args.joint,args.fields)
    model=base.joined(base.sh,False)
    pill=Pillbox(model,base.c,body,d,pair)
    emit(kind='preparation',grid=base.grid,winding=base.N,source_sha256=digest(args.state),
         joint_sha256=digest(args.joint),fields_sha256=digest(args.fields),
         curvature_drive_minus_rest=0,component_order=['amplitude_swirl_and_well','phase_and_order_inertia','order_gradient_and_well'])
    for h in (.12*2/7,.12*4/7,.12):
        for compact in (False,True):
            emit(**pill.audit(h,compact,8))
    emit(kind='surface_quadrature_replica',quad=12,readout=pill.sides(1,.12,False,12)[0].tolist(),
         receiving=pill.sides(2,.12,False,12)[0].tolist())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('state','joint','fields'):
        ap.add_argument(key,type=Path)
    run(ap.parse_args())
