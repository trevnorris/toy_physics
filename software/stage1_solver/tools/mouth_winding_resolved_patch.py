"""Resolve a stationary through-flow and free crystalline order in a finite exit patch.

Take both bath traces from the held sharp-front solution. The artificial cylindrical
boundary beyond the lip is matched to the receiver. A disk of zero thickness or a
prescribed disk current is never used in the patch Hamiltonian. This is a finite
Dirichlet subdomain test, not self-consistent iteration with the two far reservoirs.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import factorized
from scipy.special import expit

from mouth_breathing_mobile_patch import KAPPA, WELL, WIDTH, ORDER_INERTIA
from mouth_poly_shallow_release import emit
from mouth_sealed_shape_relax import Meridian
from mouth_wound_carrier import Wound, prepared


BITS = [(i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)]


class Patch:
    def __init__(self, path, radius, exit_radius, mesh, halfwidth):
        if radius not in (1, 2):
            raise ValueError('Reserved radius')
        guard = json.loads((path.parent/'outcome.json').read_bytes())
        limits = json.loads((path.parent/'limit-validation.json').read_bytes())
        if guard['exitCode'] or guard['childOutcome']['exitCode'] or not limits['verified'] or int(limits['actual']['memory.swap.max']):
            raise ValueError('Cannot adopt an unvalidated carrier')
        config, source_order, digest = prepared(radius)
        config[1] = np.log(exit_radius / (1.5 * radius))
        self.shape = Meridian(radius, config)
        self.b = exit_radius
        self.h = halfwidth
        with np.load(path, allow_pickle=False) as fields:
            if fields['material_sha256'].item() != digest:
                raise ValueError('Wrong order preparation')
            if abs(np.max(fields['disk'][:,3])/exit_radius-1) > .008:
                raise ValueError('Exit size differs from field geometry')
            winding = int(fields['winding']) if 'winding' in fields else 1
            if winding != 1:
                raise ValueError('Wrong circulation')
            model = Wound(self.shape, source_order, 1, int(fields['mesh']), int(fields['polar']),
                          float(fields['outer']), float(fields['roof']))
            self.make_grid(mesh)
            radial, angle, height = np.meshgrid(self.x, self.a, self.y,indexing='ij')
            r = self.shape.radius(height + self.shape.length)[0]*radial
            s = r*np.sin(angle)
            fref = model.vortex(s)[0]
            f = np.zeros_like(r)
            p = np.zeros_like(r)
            for part, select in ((1,height<=0), (2,height>0)):
                ids = model.grid[part,0]
                orig = (model.axial[part],model.radial[0],model.angles)
                pts = np.stack((height[select]+self.shape.length,radial[select],angle[select]),axis=-1)
                for out, source, scale in ((f, fields['amplitude_ratio'][model.did[ids]],fref[select]),
                                           (p, fields['phase'][ids],1.)):
                    out[select] = scale*RegularGridInterpolator(orig,source,bounds_error=True)(pts)
            pts = np.stack((np.clip((height+self.shape.length)/self.shape.length,0,1),radial),axis=-1)
            native_chi = RegularGridInterpolator(tuple(np.linspace(0,1,m) for m in source_order.shape),
                                                source_order)(pts)
            # The interior face of the slab has a material rim at half order.
            # Beyond it the artificial matching cylinder samples an undriven
            # receiver continuation. Varying the port spacing tests this closure.
            c = np.where(height<=0, native_chi, expit(-height/WIDTH))
            initial = np.stack([f,p,c])
            self.seed = np.stack([q.ravel()[self.unique] for q in initial])
            self.record = dict(source_path=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                               native_material_sha256=digest,entrance=radius,exit_radius=exit_radius,
                               azimuthal_winding=winding,mesh=mesh,half_port_distance=halfwidth,
                               sharp_intake=fields['disk'][:,0]@fields['disk'][:,4])
        self.set_boundaries()
        self.make_matrices()

    def make_grid(self, mesh):
        nr, na, ny = mesh
        self.x = np.linspace(0,1,nr+1)
        self.a = np.linspace(0,np.pi/2,na+1)
        self.y = np.linspace(-self.h,self.h,ny+1)
        self.mid = ny//2
        ids = np.arange((nr+1)*(na+1)*(ny+1)).reshape(nr+1,na+1,ny+1)
        ids[0] = ids[0,:1]
        _, self.unique, inv = np.unique(ids,return_index=True,return_inverse=True)
        self.ids = inv.reshape(ids.shape)
        self.size = len(self.unique)
        self.nodes = np.stack([self.ids[i:nr+i,j:na+j,k:ny+k].ravel() for i,j,k in BITS],axis=-1)
        self.dims = (nr,na,ny)

    def set_boundaries(self):
        ids = self.ids
        masks = np.ones((3,self.size),dtype=bool)
        masks[:,ids[:,:,0]] = False
        masks[:,ids[:,:,-1]] = False
        # The cylinder above the exit is an artificial cut through receiver material;
        # the wall below is impermeable with natural amplitude traction.
        masks[:,ids[-1,:,self.mid+1:]] = False
        masks[2,ids[-1,:,:self.mid+1]] = False
        # Regularity of a one-quantum amplitude on the circulation axis.
        masks[0,ids[:,0,:]] = False
        self.masks = masks
        self.slices = np.cumsum([0]+[int(sum(m)) for m in masks])

    def make_matrices(self):
        self.V, self.G = [], [[],[],[]]
        self.D = []
        nr,na,ny=self.dims
        dx,da,dy=self.x[1],self.a[1],self.y[1]-self.y[0]
        x0,a0,y0=np.meshgrid(self.x[:-1],self.a[:-1],self.y[:-1],indexing='ij')
        size = (len(self.nodes),self.size)
        rows = np.broadcast_to(np.arange(len(self.nodes))[:,None],self.nodes.shape)
        def matrix(values):
            return sparse.coo_matrix((np.broadcast_to(values,self.nodes.shape).ravel(),
                                      (rows.ravel(),self.nodes.ravel())),shape=size).tocsr()
        samples = (np.array([-1,1])/np.sqrt(3)+1)/2
        core_density_half = 1.547951
        for qx in samples:
            for qa in samples:
                for qy in samples:
                    fractions=np.array([qx,qa,qy])
                    factor=np.where(np.array(BITS),fractions,1-fractions)
                    signs=2*np.array(BITS)-1
                    b=factor.prod(axis=1)
                    dd=np.array([signs[:,i]*factor[:,(i+1)%3]*factor[:,(i+2)%3] for i in range(3)])
                    xx=(x0+qx*dx).ravel()
                    a=(a0+qa*da).ravel()
                    y=(y0+qy*dy).ravel()
                    radius,slope=self.shape.radius(self.shape.length+y)
                    rr=radius*xx
                    rr_grad=dd[0]/(dx*radius[:,None])
                    az_grad=dd[1]/(da*rr[:,None])
                    zz_grad=dd[2]/dy-(xx*slope)[:,None]*rr_grad
                    w=4*np.pi*rr**2*np.sin(a)*radius*dx*da*dy/8
                    self.V.append(matrix(b))
                    for store, grad in zip(self.G,(rr_grad,az_grad,zz_grad)):
                        store.append(matrix(grad))
                    s=rr*np.sin(a)
                    self.D.append(np.column_stack((w,xx,y,s,rr,s<core_density_half)))
        self.V=sparse.vstack(self.V).tocsr()
        self.G=[sparse.vstack(g).tocsr() for g in self.G]
        self.D=np.vstack(self.D)
        w=self.D[:,0]
        n=(self.V@self.seed[0])**2
        def metric(weights,bulk):
            S=self.V.T@self.V.multiply((w*bulk)[:,None])
            for grad in self.G:
                S+=grad.T@grad.multiply((w*weights)[:,None])
            return S
        self.solvers=[factorized(M[m][:,m].tocsc()) for m,M in zip(self.masks,
                         (metric(1.,1/self.D[:,3]**2+1.),metric(n,0.),metric(KAPPA,2*WELL)))]
        self.Nq=len(w)

    def unpack(self,x,seed):
        u=seed.copy()
        for i,m in enumerate(self.masks):
            u[i,m]=x[self.slices[i]:self.slices[i+1]]
        return u

    def vars(self,u):
        return np.concatenate([a[m] for a,m in zip(u,self.masks)])

    def norm(self,g):
        return np.concatenate([solve(g[a:b]) for solve,a,b in zip(self.solvers,self.slices,self.slices[1:])])

    def evaluate(self,u,read=False):
        w,_,_,s,_,_=self.D.T
        f,p,c=u
        f=self.V@f
        n=f*f
        chi=self.V@c
        gf,gp,gc=(np.column_stack([G@a for G in self.G]) for a in u)
        square=lambda a:np.einsum('ij,ij->i',a,a)
        gg=square(gc)
        proj=np.einsum('ij,ij->i',gp,gc)
        pi=ORDER_INERTIA*n*proj/(n+ORDER_INERTIA*gg)
        vel=gp-(ORDER_INERTIA*proj/(n+ORDER_INERTIA*gg))[:,None]*gc
        curr=n[:,None]*vel
        u0=(n**5-5*n+4)/20
        pot=WELL*(chi*(1-chi))**2
        E=.5*(square(gf)+n/s**2+KAPPA*gg+
              np.einsum('ij,ij->i',gp,curr))+u0+pot
        if read:
            t=n*vel[:,2]**2+gf[:,2]**2+KAPPA*gc[:,2]**2-E
            tx=n*vel[:,2]*vel[:,0]+gf[:,2]*gf[:,0]+KAPPA*gc[:,2]*gc[:,0]
            ta=n*vel[:,2]*vel[:,1]+gf[:,2]*gf[:,1]+KAPPA*gc[:,2]*gc[:,1]
            return curr,vel,chi,n,t,tx,ta,E
        # The derivative of minimized phase energy with respect to n is
        # (|u|²+2 pi u.grad(chi)/n)/2.
        gp2=square(gp)
        dep=gp2-ORDER_INERTIA**2*gg*proj**2/(n+ORDER_INERTIA*gg)**2
        gf0=f*(1/s**2 + (n**4-1)/2 + dep)
        gr=np.array([self.V.T@(w*gf0),np.zeros(self.size),
                     self.V.T@(w*(2*WELL*chi*(1-chi)*(1-2*chi)))])
        for i,G in enumerate(self.G):
            gr[0]+=G.T@(w*gf[:,i])
            gr[1]+=G.T@(w*curr[:,i])
            gr[2]+=G.T@(w*(KAPPA*gc[:,i]-pi*vel[:,i]))
        return w@E,gr

    def solve(self,rest):
        seed=self.seed.copy()
        if rest:
            seed[1]=0.
        x=self.vars(seed)
        E,gr=self.evaluate(seed)
        g=self.vars(gr)
        history=[]
        for iteration in range(600):
            inv=self.norm(g)
            err=max(abs(inv))
            if err<2e-7:
                break
            q=g.copy()
            a=[]
            for dx,dg,scale in reversed(history):
                coeff=scale*(dx@q)
                q-=coeff*dg
                a.append(coeff)
            dir=self.norm(q)
            for (dx,dg,scale),coeff in zip(history,reversed(a)):
                dir+=dx*(coeff-scale*(dg@dir))
            if dir@g<=0:
                dir=inv
                history=[]
            descent=dir@g
            step=min(1,.15/max(abs(dir)))
            while step>1e-9:
                d=-step*dir
                new=self.unpack(x+d,seed)
                val,gr1=self.evaluate(new)
                if val<=E-.0001*step*descent+2.e-14*max(1.,abs(E)):
                    break
                step*=.5
            else:
                raise ArithmeticError(f'Stationarity line search failed at {iteration} {err}')
            g1=self.vars(gr1)
            dg=g1-g
            curvature=d@dg
            if curvature>max(0,1.e-11*np.linalg.norm(d)*np.linalg.norm(dg)):
                history.append((d,dg,1/curvature))
                history=history[-12:]
            E,g,x=val,g1,x+d
        else:
            raise ArithmeticError(f'Unconverged stationarity {err}')
        return self.unpack(x,seed),dict(iterations=iteration,stationarity=err,energy=E)

    def measure(self,u,details):
        j,vel,c,n,T,Tx,Ta,E=self.evaluate(u,True)
        _,gr=self.evaluate(u)
        w,x,y,s,rr,core=self.D.T
        _,na,ny=self.dims
        dy=self.y[1]-self.y[0]
        flow=-(gr[1,np.unique(self.ids[:,:,0])].sum())
        receiver=gr[1,np.unique(self.ids[:,:,-1])].sum()
        side=gr[1,np.unique(self.ids[-1,:,self.mid+1:-1])].sum()
        # Vary the embedding in depth with a compact virtual displacement
        # vanishing on both ports and on the lateral boundary. The residual checks
        # integrated divergence of the full momentum flux, without face differentiation.
        radius,slope=self.shape.radius(self.shape.length+y)
        envelope=(1-x*x)**2
        dex=-4*x*(1-x*x)
        ay=(1-(y/self.h)**2)
        ayw=-2*y/self.h**2
        dr=dex/radius*ay
        dw=envelope*ayw-dex*x*slope/radius*ay
        weak_terms_x=w*Tx*dr
        weak_terms_w=w*T*dw
        scale=np.sum(abs(weak_terms_x)+abs(weak_terms_w))
        # Two separate planes, with the two Gauss layers nearest the original exit.
        # The stress is grand-canonical and includes pressure, order, quantum gradient
        # and all three meridional and azimuthal kinetic components.
        surf=[]
        for sign in (-1,1):
            choose=(sign*y>0)&(abs(y)<dy)
            area=w[choose]/dy
            jz=j[choose,2]
            m=abs(vel[choose,2])/n[choose]**2
            mt=np.sqrt(np.sum(vel[choose]**2,axis=1)+1/s[choose]**2)/n[choose]**2
            parts=[]
            for label,indices in [('all',np.ones_like(m,bool)),('normal_subsonic',m<1),
                                  ('total_subsonic',mt<1),('straight_core_interior',core[choose].astype(bool)),
                                  ('density_above_half',n[choose]>.5),('dense_total_subsonic',(mt<1)&(n[choose]>.5))]:
                a=area[indices]
                jj=jz[indices]
                tt=T[choose][indices]
                parts.append(dict(region=label,area=float(a.sum()),inward_current=float(a@jj),
                             positive_current=float(a@np.maximum(0.,jj)),negative_current=float(a@np.minimum(0.,jj)),
                             normal_stress_integral=float(a@tt),absolute_stress_integral=float(a@abs(tt))))
            surf.append(dict(sign=sign,minimum_normal_mach=m.min(),minimum_total_mach=mt.min(),
                             minimum_density=n[choose].min(),maximum_density=n[choose].max(),
                             minimum_order=c[choose].min(),maximum_order=c[choose].max(),parts=parts))
        return dict(**details,feed_inward=flow,receiver_axial=receiver,receiver_side=side,
                    port_current_error=receiver+side-flow,
                    order_range=[float(c.min()),float(c.max())],
                    compact_momentum_residual=float((weak_terms_w+weak_terms_x).sum()),
                    compact_momentum_scale=scale,
                    momentum_relative_residual=float(abs((weak_terms_w+weak_terms_x).sum())/scale),
                    disk=surf),T

    def run(self):
        emit(kind='preparation',**self.record,KAPPA=KAPPA,order_inertia=ORDER_INERTIA,well=WELL,
             matched_ports='Dirichlet density, order and phase at both finite cuts; free order at disk, fixed side rim')
        states=[]
        for rest in [False,True]:
            u,info=self.solve(rest)
            d,T=self.measure(u,info)
            emit(kind='resolved',rest=rest,**d)
            if abs(d['port_current_error'])>2.e-5 or not np.isfinite(d['compact_momentum_residual']):
                raise ArithmeticError('Port balance or stress evaluation failed')
            states.append((u,T))
        (u,T),(u0,T0)=states
        _,vel,_,n,_,Tx,_,_=self.evaluate(u,True)
        Tx0=self.evaluate(u0,True)[5]
        w,x,y,s,rr,core=self.D.T
        dy=self.y[1]-self.y[0]
        rad,slope=self.shape.radius(self.shape.length+y)
        env=(1-x*x)**2
        dex=-4*x*(1-x*x)
        win=1-(y/self.h)**2
        dw=env*(-2*y/self.h**2)-dex*x*slope/rad*win
        dr=dex/rad*win
        signed=w*(dw*(T-T0)+dr*(Tx-Tx0))
        norm=abs(w*dw*(T-T0))+abs(w*dr*(Tx-Tx0))
        emit(kind='driven_momentum_audit',drive_minus_rest_residual=float(sum(signed)),
             drive_minus_rest_scale=float(sum(norm)),
             relative_momentum_residual=float(abs(sum(signed))/sum(norm)))
        for sign in (-1,1):
            take=(y*sign>0)&(abs(y)<dy)
            ma=abs(vel[:,2])/n**2
            mt=np.sqrt(np.sum(vel**2,axis=1)+1/s**2)/n**2
            row={}
            for key,choose in [('all',take),('normal_subsonic',take&(ma<1)),('total_subsonic',take&(mt<1)),
                               ('straight_core_interior',take&core.astype(bool)),('density_above_half',take&(n>.5)),
                               ('dense_total_subsonic',take&(mt<1)&(n>.5))]:
                dF=(T-T0)[choose]
                row[key]=dict(flow_induced_normal_stress=(w[choose]/dy)@dF,
                              absolute_flow_induced_normal_stress=(w[choose]/dy)@abs(dF))
            emit(kind='matched_normal_stress',sign=sign,regions=row)


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('fields',type=Path)
    a.add_argument('radius',type=int,choices=[1,2])
    a.add_argument('exit',type=float)
    a.add_argument('--mesh',type=int,nargs=3,default=[10,8,32])
    a.add_argument('--halfwidth',type=float,default=.12)
    args=a.parse_args()
    if args.mesh[2]%2:
        raise ValueError('Exit plane must coincide with mesh interface')
    Patch(args.fields,args.radius,args.exit,args.mesh,args.halfwidth).run()
