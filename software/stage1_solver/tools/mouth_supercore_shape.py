"""Joint static shape and crystalline-order relaxation of a radius-four single-winding puncture."""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.optimize import minimize

from mouth_circulation_material import collar
from mouth_exchange_carrier import BIAS
from mouth_live_waist_fields import ORDER_INERTIA
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import CUT, Order, exterior, source
from mouth_sealed_shape_relax import DELTA, Meridian
from mouth_wound_carrier import Wound, basis


def preparation(n):
    """One native slab thickness, geometry-only interpolation and a prepared melted central pocket."""
    blank = Meridian(4.,np.zeros(4))
    q = np.array([-np.log(blank.length), np.log(4./blank.exit), 0., 0.])
    t = np.array([.25,.75])
    chart = Meridian(4.,q)
    widths = 4.*4*t*(1-t)
    q[2:] = np.linalg.solve(np.column_stack((widths,widths*(2*t-1))),4.-chart.profile(t)[0])
    shape = Meridian(4.,q)
    owner = Order(shape,n)
    # Start with an actual zero-order pocket. The initial guess has no place
    # in the final Hamiltonian. The flat feed, physical side and terminal traces persist.
    radial = owner.r/shape.neck
    depth = owner.w[:,None]/shape.length
    distance = np.sqrt((radial/.43)**2+((depth-.50)/.28)**2)
    weight = np.clip((distance-.65)/.4,0.,1.)
    chi = owner.seed*weight**2*(3-2*weight)
    chi[~owner.free] = owner.seed[~owner.free]
    registration = np.array([-np.log(shape.neck)])  # unit |dr/dS| directly at the initial rim
    return shape, chi, registration


def emit(**row):
    print(json.dumps(row,allow_nan=False),flush=True)


def snapshot(path,q,chi,registration):
    temp=path.with_suffix('.partial.npz')
    np.savez_compressed(temp,radius=4,shape_params=q,order=chi,registration=registration)
    temp.replace(path)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def order_matrices(wound,chi,amplitude,phase):
    """Assemble all dependence of the fixed-carrier Hamiltonian on nodal order."""
    rr,xw,ang,dx,dw,dt,_,_=wound.geometry[2]
    w=xw/wound.shape.length
    b,slope=wound.shape.radius(xw)
    nt,nr=np.subtract(chi.shape,1)
    i=np.minimum((rr*nr).astype(int),nr-1)
    j=np.minimum((w*nt).astype(int),nt-1)
    a,t=rr*nr-i,w*nt-j
    ids=(j*(nr+1)+i)[:,None]+np.array([0,1,nr+1,nr+2])
    nx=nr*np.stack([t-1,1-t,-t,t],-1)/b[:,None]
    nz=nt*np.stack([a-1,-a,1-a,a],-1)/wound.shape.length-(rr*slope)[:,None]*nx
    rows=np.repeat(np.arange(len(ids)),4)
    def coeff(g):
        return sparse.csr_matrix((g.ravel(),(rows,ids.ravel())),shape=(len(ids),chi.size))
    cx,cz=coeff(nx),coeff(nz)
    nodes=wound.cells[2][0]
    v2=amplitude[wound.did[nodes]].mean(axis=1)**2
    n=v2*wound.cells[2][3]
    field=phase[nodes]
    q,v=leggauss(3);q=(q+1)/2;v=v/2
    aa=np.zeros(len(ids));ab=aa.copy();bb=aa.copy()
    for a,wa in zip(q,v):
        x=rr+(a-.5)*dx
        for h,wh in zip(q,v):
            z=xw+(h-.5)*dw
            rad,slope=wound.shape.radius(z)
            r=rad*x
            for c,wc in zip(q,v):
                theta=ang+(c-.5)*dt
                f,_=wound.vortex(r*np.sin(theta))
                _,gx,gz,_=basis(a,h,c)
                pr=(field@gx)/(dx*rad)
                pw=(field@gz)/dw-x*slope*pr
                dv=4*np.pi*r*r*np.sin(theta)*rad*dx*dw*dt*wa*wh*wc*f*f
                aa+=dv*pr*pr
                ab+=dv*pr*pw
                bb+=dv*pw*pw

    def work(order):
        gr=cx@order.ravel();gw=cz@order.ravel()
        denom=n+ORDER_INERTIA*(gr*gr+gw*gw)
        r=aa*gr+ab*gw;z=ab*gr+bb*gw
        squared=r*gr+z*gw
        factor=-v2*ORDER_INERTIA/denom
        E=.5*np.sum(factor*squared)
        rem=ORDER_INERTIA*squared/denom
        g=cx.T@(factor*(r-rem*gr))+cz.T@(factor*(z-rem*gw))
        return E,g.reshape(order.shape)
    return work


def matter(shape,ref,co,n,receiver_n):
    t,x=np.linspace(0,1,2*n+1),np.linspace(0,1,n+1)
    elastic=source(shape,co,t,x,False)+core_source(shape,co,ref,co,t,x)
    outer=exterior(shape,receiver_n,CUT,side_order_work=lambda w,y:source(shape,co,w/shape.length,y,True))
    owner=Order(shape,n)
    return owner,elastic,outer['energy']+outer['static_work'],outer['residual']


def solve_order(owner,elastic,chi,flow):
    free=owner.free
    scale=np.sqrt(2*owner.stiff.diagonal()+2*owner.lump).reshape(chi.shape)[free]
    def objective(x):
        c=chi.copy()
        c[free]=x/scale
        en,g=owner.energy_gradient(c)
        phase,dp=flow(c)
        return en+np.sum(elastic*c)+phase,(g+elastic+dp)[free]/scale
    result=minimize(objective,chi[free]*scale,jac=True,method='L-BFGS-B',
                    bounds=list(zip(np.zeros_like(scale),scale)),
                    options=dict(gtol=3e-8,ftol=2e-15,maxiter=5000,maxls=60))
    station=np.max(abs(result.x-np.clip(result.x-result.jac,0,scale)))
    if station>2e-5:
        raise ArithmeticError(f'Unresolved static order {station} {result.message}')
    new=chi.copy()
    new[free]=result.x/scale
    return new,float(station)


class Relax:
    def __init__(self,n,mesh,polar,folder,flow_first=True):
        self.flow_first=flow_first
        self.ref,self.seed,self.co=preparation(n)
        self.n,self.mesh,self.polar=n,mesh,polar
        self.folder=folder
        self.q=self.ref.params.copy()
        self.chi=self.seed.copy()
        self.scale=np.array([1.,.25,.25,.25])
        self.cache=None

    def pack(self,q):
        return (q-self.ref.params)/self.scale

    def shape(self,x):
        return Meridian(4.,self.ref.params+x*self.scale)

    def model(self,shape,chi):
        return Wound(shape,chi,1,self.mesh,self.polar,28.,22.)

    def material(self,shape):
        return matter(shape,self.ref,self.co,self.n,self.n)

    def fields(self,shape,seed):
        owner,elastic,out,ores=self.material(shape)
        flow0=lambda c:(0.,np.zeros_like(c))
        chi=seed.copy() if self.flow_first else solve_order(owner,elastic,seed,flow0)[0]
        v=None
        for cycle in range(32):
            model=self.model(shape,chi)
            with contextlib.redirect_stdout(io.StringIO()):
                v,phase,info=model.solve(BIAS,v)
            flow=order_matrices(model,chi,v,phase)
            updated,station=solve_order(owner,elastic,chi,flow)
            work,g=flow(chi)
            E0,g0=owner.energy_gradient(chi)
            projected=(chi-np.clip(chi-(g0+g+elastic),0,1))[owner.free]
            err=float(max(abs(projected)))
            if err<4e-5:
                break
            chi=updated
        else:
            raise ArithmeticError(f'Unclosed order and carrier {err}')
        eg=E0+np.sum(elastic*chi)
        return (out+eg+info['excess_energy'],chi,v,phase,model,info,
                dict(material_inside=float(eg),material_outside=float(out),carrier=info['excess_energy'],
                     order_gradient=err,receiver_residual=ores,current=info['total_current'],order_cycles=cycle+1))

    def frozen(self,x,chi,v,p):
        sh=self.shape(x)
        owner,elastic,out,_=self.material(sh)
        model=self.model(sh,chi)
        return owner.energy_gradient(chi)[0]+np.sum(chi*elastic)+out+model.amplitude(v,p)[0]

    def eval(self,x):
        if self.cache is not None and np.array_equal(x,self.cache[0]):
            return self.cache[1],self.cache[2]
        sh=self.shape(x)
        if sh.wall()['minimum_radius'] <= .15 or sh.length>8 or sh.exit>8:
            return 1e8+(x@x),2*x
        E,chi,v,p,model,info,stats=self.fields(sh,self.chi)
        h=.0003
        # Distinct shape changes re-integrate both material sides and the whole
        # two-sided carrier. Hold stationary carrier coordinates during virtual work.
        I=np.eye(len(x))
        grad=np.array([(self.frozen(x+h*i,chi,v,p)-self.frozen(x-h*i,chi,v,p))/(2*h) for i in I])
        self.cache=(x.copy(),E,grad,chi,stats)
        return E,grad

    def record(self,x):
        E,g=self.eval(x)
        _,_,_,self.chi,info=self.cache
        sh=self.shape(x)
        self.q=sh.params.copy()
        tag=snapshot(self.folder/'state.npz',self.q,self.chi,self.co)
        init_pocket=(self.seed<.01)
        emit(kind='accepted',energy=E,shape_params=self.q.tolist(),depth=sh.length,exit=sh.exit,
             waist=sh.wall()['minimum_radius'],shape_forces=(-g).tolist(),
             axial_outward=-g[0]/sh.length,exit_outward=-g[1]/(sh.exit*.25),
             pocket_peak_order=float(self.chi[init_pocket].max()),
             old_pocket_mean_order=float(self.chi[init_pocket].mean()),
             shaft_mid_order=float(self.chi[self.n,0]),upper_collar=collar(sh,self.co,False),
             lower_collar=collar(sh,self.co,True),state_sha256=tag,**info)

    def run(self):
        x=self.pack(self.q)
        self.record(x)
        self.cache=None
        out=minimize(self.eval,x,jac=True,method='L-BFGS-B',callback=self.record,
                     bounds=[(-3.8,1.6),(-7.,2.7),(-2.,2.),(-2.,2.)],
                     options=dict(gtol=.0005,ftol=2e-15,maxiter=350,maxls=28))
        self.record(out.x)
        lo,hi=np.array([[-3.8,1.6],[-7.,2.7],[-2.,2.],[-2.,2.]]).T
        at=(out.x-lo<1.e-5)|(hi-out.x<1.e-5)
        emit(kind='optimization',shape_residual=float(max(abs(out.jac))),message=str(out.message),
             active_numerical_limits=at.tolist(),opt_success=bool(out.success))
        if max(abs(out.jac))>.01 or np.any(at):
            emit(kind='unresolved',reason='shape gradient or coordinate-domain boundary')
        else:
            emit(kind='conditional_stationary',scope='Finite global shapes and registered material, sharp outer reservoir matching')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    parser.add_argument('--order-mesh',type=int,default=32)
    parser.add_argument('--carrier-mesh',type=int,default=7)
    parser.add_argument('--polar',type=int,default=6)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--condition-without-flow',action='store_true',
                        help='Control that conditions the initial order before adding the through-flow')
    args=parser.parse_args()
    if args.prepare_only:
        sh,chi,co=preparation(args.order_mesh)
        sha=snapshot(args.folder/'state.npz',sh.params,chi,co)
        emit(kind='preparation',radius=4,depth=sh.length,exit=sh.exit,waist=sh.wall()['minimum_radius'],
             sha256=sha,material_order_mesh=args.order_mesh)
    else:
        Relax(args.order_mesh,args.carrier_mesh,args.polar,args.folder,not args.condition_without_flow).run()
