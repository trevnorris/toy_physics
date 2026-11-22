"""Two-sided trial relaxation from unit depth on the separately registered pocket."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss, legval, legder

from mouth_exchange_carrier import BIAS
from mouth_finite_bore_support import graph_values, quadrature
from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_fields import CoupledLoads
from mouth_live_waist_order import Order
from mouth_live_waist_reopen import reopen
from mouth_polycrystal_material import stable_energy
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import PORT, ROOT, read_case
from mouth_strained_wave import install
from mouth_waist_outer_order import exterior


CUT = 12.  # Common exterior box contains the outer elastic collar on every probe.


def emit(**items):
    print(json.dumps(items, allow_nan=False), flush=True)


def source(shape, c, t, x, outside):
    """Integrate reference W times Eulerian order shape functions on exact split cells."""
    grid = np.zeros((len(t), len(x)))
    dx = np.diff(x)[None,:]
    dt = np.diff(t)[:,None]
    points, weights = leggauss(4)
    for z, wz in zip((1+points)/2, weights/2):
        tq = t[:-1,None] + dt*z
        b, slope, _ = shape.profile(tq)
        lnp = legval(2*tq-1,c)
        p = np.exp(lnp)
        pt = 2*legval(2*tq-1,legder(c))
        if outside:
            edge = b*np.expm1(.5*p)/(CUT-b)
            if np.max(edge)>=1.:
                raise ArithmeticError("Receiving collar outside the order box")
            lo, hi = x[None,:-1], np.minimum(edge,x[None,1:])
            width, origin = CUT-b, b
        else:
            lo, hi = np.maximum(np.exp(-.5*p),x[None,:-1]), x[None,1:]
            width, origin = b, 0.
        span = np.maximum(hi-lo,0.)
        for y, wy in zip((1+points)/2, weights/2):
            xx = lo+span*y
            r = origin+xx*width
            a = -np.log(r/b)/p
            # Inactive quadrature points can lie beyond the specimen.
            a = np.clip(a,-.5,.5)
            _,rs,rz,h,J,_ = graph_values(a,tq,shape,lnp,pt)
            tr = rs*rs+rz*rz+shape.length**2+2*h*h
            density = stable_energy(tr,J)/J
            measure = (4*np.pi*r*r*width*span*shape.length*dt*wz*wy)*density
            u = (xx-x[None,:-1])/dx
            for iz,vz in ((0,1-z),(1,z)):
                for ix,vx in ((0,1-u),(1,u)):
                    grid[iz:iz+len(t)-1, ix:ix+len(x)-1] += measure*vz*vx
    # This separate integration in the reference slab checks both Jacobians.
    a,ts,w=quadrature(-1 if outside else 1,120,40)
    _,rs,rz,h,j,_=graph_values(a,ts,shape,legval(2*ts-1,c),2*legval(2*ts-1,legder(c)))
    ref=np.sum(4*np.pi*(shape.neck+a)**2*w*stable_energy(rs*rs+rz*rz+shape.length**2+2*h*h,j))
    if abs(grid.sum()-ref)>3e-4*max(1,ref) or np.min(grid)<-1e-12:
        raise ArithmeticError(f"Order/elastic quadrature mismatch {grid.sum()} vs {ref}")
    return grid


def get_input(path):
    guard=json.loads((path.parent/"outcome.json").read_text())
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    if (not guard["limitsVerified"] or guard["exitCode"] or guard["childOutcome"]["exitCode"]
            or sha not in (path.parent/"stdout").read_text()):
        raise ValueError(f"Unverified starting field {path}")
    return sha


class Trial:
    def __init__(self, radius, n, *, bonded_core=False):
        self.radius, self.n = radius, n
        prep,*_ = read_case(radius,n)
        self.mesh=prep["field_mesh"]
        kind="bounded" if radius==1 else "relaxed"
        fp=ROOT/f"_scratch/s11c/shallow_puncture_fields01/{kind}_r{radius}_m{n}/center.npz"
        gp=ROOT/f"_scratch/s11c/polycrystal_material01/r{radius}g2/coeffs.npz"
        self.fingerprints=[get_input(p) for p in (fp,gp)]
        with np.load(fp) as f, np.load(gp) as g:
            self.chi=f["order"]; self.x=f["coordinates"]
            self.co, self.cr = g["ordered"], g["receiver"]
            if (abs(Meridian(radius,self.x).length-1)>1e-12
                    or np.max(abs(self.x-g["coordinates"]))>1e-12):
                raise ValueError("Geometry mismatch at the shallow starting state")
        area=ROOT/"_scratch/s11c/poly_wave01"
        key=f"r{radius}n{n}q3"
        saved=area/"fields"/f"selected_{key}.npz"
        run=area/f"release_{key}"
        guard=json.loads((run/"outcome.json").read_text())
        record=json.loads((run/"stdout").read_text().splitlines()[0])
        sha=hashlib.sha256(saved.read_bytes()).hexdigest()
        if (guard["exitCode"] or not guard["limitsVerified"] or guard["childOutcome"]["exitCode"]
                or record["sha256"]!=sha or record["graph_sha256"]!=self.fingerprints[1]):
            raise ValueError("Stabilized material pole was not validated")
        with np.load(saved) as w:
            if not np.array_equal(w["order"],self.chi) or not np.array_equal(w["coeffs"],self.co):
                raise ValueError("Order or lattice changed at the continued wave")
            self.seed=complex(w["pole"].item()),w["vector"]
        self.fingerprints.append(sha)
        self.core_reference = (Meridian(radius,self.x), self.co) if bonded_core else None

    def evaluate(self, shifts, chi, seed=None):
        x=self.x + shifts[:4]
        shape=Meridian(self.radius,x)
        co=self.co.copy();co[0]+=shifts[4]
        cr=self.cr.copy();cr[0]+=shifts[5]
        t=np.linspace(0.,1.,2*self.n+1); s=np.linspace(0.,1.,self.n+1)
        load=source(shape,co,t,s,False)
        core_work = 0.
        if self.core_reference is not None:
            from mouth_poly_bonded_core import core_source
            core = core_source(shape,co,*self.core_reference,t,s)
            core_work = float(np.sum(core*chi))
            load += core
        def receiving(wa,xo):
            return source(shape,cr,wa/shape.length,xo,True)
        outer=exterior(shape,self.n,CUT,side_order_work=receiving)
        owner=Order(shape,self.n)
        material=CoupledLoads(shape,self.n,2*self.n,self.mesh,PORT,
                              self.seed if seed is None else seed)
        install(material.wave,shape,co,core_reference=self.core_reference)
        ActiveWave.install(material)
        material.wave.target=material.wave.seed[0].real
        ei,gi=owner.energy_gradient(chi)
        ef,gf,info=material.evaluate(chi,1.,BIAS)
        terms=dict(inner_order=ei,outer_order=outer["energy"],
                   static_inside=float(np.sum(load*chi)),static_outside=outer["static_work"],
                   two_sided_carrier=info["carrier_energy"],strained_wave=info["wave_energy"])
        return (sum(terms.values()),gi+gf+load,owner,material,
                dict(energy=terms,core_elastic_portion=core_work,receiver_residual=outer["residual"],
                     current=info["current"],outflow=info["terminal_current"],
                     depth=shape.length,waist=shape.wall()["minimum_radius"],
                     mode_overlap=info["mode_overlap"],wave_components=info["wave_components"],
                     cavity_fraction=info["wave_pocket_fraction"],pole_residual=info["pole_residual"]),
                load)


def opening_sites(test,chi,order,material,elastic):
    wave=material.wave
    om,v=wave.seed
    proj=wave.project.tocsc()
    weights=proj@chi.T.ravel()
    mass=wave.cut_m+wave.matrix(weights[:,None,None]*wave.mc)
    live=mass.diagonal()>0
    diff=(wave.port.block(om+2e-6)-wave.port.block(om-2e-6))/(4e-6)
    edge=v[wave.port.ids]
    den=2*om*(v@(mass@v))+edge@(diff@edge)
    q=wave.kc-om**2*wave.mc
    eo,go=order.energy_gradient(chi)
    ec,gc,_=material.carrier.response(chi,BIAS)
    candidates=[]
    for i,j in np.argwhere((chi==0)&order.free):
        k=j*(2*test.n+1)+i
        sl=slice(proj.indptr[k],proj.indptr[k+1])
        direct,correction,new,extension,_=reopen(proj.indices[sl],proj.data[sl],wave,q,v,live)
        slope=go[i,j]+gc[i,j]+elastic[i,j]+((direct+correction)/den).real
        candidates.append((slope,i,j,new,extension))
    if not candidates:
        return
    slope,i,j,new,extension=min(candidates,key=lambda row:row[0])
    v1=v.copy();v1[new]=extension
    checks=[]
    for h in (.001,.0005):
        trial=chi.copy();trial[i,j]+=h
        wave.seed=om,v1.copy();wave.target=om.real
        e,_,_=material.evaluate(trial,1.,BIAS)
        checks.append((order.energy_gradient(trial)[0]+e-eo-ec-om.real)/h+elastic[i,j])
    wave.seed=om,v;wave.target=om.real
    emit(kind="previously_melted",tested=len(candidates),min_slope=slope,
         negative=sum(int(row[0]<-1e-6) for row in candidates),index=[int(i),int(j)],
         fractions=[i/(2*test.n),j/test.n],direct_secants=checks)
    if abs(2*checks[1]-checks[0]-slope)>2e-5*max(1,abs(slope)):
        raise ArithmeticError("Single-site one-sided derivative unresolved")


def run(radius,n):
    test=Trial(radius,n)
    shifts=np.zeros(6); chi=test.chi
    e0,g,order,mat,initial,elastic=test.evaluate(shifts,chi)
    emit(kind="initial",radius=radius,mesh=n,source_sha256=test.fingerprints,**initial)
    opening_sites(test,chi,order,mat,elastic)
    seed=mat.wave.seed
    # A simultaneous admissible release. Restore neither a lost Jacobian nor an
    # unsigned mode crossing by assigning an arbitrary force or a new eigenmode.
    h=.0005; grad=[]
    for k in range(len(shifts)):
        offset=shifts.copy();offset[k]=h
        plus=test.evaluate(offset,chi,seed)[0]
        minus=test.evaluate(-offset,chi,seed)[0]
        grad.append((plus-minus)/(2*h))
    # Finite opening/length/shape variations vanish at the specified entrance.
    # Scale by reference units before taking a descent direction.
    gradient=np.array(grad)
    structural=-gradient/max(np.max(abs(gradient)),1e-12)
    scale=(2*order.stiff.diagonal()+1.5*order.lump).reshape(chi.shape)
    dchi=np.zeros_like(chi)
    live=order.free&(chi>0)&(chi<1)
    dchi[live]=-g[live]/scale[live]*chi[live]*(1-chi[live])
    dchi/=max(np.max(abs(dchi)),1e-12)
    pred=gradient@structural+np.sum(g*dchi)
    if pred>=0: raise ArithmeticError("No descent for the joint trial")
    emit(kind="direction",shape_and_spacing_axes=["log_depth","log_exit_radius","middle",
         "upper_minus_lower","inner_log_spacing","outer_log_spacing"],
         slopes=grad,direction=structural.tolist(),order_slope=np.sum(g*dchi),
         predicted=pred)
    for step in (.002,.001):
        cnew=chi+step*dchi
        if min(cnew.ravel())<0 or max(cnew.ravel())>1:
            raise ArithmeticError("Order trial crossed active bounds")
        e,_,_,_,out,_=test.evaluate(step*structural,cnew,seed)
        emit(kind="joint_step",step=step,change=e-e0,secant=(e-e0)/step,
             order_shift=float(np.max(abs(cnew-chi))),**out)
        if e>=e0 or out["mode_overlap"]<.95:
            raise ArithmeticError("Joint release failed its energy or wave test")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius",choices=(1,2),type=int,required=True)
    parser.add_argument("--mesh",choices=(48,72),type=int,required=True)
    args=parser.parse_args()
    run(args.radius,args.mesh)
