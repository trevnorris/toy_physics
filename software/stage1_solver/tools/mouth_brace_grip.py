"""Probe a passing neutral bilayer when conversion removes only its normal bracing."""
import argparse
import ctypes as ct
from pathlib import Path
import subprocess

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh

from mouth_middle_tag_grip import cubic_trace
from neutral_patch_open_w import OpenCell
from neutral_patch_sites import BETWEEN, GAP, TEMPERATURE, WITHIN, emit
from neutral_patch_structural_endpoints import preparation


def clock(mesh):
    trace,intake=cubic_trace(mesh)
    depths=np.array([trace.edges[0],1.175,trace.crossing(.5),trace.crossing(.1),trace.crossing(.01),1.6,2.])
    edges=np.unique(np.r_[depths,trace.edges[(trace.edges>depths[0])&(trace.edges<depths[-1])]])
    z=np.r_[np.concatenate([np.linspace(a,b,16,endpoint=False) for a,b in zip(edges[:-1],edges[1:])]),edges[-1]]
    g,w=leggauss(6)
    half=np.diff(z)[:,None]/2
    t=np.r_[0.,np.cumsum((half*(1/trace.speed(z[:-1,None]+half*(1+g)))@w))]
    tau=t[np.searchsorted(z,depths)]
    np.testing.assert_allclose(tau[-1],trace.transit(depths[0],depths[-1],8),rtol=1.e-9)
    chi=np.clip(trace.sample(trace.chi,z)[0],0.,1.)
    return depths,tau,PchipInterpolator(t,chi),intake


class Material:
    def __init__(self,args):
        self.cell=OpenCell(args.side,args.cutoff,args.power)
        self.n=self.cell.n
        self.half=self.n//2
        self.sign=np.repeat([-.5,.5],self.half)
        self.flat=np.column_stack([preparation(self.cell,0),self.sign*GAP])
        self.load=np.ascontiguousarray(self.cell.evaluate(self.flat)[1][:,2])
        self.height=self.flat[:,2].copy()
        self.brace=2*WITHIN[0]*WITHIN[1]**2 if args.brace=='elastic' else 0.
        self.mobility=args.mobility
        self.power=args.power
        self.temperature=args.temperature
        assert GAP==1.
        image_cell=OpenCell(args.side,args.cutoff+.36,args.power)
        self.indices=[np.ascontiguousarray(x,dtype=np.int32) for x in (image_cell.i,image_cell.j)]
        self.shifts=np.ascontiguousarray(image_cell.shifts,dtype=float)
        self.wells=[np.array(p,dtype=float) for p in (WITHIN,BETWEEN)]
        src=Path(__file__).with_suffix('.cpp')
        library=args.work/'brace.so'
        subprocess.run(['g++','-O3','-std=c++17','-fPIC','-shared',str(src),'-o',str(library)],check=True)
        self.advance=ct.CDLL(str(library.resolve())).evolve
        f64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        i32=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.advance.argtypes=[ct.c_int,ct.c_int,i32,i32,f64,f64,ct.c_double,ct.c_int,
                               f64,f64,f64,f64,ct.c_double,ct.c_double,ct.c_int,
                               ct.c_int,ct.c_double,f64,f64,f64]
        self.advance.restype=ct.c_int
        self.cutoff=args.cutoff

    def step(self,pos,dt,chi,noise,hold=False):
        code=self.advance(self.n,len(self.indices[0]),*self.indices,self.shifts,self.cell.box,
                          self.cutoff,self.power,*self.wells,self.height,self.load,self.brace,
                          self.mobility,int(hold),len(chi),dt,chi,noise,pos)
        if code:
            raise ArithmeticError(f'Contact path failed native step check {code}')

    def slope(self,pos,check=True):
        directions=np.zeros((3*self.n,2))
        directions[0::3,0]=self.sign
        directions[1::3,1]=self.sign
        h=self.cell.evaluate(pos,tangent=True)[0]
        k=directions.T@h@directions/self.cell.area
        if check:
            for delta in (1.e-4,.5e-4):
                for axis in range(2):
                    move=self.sign[:,None]*np.eye(3)[axis]
                    hi=self.cell.evaluate(pos+delta*move)[1]
                    lo=self.cell.evaluate(pos-delta*move)[1]
                    measured=self.sign@(hi[:,:2]-lo[:,:2])/(2*delta*self.cell.area)
                    np.testing.assert_allclose(measured,k[:,axis],atol=8e-5,rtol=3e-5)
        return k

    def validate(self):
        k0=self.slope(self.flat)
        np.testing.assert_allclose(k0,k0[0,0]*np.eye(2),atol=1.e-10)
        rng=np.random.default_rng(705)
        trial=self.flat+.018*rng.normal(size=self.flat.shape)
        _,grad,_=self.cell.evaluate(trial)
        expected=-grad
        expected[:,2]+=.45*(self.load-self.brace*(trial[:,2]-self.height))
        observed=trial.copy()
        self.step(observed,1e-5,np.array([.45]),np.zeros((1,self.n,3)))
        np.testing.assert_allclose((observed-trial)/1e-5,self.mobility*expected,atol=1e-9)
        h=self.cell.evaluate(self.flat,tangent=True)[0]
        unbraced=eigh(h,eigvals_only=True)[0]
        with_brace=h.copy()
        with_brace[2::3,2::3]+=self.brace*np.eye(self.n)
        emit(test='ordered_reference',sites=self.n,ordered_x_grip=k0[0,0],pressure=self.load[self.half:].sum()/self.cell.area,
             normal_spring=self.brace,least_at_dead_load=unbraced,least_with_chosen_brace=eigh(with_brace,eigvals_only=True)[0])
        return k0[0,0]


def evolve_interval(mat,states,t0,t1,profile,rng,bridge,refine,hold=False):
    # Identical Brownian increments over each pair of refined steps. The
    # material and matched control see the same noise and physical times.
    coarse_steps=int(np.ceil((t1-t0)/.001))
    coarse_dt=(t1-t0)/coarse_steps
    for i in range(0,coarse_steps,512):
        count=min(512,coarse_steps-i)
        delta=rng.normal(size=(count,mat.n,3))*np.sqrt(2*mat.temperature*mat.mobility*coarse_dt)
        if refine==2:
            extra=bridge.normal(size=delta.shape)*np.sqrt(2*mat.temperature*mat.mobility*coarse_dt)
            delta=np.stack(((delta+extra)/2,(delta-extra)/2),axis=1).reshape(2*count,mat.n,3)
        dt=coarse_dt/refine
        times=t0+i*coarse_dt+dt*(.5+np.arange(len(delta)))
        for pos,schedule in zip(states,profile,strict=True):
            chi=np.ascontiguousarray(schedule(times))
            mat.step(pos,dt,chi,delta,hold)


def held_probe(mat,initial,t0,order,k0,args,name):
    # Hold both original groups' x centres, allowing all other site motions.
    # Matched copies evolve under fresh, shared noise and the same order history.
    amps=np.array([0.,1.e-4,-1.e-4,.5e-4,-.5e-4])
    states=[initial.copy() for _ in amps]
    for pos,amp in zip(states,amps):
        pos[:,0]+=amp*mat.sign
    rel0=np.mean(initial[mat.half:,0])-np.mean(initial[:mat.half,0])
    rng=np.random.default_rng(args.seed+404)
    bridge=np.random.default_rng(args.seed+7404)
    old=t0
    for elapsed in (0.,.125,.5,1.,2.):
        t=t0+elapsed
        if elapsed:
            evolve_interval(mat,states,old,t,[order]*len(states),rng,bridge,args.refine,True)
        readings=np.array([mat.sign@mat.cell.evaluate(p)[1][:,0]/mat.cell.area for p in states])
        full=(readings[1]-readings[2])/(amps[1]-amps[2])
        half=(readings[3]-readings[4])/(amps[3]-amps[4])
        errors=[p[mat.half:,0].mean()-p[:mat.half,0].mean()-rel0-a for p,a in zip(states,amps)]
        emit(test='held_probe',history=name,elapsed=elapsed,order=float(order(np.array([t]))[0]),
             x_ratio_full=full/k0,x_ratio_half=half/k0,max_centroid_error=max(abs(np.asarray(errors))),
             signed_holding_tractions=readings.tolist())
        if max(abs(np.asarray(errors)))>2e-9 or abs(full-half)>1.e-3*max(k0,abs(half)):
            raise ArithmeticError('Held relative displacement or linear probe failed')
        old=t


def endpoint(mat,pos,k0):
    # A local cold basin tests whether releasing the brace has left a solid.
    # Cell shear and an opposite slip of the two original groups are distinct.
    try:
        eq,_,stress,res=mat.cell.relax(pos)
        h,b,born=mat.cell.evaluate(eq,tangent=True)
        eig,vec=eigh(h)
        if max(abs(eig[:3]))>1.e-6 or eig[3]<1.e-7:
            emit(test='cold_endpoint',resolved=False,reason='Unresolved local stability',
                 lowest_eigenvalues=eig[:5].tolist(),force_residual=res)
            return
        a=np.zeros(3*mat.n)
        a[::3]=2*mat.sign/mat.half
        soft=(vec[:,3:].T@a)
        stiffness=1/(mat.cell.area*np.sum(soft**2/eig[3:]))
        shear=(born-np.sum((vec[:,3:].T@b)**2/eig[3:]))/mat.cell.area
        emit(test='cold_endpoint',resolved=True,force_residual=res,min_internal=eig[3],
             relaxed_tagged_grip_ratio=stiffness/k0,relaxed_common_cell_shear=shear,stress=stress,
             original_group_gap=eq[mat.half:,2].mean()-eq[:mat.half,2].mean())
    except ArithmeticError as exc:
        emit(test='cold_endpoint',resolved=False,reason=str(exc))


def run(args):
    args.work.mkdir(parents=True,exist_ok=True)
    depth,times,chi,intake=clock(args.mesh)
    mat=Material(args)
    k0=mat.validate()
    unit=np.ones_like
    rng=np.random.default_rng(args.seed)
    bridge=np.random.default_rng(args.seed+7013)
    initial=mat.flat.copy()
    perturb=np.random.default_rng(args.seed+11).normal(size=initial.shape)
    initial+=.001*(perturb-perturb.mean(axis=0))
    states=[initial.copy(),initial.copy()]
    emit(test='setup',side=args.side,power=args.power,seed=args.seed,site_mobility=args.mobility,
         temperature=args.temperature,cutoff=args.cutoff,intake=intake,
         depths=depth.tolist(),times=times.tolist(),refine=args.refine,brace=args.brace)
    for j in range(len(depth)):
        if j:
            evolve_interval(mat,states,times[j-1],times[j],[chi,unit],rng,bridge,args.refine)
        for name,p in zip(('converting','full_brace'),states):
            k=mat.slope(p)
            z=p[:,2]
            gap=z[mat.half:].mean()-z[:mat.half].mean()
            rough=[float(np.std(s)) for s in np.split(z,2)]
            order=float(chi(times[j])) if name=='converting' else 1.
            emit(test='passage',history=name,depth=depth[j],time=times[j],order=order,
                 gap=gap,roughness=rough,mean_grip_ratio=np.trace(k)/(2*k0),grip_tensor_over_ordered=(k/k0).tolist(),
                 eigen_grip_over_ordered=np.linalg.eigvalsh(k/k0).tolist(),
                 label_height_overlap=max(z[:mat.half])-min(z[mat.half:]))
        if j==4:
            held_probe(mat,states[0],times[j],chi,k0,args,'converting')
            held_probe(mat,states[1],times[j],unit,k0,args,'full_brace')
    endpoint(mat,states[0],k0)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34),default=34)
    p.add_argument('--side',type=int,choices=(4,6),default=4)
    p.add_argument('--power',type=int,choices=(1,2),default=1)
    p.add_argument('--seed',type=int,default=319)
    p.add_argument('--refine',type=int,choices=(1,2),default=1)
    p.add_argument('--cutoff',type=float,default=3.5)
    p.add_argument('--temperature',type=float,default=TEMPERATURE)
    p.add_argument('--mobility',type=float,default=1.)
    p.add_argument('--brace',choices=('elastic','dead'),default='elastic')
    p.add_argument('--work',type=Path,required=True)
    run(p.parse_args())
