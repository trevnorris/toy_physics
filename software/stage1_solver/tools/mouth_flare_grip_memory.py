"""Does grip that survives conversion allow the independently calibrated coating to hold the flare?

Only the strength of a local intersheet slip spring travels with a material
label. Tangential displacements relax elastically on the frozen strength field.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.optimize import brentq

from mouth_calibrated_flare_coating import read_inputs
from mouth_flare_support_bounds import minimum
from mouth_flared_receiver import emit
from mouth_middle_tag_grip import cubic_trace
from mouth_replacing_flare_slip import assemble


class PassingGrip:
    def __init__(self, trace, points, grip, subdivisions):
        roots, weights = leggauss(6)
        width = np.diff(points)[:,None]
        z = points[:-1,None]+width*(1+roots)/2
        radius,slope = trace.geometry.radius(z)
        self.weight=(grip*4*np.pi/3*radius**2*np.hypot(1,slope)*width*weights/2).ravel()
        n = len(points)
        # Same physical two-sheet energy as the archived interlayer matrix.
        # Remove the four upstream slips and leave the flare coordinate last.
        i = np.repeat(np.arange(n-1),len(roots))
        columns = np.column_stack((4*i,4*i+4,4*i+2,4*i+6))
        basis = np.column_stack(((1-roots)/2,(1+roots)/2,-(1-roots)/2,-(1+roots)/2))
        values=np.tile(basis,(n-1,1)).ravel()
        row=np.repeat(np.arange(z.size),4)
        self.slip=[]
        for offset in (0,1):
            matrix=sparse.csr_matrix((values,(row,(columns+offset).ravel())),shape=(z.size,4*n+1))
            self.slip.append(matrix[:,4:])
        self.location=trace.crossing(.01)
        samples=np.concatenate([np.linspace(a,b,subdivisions,endpoint=False)
                                for a,b in zip(points[:-1],points[1:])])
        grid=np.unique(np.r_[samples,points,z.ravel(),1.175,self.location])
        self.quad_idx=np.searchsorted(grid,z.ravel())
        self.probe_idx=np.searchsorted(grid,[1.175,self.location])
        self.order=np.clip(trace.sample(trace.chi,grid)[0],0.,1.)
        mid=(grid[1:]+grid[:-1])/2
        half=np.diff(grid)/2
        self.dt=half*((1/trace.speed(mid[:,None]+half[:,None]*roots))@weights)
        self.time=np.r_[0.,np.cumsum(self.dt)]
        self.transit=self.time[-1]

    def strength(self,tau):
        if tau==0:
            return self.order
        if np.isinf(tau):
            return np.full_like(self.order,self.order[0])
        # Exact update for a linearly varying equilibrium strength in material time.
        u=self.dt/tau
        change=-np.expm1(-u)
        lag_weight=np.empty_like(u)
        small=u<1.e-4
        lag_weight[small]=u[small]*(.5+u[small]*(-1/6+u[small]/24))
        lag_weight[~small]=1-change[~small]/u[~small]
        increments=lag_weight*np.diff(self.order)
        m=np.empty_like(self.order)
        m[0]=self.order[0]
        for i in range(len(u)):
            m[i+1]=m[i]+change[i]*(self.order[i]-m[i])+increments[i]
        return m

    def coupling(self,m):
        w=self.weight*m[self.quad_idx]
        a,b=self.slip
        return (a.T @ a.multiply(w[:,None])+2*b.T @ b.multiply(w[:,None])).tocsc()


def run(args):
    trace,intake=cubic_trace(args.mesh)
    patch,base,areas=read_inputs(args.mesh)
    points,grip,(membrane,old_between,wall,_,_)=assemble(trace,patch,args.refine)
    process=PassingGrip(trace,points,grip,args.subdivisions)
    zero=process.coupling(process.strength(0))-old_between
    agreement=np.max(abs(zero.data),initial=0.)
    bare=float(base[2,2])
    target=bare+.02*areas[2]
    slips=np.arange(membrane.shape[0]-1)
    emit(test='preparation',mesh=args.mesh,refine=args.refine,subdivisions=args.subdivisions,intake=intake,
         duration=process.transit,probe_depths=[1.175,process.location],probe_times=process.time[process.probe_idx].tolist(),
         start_order=process.order[0],instantaneous_matrix_error=agreement,target=target,
         physical_observable='intersheet differential tangential stiffness for a fresh displacement on a passing sample')
    if agreement>1.e-9 or abs(intake-.25)>1.e-8:
        raise ArithmeticError('Instantaneous interlayer or intake changed')
    for name,free,extra in [('face_locked',slips[slips%4>=2],0*wall),('neutral_registry',slips,wall)]:
        def measure(tau):
            m=process.strength(tau)
            _,force,energy,residual=minimum(membrane+extra+process.coupling(m),free)
            np.testing.assert_allclose(force,energy,rtol=1.e-8)
            if residual>1.e-8:
                raise ArithmeticError('Unbalanced sheet')
            return bare+force,m
        times=[0.,.25,1.,4.,16.,64.,np.inf]
        readings=[]
        for tau in times:
            total,m=measure(tau)
            readings.append(total)
            emit(test='memory',wall=name,time_constant=None if np.isinf(tau) else tau,
                 infinite_memory=bool(np.isinf(tau)),total=float(total),
                 local_relative_to_ordered=m[process.probe_idx].tolist())
        if min(np.diff(readings)) < -1.e-8:
            raise ArithmeticError('Nonmonotonic sampled response; inspect all crossings')
        finite_pairs=zip(times[:-2],times[1:-1],readings[:-2],readings[1:-1])
        for lo,hi,left,right in finite_pairs:
            if left<target<=right:
                tau=brentq(lambda t:measure(t)[0]-target,lo,hi,xtol=2.e-9)
                total,m=measure(tau)
                middle,low=m[process.probe_idx]
                emit(test='weak_target_crossing',wall=name,time_constant=tau,bracket=[lo,hi],total=total,
                     relative_to_ordered_at_middle=middle,relative_to_ordered_at_one_percent=low,
                     one_percent_relative_to_middle=low/middle,one_percent_order=process.order[process.probe_idx[1]])
                break


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34),required=True)
    p.add_argument('--refine',type=int,choices=(1,2,4),required=True)
    p.add_argument('--subdivisions',type=int,choices=(8,16,32),required=True)
    run(p.parse_args())
