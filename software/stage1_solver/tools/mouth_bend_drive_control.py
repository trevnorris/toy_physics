"""Finite-offset, two-reservoir carrier controls, including momentum on two internal cuts."""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_response import prepare, trace
from mouth_empty_bend_offsets import BentCarrier
from mouth_empty_bend_stiffness import shapes
from mouth_flared_receiver import emit


class DriveControl(BentCarrier):
    def evaluate(self, state, offset):
        self.last_field=state
        return super().evaluate(state, offset)

    def fixed_phase(self, offset):
        field=np.zeros((len(self.base),self.modes),complex)
        field[:,0]=self.base
        for iteration in range(45):
            residual,energy,force=self.evaluate(field,offset)
            error=0.
            updates=np.zeros_like(field)
            for l,(factor,mask) in enumerate(self.preconditioners):
                rhs=np.r_[residual[:,l].real,residual[:,l].imag]
                error=max(error,np.max(abs(rhs[mask])))
                step=np.zeros_like(rhs)
                step[mask]=factor.solve(-rhs[mask])
                real,imag=np.split(step,2)
                updates[:,l]=real+1j*imag
            if error<3e-11:
                flux=2*(field[:,0].conj()*residual[:,0]).imag
                return dict(offset=list(offset),phase_bias=np.angle(field[self.dom.bulk,0][0]),
                            feed=-flux[self.dom.feed].sum(),bulk=flux[self.dom.bulk].sum(),
                            energy=energy,holding_force=force.tolist(),residual=error,iterations=iteration)
            field+=updates
        raise ArithmeticError(f'Fixed-phase carrier did not settle, residual={error:.4g}')


def line_rule(edges):
    root,weight=leggauss(8)
    half=np.diff(edges)[:,None]/2
    return (edges[:-1,None]+half*(root+1)).ravel(),(half*weight).ravel()


def flux(solver, offset, z, v, weights, side):
    dom=solver.dom
    r,w,_,_,_,_,f,fr,fw=trace(dom,2,z,v,solver.last_field)
    radius,slope=dom.radius(w)
    _,sr,sw=shapes(r,w,'compact',dom.far,dom.roof)
    if np.any(sr):
        raise ValueError('A cut crosses the transverse cutoff of the shape map')
    bend=sw@offset
    measure=4*np.pi*r*r*weights*(1 if side else radius)
    answer=np.zeros(5)
    for mu,aw,y,dy in zip(solver.mu,solver.weights,solver.y,solver.dy,strict=True):
        val=f@y
        radial=fr@y
        angular=f@dy/r
        lateral=mu*radial+(1-mu**2)*angular
        axial=fw@y-bend*lateral
        n=abs(val)**2
        energy=.5*(abs(radial)**2+(1-mu**2)*abs(angular)**2+abs(axial)**2)
        energy+=(n-1)**2*(n**3+2*n*n+3*n+4)/20
        normal=radial-(slope+bend*mu)*axial if side else axial
        jn=(val.conj()*normal).imag
        jx=(val.conj()*lateral).imag
        advective=jx*jn/n
        momentum=(lateral.conj()*normal).real-(mu*energy if side else 0.)
        # Explicit time harmonic gauge shift of one model energy unit per
        # constituent tests the net energy budget without a vanishing gauge current.
        energy_current=-((-1j*val).conj()*normal).real
        pressure_only=-mu*energy if side else np.zeros_like(w)
        answer+=aw*np.array([measure@a for a in (jn,advective,momentum,energy_current,pressure_only)])
    return answer


def cuts(solver,offset,lower,upper):
    depths,radial,_=solver.dom.blocks[2]
    v,vwt=line_rule(radial)
    incoming=flux(solver,offset,np.full_like(v,lower),v,vwt,False)
    outgoing=flux(solver,offset,np.full_like(v,upper),v,vwt,False)
    z,zwt=line_rule(np.r_[lower,depths[(depths>lower)&(depths<upper)],upper])
    wall=flux(solver,offset,z,np.zeros_like(z),zwt,True)
    return dict(lower_depth=lower,upper_depth=upper,low=incoming.tolist(),high=outgoing.tolist(),
                side_out=wall.tolist(),closure=(outgoing-incoming+wall).tolist())


def main(a):
    dom,field=prepare(a)
    solver=DriveControl(dom,field,a.modes)
    emit(test='settings',ensemble=a.ensemble,mesh=a.mesh,modes=a.modes,
         flux_columns=['number','sideways_advection','sideways_total','energy_at_common_mu_one','side_pressure_only'],
         original_bias=np.angle(field[dom.bulk][0]),differential_conductance=solver.conductance)
    for bend in (0.,.01,-.01,.005,-.005):
        offset=np.array([0.,bend])
        row=solver.solve(offset,.25) if a.ensemble=='intake' else solver.fixed_phase(offset)
        emit(test='native',ensemble=a.ensemble,**row,
             shell=cuts(solver,offset,1.173,2.067),inner=cuts(solver,offset,1.31,1.87))
        if abs(row['bulk']-row['feed'])>3e-9:
            raise ArithmeticError('Reservoir number budget failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34,48),required=True)
    p.add_argument('--ensemble',choices=('intake','phase'),required=True)
    p.add_argument('--modes',type=int,default=3)
    a=p.parse_args()
    a.source=Path(f'_scratch/s11c/mouth_flared_delayed_n{a.mesh}')
    main(a)
