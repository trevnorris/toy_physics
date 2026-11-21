"""Opposite finite bends of the wave-empty nonlinear draining carrier at fixed total intake.

Meridional Q2 fields multiply Legendre polynomials about the displacement axis.
The two far boundaries stay fixed. Their phase difference is retuned at each
offset; both reservoir currents and the mechanical holding forces are outputs.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import Legendre,leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_curved_phase_stability import weighted_gram
from mouth_empty_bend_stiffness import emit,prepare,quadrature,shapes


class BentCarrier:
    def __init__(self,dom,base,nmodes):
        self.dom=dom
        self.r,w,self.dr,self.dw=quadrature(dom)
        self.base=base
        _,self.sr,self.sw=shapes(self.r,w,'compact',dom.far,dom.roof)
        self.modes=nmodes
        self.mu,self.weights=leggauss(max(8,5*(nmodes-1)+2))
        self.weights/=2
        polynomials=[Legendre.basis(l) for l in range(nmodes)]
        self.y=np.array([p(self.mu) for p in polynomials]).T
        self.dy=np.array([p.deriv()(self.mu) for p in polynomials]).T
        self.preconditioners=[]
        f=dom.N@base; n=abs(f)**2
        wt=dom.weight
        gram=lambda v: weighted_gram(dom.N,wt*v)
        coupling=gram(2*n**3*f.real*f.imag)
        for l in range(nmodes):
            mask=np.tile(dom.free & ((dom.points[:,0]>0)|(l==0)),2)
            scalar=.5*dom.stiff+gram((n**4-1)/4 + .5*l*(l+1)/self.r**2)
            jac=sparse.bmat([[scalar+gram(2*n**3*f.real**2),coupling],
                             [coupling,scalar+gram(2*n**3*f.imag**2)]],format='csr')
            fac=splu(jac[mask][:,mask].tocsc())
            self.preconditioners.append((fac,mask))
            if l==0:
                unit=np.zeros_like(base); unit[dom.bulk]=1j*base[dom.bulk]
                pair=np.r_[unit.real,unit.imag]
                pair[mask]=fac.solve(-(jac@pair)[mask])
                x,y=np.split(pair,2)
                self.bias_response=x+1j*y
                x,y=np.split(jac@pair,2)
                bc=self.evaluate(np.column_stack([base]+[np.zeros_like(base)]*(nmodes-1)),[0.,0.])[0][:,0]
                self.conductance=2*np.sum((self.bias_response.conj()*bc+base.conj()*(x+1j*y))[dom.bulk].imag)
        if self.conductance<=0: raise ArithmeticError('Nonpositive differential conductance')

    def evaluate(self,state,offset):
        dom=self.dom
        nval=dom.N@state
        delta=state.copy(); delta[:,0]-=1
        qr,qw=self.dr@delta,self.dw@delta
        br=self.sr@offset; bw=self.sw@offset; bb=br**2+bw**2
        field_terms=np.zeros_like(nval); radial=np.zeros_like(nval); axial=np.zeros_like(nval)
        force=np.zeros(2); energy=0.
        for m,aw,y,dy in zip(self.mu,self.weights,self.y,self.dy,strict=True):
            f=nval@y; fr=qr@y; fw=qw@y; fu=nval@dy
            z=m*fr+(1-m*m)*fu/self.r
            metric=1+m*br
            if metric.min()<=0: raise ValueError('Folded deformation')
            h=br*fr+bw*fw
            corr=-h+bb*z/metric
            nr=abs(f)**2; rho=nr-1
            pot=rho**2*(nr**3+2*nr**2+3*nr+4)/20
            flat=.5*(abs(fr)**2+abs(fw)**2+(1-m*m)*abs(fu/self.r)**2)+pot
            radial+=aw*(metric*fr-br*z+m*corr)[:,None]*y
            axial+=aw*(metric*fw-bw*z)[:,None]*y
            field_terms+=aw*((metric*(nr**4-1)*f/2)[:,None]*y+
                        ((1-m*m)*(metric*fu/self.r**2+corr/self.r))[:,None]*dy)
            energy+=aw*(dom.weight@(metric*flat-(z.conj()*h).real+.5*bb*abs(z)**2/metric))
            p=abs(z)**2/metric
            first=m*flat-(z.conj()*fr).real+p*(br-.5*bb*m/metric)
            second=-(z.conj()*fw).real+p*bw
            force+=aw*((dom.weight*first)@self.sr+(dom.weight*second)@self.sw)
        residual=.5*(self.dr.T@(dom.weight[:,None]*radial)+self.dw.T@(dom.weight[:,None]*axial)
                     +dom.N.T@(dom.weight[:,None]*field_terms))
        residual*=2*np.arange(self.modes)+1
        return residual,energy,force

    def solve(self,offset,rate):
        f=np.zeros((len(self.base),self.modes),complex); f[:,0]=self.base
        dom=self.dom
        theta=np.angle(f[dom.bulk,0][0])
        for it in range(45):
            residual,energy,force=self.evaluate(f,offset)
            flow=2*(f[:,0].conj()*residual[:,0]).imag
            bulk=flow[dom.bulk].sum(); feed=-flow[dom.feed].sum()
            error=0.; updates=[]
            for l,(fac,mask) in enumerate(self.preconditioners):
                e=residual[:,l]; raw=np.r_[e.real,e.imag]
                error=max(error,np.max(abs(raw[mask])))
                step=np.zeros_like(raw); step[mask]=fac.solve(-raw[mask]); x,y=np.split(step,2)
                updates.append(x+1j*y)
            if error<3e-11:
                if abs(bulk-rate)<1e-10 and abs(feed-rate)<2e-9:
                    return dict(offset=list(offset),phase_bias=theta,feed=feed,bulk=bulk,energy=energy,
                                fixed_rate_potential=energy-rate*theta,holding_force=force.tolist(),
                                residual=error,iterations=it)
                shift=(rate-bulk)/self.conductance
                theta+=shift; f[:,0]+=shift*self.bias_response
                f[dom.bulk,0]=np.exp(1j*theta)
            else:
                f+=np.column_stack(updates)
        raise ArithmeticError(f'Unconverged finite bend: residual={error:g} feed={feed:g} bulk={bulk:g}')


def main(a):
    dom,base=prepare(a)
    solver=BentCarrier(dom,base,a.modes)
    direction=np.array([float(v) for v in a.direction.split(',')])
    rate=.25
    emit(test='settings',modes=a.modes,direction=direction.tolist(),differential_conductance=solver.conductance)
    rows={}
    for offset in (0.,a.step,-a.step,a.step/2,-a.step/2):
        row=solver.solve(offset*direction,rate)
        rows[offset]=row
        emit(test='finite_offset',**row)
    for size in (a.step,a.step/2):
        plus,minus=rows[size],rows[-size]
        f=(np.array(plus['holding_force'])-np.array(minus['holding_force']))/(2*size)
        centered=lambda key:(plus[key]+minus[key]-2*rows[0.][key])/size**2
        emit(test='secant',modes=a.modes,step=size,force_column=f.tolist(),
             directional_virtual_work=direction@f,retuned_energy_second_difference=centered('energy'),
             fixed_rate_second_difference=centered('fixed_rate_potential'),
             mesh=a.mesh,direction=direction.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(profile='cubic',rate=.25,source=None)
    p.add_argument('--mesh',type=int,default=18)
    p.add_argument('--modes',type=int,default=2)
    p.add_argument('--direction',required=True)
    p.add_argument('--step',type=float,default=.02)
    p.add_argument('--far',type=float,default=9.)
    p.add_argument('--roof',type=float,default=10.)
    main(p.parse_args())
