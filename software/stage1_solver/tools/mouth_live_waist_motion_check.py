"""Evaluate the carrier Routhian again at signed order velocities and changed order."""
import json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import spsolve
from mouth_live_waist_fields import BIAS, Carrier, ORDER_INERTIA
from mouth_live_waist_motion import CarrierMotion
from mouth_sealed_shape_relax import Meridian

ROOT=Path(__file__).resolve().parents[3]


def main():
    archived=json.loads((ROOT/'research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json').read_bytes())
    for radius in (1,2):
        state=archived['training'][str(radius)]['states'][-2]
        shape=Meridian(radius,state['coordinates'])
        n=48
        x,t=np.meshgrid(np.linspace(0,1,n+1),np.linspace(0,1,2*n+1))
        order=.5+.49*np.tanh((1-x)*shape.profile(t)[0]/.08)*np.tanh((1-t)*shape.length/.08)
        eta=np.sin(np.pi*t)*np.cos(np.pi*x/2)
        c=Carrier(shape,n,2*n,14)
        _,g,_=c.response(order,BIAS)
        m=CarrierMotion(c,order,eta)
        out=m.matched_mass()
        h=.002
        slope=(m.held(h)['energy_change']-m.held(-h)['energy_change'])/(2*h)
        print(json.dumps(dict(radius=radius,mode_masses=out,frozen_shape_slope=slope,
                              variational_shape_slope=np.sum(g*eta))),flush=True)
        assert abs(slope-np.sum(g*eta))<1.e-10
        dr,dw=m.grad_basis
        weight=m.coeff(m.gr,m.gw)
        terms=np.zeros_like(m.phase)
        for col in range(4):
            np.add.at(terms,m.ids[:,col],weight*m.eta*(m.gr*dr[:,col]+m.gw*dw[:,col]))
        rates=[-.125,0.,.125]
        energies=[]
        for rate in rates:
            p=np.zeros_like(m.phase);p[c.bulk]=BIAS
            p[c.free]=spsolve(m.A[c.free][:,c.free],(rate*terms-m.A@p)[c.free])
            # Direct on-shell difference by subtracting each local quadratic form.
            v=p[m.ids]
            delta=p-m.phase
            L=-.5*(delta@(m.A@delta))-delta@(m.A@m.phase)+rate*(terms@p)
            L+=.5*rate**2*np.sum(weight*m.eta**2*m.vol)
            energies.append(L)
        dL=(energies[-1]-energies[0])/.25
        d2L=(energies[-1]+energies[0]-2*energies[1])/.125**2
        print(json.dumps(dict(radius=radius,finite_velocity_momentum=dL,finite_velocity_mass=d2L)),flush=True)
        if abs(dL-out['momentum_offset'])>1e-9 or abs(d2L-out['matched_mass'])>1e-9 or not 0<d2L<=1.001*ORDER_INERTIA*np.sum(m.eta**2*m.vol):
            raise ArithmeticError('Material velocity and phase-condensed inertia checks disagree')


if __name__=='__main__':
    main()
