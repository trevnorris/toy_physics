"""Matched-intake titration of an external energy store on the stationary deep order branch.

The feed carries h times H(chi) of convertible excitation per constituent.
This trial adds no independent force on a wall spring. The carrier is held at
its archived current to isolate how the order layer responds to the new drive.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq

from mouth_bend_renewal import ROOT
from mouth_carrier_lateral_response import prepare
from mouth_empty_bend_stiffness import quadrature
from mouth_free_order_diffuse import positive_factor, read, section_geometry
from mouth_free_order_residual_check import DifferenceOrder
from mouth_order_wave_pulse import digest, emit


class SuppliedOrder(DifferenceOrder):
    def __init__(self, dom, psi, width, mobility):
        super().__init__(dom, psi, width)
        self.density = abs(dom.N @ psi)**2
        self.mobility = mobility

    def equilibrium(self, seed, head):
        x = seed.copy()
        for iteration in range(30):
            c = (self.dom.N @ x)/self.quad_sum
            d = (self.dom.N @ (1-x))/self.quad_sum
            res,jac = self.measure(x)
            drive = self.wt*self.mobility*head*self.density
            res += self.Nf.T @ (drive*6*c*d)
            jac += self.mass_weight(drive*6*(d-c))
            factor = positive_factor(jac)
            if factor is None:
                return None
            step = factor.solve(-res)
            error = max(abs(step))
            if error < 5e-12:
                return x, error, iteration
            x[self.dom.free] += step*min(1., .1/error)
        return None


def main(a):
    settings=ROOT/f'mouth_free_order_delayed_n{a.mesh}_d04'
    source=ROOT/f'mouth_flared_delayed_n{a.mesh}'
    meta=json.loads((settings/'stdout').read_text().splitlines()[0])
    archive=json.loads(Path('research/pde_ledger_v3/_measurements/mouth_free_order_depth.json').read_text())
    seed_file=settings/'order_seed1.55.npz'
    if digest(seed_file)!=archive['diffuse_orders_sha256'][str(seed_file)] or digest(source/'field.npz')!=meta['source_sha256']:
        raise ValueError('Mismatched inputs')
    dom,psi=prepare(SimpleNamespace(source=source))
    with np.load(seed_file) as data:
        base=data['chi']
    Q=meta['source_flux']
    model=SuppliedOrder(dom,psi,meta['delta'],meta['mobility'])
    _,_,dr,dw=quadrature(dom)
    f=dom.N@psi
    ur,uw=(np.imag(g@(psi-1)/f) for g in (dr,dw))
    wt=dom.weight
    def readings(x,h):
        c=(dom.N@x)/model.quad_sum
        rate=ur*(dr@x)+uw*(dw@x)
        converted=wt@(-model.density*6*c*(1-c)*rate)
        loss=wt@(rate**2)/meta['mobility']
        depth=brentq(lambda z:read(dom,x,section_geometry(dom,z)[0],z)[0]-.5,.51,2.5)
        return dict(converted_number_rate=converted,head_power=Q*h,volume_supply=h*converted,
                    passive_loss=loss,unpaid_loss=loss-Q*h,wall_half_order_depth=depth,
                    quadrature_order_range=[c.min(),c.max()])
    flat_area=4*np.pi/3
    # Independent calibration on a flat unit-radius counting section,
    # n=1, with the same width, mobility and total inward number rate.
    h_flat=Q/(flat_area*6*meta['delta']*meta['mobility'])
    emit(test='calibration',mesh=a.mesh,intake=Q,radius=1.,order_source=str(seed_file),order_sha256=digest(seed_file),
         mobility=meta['mobility'],width=meta['delta'],area_flat=flat_area,flat_head=h_flat,
         definition='unstrained constant-density planar front, equal number rate and radius')
    initial=model.equilibrium(base,0.)
    if initial is None:
        raise ValueError('Original stationary branch missing')
    x=initial[0]
    h=0.
    a.out.mkdir(parents=True,exist_ok=True)
    for multiple in a.levels:
        target=multiple*h_flat
        step=.002
        lost=False
        while target-h>1e-14:
            nxt=min(target,h+step)
            solved=model.equilibrium(x,nxt)
            if solved is None:
                step/=2
                if step<1.e-5:
                    emit(test='continuation_limit',last_locally_stable_head=h,requested_head=target,
                         next_unresolved_head=nxt,**readings(x,h))
                    lost=True
                    break
                continue
            x=solved[0]
            h=nxt
            step=min(.004,step*1.5)
        if lost:
            break
        row=readings(x,h)
        file=a.out/f'order_{multiple:.5f}.npz'
        np.savez_compressed(file,chi=x)
        emit(test='supplied_order',multiple=multiple,head=h,path=str(file),sha256=digest(file),**row)
        if abs(row['converted_number_rate']/Q-1)>1.e-3:
            raise ArithmeticError('Independent volume and endpoint conversion rates disagree')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(34,48),required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--levels',type=float,nargs='+',default=[0,.05,.1,.2,.5,.8,.9,.95,1.,1.05,1.1,1.2])
    a=p.parse_args()
    if a.levels!=sorted(a.levels) or a.levels[0]!=0:
        raise ValueError('Continue monotonically from zero')
    main(a)
