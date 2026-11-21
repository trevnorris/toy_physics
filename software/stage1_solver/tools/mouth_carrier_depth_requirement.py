"""Inverse requirement on the branch where weak intake/inertia matches the reference.

The target ratio is imposed only to falsify a proposed class of passive depth
laws. Compare the axial load at the required depth with all positive mixtures
of the two independently calibrated size-one restorations.
"""
import argparse
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve
from mouth_carrier_depth_selection import DepthDomain
from mouth_carrier_depth_work import axial_work,geo_work
from mouth_carrier_weak_radius import translation
from mouth_flared_receiver import EXIT,emit


def run(a):
    folder=Path('_scratch/s11c/carrier_depth_one_01')
    if json.loads((folder/'outcome.json').read_text())['exitCode'] or not json.loads((folder/'limit-validation.json').read_text())['verified']:
        raise ValueError('Invalid unit-radius reference')
    ref=json.loads((folder/'stdout').read_text())
    def static(h):
        roof=h*round(EXIT+16/h,9)
        dom=DepthDomain(a.radius,h,a.mesh,24.,roof,1.)
        phase=dom.bulk.astype(float)
        free=dom.free
        phase[free]=spsolve(dom.stiff[free][:,free],-(dom.stiff@phase)[free])
        current=dom.stiff@phase
        c=float(phase@current)
        np.testing.assert_allclose(c,-current[dom.feed].sum(),rtol=3e-7)
        return dom,phase,c
    @lru_cache(None)
    def ratio(h):
        dom,_,c=static(h)
        m,_=translation(dom)
        return c/m
    grid=[.25,.4,.6,.8,1.,1.25,1.6,2.,2.5]
    samples=[ratio(h) for h in grid]
    emit(test='inverse_scan',r=a.radius,depths=grid,ratios=samples)
    target=ref['ratio_per_bias']
    intervals=[(x,y) for x,y,v,w in zip(grid,grid[1:],samples,samples[1:])
               if (v-target)*(w-target)<=0]
    if len(intervals)!=1:
        raise ArithmeticError('Expected one resolved branch for the local inverse test')
    lo,hi=intervals[0]
    h=brentq(lambda x:ratio(x)-target,lo,hi,xtol=4e-6)
    dom,phase,c=static(h)
    step=.004*h
    derivative=(static(h+step)[2]-static(h-step)[2])/(4*step)
    actual=-axial_work(dom,phase,a.radius,h)
    v,area=geo_work(dom,a.radius,h)
    pressure=ref['required_pressure']*v
    tension=ref['required_tension']*area
    fraction=(actual-tension)/(pressure-tension)
    emit(test='inverse_requirement',radius=a.radius,mesh=a.mesh,depth=h,
         checked_ratio=ratio(h),reference_ratio=target,required_restoration=actual,
         energy_restoration=-derivative,pressure_only=pressure,tension_only=tension,
         required_pressure_fraction_at_reference=fraction)
    np.testing.assert_allclose(actual,-derivative,rtol=.02)
    np.testing.assert_allclose(ratio(h),target,rtol=4e-5)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--mesh',type=int,default=26)
    run(p.parse_args())
