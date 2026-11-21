"""Virtual work to change the depth of a coherent throat at a fixed outer phase.

A local volumetric work cost and constant wall tension are separate hypotheses.
Geometry beyond the exit translates with the exit, keeping bath clearance fixed.
Quantities are reported per squared weak phase difference.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.sparse.linalg import spsolve

from mouth_carrier_lateral_response import trace
from mouth_carrier_weak_radius import PhysicalDomain, translation
from mouth_flared_receiver import EXIT, ROUNDING, FLARE_START, emit


def shape(w,dom):
    if w<ROUNDING:
        s=np.sqrt(w*(2*ROUNDING-w))
        return dom.rc-s,(w-ROUNDING)/s
    return dom.radius(w)


def geo_work(dom,radius,depth):
    def integrate(kind):
        def integrand(w):
            b,bp=shape(w,dom)
            if kind=='volume':
                return 4*np.pi/3*radius**3*b**3
            return 4*np.pi*radius**2*b**2*depth/np.hypot(depth,radius*bp)
        return sum(quad(integrand,a,b,epsabs=1e-10)[0] for a,b in
                   ((0.,ROUNDING),(ROUNDING,FLARE_START),(FLARE_START,EXIT)))
    return integrate('volume'),integrate('area')


def axial_work(dom,phase,radius,depth):
    # Smooth physical transverse slices through the tube. Negative-w feed
    # and the downstream chamber undergo only translations or stay at rest.
    vals=[]
    x,wt=leggauss(5)
    pg=dom.blocks[2][0]
    levels=np.r_[np.linspace(0,ROUNDING,45)[:-1],pg]
    for a,b in zip(levels[:-1],levels[1:]):
        for w,dw in zip((a+b)/2+(b-a)*x/2,wt*(b-a)/2):
            bound,_=shape(w,dom)
            # Smooth in a single element, with independent radial quadrature.
            radial=(np.arange(48)[:,None]+(x+1)/2)/48
            rr=(bound*radial).ravel()
            weights=np.tile(wt/96,48)*bound
            if w<ROUNDING:
                p=np.arctan2(w-ROUNDING,dom.rc-rr)
                span,_=dom.reach(p)
                v=(np.hypot(w-ROUNDING,dom.rc-rr)-ROUNDING)/(span-ROUNDING)
                block=0
            else:
                v=1-rr/bound
                p=np.full_like(v,w)
                block=2
            r_actual,_,_,_,_,_,_,dr,dz=trace(dom,block,p,v,phase[:,None])
            np.testing.assert_allclose(r_actual,rr*radius,atol=1e-10)
            vals.append(.5*dw*np.dot(weights*radius*4*np.pi*r_actual**2,
                                     dr[:,0]**2-dz[:,0]**2))
    return sum(vals)


def solve(radius,depth,mesh,outer,clearance,interior):
    # Round only the remote chart coordinate, within 1e-9 of the requested
    # physical roof. This avoids ambiguous membership at a prescribed bath.
    roof=depth*round(EXIT+clearance/depth,9)
    dom=PhysicalDomain(radius,depth,mesh,outer,roof,1.)
    mask=dom.free
    phase=dom.bulk.astype(float)
    phase[mask]=spsolve(dom.stiff[mask][:,mask],-(dom.stiff@phase)[mask])
    j=dom.stiff@phase
    c=phase@j
    if abs(c+j[dom.feed].sum())>2e-7*c:
        raise ArithmeticError('Independent bath rate and energy disagree')
    if not interior:
        return .5*c
    v,a=geo_work(dom,radius,depth)
    m,_=translation(dom)
    work=axial_work(dom,phase,radius,depth)
    return .5*c,work,v,a,m


def measure(args,h):
    vals={}
    for offset in (-.04,-.02,0.,.02,.04):
        vals[offset]=solve(args.radius,h*(1+offset),args.mesh,args.outer,args.clearance,offset==0)
    energy,local,v,a,m=vals[0.]
    fine=(vals[.02]-vals[-.02])/(.04*h)
    coarse=(vals[.04]-vals[-.04])/(.08*h)
    derivative=(4*fine-coarse)/3
    curv1=(vals[.02]+vals[-.02]-2*energy)/(.02*h)**2
    curv2=(vals[.04]+vals[-.04]-2*energy)/(.04*h)**2
    curvature=(4*curv1-curv2)/3
    row=dict(radius=args.radius,depth=h,mesh=args.mesh,outer=args.outer,clearance=args.clearance,
             energy=energy,conductance=2*energy,M=m,ratio_per_bias=2*energy/m,
             gradient=derivative,local_stress_gradient=local,derivatives=[coarse,fine],
             carrier_curvature=curvature,volume_gradient=v,area_gradient=a,
             required_pressure=-derivative/v,required_tension=-derivative/a)
    emit(test='depth_work',**row)
    if not np.isclose(derivative,local,rtol=.015,atol=1e-8) or m<=0 or not np.isfinite(curvature):
        raise ArithmeticError('Virtual axial work or positive inertia check failed')
    return row


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--mesh',type=int,default=26)
    p.add_argument('--outer',type=float,default=14.)
    p.add_argument('--clearance',type=float,default=11.95)
    p.add_argument('--depths',type=float,nargs='+',required=True)
    a=p.parse_args()
    for h in a.depths:
        measure(a,h)
