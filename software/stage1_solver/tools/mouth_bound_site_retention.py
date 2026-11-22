"""Conditional first passage of a cold surface constituent in a prescribed stream.

The reversible site energy follows the local stationary adsorption branch
quasistatically. The normal guide is fixed. The single longitudinal coordinate
is wall arclength, and the drag per site mass is gamma*(u-v).
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_surface_site_transfer import BASIS, Surface
from neutral_patch_sites import emit


def first_passage(s, w, potential, gamma, speed, transit):
    vfield = PchipInterpolator(s,potential)
    force = vfield.derivative()
    end=s[-1]
    def rhs(t, y):
        x,v,loss=y
        return v,-float(force(np.clip(x,0,end)))-gamma*v,gamma*v*v
    def upstream(t,y):
        return y[0]
    def downstream(t,y):
        return y[0]-end
    # Prefix/suffix maxima of the shape-preserving interpolant occur at nodes
    # or at the current location. A barrier at the current location itself
    # cannot certify a strict trap. Only interior states with two resolved
    # higher cuts and declining mechanical energy can end as confined.
    def barriers(x):
        j=np.searchsorted(s,np.clip(x,0,end),side='right')
        left=int(np.argmax(potential[:max(j,1)]))
        right=min(j,len(s)-1)+int(np.argmax(potential[min(j,len(s)-1):]))
        return left,right
    def confinement(t,y):
        x,v,_=y
        if x<=0 or x>=end:
            return -v*v/2-1e-8
        l,r=barriers(x)
        barrier=min(potential[l],potential[r])
        bottom=np.min(potential[l:r+1])
        # Resolve 5% of the well depth at the event, up to 1e-4 specific energy.
        margin=max(1e-8,min(1e-4,.05*(barrier-bottom)))
        return barrier-float(vfield(x))-.5*v*v-margin
    for event,direction in [(upstream,-1),(downstream,1),(confinement,1)]:
        event.terminal=True
        event.direction=direction
    y=np.array([0.,speed,0.])
    energy0=.5*speed**2+potential[0]
    t=err=0.
    while True:
        sol=solve_ivp(rhs,(t,t+transit),y,events=(upstream,downstream,confinement),
                      method='DOP853',atol=1e-12,rtol=1e-9,max_step=transit/128)
        if not sol.success:
            raise ArithmeticError(sol.message)
        ledger=.5*sol.y[1]**2+vfield(sol.y[0])+sol.y[2]-energy0
        err=max(err,float(np.max(abs(ledger))))
        y,t=sol.y[:,-1],sol.t[-1]
        if sol.status==1:
            break
    h=float(vfield(y[0]))+.5*y[1]**2
    event=next(i for i,a in enumerate(sol.t_events) if len(a))
    out=dict(outcome=('returned_to_entry','passed_half_order','confined_between_cuts')[event],
             time=t,depth=np.interp(y[0],s,w),speed=y[1],effective_energy=h,
             integrated_specific_dissipation=y[2],energy_ledger_max_error=err)
    if event==2:
        l,r=barriers(y[0])
        out.update(barrier_depths=[w[l],w[r]],energy_margin=min(potential[l],potential[r])-h)
        if out['energy_margin']<=20*err:
            raise ArithmeticError('Energy margin does not resolve confinement')
    if err>5e-8:
        raise ArithmeticError('Unresolved energy ledger')
    return out


def run(args):
    tr,intake=cubic_trace(args.mesh)
    patch,_,_=read_inputs(args.mesh)
    a,b=1.175,tr.crossing(.5)
    w=np.unique(np.r_[np.linspace(a,b,args.segments+1), (a+b)/2,
                      tr.edges[(tr.edges>a)&(tr.edges<b)]])
    hi=np.where(w==(a+b)/2)[0].item()
    roots,weights=leggauss(6)
    half=np.diff(w)/2
    z=w[:-1,None]+half[:,None]*(1+roots)
    ds=half*(np.hypot(1.,tr.geometry.radius(z)[1])@weights)
    s=np.r_[0.,np.cumsum(ds)]
    theta=np.unwrap(np.angle(tr.sample(tr.psi,w)[0]))
    theta-=theta[0]
    chi=tr.sample(tr.chi,w)[0]
    speed=tr.speed(w)[0]*np.hypot(1.,tr.geometry.radius(w)[1][0])
    transit=tr.transit(a,b,8)
    sample=Surface(args.side,5.,patch['pitch'],False)
    start=np.zeros(sample.dim)
    start[:2]=patch['pitch']*BASIS.sum(axis=0)/3
    start[2]=.75
    scale=2*patch['shared_scale']
    emit(test='path',mesh=args.mesh,intake=intake,host_side=args.side,path_knots=len(w),
         segments=args.segments,entry_speed=speed,carrier_transit=transit,mid_depth=w[hi])
    if np.any(np.diff(theta)<=0) or np.any(np.diff(chi)>=0):
        raise ArithmeticError('Path must have increasing phase and decreasing order')
    for mixed,label in [(False,'weaken_all_contacts'),(True,'retain_mean_attraction')]:
        energies=[]
        state=start
        for order in chi:
            state,e,_,_=sample.equilibrate(order,mixed,state)
            energies.append(e*scale)
        mu=np.array(energies)-energies[0]
        lo=mu[hi]/theta[hi]
        hi_rate=(mu[-1]-mu[hi])/(theta[-1]-theta[hi])
        emit(test='midpoint',law=label,depths=w[[0,hi,-1]].tolist(),phase=theta[[0,hi,-1]].tolist(),
             mu=mu[[0,hi,-1]].tolist(),lower_rate=lo,upper_rate=hi_rate,
             midpoint_certifies_an_interior_minimum=bool(lo<hi_rate))
        for gamma in args.gamma:
            vpot=mu-gamma*theta
            j=np.argmin(vpot)
            result=first_passage(s,w,vpot,gamma,speed,transit)
            emit(test='launch',law=label,gamma=gamma,min_knot_depth=w[j],
                 midpoint_below_entry=vpot[0]-vpot[hi],midpoint_below_exit=vpot[-1]-vpot[hi],**result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',required=True,type=int,choices=(26,34))
    p.add_argument('--segments',type=int,required=True)
    p.add_argument('--side',type=int,default=5)
    p.add_argument('--gamma',type=float,nargs='+',default=[0.,3.,5.5,12.])
    a=p.parse_args()
    if a.segments<4 or a.segments%2 or a.side<3 or not all(np.isfinite(g) and g>=0 for g in a.gamma):
        p.error('Use an even positive path resolution, side >=3 and finite nonnegative rates')
    run(a)
