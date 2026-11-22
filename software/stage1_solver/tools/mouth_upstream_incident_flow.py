"""Sideways momentum given to an unbent smooth mouth by a transverse incident stream.

Compare direct Cartesian trajectories to radial sensitivity and conserved angular
momentum on a sphere. The stream starts at the upstream cut; its mass flow is the
normalization. No later release or extra recipient of momentum is imposed.
"""
import argparse
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from mouth_upstream_normal_contact import Crossing, coefficients
from neutral_patch_sites import BETWEEN, TEMPERATURE, emit


def cartesian(tr,cuts,scale,args):
    master=Crossing(tr,cuts,scale,SimpleNamespace(law='surface',origin=args.origin,amplitude=1.,
                                angles=args.angles,radial=args.radial))
    master.state[2]+=master.q
    master.inlet_v=master.mean(master.state[2])
    emit(test='incident_setup',mesh=args.mesh,speed=master.q.ravel().tolist(),radial_nodes=args.radial,
         mu_nodes=args.angles,dt=args.dt,origin=args.origin)
    edges=np.unique(np.r_[cuts,tr.edges[(tr.edges>cuts[0])&(tr.edges<cuts[-1])]])
    for a,b in zip(edges[:-1],edges[1:],strict=True):
        n=max(1,int(np.ceil(tr.transit(a,b,8)/args.dt)))
        nodes=np.linspace(a,b,n+1)
        mids=np.linspace(a,b,2*n+1)
        mids[0]=np.nextafter(a,b)
        mids[-1]=np.nextafter(b,a)
        band=coefficients(tr,mids,cuts[0],cuts[2],False)
        band.bend.fill(0.)
        band.bend_deriv.fill(0.)
        master.advance(nodes,band)
        if b in cuts:
            gained=master.mean(master.state[2])-master.inlet_v
            values=-gained
            emit(test='incident_crossing',depth=b,chi=band.order[-1],
                 wall_impulse_per_flux=values.tolist(),
                 captured_per_incident_momentum=(values[1]-values[2])/(2*master.q[1,0]),
                 captured_half=(values[3]-values[4])/(2*master.q[3,0]),
                 independent_virtual_force_half=-(master.impulse[3]-master.impulse[4])/(2*master.q[3,0]),
                 near_fractions=master.mean(master.active.astype(float)).tolist())


def invariant(tr,cuts,scale,args):
    # A small incident velocity U e has radial part U mu and transverse part U P_n e.
    # On the unbent guide the latter follows conservation of angular momentum.
    if args.normal_rule=="finite":
        if args.radial<1: raise ValueError("A finite-interval Gaussian needs quadrature nodes")
        nodes,weight=leggauss(args.radial)
        standard=9*nodes
        weight*=np.exp(-.5*standard**2)
        weight/=weight.sum()
        launch=standard*np.sqrt(scale*TEMPERATURE)
    elif args.radial:
        nodes,weight=hermgauss(args.radial)
        weight/=np.sqrt(np.pi)
        launch=nodes*np.sqrt(2*scale*TEMPERATURE)
    else:
        launch=np.zeros(1)
        weight=np.ones(1)
    R0,Rprime=tr.geometry.radius(cuts[0])
    speed0=tr.speed(cuts[0])
    state=np.zeros((5,len(weight)))
    state[0]=R0
    state[1]=launch+Rprime*speed0
    state[3]=1. # radial sensitivity of velocity to incident radial velocity
    D,alpha,_=BETWEEN
    edges=np.unique(np.r_[cuts,tr.edges[(tr.edges>cuts[0])&(tr.edges<cuts[-1])]])
    for a,b in zip(edges[:-1],edges[1:],strict=True):
        def ode(w,y):
            r,v,s,z,_=y.reshape(state.shape)
            depth=np.clip(w,np.nextafter(a,b),np.nextafter(b,a))
            chi=tr.sample(tr.chi,depth)[0]
            R=tr.geometry.radius(depth)[0]
            ex=np.exp(-alpha*(r-R)/PITCH)
            scale0=2*D*alpha*scale*chi/PITCH
            radial=scale0*ex*(ex-1)
            gradient=scale0*alpha/PITCH*ex*(1-2*ex)
            return (np.array([v,radial,z,gradient*s,1/r**2])/tr.speed(depth)).ravel()
        sol=solve_ivp(ode,(a,b),state.ravel(),method='DOP853',rtol=2e-10,atol=1e-12,t_eval=[b])
        if not sol.success: raise ArithmeticError(sol.message)
        state=sol.y[:,-1].reshape(state.shape)
        if b in cuts:
            r,v,_,z,I=state
            angular=R0/r+R0*v*I
            radial_gauss=float(weight@(launch*v)/(scale*TEMPERATURE)) if args.radial else None
            emit(test='angular_momentum_check',depth=b,chi=tr.sample(tr.chi,b)[0],radial_nodes=args.radial,normal_rule=args.normal_rule,
                 captured_per_incident_momentum=1-weight@(z+2*angular)/3,
                 radial_derivative=weight@z,angular_derivative=weight@angular,
                 thermal_radial_by_parts=radial_gauss,
                 captured_by_thermal_integration=None if radial_gauss is None else 1-(radial_gauss+2*weight@angular)/3)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--radial',type=int,required=True)
    p.add_argument('--method',choices=['cartesian','invariant'],required=True)
    p.add_argument('--normal-rule',choices=['hermite','finite'],default='hermite')
    p.add_argument('--angles',type=int,default=16)
    p.add_argument('--dt',type=float,default=.0006)
    p.add_argument('--origin',type=float,default=0.)
    args=p.parse_args()
    tr,rate=cubic_trace(args.mesh)
    lo,hi=1.175,tr.crossing(.5)
    cuts=np.array([lo,(lo+hi)/2,hi,tr.crossing(.25),tr.crossing(.1),tr.crossing(.01)])
    patch,_,_=read_inputs(args.mesh)
    emit(test='incident_material',mesh=args.mesh,radial_nodes=args.radial,intake=rate,scale=2*patch['shared_scale'])
    (cartesian if args.method=='cartesian' else invariant)(tr,cuts,2*patch['shared_scale'],args)
