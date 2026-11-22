"""Transport small transverse deviations from an ordered inlet to the deep flare.

For axial motion along a lattice row, the transverse force at exact alignment
vanishes by reflection. Determine the response to perturbing either inlet
position or inlet velocity without assigning fresh sites a new row downstream.
"""
import argparse

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline, PchipInterpolator

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, GAP, emit, morse


def transverse_hessian(pitch):
    axis=np.linspace(0.,pitch,513)
    curve=np.zeros_like(axis)
    depth,alpha,rest=BETWEEN
    for i in range(-10,11):
        for j in range(-10,11):
            x=axis-(i+j*.5)*pitch
            y=j*np.sqrt(3.)/2*pitch
            r=np.sqrt(x*x+y*y+GAP**2)
            _,first=morse(r,*BETWEEN)
            exp=np.exp(-alpha*(r-rest))
            second=2*depth*alpha**2*exp*(2*exp-1)
            curve+=first/r+(second-first/r)*(y/r)**2
    np.testing.assert_allclose(curve[0],curve[-1],atol=1e-12)
    curve[-1]=curve[0]
    return CubicSpline(axis,curve,bc_type='periodic',extrapolate='periodic')


def transfer(mesh,substeps,phases):
    trace,_=cubic_trace(mesh)
    patch,_,_=read_inputs(mesh)
    pitch=patch['pitch']
    curvature=transverse_hessian(pitch)
    scale=2*patch['shared_scale']/PITCH**2
    longitudes=(np.arange(phases)+.5)/phases*pitch
    state=np.repeat(np.eye(2)[None,:,:],phases,axis=0)
    elapsed=0.
    length=0.
    for lo,hi in ((1.175,1.55),(1.55,1.65)):
        edges=np.unique(np.r_[np.linspace(lo,hi,8193),trace.edges[(trace.edges>lo)&(trace.edges<hi)]])
        clock=cumulative_trapezoid(1/trace.speed(edges),edges,initial=0.)
        distance=cumulative_trapezoid(np.hypot(1.,trace.geometry.radius(edges)[1]),edges,initial=0.)
        time=np.linspace(0.,clock[-1],int(np.ceil(clock[-1]*substeps))+1)
        w=PchipInterpolator(clock,edges)(time)
        arc=PchipInterpolator(clock,distance)(time)
        chi=trace.sample(trace.chi,w)[0]
        k=scale*chi[:,None]*curvature((length+arc[:,None])/PITCH+longitudes)
        dt=time[1]
        for step in range(len(time)-1):
            state[:,1]-=.5*dt*k[step,:,None]*state[:,0]
            state[:,0]+=dt*state[:,1]
            state[:,1]-=.5*dt*k[step+1,:,None]*state[:,0]
        elapsed+=clock[-1]
        length+=distance[-1]
        det_error=np.max(abs(np.linalg.det(state)-1.))
        if det_error>1.e-7:
            raise ArithmeticError('Transverse phase-area transfer failed')
        emit(test='upstream_row_control',mesh=mesh,substeps_per_time=substeps,phases=phases,
             inlet=1.175,depth=hi,elapsed=elapsed,arc=length,order_here=float(chi[-1]),start_segment_order=float(chi[0]),
             mean_transfer=state.mean(axis=0).tolist(),rms_transfer=np.sqrt((state**2).mean(axis=0)).tolist(),
             fixed_speed_position_gain_range=[float(state[:,0,0].min()),float(state[:,0,0].max())],
             max_phase_area_error=det_error)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--substeps',type=int,required=True)
    p.add_argument('--phases',type=int,default=96)
    args=p.parse_args()
    transfer(args.mesh,args.substeps,args.phases)
