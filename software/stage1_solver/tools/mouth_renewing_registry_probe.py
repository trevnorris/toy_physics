"""Free transverse response of fresh material traversing a small developed flare patch.

Tangential flow along the patch is supplied by the stationary carrier. Only the
transverse coordinate is free. The lining is a hypothesis and the driven axial
flow can do work. No transverse drag or assigned transverse pinning is used.
"""
import argparse

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import map_coordinates, spline_filter

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, GAP, emit, morse


class PeriodicLining:
    def __init__(self, pitch, size):
        self.basis = pitch*np.array([[1.,0.],[.5,np.sqrt(3.)/2]])
        self.inverse = np.linalg.inv(self.basis)
        self.size = size
        phases = np.moveaxis(np.indices((size,size)),0,-1)/size
        xy = phases@self.basis
        g = self.gradient(xy)
        self.table = [spline_filter(g[...,i],order=3,mode='grid-wrap') for i in (0,1)]

    def gradient(self, xy):
        # A symmetric shell wide enough for accurate periodic edge matching.
        result = np.zeros_like(xy)
        for i in range(-9,11):
            for j in range(-9,11):
                delta = xy - np.array([i,j])@self.basis
                r = np.sqrt((delta*delta).sum(axis=-1) + GAP**2)
                _,v = morse(r,*BETWEEN)
                result += (v/r)[...,None]*delta
        return result

    def projected_table(self, direction):
        return self.table[0]*direction[0]+self.table[1]*direction[1]

    def interpolate(self, table, fractions):
        coords = self.size*fractions
        return map_coordinates(table,coords,order=3,mode='grid-wrap',prefilter=False)


def path(trace,steps):
    lo,hi=1.55,1.65
    w = np.unique(np.r_[np.linspace(lo,hi,2049),trace.edges[(trace.edges>lo)&(trace.edges<hi)]])
    elapsed = cumulative_trapezoid(1/trace.speed(w),w,initial=0.)
    arc = cumulative_trapezoid(np.hypot(1.,trace.geometry.radius(w)[1]),w,initial=0.)
    time = np.linspace(0,elapsed[-1],steps+1)
    depth = PchipInterpolator(elapsed,w)(time)
    s = PchipInterpolator(elapsed,arc)(time)
    chi = trace.sample(trace.chi,depth)[0]
    if chi.min()<0 or chi.max()>1: raise ValueError('Unphysical wall order')
    np.testing.assert_allclose(time[-1],trace.transit(lo,hi,12),rtol=5e-7)
    return time,s,chi


def passage(lining, patch, time, arc, coupling, ray, supply, phase_n):
    t = np.array([1.,ray])/np.hypot(1.,ray)
    n = np.array([-t[1],t[0]])
    table = lining.projected_table(n)
    phase = (np.arange(phase_n)+.5)/phase_n
    if supply=='row':
        if ray!=0: raise ValueError('The coherent row experiment uses axial row flow')
        # Same transverse row, uniform starting longitude along a primitive period.
        initial = np.column_stack((phase,np.zeros_like(phase)))
    else:
        u,v = np.meshgrid(phase,phase,indexing='ij')
        initial = np.column_stack((u.ravel(),v.ravel()))
    # y, q are in carrier length, v in carrier length/time; displacements signed.
    q0=5.e-6
    v0=q0/time[-1]
    qs=q0*np.array([0.,1.,-1.,.5,-.5,0.,0.,0.,0.,.5])
    vin=v0*np.array([0.,0.,0.,0.,0.,1.,-1.,.5,-.5,-.5])
    y=np.zeros((len(qs),len(initial)))
    velocity=np.broadcast_to(vin[:,None],y.shape).copy()
    # Assign half the unit resting carrier volume density to the facing sheet.
    # Thus energy per unit site mass is twice the pressure-calibrated factor.
    acc_scale=2*patch['shared_scale']/PITCH
    curvature_coord=n@lining.inverse/PITCH
    along=arc[:,None]*(t@lining.inverse/PITCH)
    def accel(j):
        frac=(initial.T[:,None,:]+along[j,:,None,None]
              +curvature_coord[:,None,None]*(y-qs[:,None]))
        return -acc_scale*coupling[j]*lining.interpolate(table,frac)
    a=accel(0)
    dt=time[1]-time[0]
    for j in range(1,len(time)):
        y += dt*velocity+.5*dt**2*a
        a_next=accel(j)
        velocity += .5*dt*(a+a_next)
        a=a_next
    # The flux of fresh sites per area is the instantaneous reference inventory
    # divided by residence time. Reciprocal wall force is minus their momentum gain.
    # External holding force on the wall has the sign of the gain by the sites.
    sheet_mass_per_area=GAP*PITCH
    force=sheet_mass_per_area*np.mean(velocity-vin[:,None],axis=1)/time[-1]
    force-=force[0].copy()
    ymean=y.mean(axis=1)
    gain=(ymean[1]-ymean[2])/(2*q0)
    kout=(force[1]-force[2])/(2*q0)
    zin=(force[5]-force[6])/(2*v0)
    kval_half=(force[3]-force[4])/q0
    zval_half=(force[7]-force[8])/v0
    held_out=.5*(kout*q0-zin*v0)
    emit(test='free_passage',supply=supply,ray=ray,number_of_initial_phases=len(initial),
         concentration_min=coupling.min(),concentration_max=coupling.max(),
         relative_force=force.tolist(),wall_offsets=qs.tolist(),incoming_transverse_speeds=vin.tolist(),
         offset_slope=kout,half_offset_slope=kval_half,velocity_slope=zin,half_velocity_slope=zval_half,
         wall_offset_to_mean_exit_shift=gain,remaining_inlet_position_gain=1-gain,
         mean_exit_velocity_gain=(velocity[5].mean()-velocity[6].mean())/(2*v0),
         heldout_force=force[-1],predicted_mixed_force=held_out,
         rms_unforced_exit_speed=float(np.sqrt(np.mean(velocity[0]**2))))
    if not np.isfinite(y).all() or not np.isfinite(velocity).all():
        raise ArithmeticError('Invalid free transverse motion')


def main(args):
    patch,_,_=read_inputs(args.mesh)
    tr,intake=cubic_trace(args.mesh)
    time,arc,order=path(tr,args.steps)
    wall=PeriodicLining(patch['pitch'],args.grid)
    sample=np.array([[.38,-.42],[.7,.11],[.3,1.2],[-.6,-.3]])
    for axis in np.eye(2):
        direct=wall.gradient(sample)@axis
        grid=wall.interpolate(wall.projected_table(axis),(sample@wall.inverse).T)
        np.testing.assert_allclose(direct,grid,atol=1e-6)
    emit(test='setup',mesh=args.mesh,steps=args.steps,grid=args.grid,phases_per_axis=args.phases,
         intake=intake,patch_depth=[1.55,1.65],arc=arc[-1],residence=time[-1],
         pitches=arc[-1]/(PITCH*patch['pitch']),order_range=[order.min(),order.max()],
         half_cell_mass_per_three_area=GAP*PITCH)
    for coupling_name in args.couplings:
        coupling={'order':order,'persistent':np.ones_like(order),'smooth':np.zeros_like(order)}[coupling_name]
        emit(test='lining',law=coupling_name)
        # Fully persistent alignment is a separate single coherent row experiment.
        # A uniform transverse phase ensemble has an exact zero static offset
        # response by translation symmetry, including fully nonlinear trajectories.
        cases=[(0.,'row')]
        if coupling_name=='order':
            cases.extend([(0.,'uniform'),(np.sqrt(2.)/4,'uniform')])
        for ray,supply in cases:
            passage(wall,patch,time,arc,coupling,ray,supply,args.phases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=[26,34],required=True)
    p.add_argument('--steps',type=int,required=True)
    p.add_argument('--grid',type=int,required=True)
    p.add_argument('--phases',type=int,required=True)
    p.add_argument('--couplings',nargs='+',choices=['order','persistent','smooth'],default=['order','persistent','smooth'])
    main(p.parse_args())
