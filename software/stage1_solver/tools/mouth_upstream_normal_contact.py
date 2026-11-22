"""Stationary tagged matter crossing a weak sideways bend of one draining mouth.

Sites move freely in the three brane coordinates at a prescribed fourth-coordinate
speed. A single radially smooth, order-weighted, slice-radius Morse well is attached
to the held shape. No tangential registry, drag force, or particle mass for the
mouth is added. Every reported holder force is per unit tagged mass flow.
"""
import argparse
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

from mouth_calibrated_flare_coating import read_inputs
from mouth_carrier_lateral_inertia import taper
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, TEMPERATURE, emit, morse


@dataclass
class Band:
    time_per_depth: np.ndarray
    radius: np.ndarray
    order: np.ndarray
    bend: np.ndarray
    bend_deriv: np.ndarray


def coefficients(trace, w, low, high, straight):
    radius = trace.geometry.radius(w)[0]
    if straight:
        radius = np.ones_like(w)*trace.geometry.radius(np.array([low]))[0][0]
    shape, derivative = taper(w,low,high)
    return Band(1/trace.speed(w),radius,trace.sample(trace.chi,w)[0],1-shape,-derivative)


class Crossing:
    def __init__(self, trace, cuts, scale, args):
        self.q = np.array([0.,2e-5,-2e-5,1e-5,-1e-5])[:,None]
        self.scale = scale
        self.law = args.law
        self.origin = args.origin
        self.amplitude = args.amplitude
        self.q *= args.amplitude
        self.fd = 1e-7
        r0,slope = trace.geometry.radius(cuts[:1])
        mu, weight = leggauss(args.angles)
        if args.radial:
            z, weights = hermgauss(args.radial)
            speed = z*np.sqrt(2*scale*TEMPERATURE)
            self.weight = np.outer(weight/2,weights/np.sqrt(np.pi)).ravel()
        else:
            speed = np.zeros(1)
            self.weight = weight/2
        direction = np.repeat(mu,len(speed))
        transverse = np.sqrt(1-direction**2)
        velocity = np.tile(speed,len(mu))
        if args.law!='straight': velocity += slope[0]*trace.speed(cuts[:1])[0]
        state = np.stack((self.origin+r0[0]*direction,r0[0]*transverse,velocity*direction,velocity*transverse))
        self.state = np.repeat(state[:,None,:],len(self.q),axis=1)
        self.inlet_v = self.mean(self.state[2])
        self.previous_mu = np.broadcast_to(direction,self.state[0].shape).copy()
        self.initial_mu = self.previous_mu.copy()
        self.active = np.ones(self.state[0].shape,dtype=bool)
        self.always = self.active.copy()
        self.exchange = np.zeros(len(self.q))
        self.transport = np.zeros(len(self.q))
        self.left = np.zeros(len(self.q))
        self.joined = np.zeros(len(self.q))
        self.impulse = np.zeros(len(self.q))
        self.conjugate = np.zeros(len(self.q))
        self.velocity_work = np.zeros(len(self.q))

    def mean(self, a):
        return a@self.weight

    def energy(self,radius,wall_radius,order):
        return self.scale*order*morse(BETWEEN[2]+(radius-wall_radius)/PITCH,*BETWEEN)[0]

    def rhs(self,state,band,i):
        pos,perp,vx,vt = state
        dx = pos-(self.origin+self.q*band.bend[i])
        radius = np.hypot(dx,perp)
        _,gradient = morse(BETWEEN[2]+(radius-band.radius[i])/PITCH,*BETWEEN)
        radial = -self.scale*band.order[i]*gradient/PITCH
        if self.law=='unbound':
            radial = np.zeros_like(radial)
        return band.time_per_depth[i]*np.stack((vx,vt,radial*dx/radius,radial*perp/radius))

    def independently_read(self,state,band,i):
        dx = state[0]-(self.origin+self.q*band.bend[i])
        # Translate the potential at fixed absolute particle positions; this is
        # the external force which counters material-on-wall recoil at that depth.
        left = self.energy(np.hypot(dx+self.fd,state[1]),band.radius[i],band.order[i])
        right = self.energy(np.hypot(dx-self.fd,state[1]),band.radius[i],band.order[i])
        force = self.mean(right-left)/(2*self.fd)*band.time_per_depth[i]
        if self.law=='unbound': force.fill(0.)
        return np.array((force,band.bend[i]*force,band.bend_deriv[i]*self.mean(state[2])))

    def advance(self, nodes, band):
        for j,h in enumerate(np.diff(nodes)):
            k = 2*j
            s = self.state
            k1 = self.rhs(s,band,k)
            mid1 = s+h/2*k1
            k2 = self.rhs(mid1,band,k+1)
            mid2 = s+h/2*k2
            k3 = self.rhs(mid2,band,k+1)
            k4 = self.rhs(s+h*k3,band,k+2)
            self.state = s+h/6*(k1+2*k2+2*k3+k4)
            reading = h/6*(self.independently_read(s,band,k)
                       +4*self.independently_read((mid1+mid2)/2,band,k+1)
                       +self.independently_read(self.state,band,k+2))
            self.impulse += reading[0]
            self.conjugate += reading[1]
            self.velocity_work += reading[2]
            dx=self.state[0]-(self.origin+self.q*band.bend[k+2])
            radius=np.hypot(dx,self.state[1])
            mu=dx/radius
            energy=morse(BETWEEN[2]+(radius-band.radius[k+2])/PITCH,*BETWEEN)[0]
            # Geometric half-depth surface neighborhood, also measurable on the
            # potential-off control. It defines proximity, not permanent binding.
            new=energy<-.5*BETWEEN[0]
            da=new.astype(float)-self.active
            self.exchange+=self.mean((mu+self.previous_mu)*.5*da)
            self.transport+=self.mean(.5*(new.astype(float)+self.active)*(mu-self.previous_mu))
            self.left+=self.mean(np.maximum(-da,0))
            self.joined+=self.mean(np.maximum(da,0))
            self.active=new
            self.always &= new
            self.previous_mu=mu

    def read(self,depth,order,shape):
        vx=self.mean(self.state[2])
        through=vx-self.inlet_v
        dipole=3*self.mean(self.active*self.previous_mu)
        all_dipole=3*self.mean(self.previous_mu)
        by_parts=shape*vx-self.velocity_work
        fractions=self.mean(self.active.astype(float))
        def slopes(values):
            return ((values[1]-values[2])/(2*self.q[1,0]),(values[3]-values[4])/(2*self.q[3,0]))
        emit(test='cut',depth=depth,chi=order,shape=shape,absolute_force_per_flux=through.tolist(),
             direct_force_per_flux=self.impulse.tolist(),conjugate_force_per_flux=self.conjugate.tolist(),
             net_holder_slopes=slopes(through),direct_holder_slopes=slopes(self.impulse),
             conjugate_holder_slopes=slopes(self.conjugate),conjugate_by_parts_slopes=slopes(by_parts),
             angular_flux_dipole_slopes=slopes(all_dipole),near_site_dipole_slopes=slopes(dipole),
             exchange_dipole_slopes=slopes(3*self.exchange),retained_transport_dipole_slopes=slopes(3*self.transport),
             zero_near_fraction=fractions[0],all_step_near_fraction=fractions.tolist(),
             never_left_fraction=self.mean(self.always.astype(float)).tolist(),
             exits_per_tag=self.left.tolist(),entries_per_tag=self.joined.tolist())
        np.testing.assert_allclose(dipole-3*self.mean(self.initial_mu),3*(self.exchange+self.transport),atol=5e-12)
        np.testing.assert_allclose(fractions,1+self.joined-self.left,atol=5e-12)
        if not np.isfinite(self.state).all(): raise ArithmeticError('Nonfinite tag')


def main(args):
    tr,rate=cubic_trace(args.mesh)
    patch,_,_=read_inputs(args.mesh)
    low,high=1.175,tr.crossing(.5)
    cuts=np.array([low,(low+high)/2,high,tr.crossing(.25),tr.crossing(.1),tr.crossing(.01)])
    scale=2*patch['shared_scale']
    obj=Crossing(tr,cuts,scale,args)
    emit(test='setup',mesh=args.mesh,radius=1.,intake=rate,law=args.law,scale=scale,ell=PITCH,
         depths=cuts.tolist(),amplitudes=obj.q.ravel().tolist(),max_dt=args.dt,angle_nodes=args.angles,
         radial_nodes=args.radial,origin=args.origin,fd_translation=obj.fd,radial_launch_sigma=np.sqrt(scale*TEMPERATURE))
    obj.read(low,tr.sample(tr.chi,np.array([low]))[0][0],0.)
    edges=np.unique(np.r_[cuts,tr.edges[(tr.edges>cuts[0])&(tr.edges<cuts[-1])]])
    for a,b in zip(edges[:-1],edges[1:],strict=True):
        n=max(1,int(np.ceil(tr.transit(a,b,8)/args.dt)))
        nodes=np.linspace(a,b,n+1)
        samples=np.linspace(a,b,2*n+1)
        samples[0]=np.nextafter(a,b)
        samples[-1]=np.nextafter(b,a)
        band=coefficients(tr,samples,low,high,args.law=='straight')
        obj.advance(nodes,band)
        if b in cuts: obj.read(b,band.order[-1],band.bend[-1])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=[26,34],required=True)
    p.add_argument('--dt',type=float,required=True,help='Integration step in model time')
    p.add_argument('--angles',type=int,required=True)
    p.add_argument('--radial',type=int,default=0)
    p.add_argument('--origin',type=float,default=0.)
    p.add_argument('--amplitude',type=float,default=1.)
    p.add_argument('--law',choices=['surface','straight','unbound'],default='surface')
    main(p.parse_args())
