"""Reversible held-shape time control with one local order displacement and one prepared shear mode."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from mouth_live_waist_order import Order, Meridian, load_for, state, BIAS, emit
from mouth_live_waist_motion import CarrierMotion


def experiment(a):
    saved = state(a.radius,a.point)
    shape = Meridian(a.radius,saved['coordinates'])
    order = Order(shape,a.mesh)
    data = np.load(a.field)
    chi = data['order']
    pole = complex(data['pole'])
    fields = load_for(shape,a.mesh,a.field_mesh,saved['wave_frequency'],(pole,data['vector']))
    if chi.min() == 0.:
        from mouth_live_waist_active_wave import ActiveWave
        fields = ActiveWave.install(fields)
    _,_,info=fields.evaluate(chi)
    omega=info['wave_frequency']
    t = order.w[:,None]/shape.length
    x = np.linspace(0.,1.,a.mesh+1)[None,:]
    # Extend the local waist displacement smoothly, vanish on every pin.
    # Weight by chi(1-chi) to avoid pushing a saturated site through a well.
    eta = chi*(1-chi)*(1-x**2)*np.sin(np.pi*t)*np.exp(-((t-order.t_waist)/.16)**2)
    eta[~order.free]=0.
    eta/=max(eta.ravel())
    positive=eta>1e-8
    lower=max(-.25,np.max(-chi[positive]/eta[positive]))
    upper=min(.25,np.min((1-chi[positive])/eta[positive]))
    if not (lower<0<upper):
        raise ValueError('No reversible two-sided perturbation for this preparation')
    kw,mw=fields.wave.strain_weights,fields.wave.inertia_weights
    k=float(np.sum(kw*eta))
    m=float(np.sum(mw*eta))
    near=np.broadcast_to(abs(t-order.t_waist)<.1,chi.shape)
    far=np.broadcast_to(t>order.t_waist+.1,chi.shape)
    motion=CarrierMotion(fields.carrier,chi,eta)
    masses=motion.matched_mass()
    abscissa=np.linspace(lower,upper,151)
    samples=[motion.held(e) for e in abscissa]
    inertia=CubicSpline(abscissa,[s['frozen_mass'] for s in samples])
    carrier=CubicSpline(abscissa,[s['energy_change'] for s in samples])
    y=chi.ravel()-.5
    z=eta.ravel()
    w=order.lump
    coefficients=[0., 2*z@(order.stiff@y)+w@((4*y**3-y)*z),
                  z@(order.stiff@z)+w@((6*y*y-.5)*z*z),
                  w@(4*y*z**3), w@(z**4)]
    potential=np.polynomial.Polynomial(coefficients)
    def derivs(a):
        return omega**2*(1+k*a),1+m*a
    def rhs(time,vector):
        v,p,q,z=vector
        K,M=derivs(v)
        mass=inertia(v)
        accel = (p**2*inertia(v,1)/(2*mass**2)-potential.deriv()(v)-carrier(v,1)
                 +.5*m*z**2/M**2 - .5*omega**2*k*q**2)
        return np.array([p/mass,accel,z/M,-K*q])
    def energy(row):
        a,p,q,z=row
        K,M=derivs(a)
        return p*p/(2*inertia(a))+potential(a)+carrier(a)+z*z/(2*M)+K*q*q/2
    def leave_lo(t,v):
        return v[0]-lower
    def leave_hi(t,v):
        return upper-v[0]
    leave_lo.terminal=leave_hi.terminal=True
    mode0=order.axis(eta)
    tspan=2*2*np.pi/omega
    rate2=(potential.deriv(2)(0)+carrier(0,2))/inertia(0)
    emit(kind='modal_control',radius=a.radius,point=a.point,mesh=a.mesh,field_mesh=a.field_mesh,
         initial_axis_order=order.axis(chi),mode_axis_component=mode0,
         field_min=float(chi.min()), field_max=float(chi.max()),
         fractional_stiffness_derivative=k,fractional_inertia_derivative=m,
         reference_frequency=omega, approximate_bare_order_frequency=np.sign(rate2)*np.sqrt(abs(rate2)),
         mean_force_at_rest=-float(potential.deriv()(0)+carrier(0,1)+omega*(k-m)/2),
         carrier_current=info['current'],unperturbed_wave_energy_loss_rate=info['wave_loss_rate'],
         waist_depth_fraction=order.t_waist, kinetic_near_waist=np.sum((chi*mw)[near]),
         strain_near_waist=np.sum((chi*kw)[near]),kinetic_downstream=np.sum((chi*mw)[far]),
         strain_downstream=np.sum((chi*kw)[far]), allowed_range=[lower,upper],
         **masses)
    for label,initial in (('pure_strain',[0.,0.,np.sqrt(2/omega),0.]),
                          ('pure_kinetic',[0.,0.,0.,np.sqrt(2*omega)])):
        prior=None
        for relative in (2e-8,4e-10):
            init=np.asarray(initial)
            sol=solve_ivp(rhs,(0,tspan),init,method='DOP853',rtol=relative,atol=relative*1e-3,
                          dense_output=True,events=(leave_lo,leave_hi))
            if not sol.success:
                raise ValueError(sol.message)
            grid=np.linspace(0,sol.t[-1],1401)
            ys=sol.sol(grid)
            E=energy(ys)
            drift=max(abs(E-E[0]))/E[0]
            position=ys[0]
            center=order.axis(chi)+position*mode0
            upper_extent=max(chi.max(),np.max(chi+max(position)*eta))
            lower_extent=min(chi.min(),np.min(chi+min(position)*eta))
            result=dict(kind='time',phase=label,tolerance=relative,time=sol.t[-1],goal_time=tspan,
                        initial_generalized_order_acceleration=float(rhs(0,init)[1]/inertia(0)),
                        order_coordinate_min=min(position),order_coordinate_max=max(position),
                        axis_order_min=min(center),axis_order_max=max(center),
                        field_min=lower_extent,field_max=upper_extent,
                        relative_energy_drift=drift,reached_bound=sol.status==1,
                        final_state=ys[:,-1].tolist(),steps=len(sol.t)-1)
            if prior is not None:
                result['relative_trace_error']=max(abs(prior['sol'](grid)[0]-position))/max(1e-8,max(abs(position)))
                if result['relative_trace_error']>.004 or drift>1e-7:
                    emit(**result)
                    raise ArithmeticError('Order-wave transient unresolved by time refinement')
            emit(**result)
            prior={'sol':sol.sol}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--radius',type=int,choices=(1,2),required=True)
    ap.add_argument('--point',choices=('broad','mid','tight'),required=True)
    ap.add_argument('--mesh',type=int,required=True)
    ap.add_argument('--field-mesh',type=int,required=True)
    ap.add_argument('--field',type=Path,required=True)
    experiment(ap.parse_args())
