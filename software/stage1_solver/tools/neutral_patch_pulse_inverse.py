"""Infer the backing from an isolated-mouth pulse at two internal shear gauges.

Only the early outer-gauge traction enters the fit. The measured outer-gauge
velocity is the time-dependent drive for a one-dimensional source-free
sublayer of known depth. We use the cross-sectional mean projection and a
separate modal time propagator. The other gauge, late times and half pulse
are reserved for forward checks. No boundary force or imposed law is read.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq, least_squares, minimize_scalar
from scipy.signal import lfilter


def rms(values):
    return float(np.sqrt(np.mean(np.square(values))))


def relative_error(predicted, observed):
    denominator = rms(observed)
    if denominator <= 0:
        raise ValueError('Unexcited comparison channel')
    return rms(predicted-observed)/denominator


def signed_integral(t, values):
    return float(np.sum(np.diff(t)*(values[1:]+values[:-1])/2))


def propagate(drive, dt, cells, diffusivity, density, conductance):
    """Return outer traction and inner velocity/traction, per unit area.

    g in [0,1] is the steady layer conductance times H/eta for layer H=1.
    Bottom Robin face adds its own slip length to the half-cell resistance.
    Exact modal evolution uses the linear interpolant through the outer drive.
    """
    h = 1/cells
    eta = density*diffusivity
    bottom_transfer = h*conductance/(1-conductance+h*conductance/2)
    diagonal = np.full(cells,2.)
    diagonal[0] = 1+bottom_transfer
    diagonal[-1] = 3.
    rates, vectors = eigh_tridiagonal(diagonal,-np.ones(cells-1))
    rates *= diffusivity/h**2
    forcing = (2*diffusivity/h**2)*vectors[-1]
    delta = rates*dt
    decay = np.exp(-delta)
    # Integrate exp(-rate*(dt-s)) against both linear basis functions.
    first = (-np.expm1(-delta)-delta*decay)/(rates*delta)
    second = (delta+np.expm1(-delta))/(rates*delta)
    modes = np.zeros((cells,len(drive)))
    for k in range(cells):
        increment = forcing[k]*(first[k]*drive[:-1]+second[k]*drive[1:])
        modes[k,1:] = lfilter([1.],[1.,-decay[k]],increment)
    center = cells//2
    top_stress = (2*eta/h)*(drive-vectors[-1]@modes)
    deep_velocity = ((vectors[center]+vectors[center-1])/2)@modes
    deep_stress = (eta/h)*(vectors[center]-vectors[center-1])@modes
    return top_stress,deep_velocity,deep_stress


def read_pulse(path, samples):
    run_path = path.parent/('run_'+path.stem)
    if run_path.is_dir():
        outcome = json.loads((run_path/'outcome.json').read_text())
        valid = json.loads((run_path/'limit-validation.json').read_text())['verified']
        terminal = json.loads((run_path/'stdout').read_bytes().rstrip().splitlines()[-1])
        if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not valid or
                terminal.get('sha256') != hashlib.sha256(path.read_bytes()).hexdigest()):
            raise ValueError(f'Input archive failed provenance check: {path}')
    with np.load(path,allow_pickle=False) as data:
        saved_t = data['t']
        labels = list(data['probe_cases'].astype(str))
        u = data['u']
        tau = data['tau']
        h = float(data['h'])
        nu = float(data['nu'])
        rho = float(data['rho'])
        np.testing.assert_allclose(data['z'],[-.5,-1.])
    if u.shape!=tau.shape or u.shape!=(len(saved_t),len(labels),2):
        raise ValueError(f'Invalid sensor dimensions in {path}')
    if saved_t[0]!=0 or np.any(np.diff(saved_t)<=0) or not np.isfinite([nu,rho,h]).all():
        raise ValueError(f'Invalid gauge coordinates in {path}')
    t = np.linspace(0,saved_t[-1],samples)
    def differential(i,j):
        result = []
        for source in (u,tau):
            odd = (source[:,labels.index(i)]-source[:,labels.index(j)])/2
            result.append(np.array([np.interp(t,saved_t,column) for column in odd.T]))
        return result
    v,s = differential('pulse','minus_pulse')
    half_v,half_s = differential('half_pulse','minus_half_pulse')
    if np.any(v[:,0]!=0) or np.any(half_v[:,0]!=0):
        raise ValueError('The fit requires zero initial lateral perturbation')
    return t,v,s,half_v,half_s,h,nu,rho


def fit_backing(drive, stress, dt, cells, nu, rho):
    def loss(g):
        predicted = propagate(drive,dt,cells,nu,rho,g)[0]
        return relative_error(predicted,stress)
    inside = minimize_scalar(loss,bounds=(0.,1.),method='bounded',options={'xatol':1e-10})
    candidates = [(float(inside.fun),float(inside.x)),(loss(0.),0.),(loss(1.),1.)]
    minimum,g = min(candidates)
    # One-percent RMS measurement discrimination allowance above best fit.
    # These intervals assume the input trace and other material parameters
    # exact. They are a tolerance experiment, not statistical uncertainty.
    def threshold(x):
        return loss(x)-minimum-.01
    lo = 0. if threshold(0.)<=0 else brentq(threshold,0.,g)
    hi = 1. if threshold(1.)<=0 else brentq(threshold,g,1.)
    return g,minimum,[lo,hi]


def impedance(g,nu,rho):
    return None if g==1 else float(rho*nu*g/(1-g))


def measure(path,samples):
    t,v,s,hv,hs,h,nominal_nu,rho = read_pulse(path,samples)
    dt = t[1]-t[0]
    train_len = int(np.searchsorted(t,.65*t[-1],side='right'))
    tail = slice(train_len,None)
    outer, stress = v[0,:train_len],s[0,:train_len]
    cells = round(1/h)
    if cells%2 or not np.isclose(cells*h,1.):
        raise ValueError('Sensor spacing incompatible with sublayer depth')

    g,error,band = fit_backing(outer,stress,dt,cells,nominal_nu,rho)
    prediction = propagate(v[0],dt,cells,nominal_nu,rho,g)
    half = propagate(hv[0],dt,cells,nominal_nu,rho,g)
    # Independent-mesh continuum interpretation; the first fit explicitly
    # models the finite spatial extent of the data generator's gauges.
    gc,ec,_ = fit_backing(outer,stress,dt,4*cells,nominal_nu,rho)
    static_slope = float(outer@stress/(outer@outer))
    static_tail = static_slope*v[0,tail]

    def material_residual(parameters):
        nu = float(np.exp(parameters[0]))
        predicted = propagate(outer,dt,cells,nu,rho,parameters[1])[0]
        return (predicted-stress)/rms(stress)
    starts = ((np.log(.65),.4),(np.log(1.5),.8))
    fits = [least_squares(material_residual,start,bounds=([np.log(.2),0.],[np.log(5.),1.]),
                          gtol=1e-12,xtol=1e-11,ftol=1e-12) for start in starts]
    best = min(fits,key=lambda f:f.cost)
    nu_joint,g_joint = float(np.exp(best.x[0])),float(best.x[1])
    matrix = best.jac/np.sqrt(len(stress))
    singular = np.linalg.svd(matrix,compute_uv=False)
    joint = propagate(v[0],dt,cells,nu_joint,rho,g_joint)
    momentum = np.zeros(len(t))
    flux = s[0]-s[1]
    momentum[1:] = np.cumsum(dt*(flux[1:]+flux[:-1])/2)
    pulse_peak_time = .225*t[-1]
    result = dict(test='pulse_inverse',specimen=path.stem,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  samples=samples,cells=cells,train_end=t[train_len-1],input_nu=nominal_nu,input_rho=rho,
                  gauge_conductance=g,impedance=impedance(g,nominal_nu,rho),infinite_impedance=(g==1),
                  fit_relative_rms=error,g_one_percent_rms_interval=band,
                  impedance_one_percent_rms_interval=[impedance(q,nominal_nu,rho) for q in band],
                  finer_inverse_cells=4*cells,finer_inverse_conductance=gc,finer_inverse_fit_relative_rms=ec,
                  tail_outer_traction_error=relative_error(prediction[0][tail],s[0,tail]),
                  inner_velocity_error=relative_error(prediction[1],v[1]),
                  inner_traction_error=relative_error(prediction[2],s[1]),
                  half_outer_traction_error=relative_error(half[0],hs[0]),
                  half_inner_velocity_error=relative_error(half[1],hv[1]),
                  half_inner_traction_error=relative_error(half[2],hs[1]),
                  half_top_gauge_nonlinearity=relative_error(2*hv[0],v[0]),
                  static_fitted_traction_per_speed=static_slope,static_tail_traction_error=relative_error(static_tail,s[0,tail]),
                  joint_diffusivity=nu_joint,joint_conductance=g_joint,joint_traction_fit_rms=rms(best.fun),
                  joint_relative_traction_singular_values=singular.tolist(),
                  joint_start_spread=float(np.max(np.abs(fits[0].x-fits[1].x))),
                  joint_inner_velocity_error=relative_error(joint[1],v[1]),
                  probe_velocity_peak=v.max(axis=1).tolist(),
                  probe_velocity_peak_time=t[v.argmax(axis=1)].tolist(),
                  peak_lag_after_outer_actuator=(t[v.argmax(axis=1)]-pulse_peak_time).tolist(),
                  probe_traction_max=s.max(axis=1).tolist(),probe_traction_min=s.min(axis=1).tolist(),
                  momentum_stored_between_probes_max=float(momentum.max()),momentum_stored_between_probes_final=float(momentum[-1]),
                  return_over_positive_outer_shear_impulse=signed_integral(t,np.maximum(-s[0],0))/signed_integral(t,np.maximum(s[0],0)))
    print(json.dumps(result,allow_nan=False),flush=True)
    if not best.success or not np.all(np.isfinite(singular)):
        raise ArithmeticError('Sublayer response fit did not complete')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('probes',type=Path,nargs='+')
    p.add_argument('--samples',type=int,default=1601)
    args = p.parse_args()
    if args.samples < 101:
        p.error('Insufficient pulse samples')
    for path in args.probes:
        measure(path,args.samples)
