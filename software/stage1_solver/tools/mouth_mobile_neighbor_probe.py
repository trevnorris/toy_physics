"""Measure cold same-layer response with all neighbors mobile, before running a bend."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, root

from mouth_order_wave_pulse import digest, emit
from neutral_patch_conversion import NEUTRAL
from neutral_patch_sites import COUNT, WITHIN, morse

COHORT=Path('research/pde_ledger_v3/_measurements/neutral_patch_conversion.json')
TIMES=np.array([0.,.01,.0625,.25,1.,4.,8.,28.3125,28.43102768])
I,J=np.triu_indices(COUNT,1)
INCIDENCE=np.eye(COUNT)[J]-np.eye(COUNT)[I]


def energy_gradient(flat):
    delta=INCIDENCE@flat.reshape(COUNT,2)
    r=np.linalg.norm(delta,axis=1)
    energy,tension=morse(r,*WITHIN)
    return float(energy.sum()),(INCIDENCE.T@((tension/r)[:,None]*delta)).ravel()


def incremental_gradient(flat):
    """Evaluate the nonlinear gradient change without subtracting two nearly equal forces."""
    edges = INCIDENCE @ flat.reshape(COUNT, 2)
    rest_r = np.linalg.norm(edges, axis=1)
    depth, alpha, length = WITHIN
    weight = np.exp(-alpha * (rest_r - length))
    initial_side = morse(rest_r, *WITHIN)[1] / rest_r

    def change(offset):
        shift = INCIDENCE @ offset.reshape(COUNT, 2)
        radius = np.linalg.norm(edges + shift, axis=1)
        dr = np.sum(shift * (2*edges + shift), axis=1) / (radius + rest_r)
        dw = weight * np.expm1(-alpha * dr)
        dtension = 2*depth*alpha * dw * (1 - 2*weight - dw)
        dside = (dtension - initial_side*dr) / radius
        delta = edges*dside[:, None] + shift*(initial_side + dside)[:, None]
        return (INCIDENCE.T @ delta).ravel()
    return change


def curvatures(flat):
    xy=flat.reshape(COUNT,2)
    delta=INCIDENCE@xy
    r=np.linalg.norm(delta,axis=1)
    depth,alpha,rest=WITHIN
    exp=np.exp(-alpha*(r-rest))
    _,tension=morse(r,*WITHIN)
    radial=2*depth*alpha**2*exp*(2*exp-1)
    side=tension/r
    unit=delta/r[:,None]
    blocks=side[:,None,None]*np.eye(2)+(radial-side)[:,None,None]*unit[:,:,None]*unit[:,None,:]
    H=np.einsum('ei,ej,eab->iajb',INCIDENCE,INCIDENCE,blocks,optimize=True).reshape(2*COUNT,2*COUNT)
    scalar=INCIDENCE.T@(((radial+2*side)/3)[:,None]*INCIDENCE)
    return H,scalar


def rigid_frame(flat):
    xy=flat.reshape(COUNT,2)
    rel=xy-xy.mean(axis=0)
    modes=np.column_stack((np.tile([1.,0.],COUNT),np.tile([0.,1.],COUNT),
                           np.column_stack((-rel[:,1],rel[:,0])).ravel()))
    return linalg.orth(modes)


def relax(initial):
    basis=linalg.null_space(rigid_frame(initial).T)
    def objective(q):
        e,g=energy_gradient(initial+basis@q)
        return e,basis.T@g
    opt=minimize(objective,np.zeros(basis.shape[1]),jac=True,method='L-BFGS-B',
                 options={'ftol':1e-15,'gtol':1e-11,'maxiter':5000})
    def jac(q):
        return basis.T@curvatures(initial+basis@q)[0]@basis
    refined=root(lambda q:objective(q)[1],opt.x,jac=jac,tol=1e-11)
    xy=initial+basis@refined.x
    e,g=energy_gradient(xy)
    if max(abs(g))>2e-9 or e>opt.fun+1e-10:
        raise ArithmeticError('Bare upper layer did not settle')
    return xy


def relative_error(pred,actual):
    return float(np.max(abs(pred-actual))/max(1e-15,np.max(abs(actual))))


def step_response(H,ports,times):
    """Clamp the tagged site in both planar directions. Other sites can move."""
    free=np.delete(np.arange(len(H)),ports)
    rates,V=linalg.eigh(H[np.ix_(free,free)])
    if rates[0]<-2e-7:
        raise ArithmeticError('Unstable constrained coupon')
    rates=np.maximum(rates,0.)
    cross=V.T@H[np.ix_(free,ports)]
    integ=np.broadcast_to(times[:,None],(len(times),len(rates))).copy()
    live=rates>1e-8
    integ[:,live]=-np.expm1(-times[:,None]*rates[live])/rates[live]
    return H[np.ix_(ports,ports)]-np.einsum('ij,ti,ik->tjk',cross,integ,cross)


def held_direct(flat,ports,axis,amplitude):
    free=np.delete(np.arange(len(flat)),ports)
    change=incremental_gradient(flat)
    displaced=np.zeros_like(flat)
    displaced[ports[axis]]=amplitude
    def points(q):
        state=displaced.copy()
        state[free]+=q
        return state
    def drift(t,q):
        return -change(points(q))[free]
    result=solve_ivp(drift,(0.,TIMES[-1]),np.zeros(len(free)),method='BDF',
                     jac=lambda t,q:-curvatures(flat+points(q))[0][np.ix_(free,free)],
                     t_eval=TIMES,rtol=2e-9,atol=2e-13)
    if not result.success: raise RuntimeError(result.message)
    return np.array([change(points(q))[ports] for q in result.y.T])


def pulse_validation(flat,H,tags):
    """Fresh multi-site force history, including net force and torque; no clamps."""
    rng=np.random.default_rng(7103)
    force=np.zeros((COUNT,2))
    force[tags]=rng.normal(size=(len(tags),2))
    force=force.ravel()/np.linalg.norm(force)
    def profile(t): return np.sin(np.pi*t/2)**2 if t<2. else 0.
    change=incremental_gradient(flat)
    def direct(sign):
        sol=solve_ivp(lambda t,q:sign*1e-5*profile(t)*force-change(q),
                      (0.,TIMES[-1]),np.zeros_like(flat),method='BDF',t_eval=TIMES,
                      jac=lambda t,q:-curvatures(flat+q)[0],atol=1e-13,rtol=1e-9,max_step=.25)
        if not sol.success: raise RuntimeError(sol.message)
        return sol.y.T
    rates,v=linalg.eigh(H)
    rates[np.abs(rates)<1e-7]=0
    t=np.minimum(TIMES,2.)[:,None]
    integ=np.broadcast_to(t,(len(t),len(rates))).copy()
    nz=rates!=0
    integ[:,nz]=-np.expm1(-t*rates[nz])/rates[nz]
    oscillation=np.real((np.exp(1j*np.pi*t)-np.exp(-t*rates))/(rates+1j*np.pi))
    gain=.5*(integ-oscillation)*np.exp(-(TIMES[:,None]-t)*rates)
    pred=(gain*(v.T@force))@v.T
    val=(direct(1)-direct(-1))/(2e-5)
    return relative_error(pred,val)


def run(output):
    labels=np.array(json.loads(COHORT.read_text())['independent_id_extraction']['laws']['central']['strict_surviving_ids'])
    tags=labels[:,1]-COUNT
    xy=np.array(json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']['rested_xy'])[COUNT:]
    initial=xy.ravel()
    flat=relax(initial)
    change=incremental_gradient(flat)
    delta=np.random.default_rng(1081).normal(0,.001,size=len(flat))
    np.testing.assert_allclose(change(delta),energy_gradient(flat+delta)[1]-energy_gradient(flat)[1],
                               atol=2e-12,rtol=3e-11)
    analytic,_=curvatures(flat)
    # Independent weak, opposite force readings for every coordinate.
    def meter(h):
        return np.column_stack([(energy_gradient(flat+h*d)[1]-energy_gradient(flat-h*d)[1])/(2*h) for d in np.eye(len(flat))])
    measured=meter(2e-6)
    half=meter(1e-6)
    rms=relative_error(analytic,measured)
    half_error=relative_error(half,measured)
    rigid=rigid_frame(flat)
    free=linalg.null_space(rigid.T)
    H=free@(free.T@((measured+measured.T)/2)@free)@free.T
    rates,v=linalg.eigh(free.T@H@free)
    if rms>1e-7 or half_error>1e-7 or min(rates)<=1e-8: raise ArithmeticError('Unresolved material modes')
    v=np.column_stack((rigid,free@v))
    rates=np.r_[np.zeros(rigid.shape[1]),rates]
    # Previous convention for a matched frozen/mobile control: average each
    # bond curvature over three tangents first, as in the frozen-site trial.
    _,scalar=curvatures(initial)
    legacy,vec=linalg.eigh(scalar)
    if abs(legacy[0])>1e-8 or legacy[1]<=1e-8: raise ArithmeticError('Unresolved scalar control')
    legacy[0]=0.
    emit(test='coupon',ids=labels.tolist(),archive=digest(NEUTRAL),cohort=digest(COHORT),
         relaxation_max_site_shift=float(np.linalg.norm((flat-initial).reshape(-1,2),axis=1).max()),
         residual=float(max(abs(energy_gradient(flat)[1]))),old_residual=float(max(abs(energy_gradient(initial)[1]))),
         calibration_difference=rms,half_amplitude_difference=half_error,
         first_internal_rate=rates[3],fastest_rate=rates[-1],free_rigid_modes=3,
         old_frozen_scalar=scalar.diagonal()[tags].tolist(),plane_instant=(H.diagonal().reshape(-1,2)[tags].mean(axis=1)).tolist())
    readings=[]; worst=0.
    for index in tags:
        ports=2*index+np.arange(2)
        pred=step_response(H,ports,TIMES)
        direct=np.stack([(held_direct(flat,ports,a,1e-5)-held_direct(flat,ports,a,-1e-5))/(2e-5) for a in range(2)],axis=2)
        worst=max(worst,relative_error(pred,direct))
        readings.append(np.trace(direct,axis1=1,axis2=2)/2)
    readings=np.array(readings)
    err=pulse_validation(flat,H,tags)
    emit(test='mobile_holds',times=TIMES.tolist(),force_per_displacement=readings.tolist(),
         early_min_to_max=[min(readings[:,2]),max(readings[:,2])],
         transit_min_to_max=[min(readings[:,-1]),max(readings[:,-1])],
         peak_normalized_hold_error=worst,peak_normalized_pulse_error=err)
    if worst>1e-5 or err>1e-5: raise ArithmeticError('Nonlinear coupon validation failed')
    # B[o,mode,site] maps applied forces to modes in drive direction o.
    B=np.array([v[2*tags+axis].T for axis in range(2)])
    with output.open('xb') as out:
        np.savez_compressed(out,rates=rates,B=B,legacy_rates=legacy,legacy_B=vec[tags].T[None,:,:],
                            labels=labels,upper_equilibrium=flat.reshape(-1,2))
    emit(test='frozen_coupon_response',file=str(output),sha256=digest(output))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
