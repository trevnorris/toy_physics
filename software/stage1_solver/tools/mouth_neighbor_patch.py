"""Interrupt a weak shear in a patch of the archived 0PN draining material.

Three mobile features couple by elastic and viscous neighbor forces. Order
changes the elastic coefficients and receives the conjugate elastic energy.
The calibration instrument uses only the intact far-support reaction. The
return has one controlled gap, equal preload and equal time after rebirth.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh
from scipy.sparse.linalg import splu

from mouth_contact_age import PRIOR
from mouth_free_contact import Carrier
from mouth_free_order_diffuse import NU, SIGMA, Order, emit, read
from mouth_native_order_contrast import flat_window


GAPS = (3/16, 3/4)
HOLD = .5
AMPLITUDE = .002
QUARTER = .25
LINK = 1.25


class Neighbors:
    def __init__(self):
        self.edges = np.array([[1., -1., 0.], [0., 1., -1.], [0., 0., 1.]])
        self.springs = np.array([2., .8, 3.])
        self.drags = np.array([.5, 1.4, .35])
        self.modes = {}
        for connected, drag_change in ((True, 1), (False, 1), (False, 4)):
            drag = self.drags.copy()
            drag[0] *= drag_change
            damping = self.edges.T @ (drag[:, None]*self.edges)
            spring = self.edges.T @ (self.springs[:, None]*self.edges)
            if connected:
                spring[0,0] += LINK
            rates, vec = eigh(spring, damping)
            self.modes[connected, drag_change] = spring, damping, rates, vec

    def motion(self, y, drive, connected, change, scale, dt):
        spring, _, rates, vec = self.modes[connected, change]
        source = -spring @ y
        if connected:
            source[0] += LINK*drive
        return y + vec @ ((-np.expm1(-scale*rates*dt)/rates)*(vec.T @ source))

    def energy(self, y, drive, connected):
        return .5*(self.springs @ (self.edges @ y)**2 + connected*LINK*(drive-y[0])**2)

    def reactions(self, y, drive, connected, change, scale):
        spring, damping, _, _ = self.modes[connected, change]
        force = -spring @ y
        if connected:
            force[0] += LINK*drive
        velocity = np.linalg.solve(damping, scale*force)
        far = scale*self.springs[-1]*y[-1] + self.drags[-1]*velocity[-1]
        return far, scale*connected*LINK*(drive-y[0])


class Patch:
    def __init__(self, case):
        self.case = case
        carrier = Carrier(Path(case['carrier_file']).parent)
        self.dom = carrier.domain
        self.point = case['r'], case['w']
        field = Path(case['order_file'])
        if sha256(field.read_bytes()).hexdigest() != case['order_sha256']:
            raise ValueError('Archived order hash mismatch')
        with np.load(field) as stored:
            native = stored['chi']
        psi, pr, pw = read(self.dom, carrier.field, *self.point)
        self.velocity = np.array([(pr/psi).imag, (pw/psi).imag])
        order = Order(self.dom, carrier.field, .04)
        self.interp = order.Nf
        self.mass = order.mass.tocsr()
        residual, self.jac = order.measure(native)
        if np.max(abs(residual))>3e-10:
            raise ValueError('Unforced reference is not stationary')
        self.q_base = self.dom.N @ native
        r,w = (self.dom.N @ self.dom.points).T
        shape = (flat_window(r-self.point[0],.035,.07)
                 *flat_window(w-self.point[1],.025,.055))
        physical = self.dom.weight * shape**2
        physical /= physical.sum()
        self.average = self.interp.T @ physical
        self.base_scale = (self.q_base @ physical)/QUARTER
        mobility = NU/(6*SIGMA*order.delta)
        self.shear_source = self.interp.T @ (order.wt * mobility*shape**2/QUARTER)
        self.well_coeff = order.wt*order.scale
        self.factors = {}
        self.network = Neighbors()
        self.work = np.zeros_like(native)
        emit(test='patch', mesh=case['mesh'], intake=case['intake'], state=case['state'],
             center=self.point, center_order=read(self.dom,native,*self.point)[0],
             patch_average=self.base_scale*QUARTER, mobility=mobility,
             local_dofs=len(native), source=case['carrier_file'], order=str(field),
             stationary_residual=max(abs(residual)), spring=self.network.springs.tolist(),
             drag=self.network.drags.tolist(), contact_spring=LINK)

    def from_rest(self):
        return np.zeros(3), np.zeros(self.mass.shape[0])

    def scale(self, z):
        return self.base_scale + (self.average @ z)/QUARTER

    def at_center(self, z):
        self.work[self.dom.free] = z
        value, dr, dw = read(self.dom, self.work, *self.point)
        return value, self.velocity @ [dr,dw]

    def step(self, state, drive, connected, change, dt):
        y,z = state
        if dt not in self.factors:
            self.factors[dt] = splu((self.mass+.5*dt*self.jac).tocsc())
        scale = self.scale(z)
        if scale <= 0:
            raise ValueError('Loss of positive patch-average elasticity')
        mid = self.network.motion(y,drive,connected,change,scale,.5*dt)
        energy = self.network.energy(mid,drive,connected)
        # The exact change of the native cubic order residual beyond its
        # tangent. Use an explicit remainder and an implicit centered tangent.
        dz_q = self.interp @ z
        remainder = self.interp.T @ (self.well_coeff*dz_q**2*(6*self.q_base-3+2*dz_q))
        dz = self.factors[dt].solve(-dt*(self.jac@z + remainder + self.shear_source*energy))
        new_y = self.network.motion(y,drive,connected,change,.5*(scale+self.scale(z+dz)),dt)
        return new_y,z+dz

    def advance(self, state, drive, connected, change, dt, duration, sample=False):
        count = round(duration/dt)
        if abs(count*dt-duration)>1.e-12:
            raise ValueError('Timing off grid')
        forces = []
        for _ in range(count):
            state = self.step(state,drive,connected,change,dt)
            if sample:
                forces.append(self.network.reactions(state[0],drive,connected,change,
                                                      self.scale(state[1]))[0])
        return state, forces


def predict(kernel, dt, gap):
    # This routine receives a measured intact curve, the clock and the
    # external schedule only. Release selects the effective drive that gives
    # zero absolute port traction. Freeze it at stress-free reattachment.
    sample_dt = 1/4096
    knots = np.arange(len(kernel))*dt
    end = round((2*HOLD+gap)/sample_dt)
    g = PchipInterpolator(knots,kernel)(np.arange(end+1)*sample_dt)
    jumps = np.zeros(end+1)
    jumps[0] = 1.
    first, last = round(HOLD/sample_dt), round((HOLD+gap)/sample_dt)
    for j in range(first,last+1):
        jumps[j] = -(g[j] + g[1:j-first+1][::-1] @ jumps[first:j])/g[0]
    return float(jumps @ g[::-1]), float(jumps.sum())


def predict_ramps(kernel, dt, gap, divisions):
    """Average the intact kernel over each linear motion of the free port."""
    interp = PchipInterpolator(np.arange(len(kernel))*dt, kernel)
    primitive = interp.antiderivative()
    release = -interp(HOLD)/kernel[0]
    t = np.arange(round(gap*divisions)+1)/divisions
    average = np.diff(primitive(t))*divisions
    increments = np.zeros(len(t)-1)
    for j in range(len(increments)):
        previous = average[1:j+1][::-1] @ increments[:j]
        remaining = interp(HOLD+t[j+1]) + release*interp(t[j+1]) + previous
        increments[j] = -remaining/average[0]
    end = 2*HOLD+gap
    influence = -np.diff(primitive(end-HOLD-t))*divisions
    final = interp(end) + release*interp(end-HOLD) + influence @ increments
    return float(final), float(1.+release+increments.sum())


def three_reads(patch, state, drive, h, window, pulse_fraction):
    initial, drift0 = patch.at_center(state[1])
    pulse = abs(h)*pulse_fraction
    rates, increments = [], []
    for s in (-pulse,0.,pulse):
        moved,_ = patch.advance(state,drive+s,True,1,window/2,window)
        q, drift = patch.at_center(moved[1])
        # Baseline stationary material drift is common to all three runs.
        increments.append(q-initial)
        rates.append((q-initial)/window+.5*(drift+drift0))
    neg,center,pos = rates
    curvature = pos+neg-2*center
    if abs(curvature)<1.e-14:
        raise ValueError('Unresolved order curvature')
    strain = .5*pulse*(pos-neg)/curvature
    return strain, curvature, increments, rates


def run(mesh, steps, amplitude, extra_reads, ramp_check):
    rows = json.loads(PRIOR.read_text())['runs']['mouth_mixed_order_ripple02']['stdout_lines']
    case, = [c for c in map(json.loads,rows) if c['test']=='archive_patch'
             and (c['mesh'],c['profile'],c['state'])==(mesh,'cubic','near')]
    patch = Patch(case)
    h = amplitude
    dt = 1/steps
    intact = patch.network.reactions(np.zeros(3),h,True,1,patch.base_scale)[0]
    _, curve = patch.advance(patch.from_rest(),h,True,1,dt,2*HOLD+max(GAPS),True)
    kernel = np.r_[intact,curve]/h
    predictions = {gap: predict(kernel,dt,gap) for gap in GAPS}
    if ramp_check:
        for divisions in (128,512,2048):
            predictions = {gap: predict_ramps(kernel,dt,gap,divisions) for gap in GAPS}
            emit(test='ramp_predictions', inverse_divisions=divisions, gaps=list(GAPS),
                 normalized_predictions=[predictions[g][0]/kernel[0] for g in GAPS],
                 predicted_reconnection_displacement_over_h=[predictions[g][1] for g in GAPS])
    emit(test='intact_only_prediction', steps=steps, h=h, unit_step_force=kernel.tolist(),
         normalized_predictions=[predictions[g][0]/kernel[0] for g in GAPS],
         predicted_reconnection_displacement_over_h=[predictions[g][1] for g in GAPS])
    loaded,_ = patch.advance(patch.from_rest(),h,True,1,dt,HOLD)
    for gap in GAPS:
        for change in (1,4):
            free,_ = patch.advance(loaded,0.,False,change,dt,gap)
            target = free[0][0]
            at_birth = patch.network.reactions(free[0],target,True,1,patch.scale(free[1]))
            state,_ = patch.advance(free,target,True,1,dt,HOLD)
            y,z = state
            scale = patch.scale(z)
            far,contact = patch.network.reactions(y,target,True,1,scale)
            emit(test='reconnected', steps=steps, h=h, gap=gap, detached_drag_multiple=change,
                 displacements_at_break=(loaded[0]/h).tolist(),
                 displacements_at_birth=(free[0]/h).tolist(), displacements_at_read=(y/h).tolist(),
                 force_at_birth=list(at_birth), force_far=far, force_contact=contact,
                 gain=far/intact, predicted=predictions[gap][0]/kernel[0],
                 center_order_change=patch.at_center(z)[0],
                 maximum_nodal_order_change=max(abs(z)), patch_stiffness_scale=scale)
            if max(abs(np.asarray(at_birth)))>1.e-11 or abs(far-contact)>1.e-11:
                raise ArithmeticError('Stress-free birth or instantaneous force balance failed')
            for window, fraction in ((1/512,.25),(1/2048,.25),(1/8192,.25),(1/8192,.125)) if extra_reads else ((1/8192,.25),):
                strain,curvature,changes,rates = three_reads(patch,state,target,h,window,fraction)
                # Convert the measured strain to force with the intact
                # calibration and the independently recorded change of order.
                estimated = kernel[0]*(scale/patch.base_scale)*strain
                emit(test='order_triplet', gap=gap, detached_drag_multiple=change,
                     h=h, window=window, probe=abs(h)*fraction, order_increments=changes,
                     material_rates_relative_to_stationary_drift=rates,
                     rate_curvature=curvature, inferred_strain=strain,
                     order_gain=estimated/intact, direct_gain=far/intact)
                if window==1/8192 and abs((estimated-far)/intact)>.002:
                    raise ArithmeticError('Fine readout unresolved relative to the direct force')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(34,48),required=True)
    p.add_argument('--steps',type=int,choices=(32,64,128),required=True)
    p.add_argument('--amplitude',type=float,default=AMPLITUDE)
    p.add_argument('--extra-reads',action='store_true')
    p.add_argument('--ramp-check',action='store_true')
    a=p.parse_args()
    run(a.mesh,a.steps,a.amplitude,a.extra_reads,a.ramp_check)
