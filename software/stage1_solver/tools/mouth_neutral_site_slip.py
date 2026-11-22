"""Cold longitudinal impulse on one upper-layer constituent of the supported bilayer.

The other N-1 equal-mass sites take the opposite total impulse. Both their
centre of mass and the tagged constituent are read in the laboratory frame.
Only the normal bracing and its preload depend on the uniform order. The
length and specific-energy conversion is the trial equal-site-mass assignment
used for bound-site work. This experiment contains no carrier velocity port.
"""
import argparse

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from scipy.optimize import brentq

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH as MOUTH_LENGTH
from neutral_patch_open_w import OpenCell
from neutral_patch_sites import GAP, WITHIN, emit
from neutral_patch_structural_endpoints import preparation


TRIAL_RATE = 5.5


def supported_state(cell, chi):
    signs = np.repeat([-.5, .5], cell.n//2)
    xy = preparation(cell, 0)
    ordered = np.column_stack((xy, signs*GAP))
    preload = cell.evaluate(ordered)[1][:, 2].copy()
    spring = 2*WITHIN[0]*WITHIN[1]**2

    def gradient(pos):
        grad = cell.evaluate(pos)[1]
        grad[:, 2] += chi*(spring*(pos[:, 2]-ordered[:, 2])-preload)
        return grad

    # The registered, spatially uniform sheet is an invariant symmetry sector.
    # Solve its normal force, then check its full gradient and unconstrained Hessian.
    def gap_force(gap):
        return signs @ gradient(np.column_stack((xy, signs*gap)))[:, 2]

    gap = brentq(gap_force, .6*GAP, 1.01*GAP, xtol=1.e-15)
    pos = np.column_stack((xy, signs*gap))
    hessian = cell.evaluate(pos, tangent=True)[0]
    normal = np.arange(2, 3*cell.n, 3)
    hessian[normal, normal] += chi*spring
    return pos, gradient, hessian, gap, spring


def first_reversal(frequencies, weights, end):
    def velocity(t):
        return weights @ np.cos(frequencies*t)

    grid = np.linspace(0., end, 1+int(np.ceil(end*frequencies.max()*24/np.pi)))
    values = np.cos(np.outer(grid, frequencies)) @ weights
    crossings = np.flatnonzero(values < 0.)
    last = int(crossings[0]) - 1 if len(crossings) else len(grid) - 1
    # Bound departure from each positive chord by max |v''|*dt^2/8.
    curvature = weights @ frequencies**2
    if np.any(np.minimum(values[:last], values[1:last+1]) <= curvature*np.diff(grid[:last+1])**2/8):
        raise ArithmeticError('First zero not resolved by the sampling grid')
    if not len(crossings):
        return None
    start, stop = grid[last:last+2]
    derivative = -weights @ (frequencies*np.sin(frequencies*start))
    if derivative + curvature*(stop-start) >= 0.:
        raise ArithmeticError('Crossing bracket is not strictly descending')
    return brentq(velocity, start, stop, xtol=1.e-14)


def force_trajectory(gradient, pos, readout, kick_shape, amplitude, times):
    dim = pos.size

    def rhs(t, state):
        grad = gradient(pos+state[:dim].reshape(pos.shape))
        return np.r_[state[dim:], -grad.ravel()]

    def reversal(t, state):
        return (readout @ state[dim:])/amplitude

    reversal.direction = -1
    initial = np.r_[np.zeros(dim), amplitude*kick_shape]
    sol = solve_ivp(rhs, (0., times[-1]), initial, t_eval=times, events=reversal,
                    method='DOP853', rtol=2.e-9, atol=abs(amplitude)*5.e-10, max_step=.02)
    if not sol.success:
        raise ArithmeticError(sol.message)
    displacement, velocity = np.split(sol.y, 2)
    momentum = np.max(abs(velocity[::3].sum(axis=0)))/abs(amplitude)
    return (readout @ displacement)/amplitude, (readout @ velocity)/amplitude, sol.t_events[0], momentum


def measure(cell, chi, label, clock, transit):
    pos, gradient, hessian, gap, spring = supported_state(cell, chi)
    tag = cell.n//2
    readout = np.zeros(pos.size)
    readout[::3] = -1./(cell.n-1)
    readout[3*tag] = 1.
    kick_shape = readout/(readout @ readout)
    translations = np.tile(np.eye(3)[:, :2], (cell.n, 1))/np.sqrt(cell.n)
    eigenvalues, modes = eigh(hessian)
    residual = np.max(abs(gradient(pos)))
    # Full mixed-coordinate derivative check uses the microscopic gradient.
    probe = np.random.default_rng(779).normal(size=pos.shape)
    probe /= np.max(abs(probe))
    # The |z|^3 term in the mixing law is C² at a coplanar bond. Its
    # central force difference has a linear step error even on a perfect sheet.
    differences = [(gradient(pos+step*probe)-gradient(pos-step*probe)).ravel()/(2*step)
                   for step in (6.e-6, 3.e-6)]
    finite = 2*differences[1]-differences[0]
    fd_error = np.max(abs(finite-hessian @ probe.ravel()))
    weights = (modes.T @ readout)**2/(readout @ readout)
    emit(test='supported_equilibrium', order_name=label, order=chi, sites=cell.n,
         tag=tag, gap=gap, normal_spring=chi*spring, max_force=residual,
         least_eigenvalues=eigenvalues[:5].tolist(), hessian_gradient_check_error=fd_error,
         lateral_translation_residual=np.max(abs(hessian @ translations)),
         translation_weight=weights[:2].sum(), internal_weight=weights[2:].sum(),
         projected_norm_squared=readout @ readout, total_x_momentum_per_unit_kick=kick_shape[::3].sum())
    if residual > 1.e-9 or fd_error > 2.e-5 or max(abs(eigenvalues[:2])) > 1.e-8 or eigenvalues[2] < 1.e-7:
        raise ArithmeticError('Unresolved supported equilibrium or full-site tangent')
    np.testing.assert_allclose(hessian @ translations, 0., atol=1.e-10)
    np.testing.assert_allclose([weights[:2].sum(), weights[2:].sum(), kick_shape[::3].sum(),
                                readout @ kick_shape], [0., 1., 0., 1.], atol=1.e-12, rtol=0.)
    freq, weights = np.sqrt(eigenvalues[2:]), weights[2:]

    def modal_response(native_times):
        phase = np.outer(native_times, freq)
        return np.sin(phase) @ (weights/freq), np.cos(phase) @ weights

    zero = first_reversal(freq, weights, transit/clock)
    early = np.array([0., .05, .1, .15, .2, .25])
    times = np.unique(np.r_[early*clock, .05, .1, 1/TRIAL_RATE, 2/TRIAL_RATE, transit])
    displacements, velocities = modal_response(times/clock)
    emit(test='linear_response', order_name=label, model_times=times.tolist(),
         velocity_over_initial_slip=velocities.tolist(),
         displacement_over_initial_model_velocity=(displacements*clock).tolist(),
         initial_second_velocity_derivative_model=-(weights @ freq**2)/clock**2,
         first_velocity_reversal_model_time=None if zero is None else zero*clock,
         constant_drag_comparison=np.exp(-TRIAL_RATE*times).tolist())

    ref_x, ref_v = modal_response(early)
    for amplitude in (1.e-4, -1.e-4, .5e-4, -.5e-4):
        x, v, events, momentum = force_trajectory(gradient, pos, readout, kick_shape, amplitude, early)
        emit(test='nonlinear_signed_kick', order_name=label, native_initial_relative_speed=amplitude,
             model_initial_relative_speed=amplitude*MOUTH_LENGTH/clock,
             model_times=(early*clock).tolist(), relative_x_over_initial_model_velocity=(x*clock).tolist(),
             relative_velocity_over_initial=v.tolist(), max_tangent_velocity_difference=np.max(abs(v-ref_v)),
             first_velocity_reversal_model_time=events[0]*clock if len(events) else None,
             max_net_x_momentum_over_initial_relative_speed=momentum)
        np.testing.assert_allclose(v, ref_v, atol=5.e-4, rtol=0.)
        np.testing.assert_allclose(x, ref_x, atol=5.e-5, rtol=0.)
        if momentum > 2.e-7 or bool(len(events)) != (zero is not None and zero <= early[-1]):
            raise ArithmeticError('Momentum or first-reversal check failed')
        if len(events):
            np.testing.assert_allclose(events[0], zero, atol=5.e-5, rtol=0.)


def run(args):
    patch, _, _ = read_inputs(args.mesh)
    trace, intake = cubic_trace(args.mesh)
    entry = trace.sample(trace.chi, np.array([1.175]))[0].item()
    transit = trace.transit(1.175, trace.crossing(.5), 8)
    clock = MOUTH_LENGTH/np.sqrt(2*patch['shared_scale'])
    cell = OpenCell(args.side, args.cutoff, args.power)
    emit(test='setup', side=args.side, cutoff=args.cutoff, mixing_power=args.power, mesh=args.mesh,
         intake=intake, microscopic_length_in_mouth_radii=MOUTH_LENGTH,
         specific_energy_per_native_unit=2*patch['shared_scale'], model_time_per_native_time=clock,
         trial_rate=TRIAL_RATE, comparison_relaxation_time=1/TRIAL_RATE, carrier_transit=transit,
         mass_assumption='equal site inertia with native well energy / site mass = 2*shared_scale',
         tagged_object='first upper-layer host constituent relative to the other N-1 host sites')
    for name, order in [('upstream', entry), ('half', .5)]:
        measure(cell, order, name, clock, transit)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(26, 34), default=34)
    parser.add_argument('--side', type=int, choices=(4, 6, 8), default=4)
    parser.add_argument('--cutoff', type=float, default=3.5)
    parser.add_argument('--power', type=int, choices=(1, 2), default=1)
    run(parser.parse_args())
