"""Noise-free local response of two mobile layers at fixed archived order.

Both x means are held at their archived values during preparation. The only
other constraint is a common y origin. Upper/lower translations are +/-q/2.
The material force is the sum of interplane forces on the lower layer.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from scipy.optimize import minimize, root

from neutral_patch_sites import COUNT, GAP, LAWS, PITCH_UNIT, Patch, morse


ARCHIVE = Path('research/pde_ledger_v3/_measurements/neutral_patch_conversion.json')
HISTORIES = Path('_scratch/s11c/conversion_pair_histories')
LEVELS = (.95, .50, .10, .01)


def emit(test, **fields):
    print(json.dumps(dict(test=test, **fields), allow_nan=False), flush=True)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coordinates():
    fixed = np.zeros((3, 2*COUNT, 2))
    fixed[0, :COUNT, 0] = 1/COUNT
    fixed[1, COUNT:, 0] = 1/COUNT
    fixed[2, :, 1] = 1/(2*COUNT)
    drive = np.zeros((2*COUNT, 2))
    drive[:COUNT, 0] = -.5
    drive[COUNT:, 0] = .5
    fixed = fixed.reshape(3, -1)
    return fixed, null_space(fixed), drive.ravel()


def hessian(patch, flat):
    """Pair Hessian g I + 2 (dg/d|d|²) d d^T for pair gradient g d."""
    xy = flat.reshape(-1, 2)
    delta = xy[patch.j]-xy[patch.i]
    tangent2 = np.einsum('ij,ij->i', delta, delta)
    radius2 = tangent2.copy()
    # a = d(radius²)/d(tangent²)
    a = np.ones_like(tangent2)
    da = np.zeros_like(a)
    if patch.law == 'layered':
        radius2[patch.cross] = tangent2[patch.cross]**2/PITCH_UNIT**2
        a[patch.cross] = 2*tangent2[patch.cross]/PITCH_UNIT**2
        da[patch.cross] = 2/PITCH_UNIT**2
    radius2[patch.cross] += GAP**2
    r = np.sqrt(radius2)
    exponential = np.exp(-patch.alpha*(r-patch.length))
    first = morse(r, patch.depth, patch.alpha, patch.length)[1]
    second = 2*patch.depth*patch.alpha**2*exponential*(2*exponential-1)
    g = a*first/r
    twice_dg = 2*da*first/r + (a/r)**2*(second-first/r)
    block = (g[:, None, None]*np.eye(2) +
             twice_dg[:, None, None]*delta[:, :, None]*delta[:, None, :])
    result = np.zeros((2*COUNT, 2*COUNT, 2, 2))
    for i, j, sign in ((patch.i, patch.i, 1), (patch.j, patch.j, 1),
                       (patch.i, patch.j, -1), (patch.j, patch.i, -1)):
        np.add.at(result, (i, j), sign*block)
    return result.transpose(0, 2, 1, 3).reshape(len(flat), len(flat))


def rest(patch, initial, free):
    def energy_grad(z):
        points = (initial+free@z).reshape(-1, 2)
        energy, gradient, _ = patch.evaluate(points)
        return float(energy.sum()), free.T@gradient.ravel()

    def jacobian(z):
        return free.T@hessian(patch, initial+free@z)@free

    opt = minimize(energy_grad, np.zeros(free.shape[1]), method='L-BFGS-B', jac=True,
                   options=dict(gtol=1e-10, ftol=1e-15, maxiter=30000, maxls=70))
    # The energy often stops changing in double precision before the forces
    # vanish. Solve the remaining projected gradient near the local minimum.
    polished = root(lambda z: energy_grad(z)[1], opt.x, jac=jacobian, tol=1e-11)
    energy, grad = energy_grad(polished.x)
    residual = float(max(abs(grad)))
    if residual > 1e-9 or energy > opt.fun+1e-10:
        raise ValueError(f'Unrelaxed preparation: {opt.message}, {polished.message}, {residual}')
    return initial+free@polished.x, residual


def direct_force(patch, flat, free, drive, t, amplitude, pulse_duration=None):
    """Constrained nonlinear drift, direct upper-layer x energy-gradient gauge."""
    def q(s):
        if pulse_duration is None:
            return amplitude
        return amplitude*np.sin(np.pi*s/pulse_duration)**2 if s < pulse_duration else 0.

    def points(s, z):
        return flat + free@z + drive*q(s)

    def drift(s, z):
        grad = patch.evaluate(points(s, z).reshape(-1, 2))[1]
        return -free.T@grad.ravel()

    def jac(s, z):
        return -free.T@hessian(patch, points(s, z))@free

    ode = solve_ivp(drift, (0., t[-1]), np.zeros(free.shape[1]), t_eval=t,
                    method='Radau', jac=jac, rtol=2e-8, atol=2e-12,
                    max_step=(pulse_duration or t[-1])/8)
    if not ode.success:
        raise RuntimeError(ode.message)
    return np.array([patch.evaluate(points(s, z).reshape(-1, 2))[1][COUNT:, 0].sum()
                     for s, z in zip(t, ode.y.T, strict=True)])


def compare(observed, expected, baseline, t, duration):
    signal = observed-baseline
    error = signal-expected

    def relative(mask):
        total = np.trapz(signal[mask]**2, t[mask])
        return np.sqrt(np.trapz(error[mask]**2, t[mask])/total)

    return dict(full_rel_rms=relative(np.ones(len(t), dtype=bool)),
                tail_rel_rms=relative(t >= duration),
                signal_peak=max(abs(signal)), max_abs_error=max(abs(error)))


def trial(patch, initial, t_total, level, amplitude, trace_dir):
    fixed, free, drive = coordinates()
    flat, residual = rest(patch, initial, free)
    full = hessian(patch, flat)
    restricted = free.T@full@free
    rates, vec = np.linalg.eigh(restricted)
    if rates[0] <= 0:
        raise ArithmeticError(f'Nonpositive constrained curvature: {rates[0]}')
    step0 = drive@full@drive
    coupling2 = (vec.T@free.T@full@drive)**2
    weights = coupling2/rates
    plateau = step0-sum(weights)
    t = np.linspace(0, t_total, 513)
    def step(times):
        return step0 + np.expm1(-np.outer(times, rates))@weights
    kernel = step(t)
    duration = t[128]
    # Convolution with the prescribed sin² displacement; each exponential
    # convolves with q, then contributes -b² times this integral.
    s = np.minimum(t, duration)[:, None]
    decay = np.exp(-s*rates)
    omega = 2*np.pi/duration
    phi = -np.expm1(-s*rates)/rates
    int_cos = np.real((np.exp(1j*omega*s)-decay)/(rates+1j*omega))
    memory = .5*(phi-int_cos)*np.exp(-np.maximum(t[:, None]-duration, 0)*rates)
    waveform = np.where(t < duration, np.sin(np.pi*t/duration)**2, 0.)
    prediction = amplitude*(step0*waveform-memory@coupling2)

    # Finite difference of the native Patch gradient, including every column.
    delta = 2e-6
    def gradient(p):
        return patch.evaluate(p.reshape(-1, 2))[1].ravel()
    fd = np.stack([(gradient(flat+delta*col)-gradient(flat-delta*col))/(2*delta)
                   for col in np.eye(len(flat))], axis=1)
    hess_error = float(np.max(abs(full-fd))/np.max(abs(full)))
    error_mean = max(abs(fixed@(flat-initial)))
    if hess_error > 2e-8 or error_mean > 1e-12:
        raise ArithmeticError('Potential Hessian or controlled mean mismatch')

    baseline = direct_force(patch, flat, free, drive, t, 0.)
    pulse = [direct_force(patch, flat, free, drive, t, sign*amplitude, duration)
             for sign in (1., -1.)]
    held = [direct_force(patch, flat, free, drive, t, sign*amplitude)
            for sign in (1., -1.)]
    pulse_odd = (pulse[0]-pulse[1])/2
    step_odd = (held[0]-held[1])/2
    times = np.array([0., t_total/64, t_total/16, duration, t_total])
    values = step(times)
    long = rates*t_total < 1
    emit('stationary_step', law=patch.law, chi=level, residual=residual,
         archive_to_rest_rms=np.sqrt(np.mean((flat-initial)**2)),
         archived_layer_means_x=(fixed@initial)[:2].tolist(),
         fixed_constraint_max_abs_change=error_mean, hessian_max_rel_error=hess_error,
         rest_material_force=baseline[0], zero_material_max_drift=max(abs(baseline-baseline[0])),
         min_constrained_rate=rates[0], time=times.tolist(), stiffness=values.tolist(),
         retention=(values/step0).tolist(), initial_slope=-sum(coupling2),
         initial_fractional_loss_rate=sum(coupling2)/step0,
         plateau=plateau, plateau_fraction=plateau/step0,
         slow_mode_count=int(sum(long)), slow_mode_weight_fraction=sum(weights[long])/step0,
         slow_mode_weight_remaining_fraction_at_T=sum(weights[long]*np.exp(-rates[long]*t_total))/step0,
         longest_time=1/rates[0], longest_mode_weight_fraction=weights[0]/step0,
         rate=rates.tolist(), exponential_weight=weights.tolist())
    emit('stationary_validation', law=patch.law, chi=level, duration=duration,
         amplitude=amplitude, uniform_drag_coefficient=drive@drive,
         peak_uniform_pulse_drag=(drive@drive)*amplitude*np.pi/duration,
         pulse_plus=compare(pulse[0], prediction, baseline, t, duration),
         pulse_minus=compare(pulse[1], -prediction, baseline, t, duration),
         pulse_odd=compare(pulse_odd, prediction, 0., t, duration),
         step_odd=compare(step_odd, amplitude*kernel, 0., t, duration))
    if trace_dir:
        trace_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(trace_dir/f'{patch.law}_{level:.2f}.npz',
                            time=t, baseline=baseline, kernel=kernel, predicted_pulse=prediction,
                            pulse_plus=pulse[0], pulse_minus=pulse[1], step_plus=held[0],
                            step_minus=held[1], rested_xy=flat.reshape(-1, 2))


def main(args):
    prior = json.loads(ARCHIVE.read_text())['independent_unforced_first_exits']
    for law in args.laws:
        path = HISTORIES/f'{law}.npz'
        checksum = digest(path)
        if checksum != prior[law]['history_sha256']:
            raise ValueError('Archived conversion history changed')
        with np.load(path) as data:
            order, clock, xy = data['chi'], data['t'], data['xy_zero']
        indexes = [int(np.argmin(abs(order-level))) for level in LEVELS]
        np.testing.assert_allclose(order[indexes], LEVELS, rtol=0, atol=2e-10)
        emit('stationary_preparation', law=law, history_path=str(path), history_sha256=checksum,
             seed=prior[law]['seeds'][0], gates=LEVELS, archived_indices=indexes,
             archived_times=clock[indexes].tolist(), gap=GAP, duration=clock[-1]-clock[0],
             bath_temperature=0., mobility=1., sites_per_layer=COUNT)
        for level, index in zip(LEVELS, indexes, strict=True):
            patch = Patch(law)
            patch.depth[patch.cross] *= level/.95
            trial(patch, xy[index, 0].ravel(), clock[-1]-clock[0], level,
                  args.amplitude, args.trace_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--laws', nargs='+', choices=LAWS, default=LAWS)
    parser.add_argument('--amplitude', type=float, default=.001)
    parser.add_argument('--trace-dir', type=Path)
    main(parser.parse_args())
