"""Use two archived velocity traces as the boundaries of a source-free slab.

Use the archive's cell widths, time grid and density. The viscous and elastic
solvers predict both face tractions without reading either traction. The
elastic trial law is tau = G0 chi(t) du/dz for displacement u. It assumes
uniform order across the gauge interval, zero initial strain and fixed
density. All times, lengths, and stresses below use the archived patch units.
"""
import json
from pathlib import Path

import numpy as np

from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import digest, emit
from neutral_patch_conversion import clock


ARCHIVE = Path('research/pde_ledger_v3/_measurements/neutral_patch_shear_pulse.json')
ROOT = Path('_scratch/s11c')
DRIVES = ['pulse', 'minus_pulse', 'half_pulse', 'minus_half_pulse']


def read_pulses(time, hashes):
    folder = ROOT/'pulse_probe01'
    paths = [folder/f'{label}_{sites}_fine.npz' for label in 'ABCD'
             for sites in ('central', 'layered')]
    paths.append(ROOT/'pulse_open01/deep10.npz')
    velocities, tractions, headers = [], [], []
    for path in paths:
        if digest(path) != hashes[str(path)]:
            raise ValueError(f'Pulse archive changed: {path}')
        with np.load(path) as data:
            np.testing.assert_array_equal(data['t'], time)
            np.testing.assert_array_equal(data['z'], [-.5, -1.])
            np.testing.assert_array_equal(data['probe_cases'], ['zero', *DRIVES])
            # A pre-pulse reference is already removed for each seed in the
            # archive. Keep all four signed traces, before taking any ratios.
            for name, output in (('u', velocities), ('tau', tractions)):
                values = data[name]
                np.testing.assert_array_equal(values[:, 0], 0.)
                np.testing.assert_array_equal(values[0], 0.)
                output.append(values[:, 1:, :])
            headers.append(tuple(float(data[key]) for key in ('h', 'rho', 'nu')))
    if len(set(headers)) != 1:
        raise ValueError('The reconstructions need a shared grid, density and viscosity')
    return paths, *headers[0], np.concatenate(velocities, axis=1), np.concatenate(tractions, axis=1)


def integrate_left(values, dt):
    return np.concatenate([np.zeros_like(values[:1]),
                           np.cumsum(dt[:, None]*values[:-1], axis=0)])


def grad_at_faces(cells, boundary, h):
    # Cells increase toward the near plane, while the files store near first.
    faces = np.empty((len(cells), cells.shape[1]+1))
    faces[:, 0] = 2*(cells[:, 0]-boundary[:, 1])/h
    faces[:, 1:-1] = np.diff(cells, axis=1)/h
    faces[:, -1] = 2*(boundary[:, 0]-cells[:, -1])/h
    return faces


def propagate(velocity, dt, h, gap, rho, coefficient, elastic=False):
    """Supply both face velocities. Record local stresses before each update."""
    count = int(round(gap/h))
    if count*h != gap or count < 2:
        raise ValueError('Gauge faces must delimit whole cells')
    state = np.zeros((velocity.shape[1], count))
    cell_speed = np.zeros_like(state)
    face_displacement = np.zeros_like(velocity[0])
    raw = np.empty_like(velocity)
    mean = np.empty(velocity.shape[:-1])
    for i in range(len(velocity)):
        boundary = face_displacement if elastic else velocity[i]
        faces = coefficient[i]*grad_at_faces(state, boundary, h)
        raw[i] = faces[:, [-1, 0]]
        mean[i] = (.5*(faces[:, 0]+faces[:, -1])+faces[:, 1:-1].sum(axis=1))/count
        if i == len(dt):
            break
        increment = (dt[i]/(rho*h))*np.diff(faces, axis=1)
        if elastic:
            cell_speed += increment
            state += dt[i]*cell_speed
            face_displacement += dt[i]*velocity[i]
        else:
            state += increment
    return raw, mean


def norm(values, weights):
    return np.sqrt(np.einsum('t,t...->...', weights, values**2))


def rel_rms(observed, predicted, weights):
    return norm(predicted-observed, weights)/norm(observed, weights)


def profile_mixture(rate, spring, stress, weights, discrepancy):
    """Common coefficients with equal trace weight and a relative L2 tolerance."""
    scale = norm(stress, weights)
    a, b, y = rate/scale, spring/scale, stress/scale

    def dot(x, z):
        return np.einsum('t,ti,ti->', weights, x, z)/x.shape[1]

    aa = dot(a, a)
    projection = dot(a, b)/aa
    orthogonal = b-projection*a
    sensitivity = np.sqrt(dot(orthogonal, orthogonal))
    g = dot(orthogonal, y)/sensitivity**2
    eta = dot(a, y)/aa-g*projection
    residual = y-eta*a-g*b
    floor = np.sqrt(dot(residual, residual))
    upper = g+np.sqrt(discrepancy**2-floor**2)/sensitivity
    eta_upper = dot(a, y)/aa-upper*projection
    at_upper = y-eta_upper*a-upper*b
    result = dict(tolerance=discrepancy, unconstrained_eta=eta, unconstrained_G0=g,
                  minimum_relative_rms=floor, G0_upper=upper, eta_at_upper=eta_upper,
                  reversible_term_rms_over_stress_rms_at_upper=upper*np.sqrt(dot(b, b)),
                  relative_rms_at_upper=np.sqrt(dot(at_upper, at_upper)),
                  sensitivity_after_eta_projection=sensitivity)
    if min(upper, eta_upper) <= 0:
        raise ArithmeticError('Mixed-law upper endpoint violates positivity')
    return result


def main():
    trace = WallTrace(48, .5)
    time, _, order, _, _ = clock(trace, 2)
    hashes = json.loads(ARCHIVE.read_text())['gauge_file_sha256']
    paths, h, rho, nu, velocity, measured = read_pulses(time, hashes)
    dt = np.diff(time)
    # Trapezoidal time weights for all measurements and fits. The original
    # pulse inversion reserved the last 35 percent of the transit.
    weights = np.r_[dt[0], dt[:-1]+dt[1:], dt[-1]]/2
    early = time <= time[0]+.65*(time[-1]-time[0])
    late = ~early
    gap = .5
    rate = (velocity[:, :, 0]-velocity[:, :, 1])/gap
    strain = integrate_left(rate, dt)
    spring = order[:, None]*strain
    trap_strain = integrate_left(np.vstack([.5*(rate[1:]+rate[:-1]), rate[-1]]), dt)
    mu = rho*nu
    viscous, mean = propagate(velocity, dt, h, gap, rho, np.full_like(time, mu))
    # Calibrate G0 on the A-central full positive pulse alone. Its synthetic
    # spatially averaged stress was reconstructed using the stated viscosity.
    x = spring[early, 0]
    y = mean[early, 0]
    w = weights[early]
    g0 = np.sum(w*x*y)/np.sum(w*x*x)
    if g0 <= 0 or np.max(dt)*np.sqrt(np.max(g0*order)/rho) >= h:
        raise ValueError('The early fitted modulus is nonpositive or underresolved in time')
    elastic, elastic_mean = propagate(velocity, dt, h, gap, rho, g0*order, elastic=True)
    viscous_error = rel_rms(measured, viscous, weights)
    viscous_mean_error = rel_rms(mean, mu*rate, weights)
    elastic_mean_error = rel_rms(g0*spring, elastic_mean, weights)
    emit(test='inputs', archive_sha256=digest(ARCHIVE), pulse_sha256={str(p):hashes[str(p)] for p in paths},
         order_inputs=trace.hashes, N=len(time), h=h, gap=gap, rho=rho, nu=nu, chi_start=order[0], chi_end=order[-1],
         duration=time[-1]-time[0], early_through=float(time[early][-1]), late_from=float(time[late][0]),
         drives=DRIVES, gauges=['near', 'deep'], calibration=paths[0].stem, G0_early=g0,
         rate_predicted_synthetic_mean_max_rel_rms=viscous_mean_error.max(),
         strain_predicted_elastic_mean_max_rel_rms=elastic_mean_error.max(),
         held_out_viscous_raw_max_rel_rms=viscous_error.max(),
         left_vs_trapezoid_strain_max_rel_rms=rel_rms(strain, trap_strain, weights).max())
    if max(viscous_error.max(), viscous_mean_error.max(), elastic_mean_error.max()) > 1e-8:
        raise ArithmeticError('Archive raw traction or discrete mean check failed')

    for n, path in enumerate(paths):
        block = slice(4*n, 4*n+4)
        target, basis, direct = mean[:, block], spring[:, block], measured[:, block]
        late_elastic = g0*basis[late]
        elastic_raw = elastic[:, block]
        # Positive-pulse instant carrying the largest traction incompatible in
        # sign with G0*chi*strain for every positive G0, when such a time exists.
        opposite = late & (target[:, 0]*basis[:, 0] < 0)
        witness = None
        if opposite.any():
            loc = np.flatnonzero(opposite)[np.argmax(abs(target[opposite, 0]))]
            witness = dict(t=time[loc], order=order[loc], strain=strain[loc, 4*n],
                           strain_rate=rate[loc, 4*n], synthetic_mean_stress=target[loc, 0],
                           elastic_mean_stress=g0*basis[loc, 0])
        # Even the optimal scalar for the late interval cannot repair a shape
        # mismatch. Calculate this separately from the shared predictive fit.
        ww = weights[late, None]
        best = np.maximum(0., np.sum(ww*target[late]*basis[late], axis=0)/
                          np.sum(ww*basis[late]**2, axis=0))
        emit(test='comparison', file=path.stem,
             early_mean_elastic_rel_rms=rel_rms(target[early], g0*basis[early], weights[early]).tolist(),
             late_mean_elastic_rel_rms=rel_rms(target[late], late_elastic, weights[late]).tolist(),
             late_best_positive_G0=best.tolist(),
             late_best_mean_elastic_rel_rms=rel_rms(target[late], best*basis[late], weights[late]).tolist(),
             late_raw_elastic_rel_rms=rel_rms(direct[late], elastic_raw[late], weights[late]).tolist(),
             late_raw_viscous_rel_rms=rel_rms(direct[late], viscous[late, block], weights[late]).tolist(),
             late_sign_witness=witness)
    emit(test='common_parallel_stress_resolution', **profile_mixture(rate, spring, mean, weights, .01))


if __name__ == '__main__':
    main()
