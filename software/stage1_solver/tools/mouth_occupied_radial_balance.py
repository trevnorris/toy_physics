"""Test independent stationary normal work on the occupied radius-one mouth."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import yaml
from scipy import sparse
from scipy.optimize import linprog
from scipy.sparse.linalg import eigsh
from scipy.special import expit

from mouth_carrier_lateral_response import prepare
from mouth_deep_attachment_traction import BANDS, TRAIN
from mouth_deep_support_calibration import band_work, quadrature
from mouth_flared_receiver import EXIT, ROUNDING, emit
from mouth_trapped_wave_ports import ARCHIVE
from mouth_wave_flare_work import Wave, radial_profile


def wave_readings(step, pad):
    p = yaml.safe_load(ARCHIVE.read_text())['parameters']
    wave = Wave(1., 1., step, pad, p)
    omega, u, op, _, radial_gate = wave.ground(1.)
    residual = np.linalg.norm(op @ u - omega**2*u)
    assert omega < p['mu_out'] and residual < 1e-7
    np.testing.assert_allclose(u @ u, 1.)
    z = wave.w.ravel()
    gate = wave.gate(1.)[0]
    width = p['delta_r']/2
    radial_derivative = wave.contrast*gate*(1-radial_gate)/width
    slab_force = (u.reshape(radial_gate.shape)**2 * radial_derivative).sum(axis=0)/(2*omega)
    displacements, volumes = [], []
    for lo, hi in BANDS:
        w, dw = quadrature(z, lo, hi)
        t = (w-lo)/(hi-lo)
        bump = 16*t**2*(1-t)**2
        a = radial_profile(w)[0]
        volumes.append(4*np.pi * dw @ (a**2*bump))
        # Integrate the piecewise-linear axial force density against the actual
        # band displacement, including bands narrower than three wave cells.
        i = np.searchsorted(z, w) - 1
        right = (w-z[i])/step
        weights = np.zeros_like(z)
        np.add.at(weights, i, dw*bump*(1-right)/step)
        np.add.at(weights, i+1, dw*bump*right/step)
        displacements.append(weights)
    displacements = np.array(displacements)
    wave_pressure = (displacements @ slab_force)/volumes
    # Move the radial edge of the well by the same projected test shape.
    # The axial gate and the far boundaries are fixed.
    controls = np.arange(1, len(BANDS)+1)/len(BANDS)
    shape = controls @ displacements
    a = radial_profile(wave.w)[0]
    axial_gate = (expit((wave.w-ROUNDING)/(p['delta_parallel']/2))
                  * expit((EXIT-wave.w)/(p['delta_parallel']/2)))

    def frequency(delta):
        shifted = expit((a+delta*shape-wave.r)/width)*axial_gate
        operator = wave.free + sparse.diags((p['mu_out']**2-wave.contrast*shifted).ravel())
        val = eigsh(operator, k=1, sigma=0., which='LM', v0=u, ncv=14,
                    return_eigenvectors=False, tol=2e-11)
        return np.sqrt(val[0])

    predicted = shape @ slab_force
    for eps in (.001, .0005):
        measured = (frequency(-eps)-frequency(eps))/(2*eps)
        emit(test='wave_work', eps=eps, predicted=predicted, finite_difference=measured)
        np.testing.assert_allclose(measured, predicted, rtol=2e-5)
    emit(test='wave_state', step=step, pad=pad, omega=omega, residual=residual,
         pressure_per_action=wave_pressure.tolist(), projected_volumes=volumes)
    return wave_pressure


def compare(curvature, pressure, wave_pressure, method):
    # Outward carrier and wave loads oppose uniform backing and area tension.
    design = np.column_stack((np.ones(len(BANDS)), curvature, -wave_pressure))
    scale = np.linalg.norm(design, axis=0)
    design = design/scale
    trained = np.linalg.solve(design[TRAIN], pressure[TRAIN])
    emit(test='three_band_fit', method=method, train=TRAIN,
         backing_tension_action=(trained/scale).tolist(),
         train_condition=np.linalg.cond(design[TRAIN]),
         residual_observed_minus_predicted=(pressure-design @ trained).tolist())
    for positive in (False, True):
        signs = np.vstack((design, -design))
        matrix = np.column_stack((signs, -np.ones(2*len(BANDS))))
        target = np.r_[pressure, -pressure]
        free = (None, None)
        material = (0., None) if positive else free
        fit = linprog([0, 0, 0, 1], A_ub=matrix, b_ub=target,
                      bounds=[free, material, material, (0., None)], method='highs')
        if not fit.success:
            raise ValueError(fit.message)
        dual = fit.ineqlin.marginals
        reduced = [0, 0, 0, 1] - matrix.T @ dual
        np.testing.assert_allclose(reduced[0], 0., atol=1e-9)
        assert np.min(reduced)>-1e-9 and np.max(dual)<1e-9
        if not positive:
            np.testing.assert_allclose(reduced[:3], 0., atol=1e-9)
        np.testing.assert_allclose(fit.fun, target @ dual, atol=1e-10)
        errors = pressure-design @ fit.x[:3]
        np.testing.assert_allclose(np.max(abs(errors)), fit.fun, atol=1e-9)
        emit(test='eight_band_bound', method=method, nonnegative_material=positive,
             backing_tension_action=(fit.x[:3]/scale).tolist(),
             max_error=fit.fun, residual=errors.tolist(), dual=dual.tolist())


def run(args):
    source = Path(f'_scratch/s11c/mouth_flared_cubic_n{args.mesh}')
    state = json.loads((source/'stdout').read_text().splitlines()[0])
    assert state['profile']=='cubic' and state['roof']==10. and state['feed_far']==9.
    for key in ('feed', 'bulk'):
        np.testing.assert_allclose(state[key], .25, atol=1e-8, rtol=0.)
    with np.load(source/'field.npz') as fields:
        emit(test='state_fields', keys=fields.files, intake=state['bulk'])
    dom, field = prepare(SimpleNamespace(source=source))
    curvature, surface, bulk = np.array([band_work(dom, field, *band) for band in BANDS]).T
    emit(test='carrier_work', bands=BANDS, curvature=curvature.tolist(),
         surface=surface.tolist(), volume=bulk.tolist())
    unit_wave = wave_readings(args.step, args.pad)
    for method, pressure in (('surface', surface), ('volume', bulk)):
        compare(curvature, pressure, unit_wave, method)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, required=True)
    parser.add_argument('--step', type=float, required=True)
    parser.add_argument('--pad', type=float, required=True)
    run(parser.parse_args())
