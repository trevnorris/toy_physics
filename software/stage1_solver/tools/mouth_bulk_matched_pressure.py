"""Relax a four-periodic-dimensional sample at zero applied pressure and test its stiffness."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm
from scipy.optimize import minimize

from mouth_bulk_isotropic_4d import strain_basis
from mouth_bulk_relaxed_shape import PeriodicMaterial, constitutive, read_stress
from neutral_patch_sites import WITHIN, emit


def free_cell(material, initial_y, initial_cell):
    """Minimize energy over sites and all ten strains in an upper triangular cell gauge."""
    length = np.linalg.det(initial_cell)**.25
    weight = np.sqrt(material.n)
    diagonal = np.diag_indices(4)
    upper = np.triu_indices(4, 1)
    # The variable cell uses logarithmic lengths for positive volume. Its
    # off-diagonal coordinates and the fractional site displacements are independent.
    def evaluate(flat, details=False):
        displacement = flat[:-10].reshape(-1, 4)
        params = flat[-10:] / weight
        y = initial_y + displacement / length
        cell = initial_cell.copy()
        cell[diagonal] *= np.exp(params[:4])
        cell[upper] += length * params[4:]
        energy, gradient, cell_gradient = material.read(y, cell)
        if details:
            return y, cell, energy, params
        shape_gradient = np.r_[cell[diagonal] * cell_gradient[diagonal], length * cell_gradient[upper]]
        return energy, np.r_[((gradient @ cell.T) / length).ravel(), shape_gradient / weight]

    size = 4 * material.n
    bounds = [(None, None)] * size + [(-.6*weight, .6*weight)] * 4 + [(-weight, weight)] * 6
    opt = minimize(evaluate, np.zeros(size + 10), jac=True, method='L-BFGS-B', bounds=bounds,
                   options=dict(maxiter=20000, maxls=70, ftol=2e-16, gtol=1e-9, maxcor=30))
    y, cell, _, params = evaluate(opt.x, True)
    if np.max(abs(params[:4])) >= .599 or np.max(abs(params[4:])) >= .999:
        raise ArithmeticError('Accepted cell touched the line-search bounds')
    if abs(opt.jac).max() > 1e-4:
        raise ArithmeticError(f'Unresolved joint site and cell equilibrium ({opt.message}, {abs(opt.jac).max()})')
    return y, cell


def diffraction(y, cell, q_max=12.):
    """Max normalized intensity for nonzero reciprocal vectors of physical length at most q_max."""
    extent = np.floor(q_max*np.linalg.norm(cell, axis=1)/(2*np.pi)).astype(int)
    integers = np.array([h for h in itertools.product(*(range(-k, k+1) for k in extent))
                         if h > (0, 0, 0, 0)])
    vectors = 2*np.pi*integers @ np.linalg.inv(cell).T
    magnitudes = np.linalg.norm(vectors, axis=1)
    use = magnitudes <= q_max
    h = integers[use]
    intensities = abs(np.exp(2j*np.pi*(h @ y.T)).mean(axis=1))**2
    best = np.argmax(intensities)
    return dict(q_max=q_max, reciprocal_pairs=len(h), peak=float(intensities[best]),
                indices=h[best].tolist(), q=float(magnitudes[use][best]))


def run(source, work):
    outcome = json.loads((source.parent.parent/'outcome.json').read_text())
    if not outcome['limitsVerified'] or outcome['exitCode'] or outcome['childOutcome']['exitCode']:
        raise ValueError('Expected a verified guarded shape-relaxation input')
    with np.load(source) as data:
        initial_y, initial_cell = data['y'], data['cell']
        seed, cutoff = int(data['seed']), float(data['rc'])
    if cutoff != 2.5:
        raise ValueError('The matched-pressure comparison uses cutoff 2.5 throughout')
    work.mkdir(parents=True, exist_ok=True)
    material = PeriodicMaterial(len(initial_y), cutoff)
    basis = strain_basis()
    before = read_stress(material, initial_y, initial_cell, basis)
    emit(test='input', seed=seed, source=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
         sites=material.n, cutoff=cutoff, well=list(WITHIN), temperature=0., applied_pressure=0.,
         initial_density=material.n/np.linalg.det(initial_cell), initial_pressure=-before[3][0],
         initial_energy_per_site=before[0]/material.n)
    y, cell = free_cell(material, initial_y, initial_cell)
    eigenvalues, eigenvectors = eigh(constitutive(material, y, cell))
    if abs(eigenvalues[:4]).max() > 1e-7 or eigenvalues[4] <= 1e-7:
        raise ArithmeticError(f'Internal Hessian unresolved: {eigenvalues[:6]}')
    internal = eigenvectors[:, 4:]

    def relax_sites(start_y, trial_cell):
        trial_y = start_y.copy()
        inv = np.linalg.inv(trial_cell)
        for _ in range(30):
            energy, grad, _, readings = read_stress(material, trial_y, trial_cell, basis)
            if abs(grad).max() < 2e-10:
                return trial_y, energy, grad, readings
            change = internal @ ((internal.T @ grad.ravel()) / eigenvalues[4:])
            trial_y -= change.reshape(-1, 4) @ inv
        raise ArithmeticError('Signed equilibrium failed to reach force tolerance')

    def tangent(step):
        columns = []
        for generator in basis:
            plus = relax_sites(y, cell @ expm(step*generator))[3]
            minus = relax_sites(y, cell @ expm(-step*generator))[3]
            columns.append((plus-minus)/(2*step))
        return np.array(columns).T

    # Remove the L-BFGS energy tolerance from the pressure match using the
    # relaxed-site stress Jacobian. Each signed strain below re-equilibrates the sites.
    y, _, _, residual = relax_sites(y, cell)
    if abs(residual).max() > 1e-9:
        correction = np.linalg.solve(tangent(1e-4), -residual)
        if abs(correction).max() > 1e-4:
            raise ArithmeticError('Cell-stress correction is outside the local tangent regime')
        cell = cell @ expm(np.einsum('q,qab->ab', correction, basis))
    y, energy, grad, residual = relax_sites(y, cell)
    if abs(residual).max() > 1e-8:
        raise ArithmeticError(f'Zero-load stress residual {residual}')
    if energy > before[0] + 1e-9:
        raise ArithmeticError('Releasing the volume raised the energy')
    emit(test='zero_load', seed=seed, density=material.n/np.linalg.det(cell),
         pressure=-residual[0], shear_residual_norm=np.linalg.norm(residual[1:]),
         max_site_force=abs(grad).max(), volume_ratio=np.linalg.det(cell)/np.linalg.det(initial_cell),
         energy_per_site=energy/material.n, cell=cell.tolist(), diffraction=diffraction(y, cell))
    matrices = []
    for step in (1e-4, 5e-5):
        raw = tangent(step)
        asymmetry = abs(raw-raw.T).max()
        if asymmetry > .01:
            raise ArithmeticError(f'Relaxed-site stress Jacobian asymmetry {asymmetry}')
        response = .5*(raw+raw.T)
        matrices.append(response)
        bulk, mixed, shear = response[0,0], response[1:,0], response[1:,1:]
        # A change in pressure permits every shear to relax; a shear at target
        # pressure permits volume to relax. Both Schur complements use the same tangent.
        shear_at_pressure = shear - np.outer(mixed, mixed)/bulk
        emit(test='signed_probe', step=step, bulk_fixed_shape=bulk,
             bulk_free_shape=bulk-mixed @ np.linalg.solve(shear, mixed),
             shear_eigenvalues_fixed_volume=np.linalg.eigvalsh(shear).tolist(),
             shear_eigenvalues_fixed_pressure=np.linalg.eigvalsh(shear_at_pressure).tolist(),
             min_full_eigenvalue=np.linalg.eigvalsh(response)[0], bulk_shear_coupling=mixed.tolist(),
             max_asymmetry=asymmetry)
    np.testing.assert_allclose(matrices[0], matrices[1], rtol=2e-4, atol=.005)
    emit(test='step_check', max_tangent_difference=abs(matrices[1]-matrices[0]).max())
    np.savez(work/'endpoint.npz', y=y, cell=cell, rc=cutoff, seed=seed, tangent=matrices[-1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--work', required=True, type=Path)
    args = parser.parse_args()
    run(args.source, args.work)
