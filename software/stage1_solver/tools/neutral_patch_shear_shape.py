"""Release the fixed-area cell shear on the free-w 72-site defect preparation."""
import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from neutral_patch_open_w import OpenCell, modulus
from neutral_patch_sites import GAP, emit
from neutral_patch_structural_endpoints import Cell, preparation


def prepared_positions(power):
    # Reconstruct the specified fixed-shape endpoint as the initial condition.
    planar = Cell(6, 6.)
    xy = planar.settle(preparation(planar, 1321))[0]
    heights = GAP*np.repeat([-.5, .5], planar.n//2)
    start = np.column_stack((xy, heights))
    noise = np.random.default_rng(319).normal(size=start.shape)
    cell = OpenCell(6, 6., power)
    return cell.relax(start + .012*(noise-noise.mean(axis=0)))[0]


def wide_stencil(power, cutoff):
    # For |gamma| <= 0.5 any included physical bond has unsheared
    # planar length <= (1+|gamma|)*cutoff. Keep images out to this limit.
    cell = OpenCell(6, 1.5*cutoff, power)
    cell.cutoff = cutoff
    return cell


def tangent(cell, pos, shear, step):
    """Differentiate the unsheared-coordinate gradient and conjugate cell force."""
    n = pos.size
    hessian = np.zeros((n, n))
    for k in range(n):
        offset = np.zeros_like(pos)
        offset.flat[k] = step
        gp = cell.evaluate(pos+offset, shear)[1]
        gm = cell.evaluate(pos-offset, shear)[1]
        hessian[:, k] = ((gp-gm)/(2*step)).ravel()
    np.testing.assert_allclose(hessian, hessian.T, atol=4e-6, rtol=1e-7)
    _, gp, sp = cell.evaluate(pos, shear+step)
    _, gm, sm = cell.evaluate(pos, shear-step)
    mixed = ((gp-gm)/(2*step)).ravel()
    affine = cell.area*(sp-sm)/(2*step)
    hessian = (hessian+hessian.T)/2
    slope, response, gap = modulus(hessian, mixed, affine, cell.area, 3)
    return hessian, mixed, affine, slope, response.reshape(pos.shape), gap


def release(cell, initial):
    scale = np.sqrt(cell.n)
    def unpack(state):
        return np.vstack((initial[0], state[:-1].reshape(-1, 3))), state[-1]/scale
    def objective(state):
        pos, shear = unpack(state)
        energy, grad, stress = cell.evaluate(pos, shear)
        return energy, np.r_[grad[1:].ravel(), cell.area*stress/scale]
    state = np.r_[initial[1:].ravel(), 0.]
    bounds = [(None, None)]*(len(state)-1) + [(-.49*scale, .49*scale)]
    fit = minimize(objective, state, jac=True, method='L-BFGS-B', bounds=bounds,
                   options=dict(ftol=2e-16, gtol=1e-10, maxiter=8000, maxls=60, maxcor=20))
    pos, shear = unpack(fit.x)
    if abs(shear) > .48:
        raise ArithmeticError('Shape release reached image-stencil trust boundary')
    h, b, affine, _, _, _ = tangent(cell, pos, shear, 2e-6)
    # Newton polish on the physical shear and all internal coordinates.
    block = np.empty((pos.size-2, pos.size-2))
    block[:-1, :-1] = h[3:, 3:]
    block[-1, :-1] = block[:-1, -1] = b[3:]
    block[-1, -1] = affine
    _, g, stress = cell.evaluate(pos, shear)
    delta = np.linalg.solve(block, np.r_[g[1:].ravel(), cell.area*stress])
    if np.max(abs(delta)) > 2e-4:
        raise ArithmeticError('Equilibrium unresolved after shape optimization')
    pos[1:] -= delta[:-1].reshape(-1, 3)
    shear -= delta[-1]
    return pos, shear


def run(power, output):
    cell = wide_stencil(power, 6.)
    initial = prepared_positions(power)
    e0, g0, s0 = cell.evaluate(initial)
    h0, b0, born0 = cell.evaluate(initial, tangent=True)
    initial_slope = modulus(h0, b0, born0, cell.area, 3)[0]
    pos, shear = release(cell, initial)
    e, g, stress = cell.evaluate(pos, shear)
    _, _, _, slope, correction, gap = tangent(cell, pos, shear, 1e-6)
    correction -= correction[0].copy()
    checks = []
    for step in (1e-4, 5e-5):
        trials = [cell.relax(pos+sign*step*correction, shear+sign*step) for sign in (-1, 1)]
        _, em, sm, _ = trials[0]
        _, ep, sp, _ = trials[1]
        checks.append(dict(step=step, stress_slope=(sp-sm)/(2*step),
                           energy_curvature=(ep+em-2*e)/(step**2*cell.area)))
    np.savez(output, initial=initial, final=pos, shear=shear, power=power, box=cell.box)
    emit(test='zero_shear_load', power=power, initial_stress=s0, initial_slope=initial_slope,
         initial_residual=np.max(abs(g0)), shear=shear, energy_change_per_site=(e-e0)/cell.n,
         final_stress=stress, final_slope=slope, final_residual=np.max(abs(g)), least_mode=gap,
         checks=checks)
    if np.max(abs(g)) > 2e-8 or abs(stress) > 2e-8 or e >= e0 or slope <= 0:
        raise ArithmeticError('Unresolved stress-free local minimum')
    for check in checks:
        np.testing.assert_allclose([check['stress_slope'],check['energy_curvature']],
                                   slope, rtol=.001, atol=.004)
    ext = wide_stencil(power, 7.)
    ex, gx, sx = ext.evaluate(pos, shear)
    np.testing.assert_allclose([ex,sx], [e,stress], atol=2e-8, rtol=0)
    np.testing.assert_allclose(gx, 0, atol=2e-7, rtol=0)


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('power', type=int, choices=[1, 2])
    cli.add_argument('output', type=Path)
    args = cli.parse_args()
    run(args.power, args.output)
