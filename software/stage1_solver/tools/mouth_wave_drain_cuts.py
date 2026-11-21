"""Locate the carrier energy transferred by a trapped-wave density perturbation."""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare, trace
from mouth_empty_bend_stiffness import quadrature
from mouth_flared_receiver import EXIT, ROUNDING, emit
from mouth_wave_drain_power import linear_carrier, wave_shape


def section_rule(dom):
    grid = dom.blocks[2][1]
    nodes, weights = leggauss(6)
    half = np.diff(grid)[:, None]/2
    v = (grid[:-1, None]+half*(1+nodes)).ravel()
    return v, (half*weights).ravel()


def run(a):
    dom, psi = prepare(SimpleNamespace(source=Path(f'_scratch/s11c/mouth_flared_cubic_n{a.mesh}')))
    profile, omega = wave_shape(dom, a.step, a.pad)
    hessian, skew, source = linear_carrier(dom, psi, profile)
    _, depths, grad_r, grad_w = quadrature(dom)
    regions = [depths < ROUNDING, (depths >= ROUNDING) & (depths < EXIT), depths >= EXIT]
    f = dom.N @ psi
    n = abs(f)**2
    basis = dom.N[:, dom.free]
    regional_sources = []
    for mask in regions:
        b = basis.T @ (dom.weight * mask * profile * f)
        regional_sources.append(np.r_[b.real, b.imag])
    regional_sources = np.array(regional_sources)
    np.testing.assert_allclose(source, regional_sources.sum(axis=0), atol=1e-14)
    dt = 2*np.pi/(omega*a.steps)
    solve = splu((skew-.5*dt*hessian).tocsc()).solve
    x = np.zeros_like(source)
    cut_work = np.zeros(2)
    applied_work = np.zeros(3)
    last_v = 0.
    full_increment = np.zeros_like(psi)
    full_midpoint = np.zeros_like(psi)
    vnodes, dv = section_rule(dom)
    for tick in range(1, 2*a.steps+1):
        # Scalar seed energy 0.01, trial squared-gap density slope +1.
        v = .01/omega**2 * np.sin(2*np.pi*(tick % a.steps)/a.steps)**2
        dx = solve(dt*(hessian @ x + source*(last_v+v)/2))
        midpoint = x+.5*dx
        applied_work += 2*(regional_sources @ midpoint)*(v-last_v)
        full_increment[dom.free] = dx[:dx.size//2]+1j*dx[dx.size//2:]
        full_midpoint[dom.free] = midpoint[:dx.size//2]+1j*midpoint[dx.size//2:]
        fields = np.column_stack((full_increment, full_midpoint))
        for index, (w, outward) in enumerate(((ROUNDING, -1), (EXIT, 1))):
            vals = trace(dom, 2, np.full_like(vnodes, w), vnodes, fields)
            r, width = vals[0], abs(vals[3])
            delta_field, axial_gradient = vals[-3][:, 0], vals[-1][:, 1]
            cut_work[index] -= outward*(4*np.pi*r**2*width*dv) @ (delta_field.conj()*axial_gradient).real
        x += dx
        last_v = v
        if tick % a.steps:
            continue
        end = full_midpoint+.5*full_increment
        h = dom.N @ end
        energy_density = (.5*(abs(grad_r @ end)**2+abs(grad_w @ end)**2)
                          +(n**4-1)/4*abs(h)**2+2*n**3*(f.conj()*h).real**2)
        stored = np.array([(dom.weight*mask) @ energy_density for mask in regions])
        net_from_flare = (stored-applied_work)[[0, 2]]
        emit(test='carrier_cuts', mesh=a.mesh, steps_per_cycle=a.steps, cycle=tick//a.steps,
             regions=['feeding side of waist', 'waist to outlet', 'receiving chamber'],
             applied_work=applied_work.tolist(), stored=stored.tolist(),
             outward_cumulative_from_local_balances=net_from_flare.tolist(),
             outward_cumulative_from_flux=cut_work.tolist(),
             direct_minus_balance=(cut_work-net_from_flare).tolist())
        np.testing.assert_allclose(stored.sum(), applied_work.sum(), rtol=2e-9)
        np.testing.assert_allclose(stored.sum(), x @ (hessian @ x), rtol=2e-9)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True, type=int)
    parser.add_argument('--step', required=True, type=float)
    parser.add_argument('--pad', required=True, type=float)
    parser.add_argument('--steps', required=True, type=int)
    run(parser.parse_args())
