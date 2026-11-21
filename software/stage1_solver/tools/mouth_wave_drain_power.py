"""Measure reversible density work during two cycles of a real trapped mode.

The scalar gap has a trial slope g with respect to carrier density about the
prepared mouth. The carrier sees the reciprocal potential g Phi^2/2. Measure
the leading order in g from a carrier initially in its stationary state.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import yaml
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.interpolate import RectBivariateSpline
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare
from mouth_curved_phase_stability import weighted_gram
from mouth_flared_receiver import emit
from mouth_trapped_wave_ports import ARCHIVE
from mouth_wave_flare_work import Wave


def wave_shape(dom, step, pad):
    params = yaml.safe_load(ARCHIVE.read_text())['parameters']
    wave = Wave(1., 1., step, pad, params)
    omega, u, op, _, _ = wave.ground(1.)
    residual = np.linalg.norm(op @ u - omega**2*u)
    np.testing.assert_allclose(u @ u, 1.)
    assert residual < 1e-7 and omega < params['mu_out']
    # u is the unit Euclidean eigenvector of the meridional radial operator.
    # The physical scalar u/(sqrt(4*pi)*r*step) has unit four-volume norm.
    r = np.r_[0., wave.r[:, 0], wave.r[-1, 0]+step]
    w = np.r_[wave.w[0, 0]-step, wave.w[0], wave.w[0, -1]+step]
    reduced = np.pad(u.reshape(len(r)-2, len(w)-2), 1)/step
    field = RectBivariateSpline(r, w, reduced)
    roots, _ = leggauss(4)
    xi, eta = [x.ravel() for x in np.meshgrid(roots, roots, indexing='ij')]
    shape = []
    for block, (pg, vg, cells) in dom.blocks.items():
        i, j = np.unravel_index(np.arange(len(cells)), (len(pg)-1, len(vg)-1))
        p = pg[i, None] + .5*np.diff(pg)[i, None]*(1+xi)
        v = vg[j, None] + .5*np.diff(vg)[j, None]*(1+eta)
        rq, wq = (a.ravel() for a in dom.chart(block, p, v)[:2])
        inside = (rq < r[-1]) & (wq > w[0]) & (wq < w[-1])
        values = np.zeros_like(rq)
        values[inside] = field.ev(rq[inside], wq[inside])/(np.sqrt(4*np.pi)*rq[inside])
        shape.append(values)
    profile = np.concatenate(shape)**2
    overlap = dom.weight @ profile
    emit(test='scalar', step=step, pad=pad, omega=omega, exterior_cutoff=params['mu_out'],
         relative_eigen_residual=residual/omega**2, norm_in_carrier=overlap)
    assert 0 < overlap < 1
    return profile, omega


def linear_carrier(dom, psi, profile):
    basis = dom.N[:, dom.free]
    f = dom.N @ psi
    n = abs(f)**2
    weight = dom.weight
    mass = weighted_gram(basis, weight)
    kinetic = dom.stiff[dom.free][:, dom.free]*.5
    scalar = kinetic + weighted_gram(basis, weight*(n**4-1)/4)
    slope = 2*weight*n**3
    cross = weighted_gram(basis, slope*f.real*f.imag)
    operator = sparse.bmat([[scalar+weighted_gram(basis, slope*f.real**2), cross],
                            [cross, scalar+weighted_gram(basis, slope*f.imag**2)]], format='csr')
    skew = sparse.bmat([[None, -mass], [mass, None]], format='csr')
    source = basis.T @ (weight*profile*f)
    return operator, skew, np.r_[source.real, source.imag]


def transfer(dom, field, profile, omega, steps, preparation):
    operator, skew, source = linear_carrier(dom, field, profile)
    dt = 2*np.pi/(omega*steps)
    solve = splu((skew-.5*dt*operator).tocsc()).solve
    # Independent initial scalar energies, reversed slope and original law.
    seeds = np.array([.0025, .01, .01, .01])
    slopes = np.array([1., 1., -1., 0.])
    strength = slopes*seeds/omega**2
    x = np.zeros((len(source), len(seeds)))
    probe_work = np.zeros(len(seeds))
    previous_projection = np.zeros(len(seeds))
    previous_v = np.zeros(len(seeds))
    last_work = np.zeros(len(seeds))
    peak_density = np.zeros(len(seeds))
    f = dom.N @ field
    basis = dom.N[:, dom.free]
    half = x.shape[0]//2
    for tick in range(1, 2*steps+1):
        envelope = np.sin(2*np.pi*(tick % steps)/steps)**2
        v = strength*envelope
        # S xdot=H x+b v. Midpoint treats the conservative carrier exactly in
        # its quadratic energy identity for any timestep and time-varying v.
        dx = solve(dt*(operator @ x + source[:, None]*(v+previous_v)/2))
        x += dx
        projection = 2*(source @ x)
        probe_work += (v-previous_v)*(projection+previous_projection)/2
        previous_projection, previous_v = projection, v
        if tick % (steps//8) == 0:
            density = 2*(f.real[:, None]*(basis @ x[:half]) + f.imag[:, None]*(basis @ x[half:]))
            peak_density = np.maximum(peak_density, np.max(abs(density), axis=0))
            np.testing.assert_allclose((dom.weight*profile) @ density, projection, rtol=1e-8, atol=1e-16)
        if tick % steps == 0:
            stored = np.sum(x*(operator @ x), axis=0)
            emit(test='cycle', preparation=preparation, cycle=tick//steps, steps_per_cycle=steps,
                 starting_scalar_energies=seeds.tolist(), density_gap_slopes=slopes.tolist(),
                 net_into_carrier_this_cycle=(probe_work-last_work).tolist(),
                 work_into_carrier_since_start=probe_work.tolist(),
                 final_carrier_perturbation_energy=stored.tolist(),
                 sampled_peak_abs_delta_n=peak_density.tolist())
            np.testing.assert_allclose(stored, probe_work, rtol=2e-8, atol=1e-17)
            np.testing.assert_allclose(probe_work[1:3], 16*probe_work[0], rtol=2e-8, atol=1e-17)
            assert probe_work[-1] == 0 and np.max(peak_density)<.005
            last_work = probe_work.copy()


def run(a):
    folder = Path(f'_scratch/s11c/mouth_flared_cubic_n{a.mesh}')
    record = json.loads((folder/'stdout').read_text().splitlines()[0])
    for key in ('bulk', 'feed'):
        np.testing.assert_allclose(record[key], .25, rtol=0., atol=1e-8)
    assert record['profile']=='cubic' and record['roof']==10. and record['feed_far']==9.
    dom, field = prepare(SimpleNamespace(source=folder))
    profile, omega = wave_shape(dom, a.step, a.pad)
    transfer(dom, field, profile, omega, a.steps, 'intake_0.25')
    transfer(dom, np.ones_like(field), profile, omega, a.steps, 'rest_control')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True, type=int)
    parser.add_argument('--step', required=True, type=float)
    parser.add_argument('--pad', required=True, type=float)
    parser.add_argument('--steps', required=True, type=int)
    args = parser.parse_args()
    if args.steps < 8 or args.steps % 8:
        parser.error('--steps must be a positive multiple of eight')
    run(args)
