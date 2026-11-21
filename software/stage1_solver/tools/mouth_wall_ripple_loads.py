"""Spatially projected field forces for four weak dipole motions along the neck.

The columns are prescribed displacements. The rows are the conjugate forces
needed to counter the carrier and cycle-averaged trapped-wave loads.
"""
import argparse

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, splu
from scipy.special import expit

from mouth_carrier_bend_budget import settings
from mouth_open_reservoirs import ReservoirDomain
from mouth_smooth_dipole_relaxation import AxialModes
from mouth_wave_outgoing_sidebands import gate


FREQUENCIES = (1.6, 2.)


def shapes(u):
    """Four unit-peak smooth displacement windows of width one quarter arm."""
    phase = 4 * np.asarray(u)[:, None] - np.arange(4)
    return np.where((phase > 0) & (phase < 1), np.sin(np.pi * phase)**2, 0.)


def emit_matrix(field, resolution, omega, matrix):
    for row, entries in enumerate(matrix, 1):
        re = ",".join(f"{x.real:.11g}" for x in entries)
        im = ",".join(f"{x.imag:.11g}" for x in entries)
        print(f"external field={field} resolution={resolution} omega={omega:g} "
              f"read_quarter={row} re={re} im={im}", flush=True)


def carrier_sources(domain, centers, edges):
    n = domain.sleeve_cells
    side, = np.where((domain.r_index == n-1) & (domain.w_index < n))
    inner, = np.where((domain.r_index == n-2) & (domain.w_index < n))
    assert np.array_equal(domain.w_index[side], np.arange(n))
    assert np.array_equal(domain.w_index[inner], np.arange(n))
    trace = ((9*domain.psi[side] - domain.psi[inner])/8)[:, None]
    area = 4*np.pi*domain.radius**2/3
    dh = 1.25*domain.eos_k*(domain.bath_density**4 - abs(trace)**8)
    spatial = 2 * area * domain.dw * centers * dh * trace
    edge_flux = 2 * domain.kinetic * area * edges * np.diff(trace, axis=0)/domain.dw
    spatial[:-1] += edge_flux
    spatial[1:] -= edge_flux
    temporal = -1j * domain.hbar * area * domain.dw * centers * trace
    def trace_adjoint(values):
        y = np.zeros((len(domain.psi), values.shape[1]), dtype=complex)
        y[side], y[inner] = 9*values/8, -values/8
        return np.vstack((y.real, y.imag))
    return trace_adjoint(spatial), trace_adjoint(temporal)


def carrier(cells):
    point, _, _, flow, archive = settings()
    params = archive["parameters"]
    radius = archive["phase_A"]["q_star"][0]
    domain = ReservoirDomain(radius, float(point["L"])/2, radius, cells, 2, 4, 3, params)
    domain.set_current(float(flow["reference_inward_current"])/params["m"])
    centers = (np.arange(cells) + .5)/cells
    edges = np.arange(1, cells)/cells
    # Check generalized sources against the established whole-bend operator.
    def full(u):
        return np.sin(np.pi*u[:, None])**2
    first, second = carrier_sources(domain, full(centers), full(edges))
    assert np.allclose(first[:, 0], domain.spatial, atol=1e-14)
    assert np.allclose(second[:, 0], domain.temporal, atol=1e-14)
    spatial, temporal = carrier_sources(domain, shapes(centers), shapes(edges))
    volume, _, _, _, angular, lap = domain.ac
    weight = volume/3
    p = domain.psi
    hc = 1.25*domain.eos_k
    h = lap/3 + sparse.diags(domain.kinetic*angular + hc*(abs(p)**8-domain.bath_density**4)*weight)
    g = 8 * hc * abs(p)**6 * weight
    cross = g*p.real*p.imag
    for omega in FREQUENCIES:
        skew = 1j*omega*domain.hbar*weight
        operator = sparse.bmat([[h+sparse.diags(g*p.real**2), sparse.diags(cross-skew)],
                                [sparse.diags(cross+skew), h+sparse.diags(g*p.imag**2)]], format="csc")
        source = spatial + 1j*omega*temporal
        response = splu(operator).solve(source/2)
        reaction = source.conj().T @ response
        emit_matrix("fluid", cells, omega, -reaction)
        # Include a complex superposition to check directional cross loads.
        for drive in (*np.eye(4), np.array([1., .25j, -.7, -.1j])):
            power = omega/2 * np.vdot(drive, reaction@drive).imag
            for cut in (1.25, 1.875):
                flux = domain.cuts(omega, response@drive, round(cut*cells), round((cut-.5)*cells))
                assert min(flux) > 0 and abs(flux.sum()-power) < 5e-7*power


def supporting_wave(step, extent, absorber):
    point, _, _, _, archive = settings()
    params = archive["parameters"]
    radius = archive["phase_A"]["q_star"][0]
    length = float(point["L"])
    mesh = AxialModes(radius, params, step, round(7.2/step)*step)
    potential, _, n = gate(mesh, length, step, params)
    h0 = mesh.kinetic + sparse.diags(potential.ravel())
    val, vec = eigsh(h0, k=1, sigma=0., which="LM", v0=mesh.seed, tol=1e-11)
    ws = np.sqrt(val[0])
    support = vec[:, 0].reshape(n-1, n)
    scale = params["I_wave"]/(2*ws)
    depth = shapes(2*mesh.w/length)
    axial = expit((length-2*np.sqrt(mesh.w**2+params["eps_w"]**2))/params["delta_parallel"])
    grad = mesh.depth*(2/params["delta_r"])*mesh.gr*(1-mesh.gr)
    r = step*np.arange(1, n)
    curve = (2/params["delta_r"])*(2*mesh.gr-1)*grad + 2*grad/r
    direct = scale * depth.T @ ((curve@(support**2) * axial/3)[:, None]*depth)
    source = (-support * np.outer(grad, axial)/np.sqrt(3))[:, :, None]*depth[None, :, :]
    closed = h0 + sparse.diags(2*mesh.centrifugal)
    outer = AxialModes(radius, params, step, extent)
    pout, _, m = gate(outer, length, step, params)
    def damping(x):
        return absorber*(np.maximum(x-5, 0)/(extent-5))**3
    cap = damping(step*np.arange(1, m))[:, None]+damping(outer.w)[None, :]
    continuum = outer.kinetic + sparse.diags(pout.ravel()+2*outer.centrifugal-1j*cap.ravel())
    extended = np.zeros((m-1, m, 4))
    extended[:n-1, :n] = source
    upper_sources = extended.reshape(-1, 4)
    lower_sources = source.reshape(-1, 4)
    for omega in FREQUENCIES:
        assert abs(ws-omega) < params["mu_out"] < ws+omega
        low = splu((closed - (ws-omega)**2*sparse.eye(closed.shape[0])).tocsc()).solve(lower_sources)
        high = splu((continuum - (ws+omega)**2*sparse.eye(continuum.shape[0])).tocsc()).solve(upper_sources)
        upper = upper_sources.T @ high
        absorbed = high.conj().T @ (cap.ravel()[:, None]*high)
        assert np.allclose((upper-upper.conj().T)/(2j), absorbed, atol=2e-11, rtol=2e-8)
        load = direct - scale*(lower_sources.T@low + upper)
        emit_matrix("wave", step, omega, load)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="field", required=True)
    fluid = sub.add_parser("fluid")
    fluid.add_argument("--cells", type=int, required=True)
    wave = sub.add_parser("wave")
    wave.add_argument("--step", type=float, required=True)
    wave.add_argument("--extent", type=float, required=True)
    wave.add_argument("--absorber", type=float, required=True)
    args = parser.parse_args()
    if args.field == "fluid":
        carrier(args.cells)
    else:
        supporting_wave(args.step, args.extent, args.absorber)
