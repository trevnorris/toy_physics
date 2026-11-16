"""Predict the reservoir shift of the one-mouth incident flow at leading order in intake."""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legvander
from scipy.special import roots_legendre

from incident_flow_falloff import checked_grid, check_current
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet
from mouth_feed_weak_check import dx4


ARCHIVE = Path("docs/_measurements/incident_flow_falloff.json")
GRIDS = ((48, 4, 16), (96, 4, 16), (96, 6, 24))
DISTANCES = (4, 8, 16)
RATE = 0.25


def preparation(slab, baseline):
    path, digest = baseline
    with verified_packet(path, {path: digest}) as p:
        checked_grid(slab, p)
        if int(p["d"]) != 0:
            raise ValueError("Expected the centered isolated preparation")
        target = check_current(slab, 0, p["currents"])
        density = abs(p["layers" if "layers" in p else "field"][-1])**2
    # Linearize j = J (w/n_cal)|psi|^2 / integral[(w/n_cal)|psi|^2]
    # at uniform psi, keeping the actual normal-speed calibration fixed.
    w = np.maximum(1-slab.x[:, None]**2-slab.r**2, 0)**3
    raw = w/density
    flux = RATE*raw/np.sum(slab.dv*raw)
    check_current(slab, 0, flux[None])
    weights = target*slab.dv
    return flux, weights/weights.sum()


def separated_speed(slab, flux, weights, d, radial_shift, measure, finite):
    """Top-layer velocity from separated radial, normal and axial Green functions."""
    h = slab.dx
    xs = abs(slab.x) < 1
    rs = slab.r < 1
    offset = slab.x[xs, None]/h
    j = np.arange(np.count_nonzero(rs))
    # For radius r_j=(j+1/2)h, the finite-volume radial recurrence
    # has eigenfunction P_j(mu) and eigenvalue (1-mu)/h^2 for -lap/2.
    # On the infinite grid its normalized spectral measure is dmu/2.
    modes = legvander(1-radial_shift, j[-1])
    source = (flux[xs][:, rs]*(2*j+1)) @ modes.T
    target = weights[xs][:, rs] @ modes.T
    cells = len(slab.x)
    separation = round(d/h)
    if not np.isclose(separation*h, d) or separation+offset.min()-offset.max() <= 2:
        raise ValueError("Derivative stencil overlaps source or translation is off-grid")
    normal_cells = len(slab.sqrt_z)-1
    speed = 0.
    for k in range(normal_cells+1):
        ez = 2*(np.sin(np.pi*k/(2*normal_cells))/slab.dz)**2
        alpha = 2*np.arcsinh(np.sqrt((radial_shift+h*h*ez)/2))
        # The infinite axial inverse is h^2 exp(-alpha |i-m|)/sinh(alpha).
        # Acting with minus the same four-point derivative as the measurement
        # leaves h(4-cosh(alpha))/3 times the exponential at disjoint mouths.
        src_decay = np.exp(-alpha*(separation/2-offset))
        tgt_decay = np.exp(-alpha*(separation/2+offset))
        kernel = h*(4-np.cosh(alpha))/3
        if finite:
            # Odd reflection at the half-cell Dirichlet faces factors the
            # 1D Green function into sinh at the source and sinh at the target.
            # Differentiating the target sinh gives its cosh (the plus below).
            src_decay *= -np.expm1(-2*alpha*(cells/2-separation+offset))
            tgt_decay *= 1+np.exp(-2*alpha*(cells/2-offset))
            kernel /= -np.expm1(-2*alpha*cells)
        # The top-to-top DCT-I projector divided by dz gives half weight
        # to the constant and Nyquist modes, full weight to the others.
        depth_weight = (0.5 if k in (0, normal_cells) else 1.)/(normal_cells*slab.dz)
        speed += depth_weight*np.dot(measure*kernel,
                                      (src_decay*source).sum(axis=0)*(tgt_decay*target).sum(axis=0))
    return speed


def open_spectrum(order):
    nodes, weights = roots_legendre(order)
    t = (1+nodes)/2
    # mu = 1-2 t^2 removes the square-root endpoint of the constant depth mode.
    return 2*t*t, t*weights


def predictions(grid, baseline):
    extent, cells, depth = grid
    slab = LateralSlab(extent, cells, depth)
    flux, weights = preparation(slab, baseline)
    # The axial and normal eigenvalues can be subtracted from the tensor sum
    # at its lowest mode. The first radial cell has P_0=1.
    ex0 = 2*(np.sin(np.pi/(2*len(slab.x)))/slab.dx)**2
    radial_shift = slab.dx**2*(slab.eigenvalues[0, 0]-ex0)
    measure = slab.modes[0]**2
    finite, opened = [], []
    for d in DISTANCES:
        source = np.roll(flux, -round(d/slab.dx), axis=0)
        phase = slab.elliptic(np.zeros(slab.eigenvalues.shape), source, 0.)
        direct = -np.sum(weights*dx4(phase[-1], slab.dx))
        del phase
        bounded = separated_speed(slab, flux, weights, d, radial_shift, measure, True)
        free = [separated_speed(slab, flux, weights, d, *open_spectrum(n), False)
                for n in (64, 128)]
        relative_check = abs(bounded/direct-1)
        quadrature_change = abs(free[0]/free[1]-1)
        print(f"L={extent} c={cells} d={d} weak_finite_speed={direct:.15g} "
              f"weak_open_speed={free[1]:.15g} finite_solver_relative_check={relative_check:.4g} "
              f"open_quadrature_relative_change={quadrature_change:.4g}", flush=True)
        if (not np.isfinite([*free, direct, bounded]).all() or direct == 0 or
                relative_check > 2e-10 or quadrature_change > 2e-10):
            raise ArithmeticError("Finite solver or open radial quadrature check failed")
        finite.append(direct)
        opened.append(free[1])
    return np.array([finite, opened])


def main():
    archive = json.loads(ARCHIVE.read_text())
    rows = {r["name"]: r for r in archive["runs"]}
    predicted = {}
    for grid in GRIDS:
        extent, cells, _ = grid
        runs = [rows[f"L{extent}C{cells}_d{d}"] for d in DISTANCES]
        controls = []
        for run in runs:
            if run["guard_exit"] or run["worker_exit"] or not run["limits_verified"]:
                raise ValueError("Unverified input")
            baseline, = [(Path(p), h) for p, h in run["input_sha256"].items()
                         if Path(p).name in ("d0.npz", f"L{extent}C{cells}_baseline.npz")]
            controls.append(baseline)
        if len(set(controls)) != 1:
            raise ValueError("Calibrations differ across separations")
        predicted[grid] = predictions(grid, controls[0])

    # Compare only after computing the finite and open predictions on every grid.
    ratios = {}
    for grid, speeds in predicted.items():
        extent, cells, _ = grid
        measured = []
        for d in DISTANCES:
            row = json.loads(rows[f"L{extent}C{cells}_d{d}"]["stdout"][-1])
            if tuple(row["grid"]) != grid or row["d"] != d:
                raise ValueError("Mismatched nonlinear readout")
            measured.append(row["mean_speed"])
        all_speeds = np.vstack([speeds, measured])
        ratio = all_speeds[:, 1:]/all_speeds[:, :-1]
        ratios[grid] = ratio
        for label, (finite, opened, observed) in zip(("4_to_8", "8_to_16"), ratio.T, strict=True):
            print(f"L={extent} c={cells} leg={label} weak_finite={finite:.12g} "
                  f"weak_open={opened:.12g} nonlinear_measured={observed:.12g} "
                  f"open_minus_finite_pp={100*(opened-finite):.9g} "
                  f"measured_minus_weak_finite_pp={100*(observed-finite):.9g} "
                  f"weak_open_shortfall_pp={100*(0.25-opened):.9g}")
    for label, (finite, opened, observed) in zip(("4_to_8", "8_to_16"),
            (100*(ratios[GRIDS[1]]-ratios[GRIDS[0]])).T, strict=True):
        print(f"move_48_to_96 leg={label} predicted_shift_pp={finite:.12g} "
              f"measured_shift_pp={observed:.12g} difference_pp={observed-finite:.9g} "
              f"open_calibration_shift_pp={opened:.9g}")


if __name__ == "__main__":
    main()
