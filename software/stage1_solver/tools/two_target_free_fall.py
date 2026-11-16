"""Compare two target masses with the same spherical intakes and outlet law at release."""

import argparse
from math import pi
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from pair_holding_probe import DENSITY, PROBE, basis
from pair_probe_transfer_check import wolfram_rows


def exterior(distance, degree):
    centers = [0.0, float(distance)] if distance else [0.0]
    count = len(centers)
    mu, _ = leggauss(degree + 1)
    sine = np.sqrt(1 - mu**2)
    blocks, conditions = [], []
    # Column zero drives steady suction. Each subsequent column accelerates one
    # center along +x while both radii and relative suction speeds stay fixed.
    for i, center in enumerate(centers):
        _, x_grad, r_grad = basis(center + mu, sine, centers, degree)
        blocks.append(mu[:, None] * x_grad + sine[:, None] * r_grad)
        target = np.zeros((len(mu), count + 1))
        target[:, 0], target[:, i + 1] = -1, mu
        conditions.append(target)
    coefficients = np.linalg.solve(np.vstack(blocks), np.vstack(conditions))

    mu, weights = leggauss(4 * degree + 9)
    sine = np.sqrt(1 - mu**2)
    load, added_mass = np.zeros(count), np.zeros((count, count))
    boundary_error = 0.0
    for i, center in enumerate(centers):
        potential, x_grad, r_grad = basis(center + mu, sine, centers, degree)
        ux, ur = x_grad @ coefficients[:, 0], r_grad @ coefficients[:, 0]
        normal = (mu[:, None] * x_grad + sine[:, None] * r_grad) @ coefficients
        # For the k=0 outlet the outgoing axial momentum and extra axial traction
        # vanish. Integrate incoming material momentum and Bernoulli pressure.
        load[i] = 2 * pi * DENSITY * weights @ (
            -ux * normal[:, 0] + mu * (ux**2 + ur**2) / 2
        )
        added_mass[i] = -2 * pi * DENSITY * (weights * mu) @ potential @ coefficients[:, 1:]
        normal[:, 0] += 1
        normal[:, i + 1] -= mu
        boundary_error = max(boundary_error, abs(normal).max())
    return load, added_mass, boundary_error


def run_case(distance, degree, baseline, baseline_mobility):
    load, added_mass, normal_error = exterior(distance, degree)
    count = len(load)
    masses = np.ones(count)
    masses[0] = 2
    inertia = np.diag(masses) + added_mass
    free = np.linalg.solve(inertia, load)
    probes = PROBE * np.eye(count)
    minus = np.linalg.solve(inertia, load[:, None] - probes)
    plus = np.linalg.solve(inertia, load[:, None] + probes)
    mobility = (plus - minus) / (2 * PROBE)
    reconstructed = mobility @ load
    effective_mass = 1 / mobility[0, 0]
    hold, free_a = float(baseline[0]), float(baseline[3])
    original_row = np.array(baseline_mobility, dtype=float)

    # Independent prediction from the archived exact-arithmetic, unit-mass
    # target. Adding bare mass changes only the target diagonal of its inertia.
    mass_increment = masses[0] - 1
    divisor = 1 + mass_increment * original_row[0]
    reference_row = original_row / divisor
    reference = np.array([hold, free_a / divisor, *reference_row])
    direct = np.array([-load[0], free[0], *mobility[0]])
    cross_error = abs(direct - reference).max()
    reconstruction_error = abs(reconstructed - free).max()
    reciprocity_error = abs(inertia - inertia.T).max()
    force_sum = load.sum()
    print(f"d={distance or 'isolated'} degree={degree} hold_B={-load[0]:+.13g} "
          f"a_B={free[0]:+.13g} probe_minus={minus[0,0]:+.13g} probe_plus={plus[0,0]:+.13g}")
    print(f"  mu_B={mobility[0].tolist()} calibrated_M_B={effective_mass:.13g} "
          f"calibrated_M_A={1/original_row[0]:.13g}")
    print(f"  reference_hold={hold:+.13g} reference_a_B={reference[1]:+.13g} "
          f"reference_mu_B={reference_row.tolist()} a_A={free_a:+.13g}")
    print(f"  max_reference_error={cross_error:.3g} normal_error={normal_error:.3g} "
          f"reconstruction_error={reconstruction_error:.3g} reciprocity_error={reciprocity_error:.3g} "
          f"net_held_state_load={force_sum:+.3g}", flush=True)
    if distance:
        print(f"  a_B/a_A={free[0]/free_a:.13g} d2_a_A={distance**2*free_a:.13g} "
              f"d2_a_B={distance**2*free[0]:.13g}", flush=True)
    if not np.isfinite(np.r_[direct, reference, normal_error]).all():
        raise ArithmeticError("nonfinite surface or response reading")
    if normal_error > 1e-8 or max(cross_error, reconstruction_error,
                                   reciprocity_error, abs(force_sum)) > 2e-12:
        raise ArithmeticError("surface, force, probe or reference check failed")
    return direct


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mass_one_wolfram_stdout", type=Path)
    rows, mobilities = wolfram_rows(parser.parse_args().mass_one_wolfram_stdout)
    for separation in (0, 10, 20):
        readings = [
            run_case(separation, degree, rows[separation, 0], mobilities[separation])
            for degree in (8, 14)
        ]
        difference = abs(readings[1] - readings[0]).max()
        print(f"  max_refinement_change={difference:.3g}", flush=True)
        if difference > 2e-12:
            raise ArithmeticError("surface-force or probe refinement failed")
