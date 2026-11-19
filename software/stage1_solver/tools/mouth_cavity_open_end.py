"""Rayleigh certificates for a scalar finite cylinder with an open receiver.

The well is B^3(R) x (0,L) in R^3 x R_+, with a Dirichlet plane at z=0.
A compact change of squared gap leaves the continuum starting at mu_out^2.
Use a product f(r) Y_lm g(z) with continuous piecewise linear f and g,
extended by zero beyond their mesh. Each mesh continues past the well.
The angular harmonics have unit norm. Normalizations cancel in the quotient.
SciPy selects trial vectors. All reported comparisons to threshold use exact
rational evaluation of integer-rounded vectors and archived decimal inputs.
"""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path

import numpy as np
import yaml
from scipy.linalg import eigh


ARCHIVE = Path(__file__).resolve().parents[1] / "reports/pathA_26_derrick_results.yaml"
CELLS_INSIDE = 16
CELLS_TOTAL = 4 * CELLS_INSIDE


def forms(ell):
    """Dimensionless mass, inside mass and kinetic matrices times 60.

    Grid coordinate is distance in units of one cell. For the radial
    factor the mass and derivative terms carry x^2 and the angular
    term carries ell*(ell+1). For the axial factor ell is None.
    """
    shape = (CELLS_TOTAL + 1,) * 2
    mass = np.zeros(shape, dtype=int)
    inside = np.zeros_like(mass)
    kinetic = np.zeros_like(mass)
    for i in range(CELLS_TOTAL):
        pair = np.s_[i:i+2, i:i+2]
        if ell is None:
            a, b, d = 20, 10, 20
            grad, angular = 60, 0
        else:
            a = 20*i*i + 10*i + 2
            b = 10*i*i + 10*i + 3
            d = 20*i*i + 30*i + 12
            grad, angular = 60*i*i + 60*i + 20, ell*(ell+1)
        local = np.array([[a, b], [b, d]])
        mass[pair] += local
        if i < CELLS_INSIDE:
            inside[pair] += local
        kinetic[pair] += [[grad + 20*angular, -grad + 10*angular],
                          [-grad + 10*angular, grad + 20*angular]]
    # Regular monopole at origin. Dipole and axial amplitudes vanish there.
    active = slice(0 if ell == 0 else 1, CELLS_TOTAL)
    return tuple(matrix[active, active] for matrix in (mass, inside, kinetic))


def choose_factor(matrices, cell, weighted_depth):
    mass, inside, kinetic = matrices
    hamiltonian = kinetic / float(cell)**2 - float(weighted_depth) * inside
    _, columns = eigh(hamiltonian, mass, subset_by_index=(0, 0))
    v = columns[:, 0]
    return v, float((v @ inside @ v) / (v @ mass @ v))


def exact_factor(matrices, cell, vector):
    v = np.rint(vector / np.max(np.abs(vector)) * 10**6).astype(int).astype(object)
    norm, inner, gradient = (int(v @ m.astype(object) @ v) for m in matrices)
    return Q(gradient, norm) / cell**2, Q(inner, norm)


def product(ell, radius, length, speed, gap_in, gap_out):
    depth = gap_out**2 - gap_in**2
    r_mats, z_mats = forms(ell), forms(None)
    dr, dz = radius / CELLS_INSIDE, length / CELLS_INSIDE
    probability_z = 1
    for _ in range(16):
        r_vec, probability_r = choose_factor(r_mats, dr, depth * probability_z / speed**2)
        z_vec, probability_z = choose_factor(z_mats, dz, depth * probability_r / speed**2)
    tr, pr = exact_factor(r_mats, dr, r_vec)
    tz, pz = exact_factor(z_mats, dz, z_vec)
    return gap_out**2 + speed**2 * (tr + tz) - depth*pr*pz, pr, pz


def upper_decimal(value, places=8):
    """Round toward +infinity, preserving the variational bound in print."""
    scale = 10**places
    units = -((-value.numerator * scale) // value.denominator)
    return f"{units // scale}.{units % scale:0{places}d}"


def dipole_exclusion(length, speed, depth):
    """Return a conservative radius below which the dipole form is nonnegative
    relative to the continuum. The separate one-dimensional inequalities are
    T_r >= (pi/R)^2 P_r and T_z >= (pi/(2L))^2 P_z. Each inside mass is
    at least the mass in the product well, for an arbitrary joint state.

    The axial inequality follows from the first Dirichlet-Neumann eigenvalue
    on (0,L). For the radial inequality, the least-cost exterior continuation
    is f(R)*(R/r)^2. Its slope gives the Robin condition R f'(R)=-2 f(R).
    The first interior eigenfunction j_1(kr) then obeys j_0(kR)=0, so kR=pi.
    """
    pi_lower = Q(333, 106)
    remaining = depth - (speed * pi_lower / (2*length))**2
    if remaining <= 0:
        return "all"
    radius_squared = (speed * pi_lower)**2 / remaining
    scale = 10**6
    units = isqrt(radius_squared.numerator * scale**2 // radius_squared.denominator)
    return f"{units // scale}.{units % scale:06d}"


def main():
    source = yaml.load(ARCHIVE.read_text(), Loader=yaml.BaseLoader)
    r0, length = map(Q, source["phase_A"]["q_star"])
    speed, gap_in, gap_out = (Q(source["parameters"][k]) for k in ("c_w", "mu_in", "mu_out"))
    print(f"R0={float(r0):.12f} L0={float(length):.12f} continuum_omega2={gap_out**2}")
    excluded = dipole_exclusion(length, speed, gap_out**2 - gap_in**2)
    print(f"dipole_below_continuum_excluded_for_R_at_most={excluded} (sufficient)")
    for factor in map(Q, ("0.9", "1", "1.125")):
        for ell in (0, 1):
            upper, pr, pz = product(ell, factor*r0, length, speed, gap_in, gap_out)
            print(f"R/R0={float(factor):g} ell={ell} omega2_upper={upper_decimal(upper)} "
                  f"Pr={float(pr):.8f} Pz={float(pz):.8f} "
                  f"bound_state_certified={upper < gap_out**2}")


if __name__ == "__main__":
    main()
