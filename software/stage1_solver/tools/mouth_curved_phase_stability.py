"""Check a saved Q2 carrier against a half-plane bound and the fixed-bath energy.

The Bernstein certificate is exact for the stored binary nodal values and bias.
The spectral checks use the solver's geometry, basis and four-point quadrature.
No carrier equation or front equation is solved again.
"""

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from mouth_curved_phase import Domain


def rotation_bounds(bias):
    """Enclose cos(bias/2), sin(bias/2) by dyadics, using alternating series."""
    half = Fraction.from_float(bias) / 2
    x = abs(half)
    if x > 1:
        raise ValueError("The alternating-series enclosure requires |bias/2| <= 1")
    scale = 1 << 70
    bounds = []
    for parity in (0, 1):
        terms = [(-1)**k * x ** (2*k + parity) / math.factorial(2*k + parity)
                 for k in range(10)]
        lower, upper = sum(terms), sum(terms[:-1])
        if parity and half < 0:
            lower, upper = -upper, -lower
        bounds.append((math.floor(lower * scale), math.ceil(upper * scale)))
    return bounds, scale


def half_plane_bound(domain, field, bias):
    # Map [-1,0,1] nodal values to the three quadratic Bernstein coefficients,
    # multiplied by two. In two dimensions the tensor product gives a factor four.
    transform = np.array([[2, 0, 0], [-1, 4, -1], [0, 0, 2]], dtype=object)
    pairs = [[v.as_integer_ratio() for v in values] for values in (field.real, field.imag)]
    scale = max(den for component in pairs for _, den in component)
    nodes = [np.array([num * (scale // den) for num, den in component], dtype=object)
             for component in pairs]
    dofs = np.concatenate([block[2] for block in domain.blocks.values()]).reshape(-1, 3, 3)
    trig, trig_scale = rotation_bounds(bias)
    numerator = np.zeros(dofs.shape, dtype=object)
    for component, (lower, upper) in zip(nodes, trig, strict=True):
        coeff = transform @ component[dofs] @ transform.T
        numerator += np.where(coeff >= 0, lower, upper) * coeff
    bound = Fraction(int(numerator.min()), 4 * scale * trig_scale)
    # The float field is rounded towards minus infinity for convenient quoting.
    rounded = float(bound)
    if Fraction.from_float(rounded) > bound:
        rounded = math.nextafter(rounded, -math.inf)
    return dict(q2_cells=len(dofs), coefficients=numerator.size,
                half_plane_lower=rounded, exact_rational_lower=str(bound))


def weighted_gram(basis, weights):
    return basis.T @ basis.multiply(weights[:, None])


def lowest_mode(operator, mass, floor):
    # The kinetic and radial-response terms are positive semidefinite. A shift
    # strictly below min V reaches the algebraically smallest mode from below.
    initial = np.linspace(0.75, 1.25, operator.shape[0])
    eigenvalues, vectors = eigsh(operator, M=mass, k=1, sigma=floor - 0.1,
                                which="LM", v0=initial, tol=1e-11)
    eigenvalue, vector = eigenvalues[0], vectors[:, 0]
    mv = mass @ vector
    residual = np.linalg.norm(operator @ vector - eigenvalue * mv)
    return dict(eigenvalue=float(eigenvalue),
                relative_residual=float(residual / (abs(eigenvalue) * np.linalg.norm(mv))))


def check(directory, full_hessian):
    state, = [row for row in (json.loads(s) for s in (directory / "stdout").read_text().splitlines())
              if row["test"] == "stationary"]
    domain = Domain(state["n"], state["far"])
    with np.load(directory / "field.npz") as archive:
        np.testing.assert_array_equal(archive["points"], domain.points)
        field = archive["psi"]
    bias = state["bias"]
    np.testing.assert_allclose(field[domain.feed], 1, rtol=0, atol=1e-14)
    np.testing.assert_allclose(field[domain.bulk], np.exp(1j * bias), rtol=0, atol=1e-14)
    certificate = half_plane_bound(domain, field, bias)

    free = domain.free
    basis = domain.N[:, free]
    sampled = domain.N @ field
    n = np.abs(sampled)**2
    potential = (n**4 - 1) / 4
    kinetic = domain.stiff[free][:, free] * 0.5
    mass = weighted_gram(basis, domain.weight)
    scalar = kinetic + weighted_gram(basis, domain.weight * potential)
    residual = 0.5 * (domain.stiff @ (field - 1))[free]
    residual += basis.T @ (domain.weight * potential * sampled)
    result = dict(n=state["n"], far=state["far"], rate=state["rate"], bias=bias,
                  free_scalar_dofs=int(free.sum()), **certificate,
                  projected_residual=float(np.max(abs(residual))),
                  sampled_potential_floor=float(potential.min()),
                  scalar_spectrum=lowest_mode(scalar, mass, potential.min()))
    if full_hessian:
        w = 2 * domain.weight * n**3
        mixed = weighted_gram(basis, w * sampled.real * sampled.imag)
        blocks = [[scalar + weighted_gram(basis, w * sampled.real**2), mixed],
                  [mixed, scalar + weighted_gram(basis, w * sampled.imag**2)]]
        result["half_hessian_spectrum"] = lowest_mode(
            sparse.bmat(blocks, format="csr"), sparse.block_diag([mass, mass]), potential.min())
    print(json.dumps(result), flush=True)
    for key in ("scalar_spectrum", "half_hessian_spectrum"):
        if key in result:
            mode = result[key]
            if mode["eigenvalue"] <= 0 or mode["relative_residual"] > 1e-7:
                raise ArithmeticError(f"{key} failed the fixed-bath spectral check")
    if certificate["half_plane_lower"] <= 0 or result["projected_residual"] > 3e-11:
        raise ArithmeticError("The stored branch failed its positivity or stationarity check")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", type=Path, nargs="+")
    parser.add_argument("--full-hessian", action="store_true")
    args = parser.parse_args()
    for run in args.directories:
        check(run, args.full_hessian)
