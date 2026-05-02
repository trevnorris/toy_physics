#!/usr/bin/env python3
"""Extended reduced branch-family scan with an outgoing-profile curvature coordinate."""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass

import mpmath as mp
import numpy as np
import sympy as sp


def assert_close(label: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label} failed: {actual} vs {expected} (tol={tol})")


@dataclass(frozen=True)
class BranchEval:
    residuals: tuple[float, float, float, float]
    norm: float


w = sp.symbols("w", real=True)
MODES = (0, 2, 4, 6)
mp.mp.dps = 50


def basis_mode(n: int) -> sp.Expr:
    return sp.simplify(
        sp.hermite(n, w) * sp.exp(-w**2 / 2) / (sp.pi ** sp.Rational(1, 4) * sp.sqrt(2**n * sp.factorial(n)))
    )


BASIS = [basis_mode(n) for n in MODES]
BASIS_PRIME = [sp.diff(expr, w) for expr in BASIS]
BASIS_F = [sp.lambdify(w, expr, "mpmath") for expr in BASIS]
BASIS_PRIME_F = [sp.lambdify(w, expr, "mpmath") for expr in BASIS_PRIME]


def quad_inf(func) -> float:
    return float(mp.quad(func, [-mp.inf, mp.inf]))


def symmetric_matrix(size: int, entry_fn) -> np.ndarray:
    out = np.zeros((size, size), dtype=float)
    for i in range(size):
        for j in range(i, size):
            value = entry_fn(i, j)
            out[i, j] = value
            out[j, i] = value
    return out


def source_vector(size: int, entry_fn) -> np.ndarray:
    return np.array([entry_fn(i) for i in range(size)], dtype=float)


def low_series_packet(mass_matrix: np.ndarray, stiff_matrix: np.ndarray, coeffs: np.ndarray) -> tuple[float, float, float]:
    a0 = np.linalg.solve(stiff_matrix, coeffs)
    a2 = -np.linalg.solve(stiff_matrix, mass_matrix @ a0)
    a4 = -np.linalg.solve(stiff_matrix, mass_matrix @ a2)
    return float(coeffs @ a0), float(coeffs @ a2), float(coeffs @ a4)


def evaluate_branch(log_amp: float, log_r_width: float, log_beta_width: float, outgoing_delta: float) -> BranchEval:
    amp = math.exp(log_amp)
    r_width = math.exp(log_r_width)
    beta_width = math.exp(log_beta_width)
    outgoing_weight = 0.5 + outgoing_delta

    R0 = amp * sp.exp(-r_width * w**2 / 2)
    beta = amp * sp.exp(-beta_width * w**2 / 2)
    R0p = sp.diff(R0, w)

    mu_eta = 1 + (1 + w**2 / 4) * R0
    tw_shape = 1 + w**2 / 6
    T_w = 1 + tw_shape * R0
    T_omega = (1 + (1 + w**2 / 8) * R0) / 6
    U_scale = sp.simplify((sp.diff(T_w * R0p, w) - tw_shape * R0p**2 / 2) / R0)
    K_eta = sp.simplify(U_scale - sp.diff(tw_shape * R0p, w))

    mu_f = sp.lambdify(w, mu_eta, "mpmath")
    tw_f = sp.lambdify(w, T_w, "mpmath")
    to_f = sp.lambdify(w, T_omega, "mpmath")
    keta_f = sp.lambdify(w, K_eta, "mpmath")
    beta_f = sp.lambdify(w, beta, "mpmath")
    beta_prime_f = sp.lambdify(w, sp.diff(beta, w), "mpmath")
    phi_B_f = sp.lambdify(w, beta, "mpmath")
    phi_Z_f = sp.lambdify(w, R0 * beta, "mpmath")
    phi_N_f = sp.lambdify(w, (1 + outgoing_weight * w**2) * beta, "mpmath")

    size = len(MODES)
    mass_matrix = symmetric_matrix(
        size,
        lambda i, j: quad_inf(lambda x: mu_f(x) * BASIS_F[i](x) * BASIS_F[j](x)),
    )
    stiff_matrix = symmetric_matrix(
        size,
        lambda i, j: quad_inf(
            lambda x: tw_f(x) * BASIS_PRIME_F[i](x) * BASIS_PRIME_F[j](x)
            + (keta_f(x) + 6.0 * to_f(x)) * BASIS_F[i](x) * BASIS_F[j](x)
        ),
    )

    coeffs_B = source_vector(size, lambda i: quad_inf(lambda x: BASIS_F[i](x) * phi_B_f(x)))
    coeffs_Z = source_vector(size, lambda i: quad_inf(lambda x: BASIS_F[i](x) * phi_Z_f(x)))
    coeffs_N = source_vector(size, lambda i: quad_inf(lambda x: BASIS_F[i](x) * phi_N_f(x)))

    B0, B2, B4 = low_series_packet(mass_matrix, stiff_matrix, coeffs_B)
    Z0, Z2, Z4 = low_series_packet(mass_matrix, stiff_matrix, coeffs_Z)
    N0, N2, N4 = low_series_packet(mass_matrix, stiff_matrix, coeffs_N)

    M_sigma = quad_inf(lambda x: mu_f(x) * beta_f(x) ** 2)
    K_sigma = quad_inf(lambda x: tw_f(x) * beta_prime_f(x) ** 2 + (keta_f(x) + 6.0 * to_f(x)) * beta_f(x) ** 2)

    D0 = K_sigma - B0 - Z0
    D2 = -(M_sigma + B2 + Z2)
    D4 = -(B4 + Z4)
    P0 = N0 / D0
    P2 = (D0 * N2 - 2.0 * D2 * N0) / D0**2
    P4 = (D0**2 * N4 - 2.0 * D0 * (D2 * N2 + D4 * N0) + 3.0 * D2**2 * N0) / D0**3

    R_pole = D0 * (B4 + Z4) - 3.0 * (M_sigma + B2 + Z2) ** 2
    R_norm = P0 - 54.0 / 5.0
    residuals = (R_pole, R_norm, P2, P4)
    return BranchEval(residuals=residuals, norm=float(np.linalg.norm(np.array(residuals, dtype=float))))


def central_jacobian(base: np.ndarray, h: float) -> np.ndarray:
    columns: list[np.ndarray] = []
    for idx in range(base.size):
        delta = np.zeros_like(base)
        delta[idx] = h
        plus = evaluate_branch(*(base + delta))
        minus = evaluate_branch(*(base - delta))
        columns.append((np.array(plus.residuals) - np.array(minus.residuals)) / (2.0 * h))
    return np.column_stack(columns)


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_outgoing_family_scan",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "geometry_status": "reduced_parent_background_plus_outgoing_curvature_family",
        "family_coordinates": [
            "log_R0_amplitude",
            "log_R0_inverse_width",
            "log_beta_inverse_width",
            "delta_outgoing_quadratic_weight",
        ],
        "baseline_coordinates": [0.0, 0.0, 0.0, 0.0],
        "basis_family": [0, 2, 4, 6],
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    baseline = evaluate_branch(0.0, 0.0, 0.0, 0.0)
    expected_residuals = (-13.134593938872376, -10.33719584868593, 0.37009844569768474, 0.8889149882257381)
    for idx, expected in enumerate(expected_residuals):
        assert_close(f"baseline residual {idx}", baseline.residuals[idx], expected, 1e-8)

    h = 0.02
    base_coords = np.zeros(4, dtype=float)
    jac = central_jacobian(base_coords, h)
    svals = np.linalg.svd(jac, compute_uv=False)
    rank = int(np.sum(svals > 1e-6))

    residual_vec = np.array(baseline.residuals, dtype=float)
    delta_ls, _, _, _ = np.linalg.lstsq(jac, -residual_vec, rcond=None)
    residual_best_linear = residual_vec + jac @ delta_ls
    irreducible_linear_norm = float(np.linalg.norm(residual_best_linear))

    trial_scales = (0.01, 0.02, 0.03, 0.05, 0.08, 0.10)
    trial_results = []
    for scale in trial_scales:
        trial_eval = evaluate_branch(*(base_coords + scale * delta_ls))
        trial_results.append((scale, trial_eval))
    best_scale, best_eval = min(trial_results, key=lambda item: item[1].norm)

    assert rank == 4
    assert_close("smallest singular value", float(svals[-1]), 0.108888019657, 1e-6)
    assert_close("linearized residual norm", irreducible_linear_norm, 4.4642745655749085e-13, 1e-10)
    assert best_eval.norm < baseline.norm
    assert abs(best_eval.residuals[2]) < 1e-3
    assert abs(best_eval.residuals[3]) < 1e-3
    assert abs(best_eval.residuals[1]) > 10.0

    print("STEP 21 OUTGOING FAMILY SCAN AUDIT")
    print("Extended the reduced branch-family scan by one upstream outgoing-profile curvature coordinate.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  geometry_status =", branch_metadata["geometry_status"])
    print("Frozen baseline residual packet:")
    print("  (R_pole, R_norm, R_P2, R_P4) =", baseline.residuals)
    print("  residual norm =", baseline.norm)
    print("Extended local family:")
    print("  coordinates = (log R0 amplitude, log R0 inverse-width, log beta inverse-width, delta outgoing quadratic weight)")
    print("  finite-difference step =", h)
    print("  Jacobian rows = (R_pole, R_norm, R_P2, R_P4)")
    print("  Jacobian =")
    for row in jac:
        print("   ", [round(float(x), 12) for x in row.tolist()])
    print("  singular values =", [round(float(x), 12) for x in svals.tolist()])
    print("  numerical rank =", rank)
    print("  least-squares linearized delta =", [round(float(x), 12) for x in delta_ls.tolist()])
    print("  best linearized residual =", tuple(float(x) for x in residual_best_linear))
    print("  irreducible linearized norm =", irreducible_linear_norm)
    print("Actual trial steps along the least-squares direction:")
    for scale, trial_eval in trial_results:
        print("  scale =", scale, "coords =", tuple(float(x) for x in (base_coords + scale * delta_ls)), "norm =", trial_eval.norm, "residuals =", trial_eval.residuals)
    print("  best actual trial scale =", best_scale)
    print("  best actual residual norm =", best_eval.norm)
    print("  best actual residual packet =", best_eval.residuals)
    print("Interpretation:")
    print("  The added outgoing coordinate restores full local rank and can nearly cancel R_P2 and R_P4 on an actual finite step.")
    print("  But R_norm remains O(10), so static normalization is still the dominant obstruction in this enlarged reduced family.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
