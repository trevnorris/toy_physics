#!/usr/bin/env python3
"""Reduced branch-family sensitivity scan around the frozen step-19 packet."""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass

import mpmath as mp
import numpy as np
import sympy as sp

from step_19_parent_throat_action_actual_branch_export_sympy import export_step19_four_mode_residuals


def assert_close(label: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label} failed: {actual} vs {expected} (tol={tol})")


def assert_small(label: str, value: float, tol: float) -> None:
    if abs(value) > tol:
        raise AssertionError(f"{label} failed: |value|={abs(value)} > {tol}")


def assert_greater(label: str, actual: float, floor: float) -> None:
    if actual <= floor:
        raise AssertionError(f"{label} failed: {actual} <= {floor}")


@dataclass(frozen=True)
class BranchEval:
    M_sigma: float
    K_sigma: float
    B: tuple[float, float, float]
    Z: tuple[float, float, float]
    N: tuple[float, float, float]
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


def build_numeric_functions(
    log_amp: float,
    log_r_width: float,
    log_beta_width: float,
    *,
    mu_eta_profile_sign: float = 1.0,
) -> dict[str, object]:
    amp = math.exp(log_amp)
    r_width = math.exp(log_r_width)
    beta_width = math.exp(log_beta_width)

    R0 = amp * sp.exp(-r_width * w**2 / 2)
    beta = amp * sp.exp(-beta_width * w**2 / 2)
    R0p = sp.diff(R0, w)
    tw_shape = 1 + w**2 / 6
    mu_eta = 1 + mu_eta_profile_sign * (1 + w**2 / 4) * R0
    T_w = 1 + tw_shape * R0
    T_omega = (1 + (1 + w**2 / 8) * R0) / 6
    U_scale = sp.simplify((sp.diff(T_w * R0p, w) - tw_shape * R0p**2 / 2) / R0)
    K_eta = sp.simplify(U_scale - sp.diff(tw_shape * R0p, w))

    return {
        "amp": amp,
        "r_width": r_width,
        "beta_width": beta_width,
        "mu_eta": sp.lambdify(w, mu_eta, "mpmath"),
        "T_w": sp.lambdify(w, T_w, "mpmath"),
        "T_omega": sp.lambdify(w, T_omega, "mpmath"),
        "K_eta": sp.lambdify(w, K_eta, "mpmath"),
        "beta": sp.lambdify(w, beta, "mpmath"),
        "beta_prime": sp.lambdify(w, sp.diff(beta, w), "mpmath"),
        "phi_B": sp.lambdify(w, beta, "mpmath"),
        "phi_Z": sp.lambdify(w, R0 * beta, "mpmath"),
        "phi_N": sp.lambdify(w, (1 + w**2 / 2) * beta, "mpmath"),
    }


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


def low_series_packet(
    mass_matrix: np.ndarray,
    stiff_matrix: np.ndarray,
    coeffs: np.ndarray,
) -> tuple[float, float, float]:
    a0 = np.linalg.solve(stiff_matrix, coeffs)
    a2 = -np.linalg.solve(stiff_matrix, mass_matrix @ a0)
    a4 = -np.linalg.solve(stiff_matrix, mass_matrix @ a2)
    return float(coeffs @ a0), float(coeffs @ a2), float(coeffs @ a4)


def evaluate_branch(
    log_amp: float,
    log_r_width: float,
    log_beta_width: float,
    *,
    mu_eta_profile_sign: float = 1.0,
) -> BranchEval:
    funcs = build_numeric_functions(
        log_amp,
        log_r_width,
        log_beta_width,
        mu_eta_profile_sign=mu_eta_profile_sign,
    )
    mu_f = funcs["mu_eta"]
    tw_f = funcs["T_w"]
    to_f = funcs["T_omega"]
    keta_f = funcs["K_eta"]
    beta_f = funcs["beta"]
    beta_prime_f = funcs["beta_prime"]
    phi_B_f = funcs["phi_B"]
    phi_Z_f = funcs["phi_Z"]
    phi_N_f = funcs["phi_N"]

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
    norm = float(np.linalg.norm(np.array(residuals, dtype=float)))

    return BranchEval(
        M_sigma=M_sigma,
        K_sigma=K_sigma,
        B=(B0, B2, B4),
        Z=(Z0, Z2, Z4),
        N=(N0, N2, N4),
        residuals=residuals,
        norm=norm,
    )


def central_jacobian(base: np.ndarray, h: float) -> tuple[np.ndarray, list[BranchEval]]:
    columns: list[np.ndarray] = []
    probes: list[BranchEval] = []
    for idx in range(base.size):
        delta = np.zeros_like(base)
        delta[idx] = h
        plus = evaluate_branch(*(base + delta))
        minus = evaluate_branch(*(base - delta))
        probes.extend([plus, minus])
        columns.append((np.array(plus.residuals) - np.array(minus.residuals)) / (2.0 * h))
    return np.column_stack(columns), probes


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_reduced_family_scan",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "geometry_status": "reduced_parent_background_profile_family",
        "family_coordinates": ["log_R0_amplitude", "log_R0_inverse_width", "log_beta_inverse_width"],
        "baseline_coordinates": [0.0, 0.0, 0.0],
        "basis_family": [0, 2, 4, 6],
    }

    baseline = evaluate_branch(0.0, 0.0, 0.0)
    expected_B = (0.7816273402373896, -0.6449331985854891, 0.5321531025367578)
    expected_Z = (0.5083076399936368, -0.4143227047792866, 0.3408544970748318)
    expected_N = (1.311690146078847, -1.062371646707426, 0.8725297055800806)
    expected_residuals = export_step19_four_mode_residuals()

    for idx, expected in enumerate(expected_B):
        assert_close(f"baseline B{2*idx}", baseline.B[idx], expected, 1e-8)
    for idx, expected in enumerate(expected_Z):
        assert_close(f"baseline Z{2*idx}", baseline.Z[idx], expected, 1e-8)
    for idx, expected in enumerate(expected_N):
        assert_close(f"baseline N{2*idx}", baseline.N[idx], expected, 1e-8)
    for idx, expected in enumerate(expected_residuals):
        assert_close(f"baseline residual {idx}", baseline.residuals[idx], expected, 1e-8)

    h = 0.02
    base_coords = np.zeros(3, dtype=float)
    jac, _ = central_jacobian(base_coords, h)
    svals = np.linalg.svd(jac, compute_uv=False)
    rank = int(np.sum(svals > 1e-6))

    residual_vec = np.array(baseline.residuals, dtype=float)
    delta_ls, _, _, _ = np.linalg.lstsq(jac, -residual_vec, rcond=None)
    residual_best_linear = residual_vec + jac @ delta_ls
    irreducible_linear_norm = float(np.linalg.norm(residual_best_linear))

    trial_scales = (0.05, 0.1, 0.2)
    trial_results = []
    for scale in trial_scales:
        trial_eval = evaluate_branch(*(base_coords + scale * delta_ls))
        trial_results.append((scale, trial_eval))
    best_scale, best_eval = min(trial_results, key=lambda item: item[1].norm)
    mu_eta_sign_flip = evaluate_branch(0.0, 0.0, 0.0, mu_eta_profile_sign=-1.0)
    mu_eta_sign_flip_delta_norm = float(
        np.linalg.norm(np.array(mu_eta_sign_flip.residuals, dtype=float) - np.array(baseline.residuals, dtype=float))
    )

    if rank != 3:
        raise AssertionError(f"expected full 3-column numerical rank, got {rank}")
    assert_close("smallest singular value", float(svals[-1]), 0.14131229158, 1e-9)
    assert_small("baseline-vs-step19 4-mode norm drift", baseline.norm - float(np.linalg.norm(expected_residuals)), 1e-8)
    assert_greater("mu_eta sign-flip residual-packet delta", mu_eta_sign_flip_delta_norm, 1.0)

    branch_freeze_payload = {
        "metadata": branch_metadata,
        "baseline_residuals": baseline.residuals,
        "singular_values": [float(value) for value in svals.tolist()],
        "rank": rank,
        "irreducible_linearized_norm": irreducible_linear_norm,
        "best_actual_residual_norm": best_eval.norm,
        "mu_eta_sign_flip_residual_norm": mu_eta_sign_flip.norm,
        "mu_eta_sign_flip_delta_norm": mu_eta_sign_flip_delta_norm,
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_freeze_payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    print("STEP 20 REDUCED BRANCH FAMILY SCAN AUDIT")
    print("Built a target-blind local branch-family scan around the frozen step-19 reduced packet.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  geometry_status =", branch_metadata["geometry_status"])
    print("Frozen baseline packet (matches step-19 4-mode audit):")
    print("  M_sigma, K_sigma =", baseline.M_sigma, baseline.K_sigma)
    print("  B0, B2, B4 =", baseline.B)
    print("  Z0, Z2, Z4 =", baseline.Z)
    print("  N0, N2, N4 =", baseline.N)
    print("  residual packet (R_pole, R_norm, R_P2, R_P4) =", baseline.residuals)
    print("  residual norm =", baseline.norm)
    print("Local log-parameter family:")
    print("  coordinates = (log R0 amplitude, log R0 inverse-width, log beta inverse-width)")
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
        print(
            "  scale =",
            scale,
            "coords =",
            tuple(float(x) for x in (base_coords + scale * delta_ls)),
            "norm =",
            trial_eval.norm,
            "residuals =",
            trial_eval.residuals,
        )
    print("  best actual trial scale =", best_scale)
    print("  best actual residual norm =", best_eval.norm)
    print("  best actual residual packet =", best_eval.residuals)
    print("Sign-flip mutation guard:")
    print("  mu_eta sign-flip residual norm =", mu_eta_sign_flip.norm)
    print("  mu_eta sign-flip norm separation =", mu_eta_sign_flip.norm - baseline.norm)
    print("  mu_eta sign-flip residual-packet delta =", mu_eta_sign_flip_delta_norm)
    print("Interpretation:")
    print("  This scan is diagnostic only: it varies upstream reduced branch coordinates before the target check.")
    print("  A nonzero irreducible linearized norm means this 3-parameter family cannot cancel all four isotropic residuals at first order near the frozen branch.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
