#!/usr/bin/env python3
"""Static-normalized diagnostic slice through the step-21 outgoing-family direction."""
from __future__ import annotations

import hashlib
import json
import math

import numpy as np
import sympy as sp
from scipy import integrate as scipy_integrate


def assert_close(label: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label} failed: {actual} vs {expected} (tol={tol})")


def assert_greater(label: str, actual: float, threshold: float) -> None:
    if not actual > threshold:
        raise AssertionError(f"{label} failed: {actual} <= {threshold}")


w = sp.symbols("w", real=True)
MODES = (0, 2, 4, 6)
# The basis/profile packets are Gaussian in the sampled range; finite bounds
# avoid unstable exp(+w^2) * exp(-w^2) cancellation probes at infinity.
QUAD_BOUND = 12.0


def basis_mode(n: int) -> sp.Expr:
    return sp.simplify(
        sp.hermite(n, w) * sp.exp(-w**2 / 2) / (sp.pi ** sp.Rational(1, 4) * sp.sqrt(2**n * sp.factorial(n)))
    )


BASIS = [basis_mode(n) for n in MODES]
BASIS_PRIME = [sp.diff(expr, w) for expr in BASIS]
BASIS_F = [sp.lambdify(w, expr, "numpy") for expr in BASIS]
BASIS_PRIME_F = [sp.lambdify(w, expr, "numpy") for expr in BASIS_PRIME]


def quad_inf(func, bound: float = QUAD_BOUND) -> float:
    value, _ = scipy_integrate.quad(func, -bound, bound, epsabs=1e-10, epsrel=1e-10, limit=200)
    return float(value)


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


def evaluate_branch(
    log_amp: float,
    log_r_width: float,
    log_beta_width: float,
    outgoing_delta: float,
    quad_bound: float = QUAD_BOUND,
) -> tuple[float, float, float, float]:
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

    mu_f = sp.lambdify(w, mu_eta, "numpy")
    tw_f = sp.lambdify(w, T_w, "numpy")
    to_f = sp.lambdify(w, T_omega, "numpy")
    keta_f = sp.lambdify(w, K_eta, "numpy")
    beta_f = sp.lambdify(w, beta, "numpy")
    beta_prime_f = sp.lambdify(w, sp.diff(beta, w), "numpy")
    phi_B_f = sp.lambdify(w, beta, "numpy")
    phi_Z_f = sp.lambdify(w, R0 * beta, "numpy")
    phi_N_f = sp.lambdify(w, (1 + outgoing_weight * w**2) * beta, "numpy")
    q = lambda func: quad_inf(func, bound=quad_bound)

    size = len(MODES)
    mass_matrix = symmetric_matrix(
        size,
        lambda i, j: q(lambda x: mu_f(x) * BASIS_F[i](x) * BASIS_F[j](x)),
    )
    stiff_matrix = symmetric_matrix(
        size,
        lambda i, j: q(
            lambda x: tw_f(x) * BASIS_PRIME_F[i](x) * BASIS_PRIME_F[j](x)
            + (keta_f(x) + 6.0 * to_f(x)) * BASIS_F[i](x) * BASIS_F[j](x)
        ),
    )

    coeffs_B = source_vector(size, lambda i: q(lambda x: BASIS_F[i](x) * phi_B_f(x)))
    coeffs_Z = source_vector(size, lambda i: q(lambda x: BASIS_F[i](x) * phi_Z_f(x)))
    coeffs_N = source_vector(size, lambda i: q(lambda x: BASIS_F[i](x) * phi_N_f(x)))

    B0, B2, B4 = low_series_packet(mass_matrix, stiff_matrix, coeffs_B)
    Z0, Z2, Z4 = low_series_packet(mass_matrix, stiff_matrix, coeffs_Z)
    N0, N2, N4 = low_series_packet(mass_matrix, stiff_matrix, coeffs_N)

    M_sigma = q(lambda x: mu_f(x) * beta_f(x) ** 2)
    K_sigma = q(lambda x: tw_f(x) * beta_prime_f(x) ** 2 + (keta_f(x) + 6.0 * to_f(x)) * beta_f(x) ** 2)

    D0 = K_sigma - B0 - Z0
    D2 = -(M_sigma + B2 + Z2)
    D4 = -(B4 + Z4)
    P0 = N0 / D0
    P2 = (D0 * N2 - 2.0 * D2 * N0) / D0**2
    P4 = (D0**2 * N4 - 2.0 * D0 * (D2 * N2 + D4 * N0) + 3.0 * D2**2 * N0) / D0**3
    R_pole = D0 * (B4 + Z4) - 3.0 * (M_sigma + B2 + Z2) ** 2
    return R_pole, P0, P2, P4


def central_jacobian(base: np.ndarray, h: float) -> np.ndarray:
    columns: list[np.ndarray] = []
    for idx in range(base.size):
        delta = np.zeros_like(base)
        delta[idx] = h
        plus = evaluate_branch(*(base + delta))
        minus = evaluate_branch(*(base - delta))
        plus_vec = np.array([plus[0], plus[1] - 54.0 / 5.0, plus[2], plus[3]], dtype=float)
        minus_vec = np.array([minus[0], minus[1] - 54.0 / 5.0, minus[2], minus[3]], dtype=float)
        columns.append((plus_vec - minus_vec) / (2.0 * h))
    return np.column_stack(columns)


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_static_normalized_slice",
        "pre_target_freeze": False,
        "target_blind": False,
        "no_post_residual_refit": False,
        "declared_slice": "R_norm == 0",
        "boundary_class": "open_impedance_demo",
        "geometry_status": "diagnostic_slice_through_step21_outgoing_family",
        "family_coordinates": [
            "log_R0_amplitude",
            "log_R0_inverse_width",
            "log_beta_inverse_width",
            "delta_outgoing_quadratic_weight",
        ],
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    base_coords = np.zeros(4, dtype=float)
    h = 0.02
    baseline = evaluate_branch(*base_coords)
    quad_bound_packets = np.array([evaluate_branch(*base_coords, quad_bound=bound) for bound in (10.0, 12.0, 14.0)], dtype=float)
    quad_bound_packet_spread = float(np.max(np.ptp(quad_bound_packets, axis=0)))
    if not quad_bound_packet_spread < 1.0e-9:
        raise AssertionError(f"quadrature-bound sensitivity too large: {quad_bound_packet_spread}")
    baseline_vec = np.array([baseline[0], baseline[1] - 54.0 / 5.0, baseline[2], baseline[3]], dtype=float)
    jac = central_jacobian(base_coords, h)
    delta_ls, _, jac_rank, jac_singular_values = np.linalg.lstsq(jac, -baseline_vec, rcond=None)
    if jac_rank != jac.shape[1]:
        raise AssertionError(f"linearized packet Jacobian is rank deficient: rank={jac_rank}, singular values={jac_singular_values}")
    jac_condition_number = float(jac_singular_values[0] / jac_singular_values[-1])
    if not jac_condition_number < 1.0e4:
        raise AssertionError(f"linearized packet Jacobian is poorly conditioned: {jac_condition_number}")
    linear_baseline_norm = float(np.linalg.norm(baseline_vec))
    linear_correct_norm = float(np.linalg.norm(baseline_vec + jac @ delta_ls))
    linear_sign_flip_norm = float(np.linalg.norm(baseline_vec - jac @ delta_ls))
    assert_greater("sign-flipped outgoing direction worsens linearized packet", linear_sign_flip_norm, linear_baseline_norm)

    lo = 0.09
    hi = 0.092
    rp_lo = evaluate_branch(*(base_coords + lo * delta_ls))[0]
    rp_hi = evaluate_branch(*(base_coords + hi * delta_ls))[0]
    if not (rp_lo < 0.0 < rp_hi):
        raise AssertionError("step-22 bracket does not straddle the one-pole root")

    for _ in range(12):
        mid = 0.5 * (lo + hi)
        rp_mid = evaluate_branch(*(base_coords + mid * delta_ls))[0]
        if rp_mid > 0.0:
            hi = mid
        else:
            lo = mid

    scale = 0.5 * (lo + hi)
    Rpole, P0, P2, P4 = evaluate_branch(*(base_coords + scale * delta_ls))
    mhat0_static = math.sqrt((54.0 / 5.0) / P0)
    packet = (Rpole, 0.0, P2, P4)
    packet_norm = float(np.linalg.norm(np.array(packet, dtype=float)))

    assert_close("baseline R_pole", baseline[0], -13.134593938872369, 1e-8)
    if not (0.091 < scale < 0.092):
        raise AssertionError(f"root scale left expected bracket: {scale}")
    assert abs(Rpole) < 5e-4
    assert abs(P2) < 5e-4
    assert abs(P4) < 5e-4
    assert packet_norm < 5e-4
    assert mhat0_static > 100.0

    print("STEP 22 STATIC-NORMALIZED SLICE AUDIT")
    print("Declared a non-target-blind diagnostic slice through the step-21 outgoing-family direction.")
    print("V2 diagnostic metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  declared_slice =", branch_metadata["declared_slice"])
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  geometry_status =", branch_metadata["geometry_status"])
    print("Baseline residual packet from the underlying target-blind family:")
    print("  (R_pole, R_norm, R_P2, R_P4) =", tuple(float(x) for x in baseline_vec))
    print("  quadrature-bound packet spread =", quad_bound_packet_spread)
    print("  QUAD_BOUND insensitivity guard = PASS")
    print("Step-21 least-squares outgoing-family direction:")
    print("  delta =", [round(float(x), 12) for x in delta_ls.tolist()])
    print("  jacobian condition number =", jac_condition_number)
    print("  sign-flipped linearized direction guard = PASS")
    print("One-pole root bracket on the static-normalized slice:")
    print("  bracket = (0.09, 0.092)")
    print("  root scale =", scale)
    print("Static-normalized slice data at the root:")
    print("  mhat0_static =", mhat0_static)
    print("  packet (R_pole, R_norm, R_P2, R_P4) =", packet)
    print("  packet norm =", packet_norm)
    print("Interpretation:")
    print("  On this declared static-normalized slice, the isotropic packet can be made numerically tiny inside the reduced outgoing family.")
    print("  But doing so requires a very large source normalization mhat0, so this is a residual-isolation diagnostic, not a realized-branch success verdict.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
