#!/usr/bin/env python3
"""Target-blind normalization frontier scan on the step-21 outgoing family."""
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
class SamplePoint:
    scale: float
    p0: float
    mhat0_req: float
    q_iso: float
    full_norm: float
    packet: tuple[float, float, float, float]


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


def evaluate_branch(log_amp: float, log_r_width: float, log_beta_width: float, outgoing_delta: float) -> tuple[float, float, float, float]:
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
    return R_pole, P0, P2, P4


def central_jacobian(base: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    baseline = evaluate_branch(*base)
    baseline_vec = np.array([baseline[0], baseline[1] - 54.0 / 5.0, baseline[2], baseline[3]], dtype=float)
    columns: list[np.ndarray] = []
    for idx in range(base.size):
        delta = np.zeros_like(base)
        delta[idx] = h
        plus = evaluate_branch(*(base + delta))
        minus = evaluate_branch(*(base - delta))
        plus_vec = np.array([plus[0], plus[1] - 54.0 / 5.0, plus[2], plus[3]], dtype=float)
        minus_vec = np.array([minus[0], minus[1] - 54.0 / 5.0, minus[2], minus[3]], dtype=float)
        columns.append((plus_vec - minus_vec) / (2.0 * h))
    return baseline_vec, np.column_stack(columns)


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_normalization_frontier_scan",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "geometry_status": "target_blind_frontier_scan_on_step21_ray",
        "family_coordinates": [
            "log_R0_amplitude",
            "log_R0_inverse_width",
            "log_beta_inverse_width",
            "delta_outgoing_quadratic_weight",
        ],
        "sample_scales": [0.0, 0.03, 0.05, 0.08, 0.088, 0.09, 0.092],
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    base_coords = np.zeros(4, dtype=float)
    h = 0.02
    baseline_vec, jac = central_jacobian(base_coords, h)
    delta_ls, _, _, _ = np.linalg.lstsq(jac, -baseline_vec, rcond=None)

    assert_close("baseline R_pole", baseline_vec[0], -13.134593938872369, 1e-8)
    assert_close("step21 outgoing delta[3]", float(delta_ls[3]), 74.190514652389, 1e-8)

    sample_points: list[SamplePoint] = []
    for scale in branch_metadata["sample_scales"]:
        Rpole, P0, P2, P4 = evaluate_branch(*(base_coords + scale * delta_ls))
        if P0 <= 0.0:
            raise AssertionError(f"P0 became nonpositive at scale {scale}: {P0}")
        mhat0_req = math.sqrt((54.0 / 5.0) / P0)
        q_iso = float(np.linalg.norm(np.array([Rpole, P2, P4], dtype=float)))
        full_norm = float(np.linalg.norm(np.array([Rpole, P0 - 54.0 / 5.0, P2, P4], dtype=float)))
        sample_points.append(
            SamplePoint(
                scale=scale,
                p0=P0,
                mhat0_req=mhat0_req,
                q_iso=q_iso,
                full_norm=full_norm,
                packet=(Rpole, P0 - 54.0 / 5.0, P2, P4),
            )
        )

    frontier = []
    for point in sample_points:
        dominated = False
        for other in sample_points:
            if other is point:
                continue
            if other.q_iso <= point.q_iso and other.mhat0_req <= point.mhat0_req and (
                other.q_iso < point.q_iso or other.mhat0_req < point.mhat0_req
            ):
                dominated = True
                break
        if not dominated:
            frontier.append(point)

    qiso_values = [point.q_iso for point in sample_points]
    mhat_values = [point.mhat0_req for point in sample_points]
    if any(b > a for a, b in zip(qiso_values, qiso_values[1:])):
        raise AssertionError("Q_iso was not monotone decreasing along the frozen sample ray")
    if any(b < a for a, b in zip(mhat_values, mhat_values[1:])):
        raise AssertionError("mhat0_req was not monotone increasing along the frozen sample ray")

    first_subunit = next(point for point in sample_points if point.q_iso < 1.0)
    best_qiso = min(sample_points, key=lambda point: point.q_iso)

    assert_close("first subunit-qiso scale", first_subunit.scale, 0.09, 1e-12)
    assert_close("first subunit-qiso mhat0", first_subunit.mhat0_req, 194.6081703105869, 1e-6)
    assert_close("best sampled qiso scale", best_qiso.scale, 0.092, 1e-12)
    assert_close("best sampled qiso", best_qiso.q_iso, 0.26705543084121786, 1e-9)

    print("STEP 23 NORMALIZATION FRONTIER SCAN AUDIT")
    print("Built a target-blind normalization frontier on a frozen sample ray inside the step-21 outgoing family.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  geometry_status =", branch_metadata["geometry_status"])
    print("Underlying step-21 least-squares direction:")
    print("  delta =", [round(float(x), 12) for x in delta_ls.tolist()])
    print("Frozen sample ray:")
    print("  sample_scales =", branch_metadata["sample_scales"])
    print("Sampled normalization frontier points:")
    for point in sample_points:
        print(
            "  scale =",
            point.scale,
            "mhat0_req =",
            point.mhat0_req,
            "Q_iso =",
            point.q_iso,
            "full_norm =",
            point.full_norm,
            "packet =",
            point.packet,
        )
    print("Pareto frontier on the sampled ray:")
    for point in frontier:
        print("  frontier scale =", point.scale, "mhat0_req =", point.mhat0_req, "Q_iso =", point.q_iso)
    print("Key frontier diagnostics:")
    print("  first sampled point with Q_iso < 1 =", (first_subunit.scale, first_subunit.mhat0_req, first_subunit.q_iso))
    print("  best sampled Q_iso point =", (best_qiso.scale, best_qiso.mhat0_req, best_qiso.q_iso))
    print("Interpretation:")
    print("  Along this frozen target-blind ray, reducing the isotropic packet requires rapidly increasing source normalization.")
    print("  On the sampled ray, Q_iso first drops below 1 only once mhat0_req is already about 1.95e2.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
