#!/usr/bin/env python3
"""Exact transfer/outlet amplitude interpretation for the step-24 lambda_out coordinate."""
from __future__ import annotations

import hashlib
import json

import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr != 0:
        raise AssertionError(f"{label} failed: {expr}")


def assert_close(label: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label} failed: {actual} vs {expected} (tol={tol})")


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_transfer_amplitude_interpretation",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "lambda_out compared against exact transfer-shape scaling and the stage95 hybrid outlet branch-B amplitude coordinate",
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    lam = sp.symbols("lambda_out", positive=True, real=True)
    D0, D2, D4 = sp.symbols("D0 D2 D4", nonzero=True, real=True)
    N0, N2, N4 = sp.symbols("N0 N2 N4", real=True)
    A, B, G5 = sp.symbols("A B G5", real=True)
    K, T2 = sp.symbols("K T2", positive=True, real=True)
    sigma, gamma = sp.symbols("sigma gamma", real=True)

    # Stage-005 internal prefactor packet.
    P0 = sp.simplify(N0 / D0)
    P2 = sp.simplify((D0 * N2 - 2 * D2 * N0) / D0**2)
    P4 = sp.simplify((D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3)

    P0_scaled = sp.simplify(P0.subs({N0: lam * N0, N2: lam * N2, N4: lam * N4}))
    P2_scaled = sp.simplify(P2.subs({N0: lam * N0, N2: lam * N2, N4: lam * N4}))
    P4_scaled = sp.simplify(P4.subs({N0: lam * N0, N2: lam * N2, N4: lam * N4}))

    assert_zero("uniform N-scaling gives uniform P0 scaling", P0_scaled - lam * P0)
    assert_zero("uniform N-scaling gives uniform P2 scaling", P2_scaled - lam * P2)
    assert_zero("uniform N-scaling gives uniform P4 scaling", P4_scaled - lam * P4)

    # Stage-005 grouped outgoing compiler.
    K0 = sp.simplify(P0)
    K2 = sp.simplify(P2 + A * P0)
    K4 = sp.simplify(P4 + A * P2 + B * P0)
    Gamma5 = sp.simplify(G5 * P0)
    K0_scaled = sp.simplify(P0_scaled)
    K2_scaled = sp.simplify(P2_scaled + A * P0_scaled)
    K4_scaled = sp.simplify(P4_scaled + A * P2_scaled + B * P0_scaled)
    Gamma5_scaled = sp.simplify(G5 * P0_scaled)

    assert_zero("compiler preserves uniform K0 scaling", K0_scaled - lam * K0)
    assert_zero("compiler preserves uniform K2 scaling", K2_scaled - lam * K2)
    assert_zero("compiler preserves uniform K4 scaling", K4_scaled - lam * K4)
    assert_zero("compiler preserves uniform Gamma5 scaling", Gamma5_scaled - lam * Gamma5)

    K2_shape = sp.simplify(K2 / K0)
    K4_shape = sp.simplify(K4 / K0)
    assert_zero(
        "uniform scaling leaves K2/K0 invariant",
        sp.simplify(K2_scaled / K0_scaled - K2_shape),
    )
    assert_zero(
        "uniform scaling leaves K4/K0 invariant",
        sp.simplify(K4_scaled / K0_scaled - K4_shape),
    )

    # Stage-163 transfer-shape theorem.
    transfer_shape_scaled = sp.simplify((lam * N0) / K)
    assert_zero("transfer-shape amplitude law", transfer_shape_scaled - lam * (N0 / K))

    # Stage-95 hybrid branch B.
    chi_B = sp.simplify((1 - 9 * sigma * gamma) / (1 - sigma))
    chi_B_canonical = sp.simplify(chi_B.subs(gamma, sp.Rational(1, 9)))
    scaled_identity_factor = sp.simplify(1 - sigma)
    sigma_from_lambda = sp.solve(sp.Eq(scaled_identity_factor, lam), sigma)[0]

    assert_zero("stage95 canonical branch-B keeps chi_Q = 1", chi_B_canonical - 1)
    assert_zero("lambda_out to hybrid-branch sigma map", sigma_from_lambda - (1 - lam))

    point_budget50 = {
        "label": "step24 best point under mhat0_req <= 50",
        "scale": 0.092,
        "lambda_out": 20.0,
        "mhat0_req": 47.912441136331765,
        "Q_iso": 0.2670781822626949,
    }
    point_q1 = {
        "label": "step24 best point with Q_iso <= 1",
        "scale": 0.09,
        "lambda_out": 2000.0,
        "mhat0_req": 4.351570977913287,
        "Q_iso": 0.618690285150578,
    }
    point_qhalf = {
        "label": "step24 best point with Q_iso <= 0.5",
        "scale": 0.092,
        "lambda_out": 2000.0,
        "mhat0_req": 4.7912441136331765,
        "Q_iso": 0.4394839373049669,
    }

    def sigma_branch_b(point: dict[str, float]) -> float:
        return 1.0 - point["lambda_out"]

    sigma_budget50 = sigma_branch_b(point_budget50)
    sigma_q1 = sigma_branch_b(point_q1)
    sigma_qhalf = sigma_branch_b(point_qhalf)

    assert_close("budget50 sigma", sigma_budget50, -19.0, 1e-12)
    assert_close("q<=1 sigma", sigma_q1, -1999.0, 1e-12)
    assert_close("q<=0.5 sigma", sigma_qhalf, -1999.0, 1e-12)

    print("STEP 27 TRANSFER-AMPLITUDE INTERPRETATION AUDIT")
    print("Tested whether the step-24 lambda_out coordinate has an exact PDE-side interpretation independent of N_Q.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Exact prefactor / compiler scaling laws:")
    print("  P0(lambda_out*N) =", sp.sstr(P0_scaled))
    print("  P2(lambda_out*N) =", sp.sstr(P2_scaled))
    print("  P4(lambda_out*N) =", sp.sstr(P4_scaled))
    print("  uniform scaling leaves K2/K0 invariant with K2/K0 =", sp.sstr(K2_shape))
    print("  uniform scaling leaves K4/K0 invariant with K4/K0 =", sp.sstr(K4_shape))
    print("  exact transfer-shape amplitude law = T_eff^2 -> lambda_out * T_eff^2")
    print("Stage95 hybrid branch-B amplitude coordinate:")
    print("  chi_B(sigma, gamma) =", sp.sstr(chi_B))
    print("  canonical branch-B chi_Q =", sp.sstr(chi_B_canonical))
    print("  branch-B scaling factor =", sp.sstr(scaled_identity_factor))
    print("  sigma(lambda_out) =", sp.sstr(sigma_from_lambda))
    print("Concrete step24 mappings:")
    for point, sigma_value in (
        (point_budget50, sigma_budget50),
        (point_q1, sigma_q1),
        (point_qhalf, sigma_qhalf),
    ):
        print(" ", point["label"], "=", point)
        print("    hybrid branch-B sigma =", sigma_value)
    print("Interpretation:")
    print("  Step24 lambda_out is not admissible as N_Q, but exact PDE algebra does leave room for an independent transfer/outlet amplitude coordinate.")
    print("  Uniform scaling of the outgoing N-packet scales the whole P-packet and the compiled branch coefficients without changing the normalized shape ratios.")
    print("  Stage95 branch B gives a concrete exact outlet family where chi_Q stays canonical while the overall outgoing amplitude scales by 1 - sigma.")
    print("  The remaining burden is quantitative, not algebraic: the useful step24 points correspond to sigma = -19 or sigma = -1999 in that minimal family.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
