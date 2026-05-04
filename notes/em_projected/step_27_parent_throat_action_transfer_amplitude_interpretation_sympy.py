#!/usr/bin/env python3
"""Exact transfer/outlet amplitude interpretation for the step-24 lambda_out coordinate."""
from __future__ import annotations

import hashlib
import json

import sympy as sp

from step_24_parent_throat_action_outgoing_amplitude_frontier_sympy import export_step24_frontier_points
from step_29_parent_throat_action_moderate_branchb_sector_sympy import export_step29_moderate_branchb_patch


def assert_zero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr != 0:
        raise AssertionError(f"{label} failed: {expr}")


def assert_close(label: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label} failed: {actual} vs {expected} (tol={tol})")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


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
    eps = sp.symbols("eps", real=True, nonzero=True)

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
    assert_nonzero("uniform N-scaling detects mutated lambda", P0_scaled - (lam + eps) * P0)

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
    gamma_canonical = sp.Rational(1, 9)
    # Sanity-check the imported Stage-95 value locally; this is not an upstream
    # derivation of gamma_canonical.
    chi_B_canonical = sp.simplify(chi_B.subs(gamma, gamma_canonical))
    scaled_identity_factor = sp.simplify(1 - sigma)
    sigma_from_lambda = sp.solve(sp.Eq(scaled_identity_factor, lam), sigma)[0]

    assert_zero("lambda_out to hybrid-branch sigma map", sigma_from_lambda - (1 - lam))
    assert_nonzero("lambda_out sigma map detects shifted branch factor", sigma_from_lambda - (1 - (lam + eps)))

    step24_exports = export_step24_frontier_points()
    step29_patch = export_step29_moderate_branchb_patch()
    patch_20 = step29_patch["same_scale_lambda20"]
    point_budget50 = {
        "label": "step24 best point under mhat0_req <= 50",
        "scale": float(patch_20["scale"]),
        "lambda_out": float(patch_20["lambda_out"]),
        "mhat0_req": float(patch_20["mhat0_req"]),
        "Q_iso": float(patch_20["Q_iso"]),
    }
    point_q1 = {
        "label": "step24 best point with Q_iso <= 1",
        **step24_exports["best_q_iso_le_1"],
    }
    point_qhalf = {
        "label": "step24 best point with Q_iso <= 0.5",
        **step24_exports["best_q_iso_le_half"],
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
    print("  scaling mutation guards = PASS")
    print("Stage95 hybrid branch-B amplitude coordinate:")
    print("  canonical branch-B gamma =", sp.sstr(gamma_canonical))
    print("  gamma provenance = imported Stage-95 canonical value; no independent upstream derivation in this local script")
    print("  chi_B(sigma, gamma) =", sp.sstr(chi_B))
    print("  canonical branch-B chi_Q =", sp.sstr(chi_B_canonical))
    print("  canonical chi_Q check is a local sanity check, not an upstream gamma derivation")
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
