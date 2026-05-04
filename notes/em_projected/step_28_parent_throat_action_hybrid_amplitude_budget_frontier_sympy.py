#!/usr/bin/env python3
"""Budgeted frontier for the exact Stage95 branch-B amplitude interpretation."""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass

import numpy as np

import step_23_parent_throat_action_normalization_frontier_sympy as step23


def assert_close(label: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label} failed: {actual} vs {expected} (tol={tol})")


@dataclass(frozen=True)
class HybridPoint:
    scale: float
    lambda_out: float
    sigma: float
    mhat0_req: float
    q_iso: float
    packet: tuple[float, float, float, float]


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_hybrid_amplitude_budget_frontier",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "step24 lambda_out reinterpreted as the exact Stage95 branch-B amplitude factor 1 - sigma",
        "sample_scales": [0.0, 0.03, 0.05, 0.08, 0.088, 0.09, 0.092],
        "sample_lambda_out": [1, 5, 20, 50, 100, 200, 500, 1000, 2000],
        "sigma_budgets": [4, 19, 49, 199, 1999],
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    base_coords = np.zeros(4, dtype=float)
    baseline_vec, jac = step23.central_jacobian(base_coords, 0.02)
    delta_ls, _, _, _ = np.linalg.lstsq(jac, -baseline_vec, rcond=None)

    points: list[HybridPoint] = []
    for scale in branch_metadata["sample_scales"]:
        Rpole, P0, P2, P4 = step23.evaluate_branch(*(base_coords + scale * delta_ls))
        if P0 <= 0.0:
            raise AssertionError(f"P0 became nonpositive at scale {scale}: {P0}")
        for lambda_out in branch_metadata["sample_lambda_out"]:
            p0_scaled = lambda_out * P0
            p2_scaled = lambda_out * P2
            p4_scaled = lambda_out * P4
            mhat0_req = math.sqrt((54.0 / 5.0) / p0_scaled)
            q_iso = float(np.linalg.norm(np.array([Rpole, p2_scaled, p4_scaled], dtype=float)))
            packet = (Rpole, p0_scaled - 54.0 / 5.0, p2_scaled, p4_scaled)
            points.append(
                HybridPoint(
                    scale=scale,
                    lambda_out=float(lambda_out),
                    sigma=1.0 - float(lambda_out),
                    mhat0_req=mhat0_req,
                    q_iso=q_iso,
                    packet=packet,
                )
            )

    best_q_by_budget: list[tuple[int, HybridPoint]] = []
    best_mhat_q1_by_budget: list[tuple[int, HybridPoint]] = []
    for budget in branch_metadata["sigma_budgets"]:
        allowed = [point for point in points if abs(point.sigma) <= budget]
        best_q = min(allowed, key=lambda point: point.q_iso)
        best_q_by_budget.append((budget, best_q))

        allowed_q1 = [point for point in allowed if point.q_iso <= 1.0]
        if not allowed_q1:
            raise AssertionError(f"No Q_iso <= 1 point found under sigma budget {budget}")
        best_mhat_q1 = min(allowed_q1, key=lambda point: point.mhat0_req)
        best_mhat_q1_by_budget.append((budget, best_mhat_q1))

    q_values = [point.q_iso for _, point in best_q_by_budget]
    mhat_values = [point.mhat0_req for _, point in best_mhat_q1_by_budget]
    if any(b > a for a, b in zip(q_values, q_values[1:])):
        raise AssertionError("Best Q_iso did not improve monotonically with relaxed sigma budgets")
    if any(b > a for a, b in zip(mhat_values, mhat_values[1:])):
        raise AssertionError("Best mhat0_req at Q_iso<=1 did not improve monotonically with relaxed sigma budgets")
    if any(point.lambda_out != 1.0 for _, point in best_q_by_budget):
        raise AssertionError("The best-Q frontier no longer stayed on the undeformed lambda_out = 1 branch")

    assert_close("best Q at sigma<=4", best_q_by_budget[0][1].q_iso, 0.26705543084121786, 1e-9)
    assert_close("best Q at sigma<=1999", best_q_by_budget[-1][1].q_iso, 0.26705543084121786, 1e-9)
    assert_close("best mhat q<=1 at sigma<=4", best_mhat_q1_by_budget[0][1].mhat0_req, 87.03141955817699, 1e-6)
    assert_close("best mhat q<=1 at sigma<=19", best_mhat_q1_by_budget[1][1].mhat0_req, 43.515709779088496, 1e-6)
    assert_close("best mhat q<=1 at sigma<=1999", best_mhat_q1_by_budget[-1][1].mhat0_req, 4.3515709779088505, 1e-6)
    assert_close("best q<=1 qiso at sigma<=4", best_mhat_q1_by_budget[0][1].q_iso, 0.4513189656743675, 1e-9)
    assert_close("best q<=1 qiso at sigma<=1999", best_mhat_q1_by_budget[-1][1].q_iso, 0.618690285150578, 1e-9)

    print("STEP 28 HYBRID AMPLITUDE BUDGET FRONTIER AUDIT")
    print("Reinterpreted the step-24 outgoing amplitude as the exact Stage95 branch-B factor 1 - sigma and scanned sigma budgets.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("  sample_scales =", branch_metadata["sample_scales"])
    print("  sample_lambda_out =", branch_metadata["sample_lambda_out"])
    print("  sigma_budgets =", branch_metadata["sigma_budgets"])
    print("Best sampled Q_iso by sigma budget:")
    for budget, point in best_q_by_budget:
        print(
            "  |sigma| <=",
            budget,
            "->",
            (point.scale, point.lambda_out, point.sigma, point.mhat0_req, point.q_iso),
        )
    print("Best sampled normalization among points with Q_iso <= 1 by sigma budget:")
    for budget, point in best_mhat_q1_by_budget:
        print(
            "  |sigma| <=",
            budget,
            "->",
            (point.scale, point.lambda_out, point.sigma, point.mhat0_req, point.q_iso),
        )
    print("Interpretation:")
    print("  This budgeted frontier asks whether the exact Stage95 branch-B amplitude channel helps before the deformation becomes extremely large.")
    print("  The best-Q frontier measures how much isotropic defect reduction is available at each amplitude budget.")
    print("  The Q_iso<=1 frontier measures how much normalization relief is available at each amplitude budget.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
