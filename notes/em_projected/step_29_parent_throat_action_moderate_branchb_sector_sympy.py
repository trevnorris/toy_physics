#!/usr/bin/env python3
"""Moderate Stage95 branch-B sector on the sampled reduced frontier."""
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


def assert_greater(label: str, actual: float, threshold: float) -> None:
    if not actual > threshold:
        raise AssertionError(f"{label} failed: {actual} <= {threshold}")


@dataclass(frozen=True)
class Point:
    scale: float
    lambda_out: float
    sigma: float
    mhat0_req: float
    q_iso: float
    delta_ln_t2: float


STEP29_MODERATE_BRANCHB_EXPORTS = {
    "P0_base": 0.00023523241237876055,
    "same_scale_lambda20": {
        "scale": 0.092,
        "lambda_out": 20.0,
        "sigma": -19.0,
        "mhat0_req": 47.91244113628207,
        "Q_iso": 0.2670781822533451,
        "delta_ln_t2": 2.995732273553991,
    },
    "same_scale_lambda50": {
        "scale": 0.092,
        "lambda_out": 50.0,
        "sigma": -49.0,
        "mhat0_req": 30.302488449879455,
        "Q_iso": 0.26719789464729127,
        "delta_ln_t2": 3.912023005428146,
    },
}


def export_step29_moderate_branchb_patch() -> dict[str, object]:
    return {
        "P0_base": STEP29_MODERATE_BRANCHB_EXPORTS["P0_base"],
        "same_scale_lambda20": dict(STEP29_MODERATE_BRANCHB_EXPORTS["same_scale_lambda20"]),
        "same_scale_lambda50": dict(STEP29_MODERATE_BRANCHB_EXPORTS["same_scale_lambda50"]),
    }


def pareto_frontier(points: list[Point]) -> list[Point]:
    frontier: list[Point] = []
    for point in points:
        dominated = False
        for other in points:
            if other is point:
                continue
            if other.mhat0_req <= point.mhat0_req and other.q_iso <= point.q_iso and (
                other.mhat0_req < point.mhat0_req or other.q_iso < point.q_iso
            ):
                dominated = True
                break
        if not dominated:
            frontier.append(point)
    return sorted(frontier, key=lambda point: (point.q_iso, point.mhat0_req))


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_moderate_branchb_sector",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "sampled low-defect Pareto sector under the exact Stage95 branch-B amplitude coordinate",
        "sample_scales": [0.0, 0.03, 0.05, 0.08, 0.088, 0.09, 0.092],
        "sample_lambda_out": [1, 5, 20, 50, 100, 200, 500, 1000, 2000],
        "sigma_budgets": [19, 49],
        "q_sector_cap": 1.5,
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    base_coords = np.zeros(4, dtype=float)
    baseline_vec, jac = step23.central_jacobian(base_coords, 0.02)
    delta_ls, _, _, _ = np.linalg.lstsq(jac, -baseline_vec, rcond=None)
    linear_baseline_norm = float(np.linalg.norm(baseline_vec))
    linear_correct_norm = float(np.linalg.norm(baseline_vec + jac @ delta_ls))
    linear_sign_flip_norm = float(np.linalg.norm(baseline_vec - jac @ delta_ls))
    if not linear_correct_norm < 1e-8:
        raise AssertionError(f"least-squares direction no longer closes the linearized packet: {linear_correct_norm}")
    assert_greater("sign-flipped outgoing direction worsens linearized packet", linear_sign_flip_norm, linear_baseline_norm)

    all_points: list[Point] = []
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
            all_points.append(
                Point(
                    scale=scale,
                    lambda_out=float(lambda_out),
                    sigma=1.0 - float(lambda_out),
                    mhat0_req=mhat0_req,
                    q_iso=q_iso,
                    delta_ln_t2=math.log(float(lambda_out)),
                )
            )

    frontiers: dict[int, list[Point]] = {}
    for budget in branch_metadata["sigma_budgets"]:
        sector_points = [
            point
            for point in all_points
            if abs(point.sigma) <= budget and point.q_iso <= branch_metadata["q_sector_cap"]
        ]
        frontiers[budget] = pareto_frontier(sector_points)

    if len(frontiers[19]) != 5:
        raise AssertionError(f"Unexpected low-defect frontier size for sigma<=19: {len(frontiers[19])}")
    if len(frontiers[49]) != 6:
        raise AssertionError(f"Unexpected low-defect frontier size for sigma<=49: {len(frontiers[49])}")

    baseline = next(point for point in frontiers[19] if point.scale == 0.092 and point.lambda_out == 1.0)
    lam20_same_scale = next(point for point in frontiers[19] if point.scale == 0.092 and point.lambda_out == 20.0)
    lam20_shifted = next(point for point in frontiers[19] if point.scale == 0.09 and point.lambda_out == 20.0)
    lam50_same_scale = next(point for point in frontiers[49] if point.scale == 0.092 and point.lambda_out == 50.0)
    export_20 = STEP29_MODERATE_BRANCHB_EXPORTS["same_scale_lambda20"]
    export_50 = STEP29_MODERATE_BRANCHB_EXPORTS["same_scale_lambda50"]
    p0_base_export = STEP29_MODERATE_BRANCHB_EXPORTS["P0_base"]

    q_delta_20_same = lam20_same_scale.q_iso - baseline.q_iso
    q_rel_20_same = q_delta_20_same / baseline.q_iso
    norm_drop_20_same = baseline.mhat0_req / lam20_same_scale.mhat0_req

    q_delta_50_same = lam50_same_scale.q_iso - baseline.q_iso
    q_rel_50_same = q_delta_50_same / baseline.q_iso
    norm_drop_50_same = baseline.mhat0_req / lam50_same_scale.mhat0_req
    p0_base_from_20 = (54.0 / 5.0) / (lam20_same_scale.mhat0_req**2 * lam20_same_scale.lambda_out)
    p0_base_from_50 = (54.0 / 5.0) / (lam50_same_scale.mhat0_req**2 * lam50_same_scale.lambda_out)

    assert_close("baseline q", baseline.q_iso, 0.2670554308318671, 1e-9)
    assert_close("exported P0_base from lambda20", p0_base_from_20, p0_base_export, 1e-18)
    assert_close("exported P0_base from lambda50", p0_base_from_50, p0_base_export, 1e-18)
    assert_close("same-scale lambda20 scale", lam20_same_scale.scale, export_20["scale"], 1e-12)
    assert_close("same-scale lambda20 lambda", lam20_same_scale.lambda_out, export_20["lambda_out"], 1e-12)
    assert_close("same-scale lambda20 sigma", lam20_same_scale.sigma, export_20["sigma"], 1e-12)
    assert_close("same-scale lambda20 q", lam20_same_scale.q_iso, export_20["Q_iso"], 1e-9)
    assert_close("same-scale lambda20 mhat", lam20_same_scale.mhat0_req, export_20["mhat0_req"], 1e-6)
    assert_close("same-scale lambda20 q delta", q_delta_20_same, 2.2751421478006684e-05, 1e-12)
    assert_close("same-scale lambda20 q relative", q_rel_20_same, 8.519362967881576e-05, 1e-12)
    assert_close("same-scale lambda20 normalization drop", norm_drop_20_same, math.sqrt(20.0), 1e-12)
    assert_close("same-scale lambda20 log amplitude", lam20_same_scale.delta_ln_t2, export_20["delta_ln_t2"], 1e-12)
    sigma_sign_flip_residual_20 = abs(lam20_same_scale.sigma - (lam20_same_scale.lambda_out - 1.0))
    inverse_lambda_mutation_residual_20 = abs((lam20_same_scale.mhat0_req / baseline.mhat0_req) ** 2 - lam20_same_scale.lambda_out)
    assert_close("lambda20 branch-B sigma law", lam20_same_scale.sigma, 1.0 - lam20_same_scale.lambda_out, 1e-12)
    assert_greater("lambda20 sigma sign-flip mutation residual", sigma_sign_flip_residual_20, 10.0)
    assert_greater("lambda20 inverse lambda scaling mutation residual", inverse_lambda_mutation_residual_20, 10.0)

    assert_close("shifted lambda20 q", lam20_shifted.q_iso, 0.4513375659863663, 1e-9)
    assert_close("shifted lambda20 mhat", lam20_shifted.mhat0_req, 43.51570977908849, 1e-6)

    assert_close("same-scale lambda50 scale", lam50_same_scale.scale, export_50["scale"], 1e-12)
    assert_close("same-scale lambda50 lambda", lam50_same_scale.lambda_out, export_50["lambda_out"], 1e-12)
    assert_close("same-scale lambda50 sigma", lam50_same_scale.sigma, export_50["sigma"], 1e-12)
    assert_close("same-scale lambda50 q", lam50_same_scale.q_iso, export_50["Q_iso"], 1e-9)
    assert_close("same-scale lambda50 mhat", lam50_same_scale.mhat0_req, export_50["mhat0_req"], 1e-6)
    assert_close("same-scale lambda50 q delta", q_delta_50_same, 0.0001424638154241542, 1e-12)
    assert_close("same-scale lambda50 q relative", q_rel_50_same, 0.0005334615925255107, 1e-12)
    assert_close("same-scale lambda50 normalization drop", norm_drop_50_same, math.sqrt(50.0), 1e-12)
    assert_close("same-scale lambda50 log amplitude", lam50_same_scale.delta_ln_t2, export_50["delta_ln_t2"], 1e-12)
    sigma_sign_flip_residual_50 = abs(lam50_same_scale.sigma - (lam50_same_scale.lambda_out - 1.0))
    inverse_lambda_mutation_residual_50 = abs((lam50_same_scale.mhat0_req / baseline.mhat0_req) ** 2 - lam50_same_scale.lambda_out)
    assert_close("lambda50 branch-B sigma law", lam50_same_scale.sigma, 1.0 - lam50_same_scale.lambda_out, 1e-12)
    assert_greater("lambda50 sigma sign-flip mutation residual", sigma_sign_flip_residual_50, 50.0)
    assert_greater("lambda50 inverse lambda scaling mutation residual", inverse_lambda_mutation_residual_50, 25.0)

    print("STEP 29 MODERATE BRANCH-B SECTOR AUDIT")
    print("Extracted the sampled low-defect Pareto sector for moderate exact Stage95 branch-B amplitude budgets.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("  sigma_budgets =", branch_metadata["sigma_budgets"])
    print("  q_sector_cap =", branch_metadata["q_sector_cap"])
    print("Low-defect Pareto frontier for |sigma| <= 19:")
    for point in frontiers[19]:
        print(" ", (point.scale, point.lambda_out, point.sigma, point.mhat0_req, point.q_iso, point.delta_ln_t2))
    print("Low-defect Pareto frontier for |sigma| <= 49:")
    for point in frontiers[49]:
        print(" ", (point.scale, point.lambda_out, point.sigma, point.mhat0_req, point.q_iso, point.delta_ln_t2))
    print("Direction and amplitude mutation guards:")
    print("  sign-flipped linearized direction guard = PASS")
    print("  sigma sign-flip mutation residuals =", (sigma_sign_flip_residual_20, sigma_sign_flip_residual_50))
    print("  inverse lambda scaling mutation residuals =", (inverse_lambda_mutation_residual_20, inverse_lambda_mutation_residual_50))
    print("Key moderate-sector diagnostics:")
    print("  baseline point =", (baseline.scale, baseline.lambda_out, baseline.mhat0_req, baseline.q_iso))
    print("  same-scale lambda_out = 20 point =", (lam20_same_scale.scale, lam20_same_scale.lambda_out, lam20_same_scale.mhat0_req, lam20_same_scale.q_iso))
    print("    q increase =", q_delta_20_same)
    print("    relative q increase =", q_rel_20_same)
    print("    normalization drop factor =", norm_drop_20_same)
    print("    delta ln(T_eff^2) =", lam20_same_scale.delta_ln_t2)
    print("  shifted scale lambda_out = 20 point =", (lam20_shifted.scale, lam20_shifted.lambda_out, lam20_shifted.mhat0_req, lam20_shifted.q_iso))
    print("  same-scale lambda_out = 50 point =", (lam50_same_scale.scale, lam50_same_scale.lambda_out, lam50_same_scale.mhat0_req, lam50_same_scale.q_iso))
    print("    q increase =", q_delta_50_same)
    print("    relative q increase =", q_rel_50_same)
    print("    normalization drop factor =", norm_drop_50_same)
    print("    delta ln(T_eff^2) =", lam50_same_scale.delta_ln_t2)
    print("Interpretation:")
    print("  The moderate branch-B sector is much stronger than the coarse step-28 slices suggested.")
    print("  At fixed scale 0.092, sigma = -19 cuts the normalization burden by sqrt(20) while increasing Q_iso by only about 2.28e-05.")
    print("  Even sigma = -49 cuts the normalization burden by sqrt(50) while increasing Q_iso by only about 1.42e-04.")
    print("  So the first physically interesting regime is not the extreme sigma = -1999 corner; it is the moderate sigma = -19 to -49 corridor.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
