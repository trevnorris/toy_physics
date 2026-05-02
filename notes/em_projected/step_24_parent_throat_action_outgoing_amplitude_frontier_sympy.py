#!/usr/bin/env python3
"""Target-blind outgoing-amplitude frontier on the frozen step-21 ray."""
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
class FrontierPoint:
    scale: float
    lambda_out: float
    mhat0_req: float
    q_iso: float
    packet: tuple[float, float, float, float]


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_outgoing_amplitude_frontier_scan",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "geometry_status": "target_blind_outgoing_amplitude_frontier_on_step21_ray",
        "sample_scales": [0.0, 0.03, 0.05, 0.08, 0.088, 0.09, 0.092],
        "sample_lambda_out": [1, 5, 20, 50, 100, 200, 500, 1000, 2000],
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    base_coords = np.zeros(4, dtype=float)
    baseline_vec, jac = step23.central_jacobian(base_coords, 0.02)
    delta_ls, _, _, _ = np.linalg.lstsq(jac, -baseline_vec, rcond=None)

    points: list[FrontierPoint] = []
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
            points.append(FrontierPoint(scale, float(lambda_out), mhat0_req, q_iso, packet))

    best_q_under_50 = min((point for point in points if point.mhat0_req <= 50.0), key=lambda point: point.q_iso)
    best_q_under_20 = min((point for point in points if point.mhat0_req <= 20.0), key=lambda point: point.q_iso)
    best_q_under_10 = min((point for point in points if point.mhat0_req <= 10.0), key=lambda point: point.q_iso)
    best_q_under_5 = min((point for point in points if point.mhat0_req <= 5.0), key=lambda point: point.q_iso)
    best_mhat_under_q1 = min((point for point in points if point.q_iso <= 1.0), key=lambda point: point.mhat0_req)
    best_mhat_under_qhalf = min((point for point in points if point.q_iso <= 0.5), key=lambda point: point.mhat0_req)

    assert_close("best q under 50 scale", best_q_under_50.scale, 0.092, 1e-12)
    assert_close("best q under 50 lambda", best_q_under_50.lambda_out, 20.0, 1e-12)
    assert_close("best q under 50 mhat0", best_q_under_50.mhat0_req, 47.912441136331765, 1e-6)
    assert_close("best q under 50 qiso", best_q_under_50.q_iso, 0.2670781822626949, 1e-9)

    assert_close("best q under 20 lambda", best_q_under_20.lambda_out, 200.0, 1e-12)
    assert_close("best q under 20 mhat0", best_q_under_20.mhat0_req, 15.151244224955443, 1e-6)
    assert_close("best q under 20 qiso", best_q_under_20.q_iso, 0.2693266571833183, 1e-9)

    assert_close("best q under 10 lambda", best_q_under_10.lambda_out, 500.0, 1e-12)
    assert_close("best q under 10 mhat0", best_q_under_10.mhat0_req, 9.582488227266353, 1e-6)
    assert_close("best q under 10 qiso", best_q_under_10.q_iso, 0.28094980884298526, 1e-9)

    assert_close("best q under 5 lambda", best_q_under_5.lambda_out, 2000.0, 1e-12)
    assert_close("best q under 5 mhat0", best_q_under_5.mhat0_req, 4.7912441136331765, 1e-6)
    assert_close("best q under 5 qiso", best_q_under_5.q_iso, 0.4394839373049669, 1e-9)

    assert_close("best mhat under q<=1 scale", best_mhat_under_q1.scale, 0.09, 1e-12)
    assert_close("best mhat under q<=1 lambda", best_mhat_under_q1.lambda_out, 2000.0, 1e-12)
    assert_close("best mhat under q<=1", best_mhat_under_q1.mhat0_req, 4.351570977913287, 1e-6)
    assert_close("best mhat under q<=1 qiso", best_mhat_under_q1.q_iso, 0.618690285150578, 1e-9)

    assert_close("best mhat under q<=0.5 scale", best_mhat_under_qhalf.scale, 0.092, 1e-12)
    assert_close("best mhat under q<=0.5 lambda", best_mhat_under_qhalf.lambda_out, 2000.0, 1e-12)
    assert_close("best mhat under q<=0.5", best_mhat_under_qhalf.mhat0_req, 4.7912441136331765, 1e-6)
    assert_close("best mhat under q<=0.5 qiso", best_mhat_under_qhalf.q_iso, 0.4394839373049669, 1e-9)

    print("STEP 24 OUTGOING AMPLITUDE FRONTIER AUDIT")
    print("Built a target-blind two-coordinate frontier using the frozen step-21 ray plus an outgoing-amplitude coordinate.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  geometry_status =", branch_metadata["geometry_status"])
    print("Frozen sampling grid:")
    print("  sample_scales =", branch_metadata["sample_scales"])
    print("  sample_lambda_out =", branch_metadata["sample_lambda_out"])
    print("Underlying step-21 least-squares direction:")
    print("  delta =", [round(float(x), 12) for x in delta_ls.tolist()])
    print("Best sampled Q_iso at fixed normalization budgets:")
    print("  mhat0_req <= 50  ->", (best_q_under_50.scale, best_q_under_50.lambda_out, best_q_under_50.mhat0_req, best_q_under_50.q_iso))
    print("  mhat0_req <= 20  ->", (best_q_under_20.scale, best_q_under_20.lambda_out, best_q_under_20.mhat0_req, best_q_under_20.q_iso))
    print("  mhat0_req <= 10  ->", (best_q_under_10.scale, best_q_under_10.lambda_out, best_q_under_10.mhat0_req, best_q_under_10.q_iso))
    print("  mhat0_req <= 5   ->", (best_q_under_5.scale, best_q_under_5.lambda_out, best_q_under_5.mhat0_req, best_q_under_5.q_iso))
    print("Best sampled normalization at isotropic-defect thresholds:")
    print("  Q_iso <= 1.0 ->", (best_mhat_under_q1.scale, best_mhat_under_q1.lambda_out, best_mhat_under_q1.mhat0_req, best_mhat_under_q1.q_iso))
    print("  Q_iso <= 0.5 ->", (best_mhat_under_qhalf.scale, best_mhat_under_qhalf.lambda_out, best_mhat_under_qhalf.mhat0_req, best_mhat_under_qhalf.q_iso))
    print("Interpretation:")
    print("  A target-blind outgoing amplitude coordinate changes the normalization frontier dramatically.")
    print("  On the sampled grid, Q_iso <= 1 is compatible with mhat0_req around 4.35, and Q_iso <= 0.5 with mhat0_req around 4.79.")
    print("  So static normalization no longer looks structurally fatal if this outgoing-amplitude family is physically admissible.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
