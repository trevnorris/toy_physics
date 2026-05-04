#!/usr/bin/env python3
"""Natural-source outgoing burden implied by the step-24 amplitude frontier."""
from __future__ import annotations

import hashlib
import json

import sympy as sp

from step_24_parent_throat_action_outgoing_amplitude_frontier_sympy import export_step24_frontier_points


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
        "branch_id": "v2_local_parent_background_natural_source_burden",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "step24 outgoing amplitude compared against the exact PDE natural source map and the minimal Robin outlet model",
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    P0_base, P0_target = sp.symbols("P0_base P0_target", positive=True, real=True)
    mhat0_req, lambda_out = sp.symbols("mhat0_req lambda_out", positive=True, real=True)
    chi_Q, N_Q_nat = sp.symbols("chi_Q N_Q_nat", positive=True, real=True)
    rho = sp.symbols("rho", real=True)
    eps = sp.symbols("eps", real=True, nonzero=True)

    # Reduced static-normalized frontier law from step 24.
    static_scaled = sp.Eq(mhat0_req**2 * lambda_out * P0_base, P0_target)

    # Exact PDE natural-source-map law.
    natural_source_map = sp.Eq(N_Q_nat, 1 / chi_Q)

    # Combine with the odd-closure identity at mhat0 = 1.
    chi_from_static = sp.solve(sp.expand(static_scaled.lhs - static_scaled.rhs), P0_base)[0] / P0_target
    chi_from_static = sp.simplify(chi_from_static)
    NQ_from_static = sp.simplify(1 / chi_from_static)
    static_scaled_mutated_target = sp.Eq(mhat0_req**2 * lambda_out * P0_base, P0_target + eps)
    chi_from_mutated_target = sp.solve(
        sp.expand(static_scaled_mutated_target.lhs - static_scaled_mutated_target.rhs),
        P0_base,
    )[0] / P0_target
    chi_from_mutated_target = sp.simplify(chi_from_mutated_target)
    NQ_from_mutated_target = sp.simplify(1 / chi_from_mutated_target)
    target_mutation_residual = sp.factor(sp.together(NQ_from_mutated_target - NQ_from_static))

    assert_zero("chi_Q from static-normalized scaled packet", chi_from_static - 1 / (lambda_out * mhat0_req**2))
    assert_zero("natural-source N_Q burden", NQ_from_static - lambda_out * mhat0_req**2)
    assert_zero("natural-source map identity", sp.simplify(natural_source_map.lhs - natural_source_map.rhs).subs(N_Q_nat, 1 / chi_Q))
    assert_nonzero(
        "natural-source N_Q burden detects upstream P0_target mutation",
        target_mutation_residual,
    )

    # Minimal Robin outlet model from PDE stage 93:
    #   chi_Q^R = 3 / (3 - rho)
    chi_robin = sp.simplify(3 / (3 - rho))
    rho_from_NQ = sp.solve(sp.Eq(1 / chi_robin, N_Q_nat), rho)[0]
    assert_zero("Robin outlet inversion", rho_from_NQ - (3 - 3 * N_Q_nat))
    assert_nonzero("Robin outlet inversion detects mutated slope", rho_from_NQ - (3 - (3 + eps) * N_Q_nat))

    step24_exports = export_step24_frontier_points()
    point_step23 = {
        "label": "step23 target-blind ray point",
        **step24_exports["target_blind_scale_0_09"],
    }
    point_step24_q1 = {
        "label": "step24 best point with Q_iso <= 1",
        **step24_exports["best_q_iso_le_1"],
    }
    point_step24_qhalf = {
        "label": "step24 best point with Q_iso <= 0.5",
        **step24_exports["best_q_iso_le_half"],
    }

    def natural_nq(point: dict[str, float]) -> float:
        return point["lambda_out"] * point["mhat0_req"] ** 2

    def natural_chi(point: dict[str, float]) -> float:
        return 1.0 / natural_nq(point)

    def robin_rho(point: dict[str, float]) -> float:
        return 3.0 - 3.0 * natural_nq(point)

    nq_step23 = natural_nq(point_step23)
    nq_q1 = natural_nq(point_step24_q1)
    nq_qhalf = natural_nq(point_step24_qhalf)
    chi_step23 = natural_chi(point_step23)
    chi_q1 = natural_chi(point_step24_q1)
    chi_qhalf = natural_chi(point_step24_qhalf)
    rho_step23 = robin_rho(point_step23)
    rho_q1 = robin_rho(point_step24_q1)
    rho_qhalf = robin_rho(point_step24_qhalf)
    q1_lambda_mutation_residual = natural_nq({**point_step24_q1, "lambda_out": point_step24_q1["lambda_out"] - 1.0}) - nq_q1

    assert_close("same-scale natural N_Q burden is lambda-invariant", nq_step23, nq_q1, 1e-9)
    assert_close("same-scale natural chi is lambda-invariant", chi_step23, chi_q1, 1e-15)
    assert_close("step23 natural N_Q burden", nq_step23, 37872.33995155716, 1e-6)
    assert_close("step24 q<=0.5 natural N_Q burden", nq_qhalf, 45912.04031275389, 1e-6)
    assert_close("step23 natural chi", chi_step23, 2.6404494712422554e-05, 1e-12)
    assert_close("step24 q<=0.5 natural chi", chi_qhalf, 2.1780778923914126e-05, 1e-12)
    assert_close("step23 minimal Robin rho", rho_step23, -113614.01985490319, 1e-6)
    assert_close("step24 q<=0.5 minimal Robin rho", rho_qhalf, -137733.1209385474, 1e-6)
    if abs(q1_lambda_mutation_residual) < 1.0:
        raise AssertionError("lambda_out mutation did not perturb the natural-source burden enough")

    print("STEP 26 NATURAL-SOURCE OUTGOING BURDEN AUDIT")
    print("Tested whether the step-24 outgoing-amplitude frontier improves the exact PDE natural-source outgoing burden.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Exact PDE bridge:")
    print("  chi_Q from the static-normalized scaled packet =", sp.sstr(chi_from_static))
    print("  natural-source N_Q burden =", sp.sstr(NQ_from_static))
    print("  upstream P0_target mutation residual =", sp.sstr(target_mutation_residual))
    print("  minimal Robin rho(N_Q) =", sp.sstr(rho_from_NQ))
    print("  q<=1 lambda mutation burden residual =", q1_lambda_mutation_residual)
    print("Concrete points:")
    for point, nq, chi, rho_value in (
        (point_step23, nq_step23, chi_step23, rho_step23),
        (point_step24_q1, nq_q1, chi_q1, rho_q1),
        (point_step24_qhalf, nq_qhalf, chi_qhalf, rho_qhalf),
    ):
        print(" ", point["label"], "=", point)
        print("    natural-source N_Q burden =", nq)
        print("    natural-source chi_Q =", chi)
        print("    minimal Robin rho =", rho_value)
    print("Interpretation:")
    print("  The target-blind outgoing amplitude coordinate does not reduce the exact natural-source outgoing burden.")
    print("  At fixed frozen scale it only redistributes that burden between lambda_out and mhat0_req.")
    print("  The scale-0.09 step23 point and the scale-0.09 step24 Q_iso<=1 point require the same natural-source N_Q ~ 3.787e4.")
    print("  In the minimal Robin outlet model, those chi_Q values would require rho around -1.136e5 to -1.377e5.")
    print("  So the step-24 improvement is not yet a realized outgoing branch; it still needs a genuinely independent branch-derived transfer/outlet coordinate.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
