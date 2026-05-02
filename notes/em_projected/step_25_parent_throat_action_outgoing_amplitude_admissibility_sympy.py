#!/usr/bin/env python3
"""Admissibility test for the step-24 outgoing amplitude coordinate."""
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


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    if expr == 0:
        raise AssertionError(f"{label} was unexpectedly zero")


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_outgoing_amplitude_admissibility",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "lambda_out identified with outgoing normalization scalar N_Q",
    }

    P0_base, P0_target = sp.symbols("P0_base P0_target", positive=True, real=True)
    mhat0, chi_Q, N_Q = sp.symbols("mhat0 chi_Q N_Q", positive=True, real=True)
    X = sp.symbols("X", positive=True, real=True)

    # Exact laws being combined:
    # 1. static normalization on the scaled outgoing packet
    #       mhat0^2 N_Q P0_base = P0_target
    # 2. exact odd closure
    #       mhat0^2 chi_Q N_Q = 1
    #
    # If the outgoing amplitude coordinate is identified with N_Q, eliminating
    # mhat0^2 N_Q shows that chi_Q is fixed by the unscaled branch P0_base.
    static_scaled = sp.Eq(P0_base * X, P0_target)
    odd_closure = sp.Eq(chi_Q * X, 1)

    X_from_elimination = P0_target / P0_base
    chi_from_elimination = P0_base / P0_target
    NQ_from_elimination = P0_target / (mhat0**2 * P0_base)
    jac_det = sp.Matrix(
        [
            [sp.diff(static_scaled.lhs - static_scaled.rhs, X), sp.diff(static_scaled.lhs - static_scaled.rhs, chi_Q)],
            [sp.diff(odd_closure.lhs - odd_closure.rhs, X), sp.diff(odd_closure.lhs - odd_closure.rhs, chi_Q)],
        ]
    ).det()

    assert_zero("static law solved by X elimination", (static_scaled.lhs - static_scaled.rhs).subs(X, X_from_elimination))
    assert_zero(
        "odd closure solved by chi_Q elimination",
        (odd_closure.lhs - odd_closure.rhs).subs({X: X_from_elimination, chi_Q: chi_from_elimination}),
    )
    assert_zero("chi_Q elimination law", chi_from_elimination - P0_base / P0_target)
    assert_zero("N_Q elimination law", NQ_from_elimination - P0_target / (mhat0**2 * P0_base))
    assert_zero("elimination Jacobian determinant", jac_det.subs(X, X_from_elimination) - P0_target)
    assert_nonzero(
        "mutating chi_Q away from P0_base/P0_target breaks odd closure",
        (odd_closure.lhs - odd_closure.rhs).subs({X: X_from_elimination, chi_Q: 1 / P0_target}),
    )
    assert_nonzero(
        "mutating X away from P0_target/P0_base breaks static law",
        (static_scaled.lhs - static_scaled.rhs).subs(X, P0_target / (2 * P0_base)),
    )

    # Stage-82 reduced finish-line law: every outgoing normalization residual
    # equals N_Q - 1.
    outgoing_finish_line_residual = N_Q - 1

    # Concrete step-24 frontier points.
    point_step23 = {
        "scale": 0.09,
        "lambda_out": 1.0,
        "mhat0_req": 194.6081703105869,
        "Q_iso": 0.4513177752288337,
    }
    point_step24_q1 = {
        "scale": 0.09,
        "lambda_out": 2000.0,
        "mhat0_req": 4.351570977913287,
        "Q_iso": 0.618690285150578,
    }
    point_step24_qhalf = {
        "scale": 0.092,
        "lambda_out": 2000.0,
        "mhat0_req": 4.7912441136331765,
        "Q_iso": 0.4394839373049669,
    }

    def inferred_chi(point: dict[str, float]) -> float:
        return 1.0 / (point["mhat0_req"] ** 2 * point["lambda_out"])

    def inferred_nq_defect(point: dict[str, float]) -> float:
        return float(outgoing_finish_line_residual.subs(N_Q, point["lambda_out"]))

    def odd_closure_residual(point: dict[str, float], chi: float) -> float:
        return point["mhat0_req"] ** 2 * point["lambda_out"] * chi - 1.0

    chi_step23 = inferred_chi(point_step23)
    chi_q1 = inferred_chi(point_step24_q1)
    chi_qhalf = inferred_chi(point_step24_qhalf)
    nq_defect_q1 = inferred_nq_defect(point_step24_q1)
    nq_defect_qhalf = inferred_nq_defect(point_step24_qhalf)
    static_load_ratio_q1 = (
        point_step24_q1["mhat0_req"] ** 2
        * point_step24_q1["lambda_out"]
        / (point_step23["mhat0_req"] ** 2 * point_step23["lambda_out"])
    )
    q1_mutated_lambda_residual = odd_closure_residual({**point_step24_q1, "lambda_out": 1999.0}, chi_q1)

    branch_freeze_payload = {
        "metadata": branch_metadata,
        "frontier_points": [point_step23, point_step24_q1, point_step24_qhalf],
        "symbolic_laws": {
            "chi_Q": sp.sstr(chi_from_elimination),
            "N_Q": sp.sstr(NQ_from_elimination),
            "outgoing_finish_line_residual": sp.sstr(outgoing_finish_line_residual),
        },
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_freeze_payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    assert_close("same-scale chi_Q invariance under lambda_out", chi_step23, chi_q1, 1e-15)
    assert_close("same-scale static load ratio under lambda_out", static_load_ratio_q1, 1.0, 1e-12)
    assert_close("step23-point chi_Q requirement", chi_step23, 2.6404494712422554e-05, 1e-12)
    assert_close("step24 q<=1 chi_Q requirement", chi_q1, 2.6404494712422554e-05, 1e-12)
    assert_close("step24 q<=0.5 chi_Q requirement", chi_qhalf, 2.1780778923914126e-05, 1e-12)
    assert_close("step24 q<=1 N_Q-1", nq_defect_q1, 1999.0, 1e-12)
    assert_close("step24 q<=0.5 N_Q-1", nq_defect_qhalf, 1999.0, 1e-12)
    assert_close("q<=1 odd closure residual", odd_closure_residual(point_step24_q1, chi_q1), 0.0, 1e-14)
    if abs(q1_mutated_lambda_residual) < 1e-5:
        raise AssertionError("lambda_out mutation did not perturb the odd closure enough")

    print("STEP 25 OUTGOING AMPLITUDE ADMISSIBILITY AUDIT")
    print("Tested whether the step-24 outgoing amplitude coordinate is admissible if interpreted as the PDE outgoing-normalization scalar N_Q.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Exact elimination laws:")
    print("  X = mhat0^2 N_Q =", sp.sstr(X_from_elimination))
    print("  chi_Q =", sp.sstr(chi_from_elimination))
    print("  N_Q   =", sp.sstr(NQ_from_elimination))
    print("  elimination Jacobian determinant =", sp.sstr(jac_det.subs(X, X_from_elimination)))
    print("  outgoing finish-line residual =", sp.sstr(outgoing_finish_line_residual))
    print("  Therefore chi_Q is independent of the outgoing amplitude once the static target is imposed.")
    print("Concrete step-24 / step-23 points:")
    print("  target-blind ray point with lambda_out = 1 :", point_step23)
    print("    inferred chi_Q =", chi_step23)
    print("  best sampled point with Q_iso <= 1 :", point_step24_q1)
    print("    inferred chi_Q =", chi_q1)
    print("    outgoing finish-line defect N_Q - 1 =", nq_defect_q1)
    print("    same-scale static load ratio =", static_load_ratio_q1)
    print("    lambda_out mutation residual (1999 instead of 2000) =", q1_mutated_lambda_residual)
    print("  best sampled point with Q_iso <= 0.5 :", point_step24_qhalf)
    print("    inferred chi_Q =", chi_qhalf)
    print("    outgoing finish-line defect N_Q - 1 =", nq_defect_qhalf)
    print("Interpretation:")
    print("  If lambda_out is really the PDE outgoing-normalization scalar N_Q, then the step-24 improvement does not rescue outgoing admissibility.")
    print("  It leaves chi_Q forced to O(10^-5), far from the canonical finish line chi_Q = 1, and the large-amplitude points carry N_Q - 1 = 1999.")
    print("  So step-24 should be read as an upper-bound diagnostic unless a genuinely independent outgoing amplitude coordinate is derived from the branch.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
