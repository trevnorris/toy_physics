#!/usr/bin/env python3
"""
Stage V2-16 — Branch-freeze / no-refit protocol audit.

Purpose
-------
This script turns the "no-refit" requirement into executable checks:

1. Build a topological ordering for the branch-freeze protocol and verify that
   target/residual information has no path back into branch definitions.
2. Rebuild the grouped-P2 target residuals used by the V2-13/V2-15 bridge.
3. Verify the one-pole, constant-prefactor, 2.5PN normalization, and 4PN-tail
   interface identities.
4. Demonstrate why a freeze certificate is necessary: if post-hoc fitting were
   allowed, a small set of variables has full rank against the target residuals.

The script is intentionally symbolic. It does not decide whether the PDE branch
is physically realized. It certifies the ordering and the algebraic comparison
surface used after a branch has already been frozen.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import sympy as sp


def bool_incidence(expressions: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> sp.Matrix:
    """Return a 0/1 incidence matrix: entry (i,j)=1 iff var_j appears in expr_i."""
    rows = []
    for expr in expressions:
        fs = expr.free_symbols
        rows.append([sp.Integer(1) if v in fs else sp.Integer(0) for v in variables])
    return sp.Matrix(rows)


def is_strictly_upper_triangular(A: sp.Matrix) -> bool:
    """Check that a matrix has no entries on or below the diagonal."""
    n, m = A.shape
    assert n == m
    return all(A[i, j] == 0 for i in range(n) for j in range(i + 1))


def boolean_transitive_closure(A: sp.Matrix) -> sp.Matrix:
    """Boolean transitive closure for a small adjacency matrix."""
    n, m = A.shape
    assert n == m
    reach = [[bool(A[i, j]) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    return sp.Matrix([[sp.Integer(1) if reach[i][j] else sp.Integer(0) for j in range(n)] for i in range(n)])


@dataclass(frozen=True)
class Check:
    name: str
    passed: bool
    detail: str = ""


def main() -> None:
    checks: List[Check] = []

    # ------------------------------------------------------------------
    # 1. Protocol DAG
    # ------------------------------------------------------------------
    nodes = [
        "parent_action",
        "gauge_convention",
        "wall_action_and_geometry",
        "open_boundary_protocol",
        "projection_and_source_map",
        "support_profile_family",
        "branch_solve",
        "coefficient_extraction",
        "target_residual_evaluation",
        "target_decision",
    ]
    idx = {name: i for i, name in enumerate(nodes)}

    allowed_edges = [
        ("parent_action", "gauge_convention"),
        ("parent_action", "wall_action_and_geometry"),
        ("gauge_convention", "support_profile_family"),
        ("wall_action_and_geometry", "open_boundary_protocol"),
        ("open_boundary_protocol", "support_profile_family"),
        ("projection_and_source_map", "support_profile_family"),
        ("support_profile_family", "branch_solve"),
        ("branch_solve", "coefficient_extraction"),
        ("coefficient_extraction", "target_residual_evaluation"),
        ("target_residual_evaluation", "target_decision"),
    ]

    A = sp.zeros(len(nodes))
    for source, target in allowed_edges:
        A[idx[source], idx[target]] = 1

    checks.append(Check(
        "protocol_dag_is_topologically_ordered",
        is_strictly_upper_triangular(A),
        "All allowed arrows point from frozen/earlier choices to later evaluation."
    ))

    reach = boolean_transitive_closure(A)
    forbidden_pairs = [
        ("target_decision", "branch_solve"),
        ("target_decision", "support_profile_family"),
        ("target_residual_evaluation", "branch_solve"),
        ("target_residual_evaluation", "wall_action_and_geometry"),
        ("target_residual_evaluation", "open_boundary_protocol"),
    ]
    forbidden_ok = all(reach[idx[src], idx[dst]] == 0 for src, dst in forbidden_pairs)
    checks.append(Check(
        "no_target_to_branch_feedback_path",
        forbidden_ok,
        "No residual/decision node can reach a branch-definition node."
    ))

    # A deliberately bad edge must be detected.
    A_bad = sp.Matrix(A)
    A_bad[idx["target_residual_evaluation"], idx["support_profile_family"]] = 1
    checks.append(Check(
        "bad_refit_edge_is_detected",
        not is_strictly_upper_triangular(A_bad),
        "Adding residual_evaluation -> support_profile_family breaks the freeze order."
    ))

    # ------------------------------------------------------------------
    # 2. Symbolic target residuals
    # ------------------------------------------------------------------
    K, M = sp.symbols("K M", nonzero=True)
    B0, B2, B4 = sp.symbols("B0 B2 B4")
    Z0, Z2, Z4 = sp.symbols("Z0 Z2 Z4")
    N0, N2, N4 = sp.symbols("N0 N2 N4", nonzero=True)
    mhat0, Sport = sp.symbols("mhat0 Sport", nonzero=True)
    G, c, cs, a_th = sp.symbols("G c c_s a", positive=True)
    Theta_tail, Mtot = sp.symbols("Theta_tail Mtot")

    D0 = K - B0 - Z0
    D2 = -(M + B2 + Z2)
    D4 = -(B4 + Z4)

    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2

    P0 = N0 / D0
    P2 = (D0 * N2 - 2 * D2 * N0) / D0**2
    P4 = (
        D0**2 * N4
        - 2 * D0 * (D2 * N2 + D4 * N0)
        + 3 * D2**2 * N0
    ) / D0**3

    R_one_pole = sp.simplify(D0 * (B4 + Z4) - 3 * (M + B2 + Z2)**2)

    N2_const = sp.simplify(2 * D2 * N0 / D0)
    N4_const = sp.simplify(
        (2 * D0 * (D2 * N2 + D4 * N0) - 3 * D2**2 * N0) / D0**2
    )

    P_eff = mhat0**2 * Sport * P0
    P_target = 54 * G * cs**5 / (5 * a_th**5 * c**5)
    R_norm = sp.simplify(P_eff - P_target)

    gamma_eff = sp.simplify(P_eff * a_th**5 / (27 * cs**5))
    gamma_GR = 2 * G / (5 * c**5)
    gamma_res = sp.simplify(gamma_eff - gamma_GR)

    R_tail = sp.simplify(Theta_tail * (c / cs)**3 - 1)
    C_tail_toy = sp.simplify(Theta_tail * G * Mtot * gamma_eff / (2 * cs**3))
    C_tail_GR = sp.simplify(G * Mtot * gamma_GR / (2 * c**3))

    one_pole_identity = sp.simplify((u4 - 4 * u2**2) * D0**2 - R_one_pole)
    checks.append(Check(
        "one_pole_identity",
        one_pole_identity == 0,
        "(u4 - 4u2^2)D0^2 equals D0(B4+Z4)-3(M+B2+Z2)^2."
    ))

    P2_const_check = sp.simplify(P2.subs(N2, N2_const))
    checks.append(Check(
        "constant_prefactor_P2_condition",
        P2_const_check == 0,
        "P2=0 iff N2=2D2N0/D0."
    ))

    P4_const_check = sp.simplify(P4.subs(N4, N4_const))
    checks.append(Check(
        "constant_prefactor_P4_condition",
        P4_const_check == 0,
        "P4=0 iff the displayed N4 relation holds."
    ))

    gamma_norm_identity = sp.simplify(gamma_res - R_norm * a_th**5 / (27 * cs**5))
    checks.append(Check(
        "normalization_equivalence_to_gamma_GR",
        gamma_norm_identity == 0,
        "R_norm=0 is equivalent to gamma_eff=2G/(5c^5)."
    ))

    tail_sub_check = sp.simplify((C_tail_toy - C_tail_GR).subs({
        P_eff: P_target,  # this direct substitution is not structurally useful
    }))
    # Instead substitute gamma_eff -> gamma_GR and R_tail -> 0 by setting Theta_tail.
    tail_gate_value = sp.solve(sp.Eq(R_tail, 0), Theta_tail)[0]
    tail_gate_check = sp.simplify((C_tail_toy - C_tail_GR).subs({
        gamma_eff: gamma_GR,
        Theta_tail: tail_gate_value,
    }))
    # Because C_tail_toy is already expanded in terms of P_eff, substitute P_eff target explicitly too.
    tail_gate_check = sp.simplify((Theta_tail * G * Mtot * gamma_GR / (2 * cs**3) - C_tail_GR).subs(Theta_tail, tail_gate_value))
    checks.append(Check(
        "tail_transport_gate",
        tail_gate_check == 0,
        "With gamma_eff=gamma_GR and Theta_tail*(c/c_s)^3=1, the toy tail coefficient equals the GR tail coefficient."
    ))

    # ------------------------------------------------------------------
    # 3. Incidence and post-fit rank warning
    # ------------------------------------------------------------------
    residuals = [R_one_pole, P2, P4, R_norm, R_tail]
    residual_names = ["R_one_pole", "R_P2", "R_P4", "R_norm", "R_tail"]
    branch_variables = [K, M, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4, mhat0, Sport, a_th, cs, Theta_tail]
    Inc = bool_incidence(residuals, branch_variables)

    # If post-hoc fitting were allowed, these five variables can generically solve the five residuals.
    post_fit_knobs = [K, N2, N4, N0, Theta_tail]
    J = sp.Matrix(residuals).jacobian(post_fit_knobs)
    # Use the square determinant of the natural submatrix, which is J itself here.
    det_J = sp.factor(sp.together(J.det()))
    expected_det = sp.factor(sp.together((B4 + Z4) * mhat0**2 * Sport * (c / cs)**3 / D0**3))
    checks.append(Check(
        "post_fit_rank_warning_jacobian_nonzero",
        sp.simplify(det_J - expected_det) == 0,
        "The target residuals have a full-rank post-fit sub-Jacobian under generic nonzero/stable conditions."
    ))

    post_fit_rank = J.rank()
    checks.append(Check(
        "post_fit_rank_is_five",
        post_fit_rank == 5,
        "Without branch freeze, five algebraic knobs can generically tune five residuals."
    ))

    # ------------------------------------------------------------------
    # 4. Freeze packet and deterministic certificate hash
    # ------------------------------------------------------------------
    freeze_packet: Dict[str, object] = {
        "stage": "V2-16",
        "protocol": "branch_freeze_no_refit",
        "must_freeze_before_target_evaluation": [
            "parent action and current bookkeeping",
            "gauge-fixing convention",
            "wall/throat action or effective interface action",
            "open-exit impedance boundary protocol",
            "projection/source-map convention",
            "support profile family and number of modes/ports",
            "stability acceptance gates",
            "coefficient extraction formulas",
        ],
        "may_evaluate_only_after_freeze": [
            "one-pole residual",
            "constant-prefactor residuals",
            "universal quadrupole normalization",
            "tail-transport scalar",
            "weak-axisymmetric prefactor slope Xi_1",
        ],
        "forbidden_after_evaluation": [
            "changing support-cardinality",
            "changing boundary condition class",
            "changing gauge convention",
            "changing port/source normalization convention",
            "dropping dark or unstable branches only after target miss unless this was predeclared",
            "adding compensating modes because a target residual is nonzero",
        ],
    }
    packet_json = json.dumps(freeze_packet, sort_keys=True, indent=2)
    packet_hash = hashlib.sha256(packet_json.encode("utf-8")).hexdigest()

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    print("Stage V2-16 — Branch-freeze / no-refit protocol audit")
    print("=" * 72)
    for ch in checks:
        status = "PASS" if ch.passed else "FAIL"
        print(f"{status:4s}  {ch.name}")
        if ch.detail:
            print(f"      {ch.detail}")

    print("\nProtocol nodes in frozen topological order:")
    for i, node in enumerate(nodes):
        print(f"  {i:02d}: {node}")

    print("\nForbidden feedback paths:")
    for src, dst in forbidden_pairs:
        print(f"  {src:30s} -> {dst:30s}: {'BLOCKED' if reach[idx[src], idx[dst]] == 0 else 'OPEN'}")

    print("\nResidual definitions:")
    for name, expr in zip(residual_names, residuals):
        print(f"  {name} = {sp.factor(expr)}")

    print("\nIncidence matrix rows=[R_one_pole,R_P2,R_P4,R_norm,R_tail]")
    print("columns=[K,M,B0,B2,B4,Z0,Z2,Z4,N0,N2,N4,mhat0,Sport,a,c_s,Theta_tail]")
    print(Inc)

    print("\nPost-fit warning:")
    print(f"  post_fit_knobs = {[str(v) for v in post_fit_knobs]}")
    print(f"  det(d residuals / d post_fit_knobs) = {det_J}")
    print("  Generic nonzero conditions: D0!=0, B4+Z4!=0, mhat0*Sport!=0, c/c_s!=0.")
    print("  Therefore a no-refit freeze certificate is mandatory; algebra alone gives too many ways to tune.")

    print("\nFreeze packet SHA256:")
    print(f"  {packet_hash}")

    print("\nFreeze packet JSON:")
    print(packet_json)

    failed = [ch.name for ch in checks if not ch.passed]
    print("\nSummary:")
    print(f"  checks_total = {len(checks)}")
    print(f"  checks_passed = {len(checks)-len(failed)}")
    print(f"  checks_failed = {len(failed)}")
    if failed:
        print(f"  failed = {failed}")
    else:
        print("  all checks passed")


if __name__ == "__main__":
    main()
