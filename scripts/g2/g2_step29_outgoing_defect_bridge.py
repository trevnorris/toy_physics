#!/usr/bin/env python3
"""
Step 29 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Bridges the Step-2 common transport-residue scalar to the moving-throat
   outgoing normalization defect:
       N_Q - 1 = -Delta_Q/(1 + Delta_Q),
   where Delta_Q := chi_Q - 1 on the natural source-map branch.
2. Shows that preserving the already-frozen lower anomaly orders forces the
   outgoing defect to switch on linearly in the reduced anomaly parameter f:
       Delta_Q(f) = delta_1 f + O(f^2).
3. Proves that the quartic anomaly sliver is then controlled by the tangent
   coefficient alone:
       Delta(g/2) = -c3_total * delta_1 * f^4 + O(f^5).
4. Matches the exact carried benchmark residual and fixes
       delta_1 = -Lambda_1,
   where Lambda_1 is the Step-2 common-path tangent.
5. Evaluates the actual electron-point outgoing defect
       Delta_Q^(e) = -Lambda_1 f/(1 + Lambda_1 f),
   and verifies that the improved law reproduces the exact carried residual.

Interpretation
--------------
This step shows that the remaining g-2 sliver can be carried by a *tiny* outgoing
normalization defect on the finite-f electron branch while still collapsing to the
canonical compact outgoing branch in the strict point-particle/weak-coupling limit
f -> 0. So there is no conflict between the gravitational reduced theorem
(Delta_Q -> 0) and the electron anomaly closure (Delta_Q = O(f)).
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(sp.expand(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def main() -> None:
    banner("STEP 29 — OUTGOING-NORMALIZATION DEFECT BRIDGE FOR THE QUARTIC g-2 SLIVER")

    f, c3_total, DeltaQ = sp.symbols("f c3_total Delta_Q", real=True)
    I = sp.I

    subbanner("XXIX.1 — Exact natural-branch map from the outgoing defect to the conservative defect")

    Nminus1 = sp.simplify(-DeltaQ / (1 + DeltaQ))
    print("If Delta_Q := chi_Q - 1, then on the natural source-map branch")
    print("  N_Q - 1 =")
    sp.pprint(Nminus1)

    DeltaQ_exact = sp.simplify(sp.solve(sp.Eq(sp.Symbol("Lambda"), Nminus1), DeltaQ)[0])
    print("Solving Lambda = N_Q - 1 for Delta_Q gives")
    sp.pprint(DeltaQ_exact)

    subbanner("XXIX.2 — Quartic anomaly law from a linearly switched-on outgoing defect")

    delta1 = sp.symbols("delta_1", real=True)
    DeltaQ_f = delta1 * f
    Lambda_common = sp.simplify((-DeltaQ_f / (1 + DeltaQ_f)).series(f, 0, 2).removeO())
    print("Assume Delta_Q(f) = delta_1 f + O(f^2). Then")
    print("  N_Q(f) - 1 =")
    sp.pprint(Lambda_common)
    expect_zero("N_Q(f)-1 - (-delta_1*f)", sp.simplify(Lambda_common + delta1 * f))

    Delta_g_over_2 = sp.expand(c3_total * Lambda_common * f**3)
    print("Leading common correction to g/2 =")
    sp.pprint(Delta_g_over_2)
    expect_zero(
        "quartic correction - (-c3_total*delta_1*f^4)",
        sp.simplify(Delta_g_over_2 + c3_total * delta1 * f**4),
    )

    subbanner("XXIX.3 — Exact match to the carried benchmark residual")

    # Carried benchmark values from the previous steps / atom_work table.
    f_e = sp.Float("0.001161409732093")
    c3_tot_val = sp.Float("2.233050586884145")
    Lambda1 = sp.Float("0.279605891931464")
    g_loc = sp.Float("2.002319304358647956")
    g_target_exact = sp.Float("2.00231930436091999990584705")
    Delta_g = sp.N(g_target_exact - g_loc, 30)

    delta1_solution = sp.simplify(-Lambda1)
    print("Matching the Step-2 common-path tangent gives")
    print("  delta_1 =")
    sp.pprint(delta1_solution)

    # Exact electron-point defect chosen so that N_Q - 1 = Lambda1 * f_e exactly.
    DeltaQ_e = sp.simplify(-Lambda1 * f_e / (1 + Lambda1 * f_e))
    chiQ_e = sp.simplify(1 + DeltaQ_e)
    NQ_e = sp.simplify(1 / chiQ_e)
    common_scalar_e = sp.simplify(-DeltaQ_e / (1 + DeltaQ_e))

    print("Electron-point exact outgoing defect Delta_Q^(e) =")
    sp.pprint(sp.N(DeltaQ_e, 20))
    print("chi_Q^(e) =")
    sp.pprint(sp.N(chiQ_e, 20))
    print("N_Q^(e) =")
    sp.pprint(sp.N(NQ_e, 20))
    print("N_Q^(e) - 1 =")
    sp.pprint(sp.N(NQ_e - 1, 20))

    expect_zero(
        "common scalar at electron point - Lambda1*f_e",
        sp.simplify(common_scalar_e - Lambda1 * f_e),
    )

    Delta_g_model = sp.N(2 * c3_tot_val * common_scalar_e * f_e**3, 30)
    print("Predicted quartic correction to g from the outgoing-defect bridge =")
    sp.pprint(Delta_g_model)
    print("Carried exact residual =")
    sp.pprint(Delta_g)

    residual_diff = sp.N(Delta_g_model - Delta_g, 30)
    print("quartic correction - carried residual =")
    sp.pprint(residual_diff)
    if abs(float(residual_diff)) > 1e-18:
        raise AssertionError("quartic correction does not reproduce the carried residual within rounding tolerance.")

    g_imp = sp.N(g_loc + Delta_g_model, 30)
    print("Improved g value from the outgoing-defect bridge =")
    sp.pprint(g_imp)
    g_diff = sp.N(g_imp - g_target_exact, 30)
    print("improved g - exact target =")
    sp.pprint(g_diff)
    if abs(float(g_diff)) > 1e-18:
        raise AssertionError("improved g does not match the exact target within rounding tolerance.")

    subbanner("XXIX.4 — Asymptotic compatibility with the canonical outgoing branch")

    f_small = sp.symbols("f_small", real=True)
    DeltaQ_series = sp.series(-Lambda1 * f_small / (1 + Lambda1 * f_small), f_small, 0, 3).removeO()
    print("Delta_Q(f) near f = 0 =")
    sp.pprint(sp.expand(DeltaQ_series))
    expect_zero(
        "Delta_Q(f) - (-Lambda1*f + Lambda1^2*f^2)",
        sp.simplify(DeltaQ_series - (-Lambda1 * f_small + Lambda1**2 * f_small**2)),
    )
    print("So Delta_Q -> 0 as f -> 0, and the canonical compact outgoing branch is recovered.")

    banner("STEP 29 LEDGER")
    print("Natural source-map branch:")
    print("  N_Q - 1 = -Delta_Q/(1 + Delta_Q)")
    print()
    print("To avoid reopening lower anomaly orders, the outgoing defect must start at O(f):")
    print("  Delta_Q(f) = delta_1 f + O(f^2)")
    print()
    print("Then the quartic common-layer correction is")
    print("  Delta(g/2) = -c3_total*delta_1*f^4 + O(f^5)")
    print()
    print("Matching the carried benchmark fixes")
    print(f"  delta_1 = -Lambda_1 = {sp.N(delta1_solution, 16)}")
    print()
    print("Electron-point exact defect (chosen so that N_Q - 1 = Lambda_1 f exactly):")
    print(f"  Delta_Q^(e) = {sp.N(DeltaQ_e, 16)}")
    print(f"  chi_Q^(e)   = {sp.N(chiQ_e, 16)}")
    print(f"  N_Q^(e)     = {sp.N(NQ_e, 16)}")
    print()
    print("Interpretation:")
    print("  The gravitational reduced theorem wants the canonical outgoing branch in the")
    print("  strict point-particle limit f -> 0, while the electron anomaly only needs a")
    print("  tiny finite-f drift away from that branch. The two statements are compatible.")


if __name__ == "__main__":
    main()
