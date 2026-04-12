#!/usr/bin/env python3
"""
Step 30 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Uses the Step-28 isotropic DtN deformation law
       chi_Q = 1 + eps*(5 b + a_0/3 + 9 a_5) + O(eps^2)
   and identifies the anomaly expansion parameter eps with the reduced
   electron-anomaly parameter f.
2. Imports the Step-29 outgoing-defect tangent
       Delta_Q(f) = -Lambda_1 f + O(f^2),
   so that the isotropic DtN branch-selection combination is fixed by
       5 b + a_0/3 + 9 a_5 = -Lambda_1.
3. Verifies that the quartic anomaly correction can be written directly as
       Delta(g/2) = -c3_total * (5 b + a_0/3 + 9 a_5) * f^4 + O(f^5).
4. Computes three pure-subbranch realizations (pure argument deformation,
   pure static slot, pure odd-core outlet) and one bookkeeping minimum-norm
   realization of the DtN deformation triple.

Interpretation
--------------
After this step the quartic g-2 closure is no longer an abstract “common PDE
layer.” It is one explicit isotropic DtN branch-selection target. The remaining
PDE task is therefore: determine which actual moving-throat branch produces the
required combination 5 b + a_0/3 + 9 a_5.
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
    banner("STEP 30 — THE QUARTIC g-2 SLIVER AS AN EXPLICIT DtN BRANCH-SELECTION CONSTRAINT")

    f, c3_total, Lambda1 = sp.symbols("f c3_total Lambda_1", real=True)
    b, a0, a5 = sp.symbols("b a_0 a_5", real=True)
    C = sp.simplify(5 * b + a0 / 3 + 9 * a5)

    subbanner("XXX.1 — Exact tangent constraint from the Step-28 DtN deformation law")

    # Step-28 linearized isotropic DtN deformation law with eps identified as f.
    DeltaQ_tangent = sp.expand(f * C)
    print("With eps identified as the anomaly parameter f, the Step-28 law gives")
    print("  Delta_Q(f) =")
    sp.pprint(DeltaQ_tangent)
    print("through first order.")

    target_constraint = sp.Eq(C, -Lambda1)
    print("Matching the Step-29 outgoing-defect tangent requires")
    sp.pprint(target_constraint)

    expect_zero(
        "quartic tangent constraint",
        sp.simplify((C + Lambda1).subs({b: (-Lambda1) / 5, a0: 0, a5: 0})),
    )

    subbanner("XXX.2 — Quartic anomaly law directly in the DtN branch variables")

    Delta_g_over_2 = sp.expand(-c3_total * C * f**4)
    print("Leading quartic anomaly correction =")
    sp.pprint(Delta_g_over_2)

    b_constraint = sp.simplify((-Lambda1 - a0 / 3 - 9 * a5) / 5)
    expect_zero(
        "quartic correction - c3_total*Lambda1*f^4 after imposing the DtN constraint",
        sp.simplify(Delta_g_over_2.subs({b: b_constraint}) - c3_total * Lambda1 * f**4),
    )

    subbanner("XXX.3 — Special isotropic DtN subbranches")

    pure_b = {b: -Lambda1 / 5, a0: sp.Integer(0), a5: sp.Integer(0)}
    pure_a0 = {b: sp.Integer(0), a0: -3 * Lambda1, a5: sp.Integer(0)}
    pure_a5 = {b: sp.Integer(0), a0: sp.Integer(0), a5: -Lambda1 / 9}

    print("Pure argument-deformation branch:")
    sp.pprint(pure_b)
    expect_zero("pure-b constraint", sp.simplify(C.subs(pure_b) + Lambda1))

    print("Pure static-slot branch:")
    sp.pprint(pure_a0)
    expect_zero("pure-a0 constraint", sp.simplify(C.subs(pure_a0) + Lambda1))

    print("Pure odd-core-outlet branch:")
    sp.pprint(pure_a5)
    expect_zero("pure-a5 constraint", sp.simplify(C.subs(pure_a5) + Lambda1))

    subbanner("XXX.4 — Euclidean minimum-norm bookkeeping realization")

    # Solve min ||x||^2 subject to w.x = -Lambda1 with w = (5,1/3,9).
    w = sp.Matrix([5, sp.Rational(1, 3), 9])
    x = -Lambda1 * w / (w.dot(w))
    b_min, a0_min, a5_min = [sp.simplify(v) for v in x]

    print("Minimum-norm bookkeeping realization (in the unweighted Euclidean norm) =")
    sp.pprint(sp.Matrix([b_min, a0_min, a5_min]))
    expect_zero(
        "minimum-norm realization constraint",
        sp.simplify(5 * b_min + a0_min / 3 + 9 * a5_min + Lambda1),
    )

    subbanner("XXX.5 — Numerical values on the carried benchmark")

    Lambda1_num = sp.Float("0.279605891931464")
    f_num = sp.Float("0.001161409732093")
    pure_b_num = {k: sp.N(v.subs(Lambda1, Lambda1_num), 20) for k, v in pure_b.items()}
    pure_a0_num = {k: sp.N(v.subs(Lambda1, Lambda1_num), 20) for k, v in pure_a0.items()}
    pure_a5_num = {k: sp.N(v.subs(Lambda1, Lambda1_num), 20) for k, v in pure_a5.items()}
    x_num = [sp.N(v.subs(Lambda1, Lambda1_num), 20) for v in [b_min, a0_min, a5_min]]

    print("Pure-b tangent branch:")
    sp.pprint(pure_b_num)
    print("Pure-a0 tangent branch:")
    sp.pprint(pure_a0_num)
    print("Pure-a5 tangent branch:")
    sp.pprint(pure_a5_num)
    print("Minimum-norm tangent branch:")
    sp.pprint(sp.Matrix(x_num))

    actual_linear_defect = sp.N((-Lambda1_num) * f_num, 20)
    print("Actual electron-point linearized defect f*(5b + a_0/3 + 9a_5) =")
    sp.pprint(actual_linear_defect)

    banner("STEP 30 LEDGER")
    print("Step-28 linearized isotropic DtN branch law:")
    print("  chi_Q = 1 + f*(5 b + a_0/3 + 9 a_5) + O(f^2)")
    print()
    print("Step-29 outgoing-defect tangent:")
    print("  Delta_Q(f) = -Lambda_1 f + O(f^2)")
    print()
    print("Therefore the exact DtN tangent constraint for the quartic g-2 sliver is")
    print("  5 b + a_0/3 + 9 a_5 = -Lambda_1")
    print()
    print("Equivalent quartic anomaly law:")
    print("  Delta(g/2) = -c3_total*(5 b + a_0/3 + 9 a_5)*f^4 + O(f^5)")
    print()
    print("Three pure realizations:")
    print(f"  pure b   : b = {-sp.N(Lambda1_num/5,16)},   a_0 = 0,                a_5 = 0")
    print(f"  pure a_0 : b = 0,                a_0 = {-sp.N(3*Lambda1_num,16)}, a_5 = 0")
    print(f"  pure a_5 : b = 0,                a_0 = 0,                a_5 = {-sp.N(Lambda1_num/9,16)}")
    print()
    print("Interpretation:")
    print("  The quartic anomaly layer is now an explicit isotropic DtN branch-selection")
    print("  target. The remaining PDE job is to determine which actual moving-throat")
    print("  branch produces the required combination 5 b + a_0/3 + 9 a_5.")


if __name__ == "__main__":
    main()
