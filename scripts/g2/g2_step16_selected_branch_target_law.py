#!/usr/bin/env python3
"""
Step 16 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Takes the Step-15 universal transfer-shape law
       delta ln T^2 = Lambda_1,
   together with the exact selected-branch identity
       R_target T^2 = Lambda_0 (1 - epsilon_eta).
2. Derives the carried-order demand-ratio law
       delta ln R_target = delta ln(1 - epsilon_eta) - Lambda_1.
3. Rewrites this directly in the quotient dressing coordinate q_eta.
4. Evaluates the direct and coherent tracking-rigid closures from Step 15.
5. Shows that the quartic branch-selection problem is equivalent to deciding
   whether the omitted common layer changes the selected-branch demand ratio or
   leaves it fixed.

Interpretation
--------------
Step 15 proved that the direct transfer-shape update is universal. This step
shows how the same quartic layer appears from the selected-branch side. The two
tracking-rigid laws are now separated by a single yes/no: does R_target drift?
If the physical moving-throat branch is already organized by a fixed selected-
branch spectral target, the coherent closure is the natural one; if not, the
branch may instead absorb the quartic layer by retargeting R_target itself.
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
    banner("STEP 16 — SELECTED-BRANCH DEMAND-RATIO LAW")

    eps_eta_star = sp.symbols("epsilon_eta_star", positive=True, real=True)
    Lambda = sp.symbols("Lambda_1", real=True)
    q_eta = sp.symbols("q_eta", real=True)

    # Stage 167 selected-branch identity at the carried first omitted common order:
    #   R_target * T^2 = Lambda_0 (1 - epsilon_eta),
    # with Lambda_0 inert at this order.
    dln_Tsq = Lambda
    dln_one_minus_eps = sp.simplify(-eps_eta_star / (1 - eps_eta_star) * q_eta)
    dln_Rtarget = sp.simplify(dln_one_minus_eps - dln_Tsq)

    subbanner("XVI.1 — General demand-ratio law")

    print("Using R_target T^2 = Lambda_0 (1 - epsilon_eta) and delta ln T^2 = Lambda_1,")
    print("the carried selected-branch demand-ratio drift is")
    print("  delta ln R_target =")
    sp.pprint(dln_Rtarget)

    print("Equivalently,")
    print("  delta ln R_target = -[epsilon_eta,*/(1-epsilon_eta,*)] q_eta - Lambda_1.")

    # Residual relation from Step 15.
    R1 = sp.simplify(dln_one_minus_eps - Lambda)
    expect_zero("R_1 - delta ln R_target", sp.simplify(R1 - dln_Rtarget))

    print("So the selected-branch residual is literally the demand-ratio drift:")
    print("  R_1 = delta ln R_target.")

    subbanner("XVI.2 — Tracking-rigid direct closure")

    direct_map = {q_eta: 0}
    dln_Rtarget_direct = sp.simplify(dln_Rtarget.subs(direct_map))
    expect_zero("direct closure delta ln R_target + Lambda_1", dln_Rtarget_direct + Lambda)

    print("Direct tracking-rigid closure:")
    print("  delta ln(1-epsilon_eta) = 0")
    print("  delta ln T^2            = Lambda_1")
    print("  delta ln R_target       =")
    sp.pprint(dln_Rtarget_direct)

    subbanner("XVI.3 — Tracking-rigid selected-branch coherent closure")

    q_eta_coh = sp.simplify(-(1 - eps_eta_star) / eps_eta_star * Lambda)
    coherent_map = {q_eta: q_eta_coh}
    dln_Rtarget_coh = sp.simplify(dln_Rtarget.subs(coherent_map))
    expect_zero("coherent closure delta ln R_target", dln_Rtarget_coh)

    print("Tracking-rigid coherent closure:")
    print("  q_eta =")
    sp.pprint(q_eta_coh)
    print("  delta ln(1-epsilon_eta) = Lambda_1")
    print("  delta ln T^2            = Lambda_1")
    print("  delta ln R_target       =")
    sp.pprint(dln_Rtarget_coh)

    subbanner("XVI.4 — Reduced branch-selection theorem")

    print("At the carried quartic order, the two tracking-rigid laws are now separated by")
    print("one scalar criterion:")
    print()
    print("  direct law    <=> delta ln R_target = -Lambda_1")
    print("  coherent law  <=> delta ln R_target = 0")
    print()
    print("So the physical question is no longer whether the transfer shape moves; it always")
    print("does. The real question is whether the selected-branch demand ratio R_target is")
    print("re-targeted by the omitted common layer or remains fixed while the dressing sector")
    print("co-moves to absorb the universal transfer-shape update.")


if __name__ == "__main__":
    main()
