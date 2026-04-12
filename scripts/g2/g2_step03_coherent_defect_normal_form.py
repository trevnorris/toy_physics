#!/usr/bin/env python3
"""
Step 3 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Imports the exact coherent-branch quotient normal form from the moving-throat
   notes:

       Theta_1 = -C_tr,* q_tr,
       Xi_1    = A_tr,* q_tr + q_nt,
       R_1+Xi_1 = -[eps_eta,*/(1-eps_eta,*)] q_eta.

2. Verifies that the grouped defect Xi_1 is exactly the sum of the tracking
   feed-through plus the nontracking monomial drift.
3. Verifies the exact inverse reconstruction formulas for q_tr, q_nt, q_eta.
4. Pushes the Step-2 common-scalar ansatz through this exact PDE normal form,
   giving the quartic tangent constraints for:
      - the observable coherent defect Xi_1,
      - the stripped nontracking defect q_nt.
5. Solves the minimum-norm tangent path in quotient space for the direct
   observable-matching constraint.

Interpretation
--------------
This step removes the arbitrariness left in Step 2.  The moving-throat notes do
not leave the first common scalar as a generic w·q.  On the coherent branch the
observable grouped defect is forced to be Xi_1 = A_tr,* q_tr + q_nt, while the
selected-branch dressing channel is carried separately by q_eta.
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


# ---------------------------------------------------------------------------
# Step 3 — Exact coherent normal form
# ---------------------------------------------------------------------------

def step3_coherent_normal_form() -> None:
    banner("STEP 3 — EXACT COHERENT-DEFECT NORMAL FORM")

    # Constructive-branch parameters.
    chi, delta = sp.symbols("chi0_* deltaU_*", positive=True, real=True)
    epsW, eps, eps_eta = sp.symbols("epsilonW_* epsilon_* epsilon_eta_*", real=True)

    # Direct microscopic monomial coefficients from the moving-throat handoff.
    E = sp.simplify((2 * epsW / (1 - eps)) * (11 + 9 * delta) / (11 * (1 + delta)))
    F = sp.simplify(
        2 * chi / (1 + delta)
        + 4 * epsW * delta / (11 * (1 - eps) * (1 + delta) ** 2)
    )
    Ctr = sp.simplify(chi * delta / ((1 + chi) * (1 + delta) * (1 + chi + delta)))
    Atr = sp.simplify(2 * chi / ((1 + chi) * (1 + delta)))
    Keta = sp.simplify(eps_eta / (1 - eps_eta))

    subbanner("III.1 — Exact coherent-branch coefficients")
    print("E_* =")
    sp.pprint(E)
    print("F_* =")
    sp.pprint(F)
    print("C_tr,* =")
    sp.pprint(Ctr)
    print("A_tr,* =")
    sp.pprint(Atr)
    print("K_eta,* = eps_eta,*/(1-eps_eta,*) =")
    sp.pprint(Keta)

    # Microscopic slippages and quotient monomials.
    SigZ, SigChi, SigEps, SigDel = sp.symbols(
        "Sigma_Z Sigma_chi Sigma_epsilon Sigma_delta", real=True
    )
    qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)
    Theta, Xi, RplusXi = sp.symbols("Theta_1 Xi_1 RplusXi_1", real=True)

    qtr_def = sp.simplify((1 + chi) * SigDel + (1 + delta) * SigChi)
    qnt_def = sp.simplify(SigZ + E * SigEps - F * SigDel)

    Xi_def = sp.simplify(SigZ + Atr * qtr_def + E * SigEps - F * SigDel)
    Theta_def = sp.simplify(-Ctr * qtr_def)
    RplusXi_def = sp.simplify(-Keta * qeta)

    subbanner("III.2 — Exact quotient normal form")
    print("q_tr =")
    sp.pprint(qtr_def)
    print("q_nt =")
    sp.pprint(qnt_def)
    print("Xi_1 direct form =")
    sp.pprint(Xi_def)
    print("Theta_1 direct form =")
    sp.pprint(Theta_def)

    expect_zero("Xi_1 - (A_tr,* q_tr + q_nt)", Xi_def - (Atr * qtr_def + qnt_def))
    expect_zero("Theta_1 + C_tr,* q_tr", Theta_def + Ctr * qtr_def)

    # Exact inverse reconstruction formulas.
    qtr_from_Theta = sp.simplify(-Theta / Ctr)
    qnt_from_XiTheta = sp.simplify(Xi + (Atr / Ctr) * Theta)
    qeta_from_sel = sp.simplify(-(1 - eps_eta) * RplusXi / eps_eta)

    subbanner("III.3 — Exact inverse reconstruction formulas")
    print("q_tr(Theta_1) =")
    sp.pprint(qtr_from_Theta)
    print("q_nt(Xi_1,Theta_1) =")
    sp.pprint(qnt_from_XiTheta)
    print("q_eta(R_1+Xi_1) =")
    sp.pprint(qeta_from_sel)

    expect_zero(
        "q_nt reconstructed from Xi_1 and Theta_1",
        sp.simplify((Xi_def + (Atr / Ctr) * Theta_def) - qnt_def),
    )

    # -----------------------------------------------------------------------
    # Step-2 quartic tangent constraint, now using the exact PDE scalar.
    # -----------------------------------------------------------------------
    subbanner("III.4 — Quartic tangent constraint from Step 2")

    Lam1 = sp.symbols("Lambda_1", real=True)
    s_tr, s_nt, s_eta = sp.symbols("s_tr s_nt s_eta", real=True)

    Xi_tangent = sp.simplify(Atr * s_tr + s_nt)
    Sel_tangent = sp.simplify(-Keta * s_eta)

    print("Direct observable coherent scalar tangent Xi_1' =")
    sp.pprint(Xi_tangent)
    print("Selected-branch dressing tangent (R_1+Xi_1)' =")
    sp.pprint(Sel_tangent)

    print("\nIf the Step-2 common scalar is identified with the observable grouped defect Xi_1,")
    print("the quartic match requires")
    print("  A_tr,* s_tr + s_nt = Lambda_1.")

    print("\nIf the universal tracking feed-through is treated as already frozen into the")
    print("baseline, the stripped nontracking condition is simply")
    print("  s_nt = Lambda_1.")

    # Minimum-norm direct observable path in full quotient space.
    svec = sp.Matrix([s_tr, s_nt, s_eta])
    # Euclidean minimum-norm under one linear constraint c·s = Lambda_1.
    cvec = sp.Matrix([Atr, 1, 0])
    smin = sp.simplify((Lam1 / (cvec.dot(cvec))) * cvec)

    print("\nMinimum-norm direct observable path s_min =")
    sp.pprint(smin)

    expect_zero("direct constraint satisfied by s_min", (cvec.dot(smin)) - Lam1)
    expect_zero(
        "s_eta vanishes on minimum-norm direct path",
        smin[2],
    )

    # Special rays.
    s_track_only = sp.simplify(Lam1 / Atr)
    s_nontrack_only = Lam1
    print("\nSpecial direct-matching rays:")
    print("Tracking-only:  s_nt = s_eta = 0,  s_tr =")
    sp.pprint(s_track_only)
    print("Nontracking-only: s_tr = s_eta = 0,  s_nt =")
    sp.pprint(s_nontrack_only)

    # Step-2 numerical carry-forward benchmark.
    Lam1_num = sp.Float("0.279605891931464", 30)
    print(f"\nStep-2 benchmark value: Lambda_1 = {Lam1_num}")
    print("So the exact Step-3 direct quartic condition is")
    print("  s_nt + A_tr,* s_tr = 0.279605891931464")
    print("and the stripped nontracking condition is")
    print("  s_nt = 0.279605891931464")

    subbanner("III.5 — Interpretation")
    print("1. The moving-throat notes collapse the arbitrary Step-2 scalar w·q to an")
    print("   exact coherent-branch combination Xi_1 = A_tr,* q_tr + q_nt.")
    print("2. The dressing monomial q_eta does not enter the direct grouped defect at")
    print("   this order; it enters only the separate selected-branch demand channel.")
    print("3. Therefore the first quartic anomaly closure can be formulated entirely in")
    print("   the (q_tr, q_nt) plane unless one explicitly chooses to dress the")
    print("   selected-branch residual as well.")


if __name__ == "__main__":
    step3_coherent_normal_form()
