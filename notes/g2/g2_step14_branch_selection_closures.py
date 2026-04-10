#!/usr/bin/env python3
"""
Step 14 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the exact Step-13 quotient laws
       Theta_1 = -C_tr q_tr,
       Xi_1    = A_tr q_tr + q_nt,
       R_1+Xi_1 = -(epsilon_eta,*/(1-epsilon_eta,*)) q_eta.
2. Imposes the carried quartic anomaly target Xi_1 = Lambda_1 and writes the
   full exact matching family in quotient coordinates.
3. Derives three natural closure families:
       (a) tracking-rigid + dressing-rigid,
       (b) tracking-rigid + selected-branch coherent (R_1 = 0),
       (c) dressing-rigid minimum-norm in the quotient plane.
4. Pushes each family through the canonical Step-13 quotient section to obtain
   the minimal microscopic representative in the variables
       (K_eta^(eff), mu_W, T_U).
5. Derives the finite multiplicative laws for the direct microscopic monomials
       C_tr, C_nt, epsilon_eta
   on each closure branch.
6. Computes the exact dressing-rigid minimum-microscopic-section-norm branch,
   which gives the cleanest comparison between a quotient-space optimum and a
   microscopic optimum.

Interpretation
--------------
Step 13 removed the similarity redundancy. Step 14 uses that quotient to turn
"the missing O(f^4) common layer" into explicit finite branch laws.
The main simplification is that, if tracking is kept rigid, the quartic anomaly
is either:

  * a pure nontracking monomial drift, or
  * a locked nontracking+dressing co-drift if the selected branch is also kept
    coherent (R_1 = 0).
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
    banner("STEP 14 — EXACT BRANCH-SELECTION CLOSURES IN THE QUOTIENT")

    # ------------------------------------------------------------------
    # I. Step-13 quotient laws and the general quartic matching family
    # ------------------------------------------------------------------
    subbanner("XIV.1 — Exact quotient laws and the general quartic matching family")

    chi_star, delta_star = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
    F_star = sp.symbols("F_star", real=True)
    eps_eta_star = sp.symbols("epsilon_eta_star", positive=True, real=True)
    Lambda = sp.symbols("Lambda_1", real=True)

    q_tr, q_nt, q_eta = sp.symbols("q_tr q_nt q_eta", real=True)

    C_tr = sp.simplify(
        chi_star * delta_star
        / ((1 + chi_star) * (1 + delta_star) * (1 + chi_star + delta_star))
    )
    A_tr = sp.simplify(2 * chi_star / ((1 + chi_star) * (1 + delta_star)))
    alpha = sp.simplify(1 / (1 + chi_star))
    beta = sp.simplify(F_star / (1 + chi_star))

    Theta = sp.simplify(-C_tr * q_tr)
    Xi = sp.simplify(A_tr * q_tr + q_nt)
    R = sp.simplify(-Xi - eps_eta_star / (1 - eps_eta_star) * q_eta)

    print("Theta_1 =")
    sp.pprint(Theta)
    print("Xi_1 =")
    sp.pprint(Xi)
    print("R_1 =")
    sp.pprint(R)

    q_nt_match = sp.simplify(Lambda - A_tr * q_tr)
    expect_zero("Xi_1 - Lambda on the general matching family", Xi.subs(q_nt, q_nt_match) - Lambda)

    print("General quartic matching family:")
    print("  q_nt = Lambda_1 - A_tr q_tr")
    sp.pprint(q_nt_match)

    # Canonical section from Step 13: only K_eta, mu_W, T_U move.
    section_vec = sp.Matrix([
        -q_eta,
        beta * q_tr + q_nt - q_eta,
        alpha * q_tr,
    ])
    print("Canonical microscopic section coordinates (Delta ln K_eta, Delta ln mu_W, Delta ln T_U)^T =")
    sp.pprint(section_vec)

    # ------------------------------------------------------------------
    # II. Tracking-rigid + dressing-rigid = pure nontracking law
    # ------------------------------------------------------------------
    subbanner("XIV.2 — Tracking-rigid + dressing-rigid closure")

    pure_map = {
        q_tr: 0,
        q_nt: Lambda,
        q_eta: 0,
    }
    pure_section = sp.simplify(section_vec.subs(pure_map))

    expect_zero("pure closure Xi_1 - Lambda", Xi.subs(pure_map) - Lambda)
    expect_zero("pure closure Theta_1", Theta.subs(pure_map))
    expect_zero("pure closure R_1 + Lambda", sp.simplify(R.subs(pure_map) + Lambda))

    print("Tracking-rigid + dressing-rigid quotient point:")
    print("  (q_tr, q_nt, q_eta) = (0, Lambda_1, 0)")
    print("Canonical microscopic representative:")
    sp.pprint(pure_section)

    # ------------------------------------------------------------------
    # III. Tracking-rigid + selected-branch coherence (R_1 = 0)
    # ------------------------------------------------------------------
    subbanner("XIV.3 — Tracking-rigid + selected-branch coherent closure")

    q_eta_sel = sp.simplify(-(1 - eps_eta_star) / eps_eta_star * Lambda)
    sel_map = {
        q_tr: 0,
        q_nt: Lambda,
        q_eta: q_eta_sel,
    }
    sel_section = sp.simplify(section_vec.subs(sel_map))

    expect_zero("selected-branch coherent Xi_1 - Lambda", Xi.subs(sel_map) - Lambda)
    expect_zero("selected-branch coherent Theta_1", Theta.subs(sel_map))
    expect_zero("selected-branch coherent R_1", sp.simplify(R.subs(sel_map)))

    print("Selected-branch coherent quotient point:")
    print("  q_tr = 0")
    print("  q_nt = Lambda_1")
    print("  q_eta =")
    sp.pprint(q_eta_sel)
    print("Canonical microscopic representative:")
    sp.pprint(sel_section)

    # ------------------------------------------------------------------
    # IV. Dressing-rigid minimum norm in the quotient plane
    # ------------------------------------------------------------------
    subbanner("XIV.4 — Dressing-rigid minimum norm in the quotient plane")

    qtr_var = sp.symbols("qtr_var", real=True)
    qnt_line = sp.simplify(Lambda - A_tr * qtr_var)
    Nq = sp.expand(qtr_var**2 + qnt_line**2)
    dNq = sp.diff(Nq, qtr_var)
    qtr_qmin_raw = sp.simplify(sp.solve(sp.Eq(dNq, 0), qtr_var)[0])
    qtr_qmin = sp.simplify(A_tr * Lambda / (1 + A_tr**2))
    qnt_qmin = sp.simplify(Lambda / (1 + A_tr**2))
    expect_zero("closed quotient minimum q_tr", qtr_qmin_raw - qtr_qmin)
    expect_zero("closed quotient minimum q_nt", qnt_line.subs(qtr_var, qtr_qmin) - qnt_qmin)

    qmin_map = {
        q_tr: qtr_qmin,
        q_nt: qnt_qmin,
        q_eta: 0,
    }
    qmin_section = sp.simplify(section_vec.subs(qmin_map))

    expect_zero("quotient-norm minimum Xi_1 - Lambda", Xi.subs(qmin_map) - Lambda)
    expect_zero("quotient-norm stationary condition", dNq.subs(qtr_var, qtr_qmin))

    print("Dressing-rigid minimum-norm quotient point:")
    print("  q_tr =")
    sp.pprint(qtr_qmin)
    print("  q_nt =")
    sp.pprint(qnt_qmin)
    print("Canonical microscopic representative:")
    sp.pprint(qmin_section)

    # ------------------------------------------------------------------
    # V. Dressing-rigid minimum norm in the canonical microscopic section
    # ------------------------------------------------------------------
    subbanner("XIV.5 — Dressing-rigid minimum norm in the canonical microscopic section")

    qnt_line2 = sp.simplify(Lambda - A_tr * qtr_var)
    section_sq = sp.expand((beta * qtr_var + qnt_line2)**2 + (alpha * qtr_var)**2)
    d_section_sq = sp.diff(section_sq, qtr_var)
    qtr_mmin_raw = sp.simplify(sp.solve(sp.Eq(d_section_sq, 0), qtr_var)[0])
    qtr_mmin = sp.simplify((A_tr - beta) * Lambda / ((A_tr - beta)**2 + alpha**2))
    qnt_mmin = sp.simplify(Lambda - A_tr * qtr_mmin)
    expect_zero("closed microscopic minimum q_tr", qtr_mmin_raw - qtr_mmin)
    expect_zero("closed microscopic minimum q_nt", qnt_line2.subs(qtr_var, qtr_mmin) - qnt_mmin)

    mmin_map = {
        q_tr: qtr_mmin,
        q_nt: qnt_mmin,
        q_eta: 0,
    }
    mmin_section = sp.simplify(section_vec.subs(mmin_map))

    expect_zero("microscopic-section minimum Xi_1 - Lambda", Xi.subs(mmin_map) - Lambda)
    expect_zero("microscopic-section stationary condition", d_section_sq.subs(qtr_var, qtr_mmin))

    print("Dressing-rigid minimum-microscopic-section quotient point:")
    print("  q_tr =")
    sp.pprint(qtr_mmin)
    print("  q_nt =")
    sp.pprint(qnt_mmin)
    print("Canonical microscopic representative:")
    sp.pprint(mmin_section)

    # Special condition under which the microscopic minimum collapses back to
    # the pure nontracking law.
    special_F = sp.solve(sp.Eq(beta, A_tr), F_star)[0]
    print("Special condition beta = A_tr  <=>  F_* =")
    sp.pprint(sp.simplify(special_F))
    expect_zero(
        "microscopic minimum collapses to pure nontracking when beta = A_tr",
        sp.simplify(qtr_mmin.subs(F_star, special_F)),
    )

    # ------------------------------------------------------------------
    # VI. Finite monomial laws on the main closure branches
    # ------------------------------------------------------------------
    subbanner("XIV.6 — Finite monomial laws on the main closure branches")

    Ctr, Cnt, eps_eta = sp.symbols("C_tr C_nt epsilon_eta", positive=True, real=True)

    def monomial_update(qtr_expr: sp.Expr, qnt_expr: sp.Expr, qeta_expr: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
        return (
            sp.simplify(Ctr * sp.exp(qtr_expr)),
            sp.simplify(Cnt * sp.exp(qnt_expr)),
            sp.simplify(eps_eta * sp.exp(qeta_expr)),
        )

    pure_mon = monomial_update(0, Lambda, 0)
    sel_mon = monomial_update(0, Lambda, q_eta_sel)
    qmin_mon = monomial_update(qtr_qmin, qnt_qmin, 0)
    mmin_mon = monomial_update(qtr_mmin, qnt_mmin, 0)

    print("Tracking-rigid + dressing-rigid finite monomial law:")
    print("  C_tr  ->")
    sp.pprint(pure_mon[0])
    print("  C_nt  ->")
    sp.pprint(pure_mon[1])
    print("  eps_eta ->")
    sp.pprint(pure_mon[2])

    print("Tracking-rigid + selected-branch coherent finite monomial law:")
    print("  C_tr  ->")
    sp.pprint(sel_mon[0])
    print("  C_nt  ->")
    sp.pprint(sel_mon[1])
    print("  eps_eta ->")
    sp.pprint(sel_mon[2])

    print("Dressing-rigid quotient-minimum finite monomial law:")
    print("  C_tr  ->")
    sp.pprint(qmin_mon[0])
    print("  C_nt  ->")
    sp.pprint(qmin_mon[1])

    print("Dressing-rigid microscopic-minimum finite monomial law:")
    print("  C_tr  ->")
    sp.pprint(mmin_mon[0])
    print("  C_nt  ->")
    sp.pprint(mmin_mon[1])

    # ------------------------------------------------------------------
    # VII. One carried numerical datum from Step 2
    # ------------------------------------------------------------------
    subbanner("XIV.7 — Carried numerical size of the pure nontracking finite update")

    Lambda_num = sp.Float("0.279605891931464")
    amplification = sp.N(sp.exp(Lambda_num), 16)
    print("exp(Lambda_1) with the carried Step-2 value Lambda_1 = 0.279605891931464 is")
    print(amplification)

    banner("STEP 14 COMPLETE")


if __name__ == "__main__":
    main()
