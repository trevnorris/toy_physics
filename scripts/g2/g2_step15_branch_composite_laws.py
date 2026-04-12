#!/usr/bin/env python3
"""
Step 15 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the Step-14 exact quotient matching family
       q_nt = Lambda_1 - A_tr q_tr.
2. Uses the Stage-167 branch-composite dictionary at the carried weak-
   axisymmetric order,
       delta ln R_tr      = -(1/C_*) q_tr,
       delta ln N_*       = q_nt,
       delta ln epsilon_eta = q_eta,
   where N_* = T^2 R_tr^{B_*}.
3. Proves the exact coefficient identity
       A_tr = B_*/C_*,
   which makes the direct transfer-shape drift universal.
4. Derives the order-f^4 branch-composite law
       delta ln T^2 = Lambda_1,
   independent of the choice of q_tr or q_eta.
5. Evaluates the two tracking-rigid Step-14 closures:
       (a) direct nontracking closure,
       (b) selected-branch coherent closure,
   and shows that they differ only in the dressing response.
6. Shows that the geometric minimum-norm closures change only the partition
   between R_tr and N_*, never the universal transfer-shape update.

Interpretation
--------------
After Step 14 the quotient closures were explicit, but still lived in the exact
monomial quotient coordinates. This step translates them back into the actual
moving-throat branch composites. The key result is that, once Xi_1 = Lambda_1 is
imposed, the direct transfer shape T^2 always acquires the same quartic update.
The only remaining physical branch choice is whether the dressing sector co-drifts
with that update or stays rigid.
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
    banner("STEP 15 — BRANCH-COMPOSITE LAWS AND UNIVERSAL TRANSFER-SHAPE DRIFT")

    # ------------------------------------------------------------------
    # I. Step-14 matching family and the Stage-167 composite dictionary
    # ------------------------------------------------------------------
    subbanner("XV.1 — Matching family and branch-composite dictionary")

    chi_star, delta_star = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
    eps_eta_star = sp.symbols("epsilon_eta_star", positive=True, real=True)
    Lambda = sp.symbols("Lambda_1", real=True)
    beta = sp.symbols("beta", real=True)

    q_tr, q_eta = sp.symbols("q_tr q_eta", real=True)
    q_nt = sp.simplify(Lambda - (2 * chi_star / ((1 + chi_star) * (1 + delta_star))) * q_tr)

    A_tr = sp.simplify(2 * chi_star / ((1 + chi_star) * (1 + delta_star)))
    B_star = sp.simplify(2 * (1 + chi_star + delta_star) / delta_star)
    C_star = sp.simplify((1 + chi_star) * (1 + delta_star) * (1 + chi_star + delta_star) / (chi_star * delta_star))

    print("Step-14 exact matching family:")
    print("  q_nt = Lambda_1 - A_tr q_tr")
    sp.pprint(q_nt)

    print("Stage-167 reference-branch composite coefficients:")
    print("  A_tr =")
    sp.pprint(A_tr)
    print("  B_*  =")
    sp.pprint(B_star)
    print("  C_*  =")
    sp.pprint(C_star)

    # Stage-167 branch-composite dictionary at the carried weak-axisymmetric order.
    dln_Rtr = sp.simplify(-q_tr / C_star)
    dln_Nstar = sp.simplify(q_nt)
    dln_eps = sp.simplify(q_eta)

    print("Carried branch-composite drifts:")
    print("  delta ln R_tr =")
    sp.pprint(dln_Rtr)
    print("  delta ln N_* =")
    sp.pprint(dln_Nstar)
    print("  delta ln epsilon_eta =")
    sp.pprint(dln_eps)

    # ------------------------------------------------------------------
    # II. Exact coefficient identity behind the universal transfer-shape law
    # ------------------------------------------------------------------
    subbanner("XV.2 — Exact identity A_tr = B_*/C_*")

    expect_zero("A_tr - B_*/C_*", sp.simplify(A_tr - B_star / C_star))

    # Since N_* = T^2 R_tr^{B_*}, one has at the carried order
    #   delta ln T^2 = delta ln N_* - B_* delta ln R_tr.
    dln_Tsq = sp.simplify(dln_Nstar - B_star * dln_Rtr)

    print("delta ln T^2 = delta ln N_* - B_* delta ln R_tr =")
    sp.pprint(dln_Tsq)
    expect_zero("delta ln T^2 - Lambda_1", dln_Tsq - Lambda)

    print("Result: once Xi_1 = Lambda_1 is imposed, the direct transfer shape obeys")
    print("        delta ln T^2 = Lambda_1 independent of q_tr and q_eta.")

    # ------------------------------------------------------------------
    # III. Selected-branch composite and residual bookkeeping
    # ------------------------------------------------------------------
    subbanner("XV.3 — Selected-branch composite bookkeeping")

    # Stage 167: E := 1 - epsilon_eta, so at the carried order
    #   delta ln E = -(eps_eta,*/(1-eps_eta,*)) q_eta.
    dln_E = sp.simplify(-eps_eta_star / (1 - eps_eta_star) * q_eta)
    R1 = sp.simplify(dln_E - Lambda)

    print("For E := 1 - epsilon_eta, the carried-order selected-branch drift is")
    print("  delta ln E =")
    sp.pprint(dln_E)
    print("so the residual is")
    print("  R_1 = delta ln E - Lambda_1 =")
    sp.pprint(R1)

    # ------------------------------------------------------------------
    # IV. Tracking-rigid closures from Step 14
    # ------------------------------------------------------------------
    subbanner("XV.4 — Tracking-rigid closures in actual branch variables")

    q_eta_coh = sp.simplify(-(1 - eps_eta_star) / eps_eta_star * Lambda)

    direct_map = {q_tr: 0, q_eta: 0}
    coherent_map = {q_tr: 0, q_eta: q_eta_coh}

    direct_vec = sp.Matrix([
        sp.simplify(dln_Rtr.subs(direct_map)),
        sp.simplify(dln_Nstar.subs(direct_map)),
        sp.simplify(dln_Tsq.subs(direct_map)),
        sp.simplify(dln_eps.subs(direct_map)),
        sp.simplify(dln_E.subs(direct_map)),
        sp.simplify(R1.subs(direct_map)),
    ])
    coherent_vec = sp.Matrix([
        sp.simplify(dln_Rtr.subs(coherent_map)),
        sp.simplify(dln_Nstar.subs(coherent_map)),
        sp.simplify(dln_Tsq.subs(coherent_map)),
        sp.simplify(dln_eps.subs(coherent_map)),
        sp.simplify(dln_E.subs(coherent_map)),
        sp.simplify(R1.subs(coherent_map)),
    ])

    print("Columns: (delta ln R_tr, delta ln N_*, delta ln T^2, delta ln epsilon_eta, delta ln E, R_1)^T")
    print("Direct tracking-rigid closure:")
    sp.pprint(direct_vec)
    print("Tracking-rigid + selected-branch coherent closure:")
    sp.pprint(coherent_vec)

    expect_zero("direct closure delta ln T^2 - Lambda", direct_vec[2] - Lambda)
    expect_zero("direct closure delta ln E", direct_vec[4])
    expect_zero("direct closure R_1 + Lambda", direct_vec[5] + Lambda)

    expect_zero("coherent closure delta ln T^2 - Lambda", coherent_vec[2] - Lambda)
    expect_zero("coherent closure delta ln E - Lambda", coherent_vec[4] - Lambda)
    expect_zero("coherent closure R_1", coherent_vec[5])

    print("Conclusion: the two tracking-rigid laws have the SAME direct transfer-shape update")
    print("            delta ln T^2 = Lambda_1. They differ only in the dressing response.")

    # ------------------------------------------------------------------
    # V. General dressing-rigid families: q_tr only repartitions R_tr and N_*
    # ------------------------------------------------------------------
    subbanner("XV.5 — General dressing-rigid closures only repartition R_tr and N_*")

    qtr_var = sp.symbols("qtr_var", real=True)
    qnt_var = sp.simplify(Lambda - A_tr * qtr_var)

    gen_dln_Rtr = sp.simplify((-q_tr / C_star).subs(q_tr, qtr_var))
    gen_dln_Nstar = sp.simplify(qnt_var)
    gen_dln_Tsq = sp.simplify((gen_dln_Nstar - B_star * gen_dln_Rtr))

    print("Generic dressing-rigid branch:")
    print("  delta ln R_tr =")
    sp.pprint(gen_dln_Rtr)
    print("  delta ln N_* =")
    sp.pprint(gen_dln_Nstar)
    print("  delta ln T^2 =")
    sp.pprint(gen_dln_Tsq)
    expect_zero("generic dressing-rigid delta ln T^2 - Lambda", gen_dln_Tsq - Lambda)

    # Step-14 quotient minimum-norm branch in quotient plane.
    qtr_qmin = sp.simplify(A_tr * Lambda / (1 + A_tr**2))
    qnt_qmin = sp.simplify(Lambda / (1 + A_tr**2))

    # Step-14 minimum norm in the canonical microscopic section.
    alpha = sp.simplify(1 / (1 + chi_star))
    qtr_mmin = sp.simplify((A_tr - beta) * Lambda / ((A_tr - beta)**2 + alpha**2))
    qnt_mmin = sp.simplify(Lambda - A_tr * qtr_mmin)

    expect_zero(
        "quotient-minimum branch still has delta ln T^2 = Lambda",
        sp.simplify((qnt_qmin + B_star * qtr_qmin / C_star) - Lambda),
    )
    expect_zero(
        "microscopic-minimum branch still has delta ln T^2 = Lambda",
        sp.simplify((qnt_mmin + B_star * qtr_mmin / C_star) - Lambda),
    )

    print("So any q_tr admixture changes only how the update is split between")
    print("R_tr and the corrected nontracking composite N_*. The direct transfer")
    print("shape T^2 itself always receives the same quartic drift.")

    # ------------------------------------------------------------------
    # VI. Minimal branch-selection statement after the translation
    # ------------------------------------------------------------------
    subbanner("XV.6 — Reduced physical branch-selection question")

    print("At the carried first omitted common order, the quartic anomaly layer reduces to:")
    print("  1. universal direct law  : delta ln T^2 = Lambda_1,")
    print("  2. optional tracking law : delta ln R_tr = -(1/C_*) q_tr,")
    print("  3. optional dressing law : delta ln epsilon_eta = q_eta.")
    print()
    print("Hence the actual moving-throat branch no longer needs to tell us the direct")
    print("transfer-shape update; that part is fixed. It only needs to tell us whether")
    print("the universal transfer-shape drift is accompanied by")
    print("  - no dressing response, or")
    print("  - the locked dressing co-drift that keeps the selected branch coherent.")


if __name__ == "__main__":
    main()
