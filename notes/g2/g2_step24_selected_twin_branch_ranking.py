#!/usr/bin/env python3
"""
Step 24 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Imports the Step-23 selected-branch result for the natural minimal isotropic
   passive/outgoing quadrupole branch,
       Pi_tr / C_mix = 4/3,
   and rewrites it as the exact one-parameter selected twin-support curve
       epsilon_* = 1 - 3 varrho / 2,
       sigma     = 4/(3 varrho) - 2.
2. Proves that this curve lies strictly inside the Step-22 symmetric-lowest-twin
   window for every allowed varrho in (0, 2/3), so mixed-only and non-twin
   support branches are removed from the live anomaly closure.
3. Pushes the Step-21 primitive ranking thresholds onto that selected curve and
   derives the two exact surviving varrho-thresholds:
       q_W = q_Lambda  at varrho_WL = 2(1+beta^2)/(3(2+beta^2)),
       |q_U| = q_Lambda at varrho_UL = 2(1+beta^2)/(3(1+beta+beta^2)).
4. Uses those thresholds to give the complete quartic primitive ranking on the
   selected twin-support branch.

Interpretation
--------------
Before Step 23, the coherent branch still had three support sectors and three
primitive ranking regimes. After Step 23, only one support sector survives:
exactly the symmetric-lowest-twin curve. Step 24 shows that the whole remaining
quartic ambiguity on that selected branch reduces to two exact thresholds in the
single selector varrho.
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
    banner("STEP 24 — EXACT PRIMITIVE RANKING ON THE SELECTED TWIN-SUPPORT BRANCH")

    beta = sp.symbols("beta", positive=True, real=True)
    varrho = sp.symbols("varrho", positive=True, real=True)
    eps = sp.symbols("epsilon_star", positive=True, real=True)
    sigma = sp.symbols("sigma", real=True)

    subbanner("XXIV.1 — Selected twin-support curve from Step 23")

    eps_sel = sp.simplify(1 - sp.Rational(3, 2) * varrho)
    sigma_sel = sp.simplify(2 * eps_sel / (1 - eps_sel))

    print("Selected branch relation from Pi_tr/C_mix = 4/3:")
    print("epsilon_*(varrho) =")
    sp.pprint(eps_sel)
    print("sigma(varrho) =")
    sp.pprint(sigma_sel)

    expect_zero("epsilon_* - (1 - 3 varrho/2)", sp.simplify(eps_sel - (1 - sp.Rational(3, 2) * varrho)))
    expect_zero("sigma - (4/(3 varrho) - 2)", sp.simplify(sigma_sel - (4 / (3 * varrho) - 2)))

    print("Allowed selected branch domain from 0 < epsilon_* < 1:")
    print("  0 < varrho < 2/3")

    subbanner("XXIV.2 — Exact placement inside the Step-22 support windows")

    sigma_mix_floor = sp.simplify(1 / varrho - 2)
    sigma_twin_ceiling = sp.simplify(2 / varrho - 2)

    print("Mixed-only / twin / non-twin sigma windows from Step 22:")
    print("  mixed-only enough:         0 < sigma <= 1/varrho - 2")
    print("  symmetric lowest twin:     1/varrho - 2 < sigma <= 2/varrho - 2")
    print("  non-twin asymmetry needed: sigma > 2/varrho - 2")

    delta_low = sp.simplify(sigma_sel - sigma_mix_floor)
    delta_high = sp.simplify(sigma_twin_ceiling - sigma_sel)

    print("Selected branch position relative to the twin window:")
    print("sigma_sel - (1/varrho - 2) =")
    sp.pprint(delta_low)
    print("(2/varrho - 2) - sigma_sel =")
    sp.pprint(delta_high)
    expect_zero("sigma_sel - mixed-only floor - 1/(3 varrho)", sp.simplify(delta_low - 1 / (3 * varrho)))
    expect_zero("twin ceiling - sigma_sel - 2/(3 varrho)", sp.simplify(delta_high - 2 / (3 * varrho)))

    print("So the selected branch lies strictly inside the symmetric-lowest-twin window for every varrho in (0, 2/3).")

    subbanner("XXIV.3 — Exact surviving ranking thresholds on the selected branch")

    eps_WL = sp.simplify(1 / (2 + beta**2))
    eps_UL = sp.simplify(beta / (1 + beta + beta**2))

    varrho_WL = sp.simplify(2 * (1 - eps_WL) / 3)
    varrho_UL = sp.simplify(2 * (1 - eps_UL) / 3)

    print("q_W = q_Lambda threshold on the selected branch:")
    sp.pprint(varrho_WL)
    print("|q_U| = q_Lambda threshold on the selected branch:")
    sp.pprint(varrho_UL)

    expect_zero(
        "varrho_WL - 2(1+beta^2)/(3(2+beta^2))",
        sp.simplify(varrho_WL - 2 * (1 + beta**2) / (3 * (2 + beta**2))),
    )
    expect_zero(
        "varrho_UL - 2(1+beta^2)/(3(1+beta+beta^2))",
        sp.simplify(varrho_UL - 2 * (1 + beta**2) / (3 * (1 + beta + beta**2))),
    )

    print("Ordering of the two thresholds:")
    diff_thresholds = sp.simplify(varrho_UL - varrho_WL)
    ceiling_gap = sp.simplify(sp.Rational(2, 3) - varrho_UL)
    print("varrho_UL - varrho_WL =")
    sp.pprint(diff_thresholds)
    print("2/3 - varrho_UL =")
    sp.pprint(ceiling_gap)
    expect_zero(
        "varrho_UL - varrho_WL - 2(1+beta^2)(1-beta)/(3(1+beta+beta^2)(2+beta^2))",
        sp.simplify(diff_thresholds - 2 * (1 + beta**2) * (1 - beta) / (3 * (1 + beta + beta**2) * (2 + beta**2))),
    )
    expect_zero(
        "2/3 - varrho_UL - 2 beta/(3(1+beta+beta^2))",
        sp.simplify(ceiling_gap - 2 * beta / (3 * (1 + beta + beta**2))),
    )

    print("Hence, for 0 < beta < 1,")
    print("  0 < varrho_WL < varrho_UL < 2/3.")

    subbanner("XXIV.4 — Exact primitive ranking sectors on the selected branch")

    print("The Step-21 branch-independent orderings still hold:")
    print("  q_chi > q_Z,")
    print("  q_Z   > q_W,")
    print("  q_W   > |q_U|.")
    print("Only q_Lambda moves relative to q_W and |q_U|.")
    print()
    print("Region I: 0 < varrho < varrho_WL")
    print("  -> epsilon_* > 1/(2+beta^2)")
    print("  -> q_chi > q_Z > q_Lambda > q_W > |q_U|.")
    print()
    print("Region II: varrho_WL < varrho < varrho_UL")
    print("  -> beta/(1+beta+beta^2) < epsilon_* < 1/(2+beta^2)")
    print("  -> q_chi > q_Z > q_W > q_Lambda > |q_U|.")
    print()
    print("Region III: varrho_UL < varrho < 2/3")
    print("  -> 0 < epsilon_* < beta/(1+beta+beta^2)")
    print("  -> q_chi > q_Z > q_W > |q_U| > q_Lambda.")

    subbanner("XXIV.5 — Exact numerical windows on the constructive coherent branch")

    beta_max = sp.Rational(2, 11)
    varrho_WL_min = sp.simplify(sp.limit(varrho_WL, beta, 0, dir="+"))
    varrho_WL_max = sp.simplify(varrho_WL.subs(beta, beta_max))
    varrho_UL_min = sp.simplify(varrho_UL.subs(beta, beta_max))
    varrho_UL_max = sp.simplify(sp.limit(varrho_UL, beta, 0, dir="+"))

    print("Using the coherent constructive bound 0 < beta < 2/11:")
    print("  varrho_WL in (")
    sp.pprint(varrho_WL_min)
    print(",")
    sp.pprint(varrho_WL_max)
    print(")")
    print("  varrho_UL in (")
    sp.pprint(varrho_UL_min)
    print(",")
    sp.pprint(varrho_UL_max)
    print(")")

    samples = [
        (sp.Rational(1, 10), sp.Rational(1, 4)),
        (sp.Rational(1, 10), sp.Rational(2, 5)),
        (sp.Rational(1, 10), sp.Rational(3, 5)),
    ]
    for beta_val, varrho_val in samples:
        eps_val = sp.N(eps_sel.subs(varrho, varrho_val))
        varrho_WL_val = sp.N(varrho_WL.subs(beta, beta_val))
        varrho_UL_val = sp.N(varrho_UL.subs(beta, beta_val))
        if varrho_val < varrho_WL.subs(beta, beta_val):
            verdict = "q_Lambda > q_W > |q_U|"
        elif varrho_val < varrho_UL.subs(beta, beta_val):
            verdict = "q_W > q_Lambda > |q_U|"
        else:
            verdict = "q_W > |q_U| > q_Lambda"
        print(f"beta = {beta_val}, varrho = {varrho_val}: epsilon_* = {eps_val}")
        print(f"  varrho_WL = {varrho_WL_val}, varrho_UL = {varrho_UL_val}")
        print(f"  verdict   = {verdict}\n")

    banner("STEP 24 LEDGER")
    print("Selected support branch from Step 23:")
    print("  Pi_tr/C_mix = 4/3")
    print("  epsilon_* = 1 - 3 varrho/2")
    print("  sigma     = 4/(3 varrho) - 2")
    print("  0 < varrho < 2/3")
    print()
    print("This curve lies strictly inside the Step-22 symmetric-lowest-twin window:")
    print("  1/varrho - 2 < sigma < 2/varrho - 2")
    print()
    print("Exact surviving thresholds on the selected branch:")
    print("  q_W = q_Lambda  at varrho_WL = 2(1+beta^2)/(3(2+beta^2))")
    print("  |q_U| = q_Lambda at varrho_UL = 2(1+beta^2)/(3(1+beta+beta^2))")
    print("  with 0 < varrho_WL < varrho_UL < 2/3")
    print()
    print("Exact rankings on the selected twin-support branch:")
    print("  Region I:   0 < varrho < varrho_WL      -> q_chi > q_Z > q_Lambda > q_W > |q_U|")
    print("  Region II:  varrho_WL < varrho < varrho_UL -> q_chi > q_Z > q_W > q_Lambda > |q_U|")
    print("  Region III: varrho_UL < varrho < 2/3    -> q_chi > q_Z > q_W > |q_U| > q_Lambda")
    print()
    print("Interpretation:")
    print("  The selected minimal isotropic branch collapses the old support phase diagram")
    print("  to one exact twin-support curve with only two remaining thresholds. The")
    print("  outgoing-scale drift q_Lambda beats q_W only in the low-varrho/high-blocking")
    print("  corner; across most of the selected curve the quartic layer is still carried")
    print("  primarily by q_W, while q_U only beats q_Lambda near the very weak-blocking end.")


if __name__ == "__main__":
    main()
