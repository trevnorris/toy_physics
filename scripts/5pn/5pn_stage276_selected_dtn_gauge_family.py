#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 276 — exact selected-branch DtN gauge family.

What this script does
---------------------
1. Starts from the exact isotropic deformation formulas

       N_Q   = 3 / (3 S - Sigma_0),
       chi_Q = 3 (S beta^5 + 9 Sigma_5) / (3 S - Sigma_0).

2. Solves exactly for (Sigma_0, Sigma_5) in terms of the observable pair
   (N_Q, chi_Q) and the gauge-like choices (S, beta).
3. Records three useful exact gauges:
      - core gauge      : S = 1, beta = 1
      - scale gauge     : Sigma_0 = 0, beta = 1
      - argument gauge  : Sigma_0 = 0, Sigma_5 = 0 (positive beta branch)
4. Shows explicitly that the selected moving-throat outgoing branch determines a
   two-parameter *family* of DtN deformation triples rather than one unique
   (beta, Sigma_0, Sigma_5) assignment.
"""


def main() -> None:
    banner("STAGE 276 — SELECTED-BRANCH DTN GAUGE FAMILY")

    NQ, chiQ = sp.symbols("N_Q chi_Q", positive=True, real=True)
    S, beta = sp.symbols("S beta", positive=True, real=True)
    Sigma0, Sigma5 = sp.symbols("Sigma_0 Sigma_5", real=True)

    subbanner("I. Exact isotropic deformation relations")
    NQ_expr = sp.simplify(3 / (3 * S - Sigma0))
    chiQ_expr = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))

    print("N_Q =")
    sp.pprint(NQ_expr)
    print("chi_Q =")
    sp.pprint(chiQ_expr)

    subbanner("II. Exact inverse family for Sigma_0 and Sigma_5")
    Sigma0_sol = sp.solve(sp.Eq(NQ_expr, NQ), Sigma0)[0]
    Sigma5_sol = sp.solve(sp.Eq(chiQ_expr.subs(Sigma0, Sigma0_sol), chiQ), Sigma5)[0]

    print("Sigma_0 =")
    sp.pprint(sp.factor(Sigma0_sol))
    print("Sigma_5 =")
    sp.pprint(sp.factor(Sigma5_sol))

    expect_zero("N_Q inverse check", sp.simplify(NQ_expr.subs(Sigma0, Sigma0_sol) - NQ))
    expect_zero(
        "chi_Q inverse check",
        sp.simplify(chiQ_expr.subs({Sigma0: Sigma0_sol, Sigma5: Sigma5_sol}) - chiQ),
    )

    print("So the observable pair (N_Q, chi_Q) fixes only the two combinations")
    print("  Sigma_0 = 3S - 3/N_Q,")
    print("  Sigma_5 = chi_Q/(9N_Q) - S beta^5/9,")
    print("leaving a genuine two-parameter DtN gauge family labeled by (S, beta).")

    subbanner("III. Three exact gauges")

    # Core gauge: S = 1, beta = 1
    Sigma0_core = sp.simplify(Sigma0_sol.subs({S: 1, beta: 1}))
    Sigma5_core = sp.simplify(Sigma5_sol.subs({S: 1, beta: 1}))
    print("Core gauge  (S = 1, beta = 1):")
    print("  Sigma_0 =")
    sp.pprint(sp.factor(Sigma0_core))
    print("  Sigma_5 =")
    sp.pprint(sp.factor(Sigma5_core))
    expect_zero("core gauge N_Q", sp.simplify(NQ_expr.subs({S: 1, beta: 1, Sigma0: Sigma0_core}) - NQ))
    expect_zero(
        "core gauge chi_Q",
        sp.simplify(chiQ_expr.subs({S: 1, beta: 1, Sigma0: Sigma0_core, Sigma5: Sigma5_core}) - chiQ),
    )

    # Scale gauge: Sigma0 = 0, beta = 1
    S_scale = sp.solve(sp.Eq(Sigma0_sol.subs(beta, 1), 0), S)[0]
    Sigma5_scale = sp.simplify(Sigma5_sol.subs({S: S_scale, beta: 1}))
    print("Scale gauge (Sigma_0 = 0, beta = 1):")
    print("  S =")
    sp.pprint(sp.factor(S_scale))
    print("  Sigma_5 =")
    sp.pprint(sp.factor(Sigma5_scale))
    expect_zero("scale gauge N_Q", sp.simplify(NQ_expr.subs({S: S_scale, beta: 1, Sigma0: 0}) - NQ))
    expect_zero(
        "scale gauge chi_Q",
        sp.simplify(chiQ_expr.subs({S: S_scale, beta: 1, Sigma0: 0, Sigma5: Sigma5_scale}) - chiQ),
    )

    # Argument gauge: Sigma0 = 0, Sigma5 = 0, positive beta branch
    S_arg = S_scale
    beta_arg = sp.solve(sp.Eq(Sigma5_sol.subs({S: S_arg, Sigma0: 0, Sigma5: 0}), 0), beta)[0]
    print("Argument gauge (Sigma_0 = 0, Sigma_5 = 0, positive beta branch):")
    print("  S =")
    sp.pprint(sp.factor(S_arg))
    print("  beta =")
    sp.pprint(sp.factor(beta_arg))
    expect_zero("argument gauge N_Q", sp.simplify(NQ_expr.subs({S: S_arg, Sigma0: 0}) - NQ))
    expect_zero(
        "argument gauge chi_Q",
        sp.simplify(chiQ_expr.subs({S: S_arg, beta: beta_arg, Sigma0: 0, Sigma5: 0}) - chiQ),
    )

    subbanner("IV. Useful exact differences between gauges")
    print("Core minus scale gauge static core shift =")
    sp.pprint(sp.simplify(Sigma0_core - 0))
    print("Core minus scale gauge odd core shift =")
    sp.pprint(sp.simplify(Sigma5_core - Sigma5_scale))

    print("\nInterpretation:")
    print("- The actual selected moving-throat outgoing branch does not determine a unique")
    print("  (beta, Sigma_0, Sigma_5) triple by itself.")
    print("- It determines an exact two-parameter gauge family of such triples.")
    print("- The core gauge packages the defect into a static core shift plus an odd core outlet.")
    print("- The scale gauge packages the same defect into a pure mouth normalization plus an odd core outlet.")
    print("- The argument gauge packages it into a pure effective argument deformation on the positive branch.")


if __name__ == "__main__":
    main()
