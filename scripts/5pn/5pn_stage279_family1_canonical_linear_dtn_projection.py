#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 279 — linear Family-1 canonical-even projection into the DtN deformation variables.

What this script does
---------------------
1. Starts from the exact compensated hybrid outlet linear defect formulas

       delta E2 = [deltaC - 9 sigma_* delta kappa_W] / [27(1-sigma_*)],
       delta E4 = [5 deltaC - 72 sigma_* delta kappa_W] / [243(1-sigma_*)],
       Delta_Q  = [deltaC - 27 sigma_* delta gamma_W] / [3(1-sigma_*)].

2. Imposes the canonical-even preservation gate deltaE2 = deltaE4 = 0 and proves

       deltaC = 0,
       delta kappa_W = 0.

3. Uses the Stage-141/142 Family-1 transport law to show that this forces

       delta g = 0,

   so all first-order mouth/core broadening transport is frozen out on the canonical-even branch.
4. Translates the remaining outgoing defect into the Stage-92 linear DtN variables
   (b, a_0, a_5). In the natural compensated core gauge one gets

       b = 0,
       a_0 = 0,
       a_5 = - sigma_*/(1-sigma_*) delta gamma_W,

   and therefore

       Delta_Q = 9 a_5.

This is the sharpest first-order bridge from the explicit Family-1 mouth/core system
into the isotropic DtN deformation variables.
"""


def main() -> None:
    banner("STAGE 279 — FAMILY-1 CANONICAL-EVEN DTN PROJECTION")

    sigma, deltaC, dkappa, dgamma = sp.symbols("sigma_star delta_C delta_kappa_W delta_gamma_W", positive=True, real=True)
    dR, dg = sp.symbols("delta_R delta_g", real=True)
    rF1 = sp.symbols("r_F1", positive=True, real=True)
    b, a0, a5 = sp.symbols("b a_0 a_5", real=True)

    subbanner("I. Exact compensated-hybrid linear outlet defects")
    dE2 = sp.simplify((deltaC - 9 * sigma * dkappa) / (27 * (1 - sigma)))
    dE4 = sp.simplify((5 * deltaC - 72 * sigma * dkappa) / (243 * (1 - sigma)))
    DeltaQ = sp.simplify((deltaC - 27 * sigma * dgamma) / (3 * (1 - sigma)))

    print("delta E2 =")
    sp.pprint(dE2)
    print("delta E4 =")
    sp.pprint(dE4)
    print("Delta_Q =")
    sp.pprint(DeltaQ)

    subbanner("II. Exact canonical-even gate")
    # Work with the linear numerators; the physical branch already assumes sigma_* != 0 and sigma_* != 1.
    eq1 = sp.expand((27 * (1 - sigma)) * dE2)
    eq2 = sp.expand((243 * (1 - sigma)) * dE4)
    M = sp.Matrix([[sp.diff(eq1, deltaC), sp.diff(eq1, dkappa)], [sp.diff(eq2, deltaC), sp.diff(eq2, dkappa)]])
    rhs = -sp.Matrix([eq1.subs({deltaC: 0, dkappa: 0}), eq2.subs({deltaC: 0, dkappa: 0})])
    sol_vec = sp.simplify(M.LUsolve(rhs))
    sol_even = {deltaC: sp.simplify(sol_vec[0]), dkappa: sp.simplify(sol_vec[1])}
    print("Canonical-even solution =")
    sp.pprint(sol_even)

    expect_zero("deltaC on canonical-even branch", sol_even[deltaC])
    expect_zero("delta kappa_W on canonical-even branch", sol_even[dkappa])

    DeltaQ_even = sp.simplify(DeltaQ.subs(sol_even))
    print("Remaining outgoing defect on the canonical-even branch =")
    sp.pprint(DeltaQ_even)
    expect_zero(
        "Delta_Q + 9 sigma_* delta gamma_W/(1-sigma_*)",
        sp.simplify(DeltaQ_even + 9 * sigma * dgamma / (1 - sigma)),
    )

    subbanner("III. Family-1 mouth/core transport")
    dR_from_dg = sp.simplify(-dg / sp.sqrt(1 + rF1**2))
    dC_from_dR = sp.simplify(-16 * sigma * dR)
    dC_from_dg = sp.simplify(dC_from_dR.subs(dR, dR_from_dg))

    print("delta R from delta g =")
    sp.pprint(dR_from_dg)
    print("delta C from delta g =")
    sp.pprint(dC_from_dg)

    dg_from_even = sp.solve(sp.Eq(deltaC, dC_from_dg), dg)[0].subs(deltaC, 0)
    print("On the canonical-even branch, delta g =")
    sp.pprint(sp.simplify(dg_from_even))
    expect_zero("delta g on canonical-even branch", sp.simplify(dg_from_even))

    subbanner("IV. Projection into the Stage-92 DtN deformation triple")
    DeltaQ_lin = sp.simplify(5 * b + a0 / 3 + 9 * a5)
    print("Stage-92 linear law =")
    sp.pprint(DeltaQ_lin)

    # Natural compensated core gauge: no argument shift, no static DtN shift once canonical-even is enforced.
    a5_core = sp.simplify(DeltaQ_even / 9)
    print("Natural compensated core-gauge assignment:")
    print("  b = 0")
    print("  a_0 = 0")
    print("  a_5 =")
    sp.pprint(a5_core)

    expect_zero(
        "core-gauge DtN projection check",
        sp.simplify(DeltaQ_lin.subs({b: 0, a0: 0, a5: a5_core}) - DeltaQ_even),
    )

    subbanner("V. Best current reading")
    print("1. Exact first-order preservation of the canonical even l=2 fingerprint forces")
    print("      delta C = 0,   delta kappa_W = 0.")
    print("2. Through the explicit Family-1 transport law, that also forces")
    print("      delta g = 0.")
    print("3. So on the canonical-even compensated branch the only surviving linear mouth/core")
    print("   defect is the odd mixed-channel renormalization delta gamma_W.")
    print("4. In the natural compensated core gauge this is represented by")
    print("      beta = 1,   Sigma_0 = 0,   Sigma_5 = - sigma_*/(1-sigma_*) delta gamma_W,")
    print("   so")
    print("      Delta_Q = 9 Sigma_5.")
    print("5. This makes chi_Q a direct continuum-kernel observable of the selected outgoing branch,")
    print("   not an abstract leftover parameter.")


if __name__ == "__main__":
    main()
