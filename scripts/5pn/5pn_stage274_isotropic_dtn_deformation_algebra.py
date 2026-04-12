#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 274 — isotropic DtN deformation algebra after canonical outgoing matching.

What this script does
---------------------
1. Starts from the exact outgoing DtN branch

       Lambda_2^out(z) = -3 + z^2/3 + z^4/9 + i z^5/9 + ...

   and considers the most general isotropic deformation of the form

       Lambda_2^def(z)
       = S Lambda_2^out(beta z) + Sigma0 + Sigma2 z^2 + Sigma4 z^4 + i Sigma5 z^5 + O(z^6).

2. Uses the exact relation Yhat = -3/Lambda to compute the deformed normalized response.
3. Imposes the exact condition that the normalized even moments remain canonical, i.e.

       y2/y0 = 1/9,
       y4/y0 = 4/81,

   and solves exactly for Sigma2 and Sigma4.
4. Extracts the exact outgoing-normalization factor

       chi_Q = 27 y5 / y0
             = 3 (S beta^5 + 9 Sigma5) / (3 S - Sigma0).

5. Shows that only beta, Sigma0, and Sigma5 can shift chi_Q once the even moments
   have been held fixed.
"""


def main() -> None:
    banner("STAGE 274 — ISOTROPIC DTN DEFORMATION ALGEBRA")

    z = sp.symbols("z", real=True)
    S, beta = sp.symbols("S beta", positive=True, real=True)
    Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols("Sigma_0 Sigma_2 Sigma_4 Sigma_5", real=True)

    subbanner("I. Deformed isotropic DtN branch")
    Lambda_out = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9
    Lambda_def = sp.expand(S * Lambda_out.subs(z, beta * z) + Sigma0 + Sigma2 * z**2 + Sigma4 * z**4 + sp.I * Sigma5 * z**5)
    Y_def = sp.series(-3 / Lambda_def, z, 0, 6).removeO()

    print("Lambda_2^def(z) =")
    sp.pprint(Lambda_def)
    print("Yhat_2^def(z) = -3/Lambda_2^def(z) =")
    sp.pprint(Y_def)

    y0 = sp.simplify(Y_def.subs(z, 0))
    y2 = sp.simplify(sp.expand(Y_def).coeff(z, 2))
    y4 = sp.simplify(sp.expand(Y_def).coeff(z, 4))
    y5 = sp.simplify(sp.expand(Y_def).coeff(z, 5) / sp.I)

    subbanner("II. Exact even-moment preserving completion")
    sol = sp.solve(
        [sp.Eq(y2 / y0, sp.Rational(1, 9)), sp.Eq(y4 / y0, sp.Rational(4, 81))],
        [Sigma2, Sigma4],
        dict=True,
    )
    if len(sol) != 1:
        raise AssertionError("Expected a unique isotropic even-preserving solve for Sigma2, Sigma4")
    sol = sol[0]

    Sigma2_sol = sp.simplify(sol[Sigma2])
    Sigma4_sol = sp.simplify(sol[Sigma4])

    print("Sigma2 =")
    sp.pprint(Sigma2_sol)
    print("Sigma4 =")
    sp.pprint(Sigma4_sol)

    Y_fix = sp.expand(Y_def.subs({Sigma2: Sigma2_sol, Sigma4: Sigma4_sol}))
    y0_fix = sp.simplify(Y_fix.subs(z, 0))
    y2_fix = sp.simplify(sp.expand(Y_fix).coeff(z, 2))
    y4_fix = sp.simplify(sp.expand(Y_fix).coeff(z, 4))
    y5_fix = sp.simplify(sp.expand(Y_fix).coeff(z, 5) / sp.I)

    expect_zero("(y2/y0) - 1/9", sp.simplify(y2_fix / y0_fix - sp.Rational(1, 9)))
    expect_zero("(y4/y0) - 4/81", sp.simplify(y4_fix / y0_fix - sp.Rational(4, 81)))

    subbanner("III. Exact outgoing-normalization factor")
    NQ = sp.simplify(y0_fix)
    chiQ = sp.simplify(27 * y5_fix / y0_fix)

    print("N_Q = y0 =")
    sp.pprint(NQ)
    print("chi_Q = 27 y5 / y0 =")
    sp.pprint(chiQ)

    expect_zero("N_Q - 3/(3S - Sigma0)", sp.simplify(NQ - 3 / (3 * S - Sigma0)))
    expect_zero(
        "chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0)",
        sp.simplify(chiQ - 3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0)),
    )

    subbanner("IV. Best current reading")
    print("1. Once the normalized even moments are held canonical, Sigma2 and Sigma4 are fixed exactly.")
    print("2. The remaining isotropic outgoing deformation data enter only through")
    print("      beta, Sigma0, Sigma5.")
    print("3. The conservative scalar defect is")
    print("      N_Q = 3 / (3S - Sigma0),")
    print("   while the outgoing normalization factor is")
    print("      chi_Q = 3 (S beta^5 + 9 Sigma5) / (3S - Sigma0).")
    print("4. So after the actual isotropic branch has been reduced to N_Q and chi_Q, the only")
    print("   isotropic DtN data still capable of moving the outgoing theorem are beta, Sigma0, and Sigma5.")


if __name__ == "__main__":
    main()
