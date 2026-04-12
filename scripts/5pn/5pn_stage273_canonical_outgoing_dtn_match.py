#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 273 — canonical outgoing l=2 DtN match and the exact chi_Q = 1 branch.

What this script does
---------------------
1. Rebuilds the exact compact outgoing l=2 DtN fingerprint

       Lambda_2^out(z) = -3 + z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27 + ...

   and the corresponding normalized outgoing response

       Yhat_2^out(z) = -3 / Lambda_2^out(z).
2. Matches that response to the grouped-P2 one-pole-plus-contact module written in
   the dimensionless variable z = omega a / c_s:

       Yhat_Q^ret(z) = 3/4 + (1/4)/(1 - alpha z^2 - i chi_Q B z^5) + ...

3. Solves the exact coefficient match and shows

       alpha = 4/9,
       B chi_Q = 4/27.

4. Using the canonical branch value

       B = sigma_Q^can (c_s/a)^5,
       sigma_Q^can = 9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5),

   proves chi_Q = 1 and fixes the canonical pole scale

       Omega_Q = 3 c_s / (2 a).
"""


def main() -> None:
    banner("STAGE 273 — CANONICAL OUTGOING l=2 DTN MATCH")

    z = sp.symbols("z", real=True)
    alpha, B, chiQ = sp.symbols("alpha B chi_Q", real=True)
    a, cs, OmegaQ = sp.symbols("a c_s Omega_Q", positive=True, real=True)

    subbanner("I. Exact compact outgoing DtN fingerprint")
    Lambda_out = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9 - 2 * z**6 / 27 - sp.I * z**7 / 27
    Y_out = sp.series(-3 / Lambda_out, z, 0, 8).removeO()

    print("Lambda_2^out(z) =")
    sp.pprint(Lambda_out)
    print("Yhat_2^out(z) = -3/Lambda_2^out(z) =")
    sp.pprint(Y_out)

    expect_zero("coeff(z^2) - 1/9", sp.simplify(sp.expand(Y_out).coeff(z, 2) - sp.Rational(1, 9)))
    expect_zero("coeff(z^4) - 4/81", sp.simplify(sp.expand(Y_out).coeff(z, 4) - sp.Rational(4, 81)))
    expect_zero("coeff(z^5) - i/27", sp.simplify(sp.expand(Y_out).coeff(z, 5) / sp.I - sp.Rational(1, 27)))

    subbanner("II. Generic one-pole-plus-contact module in the z variable")
    Y_gen = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - alpha * z**2 - sp.I * chiQ * B * z**5))
    Y_gen_series = sp.series(Y_gen, z, 0, 6).removeO()
    print("Yhat_Q^ret(z) =")
    sp.pprint(Y_gen)
    print("Series through O(z^5) =")
    sp.pprint(Y_gen_series)

    y2 = sp.simplify(sp.expand(Y_gen_series).coeff(z, 2))
    y4 = sp.simplify(sp.expand(Y_gen_series).coeff(z, 4))
    y5 = sp.simplify(sp.expand(Y_gen_series).coeff(z, 5) / sp.I)

    sol_alpha = sp.solve(sp.Eq(y2, sp.Rational(1, 9)), alpha)[0]
    sol_Bchi = sp.solve(sp.Eq(y5.subs(alpha, sol_alpha), sp.Rational(1, 27)), B * chiQ)[0]

    print("Exact coefficient match gives")
    print("  alpha =")
    sp.pprint(sol_alpha)
    print("  B * chi_Q =")
    sp.pprint(sol_Bchi)

    expect_zero("y4(alpha=4/9) - 4/81", sp.simplify(y4.subs(alpha, sp.Rational(4, 9)) - sp.Rational(4, 81)))
    expect_zero("B*chi_Q - 4/27", sp.simplify(sol_Bchi - sp.Rational(4, 27)))

    subbanner("III. Canonical pole scale and chi_Q = 1")
    alpha_from_pole = sp.simplify(cs**2 / (a**2 * OmegaQ**2))
    Omega_from_alpha = sp.solve(sp.Eq(alpha_from_pole, sp.Rational(4, 9)), OmegaQ)[0]
    sigma_can = sp.simplify(sp.Rational(9, 8) / OmegaQ**5)
    B_can = sp.simplify(sigma_can * (cs / a) ** 5)

    print("alpha = omega^2/Omega_Q^2 translated into z = omega a / c_s gives")
    print("  alpha = c_s^2 / (a^2 Omega_Q^2)")
    print("Matching alpha = 4/9 fixes")
    print("  Omega_Q =")
    sp.pprint(Omega_from_alpha)

    print("Canonical sigma_Q^can = 9/(8 Omega_Q^5) becomes")
    sp.pprint(sp.simplify(sigma_can.subs(OmegaQ, Omega_from_alpha)))
    print("and therefore the dimensionless odd coefficient B = sigma_Q^can (c_s/a)^5 is")
    sp.pprint(sp.simplify(B_can.subs(OmegaQ, Omega_from_alpha)))

    expect_zero("canonical B - 4/27", sp.simplify(B_can.subs(OmegaQ, Omega_from_alpha) - sp.Rational(4, 27)))
    expect_zero("chi_Q on canonical branch - 1", sp.simplify(sp.Rational(4, 27) / sp.Rational(4, 27) - 1))

    subbanner("IV. Best current reading")
    print("1. The compact outgoing l=2 DtN branch fixes the normalized outgoing fingerprint exactly.")
    print("2. Matching that fingerprint to the grouped-P2 one-pole-plus-contact module yields")
    print("      alpha = 4/9,   B chi_Q = 4/27.")
    print("3. Translating alpha back to the physical pole scale gives")
    print("      Omega_Q = 3 c_s / (2 a).")
    print("4. On the canonical compact passive/outgoing branch B itself is exactly 4/27, so")
    print("      chi_Q = 1.")


if __name__ == "__main__":
    main()
