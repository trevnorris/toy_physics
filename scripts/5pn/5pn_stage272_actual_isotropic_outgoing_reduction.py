#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 272 — actual isotropic outgoing reduction after the geometry-lane verdict.

What this script does
---------------------
1. Encodes the exact actual-isotropic grouped-P2 one-pole conservative branch

       Y_Q^cons(omega) = N_Q [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2) ]

   and shows that every conservative low-frequency coefficient scales by the
   same scalar defect N_Q.
2. Adds the leading retarded normalization factor chi_Q and proves that the odd
   coefficient scales as chi_Q N_Q.
3. Rewrites the full odd normalization condition as

       mhat_0^2 chi_Q N_Q = 1,

   and on the natural source-map branch mhat_0 -> 1 proves that the remaining
   reduced obstruction is purely outgoing:

       N_Q = 1 / chi_Q.
4. Verifies the canonical target relation

       K0_target * sigma_Q^can / 4 = 2 G / (5 c^5).
"""


def main() -> None:
    banner("STAGE 272 — ACTUAL ISOTROPIC OUTGOING REDUCTION")

    omega, OmegaQ = sp.symbols("omega Omega_Q", positive=True, real=True)
    NQ, chiQ, mhat0 = sp.symbols("N_Q chi_Q mhat_0", positive=True, real=True)
    G, c, cs, a = sp.symbols("G c c_s a", positive=True, real=True)
    sigma_can = sp.symbols("sigma_Q_can", positive=True, real=True)

    subbanner("I. Exact actual-isotropic conservative branch")
    Y_cons = sp.simplify(NQ * (sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / OmegaQ**2)))
    Y_cons_series = sp.series(Y_cons, omega, 0, 6).removeO()
    print("Y_Q^cons(omega) =")
    sp.pprint(Y_cons)
    print("Low-frequency series through O(omega^4) =")
    sp.pprint(Y_cons_series)

    K0 = sp.simplify(Y_cons_series.subs(omega, 0))
    K2 = sp.simplify(sp.expand(Y_cons_series).coeff(omega, 2))
    K4 = sp.simplify(sp.expand(Y_cons_series).coeff(omega, 4))

    print("Extracted conservative coefficients:")
    print("  K0 =")
    sp.pprint(K0)
    print("  K2 =")
    sp.pprint(K2)
    print("  K4 =")
    sp.pprint(K4)

    expect_zero("K2/K0 - 1/(4 Omega_Q^2)", sp.simplify(K2 / K0 - 1 / (4 * OmegaQ**2)))
    expect_zero("K4/K0 - 1/(4 Omega_Q^4)", sp.simplify(K4 / K0 - 1 / (4 * OmegaQ**4)))

    print("So once the actual isotropic one-pole conservative carrier is accepted,")
    print("all conservative low-frequency coefficients scale by the single scalar N_Q.")

    subbanner("II. Retarded one-pole completion and the exact odd factorization")
    Y_ret = sp.simplify(
        NQ * (sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / OmegaQ**2 - sp.I * chiQ * sigma_can * omega**5))
    )
    Y_ret_series = sp.series(Y_ret, omega, 0, 6).removeO()
    print("Y_Q^ret(omega) =")
    sp.pprint(Y_ret)
    print("Low-frequency series through O(omega^5) =")
    sp.pprint(Y_ret_series)

    Gamma5 = sp.simplify(sp.expand(Y_ret_series).coeff(omega, 5) / sp.I)
    print("Odd coefficient Gamma5 =")
    sp.pprint(Gamma5)

    expect_zero(
        "Gamma5/(K0 sigma_Q^can/4) - chi_Q",
        sp.simplify(Gamma5 / (K0 * sigma_can / 4) - chiQ),
    )

    print("Therefore the retarded grouped-P2 one-pole branch carries the exact odd scaling")
    print("  Gamma5 / Gamma5_target = chi_Q N_Q")
    print("once the conservative normalization is packaged into N_Q.")

    subbanner("III. Exact target normalization and natural source-map reduction")
    K0_target = sp.simplify(54 * G * cs**5 / (5 * a**5 * c**5))
    sigma_can_target = sp.simplify(4 * a**5 / (27 * cs**5))
    Gamma5_target = sp.simplify(K0_target * sigma_can_target / 4)

    print("K0_target =")
    sp.pprint(K0_target)
    print("sigma_Q^can =")
    sp.pprint(sigma_can_target)
    print("Gamma5_target = K0_target sigma_Q^can / 4 =")
    sp.pprint(Gamma5_target)

    expect_zero(
        "Gamma5_target - 2 G/(5 c^5)",
        sp.simplify(Gamma5_target - 2 * G / (5 * c**5)),
    )

    odd_norm = sp.simplify(mhat0**2 * chiQ * NQ)
    print("Full odd normalization condition =")
    sp.pprint(sp.Eq(odd_norm, 1))

    NQ_natural = sp.simplify(1 / chiQ)
    print("On the natural point-particle source map mhat_0 -> 1, the remaining reduced obstruction is")
    print("  N_Q =")
    sp.pprint(NQ_natural)

    subbanner("IV. Best current reduced reading")
    print("1. The actual isotropic grouped-P2 one-pole conservative branch packages every even moment")
    print("   into the single scalar normalization defect N_Q.")
    print("2. The leading retarded odd coefficient is multiplied by one extra outgoing factor chi_Q.")
    print("3. So the full odd normalization condition is exactly mhat_0^2 chi_Q N_Q = 1.")
    print("4. On the natural source-map branch, the remaining reduced obstruction is purely outgoing:")
    print("      N_Q = 1/chi_Q.")


if __name__ == "__main__":
    main()
