#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import (
    banner,
    subbanner,
    expect_zero,
    response_moments_from_D,
    prefactor_moments,
)

"""
Stage 212 — exact isotropic full-bundle target surface.

What this script does
---------------------
1. Rewrites the isotropic grouped-lane bundle in terms of the conservative operator
   moments
      D0 = K - B0 - Z0,
      D2 = -(M + B2 + Z2),
      D4 = -(B4 + Z4).
2. Compiles the normalized grouped-response moments (u2,u4) and outgoing prefactor
   moments (P0,P2,P4).
3. Proves the exact isotropic one-pole defect formula
      u4 - 4 u2^2 = [D0(B4+Z4) - 3(M+B2+Z2)^2] / D0^2.
4. Records the exact constant-prefactor branch conditions P2 = P4 = 0.
5. Packages the universal 2.5PN / 4PN normalization hit into the same isotropic
   target surface.

Interpretation
--------------
This is the sharpest isotropic full-bundle statement before any true moving-throat
branch solve: the completed branch has to land on one exact algebraic surface in
(Bn, Zn, Nn, K, M, mhat0).
"""


if __name__ == "__main__":
    banner("STAGE 212 — EXACT ISOTROPIC FULL-BUNDLE TARGET SURFACE")

    G, c, cs, a = sp.symbols("G c c_s a", positive=True, real=True)
    mhat0 = sp.symbols("mhat_0", positive=True, real=True)
    K, M = sp.symbols("K M", positive=True, real=True)
    B0, B2, B4 = sp.symbols("B0 B2 B4", real=True)
    Z0, Z2, Z4 = sp.symbols("Z0 Z2 Z4", real=True)
    N0, N2, N4 = sp.symbols("N0 N2 N4", real=True)

    D0 = sp.simplify(K - B0 - Z0)
    D2 = sp.simplify(-(M + B2 + Z2))
    D4 = sp.simplify(-(B4 + Z4))

    subbanner("I. Conservative operator moments and normalized response")
    print("D0 =")
    sp.pprint(D0)
    print("D2 =")
    sp.pprint(D2)
    print("D4 =")
    sp.pprint(D4)

    moms = response_moments_from_D(D0, D2, D4)
    print("u2 =")
    sp.pprint(moms["u2"])
    print("u4 =")
    sp.pprint(moms["u4"])

    pole_defect = sp.simplify(moms["u4"] - 4 * moms["u2"] ** 2)
    pole_defect_expected = sp.simplify((D0 * (B4 + Z4) - 3 * (M + B2 + Z2) ** 2) / D0**2)
    print("u4 - 4 u2^2 =")
    sp.pprint(pole_defect)
    expect_zero("one-pole defect formula", pole_defect - pole_defect_expected)

    subbanner("II. Exact isotropic one-pole surface")
    one_pole_surface = sp.simplify(sp.factor(sp.together(sp.numer(pole_defect_expected))))
    print("One-pole surface numerator =")
    sp.pprint(one_pole_surface)
    print("Therefore the exact isotropic one-pole condition is")
    print("  D0 (B4 + Z4) = 3 (M + B2 + Z2)^2")

    subbanner("III. Outgoing-prefactor moments and constant-prefactor branch")
    pref = prefactor_moments(D0, D2, D4, N0, N2, N4)
    print("P0 =")
    sp.pprint(pref["P0"])
    print("P2 =")
    sp.pprint(pref["P2"])
    print("P4 =")
    sp.pprint(pref["P4"])

    N2_const = sp.solve(sp.Eq(pref["P2"], 0), N2)[0]
    N4_const = sp.solve(sp.Eq(pref["P4"].subs(N2, N2_const), 0), N4)[0]
    print("N2 from P2 = 0:")
    sp.pprint(sp.simplify(N2_const))
    print("N4 from P4 = 0 after imposing P2 = 0:")
    sp.pprint(sp.simplify(N4_const))
    expect_zero("P2 on constant-prefactor branch", pref["P2"].subs(N2, N2_const))
    expect_zero("P4 on constant-prefactor branch", pref["P4"].subs({N2: N2_const, N4: N4_const}))

    subbanner("IV. Universal normalization hit")
    P0_target = sp.simplify(54 * G * cs**5 / (5 * a**5 * c**5 * mhat0**2))
    normalization_defect = sp.simplify(mhat0**2 * pref["P0"] - 54 * G * cs**5 / (5 * a**5 * c**5))
    print("mhat0^2 P0 - 54 G c_s^5/(5 a^5 c^5) =")
    sp.pprint(normalization_defect)
    print("Equivalent exact normalization target:")
    print("  mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)")
    print("or")
    print("  P0 =")
    sp.pprint(P0_target)

    banner("STAGE 212 LEDGER")
    print("1. The isotropic one-pole defect is now fully explicit:")
    print("      u4 - 4 u2^2 = [D0(B4+Z4) - 3(M+B2+Z2)^2] / D0^2.")
    print("2. So the exact isotropic one-pole surface is")
    print("      D0(B4+Z4) = 3(M+B2+Z2)^2.")
    print("3. The constant-prefactor branch is exact and gives")
    print("      N2 = 2 D2 N0 / D0,")
    print("      N4 = [2 D0(D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2.")
    print("4. The universal 2.5PN / 4PN hit is the single isotropic ratio")
    print("      mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5).")
    print("5. So the completed isotropic moving-throat branch has to land on one exact")
    print("   combined target surface in the conservative bundle plus outgoing moments.")
