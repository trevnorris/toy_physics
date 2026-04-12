#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage258_261_common import banner, subbanner, expect_zero, family1_refreshed_numbers

"""
Stage 259 — exact geometry-lane regime surface in (epsilon_2, epsilon_4).

This script propagates the exact one-pole support-regime split back into the
selected moving-throat grouped-P2 branch using the exact Stage-75 relation

    c_pole = (1 + epsilon_4) / [4 (1 + epsilon_2)^2].

The result is an exact phase portrait in geometry-contamination language:

- mixed-only boundary     : 1 + epsilon_4 = 0,
- lowest-twin boundary    : 1 + epsilon_4 = 2 (1 + epsilon_2)^2,
- non-twin region         : 1 + epsilon_4 > 2 (1 + epsilon_2)^2,

plus the refreshed Family-1 upper ceiling at eps_blk = 0.  This makes the
remaining selected-branch demand question visible directly in the geometry-lane
variables rather than only in c_pole language.
"""


def main() -> None:
    banner("STAGE 259 — EXACT GEOMETRY-CONTAMINATION REGIME SURFACE")

    eps2, eps4 = sp.symbols("epsilon_2 epsilon_4", real=True)
    c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2) ** 2))
    rho_alpha = sp.simplify(1 / (1 - c_pole))
    twin_indicator = sp.simplify((1 + eps4) - 2 * (1 + eps2) ** 2)

    subbanner("I. Exact one-pole carrier in geometry-lane variables")
    print("c_pole(epsilon_2, epsilon_4) =")
    sp.pprint(c_pole)
    print("rho_alpha(epsilon_2, epsilon_4) =")
    sp.pprint(sp.factor(rho_alpha))
    print("Pi_tr / C_mix = rho_alpha =")
    sp.pprint(sp.factor(rho_alpha))

    subbanner("II. Exact regime surfaces")
    # mixed-only boundary: c_pole = 0
    eps4_mixed = sp.solve(sp.Eq(c_pole, 0), eps4)[0]
    eps4_twin = sp.solve(sp.Eq(c_pole, sp.Rational(1, 2)), eps4)[0]
    print("Mixed-only boundary c_pole = 0  ->  epsilon_4 =")
    sp.pprint(sp.expand(eps4_mixed))
    print("Lowest-twin boundary c_pole = 1/2  ->  epsilon_4 =")
    sp.pprint(sp.expand(eps4_twin))
    expect_zero("mixed-only surface", sp.simplify(eps4_mixed + 1))
    expect_zero("lowest-twin surface", sp.simplify(eps4_twin - (2 * (1 + eps2) ** 2 - 1)))

    print()
    print("Equivalent sign test:")
    print("  c_pole - 1/2 has the same sign as")
    sp.pprint(twin_indicator)
    expect_zero(
        "sign-equivalent numerator",
        sp.simplify(c_pole - sp.Rational(1, 2) - twin_indicator / (4 * (1 + eps2) ** 2)),
    )

    subbanner("III. Refreshed Family-1 admissible strip at eps_blk = 0")
    nums = family1_refreshed_numbers()
    cmax0 = mp.mpf(nums["c_pole_max_0"])
    upper_factor = 4 * cmax0
    non_twin_headroom = upper_factor - 2
    print(f"Family-1 hard unblocked ceiling: c_pole < {cmax0}")
    print(f"Equivalent upper geometry surface: 1 + epsilon_4 < {upper_factor} (1 + epsilon_2)^2")
    print(f"Non-twin strip width factor above the lowest-twin boundary: {non_twin_headroom}")
    print()
    print("So at eps_blk = 0 the exact admissible non-twin strip is")
    print("  2 (1 + epsilon_2)^2 < 1 + epsilon_4 < 4 c_pole,max^(F1) (1 + epsilon_2)^2.")

    subbanner("IV. The actual isotropic grouped-P2 point")
    c_actual = c_pole.subs({eps2: 0, eps4: 0})
    rho_actual = rho_alpha.subs({eps2: 0, eps4: 0})
    margin_twin = sp.simplify((2 * (1 + eps2) ** 2 - (1 + eps4)).subs({eps2: 0, eps4: 0}))
    print(f"At epsilon_2 = epsilon_4 = 0: c_pole = {c_actual}, rho_alpha = {rho_actual}.")
    print(f"Exact lowest-twin margin at the isotropic point: {margin_twin}")
    print("So the actual isotropic grouped-P2 branch lies strictly below the non-twin surface.")

    banner("STAGE 259 LEDGER")
    print("1. The exact selected-branch support regimes can be read directly in geometry-contamination")
    print("   language through c_pole = (1+epsilon_4)/[4(1+epsilon_2)^2].")
    print("2. The universal lowest-twin boundary is")
    print("      1 + epsilon_4 = 2 (1 + epsilon_2)^2.")
    print("   Above it, non-twin asymmetry is required; below it, the symmetric lowest twin remains enough.")
    print("3. On the exact Lambda_EM-refreshed Family-1 branch at eps_blk = 0, the admissible non-twin")
    print("   strip is still finite, while the actual isotropic point (epsilon_2,epsilon_4)=(0,0)")
    print("   stays safely inside the twin-sufficient region.")


if __name__ == "__main__":
    main()
