#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage258_261_common import banner, subbanner, expect_zero

"""
Stage 262 — actual selected moving-throat grouped-P2 branch map in
(c_pole, epsilon_2, epsilon_4) language.

This script takes the exact Stage-75/79 carrier formula

    c_pole = (1 + epsilon_4) / [4 (1 + epsilon_2)^2]

and pushes it directly into the blocked support-demand quantity zeta_req.
The point is to stop speaking about a generic one-pole carrier and instead
use the actual moving-throat grouped-P2 branch variables themselves.

Main exact outputs:
1. the actual blocked demand map zeta_req(epsilon_2, epsilon_4; eps_blk),
2. the exact isotropic point (epsilon_2,epsilon_4)=(0,0),
3. the exact twin margin carried by the branch,
4. the monotonic sign theorems showing that positive epsilon_2 softens the
   support demand while positive epsilon_4 hardens it.
"""


def main() -> None:
    banner("STAGE 262 — ACTUAL SELECTED GROUPED-P2 BRANCH MAP")

    eps2, eps4, eps_blk = sp.symbols("epsilon_2 epsilon_4 eps_blk", real=True)

    c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2) ** 2))
    rho_alpha = sp.simplify(1 / (1 - c_pole))
    zeta_req = sp.simplify(c_pole / (1 - eps_blk - (1 - 2 * eps_blk) * c_pole))
    blocked_den = sp.simplify(4 * (1 + eps2) ** 2 * (1 - eps_blk) - (1 - 2 * eps_blk) * (1 + eps4))

    subbanner("I. Exact actual grouped-P2 branch map")
    print("c_pole(epsilon_2, epsilon_4) =")
    sp.pprint(c_pole)
    print("rho_alpha(epsilon_2, epsilon_4) =")
    sp.pprint(sp.factor(rho_alpha))
    print("zeta_req(epsilon_2, epsilon_4; eps_blk) =")
    sp.pprint(sp.factor(zeta_req))

    expect_zero(
        "zeta_req - explicit blocked denominator form",
        sp.simplify(zeta_req - (1 + eps4) / blocked_den),
    )

    subbanner("II. Exact isotropic actual branch")
    c_iso = sp.simplify(c_pole.subs({eps2: 0, eps4: 0}))
    rho_iso = sp.simplify(rho_alpha.subs({eps2: 0, eps4: 0}))
    zeta_iso = sp.simplify(zeta_req.subs({eps2: 0, eps4: 0}))
    print("At the actual isotropic point epsilon_2 = epsilon_4 = 0:")
    print("c_pole =")
    sp.pprint(c_iso)
    print("rho_alpha =")
    sp.pprint(rho_iso)
    print("zeta_req =")
    sp.pprint(zeta_iso)
    expect_zero("c_pole(actual)-1/4", sp.simplify(c_iso - sp.Rational(1, 4)))
    expect_zero("rho_alpha(actual)-4/3", sp.simplify(rho_iso - sp.Rational(4, 3)))
    expect_zero("zeta_req(actual)-1/(3-2 eps_blk)", sp.simplify(zeta_iso - 1 / (3 - 2 * eps_blk)))

    subbanner("III. Exact twin-safety numerator carried directly by the actual branch")
    M_twin = sp.expand(2 * (1 + eps2) ** 2 - (1 + eps4))
    print("Twin-safety margin M_twin := 2(1+epsilon_2)^2 - (1+epsilon_4) =")
    sp.pprint(M_twin)
    expect_zero(
        "c_pole - 1/2 + M_twin/[4(1+epsilon_2)^2]",
        sp.simplify(c_pole - sp.Rational(1, 2) + M_twin / (4 * (1 + eps2) ** 2)),
    )
    expect_zero(
        "zeta_req - 1 + 2(1-eps_blk) M_twin/blocked_den",
        sp.simplify(zeta_req - 1 + 2 * (1 - eps_blk) * M_twin / blocked_den),
    )
    print()
    print("So the actual grouped-P2 branch is:")
    print("  twin-safe      iff M_twin >= 0,")
    print("  exactly on the twin boundary iff M_twin = 0,")
    print("  non-twin       iff M_twin < 0.")

    subbanner("IV. Exact contamination monotonicity")
    dc_de2 = sp.simplify(sp.diff(c_pole, eps2))
    dc_de4 = sp.simplify(sp.diff(c_pole, eps4))
    dz_de2 = sp.simplify(sp.diff(zeta_req, eps2))
    dz_de4 = sp.simplify(sp.diff(zeta_req, eps4))
    print("d c_pole / d epsilon_2 =")
    sp.pprint(sp.factor(dc_de2))
    print("d c_pole / d epsilon_4 =")
    sp.pprint(sp.factor(dc_de4))
    print("d zeta_req / d epsilon_2 =")
    sp.pprint(sp.factor(dz_de2))
    print("d zeta_req / d epsilon_4 =")
    sp.pprint(sp.factor(dz_de4))

    expect_zero(
        "dc/de2 + (1+epsilon_4)/(2(1+epsilon_2)^3)",
        sp.simplify(dc_de2 + (1 + eps4) / (2 * (1 + eps2) ** 3)),
    )
    expect_zero(
        "dc/de4 - 1/[4(1+epsilon_2)^2]",
        sp.simplify(dc_de4 - 1 / (4 * (1 + eps2) ** 2)),
    )
    expect_zero(
        "dz/de2 + 8(1-eps_blk)(1+epsilon_2)(1+epsilon_4)/blocked_den^2",
        sp.simplify(dz_de2 + 8 * (1 - eps_blk) * (1 + eps2) * (1 + eps4) / blocked_den ** 2),
    )
    expect_zero(
        "dz/de4 - 4(1-eps_blk)(1+epsilon_2)^2/blocked_den^2",
        sp.simplify(dz_de4 - 4 * (1 - eps_blk) * (1 + eps2) ** 2 / blocked_den ** 2),
    )
    print()
    print("On every admissible blocked branch with 1+epsilon_4 > 0 and 0 <= eps_blk < 1:")
    print("  - positive epsilon_2 lowers c_pole and lowers zeta_req,")
    print("  - positive epsilon_4 raises c_pole and raises zeta_req.")

    banner("STAGE 262 LEDGER")
    print("1. The actual selected moving-throat grouped-P2 branch enters the blocked support-demand")
    print("   problem only through the exact carrier")
    print("      c_pole = (1+epsilon_4) / [4(1+epsilon_2)^2].")
    print("2. At the isotropic actual point epsilon_2 = epsilon_4 = 0 one gets exactly")
    print("      c_pole = 1/4,  rho_alpha = 4/3,  zeta_req = 1/(3-2 eps_blk).")
    print("3. The exact twin-safety numerator is")
    print("      M_twin = 2(1+epsilon_2)^2 - (1+epsilon_4),")
    print("   and the actual branch is twin-safe iff M_twin >= 0.")
    print("4. Positive epsilon_2 contamination softens the support demand, while positive epsilon_4")
    print("   contamination hardens it.  So the geometry lane affects the selected-branch support")
    print("   problem asymmetrically even before any further PDE input is used.")


if __name__ == "__main__":
    main()
