#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage258_261_common import banner, subbanner, expect_zero, family1_refreshed_numbers

"""
Stage 263 — exact blocked Family-1 corridor test written directly in the actual
selected moving-throat grouped-P2 branch variables.

This script takes the exact blocked Family-1 ceiling from Stage 261 and rewrites
it as an inequality in (epsilon_2, epsilon_4) using the actual grouped-P2 branch
carrier.  It then evaluates the exact margins at the actual isotropic point and
shows how blocking changes only the upper Family-1 margin, not the universal
lowest-twin margin.
"""


def main() -> None:
    banner("STAGE 263 — ACTUAL GROUPED-P2 BLOCKED CORRIDOR TEST")

    eps2, eps4, eps_blk, zmax = sp.symbols("epsilon_2 epsilon_4 eps_blk zeta_max", positive=True, real=True)

    c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2) ** 2))
    cmax = sp.simplify(zmax * (1 - eps_blk) / (1 + (1 - 2 * eps_blk) * zmax))

    M_twin = sp.expand(2 * (1 + eps2) ** 2 - (1 + eps4))
    M_F1 = sp.expand(4 * cmax * (1 + eps2) ** 2 - (1 + eps4))

    subbanner("I. Exact blocked corridor in the actual branch variables")
    print("Lowest-twin boundary numerator M_twin =")
    sp.pprint(M_twin)
    print("Family-1 admissibility numerator M_F1 =")
    sp.pprint(sp.factor(M_F1))

    expect_zero(
        "c_pole - 1/2 + M_twin/[4(1+epsilon_2)^2]",
        sp.simplify(c_pole - sp.Rational(1, 2) + M_twin / (4 * (1 + eps2) ** 2)),
    )
    expect_zero(
        "c_pole - cmax + M_F1/[4(1+epsilon_2)^2]",
        sp.simplify(c_pole - cmax + M_F1 / (4 * (1 + eps2) ** 2)),
    )

    print()
    print("So the exact blocked corridor is")
    print("  M_twin >= 0   -> symmetric lowest twin is enough,")
    print("  M_twin < 0    -> non-twin asymmetry is required,")
    print("  M_F1 > 0      -> the exact Lambda_EM-refreshed Family-1 side is still admissible.")

    subbanner("II. Exact actual isotropic margins")
    M_twin_iso = sp.simplify(M_twin.subs({eps2: 0, eps4: 0}))
    M_F1_iso = sp.simplify(M_F1.subs({eps2: 0, eps4: 0}))
    print("M_twin at the actual isotropic point =")
    sp.pprint(M_twin_iso)
    print("M_F1 at the actual isotropic point =")
    sp.pprint(sp.factor(M_F1_iso))
    expect_zero("M_twin(actual)-1", sp.simplify(M_twin_iso - 1))
    expect_zero("M_F1(actual) - (4 cmax - 1)", sp.simplify(M_F1_iso - (4 * cmax - 1)))

    dMF_de = sp.simplify(sp.diff(M_F1_iso, eps_blk))
    print("d M_F1(actual) / d eps_blk =")
    sp.pprint(sp.factor(dMF_de))
    expect_zero(
        "dM_F1(actual)/de - 4 zmax(zmax-1)/(1+(1-2 eps_blk)zmax)^2",
        sp.simplify(dMF_de - 4 * zmax * (zmax - 1) / (1 + (1 - 2 * eps_blk) * zmax) ** 2),
    )
    print()
    print("So the exact isotropic branch sits one full unit inside the twin-safe strip, while")
    print("its Family-1 admissibility margin is 4 cmax(eps_blk)-1 and grows monotonically with blocking")
    print("whenever zeta_max^(F1) > 1.")

    subbanner("III. Refreshed numerical margins on the exact Lambda_EM branch")
    nums = family1_refreshed_numbers()
    zmax_num = mp.mpf(nums["zeta_max"])
    cmax0_num = mp.mpf(nums["c_pole_max_0"])
    epscrit_num = mp.mpf(nums["eps_blk_crit"])

    def cmax_num(eb: mp.mpf) -> mp.mpf:
        return zmax_num * (1 - eb) / (1 + (1 - 2 * eb) * zmax_num)

    margins = {
        "unblocked": 4 * cmax0_num - 1,
        "critical-blocking": 4 * cmax_num(epscrit_num) - 1,
    }
    print(f"zeta_max^(F1) ≈ {zmax_num}")
    print(f"c_pole,max^(F1)(0) ≈ {cmax0_num}")
    print(f"eps_blk^crit ≈ {epscrit_num}")
    print(f"Exact twin margin at the actual isotropic point: 1")
    print(f"Family-1 admissibility margin at eps_blk = 0: {margins['unblocked']}")
    print(f"Family-1 admissibility margin at eps_blk = eps_crit: {margins['critical-blocking']}")
    print()
    print("So the actual isotropic branch sits at")
    print("  - unit distance below the universal twin/non-twin boundary, and")
    print("  - a strictly positive distance below the exact blocked Family-1 ceiling throughout")
    print("    the whole admissible blocked interval.")

    subbanner("IV. Exact corridor read directly as epsilon_4 bounds")
    eps4_twin = sp.expand(2 * (1 + eps2) ** 2 - 1)
    eps4_F1 = sp.expand(4 * cmax * (1 + eps2) ** 2 - 1)
    print("Twin-safe upper bound on epsilon_4 =")
    sp.pprint(eps4_twin)
    print("Family-1 admissibility upper bound on epsilon_4 =")
    sp.pprint(sp.factor(eps4_F1))
    print()
    print("At fixed epsilon_2, the actual grouped-P2 branch must satisfy")
    print("  epsilon_4 <= 2(1+epsilon_2)^2 - 1        for lowest-twin sufficiency,")
    print("  epsilon_4 <  4 cmax^(F1)(eps_blk)(1+epsilon_2)^2 - 1 for exact Family-1 admissibility.")

    banner("STAGE 263 LEDGER")
    print("1. The exact blocked Family-1 corridor can be tested directly in the actual grouped-P2")
    print("   branch variables by the two numerators")
    print("      M_twin = 2(1+epsilon_2)^2 - (1+epsilon_4),")
    print("      M_F1   = 4 cmax^(F1)(eps_blk)(1+epsilon_2)^2 - (1+epsilon_4).")
    print("2. At the actual isotropic point, M_twin = 1 exactly, so the natural branch is not merely")
    print("   twin-safe — it sits one full unit below the universal twin/non-twin boundary.")
    print("3. The exact Family-1 admissibility margin at the actual isotropic point is 4 cmax^(F1)-1,")
    print("   which is strictly positive and increases with blocking on the Lambda_EM-refreshed branch.")
    print("4. So the actual isotropic grouped-P2 branch is safely inside the exact blocked corridor;")
    print("   any future failure must come from real geometry contamination, not from the isotropic point itself.")


if __name__ == "__main__":
    main()
