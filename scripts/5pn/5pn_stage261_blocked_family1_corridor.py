#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage258_261_common import banner, subbanner, expect_zero, family1_refreshed_numbers

"""
Stage 261 — exact blocked Family-1 admissible non-twin corridor after the
Lambda_EM refresh.

This script combines the Stage-254 Family-1 support ceiling with the Stage-258
regime split and the Stage-260 blocking monotonicity theorem.

The main exact results are:

1. the hard Family-1 pole ceiling at fixed blocking is

       c_pole,max(eps_blk)
         = zeta_max (1 - eps_blk)
           / [1 + (1 - 2 eps_blk) zeta_max],

2. because zeta_max > 1, this ceiling increases monotonically with blocking,
3. the lowest-twin boundary c_pole = 1/2 stays fixed,
4. so the admissible non-twin Family-1 corridor

       1/2 < c_pole < c_pole,max(eps_blk)

   widens with blocking,
5. while the maximum required asymmetry excess stays fixed at

       Delta_zeta,max = zeta_max - 1,

   because the Family-1 ceiling is set at fixed zeta_max in zeta-language.

This is the exact blocked analogue of the unblocked Family-1 headroom statement
from Stage 257.
"""


def main() -> None:
    banner("STAGE 261 — EXACT BLOCKED FAMILY-1 NON-TWIN CORRIDOR")

    zmax, eps_blk = sp.symbols("zeta_max eps_blk", positive=True, real=True)
    cmax = sp.simplify(zmax * (1 - eps_blk) / (1 + (1 - 2 * eps_blk) * zmax))
    dcmax_de = sp.factor(sp.diff(cmax, eps_blk))
    width_nt = sp.simplify(cmax - sp.Rational(1, 2))
    dwidth_de = sp.factor(sp.diff(width_nt, eps_blk))

    subbanner("I. Exact blocked Family-1 pole ceiling")
    print("c_pole,max^(F1)(eps_blk) =")
    sp.pprint(sp.factor(cmax))
    print("d c_pole,max / d eps_blk =")
    sp.pprint(dcmax_de)
    print("non-twin corridor width = c_pole,max - 1/2 =")
    sp.pprint(sp.factor(width_nt))
    print("d(width)/d eps_blk =")
    sp.pprint(dwidth_de)

    expect_zero(
        "dcmax/de - positive form",
        sp.simplify(dcmax_de - zmax * (zmax - 1) / (1 + (1 - 2 * eps_blk) * zmax) ** 2),
    )

    subbanner("II. Exact endpoint values")
    cmax0 = sp.simplify(cmax.subs(eps_blk, 0))
    eps_crit = sp.simplify(1 / zmax)
    cmaxcrit = sp.simplify(cmax.subs(eps_blk, eps_crit))
    print("c_pole,max^(F1)(0) =")
    sp.pprint(cmax0)
    print("eps_blk^crit =")
    sp.pprint(eps_crit)
    print("c_pole,max^(F1)(eps_blk^crit) =")
    sp.pprint(cmaxcrit)
    expect_zero("cmax(0)-z/(1+z)", sp.simplify(cmax0 - zmax / (1 + zmax)))
    expect_zero("cmax(eps_crit)-1", sp.simplify(cmaxcrit - 1))

    subbanner("III. Exact non-twin corridor in geometry-contamination language")
    eps2, eps4 = sp.symbols("epsilon_2 epsilon_4", real=True)
    twin_lower = sp.simplify(2 * (1 + eps2) ** 2)
    upper = sp.simplify(4 * cmax * (1 + eps2) ** 2)
    print("Lower non-twin boundary: 1 + epsilon_4 =")
    sp.pprint(twin_lower)
    print("Upper Family-1 boundary: 1 + epsilon_4 =")
    sp.pprint(sp.factor(upper))
    print("So the blocked admissible non-twin strip is")
    print("  2(1+epsilon_2)^2 < 1+epsilon_4 < 4 c_pole,max^(F1)(eps_blk)(1+epsilon_2)^2.")

    subbanner("IV. Exact asymmetry-ceiling invariance")
    c_pole = sp.symbols("c_pole", real=True)
    zeta_req = sp.simplify(c_pole / (1 - eps_blk - (1 - 2 * eps_blk) * c_pole))
    Delta_zeta = sp.simplify(zeta_req - 1)
    Delta_max = sp.simplify(Delta_zeta.subs(c_pole, cmax))
    print("Delta_zeta,max on the blocked Family-1 ceiling =")
    sp.pprint(sp.factor(Delta_max))
    expect_zero("Delta_max-(zmax-1)", sp.simplify(Delta_max - (zmax - 1)))

    subbanner("V. Refreshed numerical window")
    nums = family1_refreshed_numbers()
    zmax_num = mp.mpf(nums["zeta_max"])
    cmax0_num = mp.mpf(nums["c_pole_max_0"])
    epscrit_num = mp.mpf(nums["eps_blk_crit"])
    width0_num = cmax0_num - mp.mpf('1')/2
    widthcrit_num = mp.mpf('1')/2
    upper0 = 4 * cmax0_num
    uppercrit = mp.mpf('4')
    print(f"zeta_max^(F1) ≈ {zmax_num}")
    print(f"c_pole,max^(F1)(0) ≈ {cmax0_num}")
    print(f"eps_blk^crit ≈ {epscrit_num}")
    print(f"Initial non-twin corridor width  c_max(0)-1/2 ≈ {width0_num}")
    print(f"Critical-limit width            c_max(eps_crit)-1/2 = {widthcrit_num}")
    print(f"Upper geometry factor at eps_blk=0      : 4 c_max(0) ≈ {upper0}")
    print(f"Upper geometry factor at eps_blk=epscrit: 4 c_max(eps_crit) = {uppercrit}")
    print(f"Maximum asymmetry excess Delta_zeta,max = zeta_max-1 ≈ {zmax_num - 1}")

    banner("STAGE 261 LEDGER")
    print("1. The exact Lambda_EM-refreshed Family-1 hard ceiling in c_pole language grows with")
    print("   blocking whenever zeta_max^(F1) > 1:")
    print("      d c_pole,max^(F1) / d eps_blk > 0.")
    print("2. The universal lowest-twin boundary c_pole = 1/2 stays fixed, so the admissible")
    print("   Family-1 non-twin corridor widens with blocking:")
    print("      1/2 < c_pole < c_pole,max^(F1)(eps_blk).")
    print("3. In geometry-contamination language this is the exact blocked strip")
    print("      2(1+epsilon_2)^2 < 1+epsilon_4 < 4 c_pole,max^(F1)(eps_blk)(1+epsilon_2)^2.")
    print("4. The maximal asymmetry demand on that blocked Family-1 ceiling stays fixed at")
    print("      Delta_zeta,max = zeta_max^(F1) - 1,")
    print("   because the ceiling itself is defined at fixed zeta_max in support language.")


if __name__ == "__main__":
    main()
