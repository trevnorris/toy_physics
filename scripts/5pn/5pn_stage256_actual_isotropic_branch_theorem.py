#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 256 — exact theorem that the actual isotropic one-pole branch is automatically
safe on the Lambda_EM-refreshed Family-1 support/source branch throughout the full
admissible blocked regime.
"""


def refreshed_zeta_max() -> mp.mpf:
    mp.mp.dps = 80
    x01 = mp.besseljzero(0, 1)
    Lambda_EM = mp.sqrt(2) * mp.pi / x01
    Lambda_ell = 20 * Lambda_EM
    eta = Lambda_ell
    kappa = mp.mpf("9") * Lambda_ell**2 / 5
    f = lambda y: y * mp.tan(y) - eta
    lo = mp.mpf("1.4")
    hi = mp.mpf("1.56")
    for _ in range(400):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    y = (lo + hi) / 2
    A_F1 = (kappa + mp.pi**2 / 4) / (kappa + y**2)
    return A_F1 * mp.pi**2 / 4


def main() -> None:
    banner("STAGE 256 — ACTUAL ISOTROPIC BRANCH IS AUTOMATIC ON FAMILY-1")

    eps_blk, zmax = sp.symbols("eps_blk zmax", positive=True, real=True)
    c_pole = sp.Rational(1, 4)
    rho_alpha = sp.simplify(1 / (1 - c_pole))
    zeta_req = sp.simplify((rho_alpha - 1) / (1 - eps_blk * (2 - rho_alpha)))

    subbanner("I. Exact blocked-regime demand of the actual isotropic one-pole branch")
    print("rho_alpha^(act) =")
    sp.pprint(rho_alpha)
    print("zeta_req^(act)(eps_blk) =")
    sp.pprint(zeta_req)
    expect_zero("zeta_req - 1/(3-2 eps_blk)", sp.simplify(zeta_req - 1 / (3 - 2 * eps_blk)))

    print("Pi_tr / C_mix =")
    sp.pprint(rho_alpha)
    print("So Pi_tr / C_mix stays exactly 4/3, independent of blocking.")

    subbanner("II. Exact worst-case demand at the admissible blocking ceiling")
    eps_crit = sp.simplify(1 / zmax)
    zeta_worst = sp.simplify(zeta_req.subs(eps_blk, eps_crit))
    print("eps_blk^crit =")
    sp.pprint(eps_crit)
    print("zeta_req^(act)(eps_blk^crit) =")
    sp.pprint(zeta_worst)
    expect_zero(
        "zeta_max - zeta_worst",
        sp.simplify(zmax - zeta_worst - 3 * zmax * (zmax - 1) / (3 * zmax - 2)),
    )
    print("zeta_max - zeta_worst = 3 zeta_max (zeta_max - 1)/(3 zeta_max - 2) > 0 whenever zeta_max > 1.")

    subbanner("III. Exact twin-safe theorem")
    print("zeta_req^(act) <= 1  iff  eps_blk <= 1.")
    print("Since the admissible Family-1 blocked regime satisfies eps_blk < 1/zeta_max and zeta_max > 1,")
    print("the actual isotropic branch stays in the symmetric-lowest-twin-safe regime throughout.")

    zeta_max_num = refreshed_zeta_max()
    eps_crit_num = 1 / zeta_max_num
    zeta_worst_num = 1 / (3 - 2 * eps_crit_num)
    margin_num = zeta_max_num - zeta_worst_num

    subbanner("IV. Exact Lambda_EM-refreshed numerical margins")
    print(f"zeta_max^(F1) ≈ {zeta_max_num}")
    print(f"eps_blk^crit ≈ {eps_crit_num}")
    print(f"zeta_req^(act)(0) = {mp.mpf('1')/3}")
    print(f"zeta_req^(act)(eps_blk^crit) ≈ {zeta_worst_num}")
    print(f"hard margin at worst admissible blocking ≈ {margin_num}")
    print(f"Pi_tr^(act) / C_mix = {mp.mpf('4')/3}")

    banner("STAGE 256 LEDGER")
    print("1. The actual isotropic one-pole conservative branch has")
    print("      rho_alpha = 4/3,   zeta_req = 1/(3 - 2 eps_blk),   Pi_tr/C_mix = 4/3.")
    print("2. Because zeta_max^(F1) > 1 on the refreshed Lambda_EM Family-1 branch,")
    print("   the whole admissible blocked interval stays safely below the hard support ceiling.")
    print("3. So the explicit Family-1 support/source side is automatic on the actual isotropic")
    print("   branch, not just at eps_blk = 0 but throughout the full admissible blocked regime.")


if __name__ == "__main__":
    main()
