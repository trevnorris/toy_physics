#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 250 — replace the shorthand Family-1 reference ratio L/a = 37/20 by the
exact carried EM-branch geometry relation

    Lambda_EM = sqrt(2) * pi / x_01,

where x_01 is the first positive J_0 zero.

This script is deliberately narrow: it refreshes the geometric lock, the healing
lock, and the resulting explicit support variables before the threshold window and
wall-depth comparison are recomputed downstream.
"""


def main() -> None:
    banner("STAGE 250 — EXACT LAMBDA_EM FAMILY-1 GEOMETRY REFRESH")

    x01 = sp.symbols("x_01", positive=True, real=True)
    eps_r = sp.Rational(1, 20)

    Lambda_EM = sp.sqrt(2) * sp.pi / x01
    Lambda_ell = sp.simplify(Lambda_EM / eps_r)
    eta = sp.simplify(Lambda_ell)
    chi_s = sp.simplify(Lambda_ell / 2)
    kappa = sp.simplify(4 * chi_s**2 + sp.Rational(4, 5) * Lambda_ell**2)

    subbanner("I. Exact carried geometry relation")
    print("Lambda_EM =")
    sp.pprint(Lambda_EM)
    print("Lambda_ell =")
    sp.pprint(Lambda_ell)
    print("chi_s =")
    sp.pprint(chi_s)
    print("kappa =")
    sp.pprint(kappa)

    expect_zero(
        "Lambda_ell - 20*sqrt(2)*pi/x_01",
        sp.simplify(Lambda_ell - 20 * sp.sqrt(2) * sp.pi / x01),
    )
    expect_zero(
        "chi_s - 10*sqrt(2)*pi/x_01",
        sp.simplify(chi_s - 10 * sp.sqrt(2) * sp.pi / x01),
    )
    expect_zero(
        "kappa - 1440*pi^2/x_01^2",
        sp.simplify(kappa - 1440 * sp.pi**2 / x01**2),
    )

    subbanner("II. Numerical carried value from the exact EM branch")
    mp.mp.dps = 80
    x01_num = mp.besseljzero(0, 1)
    Lambda_EM_num = mp.sqrt(2) * mp.pi / x01_num
    Lambda_ell_num = 20 * Lambda_EM_num
    chi_s_num = Lambda_ell_num / 2
    kappa_num = mp.mpf('9') * Lambda_ell_num**2 / 5

    Lambda_freeze = mp.mpf('37') / 20
    Lambda_ell_freeze = mp.mpf('37')
    kappa_freeze = mp.mpf(12321) / 5

    print(f"x_01 ≈ {x01_num}")
    print(f"Lambda_EM ≈ {Lambda_EM_num}")
    print(f"Lambda_ell ≈ {Lambda_ell_num}")
    print(f"chi_s ≈ {chi_s_num}")
    print(f"kappa ≈ {kappa_num}")
    print()
    print(f"relative shift in L/a from 37/20 ≈ {(Lambda_EM_num - Lambda_freeze) / Lambda_freeze}")
    print(f"relative shift in L/ell from 37 ≈ {(Lambda_ell_num - Lambda_ell_freeze) / Lambda_ell_freeze}")
    print(f"absolute shift in kappa from 12321/5 ≈ {kappa_num - kappa_freeze}")

    banner("STAGE 250 LEDGER")
    print("1. The reference Family-1 branch should use the exact carried geometry ratio")
    print("      L/a = Lambda_EM = sqrt(2) pi / x_01,")
    print("   not the shorthand freeze 37/20.")
    print("2. With ell/a = 1/20, the explicit branch variables become")
    print("      Lambda_ell = 20 Lambda_EM,")
    print("      chi_s = 10 Lambda_EM,")
    print("      kappa = (9/5) Lambda_ell^2 = 1440 pi^2 / x_01^2.")
    print("3. Numerically this is only a ~1.36e-3 relative shift from the old shorthand,")
    print("   but it should be propagated exactly in all later explicit-branch formulas.")


if __name__ == "__main__":
    main()
