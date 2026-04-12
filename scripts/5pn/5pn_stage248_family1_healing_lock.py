#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 248 — healing-length lock and exact support-scale reduction on the explicit
Family-1 reference branch.
"""


def main() -> None:
    banner("STAGE 248 — FAMILY-1 HEALING LOCK")

    hbar, m, csw = sp.symbols("hbar m csw", positive=True, real=True)
    Lambda_ell = sp.Integer(37)

    subbanner("I. Healing-width closure")
    ell_h = sp.simplify(hbar / (2 * m * csw))
    print("ell_h =")
    sp.pprint(ell_h)

    chi_s = sp.simplify(m * csw * (Lambda_ell * ell_h) / hbar)
    print()
    print("chi_s =")
    sp.pprint(chi_s)
    expect_zero("chi_s = Lambda_ell/2", sp.simplify(chi_s - sp.Rational(Lambda_ell, 2)))

    subbanner("II. Exact reference-branch support scales")
    kappa = sp.simplify(4 * chi_s**2 + sp.Rational(4, 5) * Lambda_ell**2)
    alpha = sp.simplify(sp.sqrt(kappa))

    print("kappa =")
    sp.pprint(kappa)
    print()
    print("alpha = sqrt(kappa) =")
    sp.pprint(alpha)

    expect_zero("kappa = 12321/5", sp.simplify(kappa - sp.Rational(12321, 5)))
    expect_zero("alpha = 111/sqrt(5)", sp.simplify(alpha - 111 / sp.sqrt(5)))

    banner("STAGE 248 LEDGER")
    print("1. Identifying the active shell width with the GNLS healing width gives")
    print("      ell = hbar / (2 m c_(s,w)).")
    print("2. Therefore the explicit Family-1 support scale is locked to")
    print("      chi_s = L/(2 ell) = Lambda_ell/2 = 37/2.")
    print("3. On the same branch,")
    print("      kappa = 4 chi_s^2 + (4/5) Lambda_ell^2 = 12321/5,")
    print("      alpha = sqrt(kappa) = 111/sqrt(5).")
    print("4. So after the healing lock, the explicit branch has fixed (chi_s, eta, kappa) and only the wall-loading datum remains open.")


if __name__ == "__main__":
    main()
