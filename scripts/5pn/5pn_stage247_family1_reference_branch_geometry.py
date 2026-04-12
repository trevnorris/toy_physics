#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 247 — Family-1 reference-branch geometry map.
"""


def main() -> None:
    banner("STAGE 247 — FAMILY-1 REFERENCE-BRANCH GEOMETRY MAP")

    eps_r = sp.Rational(1, 20)
    Lambda_star = sp.Rational(37, 20)

    subbanner("I. Exact reference-branch wall width")
    print("epsilon_r =")
    sp.pprint(eps_r)
    print()
    print("L/a =")
    sp.pprint(Lambda_star)

    Lambda_ell = sp.simplify(Lambda_star / eps_r)
    print()
    print("Lambda_ell = L/ell =")
    sp.pprint(Lambda_ell)
    expect_zero("Lambda_ell = 37", sp.simplify(Lambda_ell - 37))

    subbanner("II. Immediate Robin consequence")
    eta = sp.simplify(Lambda_ell)
    expect_zero("eta = 37", sp.simplify(eta - 37))
    print("eta =")
    sp.pprint(eta)

    banner("STAGE 247 LEDGER")
    print("1. On the balanced Family-1 reference branch, the thin-wall identification is")
    print("      ell/a = 1/20.")
    print("2. The carried preferred aspect ratio is")
    print("      L/a = 37/20.")
    print("3. Therefore the explicit branch fixes")
    print("      Lambda_ell = L/ell = 37,")
    print("      eta = 37.")
    print("4. So the first explicit moving-throat support branch is pinned to one concrete large-eta Robin regime.")


if __name__ == "__main__":
    main()
