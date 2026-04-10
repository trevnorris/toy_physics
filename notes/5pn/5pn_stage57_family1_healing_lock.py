#!/usr/bin/env python3
"""
5pn_stage57_family1_healing_lock.py

Stage 57 audit: healing-length lock and the actual reference-branch support scale.
"""

from __future__ import annotations

import math
import mpmath as mp
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 57 — HEALING-LENGTH LOCK ON THE FAMILY-1 BRANCH")

Lambda_ell = sp.Integer(37)
chi_s = sp.simplify(Lambda_ell / 2)
kappa = sp.simplify(4 * chi_s**2 + sp.Rational(4,5) * Lambda_ell**2)
expect_zero("kappa - (9/5) Lambda_ell^2", kappa - sp.Rational(9,5) * Lambda_ell**2)
expect_zero("kappa - 12321/5", kappa - sp.Rational(12321, 5))

alpha = sp.sqrt(kappa)
expect_zero("alpha - 111/sqrt(5)", alpha - sp.Rational(111,1) / sp.sqrt(5))

print("Lambda_ell =")
sp.pprint(Lambda_ell)
print("chi_s =")
sp.pprint(chi_s)
print("kappa =")
sp.pprint(kappa)
print("alpha =")
sp.pprint(alpha)
print("alpha numeric =", mp.sqrt(mp.mpf(12321)/5))

banner("STAGE 57 FINAL LEDGER")
print("With the healing-width closure ell = hbar/(2 m c_(s,w)), the explicit Family-1 branch fixes")
print("  chi_s = Lambda_ell/2 = 37/2,")
print("  kappa = 12321/5,")
print("  alpha = 111/sqrt(5) ≈ 49.6407091.")
