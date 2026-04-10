#!/usr/bin/env python3
"""
5pn_stage56_family1_geometry_map.py

Stage 56 audit: Family-1 reference-branch geometry map.
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

banner("STAGE 56 — FAMILY-1 REFERENCE-BRANCH GEOMETRY MAP")

eps_r = sp.Rational(1, 20)
Lambda_star = sp.Rational(37, 20)
Lambda_ell = sp.simplify(Lambda_star / eps_r)
expect_zero("Lambda_ell - 37", Lambda_ell - 37)

eta = Lambda_ell

print("epsilon_r =")
sp.pprint(eps_r)
print("Lambda_* = L/a =")
sp.pprint(Lambda_star)
print("Lambda_ell = L/ell =")
sp.pprint(Lambda_ell)
print("eta =")
sp.pprint(eta)

banner("STAGE 56 FINAL LEDGER")
print("On the balanced Family-1 reference branch with ell/a = 1/20 and L/a = 37/20,")
print("the support-geometry ratio is fixed exactly to")
print("  Lambda_ell = 37,  eta = 37.")
