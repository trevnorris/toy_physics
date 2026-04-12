#!/usr/bin/env python3
"""
5pn_stage71_loading_ratio_from_minimal_module.py

Stage 71 audit: loading-ratio extraction from the minimal isotropic quadrupole module.
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

banner("STAGE 71 — LOADING RATIO FROM THE MINIMAL ISOTROPIC MODULE")

rho_alpha, Omega_Q, omega = sp.symbols("rho_alpha Omega_Q omega", positive=True, real=True)
c0, c1, C_mix, Pi_tr = sp.symbols("c0 c1 C_mix Pi_tr", positive=True, real=True)

subbanner("71.1 — Natural contact-plus-pole identification")
Y_cons = sp.simplify(1 / rho_alpha + (rho_alpha - 1) / rho_alpha / (1 - omega**2 / Omega_Q**2))
print("Y_Q^cons(omega) =")
sp.pprint(Y_cons)

c0_expr = sp.simplify(1 / rho_alpha)
c1_expr = sp.simplify((rho_alpha - 1) / rho_alpha)
print("c0 =")
sp.pprint(c0_expr)
print("c1 =")
sp.pprint(c1_expr)
expect_zero("c0 + c1 - 1", c0_expr + c1_expr - 1)

subbanner("71.2 — Exact inverse formulas")
rho_from_c0 = sp.simplify(1 / c0)
rho_from_c1 = sp.simplify(1 / (1 - c1))
zeta_from_cc = sp.simplify(c1 / c0)
print("rho_alpha = 1/c0 = 1/(1-c1)")
print("zeta_req = c1/c0")
sp.pprint(rho_from_c0)
sp.pprint(rho_from_c1)
sp.pprint(zeta_from_cc)

subbanner("71.3 — Minimal isotropic quadrupole module")
c0_min = sp.Rational(3, 4)
c1_min = sp.Rational(1, 4)
rho_min = sp.simplify(1 / c0_min)
zeta_min = sp.simplify(c1_min / c0_min)
Pi_min = sp.simplify(rho_min * C_mix)
expect_zero("rho_alpha - 4/3", rho_min - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_min - sp.Rational(1, 3))
expect_zero("Pi_tr - (4/3) C_mix", Pi_min - sp.Rational(4, 3) * C_mix)

banner("STAGE 71 FINAL LEDGER")
print("Under the natural contact-plus-pole interpretation of the minimal isotropic conservative")
print("quadrupole precursor")
print("  Y_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2),")
print("the explicit support loading ratio is fixed exactly to")
print("  rho_alpha = 4/3,  zeta_req = 1/3,  Pi_tr = (4/3) C_mix.")
print("So the branch lies in the symmetric-lowest-twin regime and does not require non-twin asymmetry.")
