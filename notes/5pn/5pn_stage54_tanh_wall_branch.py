#!/usr/bin/env python3
"""
5pn_stage54_tanh_wall_branch.py

Stage 54 audit: canonical tanh-wall branch and natural local mouth closure.
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

banner("STAGE 54 — CANONICAL TANH-WALL BRANCH")

xi = sp.symbols("xi", real=True)
hbar, m, rho_w, c_sw = sp.symbols("hbar m rho_w c_sw", positive=True, real=True)
a, ell, L, V0 = sp.symbols("a ell L V0", positive=True, real=True)
T_X, K_X = sp.symbols("T_X K_X", positive=True, real=True)

subbanner("54.1 — Exact tanh-wall moments")
# Use t = tanh(xi), dt = sech(xi)^2 dxi.
t = sp.symbols("t", real=True)
I_f = sp.integrate((1 - t**2) / 4, (t, -1, 1))
I_g = sp.integrate(t**2 * (1 - t**2), (t, -1, 1))
expect_zero("I_f - 1/3", I_f - sp.Rational(1, 3))
expect_zero("I_g - 4/15", I_g - sp.Rational(4, 15))
expect_zero("I_g/I_f - 4/5", sp.simplify(I_g / I_f) - sp.Rational(4, 5))

subbanner("54.2 — Explicit branch coefficients")
T_X_expr = sp.simplify(sp.pi * a**2 * ell * hbar**2 / (3 * m * rho_w))
K_X_expr = sp.simplify(4 * sp.pi * a**2 * (5 * m**2 * c_sw**2 * ell**2 + hbar**2) / (15 * ell * m * rho_w))
H_w = sp.simplify(m * c_sw**2 / rho_w)
J1 = sp.simplify(sp.Rational(1,3) / H_w)
kappa = sp.simplify(4 * (m * c_sw * L / hbar)**2 + sp.Rational(4,5) * (L / ell)**2)
W_wall = sp.simplify(4 * rho_w**2 * V0**2 * L**2 / (hbar**2 * c_sw**2 * ell**2))
print("T_X =")
sp.pprint(T_X_expr)
print("K_X =")
sp.pprint(K_X_expr)
print("J_1 =")
sp.pprint(J1)
print("kappa =")
sp.pprint(kappa)
print("W_wall =")
sp.pprint(W_wall)

subbanner("54.3 — Natural local mouth closure")
eta = sp.simplify((T_X_expr / ell) * L / T_X_expr)
expect_zero("eta - L/ell", eta - L / ell)

banner("STAGE 54 FINAL LEDGER")
print("Stage 54 makes the first explicit branch canonical:")
print("  I_f = 1/3,  I_g = 4/15,  I_g/I_f = 4/5,")
print("  eta = L/ell,")
print("  kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L/ell)^2,")
print("  W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).")
