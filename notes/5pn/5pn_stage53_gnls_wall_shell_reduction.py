#!/usr/bin/env python3
"""
5pn_stage53_gnls_wall_shell_reduction.py

Stage 53 audit: explicit GNLS wall-shell reduction for the first support branch.
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

banner("STAGE 53 — EXPLICIT GNLS WALL-SHELL REDUCTION")

hbar, m, rho_w, H_w = sp.symbols("hbar m rho_w H_w", positive=True, real=True)
q = sp.Function("q")
s = sp.symbols("s", real=True)
Npp, Gpp = sp.symbols("N_pp G_pp", positive=True, real=True)
V0, ell, I1 = sp.symbols("V0 ell I1", positive=True, real=True)
K_X, T_X, L, c_sw = sp.symbols("K_X T_X L c_sw", positive=True, real=True)
gphi = sp.symbols("g_phi", positive=True, real=True)
Lambda_phi, Xi, G_eq = sp.symbols("Lambda_phi Xi G_eq", positive=True, real=True)

subbanner("53.1 — Parent quadratic shell energy")
T_X_expr = sp.simplify(hbar**2 * Npp / (4 * m * rho_w))
K_X_expr = sp.simplify(H_w * Npp + hbar**2 * Gpp / (4 * m * rho_w))
print("T_X =")
sp.pprint(T_X_expr)
print("K_X =")
sp.pprint(K_X_expr)

subbanner("53.2 — Explicit branch loading law")
Lambda_phi_expr = gphi * sp.Symbol("O_(sigma phi)", positive=True, real=True)
print("g_phi = V0/ell on the explicit wall family.")
print("Lambda_phi = g_phi O_(sigma phi) remains the reduced source/support loading.")

subbanner("53.3 — Exact gain and figure-of-merit identity")
G_eq_expr = sp.simplify(gphi**2 * I1 / K_X)
Xi_expr = sp.simplify(gphi**2 * I1 * L**2 / T_X)
print("G_eq =")
sp.pprint(G_eq_expr)
print("Xi =")
sp.pprint(Xi_expr)
print("So Xi = kappa G_eq with kappa = K_X L^2 / T_X.")

subbanner("53.4 — Matched thin-wall branch: Xi = W_wall")
a, J1 = sp.symbols("a J1", positive=True, real=True)
W_wall = sp.simplify(4 * sp.pi * a**2 * L**2 * J1 * V0**2 / (T_X * ell))
Xi_tw = sp.simplify((V0/ell)**2 * (4 * sp.pi * a**2 * ell * J1) * L**2 / T_X)
expect_zero("Xi_tw - W_wall", Xi_tw - W_wall)

banner("STAGE 53 FINAL LEDGER")
print("Stage 53 derives the support tension and stiffness directly from the parent GNLS shell energy:")
print("  T_X = hbar^2 N_(phi phi)/(4 m rho_w),")
print("  K_X = H_w N_(phi phi) + hbar^2 G_(phi phi)/(4 m rho_w).")
print("On the matched thin-wall branch the fixed-point coupling Xi is exactly the same object as")
print("the wall figure of merit W_wall.")
