#!/usr/bin/env python3
"""
5pn_stage59_n5_wall_depth_lock.py

Stage 59 audit: exact n=5 wall-depth lock for the Family-1 branch.
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

banner("STAGE 59 — EXACT n=5 WALL-DEPTH LOCK")

K, rho, m, c_s, hbar, ell, rho_w, lambda_mu = sp.symbols(
    "K rho m c_s hbar ell rho_w lambda_mu", positive=True, real=True
)

subbanner("59.1 — Exact n=5 enthalpy–sound-speed identity")
P = K * rho**5
U = K * rho**5 / 4
h = sp.diff(U, rho)
c_s_sq = sp.diff(P, rho) / m
expect_zero("h - m c_s^2 / 4", h - m * c_s_sq / 4)

subbanner("59.2 — Local wall-depth datum")
mu_star = lambda_mu * m * c_s_sq / 4
Theta_w = sp.simplify(4 * rho_w**2 * mu_star**2 / (hbar**2 * c_s_sq))
print("Theta_w =")
sp.pprint(Theta_w)

subbanner("59.3 — Healing-lock reduction")
Theta_heal = sp.simplify(Theta_w.subs(hbar, 2 * m * sp.sqrt(c_s_sq) * ell))
# Because hbar = 2 m c_s ell with c_s = sqrt(c_s_sq)
print("Theta_w under hbar = 2 m c_s ell =")
sp.pprint(Theta_heal)

# Rewrite by hand into the exact compact form.
Theta_compact = sp.simplify(lambda_mu**2 * rho_w**2 / (16 * ell**2))
expect_zero("Theta_heal - lambda_mu^2 rho_w^2/(16 ell^2)", Theta_heal - Theta_compact)

subbanner("59.4 — Reference-branch normalized form")
Theta_ref = sp.simplify(Theta_compact.subs(ell, sp.Rational(1,20)))
expect_zero("Theta_ref - 25 lambda_mu^2 rho_w^2", Theta_ref - 25 * lambda_mu**2 * rho_w**2)

banner("STAGE 59 FINAL LEDGER")
print("On the explicit Family-1 branch, the wall-depth datum is no longer opaque.")
print("Using the n=5 enthalpy relation and the healing-width lock, one gets")
print("  Theta_w = lambda_mu^2 rho_w^2 / (16 ell^2),")
print("and in normalized Family-1 wall variables with ell = 1/20,")
print("  Theta_w = 25 lambda_mu^2 rho_w^2.")
