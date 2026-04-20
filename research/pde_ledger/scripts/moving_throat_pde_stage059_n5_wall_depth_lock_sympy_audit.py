#!/usr/bin/env python3
"""
Stage 59 SymPy audit.

Checks:
1. Exact n=5 enthalpy-sound-speed identity h = m c_s^2 / 4.
2. Exact algebraic reduction of Theta_w under the local enthalpy lock mu_* = lambda_mu h_w.
3. Exact healing-lock reduction Theta_w = lambda_mu^2 rho_w^2 / (16 ell^2).
4. Reference-branch form Theta_w = (25/4) lambda_mu^2 rho_w^2 in normalized Family-1 wall units.
"""
from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 59 — EXACT n=5 WALL-DEPTH LOCK")

K, rho, m, hbar, csw = sp.symbols("K rho m hbar c_sw", positive=True, real=True)
lambda_mu, rho_w, ell, a = sp.symbols("lambda_mu rho_w ell a", positive=True, real=True)

P = K * rho**5
U = K * rho**5 / 4
cs2 = sp.diff(P, rho) / m
h = sp.diff(U, rho)

print("c_s^2(rho) =", sp.simplify(cs2))
print("h(rho)     =", sp.simplify(h))
expect_zero("n=5 enthalpy identity", h - m * cs2 / 4)

mu_star = lambda_mu * m * csw**2 / 4
Theta_w = sp.simplify(4 * rho_w**2 * mu_star**2 / (hbar**2 * csw**2))
Theta_expected = sp.simplify(lambda_mu**2 * m**2 * rho_w**2 * csw**2 / (4 * hbar**2))
print("Theta_w (enthalpy lock) =", Theta_w)
expect_zero("Theta_w - expected", Theta_w - Theta_expected)

healing_sub = {ell: hbar / (2 * m * csw)}
Theta_heal = sp.simplify(Theta_w.subs(healing_sub))
Theta_heal_expected = sp.simplify(lambda_mu**2 * rho_w**2 / (16 * ell**2))
# verify by substituting back c_sw in terms of ell
cs_sub = {csw: hbar / (2 * m * ell)}
expect_zero("healing-lock reduction", Theta_w.subs(cs_sub) - Theta_heal_expected)
print("Theta_w (healing lock) =", Theta_heal_expected)

ref_sub = {ell: a / 20}
Theta_ref = sp.simplify(Theta_heal_expected.subs(ref_sub))
print("Theta_w (reference branch, general a) =", Theta_ref)
Theta_ref_norm = sp.simplify(Theta_ref.subs(a, 1))
print("Theta_w (reference branch, normalized wall units) =", Theta_ref_norm)
expect_zero("normalized reference factor", Theta_ref_norm - 25 * lambda_mu**2 * rho_w**2)
