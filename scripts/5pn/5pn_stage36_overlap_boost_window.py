#!/usr/bin/env python3
"""
5pn_stage36_overlap_boost_window.py

SymPy audit for Moving-Throat PDE Stage 36.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 36 — EXACT OVERLAP-BOOST WINDOW")

# Symbols
alpha, L, s = sp.symbols("alpha L s", positive=True, real=True)
zeta_req = sp.symbols("zeta_req", positive=True, real=True)
pi = sp.pi

# D/N lowest mode and mixed baseline overlap
chi0 = sp.sqrt(2 / L) * sp.sin(pi * s / (2 * L))
I_W = sp.simplify(sp.integrate(chi0, (s, 0, L)))
print("chi_0(s) =", chi0)
print("I_W      =", I_W)
expect_zero("I_W - 2 sqrt(2L)/pi", sp.simplify(I_W - 2 * sp.sqrt(2 * L) / pi))

# General exponential source density
sigma_alpha = sp.simplify(alpha * sp.exp(alpha * s / L) / (sp.exp(alpha) - 1))
print("\nsigma_alpha(s) =", sigma_alpha)
expect_zero("total source strength - L", sp.simplify(sp.integrate(sigma_alpha, (s, 0, L)) - L))

I_alpha = sp.simplify(sp.integrate(sigma_alpha * chi0, (s, 0, L)))
Omega_exp = sp.simplify(I_alpha / I_W)

print("\nI_alpha =", I_alpha)
print("Omega_exp(alpha) =", Omega_exp)

expect_zero(
    "Omega_exp(alpha) - closed form",
    Omega_exp
    - pi * alpha * (2 * alpha * sp.exp(alpha) + pi)
    / ((4 * alpha**2 + pi**2) * (sp.exp(alpha) - 1)),
)

print("\nmax chi_0 =", sp.simplify(sp.sqrt(2 / L)))
print("sharp overlap ceiling =", sp.simplify((sp.sqrt(2 / L) * L) / I_W))

expect_zero("Omega_exp(0) - 1", sp.simplify(sp.limit(Omega_exp, alpha, 0) - 1))
expect_zero("Omega_exp(infty) - pi/2", sp.simplify(sp.limit(Omega_exp, alpha, sp.oo) - pi / 2))
print("small-alpha series =", sp.series(Omega_exp, alpha, 0, 3).removeO())
expect_zero(
    "large-alpha asymptotic coefficient",
    sp.simplify(sp.limit((pi / 2 - Omega_exp) * alpha**2, alpha, sp.oo) - pi**3 / 8),
)

banner("STAGE 36 THEOREM LEDGER")
print("For any nonnegative lowest-lane source density with the same total source strength as")
print("the mixed lane, the overlap boost satisfies")
print("  0 <= Omega_0 <= pi/2,")
print("so")
print("  A_I = Omega_0^2 <= pi^2/4.")
print()
print("The explicit exponential constructive family has")
print("  Omega_exp(0) = 1,")
print("  lim_(alpha->+infty) Omega_exp = pi/2,")
print("and therefore pure overlap rescue alone can beat the Stage-35 threshold only if")
print("  zeta_req <= pi^2/4.")
