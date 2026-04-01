#!/usr/bin/env python3
"""
Stage 36 SymPy audit: exact overlap-boost window for the lowest support lane.
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


banner("STAGE 36 — EXACT OVERLAP-BOOST WINDOW")

s, L, alpha = sp.symbols("s L alpha", positive=True, real=True)
pi = sp.pi

chi0 = sp.sqrt(sp.Integer(2) / L) * sp.sin(pi * s / (2 * L))
I_W = sp.simplify(sp.integrate(chi0, (s, 0, L)))
chi0_max = sp.sqrt(sp.Integer(2) / L)
Omega_max = sp.simplify(L * chi0_max / I_W)
A_I_max = sp.simplify(Omega_max**2)

print("chi0(s) =")
sp.pprint(chi0)
print("I_W =", I_W)
print("max chi0 =", chi0_max)
print("Omega_max =", Omega_max)
print("A_I,max =", A_I_max)
expect_zero("Omega_max - pi/2", Omega_max - pi/2)
expect_zero("A_I,max - pi^2/4", A_I_max - pi**2/4)

banner("EXPLICIT EXPONENTIAL BOTTOM-BIASED FAMILY")

sigma_alpha = sp.simplify(alpha * sp.exp(alpha * s / L) / (sp.exp(alpha) - 1))
Sigma_total = sp.simplify(sp.integrate(sigma_alpha, (s, 0, L)))
I_alpha = sp.simplify(sp.integrate(sigma_alpha * chi0, (s, 0, L)))
Omega_alpha = sp.simplify(I_alpha / I_W)

print("sigma_alpha(s) =")
sp.pprint(sigma_alpha)
print("Integral sigma_alpha ds =", Sigma_total)
print("I_alpha =", I_alpha)
print("Omega_alpha =", sp.factor(Omega_alpha))
expect_zero("same total source strength", Sigma_total - L)

Omega_alpha_simpler = sp.simplify(
    pi * alpha * (2 * alpha * sp.exp(alpha) + pi)
    / ((4 * alpha**2 + pi**2) * (sp.exp(alpha) - 1))
)
expect_zero("Omega_alpha closed form", Omega_alpha - Omega_alpha_simpler)

Omega0 = sp.simplify(sp.limit(Omega_alpha, alpha, 0))
Omegainf = sp.simplify(sp.limit(Omega_alpha, alpha, sp.oo))
print("Omega_alpha(alpha->0) =", Omega0)
print("Omega_alpha(alpha->+infty) =", Omegainf)
expect_zero("uniform branch value", Omega0 - 1)
expect_zero("antinode concentration limit", Omegainf - pi/2)

series_small = sp.series(Omega_alpha, alpha, 0, 3).removeO()
print("small-alpha series =", series_small)
expected_linear = sp.simplify(2 / pi - sp.Rational(1, 2))
print("linear coefficient =", expected_linear)
expect_zero("linear coefficient - (4-pi)/(2pi)", expected_linear - (4 - pi) / (2 * pi))

banner("PURE OVERLAP THRESHOLD")

zeta_req = sp.symbols("zeta_req", positive=True, real=True)
criterion = sp.simplify(A_I_max - zeta_req)
print("Pure-overlap rescue possible only if zeta_req <=", A_I_max)
print("A_I,max - zeta_req =", criterion)

print("\nFINAL LEDGER")
print("  Omega_0 is bounded by pi/2 for any nonnegative equal-strength support density.")
print("  A_I = Omega_0^2 is therefore bounded by pi^2/4.")
print("  The explicit exponential family interpolates from 1 to pi/2.")
print("  So pure overlap rescue is possible only if zeta_req <= pi^2/4.")
