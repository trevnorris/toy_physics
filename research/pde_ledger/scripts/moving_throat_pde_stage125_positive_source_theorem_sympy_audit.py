#!/usr/bin/env python3
"""
Stage 125 SymPy audit: positive local mouth-source theorem.
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


def expect_true(name: str, cond: bool) -> None:
    print(f"{name} = {cond}")
    if not cond:
        raise AssertionError(f"{name} is false")


banner("STAGE 125 — POSITIVE LOCAL MOUTH-SOURCE THEOREM")

z, L = sp.symbols("z L", positive=True, real=True)
x = sp.symbols("x", real=True)
k = sp.pi / (2 * L)
kernel = sp.cos(sp.pi * x / 2)

print("D/N half-wave derivative factor:")
print("cos(k z) with k = pi/(2L) =", sp.cos(k * z))

kernel_min = sp.calculus.util.minimum(kernel, x, sp.Interval(0, 1))
kernel_max = sp.calculus.util.maximum(kernel, x, sp.Interval(0, 1))
expect_zero("kernel minimum on [0,1]", kernel_min)
expect_zero("kernel maximum on [0,1] - 1", kernel_max - 1)

print("\nBecause the normalized kernel lies in [0,1] on x in [0,1],")
print("every positive normalized source law has its cosine moment in [0,1].")

r = sp.simplify(sp.sqrt(sp.Rational(12, 1) * sp.Rational(37, 20) ** 2 / sp.pi**2 - 1))
R = sp.sqrt(4107 - 100 * sp.pi**2)
gminus = sp.simplify((2 * R - 37 * sp.sqrt(3)) / (20 * sp.pi))
gplus = sp.simplify((2 * R + 37 * sp.sqrt(3)) / (20 * sp.pi))

expect_zero("r - sqrt(4107 - 100 pi^2)/(10 pi)", r - R / (10 * sp.pi))
expect_zero("lower branch balance relation", 1 + r**2 - 4 * (gminus - r) ** 2)
expect_zero("upper branch balance relation", 1 + r**2 - 4 * (gplus - r) ** 2)

print("\nExplicit Family-1 compensated branches:")
print("g_- =", gminus)
print("g_+ =", gplus)
print("g_- (numeric) =", sp.N(gminus, 20))
print("g_+ (numeric) =", sp.N(gplus, 20))

expect_true("g_- > 0", bool(sp.N(gminus) > 0))
expect_true("g_- < 1", bool(sp.N(gminus) < 1))
expect_true("g_+ > 1", bool(sp.N(gplus) > 1))

# Integral-bound test on a one-parameter family of positive normalized sources.
# sigma_a(z) = (a + 1) * (z/L)^a / L is nonneg on [0,L] for a >= 0 and integrates to 1.
# Endpoint values: a = 0 -> uniform (g = 2/pi); a -> oo -> peaked at z = L (g -> 0).
sigma_a = sp.symbols("sigma_param", nonnegative=True, real=True)
sigma_profile = (sigma_a + 1) * (z / L) ** sigma_a / L
norm_a = sp.simplify(sp.integrate(sigma_profile, (z, 0, L)))
expect_zero("parametric family normalization", norm_a - 1)

g_a = sp.simplify(sp.integrate(sigma_profile * sp.cos(sp.pi * z / (2 * L)), (z, 0, L)))
print("g_a (parametric moment) =", g_a)

expect_zero("moment g[uniform] - 2/pi", g_a.subs(sigma_a, 0) - 2 / sp.pi)

# Numerical check that g_a -> 0 as a -> oo (SymPy can't take a symbolic limit on the
# resulting hypergeometric form, so verify the trend at a = 100 instead).
g_a_large = sp.N(g_a.subs(sigma_a, 100))
print("g_a at sigma_param = 100 (peaked-at-L proxy) =", g_a_large)
# Genuine range check of the paper bound 0 <= g[sigma] <= 1 at the peaked endpoint.
# NOT wrapped in abs(): a small NEGATIVE moment (a sign error in the closed form)
# now FAILS the lower bound, instead of passing a mere smallness test.
expect_true("g[peaked@L proxy a=100] >= 0", bool(g_a_large >= 0))
expect_true("g[peaked@L proxy a=100] <= 1", bool(g_a_large <= 1))
# Retain the smallness fact (peaked source biases g toward the lower endpoint):
# a strictly positive value below 1/20 confirms the trend without admitting negatives.
expect_true("g[peaked@L proxy a=100] < 1/20", bool(g_a_large < sp.Rational(1, 20)))

expect_true("g[uniform] >= 0", bool(sp.N(g_a.subs(sigma_a, 0)) >= 0))
expect_true("g[uniform] <= 1", bool(sp.N(g_a.subs(sigma_a, 0)) <= 1))

print("\nConclusion: under any positive localized mouth source law,")
print("the upper compensated branch is impossible and the lower branch is unique.")
