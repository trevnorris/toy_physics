#!/usr/bin/env python3
"""
Stage 108 SymPy audit: positive local mouth-source theorem.
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


banner("STAGE 108 — POSITIVE LOCAL MOUTH-SOURCE THEOREM")

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

print("\nConclusion: under any positive localized mouth source law,")
print("the upper compensated branch is impossible and the lower branch is unique.")
