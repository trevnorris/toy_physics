#!/usr/bin/env python3
"""
Stage 29 SymPy audit.

Checks:
1. Exact tracking-branch formulas G_tr and F_tr.
2. Flat-branch and strong-split endpoint limits.
3. Exact derivatives with respect to R.
4. Exact flat-branch comparison formulas G_tr - G_flat and F_flat - F_tr.
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


def expect_positive_coefficients(name: str, expr: sp.Expr, *gens: sp.Symbol) -> None:
    poly = sp.Poly(expr, *gens)
    coeffs = poly.coeffs()
    print(f"{name} coefficients = {coeffs}")
    if any(coef <= 0 for coef in coeffs):
        raise AssertionError(f"{name} has a non-positive coefficient")


banner("STAGE 29 — TRACKING-BRANCH BOUNDS AUDIT")

xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)

G_tr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R ** 2) * xi))
F_tr = sp.simplify(
    (9 * delta + (9 + 2 * R ** 2) * xi) ** 2
    * (9 * delta + (9 + 2 * R) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta ** 2 + 18 * delta * xi + (9 + 2 * R ** 2) * xi ** 2) ** 2)
)

G_flat = sp.simplify(G_tr.subs(R, 1))
F_flat = sp.simplify(F_tr.subs(R, 1))

banner("1. Flat and strong-split endpoints")
print("G_tr =", G_tr)
print("F_tr =", F_tr)
print("G_flat =", G_flat)
print("F_flat =", F_flat)
expect_zero("strong-split endpoint for G", sp.simplify(G_tr.subs(R, 0) - xi))
expect_zero("strong-split endpoint for F", sp.simplify(F_tr.subs(R, 0) - 1 / (1 - xi)))

banner("2. Exact derivatives with respect to R")
dG_dR = sp.factor(sp.diff(G_tr, R))
dF_dR = sp.factor(sp.diff(F_tr, R))
print("dG_tr/dR =", dG_dR)
print("dF_tr/dR =", dF_dR)

P_R = (
    4 * R ** 4 * xi ** 3
    + 54 * R ** 2 * delta ** 2 * xi
    + 90 * R ** 2 * delta * xi ** 2
    + 36 * R ** 2 * xi ** 3
    + 162 * R * delta ** 3
    + 324 * R * delta ** 2 * xi
    + 162 * R * delta * xi ** 2
    + 81 * delta ** 3
    + 243 * delta ** 2 * xi
    + 243 * delta * xi ** 2
    + 81 * xi ** 3
)

dG_expected = -36 * R * xi ** 2 * (delta + xi) / (2 * R ** 2 * xi + 9 * delta + 9 * xi) ** 2
dF_expected = (
    4 * xi * (2 * R * xi + 9 * delta + 9 * xi) * (2 * R ** 2 * xi + 9 * delta + 9 * xi) * P_R
    / (81 * (1 - xi) * (2 * R ** 2 * xi ** 2 + 9 * delta ** 2 + 18 * delta * xi + 9 * xi ** 2) ** 3)
)

expect_zero("dG_tr/dR formula", sp.simplify(dG_dR - dG_expected))
expect_zero("dF_tr/dR formula", sp.simplify(dF_dR - dF_expected))
expect_positive_coefficients("P_R", P_R, R, delta, xi)

banner("3. Exact comparison with the flat branch")

delta_G = sp.factor(sp.simplify(G_tr - G_flat))
delta_F = sp.factor(sp.simplify(F_flat - F_tr))
print("G_tr - G_flat =", delta_G)
print("F_flat - F_tr =", delta_F)

P1 = (
    18 * R ** 2 * delta ** 2 * xi
    + 36 * R ** 2 * delta * xi ** 2
    + 22 * R ** 2 * xi ** 3
    + 81 * R * delta ** 3
    + 180 * R * delta ** 2 * xi
    + 99 * R * delta * xi ** 2
    + 162 * delta ** 3
    + 423 * delta ** 2 * xi
    + 360 * delta * xi ** 2
    + 99 * xi ** 3
)

P2 = (
    18 * R ** 3 * delta ** 2 * xi ** 2
    + 36 * R ** 3 * delta * xi ** 3
    + 22 * R ** 3 * xi ** 4
    + 81 * R ** 2 * delta ** 3 * xi
    + 324 * R ** 2 * delta ** 2 * xi ** 2
    + 459 * R ** 2 * delta * xi ** 3
    + 220 * R ** 2 * xi ** 4
    + 81 * R * delta ** 3 * xi
    + 243 * R * delta ** 2 * xi ** 2
    + 261 * R * delta * xi ** 3
    + 99 * R * xi ** 4
    + 729 * delta ** 4
    + 3078 * delta ** 3 * xi
    + 4959 * delta ** 2 * xi ** 2
    + 3600 * delta * xi ** 3
    + 990 * xi ** 4
)

G_diff_expected = 18 * xi ** 2 * (1 - R ** 2) * (delta + xi) / (
    (9 * delta + 11 * xi) * (2 * R ** 2 * xi + 9 * delta + 9 * xi)
)
F_diff_expected = 4 * xi * (1 - R) * P1 * P2 / (
    81 * (1 - xi) * (9 * delta ** 2 + 18 * delta * xi + 11 * xi ** 2) ** 2
    * (2 * R ** 2 * xi ** 2 + 9 * delta ** 2 + 18 * delta * xi + 9 * xi ** 2) ** 2
)

expect_zero("G_tr - G_flat formula", sp.simplify(delta_G - G_diff_expected))
expect_zero("F_flat - F_tr formula", sp.simplify(delta_F - F_diff_expected))
expect_positive_coefficients("P1", P1, R, delta, xi)
expect_positive_coefficients("P2", P2, R, delta, xi)

banner("3b. Sign verification of branch-difference factors")

# Boundary values: at R=1 the tracking and flat branches must coincide; at R=0
# the tracking branch must hit the strong-split endpoint.
expect_zero(
    "G_tr - G_flat vanishes at R=1",
    sp.simplify((G_tr - G_flat).subs(R, 1)),
)
expect_zero(
    "F_flat - F_tr vanishes at R=1",
    sp.simplify((F_flat - F_tr).subs(R, 1)),
)
expect_zero(
    "G_tr at R=0 equals xi",
    sp.simplify(G_tr.subs(R, 0) - xi),
)
expect_zero(
    "F_tr at R=0 equals 1/(1-xi)",
    sp.simplify(F_tr.subs(R, 0) - 1 / (1 - xi)),
)

# Numerical sign sampling of delta_G = G_tr - G_flat and delta_F = F_flat - F_tr
# at three representative interior points (R, xi, delta). All must be strictly
# positive for the bound G_flat <= G_tr and F_tr <= F_flat to hold on the open
# interval R in (0,1).
delta_G_value = sp.simplify(G_tr - G_flat)
delta_F_value = sp.simplify(F_flat - F_tr)
sample_points = [
    {R: sp.Rational(1, 4), xi: sp.Rational(1, 3), delta: sp.Rational(1, 2)},
    {R: sp.Rational(1, 2), xi: sp.Rational(1, 2), delta: sp.Rational(1, 4)},
    {R: sp.Rational(3, 4), xi: sp.Rational(1, 5), delta: sp.Rational(2, 3)},
]
for idx, pt in enumerate(sample_points, start=1):
    g_sample = sp.simplify(delta_G_value.subs(pt))
    f_sample = sp.simplify(delta_F_value.subs(pt))
    print(f"sample {idx}: G_tr - G_flat = {g_sample}, F_flat - F_tr = {f_sample}")
    if g_sample <= 0:
        raise AssertionError(
            f"sample {idx} violates G_tr > G_flat: got {g_sample}"
        )
    if f_sample <= 0:
        raise AssertionError(
            f"sample {idx} violates F_tr < F_flat: got {f_sample}"
        )

banner("4. Endpoint bounds")
print("For 0 <= R <= 1:")
print("  G_flat <= G_tr <= xi")
print("  1/(1-xi) <= F_tr <= F_flat")

print("\nAll Stage-29 symbolic checks passed.")
