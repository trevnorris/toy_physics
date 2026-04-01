#!/usr/bin/env python3
"""
Moving-throat PDE Stage 19 SymPy audit.

Checks:
1. Exact dimensionless support-feasibility function G(xi,delta).
2. Exact monotonicity and endpoint values of G.
3. Exact factorization of the required support loading through G - M_mix.
4. Near-onset asymptotics of the support-feasibility frontier.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr):
    simp = sp.simplify(sp.expand(expr))
    print(f"{name} = {simp}")
    if simp != 0:
        raise AssertionError(f"{name} is not zero")


A, delta = sp.symbols("A delta", positive=True, real=True)
xi = sp.symbols("xi", nonnegative=True, real=True)
Chi, OmegaU, Delta0 = sp.symbols("Chi Omega_U Delta_0", positive=True, real=True)

# Dimensionless support-feasibility variables.
Mmix = sp.symbols("M_mix", positive=True, real=True)
G = sp.simplify(9 * xi * (xi + delta) / (9 * delta + 11 * xi))
alpha_mix = sp.simplify(Chi**2 / (OmegaU**2 * Delta0))
Mmix_expr = sp.simplify(8 * alpha_mix / (sp.pi**2 * A))

a_req = sp.simplify(9 * sp.pi**2 * A * xi * (xi + delta) / (8 * (9 * delta + 11 * xi)))
gBreq_sq_over_varpi2 = sp.simplify(a_req - alpha_mix)

banner("STAGE 19.1 — EXACT SUPPORT-FEASIBILITY FUNCTION")
print("G(xi,delta) =", G)
print("M_mix =", Mmix_expr)
expect_zero(
    "g_B,req^2/varpi^2 - (pi^2 A / 8) (G - M_mix)",
    gBreq_sq_over_varpi2 - (sp.pi**2 * A / 8) * (G - Mmix_expr),
)

banner("STAGE 19.2 — EXACT MONOTONICITY AND ENDPOINTS OF G")
dG = sp.simplify(sp.diff(G, xi))
dG_target = sp.simplify(9 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) / (9 * delta + 11 * xi) ** 2)
print("dG/dxi =", dG)
expect_zero("dG/dxi - manifestly positive form", dG - dG_target)
expect_zero("G(0,delta)", sp.simplify(G.subs(xi, 0)))
Gmax = sp.simplify(sp.limit(G, xi, 1, dir='-'))
Gmax_target = sp.simplify(9 * (1 + delta) / (9 * delta + 11))
print("G_max(delta) =", Gmax)
expect_zero("G_max - closed form", Gmax - Gmax_target)

banner("STAGE 19.3 — NEAR-ONSET ASYMPTOTICS")
G_series = sp.series(G, xi, 0, 3).removeO()
G_series_target = sp.simplify(xi - 2 * xi**2 / (9 * delta))
print("G(xi,delta) near xi=0 =", G_series)
expect_zero("G near-onset series through O(xi^2)", G_series - G_series_target)

print("All Stage 19 checks passed.")
