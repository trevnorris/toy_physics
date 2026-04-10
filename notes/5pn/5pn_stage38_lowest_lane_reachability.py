#!/usr/bin/env python3
"""
5pn_stage38_lowest_lane_reachability.py

SymPy audit for Moving-Throat PDE Stage 38.
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


banner("STAGE 38 — EXPLICIT REACHABILITY OF THE NON-TWIN LOWEST SUPPORT LANE")

# Symbols
alpha, x = sp.symbols("alpha x", positive=True, real=True)
y, zeta_req = sp.symbols("y zeta_req", positive=True, real=True)
pi = sp.pi

Omega_exp = sp.simplify(
    pi * alpha * (2 * alpha * sp.exp(alpha) + pi)
    / ((4 * alpha**2 + pi**2) * (sp.exp(alpha) - 1))
)
A_K = sp.simplify(1 / (1 - x / 4 + x * y**2 / pi**2))
zeta_expR = sp.simplify(Omega_exp**2 * A_K)

print("zeta_0^(exp+R)(alpha,eta) =")
sp.pprint(zeta_expR)

expect_zero("symmetric twin point", sp.simplify(sp.limit(zeta_expR.subs(y, pi / 2), alpha, 0) - 1))
expect_zero(
    "closure ceiling",
    sp.simplify(sp.limit(zeta_expR.subs(y, 0), alpha, sp.oo) - pi**2 / (4 - x)),
)

reachability_floor = sp.simplify(4 - pi**2 / zeta_req)
KX_over_KW = sp.symbols("KX_over_KW", positive=True, real=True)
expect_zero(
    "stiffness-ratio form of the reachability floor",
    sp.simplify((1 - x / 4).subs(x, reachability_floor) - pi**2 / (4 * zeta_req)),
)

print("\nExact Stage-35 reachability criterion for this explicit family:")
print("  zeta_req <= pi^2 / (4 - x)")
print("equivalently")
print("  x >= 4 - pi^2 / zeta_req")

banner("STAGE 38 THEOREM LEDGER")
print("Three exact sub-regimes follow.")
print("  A: zeta_req <= pi^2/4")
print("     overlap enhancement alone can reach the threshold.")
print("  B: pi^2/4 < zeta_req <= pi^2/(4-x)")
print("     the explicit combined family reaches the threshold using both overlap and softening.")
print("  C: zeta_req > pi^2/(4-x)")
print("     even the explicit combined family fails; a stronger operator deformation is required.")
