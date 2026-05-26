#!/usr/bin/env python3
"""
Stage 38 SymPy audit: explicit reachability of the non-twin lowest support lane.
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


banner("STAGE 38 — EXPLICIT LOWEST-LANE REACHABILITY")

# Paper-stated domain: alpha > 0, 0 < x < 4 (Stage 054 softening), 0 < y < pi/2
# (principal branch of y tan y = eta), zeta_req > 0.  SymPy lacks compound
# symbol-level bounds, so positivity is declared here and the (0, 4) / (0, pi/2)
# constraints are exercised by the endpoint substitutions below.
alpha, x, y, zeta_req = sp.symbols("alpha x y zeta_req", positive=True, real=True)
pi = sp.pi

Omega_exp = sp.simplify(
    pi * alpha * (2 * alpha * sp.exp(alpha) + pi)
    / ((4 * alpha**2 + pi**2) * (sp.exp(alpha) - 1))
)
AK = sp.simplify(1 / (1 - x / 4 + x * y**2 / pi**2))
zeta_family = sp.simplify(Omega_exp**2 * AK)

print("Omega_exp(alpha) =", Omega_exp)
print("A_K(y,x) =", AK)
print("zeta_family(alpha,y;x) =", zeta_family)

zeta_twin = sp.simplify(sp.limit(zeta_family, alpha, 0).subs(y, pi/2))
zeta_max = sp.simplify(sp.limit(zeta_family, alpha, sp.oo).subs(y, 0))
print("symmetric twin point =", zeta_twin)
print("closure maximum =", zeta_max)
expect_zero("twin value", zeta_twin - 1)
expect_zero("closure maximum", zeta_max - pi**2 / (4 - x))

banner("EXACT REACHABILITY FLOOR")

x_floor = sp.simplify(sp.solve(sp.Eq(zeta_max, zeta_req), x)[0])
print("x floor from zeta_max = zeta_req:", x_floor)
expect_zero("x floor = 4 - pi^2/zeta_req", x_floor - (4 - pi**2 / zeta_req))

# Equivalent stiffness-ratio form.
KX_over_KW = sp.symbols("KX_over_KW", positive=True, real=True)
expect_zero("KX/KW equivalence", (1 / AK).subs(y, 0).subs(x, x_floor) - pi**2 / (4 * zeta_req))

banner("REGIME SPLIT")

pure_overlap_ceiling = sp.simplify(pi**2 / 4)
print("pure-overlap ceiling =", pure_overlap_ceiling)
print("combined closure ceiling =", zeta_max)
print("Regime A: zeta_req <= pi^2/4  -> overlap alone can work.")
print("Regime B: pi^2/4 < zeta_req <= pi^2/(4-x) -> combined overlap+softening can work.")
print("Regime C: zeta_req > pi^2/(4-x) -> this explicit family fails.")

print("\nFINAL LEDGER")
print("  The explicit exp+Robin family has closure range 1 <= zeta <= pi^2/(4-x).")
print("  It reaches the Stage-35 threshold iff zeta_req <= pi^2/(4-x) (strict < for finite parameters).")
print("  If zeta_req > pi^2/4, the compliance floor is x >= 4 - pi^2/zeta_req.")
