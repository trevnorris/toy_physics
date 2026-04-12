#!/usr/bin/env python3
"""
5pn_stage25_dimensionless_normalization_and_support_frontier.py

Twenty-fifth executable SymPy audit for the 5PN grouped-real-P2 / moving-throat
program.

What this script does
---------------------
1. Reduces the selected-branch normalization problem to the exact D/N
   dimensionless shape function F(xi,delta).
2. Proves strict monotonicity of F and its endpoint values.
3. Derives the exact support-feasibility function G(xi,delta) and its endpoint
   geometry.
4. Records the near-onset asymptotics of both branch functions.

Interpretation
--------------
The selected quadrupole branch is controlled by two universal scalar functions,
F and G, of the same branch coordinate xi.
"""

from __future__ import annotations

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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("I. EXACT D/N DIMENSIONLESS SHAPE FUNCTION F(xi,delta)")

# Exact D/N constants.
kappa0 = 2 * sp.sqrt(2) / sp.pi
kappa1 = -4 / (3 * sp.pi)
eta = sp.simplify(kappa1**2 / kappa0**2)
print("eta = kappa_1^2 / kappa_0^2 =", eta)
expect_zero("eta - 2/9", eta - sp.Rational(2, 9))

# Dimensionless branch variables.
xi, delta = sp.symbols("xi delta", positive=True, real=True)
A, beta0 = sp.symbols("A beta_0", positive=True, real=True)

x = sp.simplify(A * xi)
DeltaK_ax = sp.simplify(A * delta)

# Exact Stage-24 formulas rewritten in dimensionless form.
s_x = sp.simplify(
    (kappa0**2 * (x + DeltaK_ax) + kappa1**2 * x) ** 2
    / (kappa0**2 * (x + DeltaK_ax) ** 2 + kappa1**2 * x**2)
)
N_minus = sp.simplify(beta0 * s_x**2 / (kappa0**2 * (A - x)))
N_onset = sp.simplify(beta0 * kappa0**2 / A)
F = sp.simplify(N_minus / N_onset)

F_expected = sp.simplify(
    (9 * delta + 11 * xi) ** 4
    / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
)

print("F(xi,delta) =")
sp.pprint(F)
expect_zero("F - expected closed form", F - F_expected)

subbanner("I.1 — Exact monotonicity and endpoint values")
dF_dxi = sp.simplify(sp.diff(F_expected, xi))
dF_expected = sp.simplify(
    (9 * delta + 11 * xi) ** 3
    * (81 * delta**3 + 72 * delta**2 + 189 * delta**2 * xi + 297 * delta * xi**2 + 121 * xi**3)
    / (81 * (1 - xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 3)
)
expect_zero("dF/dxi minus expected positive form", dF_dxi - dF_expected)

print("F(0,delta) =", sp.simplify(F_expected.subs(xi, 0)))
print("lim_{xi->1^-} F(xi,delta) =", sp.limit(F_expected, xi, 1, dir='-'))

C_delta = sp.simplify((9 * delta + 11) ** 4 / (81 * (9 * delta**2 + 18 * delta + 11) ** 2))
print("Near-softening coefficient C(delta) =")
sp.pprint(C_delta)
expect_zero(
    "(1-xi) F(xi,delta) near xi->1 equals C(delta)",
    sp.simplify(sp.limit((1 - xi) * F_expected, xi, 1, dir='-') - C_delta),
)

banner("II. EXACT SUPPORT-FEASIBILITY FUNCTION G(xi,delta)")

G = sp.simplify(9 * xi * (xi + delta) / (9 * delta + 11 * xi))
print("G(xi,delta) =")
sp.pprint(G)

# Connection to required total loading.
alpha_req = sp.simplify(sp.pi**2 * A * G / 8)
print("alpha_req(xi,delta) =")
sp.pprint(alpha_req)

dG_dxi = sp.simplify(sp.diff(G, xi))
dG_expected = sp.simplify(9 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) / (9 * delta + 11 * xi) ** 2)
expect_zero("dG/dxi minus expected positive form", dG_dxi - dG_expected)

print("G(0,delta) =", sp.simplify(G.subs(xi, 0)))
G_max = sp.simplify(sp.limit(G, xi, 1, dir='-'))
print("G_max(delta) =")
sp.pprint(G_max)
expect_zero("G_max - 9(1+delta)/(9 delta + 11)", G_max - sp.simplify(9 * (1 + delta) / (9 * delta + 11)))

banner("III. COMBINED ADMISSIBILITY GEOMETRY")

R_target, M_mix = sp.symbols("R_target M_mix", positive=True, real=True)
print("Exact admissibility frontier for fixed delta:")
print("  R_target = F(xi,delta)")
print("  M_crit   = G(xi,delta)")
print("with xi in [0,1).")

# Near-onset asymptotics.
subbanner("III.1 — Near-onset asymptotics")
F_series = sp.series(F_expected, xi, 0, 3).removeO()
G_series = sp.series(G, xi, 0, 3).removeO()
print("F(xi,delta) through O(xi^2) =")
sp.pprint(sp.expand(F_series))
print("G(xi,delta) through O(xi^2) =")
sp.pprint(sp.expand(G_series))

F_series_expected = sp.simplify(1 + (1 + 8 / (9 * delta)) * xi + (1 + 8 / (9 * delta) - 28 / (27 * delta**2)) * xi**2)
G_series_expected = sp.simplify(xi - 2 * xi**2 / (9 * delta))
expect_zero("F near-onset series minus expected", sp.expand(F_series - F_series_expected))
expect_zero("G near-onset series minus expected", sp.expand(G_series - G_series_expected))

eps_R = sp.symbols("eps_R", real=True)
xi_req_near_onset = sp.simplify(eps_R / (1 + 8 / (9 * delta)))
print("If R_target = 1 + eps_R with eps_R << 1, then")
print("xi_req ~=", xi_req_near_onset)
print("M_crit ~=", xi_req_near_onset)

banner("FINAL LEDGER")
print("1. The exact normalization locus is F(xi,delta) = R_target.")
print("2. F is strictly increasing from 1 at onset to +infinity at softening.")
print("3. The exact support-feasibility frontier is G(xi,delta), also strictly increasing.")
print("4. For fixed delta, the physical defect must lie on the unique F-locus and below the G-frontier.")
