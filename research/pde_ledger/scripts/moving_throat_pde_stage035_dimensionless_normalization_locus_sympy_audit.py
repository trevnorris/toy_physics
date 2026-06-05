#!/usr/bin/env python3
"""
Moving-throat PDE Stage 35 SymPy audit.

Checks:
1. Exact D/N dimensionless shape function F(xi, delta).
2. Exact monotonicity of F on the stable branch.
3. Exact onset and softening limits.
4. Exact required total loading alpha_req(xi, delta).
5. Exact required support-loading formula and near-onset / near-softening asymptotics.
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


# Exact D/N overlap constants.
kappa0_sq = sp.Rational(8) / sp.pi**2
kappa1_sq = sp.Rational(16) / (9 * sp.pi**2)
eta = sp.simplify(kappa1_sq / kappa0_sq)

# Dimensionless variables.
A, DeltaK, beta0, NQ = sp.symbols("A DeltaK beta0 NQ", positive=True, real=True)
xi = sp.symbols("xi", nonnegative=True, real=True)
delta = sp.symbols("delta", positive=True, real=True)
Chi, OmegaU, Delta0 = sp.symbols("Chi Omega_U Delta_0", positive=True, real=True)

x = sp.simplify(A * xi)
DeltaK_sub = sp.simplify(A * delta)

alpha_x = sp.simplify(x * (x + DeltaK_sub) / (kappa0_sq * (x + DeltaK_sub) + kappa1_sq * x))
N_x = sp.simplify(
    beta0 * (kappa0_sq * (x + DeltaK_sub) + kappa1_sq * x) ** 4
    / (kappa0_sq * (A - x) * (kappa0_sq * (x + DeltaK_sub) ** 2 + kappa1_sq * x**2) ** 2)
)

banner("STAGE 35.1 — EXACT D/N DIMENSIONLESS SHAPE FUNCTION")
F = sp.simplify(N_x / (beta0 * kappa0_sq / A))
F_target = sp.simplify((9 * delta + 11 * xi) ** 4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2))
print("eta =", eta)
print("F(xi,delta) =")
sp.pprint(F)
expect_zero("F - closed D/N form", F - F_target)

R_target = sp.simplify(NQ * A / (beta0 * kappa0_sq))
print("R_target =", R_target)
print("Normalization locus: F(xi,delta) = R_target")

banner("STAGE 35.2 — EXACT MONOTONICITY, ONSET, AND SOFTENING LIMITS")
dF = sp.simplify(sp.diff(F_target, xi))
dF_target = sp.simplify(
    (9 * delta + 11 * xi) ** 3
    * (81 * delta**3 + 72 * delta**2 + 189 * delta**2 * xi + 297 * delta * xi**2 + 121 * xi**3)
    / (81 * (1 - xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 3)
)
print("dF/dxi =")
sp.pprint(dF)
expect_zero("dF/dxi - manifestly positive form", dF - dF_target)
expect_zero("F(0,delta) - 1", sp.simplify(F_target.subs(xi, 0) - 1))
soft_const = sp.simplify(sp.limit((1 - xi) * F_target, xi, 1, dir="-"))
soft_const_target = sp.simplify((9 * delta + 11) ** 4 / (81 * (9 * delta**2 + 18 * delta + 11) ** 2))
print("softening constant C(delta) =", soft_const)
expect_zero("softening constant - closed form", soft_const - soft_const_target)
eps_soft = sp.symbols("eps_soft", positive=True, real=True)
soft_scaled_series = sp.series(eps_soft * F_target.subs(xi, 1 - eps_soft), eps_soft, 0, 1).removeO()
print("eps_soft * F(1-eps_soft, delta) near eps_soft=0 =", soft_scaled_series)
expect_zero("near-softening asymptotic coefficient", soft_scaled_series - soft_const_target)

banner("STAGE 35.3 — EXACT REQUIRED TOTAL LOADING AND SUPPORT COUPLING")
alpha_req = sp.simplify(alpha_x)
alpha_req_target = sp.simplify(9 * sp.pi**2 * A * xi * (xi + delta) / (8 * (9 * delta + 11 * xi)))
print("alpha_req(xi,delta) =")
sp.pprint(alpha_req)
expect_zero("alpha_req - closed D/N form", alpha_req - alpha_req_target)

alpha_crit = sp.simplify(sp.limit(alpha_req_target, xi, 1, dir="-"))
alpha_crit_target = sp.simplify(9 * sp.pi**2 * A * (1 + delta) / (8 * (11 + 9 * delta)))
print("alpha_crit =", alpha_crit)
expect_zero("alpha_crit - closed form", alpha_crit - alpha_crit_target)

alpha_mix = sp.simplify(Chi**2 / (OmegaU**2 * Delta0))
gBreq_sq_over_varpi2 = sp.simplify(alpha_req_target - alpha_mix)
print("g_B,req^2 / varpi^2 =")
sp.pprint(gBreq_sq_over_varpi2)

banner("STAGE 35.4 — EXACT ASYMPTOTICS OF THE NORMALIZATION LOCUS")
F_series = sp.series(F_target, xi, 0, 3).removeO()
F_series_target = sp.simplify(1 + xi * (1 + sp.Rational(8, 9) / delta) + xi**2 * (1 + sp.Rational(8, 9) / delta - sp.Rational(28, 27) / delta**2))
print("F(xi,delta) near xi=0 =", F_series)
expect_zero("near-onset series through O(xi^2)", F_series - F_series_target)

alpha_series = sp.series(alpha_req_target, xi, 0, 3).removeO()
alpha_series_target = sp.simplify(sp.pi**2 * A * xi / 8 - sp.pi**2 * A * xi**2 / (36 * delta))
print("alpha_req near xi=0 =", alpha_series)
expect_zero("alpha_req near-onset series through O(xi^2)", alpha_series - alpha_series_target)

print("All Stage 35 checks passed.")
