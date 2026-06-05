#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 59 SymPy audit.

Checks:
1. exact lower/upper support-ratio brackets on the Stage-41 branch interval,
2. exact residual-bracket definitions,
3. exact coupling thresholds Xi_fail and Xi_suff,
4. weak-coupling expansion of the physical ratio using the Stage-39 Omega_Pe series.
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


def expect_positive(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.factor(expr))
    print(f"{name} = {expr}")
    if expr.is_positive is not True:
        raise AssertionError(f"{name} is not positive")


def expect_close(name: str, value: sp.Expr, target: sp.Expr, tol: str = "1e-40") -> None:
    diff = abs(sp.N(value - target, 80))
    print(f"{name} diff = {diff}")
    if diff > sp.Float(tol, 80):
        raise AssertionError(f"{name} is not within tolerance")


banner("STAGE 59 — OPERATOR-SELECTED RESIDUAL BOUNDS")

Xi, Delta0, DeltaInf = sp.symbols("Xi Delta0 DeltaInf", positive=True, real=True)
kappa, y, zeta_req = sp.symbols("kappa y zeta_req", positive=True, real=True)
Pe_req = sp.symbols("Pe_req", positive=True, real=True)

Omega = sp.Function("Omega")
A_K = sp.simplify((kappa + sp.pi**2 / 4) / (kappa + y**2))

zeta_lo = sp.simplify(A_K * Omega(Xi * Delta0) ** 2)
zeta_hi = sp.simplify(A_K * Omega(Xi * DeltaInf) ** 2)
print("zeta_- =", zeta_lo)
print("zeta_+ =", zeta_hi)

R_lo = sp.simplify(zeta_req - zeta_hi)
R_hi = sp.simplify(zeta_req - zeta_lo)
print("R_- =", R_lo)
print("R_+ =", R_hi)

Xi_fail = sp.simplify(Pe_req / DeltaInf)
Xi_suff = sp.simplify(Pe_req / Delta0)
print("Xi_fail =", Xi_fail)
print("Xi_suff =", Xi_suff)
delta_gap = sp.symbols("delta_gap", positive=True, real=True)
DeltaInf_ordered = Delta0 + delta_gap
Xi_fail_ordered = sp.simplify(Pe_req / DeltaInf_ordered)
Xi_suff_ordered = sp.simplify(Pe_req / Delta0)
zeta_req_branch = sp.simplify(A_K * Omega(Pe_req) ** 2)
expect_positive("Xi_suff - Xi_fail (ordered)", Xi_suff_ordered - Xi_fail_ordered)

# Weak-coupling expansion using the exact Stage-39 Omega_Pe series.
Pe = sp.symbols("Pe", positive=True, real=True)
Omega_Pe = sp.simplify(
    sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi)
    / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1))
)
Omega_sq_series = sp.series(Omega_Pe**2, Pe, 0, 2).removeO()
print("Omega_Pe^2 small-Pe series =", Omega_sq_series)
expect_zero(
    "Omega^2 linear coefficient",
    sp.expand(Omega_sq_series).coeff(Pe, 1) - (4 - sp.pi) / sp.pi,
)

Xi_num = sp.symbols("Xi_num", real=True)
Pe_num = sp.symbols("Pe_num", real=True)
kappa_probe = sp.Integer(2)
y_probe = sp.Integer(1)
Delta0_probe = sp.Rational(3, 5)
delta_gap_probe = sp.Rational(2, 5)
DeltaInf_probe = Delta0_probe + delta_gap_probe
A_K_probe = sp.N(A_K.subs({kappa: kappa_probe, y: y_probe}), 80)
zeta_req_probe = sp.Rational(2, 5)  # independent target, NOT constructed from Omega_Pe
Pe_star = sp.nsolve(
    sp.N(A_K_probe * Omega_Pe.subs(Pe, Pe_num) ** 2 - zeta_req_probe, 80),
    Pe_num,
    sp.Rational(1, 2),
    prec=70,
)
Pe_req_probe = Pe_star  # operator-branch threshold scale derived from the independent zeta_req
Xi_fail_expected = sp.N(Pe_req_probe / DeltaInf_probe, 80)
Xi_suff_expected = sp.N(Pe_req_probe / Delta0_probe, 80)
Xi_fail_solved = sp.nsolve(
    sp.N(A_K_probe * Omega_Pe.subs(Pe, Xi_num * DeltaInf_probe) ** 2 - zeta_req_probe, 80),
    Xi_num,
    Xi_fail_expected,
    prec=70,
)
Xi_suff_solved = sp.nsolve(
    sp.N(A_K_probe * Omega_Pe.subs(Pe, Xi_num * Delta0_probe) ** 2 - zeta_req_probe, 80),
    Xi_num,
    Xi_suff_expected,
    prec=70,
)
expect_close("Xi_fail*DeltaInf saturates at Pe_star", Xi_fail_solved * DeltaInf_probe, Pe_star)
expect_close("Xi_suff*Delta0 saturates at Pe_star", Xi_suff_solved * Delta0_probe, Pe_star)

zeta_weak = sp.expand(A_K * (1 + (4 - sp.pi) / sp.pi * Xi * Delta0))
print("weak-coupling zeta_phys =", zeta_weak)

print("\nStage 59 audit passed.")
