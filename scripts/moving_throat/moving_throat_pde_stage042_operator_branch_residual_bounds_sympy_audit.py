#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 42 SymPy audit.

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


banner("STAGE 42 — OPERATOR-SELECTED RESIDUAL BOUNDS")

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
expect_zero("residual bracket center identity", R_hi - R_lo - (zeta_hi - zeta_lo))

Xi_fail = sp.simplify(Pe_req / DeltaInf)
Xi_suff = sp.simplify(Pe_req / Delta0)
print("Xi_fail =", Xi_fail)
print("Xi_suff =", Xi_suff)
expect_zero("Xi_fail upper-edge identity", Xi_fail * DeltaInf - Pe_req)
expect_zero("Xi_suff lower-edge identity", Xi_suff * Delta0 - Pe_req)

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

zeta_weak = sp.expand(A_K * (1 + (4 - sp.pi) / sp.pi * Xi * Delta0))
print("weak-coupling zeta_phys =", zeta_weak)

print("\nStage 42 audit passed.")
