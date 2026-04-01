#!/usr/bin/env python3
"""
SymPy audit for Stage 62.

Verifies the explicit Family-1 map from the selected quadrupole demand zeta_req
into the required constructive transport bias Pe_req.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_small(name: str, expr: sp.Expr, tol: float = 1e-12) -> None:
    val = sp.N(sp.simplify(expr), 50)
    print(f"{name} = {val}")
    if abs(complex(val)) > tol:
        raise AssertionError(f"{name} is not within tolerance {tol}")


banner("STAGE 62 — FAMILY-1 QUADRUPOLE-DEMAND / PE MAP")

Pe = sp.symbols('Pe', positive=True, real=True)
y = sp.symbols('y', real=True)

kappa_F1 = sp.Rational(12321, 5)
eta_F1 = sp.Integer(37)

y_F1 = sp.nsolve(y * sp.tan(y) - eta_F1, 1.53, tol=1e-34, maxsteps=100)
A_F1 = (kappa_F1 + sp.pi**2 / 4) / (kappa_F1 + y_F1**2)

print("y_F1 =", sp.N(y_F1, 30))
print("A_F1 =", sp.N(A_F1, 30))
print("kappa_F1 =", kappa_F1)
print("eta_F1   =", eta_F1)
print("Robin residual =", sp.N(y_F1 * sp.tan(y_F1) - eta_F1, 30))

Omega = sp.simplify(sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi) / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1)))
Omega0 = sp.simplify(sp.limit(Omega, Pe, 0, dir='+'))
Omega_inf = sp.simplify(sp.limit(Omega, Pe, sp.oo))

print("Omega(Pe) =", Omega)
print("Omega(0+) =", Omega0)
print("Omega(inf) =", Omega_inf)
expect_small("Omega(0+) - 1", Omega0 - 1)
expect_small("Omega(inf) - pi/2", Omega_inf - sp.pi / 2)

zeta_F1 = sp.simplify(A_F1 * Omega**2)
zeta0 = sp.simplify(sp.limit(zeta_F1, Pe, 0, dir='+'))
zeta_inf = sp.simplify(sp.limit(zeta_F1, Pe, sp.oo))

print("zeta_F1(Pe) = A_F1 Omega(Pe)^2")
print("zeta_F1(0+) =", sp.N(zeta0, 30))
print("zeta_F1(inf) =", sp.N(zeta_inf, 30))
expect_small("zeta_F1(0+) - A_F1", zeta0 - A_F1)
expect_small("zeta_F1(inf) - A_F1 pi^2/4", zeta_inf - A_F1 * sp.pi**2 / 4)

# Small-Pe expansion check.
zeta_series = sp.series(zeta_F1, Pe, 0, 2).removeO()
expected_series = sp.expand(A_F1 * (1 + ((4 - sp.pi) / sp.pi) * Pe))
print("zeta_F1(Pe) small-Pe =", sp.expand(zeta_series))
expect_small("small-Pe coefficient check", sp.expand(zeta_series - expected_series), tol=1e-12)

banner("FINAL LEDGER")
print("zeta_F1(Pe) = A_F1 Omega_Pe^2")
print("A_F1 =", sp.N(A_F1, 20))
print("zeta_F1(0) =", sp.N(zeta0, 20))
print("zeta_max^(F1) =", sp.N(zeta_inf, 20))
print("A unique constructive Pe_req exists iff A_F1 <= zeta_req <= zeta_max^(F1).")
