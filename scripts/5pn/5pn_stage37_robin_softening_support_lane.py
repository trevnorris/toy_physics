#!/usr/bin/env python3
"""
5pn_stage37_robin_softening_support_lane.py

SymPy audit for Moving-Throat PDE Stage 37.
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


banner("STAGE 37 — ROBIN-COMPLIANCE SOFTENING OF THE LOWEST SUPPORT LANE")

# Symbols
y, eta = sp.symbols("y eta", positive=True, real=True)
x, zeta_req = sp.symbols("x zeta_req", positive=True, real=True)
pi = sp.pi

print("Robin eigenvalue equation:")
print("  y tan y = eta")

A_K = sp.simplify(1 / (1 - x / 4 + x * y**2 / pi**2))
print("\nA_K(eta) =", A_K)

expect_zero("A_K(y=pi/2) - 1", sp.simplify(A_K.subs(y, pi / 2) - 1))
expect_zero("A_K(y=0) - 4/(4-x)", sp.simplify(sp.limit(A_K, y, 0) - 4 / (4 - x)))

y_req_sq = sp.simplify((pi**2 / x) * (1 / zeta_req - 1 + x / 4))
print("\ny_req^2 =", y_req_sq)
expect_zero("A_K(y_req) - zeta_req", sp.simplify(A_K.subs(y**2, y_req_sq) - zeta_req))

eta_req = sp.simplify(sp.sqrt(y_req_sq) * sp.tan(sp.sqrt(y_req_sq)))
print("eta_req =", eta_req)

print("\nPure-softening rescue exists only if")
print("  zeta_req <= 4/(4-x)")
print("equivalently")
print("  x >= 4 - 4/zeta_req")

banner("STAGE 37 THEOREM LEDGER")
print("The Robin-deformed lowest support branch is determined by")
print("  y tan y = eta,   0 < y < pi/2.")
print("Its exact softening factor is")
print("  A_K = 1 / [ 1 - x/4 + x y^2 / pi^2 ].")
print("The exact support-softening window is")
print("  1 <= A_K <= 4/(4-x),")
print("with the upper endpoint reached only in the soft-mouth closure eta -> 0^+.")
