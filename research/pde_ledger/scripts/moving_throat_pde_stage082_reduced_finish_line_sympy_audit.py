#!/usr/bin/env python3
"""
Stage 82 SymPy audit: the reduced finish line is a single normalization defect.
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


banner("STAGE 82 — REDUCED FINISH LINE")

G, c, c_s, a, Omega_Q = sp.symbols("G c c_s a Omega_Q", positive=True, real=True)
omega = sp.symbols("omega", real=True)
N_Q = sp.symbols("N_Q", positive=True, real=True)

Yhat_Q_cons = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / Omega_Q**2))
print("Yhat_Q^cons(omega) =", Yhat_Q_cons)

K0_target = sp.simplify(64 * G * Omega_Q**5 / (45 * c**5))
K0_target_geom = sp.simplify(K0_target.subs(Omega_Q, 3 * c_s / (2 * a)))
expect_zero("K0_target geometric form", K0_target_geom - 54 * G * c_s**5 / (5 * a**5 * c**5))

K0 = sp.simplify(N_Q * K0_target)
K2 = sp.simplify(K0 / (4 * Omega_Q**2))
K4 = sp.simplify(K0 / (4 * Omega_Q**4))
Gamma5 = sp.simplify(9 * K2 ** sp.Rational(5, 2) / K0 ** sp.Rational(3, 2))

K2_target = sp.simplify(K0_target / (4 * Omega_Q**2))
K4_target = sp.simplify(K0_target / (4 * Omega_Q**4))
Gamma5_target = sp.simplify(9 * K2_target ** sp.Rational(5, 2) / K0_target ** sp.Rational(3, 2))

R0 = sp.simplify(K0 / K0_target - 1)
R2 = sp.simplify(K2 / K2_target - 1)
R4 = sp.simplify(K4 / K4_target - 1)
R5 = sp.simplify(Gamma5 / Gamma5_target - 1)

expect_zero("R0 - (N_Q - 1)", R0 - (N_Q - 1))
expect_zero("R2 - (N_Q - 1)", R2 - (N_Q - 1))
expect_zero("R4 - (N_Q - 1)", R4 - (N_Q - 1))
expect_zero("R5 - (N_Q - 1)", R5 - (N_Q - 1))

print("R0 =", sp.factor(R0))
print("R2 =", sp.factor(R2))
print("R4 =", sp.factor(R4))
print("R5 =", sp.factor(R5))
print("\nSTAGE 82 AUDIT PASSED")
