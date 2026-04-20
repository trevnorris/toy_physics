#!/usr/bin/env python3
"""
moving_throat_pde_stage156_axisymmetric_loading_mismatch_sympy_audit.py

SymPy-backed audit for the weak-axisymmetric transport of the grouped physical
slopes and the collapse to one static loading mismatch.

Checks:
1. Exact weak-axisymmetric expansions for u2, u4, and P0.
2. Hidden-even operator identity on the canonical branch.
3. Even-preserving collapse D21 = -D01/9, D41 = -D01/27.
4. Remaining defect Xi_load = N01/N0 - D01/D0.
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

banner("STAGE 156 — WEAK-AXISYMMETRIC PHYSICAL-SLOPE TRANSPORT")

eps, lam = sp.symbols("eps lam", real=True)
D0, D01, D2, D21, D4, D41 = sp.symbols("D0 D01 D2 D21 D4 D41", real=True, nonzero=True)
N0, N01 = sp.symbols("N0 N01", real=True, nonzero=True)
u2, u4, P0 = sp.symbols("u2 u4 P0", real=True)

D0A = D0 + eps*lam*D01
D2A = D2 + eps*lam*D21
D4A = D4 + eps*lam*D41
N0A = N0 + eps*lam*N01

u2A = sp.expand(sp.series(-D2A/D0A, eps, 0, 2).removeO())
u4A = sp.expand(sp.series((D2A**2 - D0A*D4A)/D0A**2, eps, 0, 2).removeO())
P0A = sp.expand(sp.series(N0A/D0A, eps, 0, 2).removeO())

u21 = sp.simplify(sp.diff(u2A, eps).subs(eps, 0)/lam)
u41 = sp.simplify(sp.diff(u4A, eps).subs(eps, 0)/lam)
P1 = sp.simplify(sp.diff(P0A, eps).subs(eps, 0)/lam)

print("u2^(1) general =", u21)
print("u4^(1) general =", u41)
print("P1 general     =", P1)

expect_zero(
    "u2 slope identity",
    u21 - (-(D21 + u2*D01)/D0).subs({u2: -D2/D0})
)

banner("Canonical branch formulas")
u21_can = sp.simplify(u21.subs(D2, -sp.Rational(1,9)*D0))
u41_can = sp.simplify(u41.subs({D2: -sp.Rational(1,9)*D0, D4: -sp.Rational(1,27)*D0}))
P1_ratio = sp.simplify(P1 / (N0/D0))

print("u2^(1) canonical =", u21_can)
print("u4^(1) canonical =", u41_can)
print("P1/P0            =", P1_ratio)

expect_zero(
    "u4 canonical formula",
    u41_can + (5*D01 + 18*D21 + 81*D41)/(81*D0)
)
expect_zero(
    "P1/P0 formula",
    P1_ratio - (N01/N0 - D01/D0)
)

banner("Hidden-even operator identity")
hidden_even_residual = sp.expand(
    u41_can - sp.Rational(8,9)*u21_can
    - (D01/(27*D0) + 2*D21/(3*D0) - D41/D0)
)
expect_zero("hidden-even residual", hidden_even_residual)

D41_hidden = sp.simplify(sp.solve(sp.Eq(u41_can, sp.Rational(8,9)*u21_can), D41)[0])
print("D41 from hidden-even relation =", D41_hidden)

banner("Even-preserving collapse")
u21_zero_D21 = sp.simplify(sp.solve(sp.Eq(u21_can, 0), D21)[0])
print("D21 from u2^(1)=0 =", u21_zero_D21)

D41_even = sp.simplify(D41_hidden.subs(D21, u21_zero_D21))
print("D41 on even-preserving branch =", D41_even)

expect_zero("D21 + D01/9", u21_zero_D21 + D01/9)
expect_zero("D41 + D01/27", D41_even + D01/27)

banner("Single static loading mismatch")
Xi_load = sp.simplify(N01/N0 - D01/D0)
print("Xi_load =", Xi_load)

lam20 = sp.Integer(1)
lam21 = sp.Rational(1,2)
lam22 = sp.Integer(-1)

Delta20 = sp.simplify(lam20 * Xi_load)
Delta21 = sp.simplify(lam21 * Xi_load)
Delta22 = sp.simplify(lam22 * Xi_load)

print("Delta_Q^(20)/eps =", Delta20)
print("Delta_Q^(21)/eps =", Delta21)
print("Delta_Q^(22)/eps =", Delta22)

print("\nCarry-forward formulas:")
print("  u2^(1) = -(D21 + u2 D01)/D0")
print("  u4^(1) = -(5 D01 + 18 D21 + 81 D41)/(81 D0) on the canonical branch")
print("  P1/P0  = N01/N0 - D01/D0")
print("  hidden-even  <=>  D41 = 2 D21/3 + D01/27")
print("  if u2^(1)=0, then D21 = -D01/9 and D41 = -D01/27")
print("  remaining defect Xi_load = N01/N0 - D01/D0")
