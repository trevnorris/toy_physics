#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 156 — AXISYMMETRIC TRANSPORT AND COLLAPSE TO ONE STATIC LOADING MISMATCH")

eps = sp.symbols('epsilon', real=True)
lambda20, lambda21, lambda22 = sp.Integer(1), sp.Rational(1,2), -sp.Integer(1)
D0, N0 = sp.symbols('D0 N0', nonzero=True, real=True)
D01, D21, D41, N01 = sp.symbols('D01 D21 D41 N01', real=True)
K1, M1, B01, B21, Z01, Z21, N1 = sp.symbols('K1 M1 B01 B21 Z01 Z21 N1', real=True)
D01_expr = sp.simplify(K1 - B01 - Z01)
D21_expr = sp.simplify(-(M1 + B21 + Z21))
N01_expr = N1

subbanner("1. Actual grouped weak-axisymmetric operator slopes")
u2_1 = sp.simplify(-(D21 + D01/9)/D0)
u4_1 = sp.simplify(-(5*D01 + 18*D21 + 81*D41)/(81*D0))
Xi_load = sp.simplify(N01/N0 - D01/D0)

print("u2^(1) =", u2_1)
print("u4^(1) =", u4_1)
print("Xi_load =", Xi_load)

subbanner("2. Hidden-even transport on the canonical branch")
expect_zero("u4^(1) - (8/9) u2^(1) + (D41 - 2 D21/3 - D01/27)/D0", u4_1 - sp.Rational(8,9)*u2_1 + (D41 - 2*D21/3 - D01/27)/D0)

D21_even = -D01/9
D41_even = -D01/27
expect_zero("u2^(1) on even-preserving branch", u2_1.subs({D21: D21_even}))
expect_zero("u4^(1) on even-preserving+hidden-even branch", u4_1.subs({D21: D21_even, D41: D41_even}))

subbanner("3. Grouped defect pattern")
DQ20 = sp.simplify(eps*lambda20*Xi_load)
DQ21 = sp.simplify(eps*lambda21*Xi_load)
DQ22 = sp.simplify(eps*lambda22*Xi_load)
Xi_bar, Xi_a, Xi_b = grouped_triplet(DQ20, DQ21, DQ22)
expect_zero("weak-axisymmetric trace vanishes", Xi_bar)
expect_zero("axisymmetric fingerprint b - 3 a", Xi_b - 3*Xi_a)
print("Delta_Q^(20) =", DQ20)
print("Delta_Q^(21) =", DQ21)
print("Delta_Q^(22) =", DQ22)

subbanner("4. Microscopic decomposition of the static slopes")
print("D01 =", D01_expr)
print("D21 =", D21_expr)
print("N01 =", N01_expr)
u2_1_micro = sp.simplify(-(D21_expr + D01_expr/9)/D0)
P1_over_P0_micro = sp.simplify(N01_expr/N0 - D01_expr/D0)
expect_zero("u2^(1) microscopic form", u2_1.subs({D01:D01_expr, D21:D21_expr}) - u2_1_micro)
expect_zero("Xi_load microscopic form", Xi_load.subs({D01:D01_expr, N01:N01_expr}) - P1_over_P0_micro)

banner("STAGE 156 LEDGER")
print("1. The weak-axisymmetric physical slopes are exact transport coefficients of")
print("   the actual grouped response moments:")
print("      u2^(1) = -(D21 + D01/9)/D0,")
print("      P1/P0  = N01/N0 - D01/D0.")
print("2. Hidden-even consistency is")
print("      D41 = (2/3) D21 + D01/27.")
print("3. On the even-preserving branch,")
print("      D21 = -D01/9,   D41 = -D01/27,")
print("   so the conservative grouped response is transported by one static slope D01.")
print("4. The remaining linear grouped 2.5PN defect is one static loading mismatch")
print("      Xi_load = N01/N0 - D01/D0.")
print("5. Its grouped-lane pattern is fixed by the weak-axisymmetric signature")
print("      (20,21,22) ~ (1, 1/2, -1).")
