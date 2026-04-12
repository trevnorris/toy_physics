#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage296_299_common import *

"""
Stage 296 — actual moving-throat static-slope compiler.

This stage rewrites the physical weak-axisymmetric slopes u2^(1) and P1/P0 in the
sharpest static grouped-bundle language:
- one conservative static slope D01 / D0,
- one outgoing-transfer static slope N1 / N0,
- and, on the even-preserving branch, only their difference Xi_1.
"""

banner("STAGE 296 — ACTUAL MOVING-THROAT STATIC-SLOPE COMPILER")

# Common grouped-bundle symbols.
D0, D2, u2 = sp.symbols("D_0 D_2 u_2", positive=True, real=True)
K1, M1 = sp.symbols("K_1 M_1", real=True)
B01, B21, Z01, Z21 = sp.symbols("B0_1 B2_1 Z0_1 Z2_1", real=True)
N0, N1 = sp.symbols("N_0 N_1", positive=True, real=True)
kappa1, nuN = sp.symbols("kappa_1 bar_nu_N", real=True)

subbanner("I. Exact microscopic decomposition of the two physical slopes")
D01 = sp.simplify(K1 - B01 - Z01)
D21 = sp.simplify(-(M1 + B21 + Z21))
u2_1 = sp.simplify(-(D21 + u2 * D01) / D0)
P1_over_P0 = sp.simplify(N1 / N0 - D01 / D0)

print("D01 =")
sp.pprint(D01)
print("D21 =")
sp.pprint(D21)
print("u2^(1) =")
sp.pprint(u2_1)
print("P1/P0 =")
sp.pprint(P1_over_P0)

subbanner("II. Even-preserving branch")
# On the even-preserving branch, D21 = -D01/9.
u2_1_even = sp.simplify(u2_1.subs({D21: -D01 / 9, u2: sp.Rational(1, 9)}))
print("u2^(1) on the even-preserving branch =")
sp.pprint(u2_1_even)
expect_zero("u2^(1)_(even)", u2_1_even)

Xi_load_even = sp.simplify(P1_over_P0.subs({D01: kappa1 * D0, N1: nuN * N0}))
print("Xi_1 = P1/P0 on the even-preserving branch =")
sp.pprint(Xi_load_even)
expect_zero("Xi_1 - (bar_nu_N - kappa_1)", Xi_load_even - (nuN - kappa1))

subbanner("III. Portwise outgoing-transfer slope compiler")
w1, w2, w3 = sp.symbols("w_1 w_2 w_3", real=True)
nu1, nu2, nu3 = sp.symbols("nu_1 nu_2 nu_3", real=True)
w_norm = sp.simplify(w1 + w2 + w3)
nu_bar = sp.simplify((w1 * nu1 + w2 * nu2 + w3 * nu3) / w_norm)

print("bar_nu_N =")
sp.pprint(nu_bar)
print("If the outgoing weights are normalized, w1 + w2 + w3 = 1, then")
sp.pprint(sp.simplify(nu_bar.subs({w_norm: 1})))

# Recover the grouped Y20 weak-axisymmetric pattern explicitly.
eps = sp.symbols("eps", real=True)
lam20 = sp.Integer(1)
lam21 = sp.Rational(1, 2)
lam22 = -sp.Integer(1)
Xi20 = sp.simplify(eps * lam20 * (nuN - kappa1))
Xi21 = sp.simplify(eps * lam21 * (nuN - kappa1))
Xi22 = sp.simplify(eps * lam22 * (nuN - kappa1))
gXi = grouped_trace_anomaly(Xi20, Xi21, Xi22)
print("Xi^(20), Xi^(21), Xi^(22) =")
sp.pprint(sp.Matrix([Xi20, Xi21, Xi22]))
print("Grouped trace/anomaly:")
sp.pprint(sp.Matrix([gXi["bar"], gXi["a"], gXi["b"]]))
expect_zero("b_Xi - 3 a_Xi", gXi["b"] - 3 * gXi["a"])

subbanner("IV. Interpretation")
print("The actual physical weak-axisymmetric continuation point is now:")
print("  1. one conservative static slope D01/D0,")
print("  2. one outgoing-transfer static slope N1/N0,")
print("  3. and, on the even-preserving branch, only their difference")
print("       Xi_1 = P1/P0 = bar_nu_N - kappa_1.")
