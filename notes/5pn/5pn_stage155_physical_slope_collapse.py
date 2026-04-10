#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 155 — COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM TO u2^(1) AND P1")

eps = sp.symbols('epsilon', real=True)
lam = sp.symbols('lambda_A', real=True)
D0, D2, D4, N0 = sp.symbols('D0 D2 D4 N0', nonzero=True, real=True)
D01, D21, D41, N01 = sp.symbols('D01 D21 D41 N01', real=True)
P0 = N0/D0

u2 = sp.Rational(1, 9)
u4 = sp.Rational(4, 81)

subbanner("1. Weak-axisymmetric grouped physical response")
D0A = D0 + eps*lam*D01
D2A = D2 + eps*lam*D21
D4A = D4 + eps*lam*D41
N0A = N0 + eps*lam*N01

u2A = sp.expand(sp.series(-D2A/D0A, eps, 0, 2).removeO()).subs({D2: -D0/9})
u4A = sp.expand(sp.series((D2A**2 - D0A*D4A)/D0A**2, eps, 0, 2).removeO()).subs({D2: -D0/9, D4: -D0/27})
P0A = sp.expand(sp.series(N0A/D0A, eps, 0, 2).removeO())

u2_1 = sp.simplify(u2A.coeff(eps,1)/lam)
u4_1 = sp.simplify(u4A.coeff(eps,1)/lam)
P1 = sp.simplify(P0A.coeff(eps,1)/lam)

print("u2^(1) =", u2_1)
print("u4^(1) =", u4_1)
print("P1     =", P1)
print("P1/P0  =", sp.simplify(P1/P0))

subbanner("2. Collapse of (frakK_1, frakG_1)")
frakK1, frakG1 = sp.symbols('frakK1 frakG1', real=True)
mathcalK = sp.simplify(D21 + D01/9)
mathcalG = sp.simplify(N01 - P0*D01)
expect_zero("frakK_1 + D0 u2^(1)", mathcalK + D0*u2_1)
expect_zero("frakG_1 - D0 P1", mathcalG - D0*P1)

subbanner("3. Hidden-even consistency relation")
expect_zero("u4^(1) - (8/9) u2^(1) + (D41 - 2 D21/3 - D01/27)/D0", u4_1 - sp.Rational(8,9)*u2_1 + (D41 - 2*D21/3 - D01/27)/D0)

subbanner("4. Direct outlet amplitudes in physical variables")
sigma_star = sp.symbols('sigma_star', nonzero=True, real=True)
kappa1 = sp.simplify(-3*(1-sigma_star)*u2_1/sigma_star)
gamma1 = sp.simplify(-(1-sigma_star)*(P1/P0)/(9*sigma_star))
print("kappa_1 =", kappa1)
print("gamma_1 =", gamma1)

subbanner("5. Even-preserving branch")
Xi_load = sp.simplify(P1/P0)
D21_even = -D01/9
D41_even = -D01/27
expect_zero("u2^(1) on even-preserving branch", u2_1.subs({D21: D21_even}))
expect_zero("u4^(1) on even+hidden-even branch", u4_1.subs({D21: D21_even, D41: D41_even}))
print("Xi_load =", Xi_load)

banner("STAGE 155 LEDGER")
print("1. On the canonical compensated branch, the microscopic obstruction pair is")
print("      frakK_1 = -D0 u2^(1),")
print("      frakG_1 =  D0 P1.")
print("2. The one-parameter hidden-even consistency condition is exactly")
print("      u4^(1) = (8/9) u2^(1),")
print("   equivalently D41 = (2/3) D21 + D01/27.")
print("3. The direct outlet amplitudes are physical:")
print("      kappa_1 = -3(1-sigma_*) u2^(1)/sigma_*,")
print("      gamma_1 = -(1-sigma_*) (P1/P0)/(9 sigma_*).")
print("4. On the even-preserving branch the remaining linear grouped defect is just")
print("      Xi_load = P1/P0.")
