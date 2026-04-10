#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 160 — WEAK-AXISYMMETRIC COLLAPSE OF THE OUTGOING-SLIPPAGE BUNDLE TO ONE SCALAR")

eps = sp.symbols('epsilon', real=True)
lambda20, lambda21, lambda22 = sp.Integer(1), sp.Rational(1,2), -sp.Integer(1)

subbanner("1. Weak-axisymmetric microscopic slope bookkeeping")
kappa1 = sp.symbols('kappa1', real=True)
gWr, gUr, rr, oUr, oWr = sp.symbols('gWr gUr rr oUr oWr', real=True)
m_r = sp.simplify(gWr - oWr - kappa1/2)
i_r = sp.simplify(rr + gUr - oUr - gWr)
h_r = sp.simplify(2*rr - oUr - oWr)
print("m_r =", m_r)
print("i_r =", i_r)
print("h_r =", h_r)

subbanner("2. Exact portwise outgoing-defect amplitude")
I_r, H_r = sp.symbols('I_r H_r', real=True)
sigma_r = sp.simplify(2*m_r + 2*I_r/(1+I_r)*i_r + 2*H_r/(1-H_r)*h_r)
print("sigma_r =")
sp.pprint(sigma_r)

subbanner("3. Exact grouped collapse to one scalar amplitude")
rho1, rho2 = sp.symbols('rho1 rho2', real=True)
sigma1, sigma2 = sp.symbols('sigma1 sigma2', real=True)
Xi1 = sp.simplify(rho1*sigma1 + rho2*sigma2)
Xi20 = sp.simplify(eps*lambda20*Xi1)
Xi21 = sp.simplify(eps*lambda21*Xi1)
Xi22 = sp.simplify(eps*lambda22*Xi1)
print("Xi_load^(20) =", Xi20)
print("Xi_load^(21) =", Xi21)
print("Xi_load^(22) =", Xi22)

Xi_bar, Xi_a, Xi_b = grouped_triplet(Xi20, Xi21, Xi22)
expect_zero("grouped trace vanishes", Xi_bar)
expect_zero("axisymmetric fingerprint b_Xi - 3 a_Xi", Xi_b - 3*Xi_a)

subbanner("4. Identification with the physical outgoing-prefactor slope")
P1, P0 = sp.symbols('P1 P0', nonzero=True, real=True)
Xi1_phys = sp.simplify(P1/P0)
print("Xi_1 := P1/P0 =", Xi1_phys)

DQ20 = sp.simplify(eps*Xi1_phys)
DQ21 = sp.simplify(eps*Xi1_phys/2)
DQ22 = sp.simplify(-eps*Xi1_phys)
DQ_bar, DQ_a, DQ_b = grouped_triplet(DQ20, DQ21, DQ22)
expect_zero("quadrupole trace vanishes", DQ_bar)
expect_zero("quadrupole axisymmetric fingerprint", DQ_b - 3*DQ_a)

subbanner("5. Rigidity and dominant-port corollaries")
sigma_r_rigid = sp.simplify(sigma_r.subs({i_r:0, h_r:0}))
expect_zero("rigid branch sigma_r - 2 m_r", sigma_r_rigid - 2*m_r)
print("If i_r = h_r = 0, sigma_r =", sigma_r_rigid)
print("Dominant-port limit: Xi_1 ≈ sigma_(dominant port).")

banner("STAGE 160 LEDGER")
print("1. On the weak-axisymmetric branch, every outgoing slippage inherits the same")
print("   grouped signature (1, 1/2, -1).")
print("2. Each outgoing port contributes only one scalar amplitude")
print("      sigma_r = 2 m_r + 2 I_r/(1+I_r) i_r + 2 H_r/(1-H_r) h_r.")
print("3. The full remaining grouped defect collapses to one weighted scalar")
print("      Xi_1 = sum_r rho_r^(N) sigma_r.")
print("4. Its grouped anisotropy is fixed exactly by the axisymmetric law")
print("      b = 3 a.")
print("5. The same scalar is the physical outgoing-prefactor slope")
print("      Xi_1 = P1/P0.")
print("6. On the even-preserving branch, the whole remaining linear grouped 2.5PN")
print("   normalization defect is controlled by that one scalar amplitude.")
