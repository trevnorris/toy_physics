#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage296_299_common import *

"""
Stage 298 — microscopic coherent-kernel slippage decomposition and exact
triangular normal form.

This stage pushes the physical-branch defect law down to the actual microscopic
coherent-kernel drifts, defines the branch-adapted slippages
(Sigma_tr, Sigma_nt, Sigma_eta), and proves the exact triangular normal form
for (Theta_1, Xi_1, R_1).
"""

banner("STAGE 298 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE NORMAL FORM")

# Microscopic weak-axisymmetric logarithmic drifts.
lam1, c1, gam1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
kapU, kapEta, kapW, mu1, tau1 = sp.symbols("kappa_U kappa_eta kappa_W mu_1 tau_1", real=True)

# Coherent reference-branch variables.
chi0, deltaU, epsW, eps, epsEta = sp.symbols("chi_0 delta_U epsilon_W epsilon epsilon_eta", positive=True, real=True)

subbanner("I. Direct physical-branch drifts from the microscopic variables")
SigmaZ = sp.simplify(2 * lam1 + mu1 - kapEta - 2 * kapW)
SigmaChi = sp.simplify(gam1 + c1 - kapU)
SigmaEta = sp.simplify(2 * c1 - kapU - kapEta)
SigmaEps = sp.simplify(2 * gam1 + 2 * lam1 - kapU - kapW)
SigmaDelta = sp.simplify(tau1 - kapU)

zetaZ = sp.simplify(SigmaZ + (kapW - mu1))  # not used later; just anchor relation
omegaW = sp.simplify(kapW - mu1)
chi1 = sp.simplify(chi0 * SigmaChi)
eta1 = sp.simplify(epsEta * SigmaEta)
epsW1 = sp.simplify(epsW * SigmaEps)
deltaU1 = sp.simplify(deltaU * SigmaDelta)

print("Sigma_Z ="); sp.pprint(SigmaZ)
print("Sigma_chi ="); sp.pprint(SigmaChi)
print("Sigma_eta ="); sp.pprint(SigmaEta)
print("Sigma_epsilon ="); sp.pprint(SigmaEps)
print("Sigma_delta ="); sp.pprint(SigmaDelta)

subbanner("II. Exact defect law in microscopic slippages")
Xi1 = sp.simplify(
    SigmaZ
    + 2 * chi0 * SigmaChi / (1 + chi0)
    + (2 * epsW / (1 - eps))
    * (
        (11 + 9 * deltaU) * SigmaEps / (11 * (1 + deltaU))
        - 2 * deltaU * SigmaDelta / (11 * (1 + deltaU) ** 2)
    )
)
R1 = sp.simplify(-epsEta * SigmaEta / (1 - epsEta) - Xi1)
print("Xi_1 =")
sp.pprint(Xi1)
print("R_1 =")
sp.pprint(R1)
expect_zero("R_1 + Xi_1 + epsilon_eta*Sigma_eta/(1-epsilon_eta)", sp.simplify(R1 + Xi1 + epsEta * SigmaEta / (1 - epsEta)))

subbanner("III. Tracking/nontracking split")
SigmaTr = sp.simplify((1 + chi0) * SigmaDelta + (1 + deltaU) * SigmaChi)
C_tr = sp.simplify(chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)))
A_tr = sp.simplify(2 * chi0 / ((1 + chi0) * (1 + deltaU)))
SigmaNt = sp.simplify(
    SigmaZ
    + (2 * epsW / (1 - eps)) * (11 + 9 * deltaU) * SigmaEps / (11 * (1 + deltaU))
    - (2 * chi0 / (1 + deltaU) + 4 * epsW * deltaU / (11 * (1 - eps) * (1 + deltaU) ** 2)) * SigmaDelta
)
Theta1 = sp.simplify(-C_tr * SigmaTr)
Xi1_tri = sp.simplify(A_tr * SigmaTr + SigmaNt)

print("Sigma_tr =")
sp.pprint(SigmaTr)
print("Sigma_nt =")
sp.pprint(SigmaNt)
print("Theta_1 =")
sp.pprint(Theta1)
print("Xi_1 (triangular form) =")
sp.pprint(Xi1_tri)

expect_zero("Xi_1 - (A_tr Sigma_tr + Sigma_nt)", sp.simplify(Xi1 - Xi1_tri))

subbanner("IV. Exact inverse reconstruction")
SigmaTr_inv = sp.simplify(-Theta1 / C_tr)
SigmaNt_inv = sp.simplify(Xi1 + A_tr * Theta1 / C_tr)
SigmaEta_inv = sp.simplify(-(1 - epsEta) * (R1 + Xi1) / epsEta)

print("Sigma_tr from Theta_1 =")
sp.pprint(SigmaTr_inv)
print("Sigma_nt from (Theta_1, Xi_1) =")
sp.pprint(SigmaNt_inv)
print("Sigma_eta from (R_1, Xi_1) =")
sp.pprint(SigmaEta_inv)

expect_zero("Sigma_tr inverse", sp.simplify(SigmaTr_inv - SigmaTr))
expect_zero("Sigma_nt inverse", sp.simplify(SigmaNt_inv - SigmaNt))
expect_zero("Sigma_eta inverse", sp.simplify(SigmaEta_inv - SigmaEta))

subbanner("V. Zero-defect theorem")
print("Theta_1 = Xi_1 = R_1 = 0  iff  Sigma_tr = Sigma_nt = Sigma_eta = 0")
