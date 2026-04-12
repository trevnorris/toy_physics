#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage300_302_common import *

"""
Stage 302 — exact first-order defect compiler from quotient residuals.

This stage takes the three quotient residual coordinates of Stage 300 / 301 and
compiles them into the physical first-order defect triplet
    (Theta_1, Xi_1, R_1),
showing that the actual moving-throat branch can be tested either in quotient
coordinates or in direct physical defect language with an exact invertible bridge.
"""

banner("STAGE 302 — SELECTED BRANCH DEFECT COMPILER")

chi0s, deltaUs, epsEtas = sp.symbols("chi0_star deltaU_star epsilon_eta_star", positive=True, real=True)
qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

subbanner("I. Exact first-order map from quotient residuals to physical defects")
C_tr = sp.simplify(chi0s * deltaUs / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs)))
A_tr = sp.simplify(2 * chi0s / ((1 + chi0s) * (1 + deltaUs)))

Theta1 = sp.simplify(-C_tr * qtr)
Xi1 = sp.simplify(A_tr * qtr + qnt)
R1 = sp.simplify(-epsEtas * qeta / (1 - epsEtas) - Xi1)

defect_matrix = sp.Matrix([
    [-C_tr, 0, 0],
    [A_tr, 1, 0],
    [-A_tr, -1, -epsEtas / (1 - epsEtas)],
])
qvec = sp.Matrix([qtr, qnt, qeta])
defects = sp.Matrix([Theta1, Xi1, R1])

print("Defect matrix D_q =")
sp.pprint(defect_matrix)
print("(Theta_1, Xi_1, R_1)^T = D_q (q_tr, q_nt, q_eta)^T =")
sp.pprint(defects)
expect_zero("D_q q - defects", sp.simplify(defect_matrix * qvec - defects))

subbanner("II. Exact inverse map")
qtr_inv = sp.simplify(-Theta1 / C_tr)
qnt_inv = sp.simplify(Xi1 + A_tr * Theta1 / C_tr)
qeta_inv = sp.simplify(-(1 - epsEtas) * (R1 + Xi1) / epsEtas)

print("q_tr =")
sp.pprint(qtr_inv)
print("q_nt =")
sp.pprint(qnt_inv)
print("q_eta =")
sp.pprint(qeta_inv)

expect_zero("q_tr inverse", sp.simplify(qtr_inv - qtr))
expect_zero("q_nt inverse", sp.simplify(qnt_inv - qnt))
expect_zero("q_eta inverse", sp.simplify(qeta_inv - qeta))

subbanner("III. Exact nondegeneracy on the physical reference branch")
detD = sp.factor(defect_matrix.det())
print("det(D_q) =")
sp.pprint(detD)
expect_zero(
    "det(D_q) - chi0*deltaU*epsilon_eta / ((1+chi0)(1+deltaU)(1+chi0+deltaU)(1-epsilon_eta))",
    sp.simplify(detD - chi0s * deltaUs * epsEtas / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs) * (1 - epsEtas))),
)

subbanner("IV. Zero-defect theorem in actual-branch quotient language")
print("Theta_1 = Xi_1 = R_1 = 0  iff  q_tr = q_nt = q_eta = 0.")
print("So an actual moving-throat branch can be tested in either of two exactly equivalent first-order languages:")
print("  1. finite quotient residuals (q_tr,q_nt,q_eta), or")
print("  2. physical defect triplet (Theta_1, Xi_1, R_1).")
