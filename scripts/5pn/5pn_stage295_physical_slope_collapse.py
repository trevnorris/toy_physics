#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage289_292_common import *

"""
Stage 295 — collapse of the weak grouped outlet problem to the physical slopes
u_2^(1) and P_1/P_0.
"""

banner("STAGE 295 — PHYSICAL SLOPE COLLAPSE")

D0, N0, P0, sigma = sp.symbols("D_0 N_0 P_0 sigma_star", positive=True, real=True)
eps = sp.symbols("eps", real=True)
K1m, G1m = sp.symbols("mathfrakK1 mathfrakG1", real=True)
nu21, P1 = sp.symbols("u2_1 P_1", real=True)

# Standard weak-axisymmetric Y20 grouped signature.
l20 = sp.Integer(1)
l21 = sp.Rational(1, 2)
l22 = -sp.Integer(1)

subbanner("I. Weak-axisymmetric grouped signature")
# Standard weak-axisymmetric Y20 pattern.

# Physical grouped slopes.
du2_20 = eps * l20 * nu21
du2_21 = eps * l21 * nu21
du2_22 = eps * l22 * nu21
dP0_20 = eps * l20 * P1
dP0_21 = eps * l21 * P1
dP0_22 = eps * l22 * P1

print("delta u2^(A) = eps lambda_A u2^(1)")
print("delta P0^(A) = eps lambda_A P1")

subbanner("II. Microscopic obstruction amplitudes are the physical slopes")
# Stage-154/155 dictionary.
K1_from_phys = sp.simplify(-D0 * nu21)
G1_from_phys = sp.simplify(D0 * P1)
print("mathfrak K_1 =")
sp.pprint(K1_from_phys)
print("mathfrak G_1 =")
sp.pprint(G1_from_phys)

# Direct outlet deformations.
dkappa20 = sp.simplify(-3 * (1 - sigma) * du2_20 / sigma)
dkappa21 = sp.simplify(-3 * (1 - sigma) * du2_21 / sigma)
dkappa22 = sp.simplify(-3 * (1 - sigma) * du2_22 / sigma)

dgamma20 = sp.simplify(-(1 - sigma) * dP0_20 / (9 * sigma * P0))
dgamma21 = sp.simplify(-(1 - sigma) * dP0_21 / (9 * sigma * P0))
dgamma22 = sp.simplify(-(1 - sigma) * dP0_22 / (9 * sigma * P0))

print("delta kappa_W^(20) =")
sp.pprint(dkappa20)
print("delta kappa_W^(21) =")
sp.pprint(dkappa21)
print("delta kappa_W^(22) =")
sp.pprint(dkappa22)
print("delta gamma_W^(20) =")
sp.pprint(dgamma20)
print("delta gamma_W^(21) =")
sp.pprint(dgamma21)
print("delta gamma_W^(22) =")
sp.pprint(dgamma22)

# Scalar amplitudes.
kappa1_amp = sp.simplify(-3 * (1 - sigma) * nu21 / sigma)
gamma1_amp = sp.simplify(-(1 - sigma) * P1 / (9 * sigma * P0))
print("kappa_1 =")
sp.pprint(kappa1_amp)
print("gamma_1 =")
sp.pprint(gamma1_amp)

subbanner("III. Grouped trace/anomaly form")
# Hidden outlet grouped anomalies inherit the same Y20 signature.
gk = grouped_trace_anomaly(dkappa20, dkappa21, dkappa22)
gg = grouped_trace_anomaly(dgamma20, dgamma21, dgamma22)
print("a_kappa =")
sp.pprint(gk["a"])
print("b_kappa =")
sp.pprint(gk["b"])
print("a_gamma =")
sp.pprint(gg["a"])
print("b_gamma =")
sp.pprint(gg["b"])
expect_zero("b_kappa - 3 a_kappa", sp.simplify(gk["b"] - 3 * gk["a"]))
expect_zero("b_gamma - 3 a_gamma", sp.simplify(gg["b"] - 3 * gg["a"]))

subbanner("IV. Interpretation")
print("On the weak axisymmetric grouped branch, the whole linear outlet problem")
print("collapses to two physical slopes only:")
print("  u2^(1) and P1/P0.")
print("Equivalently, the microscopic obstruction amplitudes are")
print("  mathfrak K_1 = -D_0 u2^(1),")
print("  mathfrak G_1 =  D_0 P_1.")
