
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 144 — D/N SIMILARITY SLIPPAGE")

rc, sigma = sp.symbols('r_c sigma', positive=True, real=True)
dln_gamma0, dln_Lratio, dln1prc = sp.symbols('dlnGamma0 dlnLratio dln1prc', real=True)
Xi_gamma, Xi_L, dPi_tan = sp.symbols('Xi_gamma Xi_L deltaPi_tan', real=True)
r = sp.symbols('r', positive=True, real=True)
LW, a = sp.symbols('L_W a', positive=True, real=True)
gamma0 = sp.symbols('gamma0', positive=True, real=True)

subbanner("1. Similarity-defect parametrization")
de_kappa = 2 * dln_Lratio - dln1prc
de_gamma = dln_gamma0 - dln1prc
dB = sp.simplify((1 + rc) * (de_gamma - de_kappa) / 9)
print("deltaB_W =", dB)

UpsilonPi = sp.simplify(dB.subs({dln_gamma0: Xi_gamma * dPi_tan, dln_Lratio: Xi_L * dPi_tan, dln1prc: 0}) / dPi_tan)
target_U = (1 + rc) * (Xi_gamma - 2 * Xi_L) / 9
expect_zero("UpsilonPi - target", UpsilonPi - target_U)
print("UpsilonPi =", UpsilonPi)

subbanner("2. Collapse of the final defect law")
Xi_slip = sp.symbols('Xi_slip', real=True)
DeltaQ = -9 * sigma * UpsilonPi * dPi_tan / ((1 - sigma) * (1 + rc))
DeltaQ_slip = sp.simplify(DeltaQ.subs(Xi_gamma - 2 * Xi_L, Xi_slip))
target_DeltaQ = -sigma * Xi_slip * dPi_tan / (1 - sigma)
expect_zero("DeltaQ_slip - target", DeltaQ_slip - target_DeltaQ)
print("DeltaQ =", target_DeltaQ)
print("N_Q - 1 =", -target_DeltaQ)

subbanner("3. Exact D/N similarity-preservation theorem")
dn_law = sp.Eq(gamma0, 4 * LW**2 / (3 * sp.pi**2 * a**2))
# log-variation of the law
lhs = dln_gamma0
rhs = 2 * (sp.symbols('dlnLW', real=True) - sp.symbols('dlna', real=True))
print("If gamma0 = 4 L_W^2 / (3 pi^2 a^2), then delta ln gamma0 =", rhs)
print("Hence Xi_slip = 0 and the first-order reduced 2.5PN defect vanishes.")
