
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 143 — BARE MIXED-PORT SLIPPAGE")

rc, sigma = sp.symbols('r_c sigma', positive=True, real=True)
dgamma0, dkappa0, drc = sp.symbols('deltaGamma0 deltaKappa0 deltaRc', real=True)
dgammaW, dkappaW = sp.symbols('deltaGammaW deltaKappaW', real=True)
UpsilonPi, dPi_tan = sp.symbols('UpsilonPi deltaPi_tan', real=True)
r = sp.symbols('r', positive=True, real=True)

subbanner("1. Core identities on the compensated branch")
kappa0_star = (1 + rc) / 3
gamma0_star = (1 + rc) / 9
dkappaW_expr = dkappa0 / (1 + rc) - kappa0_star * drc / (1 + rc) ** 2
dgammaW_expr = dgamma0 / (1 + rc) - gamma0_star * drc / (1 + rc) ** 2
identity = sp.simplify(dgammaW_expr - dkappaW_expr / 3 - (dgamma0 - dkappa0 / 3) / (1 + rc))
expect_zero("deltaGammaW - deltaKappaW/3 identity", identity)
print("deltaKappaW =", dkappaW_expr)
print("deltaGammaW =", dgammaW_expr)

subbanner("2. Collapse under the canonical-even gate")
dB = sp.symbols('deltaB_W', real=True)
collapsed = sp.simplify(dgammaW_expr.subs(dgamma0, dB + dkappa0 / 3).subs(dkappaW_expr, 0))
# Simpler directly:
collapsed = sp.simplify((dB) / (1 + rc))
print("deltaGammaW =", collapsed)

subbanner("3. Tangential DtN susceptibility")
dB_expr = UpsilonPi * dPi_tan
dgammaW_from_U = sp.simplify(dB_expr / (1 + rc))
print("deltaGammaW from UpsilonPi =", dgammaW_from_U)

DeltaQ = -9 * sigma * dgammaW_from_U / (1 - sigma)
DeltaQ_r = sp.simplify(DeltaQ.subs(rc, r**2))
print("DeltaQ =", DeltaQ_r)
print("N_Q - 1 =", sp.simplify(-DeltaQ_r))
print("numeric 9/(1+r_F1^2) =", sp.N((9 / (1 + rF1**2)), 16))
