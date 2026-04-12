#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage284_288_common import *

"""
Stage 285 — bare mixed-port slippage and exact D/N similarity decomposition.

What this script does
---------------------
1. Re-derives the compensated-branch core identity relating delta gamma_W and
   delta kappa_W to the bare mixed-port slippage scalar.
2. Introduces the tangential susceptibility delta B_W = Upsilon_Pi delta Pi_tan.
3. Rewrites the susceptibility in D/N similarity form and proves that the final
   outlet defect collapses to Delta_Q = -(sigma_*/(1-sigma_*)) Xi_slip delta Pi_tan.
"""

banner("STAGE 285 — BARE MIXED-PORT SLIPPAGE AND D/N SIMILARITY")

rcs = sp.symbols("r_c_star", positive=True, real=True)
drc, dk0, dg0, dkW = sp.symbols("drc dkappa0 dgamma0 dkappaW", real=True)
sigma_star, Upsilon_Pi, dPi_tan, Xi_gamma, Xi_L = sp.symbols(
    "sigma_star Upsilon_Pi deltaPi_tan Xi_gamma Xi_L", real=True
)

subbanner("I. Exact compensated-branch core identity")
k0_star = (1 + rcs) / 3
g0_star = (1 + rcs) / 9

dkW_expr = sp.simplify(dk0 / (1 + rcs) - k0_star * drc / (1 + rcs) ** 2)
dgW_expr = sp.simplify(dg0 / (1 + rcs) - g0_star * drc / (1 + rcs) ** 2)
identity = sp.simplify(dgW_expr - dkW_expr / 3 - (dg0 - dk0 / 3) / (1 + rcs))
print("delta kappa_W =")
sp.pprint(dkW_expr)
print("delta gamma_W =")
sp.pprint(dgW_expr)
expect_zero("compensated branch identity", identity)

dB = sp.symbols("deltaB_W", real=True)
print("If delta kappa_W = 0, then")
print("delta gamma_W =")
sp.pprint(dB / (1 + rcs))

subbanner("II. Tangential susceptibility and outlet defect")
dgW_from_tan = sp.simplify(Upsilon_Pi * dPi_tan / (1 + rcs))
DeltaQ = sp.simplify(-9 * sigma_star / (1 - sigma_star) * dgW_from_tan)
print("delta gamma_W from tangential mouth deformation =")
sp.pprint(dgW_from_tan)
print("Delta_Q =")
sp.pprint(DeltaQ)

subbanner("III. Exact D/N similarity decomposition")
dln_gamma0, dln_Lratio = sp.symbols("dln_gamma0 dln_Lratio", real=True)
# On the compensated branch: gamma0_* = (1+r_c*)/9 and kappa0_* = 4 (L_W/a)^2 / pi^2 = (1+r_c*)/3.
# Therefore deltaB_W = gamma0_* dln gamma0 - (1/3) kappa0_* * 2 dln(L_W/a).
dB_from_similarity = sp.simplify((1 + rcs) / 9 * (dln_gamma0 - 2 * dln_Lratio))
print("delta B_W from D/N similarity strain =")
sp.pprint(dB_from_similarity)

Xi_slip = sp.simplify(Xi_gamma - 2 * Xi_L)
Upsilon_from_similarity = sp.simplify((1 + rcs) / 9 * Xi_slip)
print("Upsilon_Pi =")
sp.pprint(Upsilon_from_similarity)

DeltaQ_similarity = sp.simplify(
    -9 * sigma_star / ((1 - sigma_star) * (1 + rcs)) * Upsilon_from_similarity * dPi_tan
)
print("Delta_Q after similarity decomposition =")
sp.pprint(DeltaQ_similarity)
expect_zero(
    "Delta_Q + (sigma_*/(1-sigma_*)) Xi_slip deltaPi_tan",
    sp.simplify(DeltaQ_similarity + sigma_star * Xi_slip * dPi_tan / (1 - sigma_star)),
)

subbanner("IV. Final compact defect law")
print("Define Xi_slip := Xi_gamma - 2 Xi_L. Then")
print("Delta_Q =")
sp.pprint(DeltaQ_similarity)
