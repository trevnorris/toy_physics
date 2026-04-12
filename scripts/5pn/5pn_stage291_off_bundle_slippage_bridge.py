#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage289_292_common import *

"""
Stage 291 — exact off-bundle slippage bridge.
"""

banner("STAGE 291 — OFF-BUNDLE SLIPPAGE BRIDGE")

refs = family1_reference_values()
g_star = refs["g_star"]
r_star = refs["r_F1"]
Sigma0_can = refs["Sigma0_can"]
S_can = refs["S_can"]

# Stage-147 microscopic channels.
dlnZratio, dlncsw, dlncs, dlna = sp.symbols("dln_Zq_over_rho dln_csw dln_cs dln_a", real=True)

# Lower-branch transport laws plus slippages.
epsL, epsv, epsT = sp.symbols("eps_L eps_v eps_T", real=True)
dlnLW = dlna + epsL
dlnvw0 = sp.Rational(1, 2) * dlnZratio + sp.Rational(3, 2) * dlncsw + dlncs - sp.Rational(5, 2) * dlna + epsv
dlnTm = sp.Rational(1, 2) * dlnZratio + sp.Rational(3, 2) * dlncsw - dlncs - sp.Rational(3, 2) * dlna + epsT

Astar = sp.simplify(g_star + 1 / (4 * sp.sqrt(1 + r_star**2)))
Bstar = sp.simplify(1 / (2 * sp.sqrt(1 + r_star**2)))
Cstar = sp.simplify(2 * g_star + sp.Rational(3, 1) / (4 * sp.sqrt(1 + r_star**2)))

delta_perp = sp.simplify(
    Astar * dlnZratio + 3 * Astar * dlncsw + Bstar * dlncs
    - g_star * dlnTm - (g_star + Bstar) * dlnvw0 - 2 * Astar * dlna - Cstar * dlnLW
)
eps_perp = sp.simplify(g_star * epsT + (g_star + Bstar) * epsv + Cstar * epsL)

eps_perp_sym = sp.symbols("eps_perp", real=True)

subbanner("I. Exact collapse of the off-family scalar")
print("delta_perp =")
sp.pprint(delta_perp)
print("eps_perp =")
sp.pprint(eps_perp)
expect_zero("delta_perp + eps_perp", sp.simplify(delta_perp + eps_perp))

subbanner("II. Numerical Family-1 weights")
num_g = sp.N(g_star)
num_gB = sp.N(g_star + Bstar)
num_C = sp.N(Cstar)
num_epspref = sp.N(Sigma0_can * S_can / sp.sqrt(1 + r_star**2))
print(f"coefficient of eps_T   = {-num_g}")
print(f"coefficient of eps_v   = {-num_gB}")
print(f"coefficient of eps_L   = {-num_C}")
print(f"mouth-bias coefficient = {-num_epspref}")

subbanner("III. Mouth-bias transport")
dSigma0, dScal = sp.symbols("dSigma0 dScal", real=True)
dPi_tan = sp.simplify((1 - sp.Rational(1, 4) * S_can) * dSigma0 - Sigma0_can * dScal / 4)
dPi_full = sp.simplify(dPi_tan - Sigma0_can * S_can * eps_perp / sp.sqrt(1 + r_star**2))
print("deltaPi_tan =")
sp.pprint(dPi_tan)
print("deltaPi =")
sp.pprint(dPi_full)
print("deltaPi (numerical coefficients) =")
sp.pprint(sp.expand(sp.N(dPi_full)))

subbanner("IV. Direct outlet defect ledger")
sigma_star, dkappaW, dgammaW = sp.symbols("sigma_star delta_kappa_W delta_gamma_W", real=True)
deltaC = sp.simplify(-16 * sigma_star * eps_perp_sym / sp.sqrt(1 + r_star**2))
dE2 = sp.simplify(sigma_star * (-16 * eps_perp_sym / sp.sqrt(1 + r_star**2) - 9 * dkappaW) / (27 * (1 - sigma_star)))
dE4 = sp.simplify(sigma_star * (-80 * eps_perp_sym / sp.sqrt(1 + r_star**2) - 72 * dkappaW) / (243 * (1 - sigma_star)))
DeltaQ = sp.simplify(sigma_star * (-16 * eps_perp_sym / sp.sqrt(1 + r_star**2) - 27 * dgammaW) / (3 * (1 - sigma_star)))
print("delta C =")
sp.pprint(deltaC)
print("delta E2 =")
sp.pprint(dE2)
print("delta E4 =")
sp.pprint(dE4)
print("Delta_Q =")
sp.pprint(DeltaQ)

subbanner("V. Even-preservation corollary")
DeltaQ_even_preserved = sp.simplify(DeltaQ.subs({eps_perp_sym: 0, dkappaW: 0}))
print("Delta_Q |_{eps_perp = 0, delta kappa_W = 0} =")
sp.pprint(DeltaQ_even_preserved)
expect_zero(
    "DeltaQ even-preserved form",
    sp.simplify(DeltaQ_even_preserved + 9 * sigma_star * dgammaW / (1 - sigma_star)),
)

subbanner("VI. Interpretation")
print("All first-order off-bundle normal motion is carried by one scalar eps_perp.")
print("The mouth-bias and conservative-even defects depend on the same scalar, while")
print("the remaining odd normalization defect still needs delta gamma_W unless")
print("eps_perp and delta kappa_W are both killed.")
