#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage187_192_common import *

"""
Stage 189 — exact reference-orbit transport laws in ratio form.

What this script does
---------------------
1. Rewrites the finite orbit action as exact multiplicative transport laws relative to
   a chosen reference point on the orbit.
2. Expresses the dependent-coordinate ratios purely in terms of the five free-coordinate
   ratios.
3. Verifies those laws directly against the finite orbit action.
"""

banner("STAGE 189 — EXACT REFERENCE-ORBIT TRANSPORT LAWS")

syms = base_symbols()
rat = free_ratio_symbols()
laws = orbit_ratio_laws(syms, rat)

subbanner("I. Exact dependent-ratio laws from free-coordinate ratios")
print("R_Keta^(orbit) =")
sp.pprint(laws["RKeta"])
print("R_TU^(orbit) =")
sp.pprint(laws["RTU"])
print("R_muW^(orbit) =")
sp.pprint(laws["Rmu"])

subbanner("II. Direct verification from the finite orbit action")
pars = {
    "Lam": sp.log(rat["Rlam"]),
    "C": sp.log(rat["Rc"]),
    "Gam": sp.log(rat["Rgam"]),
    "U": sp.log(rat["RU"]),
    "W": sp.log(rat["RW"]),
}
syms1 = orbit_transform(syms, pars)
expect_zero("R_Keta law", sp.simplify((syms1["Keta"] / syms["Keta"]) / laws["RKeta"] - 1))
expect_zero("R_TU law", sp.simplify((syms1["TU"] / syms["TU"]) / laws["RTU"] - 1))
expect_zero("R_muW law", sp.simplify((syms1["muW"] / syms["muW"]) / laws["Rmu"] - 1))

subbanner("III. Orbit-lock criterion relative to a reference orbit point")
print("Given a reference orbit point x_ref and an actual candidate with free-coordinate ratios")
print("(R_lambda, R_c, R_gamma, R_U, R_W), the exact single-orbit prediction is:")
print("  K_eta(actual) / K_eta(ref) = R_Keta^(orbit)")
print("  T_U(actual)   / T_U(ref)   = R_TU^(orbit)")
print("  mu_W(actual)  / mu_W(ref)  = R_muW^(orbit)")
print("Any violation of these three transport laws is exact off-orbit motion.")

banner("STAGE 189 LEDGER")
print("1. The finite orbit law is now an exact coevolution law relative to a reference point.")
print("2. The five free-coordinate ratios determine the three dependent-coordinate ratios uniquely.")
print("3. This is the direct finite version of the earlier tangent-space orbit law.")
