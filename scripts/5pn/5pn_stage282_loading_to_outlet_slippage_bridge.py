#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage280_283_common import banner, subbanner, expect_zero

banner("STAGE 282 — LOADING/PRODUCT TO OUTLET-SLIPPAGE BRIDGE")

sigma_star = sp.symbols('sigma_star', positive=True, real=True)
Xi_slip, dPi_tan = sp.symbols('Xi_slip delta_Pi_tan', real=True)
P1, P0 = sp.symbols('P_1 P_0', real=True, nonzero=True)

dgamma = sp.symbols('delta_gamma_W', real=True)

subbanner("I. Exact tangent transport from Family-1 loading to outlet odd slippage")
DeltaQ_from_load = sp.simplify(-sigma_star*Xi_slip*dPi_tan/(1 - sigma_star))
DeltaQ_from_gamma = sp.simplify(-9*sigma_star*dgamma/(1 - sigma_star))
print("Delta_Q from tangential loading =")
sp.pprint(DeltaQ_from_load)
print("Delta_Q from outlet odd slippage =")
sp.pprint(DeltaQ_from_gamma)

dgamma_from_load = sp.simplify(sp.solve(sp.Eq(DeltaQ_from_load, DeltaQ_from_gamma), dgamma)[0])
print("delta gamma_W implied by the two exact Delta_Q formulas =")
sp.pprint(dgamma_from_load)
expect_zero(
    "delta gamma transport law",
    dgamma_from_load - Xi_slip*dPi_tan/sp.Integer(9),
)

subbanner("II. Weak-axisymmetric grouped form")
gamma1 = sp.simplify(-(1 - sigma_star)/(9*sigma_star) * P1/P0)
DeltaQ_from_P = sp.simplify(-9*sigma_star*gamma1/(1 - sigma_star))
print("gamma_1 =")
sp.pprint(gamma1)
print("Delta_Q from gamma_1 =")
sp.pprint(DeltaQ_from_P)
expect_zero("Delta_Q equals P1/P0", DeltaQ_from_P - P1/P0)

P_ratio_from_load = sp.simplify(sp.solve(sp.Eq(DeltaQ_from_load, P1/P0), P1/P0)[0])
print("P1/P0 implied by the loading-side transport =")
sp.pprint(P_ratio_from_load)
expect_zero(
    "P1/P0 loading bridge",
    P_ratio_from_load + sigma_star*Xi_slip*dPi_tan/(1 - sigma_star),
)

subbanner("III. End-to-end meaning")
print("The final outgoing theorem gap is one scalar, written equivalently as:")
print("  Delta_Q")
print("  P1/P0")
print("  Xi_slip * delta Pi_tan")
print()
print("with exact bridge laws:")
print("  delta gamma_W = Xi_slip * delta Pi_tan / 9")
print("  Delta_Q       = -(sigma_*/(1-sigma_*)) Xi_slip delta Pi_tan")
print("  Delta_Q       = P1/P0")
