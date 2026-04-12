#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage280_283_common import banner, subbanner, expect_zero

banner("STAGE 283 — END-TO-END SELECTED-BRANCH OUTGOING BRIDGE")

sigma_star, gamma_star = sp.symbols('sigma_star gamma_star', positive=True, real=True)
eps_blk = sp.symbols('eps_blk', real=True)
Xi_slip, dPi_tan = sp.symbols('Xi_slip delta_Pi_tan', real=True)

rho_alpha = sp.Rational(4, 3)
zeta_req = sp.simplify((rho_alpha - 1)/(1 - eps_blk*(2 - rho_alpha)))
Pi_over_C = sp.simplify((1 + (1 - 2*eps_blk)*zeta_req)/(1 - eps_blk*zeta_req))
chi_Q = sp.simplify((1 - 9*sigma_star*gamma_star)/(1 - sigma_star))

subbanner("I. Conservative selected branch")
print("rho_alpha =")
sp.pprint(rho_alpha)
print("zeta_req(eps_blk) =")
sp.pprint(zeta_req)
print("Pi_tr / C_mix =")
sp.pprint(Pi_over_C)
expect_zero("conservative ratio remains 4/3", Pi_over_C - sp.Rational(4, 3))

subbanner("II. Outgoing outlet branch")
print("chi_Q(sigma_W,gamma_W) =")
sp.pprint(chi_Q)
chi_can = sp.simplify(chi_Q.subs(gamma_star, sp.Rational(1, 9)))
print("chi_Q on the canonical odd point gamma_W = 1/9 =")
sp.pprint(chi_can)
expect_zero("canonical odd point preserves chi_Q", chi_can - 1)

subbanner("III. Tangent transport from selected loading to outgoing defect")
dgamma_from_load = sp.simplify(Xi_slip*dPi_tan/sp.Integer(9))
DeltaQ = sp.simplify(-sigma_star*Xi_slip*dPi_tan/(1 - sigma_star))
print("delta gamma_W =")
sp.pprint(dgamma_from_load)
print("Delta_Q =")
sp.pprint(DeltaQ)

subbanner("IV. Final one-language theorem ledger")
print("Selected conservative branch:")
print("  Pi_tr / C_mix = 4/3")
print("  zeta_req = 1 / (3 - 2 eps_blk)")
print()
print("Compensated canonical-even outlet:")
print("  chi_Q = (1 - 9 sigma_W gamma_W) / (1 - sigma_W)")
print("  chi_Q = 1  iff  gamma_W = 1/9  on the nontrivial branch")
print()
print("Tangent loading-to-outlet bridge:")
print("  delta gamma_W = Xi_slip * delta Pi_tan / 9")
print("  Delta_Q = -(sigma_*/(1-sigma_*)) Xi_slip delta Pi_tan")
print()
print("So the final outgoing theorem gap is written in one language end to end:")
print("  conservative selected loading/product side  ->  tangential load slippage  ->  outlet odd slippage  ->  chi_Q")
