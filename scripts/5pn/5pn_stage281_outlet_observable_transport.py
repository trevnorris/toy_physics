#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage280_283_common import banner, subbanner, expect_zero

banner("STAGE 281 — OUTLET-OBSERVABLE TRANSPORT ON THE COMPENSATED CANONICAL-EVEN BRANCH")

sigma_star, gamma_star = sp.symbols('sigma_star gamma_star', positive=True, real=True)
dsigma, dgamma = sp.symbols('delta_sigma delta_gamma', real=True)
dC, dkappaW = sp.symbols('delta_C delta_kappa_W', real=True)

chi = sp.simplify((1 - 9*sigma_star*gamma_star)/(1 - sigma_star))
subbanner("I. General compensated-hybrid outgoing factor and its canonical differential")
print("chi_Q(sigma,gamma) =")
sp.pprint(chi)

dchi = sp.simplify(sp.diff(chi, sigma_star)*dsigma + sp.diff(chi, gamma_star)*dgamma)
print("General first differential d chi_Q =")
sp.pprint(dchi)

dchi_can = sp.simplify(dchi.subs(gamma_star, sp.Rational(1, 9)))
print("At the canonical outgoing point gamma_* = 1/9:")
sp.pprint(dchi_can)
expect_zero(
    "canonical dchi minus pure dgamma term",
    dchi_can + 9*sigma_star*dgamma/(1 - sigma_star),
)

subbanner("II. Exact canonical-even gate from the outlet projection")
dE2 = sp.simplify((dC - 9*sigma_star*dkappaW)/(27*(1 - sigma_star)))
dE4 = sp.simplify((5*dC - 72*sigma_star*dkappaW)/(243*(1 - sigma_star)))
DeltaQ = sp.simplify((dC - 27*sigma_star*dgamma)/(3*(1 - sigma_star)))
print("deltaE2 =")
sp.pprint(dE2)
print("deltaE4 =")
sp.pprint(dE4)
print("Delta_Q =")
sp.pprint(DeltaQ)

sol = sp.solve([sp.Eq(dE2, 0), sp.Eq(dE4, 0)], [dC, dkappaW], dict=True)
print("Solve deltaE2 = deltaE4 = 0 gives:")
print(sol)
sol0 = sol[0]
DeltaQ_can = sp.simplify(DeltaQ.subs(sol0))
print("Delta_Q on the canonical-even tangent =")
sp.pprint(DeltaQ_can)
expect_zero(
    "Delta_Q canonical tangent law",
    DeltaQ_can + 9*sigma_star*dgamma/(1 - sigma_star),
)

subbanner("III. Canonical branch meaning")
print("Once the canonical-even gate is imposed, sigma_W and kappa_W are frozen to first order.")
print("The only surviving isotropic outlet defect is the odd slippage d gamma_W, amplified by the static factor sigma_*/(1-sigma_*).")
