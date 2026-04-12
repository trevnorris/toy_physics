#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage280_283_common import banner, subbanner, expect_zero

banner("STAGE 280 — COMPENSATED HYBRID OUTLET AND SELECTED LOADING/PRODUCT BRIDGE")

z = sp.symbols('z', real=True)
sigma_W, gamma_W = sp.symbols('sigma_W gamma_W', real=True)
eps_blk = sp.symbols('eps_blk', real=True)

Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + sp.I*z**5/sp.Integer(9)
Lambda_hyb = Lambda_out + 4*sigma_W - sigma_W/(1 - z**2/sp.Integer(3) - sp.I*gamma_W*z**5)
Yhat_hyb = sp.expand(sp.series((Lambda_hyb.subs(z, 0))/Lambda_hyb, z, 0, 6).removeO())
chi_Q = sp.simplify((1 - 9*sigma_W*gamma_W)/(1 - sigma_W))

subbanner("I. Exact compensated-hybrid normalized response")
print("Yhat_hyb(z) =")
sp.pprint(Yhat_hyb)
print("chi_Q =")
sp.pprint(chi_Q)

Yhat_module = sp.expand(
    sp.series(
        sp.Rational(3, 4)
        + sp.Rational(1, 4)/(1 - sp.Rational(4, 9)*z**2 - sp.I*sp.Rational(4, 27)*chi_Q*z**5),
        z, 0, 6,
    ).removeO()
)
print("Equivalent minimal contact-plus-pole module =")
sp.pprint(Yhat_module)
expect_zero("compensated hybrid minus minimal module", Yhat_hyb - Yhat_module)

subbanner("II. Exact selected loading/product ratio on the canonical-even branch")
rho_alpha = sp.Rational(4, 3)
zeta_req = sp.simplify((rho_alpha - 1)/(1 - eps_blk*(2 - rho_alpha)))
Q = sp.simplify((1 + (1 - 2*eps_blk)*zeta_req)/(1 - eps_blk*zeta_req))
print("rho_alpha =")
sp.pprint(rho_alpha)
print("zeta_req(eps_blk) =")
sp.pprint(zeta_req)
print("Pi_tr / C_mix = Q(zeta_req, eps_blk) =")
sp.pprint(Q)
expect_zero("exact product ratio is fixed", Q - sp.Rational(4, 3))

subbanner("III. Conservative/odd factorization theorem")
print("On the compensated canonical-even outlet:")
print("  c0 = 3/4")
print("  c1 = 1/4")
print("  Pi_tr / C_mix = 4/3  (independent of eps_blk once zeta_req is imposed)")
print("  chi_Q = (1 - 9 sigma_W gamma_W) / (1 - sigma_W)")
print()
print("So sigma_W and gamma_W do not move the conservative selected loading ratio; they move only the outgoing normalization factor chi_Q.")
