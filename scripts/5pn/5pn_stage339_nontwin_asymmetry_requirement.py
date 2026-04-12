#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage336_339_common import *

"""
Stage 339 — exact non-twin asymmetry requirement once rho_alpha exceeds 2.

What this script does
---------------------
1. Rewrites the exact required support ratio zeta_req directly in rho_alpha.
2. Proves zeta_req > 1 iff rho_alpha > 2.
3. Derives the exact excess-above-twin formula
      zeta_req - 1 = (1-eps)(rho_alpha - 2) / (1 + eps(rho_alpha - 2)).
4. Composes this with the selected-branch loading ratio rho_alpha(xi) to obtain
   the exact non-twin asymmetry demand on the physical branch.
"""

banner("STAGE 339 — EXACT NON-TWIN ASYMMETRY REQUIREMENT")

rho_alpha, eps = sp.symbols("rho_alpha eps", positive=True, real=True)
zeta_req = zeta_req_from_rho(rho_alpha, eps)

subbanner("I. Exact support-ratio demand")
print("zeta_req(rho_alpha, eps) =")
sp.pprint(zeta_req)

subbanner("II. Exact twin/non-twin threshold")
zeta_excess = sp.simplify(zeta_req - 1)
print("zeta_req - 1 =")
sp.pprint(zeta_excess)

expect_zero(
    "zeta_req - 1 - (1-eps)(rho_alpha-2)/(1+eps(rho_alpha-2))",
    zeta_excess - sp.simplify((1 - eps) * (rho_alpha - 2) / (1 + eps * (rho_alpha - 2))),
)

dzeta_drho = sp.simplify(sp.diff(zeta_req, rho_alpha))
print("d zeta_req / d rho_alpha =")
sp.pprint(dzeta_drho)

expect_zero(
    "dzeta/drho - (1-eps)/(1+eps(rho-2))^2",
    dzeta_drho - sp.simplify((1 - eps) / (1 + eps * (rho_alpha - 2))**2),
)

print("Therefore:")
print("  rho_alpha <= 1   -> mixed-only enough")
print("  1 < rho_alpha <= 2 -> lowest twin enough (0 < zeta_req <= 1)")
print("  rho_alpha > 2    -> non-twin asymmetry required (zeta_req > 1)")

subbanner("III. Exact branch-level non-twin demand")
xi, delta, R, Mmix = sp.symbols("xi delta R M_mix", positive=True, real=True)
rho_branch = sp.simplify(G_tr(xi, delta, R) / Mmix)
zeta_branch = sp.simplify(zeta_req.subs(rho_alpha, rho_branch))

print("rho_alpha^(branch) =")
sp.pprint(rho_branch)
print("zeta_req^(branch) =")
sp.pprint(zeta_branch)

print("So the physical branch needs non-twin asymmetry exactly when")
print("  G_tr(xi,delta;R) > 2 M_mix,")
print("equivalently when")
print("  xi_phys > xi_(2x).")
