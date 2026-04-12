#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage333_335_common import *

"""
Stage 333 — exact branch-product support regimes.

What this script does
---------------------
1. Rebuilds the exact support-demand map in the branch-product language
   (Pi_tr, C_mix, eps).
2. Proves the three exact regimes:
      Pi_tr <= C_mix              -> mixed-only enough,
      C_mix < Pi_tr <= 2 C_mix    -> symmetric lowest twin enough,
      Pi_tr > 2 C_mix             -> non-twin asymmetry required.
3. Verifies the exact derivative dzeta_req / dPi_tr > 0.
"""

banner("STAGE 333 — EXACT BRANCH-PRODUCT SUPPORT REGIMES")

Pi_tr, Lambda, eps = sp.symbols("Pi_tr Lambda eps", positive=True, real=True)
Cmix = c_mix(Lambda, eps)
Sreq = s_req_from_products(Pi_tr, Cmix)
zeta_req = zeta_req_from_products(Pi_tr, Cmix, eps)

subbanner("I. Exact demand map")
print("C_mix =")
sp.pprint(Cmix)
print("S_req = Pi_tr / C_mix =")
sp.pprint(Sreq)
print("zeta_req(Pi_tr,C_mix,eps) =")
sp.pprint(zeta_req)

expect_zero(
    "zeta_req - (S_req - 1)/(1 + eps*(S_req - 2))",
    zeta_req - sp.simplify((Sreq - 1) / (1 + eps*(Sreq - 2))),
)

subbanner("II. Exact regime thresholds")
zeta_at_mixed = sp.simplify(zeta_req.subs(Pi_tr, Cmix))
zeta_at_twin = sp.simplify(zeta_req.subs(Pi_tr, 2*Cmix))
print("zeta_req at Pi_tr = C_mix ->")
sp.pprint(zeta_at_mixed)
print("zeta_req at Pi_tr = 2 C_mix ->")
sp.pprint(zeta_at_twin)
expect_zero("zeta(C_mix)", zeta_at_mixed)
expect_zero("zeta(2 C_mix) - 1", zeta_at_twin - 1)

print("Exact regime split:")
print("  Pi_tr <= C_mix              -> mixed-only already enough")
print("  C_mix < Pi_tr <= 2 C_mix    -> symmetric lowest twin enough")
print("  Pi_tr > 2 C_mix             -> non-twin asymmetry required")

subbanner("III. Strict monotonicity of the support demand")
dzeta_dPi = sp.factor(sp.diff(zeta_req, Pi_tr))
print("d zeta_req / d Pi_tr =")
sp.pprint(dzeta_dPi)

expect_zero(
    "dzeta_dPi denominator identity",
    sp.simplify(
        dzeta_dPi
        - Cmix*(1 - eps) / (Cmix - eps*(2*Cmix - Pi_tr))**2
    ),
)

print("Since C_mix > 0 and 0 < eps < 1 on the blocked branch, dzeta_req/dPi_tr > 0.")
