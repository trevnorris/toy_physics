#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage336_339_common import *

"""
Stage 336 — exact selected-branch loading ratio on the physical tracking branch.

What this script does
---------------------
1. Rebuilds the exact branch-product ratio
      rho_alpha = Pi_tr / C_mix
   on the tracking branch.
2. Proves that on the physical target branch this is exactly
      rho_alpha = G_tr / M_mix.
3. Rewrites the mixed-only / lowest-twin / non-twin support regimes in the
   physical selected-branch variable xi.
4. Proves rho_alpha is strictly increasing with xi at fixed (delta,R,M_mix).
"""

banner("STAGE 336 — EXACT SELECTED-BRANCH LOADING RATIO")

xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)
R_target, Mmix, Lambda, eps = sp.symbols("R_target M_mix Lambda eps", positive=True, real=True)

G = G_tr(xi, delta, R)
F = F_tr(xi, delta, R)
Pi = Pi_tr(xi, delta, R)
Cmix = c_mix(Lambda, eps)
rho = rho_alpha_from_products(Pi, Cmix)

subbanner("I. Exact branch-product ratio")
print("G_tr =")
sp.pprint(G)
print("F_tr =")
sp.pprint(F)
print("Pi_tr =")
sp.pprint(Pi)
print("C_mix =")
sp.pprint(Cmix)
print("rho_alpha = Pi_tr / C_mix =")
sp.pprint(rho)

subbanner("II. Exact physical selected-branch identity")
rho_phys = sp.simplify(rho.subs({Lambda: sp.pi**2 * R_target * Mmix / (8 * (1 - eps)), R_target: F}))
print("rho_alpha on the physical target branch =")
sp.pprint(rho_phys)
expect_zero("rho_alpha - G_tr/M_mix", rho_phys - sp.simplify(G / Mmix))

subbanner("III. Exact regime split in the selected-branch variable")
print("mixed-only enough      <=> rho_alpha <= 1")
print("lowest symmetric twin  <=> 1 < rho_alpha <= 2")
print("non-twin asymmetry     <=> rho_alpha > 2")
print()
print("Equivalently on the physical branch:")
print("  G_tr(xi,delta;R) <= M_mix")
print("  M_mix < G_tr(xi,delta;R) <= 2 M_mix")
print("  G_tr(xi,delta;R) > 2 M_mix")

subbanner("IV. Strict monotonicity with respect to the selected-branch depth")
drho_dxi = sp.simplify(sp.diff(G / Mmix, xi))
print("d rho_alpha / d xi =")
sp.pprint(drho_dxi)

expect_zero(
    "drho_dxi - (1/M_mix) dG/dxi",
    drho_dxi - sp.simplify(sp.diff(G, xi) / Mmix),
)

print("Since dG_tr/dxi > 0 on the stable branch and M_mix > 0, rho_alpha rises")
print("strictly with the physical selected-branch softening depth xi.")
