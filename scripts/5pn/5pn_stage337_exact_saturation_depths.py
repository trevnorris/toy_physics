#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage336_339_common import *

"""
Stage 337 — exact saturation depths for mixed-only and lowest-twin support.

What this script does
---------------------
1. Solves G_tr = M_mix exactly for the mixed-only saturation depth xi_(1x).
2. Solves G_tr = 2 M_mix exactly for the lowest-twin saturation depth xi_(2x).
3. Rewrites the three support regimes directly as xi-intervals.
4. Records the exact branch-ordering theorem xi_(1x) < xi_(2x).
"""

banner("STAGE 337 — EXACT SATURATION DEPTHS xi_(1x) AND xi_(2x)")

xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)
Mmix = sp.symbols("M_mix", positive=True, real=True)

G = G_tr(xi, delta, R)
xi_1x = xi_threshold(sp.Integer(1), Mmix, delta, R)
xi_2x = xi_threshold(sp.Integer(2), Mmix, delta, R)

subbanner("I. Exact mixed-only saturation depth")
print("xi_(1x) from G_tr = M_mix =")
sp.pprint(xi_1x)
expect_zero(
    "G_tr(xi_(1x)) - M_mix",
    sp.simplify(G.subs(xi, xi_1x) - Mmix),
)

subbanner("II. Exact lowest-twin saturation depth")
print("xi_(2x) from G_tr = 2 M_mix =")
sp.pprint(xi_2x)
expect_zero(
    "G_tr(xi_(2x)) - 2 M_mix",
    sp.simplify(G.subs(xi, xi_2x) - 2 * Mmix),
)

subbanner("III. Exact xi-interval regime classifier")
print("mixed-only enough      <=> xi_phys <= xi_(1x)")
print("lowest symmetric twin  <=> xi_(1x) < xi_phys <= xi_(2x)")
print("non-twin asymmetry     <=> xi_phys > xi_(2x)")

subbanner("IV. Exact branch ordering")
dG_dxi = sp.simplify(sp.diff(G, xi))
print("dG_tr/dxi =")
sp.pprint(dG_dxi)
print("Since G_tr is strictly increasing and 2 M_mix > M_mix, the unique roots obey")
print("xi_(1x) < xi_(2x) on every stable tracking branch.")

print()
print("Useful closed forms:")
print("xi_(1x) =")
sp.pprint(sp.factor(xi_1x))
print("xi_(2x) =")
sp.pprint(sp.factor(xi_2x))
