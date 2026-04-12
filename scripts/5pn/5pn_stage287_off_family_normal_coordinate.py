#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage284_288_common import *

"""
Stage 287 — off-family normal coordinate and microscopic compensation defect.

What this script does
---------------------
1. Defines the exact off-family normal coordinate delta_perp.
2. Proves that it is the normalized first-order compensation defect.
3. Derives the exact microscopic formula for delta_perp in the parent ratios
   (K_s, K_q, lambda, g_s, g_q).
4. Shows how delta R_q is transported by the same scalar.
"""

banner("STAGE 287 — OFF-FAMILY NORMAL COORDINATE")

r = sp.symbols("r", positive=True, real=True)
dg, dr = sp.symbols("deltag deltar", real=True)

dlgq, dlgs, dllam, dlKs, dlKq = sp.symbols(
    "dln_gq dln_gs dln_lambda dln_Ks dln_Kq", real=True
)

subbanner("I. Exact parent compensation function and normal coordinate")
gm = g_minus(r)
gprime = g_minus_prime(r)
F = sp.simplify(1 + r**2 - 4 * (sp.Symbol("g") - r) ** 2)

delta_perp = sp.simplify(dg - gprime * dr)
print("g_-(r) =")
sp.pprint(gm)
print("delta_perp =")
sp.pprint(delta_perp)

subbanner("II. delta F and delta R_q are exact images of delta_perp")
dF = sp.simplify(2*r*dr - 8*(gm - r)*(dg - dr))
print("delta F on the lower branch =")
sp.pprint(dF)
expect_zero("deltaF - 4 sqrt(1+r^2) delta_perp", sp.simplify(dF - 4*sp.sqrt(1+r**2)*delta_perp))

Rq = R_of_g(sp.Symbol("g"), r)
dRq = sp.simplify(sp.diff(Rq, sp.Symbol("g")).subs({sp.Symbol("g"): gm}) * dg + sp.diff(Rq, r).subs({sp.Symbol("g"): gm}) * dr)
print("delta R_q on the lower branch =")
sp.pprint(dRq)
expect_zero("deltaRq + delta_perp/sqrt(1+r^2)", sp.simplify(dRq + delta_perp/sp.sqrt(1+r**2)))

subbanner("III. Exact microscopic parent-variable formula")
dr_parent = sp.simplify(r * (dllam - sp.Rational(1, 2) * dlKs - sp.Rational(1, 2) * dlKq))
dg_parent = sp.simplify(gm * (dlgq - dlgs + sp.Rational(1, 2) * dlKs - sp.Rational(1, 2) * dlKq))
delta_perp_parent = sp.simplify(dg_parent - gprime * dr_parent)

log_imbalance_1 = sp.simplify(dlgq + dlKs - dlgs - dllam)
log_imbalance_2 = sp.simplify(dlKs + dlKq - 2 * dllam)
delta_perp_expected = sp.simplify(gm * log_imbalance_1 + log_imbalance_2 / (4 * sp.sqrt(1 + r**2)))

print("delta_perp from parent drifts =")
sp.pprint(delta_perp_parent)
print("Expected compact form =")
sp.pprint(delta_perp_expected)
expect_zero("delta_perp parent formula", sp.simplify(delta_perp_parent - delta_perp_expected))

subbanner("IV. Final off-family transport laws")
print("delta_perp =")
sp.pprint(delta_perp_expected)
print("delta F =")
sp.pprint(4 * sp.sqrt(1 + r**2) * delta_perp_expected)
print("delta R_q =")
sp.pprint(-delta_perp_expected / sp.sqrt(1 + r**2))
