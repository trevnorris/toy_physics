#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage289_292_common import *

"""
Stage 293 — exact linear grouped outlet map into the direct outlet coefficients
(delta kappa_W, delta gamma_W).
"""

banner("STAGE 293 — LINEAR GROUPED OUTLET MAP")

D0, N0, P0, sigma = sp.symbols("D_0 N_0 P_0 sigma_star", positive=True, real=True)
dD0, dD2, dD4, dN0 = sp.symbols("deltaD0 deltaD2 deltaD4 deltaN0", real=True)

subbanner("I. Exact grouped conservative / outgoing slopes")
du2 = sp.simplify(-(dD2 + dD0 / 9) / D0)
du4 = sp.simplify(-(dD4 + sp.Rational(2, 9) * dD2 + sp.Rational(5, 81) * dD0) / D0)
dP0 = sp.simplify((dN0 - P0 * dD0) / D0)

print("delta u2 =")
sp.pprint(du2)
print("delta u4 =")
sp.pprint(du4)
print("delta P0 =")
sp.pprint(dP0)

subbanner("II. Direct outlet coefficients")
dkappaW = sp.simplify(-3 * (1 - sigma) * du2 / sigma)
dgammaW = sp.simplify(-(1 - sigma) * dP0 / (9 * sigma * P0))
print("delta kappa_W =")
sp.pprint(dkappaW)
print("delta gamma_W =")
sp.pprint(dgammaW)
expect_zero(
    "delta kappa_W transport",
    sp.simplify(dkappaW - 3 * (1 - sigma) * (dD2 + dD0 / 9) / (sigma * D0)),
)
expect_zero(
    "delta gamma_W transport",
    sp.simplify(dgammaW + (1 - sigma) * (dN0 - P0 * dD0) / (9 * sigma * P0 * D0)),
)

subbanner("III. Exact one-parameter even consistency relation")
constraint = sp.simplify(du4 - sp.Rational(8, 9) * du2)
print("delta u4 - 8 delta u2 / 9 =")
sp.pprint(constraint)
expect_zero(
    "even consistency iff dD4 = 2 dD2/3 + dD0/27",
    sp.simplify(constraint.subs({dD4: sp.Rational(2, 3) * dD2 + dD0 / 27})),
)

constraint_solved = sp.solve(sp.Eq(constraint, 0), dD4)[0]
print("Equivalent microscopic even consistency relation: delta D4 =")
sp.pprint(constraint_solved)
expect_zero(
    "solved even consistency relation",
    sp.simplify(constraint_solved - (sp.Rational(2, 3) * dD2 + dD0 / 27)),
)

subbanner("IV. Interpretation")
print("The linear grouped-even problem collapses to")
print("  dD2 + dD0/9")
print("plus the exact hidden-even compatibility")
print("  dD4 = 2 dD2/3 + dD0/27.")
print("The linear grouped-odd problem collapses to")
print("  dN0 - P0 dD0.")
