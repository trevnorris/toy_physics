#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage284_288_common import *

"""
Stage 286 — parent compensation-family rigidity and automatic similarity preservation.

What this script does
---------------------
1. Uses the exact parent compensation family to prove D/N similarity preservation
   identically along the family.
2. Shows that on the lower compensated branch the canonical-even condition
   delta g = 0 forces delta r = 0.
3. Concludes that the first-order reduced outgoing defect vanishes automatically
   on-family.
"""

banner("STAGE 286 — PARENT COMPENSATION RIGIDITY")

r = sp.symbols("r", positive=True, real=True)
dg, dr = sp.symbols("deltag deltar", real=True)
sigma_star, dPi_tan = sp.symbols("sigma_star deltaPi_tan", real=True)

subbanner("I. Exact parent compensation family")
gm = g_minus(r)
Lratio = sp.simplify(sp.pi / 2 * sp.sqrt((1 + r**2) / 3))
gamma0 = sp.simplify((1 + r**2) / 9)

print("g_-(r) =")
sp.pprint(gm)
print("L_W/a =")
sp.pprint(Lratio)
print("gamma_0 =")
sp.pprint(gamma0)

subbanner("II. Automatic D/N similarity preservation on the parent family")
Xi_slip_family = sp.simplify(sp.diff(sp.log(gamma0), r) - 2 * sp.diff(sp.log(Lratio), r))
print("d/dr [ ln gamma0 - 2 ln(L_W/a) ] =")
sp.pprint(Xi_slip_family)
expect_zero("parent-family Xi_slip", Xi_slip_family)

subbanner("III. Lower compensated branch rigidity under delta g = 0")
gprime = sp.simplify(sp.diff(gm, r))
print("g_-'(r) =")
sp.pprint(gprime)
print("Alternative manifestly positive form =")
positive_form = sp.simplify((4 + 3*r**2) / (2*sp.sqrt(1 + r**2) * (2*sp.sqrt(1 + r**2) + r)))
sp.pprint(positive_form)
expect_zero("g_-'(r) - positive form", sp.simplify(gprime - positive_form))

print("Therefore delta g = g_-'(r_*) delta r, and on the lower branch delta g = 0 implies delta r = 0.")

subbanner("IV. Collapse of the first-order defect ledger on-family")
# Along the exact parent family, with delta r = 0, every first-order similarity drift vanishes.
dln_gamma0 = sp.simplify(sp.diff(sp.log(gamma0), r) * dr)
dln_Lratio = sp.simplify(sp.diff(sp.log(Lratio), r) * dr)
print("delta ln gamma0 =")
sp.pprint(dln_gamma0)
print("delta ln(L_W/a) =")
sp.pprint(dln_Lratio)
expect_zero("delta ln gamma0 - 2 delta ln(L_W/a)", sp.simplify(dln_gamma0 - 2*dln_Lratio))

Xi_slip = sp.symbols("Xi_slip", real=True)
DeltaQ = sp.simplify(-sigma_star * Xi_slip * dPi_tan / (1 - sigma_star))
print("If the branch stays on the exact parent family and delta g = 0, then Xi_slip = 0 and")
print("Delta_Q =")
sp.pprint(DeltaQ.subs({Xi_slip: 0}))

subbanner("V. Numerical rigidity at the Family-1 value")
r_F1 = family1_reference_values()["r_F1"]
gprime_num = sp.N(gprime.subs({r: r_F1}))
print(f"g_-'(r_F1) = {gprime_num}")
print(f"So delta r ≈ {sp.N(1/gprime_num)} delta g on the lower branch.")
