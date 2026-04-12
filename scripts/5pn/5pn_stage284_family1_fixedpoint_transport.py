#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage284_288_common import *

"""
Stage 284 — exact co-evolving Family-1 fixed-point transport.

What this script does
---------------------
1. Re-derives the linear transport from the renormalized Family-1 canonical point
   into the mouth/core load ledger (delta M_s, delta M_q, delta Pi).
2. Specializes to the canonical-even tangent delta g = 0 and extracts the exact
   tangential defect variable delta Pi_tan.
3. Verifies the numerical coefficients recorded in the note.
"""

banner("STAGE 284 — FAMILY-1 FIXED-POINT TRANSPORT")

refs = family1_reference_values()
Sigma0_can = refs["Sigma0_can"]
Tmhat_can = refs["Tmhat_can"]
g_star = refs["g_star"]
R_can = refs["R_can"]
S_can = refs["S_can"]
r_F1 = refs["r_F1"]

Sigma0, S, g, r = sp.symbols("Sigma_0 S g r", positive=True, real=True)
dSigma0, dS, dg, dTmhat = sp.symbols("deltaSigma0 deltaS deltag deltaTmhat", real=True)

subbanner("I. Exact Family-1 gain and slope relations")
R = R_of_g(g, r)
Ms = Sigma0
Mq = sp.simplify(-Sigma0 * R)
Pi = sp.simplify(Sigma0 * (1 - R * S))

print("R(g,r) =")
sp.pprint(R)
print("M_s =")
sp.pprint(Ms)
print("M_q =")
sp.pprint(Mq)
print("Pi =")
sp.pprint(Pi)

subbanner("II. Lower compensated branch transport: delta g -> delta R")
g_lower = g_minus(r)
dR_lower = sp.simplify(sp.diff(R, g).subs({g: g_lower}) * dg)
print("g_-(r) =")
sp.pprint(g_lower)
print("delta R on lower branch =")
sp.pprint(dR_lower)
expect_zero(
    "deltaR + dg/sqrt(1+r^2)",
    sp.simplify(dR_lower + dg / sp.sqrt(1 + r**2)),
)

subbanner("III. Exact load transport about the canonical point")
# Linearized mouth/core transport at fixed r, with the lower-branch deltaR law inserted.
dMq = sp.simplify(sp.diff(Mq, Sigma0) * dSigma0 + sp.diff(Mq, g).subs({g: g_lower}) * dg)
dPi = sp.simplify(
    sp.diff(Pi, Sigma0) * dSigma0
    + sp.diff(Pi, S) * dS
    + sp.diff(Pi, g).subs({g: g_lower}) * dg
)

print("delta M_s =")
sp.pprint(dSigma0)
print("delta M_q =")
sp.pprint(dMq)
print("delta Pi =")
sp.pprint(dPi)

# Canonical-even tangent is delta g = 0.
dPi_tan = sp.simplify(dPi.subs({dg: 0, g: g_star, Sigma0: Sigma0_can, S: S_can, r: r_F1}))
print("delta Pi_tan =")
sp.pprint(dPi_tan)

# Relation Sigma0 = 20/9 * Tmhat^2
Tmhat = sp.symbols('Tmhat', positive=True, real=True)
Sigma0_from_T = sp.Rational(20, 9) * Tmhat ** 2
print("Sigma0(Tmhat) =")
sp.pprint(Sigma0_from_T)

dSigma0_from_T = sp.simplify(sp.diff(Sigma0_from_T, Tmhat).subs({Tmhat: Tmhat_can}) * dTmhat)
print("delta Sigma0 from delta Tmhat =")
sp.pprint(dSigma0_from_T)

subbanner("IV. Numerical coefficients at the renormalized canonical point")
coef_dSigma0 = sp.N(sp.diff(dPi_tan, dSigma0))
coef_dS = sp.N(sp.diff(dPi_tan, dS))
dPi_tan_T = sp.expand(dPi_tan.subs({dSigma0: dSigma0_from_T / dTmhat * dTmhat}))
coef_dTmhat = sp.N(sp.diff(dPi_tan_T, dTmhat))

print(f"coefficient of deltaSigma0  = {coef_dSigma0}")
print(f"coefficient of deltaS       = {coef_dS}")
print(f"coefficient of deltaTmhat   = {coef_dTmhat}")

for name, value, target in [
    ("deltaPi_tan numerical dSigma0 coefficient", coef_dSigma0, sp.Float("0.832409471081635")),
    ("deltaPi_tan numerical dS coefficient", coef_dS, sp.Float("-1.16275838754222")),
    ("deltaPi_tan numerical dTmhat coefficient", coef_dTmhat, sp.Float("5.35223887169622")),
]:
    diff = sp.N(value - target)
    print(f"{name} residual = {diff}")
    if abs(float(diff)) > 1e-12:
        raise AssertionError(f"{name} mismatch")

subbanner("V. Final compact transport law")
print("On the canonical-even tangent delta g = 0,")
print("delta Pi_tan =")
sp.pprint(dPi_tan)
print()
print("Equivalently, using delta Sigma0 = (40/9) Tmhat_can delta Tmhat,")
print("delta Pi_tan =")
sp.pprint(dPi_tan_T)
