#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage346_349_common import *

banner("STAGE 346 — ACTUAL COHERENT TRACKING-BRANCH VALUES")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
epsW, epsEta = sp.symbols("epsilon_W epsilon_eta", positive=True, real=True)
ZW, Lambda = sp.symbols("Z_W Lambda", positive=True, real=True)
zeta = sp.symbols("zeta", positive=True, real=True)
chiQ = sp.symbols("chi_Q", positive=True, real=True)

eps = sp.simplify(epsW * (1 - sp.Rational(2,11) * deltaU / (1 + deltaU)))
Rtr = sp.simplify((1 + chi0/(1 + deltaU)) / (1 + chi0))
Rtarget = sp.simplify(Lambda * (1 - epsEta) * (1 - eps)**2 / (ZW * (1 + chi0)**2))
Mmix = sp.simplify(8 * ZW * (1 + chi0)**2 / (sp.pi**2 * (1 - epsEta) * (1 - eps)))
S = sp.simplify(1 + zeta*(1 - eps)/(1 - zeta*eps))
Mtr = sp.simplify(Mmix * S)
NQ = sp.simplify(1/chiQ)

subbanner("I. Actual coherent local D/N branch values")
print("R_tr =")
sp.pprint(Rtr)
print("epsilon =")
sp.pprint(eps)
print("R_target =")
sp.pprint(Rtarget)
print("epsilon_eta =")
sp.pprint(epsEta)
print("M_mix =")
sp.pprint(Mmix)
print("S(zeta;epsilon) =")
sp.pprint(S)
print("M_tr =")
sp.pprint(Mtr)
print("N_Q on the natural source-map branch =")
sp.pprint(NQ)

subbanner("II. Exact support-blindness of the orbit packet")
expect_zero("d R_tr / d zeta", sp.diff(Rtr, zeta))
expect_zero("d epsilon_eta / d zeta", sp.diff(epsEta, zeta))
expect_zero("d R_target / d zeta", sp.diff(Rtarget, zeta))

subbanner("III. Exact product law on the physical coherent branch")
product_law = sp.simplify(Rtarget * Mmix)
print("R_target * M_mix =")
sp.pprint(product_law)
expect_zero(
    "product law residual",
    product_law - sp.simplify(8 * Lambda * (1 - eps) / sp.pi**2),
)

print("\nInterpretation:")
print("  The actual coherent branch packet is carried by (R_tr, R_target, epsilon_eta),")
print("  while the support lane enters only through the baseline enhancement M_tr = M_mix S.")
print("  So zeta can move the steady normalization baseline but cannot move the orbit packet.")
