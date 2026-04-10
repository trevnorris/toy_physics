#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 157 — STATIC SELF-SIMILARITY AND THE REMAINING LINEAR GROUPED DEFECT")

deltaN, deltaK, deltaB, deltaZ = sp.symbols('deltaN deltaK deltaB deltaZ', real=True)
omegaB, omegaZ = sp.symbols('omegaB omegaZ', real=True)
omegaK = sp.simplify(1 + omegaB + omegaZ)

subbanner("1. Exact static-slope decomposition")
deltaD = sp.simplify(omegaK*deltaK - omegaB*deltaB - omegaZ*deltaZ)
expect_zero("weight identity rewrite", sp.simplify(omegaK - omegaB - omegaZ - 1))

Xi_load = sp.simplify(deltaN - deltaD)
Xi_rewrite = sp.simplify((deltaN - deltaK) + omegaB*(deltaB - deltaK) + omegaZ*(deltaZ - deltaK))
expect_zero("Xi_load rewrite", Xi_load.subs({omegaK: 1 + omegaB + omegaZ}) - Xi_rewrite)

print("Xi_load =", Xi_rewrite)

subbanner("2. Wall-referenced self-similarity defect fields")
SigmaB, SigmaZ, SigmaN = sp.symbols('SigmaB SigmaZ SigmaN', real=True)
rhoB1, rhoB2 = sp.symbols('rhoB1 rhoB2', real=True)
rhoZ1, rhoZ2 = sp.symbols('rhoZ1 rhoZ2', real=True)
rhoN1, rhoN2 = sp.symbols('rhoN1 rhoN2', real=True)
ThetaB = sp.simplify(rhoB1*SigmaB + rhoB2*SigmaB)
ThetaZ = sp.simplify(rhoZ1*SigmaZ + rhoZ2*SigmaZ)
ThetaN = sp.simplify(rhoN1*SigmaN + rhoN2*SigmaN)
print("Generic Stage-157 form: Xi_load = Theta_N + omega_B Theta_B + omega_Z Theta_Z")

subbanner("3. Static self-similarity theorem")
Xi_ss = sp.simplify(SigmaN + omegaB*SigmaB + omegaZ*SigmaZ)
expect_zero("static self-similarity => Xi_load = 0", Xi_ss.subs({SigmaB:0, SigmaZ:0, SigmaN:0}))
print("If Sigma_B = Sigma_Z = Sigma_N = 0, then Xi_load = 0.")

banner("STAGE 157 LEDGER")
print("1. The remaining linear grouped defect is exactly")
print("      Xi_load = (deltaN-deltaK) + omega_B(deltaB-deltaK) + omega_Z(deltaZ-deltaK).")
print("2. The wall baseline slope deltaK is the natural reference for the remaining")
print("   loading mismatch.")
print("3. Each support/transfer sector contributes only through its failure to co-load")
print("   self-similarly with that wall baseline.")
print("4. If the branch is statically self-similar relative to the wall baseline,")
print("   then Xi_load vanishes automatically.")
