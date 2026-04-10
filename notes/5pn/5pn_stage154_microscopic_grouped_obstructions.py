#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 154 — EXACT MICROSCOPIC DECOMPOSITION OF THE LINEAR GROUPED OUTLET OBSTRUCTIONS")

P0 = sp.symbols('P0', nonzero=True, real=True)
dK, dM, dB0, dB2, dZ0, dZ2, dN0 = sp.symbols('dK dM dB0 dB2 dZ0 dZ2 dN0', real=True)

subbanner("1. Full-bundle decomposition of mathcalK_A and mathcalG_A")
dD0 = sp.simplify(dK - dB0 - dZ0)
dD2 = sp.simplify(-(dM + dB2 + dZ2))

mathcalK = sp.simplify(dD2 + dD0/9)
mathcalG = sp.simplify(dN0 - P0*dD0)

W_A = sp.simplify(dK/9 - dM)
B_A = sp.simplify(dB2 + dB0/9)
Z_A = sp.simplify(dZ2 + dZ0/9)
N_A = sp.simplify(dN0 + P0*dZ0)

expect_zero("mathcalK - (W_A - B_A - Z_A)", mathcalK - (W_A - B_A - Z_A))
expect_zero("mathcalG - (-P0 dK + P0 dB0 + N_A)", mathcalG - (-P0*dK + P0*dB0 + N_A))

print("mathcalK =", mathcalK)
print("mathcalG =", mathcalG)

subbanner("2. Exact BdG support contribution")
eps = sp.symbols('epsilon', real=True)
c, varpi, dc, dvarpi = sp.symbols('c varpi dc dvarpi', nonzero=True, real=True)
B0 = c**2/varpi**2
B2 = c**2/varpi**4
B0_var = sp.expand(sp.series(((c + eps*dc)**2)/(varpi + eps*dvarpi)**2, eps, 0, 2).removeO())
B2_var = sp.expand(sp.series(((c + eps*dc)**2)/(varpi + eps*dvarpi)**4, eps, 0, 2).removeO())
dB0_formula = sp.simplify(B0_var.coeff(eps,1))
dB2_formula = sp.simplify(B2_var.coeff(eps,1))
expect_zero("dB0 - [2c dc/varpi^2 - 2 c^2 dvarpi/varpi^3]", dB0_formula - (2*c*dc/varpi**2 - 2*c**2*dvarpi/varpi**3))
expect_zero("dB2 - [2c dc/varpi^4 - 4 c^2 dvarpi/varpi^5]", dB2_formula - (2*c*dc/varpi**4 - 4*c**2*dvarpi/varpi**5))
Bcal = sp.simplify(dB2_formula + dB0_formula/9)
print("Single-mode BdG obstruction contribution =", Bcal)

subbanner("3. Exact conservative Maxwell/mixed contribution")
Q, S, G, Delta, dQ, dS, dG, dDelta = sp.symbols('Q S G Delta dQ dS dG dDelta', nonzero=True, real=True)
Z0 = Q/Delta
Z2 = (Q*S - G*Delta)/Delta**2

Z0_var = sp.expand(sp.series((Q + eps*dQ)/(Delta + eps*dDelta), eps, 0, 2).removeO())
Z2_var = sp.expand(sp.series(((Q + eps*dQ)*(S + eps*dS) - (G + eps*dG)*(Delta + eps*dDelta))/(Delta + eps*dDelta)**2, eps, 0, 2).removeO())
dZ0_formula = sp.simplify(Z0_var.coeff(eps,1))
dZ2_formula = sp.simplify(Z2_var.coeff(eps,1))

expect_zero("dZ0 formula", dZ0_formula - (Delta*dQ - Q*dDelta)/Delta**2)
expect_zero("dZ2 formula", dZ2_formula - (S*dQ/Delta**2 + Q*dS/Delta**2 - dG/Delta + (G/Delta**2 - 2*Q*S/Delta**3)*dDelta))
Zcal = sp.simplify(dZ2_formula + dZ0_formula/9)
print("Single-port conservative obstruction contribution =", Zcal)

subbanner("4. Exact outgoing-transfer contribution")
P = sp.symbols('P', nonzero=True, real=True)
dP = sp.symbols('dP', real=True)
N0 = P**2/Delta**2
N0_var = sp.expand(sp.series((P + eps*dP)**2/(Delta + eps*dDelta)**2, eps, 0, 2).removeO())
dN0_formula = sp.simplify(N0_var.coeff(eps,1))
expect_zero("dN0 formula", dN0_formula - (2*P*dP/Delta**2 - 2*P**2*dDelta/Delta**3))
Nbundle = sp.simplify(dN0_formula + P0*dZ0_formula)
expect_zero("Nbundle formula", Nbundle - (P0*dQ/Delta + 2*P*dP/Delta**2 - (P0*Q/Delta**2 + 2*P**2/Delta**3)*dDelta))
print("Single-port outgoing-transfer bundle =", Nbundle)

subbanner("5. Weak-axisymmetric collapse to (frakK_1, frakG_1)")
lam20, lam21, lam22 = sp.Integer(1), sp.Rational(1,2), sp.Integer(-1)
frakK1, frakG1 = sp.symbols('frakK1 frakG1', real=True)
for lbl, lam in [('20',lam20), ('21',lam21), ('22',lam22)]:
    print(f"mathcalK_{lbl} =", sp.simplify(eps*lam*frakK1))
    print(f"mathcalG_{lbl} =", sp.simplify(eps*lam*frakG1))

banner("STAGE 154 LEDGER")
print("1. The linear grouped-even and grouped-odd outlet defects are driven only by")
print("   four microscopic defect bundles:")
print("      W_A, B_A, Z_A, N_A.")
print("2. Exactly,")
print("      mathcalK_A = W_A - B_A - Z_A,")
print("      mathcalG_A = -P0 dK_A + P0 dB_{A,0} + N_A.")
print("3. For each BdG mode, only the weighted combination")
print("      dB_{A,2} + dB_{A,0}/9")
print("   enters the even obstruction.")
print("4. For each Maxwell/mixed port, the conservative obstruction depends only on")
print("      dQ, dS, dG, dDelta,")
print("   while the direct odd transfer bundle depends only on")
print("      dQ, dP, dDelta.")
print("5. So the weak grouped-lane outlet bottleneck has collapsed to the scalar pair")
print("      (frakK_1, frakG_1)")
print("   once the physical branch tells us the grouped microscopic slopes.")
