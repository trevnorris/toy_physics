#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 153 — EXACT LINEAR GROUPED-P2 MAP INTO DIRECT OUTLET COEFFICIENTS")

D0, N0, sigma_star = sp.symbols('D0 N0 sigma_star', nonzero=True, real=True)
u2 = sp.Rational(1, 9)
P0 = N0/D0

dD0, dD2, dD4, dN0 = sp.symbols('dD0 dD2 dD4 dN0', real=True)
mathcalK = sp.simplify(dD2 + dD0/9)
mathcalG = sp.simplify(dN0 - P0*dD0)

subbanner("1. Physical grouped response variations")
du2 = sp.simplify(-(dD2 + u2*dD0)/D0)
du4 = sp.simplify(-(5*dD0 + 18*dD2 + 81*dD4)/(81*D0))
dP0 = sp.simplify((dN0 - P0*dD0)/D0)

expect_zero("du2 + mathcalK/D0", du2 + mathcalK/D0)
expect_zero("dP0 - mathcalG/D0", dP0 - mathcalG/D0)

print("delta u2 =", du2)
print("delta u4 =", du4)
print("delta P0 =", dP0)

subbanner("2. Direct outlet map")
delta_kappa, delta_gamma = sp.symbols('delta_kappa delta_gamma', real=True)
# Stage-142/153 direct formulas on the pure linear grouped-anisotropy branch.
eq_kappa = sp.Eq(du2, -sigma_star*delta_kappa/(3*(1-sigma_star)))
eq_gamma = sp.Eq(dP0/P0, -9*sigma_star*delta_gamma/(1-sigma_star))
sol_kappa = sp.solve(eq_kappa, delta_kappa)[0]
sol_gamma = sp.solve(eq_gamma, delta_gamma)[0]

expect_zero("delta kappa formula", sol_kappa - 3*(1-sigma_star)*mathcalK/(sigma_star*D0))
expect_zero("delta gamma formula", sol_gamma + (1-sigma_star)*mathcalG/(9*sigma_star*N0))

print("delta kappa_W =", sol_kappa)
print("delta gamma_W =", sol_gamma)

subbanner("3. One-parameter hidden-even consistency")
hidden_even_residual = sp.simplify(du4 - sp.Rational(8,9)*du2)
expect_zero("hidden-even residual + (dD4 - (2 dD2)/3 - dD0/27)/D0", hidden_even_residual + (dD4 - 2*dD2/3 - dD0/27)/D0)
print("hidden-even residual =", hidden_even_residual)
print("Equivalent operator law: dD4 = 2 dD2/3 + dD0/27")

subbanner("4. Weak-axisymmetric one-parameter collapse")
eps = sp.symbols('epsilon', real=True)
lam20, lam21, lam22 = sp.Integer(1), sp.Rational(1,2), sp.Integer(-1)
K1, G1 = sp.symbols('K1 G1', real=True)

for lbl, lam in [('20',lam20), ('21',lam21), ('22',lam22)]:
    print(f"mathcalK_{lbl} =", sp.simplify(eps*lam*K1))
    print(f"mathcalG_{lbl} =", sp.simplify(eps*lam*G1))

kappa1 = sp.simplify(3*(1-sigma_star)*K1/(sigma_star*D0))
gamma1 = sp.simplify(-(1-sigma_star)*G1/(9*sigma_star*N0))
print("kappa_1 =", kappa1)
print("gamma_1 =", gamma1)

banner("STAGE 153 LEDGER")
print("1. The full linear grouped-P2 outlet problem collapses to two exact microscopic")
print("   combinations:")
print("      mathcalK_A = dD_{A,2} + dD_{A,0}/9,")
print("      mathcalG_A = dN_{A,0} - P0 dD_{A,0}.")
print("2. The direct outlet deformations are")
print("      delta kappa_W^(A) = 3(1-sigma_*) mathcalK_A / (sigma_* D0),")
print("      delta gamma_W^(A) = -(1-sigma_*) mathcalG_A / (9 sigma_* N0).")
print("3. The one-parameter hidden-even consistency condition is")
print("      dD4 = (2/3) dD2 + dD0/27,")
print("   equivalently delta u4 = (8/9) delta u2.")
print("4. So the remaining linear grouped-anisotropy bottleneck is not the whole")
print("   bundle of dD and dN data, but just the pair (mathcalK_A, mathcalG_A).")
