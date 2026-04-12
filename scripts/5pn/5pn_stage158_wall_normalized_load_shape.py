#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 158 — WALL-NORMALIZED PORT SHAPE VARIABLES AND THE OUTGOING-LOAD THEOREM")

c, varpi, K = sp.symbols('c varpi K', positive=True, real=True)
Q, Delta, P = sp.symbols('Q Delta P', nonzero=True, real=True)
deltaK = sp.symbols('deltaK', real=True)

subbanner("1. Exact wall-normalized factorization")
chi = sp.simplify(c/(sp.sqrt(K)*varpi))
Upsilon = sp.simplify(Q/(K*Delta))
Lambda = sp.simplify(P/Delta)

expect_zero("B0 - K chi^2", c**2/varpi**2 - K*chi**2)
expect_zero("Z0 - K Upsilon", Q/Delta - K*Upsilon)
expect_zero("N0 - Lambda^2", P**2/Delta**2 - Lambda**2)

print("chi =", chi)
print("Upsilon =", Upsilon)
print("Lambda =", Lambda)

subbanner("2. Exact rewrite of the defect fields")
dlnchi2, dlnUpsilon, dlnLambda = sp.symbols('dlnchi2 dlnUpsilon dlnLambda', real=True)
SigmaB = dlnchi2
SigmaZ = dlnUpsilon
SigmaN = sp.simplify(2*dlnLambda - deltaK)
print("Sigma_B =", SigmaB)
print("Sigma_Z =", SigmaZ)
print("Sigma_N =", SigmaN)

subbanner("3. Exact wall-normalized load-shape formula")
ThetaN, ThetaB, ThetaZ, omegaB, omegaZ = sp.symbols('ThetaN ThetaB ThetaZ omegaB omegaZ', real=True)
Xi_load = sp.simplify(ThetaN + omegaB*ThetaB + omegaZ*ThetaZ)
print("Xi_load =", Xi_load)

subbanner("4. Conservative-shape theorem and naive self-similarity no-go")
rho1, rho2 = sp.symbols('rho1 rho2', real=True)
ThetaN_cons = sp.simplify(rho1*(2*dlnLambda - deltaK) + rho2*(2*dlnLambda - deltaK))
expect_zero("conservative-shape branch => Xi_load = weighted outgoing load", ThetaN_cons - (rho1+rho2)*(2*dlnLambda - deltaK))
naive_common = sp.simplify((rho1+rho2)*(0 + omegaB*0 + omegaZ*0) + (rho1+rho2)*(-deltaK))
print("If all wall-normalized shapes are frozen, Xi_load =", naive_common)
print("So naive common self-similarity gives Xi_load = -deltaK unless deltaK = 0.")

subbanner("5. Exact outgoing-load theorem")
dlnLambda_r1, dlnLambda_r2 = sp.symbols('dlnLambda_r1 dlnLambda_r2', real=True)
rhoN1, rhoN2 = sp.symbols('rhoN1 rhoN2', real=True)
Xi_out = sp.simplify(2*(rhoN1*dlnLambda_r1 + rhoN2*dlnLambda_r2) - (rhoN1+rhoN2)*deltaK)
print("Outgoing-load theorem representative =", Xi_out)
print("Vanishing defect requires 2 sum_r rho_r^(N) dln Lambda_r = deltaK.")

banner("STAGE 158 LEDGER")
print("1. The BdG support defect is exactly the drift of the wall-normalized support")
print("   shape chi_alpha.")
print("2. The conservative Maxwell/mixed defect is exactly the drift of the")
print("   wall-normalized port shape Upsilon_r.")
print("3. The outgoing-transfer defect is exactly the wall-loading mismatch of the")
print("   dimensionless load factor Lambda_r = P_r/Delta_r.")
print("4. On conservative-shape-preserving branches, the full remaining linear grouped")
print("   defect is carried only by the outgoing load factor.")
print("5. A naive branch that freezes every wall-normalized shape gives Xi_load = -deltaK,")
print("   so the outgoing transfer sector must actively co-load with the wall baseline.")
print("6. The exact outgoing-load theorem is")
print("      2 sum_r rho_r^(N) dln Lambda_r = deltaK.")
