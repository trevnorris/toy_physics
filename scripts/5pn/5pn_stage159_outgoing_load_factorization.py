#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 159 — OUTGOING LOAD-FACTOR FACTORIZATION AND THE SQUARE-ROOT MIXED-LEG LAW")

GW, GU, R, OU2, OW2, K = sp.symbols('G_W G_U R O_U2 O_W2 K', positive=True, real=True)
dlnGW, dlnGU, dlnR, dlnOU2, dlnOW2, deltaK = sp.symbols('dlnGW dlnGU dlnR dlnOU2 dlnOW2 deltaK', real=True)

subbanner("1. Exact factorization of the outgoing load")
P = OU2*GW + R*GU
Delta = OU2*OW2 - R**2
Lambda = sp.simplify(P/Delta)

M = sp.simplify(GW/(OW2*sp.sqrt(K)))
I = sp.simplify(R*GU/(OU2*GW))
H = sp.simplify(R**2/(OU2*OW2))

expect_zero("Lambda^2/K factorization", sp.simplify(Lambda**2/K - M**2*(1+I)**2/(1-H)**2))
print("M_r =", M)
print("I_r =", I)
print("H_r =", H)

subbanner("2. Exact decomposition of the outgoing defect field")
SigmaN = sp.simplify(2*sp.Symbol('dlnM') + 2*sp.Symbol('dln1pI') - 2*sp.Symbol('dln1mH'))
print("Sigma_N = 2 dln M + 2 dln(1+I) - 2 dln(1-H)")

subbanner("3. First-order transport of the three microscopic slippages")
dlnM = sp.simplify(dlnGW - dlnOW2 - deltaK/2)
dlnI = sp.simplify(dlnR + dlnGU - dlnOU2 - dlnGW)
dlnH = sp.simplify(2*dlnR - dlnOU2 - dlnOW2)
SigmaN_transport = sp.simplify(2*dlnM + 2*I/(1+I)*dlnI + 2*H/(1-H)*dlnH)
print("delta ln M =", dlnM)
print("delta ln I =", dlnI)
print("delta ln H =", dlnH)
print("Sigma_N =")
sp.pprint(SigmaN_transport)

subbanner("4. Fully expanded primitive-port formula")
SigmaN_expanded = sp.expand(SigmaN_transport)
print("Expanded Sigma_N =")
sp.pprint(SigmaN_expanded)

subbanner("5. Rigidity theorem and square-root mixed-leg law")
Sigma_rigid = sp.simplify(SigmaN_transport.subs({dlnI: 0, dlnH: 0}))
expect_zero("rigid branch => Sigma_N - 2 dln M", Sigma_rigid - 2*dlnM)
expect_zero("square-root mixed-leg law residual", dlnM.subs({dlnGW: dlnOW2 + deltaK/2}))
print("If dln I = dln H = 0, then Sigma_N =", Sigma_rigid)
print("So zero defect requires GW/OW2 to load as sqrt(K).")

subbanner("6. Dominant-port corollary")
print("If one port dominates, Xi_load ≈ Sigma_(dominant port).")
print("Under interference/hybridization rigidity, Xi_load ≈ 2 dln(GW/(OW2 sqrt(K))).")

banner("STAGE 159 LEDGER")
print("1. The outgoing load factor splits exactly into")
print("      raw mixed-leg load × interference ratio × hybridization ratio.")
print("2. The defect field is")
print("      Sigma_N = 2 dln M + 2 I/(1+I) dln I + 2 H/(1-H) dln H.")
print("3. So the remaining linear grouped defect is the weighted failure of three")
print("   microscopic loading factors to track one another.")
print("4. If interference and hybridization are rigid, the whole defect collapses to")
print("   the square-root wall-loading law for the raw mixed leg:")
print("      GW/OW2 ∝ sqrt(K).")
