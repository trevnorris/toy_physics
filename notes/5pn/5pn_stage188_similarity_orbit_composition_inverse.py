#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage187_192_common import *

"""
Stage 188 — exact composition, inverse, and parameter recovery for the finite orbit.

What this script does
---------------------
1. Proves that the finite similarity-orbit action composes by addition of the five
   generator exponents.
2. Proves that the inverse action is obtained by negating those exponents.
3. Shows that, for orbit-related states, the five free exponents are recovered exactly
   from the ratios of the five free microscopic coordinates.
"""

banner("STAGE 188 — EXACT ORBIT COMPOSITION, INVERSE, AND PARAMETER RECOVERY")

syms = base_symbols()
A = {k: sp.symbols(f"{k}_A", real=True) for k in ("Lam", "C", "Gam", "U", "W")}
B = {k: sp.symbols(f"{k}_B", real=True) for k in ("Lam", "C", "Gam", "U", "W")}

symsA = orbit_transform(syms, A)
symsAB = orbit_transform(symsA, B)
ABsum = {k: sp.simplify(A[k] + B[k]) for k in A}
symsSum = orbit_transform(syms, ABsum)

subbanner("I. Group law: composition by exponent addition")
for key in ("lamW", "c_etaU", "gamma", "KU", "Keta", "KW", "muW", "TU"):
    expect_zero(f"compose[{key}] - sum[{key}]", sp.simplify(symsAB[key] / symsSum[key] - 1))

subbanner("II. Exact inverse action")
Ainv = {k: -A[k] for k in A}
symsId = orbit_transform(symsA, Ainv)
for key in ("lamW", "c_etaU", "gamma", "KU", "Keta", "KW", "muW", "TU"):
    expect_zero(f"inverse[{key}] - identity", sp.simplify(symsId[key] / syms[key] - 1))

subbanner("III. Parameter recovery from orbit-related free-coordinate ratios")
lamWp, cetaUp, gammap, KUp, KWp = sp.symbols("lambda_Wp c_etaUp gammap K_Up K_Wp", positive=True, real=True)
Lam_rec = sp.log(lamWp / syms["lamW"])
C_rec = sp.log(cetaUp / syms["c_etaU"])
Gam_rec = sp.log(gammap / syms["gamma"])
U_rec = sp.log(KUp / syms["KU"])
W_rec = sp.log(KWp / syms["KW"])

print("Recovered exponents from free-coordinate ratios:")
print("Lambda = log(lambda_W'/lambda_W)")
print("C      = log(c_etaU'/c_etaU)")
print("Gamma  = log(gamma'/gamma)")
print("U      = log(K_U'/K_U)")
print("W      = log(K_W'/K_W)")

Arec = {"Lam": Lam_rec, "C": C_rec, "Gam": Gam_rec, "U": U_rec, "W": W_rec}
symsRec = orbit_transform(syms, Arec)
expect_zero("recover lambda_W'", symsRec["lamW"] - lamWp)
expect_zero("recover c_etaU'", symsRec["c_etaU"] - cetaUp)
expect_zero("recover gamma'", symsRec["gamma"] - gammap)
expect_zero("recover K_U'", symsRec["KU"] - KUp)
expect_zero("recover K_W'", symsRec["KW"] - KWp)

print("Predicted dependent transport from recovered exponents:")
print("K_eta' =")
sp.pprint(symsRec["Keta"])
print("T_U' =")
sp.pprint(symsRec["TU"])
print("mu_W' =")
sp.pprint(symsRec["muW"])

banner("STAGE 188 LEDGER")
print("1. The finite similarity action forms an exact abelian five-parameter group.")
print("2. Its inverse is obtained by exponent negation.")
print("3. For orbit-related states, the five free generator parameters are recovered")
print("   exactly from the five free-coordinate ratios.")
print("4. Once those ratios are known, the dependent triple is predicted uniquely.")
