
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage193_198_common import *

"""
Stage 195 — pairwise quotient-to-observable compiler.

What this script does
---------------------
1. Composes the exact pairwise quotient coordinates with the Stage-170 linear observable map.
2. Produces the linear weak-axisymmetric observable triple predicted by any small pairwise
   orbit mismatch.
3. Verifies the pure-channel signatures for isolated T_U, K_eta, and mu_W residuals.
4. Gives the exact inverse linear map from the observable triple back to the quotient
   coordinates, and then to the residual mismatch ratios.
"""

banner("STAGE 195 — PAIRWISE QUOTIENT-TO-OBSERVABLE COMPILER")

star = base_symbols("")
chi0s = star["chi0s"]
deltaUs = star["deltaUs"]
Fstar = star["Fstar"]

t, k, m = sp.symbols("z_T z_K z_M", positive=True, real=True)
mT = t
mK = k
mMu = m

subbanner("I. Exact quotient coordinates from pairwise residual ratios")
q = q_from_mismatch(star, mT, mK, mMu)
print("q_tr =")
sp.pprint(q["qtr"])
print("q_nt =")
sp.pprint(q["qnt"])
print("q_eta =")
sp.pprint(q["qeta"])

subbanner("II. Linear observable compiler")
obs = linear_observable_map(star, q["qtr"], q["qnt"], q["qeta"])
print("Theta_1^(lin) =")
sp.pprint(obs["Theta"])
print("Xi_1^(lin) =")
sp.pprint(obs["Xi"])
print("R_1 + Xi_1^(lin) =")
sp.pprint(obs["Rsum"])
print("R_1^(lin) =")
sp.pprint(obs["R1"])

subbanner("III. Pure-channel signatures")

# Pure T_U mismatch
obs_T = linear_observable_map(star, *q_from_mismatch(star, t, 1, 1).values())
print("Pure T_U mismatch:")
print("  q_tr =")
sp.pprint(q_from_mismatch(star, t, 1, 1)["qtr"])
print("  q_nt =")
sp.pprint(q_from_mismatch(star, t, 1, 1)["qnt"])
print("  q_eta =")
sp.pprint(q_from_mismatch(star, t, 1, 1)["qeta"])
print("  Theta_1^(lin) =")
sp.pprint(obs_T["Theta"])
print("  Xi_1^(lin) =")
sp.pprint(obs_T["Xi"])
print("  R_1 + Xi_1^(lin) =")
sp.pprint(obs_T["Rsum"])

# Pure K_eta mismatch
obs_K = linear_observable_map(star, *q_from_mismatch(star, 1, k, 1).values())
print("\nPure K_eta mismatch:")
print("  q_tr =")
sp.pprint(q_from_mismatch(star, 1, k, 1)["qtr"])
print("  q_nt =")
sp.pprint(q_from_mismatch(star, 1, k, 1)["qnt"])
print("  q_eta =")
sp.pprint(q_from_mismatch(star, 1, k, 1)["qeta"])
print("  Theta_1^(lin) =")
sp.pprint(obs_K["Theta"])
print("  Xi_1^(lin) =")
sp.pprint(obs_K["Xi"])
print("  R_1 + Xi_1^(lin) =")
sp.pprint(obs_K["Rsum"])

# Pure mu_W mismatch
obs_M = linear_observable_map(star, *q_from_mismatch(star, 1, 1, m).values())
print("\nPure mu_W mismatch:")
print("  q_tr =")
sp.pprint(q_from_mismatch(star, 1, 1, m)["qtr"])
print("  q_nt =")
sp.pprint(q_from_mismatch(star, 1, 1, m)["qnt"])
print("  q_eta =")
sp.pprint(q_from_mismatch(star, 1, 1, m)["qeta"])
print("  Theta_1^(lin) =")
sp.pprint(obs_M["Theta"])
print("  Xi_1^(lin) =")
sp.pprint(obs_M["Xi"])
print("  R_1 + Xi_1^(lin) =")
sp.pprint(obs_M["Rsum"])

subbanner("IV. Exact inverse linear map from observables back to q and m")
Theta1, Xi1, Rsum = sp.symbols("Theta_1 Xi_1 Rsum", real=True)
qtr_inv = sp.simplify(-(1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs) / (chi0s * deltaUs) * Theta1)
qnt_inv = sp.simplify(Xi1 + 2 * (1 + chi0s + deltaUs) / deltaUs * Theta1)
qeta_inv = sp.simplify(-(1 - obs["eps_eta_star"]) / obs["eps_eta_star"] * Rsum)

print("q_tr^(inv) =")
sp.pprint(qtr_inv)
print("q_nt^(inv) =")
sp.pprint(qnt_inv)
print("q_eta^(inv) =")
sp.pprint(qeta_inv)

mT_inv = sp.simplify(sp.exp(qtr_inv / (1 + chi0s)))
mK_inv = sp.simplify(sp.exp(-qeta_inv))
mMu_inv = sp.simplify(sp.exp(qnt_inv - qeta_inv + Fstar * qtr_inv / (1 + chi0s)))

print("m_T^(inv) =")
sp.pprint(mT_inv)
print("m_K^(inv) =")
sp.pprint(mK_inv)
print("m_mu^(inv) =")
sp.pprint(mMu_inv)

banner("STAGE 195 LEDGER")
print("1. Any small pairwise orbit mismatch has a universal linear weak-axisymmetric observable")
print("   signature once written in the exact quotient coordinates.")
print("2. Pure T_U mismatch drives a coupled (Theta_1, Xi_1) pattern with q_nt = -F_* ln m_T.")
print("3. Pure K_eta mismatch drives only q_eta and q_nt, with Theta_1 = 0.")
print("4. Pure mu_W mismatch drives only q_nt, hence Xi_1 alone at linear order.")
print("5. The linear observable map is exactly invertible back to q, and hence back to the")
print("   finite residual ratios by exponentiation.")
