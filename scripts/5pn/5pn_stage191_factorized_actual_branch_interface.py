#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage187_192_common import *

"""
Stage 191 — factorized actual-branch interface relative to a reference orbit.

What this script does
---------------------
1. Packages the exact reference-orbit transport laws and the dependent residual ratios
   into a single actual-branch factorization interface.
2. Shows that restoration to the orbit at fixed free-coordinate ratios is achieved by
   dividing only by the residual mismatch ratios.
3. Produces the direct finite criterion that a candidate branch is orbit-locked iff the
   three residual mismatch ratios are unity.
"""

banner("STAGE 191 — FACTORIZED ACTUAL-BRANCH INTERFACE")

syms = base_symbols()
rat = free_ratio_symbols()
mis = mismatch_ratio_symbols()
laws = orbit_ratio_laws(syms, rat)

subbanner("I. Exact factorized actual-branch laws")
Keta_actual = sp.simplify(syms["Keta"] * laws["RKeta"] * mis["mK"])
TU_actual = sp.simplify(syms["TU"] * laws["RTU"] * mis["mT"])
mu_actual = sp.simplify(syms["muW"] * laws["Rmu"] * mis["mMu"])

print("K_eta(actual) =")
sp.pprint(Keta_actual)
print("T_U(actual) =")
sp.pprint(TU_actual)
print("mu_W(actual) =")
sp.pprint(mu_actual)

subbanner("II. Exact restoration at fixed free-coordinate ratios")
Keta_restore = sp.simplify(Keta_actual / mis["mK"])
TU_restore = sp.simplify(TU_actual / mis["mT"])
mu_restore = sp.simplify(mu_actual / mis["mMu"])

expect_zero("restore K_eta to orbit transport", Keta_restore / (syms["Keta"] * laws["RKeta"]) - 1)
expect_zero("restore T_U to orbit transport", TU_restore / (syms["TU"] * laws["RTU"]) - 1)
expect_zero("restore mu_W to orbit transport", mu_restore / (syms["muW"] * laws["Rmu"]) - 1)

print("K_eta^(restore) =")
sp.pprint(Keta_restore)
print("T_U^(restore) =")
sp.pprint(TU_restore)
print("mu_W^(restore) =")
sp.pprint(mu_restore)

subbanner("III. Direct orbit-lock criterion")
print("Given a reference orbit point and the five free-coordinate ratios, the actual branch")
print("is on the same exact orbit iff")
print("  m_T = 1,  m_K = 1,  m_mu = 1.")
print("Equivalently, the actual dependent triple equals the reference-orbit transport law:")
print("  K_eta(actual) = K_eta(ref) * R_Keta^(orbit)")
print("  T_U(actual)   = T_U(ref)   * R_TU^(orbit)")
print("  mu_W(actual)  = mu_W(ref)  * R_muW^(orbit)")

banner("STAGE 191 LEDGER")
print("1. The exact actual-branch interface now factors every candidate into")
print("   free-coordinate orbit transport × dependent residual mismatch.")
print("2. Restoration at fixed free-coordinate ratios acts only on the residual triple.")
print("3. Orbit lock is exactly the statement that the residual triple is trivial.")
