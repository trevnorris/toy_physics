
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage193_198_common import *

"""
Stage 196 — exact inversion from invariant ratios to dependent residual ratios.

What this script does
---------------------
1. Starts from the reference-independent invariant-ratio packet
      (R_Ctr, R_Cnt, R_eps)
   between two candidate states.
2. Solves exactly for the pairwise residual mismatch ratios
      (m_T, m_K, m_mu).
3. Verifies that the restored dependent triple constructed from those mismatch ratios
   reproduces the same invariant ratios exactly.
4. Shows that free-coordinate ratios are not needed once the invariant packet is known.
"""

banner("STAGE 196 — EXACT INVARIANT-RATIO INVERSION")

star = base_symbols("")
RCtr, RCnt, Reps = sp.symbols("R_Ctr R_Cnt R_eps", positive=True, real=True)

subbanner("I. Exact inversion formulas")
mis = mismatch_from_invariant_ratios(star, RCtr, RCnt, Reps)
print("m_T =")
sp.pprint(mis["mT"])
print("m_K =")
sp.pprint(mis["mK"])
print("m_mu =")
sp.pprint(mis["mMu"])

subbanner("II. Exact consistency check")
inv_rebuilt = {
    "RCtr": sp.simplify(mis["mT"] ** (1 + star["chi0s"])),
    "RCnt": sp.simplify(mis["mMu"] / (mis["mK"] * mis["mT"] ** star["Fstar"])),
    "Reps": sp.simplify(1 / mis["mK"]),
}
expect_zero("ln RCtr rebuilt - ln RCtr", sp.log(inv_rebuilt["RCtr"]) - sp.log(RCtr))
expect_zero("ln RCnt rebuilt - ln RCnt", sp.log(inv_rebuilt["RCnt"]) - sp.log(RCnt))
expect_zero("ln Reps rebuilt - ln Reps", sp.log(inv_rebuilt["Reps"]) - sp.log(Reps))

subbanner("III. Restoration formulas from invariant-ratio data alone")
Keta, TU, muW = sp.symbols("K_eta T_U mu_W", positive=True, real=True)
Keta_restore = sp.simplify(Keta / mis["mK"])
TU_restore = sp.simplify(TU / mis["mT"])
mu_restore = sp.simplify(muW / mis["mMu"])

print("K_eta^(restore) =")
sp.pprint(Keta_restore)
print("T_U^(restore) =")
sp.pprint(TU_restore)
print("mu_W^(restore) =")
sp.pprint(mu_restore)

banner("STAGE 196 LEDGER")
print("1. The invariant-ratio packet (R_Ctr,R_Cnt,R_eps) already contains the full pairwise")
print("   orbit-lock information.")
print("2. The dependent residual ratios are recovered exactly by")
print("      m_T = R_Ctr^(1/(1+chi0_*)),")
print("      m_K = 1/R_eps,")
print("      m_mu = R_Cnt * m_K * m_T^F_*.")
print("3. So once invariant ratios are available, free-coordinate ratios are no longer needed")
print("   for the orbit-lock verdict or the dependent-coordinate restoration.")
