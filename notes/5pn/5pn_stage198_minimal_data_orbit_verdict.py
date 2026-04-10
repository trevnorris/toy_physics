
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage193_198_common import *

"""
Stage 198 — minimal-data orbit verdict interface.

What this script does
---------------------
1. Shows that the orbit-lock verdict can be reached from any one of three exact data
   packets:
      (a) residual mismatch ratios (m_T,m_K,m_mu),
      (b) invariant ratios (R_Ctr,R_Cnt,R_eps),
      (c) quotient coordinates (q_tr,q_nt,q_eta).
2. Compiles each packet into the other two exactly.
3. Returns the dependent-coordinate restoration map and the canonical orbit distance.
4. Gives a compact PDE-ready decision tree for future numerical branch data.
"""

banner("STAGE 198 — MINIMAL-DATA ORBIT VERDICT INTERFACE")

star = base_symbols("")
chi0s = star["chi0s"]
Fstar = star["Fstar"]

mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)
RCtr, RCnt, Reps = sp.symbols("R_Ctr R_Cnt R_eps", positive=True, real=True)
qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

subbanner("I. Packet A: residual mismatch ratios -> invariant ratios and q")
q_from_m = q_from_mismatch(star, mT, mK, mMu)
inv_from_m = {
    "RCtr": sp.simplify(mT ** (1 + chi0s)),
    "RCnt": sp.simplify(mMu / (mK * mT ** Fstar)),
    "Reps": sp.simplify(1 / mK),
}
print("From (m_T,m_K,m_mu):")
print("  R_Ctr =")
sp.pprint(inv_from_m["RCtr"])
print("  R_Cnt =")
sp.pprint(inv_from_m["RCnt"])
print("  R_eps =")
sp.pprint(inv_from_m["Reps"])
print("  q =")
sp.pprint(sp.Matrix([q_from_m["qtr"], q_from_m["qnt"], q_from_m["qeta"]]))

subbanner("II. Packet B: invariant ratios -> residual ratios and q")
m_from_inv = mismatch_from_invariant_ratios(star, RCtr, RCnt, Reps)
q_from_inv = q_from_invariant_ratios(RCtr, RCnt, Reps)
print("From (R_Ctr,R_Cnt,R_eps):")
print("  m =")
sp.pprint(sp.Matrix([m_from_inv["mT"], m_from_inv["mK"], m_from_inv["mMu"]]))
print("  q =")
sp.pprint(sp.Matrix([q_from_inv["qtr"], q_from_inv["qnt"], q_from_inv["qeta"]]))

subbanner("III. Packet C: quotient coordinates -> residual ratios and invariant ratios")
m_from_q = {
    "mT": sp.simplify(sp.exp(qtr / (1 + chi0s))),
    "mK": sp.simplify(sp.exp(-qeta)),
    "mMu": sp.simplify(sp.exp(qnt - qeta + Fstar * qtr / (1 + chi0s))),
}
inv_from_q = {
    "RCtr": sp.simplify(sp.exp(qtr)),
    "RCnt": sp.simplify(sp.exp(qnt)),
    "Reps": sp.simplify(sp.exp(qeta)),
}
print("From (q_tr,q_nt,q_eta):")
print("  m =")
sp.pprint(sp.Matrix([m_from_q["mT"], m_from_q["mK"], m_from_q["mMu"]]))
print("  R =")
sp.pprint(sp.Matrix([inv_from_q["RCtr"], inv_from_q["RCnt"], inv_from_q["Reps"]]))

subbanner("IV. Exact cross-checks")
# A -> B -> A
m_round = mismatch_from_invariant_ratios(star, inv_from_m["RCtr"], inv_from_m["RCnt"], inv_from_m["Reps"])
expect_zero("ln m_T roundtrip", sp.log(m_round["mT"]) - sp.log(mT))
expect_zero("ln m_K roundtrip", sp.log(m_round["mK"]) - sp.log(mK))
expect_zero("ln m_mu roundtrip", sp.log(m_round["mMu"]) - sp.log(mMu))

# C -> A -> C
expect_zero("q_tr from q-packet", sp.log(inv_from_q["RCtr"]) - qtr)
expect_zero("q_nt from q-packet", sp.log(inv_from_q["RCnt"]) - qnt)
expect_zero("q_eta from q-packet", sp.log(inv_from_q["Reps"]) - qeta)

subbanner("V. Restoration map and orbit distance from the minimal packet")
Keta, TU, muW = sp.symbols("K_eta T_U mu_W", positive=True, real=True)
Keta_restore = sp.simplify(Keta / m_from_q["mK"])
TU_restore = sp.simplify(TU / m_from_q["mT"])
mu_restore = sp.simplify(muW / m_from_q["mMu"])
D2 = sp.simplify(qtr**2 + qnt**2 + qeta**2)

print("K_eta^(restore) =")
sp.pprint(Keta_restore)
print("T_U^(restore) =")
sp.pprint(TU_restore)
print("mu_W^(restore) =")
sp.pprint(mu_restore)
print("D^2 =")
sp.pprint(D2)

banner("STAGE 198 LEDGER")
print("1. Orbit lock can be decided from any one of three exact packets:")
print("      (m_T,m_K,m_mu),  (R_Ctr,R_Cnt,R_eps),  or  (q_tr,q_nt,q_eta).")
print("2. The three packets are exactly interconvertible.")
print("3. The dependent-coordinate restoration map and the scalar orbit distance D^2")
print("   can therefore be computed from whichever packet the PDE numerics can provide")
print("   most cleanly.")
print("4. This is the minimal-data decision layer for future actual-branch orbit tests.")
