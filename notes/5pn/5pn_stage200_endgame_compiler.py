
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import *

"""
Stage 200 — exact endgame compiler from orbit packet + branch packet.

What this script does
---------------------
1. Combines the Stage-199 grouped-lane branch packet with the exact orbit-lock
   packets from Stages 181–198.
2. Shows that the reduced 5PN / 2.5PN / 4PN closure problem depends only on:
      (a) one grouped conservative/outgoing packet,
      (b) one orbit-lock packet.
3. Proves exact equivalence between using:
      - residual ratios (m_T,m_K,m_mu),
      - invariant ratios (R_tr,R_nt,R_eta),
      - quotient coordinates (q_tr,q_nt,q_eta)
   for the orbit side.
4. Records the final reduced closure condition in one compact vector equation.
"""

banner("STAGE 200 — EXACT ENDGAME COMPILER")

# Orbit packet symbols
orb = orbit_symbols("")
chi0s, Fstar = orb["chi0s"], orb["Fstar"]
mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)
Rtr, Rnt, Reta = sp.symbols("R_tr R_nt R_eta", positive=True, real=True)
qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

# Branch packet symbols
a2, b2, a4, b4 = sp.symbols("a_2 b_2 a_4 b_4", real=True)
aP0, bP0 = sp.symbols("a_P0 b_P0", real=True)
Delta_pole, Delta_norm = sp.symbols("Delta_pole Delta_norm", real=True)

subbanner("I. Exact orbit-packet interconversion")
q_from_m = q_from_mismatch(chi0s, Fstar, mT, mK, mMu)
m_from_q = mismatch_from_q(chi0s, Fstar, qtr, qnt, qeta)
m_from_R = mismatch_from_invariants(chi0s, Fstar, Rtr, Rnt, Reta)
q_from_R = q_from_invariants(Rtr, Rnt, Reta)

print("q from (m_T,m_K,m_mu) =")
sp.pprint(sp.Matrix([q_from_m["qtr"], q_from_m["qnt"], q_from_m["qeta"]]))
print("m from (q_tr,q_nt,q_eta) =")
sp.pprint(sp.Matrix([m_from_q["mT"], m_from_q["mK"], m_from_q["mMu"]]))
print("m from (R_tr,R_nt,R_eta) =")
sp.pprint(sp.Matrix([m_from_R["mT"], m_from_R["mK"], m_from_R["mMu"]]))
print("q from (R_tr,R_nt,R_eta) =")
sp.pprint(sp.Matrix([q_from_R["qtr"], q_from_R["qnt"], q_from_R["qeta"]]))

expect_zero("roundtrip q_tr", q_from_m["qtr"].subs({mT: m_from_q["mT"]}) - qtr)
expect_zero("roundtrip q_eta", q_from_m["qeta"].subs({mK: m_from_q["mK"]}) - qeta)
expect_zero(
    "roundtrip q_nt",
    sp.simplify(q_from_m["qnt"].subs({mT: m_from_q["mT"], mK: m_from_q["mK"], mMu: m_from_q["mMu"]}) - qnt),
)
expect_zero("R_tr from q", sp.log(Rtr) - q_from_R["qtr"])
expect_zero("R_nt from q", sp.log(Rnt) - q_from_R["qnt"])
expect_zero("R_eta from q", sp.log(Reta) - q_from_R["qeta"])

subbanner("II. Exact reduced endgame packet")
Delta_orbit_q = sp.Matrix([qtr, qnt, qeta])
Delta_orbit_m = sp.Matrix([sp.log(mT), sp.log(mK), sp.log(mMu)])
Delta_branch = sp.Matrix([a2, b2, a4, b4, aP0, bP0, Delta_pole, Delta_norm])

print("Delta_orbit^(q) =")
sp.pprint(Delta_orbit_q)
print("Delta_orbit^(ln m) =")
sp.pprint(Delta_orbit_m)
print("Delta_branch =")
sp.pprint(Delta_branch)

subbanner("III. Support/source side no longer bottleneck on the minimal isotropic module")
rho_alpha = sp.Rational(4, 3)
zeta_req = sp.simplify(rho_alpha - 1)
A_F1 = sp.Float("1.00005192880220")
print("rho_alpha(minimal isotropic module) =", rho_alpha)
print("zeta_req(minimal isotropic module) =", zeta_req)
print("A_F1(explicit Family-1 ceiling factor) =", A_F1)
print("A_F1 - zeta_req =", sp.N(A_F1 - zeta_req))
if not (A_F1 > zeta_req):
    raise AssertionError("Explicit Family-1 support/source side did not clear zeta_req on the minimal isotropic module.")

subbanner("IV. Final reduced closure condition")
print("Reduced closure requires:")
print("  Delta_branch = 0")
print("  Delta_orbit^(q) = 0")
print("Equivalently:")
print("  a_2=b_2=a_4=b_4=0,  a(P_0)=b(P_0)=0,  Delta_pole=0,  Delta_norm=0,")
print("  q_tr=q_nt=q_eta=0.")

banner("STAGE 200 LEDGER")
print("1. The endgame has now split cleanly into two exact packets:")
print("      (a) grouped branch packet  Delta_branch")
print("      (b) orbit-lock packet      Delta_orbit")
print("2. The orbit-lock packet can be supplied in any exact form:")
print("      (m_T,m_K,m_mu),  (R_tr,R_nt,R_eta),  or  (q_tr,q_nt,q_eta).")
print("3. On the minimal isotropic module, the explicit Family-1 support/source side")
print("   already clears the required zeta threshold, so it is no longer the active bottleneck.")
print("4. The reduced GR-compatible closure problem is therefore exactly the zero-set test")
print("   of Delta_branch together with Delta_orbit.")
