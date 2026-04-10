
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import *

"""
Stage 201 — home-stretch theorem and minimal PDE data packet.

What this script does
---------------------
1. Packages the final reduced closure problem as a theorem with explicit hypotheses.
2. Identifies the smallest exact PDE data packet that is still needed.
3. Separates what is already settled inside the hierarchy from what the completed
   moving-throat PDE must still compute.
4. Produces a compact decision tree for the final numerical / analytic endgame.
"""

banner("STAGE 201 — HOME-STRETCH THEOREM AND MINIMAL PDE DATA PACKET")

G, c, cs, a, mhat0 = sp.symbols("G c c_s a mhat_0", positive=True, real=True)
chi0s, Fstar = sp.symbols("chi0_star F_star", positive=True, real=True)

subbanner("I. Minimal exact PDE data packet")
print("Packet A — grouped bundle data:")
print("  (D_A0, D_A2, D_A4, N_A0, N_A2, N_A4) for A in {20,21,22}")
print("  plus the source-map factor mhat_0")
print()
print("Packet B — orbit/invariant data:")
print("  one of")
print("    (m_T, m_K, m_mu),")
print("    (R_tr, R_nt, R_eta),")
print("    (q_tr, q_nt, q_eta)")
print()
print("Everything else in the reduced closure test is an exact compiler output of these packets.")

subbanner("II. Final reduced closure equations")
a2, b2, a4, b4, aP0, bP0, Delta_pole, Delta_norm = sp.symbols(
    "a_2 b_2 a_4 b_4 a_P0 b_P0 Delta_pole Delta_norm", real=True
)
qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

Delta_branch = sp.Matrix([a2, b2, a4, b4, aP0, bP0, Delta_pole, Delta_norm])
Delta_orbit = sp.Matrix([qtr, qnt, qeta])

print("Delta_branch =")
sp.pprint(Delta_branch)
print("Delta_orbit =")
sp.pprint(Delta_orbit)

print("Reduced closure theorem:")
print("  Delta_branch = 0  and  Delta_orbit = 0")

subbanner("III. What is already settled inside the hierarchy")
print("Settled inside the current hierarchy:")
print("  - the support/source side of the minimal isotropic module is already non-bottlenecked;")
print("  - the adiabatic-elastic branch theorem reduces weak-axisymmetric closure to orbit lock;")
print("  - the orbit-lock verdict is exactly the zero-set of (q_tr,q_nt,q_eta);")
print("  - the grouped real P2 conservative/outgoing test is exactly the zero-set of Delta_branch.")

subbanner("IV. What the completed moving-throat PDE still has to compute")
D0, N0 = sp.symbols("D_0 N_0", positive=True, real=True)
P0 = sp.simplify(N0 / D0)
Delta_norm_iso = sp.simplify(mhat0**2 * P0 - 54*G*cs**5/(5*a**5*c**5))
print("Still open on the actual branch:")
print("  1. the true grouped-lane bundle data and their isotropy defects;")
print("  2. the one-pole conservative defect Delta_pole;")
print("  3. the quadrupole-normalization hit")
print("         Delta_norm =")
sp.pprint(Delta_norm_iso)
print("  4. the dependent residual triple, equivalently the orbit packet.")

subbanner("V. Compact decision tree")
print("Step 1: compile Packet A -> Delta_branch.")
print("Step 2: compile Packet B -> Delta_orbit.")
print("Step 3: if any component of Delta_branch is nonzero, the branch fails the reduced GR test.")
print("Step 4: if Delta_branch = 0 but Delta_orbit != 0, the branch is isotropic-but-off-orbit.")
print("Step 5: only if both packets vanish is the reduced closure complete inside the current hierarchy.")

banner("STAGE 201 LEDGER")
print("1. The endgame is now fully localized to two exact finite packets.")
print("2. The smallest unresolved theorem problem is no longer diffuse.")
print("3. The completed moving-throat PDE only has to supply the data needed to evaluate")
print("      Delta_branch and Delta_orbit.")
print("4. Once those are known, the reduced verdict is immediate and exact.")
