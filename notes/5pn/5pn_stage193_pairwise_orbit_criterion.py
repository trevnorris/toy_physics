
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage193_198_common import *

"""
Stage 193 — exact pairwise orbit criterion.

What this script does
---------------------
1. Takes two arbitrary positive microscopic states x and y with the same branch constants
   (chi0_*, deltaU_*, E_*, F_*).
2. Builds the free-coordinate ratios, the exact orbit-predicted dependent ratios, and the
   actual dependent ratios.
3. Defines the exact pairwise residual mismatch ratios
      (m_T^(x->y), m_K^(x->y), m_mu^(x->y)).
4. Proves that the invariant ratios between y and x depend only on those three residual
   ratios, with the five free-coordinate ratios cancelling exactly.
5. Concludes that x and y lie on the same exact similarity orbit iff the residual triple
   is trivial, equivalently iff the invariant ratios are all unity.
"""

banner("STAGE 193 — EXACT PAIRWISE ORBIT CRITERION")

# Use a shared branch-constant dictionary and two independent microscopic states.
star = base_symbols("")
x = base_symbols("x_")
y = base_symbols("y_")

# Identify shared branch constants.
for state in (x, y):
    state["chi0s"] = star["chi0s"]
    state["deltaUs"] = star["deltaUs"]
    state["Estar"] = star["Estar"]
    state["Fstar"] = star["Fstar"]
    state["sigma"] = star["sigma"]
    state["L"] = star["L"]

subbanner("I. Pairwise free-coordinate ratios and orbit prediction")
pair = pairwise_mismatch_ratios(x, y, star)

print("Free-coordinate ratios R_(x->y) =")
for k in ("Rlam", "Rc", "Rgam", "RU", "RW"):
    print(f"  {k} =")
    sp.pprint(pair["rats"][k])

print("\nOrbit-predicted dependent ratios =")
print("  R_Keta^(orbit) =")
sp.pprint(pair["orbit"]["RKeta"])
print("  R_TU^(orbit) =")
sp.pprint(pair["orbit"]["RTU"])
print("  R_muW^(orbit) =")
sp.pprint(pair["orbit"]["Rmu"])

print("\nActual dependent ratios =")
print("  R_Keta^(act) =")
sp.pprint(pair["rats"]["RKeta_act"])
print("  R_TU^(act) =")
sp.pprint(pair["rats"]["RTU_act"])
print("  R_muW^(act) =")
sp.pprint(pair["rats"]["Rmu_act"])

subbanner("II. Exact pairwise residual mismatch ratios")
mT = pair["mT"]; mK = pair["mK"]; mMu = pair["mMu"]
print("m_T^(x->y) =")
sp.pprint(mT)
print("m_K^(x->y) =")
sp.pprint(mK)
print("m_mu^(x->y) =")
sp.pprint(mMu)

subbanner("III. Invariant ratios depend only on the residual triple")
inv_rat = pairwise_invariant_ratios(x, y)
print("R_Ctr^(x->y) =")
sp.pprint(inv_rat["RCtr"])
print("R_Cnt^(x->y) =")
sp.pprint(inv_rat["RCnt"])
print("R_eps^(x->y) =")
sp.pprint(inv_rat["Reps"])

expect_zero(
    "ln R_Ctr - (1+chi0_*) ln m_T",
    sp.log(inv_rat["RCtr"]) - (1 + star["chi0s"]) * sp.log(mT),
)
expect_zero(
    "ln R_eps + ln m_K",
    sp.log(inv_rat["Reps"]) + sp.log(mK),
)
expect_zero(
    "ln R_Cnt - (ln m_mu - ln m_K - F_* ln m_T)",
    sp.log(inv_rat["RCnt"]) - (sp.log(mMu) - sp.log(mK) - star["Fstar"] * sp.log(mT)),
)

subbanner("IV. Exact same-orbit criterion")
print("x and y lie on the same exact similarity orbit iff any one of the following equivalent")
print("pairwise tests holds:")
print("  (a) m_T^(x->y) = m_K^(x->y) = m_mu^(x->y) = 1")
print("  (b) R_Ctr^(x->y) = R_Cnt^(x->y) = R_eps^(x->y) = 1")
print("  (c) the actual dependent ratios equal the orbit-predicted dependent ratios")
print("      at the same five free-coordinate ratios.")

banner("STAGE 193 LEDGER")
print("1. The finite orbit-lock test is now reference-independent: any two candidate states")
print("   can be compared directly.")
print("2. The five free-coordinate ratios drop out of the invariant-ratio test exactly.")
print("3. The entire pairwise orbit question collapses to the residual triple")
print("      (m_T^(x->y), m_K^(x->y), m_mu^(x->y)).")
