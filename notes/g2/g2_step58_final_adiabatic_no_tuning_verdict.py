#!/usr/bin/env python3
"""
g2_step58_final_adiabatic_no_tuning_verdict.py

Step 58 of the g-2 chain.

What this script does
---------------------
1. Takes the no-tuning consequences of Steps 56–57:
       Delta_Q^(nat) = 0,
       chi_Q^(nat) = 1,
       N_Q^(nat) = 1,
       ell^(nat) = 0.
2. Identifies the corresponding g prediction with the carried local closure
   from atom_work.md.
3. Compares that no-tuning prediction to the electron target used throughout the
   write-up.
4. Compares the no-tuning branch to the explicit electron-point outgoing defect
   required by the earlier outgoing-bridge steps.
5. States the strongest final verdict available inside the current closure.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


banner("STEP 58 — FINAL ADIABATIC NO-TUNING VERDICT")

# Carried benchmark values from atom_work.md and the earlier outgoing-bridge steps.
g_loc = sp.Float("2.00231930435865")
g_e = sp.Float("2.00231930436092")
residual_g = sp.N(g_e - g_loc, 20)
residual_g_over_2 = sp.N(residual_g / 2, 20)

Delta_Q_e = sp.Float("-3.24631584151692e-4")
chi_Q_e = sp.Float("0.999675368415848")
N_Q_e = sp.Float("1.00032473700404")
ell_e = sp.N(sp.log(N_Q_e), 20)

# Natural compensated lower-branch first-order values from Steps 56-57.
Delta_Q_nat = sp.Integer(0)
chi_Q_nat = sp.Integer(1)
N_Q_nat = sp.Integer(1)
ell_nat = sp.Integer(0)

print("Natural branch:")
print("  Delta_Q^(nat) =", Delta_Q_nat)
print("  chi_Q^(nat)   =", chi_Q_nat)
print("  N_Q^(nat)     =", N_Q_nat)
print("  ell^(nat)     =", ell_nat)
print()
print("Electron-point outgoing branch values required by the earlier bridge:")
print("  Delta_Q^(e) =", Delta_Q_e)
print("  chi_Q^(e)   =", chi_Q_e)
print("  N_Q^(e)     =", N_Q_e)
print("  ell^(e)     =", ell_e)
print()
print("No-tuning g prediction =", g_loc)
print("Electron target        =", g_e)
print("Residual g_e - g_loc   =", residual_g)
print("Residual in g/2        =", residual_g_over_2)

# Consistency check of the outgoing bridge numbers.
check_NQ = sp.N(1 / (1 + Delta_Q_e), 18)
print("1/(1 + Delta_Q^(e)) =", check_NQ)
print("Difference from carried N_Q^(e) =", sp.N(check_NQ - N_Q_e, 18))

banner("FINAL LEDGER")
print("Steps 56-57 imply that the natural compensated lower branch has")
print("  Delta_Q = 0, chi_Q = 1, N_Q = 1, ell = 0")
print("at first order. So inside the current no-tuning closure there is no outgoing")
print("normalization correction on top of the carried local anomaly law.")
print()
print("Therefore the natural adiabatic no-tuning prediction is")
print(f"  g_pred^(nat) = {g_loc}")
print("which differs from the electron target by")
print(f"  g_e - g_pred^(nat) = {residual_g}")
print()
print("The exact electron-point outgoing deformation that would be needed instead is")
print(f"  Delta_Q^(e) = {Delta_Q_e}")
print(f"  chi_Q^(e)   = {chi_Q_e}")
print(f"  N_Q^(e)     = {N_Q_e}")
print(f"  ell^(e)     = {ell_e}")
print()
print("So the strongest honest conclusion inside the present closure is:")
print("  - the canonical background / no-tuning branch is naturally derived;")
print("  - the exact electron sliver is not naturally forced by that branch;")
print("  - reproducing the exact electron value still needs either")
print("      (i) a genuine off-family scalar slippage epsilon_perp != 0,")
print("      (ii) a direct odd mixed-port renormalization delta_gamma_W != 0,")
print("      or (iii) a beyond-first-order effect outside the current reduced closure.")
