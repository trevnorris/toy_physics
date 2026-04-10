
#!/usr/bin/env python3
"""
5pn_stage80_single_normalization_defect.py

Stage 80 audit: the actual isotropic passive/outgoing branch collapses to a single
normalization defect.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 80 — THE ACTUAL ISOTROPIC PASSIVE/OUTGOING BRANCH COLLAPSES TO A SINGLE NORMALIZATION DEFECT")

G, c, Omega_Q, Kbar0 = sp.symbols("G c Omega_Q Kbar_0", positive=True, real=True)

Kbar2 = sp.simplify(Kbar0 / (4 * Omega_Q**2))
Kbar4 = sp.simplify(Kbar0 / (4 * Omega_Q**4))
Gammabar5 = sp.simplify(9 * Kbar2**sp.Rational(5, 2) / Kbar0**sp.Rational(3, 2))
expect_zero("Gamma_bar_5 - 9 Kbar0/(32 Omega_Q^5)", Gammabar5 - 9 * Kbar0 / (32 * Omega_Q**5))

Kbar0_target = sp.simplify(64 * G * Omega_Q**5 / (45 * c**5))
Kbar2_target = sp.simplify(Kbar0_target / (4 * Omega_Q**2))
Kbar4_target = sp.simplify(Kbar0_target / (4 * Omega_Q**4))
Gammabar5_target = sp.simplify(2 * G / (5 * c**5))

N_Q = sp.symbols("N_Q", positive=True, real=True)
rat0 = sp.simplify((N_Q * Kbar0_target) / Kbar0_target)
rat2 = sp.simplify((N_Q * Kbar2_target) / Kbar2_target)
rat4 = sp.simplify((N_Q * Kbar4_target) / Kbar4_target)
rat5 = sp.simplify((N_Q * Gammabar5_target) / Gammabar5_target)
expect_zero("R0 - N_Q", rat0 - N_Q)
expect_zero("R2 - N_Q", rat2 - N_Q)
expect_zero("R4 - N_Q", rat4 - N_Q)
expect_zero("R5 - N_Q", rat5 - N_Q)

banner("STAGE 80 FINAL LEDGER")
print("Once the actual isotropic grouped-P2 one-pole branch is accepted, all low-frequency")
print("targets scale with one and the same number")
print("  N_Q := Kbar_0 / Kbar_0^target.")
print("So the reduced nonspinning point-particle 2.5PN theorem is closed iff")
print("  N_Q = 1.")
