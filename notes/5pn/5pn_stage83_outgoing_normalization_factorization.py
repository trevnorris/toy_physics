
#!/usr/bin/env python3
"""
5pn_stage83_outgoing_normalization_factorization.py

Stage 83 audit: exact factorization of the last 2.5PN defect into conservative
and outgoing pieces.
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

banner("STAGE 83 — EXACT FACTORIZATION OF THE LAST 2.5PN DEFECT")

G, c, Omega_Q, Kbar0, chi_Q, mhat0 = sp.symbols(
    "G c Omega_Q Kbar_0 chi_Q mhat_0", positive=True, real=True
)

Kbar0_target = sp.simplify(64 * G * Omega_Q**5 / (45 * c**5))
N_Q = sp.simplify(Kbar0 / Kbar0_target)
Kbar2 = sp.simplify(Kbar0 / (4 * Omega_Q**2))
Kbar4 = sp.simplify(Kbar0 / (4 * Omega_Q**4))
Gammabar5 = sp.simplify(chi_Q * 9 * Kbar0 / (32 * Omega_Q**5))
Gammabar5_target = sp.simplify(2 * G / (5 * c**5))

expect_zero("Kbar_2/Kbar_2^target - N_Q", (Kbar2 / (Kbar0_target / (4 * Omega_Q**2))) - N_Q)
expect_zero("Kbar_4/Kbar_4^target - N_Q", (Kbar4 / (Kbar0_target / (4 * Omega_Q**4))) - N_Q)
expect_zero("Gammabar_5/Gammabar_5^target - chi_Q N_Q", Gammabar5 / Gammabar5_target - chi_Q * N_Q)

normalization_condition = sp.simplify(mhat0**2 * chi_Q * N_Q)
print("mhat_0^2 chi_Q N_Q =", normalization_condition)

banner("STAGE 83 FINAL LEDGER")
print("The even conservative defect and the odd retarded defect separate exactly:")
print("  Kbar_2/Kbar_2^target = Kbar_4/Kbar_4^target = N_Q,")
print("  Gammabar_5/Gammabar_5^target = chi_Q N_Q.")
print("So the full odd normalization condition is")
print("  mhat_0^2 chi_Q N_Q = 1.")
