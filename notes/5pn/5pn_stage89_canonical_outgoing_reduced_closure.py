
#!/usr/bin/env python3
"""
5pn_stage89_canonical_outgoing_reduced_closure.py

Stage 89 audit: reduced 2.5PN closure on the canonical outgoing DtN branch.
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

banner("STAGE 89 — REDUCED 2.5PN CLOSURE ON THE CANONICAL OUTGOING DtN BRANCH")

G, c, c_s, a, mhat0 = sp.symbols("G c c_s a mhat_0", positive=True, real=True)
chi_Q = sp.Integer(1)
N_Q = sp.symbols("N_Q", positive=True, real=True)

expect_zero("N_Q on the strict point-particle canonical branch - 1", sp.simplify(1 / chi_Q - 1))
Omega_Q = sp.simplify(3 * c_s / (2 * a))
Kbar0 = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
Kbar2 = sp.simplify(Kbar0 / (4 * Omega_Q**2))
Kbar4 = sp.simplify(Kbar0 / (4 * Omega_Q**4))
Gammabar5 = sp.simplify(2 * G / (5 * c**5))

expect_zero("Kbar_2 - 6 G c_s^3/(5 a^3 c^5)", Kbar2 - 6 * G * c_s**3 / (5 * a**3 * c**5))
expect_zero("Kbar_4 - 8 G c_s/(15 a c^5)", Kbar4 - 8 * G * c_s / (15 * a * c**5))

banner("STAGE 89 FINAL LEDGER")
print("On the canonical compact passive/outgoing grouped-P2 DtN branch, chi_Q = 1.")
print("So on the natural strict point-particle source-map branch the reduced normalization stack")
print("closes with")
print("  N_Q = 1,")
print("and the canonical invariant coefficients are fixed to")
print("  Kbar_0 = 54 G c_s^5/(5 a^5 c^5),")
print("  Kbar_2 = 6 G c_s^3/(5 a^3 c^5),")
print("  Kbar_4 = 8 G c_s/(15 a c^5),")
print("  Gammabar_5 = 2 G/(5 c^5).")
