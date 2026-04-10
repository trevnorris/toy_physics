#!/usr/bin/env python3
"""
5pn_stage18_full_bundle_isotropic_targets.py

Eighteenth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Writes the exact isotropic full-bundle conservative coefficients
      D0 = K - B0 - Z0,
      D2 = -(M + B2 + Z2),
      D4 = -(B4 + Z4).
2. Derives the exact grouped-response moments
      u2, u4,
   and the exact outgoing-prefactor moments
      P0, P2, P4.
3. Computes the exact single-pole defect
      u4 - 4 u2^2,
   and solves it for the full-bundle 5PN support datum B4 + Z4.
4. Solves the constant-prefactor conditions
      P2 = 0,
      P4 = 0
   for the outgoing-transfer moments N2 and N4.
5. Solves the universal 2.5PN/4PN normalization condition for N0 and records
   the complete isotropic target surface.

Interpretation
--------------
This script does not touch weak-axisymmetric defects. It sharpens the isotropic
full-bundle theorem target: once the completed moving-throat PDE supplies the
coupled moments (B_n, Z_n, N_n), the 5PN conservative branch and the 2.5PN/4PN
outgoing normalization are constrained by one exact algebraic target surface.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("I. EXACT ISOTROPIC FULL-BUNDLE COEFFICIENTS")

K, M = sp.symbols("K M", positive=True, real=True)
B0, B2, B4 = sp.symbols("B0 B2 B4", real=True)
Z0, Z2, Z4 = sp.symbols("Z0 Z2 Z4", real=True)
N0, N2, N4 = sp.symbols("N0 N2 N4", real=True)

D0 = sp.simplify(K - B0 - Z0)
D2 = sp.simplify(-(M + B2 + Z2))
D4 = sp.simplify(-(B4 + Z4))

print("D0 =", D0)
print("D2 =", D2)
print("D4 =", D4)

u2 = sp.simplify(-D2 / D0)
u4 = sp.simplify((D2**2 - D0 * D4) / D0**2)
P0 = sp.simplify(N0 / D0)
P2 = sp.simplify((D0 * N2 - 2 * D2 * N0) / D0**2)
P4 = sp.simplify((D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3)

print("\nu2 =")
sp.pprint(u2)
print("\nu4 =")
sp.pprint(u4)
print("\nP0 =")
sp.pprint(P0)
print("\nP2 =")
sp.pprint(P2)
print("\nP4 =")
sp.pprint(P4)


banner("II. EXACT SINGLE-POLE DEFECT")

single_pole_defect = sp.simplify(u4 - 4 * u2**2)
print("u4 - 4 u2^2 =")
sp.pprint(single_pole_defect)

single_pole_expected = sp.simplify((D0 * (B4 + Z4) - 3 * (M + B2 + Z2) ** 2) / D0**2)
expect_zero("single-pole defect bridge", single_pole_defect - single_pole_expected)

B4Z4_target = sp.simplify(sp.solve(sp.Eq(single_pole_defect, 0), B4 + Z4)[0])
print("\nExact isotropic 5PN one-pole target:")
print("  B4 + Z4 =")
sp.pprint(B4Z4_target)

subbanner("Support interpretation")
print("The isotropic 5PN conservative branch is one-pole iff")
print("  D0 (B4 + Z4) = 3 (M + B2 + Z2)^2.")
print("So the full bundle must place the fourth even support moment on that exact surface.")


banner("III. EXACT CONSTANT-PREFACTOR CONDITIONS")

N2_const = sp.simplify(sp.solve(sp.Eq(P2, 0), N2)[0])
N4_const = sp.simplify(sp.solve(sp.Eq(P4.subs(N2, N2_const), 0), N4)[0])

print("P2 = 0 fixes N2 =")
sp.pprint(N2_const)
print("\nP4 = 0 then fixes N4 =")
sp.pprint(N4_const)

expect_zero("P2 after imposing N2_const", sp.simplify(P2.subs(N2, N2_const)))
expect_zero("P4 after imposing N2_const,N4_const", sp.simplify(P4.subs({N2: N2_const, N4: N4_const})))

subbanner("Equivalent D-language form")
expect_zero("N2_const - 2 D2 N0 / D0", sp.simplify(N2_const - 2 * D2 * N0 / D0))
expect_zero(
    "N4_const - (2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0)/D0^2",
    sp.simplify(N4_const - (2 * D0 * (D2 * N2_const + D4 * N0) - 3 * D2**2 * N0) / D0**2),
)

print("So the constant-prefactor branch does not require N2 and N4 to vanish.")
print("It requires them to be locked to the conservative moments D0, D2, D4 exactly.")


banner("IV. UNIVERSAL NORMALIZATION CONDITION")

mhat0, Ggrav, cs, a_th, c_light = sp.symbols("mhat_0 G c_s a c", positive=True, real=True)
P0_target = sp.simplify(54 * Ggrav * cs**5 / (5 * a_th**5 * c_light**5))

N0_target = sp.simplify(sp.solve(sp.Eq(mhat0**2 * P0, P0_target), N0)[0])

print("mhat0^2 P0 target =")
sp.pprint(P0_target)
print("\nEquivalent N0 target =")
sp.pprint(N0_target)

expect_zero(
    "normalization solve check",
    sp.simplify(mhat0**2 * P0.subs(N0, N0_target) - P0_target),
)


banner("V. COMPLETE ISOTROPIC TARGET SURFACE")

compatibility_surface = sp.simplify(N0_target / P0_target)
print("The static outgoing target can be written as")
print("  N0 = P0_target * D0 / mhat0^2, with D0 = K - B0 - Z0.")
print("\nIf the isotropic branch is also one-pole, then")
print("  B4 + Z4 =")
sp.pprint(B4Z4_target)
print("\nIf the isotropic branch is also constant-prefactor, then")
print("  N2 =")
sp.pprint(N2_const)
print("  N4 =")
sp.pprint(N4_const)

print("\nSo the full isotropic moving-throat bundle must land on the combined target:")
print("  D0 = K - B0 - Z0,")
print("  D0 (B4 + Z4) = 3 (M + B2 + Z2)^2,")
print("  mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5),")
print("  N2 = 2 D2 N0 / D0,")
print("  N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0]/D0^2.")

print("\nThat is the sharpest isotropic 5PN / 2.5PN / 4PN full-bundle target surface")
print("I can write before solving the actual moving-throat overlap problem.")
