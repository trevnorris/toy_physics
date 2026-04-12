#!/usr/bin/env python3
"""
g2_step49_overlap_normalization_gate.py

Step 49 of the g-2 rebuild.

What this script does
---------------------
1. Starts from the Stage-6/7 isotropic grouped-P2 overlap formulas
       D0 = K - B0 - Z0,
       D2 = -(M + B2 + Z2),
       D4 = -(B4 + Z4),
   together with
       P0 = N0 / D0.
2. Rewrites the universal quadrupole normalization condition as the exact reduced
   overlap-integral theorem gate
       mhat0^2 N0 / (K - B0 - Z0) = 54 G c_s^5 / (5 a^5 c^5).
3. Pulls the constant-prefactor branch conditions P2 = P4 = 0 down to the overlap
   level:
       N2 = -2 (M + B2 + Z2) N0 / (K - B0 - Z0),
       N4 = [ (M + B2 + Z2)^2 - 2 (K - B0 - Z0)(B4 + Z4) ] N0 / (K - B0 - Z0)^2.
4. Derives the first weak-axisymmetric anisotropy diagnostic:
   if
       D_A = D0 + eps lambda_A D1,
       N_A = N0 + eps lambda_A N1,
   with
       lambda = (1,1/2,-1),
   then the normalization defects satisfy
       b_P = 3 a_P.
   This is the direct reduced-PDE failure test for deviations from the isotropic
   branch.

Interpretation
--------------
The next honest theorem gate is no longer an abstract 'solve more PDE' request.
It is the explicit overlap-integral system above: compute the actual branch
overlaps, check isotropy, test the single normalization ratio, and then see whether
any weak anisotropy obeys the universal axisymmetric law b = 3a.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STEP 49A — EXACT ISOTROPIC OVERLAP-LEVEL BUNDLE COEFFICIENTS")

K, M = sp.symbols("K M", real=True)
B0, B2, B4 = sp.symbols("B0 B2 B4", real=True)
Z0, Z2, Z4 = sp.symbols("Z0 Z2 Z4", real=True)
N0, N2, N4 = sp.symbols("N0 N2 N4", real=True)
mhat0, G, a, c, c_s = sp.symbols("mhat0 G a c c_s", positive=True, real=True)

D0 = sp.simplify(K - B0 - Z0)
D2 = sp.simplify(-(M + B2 + Z2))
D4 = sp.simplify(-(B4 + Z4))
P0 = sp.simplify(N0 / D0)

print("D0 =", D0)
print("D2 =", D2)
print("D4 =", D4)
print("P0 =", P0)


banner("STEP 49B — THE EXACT OVERLAP-INTEGRAL NORMALIZATION GATE")

target_ratio = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
norm_gate = sp.Eq(sp.simplify(mhat0**2 * N0 / D0), target_ratio)

print("Normalization gate:")
sp.pprint(norm_gate)

# Equivalent cross-multiplied theorem gate.
norm_cross = sp.simplify(mhat0**2 * N0 - target_ratio * D0)
print("Cross-multiplied residual =", norm_cross)

print("\nInterpretation:")
print("  On the natural isotropic branch the full moving-throat theorem gate is the")
print("  single ratio mhat0^2 N0 / (K - B0 - Z0).")
print("  The Maxwell/mixed outgoing-transfer weight N0 competes directly against the")
print("  conservative dressed wall stiffness K - B0 - Z0.")


banner("STEP 49C — CONSTANT-PREFACTOR CONDITIONS IN PURE OVERLAP VARIABLES")

N2_overlap = sp.simplify(2 * D2 * N0 / D0)
N4_overlap = sp.simplify((D2**2 + 2 * D0 * D4) * N0 / D0**2)

print("N2 overlap condition =", N2_overlap)
print("N4 overlap condition =", N4_overlap)

expect_zero(
    "N2 overlap form",
    sp.simplify(N2_overlap + 2 * (M + B2 + Z2) * N0 / (K - B0 - Z0)),
)
expect_zero(
    "N4 overlap form",
    sp.simplify(
        N4_overlap
        - ((M + B2 + Z2) ** 2 - 2 * (K - B0 - Z0) * (B4 + Z4)) * N0 / (K - B0 - Z0) ** 2
    ),
)

print("\nSo if the minimal constant-prefactor outgoing branch is correct, the higher")
print("transfer moments are not independent. They are slaved to the conservative")
print("overlap data exactly by the two formulas printed above.")


banner("STEP 49D — FIRST WEAK-AXISYMMETRIC ANISOTROPY DIAGNOSTIC")

eps = sp.symbols("eps", real=True)
D1, N1 = sp.symbols("D1 N1", real=True)

lam20 = sp.Integer(1)
lam21 = sp.Rational(1, 2)
lam22 = -sp.Integer(1)

P20 = sp.simplify(sp.series((N0 + eps * lam20 * N1) / (D0 + eps * lam20 * D1), eps, 0, 2).removeO())
P21 = sp.simplify(sp.series((N0 + eps * lam21 * N1) / (D0 + eps * lam21 * D1), eps, 0, 2).removeO())
P22 = sp.simplify(sp.series((N0 + eps * lam22 * N1) / (D0 + eps * lam22 * D1), eps, 0, 2).removeO())

Pbar = sp.simplify((P20 + 2 * P21 + 2 * P22) / 5)
aP = sp.simplify((2 * P20 - P21 - P22) / 10)
bP = sp.simplify((P21 - P22) / 2)

P1 = sp.simplify((N1 * D0 - N0 * D1) / D0**2)

print("P20 =", P20)
print("P21 =", P21)
print("P22 =", P22)
print("Pbar =", Pbar)
print("aP =", aP)
print("bP =", bP)
print("P1 =", P1)

expect_zero("isotropic average stays P0", sp.simplify(Pbar - N0 / D0))
expect_zero("aP law", sp.simplify(aP - eps * P1 / 4))
expect_zero("bP law", sp.simplify(bP - 3 * eps * P1 / 4))
expect_zero("axisymmetric normalization defect law", sp.simplify(bP - 3 * aP))

print("\nInterpretation:")
print("  If a future PDE output shows weak grouped-lane normalization anisotropy,")
print("  the pure axisymmetric quadrupole branch must satisfy b_P = 3 a_P.")
print("  Any violation means the branch is not a simple weak axisymmetric deformation")
print("  of the isotropic moving-throat bundle.")


banner("FINAL STEP-49 LEDGER")
print("Exact isotropic overlap-level theorem gate:")
print("  mhat0^2 N0 / (K - B0 - Z0) = 54 G c_s^5 / (5 a^5 c^5)")
print()
print("Minimal constant-prefactor overlap constraints:")
print("  N2 = -2 (M + B2 + Z2) N0 / (K - B0 - Z0)")
print("  N4 = [ (M + B2 + Z2)^2 - 2 (K - B0 - Z0)(B4 + Z4) ] N0 / (K - B0 - Z0)^2")
print()
print("Weak axisymmetric normalization diagnostic:")
print("  b_P = 3 a_P")
print()
print("So the next honest PDE-facing falsification test is:")
print("  1. compute the actual overlap data Bn, Zn, Nn;")
print("  2. check whether the grouped anisotropy defects vanish on the natural branch;")
print("  3. test the single ratio mhat0^2 N0 / (K - B0 - Z0);")
print("  4. if weak anisotropy survives, verify whether it obeys the universal law b = 3 a.")
