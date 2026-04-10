#!/usr/bin/env python3
"""
5pn_stage4_axisymmetric_transport.py

Fourth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Introduces a weak-axisymmetric grouped-lane deformation with the exact Stage-7
   pattern
       lambda_(20) = 1,
       lambda_(21) = 1/2,
       lambda_(22) = -1.
2. Decomposes the grouped operator slopes D01,D21,D41 and transfer slope N01
   into microscopic wall / BdG / Maxwell-mixed bundles, following Stage 154.
3. Verifies the exact obstruction formulas
       K_1 := D21 + D01/9 = W_1 - B_1 - Z_1,
       G_1 := N01 - P0 D01 = -P0 K1 + P0 B01 + Nbundle_1.
4. Derives the physical weak-axisymmetric slopes of the grouped response and
   outgoing prefactor:
       u2^(1), u4^(1), P1/P0.
5. Verifies the hidden-even consistency law on the canonical compensated branch:
       u4^(1) = (8/9) u2^(1)
   iff
       D41 = (2/3) D21 + D01/27.
6. Shows that on the even-preserving branch u2^(1)=0 the whole conservative
   grouped response is transported by one static slope D01, while the remaining
   linear grouped 2.5PN defect collapses to one scalar loading mismatch
       Xi_load = N01/N0 - D01/D0 = P1/P0.

Interpretation
--------------
This script does not yet solve the true anisotropic moving-throat branch. It
turns the handoff's grouped-lane anisotropy problem into the exact quantities
that the next PDE/overlap computation must deliver: D01 and N01, or equivalently
u2^(1) and P1/P0, and on the even-preserving branch just Xi_load.
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


# ---------------------------------------------------------------------------
# I. Weak-axisymmetric grouped-lane setup
# ---------------------------------------------------------------------------

banner("I. WEAK-AXISYMMETRIC GROUPED-LANE SETUP")

eps = sp.symbols("eps", real=True)
lambda20 = sp.Integer(1)
lambda21 = sp.Rational(1, 2)
lambda22 = -sp.Integer(1)
lambdaA = sp.symbols("lambda_A", real=True)

print("lambda_(20) =", lambda20)
print("lambda_(21) =", lambda21)
print("lambda_(22) =", lambda22)

D0, D2, D4, N0 = sp.symbols("D0 D2 D4 N0", nonzero=True, real=True)
P0 = N0 / D0
u2 = -D2 / D0
u4 = (D2**2 - D0 * D4) / D0**2

print("Baseline grouped quantities:")
print("u2 =", u2)
print("u4 =", u4)
print("P0 =", P0)


# ---------------------------------------------------------------------------
# II. Microscopic slope decomposition (Stage 154 structure)
# ---------------------------------------------------------------------------

banner("II. MICROSCOPIC SLOPE DECOMPOSITION")

K1_wall, M1_wall = sp.symbols("K1_wall M1_wall", real=True)
B01, B21, B41 = sp.symbols("B01 B21 B41", real=True)
Z01, Z21, Z41 = sp.symbols("Z01 Z21 Z41", real=True)
N01 = sp.symbols("N01", real=True)

D01 = K1_wall - B01 - Z01
D21 = -(M1_wall + B21 + Z21)
D41 = -(B41 + Z41)

print("D01 =", D01)
print("D21 =", D21)
print("D41 =", D41)

W1 = sp.simplify(K1_wall / 9 - M1_wall)
Bcal1 = sp.simplify(B21 + B01 / 9)
Zcal1 = sp.simplify(Z21 + Z01 / 9)
Nbundle1 = sp.simplify(N01 + P0 * Z01)

K_obstruction_1 = sp.simplify(D21 + D01 / 9)
G_obstruction_1 = sp.simplify(N01 - P0 * D01)

expect_zero("K_obstruction_1 - (W1 - Bcal1 - Zcal1)", K_obstruction_1 - (W1 - Bcal1 - Zcal1))
expect_zero("G_obstruction_1 - (-P0*K1_wall + P0*B01 + Nbundle1)", G_obstruction_1 - (-P0 * K1_wall + P0 * B01 + Nbundle1))

print("\nW1       =", W1)
print("Bcal1    =", Bcal1)
print("Zcal1    =", Zcal1)
print("Nbundle1 =", Nbundle1)
print("K_obstruction_1 =", K_obstruction_1)
print("G_obstruction_1 =", G_obstruction_1)


# ---------------------------------------------------------------------------
# III. Physical weak-axisymmetric slopes of u2, u4, and P0
# ---------------------------------------------------------------------------

banner("III. PHYSICAL WEAK-AXISYMMETRIC SLOPES")

D0A = D0 + eps * lambdaA * D01
D2A = D2 + eps * lambdaA * D21
D4A = D4 + eps * lambdaA * D41
N0A = N0 + eps * lambdaA * N01

u2A = sp.expand(sp.series(-D2A / D0A, eps, 0, 2).removeO())
u4A = sp.expand(sp.series((D2A**2 - D0A * D4A) / D0A**2, eps, 0, 2).removeO())
P0A = sp.expand(sp.series(N0A / D0A, eps, 0, 2).removeO())

u2_1 = sp.simplify(u2A.coeff(eps, 1) / lambdaA)
u4_1 = sp.simplify(u4A.coeff(eps, 1) / lambdaA)
P1 = sp.simplify(P0A.coeff(eps, 1) / lambdaA)

print("u2^(A) =", u2A)
print("u4^(A) =", u4A)
print("P0^(A) =", P0A)

print("\nu2^(1) =", u2_1)
print("u4^(1) =", u4_1)
print("P1     =", P1)
print("P1/P0  =", sp.simplify(P1 / P0))

expect_zero("u2^(1) + (D21 + u2*D01)/D0", u2_1 + (D21 + u2 * D01) / D0)
expect_zero("P1 - (D0*N01 - N0*D01)/D0^2", P1 - (D0 * N01 - N0 * D01) / D0**2)
expect_zero("P1/P0 - (N01/N0 - D01/D0)", sp.simplify(P1 / P0) - (N01 / N0 - D01 / D0))
expect_zero("K_obstruction_1 + D0*u2^(1) - (1/9 - u2)*D01", K_obstruction_1 + D0 * u2_1 - (sp.Rational(1, 9) - u2) * D01)
expect_zero("G_obstruction_1 - D0*P1", G_obstruction_1 - D0 * P1)


# ---------------------------------------------------------------------------
# IV. Hidden-even consistency on the canonical compensated branch
# ---------------------------------------------------------------------------

banner("IV. HIDDEN-EVEN CONSISTENCY ON THE CANONICAL BRANCH")

canonical_subs = {
    D2: -D0 / 9,
    D4: -D0 / 27,
}

u2_1_can = sp.simplify(u2_1.subs(canonical_subs))
u4_1_can = sp.simplify(u4_1.subs(canonical_subs))

print("u2^(1) on the canonical branch =", u2_1_can)
print("u4^(1) on the canonical branch =", u4_1_can)

expect_zero("u2^(1)_canonical + (D21 + D01/9)/D0", u2_1_can + (D21 + D01 / 9) / D0)
expect_zero(
    "u4^(1)_canonical + (5*D01 + 18*D21 + 81*D41)/(81*D0)",
    u4_1_can + (5 * D01 + 18 * D21 + 81 * D41) / (81 * D0),
)

hidden_even_residual = sp.simplify(sp.expand(u4_1_can - sp.Rational(8, 9) * u2_1_can))
print("\nHidden-even residual u4^(1) - (8/9) u2^(1) =", hidden_even_residual)
hidden_even_law = sp.Rational(2, 3) * D21 + D01 / 27
print("Equivalent operator law: D41 =", hidden_even_law)
expect_zero(
    "hidden-even residual + (D41 - (2*D21/3 + D01/27))/D0",
    hidden_even_residual + (D41 - hidden_even_law) / D0,
)


# ---------------------------------------------------------------------------
# V. Even-preserving branch and collapse to Xi_load
# ---------------------------------------------------------------------------

banner("V. EVEN-PRESERVING BRANCH AND STATIC LOADING MISMATCH")

D21_even = -D01 / 9
print("D21 on the even-preserving branch =", D21_even)
expect_zero("D21 + D01/9", D21_even + D01 / 9)

D41_even = -D01 / 27
print("D41 on the even-preserving + hidden-even branch =", D41_even)
expect_zero("D41 + D01/27", D41_even + D01 / 27)

Xi_load = sp.simplify(N01 / N0 - D01 / D0)
expect_zero("Xi_load - P1/P0", Xi_load - sp.simplify(P1 / P0))
print("Xi_load =", Xi_load)

subbanner("Grouped-lane normalization defect pattern")
DeltaQ20 = sp.simplify(eps * lambda20 * Xi_load)
DeltaQ21 = sp.simplify(eps * lambda21 * Xi_load)
DeltaQ22 = sp.simplify(eps * lambda22 * Xi_load)
print("Delta_Q^(20) =", DeltaQ20)
print("Delta_Q^(21) =", DeltaQ21)
print("Delta_Q^(22) =", DeltaQ22)

# Check grouped defect law again in physical variables.
P20 = sp.simplify(P0 + DeltaQ20 * P0)
P21 = sp.simplify(P0 + DeltaQ21 * P0)
P22 = sp.simplify(P0 + DeltaQ22 * P0)

abarP = sp.simplify((P20 + 2 * P21 + 2 * P22) / 5 - P0)
aP = sp.simplify((2 * P20 - P21 - P22) / 10)
bP = sp.simplify((P21 - P22) / 2)

expect_zero("grouped prefactor trace defect", abarP)
expect_zero("aP - eps*P0*Xi_load/4", aP - eps * P0 * Xi_load / 4)
expect_zero("bP - 3*eps*P0*Xi_load/4", bP - 3 * eps * P0 * Xi_load / 4)
expect_zero("bP - 3 aP", bP - 3 * aP)

subbanner("Direct outlet amplitudes in physical variables")
sigma_star = sp.symbols("sigma_star", nonzero=True)
kappa1 = sp.simplify(-3 * (1 - sigma_star) * u2_1_can / sigma_star)
gamma1 = sp.simplify(-(1 - sigma_star) * Xi_load / (9 * sigma_star))
print("kappa_1 =", kappa1)
print("gamma_1 =", gamma1)


# ---------------------------------------------------------------------------
# VI. Final ledger
# ---------------------------------------------------------------------------

banner("VI. AXISYMMETRIC TRANSPORT LEDGER")
print("1. The grouped-lane obstruction pair is exactly")
print("      K_1 = D21 + D01/9,")
print("      G_1 = N01 - P0 D01,")
print("   with the exact microscopic decomposition")
print("      K_1 = W1 - Bcal1 - Zcal1,")
print("      G_1 = -P0 K1_wall + P0 B01 + Nbundle1.")
print("2. The physical weak-axisymmetric slopes are")
print("      u2^(1) = -(D21 + u2 D01)/D0,")
print("      P1/P0  = N01/N0 - D01/D0.")
print("3. On the canonical compensated branch, hidden-even consistency is")
print("      u4^(1) = (8/9) u2^(1)")
print("   iff")
print("      D41 = (2/3) D21 + D01/27.")
print("4. On the even-preserving branch u2^(1)=0, the conservative grouped response")
print("   is transported by one static slope D01 only:")
print("      D21 = -D01/9,   D41 = -D01/27.")
print("5. The remaining linear grouped 2.5PN defect collapses to one scalar loading")
print("   mismatch")
print("      Xi_load = N01/N0 - D01/D0 = P1/P0.")
print("6. Its grouped-lane signature is fixed:")
print("      (20,21,22) ~ (1, 1/2, -1),")
print("   equivalently b = 3 a.")
print("7. The next honest PDE step is now very narrow: compute D01 and N01 — or")
print("   directly Xi_load — from a primitive weak-axisymmetric deformation of the")
print("   finite-throat overlap model.")
