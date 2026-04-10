#!/usr/bin/env python3
"""
g2_step46_outgoing_bundle_target_surface.py

Step 46 of the g-2 rebuild.

What this script does
---------------------
1. Starts from the moving-throat constant-prefactor outgoing branch
       K0 = P0,
       K2 = P0 a^2 /(9 c_s^2),
       K4 = 4 P0 a^4 /(81 c_s^4),
       Gamma5 = P0 a^5 /(27 c_s^5),
   together with the universal normalization target
       P0 = 54 G c_s^5 /(5 a^5 c^5).
2. Eliminates P0 and derives the exact outgoing bundle
       K0 = 54 G s^5 /(5 c^5),
       K2 = 6 G s^3 /(5 c^5),
       K4 = 8 G s /(15 c^5),
       Gamma5 = 2 G /(5 c^5),
   where s = c_s / a.
3. Eliminates the remaining scale ratio s and derives exact algebraic target-surface
   relations among the observable bundle coefficients:
       K2^2 = (1/4) K0 K4,
       K0 = 4 K2^2 / K4,
       K2 = 81 K4^3 / (64 Gamma5^2),
       K0 = 6561 K4^5 / (1024 Gamma5^4).

Interpretation
--------------
If the true moving-throat branch selects the minimal constant-prefactor isotropic
outgoing continuation, then the full adiabatic outgoing bundle does not live in a
generic four-parameter space. It lies on a one-parameter algebraic leaf with a
universal odd quadrupole coefficient.
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


banner("STEP 46A — CONSTANT-PREFACTOR OUTGOING BUNDLE ON THE UNIVERSAL TARGET CURVE")

G, c, a, c_s = sp.symbols("G c a c_s", positive=True, real=True)

P0 = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
K0 = sp.simplify(P0)
K2 = sp.simplify(P0 * a**2 / (9 * c_s**2))
K4 = sp.simplify(4 * P0 * a**4 / (81 * c_s**4))
Gamma5 = sp.simplify(P0 * a**5 / (27 * c_s**5))

print("P0       =", P0)
print("K0       =", K0)
print("K2       =", K2)
print("K4       =", K4)
print("Gamma5   =", Gamma5)

expect_zero("Gamma5 - 2G/(5 c^5)", sp.simplify(Gamma5 - 2 * G / (5 * c**5)))

banner("STEP 46B — REDUCE TO THE SINGLE SCALE RATIO s = c_s / a")

Gamma5_sym, s = sp.symbols("Gamma_5 s", positive=True, real=True)

K0_s = sp.simplify(K0.subs({G: 5 * c**5 * Gamma5_sym / 2, c_s: s * a}))
K2_s = sp.simplify(K2.subs({G: 5 * c**5 * Gamma5_sym / 2, c_s: s * a}))
K4_s = sp.simplify(K4.subs({G: 5 * c**5 * Gamma5_sym / 2, c_s: s * a}))
Gamma5_s = sp.simplify(Gamma5.subs({G: 5 * c**5 * Gamma5_sym / 2}))

print("K0(s,Gamma5)     =", K0_s)
print("K2(s,Gamma5)     =", K2_s)
print("K4(s,Gamma5)     =", K4_s)
print("Gamma5(s)        =", Gamma5_s)

expect_zero("K0 - 27 Gamma5 s^5", K0_s - 27 * Gamma5_sym * s**5)
expect_zero("K2 - 3 Gamma5 s^3", K2_s - 3 * Gamma5_sym * s**3)
expect_zero("K4 - (4/3) Gamma5 s", K4_s - sp.Rational(4, 3) * Gamma5_sym * s)
expect_zero("Gamma5 invariant in s-parameterization", Gamma5_s - Gamma5_sym)

banner("STEP 46C — ELIMINATE s: EXACT ALGEBRAIC TARGET SURFACE")

inv1 = sp.simplify(K2_s**2 - K0_s * K4_s / 4)
inv2 = sp.simplify(K0_s - 4 * K2_s**2 / K4_s)
inv3 = sp.simplify(K2_s - 81 * K4_s**3 / (64 * Gamma5_sym**2))
inv4 = sp.simplify(K0_s - 6561 * K4_s**5 / (1024 * Gamma5_sym**4))

expect_zero("K2^2 - (1/4) K0 K4", inv1)
expect_zero("K0 - 4 K2^2 / K4", inv2)
expect_zero("K2 - 81 K4^3 / (64 Gamma5^2)", inv3)
expect_zero("K0 - 6561 K4^5 / (1024 Gamma5^4)", inv4)

print("\nReading:")
print("  The constant-prefactor outgoing bundle is a one-parameter leaf.")
print("  Once Gamma5 is fixed to its universal value, the three even slots are")
print("  forced onto the algebraic surface K2^2 = (1/4) K0 K4.")

banner("STEP 46D — OPTIONAL NUMERICAL CHECK IN TERMS OF G AND c")

G_num = sp.Symbol("G_num", positive=True)
c_num = sp.Symbol("c_num", positive=True)

# Pure symbolic check: substitute the universal Gamma5 back in.
Gamma5_univ = sp.simplify(2 * G_num / (5 * c_num**5))
K0_univ = sp.simplify(27 * Gamma5_univ * s**5)
K2_univ = sp.simplify(3 * Gamma5_univ * s**3)
K4_univ = sp.simplify(sp.Rational(4, 3) * Gamma5_univ * s)

expect_zero("bundle relation with universal Gamma5",
            sp.simplify(K2_univ**2 - K0_univ * K4_univ / 4))

print("\nFinal theorem ledger of Step 46:")
print("  K0     = 54 G s^5 / (5 c^5)")
print("  K2     = 6 G s^3 / (5 c^5)")
print("  K4     = 8 G s / (15 c^5)")
print("  Gamma5 = 2 G / (5 c^5)")
print("  with s = c_s / a,")
print("  and therefore K2^2 = (1/4) K0 K4 exactly.")
