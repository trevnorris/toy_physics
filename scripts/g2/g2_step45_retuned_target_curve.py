#!/usr/bin/env python3
"""
g2_step45_retuned_target_curve.py

Step 45 of the g-2 rebuild.

What this script does
---------------------
1. Takes the Step-44 normalization-product bridge
       P0 = exp(ell) P0_target,
   with
       P0_target(a,c_s) = 54 G c_s^5 / (5 a^5 c^5),
   and the adiabatic anomaly-only drift law from Steps 40-41
       delta ln a = 0,
       delta ln c_s = ell/5.
2. Proves that the normalization-product rescaling from Step 44 is exactly the
   same as evaluating the universal moving-throat target curve at the retuned
   adiabatic branch values:
       P0(a, c_s e^{ell/5}) = exp(ell) P0(a,c_s).
3. Records the minimal grouped-P2 continuation if one also adopts the constant-
   isotropic-prefactor branch P2 = P4 = 0:
       K2 = P0 a^2 / (9 c_s^2),
       K4 = 4 P0 a^4 / (81 c_s^4).
4. Derives the finite adiabatic scaling laws
       P0 -> e^{ell},
       K2 -> e^{3 ell / 5},
       K4 -> e^{ell / 5},
   and the corresponding tangent drifts.
5. Verifies the exact transfer-moment conditions of the constant-prefactor branch
   in terms of the conservative grouped-response moments:
       N2 = -2 u2 N0,
       N4 = (3 u2^2 - 2 u4) N0.

Interpretation
--------------
On the adiabatic anomaly track, the moving-throat normalization product does not
leave the universal target curve. It slides along the same curve because the core
sound scale retunes by exactly the amount needed. The extra constant-prefactor
branch conditions are the minimal grouped-P2 completion of that picture, not yet a
fully forced theorem of the unsolved PDE.
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


banner("STEP 45A — THE ADIABATIC ANOMALY RIDES THE SAME UNIVERSAL NORMALIZATION CURVE")

G, a, c, c_s, ell = sp.symbols("G a c c_s ell", positive=True, real=True)

P0_target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
P0_retuned = sp.simplify(P0_target.subs(c_s, c_s * sp.exp(ell / 5)))
P0_scaled = sp.simplify(sp.exp(ell) * P0_target)

print("P0_target(a,c_s) =", P0_target)
print("P0_target(a,c_s e^{ell/5}) =", P0_retuned)
print("exp(ell) P0_target         =", P0_scaled)

expect_zero(
    "retuned target curve = Step-44 scaling",
    sp.simplify(P0_retuned - P0_scaled),
)

print("\nReading:")
print("  The Step-44 prefactor enhancement does not move the branch off the")
print("  universal moving-throat normalization law. It is exactly what one gets")
print("  by evaluating the same target product at the retuned adiabatic sound")
print("  scale c_s -> c_s exp(ell/5), with the wall radius fixed.")


banner("STEP 45B — MINIMAL CONSTANT-PREFACTOR GROUPED-P2 CONTINUATION")

K2, K4 = sp.symbols("K2 K4", positive=True, real=True)

K2_min = sp.simplify(P0_target * a**2 / (9 * c_s**2))
K4_min = sp.simplify(4 * P0_target * a**4 / (81 * c_s**4))

print("K2_min =", K2_min)
print("K4_min =", K4_min)

# Retuned branch values under c_s -> c_s e^{ell/5}, a fixed.
K2_retuned = sp.simplify(K2_min.subs(c_s, c_s * sp.exp(ell / 5)))
K4_retuned = sp.simplify(K4_min.subs(c_s, c_s * sp.exp(ell / 5)))

expect_zero("K2 retuned scaling", sp.simplify(K2_retuned - sp.exp(3 * ell / 5) * K2_min))
expect_zero("K4 retuned scaling", sp.simplify(K4_retuned - sp.exp(ell / 5) * K4_min))

print("K2_retuned =", K2_retuned)
print("K4_retuned =", K4_retuned)

print("\nInterpretation:")
print("  If one chooses the simplest isotropic grouped-P2 continuation P2 = P4 = 0,")
print("  then the even coefficients are slaved to the same retuned target curve:")
print("      P0  -> exp(ell),")
print("      K2  -> exp(3 ell/5),")
print("      K4  -> exp(ell/5).")


banner("STEP 45C — TANGENT DRIFTS OF THE MINIMAL GROUPED-P2 CONTINUATION")

eps = sp.symbols("eps", real=True)
P0_eps = sp.simplify(P0_target.subs(c_s, c_s * sp.exp(eps * ell / 5)))
K2_eps = sp.simplify(K2_min.subs(c_s, c_s * sp.exp(eps * ell / 5)))
K4_eps = sp.simplify(K4_min.subs(c_s, c_s * sp.exp(eps * ell / 5)))

dlnP0 = sp.simplify(sp.diff(sp.log(P0_eps), eps).subs(eps, 0))
dlnK2 = sp.simplify(sp.diff(sp.log(K2_eps), eps).subs(eps, 0))
dlnK4 = sp.simplify(sp.diff(sp.log(K4_eps), eps).subs(eps, 0))

print("delta ln P0 =", dlnP0)
print("delta ln K2 =", dlnK2)
print("delta ln K4 =", dlnK4)

expect_zero("tangent drift P0", sp.simplify(dlnP0 - ell))
expect_zero("tangent drift K2", sp.simplify(dlnK2 - 3 * ell / 5))
expect_zero("tangent drift K4", sp.simplify(dlnK4 - ell / 5))


banner("STEP 45D — CONSTANT-PREFACTOR TRANSFER-MOMENT CONDITIONS")

D0, D2, D4, N0, N2, N4, u2, u4 = sp.symbols("D0 D2 D4 N0 N2 N4 u2 u4", nonzero=True, real=True)

u2_def = sp.Eq(u2, -D2 / D0)
u4_def = sp.Eq(u4, (D2**2 - D0 * D4) / D0**2)

P2_expr = sp.simplify((D0 * N2 - 2 * D2 * N0) / D0**2)
P4_expr = sp.simplify((D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3)

N2_sol = sp.solve(sp.Eq(P2_expr, 0), N2)[0]
N4_sol = sp.solve(sp.Eq(P4_expr.subs(N2, N2_sol), 0), N4)[0]

N2_u = sp.simplify(N2_sol.subs({D2: -u2 * D0}))
N4_u = sp.simplify(N4_sol.subs({D2: -u2 * D0, D4: (u2**2 - u4) * D0}))

print("P2 =", P2_expr)
print("P4 =", P4_expr)
print("N2 on constant-prefactor branch =", N2_sol)
print("N4 on constant-prefactor branch =", N4_sol)
print("N2 in (u2,u4) language =", N2_u)
print("N4 in (u2,u4) language =", N4_u)

expect_zero("N2 constant-prefactor identity", sp.simplify(N2_u + 2 * u2 * N0))
expect_zero("N4 constant-prefactor identity", sp.simplify(N4_u - (3 * u2**2 - 2 * u4) * N0))

print("\nReading:")
print("  On the minimal constant-prefactor isotropic branch the outgoing transfer")
print("  moments are not arbitrary. They are slaved to the conservative grouped")
print("  response moments by")
print("      N2 = -2 u2 N0,")
print("      N4 = (3 u2^2 - 2 u4) N0.")


banner("STEP 45E — ELECTRON-POINT RATIOS")

Lambda1, f = sp.symbols("Lambda1 f", positive=True, real=True)
ell_f = sp.log(1 + Lambda1 * f)

P0_ratio = sp.simplify(sp.exp(ell_f))
K2_ratio = sp.simplify(sp.exp(3 * ell_f / 5))
K4_ratio = sp.simplify(sp.exp(ell_f / 5))

print("P0 ratio =", P0_ratio)
print("K2 ratio =", K2_ratio)
print("K4 ratio =", K4_ratio)

Lambda1_num = sp.Float("0.279605891931464")
f_num = sp.Float("0.001161409732093")

P0_ratio_num = sp.N(P0_ratio.subs({Lambda1: Lambda1_num, f: f_num}), 18)
K2_ratio_num = sp.N(K2_ratio.subs({Lambda1: Lambda1_num, f: f_num}), 18)
K4_ratio_num = sp.N(K4_ratio.subs({Lambda1: Lambda1_num, f: f_num}), 18)

print("P0 electron ratio =", P0_ratio_num)
print("K2 electron ratio =", K2_ratio_num)
print("K4 electron ratio =", K4_ratio_num)
print("P0 ppm shift      =", sp.N((P0_ratio_num - 1) * 10**6, 18))
print("K2 ppm shift      =", sp.N((K2_ratio_num - 1) * 10**6, 18))
print("K4 ppm shift      =", sp.N((K4_ratio_num - 1) * 10**6, 18))

print("\nFinal reading of Step 45:")
print("  - the adiabatic anomaly-only branch keeps the moving-throat normalization")
print("    law on the same universal target curve,")
print("  - the minimal constant-prefactor grouped-P2 continuation rescales")
print("        P0 by exp(ell), K2 by exp(3 ell/5), K4 by exp(ell/5),")
print("  - and, if that minimal branch is the one selected by the PDE, the")
print("    outgoing transfer moments must obey")
print("        N2 = -2 u2 N0,")
print("        N4 = (3 u2^2 - 2 u4) N0.")
