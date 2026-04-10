#!/usr/bin/env python3
"""
g2_step48_bundle_target_pullback.py

Step 48 of the g-2 rebuild.

What this script does
---------------------
1. Starts from the grouped-P2 bundle formulas
       P0 = N0 / D0,
       P2 = (D0 N2 - 2 D2 N0) / D0^2,
       P4 = (D0^2 N4 - 2 D0(D2 N2 + D4 N0) + 3 D2^2 N0) / D0^3.
2. Solves the exact constant-prefactor conditions
       P2 = 0,
       P4 = 0,
   for N2 and N4 in terms of the conservative bundle moments D0,D2,D4 and N0.
3. Inserts the compact outgoing l=2 fingerprint
       A = a^2 / (9 c_s^2),
       B = 4 a^4 / (81 c_s^4) = 4 A^2,
       G5 = a^5 / (27 c_s^5),
   and proves that the observable target surface
       K2^2 = (1/4) K0 K4
   is automatic once the constant-prefactor microscopic closure is imposed.
4. Rewrites the adiabatic anomaly transport directly at the microscopic port level:
       P0 -> e^ell P0,
       A  -> e^{-2 ell/5} A,
       B  -> e^{-4 ell/5} B,
       G5 -> e^{-ell} G5,
   which reproduces the Step-47 observable bundle flow and explains immediately why
       Gamma5
   stays invariant.

Interpretation
--------------
The minimal outgoing observable target surface is not an extra independent condition.
It is exactly the image of:
  (i) the constant-prefactor microscopic bundle closure, and
  (ii) the compact outgoing l=2 fingerprint.
This is the cleanest pullback of the Step-46/47 outgoing bundle surface to the
microscopic grouped-P2 bundle.
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


banner("STEP 48A — EXACT BUNDLE-TO-OUTGOING MAP")

D0, D2, D4, N0, N2, N4 = sp.symbols("D0 D2 D4 N0 N2 N4", nonzero=True, real=True)
A, B, G5 = sp.symbols("A B G5", positive=True, real=True)

P0 = sp.simplify(N0 / D0)
P2 = sp.simplify((D0 * N2 - 2 * D2 * N0) / D0**2)
P4 = sp.simplify((D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3)

K0 = sp.simplify(P0)
K2 = sp.simplify(P2 + A * P0)
K4 = sp.simplify(P4 + A * P2 + B * P0)
Gamma5 = sp.simplify(G5 * P0)

print("P0 =", P0)
print("P2 =", P2)
print("P4 =", P4)
print("K0 =", K0)
print("K2 =", K2)
print("K4 =", K4)
print("Gamma5 =", Gamma5)


banner("STEP 48B — CONSTANT-PREFACTOR MICROSCOPIC CLOSURE")

N2_sol = sp.solve(sp.Eq(P2, 0), N2)[0]
N4_sol = sp.solve(sp.Eq(sp.simplify(P4.subs(N2, N2_sol)), 0), N4)[0]

print("N2 solution =", N2_sol)
print("N4 solution =", N4_sol)

expect_zero("P2 vanishes on N2_sol", sp.simplify(P2.subs(N2, N2_sol)))
expect_zero("P4 vanishes on N2_sol,N4_sol", sp.simplify(P4.subs({N2: N2_sol, N4: N4_sol})))

u2 = sp.symbols("u2", real=True)
u4 = sp.symbols("u4", real=True)

u2_def = sp.simplify(-D2 / D0)
u4_def = sp.simplify((D2**2 - D0 * D4) / D0**2)

expect_zero("u2 definition consistency", sp.simplify(u2_def + D2 / D0))
expect_zero("u4 definition consistency", sp.simplify(u4_def - (D2**2 - D0 * D4) / D0**2))

N2_u = sp.simplify(N2_sol.subs(D2, -u2 * D0))
N4_u = sp.simplify(N4_sol.subs({D2: -u2 * D0, D4: (u2**2 - u4) * D0}))

print("N2 written in (u2,N0) form =", N2_u)
print("N4 written in (u2,u4,N0) form =", N4_u)

expect_zero("N2 = -2 u2 N0", sp.simplify(N2_u + 2 * u2 * N0))
expect_zero("N4 = (3u2^2 - 2u4) N0", sp.simplify(N4_u - (3 * u2**2 - 2 * u4) * N0))


banner("STEP 48C — THE STEP-46 TARGET SURFACE IS THE IMAGE OF THE MICROSCOPIC CLOSURE")

K0_cp = sp.simplify(K0.subs({N2: N2_sol, N4: N4_sol}))
K2_cp = sp.simplify(K2.subs({N2: N2_sol, N4: N4_sol}))
K4_cp = sp.simplify(K4.subs({N2: N2_sol, N4: N4_sol}))
Gamma5_cp = sp.simplify(Gamma5.subs({N2: N2_sol, N4: N4_sol}))

print("K0 on constant-prefactor branch =", K0_cp)
print("K2 on constant-prefactor branch =", K2_cp)
print("K4 on constant-prefactor branch =", K4_cp)
print("Gamma5 on constant-prefactor branch =", Gamma5_cp)

expect_zero("K2 reduces to A K0", sp.simplify(K2_cp - A * K0_cp))
expect_zero("K4 reduces to B K0", sp.simplify(K4_cp - B * K0_cp))
expect_zero("Gamma5 reduces to G5 K0", sp.simplify(Gamma5_cp - G5 * K0_cp))

expect_zero(
    "Step-46 algebraic surface from B = 4 A^2",
    sp.simplify((K2_cp**2 - sp.Rational(1, 4) * K0_cp * K4_cp).subs(B, 4 * A**2)),
)


banner("STEP 48D — INSERTING THE COMPACT OUTGOING l=2 FINGERPRINT")

a, c_s = sp.symbols("a c_s", positive=True, real=True)

A_port = sp.simplify(a**2 / (9 * c_s**2))
B_port = sp.simplify(4 * a**4 / (81 * c_s**4))
G5_port = sp.simplify(a**5 / (27 * c_s**5))

expect_zero("B_port = 4 A_port^2", sp.simplify(B_port - 4 * A_port**2))

K0_port = sp.simplify(K0_cp)
K2_port = sp.simplify(K2_cp.subs(A, A_port))
K4_port = sp.simplify(K4_cp.subs({A: A_port, B: B_port}))
Gamma5_port = sp.simplify(Gamma5_cp.subs(G5, G5_port))

print("A_port =", A_port)
print("B_port =", B_port)
print("G5_port =", G5_port)
print("K2/K0 on compact branch =", sp.simplify(K2_port / K0_port))
print("K4/K0 on compact branch =", sp.simplify(K4_port / K0_port))
print("Gamma5/K0 on compact branch =", sp.simplify(Gamma5_port / K0_port))


banner("STEP 48E — MICROSCOPIC PORT TRANSPORT ON THE ADIABATIC ANOMALY TRACK")

ell = sp.symbols("ell", real=True)
P0_star = sp.symbols("P0_star", positive=True, real=True)

P0_flow = sp.simplify(sp.exp(ell) * P0_star)
A_flow = sp.simplify(A_port.subs(c_s, c_s * sp.exp(ell / 5)))
B_flow = sp.simplify(B_port.subs(c_s, c_s * sp.exp(ell / 5)))
G5_flow = sp.simplify(G5_port.subs(c_s, c_s * sp.exp(ell / 5)))

expect_zero("A scales as exp(-2 ell/5)", sp.simplify(A_flow - sp.exp(-2 * ell / 5) * A_port))
expect_zero("B scales as exp(-4 ell/5)", sp.simplify(B_flow - sp.exp(-4 * ell / 5) * B_port))
expect_zero("G5 scales as exp(-ell)", sp.simplify(G5_flow - sp.exp(-ell) * G5_port))

K0_flow = sp.simplify(P0_flow)
K2_flow = sp.simplify(A_flow * P0_flow)
K4_flow = sp.simplify(B_flow * P0_flow)
Gamma5_flow = sp.simplify(G5_flow * P0_flow)

expect_zero("K0 flow", sp.simplify(K0_flow - sp.exp(ell) * P0_star))
expect_zero("K2 flow", sp.simplify(K2_flow - sp.exp(3 * ell / 5) * A_port * P0_star))
expect_zero("K4 flow", sp.simplify(K4_flow - sp.exp(ell / 5) * B_port * P0_star))
expect_zero("Gamma5 invariance", sp.simplify(Gamma5_flow - G5_port * P0_star))

print("P0 ->", P0_flow)
print("A  ->", A_flow)
print("B  ->", B_flow)
print("G5 ->", G5_flow)

print("\nInterpretation:")
print("  The adiabatic anomaly acts microscopically as:")
print("    - one outgoing-normalization rescaling P0 -> exp(ell) P0,")
print("    - plus one universal compact-port retuning driven by c_s -> c_s exp(ell/5).")
print("  Because G5 rescales as exp(-ell), the odd coefficient Gamma5 = G5 P0")
print("  stays exactly invariant even while the even outgoing bundle moves.")


banner("FINAL STEP-48 LEDGER")
print("Constant-prefactor microscopic bundle closure:")
print("  N2 = 2 D2 N0 / D0 = -2 u2 N0")
print("  N4 = (D2^2 + 2 D0 D4) N0 / D0^2 = (3u2^2 - 2u4) N0")
print()
print("Compact outgoing l=2 fingerprint:")
print("  A  = a^2 / (9 c_s^2)")
print("  B  = 4 a^4 / (81 c_s^4) = 4 A^2")
print("  G5 = a^5 / (27 c_s^5)")
print()
print("Observable surface as image of microscopic closure:")
print("  K0 = P0")
print("  K2 = A K0")
print("  K4 = B K0")
print("  Gamma5 = G5 K0")
print("  => K2^2 = (1/4) K0 K4")
print()
print("Adiabatic anomaly port transport:")
print("  P0 -> exp(ell) P0")
print("  A  -> exp(-2 ell/5) A")
print("  B  -> exp(-4 ell/5) B")
print("  G5 -> exp(-ell) G5")
print("  => Gamma5 invariant, while K0,K2,K4 follow the Step-47 flow.")
