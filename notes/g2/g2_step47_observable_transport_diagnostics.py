#!/usr/bin/env python3
"""
g2_step47_observable_transport_diagnostics.py

Step 47 of the g-2 rebuild.

What this script does
---------------------
1. Starts from the Step-46 outgoing target surface
       K0 = 27 Gamma5 s^5,
       K2 = 3 Gamma5 s^3,
       K4 = (4/3) Gamma5 s,
       Gamma5 = const,
   where s = c_s / a.
2. Inserts the adiabatic anomaly-only retuning
       s -> exp(ell/5) s
   that follows from
       delta ln a = 0,
       delta ln c_s = ell/5.
3. Derives the exact finite flow
       K0 -> exp(ell) K0,
       K2 -> exp(3 ell/5) K2,
       K4 -> exp(ell/5) K4,
       Gamma5 -> Gamma5.
4. Derives inversion formulas for ell from any even outgoing coefficient.
5. Derives the differential/logarithmic slope laws among observable bundle slots.
6. Evaluates the electron-point ratios numerically.

Interpretation
--------------
If the true moving-throat branch selects the minimal isotropic constant-prefactor
outgoing continuation, then the adiabatic anomaly does not deform the outgoing
bundle arbitrarily. It moves it along one exact one-parameter flow with a fixed
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


banner("STEP 47A — EXACT FINITE FLOW ON THE OUTGOING TARGET SURFACE")

Gamma5, s, ell = sp.symbols("Gamma_5 s ell", positive=True, real=True)

K0 = sp.simplify(27 * Gamma5 * s**5)
K2 = sp.simplify(3 * Gamma5 * s**3)
K4 = sp.simplify(sp.Rational(4, 3) * Gamma5 * s)

K0_flow = sp.simplify(K0.subs(s, s * sp.exp(ell / 5)))
K2_flow = sp.simplify(K2.subs(s, s * sp.exp(ell / 5)))
K4_flow = sp.simplify(K4.subs(s, s * sp.exp(ell / 5)))
Gamma5_flow = sp.simplify(Gamma5.subs(s, s * sp.exp(ell / 5)))

expect_zero("K0 flow", sp.simplify(K0_flow - sp.exp(ell) * K0))
expect_zero("K2 flow", sp.simplify(K2_flow - sp.exp(3 * ell / 5) * K2))
expect_zero("K4 flow", sp.simplify(K4_flow - sp.exp(ell / 5) * K4))
expect_zero("Gamma5 flow", sp.simplify(Gamma5_flow - Gamma5))

print("K0 ->", K0_flow)
print("K2 ->", K2_flow)
print("K4 ->", K4_flow)
print("Gamma5 ->", Gamma5_flow)


banner("STEP 47B — INVERSION OF THE ANOMALY PARAMETER FROM ANY EVEN SLOT")

K0s, K2s, K4s = sp.symbols("K0_star K2_star K4_star", positive=True, real=True)

K0_new = sp.simplify(sp.exp(ell) * K0s)
K2_new = sp.simplify(sp.exp(3 * ell / 5) * K2s)
K4_new = sp.simplify(sp.exp(ell / 5) * K4s)

# Exact reconstruction identities (written as zero checks after substitution).
ell_from_K0 = sp.log(K0_new / K0s)
ell_from_K2 = sp.Rational(5, 3) * sp.log(K2_new / K2s)
ell_from_K4 = 5 * sp.log(K4_new / K4s)

expect_zero("ell from K0 = ell", sp.simplify(ell_from_K0 - ell))
expect_zero("ell from K2 = ell", sp.simplify(ell_from_K2 - ell))
expect_zero("ell from K4 = ell", sp.simplify(ell_from_K4 - ell))

print("ell = ln(K0/K0*)")
print("ell = (5/3) ln(K2/K2*)")
print("ell = 5 ln(K4/K4*)")


banner("STEP 47C — DIFFERENTIAL / LOG-SLOPE LAWS")

eps = sp.symbols("eps", real=True)

K0_eps = sp.simplify(K0.subs(s, s * sp.exp(eps * ell / 5)))
K2_eps = sp.simplify(K2.subs(s, s * sp.exp(eps * ell / 5)))
K4_eps = sp.simplify(K4.subs(s, s * sp.exp(eps * ell / 5)))

dlnK0 = sp.simplify(sp.diff(sp.log(K0_eps), eps).subs(eps, 0))
dlnK2 = sp.simplify(sp.diff(sp.log(K2_eps), eps).subs(eps, 0))
dlnK4 = sp.simplify(sp.diff(sp.log(K4_eps), eps).subs(eps, 0))

expect_zero("d ln K0 - ell", dlnK0 - ell)
expect_zero("d ln K2 - 3 ell/5", dlnK2 - 3 * ell / 5)
expect_zero("d ln K4 - ell/5", dlnK4 - ell / 5)

expect_zero("observable slope law 1", sp.simplify(dlnK0 - 5 * dlnK4))
expect_zero("observable slope law 2", sp.simplify(dlnK2 - 3 * dlnK4))
expect_zero("observable slope law 3", sp.simplify(dlnK0 - sp.Rational(5, 3) * dlnK2))

print("d ln K0 =", dlnK0)
print("d ln K2 =", dlnK2)
print("d ln K4 =", dlnK4)
print("d ln Gamma5 = 0")


banner("STEP 47D — ELECTRON-POINT NUMERICS")

Lambda1_num = sp.Float("0.279605891931464")
f_num = sp.Float("0.001161409732093")
ell_num = sp.N(sp.log(1 + Lambda1_num * f_num), 20)

ratio_K0 = sp.N(sp.exp(ell_num), 20)
ratio_K2 = sp.N(sp.exp(3 * ell_num / 5), 20)
ratio_K4 = sp.N(sp.exp(ell_num / 5), 20)

print("ell_electron =", ell_num)
print("K0/K0*       =", ratio_K0)
print("K2/K2*       =", ratio_K2)
print("K4/K4*       =", ratio_K4)
print("Gamma5/Gamma5* = 1")
print("K0 ppm shift =", sp.N((ratio_K0 - 1) * 10**6, 18))
print("K2 ppm shift =", sp.N((ratio_K2 - 1) * 10**6, 18))
print("K4 ppm shift =", sp.N((ratio_K4 - 1) * 10**6, 18))

banner("STEP 47E — FINAL LEDGER")

print("Finite-flow law:")
print("  K0     -> exp(ell) K0")
print("  K2     -> exp(3 ell/5) K2")
print("  K4     -> exp(ell/5) K4")
print("  Gamma5 -> Gamma5")
print()
print("Inversion law:")
print("  ell = ln(K0/K0*) = (5/3) ln(K2/K2*) = 5 ln(K4/K4*)")
print()
print("Observable slope laws:")
print("  d ln K0 = 5 d ln K4")
print("  d ln K2 = 3 d ln K4")
print("  d ln K0 = (5/3) d ln K2")
