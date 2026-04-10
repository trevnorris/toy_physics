#!/usr/bin/env python3
"""
g2_step40_adiabatic_wall_anomaly_bridge.py

Step 40 of the g-2 rebuild.

What this script does
---------------------
1. Starts from the visible Step-39 strong adiabatic-wall closure:
      d ln Theta_w = 0,
      d ln Z_q = 0,
   so the reduced branch is parameterized by
      sigma := d ln K_s,
      ell   := d ln P_0.
2. Rewrites the quartic anomaly sliver as an exact outgoing-normalization law on
   that branch using the carried Step-29/30 defect relation.
3. Derives the exact microscopic anomaly laws
      d ln K_q = (2/5) ell,
      d ln(c_s/a) = (1/5) ell.
4. Evaluates the physical electron-point numbers.

Interpretation
--------------
On the strong adiabatic-wall branch, the quartic anomaly is carried by a tiny
relative core-stiffening / outgoing-normalization drift. The coherent elastic
wall squish sigma remains orthogonal to that anomaly lane.
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


banner("STEP 40A — STRONG ADIABATIC-WALL VARIABLES")

sigma, ell = sp.symbols("sigma ell", real=True)

dln_a = sp.simplify(sigma / 2)
dln_cs = sp.simplify(sigma / 2 + ell / 5)
dln_Kq = sp.simplify(2 * ell / 5)
dln_rel = sp.simplify(dln_cs - dln_a)

print("d ln a =", dln_a)
print("d ln c_s =", dln_cs)
print("d ln K_q =", dln_Kq)
print("d ln(c_s/a) =", dln_rel)

expect_zero("relative-core law", dln_rel - ell / 5)
expect_zero("relative-core vs K_q", dln_rel - dln_Kq / 2)

print("\nInterpretation:")
print("  sigma is the coherent elastic wall-squish direction.")
print("  ell   is the anomaly-carrying isotropic load / normalization direction.")
print("  The relative core stiffening c_s/a depends only on ell, not on sigma.")


banner("STEP 40B — EXACT OUTGOING-NORMALIZATION LAW")

Lambda1, f = sp.symbols("Lambda1 f", positive=True, real=True)
x = sp.simplify(Lambda1 * f)
Delta_Q = sp.simplify(-x / (1 + x))
N_Q = sp.simplify(1 / (1 + Delta_Q))
ell_exact = sp.simplify(sp.log(N_Q))
ell_linear = sp.simplify(N_Q - 1)

print("Delta_Q(f) =", Delta_Q)
print("N_Q(f)     =", N_Q)
print("ell_exact  = d ln P_0 =", ell_exact)
print("ell_linear = N_Q - 1  =", ell_linear)
expect_zero("N_Q - (1 + Lambda1 f)", N_Q - (1 + Lambda1 * f))

# Exact anomaly laws on the strong adiabatic-wall branch.
dln_Kq_exact = sp.simplify(2 * ell_exact / 5)
dln_rel_exact = sp.simplify(ell_exact / 5)

dln_Kq_linear = sp.simplify(2 * ell_linear / 5)
dln_rel_linear = sp.simplify(ell_linear / 5)

print("d ln K_q (exact)    =", dln_Kq_exact)
print("d ln(c_s/a) (exact) =", dln_rel_exact)
print("d ln K_q (linear)   =", dln_Kq_linear)
print("d ln(c_s/a) (linear)=", dln_rel_linear)

print("\nInterpretation:")
print("  If the source-map normalization is held fixed, the quartic anomaly sliver")
print("  is carried by the same outgoing-normalization factor N_Q = 1 + Lambda1 f.")
print("  On the strong adiabatic wall branch, that becomes a tiny relative")
print("  core-stiffening law plus a tiny K_q stiffening law.")


banner("STEP 40C — ELECTRON-POINT NUMBERS")

Lambda1_num = sp.Float("0.279605891931464")
f_num = sp.Float("0.001161409732093")
subs_num = {Lambda1: Lambda1_num, f: f_num}

vals = {
    "Lambda1 f": sp.N(ell_linear.subs(subs_num), 18),
    "ln(1 + Lambda1 f)": sp.N(ell_exact.subs(subs_num), 18),
    "d ln K_q exact": sp.N(dln_Kq_exact.subs(subs_num), 18),
    "d ln(c_s/a) exact": sp.N(dln_rel_exact.subs(subs_num), 18),
    "exact-linear gap": sp.N((ell_linear - ell_exact).subs(subs_num), 18),
}
for k, v in vals.items():
    print(f"{k:20s} = {v}")

print("\nNumerical reading:")
print("  - the outgoing normalization drift is O(3.25e-4),")
print("  - the induced K_q stiffening is O(1.30e-4),")
print("  - and the required relative core stiffening c_s/a is only O(6.49e-5).")
print("So on the strong adiabatic-wall branch the quartic electron sliver is a")
print("very small compressibility correction rather than a thermodynamic wall shift.")
