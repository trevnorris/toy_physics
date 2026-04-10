#!/usr/bin/env python3
"""
g2_step44_normalization_product_bridge.py

Step 44 of the g-2 rebuild.

What this script does
---------------------
1. Takes the exact reduced factorization from Step 28,
       mhat_0^2 * chi_Q * N_Q = 1,
   together with the Step-43 pure-odd normalized DtN law
       chi_Q = exp(-ell).
2. Restricts to the natural point-particle source-map branch mhat_0 = 1 and
   proves the reciprocal anomaly law
       N_Q = exp(+ell).
3. Uses the Step-27 normalization-defect scaling and the moving-throat Stage-5
   normalization product
       Gamma_5 = P_0 a^5 / (27 c_s^5),
       P_0^target = 54 G c_s^5 / (5 a^5 c^5),
   to show
       P_0 = N_Q P_0^target = exp(ell) P_0^target
   on the natural source-map branch.
4. Rewrites this equivalently as
       P_0 * chi_Q = P_0^target,
   i.e. the pure-odd outlet deformation and the grouped-P2 prefactor rescaling are
   exact reciprocals of the same quartic branch motion.
5. Evaluates the finite electron-point scale factor.

Interpretation
--------------
Step 43 did not just produce a nice DtN picture. It also fixes how the adiabatic
pure-odd outlet must feed into the moving-throat grouped-P2 normalization product.
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


banner("STEP 44A — PURE ODD DTN DEFECT AND THE REDUCED NORMALIZATION FACTORIZATION")

ell, chi_Q, N_Q, mhat0 = sp.symbols("ell chi_Q N_Q mhat0", positive=True, real=True)

chi_pure_odd = sp.exp(-ell)
factorization = sp.Eq(mhat0**2 * chi_pure_odd * N_Q, 1)
NQ_exact = sp.solve(factorization, N_Q)[0]
NQ_source_map = sp.simplify(NQ_exact.subs(mhat0, 1))

print("chi_Q(ell) =", chi_pure_odd)
print("exact factorization =", factorization)
print("N_Q exact =", NQ_exact)
print("N_Q on natural source-map branch =", NQ_source_map)

expect_zero("source-map reciprocal law", sp.simplify(NQ_source_map - sp.exp(ell)))

print("\nReading:")
print("  On the natural point-particle source-map branch, the pure-odd outlet")
print("  deformation chi_Q = exp(-ell) forces the grouped normalization defect")
print("  to move in the opposite direction:")
print("      N_Q = exp(+ell).")


banner("STEP 44B — THE MOVING-THROAT GROUPED-P2 NORMALIZATION PRODUCT")

G, a, c, c_s, P0 = sp.symbols("G a c c_s P0", positive=True, real=True)
Gamma5_target = sp.simplify(2 * G / (5 * c**5))
P0_target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))

# Stage-27 normalization-defect scaling of the odd grouped-P2 coefficient.
Gamma5_actual = sp.simplify(N_Q * Gamma5_target)

# Stage-5 moving-throat branch relation.
Gamma5_from_P0 = sp.simplify(P0 * a**5 / (27 * c_s**5))

# Solve P0 in terms of N_Q and the target product.
P0_exact = sp.solve(sp.Eq(Gamma5_from_P0, Gamma5_actual), P0)[0]
P0_source_map = sp.simplify(P0_exact.subs({N_Q: NQ_source_map}))

print("Gamma5_target =", Gamma5_target)
print("P0_target     =", P0_target)
print("Gamma5_from_P0 =", Gamma5_from_P0)
print("P0 exact =", P0_exact)
print("P0 on natural source-map branch =", P0_source_map)

expect_zero("P0 = N_Q * P0_target", sp.simplify(P0_exact - N_Q * P0_target))
expect_zero("source-map product bridge", sp.simplify(P0_source_map - sp.exp(ell) * P0_target))

print("\nReading:")
print("  On the natural source-map branch the grouped-P2 static prefactor is")
print("  renormalized by the exact inverse of the pure-odd outlet factor:")
print("      P0 = exp(+ell) P0_target.")


banner("STEP 44C — EXACT RECIPROCAL LAW BETWEEN THE ODD DTN DEFORMATION AND THE STATIC PREFactor")

reciprocal_law = sp.simplify(P0_source_map * chi_pure_odd)
print("P0 * chi_Q =", reciprocal_law)
expect_zero("reciprocal anomaly bridge", sp.simplify(reciprocal_law - P0_target))

print("\nEquivalent statement:")
print("  The adiabatic quartic branch may be read either as")
print("      chi_Q = exp(-ell)")
print("  in normalized DtN language, or as")
print("      P0 = exp(+ell) P0_target")
print("  in grouped-P2 normalization language.")
print("  Those are exact reciprocals on the natural source-map branch.")


banner("STEP 44D — FINITE ELECTRON-POINT SCALE FACTOR")

Lambda1, f = sp.symbols("Lambda1 f", positive=True, real=True)
ell_from_f = sp.log(1 + Lambda1 * f)
scale_factor = sp.simplify(sp.exp(ell_from_f))

print("ell(f)        =", ell_from_f)
print("exp(ell(f))   =", scale_factor)
expect_zero("exp(log(1+x)) identity", sp.simplify(scale_factor - (1 + Lambda1 * f)))

Lambda1_num = sp.Float("0.279605891931464")
f_num = sp.Float("0.001161409732093")
ell_num = sp.N(ell_from_f.subs({Lambda1: Lambda1_num, f: f_num}), 18)
scale_num = sp.N(scale_factor.subs({Lambda1: Lambda1_num, f: f_num}), 18)

print("ell electron-point      =", ell_num)
print("normalization scale     =", scale_num)
print("ppm shift of P0         =", sp.N((scale_num - 1) * 10**6, 18))

print("\nFinal reading of Step 44:")
print("  - the Step-43 pure-odd outlet deformation chi_Q = exp(-ell) is not the")
print("    whole story,")
print("  - on the natural point-particle source-map branch it forces the grouped-P2")
print("    normalization defect to be N_Q = exp(+ell),")
print("  - and therefore the moving-throat normalization product rescales as")
print("        P0 = exp(+ell) P0_target.")
