#!/usr/bin/env python3
"""
g2_step52_final_notuning_verdict.py

Step 52 of the g-2 rebuild.

What this script does
---------------------
1. Takes the exact isotropic DtN deformation law from the moving-throat notes,
       chi_Q = 3 (S beta^5 + 9 Sigma5) / (3 S - Sigma0),
   with the canonical-even constraints imposed.
2. Proves the robustness classes:
       chi_Q = 1
   for pure scale deformations and, once the even fingerprint is fixed,
   for pure scale+argument deformations as well.
3. Inserts the adiabatic electron target
       chi_e = exp(-ell) = 1 / (1 + Lambda1 f),
   and solves the full exact deformation law on the adiabatic even-frozen slice
       beta = 1, Sigma0 = 0.
4. Recovers the pure odd representative
       Sigma5 = - S Lambda1 f / [9 (1 + Lambda1 f)],
   and the tangent law
       a5 = -Lambda1/9.
5. States the closing verdict: the current stack naturally derives the canonical
   no-tuning branch, and naturally compresses the electron sliver to one pure odd
   scalar, but it does not yet derive the electron-point magnitude without an
   additional branch-selection law.

Interpretation
--------------
This step is the honest finish line for the current reduced derivation chain.
The background canonical branch is not tuned, but the observed electron anomaly is
still one tiny branch-selection datum rather than a number already forced by the
frozen files alone.
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


banner("STEP 52A — EXACT ISOTROPIC DTN DEFORMATION ALGEBRA")

S, beta, Sigma0, Sigma5 = sp.symbols("S beta Sigma0 Sigma5", positive=True, real=True)
chi_Q = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))

print("chi_Q deformation law =", chi_Q)

# Robustness classes.
chi_scale = sp.simplify(chi_Q.subs({beta: 1, Sigma0: 0, Sigma5: 0}))
chi_scale_arg = sp.simplify(chi_Q.subs({Sigma0: 0, Sigma5: S * (1 - beta**5) / 9}))

print("chi_Q (pure scale) =", chi_scale)
print("chi_Q (preserving submanifold for beta deformation) =", chi_scale_arg)
expect_zero("pure scale robustness", sp.simplify(chi_scale - 1))
expect_zero("preserving submanifold robustness", sp.simplify(chi_scale_arg - 1))

print("\nReading:")
print("  The canonical outgoing branch is robust against pure overall mouth")
print("  normalization, and against pure argument rescaling once the canonical")
print("  even fingerprint is preserved. A nontrivial chi_Q therefore requires a")
print("  genuine throat-core deformation.")


banner("STEP 52B — INSERT THE ADIABATIC ELECTRON TARGET")

ell, Lambda1, f = sp.symbols("ell Lambda1 f", positive=True, real=True)
chi_e = sp.exp(-ell)
chi_e_f = sp.simplify(1 / (1 + Lambda1 * f))

print("chi_e(ell) =", chi_e)
print("chi_e(f)   =", chi_e_f)
expect_zero("ell to f rewrite", sp.simplify(chi_e - chi_e_f.subs(Lambda1 * f, sp.exp(ell) - 1)))

Delta_Q = sp.simplify(chi_e_f - 1)
print("Delta_Q exact =", Delta_Q)

Lambda1_num = sp.Float("0.279605891931464")
f_num = sp.Float("0.001161409732093")
Delta_num = sp.N(Delta_Q.subs({Lambda1: Lambda1_num, f: f_num}), 18)
print("Delta_Q electron-point =", Delta_num)
print("Delta_Q in ppm         =", sp.N(Delta_num * 10**6, 18))


banner("STEP 52C — ADIABATIC EVEN-FROZEN SLICE: THE ANOMALY IS ONE PURE ODD SCALAR")

chi_even_frozen = sp.simplify(chi_Q.subs({beta: 1, Sigma0: 0}))
Sigma5_exact = sp.simplify(S * (chi_e_f - 1) / 9)

print("chi_even_frozen =", chi_even_frozen)
print("Sigma5 exact    =", Sigma5_exact)
expect_zero(
    "pure odd finite-f law",
    sp.simplify(Sigma5_exact + S * Lambda1 * f / (9 * (1 + Lambda1 * f))),
)
expect_zero(
    "chi_even_frozen consistency",
    sp.simplify(chi_even_frozen.subs(Sigma5, Sigma5_exact) - chi_e_f),
)

Sigma5_num = sp.N(Sigma5_exact.subs({S: 1, Lambda1: Lambda1_num, f: f_num}), 18)
print("Sigma5 electron-point (S=1) =", Sigma5_num)

print("\nReading:")
print("  Once the adiabatic wall track freezes the even branch, the whole electron")
print("  anomaly is carried by one pure odd isotropic throat-core outlet scalar.")


banner("STEP 52D — LINEARIZED BRANCH-SELECTION LAW")

eps, b, a0, a5 = sp.symbols("eps b a0 a5", real=True)
chi_lin = sp.expand(sp.series(chi_Q.subs({S: 1, beta: 1 + eps*b, Sigma0: eps*a0, Sigma5: eps*a5}), eps, 0, 2).removeO())
chi_target_lin = sp.expand(sp.series(1 / (1 + eps * Lambda1), eps, 0, 2).removeO())
constraint = sp.expand((chi_lin - chi_target_lin).coeff(eps, 1))

print("linear constraint =", constraint)
expect_zero(
    "general tangent deformation law",
    sp.simplify(constraint - (5 * b + a0 / 3 + 9 * a5 + Lambda1)),
)

a5_sol = sp.solve(sp.Eq(constraint.subs({b: 0, a0: 0}), 0), a5)[0]
print("a5 on the even-frozen adiabatic slice =", a5_sol)
expect_zero("pure odd tangent value", sp.simplify(a5_sol + Lambda1 / 9))

print("\nReading:")
print("  At tangent level the adiabatic anomaly branch is the pure odd selection")
print("      a5 = -Lambda1/9,")
print("  with no need to move the even DtN slots at all.")


banner("STEP 52E — FINAL VERDICT ON 'NO TUNING'")

print("1. The canonical compact outgoing branch itself is a genuine no-tuning result:")
print("      chi_Q = 1, N_Q = 1.")
print("2. The adiabatic electron anomaly does not reopen a broad fit space.")
print("   It collapses to one tiny pure odd scalar Sigma5 (or a5 at tangent level).")
print("3. But the current frozen stack still does NOT derive the numerical electron")
print("   value of that scalar from first principles alone.")
print("   It needs an additional branch-selection law for the actual isotropic DtN")
print("   deformation of the electron branch.")
print("4. So the strongest honest status is:")
print("      - the background canonical branch is naturally derived with no tuning;")
print("      - the observed electron sliver is one constrained branch datum, not a")
print("        multi-parameter tune;")
print("      - its actual numerical value is still open.")
