#!/usr/bin/env python3
"""
g2_step43_pure_odd_dtn_representative.py

Step 43 of the g-2 rebuild.

What this script does
---------------------
1. Takes the exact isotropic DtN deformation surface from Step 31,
       chi_Q = 3 (S beta^5 + 9 Sigma5) / (3 S - Sigma0),
   and the Step-42 anomaly target chi_Q = exp(-ell) = 1/(1+x).
2. Imposes the adiabatic anomaly-only closure in the normalized DtN language:
   the even branch is frozen (beta = 1, Sigma0 = 0), so only the odd slot is
   allowed to move.
3. Solves the exact finite-f odd-slot law
       Sigma5 = S (chi_Q - 1) / 9 = - S x / [9 (1 + x)].
4. Shows that this pure-odd normalized representative is independent of the
   microscopic loading-share drift sigma_W: once the hybrid/core outlet is
   reduced to the normalized DtN branch, the whole anomaly is carried by one odd
   coefficient.
5. Expands the result to tangent level and recovers
       b = 0, a0 = 0, a5 = -Lambda1/9
   in canonical normalized gauge S = 1.

Interpretation
--------------
The adiabatic anomaly-only branch admits a very clean normalized DtN description:
all of the observable branch motion can be represented by a single pure odd
isotropic outlet coefficient, with the even branch frozen exactly.
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


banner("STEP 43A — EXACT ISOTROPIC DTN SURFACE ON THE ADIABATIC ANOMALY BRANCH")

ell, x, S, beta, Sigma0, Sigma5 = sp.symbols("ell x S beta Sigma0 Sigma5", real=True)
chi_target = sp.exp(-ell)
chi_target_x = 1 / (1 + x)
chi_surface = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))

print("chi_target(ell) =", chi_target)
print("chi_target(x)   =", chi_target_x)
print("chi_surface     =", chi_surface)
expect_zero("x = exp(ell)-1 rewrite", sp.simplify(chi_target - chi_target_x.subs(x, sp.exp(ell) - 1)))

print("\nReading:")
print("  This is the exact finite-f isotropic DtN branch-selection surface.")
print("  The adiabatic anomaly-only branch now asks how that surface is realized")
print("  once the parent compensation sheet and even D/N geometry are frozen.")


banner("STEP 43B — PURE ODD REPRESENTATIVE OF THE ADIABATIC ANOMALY")

# Fixed-even slice in normalized DtN language:
#   beta = 1   (no argument deformation)
#   Sigma0 = 0 (no static additive isotropic defect)
chi_pure_odd = sp.simplify(chi_surface.subs({beta: 1, Sigma0: 0}))
Sigma5_exact = sp.solve(sp.Eq(chi_pure_odd, chi_target), Sigma5)[0]
Sigma5_x = sp.simplify(Sigma5_exact.subs(sp.exp(ell), 1 + x))

print("chi_pure_odd =", chi_pure_odd)
print("Sigma5(ell)  =", Sigma5_exact)
print("Sigma5(x)    =", Sigma5_x)

expect_zero("pure-odd exact law", sp.simplify(chi_pure_odd.subs(Sigma5, Sigma5_exact) - chi_target))
expect_zero("pure-odd x-law", sp.simplify(Sigma5_x - (-S * x / (9 * (1 + x)))))

print("\nInterpretation:")
print("  Once the even branch is frozen, the exact electron anomaly is represented")
print("  by one pure odd isotropic outlet coefficient Sigma5 and nothing else.")


banner("STEP 43C — THE PURE ODD REPRESENTATIVE IS INDEPENDENT OF THE MICROSCOPIC LOADING SHARE")

sigma = sp.symbols("sigma", positive=True, real=True)
gamma = sp.symbols("gamma", real=True)
chi_hyb = sp.simplify((1 - 9 * sigma * gamma) / (1 - sigma))

# Solve the compensated hybrid/core outlet for gamma at the same chi target.
gamma_from_target = sp.solve(sp.Eq(chi_hyb, chi_target), gamma)[0]
Sigma5_from_hybrid = sp.simplify(S * (chi_hyb - 1) / 9).subs(gamma, gamma_from_target)

print("chi_hyb =", chi_hyb)
print("gamma_from_target =", gamma_from_target)
print("Sigma5 from hybrid normalization =", Sigma5_from_hybrid)

expect_zero(
    "hybrid -> pure odd representative",
    sp.simplify(Sigma5_from_hybrid - Sigma5_exact),
)

print("\nReading:")
print("  The microscopic loading share sigma_W can drift on the compensated hybrid")
print("  outlet, but after reduction to the normalized DtN branch that dependence")
print("  disappears. The observable anomaly is carried entirely by the pure odd")
print("  coefficient Sigma5.")


banner("STEP 43D — TANGENT LAW IN CANONICAL NORMALIZED GAUGE")

eps, b, a0, a5, Lambda1, f = sp.symbols("eps b a0 a5 Lambda1 f", real=True)
chi_eps = chi_surface.subs({S: 1, beta: 1 + eps * b, Sigma0: eps * a0, Sigma5: eps * a5})
chi_eps_series = sp.expand(sp.series(chi_eps, eps, 0, 2).removeO())
chi_target_eps = sp.expand(sp.series(1 / (1 + eps * Lambda1), eps, 0, 2).removeO())

print("chi_eps series    =", chi_eps_series)
print("chi_target series =", chi_target_eps)

# General tangent law first.
constraint = sp.expand(sp.series(chi_eps - chi_target_eps, eps, 0, 2).removeO().coeff(eps, 1))
print("general tangent constraint =", constraint)
expect_zero(
    "Step-30 tangent law",
    sp.simplify(constraint - (5 * b + a0 / 3 + 9 * a5 + Lambda1)),
)

# Adiabatic anomaly-only even-frozen slice: b = 0, a0 = 0.
a5_sol = sp.solve(sp.Eq(constraint.subs({b: 0, a0: 0}), 0), a5)[0]
print("a5 (even-frozen slice) =", a5_sol)
expect_zero("pure odd tangent value", sp.simplify(a5_sol + Lambda1 / 9))

# Canonical finite-f gauge S = 1.
Sigma5_canonical = sp.simplify(Sigma5_x.subs(S, 1).subs(x, Lambda1 * f))
print("Sigma5 canonical finite-f =", Sigma5_canonical)

Lambda1_num = sp.Float("0.279605891931464")
f_num = sp.Float("0.001161409732093")
Sigma5_num = sp.N(Sigma5_canonical.subs({Lambda1: Lambda1_num, f: f_num}), 18)
a5_num = sp.N(a5_sol.subs(Lambda1, Lambda1_num), 18)
print("Sigma5 canonical electron-point =", Sigma5_num)
print("a5 canonical tangent           =", a5_num)

print("\nFinal reading of Step 43:")
print("  - the adiabatic anomaly-only branch admits an exact pure-odd normalized")
print("    DtN representative,")
print("  - in canonical normalized gauge this is Sigma5 = -x/[9(1+x)] with")
print("    x = Lambda1 f,")
print("  - and at tangent level the branch simply picks")
print("        b = 0, a0 = 0, a5 = -Lambda1/9.")
