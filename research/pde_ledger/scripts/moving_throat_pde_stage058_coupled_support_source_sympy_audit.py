#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 41 SymPy audit.

Checks:
1. exact support-drop kernel and its positive derivative identity,
2. normalization of the constructive source family Sigma_Pe,
3. exact antiderivative formulas for Ic and Is,
4. exact uniform-source drop Delta_0,
5. exact sharp-bottom endpoint Delta_inf,
6. fixed-point branch bracket Pe_* in [Xi Delta_0, Xi Delta_inf].
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


banner("STAGE 41 — COUPLED SUPPORT/SOURCE OPERATOR")

x, s = sp.symbols("x s", real=True)
Pe, alpha, eta, Xi = sp.symbols("Pe alpha eta Xi", positive=True, real=True)

W = sp.simplify(alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))
K = sp.simplify(
    (sp.cosh(alpha * x) + (eta / alpha) * sp.sinh(alpha * x) - sp.cosh(alpha * (1 - x))) / W
)
Kprime = sp.simplify(sp.diff(K, x))
print("K_(alpha,eta)(x) =", K)
print("dK/dx =", Kprime)
expect_zero(
    "Kprime identity",
    Kprime
    - (alpha * sp.sinh(alpha * x) + eta * sp.cosh(alpha * x) + alpha * sp.sinh(alpha * (1 - x))) / W,
)

Sigma = sp.simplify(Pe * sp.exp(Pe * x) / (sp.exp(Pe) - 1))
print("Sigma_Pe(x) =", Sigma)
expect_zero("Sigma normalization", sp.integrate(Sigma, (x, 0, 1)) - 1)

# Exact auxiliary integrals for the closed-form support drop.
Fc = sp.exp(Pe * x) * (Pe * sp.cosh(alpha * x) - alpha * sp.sinh(alpha * x)) / (Pe**2 - alpha**2)
Fs = sp.exp(Pe * x) * (Pe * sp.sinh(alpha * x) - alpha * sp.cosh(alpha * x)) / (Pe**2 - alpha**2)
expect_zero("Ic antiderivative check", sp.diff(Fc, x) - sp.exp(Pe * x) * sp.cosh(alpha * x))
expect_zero("Is antiderivative check", sp.diff(Fs, x) - sp.exp(Pe * x) * sp.sinh(alpha * x))

Ic = sp.simplify(Fc.subs(x, 1) - Fc.subs(x, 0))
Is = sp.simplify(Fs.subs(x, 1) - Fs.subs(x, 0))
print("Ic(Pe,alpha) =", Ic)
print("Is(Pe,alpha) =", Is)

Delta = sp.simplify(
    Pe / (sp.exp(Pe) - 1) * ((1 - sp.cosh(alpha)) * Ic + (eta / alpha + sp.sinh(alpha)) * Is) / W
)
print("Delta(Pe;alpha,eta) =", Delta)

Delta0 = sp.simplify(sp.limit(Delta, Pe, 0))
Delta0_expected = sp.simplify(eta * (sp.cosh(alpha) - 1) / (alpha**2 * W))
print("Delta_0 =", Delta0)
expect_zero("Delta0 formula", Delta0 - Delta0_expected)
expect_zero("Delta0 integral identity", Delta0 - sp.integrate(K, (x, 0, 1)))

Delta_inf = sp.simplify(K.subs(x, 1))
Delta_inf_expected = sp.simplify((sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) / W)
print("Delta_inf =", Delta_inf)
expect_zero("Delta_inf direct substitution (sanity)", Delta_inf - Delta_inf_expected)

# Fixed-point branch bracket.
Pe_lo = sp.simplify(Xi * Delta0)
Pe_hi = sp.simplify(Xi * Delta_inf)
print("Pe_lo = Xi Delta_0 =", Pe_lo)
print("Pe_hi = Xi Delta_inf =", Pe_hi)

# Bracket non-emptiness: Delta_inf >= Delta_0 for all alpha, eta > 0.
bracket_gap = sp.simplify(sp.together(Delta_inf - Delta0))
bracket_gap_expected = sp.simplify(
    ((alpha**2 - eta) * (sp.cosh(alpha) - 1) + alpha * eta * sp.sinh(alpha))
    / (alpha**2 * W)
)
expect_zero("bracket gap closed form", bracket_gap - bracket_gap_expected)
for a_val in [sp.Rational(1, 10), sp.Integer(1), sp.Integer(3)]:
    for e_val in [sp.Rational(1, 10), sp.Integer(1), sp.Integer(10)]:
        val = float(sp.N(bracket_gap.subs({alpha: a_val, eta: e_val})))
        if val <= 0:
            raise AssertionError(
                f"bracket gap non-positive at alpha={a_val}, eta={e_val}: {val}"
            )
print("bracket gap positivity sweep = PASS")

# Delta_inf is the sharp-bottom (Pe -> oo) endpoint of Delta(Pe).
Delta_inf_limit = sp.simplify(sp.limit(Delta, Pe, sp.oo))
print("Delta(Pe -> oo) =", Delta_inf_limit)
expect_zero("Delta_inf as Pe -> oo limit", sp.simplify((Delta_inf_limit - Delta_inf_expected).rewrite(sp.exp)))

# Weak-coupling branch law.
Delta_series = sp.series(Delta, Pe, 0, 2).removeO()
print("Delta(Pe) small-Pe series =", Delta_series)
expect_zero("weak-coupling constant term", sp.expand(Delta_series).coeff(Pe, 0) - Delta0)
Pe1_coeff = sp.simplify(sp.expand(Delta_series).coeff(Pe, 1))
print("Delta(Pe) first-order coefficient =", Pe1_coeff)
Pe1_val = sp.N(Pe1_coeff.subs({alpha: sp.Integer(1), eta: sp.Integer(1)}))
if sp.Abs(Pe1_val) < sp.Rational(1, 10**8):
    raise AssertionError(
        f"weak-coupling first-order coefficient vanishes at alpha=eta=1: {Pe1_val}"
    )
print("weak-coupling first-order coefficient nonvanishing at alpha=eta=1: PASS")

print("\nStage 41 audit passed.")
