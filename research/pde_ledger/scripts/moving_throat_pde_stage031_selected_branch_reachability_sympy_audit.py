#!/usr/bin/env python3
"""
Stage 14 SymPy audit.

Checks:
1. Exact derivative of the selected overlap s_-(alpha).
2. Exact derivative formula for the selected prefactor P0_-(alpha).
3. Initial-value formulas at alpha = 0.
4. Exact softening threshold and determinant factorization.
5. Stable-side monotonicity structure in exact symbolic form.
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


banner("PART I — EXACT SELECTED OVERLAP DERIVATIVE")

A, DK = sp.symbols("A DK", positive=True, real=True)
alpha = sp.symbols("alpha", nonnegative=True, real=True)
x0, x1 = sp.symbols("x0 x1", positive=True, real=True)
beta0 = sp.symbols("beta0", positive=True, real=True)

sigma = x0 + x1
delta_kappa = x0 - x1
KappaProd = x0 * x1
R = sp.sqrt((DK + alpha * delta_kappa) ** 2 + 4 * alpha**2 * KappaProd)

lam_minus = sp.simplify((2 * A + DK - alpha * sigma - R) / 2)
lam_plus = sp.simplify((2 * A + DK - alpha * sigma + R) / 2)

s_minus = sp.simplify(-sp.diff(lam_minus, alpha))
ds_exact = sp.simplify(sp.diff(s_minus, alpha))
ds_expected = sp.simplify(2 * DK**2 * KappaProd / R**3)
expect_zero("ds_-/dalpha exact formula", ds_exact - ds_expected)
print("ds_-/dalpha =")
sp.pprint(ds_expected)

banner("PART II — EXACT PREFATOR MONOTONICITY IDENTITY")

P0_sel = sp.simplify(beta0 * s_minus / lam_minus)

dP_direct = sp.diff(P0_sel, alpha)
dP_physical = beta0 * (ds_expected * lam_minus + s_minus**2) / lam_minus**2
expect_zero("dP0_-/dalpha direct identity", sp.simplify(dP_direct - dP_physical))

print("dP0_-/dalpha =")
sp.pprint(sp.simplify(dP_physical))

banner("PART III — INITIAL VALUES")

expect_zero("lambda_-(0)", sp.simplify(lam_minus.subs(alpha, 0) - A))
expect_zero("s_-(0)", sp.simplify(s_minus.subs(alpha, 0) - x0))
expect_zero("P0_-(0)", sp.simplify(P0_sel.subs(alpha, 0) - beta0 * x0 / A))

banner("PART IV — EXACT SOFTENING THRESHOLD")

T0 = sp.simplify((A + DK) * x0 + A * x1)
alpha_crit = sp.simplify(A * (A + DK) / T0)
expect_zero("det factorization", sp.expand(lam_minus * lam_plus - (A * (A + DK) - alpha * T0)))
root_crit = A**2 * x1 + (A + DK) ** 2 * x0
lam_minus_crit = sp.together(lam_minus.subs(alpha, alpha_crit)).replace(
    lambda expr: isinstance(expr, sp.Pow) and expr.exp == sp.Rational(1, 2) and sp.simplify(expr.base - root_crit**2) == 0,
    lambda expr: root_crit,
)
expect_zero("lam_-(alpha_crit)", sp.simplify(lam_minus_crit))

lambda_plus_crit = sp.simplify(lam_plus.subs(alpha, alpha_crit))
radcrit_derived = sp.expand(T0**2 * (R**2).subs(alpha, alpha_crit))
expect_zero(
    "threshold radical square identity",
    sp.expand(radcrit_derived - (A**2 * x1 + (A + DK) ** 2 * x0) ** 2),
)
print("alpha_crit =")
sp.pprint(alpha_crit)
print("lambda_+^(crit) =")
sp.pprint(lambda_plus_crit)

banner("PART V — STABLE-SIDE DIVERGENCE STRUCTURE")

# Exact identity showing linear vanishing of lambda_- through the determinant factor.
expect_zero(
    "lambda_- * lambda_+ - T0*(alpha_crit-alpha)",
    sp.expand(lam_minus * lam_plus - T0 * (alpha_crit - alpha)),
)

# The prefactor rewritten to expose the stable-side divergence.
P0_factored = sp.simplify(beta0 * s_minus * lam_plus / (T0 * (alpha_crit - alpha)))
expect_zero("P0_- factorization", sp.simplify(P0_sel - P0_factored))
print("P0_- =")
sp.pprint(P0_factored)

banner("STAGE 14 AUDIT COMPLETE")
print("Verified:")
print("  • exact overlap derivative ds_-/dalpha = 2 DK^2 kappa0^2 kappa1^2 / R^3")
print("  • exact monotonicity identity for dP0_-/dalpha")
print("  • exact initial-value formulas at alpha = 0")
print("  • exact softening threshold alpha_crit and determinant factorization")
print("  • exact stable-side factorization P0_- = beta0 s_- lambda_+ / [T0 (alpha_crit-alpha)]")
