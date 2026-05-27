#!/usr/bin/env python3
"""Stage 109 SymPy audit — linearized branch-selection law.

Paper card's secondary `Checks` items (Robin and standalone mixed-pole limits;
even-coefficient preservation under the compensated branch) are exercised at
downstream stages: stage 110 (Robin outlet model), stage 111 (mixed side-channel
pole), stage 112 (hybrid Robin+mixed compensation, which verifies the
compensated branch's even-coefficient preservation alongside the odd-norm
match). This stage establishes the linearized framework that those downstream
verifications consume.
"""
from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


eps = sp.symbols('eps', real=True)
s, b, a0, a5 = sp.symbols('s b a0 a5', real=True)

banner('STAGE 109 — LINEARIZED BRANCH-SELECTION LAW')

S = 1 + eps*s
beta = 1 + eps*b
Sigma0 = eps*a0
Sigma5 = eps*a5
chi = sp.simplify(3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0))
chi_series = sp.series(chi, eps, 0, 2).removeO()
print('chi_Q series =', sp.expand(chi_series))

expected = 1 + eps*(5*b + a0/sp.Integer(3) + 9*a5)
expect_zero('linearized chi law', chi_series - expected)

coeff = sp.expand((chi_series - 1)/eps)
print('first-order defect coefficient =', coeff)
expect_zero('overall scale cancels', sp.diff(coeff, s))

# Linearized preservation condition.
a5_sol = sp.solve(sp.Eq(coeff, 0), a5)[0]
print('a5 preservation condition =', sp.simplify(a5_sol))
expected_a5_sol = -sp.Rational(5, 9)*b - sp.Rational(1, 27)*a0
expect_zero('a5 preservation closed-form', sp.simplify(a5_sol - expected_a5_sol))
expect_zero('preservation substitution', coeff.subs(a5, a5_sol))
