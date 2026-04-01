#!/usr/bin/env python3
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


z = sp.symbols('z', real=True)
S, beta = sp.symbols('S beta', nonzero=True, real=True)
Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols('Sigma0 Sigma2 Sigma4 Sigma5', real=True)
I = sp.I

banner('STAGE 91 — ROBUSTNESS CLASSES FOR chi_Q')

Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)

# Class A: pure scale
Y_scale = sp.series((-3*S)/(S*Lambda_out), z, 0, 6).removeO()
Y_can = sp.series((-3)/Lambda_out, z, 0, 6).removeO()
expect_zero('pure scale invariance', Y_scale - Y_can)

# Class B: pure scale+argument
Y_arg = sp.series((-3*S)/(S*Lambda_out.subs(z, beta*z)), z, 0, 6).removeO()
print('Y_scale+arg(z) =', sp.expand(Y_arg))
m2_arg = sp.expand(Y_arg).coeff(z, 2)
m4_arg = sp.expand(Y_arg).coeff(z, 4)
chi_arg = sp.simplify(sp.im(sp.expand(Y_arg).coeff(z, 5)) / sp.Rational(1, 27))
print('m2_arg =', m2_arg)
print('m4_arg =', m4_arg)
print('chi_arg =', chi_arg)

sol_beta = sp.solve([sp.Eq(m2_arg, sp.Rational(1, 9)), sp.Eq(m4_arg, sp.Rational(4, 81))], [beta], dict=True)
print('solutions preserving canonical even fingerprint =', sol_beta)
if {sol[beta] for sol in sol_beta} != {sp.Integer(-1), sp.Integer(1)}:
    raise AssertionError('Expected the two algebraic roots beta=±1.')
expect_zero('chi_arg(beta=1) - 1', chi_arg.subs(beta, 1) - 1)

# Class C: additive core, beta=1
Lambda_add = sp.expand(S*Lambda_out + Sigma0 + Sigma2*z**2 + Sigma4*z**4 + I*Sigma5*z**5)
L0 = sp.simplify(Lambda_add.subs(z, 0))
L2 = sp.simplify(sp.expand(Lambda_add).coeff(z, 2))
L4 = sp.simplify(sp.expand(Lambda_add).coeff(z, 4))
L5 = sp.simplify(sp.im(sp.expand(Lambda_add).coeff(z, 5)))
m2 = sp.simplify(-L2/L0)
m4 = sp.simplify(L2**2/L0**2 - L4/L0)
sol = sp.solve([sp.Eq(m2, sp.Rational(1, 9)), sp.Eq(m4, sp.Rational(4, 81))], [Sigma2, Sigma4], dict=True)
if len(sol) != 1:
    raise AssertionError('Expected unique additive even-match solution.')
sol = sol[0]
print('Sigma2(beta=1) =', sp.simplify(sol[Sigma2]))
print('Sigma4(beta=1) =', sp.simplify(sol[Sigma4]))
chi_add = sp.simplify((-(L5/L0)).subs(sol) / sp.Rational(1, 27))
print('chi_add =', sp.factor(chi_add))
expect_zero('chi_add - 3(S+9Sigma5)/(3S-Sigma0)', chi_add - 3*(S + 9*Sigma5)/(3*S - Sigma0))

chi_pres = sp.solve(sp.Eq(chi_add, 1), Sigma5)[0]
print('Sigma5 preservation locus =', sp.simplify(chi_pres))
expect_zero('preservation locus check', chi_add.subs(Sigma5, chi_pres) - 1)
