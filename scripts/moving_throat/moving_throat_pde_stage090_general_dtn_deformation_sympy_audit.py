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

banner('STAGE 90 — GENERAL ISOTROPIC DTN DEFORMATION ALGEBRA')

Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)
print('Lambda_out(z) =', Lambda_out)

Lambda_def = sp.expand(S*Lambda_out.subs(z, beta*z) + Sigma0 + Sigma2*z**2 + Sigma4*z**4 + I*Sigma5*z**5)
print('Lambda_def(z) =', Lambda_def)

L0 = sp.simplify(Lambda_def.subs(z, 0))
L2 = sp.simplify(sp.expand(Lambda_def).coeff(z, 2))
L4 = sp.simplify(sp.expand(Lambda_def).coeff(z, 4))
L5 = sp.simplify(sp.im(sp.expand(Lambda_def).coeff(z, 5)))
print('L0 =', L0)
print('L2 =', L2)
print('L4 =', L4)
print('L5 =', L5)

Y_direct = sp.series(L0 / Lambda_def, z, 0, 6).removeO()
Y_formula = 1 - (L2/L0)*z**2 + (L2**2/L0**2 - L4/L0)*z**4 - I*(L5/L0)*z**5
expect_zero('normalized expansion direct-formula', sp.expand(Y_direct - Y_formula))
print('Y_def(z) =', sp.expand(Y_formula))

m2 = sp.simplify(-L2 / L0)
m4 = sp.simplify(L2**2 / L0**2 - L4 / L0)
chiQ = sp.simplify((-L5 / L0) / sp.Rational(1, 27))
print('m2 =', m2)
print('m4 =', m4)
print('chi_Q =', chiQ)

sol = sp.solve([sp.Eq(m2, sp.Rational(1, 9)), sp.Eq(m4, sp.Rational(4, 81))], [Sigma2, Sigma4], dict=True)
if len(sol) != 1:
    raise AssertionError('Expected a unique solution for (Sigma2,Sigma4).')
sol = sol[0]
print('Sigma2_evenmatch =', sp.simplify(sol[Sigma2]))
print('Sigma4_evenmatch =', sp.simplify(sol[Sigma4]))

chi_even = sp.simplify(chiQ.subs(sol))
print('chi_Q under canonical-even matching =', sp.factor(chi_even))
expect_zero(
    'chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0)',
    chi_even - 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0)
)

print('\nFinal exact formulas:')
print('  Sigma2 =', sp.simplify(sol[Sigma2]))
print('  Sigma4 =', sp.simplify(sol[Sigma4]))
print('  chi_Q  =', sp.factor(chi_even))
