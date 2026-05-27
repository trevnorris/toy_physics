#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

z, sigma, kappa, gamma = sp.symbols('z sigma kappa gamma', real=True)
I = sp.I

Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)
Lambda_mix = sp.expand(sp.series(Lambda_out - sigma/(1 - kappa*z**2 - I*gamma*z**5), z, 0, 6).removeO())

L0 = sp.simplify(Lambda_mix.coeff(z, 0))
L2 = sp.simplify(Lambda_mix.coeff(z, 2))
L4 = sp.simplify(Lambda_mix.coeff(z, 4))
L5 = sp.simplify(Lambda_mix.coeff(z, 5) / I)

print('Lambda_mix(z) =', Lambda_mix)
print('L0 =', L0)
print('L2 =', L2)
print('L4 =', L4)
print('L5 =', L5)

kappa_match = sp.solve(sp.Eq(-L2/L0, sp.Rational(1, 9)), kappa)[0]
sigma_match = sp.solve(sp.Eq((L2**2/L0**2 - L4/L0).subs(kappa, kappa_match), sp.Rational(4, 81)), sigma)[0]
chi_mix = sp.simplify((-L5/L0) / sp.Rational(1, 27))
chi_mix_lin = sp.expand(sp.series(chi_mix, sigma, 0, 2).removeO())

print('kappa from z^2 matching =', kappa_match)
print('sigma from z^4 matching =', sigma_match)
print('chi_Q^mix =', chi_mix)
print('chi_Q^mix linearized =', chi_mix_lin)

assert sp.simplify(kappa_match + sp.Rational(1, 9)) == 0
assert sp.simplify(sigma_match) == 0
assert sp.simplify(chi_mix - 3*(1 - 9*sigma*gamma)/(3 + sigma)) == 0
assert sp.simplify(chi_mix_lin - (1 - sigma*(sp.Rational(1,3) + 9*gamma))) == 0
print('stage111: PASS')
