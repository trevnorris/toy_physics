#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

z, rho = sp.symbols('z rho', real=True)
I = sp.I

Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)
Lambda_R = Lambda_out + rho
Y_R = sp.simplify((-3 + rho) / Lambda_R)
Y_R_series = sp.expand(sp.series(Y_R, z, 0, 6).removeO())

c2 = sp.simplify(Y_R_series.coeff(z, 2))
c4 = sp.simplify(Y_R_series.coeff(z, 4))
c5 = sp.simplify(Y_R_series.coeff(z, 5) / I)
chi_R = sp.simplify(c5 / sp.Rational(1, 27))
chi_R_lin = sp.expand(sp.series(chi_R, rho, 0, 3).removeO())

print('Y_R(z) =', Y_R_series)
print('c2 =', c2)
print('c4 =', c4)
print('c5 =', c5)
print('chi_Q^R =', chi_R)
print('chi_Q^R linearized =', chi_R_lin)

assert sp.simplify(c2 - 1/(9 - 3*rho)) == 0
assert sp.simplify(c4 - (4 - rho)/(9*(3 - rho)**2)) == 0
assert sp.simplify(c5 - 1/(27 - 9*rho)) == 0
assert sp.simplify(chi_R - 3/(3 - rho)) == 0
assert sp.simplify(chi_R_lin - (1 + rho/3 + rho**2/9)) == 0
print('stage110: PASS')
