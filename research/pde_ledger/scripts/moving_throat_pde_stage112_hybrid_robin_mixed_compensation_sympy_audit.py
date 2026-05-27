#!/usr/bin/env python3
"""
Stage 112 SymPy audit.

Provenance notes
----------------
- `Lambda_out` is the exact canonical outgoing DtN branch from Stages 104/105.
- The Robin denominator and the two-branch hybrid solve are the same Stage 110
  / Stage 111 compensation seam; this script only checks the exact `chi_Q`
  consequences of those carried solves.
- Branches A/B are the two Stage 111 hybrid Robin-mixed compensation solves;
  this script keeps that two-branch naming and only checks the exact `chi_Q`
  consequences of each solve.
"""

from __future__ import annotations
import sympy as sp

z, rho, sigma, kappa, gamma = sp.symbols('z rho sigma kappa gamma', real=True)
I = sp.I

# Stage 104 canonical outgoing DtN branch carried into the hybrid solve.
Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)
Lambda_hyb = sp.expand(sp.series(Lambda_out + rho - sigma/(1 - kappa*z**2 - I*gamma*z**5), z, 0, 6).removeO())

L0 = sp.simplify(Lambda_hyb.coeff(z, 0))
L2 = sp.simplify(Lambda_hyb.coeff(z, 2))
L4 = sp.simplify(Lambda_hyb.coeff(z, 4))
L5 = sp.simplify(Lambda_hyb.coeff(z, 5) / I)

sols = sp.solve([
    sp.Eq(-L2/L0, sp.Rational(1, 9)),
    sp.Eq(L2**2/L0**2 - L4/L0, sp.Rational(4, 81))
], [rho, kappa], dict=True)

print('Lambda_hyb(z) =', Lambda_hyb)
print('canonical-even solutions =', sols)
assert sols == [{kappa: 0, rho: sigma}, {kappa: sp.Rational(1, 3), rho: 4*sigma}]

# Branch A
chi_A = sp.simplify(((-L5/L0) / sp.Rational(1, 27)).subs(sols[0]))
print('chi_Q branch A =', chi_A)
assert sp.simplify(chi_A - (1 - 9*sigma*gamma)) == 0

# Branch B
chi_B = sp.simplify(((-L5/L0) / sp.Rational(1, 27)).subs(sols[1]))
print('chi_Q branch B =', chi_B)
assert sp.simplify(chi_B - (1 - 9*sigma*gamma)/(1 - sigma)) == 0
assert sp.simplify(chi_B.subs(gamma, sp.Rational(1, 9)) - 1) == 0

scaled_identity = sp.expand((Lambda_hyb.subs(sols[1]) - (1 - sigma)*Lambda_out).subs(gamma, sp.Rational(1, 9)))
print('scaled identity on branch B =', scaled_identity)
assert sp.simplify(scaled_identity) == 0
print('stage112: PASS')
