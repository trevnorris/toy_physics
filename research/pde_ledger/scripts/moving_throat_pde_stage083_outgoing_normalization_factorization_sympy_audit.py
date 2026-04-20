#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

G, c = sp.symbols('G c', positive=True, real=True)
Omega = sp.symbols('Omega', positive=True, real=True)
K0, mhat0, chiQ = sp.symbols('K0 mhat0 chiQ', positive=True, real=True)
omega = sp.symbols('omega', real=True)

sigma_can = sp.Rational(9, 8) / Omega**5
Y = sp.Rational(3,4) + sp.Rational(1,4) / (1 - omega**2/Omega**2 - sp.I*chiQ*sigma_can*omega**5)
Yser = sp.expand(sp.series(Y, omega, 0, 6).removeO())

K2 = sp.simplify(K0 * Yser.coeff(omega, 2))
K4 = sp.simplify(K0 * Yser.coeff(omega, 4))
Gamma5 = sp.simplify(sp.im(Yser.coeff(omega, 5)) * K0)

K0_t = sp.simplify(64*G*Omega**5/(45*c**5))
K2_t = sp.simplify(K0_t/(4*Omega**2))
K4_t = sp.simplify(K0_t/(4*Omega**4))
Gamma5_t = sp.simplify(2*G/(5*c**5))
NQ = sp.simplify(K0 / K0_t)

print('Yhat_Q^ret series =', Yser)
print('K2 =', K2)
print('K4 =', K4)
print('Gamma5 =', Gamma5)
print('NQ =', NQ)
print('K2/K2_target - NQ =', sp.simplify(K2/K2_t - NQ))
print('K4/K4_target - NQ =', sp.simplify(K4/K4_t - NQ))
print('Gamma5/Gamma5_target - chiQ*NQ =', sp.simplify(Gamma5/Gamma5_t - chiQ*NQ))
print('mhat0^2*Gamma5/Gamma5_target - mhat0^2*chiQ*NQ =', sp.simplify(mhat0**2*Gamma5/Gamma5_t - mhat0**2*chiQ*NQ))
print('NQ - 1/(mhat0^2*chiQ) after odd normalization =', sp.simplify(sp.solve(sp.Eq(mhat0**2*chiQ*NQ, 1), NQ)[0] - 1/(mhat0**2*chiQ)))
