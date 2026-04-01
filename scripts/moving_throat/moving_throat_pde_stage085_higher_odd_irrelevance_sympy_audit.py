#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

omega = sp.symbols('omega', real=True)
Omega, chiQ, tauQ = sp.symbols('Omega chiQ tauQ', positive=True, real=True)

sigma_can = sp.Rational(9, 8) / Omega**5
Y = sp.Rational(3,4) + sp.Rational(1,4) / (1 - omega**2/Omega**2 - sp.I*chiQ*sigma_can*omega**5 - sp.I*tauQ*omega**7)
Yser5 = sp.expand(sp.series(Y, omega, 0, 6).removeO())
Yser8 = sp.expand(sp.series(Y, omega, 0, 8).removeO())

print('series through O(omega^5) =', Yser5)
print('series through O(omega^7) =', Yser8)
print('tauQ coefficient in omega^5 term =', sp.simplify(sp.diff(sp.im(Yser5.coeff(omega,5)), tauQ)))
print('tauQ coefficient in omega^7 term =', sp.simplify(sp.diff(sp.im(Yser8.coeff(omega,7)), tauQ)))
print('check canonical odd coefficient =', sp.simplify(sp.im(Yser5.coeff(omega,5)) - chiQ*sp.Rational(9,32)/Omega**5))
