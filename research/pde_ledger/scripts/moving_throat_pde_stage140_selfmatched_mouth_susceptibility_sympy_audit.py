#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

L, a, ell, rho, Tm, hbar, cs, m, pi = sp.symbols('L a ell rho T_m hbar c_s m pi', positive=True, real=True)
Js = 4*sp.pi*a**2*ell/3
Hw = m*cs**2/rho
Theta = Hw*Js
Ks = 3*sp.pi*a**2*hbar**2/(5*m*rho*ell)
gs = Tm*Js
Sigma0 = sp.simplify(L*gs**2/(Ks*Theta))
print('Sigma_0 =', Sigma0)
That = sp.symbols('That', positive=True, real=True)
Sigma0_hat = sp.simplify(Sigma0.subs(Tm, That*hbar*cs/(rho*ell*sp.sqrt(L))))
print('Sigma_0 in terms of That =', Sigma0_hat)
assert sp.simplify(Sigma0_hat - sp.Rational(20,9)*That**2) == 0

Ms_nat = sp.N('1.6685425296562397', 30)
Ms_comp = sp.N('1.80594111095636', 30)
That_nat = sp.N(sp.sqrt(9*Ms_nat/20), 30)
That_comp = sp.N(sp.sqrt(9*Ms_comp/20), 30)
print('That_nat =', That_nat)
print('That_comp =', That_comp)
print('fractional traction enhancement =', sp.N(That_comp/That_nat - 1, 20))
assert sp.Abs(That_nat - sp.Float('0.866512630228382', 30)) < sp.Float('1e-12', 30)
assert sp.Abs(That_comp - sp.Float('0.901484054174206', 30)) < sp.Float('1e-12', 30)
assert sp.Abs(sp.N(That_comp/That_nat - 1, 30) - sp.Float('0.0403588161624', 30)) < sp.Float('1e-12', 30)
print('numeric fixed-point checks PASS')
