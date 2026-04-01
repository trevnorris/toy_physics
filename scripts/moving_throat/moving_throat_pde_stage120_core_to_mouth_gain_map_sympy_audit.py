#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

L, Theta = sp.symbols('L Theta', positive=True, real=True)
Ks, Kq, lam, gs, gq = sp.symbols('K_s K_q lam g_s g_q', positive=True, real=True)
rc = sp.symbols('r_c', positive=True, real=True)

rho_c = sp.simplify(gs**2 / Ks)
sigma_c = sp.simplify((Ks*gq - lam*gs)**2 / (Ks*(Ks*Kq + lam**2)))

Ms = sp.simplify(L * rho_c / Theta)
Mq = sp.simplify(-L * sigma_c / Theta)

print('rho_c =', rho_c)
print('sigma_c =', sigma_c)
print('M_s =', Ms)
print('M_q =', Mq)

# Check equivalence with r_c notation.
expr_rc = sp.simplify((Ks*gq - lam*gs)**2 / (Ks**2*Kq*(1 + lam**2/(Ks*Kq))))
print('sigma_c (r_c form) =', expr_rc)
assert sp.simplify(sigma_c - expr_rc) == 0

print('\nFinal explicit gain map verified.')
