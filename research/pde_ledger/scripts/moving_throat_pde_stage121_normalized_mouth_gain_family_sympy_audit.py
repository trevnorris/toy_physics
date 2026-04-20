#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

L, Theta = sp.symbols('L Theta', positive=True, real=True)
Ks, Kq, lam, gs, gq = sp.symbols('K_s K_q lam g_s g_q', positive=True, real=True)
r, gc = sp.symbols('r g_c', real=True)
Sigma0 = sp.symbols('Sigma_0', positive=True, real=True)

Ms = sp.simplify(L*gs**2/(Ks*Theta))
Mq = sp.simplify(-L*(Ks*gq - lam*gs)**2/(Ks*(Ks*Kq + lam**2)*Theta))

subs = {
    lam: r*sp.sqrt(Ks*Kq),
    gq: gc*gs*sp.sqrt(Kq/Ks),
}
Ms_norm = sp.simplify(Ms.subs(subs))
Mq_norm = sp.simplify(Mq.subs(subs))
Rq_raw = sp.simplify(-Mq_norm / Ms_norm)
Rq = sp.simplify(Rq_raw.subs(Ms_norm, Sigma0))

print('M_s normalized =', Sigma0)
print('M_q normalized =', sp.simplify(Mq_norm.subs(Ms_norm, Sigma0)))
print('R_q =', Rq)

# Compensation family check
sol_plus = sp.simplify(Rq.subs(gc, r + sp.sqrt(1 + r**2)/2))
sol_minus = sp.simplify(Rq.subs(gc, r - sp.sqrt(1 + r**2)/2))
print('R_q on + branch =', sol_plus)
print('R_q on - branch =', sol_minus)
assert sp.simplify(sol_plus - sp.Rational(1,4)) == 0
assert sp.simplify(sol_minus - sp.Rational(1,4)) == 0

print('\nNormalized mouth-gain family verified.')
