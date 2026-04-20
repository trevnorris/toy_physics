#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

rF = sp.N(sp.sqrt(sp.Rational(12,1)/sp.pi**2 * (sp.Rational(37,20))**2 - 1), 30)
Pi_star = sp.N('1.50882951349316', 30)
Sq_star = sp.N('0.658075937605429', 30)

Rq_nat = sp.N((1 - rF)**2 / (1 + rF**2), 30)
Ms_nat = sp.N(Pi_star / (1 - Rq_nat * Sq_star), 30)
Mq_nat = sp.N(-Rq_nat * Ms_nat, 30)

g_minus = sp.N(rF - sp.sqrt(1 + rF**2)/2, 30)
Rq_comp = sp.N((g_minus - rF)**2 / (1 + rF**2), 30)
Ms_comp = sp.N(Pi_star / (1 - Rq_comp * Sq_star), 30)
Mq_comp = sp.N(-Rq_comp * Ms_comp, 30)

print('r_F1 =', rF)
print('S_q(Pi_*) =', Sq_star)
print('R_q^nat =', Rq_nat)
print('M_s^nat,* =', Ms_nat)
print('M_q^nat,* =', Mq_nat)
print('g_-^F1 =', g_minus)
print('R_q^comp =', Rq_comp)
print('M_s^comp,* =', Ms_comp)
print('M_q^comp,* =', Mq_comp)
print('shell gain fractional shift =', sp.N(Ms_comp/Ms_nat - 1, 20))
print('mixed gain magnitude ratio =', sp.N(abs(Mq_comp)/abs(Mq_nat), 20))

assert abs(Rq_comp - sp.Rational(1,4)) < sp.Float('1e-25')
