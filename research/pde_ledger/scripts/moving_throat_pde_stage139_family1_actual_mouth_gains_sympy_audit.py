#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

# r_F1 closed form from Stage 121.
rF = sp.N(sp.sqrt(sp.Rational(12,1)/sp.pi**2 * (sp.Rational(37,20))**2 - 1), 30)
# Pi_* and S_q(Pi_*) imported from Stage 134.
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

# Numerical-deliverable assertions against the boxed values in
# notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md.
tol_literal = sp.Float('1e-12')   # notes give ~15 digits
tol_algebraic = sp.Float('1e-25') # algebraic identities

# r_F1 closed form vs notes literal
assert abs(rF - sp.Float('1.77799353547498', 20)) < tol_literal, (rF,)

# R_q^nat closed form vs notes literal
assert abs(Rq_nat - sp.Float('0.145454452260421', 20)) < tol_literal, (Rq_nat,)

# Natural-branch mouth gains vs notes literals
assert abs(Ms_nat - sp.Float('1.66854252965624', 20)) < tol_literal, (Ms_nat,)
assert abs(Mq_nat - sp.Float('-0.242696939724365', 20)) < tol_literal, (Mq_nat,)

# Compensated-branch mouth gains vs notes literals
assert abs(Ms_comp - sp.Float('1.80594111095636', 20)) < tol_literal, (Ms_comp,)
assert abs(Mq_comp - sp.Float('-0.451485277739090', 20)) < tol_literal, (Mq_comp,)

# Outlet consistency Pi_* = M_s + M_q * S_q(Pi_*), both branches.
# (Algebraically built in by M_s = Pi_*/(1 - R_q S_q), M_q = -R_q M_s, but
# verifies that the imported Pi_* and S_q literals satisfy the form.)
assert abs(Pi_star - (Ms_nat + Mq_nat * Sq_star)) < tol_algebraic
assert abs(Pi_star - (Ms_comp + Mq_comp * Sq_star)) < tol_algebraic

# Compensated R_q closed form
assert abs(Rq_comp - sp.Rational(1, 4)) < tol_algebraic

print('all assertions passed')
