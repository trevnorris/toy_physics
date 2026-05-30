#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

# r_F1 closed form from Stage 121.
rF = sp.N(sp.sqrt(sp.Rational(12,1)/sp.pi**2 * (sp.Rational(37,20))**2 - 1), 30)
# Pi_* and S_q(Pi_*) imported from Stage 134.
Pi_star = sp.N('1.50882951349316', 30)
Sq_star = sp.N('0.658075937605429', 30)
# (Self-matched susceptibility closure Sigma_0 = M_s is established at Stage 140, not here;
#  Stage 139 only evaluates the gain pair on the Family-1 branch. See P4-42.)

Rq_nat = sp.N((1 - rF)**2 / (1 + rF**2), 30)
Ms_nat = sp.N(Pi_star / (1 - Rq_nat * Sq_star), 30)
Mq_nat = sp.N(-Rq_nat * Ms_nat, 30)

# g_- is the LOWER compensated branch g_c = rF - (1/2) sqrt(1 + rF^2) (notes stage139 lines
# 91-100). R_q^comp = 1/4 is DEFINITIONAL on this branch (true for any rF) AND branch-blind
# (g_+ gives 1/4 too). The falsifiable content is the g_-^F1 VALUE (discriminates the lower
# branch, below) and the natural-branch R_q^nat != 1/4 (line 41), NOT R_q^comp = 1/4.
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

# (R2 anchor) g_-^F1 VALUE vs the canonical cross-stage literal (owned at 127/142/144/164/169).
# This DISCRIMINATES the lower compensated branch (g_- ~ 0.758) from the upper (g_+ ~ 2.79),
# which R_q = 1/4 cannot. Falsifiable: a sign/branch/rF typo gives ~2.79 and FAILS.
assert abs(g_minus - sp.Float('0.758035078944662826919680890414', 30)) < tol_literal, (g_minus,)

# Natural-branch mouth gains vs notes literals
assert abs(Ms_nat - sp.Float('1.66854252965624', 20)) < tol_literal, (Ms_nat,)
assert abs(Mq_nat - sp.Float('-0.242696939724365', 20)) < tol_literal, (Mq_nat,)

# Compensated-branch mouth gains vs notes literals
assert abs(Ms_comp - sp.Float('1.80594111095636', 20)) < tol_literal, (Ms_comp,)
assert abs(Mq_comp - sp.Float('-0.451485277739090', 20)) < tol_literal, (Mq_comp,)

# R1 (was tautological outlet residual): independently reconstruct S_q(Pi_*) from
# the Stage 134 closed-form mouth-response kernel S(Pi, kappa) at kappa = pi/2, and
# confirm it equals the imported Stage 134 literal Sq_star. This exercises the
# imported value via a route that does NOT reuse M_s = Pi_*/(1 - R_q S_q).
kappa_q = sp.pi/2
S_kernel = (Pi_star * (kappa_q*sp.tanh(kappa_q)
            + Pi_star*(sp.exp(-Pi_star)/sp.cosh(kappa_q) - 1))
            / ((1 - sp.exp(-Pi_star)) * (kappa_q**2 - Pi_star**2)))
Sq_recon = sp.N(S_kernel, 30)
assert abs(Sq_recon - Sq_star) < tol_literal, (Sq_recon, Sq_star)

# Documented identity (NOT independent): with M_s = Pi_*/(1 - R_q S_q) and M_q = -R_q M_s
# the outlet form Pi_* = M_s + M_q S_q is true by construction; kept only as a structural
# sanity print, NOT as a verification of the imported literals (see R1 in codex_review).
print('outlet form residual (nat, structural) =',
      sp.N(Pi_star - (Ms_nat + Mq_nat * Sq_star), 5))
print('outlet form residual (comp, structural) =',
      sp.N(Pi_star - (Ms_comp + Mq_comp * Sq_star), 5))

# (R2 definitional-consistency, NOT the anchor) R_q^comp = 1/4 holds by construction of the
# compensated branch (g_- = rF - sqrt(1+rF^2)/2 ⇒ (g_- - rF)^2/(1+rF^2) = 1/4 for any rF, both
# branches). The falsifiable checks are the g_-^F1 value (line ~42) and R_q^nat != 1/4 (line 41).
assert abs(Rq_comp - sp.Rational(1, 4)) < tol_algebraic, (Rq_comp,)

print('all assertions passed')
