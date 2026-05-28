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

# Anchor M_s, M_q against the literal forms quoted in paper/stages/stage_137.tex.
# Ms and Mq above were built via the route L*rho_c/Theta and -L*sigma_c/Theta;
# the *_paper forms are constructed directly from the paper-card primitives.
# Any sign/factor error in rho_c, sigma_c, Ms, or Mq propagation will fail this.
Ms_paper = L * gs**2 / (Ks * Theta)
Mq_paper = -L * (Ks*gq - lam*gs)**2 / (Ks * (Ks*Kq + lam**2) * Theta)
assert sp.simplify(Ms - Ms_paper) == 0
assert sp.simplify(Mq - Mq_paper) == 0
print('M_s matches paper card closed form.')
print('M_q matches paper card closed form.')

# Schur-complement static-limit anchor (notes Sec. on Stage 97 form).
# delta_Lambda_core(z) = rho_c - sigma_c / (1 - kappa_c z^2 - I gamma_c z^5)
# at z -> 0 must reduce to rho_c - sigma_c.
kappa_c, gamma_c, z_var = sp.symbols('kappa_c gamma_c z_var', positive=True, real=True)
delta_Lambda_core = rho_c - sigma_c / (1 - kappa_c*z_var**2 - sp.I*gamma_c*z_var**5)
static_limit = sp.limit(delta_Lambda_core, z_var, 0)
assert sp.simplify(static_limit - (rho_c - sigma_c)) == 0
print('Schur-complement static limit matches rho_c - sigma_c.')

# Card Check 1 (downgraded to carry-forward at stage 134, but a self-contained
# Sq=0 limit of the gain pair fixed-point law is still well-defined here):
# Pi = M_s + M_q * S_q(Pi). At Sq = 0, recover Pi = M_s.
Pi_var, Sq_var = sp.symbols('Pi_var S_q_var', real=True)
outlet_residual = Pi_var - (Ms + Mq * Sq_var)
assert sp.simplify(outlet_residual.subs(Sq_var, 0).subs(Pi_var, Ms)) == 0
print('Outlet consistency reduces to Pi = M_s at S_q = 0.')

# Check equivalence with r_c notation.
expr_rc = sp.simplify((Ks*gq - lam*gs)**2 / (Ks**2*Kq*(1 + lam**2/(Ks*Kq))))
print('sigma_c (r_c form) =', expr_rc)
assert sp.simplify(sigma_c - expr_rc) == 0

print('\nFinal explicit gain map verified.')
