#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage121_140_common import banner

banner("STAGE 121 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO")

r, gc, Sigma0 = sp.symbols('r gc Sigma0', positive=True, real=True)
Ks, Kq, lam, gs, gq, L, Theta = sp.symbols('K_s K_q lambda g_s g_q L Theta_sigma', positive=True, real=True)

Ms = sp.simplify(L * gs**2 / (Ks * Theta))
Mq = sp.simplify(-L * (Ks * gq - lam * gs)**2 / (Ks * (Ks * Kq + lam**2) * Theta))

r_def = lam / sp.sqrt(Ks * Kq)
gc_def = gq * sp.sqrt(Ks) / (gs * sp.sqrt(Kq))
Sigma0_def = L * gs**2 / (Ks * Theta)

expr_Ms = Sigma0
expr_Mq = -Sigma0 * (gc - r)**2 / (1 + r**2)
Rq = sp.simplify(-expr_Mq / expr_Ms)

print("M_s =", expr_Ms)
print("M_q =", expr_Mq)
print("R_q =", Rq)
print("Pi = Sigma_0 * (1 - R_q * S_q(Pi))")

# Exact compensation family
branch = r + sp.Rational(1, 2) * sp.sqrt(1 + r**2)
branch2 = r - sp.Rational(1, 2) * sp.sqrt(1 + r**2)
Rq_plus = sp.simplify(Rq.subs(gc, branch))
Rq_minus = sp.simplify(Rq.subs(gc, branch2))
print("R_q on upper compensation family =", Rq_plus)
print("R_q on lower compensation family =", Rq_minus)

assert sp.simplify(Rq_plus - sp.Rational(1, 4)) == 0
assert sp.simplify(Rq_minus - sp.Rational(1, 4)) == 0

print("\nResult: M_s = Sigma_0,  M_q = -Sigma_0 * ((g_c-r)^2/(1+r^2)),  and on the exact compensation family R_q = 1/4.")
