#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
import mpmath as mp
from fivepn_stage121_140_common import banner, Pi_star, rF1, S_star

banner("STAGE 123 — SELF-MATCHED MOUTH SUSCEPTIBILITY CLOSURE")

L, ell, rho_w, Tm, hbar, cs = sp.symbols('L ell rho_w T_m hbar c_s_w', positive=True, real=True)
J_s, H_w, g_s, K_s, Theta = sp.symbols('J_s H_w g_s K_s Theta_sigma', positive=True, real=True)

J_s_expr = 4 * sp.pi * sp.Symbol('a', positive=True)**2 * ell / 3
H_w_expr = sp.Symbol('m_psi', positive=True) * cs**2 / rho_w
Theta_expr = sp.simplify(H_w_expr * J_s_expr)
g_s_expr = Tm * J_s_expr
K_s_expr = 3 * sp.pi * sp.Symbol('a', positive=True)**2 * hbar**2 / (5 * sp.Symbol('m_psi', positive=True) * rho_w * ell)

Sigma0_expr = sp.simplify(sp.Symbol('L', positive=True) * g_s_expr**2 / (K_s_expr * Theta_expr))
Sigma0_simplified = sp.simplify(Sigma0_expr.subs(sp.Symbol('L', positive=True), L))
That = sp.symbols('That_m', positive=True, real=True)
That_def = rho_w * ell * sp.sqrt(L) * Tm / (hbar * cs)
Sigma0_That = sp.simplify(Sigma0_simplified.subs(Tm, That * hbar * cs / (rho_w * ell * sp.sqrt(L))))

print("Sigma_0 =", Sigma0_simplified)
print("With T_hat := rho_w ell sqrt(L) T_m / (hbar c_s,w),  Sigma_0 =", Sigma0_That)

# Numerical targets from Stage 122
R_nat = (1 - rF1) ** 2 / (1 + rF1**2)
Ms_nat = Pi_star / (1 - R_nat * S_star)
Ms_comp = Pi_star / (1 - S_star / 4)
T_nat = mp.sqrt(9 * Ms_nat / 20)
T_comp = mp.sqrt(9 * Ms_comp / 20)

print("\nRequired normalized tractions on Family-1:")
print("T_hat_nat,* =", mp.nstr(T_nat, 18))
print("T_hat_comp,* =", mp.nstr(T_comp, 18))
print("fractional enhancement =", mp.nstr(T_comp / T_nat - 1, 18))
