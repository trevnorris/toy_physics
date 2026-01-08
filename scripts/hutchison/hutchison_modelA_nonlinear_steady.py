
"""
Hutchison Toy Model (Model A) -- nonlinear upgrade (steady-state)
----------------------------------------------------------------
Implements:
  - Duffing amplitude-dependent frequency shift
  - Saturating coupling in the drive
  - Rigidity proxy R(eps) = 1/(1+(eps/epscrit)^2)

This provides fast (f, pA) response maps and threshold contours.

NOT a blueprint for experiments; theory-only thought model.
"""

import math
import numpy as np

def f0(L, c_star, chi=1.0):
    return chi*c_star/(2*L)

def duffing_beta_for_linewidth(L, Q, eps_ref, c_star, chi=1.0, sign=+1.0):
    w0 = 2*math.pi*f0(L, c_star, chi)
    a = math.sqrt(2.0)*eps_ref*L
    beta_mag = (4.0*w0*w0) / (3.0*Q*a*a)
    return sign*beta_mag

def solve_duffing_amplitude(
    L, Q, pA, f_drive,
    beta, a_sat,
    rho_m, alpha_over_beta,
    c_star, chi=1.0,
    max_iter=30, tol=1e-8
):
    w0 = 2*math.pi*f0(L, c_star, chi)
    w  = 2*math.pi*f_drive
    gamma = w0/(2*Q)
    f_base = (1.0/alpha_over_beta) * (pA/(rho_m*L))
    a = 0.0
    f_eff = f_base

    for _ in range(max_iter):
        A = (3.0/4.0)*beta
        B = (w0*w0 - w*w)
        C = (2.0*gamma*w)
        coef = [A*A, 2*A*B, (B*B + C*C), -f_eff*f_eff]
        roots = np.roots(coef)
        real_pos = [r.real for r in roots if abs(r.imag) < 1e-8 and r.real > 0]
        a_new = math.sqrt(max(real_pos)) if real_pos else 0.0
        f_new = f_base / (1.0 + (a_new/a_sat)**2)
        if abs(a_new - a) < tol*max(1.0, abs(a_new)):
            a, f_eff = a_new, f_new
            break
        a, f_eff = a_new, f_new

    return a

def eps_from_amplitude(a, L):
    return (a/math.sqrt(2.0))/L

def rigidity_proxy(eps, epscrit, power=2.0):
    return 1.0/(1.0 + (eps/epscrit)**power)
