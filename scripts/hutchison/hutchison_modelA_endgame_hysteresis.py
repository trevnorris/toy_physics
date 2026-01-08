"""
Hutchison Toy Model (Model A) -- nonlinear endgame with hysteresis
------------------------------------------------------------------
Includes:
  - Duffing nonlinearity (amplitude-dependent resonance shift)
  - Saturating drive coupling
  - Hysteresis via branch selection + continuation (up vs down sweep)

Theory-only thought model; not an experimental blueprint.
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

def duffing_roots_y_for_feff(L, Q, f_drive, beta, f_eff, c_star, chi=1.0):
    w0 = 2*math.pi*f0(L, c_star, chi)
    w  = 2*math.pi*f_drive
    gamma = w0/(2*Q)
    A = (3.0/4.0)*beta
    B = (w0*w0 - w*w)
    C = (2.0*gamma*w)
    coef = [A*A, 2*A*B, (B*B + C*C), -f_eff*f_eff]
    roots = np.roots(coef)
    return sorted([r.real for r in roots if abs(r.imag) < 1e-8 and r.real > 0])

def solve_continuation(
    L, Q, pA, f_grid, beta, a_sat,
    rho_m, alpha_over_beta,
    c_star, chi=1.0,
    direction="up",
    start_branch="low",
    max_iter=25
):
    f_base = (1.0/alpha_over_beta) * (pA/(rho_m*L))
    seq = f_grid[::-1] if direction=="down" else f_grid
    a_prev = 0.0
    out = []

    for idx, fdrv in enumerate(seq):
        a_guess = a_prev
        f_eff = f_base/(1.0 + (a_guess/a_sat)**2)
        a_sel = a_guess
        for _ in range(max_iter):
            y_roots = duffing_roots_y_for_feff(L, Q, fdrv, beta, f_eff, c_star, chi)
            if not y_roots:
                a_cand = 0.0
            else:
                a_cands = [math.sqrt(y) for y in y_roots]
                if idx == 0:
                    a_cand = min(a_cands) if start_branch=="low" else max(a_cands)
                else:
                    a_cand = min(a_cands, key=lambda a: abs(a-a_prev))
            f_new = f_base/(1.0 + (a_cand/a_sat)**2)
            if abs(f_new - f_eff) < 1e-10:
                a_sel = a_cand
                break
            a_sel = a_cand
            f_eff = f_new
        a_prev = a_sel
        out.append(a_sel)

    out = out[::-1] if direction=="down" else out
    return np.array(out)

def eps_rms_from_amplitude(a, L):
    return (a/math.sqrt(2.0))/L

def rigidity_proxy(eps, epscrit, power=2.0):
    return 1.0/(1.0 + (eps/epscrit)**power)
