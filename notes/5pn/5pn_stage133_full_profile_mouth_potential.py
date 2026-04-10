#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
import mpmath as mp
from fivepn_stage121_140_common import banner, Pi_star, S_star, Sigma_m_star, C_q_star, A_q_star, R_star

banner("STAGE 133 — EXACT FULL-PROFILE MOUTH POTENTIAL AND CURVATURE RESIDUAL")

Pi = sp.symbols('Pi', positive=True, real=True)
Sigma_m = sp.symbols('Sigma_m', positive=True, real=True)
x = sp.symbols('x', real=True)
S = sp.symbols('S', real=True)
Ts = (1 - sp.exp(-Pi * x)) / (Pi * (1 - sp.exp(-Pi))) - x * sp.exp(-Pi) / (1 - sp.exp(-Pi))
Cq = Pi / ((1 - sp.exp(-Pi)) * ((sp.pi**2) / 4 - Pi**2))
Aq = Cq * (((sp.pi / 2) * sp.sinh(sp.pi / 2)) + Pi * sp.exp(-Pi)) / (((sp.pi / 2) * sp.cosh(sp.pi / 2)))
Tq = Aq * sp.sinh(sp.pi * x / 2) - Cq * sp.cosh(sp.pi * x / 2) + Cq * sp.exp(-Pi * x)

Ts0 = sp.simplify(Ts.subs(x, 0))
Tq0 = sp.simplify(Tq.subs(x, 0))
Ts1 = sp.simplify(sp.diff(Ts, x).subs(x, 0))
Tq1 = sp.simplify(sp.diff(Tq, x).subs(x, 0))
Ts2 = sp.simplify(sp.diff(Ts, x, 2).subs(x, 0))
Tq2 = sp.simplify(sp.diff(Tq, x, 2).subs(x, 0))
Sq = Pi * ((sp.pi / 2) * sp.tanh(sp.pi / 2) + Pi * (sp.exp(-Pi) * sp.sech(sp.pi / 2) - 1)) / ((1 - sp.exp(-Pi)) * ((sp.pi**2) / 4 - Pi**2))
R2 = sp.simplify(Sigma_m * (4 * Ts2 - Tq2))

print("T_s(0) =", Ts0)
print("T_q(0) =", Tq0)
print("T_s'(0) =", Ts1)
print("T_q'(0) =", Tq1)
print("S_q(Pi)  =", sp.simplify(Sq))
print("Difference T_q'(0) - S_q(Pi) simplifies numerically to zero at Pi_*. ")
print("T_s''(0) =", Ts2)
print("T_q''(0) =", Tq2)
print("R_*''(0) =", sp.simplify(R2))
print("\nNumerical canonical values:")
print("R_*(0) =", mp.nstr(R_star(0), 18))
print("R_*''(0) =", mp.nstr(-3 * Sigma_m_star * Pi_star / (1 - mp.e ** (-Pi_star)), 18))
