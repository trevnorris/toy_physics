#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, Sigma0_star, Tm_star, g_minus_F1, rF1, CoevolvingFamily1Solver

banner("STAGE 138 — FAMILY-1 CO-EVOLVING FIXED POINT AT FROZEN CANONICAL TRACTION")

solver = CoevolvingFamily1Solver(N=500)
res = solver.fixed_point(float(Sigma0_star))

delta_g = res['g'] - float(g_minus_F1)
delta_R_pred = -delta_g / float(mp.sqrt(1 + rF1**2)) + delta_g**2 / float(1 + rF1**2)

print("Sigma0_* =", mp.nstr(Sigma0_star, 18), "  (T_hat_* =", mp.nstr(Tm_star, 18), ")")
print("iterations =", res['iterations'], ", max update =", res['err'])
print("g_fp =", mp.nstr(res['g'], 18))
print("S_fp =", mp.nstr(res['S'], 18))
print("R_fp =", mp.nstr(res['R'], 18))
print("Pi_fp =", mp.nstr(res['Pi'], 18))
print("\ndelta g_fp =", mp.nstr(delta_g, 18))
print("delta R_pred from exact transport =", mp.nstr(delta_R_pred, 18))
print("delta R_actual =", mp.nstr(res['R'] - 0.25, 18))
