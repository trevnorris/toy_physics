#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, Sigma0_star, Tm_star, Pi_star, CoevolvingFamily1Solver

banner("STAGE 140 — CORE–MOUTH CO-EVOLUTION STATUS")

solver = CoevolvingFamily1Solver(N=500)
frozen = solver.fixed_point(float(Sigma0_star))
root = solver.canonical_sigma0()

print("At the old canonical traction:")
print("Sigma0_* =", mp.nstr(Sigma0_star, 18), ", T_hat_* =", mp.nstr(Tm_star, 18), ", Pi_* =", mp.nstr(Pi_star, 18))
print("g_fp =", mp.nstr(frozen['g'], 18))
print("R_fp =", mp.nstr(frozen['R'], 18))
print("Pi_fp =", mp.nstr(frozen['Pi'], 18))
print("\nRestoring exact compensation requires the renormalized traction:")
print("Sigma0_can =", mp.nstr(root['Sigma0'], 18))
print("T_hat_can  =", mp.nstr(root['T_hat'], 18))
print("Pi_can     =", mp.nstr(root['Pi'], 18))
print("\nRelative renormalizations:")
print("Sigma0_can/Sigma0_* - 1 =", mp.nstr(root['Sigma0'] / float(Sigma0_star) - 1, 18))
print("T_hat_can/T_hat_* - 1   =", mp.nstr(root['T_hat'] / float(Tm_star) - 1, 18))
print("Pi_can/Pi_* - 1         =", mp.nstr(root['Pi'] / float(Pi_star) - 1, 18))
