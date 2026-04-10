#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, g_minus_F1, CoevolvingFamily1Solver

banner("STAGE 139 — RENORMALIZED CANONICAL BRANCH UNDER FULL CORE–MOUTH CO-EVOLUTION")

solver = CoevolvingFamily1Solver(N=500)
root = solver.canonical_sigma0()

print("Unique renormalized canonical gain:")
print("Sigma0_can =", mp.nstr(root['Sigma0'], 18))
print("T_hat_can  =", mp.nstr(root['T_hat'], 18))
print("\nRestored canonical fixed point:")
print("g_can =", mp.nstr(root['g'], 18), "  (target g_* =", mp.nstr(g_minus_F1, 18), ")")
print("R_can =", mp.nstr(root['R'], 18))
print("S_can =", mp.nstr(root['S'], 18))
print("Pi_can =", mp.nstr(root['Pi'], 18))
print("iterations at root =", root['iterations'])
