#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, rF1, g_minus_F1, Pi_star, S_star, g_pi

banner("STAGE 122 — ACTUAL FAMILY-1 MOUTH GAINS")

R_nat = (1 - rF1) ** 2 / (1 + rF1**2)
Ms_nat = Pi_star / (1 - R_nat * S_star)
Mq_nat = -R_nat * Ms_nat

R_comp = mp.mpf(1) / 4
Ms_comp = Pi_star / (1 - S_star / 4)
Mq_comp = -Ms_comp / 4

gc_comp = g_minus_F1

print("r_F1 =", mp.nstr(rF1, 18))
print("Pi_* =", mp.nstr(Pi_star, 18))
print("S_q(Pi_*) =", mp.nstr(S_star, 18))
print("\nNatural equal-normalized branch:")
print("R_q^nat =", mp.nstr(R_nat, 18))
print("M_s^nat,* =", mp.nstr(Ms_nat, 18))
print("M_q^nat,* =", mp.nstr(Mq_nat, 18))
print("\nExact compensated branch:")
print("g_c^- =", mp.nstr(gc_comp, 18))
print("R_q =", mp.nstr(R_comp, 18))
print("M_s^comp,* =", mp.nstr(Ms_comp, 18))
print("M_q^comp,* =", mp.nstr(Mq_comp, 18))
print("\nRelative shifts:")
print("shell gain fractional change =", mp.nstr(Ms_comp / Ms_nat - 1, 18))
print("mixed gain magnitude ratio =", mp.nstr(abs(Mq_comp) / abs(Mq_nat), 18))
