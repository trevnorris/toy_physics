#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, rF1, Pi_star, S_star

banner("STAGE 124 — MOUTH GAIN STATUS UPDATE")

R_nat = (1 - rF1) ** 2 / (1 + rF1**2)
Ms_nat = Pi_star / (1 - R_nat * S_star)
Mq_nat = -R_nat * Ms_nat
Ms_comp = Pi_star / (1 - S_star / 4)
Mq_comp = -Ms_comp / 4
T_nat = mp.sqrt(9 * Ms_nat / 20)
T_comp = mp.sqrt(9 * Ms_comp / 20)

print("Derived gain formulas:")
print("M_s = L g_s^2 / (K_s Theta_sigma)")
print("M_q = -L (K_s g_q - lambda g_s)^2 / [K_s(K_s K_q + lambda^2) Theta_sigma]")
print("\nExact compensated ratio: M_q = - M_s / 4")
print("\nFamily-1 numbers:")
print("natural branch:    M_s =", mp.nstr(Ms_nat, 12), ", M_q =", mp.nstr(Mq_nat, 12))
print("compensated branch M_s =", mp.nstr(Ms_comp, 12), ", M_q =", mp.nstr(Mq_comp, 12))
print("\nSelf-matched traction law: M_s = (20/9) T_hat_m^2")
print("natural vs canonical traction difference =", mp.nstr(T_comp / T_nat - 1, 12))
