#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, Pi_star, Tm_star, g_minus_F1, g_pi

banner("STAGE 128 — MOUTH BRANCH SELECTION STATUS")

print("Positive mouth-layer family:")
print("g_Pi = 2 Pi (2 Pi e^Pi + pi) / ((4 Pi^2 + pi^2)(e^Pi - 1))")
print("range: 2/pi < g_Pi < 1 for Pi > 0")
print("\nSelf-consistent Family-1 law:")
print("Pi = Sigma_0 [1 - R_q(Pi) S_q(Pi)]")
print("\nCanonical lower branch:")
print("g_* =", mp.nstr(g_minus_F1, 18))
print("Pi_* =", mp.nstr(Pi_star, 18))
print("T_hat_* =", mp.nstr(Tm_star, 18))
print("\nWhat remains open is no longer branch choice but the size of finite corrections to the unique regular canonical point.")
