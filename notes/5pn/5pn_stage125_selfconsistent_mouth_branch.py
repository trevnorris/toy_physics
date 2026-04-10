#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, rF1, Pi_star, S_star, g_pi, R_q, Sigma0_of_Pi, T_hat_of_Pi, Tm_star

banner("STAGE 125 — SELF-CONSISTENT MOUTH-BRANCH LAW")

print("Self-consistent gain ratio:")
print("R_q(Pi) = ((g_Pi - r_F1)^2)/(1 + r_F1^2)")
print("with g_Pi =", mp.nstr(g_pi(Pi_star), 18), "at Pi_*")
print("\nExact shell-gain branch:")
print("Sigma_0(Pi) = Pi / (1 - R_q(Pi) S_q(Pi))")
print("T_hat_m(Pi) = sqrt(9 Pi / [20 (1 - R_q(Pi) S_q(Pi))])")
print("\nCanonical point:")
print("Pi_* =", mp.nstr(Pi_star, 18))
print("R_q(Pi_*) =", mp.nstr(R_q(Pi_star), 18))
print("Sigma_0(Pi_*) =", mp.nstr(Sigma0_of_Pi(Pi_star), 18))
print("T_hat_m(Pi_*) =", mp.nstr(T_hat_of_Pi(Pi_star), 18))
print("(Matches earlier canonical traction =", mp.nstr(Tm_star, 18), ")")
