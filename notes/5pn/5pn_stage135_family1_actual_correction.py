#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import (
    banner,
    cov_c_R_star,
    cov_Kq_R_star,
    delta_g_act,
    delta_S_act,
    delta_Pi_act,
    delta_T_act,
    Pi_corr,
    T_corr,
    g1,
    S1,
    Pi1,
    T1,
    lambda_eff_Pi,
    lambda_eff_T,
)

banner("STAGE 135 — ACTUAL FAMILY-1 MOUTH CORRECTION AND ONE-STEP NONLINEAR CHECK")

print("Residual covariances:")
print("Cov_*(c, R_*) =", mp.nstr(cov_c_R_star, 18))
print("Cov_*(K_q, R_*) =", mp.nstr(cov_Kq_R_star, 18))
print("\nActual first-order moment shifts:")
print("delta g_act =", mp.nstr(delta_g_act, 18))
print("delta S_act =", mp.nstr(delta_S_act, 18))
print("delta Pi_act =", mp.nstr(delta_Pi_act, 18))
print("delta T_act  =", mp.nstr(delta_T_act, 18))
print("Pi_corr =", mp.nstr(Pi_corr, 18))
print("T_corr  =", mp.nstr(T_corr, 18))
print("\nOne-step nonlinear iterate:")
print("g_1 =", mp.nstr(g1, 18))
print("S_1 =", mp.nstr(S1, 18))
print("Pi_1 =", mp.nstr(Pi1, 18))
print("T_1  =", mp.nstr(T1, 18))
print("\nEffective positive-family interpolation parameters:")
print("lambda_eff^(Pi) =", mp.nstr(lambda_eff_Pi, 18))
print("lambda_eff^(T)  =", mp.nstr(lambda_eff_T, 18))
