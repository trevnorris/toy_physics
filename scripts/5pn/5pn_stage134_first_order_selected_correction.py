#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, cov_c_R_star, cov_Kq_R_star, gprime_star, A_T, B_T

banner("STAGE 134 — FIRST-ORDER SOURCE CORRECTION SELECTED BY THE FULL MOUTH PROFILE")

print("Exact first-order selected correction:")
print("Sigma_act(x) = Sigma_*(x) [1 - (R_*(x) - <R_*>_*)] + O(R_*^2)")
print("\nSelected moment shifts:")
print("delta g_act = - Cov_*(c, R_*) =", mp.nstr(-cov_c_R_star, 18))
print("delta S_act = - Cov_*(K_q, R_*) =", mp.nstr(-cov_Kq_R_star, 18))
print("\nRetuning and traction transport:")
print("delta Pi_act = Cov_*(c, R_*) / g_*' =", mp.nstr(cov_c_R_star / gprime_star, 18))
print("delta T_act  = -A_T Cov_*(c,R_*) - B_T Cov_*(K_q,R_*)")
print("            =", mp.nstr(-A_T * cov_c_R_star - B_T * cov_Kq_R_star, 18))
