#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, Pi_star, g_minus_F1, S_star, gprime_star, Sprime_star, c_kernel, Kq_kernel, g_u, S_u, g_d, S_d

banner("STAGE 129 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS")

print("Canonical point:")
print("Pi_* =", mp.nstr(Pi_star, 18))
print("g_* =", mp.nstr(g_minus_F1, 18))
print("S_* =", mp.nstr(S_star, 18))
print("g_*' =", mp.nstr(gprime_star, 18))
print("S_*' =", mp.nstr(Sprime_star, 18))
print("\nFor any positive normalized deformation family Sigma_eps = (1-eps) Sigma_* + eps varsigma,")
print("the only first-order moments are ḡ_varsigma and S̄_varsigma.")
print("Uniform family:")
print("g_u =", mp.nstr(g_u, 18))
print("S_u =", mp.nstr(S_u, 18))
print("Derivative family:")
print("g_d =", mp.nstr(g_d, 18))
print("S_d =", mp.nstr(S_d, 18))
print("\nExact first-order retuning law:")
print("delta Pi = - eps * (ḡ_varsigma - g_*) / g_*'")
print("delta S_q = eps * [ (S̄_varsigma - S_*) - (S_*'/g_*') (ḡ_varsigma - g_*) ]")
