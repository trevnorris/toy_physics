#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, rF1, g_minus_F1

banner("STAGE 137 — EXACT CO-EVOLVING CORE–MOUTH FIXED-POINT MAP")

sqrtfac = mp.sqrt(1 + rF1**2)
print("Exact co-evolving map:")
print("g[Sigma] = ∫ Sigma(x) cos(pi x/2) dx")
print("S[Sigma] = ∫ Sigma(x) cosh(pi(1-x)/2)/cosh(pi/2) dx")
print("R[Sigma] = ((g[Sigma] - r_F1)^2)/(1 + r_F1^2)")
print("Phi_{Sigma0}[Sigma](x) = Sigma0 [T_s[Sigma](x) - R[Sigma] T_q[Sigma](x)]")
print("Sigma(x) = exp(-Phi)/∫exp(-Phi)")
print("\nCanonical compensation inside the co-evolving map is g[Sigma] = g_* <=> R[Sigma] = 1/4.")
print("\nExact first-order defect transport:")
print("delta R = - delta g / sqrt(1+r_F1^2) + O(delta g^2)")
print("sqrt(1+r_F1^2) =", mp.nstr(sqrtfac, 18))
print("linear coefficient =", mp.nstr(-1 / sqrtfac, 18))
print("\nLocal slope identity: Pi[Sigma] = Sigma0 [1 - R[Sigma] S[Sigma]]")
print("Under self-matched susceptibility: Sigma0 = (20/9) T_hat_m^2")
