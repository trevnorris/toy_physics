#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, A_T, B_T, Tm_star, g_minus_F1, S_star, Pi_star

banner("STAGE 130 — FIRST-ORDER RIGIDITY KERNEL")

print("Canonical branch traction law:")
print("T_hat_m = sqrt(9 Sigma_0 / 20),  Sigma_0 = Pi / (1 - S_q/4)")
print("\nRigidity coefficients:")
print("A_T =", mp.nstr(A_T, 18))
print("B_T =", mp.nstr(B_T, 18))
print("|A_T|/B_T =", mp.nstr(abs(A_T) / B_T, 18))
print("\nFirst-order traction law:")
print("delta T_hat_m = eps [ A_T (ḡ_varsigma - g_*) + B_T (S̄_varsigma - S_*) ]")
print("\nEquivalent centered kernel form:")
print("W_*(x) = A_T (c(x) - g_*) + B_T (K_q(x) - S_*)")
print("delta T_hat_m = eps ∫ W_*(x) [varsigma(x) - Sigma_*(x)] dx")
