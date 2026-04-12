#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, A_T, B_T, lambda_Pi_zero, lambda_T_zero

banner("STAGE 132 — NON-EXPONENTIAL MOUTH-RIGIDITY STATUS")

print("At first order, every positive normalized deformation enters only through two moments:")
print("ḡ_varsigma = ∫ varsigma c,   S̄_varsigma = ∫ varsigma K_q")
print("\nRetuning and traction transport:")
print("delta Pi = - eps (ḡ_varsigma - g_*) / g_*'")
print("delta T_hat = eps [A_T (ḡ_varsigma - g_*) + B_T (S̄_varsigma - S_*)]")
print("A_T =", mp.nstr(A_T, 12), ", B_T =", mp.nstr(B_T, 12), ", |A_T|/B_T =", mp.nstr(abs(A_T) / B_T, 12))
print("\nRepresentative zero-shift mixtures:")
print("lambda_(Pi,0) =", mp.nstr(lambda_Pi_zero, 12))
print("lambda_(T,0) =", mp.nstr(lambda_T_zero, 12))
print("\nInterpretation: the canonical branch is rigid but not brittle; overlap changes dominate mixed-kernel changes by ~31.68x.")
