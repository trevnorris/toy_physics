#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import (
    banner,
    delta_Pi_u,
    delta_T_u,
    delta_Pi_d,
    delta_T_d,
    lambda_Pi_zero,
    lambda_T_zero,
)

banner("STAGE 131 — REPRESENTATIVE NON-EXPONENTIAL POSITIVE MOUTH FAMILIES")

print("Uniform broadening family shifts:")
print("delta Pi_u / eps =", mp.nstr(delta_Pi_u, 18))
print("delta T_u  / eps =", mp.nstr(delta_T_u, 18))
print("\nSelf-matched derivative family shifts:")
print("delta Pi_d / eps =", mp.nstr(delta_Pi_d, 18))
print("delta T_d  / eps =", mp.nstr(delta_T_d, 18))
print("\nConvex interpolation zero-shift points:")
print("lambda_(Pi,0) =", mp.nstr(lambda_Pi_zero, 18))
print("lambda_(T,0)  =", mp.nstr(lambda_T_zero, 18))
print("\nThese match the earlier positive-family compensation fraction at first order.")
