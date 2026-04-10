#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, Pi_corr, T_corr, Pi1, T1

banner("STAGE 136 — FULL MOUTH CORRECTION STATUS")

print("The full explicit GNLS + localized-Maxwell Family-1 mouth boundary layer selects a definite non-exponential correction.")
print("Corrected canonical point from the linear projection:")
print("Pi_corr =", mp.nstr(Pi_corr, 18))
print("T_corr  =", mp.nstr(T_corr, 18))
print("\nOne-step nonlinear iterate:")
print("Pi_1 =", mp.nstr(Pi1, 18))
print("T_1  =", mp.nstr(T1, 18))
print("\nInterpretation: branch selection is finished; what remains is a finite correction problem around the unique regular branch.")
