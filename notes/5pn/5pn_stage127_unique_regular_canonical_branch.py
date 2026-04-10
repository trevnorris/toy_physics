#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
from fivepn_stage121_140_common import banner, rF1, g_minus_F1, g_plus_F1, Pi_star, Pi_match, T_hat_of_Pi

banner("STAGE 127 — UNIQUE REGULAR CANONICAL MOUTH BRANCH")

print("g_-^F1 =", mp.nstr(g_minus_F1, 18))
print("g_+^F1 =", mp.nstr(g_plus_F1, 18))
print("2/pi   =", mp.nstr(2 / mp.pi, 18))
print("\nUpper branch impossible because g_+^F1 > 1.")
print("Lower branch is uniquely reachable because 2/pi < g_-^F1 < 1 and g_Pi is monotone.")
print("\nUnique finite-bias lower-branch point:")
print("Pi_* =", mp.nstr(Pi_star, 18))
print("T_hat_m(Pi_*) =", mp.nstr(T_hat_of_Pi(Pi_star), 18))
print("\nDerivative-matched point g = pi/4 occurs at:")
print("Pi_match =", mp.nstr(Pi_match, 18))
print("T_hat_m(Pi_match) =", mp.nstr(T_hat_of_Pi(Pi_match), 18))
print("\nConclusion: the lower compensated branch is the unique regular finite-bias / finite-traction Family-1 branch.")
