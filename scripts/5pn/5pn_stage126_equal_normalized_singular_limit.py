#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
import mpmath as mp
from fivepn_stage121_140_common import banner, rF1, g_pi, S_q, R_q, T_hat_of_Pi

banner("STAGE 126 — EQUAL-NORMALIZED BRANCH IS A SINGULAR LIMIT")

Pi = sp.symbols('Pi', positive=True, real=True)
g_expr = 2 * Pi * (2 * Pi * sp.exp(Pi) + sp.pi) / ((4 * Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1))
num = sp.expand(sp.pi**2 * (sp.exp(Pi) - 1) - 2 * sp.pi * Pi - 4 * Pi**2)
decomp = sp.expand(sp.pi**2 * (sp.exp(Pi) - 1 - Pi - Pi**2 / 2) + Pi * (sp.pi**2 - 2 * sp.pi) + Pi**2 * (sp.pi**2 / 2 - 4))
assert sp.simplify(num - decomp) == 0

print("1 - g_Pi = numerator / denominator with")
print("numerator =", num)
print("decomposition =", decomp)
print("denominator = (4 Pi^2 + pi^2)(e^Pi - 1) > 0")
print("Hence 0 < g_Pi < 1 for all finite Pi > 0.")

for P in [2, 5, 10, 20, 50]:
    print(f"Pi={P:>2}  g_Pi={mp.nstr(g_pi(P), 18)}  T_hat(Pi)={mp.nstr(T_hat_of_Pi(P), 18)}")

print("\nLimit Pi -> infinity: g_Pi ->", mp.nstr(mp.limit(lambda z: g_pi(z), mp.inf), 18))
print("Equal-normalized branch g_c = 1 occurs only in the singular point-source limit Pi -> +infty.")
