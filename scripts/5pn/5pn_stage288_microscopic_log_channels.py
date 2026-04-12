#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage284_288_common import *

"""
Stage 288 — explicit microscopic log channels and lower-branch cancellation.

What this script does
---------------------
1. Rewrites the off-family scalar delta_perp in the Stage-147 microscopic
   log-channel variables (Z_q, rho_w, c_{s,w}, c_s, T_m, v_{w0}, a, L_W).
2. Inserts the exact lower compensated branch drift laws from Stage 148.
3. Verifies that both imbalance channels vanish identically on the exact lower
   branch, hence delta_perp = 0 on-family.
"""

banner("STAGE 288 — MICROSCOPIC LOG CHANNELS AND LOWER-BRANCH CANCELLATION")

refs = family1_reference_values()
g_star = refs["g_star"]
r_F1 = refs["r_F1"]

ch = tangent_log_channels()
Zq = ch["dln_Zq"]
rho_w = ch["dln_rho_w"]
csw = ch["dln_csw"]
cs = ch["dln_cs"]
Tm = ch["dln_Tm"]
vw0 = ch["dln_vw0"]
a = ch["dln_a"]
LW = ch["dln_LW"]

subbanner("I. Stage-147 microscopic imbalance channels")
imbalance_1 = sp.simplify(Zq + 3 * csw - rho_w - Tm - vw0 - 2 * a - 2 * LW)
imbalance_2 = sp.simplify(Zq + 2 * cs + 3 * csw - rho_w - 2 * vw0 - 2 * a - 3 * LW)

print("delta ln(g_q K_s / (g_s lambda)) =")
sp.pprint(imbalance_1)
print("delta ln(K_s K_q / lambda^2) =")
sp.pprint(imbalance_2)

coef = sp.N(1 / (4 * sp.sqrt(1 + r_F1**2)))
delta_perp = sp.simplify(g_star * imbalance_1 + imbalance_2 / (4 * sp.sqrt(1 + r_F1**2)))
print("delta_perp =")
sp.pprint(delta_perp)
print(f"Numerical coefficient 1/(4 sqrt(1+r_F1^2)) = {coef}")

subbanner("II. Exact lower-branch transport laws from Stage 148")
lLW = a
lvw0 = sp.Rational(1, 2) * (Zq - rho_w + 3 * csw + 2 * cs - 5 * a)
lTm = sp.Rational(1, 2) * (Zq - rho_w + 3 * csw - 2 * cs - 3 * a)

print("delta ln L_W =")
sp.pprint(lLW)
print("delta ln v_{w0} =")
sp.pprint(lvw0)
print("delta ln T_m =")
sp.pprint(lTm)

subbanner("III. Exact cancellation on the lower compensated branch")
imbalance_1_lower = sp.simplify(imbalance_1.subs({LW: lLW, vw0: lvw0, Tm: lTm}))
imbalance_2_lower = sp.simplify(imbalance_2.subs({LW: lLW, vw0: lvw0}))
delta_perp_lower = sp.simplify(delta_perp.subs({LW: lLW, vw0: lvw0, Tm: lTm}))

expect_zero("lower-branch imbalance channel 1", imbalance_1_lower)
expect_zero("lower-branch imbalance channel 2", imbalance_2_lower)
expect_zero("lower-branch delta_perp", delta_perp_lower)

subbanner("IV. Final interpretation")
print("On the exact lower compensated branch, the co-transport laws force")
print("  delta ln(g_q K_s/(g_s lambda)) = 0,")
print("  delta ln(K_s K_q/lambda^2) = 0,")
print("hence delta_perp = 0 identically.")
