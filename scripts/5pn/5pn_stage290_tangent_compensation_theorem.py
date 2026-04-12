#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage289_292_common import *

"""
Stage 290 — exact tangent-compensation theorem for arbitrary first-order isotropic
bundle drift.
"""

banner("STAGE 290 — TANGENT-COMPENSATION THEOREM")

refs = family1_reference_values()
g_star = refs["g_star"]
r_star = refs["r_F1"]
Sigma0_can = refs["Sigma0_can"]
S_can = refs["S_can"]

th, Ks, Kq, P0 = sp.symbols("dln_Theta_w dln_K_s dln_K_q dln_P_0", real=True)

rho_w = th / 2
csw = th
ell = -th
a = Ks / 2 - th / 4
LW = a
cs = Ks / 2 - th / 4 + P0 / 5
Zq = Kq - 2 * P0 / 5
vw0 = -sp.Rational(3, 4) * Ks + sp.Rational(1, 2) * Kq + sp.Rational(13, 8) * th
Tm = -sp.Rational(5, 4) * Ks + sp.Rational(1, 2) * Kq + sp.Rational(15, 8) * th - sp.Rational(2, 5) * P0
gs = -sp.Rational(1, 4) * Ks + sp.Rational(1, 2) * Kq + sp.Rational(3, 8) * th - sp.Rational(2, 5) * P0
gq = -sp.Rational(3, 4) * Ks + Kq + sp.Rational(3, 8) * th - sp.Rational(2, 5) * P0
lam = sp.Rational(1, 2) * (Ks + Kq)

subbanner("I. Parent-invariant drift ledger")
dln_rc = sp.simplify(2 * lam - Ks - Kq)
dln_r = sp.simplify(lam - (Ks + Kq) / 2)
dln_g = sp.simplify(gq + Ks / 2 - gs - Kq / 2)

print("dln r_c =")
sp.pprint(dln_rc)
print("dln r =")
sp.pprint(dln_r)
print("dln g =")
sp.pprint(dln_g)

expect_zero("hybridization ratio invariant", dln_rc)
expect_zero("normalized mixed flow invariant", dln_r)
expect_zero("normalized mouth-coupling invariant", dln_g)

subbanner("II. Stage-147 logarithmic imbalance channels")
imbalance_1 = sp.simplify(gq + Ks - gs - lam)
imbalance_2 = sp.simplify(Ks + Kq - 2 * lam)
print("dln(g_q K_s/(g_s lambda)) =")
sp.pprint(imbalance_1)
print("dln(K_s K_q/lambda^2) =")
sp.pprint(imbalance_2)
expect_zero("imbalance channel 1", imbalance_1)
expect_zero("imbalance channel 2", imbalance_2)

subbanner("III. Exact off-family normal coordinate")
delta_perp = sp.simplify(g_star * imbalance_1 + imbalance_2 / (4 * sp.sqrt(1 + r_star**2)))
print("delta_perp =")
sp.pprint(delta_perp)
expect_zero("delta_perp vanishes on isotropic bundle transport", delta_perp)

subbanner("IV. Mouth-bias transport collapses to the tangential piece")
dSigma0, dScal = sp.symbols("dSigma0 dScal", real=True)
dPi_tan = sp.simplify((1 - sp.Rational(1, 4) * S_can) * dSigma0 - Sigma0_can * dScal / 4)
dPi_full = sp.simplify(dPi_tan + Sigma0_can * S_can * delta_perp / sp.sqrt(1 + r_star**2))
print("deltaPi_tan =")
sp.pprint(dPi_tan)
print("deltaPi_full =")
sp.pprint(dPi_full)
expect_zero("deltaPi_full - deltaPi_tan", sp.simplify(dPi_full - dPi_tan))

subbanner("V. Interpretation")
print("Arbitrary first-order isotropic bundle drift in")
print("  (Theta_w, K_s, K_q, P_0)")
print("preserves r_c, r, and g exactly, so it moves the co-evolving Family-1 branch")
print("tangentially inside the exact parent compensation family.")
print("Therefore the first genuine first-order danger cannot come from isotropic bundle")
print("transport itself. It must come from off-bundle structure.")
