#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage289_292_common import *

"""
Stage 289 — exact bundle transport compiler from (Theta_w, K_s, K_q, P_0)
into the remaining lower-branch microscopic drifts.
"""

banner("STAGE 289 — BUNDLE TRANSPORT COMPILER")

th, Ks, Kq, P0 = sp.symbols("dln_Theta_w dln_K_s dln_K_q dln_P_0", real=True)
N0, D0 = sp.symbols("dln_N_0 dln_D_0", real=True)

rho_w, a, cs, Zq = sp.symbols("dln_rho_w dln_a dln_c_s dln_Z_q", real=True)
csw, ell, LW = sp.symbols("dln_c_sw dln_ell dln_L_W", real=True)
vw0, Tm, gs, gq, lam = sp.symbols("dln_v_w0 dln_T_m dln_g_s dln_g_q dln_lambda", real=True)

subbanner("I. Exact inversion of the last four irreducible drifts")
eqs = [
    sp.Eq(th, 2 * rho_w),
    sp.Eq(Ks, 2 * a + rho_w),
    sp.Eq(Kq, Zq + 2 * cs - 2 * a),
    sp.Eq(P0, 5 * (cs - a)),
]
sol = sp.solve(eqs, [rho_w, a, cs, Zq], dict=True)[0]
for key in [rho_w, a, cs, Zq]:
    print(f"{key} =")
    sp.pprint(sol[key])

expect_zero("rho_w - Theta_w/2", sol[rho_w] - th / 2)
expect_zero("a - (K_s/2 - Theta_w/4)", sol[a] - (Ks / 2 - th / 4))
expect_zero("c_s - (K_s/2 - Theta_w/4 + P_0/5)", sol[cs] - (Ks / 2 - th / 4 + P0 / 5))
expect_zero("Z_q - (K_q - 2 P_0/5)", sol[Zq] - (Kq - 2 * P0 / 5))

subbanner("II. Exact co-transport of the remaining mouth/background variables")
transport = {
    csw: th,
    ell: -th,
    LW: Ks / 2 - th / 4,
    vw0: -sp.Rational(3, 4) * Ks + sp.Rational(1, 2) * Kq + sp.Rational(13, 8) * th,
    Tm: -sp.Rational(5, 4) * Ks + sp.Rational(1, 2) * Kq + sp.Rational(15, 8) * th - sp.Rational(2, 5) * P0,
    gs: -sp.Rational(1, 4) * Ks + sp.Rational(1, 2) * Kq + sp.Rational(3, 8) * th - sp.Rational(2, 5) * P0,
    gq: -sp.Rational(3, 4) * Ks + Kq + sp.Rational(3, 8) * th - sp.Rational(2, 5) * P0,
    lam: sp.Rational(1, 2) * (Ks + Kq),
}
for key in [csw, ell, LW, vw0, Tm, gs, gq, lam]:
    print(f"{key} =")
    sp.pprint(transport[key])

subbanner("III. Full-bundle rewrite using P_0 = N_0 / D_0")
P0_bundle = N0 - D0
cs_bundle = sp.simplify(sol[cs].subs({P0: P0_bundle}))
Zq_bundle = sp.simplify(sol[Zq].subs({P0: P0_bundle}))
print("dln c_s =")
sp.pprint(cs_bundle)
print("dln Z_q =")
sp.pprint(Zq_bundle)
expect_zero(
    "bundle c_s rewrite",
    cs_bundle - (Ks / 2 - th / 4 + (N0 - D0) / 5),
)
expect_zero(
    "bundle Z_q rewrite",
    Zq_bundle - (Kq - 2 * (N0 - D0) / 5),
)

subbanner("IV. Frozen-wall corollary")
print("If dln Theta_w = 0, then")
print("dln rho_w =")
sp.pprint(sp.simplify(sol[rho_w].subs({th: 0})))
print("dln a =")
sp.pprint(sp.simplify(sol[a].subs({th: 0})))
print("dln c_s =")
sp.pprint(sp.simplify(sol[cs].subs({th: 0})))
print("dln Z_q =")
sp.pprint(sp.simplify(sol[Zq].subs({th: 0})))

subbanner("V. Interpretation")
print("The last four irreducible lower-branch drifts are exact algebraic images of")
print("  (dln Theta_w, dln K_s, dln K_q, dln P_0).")
print("Everything else is co-transported by the same bundle data.")
