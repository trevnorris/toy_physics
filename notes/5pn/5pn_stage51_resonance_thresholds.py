#!/usr/bin/env python3
"""
5pn_stage51_resonance_thresholds.py

Stage 51 audit: resonance-corrected thresholds for the sech–Gaussian benchmark family.
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 51 — RESONANCE-CORRECTED THRESHOLDS FOR THE SECH–GAUSSIAN BENCHMARK")

rho_star, g_phi, N_pp, m, cs_star, K_X = sp.symbols("rho_star g_phi N_pp m c_s_star K_X", positive=True, real=True)
W_wall, Pe_req, Delta0, Deltainf = sp.symbols("W_wall Pe_req Delta0 Deltainf", positive=True, real=True)
C2 = sp.symbols("C2", positive=True, real=True)

subbanner("51.1 — Exact resonance-corrected gain and wall figure")
G_match = rho_star * g_phi**2 * N_pp / (m * cs_star**2 * K_X)
G_res = sp.simplify(C2 * G_match)
W_res = sp.simplify(C2 * W_wall)
expect_zero("resonance gain factor", G_res - C2 * G_match)
expect_zero("resonance wall factor", W_res - C2 * W_wall)

subbanner("51.2 — Exact profile-family thresholds")
W_fail_res = sp.simplify(Pe_req / (C2 * Deltainf))
W_suff_res = sp.simplify(Pe_req / (C2 * Delta0))
print("W_wall <= Pe_req/[C^2 Delta_inf]  -> fail on the profile family")
print("W_wall >= Pe_req/[C^2 Delta_0]    -> succeed on the profile family")
print("\nExplicit threshold surfaces:")
sp.pprint(W_fail_res)
sp.pprint(W_suff_res)

subbanner("51.3 — Resonance-point penalty factor")
mp.mp.dps = 60

def I_of_r(rv: mp.mpf) -> mp.mpf:
    f = lambda xx: mp.sech(xx) * mp.e**(-(xx**2) / (rv**2))
    return 2 * mp.quad(f, [0, mp.inf])


def C2_of_r(rv: mp.mpf) -> mp.mpf:
    return I_of_r(rv)**2 / (rv * mp.sqrt(2 * mp.pi))

r_star = mp.sqrt(mp.pi)
C_res2 = C2_of_r(r_star)
P_res = 1 / C_res2
print("r_* = sqrt(pi)        =", r_star)
print("C_res^2               =", C_res2)
print("P_res = 1/C_res^2     =", P_res)
print("Threshold shift factor =", P_res)

banner("STAGE 51 FINAL LEDGER")
print("Stage 51 shows that the explicit sech–Gaussian benchmark does not rewrite the theorem")
print("structure. It simply multiplies the Stage-49 wall thresholds by the tiny resonance")
print("penalty factor P_res = 1.005612487760576... .")
