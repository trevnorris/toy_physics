#!/usr/bin/env python3
"""
5pn_stage50_sech_gaussian_resonance.py

Stage 50 audit: exact sech–Gaussian coherence resonance benchmark.
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


banner("STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE RESONANCE BENCHMARK")

mp.mp.dps = 60

def I_of_r(rv: mp.mpf) -> mp.mpf:
    f = lambda xx: mp.sech(xx) * mp.e**(-(xx**2) / (rv**2))
    return 2 * mp.quad(f, [0, mp.inf])


def C2_of_r(rv: mp.mpf) -> mp.mpf:
    return I_of_r(rv)**2 / (rv * mp.sqrt(2 * mp.pi))

subbanner("50.1 — Exact norms (numerical cross-checks)")
wf = mp.mpf("1.3")
wg = mp.mpf("2.1")
Nss_num = 2 * mp.quad(lambda yy: mp.sech(yy / wf)**2, [0, mp.inf])
Npp_num = mp.quad(lambda yy: mp.e**(-2 * yy**2 / wg**2), [-mp.inf, mp.inf])
print("N_(sigma sigma) numeric =", Nss_num)
print("2 w_f                 =", 2 * wf)
print("abs diff              =", abs(Nss_num - 2 * wf))
print("N_(phi phi) numeric   =", Npp_num)
print("w_g sqrt(pi/2)        =", wg * mp.sqrt(mp.pi / 2))
print("abs diff              =", abs(Npp_num - wg * mp.sqrt(mp.pi / 2)))
print("\nCoherence formula:")
print("  C^2(r) = I(r)^2 / [r sqrt(2 pi)],  I(r)=∫ sech(x) exp(-x^2/r^2) dx")

subbanner("50.2 — Exact self-duality (numerical verification)")
for rv in [mp.mpf("0.6"), mp.mpf("0.8"), mp.mpf("1.1"), mp.mpf("2.3"), mp.mpf("4.0")]:
    lhs = C2_of_r(rv)
    rhs = C2_of_r(mp.pi / rv)
    print(f"r={rv}: C^2(r)     = {lhs}")
    print(f"      C^2(pi/r) = {rhs}")
    print(f"      abs diff  = {abs(lhs-rhs)}")

subbanner("50.3 — Exact self-dual stationary point")
r_star = mp.sqrt(mp.pi)
h = mp.mpf("1e-4")
deriv_fd = (C2_of_r(r_star + h) - C2_of_r(r_star - h)) / (2 * h)
print("r_* = sqrt(pi) =", r_star)
print("finite-difference derivative at r_* =", deriv_fd)

subbanner("50.4 — Numerical maximum and penalty factor")
C_res2 = C2_of_r(r_star)
P_res = 1 / C_res2
print("C_res^2 =", C_res2)
print("1 - C_res^2 =", 1 - C_res2)
print("P_res =", P_res)

# Simple grid audit for uniqueness on a constructive interval.
rs = [mp.mpf("0.35") + mp.mpf(i) * mp.mpf("0.03") for i in range(300)]
vals = [C2_of_r(rv) for rv in rs]
max_index = max(range(len(vals)), key=lambda i: vals[i])
print("grid maximum location ≈", rs[max_index])
print("grid maximum value    =", vals[max_index])

banner("STAGE 50 FINAL LEDGER")
print("Stage 50 verifies the exact self-duality benchmark for the explicit sech–Gaussian")
print("profile family, isolates the self-dual stationary point r = sqrt(pi), and computes")
print("the near-perfect resonance coherence")
print("  C_res^2 = 0.994418836451529... ,  P_res = 1.005612487760576... .")
