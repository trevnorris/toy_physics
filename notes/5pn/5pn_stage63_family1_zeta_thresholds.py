#!/usr/bin/env python3
"""
5pn_stage63_family1_zeta_thresholds.py

Stage 63 audit: explicit Family-1 conversion of the Pe_req window into zeta_req thresholds.
"""

from __future__ import annotations

import math
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


def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 63 — FAMILY-1 zeta THRESHOLDS")

mp.mp.dps = 80

kappa_F1 = mp.mpf(12321) / 5
eta_F1 = mp.mpf(37)
f = lambda y: y * mp.tan(y) - eta_F1
y_F1 = mp.findroot(f, mp.mpf("1.53"))
A_F1 = (kappa_F1 + mp.pi**2/4) / (kappa_F1 + y_F1**2)

def Omega(Pe: mp.mpf) -> mp.mpf:
    if abs(Pe) < mp.mpf("1e-20"):
        return mp.mpf(1)
    return mp.pi * Pe * (2 * Pe * mp.e**Pe + mp.pi) / ((4 * Pe**2 + mp.pi**2) * (mp.e**Pe - 1))

def zeta_F1(Pe: mp.mpf) -> mp.mpf:
    return A_F1 * Omega(Pe)**2

Pe_suff_chi = mp.mpf("96.5285247264386")
Pe_fail_chi = mp.mpf("11220.5441626259")
Pe_suff_J = mp.mpf("22.0062226330754")
Pe_fail_J = mp.mpf("2558.01892349205")

zeta_suff_chi = zeta_F1(Pe_suff_chi)
zeta_fail_chi = zeta_F1(Pe_fail_chi)
zeta_suff_J = zeta_F1(Pe_suff_J)
zeta_fail_J = zeta_F1(Pe_fail_J)
zeta_max = A_F1 * mp.pi**2 / 4

print("zeta_suff^(chi)(1) =", zeta_suff_chi)
print("zeta_fail^(chi)(1) =", zeta_fail_chi)
print("zeta_suff^(J)(1)   =", zeta_suff_J)
print("zeta_fail^(J)(1)   =", zeta_fail_J)
print("zeta_max^(F1)      =", zeta_max)

banner("STAGE 63 FINAL LEDGER")
print("The Stage-61 Family-1 wall-depth verdict is now written directly in the quadrupole-demand")
print("variable zeta_req. At lambda_mu = 1:")
print("  zeta_req <= 2.46622291347846  -> guaranteed success (natural shell-weighted branch),")
print("  zeta_req >= 2.46752913273870  -> guaranteed failure,")
print("with the hard explicit ceiling")
print("  zeta_req <  2.46752922945601.")
