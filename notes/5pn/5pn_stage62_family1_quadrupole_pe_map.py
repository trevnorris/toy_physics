#!/usr/bin/env python3
"""
5pn_stage62_family1_quadrupole_pe_map.py

Stage 62 audit: exact Family-1 map from quadrupole demand zeta_req to Pe_req.
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

banner("STAGE 62 — FAMILY-1 MAP FROM zeta_req TO Pe_req")

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

zeta_max = A_F1 * mp.pi**2 / 4

print("y_F1 =", y_F1)
print("A_F1 =", A_F1)
print("zeta_F1(0) =", A_F1)
print("zeta_max^(F1) =", zeta_max)

subbanner("62.1 — Small-demand expansion")
eps_z = sp.symbols("eps_z", real=True)
Pe = sp.symbols("Pe", real=True)
coeff = sp.simplify((4 - sp.pi) / sp.pi)
print("Omega_Pe^2 = 1 + ((4-pi)/pi) Pe + O(Pe^2)")
print("So for zeta_req = A_F1 (1 + eps_z),")
print("  Pe_req ~= (pi/(4-pi)) eps_z.")
print("pi/(4-pi) =", mp.pi / (4 - mp.pi))

banner("STAGE 62 FINAL LEDGER")
print("The explicit Family-1 constructive branch spans the exact interval")
print("  A_F1 <= zeta_req <= zeta_max^(F1)")
print("with")
print("  A_F1 ≈ 1.00005192880220,")
print("  zeta_max^(F1) ≈ 2.46752922945601.")
print("So the remaining demand-side question is whether the final selected quadrupole branch")
print("asks for a support ratio below this explicit ceiling.")
