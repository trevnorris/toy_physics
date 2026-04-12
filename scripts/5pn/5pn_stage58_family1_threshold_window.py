#!/usr/bin/env python3
"""
5pn_stage58_family1_threshold_window.py

Stage 58 audit: explicit Family-1 threshold window and the last remaining wall-amplitude datum.
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

banner("STAGE 58 — EXPLICIT FAMILY-1 THRESHOLD WINDOW")

mp.mp.dps = 80

Lambda_ell = mp.mpf(37)
eta = mp.mpf(37)
kappa = mp.mpf(12321) / 5
alpha = mp.sqrt(kappa)

f = lambda y: y * mp.tan(y) - eta
y = mp.findroot(f, mp.mpf("1.53"))

Delta0 = eta * (mp.cosh(alpha) - 1) / (alpha**2 * (alpha * mp.sinh(alpha) + eta * mp.cosh(alpha)))
Deltainf = (mp.cosh(alpha) + (eta / alpha) * mp.sinh(alpha) - 1) / (alpha * mp.sinh(alpha) + eta * mp.cosh(alpha))

Upsilon_fail_coeff = 1 / (Lambda_ell**2 * Deltainf)
Upsilon_suff_coeff = 1 / (Lambda_ell**2 * Delta0)
Xi_fail_coeff = 1 / Deltainf
Xi_suff_coeff = 1 / Delta0

print("y_F1 =", y)
print("Delta_0(kappa_F1,eta_F1)   =", Delta0)
print("Delta_inf(kappa_F1,eta_F1) =", Deltainf)
print("Upsilon_fail =", mp.nstr(Upsilon_fail_coeff, 18), "* Pe_req")
print("Upsilon_suff =", mp.nstr(Upsilon_suff_coeff, 18), "* Pe_req")
print("Xi_fail =", mp.nstr(Xi_fail_coeff, 18), "* Pe_req")
print("Xi_suff =", mp.nstr(Xi_suff_coeff, 18), "* Pe_req")

alpha_r = mp.mpf(10)
Theta_fail_coeff = Upsilon_fail_coeff / alpha_r**2
Theta_suff_coeff = Upsilon_suff_coeff / alpha_r**2
print("Theta_fail =", mp.nstr(Theta_fail_coeff, 18), "* Pe_req")
print("Theta_suff =", mp.nstr(Theta_suff_coeff, 18), "* Pe_req")

banner("STAGE 58 FINAL LEDGER")
print("The explicit Family-1/healing-locked branch fixes the operator thresholds numerically:")
print("  Upsilon_fail ≈ 0.0362605617972939 Pe_req,")
print("  Upsilon_suff ≈ 4.21495341569977 Pe_req,")
print("or, after Upsilon_w = 100 Theta_w,")
print("  Theta_fail ≈ 3.62605617972939e-4 Pe_req,")
print("  Theta_suff ≈ 4.21495341569977e-2 Pe_req.")
