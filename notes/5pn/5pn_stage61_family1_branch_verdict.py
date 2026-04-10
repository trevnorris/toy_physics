#!/usr/bin/env python3
"""
5pn_stage61_family1_branch_verdict.py

Stage 61 audit: explicit Family-1 branch comparison and closing verdict for the subprogram.
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

banner("STAGE 61 — EXPLICIT FAMILY-1 BRANCH COMPARISON")

mp.mp.dps = 80

Theta_chi_coeff = mp.mpf("4.06863235008162")
Theta_J_coeff = mp.mpf("0.927552032539308")
Theta_fail_coeff = mp.mpf("3.62605617972939e-4")
Theta_suff_coeff = mp.mpf("4.21495341569977e-2")

Pe_suff_chi = Theta_chi_coeff / Theta_suff_coeff
Pe_fail_chi = Theta_chi_coeff / Theta_fail_coeff
Pe_suff_J = Theta_J_coeff / Theta_suff_coeff
Pe_fail_J = Theta_J_coeff / Theta_fail_coeff

print("Pe_suff^(chi) =", mp.nstr(Pe_suff_chi, 18), "* lambda_mu^2")
print("Pe_fail^(chi) =", mp.nstr(Pe_fail_chi, 18), "* lambda_mu^2")
print("Pe_suff^(J)   =", mp.nstr(Pe_suff_J, 18), "* lambda_mu^2")
print("Pe_fail^(J)   =", mp.nstr(Pe_fail_J, 18), "* lambda_mu^2")

banner("STAGE 61 FINAL LEDGER")
print("On the explicit Family-1 branch, the wall-depth supply is not the natural bottleneck.")
print("Using the shell-weighted datum:")
print("  Pe_req <= 96.5285247264386 lambda_mu^2  -> guaranteed success,")
print("  Pe_req >= 11220.5441626259 lambda_mu^2 -> guaranteed failure.")
print("Even the conservative floor gives")
print("  Pe_req <= 22.0062226330754 lambda_mu^2 -> guaranteed success,")
print("  Pe_req >= 2558.01892349205 lambda_mu^2 -> guaranteed failure.")
