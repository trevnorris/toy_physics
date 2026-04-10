#!/usr/bin/env python3
"""
5pn_stage60_family1_theta_extraction.py

Stage 60 audit: shell-weighted extraction of Theta_w on the explicit Family-1 wall.
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

banner("STAGE 60 — SHELL-WEIGHTED EXTRACTION OF THETA_w")

mp.mp.dps = 80

alpha_r = mp.mpf(10)

def S(xi: mp.mpf) -> mp.mpf:
    return (1 + mp.tanh(xi)) / 2

def chi_phi(xi: mp.mpf) -> mp.mpf:
    return mp.mpf('0.5') * mp.sech(xi)**2

xi_star = mp.atanh(2 / mp.sqrt(alpha_r) - 1)

def rho_r(xi: mp.mpf) -> mp.mpf:
    val = 1 - alpha_r * S(xi)**2
    return val**mp.mpf('0.25') if val > 0 else mp.mpf('0')

I_f = mp.quad(lambda xx: chi_phi(xx)**2, [-mp.inf, mp.inf])

num_rho = mp.quad(lambda xx: chi_phi(xx)**2 * rho_r(xx), [-mp.inf, xi_star])
num_rho2 = mp.quad(lambda xx: chi_phi(xx)**2 * rho_r(xx)**2, [-mp.inf, xi_star])

rho_avg = num_rho / I_f
rho2_avg = num_rho2 / I_f

Theta_chi_coeff = 25 * rho2_avg
Theta_J_coeff = 25 * rho_avg**2

print("xi_* =", xi_star)
print("I_f =", I_f)
print("<rho_r>_chi =", rho_avg)
print("<rho_r^2>_chi =", rho2_avg)
print("Theta_w^(chi) =", mp.nstr(Theta_chi_coeff, 18), "* lambda_mu^2")
print("Theta_w^(J)   =", mp.nstr(Theta_J_coeff, 18), "* lambda_mu^2")

banner("STAGE 60 FINAL LEDGER")
print("The explicit Family-1 wall-depth datum is now extracted directly from the canonical wall")
print("profile and canonical support weight:")
print("  Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2,")
print("  Theta_w^(J)   ≈ 0.927552032539308 lambda_mu^2.")
