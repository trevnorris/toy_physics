#!/usr/bin/env python3
"""
5pn_stage45_parent_action_gain.py

Stage 45 audit: parent-action projection of the microscopic support/source gain.
"""

from __future__ import annotations

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


banner("STAGE 45 — PARENT-ACTION PROJECTION OF THE MICROSCOPIC SUPPORT/SOURCE GAIN")

rho, rho_star, K, m = sp.symbols("rho rho_star K m", positive=True, real=True)
g_phi, K_X, T_X, L = sp.symbols("g_phi K_X T_X L", positive=True, real=True)
N_ss, N_pp, O_sp = sp.symbols("N_ss N_pp O_sp", positive=True, real=True)

subbanner("45.1 — Exact n=5 compressional stiffness")
U = K * rho**5 / 4
h = sp.diff(U, rho)
cs2 = sp.diff(K * rho**5, rho) / m
hprime = sp.diff(h, rho)
expect_zero("h(rho) - 5K rho^4 / 4", h - 5 * K * rho**4 / 4)
expect_zero("h'(rho) - 5K rho^3", hprime - 5 * K * rho**3)
expect_zero("h'(rho) - m c_s^2 / rho", hprime - m * cs2 / rho)

subbanner("45.2 — One-channel projection of the parent matter energy")
Theta_sigma = sp.simplify(hprime.subs(rho, rho_star) * N_ss)
Lambda_phi = g_phi * O_sp
print("Theta_sigma =")
sp.pprint(Theta_sigma)
print("Lambda_phi =")
sp.pprint(Lambda_phi)

chi_sigma_eff = sp.simplify(1 / Theta_sigma)
chi_sigma_formula = sp.simplify(rho_star / (m * cs2.subs(rho, rho_star) * N_ss))
expect_zero("effective susceptibility", chi_sigma_eff - chi_sigma_formula)

subbanner("45.3 — Exact parent formula for G_micro")
G_micro = sp.simplify(chi_sigma_eff * Lambda_phi**2 / K_X)
G_formula = sp.simplify(rho_star * g_phi**2 * O_sp**2 / (m * cs2.subs(rho, rho_star) * K_X * N_ss))
expect_zero("parent gain formula", G_micro - G_formula)
print("G_micro =")
sp.pprint(G_formula)

C2 = sp.symbols("C2", nonnegative=True, real=True)
C2_formula = O_sp**2 / (N_ss * N_pp)
G_factorized = sp.simplify(rho_star * g_phi**2 * N_pp / (m * cs2.subs(rho, rho_star) * K_X) * C2_formula)
expect_zero("coherence-factor decomposition", G_formula - G_factorized)
print("C_(sigma phi)^2 = O_(sigma phi)^2 / [N_(sigma sigma) N_(phi phi)]")

subbanner("45.4 — Exact microscopic fixed-point strength")
Xi_micro = sp.simplify((K_X * L**2 / T_X) * G_formula)
Xi_formula = sp.simplify(rho_star * g_phi**2 * O_sp**2 * L**2 / (m * cs2.subs(rho, rho_star) * T_X * N_ss))
expect_zero("Xi_micro parent formula", Xi_micro - Xi_formula)

banner("STAGE 45 FINAL LEDGER")
print("Stage 45 pushes the support/source gain all the way back to the parent 4D action.")
print("The exact n=5 GNLS compressional susceptibility and one support/source overlap now give")
print("  G_micro = rho_* g_phi^2 O_(sigma phi)^2 / [m c_(s,*)^2 K_X N_(sigma sigma)],")
print("with the coherence-factor decomposition and Xi_micro following immediately.")
