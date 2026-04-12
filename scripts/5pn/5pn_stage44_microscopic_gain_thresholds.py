#!/usr/bin/env python3
"""
5pn_stage44_microscopic_gain_thresholds.py

Stage 44 audit: microscopic gain thresholds and the exact operator phase diagram.
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


banner("STAGE 44 — MICROSCOPIC GAIN THRESHOLDS AND THE EXACT OPERATOR PHASE DIAGRAM")

Pe_req = sp.symbols("Pe_req", positive=True, real=True)
kappa, eta = sp.symbols("kappa eta", positive=True, real=True)
alpha = sp.sqrt(kappa)
chi_sigma, Lambda_phi, T_X, K_X, L = sp.symbols("chi_sigma Lambda_phi T_X K_X L", positive=True, real=True)

Delta0 = eta * (sp.cosh(alpha) - 1) / (kappa * (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)))
Deltainf = (sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) / (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))

subbanner("44.1 — Exact microscopic control variable")
Xi_micro = chi_sigma * Lambda_phi**2 * L**2 / T_X
G_micro = chi_sigma * Lambda_phi**2 / K_X
expect_zero("Xi_micro - kappa*G_micro", (Xi_micro - kappa * G_micro).subs(K_X * L**2 / T_X, kappa).subs(K_X, kappa * T_X / L**2))

G_fail = sp.simplify(Pe_req / (kappa * Deltainf))
G_suff = sp.simplify(Pe_req / (kappa * Delta0))
print("G_fail(kappa,eta) =")
sp.pprint(G_fail)
print("G_suff(kappa,eta) =")
sp.pprint(G_suff)

subbanner("44.2 — Threshold surfaces in microscopic variables")
chi_fail = sp.simplify(T_X * Pe_req / (Lambda_phi**2 * L**2 * Deltainf))
chi_suff = sp.simplify(T_X * Pe_req / (Lambda_phi**2 * L**2 * Delta0))
Lam_fail = sp.simplify(T_X * Pe_req / (chi_sigma * L**2 * Deltainf))
Lam_suff = sp.simplify(T_X * Pe_req / (chi_sigma * L**2 * Delta0))
print("chi_sigma fail/succeed thresholds =")
sp.pprint(chi_fail)
sp.pprint(chi_suff)
print("Lambda_phi^2 fail/succeed thresholds =")
sp.pprint(Lam_fail)
sp.pprint(Lam_suff)

subbanner("44.3 — Exact soft-support limit kappa -> 0+")
Delta0_soft = sp.limit(Delta0, kappa, 0, dir='+')
Deltainf_soft = sp.limit(Deltainf, kappa, 0, dir='+')
expect_zero("Delta_0 soft-support limit", Delta0_soft - sp.Rational(1, 2))
expect_zero("Delta_inf soft-support limit", Deltainf_soft - 1)
print("Therefore G_fail ~ Pe_req / kappa and G_suff ~ 2 Pe_req / kappa.")

subbanner("44.4 — Exact highly compliant-mouth limit eta -> +infinity")
Delta0_inf = sp.simplify(sp.limit(Delta0, eta, sp.oo))
Deltainf_inf = sp.simplify(sp.limit(Deltainf, eta, sp.oo))
expect_zero("Delta_0^(inf)", Delta0_inf - (1 - sp.sech(alpha)) / kappa)
expect_zero("Delta_inf^(inf)", Deltainf_inf - sp.tanh(alpha) / alpha)
G_fail_inf = sp.simplify(Pe_req / (kappa * Deltainf_inf))
G_suff_inf = sp.simplify(Pe_req / (kappa * Delta0_inf))
expect_zero("G_fail^(inf)", G_fail_inf - Pe_req / (alpha * sp.tanh(alpha)))
expect_zero("G_suff^(inf)", G_suff_inf - Pe_req / (1 - sp.sech(alpha)))
print("\nLarge-kappa compliant-mouth asymptotics:")
print("  G_fail^(inf) ~ Pe_req / sqrt(kappa)")
print("  G_suff^(inf) ~ Pe_req")

banner("STAGE 44 FINAL LEDGER")
print("Stage 44 reorganizes the theorem gap into a phase diagram for the exact dimensionless")
print("microscopic gain G_micro = chi_sigma Lambda_phi^2 / K_X against the geometry-controlled")
print("thresholds G_fail(kappa,eta) and G_suff(kappa,eta).")
