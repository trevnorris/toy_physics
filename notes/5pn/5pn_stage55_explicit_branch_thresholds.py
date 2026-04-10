#!/usr/bin/env python3
"""
5pn_stage55_explicit_branch_thresholds.py

Stage 55 audit: explicit branch placement map and threshold surfaces.
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

banner("STAGE 55 — EXPLICIT BRANCH PLACEMENT MAP AND THRESHOLD SURFACES")

chi_s, Lambda_ell, Upsilon_w = sp.symbols("chi_s Lambda_ell Upsilon_w", positive=True, real=True)
Pe_req, Delta0, Deltainf = sp.symbols("Pe_req Delta0 Delta_inf", positive=True, real=True)
chi = sp.symbols("chi", positive=True, real=True)
Lam = sp.symbols("Lam", positive=True, real=True)

subbanner("55.1 — Exact explicit-branch placement map")
kappa = sp.simplify(4 * chi_s**2 + sp.Rational(4,5) * Lambda_ell**2)
eta = Lambda_ell
W_wall = Upsilon_w * Lambda_ell**2
print("kappa =")
sp.pprint(kappa)
print("eta =")
sp.pprint(eta)
print("W_wall =")
sp.pprint(W_wall)

Upsilon_fail = sp.simplify(Pe_req / (Lambda_ell**2 * Deltainf))
Upsilon_suff = sp.simplify(Pe_req / (Lambda_ell**2 * Delta0))
print("Upsilon_fail =")
sp.pprint(Upsilon_fail)
print("Upsilon_suff =")
sp.pprint(Upsilon_suff)

subbanner("55.2 — Shell-gradient dominated asymptotics")
alpha_grad = sp.simplify(2 * Lam / sp.sqrt(5))
Delta_inf_grad = sp.simplify(1 / alpha_grad)
Delta_0_grad = sp.simplify(Lam / (alpha_grad**2 * (alpha_grad + Lam)))
U_fail_grad = sp.simplify(Pe_req / (Lam**2 * Delta_inf_grad))
U_suff_grad = sp.simplify(Pe_req / (Lam**2 * Delta_0_grad))
print("Upsilon_fail^(grad) =")
sp.pprint(U_fail_grad)
print("Upsilon_suff^(grad) =")
sp.pprint(U_suff_grad)
expect_zero("U_fail^(grad) - 2 Pe_req/(sqrt(5) Lam)", U_fail_grad - 2 * Pe_req / (sp.sqrt(5) * Lam))
expect_zero(
    "U_suff^(grad) - (4/5)(1+2/sqrt(5)) Pe_req",
    U_suff_grad - sp.Rational(4,5) * (1 + 2 / sp.sqrt(5)) * Pe_req
)

subbanner("55.3 — Compression-dominated asymptotics")
alpha_comp = sp.simplify(2 * chi)
Delta_inf_comp = sp.simplify(1 / alpha_comp)
Delta_0_comp = sp.simplify(Lam / (alpha_comp**2 * (alpha_comp + Lam)))
U_fail_comp = sp.simplify(Pe_req / (Lam**2 * Delta_inf_comp))
U_suff_comp = sp.simplify(Pe_req / (Lam**2 * Delta_0_comp))
print("Upsilon_fail^(comp) =")
sp.pprint(U_fail_comp)
print("Upsilon_suff^(comp) =")
sp.pprint(U_suff_comp)
expect_zero("U_fail^(comp) - 2 Pe_req chi / Lam^2", U_fail_comp - 2 * Pe_req * chi / Lam**2)
expect_zero(
    "U_suff^(comp) - 4 Pe_req chi^2 (Lam + 2 chi) / Lam^3",
    U_suff_comp - 4 * Pe_req * chi**2 * (Lam + 2 * chi) / Lam**3
)

banner("STAGE 55 FINAL LEDGER")
print("Stage 55 writes the first explicit moving-throat support/source theorem directly in the")
print("parent branch variables (chi_s, Lambda_ell, Upsilon_w), with exact fail/succeed surfaces")
print("Upsilon_fail and Upsilon_suff and clean shell-gradient/compression asymptotics.")
