#!/usr/bin/env python3
"""
5pn_stage46_parent_threshold_theorem.py

Stage 46 audit: parent-overlap threshold theorem and exact microscopic success/no-go tests.
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


banner("STAGE 46 — PARENT-OVERLAP THRESHOLD THEOREM")

rho_star, g_phi, m, cs_star, K_X = sp.symbols("rho_star g_phi m c_s_star K_X", positive=True, real=True)
N_ss, N_pp, O_sp = sp.symbols("N_ss N_pp O_sp", positive=True, real=True)
C2 = sp.symbols("C2", nonnegative=True, real=True)
Pe_req, Delta0, Deltainf, kappa = sp.symbols("Pe_req Delta0 Deltainf kappa", positive=True, real=True)
T_X, L = sp.symbols("T_X L", positive=True, real=True)

G_micro = rho_star * g_phi**2 * O_sp**2 / (m * cs_star**2 * K_X * N_ss)
G_fail = Pe_req / (kappa * Deltainf)
G_suff = Pe_req / (kappa * Delta0)

subbanner("46.1 — Exact thresholds on the parent loading amplitude")
g_fail_sq = sp.simplify(m * cs_star**2 * K_X * N_ss * G_fail / (rho_star * O_sp**2))
g_suff_sq = sp.simplify(m * cs_star**2 * K_X * N_ss * G_suff / (rho_star * O_sp**2))
print("g_(phi,fail)^2 =")
sp.pprint(g_fail_sq)
print("g_(phi,suff)^2 =")
sp.pprint(g_suff_sq)

subbanner("46.2 — Exact coherence thresholds")
G_micro_C = rho_star * g_phi**2 * N_pp * C2 / (m * cs_star**2 * K_X)
expect_zero("coherence-form gain", G_micro_C - G_micro.subs(O_sp**2, N_ss * N_pp * C2))
C_fail_sq = sp.simplify(m * cs_star**2 * K_X * G_fail / (rho_star * g_phi**2 * N_pp))
C_suff_sq = sp.simplify(m * cs_star**2 * K_X * G_suff / (rho_star * g_phi**2 * N_pp))
print("C_fail^2 =")
sp.pprint(C_fail_sq)
print("C_suff^2 =")
sp.pprint(C_suff_sq)

subbanner("46.3 — Exact Cauchy no-go theorem")
G_max = sp.simplify(rho_star * g_phi**2 * N_pp / (m * cs_star**2 * K_X))
expect_zero("G_max from perfect alignment", G_micro_C.subs(C2, 1) - G_max)
print("If G_max(g_phi) < G_fail(kappa,eta), then no profile engineering can rescue the branch.")

subbanner("46.4 — Thresholds in terms of Pe_req, Delta_0, Delta_inf")
# Insert kappa = K_X L^2 / T_X and verify K_X cancellation in the prefactor.
g_fail_explicit = sp.simplify(g_fail_sq.subs(K_X, kappa * T_X / L**2))
g_suff_explicit = sp.simplify(g_suff_sq.subs(K_X, kappa * T_X / L**2))
expect_zero(
    "g_fail explicit threshold",
    g_fail_explicit - m * cs_star**2 * T_X * N_ss * Pe_req / (rho_star * L**2 * O_sp**2 * Deltainf),
)
expect_zero(
    "g_suff explicit threshold",
    g_suff_explicit - m * cs_star**2 * T_X * N_ss * Pe_req / (rho_star * L**2 * O_sp**2 * Delta0),
)

subbanner("46.5 — Matched-profile special case")
g_fail_matched = sp.simplify(g_fail_sq.subs({N_ss: 1, N_pp: 1, O_sp: 1}))
g_suff_matched = sp.simplify(g_suff_sq.subs({N_ss: 1, N_pp: 1, O_sp: 1}))
print("Matched/normalized-lane thresholds =")
sp.pprint(g_fail_matched)
sp.pprint(g_suff_matched)

banner("STAGE 46 FINAL LEDGER")
print("Stage 46 turns the support/source theorem gap into exact fail/succeed surfaces for")
print("the parent confinement-loading amplitude g_phi and the profile coherence C_(sigma phi)^2.")
print("It also gives the first exact Cauchy no-go theorem at fixed support channel.")
