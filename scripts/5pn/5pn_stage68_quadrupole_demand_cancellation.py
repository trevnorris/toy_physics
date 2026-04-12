#!/usr/bin/env python3
"""
5pn_stage68_quadrupole_demand_cancellation.py

Stage 68 audit: exact cancellation of the outgoing-normalization factors in the selected quadrupole-demand product.
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

banner("STAGE 68 — QUADRUPOLE-DEMAND CANCELLATION")

NQ_t, beta0, alpha_req, alpha_mix = sp.symbols("N_Q^(target) beta_0 alpha_req alpha_mix", positive=True, real=True)
A, kappa0_sq = sp.symbols("A kappa_0^2", positive=True, real=True)
Pi_tr, C_mix = sp.symbols("Pi_tr C_mix", positive=True, real=True)
mhat, s_minus, lam_minus = sp.symbols("mhat_- s_- lambda_-", positive=True, real=True)
eps_blk = sp.symbols("eps_blk", nonnegative=True, real=True)

subbanner("68.1 — Exact product identities")
R_target = sp.simplify(NQ_t * A / (beta0 * kappa0_sq))
G_tr = sp.simplify(8 * alpha_req / (sp.pi**2 * A))
M_mix = sp.simplify(8 * alpha_mix / (sp.pi**2 * A))

Pi_expr = sp.simplify((R_target * G_tr).subs(kappa0_sq, 8 / sp.pi**2))
C_expr = sp.simplify((R_target * M_mix).subs(kappa0_sq, 8 / sp.pi**2))
print("Pi_tr =")
sp.pprint(Pi_expr)
print("C_mix =")
sp.pprint(C_expr)
expect_zero("Pi_tr/C_mix - alpha_req/alpha_mix", sp.simplify(Pi_expr / C_expr) - alpha_req / alpha_mix)

subbanner("68.2 — Spectral form from selected-mode normalization")
NQ_rel = sp.Eq(NQ_t / beta0, mhat**2 * s_minus / lam_minus)
print("Selected-mode normalization relation:")
sp.pprint(NQ_rel)
Pi_spec = sp.simplify(Pi_expr.subs(NQ_t / beta0, mhat**2 * s_minus / lam_minus))
C_spec = sp.simplify(C_expr.subs(NQ_t / beta0, mhat**2 * s_minus / lam_minus))
print("Pi_tr =")
sp.pprint(Pi_spec)
print("C_mix =")
sp.pprint(C_spec)

subbanner("68.3 — Pure loading-ratio support demand")
rho_alpha = sp.symbols("rho_alpha", positive=True, real=True)
zeta_rho = sp.simplify((rho_alpha - 1) / (1 - eps_blk * (2 - rho_alpha)))
print("zeta_req(rho_alpha,eps_blk) =")
sp.pprint(zeta_rho)
expect_zero(
    "unblocked law zeta_req - (rho_alpha - 1)",
    sp.simplify(zeta_rho.subs(eps_blk, 0) - (rho_alpha - 1))
)

banner("STAGE 68 FINAL LEDGER")
print("Once the outgoing quadrupole branch is normalized, the explicit support theorem loses all")
print("separate dependence on (mhat_-, beta_0, s_-, lambda_-).")
print("It depends only on the loading ratio")
print("  rho_alpha = alpha_req / alpha_mix,")
print("with support demand")
print("  zeta_req = (rho_alpha - 1) / [1 - eps_blk (2 - rho_alpha)].")
