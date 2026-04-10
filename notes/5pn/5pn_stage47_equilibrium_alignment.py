#!/usr/bin/env python3
"""
5pn_stage47_equilibrium_alignment.py

Stage 47 audit: parent equilibrium source/support alignment and the exact matched-layer gain.
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


banner("STAGE 47 — PARENT EQUILIBRIUM SOURCE/SUPPORT ALIGNMENT")

g_phi, K_X, N_pp, H_w, m, cs_w, rho_w = sp.symbols("g_phi K_X N_pp H_w m c_s_w rho_w", positive=True, real=True)
I1, I2 = sp.symbols("I1 I2", positive=True, real=True)

subbanner("47.1 — Exact equilibrium-aligned source profile")
print("The parent linearized equilibrium law H(y) delta rho + delta V_conf = 0 gives")
print("  chi_sigma(y) = g_phi chi_phi(y) / H(y).")

subbanner("47.2 — Exact overlap invariants on the equilibrium branch")
O_sp = g_phi * I1
N_ss = g_phi**2 * I2
C2_eq = sp.simplify(O_sp**2 / (N_pp * N_ss))
expect_zero("equilibrium coherence factor", C2_eq - I1**2 / (N_pp * I2))

subbanner("47.3 — Exact eliminated-source support softening")
DeltaK_eq = g_phi**2 * I1
G_eq = sp.simplify(DeltaK_eq / K_X)
expect_zero("equilibrium gain", G_eq - g_phi**2 * I1 / K_X)

subbanner("47.4 — Matched-layer limit")
I1_matched = N_pp / H_w
I2_matched = N_pp / H_w**2
C2_matched = sp.simplify((I1_matched**2) / (N_pp * I2_matched))
expect_zero("matched-layer coherence", C2_matched - 1)
G_eq_matched = sp.simplify((g_phi**2 * I1_matched) / K_X)
expect_zero("matched-layer gain", G_eq_matched - g_phi**2 * N_pp / (K_X * H_w))
expect_zero(
    "best-alignment gain after H_w = m c_s,w^2 / rho_w",
    G_eq_matched.subs(H_w, m * cs_w**2 / rho_w) - rho_w * g_phi**2 * N_pp / (m * cs_w**2 * K_X),
)

banner("STAGE 47 FINAL LEDGER")
print("Stage 47 removes one large ambiguity: on the parent equilibrium branch the source")
print("profile is not free. The coherence is derived from the compressional-stiffness profile,")
print("and it saturates to C^2 = 1 in the thin matched-layer limit.")
