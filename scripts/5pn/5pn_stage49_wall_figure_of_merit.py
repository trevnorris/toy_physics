#!/usr/bin/env python3
"""
5pn_stage49_wall_figure_of_merit.py

Stage 49 audit: dimensionless wall figure of merit for the first explicit parent branch.
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


banner("STAGE 49 — DIMENSIONLESS WALL FIGURE OF MERIT")

a, L, V0, T_X, ell, J1, I_f, H_w = sp.symbols("a L V0 T_X ell J1 I_f H_w", positive=True, real=True)
K_X, kappa = sp.symbols("K_X kappa", positive=True, real=True)
Pe_req, Delta0, Deltainf = sp.symbols("Pe_req Delta0 Deltainf", positive=True, real=True)

subbanner("49.1 — Exact wall figure of merit")
G_eq_tw = 4 * sp.pi * a**2 * V0**2 * J1 / (K_X * ell)
W_wall = 4 * sp.pi * a**2 * L**2 * J1 * V0**2 / (T_X * ell)
expect_zero("W_wall - kappa G_eq^(tw)", W_wall - (kappa * G_eq_tw).subs(K_X, kappa * T_X / L**2))

subbanner("49.2 — Exact fail/succeed thresholds in wall form")
W_fail = Pe_req / Deltainf
W_suff = Pe_req / Delta0
print("W_fail =")
sp.pprint(W_fail)
print("W_suff =")
sp.pprint(W_suff)
print("\nSo the thin-wall matched branch obeys")
print("  W_wall <= W_fail  -> fail")
print("  W_wall >= W_suff  -> succeed")

subbanner("49.3 — Monotonicity")
print("W_wall = 4 pi a^2 L^2 J1 V0^2 / (T_X ell)")
print("It is monotone increasing in V0, a, L, and J1, and monotone decreasing in T_X and ell.")

subbanner("49.4 — Constant-compressibility wall form")
W_H = sp.simplify(W_wall.subs(J1, I_f / H_w))
print("W_H =")
sp.pprint(W_H)

banner("STAGE 49 FINAL LEDGER")
print("Stage 49 compresses the first explicit parent branch to one number:")
print("  W_wall = 4 pi a^2 L^2 J1 V0^2 / (T_X ell).")
print("The support/source question is now whether the real branch lands below, within, or")
print("above the exact window [Pe_req/Delta_inf, Pe_req/Delta_0].")
