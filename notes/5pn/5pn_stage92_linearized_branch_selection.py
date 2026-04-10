#!/usr/bin/env python3
"""
5pn_stage92_linearized_branch_selection.py

Stage 92 audit: linearized branch-selection law near the canonical outgoing branch.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 92 — LINEARIZED BRANCH-SELECTION LAW")

eps = sp.symbols("eps", real=True)
s, b, a0, a5 = sp.symbols("s b a_0 a_5", real=True)

S = 1 + eps * s
beta = 1 + eps * b
Sigma0 = eps * a0
Sigma5 = eps * a5

chi_Q = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
chi_series = sp.series(chi_Q, eps, 0, 2).removeO()
print("chi_Q =")
sp.pprint(chi_series)

lin = sp.simplify(sp.expand(chi_series).coeff(eps, 1))
print("linear coefficient =", lin)

expect_zero("linear coefficient - (5b + a0/3 + 9a5)", lin - (5 * b + a0 / 3 + 9 * a5))
expect_zero("overall scale s cancels", sp.diff(lin, s))

a5_pres = sp.solve(sp.Eq(lin, 0), a5)[0]
print("first-order preservation law a5 =", sp.factor(a5_pres))
expect_zero("a5 preservation formula", a5_pres + 5 * b / 9 + a0 / 27)

banner("STAGE 92 FINAL LEDGER")
print("Linearized outgoing-normalization shift:")
print("  delta chi_Q = eps * (5 b + a_0/3 + 9 a_5) + O(eps^2).")
print("Overall mouth scaling drops out to first order; the minimal isotropic branch-selection data are (b, a_0, a_5).")
