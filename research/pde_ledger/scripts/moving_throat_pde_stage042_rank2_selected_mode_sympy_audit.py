#!/usr/bin/env python3
"""
Moving-throat PDE — Stage 25 SymPy audit.

What this audit verifies
------------------------
1. After inserting the exact Stage-24 support loading, the selected lower-mode
   eigenvector ratio has a closed-form rank-2 expression.
2. The outgoing/mixed and source overlaps are exact and lead to a generalized
   normalization function F_(q,r,t)(xi,delta;m).
3. If the support tracks the mixed direction, the exact rank-2 normalization law
   collapses back to the Stage-23 two-vector law.
4. The source-tied split-U specialization yields the exact function F_src.
5. Setting R_U = 1 recovers the flat-U branch exactly.
6. The first-order source-tied deformation about R_U = 1 is exact.
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
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


A0, delta = sp.symbols("A0 delta", positive=True, real=True)
xi = sp.symbols("xi", positive=True, real=True)
m = sp.symbols("m", real=True)
q, r, t = sp.symbols("q r t", real=True)
lam0 = sp.Rational(2, 9)
R_U, eps = sp.symbols("R_U eps", positive=True, real=True)

banner("STAGE 25 — SELECTED-MODE NORMALIZATION UNDER RANK-2 SUPPORT COMPLETION")

subbanner("25.1 — Exact selected-mode eigenvector ratio after inserting the Stage-24 support loading")

n_req = sp.simplify(
    (xi * (delta + xi) - m * (delta + (1 + q**2) * xi))
    / (delta + (1 + r**2) * xi - m * (q - r)**2)
)

# First-row and second-row versions of the eigenvector ratio e1/e0.
ratio_row1 = sp.simplify((xi - m - n_req) / (m * q + n_req * r))
ratio_row2 = sp.simplify((m * q + n_req * r) / (delta + xi - m * q**2 - n_req * r**2))
ratio_expected = sp.simplify((m * (q - r) + r * xi) / (delta + xi - m * q * (q - r)))

print("e1/e0 from row 1 =")
sp.pprint(sp.factor(ratio_row1))
print("e1/e0 from row 2 =")
sp.pprint(sp.factor(ratio_row2))
print("expected e1/e0 =")
sp.pprint(sp.factor(ratio_expected))

expect_zero("row1 - expected", ratio_row1 - ratio_expected)
expect_zero("row2 - expected", ratio_row2 - ratio_expected)

subbanner("25.2 — Exact overlap formulas and generalized rank-2 normalization function")

D_qr = sp.simplify((delta + xi - m * q * (q - r))**2 + (m * (q - r) + r * xi)**2)
Z_overlap = sp.simplify((1 + q * ratio_expected)**2 / (1 + ratio_expected**2))
S_overlap = sp.simplify((1 + t * ratio_expected)**2 / (1 + ratio_expected**2))

Z_expected = sp.simplify((delta + (1 + q * r) * xi)**2 / D_qr)
S_expected = sp.simplify((delta + (1 + r * t) * xi - m * (q - r) * (q - t))**2 / D_qr)

print("(z.e_-)^2 / z0^2 =")
sp.pprint(sp.factor(Z_overlap))
print("(s.e_-)^2 / s0^2 =")
sp.pprint(sp.factor(S_overlap))

expect_zero("Z_overlap - expected", Z_overlap - Z_expected)
expect_zero("S_overlap - expected", S_overlap - S_expected)

F_general = sp.simplify(Z_overlap * S_overlap / (1 - xi))
print("F_(q,r,t) =")
sp.pprint(sp.factor(F_general))

F_expected = sp.simplify(
    (delta + (1 + q * r) * xi)**2
    * (delta + (1 + r * t) * xi - m * (q - r) * (q - t))**2
    / ((1 - xi) * D_qr**2)
)
expect_zero("F_general - expected", F_general - F_expected)

subbanner("25.3 — Tracking-support collapse back to Stage 23")

F_track = sp.simplify(F_expected.subs(r, q))
F_stage23 = sp.simplify(
    (delta + (1 + q**2) * xi)**2 * (delta + (1 + q * t) * xi)**2
    / ((1 - xi) * ((delta + xi)**2 + q**2 * xi**2)**2)
)
print("F_track =")
sp.pprint(sp.factor(F_track))
expect_zero("tracking collapse", F_track - F_stage23)

subbanner("25.4 — Source-tied split-U specialization")

# Source-tied support: r = t, q = t R_U, with t^2 = lam0.
q_tied_qr = lam0 * R_U
a1 = delta + (1 + q_tied_qr) * xi
b1 = delta + (1 + lam0) * xi - m * lam0 * (R_U - 1)**2
D_src = sp.simplify((delta + xi - m * lam0 * R_U * (R_U - 1))**2 + lam0 * (xi + m * (R_U - 1))**2)
F_src = sp.simplify(a1**2 * b1**2 / ((1 - xi) * D_src**2))

print("F_src =")
sp.pprint(sp.factor(F_src))

# Verify by direct substitution into the general formula.
F_src_direct = sp.simplify(
    F_expected.subs({q: sp.sqrt(lam0) * R_U, r: sp.sqrt(lam0), t: sp.sqrt(lam0)})
)
expect_zero("source-tied specialization", sp.simplify(F_src_direct - F_src))

subbanner("25.5 — Exact flat-U recovery")

F_flat = sp.simplify(
    (delta + (1 + lam0) * xi)**4
    / ((1 - xi) * ((delta + xi)**2 + lam0 * xi**2)**2)
)
expect_zero("F_src(R_U=1) - F_flat", sp.simplify(F_src.subs(R_U, 1) - F_flat))

subbanner("25.6 — First-order source-tied deformation about R_U = 1")

n_src = sp.simplify(
    (xi * (delta + xi) - m * (delta + (1 + lam0 * R_U**2) * xi))
    / (delta + (1 + lam0) * xi - m * lam0 * (R_U - 1)**2)
)
G_flat = sp.simplify(xi * (delta + xi) / (delta + (1 + lam0) * xi))

H_n_src = sp.simplify(sp.diff(n_src, R_U).subs(R_U, 1))
H_n_expected = sp.simplify(-2 * lam0 * m * xi / (delta + (1 + lam0) * xi))
print("H_n^(src) =")
sp.pprint(sp.factor(H_n_src))
expect_zero("H_n^(src) - expected", H_n_src - H_n_expected)

F_ratio = sp.simplify(F_src / F_flat)
H_F_src = sp.simplify(sp.diff(F_ratio, R_U).subs(R_U, 1))
H_F_expected = sp.simplify(
    2 * lam0 * (
        xi * ((delta + xi)**2 + lam0 * xi**2)
        + 2 * m * delta * (delta + (1 + lam0) * xi)
    )
    / ((delta + (1 + lam0) * xi) * ((delta + xi)**2 + lam0 * xi**2))
)
print("H_F^(src) =")
sp.pprint(sp.factor(H_F_src))
expect_zero("H_F^(src) - expected", H_F_src - H_F_expected)

# Linearized expansions.
R_series = 1 + eps
n_linear = sp.simplify((G_flat - m) + eps * H_n_expected)
F_linear = sp.simplify(1 + eps * H_F_expected)
expect_zero(
    "n_src - linear expansion",
    sp.expand(sp.series(n_src.subs(R_U, R_series), eps, 0, 2).removeO() - n_linear),
)
expect_zero(
    "F_src/F_flat - linear expansion",
    sp.expand(sp.series(F_ratio.subs(R_U, R_series), eps, 0, 2).removeO() - F_linear),
)

banner("STAGE 25 THEOREM LEDGER")
print("1. With two loading directions, the selected-mode normalization remains exact and is")
print("   governed by the three-direction function F_(q,r,t)(xi,delta;m).")
print("2. If the support tracks the mixed direction, the entire rank-2 normalization law")
print("   collapses exactly to the Stage-23 two-vector function.")
print("3. If the support stays tied to the original source direction, the split-U continuum")
print("   yields the exact source-tied normalization function F_src(xi,delta;m,R_U).")
print("4. Setting R_U = 1 recovers the flat-U branch exactly.")
print("5. On the natural constructive split-U branch with R_U < 1, the source-tied closure")
print("   raises the required support loading and lowers the selected-mode normalization")
print("   relative to the flat-U subtraction picture.")
