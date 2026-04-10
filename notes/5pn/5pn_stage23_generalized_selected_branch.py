#!/usr/bin/env python3
"""
5pn_stage23_generalized_selected_branch.py

SymPy audit for Moving-Throat PDE Stage 23.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 23 — GENERALIZED SELECTED-BRANCH NORMALIZATION")

# Core symbols
A0, alpha, z0 = sp.symbols("A0 alpha z0", positive=True, real=True)
q, eta, delta, xi, m = sp.symbols("q eta delta xi m", real=True)

# Exact rank-1 loaded 2x2 branch solve
K_loaded = sp.Matrix([
    [A0 - alpha * z0**2, -alpha * z0**2 * q],
    [-alpha * z0**2 * q, A0 * (1 + delta) - alpha * z0**2 * q**2],
])
lambda_minus = A0 * (1 - xi)
char_residual = sp.expand((K_loaded - lambda_minus * sp.eye(2)).det())
alpha_req = sp.simplify(sp.solve(sp.Eq(char_residual, 0), alpha)[0])
alpha_target = sp.simplify(A0 * xi * (xi + delta) / (z0**2 * (delta + (1 + q**2) * xi)))

print("alpha_req    =", alpha_req)
print("alpha_target =", alpha_target)
expect_zero("alpha_req - alpha_target", alpha_req - alpha_target)

# Eigenvector ratio
r = sp.symbols("r", real=True)
eig_ratio = sp.simplify(sp.solve(sp.Eq((K_loaded - lambda_minus * sp.eye(2))[0, 0] + (K_loaded - lambda_minus * sp.eye(2))[0, 1] * r, 0), r)[0].subs(alpha, alpha_req))
print("e_1/e_0 =", eig_ratio)
expect_zero("eigenvector ratio theorem", eig_ratio - q * xi / (delta + xi))

# Exact overlap formulas with distinct source and loading directions
# Choose the unnormalized lower-mode eigenvector e_- proportional to (delta+xi, q xi).
den = (delta + xi) ** 2 + q**2 * xi**2
z_overlap_sq = sp.simplify((delta + (1 + q**2) * xi) ** 2 / den)
s_overlap_sq = sp.simplify((delta + (1 + eta) * xi) ** 2 / den)

print("\n(z.e_-)^2 / z_0^2 =", z_overlap_sq)
print("(s.e_-)^2 / s_0^2 =", s_overlap_sq)

# Stage-23 shape functions
F_qeta = sp.simplify(
    (delta + (1 + q**2) * xi) ** 2 * (delta + (1 + eta) * xi) ** 2
    / ((1 - xi) * ((delta + xi) ** 2 + q**2 * xi**2) ** 2)
)
G_q = sp.simplify(xi * (xi + delta) / (delta + (1 + q**2) * xi))

print("\nF_(q,eta)(xi,delta) =", F_qeta)
print("G_q(xi,delta)       =", G_q)

# Specialization to the split-U continuum
R = sp.symbols("R", real=True)
lambda_0 = sp.Rational(2, 9)
q_U = -sp.sqrt(lambda_0) * R
eta_U = lambda_0 * R

F_U = sp.simplify(F_qeta.subs({q: q_U, eta: eta_U}))
G_U = sp.simplify(G_q.subs(q**2, lambda_0 * R**2))
F_flat = sp.simplify(F_U.subs(R, 1))
G_flat = sp.simplify(G_U.subs(R, 1))

print("\nF_U(xi,delta;R) =", F_U)
print("G_U(xi,delta;R) =", G_U)
print("F_flat          =", F_flat)
print("G_flat          =", G_flat)

# Exact recovery of the flat-U branch at R = 1
expect_zero(
    "F_U(R=1) - Stage-18/19 flat branch",
    F_flat - (9 * delta + 11 * xi) ** 4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2),
)
expect_zero(
    "G_U(R=1) - Stage-18/19 flat branch",
    G_flat - 9 * xi * (xi + delta) / (9 * delta + 11 * xi),
)

# Exact small-deformation expansion about the flat branch
# R = 1 + eps

eps = sp.symbols("eps", real=True)
HF_series = sp.series((F_U / F_flat).subs(R, 1 + eps), eps, 0, 2).removeO()
HG_series = sp.series((G_U / G_flat).subs(R, 1 + eps), eps, 0, 2).removeO()
H_F = sp.simplify((HF_series - 1) / eps)
H_G = sp.simplify((HG_series - 1) / eps)

H_F_target = sp.simplify(
    4 * xi * (27 * delta**2 + 36 * delta * xi + 11 * xi**2)
    / ((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2))
)
H_G_target = sp.simplify(-4 * xi / (9 * delta + 11 * xi))

print("\nH_F =", H_F)
print("H_G =", H_G)
expect_zero("H_F exact formula", H_F - H_F_target)
expect_zero("H_G exact formula", H_G - H_G_target)

banner("STAGE 23 THEOREM LEDGER")
print("The first non-flat U structure does not destroy the selected-branch geometry.")
print("It deforms it from the old one-vector branch to the exact two-vector branch")
print("  F_(q,eta)(xi,delta),   G_q(xi,delta),")
print("with the split-U continuum collapsing to the one-parameter family")
print("  F_U(xi,delta;R_U),   G_U(xi,delta;R_U).")
print()
print("Setting R_U = 1 recovers the old flat branch exactly.")
print("For R_U = 1 + eps, the exact first deformation is")
print("  F_U/F_flat = 1 + eps H_F + O(eps^2),")
print("  G_U/G_flat = 1 + eps H_G + O(eps^2),")
print("with the verified closed forms for H_F and H_G printed above.")
