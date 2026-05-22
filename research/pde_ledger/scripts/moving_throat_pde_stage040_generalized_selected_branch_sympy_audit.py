#!/usr/bin/env python3
"""
Moving-throat PDE — Stage 23 SymPy audit.

What this audit verifies
------------------------
1. The selected lower wall branch for a diagonal 2x2 baseline plus a rank-1
   loading vector z has exact closed-form eigenvalue and eigenvector formulas.
2. If the orbital/worldtube source vector s is not collinear with z, the old
   Stage-18 normalization function deforms to an exact two-vector shape
   function F_{q,eta}.
3. Specializing to the split-U continuum of Stage 22 yields a one-parameter
   family F_U(xi,delta;R_U) and G_U(xi,delta;R_U).
4. Setting R_U = 1 recovers the Stage-18/19 functions exactly.
5. The first-order deformation around the flat-U limit is exact.
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
alpha, z0 = sp.symbols("alpha z0", positive=True, real=True)
q = sp.symbols("q", nonzero=True, real=True)          # q := z1/z0
eta = sp.symbols("eta", real=True)                   # eta := (s1/s0) q
xi = sp.symbols("xi", positive=True, real=True)
R_U, eps = sp.symbols("R_U eps", positive=True, real=True)
pi = sp.pi
lam0 = sp.Rational(2, 9)

banner("STAGE 23 — GENERALIZED SELECTED-BRANCH NORMALIZATION WITH SOURCE/LOADING MISMATCH")

subbanner("23.1 — Exact 2x2 selected-branch solve with a generic loading ratio q = z1/z0")

# Baseline wall matrix diag(A0, A1) with A1 = A0 (1 + delta), rank-1 loading alpha z z^T.
lam_minus = A0 * (1 - xi)
alpha_req = sp.simplify(A0 * xi * (delta + xi) / (z0**2 * (delta + xi) + (q * z0)**2 * xi))
alpha_req_q = sp.simplify(alpha_req)
print("alpha_req =", alpha_req_q)

# Lower eigenvector ratio e1/e0.
r = sp.simplify((A0 * xi - alpha_req_q * z0**2) / (alpha_req_q * z0 * (q * z0)))
print("e1/e0 =", r)
expect_zero("e1/e0 closed form", r - xi * q / (delta + xi))

# Explicit eigenvalue/eigenvector residual check.
# Convention (consistent with the sign chosen for r above): the perturbed
# matrix is M - alpha z z^T, where M = diag(A0, A0*(1+delta)) and
# z = (z0, q*z0)^T. The selected-branch claim is that lam_minus = A0*(1-xi)
# is an eigenvalue and (1, r)^T with r = q*xi/(delta+xi) is the eigenvector.
M_base = sp.Matrix([[A0, 0], [0, A0 * (1 + delta)]])
z_vec = sp.Matrix([z0, q * z0])
M_perturbed = M_base - alpha_req_q * (z_vec * z_vec.T)
e_minus = sp.Matrix([1, q * xi / (delta + xi)])
eig_residual = sp.simplify(M_perturbed * e_minus - lam_minus * e_minus)
expect_zero("eigenvector residual row 0", eig_residual[0])
expect_zero("eigenvector residual row 1", eig_residual[1])

subbanner("23.2 — Exact overlap formulas and generalized F,G functions")

z_overlap_sq = sp.simplify((1 + q * r)**2 / (1 + r**2))
s_overlap_sq = sp.simplify((1 + eta * xi / (delta + xi))**2 / (1 + r**2))

print("(z.e_-)^2 / z0^2 =", sp.factor(z_overlap_sq))
print("(s.e_-)^2 / s0^2 =", sp.factor(s_overlap_sq))

F_general = sp.simplify((A0 / lam_minus) * z_overlap_sq * s_overlap_sq)
G_general = sp.simplify((z0**2 / A0) * alpha_req_q)

print("F_(q,eta) =", sp.factor(F_general))
print("G_q =", sp.factor(G_general))

F_expected = sp.simplify((delta + (1 + q**2) * xi)**2 * (delta + (1 + eta) * xi)**2 / ((1 - xi) * ((delta + xi)**2 + q**2 * xi**2)**2))
G_expected = sp.simplify(xi * (delta + xi) / (delta + (1 + q**2) * xi))

expect_zero("F_general - expected", F_general - F_expected)
expect_zero("G_general - expected", G_general - G_expected)

subbanner("23.3 — Specialization to the split-U continuum of Stage 22")

# For the physical split-U case, the source vector is the original v direction,
# while the loading vector obeys z1/z0 = (kappa1/kappa0) R_U = -sqrt(lambda0) R_U.
q_U = sp.simplify(-sp.sqrt(lam0) * R_U)
eta_U = sp.simplify(lam0 * R_U)
F_U = sp.simplify(F_expected.subs({q: q_U, eta: eta_U}))
G_U = sp.simplify(G_expected.subs({q: q_U}))

print("F_U(xi,delta;R_U) =", sp.factor(F_U))
print("G_U(xi,delta;R_U) =", sp.factor(G_U))

# F_stage18 reproduces the Stage-18 closed-form normalization F(xi, delta)
# verified in scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py lines 46-58.
# G_stage19 reproduces the Stage-19 closed-form loading G(xi, delta)
# verified in scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py lines 53-70.
# Keep these literals in sync with the upstream source of truth.
F_stage18 = sp.simplify((9 * delta + 11 * xi)**4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2)**2))
G_stage19 = sp.simplify(9 * xi * (delta + xi) / (9 * delta + 11 * xi))

expect_zero("F_U(R_U=1) - Stage18 F", sp.simplify(F_U.subs(R_U, 1) - F_stage18))
expect_zero("G_U(R_U=1) - Stage19 G", sp.simplify(G_U.subs(R_U, 1) - G_stage19))

subbanner("23.4 — Independent cross-check of first-order deformation about flat-U limit")

# Define R_U-dependent q and eta substitutions parametrized by eps,
# so that R_U = 1 + eps.
q_U_eps = sp.simplify(-sp.sqrt(lam0) * (1 + eps))
eta_U_eps = sp.simplify(lam0 * (1 + eps))

# Path 1 (script's original route): differentiate F_U with respect to R_U.
HF = sp.simplify(sp.diff(F_U.subs(R_U, 1 + eps), eps).subs(eps, 0) / F_stage18)
HG = sp.simplify(sp.diff(G_U.subs(R_U, 1 + eps), eps).subs(eps, 0) / G_stage19)

# Path 2 (independent): substitute the eps-parametrized (q, eta) into
# the general two-vector functions F_general, G_general (section 23.2),
# then expand to first order in eps and read off the coefficient.
F_general_eps = F_general.subs({q: q_U_eps, eta: eta_U_eps})
G_general_eps = G_general.subs({q: q_U_eps})
HF_direct = sp.simplify(sp.diff(F_general_eps, eps).subs(eps, 0) / F_stage18)
HG_direct = sp.simplify(sp.diff(G_general_eps, eps).subs(eps, 0) / G_stage19)

print("H_F (via F_U)        =", sp.factor(HF))
print("H_F (via F_general)  =", sp.factor(HF_direct))
print("H_G (via G_U)        =", sp.factor(HG))
print("H_G (via G_general)  =", sp.factor(HG_direct))

expect_zero("H_F cross-check (F_U vs F_general)", sp.simplify(HF - HF_direct))
expect_zero("H_G cross-check (G_U vs G_general)", sp.simplify(HG - HG_direct))

banner("STAGE 23 THEOREM LEDGER")
print("1. Once source and loading vectors differ, the Stage-18 normalization function deforms")
print("   from the single-vector form F(xi,delta) to the exact two-vector function")
print("      F_(q,eta)(xi,delta).")
print("2. The required baseline loading is still exact and one-dimensional:")
print("      G_q(xi,delta) = xi (xi+delta) / [delta + (1+q^2) xi].")
print("3. For the split-U continuum, the source/loading mismatch collapses to one exact")
print("   parameter R_U through")
print("      q = -(sqrt(2)/3) R_U,  eta = (2/9) R_U.")
print("4. Setting R_U = 1 recovers the Stage-18/19 flat-U branch exactly.")
print("5. Therefore the first non-flat U structure does not destroy the selected-branch")
print("   theorem geometry, but it deforms it into a one-parameter family indexed by R_U.")
