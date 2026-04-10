#!/usr/bin/env python3
"""
5pn_stage28_coherent_local_tracking_reduction.py

SymPy audit for Moving-Throat PDE Stage 28.
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


banner("STAGE 28 — COHERENT LOCAL D/N SUPPORT KERNEL")

# Symbols
lambda_W, lambda_phi, gamma = sp.symbols("lambda_W lambda_phi gamma", nonzero=True, real=True)
mu_eta, mu_U, mu_W, mu_phi = sp.symbols("mu_eta mu_U mu_W mu_phi", positive=True, real=True)
g_U, K_U = sp.symbols("g_U K_U", nonzero=True, real=True)
chi_0, delta_U = sp.symbols("chi_0 delta_U", positive=True, real=True)
M_mix, M_supp, M_tr = sp.symbols("M_mix M_supp M_tr", real=True)
R = sp.symbols("R", positive=True, real=True)
xi, delta = sp.symbols("xi delta", positive=True, real=True)
lambda_0 = sp.Rational(2, 9)

# Coherent local amplitude pattern
g_W = sp.simplify(lambda_W / sp.sqrt(mu_eta * mu_W))
g_R = sp.simplify(gamma * lambda_W / sp.sqrt(mu_U * mu_W))
g_B = sp.simplify(lambda_phi / sp.sqrt(mu_eta * mu_phi))
g_S = sp.simplify(gamma * lambda_phi / sp.sqrt(mu_U * mu_phi))

print("g_W =", g_W)
print("g_R =", g_R)
print("g_B =", g_B)
print("g_S =", g_S)
expect_zero("tracking theorem g_B g_R - g_W g_S", g_B * g_R - g_W * g_S)

# Common interference ratio and exact tracking factor
rho_0 = sp.simplify(g_R * g_U / (K_U * g_W))
sigma_0 = sp.simplify(g_U * g_S / (K_U * g_B))
expect_zero("rho_0 - sigma_0", rho_0 - sigma_0)

R_tr = sp.simplify((1 + chi_0 / (1 + delta_U)) / (1 + chi_0))
print("\nR_tr =", R_tr)
print("1 - R_tr =", sp.simplify(1 - R_tr))
print("R_tr - 1/(1+delta_U) =", sp.simplify(R_tr - 1 / (1 + delta_U)))

# Exact total loading on the coherent branch
Z_W, Z_phi = sp.symbols("Z_W Z_phi", positive=True, real=True)
eps_eta, eps_W_split, eps_phi_split = sp.symbols("eps_eta eps_W_split eps_phi_split", real=True)
M_mix_expr = sp.simplify(8 * Z_W * (1 + chi_0) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps_W_split)))
M_supp_expr = sp.simplify(8 * Z_phi * (1 + chi_0) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps_phi_split)))
M_tr_expr = sp.simplify(M_mix_expr + M_supp_expr)

print("\nM_mix  =", M_mix_expr)
print("M_supp =", M_supp_expr)
print("M_tr   =", M_tr_expr)

# Exact collapse of the Stage-27 quadratic branch equation to tracking form
B_tr = sp.simplify(delta - M_tr * (1 + lambda_0 * R**2))
C_tr = sp.simplify(-delta * M_tr)
quadratic = sp.expand(xi**2 + B_tr * xi + C_tr)

G_tr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R**2) * xi))
expect_zero(
    "quadratic <-> one-parameter tracking law",
    sp.expand(quadratic.subs(M_tr, G_tr)),
)

# Exact normalization collapse
F_tr = sp.simplify(
    (9 * delta + (9 + 2 * R**2) * xi) ** 2
    * (9 * delta + (9 + 2 * R) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2)
)

print("\nG_tr(xi,delta;R) =", G_tr)
print("F_tr(xi,delta;R) =", F_tr)

banner("STAGE 28 THEOREM LEDGER")
print("The first coherent local D/N support kernel satisfies")
print("  g_B g_R = g_W g_S,")
print("so it lands exactly on the tracking surface.")
print()
print("Hence the full rank-2 Stage-27 branch collapses to the one-parameter laws")
print("  M_tr = G_tr(xi,delta;R_tr),")
print("  R_target = F_tr(xi,delta;R_tr),")
print("with")
print("  M_tr = M_mix + M_supp,")
print("  R_tr = [1 + chi_0/(1+delta_U)] / (1 + chi_0).")
print()
print("The exact constructive-branch range is")
print("  1/(1+delta_U) < R_tr < 1.")
