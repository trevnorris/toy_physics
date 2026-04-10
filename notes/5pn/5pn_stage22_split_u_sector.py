#!/usr/bin/env python3
"""
5pn_stage22_split_u_sector.py

SymPy audit for Moving-Throat PDE Stage 22.
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


banner("STAGE 22 — FIRST NON-FLAT U DOUBLET")

# Symbols
K_eta_eff, K_U, K_W_eff = sp.symbols("K_eta_eff K_U K_W_eff", positive=True, real=True)
mu_eta, mu_W = sp.symbols("mu_eta mu_W", positive=True, real=True)
T_U, T_w, L = sp.symbols("T_U T_w L", positive=True, real=True)
c_etaU, c_UW, c_etaW = sp.symbols("c_etaU c_UW c_etaW", nonzero=True, real=True)
rho_0, g_W, eps_eta, eps_W, delta_0 = sp.symbols("rho_0 g_W eps_eta eps_W delta_0", real=True)
G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
pi = sp.pi

# Exact D/N overlap data
kappa_0 = sp.simplify(2 * sp.sqrt(2) / pi)
kappa_1 = sp.simplify(-4 / (3 * pi))
sigma = sp.simplify(kappa_0**2 + kappa_1**2)
lambda_0 = sp.simplify(kappa_1**2 / kappa_0**2)

print("kappa_0 =", kappa_0)
print("kappa_1 =", kappa_1)
print("sigma   =", sigma)
print("lambda_0=", lambda_0)

# First non-flat U doublet

delta_U = sp.simplify(pi**2 * T_U / (L**2 * K_U))
K_U0 = K_U
K_U1 = sp.simplify(K_U * (1 + delta_U))

print("\ndelta_U =", delta_U)
print("K_(U0)  =", K_U0)
print("K_(U1)  =", K_U1)

# Split direct softening and shifted anisotropy
A_0 = sp.simplify((K_eta_eff - c_etaU**2 / K_U) / mu_eta)
A_1 = sp.simplify((K_eta_eff * (1 + delta_0) - c_etaU**2 / K_U1) / mu_eta)

eps_eta_def = sp.simplify(c_etaU**2 / (K_U * K_eta_eff))
delta_split = sp.simplify((delta_0 + eps_eta * delta_U / (1 + delta_U)) / (1 - eps_eta))

print("\nA_0         =", A_0)
print("A_1         =", A_1)
print("delta_split =", delta_split)
expect_zero(
    "A_1 - A_0 (1 + delta_split)",
    (A_1 - A_0 * (1 + delta_split)).subs(eps_eta, eps_eta_def),
)

# Exact split mixed blocking ratio
S_U = sp.simplify(kappa_0**2 / K_U + kappa_1**2 / K_U1)
eps_W_split_direct = sp.simplify(c_UW**2 * S_U / K_W_eff)
eps_W_direct = sp.simplify(c_UW**2 * sigma / (K_U * K_W_eff))
eps_W_split = sp.simplify(eps_W * (1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)))

print("\nS_U                 =", S_U)
print("eps_W^(split,direct)=", eps_W_split_direct)
print("eps_W^(flat)        =", eps_W_direct)
expect_zero(
    "eps_W_split direct - exact formula",
    eps_W_split_direct.subs(c_UW**2 / K_W_eff, eps_W * K_U / sigma) - eps_W_split,
)

# Mixed loading vector and direction splitting theorem
z_0 = sp.simplify(kappa_0 * g_W * (1 + rho_0))
z_1 = sp.simplify(kappa_1 * g_W * (1 + rho_0 / (1 + delta_U)))
R_U = sp.simplify((1 + rho_0 / (1 + delta_U)) / (1 + rho_0))
D_dir = sp.simplify(kappa_0 * z_1 - kappa_1 * z_0)

print("\nz_0  =", z_0)
print("z_1  =", z_1)
print("R_U  =", R_U)
print("D_dir=", D_dir)
expect_zero("z_1/z_0 - (kappa_1/kappa_0) R_U", sp.simplify(z_1 / z_0 - (kappa_1 / kappa_0) * R_U))
expect_zero(
    "D_dir factorization",
    D_dir + sp.simplify(kappa_0 * kappa_1 * g_W * rho_0 * delta_U / (1 + delta_U)),
)

# Exact split continuum placement map
Z_W, Lambda = sp.symbols("Z_W Lambda", positive=True, real=True)
eps_W_split_sym = sp.symbols("eps_W_split", real=True)

delta_split_sym = sp.symbols("delta_split", real=True)
M_mix_split = sp.simplify(8 * Z_W * (1 + rho_0)**2 / (pi**2 * (1 - eps_eta) * (1 - eps_W_split_sym)))
R_target_split = sp.simplify(Lambda * (1 - eps_eta) * (1 - eps_W_split_sym)**2 / (Z_W * (1 + rho_0)**2))
expect_zero(
    "split-U product relation",
    sp.simplify(M_mix_split * R_target_split - 8 * Lambda * (1 - eps_W_split_sym) / pi**2),
)

# Small-delta_U expansions
small_R_U = sp.series(R_U, delta_U, 0, 2).removeO()
small_eps_W_split = sp.series((1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)), delta_U, 0, 2).removeO()

print("\nR_U expansion at small delta_U        =", small_R_U)
print("eps_W^(split)/eps_W expansion        =", small_eps_W_split)
print("delta_split expansion at small delta_U =", sp.series(delta_split, delta_U, 0, 2).removeO())

banner("STAGE 22 THEOREM LEDGER")
print("Turning on the first axial U structure preserves the scalar placement map but breaks")
print("the old source/loading collinearity except on the exact surfaces")
print("  delta_U = 0  or  rho_0 = 0.")
print()
print("The exact new objects are")
print("  delta_split   = [delta_0 + eps_eta delta_U/(1+delta_U)] / (1-eps_eta)")
print("  eps_W_split   = eps_W [1 - (2/11) delta_U/(1+delta_U)]")
print("  R_U           = [1 + rho_0/(1+delta_U)] / (1 + rho_0)")
print("  D_dir         = - kappa_0 kappa_1 g_W rho_0 delta_U / (1+delta_U)")
print()
print("So Stage 22 preserves the scalar continuum placement law but introduces the exact")
print("direction-splitting invariant that Stage 23 must resolve in the selected-branch normalization.")
