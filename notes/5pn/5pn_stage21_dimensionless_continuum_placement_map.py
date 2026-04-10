#!/usr/bin/env python3
"""
5pn_stage21_dimensionless_continuum_placement_map.py

SymPy audit for Moving-Throat PDE Stage 21.
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


banner("STAGE 21 — DIMENSIONLESS CONTINUUM PLACEMENT MAP")

# Microscopic definitions
G, c, c_s = sp.symbols("G c c_s", positive=True, real=True)
a, L = sp.symbols("a L", positive=True, real=True)
K_eta_eff, K_W_eff, K_U = sp.symbols("K_eta_eff K_W_eff K_U", positive=True, real=True)
mu_W, mu_eta = sp.symbols("mu_W mu_eta", positive=True, real=True)
T_w = sp.symbols("T_w", positive=True, real=True)
c_etaU, c_UW, c_etaW = sp.symbols("c_etaU c_UW c_etaW", nonzero=True, real=True)
sigma = sp.symbols("sigma", positive=True, real=True)
pi = sp.pi

# Dimensionless kernel ratios from microscopic data

eps_eta_micro = sp.simplify(c_etaU**2 / (K_U * K_eta_eff))
eps_W_micro = sp.simplify(c_UW**2 * sigma / (K_U * K_W_eff))
rho_micro = sp.simplify(c_UW * c_etaU / (K_U * c_etaW))
Z_W_micro = sp.simplify(c_etaW**2 / (K_eta_eff * K_W_eff))
delta_0_micro = sp.simplify(pi**2 * T_w / (L**2 * K_eta_eff))
Lambda_micro = sp.simplify(27 * pi**2 * G * c_s**5 * K_W_eff / (20 * a**5 * c**5 * mu_W))

print("eps_eta (micro) =", eps_eta_micro)
print("eps_W   (micro) =", eps_W_micro)
print("rho     (micro) =", rho_micro)
print("Z_W     (micro) =", Z_W_micro)
print("delta_0 (micro) =", delta_0_micro)
print("Lambda  (micro) =", Lambda_micro)

# Abstract Stage-21 ratios for derivative algebra
eps_eta, eps_W, rho, Z_W, delta_0, Lambda = sp.symbols(
    "eps_eta eps_W rho Z_W delta_0 Lambda", real=True
)

delta = sp.simplify(delta_0 / (1 - eps_eta))
M_mix = sp.simplify(8 * Z_W * (1 + rho)**2 / (pi**2 * (1 - eps_eta) * (1 - eps_W)))
R_target = sp.simplify(Lambda * (1 - eps_eta) * (1 - eps_W)**2 / (Z_W * (1 + rho)**2))

print("\ndelta    =", delta)
print("M_mix    =", M_mix)
print("R_target =", R_target)

# Useful corollary for beta_0
beta_0 = sp.simplify((mu_W / mu_eta) * (K_eta_eff / K_W_eff) * Z_W * (1 + rho)**2 / (1 - eps_W)**2)
print("beta_0   =", beta_0)

# Exact product relation
product_relation = sp.simplify(R_target * M_mix)
expected_product = sp.simplify(8 * Lambda * (1 - eps_W) / pi**2)
expect_zero("R_target*M_mix - 8 Lambda (1-eps_W)/pi^2", product_relation - expected_product)

# Restore microscopic Lambda and verify the physical product formula
to_micro = {Lambda: Lambda_micro, eps_W: eps_W_micro}
expect_zero(
    "product relation in microscopic variables",
    product_relation.subs(to_micro)
    - sp.simplify(54 * G * c_s**5 * K_W_eff * (1 - eps_W_micro) / (5 * a**5 * c**5 * mu_W)),
)

# Exact one-way derivative tendencies

d_delta_deps_eta = sp.simplify(sp.diff(delta, eps_eta))
d_Mmix_deps_eta = sp.simplify(sp.diff(M_mix, eps_eta))
d_Rtarget_deps_eta = sp.simplify(sp.diff(R_target, eps_eta))

d_Mmix_deps_W = sp.simplify(sp.diff(M_mix, eps_W))
d_Rtarget_deps_W = sp.simplify(sp.diff(R_target, eps_W))

d_Mmix_dZ_W = sp.simplify(sp.diff(M_mix, Z_W))
d_Rtarget_dZ_W = sp.simplify(sp.diff(R_target, Z_W))

d_Mmix_drho = sp.simplify(sp.diff(M_mix, rho))
d_Rtarget_drho = sp.simplify(sp.diff(R_target, rho))

print("\nd delta / d eps_eta    =", d_delta_deps_eta)
print("d M_mix / d eps_eta    =", d_Mmix_deps_eta)
print("d R_target / d eps_eta =", d_Rtarget_deps_eta)
print("\nd M_mix / d eps_W      =", d_Mmix_deps_W)
print("d R_target / d eps_W   =", d_Rtarget_deps_W)
print("\nd M_mix / d Z_W        =", d_Mmix_dZ_W)
print("d R_target / d Z_W     =", d_Rtarget_dZ_W)
print("\nd M_mix / d rho        =", d_Mmix_drho)
print("d R_target / d rho     =", d_Rtarget_drho)

banner("STAGE 21 THEOREM LEDGER")
print("The exact placement coordinates are")
print("  delta    = delta_0 / (1 - eps_eta)")
print("  M_mix    = 8 Z_W (1+rho)^2 / [pi^2 (1-eps_eta)(1-eps_W)]")
print("  R_target = Lambda (1-eps_eta)(1-eps_W)^2 / [Z_W (1+rho)^2]")
print()
print("The exact product law is")
print("  R_target M_mix = 8 Lambda (1-eps_W)/pi^2")
print("                  = 54 G c_s^5 K_W^(eff) (1-eps_W) / (5 a^5 c^5 mu_W)")
print()
print("So Stage 21 separates the continuum kernel into")
print("  - geometry lane:        delta(eps_eta)")
print("  - mixed-stability lane: product set by (eps_W, Lambda)")
print("  - redistribution lane:  (eps_eta, Z_W, rho)")
