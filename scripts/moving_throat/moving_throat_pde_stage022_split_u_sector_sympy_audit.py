#!/usr/bin/env python3
"""
Moving-throat PDE — Stage 22 SymPy audit.

What this audit verifies
------------------------
1. Turning on the first axial stiffness in the internal U sector splits the flat
   U-doublet exactly.
2. The direct wall softening becomes mode-dependent and yields an exact shifted
   anisotropy ratio delta_split.
3. The mixed U/W blocking ratio becomes an exact split effective quantity
   eps_W_split.
4. The mixed loading vector z is no longer generically collinear with the source
   vector v, and the exact direction-splitting invariant factors cleanly.
5. The Stage-21 continuum placement map survives at the scalar level with
   eps_W -> eps_W_split and delta -> delta_split.
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


pi = sp.pi

# Basic positive symbols.
K_eta_eff, K_U, K_W_eff = sp.symbols("K_eta_eff K_U K_W_eff", positive=True, real=True)
mu_eta, mu_W = sp.symbols("mu_eta mu_W", positive=True, real=True)
c_etaU, c_etaW, c_UW = sp.symbols("c_etaU c_etaW c_UW", positive=True, real=True)
Tw, L = sp.symbols("T_w L", positive=True, real=True)
G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)

delta0, deltaU = sp.symbols("delta0 deltaU", positive=True, real=True)
eps_eta, eps_W, rho0, Z_W, Lambda = sp.symbols("eps_eta eps_W rho0 Z_W Lambda", positive=True, real=True)

kappa0 = 2 * sp.sqrt(2) / pi
kappa1 = -4 / (3 * pi)
sigma = sp.simplify(kappa0**2 + kappa1**2)
lam0 = sp.simplify(kappa1**2 / kappa0**2)

banner("STAGE 22 — SPLIT-U SECTOR: EXACT CONTINUUM REDUCTION")
print("kappa0 =", sp.simplify(kappa0))
print("kappa1 =", sp.simplify(kappa1))
print("sigma  =", sigma)
print("lambda0 = kappa1^2/kappa0^2 =", lam0)

subbanner("22.1 — Exact U-mode split and direct wall softening")

# First nontrivial U-sector axial structure.
K_U1 = sp.simplify(K_U * (1 + deltaU))
A0 = sp.simplify((K_eta_eff - c_etaU**2 / K_U) / mu_eta)
A1 = sp.simplify((K_eta_eff * (1 + delta0) - c_etaU**2 / K_U1) / mu_eta)

eps_eta_def = sp.simplify(c_etaU**2 / (K_U * K_eta_eff))
delta0_def = sp.symbols("delta0_def", positive=True)  # placeholder not used algebraically

delta_split = sp.simplify((delta0 + eps_eta * deltaU / (1 + deltaU)) / (1 - eps_eta))
A0_expected = sp.simplify(K_eta_eff * (1 - eps_eta) / mu_eta)
A1_expected = sp.simplify(A0_expected * (1 + delta_split))

print("A0 =", A0_expected)
print("A1 =", sp.expand(A1_expected))
print("delta_split =", delta_split)

expect_zero("A0 direct - expected", A0.subs({c_etaU**2: eps_eta * K_U * K_eta_eff}) - A0_expected)
expect_zero("A1 direct - expected", A1.subs({c_etaU**2: eps_eta * K_U * K_eta_eff}) - A1_expected)

subbanner("22.2 — Exact mixed blocking ratio with split U sector")

# Split-U mixed blocking.
S_U = sp.simplify(kappa0**2 / K_U + kappa1**2 / K_U1)
eps_W_direct = sp.simplify(c_UW**2 * S_U / K_W_eff)
eps_W_split = sp.simplify(eps_W * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))

print("S_U =", S_U)
print("eps_W_split =", eps_W_split)
expect_zero(
    "eps_W direct - split formula",
    eps_W_direct.subs({c_UW**2: eps_W * K_U * K_W_eff / sigma}) - eps_W_split,
)

subbanner("22.3 — Mixed loading vector and exact direction-splitting invariant")

g_W = sp.simplify(c_etaW / sp.sqrt(mu_eta * mu_W))

z0 = sp.simplify(kappa0 * g_W * (1 + rho0))
z1 = sp.simplify(kappa1 * g_W * (1 + rho0 / (1 + deltaU)))
R_U = sp.simplify((1 + rho0 / (1 + deltaU)) / (1 + rho0))

print("z0 =", z0)
print("z1 =", z1)
print("R_U =", R_U)
expect_zero("z1/z0 - (kappa1/kappa0) R_U", sp.simplify(z1 / z0 - (kappa1 / kappa0) * R_U))

D_dir = sp.simplify(kappa0 * z1 - kappa1 * z0)
D_dir_expected = sp.simplify(-kappa0 * kappa1 * g_W * rho0 * deltaU / (1 + deltaU))
print("D_dir =", D_dir)
expect_zero("direction-splitting invariant", D_dir - D_dir_expected)

print("Collinearity theorem: D_dir = 0 iff deltaU = 0 or rho0 = 0 (equivalently g_R g_U = 0).")

subbanner("22.4 — Split-U continuum placement map")

# Stage-21 dimensionless kernel ratios, now with split effective blocking.
M_mix_split = sp.simplify(8 * Z_W * (1 + rho0)**2 / (pi**2 * (1 - eps_eta) * (1 - eps_W_split)))
R_target_split = sp.simplify(Lambda * (1 - eps_eta) * (1 - eps_W_split)**2 / (Z_W * (1 + rho0)**2))
product = sp.simplify(M_mix_split * R_target_split)

print("M_mix^(split U) =", M_mix_split)
print("R_target^(split U) =", R_target_split)
print("product =", product)
expect_zero("product law", product - 8 * Lambda * (1 - eps_W_split) / pi**2)

subbanner("22.5 — Small-splitting expansions")

delta_split_series = sp.simplify(sp.series(delta_split, deltaU, 0, 2).removeO())
eps_W_series = sp.simplify(sp.series(eps_W_split, deltaU, 0, 2).removeO())
R_U_series = sp.simplify(sp.series(R_U, deltaU, 0, 2).removeO())
M0 = sp.simplify(M_mix_split.subs(deltaU, 0))
R0 = sp.simplify(R_target_split.subs(deltaU, 0))
M_ratio = sp.simplify(sp.series(M_mix_split / M0, deltaU, 0, 2).removeO())
R_ratio = sp.simplify(sp.series(R_target_split / R0, deltaU, 0, 2).removeO())

print("delta_split =", delta_split_series, "+ O(deltaU^2)")
print("eps_W_split =", eps_W_series, "+ O(deltaU^2)")
print("R_U =", R_U_series, "+ O(deltaU^2)")
print("M_mix_split / M_mix_flat =", M_ratio, "+ O(deltaU^2)")
print("R_target_split / R_target_flat =", R_ratio, "+ O(deltaU^2)")

banner("STAGE 22 THEOREM LEDGER")
print("1. Turning on the first axial U stiffness keeps the scalar placement map exact, but")
print("   renormalizes the direct wall anisotropy and the mixed blocking ratio.")
print("2. The exact shifted anisotropy ratio is")
print("      delta_split = [delta0 + eps_eta deltaU/(1+deltaU)]/(1-eps_eta).")
print("3. The exact mixed blocking ratio becomes")
print("      eps_W_split = eps_W [1 - 2 deltaU/(11(1+deltaU))].")
print("4. The mixed loading vector rotates away from the source/support direction by the")
print("   exact invariant D_dir = -kappa0*kappa1*g_W*rho0*deltaU/(1+deltaU).")
print("5. Therefore the Stage-21 factorization survives at the scalar-placement level, but")
print("   source/loading collinearity is no longer generic once deltaU != 0 and rho0 != 0.")
