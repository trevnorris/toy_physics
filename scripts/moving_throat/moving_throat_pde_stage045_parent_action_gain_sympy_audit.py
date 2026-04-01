#!/usr/bin/env python3
"""
Stage 45 SymPy audit.

Checks:
1. n=5 EOS identities h'(rho)=m c_s^2/rho.
2. Parent-channel projection formulas Theta_sigma and Lambda_phi.
3. Exact parent formula for G_micro.
4. Coherence-factor decomposition.
5. Xi_micro = kappa G_micro reduction.
"""

from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 80
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 45 — PARENT-ACTION PROJECTION OF THE MICROSCOPIC GAIN")

# EOS identities
rho, K, m = sp.symbols("rho K m", positive=True, real=True)
U = K * rho**5 / 4
h = sp.diff(U, rho)
hprime = sp.diff(h, rho)
cs_sq = (1 / m) * sp.diff(K * rho**5, rho)

print("U(rho) =", sp.simplify(U))
print("h(rho) =", sp.simplify(h))
print("h'(rho) =", sp.simplify(hprime))
print("c_s^2(rho) =", sp.simplify(cs_sq))
expect_zero("h'(rho) - m c_s^2 / rho", hprime - m * cs_sq / rho)

# Projected source/support coefficients
rho_star, cs_star_sq = sp.symbols("rho_star cs_star_sq", positive=True, real=True)
Nss, Npp, Osp = sp.symbols("N_ss N_pp O_sp", positive=True, real=True)
g_phi, KX, TX, L = sp.symbols("g_phi K_X T_X L", positive=True, real=True)

Theta_sigma = (m * cs_star_sq / rho_star) * Nss
Lambda_phi = g_phi * Osp
chi_eff = sp.simplify(1 / Theta_sigma)
G_micro = sp.simplify(chi_eff * Lambda_phi**2 / KX)
G_expected = sp.simplify(rho_star * g_phi**2 * Osp**2 / (m * cs_star_sq * KX * Nss))

print("Theta_sigma =", Theta_sigma)
print("Lambda_phi =", Lambda_phi)
print("chi_sigma^(eff) =", chi_eff)
print("G_micro =", G_micro)
expect_zero("G_micro - expected parent formula", G_micro - G_expected)

# Coherence factor decomposition
C2 = sp.symbols("C_sp_sq", nonnegative=True, real=True)
C2_def = sp.simplify(Osp**2 / (Nss * Npp))
G_coh = sp.simplify(rho_star * g_phi**2 * Npp * C2 / (m * cs_star_sq * KX))
print("C_(sigma phi)^2 definition =", C2_def)
expect_zero(
    "coherence-factor decomposition",
    G_expected.subs(Osp**2, C2 * Nss * Npp) - G_coh,
)

# Xi_micro relation
kappa = sp.symbols("kappa", positive=True, real=True)
Xi_micro = sp.simplify(kappa * G_micro)
Xi_expected = sp.simplify(rho_star * g_phi**2 * Osp**2 * L**2 / (m * cs_star_sq * TX * Nss))
print("Xi_micro =", Xi_micro)
expect_zero(
    "Xi_micro - parent projected formula",
    Xi_micro.subs(kappa, KX * L**2 / TX) - Xi_expected,
)

print("\nAll Stage 45 symbolic checks passed.")
