#!/usr/bin/env python3
"""
moving_throat_pde_stage21_dimensionless_continuum_placement_sympy_audit.py

Stage 21 SymPy audit:
compress the Stage-20 continuum formulas into the exact dimensionless kernel map,
verify the product relation, and factor the one-way parameter tendencies.
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


def expect_zero(name: str, expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 21 — DIMENSIONLESS CONTINUUM PLACEMENT AUDIT")

# ---------------------------------------------------------------------------
# 1. Stage-20 continuum formulas
# ---------------------------------------------------------------------------

subbanner("1. Stage-20 continuum formulas")

G, c_light, c_s, a, L = sp.symbols("G c_light c_s a L", positive=True, real=True)
mu_eta, mu_W = sp.symbols("mu_eta mu_W", positive=True, real=True)
K_eta, T_Omega, T_w = sp.symbols("K_eta T_Omega T_w", positive=True, real=True)
K_U, K_W, T_W = sp.symbols("K_U K_W T_W", positive=True, real=True)
c_etaU, c_etaW, c_UW = sp.symbols("c_etaU c_etaW c_UW", real=True)

sigma = sp.Rational(88, 9) / sp.pi**2
Keta_eff = sp.simplify(K_eta + 6 * T_Omega)
KWeff = sp.simplify(K_W + sp.pi**2 * T_W / (4 * L**2))

A = sp.simplify((K_U * Keta_eff - c_etaU**2) / (mu_eta * K_U))
delta = sp.simplify(sp.pi**2 * T_w * K_U / (L**2 * (K_U * Keta_eff - c_etaU**2)))
M_mix = sp.simplify(
    8 * (K_U * c_etaW + c_UW * c_etaU) ** 2
    / (sp.pi**2 * (K_U * Keta_eff - c_etaU**2) * (K_U * KWeff - c_UW**2 * sigma))
)
beta0 = sp.simplify(
    (mu_W / mu_eta) * (K_U * c_etaW + c_UW * c_etaU) ** 2
    / ((K_U * KWeff - c_UW**2 * sigma) ** 2)
)
R_target = sp.simplify(
    (54 * G * c_s**5 / (5 * a**5 * c_light**5)) * A / ((8 / sp.pi**2) * beta0)
)

print(f"A         = {A}")
print(f"delta     = {delta}")
print(f"M_mix     = {M_mix}")
print(f"beta0     = {beta0}")
print(f"R_target  = {R_target}")

# ---------------------------------------------------------------------------
# 2. Dimensionless kernel ratios
# ---------------------------------------------------------------------------

subbanner("2. Dimensionless kernel substitutions")

eps_eta, eps_W, rho, Z_W, delta0, Lambda = sp.symbols(
    "eps_eta eps_W rho Z_W delta0 Lambda", positive=True, real=True
)

subs_dimless = {
    c_etaU**2: eps_eta * K_U * Keta_eff,
    c_UW**2: eps_W * K_U * KWeff / sigma,
    c_UW * c_etaU: rho * K_U * c_etaW,
    c_etaW**2: Z_W * Keta_eff * KWeff,
    T_w: delta0 * L**2 * Keta_eff / sp.pi**2,
    G: 20 * Lambda * a**5 * c_light**5 * mu_W / (27 * sp.pi**2 * c_s**5 * KWeff),
}

def apply_dimless(expr):
    expr = sp.expand(expr.subs({c_UW * c_etaU: rho * K_U * c_etaW}))
    expr = sp.expand(expr.subs(c_etaW**2, Z_W * Keta_eff * KWeff))
    expr = sp.expand(expr.subs({
        c_etaU**2: eps_eta * K_U * Keta_eff,
        c_UW**2: eps_W * K_U * KWeff / sigma,
        T_w: delta0 * L**2 * Keta_eff / sp.pi**2,
        G: 20 * Lambda * a**5 * c_light**5 * mu_W / (27 * sp.pi**2 * c_s**5 * KWeff),
    }))
    return sp.simplify(expr)

delta_dimless = apply_dimless(delta)
M_dimless = apply_dimless(M_mix)
R_dimless = apply_dimless(R_target)
beta_dimless = apply_dimless(beta0)

expect_zero("delta - delta0/(1-eps_eta)", delta_dimless - delta0 / (1 - eps_eta))
expect_zero(
    "M_mix - 8 Z_W (1+rho)^2 / [pi^2 (1-eps_eta)(1-eps_W)]",
    M_dimless - 8 * Z_W * (1 + rho) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps_W)),
)
expect_zero(
    "R_target - Lambda (1-eps_eta)(1-eps_W)^2 / [Z_W (1+rho)^2]",
    R_dimless - Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W * (1 + rho) ** 2),
)
expect_zero(
    "beta0 - (mu_W/mu_eta)(Keta_eff/KWeff) Z_W (1+rho)^2/(1-eps_W)^2",
    beta_dimless - (mu_W / mu_eta) * (Keta_eff / KWeff) * Z_W * (1 + rho) ** 2 / (1 - eps_W) ** 2,
)

print(f"delta(dimless)    = {sp.simplify(delta_dimless)}")
print(f"M_mix(dimless)    = {sp.simplify(M_dimless)}")
print(f"R_target(dimless) = {sp.simplify(R_dimless)}")
print(f"beta0(dimless)    = {sp.simplify(beta_dimless)}")

# ---------------------------------------------------------------------------
# 3. Exact product relation
# ---------------------------------------------------------------------------

subbanner("3. Exact product relation")

product = sp.simplify(R_dimless * M_dimless)
expect_zero("R_target M_mix - 8 Lambda (1-eps_W)/pi^2", product - 8 * Lambda * (1 - eps_W) / sp.pi**2)

# Also verify the Stage-20 equivalent form.
OmegaU2 = sp.symbols("OmegaU2", positive=True, real=True)
Delta0 = sp.symbols("Delta0", positive=True, real=True)
NQ = 54 * G * c_s**5 / (5 * a**5 * c_light**5)
expect_zero("8 Lambda (1-eps_W)/pi^2 - NQ * KWeff(1-eps_W)/mu_W",
            8 * Lambda * (1 - eps_W) / sp.pi**2
            - (54 * G * c_s**5 / (5 * a**5 * c_light**5) * KWeff * (1 - eps_W) / mu_W).subs(subs_dimless))

print(f"R_target M_mix = {product}")

# ---------------------------------------------------------------------------
# 4. One-way parameter tendencies
# ---------------------------------------------------------------------------

subbanner("4. Exact derivative factors")

delta_expr = delta0 / (1 - eps_eta)
M_expr = 8 * Z_W * (1 + rho) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps_W))
R_expr = Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W * (1 + rho) ** 2)

d_delta_deps_eta = sp.simplify(sp.diff(delta_expr, eps_eta))
dM_dZ = sp.simplify(sp.diff(M_expr, Z_W))
dR_dZ = sp.simplify(sp.diff(R_expr, Z_W))
dM_deps_eta = sp.simplify(sp.diff(M_expr, eps_eta))
dR_deps_eta = sp.simplify(sp.diff(R_expr, eps_eta))
dM_deps_W = sp.simplify(sp.diff(M_expr, eps_W))
dR_deps_W = sp.simplify(sp.diff(R_expr, eps_W))
dM_drho = sp.simplify(sp.diff(M_expr, rho))
dR_drho = sp.simplify(sp.diff(R_expr, rho))

expect_zero("d delta / d eps_eta factorization", d_delta_deps_eta - delta0 / (1 - eps_eta) ** 2)
expect_zero("d M / d Z_W factorization", dM_dZ - 8 * (1 + rho) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps_W)))
expect_zero("d R / d Z_W factorization", dR_dZ + Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W**2 * (1 + rho) ** 2))
expect_zero("d M / d eps_eta factorization", dM_deps_eta - 8 * Z_W * (1 + rho) ** 2 / (sp.pi**2 * (1 - eps_eta) ** 2 * (1 - eps_W)))
expect_zero("d R / d eps_eta factorization", dR_deps_eta + Lambda * (1 - eps_W) ** 2 / (Z_W * (1 + rho) ** 2))
expect_zero("d M / d eps_W factorization", dM_deps_W - 8 * Z_W * (1 + rho) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps_W) ** 2))
expect_zero("d R / d eps_W factorization", dR_deps_W + 2 * Lambda * (1 - eps_eta) * (1 - eps_W) / (Z_W * (1 + rho) ** 2))
expect_zero("d M / d rho factorization", dM_drho - 16 * Z_W * (1 + rho) / (sp.pi**2 * (1 - eps_eta) * (1 - eps_W)))
expect_zero("d R / d rho factorization", dR_drho + 2 * Lambda * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W * (1 + rho) ** 3))

# Sign assertions under the natural transfer-branch assumption
# (0 < eps_eta < 1, 0 < eps_W < 1, 1 + rho > 0, Z_W > 0, Lambda > 0, delta0 > 0)
# Each derivative factors as (sign) * (manifestly positive template); we assert the sign.

# delta = delta0/(1 - eps_eta), d delta / d eps_eta = + delta0/(1 - eps_eta)^2
expect_zero(
    "sign(d delta / d eps_eta) - 1",
    d_delta_deps_eta * (1 - eps_eta)**2 / delta0 - 1,
)
expect_zero(
    "sign(d M / d Z_W) - 1",
    dM_dZ * sp.pi**2 * (1 - eps_eta) * (1 - eps_W) / (8 * (1 + rho)**2) - 1,
)
expect_zero(
    "sign(d R / d Z_W) + 1",
    dR_dZ * Z_W**2 * (1 + rho)**2 / (Lambda * (1 - eps_eta) * (1 - eps_W)**2) + 1,
)
expect_zero(
    "sign(d M / d eps_eta) - 1",
    dM_deps_eta * sp.pi**2 * (1 - eps_eta)**2 * (1 - eps_W) / (8 * Z_W * (1 + rho)**2) - 1,
)
expect_zero(
    "sign(d R / d eps_eta) + 1",
    dR_deps_eta * Z_W * (1 + rho)**2 / (Lambda * (1 - eps_W)**2) + 1,
)
expect_zero(
    "sign(d M / d eps_W) - 1",
    dM_deps_W * sp.pi**2 * (1 - eps_eta) * (1 - eps_W)**2 / (8 * Z_W * (1 + rho)**2) - 1,
)
expect_zero(
    "sign(d R / d eps_W) + 1",
    dR_deps_W * Z_W * (1 + rho)**2 / (2 * Lambda * (1 - eps_eta) * (1 - eps_W)) + 1,
)
expect_zero(
    "sign(d M / d rho) - 1",
    dM_drho * sp.pi**2 * (1 - eps_eta) * (1 - eps_W) / (16 * Z_W * (1 + rho)) - 1,
)
expect_zero(
    "sign(d R / d rho) + 1",
    dR_drho * Z_W * (1 + rho)**3 / (2 * Lambda * (1 - eps_eta) * (1 - eps_W)**2) + 1,
)

print("On the natural nonvanishing transfer branch 1+rho > 0:")
print("  delta increases with eps_eta")
print("  M_mix increases with eps_eta, eps_W, Z_W, rho")
print("  R_target decreases with eps_eta, eps_W, Z_W, rho")

banner("STAGE 21 AUDIT COMPLETE")
print("The Stage-20 continuum operator collapses to an exact five-ratio placement map")
print("together with the product relation R_target M_mix = 8 Lambda (1-eps_W)/pi^2.")
