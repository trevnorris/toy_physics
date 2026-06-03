#!/usr/bin/env python3
"""
moving_throat_pde_stage20_continuum_kernel_sympy_audit.py

Stage 20 SymPy audit:
derive the Stage-17/19 reduced branch data from the first explicit finite-throat
continuum kernel and verify the exact closed formulas for

    A, DeltaK_ax, alpha_mix, beta_0, M_mix, delta.

The script also verifies the exact Schur-complement factorization

    Sigma_wall = Xi I_2 + alpha v v^T

for the reduced continuum operator.
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
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(sp.expand(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 37 — CONTINUUM-KERNEL EXTRACTION AUDIT")

# ---------------------------------------------------------------------------
# 1. Exact finite-throat mode system
# ---------------------------------------------------------------------------

subbanner("1. Exact N/N and D/N modes")

s, L = sp.symbols("s L", positive=True, real=True)

u0 = 1 / sp.sqrt(L)
u1 = sp.sqrt(2 / L) * sp.cos(sp.pi * s / L)
f0 = sp.sqrt(2 / L) * sp.sin(sp.pi * s / (2 * L))

expect_zero("u0 normalization - 1", sp.integrate(u0 * u0, (s, 0, L)) - 1)
expect_zero("u1 normalization - 1", sp.integrate(u1 * u1, (s, 0, L)) - 1)
expect_zero("f0 normalization - 1", sp.integrate(f0 * f0, (s, 0, L)) - 1)
expect_zero("u0-u1 orthogonality", sp.integrate(u0 * u1, (s, 0, L)))
expect_zero("u0 N/N boundary slope @0", sp.diff(u0, s).subs(s, 0))
expect_zero("u1 N/N boundary slope @0", sp.diff(u1, s).subs(s, 0))
expect_zero("u1 N/N boundary slope @L", sp.diff(u1, s).subs(s, L))
expect_zero("f0 D/N boundary value @0", f0.subs(s, 0))
expect_zero("f0 D/N boundary slope @L", sp.diff(f0, s).subs(s, L))

kappa0 = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
kappa1 = sp.simplify(sp.integrate(u1 * f0, (s, 0, L)))
sigma = sp.simplify(kappa0**2 + kappa1**2)

print(f"kappa0 = {kappa0}")
print(f"kappa1 = {kappa1}")
print(f"sigma  = {sigma}")

expect_zero("kappa0 - 2 sqrt(2)/pi", kappa0 - 2 * sp.sqrt(2) / sp.pi)
expect_zero("kappa1 + 4/(3 pi)", kappa1 + 4 / (3 * sp.pi))
expect_zero("sigma - 88/(9 pi^2)", sigma - sp.Rational(88, 9) / sp.pi**2)

# ---------------------------------------------------------------------------
# 2. Mass-normalized modal kernels
# ---------------------------------------------------------------------------

subbanner("2. Mass-normalized projected kernels")

mu_eta, mu_U, mu_phi, mu_W = sp.symbols("mu_eta mu_U mu_phi mu_W", positive=True, real=True)
K_eta, T_Omega, T_w = sp.symbols("K_eta T_Omega T_w", positive=True, real=True)
K_U, K_phi, T_phi = sp.symbols("K_U K_phi T_phi", positive=True, real=True)
K_W, T_W = sp.symbols("K_W T_W", positive=True, real=True)

Keta_eff = sp.simplify(K_eta + 6 * T_Omega)
KWeff = sp.simplify(K_W + sp.pi**2 * T_W / (4 * L**2))

K0 = sp.simplify(Keta_eff / mu_eta)
DeltaK_ax = sp.simplify(sp.pi**2 * T_w / (mu_eta * L**2))
varpi2 = sp.simplify((K_phi + sp.pi**2 * T_phi / (4 * L**2)) / mu_phi)
OmegaU2 = sp.simplify(K_U / mu_U)
OmegaW2 = sp.simplify(KWeff / mu_W)

print(f"K0          = {K0}")
print(f"DeltaK_ax   = {DeltaK_ax}")
print(f"varpi^2     = {varpi2}")
print(f"Omega_U^2   = {OmegaU2}")
print(f"Omega_W^2   = {OmegaW2}")

# ---------------------------------------------------------------------------
# 3. Exact Schur complement of the internal block
# ---------------------------------------------------------------------------

subbanner("3. Schur-complement factorization")

A_U, A_phi, A_W = sp.symbols("A_U A_phi A_W", real=True)
gU, gB, gW, gR = sp.symbols("gU gB gW gR", real=True)

B = sp.Matrix([
    [A_U, 0,   0,         -gR * kappa0],
    [0,   A_U, 0,         -gR * kappa1],
    [0,   0,   A_phi,      0],
    [-gR * kappa0, -gR * kappa1, 0, A_W],
])

C = sp.Matrix([
    [gU, 0,  gB * kappa0, gW * kappa0],
    [0,  gU, gB * kappa1, gW * kappa1],
])

Sigma = sp.simplify(C * B.inv() * C.T)

Delta_UW = sp.simplify(A_U * A_W - gR**2 * sigma)
Xi = sp.simplify(gU**2 / A_U)
alpha = sp.simplify(gB**2 / A_phi + (A_U * gW + gR * gU) ** 2 / (A_U * Delta_UW))
Sigma_expected = sp.simplify(Xi * sp.eye(2) + alpha * sp.Matrix([
    [kappa0**2, kappa0 * kappa1],
    [kappa0 * kappa1, kappa1**2],
]))

expect_zero("Sigma_wall - [Xi I + alpha v v^T]", Sigma - Sigma_expected)

# ---------------------------------------------------------------------------
# 4. Continuum couplings and closed reduced formulas
# ---------------------------------------------------------------------------

subbanner("4. Continuum extraction of A, alpha_mix, beta_0, M_mix, delta")

c_etaU, c_etaPhi, c_etaW, c_UW = sp.symbols("c_etaU c_etaPhi c_etaW c_UW", real=True)

gU_cont = sp.simplify(c_etaU / sp.sqrt(mu_eta * mu_U))
gB_cont = sp.simplify(c_etaPhi / sp.sqrt(mu_eta * mu_phi))
gW_cont = sp.simplify(c_etaW / sp.sqrt(mu_eta * mu_W))
gR_cont = sp.simplify(c_UW / sp.sqrt(mu_U * mu_W))

A = sp.simplify(K0 - gU_cont**2 / OmegaU2)
Delta0 = sp.simplify(OmegaU2 * OmegaW2 - gR_cont**2 * sigma)
Chi = sp.simplify(OmegaU2 * gW_cont + gR_cont * gU_cont)
beta0 = sp.simplify(Chi**2 / Delta0**2)
alpha_mix = sp.simplify(Chi**2 / (OmegaU2 * Delta0))
M_mix = sp.simplify(8 * alpha_mix / (sp.pi**2 * A))
delta = sp.simplify(DeltaK_ax / A)

A_expected = sp.simplify((K_U * Keta_eff - c_etaU**2) / (mu_eta * K_U))
Delta0_expected = sp.simplify((K_U * KWeff - c_UW**2 * sigma) / (mu_U * mu_W))
Chi_expected = sp.simplify((K_U * c_etaW + c_UW * c_etaU) / (mu_U * sp.sqrt(mu_eta * mu_W)))
beta0_expected = sp.simplify(
    (mu_W / mu_eta) * (K_U * c_etaW + c_UW * c_etaU) ** 2
    / ((K_U * KWeff - c_UW**2 * sigma) ** 2)
)
alpha_mix_expected = sp.simplify(
    (K_U * c_etaW + c_UW * c_etaU) ** 2
    / (mu_eta * K_U * (K_U * KWeff - c_UW**2 * sigma))
)
M_mix_expected = sp.simplify(
    8 * (K_U * c_etaW + c_UW * c_etaU) ** 2
    / (sp.pi**2 * (K_U * Keta_eff - c_etaU**2) * (K_U * KWeff - c_UW**2 * sigma))
)
delta_expected = sp.simplify(
    sp.pi**2 * T_w * K_U / (L**2 * (K_U * Keta_eff - c_etaU**2))
)

expect_zero("A continuum formula", A - A_expected)
expect_zero("Delta0 continuum formula", Delta0 - Delta0_expected)
expect_zero("Chi continuum formula", Chi - Chi_expected)
expect_zero("beta0 continuum formula", beta0 - beta0_expected)
expect_zero("alpha_mix continuum formula", alpha_mix - alpha_mix_expected)
expect_zero("M_mix continuum formula", M_mix - M_mix_expected)
expect_zero("delta continuum formula", delta - delta_expected)

print(f"A         = {A_expected}")
print(f"Delta0    = {Delta0_expected}")
print(f"Chi       = {Chi_expected}")
print(f"beta0     = {beta0_expected}")
print(f"alpha_mix = {alpha_mix_expected}")
print(f"M_mix     = {M_mix_expected}")
print(f"delta     = {delta_expected}")

# ---------------------------------------------------------------------------
# 5. Exact continuum stability surfaces
# ---------------------------------------------------------------------------

subbanner("5. Exact continuum stability inequalities")

expect_zero(
    "A - [K_U K_eta^(eff) - c_(etaU)^2]/(mu_eta K_U)",
    A - (K_U * Keta_eff - c_etaU**2) / (mu_eta * K_U),
)
expect_zero(
    "Delta0 - [K_U K_W^(eff) - c_(UW)^2 sigma]/(mu_U mu_W)",
    Delta0 - (K_U * KWeff - c_UW**2 * sigma) / (mu_U * mu_W),
)

print("Stable selected branch requires:")
print("  K_U K_eta^(eff) > c_(etaU)^2")
print("  K_U K_W^(eff)   > c_(UW)^2 sigma")
print("with sigma =", sigma)

banner("STAGE 37 AUDIT COMPLETE")
print("The first explicit finite-throat continuum kernel reproduces the full reduced")
print("Stage-17/19 branch data exactly.")
