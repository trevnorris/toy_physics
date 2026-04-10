#!/usr/bin/env python3
"""
5pn_stage26_continuum_kernel_extraction_and_placement_map.py

Twenty-sixth executable SymPy audit for the 5PN grouped-real-P2 / moving-throat
program.

What this script does
---------------------
1. Extracts the reduced branch data A, DeltaK_ax, alpha_mix, beta_0, M_mix, and
   delta from the first explicit finite-throat continuum kernel.
2. Compresses the microscopic operator into the exact dimensionless placement
   ratios (eps_eta, eps_W, rho, Z_W, delta_0, Lambda).
3. Proves the exact placement formulas for delta, M_mix, and R_target.
4. Proves the exact product law and first derivative tendencies.

Interpretation
--------------
The continuum selected-branch placement problem factorizes into a geometry lane,
a mixed-stability/product lane, and a redistribution lane.
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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("I. CONTINUUM-KERNEL EXTRACTION OF THE REDUCED BRANCH DATA")

# Exact D/N overlap norm.
sigma = sp.Rational(88, 1) / (9 * sp.pi**2)
print("sigma =", sigma)

# Continuum-kernel parameters.
K_U = sp.symbols("K_U", positive=True, real=True)
K_eta_eff, K_W_eff = sp.symbols("K_eta_eff K_W_eff", positive=True, real=True)
mu_eta, mu_U, mu_W = sp.symbols("mu_eta mu_U mu_W", positive=True, real=True)
c_etaU, c_etaW, c_UW = sp.symbols("c_etaU c_etaW c_UW", real=True)
T_w, L = sp.symbols("T_w L", positive=True, real=True)
G_N, c_s, a_th, c_light = sp.symbols("G c_s a_th c", positive=True, real=True)

A = sp.simplify((K_U * K_eta_eff - c_etaU**2) / (mu_eta * K_U))
DeltaK_ax = sp.simplify(sp.pi**2 * T_w / (mu_eta * L**2))
Delta0 = sp.simplify((K_U * K_W_eff - c_UW**2 * sigma) / (mu_U * mu_W))
Chi = sp.simplify((K_U * c_etaW + c_UW * c_etaU) / (mu_U * sp.sqrt(mu_eta * mu_W)))
alpha_mix = sp.simplify((K_U * c_etaW + c_UW * c_etaU) ** 2 / (mu_eta * K_U * (K_U * K_W_eff - c_UW**2 * sigma)))
beta0 = sp.simplify((mu_W / mu_eta) * (K_U * c_etaW + c_UW * c_etaU) ** 2 / (K_U * K_W_eff - c_UW**2 * sigma) ** 2)
M_mix = sp.simplify(8 * alpha_mix / (sp.pi**2 * A))
delta = sp.simplify(DeltaK_ax / A)

print("A =")
sp.pprint(A)
print("DeltaK_ax =")
sp.pprint(DeltaK_ax)
print("Delta_0 =")
sp.pprint(Delta0)
print("Chi =")
sp.pprint(Chi)
print("alpha_mix =")
sp.pprint(alpha_mix)
print("beta_0 =")
sp.pprint(beta0)
print("M_mix =")
sp.pprint(M_mix)
print("delta =")
sp.pprint(delta)

subbanner("I.1 — Exact continuum stability inequalities")
print("A > 0  <=>  K_U K_eta^(eff) > c_(etaU)^2")
print("Delta_0 > 0  <=>  K_U K_W^(eff) > c_(UW)^2 sigma")

banner("II. EXACT DIMENSIONLESS PLACEMENT RATIOS")

eps_eta, eps_W, rho, Z_W, delta0, Lambda = sp.symbols(
    "eps_eta eps_W rho Z_W delta_0 Lambda", positive=True, real=True
)

# Define the exact placement ratios from the continuum kernel.
eps_eta_expr = sp.simplify(c_etaU**2 / (K_U * K_eta_eff))
eps_W_expr = sp.simplify(c_UW**2 * sigma / (K_U * K_W_eff))
rho_expr = sp.simplify(c_UW * c_etaU / (K_U * c_etaW))
Z_W_expr = sp.simplify(c_etaW**2 / (K_eta_eff * K_W_eff))
delta0_expr = sp.simplify(sp.pi**2 * T_w / (L**2 * K_eta_eff))
Lambda_expr = sp.simplify(27 * sp.pi**2 * G_N * c_s**5 * K_W_eff / (20 * a_th**5 * c_light**5 * mu_W))

print("eps_eta =")
sp.pprint(eps_eta_expr)
print("eps_W =")
sp.pprint(eps_W_expr)
print("rho =")
sp.pprint(rho_expr)
print("Z_W =")
sp.pprint(Z_W_expr)
print("delta_0 =")
sp.pprint(delta0_expr)
print("Lambda =")
sp.pprint(Lambda_expr)

# Exact placement formulas.
delta_expected = sp.simplify(delta0_expr / (1 - eps_eta_expr))
M_mix_expected = sp.simplify(8 * Z_W_expr * (1 + rho_expr) ** 2 / (sp.pi**2 * (1 - eps_eta_expr) * (1 - eps_W_expr)))
R_target_expr = sp.simplify((54 * G_N * c_s**5 / (5 * a_th**5 * c_light**5)) * A / (beta0 * (8 / sp.pi**2)))
R_target_expected = sp.simplify(Lambda_expr * (1 - eps_eta_expr) * (1 - eps_W_expr) ** 2 / (Z_W_expr * (1 + rho_expr) ** 2))

expect_zero("delta - delta_0/(1-eps_eta)", delta - delta_expected)
expect_zero("M_mix - expected placement formula", M_mix - M_mix_expected)
expect_zero("R_target - expected placement formula", R_target_expr - R_target_expected)

print("R_target =")
sp.pprint(R_target_expr)

subbanner("II.1 — Outgoing transfer factor in dimensionless form")
mu_ratio = sp.symbols("mu_ratio", positive=True, real=True)
beta0_expected = sp.simplify((mu_W / mu_eta) * (K_eta_eff / K_W_eff) * Z_W_expr * (1 + rho_expr) ** 2 / (1 - eps_W_expr) ** 2)
expect_zero("beta_0 - expected dimensionless form", beta0 - beta0_expected)
print("beta_0 =")
sp.pprint(beta0_expected)

banner("III. EXACT PRODUCT LAW")

product_expr = sp.simplify(R_target_expr * M_mix)
product_expected_1 = sp.simplify(8 * Lambda_expr * (1 - eps_W_expr) / sp.pi**2)
product_expected_2 = sp.simplify(54 * G_N * c_s**5 * K_W_eff * (1 - eps_W_expr) / (5 * a_th**5 * c_light**5 * mu_W))
product_expected_3 = sp.simplify((54 * G_N * c_s**5 / (5 * a_th**5 * c_light**5)) * Delta0 / (K_U / mu_U))

expect_zero("R_target M_mix - 8 Lambda (1-eps_W)/pi^2", product_expr - product_expected_1)
expect_zero("R_target M_mix - 54 G c_s^5 K_W^(eff)(1-eps_W)/(5 a^5 c^5 mu_W)", product_expr - product_expected_2)
expect_zero("R_target M_mix - N_Q_target Delta_0 / Omega_U^2", product_expr - product_expected_3)

print("R_target M_mix =")
sp.pprint(product_expr)

banner("IV. EXACT ONE-WAY PARAMETER TENDENCIES")

# Treat the dimensionless ratios as independent variables for derivative sign structure.
eps_eta_s, eps_W_s, rho_s, Z_W_s, delta0_s, Lambda_s = sp.symbols(
    "eps_eta_s eps_W_s rho_s Z_W_s delta0_s Lambda_s",
    positive=True,
    real=True,
)

delta_sym = sp.simplify(delta0_s / (1 - eps_eta_s))
M_mix_sym = sp.simplify(8 * Z_W_s * (1 + rho_s) ** 2 / (sp.pi**2 * (1 - eps_eta_s) * (1 - eps_W_s)))
R_target_sym = sp.simplify(Lambda_s * (1 - eps_eta_s) * (1 - eps_W_s) ** 2 / (Z_W_s * (1 + rho_s) ** 2))

print("d delta / d eps_eta =")
sp.pprint(sp.simplify(sp.diff(delta_sym, eps_eta_s)))
print("d M_mix / d eps_eta =")
sp.pprint(sp.simplify(sp.diff(M_mix_sym, eps_eta_s)))
print("d R_target / d eps_eta =")
sp.pprint(sp.simplify(sp.diff(R_target_sym, eps_eta_s)))

print("d M_mix / d eps_W =")
sp.pprint(sp.simplify(sp.diff(M_mix_sym, eps_W_s)))
print("d R_target / d eps_W =")
sp.pprint(sp.simplify(sp.diff(R_target_sym, eps_W_s)))

print("d M_mix / d Z_W =")
sp.pprint(sp.simplify(sp.diff(M_mix_sym, Z_W_s)))
print("d R_target / d Z_W =")
sp.pprint(sp.simplify(sp.diff(R_target_sym, Z_W_s)))

print("d M_mix / d rho =")
sp.pprint(sp.simplify(sp.diff(M_mix_sym, rho_s)))
print("d R_target / d rho =")
sp.pprint(sp.simplify(sp.diff(R_target_sym, rho_s)))

banner("V. THREE-LANE FACTORIZATION")
print("Geometry lane:")
print("  delta = delta_0 / (1 - eps_eta)")
print("Mixed-stability/product lane:")
print("  R_target M_mix = 8 Lambda (1 - eps_W) / pi^2")
print("Redistribution lane:")
print("  (eps_eta, Z_W, rho) move the branch along the fixed product curve.")

banner("FINAL LEDGER")
print("1. The first explicit continuum kernel already determines A, DeltaK_ax, beta_0, alpha_mix, M_mix, and delta exactly.")
print("2. The continuum placement problem compresses to the five ratios (eps_eta, eps_W, rho, Z_W, delta_0) plus Lambda.")
print("3. The exact product law shows that eps_eta, Z_W, and rho only redistribute the defect along a fixed product curve.")
print("4. The remaining theorem gate is whether the true moving-throat branch lands on the admissible F/G geometry with these continuum ratios.")
