#!/usr/bin/env python3
"""
moving_throat_pde_stage161_outgoing_port_coloading_sympy_audit.py

SymPy-backed audit for Stage 161.

Checks:
1. Exact weak-axisymmetric slope of the static outgoing-transfer coefficient
   N_{A,0}^{(r)} = P_{A,r}^2 / Delta_{A,r}^2.
2. Exact collapse Xi_1 = sum_r rho_r^(N) (nu_r - kappa_1) = nubar_N - kappa_1.
3. Exact formulas for the actual port slopes p_r and d_r.
4. Exact equivalence between the direct port-slope formula and the Stage-160
   slippage formula sigma_r = nu_r - kappa_1.
"""

from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 161 — OUTGOING-PORT CO-LOADING THEOREM")

# ---------------------------------------------------------------------------
# 1. Static outgoing-transfer slope
# ---------------------------------------------------------------------------

banner("1. Exact weak-axisymmetric slope of N_{A,0}^{(r)}")

eps, lam = sp.symbols('eps lam', real=True)
P0r, D0r, p_r, d_r = sp.symbols('P0r D0r p_r d_r', positive=True, real=True)

PAr = P0r * (1 + eps * lam * p_r)
DAr = D0r * (1 + eps * lam * d_r)
N0r = P0r**2 / D0r**2
NAr = sp.expand(PAr**2 / DAr**2)

NAr_lin = sp.expand(sp.series(NAr, eps, 0, 2).removeO())
nu_from_series = sp.simplify((NAr_lin / N0r - 1) / (eps * lam))
expect_zero("nu_r - 2(p_r-d_r)", nu_from_series - 2 * (p_r - d_r))

print("N_{A,0}^{(r)} =")
sp.pprint(NAr_lin)
print("nu_r =", sp.simplify(2 * (p_r - d_r)))

# ---------------------------------------------------------------------------
# 2. Weighted collapse of Xi_1
# ---------------------------------------------------------------------------

banner("2. Weighted collapse Xi_1 = nubar_N - kappa_1")

rho1, rho2, rho3 = sp.symbols('rho1 rho2 rho3', real=True)
nu1, nu2, nu3, kappa1 = sp.symbols('nu1 nu2 nu3 kappa1', real=True)

Xi = sp.expand(rho1 * (nu1 - kappa1) + rho2 * (nu2 - kappa1) + rho3 * (nu3 - kappa1))
nu_bar = sp.expand(rho1 * nu1 + rho2 * nu2 + rho3 * nu3)
expect_zero(
    "Xi_1 - (nubar_N - kappa_1)",
    Xi.subs(rho3, 1 - rho1 - rho2) - (nu_bar.subs(rho3, 1 - rho1 - rho2) - kappa1),
)

print("Xi_1 =", Xi)
print("nubar_N =", nu_bar)

# ---------------------------------------------------------------------------
# 3. Exact formulas for p_r and d_r from actual port data
# ---------------------------------------------------------------------------

banner("3. Exact formulas for p_r and d_r from actual port data")

OU2, OW2, R, GU, GW = sp.symbols('OU2 OW2 R GU GW', positive=True, real=True)
oU, oW, rr, gU, gW = sp.symbols('oU oW rr gU gW', real=True)

P = OU2 * GW + R * GU
Delta = OU2 * OW2 - R**2

P_A = OU2 * GW * (1 + eps * lam * (oU + gW)) + R * GU * (1 + eps * lam * (rr + gU))
D_A = OU2 * OW2 * (1 + eps * lam * (oU + oW)) - R**2 * (1 + 2 * eps * lam * rr)

p_from_series = sp.simplify((sp.series(P_A / P, eps, 0, 2).removeO() - 1) / (eps * lam))
d_from_series = sp.simplify((sp.series(D_A / Delta, eps, 0, 2).removeO() - 1) / (eps * lam))

alpha = sp.simplify(OU2 * GW / P)
beta = sp.simplify(R * GU / P)
chi = sp.simplify(OU2 * OW2 / Delta)
zeta = sp.simplify(R**2 / Delta)

p_expected = sp.simplify(alpha * (oU + gW) + beta * (rr + gU))
d_expected = sp.simplify(chi * (oU + oW) - 2 * zeta * rr)

expect_zero("p_r formula", p_from_series - p_expected)
expect_zero("d_r formula", d_from_series - d_expected)
expect_zero("alpha+beta-1", sp.simplify(alpha + beta - 1))
expect_zero("chi-zeta-1", sp.simplify(chi - zeta - 1))

print("alpha_r =", alpha)
print("beta_r  =", beta)
print("chi_r   =", chi)
print("zeta_r  =", zeta)
print("p_r     =", p_expected)
print("d_r     =", d_expected)

# ---------------------------------------------------------------------------
# 4. Equivalence with the Stage-160 slippage formula
# ---------------------------------------------------------------------------

banner("4. Equivalence with the Stage-160 slippage formula")

Ir_expr = sp.simplify(R * GU / (OU2 * GW))
Hr_expr = sp.simplify(R**2 / (OU2 * OW2))

m_r = gW - oW - sp.Rational(1, 2) * kappa1
i_r = rr + gU - oU - gW
h_r = 2 * rr - oU - oW

nu_direct = sp.simplify(2 * p_expected - 2 * d_expected)
nu_expected = sp.simplify(
    kappa1
    + 2 * m_r
    + (2 * Ir_expr / (1 + Ir_expr)) * i_r
    + (2 * Hr_expr / (1 - Hr_expr)) * h_r
)
expect_zero("nu_r - [kappa1 + sigma_r]", nu_direct - nu_expected)

sigma_r = sp.simplify(nu_expected - kappa1)
print("Ir =", Ir_expr)
print("Hr =", Hr_expr)
print("nu_r =", nu_expected)
print("sigma_r =", sigma_r)
print("Hence Xi_1 = sum_r rho_r^(N) sigma_r = sum_r rho_r^(N)(nu_r-kappa1).")
