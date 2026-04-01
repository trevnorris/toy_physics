#!/usr/bin/env python3
"""
moving_throat_pde_stage151_off_bundle_slippage_sympy_audit.py

SymPy-backed audit for the exact first-order off-bundle slippage decomposition.

Checks:
1. Stage-147 normal coordinate expressed in microscopic log drifts.
2. Exact cancellation of the lower-branch transport laws from Stage 148.
3. Exact reduction of the normal coordinate to the three slippage variables
   eps_L, eps_v, eps_T.
4. Exact transport of the mouth-bias variation and the outlet defects.
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

banner("STAGE 151 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION")

# Symbols
g, r = sp.symbols("g r", positive=True, real=True)
s = sp.sqrt(1 + r**2)

Xi, csw, cs, a = sp.symbols("Xi csw cs a", real=True)
LW, vw0, Tm = sp.symbols("LW vw0 Tm", real=True)
epsL, epsv, epsT = sp.symbols("epsL epsv epsT", real=True)

A = g + sp.Rational(1, 4) / s
B = sp.Rational(1, 2) / s
C = 2 * g + sp.Rational(3, 4) / s

delta_perp = (
    A * Xi
    + 3 * A * csw
    + B * cs
    - g * Tm
    - (g + B) * vw0
    - 2 * A * a
    - C * LW
)

print("delta_perp =", sp.simplify(delta_perp))

banner("Stage-148 lower-branch transport laws")
Tm_br = sp.Rational(1, 2) * Xi + sp.Rational(3, 2) * csw - cs - sp.Rational(3, 2) * a
vw0_br = sp.Rational(1, 2) * Xi + sp.Rational(3, 2) * csw + cs - sp.Rational(5, 2) * a
LW_br = a

print("Tm_br  =", sp.simplify(Tm_br))
print("vw0_br =", sp.simplify(vw0_br))
print("LW_br  =", sp.simplify(LW_br))

expect_zero(
    "bundle tangency (delta_perp on exact lower branch)",
    delta_perp.subs({Tm: Tm_br, vw0: vw0_br, LW: LW_br}),
)

banner("Introduce the three off-bundle slippages")
subs_slip = {
    Tm: Tm_br + epsT,
    vw0: vw0_br + epsv,
    LW: LW_br + epsL,
}
delta_perp_slip = sp.expand(delta_perp.subs(subs_slip))
eps_perp = g * epsT + (g + B) * epsv + C * epsL

print("delta_perp with slippages =", sp.simplify(delta_perp_slip))
expect_zero("delta_perp + eps_perp", delta_perp_slip + eps_perp)

banner("Mouth-bias transport")
Sigma0, S_can, dSigma0, dS = sp.symbols("Sigma0 S_can dSigma0 dS", real=True)
deltaPi = (1 - sp.Rational(1, 4) * S_can) * dSigma0 - Sigma0 * dS / 4 + Sigma0 * S_can * delta_perp_slip / s
deltaPi_expected = (
    (1 - sp.Rational(1, 4) * S_can) * dSigma0
    - Sigma0 * dS / 4
    - Sigma0 * S_can * eps_perp / s
)
expect_zero("deltaPi transport identity", deltaPi - deltaPi_expected)

banner("Outlet defects")
sigma, dkappaW, dgammaW = sp.symbols("sigma dkappaW dgammaW", real=True)
dC = 16 * sigma * delta_perp_slip / s
dE2 = sigma * (16 * delta_perp_slip / s - 9 * dkappaW) / (27 * (1 - sigma))
dE4 = sigma * (80 * delta_perp_slip / s - 72 * dkappaW) / (243 * (1 - sigma))
DeltaQ = sigma * (16 * delta_perp_slip / s - 27 * dgammaW) / (3 * (1 - sigma))

expect_zero("dC + 16 sigma eps_perp / s", dC + 16 * sigma * eps_perp / s)
expect_zero(
    "dE2 - sigma(-16 eps_perp/s - 9 dkappaW)/(27(1-sigma))",
    dE2 - sigma * (-16 * eps_perp / s - 9 * dkappaW) / (27 * (1 - sigma)),
)
expect_zero(
    "dE4 - sigma(-80 eps_perp/s - 72 dkappaW)/(243(1-sigma))",
    dE4 - sigma * (-80 * eps_perp / s - 72 * dkappaW) / (243 * (1 - sigma)),
)
expect_zero(
    "DeltaQ - sigma(-16 eps_perp/s - 27 dgammaW)/(3(1-sigma))",
    DeltaQ - sigma * (-16 * eps_perp / s - 27 * dgammaW) / (3 * (1 - sigma)),
)

banner("Numeric Family-1 coefficients")
r_num = sp.Float("1.77799353547498", 30)
g_num = sp.Float("0.758035078944663", 30)
s_num = sp.sqrt(1 + r_num**2)
B_num = sp.Rational(1, 2) / s_num
C_num = 2 * g_num + sp.Rational(3, 4) / s_num
Sigma0_num = sp.Float("4.651033550168876", 30)
S_num = sp.Float("0.6703621156734617", 30)

print("coeff epsT in delta_perp =", sp.N(-g_num, 20))
print("coeff epsv in delta_perp =", sp.N(-(g_num + B_num), 20))
print("coeff epsL in delta_perp =", sp.N(-C_num, 20))

Pi_epsT = -Sigma0_num * S_num * g_num / s_num
Pi_epsv = -Sigma0_num * S_num * (g_num + B_num) / s_num
Pi_epsL = -Sigma0_num * S_num * C_num / s_num

print("coeff epsT in deltaPi =", sp.N(Pi_epsT, 20))
print("coeff epsv in deltaPi =", sp.N(Pi_epsv, 20))
print("coeff epsL in deltaPi =", sp.N(Pi_epsL, 20))

print("\nCarry-forward formulas:")
print("  eps_L := d ln L_W - d ln a")
print("  eps_v := d ln v_w0 - [1/2 d ln(Z_q/rho_w) + 3/2 d ln c_sw + d ln c_s - 5/2 d ln a]")
print("  eps_T := d ln T_m  - [1/2 d ln(Z_q/rho_w) + 3/2 d ln c_sw - d ln c_s - 3/2 d ln a]")
print("  delta_perp = -[ g_* eps_T + (g_* + 1/(2 sqrt(1+r_*^2))) eps_v + (2 g_* + 3/(4 sqrt(1+r_*^2))) eps_L ]")
print("  All first-order off-bundle conservative-even defects depend only on eps_perp.")
