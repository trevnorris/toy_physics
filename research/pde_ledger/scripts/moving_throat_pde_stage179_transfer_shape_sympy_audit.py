#!/usr/bin/env python3
"""
moving_throat_pde_stage179_transfer_shape_sympy_audit.py

SymPy-backed audit for the wall-normalized transfer-shape theorem.

Checks:
1. Exact factorization N0/K = T^2 in terms of wall-normalized port variables.
2. Exact weak-axisymmetric slope identity nu_r = kappa_1 + 2 tau_r.
3. Exact equivalence to the Stage 176/160/161 slippage formulas.
4. Weighted defect identity Xi_1 = 2 sum_r rho_r tau_r when sum_r rho_r = 1.
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


banner("STAGE 179 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM")

# Primitive static port data.
K = sp.symbols("K", positive=True, real=True)
OU2, OW2 = sp.symbols("OU2 OW2", positive=True, real=True)
GW, GU = sp.symbols("GW GU", positive=True, real=True)
R = sp.symbols("R", real=True)

P = OU2 * GW + R * GU
Delta = OU2 * OW2 - R**2
N0 = sp.simplify(P**2 / Delta**2)

# Wall-normalized variables.
GWh = GW / (sp.sqrt(K) * OW2)
GUh = GU / (sp.sqrt(K) * sp.sqrt(OU2) * sp.sqrt(OW2))
Rh = R / (sp.sqrt(OU2) * sp.sqrt(OW2))
T = sp.simplify((GWh + Rh * GUh) / (1 - Rh**2))

print("P   =", sp.simplify(P))
print("Delta =", sp.simplify(Delta))
print("T   =", sp.simplify(T))
expect_zero("N0/K - T^2", sp.simplify(N0 / K - T**2))

banner("Weak-axisymmetric slope identity")

eps, lam = sp.symbols("eps lam", real=True)
kappa1 = sp.symbols("kappa1", real=True)
gW, gU, oU, oW, rr = sp.symbols("gW gU oU oW rr", real=True)

KA = K * (1 + eps * lam * kappa1)
OU2A = OU2 * (1 + eps * lam * oU)
OW2A = OW2 * (1 + eps * lam * oW)
GWA = GW * (1 + eps * lam * gW)
GUA = GU * (1 + eps * lam * gU)
RA = R * (1 + eps * lam * rr)

PA = sp.expand(OU2A * GWA + RA * GUA)
DeltaA = sp.expand(OU2A * OW2A - RA**2)
N0A = sp.expand(PA**2 / DeltaA**2)

# Direct weak-axisymmetric logarithmic slope of N0.
nu_direct = sp.simplify(sp.expand(sp.series(sp.log(N0A), eps, 0, 2).removeO()).coeff(eps, 1) / lam)
print("nu_direct =")
sp.pprint(nu_direct)

# Wall-normalized logarithmic slopes.
w = sp.simplify(gW - oW - sp.Rational(1, 2) * kappa1)
u = sp.simplify(gU - sp.Rational(1, 2) * oU - sp.Rational(1, 2) * oW - sp.Rational(1, 2) * kappa1)
c = sp.simplify(rr - sp.Rational(1, 2) * oU - sp.Rational(1, 2) * oW)

alpha = sp.simplify(GWh / (GWh + Rh * GUh))
beta = sp.simplify(Rh * GUh / (GWh + Rh * GUh))
tau = sp.simplify(alpha * w + beta * (u + c) + 2 * Rh**2 / (1 - Rh**2) * c)
nu_expected = sp.simplify(kappa1 + 2 * tau)

print("tau =")
sp.pprint(tau)
print("nu_expected =")
sp.pprint(nu_expected)
expect_zero("nu_direct - (kappa1 + 2 tau)", nu_direct - nu_expected)

banner("Exact equivalence to the Stage 176/160/161 slippage formulas")

Ih = sp.simplify(Rh * GUh / GWh)
Hh = sp.simplify(Rh**2)
m = w
i = sp.simplify((u + c) - w)
h = sp.simplify(2 * c)

tau_slippage = sp.simplify(m + Ih / (1 + Ih) * i + Hh / (1 - Hh) * h)
expect_zero("tau - slippage form", tau - tau_slippage)
expect_zero("(nu-kappa1) - 2*tau_slippage", (nu_expected - kappa1) - 2 * tau_slippage)

banner("Weighted defect identity")

rho1, rho2 = sp.symbols("rho1 rho2", real=True)
GW1, GW2, GW3 = sp.symbols("GW1 GW2 GW3", positive=True, real=True)
GU1, GU2, GU3 = sp.symbols("GU1 GU2 GU3", positive=True, real=True)
OU21, OU22, OU23 = sp.symbols("OU21 OU22 OU23", positive=True, real=True)
OW21, OW22, OW23 = sp.symbols("OW21 OW22 OW23", positive=True, real=True)
R1, R2, R3 = sp.symbols("R1 R2 R3", real=True)
gW1, gW2, gW3 = sp.symbols("gW1 gW2 gW3", real=True)
gU1, gU2, gU3 = sp.symbols("gU1 gU2 gU3", real=True)
oU1, oU2, oU3 = sp.symbols("oU1 oU2 oU3", real=True)
oW1, oW2, oW3 = sp.symbols("oW1 oW2 oW3", real=True)
rr1, rr2, rr3 = sp.symbols("rr1 rr2 rr3", real=True)
rho3 = 1 - rho1 - rho2

port1 = {GW: GW1, GU: GU1, OU2: OU21, OW2: OW21, R: R1, gW: gW1, gU: gU1, oU: oU1, oW: oW1, rr: rr1}
port2 = {GW: GW2, GU: GU2, OU2: OU22, OW2: OW22, R: R2, gW: gW2, gU: gU2, oU: oU2, oW: oW2, rr: rr2}
port3 = {GW: GW3, GU: GU3, OU2: OU23, OW2: OW23, R: R3, gW: gW3, gU: gU3, oU: oU3, oW: oW3, rr: rr3}

tau1 = sp.simplify(tau.subs(port1))
tau2 = sp.simplify(tau.subs(port2))
tau3 = sp.simplify(tau.subs(port3))
nu1 = sp.simplify(nu_direct.subs(port1))
nu2 = sp.simplify(nu_direct.subs(port2))
nu3 = sp.simplify(nu_direct.subs(port3))

expect_zero("nu_1 - (kappa1 + 2 tau_1)", nu1 - (kappa1 + 2 * tau1))
expect_zero("nu_2 - (kappa1 + 2 tau_2)", nu2 - (kappa1 + 2 * tau2))
expect_zero("nu_3 - (kappa1 + 2 tau_3)", nu3 - (kappa1 + 2 * tau3))

Xi = rho1 * (nu1 - kappa1) + rho2 * (nu2 - kappa1) + rho3 * (nu3 - kappa1)
Xi_expected = 2 * (rho1 * tau1 + rho2 * tau2 + rho3 * tau3)
expect_zero("Xi_1 - 2 weighted tau", Xi - Xi_expected)

print("\nCarry-forward formulas:")
print("  N0^(r) = K * T_r^2")
print("  T_r = (Ghat_W + Rhat Ghat_U)/(1-Rhat^2)")
print("  nu_r = kappa1 + 2 tau_r")
print("  tau_r = m_r + I_r/(1+I_r) i_r + H_r/(1-H_r) h_r")
print("  Xi_1 = 2 sum_r rho_r^(N) tau_r")
