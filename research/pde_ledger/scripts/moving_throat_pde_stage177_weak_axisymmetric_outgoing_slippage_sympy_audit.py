#!/usr/bin/env python3
"""
moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py

SymPy audit for Stage 177.

Checks:
1. Weak-axisymmetric grouped slopes of M_r, I_r, H_r.
2. Exact collapse of the portwise defect to Sigma_{A,r} = eps lambda_A sigma_r.
3. Exact grouped trace/anomaly law for Xi_A = eps lambda_A Xi_1.
4. Exact identification of the weak-axisymmetric outgoing-prefactor slope.
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


banner("STAGE 160 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE")

eps, lam = sp.symbols("eps lam", real=True)
K, OU2, OW2, R, GU, GW = sp.symbols("K OU2 OW2 R GU GW", positive=True, real=True)
kappa1, gW1, gU1, r1, oU1, oW1 = sp.symbols("kappa1 gW1 gU1 r1 oU1 oW1", real=True)

# Background microscopic port invariants.
Mcal = GW / (sp.sqrt(K) * OW2)
Ical = R * GU / (OU2 * GW)
Hcal = R**2 / (OU2 * OW2)
Lambda = (OU2 * GW + R * GU) / (OU2 * OW2 - R**2)

# Weak-axisymmetric grouped perturbation in one lane with generic lambda_A = lam.
Kp = K * sp.exp(eps * lam * kappa1)
GWp = GW * sp.exp(eps * lam * gW1)
GUp = GU * sp.exp(eps * lam * gU1)
Rp = R * sp.exp(eps * lam * r1)
OU2p = OU2 * sp.exp(eps * lam * oU1)
OW2p = OW2 * sp.exp(eps * lam * oW1)

Mcal_p = GWp / (sp.sqrt(Kp) * OW2p)
Ical_p = Rp * GUp / (OU2p * GWp)
Hcal_p = Rp**2 / (OU2p * OW2p)
Lambda_p = (OU2p * GWp + Rp * GUp) / (OU2p * OW2p - Rp**2)

# Exact weak-axisymmetric logarithmic slopes.
m_r = gW1 - oW1 - sp.Rational(1, 2) * kappa1
i_r = r1 + gU1 - oU1 - gW1
h_r = 2 * r1 - oU1 - oW1

# Verify each microscopic slippage inherits the same lambda_A pattern.
dlnM_exact = sp.simplify(sp.expand(sp.series(sp.log(Mcal_p / Mcal), eps, 0, 2).removeO() / eps))
dlnI_exact = sp.simplify(sp.expand(sp.series(sp.log(Ical_p / Ical), eps, 0, 2).removeO() / eps))
dlnH_exact = sp.simplify(sp.expand(sp.series(sp.log(Hcal_p / Hcal), eps, 0, 2).removeO() / eps))

expect_zero("weak-axisymmetric d ln M", dlnM_exact - lam * m_r)
expect_zero("weak-axisymmetric d ln I", dlnI_exact - lam * i_r)
expect_zero("weak-axisymmetric d ln H", dlnH_exact - lam * h_r)

banner("Portwise outgoing-defect amplitude")

Sigma_exact = sp.simplify(
    sp.expand(
        sp.series(
            sp.log((Lambda_p**2 / Kp) / (Lambda**2 / K)),
            eps,
            0,
            2,
        ).removeO() / eps
    )
)

sigma_r = sp.simplify(
    2 * m_r
    + 2 * Ical / (1 + Ical) * i_r
    + 2 * Hcal / (1 - Hcal) * h_r
)
expect_zero("Sigma_{A,r} = lambda_A sigma_r", Sigma_exact - lam * sigma_r)

banner("Grouped weak-axisymmetric trace/anomaly law")

Xi1, epsilon = sp.symbols("Xi1 epsilon", real=True)
lam20 = sp.Integer(1)
lam21 = sp.Rational(1, 2)
lam22 = -sp.Integer(1)

Xi20 = epsilon * lam20 * Xi1
Xi21 = epsilon * lam21 * Xi1
Xi22 = epsilon * lam22 * Xi1

Xi_bar = sp.simplify((Xi20 + 2 * Xi21 + 2 * Xi22) / 5)
a_Xi = sp.simplify((2 * Xi20 - Xi21 - Xi22) / 10)
b_Xi = sp.simplify((Xi21 - Xi22) / 2)

expect_zero("grouped trace vanishes", Xi_bar)
expect_zero("a_Xi - eps Xi1/4", a_Xi - epsilon * Xi1 / 4)
expect_zero("b_Xi - 3 eps Xi1/4", b_Xi - 3 * epsilon * Xi1 / 4)
expect_zero("axisymmetric defect law b_Xi - 3 a_Xi", b_Xi - 3 * a_Xi)

banner("Outgoing-prefactor slope identification")

P0, D0, N0, d1, n1 = sp.symbols("P0 D0 N0 d1 n1", positive=True, real=True)
P0_exact = N0 / D0
DA = D0 * sp.exp(eps * lam * d1)
NA = N0 * sp.exp(eps * lam * n1)
PA = sp.simplify(NA / DA)
P_slope_series = sp.expand(sp.series(PA, eps, 0, 2).removeO())
P_slope_exact = sp.simplify(P_slope_series.coeff(eps, 1) / lam)
Xi_load = sp.simplify(n1 - d1)
expect_zero("P_A slope = P0 * Xi_load", P_slope_exact - P0_exact * Xi_load)
expect_zero("(P1/P0) - Xi_load", sp.simplify(P_slope_exact / P0_exact - Xi_load))

print("\nCarry-forward formulas:")
print("  m_r = gW1 - oW1 - kappa1/2")
print("  i_r = r1 + gU1 - oU1 - gW1")
print("  h_r = 2 r1 - oU1 - oW1")
print("  sigma_r = 2 m_r + 2 I/(1+I) i_r + 2 H/(1-H) h_r")
print("  Xi_A = eps lambda_A Xi1  =>  abar=0, a=eps Xi1/4, b=3 eps Xi1/4")
print("  On the conservative-shape-preserving branch, Xi1 = P1/P0.")
