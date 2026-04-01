#!/usr/bin/env python3
"""
moving_throat_pde_stage159_outgoing_load_factorization_sympy_audit.py

SymPy audit for Stage 159.

Checks:
1. Exact factorization of Lambda_r^2 / K into (M_r, I_r, H_r).
2. Exact first-order defect decomposition.
3. Expanded primitive-variable transport law.
4. Rigidity corollary reducing the defect to the square-root mixed-leg law.
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


banner("STAGE 159 — OUTGOING LOAD-FACTOR FACTORIZATION")

K, OU2, OW2, R, GU, GW = sp.symbols("K OU2 OW2 R GU GW", positive=True, real=True)

Lambda = (OU2 * GW + R * GU) / (OU2 * OW2 - R**2)
Mcal = GW / (sp.sqrt(K) * OW2)
Ical = R * GU / (OU2 * GW)
Hcal = R**2 / (OU2 * OW2)

expect_zero(
    "exact factorization of Lambda^2/K",
    Lambda**2 / K - Mcal**2 * (1 + Ical) ** 2 / (1 - Hcal) ** 2,
)

print("\nCarry-forward exact identity:")
print("  Lambda^2/K = M^2 (1+I)^2 / (1-H)^2")

banner("First-order logarithmic transport")

eps = sp.symbols("eps", real=True)
dK, dOU, dOW, dR, dGU, dGW = sp.symbols("dK dOU dOW dR dGU dGW", real=True)

Kp = K * sp.exp(eps * dK)
OU2p = OU2 * sp.exp(eps * dOU)
OW2p = OW2 * sp.exp(eps * dOW)
Rp = R * sp.exp(eps * dR)
GUp = GU * sp.exp(eps * dGU)
GWp = GW * sp.exp(eps * dGW)

Lambdap = (OU2p * GWp + Rp * GUp) / (OU2p * OW2p - Rp**2)
Sigma_exact = sp.simplify(
    sp.expand(
        sp.series(
            sp.log((Lambdap**2 / Kp) / (Lambda**2 / K)),
            eps,
            0,
            2,
        ).removeO()
        / eps
    )
)

# Microscopic logarithmic drifts.
dlnM = dGW - dOW - sp.Rational(1, 2) * dK
dlnI = dR + dGU - dOU - dGW
dlnH = 2 * dR - dOU - dOW

Sigma_factored = sp.simplify(
    2 * dlnM + 2 * Ical / (1 + Ical) * dlnI + 2 * Hcal / (1 - Hcal) * dlnH
)
expect_zero("factored first-order defect formula", Sigma_exact - Sigma_factored)

Sigma_expanded = sp.simplify(
    -dK
    + 2 / (1 + Ical) * dGW
    + 2 * Ical / (1 + Ical) * dGU
    + 2 * (Ical / (1 + Ical) + 2 * Hcal / (1 - Hcal)) * dR
    - 2 * (Ical / (1 + Ical) + Hcal / (1 - Hcal)) * dOU
    - 2 / (1 - Hcal) * dOW
)
expect_zero("expanded primitive-variable transport", Sigma_exact - Sigma_expanded)

print("\nFirst-order formulas:")
print("  d ln M = d ln G_W - d ln Omega_W^2 - 1/2 dK")
print("  d ln I = d ln R + d ln G_U - d ln Omega_U^2 - d ln G_W")
print("  d ln H = 2 d ln R - d ln Omega_U^2 - d ln Omega_W^2")
print("  Sigma^(N) = 2 d ln M + 2 I/(1+I) d ln I + 2 H/(1-H) d ln H")

banner("Rigidity corollary")

Sigma_rigid = sp.simplify(Sigma_factored.subs({dlnI: 0, dlnH: 0}))
expect_zero("rigidity reduction to 2 d ln M", Sigma_rigid - 2 * dlnM)

print("\nConclusion:")
print("  If the interference and hybridization ratios are rigid, then")
print("      Sigma^(N) = 2 d ln[ G_W / (Omega_W^2 sqrt(K)) ].")
print("  So the remaining linear grouped defect vanishes when the raw mixed leg")
print("  obeys the square-root wall-loading law G_W / Omega_W^2 ∝ sqrt(K).")
