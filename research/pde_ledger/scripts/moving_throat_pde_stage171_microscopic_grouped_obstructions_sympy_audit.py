#!/usr/bin/env python3
"""
moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py

SymPy-backed audit for the microscopic decomposition of the linear grouped outlet
obstructions for Stage 171.
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


banner("STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS")

# ---------------------------------------------------------------------------
# 1. Exact decomposition of K_A and G_A.
# ---------------------------------------------------------------------------

P0 = sp.symbols('P0', nonzero=True)
dK, dM = sp.symbols('dK dM')
dB0, dB2 = sp.symbols('dB0 dB2')
dZ0, dZ2 = sp.symbols('dZ0 dZ2')
dN0 = sp.symbols('dN0')

K_exact = sp.expand((-dM - dB2 - dZ2) + sp.Rational(1, 9) * (dK - dB0 - dZ0))
K_split = sp.expand((sp.Rational(1, 9) * dK - dM) - (dB2 + sp.Rational(1, 9) * dB0) - (dZ2 + sp.Rational(1, 9) * dZ0))
expect_zero("K_A split", K_exact - K_split)

G_exact = sp.expand(dN0 - P0 * (dK - dB0 - dZ0))
G_split = sp.expand(-P0 * dK + P0 * dB0 + (dN0 + P0 * dZ0))
expect_zero("G_A split", G_exact - G_split)

# ---------------------------------------------------------------------------
# 2. Exact BdG linearization.
# ---------------------------------------------------------------------------

c, w = sp.symbols('c w', nonzero=True)
dc, dw = sp.symbols('dc dw')
B0 = c**2 / w**2
B2 = c**2 / w**4

dB0_exact = sp.diff(B0, c) * dc + sp.diff(B0, w) * dw
dB2_exact = sp.diff(B2, c) * dc + sp.diff(B2, w) * dw

dB0_formula = 2 * c * dc / w**2 - 2 * c**2 * dw / w**3
dB2_formula = 2 * c * dc / w**4 - 4 * c**2 * dw / w**5
expect_zero("delta B0 formula", dB0_exact - dB0_formula)
expect_zero("delta B2 formula", dB2_exact - dB2_formula)

Bcomb_exact = sp.expand(dB2_exact + sp.Rational(1, 9) * dB0_exact)
Bcomb_formula = sp.expand(
    2 * c * (1 / w**4 + 1 / (9 * w**2)) * dc
    - 2 * c**2 * (2 / w**5 + 1 / (9 * w**3)) * dw
)
expect_zero("BdG obstruction bundle", Bcomb_exact - Bcomb_formula)

# ---------------------------------------------------------------------------
# 3. Exact portwise Maxwell/mixed linearization.
# ---------------------------------------------------------------------------

U, W, R = sp.symbols('U W R', nonzero=True)
gu, gw = sp.symbols('gu gw')
dU, dW, dR = sp.symbols('dU dW dR')
dgu, dgw = sp.symbols('dgu dgw')

Delta_expr = U * W - R**2
S_expr = U + W
Q_expr = gu**2 * W + 2 * gu * gw * R + gw**2 * U
G_expr = gu**2 + gw**2
P_expr = U * gw + R * gu

# Primitive variations.
dDelta_expr = sp.diff(Delta_expr, U) * dU + sp.diff(Delta_expr, W) * dW + sp.diff(Delta_expr, R) * dR
dS_expr = sp.diff(S_expr, U) * dU + sp.diff(S_expr, W) * dW
dQ_expr = (
    sp.diff(Q_expr, U) * dU + sp.diff(Q_expr, W) * dW + sp.diff(Q_expr, R) * dR
    + sp.diff(Q_expr, gu) * dgu + sp.diff(Q_expr, gw) * dgw
)
dG_expr = sp.diff(G_expr, gu) * dgu + sp.diff(G_expr, gw) * dgw
dP_expr = sp.diff(P_expr, U) * dU + sp.diff(P_expr, R) * dR + sp.diff(P_expr, gu) * dgu + sp.diff(P_expr, gw) * dgw

expect_zero("delta Delta formula", dDelta_expr - (W * dU + U * dW - 2 * R * dR))
expect_zero("delta S formula", dS_expr - (dU + dW))
expect_zero("delta G formula", dG_expr - (2 * gu * dgu + 2 * gw * dgw))
expect_zero("delta P formula", dP_expr - (gw * dU + gu * dR + R * dgu + U * dgw))
expect_zero(
    "delta Q formula",
    dQ_expr - (
        gw**2 * dU + gu**2 * dW + 2 * gu * gw * dR
        + 2 * (gu * W + gw * R) * dgu
        + 2 * (gw * U + gu * R) * dgw
    )
)

Delta, S, Q, Gsym, P = sp.symbols('Delta S Q G P', nonzero=True)
dDelta, dS, dQ, dG, dP = sp.symbols('dDelta dS dQ dG dP')
Z0 = Q / Delta
Z2 = (Q * S - Gsym * Delta) / Delta**2
N0expr = P**2 / Delta**2

dZ0_exact = sp.diff(Z0, Q) * dQ + sp.diff(Z0, Delta) * dDelta
dZ2_exact = (
    sp.diff(Z2, Q) * dQ + sp.diff(Z2, S) * dS
    + sp.diff(Z2, Gsym) * dG + sp.diff(Z2, Delta) * dDelta
)
dN0_exact = sp.diff(N0expr, P) * dP + sp.diff(N0expr, Delta) * dDelta

expect_zero("delta Z0 formula", dZ0_exact - (Delta * dQ - Q * dDelta) / Delta**2)
expect_zero(
    "delta Z2 formula",
    dZ2_exact - (
        S * dQ / Delta**2
        + Q * dS / Delta**2
        - dG / Delta
        + (Gsym / Delta**2 - 2 * Q * S / Delta**3) * dDelta
    )
)
expect_zero("delta N0 formula", dN0_exact - (2 * P * dP / Delta**2 - 2 * P**2 * dDelta / Delta**3))

Zcomb_exact = sp.expand(dZ2_exact + sp.Rational(1, 9) * dZ0_exact)
Zcomb_formula = sp.expand(
    (S / Delta**2 + sp.Rational(1, 9) / Delta) * dQ
    + (Q / Delta**2) * dS
    - dG / Delta
    + (Gsym / Delta**2 - Q / (9 * Delta**2) - 2 * Q * S / Delta**3) * dDelta
)
expect_zero("Z obstruction bundle", Zcomb_exact - Zcomb_formula)

Ncomb_exact = sp.expand(dN0_exact + P0 * dZ0_exact)
Ncomb_formula = sp.expand(
    (P0 / Delta) * dQ + (2 * P / Delta**2) * dP
    - (P0 * Q / Delta**2 + 2 * P**2 / Delta**3) * dDelta
)
expect_zero("N obstruction bundle", Ncomb_exact - Ncomb_formula)

# ---------------------------------------------------------------------------
# 4. Weak-axisymmetric reduction.
# ---------------------------------------------------------------------------

eps = sp.symbols('eps')
l20, l21, l22 = sp.Integer(1), sp.Rational(1, 2), sp.Integer(-1)
K1, M1, B01, B21, Z01, Z21, N01 = sp.symbols('K1 M1 B01 B21 Z01 Z21 N01')

for lam in (l20, l21, l22):
    K_micro = eps * lam * (sp.Rational(1, 9) * K1 - M1 - B21 - sp.Rational(1, 9) * B01 - Z21 - sp.Rational(1, 9) * Z01)
    G_micro = eps * lam * (N01 - P0 * K1 + P0 * B01 + P0 * Z01)
    K_rebuilt = eps * lam * ((-M1 - B21 - Z21) + sp.Rational(1, 9) * (K1 - B01 - Z01))
    G_rebuilt = eps * lam * (N01 - P0 * (K1 - B01 - Z01))
    expect_zero(f"weak-axisymmetric K obstruction lambda={lam}", K_micro - K_rebuilt)
    expect_zero(f"weak-axisymmetric G obstruction lambda={lam}", G_micro - G_rebuilt)

print("\nCarry-forward formulas:")
print("  K_A = (1/9 delta K_A - delta M_A) - (delta B_A2 + 1/9 delta B_A0) - (delta Z_A2 + 1/9 delta Z_A0)")
print("  G_A = -P0 delta K_A + P0 delta B_A0 + (delta N_A0 + P0 delta Z_A0)")
print("  The weak grouped branch collapses to the scalar pair (K_1, G_1).")
