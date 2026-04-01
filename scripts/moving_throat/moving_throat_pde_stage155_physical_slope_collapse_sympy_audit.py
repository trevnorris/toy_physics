#!/usr/bin/env python3
"""
moving_throat_pde_stage155_physical_slope_collapse_sympy_audit.py

SymPy-backed audit for Stage 155.

Checks:
1. Exact first-order variation formulas for the physical grouped variables u_2^(A) and P_0^(A).
2. Exact collapse of the Stage-154 obstruction pair (K_A,G_A) to the physical slopes.
3. Canonical-even collapse K_A = -D_0 delta u_2^(A) at u_2 = 1/9.
4. Equivalence of the microscopic even-consistency relation to delta u_4 = (8/9) delta u_2.
5. Direct outlet formulas in physical grouped variables.
6. Exact identity Delta_Q = delta P_0 / P_0 on the even-preserving branch.
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

banner("STAGE 155 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM")

# ---------------------------------------------------------------------------
# First-order grouped physical variables
# ---------------------------------------------------------------------------

D0, N0, u2 = sp.symbols("D0 N0 u2", real=True, nonzero=True)
dD0, dD2, dD4, dN0 = sp.symbols("dD0 dD2 dD4 dN0", real=True)
eps, lam = sp.symbols("eps lam", real=True)
sigma = sp.symbols("sigma", real=True, nonzero=True)

D2 = -u2 * D0
P0 = N0 / D0

DA0 = D0 + eps * lam * dD0
DA2 = D2 + eps * lam * dD2
NA0 = N0 + eps * lam * dN0

u2A = -DA2 / DA0
P0A = NA0 / DA0

delta_u2 = sp.simplify((sp.series(u2A, eps, 0, 2).removeO() - u2) / (eps * lam))
delta_P0 = sp.simplify((sp.series(P0A, eps, 0, 2).removeO() - P0) / (eps * lam))

print("delta u_2^(A) =", delta_u2)
print("delta P_0^(A) =", delta_P0)

KA = dD2 + sp.Rational(1, 9) * dD0
GA = dN0 - P0 * dD0

banner("Exact obstruction-pair rewrite")
expect_zero(
    "K_A + D0*delta_u2 - (1/9-u2)*dD0",
    KA + D0 * delta_u2 - (sp.Rational(1, 9) - u2) * dD0,
)
expect_zero(
    "G_A - D0*delta_P0",
    GA - D0 * delta_P0,
)

banner("Canonical-even collapse at u2 = 1/9")
expect_zero(
    "K_A + D0*delta_u2 on canonical branch",
    (KA + D0 * delta_u2).subs(u2, sp.Rational(1, 9)),
)

# ---------------------------------------------------------------------------
# Hidden-even consistency relation in physical variables
# ---------------------------------------------------------------------------

banner("Physical form of the hidden-even consistency relation")

u2_star = sp.Rational(1, 9)
u4_star = sp.Rational(4, 81)
D2_star = -u2_star * D0
D4_star = D0 * (u2_star**2 - u4_star)

DA2_star = D2_star + eps * lam * dD2
DA4_star = D4_star + eps * lam * dD4

u2A_star = -DA2_star / DA0
u4A_star = (DA2_star**2 - DA0 * DA4_star) / DA0**2

delta_u2_star = sp.simplify((sp.series(u2A_star, eps, 0, 2).removeO() - u2_star) / (eps * lam))
delta_u4_star = sp.simplify((sp.series(u4A_star, eps, 0, 2).removeO() - u4_star) / (eps * lam))

print("delta u_2^(A) on canonical branch =", delta_u2_star)
print("delta u_4^(A) on canonical branch =", delta_u4_star)

even_relation = sp.Rational(2, 3) * dD2 + sp.Rational(1, 27) * dD0
expect_zero(
    "delta u_4 - (8/9) delta u_2 under microscopic even-consistency relation",
    (delta_u4_star - sp.Rational(8, 9) * delta_u2_star).subs(dD4, even_relation),
)

# ---------------------------------------------------------------------------
# Direct outlet formulas in physical variables
# ---------------------------------------------------------------------------

banner("Direct outlet formulas in physical grouped variables")

delta_kappa = sp.simplify(3 * (1 - sigma) * KA / (sigma * D0))
delta_gamma = sp.simplify(-(1 - sigma) * GA / (9 * sigma * N0))

expect_zero(
    "delta kappa_W^(A) + 3(1-sigma)/sigma * delta u_2^(A)",
    (delta_kappa + 3 * (1 - sigma) * delta_u2 / sigma).subs(u2, sp.Rational(1, 9)),
)

expect_zero(
    "delta gamma_W^(A) + (1-sigma)/(9 sigma P0) * delta P_0^(A)",
    delta_gamma + (1 - sigma) * delta_P0 / (9 * sigma * P0),
)

DeltaQ = sp.simplify(-9 * sigma * delta_gamma / (1 - sigma))
expect_zero(
    "Delta_Q^(A) - delta P_0^(A)/P0",
    DeltaQ - delta_P0 / P0,
)

print("\nCarry-forward formulas:")
print("  K_A = -D_0 delta u_2^(A) on the canonical even branch")
print("  G_A = D_0 delta P_0^(A) exactly")
print("  delta u_4^(A) = (8/9) delta u_2^(A) iff dD_4 = (2/3) dD_2 + (1/27) dD_0")
print("  delta kappa_W^(A) = -3(1-sigma_*)/sigma_* * delta u_2^(A)")
print("  delta gamma_W^(A) = -(1-sigma_*)/(9 sigma_*) * delta P_0^(A)/P_0")
print("  Delta_Q^(A) = delta P_0^(A)/P_0 on the even-preserving branch")
