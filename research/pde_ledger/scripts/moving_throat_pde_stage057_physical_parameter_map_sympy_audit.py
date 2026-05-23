#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 40 SymPy audit.

Checks:
1. x -> kappa substitution,
2. Robin softening simplification A_K = (kappa + pi^2/4)/(kappa + y^2),
3. exact physical lowest-lane family zeta_0^(Pe+R),
4. monotonic derivative identities in kappa and y,
5. closure ceiling zeta_max(kappa),
6. exact kappa ceiling and parameter-threshold formulas.
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


banner("STAGE 40 — PHYSICAL (Pe, kappa, eta) PLACEMENT MAP")

Pe, kappa, y, zeta_req = sp.symbols("Pe kappa y zeta_req", positive=True, real=True)

Omega_Pe = sp.simplify(
    sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi)
    / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1))
)

x = sp.symbols("x", positive=True, real=True)
A_K_x = sp.simplify(1 / (1 - x / 4 + x * y**2 / sp.pi**2))

x_sub = sp.simplify(sp.pi**2 / (kappa + sp.pi**2 / 4))
A_K_kappa = sp.simplify(A_K_x.subs(x, x_sub))
print("x(kappa) =", x_sub)
print("A_K(eta;kappa) =", A_K_kappa)
expect_zero(
    "A_K - (kappa+pi^2/4)/(kappa+y^2)",
    A_K_kappa - (kappa + sp.pi**2 / 4) / (kappa + y**2),
)

zeta_phys = sp.simplify(Omega_Pe**2 * A_K_kappa)
print("zeta_0^(Pe+R) =", zeta_phys)
expect_zero(
    "zeta_phys - Omega_Pe^2*(kappa+pi^2/4)/(kappa+y^2)",
    zeta_phys - Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y**2),
)

# Monotonic derivatives.
dkappa = sp.simplify(sp.diff(zeta_phys, kappa))
dy = sp.simplify(sp.diff(zeta_phys, y))
print("partial_kappa zeta =", dkappa)
print("partial_y zeta =", dy)
expect_zero(
    "partial_kappa identity",
    dkappa - Omega_Pe**2 * (y**2 - sp.pi**2 / 4) / (kappa + y**2)**2,
)
expect_zero(
    "partial_y identity",
    dy + 2 * Omega_Pe**2 * y * (kappa + sp.pi**2 / 4) / (kappa + y**2)**2,
)

# Constructive-branch closure ceiling.
zeta_max = sp.simplify((sp.pi**2 / 4) * (kappa + sp.pi**2 / 4) / kappa)
print("zeta_max(kappa) =", zeta_max)
expect_zero(
    "zeta_max - limit(Pe->inf, y->0)",
    zeta_max - sp.simplify(sp.limit(sp.limit(zeta_phys, Pe, sp.oo), y, 0, dir="+")),
)

# Exact stiffness ceiling from zeta_req <= zeta_max.
kappa_max = sp.simplify(sp.solve(sp.Eq(zeta_req, zeta_max), kappa)[0])
print("kappa_max(zeta_req) =", kappa_max)
expect_zero(
    "kappa_max identity",
    kappa_max - sp.pi**4 / (4 * (4 * zeta_req - sp.pi**2)),
)

# Exact parameter-threshold formulas.
Omega_req_sq = sp.simplify(zeta_req * (kappa + y**2) / (kappa + sp.pi**2 / 4))
y_req_sq = sp.simplify((Omega_Pe**2 / zeta_req) * (kappa + sp.pi**2 / 4) - kappa)
kappa_req = sp.simplify(
    sp.solve(
        sp.Eq(zeta_req, Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y**2)),
        kappa,
    )[0]
)
print("Omega_req^2 =", Omega_req_sq)
print("y_req^2 =", y_req_sq)
print("kappa_req =", kappa_req)
expect_zero(
    "kappa_req identity",
    kappa_req - (Omega_Pe**2 * sp.pi**2 / 4 - zeta_req * y**2) / (zeta_req - Omega_Pe**2),
)
expect_zero(
    "y_req defining equation",
    zeta_req - sp.simplify(
        Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y_req_sq)
    ),
)

print("\nStage 40 audit passed.")
