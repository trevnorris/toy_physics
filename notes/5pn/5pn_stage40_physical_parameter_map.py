#!/usr/bin/env python3
"""
5pn_stage40_physical_parameter_map.py

SymPy audit for Moving-Throat PDE Stage 40.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 40 — PHYSICAL (Pe, kappa, eta) PLACEMENT MAP")

# Symbols
Pe, kappa = sp.symbols("Pe kappa", positive=True, real=True)
eta, y = sp.symbols("eta y", positive=True, real=True)
zeta_req = sp.symbols("zeta_req", positive=True, real=True)
pi = sp.pi

# Physical parameter map
Omega_Pe = sp.simplify(
    pi * Pe * (2 * Pe * sp.exp(Pe) + pi)
    / ((4 * Pe**2 + pi**2) * (sp.exp(Pe) - 1))
)
x = sp.simplify(pi**2 / (kappa + pi**2 / 4))
A_K_phys = sp.simplify((kappa + pi**2 / 4) / (kappa + y**2))
A_K_abstract = sp.simplify(1 / (1 - x / 4 + x * y**2 / pi**2))
expect_zero("physical vs abstract A_K", A_K_phys - A_K_abstract)

zeta_phys = sp.simplify(Omega_Pe**2 * A_K_phys)
zeta_max = sp.simplify((pi**2 / 4) * (kappa + pi**2 / 4) / kappa)

print("x(kappa) =", x)
print("A_K(eta;kappa) =", A_K_phys)
print("zeta_0^(Pe+R)(Pe,eta;kappa) =")
sp.pprint(zeta_phys)
print("zeta_max(kappa) =", zeta_max)

# Monotone placement derivatives
expect_zero(
    "partial_kappa zeta - expected",
    sp.simplify(sp.diff(zeta_phys, kappa) - Omega_Pe**2 * (y**2 - pi**2 / 4) / (kappa + y**2) ** 2),
)

dydeta = sp.simplify(1 / (sp.tan(y) + y / sp.cos(y) ** 2))
dA_deta = sp.simplify(sp.diff(A_K_phys, y) * dydeta)
print("\ndy/deta =", dydeta)
print("dA_K/deta =", dA_deta)

expect_zero(
    "closure ceiling zeta_max",
    sp.simplify(sp.limit(zeta_phys.subs(y, 0), Pe, sp.oo) - zeta_max),
)

kappa_max = sp.solve(sp.Eq(zeta_req, zeta_max), kappa)[0]
print("\nkappa_max(zeta_req) =", kappa_max)
expect_zero(
    "zeta_max(kappa_max) - zeta_req",
    sp.simplify(zeta_max.subs(kappa, kappa_max) - zeta_req),
)

# Exact threshold surfaces
Omega_req_sq = sp.simplify(zeta_req * (kappa + y**2) / (kappa + pi**2 / 4))
y_req_sq = sp.simplify((Omega_Pe**2 / zeta_req) * (kappa + pi**2 / 4) - kappa)
kappa_req = sp.simplify((Omega_Pe**2 * pi**2 / 4 - zeta_req * y**2) / (zeta_req - Omega_Pe**2))

print("\nOmega_req^2 =", Omega_req_sq)
print("y_req^2     =", y_req_sq)
print("kappa_req   =", kappa_req)

expect_zero(
    "support equality at Omega_req",
    sp.simplify((Omega_req_sq * A_K_phys) - zeta_req),
)
expect_zero(
    "support equality at y_req",
    sp.simplify((Omega_Pe**2 * (kappa + pi**2 / 4) / (kappa + y_req_sq)) - zeta_req),
)
expect_zero(
    "support equality at kappa_req",
    sp.simplify((Omega_Pe**2 * (kappa_req + pi**2 / 4) / (kappa_req + y**2)) - zeta_req),
)

banner("STAGE 40 THEOREM LEDGER")
print("The explicit lowest-lane family is now written entirely in physical operator variables:")
print("  Pe     = v_sigma L / D_sigma")
print("  kappa  = K_X L^2 / T_X")
print("  eta    = hL = K_m L / T_X")
print()
print("The exact constructive-branch ceiling is")
print("  zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa.")
print("So the Stage-35 demand is reachable on this first physical family iff")
print("  zeta_req <= zeta_max(kappa),")
print("equivalently, whenever zeta_req > pi^2/4,")
print("  kappa <= pi^4 / [4(4 zeta_req - pi^2)].")
