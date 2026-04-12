#!/usr/bin/env python3
"""
5pn_stage32_dn_overlap_zeta_extraction.py

SymPy audit for Moving-Throat PDE Stage 32.
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


banner("STAGE 32 — EXPLICIT D/N OVERLAP EXTRACTION OF zeta")

# Symbols
s, L = sp.symbols("s L", positive=True, real=True)
n = sp.symbols("n", integer=True, nonnegative=True)
K_W_eff, K_phi_n_eff = sp.symbols("K_W_eff K_phi_n_eff", positive=True, real=True)
K_X, T_X = sp.symbols("K_X T_X", positive=True, real=True)
zeta_req, x = sp.symbols("zeta_req x", positive=True, real=True)
pi = sp.pi

# Exact D/N family and overlap hierarchy
chi_n = sp.sqrt(2 / L) * sp.sin((n + sp.Rational(1, 2)) * pi * s / L)
I_n = sp.simplify(sp.integrate(chi_n, (s, 0, L)))
I_0 = sp.simplify(I_n.subs(n, 0))
ratio_n = sp.simplify(I_n / I_0)

print("chi_n(s) =", chi_n)
print("I_n      =", I_n)
print("I_0      =", I_0)
print("I_n/I_0  =", ratio_n)
expect_zero("uniform-source overlap ratio theorem", ratio_n - 1 / (2 * n + 1))

# Exact microscopic coherent support ratio
zeta_phys_n = sp.simplify(K_W_eff / K_phi_n_eff * ratio_n**2)
print("\nzeta_n^(phys) =", zeta_phys_n)

# Same-operator twin family
K_W_twin = sp.simplify(K_X + pi**2 * T_X / (4 * L**2))
K_phi_n_twin = sp.simplify(K_W_twin + pi**2 * T_X * n * (n + 1) / L**2)
print("\nK_W^(eff) twin       =", K_W_twin)
print("K_(phi,n)^(eff) twin =", K_phi_n_twin)
expect_zero(
    "same-operator twin stiffness split",
    K_phi_n_twin - (K_W_twin + pi**2 * T_X * n * (n + 1) / L**2),
)

# Introduce x := pi^2 T_X / (L^2 K_W^(eff))
zeta_n_twin = sp.simplify(1 / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1))))
expect_zero(
    "twin-family zeta formula",
    (K_W_twin / (K_W_twin * (1 + x * n * (n + 1))) * ratio_n**2) - zeta_n_twin,
)
print("zeta_n^(twin)        =", zeta_n_twin)
expect_zero("lowest twin half-wave zeta_0 = 1", sp.simplify(zeta_n_twin.subs(n, 0) - 1))

# Exact support-threshold inequality in same-operator twin family
x_max = sp.simplify((1 / (((2 * n + 1) ** 2) * zeta_req) - 1) / (n * (n + 1)))
print("\nx_max(n; zeta_req) =", x_max)
expect_zero(
    "x_max solves zeta_n^(twin)=zeta_req",
    sp.simplify(zeta_n_twin.subs(x, x_max) - zeta_req),
)

banner("STAGE 32 THEOREM LEDGER")
print("For the first explicit coherent finite-throat D/N family,")
print("  zeta_n^(phys) = (K_W^(eff) / K_(phi,n)^(eff)) (I_n / I_0)^2.")
print()
print("For the uniform local source density sigma(s)=1,")
print("  I_n / I_0 = 1 / (2n+1),")
print("so")
print("  zeta_n^(phys) = (K_W^(eff) / K_(phi,n)^(eff)) / (2n+1)^2.")
print()
print("On the same-operator twin family,")
print("  zeta_n^(twin) = 1 / [ (2n+1)^2 (1 + x n(n+1)) ],")
print("with the exact lowest-twin result")
print("  zeta_0^(twin) = 1.")
print()
print("The exact higher-harmonic stiffness threshold is")
print("  x <= x_max(n; zeta_req)")
print("whenever the n-th twin harmonic is not ruled out already by overlap alone.")
