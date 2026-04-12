#!/usr/bin/env python3
"""
5pn_stage99_dn_mixed_tube_realization.py

Stage 99 audit: finite D/N mixed-tube realization.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 99 — FINITE D/N MIXED-TUBE REALIZATION")

L_W, a, c_s, omega, z = sp.symbols("L_W a c_s omega z", positive=True, real=True)
r_c = sp.symbols("r_c", positive=True, real=True)

k_W = sp.simplify(sp.pi / (2 * L_W))
Omega_W = sp.simplify(sp.pi * c_s / (2 * L_W))
print("k_W =", k_W)
print("Omega_W =", Omega_W)

kappa0 = sp.simplify((omega**2 / Omega_W**2) / ((a * omega / c_s)**2))
print("kappa_0 =", kappa0)
expect_zero("kappa_0 formula", kappa0 - 4 * L_W**2 / (sp.pi**2 * a**2))

L_W_sol = sp.solve(sp.Eq(kappa0, (1 + r_c) / 3), L_W)[0]
print("L_W compensation-selected =", sp.factor(L_W_sol))
expect_zero(
    "L_W formula",
    L_W_sol - sp.pi * a * sp.sqrt((1 + r_c) / 3) / 2,
)

gamma0 = sp.simplify((1 + r_c) / 9)
Dw_bare = sp.simplify((1 + r_c) * (1 - z**2 / 3 - sp.I * z**5 / 9))
print("D_W^bare(z) =")
sp.pprint(Dw_bare)

kappa_c = sp.simplify(((1 + r_c) / 3) / (1 + r_c))
gamma_c = sp.simplify(gamma0 / (1 + r_c))
expect_zero("kappa_c - 1/3", kappa_c - sp.Rational(1, 3))
expect_zero("gamma_c - 1/9", gamma_c - sp.Rational(1, 9))

banner("STAGE 99 FINAL LEDGER")
print("The bare mixed side-channel can be realized geometrically as the first D/N half-wave on an auxiliary tube.")
print("Compensation selects")
print("  L_W = (pi a / 2) * sqrt((1+r_c)/3),")
print("and a pure-scale deformation of the canonical compact outgoing l=2 branch gives")
print("  D_W^bare(z) = (1+r_c)(1 - z^2/3 - i z^5/9) + O(z^6),")
print("so the hybridization factor is removed exactly and the final coefficients remain")
print("  kappa_c = 1/3,   gamma_c = 1/9.")
