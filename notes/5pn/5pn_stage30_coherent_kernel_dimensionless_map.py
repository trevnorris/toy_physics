#!/usr/bin/env python3
"""
5pn_stage30_coherent_kernel_dimensionless_map.py

SymPy audit for Moving-Throat PDE Stage 30.
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


banner("STAGE 30 — COHERENT-KERNEL DIMENSIONLESS MAP")

# Symbols
chi_0, delta_U = sp.symbols("chi_0 delta_U", positive=True, real=True)
eps_eta, eps_W, zeta = sp.symbols("eps_eta eps_W zeta", real=True)
Z_W, delta_0, Lambda = sp.symbols("Z_W delta_0 Lambda", positive=True, real=True)
pi = sp.pi

# Exact coherent branch data
R_tr = sp.simplify((1 + chi_0 / (1 + delta_U)) / (1 + chi_0))
eps = sp.simplify(eps_W * (1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)))
delta = sp.simplify((delta_0 + eps_eta * delta_U / (1 + delta_U)) / (1 - eps_eta))

M_mix = sp.simplify(8 * Z_W * (1 + chi_0) ** 2 / (pi**2 * (1 - eps_eta) * (1 - eps)))
M_supp = sp.simplify(8 * zeta * Z_W * (1 + chi_0) ** 2 / (pi**2 * (1 - eps_eta) * (1 - zeta * eps)))
R_target = sp.simplify(Lambda * (1 - eps_eta) * (1 - eps) ** 2 / (Z_W * (1 + chi_0) ** 2))

print("R_tr    =", R_tr)
print("eps     =", eps)
print("delta   =", delta)
print("M_mix   =", M_mix)
print("M_supp  =", M_supp)
print("R_target=", R_target)

# Exact support-enhancement factor
S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
M_tr = sp.simplify(M_mix + M_supp)

print("\nS(zeta;eps) =", S)
expect_zero("M_tr - M_mix * S", M_tr - sp.simplify(M_mix * S))
print("dS/dzeta    =", sp.simplify(sp.diff(S, zeta)))
print("S(0;eps)    =", sp.simplify(S.subs(zeta, 0)))
print("S(1;eps)    =", sp.simplify(S.subs(zeta, 1)))

# Exact target product law on coherent branch
expect_zero(
    "R_target * M_tr - 8 Lambda (1-eps)/pi^2 * S",
    sp.simplify(R_target * M_tr - 8 * Lambda * (1 - eps) / pi**2 * S),
)

banner("STAGE 30 THEOREM LEDGER")
print("The coherent local D/N kernel collapses to the exact dimensionless map")
print("  R_tr    = [1 + chi_0/(1+delta_U)] / (1 + chi_0)")
print("  eps     = eps_W [1 - (2/11) delta_U/(1+delta_U)]")
print("  delta   = [delta_0 + eps_eta delta_U/(1+delta_U)] / (1 - eps_eta)")
print("  M_mix   = 8 Z_W (1+chi_0)^2 / [pi^2 (1-eps_eta)(1-eps)]")
print("  M_supp  = 8 zeta Z_W (1+chi_0)^2 / [pi^2 (1-eps_eta)(1-zeta eps)]")
print("  R_target= Lambda (1-eps_eta)(1-eps)^2 / [Z_W (1+chi_0)^2]")
print()
print("The support lane enters only through the exact enhancement factor")
print("  S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps),")
print("with M_tr = M_mix S(zeta;eps).")
print()
print("So the support problem is now a one-parameter enhancement problem rather than a")
print("diffuse multidimensional closure question.")
