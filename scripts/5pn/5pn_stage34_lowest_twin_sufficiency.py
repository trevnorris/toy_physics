#!/usr/bin/env python3
"""
5pn_stage34_lowest_twin_sufficiency.py

SymPy audit for Moving-Throat PDE Stage 34.
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


banner("STAGE 34 — EXACT LOWEST-TWIN SUFFICIENCY CRITERION")

# Symbols
xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)
Lambda, eps = sp.symbols("Lambda eps", positive=True, real=True)
M_mix = sp.symbols("M_mix", positive=True, real=True)
chi0, eps_eta, Z_W = sp.symbols("chi0 eps_eta Z_W", positive=True, real=True)
pi = sp.pi

# Tracking-branch functions
G_tr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R**2) * xi))
F_tr = sp.simplify(
    (9 * delta + (9 + 2 * R**2) * xi) ** 2
    * (9 * delta + (9 + 2 * R) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2)
)
Pi_tr = sp.simplify(F_tr * G_tr)
Pi_tr_expected = sp.simplify(
    xi * (xi + delta) * (9 * delta + (9 + 2 * R) * xi) ** 2 * (9 * delta + (9 + 2 * R**2) * xi)
    / (9 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2)
)

print("G_tr(xi,delta;R) =")
sp.pprint(G_tr)
print("\nF_tr(xi,delta;R) =")
sp.pprint(F_tr)
print("\nPi_tr(xi,delta;R) =")
sp.pprint(Pi_tr)

expect_zero("Pi_tr - closed product form", Pi_tr - Pi_tr_expected)
expect_zero("Pi_tr(0,delta;R)", sp.simplify(Pi_tr.subs(xi, 0)))
print("lim_(xi->1^-) Pi_tr =", sp.simplify(sp.limit(Pi_tr, xi, 1, dir="-")))

# Eliminate zeta_req
C_mix = sp.simplify(8 * Lambda * (1 - eps) / pi**2)
S_req = sp.simplify(Pi_tr / C_mix)
S = sp.symbols("S", positive=True, real=True)
zeta_req = sp.simplify((S_req - 1) / (1 + eps * (S_req - 2)))
zeta_of_S = sp.simplify((S - 1) / (1 + eps * (S - 2)))

print("\nC_mix =", C_mix)
print("S_req =", S_req)
print("zeta_req =", zeta_req)

expect_zero("zeta(S=2) - 1", sp.simplify(zeta_of_S.subs(S, 2) - 1))

Lambda_twin_req = sp.simplify(pi**2 * Pi_tr / (16 * (1 - eps)))
print("\nLambda_twin,req =", Lambda_twin_req)
expect_zero(
    "S_req(Lambda = Lambda_twin,req) - 2",
    sp.simplify(S_req.subs(Lambda, Lambda_twin_req) - 2),
)

# Mixed baseline and wall/mixed overlap threshold
M_mix_twin_req = sp.simplify(G_tr / 2)
print("\nM_mix^(twin,req) =", M_mix_twin_req)

M_mix_formula = sp.simplify(8 * Z_W * (1 + chi0) ** 2 / (pi**2 * (1 - eps_eta) * (1 - eps)))
Z_W_twin_req = sp.simplify(
    pi**2 * (1 - eps_eta) * (1 - eps) * G_tr / (16 * (1 + chi0) ** 2)
)
print("M_mix(Z_W) =", M_mix_formula)
print("Z_W^(twin,req) =", Z_W_twin_req)
expect_zero(
    "M_mix(Z_W_twin,req) - G_tr/2",
    sp.simplify(M_mix_formula.subs(Z_W, Z_W_twin_req) - G_tr / 2),
)

# Exact twin-saturation depth
xi_2x = sp.simplify(
    (
        2 * M_mix * (9 + 2 * R**2)
        - 9 * delta
        + sp.sqrt((2 * M_mix * (9 + 2 * R**2) - 9 * delta) ** 2 + 648 * M_mix * delta)
    )
    / 18
)
print("\nxi_(2x)(M_mix,delta;R) =")
sp.pprint(xi_2x)
expect_zero("G_tr(xi_2x) - 2 M_mix", sp.simplify(G_tr.subs(xi, xi_2x) - 2 * M_mix))

banner("STAGE 34 THEOREM LEDGER")
print("The lowest symmetric twin lane is sufficient iff")
print("  Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1-eps) / pi^2.")
print()
print("Equivalent exact thresholds:")
print("  Lambda_twin,req = pi^2 Pi_tr / [16(1-eps)]")
print("  M_mix^(twin,req) = G_tr / 2")
print("  Z_W^(twin,req) = pi^2 (1-eps_eta)(1-eps) G_tr / [16(1+chi0)^2]")
print()
print("The exact twin-saturation depth at fixed mixed baseline is the quadratic root xi_(2x)")
print("printed above, defined by G_tr(xi_(2x)) = 2 M_mix.")
