#!/usr/bin/env python3
"""
5pn_stage27_continuum_selected_rank2_closure.py

SymPy audit for Moving-Throat PDE Stage 27.
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


banner("STAGE 27 — CONTINUUM-SELECTED RANK-2 CLOSURE")

# Symbols
xi, delta = sp.symbols("xi delta", real=True)
M_mix, M_supp = sp.symbols("M_mix M_supp", real=True)
R_U, R_phi = sp.symbols("R_U R_phi", real=True)
rho_0, sigma_0, delta_U = sp.symbols("rho_0 sigma_0 delta_U", real=True)
t = sp.symbols("t", real=True)
lambda_0 = sp.Rational(2, 9)

# Exact continuum-selected branch data
q = sp.simplify(t * R_U)
r = sp.simplify(t * R_phi)
print("q =", q)
print("r =", r)
print("lambda_0 =", lambda_0)

# Exact continuum-selected branch equation
n_req = sp.simplify(
    (xi * (delta + xi) - M_mix * (delta + (1 + q**2) * xi))
    / (delta + (1 + r**2) * xi - M_mix * (q - r) ** 2)
)
print("\nn_req(xi,delta;M_mix,q,r) =", n_req)

# Rearrangement to exact quadratic
poly = sp.expand(sp.together(M_supp - n_req).as_numer_denom()[0])
B_cont = sp.simplify(delta - M_mix * (1 + t**2 * R_U**2) - M_supp * (1 + t**2 * R_phi**2))
C_cont = sp.simplify(-delta * (M_mix + M_supp) + t**2 * M_mix * M_supp * (R_U - R_phi) ** 2)
quadratic = sp.expand(xi**2 + B_cont * xi + C_cont)

expect_zero("selected-branch numerator + quadratic", poly + quadratic)
print("\nB_cont =", B_cont)
print("C_cont =", C_cont)

xi_phys = sp.simplify((-B_cont + sp.sqrt(B_cont**2 - 4 * C_cont)) / 2)
print("xi_phys =", xi_phys)
print("xi_phys(M_mix=0, M_supp=0) =", sp.simplify(xi_phys.subs({M_mix: 0, M_supp: 0})))

# Exact continuum-selected normalization function
lambda0 = lambda_0
F_cont = sp.simplify(
    (delta + (1 + lambda0 * R_U * R_phi) * xi) ** 2
    * (delta + (1 + lambda0 * R_phi) * xi - lambda0 * M_mix * (R_U - R_phi) * (R_U - 1)) ** 2
    / (
        (1 - xi)
        * (
            (delta + xi - lambda0 * M_mix * R_U * (R_U - R_phi)) ** 2
            + lambda0 * (M_mix * (R_U - R_phi) + R_phi * xi) ** 2
        ) ** 2
    )
)
D_cont = sp.simplify(
    (delta + xi - lambda0 * M_mix * R_U * (R_U - R_phi)) ** 2
    + lambda0 * (M_mix * (R_U - R_phi) + R_phi * xi) ** 2
)

print("\nF_cont(xi) =", F_cont)
print("D_cont(xi) =", D_cont)

# Exact special surfaces
R_U_expr = sp.simplify((1 + rho_0 / (1 + delta_U)) / (1 + rho_0))
R_phi_expr = sp.simplify((1 + sigma_0 / (1 + delta_U)) / (1 + sigma_0))
Delta_R = sp.simplify(R_U_expr - R_phi_expr)
Delta_R_target = sp.simplify(delta_U * (sigma_0 - rho_0) / ((1 + delta_U) * (1 + rho_0) * (1 + sigma_0)))
expect_zero("Delta_R exact factorization", Delta_R - Delta_R_target)
expect_zero("minimal-kernel source-tied surface (sigma_0=0)", sp.simplify(R_phi_expr.subs(sigma_0, 0) - 1))
expect_zero("interference-matched tracking surface (sigma_0=rho_0)", sp.simplify(R_phi_expr.subs(sigma_0, rho_0) - R_U_expr))

banner("STAGE 27 THEOREM LEDGER")
print("The continuum-selected rank-2 branch is pinned by the exact quadratic")
print("  xi^2 + B_cont xi + C_cont = 0,")
print("with")
print("  B_cont = delta - M_mix (1+t^2 R_U^2) - M_supp (1+t^2 R_phi^2),")
print("  C_cont = - delta (M_mix + M_supp) + t^2 M_mix M_supp (R_U - R_phi)^2.")
print()
print("The exact normalization test is")
print("  R_target = F_cont(xi_phys).")
print()
print("The source-tied surface is sigma_0 = 0, while the tracking surface is sigma_0 = rho_0.")
print("The exact direction mismatch is")
print("  Delta_R = delta_U (sigma_0 - rho_0) / [(1+delta_U)(1+rho_0)(1+sigma_0)].")
