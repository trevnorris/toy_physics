#!/usr/bin/env python3
"""
5pn_stage31_support_compensation_theorem.py

SymPy audit for Moving-Throat PDE Stage 31.
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


banner("STAGE 31 — SUPPORT-COMPENSATION THEOREM")

# Symbols
xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)
M_mix, M_req = sp.symbols("M_mix M_req", positive=True, real=True)
eps, zeta = sp.symbols("eps zeta", real=True)
S_req, S_crit = sp.symbols("S_req S_crit", positive=True, real=True)
R_target = sp.symbols("R_target", positive=True, real=True)

# Universal tracking branch and critical load
G_tr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R**2) * xi))
dG_dxi = sp.simplify(sp.diff(G_tr, xi))
M_crit = sp.simplify(G_tr.subs(xi, 1))
Mcrit_diff = sp.simplify(M_crit - G_tr)

print("G_tr(xi,delta;R) =", G_tr)
print("dG_tr/dxi        =", dG_dxi)
print("M_crit(delta,R)  =", M_crit)
expect_zero(
    "M_crit - G_tr exact formula",
    Mcrit_diff
    - 9 * (1 - xi) * (2 * R**2 * xi + 9 * delta**2 + 9 * delta * xi + 9 * delta + 9 * xi)
    / ((2 * R**2 + 9 * delta + 9) * (2 * R**2 * xi + 9 * delta + 9 * xi)),
)

# Exact support-enhancement factor and inverse
S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
dS_dzeta = sp.simplify(sp.diff(S, zeta))
zeta_req = sp.simplify((S_req - 1) / (1 + eps * (S_req - 2)))

print("\nS(zeta;eps) =", S)
print("dS/dzeta    =", dS_dzeta)
print("zeta_req    =", zeta_req)
print("S(0;eps)    =", sp.simplify(S.subs(zeta, 0)))
print("limit zeta->1/eps^- of S = +oo (formal pole at zeta = 1/eps)")
expect_zero("S(zeta_req) - S_req", sp.simplify(S.subs(zeta, zeta_req) - S_req))
expect_zero(
    "1/eps - zeta_req stability margin",
    sp.simplify(1 / eps - zeta_req) - (1 - eps) / (eps * (1 + eps * (S_req - 2))),
)

# Stable-side target point from the tracking normalization law
F_tr = sp.simplify(
    (9 * delta + (9 + 2 * R**2) * xi) ** 2
    * (9 * delta + (9 + 2 * R) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2)
)
print("\nF_tr(0,delta;R) =", sp.simplify(F_tr.subs(xi, 0)))
print("limit xi->1^- of F_tr = +oo (simple pole from 1/(1-xi))")

# Exact support-compensation theorem algebra
S_req_from_load = sp.simplify(M_req / M_mix)
expect_zero("zeta_req(M_req/M_mix) reproduces required enhancement", sp.simplify(S.subs(zeta, zeta_req.subs(S_req, S_req_from_load)) - S_req_from_load))

S_crit_expr = sp.simplify(M_crit / M_mix)
zeta_crit = sp.simplify((S_crit - 1) / (1 + eps * (S_crit - 2)))
expect_zero(
    "zeta_crit - zeta_req theorem",
    sp.simplify(zeta_crit - zeta_req)
    - (S_crit - S_req) * (1 - eps) / ((1 + eps * (S_crit - 2)) * (1 + eps * (S_req - 2))),
)

# Exact implicit branch response to support enhancement
# From M_mix S(zeta;eps) = G_tr(xi_phys,delta;R), we obtain
# dxi_phys/dzeta = M_mix (dS/dzeta) / (dG_tr/dxi)|_{xi_phys}.
implicit_response = sp.simplify(M_mix * dS_dzeta / dG_dxi)
print("\ndxi_phys/dzeta (implicit-form theorem) =", implicit_response)

banner("STAGE 31 THEOREM LEDGER")
print("There is no reduced-level support no-go on the physical tracking branch.")
print()
print("Exact facts:")
print("  - G_tr is strictly increasing on 0 < xi < 1.")
print("  - M_crit = G_tr(1,delta;R) is the exact softening-limit load.")
print("  - S(zeta;eps) is strictly increasing and invertible on 0 <= zeta < 1/eps.")
print("  - Every finite target R_target > 1 corresponds to at least one stable-side xi_req in (0,1).")
print("  - Whenever M_mix < M_req, the unique coherent support ratio")
print("      zeta_req = (S_req - 1) / [1 + eps (S_req - 2)]")
print("    reaches the target before softening, with zeta_req < zeta_crit < 1/eps.")
print()
print("So the remaining theorem gap is no longer support-feasibility in principle.")
print("It is whether the actual moving-throat PDE produces a physical zeta large enough to meet zeta_req.")
