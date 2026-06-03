#!/usr/bin/env python3
"""
Stage 31 SymPy audit.

Checks:
1. Exact tracking-branch critical load and monotonicity in xi.
2. Exact support-enhancement inverse map zeta(S).
3. Exact support-feasibility formulas zeta_req and zeta_crit.
4. Exact stability margins zeta_req < 1/eps and zeta_req < zeta_crit when S_req < S_crit.
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


banner("STAGE 48 — SUPPORT COMPENSATION THEOREM AUDIT")

xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)
eps, zeta = sp.symbols("eps zeta", positive=True, real=True)
Mmix = sp.symbols("M_mix", positive=True, real=True)
Sreq, Scrit = sp.symbols("S_req S_crit", positive=True, real=True)

# Tracking branch.
G_tr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R ** 2) * xi))
F_tr = sp.simplify(
    (9 * delta + (9 + 2 * R ** 2) * xi) ** 2
    * (9 * delta + (9 + 2 * R) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta ** 2 + 18 * delta * xi + (9 + 2 * R ** 2) * xi ** 2) ** 2)
)

banner("1. Tracking branch critical load and monotonicity")
Mcrit = sp.simplify(sp.limit(G_tr, xi, 1, dir="-"))
dG_dxi = sp.factor(sp.diff(G_tr, xi))
print("G_tr =", G_tr)
print("M_crit =", Mcrit)
print("dG_tr/dxi =", dG_dxi)
expect_zero(
    "dG_tr/dxi formula",
    dG_dxi - 9 * (2 * R ** 2 * xi ** 2 + 9 * delta ** 2 + 18 * delta * xi + 9 * xi ** 2)
    / (2 * R ** 2 * xi + 9 * delta + 9 * xi) ** 2,
)
expect_zero("F_tr(xi=0)-1", sp.simplify(F_tr.subs(xi, 0) - 1))
print("limit xi->1^- of F_tr =", sp.limit(F_tr, xi, 1, dir="-"))
soft_coeff = sp.simplify(sp.limit((1 - xi) * F_tr, xi, 1, dir="-"))
print("softening coefficient for F_tr =", soft_coeff)
expect_zero(
    "(1-xi) F_tr softening coefficient",
    soft_coeff
    - (9 * delta + 9 + 2 * R ** 2) ** 2 * (9 * delta + 9 + 2 * R) ** 2
    / (81 * (9 * delta ** 2 + 18 * delta + 9 + 2 * R ** 2) ** 2),
)
Mgap = sp.factor(sp.simplify(Mcrit - G_tr))
print("M_crit - G_tr =", Mgap)
expect_zero(
    "M_crit - G_tr formula",
    Mgap
    - 9 * (1 - xi) * (2 * R ** 2 * xi + 9 * delta ** 2 + 9 * delta * xi + 9 * delta + 9 * xi)
    / ((2 * R ** 2 + 9 * delta + 9) * (2 * R ** 2 * xi + 9 * delta + 9 * xi)),
)

# Support-enhancement factor.
S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))

banner("2. Support-enhancement inverse map")
print("S(zeta;eps) =", S)
expect_zero("S(zeta=0)-1", sp.simplify(S.subs(zeta, 0) - 1))
expect_zero("dS/dzeta - (1-eps)/(1-zeta eps)^2", sp.diff(S, zeta) - (1 - eps) / (1 - zeta * eps) ** 2)
print("limit zeta->(1/eps)^- of S =", sp.limit(S, zeta, 1 / eps, dir="-"))
nu = sp.symbols("nu", positive=True, real=True)
eps_phys = sp.simplify(1 / (1 + nu))
S_phys = sp.simplify(S.subs(eps, eps_phys))
limit_phys = sp.limit(S_phys, zeta, 1 / eps_phys, dir="-")
print("limit zeta->(1/eps)^- of S under 0<eps<1 =", limit_phys)
if limit_phys != sp.oo:
    raise AssertionError("physical-branch support enhancement does not diverge to +oo")

zeta_req = sp.simplify((Sreq - 1) / (1 + eps * (Sreq - 2)))
zeta_crit = sp.simplify((Scrit - 1) / (1 + eps * (Scrit - 2)))

expect_zero("inverse map S(zeta_req)-S_req", sp.simplify(S.subs(zeta, zeta_req) - Sreq))
expect_zero("inverse map S(zeta_crit)-S_crit", sp.simplify(S.subs(zeta, zeta_crit) - Scrit))

banner("3. Exact support-feasibility margins")
print("zeta_req =", zeta_req)
print("zeta_crit =", zeta_crit)

margin_pole = sp.factor(sp.simplify(1 / eps - zeta_req))
margin_branch = sp.factor(sp.simplify(zeta_crit - zeta_req))
print("1/eps - zeta_req =", margin_pole)
print("zeta_crit - zeta_req =", margin_branch)

expect_zero(
    "pole margin formula",
    margin_pole - (1 - eps) / (eps * (1 + eps * (Sreq - 2))),
)
expect_zero(
    "branch margin formula",
    margin_branch - (Scrit - Sreq) * (1 - eps) / ((1 + eps * (Scrit - 2)) * (1 + eps * (Sreq - 2))),
)

banner("4. Exact monotone softening response to support enhancement")
# Implicit derivative using M_mix * S = G_tr.
dxi_dzeta = sp.simplify(Mmix * sp.diff(S, zeta) / sp.diff(G_tr, xi))
print("dxi_phys/dzeta =", sp.factor(dxi_dzeta))
expect_zero(
    "dxi_phys/dzeta formula",
    dxi_dzeta
    - Mmix * (1 - eps) * (2 * R ** 2 * xi + 9 * delta + 9 * xi) ** 2
    / ((1 - zeta * eps) ** 2 * 9 * (2 * R ** 2 * xi ** 2 + 9 * delta ** 2 + 18 * delta * xi + 9 * xi ** 2)),
)

print("\nAll Stage-048 symbolic checks passed.")
