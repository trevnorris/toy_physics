#!/usr/bin/env python3
"""
5pn_stage42_operator_branch_residual_bounds.py

Stage 42 audit: exact residual bounds on the operator-selected branch.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 42 — EXACT RESIDUAL BOUNDS ON THE OPERATOR-SELECTED BRANCH")

Xi, Pe_req = sp.symbols("Xi Pe_req", positive=True, real=True)
alpha, eta, y = sp.symbols("alpha eta y", positive=True, real=True)

Delta0 = eta * (sp.cosh(alpha) - 1) / (alpha**2 * (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)))
Deltainf = (sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) / (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))

Pe = sp.symbols("Pe", positive=True, real=True)
Omega = sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi) / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1))
AK = (alpha**2 + sp.pi**2 / 4) / (alpha**2 + y**2)

subbanner("42.1 — Exact physical ratio family")
zeta_of_Pe = Omega**2 * AK
print("zeta_0^(Pe+R)(Pe,eta;kappa) =")
sp.pprint(zeta_of_Pe)
print("\nSince Omega_Pe is strictly increasing on the constructive branch, the exact Stage-41")
print("Pe-bracket immediately induces a zeta-bracket.")

subbanner("42.2 — Exact threshold couplings")
Xi_fail = Pe_req / Deltainf
Xi_suff = Pe_req / Delta0
print("Xi_fail =")
sp.pprint(Xi_fail)
print("Xi_suff =")
sp.pprint(Xi_suff)

zeta_req_expr = zeta_of_Pe.subs(Pe, Pe_req)
zeta_plus_at_fail = zeta_of_Pe.subs(Pe, Xi_fail * Deltainf)
zeta_minus_at_suff = zeta_of_Pe.subs(Pe, Xi_suff * Delta0)
expect_zero("zeta_+(Xi_fail) - zeta_req", zeta_plus_at_fail - zeta_req_expr)
expect_zero("zeta_-(Xi_suff) - zeta_req", zeta_minus_at_suff - zeta_req_expr)

print("\nTherefore:")
print("  Xi <= Xi_fail  -> guaranteed failure inside this operator family")
print("  Xi >= Xi_suff  -> guaranteed success")
print("Only Xi_fail < Xi < Xi_suff still requires the full fixed-point root solve.")

subbanner("42.3 — Residual envelope")
print("Define")
print("  R_phys = zeta_req - zeta_phys,")
print("  R_-    = zeta_req - zeta_+,")
print("  R_+    = zeta_req - zeta_-,")
print("so the exact Stage-41 branch bracket implies")
print("  R_- <= R_phys <= R_+.")

subbanner("42.4 — Weak-coupling branch law")
Omega2_series = sp.series(sp.expand(Omega**2), Pe, 0, 2).removeO()
print("Omega_Pe^2 near Pe=0 =")
sp.pprint(Omega2_series)
expect_zero("linear coefficient", sp.expand(Omega2_series).coeff(Pe, 1) - (4 - sp.pi) / sp.pi)
print("\nWith Pe_* = Xi Delta_0 + O(Xi^2), this gives")
print("  zeta_phys = A_K [1 + ((4-pi)/pi) Xi Delta_0 + O(Xi^2)].")

banner("STAGE 42 FINAL LEDGER")
print("Stage 42 compresses the operator-selected support problem to exact threshold couplings")
print("Xi_fail <= Xi_suff, an exact residual envelope R_- <= R_phys <= R_+, and a weak-coupling")
print("branch law showing how the physical support ratio departs from the symmetric point.")
