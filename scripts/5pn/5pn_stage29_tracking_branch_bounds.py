#!/usr/bin/env python3
"""
5pn_stage29_tracking_branch_bounds.py

SymPy audit for Moving-Throat PDE Stage 29.
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


banner("STAGE 29 — TRACKING-BRANCH BOUNDS AND RESIDUAL COMPARISON")

xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)

G_tr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R**2) * xi))
F_tr = sp.simplify(
    (9 * delta + (9 + 2 * R**2) * xi) ** 2
    * (9 * delta + (9 + 2 * R) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2)
)
G_flat = sp.simplify(G_tr.subs(R, 1))
F_flat = sp.simplify(F_tr.subs(R, 1))

print("G_tr   =", G_tr)
print("F_tr   =", F_tr)
print("G_flat =", G_flat)
print("F_flat =", F_flat)

# Exact monotonicity in R

dG = sp.simplify(sp.diff(G_tr, R))
print("\ndG_tr/dR =", dG)
expect_zero(
    "dG_tr/dR theorem",
    dG + 36 * R * xi**2 * (delta + xi) / (2 * R**2 * xi + 9 * delta + 9 * xi) ** 2,
)

# dF/dR factorization
P_R = sp.expand(
    4 * R**4 * xi**3
    + 54 * R**2 * delta**2 * xi + 90 * R**2 * delta * xi**2 + 36 * R**2 * xi**3
    + 162 * R * delta**3 + 324 * R * delta**2 * xi + 162 * R * delta * xi**2
    + 81 * delta**3 + 243 * delta**2 * xi + 243 * delta * xi**2 + 81 * xi**3
)

dF = sp.simplify(sp.diff(F_tr, R))
dF_target = sp.simplify(
    4 * xi * (2 * R * xi + 9 * delta + 9 * xi) * (2 * R**2 * xi + 9 * delta + 9 * xi) * P_R
    / (81 * (1 - xi) * (2 * R**2 * xi**2 + 9 * delta**2 + 18 * delta * xi + 9 * xi**2) ** 3)
)
print("\ndF_tr/dR =", dF)
expect_zero("dF_tr/dR factorization", dF - dF_target)

# Exact comparison with the flat branch
G_excess = sp.simplify(G_tr - G_flat)
G_excess_target = sp.simplify(
    18 * xi**2 * (1 - R**2) * (delta + xi)
    / ((9 * delta + 11 * xi) * (2 * R**2 * xi + 9 * delta + 9 * xi))
)
print("\nG_tr - G_flat =", G_excess)
expect_zero("G_tr - G_flat theorem", G_excess - G_excess_target)

P1 = sp.expand(
    18 * R**2 * delta**2 * xi + 36 * R**2 * delta * xi**2 + 22 * R**2 * xi**3
    + 81 * R * delta**3 + 180 * R * delta**2 * xi + 99 * R * delta * xi**2
    + 162 * delta**3 + 423 * delta**2 * xi + 360 * delta * xi**2 + 99 * xi**3
)
P2 = sp.expand(
    18 * R**3 * delta**2 * xi**2 + 36 * R**3 * delta * xi**3 + 22 * R**3 * xi**4
    + 81 * R**2 * delta**3 * xi + 324 * R**2 * delta**2 * xi**2 + 459 * R**2 * delta * xi**3 + 220 * R**2 * xi**4
    + 81 * R * delta**3 * xi + 243 * R * delta**2 * xi**2 + 261 * R * delta * xi**3 + 99 * R * xi**4
    + 729 * delta**4 + 3078 * delta**3 * xi + 4959 * delta**2 * xi**2 + 3600 * delta * xi**3 + 990 * xi**4
)
F_deficit = sp.simplify(F_flat - F_tr)
F_deficit_target = sp.simplify(
    4 * xi * (1 - R) * P1 * P2
    / (
        81 * (1 - xi)
        * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2
        * (2 * R**2 * xi**2 + 9 * delta**2 + 18 * delta * xi + 9 * xi**2) ** 2
    )
)
print("\nF_flat - F_tr =", F_deficit)
expect_zero("F_flat - F_tr theorem", F_deficit - F_deficit_target)

# Endpoint formulas and bounds
expect_zero("G_tr(R=0) - xi", sp.simplify(G_tr.subs(R, 0) - xi))
expect_zero("F_tr(R=0) - 1/(1-xi)", sp.simplify(F_tr.subs(R, 0) - 1 / (1 - xi)))

# Residual comparison theorem
R_target = sp.symbols("R_target", real=True)
E_tr = sp.simplify(R_target - F_tr)
E_flat = sp.simplify(R_target - F_flat)
expect_zero("E_tr - E_flat - (F_flat - F_tr)", sp.simplify(E_tr - E_flat - (F_flat - F_tr)))

banner("STAGE 29 THEOREM LEDGER")
print("At fixed (xi,delta) the coherent tracking branch is ordered by R:")
print("  dG_tr/dR < 0,   dF_tr/dR > 0  for the physical branch.")
print()
print("Therefore, for 0 < R < 1,")
print("  G_tr > G_flat,   F_tr < F_flat.")
print()
print("The exact endpoint formulas are")
print("  G_tr(R=0) = xi,")
print("  F_tr(R=0) = 1/(1-xi).")
print()
print("So the constructive split-U deformation makes the normalization target harder, not easier:")
print("it raises the required loading and lowers the normalized response relative to the flat branch.")
