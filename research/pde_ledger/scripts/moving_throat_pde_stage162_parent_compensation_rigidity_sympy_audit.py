#!/usr/bin/env python3
"""
moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py

SymPy-backed audit for Stage 162.

Checks:
1. Exact parent-family similarity identity
       d ln(gamma0) - 2 d ln(L_W/a) = 0.
2. Exact lower-branch differential law
       delta g = (1 - r/(2 sqrt(1+r^2))) delta r.
3. Exact positivity decomposition of the lower-branch slope.
4. Numerical rigidity factor on the Family-1 branch.
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


banner("STAGE 145 — PARENT COMPENSATION-SURFACE RIGIDITY")

r, dr = sp.symbols("r dr", real=True)
pi = sp.pi

# Exact parent family formulas from Stages 99 and 102.
gamma0 = (1 + r**2) / 9
Lratio = pi / 2 * sp.sqrt((1 + r**2) / 3)
g_lower = r - sp.sqrt(1 + r**2) / 2

print("gamma0(r) =", gamma0)
print("L_W/a (r) =", Lratio)
print("g_lower(r) =", g_lower)

# 1. Automatic D/N similarity preservation on the exact parent family.
dlog_gamma = sp.simplify(sp.diff(sp.log(gamma0), r) * dr)
dlog_L = sp.simplify(sp.diff(sp.log(Lratio), r) * dr)
print("d ln gamma0 =", dlog_gamma)
print("d ln(L_W/a) =", dlog_L)
expect_zero("similarity identity", dlog_gamma - 2 * dlog_L)

# 2. Lower-branch differential law.
slope = sp.simplify(sp.diff(g_lower, r))
print("dg_lower/dr =", slope)
expect_zero(
    "lower-branch differential law",
    slope - (1 - r / (2 * sp.sqrt(1 + r**2))),
)

# 3. Exact positivity decomposition.
slope_pos = sp.simplify(
    (4 + 3 * r**2) / (2 * sp.sqrt(1 + r**2) * (2 * sp.sqrt(1 + r**2) + r))
)
expect_zero("positive slope decomposition", slope - slope_pos)
print("The lower-branch slope is manifestly positive for all real r.")

# 4. Numerical Family-1 rigidity factor.
rF1 = sp.Float("1.77799353547498", 30)
slope_num = sp.N(slope.subs(r, rF1), 30)
inv_slope_num = sp.N(1 / slope_num, 30)
print("r_F1 =", rF1)
print("(dg/dr)|_F1 =", slope_num)
print("dr/dg |_F1 =", inv_slope_num)

print("\nCarry-forward formulas:")
print("  On the exact parent family: d ln gamma0 = 2 d ln(L_W/a), so Xi_slip = 0.")
print("  On the lower branch: delta g = (dg/dr) delta r with dg/dr > 0.")
print("  Therefore delta g = 0 implies delta r = 0, hence the first-order")
print("  D/N similarity defect and the reduced 2.5PN normalization defect vanish.")
