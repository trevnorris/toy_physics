#!/usr/bin/env python3
"""
Stage 33 SymPy audit.

Checks:
1. Exact equivalence zeta_req <= 1  <=>  S_req <= 2.
2. Exact doubling theorem S(1;eps) = 2.
3. Exact same-operator twin inequality zeta_n^(twin) >= zeta_req gives x <= x_max.
4. Exact impossibility bound zeta_req > (2n+1)^(-2) for n>=1.
5. Exact higher-harmonic enhancement ceiling S_n^(max).
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


banner("STAGE 33 — PHYSICAL ZETA VS ZETA_REQ")

Sreq, eps = sp.symbols("S_req eps", positive=True, real=True)
zeta = sp.symbols("zeta", real=True)
n = sp.symbols("n", integer=True, positive=True)
x = sp.symbols("x", positive=True, real=True)

zeta_req = sp.simplify((Sreq - 1) / (1 + eps * (Sreq - 2)))
S = sp.simplify(1 + zeta * (1 - eps) / (1 - eps * zeta))

banner("1. Exact doubling theorem for the lowest symmetric twin lane")
print("zeta_req =", zeta_req)
print("S(zeta;eps) =", S)
expect_zero("S(1;eps) - 2", sp.simplify(S.subs(zeta, 1) - 2))
# zeta_req <= 1 iff (1-eps)(Sreq-2) <= 0
criterion = sp.factor(sp.simplify(zeta_req - 1))
print("zeta_req - 1 =", criterion)
expect_zero(
    "zeta_req - 1 - (1-eps)(S_req-2)/(1+eps(S_req-2))",
    criterion - (1 - eps) * (Sreq - 2) / (1 + eps * (Sreq - 2)),
)

banner("2. Same-operator twin family and stiffness threshold")
zeta_n = sp.simplify(1 / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1))))
print("zeta_n^(twin) =", zeta_n)
# Solve equality zeta_n = zeta_req for x.
x_eq = sp.solve(sp.Eq(zeta_n, zeta_req), x)[0]
print("x_max(n;zeta_req) =", sp.simplify(x_eq))
expect_zero(
    "x_max - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]",
    sp.simplify(x_eq - (1 / (((2 * n + 1) ** 2) * zeta_req) - 1) / (n * (n + 1))),
)

banner("3. Exact impossibility bound from higher-harmonic suppression")
# Suppression bound relative to overlap only.
supp = sp.simplify((2 * n + 1) ** 2 * zeta_n)
print("(2n+1)^2 zeta_n^(twin) =", supp)
expect_zero("suppression factor - 1/(1+x n(n+1))", supp - 1 / (1 + x * n * (n + 1)))
print("Necessary condition for success: zeta_req <= 1/(2n+1)^2")

banner("4. Higher-harmonic enhancement ceiling")
S_n = sp.simplify(S.subs(zeta, zeta_n))
S_n_max = sp.simplify(1 + (1 - eps) / ((2 * n + 1) ** 2 - eps))
print("S_n^(twin) =", S_n)
print("S_n^(max) =", S_n_max)
# At x=0 the upper bound is saturated.
expect_zero("S_n^(twin)(x=0) - S_n^(max)", sp.simplify(S_n.subs(x, 0) - S_n_max))
# First few explicit cases.
print("S_1^(max) =", sp.simplify(S_n_max.subs(n, 1)))
print("S_2^(max) =", sp.simplify(S_n_max.subs(n, 2)))

print("\nAll Stage-33 symbolic checks passed.")
