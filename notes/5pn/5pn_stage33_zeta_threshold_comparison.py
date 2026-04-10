#!/usr/bin/env python3
"""
5pn_stage33_zeta_threshold_comparison.py

SymPy audit for Moving-Throat PDE Stage 33.
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


banner("STAGE 33 — EXACT COMPARISON OF PHYSICAL zeta AGAINST zeta_req")

# Symbols
eps, zeta_req, S_req = sp.symbols("eps zeta_req S_req", positive=True, real=True)
n, x = sp.symbols("n x", integer=True, positive=True)
pi = sp.pi

# Stage-31 support requirement map
zeta_req_expr = sp.simplify((S_req - 1) / (1 + eps * (S_req - 2)))
S = lambda z: sp.simplify(1 + z * (1 - eps) / (1 - eps * z))

print("zeta_req =", zeta_req_expr)
expect_zero(
    "zeta_req - 1 factorization",
    sp.simplify(zeta_req_expr - 1)
    + (S_req - 2) * (eps - 1) / (S_req * eps - 2 * eps + 1),
)

# Lowest symmetric twin lane
zeta0_twin = sp.Integer(1)
S0 = sp.simplify(S(zeta0_twin))
print("\nzeta_0^(twin) =", zeta0_twin)
print("S_0           =", S0)
expect_zero("lowest twin doubling theorem", S0 - 2)

# Exact equivalence zeta_req <= 1 <=> S_req <= 2 encoded algebraically
expect_zero(
    "zeta_req(S_req=2) - 1",
    sp.simplify(zeta_req_expr.subs(S_req, 2) - 1),
)

# Higher D/N support harmonics
zeta_n_twin = sp.simplify(1 / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1))))
upper_overlap_bound = sp.simplify(1 / (2 * n + 1) ** 2)
print("\nzeta_n^(twin) =", zeta_n_twin)
print("upper overlap-only bound =", upper_overlap_bound)
print("zeta_1 overlap-only ceiling =", upper_overlap_bound.subs(n, 1))
print("zeta_2 overlap-only ceiling =", upper_overlap_bound.subs(n, 2))
print("zeta_3 overlap-only ceiling =", upper_overlap_bound.subs(n, 3))

x_max = sp.simplify((1 / (((2 * n + 1) ** 2) * zeta_req) - 1) / (n * (n + 1)))
print("\nx_max(n; zeta_req) =", x_max)
expect_zero(
    "x_max solves zeta_n^(twin)=zeta_req",
    sp.simplify(zeta_n_twin.subs(x, x_max) - zeta_req),
)

# Exact enhancement ceiling for higher harmonics
S_n_twin = sp.simplify(S(zeta_n_twin))
S_n_max = sp.simplify(1 + (1 - eps) / ((2 * n + 1) ** 2 - eps))
expect_zero(
    "S_n^(max) = S(zeta_max)",
    sp.simplify(S(upper_overlap_bound) - S_n_max),
)

print("\nS_n^(twin) =", S_n_twin)
print("S_n^(max)  =", S_n_max)

banner("STAGE 33 THEOREM LEDGER")
print("The explicit coherent support comparison is now exact.")
print()
print("Lowest symmetric twin lane:")
print("  zeta_0^(twin) = 1,   S_0 = 2.")
print("So the lowest symmetric twin branch succeeds iff")
print("  zeta_req <= 1,   equivalently   S_req <= 2.")
print()
print("Higher D/N harmonics:")
print("  zeta_n^(twin) = 1 / [ (2n+1)^2 (1 + x n(n+1)) ] < 1/(2n+1)^2.")
print("So the exact impossibility tower begins with")
print("  n=1 impossible if zeta_req > 1/9,")
print("  n=2 impossible if zeta_req > 1/25,")
print("  n=3 impossible if zeta_req > 1/49, ...")
print()
print("When a higher harmonic is not ruled out immediately, the exact softness threshold is")
print("  x <= x_max(n; zeta_req).")
print()
print("Thus the explicit finite-throat D/N support tower is strongly biased toward the")
print("lowest coherent twin lane.")
