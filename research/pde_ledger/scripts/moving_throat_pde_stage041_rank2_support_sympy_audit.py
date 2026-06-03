#!/usr/bin/env python3
"""
Moving-throat PDE — Stage 24 SymPy audit.

What this audit verifies
------------------------
1. The selected lower branch for a diagonal 2x2 baseline plus two rank-1
   loadings remains analytically solvable because the determinant is linear
   in the support loading.
2. The exact required support loading n_req(xi,delta;m,q,r) is correct.
3. The mixed-baseline monotonicity theorem dn_req/dm < 0 is exact.
4. If the support tracks the mixed direction (r=q), the rank-2 completion
   collapses exactly to the Stage-23 one-direction geometry.
5. If the support remains tied to the original source direction, the exact
   split-U source-tied support formula follows.
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
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


A0, delta = sp.symbols("A0 delta", positive=True, real=True)
xi = sp.symbols("xi", positive=True, real=True)
m, n = sp.symbols("m n", real=True)
q, r = sp.symbols("q r", real=True)
lam0 = sp.Rational(2, 9)
R_U = sp.symbols("R_U", positive=True, real=True)

banner("STAGE 41 — RANK-2 SUPPORT COMPLETION")

subbanner("24.1 — Exact rank-2 determinant and support-loading theorem")

# Reduced dimensionless loaded matrix after dividing by A0.
M = sp.Matrix([
    [1 - m - n, -(m*q + n*r)],
    [-(m*q + n*r), 1 + delta - m*q**2 - n*r**2],
])

lam = 1 - xi
Det = sp.expand((M - lam * sp.eye(2)).det())
print("det(M - (1-xi)I) =")
sp.pprint(Det)

Det_expected = sp.expand(
    xi * (delta + xi)
    - m * (delta + (1 + q**2) * xi)
    - n * (delta + (1 + r**2) * xi)
    + m * n * (q - r)**2
)
expect_zero("determinant decomposition", Det - Det_expected)

n_req = sp.simplify(sp.solve(sp.Eq(Det_expected, 0), n)[0])
print("n_req =")
sp.pprint(sp.factor(n_req))

n_expected = sp.simplify(
    (xi * (delta + xi) - m * (delta + (1 + q**2) * xi))
    / (delta + (1 + r**2) * xi - m * (q - r)**2)
)
expect_zero("n_req - expected", n_req - n_expected)

subbanner("24.2 — Exact monotonicity with respect to mixed baseline loading")

dn_dm = sp.simplify(sp.diff(n_expected, m))
print("d n_req / d m =")
sp.pprint(sp.factor(dn_dm))

monotone_expected = sp.simplify(
    - (delta + (1 + q * r) * xi)**2
    / (delta + (1 + r**2) * xi - m * (q - r)**2)**2
)
expect_zero("dn/dm - expected", dn_dm - monotone_expected)

subbanner("24.3 — Tracking theorem: support follows the mixed direction")

n_track = sp.simplify(n_expected.subs(r, q))
G_q = sp.simplify(xi * (delta + xi) / (delta + (1 + q**2) * xi))
print("n_req(r=q) =")
sp.pprint(sp.factor(n_track))
expect_zero("tracking collapse", n_track - (G_q - m))

subbanner("24.4 — Source-tied theorem: support remains aligned with the original source vector")

# The physical source-tied specialization uses q = t R_U, r = t with t^2 = lam0.
# Derive n_src by substituting into the general n_expected from section 24.1,
# not by re-stating the substituted formula.
t = sp.symbols("t", real=True)
n_src = sp.simplify(
    n_expected.subs({q: t * R_U, r: t})
)
n_src = sp.simplify(sp.expand(n_src).subs(t**2, lam0))
print("n_req^(src) =")
sp.pprint(sp.factor(n_src))

n_src_expected = sp.simplify(
    (xi * (delta + xi) - m * (delta + (1 + lam0 * R_U**2) * xi))
    / (delta + (1 + lam0) * xi - m * lam0 * (R_U - 1)**2)
)
expect_zero("source-tied formula", n_src - n_src_expected)

# Exact regularity and positivity thresholds.
reg_threshold = sp.simplify((delta + (1 + lam0) * xi) / (lam0 * (R_U - 1)**2))
num_zero_threshold = sp.simplify(xi * (delta + xi) / (delta + (1 + lam0 * R_U**2) * xi))
print("mixed-loading pole threshold =")
sp.pprint(reg_threshold)
print("mixed-loading ceiling from n_req >= 0 =")
sp.pprint(num_zero_threshold)

# Exact derivative in the source-tied case.
dn_dm_src = sp.simplify(sp.diff(n_src_expected, m))
print("d n_req^(src) / d m =")
sp.pprint(sp.factor(dn_dm_src))

dn_dm_src_expected = sp.simplify(
    - (delta + (1 + lam0 * R_U) * xi)**2
    / (delta + (1 + lam0) * xi - m * lam0 * (R_U - 1)**2)**2
)
expect_zero("source-tied dn/dm", dn_dm_src - dn_dm_src_expected)

banner("STAGE 41 THEOREM LEDGER")
print("1. The selected rank-2 wall problem remains analytically solvable because the")
print("   determinant is linear in the support loading n.")
print("2. The exact support-loading theorem is")
print("      n_req = [xi(delta+xi) - m(delta + (1+q^2)xi)]")
print("              / [delta + (1+r^2)xi - m(q-r)^2].")
print("3. The mixed baseline always lowers the support needed to hit the same branch:")
print("      d n_req / d m < 0.")
print("4. If support tracks the mixed direction (r=q), the rank-2 completion collapses")
print("   exactly to the Stage-23 one-direction geometry:")
print("      n_req = G_q - m.")
print("5. If support stays tied to the original source direction, the selected branch")
print("   acquires a genuine new source-tied support-feasibility window.")
