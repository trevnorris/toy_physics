#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage012_selected_branch_signature_sympy_audit.py

Stage 012 — compare the Stage-011 numerator/denominator split against the actual
moving-throat selected branch.

What this script does
---------------------
1. Starts from the exact selected-branch normalization product
       N_-(x) = beta0 * s_-(x)^2 / (kappa0^2 * (A - x))
   and the exact overlap formula from the moving-throat selected-mode chain.
2. Introduces the dimensionless softening variables
       xi = x/A, delta = DeltaK_ax/A
   and verifies the universal factorization
       N_-(x) = (beta0/A) * F_num(xi,delta) / (1 - xi).
3. Computes the exact numerator-like and denominator-like log-slopes
       L_num = d_xi ln F_num,
       L_den = d_xi ln[(1-xi)^(-1)].
4. Builds the exact classifier ratio
       R_ND = L_num / L_den,
   whose sign relative to 1 decides whether the selected branch is
   numerator-like or denominator-like.
5. Proves the exact onset threshold delta = 8/9.
6. Proves the exact crossover theorem using a strictly increasing cubic.
7. Reports sample crossover depths xi_*(delta) for representative delta values.

Main structural result
----------------------
The actual selected branch is never literally numerator-rigid or denominator-rigid.
It is an exact co-loading product.

Its signature is nonetheless universal:

- if delta >= 8/9, it is denominator-like on the entire stable branch;
- if 0 < delta < 8/9, it begins numerator-like but crosses exactly once to
  denominator-like;
- near softening it is always denominator-like.
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

def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.factor(sp.simplify(sp.expand(z))))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.factor(sp.simplify(sp.expand(expr)))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")

banner("STAGE 012 — SELECTED-BRANCH NUMERATOR / DENOMINATOR SIGNATURE")

# ---------------------------------------------------------------------------
# Exact selected-branch normalization product
# ---------------------------------------------------------------------------

subbanner("I. Exact selected-branch normalization product")

A, x, DeltaK_ax, beta0 = sp.symbols("A x DeltaK_ax beta_0", positive=True, real=True)
pi = sp.pi
kappa0_sq = sp.Rational(8, 1) / pi**2
kappa1_sq = sp.Rational(16, 9) / pi**2

s_x = sp.simplify(
    (kappa0_sq * (x + DeltaK_ax) + kappa1_sq * x) ** 2
    / (kappa0_sq * (x + DeltaK_ax) ** 2 + kappa1_sq * x**2)
)

N_x = sp.simplify(beta0 * s_x**2 / (kappa0_sq * (A - x)))

print("kappa0^2 =")
sp.pprint(kappa0_sq)
print("kappa1^2 =")
sp.pprint(kappa1_sq)
print("Exact s_-(x) =")
sp.pprint(sp.factor(s_x))
print("Exact N_-(x) =")
sp.pprint(sp.factor(N_x))

# ---------------------------------------------------------------------------
# Dimensionless factorization
# ---------------------------------------------------------------------------

subbanner("II. Dimensionless factorization")

xi, delta = sp.symbols("xi delta", positive=True, real=True)
subs_dim = {
    x: A * xi,
    DeltaK_ax: A * delta,
}

N_dim = sp.simplify((N_x.subs(subs_dim)) / (beta0 / A))

F = sp.simplify(
    (9 * delta + 11 * xi) ** 4
    / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
)

expect_zero("N_-(x)/(beta0/A) - (8/pi^2) F(xi,delta)", sp.simplify(N_dim - (sp.Rational(8,1)/sp.pi**2) * F))

F_num = sp.simplify(
    (9 * delta + 11 * xi) ** 4
    / (81 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
)
F_den = sp.simplify(1 / (1 - xi))
expect_zero("F - F_num * F_den", sp.simplify(F - F_num * F_den))

print("F_num(xi,delta) =")
sp.pprint(F_num)
print("F_den(xi) =")
sp.pprint(F_den)

# ---------------------------------------------------------------------------
# Exact log-slope classifier
# ---------------------------------------------------------------------------

subbanner("III. Exact log-slope classifier")

L_num = sp.factor(sp.simplify(sp.diff(sp.log(F_num), xi)))
L_den = sp.factor(sp.simplify(sp.diff(sp.log(F_den), xi)))
R_ND = sp.factor(sp.simplify(L_num / L_den))

print("L_num = d_xi ln(F_num) =")
sp.pprint(L_num)
print("L_den = d_xi ln(F_den) =")
sp.pprint(L_den)
print("R_ND = L_num / L_den =")
sp.pprint(R_ND)

expect_zero(
    "R_ND - 72 delta^2 (1-xi) / ((9 delta + 11 xi)(9 delta^2 + 18 delta xi + 11 xi^2))",
    sp.simplify(
        R_ND
        - 72 * delta**2 * (1 - xi)
        / ((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2))
    ),
)

print("Onset classifier R_ND(0,delta) =")
sp.pprint(sp.simplify(R_ND.subs({xi: 0})))
print("Near-softening limit lim_{xi->1^-} R_ND =")
sp.pprint(sp.simplify(sp.limit(R_ND, xi, 1, dir="-")))

# ---------------------------------------------------------------------------
# Exact crossover theorem
# ---------------------------------------------------------------------------

subbanner("IV. Exact crossover theorem")

P_poly = sp.expand((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) - 72 * delta**2 * (1 - xi))
dP_dxi = sp.factor(sp.diff(P_poly, xi))

print("Crossover polynomial P(xi,delta) =")
sp.pprint(P_poly)
print("dP/dxi =")
sp.pprint(dP_dxi)

expect_zero(
    "P(0,delta) - 9 delta^2 (9 delta - 8)",
    sp.simplify(P_poly.subs({xi: 0}) - 9 * delta**2 * (9 * delta - 8)),
)

print("Strict positivity of dP/dxi for xi>=0, delta>0 is manifest.")
print("Therefore P(xi,delta) is strictly increasing in xi.")

print("\nTheorem split:")
print("  - if delta >= 8/9, then P(0,delta) >= 0 and P stays positive,")
print("    so the branch is denominator-like for the whole stable interval;")
print("  - if 0 < delta < 8/9, then P(0,delta) < 0 but P(1,delta) > 0,")
print("    so there is a unique crossover xi_*(delta) in (0,1).")

# Verify P(1,delta) > 0 explicitly
P_at_1 = sp.factor(sp.simplify(P_poly.subs({xi: 1})))
print("P(1,delta) =")
sp.pprint(P_at_1)

# ---------------------------------------------------------------------------
# Sample crossover depths
# ---------------------------------------------------------------------------

subbanner("V. Sample crossover depths")

def sample_root(delta_value: sp.Rational):
    poly = sp.Poly(sp.expand(P_poly.subs({delta: delta_value})), xi)
    roots = [complex(r) for r in sp.nroots(poly)]
    real_roots = sorted(float(r.real) for r in roots if abs(r.imag) < 1e-12 and 0.0 < r.real < 1.0)
    if len(real_roots) != 1:
        raise AssertionError(f"Expected exactly one physical crossover root for delta={delta_value}, got {real_roots}")
    return real_roots[0]

for dval in [sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4)]:
    root = sample_root(dval)
    print(f"delta = {sp.nsimplify(dval)}  ->  xi_* ≈ {root:.15f}")

subbanner("VI. Final verdict")
print("The selected branch is an exact numerator/denominator co-loading product.")
print("It is never literally numerator-rigid or denominator-rigid.")
print("But its signature is universal:")
print("  * delta >= 8/9      -> denominator-like everywhere")
print("  * 0 < delta < 8/9   -> numerator-like near onset, unique crossover, denominator-like later")
print("  * xi -> 1^-         -> always denominator-like")
