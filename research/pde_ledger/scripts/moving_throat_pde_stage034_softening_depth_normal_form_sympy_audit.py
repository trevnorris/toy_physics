#!/usr/bin/env python3
"""
Moving-throat PDE Stage 17 SymPy audit.

Checks:
1. Exact secular elimination from selected eigenvalue lambda_- to softening depth x = A - lambda_-.
2. Exact formulas for alpha(x), s_-(x), and N_-(x).
3. Exact monotonicity of alpha(x).
4. Exact required support-loading formula g_B^2/varpi^2 in the softening-depth variable.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr):
    simp = sp.simplify(sp.expand(expr))
    print(f"{name} = {simp}")
    if simp != 0:
        raise AssertionError(f"{name} is not zero")


# Generic positive constants for the rank-1 selected-branch model.
ks0, ks1 = sp.symbols("kappa0_sq kappa1_sq", positive=True, real=True)
A, DeltaK = sp.symbols("A DeltaK", positive=True, real=True)
beta0 = sp.symbols("beta0", positive=True, real=True)
lam = sp.symbols("lambda", positive=True, real=True)
x = sp.symbols("x", nonnegative=True, real=True)

banner("STAGE 34.1 — EXACT SOFTENING-DEPTH SECULAR REDUCTION")

# Selected eigenvalue lambda_- lies below the flat branch A.
alpha_lam = sp.simplify(1 / (ks0 / (A - lam) + ks1 / (A + DeltaK - lam)))
S1 = sp.simplify(ks0 / (A - lam) + ks1 / (A + DeltaK - lam))
S1p = sp.simplify(sp.diff(S1, lam))
s_lam = sp.simplify(S1**2 / S1p)
N_lam = sp.simplify(beta0 * s_lam**2 / (ks0 * lam))

# Softening depth x := A - lambda_-.
alpha_x = sp.simplify(x * (x + DeltaK) / (ks0 * (x + DeltaK) + ks1 * x))
s_x = sp.simplify((ks0 * (x + DeltaK) + ks1 * x) ** 2 / (ks0 * (x + DeltaK) ** 2 + ks1 * x**2))
N_x = sp.simplify(beta0 * (ks0 * (x + DeltaK) + ks1 * x) ** 4 / (ks0 * (A - x) * (ks0 * (x + DeltaK) ** 2 + ks1 * x**2) ** 2))

print("alpha(lambda) =")
sp.pprint(alpha_lam)
print("alpha(x) =")
sp.pprint(alpha_x)
print("s_-(x) =")
sp.pprint(s_x)
print("N_-(x) =")
sp.pprint(N_x)

expect_zero("alpha(lambda=A-x) - alpha(x)", alpha_lam.subs(lam, A - x) - alpha_x)
expect_zero("s(lambda=A-x) - s(x)", s_lam.subs(lam, A - x) - s_x)
expect_zero("N(lambda=A-x) - N(x)", N_lam.subs(lam, A - x) - N_x)

banner("STAGE 34.2 — EXACT MONOTONICITY OF THE SOFTENING MAP")
dalpha_dx = sp.simplify(sp.diff(alpha_x, x))
dalpha_target = sp.simplify((ks0 * (x + DeltaK) ** 2 + ks1 * x**2) / (ks0 * (x + DeltaK) + ks1 * x) ** 2)
print("d alpha / dx =")
sp.pprint(dalpha_dx)
expect_zero("dalpha/dx - manifestly positive form", dalpha_dx - dalpha_target)
expect_zero("s_x * d alpha / dx - 1", sp.simplify(s_x * dalpha_dx - 1))

banner("STAGE 34.3 — REQUIRED SUPPORT LOADING IN SOFTENING-DEPTH FORM")
gB_sq_over_varpi2 = sp.symbols("gB_sq_over_varpi2", real=True)
Chi, OmegaU, Delta0 = sp.symbols("Chi Omega_U Delta_0", positive=True, real=True)
alpha_mix = sp.simplify(Chi**2 / (OmegaU**2 * Delta0))
alpha_total = sp.simplify(gB_sq_over_varpi2 + alpha_mix)
gBreq_sq_over_varpi2 = sp.simplify(sp.solve(sp.Eq(alpha_total, alpha_x), gB_sq_over_varpi2)[0])
print("alpha_mix =", alpha_mix)
print("g_B,req^2 / varpi^2 =")
sp.pprint(gBreq_sq_over_varpi2)

# Independent solve of the loading equation in the original lambda variable,
# checked against the x-form solution under lambda = A - x. This is non-trivial:
# it exercises the closed-form alpha_lam together with the substitution.
gBreq_lambda = sp.simplify(
    sp.solve(sp.Eq(gB_sq_over_varpi2 + alpha_mix, alpha_lam), gB_sq_over_varpi2)[0]
)
expect_zero(
    "lambda-form vs x-form support loading",
    gBreq_lambda.subs(lam, A - x) - gBreq_sq_over_varpi2,
)

print("All Stage 17 checks passed.")
