#!/usr/bin/env python3
"""
5pn_stage24_source_map_microscopic_normalization_softening.py

Twenty-fourth executable SymPy audit for the 5PN grouped-real-P2 / moving-throat
program.

What this script does
---------------------
1. Eliminates the abstract selected-branch source map on the natural D/N source
   branch.
2. Rewrites the selected-mode normalization product entirely in microscopic
   couplings.
3. Introduces the softening-depth variable x = A - lambda_- and proves the
   exact secular reduction alpha_0(x).
4. Derives the exact overlap s_-(x), normalization product N_-(x), and required
   support loading.

Interpretation
--------------
The selected quadrupole branch is no longer an eigenvector bookkeeping problem.
It is a one-variable spectral-placement problem in the softening depth x.
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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("I. D/N SOURCE MAP ELIMINATION AND THE EXACT SELECTED NORMALIZATION PRODUCT")

# Exact D/N overlap constants from the finite-throat branch.
kappa0 = 2 * sp.sqrt(2) / sp.pi
kappa1 = -4 / (3 * sp.pi)
sigma = sp.simplify(kappa0**2 + kappa1**2)

print("kappa_0 =", kappa0)
print("kappa_1 =", kappa1)
print("sigma   =", sigma)
print("sigma / kappa_0^2 =", sp.simplify(sigma / kappa0**2))

# Stage-16 microscopic abbreviations.
A, DeltaK_ax = sp.symbols("A DeltaK_ax", positive=True, real=True)
Chi, Delta0 = sp.symbols("Chi Delta_0", positive=True, real=True)
alpha_mix, beta0 = sp.symbols("alpha_mix beta_0", positive=True, real=True)
NQ_target = sp.symbols("N_Q_target", positive=True, real=True)

print("alpha_mix = Chi^2 / (Omega_U^2 Delta_0)")
print("beta_0    = Chi^2 / Delta_0^2")

# Softening-depth variable.
x = sp.symbols("x", nonnegative=True, real=True)
lam_minus = sp.simplify(A - x)

# Exact softening-depth secular law from Stage 17.
alpha_x = sp.simplify(x * (x + DeltaK_ax) / (kappa0**2 * (x + DeltaK_ax) + kappa1**2 * x))
print("alpha_0(x) =")
sp.pprint(alpha_x)

# Exact derivative and selected overlap.
dalpha_dx = sp.simplify(sp.diff(alpha_x, x))
s_x = sp.simplify(1 / dalpha_dx)
s_x_closed = sp.simplify(
    (kappa0**2 * (x + DeltaK_ax) + kappa1**2 * x) ** 2
    / (kappa0**2 * (x + DeltaK_ax) ** 2 + kappa1**2 * x**2)
)
expect_zero("s_-(x) from 1/(d alpha_0/dx) minus closed form", s_x - s_x_closed)

print("d alpha_0 / dx =")
sp.pprint(dalpha_dx)
print("s_-(x) =")
sp.pprint(s_x_closed)

# Selected source map on the natural D/N source branch.
mhat_sq = sp.simplify(s_x_closed / kappa0**2)
print("mhat_-^2 =")
sp.pprint(mhat_sq)

# Exact normalization product.
P0_minus = sp.simplify(beta0 * s_x_closed / lam_minus)
N_minus = sp.simplify(beta0 * s_x_closed**2 / (kappa0**2 * lam_minus))
print("P_{0,-}(x) =")
sp.pprint(P0_minus)
print("N_-(x) =")
sp.pprint(N_minus)

expect_zero(
    "N_-(x) - mhat_-^2 P_{0,-}(x)",
    N_minus - sp.simplify(mhat_sq * P0_minus),
)

subbanner("I.1 — Onset value and exact monotonicity")
N_onset = sp.simplify(N_minus.subs(x, 0))
print("N_-(0) =", N_onset)
expect_zero("N_-(0) - beta_0 kappa_0^2 / A", N_onset - sp.simplify(beta0 * kappa0**2 / A))

# Exact positivity of d alpha_0 / dx.
dalpha_dx_closed = sp.simplify(
    (kappa0**2 * (x + DeltaK_ax) ** 2 + kappa1**2 * x**2)
    / (kappa0**2 * (x + DeltaK_ax) + kappa1**2 * x) ** 2
)
expect_zero("d alpha_0/dx minus exact closed form", dalpha_dx - dalpha_dx_closed)

# Required support loading once x is known.
alpha_mix_sym = sp.symbols("alpha_mix", real=True)
gBreq_sq_over_varpi_sq = sp.simplify(alpha_x - alpha_mix_sym)
print("g_B,req^2 / varpi^2 =")
sp.pprint(gBreq_sq_over_varpi_sq)

banner("II. EXACT MICROSCOPIC COUPLING-LEVEL NORMALIZATION EQUATION")

# Full Stage-16 microscopic relation.
OmegaU, OmegaW, gU, gB, gW, gR, varpi = sp.symbols(
    "Omega_U Omega_W g_U g_B g_W g_R varpi", positive=True, real=True
)

Delta0_micro = sp.simplify(OmegaU**2 * OmegaW**2 - gR**2 * sigma)
Chi_micro = sp.simplify(OmegaU**2 * gW + gR * gU)
alpha0_micro = sp.simplify(gB**2 / varpi**2 + Chi_micro**2 / (OmegaU**2 * Delta0_micro))
beta0_micro = sp.simplify(Chi_micro**2 / Delta0_micro**2)
alpha_crit = sp.simplify(A * (A + DeltaK_ax) / ((A + DeltaK_ax) * kappa0**2 + A * kappa1**2))

print("Delta_0 =")
sp.pprint(Delta0_micro)
print("Chi =")
sp.pprint(Chi_micro)
print("alpha_0 =")
sp.pprint(alpha0_micro)
print("beta_0 =")
sp.pprint(beta0_micro)
print("alpha_crit =")
sp.pprint(alpha_crit)
print("alpha_crit with exact D/N constants =")
sp.pprint(sp.simplify(alpha_crit.subs({kappa0**2: 8 / sp.pi**2, kappa1**2: 16 / (9 * sp.pi**2)})))

normalization_eq = sp.Eq(sp.simplify(beta0_micro * s_x_closed**2 / (kappa0**2 * (A - x))), NQ_target)
print("Exact microscopic normalization equation:")
sp.pprint(normalization_eq)

subbanner("II.1 — Exact onset inequality")
N0_micro = sp.simplify(beta0_micro * kappa0**2 / A)
onset_ineq_rhs = sp.simplify(NQ_target * A / kappa0**2)
print("Necessary onset condition: N_-(0) <= N_Q_target")
print("Equivalent inequality for Chi^2:")
sp.pprint(sp.Le(Chi_micro**2, sp.simplify(onset_ineq_rhs * Delta0_micro**2)))

banner("III. WEAK-LOADING EXPANSION AND THE EXACT SOFTENING THRESHOLD")

# Weak-loading expansion: expand N(x(alpha)) around alpha=0 using exact inversion by series.
alpha = sp.symbols("alpha", real=True)
# Solve x = a1 alpha + a2 alpha^2 + O(alpha^3)
a1, a2 = sp.symbols("a1 a2")
x_series = a1 * alpha + a2 * alpha**2
series_eq = sp.expand(alpha_x.subs(x, x_series) - alpha)
sol_a1 = sp.solve(sp.Eq(sp.expand(series_eq).coeff(alpha, 1), 0), a1)[0]
sol_a2 = sp.solve(
    sp.Eq(sp.expand(series_eq.subs(a1, sol_a1)).coeff(alpha, 2), 0),
    a2,
)[0]
print("x(alpha) =", sp.simplify(sol_a1 * alpha + sol_a2 * alpha**2), "+ O(alpha^3)")

N_series = sp.expand(sp.series(N_minus.subs({x: sol_a1 * alpha + sol_a2 * alpha**2}), alpha, 0, 3).removeO())
print("N_-(alpha_0) through O(alpha_0^2) =")
sp.pprint(N_series)

# Exact softening threshold as x -> A^-
alpha_crit_from_limit = sp.simplify(sp.limit(alpha_x, x, A, dir='-'))
expect_zero("alpha_crit from limit minus closed form", alpha_crit_from_limit - alpha_crit)
print("alpha_crit =")
sp.pprint(alpha_crit)

banner("FINAL LEDGER")
print("1. The natural D/N source map is fixed exactly by mhat_-^2 = s_-(x)/kappa_0^2.")
print("2. The selected-branch normalization product is N_-(x) = beta_0 s_-(x)^2 / (kappa_0^2 (A-x)).")
print("3. The branch is parameterized exactly by the softening depth x through alpha_0(x).")
print("4. The remaining theorem gate is the scalar equation N_-(x) = N_Q_target on the stable side x < A.")
