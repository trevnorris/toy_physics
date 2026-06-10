#!/usr/bin/env python3
"""
Moving-throat PDE Stage 36 SymPy audit.

Checks:
1. Exact dimensionless support-feasibility function G(xi,delta).
2. Exact monotonicity and endpoint values of G.
3. Exact factorization of the required support loading through G - M_mix.
4. Near-onset asymptotics of the support-feasibility frontier.

Provenance notes:
- `xi`, `delta`, and `M_mix` are the same Stage 035/036 dimensionless support
  variables; this audit keeps the onset boundary `xi >= 0` and the support-gap
  notation unchanged instead of renaming the frontier data.
"""

from __future__ import annotations

import sympy as sp


KAPPA0_SQ = sp.Rational(8) / sp.pi**2
KAPPA1_SQ = sp.Rational(16) / (9 * sp.pi**2)


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


def expect_true(name: str, cond: bool, detail: str) -> None:
    print(f"{name} = {detail}")
    if not cond:
        raise AssertionError(f"{name} is not true")


A, delta = sp.symbols("A delta", positive=True, real=True)
xi = sp.symbols("xi", nonnegative=True, real=True)
Chi, OmegaU, Delta0 = sp.symbols("Chi Omega_U Delta_0", positive=True, real=True)
beta0, NQ = sp.symbols("beta0 NQ", positive=True, real=True)

# Dimensionless support-feasibility variables.
Mmix = sp.symbols("M_mix", positive=True, real=True)
F = sp.simplify((9 * delta + 11 * xi) ** 4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2))
G = sp.simplify(9 * xi * (xi + delta) / (9 * delta + 11 * xi))
alpha_mix = sp.simplify(Chi**2 / (OmegaU**2 * Delta0))
Mmix_expr = sp.simplify(8 * alpha_mix / (sp.pi**2 * A))
R_target = sp.simplify(NQ * A / (beta0 * (sp.Rational(8) / sp.pi**2)))

# Stage 035 eq.(alpha-req): required support amplitude carried into Stage 036.
a_req = sp.simplify(9 * sp.pi**2 * A * xi * (xi + delta) / (8 * (9 * delta + 11 * xi)))
gBreq_sq_over_varpi2 = sp.simplify(a_req - alpha_mix)

banner("STAGE 36.1 — EXACT SUPPORT-FEASIBILITY FUNCTION")
print("G(xi,delta) =", G)
print("F(xi,delta) =", F)
print("M_mix =", Mmix_expr)
print("R_target =", R_target)
# Definitional self-consistency: a_req and G are both hardcoded closed forms,
# so this just confirms a_req = (pi^2 A / 8) G after Mmix_expr cancels.
# The genuine anchor for F (and hence the closed form of G via the same algebra)
# is the symbolic kappa derivation below ("symbolic kappa derivation: F(xi,delta) - R_target_sym").
expect_zero(
    "g_B,req^2/varpi^2 - (pi^2 A / 8) (G - M_mix)",
    gBreq_sq_over_varpi2 - (sp.pi**2 * A / 8) * (G - Mmix_expr),
)

banner("STAGE 36.2 — EXACT MONOTONICITY AND ENDPOINTS OF G")
dG = sp.simplify(sp.diff(G, xi))
dG_target = sp.simplify(9 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) / (9 * delta + 11 * xi) ** 2)
print("dG/dxi =", dG)
expect_zero("dG/dxi - manifestly positive form", dG - dG_target)
expect_zero("G(0,delta)", sp.simplify(G.subs(xi, 0)))
Gmax = sp.simplify(sp.limit(G, xi, 1, dir='-'))
Gmax_target = sp.simplify(9 * (1 + delta) / (9 * delta + 11))
print("G_max(delta) =", Gmax)
expect_zero("G_max - closed form", Gmax - Gmax_target)

banner("STAGE 36.3 — PARAMETRIC FRONTIER AND FINAL ADMISSIBILITY TEST")
xi_req = sp.symbols("xi_req", nonnegative=True, real=True)
print("Parametric frontier: xi -> (F(xi,delta), G(xi,delta))")
# Same definitional identity as above, with xi -> xi_req. Definitional, not load-bearing.
expect_zero(
    "final-test support inequality <-> nonnegative required support loading",
    (sp.pi**2 * A / 8) * (G.subs(xi, xi_req) - Mmix_expr) - gBreq_sq_over_varpi2.subs(xi, xi_req),
)

delta_sample = sp.Integer(1)
xi_sample = sp.Rational(1, 2)
F_sample = sp.simplify(F.subs({delta: delta_sample, xi: xi_sample}))
G_sample = sp.simplify(G.subs({delta: delta_sample, xi: xi_sample}))
Mmix_admissible = sp.N(Mmix_expr.subs({Chi: 1, OmegaU: 1, Delta0: 1, A: 29}))
Mmix_inadmissible = sp.N(Mmix_expr.subs({Chi: 1, OmegaU: 1, Delta0: 1, A: 1}))
G_sample_n = sp.N(G_sample)
expect_true("admissible sample: R_target >= 1", bool(F_sample >= 1), f"R_target={F_sample}")

# Stronger middle-conjunct witness: derive R_target from the Stage-034
# microscopic normalization product N_-(x), not by defining it to equal F.
A_host = sp.Integer(3)
beta0_host = sp.Integer(5)
x_host = sp.simplify(A_host * xi_sample)
deltaK_host = sp.simplify(A_host * delta_sample)
N_host = sp.simplify(
    beta0_host
    * (KAPPA0_SQ * (x_host + deltaK_host) + KAPPA1_SQ * x_host) ** 4
    / (
        KAPPA0_SQ
        * (A_host - x_host)
        * (KAPPA0_SQ * (x_host + deltaK_host) ** 2 + KAPPA1_SQ * x_host**2) ** 2
    )
)
R_target_host = sp.simplify(N_host * A_host / (beta0_host * KAPPA0_SQ))
expect_zero(
    "admissible sample: F(xi_req,delta) - R_target(host)",
    F.subs({delta: delta_sample, xi: xi_sample}) - R_target_host,
)
# Symbolic kappa-based cross-check of the support-feasibility content.
# Builds R_target_sym from the Stage-034 microscopic kappa expansion
# symbolically in (xi, delta, A_sym, beta0_sym) and confirms it equals F.
A_sym, beta0_sym = sp.symbols("A_sym beta0_sym", positive=True, real=True)
x_sym = A_sym * xi
deltaK_sym = A_sym * delta
N_sym = (
    beta0_sym
    * (KAPPA0_SQ * (x_sym + deltaK_sym) + KAPPA1_SQ * x_sym) ** 4
    / (
        KAPPA0_SQ
        * (A_sym - x_sym)
        * (KAPPA0_SQ * (x_sym + deltaK_sym) ** 2 + KAPPA1_SQ * x_sym ** 2) ** 2
    )
)
R_target_sym = sp.simplify(N_sym * A_sym / (beta0_sym * KAPPA0_SQ))
expect_zero(
    "symbolic kappa derivation: F(xi,delta) - R_target_sym",
    sp.simplify(F - R_target_sym),
)
expect_true(
    "admissible sample: M_mix < G(xi_req,delta)",
    bool(Mmix_admissible < G_sample_n),
    f"M_mix={Mmix_admissible}, G={G_sample_n}",
)
expect_true(
    "inadmissible sample: support deficit blocks the branch",
    bool(Mmix_inadmissible > G_sample_n),
    f"M_mix={Mmix_inadmissible}, G={G_sample_n}",
)

banner("STAGE 36.4 — NEAR-ONSET ASYMPTOTICS")
G_series = sp.series(G, xi, 0, 3).removeO()
G_series_target = sp.simplify(xi - 2 * xi**2 / (9 * delta))
print("G(xi,delta) near xi=0 =", G_series)
expect_zero("G near-onset series through O(xi^2)", G_series - G_series_target)

print("All Stage 36 checks passed.")
