#!/usr/bin/env python3
"""
Stage 46 SymPy audit.

Checks:
1. Parent fail/succeed thresholds from G_micro.
2. Coherence-threshold equivalence.
3. Cauchy upper-gain factorization.
4. Exact insertion of G_fail/G_suff = Pe_req/(kappa Delta).
5. Cancellation of explicit K_X in the prefactor after inserting kappa.
"""

from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 80
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 46 — PARENT-OVERLAP THRESHOLD THEOREM")

rho_star, cs_star_sq = sp.symbols("rho_star cs_star_sq", positive=True, real=True)
g_phi, KX, TX, L = sp.symbols("g_phi K_X T_X L", positive=True, real=True)
Nss, Npp, Osp = sp.symbols("N_ss N_pp O_sp", positive=True, real=True)
C2 = sp.symbols("C_sp_sq", positive=True, real=True)
G_fail, G_suff = sp.symbols("G_fail G_suff", positive=True, real=True)
Pe_req, Delta0, Deltainf = sp.symbols("Pe_req Delta0 Deltainf", positive=True, real=True)
kappa = sp.symbols("kappa", positive=True, real=True)

G_micro = sp.simplify(rho_star * g_phi**2 * Osp**2 / (m := sp.symbols('m', positive=True, real=True)) / cs_star_sq / KX / Nss)

# Thresholds in overlap form
# g_phi^2 thresholds
g_fail_sq = sp.simplify(m * cs_star_sq * KX * Nss * G_fail / (rho_star * Osp**2))
g_suff_sq = sp.simplify(m * cs_star_sq * KX * Nss * G_suff / (rho_star * Osp**2))
print("g_(phi,fail)^2 =", g_fail_sq)
print("g_(phi,suff)^2 =", g_suff_sq)

# Coherence-factor substitution
expect_zero(
    "coherence substitution in g_fail^2",
    g_fail_sq.subs(Osp**2, C2 * Nss * Npp)
    - sp.simplify(m * cs_star_sq * KX * G_fail / (rho_star * Npp * C2)),
)
expect_zero(
    "coherence substitution in g_suff^2",
    g_suff_sq.subs(Osp**2, C2 * Nss * Npp)
    - sp.simplify(m * cs_star_sq * KX * G_suff / (rho_star * Npp * C2)),
)

# Coherence thresholds from fixed g_phi
C_fail_sq = sp.simplify(m * cs_star_sq * KX * G_fail / (rho_star * g_phi**2 * Npp))
C_suff_sq = sp.simplify(m * cs_star_sq * KX * G_suff / (rho_star * g_phi**2 * Npp))
print("C_fail^2 =", C_fail_sq)
print("C_suff^2 =", C_suff_sq)
expect_zero("C_suff^2/C_fail^2 - G_suff/G_fail", sp.simplify(C_suff_sq / C_fail_sq - G_suff / G_fail))

# Best-case gain at perfect alignment
G_max = sp.simplify(rho_star * g_phi**2 * Npp / (m * cs_star_sq * KX))
expect_zero(
    "G_micro - G_max*C^2",
    G_micro.subs(Osp**2, C2 * Nss * Npp) - sp.simplify(G_max * C2),
)

# Insert Stage-44 threshold formulas
G_fail_sub = Pe_req / (kappa * Deltainf)
G_suff_sub = Pe_req / (kappa * Delta0)

expect_zero(
    "KX*g_fail threshold with kappa inserted",
    g_fail_sq.subs({G_fail: G_fail_sub, kappa: KX * L**2 / TX})
    - sp.simplify(m * cs_star_sq * TX * Nss * Pe_req / (rho_star * Osp**2 * L**2 * Deltainf)),
)
expect_zero(
    "KX*g_suff threshold with kappa inserted",
    g_suff_sq.subs({G_suff: G_suff_sub, kappa: KX * L**2 / TX})
    - sp.simplify(m * cs_star_sq * TX * Nss * Pe_req / (rho_star * Osp**2 * L**2 * Delta0)),
)

expect_zero(
    "coherence-form fail threshold with kappa inserted",
    g_fail_sq.subs({Osp**2: C2 * Nss * Npp, G_fail: G_fail_sub, kappa: KX * L**2 / TX})
    - sp.simplify(m * cs_star_sq * TX * Pe_req / (rho_star * Npp * C2 * L**2 * Deltainf)),
)
expect_zero(
    "coherence-form suff threshold with kappa inserted",
    g_suff_sq.subs({Osp**2: C2 * Nss * Npp, G_suff: G_suff_sub, kappa: KX * L**2 / TX})
    - sp.simplify(m * cs_star_sq * TX * Pe_req / (rho_star * Npp * C2 * L**2 * Delta0)),
)

print("\nAll Stage 46 symbolic checks passed.")
