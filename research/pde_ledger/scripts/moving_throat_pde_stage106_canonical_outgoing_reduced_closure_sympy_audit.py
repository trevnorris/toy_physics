#!/usr/bin/env python3
"""
Stage 106 SymPy audit: reduced 2.5PN closure on the canonical outgoing DtN branch.

Paper-card `Checks` items (ii) "higher odd terms beyond the point-particle 2.5PN
coefficient" and (iii) "outgoing l=2 DtN fingerprint against the normalized
z = omega a / c_s expansion" are verified upstream. Specifically:
  * Stage 102 (higher_odd_irrelevance) covers item (ii).
  * Stage 104 (outgoing_dtn_fingerprint) covers item (iii) and is the source of
    chi_Q = 1 on the canonical compact passive/outgoing grouped-P_2 branch.
This stage uses chi_Q = 1 as a carry-in and focuses on item (i): factorization
separability of m0hat^2 * chi_Q * N_Q and the closure to N_Q = 1.
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


banner("STAGE 106 — REDUCED 2.5PN CLOSURE ON CANONICAL OUTGOING DtN BRANCH")

G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
m0hat, chi_Q, N_Q = sp.symbols("m0hat chi_Q N_Q", positive=True, real=True)

constraint = sp.Eq(m0hat**2 * chi_Q * N_Q, 1)
print("Normalization constraint =", constraint)

NQ_general = sp.simplify(sp.solve(constraint, N_Q)[0])
print("N_Q on the general outgoing branch =", NQ_general)

NQ_canonical = sp.simplify(NQ_general.subs(chi_Q, 1))
print("N_Q on the canonical outgoing branch =", NQ_canonical)
expect_zero("point-particle canonical branch gives N_Q = 1", NQ_canonical.subs(m0hat, 1) - 1)

K0_target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
K2_target = sp.simplify(6 * G * c_s**3 / (5 * a**3 * c**5))
K4_target = sp.simplify(8 * G * c_s / (15 * a * c**5))
Gamma5_target = sp.simplify(2 * G / (5 * c**5))

K0 = sp.simplify(NQ_general * K0_target)
K2 = sp.simplify(K0 / (4 * (3 * c_s / (2 * a)) ** 2))
K4 = sp.simplify(K0 / (4 * (3 * c_s / (2 * a)) ** 4))
Gamma5 = sp.simplify(NQ_general * Gamma5_target)

print("K0 =", K0)
print("K2 =", K2)
print("K4 =", K4)
print("Gamma5 =", Gamma5)

expect_zero(
    "target identity K0_target K4_target - 4 K2_target^2",
    K0_target * K4_target - 4 * K2_target**2,
)
expect_zero(
    "target identity Gamma5_target - 9 sqrt(K2_target^5 / K0_target^3)",
    Gamma5_target - 9 * sp.sqrt(K2_target**5 / K0_target**3),
)

gamma_eff_canonical = sp.simplify((m0hat**2 * Gamma5).subs(chi_Q, 1))
expect_zero(
    "canonical gamma_eff - target",
    gamma_eff_canonical - Gamma5_target,
)

# F4: Delta_Q first-order sensitivity. On the natural branch (m0hat=1) we have
# N_Q = 1/chi_Q, so gamma_eff = m0hat^2 * N_Q * Gamma5_target collapses to
# Gamma5_target / chi_Q. Expanding around chi_Q = 1 + Delta_Q gives a real
# first-order slope of -Gamma5_target. This exercises the algebra non-trivially
# (a sign flip or factor error in N_Q would change the slope).
Delta_Q = sp.symbols("Delta_Q", real=True)
gamma_eff_off = (m0hat**2 * Gamma5).subs([(m0hat, 1), (chi_Q, 1 + Delta_Q)])
gamma_eff_series = sp.series(gamma_eff_off, Delta_Q, 0, 2).removeO()
linear_coeff = sp.simplify(gamma_eff_series.coeff(Delta_Q, 1))
expect_zero(
    "Delta_Q first-order sensitivity coefficient",
    linear_coeff - (-2 * G / (5 * c**5)),
)
zeroth_coeff = sp.simplify(gamma_eff_series.coeff(Delta_Q, 0))
expect_zero(
    "Delta_Q zeroth-order coefficient equals Gamma5_target",
    zeroth_coeff - Gamma5_target,
)

print("\nRESULT:")
print("  On the canonical outgoing branch, N_Q = 1/m0hat^2 and therefore")
print("  the strict point-particle limit m0hat = 1 collapses the reduced closure to N_Q = 1.")
print("  The effective odd coefficient m0hat^2 Gamma5 then reproduces 2G/(5 c^5) exactly.")
