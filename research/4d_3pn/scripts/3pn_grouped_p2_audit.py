#!/usr/bin/env python3
"""
3pn_grouped_p2_audit.py

Kickoff SymPy audit for the grouped-P2 conservative gate that should be checked first
on the road to a full 3PN clean derivation.

What this script does
---------------------
1. Defines the grouped real P2 coefficients through O(omega^2) and O(omega^4).
2. Verifies the exact inverse map between grouped coefficients and the axisymmetric
   trace/anisotropy variables (ubar, a, b).
3. States the 3PN first-pass isotropy gate in the raw normalized-support convention.
4. Records the minimal-branch formulas that become available if isotropy passes.

Interpretation:
- 3PN start: extract only u2^(20), u2^(21), u2^(22) and test isotropy.
- 5PN / O(omega^4) later: extract u4^(20), u4^(21), u4^(22) and test the branch identity.
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


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

G, c = sp.symbols("G c", positive=True, real=True)
omega = sp.symbols("omega", real=True)

u220, u221, u222 = sp.symbols("u2_20 u2_21 u2_22", real=True)
u420, u421, u422 = sp.symbols("u4_20 u4_21 u4_22", real=True)

ubar2, a2, b2 = sp.symbols("ubar2 a2 b2", real=True)
ubar4, a4, b4 = sp.symbols("ubar4 a4 b4", real=True)


# ---------------------------------------------------------------------------
# Grouped P2 expansions in normalized-support convention
# ---------------------------------------------------------------------------

banner("GROUPED REAL P2 EXPANSIONS")

Y20 = 1 + u220 * omega**2 + u420 * omega**4
Y21 = 1 + u221 * omega**2 + u421 * omega**4
Y22 = 1 + u222 * omega**2 + u422 * omega**4

print("Y20(omega) =", Y20)
print("Y21(omega) =", Y21)
print("Y22(omega) =", Y22)

print("\n3PN first-pass data:")
print("  u2^(20), u2^(21), u2^(22)")
print("5PN / O(omega^4) follow-up data:")
print("  u4^(20), u4^(21), u4^(22)")


# ---------------------------------------------------------------------------
# Exact grouped -> axisymmetric inverse map
# ---------------------------------------------------------------------------

banner("EXACT GROUPED -> AXISYMMETRIC INVERSE MAP")

ubar2_expr = (u220 + 2 * u221 + 2 * u222) / 5
a2_expr = (2 * u220 - u221 - u222) / 10
b2_expr = (u221 - u222) / 2

ubar4_expr = (u420 + 2 * u421 + 2 * u422) / 5
a4_expr = (2 * u420 - u421 - u422) / 10
b4_expr = (u421 - u422) / 2

print("ubar2 =", ubar2_expr)
print("a2    =", a2_expr)
print("b2    =", b2_expr)
print("ubar4 =", ubar4_expr)
print("a4    =", a4_expr)
print("b4    =", b4_expr)

# Forward map from (ubar, a, b) back to grouped coefficients.
u220_fwd = ubar2 + 4 * a2
u221_fwd = ubar2 - a2 + b2
u222_fwd = ubar2 - a2 - b2

u420_fwd = ubar4 + 4 * a4
u421_fwd = ubar4 - a4 + b4
u422_fwd = ubar4 - a4 - b4

# Verify exact inverse relations.
expect_zero("u2^(20) recovered", u220_fwd.subs({ubar2: ubar2_expr, a2: a2_expr, b2: b2_expr}) - u220)
expect_zero("u2^(21) recovered", u221_fwd.subs({ubar2: ubar2_expr, a2: a2_expr, b2: b2_expr}) - u221)
expect_zero("u2^(22) recovered", u222_fwd.subs({ubar2: ubar2_expr, a2: a2_expr, b2: b2_expr}) - u222)

expect_zero("u4^(20) recovered", u420_fwd.subs({ubar4: ubar4_expr, a4: a4_expr, b4: b4_expr}) - u420)
expect_zero("u4^(21) recovered", u421_fwd.subs({ubar4: ubar4_expr, a4: a4_expr, b4: b4_expr}) - u421)
expect_zero("u4^(22) recovered", u422_fwd.subs({ubar4: ubar4_expr, a4: a4_expr, b4: b4_expr}) - u422)


# ---------------------------------------------------------------------------
# Weighted-sum constraints and anisotropy norms
# ---------------------------------------------------------------------------

banner("WEIGHTED-SUM CONSTRAINTS AND ANISOTROPY NORMS")

expect_zero(
    "weighted sum constraint at O(omega^2)",
    (u220 - ubar2_expr) + 2 * (u221 - ubar2_expr) + 2 * (u222 - ubar2_expr),
)
expect_zero(
    "weighted sum constraint at O(omega^4)",
    (u420 - ubar4_expr) + 2 * (u421 - ubar4_expr) + 2 * (u422 - ubar4_expr),
)

A2_sq = sp.simplify(4 * a2_expr**2 + sp.Rational(4, 5) * b2_expr**2)
A4_sq = sp.simplify(4 * a4_expr**2 + sp.Rational(4, 5) * b4_expr**2)

print("A2^2 =", A2_sq)
print("A4^2 =", A4_sq)

print("\n3PN isotropy gate:")
print("  a2 = 0 and b2 = 0")
print("Equivalently: A2^2 = 0")


# ---------------------------------------------------------------------------
# Minimal isotropic branch formulas
# ---------------------------------------------------------------------------

banner("MINIMAL ISOTROPIC BRANCH FORMULAS")

u2 = sp.symbols("u2", positive=True, real=True)
u4 = sp.symbols("u4", positive=True, real=True)
K0bar = sp.symbols("K0bar", positive=True, real=True)

OmegaQ_sq = sp.simplify(1 / (4 * u2))
Gamma5_norm = sp.simplify(9 * u2**sp.Rational(5, 2))
K0bar_target = sp.simplify(2 * G / (45 * c**5 * u2**sp.Rational(5, 2)))

print("If isotropy passes and the single-pole branch is assumed:")
print("  Omega_Q^2      =", OmegaQ_sq)
print("  Gamma5_norm    =", Gamma5_norm)
print("  K0bar_target   =", K0bar_target)
print()
print("Full 2.5PN closure still requires the 5PN / O(omega^4) branch identity:")
print("  u4 = 4 u2^2")
print()
print("So the clean division of labor is:")
print("  - 3PN: determine (ubar2, a2, b2) and test isotropy.")
print("  - 5PN: determine (ubar4, a4, b4) and test the single-pole identity u4 = 4 u2^2.")


# ---------------------------------------------------------------------------
# Final kickoff ledger
# ---------------------------------------------------------------------------

banner("FINAL GROUPED-P2 KICKOFF LEDGER")
print("Exact grouped trace/anomaly variables:")
print("  ubar2 = (u2^(20) + 2 u2^(21) + 2 u2^(22)) / 5")
print("  a2    = (2 u2^(20) - u2^(21) - u2^(22)) / 10")
print("  b2    = (u2^(21) - u2^(22)) / 2")
print()
print("First 3PN pass/fail test:")
print("  a2 = 0 and b2 = 0")
print()
print("If that test fails, the minimal isotropic quadrupole branch is already ruled out at 3PN.")
print("If it passes, carry ubar2 forward as the candidate pole datum for the 5PN / O(omega^4) test.")
