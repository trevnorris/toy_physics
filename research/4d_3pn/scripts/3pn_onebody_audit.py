#!/usr/bin/env python3
"""
3pn_onebody_audit.py

Kickoff SymPy audit for the 3PN one-body gate.

What this script does
---------------------
1. Expands the exact isotropic Schwarzschild test-mass Lagrangian through 3PN.
2. Extends the carried 2PN denominator-style self sector to cubic order.
3. Shows that a cubic denominator repair alone is not enough at 3PN.
4. Solves the minimal one-body 3PN ledger:
      - cubic static gate mu_rho3,
      - cubic denominator coefficient d3,
      - one extra self slot s24 for U^2 v^4 / c^6.
5. Computes the cubic term predicted by the unextended 2PN invariant denominator,
   and the minimal cubic geometry-invariant correction needed to repair it.

This is not yet a full 3PN paper derivation. It is the exact kickoff audit for the
strict one-body gate.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def show_coeff(name: str, expr: sp.Expr) -> None:
    print(f"{name} = {sp.simplify(sp.expand(expr))}")


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

c, U, v = sp.symbols("c U v", positive=True, real=True)
u = U / c**2

d3, mu_rho3, s24 = sp.symbols("d3 mu_rho3 s24", real=True)
nu = sp.symbols("nu", real=True)


# ---------------------------------------------------------------------------
# Exact isotropic Schwarzschild target
# ---------------------------------------------------------------------------

banner("EXACT ISOTROPIC SCHWARZSCHILD TARGET THROUGH 3PN")

L_exact = -c**2 * sp.sqrt(((1 - u / 2) / (1 + u / 2))**2 - (1 + u / 2)**4 * v**2 / c**2)
L_exact_series = sp.expand(sp.series(L_exact, c, sp.oo, 7).removeO())

print("L_exact/m =")
sp.pprint(L_exact_series)

show_coeff("c^-6 coefficient of v^8", sp.expand(L_exact_series).coeff(c, -6).coeff(v, 8))
show_coeff("c^-6 coefficient of U v^6", sp.expand(L_exact_series).coeff(c, -6).coeff(U, 1).coeff(v, 6))
show_coeff("c^-6 coefficient of U^2 v^4", sp.expand(L_exact_series).coeff(c, -6).coeff(U, 2).coeff(v, 4))
show_coeff("c^-6 coefficient of U^3 v^2", sp.expand(L_exact_series).coeff(c, -6).coeff(U, 3).coeff(v, 2))
show_coeff("c^-6 coefficient of U^4", sp.expand(L_exact_series).coeff(c, -6).coeff(U, 4).coeff(v, 0))


# ---------------------------------------------------------------------------
# Carried 2PN denominator closure extended minimally to cubic order
# ---------------------------------------------------------------------------

banner("CARRIED DENOMINATOR-STYLE SELF SECTOR EXTENDED TO CUBIC ORDER")

D = 1 - 4 * u + 8 * u**2 + d3 * u**3
L_red = -c**2 * (1 - u) * sp.sqrt(1 - (v**2 / c**2) / D)
L_red_series = sp.expand(sp.series(L_red, c, sp.oo, 7).removeO())

print("L_red/m =")
sp.pprint(L_red_series)

show_coeff("red c^-6 coefficient of v^8", sp.expand(L_red_series).coeff(c, -6).coeff(v, 8))
show_coeff("red c^-6 coefficient of U v^6", sp.expand(L_red_series).coeff(c, -6).coeff(U, 1).coeff(v, 6))
show_coeff("red c^-6 coefficient of U^2 v^4", sp.expand(L_red_series).coeff(c, -6).coeff(U, 2).coeff(v, 4))
show_coeff("red c^-6 coefficient of U^3 v^2", sp.expand(L_red_series).coeff(c, -6).coeff(U, 3).coeff(v, 2))

print("\nObservation: the carried cubic denominator reproduces v^8 and U v^6 automatically,")
print("but leaves both the U^3 v^2 slot and the U^2 v^4 slot nontrivial.")


# ---------------------------------------------------------------------------
# Minimal 3PN one-body repair ledger
# ---------------------------------------------------------------------------

banner("MINIMAL 3PN ONE-BODY REPAIR LEDGER")

# Keep the carried lower-order static sector and extend it by a cubic static coefficient.
# Add one explicit new self slot s24 * U^2 v^4 / c^6.
L_candidate = (
    L_red_series
    - U**2 / (2 * c**2)
    + U**3 / (4 * c**4)
    - mu_rho3 * U**4 / (2 * c**6)
    + s24 * U**2 * v**4 / c**6
)

residual = sp.expand(L_exact_series - L_candidate)
print("Exact target minus candidate =")
sp.pprint(residual)

# Solve coefficient-by-coefficient.
sol_mu = sp.solve(sp.Eq(residual.coeff(U, 4).coeff(v, 0).coeff(c, -6), 0), mu_rho3)[0]
sol_d3 = sp.solve(sp.Eq(residual.coeff(U, 3).coeff(v, 2).coeff(c, -6), 0), d3)[0]
sol_s24 = sp.solve(sp.Eq(residual.coeff(U, 2).coeff(v, 4).coeff(c, -6), 0), s24)[0]

show_coeff("mu_rho3", sol_mu)
show_coeff("d3", sol_d3)
show_coeff("s24", sol_s24)

residual_solved = sp.expand(residual.subs({mu_rho3: sol_mu, d3: sol_d3, s24: sol_s24}))
show_coeff("residual after solving", residual_solved)
if sp.simplify(residual_solved) != 0:
    raise AssertionError("Minimal 3PN one-body repair did not match the exact target.")


# ---------------------------------------------------------------------------
# What the unextended 2PN invariant denominator predicts at cubic order
# ---------------------------------------------------------------------------

banner("UNEXTENDED 2PN INVARIANT DENOMINATOR PREDICTION AT CUBIC ORDER")

g1 = sp.Rational(57, 64)
g2 = sp.Rational(298821, 131072)
mu = sp.Rational(32768, 3249)

G_series = 1 + g1 * u + g2 * u**2
D_carry = sp.expand((1 - 4 * u) * (1 + mu * (G_series - 1) ** 2))
D_carry_series = sp.expand(sp.series(D_carry, u, 0, 4).removeO())

print("D_carry(u) =")
sp.pprint(D_carry_series)

show_coeff("carried cubic coefficient d3_carry", D_carry_series.coeff(u, 3))
show_coeff("target cubic coefficient d3_target", sol_d3)


# ---------------------------------------------------------------------------
# Minimal cubic invariant correction if the DtN denominator philosophy is retained
# ---------------------------------------------------------------------------

banner("MINIMAL CUBIC GEOMETRY-INVARIANT CORRECTION")

D_repaired = sp.expand((1 - 4 * u) * (1 + mu * (G_series - 1) ** 2 + nu * (g1 * u) ** 3))
D_repaired_series = sp.expand(sp.series(D_repaired, u, 0, 4).removeO())
print("D_repaired(u) =")
sp.pprint(D_repaired_series)

nu_sol = sp.solve(sp.Eq(D_repaired_series.coeff(u, 3), sol_d3), nu)[0]
show_coeff("nu", nu_sol)

show_coeff("repaired cubic coefficient", sp.simplify(D_repaired_series.coeff(u, 3).subs(nu, nu_sol)))

print("\nImportant: this cubic invariant repair can fix d3, but it still does not fix")
print("the extra one-body self mismatch in the U^2 v^4 slot. That slot requires")
print("one additional 3PN self datum beyond the simple denominator extension.")


# ---------------------------------------------------------------------------
# Final theorem ledger
# ---------------------------------------------------------------------------

banner("FINAL 3PN ONE-BODY KICKOFF LEDGER")
print("Exact isotropic Schwarzschild 3PN target:")
print("  v^8      coefficient = 5/128")
print("  U v^6    coefficient = 11/16")
print("  U^2 v^4  coefficient = 47/16")
print("  U^3 v^2  coefficient = 13/8")
print("  U^4      coefficient = -1/8")
print()
print("Minimal carried-to-exact repair:")
print(f"  mu_rho3 = {sol_mu}")
print(f"  d3      = {sol_d3}")
print(f"  s24     = {sol_s24}")
print()
print("Interpretation:")
print("  - 3PN needs a new cubic static gate (mu_rho3 = 1/4).")
print("  - 3PN needs a new cubic denominator datum (d3 = -45/4).")
print("  - 3PN also opens one genuinely new self slot: U^2 v^4 / c^6 with")
print("    coefficient s24 = -1/16 relative to the carried denominator-style self sector.")
print("  - So 3PN is not just a one-parameter extension of the 2PN one-body closure.")
