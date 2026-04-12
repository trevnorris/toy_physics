#!/usr/bin/env python3
"""
4pn_onebody_audit.py

Kickoff SymPy audit for the strict 4PN one-body gate.

What this script does
---------------------
1. Expands the exact isotropic Schwarzschild test-mass Lagrangian through 4PN.
2. Continues the carried 3PN denominator-style self sector to quartic order.
3. Shows explicitly that a quartic denominator plus a quartic static gate is still
   not enough at 4PN.
4. Solves the minimal direct continuation of the 3PN one-body ledger:
      - quartic static gate mu_rho4,
      - quartic denominator coefficient d4,
      - two genuinely new self slots s34 and s26,
        corresponding to U^3 v^4 / c^8 and U^2 v^6 / c^8.

This is not yet a full 4PN paper derivation. It is the exact first gate for the
strict one-body sector.
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

d3, d4 = sp.symbols("d3 d4", real=True)
mu_rho3, mu_rho4 = sp.symbols("mu_rho3 mu_rho4", real=True)
s24, s34, s26 = sp.symbols("s24 s34 s26", real=True)


# ---------------------------------------------------------------------------
# Exact isotropic Schwarzschild target through 4PN
# ---------------------------------------------------------------------------

banner("EXACT ISOTROPIC SCHWARZSCHILD TARGET THROUGH 4PN")

L_exact = -c**2 * sp.sqrt(((1 - u / 2) / (1 + u / 2))**2 - (1 + u / 2)**4 * v**2 / c**2)
L_exact_series = sp.expand(sp.series(L_exact, c, sp.oo, 9).removeO())

print("L_exact/m =")
sp.pprint(L_exact_series)

show_coeff("c^-8 coefficient of v^10", sp.expand(L_exact_series).coeff(c, -8).coeff(v, 10))
show_coeff("c^-8 coefficient of U v^8", sp.expand(L_exact_series).coeff(c, -8).coeff(U, 1).coeff(v, 8))
show_coeff("c^-8 coefficient of U^2 v^6", sp.expand(L_exact_series).coeff(c, -8).coeff(U, 2).coeff(v, 6))
show_coeff("c^-8 coefficient of U^3 v^4", sp.expand(L_exact_series).coeff(c, -8).coeff(U, 3).coeff(v, 4))
show_coeff("c^-8 coefficient of U^4 v^2", sp.expand(L_exact_series).coeff(c, -8).coeff(U, 4).coeff(v, 2))
show_coeff("c^-8 coefficient of U^5", sp.expand(L_exact_series).coeff(c, -8).coeff(U, 5).coeff(v, 0))


# ---------------------------------------------------------------------------
# Carried 3PN denominator-style self sector, extended to quartic order
# ---------------------------------------------------------------------------

banner("CARRIED 3PN PACKAGING EXTENDED TO QUARTIC ORDER")

# Direct continuation of the 3PN one-body package.
# Carried lower-order data:
#   d3      = -45/4
#   mu_rho3 = 1/4  via  -mu_rho3 * U^4 / (2 c^6)
#   s24     = -1/16 via +s24 * U^2 v^4 / c^6
D = 1 - 4 * u + 8 * u**2 + d3 * u**3 + d4 * u**4
L_red = -c**2 * (1 - u) * sp.sqrt(1 - (v**2 / c**2) / D)
L_red_series = sp.expand(sp.series(L_red, c, sp.oo, 9).removeO())

L_candidate_min = (
    L_red_series
    - U**2 / (2 * c**2)
    + U**3 / (4 * c**4)
    - mu_rho3 * U**4 / (2 * c**6)
    + mu_rho4 * U**5 / (2 * c**8)
    + s24 * U**2 * v**4 / c**6
)

carried_vals = {
    d3: -sp.Rational(45, 4),
    mu_rho3: sp.Rational(1, 4),
    s24: -sp.Rational(1, 16),
}

residual_min = sp.expand(L_exact_series - L_candidate_min.subs(carried_vals))
print("Residual with only d4 and mu_rho4 still open =")
sp.pprint(residual_min)

show_coeff("residual U^5 / c^8 slot", residual_min.coeff(U, 5).coeff(v, 0).coeff(c, -8))
show_coeff("residual U^4 v^2 / c^8 slot", residual_min.coeff(U, 4).coeff(v, 2).coeff(c, -8))
show_coeff("residual U^3 v^4 / c^8 slot", residual_min.coeff(U, 3).coeff(v, 4).coeff(c, -8))
show_coeff("residual U^2 v^6 / c^8 slot", residual_min.coeff(U, 2).coeff(v, 6).coeff(c, -8))

print("\nObservation:")
print("  A quartic denominator plus a quartic static gate is NOT enough.")
print("  Two genuinely new self-sector slots remain at 4PN:")
print("    U^3 v^4 / c^8 and U^2 v^6 / c^8.")


# ---------------------------------------------------------------------------
# Minimal 4PN one-body repair ledger
# ---------------------------------------------------------------------------

banner("MINIMAL 4PN ONE-BODY REPAIR LEDGER")

L_candidate = (
    L_red_series
    - U**2 / (2 * c**2)
    + U**3 / (4 * c**4)
    - mu_rho3 * U**4 / (2 * c**6)
    + mu_rho4 * U**5 / (2 * c**8)
    + s24 * U**2 * v**4 / c**6
    + s34 * U**3 * v**4 / c**8
    + s26 * U**2 * v**6 / c**8
)

residual = sp.expand(L_exact_series - L_candidate.subs(carried_vals))
print("Exact target minus candidate =")
sp.pprint(residual)

sol_mu4 = sp.solve(sp.Eq(residual.coeff(U, 5).coeff(v, 0).coeff(c, -8), 0), mu_rho4)[0]
sol_d4 = sp.solve(sp.Eq(residual.coeff(U, 4).coeff(v, 2).coeff(c, -8), 0), d4)[0]
sol_s34 = sp.solve(sp.Eq(residual.coeff(U, 3).coeff(v, 4).coeff(c, -8), 0), s34)[0]
sol_s26 = sp.solve(sp.Eq(residual.coeff(U, 2).coeff(v, 6).coeff(c, -8), 0), s26)[0]

show_coeff("mu_rho4", sol_mu4)
show_coeff("d4", sol_d4)
show_coeff("s34", sol_s34)
show_coeff("s26", sol_s26)

residual_solved = sp.expand(
    residual.subs({mu_rho4: sol_mu4, d4: sol_d4, s34: sol_s34, s26: sol_s26})
)
show_coeff("residual after solving", residual_solved)
if sp.simplify(residual_solved) != 0:
    raise AssertionError("Minimal 4PN one-body repair did not match the exact target.")


# ---------------------------------------------------------------------------
# Final theorem ledger
# ---------------------------------------------------------------------------

banner("FINAL 4PN ONE-BODY KICKOFF LEDGER")
print("Exact isotropic Schwarzschild 4PN target:")
print("  v^10      coefficient = 7/256")
print("  U v^8     coefficient = 75/128")
print("  U^2 v^6   coefficient = 59/16")
print("  U^3 v^4   coefficient = 203/32")
print("  U^4 v^2   coefficient = 31/32")
print("  U^5       coefficient = 1/16")
print()
print("Minimal direct continuation of the 3PN one-body package:")
print(f"  mu_rho4 = {sol_mu4}")
print(f"  d4      = {sol_d4}")
print(f"  s34     = {sol_s34}")
print(f"  s26     = {sol_s26}")
print()
print("Interpretation:")
print("  - 4PN needs a new quartic static gate (mu_rho4 = 1/8).")
print("  - 4PN needs a new quartic denominator datum (d4 = 205/16).")
print("  - 4PN also opens two genuinely new self slots:")
print("      s34 for U^3 v^4 / c^8 with coefficient -15/32,")
print("      s26 for U^2 v^6 / c^8 with coefficient -1/16.")
print("  - So, on the clean continuation of the 3PN ansatz, the 4PN one-body repair ledger is")
print("    strictly larger than the 3PN ledger: quartic denominator + quartic static gate are")
print("    still not enough by themselves.")
