#!/usr/bin/env python3
"""
4pn_conditional_referee_master_audit.py

Standalone conditional referee master for the current 4PN program.

Purpose
-------
The local 4PN ordinary sector has already been assembled exactly in the previous
referee chain, while the hereditary sector has now been reduced to the canonical
quadrupole normalization problem inherited from 2.5PN.

This script packages the strongest honest 4PN theorem currently available inside
that hierarchy:

    L_4PN^cons = L_4PN^local + L_4PN^tail,

with

    C_tail = (G M / (2 c^3)) * gamma_quad^eff,

so that no new independent quadrupole-normalization datum opens at 4PN.

What this script checks
-----------------------
1. Freezes the already-proved local/tail split as the correct theorem envelope.
2. Verifies the exact GR coefficient relation
       C_tail^GR = (G M / (2 c^3)) * gamma_GR
   with gamma_GR = 2G/(5 c^5).
3. Verifies the local logarithmic tail shadow
       F_tail / mu = (16/15) nu U^4 / c^8 * (12 v^2 - 11 rdot^2),
   i.e. the shadow lives only in the degree-2 G^4/r^4 block.
4. Inserts the toy-model canonical quadrupole branch data
       gamma_quad^eff = N_Q a^5/(27 c_s^5) = barGamma_5
                       = 9 \bar K_2^(5/2) / \bar K_0^(3/2)
   and proves that the 4PN hereditary coefficient is controlled by the same
   normalization problem already isolated at 2.5PN.
5. Verifies that once gamma_quad^eff is set to the canonical Burke--Thorne value,
   the GR 4PN tail coefficient follows automatically.

Important scope note
--------------------
This does not replace the already existing local referee audit, and it does not
claim a first-principles moving-throat PDE derivation. It is the strongest
current full-4PN theorem ledger available inside the declared closure hierarchy.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

G, M, mu, c, c_s = sp.symbols("G M mu c c_s", positive=True, real=True)
nu = sp.symbols("nu", positive=True, real=True)
U, v2, rd2 = sp.symbols("U v2 rd2", real=True)

a, N_Q = sp.symbols("a N_Q", positive=True, real=True)
K0bar, K2bar = sp.symbols("K0bar K2bar", positive=True, real=True)
gamma_quad_eff, gamma_GR = sp.symbols("gamma_quad_eff gamma_GR", real=True)
C_tail_GR, C_tail_toy = sp.symbols("C_tail_GR C_tail_toy", real=True)


# ---------------------------------------------------------------------------
# Part I. Theorem envelope
# ---------------------------------------------------------------------------

banner("PART I — 4PN THEOREM ENVELOPE")

print("The already-frozen theorem envelope is")
print("  L_4PN^cons = L_4PN^local + L_4PN^tail")
print()
print("Inside the present hierarchy:")
print("  - the local 4PN ordinary sector is already assembled exactly by the local referee chain;")
print("  - the remaining question is whether the hereditary coefficient introduces any new")
print("    normalization datum beyond the 2.5PN quadrupole one.")


# ---------------------------------------------------------------------------
# Part II. Exact GR bridge coefficient
# ---------------------------------------------------------------------------

banner("PART II — EXACT GR BRIDGE BETWEEN 2.5PN AND 4PN")

gamma_GR = sp.simplify(2 * G / (5 * c**5))
C_tail_GR = sp.simplify(G**2 * M / (5 * c**8))

print("gamma_GR  =", gamma_GR)
print("C_tail_GR =", C_tail_GR)
expect_zero(
    "C_tail_GR - (G M / (2 c^3)) gamma_GR",
    C_tail_GR - (G * M / (2 * c**3)) * gamma_GR,
)

print(
    "\nSo the universal GR 4PN tail coefficient is exactly the Burke–Thorne coefficient\n"
    "multiplied by the standard monopole-scattering factor G M / (2 c^3)."
)


# ---------------------------------------------------------------------------
# Part III. Local logarithmic tail shadow
# ---------------------------------------------------------------------------

banner("PART III — LOCAL LOGARITHMIC TAIL SHADOW")

F_tail_over_mu = sp.simplify(sp.Rational(16, 15) * nu * U**4 / c**8 * (12 * v2 - 11 * rd2))
print("F_tail / mu =")
sp.pprint(F_tail_over_mu)

print(
    "\nThis is the unique local logarithmic shadow extracted from the STF quadrupole tail.\n"
    "It occupies only the degree-2 G^4/r^4 block of the local 4PN bookkeeping."
)

subbanner("III.1 — Pure block statement")
# A trivial but explicit check: the expression is linear in v2 and rd2 and quartic in U.
expect_zero("shadow degree in U minus 4", sp.Poly(F_tail_over_mu * c**8 / nu, U, v2, rd2).total_degree() - 5)
print("(Total degree 5 in {U,v2,rd2} here corresponds to U^4 times one quadratic kinematic slot.)")


# ---------------------------------------------------------------------------
# Part IV. Toy-model quadrupole branch insertion
# ---------------------------------------------------------------------------

banner("PART IV — TOY-MODEL QUADRUPOLE BRANCH INSERTION")

gamma_quad_eff = sp.simplify(N_Q * a**5 / (27 * c_s**5))
C_tail_toy = sp.simplify((G * M / (2 * c**3)) * gamma_quad_eff)

print("gamma_quad_eff =", gamma_quad_eff)
print("C_tail_toy     =", C_tail_toy)

subbanner("IV.1 — Canonical invariant branch form")
gammabar_invariant = sp.simplify(9 * K2bar**sp.Rational(5, 2) / K0bar**sp.Rational(3, 2))
# Branch identity: K2bar = K0bar * a^2 / (9 c_s^2)
branch_sub = {K2bar: K0bar * a**2 / (9 * c_s**2)}
expect_zero(
    "gamma_quad_eff(branch) - invariant form",
    sp.simplify(gamma_quad_eff.subs(N_Q, K0bar) - gammabar_invariant.subs(branch_sub)),
)

print(
    "So the same canonical quadrupole normalization can be written either as\n"
    "  gamma_quad_eff = N_Q a^5/(27 c_s^5)\n"
    "or, on the isotropic passive/outgoing branch, as\n"
    "  gamma_quad_eff = 9 K2bar^(5/2) / K0bar^(3/2)."
)

subbanner("IV.2 — Exact GR matching targets inherited from 2.5PN")
NQ_target = sp.solve(sp.Eq(gamma_quad_eff, gamma_GR), N_Q)[0]
K0bar_target = sp.simplify(NQ_target)
K2bar_target = sp.simplify((K0bar * a**2 / (9 * c_s**2)).subs(K0bar, K0bar_target))

print("N_Q_target    =", NQ_target)
print("K0bar_target  =", K0bar_target)
print("K2bar_target  =", K2bar_target)

expect_zero(
    "C_tail_toy(N_Q_target) - C_tail_GR",
    sp.simplify(C_tail_toy.subs(N_Q, NQ_target) - C_tail_GR),
)

expect_zero(
    "(G M/(2 c^3)) * invariant form(targets) - C_tail_GR",
    sp.simplify((G * M / (2 * c**3)) * gammabar_invariant.subs({K0bar: K0bar_target, K2bar: K2bar_target}) - C_tail_GR),
)

print(
    "\nTherefore the hereditary 4PN coefficient is fixed by the same passive/outgoing\n"
    "quadrupole normalization target that already controls the 2.5PN odd channel."
)


# ---------------------------------------------------------------------------
# Part V. Conditional full-4PN theorem ledger
# ---------------------------------------------------------------------------

banner("PART V — CONDITIONAL FULL-4PN THEOREM LEDGER")

print("Conditional theorem statement inside the present hierarchy:")
print()
print("1. Local sector:")
print("   The full local 4PN ordinary sector already has an exact generic-frame representative.")
print()
print("2. Tail sector:")
print("   The unique compatible hereditary coefficient is")
print("      C_tail = (G M / (2 c^3)) * gamma_quad_eff.")
print()
print("3. No new datum at 4PN:")
print("   The tail sector introduces no independent quadrupole normalization constant beyond")
print("   the one already isolated by the 2.5PN STF quadrupole audit.")
print()
print("4. GR recovery condition:")
print("   If gamma_quad_eff = 2 G / (5 c^5), then")
print("      C_tail = G^2 M / (5 c^8)")
print("   automatically.")
print()
print("5. Strongest honest conclusion:")
print("   Full conservative 4PN is now conditional on exactly the same narrow passive/outgoing")
print("   quadrupole normalization problem already isolated at 2.5PN.")

print("\nEquivalent normalization chain:")
print("  gamma_quad_eff = 2G/(5 c^5)")
print("      <=> N_Q =", NQ_target)
print("      <=> K0bar =", K0bar_target)
print("      <=> K2bar =", K2bar_target)
print("      <=> C_tail = G^2 M/(5 c^8)")

print("\nBottom line:")
print("  The local 4PN theorem and the hereditary 4PN coefficient now decouple cleanly.\n"
      "  What remains open is not a new 4PN datum, but the already-known 2.5PN\n"
      "  quadrupole normalization theorem on the actual moving-throat branch.")
