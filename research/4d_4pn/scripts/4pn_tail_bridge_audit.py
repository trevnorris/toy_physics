#!/usr/bin/env python3
"""
4pn_tail_bridge_audit.py

Tail / hereditary bridge audit for the 4PN program.

Purpose
-------
The local 4PN theorem is already closed. What remains is the separate conservative
nonlocal tail bridge. This audit freezes the exact GR tail functional, derives its
order-reduced Newtonian quadrupole kernel, and identifies the minimal toy-model
bridge datum that still remains open.

What this script does
---------------------
1. Records the exact GR conservative 4PN tail action / Hamiltonian coefficient.
2. Verifies the exact arithmetic relation between the GR 4PN tail coefficient and
   the already-frozen 2.5PN Burke--Thorne coefficient.
3. Derives the Newtonian order-reduced STF quadrupole third derivative
      I_ij^{(3)}
   and its invariant square.
4. Computes the universal local logarithmic tail shadow
      F = (2/5) (G^2 M / c^8) (I_ij^{(3)})^2,
   and shows that it occupies only the 4PN G^4/r^4 degree-2 block.
5. Encodes the minimal toy-model tail bridge as a single scalar transport factor
      Theta_tail,
   multiplying the same canonical quadrupole channel already isolated by the
   2.5PN notes.

Interpretation
--------------
- The 4PN tail is not a new representation problem.
- The source channel is exactly the same canonical STF quadrupole that survived the
  2.5PN audit.
- After the local 4PN theorem, the remaining bridge is one scalar transport question,
  not another large generic-frame existence problem.
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


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


def stf_outer(a: sp.Matrix, b: sp.Matrix | None = None) -> sp.Matrix:
    if b is None:
        b = a
    mat = (a * b.T + b * a.T) / 2
    return sp.simplify(mat - sp.trace(mat) * sp.eye(3) / 3)


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

G, M, mu, c, c_s, s = sp.symbols("G M mu c c_s s", positive=True, real=True)
r, rd, vt, v2 = sp.symbols("r rd vt v2", positive=True, real=True)
nu = sp.symbols("nu", positive=True, real=True)
gamma_toy, gamma_GR = sp.symbols("gamma_toy gamma_GR", real=True)
Theta_tail = sp.symbols("Theta_tail", real=True)
alpha_tail_GR, alpha_tail_toy = sp.symbols("alpha_tail_GR alpha_tail_toy", real=True)


# ---------------------------------------------------------------------------
# Part I. Exact GR tail functional coefficient and its 2.5PN relation
# ---------------------------------------------------------------------------

banner("PART I — EXACT GR 4PN TAIL COEFFICIENT")

gamma_GR = sp.Rational(2, 5) * G / c**5
alpha_tail_GR = G**2 * M / (5 * c**8)

print("Canonical 2.5PN Burke–Thorne coefficient:")
print("  gamma_GR      =", gamma_GR)
print("Exact 4PN conservative tail coefficient:")
print("  alpha_tail_GR =", alpha_tail_GR)

subbanner("I.1 — Exact arithmetic bridge between 2.5PN and 4PN coefficients")
expect_zero(
    "alpha_tail_GR - (G M / (2 c^3)) gamma_GR",
    alpha_tail_GR - (G * M / (2 * c**3)) * gamma_GR,
)

print(
    "\nThis is the exact GR coefficient relation:\n"
    "  alpha_tail_GR = (G M / (2 c^3)) * gamma_GR.\n"
    "So once the canonical quadrupole coefficient gamma is fixed, the 4PN tail\n"
    "coefficient is determined up to the monopole-scattering transport factor."
)


# ---------------------------------------------------------------------------
# Part II. Newtonian order-reduced quadrupole kernel
# ---------------------------------------------------------------------------

banner("PART II — NEWTONIAN ORDER-REDUCED STF QUADRUPOLE KERNEL")

# Choose a representative orbital frame with x along the x-axis and tangential
# velocity along y. This is enough because the final contraction is rotationally
# invariant.
x = sp.Matrix([r, 0, 0])
v = sp.Matrix([rd, vt, 0])

# Newtonian acceleration and jerk.
a = -G * M * x / r**3
j = -G * M * (v / r**3 - 3 * rd * x / r**4)

# I_ij = mu STF(x_i x_j)
# I_ij^(1) = 2 mu STF(x_i v_j)
# I_ij^(2) = 2 mu STF(v_i v_j + x_i a_j)
# I_ij^(3) = 2 mu STF(3 a_i v_j + x_i j_j)
I3_direct = sp.simplify(2 * mu * (3 * stf_outer(a, v) + stf_outer(x, j)))

A = stf_outer(x, v)
B = stf_outer(x, x)
I3_formula = sp.simplify(-2 * G * M * mu / r**3 * (4 * A - 3 * rd / r * B))

subbanner("II.1 — Order-reduced STF formula")
for i in range(3):
    for j_idx in range(3):
        expect_zero(f"I3[{i},{j_idx}] direct - formula", I3_direct[i, j_idx] - I3_formula[i, j_idx])

print("\nCompact order-reduced formula:")
print("  I_ij^(3) = -2 G M mu / r^3 * ( 4 x_<i v_j> - 3 (rd/r) x_<i x_j> )")

subbanner("II.2 — Invariant square of I_ij^(3)")
I3_sq = sp.factor(sp.expand(sum(I3_formula[i, j_idx] ** 2 for i in range(3) for j_idx in range(3))))
print("I3_sq =")
sp.pprint(I3_sq)

I3_sq_expected = sp.Rational(8, 3) * G**2 * M**2 * mu**2 / r**4 * (12 * v2 - 11 * rd**2)
expect_zero("I3_sq - expected", I3_sq.subs(vt**2, v2 - rd**2) - I3_sq_expected)

print(
    "\nUniversal invariant form:\n"
    "  (I_ij^(3))^2 = (8/3) G^2 M^2 mu^2 / r^4 * (12 v^2 - 11 rdot^2)."
)


# ---------------------------------------------------------------------------
# Part III. Universal local logarithmic tail shadow
# ---------------------------------------------------------------------------

banner("PART III — UNIVERSAL LOCAL LOGARITHMIC TAIL SHADOW")

F_tail = sp.simplify(sp.Rational(2, 5) * G**2 * M / c**8 * I3_sq_expected)
print("F_tail =")
sp.pprint(F_tail)

F_tail_per_mu = sp.simplify((F_tail / mu).subs(mu, nu * M).subs(G * M / r, sp.Symbol("U")))
print("\nPer reduced mass (mu = nu M, U = G M / r):")
sp.pprint(F_tail_per_mu)

print(
    "\nBlock placement:\n"
    "  F_tail / mu = (16/15) * nu * U^4 / c^8 * (12 v^2 - 11 rdot^2).\n"
    "So the local logarithmic shadow of the 4PN tail occupies only the G^4/r^4\n"
    "degree-2 block (the U-block in the local 4PN bookkeeping)."
)

subbanner("III.1 — Reduced COM Hamiltonian form")
# In the reduced COM Hamiltonian variables one may set p^2 ~ v^2 and p_r^2 ~ rdot^2
# at Newtonian order, which is sufficient for the tail shadow.
p2, pr2 = sp.symbols("p2 pr2", real=True)
F_tail_red = sp.simplify(sp.Rational(16, 15) * nu * sp.Symbol("u")**4 * (12 * p2 - 11 * pr2))
print("F_tail_reduced =")
sp.pprint(F_tail_red)


# ---------------------------------------------------------------------------
# Part IV. Minimal toy-model tail bridge ansatz
# ---------------------------------------------------------------------------

banner("PART IV — MINIMAL TOY-MODEL TAIL BRIDGE ANSATZ")

alpha_tail_toy = sp.simplify(Theta_tail * (G * M / (2 * c_s**3)) * gamma_toy)
print("alpha_tail_toy =")
sp.pprint(alpha_tail_toy)

print(
    "\nMinimal transport parameterization:\n"
    "  alpha_tail_toy = Theta_tail * (G M / (2 c_s^3)) * gamma_toy.\n"
    "Here gamma_toy is the canonically normalized STF quadrupole coefficient already\n"
    "isolated by the 2.5PN audit, and Theta_tail is the single extra mass-scattering\n"
    "transport factor not yet derived from the moving-throat PDE."
)

subbanner("IV.1 — GR recovery check on the c_s = c branch")
expect_zero(
    "alpha_tail_toy(c_s=c,Theta=1,gamma=gamma_GR) - alpha_tail_GR",
    alpha_tail_toy.subs({c_s: c, Theta_tail: 1, gamma_toy: gamma_GR}) - alpha_tail_GR,
)

print(
    "\nSo the tail bridge is maximally narrow:\n"
    "  - the source representation is already fixed by the 2.5PN quadrupole branch,\n"
    "  - the exact GR coefficient is recovered when c_s = c, gamma_toy = gamma_GR,\n"
    "    and Theta_tail = 1,\n"
    "  - therefore the remaining 4PN hereditary gap is one scalar transport theorem,\n"
    "    not another large comparable-mass existence problem."
)


# ---------------------------------------------------------------------------
# Part V. Final theorem ledger
# ---------------------------------------------------------------------------

banner("FINAL 4PN TAIL-BRIDGE LEDGER")
print("1. Exact GR conservative 4PN tail Hamiltonian coefficient:")
print("     alpha_tail_GR = G^2 M / (5 c^8)")
print()
print("2. Exact GR bridge to the 2.5PN Burke–Thorne coefficient:")
print("     alpha_tail_GR = (G M / (2 c^3)) * gamma_GR")
print("     with gamma_GR = 2 G / (5 c^5)")
print()
print("3. Universal Newtonian quadrupole shadow:")
print("     (I_ij^(3))^2 = (8/3) G^2 M^2 mu^2 / r^4 * (12 v^2 - 11 rdot^2)")
print()
print("4. Therefore the local logarithmic tail shadow is")
print("     F_tail / mu = (16/15) * nu * U^4 / c^8 * (12 v^2 - 11 rdot^2)")
print("   so it lives only in the 4PN G^4/r^4 degree-2 block.")
print()
print("5. Minimal toy-model tail bridge ansatz:")
print("     alpha_tail_toy = Theta_tail * (G M / (2 c_s^3)) * gamma_toy")
print()
print("6. Hence, after the local 4PN theorem, the remaining hereditary bridge is one")
print("   scalar transport question: does the moving-throat model give Theta_tail = 1")
print("   on the same canonical STF quadrupole branch already isolated by the 2.5PN work?")
