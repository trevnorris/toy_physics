#!/usr/bin/env python3
"""
5pn_stage17_lower_bound_even_gate_solution.py

Seventeenth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Starts from the exact Stage-13 normalized zero-defect kernel with the two
   Stage-14 support-blind BdG directions added.
2. Adds the minimal conservative lower-bound even gates obtained by keeping the
   explicit wall-only and pure-BdG pieces of
      K1, H_even.
3. Builds the exact 2 x 7 gate matrix acting on the seven free parameters of the
   extended monomial-rigid branch and proves it has rank 2 generically.
4. Solves the lower-bound gates algebraically for
      dln K and dln lambda_B
   in terms of the remaining five free parameters.
5. Constructs one exact five-direction null basis for the lower-bound system.

Interpretation
--------------
This is not the full 5PN theorem: the conservative Maxwell/mixed Z2/Z4 pieces
are still omitted here on purpose. But it is the first exact intersection
calculation showing how the monomial-rigid kernel is reduced further once the
explicit even gates are imposed. It gives a clean lower-bound picture of the
remaining 5PN tangent freedom before the Z-sector is reinstated.
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


banner("I. EXTENDED MONOMIAL-RIGID KERNEL WITH EXPLICIT LOWER-BOUND EVEN GATES")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)
K, M = sp.symbols("K M", positive=True, real=True)
B0, varpi = sp.symbols("B0 varpi", positive=True, real=True)

# Seven free parameters on the extended monomial-rigid kernel:
#   five normalized similarity directions + two support-blind BdG directions.
aK, aW, aU, aR, aOU = sp.symbols("alpha_K alpha_W alpha_U alpha_R alpha_OmegaU", real=True)
bB, bv = sp.symbols("beta_B beta_varpi", real=True)

# Stage-13 kernel formulas.
dw = aW
du = aU
dr = aR
dK = aK
dOU = aOU
dDelta = sp.simplify(-((1 + deltaU) / (1 + chi0)) * (dr + du - dw - 2 * dOU))
dM = sp.simplify(dK - 2 * du + 2 * dOU)
dOW = sp.simplify((dw - du + (1 - Estar) * dOU + Estar * dr - sp.Rational(1, 2) * Fstar * dDelta) / (Estar + 2))

dln_lamB = bB
dln_varpi = bv

print("The extended monomial-rigid kernel is parameterized by")
print("  (alpha_K, alpha_W, alpha_U, alpha_R, alpha_OmegaU, beta_B, beta_varpi).")

# Lower-bound conservative even gates: wall-only + pure-BdG pieces.
K1_wall = sp.simplify(K * dK / 9 - M * dM)
He_wall = sp.simplify(-K * dK / 27 + sp.Rational(2, 3) * M * dM)

D01_B = sp.simplify(-2 * B0 * (dln_lamB - dln_varpi))
D21_B = sp.simplify(-2 * B0 * (dln_lamB - 2 * dln_varpi) / varpi**2)
D41_B = sp.simplify(-2 * B0 * (dln_lamB - 3 * dln_varpi) / varpi**4)

K1_B = sp.simplify(D21_B + D01_B / 9)
He_B = sp.simplify(D41_B - sp.Rational(2, 3) * D21_B - D01_B / 27)

K1_min = sp.expand(K1_wall + K1_B)
He_min = sp.expand(He_wall + He_B)

print("K1_min =")
sp.pprint(K1_min)
print("\nHe_min =")
sp.pprint(He_min)


banner("II. EXACT 2 x 7 LOWER-BOUND GATE MATRIX")

free_vars = (aK, aW, aU, aR, aOU, bB, bv)
Gate = sp.Matrix([
    [sp.diff(K1_min, v) for v in free_vars],
    [sp.diff(He_min, v) for v in free_vars],
])

print("Gate =")
sp.pprint(Gate)
print("rank Gate =", Gate.rank())

subbanner("Immediate structural reading")
print("The columns for alpha_W and alpha_R vanish in the lower-bound wall+BdG model.")
print("So those two normalized similarity directions are untouched until the omitted")
print("Maxwell/mixed Z2,Z4 sector is reinstated.")
expect_zero("alpha_W column vanishes", Gate[:, 1])
expect_zero("alpha_R column vanishes", Gate[:, 3])


banner("III. EXACT LOWER-BOUND SOLVE FOR alpha_K AND beta_B")

solution = sp.solve([sp.Eq(K1_min, 0), sp.Eq(He_min, 0)], [aK, bB], dict=True)
if len(solution) != 1:
    raise AssertionError("Expected a unique algebraic lower-bound solve for alpha_K and beta_B.")
solution = solution[0]

alphaK_sol = sp.simplify(sp.factor(solution[aK]))
betaB_sol = sp.simplify(sp.factor(solution[bB]))

print("alpha_K =")
sp.pprint(alphaK_sol)
print("\nbeta_B =")
sp.pprint(betaB_sol)

expect_zero("K1_min after solve", K1_min.subs(solution))
expect_zero("He_min after solve", He_min.subs(solution))

print(
    "So, at this lower-bound level, imposing the even gates determines the wall"
    " amplitude scaling alpha_K and the BdG amplitude drift beta_B in terms of"
    " the remaining five free directions."
)


banner("IV. EXACT FIVE-DIRECTION NULL BASIS OF THE LOWER-BOUND SYSTEM")

null_basis = Gate.nullspace()
print("nullity =", len(null_basis))
if len(null_basis) != 5:
    raise AssertionError("Expected a 5-dimensional nullspace for the lower-bound 2 x 7 gate matrix.")

for i, vec in enumerate(null_basis, start=1):
    print(f"\nBasis vector v_{i} =")
    sp.pprint(vec)
    expect_zero(f"Gate * v_{i}", Gate * vec)

subbanner("Interpretation")
print("v_1: pure alpha_W direction (survives because the lower-bound wall+BdG gates")
print("     do not yet see the omitted Z2,Z4 mixed-sector moments).")
print("v_2: pure alpha_R direction (same caveat).")
print("v_3: matched alpha_U/alpha_OmegaU direction.")
print("v_4: BdG-amplitude deformation beta_B with exact compensating alpha_K, alpha_U.")
print("v_5: BdG-frequency deformation beta_varpi with exact compensating alpha_K, alpha_U.")


banner("V. FINAL THEOREM LEDGER")
print("1. Start with the 7-dimensional extended monomial-rigid kernel:")
print("      5 normalized similarity directions + 2 support-blind BdG directions.")
print("2. Add the explicit wall-only and pure-BdG lower-bound even gates K1 = H_even = 0.")
print("3. The resulting 2 x 7 gate matrix has generic rank 2, so the lower-bound")
print("   even-gate intersection has dimension 5.")
print("4. A convenient algebraic solve is to determine")
print("      alpha_K and beta_B")
print("   in terms of the remaining five directions.")
print("5. Two directions, alpha_W and alpha_R, remain completely untouched at this")
print("   lower-bound stage. Their eventual fate is deferred to the omitted")
print("   Maxwell/mixed Z2,Z4 sector.")
print("6. So the next honest 5PN task is now sharply localized:")
print("   reinstate the Z2,Z4 conservative mixed-sector moments and see how they cut")
print("   the remaining 5-dimensional lower-bound tangent family.")
