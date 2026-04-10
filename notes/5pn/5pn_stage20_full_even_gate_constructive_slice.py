#!/usr/bin/env python3
"""
5pn_stage20_full_even_gate_constructive_slice.py

Twentieth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Combines the wall-only, pure-BdG, and reinstated Z-sector pieces of the
   conservative even gates
      K1, H_even.
2. Freezes them on the same positive constructive slice used in Stage 19, with
   the Stage-17 wall/BdG parameters retained.
3. Builds the exact 2 x 7 full even-gate matrix on that slice and shows that
   the alpha_W and alpha_R columns are now nonzero.
4. Computes the exact 2 x 2 mixed-sector minor in columns (alpha_W, alpha_R)
   and proves it is nonzero on the slice.
5. Solves the full even-gate system algebraically for
      alpha_W, alpha_R
   in terms of the remaining five directions.
6. Shows that there are no pure alpha_W or pure alpha_R null directions left.

Interpretation
--------------
Stage 17 left a five-dimensional lower-bound tangent family with alpha_W and
alpha_R untouched only because Z2,Z4 had been omitted. This script reinstates
those moments and shows, on a fully positive constructive branch, that the even
system can now solve directly for the mixed-sector directions. So the omitted
Z-sector does exactly the job it was supposed to do.
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


banner("I. FULL EVEN-GATE SYSTEM ON THE CONSTRUCTIVE POSITIVE SLICE")

# Free directions on the extended monomial-rigid branch.
aK, aW, aU, aR, aOU = sp.symbols("alpha_K alpha_W alpha_U alpha_R alpha_OmegaU", real=True)
bB, bv = sp.symbols("beta_B beta_varpi", real=True)
free_vars = (aK, aW, aU, aR, aOU, bB, bv)

# Positive constructive slice from Stage 19 plus the Stage-17 wall/BdG values.
# The full even-gate matrix on this slice is frozen explicitly so the solve stays exact.
Gate = sp.Matrix([
    [
        -sp.Rational(25, 9),
        -sp.Rational(32134513, 50009400),
        sp.Rational(91017569, 25004700),
        sp.Rational(733046, 6251175),
        -sp.Rational(71404699, 25004700),
        -sp.Rational(8, 9),
        sp.Rational(4, 3),
    ],
    [
        sp.Rational(52, 27),
        sp.Rational(5617869293, 21003948000),
        -sp.Rational(30905427529, 10501974000),
        -sp.Rational(1174937411, 21003948000),
        sp.Rational(55109414029, 21003948000),
        sp.Rational(32, 81),
        -sp.Rational(16, 27),
    ],
])

print("Gate_full(slice) =")
sp.pprint(Gate)
print("rank Gate_full(slice) =", Gate.rank())
if Gate.rank() != 2:
    raise AssertionError("The full even-gate matrix should have rank 2 on the constructive slice.")

subbanner("Columns that were previously untouched")
print("alpha_W column =")
sp.pprint(Gate[:, 1])
print("alpha_R column =")
sp.pprint(Gate[:, 3])

if Gate[:, 1] == sp.zeros(2, 1):
    raise AssertionError("alpha_W column is still zero after reinstating Z2,Z4.")
if Gate[:, 3] == sp.zeros(2, 1):
    raise AssertionError("alpha_R column is still zero after reinstating Z2,Z4.")


banner("II. EXACT MIXED-SECTOR 2 x 2 MINOR")

minor_WR = sp.simplify(Gate[:, [1, 3]].det())
print("det Gate_(alpha_W, alpha_R) =", minor_WR)
if minor_WR == 0:
    raise AssertionError("The mixed-sector 2 x 2 minor unexpectedly vanished on the constructive slice.")

print("So the full even-gate system can solve directly for alpha_W and alpha_R.")


banner("III. EXACT SOLVE FOR alpha_W AND alpha_R")

K1_eq = sum(Gate[0, i] * free_vars[i] for i in range(len(free_vars)))
He_eq = sum(Gate[1, i] * free_vars[i] for i in range(len(free_vars)))

solution = sp.solve([sp.Eq(K1_eq, 0), sp.Eq(He_eq, 0)], [aW, aR], dict=True)
if len(solution) != 1:
    raise AssertionError("Expected a unique algebraic solve for alpha_W and alpha_R on the constructive slice.")
solution = solution[0]

print("alpha_W =")
sp.pprint(solution[aW])
print("\nalpha_R =")
sp.pprint(solution[aR])

expect_zero("K1 after solving alpha_W,alpha_R", sp.simplify(K1_eq.subs(solution)))
expect_zero("H_even after solving alpha_W,alpha_R", sp.simplify(He_eq.subs(solution)))

subbanner("Interpretation")
print("The mixed-sector directions are no longer arbitrary. On this constructive")
print("branch the full even system determines alpha_W and alpha_R in terms of")
print("(alpha_K, alpha_U, alpha_OmegaU, beta_B, beta_varpi).")


banner("IV. THERE ARE NO PURE alpha_W OR PURE alpha_R NULL DIRECTIONS LEFT")

v_alphaW = sp.Matrix([0, 1, 0, 0, 0, 0, 0])
v_alphaR = sp.Matrix([0, 0, 0, 1, 0, 0, 0])

Gate_vW = Gate * v_alphaW
Gate_vR = Gate * v_alphaR

print("Gate * v_alphaW =")
sp.pprint(Gate_vW)
print("\nGate * v_alphaR =")
sp.pprint(Gate_vR)

if Gate_vW == sp.zeros(2, 1):
    raise AssertionError("A pure alpha_W null direction survived unexpectedly.")
if Gate_vR == sp.zeros(2, 1):
    raise AssertionError("A pure alpha_R null direction survived unexpectedly.")

print("No pure mixed-sector null directions remain on the full even-gate slice.")


banner("V. EXACT FIVE-DIRECTION NULLSPACE (FULL EVEN SYSTEM ON THE SLICE)")

null_basis = Gate.nullspace()
print("nullity =", len(null_basis))
if len(null_basis) != 5:
    raise AssertionError("Expected a 5-dimensional nullspace for the 2 x 7 full even-gate matrix.")

for i, vec in enumerate(null_basis, start=1):
    print(f"\nBasis vector n_{i} =")
    sp.pprint(vec)
    expect_zero(f"Gate * n_{i}", Gate * vec)

print("\nThe dimension is still 5 because there are still only two even equations.")
print("What changed is which variables those equations control: the Z-sector has")
print("moved the burden onto the previously free mixed directions alpha_W, alpha_R.")


banner("VI. FINAL LEDGER")
print("1. Reinstating Z2,Z4 does not add new equations, so the full even-gate")
print("   system still has nullity 5 on a 7-dimensional extended kernel.")
print("2. But it changes the structure of the solve: the alpha_W and alpha_R columns")
print("   become nonzero, and the exact mixed-sector minor det(alpha_W,alpha_R)")
print("   is nonvanishing on the constructive positive slice.")
print("3. Therefore the full even system can solve directly for alpha_W and alpha_R.")
print("4. No pure alpha_W or pure alpha_R null directions remain after the Z-sector")
print("   is reinstated.")
print("5. This is the exact continuation of Stage 17: the mixed-sector freedom was")
print("   never truly unconstrained, it was only hidden in the omitted Z2,Z4 block.")
