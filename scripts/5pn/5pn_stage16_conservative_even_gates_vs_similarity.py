#!/usr/bin/env python3
"""
5pn_stage16_conservative_even_gates_vs_similarity.py

Sixteenth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Builds the exact Stage-10/13 weak-axisymmetric normalization defect
      Xi = A_tr Sigma_tr + Sigma_nt,
   together with the Stage-6 wall-only and pure-BdG conservative even gates
      K1, H_even.
2. Evaluates all three objects on representative directions that already lie in
   the Stage-13/14 zero-defect similarity kernel:
      - matched wall-only co-scaling,
      - BdG amplitude direction,
      - BdG frequency direction,
      - pure BdG self-similar branch x_c = x_varpi.
3. Proves that Xi vanishes on all of those directions while K1 and H_even remain
   generically nonzero.

Interpretation
--------------
This is the cleanest executable demonstration yet that the similarity-orbit /
monomial-rigidity theorem and the conservative 5PN even gates are genuinely
complementary. The orbit theorem kills the normalization defect, but it does not
close the even-preserving or hidden-even conservative slots.
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


banner("I. THE NORMALIZATION DEFECT Xi AND THE CONSERVATIVE EVEN GATES")

# Weak-axisymmetric normalization sector from Stages 10-13.
chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
Estar, Fstar, Atr = sp.symbols("E_star F_star A_tr", real=True)

dw, du, dr = sp.symbols("dln_GW dln_GU dln_R", real=True)
dK, dM = sp.symbols("dln_K dln_M", real=True)
dOU, dOW = sp.symbols("dln_Omega_U dln_Omega_W", real=True)
dDelta = sp.symbols("dln_delta_U", real=True)

dln_lamB, dln_varpi = sp.symbols("dln_lambda_B dln_varpi", real=True)

Sigma_tr = sp.expand((1 + deltaU) * (dr + du - dw - 2 * dOU) + (1 + chi0) * dDelta)
Sigma_nt = sp.expand(dM + 2 * dw - dK - 4 * dOW + 2 * Estar * (dr - dOU - dOW) - Fstar * dDelta)
Xi = sp.expand(Atr * Sigma_tr + Sigma_nt)

# Conservative even-gate model pieces.
K, M = sp.symbols("K M", positive=True, real=True)
B0, varpi = sp.symbols("B0 varpi", positive=True, real=True)

# Wall-only pieces.
K1_wall = sp.simplify(K * dK / 9 - M * dM)
He_wall = sp.simplify(-K * dK / 27 + sp.Rational(2, 3) * M * dM)

# Pure BdG pieces from Stage 6: x_c = dln lambda_B, x_varpi = dln varpi.
D01_B = sp.simplify(-2 * B0 * (dln_lamB - dln_varpi))
D21_B = sp.simplify(-2 * B0 * (dln_lamB - 2 * dln_varpi) / varpi**2)
D41_B = sp.simplify(-2 * B0 * (dln_lamB - 3 * dln_varpi) / varpi**4)

K1_B = sp.simplify(D21_B + D01_B / 9)
He_B = sp.simplify(D41_B - sp.Rational(2, 3) * D21_B - D01_B / 27)

print("Xi      =")
sp.pprint(Xi)
print("\nK1_wall =", K1_wall)
print("He_wall =", He_wall)
print("\nK1_B    =")
sp.pprint(K1_B)
print("\nHe_B    =")
sp.pprint(He_B)


banner("II. REPRESENTATIVE DIRECTIONS INSIDE THE ZERO-DEFECT SIMILARITY KERNEL")

# 1) Matched wall-only co-scaling direction from Stage 13.
wall_match = {
    dw: 0, du: 0, dr: 0, dOU: 0, dOW: 0, dDelta: 0,
    dK: 1, dM: 1,
    dln_lamB: 0, dln_varpi: 0,
}

subbanner("Matched wall-only co-scaling")
expect_zero("Xi on matched wall co-scaling", Xi.subs(wall_match))
print("K1_wall on matched wall co-scaling =", sp.simplify(K1_wall.subs(wall_match)))
print("He_wall on matched wall co-scaling =", sp.simplify(He_wall.subs(wall_match)))

# 2) Explicit support-blind directions from Stage 14.
e_B = {
    dw: 0, du: 0, dr: 0, dK: 0, dM: 0, dOU: 0, dOW: 0, dDelta: 0,
    dln_lamB: 1, dln_varpi: 0,
}
e_varpi = {
    dw: 0, du: 0, dr: 0, dK: 0, dM: 0, dOU: 0, dOW: 0, dDelta: 0,
    dln_lamB: 0, dln_varpi: 1,
}

subbanner("BdG amplitude direction e_B")
expect_zero("Xi on e_B", Xi.subs(e_B))
print("K1_B on e_B =")
sp.pprint(sp.simplify(K1_B.subs(e_B)))
print("He_B on e_B =")
sp.pprint(sp.simplify(He_B.subs(e_B)))

subbanner("BdG frequency direction e_varpi")
expect_zero("Xi on e_varpi", Xi.subs(e_varpi))
print("K1_B on e_varpi =")
sp.pprint(sp.simplify(K1_B.subs(e_varpi)))
print("He_B on e_varpi =")
sp.pprint(sp.simplify(He_B.subs(e_varpi)))

# 3) Pure BdG self-similar branch from Stage 6.
delta = sp.symbols("delta", real=True)
bdg_selfsimilar = {
    dw: 0, du: 0, dr: 0, dK: 0, dM: 0, dOU: 0, dOW: 0, dDelta: 0,
    dln_lamB: delta, dln_varpi: delta,
}

subbanner("Pure BdG self-similar branch dln lambda_B = dln varpi = delta")
expect_zero("Xi on the BdG self-similar branch", Xi.subs(bdg_selfsimilar))
print("K1_B on the BdG self-similar branch =")
sp.pprint(sp.simplify(K1_B.subs(bdg_selfsimilar)))
print("He_B on the BdG self-similar branch =")
sp.pprint(sp.simplify(He_B.subs(bdg_selfsimilar)))

print(
    "These formulas vanish only at delta = 0 unless a separate tuning is imposed,"
    " exactly as in the Stage-6 pure-BdG no-go."
)


banner("III. FINAL THEOREM LEDGER")
print("1. The Stage-10/13 similarity-orbit theorem annihilates the weak-axisymmetric")
print("   normalization defect Xi on the exact normalized kernel.")
print("2. The matched wall-only co-scaling direction already lies in that kernel, yet")
print("   it produces")
print("      K1_wall = K/9 - M,   He_wall = -K/27 + 2M/3,")
print("   which are generically nonzero.")
print("3. The explicit support-blind BdG directions dln lambda_B and dln varpi also")
print("   lie in the zero-defect normalization kernel, but they excite K1_B and He_B.")
print("4. Even the pure-BdG self-similar branch dln lambda_B = dln varpi kills Xi")
print("   while leaving the even-preserving and hidden-even gates generically open.")
print("5. Therefore the orbit theorem and the conservative 5PN even gates are")
print("   complementary, not redundant: monomial-rigidity kills Xi, but it does not")
print("   close K1 or H_even.")
