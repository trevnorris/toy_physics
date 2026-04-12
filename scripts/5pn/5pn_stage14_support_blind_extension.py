#!/usr/bin/env python3
"""
5pn_stage14_support_blind_extension.py

Fourteenth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Extends the Stage-13 normalized monomial-drift matrix by the explicit Stage-5
   BdG-support drifts
      dln lambda_B, dln varpi.
2. Proves that the Stage-11 weak-axisymmetric monomials are exactly support-blind:
   the BdG columns of the extended matrix vanish identically.
3. Shows that the extended primitive drift space therefore splits as
      ker(M_ext) = span{BdG amplitude, BdG frequency} (+) ker(M_norm),
   with dimension 7.
4. States the practical consequence: monomial-rigidity / similarity-orbit tests do
   not by themselves constrain the BdG primitive drifts. Those directions must be
   fixed, if at all, by the separate conservative 5PN front-end conditions
   (K1, hidden-even consistency, single-pole test, etc.).

Interpretation
--------------
This is the missing conceptual bridge between the Stage-5 primitive overlap model
and the Stage-10/11 similarity-orbit theorem. The latter governs the mixed/
wall/U normalization problem, but it deliberately does not close the BdG support
sector. That is why the conservative grouped-P2 extraction still matters.
"""

from __future__ import annotations

import sympy as sp



def banner(title: str) -> None:
    line = "=" * 88
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


banner("I. EXTENDED PRIMITIVE DRIFT MATRIX INCLUDING THE BDG SUPPORT DRIFTS")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)

# Full Stage-5-style primitive log drifts.
dln_lamB, dln_varpi = sp.symbols("dln_lambda_B dln_varpi", real=True)
dw, du, dr = sp.symbols("dln_GW dln_GU dln_R", real=True)
dK, dM = sp.symbols("dln_K dln_M", real=True)
dOU, dOW = sp.symbols("dln_Omega_U dln_Omega_W", real=True)
dDelta = sp.symbols("dln_delta_U", real=True)

# The direct Stage-11 monomials do not depend on the BdG support variables.
dlog_Ctr = sp.expand((1 + deltaU) * (dr + du - dw - 2 * dOU) + (1 + chi0) * dDelta)
dlog_Cnt = sp.expand(dM + 2 * dw - dK - 4 * dOW + 2 * Estar * (dr - dOU - dOW) - Fstar * dDelta)
dlog_eta = sp.expand(dM + 2 * du - dK - 2 * dOU)

vars_ext = (dln_lamB, dln_varpi, dw, du, dr, dK, dM, dOU, dOW, dDelta)
Mext = sp.Matrix([
    [sp.expand(dlog_Ctr).coeff(v) for v in vars_ext],
    [sp.expand(dlog_Cnt).coeff(v) for v in vars_ext],
    [sp.expand(dlog_eta).coeff(v) for v in vars_ext],
])

print("M_ext =")
sp.pprint(Mext)

expect_zero("BdG amplitude column vanishes", Mext[:, 0])
expect_zero("BdG frequency column vanishes", Mext[:, 1])
print("rank M_ext =", Mext.rank())
print("dim ker(M_ext) =", Mext.shape[1] - Mext.rank())


banner("II. DIRECT-SUM DECOMPOSITION OF THE EXTENDED ZERO-DEFECT TANGENT SPACE")

# The normalized Stage-13 kernel matrix is the right 3x8 block.
Mnorm = Mext[:, 2:]
print("M_norm =")
sp.pprint(Mnorm)
print("rank M_norm =", Mnorm.rank())
print("dim ker(M_norm) =", Mnorm.shape[1] - Mnorm.rank())

# Two exact support-blind basis directions.
e_B = sp.Matrix([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
e_varpi = sp.Matrix([0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
expect_zero("M_ext e_B", Mext * e_B)
expect_zero("M_ext e_varpi", Mext * e_varpi)

print("So dln lambda_B and dln varpi are exact support-blind weak-axisymmetric directions.")
print("They do not enter the direct monomial drifts at all.")

# Build the obvious injection of the Stage-13 normalized kernel into the extended space.
alphaK, alphaW, alphaU, alphaR, alphaOU = sp.symbols("alpha_K alpha_W alpha_U alpha_R alpha_OmegaU", real=True)

# Stage-13 normalized compatibility formulas.
dDelta_expr = sp.simplify(-((1 + deltaU) / (1 + chi0)) * (dr + du - dw - 2 * dOU))
dM_expr = sp.simplify(dK - 2 * du + 2 * dOU)
dOW_expr = sp.simplify((dw - du + (1 - Estar) * dOU + Estar * dr - sp.Rational(1, 2) * Fstar * dDelta_expr) / (Estar + 2))

kernel_norm = sp.Matrix([
    dw,
    du,
    dr,
    dK,
    dM_expr,
    dOU,
    dOW_expr,
    dDelta_expr,
])

kernel_ext = sp.Matrix([
    0,
    0,
    *list(kernel_norm),
])
expect_zero("M_ext * injected normalized kernel", sp.simplify(Mext * kernel_ext))

print("Therefore")
print("  ker(M_ext) = span{e_B, e_varpi} (+) ker(M_norm),")
print("with dimension 2 + 5 = 7.")


banner("III. PRACTICAL CONSEQUENCE FOR THE 5PN PROGRAM")

print("The Stage-10/11 similarity-orbit / monomial-rigidity theorem constrains only")
print("the mixed/wall/U normalization sector:")
print("  (G_W, G_U, R, K, M, Omega_U, Omega_W, delta_U).")
print("It does NOT constrain the explicit BdG support drifts")
print("  dln lambda_B, dln varpi.")
print()
print("So monomial-rigidity alone cannot finish the full 5PN conservative problem.")
print("The BdG support directions remain to be fixed, if at all, by the separate")
print("conservative grouped-P2 front-end conditions such as")
print("  - K1 = 0,")
print("  - the hidden-even consistency slot, and")
print("  - the O(omega^4) single-pole / grouped-response test.")
print()
print("This is why the Stage-5/6 conservative extraction is still essential even")
print("after the Stage-10/11 similarity-orbit collapse.")


banner("IV. FINAL THEOREM LEDGER")
print("1. Extending the primitive Stage-5 drift space by the BdG drifts")
print("      dln lambda_B, dln varpi")
print("   adds two exact zero columns to the Stage-11 monomial-drift matrix.")
print("2. The direct monomial drifts are therefore exactly support-blind.")
print("3. The extended primitive zero-defect tangent space has dimension 7:")
print("      2 support-blind BdG directions + 5 normalized similarity directions.")
print("4. Consequently, the weak-axisymmetric similarity-orbit theorem and the")
print("   conservative grouped-P2 extraction theorem are complementary, not redundant.")
print("5. The normalization problem lives in the mixed/wall/U sector; the remaining")
print("   conservative 5PN closure still has to control the BdG support drifts separately.")
