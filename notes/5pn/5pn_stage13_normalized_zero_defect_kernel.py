#!/usr/bin/env python3
"""
5pn_stage13_normalized_zero_defect_kernel.py

Thirteenth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Builds the exact Stage-168/169 monomial-drift matrix directly in the
   Stage-5-style normalized variables
      (G_W, G_U, R, K, M, Omega_U, Omega_W, delta_U).
2. Proves that the normalized matrix has exact rank 3, so the normalized
   zero-defect branch has dimension 5.
3. Derives the exact normalized compatibility ledger in a triangular form:
      - tracking fixes d ln delta_U,
      - dressing fixes d ln M,
      - nontracking fixes d ln Omega_W.
4. Builds one exact five-direction normalized kernel basis and proves it spans
   ker(M_norm).
5. Restricts to the original Stage-5 unsplit slice d ln delta_U = 0 and shows
   that the kernel drops to a 4-dimensional slice obeying an extra tracking
   constraint.
6. Translates the normalized log-drift ledger back into practical absolute-slope
   formulas for the Stage-5 primitive variables
      dK, dM, d(lambda_U), d(lambda_W), d(lambda_R), d(Omega_U), d(Omega_W), d(delta_U).

Interpretation
--------------
After this stage the Stage-10/11 similarity-orbit theorem is usable directly in
Stage-5 language. The normalized primitive prototype, once augmented by the
split-U variable delta_U, already has the exact codimension-3 zero-defect law.
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


# ---------------------------------------------------------------------------
# I. Exact normalized monomial-drift matrix
# ---------------------------------------------------------------------------

banner("I. EXACT NORMALIZED MONOMIAL-DRIFT MATRIX")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)

# Normalized log drifts.
dw, du, dr = sp.symbols("dln_GW dln_GU dln_R", real=True)
dK, dM = sp.symbols("dln_K dln_M", real=True)
dOU, dOW = sp.symbols("dln_Omega_U dln_Omega_W", real=True)
dDelta = sp.symbols("dln_delta_U", real=True)

# Direct monomial drifts in normalized variables.
dlog_Ctr = sp.expand((1 + deltaU) * (dr + du - dw - 2 * dOU) + (1 + chi0) * dDelta)
dlog_Cnt = sp.expand(dM + 2 * dw - dK - 4 * dOW + 2 * Estar * (dr - dOU - dOW) - Fstar * dDelta)
dlog_eta = sp.expand(dM + 2 * du - dK - 2 * dOU)

vars_vec = (dw, du, dr, dK, dM, dOU, dOW, dDelta)
Mnorm = sp.Matrix([
    [sp.expand(dlog_Ctr).coeff(v) for v in vars_vec],
    [sp.expand(dlog_Cnt).coeff(v) for v in vars_vec],
    [sp.expand(dlog_eta).coeff(v) for v in vars_vec],
])

print("M_norm =")
sp.pprint(Mnorm)

Mnorm_expected = sp.Matrix([
    [-(1 + deltaU), 1 + deltaU, 1 + deltaU, 0, 0, -2 * (1 + deltaU), 0, 1 + chi0],
    [2, 0, 2 * Estar, -1, 1, -2 * Estar, -(4 + 2 * Estar), -Fstar],
    [0, 2, 0, -1, 1, -2, 0, 0],
])
expect_zero("M_norm - expected", Mnorm - Mnorm_expected)

subbanner("Rank test")
minor = Mnorm[:, [7, 4, 6]]  # columns (dln delta_U, dln M, dln Omega_W)
minor_det = sp.simplify(minor.det())
print("det M_norm^(delta_U,M,Omega_W) =", minor_det)
expect_zero("minor determinant - 2 (E_* + 2) (1 + chi_0)", minor_det - 2 * (Estar + 2) * (1 + chi0))
print("Therefore rank M_norm = 3 and dim ker(M_norm) = 5 on the constructive branch E_* != -2.")


# ---------------------------------------------------------------------------
# II. Exact triangular compatibility solve
# ---------------------------------------------------------------------------

banner("II. EXACT TRIANGULAR COMPATIBILITY SOLVE")

# Tracking and dressing solve immediately.
dDelta_expr = sp.simplify(-((1 + deltaU) / (1 + chi0)) * (dr + du - dw - 2 * dOU))
dM_expr = sp.simplify(dK - 2 * du + 2 * dOU)

print("d ln delta_U =", dDelta_expr)
print("d ln M       =", dM_expr)

# Nontracking solve after substituting the first two equations.
eq_nt = sp.simplify(dlog_Cnt.subs({dDelta: dDelta_expr, dM: dM_expr}))
dOW_expr = sp.simplify(sp.solve(sp.Eq(eq_nt, 0), dOW)[0])

print("d ln Omega_W =", dOW_expr)
expect_zero("tracking equation after substitution", dlog_Ctr.subs(dDelta, dDelta_expr))
expect_zero("dressing equation after substitution", dlog_eta.subs(dM, dM_expr))
expect_zero("nontracking equation after substitution", dlog_Cnt.subs({dDelta: dDelta_expr, dM: dM_expr, dOW: dOW_expr}))

subbanner("More readable nontracking form")
dOW_compact = sp.simplify((dw - du + (1 - Estar) * dOU + Estar * dr - sp.Rational(1, 2) * Fstar * dDelta_expr) / (Estar + 2))
print("compact d ln Omega_W =", dOW_compact)
expect_zero("compact Omega_W formula", dOW_expr - dOW_compact)


# ---------------------------------------------------------------------------
# III. Exact five-direction normalized kernel basis
# ---------------------------------------------------------------------------

banner("III. EXACT FIVE-DIRECTION NORMALIZED KERNEL BASIS")

alphaK, alphaW, alphaU, alphaR, alphaOU = sp.symbols("alpha_K alpha_W alpha_U alpha_R alpha_OmegaU", real=True)

free_assign = {
    dK: alphaK,
    dw: alphaW,
    du: alphaU,
    dr: alphaR,
    dOU: alphaOU,
}

kernel_vector = sp.Matrix([
    free_assign[dw],
    free_assign[du],
    free_assign[dr],
    free_assign[dK],
    sp.simplify(dM_expr.subs(free_assign)),
    free_assign[dOU],
    sp.simplify(dOW_expr.subs(free_assign)),
    sp.simplify(dDelta_expr.subs(free_assign)),
])

print("General kernel vector in free parameters (alpha_K, alpha_W, alpha_U, alpha_R, alpha_OmegaU) =")
sp.pprint(kernel_vector)
expect_zero("M_norm * kernel_vector", sp.simplify(Mnorm * kernel_vector))

# Extract an explicit basis by switching on one free parameter at a time.
basis_vectors = []
for i, a in enumerate((alphaK, alphaW, alphaU, alphaR, alphaOU)):
    subs = {alphaK: 0, alphaW: 0, alphaU: 0, alphaR: 0, alphaOU: 0}
    subs[a] = 1
    vec = sp.simplify(kernel_vector.subs(subs))
    basis_vectors.append(vec)
    print(f"\nBasis vector b_{i+1} =")
    sp.pprint(vec)

Onorm = sp.Matrix.hstack(*basis_vectors)
print("\nO_norm =")
sp.pprint(Onorm)
expect_zero("M_norm O_norm", sp.simplify(Mnorm * Onorm))
print("rank O_norm =", Onorm.rank())

subbanner("Interpretation of the first basis direction")
print("b_1 has only dln K = dln M = 1 and all other drifts zero.")
print("So the matched wall-only co-scaling (K,M) -> e^s (K,M) is an exact similarity direction.")


# ---------------------------------------------------------------------------
# IV. Restriction to the original Stage-5 unsplit slice
# ---------------------------------------------------------------------------

banner("IV. RESTRICTION TO THE UNSPLIT STAGE-5 SLICE dln delta_U = 0")

M_unsplit = Mnorm[:, :-1]  # remove dln delta_U column
print("M_unsplit =")
sp.pprint(M_unsplit)
print("rank M_unsplit =", M_unsplit.rank())
print("dim ker(M_unsplit) =", M_unsplit.shape[1] - M_unsplit.rank())

tracking_unsplit = sp.simplify(sp.solve(sp.Eq(dlog_Ctr.subs(dDelta, 0), 0), dr)[0])
print("Tracking condition on the unsplit slice:")
print("  dln R =", tracking_unsplit)

# Solve the remaining equations on the unsplit slice.
dM_unsplit = sp.simplify(dM_expr)
dOW_unsplit = sp.simplify(dOW_expr.subs(dDelta, 0))
print("d ln M on the unsplit slice =", dM_unsplit)
print("d ln Omega_W on the unsplit slice =", dOW_unsplit)
expect_zero("unsplit kernel conditions", M_unsplit * sp.Matrix([dw, du, dr, dK, dM_unsplit, dOU, dOW_unsplit]).subs(dr, tracking_unsplit))

print("So the original Stage-5 prototype is the dln delta_U = 0 slice of the full")
print("normalized kernel. It still carries a 4-dimensional zero-defect tangent space,")
print("but only after imposing the extra tracking condition dln R + dln G_U - dln G_W - 2 dln Omega_U = 0.")


# ---------------------------------------------------------------------------
# V. Absolute-slope formulas in Stage-5 variables
# ---------------------------------------------------------------------------

banner("V. ABSOLUTE-SLOPE FORMULAS IN STAGE-5 VARIABLES")

lamU, lamW, lamR = sp.symbols("lambda_U lambda_W lambda_R", positive=True, real=True)
Kabs, Mabs = sp.symbols("K M", positive=True, real=True)
OmegaU_abs, OmegaW_abs = sp.symbols("Omega_U Omega_W", positive=True, real=True)
deltaU_abs = sp.symbols("delta_U", positive=True, real=True)

dlamU, dlamW, dlamR = sp.symbols("d_lambda_U d_lambda_W d_lambda_R", real=True)
dK_abs, dM_abs = sp.symbols("dK dM", real=True)
dOmegaU_abs, dOmegaW_abs = sp.symbols("d_Omega_U d_Omega_W", real=True)
ddeltaU_abs = sp.symbols("d_delta_U", real=True)

abs_subs = {
    dw: dlamW / lamW,
    du: dlamU / lamU,
    dr: dlamR / lamR,
    dK: dK_abs / Kabs,
    dM: dM_abs / Mabs,
    dOU: dOmegaU_abs / OmegaU_abs,
    dOW: dOmegaW_abs / OmegaW_abs,
    dDelta: ddeltaU_abs / deltaU_abs,
}

tracking_abs = sp.simplify((dlog_Ctr).subs(abs_subs))
dressing_abs = sp.simplify((dlog_eta).subs(abs_subs))
nontracking_abs = sp.simplify((dlog_Cnt).subs(abs_subs))

print("Tracking equation =", tracking_abs)
print("Dressing equation =", dressing_abs)
print("Nontracking equation =", nontracking_abs)

# Solve the triangular system in absolute slopes.
ddeltaU_expr = sp.simplify(ddeltaU_abs * dDelta_expr.subs({dw: dlamW / lamW, du: dlamU / lamU, dr: dlamR / lamR, dOU: dOmegaU_abs / OmegaU_abs, dDelta: 1}) )
# Simpler and exact: ddeltaU = delta_U * dln delta_U.
ddeltaU_expr = sp.simplify(deltaU_abs * dDelta_expr.subs({dw: dlamW / lamW, du: dlamU / lamU, dr: dlamR / lamR, dOU: dOmegaU_abs / OmegaU_abs}))
dM_expr_abs = sp.simplify(Mabs * dM_expr.subs({dK: dK_abs / Kabs, du: dlamU / lamU, dOU: dOmegaU_abs / OmegaU_abs}))
dOmegaW_expr_abs = sp.simplify(
    OmegaW_abs * dOW_expr.subs({
        dw: dlamW / lamW,
        du: dlamU / lamU,
        dr: dlamR / lamR,
        dK: dK_abs / Kabs,
        dOU: dOmegaU_abs / OmegaU_abs,
    })
)

print("d(delta_U) =", ddeltaU_expr)
print("dM         =", dM_expr_abs)
print("d(Omega_W) =", dOmegaW_expr_abs)

subbanner("Implicit triangular form")
dOmegaW_compact_abs = sp.simplify(
    OmegaW_abs
    * (
        dlamW / lamW
        - dlamU / lamU
        + (1 - Estar) * dOmegaU_abs / OmegaU_abs
        + Estar * dlamR / lamR
        - sp.Rational(1, 2) * Fstar * ddeltaU_abs / deltaU_abs
    ) / (Estar + 2)
)
print("d(Omega_W) compact =", dOmegaW_compact_abs)


# ---------------------------------------------------------------------------
# VI. Final ledger
# ---------------------------------------------------------------------------

banner("VI. FINAL THEOREM LEDGER")
print("1. The exact zero-defect equations in normalized Stage-5 variables are")
print("      (1+delta_U)(dln R + dln G_U - dln G_W - 2 dln Omega_U) + (1+chi_0) dln delta_U = 0,")
print("      dln M + 2 dln G_U - dln K - 2 dln Omega_U = 0,")
print("      dln M + 2 dln G_W - dln K - 4 dln Omega_W")
print("        + 2 E_* (dln R - dln Omega_U - dln Omega_W) - F_* dln delta_U = 0.")
print("2. The normalized monomial-drift matrix has exact rank 3, so the normalized")
print("   zero-defect branch has dimension 5.")
print("3. A triangular solve is available:")
print("      tracking fixes dln delta_U,")
print("      dressing fixes dln M,")
print("      nontracking fixes dln Omega_W.")
print("4. The matched wall-only co-scaling dln K = dln M is one exact similarity direction.")
print("5. The original Stage-5 primitive model is the dln delta_U = 0 slice of the full")
print("   normalized kernel, leaving a 4-dimensional zero-defect tangent space once the")
print("   extra tracking condition dln R + dln G_U - dln G_W - 2 dln Omega_U = 0 is imposed.")
print("6. In absolute Stage-5 slopes the exact triangular compatibility formulas are")
print("      d(delta_U) = -delta_U (1+delta_U)/(1+chi_0) [ d(lambda_R)/lambda_R + d(lambda_U)/lambda_U")
print("                     - d(lambda_W)/lambda_W - 2 d(Omega_U)/Omega_U ],")
print("      dM = M [ dK/K - 2 d(lambda_U)/lambda_U + 2 d(Omega_U)/Omega_U ],")
print("      d(Omega_W) = Omega_W/(E_*+2) [ d(lambda_W)/lambda_W - d(lambda_U)/lambda_U")
print("                     + (1-E_*) d(Omega_U)/Omega_U + E_* d(lambda_R)/lambda_R")
print("                     - (F_*/2) d(delta_U)/delta_U ].")
print("7. So the Stage-10/11 similarity-orbit theorem is now directly usable in the")
print("   Stage-5 primitive deformation language.")
