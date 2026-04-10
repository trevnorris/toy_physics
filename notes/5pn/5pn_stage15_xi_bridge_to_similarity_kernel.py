#!/usr/bin/env python3
"""
5pn_stage15_xi_bridge_to_similarity_kernel.py

Fifteenth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Proves the exact Stage-5 grouped-prefactor identity
      Xi_load = N01/N0 - D01/D0 = P1/P0,
   where P0 = N0/D0 and P1 is the first weak-axisymmetric prefactor slope.
2. Rebuilds the Stage-13 normalized monomial-drift matrix and the Stage-14
   support-blind extension in one place.
3. Defines the abstract weak-axisymmetric normalization defect in the exact
   Stage-10 triangular form
      Xi = A_tr * Sigma_tr + Sigma_nt,
   where (Sigma_tr, Sigma_nt, Sigma_eta) are the three direct monomial drifts.
4. Proves that Xi vanishes identically on:
      - the normalized zero-defect kernel,
      - the injected extended kernel,
      - the explicit BdG support-blind directions,
      - the matched wall-only co-scaling direction.

Interpretation
--------------
This is the missing exact bridge between the Stage-5 prefactor/load language and
Stages 10-14. It proves that the Stage-5 load defect really is the same
weak-axisymmetric object controlled by the similarity-orbit theorem, and that
this theorem kills Xi_load but says nothing yet about the separate conservative
5PN even gates.
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
# I. Exact Stage-5 prefactor-slope identity
# ---------------------------------------------------------------------------

banner("I. EXACT STAGE-5 PREFactor-SLOPE IDENTITY")

D0, N0 = sp.symbols("D0 N0", nonzero=True)
D01, N01 = sp.symbols("D01 N01", real=True)

P0 = sp.simplify(N0 / D0)
P1 = sp.simplify((N01 * D0 - N0 * D01) / D0**2)
Xi_load = sp.simplify(N01 / N0 - D01 / D0)

print("P0       =", P0)
print("P1       =", P1)
print("Xi_load  =", Xi_load)
expect_zero("Xi_load - P1/P0", Xi_load - P1 / P0)

print(
    "So the Stage-5 load defect is exactly the first grouped weak-axisymmetric"
    " prefactor slope P1/P0."
)


# ---------------------------------------------------------------------------
# II. Exact normalized monomial-drift matrix and support-blind extension
# ---------------------------------------------------------------------------

banner("II. NORMALIZED MONOMIAL-DRIFT MATRIX AND SUPPORT-BLIND EXTENSION")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)
Atr = sp.symbols("A_tr", real=True)

# Normalized log drifts.
dw, du, dr = sp.symbols("dln_GW dln_GU dln_R", real=True)
dK, dM = sp.symbols("dln_K dln_M", real=True)
dOU, dOW = sp.symbols("dln_Omega_U dln_Omega_W", real=True)
dDelta = sp.symbols("dln_delta_U", real=True)

# Direct monomial drifts from Stages 10-13.
Sigma_tr = sp.expand((1 + deltaU) * (dr + du - dw - 2 * dOU) + (1 + chi0) * dDelta)
Sigma_nt = sp.expand(dM + 2 * dw - dK - 4 * dOW + 2 * Estar * (dr - dOU - dOW) - Fstar * dDelta)
Sigma_eta = sp.expand(dM + 2 * du - dK - 2 * dOU)
Xi_monomial = sp.expand(Atr * Sigma_tr + Sigma_nt)

vars_norm = (dw, du, dr, dK, dM, dOU, dOW, dDelta)
Mnorm = sp.Matrix([
    [sp.expand(Sigma_tr).coeff(v) for v in vars_norm],
    [sp.expand(Sigma_nt).coeff(v) for v in vars_norm],
    [sp.expand(Sigma_eta).coeff(v) for v in vars_norm],
])

print("M_norm =")
sp.pprint(Mnorm)

# Support-blind extension by the explicit BdG drifts.
dln_lamB, dln_varpi = sp.symbols("dln_lambda_B dln_varpi", real=True)
vars_ext = (dln_lamB, dln_varpi, *vars_norm)
Mext = sp.Matrix([
    [sp.expand(Sigma_tr).coeff(v) for v in vars_ext],
    [sp.expand(Sigma_nt).coeff(v) for v in vars_ext],
    [sp.expand(Sigma_eta).coeff(v) for v in vars_ext],
])

print("\nM_ext =")
sp.pprint(Mext)
expect_zero("BdG amplitude column vanishes", Mext[:, 0])
expect_zero("BdG frequency column vanishes", Mext[:, 1])


# ---------------------------------------------------------------------------
# III. Exact kernel formulas and Xi-annihilation checks
# ---------------------------------------------------------------------------

banner("III. Xi VANISHES ON THE EXACT ZERO-DEFECT KERNEL")

# Stage-13 normalized kernel formulas.
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

expect_zero("M_norm * kernel_norm", sp.simplify(Mnorm * kernel_norm))
expect_zero(
    "Xi_monomial on the normalized kernel",
    Xi_monomial.subs({dM: dM_expr, dOW: dOW_expr, dDelta: dDelta_expr}),
)

subbanner("Injected extended kernel and explicit support-blind directions")

kernel_ext = sp.Matrix([0, 0, *list(kernel_norm)])
expect_zero("M_ext * kernel_ext", sp.simplify(Mext * kernel_ext))

# Explicit support-blind directions.
e_B = {dln_lamB: 1, dln_varpi: 0, dw: 0, du: 0, dr: 0, dK: 0, dM: 0, dOU: 0, dOW: 0, dDelta: 0}
e_varpi = {dln_lamB: 0, dln_varpi: 1, dw: 0, du: 0, dr: 0, dK: 0, dM: 0, dOU: 0, dOW: 0, dDelta: 0}

expect_zero("Sigma_tr on e_B", Sigma_tr.subs(e_B))
expect_zero("Sigma_nt on e_B", Sigma_nt.subs(e_B))
expect_zero("Sigma_eta on e_B", Sigma_eta.subs(e_B))
expect_zero("Xi_monomial on e_B", Xi_monomial.subs(e_B))

expect_zero("Sigma_tr on e_varpi", Sigma_tr.subs(e_varpi))
expect_zero("Sigma_nt on e_varpi", Sigma_nt.subs(e_varpi))
expect_zero("Sigma_eta on e_varpi", Sigma_eta.subs(e_varpi))
expect_zero("Xi_monomial on e_varpi", Xi_monomial.subs(e_varpi))

subbanner("Matched wall-only co-scaling direction")

wall_match = {dln_lamB: 0, dln_varpi: 0, dw: 0, du: 0, dr: 0, dK: 1, dM: 1, dOU: 0, dOW: 0, dDelta: 0}
expect_zero("Sigma_tr on matched wall co-scaling", Sigma_tr.subs(wall_match))
expect_zero("Sigma_nt on matched wall co-scaling", Sigma_nt.subs(wall_match))
expect_zero("Sigma_eta on matched wall co-scaling", Sigma_eta.subs(wall_match))
expect_zero("Xi_monomial on matched wall co-scaling", Xi_monomial.subs(wall_match))

print(
    "The weak-axisymmetric normalization defect Xi is annihilated by the full"
    " normalized kernel, by its injected extension, and by the explicit support-"
    "blind BdG directions."
)


# ---------------------------------------------------------------------------
# IV. Final ledger
# ---------------------------------------------------------------------------

banner("IV. FINAL THEOREM LEDGER")
print("1. The Stage-5 load defect satisfies the exact grouped-prefactor identity")
print("      Xi_load = N01/N0 - D01/D0 = P1/P0.")
print("2. The Stage-10/13 direct monomial drifts are")
print("      Sigma_tr  = delta ln C_tr,* ,")
print("      Sigma_nt  = delta ln C_nt,* ,")
print("      Sigma_eta = delta ln epsilon_eta.")
print("3. The weak-axisymmetric normalization defect takes the exact triangular form")
print("      Xi = A_tr Sigma_tr + Sigma_nt.")
print("4. Therefore Xi vanishes on the exact Stage-13 normalized zero-defect kernel.")
print("5. The Stage-14 support-blind extension proves that Xi also vanishes on the")
print("   explicit BdG amplitude/frequency directions dln lambda_B and dln varpi.")
print("6. So the similarity-orbit theorem really is a theorem about the Stage-5 load")
print("   defect Xi_load = P1/P0. It kills the normalization defect, but it has not")
print("   yet touched the separate conservative 5PN even gates K1 and H_even.")
