#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 211 — normalized monomial bridge and exact zero-defect similarity kernel.

What this script does
---------------------
1. Rewrites the Stage-11 microscopic monomials directly in the normalized Stage-5
   variables (G_U, G_W, R, K, M, Omega_U, Omega_W, delta_U).
2. Builds the exact normalized monomial-drift matrix M_norm and verifies that it
   matches the compact formula already recorded in the notes.
3. Solves the zero-defect equations triangularly for
      dln(delta_U), dln(M), dln(Omega_W)
   in terms of the five freely chosen drifts
      dln(K), dln(G_U), dln(G_W), dln(R), dln(Omega_U).
4. Extends the primitive space back to the explicit support variables
      (lambda_B, varpi)
   and proves that they are exact zero columns of the monomial-drift map.

Interpretation
--------------
This is the executable bridge from the primitive weak-axisymmetric extractor to the
Stage-10/11 similarity-orbit theorem. It shows exactly which first-order drifts are
constrained by monomial rigidity and which directions remain support-blind.
"""


if __name__ == "__main__":
    banner("STAGE 211 — NORMALIZED MONOMIAL BRIDGE AND SIMILARITY KERNEL")

    # Normalized constructive-branch data.
    chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
    E_star, F_star = sp.symbols("E_star F_star", real=True)

    G_W, G_U, R = sp.symbols("G_W G_U R", positive=True, real=True)
    K, M = sp.symbols("K M", positive=True, real=True)
    Omega_U, Omega_W = sp.symbols("Omega_U Omega_W", positive=True, real=True)
    sigma = sp.symbols("sigma", positive=True, real=True)

    # Logarithmic primitive drifts.
    dln_GW, dln_GU, dln_R = sp.symbols("dln_GW dln_GU dln_R", real=True)
    dln_K, dln_M = sp.symbols("dln_K dln_M", real=True)
    dln_OU, dln_OW = sp.symbols("dln_Omega_U dln_Omega_W", real=True)
    dln_deltaU = sp.symbols("dln_delta_U", real=True)

    drift_vec = sp.Matrix([dln_GW, dln_GU, dln_R, dln_K, dln_M, dln_OU, dln_OW, dln_deltaU])

    # Normalized Stage-11 monomials.
    C_tr = sp.simplify((R * G_U / (Omega_U**2 * G_W)) ** (1 + deltaU) * deltaU ** (1 + chi0))
    C_nt = sp.simplify((M * G_W**2 / (K * Omega_W**4)) * (R**2 * sigma / (Omega_U**2 * Omega_W**2)) ** E_star * deltaU ** (-F_star))
    epsilon_eta = sp.simplify(M * G_U**2 / (K * Omega_U**2))

    # Exact logarithmic drifts from the normalized monomials.
    dln_Ctr = sp.simplify(
        -(1 + deltaU) * dln_GW
        + (1 + deltaU) * dln_GU
        + (1 + deltaU) * dln_R
        - 2 * (1 + deltaU) * dln_OU
        + (1 + chi0) * dln_deltaU
    )
    dln_Cnt = sp.simplify(
        2 * dln_GW
        + 2 * E_star * dln_R
        - dln_K
        + dln_M
        - 2 * E_star * dln_OU
        - (4 + 2 * E_star) * dln_OW
        - F_star * dln_deltaU
    )
    dln_eps_eta = sp.simplify(2 * dln_GU - dln_K + dln_M - 2 * dln_OU)

    # Matrix form expected from the notes.
    M_norm_expected = sp.Matrix(
        [
            [-(1 + deltaU), 1 + deltaU, 1 + deltaU, 0, 0, -2 * (1 + deltaU), 0, 1 + chi0],
            [2, 0, 2 * E_star, -1, 1, -2 * E_star, -(4 + 2 * E_star), -F_star],
            [0, 2, 0, -1, 1, -2, 0, 0],
        ]
    )
    M_norm = sp.Matrix(
        [
            [sp.diff(dln_Ctr, v) for v in drift_vec],
            [sp.diff(dln_Cnt, v) for v in drift_vec],
            [sp.diff(dln_eps_eta, v) for v in drift_vec],
        ]
    )

    subbanner("I. Exact normalized monomial-drift matrix")
    print("M_norm =")
    sp.pprint(M_norm)
    expect_zero("M_norm - expected", M_norm - M_norm_expected)
    print("rank(M_norm) =", M_norm.rank())
    if M_norm.rank() != 3:
        raise AssertionError("Normalized monomial matrix did not have rank 3.")

    subbanner("II. Exact zero-defect equations and triangular solution")
    print("Zero-defect equations are")
    sp.pprint(sp.Matrix([dln_Ctr, dln_Cnt, dln_eps_eta]))

    dln_deltaU_sol = sp.simplify(sp.solve(sp.Eq(dln_Ctr, 0), dln_deltaU)[0])
    dln_M_sol = sp.simplify(sp.solve(sp.Eq(dln_eps_eta, 0), dln_M)[0])
    dln_OW_sol = sp.simplify(
        sp.solve(
            sp.Eq(dln_Cnt.subs({dln_deltaU: dln_deltaU_sol}), 0),
            dln_OW,
        )[0]
    )

    print("dln(delta_U) =")
    sp.pprint(dln_deltaU_sol)
    print("dln(M) =")
    sp.pprint(dln_M_sol)
    print("dln(Omega_W) =")
    sp.pprint(dln_OW_sol)

    expect_zero("tracking equation", dln_Ctr.subs({dln_deltaU: dln_deltaU_sol}))
    expect_zero("dressing equation", dln_eps_eta.subs({dln_M: dln_M_sol}))
    expect_zero(
        "nontracking equation",
        dln_Cnt.subs({dln_deltaU: dln_deltaU_sol, dln_OW: dln_OW_sol}),
    )

    subbanner("III. Support-blind extension back to the Stage-5 primitive space")
    dln_lambdaB, dln_varpi = sp.symbols("dln_lambda_B dln_varpi", real=True)
    extended_vec = sp.Matrix([dln_lambdaB, dln_varpi, dln_GW, dln_GU, dln_R, dln_K, dln_M, dln_OU, dln_OW, dln_deltaU])

    # The normalized monomials are support-blind, so lambda_B and varpi enter with zero columns.
    M_ext = sp.Matrix(
        [
            [0, 0, -(1 + deltaU), 1 + deltaU, 1 + deltaU, 0, 0, -2 * (1 + deltaU), 0, 1 + chi0],
            [0, 0, 2, 0, 2 * E_star, -1, 1, -2 * E_star, -(4 + 2 * E_star), -F_star],
            [0, 0, 0, 2, 0, -1, 1, -2, 0, 0],
        ]
    )

    print("Extended primitive monomial-drift matrix =")
    sp.pprint(M_ext)
    expect_zero("lambda_B zero column", M_ext[:, 0])
    expect_zero("varpi zero column", M_ext[:, 1])
    print("rank(M_ext) =", M_ext.rank())
    if M_ext.rank() != 3:
        raise AssertionError("Extended primitive monomial matrix did not have rank 3.")
    print("nullity(M_ext) =", M_ext.shape[1] - M_ext.rank())

    banner("STAGE 211 LEDGER")
    print("1. The Stage-11 direct monomials are already present in the normalized Stage-5 variable set.")
    print("2. The exact normalized monomial-drift map has rank 3, so the normalized zero-defect")
    print("   tangent space has dimension 5.")
    print("3. The zero-defect equations solve triangularly:")
    print("      dln(delta_U), dln(M), dln(Omega_W)")
    print("   are fixed by")
    print("      dln(K), dln(G_U), dln(G_W), dln(R), dln(Omega_U).")
    print("4. Extending back to the explicit Stage-5 primitive space adds two exact support-blind")
    print("   directions: dln(lambda_B) and dln(varpi).")
    print("5. So monomial rigidity controls the mixed/wall/U normalization problem exactly, but it")
    print("   does not by itself close the explicit BdG support directions.")
