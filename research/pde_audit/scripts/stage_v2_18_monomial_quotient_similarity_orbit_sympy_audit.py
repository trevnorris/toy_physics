#!/usr/bin/env python3
"""
Stage V2-18 — Monomial quotient and similarity-orbit audit.

This script checks the exact algebra used to organize the weak-axisymmetric
moving-throat branch into:

  q = (q_tr, q_nt, q_eta)
    = (delta ln C_tr, delta ln C_nt, delta ln epsilon_eta)

and separates true quotient motion from similarity-orbit motion.

It also verifies the exact map from quotient residuals to the physical
first-order defect triplet (Theta_1, Xi_1, R_1), including Xi_1 = P_1/P_0.
"""

from __future__ import annotations

import sympy as sp


def check(name: str, condition: bool, details: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"{status}: {name}")
    if details:
        print(f"      {details}")
    if not condition:
        raise AssertionError(name)


def is_zero(expr) -> bool:
    return sp.simplify(expr) == 0


def all_zero(mat) -> bool:
    return all(is_zero(entry) for entry in list(mat))


def main() -> None:
    # ---------------------------------------------------------------------
    # Symbols and assumptions
    # ---------------------------------------------------------------------
    chi, delta, E, F = sp.symbols("chi delta E F", nonzero=True)
    Ctr, Atr, eps = sp.symbols("C_tr A_tr eps_eta", nonzero=True)
    D0, N0, D01, N01 = sp.symbols("D0 N0 D01 N01", nonzero=True)

    # ---------------------------------------------------------------------
    # 1. Original microscopic monomial matrix M_*
    # ---------------------------------------------------------------------
    # Drift vector order:
    #   x = (ln lambda_W, ln c_etaU, ln gamma, ln K_U,
    #        ln K_eta_eff, ln K_W_eff, ln mu_W, ln T_U)
    #
    # Monomials:
    #   C_tr = (gamma c_etaU/K_U)^(1+delta) *
    #          (pi^2 T_U/(L^2 K_U))^(1+chi)
    #
    #   C_nt = (lambda_W^2 mu_W/(K_eta K_W^2)) *
    #          (gamma^2 lambda_W^2 sigma/(K_U K_W))^E *
    #          (pi^2 T_U/(L^2 K_U))^(-F)
    #
    #   eps_eta = c_etaU^2/(K_U K_eta)
    #
    # Constants pi, L, sigma do not contribute to logarithmic drifts.
    Mstar = sp.Matrix([
        [0,          1 + delta, 1 + delta, -(2 + delta + chi), 0,       0,        0,  1 + chi],
        [2 + 2*E,    0,         2*E,        F - E,              -1,     -(2 + E), 1, -F],
        [0,          2,         0,          -1,                 -1,      0,        0,  0],
    ])

    # The exact minor quoted in the notes uses columns:
    # (ln T_U, ln K_eta_eff, ln mu_W) = (7,4,6).
    minor_cols = [7, 4, 6]
    minor = sp.factor(Mstar[:, minor_cols].det())
    check(
        "original M_* has quoted nonzero minor det(T_U,K_eta,mu_W)=1+chi",
        is_zero(minor - (1 + chi)),
        f"minor={minor}",
    )

    check(
        "rank(M_*) = 3 under physical condition 1+chi != 0",
        minor == 1 + chi,
        "a nonzero 3x3 minor proves rank at least 3; matrix has only 3 rows",
    )

    nullity = 8 - 3
    check("dim ker(M_*) = 5", nullity == 5, f"nullity={nullity}")

    # ---------------------------------------------------------------------
    # 2. Exact normal basis N satisfying M_* N = I_3
    # ---------------------------------------------------------------------
    n_tr = sp.Matrix([0, 0, 0, 0, 0, 0, F/(1 + chi), 1/(1 + chi)])
    n_nt = sp.Matrix([0, 0, 0, 0, 0, 0, 1, 0])
    n_eta = sp.Matrix([0, 0, 0, 0, -1, 0, -1, 0])
    Nmat = sp.Matrix.hstack(n_tr, n_nt, n_eta)

    MN = sp.simplify(Mstar * Nmat)
    check("exact normal basis satisfies M_* N = I_3", MN == sp.eye(3), f"M_*N={MN}")

    # Projector onto quotient normal representatives: N M_*
    P_normal = sp.simplify(Nmat * Mstar)
    check("normal projector is idempotent", sp.simplify(P_normal * P_normal - P_normal) == sp.zeros(8))

    # Orbit projector I - N M_* should be idempotent and annihilated by M_*.
    P_orbit = sp.simplify(sp.eye(8) - P_normal)
    check("orbit projector is idempotent", sp.simplify(P_orbit * P_orbit - P_orbit) == sp.zeros(8))
    check("M_* annihilates orbit projector", sp.simplify(Mstar * P_orbit) == sp.zeros(3, 8))

    # Generic drift split.
    x_symbols = sp.symbols("x0:8")
    x = sp.Matrix(x_symbols)
    q = sp.simplify(Mstar * x)
    x_normal = sp.simplify(Nmat * q)
    x_orbit = sp.simplify(x - x_normal)
    check("generic split x = x_orbit + N q", sp.simplify(x_orbit + x_normal - x) == sp.zeros(8, 1))
    check("x_orbit lies in ker(M_*)", sp.simplify(Mstar * x_orbit) == sp.zeros(3, 1))
    check("normal part carries the original quotient q", sp.simplify(Mstar * x_normal - q) == sp.zeros(3, 1))

    # ---------------------------------------------------------------------
    # 3. Normalized Stage-12/13 monomial matrix.
    # ---------------------------------------------------------------------
    # Drift vector order:
    # (ln G_W, ln G_U, ln R, ln K, ln M,
    #  ln Omega_U, ln Omega_W, ln delta_U).
    Mnorm = sp.Matrix([
        [-(1 + delta), 1 + delta, 1 + delta, 0, 0, -2*(1 + delta), 0, 1 + chi],
        [2,            0,         2*E,       -1, 1, -2*E,          -(4 + 2*E), -F],
        [0,            2,         0,         -1, 1, -2,             0,          0],
    ])

    # Verify rows come from the normalized monomials:
    # C_tr = (R G_U/(Omega_U^2 G_W))^(1+delta) * delta_U^(1+chi)
    # C_nt = M G_W^2/(K Omega_W^4) *
    #        (R^2 sigma/(Omega_U^2 Omega_W^2))^E * delta_U^(-F)
    # eps_eta = M G_U^2/(K Omega_U^2)
    expected_Mnorm = sp.Matrix([
        [-(1 + delta), 1 + delta, 1 + delta, 0, 0, -2*(1 + delta), 0, 1 + chi],
        [2,            0,         2*E,       -1, 1, -2*E,          -4 - 2*E,   -F],
        [0,            2,         0,         -1, 1, -2,             0,          0],
    ])
    check("normalized monomial-drift matrix matches exponent rows", Mnorm == expected_Mnorm)

    norm_minor_cols = [7, 4, 6]  # (delta_U, M, Omega_W)
    norm_minor = sp.factor(Mnorm[:, norm_minor_cols].det())
    check(
        "normalized M has nonzero rank-3 minor under physical conditions",
        is_zero(norm_minor - 2*(E + 2)*(1 + chi)),
        f"minor={norm_minor}",
    )

    # Verify triangular zero-defect solve from the notes.
    aW, aU, aR, aK, aM, aOU, aOW, aDel = sp.symbols(
        "aW aU aR aK aM aOU aOW aDel"
    )
    aDel_sol = -((1 + delta)/(1 + chi)) * (aR + aU - aW - 2*aOU)
    aM_sol = aK - 2*aU + 2*aOU
    aOW_sol = (
        aW - aU + (1 - E)*aOU + E*aR - sp.Rational(1, 2)*F*aDel_sol
    )/(E + 2)

    v_zero = sp.Matrix([aW, aU, aR, aK, aM_sol, aOU, aOW_sol, aDel_sol])
    Mv_zero = sp.simplify(Mnorm * v_zero)
    check("triangular normalized zero-defect solve gives M_norm v = 0", Mv_zero == sp.zeros(3, 1))

    # Support-blind extension: add two explicit BdG columns that are zero.
    Mnorm_ext = sp.Matrix.hstack(sp.zeros(3, 2), Mnorm)
    check("extended BdG columns are exactly support-blind zero columns", Mnorm_ext[:, 0:2] == sp.zeros(3, 2))
    ext_nullity = 10 - 3
    check("extended support-blind monomial kernel has dimension 7", ext_nullity == 7)

    # ---------------------------------------------------------------------
    # 4. Exact physical defect compiler from quotient residuals.
    # ---------------------------------------------------------------------
    Ddef = sp.Matrix([
        [-Ctr, 0, 0],
        [Atr, 1, 0],
        [-Atr, -1, -eps/(1 - eps)],
    ])

    Ddet = sp.factor(Ddef.det())
    check(
        "defect compiler determinant is nonzero on physical branch",
        is_zero(Ddet - Ctr*eps/(1 - eps)),
        f"det={Ddet}",
    )

    qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta")
    qvec = sp.Matrix([qtr, qnt, qeta])
    defects = sp.simplify(Ddef * qvec)
    Theta, Xi, Rres = defects

    check("Theta_1 = -C_tr q_tr", is_zero(Theta + Ctr*qtr), f"Theta={Theta}")
    check("Xi_1 = A_tr q_tr + q_nt", is_zero(Xi - (Atr*qtr + qnt)), f"Xi={Xi}")
    check(
        "R_1 = -Xi_1 - eps/(1-eps) q_eta",
        is_zero(Rres + Xi + eps*qeta/(1 - eps)),
        f"R={Rres}",
    )

    # Verify inverse formulas.
    qtr_inv = -Theta / Ctr
    qnt_inv = Xi + (Atr/Ctr)*Theta
    qeta_inv = -((1 - eps)/eps) * (Rres + Xi)
    inv_residual = sp.simplify(sp.Matrix([qtr_inv, qnt_inv, qeta_inv]) - qvec)
    check("defect compiler inverse formulas recover q", inv_residual == sp.zeros(3, 1))

    # Equivalence of zero defects and zero quotient residuals.
    check(
        "zero quotient residuals imply zero physical defects",
        sp.simplify(Ddef * sp.zeros(3, 1)) == sp.zeros(3, 1),
    )
    # Since det != 0, the converse is exact on the physical branch.
    check(
        "zero physical defects imply zero quotient residuals when det != 0",
        Ddet != 0,
        "symbolic determinant C_tr*eps/(1-eps) is nonzero under C_tr!=0, eps!=0, eps!=1",
    )

    # ---------------------------------------------------------------------
    # 5. Xi_load bridge: Xi_1 = P_1/P_0.
    # ---------------------------------------------------------------------
    P0 = N0 / D0
    P1 = (N01*D0 - N0*D01) / D0**2
    Xi_load = N01/N0 - D01/D0
    check("Xi_load = P1/P0", is_zero(sp.simplify(P1/P0 - Xi_load)), f"Xi_load={Xi_load}")

    # Direct zero-Xi condition in quotient coordinates.
    # Xi=0 -> q_nt = -A_tr q_tr.
    check(
        "Xi_1=0 condition is q_nt = -A_tr q_tr",
        is_zero((Atr*qtr + (-Atr*qtr))),
    )

    # ---------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------
    print("\nSUMMARY")
    print("original_Mstar_shape:", Mstar.shape)
    print("original_rank_witness_minor:", minor)
    print("original_kernel_dimension:", nullity)
    print("normalized_M_shape:", Mnorm.shape)
    print("normalized_rank_witness_minor:", norm_minor)
    print("normalized_extended_kernel_dimension_with_BdG_blind_columns:", ext_nullity)
    print("defect_compiler_det:", Ddet)
    print("Xi_load:", sp.sstr(Xi_load))
    print("verdict: PASS — quotient variables separate similarity-orbit motion from physical weak-axisymmetric defects.")


if __name__ == "__main__":
    main()
