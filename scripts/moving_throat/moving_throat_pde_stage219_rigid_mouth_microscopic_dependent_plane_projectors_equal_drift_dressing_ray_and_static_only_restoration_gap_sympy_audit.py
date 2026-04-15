import sympy as sp


def assert_zero(expr: sp.Expr) -> None:
    simplified = sp.simplify(expr)
    if simplified != 0:
        raise AssertionError(f"Expression is not zero: {simplified}")


def assert_matrix_zero(mat: sp.Matrix) -> None:
    for entry in mat:
        simplified = sp.simplify(entry)
        if simplified != 0:
            raise AssertionError(f"Matrix entry is not zero: {simplified}")


def main() -> None:
    sp.init_printing()
    print(
        "=== Stage 219 SymPy audit: rigid-mouth microscopic dependent-plane projectors, equal-drift dressing ray, and static-only restoration gap ==="
    )

    # ------------------------------------------------------------------
    # 1. Carry-forward rigid-mouth packet and dependent section
    # ------------------------------------------------------------------
    chi0_star, F_star = sp.symbols("chi0_star F_star", positive=True, real=True)
    eps_eta_star = sp.symbols("eps_eta_star", positive=True, real=True)
    c_eta = eps_eta_star / (1 - eps_eta_star)

    q_tr, q_nt, q_eta = sp.symbols("q_tr q_nt q_eta", real=True)
    R1, E1 = sp.symbols("R1 E1", real=True)

    M_rm = sp.Matrix([[-1, -c_eta], [0, 1]])
    x_rm = sp.Matrix([R1, E1])
    q_rm = M_rm * x_rm
    assert_matrix_zero(q_rm - sp.Matrix([-R1 - c_eta * E1, E1]))

    Delta_T_q = q_tr / (1 + chi0_star)
    Delta_Keta_q = -q_eta
    Delta_mu_q = F_star * q_tr / (1 + chi0_star) + q_nt - q_eta

    y_general = sp.Matrix([Delta_T_q, Delta_Keta_q, Delta_mu_q])
    y_rm = y_general.subs(q_tr, 0)
    S_rm_dep = sp.Matrix([[0, 0], [0, -1], [1, -1]])

    assert_matrix_zero(y_rm - S_rm_dep * sp.Matrix([q_nt, q_eta]))

    print("\nRigid-mouth dependent section:")
    print("S_rm_dep =")
    sp.pprint(S_rm_dep)
    print("y_rm(q_nt,q_eta) =")
    sp.pprint(y_rm)

    # ------------------------------------------------------------------
    # 2. Direct-observable-to-dependent compiler
    # ------------------------------------------------------------------
    C_rm_dep = sp.simplify(S_rm_dep * M_rm)
    C_expected = sp.Matrix(
        [
            [0, 0],
            [0, -1],
            [-1, -1 / (1 - eps_eta_star)],
        ]
    )

    assert_matrix_zero(C_rm_dep - C_expected)

    y_from_x = sp.simplify(C_rm_dep * x_rm)
    assert_matrix_zero(
        y_from_x
        - sp.Matrix(
            [
                0,
                -E1,
                -R1 - E1 / (1 - eps_eta_star),
            ]
        )
    )

    print("\nDirect-observable-to-dependent compiler:")
    print("C_rm_dep =")
    sp.pprint(C_rm_dep)
    print("y_rm = C_rm_dep x_rm =")
    sp.pprint(y_from_x)

    # ------------------------------------------------------------------
    # 3. Exact dependent-plane packet projectors
    # ------------------------------------------------------------------
    L_rm_dep = sp.Matrix([[0, -1, 1], [0, -1, 0]])
    assert_matrix_zero(L_rm_dep * S_rm_dep - sp.eye(2))

    q_from_y = sp.simplify(L_rm_dep * y_rm)
    assert_matrix_zero(q_from_y - sp.Matrix([q_nt, q_eta]))

    Q_nt = sp.diag(1, 0)
    Q_eta = sp.diag(0, 1)

    P_nt_dep = sp.simplify(S_rm_dep * Q_nt * L_rm_dep)
    P_eta_dep = sp.simplify(S_rm_dep * Q_eta * L_rm_dep)

    P_nt_expected = sp.Matrix([[0, 0, 0], [0, 0, 0], [0, -1, 1]])
    P_eta_expected = sp.Matrix([[0, 0, 0], [0, 1, 0], [0, 1, 0]])

    assert_matrix_zero(P_nt_dep - P_nt_expected)
    assert_matrix_zero(P_eta_dep - P_eta_expected)
    assert_matrix_zero(P_nt_dep * P_nt_dep - P_nt_dep)
    assert_matrix_zero(P_eta_dep * P_eta_dep - P_eta_dep)
    assert_matrix_zero(P_nt_dep * P_eta_dep)
    assert_matrix_zero(P_eta_dep * P_nt_dep)

    plane_identity = sp.Matrix([[0, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert_matrix_zero(P_nt_dep + P_eta_dep - plane_identity)

    y_nt = sp.simplify(P_nt_dep * y_rm)
    y_eta = sp.simplify(P_eta_dep * y_rm)

    assert_matrix_zero(y_nt - sp.Matrix([0, 0, q_nt]))
    assert_matrix_zero(y_eta + q_eta * sp.Matrix([0, 1, 1]))
    assert_matrix_zero(y_rm - y_nt - y_eta)

    print("\nDependent-plane packet projectors:")
    print("L_rm_dep =")
    sp.pprint(L_rm_dep)
    print("P_nt_dep =")
    sp.pprint(P_nt_dep)
    print("P_eta_dep =")
    sp.pprint(P_eta_dep)
    print("y_nt =")
    sp.pprint(y_nt)
    print("y_eta =")
    sp.pprint(y_eta)

    # ------------------------------------------------------------------
    # 4. Static-strip equivalence and equal-drift dressing ray
    # ------------------------------------------------------------------
    # Static strip q_nt = 0 iff Delta_mu = Delta_Keta.
    assert_zero(sp.simplify((y_rm[2] - y_rm[1]) - q_nt))

    R1_static = -c_eta * E1
    y_static_strip = sp.simplify(y_from_x.subs(R1, R1_static))
    assert_matrix_zero(y_static_strip - sp.Matrix([0, -E1, -E1]))
    assert_matrix_zero(y_static_strip - (R1_static / c_eta) * sp.Matrix([0, 1, 1]))

    norm_sq_qeta = sp.simplify((y_eta.T * y_eta)[0])
    assert_zero(norm_sq_qeta - 2 * q_eta**2)

    norm_sq_E = sp.simplify((y_static_strip.T * y_static_strip)[0])
    assert_zero(norm_sq_E - 2 * E1**2)
    assert_zero(norm_sq_E - 2 * R1_static**2 / c_eta**2)

    print("\nStatic-strip and dressing-ray checks:")
    print("Delta_mu - Delta_Keta =", sp.simplify(y_rm[2] - y_rm[1]))
    print("On q_nt = 0, y_rm =")
    sp.pprint(y_static_strip)
    print("||y_eta||^2 =", norm_sq_qeta)

    # ------------------------------------------------------------------
    # 5. Microscopic correction compilers
    # ------------------------------------------------------------------
    Delta_y_static = -y_nt
    Delta_y_orbit = -y_rm

    assert_matrix_zero(Delta_y_static - sp.Matrix([0, 0, -q_nt]))
    assert_matrix_zero(y_rm + Delta_y_static - y_eta)
    assert_matrix_zero(
        Delta_y_orbit - (Delta_y_static + q_eta * sp.Matrix([0, 1, 1]))
    )
    assert_matrix_zero(y_rm + Delta_y_orbit)

    print("\nMicroscopic correction compilers:")
    print("Delta y_static =")
    sp.pprint(Delta_y_static)
    print("Delta y_orbit =")
    sp.pprint(Delta_y_orbit)
    print("Static-only and full orbit-lock correction identities verified exactly.")

    print("\nAll Stage 219 symbolic checks passed.")


if __name__ == "__main__":
    main()
