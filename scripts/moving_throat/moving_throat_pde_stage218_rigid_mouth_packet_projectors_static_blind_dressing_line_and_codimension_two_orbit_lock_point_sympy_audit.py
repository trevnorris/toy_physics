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
        "=== Stage 218 SymPy audit: rigid-mouth packet projectors, static-blind dressing line, and codimension-two orbit-lock point ==="
    )

    # ------------------------------------------------------------------
    # 1. Exact rigid-mouth packet compiler and involution
    # ------------------------------------------------------------------
    c_eta = sp.symbols("c_eta", positive=True, real=True)
    R1, E1 = sp.symbols("R1 E1", real=True)
    q_nt, q_eta = sp.symbols("q_nt q_eta", real=True)

    x_rm = sp.Matrix([R1, E1])
    M_rm = sp.Matrix([[-1, -c_eta], [0, 1]])
    q_rm = M_rm * x_rm
    Xi1 = -R1 - c_eta * E1

    assert_matrix_zero(q_rm - sp.Matrix([Xi1, E1]))
    assert_matrix_zero(M_rm * M_rm - sp.eye(2))

    x_from_q = M_rm * sp.Matrix([q_nt, q_eta])
    assert_matrix_zero(x_from_q - sp.Matrix([-q_nt - c_eta * q_eta, q_eta]))

    print("\nRigid-mouth packet compiler:")
    print("M_rm =")
    sp.pprint(M_rm)
    print("q_rm = M_rm x_rm =")
    sp.pprint(q_rm)
    print("M_rm^2 = I_2 verified exactly.")

    # ------------------------------------------------------------------
    # 2. Canonical packet projectors and direct-space decomposition
    # ------------------------------------------------------------------
    Q_nt = sp.diag(1, 0)
    Q_eta = sp.diag(0, 1)

    P_nt = M_rm * Q_nt * M_rm
    P_eta = M_rm * Q_eta * M_rm

    P_nt_expected = sp.Matrix([[1, c_eta], [0, 0]])
    P_eta_expected = sp.Matrix([[0, -c_eta], [0, 1]])

    assert_matrix_zero(P_nt - P_nt_expected)
    assert_matrix_zero(P_eta - P_eta_expected)
    assert_matrix_zero(P_nt * P_nt - P_nt)
    assert_matrix_zero(P_eta * P_eta - P_eta)
    assert_matrix_zero(P_nt * P_eta)
    assert_matrix_zero(P_eta * P_nt)
    assert_matrix_zero(P_nt + P_eta - sp.eye(2))

    x_nt = P_nt * x_rm
    x_eta = P_eta * x_rm

    assert_matrix_zero(x_nt - sp.Matrix([R1 + c_eta * E1, 0]))
    assert_matrix_zero(x_nt - sp.Matrix([-Xi1, 0]))
    assert_matrix_zero(x_eta - sp.Matrix([-c_eta * E1, E1]))
    assert_matrix_zero(x_rm - x_nt - x_eta)

    print("\nDirect-space packet projectors:")
    print("P_nt =")
    sp.pprint(P_nt)
    print("P_eta =")
    sp.pprint(P_eta)
    print("x_nt =")
    sp.pprint(x_nt)
    print("x_eta =")
    sp.pprint(x_eta)

    # ------------------------------------------------------------------
    # 3. Codimension-two orbit-lock equivalences
    # ------------------------------------------------------------------
    # q_rm = 0 iff x_rm = 0 because M_rm is involutive/invertible.
    assert_zero(M_rm.det() + 1)  # det = -1

    sol_xi_E = sp.solve([sp.Eq(Xi1, 0), sp.Eq(E1, 0)], [R1, E1], dict=True)
    if sol_xi_E != [{R1: 0, E1: 0}]:
        raise AssertionError(f"Unexpected solution for Xi1=0 and E1=0: {sol_xi_E}")

    sol_xi_R = sp.solve([sp.Eq(Xi1, 0), sp.Eq(R1, 0)], [R1, E1], dict=True)
    if sol_xi_R != [{R1: 0, E1: 0}]:
        raise AssertionError(f"Unexpected solution for Xi1=0 and R1=0: {sol_xi_R}")

    print("\nCodimension-two orbit-lock checks:")
    print("det(M_rm) =", sp.simplify(M_rm.det()))
    print("q_rm = 0 iff R1 = 0 and E1 = 0.")
    print("Also Xi1 = 0 and E1 = 0 -> (R1,E1) = (0,0).")
    print("And Xi1 = 0 and R1 = 0 -> (R1,E1) = (0,0).")

    # ------------------------------------------------------------------
    # 4. Static-blind dressing line and exact norm
    # ------------------------------------------------------------------
    qeta_line = sp.symbols("qeta_line", real=True)
    x_blind = sp.Matrix([-c_eta * qeta_line, qeta_line])
    Xi1_blind = sp.simplify((-x_blind[0] - c_eta * x_blind[1]))
    assert_zero(Xi1_blind)

    norm_sq = sp.simplify((x_blind.T * x_blind)[0])
    assert_zero(norm_sq - (1 + c_eta**2) * qeta_line**2)

    L = sp.symbols("L", positive=True, real=True)
    qeta_choice = L / sp.sqrt(1 + c_eta**2)
    norm_sq_choice = sp.simplify(norm_sq.subs(qeta_line, qeta_choice))
    assert_zero(norm_sq_choice - L**2)

    print("\nStatic-blind dressing line:")
    print("x_eta(q_eta) =")
    sp.pprint(x_blind)
    print("Xi1 on the line =", Xi1_blind)
    print("||x_eta||^2 =", norm_sq)
    print("Choosing q_eta = L/sqrt(1+c_eta^2) gives ||x_eta|| = L exactly.")

    # ------------------------------------------------------------------
    # 5. Canonical correction compilers
    # ------------------------------------------------------------------
    q_vec = sp.Matrix([q_nt, q_eta])
    x_q = M_rm * q_vec
    x_nt_from_q = M_rm * Q_nt * q_vec
    x_eta_from_q = M_rm * Q_eta * q_vec

    Delta_x_static = -x_nt_from_q
    Delta_x_orbit = -x_q

    assert_matrix_zero(x_q - sp.Matrix([-q_nt - c_eta * q_eta, q_eta]))
    assert_matrix_zero(x_nt_from_q - sp.Matrix([-q_nt, 0]))
    assert_matrix_zero(x_eta_from_q - sp.Matrix([-c_eta * q_eta, q_eta]))
    assert_matrix_zero(Delta_x_static - sp.Matrix([q_nt, 0]))
    assert_matrix_zero(x_q + Delta_x_static - x_eta_from_q)
    assert_matrix_zero(Delta_x_orbit - sp.Matrix([q_nt + c_eta * q_eta, -q_eta]))
    assert_matrix_zero(
        Delta_x_orbit - (Delta_x_static + sp.Matrix([c_eta * q_eta, -q_eta]))
    )
    assert_matrix_zero(x_q + Delta_x_orbit)

    print("\nCanonical correction compilers:")
    print("Delta x_static =")
    sp.pprint(Delta_x_static)
    print("Delta x_orbit =")
    sp.pprint(Delta_x_orbit)
    print("Static-only and full orbit-lock correction identities verified exactly.")

    print("\nAll Stage 218 symbolic checks passed.")


if __name__ == "__main__":
    main()
