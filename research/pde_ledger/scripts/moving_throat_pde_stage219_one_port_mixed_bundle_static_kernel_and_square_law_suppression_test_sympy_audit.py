import sympy as sp


def main() -> None:
    sp.init_printing()

    # ------------------------------------------------------------------
    # Symbols and bundle invariants
    # ------------------------------------------------------------------
    K_star, Omega_U, Omega_W, R = sp.symbols(
        'K_star Omega_U Omega_W R', nonzero=True
    )
    G_U, G_W = sp.symbols('G_U G_W', nonzero=True)
    J_q, J_U, J_W = sp.symbols('J_q J_U J_W')
    s_q, s_U, s_W, S = sp.symbols('s_q s_U s_W S')
    x, kappa = sp.symbols('x kappa', positive=True)
    beta_Q, beta_U, beta_W = sp.symbols('beta_Q beta_U beta_W')

    Delta = Omega_U**2 * Omega_W**2 - R**2
    Q = G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2
    P = Omega_U**2 * G_W + R * G_U
    P_U = G_U * Omega_W**2 + R * G_W
    D0 = sp.simplify(K_star - Q / Delta)
    D0_num = sp.expand(K_star * Delta - Q)

    print('=== Stage 219 SymPy audit: one-port mixed-bundle static kernel ===')
    print('\nBundle invariants:')
    print('Delta =', Delta)
    print('Q     =', Q)
    print('P     =', P)
    print('P_U   =', P_U)
    print('D0    =', D0)

    # ------------------------------------------------------------------
    # Reduced static bundle and determinant identity
    # ------------------------------------------------------------------
    K_red = sp.Matrix([
        [K_star, -G_U, -G_W],
        [-G_U, Omega_U**2, -R],
        [-G_W, -R, Omega_W**2],
    ])

    det_identity = sp.factor(K_red.det())
    print('\nDeterminant of K_red:')
    print(det_identity)
    assert sp.simplify(det_identity - Delta * D0) == 0

    # Internal block and Schur complement on the admissible branch
    K_int = sp.Matrix([
        [Omega_U**2, -R],
        [-R, Omega_W**2],
    ])
    coupling = sp.Matrix([[-G_U, -G_W]])
    schur = sp.simplify(K_star - (coupling * K_int.inv() * coupling.T)[0])
    print('\nInternal block determinant:')
    print(sp.factor(K_int.det()))
    print('Schur complement of the internal (U,W) block:')
    print(sp.factor(schur))
    assert sp.simplify(schur - D0) == 0

    # ------------------------------------------------------------------
    # Exact inverse entries
    # ------------------------------------------------------------------
    Kinv = sp.simplify(K_red.inv())

    chi_qq_expected = 1 / D0
    chi_qU_expected = P_U / (Delta * D0)
    chi_qW_expected = P / (Delta * D0)
    chi_UU_expected = (K_star * Omega_W**2 - G_W**2) / (Delta * D0)
    chi_UW_expected = (K_star * R + G_U * G_W) / (Delta * D0)
    chi_WW_expected = (K_star * Omega_U**2 - G_U**2) / (Delta * D0)

    inverse_expectations = {
        'chi_qq': (Kinv[0, 0], chi_qq_expected),
        'chi_qU': (Kinv[0, 1], chi_qU_expected),
        'chi_qW': (Kinv[0, 2], chi_qW_expected),
        'chi_UU': (Kinv[1, 1], chi_UU_expected),
        'chi_UW': (Kinv[1, 2], chi_UW_expected),
        'chi_WW': (Kinv[2, 2], chi_WW_expected),
    }

    print('\nExact inverse entries:')
    for label, (actual, expected) in inverse_expectations.items():
        diff = sp.simplify(sp.together(actual - expected))
        assert diff == 0
        print(f'{label} = {sp.simplify(expected)}')

    # ------------------------------------------------------------------
    # Exact static susceptibility kernel
    # ------------------------------------------------------------------
    J = sp.Matrix([J_q, J_U, J_W])
    delta_V_mix = sp.simplify(-sp.Rational(1, 2) * (J.T * Kinv * J)[0])
    print('\nGeneral on-shell static energy shift delta_V_mix:')
    print(delta_V_mix)

    # Collinear-source theorem
    J_col = S * sp.Matrix([s_q, s_U, s_W])
    delta_V_col = sp.simplify(-sp.Rational(1, 2) * (J_col.T * Kinv * J_col)[0])

    N_s = sp.expand(
        Delta * s_q**2
        + 2 * P_U * s_q * s_U
        + 2 * P * s_q * s_W
        + (K_star * Omega_W**2 - G_W**2) * s_U**2
        + 2 * (K_star * R + G_U * G_W) * s_U * s_W
        + (K_star * Omega_U**2 - G_U**2) * s_W**2
    )
    chi_s = sp.simplify(N_s / (Delta * D0))
    delta_V_col_expected = sp.simplify(-sp.Rational(1, 2) * chi_s * S**2)
    assert sp.simplify(sp.together(delta_V_col - delta_V_col_expected)) == 0

    print('\nCollinear-source theorem:')
    print('chi_s =', chi_s)
    print('delta_V_mix(collinear) =', delta_V_col_expected)

    # ------------------------------------------------------------------
    # Bridge to outgoing prefactor data
    # ------------------------------------------------------------------
    Lambda = sp.simplify(P / Delta)
    N0 = sp.simplify(Lambda**2)
    P0 = sp.simplify(N0 / D0)
    chi_qW = sp.simplify(Kinv[0, 2])

    print('\nOutgoing-prefactor bridge:')
    print('Lambda =', Lambda)
    print('N0     =', N0)
    print('P0     =', P0)
    print('chi_qW =', chi_qW)

    assert sp.simplify(sp.together(chi_qW - Lambda / D0)) == 0
    assert sp.simplify(sp.together(chi_qW**2 - P0 / D0)) == 0
    print('Verified: chi_qW = Lambda / D0 and chi_qW^2 = P0 / D0')

    # ------------------------------------------------------------------
    # Product-kernel theorem for primitive source families
    # ------------------------------------------------------------------
    S_Q = x**-3
    S_Y = sp.exp(-2 * kappa * x) / x
    J_primitive = sp.Matrix([
        beta_Q * S_Q,
        beta_U * S_Y,
        beta_W * S_Y,
    ])

    delta_V_primitive = sp.simplify(-sp.Rational(1, 2) * (J_primitive.T * Kinv * J_primitive)[0])

    C6 = sp.simplify(chi_qq_expected * beta_Q**2)
    C4 = sp.simplify(chi_qU_expected * beta_Q * beta_U + chi_qW_expected * beta_Q * beta_W)
    C2 = sp.simplify(
        chi_UU_expected * beta_U**2
        + 2 * chi_UW_expected * beta_U * beta_W
        + chi_WW_expected * beta_W**2
    )

    delta_V_primitive_expected = sp.simplify(
        -sp.Rational(1, 2)
        * (
            C6 / x**6
            + 2 * C4 * sp.exp(-2 * kappa * x) / x**4
            + C2 * sp.exp(-4 * kappa * x) / x**2
        )
    )

    assert sp.simplify(sp.together(delta_V_primitive - delta_V_primitive_expected)) == 0

    print('\nPrimitive-source product-kernel theorem:')
    print('C6 =', C6)
    print('C4 =', C4)
    print('C2 =', C2)
    print('delta_V_mix(primitive) =')
    print(delta_V_primitive_expected)

    # ------------------------------------------------------------------
    # A simple constructive admissible slice
    # ------------------------------------------------------------------
    sample_subs = {
        K_star: sp.Integer(11),
        Omega_U: sp.Integer(3),
        Omega_W: sp.Integer(4),
        R: sp.Integer(2),
        G_U: sp.Integer(1),
        G_W: sp.Integer(2),
    }
    Delta_sample = sp.simplify(Delta.subs(sample_subs))
    D0_sample = sp.simplify(D0.subs(sample_subs))
    det_sample = sp.simplify(det_identity.subs(sample_subs))
    K_red_sample = K_red.subs(sample_subs)
    minor_11 = sp.simplify(K_red_sample[0, 0])
    minor_22 = sp.simplify(K_red_sample[:2, :2].det())
    eigenvalues_numeric = [complex(sp.N(ev)) for ev in K_red_sample.eigenvals().keys()]

    print('\nConstructive admissible slice (for sanity check):')
    print('Delta(sample) =', Delta_sample)
    print('D0(sample)    =', D0_sample)
    print('det(sample)   =', det_sample)
    print('leading principal minors(sample) =', minor_11, minor_22, det_sample)
    print('eigenvalues(sample, numeric) =', eigenvalues_numeric)

    assert Delta_sample > 0
    assert D0_sample > 0
    assert minor_11 > 0 and minor_22 > 0 and det_sample > 0
    assert all(ev.real > 0 and abs(ev.imag) < 1e-12 for ev in eigenvalues_numeric)

    print('\nAll Stage 219 symbolic checks passed.')


if __name__ == '__main__':
    main()
