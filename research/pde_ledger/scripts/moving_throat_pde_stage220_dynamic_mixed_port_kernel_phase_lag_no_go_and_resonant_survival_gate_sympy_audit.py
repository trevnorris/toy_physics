import sympy as sp


def main() -> None:
    sp.init_printing()

    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    K, M, C, varpi = sp.symbols('K M C varpi', real=True, nonzero=True)
    Omega_U, Omega_W, R = sp.symbols('Omega_U Omega_W R', real=True, nonzero=True)
    G_U, G_W = sp.symbols('G_U G_W', real=True, nonzero=True)
    omega = sp.symbols('omega', real=True)
    Pi = sp.symbols('Pi')
    Gamma = sp.symbols('Gamma', real=True, nonnegative=True)

    J_q, J_U, J_W = sp.symbols('J_q J_U J_W', real=True)
    s_q, s_U, s_W, S = sp.symbols('s_q s_U s_W S', real=True)
    x, kappa = sp.symbols('x kappa', positive=True)
    beta_Q, beta_U, beta_W = sp.symbols('beta_Q beta_U beta_W', real=True)

    # ------------------------------------------------------------------
    # Carried static bundle data
    # ------------------------------------------------------------------
    K_star = K - C**2 / varpi**2
    Delta0 = Omega_U**2 * Omega_W**2 - R**2
    Q0 = G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2
    P_static = Omega_U**2 * G_W + R * G_U
    D0_static = sp.together(K_star - Q0 / Delta0)

    print('=== Stage 220 SymPy audit: dynamic mixed-port kernel ===')
    print('\nStatic carried data:')
    print('K_*   =', K_star)
    print('Delta =', Delta0)
    print('Q     =', Q0)
    print('P     =', P_static)
    print('D0    =', D0_static)

    # ------------------------------------------------------------------
    # Dynamic bundle
    # ------------------------------------------------------------------
    K_B = K - M * omega**2 - C**2 / (varpi**2 - omega**2)
    A = Omega_U**2 - omega**2
    Wfun = Omega_W**2 - omega**2 - Pi

    Delta_Pi = A * Wfun - R**2
    Q_Pi = G_U**2 * Wfun + 2 * G_U * G_W * R + G_W**2 * A
    D_Pi = sp.together(K_B - Q_Pi / Delta_Pi)

    print('\nDynamic bundle invariants:')
    print('K_B(omega)   =', K_B)
    print('A(omega)     =', A)
    print('W(omega)     =', Wfun)
    print('Delta_Pi     =', Delta_Pi)
    print('Q_Pi         =', Q_Pi)
    print('D_Pi         =', D_Pi)

    K_dyn = sp.Matrix([
        [K_B, -G_U, -G_W],
        [-G_U, A, -R],
        [-G_W, -R, Wfun],
    ])

    det_identity = K_dyn.det()
    assert sp.simplify(det_identity - Delta_Pi * D_Pi) == 0
    print('\nVerified determinant identity: det(K_dyn) = Delta_Pi * D_Pi')

    # Static reduction
    assert sp.simplify(K_B.subs({omega: 0}) - K_star) == 0
    assert sp.simplify(A.subs({omega: 0}) - Omega_U**2) == 0
    assert sp.simplify(Wfun.subs({omega: 0, Pi: 0}) - Omega_W**2) == 0
    assert sp.simplify(Delta_Pi.subs({omega: 0, Pi: 0}) - Delta0) == 0
    assert sp.simplify(Q_Pi.subs({omega: 0, Pi: 0}) - Q0) == 0
    assert sp.simplify(D_Pi.subs({omega: 0, Pi: 0}) - D0_static) == 0
    print('Verified static reduction back to the Stage 219 one-port bundle')

    # ------------------------------------------------------------------
    # Exact inverse entries
    # ------------------------------------------------------------------
    Kinv = K_dyn.inv()

    P_U = G_U * Wfun + R * G_W
    P = A * G_W + R * G_U

    chi_qq = 1 / D_Pi
    chi_qU = P_U / (Delta_Pi * D_Pi)
    chi_qW = P / (Delta_Pi * D_Pi)
    chi_UU = (K_B * Wfun - G_W**2) / (Delta_Pi * D_Pi)
    chi_UW = (K_B * R + G_U * G_W) / (Delta_Pi * D_Pi)
    chi_WW = (K_B * A - G_U**2) / (Delta_Pi * D_Pi)

    inverse_checks = {
        'chi_qq': sp.simplify(Kinv[0, 0] - chi_qq),
        'chi_qU': sp.simplify(Kinv[0, 1] - chi_qU),
        'chi_qW': sp.simplify(Kinv[0, 2] - chi_qW),
        'chi_UU': sp.simplify(Kinv[1, 1] - chi_UU),
        'chi_UW': sp.simplify(Kinv[1, 2] - chi_UW),
        'chi_WW': sp.simplify(Kinv[2, 2] - chi_WW),
    }
    assert all(v == 0 for v in inverse_checks.values())

    print('\nExact inverse entries:')
    print('chi_qq =', chi_qq)
    print('chi_qU =', chi_qU)
    print('chi_qW =', chi_qW)
    print('chi_UU =', chi_UU)
    print('chi_UW =', chi_UW)
    print('chi_WW =', chi_WW)

    # ------------------------------------------------------------------
    # Exact dynamic susceptibility law
    # ------------------------------------------------------------------
    J = sp.Matrix([J_q, J_U, J_W])
    deltaV_matrix = sp.together(-sp.Rational(1, 2) * (J.T * Kinv * J)[0])
    deltaV_formula = sp.together(-sp.Rational(1, 2) * (
        chi_qq * J_q**2
        + 2 * chi_qU * J_q * J_U
        + 2 * chi_qW * J_q * J_W
        + chi_UU * J_U**2
        + 2 * chi_UW * J_U * J_W
        + chi_WW * J_W**2
    ))
    assert sp.simplify(deltaV_matrix - deltaV_formula) == 0
    print('\nVerified dynamic susceptibility law:')
    print('V_mix(x,omega) = -1/2 J^T K_dyn^{-1} J')

    # Collinear-source theorem
    J_col = S * sp.Matrix([s_q, s_U, s_W])
    N_s = (
        Delta_Pi * s_q**2
        + 2 * P_U * s_q * s_U
        + 2 * P * s_q * s_W
        + (K_B * Wfun - G_W**2) * s_U**2
        + 2 * (K_B * R + G_U * G_W) * s_U * s_W
        + (K_B * A - G_U**2) * s_W**2
    )
    chi_s = sp.together(N_s / (Delta_Pi * D_Pi))
    deltaV_col = sp.together(-sp.Rational(1, 2) * (J_col.T * Kinv * J_col)[0])
    assert sp.simplify(deltaV_col + sp.Rational(1, 2) * chi_s * S**2) == 0
    print('\nVerified collinear-source theorem:')
    print('chi_s =', chi_s)

    # ------------------------------------------------------------------
    # Primitive-source product-family theorem
    # ------------------------------------------------------------------
    S_Q = x**-3
    S_Y = sp.exp(-2 * kappa * x) / x
    J_primitive = sp.Matrix([
        beta_Q * S_Q,
        beta_U * S_Y,
        beta_W * S_Y,
    ])

    deltaV_primitive = sp.together(-sp.Rational(1, 2) * (J_primitive.T * Kinv * J_primitive)[0])

    C6 = sp.together(chi_qq * beta_Q**2)
    C4 = sp.together(chi_qU * beta_Q * beta_U + chi_qW * beta_Q * beta_W)
    C2 = sp.together(chi_UU * beta_U**2 + 2 * chi_UW * beta_U * beta_W + chi_WW * beta_W**2)

    deltaV_primitive_expected = sp.together(
        -sp.Rational(1, 2) * (
            C6 / x**6
            + 2 * C4 * sp.exp(-2 * kappa * x) / x**4
            + C2 * sp.exp(-4 * kappa * x) / x**2
        )
    )
    assert sp.simplify(deltaV_primitive - deltaV_primitive_expected) == 0

    print('\nVerified primitive-source product-family theorem:')
    print('C6(omega) =', C6)
    print('C4(omega) =', C4)
    print('C2(omega) =', C2)

    # ------------------------------------------------------------------
    # Outgoing-port derivative identity
    # ------------------------------------------------------------------
    T_J = sp.together(chi_qW * J_q + chi_UW * J_U + chi_WW * J_W)
    dVdPi = sp.together(sp.diff(deltaV_formula, Pi))
    assert sp.simplify(dVdPi + sp.Rational(1, 2) * T_J**2) == 0

    print('\nVerified outgoing-port derivative identity:')
    print('dV_mix/dPi = -1/2 * T_J^2')
    print('T_J =', T_J)

    # ------------------------------------------------------------------
    # Linear outgoing correction and phase-lag theorem
    # ------------------------------------------------------------------
    T_J0 = sp.together(T_J.subs(Pi, 0))
    deltaV_linear_out = sp.together(dVdPi.subs(Pi, 0) * sp.I * Gamma)
    deltaV_linear_expected = sp.together(-sp.Rational(1, 2) * sp.I * Gamma * T_J0**2)
    assert sp.simplify(deltaV_linear_out - deltaV_linear_expected) == 0

    real_part, imag_part = sp.expand_complex(deltaV_linear_expected).as_real_imag()
    assert real_part == 0
    P_abs = sp.factor(-omega * imag_part)

    print('\nPhase-lag / pumping theorem:')
    print('deltaV_linear =', deltaV_linear_expected)
    print('Re(deltaV_linear) =', real_part)
    print('Im(deltaV_linear) =', imag_part)
    print('P_abs = -omega * Im(deltaV_linear) =', P_abs)

    # ------------------------------------------------------------------
    # Constructive off-pole slice
    # ------------------------------------------------------------------
    sample_subs = {
        K: sp.Integer(11),
        M: sp.Integer(2),
        C: sp.Integer(1),
        varpi: sp.Integer(5),
        Omega_U: sp.Integer(3),
        Omega_W: sp.Integer(4),
        R: sp.Integer(2),
        G_U: sp.Integer(1),
        G_W: sp.Integer(2),
        omega: sp.Rational(1, 2),
        J_q: sp.Integer(1),
        J_U: sp.Integer(2),
        J_W: sp.Integer(1),
    }

    Delta_sample = sp.N(Delta_Pi.subs({**sample_subs, Pi: 0}))
    D_sample = sp.N(D_Pi.subs({**sample_subs, Pi: 0}))
    V_cons_sample = sp.N(deltaV_formula.subs({**sample_subs, Pi: 0}))
    V_phase_sample = sp.N(deltaV_linear_expected.subs({**sample_subs, Gamma: sp.Rational(1, 10)}))
    P_abs_sample = sp.N(P_abs.subs({**sample_subs, Gamma: sp.Rational(1, 10)}))

    print('\nConstructive off-pole slice:')
    print('Delta_Pi(sample, Pi=0) =', Delta_sample)
    print('D_Pi(sample, Pi=0)     =', D_sample)
    print('V_mix(sample, Pi=0)    =', V_cons_sample)
    print('deltaV_linear(sample)  =', V_phase_sample)
    print('P_abs(sample)          =', P_abs_sample)

    assert Delta_sample != 0
    assert D_sample != 0
    assert abs(float(sp.re(V_phase_sample))) < 1e-12
    assert abs(float(sp.im(V_phase_sample))) > 0
    assert P_abs_sample > 0

    print('\nAll Stage 220 symbolic checks passed.')


if __name__ == '__main__':
    main()
