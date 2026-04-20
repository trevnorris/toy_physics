import sympy as sp


def main() -> None:
    sp.init_printing()

    print('=== Stage 204 SymPy audit: resonance / linewidth tradeoff ===')

    # ------------------------------------------------------------------
    # Generic simple-pole normal form
    # ------------------------------------------------------------------
    delta = sp.symbols('delta', real=True)
    Gamma_out, Z_star, Fprime, Num_star = sp.symbols(
        'Gamma_out Z_star Fprime Num_star', positive=True
    )
    Pi = sp.symbols('Pi')

    F_lin = Fprime * delta - Pi * Z_star
    chi_lin = sp.together(Num_star / F_lin)

    A_star = Num_star / Fprime
    gamma_star = Gamma_out * Z_star / Fprime
    chi_bw = sp.together(A_star / (delta - Pi * Z_star / Fprime))
    assert sp.simplify(chi_lin - chi_bw) == 0

    chi_passive = sp.together(chi_lin.subs(Pi, sp.I * Gamma_out))
    chi_passive_expected = sp.together(A_star / (delta - sp.I * gamma_star))
    assert sp.simplify(chi_passive - chi_passive_expected) == 0

    print('\nVerified generic simple-pole normal form:')
    print('chi(omega) =', chi_passive_expected)
    print('A_*        =', sp.simplify(A_star))
    print('gamma_*    =', sp.simplify(gamma_star))

    # ------------------------------------------------------------------
    # Exact Stage-203 wall derivative identity
    # ------------------------------------------------------------------
    K, M, C, varpi = sp.symbols('K M C varpi', real=True, nonzero=True)
    Omega_U, Omega_W, R = sp.symbols('Omega_U Omega_W R', real=True, nonzero=True)
    G_U, G_W = sp.symbols('G_U G_W', real=True, nonzero=True)
    omega = sp.symbols('omega', real=True)

    K_B = K - M * omega**2 - C**2 / (varpi**2 - omega**2)
    A = Omega_U**2 - omega**2
    Wfun = Omega_W**2 - omega**2 - Pi

    Delta_Pi = A * Wfun - R**2
    Q_Pi = G_U**2 * Wfun + 2 * G_U * G_W * R + G_W**2 * A
    D_Pi = sp.together(K_B - Q_Pi / Delta_Pi)
    Nfun = sp.together((A * G_W + R * G_U) ** 2 / Delta_Pi**2)

    assert sp.simplify(sp.diff(D_Pi, Pi) + Nfun) == 0
    print('\nVerified exact Stage-203 derivative identity:')
    print('dD_Pi/dPi = -N(omega)')

    # Wall-like pole normal form on the sign-fixed slice D0prime > 0.
    D0prime, N_star = sp.symbols('D0prime N_star', positive=True)
    D_lin = D0prime * delta - Pi * N_star
    chiqq_lin = sp.together(1 / D_lin)
    gamma_wall = Gamma_out * N_star / D0prime
    chiqq_expected = sp.together((1 / D0prime) / (delta - sp.I * gamma_wall))
    assert sp.simplify(chiqq_lin.subs(Pi, sp.I * Gamma_out) - chiqq_expected) == 0

    print('\nVerified wall-like pole specialization:')
    print('chi_qq(omega) =', chiqq_expected)
    print('gamma_wall    =', sp.simplify(gamma_wall))

    # ------------------------------------------------------------------
    # Exact line-shape algebra on the positive-detuning slice delta = r*gamma
    # ------------------------------------------------------------------
    Aabs, gamma, r, eta = sp.symbols('Aabs gamma r eta', positive=True)
    chi_r = sp.together(Aabs / (r * gamma - sp.I * gamma))
    re_r, im_r = sp.expand_complex(chi_r).as_real_imag()

    re_expected = sp.simplify(Aabs * r / (gamma * (1 + r**2)))
    im_expected = sp.simplify(Aabs / (gamma * (1 + r**2)))
    assert sp.simplify(re_r - re_expected) == 0
    assert sp.simplify(im_r - im_expected) == 0
    assert sp.simplify(re_r / im_r - r) == 0

    print('\nExact line-shape components on the positive-detuning slice:')
    print('Re chi =', re_expected)
    print('Im chi =', im_expected)
    print('|Re|/|Im| =', sp.simplify(re_expected / im_expected))

    # ------------------------------------------------------------------
    # Exact resonance optimum and low-loss bound
    # ------------------------------------------------------------------
    f = sp.simplify(r / (1 + r**2))
    max_identity = sp.factor(sp.together(f - sp.Rational(1, 2)))
    assert sp.simplify(max_identity + (r - 1)**2 / (2 * (1 + r**2))) == 0

    low_loss_identity = sp.factor(
        sp.together(f - eta / (1 + eta**2))
    )
    assert sp.simplify(low_loss_identity + (r - eta) * (eta * r - 1) / ((1 + r**2) * (1 + eta**2))) == 0

    print('\nExact conservative-shape factor: f(r) =', f)
    print('Maximum identity: f(r) - 1/2 =', max_identity)
    print('Low-loss identity: f(r) - eta/(1+eta^2) =', low_loss_identity)

    # ------------------------------------------------------------------
    # Barrier / absorbed-power ratio
    # ------------------------------------------------------------------
    Sfam, omega_drive = sp.symbols('Sfam omega_drive', positive=True)
    U_disp = sp.simplify(-sp.Rational(1, 2) * re_expected * Sfam**2)
    P_abs = sp.simplify(sp.Rational(1, 2) * omega_drive * im_expected * Sfam**2)
    ratio_barrier = sp.simplify(P_abs / (omega_drive * (-U_disp)))
    assert sp.simplify(ratio_barrier - 1 / r) == 0

    print('\nBarrier / absorbed-power ratio:')
    print('U_disp =', U_disp)
    print('P_abs  =', P_abs)
    print('P_abs / (omega |U_disp|) =', ratio_barrier)

    # ------------------------------------------------------------------
    # Quality-factor detuning bound
    # ------------------------------------------------------------------
    omega_star, Qfac = sp.symbols('omega_star Qfac', positive=True)
    detuning_fraction = sp.simplify((r * gamma) / omega_star)
    detuning_q_form = sp.simplify(
        detuning_fraction.subs(gamma, omega_star / (2 * Qfac))
    )
    assert sp.simplify(detuning_q_form - r / (2 * Qfac)) == 0

    low_loss_detuning = sp.simplify(detuning_q_form.subs(r, 1 / eta))
    assert sp.simplify(low_loss_detuning - 1 / (2 * Qfac * eta)) == 0

    print('\nQuality-factor detuning relation:')
    print('|omega - omega_*| / omega_* =', detuning_q_form)
    print('Low-loss boundary r = 1/eta gives', low_loss_detuning)

    # ------------------------------------------------------------------
    # Linear survival window
    # ------------------------------------------------------------------
    DeltaV_req = sp.symbols('DeltaV_req', positive=True)
    Udisp_lowloss_max = sp.simplify(
        sp.Rational(1, 2) * Aabs / gamma * eta / (1 + eta**2) * Sfam**2
    )
    residue_requirement = sp.simplify(
        2 * DeltaV_req / Sfam**2 * (1 + eta**2) / eta
    )

    print('\nLinear survival-window formulas:')
    print('|U_disp|_max in the low-loss window =', Udisp_lowloss_max)
    print('Necessary residue/linewidth ratio >=', residue_requirement)

    # ------------------------------------------------------------------
    # Constructive numerical slice
    # ------------------------------------------------------------------
    # Probe-only numerical slice for sign/scale sanity checks. These values are
    # not part of the symbolic proof path and should not be cited as derived
    # constants.
    probe_values = {
        Aabs: sp.Integer(7),
        gamma: sp.Integer(2),
        r: sp.Integer(3),
        eta: sp.Rational(1, 4),
        Sfam: sp.Integer(5),
        omega_drive: sp.Integer(11),
        Qfac: sp.Integer(40),
        omega_star: sp.Integer(80),
        DeltaV_req: sp.Rational(3, 5),
    }

    print('\nConstructive numerical slice:')
    print('Re chi            =', sp.N(re_expected.subs(probe_values)))
    print('Im chi            =', sp.N(im_expected.subs(probe_values)))
    print('|Re|/|Im|         =', sp.N((re_expected / im_expected).subs(probe_values)))
    print('|U_disp|          =', sp.N((-U_disp).subs(probe_values)))
    print('P_abs             =', sp.N(P_abs.subs(probe_values)))
    print('P_abs/(omega|U|)  =', sp.N(ratio_barrier.subs(probe_values)))
    print('low-loss |U|_max  =', sp.N(Udisp_lowloss_max.subs(probe_values)))
    print('required |A|/gamma=', sp.N(residue_requirement.subs(probe_values)))


if __name__ == '__main__':
    main()
