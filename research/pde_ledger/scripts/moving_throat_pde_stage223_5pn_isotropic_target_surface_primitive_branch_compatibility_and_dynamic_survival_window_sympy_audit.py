
import math
from typing import Iterable, Sequence

import numpy as np
import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-11) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f'{actual} !~= {expected} (tol={tol})')


def assert_monotone_increasing(values: Iterable[float]) -> None:
    seq = list(values)
    for left, right in zip(seq, seq[1:]):
        if not right > left:
            raise AssertionError(f'Expected a strictly increasing sequence, got {seq}')


def assert_monotone_decreasing(values: Iterable[float]) -> None:
    seq = list(values)
    for left, right in zip(seq, seq[1:]):
        if not right < left:
            raise AssertionError(f'Expected a strictly decreasing sequence, got {seq}')


def positive_frequency_roots_from_coeffs(coeffs: Sequence[complex]) -> list[float]:
    """Return positive real omega roots from quartic coefficients in y = omega^2."""
    roots = np.roots(np.asarray(coeffs, dtype=np.complex128))
    positive_y: list[float] = []
    for root in roots:
        if abs(root.imag) < 1e-10 and root.real > 0:
            positive_y.append(float(root.real))
    positive_y.sort()
    return [math.sqrt(val) for val in positive_y]


def bisect_monotone(
    func,
    left: float,
    right: float,
    target: float,
    index: int,
    iters: int = 70,
) -> tuple[float, float]:
    """
    Solve func(lw)[index] = target on a monotone interval [left, right].

    Returns
    -------
    (lambda_W_root, P0_target_compat_at_root)
    """
    f_left = func(left)[index] - target
    f_right = func(right)[index] - target
    if f_left == 0:
        return left, func(left)[0]
    if f_right == 0:
        return right, func(right)[0]
    if f_left * f_right > 0:
        raise AssertionError(
            f'Bisection bracket does not straddle the target: '
            f'left={left}, right={right}, f_left={f_left}, f_right={f_right}'
        )

    a = left
    b = right
    fa = f_left
    for _ in range(iters):
        mid = 0.5 * (a + b)
        fm = func(mid)[index] - target
        if fa * fm > 0:
            a = mid
            fa = fm
        else:
            b = mid
    root = 0.5 * (a + b)
    return root, func(root)[0]


def main() -> None:
    sp.init_printing()

    print('=== Stage 223 SymPy audit: isotropic target-surface compatibility window ===')

    # ------------------------------------------------------------------
    # Exact primitive overlap constant
    # ------------------------------------------------------------------
    s, L = sp.symbols('s L', positive=True)
    u0 = 1 / sp.sqrt(L)
    f0 = sp.sqrt(sp.Integer(2) / L) * sp.sin(sp.pi * s / (2 * L))
    kappa_primitive = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
    assert sp.simplify(kappa_primitive - 2 * sp.sqrt(2) / sp.pi) == 0

    print('\nVerified primitive finite-throat overlap constant:')
    print('kappa =', kappa_primitive)

    # ------------------------------------------------------------------
    # Exact primitive isotropic bundle
    # ------------------------------------------------------------------
    omega, y = sp.symbols('omega y', real=True)
    lam_B, lam_U, lam_W, lam_R = sp.symbols('lam_B lam_U lam_W lam_R', positive=True)
    Omega_U, Omega_W, varpi = sp.symbols('Omega_U Omega_W varpi', positive=True)
    K, M = sp.symbols('K M', positive=True)
    P0_target = sp.symbols('P0_target', positive=True)
    a, c_s = sp.symbols('a c_s', positive=True)
    kap = sp.symbols('kap', positive=True)

    C = kap * lam_B
    G_U = lam_U
    G_W = kap * lam_W
    R = kap * lam_R

    Delta = Omega_U**2 * Omega_W**2 - R**2
    S2 = Omega_U**2 + Omega_W**2
    H = G_U**2 + G_W**2
    Q = G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2
    P = Omega_U**2 * G_W + R * G_U

    B0 = C**2 / varpi**2
    B2 = C**2 / varpi**4
    B4 = C**2 / varpi**6

    Z0 = Q / Delta
    Z2 = (Q * S2 - H * Delta) / Delta**2
    Z4 = (Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3
    N0 = P**2 / Delta**2

    D0 = K - B0 - Z0
    D2 = -(M + B2 + Z2)
    D4 = -(B4 + Z4)

    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2
    P0 = N0 / D0

    print('\nPrimitive isotropic one-port data:')
    print('B0 =', B0)
    print('B2 =', B2)
    print('B4 =', B4)
    print('Z0 =', Z0)
    print('Z2 =', Z2)
    print('Z4 =', Z4)
    print('N0 =', N0)

    # ------------------------------------------------------------------
    # Exact isotropic one-pole condition and normalization compatibility
    # ------------------------------------------------------------------
    one_pole_identity = sp.cancel(
        u4 - 4 * u2**2 - (D0 * (B4 + Z4) - 3 * (M + B2 + Z2)**2) / D0**2
    )
    assert one_pole_identity == 0

    K_pole = sp.simplify(3 * (M + B2 + Z2)**2 / (B4 + Z4) + B0 + Z0)
    K_norm = sp.simplify(N0 / P0_target + B0 + Z0)
    P0_target_compat = sp.simplify(N0 * (B4 + Z4) / (3 * (M + B2 + Z2)**2))
    D0_compat = sp.simplify(K_pole - B0 - Z0)

    compatibility_identity = sp.cancel(N0 / P0_target_compat - 3 * (M + B2 + Z2)**2 / (B4 + Z4))
    assert compatibility_identity == 0
    assert sp.cancel(D0_compat - N0 / P0_target_compat) == 0

    primitive_specialization = sp.simplify(
        (P**2 / Delta**2) * (C**2 / varpi**6 + Z4) / (3 * (M + C**2 / varpi**4 + Z2)**2)
    )
    assert sp.cancel(P0_target_compat - primitive_specialization) == 0

    print('\nVerified isotropic one-pole compatibility identities:')
    print('K_pole =', K_pole)
    print('K_norm =', K_norm)
    print('P0_target_compat =', P0_target_compat)
    print('D0_compat =', D0_compat)

    # ------------------------------------------------------------------
    # Exact quartic pole polynomial
    # ------------------------------------------------------------------
    F_y = sp.expand(
        (((K - M * y) * (varpi**2 - y) - C**2) * ((Omega_U**2 - y) * (Omega_W**2 - y) - R**2)
         - (varpi**2 - y) * (G_U**2 * (Omega_W**2 - y) + 2 * G_U * G_W * R + G_W**2 * (Omega_U**2 - y)))
    )
    assert sp.Poly(F_y, y).degree() == 4

    print('\nExact quartic pole polynomial:')
    print('F(y) =', F_y)

    # ------------------------------------------------------------------
    # Numerical primitive slice from Stage 222/206
    # ------------------------------------------------------------------
    sample_slice = {
        kap: kappa_primitive,
        lam_B: sp.Rational(1, 2),
        lam_U: sp.Rational(3, 10),
        lam_W: sp.Rational(2, 5),
        lam_R: sp.Rational(1, 4),
        Omega_U: sp.Integer(1),
        Omega_W: sp.Rational(7, 5),
        varpi: sp.Integer(2),
        M: sp.Integer(1),
        a: sp.Integer(1),
        c_s: sp.Integer(1),
    }

    C_num = float(sp.N(C.subs(sample_slice), 30))
    GW_num = float(sp.N(G_W.subs(sample_slice), 30))
    R_num = float(sp.N(R.subs(sample_slice), 30))
    Delta_num = float(sp.N(Delta.subs(sample_slice), 30))
    Q_num = float(sp.N(Q.subs(sample_slice), 30))
    P_num = float(sp.N(P.subs(sample_slice), 30))
    B0_num = float(sp.N(B0.subs(sample_slice), 30))
    B2_num = float(sp.N(B2.subs(sample_slice), 30))
    B4_num = float(sp.N(B4.subs(sample_slice), 30))
    Z0_num = float(sp.N(Z0.subs(sample_slice), 30))
    Z2_num = float(sp.N(Z2.subs(sample_slice), 30))
    Z4_num = float(sp.N(Z4.subs(sample_slice), 30))
    N0_num = float(sp.N(N0.subs(sample_slice), 30))
    P0_target_compat_num = float(sp.N(P0_target_compat.subs(sample_slice), 30))
    K_compat_num = float(sp.N(K_pole.subs(sample_slice), 30))
    D0_compat_num = float(sp.N(D0_compat.subs(sample_slice), 30))

    print('\nConcrete sample compatibility branch:')
    print('C                 =', C_num)
    print('G_U               =', sample_slice[lam_U])
    print('G_W               =', GW_num)
    print('R                 =', R_num)
    print('Delta             =', Delta_num)
    print('Q                 =', Q_num)
    print('P                 =', P_num)
    print('B0                =', B0_num)
    print('B2                =', B2_num)
    print('B4                =', B4_num)
    print('Z0                =', Z0_num)
    print('Z2                =', Z2_num)
    print('Z4                =', Z4_num)
    print('N0                =', N0_num)
    print('P0_target_compat  =', P0_target_compat_num)
    print('K_compat          =', K_compat_num)
    print('D0_compat         =', D0_compat_num)

    assert_close(C_num, 0.45015815807855303)
    assert_close(GW_num, 0.36012652646284243)
    assert_close(R_num, 0.22507907903927652)
    assert_close(Delta_num, 1.9093394081788311)
    assert_close(Q_num, 0.3547252832105145)
    assert_close(P_num, 0.4276502501746254)
    assert_close(B0_num, 0.050660591821168886)
    assert_close(B2_num, 0.012665147955292221)
    assert_close(B4_num, 0.0031662869888230554)
    assert_close(Z0_num, 0.18578429884755747)
    assert_close(Z2_num, 0.1729553206266028)
    assert_close(Z4_num, 0.17082528586066765)
    assert_close(N0_num, 0.05016619802495911)
    assert_close(P0_target_compat_num, 0.0020697923180628827)
    assert_close(K_compat_num, 24.473754879290977)
    assert_close(D0_compat_num, 24.23730998862225)

    # ------------------------------------------------------------------
    # Pole census on the compatibility branch
    # ------------------------------------------------------------------
    F_y_sample = sp.expand(F_y.subs(sample_slice).subs(K, K_compat_num))
    sample_coeffs = [complex(sp.N(coeff, 30)) for coeff in sp.Poly(F_y_sample, y).all_coeffs()]
    coupled_poles = positive_frequency_roots_from_coeffs(sample_coeffs)

    expected_coupled = [
        0.971575315129468,
        1.41651290122561,
        1.99753567893361,
        4.94905432364313,
    ]
    for got, expected in zip(coupled_poles, expected_coupled):
        assert_close(got, expected, tol=2e-10)

    N_omega = sp.simplify(
        ((Omega_U**2 - omega**2) * G_W + R * G_U)**2
        / (((Omega_U**2 - omega**2) * (Omega_W**2 - omega**2) - R**2)**2)
    )
    RQ_expr = sp.simplify(27 * c_s**5 / (a**5 * omega**5 * N_omega))
    rq_values = [float(sp.N(RQ_expr.subs(sample_slice).subs(omega, pole), 30)) for pole in coupled_poles]

    expected_rq_values = [
        0.159888393135835,
        0.000806281535937178,
        30.1999075602499,
        36.1711864832695,
    ]
    for got, expected in zip(rq_values, expected_rq_values):
        assert_close(got, expected, tol=2e-10)

    print('\nCompatibility-branch pole census:')
    for pole, rq in zip(coupled_poles, rq_values):
        print(f'omega_* = {pole:.15f},   R_Q,* = {rq:.15f}')

    # ------------------------------------------------------------------
    # Barrier benchmark thresholds
    # ------------------------------------------------------------------
    V_known = sp.Float('1.181909222592')
    barrier_floor = sp.Float('0.1')
    DeltaV_req = V_known - barrier_floor

    threshold_eta_01 = float(sp.N(2 * DeltaV_req * (1 + sp.Rational(1, 10) ** 2) / sp.Rational(1, 10), 30))
    threshold_eta_03 = float(sp.N(2 * DeltaV_req * (1 + sp.Rational(3, 10) ** 2) / sp.Rational(3, 10), 30))

    assert_close(threshold_eta_01, 21.854566296358396)
    assert_close(threshold_eta_03, 7.8618736841685335)

    print('\nIllustrative compatibility-branch thresholds at x = 1:')
    print('R_Q,* required for eta = 0.1 :', threshold_eta_01)
    print('R_Q,* required for eta = 0.3 :', threshold_eta_03)

    # ------------------------------------------------------------------
    # Fast lambdified compatibility-family scan in lambda_W
    # ------------------------------------------------------------------
    lw_scan_symbol = sp.symbols('lw_scan_symbol', positive=True)
    C_scan = kappa_primitive * sp.Rational(1, 2)
    GU_scan = sp.Rational(3, 10)
    GW_scan = kappa_primitive * lw_scan_symbol
    R_scan = kappa_primitive * sp.Rational(1, 4)

    Delta_scan = sp.Integer(1) * sp.Rational(49, 25) - R_scan**2
    S2_scan = sp.Integer(1) + sp.Rational(49, 25)
    H_scan = GU_scan**2 + GW_scan**2
    Q_scan = GU_scan**2 * sp.Rational(49, 25) + 2 * GU_scan * GW_scan * R_scan + GW_scan**2
    P_scan = GW_scan + R_scan * GU_scan

    B0_scan = C_scan**2 / 4
    B2_scan = C_scan**2 / 16
    B4_scan = C_scan**2 / 64
    Z0_scan = Q_scan / Delta_scan
    Z2_scan = (Q_scan * S2_scan - H_scan * Delta_scan) / Delta_scan**2
    Z4_scan = (Q_scan * (S2_scan**2 - Delta_scan) - S2_scan * H_scan * Delta_scan) / Delta_scan**3
    N0_scan = P_scan**2 / Delta_scan**2

    P0_target_compat_scan = sp.simplify(N0_scan * (B4_scan + Z4_scan) / (3 * (1 + B2_scan + Z2_scan)**2))
    K_compat_scan = sp.simplify(3 * (1 + B2_scan + Z2_scan)**2 / (B4_scan + Z4_scan) + B0_scan + Z0_scan)

    F_y_scan = sp.expand(
        (((K_compat_scan - y) * (4 - y) - C_scan**2) * ((1 - y) * (sp.Rational(49, 25) - y) - R_scan**2)
         - (4 - y) * (GU_scan**2 * (sp.Rational(49, 25) - y) + 2 * GU_scan * GW_scan * R_scan + GW_scan**2 * (1 - y)))
    )
    coeff_funcs = [sp.lambdify(lw_scan_symbol, coeff, 'numpy') for coeff in sp.Poly(F_y_scan, y).all_coeffs()]
    pcompat_func = sp.lambdify(lw_scan_symbol, P0_target_compat_scan, 'numpy')
    kcompat_func = sp.lambdify(lw_scan_symbol, K_compat_scan, 'numpy')
    rq_func = sp.lambdify(
        (lw_scan_symbol, omega),
        sp.simplify(
            27 / (
                omega**5
                * (((1 - omega**2) * GW_scan + R_scan * GU_scan)**2 / (((1 - omega**2) * (sp.Rational(49, 25) - omega**2) - R_scan**2)**2))
            )
        ),
        'numpy',
    )

    def compatibility_family_data(lw_value: float) -> tuple[float, float, float, float]:
        coeffs_num = [complex(func(lw_value)) for func in coeff_funcs]
        roots = positive_frequency_roots_from_coeffs(coeffs_num)
        wall_like = roots[-2:]
        lower_wall, upper_wall = wall_like
        lower_rq = float(rq_func(lw_value, lower_wall))
        upper_rq = float(rq_func(lw_value, upper_wall))
        return float(pcompat_func(lw_value)), float(kcompat_func(lw_value)), lower_rq, upper_rq

    lambda_W_scan = [0.2, 0.4, 0.6, 0.8, 1.0]
    compatibility_scan = [compatibility_family_data(val) for val in lambda_W_scan]

    print('\nCompatibility-family scan:')
    print('(lambda_W, P0_target_compat, K_compat, lower wall R_Q, upper wall R_Q)')
    for lw_value, row in zip(lambda_W_scan, compatibility_scan):
        print((lw_value,) + row)

    expected_scan = [
        (0.000576970879843045, 29.3158464872314, 138.814136942081, 137.502546600713),
        (0.002069792318062883, 24.4737548792910, 30.1999075602499, 36.1711864832695),
        (0.00486568120048608, 21.1544287401845, 12.8348600273988, 16.7575510327116),
        (0.00916991368157304, 19.0298300900561, 7.06074242207991, 9.69035785242054),
        (0.0149811903240906, 17.7824591822917, 4.45922850098679, 6.30111094469551),
    ]
    for got_row, expected_row in zip(compatibility_scan, expected_scan):
        for got_value, expected_value in zip(got_row, expected_row):
            assert_close(got_value, expected_value, tol=3e-11)

    assert_monotone_increasing([row[0] for row in compatibility_scan])
    assert_monotone_decreasing([row[1] for row in compatibility_scan])
    assert_monotone_decreasing([row[2] for row in compatibility_scan])
    assert_monotone_decreasing([row[3] for row in compatibility_scan])

    # ------------------------------------------------------------------
    # Finite dynamic survival windows in the branch-compatible target
    # ------------------------------------------------------------------
    lower_strict_root, lower_strict_pcompat = bisect_monotone(
        compatibility_family_data, 0.4, 0.6, threshold_eta_01, index=2
    )
    upper_strict_root, upper_strict_pcompat = bisect_monotone(
        compatibility_family_data, 0.4, 0.6, threshold_eta_01, index=3
    )
    lower_loose_root, lower_loose_pcompat = bisect_monotone(
        compatibility_family_data, 0.6, 0.8, threshold_eta_03, index=2
    )
    upper_loose_root, upper_loose_pcompat = bisect_monotone(
        compatibility_family_data, 0.8, 1.0, threshold_eta_03, index=3
    )

    print('\nFinite dynamic survival windows in P0_target_compat:')
    print('10% loss benchmark:')
    print('  lower wall root lambda_W =', lower_strict_root)
    print('  lower wall P0_target_compat <=', lower_strict_pcompat)
    print('  upper wall root lambda_W =', upper_strict_root)
    print('  upper wall P0_target_compat <=', upper_strict_pcompat)
    print('30% loss benchmark:')
    print('  lower wall root lambda_W =', lower_loose_root)
    print('  lower wall P0_target_compat <=', lower_loose_pcompat)
    print('  upper wall root lambda_W =', upper_loose_root)
    print('  upper wall P0_target_compat <=', upper_loose_pcompat)

    assert_close(lower_strict_pcompat, 0.00283133168555932, tol=5e-13)
    assert_close(upper_strict_pcompat, 0.00359651058968466, tol=5e-13)
    assert_close(lower_loose_pcompat, 0.00817339430971383, tol=5e-13)
    assert_close(upper_loose_pcompat, 0.0116633929790174, tol=5e-13)

    print('\nStage 223 audit completed successfully.')


if __name__ == '__main__':
    main()
