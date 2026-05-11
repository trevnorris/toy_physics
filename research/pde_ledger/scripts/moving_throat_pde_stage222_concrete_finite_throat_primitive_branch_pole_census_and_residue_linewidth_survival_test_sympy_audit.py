import math
from typing import Iterable, Sequence

import sympy as sp


def positive_frequency_roots(poly_expr: sp.Expr, y: sp.Symbol) -> list[float]:
    """Return positive real frequencies omega = sqrt(y) from a polynomial in y."""
    roots = sp.nroots(sp.expand(poly_expr), n=30, maxsteps=100)
    positive_y: list[float] = []
    for root in roots:
        croot = complex(root)
        if abs(croot.imag) < 1e-20 and croot.real > 0:
            positive_y.append(croot.real)
    positive_y.sort()
    return [math.sqrt(val) for val in positive_y]


def nearest_label(omega: float, wall_uncoupled: Sequence[float], internal_uncoupled: Sequence[float]) -> str:
    dw = min(abs(omega - root) for root in wall_uncoupled)
    di = min(abs(omega - root) for root in internal_uncoupled)
    return 'wall-like' if dw < di else 'internal-like'


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
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


def main() -> None:
    sp.init_printing()

    print('=== Stage 222 SymPy audit: primitive finite-throat pole census ===')

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
    # Exact primitive one-port branch
    # ------------------------------------------------------------------
    y, omega = sp.symbols('y omega', real=True)
    lam_B, lam_U, lam_W, lam_R = sp.symbols('lam_B lam_U lam_W lam_R', positive=True)
    Omega_U, Omega_W, varpi = sp.symbols('Omega_U Omega_W varpi', positive=True)
    K, M = sp.symbols('K M', positive=True)
    a, c_s = sp.symbols('a c_s', positive=True)

    kap = sp.symbols('kap', positive=True)
    C = kap * lam_B
    G_U = lam_U
    G_W = kap * lam_W
    R = kap * lam_R

    K_B = K - M * omega**2 - C**2 / (varpi**2 - omega**2)
    A = Omega_U**2 - omega**2
    W = Omega_W**2 - omega**2
    Delta = A * W - R**2
    Q = G_U**2 * W + 2 * G_U * G_W * R + G_W**2 * A
    D = sp.together(K_B - Q / Delta)

    print('\nPrimitive branch invariants:')
    print('K_B(omega) =', K_B)
    print('Delta(omega) =', Delta)
    print('Q(omega) =', Q)
    print('D(omega) =', D)

    # ------------------------------------------------------------------
    # Exact quartic pole polynomial
    # ------------------------------------------------------------------
    F_y = sp.expand(
        (((K - M * y) * (varpi**2 - y) - C**2) * ((Omega_U**2 - y) * (Omega_W**2 - y) - R**2)
         - (varpi**2 - y) * (G_U**2 * (Omega_W**2 - y) + 2 * G_U * G_W * R + G_W**2 * (Omega_U**2 - y)))
    )

    quartic_relation = sp.cancel(
        D - F_y.subs(y, omega**2) / ((varpi**2 - omega**2) * Delta)
    )
    assert quartic_relation == 0
    assert sp.Poly(F_y, y).degree() == 4

    print('\nVerified exact quartic pole polynomial:')
    print('F(y) =', F_y)
    print('with D(omega) = F(omega^2)/((varpi^2-omega^2) Delta(omega))')

    # ------------------------------------------------------------------
    # Exact residue / linewidth cancellation
    # ------------------------------------------------------------------
    D0prime, omega_star, N_star = sp.symbols('D0prime omega_star N_star', positive=True)
    eta, x, DeltaV_req = sp.symbols('eta x DeltaV_req', positive=True)

    Gamma5 = a**5 / (27 * c_s**5)
    Aqq_star = 1 / D0prime
    gamma_star = Gamma5 * omega_star**5 * N_star / D0prime
    RQ_star = sp.simplify(Aqq_star / gamma_star)
    assert sp.simplify(RQ_star - 27 * c_s**5 / (a**5 * omega_star**5 * N_star)) == 0

    low_loss_peak = sp.simplify(
        sp.Rational(1, 2) * RQ_star * eta / (1 + eta**2) / x**6
    )
    survival_threshold = sp.simplify(
        2 * DeltaV_req * (1 + eta**2) * x**6 / eta
    )

    print('\nVerified residue / linewidth cancellation:')
    print('R_Q,* =', RQ_star)
    print('\nLow-loss survival formulas:')
    print('|Re V_Q|_max =', low_loss_peak)
    print('Required R_Q,* >=', survival_threshold)

    # ------------------------------------------------------------------
    # Explicit numerical primitive slice
    # ------------------------------------------------------------------
    numerical_slice = {
        kap: kappa_primitive,
        lam_B: sp.Rational(1, 2),
        lam_U: sp.Rational(3, 10),
        lam_W: sp.Rational(2, 5),
        lam_R: sp.Rational(1, 4),
        Omega_U: sp.Integer(1),
        Omega_W: sp.Rational(7, 5),
        varpi: sp.Integer(2),
        K: sp.Integer(3),
        M: sp.Integer(1),
        a: sp.Integer(1),
        c_s: sp.Integer(1),
    }

    Delta0 = sp.simplify(Delta.subs(omega, 0))
    Q0 = sp.simplify(Q.subs(omega, 0))
    D0 = sp.simplify(D.subs(omega, 0))
    N_omega = sp.simplify((A * G_W + R * G_U) ** 2 / Delta**2)
    N0 = sp.simplify(N_omega.subs(omega, 0))
    P0 = sp.simplify(N0 / D0)

    print('\nStatic sample-slice data:')
    print('C      =', sp.N(C.subs(numerical_slice), 18))
    print('G_U    =', sp.N(G_U.subs(numerical_slice), 18))
    print('G_W    =', sp.N(G_W.subs(numerical_slice), 18))
    print('R      =', sp.N(R.subs(numerical_slice), 18))
    print('Delta0 =', sp.N(Delta0.subs(numerical_slice), 18))
    print('D0     =', sp.N(D0.subs(numerical_slice), 18))
    print('N0     =', sp.N(N0.subs(numerical_slice), 18))
    print('P0     =', sp.N(P0.subs(numerical_slice), 18))

    assert_close(float(sp.N(Delta0.subs(numerical_slice), 30)), 1.9093394081788311)
    assert_close(float(sp.N(D0.subs(numerical_slice), 30)), 2.7635551093312736)
    assert_close(float(sp.N(N0.subs(numerical_slice), 30)), 0.05016619802495911)
    assert_close(float(sp.N(P0.subs(numerical_slice), 30)), 0.018152776420332848)

    # Uncoupled and coupled pole census.
    KB_num_y = sp.expand(((K - M * y) * (varpi**2 - y) - C**2).subs(numerical_slice))
    Delta_y = sp.expand(((Omega_U**2 - y) * (Omega_W**2 - y) - R**2).subs(numerical_slice))
    F_y_num = sp.expand(F_y.subs(numerical_slice))

    wall_uncoupled = positive_frequency_roots(KB_num_y, y)
    internal_uncoupled = positive_frequency_roots(Delta_y, y)
    coupled_poles = positive_frequency_roots(F_y_num, y)

    print('\nUncoupled wall/BdG roots:', wall_uncoupled)
    print('Uncoupled internal U/W roots:', internal_uncoupled)
    print('Full coupled pole census:', coupled_poles)

    expected_wall_uncoupled = [1.6814318259147836, 2.0427400751933362]
    expected_internal_uncoupled = [0.9746017237463136, 1.417798109771174]
    expected_coupled = [0.9382727417467537, 1.3914108765380409, 1.7204537104800286, 2.045399487836587]

    for got, expected in zip(wall_uncoupled, expected_wall_uncoupled):
        assert_close(got, expected)
    for got, expected in zip(internal_uncoupled, expected_internal_uncoupled):
        assert_close(got, expected)
    for got, expected in zip(coupled_poles, expected_coupled):
        assert_close(got, expected)

    print('\nPure-Q residue / linewidth figures:')
    RQ_expr = sp.simplify(27 * c_s**5 / (a**5 * omega**5 * N_omega))
    expected_RQ = [
        18.7069287828307,
        0.380740659074003,
        16.0250330226177,
        32.0025481088465,
    ]
    for pole, expected in zip(coupled_poles, expected_RQ):
        rq_value = float(sp.N(RQ_expr.subs(numerical_slice).subs(omega, pole), 30))
        label = nearest_label(pole, wall_uncoupled, internal_uncoupled)
        print(f'omega_* = {pole:.15f} ({label}),   R_Q,* = {rq_value:.15f}')
        assert_close(rq_value, expected)

    # ------------------------------------------------------------------
    # Illustrative barrier benchmark at x = 1
    # ------------------------------------------------------------------
    V_known = sp.Float('1.181909222592')
    epsilon = sp.Float('0.1')
    DeltaV_req_num = V_known - epsilon

    threshold_eta_01 = float(sp.N((2 * DeltaV_req_num * (1 + sp.Rational(1, 10) ** 2) / sp.Rational(1, 10)), 30))
    threshold_eta_03 = float(sp.N((2 * DeltaV_req_num * (1 + sp.Rational(3, 10) ** 2) / sp.Rational(3, 10)), 30))
    assert_close(threshold_eta_01, 21.854566296358396)
    assert_close(threshold_eta_03, 7.8618736841685335)

    print('\nIllustrative low-loss thresholds at x = 1:')
    print('eta = 0.1 -> R_Q,* >=', threshold_eta_01)
    print('eta = 0.3 -> R_Q,* >=', threshold_eta_03)

    # ------------------------------------------------------------------
    # Static / dynamic tension under a lambda_W scan
    # ------------------------------------------------------------------
    lambda_W_scan = [
        sp.Rational(1, 5),
        sp.Rational(2, 5),
        sp.Rational(3, 5),
        sp.Rational(4, 5),
        sp.Integer(1),
    ]
    scan_rows: list[tuple[float, float, float, float, float]] = []

    print('\nOutgoing-leg scan (lambda_W, P0, D0, upper wall omega_*, upper wall R_Q,*):')
    for lw in lambda_W_scan:
        scan_subs = dict(numerical_slice)
        scan_subs[lam_W] = lw

        p0_value = float(sp.N(P0.subs(scan_subs), 30))
        d0_value = float(sp.N(D0.subs(scan_subs), 30))
        coupled_scan = positive_frequency_roots(sp.expand(F_y.subs(scan_subs)), y)
        upper_pole = max(coupled_scan)
        rq_upper = float(sp.N(RQ_expr.subs(scan_subs).subs(omega, upper_pole), 30))

        scan_rows.append((float(lw), p0_value, d0_value, upper_pole, rq_upper))
        print((float(lw), p0_value, d0_value, upper_pole, rq_upper))

    expected_scan = [
        (0.2, 0.005947405317693074, 2.8272344215844973, 2.04402272302752, 145.4838586578641),
        (0.4, 0.018152776420332847, 2.7635551093312736, 2.045399487836587, 32.00254810884649),
        (0.6, 0.03800016314040507, 2.665913497209664, 2.04793277506821, 13.688535635680836),
        (0.8, 0.06717078268091271, 2.5343095852196678, 2.0519066888921147, 7.580971267465815),
        (1.0, 0.10847330811047598, 2.368743373361286, 2.0577833903550973, 4.827389255647022),
    ]
    for got_row, expected_row in zip(scan_rows, expected_scan):
        for got_value, expected_value in zip(got_row, expected_row):
            assert_close(got_value, expected_value)

    assert_monotone_increasing(row[1] for row in scan_rows)
    assert_monotone_decreasing(row[4] for row in scan_rows)

    print('\nVerified static/dynamic tension on the explicit primitive family:')
    print('P0 increases monotonically with lambda_W, while the upper-wall R_Q,* decreases monotonically.')

    print('\nAll Stage 222 symbolic and numerical checks passed.')


if __name__ == '__main__':
    main()
