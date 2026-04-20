import math
from typing import Iterable, Sequence

import numpy as np
import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-10) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def vector_to_float_list(v: sp.Matrix, digits: int = 15) -> list[float]:
    return [float(sp.N(v[i], digits)) for i in range(v.rows)]


def coeff_vector(expr: sp.Expr, vars_: list[sp.Symbol]) -> sp.Matrix:
    expr = sp.expand(expr)
    return sp.Matrix([sp.simplify(expr.coeff(v)) for v in vars_])


def unit_vector(v: sp.Matrix) -> sp.Matrix:
    return sp.simplify(v / sp.sqrt((v.T * v)[0]))


def positive_frequency_roots_from_coeffs(coeffs: Sequence[complex]) -> list[float]:
    """Return positive real omega roots from quartic coefficients in y = omega^2."""
    roots = np.roots(np.asarray(coeffs, dtype=np.complex128))
    positive_y: list[float] = []
    for root in roots:
        if abs(root.imag) < 1e-10 and root.real > 0:
            positive_y.append(float(root.real))
    positive_y.sort()
    return [math.sqrt(val) for val in positive_y]


def finite_dynamic_ceilings(
    rq_values: Sequence[float],
    slopes: Sequence[float],
    threshold: float,
) -> tuple[float, float]:
    """
    Return (both_poles_ceiling, nonempty_ceiling).

    both_poles_ceiling = earliest loss among degrading poles.
    nonempty_ceiling = latest loss among degrading poles, unless one pole improves,
    in which case the nonempty ceiling is infinite at first order.
    """
    finite_thresholds: list[float] = []
    improving_present = False
    for rq, slope in zip(rq_values, slopes):
        if slope < 0:
            finite_thresholds.append(math.log(rq / threshold) / (-slope))
        else:
            improving_present = True
    if not finite_thresholds:
        return math.inf, math.inf
    both = min(finite_thresholds)
    nonempty = math.inf if improving_present else max(finite_thresholds)
    return both, nonempty


def main() -> None:
    sp.init_printing()
    print("=== Stage 211 SymPy audit: numerator/denominator split and first actual dynamic-window test ===")

    # ------------------------------------------------------------------
    # 1. Primitive one-port compiler on the explicit Stage-206 branch
    # ------------------------------------------------------------------
    eps = sp.symbols("eps", real=True)
    y, omega = sp.symbols("y omega", real=True)

    kappa = sp.symbols("kappa", positive=True)
    K, M = sp.symbols("K M", positive=True)
    lamB, varpi = sp.symbols("lamB varpi", positive=True)
    lamU, lamW, lamR = sp.symbols("lamU lamW lamR", positive=True)
    OmU, OmW = sp.symbols("OmU OmW", positive=True)
    a, c_s = sp.symbols("a c_s", positive=True)

    xK, xM = sp.symbols("xK xM", real=True)
    xLB, xV = sp.symbols("xLB xV", real=True)
    xLU, xLW, xLR = sp.symbols("xLU xLW xLR", real=True)
    xOU, xOW = sp.symbols("xOU xOW", real=True)

    mixed_vars = [xLU, xLW, xLR, xOU, xOW]

    C = kappa * lamB
    GU = lamU
    GW = kappa * lamW
    R = kappa * lamR

    Delta = OmU**2 * OmW**2 - R**2
    S2 = OmU**2 + OmW**2
    H_mix = GU**2 + GW**2
    Q = GU**2 * OmW**2 + 2 * GU * GW * R + GW**2 * OmU**2
    P = OmU**2 * GW + R * GU

    B0 = C**2 / varpi**2
    B2 = C**2 / varpi**4
    B4 = C**2 / varpi**6
    Z0 = Q / Delta
    Z2 = (Q * S2 - H_mix * Delta) / Delta**2
    Z4 = (Q * (S2**2 - Delta) - S2 * H_mix * Delta) / Delta**3
    N0 = P**2 / Delta**2

    # Dressed primitive parameters.
    Kd = K * sp.exp(eps * xK)
    Md = M * sp.exp(eps * xM)
    lamBd = lamB * sp.exp(eps * xLB)
    vard = varpi * sp.exp(eps * xV)
    lamUd = lamU * sp.exp(eps * xLU)
    lamWd = lamW * sp.exp(eps * xLW)
    lamRd = lamR * sp.exp(eps * xLR)
    OmUd = OmU * sp.exp(eps * xOU)
    OmWd = OmW * sp.exp(eps * xOW)

    Cd = kappa * lamBd
    GUd = lamUd
    GWd = kappa * lamWd
    Rd = kappa * lamRd

    Deltad = OmUd**2 * OmWd**2 - Rd**2
    S2d = OmUd**2 + OmWd**2
    H_mix_d = GUd**2 + GWd**2
    Qd = GUd**2 * OmWd**2 + 2 * GUd * GWd * Rd + GWd**2 * OmUd**2
    Pd = OmUd**2 * GWd + Rd * GUd

    B0d = Cd**2 / vard**2
    B2d = Cd**2 / vard**4
    B4d = Cd**2 / vard**6
    Z0d = Qd / Deltad
    Z2d = (Qd * S2d - H_mix_d * Deltad) / Deltad**2
    Z4d = (Qd * (S2d**2 - Deltad) - S2d * H_mix_d * Deltad) / Deltad**3
    N0d = Pd**2 / Deltad**2

    D01 = sp.simplify(sp.diff(Kd - B0d - Z0d, eps).subs(eps, 0))
    D21 = sp.simplify(sp.diff(-(Md + B2d + Z2d), eps).subs(eps, 0))
    D41 = sp.simplify(sp.diff(-(B4d + Z4d), eps).subs(eps, 0))
    N01 = sp.simplify(sp.diff(N0d, eps).subs(eps, 0))
    P01 = sp.simplify(sp.diff(Pd, eps).subs(eps, 0))
    Delta01 = sp.simplify(sp.diff(Deltad, eps).subs(eps, 0))

    sample = {
        kappa: 2 * sp.sqrt(2) / sp.pi,
        lamB: sp.Rational(1, 2),
        lamU: sp.Rational(3, 10),
        lamW: sp.Rational(2, 5),
        lamR: sp.Rational(1, 4),
        OmU: sp.Integer(1),
        OmW: sp.Rational(7, 5),
        varpi: sp.Integer(2),
        M: sp.Integer(1),
        a: sp.Integer(1),
        c_s: sp.Integer(1),
    }

    # Stage-206 compatibility wall stiffness.
    B0_s = sp.simplify(B0.subs(sample))
    B2_s = sp.simplify(B2.subs(sample))
    B4_s = sp.simplify(B4.subs(sample))
    Z0_s = sp.simplify(Z0.subs(sample))
    Z2_s = sp.simplify(Z2.subs(sample))
    Z4_s = sp.simplify(Z4.subs(sample))

    D0_compat = sp.simplify(3 * (sp.Integer(1) + B2_s + Z2_s) ** 2 / (B4_s + Z4_s))
    K_compat = sp.simplify(B0_s + Z0_s + D0_compat)

    sample_full = dict(sample)
    sample_full[K] = K_compat

    zero_nonmixed = {xK: 0, xM: 0, xLB: 0, xV: 0}

    D01_mixed = sp.expand(sp.simplify(D01.subs(sample_full))).subs(zero_nonmixed)
    D21_mixed = sp.expand(sp.simplify(D21.subs(sample_full))).subs(zero_nonmixed)
    D41_mixed = sp.expand(sp.simplify(D41.subs(sample_full))).subs(zero_nonmixed)
    N01_mixed = sp.expand(sp.simplify(N01.subs(sample_full))).subs(zero_nonmixed)
    P01_mixed = sp.expand(sp.simplify(P01.subs(sample_full))).subs(zero_nonmixed)
    Delta01_mixed = sp.expand(sp.simplify(Delta01.subs(sample_full))).subs(zero_nonmixed)

    P_s = sp.simplify(P.subs(sample_full))
    Delta_s = sp.simplify(Delta.subs(sample_full))
    N0_s = sp.simplify(N0.subs(sample_full))

    Xi_transfer = sp.simplify(N01_mixed / N0_s)
    pi1 = sp.expand(sp.simplify(P01_mixed / P_s))
    delta1 = sp.expand(sp.simplify(Delta01_mixed / Delta_s))

    # ------------------------------------------------------------------
    # 2. Exact numerator/denominator split on the pure-transfer corridor
    # ------------------------------------------------------------------
    assert sp.simplify(Xi_transfer - 2 * (pi1 - delta1)) == 0

    expected_pi = [sp.Rational(3, 19), sp.Rational(16, 19), sp.Rational(3, 19), sp.Rational(32, 19), 0]
    expected_delta = [0, 0, 50 / (25 - 98 * sp.pi**2), 196 * sp.pi**2 / (98 * sp.pi**2 - 25), 196 * sp.pi**2 / (98 * sp.pi**2 - 25)]

    pi_coeff = coeff_vector(pi1, mixed_vars)
    delta_coeff = coeff_vector(delta1, mixed_vars)

    for got, expected in zip(list(pi_coeff), expected_pi):
        assert sp.simplify(got - expected) == 0
    for got, expected in zip(list(delta_coeff), expected_delta):
        assert sp.simplify(got - expected) == 0

    print("\nVerified exact numerator/denominator split theorem:")
    print("Xi_1   =", Xi_transfer)
    print("pi_1   =", pi1)
    print("delta_1=", delta1)

    # ------------------------------------------------------------------
    # 3. Pure-transfer corridor and exact rigid subcorridor counts
    # ------------------------------------------------------------------
    M_transfer = sp.Matrix([
        [sp.simplify(D01_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D21_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D41_mixed.coeff(v)) for v in mixed_vars],
    ])
    assert M_transfer.rank() == 3
    transfer_null = M_transfer.nullspace()
    assert len(transfer_null) == 2
    T = sp.Matrix.hstack(*transfer_null)

    M_num = sp.Matrix([
        [sp.simplify(D01_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D21_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D41_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(pi1.coeff(v)) for v in mixed_vars],
    ])
    M_den = sp.Matrix([
        [sp.simplify(D01_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D21_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D41_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(delta1.coeff(v)) for v in mixed_vars],
    ])
    M_both = sp.Matrix([
        [sp.simplify(D01_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D21_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(D41_mixed.coeff(v)) for v in mixed_vars],
        [sp.simplify(pi1.coeff(v)) for v in mixed_vars],
        [sp.simplify(delta1.coeff(v)) for v in mixed_vars],
    ])

    assert M_num.rank() == 4
    assert len(M_num.nullspace()) == 1
    assert M_den.rank() == 4
    assert len(M_den.nullspace()) == 1
    assert M_both.rank() == 5
    assert len(M_both.nullspace()) == 0

    det_reduced = sp.factor(sp.simplify(sp.Matrix.vstack(pi_coeff.T * T, delta_coeff.T * T).det()))
    expected_det = sp.factor(
        196 * (200 + 147 * sp.pi**2) * (80000 + 343225 * sp.pi**2 + 43218 * sp.pi**4)
        / (475 * (8670000 + 14894275 * sp.pi**2 + 2117682 * sp.pi**4))
    )
    assert sp.simplify(det_reduced - expected_det) == 0

    print("\nVerified exact subcorridor counts:")
    print("rank(pure-transfer) =", M_transfer.rank(), " nullity =", len(transfer_null))
    print("rank(+ numerator rigidity) =", M_num.rank(), " nullity =", len(M_num.nullspace()))
    print("rank(+ denominator rigidity) =", M_den.rank(), " nullity =", len(M_den.nullspace()))
    print("rank(+ both rigidities) =", M_both.rank(), " nullity =", len(M_both.nullspace()))
    print("det[(pi_1, delta_1)|_pure transfer] =", det_reduced)

    # ------------------------------------------------------------------
    # 4. Positive-Xi_1 Euclidean unit generators on the concrete branch
    # ------------------------------------------------------------------
    def orient_positive_xi(v: sp.Matrix) -> sp.Matrix:
        uv = unit_vector(v)
        xi_val = sp.N((pi_coeff.T * uv)[0] * 0 + (sp.Matrix([sp.N(sp.expand(Xi_transfer).coeff(var), 50) for var in mixed_vars]).T * uv)[0], 40)
        return sp.simplify(-uv if xi_val < 0 else uv)

    Xi_coeff_num = sp.Matrix([sp.N(sp.expand(Xi_transfer).coeff(v), 50) for v in mixed_vars])

    v_num = orient_positive_xi(M_num.nullspace()[0])
    v_den = orient_positive_xi(M_den.nullspace()[0])

    expected_v_num = [-0.55551149, 0.31814576, -0.65766801, -0.04533730, -0.39447126]
    expected_v_den = [-0.26583993, 0.18448137, 0.94454459, 0.04984499, -0.02543112]

    for actual, expected in zip(vector_to_float_list(v_num), expected_v_num):
        assert_close(actual, expected, tol=1e-8)
    for actual, expected in zip(vector_to_float_list(v_den), expected_v_den):
        assert_close(actual, expected, tol=1e-8)

    pi_num = float(sp.N((pi_coeff.T * v_num)[0], 30))
    delta_num = float(sp.N((delta_coeff.T * v_num)[0], 30))
    xi_num = float(sp.N((Xi_coeff_num.T * v_num)[0], 30))

    pi_den = float(sp.N((pi_coeff.T * v_den)[0], 30))
    delta_den = float(sp.N((delta_coeff.T * v_den)[0], 30))
    xi_den = float(sp.N((Xi_coeff_num.T * v_den)[0], 30))

    assert_close(pi_num, 0.0, tol=1e-12)
    assert_close(delta_num, -0.868056174838381, tol=1e-10)
    assert_close(xi_num, 1.73611234967676, tol=1e-10)

    assert_close(delta_den, 0.0, tol=1e-12)
    assert_close(pi_den, 0.346466075906019, tol=1e-10)
    assert_close(xi_den, 0.692932151812037, tol=1e-10)

    assert_close(xi_num, -2 * delta_num, tol=1e-10)
    assert_close(xi_den, 2 * pi_den, tol=1e-10)

    print("\nVerified positive-Xi_1 unit rigid survivors:")
    print("v_num =", vector_to_float_list(v_num, 15))
    print("pi_1(v_num)    =", pi_num)
    print("delta_1(v_num) =", delta_num)
    print("Xi_1(v_num)    =", xi_num)
    print("v_den =", vector_to_float_list(v_den, 15))
    print("pi_1(v_den)    =", pi_den)
    print("delta_1(v_den) =", delta_den)
    print("Xi_1(v_den)    =", xi_den)

    # ------------------------------------------------------------------
    # 5. Full pole census on the Stage-206 compatibility branch
    # ------------------------------------------------------------------
    F_y = sp.expand(
        (((K - M * y) * (varpi**2 - y) - C**2) * ((OmU**2 - y) * (OmW**2 - y) - R**2)
         - (varpi**2 - y) * (GU**2 * (OmW**2 - y) + 2 * GU * GW * R + GW**2 * (OmU**2 - y)))
    )
    F_poly = sp.Poly(F_y, y)
    coeff_exprs = [sp.lambdify((kappa, lamB, lamU, lamW, lamR, OmU, OmW, varpi, M, K), coeff, 'numpy') for coeff in F_poly.all_coeffs()]

    N_omega = sp.simplify(((OmU**2 - omega**2) * GW + R * GU) ** 2 / (((OmU**2 - omega**2) * (OmW**2 - omega**2) - R**2) ** 2))
    RQ_expr = sp.simplify(27 * c_s**5 / (a**5 * omega**5 * N_omega))
    rq_func = sp.lambdify((omega, kappa, lamU, lamW, lamR, OmU, OmW), RQ_expr.subs({a: 1, c_s: 1}), 'numpy')

    base_vals = {
        kappa: float(sp.N(sample[kappa], 30)),
        lamB: float(sp.N(sample[lamB], 30)),
        lamU: float(sp.N(sample[lamU], 30)),
        lamW: float(sp.N(sample[lamW], 30)),
        lamR: float(sp.N(sample[lamR], 30)),
        OmU: float(sp.N(sample[OmU], 30)),
        OmW: float(sp.N(sample[OmW], 30)),
        varpi: float(sp.N(sample[varpi], 30)),
        M: float(sp.N(sample[M], 30)),
        K: float(sp.N(K_compat, 30)),
    }

    def poles_rq_p0(direction: Sequence[float], e: float) -> tuple[list[float], list[float], float]:
        vals = dict(base_vals)
        keys = [lamU, lamW, lamR, OmU, OmW]
        for key, drift in zip(keys, direction):
            vals[key] = vals[key] * math.exp(e * drift)

        coeffs = [
            complex(func(vals[kappa], vals[lamB], vals[lamU], vals[lamW], vals[lamR], vals[OmU], vals[OmW], vals[varpi], vals[M], vals[K]))
            for func in coeff_exprs
        ]
        poles = positive_frequency_roots_from_coeffs(coeffs)
        wall_like = poles[-2:]
        rq_vals = [float(rq_func(pole, vals[kappa], vals[lamU], vals[lamW], vals[lamR], vals[OmU], vals[OmW])) for pole in wall_like]

        C_val = vals[kappa] * vals[lamB]
        GU_val = vals[lamU]
        GW_val = vals[kappa] * vals[lamW]
        R_val = vals[kappa] * vals[lamR]
        Delta_val = vals[OmU] ** 2 * vals[OmW] ** 2 - R_val ** 2
        P_val = vals[OmU] ** 2 * GW_val + R_val * GU_val
        N0_val = P_val ** 2 / Delta_val ** 2
        Z0_val = (GU_val ** 2 * vals[OmW] ** 2 + 2 * GU_val * GW_val * R_val + GW_val ** 2 * vals[OmU] ** 2) / Delta_val
        D0_val = vals[K] - C_val ** 2 / vals[varpi] ** 2 - Z0_val
        P0_val = N0_val / D0_val
        return wall_like, rq_vals, P0_val

    wall_base, rq_base, p0_base = poles_rq_p0(vector_to_float_list(v_num), 0.0)
    expected_wall_base = [1.99753567893361, 4.94905432364313]
    expected_rq_base = [30.1999075602499, 36.1711864832695]
    for got, expected in zip(wall_base, expected_wall_base):
        assert_close(got, expected, tol=2e-11)
    for got, expected in zip(rq_base, expected_rq_base):
        assert_close(got, expected, tol=2e-10)
    assert_close(p0_base, 0.002069792318062883, tol=2e-15)

    print("\nVerified undeformed compatibility-branch wall-like poles:")
    print("omega_- =", wall_base[0], "   R_Q,- =", rq_base[0])
    print("omega_+ =", wall_base[1], "   R_Q,+ =", rq_base[1])
    print("P_0     =", p0_base)

    # ------------------------------------------------------------------
    # 6. First actual dynamic-window audit on the two rigid survivors
    # ------------------------------------------------------------------
    step = 1.0e-8

    def slope_data(v: sp.Matrix) -> dict[str, float | list[float]]:
        direction = vector_to_float_list(v, 20)
        wall_m, rq_m, p0_m = poles_rq_p0(direction, -step)
        wall_p, rq_p, p0_p = poles_rq_p0(direction, +step)
        wall_0, rq_0, p0_0 = poles_rq_p0(direction, 0.0)

        wall_slopes = [(math.log(wp) - math.log(wm)) / (2 * step) for wp, wm in zip(wall_p, wall_m)]
        rq_slopes = [(math.log(rp) - math.log(rm)) / (2 * step) for rp, rm in zip(rq_p, rq_m)]
        p0_slope = (math.log(p0_p) - math.log(p0_m)) / (2 * step)

        return {
            "wall_0": wall_0,
            "rq_0": rq_0,
            "wall_slopes": wall_slopes,
            "rq_slopes": rq_slopes,
            "p0_slope": p0_slope,
        }

    num_data = slope_data(v_num)
    den_data = slope_data(v_den)

    assert_close(float(num_data["p0_slope"]), xi_num, tol=5e-7)
    assert_close(float(den_data["p0_slope"]), xi_den, tol=5e-7)

    # Negligible pole-frequency drift.
    for slope in num_data["wall_slopes"]:  # type: ignore[index]
        assert abs(float(slope)) < 5e-5
    for slope in den_data["wall_slopes"]:  # type: ignore[index]
        assert abs(float(slope)) < 5e-5

    num_rq_slopes = [float(s) for s in num_data["rq_slopes"]]  # lower, upper
    den_rq_slopes = [float(s) for s in den_data["rq_slopes"]]

    assert_close(num_rq_slopes[1], -0.52346582, tol=1e-6)
    assert_close(num_rq_slopes[0], 0.71358484, tol=1e-6)
    assert_close(den_rq_slopes[1], -0.35245541, tol=1e-6)
    assert_close(den_rq_slopes[0], -0.23169484, tol=1e-6)

    print("\nVerified first-order dynamic slopes on the rigid subcorridors:")
    print("numerator-rigid:   dln P0 =", float(num_data["p0_slope"]))
    print("                   dln R_Q,+ =", num_rq_slopes[1], "   dln R_Q,- =", num_rq_slopes[0])
    print("denominator-rigid: dln P0 =", float(den_data["p0_slope"]))
    print("                   dln R_Q,+ =", den_rq_slopes[1], "   dln R_Q,- =", den_rq_slopes[0])

    # ------------------------------------------------------------------
    # 7. Dynamic ceilings versus transported static ceilings
    # ------------------------------------------------------------------
    threshold_eta_01 = 21.854566296358396
    num_both_dyn, num_nonempty_dyn = finite_dynamic_ceilings(rq_base, num_rq_slopes, threshold_eta_01)
    den_both_dyn, den_nonempty_dyn = finite_dynamic_ceilings(rq_base, den_rq_slopes, threshold_eta_01)

    assert_close(num_both_dyn, 0.96253269, tol=2e-7)
    assert math.isinf(num_nonempty_dyn)
    assert_close(den_both_dyn, 1.39592653, tol=5e-7)
    assert_close(den_nonempty_dyn, 1.42955095, tol=5e-7)

    budget_both = 0.367930328492646
    budget_nonempty = 0.737619063660757

    num_both_stat = budget_both / xi_num
    num_nonempty_stat = budget_nonempty / xi_num
    den_both_stat = budget_both / xi_den
    den_nonempty_stat = budget_nonempty / xi_den

    assert_close(num_both_stat, 0.21192772, tol=5e-8)
    assert_close(num_nonempty_stat, 0.42486828, tol=5e-8)
    assert_close(den_both_stat, 0.53097598, tol=5e-8)
    assert_close(den_nonempty_stat, 1.06448959, tol=5e-8)

    assert num_both_dyn > num_both_stat
    assert den_both_dyn > den_both_stat
    assert den_nonempty_dyn > den_nonempty_stat

    print("\nVerified first actual dynamic-window ceilings:")
    print("numerator-rigid dynamic ceilings   =", (num_both_dyn, num_nonempty_dyn))
    print("denominator-rigid dynamic ceilings =", (den_both_dyn, den_nonempty_dyn))
    print("numerator-rigid static ceilings    =", (num_both_stat, num_nonempty_stat))
    print("denominator-rigid static ceilings  =", (den_both_stat, den_nonempty_stat))

    print("\nStage 211 audit completed successfully.")


if __name__ == "__main__":
    main()
