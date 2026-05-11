import math
from typing import Sequence

import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def select_stable_real_root(poly: sp.Expr, var: sp.Symbol) -> float:
    roots = sp.nroots(poly)
    stable: list[float] = []
    for root in roots:
        real_part = float(sp.re(root))
        imag_part = float(sp.im(root))
        if abs(imag_part) < 1e-12 and 0.0 < real_part < 1.0:
            stable.append(real_part)
    if len(stable) != 1:
        raise AssertionError(f"Expected one stable real root in (0,1), found {stable}")
    return stable[0]


def main() -> None:
    sp.init_printing()
    print("=== Stage 229 SymPy audit: selected-branch numerator/denominator signature ===")

    # ------------------------------------------------------------------
    # 1. Exact selected-branch product and dimensionless factorization
    # ------------------------------------------------------------------
    x, A, DeltaK_ax, beta0 = sp.symbols("x A DeltaK_ax beta0", positive=True)
    xi, delta = sp.symbols("xi delta", positive=True)

    kappa0_sq = sp.Rational(8, 1) / sp.pi**2
    kappa1_sq = sp.Rational(16, 1) / (9 * sp.pi**2)

    s_minus = sp.simplify(
        (kappa0_sq * (x + DeltaK_ax) + kappa1_sq * x) ** 2
        / (kappa0_sq * (x + DeltaK_ax) ** 2 + kappa1_sq * x**2)
    )
    N_minus = sp.simplify(beta0 * s_minus**2 / (kappa0_sq * (A - x)))

    F = (9 * delta + 11 * xi) ** 4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
    reduced = sp.simplify(
        N_minus.subs({x: A * xi, DeltaK_ax: A * delta})
        - (8 * beta0 / (sp.pi**2 * A)) * F
    )
    assert reduced == 0

    F_num = (9 * delta + 11 * xi) ** 4 / (81 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
    F_den = 1 / (1 - xi)
    assert sp.simplify(F - F_num * F_den) == 0

    print("Verified exact selected-branch reduction:")
    print("s_-(x) =", s_minus)
    print("N_-(x) =", sp.factor(N_minus))
    print("F(xi,delta) =", F)

    # ------------------------------------------------------------------
    # 2. Exact log-slope classifier and onset / softening limits
    # ------------------------------------------------------------------
    L_num = sp.simplify(sp.diff(sp.log(F_num), xi))
    L_den = sp.simplify(sp.diff(sp.log(F_den), xi))
    R_ND = sp.simplify(L_num / L_den)

    expected_L_num = 72 * delta**2 / ((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2))
    expected_L_den = 1 / (1 - xi)
    expected_R_ND = 72 * delta**2 * (1 - xi) / ((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2))

    assert sp.simplify(L_num - expected_L_num) == 0
    assert sp.simplify(L_den - expected_L_den) == 0
    assert sp.simplify(R_ND - expected_R_ND) == 0

    onset = sp.simplify(R_ND.subs(xi, 0))
    assert sp.simplify(onset - 8 / (9 * delta)) == 0

    soft_limit = sp.limit(R_ND, xi, 1, dir='-')
    assert soft_limit == 0

    L_num_soft = sp.simplify(sp.limit(L_num, xi, 1, dir='-'))
    expected_L_num_soft = 72 * delta**2 / ((9 * delta + 11) * (9 * delta**2 + 18 * delta + 11))
    assert sp.simplify(L_num_soft - expected_L_num_soft) == 0

    print("\nVerified exact classifier data:")
    print("L_num =", L_num)
    print("L_den =", L_den)
    print("R_ND  =", R_ND)
    print("R_ND(0,delta) =", onset)
    print("lim_{xi->1^-} R_ND =", soft_limit)

    # ------------------------------------------------------------------
    # 3. Exact crossover cubic and monotonicity
    # ------------------------------------------------------------------
    numerator = sp.expand(-sp.together(R_ND - 1).as_numer_denom()[0])
    expected_P = 121 * xi**3 + 297 * delta * xi**2 + 333 * delta**2 * xi + 81 * delta**3 - 72 * delta**2
    assert sp.simplify(numerator - expected_P) == 0

    dP = sp.expand(sp.diff(expected_P, xi))
    expected_dP = 363 * xi**2 + 594 * delta * xi + 333 * delta**2
    assert sp.simplify(dP - expected_dP) == 0

    P0 = sp.simplify(expected_P.subs(xi, 0))
    assert sp.simplify(P0 - 9 * delta**2 * (9 * delta - 8)) == 0

    print("\nVerified exact crossover theorem data:")
    print("P(xi,delta)   =", expected_P)
    print("dP/dxi        =", expected_dP)
    print("P(0,delta)    =", P0)
    print("delta = 8/9 threshold gives onset balance exactly.")

    # ------------------------------------------------------------------
    # 4. Sample crossover depths on the mixed regime
    # ------------------------------------------------------------------
    xi_var = sp.symbols("xi_var", real=True)
    samples: Sequence[tuple[sp.Rational, float]] = (
        (sp.Rational(1, 4), 0.107223051105697),
        (sp.Rational(1, 2), 0.081847937860074),
        (sp.Rational(3, 4), 0.032505121082825),
    )

    for delta_value, expected_root in samples:
        poly = sp.expand(expected_P.subs({delta: delta_value, xi: xi_var}))
        root = select_stable_real_root(poly, xi_var)
        assert_close(root, expected_root, tol=5e-13)
        # Check numerator-like before the root and denominator-like after the root.
        left_probe = max(1e-6, 0.5 * root)
        right_probe = min(1 - 1e-6, 0.5 * (root + 1.0))
        left_val = float(sp.N(R_ND.subs({xi: left_probe, delta: delta_value}), 30))
        right_val = float(sp.N(R_ND.subs({xi: right_probe, delta: delta_value}), 30))
        if not left_val > 1.0:
            raise AssertionError(f"Expected numerator-like left of crossover for delta={delta_value}, got {left_val}")
        if not right_val < 1.0:
            raise AssertionError(f"Expected denominator-like right of crossover for delta={delta_value}, got {right_val}")
        print(f"delta = {sp.sstr(delta_value):>3}  ->  xi_* = {root:.15f}")

    # ------------------------------------------------------------------
    # 5. Representative always-denominator check above the threshold
    # ------------------------------------------------------------------
    delta_test = sp.Rational(1, 1)
    for probe in (sp.Rational(1, 100), sp.Rational(1, 5), sp.Rational(3, 5), sp.Rational(9, 10)):
        value = sp.N(R_ND.subs({xi: probe, delta: delta_test}), 30)
        if not float(value) < 1.0:
            raise AssertionError(f"Expected denominator-like value for delta=1 at xi={probe}, got {value}")
    print("\nVerified a representative always-denominator slice for delta = 1.")

    print("\nAll Stage 229 symbolic and numerical audits passed.")


if __name__ == "__main__":
    main()
