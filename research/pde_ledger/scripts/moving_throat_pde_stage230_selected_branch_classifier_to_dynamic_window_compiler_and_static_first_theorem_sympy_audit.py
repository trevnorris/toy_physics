import math
from typing import Sequence

import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def dynamic_ceilings(
    s_plus: float,
    s_minus: float,
    ell_plus: float,
    ell_minus: float,
) -> tuple[float, float]:
    """Return (both_poles_ceiling, nonempty_ceiling) in |eps Xi_1|."""
    finite_thresholds: list[float] = []
    improving_present = False

    if s_plus < 0:
        finite_thresholds.append(ell_plus / (-s_plus))
    else:
        improving_present = True

    if s_minus < 0:
        finite_thresholds.append(ell_minus / (-s_minus))
    else:
        improving_present = True

    if not finite_thresholds:
        return math.inf, math.inf

    both = min(finite_thresholds)
    nonempty = math.inf if improving_present else max(finite_thresholds)
    return both, nonempty


def main() -> None:
    sp.init_printing()
    print("=== Stage 230 SymPy audit: selected-branch classifier -> dynamic window ===")

    # ------------------------------------------------------------------
    # 1. Carry the exact Stage 229 classifier and its monotonicity
    # ------------------------------------------------------------------
    xi, delta = sp.symbols("xi delta", positive=True, real=True)
    R = sp.symbols("R", nonnegative=True, real=True)

    R_ND = sp.simplify(
        72 * delta**2 * (1 - xi)
        / ((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2))
    )
    dR_dxi = sp.simplify(sp.diff(R_ND, xi))
    expected_dR_dxi = -72 * delta**2 * (
        81 * delta**3
        + 261 * delta**2
        + 297 * delta * xi * (2 - xi)
        + xi**2 * (363 - 242 * xi)
    ) / ((9 * delta + 11 * xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
    assert sp.simplify(dR_dxi - expected_dR_dxi) == 0

    onset = sp.simplify(R_ND.subs(xi, 0))
    assert sp.simplify(onset - 8 / (9 * delta)) == 0
    soft_limit = sp.limit(R_ND, xi, 1, dir="-")
    assert soft_limit == 0

    print("Verified Stage 229 selected-branch classifier:")
    print("R_ND(xi,delta) =", R_ND)
    print("dR_ND/dxi      =", dR_dxi)
    print("R_ND(0,delta)  =", onset)
    print("lim_{xi->1^-} R_ND =", soft_limit)

    # ------------------------------------------------------------------
    # 2. Carry the exact Stage 228 per-unit-Xi rigid dynamic slopes
    # ------------------------------------------------------------------
    s_plus_den = sp.Rational("-0.301516097158113")
    s_minus_den = sp.Rational("0.411024574532864")
    s_plus_num = sp.Rational("-0.508643465308977")
    s_minus_num = sp.Rational("-0.334368725711457")

    if not (s_plus_den < 0 and s_plus_num < 0 and s_minus_num < 0 < s_minus_den):
        raise AssertionError("Rigid dynamic slope signs do not match the Stage 228 carry-forward data.")

    w_num = sp.simplify(R / (1 + R))
    w_den = sp.simplify(1 / (1 + R))
    S_plus = sp.simplify(w_num * s_plus_num + w_den * s_plus_den)
    S_minus = sp.simplify(w_num * s_minus_num + w_den * s_minus_den)

    expected_S_plus = sp.simplify((R * s_plus_num + s_plus_den) / (1 + R))
    expected_S_minus = sp.simplify((R * s_minus_num + s_minus_den) / (1 + R))
    assert sp.simplify(S_plus - expected_S_plus) == 0
    assert sp.simplify(S_minus - expected_S_minus) == 0

    dS_plus = sp.simplify(sp.diff(S_plus, R))
    dS_minus = sp.simplify(sp.diff(S_minus, R))
    assert dS_plus < 0
    assert dS_minus < 0

    print("\nVerified rigid-split share compiler:")
    print("w_num(R) =", w_num)
    print("w_den(R) =", w_den)
    print("S_+(R)   =", S_plus)
    print("S_-(R)   =", S_minus)
    print("dS_+/dR  =", dS_plus)
    print("dS_-/dR  =", dS_minus)

    # ------------------------------------------------------------------
    # 3. Exact sign-flip threshold and onset threshold
    # ------------------------------------------------------------------
    R_star = sp.simplify(s_minus_den / (-s_minus_num))
    assert_close(float(sp.N(R_star, 40)), 1.229255438463336, tol=1e-15)

    delta_dyn_star = sp.simplify(8 / (9 * R_star))
    assert_close(float(sp.N(delta_dyn_star, 40)), 0.723111617875019, tol=1e-15)

    # Exact relations:
    assert sp.simplify(S_minus.subs(R, R_star)) == 0
    delta_solutions = sp.solve(sp.Eq(onset, R_star), delta)
    assert len(delta_solutions) == 1
    assert sp.simplify(delta_solutions[0] - delta_dyn_star) == 0

    # Representative sign checks.
    if not float(sp.N(S_plus.subs(R, sp.Rational(0)), 30)) < 0:
        raise AssertionError("Expected S_+(0) < 0.")
    if not float(sp.N(S_minus.subs(R, sp.Rational(1, 2)), 30)) > 0:
        raise AssertionError("Expected lower wall-like pole to improve below R_*.")
    if not float(sp.N(S_minus.subs(R, sp.Rational(2, 1)), 30)) < 0:
        raise AssertionError("Expected lower wall-like pole to worsen above R_*.")

    # Denominator-like whole-branch check for a representative slice delta = 1.
    for probe in (sp.Rational(1, 100), sp.Rational(1, 5), sp.Rational(3, 5), sp.Rational(9, 10)):
        value = float(sp.N(R_ND.subs({delta: 1, xi: probe}), 40))
        if not value < float(sp.N(R_star, 30)):
            raise AssertionError(f"Expected R_ND < R_* on delta=1 slice, got {value} at xi={probe}")

    print("\nVerified sign-flip thresholds:")
    print("R_*          =", sp.N(R_star, 30))
    print("delta_dyn_*  =", sp.N(delta_dyn_star, 30))

    # ------------------------------------------------------------------
    # 4. Dynamic window formulas in |eps Xi_1|
    # ------------------------------------------------------------------
    RQ_minus = sp.Float("30.199907560250075", 50)
    RQ_plus = sp.Float("36.171186483269487", 50)
    RQ_req = sp.Float("21.854566296358396", 50)

    ell_minus = sp.log(RQ_minus / RQ_req)
    ell_plus = sp.log(RQ_plus / RQ_req)

    assert_close(float(sp.N(ell_minus, 40)), 0.323428979934714, tol=2e-15)
    assert_close(float(sp.N(ell_plus, 40)), 0.503852964869151, tol=2e-15)

    B_plus_expr = sp.simplify(ell_plus / (-S_plus))
    B_minus_expr = sp.simplify(ell_minus / (-S_minus))

    B_both_0, B_nonempty_0 = dynamic_ceilings(
        float(sp.N(S_plus.subs(R, 0), 40)),
        float(sp.N(S_minus.subs(R, 0), 40)),
        float(sp.N(ell_plus, 40)),
        float(sp.N(ell_minus, 40)),
    )
    assert_close(B_both_0, 1.671064893775584, tol=2e-15)
    assert math.isinf(B_nonempty_0)

    B_plus_inf = float(sp.N(sp.limit(B_plus_expr, R, sp.oo), 40))
    B_minus_inf = float(sp.N(sp.limit(B_minus_expr, R, sp.oo), 40))
    B_both_inf = min(B_plus_inf, B_minus_inf)
    B_nonempty_inf = max(B_plus_inf, B_minus_inf)

    assert_close(B_both_inf, 0.967282389363822, tol=2e-15)
    assert_close(B_nonempty_inf, 0.990581810705233, tol=2e-15)

    sample_points: Sequence[tuple[sp.Expr, float, float, float, float | None]] = (
        (sp.Integer(0), -0.301516097158113, 0.411024574532864, 1.671064893775584, None),
        (sp.Integer(1), -0.405079781233545, 0.0383279244107035, 1.243836370541187, None),
        (R_star, -0.415730215182002, 0.0, 1.211971000588856, None),
        (sp.Integer(10), -0.489813704567990, -0.266605698416519, 1.028662448947899, 1.213136035184892),
    )

    print("\nSample classifier points:")
    for sample_R, expected_Sp, expected_Sm, expected_both, expected_nonempty in sample_points:
        s_plus_val = float(sp.N(S_plus.subs(R, sample_R), 40))
        s_minus_val = float(sp.N(S_minus.subs(R, sample_R), 40))
        both, nonempty = dynamic_ceilings(
            s_plus_val,
            s_minus_val,
            float(sp.N(ell_plus, 40)),
            float(sp.N(ell_minus, 40)),
        )
        assert_close(s_plus_val, expected_Sp, tol=5e-13)
        assert_close(s_minus_val, expected_Sm, tol=5e-13)
        assert_close(both, expected_both, tol=5e-13)
        if expected_nonempty is None:
            if not math.isinf(nonempty):
                raise AssertionError(f"Expected infinite nonempty ceiling at R={sample_R}, got {nonempty}")
        else:
            assert_close(nonempty, expected_nonempty, tol=5e-13)
        print(
            f"R_ND = {str(sp.N(sample_R, 16)):>18}  "
            f"S_+ = {s_plus_val:+.15f}  S_- = {s_minus_val:+.15f}  "
            f"B_both = {both:.15f}  "
            + ("B_nonempty = +inf" if math.isinf(nonempty) else f"B_nonempty = {nonempty:.15f}")
        )

    # ------------------------------------------------------------------
    # 5. Universal transported static ceilings and the static-first theorem
    # ------------------------------------------------------------------
    B_stat_both = sp.Float("0.367930328492646", 50)
    B_stat_nonempty = sp.Float("0.737619063660757", 50)

    assert B_both_inf > float(B_stat_both)
    assert B_nonempty_inf > float(B_stat_nonempty)

    print("\nVerified endpoint dynamic ceilings:")
    print("ell_- =", sp.N(ell_minus, 30))
    print("ell_+ =", sp.N(ell_plus, 30))
    print("B_dyn^(both)(0) =", B_both_0)
    print("inf B_dyn^(both) =", B_both_inf)
    print("inf B_dyn^(nonempty) over finite region =", B_nonempty_inf)

    print("\nVerified universal transported static budgets:")
    print("B_stat^(both)    =", sp.N(B_stat_both, 30))
    print("B_stat^(nonempty)=", sp.N(B_stat_nonempty, 30))

    print("\nStatic-first theorem verified:")
    print("inf B_dyn^(both)     > B_stat^(both)")
    print("inf B_dyn^(nonempty) > B_stat^(nonempty)")

    print("\nAll Stage 230 symbolic and numerical audits passed.")


if __name__ == "__main__":
    main()
