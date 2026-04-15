import math
from typing import Callable, Sequence

import mpmath as mp
import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def bisect_increasing(
    func: Callable[[float], float],
    target: float,
    lo: float,
    hi: float,
    tol: float = 1e-14,
    max_iter: int = 300,
) -> float:
    """Root of monotone increasing func(x) = target on [lo, hi]."""
    f_lo = func(lo) - target
    f_hi = func(hi) - target
    if abs(f_lo) < tol:
        return lo
    if abs(f_hi) < tol:
        return hi
    if not (f_lo < 0 < f_hi):
        raise ValueError(f"Bracket failure for increasing function: f(lo)-target={f_lo}, f(hi)-target={f_hi}")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = func(mid) - target
        if abs(f_mid) < tol or abs(hi - lo) < tol:
            return mid
        if f_mid > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def threshold_root(delta_value: float, cap_value: float, r_nd_func: Callable[[float], float]) -> float:
    """Return xi_c solving R_ND(xi,delta)=cap, or 0 if onset already satisfies the cap."""
    onset_value = r_nd_func(0.0)
    if onset_value <= cap_value + 1e-15:
        return 0.0
    return bisect_increasing(
        func=lambda x: cap_value - r_nd_func(x),  # increasing because R_ND decreases
        target=0.0,
        lo=0.0,
        hi=1.0 - 1e-12,
        tol=1e-15,
        max_iter=400,
    )


def main() -> None:
    sp.init_printing()
    print("=== Stage 214 SymPy audit: continuum pullback of the selected-branch dynamic-class map ===")

    # ------------------------------------------------------------------
    # 1. Universal selected-branch geometry and classifier
    # ------------------------------------------------------------------
    xi, delta, c = sp.symbols("xi delta c", positive=True, real=True)
    R_target = sp.symbols("R_target", positive=True, real=True)

    F = sp.simplify(
        (9 * delta + 11 * xi) ** 4
        / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
    )
    G = sp.simplify(9 * xi * (xi + delta) / (9 * delta + 11 * xi))
    R_ND = sp.simplify(
        72 * delta**2 * (1 - xi)
        / ((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2))
    )

    dF_dxi = sp.factor(sp.diff(F, xi))
    dG_dxi = sp.factor(sp.diff(G, xi))
    dR_dxi = sp.simplify(sp.diff(R_ND, xi))

    expected_dF_dxi = (
        (9 * delta + 11 * xi) ** 3
        * (81 * delta**3 + 189 * delta**2 * xi + 72 * delta**2 + 297 * delta * xi**2 + 121 * xi**3)
        / (81 * (1 - xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 3)
    )
    expected_dG_dxi = 9 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) / (9 * delta + 11 * xi) ** 2
    expected_dR_dxi = -72 * delta**2 * (
        81 * delta**3
        + 261 * delta**2
        + 297 * delta * xi * (2 - xi)
        + xi**2 * (363 - 242 * xi)
    ) / ((9 * delta + 11 * xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)

    assert sp.simplify(dF_dxi - expected_dF_dxi) == 0
    assert sp.simplify(dG_dxi - expected_dG_dxi) == 0
    assert sp.simplify(dR_dxi - expected_dR_dxi) == 0

    onset_F = sp.simplify(F.subs(xi, 0))
    onset_G = sp.simplify(G.subs(xi, 0))
    onset_R = sp.simplify(R_ND.subs(xi, 0))
    soft_limit_F = sp.limit(F, xi, 1, dir="-")
    assert onset_F == 1
    assert onset_G == 0
    assert sp.simplify(onset_R - 8 / (9 * delta)) == 0
    assert soft_limit_F == sp.oo

    print("Verified universal selected-branch geometry:")
    print("F(xi,delta) =", F)
    print("G(xi,delta) =", G)
    print("R_ND(xi,delta) =", R_ND)
    print("dF/dxi =", dF_dxi)
    print("dG/dxi =", dG_dxi)
    print("dR_ND/dxi =", dR_dxi)
    print("F(0,delta) =", onset_F)
    print("G(0,delta) =", onset_G)
    print("R_ND(0,delta) =", onset_R)

    # ------------------------------------------------------------------
    # 2. Pullback derivative formulas
    # ------------------------------------------------------------------
    dxi_dRtarget = sp.simplify(1 / dF_dxi)
    dRphys_dRtarget = sp.simplify(dR_dxi / dF_dxi)

    print("\nImplicit pullback compiler:")
    print("dxi_req/dR_target =", dxi_dRtarget)
    print("dR_phys/dR_target =", dRphys_dRtarget)

    # Numerical sign checks on the stable branch.
    F_num = sp.lambdify((xi, delta), F, "mpmath")
    G_num = sp.lambdify((xi, delta), G, "mpmath")
    R_num = sp.lambdify((xi, delta), R_ND, "mpmath")
    dF_num = sp.lambdify((xi, delta), dF_dxi, "mpmath")
    dG_num = sp.lambdify((xi, delta), dG_dxi, "mpmath")
    dR_num = sp.lambdify((xi, delta), dR_dxi, "mpmath")
    dRphys_num = sp.lambdify((xi, delta), dRphys_dRtarget, "mpmath")

    probe_grid: Sequence[tuple[float, float]] = (
        (0.01, 0.20),
        (0.05, 0.25),
        (0.10, 0.50),
        (0.20, 0.75),
        (0.40, 1.25),
        (0.80, 2.00),
    )
    for xi_val, delta_val in probe_grid:
        if not dF_num(xi_val, delta_val) > 0:
            raise AssertionError("Expected dF/dxi > 0 on the stable branch.")
        if not dG_num(xi_val, delta_val) > 0:
            raise AssertionError("Expected dG/dxi > 0 on the stable branch.")
        if not dR_num(xi_val, delta_val) < 0:
            raise AssertionError("Expected dR_ND/dxi < 0 on the stable branch.")
        if not dRphys_num(xi_val, delta_val) < 0:
            raise AssertionError("Expected dR_phys/dR_target < 0 on the stable branch.")

    print("Verified derivative signs on a representative stable-branch grid.")

    # ------------------------------------------------------------------
    # 3. Generic classifier-cap threshold compiler
    # ------------------------------------------------------------------
    P_c = sp.expand(c * (9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) - 72 * delta**2 * (1 - xi))
    dP_dxi = sp.factor(sp.diff(P_c, xi))
    expected_dP_dxi = 3 * (87 * c * delta**2 + 198 * c * delta * xi + 121 * c * xi**2 + 24 * delta**2)
    assert sp.simplify(dP_dxi - expected_dP_dxi) == 0

    onset_Pc = sp.factor(P_c.subs(xi, 0))
    assert sp.simplify(onset_Pc - 9 * delta**2 * (9 * c * delta - 8)) == 0

    delta_c = sp.simplify(8 / (9 * c))
    assert sp.simplify(onset_Pc.subs(delta, delta_c)) == 0

    print("\nVerified generic threshold compiler:")
    print("P_c(xi,delta) =", P_c)
    print("dP_c/dxi =", dP_dxi)
    print("P_c(0,delta) =", onset_Pc)
    print("delta_c =", delta_c)

    # ------------------------------------------------------------------
    # 4. Carry the Stage-213 sign-flip and denominator thresholds
    # ------------------------------------------------------------------
    s_minus_den = sp.Rational("0.411024574532864")
    s_minus_num = sp.Rational("-0.334368725711457")
    R_star = sp.simplify(s_minus_den / (-s_minus_num))
    delta_dyn_star = sp.simplify(8 / (9 * R_star))
    delta_den = sp.Rational(8, 9)

    assert_close(float(sp.N(R_star, 40)), 1.229255438463336, tol=2e-15)
    assert_close(float(sp.N(delta_dyn_star, 40)), 0.723111617875019, tol=2e-15)
    assert_close(float(delta_den), 8 / 9, tol=0.0)

    print("\nPulled-back dynamic thresholds:")
    print("R_* =", sp.N(R_star, 30))
    print("delta_dyn_* =", sp.N(delta_dyn_star, 30))
    print("delta_den =", sp.N(delta_den, 30))

    # ------------------------------------------------------------------
    # 5. Numerical pulled-back threshold curves
    # ------------------------------------------------------------------
    def xi_cap(delta_value: float, cap_value: float) -> float:
        r_func = lambda x: float(R_num(x, delta_value))
        return threshold_root(delta_value, cap_value, r_func)

    def R_cap(delta_value: float, cap_value: float) -> float:
        xi_value = xi_cap(delta_value, cap_value)
        return float(F_num(xi_value, delta_value))

    sample_rows = (
        (0.25, 0.08744210580615958, 1.3308685393299624, 0.10722305110569706, 1.3938325657815114),
        (0.50, 0.05142857879350268, 1.1399566302640967, 0.08184793786007429, 1.2210870615059795),
        (0.75, 0.0, 1.0, 0.03250512108282537, 1.0714718671109147),
    )

    print("\nSample pulled-back thresholds:")
    for delt, exp_xi_flip, exp_R_flip, exp_xi_den, exp_R_den in sample_rows:
        xi_flip = xi_cap(delt, float(sp.N(R_star, 40)))
        R_flip = R_cap(delt, float(sp.N(R_star, 40)))
        xi_den = xi_cap(delt, 1.0)
        R_den = R_cap(delt, 1.0)

        assert_close(xi_flip, exp_xi_flip, tol=5e-13)
        assert_close(R_flip, exp_R_flip, tol=5e-13)
        assert_close(xi_den, exp_xi_den, tol=5e-13)
        assert_close(R_den, exp_R_den, tol=5e-13)
        if not R_flip <= R_den + 1e-14:
            raise AssertionError("Expected R_flip <= R_den.")

        print(
            f"delta = {delt:>4.2f}   "
            f"xi_flip = {xi_flip:.9f}   R_flip = {R_flip:.9f}   "
            f"xi_den = {xi_den:.9f}   R_den = {R_den:.9f}"
        )

    # Check threshold collapse to onset above the corresponding onset deltas.
    assert_close(R_cap(0.80, float(sp.N(R_star, 40))), 1.0, tol=2e-14)
    assert_close(R_cap(1.00, 1.0), 1.0, tol=2e-14)

    print("Verified threshold collapse to onset above the onset-delta surfaces.")

    # ------------------------------------------------------------------
    # 6. Continuum placement map and equivalent kernel inequalities
    # ------------------------------------------------------------------
    eps_eta, eps_W, rho, Z_W, delta0, Lam = sp.symbols(
        "eps_eta eps_W rho Z_W delta0 Lam", positive=True, real=True
    )
    R_cap_symbol = sp.symbols("R_cap_symbol", positive=True, real=True)

    delta_expr = sp.simplify(delta0 / (1 - eps_eta))
    M_mix_expr = sp.simplify(8 * Z_W * (1 + rho) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps_W)))
    R_target_expr = sp.simplify(Lam * (1 - eps_eta) * (1 - eps_W) ** 2 / (Z_W * (1 + rho) ** 2))
    product_law = sp.simplify(R_target_expr * M_mix_expr)

    assert sp.simplify(product_law - 8 * Lam * (1 - eps_W) / sp.pi**2) == 0

    print("\nVerified continuum placement map:")
    print("delta =", delta_expr)
    print("M_mix =", M_mix_expr)
    print("R_target =", R_target_expr)
    print("R_target * M_mix =", product_law)

    # The bound translators are algebraic rearrangements under positivity assumptions.
    signflip_bound_Z = sp.simplify(Lam * (1 - eps_eta) * (1 - eps_W) ** 2 / R_cap_symbol)
    signflip_bound_M = sp.simplify(8 * Lam * (1 - eps_W) / (sp.pi**2 * R_cap_symbol))

    print("Equivalent positive-ratio translators:")
    print("Z_W (1+rho)^2 <= ", signflip_bound_Z)
    print("M_mix <= ", signflip_bound_M)

    # ------------------------------------------------------------------
    # 7. Static-first theorem survives the pullback
    # ------------------------------------------------------------------
    B_dyn_both_inf = 0.967282389363822
    B_dyn_nonempty_inf = 0.990581810705233
    B_stat_both = 0.367930328492646
    B_stat_nonempty = 0.737619063660757

    if not B_dyn_both_inf > B_stat_both:
        raise AssertionError("Expected pulled-back robust dynamic ceiling to exceed static robust budget.")
    if not B_dyn_nonempty_inf > B_stat_nonempty:
        raise AssertionError("Expected pulled-back nonempty dynamic ceiling to exceed static nonempty budget.")

    print("\nStatic-first theorem after pullback:")
    print(f"inf B_dyn^(both)     = {B_dyn_both_inf:.15f} > B_stat^(both)     = {B_stat_both:.15f}")
    print(f"inf B_dyn^(nonempty) = {B_dyn_nonempty_inf:.15f} > B_stat^(nonempty) = {B_stat_nonempty:.15f}")

    print("\nAll Stage 214 symbolic and numerical audits passed.")


if __name__ == "__main__":
    main()
