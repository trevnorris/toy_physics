#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 246
SymPy audit for the compensated multimode source compiler beyond positive Family-1.

This script verifies:
1. the exact mean-preserving compensated two-mode source family,
2. the exact candidate minima and the piecewise sign-change law,
3. the exact mouth-bias and shell-loading functionals,
4. the exact invertibility of the two-moment source map,
5. the exact mixed-to-shell loading ratio compiler,
6. the exact quarter-ratio test on the compensated Family-1 branch,
7. the transported radial source closure,
8. the exact sign-change threshold on the Session-I orientation,
9. and the Session-I readback of g[sigma], R[sigma], and sigma_min.
"""

from __future__ import annotations

import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Mean-preserving compensated two-mode source family.
    # ------------------------------------------------------------------
    section("1. Mean-preserving compensated two-mode source family")
    x = sp.symbols("x", real=True)
    a, b = sp.symbols("a b", real=True)
    y = sp.symbols("y", real=True)

    sigma = 1 + a * sp.cos(sp.pi * x) + b * sp.cos(2 * sp.pi * x)
    sigma_mean = sp.simplify(sp.integrate(sigma, (x, 0, 1)))
    sigma_y = sp.expand(1 - b + a * y + 2 * b * y**2)
    y_star = sp.simplify(-a / (4 * b))
    sigma_vertex = sp.simplify(sigma_y.subs(y, y_star))
    sigma_bdy_p = sp.simplify(sigma.subs(x, 0))
    sigma_bdy_m = sp.simplify(sigma.subs(x, 1))

    print("sigma(x)              =", sigma)
    print("Integral over [0,1]   =", sigma_mean)
    print("Quadratic rewrite     =", sigma_y)
    print("Stationary y*         =", y_star)
    print("Vertex value          =", sigma_vertex)
    print("Boundary values       =", sigma_bdy_p, ",", sigma_bdy_m)

    assert sigma_mean == 1
    assert sp.simplify(sigma_y - (1 - b + a * y + 2 * b * y**2)) == 0
    # Non-tautological: the actual substitution y = cos(pi*x) must turn sigma into the quadratic.
    assert sp.simplify(sigma - (1 - b + a * sp.cos(sp.pi * x) + 2 * b * sp.cos(sp.pi * x) ** 2)) == 0
    assert sp.simplify(sigma_vertex - (1 - b - a**2 / (8 * b))) == 0
    assert sp.simplify(sigma.subs({a: 0, b: 0}) - 1) == 0

    # ------------------------------------------------------------------
    # 2. Exact sign-change law.
    # ------------------------------------------------------------------
    section("2. Exact sign-change law")
    sigma_min_piece = sp.Piecewise(
        (1 + b - sp.Abs(a), sp.Or(sp.Le(b, 0), sp.Ge(sp.Abs(a), 4 * b))),
        (1 - b - a**2 / (8 * b), True),
    )

    print("sigma_min(a,b)        =")
    sp.pprint(sigma_min_piece)

    # Check three representative regions.
    test_points = [
        (sp.Rational(1, 2), sp.Rational(-1, 5)),  # b<0 boundary case
        (sp.Rational(5, 2), sp.Rational(1, 4)),   # b>0, boundary case |a|>4b
        (sp.Rational(1, 2), sp.Rational(1, 2)),   # b>0, interior minimum case
    ]
    for aval, bval in test_points:
        if (bval <= 0) or (abs(aval) >= 4 * bval):
            sigma_min_expected = sp.simplify(1 + bval - abs(aval))
        else:
            sigma_min_expected = sp.simplify(1 - bval - aval**2 / (8 * bval))
        sigma_min_test = sp.simplify(sigma_min_piece.subs({a: aval, b: bval}))
        print(f"Test point (a,b)=({aval},{bval}) -> sigma_min =", sigma_min_test)
        # True minimum of the quadratic sigma_y on y in [-1,1], computed independently
        # of the piecewise branch logic (boundary candidates always; vertex only when
        # the parabola opens upward, 2b>0, and the vertex lies in [-1,1]).
        cand = [sp.simplify(sigma_y.subs({a: aval, b: bval, y: 1})),
                sp.simplify(sigma_y.subs({a: aval, b: bval, y: -1}))]
        if bval > 0:
            ystar_val = sp.Rational(-aval, 4 * bval)
            if -1 <= ystar_val <= 1:
                cand.append(sp.simplify(sigma_y.subs({a: aval, b: bval, y: ystar_val})))
        sigma_min_true = sp.simplify(sp.Min(*cand))
        assert sp.simplify(sigma_min_test - sigma_min_true) == 0
        assert sp.simplify(sigma_min_test - sigma_min_expected) == 0

    # ------------------------------------------------------------------
    # 3. Exact mouth-bias and shell-loading functionals.
    # ------------------------------------------------------------------
    section("3. Exact mouth-bias and shell-loading functionals")
    c = sp.cos(sp.pi * x / 2)
    K_q = sp.cosh(sp.pi * (1 - x) / 2) / sp.cosh(sp.pi / 2)

    g_sigma = sp.simplify(sp.integrate(sigma * c, (x, 0, 1)))
    S_sigma = sp.simplify(sp.integrate(sigma * K_q, (x, 0, 1)))

    g_expected = sp.simplify(2 * (1 + a / 3 - b / 15) / sp.pi)
    S_expected = sp.simplify(2 * sp.tanh(sp.pi / 2) * (1 + a / 5 + b / 17) / sp.pi)

    print("g[sigma]              =", g_sigma)
    print("S[sigma]              =", S_sigma)

    assert sp.simplify(g_sigma - g_expected) == 0
    assert sp.simplify(S_sigma - S_expected) == 0

    # ------------------------------------------------------------------
    # 4. Exact two-moment source map and inverse.
    # ------------------------------------------------------------------
    section("4. Exact two-moment source map and inverse")
    g_tilde, S_tilde = sp.symbols("g_tilde S_tilde", real=True)

    M_src = sp.Matrix([[sp.Rational(1, 3), -sp.Rational(1, 15)],
                       [sp.Rational(1, 5),  sp.Rational(1, 17)]])
    det_M_src = sp.simplify(M_src.det())
    vec_tilde = sp.Matrix([g_tilde, S_tilde])

    # Here g_tilde := pi*g/2 - 1 and S_tilde := pi*S/(2 tanh(pi/2)) - 1.
    a_inv = sp.simplify(sp.Rational(85, 42) * S_tilde + sp.Rational(25, 14) * g_tilde)
    b_inv = sp.simplify(sp.Rational(425, 42) * S_tilde - sp.Rational(85, 14) * g_tilde)
    inv_check = sp.simplify(M_src * sp.Matrix([a_inv, b_inv]) - vec_tilde)

    print("M_src                 =")
    sp.pprint(M_src)
    print("det(M_src)            =", det_M_src)
    print("Inverse formulas      =", a_inv, ",", b_inv)

    assert det_M_src == sp.Rational(14, 425)
    assert inv_check == sp.Matrix([0, 0])

    # ------------------------------------------------------------------
    # 5. Exact mixed-to-shell loading ratio.
    # ------------------------------------------------------------------
    section("5. Exact mixed-to-shell loading ratio")
    rF1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
    R_sigma = sp.simplify((g_sigma - rF1)**2 / (1 + rF1**2))
    g_sym = sp.symbols("g_sym", real=True)
    R_of_g = sp.simplify((g_sym - rF1)**2 / (1 + rF1**2))

    g_minus = sp.simplify(rF1 - sp.sqrt(1 + rF1**2) / 2)
    g_plus = sp.simplify(rF1 + sp.sqrt(1 + rF1**2) / 2)

    R_minus = sp.simplify(R_of_g.subs(g_sym, g_minus))
    R_plus = sp.simplify(R_of_g.subs(g_sym, g_plus))

    g_c = sp.symbols("g_c", real=True)
    b_from_gc = sp.solve(sp.Eq(g_sigma, g_c), b)[0]

    print("r_F1                  =", sp.N(rF1, 16))
    print("R[sigma]              =", R_sigma)
    print("g_minus(r_F1)         =", sp.N(g_minus, 16))
    print("g_plus(r_F1)          =", sp.N(g_plus, 16))
    print("R(g_minus)            =", R_minus)
    print("R(g_plus)             =", R_plus)
    print("Compensation line b(a;g_c) =", b_from_gc)

    assert sp.simplify(R_minus - sp.Rational(1, 4)) == 0
    assert sp.simplify(R_plus - sp.Rational(1, 4)) == 0
    assert sp.simplify(b_from_gc - (5 * a + 15 - 15 * sp.pi * g_c / 2)) == 0

    # ------------------------------------------------------------------
    # 6. Stationary mouth law packet.
    # ------------------------------------------------------------------
    section("6. Stationary mouth law packet")
    Sigma0 = sp.symbols("Sigma0", positive=True, real=True)
    Pi_sigma = sp.simplify(Sigma0 * (1 - R_sigma * S_sigma))
    print("Pi[sigma]             =", Pi_sigma)

    # ------------------------------------------------------------------
    # 7. Transported radial source closure.
    # ------------------------------------------------------------------
    section("7. Transported radial source closure")
    r, r_sigma = sp.symbols("r r_sigma", positive=True, real=True)
    a0, b0 = sp.symbols("a0 b0", real=True)

    s_r = sp.simplify(r_sigma**2 / (r**2 + r_sigma**2))
    a_r = sp.simplify(a0 * s_r)
    b_r = sp.simplify(b0 * s_r)

    g_r = sp.simplify(g_sigma.subs({a: a_r, b: b_r}))
    S_r = sp.simplify(S_sigma.subs({a: a_r, b: b_r}))
    R_r = sp.simplify(R_sigma.subs({a: a_r, b: b_r}))

    print("s(r)                  =", s_r)
    print("g(r)                  =", g_r)
    print("S(r)                  =", S_r)
    print("R(r)                  =", R_r)

    # ------------------------------------------------------------------
    # 8. Exact sign-change threshold on the Session-I orientation.
    # ------------------------------------------------------------------
    section("8. Exact sign-change threshold for a0>0, b0<0")
    sigma_min_transport = sp.simplify(1 + b_r - a_r)
    r_thr = sp.simplify(sp.sqrt((a0 - b0 - 1) * r_sigma**2))

    print("sigma_min(r)          =", sigma_min_transport)
    print("r_signchange          =", r_thr)

    assert sp.simplify(sigma_min_transport - (1 - (a0 - b0) * s_r)) == 0
    # Verify the Session-I orientation (a0>0, b0<0) selects the boundary-minimum branch
    # of the Section-2 piecewise, rather than asserting the boundary form by hand.
    a0p = sp.symbols("a0p", positive=True)   # stands in for a0 > 0
    b0n = sp.symbols("b0n", negative=True)   # stands in for b0 < 0
    a_r_or = a0p * s_r
    b_r_or = b0n * s_r
    sigma_min_branch = sp.piecewise_fold(
        sigma_min_piece.subs({a: a_r_or, b: b_r_or})
    )
    sigma_min_branch = sp.simplify(sigma_min_branch)
    assert sp.simplify(sigma_min_branch - (1 - (a0p - b0n) * s_r)) == 0

    # ------------------------------------------------------------------
    # 9. Session-I readback.
    # ------------------------------------------------------------------
    section("9. Session-I readback")
    subs_num = {
        a0: sp.Float("2.2"),
        b0: sp.Float("-0.6"),
        r_sigma: sp.Float("0.8"),
        r: sp.Float("1.00217028"),
    }

    s_eval = sp.N(s_r.subs(subs_num), 16)
    a_eval = sp.N(a_r.subs(subs_num), 16)
    b_eval = sp.N(b_r.subs(subs_num), 16)
    g_eval = sp.N(g_r.subs(subs_num), 16)
    S_eval = sp.N(S_r.subs(subs_num), 16)
    R_eval = sp.N(R_r.subs(subs_num), 16)
    sigma_min_eval = sp.N(sigma_min_transport.subs(subs_num), 16)
    r_thr_eval = sp.N(r_thr.subs(subs_num), 16)
    g_zero_eval = sp.N(g_r.subs({**subs_num, r: sp.Integer(0)}), 16)

    print("s(r_eval)             =", s_eval)
    print("a(r_eval)             =", a_eval)
    print("b(r_eval)             =", b_eval)
    print("g[sigma](r_eval)      =", g_eval)
    print("S[sigma](r_eval)      =", S_eval)
    print("R[sigma](r_eval)      =", R_eval)
    print("sigma_min(r_eval)     =", sigma_min_eval)
    print("r_signchange          =", r_thr_eval)
    print("g[sigma](r=0)         =", g_zero_eval)
    print("sign_change?          =", bool(float(sigma_min_eval) < 0))

    assert abs(float(g_eval) - 0.82823667) < 5e-9
    assert abs(float(R_eval) - 0.21677037) < 5e-9
    assert abs(float(sigma_min_eval) - (-0.08979545)) < 5e-9
    assert abs(float(S_eval) - 0.67584771) < 5e-9
    assert float(r_thr_eval) > float(subs_num[r])
    assert float(g_zero_eval) > 1.0

    # ------------------------------------------------------------------
    # 10. Final success banner.
    # ------------------------------------------------------------------
    section("Stage 246 audit result")
    print("All symbolic checks passed.")
    print("Verified objects:")
    print("- exact mean-preserving compensated two-mode source family,")
    print("- exact sign-change law and transported sign-change threshold,")
    print("- exact mouth-bias and shell-loading functionals,")
    print("- exact invertibility of the two-moment source compiler,")
    print("- exact mixed-to-shell loading ratio and quarter-ratio test,")
    print("- exact stationary mouth-law packet Pi = Sigma0[1 - R S],")
    print("- and an exact readback of the Session-I sampled source diagnostics.")


if __name__ == "__main__":
    main()
