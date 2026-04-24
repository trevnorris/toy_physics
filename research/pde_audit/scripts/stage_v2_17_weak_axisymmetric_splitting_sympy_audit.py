#!/usr/bin/env python3
"""
Stage V2-17 — Weak-axisymmetric grouped-P2 splitting audit.

This script verifies the exact algebra behind the weak-axisymmetric grouped-P2
signature
    (20,21,22) ~ (1, 1/2, -1),
its grouped anisotropy law b = 3a, and the transported first-order response and
outgoing-prefactor slopes:
    u2^(1) = -(D21 + u2 D01)/D0,
    Xi_1 = P1/P0 = N01/N0 - D01/D0.

It also checks the canonical hidden-even relation on the outgoing l=2 branch:
    u4^(1) = (8/9) u2^(1)
iff
    D41 = (2/3) D21 + D01/27,
when the isotropic baseline has u2 = 1/9 and u4 = 4/81.
"""

from __future__ import annotations

import sympy as sp


def simp(expr):
    return sp.factor(sp.cancel(sp.together(sp.simplify(expr))))


def coeff_eps(expr, eps):
    """Coefficient of eps in a first-order expansion."""
    return simp(sp.diff(expr, eps).subs(eps, 0))


def grouped_coords(vec):
    x20, x21, x22 = vec
    xbar = simp((x20 + 2*x21 + 2*x22) / 5)
    a = simp((2*x20 - x21 - x22) / 10)
    b = simp((x21 - x22) / 2)
    return xbar, a, b


def check(name: str, expr) -> tuple[str, bool, sp.Expr]:
    val = simp(expr)
    ok = (val == 0)
    return name, ok, val


def odd_double_factorial(n: int) -> sp.Integer:
    """Return n!! for odd n, including (-1)!! = 1."""
    if n <= 0:
        return sp.Integer(1)
    out = sp.Integer(1)
    for k in range(n, 0, -2):
        out *= k
    return out


def sphere_monomial_integral(powers: tuple[int, int, int]) -> sp.Expr:
    """Integral over S^2 of x^i y^j z^k dOmega."""
    i, j, k = powers
    if i % 2 or j % 2 or k % 2:
        return sp.Integer(0)
    a, b, c = i // 2, j // 2, k // 2
    n = a + b + c
    numerator = odd_double_factorial(2*a - 1) * odd_double_factorial(2*b - 1) * odd_double_factorial(2*c - 1)
    denominator = odd_double_factorial(2*n + 1)
    return simp(4 * sp.pi * numerator / denominator)


def sphere_poly_integral(expr, x, y, z) -> sp.Expr:
    """Integral over S^2 of a polynomial in x,y,z."""
    expr = sp.expand(expr)
    poly = sp.Poly(expr, x, y, z, domain="EX")
    total = sp.Integer(0)
    for powers, coeff in poly.terms():
        total += coeff * sphere_monomial_integral(powers)
    return simp(total)


def main() -> None:
    eps = sp.symbols('eps')
    x0, x1 = sp.symbols('x0 x1')
    D0, D2, D4 = sp.symbols('D0 D2 D4', nonzero=True)
    D01, D21, D41 = sp.symbols('D01 D21 D41')
    N0, N01 = sp.symbols('N0 N01', nonzero=True)
    P0sym, P1sym = sp.symbols('P0 P1')

    # ------------------------------------------------------------------
    # 1. Exact angular triple-overlap check in the real normalized l=2 basis.
    # Use polynomial harmonics on the unit sphere; this avoids slow trig integration.
    # ------------------------------------------------------------------
    x, y, z = sp.symbols('x y z', real=True)
    Y20 = sp.sqrt(sp.Rational(5, 16) / sp.pi) * (2*z**2 - x**2 - y**2)
    Y21c = sp.sqrt(sp.Rational(15, 4) / sp.pi) * x*z
    Y21s = sp.sqrt(sp.Rational(15, 4) / sp.pi) * y*z
    Y22c = sp.sqrt(sp.Rational(15, 16) / sp.pi) * (x**2 - y**2)
    Y22s = sp.sqrt(sp.Rational(15, 4) / sp.pi) * x*y
    Ys = [Y20, Y21c, Y21s, Y22c, Y22s]
    labels5 = ['20', '21c', '21s', '22c', '22s']

    norm_checks = []
    for i, Yi in enumerate(Ys):
        for j, Yj in enumerate(Ys):
            norm_checks.append(check(f'orthonormal_{labels5[i]}_{labels5[j]}', sphere_poly_integral(Yi*Yj, x, y, z) - (1 if i == j else 0)))

    kappa_star = sp.sqrt(5) / (7*sp.sqrt(sp.pi))
    target_diag = [1, sp.Rational(1, 2), sp.Rational(1, 2), -1, -1]
    triple_checks = []
    triple_values = []
    for i, Yi in enumerate(Ys):
        val = sphere_poly_integral(Yi * Y20 * Yi, x, y, z)
        triple_values.append(val)
        triple_checks.append(check(f'triple_Y20_{labels5[i]}', val - kappa_star*target_diag[i]))

    # ------------------------------------------------------------------
    # 2. Grouped weak-axisymmetric signature.
    # ------------------------------------------------------------------
    lam20 = sp.Integer(1)
    lam21 = sp.Rational(1, 2)
    lam22 = -sp.Integer(1)
    lambdas = [lam20, lam21, lam22]

    lam_bar, lam_a, lam_b = grouped_coords(lambdas)
    lambda_checks = [
        check('lambda_trace_zero', lam_bar),
        check('lambda_a_is_1_over_4', lam_a - sp.Rational(1, 4)),
        check('lambda_b_is_3_over_4', lam_b - sp.Rational(3, 4)),
        check('lambda_b_equals_3_lambda_a', lam_b - 3*lam_a),
    ]

    weak_vec = [x0 + eps*la*x1 for la in lambdas]
    weak_bar, weak_a, weak_b = grouped_coords(weak_vec)
    weak_checks = [
        check('weak_trace_unchanged', weak_bar - x0),
        check('weak_a_equals_eps_x1_over_4', weak_a - eps*x1/4),
        check('weak_b_equals_3_eps_x1_over_4', weak_b - 3*eps*x1/4),
        check('weak_b_equals_3a', weak_b - 3*weak_a),
    ]

    # ------------------------------------------------------------------
    # 3. Transport of conservative response u2 and outgoing prefactor P0.
    # ------------------------------------------------------------------
    D0A = [D0 + eps*la*D01 for la in lambdas]
    D2A = [D2 + eps*la*D21 for la in lambdas]
    D4A = [D4 + eps*la*D41 for la in lambdas]
    N0A = [N0 + eps*la*N01 for la in lambdas]

    u2_base = simp(-D2/D0)
    u4_base = simp((D2**2 - D0*D4)/D0**2)

    u2A = [simp(-D2A[i]/D0A[i]) for i in range(3)]
    u2_slope_A = [simp(coeff_eps(u2A[i], eps)/lambdas[i]) for i in range(3)]
    u2_slope = simp(u2_slope_A[0])
    u2_expected = simp(-(D21 + u2_base*D01)/D0)

    P0_base = simp(N0/D0)
    PA = [simp(N0A[i]/D0A[i]) for i in range(3)]
    P_slope_A = [simp(coeff_eps(PA[i], eps)/lambdas[i]) for i in range(3)]
    P_slope = simp(P_slope_A[0])
    P_expected = simp((N01*D0 - N0*D01) / D0**2)
    Xi_1 = simp(P_slope / P0_base)
    Xi_expected = simp(N01/N0 - D01/D0)

    transport_checks = [
        check('u2_slope_all_lanes_equal', u2_slope_A[1] - u2_slope_A[0]),
        check('u2_slope_22_lane_equal', u2_slope_A[2] - u2_slope_A[0]),
        check('u2_slope_formula', u2_slope - u2_expected),
        check('P_slope_all_lanes_equal', P_slope_A[1] - P_slope_A[0]),
        check('P_slope_22_lane_equal', P_slope_A[2] - P_slope_A[0]),
        check('P_slope_formula', P_slope - P_expected),
        check('Xi_formula', Xi_1 - Xi_expected),
    ]

    P_weak_vec = [P0sym + eps*la*P1sym for la in lambdas]
    Pbar, Pa, Pb = grouped_coords(P_weak_vec)
    P_group_checks = [
        check('P_group_trace_unchanged', Pbar - P0sym),
        check('P_group_a', Pa - eps*P1sym/4),
        check('P_group_b', Pb - 3*eps*P1sym/4),
        check('P_group_b_equals_3a', Pb - 3*Pa),
    ]

    # ------------------------------------------------------------------
    # 4. Hidden-even relation on the canonical outgoing l=2 isotropic branch.
    # ------------------------------------------------------------------
    u4A = [simp((D2A[i]**2 - D0A[i]*D4A[i]) / D0A[i]**2) for i in range(3)]
    u4_slope_A = [simp(coeff_eps(u4A[i], eps)/lambdas[i]) for i in range(3)]
    u4_slope = simp(u4_slope_A[0])

    canonical_subs = {D2: -D0/9, D4: -D0/27}
    u2_slope_can = simp(u2_slope.subs(canonical_subs))
    u4_slope_can = simp(u4_slope.subs(canonical_subs))
    hidden_even_residual = simp(u4_slope_can - sp.Rational(8, 9)*u2_slope_can)
    hidden_even_factor = simp(hidden_even_residual * D0)
    hidden_condition = simp(D41 - sp.Rational(2, 3)*D21 - D01/27)

    hidden_checks = [
        check('canonical_u2_base', u2_base.subs(canonical_subs) - sp.Rational(1, 9)),
        check('canonical_u4_base', u4_base.subs(canonical_subs) - sp.Rational(4, 81)),
        # hidden_even_residual = -hidden_condition / D0 in this convention.
        check('hidden_even_factor', hidden_even_factor + hidden_condition),
        check('even_preserving_u2', u2_slope_can.subs({D21: -D01/9})),
        check('even_preserving_u4', u4_slope_can.subs({D21: -D01/9, D41: -D01/27})),
    ]

    # ------------------------------------------------------------------
    # 5. Summary and output.
    # ------------------------------------------------------------------
    checks = norm_checks + triple_checks + lambda_checks + weak_checks + transport_checks + P_group_checks + hidden_checks
    failed = [item for item in checks if not item[1]]

    print('STAGE V2-17 WEAK-AXISYMMETRIC SPLITTING AUDIT')
    print('=' * 72)
    print(f'total_checks: {len(checks)}')
    print(f'passed_checks: {len(checks) - len(failed)}')
    print(f'failed_checks: {len(failed)}')
    print()

    print('Angular triple-overlap values int Y_A Y20 Y_A dOmega:')
    for lab, val, mult in zip(labels5, triple_values, target_diag):
        print(f'  {lab:>3}: {sp.sstr(val)} = ({sp.sstr(mult)}) * sqrt(5)/(7*sqrt(pi))')
    print()

    print('Grouped weak-axisymmetric signature:')
    print(f'  lambda = (1, 1/2, -1)')
    print(f'  trace(lambda) = {sp.sstr(lam_bar)}')
    print(f'  a_lambda = {sp.sstr(lam_a)}')
    print(f'  b_lambda = {sp.sstr(lam_b)}')
    print(f'  b_lambda / a_lambda = {sp.sstr(simp(lam_b/lam_a))}')
    print()

    print('Generic weak split x_A = x0 + eps lambda_A x1:')
    print(f'  xbar = {sp.sstr(weak_bar)}')
    print(f'  a_x  = {sp.sstr(weak_a)}')
    print(f'  b_x  = {sp.sstr(weak_b)}')
    print()

    print('Conservative-response transport:')
    print(f'  u2 = {sp.sstr(u2_base)}')
    print(f'  u2^(1) = {sp.sstr(u2_slope)}')
    print()

    print('Outgoing-prefactor transport:')
    print(f'  P0 = {sp.sstr(P0_base)}')
    print(f'  P1 = {sp.sstr(P_slope)}')
    print(f'  Xi_1 = P1/P0 = {sp.sstr(Xi_1)}')
    print()

    print('Canonical hidden-even branch:')
    print('  Baseline constraints: D2 = -D0/9, D4 = -D0/27')
    print(f'  u2_base = {sp.sstr(u2_base.subs(canonical_subs))}')
    print(f'  u4_base = {sp.sstr(u4_base.subs(canonical_subs))}')
    print(f'  u2^(1) = {sp.sstr(u2_slope_can)}')
    print(f'  u4^(1) = {sp.sstr(u4_slope_can)}')
    print(f'  u4^(1) - (8/9)u2^(1) = {sp.sstr(hidden_even_residual)}')
    print(f'  hidden-even condition = {sp.sstr(hidden_condition)}')
    print('  so hidden-even is equivalent to D41 = (2/3)D21 + D01/27')
    print()

    print('Even-preserving compensated branch:')
    print('  u2^(1)=0  <=>  D21 = -D01/9')
    print('  hidden-even then requires D41 = -D01/27')
    print(f'  remaining prefactor defect Xi_1 = {sp.sstr(Xi_1)}')
    print()

    if failed:
        print('FAILED CHECKS:')
        for name, ok, val in failed:
            print(f'  {name}: residual = {sp.sstr(val)}')
        raise SystemExit(1)

    print('FINAL STATUS: PASS')
    print('Interpretation: weak axisymmetric l=2 splitting is one-dimensional, obeys b=3a,')
    print('and leaves one transported outgoing-prefactor scalar Xi_1=P1/P0 after even compensation.')


if __name__ == '__main__':
    main()
