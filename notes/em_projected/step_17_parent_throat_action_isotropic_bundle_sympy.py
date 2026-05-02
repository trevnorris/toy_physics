#!/usr/bin/env python3
"""Parent throat action: isotropic full-bundle bridge audit.

This script verifies the symbolic bridge from the promoted parent wall block
(KSigma, MSigma) to the isotropic grouped-P2 bundle moments and target surface.
"""
from __future__ import annotations

import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def main() -> None:
    KSigma, MSigma = sp.symbols('KSigma MSigma', nonzero=True)
    B0, B2, B4, Z0, Z2, Z4 = sp.symbols('B0 B2 B4 Z0 Z2 Z4', nonzero=True)
    N0, N2, N4 = sp.symbols('N0 N2 N4', nonzero=True)
    mhat0, G, cs, a, c = sp.symbols('mhat0 G cs a c', nonzero=True)

    D0 = KSigma - B0 - Z0
    D2 = -(MSigma + B2 + Z2)
    D4 = -(B4 + Z4)

    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2
    one_pole_defect = sp.factor(sp.together(u4 - 4 * u2**2))
    one_pole_numerator = D0 * (B4 + Z4) - 3 * (MSigma + B2 + Z2)**2
    assert_zero('one-pole numerator equivalence', one_pole_defect - one_pole_numerator / D0**2)

    P0 = N0 / D0
    P2 = sp.factor(sp.together((D0 * N2 - 2 * D2 * N0) / D0**2))
    P4 = sp.factor(
        sp.together((D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3)
    )

    P0_target = sp.Rational(54) * G * cs**5 / (sp.Rational(5) * a**5 * c**5 * mhat0**2)

    K_from_one_pole = sp.simplify(sp.solve(sp.Eq(one_pole_defect, 0), KSigma)[0])
    K_from_norm = sp.simplify(sp.solve(sp.Eq(P0, P0_target), KSigma)[0])
    assert_zero('K from one-pole closed form', K_from_one_pole - (B0 + Z0 + 3*(MSigma + B2 + Z2)**2/(B4 + Z4)))
    assert_zero('K from normalization closed form', K_from_norm - (B0 + Z0 + N0/P0_target))
    compatibility = sp.factor(
        sp.together((K_from_one_pole - (B0 + Z0)) - (K_from_norm - (B0 + Z0)))
    )
    assert_zero('compatibility equation', compatibility - (3*(MSigma + B2 + Z2)**2/(B4 + Z4) - N0/P0_target))

    N2_const = sp.simplify(sp.solve(sp.Eq(P2, 0), N2)[0])
    N2_const_closed = 2 * N0 * (B2 + MSigma + Z2) / (B0 - KSigma + Z0)
    N4_const = sp.simplify(sp.solve(sp.Eq(P4.subs({N2: N2_const}), 0), N4)[0])
    N4_const_closed = N0 * (
        2 * B0 * B4
        + 2 * B0 * Z4
        + B2**2
        + 2 * B2 * MSigma
        + 2 * B2 * Z2
        - 2 * B4 * KSigma
        + 2 * B4 * Z0
        - 2 * KSigma * Z4
        + MSigma**2
        + 2 * MSigma * Z2
        + 2 * Z0 * Z4
        + Z2**2
    ) / (B0 - KSigma + Z0) ** 2
    P2_zero_eq = sp.expand(D0**2 * P2)
    P4_zero_eq = sp.expand(D0**3 * P4)
    const_prefactor_matrix = sp.Matrix([
        [sp.diff(P2_zero_eq, N2), sp.diff(P2_zero_eq, N4)],
        [sp.diff(P4_zero_eq, N2), sp.diff(P4_zero_eq, N4)],
    ])
    const_prefactor_solution = sp.solve([sp.Eq(P2_zero_eq, 0), sp.Eq(P4_zero_eq, 0)], [N2, N4], dict=True)[0]
    N4_one_pole = sp.simplify(
        N4_const.subs({KSigma: K_from_one_pole})
    )
    assert_zero('N2 constant-prefactor closed form', N2_const - N2_const_closed)
    assert_zero('N4 constant-prefactor closed form', N4_const - N4_const_closed)
    assert_zero('constant-prefactor solve determinant', const_prefactor_matrix.det() - D0**3)
    assert_zero('constant-prefactor P2 factorization', P2 - (N2 - N2_const_closed) / D0)
    assert_zero('constant-prefactor P4 factorization', P4.subs(N2, N2_const_closed) - (N4 - N4_const_closed) / D0)
    assert_zero('constant-prefactor N2 solve', const_prefactor_solution[N2] - N2_const_closed)
    assert_zero('constant-prefactor N4 solve', const_prefactor_solution[N4] - N4_const_closed)
    N4_md_one_pole = -5 * (MSigma + B2 + Z2)**2 * N0 / (KSigma - B0 - Z0)**2
    assert_zero('N4 one-pole md form', N4_const.subs(KSigma, K_from_one_pole) - N4_md_one_pole.subs(KSigma, K_from_one_pole))

    M_roots = sp.solve(sp.Eq(one_pole_numerator, 0), MSigma)
    M_root_positive = sp.sqrt((KSigma - B0 - Z0) * (B4 + Z4) / 3) - (B2 + Z2)
    M_root_negative = -sp.sqrt((KSigma - B0 - Z0) * (B4 + Z4) / 3) - (B2 + Z2)
    root_gap = sp.sqrt((KSigma - B0 - Z0) * (B4 + Z4) / 3)
    assert_zero(
        'one-pole numerator factorization in MSigma',
        one_pole_numerator + 3 * (MSigma - M_root_positive) * (MSigma - M_root_negative),
    )
    assert_zero('M-root Vieta sum', M_root_positive + M_root_negative + 2 * (B2 + Z2))
    assert_zero(
        'M-root Vieta product',
        M_root_positive * M_root_negative - ((B2 + Z2) ** 2 - root_gap**2),
    )
    u2_root_positive = sp.simplify(u2.subs(MSigma, M_root_positive))
    u2_root_negative = sp.simplify(u2.subs(MSigma, M_root_negative))
    assert_zero('u2 on positive M root', u2_root_positive - root_gap / D0)
    assert_zero('u2 on negative M root', u2_root_negative + root_gap / D0)
    D0stable, Tstable = sp.symbols('D0stable Tstable', positive=True)
    stable_u2_positive = sp.simplify(sp.sqrt(D0stable * Tstable / 3) / D0stable)
    stable_u2_negative = sp.simplify(-sp.sqrt(D0stable * Tstable / 3) / D0stable)
    if stable_u2_positive.is_positive is not True:
        raise AssertionError(f'Positive one-pole root should give positive u2 under D0>0 and B4+Z4>0: {stable_u2_positive}')
    if stable_u2_negative.is_negative is not True:
        raise AssertionError(f'Negative one-pole root should give negative u2 under D0>0 and B4+Z4>0: {stable_u2_negative}')
    assert_nonzero('stable response discriminates the two M roots', u2_root_positive - u2_root_negative)
    if len(M_roots) != 2:
        raise AssertionError(f'Expected two algebraic M roots, got {M_roots}')

    w = sp.symbols('w', real=True)
    beta = sp.exp(-w**2 / 2)
    mu_eta = sp.Integer(1)
    T_w = sp.Integer(1)
    T_omega = sp.Rational(1, 6)
    K_eta = sp.Integer(0)
    MSigma_example = sp.integrate(mu_eta * beta**2, (w, -sp.oo, sp.oo))
    KSigma_example = sp.integrate(T_w * sp.diff(beta, w)**2 + (K_eta + 6*T_omega) * beta**2, (w, -sp.oo, sp.oo))
    assert_zero('concrete wall inertia integral', MSigma_example - sp.sqrt(sp.pi))
    assert_zero('concrete wall stiffness integral', KSigma_example - 3*sp.sqrt(sp.pi)/2)

    lines = []
    lines.append('PARENT THROAT ACTION — ISOTROPIC FULL-BUNDLE BRIDGE AUDIT')
    lines.append('')
    lines.append(f'D0 = {D0}')
    lines.append(f'D2 = {D2}')
    lines.append(f'D4 = {D4}')
    lines.append('')
    lines.append(f'u2 = {sp.simplify(u2)}')
    lines.append(f'u4 = {sp.simplify(u4)}')
    lines.append(f'u4 - 4*u2^2 = {one_pole_defect}')
    lines.append('')
    lines.append(f'K_from_one_pole = {K_from_one_pole}')
    lines.append(f'K_from_norm = {K_from_norm}')
    lines.append(f'P0_target = {P0_target}')
    lines.append(f'Compatibility_equation = {compatibility}')
    lines.append('')
    lines.append(f'P2 = {P2}')
    lines.append(f'P4 = {P4}')
    lines.append(f'N2_on_constant_prefactor_branch = {N2_const}')
    lines.append(f'N4_on_constant_prefactor_branch = {N4_const}')
    lines.append(f'N4_on_one_pole_plus_constant_prefactor = {sp.factor(sp.together(N4_one_pole))}')
    lines.append(f'N4_md_equivalent_on_one_pole = {sp.factor(sp.together(N4_md_one_pole.subs(KSigma, K_from_one_pole)))}')
    lines.append(f'u2_on_positive_root = {sp.sstr(u2_root_positive)}')
    lines.append(f'u2_on_negative_root = {sp.sstr(u2_root_negative)}')
    lines.append(f'Concrete wall-integral example: MSigma={sp.sstr(MSigma_example)}, KSigma={sp.sstr(KSigma_example)}')
    lines.append('The one-pole algebra leaves two M roots, but under the stable-pole sign convention D0>0 and B4+Z4>0')
    lines.append('only the positive root gives u2>0; the negative root gives u2<0.')
    lines.append('')
    lines.append('STATUS: PASS')
    lines.append('The promoted parent wall block enters the isotropic full bundle only through KSigma and MSigma,')
    lines.append('and the one-pole / normalization / constant-prefactor conditions reduce to exact algebraic equations.')

    text = '\n'.join(lines)
    print(text)


if __name__ == '__main__':
    main()
