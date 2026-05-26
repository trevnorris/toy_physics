#!/usr/bin/env python3
"""Parent throat action: weak-axisymmetric packet audit.

This script verifies the exact wall-slope solve of the even gates once the
parent throat action is combined with the already-derived support and mixed
slope packets.
"""
from __future__ import annotations

import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def main() -> None:
    dKSigma, dMSigma = sp.symbols('dKSigma dMSigma')
    B01, B21, B41 = sp.symbols('B01 B21 B41')
    Z01, Z21, Z41 = sp.symbols('Z01 Z21 Z41')
    N01, N0 = sp.symbols('N01 N0', nonzero=True)
    KSigma, B0, Z0 = sp.symbols('KSigma B0 Z0', nonzero=True)

    D0 = KSigma - B0 - Z0
    D01 = dKSigma - B01 - Z01
    D21 = -(dMSigma + B21 + Z21)
    D41 = -(B41 + Z41)

    K1 = sp.expand(D21 + D01 / sp.Integer(9))
    H_even = sp.expand(D41 - sp.Rational(2, 3) * D21 - D01 / sp.Integer(27))
    Xi1 = sp.expand(N01 / N0 - D01 / D0)

    coeff_matrix = sp.Matrix([
        [sp.diff(K1, dKSigma), sp.diff(K1, dMSigma)],
        [sp.diff(H_even, dKSigma), sp.diff(H_even, dMSigma)],
    ])
    assert_zero('even-gate solve determinant', coeff_matrix.det() - sp.Rational(1, 27))

    sol = sp.solve([sp.Eq(K1, 0), sp.Eq(H_even, 0)], [dKSigma, dMSigma], dict=True)[0]
    D01_comp = sp.simplify(D01.subs(sol))
    D21_comp = sp.simplify(D21.subs(sol))
    D41_comp = sp.simplify(D41.subs(sol))
    Xi1_comp = sp.simplify(Xi1.subs(sol).subs({D0: KSigma - B0 - Z0}))
    closed_dK = B01 + Z01 + 27 * (B41 + Z41)
    closed_dM = -(B21 + Z21) + 3 * (B41 + Z41)
    assert_zero('dKSigma closed form', sol[dKSigma] - (B01 + Z01 + 27*(B41 + Z41)))
    assert_zero('dMSigma closed form', sol[dMSigma] - (-(B21 + Z21) + 3*(B41 + Z41)))
    assert_zero('compensated D01', D01_comp - 27*(B41 + Z41))
    assert_zero('compensated D21', D21_comp + 3*(B41 + Z41))
    assert_zero('compensated D41', D41_comp + B41 + Z41)
    assert_zero('compensated Xi1', Xi1_comp - (N01/N0 - 27*(B41 + Z41)/(KSigma - B0 - Z0)))

    lines = []
    lines.append('PARENT THROAT ACTION — WEAK-AXISYMMETRIC PACKET AUDIT')
    lines.append('')
    lines.append(f'D01 = {D01}')
    lines.append(f'D21 = {D21}')
    lines.append(f'D41 = {D41}')
    lines.append(f'K1 = {K1}')
    lines.append(f'H_even = {H_even}')
    lines.append(f'Xi1 = {Xi1}')
    lines.append('')
    lines.append(f'dKSigma_from_even_gates = {sol[dKSigma]}')
    lines.append(f'dMSigma_from_even_gates = {sol[dMSigma]}')
    lines.append('')
    lines.append(f'D01_on_compensated_branch = {D01_comp}')
    lines.append(f'D21_on_compensated_branch = {D21_comp}')
    lines.append(f'D41_on_compensated_branch = {D41_comp}')
    lines.append(f'Check_D21_plus_D01_over_9 = {sp.simplify(D21_comp + D01_comp/sp.Integer(9))}')
    lines.append(f'Check_D41_minus_2D21_over_3_minus_D01_over_27 = {sp.simplify(D41_comp - sp.Rational(2,3)*D21_comp - D01_comp/sp.Integer(27))}')
    lines.append(f'Xi1_on_compensated_branch = {Xi1_comp}')
    lines.append('')
    lines.append('STATUS: PASS')
    lines.append('The parent wall slopes are fixed uniquely by the even gates, and the residual first-order')
    lines.append('normalization defect collapses to one scalar Xi1.')

    text = '\n'.join(lines)
    print(text)


if __name__ == '__main__':
    main()
