#!/usr/bin/env python3
"""
Stage V2-13: Grouped normalization ratio audit.

Symbolic audit of the bridge
  (D0,D2,D4; N0,N2,N4) -> (u2,u4; P0,P2,P4),
plus the one-pole, constant-prefactor, and universal normalization surfaces.

Implementation note: in this sandbox, running Python with normal site startup can
hang. The script is compatible with `python -S` by appending the sandbox
site-packages path before importing SymPy.
"""
from __future__ import annotations

import sys
if '/opt/pyvenv/lib/python3.13/site-packages' not in sys.path:
    sys.path.append('/opt/pyvenv/lib/python3.13/site-packages')

import sympy as sp


def record(label: str, residue: sp.Expr, results: list[tuple[str, bool, sp.Expr]]) -> None:
    """Record a polynomial-style identity check without heavy simplification."""
    r = sp.expand(residue)
    results.append((label, r == 0, r))


def main() -> None:
    x = sp.symbols('x')  # x = omega^2 for the even low-frequency expansion
    D0, D2, D4 = sp.symbols('D0 D2 D4', nonzero=True)
    N0, N2, N4 = sp.symbols('N0 N2 N4')
    results: list[tuple[str, bool, sp.Expr]] = []

    D = D0 + D2*x + D4*x**2
    N = N0 + N2*x + N4*x**2

    # Series through x^2. This is exact formal algebra.
    Y_series = sp.series(D0/D, x, 0, 3).removeO().expand()
    u2 = Y_series.coeff(x, 1)
    u4 = Y_series.coeff(x, 2)
    u2_expected = -D2/D0
    u4_expected = (D2**2 - D0*D4)/D0**2
    record('u2 = -D2/D0', (u2 - u2_expected)*D0, results)
    record('u4 = (D2^2-D0 D4)/D0^2', (u4 - u4_expected)*D0**2, results)

    Pref_series = sp.series(D0*N/D**2, x, 0, 3).removeO().expand()
    P0 = Pref_series.coeff(x, 0)
    P2 = Pref_series.coeff(x, 1)
    P4 = Pref_series.coeff(x, 2)
    P0_expected = N0/D0
    P2_expected = (D0*N2 - 2*D2*N0)/D0**2
    P4_expected = (D0**2*N4 - 2*D0*(D2*N2 + D4*N0) + 3*D2**2*N0)/D0**3
    record('P0 = N0/D0', (P0 - P0_expected)*D0, results)
    record('P2 formula', (P2 - P2_expected)*D0**2, results)
    record('P4 formula', (P4 - P4_expected)*D0**3, results)

    # Constant-prefactor branch: P2=P4=0.
    N2_const = 2*D2*N0/D0
    N4_const = N0*(D2**2 + 2*D0*D4)/D0**2
    P2_num = D0*N2 - 2*D2*N0
    P4_num = D0**2*N4 - 2*D0*(D2*N2 + D4*N0) + 3*D2**2*N0
    record('P2 numerator vanishes on constant-prefactor branch', P2_num.subs(N2, N2_const), results)
    record('P4 numerator vanishes on constant-prefactor branch', P4_num.subs({N2: N2_const, N4: N4_const}), results)

    # Outgoing l=2 fingerprint multiplication.
    Aout, Bout, G5out = sp.symbols('Aout Bout G5out')
    K0_formula = P0_expected
    K2_formula = P2_expected + Aout*P0_expected
    K4_formula = P4_expected + Aout*P2_expected + Bout*P0_expected
    Gamma5_formula = G5out*P0_expected
    record('K0 outgoing product', (K0_formula - P0_expected)*D0, results)
    record('K2 outgoing product', (K2_formula - (P2_expected + Aout*P0_expected))*D0**2, results)
    record('K4 outgoing product', (K4_formula - (P4_expected + Aout*P2_expected + Bout*P0_expected))*D0**3, results)
    record('Gamma5 outgoing product', (Gamma5_formula - G5out*P0_expected)*D0, results)

    # Isotropic full-bundle one-pole condition.
    K, M = sp.symbols('K M')
    B0, B2, B4 = sp.symbols('B0 B2 B4')
    Z0, Z2, Z4 = sp.symbols('Z0 Z2 Z4')
    D0_iso = K - B0 - Z0
    D2_iso = -(M + B2 + Z2)
    D4_iso = -(B4 + Z4)
    # D0^2*(u4-4u2^2) after substituting D2,D4.
    one_pole_computed = sp.expand((D2_iso**2 - D0_iso*D4_iso) - 4*D2_iso**2)
    one_pole_expected = sp.expand(D0_iso*(B4 + Z4) - 3*(M + B2 + Z2)**2)
    record('one-pole surface', one_pole_computed - one_pole_expected, results)

    # Universal normalization with explicit port-normalization scale.
    G, cs, ath, c_light, mhat0, S_port = sp.symbols('G c_s a c mhat0 S_port', nonzero=True)
    T_GR = 54*G*cs**5/(5*ath**5*c_light**5)
    D0_norm = mhat0**2*S_port*N0/T_GR
    D0_onepole = 3*(M + B2 + Z2)**2/(B4 + Z4)
    compat_surface = D0_norm - D0_onepole
    # If mhat0^2*S_port*P0 = T_GR, then gamma = T_GR*a^5/(27 c_s^5) = 2G/(5c^5).
    gamma_norm = T_GR*ath**5/(27*cs**5)
    gamma_GR = 2*G/(5*c_light**5)
    record('normalization implies GR gamma', (gamma_norm - gamma_GR) * (135*c_light**5/G), results)

    # Weak-axisymmetric grouped transport signature.
    eps, Pbase, P1 = sp.symbols('eps Pbase P1')
    lam20, lam21, lam22 = sp.Integer(1), sp.Rational(1, 2), -sp.Integer(1)
    P20 = Pbase + eps*lam20*P1
    P21 = Pbase + eps*lam21*P1
    P22 = Pbase + eps*lam22*P1
    Pbar = sp.expand((P20 + 2*P21 + 2*P22)/5)
    aP = sp.expand((2*P20 - P21 - P22)/10)
    bP = sp.expand((P21 - P22)/2)
    record('weak-axisymmetric grouped trace unchanged', Pbar - Pbase, results)
    record('weak-axisymmetric b=3a', bP - 3*aP, results)

    passed = sum(1 for _, ok, _ in results if ok)
    failed = len(results) - passed

    print('Stage V2-13 grouped normalization ratio audit')
    print('=' * 72)
    print(f'checks_total: {len(results)}')
    print(f'checks_passed: {passed}')
    print(f'checks_failed: {failed}')
    print()
    for label, ok, residue in results:
        print(f"[{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            print('    residue =', residue)

    print('\nCore response and prefactor formulas')
    print('-' * 72)
    print('u2 =', u2_expected)
    print('u4 =', u4_expected)
    print('P0 =', P0_expected)
    print('P2 =', P2_expected)
    print('P4 =', P4_expected)

    print('\nConstant-prefactor branch')
    print('-' * 72)
    print('N2 =', N2_const)
    print('N4 =', N4_const)

    print('\nOutgoing l=2 multiplication')
    print('-' * 72)
    print('K0 =', K0_formula)
    print('K2 =', K2_formula)
    print('K4 =', K4_formula)
    print('Gamma5 =', Gamma5_formula)

    print('\nIsotropic full-bundle target surface')
    print('-' * 72)
    print('D0 =', D0_iso)
    print('D2 =', D2_iso)
    print('D4 =', D4_iso)
    print('one_pole_surface: 0 =', one_pole_expected)
    print('T_GR =', T_GR)
    print('normalization: mhat0^2*S_port*N0/D0 = T_GR')
    print('D0_from_normalization =', D0_norm)
    print('D0_from_one_pole =', D0_onepole)
    print('compatibility_surface D0_norm - D0_onepole =', compat_surface)

    print('\nWeak-axisymmetric grouped transport')
    print('-' * 72)
    print('lambda = (1, 1/2, -1)')
    print('Pbar =', Pbar)
    print('a_P =', aP)
    print('b_P =', bP)
    print('b_P - 3 a_P =', sp.expand(bP - 3*aP))

    if failed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
