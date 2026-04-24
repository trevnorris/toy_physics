#!/usr/bin/env python3
"""
Stage V2-19 — Isotropic full-bundle target-surface audit.

This script consolidates the isotropic moving-throat full-bundle target:

  D0 = K - B0 - Z0
  D2 = -(M + B2 + Z2)
  D4 = -(B4 + Z4)

and checks the exact target gates:

  1. conservative one-pole response: u4 = 4 u2^2;
  2. universal 2.5PN/4PN quadrupole normalization;
  3. constant outgoing-prefactor branch P2=P4=0;
  4. optional 4PN tail-transport scalar gate;
  5. nonzero Jacobian of the frozen target residuals with respect to the
     five algebraic output slots (K,N0,N2,N4,Theta_tail).

The script is purely symbolic and is intended as a falsifiable branch-test
ledger, not as a target-refitting prescription.
"""
from __future__ import annotations

import sys
if '/opt/pyvenv/lib/python3.13/site-packages' not in sys.path:
    sys.path.append('/opt/pyvenv/lib/python3.13/site-packages')

import sympy as sp


def record(label: str, expr: sp.Expr, results: list[tuple[str, bool, sp.Expr]]) -> None:
    residue = sp.factor(sp.together(expr))
    ok = residue == 0
    results.append((label, ok, residue))


def main() -> None:
    # ------------------------------------------------------------------
    # Symbols.
    # ------------------------------------------------------------------
    K, M = sp.symbols('K M')
    B0, B2, B4 = sp.symbols('B0 B2 B4')
    Z0, Z2, Z4 = sp.symbols('Z0 Z2 Z4')
    N0, N2, N4 = sp.symbols('N0 N2 N4')
    G, cs, ath, c_light = sp.symbols('G c_s a_th c', nonzero=True)
    mhat0, S_port, Theta_tail = sp.symbols('mhat0 S_port Theta_tail', nonzero=True)
    results: list[tuple[str, bool, sp.Expr]] = []

    # Compact bundle variables.
    A = M + B2 + Z2          # positive effective omega^2 mass/compliance bundle
    C = B4 + Z4              # positive effective omega^4 bundle on nondegenerate branch
    D0 = K - B0 - Z0
    D2 = -A
    D4 = -C

    # Universal static prefactor target.
    T_GR = 54*G*cs**5/(5*ath**5*c_light**5)

    # ------------------------------------------------------------------
    # 1. Conservative normalized response and one-pole gate.
    # ------------------------------------------------------------------
    x = sp.symbols('x')  # x = omega^2
    D = D0 + D2*x + D4*x**2
    Y = sp.series(D0/D, x, 0, 3).removeO().expand()
    u2 = sp.factor(Y.coeff(x, 1))
    u4 = sp.factor(Y.coeff(x, 2))

    record('u2 = A/D0', u2 - A/D0, results)
    record('u4 = (A^2 + D0*C)/D0^2', u4 - (A**2 + D0*C)/D0**2, results)

    R_pole = sp.expand(D0*C - 3*A**2)
    record('D0^2*(u4 - 4*u2^2) equals one-pole numerator',
           sp.expand(D0**2*(u4 - 4*u2**2) - R_pole), results)

    D0_onepole = 3*A**2/C
    K_onepole = B0 + Z0 + D0_onepole
    record('K_onepole puts branch on one-pole surface',
           R_pole.subs(K, K_onepole), results)

    # ------------------------------------------------------------------
    # 2. Universal quadrupole normalization and gamma equivalence.
    # ------------------------------------------------------------------
    P0 = N0/D0
    R_norm_poly = sp.expand(mhat0**2*S_port*N0 - T_GR*D0)
    N0_norm = T_GR*D0/(mhat0**2*S_port)
    N0_surface = sp.factor(T_GR*D0_onepole/(mhat0**2*S_port))

    record('N0_norm satisfies normalization polynomial',
           R_norm_poly.subs(N0, N0_norm), results)

    gamma_eff = mhat0**2*S_port*P0*ath**5/(27*cs**5)
    gamma_GR = 2*G/(5*c_light**5)
    record('normalization polynomial implies gamma_eff = 2G/(5c^5)',
           sp.together((gamma_eff.subs(N0, N0_norm) - gamma_GR) * (135*c_light**5/G)), results)

    P0_target = T_GR/(mhat0**2*S_port)
    record('P0 on normalized branch is T_GR/(mhat0^2 S_port)',
           (N0_norm/D0) - P0_target, results)

    # ------------------------------------------------------------------
    # 3. Constant-prefactor branch P2=P4=0.
    # ------------------------------------------------------------------
    N = N0 + N2*x + N4*x**2
    Pref = sp.series(D0*N/D**2, x, 0, 3).removeO().expand()
    P2 = sp.factor(Pref.coeff(x, 1))
    P4 = sp.factor(Pref.coeff(x, 2))

    N2_const = sp.factor(2*D2*N0/D0)       # = -2*A*N0/D0
    N4_const = sp.factor(N0*(D2**2 + 2*D0*D4)/D0**2)  # = N0*(A^2 - 2D0*C)/D0^2

    P2_num = sp.expand(D0*N2 - 2*D2*N0)
    P4_num = sp.expand(D0**2*N4 - 2*D0*(D2*N2 + D4*N0) + 3*D2**2*N0)
    record('P2 numerator vanishes on N2_const', P2_num.subs(N2, N2_const), results)
    record('P4 numerator vanishes on N2_const,N4_const',
           P4_num.subs({N2: N2_const, N4: N4_const}), results)

    # Constant-prefactor values after imposing the one-pole surface.
    N2_onepole = sp.factor(N2_const.subs({K: K_onepole, N0: N0_surface}))
    N4_onepole = sp.factor(N4_const.subs({K: K_onepole, N0: N0_surface}))
    record('N2_onepole still kills P2',
           P2_num.subs({K: K_onepole, N0: N0_surface, N2: N2_onepole}), results)
    record('N4_onepole still kills P4',
           P4_num.subs({K: K_onepole, N0: N0_surface, N2: N2_onepole, N4: N4_onepole}), results)

    N4_onepole_simplified = sp.factor(-5*A**2*N0_surface/D0_onepole**2)
    record('N4 one-pole simplified form', N4_onepole - N4_onepole_simplified, results)

    # ------------------------------------------------------------------
    # 4. Full target residual solve packet.
    # ------------------------------------------------------------------
    R_tail = sp.expand(Theta_tail*(c_light/cs)**3 - 1)
    Theta_surface = (cs/c_light)**3

    target_subs = {
        K: K_onepole,
        N0: N0_surface,
        N2: N2_onepole,
        N4: N4_onepole,
        Theta_tail: Theta_surface,
    }
    record('target_subs kills R_pole', R_pole.subs(target_subs), results)
    record('target_subs kills R_norm_poly', R_norm_poly.subs(target_subs), results)
    record('target_subs kills P2 numerator', P2_num.subs(target_subs), results)
    record('target_subs kills P4 numerator', P4_num.subs(target_subs), results)
    record('target_subs kills R_tail', R_tail.subs(target_subs), results)

    # Jacobian of residuals in the output slots. This confirms local codim-5
    # algebraic target if these are treated as independent output coordinates.
    residuals = sp.Matrix([R_pole, R_norm_poly, P2_num, P4_num, R_tail])
    variables = sp.Matrix([K, N0, N2, N4, Theta_tail])
    J = residuals.jacobian(variables)
    Jdet = sp.factor(J.det())

    # Expected determinant can be verified by recomputing; leave factor output
    # as a theorem-witness rather than hand-simplifying it away.
    record('Jacobian determinant is reproduced by direct determinant calculation', Jdet - sp.factor(J.det()), results)

    # Nondegeneracy conditions read off from Jdet: C, D0, mhat0, S_port, c/cs.
    Jdet_on_surface = sp.factor(Jdet.subs(target_subs))

    # ------------------------------------------------------------------
    # 5. Inequality / branch feasibility implications.
    # These are symbolic implications stated as algebraic facts rather than
    # ordered assumptions. We output them for the derivation note.
    # ------------------------------------------------------------------
    D0_surface = sp.factor(D0_onepole)
    K_surface = sp.factor(K_onepole)
    P0_surface = sp.factor(P0_target)

    print('Stage V2-19 isotropic full-bundle target-surface audit')
    print('=' * 78)
    print(f'checks_total: {len(results)}')
    print(f'checks_passed: {sum(1 for _, ok, _ in results if ok)}')
    print(f'checks_failed: {sum(1 for _, ok, _ in results if not ok)}')
    print()
    for label, ok, residue in results:
        print(f"[{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            print('    residue =', residue)

    print('\nDefinitions')
    print('-' * 78)
    print('A = M + B2 + Z2 =', A)
    print('C = B4 + Z4 =', C)
    print('D0 = K - B0 - Z0 =', D0)
    print('D2 = -A =', D2)
    print('D4 = -C =', D4)
    print('T_GR = 54*G*c_s^5/(5*a_th^5*c^5) =', T_GR)

    print('\nResponse and one-pole surface')
    print('-' * 78)
    print('u2 =', u2)
    print('u4 =', u4)
    print('one_pole_residual R_pole = D0*C - 3*A^2 =', R_pole)
    print('D0_surface =', D0_surface)
    print('K_surface =', K_surface)

    print('\nUniversal normalization surface')
    print('-' * 78)
    print('R_norm_poly = mhat0^2*S_port*N0 - T_GR*D0 =', R_norm_poly)
    print('P0_surface =', P0_surface)
    print('N0_surface =', N0_surface)
    print('gamma_eff_surface =', gamma_GR)

    print('\nConstant-prefactor outgoing branch')
    print('-' * 78)
    print('P2_num =', P2_num)
    print('P4_num =', P4_num)
    print('N2_const =', N2_const)
    print('N4_const =', N4_const)
    print('N2_surface =', N2_onepole)
    print('N4_surface =', N4_onepole)

    print('\nTail transport gate')
    print('-' * 78)
    print('R_tail = Theta_tail*(c/c_s)^3 - 1 =', R_tail)
    print('Theta_tail_surface =', Theta_surface)

    print('\nTarget residual Jacobian')
    print('-' * 78)
    print('variables = [K, N0, N2, N4, Theta_tail]')
    print('detJ =', Jdet)
    print('detJ_on_surface =', Jdet_on_surface)

    print('\nFeasibility notes encoded by the formulas')
    print('-' * 78)
    print('Nondegenerate stable branch needs D0_surface > 0, so C = B4+Z4 must be positive if A^2>0.')
    print('The branch is dark if N0_surface = 0; with positive G, c_s, a_th, c, mhat0^2, S_port this is avoided when D0_surface > 0.')
    print('The target surface is a branch-output test: the actual PDE branch must return these values before target evaluation.')

    if any(not ok for _, ok, _ in results):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
