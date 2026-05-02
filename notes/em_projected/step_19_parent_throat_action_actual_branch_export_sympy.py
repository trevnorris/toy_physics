#!/usr/bin/env python3
"""Checklist audit for step_19_parent_throat_action_actual_branch_export_notes.md."""
from __future__ import annotations

import hashlib
import json
import numpy as np
import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    value = sp.factor(sp.together(sp.simplify(expr)))
    if value == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def assert_matrix_zero(label: str, matrix: sp.Matrix) -> None:
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            assert_zero(f"{label}[{i},{j}]", matrix[i, j])


def assert_small_nonzero(label: str, expr: sp.Expr, upper: float) -> None:
    value = float(sp.N(sp.Abs(expr), 30))
    if not (0.0 < value < upper):
        raise AssertionError(f"{label} failed: |value|={value} not in (0, {upper})")


def assert_numeric_small(label: str, value: float, upper: float) -> None:
    if abs(value) >= upper:
        raise AssertionError(f"{label} failed: |value|={value} not below {upper}")


def assert_numeric_small_nonzero(label: str, value: float, upper: float) -> None:
    if not (0.0 < abs(value) < upper):
        raise AssertionError(f"{label} failed: |value|={value} not in (0, {upper})")


def main() -> None:
    KSigma, MSigma = sp.symbols("KSigma MSigma", nonzero=True)
    B0, B2, B4, Z0, Z2, Z4 = sp.symbols("B0 B2 B4 Z0 Z2 Z4", nonzero=True)
    N0, N2, N4 = sp.symbols("N0 N2 N4", nonzero=True)
    mhat0, G, cs, a, c = sp.symbols("mhat0 G cs a c", nonzero=True)

    D0 = KSigma - B0 - Z0
    D2 = -(MSigma + B2 + Z2)
    D4 = -(B4 + Z4)
    P0 = N0 / D0
    P2 = (D0 * N2 - 2 * D2 * N0) / D0**2
    P4 = (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3

    one_pole_residue = D0 * (B4 + Z4) - 3 * (MSigma + B2 + Z2) ** 2
    K_one_pole = B0 + Z0 + 3 * (MSigma + B2 + Z2) ** 2 / (B4 + Z4)
    assert_zero("isotropic one-pole checklist", one_pole_residue.subs(KSigma, K_one_pole))

    static_target = 54 * G * cs**5 / (5 * a**5 * c**5)
    K_norm = B0 + Z0 + mhat0**2 * N0 / static_target
    assert_zero("static normalization checklist", (mhat0**2 * P0 - static_target).subs(KSigma, K_norm))

    N2_const = sp.solve(sp.Eq(P2, 0), N2, dict=True)[0][N2]
    N2_const_closed = 2 * N0 * (B2 + MSigma + Z2) / (B0 - KSigma + Z0)
    N4_const = sp.solve(sp.Eq(P4.subs(N2, N2_const), 0), N4, dict=True)[0][N4]
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
    assert_zero("constant-prefactor N2 closed form", N2_const - N2_const_closed)
    assert_zero("constant-prefactor N4 closed form", N4_const - N4_const_closed)
    assert_zero("constant-prefactor solve determinant", const_prefactor_matrix.det() - D0**3)
    assert_zero("constant-prefactor P2 factorization", P2 - (N2 - N2_const_closed) / D0)
    assert_zero("constant-prefactor P4 factorization", P4.subs(N2, N2_const_closed) - (N4 - N4_const_closed) / D0)
    assert_zero("constant-prefactor N2 solve", const_prefactor_solution[N2] - N2_const_closed)
    assert_zero("constant-prefactor N4 solve", const_prefactor_solution[N4] - N4_const_closed)

    dKSigma, dMSigma = sp.symbols("dKSigma dMSigma")
    B01, B21, B41 = sp.symbols("B01 B21 B41")
    Z01, Z21, Z41 = sp.symbols("Z01 Z21 Z41")
    N01 = sp.symbols("N01")

    D01 = dKSigma - B01 - Z01
    D21 = -(dMSigma + B21 + Z21)
    D41 = -(B41 + Z41)
    K1 = D21 + D01 / 9
    H_even = D41 - sp.Rational(2, 3) * D21 - D01 / 27
    Xi1 = N01 / N0 - D01 / D0

    canonical = {
        dKSigma: B01 + Z01 + 27 * (B41 + Z41),
        dMSigma: -(B21 + Z21) + 3 * (B41 + Z41),
    }
    assert_zero("canonical weak K1 checklist", K1.subs(canonical))
    assert_zero("canonical weak H_even checklist", H_even.subs(canonical))
    assert_zero("canonical D01 checklist", D01.subs(canonical) - 27 * (B41 + Z41))
    assert_zero("canonical D21 checklist", D21.subs(canonical) + 3 * (B41 + Z41))
    assert_zero("canonical D41 checklist", D41 + B41 + Z41)

    N01_for_zero_Xi = sp.solve(sp.Eq(Xi1.subs(canonical), 0), N01, dict=True)[0][N01]
    assert_zero("outgoing slope needed for Xi1=0", N01_for_zero_Xi - N0 * 27 * (B41 + Z41) / D0)

    # Concrete export packet: a branch-derived wall block plus independent
    # support/mixed/outgoing toy port exports. This is intentionally not tuned
    # to the target surfaces; it checks what an honest first export does.
    w, omega = sp.symbols("w omega", real=True)
    R = sp.symbols("R", real=True)
    R0 = sp.exp(-w**2 / 2)
    R0p = sp.diff(R0, w)
    R0pp = sp.diff(R0p, w)
    beta = sp.exp(-w**2 / 2)

    mu_shape = 1 + w**2 / 4
    tw_shape = 1 + w**2 / 6
    to_shape = 1 + w**2 / 8
    mu_eta = 1 + mu_shape * R0
    T_w = 1 + tw_shape * R0
    T_omega = (1 + to_shape * R0) / 6

    # Parent action choice:
    #   muSigma(R,w) = 1 + (1 + w^2/4) R
    #   TwSigma(R,w) = 1 + (1 + w^2/6) R
    #   TOmegaSigma(R,w) = (1 + (1 + w^2/8) R)/6
    # and USigma(R,w) = 1/2 m(w) R^2 with m(w) chosen so that R0(w) is a
    # stationary isotropic background of the wall sector.
    U_scale = sp.simplify((sp.diff(T_w * R0p, w) - tw_shape * R0p**2 / 2) / R0)
    background_residue = sp.simplify(sp.diff(T_w * R0p, w) - tw_shape * R0p**2 / 2 - U_scale * R0)
    K_eta = sp.simplify(U_scale - sp.diff(tw_shape * R0p, w))
    assert_zero("concrete branch background equation", background_residue)
    assert_zero(
        "concrete branch K_eta",
        K_eta - (w**4 * sp.exp(-w**2 / 2) / 12 + w**2 + w**2 * sp.exp(-w**2 / 2) / 2 - 1),
    )

    M_export = sp.integrate(mu_eta * beta**2, (w, -sp.oo, sp.oo))
    K_export = sp.integrate(T_w * sp.diff(beta, w)**2 + (K_eta + 6*T_omega) * beta**2, (w, -sp.oo, sp.oo))
    assert_zero("branch-derived wall inertia integral", M_export - sp.sqrt(sp.pi) * (13 * sp.sqrt(6) + 36) / 36)
    assert_zero("branch-derived wall stiffness integral", K_export - sp.sqrt(sp.pi) * (24 + 13 * sp.sqrt(6)) / 24)

    # Concrete parent-background port model: export each sector from a finite
    # Galerkin resolvent of the actual l=2 parent-background wall operator
    #
    #   L_parent[f] = -d_w(T_w d_w f) + (K_eta + 6 T_omega) f
    #
    # on the explicit branch chosen above.  The basis is a three-mode even
    # Hermite Galerkin packet; B and N lie exactly in that packet, while Z is
    # exported with a measured projection residual.
    phi_B = beta
    phi_Z = R0 * beta
    phi_N = (1 + w**2 / 2) * beta

    def basis_mode(n: int) -> sp.Expr:
        return sp.simplify(
            sp.hermite(n, w) * sp.exp(-w**2 / 2) / (sp.pi ** sp.Rational(1, 4) * sp.sqrt(2**n * sp.factorial(n)))
        )

    def parent_apply(expr: sp.Expr) -> sp.Expr:
        return sp.simplify(-sp.diff(T_w * sp.diff(expr, w), w) + (K_eta + 6 * T_omega) * expr)

    modes = (0, 2, 4)
    basis = [basis_mode(n) for n in modes]
    for i, psi_i in enumerate(basis):
        assert_zero(f"basis normalization n={modes[i]}", sp.integrate(psi_i**2, (w, -sp.oo, sp.oo)) - 1)
        for j in range(i + 1, len(basis)):
            assert_zero(
                f"basis orthogonality {modes[i]}/{modes[j]}",
                sp.integrate(psi_i * basis[j], (w, -sp.oo, sp.oo)),
            )

    mass_matrix = sp.Matrix(
        [
            [
                sp.simplify(sp.integrate(mu_eta * basis[i] * basis[j], (w, -sp.oo, sp.oo)))
                for j in range(len(basis))
            ]
            for i in range(len(basis))
        ]
    )
    stiff_matrix = sp.Matrix(
        [
            [
                sp.simplify(
                    sp.integrate(
                        T_w * sp.diff(basis[i], w) * sp.diff(basis[j], w)
                        + (K_eta + 6 * T_omega) * basis[i] * basis[j],
                        (w, -sp.oo, sp.oo),
                    )
                )
                for j in range(len(basis))
            ]
            for i in range(len(basis))
        ]
    )
    assert_matrix_zero("parent mass matrix symmetry", mass_matrix - mass_matrix.T)
    assert_matrix_zero("parent stiffness matrix symmetry", stiff_matrix - stiff_matrix.T)
    assert_nonzero("parent mass matrix determinant", mass_matrix.det())
    assert_nonzero("parent stiffness matrix determinant", stiff_matrix.det())

    def project_source(phi: sp.Expr) -> tuple[sp.Matrix, sp.Expr, sp.Expr]:
        coeffs = sp.Matrix([sp.simplify(sp.integrate(basis[i] * phi, (w, -sp.oo, sp.oo))) for i in range(len(basis))])
        proj = sp.simplify(sum(coeffs[i] * basis[i] for i in range(len(basis))))
        residual = sp.simplify(sp.integrate((phi - proj) ** 2, (w, -sp.oo, sp.oo)))
        return coeffs, proj, residual

    def low_series_packet(label: str, phi: sp.Expr) -> tuple[sp.Matrix, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
        coeffs, proj, residual = project_source(phi)
        a0 = sp.Matrix([sp.simplify(x) for x in stiff_matrix.LUsolve(coeffs)])
        a2 = sp.Matrix([sp.simplify(-x) for x in stiff_matrix.LUsolve(mass_matrix * a0)])
        a4 = sp.Matrix([sp.simplify(-x) for x in stiff_matrix.LUsolve(mass_matrix * a2)])
        assert_matrix_zero(f"{label} order-omega0 matrix equation", stiff_matrix * a0 - coeffs)
        assert_matrix_zero(f"{label} order-omega2 matrix equation", stiff_matrix * a2 + mass_matrix * a0)
        assert_matrix_zero(f"{label} order-omega4 matrix equation", stiff_matrix * a4 + mass_matrix * a2)
        resp0 = sp.simplify(sum(a0[i] * basis[i] for i in range(len(basis))))
        resp2 = sp.simplify(sum(a2[i] * basis[i] for i in range(len(basis))))
        resp4 = sp.simplify(sum(a4[i] * basis[i] for i in range(len(basis))))
        for i, psi_i in enumerate(basis):
            assert_zero(
                f"{label} order-omega0 Galerkin equation mode {modes[i]}",
                sp.integrate(psi_i * (parent_apply(resp0) - proj), (w, -sp.oo, sp.oo)),
            )
            assert_zero(
                f"{label} order-omega2 Galerkin equation mode {modes[i]}",
                sp.integrate(psi_i * (parent_apply(resp2) + mu_eta * resp0), (w, -sp.oo, sp.oo)),
            )
            assert_zero(
                f"{label} order-omega4 Galerkin equation mode {modes[i]}",
                sp.integrate(psi_i * (parent_apply(resp4) + mu_eta * resp2), (w, -sp.oo, sp.oo)),
            )
        p0 = sp.simplify((coeffs.T * a0)[0])
        p2 = sp.simplify((coeffs.T * a2)[0])
        p4 = sp.simplify((coeffs.T * a4)[0])
        return coeffs, residual, p0, p2, p4

    coeffs_B, residual_B, B0_export, B2_export, B4_export = low_series_packet("support", phi_B)
    coeffs_Z, residual_Z, Z0_export, Z2_export, Z4_export = low_series_packet("mixed", phi_Z)
    coeffs_N, residual_N, N0_export, N2_export, N4_export = low_series_packet("outgoing", phi_N)
    assert_zero("support basis projection residual", residual_B)
    assert_zero("outgoing basis projection residual", residual_N)
    assert_nonzero("mixed basis projection residual", residual_Z)
    assert_nonzero("mixed port detects dropped mode-4 contribution", coeffs_Z[2])

    # Extend the Galerkin basis by one even mode and check that the mixed-source
    # projection residual drops exactly while the exported low-frequency
    # coefficients move only slightly.
    extended_modes = (0, 2, 4, 6)
    extended_basis = [basis_mode(n) for n in extended_modes]
    for i, psi_i in enumerate(extended_basis):
        assert_zero(f"extended basis normalization n={extended_modes[i]}", sp.integrate(psi_i**2, (w, -sp.oo, sp.oo)) - 1)
        for j in range(i + 1, len(extended_basis)):
            assert_zero(
                f"extended basis orthogonality {extended_modes[i]}/{extended_modes[j]}",
                sp.integrate(psi_i * extended_basis[j], (w, -sp.oo, sp.oo)),
            )

    def extended_project(phi: sp.Expr) -> tuple[sp.Matrix, sp.Expr]:
        coeffs = sp.Matrix(
            [sp.simplify(sp.integrate(extended_basis[i] * phi, (w, -sp.oo, sp.oo))) for i in range(len(extended_basis))]
        )
        proj = sp.simplify(sum(coeffs[i] * extended_basis[i] for i in range(len(extended_basis))))
        residual = sp.simplify(sp.integrate((phi - proj) ** 2, (w, -sp.oo, sp.oo)))
        return coeffs, residual

    extended_mass_matrix = sp.Matrix(
        [
            [
                sp.integrate(mu_eta * extended_basis[i] * extended_basis[j], (w, -sp.oo, sp.oo))
                for j in range(len(extended_basis))
            ]
            for i in range(len(extended_basis))
        ]
    ).evalf(60)
    extended_stiff_matrix = sp.Matrix(
        [
            [
                sp.integrate(
                    T_w * sp.diff(extended_basis[i], w) * sp.diff(extended_basis[j], w)
                    + (K_eta + 6 * T_omega) * extended_basis[i] * extended_basis[j],
                    (w, -sp.oo, sp.oo),
                )
                for j in range(len(extended_basis))
            ]
            for i in range(len(extended_basis))
        ]
    ).evalf(60)

    def extended_low_series(phi: sp.Expr) -> tuple[sp.Matrix, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
        coeffs_exact, residual = extended_project(phi)
        coeffs_num = coeffs_exact.evalf(60)
        a0 = extended_stiff_matrix.LUsolve(coeffs_num)
        a2 = -extended_stiff_matrix.LUsolve(extended_mass_matrix * a0)
        a4 = -extended_stiff_matrix.LUsolve(extended_mass_matrix * a2)
        p0 = sp.N((coeffs_num.T * a0)[0], 40)
        p2 = sp.N((coeffs_num.T * a2)[0], 40)
        p4 = sp.N((coeffs_num.T * a4)[0], 40)
        return coeffs_exact, residual, p0, p2, p4

    coeffs_B_ext, residual_B_ext, B0_ext, B2_ext, B4_ext = extended_low_series(phi_B)
    coeffs_Z_ext, residual_Z_ext, Z0_ext, Z2_ext, Z4_ext = extended_low_series(phi_Z)
    coeffs_N_ext, residual_N_ext, N0_ext, N2_ext, N4_ext = extended_low_series(phi_N)
    assert_zero("extended support projection residual", residual_B_ext)
    assert_zero("extended outgoing projection residual", residual_N_ext)
    assert_zero("mixed projection residual improvement", residual_Z_ext - residual_Z + 5 * sp.sqrt(sp.pi) / 17496)
    assert_small_nonzero("support B0 basis extension drift", B0_ext - sp.N(B0_export, 40), 2e-4)
    assert_small_nonzero("support B2 basis extension drift", B2_ext - sp.N(B2_export, 40), 2e-4)
    assert_small_nonzero("support B4 basis extension drift", B4_ext - sp.N(B4_export, 40), 2e-4)
    assert_small_nonzero("mixed Z0 basis extension drift", Z0_ext - sp.N(Z0_export, 40), 2e-4)
    assert_small_nonzero("mixed Z2 basis extension drift", Z2_ext - sp.N(Z2_export, 40), 2e-4)
    assert_small_nonzero("mixed Z4 basis extension drift", Z4_ext - sp.N(Z4_export, 40), 2e-4)
    assert_small_nonzero("outgoing N0 basis extension drift", N0_ext - sp.N(N0_export, 40), 2e-4)
    assert_small_nonzero("outgoing N2 basis extension drift", N2_ext - sp.N(N2_export, 40), 2e-4)
    assert_small_nonzero("outgoing N4 basis extension drift", N4_ext - sp.N(N4_export, 40), 2e-4)

    def adapted_basis_growth(mode_tuple: tuple[int, ...], phi: sp.Expr) -> tuple[np.ndarray, sp.Expr, float, float, float]:
        local_basis = [basis_mode(n) for n in mode_tuple]
        mass_exact = sp.Matrix(
            [
                [
                    sp.integrate(mu_eta * local_basis[i] * local_basis[j], (w, -sp.oo, sp.oo))
                    for j in range(len(local_basis))
                ]
                for i in range(len(local_basis))
            ]
        )
        stiff_exact = sp.Matrix(
            [
                [
                    sp.integrate(
                        T_w * sp.diff(local_basis[i], w) * sp.diff(local_basis[j], w)
                        + (K_eta + 6 * T_omega) * local_basis[i] * local_basis[j],
                        (w, -sp.oo, sp.oo),
                    )
                    for j in range(len(local_basis))
                ]
                for i in range(len(local_basis))
            ]
        )
        source_exact = sp.Matrix([sp.integrate(local_basis[i] * phi, (w, -sp.oo, sp.oo)) for i in range(len(local_basis))])
        proj = sp.simplify(sum(source_exact[i] * local_basis[i] for i in range(len(local_basis))))
        residual = sp.simplify(sp.integrate((phi - proj) ** 2, (w, -sp.oo, sp.oo)))

        mass_num = np.array(mass_exact.evalf(60).tolist(), dtype=float)
        stiff_num = np.array(stiff_exact.evalf(60).tolist(), dtype=float)
        source_num = np.array(source_exact.evalf(60), dtype=float).reshape(-1)

        L = np.linalg.cholesky(mass_num)
        Linv = np.linalg.inv(L)
        adapted_operator = Linv @ stiff_num @ Linv.T
        evals, U = np.linalg.eigh(adapted_operator)
        modes_num = np.linalg.solve(L.T, U)
        orth_error = float(np.max(np.abs(modes_num.T @ mass_num @ modes_num - np.eye(len(mode_tuple)))))
        diag_error = float(np.max(np.abs(modes_num.T @ stiff_num @ modes_num - np.diag(evals))))
        assert_numeric_small(f"adapted-basis M-orthonormality {mode_tuple}", orth_error, 1e-10)
        assert_numeric_small(f"adapted-basis diagonalization {mode_tuple}", diag_error, 1e-10)

        gamma = modes_num.T @ source_num
        p0 = float(np.sum(gamma**2 / evals))
        p2 = float(-np.sum(gamma**2 / evals**2))
        p4 = float(np.sum(gamma**2 / evals**3))
        return evals, residual, p0, p2, p4

    lam3_Z, residual_Z_num3, Z0_adapt3, Z2_adapt3, Z4_adapt3 = adapted_basis_growth((0, 2, 4), phi_Z)
    lam4_Z, residual_Z_num4, Z0_adapt4, Z2_adapt4, Z4_adapt4 = adapted_basis_growth((0, 2, 4, 6), phi_Z)
    lam5_Z, residual_Z_num5, Z0_adapt5, Z2_adapt5, Z4_adapt5 = adapted_basis_growth((0, 2, 4, 6, 8), phi_Z)
    lam5_B, residual_B_num5, B0_adapt5, B2_adapt5, B4_adapt5 = adapted_basis_growth((0, 2, 4, 6, 8), phi_B)
    lam5_N, residual_N_num5, N0_adapt5, N2_adapt5, N4_adapt5 = adapted_basis_growth((0, 2, 4, 6, 8), phi_N)
    assert_zero("adapted-basis 3-mode mixed residual", residual_Z_num3 - residual_Z)
    assert_zero("adapted-basis 4-mode mixed residual", residual_Z_num4 - residual_Z_ext)
    assert_zero(
        "adapted-basis 5-mode mixed residual",
        residual_Z_num5 - sp.sqrt(sp.pi) * (-890747 + 629856 * sp.sqrt(2)) / 1259712,
    )
    assert_zero("adapted-basis 5-mode support residual", residual_B_num5)
    assert_zero("adapted-basis 5-mode outgoing residual", residual_N_num5)
    assert_zero("mixed residual improvement 4-to-5", residual_Z_num5 - residual_Z_num4 + 35 * sp.sqrt(sp.pi) / 1259712)
    assert_small_nonzero("adapted Z0 drift 3-to-4", Z0_adapt4 - Z0_adapt3, 2e-4)
    assert_small_nonzero("adapted Z2 drift 3-to-4", Z2_adapt4 - Z2_adapt3, 2e-4)
    assert_small_nonzero("adapted Z4 drift 3-to-4", Z4_adapt4 - Z4_adapt3, 2e-4)
    assert_small_nonzero("adapted B0 drift 4-to-5", B0_adapt5 - B0_ext, 3e-5)
    assert_small_nonzero("adapted B2 drift 4-to-5", B2_adapt5 - B2_ext, 3e-5)
    assert_small_nonzero("adapted B4 drift 4-to-5", B4_adapt5 - B4_ext, 3e-5)
    assert_small_nonzero("adapted Z0 drift 4-to-5", Z0_adapt5 - Z0_adapt4, 3e-5)
    assert_small_nonzero("adapted Z2 drift 4-to-5", Z2_adapt5 - Z2_adapt4, 3e-5)
    assert_small_nonzero("adapted Z4 drift 4-to-5", Z4_adapt5 - Z4_adapt4, 3e-5)
    assert_small_nonzero("adapted N0 drift 4-to-5", N0_adapt5 - N0_ext, 3e-5)
    assert_small_nonzero("adapted N2 drift 4-to-5", N2_adapt5 - N2_ext, 3e-5)
    assert_small_nonzero("adapted N4 drift 4-to-5", N4_adapt5 - N4_ext, 3e-5)
    assert_numeric_small("exact-vs-adapted mixed Z0", Z0_adapt3 - float(sp.N(Z0_export, 40)), 1e-10)
    assert_numeric_small("exact-vs-adapted mixed Z2", Z2_adapt3 - float(sp.N(Z2_export, 40)), 1e-10)
    assert_numeric_small("exact-vs-adapted mixed Z4", Z4_adapt3 - float(sp.N(Z4_export, 40)), 1e-10)

    export_packet = {
        KSigma: sp.simplify(K_export),
        MSigma: sp.simplify(M_export),
        B0: sp.simplify(B0_export),
        B2: sp.simplify(B2_export),
        B4: sp.simplify(B4_export),
        Z0: sp.simplify(Z0_export),
        Z2: sp.simplify(Z2_export),
        Z4: sp.simplify(Z4_export),
        N0: sp.simplify(N0_export),
        N2: sp.simplify(N2_export),
        N4: sp.simplify(N4_export),
        mhat0: 1,
        G: 1,
        cs: 1,
        a: 1,
        c: 1,
    }
    branch_metadata = {
        "branch_id": "v2_local_parent_background_galerkin_demo",
        "pre_target_freeze": True,
        "target_blind": True,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "geometry_status": "reduced_parent_background_profile_demo",
        "basis_family": [0, 2, 4, 6, 8],
        "frozen_branch_parameters": {
            "R0_amplitude": 1.0,
            "R0_inverse_width": 1.0,
            "beta_inverse_width": 1.0,
            "mu_shape_quadratic": 0.25,
            "Tw_shape_quadratic": 1.0 / 6.0,
            "TOmega_shape_quadratic": 0.125,
        },
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]
    export_packet_4 = {
        **export_packet,
        B0: sp.N(B0_ext, 40),
        B2: sp.N(B2_ext, 40),
        B4: sp.N(B4_ext, 40),
        Z0: sp.N(Z0_ext, 40),
        Z2: sp.N(Z2_ext, 40),
        Z4: sp.N(Z4_ext, 40),
        N0: sp.N(N0_ext, 40),
        N2: sp.N(N2_ext, 40),
        N4: sp.N(N4_ext, 40),
    }
    export_packet_5 = {
        **export_packet,
        B0: sp.Float(B0_adapt5, 30),
        B2: sp.Float(B2_adapt5, 30),
        B4: sp.Float(B4_adapt5, 30),
        Z0: sp.Float(Z0_adapt5, 30),
        Z2: sp.Float(Z2_adapt5, 30),
        Z4: sp.Float(Z4_adapt5, 30),
        N0: sp.Float(N0_adapt5, 30),
        N2: sp.Float(N2_adapt5, 30),
        N4: sp.Float(N4_adapt5, 30),
    }

    def numeric_target_residues(packet: dict[sp.Symbol, sp.Expr]) -> tuple[float, float, float, float]:
        return (
            float(sp.N(one_pole_residue.subs(packet), 40)),
            float(sp.N((mhat0**2 * P0 - static_target).subs(packet), 40)),
            float(sp.N(P2.subs(packet), 40)),
            float(sp.N(P4.subs(packet), 40)),
        )

    one_pole_residue_export = sp.simplify(one_pole_residue.subs(export_packet))
    static_residue_export = sp.simplify((mhat0**2 * P0 - static_target).subs(export_packet))
    P2_residue_export = sp.simplify(P2.subs(export_packet))
    P4_residue_export = sp.simplify(P4.subs(export_packet))
    assert_nonzero("untuned branch export one-pole residue", one_pole_residue_export)
    assert_nonzero("untuned branch export static-normalization residue", static_residue_export)
    assert_nonzero("untuned branch export P2 residue", P2_residue_export)
    assert_nonzero("untuned branch export P4 residue", P4_residue_export)
    residues_3 = numeric_target_residues(export_packet)
    residues_4 = numeric_target_residues(export_packet_4)
    residues_5 = numeric_target_residues(export_packet_5)
    assert_numeric_small_nonzero("one-pole residue drift 3-to-4", residues_4[0] - residues_3[0], 2e-4)
    assert_numeric_small_nonzero("static residue drift 3-to-4", residues_4[1] - residues_3[1], 2e-4)
    assert_numeric_small_nonzero("P2 residue drift 3-to-4", residues_4[2] - residues_3[2], 2e-4)
    assert_numeric_small_nonzero("P4 residue drift 3-to-4", residues_4[3] - residues_3[3], 2e-4)
    assert_numeric_small_nonzero("one-pole residue drift 4-to-5", residues_5[0] - residues_4[0], 2e-4)
    assert_numeric_small_nonzero("static residue drift 4-to-5", residues_5[1] - residues_4[1], 2e-4)
    assert_numeric_small_nonzero("P2 residue drift 4-to-5", residues_5[2] - residues_4[2], 2e-4)
    assert_numeric_small_nonzero("P4 residue drift 4-to-5", residues_5[3] - residues_4[3], 2e-4)

    weak_packet = {
        **export_packet,
        B01: 0,
        Z01: 0,
        B21: 0,
        Z21: 0,
        B41: 1,
        Z41: 2,
    }
    weak_packet.update({
        dKSigma: canonical[dKSigma].subs(weak_packet),
        dMSigma: canonical[dMSigma].subs(weak_packet),
    })
    weak_packet[N01] = sp.simplify(N01_for_zero_Xi.subs(weak_packet))
    assert_zero("concrete weak K1", K1.subs(weak_packet))
    assert_zero("concrete weak H_even", H_even.subs(weak_packet))
    assert_zero("concrete weak Xi1", Xi1.subs(weak_packet))

    print("STEP 19 ACTUAL BRANCH EXPORT CHECKLIST AUDIT")
    print("Checked symbolic packet tests plus one branch-derived wall export and one untuned port-derived packet.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  geometry_status =", branch_metadata["geometry_status"])
    print("Branch-derived wall export:")
    print("  R0(w)   = exp(-w^2/2)")
    print("  mu_eta  =", sp.sstr(mu_eta))
    print("  T_w     =", sp.sstr(T_w))
    print("  T_omega =", sp.sstr(T_omega))
    print("  K_eta   =", sp.sstr(K_eta))
    print("  M_Sigma =", sp.sstr(M_export))
    print("  K_Sigma =", sp.sstr(K_export))
    print("Parent-background Galerkin port export on basis n=(0,2,4):")
    print("  det(M_parent) =", sp.sstr(sp.factor(mass_matrix.det())))
    print("  det(K_parent) =", sp.sstr(sp.factor(stiff_matrix.det())))
    print("  source coeffs B =", sp.sstr(coeffs_B.T))
    print("  source coeffs Z =", sp.sstr(coeffs_Z.T))
    print("  source coeffs N =", sp.sstr(coeffs_N.T))
    print("  projection residuals =", sp.sstr(residual_B), sp.sstr(residual_Z), sp.sstr(residual_N))
    print("  B0, B2, B4 =", sp.N(B0_export, 16), sp.N(B2_export, 16), sp.N(B4_export, 16))
    print("  Z0, Z2, Z4 =", sp.N(Z0_export, 16), sp.N(Z2_export, 16), sp.N(Z4_export, 16))
    print("  N0, N2, N4 =", sp.N(N0_export, 16), sp.N(N2_export, 16), sp.N(N4_export, 16))
    print("  extended-basis mixed residual =", sp.sstr(residual_Z_ext))
    print("  4-mode B0, B2, B4 =", sp.N(B0_ext, 16), sp.N(B2_ext, 16), sp.N(B4_ext, 16))
    print("  4-mode Z0, Z2, Z4 =", sp.N(Z0_ext, 16), sp.N(Z2_ext, 16), sp.N(Z4_ext, 16))
    print("  4-mode N0, N2, N4 =", sp.N(N0_ext, 16), sp.N(N2_ext, 16), sp.N(N4_ext, 16))
    print("Operator-adapted generalized eigenvalues:")
    print("  3-mode lambdas =", [round(x, 12) for x in lam3_Z.tolist()])
    print("  4-mode lambdas =", [round(x, 12) for x in lam4_Z.tolist()])
    print("  5-mode lambdas =", [round(x, 12) for x in lam5_Z.tolist()])
    print("  5-mode mixed residual =", sp.sstr(residual_Z_num5))
    print("  adapted 5-mode B0, B2, B4 =", round(B0_adapt5, 16), round(B2_adapt5, 16), round(B4_adapt5, 16))
    print("  adapted 5-mode Z0, Z2, Z4 =", round(Z0_adapt5, 16), round(Z2_adapt5, 16), round(Z4_adapt5, 16))
    print("  adapted 5-mode N0, N2, N4 =", round(N0_adapt5, 16), round(N2_adapt5, 16), round(N4_adapt5, 16))
    print("  4-mode isotropic residues =", residues_4)
    print("  adapted 5-mode isotropic residues =", residues_5)
    print("Untuned parent-background packet:")
    print("  residual packet labels = (R_pole, R_norm, R_P2, R_P4)")
    print("  R_pole =", sp.N(one_pole_residue_export, 16))
    print("  R_norm =", sp.N(static_residue_export, 16))
    print("  R_P2   =", sp.N(P2_residue_export, 16))
    print("  R_P4   =", sp.N(P4_residue_export, 16))
    print("  one-pole residue =", sp.N(one_pole_residue_export, 16))
    print("  static residue   =", sp.N(static_residue_export, 16))
    print("  P2 residue       =", sp.N(P2_residue_export, 16))
    print("  P4 residue       =", sp.N(P4_residue_export, 16))
    print("  weak dK, dM, N01 =", sp.sstr(weak_packet[dKSigma]), sp.sstr(weak_packet[dMSigma]), sp.sstr(weak_packet[N01]))
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
