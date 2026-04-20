#!/usr/bin/env python3
"""
moving_throat_pde_stage002_breathing_reduction_sympy_audit.py

Stage 002 SymPy audit for the breathing reduction back to the old (a, L)
closure and the uncoupled grouped real P2 reduction.

Scope
-----
This script now verifies the actual reduction steps claimed in Stage 002:

  • the Y_00 mouth-average bridge behind q_00 = 2 sqrt(pi) delta a,
  • insertion of the two-mode axisymmetric ansatz into the Stage 001 wall action,
  • the resulting overlap-integral formulas for M_AB and K_AB,
  • the conservative Euler-Lagrange matrix form built from those overlap integrals,
  • the grouped-real P2 orthonormality / angular stiffness matrix,
  • and the resulting isotropic degeneracy before coupling.
"""

from __future__ import annotations

import sympy as sp
from sympy.calculus.euler import euler_equations


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.trigsimp(sp.expand(z))))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
        return

    simplified = sp.simplify(sp.trigsimp(sp.expand(expr)))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


def lap_s2(expr: sp.Expr, theta: sp.Symbol, phi: sp.Symbol) -> sp.Expr:
    return sp.simplify(
        (1 / sp.sin(theta)) * sp.diff(sp.sin(theta) * sp.diff(expr, theta), theta)
        + sp.diff(expr, phi, 2) / sp.sin(theta) ** 2
    )


def grad_s2_inner(y1: sp.Expr, y2: sp.Expr, theta: sp.Symbol, phi: sp.Symbol) -> sp.Expr:
    return sp.simplify(
        sp.diff(y1, theta) * sp.diff(y2, theta)
        + sp.diff(y1, phi) * sp.diff(y2, phi) / sp.sin(theta) ** 2
    )


def normalization_bridge_audit() -> None:
    banner("SECTION I — MONOPOLE NORMALIZATION BRIDGE")

    theta, phi = sp.symbols("theta phi", real=True)
    dOmega = sp.sin(theta)
    y00 = sp.Rational(1, 2) / sp.sqrt(sp.pi)

    avg_y00 = sp.simplify(
        sp.integrate(sp.integrate(y00 * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)) / (4 * sp.pi)
    )
    norm_y00 = sp.simplify(
        sp.integrate(sp.integrate(y00 * y00 * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi))
    )

    expect_zero("Y00 mouth average - 1/(2 sqrt(pi))", avg_y00 - 1 / (2 * sp.sqrt(sp.pi)))
    expect_zero("norm(Y00) - 1", norm_y00 - 1)

    delta_a, q00 = sp.symbols("delta_a q00", real=True)
    mouth_avg = sp.simplify(q00 * avg_y00)
    expect_zero("mouth average from q00 Y00 - q00/(2 sqrt(pi))", mouth_avg - q00 / (2 * sp.sqrt(sp.pi)))
    print("Identifying the physical mouth-average shift with delta_a then gives q00 = 2 sqrt(pi) delta_a.")

    angular_prefactor = sp.simplify(
        sp.integrate(
            sp.integrate((2 * sp.sqrt(sp.pi) * y00) ** 2 * dOmega, (phi, 0, 2 * sp.pi)),
            (theta, 0, sp.pi),
        )
    )
    expect_zero("angular prefactor from (2 sqrt(pi) Y00)^2 - 4 pi", angular_prefactor - 4 * sp.pi)


def breathing_matrix_reduction_audit() -> None:
    banner("SECTION II — AXISYMMETRIC TWO-MODE REDUCTION FROM THE STAGE 001 ACTION")

    theta, phi, w, wL, wR = sp.symbols("theta phi w w_L w_R", real=True)
    dOmega = sp.sin(theta)
    y00 = sp.Rational(1, 2) / sp.sqrt(sp.pi)

    mu_eta = sp.Function("mu_eta")(w)
    T_w = sp.Function("T_w")(w)
    K0 = sp.Function("K_0")(w)
    alpha_a = sp.Function("alpha_a")(w)
    alpha_L = sp.Function("alpha_L")(w)

    da, dL, dadt, dLdt = sp.symbols("delta_a delta_L delta_a_dot delta_L_dot", real=True)

    axisym = alpha_a * da + alpha_L * dL
    axisym_t = alpha_a * dadt + alpha_L * dLdt
    axisym_w = sp.diff(alpha_a, w) * da + sp.diff(alpha_L, w) * dL

    eta = 2 * sp.sqrt(sp.pi) * axisym * y00
    eta_t = 2 * sp.sqrt(sp.pi) * axisym_t * y00
    eta_w = 2 * sp.sqrt(sp.pi) * axisym_w * y00

    lw = sp.simplify(
        sp.Rational(1, 2)
        * sp.integrate(
            sp.integrate(
                (mu_eta * eta_t**2 - T_w * eta_w**2 - K0 * eta**2) * dOmega,
                (phi, 0, 2 * sp.pi),
            ),
            (theta, 0, sp.pi),
        )
    )

    Q = sp.Matrix([da, dL])
    Qdot = sp.Matrix([dadt, dLdt])

    M_integrand = 4 * sp.pi * sp.Matrix(
        [
            [mu_eta * alpha_a**2, mu_eta * alpha_a * alpha_L],
            [mu_eta * alpha_a * alpha_L, mu_eta * alpha_L**2],
        ]
    )
    K_integrand = 4 * sp.pi * sp.Matrix(
        [
            [
                T_w * sp.diff(alpha_a, w) ** 2 + K0 * alpha_a**2,
                T_w * sp.diff(alpha_a, w) * sp.diff(alpha_L, w) + K0 * alpha_a * alpha_L,
            ],
            [
                T_w * sp.diff(alpha_a, w) * sp.diff(alpha_L, w) + K0 * alpha_a * alpha_L,
                T_w * sp.diff(alpha_L, w) ** 2 + K0 * alpha_L**2,
            ],
        ]
    )

    lw_target = sp.simplify(sp.Rational(1, 2) * ((Qdot.T * M_integrand * Qdot)[0] - (Q.T * K_integrand * Q)[0]))
    expect_zero("reduced Lagrangian density from the action - target overlap form", lw - lw_target)

    M_matrix = sp.Matrix(
        [
            [sp.Integral(M_integrand[0, 0], (w, wL, wR)), sp.Integral(M_integrand[0, 1], (w, wL, wR))],
            [sp.Integral(M_integrand[1, 0], (w, wL, wR)), sp.Integral(M_integrand[1, 1], (w, wL, wR))],
        ]
    )
    K_matrix = sp.Matrix(
        [
            [sp.Integral(K_integrand[0, 0], (w, wL, wR)), sp.Integral(K_integrand[0, 1], (w, wL, wR))],
            [sp.Integral(K_integrand[1, 0], (w, wL, wR)), sp.Integral(K_integrand[1, 1], (w, wL, wR))],
        ]
    )
    lred = sp.Integral(lw, (w, wL, wR))
    lred_target = sp.simplify(sp.Rational(1, 2) * ((Qdot.T * M_matrix * Qdot)[0] - (Q.T * K_matrix * Q)[0]))
    expect_zero("formal reduced Lagrangian - boxed matrix form", lred - lred_target)

    print("Recovered overlap matrices:")
    sp.pprint(M_matrix)
    sp.pprint(K_matrix)

    t = sp.symbols("t", real=True)
    qa = sp.Function("q_a")
    qL = sp.Function("q_L")
    Maa, MaL, MLL = M_matrix[0, 0], M_matrix[0, 1], M_matrix[1, 1]
    Kaa, KaL, KLL = K_matrix[0, 0], K_matrix[0, 1], K_matrix[1, 1]

    lred_time = sp.Rational(1, 2) * (
        Maa * sp.diff(qa(t), t) ** 2
        + 2 * MaL * sp.diff(qa(t), t) * sp.diff(qL(t), t)
        + MLL * sp.diff(qL(t), t) ** 2
        - Kaa * qa(t) ** 2
        - 2 * KaL * qa(t) * qL(t)
        - KLL * qL(t) ** 2
    )

    el_a = euler_equations(lred_time, qa(t), [t])[0]
    el_L = euler_equations(lred_time, qL(t), [t])[0]

    expect_zero(
        "Euler-Lagrange equation for q_a",
        el_a.lhs + Maa * sp.diff(qa(t), t, 2) + MaL * sp.diff(qL(t), t, 2) + Kaa * qa(t) + KaL * qL(t),
    )
    expect_zero(
        "Euler-Lagrange equation for q_L",
        el_L.lhs + MaL * sp.diff(qa(t), t, 2) + MLL * sp.diff(qL(t), t, 2) + KaL * qa(t) + KLL * qL(t),
    )


def grouped_p2_reduction_audit() -> None:
    banner("SECTION III — GROUPED REAL P2 ORTHONORMALITY AND DEGENERACY")

    theta, phi, w, wL, wR = sp.symbols("theta phi w w_L w_R", real=True)
    dOmega = sp.sin(theta)

    y20 = sp.sqrt(5) / (4 * sp.sqrt(sp.pi)) * (3 * sp.cos(theta) ** 2 - 1)
    y21c = sp.sqrt(15) / (2 * sp.sqrt(sp.pi)) * sp.sin(theta) * sp.cos(theta) * sp.cos(phi)
    y21s = sp.sqrt(15) / (2 * sp.sqrt(sp.pi)) * sp.sin(theta) * sp.cos(theta) * sp.sin(phi)
    y22c = sp.sqrt(15) / (4 * sp.sqrt(sp.pi)) * sp.sin(theta) ** 2 * sp.cos(2 * phi)
    y22s = sp.sqrt(15) / (4 * sp.sqrt(sp.pi)) * sp.sin(theta) ** 2 * sp.sin(2 * phi)

    basis = [y20, y21c, y21s, y22c, y22s]

    norm_matrix = sp.Matrix(
        [
            [
                sp.simplify(
                    sp.integrate(sp.integrate(basis[i] * basis[j] * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi))
                )
                for j in range(5)
            ]
            for i in range(5)
        ]
    )
    grad_matrix = sp.Matrix(
        [
            [
                sp.simplify(
                    sp.integrate(
                        sp.integrate(grad_s2_inner(basis[i], basis[j], theta, phi) * dOmega, (phi, 0, 2 * sp.pi)),
                        (theta, 0, sp.pi),
                    )
                )
                for j in range(5)
            ]
            for i in range(5)
        ]
    )

    expect_zero("real P2 norm matrix - identity", norm_matrix - sp.eye(5))
    expect_zero("real P2 angular stiffness matrix - 6 identity", grad_matrix - 6 * sp.eye(5))
    for idx, y in enumerate(basis):
        expect_zero(f"-Delta_S2 basis[{idx}] - 6 basis[{idx}]", -lap_s2(y, theta, phi) - 6 * y)

    mu_eta = sp.Function("mu_eta")(w)
    T_w = sp.Function("T_w")(w)
    T_omega = sp.Function("T_Omega")(w)
    K_eta = sp.Function("K_eta")(w)
    beta2 = sp.Function("beta_2")(w)

    q_syms = sp.Matrix(sp.symbols("q20 q21c q21s q22c q22s", real=True))
    qdot_syms = sp.Matrix(sp.symbols("q20_dot q21c_dot q21s_dot q22c_dot q22s_dot", real=True))

    kinetic = sp.simplify((qdot_syms.T * norm_matrix * qdot_syms)[0])
    radial = sp.simplify((q_syms.T * norm_matrix * q_syms)[0])
    angular = sp.simplify((q_syms.T * grad_matrix * q_syms)[0])

    lw_p2 = sp.simplify(
        sp.Rational(1, 2)
        * (
            mu_eta * beta2**2 * kinetic
            - T_w * sp.diff(beta2, w) ** 2 * radial
            - K_eta * beta2**2 * radial
            - T_omega * beta2**2 * angular
        )
    )
    lw_p2_target = sp.simplify(
        sp.Rational(1, 2)
        * sum(
            mu_eta * beta2**2 * qdot_syms[i] ** 2
            - (T_w * sp.diff(beta2, w) ** 2 + (K_eta + 6 * T_omega) * beta2**2) * q_syms[i] ** 2
            for i in range(5)
        )
    )
    expect_zero("grouped real P2 reduced density - degenerate five-component form", lw_p2 - lw_p2_target)

    t = sp.symbols("t", real=True)
    q = sp.Function("q")
    M2 = sp.Integral(mu_eta * beta2**2, (w, wL, wR))
    K2 = sp.Integral(T_w * sp.diff(beta2, w) ** 2 + (K_eta + 6 * T_omega) * beta2**2, (w, wL, wR))
    l2 = sp.Rational(1, 2) * (M2 * sp.diff(q(t), t) ** 2 - K2 * q(t) ** 2)
    el2 = euler_equations(l2, q(t), [t])[0]
    expect_zero("single-component P2 Euler-Lagrange equation", el2.lhs + M2 * sp.diff(q(t), t, 2) + K2 * q(t))


def main() -> None:
    normalization_bridge_audit()
    breathing_matrix_reduction_audit()
    grouped_p2_reduction_audit()

    banner("STAGE 002 SUMMARY")
    print("Verified with SymPy:")
    print("  • the Y00 normalization bridge behind q00 = 2 sqrt(pi) delta a;")
    print("  • insertion of the two-mode ansatz into the Stage 001 wall action, producing")
    print("    the boxed overlap-integral formulas for M_AB and K_AB;")
    print("  • the conservative breathing-reduction matrix system M_AB Qdd^B + K_AB Q^B")
    print("    using those overlap integrals as coefficients;")
    print("  • real P2 orthonormality, common angular stiffness 6, and the resulting")
    print("    isotropic grouped-real P2 degeneracy before coupling.")


if __name__ == "__main__":
    main()
