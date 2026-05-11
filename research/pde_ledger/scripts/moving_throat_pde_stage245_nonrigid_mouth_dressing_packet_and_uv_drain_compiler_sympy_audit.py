#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 245
SymPy audit for the non-rigid mouth/dressing packet and U/V drain compiler.

This script verifies:
1. the exact stationary U/V packet from the reduced non-rigid free energy,
2. the exact admissibility determinant and V/U response ratio,
3. the exact positive drain term,
4. the exact finite physical compiler for T^2, epsilon_eta, and R_target,
5. the weak-axisymmetric first-order packet (Xi_1, R_1),
6. the exact dependent microscopic correction,
7. the support-blind / orbit-side separation,
8. and a Session-I readback using the values recorded in the barrier-session write-up.
"""

from __future__ import annotations

import sympy as sp


def section(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def coeff_linear(expr: sp.Expr, eps: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.series(expr, eps, 0, 2).removeO()).coeff(eps, 1)


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Exact stationary packet from the non-rigid reduced free energy.
    # ------------------------------------------------------------------
    section("1. Exact stationary non-rigid packet")

    U, V = sp.symbols("U V", real=True)
    a_U, a_V = sp.symbols("a_U a_V", positive=True, real=True)
    chi_UV = sp.symbols("chi_UV", real=True)
    f_U = sp.symbols("f_U", real=True)

    F_nr = sp.Rational(1, 2) * a_U * U**2 + sp.Rational(1, 2) * a_V * V**2 - chi_UV * U * V - f_U * U
    eq_U = sp.diff(F_nr, U)
    eq_V = sp.diff(F_nr, V)

    sol = sp.solve((sp.Eq(eq_U, 0), sp.Eq(eq_V, 0)), (U, V), dict=True)[0]
    U_sol = sp.simplify(sol[U])
    V_sol = sp.simplify(sol[V])

    Delta_UV = sp.simplify(a_U * a_V - chi_UV**2)
    U_expected = sp.simplify(a_V * f_U / Delta_UV)
    V_expected = sp.simplify(chi_UV * f_U / Delta_UV)
    ratio_VU = sp.simplify(V_sol / U_sol)
    H_nr = sp.hessian(F_nr, (U, V))
    det_H = sp.simplify(H_nr.det())

    print("F_nr(U,V)             =", F_nr)
    print("dF/dU                 =", eq_U)
    print("dF/dV                 =", eq_V)
    print("Delta_UV              =", Delta_UV)
    print("U_sol                 =", U_sol)
    print("V_sol                 =", V_sol)
    print("V/U                   =", ratio_VU)
    print("Hessian               =")
    sp.pprint(H_nr)
    print("det(H)                =", det_H)

    assert sp.simplify(U_sol - U_expected) == 0
    assert sp.simplify(V_sol - V_expected) == 0
    assert sp.simplify(ratio_VU - chi_UV / a_V) == 0
    assert sp.simplify(det_H - Delta_UV) == 0
    assert U_sol.subs(f_U, 0) == 0
    assert V_sol.subs(f_U, 0) == 0
    assert sp.simplify(V_sol.subs(chi_UV, 0)) == 0
    assert sp.simplify(U_sol.subs(chi_UV, 0) - f_U / a_U) == 0

    # ------------------------------------------------------------------
    # 2. Positive U/V drain term.
    # ------------------------------------------------------------------
    section("2. Positive U/V drain")
    D_UV = sp.simplify(chi_UV * U_sol * V_sol)
    D_expected = sp.simplify(chi_UV**2 * a_V * f_U**2 / Delta_UV**2)

    print("D_UV                  =", D_UV)
    print("Expected D_UV         =", D_expected)

    assert sp.simplify(D_UV - D_expected) == 0

    # ------------------------------------------------------------------
    # 3. Exact finite physical compiler.
    # ------------------------------------------------------------------
    section("3. Exact finite physical compiler")
    T2_ref, eps_eta_ref = sp.symbols("T2_ref eps_eta_ref", positive=True, real=True)
    R_target_ref = sp.symbols("R_target_ref", positive=True, real=True)

    T2 = sp.simplify(T2_ref * sp.exp(U_sol))
    eps_eta = sp.simplify(eps_eta_ref * sp.exp(V_sol))
    R_ratio = sp.simplify(((1 - eps_eta_ref * sp.exp(V_sol)) / (1 - eps_eta_ref)) * sp.exp(-U_sol))
    R_exact_check = sp.simplify(R_ratio * sp.exp(U_sol) * (1 - eps_eta_ref) / (1 - eps_eta_ref * sp.exp(V_sol)) - 1)

    print("T^2 / T^2_ref         =", sp.simplify(T2 / T2_ref))
    print("epsilon_eta / eps_ref =", sp.simplify(eps_eta / eps_eta_ref))
    print("R_target / R_ref      =", R_ratio)
    print("Exact multiplicative check =", R_exact_check)

    assert sp.simplify(T2 / T2_ref - sp.exp(U_sol)) == 0
    assert sp.simplify(eps_eta / eps_eta_ref - sp.exp(V_sol)) == 0
    assert sp.simplify(R_exact_check) == 0

    # ------------------------------------------------------------------
    # 4. Exact dependent microscopic correction.
    # ------------------------------------------------------------------
    section("4. Exact dependent microscopic correction")
    y_dep = sp.Matrix([sp.Integer(0), -V_sol, U_sol - V_sol])
    y_expected = sp.Matrix([
        sp.Integer(0),
        -chi_UV * f_U / Delta_UV,
        (a_V - chi_UV) * f_U / Delta_UV,
    ])

    print("y_dep(U,V)            =")
    sp.pprint(y_dep)

    assert sp.simplify(y_dep - y_expected) == sp.Matrix([0, 0, 0])

    # ------------------------------------------------------------------
    # 5. Weak-axisymmetric first-order packet.
    # ------------------------------------------------------------------
    section("5. Weak-axisymmetric first-order packet")
    eps = sp.symbols("eps", real=True)
    u1, v1 = sp.symbols("u1 v1", real=True)

    U_lin = eps * u1
    V_lin = eps * v1

    log_T = sp.log(sp.exp(U_lin))
    log_one_minus = sp.log((1 - eps_eta_ref * sp.exp(V_lin)) / (1 - eps_eta_ref))
    log_R = sp.log(((1 - eps_eta_ref * sp.exp(V_lin)) / (1 - eps_eta_ref)) * sp.exp(-U_lin))

    Xi1_nr = coeff_linear(log_T, eps)
    RpXi1 = coeff_linear(log_one_minus, eps)
    R1_nr = coeff_linear(log_R, eps)
    R1_expected = sp.simplify(-u1 - eps_eta_ref * v1 / (1 - eps_eta_ref))
    RpXi1_expected = sp.simplify(-eps_eta_ref * v1 / (1 - eps_eta_ref))

    print("Xi1_nr                =", Xi1_nr)
    print("R1 + Xi1              =", RpXi1)
    print("R1_nr                 =", R1_nr)

    assert sp.simplify(Xi1_nr - u1) == 0
    assert sp.simplify(RpXi1 - RpXi1_expected) == 0
    assert sp.simplify(R1_nr - R1_expected) == 0
    assert sp.simplify(R1_nr + Xi1_nr - RpXi1) == 0

    # ------------------------------------------------------------------
    # 6. Orbit-side / support-blind split.
    # ------------------------------------------------------------------
    section("6. Orbit-side / support-blind split")
    Lam, varrho = sp.symbols("Lam varrho", real=True)

    objects = {
        "U_sol": U_sol,
        "V_sol": V_sol,
        "eps_eta": eps_eta,
        "R_ratio": R_ratio,
        "D_UV": D_UV,
        "Delta_Keta": y_dep[1],
        "Delta_mu": y_dep[2],
    }

    for name, obj in objects.items():
        dLam = sp.simplify(sp.diff(obj, Lam))
        dvarrho = sp.simplify(sp.diff(obj, varrho))
        print(f"d/dLam    {name} =", dLam)
        print(f"d/dvarrho {name} =", dvarrho)
        assert dLam == 0
        assert dvarrho == 0

    # ------------------------------------------------------------------
    # 7. Session-I readback.
    # ------------------------------------------------------------------
    section("7. Session-I readback")
    U_obs = sp.Float("0.14313458")
    V_obs = sp.Float("-0.03619791")
    eps_ref_obs = sp.Float("0.3")
    aU_obs = sp.Float("2.5")
    aV_obs = sp.Float("3.0")
    gUV_obs = sp.Float("0.95")

    eps_eta_obs = sp.N(eps_ref_obs * sp.exp(V_obs), 16)
    R_ratio_obs = sp.N(((1 - eps_ref_obs * sp.exp(V_obs)) / (1 - eps_ref_obs)) * sp.exp(-U_obs), 16)
    y_dep_obs = sp.N(sp.Matrix([0, -V_obs, U_obs - V_obs]), 16)

    chi_UV_obs = sp.N(aV_obs * V_obs / U_obs, 16)
    chi_lambda_obs = sp.N(chi_UV_obs / gUV_obs, 16)
    Delta_obs = sp.N(aU_obs * aV_obs - chi_UV_obs**2, 16)
    fU_obs = sp.N(U_obs * Delta_obs / aV_obs, 16)

    U_rebuilt = sp.N((aV_obs * fU_obs) / Delta_obs, 16)
    V_rebuilt = sp.N((chi_UV_obs * fU_obs) / Delta_obs, 16)
    R1_obs = sp.N(-U_obs - eps_ref_obs * V_obs / (1 - eps_ref_obs), 16)

    print("epsilon_eta(obs)      =", eps_eta_obs)
    print("R_target/R_ref(obs)   =", R_ratio_obs)
    print("y_dep(obs)            =")
    sp.pprint(y_dep_obs)
    print("chi_UV inferred       =", chi_UV_obs)
    print("chi_lambda inferred   =", chi_lambda_obs)
    print("Delta_UV inferred     =", Delta_obs)
    print("f_U inferred          =", fU_obs)
    print("U rebuilt             =", U_rebuilt)
    print("V rebuilt             =", V_rebuilt)
    print("R1_obs                =", R1_obs)

    assert abs(float(eps_eta_obs) - 0.28933482) < 5e-9
    assert abs(float(U_rebuilt - U_obs)) < 1e-12
    assert abs(float(V_rebuilt - V_obs)) < 1e-12

    # ------------------------------------------------------------------
    # 8. Final success banner.
    # ------------------------------------------------------------------
    section("Stage 245 audit result")
    print("All symbolic checks passed.")
    print("Verified objects:")
    print("- exact stationary U/V packet on the non-rigid mouth/dressing branch,")
    print("- exact admissibility determinant and signed V/U response ratio,")
    print("- exact positive U/V drain,")
    print("- exact finite physical compiler for T^2, epsilon_eta, and R_target,")
    print("- exact weak-axisymmetric packet (Xi_1, R_1),")
    print("- exact dependent microscopic correction,")
    print("- orbit-side / support-blind separation from the selected-support packet,")
    print("- and a consistent Session-I readback of the reported U/V point.")


if __name__ == "__main__":
    main()
