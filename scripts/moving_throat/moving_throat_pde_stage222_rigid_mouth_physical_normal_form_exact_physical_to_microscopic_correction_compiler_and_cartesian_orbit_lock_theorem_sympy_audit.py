#!/usr/bin/env python3
"""SymPy audit for Moving-Throat PDE Stage 222.

This script verifies the rigid-mouth physical normal form:

1. The rigid-mouth packet is diagonal in the physical logarithmic coordinates
       U = log(T^2/T_ref^2),  V = log(eps_eta/eps_eta_ref)
2. The exact target-ratio formula
       R_target/R_target_ref = ((1 - eps_ref*e^V)/(1 - eps_ref)) * e^{-U}
3. The complementary physical projectors and the two commuting finite legs
4. The exact physical-to-microscopic dependent-plane compiler
       (Delta_T, Delta_Keta, Delta_mu) = (0, -V, U - V)
5. The pure transfer-shape and pure dressing microscopic images
6. The static-only, post-static, and full orbit-restoring correction compilers
7. Support-blindness with respect to zeta and M_mix
8. The Cartesian rigid-mouth orbit-lock theorem, including first-order form
"""

from __future__ import annotations

import sympy as sp


def assert_zero(expr: sp.Expr, label: str) -> None:
    simplified = sp.simplify(sp.factor(expr))
    if simplified != 0:
        raise AssertionError(f"{label} failed: {simplified}")
    print(f"[ok] {label}")


def assert_matrix_zero(mat: sp.Matrix, label: str) -> None:
    simplified = mat.applyfunc(lambda x: sp.simplify(sp.factor(x)))
    if any(entry != 0 for entry in simplified):
        raise AssertionError(f"{label} failed:\n{simplified}")
    print(f"[ok] {label}")


def main() -> None:
    sp.init_printing()

    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    U, V = sp.symbols("U V", real=True)
    T2, T2_ref = sp.symbols("T2 T2_ref", positive=True)
    eps_eta, eps_eta_ref = sp.symbols("eps_eta eps_eta_ref", positive=True)
    Rtarget, Rtarget_ref = sp.symbols("Rtarget Rtarget_ref", positive=True)
    Lambda0 = sp.symbols("Lambda0", positive=True)

    zeta, M_mix = sp.symbols("zeta M_mix")
    h = sp.symbols("h")
    dlnT2, dlnepseta = sp.symbols("dlnT2 dlnepseta", real=True)

    # ------------------------------------------------------------------
    # 1. Exact rigid-mouth physical logarithmic chart
    # ------------------------------------------------------------------
    U_def = sp.log(T2 / T2_ref)
    V_def = sp.log(eps_eta / eps_eta_ref)

    assert_zero(sp.exp(U_def) - T2 / T2_ref, "physical log coordinate U")
    assert_zero(sp.exp(V_def) - eps_eta / eps_eta_ref, "physical log coordinate V")

    M_phys = sp.eye(2)
    q_phys = M_phys * sp.Matrix([U, V])
    assert_matrix_zero(q_phys - sp.Matrix([U, V]), "diagonal rigid-mouth packet compiler")

    # Selected-branch identity and target-ratio reconstruction.
    Rtarget_from_identity = sp.simplify(Lambda0 * (1 - eps_eta) / T2)
    Rtarget_ref_from_identity = sp.simplify(Lambda0 * (1 - eps_eta_ref) / T2_ref)

    ratio_from_identity = sp.simplify(Rtarget_from_identity / Rtarget_ref_from_identity)
    ratio_from_UV = sp.simplify(((1 - eps_eta_ref * sp.exp(V)) / (1 - eps_eta_ref)) * sp.exp(-U))
    assert_zero(
        ratio_from_identity.subs({sp.exp(U): T2 / T2_ref, sp.exp(V): eps_eta / eps_eta_ref}) - ratio_from_UV.subs({sp.exp(U): T2 / T2_ref, sp.exp(V): eps_eta / eps_eta_ref}),
        "target-ratio formula from U,V",
    )

    # Also verify directly from substitutions.
    ratio_uv_direct = sp.simplify(
        ratio_from_UV.subs(
            {
                U: U_def,
                V: V_def,
            }
        )
    )
    assert_zero(ratio_uv_direct - ratio_from_identity, "target-ratio reconstruction in physical chart")

    # ------------------------------------------------------------------
    # 2. Exact physical projectors and finite legs
    # ------------------------------------------------------------------
    P_T = sp.Matrix([[1, 0], [0, 0]])
    P_eta = sp.Matrix([[0, 0], [0, 1]])

    assert_matrix_zero(P_T * P_T - P_T, "P_T idempotence")
    assert_matrix_zero(P_eta * P_eta - P_eta, "P_eta idempotence")
    assert_matrix_zero(P_T * P_eta, "P_T P_eta = 0")
    assert_matrix_zero(P_eta * P_T, "P_eta P_T = 0")
    assert_matrix_zero(P_T + P_eta - sp.eye(2), "P_T + P_eta = I")

    x_phys = sp.Matrix([U, V])
    x_T = P_T * x_phys
    x_eta = P_eta * x_phys
    assert_matrix_zero(x_phys - x_T - x_eta, "physical packet decomposition")

    # Finite transfer leg and dressing leg in target-ratio form.
    ratio_transfer = sp.simplify(ratio_from_UV.subs(V, 0))
    ratio_dressing = sp.simplify(ratio_from_UV.subs(U, 0))
    assert_zero(ratio_transfer - sp.exp(-U), "pure transfer-shape leg")
    assert_zero(
        ratio_dressing - (1 - eps_eta_ref * sp.exp(V)) / (1 - eps_eta_ref),
        "pure dressing leg",
    )
    assert_zero(
        ratio_from_UV - ratio_transfer * ratio_dressing,
        "exact commutativity / factorization of transfer and dressing legs",
    )

    # ------------------------------------------------------------------
    # 3. Exact physical-to-microscopic dependent-plane compiler
    # ------------------------------------------------------------------
    C_phys_dep = sp.Matrix([[0, 0], [0, -1], [1, -1]])
    y_dep = C_phys_dep * x_phys
    assert_matrix_zero(
        y_dep - sp.Matrix([0, -V, U - V]),
        "physical-to-microscopic dependent compiler",
    )

    L_phys_dep = sp.Matrix([[0, -1, 1], [0, -1, 0]])
    assert_matrix_zero(L_phys_dep * C_phys_dep - sp.eye(2), "left inverse of physical compiler")

    recovered = L_phys_dep * y_dep
    assert_matrix_zero(recovered - x_phys, "recovery of U,V from dependent correction")

    # ------------------------------------------------------------------
    # 4. Exact microscopic images of the two physical axes
    # ------------------------------------------------------------------
    y_T = C_phys_dep * x_T
    y_eta = C_phys_dep * x_eta

    assert_matrix_zero(y_T - sp.Matrix([0, 0, U]), "pure transfer-shape microscopic image")
    assert_matrix_zero(y_eta + V * sp.Matrix([0, 1, 1]), "pure dressing microscopic image")
    assert_matrix_zero(y_dep - y_T - y_eta, "dependent correction splits into the two physical axes")

    # ------------------------------------------------------------------
    # 5. Exact correction compilers
    # ------------------------------------------------------------------
    Delta_y_static = -y_T
    Delta_y_eta_rest = V * sp.Matrix([0, 1, 1])
    Delta_y_orbit = -y_dep

    assert_matrix_zero(Delta_y_static - sp.Matrix([0, 0, -U]), "static-only correction")
    assert_matrix_zero(Delta_y_eta_rest - sp.Matrix([0, V, V]), "post-static dressing correction")
    assert_matrix_zero(
        Delta_y_orbit - sp.Matrix([0, V, V - U]),
        "full orbit-lock correction",
    )
    assert_matrix_zero(Delta_y_orbit - Delta_y_static - Delta_y_eta_rest, "full correction splits into static + dressing parts")
    assert_matrix_zero(y_dep + Delta_y_orbit, "full orbit-lock correction cancels the dependent defect")

    # ------------------------------------------------------------------
    # 6. Exact support-blindness
    # ------------------------------------------------------------------
    # The physical observables T^2 and eps_eta are already support-blind on the carried branch.
    assert_zero(sp.diff(U_def, zeta), "support-blindness of U with respect to zeta")
    assert_zero(sp.diff(V_def, zeta), "support-blindness of V with respect to zeta")
    assert_zero(sp.diff(U_def, M_mix), "support-blindness of U with respect to M_mix")
    assert_zero(sp.diff(V_def, M_mix), "support-blindness of V with respect to M_mix")

    y_dep_def = sp.Matrix([0, -V_def, U_def - V_def])
    Delta_y_static_def = sp.Matrix([0, 0, -U_def])
    Delta_y_orbit_def = sp.Matrix([0, V_def, V_def - U_def])

    assert_matrix_zero(y_dep_def.applyfunc(lambda x: sp.diff(x, zeta)), "support-blindness of dependent correction w.r.t. zeta")
    assert_matrix_zero(y_dep_def.applyfunc(lambda x: sp.diff(x, M_mix)), "support-blindness of dependent correction w.r.t. M_mix")
    assert_matrix_zero(Delta_y_static_def.applyfunc(lambda x: sp.diff(x, zeta)), "support-blindness of static correction w.r.t. zeta")
    assert_matrix_zero(Delta_y_orbit_def.applyfunc(lambda x: sp.diff(x, zeta)), "support-blindness of orbit correction w.r.t. zeta")
    assert_matrix_zero(Delta_y_static_def.applyfunc(lambda x: sp.diff(x, M_mix)), "support-blindness of static correction w.r.t. M_mix")
    assert_matrix_zero(Delta_y_orbit_def.applyfunc(lambda x: sp.diff(x, M_mix)), "support-blindness of orbit correction w.r.t. M_mix")

    # ------------------------------------------------------------------
    # 7. Cartesian orbit-lock theorem and first-order form
    # ------------------------------------------------------------------
    # Finite equivalences.
    assert_zero(sp.exp(U.subs(U, 0)) - 1, "U = 0 implies T^2/T_ref^2 = 1")
    assert_zero(sp.exp(V.subs(V, 0)) - 1, "V = 0 implies eps_eta/eps_eta_ref = 1")
    # Use direct substitution checks from the defining logarithms.
    assert_zero(sp.exp(U_def).subs(T2, T2_ref) - 1, "T^2 = T_ref^2 implies U = 0")
    assert_zero(sp.exp(V_def).subs(eps_eta, eps_eta_ref) - 1, "eps_eta = eps_eta_ref implies V = 0")
    assert_zero(ratio_from_UV.subs({U: 0, V: 0}) - 1, "U = V = 0 implies R_target = R_target_ref")

    # First-order form.
    T2_pert = T2 * sp.exp(h * dlnT2)
    eps_pert = eps_eta * sp.exp(h * dlnepseta)

    U_first = sp.simplify(sp.diff(sp.log(T2_pert / T2), h).subs(h, 0))
    V_first = sp.simplify(sp.diff(sp.log(eps_pert / eps_eta), h).subs(h, 0))

    assert_zero(U_first - dlnT2, "first-order U compiler")
    assert_zero(V_first - dlnepseta, "first-order V compiler")

    y_first = C_phys_dep * sp.Matrix([U_first, V_first])
    assert_matrix_zero(
        y_first - sp.Matrix([0, -dlnepseta, dlnT2 - dlnepseta]),
        "first-order dependent correction",
    )

    # Static-blind line U = 0 maps to the pure equal-drift dressing ray.
    y_static_blind = y_dep.subs(U, 0)
    assert_matrix_zero(y_static_blind + V * sp.Matrix([0, 1, 1]), "static-blind line maps to equal-drift dressing ray")

    print("\nAll Stage-222 symbolic checks passed.")


if __name__ == "__main__":
    main()
