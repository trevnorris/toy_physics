#!/usr/bin/env python3
"""SymPy audit for Moving-Throat PDE Stage 222.

This script verifies the rigid-mouth physical normal form:

1. The rigid-mouth packet is diagonal in the physical logarithmic coordinates
       U = log(T^2/T_ref^2),  V = log(eps_eta/eps_eta_ref)
2. The exact target-ratio formula
       R_target/R_target_ref = ((1 - eps_ref*e^V)/(1 - eps_ref)) * e^{-U}
3. The complementary physical projectors and the two commuting finite legs
4. The Stage-219/221 dependent-plane compiler propagated into the physical chart
       (Delta_T, Delta_Keta, Delta_mu) = (0, -V, U - V)
5. The pure transfer-shape and pure dressing microscopic images
6. The static-only, post-static, and full orbit-restoring correction compilers
7. Propagation of Stage-221 support-blindness through U, V, and the correction packet
8. The Cartesian orbit-lock equivalence inherited from the carried rigid-mouth packet,
   including first-order form
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


def expand_log_simplify(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.factor(sp.expand_log(expr, force=True)))


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

    chi0, deltaU = sp.symbols("chi0 deltaU", positive=True)
    ZW, OmegaW2 = sp.symbols("ZW OmegaW2", positive=True)
    eps = sp.symbols("eps", positive=True)
    B, C = sp.symbols("B C", nonzero=True)
    Rtr, Rtr_ref = sp.symbols("Rtr Rtr_ref", positive=True)

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
    # 3. Stage-219/221 dependent-plane compiler in the physical chart
    # ------------------------------------------------------------------
    # Stage 221 rigid-mouth branch observables.
    T2_stage221 = ZW * (1 + chi0) ** 2 / (OmegaW2 * (1 - eps) ** 2)
    Rtarget_stage221 = (
        Lambda0 * OmegaW2 * (1 - eps_eta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2)
    )
    assert_zero(
        Rtarget_stage221 * T2_stage221 - Lambda0 * (1 - eps_eta),
        "Stage-221 transfer-shape identity",
    )

    q_nt_general = (
        B * sp.log(Rtr / Rtr_ref)
        + sp.log((1 - eps_eta) / (1 - eps_eta_ref))
        - sp.log(Rtarget / Rtarget_ref)
    )
    q_nt_rigid_stage221 = expand_log_simplify(
        q_nt_general.subs(
            {
                Rtr: Rtr_ref,
                Rtarget: Rtarget_stage221,
                Rtarget_ref: Lambda0 * (1 - eps_eta_ref) / T2_ref,
            }
        )
    )
    U_stage221 = sp.log(T2_stage221 / T2_ref)
    V_stage221 = sp.log(eps_eta / eps_eta_ref)

    assert_zero(
        sp.exp(q_nt_rigid_stage221) - T2_stage221 / T2_ref,
        "Stage-221 rigid-mouth packet factorization for q_nt",
    )
    assert_zero(
        expand_log_simplify(q_nt_rigid_stage221 - U_stage221),
        "Stage-221 rigid-mouth identification q_nt = U",
    )
    assert_zero(
        expand_log_simplify(V_stage221 - V_def),
        "Stage-221 rigid-mouth identification q_eta = V",
    )

    # Stage 219 rigid-mouth dependent-plane compiler.
    S_rm_dep = sp.Matrix([[0, 0], [0, -1], [1, -1]])
    C_phys_dep = S_rm_dep * M_phys
    y_dep = C_phys_dep * x_phys
    assert_matrix_zero(C_phys_dep - sp.Matrix([[0, 0], [0, -1], [1, -1]]), "physical-to-microscopic compiler inherited from Stage 219")
    assert_matrix_zero(
        y_dep - sp.Matrix([0, -V, U - V]),
        "physical-to-microscopic dependent compiler",
    )

    x_phys_stage221 = sp.Matrix([U_stage221, V_stage221])
    y_dep_stage221 = S_rm_dep * sp.Matrix([q_nt_rigid_stage221, V_stage221])
    propagated_residual = (y_dep_stage221 - C_phys_dep * x_phys_stage221).applyfunc(
        expand_log_simplify
    )
    assert_matrix_zero(
        propagated_residual,
        "Stage-219/221 dependent compiler propagated into the physical chart",
    )

    # Recover q_nt, q_eta from the Stage-219 dependent-plane equations:
    #   Delta_Keta = -q_eta,   Delta_mu = q_nt - q_eta.
    q_nt_var, q_eta_var = sp.symbols("q_nt_var q_eta_var", real=True)
    yK_var, yMu_var = sp.symbols("yK_var yMu_var", real=True)
    inverse_solution = sp.solve(
        [sp.Eq(yK_var, -q_eta_var), sp.Eq(yMu_var, q_nt_var - q_eta_var)],
        [q_nt_var, q_eta_var],
        dict=True,
    )
    expected_solution = [{q_nt_var: yMu_var - yK_var, q_eta_var: -yK_var}]
    if inverse_solution != expected_solution:
        raise AssertionError(f"Unexpected dependent-plane inverse: {inverse_solution}")

    L_phys_dep = sp.Matrix([[0, -1, 1], [0, -1, 0]])
    assert_matrix_zero(
        L_phys_dep * sp.Matrix([0, yK_var, yMu_var])
        - sp.Matrix([expected_solution[0][q_nt_var], expected_solution[0][q_eta_var]]),
        "left inverse reconstructed from Stage-219 dependent-plane equations",
    )
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
    # 6. Propagation of Stage-221 support-blindness
    # ------------------------------------------------------------------
    # Stage 221 provides support-blind direct observables. Propagate that
    # hypothesis through the physical chart and the Stage-219 dependent packet.
    T2_sb = sp.Function("T2_sb")(zeta, M_mix)
    eps_eta_sb = sp.Function("eps_eta_sb")(zeta, M_mix)
    U_sb = sp.log(T2_sb / T2_ref)
    V_sb = sp.log(eps_eta_sb / eps_eta_ref)
    support_blind_subs = {
        sp.diff(T2_sb, zeta): 0,
        sp.diff(T2_sb, M_mix): 0,
        sp.diff(eps_eta_sb, zeta): 0,
        sp.diff(eps_eta_sb, M_mix): 0,
    }

    assert_zero(sp.diff(T2_stage221, zeta), "Stage-221 branch formula for T^2 is support-blind w.r.t. zeta")
    assert_zero(sp.diff(T2_stage221, M_mix), "Stage-221 branch formula for T^2 is support-blind w.r.t. M_mix")
    assert_zero(
        sp.diff(U_sb, zeta).subs(support_blind_subs),
        "Stage-221 support-blind T^2 propagates to U w.r.t. zeta",
    )
    assert_zero(
        sp.diff(V_sb, zeta).subs(support_blind_subs),
        "Stage-221 support-blind eps_eta propagates to V w.r.t. zeta",
    )
    assert_zero(
        sp.diff(U_sb, M_mix).subs(support_blind_subs),
        "Stage-221 support-blind T^2 propagates to U w.r.t. M_mix",
    )
    assert_zero(
        sp.diff(V_sb, M_mix).subs(support_blind_subs),
        "Stage-221 support-blind eps_eta propagates to V w.r.t. M_mix",
    )

    y_dep_sb = C_phys_dep * sp.Matrix([U_sb, V_sb])
    Delta_y_static_sb = -C_phys_dep * sp.Matrix([U_sb, 0])
    Delta_y_orbit_sb = -y_dep_sb

    assert_matrix_zero(
        y_dep_sb.applyfunc(lambda x: sp.simplify(sp.diff(x, zeta).subs(support_blind_subs))),
        "Stage-221 support-blindness propagates to the dependent correction w.r.t. zeta",
    )
    assert_matrix_zero(
        y_dep_sb.applyfunc(lambda x: sp.simplify(sp.diff(x, M_mix).subs(support_blind_subs))),
        "Stage-221 support-blindness propagates to the dependent correction w.r.t. M_mix",
    )
    assert_matrix_zero(
        Delta_y_static_sb.applyfunc(lambda x: sp.simplify(sp.diff(x, zeta).subs(support_blind_subs))),
        "Stage-221 support-blindness propagates to the static correction w.r.t. zeta",
    )
    assert_matrix_zero(
        Delta_y_orbit_sb.applyfunc(lambda x: sp.simplify(sp.diff(x, zeta).subs(support_blind_subs))),
        "Stage-221 support-blindness propagates to the orbit correction w.r.t. zeta",
    )
    assert_matrix_zero(
        Delta_y_static_sb.applyfunc(lambda x: sp.simplify(sp.diff(x, M_mix).subs(support_blind_subs))),
        "Stage-221 support-blindness propagates to the static correction w.r.t. M_mix",
    )
    assert_matrix_zero(
        Delta_y_orbit_sb.applyfunc(lambda x: sp.simplify(sp.diff(x, M_mix).subs(support_blind_subs))),
        "Stage-221 support-blindness propagates to the orbit correction w.r.t. M_mix",
    )

    # ------------------------------------------------------------------
    # 7. Cartesian orbit-lock equivalence and first-order form
    # ------------------------------------------------------------------
    orbit_lock_solution = sp.solve(
        [sp.Eq(y_dep[1], 0), sp.Eq(y_dep[2], 0)],
        [U, V],
        dict=True,
    )
    if orbit_lock_solution != [{U: 0, V: 0}]:
        raise AssertionError(f"Unexpected orbit-lock solution set: {orbit_lock_solution}")

    assert_matrix_zero(y_dep.subs({U: 0, V: 0}), "U = V = 0 cancels the dependent defect")
    assert_zero(U_def.subs(T2, T2_ref), "T^2 = T_ref^2 implies U = 0")
    assert_zero(V_def.subs(eps_eta, eps_eta_ref), "eps_eta = eps_eta_ref implies V = 0")

    transfer_equilibrium_solution = sp.solve(sp.Eq(T2 / T2_ref, 1), T2, dict=True)
    dressing_equilibrium_solution = sp.solve(sp.Eq(eps_eta / eps_eta_ref, 1), eps_eta, dict=True)
    if transfer_equilibrium_solution != [{T2: T2_ref}]:
        raise AssertionError(f"Unexpected transfer-equilibrium solve: {transfer_equilibrium_solution}")
    if dressing_equilibrium_solution != [{eps_eta: eps_eta_ref}]:
        raise AssertionError(f"Unexpected dressing-equilibrium solve: {dressing_equilibrium_solution}")

    assert_zero(
        ratio_from_UV.subs({U: 0, V: 0}) - 1,
        "U = V = 0 restores the target ratio",
    )

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
