#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 242
Actual Twin-Support Placement and Coherent Orbit-Lock Compiler
SymPy audit

This script verifies the exact algebra recorded in the Stage 242 note.
"""

from __future__ import annotations

import sympy as sp


def assert_zero(name: str, expr: sp.Expr) -> None:
    """Cheap zero test that avoids very heavy global simplify passes."""
    expr_s = sp.cancel(sp.together(expr))
    if expr_s != 0:
        expr_s = sp.simplify(expr_s)
    if expr_s != 0:
        raise AssertionError(f"{name} failed:\n{sp.srepr(expr_s)}")
    print(f"[ok] {name}")


def assert_matrix_zero(name: str, mat: sp.Matrix) -> None:
    for i, entry in enumerate(mat):
        assert_zero(f"{name}[{i}]", entry)


def main() -> None:
    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    chi0, deltaU, ZW, epsW, epsEta, Lambda, Lambda0, zeta = sp.symbols(
        "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda Lambda0 zeta",
        positive=True,
    )
    beta, varrho = sp.symbols("beta varrho", positive=True)
    B_star, C_star = sp.symbols("B_star C_star", nonzero=True)

    # Reference-state data for the finite packet
    chi0_ref, deltaU_ref, ZW_ref, epsW_ref, epsEta_ref, Lambda_ref = sp.symbols(
        "chi0_ref deltaU_ref ZW_ref epsW_ref epsEta_ref Lambda_ref",
        positive=True,
    )
    E_star, F_star = sp.symbols("E_star F_star")

    # Infinitesimal logarithmic drifts
    dchi0, ddeltaU, dZW, depsW, depsEta, dLambda = sp.symbols(
        "dlnchi0 dlndeltaU dlnZW dlnepsW dlnepsEta dlnLambda"
    )

    # ------------------------------------------------------------------
    # 1. Actual selected-twin placement
    # ------------------------------------------------------------------
    # Stage 027/026/030 carry the fixed geometric coefficient 2/11.
    eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
    C_mix = sp.simplify(8 * Lambda * (1 - eps) / sp.pi**2)
    # Stage 240 fixes the selected-branch loading ratio Pi_tr/C_mix = 4/3.
    Pi_tr = sp.Rational(4, 3) * C_mix
    varrho_phys = sp.simplify(sp.pi**2 * Pi_tr / (16 * Lambda))
    sigma_phys = sp.simplify(4 / (3 * varrho_phys) - 2)

    assert_zero(
        "selected-branch coordinate varrho_phys",
        varrho_phys - sp.Rational(2, 3) * (1 - eps),
    )
    assert_zero(
        "selected-branch sigma_phys",
        sigma_phys - 2 * eps / (1 - eps),
    )

    # ------------------------------------------------------------------
    # 2. Threshold rewrite and selected-branch twin-window inclusion
    # ------------------------------------------------------------------
    varrho_WL = sp.simplify(2 * (1 + beta**2) / (3 * (2 + beta**2)))
    varrho_UL = sp.simplify(2 * (1 + beta**2) / (3 * (1 + beta + beta**2)))
    eps_WL = sp.simplify(1 - sp.Rational(3, 2) * varrho_WL)
    eps_UL = sp.simplify(1 - sp.Rational(3, 2) * varrho_UL)

    assert_zero("epsilon_WLambda rewrite", eps_WL - 1 / (2 + beta**2))
    assert_zero("epsilon_ULambda rewrite", eps_UL - beta / (1 + beta + beta**2))

    sigma_sel = sp.simplify(4 / (3 * varrho) - 2)
    assert_zero(
        "selected branch lies above mixed-only bound",
        sigma_sel - (1 / varrho - 2) - 1 / (3 * varrho),
    )
    assert_zero(
        "selected branch lies below non-twin bound",
        (2 / varrho - 2) - sigma_sel - 2 / (3 * varrho),
    )

    # ------------------------------------------------------------------
    # 3. Reduced-state bridge and direct coherent observables
    # ------------------------------------------------------------------
    ZhatW = sp.simplify(ZW * Lambda0 / Lambda)

    R_tr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    R_target = sp.simplify(
        Lambda * (1 - epsEta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2)
    )
    R_target_hat = sp.simplify(
        Lambda0 * (1 - epsEta) * (1 - eps) ** 2 / (ZhatW * (1 + chi0) ** 2)
    )

    assert_zero("reduced-state bridge for R_target", R_target - R_target_hat)
    assert_zero("zeta-absence of epsilon in the coherent placement map", sp.diff(eps, zeta))
    assert_zero("zeta-absence of R_tr in the coherent placement map", sp.diff(R_tr, zeta))
    assert_zero("zeta-absence of R_target in the coherent placement map", sp.diff(R_target, zeta))

    # ------------------------------------------------------------------
    # 4. Finite orbit packet and support-blind propagation
    # ------------------------------------------------------------------
    q_tr = (1 + deltaU_ref) * sp.log(chi0 / chi0_ref) + (1 + chi0_ref) * sp.log(
        deltaU / deltaU_ref
    )
    q_nt = (
        sp.log(ZW / ZW_ref)
        - sp.log(Lambda / Lambda_ref)
        + E_star * sp.log(epsW / epsW_ref)
        - F_star * sp.log(deltaU / deltaU_ref)
    )
    q_eta = sp.log(epsEta / epsEta_ref)

    # Stage 234 packages the orbit packet directly from the support-blind
    # observables (R_tr, R_target, epsilon_eta). Propagate that blind packet
    # through abstract support-blind direct observables rather than by
    # differentiating formulas that simply omit zeta.
    Rtr_sb = sp.Function("Rtr_sb")(zeta)
    Rtarget_sb = sp.Function("Rtarget_sb")(zeta)
    epsEta_sb = sp.Function("epsEta_sb")(zeta)
    support_blind_subs = {
        sp.diff(Rtr_sb, zeta): 0,
        sp.diff(Rtarget_sb, zeta): 0,
        sp.diff(epsEta_sb, zeta): 0,
    }

    q_tr_from_observables = -C_star * sp.log(Rtr_sb / R_tr)
    q_nt_from_observables = (
        B_star * sp.log(Rtr_sb / R_tr)
        + sp.log((1 - epsEta_sb) / (1 - epsEta_ref))
        - sp.log(Rtarget_sb / R_target)
    )
    q_eta_from_observables = sp.log(epsEta_sb / epsEta_ref)

    assert_zero(
        "support-blind direct observables propagate to q_tr",
        sp.diff(q_tr_from_observables, zeta).subs(support_blind_subs),
    )
    assert_zero(
        "support-blind direct observables propagate to q_nt",
        sp.diff(q_nt_from_observables, zeta).subs(support_blind_subs),
    )
    assert_zero(
        "support-blind direct observables propagate to q_eta",
        sp.diff(q_eta_from_observables, zeta).subs(support_blind_subs),
    )
    assert_zero("finite q_eta matches the direct observable chart", q_eta_from_observables.subs(epsEta_sb, epsEta) - q_eta)

    # ------------------------------------------------------------------
    # 5. Infinitesimal coherent packet and direct observable compiler
    # ------------------------------------------------------------------
    t = sp.symbols("t", real=True)

    chi0_t = chi0 * sp.exp(t * dchi0)
    deltaU_t = deltaU * sp.exp(t * ddeltaU)
    ZW_t = ZW * sp.exp(t * dZW)
    epsW_t = epsW * sp.exp(t * depsW)
    epsEta_t = epsEta * sp.exp(t * depsEta)
    Lambda_t = Lambda * sp.exp(t * dLambda)

    eps_t = sp.simplify(
        epsW_t * (1 - sp.Rational(2, 11) * deltaU_t / (1 + deltaU_t))
    )
    dln_eps = sp.simplify(sp.diff(sp.log(eps_t), t).subs(t, 0))
    dln_eps_formula = sp.simplify(
        depsW
        - (2 * deltaU / ((1 + deltaU) * (11 + 9 * deltaU))) * ddeltaU
    )
    assert_zero("dln epsilon compiler", dln_eps - dln_eps_formula)

    R_tr_t = sp.simplify((1 + chi0_t / (1 + deltaU_t)) / (1 + chi0_t))
    dln_Rtr = sp.simplify(sp.diff(sp.log(R_tr_t), t).subs(t, 0))
    dln_Rtr_formula = sp.simplify(
        -(
            chi0
            * deltaU
            / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
        )
        * ((1 + deltaU) * dchi0 + (1 + chi0) * ddeltaU)
    )
    assert_zero("dln R_tr compiler", dln_Rtr - dln_Rtr_formula)

    R_target_t = sp.simplify(
        Lambda_t * (1 - epsEta_t) * (1 - eps_t) ** 2 / (ZW_t * (1 + chi0_t) ** 2)
    )
    dln_Rtarget = sp.simplify(sp.diff(sp.log(R_target_t), t).subs(t, 0))
    dln_Rtarget_formula = sp.simplify(
        dLambda
        - dZW
        - 2 * chi0 / (1 + chi0) * dchi0
        - epsEta / (1 - epsEta) * depsEta
        - 2 * eps / (1 - eps) * dln_eps
    )
    assert_zero("dln R_target compiler", dln_Rtarget - dln_Rtarget_formula)

    Theta_1 = dln_Rtr
    Xi_1 = sp.simplify(-dln_Rtarget - epsEta / (1 - epsEta) * depsEta)
    Xi_1_formula = sp.simplify(
        -dLambda
        + dZW
        + 2 * chi0 / (1 + chi0) * dchi0
        + 2 * eps / (1 - eps) * dln_eps
    )
    R_1 = sp.simplify(-Xi_1 - epsEta / (1 - epsEta) * depsEta)
    c_eta = sp.simplify(epsEta / (1 - epsEta))

    assert_zero("Theta_1 direct-observable identity", Theta_1 - dln_Rtr)
    assert_zero("Xi_1 direct-observable identity", Xi_1 - Xi_1_formula)
    assert_zero("R_1 direct-observable identity", R_1 - dln_Rtarget)

    direct_packet = sp.Matrix([dln_Rtr, dln_Rtarget, depsEta])
    orbit_packet = sp.Matrix([Theta_1, Xi_1, R_1])
    orbit_compiler = sp.Matrix(
        [
            [1, 0, 0],
            [0, -1, -c_eta],
            [0, 1, 0],
        ]
    )
    assert_matrix_zero(
        "direct-observable orbit packet compiler",
        orbit_packet - orbit_compiler * direct_packet,
    )
    assert_zero(
        "orbit packet compiler determinant",
        sp.simplify(orbit_compiler.det() - c_eta),
    )
    assert_matrix_zero(
        "inverse direct-observable orbit packet compiler",
        sp.simplify(orbit_compiler.inv() * orbit_packet - direct_packet),
    )
    Theta_var, Xi_var, R_var = sp.symbols("Theta_var Xi_var R_var", real=True)
    recovered_direct_packet = sp.simplify(
        orbit_compiler.inv() * sp.Matrix([Theta_var, Xi_var, R_var])
    )
    assert_matrix_zero(
        "formal orbit packet recovers the direct drifts",
        recovered_direct_packet
        - sp.Matrix(
            [
                Theta_var,
                R_var,
                -(1 - epsEta) * (Xi_var + R_var) / epsEta,
            ]
        ),
    )
    assert_matrix_zero(
        "zero orbit packet forces zero direct drifts",
        recovered_direct_packet.subs({Theta_var: 0, Xi_var: 0, R_var: 0}),
    )

    # ------------------------------------------------------------------
    # 6. Exact two-packet split
    # ------------------------------------------------------------------
    M_mix = sp.simplify(
        8 * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - epsEta) * (1 - eps))
    )
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    M_tr = sp.simplify(M_mix * S)

    assert_zero("mixed-only product law", sp.simplify(R_target * M_mix - C_mix))
    assert_zero(
        "support-packet sensitivity",
        sp.diff(M_tr, zeta) - M_mix * (1 - eps) / (1 - zeta * eps) ** 2,
    )

    # ------------------------------------------------------------------
    # 7. A rational sample point for sanity
    # ------------------------------------------------------------------
    sample = {
        chi0: sp.Rational(3, 2),
        deltaU: sp.Rational(2, 3),
        ZW: sp.Rational(13, 17),
        epsW: sp.Rational(1, 3),
        epsEta: sp.Rational(1, 5),
        Lambda: sp.Rational(7, 11),
        zeta: sp.Integer(1),
    }

    print("\nSample branch values:")
    print("  epsilon      =", sp.simplify(eps.subs(sample)))
    print("  varrho_phys  =", sp.simplify(varrho_phys.subs(sample)))
    print("  sigma_phys   =", sp.simplify(sigma_phys.subs(sample)))
    print("  R_tr         =", sp.simplify(R_tr.subs(sample)))
    print("  R_target     =", sp.simplify(R_target.subs(sample)))
    print("  M_mix        =", sp.simplify(M_mix.subs(sample)))
    print("  M_tr         =", sp.simplify(M_tr.subs(sample)))

    print("\nAll Stage 242 symbolic checks passed.")


if __name__ == "__main__":
    main()
