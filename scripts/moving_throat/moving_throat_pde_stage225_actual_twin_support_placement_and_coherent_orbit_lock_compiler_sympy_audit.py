#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 225
Actual Twin-Support Placement and Coherent Orbit-Lock Compiler
SymPy audit

This script verifies the exact algebra recorded in the Stage-225 note.
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


def main() -> None:
    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    chi0, deltaU, ZW, epsW, epsEta, Lambda, Lambda0, zeta = sp.symbols(
        "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda Lambda0 zeta",
        positive=True,
    )
    beta, varrho = sp.symbols("beta varrho", positive=True)

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
    eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
    C_mix = sp.simplify(8 * Lambda * (1 - eps) / sp.pi**2)
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
    assert_zero("support-blind epsilon", sp.diff(eps, zeta))
    assert_zero("support-blind R_tr", sp.diff(R_tr, zeta))
    assert_zero("support-blind R_target", sp.diff(R_target, zeta))

    # ------------------------------------------------------------------
    # 4. Finite orbit packet
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

    assert_zero("support-blind q_tr", sp.diff(q_tr, zeta))
    assert_zero("support-blind q_nt", sp.diff(q_nt, zeta))
    assert_zero("support-blind q_eta", sp.diff(q_eta, zeta))

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

    Xi_1 = sp.simplify(-dln_Rtarget - epsEta / (1 - epsEta) * depsEta)
    Xi_1_formula = sp.simplify(
        -dLambda
        + dZW
        + 2 * chi0 / (1 + chi0) * dchi0
        + 2 * eps / (1 - eps) * dln_eps
    )
    assert_zero("Xi_1 direct-observable identity", Xi_1 - Xi_1_formula)

    # Show the exact logical implication used by the orbit-lock theorem:
    # Xi_1 = 0 and R_1 = dln_Rtarget = 0 force dln epsEta = 0.
    reconstructed_depsEta = sp.simplify(-(1 - epsEta) * (Xi_1 + dln_Rtarget) / epsEta)
    assert_zero("orbit-lock eta reconstruction", reconstructed_depsEta - depsEta)

    # ------------------------------------------------------------------
    # 6. Exact two-packet split
    # ------------------------------------------------------------------
    M_mix = sp.simplify(
        8 * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - epsEta) * (1 - eps))
    )
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    M_tr = sp.simplify(M_mix * S)

    assert_zero("support-blind q_tr/zeta", sp.diff(q_tr, zeta))
    assert_zero("support-blind q_nt/zeta", sp.diff(q_nt, zeta))
    assert_zero("support-blind q_eta/zeta", sp.diff(q_eta, zeta))
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

    print("\nAll Stage-225 symbolic checks passed.")


if __name__ == "__main__":
    main()
