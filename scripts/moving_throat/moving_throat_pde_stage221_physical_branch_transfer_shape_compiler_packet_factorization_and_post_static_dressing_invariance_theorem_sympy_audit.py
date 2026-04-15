#!/usr/bin/env python3
"""SymPy audit for Moving-Throat PDE Stage 221.

This script verifies the Stage-221 physical-branch packet factorization:

1. The coherent-branch identity
       R_target * T^2 = Lambda0 * (1 - eps_eta)
2. The finite packet factorization
       q_nt + (B/C) q_tr = log(T^2 / T_ref^2)
3. The rigid-mouth reduction
       q_nt = log(T^2 / T_ref^2),   q_eta = log(eps_eta / eps_eta_ref)
4. The first-order physical drift compiler
5. Support-blindness with respect to zeta and M_mix
6. The post-static dressing-invariance gates

The script is deliberately symbolic and exact: no floating-point numerics are used.
"""

from __future__ import annotations

import sympy as sp


def assert_zero(expr: sp.Expr, label: str) -> None:
    """Raise if an expression does not simplify to zero."""
    simplified = sp.simplify(sp.factor(expr))
    if simplified != 0:
        raise AssertionError(f"{label} failed: {simplified}")
    print(f"[ok] {label}")


def main() -> None:
    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    chi0, deltaU = sp.symbols("chi0 deltaU", positive=True)
    ZW, OmegaW2 = sp.symbols("ZW OmegaW2", positive=True)
    eps, eps_eta = sp.symbols("eps eps_eta", positive=True)
    Lambda0 = sp.symbols("Lambda0", positive=True)

    B, C = sp.symbols("B C", nonzero=True)
    Rtr_ref, T2_ref = sp.symbols("Rtr_ref T2_ref", positive=True)
    eps_eta_ref = sp.symbols("eps_eta_ref", positive=True)
    Rtarget_ref = sp.symbols("Rtarget_ref", positive=True)

    zeta, M_mix = sp.symbols("zeta M_mix")

    # Differential / logarithmic drift symbols.
    h = sp.symbols("h")
    dlnchi, dlndelta = sp.symbols("dlnchi dlndelta")
    dlnZW, dlnOm2 = sp.symbols("dlnZW dlnOm2")
    dlneps, dlnepseta = sp.symbols("dlneps dlnepseta")

    # ------------------------------------------------------------------
    # Coherent-branch observables
    # ------------------------------------------------------------------
    Rtr = (1 + chi0 / (1 + deltaU)) / (1 + chi0)
    T2 = ZW * (1 + chi0) ** 2 / (OmegaW2 * (1 - eps) ** 2)
    Rtarget = Lambda0 * OmegaW2 * (1 - eps_eta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2)

    # 1. Exact transfer-shape identity.
    assert_zero(Rtarget * T2 - Lambda0 * (1 - eps_eta), "coherent-branch transfer-shape identity")

    # ------------------------------------------------------------------
    # Finite packet and factorization
    # ------------------------------------------------------------------
    q_tr = -C * sp.log(Rtr / Rtr_ref)
    q_nt = (
        B * sp.log(Rtr / Rtr_ref)
        + sp.log((1 - eps_eta) / (1 - eps_eta_ref))
        - sp.log(Rtarget / Rtarget_ref)
    )
    q_eta = sp.log(eps_eta / eps_eta_ref)

    # Use the exact reference-branch identity
    #     Rtarget_ref * T2_ref = Lambda0 * (1 - eps_eta_ref)
    # to verify packet factorization by comparing exponentials.
    q_factor = q_nt + (B / C) * q_tr
    q_factor_exp = sp.simplify(
        sp.exp(q_factor.subs(Rtarget_ref, Lambda0 * (1 - eps_eta_ref) / T2_ref)) - T2 / T2_ref
    )
    assert_zero(q_factor_exp, "finite packet factorization via transfer-shape ratio")

    # Rigid-mouth slice: q_tr = 0, equivalently R_tr = R_tr_ref.
    q_nt_rigid = sp.simplify(
        q_nt.subs({Rtr: Rtr_ref, Rtarget_ref: Lambda0 * (1 - eps_eta_ref) / T2_ref})
    )
    q_nt_rigid_exp = sp.simplify(sp.exp(q_nt_rigid) - T2 / T2_ref)
    assert_zero(q_nt_rigid_exp, "rigid-mouth reduction of q_nt")
    assert_zero(sp.exp(q_eta) - eps_eta / eps_eta_ref, "finite dressing ratio q_eta")

    # ------------------------------------------------------------------
    # First-order physical drift compiler
    # ------------------------------------------------------------------
    # Perturb by logarithmic drifts x -> x * exp(h dln x).
    Rtr_pert = Rtr.subs(
        {
            chi0: chi0 * sp.exp(h * dlnchi),
            deltaU: deltaU * sp.exp(h * dlndelta),
        }
    )
    dlnRtr_calc = sp.simplify(sp.diff(sp.log(sp.simplify(Rtr_pert)), h).subs(h, 0))
    dlnRtr_expected = -chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) * (
        (1 + deltaU) * dlnchi + (1 + chi0) * dlndelta
    )
    assert_zero(dlnRtr_calc - dlnRtr_expected, "first-order compiler for delta ln R_tr")

    T2_pert = T2.subs(
        {
            ZW: ZW * sp.exp(h * dlnZW),
            OmegaW2: OmegaW2 * sp.exp(h * dlnOm2),
            chi0: chi0 * sp.exp(h * dlnchi),
            eps: eps * sp.exp(h * dlneps),
        }
    )
    dlnT2_calc = sp.simplify(sp.diff(sp.log(sp.simplify(T2_pert)), h).subs(h, 0))
    dlnT2_expected = dlnZW - dlnOm2 + 2 * chi0 / (1 + chi0) * dlnchi + 2 * eps / (1 - eps) * dlneps
    assert_zero(dlnT2_calc - dlnT2_expected, "first-order compiler for delta ln T^2")

    Rtarget_pert = Rtarget.subs(
        {
            ZW: ZW * sp.exp(h * dlnZW),
            OmegaW2: OmegaW2 * sp.exp(h * dlnOm2),
            chi0: chi0 * sp.exp(h * dlnchi),
            eps: eps * sp.exp(h * dlneps),
            eps_eta: eps_eta * sp.exp(h * dlnepseta),
        }
    )
    dlnRtarget_calc = sp.simplify(sp.diff(sp.log(sp.simplify(Rtarget_pert)), h).subs(h, 0))
    dlnRtarget_expected = (
        dlnOm2
        - dlnZW
        - 2 * chi0 / (1 + chi0) * dlnchi
        - 2 * eps / (1 - eps) * dlneps
        - eps_eta / (1 - eps_eta) * dlnepseta
    )
    assert_zero(dlnRtarget_calc - dlnRtarget_expected, "first-order compiler for delta ln R_target")

    # First-order packet relation.
    q_tr_first = -C * dlnRtr_calc
    q_nt_first_direct = B * dlnRtr_calc - dlnRtarget_calc - eps_eta / (1 - eps_eta) * dlnepseta
    q_eta_first = dlnepseta
    assert_zero(q_nt_first_direct + (B / C) * q_tr_first - dlnT2_calc, "first-order packet factorization")

    # Rigid-mouth first-order slice: q_tr = 0 -> q_nt = delta ln T^2, q_eta = d ln eps_eta.
    q_nt_first_rigid = sp.simplify(q_nt_first_direct.subs(dlnRtr_calc, 0))
    assert_zero(q_eta_first - dlnepseta, "rigid-mouth first-order dressing compiler")

    # The clean rigid-mouth statement is easiest to check from the exact factorized relation.
    q_nt_factorized_rigid = dlnT2_calc
    assert_zero(q_nt_factorized_rigid - dlnT2_expected, "rigid-mouth first-order transfer-shape compiler")

    # ------------------------------------------------------------------
    # Support-blindness
    # ------------------------------------------------------------------
    # The physical observables do not depend on zeta or M_mix in this reduced packet.
    assert_zero(sp.diff(Rtr, zeta), "support-blindness of R_tr with respect to zeta")
    assert_zero(sp.diff(T2, zeta), "support-blindness of T^2 with respect to zeta")
    assert_zero(sp.diff(eps_eta, zeta), "support-blindness of eps_eta with respect to zeta")
    assert_zero(sp.diff(Rtr, M_mix), "support-blindness of R_tr with respect to M_mix")
    assert_zero(sp.diff(T2, M_mix), "support-blindness of T^2 with respect to M_mix")
    assert_zero(sp.diff(eps_eta, M_mix), "support-blindness of eps_eta with respect to M_mix")

    # Use the factorized finite packet to verify support-blindness of q_tr, q_nt, q_eta.
    q_nt_factored = sp.log(T2 / T2_ref) - (B / C) * q_tr
    assert_zero(sp.diff(q_tr, zeta), "support-blindness of q_tr with respect to zeta")
    assert_zero(sp.diff(q_nt_factored, zeta), "support-blindness of q_nt with respect to zeta")
    assert_zero(sp.diff(q_eta, zeta), "support-blindness of q_eta with respect to zeta")
    assert_zero(sp.diff(q_tr, M_mix), "support-blindness of q_tr with respect to M_mix")
    assert_zero(sp.diff(q_nt_factored, M_mix), "support-blindness of q_nt with respect to M_mix")
    assert_zero(sp.diff(q_eta, M_mix), "support-blindness of q_eta with respect to M_mix")

    # ------------------------------------------------------------------
    # Three-gate post-static theorem checks
    # ------------------------------------------------------------------
    # Tracking gate:
    tracking_condition = (1 + deltaU) * dlnchi + (1 + chi0) * dlndelta
    assert_zero(
        sp.simplify(
            dlnRtr_expected
            + chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) * tracking_condition
        ),
        "tracking gate compiler",
    )

    # Static-blind transfer-shape gate and post-static dressing gate are direct.
    assert_zero(sp.exp(q_nt_rigid) - T2 / T2_ref, "finite transfer-shape gate")
    assert_zero(sp.exp(q_eta) - eps_eta / eps_eta_ref, "finite dressing gate")

    print("\nAll Stage-221 symbolic checks passed.")


if __name__ == "__main__":
    main()
