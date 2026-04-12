from __future__ import annotations

import sympy as sp

from fivepn_stage323_325_common import (
    banner,
    subbanner,
    expect_zero,
    lambda0_constant,
    split_epsilon,
    placement_finite_packet,
    placement_monomial_drifts,
    coherent_branch_observables,
    support_packet,
)
from fivepn_stage319_322_common import (
    microscopic_reduced_drifts,
)
from fivepn_stage312_314_common import (
    direct_monomials,
    branch_adapted_slippages,
)
from fivepn_stage315_317_common import defect_map_from_sigmas


def microscopic_coherent_placement_state(
    lambdaW,
    cetaU,
    gamma,
    KU,
    Ketaeff,
    KWeff,
    muW,
    TU,
    L,
    sigma,
    Lambda0,
    lambdaPhi,
    KPhieff,
):
    """
    Exact microscopic -> actual coherent placement-state map.

    The actual coherent branch variables used in Stages 323–325 are
        (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda, zeta),
    where
        ZhatW = Z_W * Lambda0 / Lambda.
    """
    chi0 = sp.simplify(gamma * cetaU / KU)
    deltaU = sp.simplify(sp.pi**2 * TU / (L**2 * KU))
    ZW = sp.simplify(lambdaW**2 / (Ketaeff * KWeff))
    epsW = sp.simplify(gamma**2 * lambdaW**2 * sigma / (KU * KWeff))
    epsEta = sp.simplify(cetaU**2 / (KU * Ketaeff))
    Lambda = sp.simplify(Lambda0 * KWeff / muW)
    zeta = sp.simplify(lambdaPhi**2 * KWeff / (lambdaW**2 * KPhieff))
    ZhatW = sp.simplify(ZW * Lambda0 / Lambda)
    return {
        "chi0": chi0,
        "deltaU": deltaU,
        "ZW": ZW,
        "epsW": epsW,
        "epsEta": epsEta,
        "Lambda": Lambda,
        "zeta": zeta,
        "ZhatW": ZhatW,
    }



def microscopic_coherent_placement_drifts(
    lambda1,
    c1,
    gamma1,
    kappaU,
    kappaEta,
    kappaW,
    mu1,
    tau1,
    phi1,
    kappaPhi,
):
    """
    Exact microscopic grouped weak-axisymmetric drifts -> actual coherent placement drifts.

    Adds the split variables absent from Stage 319:
        dln Z_W, dln Lambda, dln zeta,
    while keeping the earlier reduced drifts.
    """
    dln_chi0 = sp.simplify(gamma1 + c1 - kappaU)
    dln_deltaU = sp.simplify(tau1 - kappaU)
    dln_ZW = sp.simplify(2 * lambda1 - kappaEta - kappaW)
    dln_epsW = sp.simplify(2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
    dln_epsEta = sp.simplify(2 * c1 - kappaU - kappaEta)
    dln_Lambda = sp.simplify(kappaW - mu1)
    dln_zeta = sp.simplify(2 * phi1 - 2 * lambda1 + kappaW - kappaPhi)
    dln_ZhatW = sp.simplify(dln_ZW - dln_Lambda)
    return {
        "dln_chi0": dln_chi0,
        "dln_deltaU": dln_deltaU,
        "dln_ZW": dln_ZW,
        "dln_epsW": dln_epsW,
        "dln_epsEta": dln_epsEta,
        "dln_Lambda": dln_Lambda,
        "dln_zeta": dln_zeta,
        "dln_ZhatW": dln_ZhatW,
    }



def microscopic_packet_from_placement_drifts(
    drifts: dict[str, sp.Expr],
    chi0_star,
    deltaU_star,
    epsW_star,
    epsEta_star,
):
    eps_star = split_epsilon(epsW_star, deltaU_star)
    sigmas = placement_monomial_drifts(
        drifts["dln_chi0"],
        drifts["dln_deltaU"],
        drifts["dln_ZW"],
        drifts["dln_epsW"],
        drifts["dln_epsEta"],
        drifts["dln_Lambda"],
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )
    defects = defect_map_from_sigmas(
        sigmas["Sigma_tr"],
        sigmas["Sigma_nt"],
        sigmas["Sigma_eta"],
        chi0_star,
        deltaU_star,
        epsEta_star,
    )
    return {
        **sigmas,
        **defects,
        "eps_star": eps_star,
    }



def microscopic_direct_sigmas(
    lambda1,
    c1,
    gamma1,
    kappaU,
    kappaEta,
    kappaW,
    mu1,
    tau1,
    chi0_star,
    deltaU_star,
    epsW_star,
    epsEta_star,
):
    """Direct Stage-313 branch-adapted slippages for comparison."""
    eps_star = split_epsilon(epsW_star, deltaU_star)
    slips = microscopic_reduced_drifts(lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1)
    return branch_adapted_slippages(
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
        slips["dln_chi0"],
        slips["dln_deltaU"],
        slips["dln_epsEta"],
        slips["dln_ZhatW"],
        slips["dln_epsW"],
    )



def microscopic_finite_packet(
    state: dict[str, sp.Expr],
    state_ref: dict[str, sp.Expr],
    chi0_star,
    deltaU_star,
    epsW_star,
):
    eps_star = split_epsilon(epsW_star, deltaU_star)
    return placement_finite_packet(
        state["chi0"],
        state["deltaU"],
        state["ZW"],
        state["epsW"],
        state["epsEta"],
        state["Lambda"],
        state_ref["chi0"],
        state_ref["deltaU"],
        state_ref["ZW"],
        state_ref["epsW"],
        state_ref["epsEta"],
        state_ref["Lambda"],
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )



def direct_microscopic_placement_monomials(
    lambdaW,
    cetaU,
    gamma,
    KU,
    Ketaeff,
    KWeff,
    muW,
    TU,
    L,
    sigma,
    Lambda0,
    chi0_star,
    deltaU_star,
    epsW_star,
):
    """
    Direct microscopic monomials in the actual placement variables.

    This keeps Z_W and Lambda separate rather than collapsing immediately to ZhatW.
    """
    eps_star = split_epsilon(epsW_star, deltaU_star)
    state = microscopic_coherent_placement_state(
        lambdaW,
        cetaU,
        gamma,
        KU,
        Ketaeff,
        KWeff,
        muW,
        TU,
        L,
        sigma,
        Lambda0,
        sp.Symbol("lambda_phi_dummy", positive=True),
        sp.Symbol("K_phi_eff_dummy", positive=True),
    )
    E_star, F_star = direct_monomials(
        lambdaW,
        cetaU,
        gamma,
        KU,
        Ketaeff,
        KWeff,
        muW,
        TU,
        L,
        sigma,
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )["E_star"], direct_monomials(
        lambdaW,
        cetaU,
        gamma,
        KU,
        Ketaeff,
        KWeff,
        muW,
        TU,
        L,
        sigma,
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )["F_star"]
    Ctr = sp.simplify(state["chi0"] ** (1 + deltaU_star) * state["deltaU"] ** (1 + chi0_star))
    Cnt = sp.simplify(
        (state["ZW"] * Lambda0 / state["Lambda"])
        * state["epsW"] ** E_star
        * state["deltaU"] ** (-F_star)
    )
    return {
        "Ctr": Ctr,
        "Cnt": Cnt,
        "epsEta": state["epsEta"],
        "E_star": E_star,
        "F_star": F_star,
    }
