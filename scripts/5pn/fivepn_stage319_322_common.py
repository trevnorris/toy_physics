from __future__ import annotations

import sympy as sp

from fivepn_stage315_317_common import (
    banner,
    subbanner,
    expect_zero,
    monomial_exponents,
    continuum_monomial_quotients,
    continuum_monomial_drifts,
    defect_map_from_sigmas,
)


def microscopic_reduced_drifts(
    lambda1,
    c1,
    gamma1,
    kappaU,
    kappaEta,
    kappaW,
    mu1,
    tau1,
):
    dln_chi0 = sp.simplify(gamma1 + c1 - kappaU)
    dln_deltaU = sp.simplify(tau1 - kappaU)
    dln_ZhatW = sp.simplify(2 * lambda1 + mu1 - kappaEta - 2 * kappaW)
    dln_epsW = sp.simplify(2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
    dln_epsEta = sp.simplify(2 * c1 - kappaU - kappaEta)
    return {
        "dln_chi0": dln_chi0,
        "dln_deltaU": dln_deltaU,
        "dln_ZhatW": dln_ZhatW,
        "dln_epsW": dln_epsW,
        "dln_epsEta": dln_epsEta,
    }


def microscopic_state_to_reduced(
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
):
    chi0 = sp.simplify(gamma * cetaU / KU)
    deltaU = sp.simplify(sp.pi**2 * TU / (L**2 * KU))
    ZhatW = sp.simplify(lambdaW**2 * muW / (Ketaeff * KWeff**2))
    epsW = sp.simplify(gamma**2 * lambdaW**2 * sigma / (KU * KWeff))
    epsEta = sp.simplify(cetaU**2 / (KU * Ketaeff))
    return {
        "chi0": chi0,
        "deltaU": deltaU,
        "ZhatW": ZhatW,
        "epsW": epsW,
        "epsEta": epsEta,
    }


def physical_branch_observables_from_microscopic(
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
):
    red = microscopic_state_to_reduced(lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma)
    chi0 = red["chi0"]
    deltaU = red["deltaU"]
    ZhatW = red["ZhatW"]
    epsW = red["epsW"]
    epsEta = red["epsEta"]

    eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
    Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    Rtarget = sp.simplify(Lambda0 * (1 - epsEta) * (1 - eps) ** 2 / (ZhatW * (1 + chi0) ** 2))
    return {
        **red,
        "eps": eps,
        "Rtr": Rtr,
        "Rtarget": Rtarget,
    }


def direct_microscopic_monomials(
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
):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    Ctr = sp.simplify(
        (gamma * cetaU / KU) ** (1 + deltaU_star)
        * (sp.pi**2 * TU / (L**2 * KU)) ** (1 + chi0_star)
    )
    Cnt = sp.simplify(
        (lambdaW**2 * muW / (Ketaeff * KWeff**2))
        * (gamma**2 * lambdaW**2 * sigma / (KU * KWeff)) ** E_star
        * (sp.pi**2 * TU / (L**2 * KU)) ** (-F_star)
    )
    epsEta = sp.simplify(cetaU**2 / (KU * Ketaeff))
    return {"Ctr": Ctr, "Cnt": Cnt, "epsEta": epsEta, "E_star": E_star, "F_star": F_star}


def coherent_support_ratio(lambdaPhi, KPhieff, lambdaW, KWeff):
    return sp.simplify(lambdaPhi**2 * KWeff / (lambdaW**2 * KPhieff))


def coherent_support_enhancement(zeta, eps):
    return sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
