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


def lambda0_constant(G, cs, a, c):
    return sp.simplify(27 * sp.pi**2 * G * cs**5 / (20 * a**5 * c**5))


def split_epsilon(epsW, deltaU):
    return sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))


def placement_to_reduced_state(ZW, Lambda, Lambda0):
    return {
        "ZhatW": sp.simplify(ZW * Lambda0 / Lambda),
    }


def coherent_branch_observables(chi0, deltaU, ZW, epsW, epsEta, Lambda):
    eps = split_epsilon(epsW, deltaU)
    Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    Rtarget = sp.simplify(Lambda * (1 - epsEta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))
    return {"eps": eps, "Rtr": Rtr, "Rtarget": Rtarget}


def placement_finite_packet(
    chi0,
    deltaU,
    ZW,
    epsW,
    epsEta,
    Lambda,
    chi0_ref,
    deltaU_ref,
    ZW_ref,
    epsW_ref,
    epsEta_ref,
    Lambda_ref,
    chi0_star,
    deltaU_star,
    epsW_star,
    eps_star,
):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    qtr = sp.simplify(
        (1 + deltaU_star) * sp.log(chi0 / chi0_ref)
        + (1 + chi0_star) * sp.log(deltaU / deltaU_ref)
    )
    qnt = sp.simplify(
        sp.log(ZW / ZW_ref)
        - sp.log(Lambda / Lambda_ref)
        + E_star * sp.log(epsW / epsW_ref)
        - F_star * sp.log(deltaU / deltaU_ref)
    )
    qeta = sp.simplify(sp.log(epsEta / epsEta_ref))
    return {"qtr": qtr, "qnt": qnt, "qeta": qeta, "E_star": E_star, "F_star": F_star}


def placement_monomial_drifts(
    dln_chi0,
    dln_deltaU,
    dln_ZW,
    dln_epsW,
    dln_epsEta,
    dln_Lambda,
    chi0_star,
    deltaU_star,
    epsW_star,
    eps_star,
):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    Sigma_tr = sp.simplify((1 + deltaU_star) * dln_chi0 + (1 + chi0_star) * dln_deltaU)
    Sigma_nt = sp.simplify(dln_ZW - dln_Lambda + E_star * dln_epsW - F_star * dln_deltaU)
    Sigma_eta = sp.simplify(dln_epsEta)
    return {"Sigma_tr": Sigma_tr, "Sigma_nt": Sigma_nt, "Sigma_eta": Sigma_eta, "E_star": E_star, "F_star": F_star}


def log_drift(expr, var_to_dln: dict[sp.Symbol, sp.Symbol]):
    expr = sp.simplify(expr)
    total = 0
    for var, dln in var_to_dln.items():
        total += sp.simplify(var / expr * sp.diff(expr, var) * dln)
    return sp.simplify(sp.expand(total))


def support_packet(chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta):
    eps = split_epsilon(epsW, deltaU)
    Mmix = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - epsEta) * (1 - eps)))
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    Mtr = sp.simplify(Mmix * S)
    Rtarget = sp.simplify(Lambda * (1 - epsEta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))
    product = sp.simplify(Mmix * Rtarget)
    return {
        "eps": eps,
        "Mmix": Mmix,
        "S": S,
        "Mtr": Mtr,
        "Rtarget": Rtarget,
        "product": product,
    }
