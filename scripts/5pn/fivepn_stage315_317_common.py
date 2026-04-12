from __future__ import annotations

import sympy as sp


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


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


def monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star):
    E_star = sp.simplify(
        2 * epsW_star * (11 + 9 * deltaU_star)
        / (11 * (1 - eps_star) * (1 + deltaU_star))
    )
    F_star = sp.simplify(
        2 * chi0_star / (1 + deltaU_star)
        + 4 * epsW_star * deltaU_star / (11 * (1 - eps_star) * (1 + deltaU_star) ** 2)
    )
    return E_star, F_star


def continuum_monomials(chi0, deltaU, ZhatW, epsW, epsEta, chi0_star, deltaU_star, epsW_star, eps_star):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    Ctr = sp.simplify(chi0 ** (1 + deltaU_star) * deltaU ** (1 + chi0_star))
    Cnt = sp.simplify(ZhatW * epsW ** E_star * deltaU ** (-F_star))
    return {"Ctr": Ctr, "Cnt": Cnt, "epsEta": sp.simplify(epsEta), "E_star": E_star, "F_star": F_star}


def continuum_monomial_quotients(
    chi0,
    deltaU,
    ZhatW,
    epsW,
    epsEta,
    chi0_ref,
    deltaU_ref,
    ZhatW_ref,
    epsW_ref,
    epsEta_ref,
    chi0_star,
    deltaU_star,
    epsW_star,
    eps_star,
):
    mono = continuum_monomials(chi0, deltaU, ZhatW, epsW, epsEta, chi0_star, deltaU_star, epsW_star, eps_star)
    mono_ref = continuum_monomials(chi0_ref, deltaU_ref, ZhatW_ref, epsW_ref, epsEta_ref, chi0_star, deltaU_star, epsW_star, eps_star)
    qtr = sp.simplify(sp.log(mono["Ctr"] / mono_ref["Ctr"]))
    qnt = sp.simplify(sp.log(mono["Cnt"] / mono_ref["Cnt"]))
    qeta = sp.simplify(sp.log(mono["epsEta"] / mono_ref["epsEta"]))
    return {"qtr": qtr, "qnt": qnt, "qeta": qeta, "E_star": mono["E_star"], "F_star": mono["F_star"]}


def continuum_monomial_drifts(
    dln_chi0,
    dln_deltaU,
    dln_ZhatW,
    dln_epsW,
    dln_epsEta,
    chi0_star,
    deltaU_star,
    epsW_star,
    eps_star,
):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    Sigma_tr = sp.simplify((1 + deltaU_star) * dln_chi0 + (1 + chi0_star) * dln_deltaU)
    Sigma_nt = sp.simplify(dln_ZhatW + E_star * dln_epsW - F_star * dln_deltaU)
    Sigma_eta = sp.simplify(dln_epsEta)
    return {"Sigma_tr": Sigma_tr, "Sigma_nt": Sigma_nt, "Sigma_eta": Sigma_eta, "E_star": E_star, "F_star": F_star}


def defect_map_from_sigmas(Sigma_tr, Sigma_nt, Sigma_eta, chi0_star, deltaU_star, epsEta_star):
    Ctr = sp.simplify(
        chi0_star * deltaU_star / ((1 + chi0_star) * (1 + deltaU_star) * (1 + chi0_star + deltaU_star))
    )
    Atr = sp.simplify(2 * chi0_star / ((1 + chi0_star) * (1 + deltaU_star)))
    Theta1 = sp.simplify(-Ctr * Sigma_tr)
    Xi1 = sp.simplify(Atr * Sigma_tr + Sigma_nt)
    R1 = sp.simplify(-Xi1 - epsEta_star * Sigma_eta / (1 - epsEta_star))
    return {"Theta1": Theta1, "Xi1": Xi1, "R1": R1, "Ctr": Ctr, "Atr": Atr}


def required_codrifts_for_orbit_lock(
    dln_chi0,
    dln_epsW,
    chi0_star,
    deltaU_star,
    epsW_star,
    eps_star,
):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    dln_deltaU = sp.simplify(-(1 + deltaU_star) * dln_chi0 / (1 + chi0_star))
    dln_ZhatW = sp.simplify(-E_star * dln_epsW + F_star * dln_deltaU)
    dln_epsEta = sp.Integer(0)
    return {
        "dln_deltaU": dln_deltaU,
        "dln_ZhatW": dln_ZhatW,
        "dln_epsEta": dln_epsEta,
        "E_star": E_star,
        "F_star": F_star,
    }


def reduced_zero_defect_basis(chi0_star, deltaU_star, epsW_star, eps_star):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    A = sp.simplify((1 + deltaU_star) / (1 + chi0_star))
    v_chi = sp.Matrix([
        1,
        -A,
        -F_star * A,
        0,
        0,
    ])
    v_epsW = sp.Matrix([
        0,
        0,
        -E_star,
        1,
        0,
    ])
    return {"v_chi": v_chi, "v_epsW": v_epsW, "E_star": E_star, "F_star": F_star, "A": A}
