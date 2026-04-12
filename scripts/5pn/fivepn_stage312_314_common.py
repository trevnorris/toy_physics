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


def direct_monomials(
    lambda_W,
    c_etaU,
    gamma,
    K_U,
    K_etaeff,
    K_Weff,
    mu_W,
    T_U,
    L,
    sigma,
    chi0_star,
    deltaU_star,
    epsW_star,
    eps_star,
):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)

    chi0 = sp.simplify(gamma * c_etaU / K_U)
    deltaU = sp.simplify(sp.pi**2 * T_U / (L**2 * K_U))
    epsEta = sp.simplify(c_etaU**2 / (K_U * K_etaeff))
    epsW = sp.simplify(gamma**2 * lambda_W**2 * sigma / (K_U * K_Weff))
    ZW_over_OmW2 = sp.simplify(lambda_W**2 * mu_W / (K_etaeff * K_Weff**2))

    Ctr = sp.simplify(chi0 ** (1 + deltaU_star) * deltaU ** (1 + chi0_star))
    Cnt = sp.simplify(ZW_over_OmW2 * epsW**E_star * deltaU ** (-F_star))

    return {
        "chi0": chi0,
        "deltaU": deltaU,
        "epsEta": epsEta,
        "epsW": epsW,
        "ZW_over_OmW2": ZW_over_OmW2,
        "Ctr": Ctr,
        "Cnt": Cnt,
        "E_star": E_star,
        "F_star": F_star,
    }


def microscopic_slippages(lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1):
    Sigma_chi = sp.simplify(gamma_1 + c_1 - kappa_U)
    Sigma_delta = sp.simplify(tau_1 - kappa_U)
    Sigma_eta = sp.simplify(2 * c_1 - kappa_U - kappa_eta)
    Sigma_Z = sp.simplify(2 * lambda_1 + mu_1 - kappa_eta - 2 * kappa_W)
    Sigma_eps = sp.simplify(2 * gamma_1 + 2 * lambda_1 - kappa_U - kappa_W)
    return {
        "Sigma_chi": Sigma_chi,
        "Sigma_delta": Sigma_delta,
        "Sigma_eta": Sigma_eta,
        "Sigma_Z": Sigma_Z,
        "Sigma_eps": Sigma_eps,
    }


def branch_adapted_slippages(
    chi0_star,
    deltaU_star,
    epsW_star,
    eps_star,
    Sigma_chi,
    Sigma_delta,
    Sigma_eta,
    Sigma_Z,
    Sigma_eps,
):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    Sigma_tr = sp.simplify((1 + chi0_star) * Sigma_delta + (1 + deltaU_star) * Sigma_chi)
    Sigma_nt = sp.simplify(Sigma_Z + E_star * Sigma_eps - F_star * Sigma_delta)
    return {"Sigma_tr": Sigma_tr, "Sigma_nt": Sigma_nt, "Sigma_eta": Sigma_eta, "E_star": E_star, "F_star": F_star}


def compatibility_ledger(
    chi0_star,
    deltaU_star,
    epsW_star,
    eps_star,
    lambda_1,
    c_1,
    gamma_1,
    kappa_U,
    kappa_W,
):
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    tau_1 = sp.simplify(
        kappa_U - (1 + deltaU_star) * (gamma_1 + c_1 - kappa_U) / (1 + chi0_star)
    )
    kappa_eta = sp.simplify(2 * c_1 - kappa_U)
    mu_1 = sp.simplify(
        2 * c_1
        - kappa_U
        + 2 * kappa_W
        - 2 * lambda_1
        - E_star * (2 * gamma_1 + 2 * lambda_1 - kappa_U - kappa_W)
        - F_star * (1 + deltaU_star) * (gamma_1 + c_1 - kappa_U) / (1 + chi0_star)
    )
    return {"tau_1": tau_1, "kappa_eta": kappa_eta, "mu_1": mu_1, "E_star": E_star, "F_star": F_star}
