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


def monomial_exponents(chi0s, deltaUs, epsWs, epss):
    E_star = sp.simplify(2 * epsWs * (11 + 9 * deltaUs) / (11 * (1 - epss) * (1 + deltaUs)))
    F_star = sp.simplify(2 * chi0s / (1 + deltaUs) + 4 * epsWs * deltaUs / (11 * (1 - epss) * (1 + deltaUs) ** 2))
    return E_star, F_star


def microscopic_monomials(lamW, cetaU, gam, KU, Keta, KW, muW, TU, L, sigma, chi0s, deltaUs, epsWs, epss):
    chi0 = sp.simplify(gam * cetaU / KU)
    deltaU = sp.simplify(sp.pi**2 * TU / (L**2 * KU))
    epsEta = sp.simplify(cetaU**2 / (KU * Keta))
    epsW = sp.simplify(gam**2 * lamW**2 * sigma / (KU * KW))
    ZW_over_OmW2 = sp.simplify(lamW**2 * muW / (Keta * KW**2))
    E_star, F_star = monomial_exponents(chi0s, deltaUs, epsWs, epss)
    Ctr = sp.simplify(chi0 ** (1 + deltaUs) * deltaU ** (1 + chi0s))
    Cnt = sp.simplify(ZW_over_OmW2 * epsW ** E_star * deltaU ** (-F_star))
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


def monomial_drift_matrix(chi0s, deltaUs, epsWs, epss):
    E_star, F_star = monomial_exponents(chi0s, deltaUs, epsWs, epss)
    Mstar = sp.Matrix([
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + E_star), 0, 2 * E_star, F_star - E_star, -1, -(2 + E_star), 1, -F_star],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ])
    return Mstar, E_star, F_star
