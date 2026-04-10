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


def expect_zero(name: str, expr):
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


def base_symbols():
    lamW, c_etaU, gamma, KU = sp.symbols("lambda_W c_etaU gamma K_U", positive=True, real=True)
    Keta, KW, muW, TU = sp.symbols("K_eta K_W mu_W T_U", positive=True, real=True)
    L, sigma = sp.symbols("L sigma", positive=True, real=True)
    chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
    Estar, Fstar = sp.symbols("E_star F_star", real=True)
    Ctrs, Cnts, epss = sp.symbols("C_tr_star C_nt_star epsilon_eta_star", positive=True, real=True)
    return {
        "lamW": lamW,
        "c_etaU": c_etaU,
        "gamma": gamma,
        "KU": KU,
        "Keta": Keta,
        "KW": KW,
        "muW": muW,
        "TU": TU,
        "L": L,
        "sigma": sigma,
        "chi0s": chi0s,
        "deltaUs": deltaUs,
        "Estar": Estar,
        "Fstar": Fstar,
        "Ctrs": Ctrs,
        "Cnts": Cnts,
        "epss": epss,
    }


def invariants_expr(syms):
    lamW = syms["lamW"]
    c_etaU = syms["c_etaU"]
    gamma = syms["gamma"]
    KU = syms["KU"]
    Keta = syms["Keta"]
    KW = syms["KW"]
    muW = syms["muW"]
    TU = syms["TU"]
    L = syms["L"]
    sigma = syms["sigma"]
    chi0s = syms["chi0s"]
    deltaUs = syms["deltaUs"]
    Estar = syms["Estar"]
    Fstar = syms["Fstar"]

    Ctr = sp.simplify((gamma * c_etaU / KU) ** (1 + deltaUs) * (sp.pi**2 * TU / (L**2 * KU)) ** (1 + chi0s))
    Cnt = sp.simplify(
        (lamW**2 * muW / (Keta * KW**2))
        * (gamma**2 * lamW**2 * sigma / (KU * KW)) ** Estar
        * (sp.pi**2 * TU / (L**2 * KU)) ** (-Fstar)
    )
    eps_eta = sp.simplify(c_etaU**2 / (KU * Keta))
    return {"Ctr": Ctr, "Cnt": Cnt, "eps_eta": eps_eta}


def finite_orbit_dependent_triple(syms):
    lamW = syms["lamW"]
    c_etaU = syms["c_etaU"]
    gamma = syms["gamma"]
    KU = syms["KU"]
    KW = syms["KW"]
    L = syms["L"]
    sigma = syms["sigma"]
    chi0s = syms["chi0s"]
    deltaUs = syms["deltaUs"]
    Estar = syms["Estar"]
    Fstar = syms["Fstar"]
    Ctrs = syms["Ctrs"]
    Cnts = syms["Cnts"]
    epss = syms["epss"]

    Keta_orbit = sp.simplify(c_etaU**2 / (KU * epss))
    TU_orbit = sp.simplify(
        (L**2 * KU / sp.pi**2)
        * (Ctrs / (gamma * c_etaU / KU) ** (1 + deltaUs)) ** (sp.Integer(1) / (1 + chi0s))
    )
    muW_orbit = sp.simplify(
        (Cnts * Keta_orbit * KW**2 / lamW**2)
        * (gamma**2 * lamW**2 * sigma / (KU * KW)) ** (-Estar)
        * (sp.pi**2 * TU_orbit / (L**2 * KU)) ** Fstar
    )
    return {"Keta_orbit": Keta_orbit, "TU_orbit": TU_orbit, "muW_orbit": muW_orbit}
