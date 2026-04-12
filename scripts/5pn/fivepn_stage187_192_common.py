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


def base_symbols(prefix: str = ""):
    lamW, c_etaU, gamma, KU = sp.symbols(f"{prefix}lambda_W {prefix}c_etaU {prefix}gamma {prefix}K_U", positive=True, real=True)
    Keta, KW, muW, TU = sp.symbols(f"{prefix}K_eta {prefix}K_W {prefix}mu_W {prefix}T_U", positive=True, real=True)
    L, sigma = sp.symbols(f"{prefix}L {prefix}sigma", positive=True, real=True)
    chi0s, deltaUs = sp.symbols(f"{prefix}chi0_star {prefix}deltaU_star", positive=True, real=True)
    Estar, Fstar = sp.symbols(f"{prefix}E_star {prefix}F_star", real=True)
    Ctrs, Cnts, epss = sp.symbols(f"{prefix}C_tr_star {prefix}C_nt_star {prefix}epsilon_eta_star", positive=True, real=True)
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


def orbit_generators():
    Lam, C, Gam, U, W = sp.symbols("Lambda C Gamma U W", real=True)
    return {"Lam": Lam, "C": C, "Gam": Gam, "U": U, "W": W}


def orbit_exponents(syms, pars):
    chi0s = syms["chi0s"]
    deltaUs = syms["deltaUs"]
    Estar = syms["Estar"]
    Fstar = syms["Fstar"]
    Lam, C, Gam, U, W = pars["Lam"], pars["C"], pars["Gam"], pars["U"], pars["W"]
    alpha = sp.simplify((1 + deltaUs) / (1 + chi0s))
    exps = {
        "lamW": Lam,
        "c_etaU": C,
        "gamma": Gam,
        "KU": U,
        "Keta": sp.simplify(2*C - U),
        "KW": W,
        "muW": sp.simplify(2*C - U + 2*W - 2*Lam - Estar*(2*Gam + 2*Lam - U - W) - Fstar*alpha*(Gam + C - U)),
        "TU": sp.simplify(U - alpha*(Gam + C - U)),
    }
    return exps


def orbit_transform(syms, pars):
    exps = orbit_exponents(syms, pars)
    out = dict(syms)
    out["lamW"] = sp.simplify(syms["lamW"] * sp.exp(exps["lamW"]))
    out["c_etaU"] = sp.simplify(syms["c_etaU"] * sp.exp(exps["c_etaU"]))
    out["gamma"] = sp.simplify(syms["gamma"] * sp.exp(exps["gamma"]))
    out["KU"] = sp.simplify(syms["KU"] * sp.exp(exps["KU"]))
    out["Keta"] = sp.simplify(syms["Keta"] * sp.exp(exps["Keta"]))
    out["KW"] = sp.simplify(syms["KW"] * sp.exp(exps["KW"]))
    out["muW"] = sp.simplify(syms["muW"] * sp.exp(exps["muW"]))
    out["TU"] = sp.simplify(syms["TU"] * sp.exp(exps["TU"]))
    return out


def free_ratio_symbols():
    Rlam, Rc, Rgam, RU, RW = sp.symbols("R_lambda R_c R_gamma R_U R_W", positive=True, real=True)
    return {"Rlam": Rlam, "Rc": Rc, "Rgam": Rgam, "RU": RU, "RW": RW}


def orbit_ratio_laws(syms, ratios):
    chi0s = syms["chi0s"]
    deltaUs = syms["deltaUs"]
    Estar = syms["Estar"]
    Fstar = syms["Fstar"]
    Rlam, Rc, Rgam, RU, RW = ratios["Rlam"], ratios["Rc"], ratios["Rgam"], ratios["RU"], ratios["RW"]
    alpha = sp.simplify((1 + deltaUs) / (1 + chi0s))
    RKeta = sp.simplify(Rc**2 / RU)
    RTU = sp.simplify(RU * (RU / (Rgam * Rc)) ** alpha)
    Rmu = sp.simplify((RKeta * RW**2 / Rlam**2) * (Rgam**2 * Rlam**2 / (RU * RW)) ** (-Estar) * (RTU / RU) ** Fstar)
    return {"RKeta": RKeta, "RTU": RTU, "Rmu": Rmu}


def mismatch_ratio_symbols():
    mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)
    return {"mT": mT, "mK": mK, "mMu": mMu}
