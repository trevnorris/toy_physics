
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
    return {
        "lamW": lamW, "c_etaU": c_etaU, "gamma": gamma, "KU": KU,
        "Keta": Keta, "KW": KW, "muW": muW, "TU": TU,
        "L": L, "sigma": sigma, "chi0s": chi0s, "deltaUs": deltaUs,
        "Estar": Estar, "Fstar": Fstar,
    }

def invariants_expr(syms):
    lamW = syms["lamW"]; c_etaU = syms["c_etaU"]; gamma = syms["gamma"]; KU = syms["KU"]
    Keta = syms["Keta"]; KW = syms["KW"]; muW = syms["muW"]; TU = syms["TU"]
    L = syms["L"]; sigma = syms["sigma"]; chi0s = syms["chi0s"]; deltaUs = syms["deltaUs"]
    Estar = syms["Estar"]; Fstar = syms["Fstar"]
    Ctr = sp.simplify((gamma * c_etaU / KU) ** (1 + deltaUs) * (sp.pi**2 * TU / (L**2 * KU)) ** (1 + chi0s))
    Cnt = sp.simplify(
        (lamW**2 * muW / (Keta * KW**2))
        * (gamma**2 * lamW**2 * sigma / (KU * KW)) ** Estar
        * (sp.pi**2 * TU / (L**2 * KU)) ** (-Fstar)
    )
    eps_eta = sp.simplify(c_etaU**2 / (KU * Keta))
    return {"Ctr": Ctr, "Cnt": Cnt, "eps_eta": eps_eta}

def free_ratio_symbols(prefix: str = ""):
    Rlam, Rc, Rgam, RU, RW = sp.symbols(f"{prefix}R_lambda {prefix}R_c {prefix}R_gamma {prefix}R_U {prefix}R_W", positive=True, real=True)
    return {"Rlam": Rlam, "Rc": Rc, "Rgam": Rgam, "RU": RU, "RW": RW}

def orbit_ratio_laws(syms, ratios):
    chi0s = syms["chi0s"]; deltaUs = syms["deltaUs"]; Estar = syms["Estar"]; Fstar = syms["Fstar"]
    Rlam, Rc, Rgam, RU, RW = ratios["Rlam"], ratios["Rc"], ratios["Rgam"], ratios["RU"], ratios["RW"]
    alpha = sp.simplify((1 + deltaUs) / (1 + chi0s))
    RKeta = sp.simplify(Rc**2 / RU)
    RTU = sp.simplify(RU * (RU / (Rgam * Rc)) ** alpha)
    Rmu = sp.simplify((RKeta * RW**2 / Rlam**2) * (Rgam**2 * Rlam**2 / (RU * RW)) ** (-Estar) * (RTU / RU) ** Fstar)
    return {"RKeta": RKeta, "RTU": RTU, "Rmu": Rmu}

def actual_pairwise_ratios(xs, ys):
    return {
        "Rlam": sp.simplify(ys["lamW"] / xs["lamW"]),
        "Rc": sp.simplify(ys["c_etaU"] / xs["c_etaU"]),
        "Rgam": sp.simplify(ys["gamma"] / xs["gamma"]),
        "RU": sp.simplify(ys["KU"] / xs["KU"]),
        "RW": sp.simplify(ys["KW"] / xs["KW"]),
        "RKeta_act": sp.simplify(ys["Keta"] / xs["Keta"]),
        "RTU_act": sp.simplify(ys["TU"] / xs["TU"]),
        "Rmu_act": sp.simplify(ys["muW"] / xs["muW"]),
    }

def pairwise_mismatch_ratios(xs, ys, syms):
    rats = actual_pairwise_ratios(xs, ys)
    free = {"Rlam": rats["Rlam"], "Rc": rats["Rc"], "Rgam": rats["Rgam"], "RU": rats["RU"], "RW": rats["RW"]}
    orbit = orbit_ratio_laws(syms, free)
    mT = sp.simplify(rats["RTU_act"] / orbit["RTU"])
    mK = sp.simplify(rats["RKeta_act"] / orbit["RKeta"])
    mMu = sp.simplify(rats["Rmu_act"] / orbit["Rmu"])
    return {"rats": rats, "orbit": orbit, "mT": mT, "mK": mK, "mMu": mMu}

def pairwise_invariant_ratios(xs, ys):
    invx = invariants_expr(xs); invy = invariants_expr(ys)
    return {
        "RCtr": sp.simplify(invy["Ctr"] / invx["Ctr"]),
        "RCnt": sp.simplify(invy["Cnt"] / invx["Cnt"]),
        "Reps": sp.simplify(invy["eps_eta"] / invx["eps_eta"]),
    }

def q_from_mismatch(syms, mT, mK, mMu):
    chi0s = syms["chi0s"]; Fstar = syms["Fstar"]
    qtr = sp.simplify((1 + chi0s) * sp.log(mT))
    qeta = sp.simplify(-sp.log(mK))
    qnt = sp.simplify(sp.log(mMu) - sp.log(mK) - Fstar * sp.log(mT))
    return {"qtr": qtr, "qnt": qnt, "qeta": qeta}

def mismatch_from_invariant_ratios(syms, RCtr, RCnt, Reps):
    chi0s = syms["chi0s"]; Fstar = syms["Fstar"]
    mT = sp.simplify(RCtr ** (sp.Integer(1) / (1 + chi0s)))
    mK = sp.simplify(1 / Reps)
    mMu = sp.simplify(RCnt * mK * mT**Fstar)
    return {"mT": mT, "mK": mK, "mMu": mMu}

def q_from_invariant_ratios(RCtr, RCnt, Reps):
    return {"qtr": sp.simplify(sp.log(RCtr)), "qnt": sp.simplify(sp.log(RCnt)), "qeta": sp.simplify(sp.log(Reps))}

def linear_observable_map(syms, qtr, qnt, qeta):
    chi0s = syms["chi0s"]; deltaUs = syms["deltaUs"]; eps_eta_star = sp.symbols("epsilon_eta_star", positive=True, real=True)
    coeff_theta = sp.simplify(chi0s * deltaUs / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs)))
    coeff_xi = sp.simplify(2 * chi0s / ((1 + chi0s) * (1 + deltaUs)))
    Theta = sp.simplify(-coeff_theta * qtr)
    Xi = sp.simplify(coeff_xi * qtr + qnt)
    Rsum = sp.simplify(-eps_eta_star / (1 - eps_eta_star) * qeta)
    R1 = sp.simplify(Rsum - Xi)
    return {"Theta": Theta, "Xi": Xi, "Rsum": Rsum, "R1": R1, "eps_eta_star": eps_eta_star}

def residual_log_symbols():
    t, k, m = sp.symbols("t_res k_res mu_res", real=True)
    return {"t": t, "k": k, "m": m}

def orbit_distance_quadratic_form(syms):
    chi0s = syms["chi0s"]; Fstar = syms["Fstar"]
    A = sp.Matrix([
        [1 + chi0s, 0, 0],
        [-Fstar, -1, 1],
        [0, -1, 0],
    ])
    Q = sp.simplify(A.T * A)
    return {"A": A, "Q": Q}
