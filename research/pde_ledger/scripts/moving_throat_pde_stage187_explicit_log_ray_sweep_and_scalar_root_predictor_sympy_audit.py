#!/usr/bin/env python3
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


banner("STAGE 187 — EXPLICIT LOG-RAY SWEEP AND SCALAR ROOT PREDICTOR")

# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------
chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)

tau = sp.symbols("tau", real=True)
slam, sc, sgam, sU, sW = sp.symbols(
    "s_lambda s_c s_gamma s_U s_W", real=True
)

lam0, c0, gam0, KU0, KW0 = sp.symbols(
    "lambda0 c0 gamma0 KU0 KW0", positive=True, real=True
)
CtrT, CntT, epsEtaT = sp.symbols(
    "Ctr_target Cnt_target epsEta_target", positive=True, real=True
)
Lgeom, sig = sp.symbols("L sigma", positive=True, real=True)

aStar = sp.simplify((1 + deltaUs) / (1 + chi0s))
sigma_delta = sp.simplify(-aStar * (sgam + sc - sU))
sigma_T = sp.simplify(sU + sigma_delta)
sigma_Keta = sp.simplify(2 * sc - sU)
sigma_mu = sp.simplify(
    2 * sc
    - sU
    + 2 * sW
    - 2 * slam
    - Estar * (2 * sgam + 2 * slam - sU - sW)
    + Fstar * sigma_delta
)

subbanner("I. Exact free-quintuple log ray")
lam_tau = sp.simplify(lam0 * sp.exp(slam * tau))
c_tau = sp.simplify(c0 * sp.exp(sc * tau))
gam_tau = sp.simplify(gam0 * sp.exp(sgam * tau))
KU_tau = sp.simplify(KU0 * sp.exp(sU * tau))
KW_tau = sp.simplify(KW0 * sp.exp(sW * tau))

print("lambda_W(tau) =")
sp.pprint(lam_tau)
print("c_etaU(tau) =")
sp.pprint(c_tau)
print("gamma(tau) =")
sp.pprint(gam_tau)
print("K_U(tau) =")
sp.pprint(KU_tau)
print("K_W(tau) =")
sp.pprint(KW_tau)
expect_zero("d/dtau ln lambda_W - s_lambda", sp.diff(sp.log(lam_tau), tau) - slam)
expect_zero("d/dtau ln c_etaU - s_c", sp.diff(sp.log(c_tau), tau) - sc)
expect_zero("d/dtau ln gamma - s_gamma", sp.diff(sp.log(gam_tau), tau) - sgam)
expect_zero("d/dtau ln K_U - s_U", sp.diff(sp.log(KU_tau), tau) - sU)
expect_zero("d/dtau ln K_W - s_W", sp.diff(sp.log(KW_tau), tau) - sW)

subbanner("II. Exact finite graph lift along the log ray")
delta_graph_0 = sp.simplify(
    (CtrT / ((gam0 * c0 / KU0) ** (1 + deltaUs))) ** (sp.Rational(1, 1) / (1 + chi0s))
)
delta_graph_tau_direct = sp.simplify(
    (CtrT / ((gam_tau * c_tau / KU_tau) ** (1 + deltaUs))) ** (sp.Rational(1, 1) / (1 + chi0s))
)
delta_graph_tau_exp = sp.simplify(delta_graph_0 * sp.exp(sigma_delta * tau))
print("delta_U^graph(0) =")
sp.pprint(delta_graph_0)
print("sigma_delta =")
sp.pprint(sigma_delta)
print("delta_U^graph(tau) [direct] =")
sp.pprint(delta_graph_tau_direct)
print("delta_U^graph(tau) [exp form] =")
sp.pprint(delta_graph_tau_exp)
expect_zero("delta_U^graph direct - exp form", delta_graph_tau_direct - delta_graph_tau_exp)

T_graph_0 = sp.simplify(Lgeom**2 * KU0 * delta_graph_0 / sp.pi**2)
T_graph_tau_direct = sp.simplify(Lgeom**2 * KU_tau * delta_graph_tau_direct / sp.pi**2)
T_graph_tau_exp = sp.simplify(T_graph_0 * sp.exp(sigma_T * tau))
print("sigma_T =")
sp.pprint(sigma_T)
expect_zero("T_U^graph direct - exp form", T_graph_tau_direct - T_graph_tau_exp)

Keta_graph_0 = sp.simplify(c0**2 / (KU0 * epsEtaT))
Keta_graph_tau_direct = sp.simplify(c_tau**2 / (KU_tau * epsEtaT))
Keta_graph_tau_exp = sp.simplify(Keta_graph_0 * sp.exp(sigma_Keta * tau))
print("sigma_Keta =")
sp.pprint(sigma_Keta)
expect_zero("K_eta^graph direct - exp form", Keta_graph_tau_direct - Keta_graph_tau_exp)

mu_graph_0 = sp.simplify(
    CntT * c0**2 * KW0**2 / (epsEtaT * KU0 * lam0**2)
    * ((gam0**2 * lam0**2 * sig) / (KU0 * KW0)) ** (-Estar)
    * delta_graph_0**Fstar
)
mu_graph_tau_direct = sp.simplify(
    CntT * c_tau**2 * KW_tau**2 / (epsEtaT * KU_tau * lam_tau**2)
    * ((gam_tau**2 * lam_tau**2 * sig) / (KU_tau * KW_tau)) ** (-Estar)
    * delta_graph_tau_direct**Fstar
)
mu_graph_tau_exp = sp.simplify(mu_graph_0 * sp.exp(sigma_mu * tau))
print("sigma_mu =")
sp.pprint(sigma_mu)
expect_zero("mu_W^graph direct - exp form", mu_graph_tau_direct - mu_graph_tau_exp)

subbanner("III. Exact finite monomial invariance on the graph ray")
Ctr_tau = sp.simplify((gam_tau * c_tau / KU_tau) ** (1 + deltaUs) * delta_graph_tau_direct ** (1 + chi0s))
Cnt_tau = sp.simplify(
    lam_tau**2 * mu_graph_tau_direct / (Keta_graph_tau_direct * KW_tau**2)
    * ((gam_tau**2 * lam_tau**2 * sig) / (KU_tau * KW_tau)) ** Estar
    * delta_graph_tau_direct ** (-Fstar)
)
epsEta_tau = sp.simplify(c_tau**2 / (KU_tau * Keta_graph_tau_direct))
print("Ctr(tau) =")
sp.pprint(Ctr_tau)
print("Cnt(tau) =")
sp.pprint(Cnt_tau)
print("eps_eta(tau) =")
sp.pprint(epsEta_tau)
expect_zero("Ctr(tau) - Ctr_target", Ctr_tau - CtrT)
expect_zero("Cnt(tau) - Cnt_target", Cnt_tau - CntT)
expect_zero("eps_eta(tau) - epsEta_target", epsEta_tau - epsEtaT)

subbanner("IV. Stage-175 quotient map annihilates the constant ray tangent")
Mstar = sp.Matrix(
    [
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)
dx_ray = sp.Matrix([slam, sc, sgam, sU, sigma_Keta, sW, sigma_mu, sigma_T])
print("M_* =")
sp.pprint(Mstar)
print("dot(Delta x)_ray =")
sp.pprint(dx_ray)
expect_zero("M_* dot(Delta x)_ray", Mstar * dx_ray)

subbanner("V. Primitive direction table")
primitive_data = {
    "e_lambda": {slam: 1, sc: 0, sgam: 0, sU: 0, sW: 0},
    "e_c": {slam: 0, sc: 1, sgam: 0, sU: 0, sW: 0},
    "e_gamma": {slam: 0, sc: 0, sgam: 1, sU: 0, sW: 0},
    "e_U": {slam: 0, sc: 0, sgam: 0, sU: 1, sW: 0},
    "e_W": {slam: 0, sc: 0, sgam: 0, sU: 0, sW: 1},
}
expected_table = {
    "e_lambda": (0, 0, 0, -2 - 2 * Estar),
    "e_c": (-aStar, -aStar, 2, 2 - Fstar * aStar),
    "e_gamma": (-aStar, -aStar, 0, -2 * Estar - Fstar * aStar),
    "e_U": (aStar, 1 + aStar, -1, -1 + Estar + Fstar * aStar),
    "e_W": (0, 0, 0, 2 + Estar),
}
for name, subs in primitive_data.items():
    actual = (
        sp.simplify(sigma_delta.subs(subs)),
        sp.simplify(sigma_T.subs(subs)),
        sp.simplify(sigma_Keta.subs(subs)),
        sp.simplify(sigma_mu.subs(subs)),
    )
    print(f"{name} exponents = {actual}")
    for idx, (act, exp) in enumerate(zip(actual, expected_table[name]), start=1):
        expect_zero(f"{name} exponent {idx}", act - exp)

subbanner("VI. Log-ray library is first-order complete")
Y0, Y1, eps_step = sp.symbols("Y0 Y1 eps_step", positive=True, real=True)
sY = sp.simplify(Y1 / Y0)
ray_component = sp.simplify(Y0 * sp.exp(sY * eps_step))
series_match = sp.series(ray_component, eps_step, 0, 2).removeO()
expect_zero("first-order local ray match", series_match - (Y0 + Y1 * eps_step))

subbanner("VII. Exact affine and log-linear root predictors")
Phi0, L0, eps = sp.symbols("Phi0 L0 eps", positive=True, real=True)
Phi0_eps = sp.simplify(1 + eps)
Phi1 = sp.simplify(L0 * Phi0_eps)
tau_aff = sp.simplify((1 - Phi0_eps) / Phi1)
tau_log = sp.simplify(-sp.log(Phi0_eps) / L0)
print("tau_aff =")
sp.pprint(tau_aff)
print("tau_log =")
sp.pprint(tau_log)

diff_series = sp.series(sp.expand(tau_log - tau_aff), eps, 0, 4).removeO()
print("series[tau_log - tau_aff] =")
sp.pprint(diff_series)
expect_zero(
    "predictor difference starts at quadratic order",
    diff_series + eps**2 / (2 * L0) - 2 * eps**3 / (3 * L0),
)

banner("STAGE 187 SYMPY AUDIT PASSED")
