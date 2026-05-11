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


def _normalize(expr):
    expr = sp.expand(expr)
    expr = sp.expand_power_base(expr, force=True)
    expr = sp.expand_power_exp(expr)
    expr = sp.expand_log(expr, force=True)
    expr = sp.powsimp(expr, force=True)
    expr = sp.cancel(expr)
    return sp.simplify(expr)


def normalize(expr):
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(_normalize)
    return _normalize(expr)


def expect_zero(name: str, expr) -> None:
    simplified = normalize(expr)
    if isinstance(simplified, sp.MatrixBase):
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
        return

    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


def expect_positive(name: str, expr) -> None:
    simplified = normalize(sp.factor(expr))
    print(f"{name} > 0 = {simplified}")
    if simplified.is_positive is not True:
        raise AssertionError(f"{name} is not provably positive")


def expect_negative(name: str, expr) -> None:
    simplified = normalize(sp.factor(expr))
    print(f"{name} < 0 = {simplified}")
    if simplified.is_negative is not True:
        raise AssertionError(f"{name} is not provably negative")


banner("STAGE 186 — FREE-QUINTUPLE SCALAR CLOSURE SLICE AND CROSSING THEOREM")

chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)

# Ordered drift basis:
# (Delta_lambda, Delta_c, Delta_gamma, Delta_U, Delta_Keta, Delta_W, Delta_mu, Delta_T)
Mstar = sp.Matrix(
    [
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)

subbanner("I. Stage 192 dependent pivot block and exact canonical section")

T_idx = 7
Keta_idx = 4
mu_idx = 6

Pdep = sp.Matrix.hstack(Mstar[:, T_idx], Mstar[:, Keta_idx], Mstar[:, mu_idx])
Pdep_inv = sp.simplify(Pdep.inv())

Edep = sp.zeros(8, 3)
Edep[T_idx, 0] = 1
Edep[Keta_idx, 1] = 1
Edep[mu_idx, 2] = 1

Sdep = sp.simplify(Edep * Pdep_inv)

print("M_* =")
sp.pprint(Mstar)
print("P_(T,K_eta,mu) =")
sp.pprint(Pdep)
print("P_(T,K_eta,mu)^(-1) =")
sp.pprint(Pdep_inv)
print("S_(T,K_eta,mu) =")
sp.pprint(Sdep)

expect_zero("M_* S - I_3", Mstar * Sdep - sp.eye(3))

subbanner("II. Exact graph-family tangent formulas from the Stage 202 target graph")

Ctr_tgt, Cnt_tgt, epsEta_tgt = sp.symbols(
    "Ctr_target Cnt_target epsEta_target", positive=True, real=True
)
L, sigma = sp.symbols("L sigma", positive=True, real=True)
lam0, cetaU0, gamma0, KU0, KW0 = sp.symbols(
    "lambda0 cetaU0 gamma0 KU0 KW0", positive=True, real=True
)
t = sp.symbols("t", real=True)
lam1, c1, gam1, kapU, kapW = sp.symbols(
    "lambda1 c1 gamma1 kappaU kappaW", real=True
)

lam_t = lam0 * sp.exp(t * lam1)
cetaU_t = cetaU0 * sp.exp(t * c1)
gamma_t = gamma0 * sp.exp(t * gam1)
KU_t = KU0 * sp.exp(t * kapU)
KW_t = KW0 * sp.exp(t * kapW)

alpha = sp.simplify((1 + deltaUs) / (1 + chi0s))
deltaU_graph_t = sp.simplify(
    (Ctr_tgt / (gamma_t * cetaU_t / KU_t) ** (1 + deltaUs))
    ** (sp.Rational(1, 1) / (1 + chi0s))
)
T_graph_t = sp.simplify(L**2 * KU_t * deltaU_graph_t / sp.pi**2)
Keta_graph_t = sp.simplify(cetaU_t**2 / (KU_t * epsEta_tgt))
mu_graph_t = sp.simplify(
    Cnt_tgt
    * Keta_graph_t
    * KW_t**2
    / lam_t**2
    * ((gamma_t**2 * lam_t**2 * sigma) / (KU_t * KW_t)) ** (-Estar)
    * deltaU_graph_t**Fstar
)

dln_delta_graph = normalize(sp.diff(sp.log(deltaU_graph_t), t).subs(t, 0))
tau1_graph = normalize(sp.diff(sp.log(T_graph_t), t).subs(t, 0))
keta1_graph = normalize(sp.diff(sp.log(Keta_graph_t), t).subs(t, 0))
mu1_graph = normalize(sp.diff(sp.log(mu_graph_t), t).subs(t, 0))

print("d ln delta_U^graph / d tau =")
sp.pprint(dln_delta_graph)
print("tau_1^graph =")
sp.pprint(tau1_graph)
print("kappa_eta^graph =")
sp.pprint(keta1_graph)
print("mu_1^graph =")
sp.pprint(mu1_graph)

expect_zero(
    "dln delta_U^graph - carried formula",
    dln_delta_graph + alpha * (gam1 + c1 - kapU),
)
expect_zero(
    "tau_1^graph - carried formula",
    tau1_graph - (kapU - alpha * (gam1 + c1 - kapU)),
)
expect_zero("kappa_eta^graph - (2 c1 - kappa_U)", keta1_graph - (2 * c1 - kapU))
expect_zero(
    "mu_1^graph - carried formula",
    mu1_graph
    - (
        2 * c1
        - kapU
        + 2 * kapW
        - 2 * lam1
        - Estar * (2 * gam1 + 2 * lam1 - kapU - kapW)
        - Fstar * alpha * (gam1 + c1 - kapU)
    ),
)

subbanner("III. Exact graph-family orbit-kernel theorem")

dx_graph = sp.Matrix([lam1, c1, gam1, kapU, keta1_graph, kapW, mu1_graph, tau1_graph])
print("dot(Delta x)_graph =")
sp.pprint(dx_graph)
expect_zero("M_* dot(Delta x)_graph", Mstar * dx_graph)

subbanner("IV. Exact graph-error packet from the direct monomials")

lam, cetaU, gamma, KU, KW = sp.symbols(
    "lambda_W c_etaU gamma K_U K_W", positive=True, real=True
)
E_T, E_K, E_mu = sp.symbols("E_T E_K E_mu", real=True)

deltaU_graph = sp.simplify(
    (Ctr_tgt / (gamma * cetaU / KU) ** (1 + deltaUs))
    ** (sp.Rational(1, 1) / (1 + chi0s))
)
T_graph = sp.simplify(L**2 * KU * deltaU_graph / sp.pi**2)
Keta_graph = sp.simplify(cetaU**2 / (KU * epsEta_tgt))
mu_graph = sp.simplify(
    Cnt_tgt
    * Keta_graph
    * KW**2
    / lam**2
    * ((gamma**2 * lam**2 * sigma) / (KU * KW)) ** (-Estar)
    * deltaU_graph**Fstar
)

TU = sp.simplify(T_graph * sp.exp(E_T))
Keta = sp.simplify(Keta_graph * sp.exp(E_K))
muW = sp.simplify(mu_graph * sp.exp(E_mu))

Ctr = sp.simplify(
    (gamma * cetaU / KU) ** (1 + deltaUs)
    * (sp.pi**2 * TU / (L**2 * KU)) ** (1 + chi0s)
)
Cnt = sp.simplify(
    (lam**2 * muW / (Keta * KW**2))
    * ((gamma**2 * lam**2 * sigma) / (KU * KW)) ** Estar
    * (sp.pi**2 * TU / (L**2 * KU)) ** (-Fstar)
)
epsEta = sp.simplify(cetaU**2 / (KU * Keta))

qtr = normalize(sp.log(Ctr / Ctr_tgt))
qnt = normalize(sp.log(Cnt / Cnt_tgt))
qeta = normalize(sp.log(epsEta / epsEta_tgt))

print("q_tr =")
sp.pprint(qtr)
print("q_nt =")
sp.pprint(qnt)
print("q_eta =")
sp.pprint(qeta)

expect_zero("q_tr - (1+chi0_*) E_T", qtr - (1 + chi0s) * E_T)
expect_zero("q_nt - (E_mu - E_K - F_* E_T)", qnt - (E_mu - E_K - Fstar * E_T))
expect_zero("q_eta + E_K", qeta + E_K)

subbanner("V. Exact inverse compiler and repair vector")

qtr_s, qnt_s, qeta_s = sp.symbols("q_tr q_nt q_eta", real=True)
ET_s, EK_s, Emu_s = sp.symbols("E_T_sol E_K_sol E_mu_sol", real=True)

inverse_solution = sp.solve(
    [
        sp.Eq((1 + chi0s) * ET_s, qtr_s),
        sp.Eq(Emu_s - EK_s - Fstar * ET_s, qnt_s),
        sp.Eq(-EK_s, qeta_s),
    ],
    [ET_s, EK_s, Emu_s],
    dict=True,
)[0]

print("Inverse compiler solution =")
sp.pprint(inverse_solution)

expect_zero("inverse E_T", inverse_solution[ET_s] - qtr_s / (1 + chi0s))
expect_zero("inverse E_K", inverse_solution[EK_s] + qeta_s)
expect_zero(
    "inverse E_mu",
    inverse_solution[Emu_s] - (qnt_s - qeta_s + Fstar * qtr_s / (1 + chi0s)),
)

q_from_errors = sp.Matrix([(1 + chi0s) * E_T, E_mu - E_K - Fstar * E_T, -E_K])
dx_rep = normalize(-Sdep * q_from_errors)
dx_rep_expected = sp.Matrix([0, 0, 0, 0, -E_K, 0, -E_mu, -E_T])

print("Delta x_rep =")
sp.pprint(dx_rep)
expect_zero("repair vector from graph errors", dx_rep - dx_rep_expected)
expect_zero("M_* Delta x_rep + q(E)", Mstar * dx_rep + q_from_errors)

subbanner("VI. Stage 197 scalar closure composed with an explicit Stage 202 graph path")

Siso, beta, Sigma0, Sigma5 = sp.symbols("S beta Sigma_0 Sigma_5", positive=True, real=True)
tau = sp.symbols("tau", real=True)
rho = sp.symbols("rho", positive=True, real=True)
lam_bar, cetaU_bar, gamma_bar, KU_bar, KW_bar = sp.symbols(
    "lambda_bar cetaU_bar gamma_bar KU_bar KW_bar", positive=True, real=True
)

beta_path = sp.simplify(1 + rho * (2 * tau - 1) / (1 + rho))
gamma_path = sp.simplify(gamma_bar * beta_path)
y_graph_path = sp.Matrix([lam_bar, cetaU_bar, gamma_path, KU_bar, KW_bar])

chi_from_stage197 = sp.simplify(3 * (Siso * beta**5 + 9 * Sigma5) / (3 * Siso - Sigma0))
closure_num_stage197 = sp.simplify(3 * Siso * (beta**5 - 1) + Sigma0 + 27 * Sigma5)

hat_chi_graph = sp.simplify(
    chi_from_stage197.subs({beta: beta_path, Sigma0: 0, Sigma5: 0})
)
hat_delta_graph = sp.simplify(hat_chi_graph - 1)
closure_num_graph = sp.simplify(
    closure_num_stage197.subs({beta: beta_path, Sigma0: 0, Sigma5: 0})
)
hat_delta_den = sp.denom(sp.together(hat_delta_graph))

print("Explicit free-quintuple graph path y(tau) =")
sp.pprint(y_graph_path)
print("beta_path(tau) =")
sp.pprint(beta_path)
print("widehat chi_Q(y(tau)) from the carried Stage 197 closure algebra =")
sp.pprint(hat_chi_graph)
print("widehat Delta_Q(tau) =")
sp.pprint(hat_delta_graph)

expect_zero("beta_path - gamma(tau)/gamma_bar", y_graph_path[2] / gamma_bar - beta_path)
expect_zero(
    "Stage 197 closure numerator identity on the graph path",
    sp.simplify(3 * Siso * hat_delta_graph - closure_num_graph),
)
expect_positive("graph residual denominator", hat_delta_den)

hat_delta_minus = sp.factor(hat_delta_graph.subs(tau, 0))
hat_delta_plus = sp.factor(hat_delta_graph.subs(tau, 1))
hat_delta_root = sp.simplify(hat_delta_graph.subs(tau, sp.Rational(1, 2)))
real_crossings = sp.solveset(sp.Eq(hat_delta_graph, 0), tau, domain=sp.S.Reals)

print("widehat Delta_Q(0) =")
sp.pprint(hat_delta_minus)
print("widehat Delta_Q(1) =")
sp.pprint(hat_delta_plus)
print("real crossing set =")
sp.pprint(real_crossings)

expect_negative("graph residual at tau_- = 0", hat_delta_minus)
expect_positive("graph residual at tau_+ = 1", hat_delta_plus)
expect_zero("graph residual at tau_* = 1/2", hat_delta_root)
if real_crossings != sp.FiniteSet(sp.Rational(1, 2)):
    raise AssertionError(f"unexpected real crossing set: {real_crossings}")

banner("STAGE 186 SYMPY AUDIT PASSED")
