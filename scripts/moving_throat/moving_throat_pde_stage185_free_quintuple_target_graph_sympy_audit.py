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


def simplify_log(expr):
    def _one(z):
        z = sp.expand(z)
        z = sp.expand_power_base(z, force=True)
        z = sp.expand_power_exp(z)
        z = sp.expand_log(z, force=True)
        z = sp.powsimp(z, force=True)
        z = sp.cancel(z)
        return sp.simplify(z)
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(_one)
    return _one(expr)


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = simplify_log(expr)
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = simplify_log(expr)
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 185 — EXACT FREE-QUINTUPLE TARGET GRAPH AND THE FIRST REDUCED-FAMILY TEST")

# ---------------------------------------------------------------------------
# I. Carry-forward direct microscopic monomials
# ---------------------------------------------------------------------------
subbanner("I. Direct microscopic monomials and target data")

chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)

Ctr_tgt, Cnt_tgt, epsEta_tgt = sp.symbols(
    "Ctr_target Cnt_target epsEta_target", positive=True, real=True
)

lam, cetaU, gamma, KU, Keta, KW, muW, TU = sp.symbols(
    "lambda_W c_etaU gamma K_U K_eta K_W mu_W T_U", positive=True, real=True
)
L, sigma = sp.symbols("L sigma", positive=True, real=True)

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

print("Ctr_* =")
sp.pprint(Ctr)
print("Cnt_* =")
sp.pprint(Cnt)
print("eps_eta =")
sp.pprint(epsEta)

# ---------------------------------------------------------------------------
# II. Exact target-graph solve over the free quintuple
# ---------------------------------------------------------------------------
subbanner("II. Exact target-graph solve")

deltaU_graph = sp.simplify(
    (Ctr_tgt / (gamma * cetaU / KU) ** (1 + deltaUs)) ** (sp.Rational(1, 1) / (1 + chi0s))
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

print("deltaU_graph =")
sp.pprint(deltaU_graph)
print("T_graph =")
sp.pprint(T_graph)
print("Keta_graph =")
sp.pprint(Keta_graph)
print("mu_graph =")
sp.pprint(mu_graph)

expect_zero(
    "log(target Ctr reconstructed by graph / Ctr_target)",
    sp.log(Ctr.subs({TU: T_graph, Keta: Keta_graph, muW: mu_graph}) / Ctr_tgt),
)
expect_zero(
    "log(target Cnt reconstructed by graph / Cnt_target)",
    sp.log(Cnt.subs({TU: T_graph, Keta: Keta_graph, muW: mu_graph}) / Cnt_tgt),
)
expect_zero(
    "log(target eps_eta reconstructed by graph / epsEta_target)",
    sp.log(epsEta.subs({TU: T_graph, Keta: Keta_graph, muW: mu_graph}) / epsEta_tgt),
)

# Verify mu_graph is independent of L once the tracking monomial is eliminated
print(f"L in mu_graph.free_symbols? {'yes' if L in mu_graph.free_symbols else 'no'}")
if L in mu_graph.free_symbols:
    raise AssertionError("mu_graph should not depend on L")

# ---------------------------------------------------------------------------
# III. Exact equality of graph projection and Stage-184 canonical projection
# ---------------------------------------------------------------------------
subbanner("III. Graph projection equals the canonical orbit projection")

Rtr_expr = sp.simplify(Ctr / Ctr_tgt)
Rnt_expr = sp.simplify(Cnt / Cnt_tgt)
Reta_expr = sp.simplify(epsEta / epsEta_tgt)

x_graph = sp.Matrix([lam, cetaU, gamma, KU, Keta_graph, KW, mu_graph, T_graph])

x_proj_can = sp.Matrix(
    [
        lam,
        cetaU,
        gamma,
        KU,
        Keta * Reta_expr,
        KW,
        muW * Rnt_expr ** (-1) * Reta_expr * Rtr_expr ** (-Fstar / (1 + chi0s)),
        TU * Rtr_expr ** (-sp.Rational(1, 1) / (1 + chi0s)),
    ]
)

print("x_graph =")
sp.pprint(x_graph)
print("x_proj^can =")
sp.pprint(x_proj_can)

expect_zero(
    "graph projection equals canonical orbit projection",
    sp.Matrix([sp.log(x_proj_can[i] / x_graph[i]) for i in range(8)]),
)

# ---------------------------------------------------------------------------
# IV. Exact graph-error packet and equivalence with the Stage-184 mismatch chart
# ---------------------------------------------------------------------------
subbanner("IV. Exact graph-error packet")

qtr = sp.log(Rtr_expr)
qnt = sp.log(Rnt_expr)
qeta = sp.log(Reta_expr)

E_T = sp.simplify(sp.log(TU / T_graph))
E_K = sp.simplify(sp.log(Keta / Keta_graph))
E_mu = sp.simplify(sp.log(muW / mu_graph))

print("E_T =")
sp.pprint(E_T)
print("E_K =")
sp.pprint(E_K)
print("E_mu =")
sp.pprint(E_mu)

expect_zero("E_T - q_tr/(1+chi0_*)", E_T - qtr / (1 + chi0s))
expect_zero("E_K + q_eta", E_K + qeta)
expect_zero(
    "E_mu - (q_nt - q_eta + F_* q_tr/(1+chi0_*))",
    E_mu - (qnt - qeta + Fstar * qtr / (1 + chi0s)),
)

mT = sp.exp(qtr / (1 + chi0s))
mK = sp.exp(-qeta)
mMu = sp.exp(qnt - qeta + Fstar * qtr / (1 + chi0s))

expect_zero("E_T - log(m_T)", E_T - sp.log(mT))
expect_zero("E_K - log(m_K)", E_K - sp.log(mK))
expect_zero("E_mu - log(m_mu)", E_mu - sp.log(mMu))

# ---------------------------------------------------------------------------
# V. Exact repair-vector rewrite in graph-error variables
# ---------------------------------------------------------------------------
subbanner("V. Exact repair vector in graph-error coordinates")

dx_rep = sp.Matrix(
    [
        0,
        0,
        0,
        0,
        qeta,
        0,
        -qnt + qeta - Fstar * qtr / (1 + chi0s),
        -qtr / (1 + chi0s),
    ]
)

dx_rep_graph = sp.Matrix([0, 0, 0, 0, -E_K, 0, -E_mu, -E_T])

print("Delta x_rep (q-chart) =")
sp.pprint(dx_rep)
print("Delta x_rep (graph-error chart) =")
sp.pprint(dx_rep_graph)

expect_zero("repair vector rewrite", dx_rep - dx_rep_graph)

# ---------------------------------------------------------------------------
# VI. First reduced-family packet and its zero set
# ---------------------------------------------------------------------------
subbanner("VI. First reduced-family packet")

chiQ = sp.symbols("chi_Q", real=True)
Delta_family_graph = sp.Matrix([chiQ - 1, E_T, E_K, E_mu])

print("Delta_family^graph =")
sp.pprint(Delta_family_graph)

# On the target graph and chi_Q = 1, the family packet vanishes exactly.
expect_zero(
    "family packet on target graph",
    Delta_family_graph.subs({TU: T_graph, Keta: Keta_graph, muW: mu_graph, chiQ: 1}),
)

# The multiplicative family chart is exactly the ratio to the target graph.
M_family = sp.Matrix([chiQ, sp.exp(E_T), sp.exp(E_K), sp.exp(E_mu)])
print("M_family =")
sp.pprint(M_family)
expect_zero(
    "multiplicative family chart on target graph",
    M_family.subs({TU: T_graph, Keta: Keta_graph, muW: mu_graph, chiQ: 1}) - sp.Matrix([1, 1, 1, 1]),
)

banner("STAGE 185 AUDIT COMPLETE")
print("Result: the target orbit is an exact graph over the five free microscopic coordinates.")
print("Result: the four-scalar reduced-family packet is (chi_Q-1, E_T, E_K, E_mu).")
print("Result: the canonical repair vector is exactly minus the graph-error triple on the dependent coordinates.")
