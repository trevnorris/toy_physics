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


banner("STAGE 186 — FREE-QUINTUPLE SCALAR CLOSURE SLICE AND CROSSING THEOREM")

chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)
lam1, c1, gam1, kapU, kapW = sp.symbols(
    "lambda1 c1 gamma1 kappaU kappaW", real=True
)
ET, EK, Emu = sp.symbols("E_T E_K E_mu", real=True)
qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

# Ordered drift basis:
# (Delta_lambda, Delta_c, Delta_gamma, Delta_U, Delta_Keta, Delta_W, Delta_mu, Delta_T)
Mstar = sp.Matrix(
    [
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)

subbanner("I. Carry-forward Stage-175 monomial-drift map")
print("M_* =")
sp.pprint(Mstar)

subbanner("II. Exact graph-family dependent log-tangent formulas")
A = sp.simplify((1 + deltaUs) / (1 + chi0s))
dln_delta_graph = sp.simplify(-A * (gam1 + c1 - kapU))
tau1_graph = sp.simplify(kapU + dln_delta_graph)
keta1_graph = sp.simplify(2 * c1 - kapU)
mu1_graph = sp.simplify(
    2 * c1
    - kapU
    + 2 * kapW
    - 2 * lam1
    - Estar * (2 * gam1 + 2 * lam1 - kapU - kapW)
    + Fstar * dln_delta_graph
)

print("d ln delta_U^graph / d tau =")
sp.pprint(dln_delta_graph)
print("tau_1^graph = d ln T_U^graph / d tau =")
sp.pprint(tau1_graph)
print("kappa_eta^graph = d ln K_eta^graph / d tau =")
sp.pprint(keta1_graph)
print("mu_1^graph = d ln mu_W^graph / d tau =")
sp.pprint(mu1_graph)

expect_zero(
    "tau_1^graph - [kappa_U - ((1+deltaU_*)/(1+chi0_*))(gamma1+c1-kappa_U)]",
    tau1_graph - (kapU - ((1 + deltaUs) / (1 + chi0s)) * (gam1 + c1 - kapU)),
)
expect_zero("kappa_eta^graph - (2 c1 - kappa_U)", keta1_graph - (2 * c1 - kapU))
expect_zero(
    "mu_1^graph - Stage-175 q=0 formula",
    mu1_graph
    - (
        2 * c1
        - kapU
        + 2 * kapW
        - 2 * lam1
        - Estar * (2 * gam1 + 2 * lam1 - kapU - kapW)
        - Fstar * ((1 + deltaUs) / (1 + chi0s)) * (gam1 + c1 - kapU)
    ),
)

subbanner("III. Exact graph-family orbit-kernel theorem")
dx_graph = sp.Matrix(
    [lam1, c1, gam1, kapU, keta1_graph, kapW, mu1_graph, tau1_graph]
)
print("dot(Delta x)_graph =")
sp.pprint(dx_graph)
expect_zero("M_* dot(Delta x)_graph", Mstar * dx_graph)

subbanner("IV. Same-free-quintuple graph-error packet -> quotient packet")
err_vec = sp.Matrix([0, 0, 0, 0, EK, 0, Emu, ET])
q_from_E = sp.simplify(Mstar * err_vec)
q_expected = sp.Matrix([(1 + chi0s) * ET, Emu - EK - Fstar * ET, -EK])
print("Delta x_err =")
sp.pprint(err_vec)
print("q_from_E = M_* Delta x_err =")
sp.pprint(q_from_E)
expect_zero("quotient packet from graph errors", q_from_E - q_expected)

subbanner("V. Exact inverse compiler from quotient packet back to graph errors")
ET_inv = sp.simplify(qtr / (1 + chi0s))
EK_inv = sp.simplify(-qeta)
Emu_inv = sp.simplify(qnt - qeta + Fstar * qtr / (1 + chi0s))
print("E_T(q) =")
sp.pprint(ET_inv)
print("E_K(q) =")
sp.pprint(EK_inv)
print("E_mu(q) =")
sp.pprint(Emu_inv)

q_subs = {qtr: q_expected[0], qnt: q_expected[1], qeta: q_expected[2]}
expect_zero("inverse compiler for E_T", ET_inv.subs(q_subs) - ET)
expect_zero("inverse compiler for E_K", EK_inv.subs(q_subs) - EK)
expect_zero("inverse compiler for E_mu", Emu_inv.subs(q_subs) - Emu)

subbanner("VI. Exact repair vector kills the quotient packet")
repair_from_q = sp.Matrix(
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
print("Delta x_rep(q) =")
sp.pprint(repair_from_q)
expect_zero(
    "M_* Delta x_rep(q) + q",
    Mstar * repair_from_q + sp.Matrix([qtr, qnt, qeta]),
)

subbanner("VII. Exact same-free-quintuple decomposition")
candidate_minus_graph = sp.simplify(dx_graph + err_vec - dx_graph)
expect_zero("candidate - graph - error packet", candidate_minus_graph - err_vec)
expect_zero(
    "quotient packet of (graph tangent + error packet) - q_from_E",
    Mstar * (dx_graph + err_vec) - q_expected,
)

banner("STAGE 186 SYMPY AUDIT PASSED")
