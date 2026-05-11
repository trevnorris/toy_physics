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
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand_log(sp.expand(z), force=True)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand_log(sp.expand(expr), force=True))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


def simplify_log(expr):
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(lambda z: sp.simplify(sp.expand_log(z, force=True)))
    return sp.simplify(sp.expand_log(expr, force=True))


banner("STAGE 184 — EXPLICIT REALIZATION COMPILER AND CANONICAL ORBIT PROJECTION")

# ---------------------------------------------------------------------------
# I. Carry-forward quotient map and canonical section
# ---------------------------------------------------------------------------
subbanner("I. Exact Stage 192 quotient map and canonical section")

chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)

# Ordered microscopic drift basis:
# (Delta_lambda, Delta_c, Delta_gamma, Delta_U, Delta_Keta, Delta_W, Delta_mu, Delta_T)
Mstar = sp.Matrix(
    [
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)

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

# ---------------------------------------------------------------------------
# II. Intrinsic multiplicative and additive realization packet
# ---------------------------------------------------------------------------
subbanner("II. Intrinsic ratio packet and exact mismatch chart")

chiQ = sp.symbols("chi_Q", real=True)
Rtr, Rnt, Reta = sp.symbols("R_tr R_nt R_eta", positive=True, real=True)
qtr = sp.log(Rtr)
qnt = sp.log(Rnt)
qeta = sp.log(Reta)

Delta_real_int = sp.Matrix([chiQ - 1, qtr, qnt, qeta])
V_real_int = sp.Matrix([chiQ, Rtr, Rnt, Reta])

mT = sp.simplify(sp.exp(qtr / (1 + chi0s)))
mK = sp.simplify(sp.exp(-qeta))
mMu = sp.simplify(sp.exp(qnt - qeta + Fstar * qtr / (1 + chi0s)))
M_real_int = sp.Matrix([chiQ, mT, mK, mMu])

print("Delta_real^int =")
sp.pprint(Delta_real_int)
print("V_real^int =")
sp.pprint(V_real_int)
print("M_real^int =")
sp.pprint(M_real_int)

expect_zero("exp(q_tr) - R_tr", sp.exp(qtr) - Rtr)
expect_zero("exp(q_nt) - R_nt", sp.exp(qnt) - Rnt)
expect_zero("exp(q_eta) - R_eta", sp.exp(qeta) - Reta)
expect_zero("m_T - R_tr^(1/(1+chi0_*))", mT - Rtr ** (sp.Rational(1, 1) / (1 + chi0s)))
expect_zero("m_K - R_eta^(-1)", mK - Reta ** (-1))
expect_zero(
    "m_mu - R_nt R_eta^(-1) R_tr^(F_*/(1+chi0_*))",
    mMu - Rnt * Reta ** (-1) * Rtr ** (Fstar / (1 + chi0s)),
)

# ---------------------------------------------------------------------------
# III. Canonical dependent-triple repair vector
# ---------------------------------------------------------------------------
subbanner("III. Exact canonical dependent-triple repair vector")

qvec = sp.Matrix([qtr, qnt, qeta])
dx_rep = sp.simplify(-Sdep * qvec)

print("Delta x_rep =")
sp.pprint(dx_rep)

expected_dx_rep = sp.Matrix(
    [
        [0],
        [0],
        [0],
        [0],
        [qeta],
        [0],
        [-qnt + qeta - Fstar * qtr / (1 + chi0s)],
        [-qtr / (1 + chi0s)],
    ]
)

expect_zero("repair vector - expected vector", simplify_log(dx_rep - expected_dx_rep))
expect_zero("M_* Delta x_rep + q", simplify_log(Mstar * dx_rep + qvec))

# ---------------------------------------------------------------------------
# IV. Exact fixed-free-quintuple uniqueness solve
# ---------------------------------------------------------------------------
subbanner("IV. Exact same-free-quintuple uniqueness")

dT, dKeta, dmu = sp.symbols("dT dKeta dmu", real=True)
dx_dep = sp.Matrix([0, 0, 0, 0, dKeta, 0, dmu, dT])

dep_equations = sp.simplify(Mstar * dx_dep + qvec)
print("Dependent-triple solve equations M_* dx_dep = -q:")
sp.pprint(dep_equations)

solution = sp.solve(
    [
        sp.Eq(dep_equations[0], 0),
        sp.Eq(dep_equations[1], 0),
        sp.Eq(dep_equations[2], 0),
    ],
    [dT, dKeta, dmu],
    dict=True,
)[0]

print("Unique dependent-triple solution =")
sp.pprint(solution)

expect_zero(
    "unique dT + q_tr/(1+chi0_*)",
    sp.simplify(solution[dT] + qtr / (1 + chi0s)),
)
expect_zero(
    "unique dKeta - q_eta",
    sp.simplify(solution[dKeta] - qeta),
)
expect_zero(
    "unique dmu - (-q_nt + q_eta - F_* q_tr/(1+chi0_*))",
    sp.simplify(solution[dmu] - (-qnt + qeta - Fstar * qtr / (1 + chi0s))),
)

# ---------------------------------------------------------------------------
# V. Exact canonical orbit projection in multiplicative variables
# ---------------------------------------------------------------------------
subbanner("V. Exact canonical orbit projection")

lam, cetaU, gamma, KU, Keta, KW, muW, TU = sp.symbols(
    "lambda_W c_etaU gamma K_U K_eta K_W mu_W T_U", positive=True, real=True
)

x = sp.Matrix([lam, cetaU, gamma, KU, Keta, KW, muW, TU])

x_proj = sp.Matrix(
    [
        lam,
        cetaU,
        gamma,
        KU,
        Keta * Reta,
        KW,
        muW * Rnt ** (-1) * Reta * Rtr ** (-Fstar / (1 + chi0s)),
        TU * Rtr ** (-sp.Rational(1, 1) / (1 + chi0s)),
    ]
)

print("x =")
sp.pprint(x)
print("Pi_O^can(x) =")
sp.pprint(x_proj)

# Compare the log-ratio x_proj / x to the canonical repair vector
dx_from_projection = sp.Matrix([sp.log(x_proj[i] / x[i]) for i in range(8)])
expect_zero("log(x_proj/x) - Delta x_rep", simplify_log(dx_from_projection - dx_rep))

# ---------------------------------------------------------------------------
# VI. Exact intrinsic packet equals the witness packet on the target orbit
# ---------------------------------------------------------------------------
subbanner("VI. Intrinsic packet equals pairwise witness packet")

Ctr_star, Cnt_star, epsEta_star = sp.symbols(
    "Ctr_star Cnt_star epsEta_star", positive=True, real=True
)
Ctr1, Cnt1, eps1 = sp.symbols("Ctr_1 Cnt_1 epsEta_1", positive=True, real=True)
Ctr2, Cnt2, eps2 = sp.symbols("Ctr_2 Cnt_2 epsEta_2", positive=True, real=True)

Delta_int = sp.Matrix([chiQ - 1, sp.log(Ctr2 / Ctr_star), sp.log(Cnt2 / Cnt_star), sp.log(eps2 / epsEta_star)])
Delta_pair = sp.Matrix([chiQ - 1, sp.log(Ctr2 / Ctr1), sp.log(Cnt2 / Cnt1), sp.log(eps2 / eps1)])

witness_subs = {Ctr1: Ctr_star, Cnt1: Cnt_star, eps1: epsEta_star}
expect_zero(
    "pairwise witness packet - intrinsic packet",
    simplify_log(Delta_pair.subs(witness_subs) - Delta_int),
)

# ---------------------------------------------------------------------------
# VII. Exact fixed-point criterion and linearized form
# ---------------------------------------------------------------------------
subbanner("VII. Fixed-point criterion and linearized realization compiler")

expect_zero("Pi^can(x) = x when R_tr=R_nt=R_eta=1", simplify_log(x_proj.subs({Rtr: 1, Rnt: 1, Reta: 1}) - x))

Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_W Delta_mu Delta_T",
    real=True,
)
Dx = sp.Matrix([Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT])
q_from_Dx = sp.simplify(Mstar * Dx)
dx_rep_lin = sp.simplify(-Sdep * q_from_Dx)

print("q^lin = M_* Delta x =")
sp.pprint(q_from_Dx)
print("Delta x_rep^lin = -S M_* Delta x =")
sp.pprint(dx_rep_lin)

expect_zero("M_*(Delta x + Delta x_rep^lin)", sp.simplify(Mstar * (Dx + dx_rep_lin)))

banner("STAGE 184 LEDGER")
print("1. The intrinsic target packet requires only four quantities from a candidate branch:")
print("      chi_Q(x), C_tr(x)/C_tr,*, C_nt(x)/C_nt,*, eps_eta(x)/eps_eta,*.")
print("2. The exact additive realization packet is")
print("      Delta_real^int = ( chi_Q-1 , ln R_tr , ln R_nt , ln R_eta ).")
print("3. The exact mismatch chart is")
print("      m_T = R_tr^(1/(1+chi0_*)),")
print("      m_K = R_eta^(-1),")
print("      m_mu = R_nt R_eta^(-1) R_tr^(F_*/(1+chi0_*)).")
print("4. The exact canonical dependent-triple repair vector is")
print("      Delta x_rep = -S q")
print("   and is supported only on")
print("      (T_U, K_eta^(eff), mu_W).")
print("5. The exact canonical orbit projection is")
print("      Pi_O^can(x) = (lambda_W, c_etaU, gamma, K_U, K_eta R_eta, K_W, mu_W R_nt^(-1) R_eta R_tr^(-F_*/(1+chi0_*)), T_U R_tr^(-1/(1+chi0_*))).")
print("6. It is the unique point on the target orbit with the same free quintuple")
print("      (lambda_W, c_etaU, gamma, K_U, K_W^(eff)).")
print("7. The explicit realization audit is therefore")
print("      x in Z_*  <=>  chi_Q(x)=1 and Pi_O^can(x)=x")
print("   equivalently")
print("      Delta_real^int(x|Z_*) = 0.")
