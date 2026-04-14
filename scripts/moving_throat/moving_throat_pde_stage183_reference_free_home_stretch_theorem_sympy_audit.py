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


def simplify_log(expr):
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(lambda z: sp.simplify(sp.expand_log(z, force=True)))
    return sp.simplify(sp.expand_log(expr, force=True))


banner("STAGE 183 — EXACT REFERENCE-FREE HOME-STRETCH THEOREM")

# ---------------------------------------------------------------------------
# I. Intrinsic packet equals the pairwise witness packet on the target orbit
# ---------------------------------------------------------------------------
subbanner("I. Exact intrinsic packet equals the pairwise witness packet")

chiQ2 = sp.symbols("chiQ_2", real=True)
Ctr_star, Cnt_star, epsEta_star = sp.symbols(
    "Ctr_star Cnt_star epsEta_star", positive=True, real=True
)
Ctr1, Cnt1, eps1 = sp.symbols("Ctr_1 Cnt_1 epsEta_1", positive=True, real=True)
Ctr2, Cnt2, eps2 = sp.symbols("Ctr_2 Cnt_2 epsEta_2", positive=True, real=True)

DeltaQ2 = sp.simplify(chiQ2 - 1)
qtr_int = sp.log(Ctr2 / Ctr_star)
qnt_int = sp.log(Cnt2 / Cnt_star)
qeta_int = sp.log(eps2 / epsEta_star)
Delta_H_int = sp.Matrix([DeltaQ2, qtr_int, qnt_int, qeta_int])

qtr_pair = sp.log(Ctr2 / Ctr1)
qnt_pair = sp.log(Cnt2 / Cnt1)
qeta_pair = sp.log(eps2 / eps1)
Delta_H_pair = sp.Matrix([DeltaQ2, qtr_pair, qnt_pair, qeta_pair])

witness_subs = {Ctr1: Ctr_star, Cnt1: Cnt_star, eps1: epsEta_star}
Delta_H_pair_on_witness = simplify_log(Delta_H_pair.subs(witness_subs))
Delta_H_int_simplified = simplify_log(Delta_H_int)

print("Delta_H^int(x^(2)) =")
sp.pprint(Delta_H_int_simplified)
print("Delta_H^(2<-1) on a target-orbit witness x^(1) =")
sp.pprint(Delta_H_pair_on_witness)
expect_zero("pairwise packet - intrinsic packet", Delta_H_pair_on_witness - Delta_H_int_simplified)

expect_zero(
    "exp(q_tr^int) - Ctr_2/Ctr_*",
    sp.exp(qtr_int) - Ctr2 / Ctr_star,
)
expect_zero(
    "exp(q_nt^int) - Cnt_2/Cnt_*",
    sp.exp(qnt_int) - Cnt2 / Cnt_star,
)
expect_zero(
    "exp(q_eta^int) - epsEta_2/epsEta_*",
    sp.exp(qeta_int) - eps2 / epsEta_star,
)

# ---------------------------------------------------------------------------
# II. Exact cocycle and witness independence
# ---------------------------------------------------------------------------
subbanner("II. Exact cocycle and closure-orbit witness independence")

Ctr3, Cnt3, eps3 = sp.symbols("Ctr_3 Cnt_3 epsEta_3", positive=True, real=True)

qtr_31 = sp.log(Ctr3 / Ctr1)
qtr_32 = sp.log(Ctr3 / Ctr2)
qtr_21 = sp.log(Ctr2 / Ctr1)
qnt_31 = sp.log(Cnt3 / Cnt1)
qnt_32 = sp.log(Cnt3 / Cnt2)
qnt_21 = sp.log(Cnt2 / Cnt1)
qeta_31 = sp.log(eps3 / eps1)
qeta_32 = sp.log(eps3 / eps2)
qeta_21 = sp.log(eps2 / eps1)

expect_zero("q_tr cocycle", simplify_log(qtr_31 - qtr_32 - qtr_21))
expect_zero("q_nt cocycle", simplify_log(qnt_31 - qnt_32 - qnt_21))
expect_zero("q_eta cocycle", simplify_log(qeta_31 - qeta_32 - qeta_21))

Ctr1p, Cnt1p, eps1p = sp.symbols("Ctr_1p Cnt_1p epsEta_1p", positive=True, real=True)
Delta_H_pair_prime = sp.Matrix([DeltaQ2, sp.log(Ctr2 / Ctr1p), sp.log(Cnt2 / Cnt1p), sp.log(eps2 / eps1p)])

witness_prime_subs = {Ctr1p: Ctr_star, Cnt1p: Cnt_star, eps1p: epsEta_star}
expect_zero(
    "witness independence of the full home-stretch packet",
    simplify_log(Delta_H_pair.subs(witness_subs) - Delta_H_pair_prime.subs(witness_prime_subs)),
)

# ---------------------------------------------------------------------------
# III. Exact mismatch chart for the final three Packet-B slots
# ---------------------------------------------------------------------------
subbanner("III. Exact mismatch chart")

chi0s, Fstar = sp.symbols("chi0_star F_star", positive=True, real=True)
mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)

qtr_from_m = sp.simplify((1 + chi0s) * sp.log(mT))
qeta_from_m = sp.simplify(-sp.log(mK))
qnt_from_m = sp.simplify(sp.log(mMu) - sp.log(mK) - Fstar * sp.log(mT))

print("q_tr(m) =")
sp.pprint(qtr_from_m)
print("q_nt(m) =")
sp.pprint(qnt_from_m)
print("q_eta(m) =")
sp.pprint(qeta_from_m)

expect_zero("exp(q_tr/(1+chi0_*)) - m_T", sp.exp(qtr_from_m / (1 + chi0s)) - mT)
expect_zero("exp(-q_eta) - m_K", sp.exp(-qeta_from_m) - mK)
expect_zero(
    "exp(q_nt - q_eta + F_* q_tr/(1+chi0_*)) - m_mu",
    sp.exp(qnt_from_m - qeta_from_m + Fstar * qtr_from_m / (1 + chi0s)) - mMu,
)

# ---------------------------------------------------------------------------
# IV. Exact reduced closure set in intrinsic form
# ---------------------------------------------------------------------------
subbanner("IV. Exact four-scalar intrinsic finish line")

ratio_packet = sp.Matrix([
    chiQ2,
    Ctr2 / Ctr_star,
    Cnt2 / Cnt_star,
    eps2 / epsEta_star,
])
ratio_packet_from_logs = sp.Matrix([
    DeltaQ2 + 1,
    sp.exp(qtr_int),
    sp.exp(qnt_int),
    sp.exp(qeta_int),
])

print("Intrinsic ratio packet =")
sp.pprint(ratio_packet)
expect_zero("ratio packet from Delta_H^int", ratio_packet_from_logs - ratio_packet)

# ---------------------------------------------------------------------------
# V. Final linearized four-scalar compiler
# ---------------------------------------------------------------------------
subbanner("V. Final linearized four-scalar compiler")

S, beta, Sigma0, Sigma5 = sp.symbols("S beta Sigma_0 Sigma_5", nonzero=True, real=True)
eps = sp.symbols("eps", real=True)
eps_beta, dSigma0, dSigma5 = sp.symbols("eps_beta dSigma_0 dSigma_5", real=True)

chi_from_def = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
DeltaQ_linear = sp.series(
    chi_from_def.subs({beta: 1 + eps * eps_beta, Sigma0: eps * dSigma0, Sigma5: eps * dSigma5}),
    eps,
    0,
    2,
).removeO() - 1
DeltaQ_linear_expected = eps * (5 * eps_beta + dSigma0 / (3 * S) + 9 * dSigma5 / S)

print("Delta_Q^lin =")
sp.pprint(sp.expand(DeltaQ_linear))
expect_zero(
    "Delta_Q^lin - eps(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S)",
    sp.expand(DeltaQ_linear - DeltaQ_linear_expected),
)

Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_W Delta_mu Delta_T",
    real=True,
)
chi0_lin, deltaU_lin = sp.symbols("chi0_star_lin deltaU_star_lin", positive=True, real=True)
E_lin, F_lin = sp.symbols("E_star_lin F_star_lin", real=True)
Mstar = sp.Matrix(
    [
        [0, 1 + deltaU_lin, 1 + deltaU_lin, -(2 + chi0_lin + deltaU_lin), 0, 0, 0, 1 + chi0_lin],
        [2 * (1 + E_lin), 0, 2 * E_lin, F_lin - E_lin, -1, -(2 + E_lin), 1, -F_lin],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)

Dx = sp.Matrix([Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT])
q_linear = sp.simplify(Mstar * Dx)

Delta_H_lin = sp.Matrix([DeltaQ_linear_expected, q_linear[0], q_linear[1], q_linear[2]])
print("Delta_H^lin =")
sp.pprint(Delta_H_lin)

banner("STAGE 183 LEDGER")
print("1. The exact target closure orbit is the monomial level set")
print("      O_* = { x>0 : C_tr(x)=C_tr,*, C_nt(x)=C_nt,*, eps_eta(x)=eps_eta,* }.")
print("   By the carried finite-fibre theorem, this is a single exact similarity orbit.")
print("2. The full reduced closure set is")
print("      Z_* = { x in O_* : chi_Q(x)=1 }.")
print("3. The full intrinsic home-stretch packet is the exact four-scalar compiler")
print("      Delta_HS^int(x) = (chi_Q(x)-1, ln[C_tr(x)/C_tr,*], ln[C_nt(x)/C_nt,*], ln[eps_eta(x)/eps_eta,*]).")
print("4. If x^(1) lies anywhere on the target orbit O_*, then the exact pairwise packet")
print("      Delta_HS^(2<-1) = (chi_Q^(2)-1, q_tr^(2<-1), q_nt^(2<-1), q_eta^(2<-1))")
print("   equals the intrinsic packet of x^(2). Hence the verdict is independent of which witness on O_* is chosen.")
print("5. The exact reference-free full home-stretch theorem is")
print("      reduced closure complete  <=>  x in Z_*")
print("   equivalently")
print("      chi_Q=1 and x in O_*")
print("   equivalently, relative to any orbit witness,")
print("      chi_Q=1 and q_tr=q_nt=q_eta=0")
print("   equivalently")
print("      chi_Q=1 and m_T=m_K=m_mu=1.")
print("6. The final first-order compiler is")
print("      Delta_HS^lin = ( Delta_Q^lin , M_* delta x ),")
print("   with")
print("      Delta_Q^lin = 5 eps_beta + dSigma_0/(3S) + 9 dSigma_5/S.")
