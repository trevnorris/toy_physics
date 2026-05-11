#!/usr/bin/env python3
"""
Stage 200 SymPy audit.

Provenance notes
----------------
- `chi0s`, `deltaUs`, `Estar`, and `Fstar` are the carried Stage 192 Packet-B
  exponents. The thermal ratio `(pi^2 T)/(L^2 K_U)` remains explicit in the
  monomials; `deltaUs` is not a silent replacement for that factor.
- The same four-scalar compiler is reused downstream in Stages 185--186, so the
  log-ratio matrix is kept in the Stage 192 notation without introducing a new
  alphabet here.
"""

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


def simplify_expr(expr):
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(
            lambda z: sp.simplify(sp.expand_log(sp.expand_power_exp(z), force=True))
        )
    return sp.simplify(sp.expand_log(sp.expand_power_exp(expr), force=True))


def expect_zero(name: str, expr) -> None:
    expr = simplify_expr(expr)
    if isinstance(expr, sp.MatrixBase):
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


def ctr_monomial(gamma, c_eta_u, k_u, t_u, chi0s, delta_us, L):
    # Stage 192 control monomial: the thermal ratio pi^2 T/(L^2 K_U) stays
    # explicit; `delta_us` is only the exponent on the gamma*c_eta/K_U factor.
    return sp.simplify(
        (gamma * c_eta_u / k_u) ** (1 + delta_us)
        * (sp.pi**2 * t_u / (L**2 * k_u)) ** (1 + chi0s)
    )


def cnt_monomial(lambda_w, gamma, k_u, k_eta, k_w, mu_w, t_u, Estar, Fstar, L, sigma):
    return sp.simplify(
        (lambda_w**2 * mu_w / (k_eta * k_w**2))
        * (gamma**2 * lambda_w**2 * sigma / (k_u * k_w)) ** Estar
        * (sp.pi**2 * t_u / (L**2 * k_u)) ** (-Fstar)
    )


def eps_eta_monomial(c_eta_u, k_u, k_eta):
    return sp.simplify(c_eta_u**2 / (k_u * k_eta))


banner("STAGE 183 — EXACT REFERENCE-FREE HOME-STRETCH THEOREM")

chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)
L, sigma = sp.symbols("L sigma", positive=True, real=True)

# ---------------------------------------------------------------------------
# I. Derive the exact Packet-B compiler from the primitive monomial ratios
# ---------------------------------------------------------------------------
subbanner("I. Exact four-scalar Packet-B compiler from primitive monomial ratios")

lam1, c1, gam1 = sp.symbols("lambda1 cEtaU1 gamma1", positive=True, real=True)
KU1, Keta1, KW1, mu1, T1 = sp.symbols(
    "KU1 Keta1 KW1 mu1 T1", positive=True, real=True
)
lam2, c2, gam2 = sp.symbols("lambda2 cEtaU2 gamma2", positive=True, real=True)
KU2, Keta2, KW2, mu2, T2 = sp.symbols(
    "KU2 Keta2 KW2 mu2 T2", positive=True, real=True
)
rla, rc, rg, rU, rK, rW, rmu, rT = sp.symbols(
    "r_lambda r_c r_gamma r_U r_K r_W r_mu r_T", positive=True, real=True
)

Ctr1 = ctr_monomial(gam1, c1, KU1, T1, chi0s, deltaUs, L)
Ctr2 = ctr_monomial(gam2, c2, KU2, T2, chi0s, deltaUs, L)
Cnt1 = cnt_monomial(lam1, gam1, KU1, Keta1, KW1, mu1, T1, Estar, Fstar, L, sigma)
Cnt2 = cnt_monomial(lam2, gam2, KU2, Keta2, KW2, mu2, T2, Estar, Fstar, L, sigma)
eps1 = eps_eta_monomial(c1, KU1, Keta1)
eps2 = eps_eta_monomial(c2, KU2, Keta2)

ratio_subs = {
    lam2: rla * lam1,
    c2: rc * c1,
    gam2: rg * gam1,
    KU2: rU * KU1,
    Keta2: rK * Keta1,
    KW2: rW * KW1,
    mu2: rmu * mu1,
    T2: rT * T1,
}

Ctr_ratio = simplify_expr((Ctr2 / Ctr1).subs(ratio_subs))
Cnt_ratio = simplify_expr((Cnt2 / Cnt1).subs(ratio_subs))
eps_ratio = simplify_expr((eps2 / eps1).subs(ratio_subs))

print("Ctr_2 / Ctr_1 =")
sp.pprint(Ctr_ratio)
print("Cnt_2 / Cnt_1 =")
sp.pprint(Cnt_ratio)
print("epsEta_2 / epsEta_1 =")
sp.pprint(eps_ratio)

Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_W Delta_mu Delta_T",
    real=True,
)
ratio_to_logs = {
    rla: sp.exp(Dl),
    rc: sp.exp(Dc),
    rg: sp.exp(Dg),
    rU: sp.exp(DU),
    rK: sp.exp(DKeta),
    rW: sp.exp(DW),
    rmu: sp.exp(Dmu),
    rT: sp.exp(DT),
}

q_pair = simplify_expr(
    sp.Matrix(
        [
            sp.log(Ctr_ratio.subs(ratio_to_logs)),
            sp.log(Cnt_ratio.subs(ratio_to_logs)),
            sp.log(eps_ratio.subs(ratio_to_logs)),
        ]
    )
)
Dvec = sp.Matrix([Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT])
Mderived = simplify_expr(q_pair.jacobian(Dvec))
Mexpected = sp.Matrix(
    [
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)

print("q^(2<-1) from primitive monomial ratios =")
sp.pprint(q_pair)
print("M_* derived from the monomial log-ratio Jacobian =")
sp.pprint(Mderived)
expect_zero("derived M_* - carried Stage 192 matrix", Mderived - Mexpected)
expect_zero("q^(2<-1) - M_* Delta x^(2<-1)", q_pair - Mderived * Dvec)

# ---------------------------------------------------------------------------
# II. Exact witness invariance from the finite similarity-orbit law
# ---------------------------------------------------------------------------
subbanner("II. Exact target-orbit witness invariance from the carried orbit law")

lams, cs, gs = sp.symbols("lambda_star cEtaU_star gamma_star", positive=True, real=True)
KUs, Ketas, KWs, mus, Ts = sp.symbols(
    "KU_star Keta_star KW_star mu_star T_star", positive=True, real=True
)
rla_w, rc_w, rg_w, rU_w, rW_w = sp.symbols(
    "r_lambda_w r_c_w r_gamma_w r_U_w r_W_w", positive=True, real=True
)
lamx, cx, gx = sp.symbols("lambda_x cEtaU_x gamma_x", positive=True, real=True)
KUx, Ketax, KWx, mux, Tx = sp.symbols(
    "KU_x Keta_x KW_x mu_x T_x", positive=True, real=True
)
chiQx = sp.symbols("chiQ_x", real=True)

alpha_star = sp.simplify((1 + deltaUs) / (1 + chi0s))
PhiT_w = sp.simplify(rU_w * (rU_w / (rg_w * rc_w)) ** alpha_star)
PhiK_w = sp.simplify(rc_w**2 / rU_w)
PhiMu_w = sp.simplify(
    PhiK_w
    * rW_w**2
    / rla_w**2
    * (rg_w**2 * rla_w**2 / (rU_w * rW_w)) ** (-Estar)
    * (PhiT_w / rU_w) ** Fstar
)

Ctr_star = ctr_monomial(gs, cs, KUs, Ts, chi0s, deltaUs, L)
Cnt_star = cnt_monomial(lams, gs, KUs, Ketas, KWs, mus, Ts, Estar, Fstar, L, sigma)
epsEta_star = eps_eta_monomial(cs, KUs, Ketas)

Ctr_witness = ctr_monomial(
    rg_w * gs,
    rc_w * cs,
    rU_w * KUs,
    PhiT_w * Ts,
    chi0s,
    deltaUs,
    L,
)
Cnt_witness = cnt_monomial(
    rla_w * lams,
    rg_w * gs,
    rU_w * KUs,
    PhiK_w * Ketas,
    rW_w * KWs,
    PhiMu_w * mus,
    PhiT_w * Ts,
    Estar,
    Fstar,
    L,
    sigma,
)
epsEta_witness = eps_eta_monomial(rc_w * cs, rU_w * KUs, PhiK_w * Ketas)

expect_zero(
    "ln[Ctr(witness) / Ctr(*)]",
    sp.expand_log(sp.log(Ctr_witness / Ctr_star), force=True),
)
expect_zero(
    "ln[Cnt(witness) / Cnt(*)]",
    sp.expand_log(sp.log(Cnt_witness / Cnt_star), force=True),
)
expect_zero(
    "ln[epsEta(witness) / epsEta(*)]",
    sp.expand_log(sp.log(epsEta_witness / epsEta_star), force=True),
)

Ctr_x = ctr_monomial(gx, cx, KUx, Tx, chi0s, deltaUs, L)
Cnt_x = cnt_monomial(lamx, gx, KUx, Ketax, KWx, mux, Tx, Estar, Fstar, L, sigma)
epsEta_x = eps_eta_monomial(cx, KUx, Ketax)

Delta_H_intrinsic = simplify_expr(
    sp.Matrix(
        [
            chiQx - 1,
            sp.log(Ctr_x / Ctr_star),
            sp.log(Cnt_x / Cnt_star),
            sp.log(epsEta_x / epsEta_star),
        ]
    )
)
Delta_H_pair_witness = simplify_expr(
    sp.Matrix(
        [
            chiQx - 1,
            sp.log(Ctr_x / Ctr_witness),
            sp.log(Cnt_x / Cnt_witness),
            sp.log(epsEta_x / epsEta_witness),
        ]
    )
)

print("Delta_H^int(x | O_*) =")
sp.pprint(Delta_H_intrinsic)
print("Delta_H^(x<-w) with an arbitrary target-orbit witness w =")
sp.pprint(Delta_H_pair_witness)
expect_zero(
    "pairwise witness packet - intrinsic orbit packet",
    Delta_H_pair_witness - Delta_H_intrinsic,
)

# ---------------------------------------------------------------------------
# III. Exact mismatch chart from the carried dependent-triple orbit solve
# ---------------------------------------------------------------------------
subbanner("III. Exact mismatch chart from the dependent-triple orbit solve")

lamf, cf, gf = sp.symbols("lambda_f cEtaU_f gamma_f", positive=True, real=True)
KUf, KWf = sp.symbols("KU_f KW_f", positive=True, real=True)
Ctr_target, Cnt_target, epsEta_target = sp.symbols(
    "Ctr_target Cnt_target epsEta_target", positive=True, real=True
)
mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)

Keta_orbit = sp.simplify(cf**2 / (KUf * epsEta_target))
T_orbit = sp.simplify(
    (L**2 * KUf / sp.pi**2)
    * (Ctr_target / (gf * cf / KUf) ** (1 + deltaUs)) ** (1 / (1 + chi0s))
)
mu_orbit = sp.simplify(
    Cnt_target
    * Keta_orbit
    * KWf**2
    / lamf**2
    * (gf**2 * lamf**2 * sigma / (KUf * KWf)) ** (-Estar)
    * (sp.pi**2 * T_orbit / (L**2 * KUf)) ** Fstar
)

T_actual = sp.simplify(mT * T_orbit)
Keta_actual = sp.simplify(mK * Keta_orbit)
mu_actual = sp.simplify(mMu * mu_orbit)

Ctr_actual_ratio = simplify_expr(
    ctr_monomial(gf, cf, KUf, T_actual, chi0s, deltaUs, L) / Ctr_target
)
Cnt_actual_ratio = simplify_expr(
    cnt_monomial(lamf, gf, KUf, Keta_actual, KWf, mu_actual, T_actual, Estar, Fstar, L, sigma)
    / Cnt_target
)
epsEta_actual_ratio = simplify_expr(eps_eta_monomial(cf, KUf, Keta_actual) / epsEta_target)

expect_zero("Ctr(actual) / Ctr_* - m_T^(1+chi0_*)", Ctr_actual_ratio - mT ** (1 + chi0s))
expect_zero("Cnt(actual) / Cnt_* - m_mu/(m_K m_T^F_*)", Cnt_actual_ratio - mMu / (mK * mT**Fstar))
expect_zero("epsEta(actual) / epsEta_* - 1/m_K", epsEta_actual_ratio - 1 / mK)

q_mismatch = simplify_expr(
    sp.Matrix(
        [
            sp.log(Ctr_actual_ratio),
            sp.log(Cnt_actual_ratio),
            sp.log(epsEta_actual_ratio),
        ]
    )
)
Dmis = sp.Matrix([0, 0, 0, 0, sp.log(mK), 0, sp.log(mMu), sp.log(mT)])
q_mismatch_expected = simplify_expr(
    sp.Matrix(
        [
            (1 + chi0s) * sp.log(mT),
            sp.log(mMu) - sp.log(mK) - Fstar * sp.log(mT),
            -sp.log(mK),
        ]
    )
)

print("q(mismatch) =")
sp.pprint(q_mismatch)
expect_zero("q(mismatch) - carried chart", q_mismatch - q_mismatch_expected)
expect_zero("M_* Delta x_mis - q(mismatch)", Mderived * Dmis - q_mismatch)

# ---------------------------------------------------------------------------
# IV. Exact cocycle law from microscopic ratio composition
# ---------------------------------------------------------------------------
subbanner("IV. Exact cocycle law from microscopic ratio composition")

rla21, rc21, rg21, rU21, rK21, rW21, rmu21, rT21 = sp.symbols(
    "rla21 rc21 rg21 rU21 rK21 rW21 rmu21 rT21", positive=True, real=True
)
rla32, rc32, rg32, rU32, rK32, rW32, rmu32, rT32 = sp.symbols(
    "rla32 rc32 rg32 rU32 rK32 rW32 rmu32 rT32", positive=True, real=True
)

D21 = sp.Matrix(
    [
        sp.log(rla21),
        sp.log(rc21),
        sp.log(rg21),
        sp.log(rU21),
        sp.log(rK21),
        sp.log(rW21),
        sp.log(rmu21),
        sp.log(rT21),
    ]
)
D32 = sp.Matrix(
    [
        sp.log(rla32),
        sp.log(rc32),
        sp.log(rg32),
        sp.log(rU32),
        sp.log(rK32),
        sp.log(rW32),
        sp.log(rmu32),
        sp.log(rT32),
    ]
)
D31 = simplify_expr(
    sp.Matrix(
        [
            sp.log(rla32 * rla21),
            sp.log(rc32 * rc21),
            sp.log(rg32 * rg21),
            sp.log(rU32 * rU21),
            sp.log(rK32 * rK21),
            sp.log(rW32 * rW21),
            sp.log(rmu32 * rmu21),
            sp.log(rT32 * rT21),
        ]
    )
)
q21 = simplify_expr(Mderived * D21)
q32 = simplify_expr(Mderived * D32)
q31 = simplify_expr(Mderived * D31)

expect_zero("Delta x^31 - Delta x^32 - Delta x^21", D31 - D32 - D21)
expect_zero("q^(31) - q^(32) - q^(21)", q31 - q32 - q21)

# ---------------------------------------------------------------------------
# V. Packet-A linearization and the assembled full four-scalar compiler
# ---------------------------------------------------------------------------
subbanner("V. Exact Packet-A linearization and assembled full compiler")

S, beta, Sigma0, Sigma5 = sp.symbols("S beta Sigma_0 Sigma_5", nonzero=True, real=True)
eps = sp.symbols("eps", real=True)
eps_beta, dSigma0, dSigma5 = sp.symbols("eps_beta dSigma_0 dSigma_5", real=True)

chi_from_def = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
DeltaQ_linear = sp.series(
    chi_from_def.subs(
        {
            beta: 1 + eps * eps_beta,
            Sigma0: eps * dSigma0,
            Sigma5: eps * dSigma5,
        }
    ),
    eps,
    0,
    2,
).removeO() - 1
DeltaQ_linear_expected = eps * (5 * eps_beta + dSigma0 / (3 * S) + 9 * dSigma5 / S)

expect_zero(
    "Delta_Q^lin - eps(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S)",
    DeltaQ_linear - DeltaQ_linear_expected,
)

Delta_H_linear = simplify_expr(
    sp.Matrix(
        [
            DeltaQ_linear_expected,
            *(Mderived * Dvec),
        ]
    )
)
print("Delta_H^lin =")
sp.pprint(Delta_H_linear)

banner("STAGE 183 LEDGER")
print("1. The Packet-B quotient coordinates are derived from the primitive coherent monomials")
print("      (C_tr, C_nt, eps_eta)")
print("   and their logarithmic Jacobian reproduces the carried Stage 192 matrix M_*.")
print("2. The target-orbit witness packet is genuinely witness-independent because the")
print("   finite similarity-orbit law preserves the three coherent monomials exactly.")
print("3. The mismatch chart q = ((1+chi0_*) ln m_T, ln m_mu - ln m_K - F_* ln m_T, -ln m_K)")
print("   is re-derived from the exact dependent-triple orbit solve rather than posited.")
print("4. The pairwise cocycle law is inherited from microscopic ratio composition")
print("      Delta x^(31) = Delta x^(32) + Delta x^(21),")
print("   followed by the monomial compiler q = M_* Delta x.")
print("5. The full linearized home-stretch compiler is")
print("      Delta_HS^lin = (Delta_Q^lin, M_* delta x),")
print("   where Delta_Q^lin is carried from the exact Packet-A closure law and the")
print("   remaining three rows come from the monomial log-ratio Jacobian.")
