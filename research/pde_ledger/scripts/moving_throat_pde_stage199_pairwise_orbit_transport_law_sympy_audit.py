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


def simplify_expr(expr):
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(lambda z: sp.simplify(sp.expand_log(sp.expand_power_exp(z), force=True)))
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


banner("STAGE 182 — EXACT PAIRWISE ORBIT-TRANSPORT LAW AND THE TWO-POINT ORBIT-LOCK THEOREM")

# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------
chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)
L, sigma = sp.symbols("L sigma", positive=True, real=True)
pi = sp.pi

# State 1 and State 2 microscopic coordinates
lam1, c1, gam1 = sp.symbols("lambda1 cEtaU1 gamma1", positive=True, real=True)
KU1, Keta1, KW1, mu1, T1 = sp.symbols(
    "KU1 Keta1 KW1 mu1 T1", positive=True, real=True
)
lam2, c2, gam2 = sp.symbols("lambda2 cEtaU2 gamma2", positive=True, real=True)
KU2, Keta2, KW2, mu2, T2 = sp.symbols(
    "KU2 Keta2 KW2 mu2 T2", positive=True, real=True
)

# Pairwise ratios
rla, rc, rg, rU, rK, rW, rmu, rT = sp.symbols(
    "r_lambda r_c r_gamma r_U r_K r_W r_mu r_T", positive=True, real=True
)

# ---------------------------------------------------------------------------
# Exact coherent monomials and pairwise ratio formulas
# ---------------------------------------------------------------------------
subbanner("I. Exact pairwise monomial ratios")

Ctr1 = sp.simplify(
    (gam1 * c1 / KU1) ** (1 + deltaUs) * (pi**2 * T1 / (L**2 * KU1)) ** (1 + chi0s)
)
Ctr2 = sp.simplify(
    (gam2 * c2 / KU2) ** (1 + deltaUs) * (pi**2 * T2 / (L**2 * KU2)) ** (1 + chi0s)
)
Cnt1 = sp.simplify(
    (lam1**2 * mu1 / (Keta1 * KW1**2))
    * (gam1**2 * lam1**2 * sigma / (KU1 * KW1)) ** Estar
    * (pi**2 * T1 / (L**2 * KU1)) ** (-Fstar)
)
Cnt2 = sp.simplify(
    (lam2**2 * mu2 / (Keta2 * KW2**2))
    * (gam2**2 * lam2**2 * sigma / (KU2 * KW2)) ** Estar
    * (pi**2 * T2 / (L**2 * KU2)) ** (-Fstar)
)
eps1 = sp.simplify(c1**2 / (KU1 * Keta1))
eps2 = sp.simplify(c2**2 / (KU2 * Keta2))

Ctr_ratio = sp.simplify(Ctr2 / Ctr1)
Cnt_ratio = sp.simplify(Cnt2 / Cnt1)
eps_ratio = sp.simplify(eps2 / eps1)

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
Ctr_ratio_pairs = sp.simplify(Ctr_ratio.subs(ratio_subs))
Cnt_ratio_pairs = sp.simplify(Cnt_ratio.subs(ratio_subs))
eps_ratio_pairs = sp.simplify(eps_ratio.subs(ratio_subs))

Ctr_ratio_expected = sp.simplify((rg * rc / rU) ** (1 + deltaUs) * (rT / rU) ** (1 + chi0s))
Cnt_ratio_expected = sp.simplify(
    (rla**2 * rmu / (rK * rW**2))
    * (rg**2 * rla**2 / (rU * rW)) ** Estar
    * (rT / rU) ** (-Fstar)
)
eps_ratio_expected = sp.simplify(rc**2 / (rU * rK))

print("Ctr_2 / Ctr_1 =")
sp.pprint(Ctr_ratio_pairs)
print("Cnt_2 / Cnt_1 =")
sp.pprint(Cnt_ratio_pairs)
print("epsEta_2 / epsEta_1 =")
sp.pprint(eps_ratio_pairs)

expect_zero("Ctr pairwise ratio - expected", Ctr_ratio_pairs - Ctr_ratio_expected)
expect_zero("Cnt pairwise ratio - expected", Cnt_ratio_pairs - Cnt_ratio_expected)
expect_zero("epsEta pairwise ratio - expected", eps_ratio_pairs - eps_ratio_expected)

# ---------------------------------------------------------------------------
# Exact pairwise orbit transport factors
# ---------------------------------------------------------------------------
subbanner("II. Exact pairwise orbit transport factors")

alpha_star = sp.simplify((1 + deltaUs) / (1 + chi0s))
PhiT = sp.simplify(rU * (rU / (rg * rc)) ** alpha_star)
PhiK = sp.simplify(rc**2 / rU)
PhiMu = sp.simplify(
    PhiK * rW**2 / rla**2
    * (rg**2 * rla**2 / (rU * rW)) ** (-Estar)
    * (PhiT / rU) ** Fstar
)
PhiMu_expanded = sp.simplify(
    rc ** (2 - Fstar * alpha_star)
    * rg ** (-2 * Estar - Fstar * alpha_star)
    * rla ** (-2 * (1 + Estar))
    * rU ** (-1 + Estar + Fstar * alpha_star)
    * rW ** (2 + Estar)
)

print("Phi_T =")
sp.pprint(PhiT)
print("Phi_K =")
sp.pprint(PhiK)
print("Phi_mu =")
sp.pprint(PhiMu)
expect_zero(
    "ln Phi_mu - ln expanded monomial form",
    sp.expand_log(sp.log(PhiMu) - sp.log(PhiMu_expanded), force=True),
)

Ctr_same_orbit = sp.simplify(Ctr_ratio_expected.subs({rT: PhiT}))
eps_same_orbit = sp.simplify(eps_ratio_expected.subs({rK: PhiK}))
Cnt_same_orbit = sp.simplify(Cnt_ratio_expected.subs({rT: PhiT, rK: PhiK, rmu: PhiMu}))

expect_zero(
    "ln[Ctr ratio on same orbit]",
    sp.expand_log(sp.log(Ctr_same_orbit), force=True),
)
expect_zero(
    "ln[epsEta ratio on same orbit]",
    sp.expand_log(sp.log(eps_same_orbit), force=True),
)
expect_zero(
    "ln[Cnt ratio on same orbit]",
    sp.expand_log(sp.log(Cnt_same_orbit), force=True),
)

# ---------------------------------------------------------------------------
# Exact mismatch packet and pairwise invariant-ratio collapse
# ---------------------------------------------------------------------------
subbanner("III. Exact reference-independent mismatch triple")

mT = sp.simplify(rT / PhiT)
mK = sp.simplify(rK / PhiK)
mMu = sp.simplify(rmu / PhiMu)

print("m_T =")
sp.pprint(mT)
print("m_K =")
sp.pprint(mK)
print("m_mu =")
sp.pprint(mMu)

expect_zero(
    "Ctr pairwise ratio - m_T^(1+chi0_*)",
    sp.expand_log(sp.log(Ctr_ratio_expected) - (1 + chi0s) * sp.log(mT), force=True),
)
expect_zero(
    "epsEta pairwise ratio - 1/m_K",
    sp.expand_log(sp.log(eps_ratio_expected) - sp.log(1 / mK), force=True),
)
expect_zero(
    "Cnt pairwise ratio - m_mu/(m_K m_T^F_*)",
    sp.expand_log(sp.log(Cnt_ratio_expected) - (sp.log(mMu) - sp.log(mK) - Fstar * sp.log(mT)), force=True),
)

# ---------------------------------------------------------------------------
# Exact q-chart and projector split on the finite pairwise log-ratio vector
# ---------------------------------------------------------------------------
subbanner("IV. Exact q-chart and Stage 192 projectors on finite two-point data")

# Ordered basis: (lambda, c, gamma, U, K_eta, W, mu, T)
Mstar = sp.Matrix(
    [
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)

Dpair = sp.Matrix([
    sp.log(rla), sp.log(rc), sp.log(rg), sp.log(rU), sp.log(rK), sp.log(rW), sp.log(rmu), sp.log(rT)
])
qpair = simplify_expr(Mstar * Dpair)
qtr, qnt, qeta = qpair

print("q^(2<-1) = M_* Delta x^(2<-1) =")
sp.pprint(qpair)
expect_zero(
    "q_tr - (1+chi0_*) ln m_T",
    qtr - (1 + chi0s) * sp.log(mT),
)
expect_zero(
    "q_eta + ln m_K",
    qeta + sp.log(mK),
)
expect_zero(
    "q_nt - [ln m_mu - ln m_K - F_* ln m_T]",
    qnt - (sp.log(mMu) - sp.log(mK) - Fstar * sp.log(mT)),
)

# Stage 192 exact projectors
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
Qquot = sp.simplify(Sdep * Mstar)
Oorb = sp.simplify(sp.eye(8) - Qquot)

Qpair = simplify_expr(Qquot * Dpair)
Opair = simplify_expr(Oorb * Dpair)

print("Q_quot Delta x^(2<-1) =")
sp.pprint(Qpair)
print("O_orb Delta x^(2<-1) =")
sp.pprint(Opair)

expected_Qpair = sp.Matrix([
    0,
    0,
    0,
    0,
    sp.log(mK),
    0,
    sp.log(mMu),
    sp.log(mT),
])
expected_Opair = sp.Matrix([
    sp.log(rla),
    sp.log(rc),
    sp.log(rg),
    sp.log(rU),
    sp.log(PhiK),
    sp.log(rW),
    sp.log(PhiMu),
    sp.log(PhiT),
])

expect_zero("Q_quot Delta x - expected", Qpair - expected_Qpair)
expect_zero("O_orb Delta x - expected", Opair - expected_Opair)
expect_zero("Q + O - Delta x", Qpair + Opair - Dpair)
expect_zero("M_* O_orb Delta x", Mstar * Opair)
expect_zero("M_* Q_quot Delta x - qpair", Mstar * Qpair - qpair)

# ---------------------------------------------------------------------------
# Exact pairwise restoration map
# ---------------------------------------------------------------------------
subbanner("V. Exact pairwise restoration map")

T2_restore = simplify_expr(T2 * sp.exp(-qtr / (1 + chi0s))).subs(ratio_subs)
K2_restore = simplify_expr(Keta2 * sp.exp(qeta)).subs(ratio_subs)
mu2_restore = simplify_expr(mu2 * sp.exp(-qnt + qeta - Fstar * qtr / (1 + chi0s))).subs(ratio_subs)

expect_zero(
    "ln T2_restore - ln(Phi_T T1)",
    sp.expand_log(sp.log(T2_restore) - sp.log(PhiT * T1), force=True),
)
expect_zero(
    "ln K2_restore - ln(Phi_K Keta1)",
    sp.expand_log(sp.log(K2_restore) - sp.log(PhiK * Keta1), force=True),
)
expect_zero(
    "ln mu2_restore - ln(Phi_mu mu1)",
    sp.expand_log(sp.log(mu2_restore) - sp.log(PhiMu * mu1), force=True),
)

# ---------------------------------------------------------------------------
# Exact composition laws
# ---------------------------------------------------------------------------
subbanner("VI. Exact cocycle/composition laws")

rla21, rc21, rg21, rU21, rK21, rW21, rmu21, rT21 = sp.symbols(
    "rla21 rc21 rg21 rU21 rK21 rW21 rmu21 rT21", positive=True, real=True
)
rla32, rc32, rg32, rU32, rK32, rW32, rmu32, rT32 = sp.symbols(
    "rla32 rc32 rg32 rU32 rK32 rW32 rmu32 rT32", positive=True, real=True
)

def PhiT_of(rl, rc_, rg_, rU_, rW_):
    return sp.simplify(rU_ * (rU_ / (rg_ * rc_)) ** alpha_star)


def PhiK_of(rl, rc_, rg_, rU_, rW_):
    return sp.simplify(rc_**2 / rU_)


def PhiMu_of(rl, rc_, rg_, rU_, rW_):
    PhiT_local = PhiT_of(rl, rc_, rg_, rU_, rW_)
    PhiK_local = PhiK_of(rl, rc_, rg_, rU_, rW_)
    return sp.simplify(
        PhiK_local * rW_**2 / rl**2
        * (rg_**2 * rl**2 / (rU_ * rW_)) ** (-Estar)
        * (PhiT_local / rU_) ** Fstar
    )

PhiT21 = PhiT_of(rla21, rc21, rg21, rU21, rW21)
PhiK21 = PhiK_of(rla21, rc21, rg21, rU21, rW21)
PhiMu21 = PhiMu_of(rla21, rc21, rg21, rU21, rW21)
PhiT32 = PhiT_of(rla32, rc32, rg32, rU32, rW32)
PhiK32 = PhiK_of(rla32, rc32, rg32, rU32, rW32)
PhiMu32 = PhiMu_of(rla32, rc32, rg32, rU32, rW32)
PhiT31 = PhiT_of(rla32 * rla21, rc32 * rc21, rg32 * rg21, rU32 * rU21, rW32 * rW21)
PhiK31 = PhiK_of(rla32 * rla21, rc32 * rc21, rg32 * rg21, rU32 * rU21, rW32 * rW21)
PhiMu31 = PhiMu_of(rla32 * rla21, rc32 * rc21, rg32 * rg21, rU32 * rU21, rW32 * rW21)

expect_zero("Phi_T^31 - Phi_T^32 Phi_T^21", PhiT31 - PhiT32 * PhiT21)
expect_zero("Phi_K^31 - Phi_K^32 Phi_K^21", PhiK31 - PhiK32 * PhiK21)
expect_zero("Phi_mu^31 - Phi_mu^32 Phi_mu^21", PhiMu31 - PhiMu32 * PhiMu21)

mT21 = rT21 / PhiT21
mT32 = rT32 / PhiT32
mT31 = (rT32 * rT21) / PhiT31
mK21 = rK21 / PhiK21
mK32 = rK32 / PhiK32
mK31 = (rK32 * rK21) / PhiK31
mMu21 = rmu21 / PhiMu21
mMu32 = rmu32 / PhiMu32
mMu31 = (rmu32 * rmu21) / PhiMu31

expect_zero(
    "ln m_T^31 - ln m_T^32 - ln m_T^21",
    sp.expand_log(sp.log(mT31) - sp.log(mT32) - sp.log(mT21), force=True),
)
expect_zero(
    "ln m_K^31 - ln m_K^32 - ln m_K^21",
    sp.expand_log(sp.log(mK31) - sp.log(mK32) - sp.log(mK21), force=True),
)
expect_zero(
    "ln m_mu^31 - ln m_mu^32 - ln m_mu^21",
    sp.expand_log(sp.log(mMu31) - sp.log(mMu32) - sp.log(mMu21), force=True),
)

D21 = sp.Matrix([
    sp.log(rla21), sp.log(rc21), sp.log(rg21), sp.log(rU21), sp.log(rK21), sp.log(rW21), sp.log(rmu21), sp.log(rT21)
])
D32 = sp.Matrix([
    sp.log(rla32), sp.log(rc32), sp.log(rg32), sp.log(rU32), sp.log(rK32), sp.log(rW32), sp.log(rmu32), sp.log(rT32)
])
D31 = sp.Matrix([
    sp.log(rla32 * rla21), sp.log(rc32 * rc21), sp.log(rg32 * rg21), sp.log(rU32 * rU21),
    sp.log(rK32 * rK21), sp.log(rW32 * rW21), sp.log(rmu32 * rmu21), sp.log(rT32 * rT21)
])
expect_zero("Delta x^31 - Delta x^32 - Delta x^21", D31 - D32 - D21)
expect_zero("q^31 - q^32 - q^21", Mstar * D31 - Mstar * D32 - Mstar * D21)

# ---------------------------------------------------------------------------
# Reduction to Stage 198
# ---------------------------------------------------------------------------
subbanner("VII. Exact reduction to Stage 198")

special_free_equal = {rla: 1, rc: 1, rg: 1, rU: 1, rW: 1}
expect_zero("Phi_T at equal free coordinates - 1", PhiT.subs(special_free_equal) - 1)
expect_zero("Phi_K at equal free coordinates - 1", PhiK.subs(special_free_equal) - 1)
expect_zero("Phi_mu at equal free coordinates - 1", PhiMu.subs(special_free_equal) - 1)
expect_zero("m_T at equal free coordinates - r_T", mT.subs(special_free_equal) - rT)
expect_zero("m_K at equal free coordinates - r_K", mK.subs(special_free_equal) - rK)
expect_zero("m_mu at equal free coordinates - r_mu", mMu.subs(special_free_equal) - rmu)

banner("STAGE 182 LEDGER")
print("1. For two positive microscopic states x^(1), x^(2), the exact finite pairwise log-ratio vector")
print("      Delta x^(2<-1) = ln(x^(2)/x^(1))")
print("   is acted on by the same Stage 192 map M_*.")
print("   Because the coherent monomials are multiplicative, this gives the exact two-point quotient packet")
print("      q^(2<-1) = M_* Delta x^(2<-1)")
print("   with")
print("      q_tr = ln(Ctr_2/Ctr_1),   q_nt = ln(Cnt_2/Cnt_1),   q_eta = ln(epsEta_2/epsEta_1).")
print("2. If the two states are on the same exact similarity orbit, then the dependent pairwise ratios are determined uniquely by the five free pairwise ratios through")
print("      Phi_T = r_U (r_U/(r_gamma r_c))^((1+deltaU_*)/(1+chi0_*)),")
print("      Phi_K = r_c^2 / r_U,")
print("      Phi_mu = Phi_K r_W^2 / r_lambda^2 * (r_gamma^2 r_lambda^2/(r_U r_W))^(-E_*) * (Phi_T/r_U)^(F_*).")
print("3. For an arbitrary pair, the exact reference-independent mismatch triple is")
print("      m_T = r_T/Phi_T,   m_K = r_K/Phi_K,   m_mu = r_mu/Phi_mu.")
print("   The pairwise invariant ratios collapse exactly to")
print("      Ctr_2/Ctr_1   = m_T^(1+chi0_*),")
print("      epsEta_2/epsEta_1 = 1/m_K,")
print("      Cnt_2/Cnt_1   = m_mu/(m_K m_T^(F_*)).")
print("4. Therefore the exact Packet-B logarithmic chart is again")
print("      q_tr = (1+chi0_*) ln m_T,   q_eta = - ln m_K,   q_nt = ln m_mu - ln m_K - F_* ln m_T.")
print("5. Applying the Stage 192 projectors to the finite pairwise log-ratio vector gives the exact two-point decomposition")
print("      Delta x^(2<-1) = O_orb Delta x^(2<-1) + Q_quot Delta x^(2<-1),")
print("   where O_orb carries the free ratios and the transported dependent ratios")
print("      (Phi_T, Phi_K, Phi_mu),")
print("   while Q_quot has support only on the dependent triple and carries")
print("      (ln m_T, ln m_K, ln m_mu).")
print("6. The pairwise transport, mismatch, and quotient packets compose exactly:")
print("      Phi^31 = Phi^32 Phi^21,   m^31 = m^32 m^21,   q^31 = q^32 + q^21.")
print("7. Hence the exact two-point orbit-lock criterion is")
print("      x^(2) in G_* . x^(1)")
print("   iff")
print("      m_T = m_K = m_mu = 1")
print("   iff")
print("      q_tr = q_nt = q_eta = 0.")
print("8. If the two states share the same five free microscopic coordinates, then")
print("      Phi_T = Phi_K = Phi_mu = 1,")
print("   so Stage 199 reduces exactly to the Stage 198 fixed-base formulas.")
