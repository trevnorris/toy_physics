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
        expr = expr.applyfunc(lambda z: sp.simplify(z))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(expr)
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 198 — EXACT FINITE ORBIT LAW FOR THE DEPENDENT TRIPLE")

# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------
lamW, cEtaU, gamma = sp.symbols("lambda_W c_etaU gamma", positive=True, real=True)
KU, Keta, KW, muW, TU = sp.symbols(
    "K_U K_eta_eff K_W_eff mu_W T_U", positive=True, real=True
)
L, sigma = sp.symbols("L sigma", positive=True, real=True)

chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)
Ctr_star, Cnt_star, epsEta_star = sp.symbols(
    "Ctr_star Cnt_star epsEta_star", positive=True, real=True
)

pi = sp.pi

# Exact coherent monomials
Ctr = sp.simplify(
    (gamma * cEtaU / KU) ** (1 + deltaUs) * (pi**2 * TU / (L**2 * KU)) ** (1 + chi0s)
)
Cnt = sp.simplify(
    (lamW**2 * muW / (Keta * KW**2))
    * (gamma**2 * lamW**2 * sigma / (KU * KW)) ** Estar
    * (pi**2 * TU / (L**2 * KU)) ** (-Fstar)
)
epsEta = sp.simplify(cEtaU**2 / (KU * Keta))

subbanner("I. Exact coherent monomials")
print("Ctr =")
sp.pprint(Ctr)
print("Cnt =")
sp.pprint(Cnt)
print("epsEta =")
sp.pprint(epsEta)

# ---------------------------------------------------------------------------
# Exact orbit solve for dependent triple
# ---------------------------------------------------------------------------
subbanner("II. Exact finite orbit solve for (T_U, K_eta, mu_W)")

Keta_orbit = sp.simplify(cEtaU**2 / (KU * epsEta_star))
TU_orbit = sp.simplify(
    (L**2 * KU / pi**2)
    * (Ctr_star / (gamma * cEtaU / KU) ** (1 + deltaUs)) ** (sp.Rational(1, 1) / (1 + chi0s))
)
muW_orbit = sp.simplify(
    Cnt_star * Keta_orbit * KW**2 / lamW**2
    * (gamma**2 * lamW**2 * sigma / (KU * KW)) ** (-Estar)
    * (pi**2 * TU_orbit / (L**2 * KU)) ** Fstar
)

print("K_eta^(orbit) =")
sp.pprint(Keta_orbit)
print("T_U^(orbit) =")
sp.pprint(TU_orbit)
print("mu_W^(orbit) =")
sp.pprint(muW_orbit)

expect_zero(
    "epsEta(Keta_orbit) / epsEta_star - 1",
    epsEta.subs({Keta: Keta_orbit}) / epsEta_star - 1,
)
expect_zero(
    "Ctr(TU_orbit) / Ctr_star - 1",
    Ctr.subs({TU: TU_orbit}) / Ctr_star - 1,
)
expect_zero(
    "Cnt(muW_orbit,Keta_orbit,TU_orbit) / Cnt_star - 1",
    Cnt.subs({muW: muW_orbit, Keta: Keta_orbit, TU: TU_orbit}) / Cnt_star - 1,
)

# ---------------------------------------------------------------------------
# Exact mismatch ratios and invariant ratios
# ---------------------------------------------------------------------------
subbanner("III. Exact dependent residual mismatch ratios")

mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)

TU_actual = sp.simplify(mT * TU_orbit)
Keta_actual = sp.simplify(mK * Keta_orbit)
muW_actual = sp.simplify(mMu * muW_orbit)

Ctr_ratio = sp.simplify(Ctr.subs({TU: TU_actual}) / Ctr_star)
epsEta_ratio = sp.simplify(epsEta.subs({Keta: Keta_actual}) / epsEta_star)
Cnt_ratio = sp.simplify(
    Cnt.subs({TU: TU_actual, Keta: Keta_actual, muW: muW_actual}) / Cnt_star
)

print("Ctr(actual) / Ctr(orbit) =")
sp.pprint(Ctr_ratio)
print("epsEta(actual) / epsEta_star =")
sp.pprint(epsEta_ratio)
print("Cnt(actual) / Cnt(orbit) =")
sp.pprint(Cnt_ratio)

expect_zero("Ctr ratio - m_T^(1+chi0_*)", Ctr_ratio - mT ** (1 + chi0s))
expect_zero("epsEta ratio - 1/m_K", epsEta_ratio - 1 / mK)
expect_zero("Cnt ratio - m_mu/(m_K m_T^F_*)", Cnt_ratio - mMu / (mK * mT**Fstar))

# ---------------------------------------------------------------------------
# Exact logarithmic chart and Stage 192 consistency
# ---------------------------------------------------------------------------
subbanner("IV. Exact logarithmic chart and agreement with the Stage 192 drift map")

tau, kappa, mu = sp.symbols("tau kappa mu", real=True)
qtr_expected = sp.simplify((1 + chi0s) * tau)
qeta_expected = sp.simplify(-kappa)
qnt_expected = sp.simplify(mu - kappa - Fstar * tau)

print("q_tr =")
sp.pprint(qtr_expected)
print("q_nt =")
sp.pprint(qnt_expected)
print("q_eta =")
sp.pprint(qeta_expected)

# Ordered drift basis: (Delta_lambda, Delta_c, Delta_gamma, Delta_U, Delta_Keta, Delta_W, Delta_mu, Delta_T)
Mstar = sp.Matrix(
    [
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)
Dx_mis = sp.Matrix([0, 0, 0, 0, kappa, 0, mu, tau])
q_from_M = sp.simplify(Mstar * Dx_mis)
print("M_* Delta x_mis =")
sp.pprint(q_from_M)
expect_zero(
    "M_* Delta x_mis - (q_tr,q_nt,q_eta)",
    q_from_M - sp.Matrix([qtr_expected, qnt_expected, qeta_expected]),
)

# ---------------------------------------------------------------------------
# Exact restoration map
# ---------------------------------------------------------------------------
subbanner("V. Exact restoration map")

TU_actual_log = sp.simplify(sp.exp(tau) * TU_orbit)
Keta_actual_log = sp.simplify(sp.exp(kappa) * Keta_orbit)
muW_actual_log = sp.simplify(sp.exp(mu) * muW_orbit)

T_restore = sp.simplify(TU_actual_log * sp.exp(-qtr_expected / (1 + chi0s)))
K_restore = sp.simplify(Keta_actual_log * sp.exp(qeta_expected))
mu_restore = sp.simplify(
    muW_actual_log
    * sp.exp(-qnt_expected + qeta_expected - Fstar * qtr_expected / (1 + chi0s))
)

print("T_U^(restore) =")
sp.pprint(T_restore)
print("K_eta^(restore) =")
sp.pprint(K_restore)
print("mu_W^(restore) =")
sp.pprint(mu_restore)
expect_zero("T_restore - T_orbit", T_restore - TU_orbit)
expect_zero("K_restore - K_orbit", K_restore - Keta_orbit)
expect_zero("mu_restore - mu_orbit", mu_restore - muW_orbit)

# ---------------------------------------------------------------------------
# Exact inverse chart and orbit-lock criterion
# ---------------------------------------------------------------------------
subbanner("VI. Exact inverse chart and finite orbit-lock criterion")

qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)
tau_inv = sp.simplify(qtr / (1 + chi0s))
kappa_inv = sp.simplify(-qeta)
mu_inv = sp.simplify(qnt - qeta + Fstar * qtr / (1 + chi0s))

print("tau(q) =")
sp.pprint(tau_inv)
print("kappa(q) =")
sp.pprint(kappa_inv)
print("mu(q) =")
sp.pprint(mu_inv)
expect_zero(
    "direct/inverse chart consistency for q_tr",
    sp.simplify((1 + chi0s) * tau_inv - qtr),
)
expect_zero(
    "direct/inverse chart consistency for q_eta",
    sp.simplify(-kappa_inv - qeta),
)
expect_zero(
    "direct/inverse chart consistency for q_nt",
    sp.simplify(mu_inv - kappa_inv - Fstar * tau_inv - qnt),
)

mT_from_q = sp.simplify(sp.exp(tau_inv))
mK_from_q = sp.simplify(sp.exp(kappa_inv))
mMu_from_q = sp.simplify(sp.exp(mu_inv))
print("m_T(q) =")
sp.pprint(mT_from_q)
print("m_K(q) =")
sp.pprint(mK_from_q)
print("m_mu(q) =")
sp.pprint(mMu_from_q)
expect_zero("m_T at q=0 - 1", mT_from_q.subs({qtr: 0}) - 1)
expect_zero("m_K at q=0 - 1", mK_from_q.subs({qeta: 0}) - 1)
expect_zero("m_mu at q=0 - 1", mMu_from_q.subs({qtr: 0, qnt: 0, qeta: 0}) - 1)

banner("STAGE 198 LEDGER")
print("1. The exact coherent monomials are triangular in the dependent microscopic triple")
print("      (T_U, K_eta^(eff), mu_W),")
print("   so for fixed free microscopic coordinates and fixed invariant triple")
print("      (Ctr_*, Cnt_*, epsEta_*),")
print("   the dependent triple is determined uniquely by")
print("      K_eta^(orbit) = c_etaU^2 / (K_U epsEta_*),")
print("      T_U^(orbit)   = (L^2 K_U / pi^2) [Ctr_* / (gamma c_etaU / K_U)^(1+deltaU_*)]^(1/(1+chi0_*)),")
print("      mu_W^(orbit)  = Cnt_* K_eta^(orbit) K_W^2 / lambda_W^2 * (gamma^2 lambda_W^2 sigma / (K_U K_W))^(-E_*) * (pi^2 T_U^(orbit)/(L^2 K_U))^(F_*).")
print("2. A candidate branch with the same five free microscopic coordinates is therefore compared to the orbit point by the exact dependent residual mismatch ratios")
print("      T_U = m_T T_U^(orbit),   K_eta^(eff) = m_K K_eta^(orbit),   mu_W = m_mu mu_W^(orbit).")
print("3. The invariant ratios are exactly")
print("      Ctr(actual)/Ctr(orbit) = m_T^(1+chi0_*),")
print("      epsEta(actual)/epsEta_* = 1/m_K,")
print("      Cnt(actual)/Cnt(orbit) = m_mu / (m_K m_T^(F_*)).")
print("4. Writing tau = ln m_T, kappa = ln m_K, mu = ln m_mu, the Packet-B quotient coordinates are exactly")
print("      q_tr = (1+chi0_*) tau,   q_eta = -kappa,   q_nt = mu - kappa - F_* tau.")
print("   Applying the Stage 192 drift compiler M_* to a pure dependent mismatch vector reproduces this chart exactly.")
print("5. Restoration to the same orbit changes only the dependent triple:")
print("      T_U^(restore)   = T_U exp[-q_tr/(1+chi0_*)],")
print("      K_eta^(restore) = K_eta exp[q_eta],")
print("      mu_W^(restore)  = mu_W exp[-q_nt + q_eta - F_* q_tr/(1+chi0_*)].")
print("6. Therefore the finite orbit-lock condition is exactly")
print("      m_T = m_K = m_mu = 1  <=>  q_tr = q_nt = q_eta = 0.")
print("   This is the exact Packet-B complement of the Stage 197 Packet-A finish-line theorem.")
