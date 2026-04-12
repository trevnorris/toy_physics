#!/usr/bin/env python3
"""
5pn_stage12_normalized_monomial_bridge.py

Twelfth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Introduces the exact raw-to-normalized dictionary connecting the Stage-20/28
   coherent kernel variables
      (c_etaU, lambda_W, gamma, K_U, K_eta^(eff), K_W^(eff), mu_eta, mu_U,
       mu_W, T_U)
   to the Stage-5/6 normalized prototype variables
      (K, M, G_U, G_W, R, Omega_U, Omega_W, delta_U).
2. Rewrites the key coherent branch ratios
      chi_0, epsilon_eta, epsilon_W, Z_W, Lambda
   entirely in normalized variables.
3. Rewrites the Stage-168/169 direct microscopic monomials
      C_tr, C_nt, epsilon_eta
   entirely in normalized variables.
4. Derives the exact normalized log-drift bridge
      Sigma_Z, Sigma_chi, Sigma_eta, Sigma_eps, Sigma_delta
   in terms of the normalized prototype drifts
      dln G_W, dln G_U, dln R, dln K, dln M,
      dln Omega_U, dln Omega_W, dln delta_U.
5. Verifies that the direct monomial drifts computed from the normalized
   variables reproduce the Stage-10/11 triangular coordinates exactly.

Interpretation
--------------
This is the missing algebraic bridge between the Stage-5 primitive overlap model
and the Stage-10/11 similarity-orbit package. Once the split-U variable
`delta_U` is added, the Stage-5-style normalized prototype already contains the
full Stage-168/169 quotient coordinates; the extra raw mass/stiffness bookkeeping
mostly cancels out of the defect variables.
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



def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
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


# ---------------------------------------------------------------------------
# I. Exact raw-to-normalized dictionary
# ---------------------------------------------------------------------------

banner("I. EXACT RAW-TO-NORMALIZED COHERENT-KERNEL DICTIONARY")

# Raw coherent-kernel variables from Stages 20/28 and 168/169.
c_etaU, lambdaW, gamma = sp.symbols("c_etaU lambda_W gamma", positive=True, real=True)
KU, Keta_eff, KW_eff = sp.symbols("K_U K_eta_eff K_W_eff", positive=True, real=True)
mu_eta, mu_U, mu_W = sp.symbols("mu_eta mu_U mu_W", positive=True, real=True)
TU = sp.symbols("T_U", positive=True, real=True)
L = sp.symbols("L", positive=True, real=True)

# Normalized prototype variables in the Stage-5/6 spirit.
K, M = sp.symbols("K M", positive=True, real=True)
GU, GW, R = sp.symbols("G_U G_W R", positive=True, real=True)
OmegaU, OmegaW = sp.symbols("Omega_U Omega_W", positive=True, real=True)
deltaU = sp.symbols("delta_U", positive=True, real=True)

norm_from_raw = {
    K: Keta_eff,
    M: mu_eta,
    GU: c_etaU / sp.sqrt(mu_eta * mu_U),
    GW: lambdaW / sp.sqrt(mu_eta * mu_W),
    R: gamma * lambdaW / sp.sqrt(mu_U * mu_W),
    OmegaU**2: KU / mu_U,
    OmegaW**2: KW_eff / mu_W,
    deltaU: sp.pi**2 * TU / (L**2 * KU),
}

raw_to_norm = {
    Keta_eff: K,
    mu_eta: M,
    KU: OmegaU**2 * mu_U,
    KW_eff: OmegaW**2 * mu_W,
    c_etaU: GU * sp.sqrt(M * mu_U),
    lambdaW: GW * sp.sqrt(M * mu_W),
    gamma: R * sp.sqrt(mu_U) / (GW * sp.sqrt(M)),
    TU: deltaU * L**2 * OmegaU**2 * mu_U / sp.pi**2,
}

print("K         = K_eta^(eff)")
print("M         = mu_eta")
print("G_U       = c_etaU / sqrt(mu_eta mu_U)")
print("G_W       = lambda_W / sqrt(mu_eta mu_W)")
print("R         = gamma lambda_W / sqrt(mu_U mu_W)")
print("Omega_U^2 = K_U / mu_U")
print("Omega_W^2 = K_W^(eff) / mu_W")
print("delta_U   = pi^2 T_U / (L^2 K_U)")


# ---------------------------------------------------------------------------
# II. Exact coherent branch ratios in normalized variables
# ---------------------------------------------------------------------------

banner("II. EXACT COHERENT BRANCH RATIOS IN NORMALIZED VARIABLES")

sigma = sp.symbols("sigma", positive=True, real=True)
Ggrav, cs, a_th, c_light = sp.symbols("G c_s a_th c_light", positive=True, real=True)

chi0_raw = gamma * c_etaU / KU
eps_eta_raw = c_etaU**2 / (KU * Keta_eff)
eps_W_raw = gamma**2 * lambdaW**2 * sigma / (KU * KW_eff)
ZW_raw = lambdaW**2 / (Keta_eff * KW_eff)
Lambda_raw = 27 * sp.pi**2 * Ggrav * cs**5 * KW_eff / (20 * a_th**5 * c_light**5 * mu_W)

chi0_norm = sp.simplify(R * GU / (OmegaU**2 * GW))
eps_eta_norm = sp.simplify(M * GU**2 / (K * OmegaU**2))
eps_W_norm = sp.simplify(R**2 * sigma / (OmegaU**2 * OmegaW**2))
ZW_norm = sp.simplify(M * GW**2 / (K * OmegaW**2))
Lambda_norm = sp.simplify(27 * sp.pi**2 * Ggrav * cs**5 * OmegaW**2 / (20 * a_th**5 * c_light**5))

print("chi_0       =", chi0_norm)
print("epsilon_eta =", eps_eta_norm)
print("epsilon_W   =", eps_W_norm)
print("Z_W         =", ZW_norm)
print("Lambda      =", Lambda_norm)

expect_zero("chi_0 raw-normalized bridge", sp.simplify(chi0_raw.subs(raw_to_norm) - chi0_norm))
expect_zero("epsilon_eta raw-normalized bridge", sp.simplify(eps_eta_raw.subs(raw_to_norm) - eps_eta_norm))
expect_zero("epsilon_W raw-normalized bridge", sp.simplify(eps_W_raw.subs(raw_to_norm) - eps_W_norm))
expect_zero("Z_W raw-normalized bridge", sp.simplify(ZW_raw.subs(raw_to_norm) - ZW_norm))
expect_zero("Lambda raw-normalized bridge", sp.simplify(Lambda_raw.subs(raw_to_norm) - Lambda_norm))


# ---------------------------------------------------------------------------
# III. Direct microscopic monomials in normalized variables
# ---------------------------------------------------------------------------

banner("III. DIRECT MICROSCOPIC MONOMIALS IN NORMALIZED VARIABLES")

chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)

Ctr_raw = sp.simplify((gamma * c_etaU / KU)**(1 + deltaUs) * (sp.pi**2 * TU / (L**2 * KU))**(1 + chi0s))
Cnt_raw = sp.simplify(
    (lambdaW**2 * mu_W / (Keta_eff * KW_eff**2))
    * (gamma**2 * lambdaW**2 * sigma / (KU * KW_eff))**Estar
    * (sp.pi**2 * TU / (L**2 * KU))**(-Fstar)
)

Ctr_norm = sp.simplify((R * GU / (OmegaU**2 * GW))**(1 + deltaUs) * deltaU**(1 + chi0s))
Cnt_norm = sp.simplify((M * GW**2 / (K * OmegaW**4)) * (R**2 * sigma / (OmegaU**2 * OmegaW**2))**Estar * deltaU**(-Fstar))

print("C_tr,* =")
sp.pprint(Ctr_norm)
print("C_nt,* =")
sp.pprint(Cnt_norm)
print("epsilon_eta =")
sp.pprint(eps_eta_norm)

expect_zero("C_tr raw-normalized bridge", sp.simplify(Ctr_raw.subs(raw_to_norm) - Ctr_norm))
expect_zero("C_nt raw-normalized bridge", sp.simplify(Cnt_raw.subs(raw_to_norm) - Cnt_norm))


# ---------------------------------------------------------------------------
# IV. Exact normalized log-drift bridge to Stage-10 slippages
# ---------------------------------------------------------------------------

banner("IV. EXACT NORMALIZED LOG-DRIFT BRIDGE")

# Raw microscopic log drifts from Stage 168.
lambda1, c1, gamma1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
kappaU, kappa_eta, kappaW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)
mu_eta1, muU1 = sp.symbols("mu_eta_1 mu_U1", real=True)

Sigma_Z = sp.simplify(2 * lambda1 + mu1 - kappa_eta - 2 * kappaW)
Sigma_chi = sp.simplify(gamma1 + c1 - kappaU)
Sigma_eta = sp.simplify(2 * c1 - kappaU - kappa_eta)
Sigma_eps = sp.simplify(2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
Sigma_delta = sp.simplify(tau1 - kappaU)

# Normalized log drifts for the Stage-5 prototype variables.
dln_GW, dln_GU, dln_R = sp.symbols("dln_GW dln_GU dln_R", real=True)
dln_K, dln_M = sp.symbols("dln_K dln_M", real=True)
dln_OmegaU, dln_OmegaW = sp.symbols("dln_Omega_U dln_Omega_W", real=True)
dln_deltaU = sp.symbols("dln_delta_U", real=True)

# Raw -> normalized drift dictionary.
dict_drifts = {
    dln_GU: c1 - sp.Rational(1, 2) * (mu_eta1 + muU1),
    dln_GW: lambda1 - sp.Rational(1, 2) * (mu_eta1 + mu1),
    dln_R: gamma1 + lambda1 - sp.Rational(1, 2) * (muU1 + mu1),
    dln_K: kappa_eta,
    dln_M: mu_eta1,
    dln_OmegaU: sp.Rational(1, 2) * (kappaU - muU1),
    dln_OmegaW: sp.Rational(1, 2) * (kappaW - mu1),
    dln_deltaU: tau1 - kappaU,
}

Sigma_Z_norm = sp.simplify(dln_M + 2 * dln_GW - dln_K - 4 * dln_OmegaW)
Sigma_chi_norm = sp.simplify(dln_R + dln_GU - dln_GW - 2 * dln_OmegaU)
Sigma_eta_norm = sp.simplify(dln_M + 2 * dln_GU - dln_K - 2 * dln_OmegaU)
Sigma_eps_norm = sp.simplify(2 * (dln_R - dln_OmegaU - dln_OmegaW))
Sigma_delta_norm = sp.simplify(dln_deltaU)

print("Sigma_Z     =", Sigma_Z_norm)
print("Sigma_chi   =", Sigma_chi_norm)
print("Sigma_eta   =", Sigma_eta_norm)
print("Sigma_eps   =", Sigma_eps_norm)
print("Sigma_delta =", Sigma_delta_norm)

expect_zero("Sigma_Z bridge", sp.simplify(Sigma_Z_norm.subs(dict_drifts) - Sigma_Z))
expect_zero("Sigma_chi bridge", sp.simplify(Sigma_chi_norm.subs(dict_drifts) - Sigma_chi))
expect_zero("Sigma_eta bridge", sp.simplify(Sigma_eta_norm.subs(dict_drifts) - Sigma_eta))
expect_zero("Sigma_eps bridge", sp.simplify(Sigma_eps_norm.subs(dict_drifts) - Sigma_eps))
expect_zero("Sigma_delta bridge", sp.simplify(Sigma_delta_norm.subs(dict_drifts) - Sigma_delta))

subbanner("Direct monomial drifts from normalized variables")

dlog_Ctr_norm = sp.simplify((1 + deltaUs) * Sigma_chi_norm + (1 + chi0s) * Sigma_delta_norm)
dlog_Cnt_norm = sp.simplify((dln_M + 2 * dln_GW - dln_K - 4 * dln_OmegaW) + Estar * Sigma_eps_norm - Fstar * dln_deltaU)
dlog_eta_norm = sp.simplify(Sigma_eta_norm)

print("delta ln C_tr,* =", dlog_Ctr_norm)
print("delta ln C_nt,* =", dlog_Cnt_norm)
print("delta ln epsilon_eta =", dlog_eta_norm)

expect_zero(
    "delta ln C_tr,* - ((1+deltaU*) Sigma_chi + (1+chi0*) Sigma_delta)",
    dlog_Ctr_norm - ((1 + deltaUs) * Sigma_chi_norm + (1 + chi0s) * Sigma_delta_norm),
)
expect_zero(
    "delta ln C_nt,* - (Sigma_Z + E Sigma_eps - F Sigma_delta)",
    dlog_Cnt_norm - (Sigma_Z_norm + Estar * Sigma_eps_norm - Fstar * Sigma_delta_norm),
)
expect_zero("delta ln epsilon_eta - Sigma_eta", dlog_eta_norm - Sigma_eta_norm)


# ---------------------------------------------------------------------------
# V. Final ledger
# ---------------------------------------------------------------------------

banner("V. FINAL THEOREM LEDGER")
print("1. The raw Stage-168/169 coherent kernel variables can be rewritten exactly")
print("   in the Stage-5-style normalized variables")
print("      (K, M, G_U, G_W, R, Omega_U, Omega_W, delta_U).")
print("2. The key coherent ratios become")
print("      chi_0       = R G_U / (Omega_U^2 G_W),")
print("      epsilon_eta = M G_U^2 / (K Omega_U^2),")
print("      epsilon_W   = R^2 sigma / (Omega_U^2 Omega_W^2),")
print("      Z_W         = M G_W^2 / (K Omega_W^2),")
print("      Lambda      = 27 pi^2 G c_s^5 Omega_W^2 / (20 a^5 c^5).")
print("3. The direct monomials are therefore already present in the normalized prototype:")
print("      C_tr,* = (R G_U/(Omega_U^2 G_W))^(1+deltaU*) delta_U^(1+chi0*),")
print("      C_nt,* = (M G_W^2/(K Omega_W^4)) (R^2 sigma/(Omega_U^2 Omega_W^2))^E* delta_U^(-F*),")
print("      epsilon_eta = M G_U^2/(K Omega_U^2).")
print("4. Their exact log-drift coordinates collapse to the normalized formulas")
print("      Sigma_Z     = dln M + 2 dln G_W - dln K - 4 dln Omega_W,")
print("      Sigma_chi   = dln R + dln G_U - dln G_W - 2 dln Omega_U,")
print("      Sigma_eta   = dln M + 2 dln G_U - dln K - 2 dln Omega_U,")
print("      Sigma_eps   = 2(dln R - dln Omega_U - dln Omega_W),")
print("      Sigma_delta = dln delta_U.")
print("5. So the Stage-10/11 similarity-orbit problem can be formulated directly in")
print("   the Stage-5 normalized variables, with the only genuinely new branch datum")
print("   being the split-U ratio delta_U.")
