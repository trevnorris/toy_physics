#!/usr/bin/env python3
"""
moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py

SymPy-backed audit for Stage 161.

Checks:
1. Exact similarity-defect decomposition
     B_W = ((1+r_c)/9) (eps_gamma - eps_kappa).
2. Linearized slippage law
     dB_W = ((1+r_c*)/9) (d eps_gamma - d eps_kappa).
3. Exact D/N-tube expression for eps_kappa and its first variation.
4. Exact cancellation
     d eps_gamma - d eps_kappa = d ln gamma_0 - 2 d ln(LW/a).
5. Collapse of the final defect law to the single similarity-slippage scalar.
6. D/N similarity preservation theorem Xi_gamma = 2 Xi_L => Delta_Q = 0.
"""

from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 161 — D/N SIMILARITY SLIPPAGE DECOMPOSITION")

rc, eps_kappa, eps_gamma = sp.symbols("r_c eps_kappa eps_gamma", real=True)
kappa0 = (1 + rc) * (1 + eps_kappa) / 3
gamma0 = (1 + rc) * (1 + eps_gamma) / 9
BW = sp.simplify(gamma0 - kappa0 / 3)
print("B_W =", BW)
expect_zero("exact similarity-defect decomposition", BW - (1 + rc) * (eps_gamma - eps_kappa) / 9)

# Linearized law about eps_kappa = eps_gamma = 0 derived from the exact BW above.
eps = sp.symbols("eps", real=True)
rc_star, drc, deps_k, deps_g = sp.symbols("r_c_star dr_c d_eps_k d_eps_g", real=True)
BW_pert = BW.subs({eps_kappa: eps * deps_k, eps_gamma: eps * deps_g, rc: rc_star + eps * drc})
dBW = sp.simplify(sp.diff(BW_pert, eps).subs(eps, 0))
print("dB_W =", dBW)
expect_zero("linearized slippage law", dBW - (1 + rc_star) * (deps_g - deps_k) / 9)

banner("D/N-TUBE EVEN DEFECT AND THE EXACT HYBRIDIZATION CANCELLATION")

LW, a, dLW, da = sp.symbols("LW a dLW da", positive=True, real=True)
# Exact even similarity-defect variable from kappa0 = 4 LW^2 / (pi^2 a^2).
epsk_exact = sp.simplify(12 * LW**2 / (sp.pi**2 * a**2 * (1 + rc)) - 1)
print("eps_kappa =", epsk_exact)

# First variation on the compensated branch 12 LW^2 = pi^2 a^2 (1+r_c).
depsk_direct = sp.diff(epsk_exact, LW) * dLW + sp.diff(epsk_exact, a) * da + sp.diff(epsk_exact, rc) * drc
depsk_branch = sp.simplify(depsk_direct.subs(12 * LW**2, sp.pi**2 * a**2 * (1 + rc)))
print("d eps_kappa =", depsk_branch)
depsk_target = 2 * dLW / LW - 2 * da / a - drc / (1 + rc)
depsk_diff = sp.simplify(sp.expand((depsk_branch - depsk_target) * sp.pi**2 * a**2 * (1 + rc) * LW))
depsk_diff = sp.simplify(depsk_diff.subs(12 * LW**2 - sp.pi**2 * a**2 * rc - sp.pi**2 * a**2, 0))
expect_zero("d eps_kappa identity", depsk_diff)

# Odd similarity-defect variable from gamma0 = (1+r_c)(1+eps_gamma)/9.
# Solve for eps_gamma exactly, then linearize about the branch.
gamma0_sym, dgamma0 = sp.symbols("gamma0 d_gamma0", positive=True, real=True)
epsg_exact = sp.simplify(9 * gamma0_sym / (1 + rc) - 1)
print("eps_gamma =", epsg_exact)
depsg_direct = sp.diff(epsg_exact, gamma0_sym) * dgamma0 + sp.diff(epsg_exact, rc) * drc
# On the compensated branch gamma0_* = (1+r_c)/9 (notes section 1).
depsg_branch = sp.simplify(depsg_direct.subs(gamma0_sym, (1 + rc) / 9))
print("d eps_gamma =", depsg_branch)

dln_gamma0 = sp.symbols("dln_gamma0", real=True)
# On the branch, d gamma_0 = gamma_{0,*} d ln gamma_0 = (1+r_c) d ln gamma_0 / 9.
expect_zero(
    "d eps_gamma = d ln gamma0 - d ln(1+r_c)",
    depsg_branch.subs(dgamma0, (1 + rc) * dln_gamma0 / 9) - (dln_gamma0 - drc / (1 + rc)),
)
diff_identity = sp.simplify(
    sp.expand(
        ((depsg_branch - depsk_branch) - ((9 * dgamma0 / (1 + rc)) - 2 * (dLW / LW - da / a)))
        * sp.pi**2 * LW * a**2 * (1 + rc)
    )
)
diff_identity = sp.simplify(diff_identity.subs(12 * LW**2 - sp.pi**2 * a**2 * rc - sp.pi**2 * a**2, 0))
expect_zero("difference identity", diff_identity)

banner("TANGENTIAL SUSCEPTIBILITY AND FINAL DEFECT LAW")

Xi_gamma, Xi_L, sigma_star = sp.symbols("Xi_gamma Xi_L sigma_*", real=True)
dPi_tan, dSigma0, dS, dThat = sp.symbols("dPi_tan dSigma0 dS dThat", real=True)

UpsilonPi = sp.simplify((1 + rc_star) * (Xi_gamma - 2 * Xi_L) / 9)
print("Upsilon_Pi =", UpsilonPi)

DeltaQ = sp.simplify(-9 * sigma_star * UpsilonPi * dPi_tan / ((1 - sigma_star) * (1 + rc_star)))
NQm1 = sp.simplify(9 * sigma_star * UpsilonPi * dPi_tan / ((1 - sigma_star) * (1 + rc_star)))
print("Delta_Q =", DeltaQ)
print("N_Q - 1 =", NQm1)
expect_zero(
    "collapsed Delta_Q law",
    DeltaQ + sigma_star * (Xi_gamma - 2 * Xi_L) * dPi_tan / (1 - sigma_star),
)
expect_zero(
    "collapsed N_Q-1 law",
    NQm1 - sigma_star * (Xi_gamma - 2 * Xi_L) * dPi_tan / (1 - sigma_star),
)

# Insert the Stage 159 tangential mouth transport.
dPi_tan_expr = sp.Float("0.832409471081635") * dSigma0 - sp.Float("1.16275838754222") * dS
DeltaQ_mouth = sp.expand(DeltaQ.subs(dPi_tan, dPi_tan_expr))
print("Delta_Q in (dSigma0,dS) =", DeltaQ_mouth)

# Use dSigma0 ≈ 6.42981496203006 dThat.
DeltaQ_T = sp.expand(DeltaQ_mouth.subs(dSigma0, sp.Float("6.42981496203006") * dThat))
print("Delta_Q in (dThat,dS) =", DeltaQ_T)

banner("D/N SIMILARITY PRESERVATION")
expect_zero(
    "Xi_gamma = 2 Xi_L => Delta_Q = 0",
    DeltaQ.subs(Xi_gamma, 2 * Xi_L),
)
expect_zero(
    "Xi_gamma = 2 Xi_L => N_Q - 1 = 0",
    NQm1.subs(Xi_gamma, 2 * Xi_L),
)

rF1 = sp.Float("1.77799353547498")
Upsilon_prefactor = sp.N((1 + rF1**2) / 9, 18)
print("(1+r_F1^2)/9 =", Upsilon_prefactor)

banner("Carry-forward formulas")
print("1) B_W = ((1+r_c)/9) (eps_gamma - eps_kappa)")
print("2) dB_W = ((1+r_c*)/9) (d eps_gamma - d eps_kappa)")
print("3) d eps_gamma - d eps_kappa = d ln gamma_0 - 2 d ln(LW/a)")
print("4) Upsilon_Pi = ((1+r_c*)/9) (Xi_gamma - 2 Xi_L)")
print("5) Delta_Q = -(sigma_*/(1-sigma_*)) (Xi_gamma - 2 Xi_L) dPi_tan")
print("6) If Xi_gamma = 2 Xi_L, then the full first-order defect vanishes.")
