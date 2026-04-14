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
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("SAME-CHARGE BARRIER AUDIT — STAGE 025")
banner("ACTUAL TWIN-SUPPORT PLACEMENT AND DIRECT-OBSERVABLE ORBIT-LOCK COMPILER")

chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta = sp.symbols(
    "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda zeta", positive=True, real=True
)
beta = sp.symbols("beta", positive=True, real=True)
pi = sp.pi

subbanner("I. Actual coherent placement state -> selected twin-support position")

epsilon = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
Cmix = sp.simplify(8 * Lambda * (1 - epsilon) / pi**2)
Pi_tr = sp.simplify(sp.Rational(4, 3) * Cmix)
varrho = sp.simplify(pi**2 * Pi_tr / (16 * Lambda))
sigma = sp.simplify(4 / (3 * varrho) - 2)

print(f"epsilon = {epsilon}")
print(f"C_mix = {Cmix}")
print(f"Pi_tr = {Pi_tr}")
print(f"varrho = {varrho}")
print(f"sigma = {sigma}")

expect_zero("varrho - 2(1-epsilon)/3", varrho - sp.Rational(2, 3) * (1 - epsilon))
expect_zero("sigma - 2 epsilon/(1-epsilon)", sigma - 2 * epsilon / (1 - epsilon))
expect_zero("Pi_tr/C_mix - 4/3", sp.simplify(Pi_tr / Cmix) - sp.Rational(4, 3))

subbanner("II. Support-regime meaning on the selected branch")
print(f"Pi_tr - C_mix = {sp.simplify(Pi_tr - Cmix)}")
print(f"2 C_mix - Pi_tr = {sp.simplify(2 * Cmix - Pi_tr)}")
expect_zero("Pi_tr - C_mix - C_mix/3", (Pi_tr - Cmix) - Cmix / 3)
expect_zero("2C_mix - Pi_tr - 2C_mix/3", (2 * Cmix - Pi_tr) - 2 * Cmix / 3)
print("Selected branch sits strictly inside the lowest-symmetric-twin window:")
print("    C_mix < Pi_tr < 2 C_mix")

subbanner("III. Exact ranking thresholds rewritten in the physical epsilon variable")
varrho_WL = sp.simplify(2 * (1 + beta**2) / (3 * (2 + beta**2)))
varrho_UL = sp.simplify(2 * (1 + beta**2) / (3 * (1 + beta + beta**2)))
epsilon_WL = sp.simplify(1 - sp.Rational(3, 2) * varrho_WL)
epsilon_UL = sp.simplify(1 - sp.Rational(3, 2) * varrho_UL)

print(f"epsilon_(WΛ) = {epsilon_WL}")
print(f"epsilon_(UΛ) = {epsilon_UL}")
expect_zero("epsilon_(WΛ) - 1/(2+beta^2)", epsilon_WL - 1 / (2 + beta**2))
expect_zero(
    "epsilon_(UΛ) - beta/(1+beta+beta^2)",
    epsilon_UL - beta / (1 + beta + beta**2),
)

subbanner("IV. Direct orbit observables on the actual coherent branch")
R_tr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
R_target = sp.simplify(Lambda * (1 - epsEta) * (1 - epsilon) ** 2 / (ZW * (1 + chi0) ** 2))
M_mix = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (pi**2 * (1 - epsEta) * (1 - epsilon)))
S = sp.simplify(1 + zeta * (1 - epsilon) / (1 - zeta * epsilon))
M_tr = sp.simplify(M_mix * S)

print(f"R_tr = {R_tr}")
print(f"R_target = {R_target}")
print(f"M_mix = {M_mix}")
print(f"S(zeta; epsilon) = {S}")
print(f"M_tr = {M_tr}")
expect_zero(
    "R_target*M_mix - 8 Lambda (1-epsilon)/pi^2",
    sp.simplify(R_target * M_mix) - 8 * Lambda * (1 - epsilon) / pi**2,
)

subbanner("V. Support-blindness of the orbit packet")
expect_zero("dR_tr/dzeta", sp.diff(R_tr, zeta))
expect_zero("dR_target/dzeta", sp.diff(R_target, zeta))
expect_zero("depsilon/dzeta", sp.diff(epsilon, zeta))
expect_zero("d epsEta / dzeta", sp.diff(epsEta, zeta))
print("zeta enters only the support packet through S(zeta; epsilon) and M_tr.")

subbanner("VI. Exact logarithmic drift compilers")
dchi0, ddeltaU, dZW, depsW, depsEta, dLambda = sp.symbols(
    "dln_chi0 dln_deltaU dln_ZW dln_epsilonW dln_epsilonEta dln_Lambda", real=True
)

dln_epsilon = sp.simplify(
    epsW * sp.diff(sp.log(epsilon), epsW) * depsW
    + deltaU * sp.diff(sp.log(epsilon), deltaU) * ddeltaU
)
dln_Rtr = sp.simplify(
    chi0 * sp.diff(sp.log(R_tr), chi0) * dchi0
    + deltaU * sp.diff(sp.log(R_tr), deltaU) * ddeltaU
)
dln_Rtarget = sp.simplify(
    Lambda * sp.diff(sp.log(R_target), Lambda) * dLambda
    + ZW * sp.diff(sp.log(R_target), ZW) * dZW
    + chi0 * sp.diff(sp.log(R_target), chi0) * dchi0
    + epsW * sp.diff(sp.log(R_target), epsW) * depsW
    + epsEta * sp.diff(sp.log(R_target), epsEta) * depsEta
    + deltaU * sp.diff(sp.log(R_target), deltaU) * ddeltaU
)

C_tr = sp.simplify(chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)))
Sigma_tr = sp.simplify((1 + deltaU) * dchi0 + (1 + chi0) * ddeltaU)

print(f"d ln epsilon = {dln_epsilon}")
print(f"d ln R_tr = {dln_Rtr}")
print(f"d ln R_target = {dln_Rtarget}")

expect_zero(
    "dln epsilon - [dln epsW - 2 deltaU dln deltaU /((1+deltaU)(11+9 deltaU))]",
    dln_epsilon - (depsW - 2 * deltaU * ddeltaU / ((1 + deltaU) * (11 + 9 * deltaU))),
)
expect_zero("dln R_tr + C_tr Sigma_tr", dln_Rtr + C_tr * Sigma_tr)
expect_zero(
    "dln R_target - [dlnLambda - dlnZW - 2 chi0 dlnchi0/(1+chi0) - epsEta dln epsEta/(1-epsEta) - 2 epsilon dln epsilon/(1-epsilon)]",
    dln_Rtarget - (
        dLambda
        - dZW
        - 2 * chi0 * dchi0 / (1 + chi0)
        - epsEta * depsEta / (1 - epsEta)
        - 2 * epsilon * dln_epsilon / (1 - epsilon)
    ),
)

subbanner("VII. Direct-observable defect packet")
Theta1 = sp.simplify(dln_Rtr)
R1 = sp.simplify(dln_Rtarget)
Xi1 = sp.simplify(-dln_Rtarget - epsEta * depsEta / (1 - epsEta))

print(f"Theta_1 = {Theta1}")
print(f"Xi_1 = {Xi1}")
print(f"R_1 = {R1}")
expect_zero(
    "Xi_1 - [-dlnLambda + dlnZW + 2 chi0 dlnchi0/(1+chi0) + 2 epsilon dln epsilon/(1-epsilon)]",
    Xi1 - (-dLambda + dZW + 2 * chi0 * dchi0 / (1 + chi0) + 2 * epsilon * dln_epsilon / (1 - epsilon)),
)
print("Zero-defect orbit lock is therefore exactly:")
print("    d ln R_tr = 0,   d ln R_target = 0,   d ln epsilon_eta = 0")
print("with the outgoing normalization finish line carried separately by N_Q = 1 (or chi_Q = 1 on the natural source-map branch).")

print("\nAll Stage 025 checks passed.")
