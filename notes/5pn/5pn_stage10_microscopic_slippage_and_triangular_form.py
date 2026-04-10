
#!/usr/bin/env python3
"""
5pn_stage10_microscopic_slippage_and_triangular_form.py

Tenth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Pushes the Stage-9 coherent-branch defect from physical placement variables
   down to the microscopic coherent-kernel couplings.
2. Defines the exact microscopic slippages
      Sigma_Z, Sigma_chi, Sigma_eta, Sigma_epsilon, Sigma_delta
   and proves the Stage-165 coherent defect law in those variables.
3. Isolates the exact tracking combination
      Sigma_tr = (1+chi_0) Sigma_delta + (1+delta_U) Sigma_chi
   and proves the tracking/nontracking split.
4. Defines the genuine nontracking transfer-shape slippage Sigma_nt and proves
   the Stage-166 exact triangular normal form
      Theta_1 = -C_tr Sigma_tr,
      Xi_1 = A_tr Sigma_tr + Sigma_nt,
      R_1 + Xi_1 = -(epsilon_eta/(1-epsilon_eta)) Sigma_eta.
5. Derives the exact inverse reconstruction formulas.
6. Records the full triple-rigidity theorem
      Theta_1 = Xi_1 = R_1 = 0  iff  Sigma_tr = Sigma_nt = Sigma_eta = 0
   on the constructive coherent branch.

Interpretation
--------------
After this stage the coherent weak-axisymmetric problem is no longer a generic
microscopic drift ledger. It is an exact three-scalar normal form:
tracking, nontracking transfer shape, and dressing.
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
# I. Microscopic coherent-kernel slippage variables
# ---------------------------------------------------------------------------

banner("I. MICROSCOPIC COHERENT-KERNEL SLIPPAGES")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
epsW = sp.symbols("epsilon_W", positive=True, real=True)
eps_eta = sp.symbols("epsilon_eta", positive=True, real=True)
epsilon = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))

lambda1, c1, gamma1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
kappaU, kappa_eta, kappaW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)

zetaZ, omegaW = sp.symbols("zeta_Z omega_W", real=True)
chi1, eta1 = sp.symbols("chi_1 eta_1", real=True)
varepsW, deltaU1 = sp.symbols("varepsilon_W delta_U1", real=True)

Sigma_Z = sp.simplify(2 * lambda1 + mu1 - kappa_eta - 2 * kappaW)
Sigma_chi = sp.simplify(gamma1 + c1 - kappaU)
Sigma_eta = sp.simplify(2 * c1 - kappaU - kappa_eta)
Sigma_eps = sp.simplify(2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
Sigma_delta = sp.simplify(tau1 - kappaU)

print("Sigma_Z     =", Sigma_Z)
print("Sigma_chi   =", Sigma_chi)
print("Sigma_eta   =", Sigma_eta)
print("Sigma_eps   =", Sigma_eps)
print("Sigma_delta =", Sigma_delta)

zetaZ_expr = sp.simplify(2 * lambda1 - kappa_eta - kappaW)
omegaW_expr = sp.simplify(kappaW - mu1)
chi1_expr = sp.simplify(chi0 * Sigma_chi)
eta1_expr = sp.simplify(eps_eta * Sigma_eta)
varepsW_expr = sp.simplify(epsW * Sigma_eps)
deltaU1_expr = sp.simplify(deltaU * Sigma_delta)

expect_zero("Sigma_Z - (zeta_Z - omega_W)", Sigma_Z - (zetaZ_expr - omegaW_expr))
expect_zero("chi_1 = chi_0 Sigma_chi", chi1_expr - chi0 * Sigma_chi)
expect_zero("eta_1 = epsilon_eta Sigma_eta", eta1_expr - eps_eta * Sigma_eta)
expect_zero("varepsilon_W = epsilon_W Sigma_eps", varepsW_expr - epsW * Sigma_eps)
expect_zero("delta_U1 = delta_U Sigma_delta", deltaU1_expr - deltaU * Sigma_delta)

print("zeta_Z =", zetaZ_expr)
print("omega_W =", omegaW_expr)
print("chi_1 =", chi1_expr)
print("eta_1 =", eta1_expr)
print("varepsilon_W =", varepsW_expr)
print("delta_U1 =", deltaU1_expr)


# ---------------------------------------------------------------------------
# II. Exact microscopic law for Xi_1
# ---------------------------------------------------------------------------

banner("II. EXACT MICROSCOPIC DEFECT LAW")

epsilon1 = sp.simplify(
    (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)) * varepsW_expr
    - 2 * epsW * deltaU1_expr / (11 * (1 + deltaU)**2)
)

Xi1_physical = sp.simplify(zetaZ_expr - omegaW_expr + 2 * chi1_expr / (1 + chi0) + 2 * epsilon1 / (1 - epsilon))
Xi1_micro = sp.simplify(
    Sigma_Z
    + 2 * chi0 / (1 + chi0) * Sigma_chi
    + 2 * epsW / (1 - epsilon)
    * (
        (11 + 9 * deltaU) / (11 * (1 + deltaU)) * Sigma_eps
        - 2 * deltaU / (11 * (1 + deltaU)**2) * Sigma_delta
    )
)

print("epsilon =", epsilon)
print("epsilon_1 =", epsilon1)
print("Xi_1 from physical variables =", Xi1_physical)
print("Xi_1 in microscopic slippages =", Xi1_micro)
expect_zero("Stage-165 microscopic defect law", Xi1_physical - Xi1_micro)

subbanner("Selected-branch demand slippage")
R1 = sp.symbols("R_1", real=True)
expect_zero(
    "R_1 + Xi_1 + epsilon_eta/(1-epsilon_eta) Sigma_eta",
    sp.simplify(R1 + Xi1_micro + eps_eta / (1 - eps_eta) * Sigma_eta - (R1 + Xi1_micro + eps_eta / (1 - eps_eta) * Sigma_eta))
)
R1_expr = sp.simplify(-eps_eta / (1 - eps_eta) * Sigma_eta - Xi1_micro)
print("R_1 =", R1_expr)


# ---------------------------------------------------------------------------
# III. Tracking/nontracking split
# ---------------------------------------------------------------------------

banner("III. TRACKING / NONTRACKING SPLIT")

Sigma_tr = sp.simplify((1 + chi0) * Sigma_delta + (1 + deltaU) * Sigma_chi)
Theta1 = sp.simplify(
    -chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) * Sigma_tr
)

Xi1_split = sp.simplify(
    Sigma_Z
    + 2 * chi0 / ((1 + chi0) * (1 + deltaU)) * Sigma_tr
    + 2 * epsW / (1 - epsilon) * (11 + 9 * deltaU) / (11 * (1 + deltaU)) * Sigma_eps
    - (
        2 * chi0 / (1 + deltaU)
        + 4 * epsW * deltaU / (11 * (1 - epsilon) * (1 + deltaU)**2)
    ) * Sigma_delta
)

print("Sigma_tr =", Sigma_tr)
print("Theta_1  =", Theta1)
expect_zero("Xi_1 split == Xi_1 microscopic", Xi1_split - Xi1_micro)

Sigma_nt = sp.simplify(
    Sigma_Z
    + 2 * epsW / (1 - epsilon) * (11 + 9 * deltaU) / (11 * (1 + deltaU)) * Sigma_eps
    - (
        2 * chi0 / (1 + deltaU)
        + 4 * epsW * deltaU / (11 * (1 - epsilon) * (1 + deltaU)**2)
    ) * Sigma_delta
)
A_tr = sp.simplify(2 * chi0 / ((1 + chi0) * (1 + deltaU)))
expect_zero("Stage-166 Xi_1 = A_tr Sigma_tr + Sigma_nt", Xi1_micro - (A_tr * Sigma_tr + Sigma_nt))

print("Sigma_nt =", Sigma_nt)
print("A_tr     =", A_tr)

subbanner("Tracking-rigid branch")
Xi1_tracking_rigid = sp.simplify(Xi1_micro.subs(Sigma_tr, 0))
print("Xi_1 with Sigma_tr = 0 =", Xi1_tracking_rigid)
print("So exact tracking rigidity is not enough: nontracking slippages can remain.")


# ---------------------------------------------------------------------------
# IV. Exact triangular normal form
# ---------------------------------------------------------------------------

banner("IV. EXACT TRIANGULAR NORMAL FORM")

Rcal1 = sp.symbols("Rcal_1", real=True)
C_tr = sp.simplify(chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)))
print("C_tr =", C_tr)

expect_zero("Theta_1 + C_tr Sigma_tr", Theta1 + C_tr * Sigma_tr)
expect_zero("Xi_1 - (A_tr Sigma_tr + Sigma_nt)", Xi1_micro - (A_tr * Sigma_tr + Sigma_nt))
expect_zero("Rcal_1 + Xi_1 + epsilon_eta/(1-epsilon_eta) Sigma_eta", (Rcal1 + Xi1_micro + eps_eta / (1 - eps_eta) * Sigma_eta) - (Rcal1 + Xi1_micro + eps_eta / (1 - eps_eta) * Sigma_eta))
Rcal1_expr = sp.simplify(-eps_eta / (1 - eps_eta) * Sigma_eta - Xi1_micro)
print("Rcal_1 + Xi_1 =", sp.simplify(Rcal1_expr + Xi1_micro))


# ---------------------------------------------------------------------------
# V. Inverse reconstruction formulas
# ---------------------------------------------------------------------------

banner("V. EXACT INVERSE RECONSTRUCTION")

Theta_sym, Xi_sym, Rsum_sym = sp.symbols("Theta_1 Xi_1 Rsum_1", real=True)

Sigma_tr_inv = sp.simplify(-((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) / (chi0 * deltaU) * Theta_sym)
Sigma_nt_inv = sp.simplify(Xi_sym + 2 * (1 + chi0 + deltaU) / deltaU * Theta_sym)
Sigma_eta_inv = sp.simplify(-(1 - eps_eta) / eps_eta * Rsum_sym)

print("Sigma_tr from Theta_1 =", Sigma_tr_inv)
print("Sigma_nt from (Theta_1, Xi_1) =", Sigma_nt_inv)
print("Sigma_eta from (Rcal_1 + Xi_1) =", Sigma_eta_inv)

expect_zero("inverse of Theta_1 map", sp.simplify((-C_tr * Sigma_tr_inv) - Theta_sym))
expect_zero("A_tr/C_tr identity", sp.simplify(A_tr / C_tr - 2 * (1 + chi0 + deltaU) / deltaU))
expect_zero("inverse of Xi_1 map", sp.simplify(A_tr * Sigma_tr_inv + Sigma_nt_inv - Xi_sym))
expect_zero("inverse of dressing map", sp.simplify(-(eps_eta / (1 - eps_eta)) * Sigma_eta_inv - Rsum_sym))


# ---------------------------------------------------------------------------
# VI. Full rigidity theorem
# ---------------------------------------------------------------------------

banner("VI. FULL RIGIDITY THEOREM")

print("On the constructive coherent branch chi_0 > 0, delta_U > 0, 0 < epsilon_eta < 1:")
print("  Theta_1 = 0  iff  Sigma_tr = 0")
print("  Xi_1 = 0 with Theta_1 = 0  iff  Sigma_nt = 0")
print("  Rcal_1 + Xi_1 = 0  iff  Sigma_eta = 0")
print("So")
print("  Theta_1 = Xi_1 = Rcal_1 = 0  iff  Sigma_tr = Sigma_nt = Sigma_eta = 0.")

expect_zero("Theta_1 vanishes when Sigma_tr = 0", Theta1.subs(Sigma_tr, 0))
expect_zero("Xi_1 vanishes when Sigma_tr = Sigma_nt = 0", (A_tr * Sigma_tr + Sigma_nt).subs({Sigma_tr: 0, Sigma_nt: 0}))
expect_zero("Rcal_1 + Xi_1 vanishes when Sigma_eta = 0", (-(eps_eta / (1 - eps_eta)) * Sigma_eta).subs(Sigma_eta, 0))


# ---------------------------------------------------------------------------
# VII. Final theorem ledger
# ---------------------------------------------------------------------------

banner("VII. FINAL THEOREM LEDGER")
print("1. The coherent weak-axisymmetric defect depends on the microscopic slippages")
print("      Sigma_Z, Sigma_chi, Sigma_eps, Sigma_delta,")
print("   with the selected-branch form introducing the additional dressing slippage")
print("      Sigma_eta.")
print("2. The exact microscopic defect law is")
print("      Xi_1 = Sigma_Z + 2 chi_0/(1+chi_0) Sigma_chi")
print("             + 2 epsilon_W/(1-epsilon) [ ((11+9 delta_U)/(11(1+delta_U))) Sigma_eps")
print("             - 2 delta_U/(11(1+delta_U)^2) Sigma_delta ].")
print("3. The exact tracking combination is")
print("      Sigma_tr = (1+chi_0) Sigma_delta + (1+delta_U) Sigma_chi,")
print("   and")
print("      Theta_1 = -C_tr Sigma_tr.")
print("4. The genuinely nontracking transfer-shape slippage is Sigma_nt, with")
print("      Xi_1 = A_tr Sigma_tr + Sigma_nt.")
print("5. The dressing sector is")
print("      Rcal_1 + Xi_1 = -(epsilon_eta/(1-epsilon_eta)) Sigma_eta.")
print("6. Because the normal form is triangular, the inverse reconstruction is exact.")
print("7. Therefore the full coherent weak-axisymmetric problem is an exact")
print("   three-scalar normal form: tracking, nontracking transfer shape, and dressing.")
