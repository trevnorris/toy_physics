
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

def expect_zero(name: str, expr):
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


"""
5pn_stage165_microscopic_coherent_slippage.py

Executable SymPy audit for the moving-throat Stage 165 theorem.

What this script does
---------------------
1. Pushes the coherent-branch defect law down to microscopic coherent-kernel
   couplings.
2. Defines the exact microscopic slippages
      Sigma_Z, Sigma_chi, Sigma_eta, Sigma_epsilon, Sigma_delta.
3. Proves the Stage-165 microscopic defect law for Xi_1.
4. Extracts the exact tracking combination
      Sigma_tr = (1+chi_0) Sigma_delta + (1+delta_U) Sigma_chi
   and the exact tracking/nontracking split.
"""

banner("STAGE 165 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE DECOMPOSITION")

# ---------------------------------------------------------------------------
# I. Microscopic slippage variables
# ---------------------------------------------------------------------------

subbanner("I. Microscopic coherent-kernel slippage variables")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
epsW = sp.symbols("epsilon_W", positive=True, real=True)
eps_eta = sp.symbols("epsilon_eta", positive=True, real=True)
epsilon = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))

lambda1, c1, gamma1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
kappaU, kappa_eta, kappaW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)

Sigma_Z = sp.simplify(2 * lambda1 + mu1 - kappa_eta - 2 * kappaW)
Sigma_chi = sp.simplify(gamma1 + c1 - kappaU)
Sigma_eta = sp.simplify(2 * c1 - kappaU - kappa_eta)
Sigma_epsilon = sp.simplify(2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
Sigma_delta = sp.simplify(tau1 - kappaU)

print("Sigma_Z =", Sigma_Z)
print("Sigma_chi =", Sigma_chi)
print("Sigma_eta =", Sigma_eta)
print("Sigma_epsilon =", Sigma_epsilon)
print("Sigma_delta =", Sigma_delta)

# Direct physical drifts
zetaZ = sp.simplify(2 * lambda1 - kappa_eta - kappaW)
omegaW = sp.simplify(kappaW - mu1)
chi1 = sp.simplify(chi0 * Sigma_chi)
eta1 = sp.simplify(eps_eta * Sigma_eta)
varepsW = sp.simplify(epsW * Sigma_epsilon)
deltaU1 = sp.simplify(deltaU * Sigma_delta)

expect_zero("Sigma_Z - (zeta_Z - omega_W)", Sigma_Z - (zetaZ - omegaW))
expect_zero("chi_1 - chi_0 Sigma_chi", chi1 - chi0 * Sigma_chi)
expect_zero("eta_1 - epsilon_eta Sigma_eta", eta1 - eps_eta * Sigma_eta)
expect_zero("varepsilon_W - epsilon_W Sigma_epsilon", varepsW - epsW * Sigma_epsilon)
expect_zero("delta_U1 - delta_U Sigma_delta", deltaU1 - deltaU * Sigma_delta)

# ---------------------------------------------------------------------------
# II. Exact microscopic defect law
# ---------------------------------------------------------------------------

subbanner("II. Exact microscopic defect law for Xi_1")

epsilon1 = sp.simplify(
    (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)) * varepsW
    - 2 * epsW * deltaU1 / (11 * (1 + deltaU)**2)
)

Xi1_physical = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * epsilon1 / (1 - epsilon))
Xi1_micro = sp.simplify(
    Sigma_Z
    + 2 * chi0 / (1 + chi0) * Sigma_chi
    + 2 * epsW / (1 - epsilon) * (
        (11 + 9 * deltaU) / (11 * (1 + deltaU)) * Sigma_epsilon
        - 2 * deltaU / (11 * (1 + deltaU)**2) * Sigma_delta
    )
)

print("epsilon =", epsilon)
print("epsilon_1 =", epsilon1)
print("Xi_1 from physical branch variables =", Xi1_physical)
print("Xi_1 in microscopic slippages =", Xi1_micro)
expect_zero("Stage-165 microscopic defect law", Xi1_physical - Xi1_micro)

R1 = sp.symbols("R_1", real=True)
R1_expr = sp.simplify(-eps_eta / (1 - eps_eta) * Sigma_eta - Xi1_micro)
print("R_1 =", R1_expr)

# ---------------------------------------------------------------------------
# III. Tracking/nontracking split
# ---------------------------------------------------------------------------

subbanner("III. Exact tracking/nontracking split")

Sigma_tr = sp.simplify((1 + chi0) * Sigma_delta + (1 + deltaU) * Sigma_chi)
Theta1 = sp.simplify(
    -chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) * Sigma_tr
)
Xi1_split = sp.simplify(
    Sigma_Z
    + 2 * chi0 / ((1 + chi0) * (1 + deltaU)) * Sigma_tr
    + 2 * epsW / (1 - epsilon) * (11 + 9 * deltaU) / (11 * (1 + deltaU)) * Sigma_epsilon
    - (
        2 * chi0 / (1 + deltaU)
        + 4 * epsW * deltaU / (11 * (1 - epsilon) * (1 + deltaU)**2)
    ) * Sigma_delta
)

print("Sigma_tr =", Sigma_tr)
print("Theta_1 =", Theta1)
expect_zero("Xi_1 split == Xi_1 microscopic", Xi1_split - Xi1_micro)

print("With Sigma_tr = 0, the grouped defect survives as a nontracking combination of")
print("Sigma_Z, Sigma_epsilon, and Sigma_delta.")

banner("FINAL STAGE-165 LEDGER")
print("1. The coherent weak-axisymmetric defect depends on five microscopic slippages:")
print("      Sigma_Z, Sigma_chi, Sigma_eta, Sigma_epsilon, Sigma_delta.")
print("2. The exact defect law is")
print("      Xi_1 = Sigma_Z + 2 chi_0/(1+chi_0) Sigma_chi")
print("             + 2 epsilon_W/(1-epsilon)[((11+9 delta_U)/(11(1+delta_U))) Sigma_epsilon")
print("             - 2 delta_U/(11(1+delta_U)^2) Sigma_delta].")
print("3. The exact tracking combination is")
print("      Sigma_tr = (1+chi_0) Sigma_delta + (1+delta_U) Sigma_chi.")
print("4. The tracking factor drift is Theta_1 = -C_tr Sigma_tr.")
print("5. Exact tracking rigidity is necessary but not sufficient: the grouped defect also carries nontracking microscopic slippages.")
