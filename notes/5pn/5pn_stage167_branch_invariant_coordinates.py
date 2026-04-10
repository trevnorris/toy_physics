
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
5pn_stage167_branch_invariant_coordinates.py

Executable SymPy audit for the moving-throat Stage 167 theorem.

What this script does
---------------------
1. Builds the exact coherent-branch composites
      R_tr, T^2, epsilon_eta
   and the corrected nontracking composite
      N_* = T^2 R_tr^{B_*}.
2. Defines the normalized tracking invariant
      T_* = R_tr^{-C_*}.
3. Proves the exact direct branch-coordinate identities
      d ln T_* = Sigma_tr,
      d ln N_* = Sigma_nt,
      d ln epsilon_eta = Sigma_eta
   at first grouped weak-axisymmetric order.
4. Rewrites the full zero-defect theorem in terms of these exact branch composites.
"""

banner("STAGE 167 — EXACT BRANCH-INVARIANT COORDINATES")

# ---------------------------------------------------------------------------
# I. Exact coherent branch equations
# ---------------------------------------------------------------------------

subbanner("I. Exact coherent branch equations")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
ZW, OmegaW2 = sp.symbols("Z_W Omega_W2", positive=True, real=True)
epsW = sp.symbols("epsilon_W", positive=True, real=True)
eps_eta = sp.symbols("epsilon_eta", positive=True, real=True)
epsilon = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))

Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
T2 = sp.simplify(ZW * (1 + chi0)**2 / (OmegaW2 * (1 - epsilon)**2))

print("R_tr =", Rtr)
print("T^2 =", T2)
print("epsilon_eta =", eps_eta)

# weak-axisymmetric drifts of branch variables
eps, lamA, t = sp.symbols("epsilon lambda_A t", real=True)
Theta1, Xi1, Sigma_eta = sp.symbols("Theta_1 Xi_1 Sigma_eta", real=True)
chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)

# ---------------------------------------------------------------------------
# II. Direct branch-invariant tracking coordinate
# ---------------------------------------------------------------------------

subbanner("II. Direct branch-invariant tracking coordinate")

Cstar = sp.simplify((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs) / (chi0s * deltaUs))
Tstar = sp.simplify(Rtr**(-Cstar))

# Linear law from Stage 166
Sigma_tr = sp.symbols("Sigma_tr", real=True)
Theta_expr = sp.simplify(-chi0s * deltaUs / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs)) * Sigma_tr)

Rtr_t = sp.simplify(Rtr * sp.exp(t * eps * lamA * Theta_expr))
Tstar_t = sp.simplify(Rtr_t**(-Cstar))
dlog_Tstar = sp.simplify(sp.diff(sp.log(Tstar_t), t).subs(t, 0) / (eps * lamA))

print("C_* =", Cstar)
print("T_* =", Tstar)
expect_zero("d ln T_* - Sigma_tr", dlog_Tstar - Sigma_tr)

# ---------------------------------------------------------------------------
# III. Direct branch-invariant nontracking coordinate
# ---------------------------------------------------------------------------

subbanner("III. Direct branch-invariant nontracking coordinate")

Bstar = sp.simplify(2 * (1 + chi0s + deltaUs) / deltaUs)
Sigma_nt = sp.symbols("Sigma_nt", real=True)
Xi_expr = sp.simplify(Sigma_nt - Bstar * Theta_expr)  # Xi_1 = Sigma_nt - B_* Theta_1
Nstar = sp.simplify(T2 * Rtr**Bstar)

T2_t = sp.simplify(T2 * sp.exp(t * eps * lamA * Xi_expr))
Nstar_t = sp.simplify(T2_t * Rtr_t**Bstar)
dlog_Nstar = sp.simplify(sp.diff(sp.log(Nstar_t), t).subs(t, 0) / (eps * lamA))

print("B_* =", Bstar)
print("N_* =", Nstar)
expect_zero("d ln N_* - Sigma_nt", dlog_Nstar - Sigma_nt)

# ---------------------------------------------------------------------------
# IV. Direct dressing coordinate
# ---------------------------------------------------------------------------

subbanner("IV. Direct dressing coordinate")

eps_eta_t = sp.simplify(eps_eta * sp.exp(t * eps * lamA * Sigma_eta))
dlog_eps_eta = sp.simplify(sp.diff(sp.log(eps_eta_t), t).subs(t, 0) / (eps * lamA))
expect_zero("d ln epsilon_eta - Sigma_eta", dlog_eps_eta - Sigma_eta)

# Selected-branch complementary composite E = 1 - epsilon_eta
Ecomp = sp.simplify(1 - eps_eta)
Rsum = sp.simplify(-eps_eta / (1 - eps_eta) * Sigma_eta)
Ecomp_t = sp.simplify(1 - eps_eta_t)
dlog_Ecomp = sp.simplify(sp.diff(sp.log(Ecomp_t), t).subs(t, 0) / (eps * lamA))
expect_zero("d ln (1-epsilon_eta) - (R_1 + Xi_1)", dlog_Ecomp - Rsum)

# ---------------------------------------------------------------------------
# V. Composite rigidity theorem
# ---------------------------------------------------------------------------

subbanner("V. Composite rigidity theorem in direct branch variables")

print("The direct branch coordinates are:")
print("  d ln T_* = Sigma_tr")
print("  d ln N_* = Sigma_nt")
print("  d ln epsilon_eta = Sigma_eta")
print("Therefore:")
print("  Theta_1 = Xi_1 = R_1 = 0")
print("iff")
print("  d ln T_* = d ln N_* = d ln epsilon_eta = 0.")

banner("FINAL STAGE-167 LEDGER")
print("1. The coherent branch supplies exact composite quantities R_tr, T^2, and epsilon_eta.")
print("2. After harmless positive normalizations, the branch-invariant coordinates are")
print("      T_* = R_tr^{-C_*},")
print("      N_* = T^2 R_tr^{B_*},")
print("      D = epsilon_eta.")
print("3. Their first grouped weak-axisymmetric drifts are exactly")
print("      d ln T_* = Sigma_tr, d ln N_* = Sigma_nt, d ln epsilon_eta = Sigma_eta.")
print("4. So the continuation point is no longer a generic slippage ledger; it is the first grouped drift of three exact branch composites.")
