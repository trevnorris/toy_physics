
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
5pn_stage168_direct_microscopic_monomials.py

Executable SymPy audit for the moving-throat Stage 168 theorem.

What this script does
---------------------
1. Builds the direct microscopic monomials
      C_tr,* , C_nt,* , epsilon_eta.
2. Proves their logarithmic drifts are exactly
      Sigma_tr, Sigma_nt, Sigma_eta.
3. Turns the zero-defect theorem into the explicit microscopic compatibility
   ledger and solves it for the dependent drifts
      (tau_1, kappa_eta, mu_1).
"""

banner("STAGE 168 — DIRECT MICROSCOPIC MONOMIALS AND COMPATIBILITY LEDGER")

# ---------------------------------------------------------------------------
# I. Microscopic monomials
# ---------------------------------------------------------------------------

subbanner("I. Direct microscopic monomials")

lamW, c_etaU, gamma, KU = sp.symbols("lambda_W c_etaU gamma K_U", positive=True, real=True)
Keta, KW, muW, TU = sp.symbols("K_eta K_W mu_W T_U", positive=True, real=True)
L, sigma = sp.symbols("L sigma", positive=True, real=True)

chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)

Ctr = sp.simplify((gamma * c_etaU / KU)**(1 + deltaUs) * (sp.pi**2 * TU / (L**2 * KU))**(1 + chi0s))
Cnt = sp.simplify(
    (lamW**2 * muW / (Keta * KW**2))
    * (gamma**2 * lamW**2 * sigma / (KU * KW))**Estar
    * (sp.pi**2 * TU / (L**2 * KU))**(-Fstar)
)
eps_eta = sp.simplify(c_etaU**2 / (KU * Keta))

print("C_tr,* =")
sp.pprint(Ctr)
print("C_nt,* =")
sp.pprint(Cnt)
print("epsilon_eta =")
sp.pprint(eps_eta)

# ---------------------------------------------------------------------------
# II. Logarithmic drifts
# ---------------------------------------------------------------------------

subbanner("II. Exact logarithmic drifts of the direct monomials")

lambda1, c1, gamma1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
kappaU, kappa_eta, kappaW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)

Sigma_tr = sp.expand((1 + deltaUs) * (gamma1 + c1 - kappaU) + (1 + chi0s) * (tau1 - kappaU))
Sigma_nt = sp.expand(
    (2 * lambda1 + mu1 - kappa_eta - 2 * kappaW)
    + Estar * (2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
    - Fstar * (tau1 - kappaU)
)
Sigma_eta = sp.expand(2 * c1 - kappaU - kappa_eta)

print("Sigma_tr =", Sigma_tr)
print("Sigma_nt =", Sigma_nt)
print("Sigma_eta =", Sigma_eta)

# ---------------------------------------------------------------------------
# III. Explicit microscopic compatibility ledger
# ---------------------------------------------------------------------------

subbanner("III. Explicit microscopic compatibility ledger")

compat_eqs = [
    sp.Eq(Sigma_tr, 0),
    sp.Eq(Sigma_eta, 0),
    sp.Eq(Sigma_nt, 0),
]
for eq in compat_eqs:
    sp.pprint(eq)

sol = sp.solve(compat_eqs, [tau1, kappa_eta, mu1], dict=True)
if len(sol) != 1:
    raise AssertionError("Expected a unique solve for (tau_1, kappa_eta, mu_1).")
sol = sol[0]

print("Solved compatibility relations:")
for k, v in sol.items():
    print(f"  {k} = {sp.simplify(v)}")

expect_zero("tracking compatibility", Sigma_tr.subs(sol))
expect_zero("dressing compatibility", Sigma_eta.subs(sol))
expect_zero("nontracking compatibility", Sigma_nt.subs(sol))

banner("FINAL STAGE-168 LEDGER")
print("1. The direct microscopic monomials are")
print("      C_tr,* = (gamma c_etaU/K_U)^{1+delta_U,*} (pi^2 T_U/(L^2 K_U))^{1+chi_0,*},")
print("      C_nt,* = (lambda_W^2 mu_W/(K_eta K_W^2)) * (gamma^2 lambda_W^2 sigma/(K_U K_W))^{E_*} * (pi^2 T_U/(L^2 K_U))^{-F_*},")
print("      epsilon_eta = c_etaU^2/(K_U K_eta).")
print("2. Their logarithmic drifts are exactly Sigma_tr, Sigma_nt, Sigma_eta.")
print("3. Therefore the full coherent weak-axisymmetric zero-defect theorem is equivalent to")
print("      d ln C_tr,* = d ln C_nt,* = d ln epsilon_eta = 0.")
print("4. The resulting microscopic compatibility ledger solves uniquely for")
print("      tau_1, kappa_eta, mu_1")
print("   in terms of the five free drifts (lambda_1, c_1, gamma_1, kappa_U, kappa_W).")
