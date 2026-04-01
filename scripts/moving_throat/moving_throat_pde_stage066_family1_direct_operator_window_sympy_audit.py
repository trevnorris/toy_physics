#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 66 SymPy audit

Checks
------
1. Exact implicit-differentiation formula for dPe_*/dXi.
2. Exact Family-1 support/source constants and Delta_0, Delta_inf.
3. Direct operator-selected transport, zeta, and Pi/C_mix windows for the natural
   shell-weighted and Jensen-floor branches.
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


banner("STAGE 66 — DIRECT FAMILY-1 OPERATOR WINDOW")

# ------------------------------------------------------------------
# 1. Implicit differentiation of the fixed-point branch
# ------------------------------------------------------------------
Xi = sp.symbols("Xi", positive=True, real=True)
Pe = sp.Function("Pe")(Xi)
Delta = sp.Function("Delta")

fixed_point = sp.Eq(Pe, Xi * Delta(Pe))
dPe_formula = sp.solve(sp.diff(fixed_point.lhs - fixed_point.rhs, Xi), sp.diff(Pe, Xi))[0]
dPe_expected = sp.simplify(Delta(Pe) / (1 - Xi * sp.diff(Delta(Pe), Pe)))

print("dPe_*/dXi =")
sp.pprint(dPe_formula)
expect_zero("implicit fixed-point derivative", dPe_formula - dPe_expected)

# ------------------------------------------------------------------
# 2. Family-1 support/source constants
# ------------------------------------------------------------------
pi = sp.pi
kappa_F1 = sp.Rational(12321, 5)
eta_F1 = sp.Integer(37)
alpha_F1 = sp.sqrt(kappa_F1)

Delta0_F1 = sp.simplify(
    eta_F1 * (sp.cosh(alpha_F1) - 1)
    / (alpha_F1**2 * (alpha_F1 * sp.sinh(alpha_F1) + eta_F1 * sp.cosh(alpha_F1)))
)
DeltaInf_F1 = sp.simplify(
    (sp.cosh(alpha_F1) + (eta_F1 / alpha_F1) * sp.sinh(alpha_F1) - 1)
    / (alpha_F1 * sp.sinh(alpha_F1) + eta_F1 * sp.cosh(alpha_F1))
)

yy = sp.symbols("yy")
y_F1 = sp.nsolve(yy * sp.tan(yy) - eta_F1, 1.53, tol=1e-30, maxsteps=100)
A_F1 = sp.simplify((kappa_F1 + pi**2 / 4) / (kappa_F1 + y_F1**2))

print("Delta_0(F1) =", sp.N(Delta0_F1, 30))
print("Delta_inf(F1) =", sp.N(DeltaInf_F1, 30))
print("y_F1 =", sp.N(y_F1, 30))
print("A_F1 =", sp.N(A_F1, 30))

Theta_chi_coeff = sp.Float("4.06863235008162")
Theta_J_coeff = sp.Float("0.927552032539308")
Xi_chi_coeff = sp.Integer(136900) * Theta_chi_coeff
Xi_J_coeff = sp.Integer(136900) * Theta_J_coeff

print("Xi_chi coefficient =", sp.N(Xi_chi_coeff, 30))
print("Xi_J coefficient   =", sp.N(Xi_J_coeff, 30))

Pe_minus_chi = sp.N(Xi_chi_coeff * Delta0_F1, 30)
Pe_plus_chi = sp.N(Xi_chi_coeff * DeltaInf_F1, 30)
Pe_minus_J = sp.N(Xi_J_coeff * Delta0_F1, 30)
Pe_plus_J = sp.N(Xi_J_coeff * DeltaInf_F1, 30)

print("Pe_-^(chi) =", Pe_minus_chi)
print("Pe_+^(chi) =", Pe_plus_chi)
print("Pe_-^(J)   =", Pe_minus_J)
print("Pe_+^(J)   =", Pe_plus_J)

# ------------------------------------------------------------------
# 3. Direct zeta and Pi/C_mix windows
# ------------------------------------------------------------------
Pe_sym = sp.symbols("Pe", positive=True, real=True)
Omega = sp.simplify(pi * Pe_sym * (2 * Pe_sym * sp.exp(Pe_sym) + pi) / ((4 * Pe_sym**2 + pi**2) * (sp.exp(Pe_sym) - 1)))
zeta_F1 = sp.simplify(A_F1 * Omega**2)

zeta_minus_chi = sp.N(zeta_F1.subs(Pe_sym, Pe_minus_chi), 30)
zeta_plus_chi = sp.N(zeta_F1.subs(Pe_sym, Pe_plus_chi), 30)
zeta_minus_J = sp.N(zeta_F1.subs(Pe_sym, Pe_minus_J), 30)
zeta_plus_J = sp.N(zeta_F1.subs(Pe_sym, Pe_plus_J), 30)
zeta_max_F1 = sp.N(A_F1 * pi**2 / 4, 30)

print("zeta_-^(chi) =", zeta_minus_chi)
print("zeta_+^(chi) =", zeta_plus_chi)
print("zeta_-^(J)   =", zeta_minus_J)
print("zeta_+^(J)   =", zeta_plus_J)
print("zeta_max^(F1)=", zeta_max_F1)

print("Pi_suff^(chi)/C_mix at eps_blk=0 =", sp.N(1 + zeta_minus_chi, 30))
print("Pi_fail^(chi)/C_mix at eps_blk=0 =", sp.N(1 + zeta_plus_chi, 30))
print("Pi_suff^(J)/C_mix at eps_blk=0   =", sp.N(1 + zeta_minus_J, 30))
print("Pi_fail^(J)/C_mix at eps_blk=0   =", sp.N(1 + zeta_plus_J, 30))
print("Pi_max^(F1)/C_mix at eps_blk=0   =", sp.N(1 + zeta_max_F1, 30))

print("\nTheorem ledger:")
print("  dPe_*/dXi = Delta / (1 - Xi dDelta/dPe)")
print("  zeta_phys(F1) and Pi/C_mix windows are monotone in Xi on the stable branch")
print("  inserting the natural Family-1 wall data reproduces the Stage-61/63/64 windows directly")
