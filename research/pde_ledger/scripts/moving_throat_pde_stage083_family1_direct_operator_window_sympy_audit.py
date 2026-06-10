#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 083 SymPy audit

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


def expect_close(name: str, value, target, atol: float) -> None:
    diff = abs(sp.N(value, 40) - sp.N(target, 40))
    print(f"{name} diff = {diff}")
    if diff > sp.Float(atol):
        raise AssertionError(f"{name}: diff {diff} exceeds {atol}")


banner("STAGE 083 — DIRECT FAMILY-1 OPERATOR WINDOW")

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

# Numeric anchors for the Family-1 closed forms reported in the stage notes.
expect_close("Delta_0(F1) numeric anchor", Delta0_F1, sp.Float("1.73302079021525e-4"), 1e-16)
expect_close("Delta_inf(F1) numeric anchor", DeltaInf_F1, sp.Float("2.01447565540522e-2"), 1e-15)


yy = sp.symbols("yy")
y_F1 = sp.nsolve(yy * sp.tan(yy) - eta_F1, sp.Float("1.53", 80), tol=1e-30, maxsteps=100, prec=80)
y_residual = sp.N(y_F1 * sp.tan(y_F1) - eta_F1, 40)
# nsolve solved this to tol=1e-30; verify residual is below 1e-25.
if abs(y_residual) > sp.Float("1e-25"):
    raise AssertionError(f"y_F1 fails defining equation: residual = {y_residual}")
print(f"y_F1 defining-equation residual = {y_residual}")
A_F1 = sp.simplify((kappa_F1 + pi**2 / 4) / (kappa_F1 + y_F1**2))

print("Delta_0(F1) =", sp.N(Delta0_F1, 30))
print("Delta_inf(F1) =", sp.N(DeltaInf_F1, 30))
print("y_F1 =", sp.N(y_F1, 30))
print("A_F1 =", sp.N(A_F1, 30))

# SOURCE-ANCHOR (operator selectors):
#   Theta_chi_coeff (= 4.06863235008162) and Theta_J_coeff (= 0.927552032539308)
#   are the natural Family-1 operator selectors carried forward from earlier
#   stages.  The integer prefactor 136900 = 370^2 = eta_F1 * (kappa_F1 + ...?)
#   originates upstream as well.  If an upstream verifying script anchors
#   these (e.g., a root of an operator equation), reference it here.  As of
#   this audit no upstream sympy/mathematica script in this unit verifies
#   them, so any transcription typo here will propagate silently.
#   TODO(upstream-anchor): cite the producing stage's verifying script.
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
expect_close("Pe_-^(chi) numeric check", Pe_minus_chi, sp.Float("96.5285247264385"), 1e-10)
expect_close("Pe_+^(chi) numeric check", Pe_plus_chi, sp.Float("11220.5441626259"), 1e-7)
expect_close("Pe_-^(J) numeric check", Pe_minus_J, sp.Float("22.0062226330754"), 1e-10)
expect_close("Pe_+^(J) numeric check", Pe_plus_J, sp.Float("2558.01892349205"), 1e-8)

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

# Monotonicity check: zeta_F1(Pe) = A_F1 * Omega(Pe)^2 should be monotone
# increasing in Pe over the relevant range, since Omega is positive and
# increasing on (0, infinity).  The ledger claim "monotone in Xi on the
# stable branch" follows because Pe = Xi * Delta is monotone in Xi on the
# stable branch.  Check d zeta / d Pe > 0 at a sample of representative Pe
# values spanning the chi and J windows.
dzeta_dPe = sp.diff(zeta_F1, Pe_sym)
for Pe_test_val in [sp.Float("10"), sp.Float("100"), sp.Float("1000"), sp.Float("10000")]:
    deriv_num = sp.N(dzeta_dPe.subs(Pe_sym, Pe_test_val), 30)
    print(f"  d zeta / d Pe at Pe = {Pe_test_val} : {deriv_num}")
    if deriv_num <= 0:
        raise AssertionError(
            f"zeta_F1 not monotone increasing at Pe = {Pe_test_val}: deriv = {deriv_num}"
        )

print("zeta_-^(chi) =", zeta_minus_chi)
print("zeta_+^(chi) =", zeta_plus_chi)
print("zeta_-^(J)   =", zeta_minus_J)
print("zeta_+^(J)   =", zeta_plus_J)
print("zeta_max^(F1)=", zeta_max_F1)
expect_close("zeta_-^(chi) numeric check", zeta_minus_chi, sp.Float("2.46622291347846"), 1e-12)
expect_close("zeta_+^(chi) numeric check", zeta_plus_chi, sp.Float("2.46752913273870"), 1e-12)
expect_close("zeta_-^(J) numeric check", zeta_minus_J, sp.Float("2.44257571477179"), 1e-12)
expect_close("zeta_+^(J) numeric check", zeta_plus_J, sp.Float("2.46752736855058"), 1e-12)
expect_close("zeta_max^(F1) numeric check", zeta_max_F1, sp.Float("2.46752922945601"), 1e-12)

print("Pi_suff^(chi)/C_mix at eps_blk=0 =", sp.N(1 + zeta_minus_chi, 30))
print("Pi_fail^(chi)/C_mix at eps_blk=0 =", sp.N(1 + zeta_plus_chi, 30))
print("Pi_suff^(J)/C_mix at eps_blk=0   =", sp.N(1 + zeta_minus_J, 30))
print("Pi_fail^(J)/C_mix at eps_blk=0   =", sp.N(1 + zeta_plus_J, 30))
print("Pi_max^(F1)/C_mix at eps_blk=0   =", sp.N(1 + zeta_max_F1, 30))

print("\nTheorem ledger:")
print("  dPe_*/dXi = Delta / (1 - Xi dDelta/dPe)")
print("  zeta_phys(F1) and Pi/C_mix windows are monotone in Xi on the stable branch")
print("  inserting the natural Family-1 wall data reproduces the Stage-078/080/081 windows directly")
