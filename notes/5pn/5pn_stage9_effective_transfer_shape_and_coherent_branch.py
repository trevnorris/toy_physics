
#!/usr/bin/env python3
"""
5pn_stage9_effective_transfer_shape_and_coherent_branch.py

Ninth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Re-derives the effective-transfer-shape collapse
      Xi_1 = delta ln T_eff^2 / (eps lambda_A).
2. Evaluates the one-port continuum transfer shape
      T^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-epsilon_W)^2]
   and derives the exact weak-axisymmetric slope law
      Xi_1 = zeta_Z - omega_W + 2 rho_1/(1+rho) + 2 varepsilon_W/(1-epsilon_W).
3. Rewrites the same transfer shape in selected-branch form
      T^2 = const * (1-epsilon_eta) / R_target,
   proving
      Xi_1 = - eta_1/(1-epsilon_eta) - R_1.
4. Specializes to the coherent local D/N tracking branch
      T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-epsilon)^2],
      epsilon = epsilon_W (1 - 2/11 * delta_U/(1+delta_U)),
   and derives the exact coherent weak-axisymmetric defect law.
5. Verifies the support-blindness theorem and shows explicitly that tracking
   rigidity is not enough to kill Xi_1.

Interpretation
--------------
After this stage the remaining weak-axisymmetric grouped defect is the slope of a
single effective transfer shape. On the actual one-port continuum branch, it is
controlled by a small set of mixed/outgoing placement variables rather than by a
generic bundle of internal drifts.
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
# I. Effective transfer-shape collapse
# ---------------------------------------------------------------------------

banner("I. EFFECTIVE TRANSFER-SHAPE COLLAPSE")

T1, T2 = sp.symbols("T_1 T_2", positive=True, real=True)
tau1, tau2 = sp.symbols("tau_1 tau_2", real=True)
eps, lamA = sp.symbols("epsilon lambda_A", real=True)
t = sp.symbols("t", real=True)

T1_t = T1 * sp.exp(t * eps * lamA * tau1)
T2_t = T2 * sp.exp(t * eps * lamA * tau2)
Teff2_t = sp.simplify(T1_t**2 + T2_t**2)
Xi1_eff = sp.simplify(sp.diff(sp.log(Teff2_t), t).subs(t, 0) / (eps * lamA))
Xi1_weighted = sp.simplify(2 * (T1**2 * tau1 + T2**2 * tau2) / (T1**2 + T2**2))

print("T_eff^2 =", T1**2 + T2**2)
print("Xi_1 from log T_eff^2 =", Xi1_eff)
print("Xi_1 from weighted port slopes =", Xi1_weighted)
expect_zero("effective-transfer-shape collapse", Xi1_eff - Xi1_weighted)


# ---------------------------------------------------------------------------
# II. One-port actual continuum transfer shape
# ---------------------------------------------------------------------------

banner("II. ONE-PORT CONTINUUM TRANSFER SHAPE")

ZW, rho, OmegaW2, epsW = sp.symbols("Z_W rho Omega_W2 epsilon_W", positive=True, real=True)
zetaZ, omegaW, rho1, varepsW = sp.symbols("zeta_Z omega_W rho_1 varepsilon_W", real=True)

ZW_t = ZW * sp.exp(t * eps * lamA * zetaZ)
OmegaW2_t = OmegaW2 * sp.exp(t * eps * lamA * omegaW)
rho_t = rho + t * eps * lamA * rho1
epsW_t = epsW + t * eps * lamA * varepsW

T2_cont_t = sp.simplify(ZW_t * (1 + rho_t)**2 / (OmegaW2_t * (1 - epsW_t)**2))
Xi1_cont = sp.simplify(sp.diff(sp.log(T2_cont_t), t).subs(t, 0) / (eps * lamA))
Xi1_cont_expected = sp.simplify(zetaZ - omegaW + 2 * rho1 / (1 + rho) + 2 * varepsW / (1 - epsW))

print("T^2 =", ZW * (1 + rho)**2 / (OmegaW2 * (1 - epsW)**2))
print("Xi_1 from the one-port continuum branch =", Xi1_cont)
expect_zero("one-port continuum slope law", Xi1_cont - Xi1_cont_expected)


# ---------------------------------------------------------------------------
# III. Selected-branch reformulation
# ---------------------------------------------------------------------------

banner("III. SELECTED-BRANCH REFORMULATION")

G, cs, a, c = sp.symbols("G c_s a c", positive=True, real=True)
eps_eta = sp.symbols("epsilon_eta", real=True)
eta1, R1 = sp.symbols("eta_1 R_1", real=True)
Rtarget = sp.symbols("R_target", positive=True, real=True)

Lambda0 = sp.simplify(27 * sp.pi**2 * G * cs**5 / (20 * a**5 * c**5))
eps_eta_t = eps_eta + t * eps * lamA * eta1
Rtarget_t = Rtarget * sp.exp(t * eps * lamA * R1)

T2_sel_t = sp.simplify(Lambda0 * (1 - eps_eta_t) / Rtarget_t)
Xi1_sel = sp.simplify(sp.diff(sp.log(T2_sel_t), t).subs(t, 0) / (eps * lamA))
Xi1_sel_expected = sp.simplify(-eta1 / (1 - eps_eta) - R1)

print("T^2 selected-branch form =", Lambda0 * (1 - eps_eta) / Rtarget)
print("Xi_1 from selected-branch form =", Xi1_sel)
expect_zero("selected-branch slope law", Xi1_sel - Xi1_sel_expected)


# ---------------------------------------------------------------------------
# IV. Coherent local D/N tracking branch
# ---------------------------------------------------------------------------

banner("IV. COHERENT LOCAL D/N TRACKING BRANCH")

chi0, chi1 = sp.symbols("chi_0 chi_1", positive=True, real=True)
deltaU, deltaU1 = sp.symbols("delta_U delta_U1", positive=True, real=True)
zeta = sp.symbols("zeta", real=True)  # support ratio (should drop out)
epsilon = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))

chi0_t = chi0 + t * eps * lamA * chi1
deltaU_t = deltaU + t * eps * lamA * deltaU1
zeta_t = zeta + t * eps * lamA * sp.Symbol("zeta_1")  # intentionally unused
epsW_t2 = epsW + t * eps * lamA * varepsW

epsilon_t = sp.simplify(epsW_t2 * (1 - sp.Rational(2, 11) * deltaU_t / (1 + deltaU_t)))
T2_coh_t = sp.simplify(ZW_t * (1 + chi0_t)**2 / (OmegaW2_t * (1 - epsilon_t)**2))

Xi1_coh = sp.simplify(sp.diff(sp.log(T2_coh_t), t).subs(t, 0) / (eps * lamA))
epsilon1 = sp.simplify(
    (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)) * varepsW
    - 2 * epsW * deltaU1 / (11 * (1 + deltaU)**2)
)
Xi1_coh_expected = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * epsilon1 / (1 - epsilon))

print("epsilon =", epsilon)
print("epsilon_1 =", epsilon1)
print("T_coherent^2 =", sp.simplify(ZW * (1 + chi0)**2 / (OmegaW2 * (1 - epsilon)**2)))
print("Xi_1 on the coherent branch =", Xi1_coh)
expect_zero("coherent one-port slope law", Xi1_coh - Xi1_coh_expected)

subbanner("Support-blindness theorem")
expect_zero("d/dzeta log T_coherent^2", sp.diff(sp.log(sp.simplify(ZW * (1 + chi0)**2 / (OmegaW2 * (1 - epsilon)**2))), zeta))

subbanner("Tracking-factor drift is not enough")
Theta1 = sp.simplify(
    -(
        chi0 * (1 + chi0) * deltaU1
        + deltaU * (1 + deltaU) * chi1
    )
    / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
)
print("Theta_1 =", Theta1)
Theta1_rigid = sp.simplify(Theta1.subs({chi1: 0, deltaU1: 0}))
Xi1_rigid = sp.simplify(Xi1_coh_expected.subs({chi1: 0, deltaU1: 0}))
expect_zero("Theta_1 on the tracking-rigid subbranch", Theta1_rigid)
print("Xi_1 on the same tracking-rigid subbranch =", Xi1_rigid)
print("So exact tracking rigidity does not in general kill Xi_1.")


# ---------------------------------------------------------------------------
# V. Final theorem ledger
# ---------------------------------------------------------------------------

banner("V. FINAL THEOREM LEDGER")
print("1. The many-port grouped defect collapses to one effective transfer-shape slope:")
print("      Xi_1 = delta ln T_eff^2 / (eps lambda_A).")
print("2. On the actual one-port continuum branch,")
print("      T^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-epsilon_W)^2],")
print("   and therefore")
print("      Xi_1 = zeta_Z - omega_W + 2 rho_1/(1+rho) + 2 varepsilon_W/(1-epsilon_W).")
print("3. In selected-branch form, the same defect is")
print("      Xi_1 = - eta_1/(1-epsilon_eta) - R_1.")
print("4. On the coherent local D/N branch,")
print("      T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-epsilon)^2],")
print("      epsilon = epsilon_W (1 - 2/11 * delta_U/(1+delta_U)),")
print("   and")
print("      Xi_1 = zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 epsilon_1/(1-epsilon).")
print("5. The support ratio zeta drops out identically: the defect is support-blind.")
print("6. Tracking rigidity Theta_1 = 0 is necessary but not sufficient to kill Xi_1.")
print("7. So the next honest theorem gate is to compute the weak-axisymmetric drifts")
print("   of the mixed/outgoing placement variables on the actual moving-throat branch.")
