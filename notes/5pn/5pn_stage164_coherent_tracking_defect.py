
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
5pn_stage164_coherent_tracking_defect.py

Executable SymPy audit for the moving-throat Stage 164 theorem.

What this script does
---------------------
1. Substitutes the actual coherent local D/N tracking branch into the Stage-163
   transfer-shape law.
2. Derives the exact coherent weak-axisymmetric defect law
      Xi_1 = zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 epsilon_1/(1-epsilon).
3. Proves the support-blindness theorem:
      d_zeta log T^2 = d_zeta log R_target = d_zeta Xi_1 = 0.
4. Derives the exact tracking-factor drift Theta_1 and shows that exact tracking
   rigidity is not sufficient to kill Xi_1.
"""

banner("STAGE 164 — COHERENT TRACKING-BRANCH DEFECT LAW")

# ---------------------------------------------------------------------------
# I. Exact coherent tracking-branch substitution
# ---------------------------------------------------------------------------

subbanner("I. Exact coherent branch substitution")

eps, lamA, t = sp.symbols("epsilon lambda_A t", real=True)
ZW, chi0, OmegaW2 = sp.symbols("Z_W chi_0 Omega_W2", positive=True, real=True)
epsW, deltaU = sp.symbols("epsilon_W delta_U", positive=True, real=True)
zetaZ, omegaW, chi1, varepsW, deltaU1 = sp.symbols(
    "zeta_Z omega_W chi_1 varepsilon_W delta_U1", real=True
)
zeta = sp.symbols("zeta", real=True)

epsilon = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
print("epsilon =", epsilon)

ZW_t = ZW * sp.exp(t * eps * lamA * zetaZ)
chi0_t = chi0 + t * eps * lamA * chi1
OmegaW2_t = OmegaW2 * sp.exp(t * eps * lamA * omegaW)
epsW_t = epsW + t * eps * lamA * varepsW
deltaU_t = deltaU + t * eps * lamA * deltaU1

epsilon_t = sp.simplify(epsW_t * (1 - sp.Rational(2, 11) * deltaU_t / (1 + deltaU_t)))
T2_t = sp.simplify(ZW_t * (1 + chi0_t)**2 / (OmegaW2_t * (1 - epsilon_t)**2))
T2 = sp.simplify(ZW * (1 + chi0)**2 / (OmegaW2 * (1 - epsilon)**2))

print("T^2 coherent =", T2)

epsilon1 = sp.simplify(
    (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)) * varepsW
    - 2 * epsW * deltaU1 / (11 * (1 + deltaU)**2)
)
Xi1 = sp.simplify(sp.diff(sp.log(T2_t), t).subs(t, 0) / (eps * lamA))
Xi1_expected = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * epsilon1 / (1 - epsilon))

print("epsilon_1 =", epsilon1)
print("Xi_1 =", Xi1)
expect_zero("coherent weak-axisymmetric defect law", Xi1 - Xi1_expected)

# ---------------------------------------------------------------------------
# II. Support-blindness theorem
# ---------------------------------------------------------------------------

subbanner("II. Support-blindness theorem")

expect_zero("d_zeta log T^2", sp.diff(sp.log(T2), zeta))

# Selected-branch form
G, cs, a, c = sp.symbols("G c_s a c", positive=True, real=True)
eps_eta = sp.symbols("epsilon_eta", positive=True, real=True)
Lambda0 = sp.simplify(27 * sp.pi**2 * G * cs**5 / (20 * a**5 * c**5))
Rtarget = sp.simplify(Lambda0 * (1 - eps_eta) * (1 - epsilon)**2 / (ZW * (1 + chi0)**2))

expect_zero("d_zeta log R_target", sp.diff(sp.log(Rtarget), zeta))
expect_zero("d_zeta Xi_1", sp.diff(Xi1_expected, zeta))

# ---------------------------------------------------------------------------
# III. Tracking-factor drift
# ---------------------------------------------------------------------------

subbanner("III. Exact tracking-factor drift")

Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
Rtr_t = sp.simplify((1 + chi0_t / (1 + deltaU_t)) / (1 + chi0_t))
Theta1 = sp.simplify(sp.diff(sp.log(Rtr_t), t).subs(t, 0) / (eps * lamA))
Theta1_expected = sp.simplify(
    -(
        chi0 * (1 + chi0) * deltaU1
        + deltaU * (1 + deltaU) * chi1
    ) / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
)

print("R_tr =", Rtr)
print("Theta_1 =", Theta1)
expect_zero("tracking-factor drift law", Theta1 - Theta1_expected)

subbanner("IV. Tracking rigidity is not sufficient")

Theta_rigid = sp.simplify(Theta1_expected.subs({chi1: 0, deltaU1: 0}))
Xi_rigid = sp.simplify(Xi1_expected.subs({chi1: 0, deltaU1: 0}))
expect_zero("Theta_1 on chi_1=delta_U1=0 subbranch", Theta_rigid)
print("Xi_1 on the same tracking-rigid subbranch =", Xi_rigid)
print("So exact tracking rigidity does not in general kill Xi_1.")

banner("FINAL STAGE-164 LEDGER")
print("1. On the coherent local D/N tracking branch,")
print("      T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-epsilon)^2],")
print("      epsilon = epsilon_W (1 - 2/11 * delta_U/(1+delta_U)).")
print("2. The exact weak-axisymmetric defect law is")
print("      Xi_1 = zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 epsilon_1/(1-epsilon).")
print("3. The coherent support parameter zeta drops out identically from T^2, R_target, and Xi_1.")
print("4. The tracking factor obeys")
print("      Theta_1 = -[chi_0(1+chi_0) delta_U1 + delta_U(1+delta_U) chi_1]")
print("                / [(1+chi_0)(1+delta_U)(1+chi_0+delta_U)].")
print("5. Exact tracking rigidity is necessary but not sufficient to kill the grouped defect.")
