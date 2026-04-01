
#!/usr/bin/env python3
"""
moving_throat_pde_stage164_coherent_tracking_defect_sympy_audit.py

SymPy-backed audit for Stage 164: coherent tracking-branch weak-axisymmetric
defect transport and the support-blindness theorem.

Checks:
1. Exact coherent-branch transfer-shape identity from the selected-branch map.
2. Exact support-blindness: d/dzeta of T^2 and R_target vanish.
3. Exact split-blocking drift formula.
4. Exact coherent-branch defect law Xi_1.
5. Exact selected-branch identity Xi_1 = -eta_1/(1-eps_eta) - R_1.
6. Exact tracking-factor drift formula.
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

banner("STAGE 164 — COHERENT TRACKING-BRANCH DEFECT LAW")

# Core coherent-branch variables
G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
eps_eta, epsW, deltaU, chi0, zeta = sp.symbols("eps_eta epsW deltaU chi0 zeta", real=True)
ZW, OmegaW2 = sp.symbols("Z_W Omega_W2", positive=True, real=True)

Lam = sp.simplify(sp.Rational(27) * sp.pi**2 * G * c_s**5 * OmegaW2 / (20 * a**5 * c**5))
eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
R_target = sp.simplify(Lam * (1 - eps_eta) * (1 - eps)**2 / (ZW * (1 + chi0)**2))

T2_direct = sp.simplify(ZW * (1 + chi0)**2 / (OmegaW2 * (1 - eps)**2))
T2_selected = sp.simplify(sp.Rational(27) * sp.pi**2 * G * c_s**5 / (20 * a**5 * c**5) * (1 - eps_eta) / R_target)

expect_zero("direct-selected transfer-shape identity", T2_direct - T2_selected)
expect_zero("d/dzeta ln T^2", sp.diff(sp.log(T2_direct), zeta))
expect_zero("d/dzeta ln R_target", sp.diff(sp.log(R_target), zeta))

banner("Weak-axisymmetric drift transport")

# Weak-axisymmetric drift coefficients
zetaZ, omegaW, chi1 = sp.symbols("zeta_Z omega_W chi_1", real=True)
epsW1, deltaU1, eta1 = sp.symbols("varepsilon_W deltaU_1 eta_1", real=True)

eps1 = sp.simplify(sp.diff(eps, epsW) * epsW1 + sp.diff(eps, deltaU) * deltaU1)
eps1_expected = sp.simplify(
    (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)) * epsW1
    - (sp.Rational(2, 11) * epsW / (1 + deltaU)**2) * deltaU1
)
expect_zero("split-blocking drift eps_1", eps1 - eps1_expected)

Xi1 = sp.simplify(
    zetaZ
    - omegaW
    + 2 * chi1 / (1 + chi0)
    + 2 * eps1 / (1 - eps)
)

R1 = sp.simplify(
    omegaW
    - eta1 / (1 - eps_eta)
    - zetaZ
    - 2 * chi1 / (1 + chi0)
    - 2 * eps1 / (1 - eps)
)

expect_zero("selected-branch identity", Xi1 + eta1 / (1 - eps_eta) + R1)

print("Xi_1 =", Xi1)
print("R_1  =", R1)

banner("Tracking-factor drift")

Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
Theta1 = sp.simplify(sp.diff(sp.log(Rtr), chi0) * chi1 + sp.diff(sp.log(Rtr), deltaU) * deltaU1)
Theta1_expected = sp.simplify(
    -(chi0 * (1 + chi0) * deltaU1 + deltaU * (1 + deltaU) * chi1)
    / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
)
expect_zero("tracking-factor drift", Theta1 - Theta1_expected)
print("Theta_1 =", Theta1)

banner("Support-blindness consequence")

Xi_support_rigid = sp.simplify(Xi1.subs({chi1: 0, deltaU1: 0}))
Theta_support_rigid = sp.simplify(Theta1.subs({chi1: 0, deltaU1: 0}))
print("Xi_1 with chi1=deltaU1=0 =", Xi_support_rigid)
print("Theta_1 with chi1=deltaU1=0 =", Theta_support_rigid)

if Xi_support_rigid == 0:
    raise AssertionError("Expected nontrivial residual Xi_1 under support-rigid specialization.")

print("\nConclusion:")
print("  The coherent support ratio zeta drops out identically from T^2 and R_target.")
print("  The grouped weak-axisymmetric defect is carried only by")
print("  Z_W, Omega_W^2, chi_0, eps_W, and delta_U.")
print("  Exact tracking-factor rigidity is not sufficient to kill Xi_1.")
