
#!/usr/bin/env python3
"""
moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py

SymPy-backed audit for Stage 181: coherent tracking-branch weak-axisymmetric
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

banner("STAGE 181 — COHERENT TRACKING-BRANCH DEFECT LAW")

# Core coherent-branch variables
G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
eps_eta, epsW, deltaU, chi0, zeta = sp.symbols("eps_eta epsW deltaU chi0 zeta", real=True)
ZW, OmegaW2 = sp.symbols("Z_W Omega_W2", positive=True, real=True)

Lam = sp.simplify(sp.Rational(27) * sp.pi**2 * G * c_s**5 * OmegaW2 / (20 * a**5 * c**5))
eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
R_target = sp.simplify(Lam * (1 - eps_eta) * (1 - eps)**2 / (ZW * (1 + chi0)**2))

T2_direct = sp.simplify(ZW * (1 + chi0)**2 / (OmegaW2 * (1 - eps)**2))
T2_selected = sp.simplify(sp.Rational(27) * sp.pi**2 * G * c_s**5 / (20 * a**5 * c**5) * (1 - eps_eta) / R_target)
Mmix = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps)))
Msupp = sp.simplify(8 * zeta * ZW * (1 + chi0) ** 2 / (sp.pi**2 * (1 - eps_eta) * (1 - zeta * eps)))
Ssupport = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
Mtr = sp.simplify(Mmix + Msupp)
product_loaded = sp.simplify(8 * Lam * (1 - eps) / sp.pi**2 * Ssupport)
R_target_loaded = sp.simplify(product_loaded / Mtr)
T2_loaded = sp.simplify((Lam / OmegaW2) * (1 - eps_eta) / R_target_loaded)

expect_zero("direct-selected transfer-shape identity", T2_direct - T2_selected)
expect_zero("support-loaded R_target reconstruction", R_target_loaded - R_target)
expect_zero("support-loaded T^2 reconstruction", T2_loaded - T2_direct)
expect_zero("d/dzeta ln T^2 (support-loaded route)", sp.diff(sp.log(T2_loaded), zeta))
expect_zero("d/dzeta ln R_target (support-loaded route)", sp.diff(sp.log(R_target_loaded), zeta))

bad = sp.symbols("bad", nonzero=True, real=True)
Msupp_spoiled = sp.simplify(Msupp + bad * zeta * Mmix)
R_target_spoiled = sp.simplify(product_loaded / (Mmix + Msupp_spoiled))
dlnR_spoiled = sp.simplify(sp.diff(sp.log(R_target_spoiled), zeta).subs(bad, 1))
print("spoiled d/dzeta ln R_target =", dlnR_spoiled)
if dlnR_spoiled == 0:
    raise AssertionError("Expected a spoiled support packet to break R_target blindness.")

banner("Weak-axisymmetric drift transport")

# Weak-axisymmetric drift coefficients
zetaZ, omegaW, chi1 = sp.symbols("zeta_Z omega_W chi_1", real=True)
epsW1, deltaU1, eta1 = sp.symbols("varepsilon_W deltaU_1 eta_1", real=True)
s = sp.symbols("s", real=True)

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

Lam_pert = Lam.subs(OmegaW2, OmegaW2*(1 + s*omegaW))
eps_pert = (epsW + s*epsW1) * (1 - sp.Rational(2, 11) * (deltaU + s*deltaU1) / (1 + deltaU + s*deltaU1))
R_target_pert = Lam_pert * (1 - (eps_eta + s*eta1)) * (1 - eps_pert)**2 \
                / ((ZW*(1 + s*zetaZ)) * (1 + (chi0 + s*chi1))**2)
R1_derived = sp.simplify(sp.diff(sp.log(R_target_pert), s).subs(s, 0))
expect_zero("R_1 derived from R_target matches closed form", R1_derived - R1)

expect_zero("selected-branch identity", Xi1 + eta1 / (1 - eps_eta) + R1)

print("Xi_1 =", Xi1)
print("R_1  =", R1)

eps_pert = (epsW + s*epsW1) * (1 - sp.Rational(2, 11) * (deltaU + s*deltaU1) / (1 + deltaU + s*deltaU1))
T2_direct_pert = (ZW*(1 + s*zetaZ)) * (1 + (chi0 + s*chi1))**2 \
                 / ((OmegaW2*(1 + s*omegaW)) * (1 - eps_pert)**2)
Xi1_derived = sp.simplify(sp.diff(sp.log(T2_direct_pert), s).subs(s, 0))
expect_zero("Xi_1 derived from T^2 matches defect law", Xi1_derived - Xi1)

# Support-blindness of the defect itself: derive Xi_1 through the zeta-bearing
# support-loaded shape (line 56 proves T2_loaded == T2_direct) and confirm no zeta.
T2_loaded_pert = T2_loaded.subs({
    ZW: ZW*(1 + s*zetaZ), OmegaW2: OmegaW2*(1 + s*omegaW),
    chi0: chi0 + s*chi1, epsW: epsW + s*epsW1, deltaU: deltaU + s*deltaU1,
})
Xi1_loaded = sp.simplify(sp.diff(sp.log(T2_loaded_pert), s).subs(s, 0))
expect_zero("d/dzeta Xi_1 (support-loaded route)", sp.diff(Xi1_loaded, zeta))

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
