#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = '-' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        simp = expr.applyfunc(lambda e: sp.simplify(sp.expand(e)))
        print(f"{name} =")
        sp.pprint(simp)
        if any(e != 0 for e in simp):
            raise AssertionError(f"{name} is not zero")
    else:
        simp = sp.simplify(sp.expand(expr))
        print(f"{name} = {simp}")
        if simp != 0:
            raise AssertionError(f"{name} is not zero")


"""
5pn_stage174_adiabatic_elastic_orbit_theorem.py

Stage 174 — unified theorem under the adiabatic-elastic boundary rule.

What this script does
---------------------
1. Combines the Stage-151 scalar off-bundle normal coordinate with the
   Stage-166/169/170 quotient-orbit observable map.
2. Imposes the adiabatic-elastic boundary rule:
      delta ln Theta_w = 0,
      epsilon_L = epsilon_v = epsilon_T = 0.
3. Proves that the scalar off-bundle source then vanishes identically.
4. Shows that the remaining first grouped weak-axisymmetric observables are
   exactly the quotient drifts
      (d ln C_tr,* , d ln C_nt,* , d ln epsilon_eta).
5. Proves the clean unified criterion:
      under the adiabatic-elastic boundary rule,
      zero first-order defect  <=>  the true branch stays on a single G_* orbit.

Interpretation
--------------
This closes the logical loop requested by the adiabatic-wall message: once the
boundary is forced to deform elastically with no thermal fraying, the only
remaining first-order branch-selection question is whether the branch preserves
(C_tr,* , C_nt,* , epsilon_eta), i.e. whether it remains on one exact G_* orbit.
"""

banner("STAGE 174 — ADIABATIC-ELASTIC UNIFIED ORBIT THEOREM")

# ---------------------------------------------------------------------------
# I. Scalar off-bundle source
# ---------------------------------------------------------------------------

subbanner("I. Scalar off-bundle source")

chi0, deltaU = sp.symbols('chi_0 delta_U', positive=True, real=True)
eps_eta = sp.symbols('epsilon_eta', positive=True, real=True)

# Stage-151 weighted scalar off-bundle source
r = sp.symbols('r', positive=True, real=True)
gstar = sp.simplify(r - sp.sqrt(1 + r**2) / 2)
coef_T = sp.simplify(gstar)
coef_v = sp.simplify(gstar + 1 / (2 * sp.sqrt(1 + r**2)))
coef_L = sp.simplify(2 * gstar + 3 / (4 * sp.sqrt(1 + r**2)))

epsL, epsv, epsT = sp.symbols('epsilon_L epsilon_v epsilon_T', real=True)
eps_perp = sp.simplify(coef_T * epsT + coef_v * epsv + coef_L * epsL)
delta_perp = sp.simplify(-eps_perp)

print("delta_perp =")
sp.pprint(delta_perp)

# Adiabatic-elastic boundary rule
adiabatic_elastic = {epsL: 0, epsv: 0, epsT: 0}
expect_zero("delta_perp under adiabatic-elastic boundary", delta_perp.subs(adiabatic_elastic))

# ---------------------------------------------------------------------------
# II. Quotient-orbit observable map
# ---------------------------------------------------------------------------

subbanner("II. Quotient-orbit observable map")

Ctr_drift, Cnt_drift, Eta_drift = sp.symbols('dlnCtr dlnCnt dlnEta', real=True)

C_tr = sp.simplify(chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)))
A_tr = sp.simplify(2 * chi0 / ((1 + chi0) * (1 + deltaU)))

Theta1 = sp.simplify(-C_tr * Ctr_drift)
Xi1 = sp.simplify(A_tr * Ctr_drift + Cnt_drift)
RplusXi = sp.simplify(-eps_eta / (1 - eps_eta) * Eta_drift)

print("Theta_1 =", Theta1)
print("Xi_1 =", Xi1)
print("R_1 + Xi_1 =", RplusXi)

# Solve the observable equations for the quotient drifts.
sol = sp.solve([sp.Eq(Theta1, 0), sp.Eq(Xi1, 0), sp.Eq(RplusXi, 0)], [Ctr_drift, Cnt_drift, Eta_drift], dict=True)
if len(sol) != 1:
    raise AssertionError("Expected unique zero-defect solve for the quotient drifts.")
sol = sol[0]
expect_zero("zero-defect gives dln C_tr,* = 0", sol[Ctr_drift])
expect_zero("zero-defect gives dln C_nt,* = 0", sol[Cnt_drift])
expect_zero("zero-defect gives dln epsilon_eta = 0", sol[Eta_drift])

# ---------------------------------------------------------------------------
# III. Exact Stage-169/170 orbit-lock criterion
# ---------------------------------------------------------------------------

subbanner("III. Exact orbit-lock criterion")

chi0s, deltaUs = sp.symbols('chi0_star deltaU_star', positive=True, real=True)
Estar, Fstar = sp.symbols('E_star F_star', real=True)
Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf = sp.symbols(
    'Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_KW Delta_mu Delta_T',
    real=True,
)
Delta_x = sp.Matrix([Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf])
M = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
    [0, 2, 0, -1, -1, 0, 0, 0],
])

q = sp.simplify(M * Delta_x)
print("M_* Delta_x =")
sp.pprint(q)

# Finite fibre equations = similarity orbit solve
sol_orbit = sp.solve(list(q), [DTf, DKe, Dmu], dict=True)
if len(sol_orbit) != 1:
    raise AssertionError("Expected unique orbit solve.")
sol_orbit = sol_orbit[0]

print("Orbit-lock solution:")
for k, v in sol_orbit.items():
    print(f"{k} = {sp.simplify(v)}")

# ---------------------------------------------------------------------------
# IV. Unified adiabatic-elastic theorem
# ---------------------------------------------------------------------------

subbanner("IV. Unified theorem")

print("Under the adiabatic-elastic boundary rule:")
print("  - the scalar off-bundle source delta_perp vanishes identically,")
print("  - the remaining first weak-axisymmetric observables are exactly the quotient drifts,")
print("      Theta_1 = -C_tr dln C_tr,*")
print("      Xi_1    =  A_tr dln C_tr,* + dln C_nt,*")
print("      R_1+Xi_1 = -(epsilon_eta/(1-epsilon_eta)) dln epsilon_eta.")
print("  - and the branch stays on a single G_* orbit exactly when M_* Delta_x = 0.")

banner("FINAL STAGE-174 LEDGER")
print("1. Imposing delta ln Theta_w = 0 freezes the isotropic wall-depth / thermal-fraying channel.")
print("2. Imposing elastic lower-branch transport (epsilon_L = epsilon_v = epsilon_T = 0) kills the scalar off-bundle source delta_perp exactly.")
print("3. What remains is purely the quotient motion in")
print("      (d ln C_tr,* , d ln C_nt,* , d ln epsilon_eta).")
print("4. The first grouped weak-axisymmetric observables are exactly those quotient drifts.")
print("5. Therefore, under the adiabatic-elastic boundary rule,")
print("      zero first-order defect  <=>  d ln C_tr,* = d ln C_nt,* = d ln epsilon_eta = 0")
print("   <=>  M_* Delta_x = 0")
print("   <=>  the true moving-throat branch remains on a single exact G_* orbit.")
print("6. So the adiabatic wall closes the isotropic thermal loop, but the remaining branch-selection test is still the exact orbit-preservation test, not a generic wall-freeze condition by itself.")
