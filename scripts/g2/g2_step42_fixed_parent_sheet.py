#!/usr/bin/env python3
"""
g2_step42_fixed_parent_sheet.py

Step 42 of the g-2 rebuild.

What this script does
---------------------
1. Starts from the Step-41 anomaly-only adiabatic branch:
      d ln K_s = 0,
      d ln P_0 = ell,
   so the microscopic drift vector is the pure core/outgoing retuning already
   derived there.
2. Pushes that drift vector back into the parent overlap ratios
      r := lambda / sqrt(K_s K_q),
      g := g_q sqrt(K_s) / (g_s sqrt(K_q)),
      r_c := lambda^2 / (K_s K_q),
   and proves that all three remain exactly fixed.
3. Shows that the auxiliary D/N geometry L_W and the even-preserving outlet data
   stay frozen.
4. Computes the induced drift of the balanced-core loading share
      rho_c = g_s^2 / K_s,
      sigma_c = rho_c / 4,
   and derives the exact odd-outlet law gamma_c(ell) required by the electron
   target chi_Q = exp(-ell).
5. Rewrites the same result in the compensated Robin-mixed outlet variables
      sigma_W, rho_R, kappa_W, gamma_W.

Interpretation
--------------
On the anomaly-only adiabatic branch, the quartic electron sliver does not come
from motion in the parent compensation ratios. It keeps the system on the same
lower compensated family and the same even D/N geometry, while forcing only a
small odd-outlet detuning on top of a slightly reduced loading share.
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


banner("STEP 42A — ANOMALY-ONLY ADIABATIC DRIFT VECTOR")

ell = sp.symbols("ell", real=True)

# Step-41 anomaly-only branch: the thermodynamic ground-state squish is frozen.
dln_Ks = sp.Integer(0)
dln_P0 = ell

dln_a = sp.Integer(0)
dln_cs = ell / 5
dln_Kq = 2 * ell / 5
dln_vw0 = ell / 5
dln_Tm = -ell / 5
dln_gs = -ell / 5
dln_gq = sp.Integer(0)
dln_lambda = ell / 5

drift_vector = {
    "d ln a": dln_a,
    "d ln c_s": dln_cs,
    "d ln K_q": dln_Kq,
    "d ln v_w0": dln_vw0,
    "d ln T_m": dln_Tm,
    "d ln g_s": dln_gs,
    "d ln g_q": dln_gq,
    "d ln lambda": dln_lambda,
}
for k, v in drift_vector.items():
    print(f"{k:12s} = {sp.simplify(v)}")

print("\nReading:")
print("  - the wall squish lane is frozen,")
print("  - the mouth radius is frozen,")
print("  - and only the core/outgoing block retunes.")


banner("STEP 42B — THE PARENT COMPENSATION SHEET IS FIXED")

rc = sp.symbols("r_c", positive=True, real=True)

# Parent overlap ratios from the lower compensated branch.
dln_r = sp.simplify(dln_lambda - (dln_Ks + dln_Kq) / 2)
dln_g = sp.simplify(dln_gq + dln_Ks / 2 - dln_gs - dln_Kq / 2)
dln_rc = sp.simplify(2 * dln_lambda - dln_Ks - dln_Kq)

# Auxiliary D/N tube law: L_W = (pi a / 2) * sqrt((1 + r_c)/3)
dln_LW = sp.simplify(dln_a + dln_rc / (2 * (1 + rc)))

print("d ln r   =", dln_r)
print("d ln g   =", dln_g)
print("d ln r_c =", dln_rc)
print("d ln L_W =", dln_LW)

expect_zero("fixed parent ratio r", dln_r)
expect_zero("fixed parent ratio g", dln_g)
expect_zero("fixed parent ratio r_c", dln_rc)
expect_zero("fixed auxiliary D/N length", dln_LW)

print("\nInterpretation:")
print("  The anomaly-only adiabatic increment stays exactly on the same lower")
print("  compensated parent sheet. It does not move the branch in g, r, or r_c,")
print("  and it does not change the side-tube D/N geometry L_W.")


banner("STEP 42C — THE LOADING SHARE DRIFTS, THE EVEN OUTLET DOES NOT")

rho_c_star, sigma_c_star = sp.symbols("rho_c_star sigma_c_star", positive=True, real=True)

# Balanced-core loading share on the fixed parent sheet.
dln_rho_c = sp.simplify(2 * dln_gs - dln_Ks)
dln_sigma_c = sp.simplify(dln_rho_c)  # because sigma_c = rho_c / 4 on the balanced family

rho_c = sp.simplify(rho_c_star * sp.exp(dln_rho_c))
sigma_c = sp.simplify(sigma_c_star * sp.exp(dln_sigma_c))

print("d ln rho_c   =", dln_rho_c)
print("d ln sigma_c =", dln_sigma_c)
print("rho_c(ell)   =", rho_c)
print("sigma_c(ell) =", sigma_c)

expect_zero("balanced-family loading share law", dln_sigma_c - dln_rho_c)

# The even-preserving hybrid / core outlet keeps kappa fixed.
kappa_c = sp.Rational(1, 3)
print("kappa_c =", kappa_c)

print("\nInterpretation:")
print("  The anomaly does not reopen the even outlet geometry. It leaves the")
print("  even D/N side-channel fixed, but it lowers the balanced loading share")
print("  by exp(-2 ell / 5) because g_s drifts while K_s stays fixed.")


banner("STEP 42D — EXACT ODD DETUNING LAW ON THE FIXED PARENT SHEET")

# Exact electron target from the outgoing-normalization branch.
x = sp.symbols("x", positive=True, real=True)
chi_Q = sp.simplify(sp.exp(-ell))
chi_from_x = sp.simplify(1 / (1 + x))

# Hybrid / balanced-core normalization law.
gamma_c = sp.simplify((1 - chi_Q * (1 - sigma_c)) / (9 * sigma_c))
gamma_c_x = sp.simplify(gamma_c.subs(sp.exp(ell), 1 + x))
gamma_c_compact = sp.simplify((sigma_c + x) / (9 * sigma_c * (1 + x)))

print("chi_Q(ell) =", chi_Q)
print("chi_Q(x)   =", chi_from_x)
print("gamma_c(ell) =", gamma_c)
print("gamma_c(x)   =", gamma_c_x)
print("gamma_c compact =", gamma_c_compact)

expect_zero("x = exp(ell) - 1 rewrite", sp.simplify(chi_Q - chi_from_x.subs(x, sp.exp(ell) - 1)))
expect_zero("compact odd-detuning law", sp.simplify(gamma_c_x.subs(x, sp.exp(ell) - 1) - gamma_c))
expect_zero("gamma_c compact form", sp.simplify(gamma_c_compact.subs(x, sp.exp(ell) - 1) - gamma_c))

# Bare odd coefficient on the same fixed parent sheet.
gamma0 = sp.simplify((1 + rc) * gamma_c)
print("gamma_0(ell) =", gamma0)

# First-order tangent around ell = 0.
gamma_c_series = sp.series(gamma_c, ell, 0, 2).removeO()
gamma0_series = sp.series(gamma0, ell, 0, 2).removeO()
print("gamma_c series =", gamma_c_series)
print("gamma_0 series =", gamma0_series)

print("\nReading:")
print("  On the fixed parent compensation sheet, the anomaly target chi_Q = exp(-ell)")
print("  is realized entirely by an odd-outlet detuning gamma_c (or gamma_0) on top")
print("  of the softened loading share sigma_c(ell).")


banner("STEP 42E — COMPENSATED ROBIN-MIXED OUTLET VARIABLES")

# Step-32/33 hybrid notation.
sigma_W = sp.simplify(sigma_c)
rho_R = sp.simplify(4 * sigma_W)
kappa_W = sp.Rational(1, 3)
gamma_W = sp.simplify((sigma_W + x) / (9 * sigma_W * (1 + x)))

print("sigma_W(ell) =", sigma_W)
print("rho_R(ell)   =", rho_R)
print("kappa_W      =", kappa_W)
print("gamma_W      =", gamma_W)

rho_R_series = sp.series(rho_R, ell, 0, 2).removeO()
gamma_W_series = sp.series(gamma_W.subs(x, sp.exp(ell) - 1), ell, 0, 2).removeO()
print("rho_R series   =", rho_R_series)
print("gamma_W series =", gamma_W_series)

print("\nFinal reading of Step 42:")
print("  - the anomaly-only adiabatic branch keeps the parent compensation ratios fixed,")
print("  - keeps the even side-channel geometry fixed,")
print("  - lowers the balanced loading share by exp(-2 ell / 5),")
print("  - and forces only a small odd-outlet detuning to hit the electron target.")
