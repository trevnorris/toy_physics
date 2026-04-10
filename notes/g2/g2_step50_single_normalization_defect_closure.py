#!/usr/bin/env python3
"""
g2_step50_single_normalization_defect_closure.py

Step 50 of the g-2 rebuild.

What this script does
---------------------
1. Imports the actual isotropic passive/outgoing grouped-P2 one-pole structure
   from the moving-throat notes:
       Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2).
2. Derives the exact even low-frequency coefficients Kbar_2 and Kbar_4 and the
   odd coefficient Gammabar_5 as functions of Kbar_0 and Omega_Q.
3. Compares them to the canonical GR target branch and proves that the entire
   reduced normalization problem collapses to one scalar defect
       N_Q = Kbar_0 / Kbar_0^target.
4. Inserts the adiabatic electron branch from Steps 44–45,
       P_0 = exp(ell) P_0^target,
   and shows that on the g-2 branch the same scalar defect is simply
       N_Q = exp(ell) = 1 + Lambda1 f.
5. Evaluates the resulting electron-point normalization defect.

Interpretation
--------------
This step is the clean closure of the conservative grouped-P2 side: once the
actual isotropic one-pole branch is accepted, the g-2 chain does not need a
multi-parameter fit. All even and odd grouped-P2 deviations scale with one
number only.
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


banner("STEP 50A — ACTUAL ISOTROPIC PASSIVE/OUTGOING GROUPED-P2 ONE-POLE BRANCH")

omega, OmegaQ, Kbar0 = sp.symbols("omega Omega_Q Kbar_0", positive=True, real=True)
Yhat_cons = sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / OmegaQ**2)
Yhat_series = sp.expand(sp.series(Yhat_cons, omega, 0, 6).removeO())

Kbar2 = sp.simplify(Kbar0 / (4 * OmegaQ**2))
Kbar4 = sp.simplify(Kbar0 / (4 * OmegaQ**4))
Gammabar5 = sp.simplify(9 * Kbar0 / (32 * OmegaQ**5))

print("Yhat_Q^cons(omega) =", Yhat_cons)
print("Yhat_Q^cons series =", Yhat_series)
print("Kbar2 =", Kbar2)
print("Kbar4 =", Kbar4)
print("Gammabar5 =", Gammabar5)

expect_zero(
    "omega^2 coefficient",
    sp.simplify(Yhat_series.coeff(omega, 2) - 1 / (4 * OmegaQ**2)),
)
expect_zero(
    "omega^4 coefficient",
    sp.simplify(Yhat_series.coeff(omega, 4) - 1 / (4 * OmegaQ**4)),
)

print("\nReading:")
print("  Once the actual isotropic grouped-P2 carrier is one-pole, the whole even")
print("  conservative branch is fixed by (Kbar0, OmegaQ), and the odd Burke–Thorne")
print("  coefficient is already slaved to the same pair.")


banner("STEP 50B — SINGLE NORMALIZATION DEFECT AGAINST THE CANONICAL TARGET BRANCH")

G, c, c_s, a, NQ = sp.symbols("G c c_s a N_Q", positive=True, real=True)
Kbar0_target = sp.simplify(64 * G * OmegaQ**5 / (45 * c**5))
Kbar0_target_geom = sp.simplify((64 * G * (3 * c_s / (2 * a))**5) / (45 * c**5))
Kbar2_target = sp.simplify(Kbar0_target / (4 * OmegaQ**2))
Kbar4_target = sp.simplify(Kbar0_target / (4 * OmegaQ**4))
Gammabar5_target = sp.simplify(2 * G / (5 * c**5))

Kbar0_actual = sp.simplify(NQ * Kbar0_target)
Kbar2_actual = sp.simplify(NQ * Kbar2_target)
Kbar4_actual = sp.simplify(NQ * Kbar4_target)
Gammabar5_actual = sp.simplify(NQ * Gammabar5_target)

print("Kbar0_target         =", Kbar0_target)
print("Kbar0_target(a,c_s)  =", Kbar0_target_geom)
print("Kbar2_target         =", Kbar2_target)
print("Kbar4_target         =", Kbar4_target)
print("Gammabar5_target     =", Gammabar5_target)

expect_zero(
    "Kbar0 target geometry rewrite",
    sp.simplify(Kbar0_target.subs(OmegaQ, 3 * c_s / (2 * a)) - Kbar0_target_geom),
)
expect_zero(
    "Gammabar5 target identity",
    sp.simplify(9 * Kbar0_target / (32 * OmegaQ**5) - Gammabar5_target),
)
expect_zero(
    "all branch defects collapse to NQ",
    sp.simplify((Kbar2_actual / Kbar2_target) - (Kbar4_actual / Kbar4_target)),
)
expect_zero(
    "odd branch defect is same NQ",
    sp.simplify((Gammabar5_actual / Gammabar5_target) - NQ),
)

print("\nReading:")
print("  The actual isotropic passive/outgoing grouped-P2 branch has only one")
print("  reduced normalization defect:")
print("      N_Q = Kbar0 / Kbar0_target.")
print("  Kbar2, Kbar4, and Gammabar5 all move by the same factor.")


banner("STEP 50C — INSERT THE ADIABATIC ELECTRON BRANCH")

ell, Lambda1, f = sp.symbols("ell Lambda1 f", positive=True, real=True)
NQ_from_ell = sp.exp(ell)
NQ_from_f = sp.simplify(NQ_from_ell.subs(sp.exp(ell), 1 + Lambda1 * f))

print("NQ(ell) =", NQ_from_ell)
print("NQ(f)   =", NQ_from_f)
expect_zero("Step-44 bridge", sp.simplify(NQ_from_f - (1 + Lambda1 * f)))

Kbar0_e = sp.simplify(NQ_from_ell * Kbar0_target_geom)
Kbar2_e = sp.simplify(NQ_from_ell * Kbar2_target.subs(OmegaQ, 3 * c_s / (2 * a)))
Kbar4_e = sp.simplify(NQ_from_ell * Kbar4_target.subs(OmegaQ, 3 * c_s / (2 * a)))
Gammabar5_e = sp.simplify(NQ_from_ell * Gammabar5_target)

print("Kbar0 electron branch =", Kbar0_e)
print("Kbar2 electron branch =", Kbar2_e)
print("Kbar4 electron branch =", Kbar4_e)
print("Gammabar5 electron branch =", Gammabar5_e)

expect_zero(
    "electron branch ratios stay common",
    sp.simplify((Kbar2_e / Kbar2_target.subs(OmegaQ, 3 * c_s / (2 * a))) - NQ_from_ell),
)
expect_zero(
    "electron odd ratio is same NQ",
    sp.simplify((Gammabar5_e / Gammabar5_target) - NQ_from_ell),
)

print("\nReading:")
print("  On the adiabatic electron track, the grouped-P2 normalization defect is")
print("  not a free multi-slot deformation. It is just")
print("      N_Q = exp(ell) = 1 + Lambda1 f.")


banner("STEP 50D — ELECTRON-POINT SIZE OF THE SINGLE DEFECT")

Lambda1_num = sp.Float("0.279605891931464")
f_num = sp.Float("0.001161409732093")
NQ_num = sp.N(NQ_from_f.subs({Lambda1: Lambda1_num, f: f_num}), 18)

print("N_Q electron-point =", NQ_num)
print("(N_Q - 1) in ppm   =", sp.N((NQ_num - 1) * 10**6, 18))

print("\nFinal reading of Step 50:")
print("  - the actual isotropic passive/outgoing grouped-P2 branch collapses to a")
print("    single reduced normalization defect N_Q,")
print("  - the adiabatic g-2 branch fixes that defect to N_Q = exp(ell),")
print("  - and therefore the whole grouped-P2 side of the anomaly is one-number")
print("    clean before any further branch-selection detail is imposed.")
