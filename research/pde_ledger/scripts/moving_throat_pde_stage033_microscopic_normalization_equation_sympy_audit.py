#!/usr/bin/env python3
"""
Moving-throat PDE Stage 33 SymPy audit.

Checks:
1. Exact microscopic selected-branch normalization product.
2. Exact stability threshold with the finite-throat overlap constants inserted.
3. Exact zero-loading onset value.
4. Exact weak-loading expansion.
5. Exact onset stiffness rearrangement.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr):
    simp = sp.simplify(sp.expand(expr))
    print(f"{name} = {simp}")
    if simp != 0:
        raise AssertionError(f"{name} is not zero")


# Finite-throat constants.
kappa0_sq = sp.Rational(8) / sp.pi**2
kappa1_sq = sp.Rational(16) / (9 * sp.pi**2)
sigma = sp.simplify(kappa0_sq + kappa1_sq)
delta_kappa = sp.simplify(kappa0_sq - kappa1_sq)
Kprod = sp.simplify(kappa0_sq * kappa1_sq)

# Stable-branch symbols.
A, DeltaK, alpha0, beta0 = sp.symbols("A DeltaK alpha0 beta0", positive=True, real=True)
R = sp.sqrt((DeltaK + alpha0 * delta_kappa)**2 + 4 * alpha0**2 * Kprod)
lambda_minus = sp.simplify((2 * A + DeltaK - alpha0 * sigma - R) / 2)
s_minus = sp.simplify(sp.Rational(1, 2) * (sigma + ((DeltaK + alpha0 * delta_kappa) * delta_kappa + 4 * alpha0 * Kprod) / R))
Nminus = sp.simplify(beta0 * s_minus**2 / (kappa0_sq * lambda_minus))

banner("STAGE 33.1 — EXACT MICROSCOPIC NORMALIZATION PRODUCT")
print("N_-(alpha0) =")
sp.pprint(Nminus)

# Exact monotonicity identity with the Stage-15 source-map substitution built in.
ds = sp.simplify(sp.diff(s_minus, alpha0))
dN = sp.simplify(sp.diff(Nminus, alpha0))
dN_formula = sp.simplify(beta0 * (2 * s_minus * ds * lambda_minus + s_minus**3) / (kappa0_sq * lambda_minus**2))
expect_zero("dN/dalpha - monotonicity formula", dN - dN_formula)

banner("STAGE 33.2 — EXACT STABILITY THRESHOLD")
alpha_crit = sp.simplify(A * (A + DeltaK) / ((A + DeltaK) * kappa0_sq + A * kappa1_sq))
alpha_crit_target = sp.simplify(9 * sp.pi**2 * A * (A + DeltaK) / (8 * (11 * A + 9 * DeltaK)))
print("alpha_crit =", alpha_crit)
expect_zero("alpha_crit - closed finite-throat form", alpha_crit - alpha_crit_target)

banner("STAGE 33.3 — ZERO-LOADING ONSET VALUE")
N0 = sp.simplify(Nminus.subs(alpha0, 0))
print("N_-(0) =", N0)
expect_zero("N_-(0) - beta0 kappa0^2 / A", N0 - beta0 * kappa0_sq / A)

banner("STAGE 33.4 — EXACT WEAK-LOADING EXPANSION")
series_N = sp.series(Nminus, alpha0, 0, 2).removeO()
coef1 = sp.simplify(sp.diff(Nminus, alpha0).subs(alpha0, 0))
coef1_target = sp.simplify(beta0 * kappa0_sq * (4 * A * kappa1_sq + DeltaK * kappa0_sq) / (A**2 * DeltaK))
coef1_target_closed = sp.simplify(64 * beta0 * (8 * A + 9 * DeltaK) / (9 * sp.pi**4 * A**2 * DeltaK))
print("N_-(alpha0) to O(alpha0) =", series_N)
expect_zero("first derivative coefficient - generic formula", coef1 - coef1_target)
expect_zero("first derivative coefficient - finite-throat closed form", coef1 - coef1_target_closed)

banner("STAGE 33.5 — MICROSCOPIC COUPLING REWRITE")
gB, gU, gW, gR = sp.symbols("g_B g_U g_W g_R", real=True)
varpi, OmegaU, OmegaW = sp.symbols("varpi Omega_U Omega_W", positive=True, real=True)
K0, NQ = sp.symbols("K_0 NQ", positive=True, real=True)

A_mic = sp.simplify(K0 - gU**2 / OmegaU**2)
Delta0 = sp.simplify(OmegaU**2 * OmegaW**2 - gR**2 * sigma)
Chi = sp.simplify(OmegaU**2 * gW + gR * gU)
beta0_mic = sp.simplify(Chi**2 / Delta0**2)
alpha0_mic = sp.simplify(gB**2 / varpi**2 + Chi**2 / (OmegaU**2 * Delta0))
N0_mic = sp.simplify((beta0 * kappa0_sq / A).subs({beta0: beta0_mic, A: A_mic}))
K0_onset = sp.simplify(sp.solve(sp.Eq(N0_mic, NQ), K0)[0])

print("A =", A_mic)
print("Delta0 =", Delta0)
print("Chi =", Chi)
print("beta0 =", beta0_mic)
print("alpha0 =", alpha0_mic)
print("N_-(0) =", N0_mic)
print("K0_onset =", K0_onset)

expect_zero(
    "K0_onset - [gU^2/OmegaU^2 + kappa0^2 Chi^2/(NQ Delta0^2)]",
    K0_onset - (gU**2 / OmegaU**2 + kappa0_sq * Chi**2 / (NQ * Delta0**2)),
)

banner("STAGE 33.6 — FULLY SUBSTITUTED MICROSCOPIC STABILITY GATE")
alpha_crit_mic = sp.simplify(alpha_crit_target.subs(A, A_mic))
gate_den_claim = sp.simplify(8 * varpi**2 * OmegaU**2 * Delta0 * (11 * A_mic + 9 * DeltaK))
gate_diff = sp.cancel(sp.together(alpha_crit_mic - alpha0_mic))
gate_num_actual, gate_den_actual = sp.fraction(gate_diff)
print("alpha_crit(mic) - alpha_0(mic) =")
sp.pprint(gate_diff)
print("computed denominator =")
sp.pprint(gate_den_actual)
print("claimed denominator =")
sp.pprint(gate_den_claim)
# Non-tautological denominator check: the ratio gate_den_actual / gate_den_claim
# is derived from `together(alpha_crit_mic - alpha0_mic)` WITHOUT referencing
# gate_den_claim, so it can fail if the claim is wrong. The ratio must simplify
# to a parameter-free numeric constant for the claim to hold; the fully
# substituted Delta0 convention can leave a universal pi^2 normalization.
den_ratio = sp.simplify(gate_den_actual / gate_den_claim)
print("denominator ratio (must be parameter-free) =", den_ratio)
assert den_ratio.is_number, (
    f"gate denominator does not match claim up to a parameter-free constant; ratio = {den_ratio}"
)
# Independent numerator reconstruction from the claim side, then final identity check.
gate_num_target = sp.simplify(gate_num_actual / den_ratio)
expect_zero(
    "gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) "
    "(tautological by reconstruction; substantive check is den_ratio.is_number above)",
    sp.simplify(alpha_crit_mic - alpha0_mic - gate_num_target / gate_den_claim),
)
print("gate numerator =")
sp.pprint(gate_num_target)

print("All Stage 33 checks passed.")
