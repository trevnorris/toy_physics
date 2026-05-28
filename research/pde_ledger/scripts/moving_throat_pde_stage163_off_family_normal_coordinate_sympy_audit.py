#!/usr/bin/env python3
"""
moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.py

SymPy-backed audit for Stage 163.

Checks:
1. Exact parent compensation-defect scalar linearization.
2. Exact compensation-ratio transport through the off-family normal coordinate.
3. Exact microscopic parent-variable formula for delta_perp.
4. Exact outlet-defect transport formulas.
5. Tangent/normal split of the mouth-bias transport.

Provenance notes
----------------
- `gminus` is the same lower compensated branch called `g_lower` in Stage 162;
  the script renames it only to keep the calculus notation compact.
- The Family-1 numbers at the end are carried readbacks from the Stage 162/146
  benchmark branch and are used only as numeric coefficient summaries.
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


banner("STAGE 163 — OFF-FAMILY NORMAL COORDINATE")

# ------------------------------------------------------------------
# 1. Parent compensation defect and compensation-ratio transport
# ------------------------------------------------------------------
r, g = sp.symbols("r g", positive=True, real=True)
s = sp.sqrt(1 + r**2)
# Stage 162 lower compensated branch `g_lower`, renamed here to `gminus`.
gminus = r - s / 2

dg, dr = sp.symbols("dg dr", real=True)
gp = sp.diff(gminus, r)
F = 1 + r**2 - 4 * (g - r) ** 2
R = (g - r) ** 2 / (1 + r**2)

delta_perp = dg - gp * dr

dF = sp.diff(F, g).subs(g, gminus) * dg + sp.diff(F, r).subs(g, gminus) * dr
dR = sp.diff(R, g).subs(g, gminus) * dg + sp.diff(R, r).subs(g, gminus) * dr

expect_zero("delta F - 4 s delta_perp", dF - 4 * s * delta_perp)
expect_zero("delta R + delta_perp/s", dR + delta_perp / s)
expect_zero("even-preserving tangent motion keeps delta F = 0", dF.subs(dg, gp * dr))
expect_zero("even-preserving tangent motion keeps delta R = 0", dR.subs(dg, gp * dr))
print("g_-'(r) =", sp.simplify(gp))

# ------------------------------------------------------------------
# 2. Microscopic parent-variable formula for delta_perp
# ------------------------------------------------------------------
Ks, Kq, lam, gs, gq = sp.symbols("K_s K_q lam g_s g_q", positive=True, real=True)
dln_Ks, dln_Kq, dln_lam, dln_gs, dln_gq = sp.symbols(
    "dln_Ks dln_Kq dln_lam dln_gs dln_gq", real=True
)

# Linearized parent ratios on the lower compensated branch.
delta_r = r * (dln_lam - sp.Rational(1, 2) * dln_Ks - sp.Rational(1, 2) * dln_Kq)
delta_g = gminus * (dln_gq - dln_gs + sp.Rational(1, 2) * dln_Ks - sp.Rational(1, 2) * dln_Kq)
delta_perp_micro = sp.simplify(delta_g - gp * delta_r)

delta_perp_expected = sp.simplify(
    gminus * (dln_gq - dln_gs - dln_lam + dln_Ks)
    + (dln_Ks + dln_Kq - 2 * dln_lam) / (4 * s)
)
expect_zero("microscopic delta_perp identity", delta_perp_micro - delta_perp_expected)
delta_perp_nice = gminus*(dln_gq - dln_gs - dln_lam + dln_Ks) + (dln_Ks + dln_Kq - 2*dln_lam)/(4*s)
print("delta_perp microscopic form =", delta_perp_nice)

# ------------------------------------------------------------------
# 3. Outlet-defect transport
# ------------------------------------------------------------------
sigma_star, dkappaW, dgammaW = sp.symbols("sigma_star dkapW dgamW", real=True)
deltaC = 16 * sigma_star * delta_perp / s

# Stage 159 outlet transport packet carried into the off-family normal frame.
dE2 = (deltaC - 9 * sigma_star * dkappaW) / (27 * (1 - sigma_star))
dE4 = (5 * deltaC - 72 * sigma_star * dkappaW) / (243 * (1 - sigma_star))
DeltaQ = (deltaC - 27 * sigma_star * dgammaW) / (3 * (1 - sigma_star))

print("delta C =", sp.simplify(deltaC))
print("delta E2 =", sp.simplify(dE2))
print("delta E4 =", sp.simplify(dE4))
print("Delta_Q =", sp.simplify(DeltaQ))

expect_zero(
    "delta C - 4 sigma_star deltaF/(1+r^2)",
    deltaC - sigma_star * dF * 4 / (1 + r**2),
)

# ------------------------------------------------------------------
# 4. Tangent/normal split of the mouth-bias transport
# ------------------------------------------------------------------
Sigma0, Sstar, dSigma0, dS = sp.symbols("Sigma0 Sstar dSigma0 dS", real=True)
Rstar = sp.Rational(1, 4)
dR_from_perp = -delta_perp / s

dPi = sp.expand((1 - Rstar * Sstar) * dSigma0 - Sigma0 * (Rstar * dS + Sstar * dR_from_perp))
dPi_expected = sp.expand((1 - Sstar / 4) * dSigma0 - Sigma0 * dS / 4 + Sigma0 * Sstar * delta_perp / s)
expect_zero("delta Pi tangent/normal split", dPi - dPi_expected)
print("delta Pi =", sp.simplify(dPi_expected))

# ------------------------------------------------------------------
# 5. Numerical carry-forward values at the Family-1 point
# ------------------------------------------------------------------
rf1 = sp.Float("1.77799353547498", 30)
gf1 = sp.Float("0.758035078944663", 30)
Sigma0_can = sp.Float("4.651033550168876", 30)
S_can = sp.Float("0.6703621156734617", 30)
sf1 = sp.sqrt(1 + rf1**2)

coeff_F = sp.N(4 * sf1, 20)
coeff_R = sp.N(-1 / sf1, 20)
coeff_micro2 = sp.N(1 / (4 * sf1), 20)
coeff_Pi = sp.N(Sigma0_can * S_can / sf1, 20)
coeff_C = sp.N(16 / sf1, 20)

banner("Family-1 numerical coefficients")
print("4 sqrt(1+r_*^2) =", coeff_F)
print("-1/sqrt(1+r_*^2) =", coeff_R)
print("g_* =", gf1)
print("1/(4 sqrt(1+r_*^2)) =", coeff_micro2)
print("Sigma0_can S_can / sqrt(1+r_*^2) =", coeff_Pi)
print("16 / sqrt(1+r_*^2) =", coeff_C)

print("\nCarry-forward formulas:")
print("  delta_perp = delta g - g_-'(r_*) delta r")
print("  delta F = 4 sqrt(1+r_*^2) delta_perp")
print("  delta R_q = -delta_perp/sqrt(1+r_*^2)")
print("  delta_perp = g_* dln[(g_q K_s)/(g_s lam)] + [4 sqrt(1+r_*^2)]^{-1} dln[(K_s K_q)/lam^2]")
print("  delta C, delta E2, delta E4, Delta_Q are all linear in delta_perp")
