#!/usr/bin/env python3
"""
Stage 140 SymPy audit.

Checks:
1. Exact Family-1 compensation equivalence g = g_* <=> R = 1/4.
2. Exact self-matched traction law Sigma0 = 20 T_hat^2 / 9.
3. The carried Stage-139 renormalized canonical tuple satisfies the exact
   branch identities.
4. The renormalized point is genuinely above the original canonical point.
5. Tangent motion on the lower compensated family keeps delta R = 0.
6. Tangent motion kills delta C and forces delta kappa_W = 0 under
   canonical-even preservation.
7. The Stage-141 expansion point is the renormalized canonical tuple.
"""

from __future__ import annotations

import json
from pathlib import Path

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


def expect_close(name: str, value: sp.Expr, target: sp.Expr, tol: float = 1e-12) -> None:
    diff = abs(float(sp.N(value - target, 30)))
    print(f"{name} diff = {diff}")
    if diff > tol:
        raise AssertionError(f"{name} mismatch: diff={diff}")


def load_constants() -> dict:
    root = Path(__file__).resolve().parents[1]
    data = json.loads((root / "scripts/numerical/stage138_139_fixedpoint_samples.json").read_text())
    return data["constants"]


banner("STAGE 140 — CORE-MOUTH COEVOLUTION STATUS")

constants = load_constants()

banner("1. Exact Family-1 compensation identities")
r = sp.symbols("r", positive=True, real=True)
g = sp.symbols("g", real=True)
R = sp.simplify((g - r) ** 2 / (1 + r**2))
g_star = sp.simplify(r - sp.sqrt(1 + r**2) / 2)
expect_zero("R(g_*) - 1/4", R.subs(g, g_star) - sp.Rational(1, 4))

banner("2. Carry-forward numerical basepoint from Stages 138-139")
rF1_exact = sp.sqrt(4107 - 100 * sp.pi**2) / (10 * sp.pi)
g_star_exact = sp.simplify(rF1_exact - sp.sqrt(1 + rF1_exact**2) / 2)
rF1 = sp.Float(str(constants["rF1"]), 30)
g_star_num = sp.Float(str(constants["g_star"]), 30)
Sigma0_star = sp.Float(str(constants["Sigma0_star"]), 30)
T_hat_star = sp.Float(str(constants["T_hat_star"]), 30)
Pi_star = sp.Float(str(constants["Pi_star"]), 30)
Sigma0_can = sp.Float(str(constants["Sigma0_can_expected"]), 30)
S_can = sp.Float(str(constants["S_can_expected"]), 30)
Pi_can = sp.Float(str(constants["Pi_can_expected"]), 30)
T_hat_can = sp.Float(str(constants["T_hat_can_expected"]), 30)

expect_close("r_F1 radical check", sp.N(rF1_exact, 25), rF1)
expect_close("g_* lower-branch check", sp.N(g_star_exact, 25), g_star_num)
expect_close("self-matched traction law", Sigma0_can, sp.Rational(20, 9) * T_hat_can**2)

R_can = sp.N(((g_star_num - rF1) ** 2) / (1 + rF1**2), 25)
expect_close("R_can = 1/4", R_can, sp.Rational(1, 4))
expect_close("Pi_can identity", Pi_can, Sigma0_can * (1 - sp.Rational(1, 4) * S_can))

print("Sigma0_star =", Sigma0_star)
print("T_hat_star  =", T_hat_star)
print("Pi_star     =", Pi_star)
print("Sigma0_can  =", Sigma0_can)
print("T_hat_can   =", T_hat_can)
print("S_can       =", S_can)
print("Pi_can      =", Pi_can)

if not (float(Sigma0_can) > float(Sigma0_star) and float(T_hat_can) > float(T_hat_star) and float(Pi_can) > float(Pi_star)):
    raise AssertionError("Renormalized canonical tuple should sit strictly above the original canonical point.")

banner("3. Tangent-on-family and even-preservation handoff")
dg, dr = sp.symbols("delta_g delta_r", real=True)
s = sp.sqrt(1 + r**2)
gminus = sp.simplify(r - s / 2)
gp = sp.diff(gminus, r)
dR = sp.diff(R, g).subs(g, gminus) * dg + sp.diff(R, r).subs(g, gminus) * dr
expect_zero("tangent motion keeps delta R = 0", dR.subs(dg, gp * dr))

sigma_star = sp.symbols("sigma_star", real=True)
deltaC, dkappa = sp.symbols("deltaC delta_kappa", real=True)
dE2 = (deltaC - 9 * sigma_star * dkappa) / (27 * (1 - sigma_star))
dE4 = (5 * deltaC - 72 * sigma_star * dkappa) / (243 * (1 - sigma_star))
even_preservation = sp.solve([sp.Eq(dE2, 0), sp.Eq(dE4, 0)], [deltaC, dkappa], dict=True)
print("canonical-even preservation solutions =", even_preservation)
if even_preservation != [{deltaC: 0, dkappa: 0}]:
    raise AssertionError("Canonical-even preservation should pin deltaC = delta kappa_W = 0.")

deltaC_from_tangent = sp.simplify(-16 * sigma_star * dR.subs(dg, gp * dr))
expect_zero("tangent motion kills delta C", deltaC_from_tangent)

banner("4. Stage-141 expansion point")
Sigma0, dSigma0, Sstar, dS = sp.symbols("Sigma0 delta_Sigma0 Sstar delta_S", real=True)
Pi = (Sigma0 + dSigma0) * (1 - (sp.Rational(1, 4) + 0) * (Sstar + dS))
Pi_lin = sp.expand(Pi).subs({dSigma0 * dS: 0})
Pi0 = Sigma0 * (1 - Sstar / 4)
dPi_expected = sp.expand((1 - Sstar / 4) * dSigma0 - Sigma0 * dS / 4)
expect_zero("Stage-141 tangent expansion packet", (Pi_lin - Pi0) - dPi_expected)
print("dPi_tan at the renormalized canonical point =", dPi_expected.subs({Sigma0: Sigma0_can, Sstar: S_can}))

print("\nOpen branch note:")
print("  This capstone closes the reduced co-evolving classification and basepoint data.")
print("  It does not assert that the full moving-throat PDE has already realized the")
print("  balance-selected branch; that remains the open microscopic question carried forward.")
