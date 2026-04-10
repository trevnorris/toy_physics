#!/usr/bin/env python3
"""
5pn_stage94_mixed_sidechannel_pole.py

Stage 94 audit: explicit mixed A_w/F_{mu w}-type side-channel pole.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 94 — EXPLICIT MIXED SIDE-CHANNEL POLE")

sigma_W, kappa_W, gamma_W, z = sp.symbols("sigma_W kappa_W gamma_W z", real=True)

Lambda_out = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9
hidden_inv = sp.series(1 / (1 - kappa_W * z**2 - sp.I * gamma_W * z**5), z, 0, 6).removeO()
print("1/(1-kappa_W z^2 - i gamma_W z^5) =")
sp.pprint(hidden_inv)

Lambda_mix = sp.expand(Lambda_out - sigma_W * hidden_inv)
print("Lambda_2^mix(z) =")
sp.pprint(Lambda_mix)

L0 = sp.expand(Lambda_mix).coeff(z, 0)
L2 = sp.expand(Lambda_mix).coeff(z, 2)
L4 = sp.expand(Lambda_mix).coeff(z, 4)
L5 = sp.simplify(sp.expand(Lambda_mix).coeff(z, 5) / sp.I)

expect_zero("L0 + (3 + sigma_W)", L0 + 3 + sigma_W)
expect_zero("L2 - (1/3 - sigma_W kappa_W)", L2 - (sp.Rational(1, 3) - sigma_W * kappa_W))
expect_zero("L4 - (1/9 - sigma_W kappa_W^2)", L4 - (sp.Rational(1, 9) - sigma_W * kappa_W**2))
expect_zero("L5 - (1/9 - sigma_W gamma_W)", L5 - (sp.Rational(1, 9) - sigma_W * gamma_W))

eq2 = sp.simplify(-L2 / L0 - sp.Rational(1, 9))
eq4 = sp.simplify(L2**2 / L0**2 - L4 / L0 - sp.Rational(4, 81))
print("canonical-even z^2 condition =", sp.factor(eq2))
print("canonical-even z^4 condition =", sp.factor(eq4))

# The z^2 condition factors to sigma_W*(9 kappa_W + 1)=0.
expect_zero("z^2 condition factor", sp.factor(sp.together(eq2)) + sigma_W * (9 * kappa_W + 1) / (9 * (sigma_W + 3)))

eq4_after_kappa = sp.simplify(eq4.subs(kappa_W, -sp.Rational(1, 9)))
print("z^4 condition after kappa_W=-1/9 =", sp.factor(eq4_after_kappa))
expect_zero("z^4 condition then forces sigma_W=0", sp.factor(sp.together(eq4_after_kappa)) + 4 * sigma_W / (81 * (sigma_W + 3)))

chi_mix = sp.simplify((-L5 / L0) / sp.Rational(1, 27))
print("chi_Q^mix =", sp.factor(chi_mix))
expect_zero("chi_Q^mix formula", chi_mix - 3 * (1 - 9 * sigma_W * gamma_W) / (3 + sigma_W))

chi_mix_series = sp.series(chi_mix, sigma_W, 0, 2).removeO()
print("small-sigma_W expansion =", chi_mix_series)
expect_zero(
    "linear mixed expansion",
    chi_mix_series - (1 - sigma_W * (sp.Rational(1, 3) + 9 * gamma_W)),
)

banner("STAGE 94 FINAL LEDGER")
print("A standalone isotropic hidden mixed pole cannot preserve the canonical even l=2 branch unless it is absent.")
print("Its linearized branch-selection triple is")
print("  (b, a_0, a_5) = (0, -sigma_W, -sigma_W gamma_W).")
