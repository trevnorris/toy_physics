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


def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("SAME-CHARGE BARRIER AUDIT — STAGE 023")
banner("SELECTED-BRANCH LOADING RATIO FROM THE MINIMAL ISOTROPIC QUADRUPOLE PRECURSOR")

subbanner("I. Exact contact-plus-pole inverse formulas")
rho_alpha, c0, c1 = sp.symbols("rho_alpha c0 c1", positive=True, real=True)
c0_expr = sp.simplify(1 / rho_alpha)
c1_expr = sp.simplify((rho_alpha - 1) / rho_alpha)
print(f"c0(rho_alpha) = {c0_expr}")
print(f"c1(rho_alpha) = {c1_expr}")

rho_from_c0 = sp.simplify(1 / c0)
rho_from_c1 = sp.simplify(1 / (1 - c1))
zeta_req = sp.simplify(c1 / c0)
print(f"rho_alpha from c0 = {rho_from_c0}")
print(f"rho_alpha from c1 = {rho_from_c1}")
print(f"zeta_req = c1/c0 = {zeta_req}")

expect_zero(
    "inverse check: c0 -> rho_alpha",
    sp.simplify(rho_from_c0.subs(c0, c0_expr) - rho_alpha),
)
expect_zero(
    "inverse check: c1 -> rho_alpha",
    sp.simplify(rho_from_c1.subs(c1, c1_expr) - rho_alpha),
)
expect_zero(
    "zeta_req from c0,c1",
    sp.simplify(zeta_req.subs({c0: c0_expr, c1: c1_expr}) - (rho_alpha - 1)),
)

subbanner("II. Matching to the minimal isotropic quadrupole module")
vals = {c0: sp.Rational(3, 4), c1: sp.Rational(1, 4)}
rho_selected = sp.simplify(rho_from_c0.subs(vals))
zeta_selected = sp.simplify(zeta_req.subs(vals))
print(f"rho_alpha on minimal isotropic module = {rho_selected}")
print(f"zeta_req on minimal isotropic module = {zeta_selected}")
expect_zero("rho_alpha - 4/3", rho_selected - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_selected - sp.Rational(1, 3))

subbanner("III. Exact support-selector form of the selected branch")
Lambda, epsilon_star, C_mix = sp.symbols("Lambda epsilon_star C_mix", positive=True, real=True)
Pi_tr = sp.simplify(sp.Rational(4, 3) * C_mix)
varrho = sp.simplify(sp.pi**2 * Pi_tr / (16 * Lambda))
C_mix_expr = sp.simplify(8 * Lambda * (1 - epsilon_star) / sp.pi**2)
varrho_expr = sp.simplify(varrho.subs(C_mix, C_mix_expr))
S_req = sp.simplify(Pi_tr / C_mix)
print(f"Pi_tr selected = {Pi_tr}")
print(f"C_mix = {C_mix_expr}")
print(f"varrho selected = {varrho_expr}")
print(f"S_req = {S_req}")
expect_zero("varrho - 2(1-epsilon_*)/3", varrho_expr - sp.Rational(2, 3) * (1 - epsilon_star))
expect_zero("S_req - 4/3", S_req - sp.Rational(4, 3))

subbanner("IV. Regime meaning")
regime_left = sp.simplify(Pi_tr - C_mix)
regime_right = sp.simplify(2 * C_mix - Pi_tr)
print(f"Pi_tr - C_mix = {regime_left}")
print(f"2 C_mix - Pi_tr = {regime_right}")
expect_zero("Pi_tr - C_mix - C_mix/3", regime_left - C_mix / 3)
expect_zero("2C_mix - Pi_tr - 2C_mix/3", regime_right - 2 * C_mix / 3)

subbanner("V. Equivalent varrho form")
Pi_from_varrho = sp.simplify(16 * Lambda * varrho_expr / sp.pi**2)
expect_zero(
    "Pi_tr recovered from varrho",
    sp.simplify(Pi_from_varrho - sp.Rational(4, 3) * C_mix_expr),
)

print("\nAll Stage 023 checks passed.")
