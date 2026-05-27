#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(s):
    print("\n" + "="*88)
    print(s)
    print("="*88)

def expect_zero(name, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 121 — GEOMETRIC SELECTION OF r")

L,a,LW,r = sp.symbols('L a L_W r', positive=True, real=True)

LW_formula = sp.pi*a/2 * sp.sqrt((1+r**2)/3)
r_geom = sp.solve(sp.Eq(LW, LW_formula), r)[0]
r_geom_simplified = sp.simplify(r_geom.subs(LW, L))
print("r_geom(L/a) =", r_geom_simplified)

expect_zero(
    "r_geom closed-form (explicit)",
    r_geom_simplified**2 - (12*L**2/(sp.pi**2 * a**2) - 1)
)

expect_zero(
    "tube-length relation",
    LW_formula.subs(r, r_geom_simplified) - L
)

R = sp.Rational(37,20)
r_F1 = sp.simplify(r_geom_simplified.subs({L: R*a}))
print("r_F1 =", r_F1)
print("r_F1 numeric =", sp.N(r_F1, 20))

r_F1_target = sp.sqrt(4107 - 100*sp.pi**2)/(10*sp.pi)
expect_zero(
    "r_F1 symbolic value at L/a = 37/20",
    sp.simplify(r_F1 - r_F1_target)
)

rc_F1 = sp.simplify(r_F1**2)
print("r_c(F1) =", rc_F1)
print("r_c(F1) numeric =", sp.N(rc_F1, 20))

rc_F1_target = sp.Rational(4107, 100)/sp.pi**2 - 1
expect_zero(
    "r_c(F1) symbolic value",
    sp.simplify(rc_F1 - rc_F1_target)
)

# Notes (boxed): with L_W = L, the mixed-tube pole is Omega_W = pi c_s / (2 L).
c_s = sp.symbols('c_s', positive=True, real=True)
OmegaW_formula_LW = sp.pi * c_s / (2 * LW)
OmegaW_at_LW_eq_L = OmegaW_formula_LW.subs(LW, L)
print("Omega_W(L_W = L) =", OmegaW_at_LW_eq_L)
expect_zero(
    "Omega_W identification at L_W = L",
    OmegaW_at_LW_eq_L - sp.pi*c_s/(2*L)
)

threshold = sp.simplify(sp.pi/(2*sp.sqrt(3)))
print("existence threshold L/a >=", threshold, "≈", sp.N(threshold, 20))

expect_zero(
    "r_geom vanishes at existence threshold",
    sp.simplify(r_geom_simplified.subs(L, a*threshold))
)
