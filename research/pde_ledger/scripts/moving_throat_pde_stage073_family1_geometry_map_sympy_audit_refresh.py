#!/usr/bin/env python3
"""
moving_throat_pde_stage56_family1_geometry_map_sympy_audit.py

SymPy audit for Stage 56.

Verifies the explicit Family-1 geometry map:
- epsilon_r = 1/20,
- carried aspect ratio L/a = 37/20,
- hence Lambda_ell = (L/a)/(ell/a) = 37,
- and under the local Robin mouth closure K_m = T_X/ell,
  the Robin variable eta = K_m L / T_X also equals 37.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 56 — FAMILY-1 GEOMETRY MAP")

# Symbolic identity: Lambda_ell = (L/a) / (ell/a) = L/ell, independent of the
# specific reference-branch values chosen below.
L_sym, a_sym, ell_sym = sp.symbols('L a ell', positive=True)
Lambda_star_sym = L_sym / a_sym
ell_over_a_sym = ell_sym / a_sym
Lambda_ell_sym = sp.simplify(Lambda_star_sym / ell_over_a_sym)
expect_zero("Lambda_ell - L/ell (symbolic)", Lambda_ell_sym - L_sym / ell_sym)

epsilon_r = sp.Rational(1, 20)
Lambda_star = sp.Rational(37, 20)  # carried reference branch L/a

ell_over_a = epsilon_r
Lambda_ell = sp.simplify(Lambda_star / ell_over_a)

print("epsilon_r =", epsilon_r)
print("ell/a     =", ell_over_a)
print("L/a       =", Lambda_star)
print("Lambda_ell = L/ell =", Lambda_ell)

expect_zero("Lambda_ell - 37", Lambda_ell - 37)

# Robin closure K_m = T_X/ell. Build eta symbolically in K_m first so the
# assertion eta - L/ell == 0 actually exercises the closure substitution and
# not a trivial T_X cancellation.
K_m, T_X, L, ell = sp.symbols('K_m T_X L ell', positive=True, real=True)
eta_sym = K_m * L / T_X
eta = sp.simplify(eta_sym.subs(K_m, T_X / ell))
print("eta under K_m = T_X/ell ->", eta)
expect_zero("eta - L/ell", eta - L / ell)
expect_zero("eta(reference) - 37", eta.subs({L / ell: 37}) - 37)

banner("FINAL LEDGER")
print("Lambda_ell = 37")
print("eta        = 37")
