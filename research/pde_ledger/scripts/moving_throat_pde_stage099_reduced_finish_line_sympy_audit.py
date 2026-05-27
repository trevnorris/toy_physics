#!/usr/bin/env python3
"""
Stage 099 SymPy audit: the reduced finish line is a single normalization defect.

Carry-forward annotations (per [[batch-IV1-paper-alignment]] Cluster A direction
(a)): paper card `\\stagefield{Checks}` items
  (i) "static limit eps_2=eps_4=0 returns c_pole=1/4"
  (ii) "l=0 and l=2 orthogonality before applying the geometry firewall"
  (iii) "minimal-module hypothesis on any support/source success statement"
are upstream carry-ins from Part III + IV.1 stages:
  (i)  derives from stage 091 (grouped P2 + static geometry) and stage 092
       (dynamic-geometry obstruction); the reduction is asserted explicitly at
       stage 094 (after orthogonality) and stage 096 (verdict).
  (ii) is the full content of stage 094 (15 angular integrals + Laplace
       eigenvalue) in both engines.
  (iii) is anchored by the Part III chain stages 088 / 089 / 090.
This stage exercises the static-slot value and pole residue of Yhat_Q^cons
locally as a sanity anchor on the partial-fraction form (F4 substantive check).
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


banner("STAGE 099 — REDUCED FINISH LINE")

G, c, c_s, a, Omega_Q = sp.symbols("G c c_s a Omega_Q", positive=True, real=True)
omega = sp.symbols("omega", real=True)
N_Q = sp.symbols("N_Q", positive=True, real=True)

Yhat_Q_cons = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / Omega_Q**2))
print("Yhat_Q^cons(omega) =", Yhat_Q_cons)

# Static slot: Yhat_Q^cons(omega=0) = 1
expect_zero("Yhat_Q^cons static slot equals 1",
            Yhat_Q_cons.subs(omega, 0) - 1)

# Pole structure: residue of Yhat_Q^cons at omega = Omega_Q is -Omega_Q/8.
# The (1/4)/(1 - omega^2/Omega_Q^2) part has residue -Omega_Q/8 at omega=Omega_Q,
# confirming c_pole = 1/4 in the partial-fraction sense (card check (a)).
expect_zero("Yhat_Q^cons pole residue at omega=Omega_Q is -Omega_Q/8",
            sp.residue(Yhat_Q_cons, omega, Omega_Q) - (-Omega_Q / 8))

K0_target = sp.simplify(64 * G * Omega_Q**5 / (45 * c**5))
K0_target_geom = sp.simplify(K0_target.subs(Omega_Q, 3 * c_s / (2 * a)))
expect_zero("K0_target geometric form", K0_target_geom - 54 * G * c_s**5 / (5 * a**5 * c**5))

# Conservative quadrupole structural relations (forced by Yhat_Q^cons):
# K0_sym is a *free* positive symbol — assertions are not forced by an N_Q recipe.
K0_sym = sp.symbols("K0_sym", positive=True, real=True)
K2_struct = K0_sym / (4 * Omega_Q**2)
K4_struct = K0_sym / (4 * Omega_Q**4)
Gamma5_struct = 9 * K0_sym / (32 * Omega_Q**5)

# Branch identity input (card stagefield Inputs at stage_099.tex:9):
expect_zero("branch identity K0 K4 = 4 K2^2", K0_sym * K4_struct - 4 * K2_struct**2)

# Equivalent Gamma_5 forms: derived from K_n series, and from canonical odd coeff.
expect_zero("Gamma_5 sqrt form matches canonical odd-coeff form",
            9 * K2_struct**sp.Rational(5, 2) / K0_sym**sp.Rational(3, 2) - Gamma5_struct)

# Appendix factorization eq:app-part04-factorized-defect-again (chi_Q = 1 branch):
# Gammabar_5 / (2 G / (5 c^5)) = N_Q  (with K0 = N_Q * K0_target).
expect_zero("Gamma_5 normalization equals N_Q on chi_Q=1 branch",
            Gamma5_struct.subs(K0_sym, N_Q * K0_target) / (2 * G / (5 * c**5)) - N_Q)

print("K2_struct =", sp.factor(K2_struct))
print("K4_struct =", sp.factor(K4_struct))
print("Gamma5_struct =", sp.factor(Gamma5_struct))
print("\nSTAGE 099 AUDIT PASSED")
