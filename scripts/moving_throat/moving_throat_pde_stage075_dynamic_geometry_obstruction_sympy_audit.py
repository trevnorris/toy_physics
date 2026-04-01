#!/usr/bin/env python3
"""
Stage 75 SymPy audit.

Derive the exact obstruction formula for the grouped-P2 + geometry split when the
geometry lane carries O(omega^2) and O(omega^4) moments.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    s = sp.simplify(sp.expand(expr))
    print(f"{name} = {s}")
    if s != 0:
        raise AssertionError(f"{name} is not zero")


omega = sp.symbols("omega", real=True)
Kg0, Kg2, Kg4 = sp.symbols("K_g0 K_g2 K_g4", real=True)
Kp, OmegaQ = sp.symbols("K_p Omega_Q", positive=True, real=True)
eps2, eps4 = sp.symbols("eps_2 eps_4", real=True)

banner("STAGE 75 — DYNAMIC-GEOMETRY OBSTRUCTION")

K = sp.simplify(Kg0 + Kg2 * omega**2 + Kg4 * omega**4 + Kp / (1 - omega**2 / OmegaQ**2))
series = sp.expand(sp.series(K, omega, 0, 6).removeO())
K0 = sp.simplify(series.coeff(omega, 0))
K2 = sp.simplify(series.coeff(omega, 2))
K4 = sp.simplify(series.coeff(omega, 4))

print("K_Q^cons(omega) =", K)
print("Series =", series)
print("K0 =", K0)
print("K2 =", K2)
print("K4 =", K4)

branch = sp.simplify(K0 * K4 - 4 * K2**2)
print("Branch obstruction =", branch)

Kg0_sol = sp.solve(sp.Eq(branch, 0), Kg0)[0]
print("K_g0 on branch =", sp.factor(Kg0_sol))

expect_zero(
    "static-geometry limit gives 3 K_pole",
    Kg0_sol.subs({Kg2: 0, Kg4: 0}) - 3 * Kp,
)

cpole = sp.simplify(Kp / (Kg0_sol + Kp))
print("c_pole =", sp.factor(cpole))

cpole_dimless = sp.simplify(
    cpole.subs({Kg2: eps2 * Kp / OmegaQ**2, Kg4: eps4 * Kp / OmegaQ**4})
)
print("c_pole in (eps2, eps4) variables =", sp.factor(cpole_dimless))

cpole_expected = sp.simplify((1 + eps4) / (4 * (1 + eps2)**2))
expect_zero("c_pole - (1+eps4)/(4(1+eps2)^2)", cpole_dimless - cpole_expected)

cgeom_dimless = sp.simplify(1 - cpole_dimless)
print("c_geom in (eps2, eps4) variables =", sp.factor(cgeom_dimless))

small_series = sp.expand(sp.series(cpole_expected, eps2, 0, 2).removeO())
small_series = sp.expand(sp.series(small_series, eps4, 0, 2).removeO())
linear_part = sp.Rational(1, 4) * (1 + eps4 - 2 * eps2)
remainder = sp.expand(small_series - linear_part)
print("First-order expansion of c_pole =", small_series)
print("Linear part                =", linear_part)
print("Dropped higher-order tail  =", remainder)

banner("FINAL LEDGER")
print("With dynamic geometry contamination,")
print("  c_pole = (1 + eps4) / [4 (1 + eps2)^2].")
print("So the 3/4 + 1/4 split is recovered iff eps2 = eps4 = 0.")
