
#!/usr/bin/env python3
"""
5pn_stage75_dynamic_geometry_obstruction.py

Stage 75 audit: exact obstruction formula if the geometry lane carries dynamic even moments.
"""

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
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 75 — EXACT OBSTRUCTION FORMULA WITH DYNAMIC GEOMETRY")

omega, Omega_Q = sp.symbols("omega Omega_Q", positive=True, real=True)
Kp, Kg0, Kg2, Kg4 = sp.symbols("K_pole K_g0 K_g2 K_g4", positive=True, real=True)
eps2, eps4 = sp.symbols("eps_2 eps_4", real=True)

subbanner("75.1 — Generalized isotropic grouped-P2 + geometry ansatz")
K0 = sp.simplify(Kg0 + Kp)
K2 = sp.simplify(Kg2 + Kp / Omega_Q**2)
K4 = sp.simplify(Kg4 + Kp / Omega_Q**4)
print("K0 =", K0)
print("K2 =", K2)
print("K4 =", K4)

subbanner("75.2 — Exact branch identity with dynamic geometry")
residual = sp.simplify(K0 * K4 - 4 * K2**2)
Kg0_sol = sp.solve(sp.Eq(residual, 0), Kg0)[0]
print("K_g0 =", sp.factor(Kg0_sol))

subbanner("75.3 — Dimensionless contamination variables")
subs_eps = {
    Kg2: eps2 * Kp / Omega_Q**2,
    Kg4: eps4 * Kp / Omega_Q**4,
}
K0_eps = sp.simplify(K0.subs(Kg0, Kg0_sol).subs(subs_eps))
print("K0 =", sp.factor(K0_eps))
expected_K0 = sp.simplify(4 * Kp * (1 + eps2)**2 / (1 + eps4))
expect_zero("K0 - expected", K0_eps - expected_K0)

c_pole = sp.simplify(Kp / K0_eps)
c_geom = sp.simplify(1 - c_pole)
expected_c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2)**2))
expect_zero("c_pole - expected", c_pole - expected_c_pole)
print("c_geom =", c_geom)

subbanner("75.4 — Static-geometry limit and small-contamination expansion")
expect_zero("static-geometry c_pole - 1/4", c_pole.subs({eps2: 0, eps4: 0}) - sp.Rational(1, 4))
expect_zero("static-geometry c_geom - 3/4", c_geom.subs({eps2: 0, eps4: 0}) - sp.Rational(3, 4))

chi = sp.symbols("chi", real=True)
e2hat, e4hat = sp.symbols("e2hat e4hat", real=True)
series_c = sp.series(c_pole.subs({eps2: chi**2 * e2hat, eps4: chi**2 * e4hat}), chi, 0, 3).removeO()
expected_series = sp.simplify(sp.Rational(1, 4) * (1 + chi**2 * (e4hat - 2 * e2hat)))
expect_zero("small-contamination series", sp.expand(series_c - expected_series))

banner("STAGE 75 FINAL LEDGER")
print("If the geometry lane carries dynamic even moments, the pole fraction is")
print("  c_pole = (1 + eps_4) / [4 (1 + eps_2)^2].")
print("So the 3/4 + 1/4 split is recovered iff the geometry lane is static through")
print("O(omega^4).")
