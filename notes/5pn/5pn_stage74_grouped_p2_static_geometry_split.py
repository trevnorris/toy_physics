
#!/usr/bin/env python3
"""
5pn_stage74_grouped_p2_static_geometry_split.py

Stage 74 audit: deriving the 3/4 + 1/4 conservative quadrupole module from the
grouped-P2 + static-geometry split.
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

banner("STAGE 74 — GROUPED-P2 + STATIC-GEOMETRY DERIVATION OF THE 3/4 + 1/4 SPLIT")

omega, Omega_Q = sp.symbols("omega Omega_Q", positive=True, real=True)
K_geom, K_pole = sp.symbols("K_geom K_pole", positive=True, real=True)

subbanner("74.1 — Minimal grouped-P2 + geometry realization")
K0 = sp.simplify(K_geom + K_pole)
K2 = sp.simplify(K_pole / Omega_Q**2)
K4 = sp.simplify(K_pole / Omega_Q**4)
print("K0 =", K0)
print("K2 =", K2)
print("K4 =", K4)

subbanner("74.2 — Branch identity forces the split")
branch_residual = sp.simplify(K0 * K4 - 4 * K2**2)
print("K0 K4 - 4 K2^2 =", branch_residual)
K_geom_sol = sp.solve(sp.Eq(branch_residual, 0), K_geom)[0]
print("K_geom =", K_geom_sol)
expect_zero("K_geom - 3 K_pole", K_geom_sol - 3 * K_pole)

c_pole = sp.simplify(K_pole / K0.subs(K_geom, K_geom_sol))
c_geom = sp.simplify(K_geom_sol / K0.subs(K_geom, K_geom_sol))
expect_zero("c_pole - 1/4", c_pole - sp.Rational(1, 4))
expect_zero("c_geom - 3/4", c_geom - sp.Rational(3, 4))

subbanner("74.3 — Normalized conservative quadrupole module")
K_cons = sp.simplify(K_geom + K_pole / (1 - omega**2 / Omega_Q**2))
Yhat = sp.simplify((K_cons / K0).subs(K_geom, K_geom_sol))
Yhat_target = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / Omega_Q**2))
expect_zero("Yhat_Q^cons - target", sp.simplify(Yhat - Yhat_target))

subbanner("74.4 — Support/source corollary")
rho_alpha = sp.simplify(1 / sp.Rational(3, 4))
zeta_req = sp.simplify((sp.Rational(1, 4)) / (sp.Rational(3, 4)))
expect_zero("rho_alpha - 4/3", rho_alpha - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_req - sp.Rational(1, 3))

banner("STAGE 74 FINAL LEDGER")
print("If the isotropic grouped-P2 conservative branch is carried by one effective pole and")
print("the geometry lane is genuinely static through O(omega^4), then the minimal isotropic")
print("branch identity forces")
print("  K_pole = K0/4,   K_geom = 3 K0/4,")
print("and therefore")
print("  Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2).")
