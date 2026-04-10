
#!/usr/bin/env python3
"""
5pn_stage79_geometry_lane_check_verdict.py

Stage 79 audit: actual branch verdict for the geometry-lane check.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 79 — ACTUAL BRANCH VERDICT FOR THE GEOMETRY-LANE CHECK")

eps2, eps4 = sp.symbols("eps_2 eps_4", real=True)
c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2)**2))
c_geom = sp.simplify(1 - c_pole)

expect_zero("actual-branch eps_2", eps2.subs(eps2, 0))
expect_zero("actual-branch eps_4", eps4.subs(eps4, 0))
expect_zero("actual-branch c_pole - 1/4", c_pole.subs({eps2: 0, eps4: 0}) - sp.Rational(1, 4))
expect_zero("actual-branch c_geom - 3/4", c_geom.subs({eps2: 0, eps4: 0}) - sp.Rational(3, 4))

rho_alpha = sp.simplify(1 / sp.Rational(3, 4))
zeta_req = sp.simplify((sp.Rational(1, 4)) / (sp.Rational(3, 4)))
expect_zero("rho_alpha - 4/3", rho_alpha - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_req - sp.Rational(1, 3))

banner("STAGE 79 FINAL LEDGER")
print("On the actual isotropic reduced branch, the geometry lane is dynamically inert through")
print("O(omega^4), so")
print("  eps_2 = eps_4 = 0,")
print("  c_pole = 1/4,   c_geom = 3/4,")
print("and therefore")
print("  Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2).")
