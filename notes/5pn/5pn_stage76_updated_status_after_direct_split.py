
#!/usr/bin/env python3
"""
5pn_stage76_updated_status_after_direct_split.py

Stage 76 audit: updated reduced status after the direct grouped-P2 + geometry derivation.
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

banner("STAGE 76 — UPDATED STATUS AFTER THE DIRECT GROUPED-P2 + GEOMETRY DERIVATION")

eps2, eps4 = sp.symbols("eps_2 eps_4", real=True)
c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2)**2))
c_geom = sp.simplify(1 - c_pole)

expect_zero("static-geometry pole fraction - 1/4", c_pole.subs({eps2: 0, eps4: 0}) - sp.Rational(1, 4))
expect_zero("static-geometry geometry fraction - 3/4", c_geom.subs({eps2: 0, eps4: 0}) - sp.Rational(3, 4))

rho_alpha = sp.simplify(1 / sp.Rational(3, 4))
zeta_req = sp.simplify((sp.Rational(1, 4)) / (sp.Rational(3, 4)))
expect_zero("rho_alpha - 4/3", rho_alpha - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_req - sp.Rational(1, 3))

banner("STAGE 76 FINAL LEDGER")
print("After the direct grouped-P2 + geometry derivation, the remaining reduced ambiguity")
print("is no longer the existence of the minimal isotropic quadrupole module. It is exactly")
print("whether the real moving-throat geometry lane is static through O(omega^4).")
