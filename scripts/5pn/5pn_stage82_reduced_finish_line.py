
#!/usr/bin/env python3
"""
5pn_stage82_reduced_finish_line.py

Stage 82 audit: the reduced finish line after the geometry-lane check.
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

banner("STAGE 82 — THE REDUCED FINISH LINE AFTER THE GEOMETRY-LANE CHECK")

G, c, c_s, a, Omega_Q = sp.symbols("G c c_s a Omega_Q", positive=True, real=True)
Kbar0_target = sp.simplify(64 * G * Omega_Q**5 / (45 * c**5))
Kbar0_target_geom = sp.simplify(Kbar0_target.subs(Omega_Q, 3 * c_s / (2 * a)))
expect_zero("Kbar0^target - 54 G c_s^5/(5 a^5 c^5)", Kbar0_target_geom - 54 * G * c_s**5 / (5 * a**5 * c**5))

N_Q = sp.symbols("N_Q", positive=True, real=True)
print("N_Q := Kbar_0 / Kbar_0^target")
print("The full reduced GR-like point-particle 2.5PN closure on the actual isotropic branch")
print("is equivalent to N_Q = 1.")

banner("STAGE 82 FINAL LEDGER")
print("On the actual natural isotropic branch we now have:")
print("  eps_2 = eps_4 = 0,")
print("  Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2),")
print("  rho_alpha = 4/3,")
print("  zeta_req = 1/3,")
print("and the explicit Family-1 support/source side is automatic.")
print("So the only remaining reduced theorem gap is the single normalization defect N_Q - 1.")
