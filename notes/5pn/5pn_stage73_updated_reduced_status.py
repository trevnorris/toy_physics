
#!/usr/bin/env python3
"""
5pn_stage73_updated_reduced_status.py

Stage 73 audit: updated reduced theorem status after the loading-ratio extraction.
"""

from __future__ import annotations
import mpmath as mp
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

banner("STAGE 73 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION")

rho_alpha, zeta_req = sp.symbols("rho_alpha zeta_req", positive=True, real=True)
c0 = sp.Rational(3, 4)
c1 = sp.Rational(1, 4)

subbanner("73.1 — Minimal isotropic contact-plus-pole data")
rho_expr = sp.simplify(1 / c0)
zeta_expr = sp.simplify(c1 / c0)
expect_zero("rho_alpha - 4/3", rho_expr - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_expr - sp.Rational(1, 3))

subbanner("73.2 — Explicit Family-1 margin")
mp.mp.dps = 50
zeta_max_F1 = mp.mpf("2.46752922945601")
A_F1 = mp.mpf("1.0000519288021953286593340837139871812220538900677990")
zeta_req_min = mp.mpf(1) / 3
print("zeta_req^(min) =", zeta_req_min)
print("A_F1 - zeta_req^(min) =", A_F1 - zeta_req_min)
print("zeta_max^(F1) - zeta_req^(min) =", zeta_max_F1 - zeta_req_min)
print("Conclusion: the explicit Family-1 branch already succeeds at Pe_req = 0.")

banner("STAGE 73 FINAL LEDGER")
print("After the loading-ratio extraction, the explicit support/source side is no longer the")
print("live reduced bottleneck:")
print("  rho_alpha = 4/3,")
print("  zeta_req = 1/3,")
print("  Pe_req = 0 on the explicit Family-1 branch.")
print("So the next theorem gate is the grouped-P2 / geometry realization of the minimal")
print("isotropic conservative quadrupole module.")
