
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

"""
5pn_stage112_explicit_mouth_boundary_layer.py

Stage 112 audit: explicit GNLS + localized-Maxwell mouth boundary layer.
"""

banner("STAGE 112 — EXPLICIT MOUTH BOUNDARY LAYER")

z, L = sp.symbols("z L", positive=True, real=True)
Theta_sigma, V1, sigma_star, Pi = sp.symbols("Theta_sigma V_1 sigma_* Pi", positive=True, real=True)

sigma = sp.Function("sigma")
mu_chem = Theta_sigma * sp.log(sigma(z) / sigma_star) + V1 * z
print("mu_sigma^chem(z) =")
sp.pprint(mu_chem)

# Zero-flux branch: Theta sigma' + V1 sigma = 0
sigma_sol = sp.simplify(Pi * sp.exp(-Pi * z / L) / (L * (1 - sp.exp(-Pi))))
print("sigma_Pi(z) =")
sp.pprint(sigma_sol)

expect_zero(
    "normalization",
    sp.integrate(sigma_sol, (z, 0, L)) - 1,
)
expect_zero(
    "zero-flux ODE",
    sp.simplify(Theta_sigma * sp.diff(sigma_sol, z) + V1 * sigma_sol).subs(Pi, V1 * L / Theta_sigma),
)

banner("STAGE 112 FINAL LEDGER")
print("The explicit positive mouth-source branch is")
print("  sigma_Pi(z) = Pi e^{-Pi z/L} / [L(1 - e^{-Pi})],")
print("with")
print("  Pi = V_1 L / Theta_sigma.")
print("So the earlier truncated exponential family is the exact zero-flux branch of")
print("a GNLS + localized-Maxwell mouth boundary layer.")
