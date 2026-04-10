
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
5pn_stage103_core_parameter_status.py

Stage 103 audit: status reduction of the concrete core parameters to the
dimensionless parent controls.
"""

banner("STAGE 103 — CORE-PARAMETER EXTRACTION STATUS")

alpha = sp.symbols("alpha", positive=True, real=True)
K_s, K_q, lam = sp.symbols("K_s K_q lambda", positive=True, real=True)
g_s, g_q = sp.symbols("g_s g_q", positive=True, real=True)
a = sp.symbols("a", positive=True, real=True)

r = sp.simplify(lam / sp.sqrt(K_s * K_q))
g = sp.simplify(g_q * sp.sqrt(K_s) / (g_s * sp.sqrt(K_q)))
L_W = sp.simplify(sp.pi * a * sp.sqrt((1 + r**2) / 3) / 2)

# The canonical scale transformation of the core leaves the normalized branch data invariant.
subs_scale = {
    K_s: alpha * K_s,
    K_q: alpha * K_q,
    lam: alpha * lam,
    g_s: sp.sqrt(alpha) * g_s,
    g_q: sp.sqrt(alpha) * g_q,
}
expect_zero("mathfrak_r scale invariance", r.subs(subs_scale) - r)
expect_zero("mathfrak_g scale invariance", g.subs(subs_scale) - g)
expect_zero("L_W/a depends only on mathfrak_r", (L_W / a).subs(subs_scale) - (L_W / a))

balance = sp.simplify(1 + r**2 - 4 * (g - r)**2)
print("balance condition =", balance)

banner("STAGE 103 FINAL LEDGER")
print("After explicit parent extraction, the surviving microscopic branch controls are:")
print("  K_s, K_q, mathfrak_r, mathfrak_g, L_W.")
print("But the compensated canonical outlet depends only on")
print("  mathfrak_r = lambda/sqrt(K_s K_q)")
print("  mathfrak_g = g_q sqrt(K_s)/(g_s sqrt(K_q))")
print("together with")
print("  L_W/a = (pi/2) sqrt((1 + mathfrak_r^2)/3).")
print("So the live microscopic theorem gate is no longer 'what are the outlet coefficients?'")
print("It is only: which branch point (mathfrak_r, mathfrak_g) does the real GNLS + localized-Maxwell core select?")
