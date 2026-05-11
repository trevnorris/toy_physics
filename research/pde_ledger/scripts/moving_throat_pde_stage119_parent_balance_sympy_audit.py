#!/usr/bin/env python3
"""
moving_throat_pde_stage119_parent_balance_sympy_audit.py

SymPy audit for the one-parameter parent compensation family.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("I. Dimensionless reduction of the compensation law")

K_s, K_q, lam, g_s, g_q = sp.symbols("K_s K_q lam g_s g_q", positive=True, real=True)
rhat, ghat = sp.symbols("rhat ghat", real=True)

law = sp.expand(g_s**2*(K_s*K_q + lam**2) - 4*(K_s*g_q - lam*g_s)**2)
law_red = sp.simplify(
    law.subs({
        lam: rhat*sp.sqrt(K_s*K_q),
        g_q: ghat*g_s*sp.sqrt(K_q)/sp.sqrt(K_s)
    }) / (g_s**2*K_s*K_q)
)
print("Reduced law =", law_red)
expect_zero("dimensionless law", law_red - (1 + rhat**2 - 4*(ghat-rhat)**2))

banner("II. Exact solution for the mouth-coupling ratio")

ghat_sol = sp.solve(sp.Eq(1 + rhat**2, 4*(ghat-rhat)**2), ghat)
print("ghat solutions =", ghat_sol)
expect_zero("positive branch check", 1 + rhat**2 - 4*(ghat_sol[0]-rhat)**2)
expect_zero("negative branch check", 1 + rhat**2 - 4*(ghat_sol[1]-rhat)**2)

banner("III. D/N-tube length in terms of the same normalized hybridization")

a = sp.symbols("a", positive=True, real=True)
L_W, rc = sp.symbols("L_W rc", positive=True, real=True)
kappa0 = 4*L_W**2/(sp.pi**2*a**2)
L_sel = sp.solve(sp.Eq(kappa0, (1+rc)/3), L_W)[0]
print("L_W =", sp.simplify(L_sel))
expect_zero("tube-length law", L_sel - sp.pi*a*sp.sqrt((1+rc)/3)/2)

banner("IV. Explicit traction law")

Zq, mu0, c_s, Tm, J_s = sp.symbols("Zq mu0 c_s Tm J_s", positive=True, real=True)
Kq = Zq*sp.pi**2*c_s**2/(4*mu0*L_W**2)
gq = Zq*sp.pi/(sp.sqrt(2)*mu0*L_W**sp.Rational(3,2))
ghat_expr = sp.simplify(gq*sp.sqrt(K_s)/(Tm*J_s*sp.sqrt(Kq)))
print("ghat explicit =", ghat_expr)
expect_zero(
    "ghat explicit simplification",
    ghat_expr - sp.sqrt(2*Zq*K_s)/(Tm*J_s*c_s*sp.sqrt(mu0*L_W))
)

Tm_sol_plus = sp.simplify(sp.solve(sp.Eq(ghat_expr, rhat + sp.sqrt(1+rhat**2)/2), Tm)[0])
Tm_sol_minus = sp.simplify(sp.solve(sp.Eq(ghat_expr, rhat - sp.sqrt(1+rhat**2)/2), Tm)[0])
print("T_m (+ branch) =", Tm_sol_plus)
print("T_m (- branch) =", Tm_sol_minus)

print("\nAll Stage 119 symbolic checks passed.")
