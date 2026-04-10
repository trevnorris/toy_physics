#!/usr/bin/env python3
"""
5pn_stage97_concrete_two_channel_core_outlet.py

Stage 97 audit: concrete two-channel core outlet model.
"""

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

banner("STAGE 97 — CONCRETE TWO-CHANNEL CORE OUTLET MODEL")

Ks, Kq, lam = sp.symbols("K_s K_q lambda", positive=True, real=True)
gs, gq = sp.symbols("g_s g_q", real=True)
kappa0, gamma0, z = sp.symbols("kappa_0 gamma_0 z", real=True)
u = sp.symbols("u", real=True)

Dw = 1 - kappa0 * z**2 - sp.I * gamma0 * z**5
M = sp.Matrix([[Ks, lam], [lam, -Kq * Dw]])
rhs = u * sp.Matrix([gs, gq])
sol = sp.simplify(M.LUsolve(rhs))
s_sol, q_sol = sol
delta = sp.simplify((gs * s_sol + gq * q_sol) / u)

print("s(u,z) =")
sp.pprint(s_sol / u)
print("q(u,z) =")
sp.pprint(q_sol / u)
print("delta Lambda_core(z) =")
sp.pprint(delta)

delta_formula = sp.simplify(
    gs**2 / Ks
    - (Ks * gq - lam * gs)**2 / (Ks * (Ks * Kq * Dw + lam**2))
)
expect_zero("exact Schur-complement formula", delta - delta_formula)

rc = sp.simplify(lam**2 / (Ks * Kq))
rho_c = sp.simplify(gs**2 / Ks)
sigma_c = sp.simplify((Ks * gq - lam * gs)**2 / (Ks**2 * Kq * (1 + rc)))
kappa_c = sp.simplify(kappa0 / (1 + rc))
gamma_c = sp.simplify(gamma0 / (1 + rc))

delta_reduced = sp.simplify(rho_c - sigma_c / (1 - kappa_c * z**2 - sp.I * gamma_c * z**5))
print("rho_c   =", rho_c)
print("sigma_c =", sp.factor(sigma_c))
print("kappa_c =", kappa_c)
print("gamma_c =", gamma_c)
expect_zero("reduced Robin-mixed form", delta - delta_reduced)

banner("STAGE 97 FINAL LEDGER")
print("The concrete two-channel core model reproduces the full reduced Robin–mixed hybrid outlet exactly.")
print("Core-level identifications:")
print("  rho_c   = g_s^2/K_s")
print("  sigma_c = (K_s g_q - lambda g_s)^2 / [K_s^2 K_q (1+r_c)]")
print("  kappa_c = kappa_0/(1+r_c)")
print("  gamma_c = gamma_0/(1+r_c)")
print("with r_c = lambda^2/(K_s K_q).")
