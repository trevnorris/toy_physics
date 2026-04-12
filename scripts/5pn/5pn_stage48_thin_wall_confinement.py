#!/usr/bin/env python3
"""
5pn_stage48_thin_wall_confinement.py

Stage 48 audit: explicit thin-wall confinement branch and parent thresholds for the wall amplitude.
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


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 48 — EXPLICIT THIN-WALL CONFINEMENT BRANCH")

a, ell, V0, K_X, T_X, L = sp.symbols("a ell V0 K_X T_X L", positive=True, real=True)
Pe_req, Delta0, Deltainf, kappa = sp.symbols("Pe_req Delta0 Deltainf kappa", positive=True, real=True)
J1, J2, J3, I_f, H_w = sp.symbols("J1 J2 J3 I_f H_w", positive=True, real=True)
r = sp.symbols("r", real=True)
z = sp.symbols("z", real=True)
f = sp.Function("f")

subbanner("48.1 — Exact wall loading amplitude")
xi = (r - a) / ell
V_conf = V0 * f(z).subs(z, xi)
chain_rule_target = V0 * sp.diff(f(z), z).subs(z, xi) / ell
expect_zero("-dV_conf/da - (V0/ell) f'(xi)", -sp.diff(V_conf, a) - chain_rule_target)
print("So the branch loading amplitude is g_phi = V0 / ell.")

subbanner("48.2 — Exact shell integral entering I_1")
I1_exact = 4 * sp.pi * ell * (a**2 * J1 + 2 * a * ell * J2 + ell**2 * J3)
print("I_1 =")
sp.pprint(I1_exact)
print("Centered symmetric wall layer (J_2=0):")
sp.pprint(sp.simplify(I1_exact.subs(J2, 0)))

subbanner("48.3 — Exact equilibrium gain and thin-wall limit")
G_eq = sp.simplify((V0 / ell)**2 * I1_exact / K_X)
expect_zero(
    "equilibrium gain formula",
    G_eq - 4 * sp.pi * V0**2 * (a**2 * J1 / ell + 2 * a * J2 + ell * J3) / K_X,
)
G_eq_tw = sp.simplify(4 * sp.pi * a**2 * V0**2 * J1 / (K_X * ell))
print("Leading thin-wall gain =")
sp.pprint(G_eq_tw)

subbanner("48.4 — Exact wall-amplitude fail/succeed thresholds")
G_fail = Pe_req / (kappa * Deltainf)
G_suff = Pe_req / (kappa * Delta0)
V0_fail_sq = sp.simplify(K_X * ell * G_fail / (4 * sp.pi * a**2 * J1))
V0_suff_sq = sp.simplify(K_X * ell * G_suff / (4 * sp.pi * a**2 * J1))

V0_fail_explicit = sp.simplify(V0_fail_sq.subs(K_X, kappa * T_X / L**2))
V0_suff_explicit = sp.simplify(V0_suff_sq.subs(K_X, kappa * T_X / L**2))
expect_zero(
    "V0_fail explicit",
    V0_fail_explicit - T_X * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * Deltainf),
)
expect_zero(
    "V0_suff explicit",
    V0_suff_explicit - T_X * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * Delta0),
)

subbanner("48.5 — Constant-compressibility wall layer")
J1_const = I_f / H_w
print("J_1 = I_f / H_w gives thresholds")
sp.pprint(sp.simplify(V0_fail_explicit.subs(J1, J1_const)))
sp.pprint(sp.simplify(V0_suff_explicit.subs(J1, J1_const)))

banner("STAGE 48 FINAL LEDGER")
print("Stage 48 converts the parent-overlap theorem into a physical wall-amplitude test on")
print("the first explicit thin-wall confinement branch. The remaining branch question is now")
print("whether the actual wall amplitude V0 clears the exact fail/succeed surfaces above.")
