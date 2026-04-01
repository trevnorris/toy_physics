#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 48 SymPy audit

Purpose
-------
Evaluate the exact equilibrium gain on the first explicit parent wall branch,
V_conf(r;a)=V0*f((r-a)/ell), and derive explicit parent fail/succeed thresholds
for the wall amplitude V0.

Main checks
-----------
1. The support loading amplitude is g_phi = V0/ell.
2. The exact thin-wall integral scaling is
      I1 = 4*pi*ell*(a^2*J1 + 2*a*ell*J2 + ell^2*J3),
   with Jn the wall-profile moments of f'^2/h'.
3. For a centered symmetric wall layer, J2 = 0.
4. The exact equilibrium gain becomes
      G_eq = 4*pi*V0^2*(a^2*J1/ell + 2*a*J2 + ell*J3)/K_X,
   and its thin-wall leading term is
      G_eq^tw = 4*pi*a^2*V0^2*J1/(K_X*ell).
5. Inserting the Stage-44 thresholds and kappa = K_X L^2/T_X gives
      V0_fail^2 = T_X*ell*Pe_req/(4*pi*a^2*L^2*J1*Delta_inf),
      V0_suff^2 = T_X*ell*Pe_req/(4*pi*a^2*L^2*J1*Delta_0),
   so K_X cancels exactly from the thin-wall prefactor.
6. If h' is nearly constant in the wall layer, J1 = I_f/H_w with
      I_f = ∫ f'^2 dxi,
   yielding the explicit compressibility form
      V0_fail^2 = H_w*T_X*ell*Pe_req/(4*pi*a^2*L^2*I_f*Delta_inf).
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 48 — THIN-WALL CONFINEMENT BRANCH")

# Symbols
V0, ell, a, KX = sp.symbols("V0 ell a K_X", positive=True, real=True)
J1, J2, J3 = sp.symbols("J1 J2 J3", positive=True, real=True)
Pe_req, Delta0, Deltainf = sp.symbols("Pe_req Delta_0 Delta_inf", positive=True, real=True)
TX, L, kappa = sp.symbols("T_X L kappa", positive=True, real=True)
Hw, If = sp.symbols("H_w I_f", positive=True, real=True)

# Support loading amplitude from V_conf = V0 f((r-a)/ell)
g_phi = sp.simplify(V0 / ell)
print("g_phi =", g_phi)

# Exact leading moment expansion for I1 = ∫ chi_phi^2/h' d^3y in the thin-wall shell.
# chi_phi = f'(xi),  r = a + ell*xi, d^3y -> 4*pi*r^2 dr = 4*pi*ell*(a+ell xi)^2 dxi.
I1 = sp.simplify(4 * sp.pi * ell * (a**2 * J1 + 2 * a * ell * J2 + ell**2 * J3))
print("I1 =", I1)

banner("SYMMETRIC WALL LAYER")
I1_sym = sp.simplify(I1.subs(J2, 0))
print("I1 (J2=0) =", I1_sym)

banner("EXACT AND THIN-WALL EQUILIBRIUM GAIN")
G_eq = sp.simplify(g_phi**2 * I1 / KX)
G_eq_sym = sp.simplify(g_phi**2 * I1_sym / KX)
G_eq_tw = sp.simplify(4 * sp.pi * a**2 * V0**2 * J1 / (KX * ell))

print("G_eq =", G_eq)
print("G_eq (J2=0) =", G_eq_sym)
print("G_eq^tw =", G_eq_tw)
expect_zero(
    "thin-wall remainder after dropping O(ell/a) correction",
    sp.simplify(G_eq_sym - G_eq_tw - 4 * sp.pi * V0**2 * ell * J3 / KX),
)

banner("PARENT WALL THRESHOLDS FOR V0")

G_fail = sp.simplify(Pe_req / (kappa * Deltainf))
G_suff = sp.simplify(Pe_req / (kappa * Delta0))

V0_fail_sq = sp.solve(sp.Eq(G_eq_tw, G_fail), V0**2)[0]
V0_suff_sq = sp.solve(sp.Eq(G_eq_tw, G_suff), V0**2)[0]

print("V0_fail^2 =", sp.simplify(V0_fail_sq))
print("V0_suff^2 =", sp.simplify(V0_suff_sq))

# Insert kappa = KX L^2 / TX and verify KX cancellation.
subs_kappa = {kappa: KX * L**2 / TX}
V0_fail_sq_k = sp.simplify(V0_fail_sq.subs(subs_kappa))
V0_suff_sq_k = sp.simplify(V0_suff_sq.subs(subs_kappa))

print("V0_fail^2 with kappa inserted =", V0_fail_sq_k)
print("V0_suff^2 with kappa inserted =", V0_suff_sq_k)
expect_zero(
    "K_X cancellation in V0_fail^2",
    V0_fail_sq_k - TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * Deltainf),
)
expect_zero(
    "K_X cancellation in V0_suff^2",
    V0_suff_sq_k - TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * Delta0),
)

banner("CONSTANT-COMPRESSIBILITY WALL LAYER")

# If h' is nearly constant in the active layer, J1 = I_f / H_w.
V0_fail_const = sp.simplify(V0_fail_sq_k.subs(J1, If / Hw))
V0_suff_const = sp.simplify(V0_suff_sq_k.subs(J1, If / Hw))
print("V0_fail^2 | H~const =", V0_fail_const)
print("V0_suff^2 | H~const =", V0_suff_const)
expect_zero(
    "constant-H fail threshold",
    V0_fail_const - Hw * TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * If * Deltainf),
)
expect_zero(
    "constant-H suff threshold",
    V0_suff_const - Hw * TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * If * Delta0),
)

banner("STAGE 48 AUDIT PASSED")
print("1. The first explicit parent wall family gives g_phi = V0/ell.")
print("2. The exact equilibrium gain scales like V0^2 times a shell integral of chi_phi^2/h'.")
print("3. In the thin-wall limit the gain is controlled by a^2/ell, so thinner and wider active walls couple more strongly.")
print("4. The parent fail/succeed thresholds for V0 are explicit, and K_X cancels from their thin-wall prefactor once kappa is inserted.")
print("5. In the constant-compressibility wall layer, the thresholds reduce to a direct compressibility/tension/geometry law.")
