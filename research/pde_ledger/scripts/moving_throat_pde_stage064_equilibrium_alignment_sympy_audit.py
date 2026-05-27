#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 47 SymPy audit

Purpose
-------
Verify the exact algebra behind the parent equilibrium source/support alignment law.

Main checks
-----------
1. A local static linear-response closure implies
      chi_sigma(y) = g_phi * chi_phi(y) / H(y),
   with H(y) = h'(rho_*(y)).
2. The overlap invariants become
      O = g_phi * I1,
      N_ss = g_phi**2 * I2,
   where
      I1 = ∫ chi_phi^2 / H,
      I2 = ∫ chi_phi^2 / H^2.
3. The coherence factor is
      C^2 = I1^2 / (N_pp * I2).
4. In the constant-compressibility limit H(y)=H_w, the branch is exactly matched:
      C^2 = 1,
      G_eq = g_phi^2 * I1 / K_X = g_phi^2 N_pp / (K_X H_w),
   reproducing the Stage-45/46 best-alignment formulas.
5. A discrete two-point model verifies the Cauchy gap exactly:
      N_pp I2 - I1^2 = w1*w2*(H1 - H2)^2/(H1^2 H2^2) >= 0.
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


banner("STAGE 47 — PARENT EQUILIBRIUM SOURCE/SUPPORT ALIGNMENT")

# Exact symbolic parent overlap variables.
g_phi, KX, Npp = sp.symbols("g_phi K_X N_pp", positive=True, real=True)
I1, I2 = sp.symbols("I1 I2", positive=True, real=True)
Hw = sp.symbols("H_w", positive=True, real=True)

# --- Local static linear-response closure ---
# At fixed y, the local sigma-energy is (1/2) H(y) sigma^2 - g_phi chi_phi(y) sigma.
# Minimising over sigma yields the closure law sigma(y) = g_phi chi_phi(y) / H(y).
y_loc = sp.symbols("y_loc", real=True)
chi_phi_loc = sp.Function("chi_phi")(y_loc)
H_loc = sp.Function("H")(y_loc)
sigma_loc = sp.symbols("sigma_loc", real=True)
F_loc = sp.Rational(1, 2) * H_loc * sigma_loc**2 - g_phi * chi_phi_loc * sigma_loc
closure_solutions = sp.solve(sp.diff(F_loc, sigma_loc), sigma_loc)
assert len(closure_solutions) == 1, "linear-response minimiser must be unique"
chi_sigma_closure = closure_solutions[0]
print("closure chi_sigma =", chi_sigma_closure)
expect_zero("closure law chi_sigma = g_phi chi_phi/H", chi_sigma_closure - g_phi * chi_phi_loc / H_loc)

# --- Integral-level overlap invariants (concrete Gaussian profile) ---
y_int, L_int = sp.symbols("y_int L_int", positive=True, real=True)
chi_phi_g = sp.exp(-y_int**2 / (2 * L_int**2))
H_g = Hw
Npp_int_check = sp.integrate(chi_phi_g**2, (y_int, -sp.oo, sp.oo))
I1_int_check = sp.integrate(chi_phi_g**2 / H_g, (y_int, -sp.oo, sp.oo))
I2_int_check = sp.integrate(chi_phi_g**2 / H_g**2, (y_int, -sp.oo, sp.oo))
chi_sigma_g = g_phi * chi_phi_g / H_g
Osp_int_check = sp.integrate(chi_sigma_g * chi_phi_g, (y_int, -sp.oo, sp.oo))
Nss_int_check = sp.integrate(chi_sigma_g**2, (y_int, -sp.oo, sp.oo))
print("Npp_int =", Npp_int_check, "  I1_int =", I1_int_check, "  I2_int =", I2_int_check)
expect_zero("overlap O = g_phi * I1 (integral form)", Osp_int_check - g_phi * I1_int_check)
expect_zero("overlap N_ss = g_phi^2 * I2 (integral form)", Nss_int_check - g_phi**2 * I2_int_check)

# Derived overlap invariants from chi_sigma = g_phi * chi_phi / H.
Osp = sp.simplify(g_phi * I1)
Nss = sp.simplify(g_phi**2 * I2)
C2 = sp.simplify(Osp**2 / (Npp * Nss))
Geq = sp.simplify(g_phi**2 * I1 / KX)

print("O_(sigma phi) =", Osp)
print("N_(sigma sigma) =", Nss)
print("C_(sigma phi)^2 =", C2)
print("G_eq =", Geq)

banner("MATCHED-LAYER INTEGRAL REDUCTION (concrete Gaussian profile)")

y, L = sp.symbols("y L", positive=True, real=True)
chi_phi_y = sp.exp(-y**2 / (2 * L**2))
H_const = Hw  # constant compressibility on the active layer

Npp_int = sp.integrate(chi_phi_y**2, (y, -sp.oo, sp.oo))
I1_int = sp.integrate(chi_phi_y**2 / H_const, (y, -sp.oo, sp.oo))
I2_int = sp.integrate(chi_phi_y**2 / H_const**2, (y, -sp.oo, sp.oo))

print("Npp_int =", Npp_int)
print("I1_int  =", I1_int)
print("I2_int  =", I2_int)
expect_zero("matched-layer I1 reduction", I1_int - Npp_int / Hw)
expect_zero("matched-layer I2 reduction", I2_int - Npp_int / Hw**2)

banner("CONSTANT-COMPRESSIBILITY / MATCHED-LAYER LIMIT")

# In the thin active layer with nearly constant H(y)=H_w:
# I1 = Npp/Hw, I2 = Npp/Hw^2.
const_subs = {I1: Npp / Hw, I2: Npp / Hw**2}
C2_const = sp.simplify(C2.subs(const_subs))
Geq_const = sp.simplify(Geq.subs(const_subs))

print("C^2 | H=const =", C2_const)
print("G_eq | H=const =", Geq_const)
expect_zero("matched-layer coherence", C2_const - 1)
expect_zero("matched-layer gain vs Stage-45 best-alignment formula", Geq_const - g_phi**2 * Npp / (KX * Hw))

banner("DISCRETE TWO-POINT CAUCHY GAP CHECK")

w1, w2 = sp.symbols("w1 w2", positive=True, real=True)
H1, H2 = sp.symbols("H1 H2", positive=True, real=True)

Npp_disc = sp.simplify(w1 + w2)
I1_disc = sp.simplify(w1 / H1 + w2 / H2)
I2_disc = sp.simplify(w1 / H1**2 + w2 / H2**2)

gap_disc = sp.simplify(Npp_disc * I2_disc - I1_disc**2)

gap_expected = sp.simplify(w1 * w2 * (H1 - H2)**2 / (H1**2 * H2**2))

print("N_pp I2 - I1^2 =", gap_disc)
print("expected gap     =", gap_expected)
expect_zero("two-point Cauchy gap identity", gap_disc - gap_expected)

banner("GENERAL-H EQUILIBRIUM SOFTENING CHECK")

# General-H equilibrium softening on a two-point parent-equilibrium-aligned branch.
# On the aligned branch chi_sigma_k = g_phi chi_phi_k / H_k, so
#   Theta   = sum_k H_k chi_sigma_k^2     = g_phi^2 sum_k chi_phi_k^2 / H_k = g_phi^2 I_1
#   Lambda  = g_phi sum_k chi_phi_k chi_sigma_k = g_phi^2 sum_k chi_phi_k^2 / H_k = g_phi^2 I_1
# Therefore Lambda^2/Theta = g_phi^2 I_1 for ANY H(y), with no matched-layer assumption.
chi_sig1 = g_phi * sp.sqrt(w1) / H1  # chi_sigma_k for w_k = chi_phi_k^2
chi_sig2 = g_phi * sp.sqrt(w2) / H2
Theta_general = sp.simplify(H1 * chi_sig1**2 + H2 * chi_sig2**2)
Lambda_general = sp.simplify(g_phi * (sp.sqrt(w1) * chi_sig1 + sp.sqrt(w2) * chi_sig2))
soft_general = sp.simplify(Lambda_general**2 / Theta_general)
print("Theta (general, two-point) =", Theta_general)
print("Lambda (general, two-point) =", Lambda_general)
print("Lambda^2/Theta (general, two-point) =", soft_general)
expect_zero(
    "general equilibrium softening equals g_phi^2 I_1",
    soft_general - g_phi**2 * I1_disc,
)

banner("ELIMINATED-SOURCE SOFTENING CHECK")

# Reduced static energy F = (Theta/2) sigma^2 - Lambda sigma phi + (KX/2) phi^2
Theta, Lambda, phi = sp.symbols("Theta Lambda phi", positive=True, real=True)
sigma = sp.symbols("sigma", real=True)
F = sp.Rational(1, 2) * Theta * sigma**2 - Lambda * sigma * phi + sp.Rational(1, 2) * KX * phi**2
sigma_stat = sp.solve(sp.Eq(sp.diff(F, sigma), 0), sigma)[0]
F_eff = sp.simplify(F.subs(sigma, sigma_stat))
print("sigma_stat =", sigma_stat)
print("F_eff =", F_eff)
expect_zero(
    "effective support softening",
    F_eff - sp.Rational(1, 2) * (KX - Lambda**2 / Theta) * phi**2,
)

banner("STAGE 47 AUDIT PASSED")
print("1. The equilibrium-induced source channel is aligned with the support loading through 1/H.")
print("2. The coherence factor is exactly C^2 = I1^2/(N_pp I2) and reaches 1 when H is constant on the active layer.")
print("3. The exact eliminated-source gain is G_eq = g_phi^2 I1/K_X.")
print("4. In the matched-layer limit this reproduces the Stage-45/46 best-alignment formulas.")
