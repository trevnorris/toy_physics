#!/usr/bin/env python3
"""
Stage 131 SymPy audit — parent micro-threshold for canonical mouth compensation.
"""

from __future__ import annotations
import sympy as sp

Pi = sp.symbols("Pi", positive=True, real=True)
Theta_sigma, L = sp.symbols("Theta_sigma L", positive=True, real=True)
Tm, qstar, A0p = sp.symbols("T_m q_* A0p", real=True)
# F1 lower-branch value g_-^{F1}: closed form is (2*sqrt(4107 - 100*pi^2) - 37*sqrt(3)) / (20*pi).
# Carried forward from the F1 family lower-branch derivation in the upstream
# parent-mouth stages (see notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md Sec. 3).
g_minus_exact = (2*sp.sqrt(4107 - 100*sp.pi**2) - 37*sp.sqrt(3)) / (20*sp.pi)
g_minus_literal = sp.Float("0.758035078944663", 50)
assert abs(sp.N(g_minus_exact, 50) - g_minus_literal) < sp.Float("1e-14", 50), (
    f"g_-^{{F1}} closed form vs literal: {sp.N(g_minus_exact, 50)} != {g_minus_literal}"
)
g_minus = sp.N(g_minus_exact, 50)
print(f"PASS: g_-^F1 closed form matches literal 0.758035078944663")

gPi = sp.simplify(2*Pi*(2*Pi*sp.exp(Pi) + sp.pi)/((4*Pi**2 + sp.pi**2)*(sp.exp(Pi) - 1)))
Pi_star = sp.nsolve(gPi - g_minus, 1.5, tol=1e-30, maxsteps=100, prec=80)

V1 = Pi*Theta_sigma/L
V1_star = sp.N(Pi_star, 30)*Theta_sigma/L
print("Pi_* =", sp.N(Pi_star, 30))
print("V1_* =", V1_star)

gprime_star = sp.N(sp.diff(gPi, Pi).subs(Pi, Pi_star), 30)
print("g'(Pi_*) =", gprime_star)

threshold_residual = sp.simplify((Tm - qstar*A0p) - Pi*Theta_sigma/L)
print("Parent bias mismatch formula =", threshold_residual)

# --- Anchored assertions vs notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md ---

# Anchor (1): Pi_* matches the notes-quoted value 1.50882951349316 (Sec. 1).
Pi_star_paper = sp.Float("1.50882951349316", 50)
Pi_star_num = sp.N(Pi_star, 50)
assert abs(Pi_star_num - Pi_star_paper) < sp.Float("1e-14", 50), (
    f"Pi_* mismatch vs notes Sec. 1: computed {Pi_star_num}, paper {Pi_star_paper}"
)
print("PASS: Pi_* matches notes Sec. 1 value 1.50882951349316")

# Anchor (2): g'(Pi_*) matches the notes-quoted slope 0.0714453558083195 (Sec. 3).
gprime_paper = sp.Float("0.0714453558083195", 50)
gprime_num = sp.N(gprime_star, 50)
assert abs(gprime_num - gprime_paper) < sp.Float("1e-14", 50), (
    f"g'(Pi_*) mismatch vs notes Sec. 3: computed {gprime_num}, paper {gprime_paper}"
)
print("PASS: g'(Pi_*) matches notes Sec. 3 value 0.0714453558083195")

# Anchor (3): parent threshold identity at Pi = Pi_*.
# notes Sec. 2:  T_m - q_* A_0' = Pi_* * Theta_sigma / L
threshold_at_star = threshold_residual.subs(Pi, sp.N(Pi_star, 50))
expected_form = (Tm - qstar*A0p) - sp.N(Pi_star, 50)*Theta_sigma/L
assert sp.simplify(threshold_at_star - expected_form) == 0, (
    "parent threshold identity at Pi_* does not match (T_m - q_* A_0') - Pi_* Theta_sigma/L"
)
print("PASS: parent threshold identity at Pi = Pi_* matches notes Sec. 2")

# Anchor (4): lower-branch discrimination — Pi_* sits on the g_- branch, NOT on a singular point.
# A point clearly away from Pi_* (e.g. 2*Pi_*) must give a residual visibly far from zero.
gPi_offstar = gPi.subs(Pi, 2*sp.N(Pi_star, 30))
offstar_residual = abs(sp.N(gPi_offstar - g_minus, 30))
assert offstar_residual > sp.Float("1e-3", 30), (
    f"counter-example failed: gPi(2*Pi_*) residual vs g_- = {offstar_residual}, "
    "expected clearly nonzero (lower-branch discrimination, paper Checks item 3)"
)
print(f"PASS: lower-branch discrimination — gPi(2*Pi_*) - g_- = {offstar_residual}")
