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

# Anchor (4): lower-vs-singular branch discrimination — paper Checks item 3
# (paper/stages/stage_131.tex:24). Pi_* must sit on the LOWER compensated branch
# g_-^{F1}, NOT on the singular equal-normalized branch (g_nat = 1) nor the upper
# branch g_+^{F1}. g_Pi(Pi) rises monotonically from 2/pi to a supremum of 1, so
# g_nat = 1 is the UNREACHABLE supremum ("singular") and g_+^{F1} > 1 is never attained.
g_nat = sp.Integer(1)                                    # equal-normalized branch, notes Sec. 1
g_plus_exact = (2*sp.sqrt(4107 - 100*sp.pi**2) + 37*sp.sqrt(3)) / (20*sp.pi)  # upper branch
gPi_at_star = sp.N(gPi.subs(Pi, Pi_star), 40)

# (4a) lower-branch MEMBERSHIP: Pi_* solves g_Pi = g_-^{F1} to high precision.
lower_residual = abs(sp.N(gPi_at_star - g_minus, 40))
assert lower_residual < sp.Float("1e-30", 40), (
    f"lower-branch membership failed: g_Pi(Pi_*) - g_-^F1 = {lower_residual}"
)
print(f"PASS: Pi_* on lower branch — |g_Pi(Pi_*) - g_-^F1| = {lower_residual}")

# (4b) SINGULAR equal-normalized branch EXCLUDED: separation matches notes Delta g_-.
sing_sep = sp.N(g_nat - gPi_at_star, 30)
delta_g_minus_notes = sp.Float("0.241964921055337", 30)   # notes stage122 Sec. 4 line 104
assert abs(sp.N(sing_sep - delta_g_minus_notes, 30)) < sp.Float("1e-12", 30), (
    f"singular-branch separation mismatch: g_nat - g_Pi(Pi_*) = {sing_sep}, "
    f"notes Delta g_- = {delta_g_minus_notes}"
)
assert sing_sep > sp.Float("1e-3", 30), (
    f"singular equal-normalized branch NOT excluded: g_nat - g_Pi(Pi_*) = {sing_sep}"
)
print(f"PASS: singular equal-normalized branch excluded — g_nat - g_Pi(Pi_*) = {sing_sep}")

# (4c) UPPER branch EXCLUDED: g_Pi(Pi_*) is far below g_+^{F1}.
upper_sep = abs(sp.N(gPi_at_star - g_plus_exact, 30))
assert upper_sep > sp.Float("1", 30), (
    f"upper branch NOT excluded: |g_Pi(Pi_*) - g_+^F1| = {upper_sep}"
)
print(f"PASS: upper branch excluded — |g_Pi(Pi_*) - g_+^F1| = {upper_sep}")
