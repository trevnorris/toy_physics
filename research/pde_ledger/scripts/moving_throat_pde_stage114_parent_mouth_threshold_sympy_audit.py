#!/usr/bin/env python3
"""
Stage 114 SymPy audit — parent micro-threshold for canonical mouth compensation.
"""

from __future__ import annotations
import sympy as sp

Pi = sp.symbols("Pi", positive=True, real=True)
Theta_sigma, L = sp.symbols("Theta_sigma L", positive=True, real=True)
Tm, qstar, A0p = sp.symbols("T_m q_* A0p", real=True)
g_minus = sp.Float("0.758035078944663")

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
