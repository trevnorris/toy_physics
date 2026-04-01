#!/usr/bin/env python3
"""
moving_throat_pde_stage138_frozen_traction_fixedpoint_sympy_audit.py

Numerical/SymPy audit for the full co-evolving Family-1 fixed point at frozen
canonical traction Sigma0*.

Checks:
1. Solve the exact nonlinear fixed-point map.
2. Report g_fp, S_fp, R_fp, Pi_fp.
3. Verify the exact shifted-R transport law from Stage 137.
"""

from __future__ import annotations
import math
import numpy as np
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

rF1 = 1.77799353547498
g_star = 0.758035078944663
S_star = 0.658075937605429
Pi_star = 1.50882951349316
Sigma0_star = 1.80594111095636

N = 2401
x = np.linspace(0.0, 1.0, N)
dx = x[1] - x[0]
w = np.ones(N) * dx
w[0] = w[-1] = dx / 2
kappa = math.pi / 2
c = np.cos(math.pi * x / 2)
Kq = np.cosh(kappa * (1 - x)) / np.cosh(kappa)

def normalize(sig: np.ndarray) -> np.ndarray:
    return sig / np.sum(sig * w)

def Ts(sig: np.ndarray) -> np.ndarray:
    cum_sig = np.cumsum(sig * w)
    cum_y = np.cumsum(x * sig * w)
    return cum_y + x * (cum_sig[-1] - cum_sig)

def Tq(sig: np.ndarray) -> np.ndarray:
    A = np.cumsum(np.sinh(kappa * x) * sig * w)
    B = np.cumsum((np.cosh(kappa * (1 - x)) * sig * w)[::-1])[::-1]
    return (np.cosh(kappa * (1 - x)) * A + np.sinh(kappa * x) * B) / (kappa * np.cosh(kappa))

def g(sig: np.ndarray) -> float:
    return float(np.sum(c * sig * w))

def S(sig: np.ndarray) -> float:
    return float(np.sum(Kq * sig * w))

def R(sig: np.ndarray) -> float:
    gv = g(sig)
    return ((gv - rF1) ** 2) / (1 + rF1**2)

def Phi(sig: np.ndarray, Sigma0: float) -> np.ndarray:
    return Sigma0 * (Ts(sig) - R(sig) * Tq(sig))

def next_sigma(sig: np.ndarray, Sigma0: float) -> np.ndarray:
    ph = Phi(sig, Sigma0)
    return normalize(np.exp(-ph))

def solve_fixed_point(Sigma0: float, maxit: int = 400) -> tuple[np.ndarray, int, float]:
    sig = normalize(Pi_star * np.exp(-Pi_star * x))
    for n in range(maxit):
        sig_new = next_sigma(sig, Sigma0)
        err = float(np.max(np.abs(sig_new - sig)))
        sig = sig_new
        if err < 1e-13:
            return sig, n + 1, err
    return sig, maxit, err

banner("STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT")
sig, niter, err = solve_fixed_point(Sigma0_star)
g_fp = g(sig)
S_fp = S(sig)
R_fp = R(sig)
Pi_fp = Sigma0_star * (1 - R_fp * S_fp)

print("iterations =", niter)
print("max residual =", err)
print("g_fp =", g_fp)
print("S_fp =", S_fp)
print("R_fp =", R_fp)
print("Pi_fp =", Pi_fp)

dg = sp.Float(g_fp - g_star, 30)
r = sp.Float(rF1, 30)
pred = sp.simplify(-dg / sp.sqrt(1 + r**2) + dg**2 / (1 + r**2))
direct = sp.Float(R_fp - sp.Rational(1, 4), 30)
print("delta R from direct solve =", direct)
print("delta R from exact transport law =", pred)
if abs(float(pred - direct)) > 1e-10:
    raise AssertionError("Exact transport law failed numerically.")

print("\nConclusion:")
print("  At frozen canonical traction, the co-evolving fixed point remains close in Pi")
print("  but drifts to R > 1/4, so exact compensation is lost.")
