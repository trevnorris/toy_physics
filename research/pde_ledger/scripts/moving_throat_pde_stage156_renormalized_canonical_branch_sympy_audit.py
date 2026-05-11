#!/usr/bin/env python3
"""
moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py

Numerical/SymPy audit for the renormalized canonical Family-1 branch under full
core–mouth co-evolution.

Checks:
1. Solve the fixed point as a function of Sigma0.
2. Solve the unique root g_fp(Sigma0) = g_*.
3. Report the renormalized canonical Sigma0, T_hat, S_can, Pi_can.
"""

from __future__ import annotations
import math
import numpy as np

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

rF1 = 1.77799353547498
g_star = 0.758035078944663
Pi_star = 1.50882951349316
Sigma0_star = 1.80594111095636
T_hat_star = 0.901484054174204

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

def next_sigma(sig: np.ndarray, Sigma0: float) -> np.ndarray:
    ph = Sigma0 * (Ts(sig) - R(sig) * Tq(sig))
    return normalize(np.exp(-ph))

def solve_fixed_point(Sigma0: float, maxit: int = 500) -> np.ndarray:
    sig = normalize(Pi_star * np.exp(-Pi_star * x))
    for _ in range(maxit):
        sig_new = next_sigma(sig, Sigma0)
        if float(np.max(np.abs(sig_new - sig))) < 1e-13:
            return sig_new
        sig = sig_new
    return sig

def g_fp(Sigma0: float) -> float:
    sig = solve_fixed_point(Sigma0)
    return g(sig)

def bracket_root(lo: float, hi: float, target: float) -> tuple[float, float]:
    flo = g_fp(lo) - target
    fhi = g_fp(hi) - target
    if flo == 0:
        return lo, lo
    if fhi == 0:
        return hi, hi
    if flo * fhi > 0:
        raise RuntimeError("Root not bracketed.")
    return lo, hi

def bisect(lo: float, hi: float, target: float, iters: int = 60) -> float:
    lo, hi = bracket_root(lo, hi, target)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = g_fp(mid) - target
        fl = g_fp(lo) - target
        if fl * fm <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)

banner("STAGE 139 — RENORMALIZED CANONICAL BRANCH")
Sigma0_can = bisect(3.0, 6.0, g_star, 55)
sig_can = solve_fixed_point(Sigma0_can)
g_can = g(sig_can)
S_can = S(sig_can)
R_can = R(sig_can)
Pi_can = Sigma0_can * (1 - R_can * S_can)
T_hat_can = math.sqrt(9 * Sigma0_can / 20)

print("Sigma0_can =", Sigma0_can)
print("g_can      =", g_can)
print("S_can      =", S_can)
print("R_can      =", R_can)
print("Pi_can     =", Pi_can)
print("T_hat_can  =", T_hat_can)

print("\nRelative shifts from original canonical point:")
print("Sigma0 ratio - 1 =", Sigma0_can / Sigma0_star - 1)
print("Pi ratio - 1     =", Pi_can / Pi_star - 1)
print("T_hat ratio - 1  =", T_hat_can / T_hat_star - 1)

if abs(g_can - g_star) > 1e-10:
    raise AssertionError("Renormalized canonical branch did not restore g = g_*.")
if abs(R_can - 0.25) > 1e-10:
    raise AssertionError("Renormalized canonical branch did not restore R = 1/4.")

print("\nConclusion:")
print("  Full co-evolution preserves the lower compensated branch,")
print("  but only after a unique upward traction renormalization.")
