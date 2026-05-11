#!/usr/bin/env python3
"""
moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py

Numerical audit of the actual Family-1 mouth correction induced by the full compensated
mouth potential.
"""

from __future__ import annotations
import sympy as sp
import mpmath as mp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_close(name: str, val: float, target: float, tol: float = 1e-10) -> None:
    diff = abs(val - target)
    print(f"{name} = {val}   (target {target}, diff {diff})")
    if diff > tol:
        raise AssertionError(f"{name} mismatch")

banner("ACTUAL FAMILY-1 MOUTH CORRECTION")

Pi_star = mp.mpf("1.50882951349316")
Sigma_m_star = mp.mpf("0.451485277739090")
g_star = mp.mpf("0.758035078944663")
S_star = mp.mpf("0.658075937605429")
gprime_star = mp.mpf("0.0714453558083195")
AT = mp.mpf("-4.27263956256927")
BT = mp.mpf("0.134875005736706")
Tstar = mp.mpf("0.901484054174205")

k = mp.pi/2

def Sigma_star(x):
    return Pi_star * mp.e**(-Pi_star * x) / (1 - mp.e**(-Pi_star))

def T_s(x):
    return (1 - mp.e**(-Pi_star * x)) / (Pi_star * (1 - mp.e**(-Pi_star))) - x * mp.e**(-Pi_star)/(1 - mp.e**(-Pi_star))

def T_q(x):
    Cq = Pi_star / ((1 - mp.e**(-Pi_star)) * (k**2 - Pi_star**2))
    Aq = Cq * (k * mp.sinh(k) + Pi_star * mp.e**(-Pi_star)) / (k * mp.cosh(k))
    return Aq * mp.sinh(k*x) - Cq * mp.cosh(k*x) + Cq * mp.e**(-Pi_star*x)

def S_q_kernel(x):
    return mp.cosh((mp.pi/2)*(1-x)) / mp.cosh(mp.pi/2)

def c_kernel(x):
    return mp.cos(mp.pi*x/2)

# Canonical compensated residual profile
def R_star(x):
    return Sigma_m_star * (4*T_s(x) - T_q(x)) - Pi_star*x

# Check the tangent constraints numerically.
for h in [mp.mpf("1e-6"), mp.mpf("1e-5")]:
    deriv = (R_star(h) - R_star(0)) / h
    print(f"R'(0) forward diff at h={h} :", deriv)
print("R(0) =", R_star(0))
print("R(1) =", R_star(1))

# Weighted averages
def avg(func):
    return mp.quad(lambda t: Sigma_star(t) * func(t), [0, 1])

def cov(f, h):
    af = avg(f)
    ah = avg(h)
    return mp.quad(lambda t: Sigma_star(t) * (f(t) - af) * (h(t) - ah), [0, 1])

Cov_cR = cov(c_kernel, R_star)
Cov_SR = cov(S_q_kernel, R_star)

delta_g = -Cov_cR
delta_S = -Cov_SR
deltaPi = Cov_cR / gprime_star
deltaT = AT * delta_g + BT * delta_S

print("Cov_*(c,R_*)      =", Cov_cR)
print("Cov_*(K_q,R_*)    =", Cov_SR)
print("delta g_act       =", delta_g)
print("delta S_act       =", delta_S)
print("delta Pi_act      =", deltaPi)
print("delta Tm_act      =", deltaT)
print("Pi_corr           =", Pi_star + deltaPi)
print("T_corr            =", Tstar + deltaT)

# One-step nonlinear Picard iterate
def sigma1_unnorm(x):
    return mp.e**(-Pi_star * x - R_star(x))

Z1 = mp.quad(lambda t: sigma1_unnorm(t), [0, 1])
def sigma1(x):
    return sigma1_unnorm(x) / Z1

g1 = mp.quad(lambda t: sigma1(t) * c_kernel(t), [0, 1])
S1 = mp.quad(lambda t: sigma1(t) * S_q_kernel(t), [0, 1])

deltaPi1 = -(g1 - g_star) / gprime_star
deltaT1 = AT * (g1 - g_star) + BT * (S1 - S_star)

print("g_1               =", g1)
print("S_1               =", S1)
print("Pi_1              =", Pi_star + deltaPi1)
print("T_1               =", Tstar + deltaT1)

# Effective interpolation parameter from Stage 148 affine laws
lam_Pi = (mp.mpf("1.69941496131430") - deltaPi) / mp.mpf("2.08240814741023")
lam_T = (mp.mpf("0.508756302215085") - deltaT) / mp.mpf("0.625700104366894")
print("lambda_eff^(Pi)   =", lam_Pi)
print("lambda_eff^(T)    =", lam_T)

# Internal consistency checks
expect_close("delta g from sigma1-linearized covariance check",
             mp.quad(lambda t: Sigma_star(t)*(1 - (R_star(t)-avg(R_star))) * c_kernel(t), [0,1]) - g_star,
             delta_g, tol=1e-9)

print("\nConclusion:")
print("  The full compensated mouth potential broadens the source relative to the tangent exponential,")
print("  shifts the canonical Family-1 point upward, and is well approximated by a moderate positive")
print("  broadening toward the uniform family (lambda_eff ≈ 0.38).")
