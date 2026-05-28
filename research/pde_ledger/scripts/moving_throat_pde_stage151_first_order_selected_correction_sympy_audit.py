#!/usr/bin/env python3
"""
moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py

Audit for the first-order source correction selected by the full mouth profile.

The paper-claimed moment-shift identities are
    delta_g = - Cov_*(c, R),   delta_S = - Cov_*(K_q, R),
where delta_Sigma = -Sigma_* (R - <R>_*) is the linearized correction to the
canonical exponential source.  The SymPy script verifies these by:

  (1) symbolically defining delta_Sigma = -Sigma_*(R - <R>_*) (hand form), and
  (2) numerically integrating delta_g_int = int c(x) delta_Sigma dx and
      checking it matches the abstract covariance form to high precision.

Numeric Pi_star, r1, r2 anchors (Pi_star = canonical Family-1 value from notes;
r1, r2 are arbitrary non-zero example residual coefficients) make the integral
identities computable yet non-trivial:  a sign flip on the centering, a missing
factor in the moment, or a misordered Sigma_* would change the numeric value
and trip the assertion.
"""

from __future__ import annotations
import mpmath as mp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_close(name: str, val, target, tol: float = 1e-15) -> None:
    diff = abs(mp.mpf(val) - mp.mpf(target))
    print(f"{name} = {val}   (target {target}, diff {diff})")
    if float(diff) > tol:
        raise AssertionError(f"{name} mismatch")


banner("FIRST-ORDER SELF-CONSISTENT SOURCE CORRECTION")

mp.mp.dps = 40

# Numeric anchors.  Pi_star = canonical Family-1 mouth bias (notes / stage 156).
# r1, r2 = arbitrary nonzero example coefficients in the residual R(x) = r1 x + r2 x^2.
# gprime, AT, BT = arbitrary nonzero example parameters; we verify the algebraic
# relations delta_Pi = -delta_g/gprime and delta_T = AT*delta_g + BT*delta_S
# without needing the specific Family-1 values (those live in stage 152).
Pi_star = mp.mpf("1.50882951349316")
r1 = mp.mpf("1.7")
r2 = mp.mpf("-0.9")
gprime = mp.mpf("0.0714453558083195")
AT = mp.mpf("-4.27263956256927")
BT = mp.mpf("0.134875005736706")
k = mp.pi / 2


def Sigma_star(x):
    return Pi_star * mp.e ** (-Pi_star * x) / (1 - mp.e ** (-Pi_star))


def c_kernel(x):
    return mp.cos(k * x)


def K_kernel(x):
    return mp.cosh(k * (1 - x)) / mp.cosh(k)


def R_residual(x):
    return r1 * x + r2 * x ** 2


def mean(f):
    return mp.quad(lambda t: Sigma_star(t) * f(t), [0, 1])


# Canonical normalization
norm = mp.quad(Sigma_star, [0, 1])
expect_close("<1>_* = 1 (canonical normalization)", norm, mp.mpf(1), tol=1e-30)

# Moments
Rbar = mean(R_residual)
cbar = mean(c_kernel)
Kbar = mean(K_kernel)
cRbar = mean(lambda t: c_kernel(t) * R_residual(t))
KRbar = mean(lambda t: K_kernel(t) * R_residual(t))
CovcR = cRbar - cbar * Rbar
CovKR = KRbar - Kbar * Rbar

print(f"<R>_*       = {Rbar}")
print(f"<c>_*       = {cbar}")
print(f"<K>_*       = {Kbar}")
print(f"Cov(c,R)    = {CovcR}")
print(f"Cov(K,R)    = {CovKR}")

# Linearized correction (hand form): delta_Sigma = -Sigma_* (R - <R>_*)
def delta_Sigma(x):
    return -Sigma_star(x) * (R_residual(x) - Rbar)


# Centering check: <delta_Sigma>_* = 0 (proves the correction is mass-preserving)
centering = mp.quad(delta_Sigma, [0, 1])
expect_close("<delta_Sigma>_*  (centering)", centering, mp.mpf(0), tol=1e-30)

# Moment-shift identities via integration
delta_g_int = mp.quad(lambda t: c_kernel(t) * delta_Sigma(t), [0, 1])
delta_S_int = mp.quad(lambda t: K_kernel(t) * delta_Sigma(t), [0, 1])
expect_close("delta_g_int = -Cov(c,R)", delta_g_int, -CovcR, tol=1e-30)
expect_close("delta_S_int = -Cov(K,R)", delta_S_int, -CovKR, tol=1e-30)

# Bias and traction retunings
deltaPi = -delta_g_int / gprime
deltaT = AT * delta_g_int + BT * delta_S_int
expect_close("deltaPi = Cov(c,R)/gprime", deltaPi, CovcR / gprime, tol=1e-30)
expect_close(
    "deltaT = -AT*Cov(c,R) - BT*Cov(K,R)",
    deltaT,
    -AT * CovcR - BT * CovKR,
    tol=1e-30,
)

print(f"\ndeltaPi  = {deltaPi}")
print(f"deltaT   = {deltaT}")

print("\nTheorem:")
print("  Once the full mouth residual R_*(x) is known, the selected first-order")
print("  source correction is completely determined by Cov_*(c,R_*) and Cov_*(K_q,R_*).")
