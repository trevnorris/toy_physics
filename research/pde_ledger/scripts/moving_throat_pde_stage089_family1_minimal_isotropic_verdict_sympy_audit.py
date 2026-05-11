#!/usr/bin/env python3
"""
moving_throat_pde_stage72_family1_minimal_isotropic_verdict_sympy_audit.py

SymPy / arithmetic audit for Stage 72.

Checks:
1. compare rho_alpha = 4/3 against the explicit Family-1 ratio window;
2. compare zeta_req = 1/3 against the explicit support ceiling;
3. verify that zeta_req < A_F1, so the explicit Family-1 transport map requires Pe_req = 0.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr, tol: float = 1e-12) -> None:
    val = sp.N(sp.simplify(expr), 50)
    print(f"{name} = {val}")
    if abs(complex(val)) > tol:
        raise AssertionError(f"{name} is not within tolerance {tol}")

banner("STAGE 72 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH")

rho_min = sp.Rational(4,3)
zeta_min = sp.Rational(1,3)

# Stage-62 Family-1 demand map data.
Pe = sp.symbols("Pe", positive=True, real=True)
y = sp.symbols("y", real=True)
kappa_F1 = sp.Rational(12321, 5)
eta_F1 = sp.Integer(37)
y_F1 = sp.nsolve(y * sp.tan(y) - eta_F1, sp.Float("1.53", 80), tol=1e-40, maxsteps=100, prec=100)
A_F1 = sp.simplify((kappa_F1 + sp.pi**2 / 4) / (kappa_F1 + y_F1**2))
Omega = sp.simplify(
    sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi)
    / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1))
)
zeta_F1 = sp.simplify(A_F1 * Omega**2)
zeta_max = sp.simplify(sp.limit(zeta_F1, Pe, sp.oo))

# Stage-63/69 thresholds evaluated at lambda_mu = 1.
Pe_suff_chi = sp.Float("96.5285247264386")
Pe_fail_chi = sp.Float("11220.5441626259")
eps_blk = sp.Integer(0)
zeta = sp.symbols("zeta", positive=True, real=True)
Q = sp.simplify((1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta))
zeta_suff = sp.simplify(zeta_F1.subs(Pe, Pe_suff_chi))
zeta_fail = sp.simplify(zeta_F1.subs(Pe, Pe_fail_chi))
rho_suff = sp.simplify(Q.subs(zeta, zeta_suff))
rho_fail = sp.simplify(Q.subs(zeta, zeta_fail))
rho_max = sp.simplify(Q.subs(zeta, zeta_max))

expect_zero("Stage-62 zeta_max = A_F1 pi^2/4", zeta_max - A_F1 * sp.pi**2 / 4, tol=1e-30)
expect_zero("Stage-69 Q(zeta;0) = 1 + zeta", Q - (1 + zeta), tol=1e-30)
expect_zero("rho_suff anchor", rho_suff - (1 + zeta_suff), tol=1e-30)
expect_zero("rho_fail anchor", rho_fail - (1 + zeta_fail), tol=1e-30)
expect_zero("rho_max anchor", rho_max - (1 + zeta_max), tol=1e-30)

Delta_suff = sp.N(rho_suff - rho_min, 25)
Delta_fail = sp.N(rho_fail - rho_min, 25)
Delta_max  = sp.N(rho_max  - rho_min, 25)
Delta_zeta = sp.N(zeta_max - zeta_min, 25)
Delta_AF1  = sp.N(A_F1 - zeta_min, 25)

print("rho_min   =", sp.N(rho_min, 25))
print("rho_suff  =", sp.N(rho_suff, 25))
print("rho_fail  =", sp.N(rho_fail, 25))
print("rho_max   =", sp.N(rho_max, 25))
print("zeta_min  =", sp.N(zeta_min, 25))
print("zeta_max  =", sp.N(zeta_max, 25))
print("A_F1      =", sp.N(A_F1, 25))

print("\nMargins:")
print("Delta_suff =", Delta_suff)
print("Delta_fail =", Delta_fail)
print("Delta_max  =", Delta_max)
print("Delta_zeta =", Delta_zeta)
print("Delta_AF1  =", Delta_AF1)

if not (rho_min < rho_suff < rho_fail < rho_max):
    raise AssertionError("Family-1 loading-ratio ordering failed.")
if not (zeta_min < 1):
    raise AssertionError("Minimal isotropic branch left the symmetric-lowest-twin regime.")
if not (zeta_min < A_F1):
    raise AssertionError("Minimal isotropic branch no longer succeeds at zero transport bias.")
if not (zeta_min < zeta_max):
    raise AssertionError("Minimal isotropic branch exceeded the Family-1 support ceiling.")

print("\nRegime checks:")
print("  rho_min < rho_suff   -> guaranteed success")
print("  zeta_min < 1         -> symmetric lowest twin already enough")
print("  zeta_min < A_F1      -> Pe_req = 0 on the explicit Family-1 transport map")
