#!/usr/bin/env python3
"""
moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py

SymPy-backed audit for Stage 167.

Checks:
1. Exact bundle transport formulas for c_{s,w}, ell, L_W, v_{w0}, T_m.
2. Exact bundle transport formulas for g_s, g_q, lambda.
3. Exact first-order invariance of r_c, frak r, frak g.
4. Exact vanishing of the Stage 163 off-family channels and delta_perp.
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

banner("STAGE 167 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION")

dTheta, dKs, dKq, dP = sp.symbols("dTheta dKs dKq dP", real=True)

# Stage 166 bundle inversion
drho = sp.Rational(1, 2) * dTheta
da = sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta
dcs = sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta + sp.Rational(1, 5) * dP
dZ = dKq - sp.Rational(2, 5) * dP

# Frozen n=5 wall-EOS and healing lock
dcsw = 2 * drho
dell = -dcsw
dLW = da

print("delta ln rho_w  =", drho)
print("delta ln a      =", da)
print("delta ln c_s    =", dcs)
print("delta ln Z_q    =", dZ)
print("delta ln c_s,w  =", dcsw)
print("delta ln ell    =", dell)
print("delta ln L_W    =", dLW)

# Stage 165 exact branch drifts
dv = sp.simplify(sp.Rational(1, 2) * (dZ - drho) + sp.Rational(3, 2) * dcsw + dcs - sp.Rational(5, 2) * da)
dT = sp.simplify(sp.Rational(1, 2) * (dZ - drho) + sp.Rational(3, 2) * dcsw - dcs - sp.Rational(3, 2) * da)

print("\ndelta ln v_w0   =", dv)
print("delta ln T_m    =", dT)
print("delta ln(v/T)   =", sp.simplify(dv - dT))
print("delta ln(v*T)   =", sp.simplify(dv + dT))

# Parent couplings
dgq = sp.simplify(dZ - sp.Rational(3, 2) * dLW)                  # g_q ∝ Z_q L_W^{-3/2}
dgs = sp.simplify(dT + 2 * da + dell)                            # g_s ∝ T_m a^2 ell
dIsq = sp.simplify(2 * da + dell + sp.Rational(1, 2) * dLW)     # I_sq ∝ a^2 ell L_W^{1/2}
dlam = sp.simplify(dv + dIsq)                                    # lambda ∝ v_w0 I_sq

print("\ndelta ln g_q    =", dgq)
print("delta ln g_s    =", dgs)
print("delta ln I_sq   =", dIsq)
print("delta ln lambda =", dlam)

banner("Parent invariants")
drc = sp.simplify(2 * dlam - dKs - dKq)
dr = sp.simplify(dlam - sp.Rational(1, 2) * (dKs + dKq))
dg = sp.simplify(dgq + sp.Rational(1, 2) * dKs - dgs - sp.Rational(1, 2) * dKq)

expect_zero("delta ln r_c", drc)
expect_zero("delta ln frak r", dr)
expect_zero("delta ln frak g", dg)

banner("Stage 163 off-family channels")
chan1 = sp.simplify(dgq + dKs - dgs - dlam)
chan2 = sp.simplify(dKs + dKq - 2 * dlam)
expect_zero("delta ln(g_q K_s/(g_s lambda))", chan1)
expect_zero("delta ln(K_s K_q/lambda^2)", chan2)

gstar, rstar = sp.symbols("gstar rstar", positive=True, real=True)
delta_perp = sp.simplify(gstar * chan1 + chan2 / (4 * sp.sqrt(1 + rstar**2)))
expect_zero("delta_perp", delta_perp)

print("\nCarry-forward formulas:")
print("  delta ln v_w0   = -3/4 delta ln K_s + 1/2 delta ln K_q + 13/8 delta ln Theta_w")
print("  delta ln T_m    = -5/4 delta ln K_s + 1/2 delta ln K_q + 15/8 delta ln Theta_w - 2/5 delta ln P_0")
print("  delta ln g_s    = -1/4 delta ln K_s + 1/2 delta ln K_q + 3/8 delta ln Theta_w - 2/5 delta ln P_0")
print("  delta ln g_q    = -3/4 delta ln K_s + delta ln K_q + 3/8 delta ln Theta_w - 2/5 delta ln P_0")
print("  delta ln lambda = 1/2 (delta ln K_s + delta ln K_q)")
print("  delta ln r_c    = delta ln frak r = delta ln frak g = 0")
print("  delta_perp      = 0")
