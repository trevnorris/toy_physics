#!/usr/bin/env python3
"""0PN stationary-mouth interaction energy flux witness from primitive fields.

Benchmark primitives (n=3):
  v_i = q_i (x-X_i) / (rho * 4*pi * |x-X_i|^3)
  v = grad(phi), p = p_inf - rho*(dphi/dt + |v|^2/2)
  E = rho*|v|^2/2
  stationary cut: J = (E+p) v  =>  (p_inf - rho*dphi/dt) v
  four-state: dQ = Q12 - Q1 - Q2 + Q0
  A = q1*q2 / (rho * 4*pi * d^2)   # attraction magnitude, same-sign q
  prediction: delta J_1 = A * s  to first order (fixed sphere about mouth 1)

Checks: speed reversal, half strength, and a control that drops unsteady
Bernoulli (must break the flux prediction).
Force is NOT used to compute the flux; A is only a reference scale.
"""
from __future__ import annotations

import math
import resource
import sys
import traceback

# ---- memory guard (no wall-clock deadline) ----
SOFT_MB, HARD_MB = 512, 512

def install_memory_guard(soft_mb=SOFT_MB, hard_mb=HARD_MB):
    soft = soft_mb * 1024 * 1024
    hard = hard_mb * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (soft, hard))

# ---- primitives ----
PI = math.pi
OMEGA_3 = 4.0 * PI

def v_monopole(q, X, x, rho):
    """Lab velocity of one radial monopole at x."""
    dx = x[0] - X[0]
    dy = x[1] - X[1]
    dz = x[2] - X[2]
    r2 = dx*dx + dy*dy + dz*dz
    r = math.sqrt(r2)
    if r < 1e-30:
        raise ValueError("hit singularity")
    pref = q / (rho * OMEGA_3 * r2 * r)  # / r^3
    return (pref*dx, pref*dy, pref*dz)

def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def fibonacci_sphere(n):
    """Approx. equal-area unit directions + weights (4pi/n)."""
    pts = []
    w = OMEGA_3 / n
    offset = 2.0 / n
    increment = PI * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = ((i * offset) - 1.0) + (offset / 2.0)
        r = math.sqrt(max(0.0, 1.0 - y*y))
        phi = (i % n) * increment
        x = math.cos(phi) * r
        z = math.sin(phi) * r
        pts.append(((x, y, z), w))
    return pts

def energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                           include_mouth1, include_mouth2, use_unsteady,
                           n_quad=4000):
    """Outward flux ∮ (p_inf - rho*dphi/dt) v·n dA on sphere about `center`.

    Mouth 1 fixed; mouth 2 translates at Xdot2 = (-s, 0, 0).
    dphi/dt = s * (v2)_x  when mouth 2 present and use_unsteady.
    """
    flux = 0.0
    for nhat, w_ang in fibonacci_sphere(n_quad):
        x = (center[0] + radius*nhat[0],
             center[1] + radius*nhat[1],
             center[2] + radius*nhat[2])
        area_weight = w_ang * radius * radius  # dA = R^2 dΩ

        v = (0.0, 0.0, 0.0)
        dphi_dt = 0.0
        if include_mouth1:
            v1 = v_monopole(q1, X1, x, rho)
            v = add(v, v1)
        if include_mouth2:
            v2 = v_monopole(q2, X2, x, rho)
            v = add(v, v2)
            if use_unsteady:
                # dphi/dt = s * (v2)_x
                dphi_dt += s * v2[0]

        # (E+p) = p_inf - rho*dphi/dt   (E cancels)
        Ep = p_inf - rho * dphi_dt
        vn = dot(v, nhat)  # n outward from enclosed region
        flux += Ep * vn * area_weight
    return flux

def four_state_J1(q1, q2, d, rho, p_inf, s, radius, use_unsteady, n_quad=4000):
    X1 = (0.0, 0.0, 0.0)
    X2 = (d, 0.0, 0.0)
    center = X1
    Q12 = energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                                 True, True, use_unsteady, n_quad)
    Q1  = energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                                 True, False, use_unsteady, n_quad)
    Q2  = energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                                 False, True, use_unsteady, n_quad)
    Q0  = energy_flux_stationary(q1, q2, X1, X2, rho, p_inf, s, center, radius,
                                 False, False, use_unsteady, n_quad)
    return Q12 - Q1 - Q2 + Q0

def A_of_d(q1, q2, rho, d):
    return q1 * q2 / (rho * OMEGA_3 * d * d)

def run_checks():
    # Drains: outward rates negative, same channel sign
    q1 = -2.0
    q2 = -3.0
    rho = 1.25
    p_inf = 0.0
    d = 4.0
    R_enc = 0.7   # companion outside: d - R_enc > 0
    assert d - R_enc > 1.0
    s = 0.4
    A = A_of_d(q1, q2, rho, d)
    assert A > 0.0

    n_quad = 6000
    # Proper prediction: unsteady Bernoulli on
    J_p = four_state_J1(q1, q2, d, rho, p_inf, +s, R_enc, True, n_quad)
    J_m = four_state_J1(q1, q2, d, rho, p_inf, -s, R_enc, True, n_quad)
    J_h = four_state_J1(q1, q2, d, rho, p_inf, +0.5*s, R_enc, True, n_quad)
    J_0 = four_state_J1(q1, q2, d, rho, p_inf, 0.0, R_enc, True, n_quad)

    pred_p = A * s
    pred_h = A * (0.5 * s)
    odd_slope = (J_p - J_m) / (2.0 * s)

    # Breaking control: drop dphi/dt in Bernoulli (steady pressure only)
    J_break = four_state_J1(q1, q2, d, rho, p_inf, +s, R_enc, False, n_quad)

    # Second radius (convergence / cut independence)
    J_p_R2 = four_state_J1(q1, q2, d, rho, p_inf, +s, 1.1, True, n_quad)

    rel = lambda a, b: abs(a - b) / max(1e-30, abs(b))

    results = {
        "A": A,
        "J(+s)": J_p,
        "pred A*s": pred_p,
        "rel_err (+s)": rel(J_p, pred_p),
        "J(-s)": J_m,
        "reversal_sum": J_p + J_m,          # ~0
        "J(+s/2)": J_h,
        "pred A*(s/2)": pred_h,
        "rel_err half": rel(J_h, pred_h),
        "J(0)": J_0,
        "odd_slope": odd_slope,
        "rel_err slope vs A": rel(odd_slope, A),
        "J_break_drop_unsteady": J_break,
        "break_rel_to_As": rel(J_break, pred_p),  # should be O(1), not <<1
        "J(+s) R=1.1": J_p_R2,
        "rel_err second R": rel(J_p_R2, pred_p),
    }

    # Pass criteria
    tol = 0.03  # quadrature tolerance for fibonacci n~6e3
    passes = (
        results["rel_err (+s)"] < tol
        and abs(results["reversal_sum"]) < tol * abs(pred_p)
        and results["rel_err half"] < tol
        and abs(results["J(0)"]) < tol * abs(pred_p)
        and results["rel_err slope vs A"] < tol
        and results["rel_err second R"] < tol
        and results["break_rel_to_As"] > 0.5  # control must NOT match As
    )

    print("=== 0PN fixed-mouth energy-flux witness ===")
    for k, v in results.items():
        print(f"{k}: {v!r}")
    print("PASS" if passes else "FAIL")
    return 0 if passes else 1

if __name__ == "__main__":
    try:
        install_memory_guard()
        sys.exit(run_checks())
    except MemoryError:
        print("Memory guard triggered", file=sys.stderr)
        sys.exit(2)
    except Exception:
        traceback.print_exc()
        sys.exit(3)
