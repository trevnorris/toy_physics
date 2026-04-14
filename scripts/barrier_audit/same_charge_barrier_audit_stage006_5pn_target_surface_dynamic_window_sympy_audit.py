#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage006_5pn_target_surface_dynamic_window_sympy_audit.py

Stage 006 — exact 5PN isotropic target-surface compatibility and the dynamic
same-charge survival window on the primitive finite-throat one-port family.

What this script does
---------------------
1. Rebuilds the exact isotropic one-port / one-pole / normalization surface used
   by the 5PN / 2.5PN / 4PN moving-throat endgame.
2. Verifies symbolically that simultaneous isotropic one-pole success and
   outgoing-normalization success reduce to one exact compatibility equation.
3. Specializes that compatibility surface to the explicit primitive finite-throat
   same-charge branch already introduced in Stage 005.
4. Evaluates one concrete compatibility branch numerically, computes its pole
   census, and checks whether the wall-like poles still survive the Stage-001
   barrier benchmark.
5. Scans the outgoing-leg coupling lambda_W *along the 5PN-compatible surface*
   and exposes the remaining static/dynamic tradeoff as a true survival window.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.factor(sp.simplify(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.factor(sp.simplify(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def monotone_decreasing(vals: list[float]) -> bool:
    return all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))


def monotone_increasing(vals: list[float]) -> bool:
    return all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))


KAPPA = 2.0 * math.sqrt(2.0) / math.pi


# ---------------------------------------------------------------------------
# Primitive branch data containers
# ---------------------------------------------------------------------------

@dataclass
class PrimitiveParams:
    lamB: float
    lamU: float
    lamW: float
    lamR: float
    OmegaU: float
    OmegaW: float
    varpi: float
    M: float
    a: float = 1.0
    cs: float = 1.0


@dataclass
class CompatibilityData:
    C: float
    GU: float
    GW: float
    R: float
    Delta: float
    S2: float
    H: float
    Q: float
    P: float
    B0: float
    B2: float
    B4: float
    Z0: float
    Z2: float
    Z4: float
    N0: float
    P0_target_compat: float
    K_compat: float
    D0_compat: float


# ---------------------------------------------------------------------------
# Part I. Exact isotropic target-surface algebra
# ---------------------------------------------------------------------------

def verify_exact_target_surface() -> None:
    banner("PART I — EXACT 5PN ISOTROPIC TARGET-SURFACE COMPATIBILITY")

    K, M = sp.symbols("K M", positive=True, real=True)
    B0, B2, B4 = sp.symbols("B0 B2 B4", positive=True, real=True)
    Z0, Z2, Z4 = sp.symbols("Z0 Z2 Z4", real=True)
    N0 = sp.symbols("N0", positive=True, real=True)
    P0t = sp.symbols("P0_target", positive=True, real=True)

    D0 = sp.simplify(K - B0 - Z0)
    D2 = sp.simplify(-(M + B2 + Z2))
    D4 = sp.simplify(-(B4 + Z4))

    u2 = sp.simplify(-D2 / D0)
    u4 = sp.simplify((D2**2 - D0 * D4) / D0**2)
    pole_defect = sp.simplify(u4 - 4 * u2**2)

    print("D0 =")
    sp.pprint(D0)
    print("D2 =")
    sp.pprint(D2)
    print("D4 =")
    sp.pprint(D4)
    print("u2 =")
    sp.pprint(u2)
    print("u4 =")
    sp.pprint(u4)
    print("u4 - 4 u2^2 =")
    sp.pprint(pole_defect)

    expect_zero(
        "one-pole numerator identity",
        sp.simplify(pole_defect - (D0 * (B4 + Z4) - 3 * (M + B2 + Z2) ** 2) / D0**2),
    )

    K_from_pole = sp.simplify(3 * (M + B2 + Z2) ** 2 / (B4 + Z4) + B0 + Z0)
    K_from_norm = sp.simplify(N0 / P0t + B0 + Z0)
    compat = sp.simplify(N0 / P0t - 3 * (M + B2 + Z2) ** 2 / (B4 + Z4))
    P0t_compat = sp.simplify(N0 * (B4 + Z4) / (3 * (M + B2 + Z2) ** 2))

    print("K from one-pole surface =")
    sp.pprint(K_from_pole)
    print("K from normalization target =")
    sp.pprint(K_from_norm)
    print("Compatibility equation =")
    sp.pprint(compat)
    print("Branch-compatible target P0_target,compat =")
    sp.pprint(P0t_compat)

    expect_zero(
        "K_from_pole - K_from_norm on compatibility surface",
        sp.simplify((K_from_pole - K_from_norm).subs(P0t, P0t_compat)),
    )


# ---------------------------------------------------------------------------
# Part II. Primitive finite-throat specialization
# ---------------------------------------------------------------------------

def primitive_compatibility_data(p: PrimitiveParams) -> CompatibilityData:
    C = KAPPA * p.lamB
    GU = p.lamU
    GW = KAPPA * p.lamW
    R = KAPPA * p.lamR

    Delta = p.OmegaU**2 * p.OmegaW**2 - R**2
    S2 = p.OmegaU**2 + p.OmegaW**2
    H = GU**2 + GW**2
    Q = GU**2 * p.OmegaW**2 + 2.0 * GU * GW * R + GW**2 * p.OmegaU**2
    P = p.OmegaU**2 * GW + R * GU

    B0 = C**2 / p.varpi**2
    B2 = C**2 / p.varpi**4
    B4 = C**2 / p.varpi**6

    Z0 = Q / Delta
    Z2 = (Q * S2 - H * Delta) / Delta**2
    Z4 = (Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3

    N0 = P**2 / Delta**2

    P0_target_compat = N0 * (B4 + Z4) / (3.0 * (p.M + B2 + Z2) ** 2)
    K_compat = 3.0 * (p.M + B2 + Z2) ** 2 / (B4 + Z4) + B0 + Z0
    D0_compat = K_compat - B0 - Z0

    return CompatibilityData(
        C=C,
        GU=GU,
        GW=GW,
        R=R,
        Delta=Delta,
        S2=S2,
        H=H,
        Q=Q,
        P=P,
        B0=B0,
        B2=B2,
        B4=B4,
        Z0=Z0,
        Z2=Z2,
        Z4=Z4,
        N0=N0,
        P0_target_compat=P0_target_compat,
        K_compat=K_compat,
        D0_compat=D0_compat,
    )


def verify_exact_primitive_specialization() -> None:
    banner("PART II — EXACT PRIMITIVE ONE-PORT SPECIALIZATION")

    kappa = sp.simplify(2 * sp.sqrt(2) / sp.pi)
    lamB, lamU, lamW, lamR = sp.symbols("lambda_B lambda_U lambda_W lambda_R", positive=True, real=True)
    OmegaU, OmegaW, varpi, M = sp.symbols("Omega_U Omega_W varpi M", positive=True, real=True)

    C = sp.simplify(kappa * lamB)
    GU = lamU
    GW = sp.simplify(kappa * lamW)
    R = sp.simplify(kappa * lamR)

    Delta = sp.simplify(OmegaU**2 * OmegaW**2 - R**2)
    S2 = sp.simplify(OmegaU**2 + OmegaW**2)
    H = sp.simplify(GU**2 + GW**2)
    Q = sp.simplify(GU**2 * OmegaW**2 + 2 * GU * GW * R + GW**2 * OmegaU**2)
    P = sp.simplify(OmegaU**2 * GW + R * GU)

    B0 = sp.simplify(C**2 / varpi**2)
    B2 = sp.simplify(C**2 / varpi**4)
    B4 = sp.simplify(C**2 / varpi**6)
    Z0 = sp.simplify(Q / Delta)
    Z2 = sp.simplify((Q * S2 - H * Delta) / Delta**2)
    Z4 = sp.simplify((Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3)
    N0 = sp.simplify(P**2 / Delta**2)

    P0t_compat = sp.simplify(N0 * (B4 + Z4) / (3 * (M + B2 + Z2) ** 2))
    K_compat = sp.simplify(3 * (M + B2 + Z2) ** 2 / (B4 + Z4) + B0 + Z0)

    print("kappa =")
    sp.pprint(kappa)
    print("C, G_U, G_W, R =")
    sp.pprint(sp.Matrix([C, GU, GW, R]))
    print("Delta, S2, H, Q, P =")
    sp.pprint(sp.Matrix([Delta, S2, H, Q, P]))
    print("B0, B2, B4 =")
    sp.pprint(sp.Matrix([B0, B2, B4]))
    print("Z0, Z2, Z4 =")
    sp.pprint(sp.Matrix([Z0, Z2, Z4]))
    print("N0 =")
    sp.pprint(N0)
    print("P0_target,compat =")
    sp.pprint(P0t_compat)
    print("K_compat =")
    sp.pprint(K_compat)

    expect_zero(
        "primitive compatibility identity",
        sp.simplify(N0 / P0t_compat - 3 * (M + B2 + Z2) ** 2 / (B4 + Z4)),
    )


# ---------------------------------------------------------------------------
# Part III. Pole census on the compatibility branch
# ---------------------------------------------------------------------------

def quartic_coefficients(p: PrimitiveParams, K_value: float) -> list[float]:
    C = KAPPA * p.lamB
    GU = p.lamU
    GW = KAPPA * p.lamW
    R = KAPPA * p.lamR

    y = sp.symbols("y", real=True)
    F = sp.expand(
        ((K_value - p.M * y) * (p.varpi**2 - y) - C**2)
        * ((p.OmegaU**2 - y) * (p.OmegaW**2 - y) - R**2)
        - (p.varpi**2 - y)
        * (GU**2 * (p.OmegaW**2 - y) + 2.0 * GU * GW * R + GW**2 * (p.OmegaU**2 - y))
    )
    return [float(cc) for cc in sp.Poly(F, y).all_coeffs()]


def uncoupled_wall_roots(p: PrimitiveParams, K_value: float) -> list[float]:
    C = KAPPA * p.lamB
    coeff = [p.M, -(K_value + p.M * p.varpi**2), K_value * p.varpi**2 - C**2]
    roots = np.roots(coeff)
    return sorted(math.sqrt(r.real) for r in roots if abs(r.imag) < 1e-10 and r.real > 1e-12)


def uncoupled_internal_roots(p: PrimitiveParams) -> list[float]:
    R = KAPPA * p.lamR
    coeff = [1.0, -(p.OmegaU**2 + p.OmegaW**2), p.OmegaU**2 * p.OmegaW**2 - R**2]
    roots = np.roots(coeff)
    return sorted(math.sqrt(r.real) for r in roots if abs(r.imag) < 1e-10 and r.real > 1e-12)


def transfer_factor_at_pole(omega: float, p: PrimitiveParams) -> tuple[float, float, float]:
    GU = p.lamU
    GW = KAPPA * p.lamW
    R = KAPPA * p.lamR
    A = p.OmegaU**2 - omega**2
    W = p.OmegaW**2 - omega**2
    Delta = A * W - R**2
    P = A * GW + R * GU
    N = P**2 / Delta**2
    return N, Delta, P


def pole_census(p: PrimitiveParams, K_value: float) -> list[dict[str, float | str]]:
    coeff = quartic_coefficients(p, K_value)
    roots = np.roots(coeff)
    ys = sorted(r.real for r in roots if abs(r.imag) < 1e-8 and r.real > 1e-10)
    poles = [math.sqrt(y) for y in ys]
    walls = uncoupled_wall_roots(p, K_value)
    ints = uncoupled_internal_roots(p)

    out: list[dict[str, float | str]] = []
    for omega in poles:
        nearest_wall = min(abs(omega - z) for z in walls)
        nearest_int = min(abs(omega - z) for z in ints)
        family = "wall-like" if nearest_wall <= nearest_int else "internal-like"
        Nstar, Deltastar, Pstar = transfer_factor_at_pole(omega, p)
        RQ = 27.0 * p.cs**5 / (p.a**5 * omega**5 * Nstar)
        out.append(
            {
                "omega": omega,
                "family": family,
                "Nstar": Nstar,
                "Delta": Deltastar,
                "Pstar": Pstar,
                "RQ": RQ,
            }
        )
    return out


def verify_sample_compatibility_branch() -> None:
    banner("PART III — CONCRETE COMPATIBILITY BRANCH AND POLE CENSUS")

    p = PrimitiveParams(
        lamB=0.5,
        lamU=0.3,
        lamW=0.4,
        lamR=0.25,
        OmegaU=1.0,
        OmegaW=1.4,
        varpi=2.0,
        M=1.0,
        a=1.0,
        cs=1.0,
    )
    data = primitive_compatibility_data(p)
    poles = pole_census(p, data.K_compat)

    print("Primitive sample branch parameters:")
    print(p)
    print("\nCompatibility data:")
    for name in [
        "C", "GU", "GW", "R", "Delta", "Q", "H", "P",
        "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0",
        "P0_target_compat", "K_compat", "D0_compat",
    ]:
        print(f"{name} = {getattr(data, name)}")

    print("\nPole census on the compatibility branch:")
    for row in poles:
        print(row)

    wall_like = [row for row in poles if row["family"] == "wall-like"]
    if len(wall_like) != 2:
        raise AssertionError("Expected exactly two wall-like poles on the sample compatibility branch.")


# ---------------------------------------------------------------------------
# Part IV. Dynamic survival window along the compatibility surface
# ---------------------------------------------------------------------------

def stage005_required_ratio(delta_v_req: float, eta: float, x: float = 1.0) -> float:
    return 2.0 * delta_v_req * (1.0 + eta**2) * x**6 / eta


def compatibility_scan_row(lamW: float) -> dict[str, float]:
    p = PrimitiveParams(
        lamB=0.5,
        lamU=0.3,
        lamW=lamW,
        lamR=0.25,
        OmegaU=1.0,
        OmegaW=1.4,
        varpi=2.0,
        M=1.0,
        a=1.0,
        cs=1.0,
    )
    d = primitive_compatibility_data(p)
    poles = pole_census(p, d.K_compat)
    walls = [row for row in poles if row["family"] == "wall-like"]
    walls = sorted(walls, key=lambda z: z["omega"])
    return {
        "lamW": lamW,
        "P0_target_compat": d.P0_target_compat,
        "K_compat": d.K_compat,
        "omega_wall_lo": walls[0]["omega"],
        "RQ_wall_lo": walls[0]["RQ"],
        "omega_wall_hi": walls[1]["omega"],
        "RQ_wall_hi": walls[1]["RQ"],
    }


def bisect_root(func: Callable[[float], float], a: float, b: float, steps: int = 70) -> float:
    fa = func(a)
    fb = func(b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if fa * fb > 0:
        raise ValueError("Bisection interval does not bracket a root.")
    lo, hi = a, b
    flo, fhi = fa, fb
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        fm = func(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def verify_dynamic_window() -> None:
    banner("PART IV — DYNAMIC SURVIVAL WINDOW ALONG THE 5PN-COMPATIBLE SURFACE")

    # Same illustrative barrier benchmark carried from Stage 005.
    V_known_1 = 1.181909222592
    eps_energy = 0.1
    delta_v_req = V_known_1 - eps_energy

    req_10 = stage005_required_ratio(delta_v_req, eta=0.1, x=1.0)
    req_30 = stage005_required_ratio(delta_v_req, eta=0.3, x=1.0)

    print("Carried local barrier benchmark:")
    print(f"V_known(1) = {V_known_1}")
    print(f"epsilon = {eps_energy}")
    print(f"DeltaV_req(1) = {delta_v_req}")
    print(f"Required R_Q,* at eta = 0.1 = {req_10}")
    print(f"Required R_Q,* at eta = 0.3 = {req_30}")

    scan_lamW = [0.2, 0.4, 0.6, 0.8, 1.0]
    rows = [compatibility_scan_row(x) for x in scan_lamW]

    print("\nCompatibility scan in lambda_W:")
    for row in rows:
        print(row)

    if not monotone_increasing([row["P0_target_compat"] for row in rows]):
        raise AssertionError("P0_target,compat should increase monotonically on the chosen scan.")
    if not monotone_decreasing([row["K_compat"] for row in rows]):
        raise AssertionError("K_compat should decrease monotonically on the chosen scan.")
    if not monotone_decreasing([row["RQ_wall_lo"] for row in rows]):
        raise AssertionError("Lower-wall R_Q should decrease monotonically on the chosen scan.")
    if not monotone_decreasing([row["RQ_wall_hi"] for row in rows]):
        raise AssertionError("Upper-wall R_Q should decrease monotonically on the chosen scan.")

    # Threshold crossings for the lower and upper wall poles.
    def lower_minus_req10(x: float) -> float:
        return compatibility_scan_row(x)["RQ_wall_lo"] - req_10

    def upper_minus_req10(x: float) -> float:
        return compatibility_scan_row(x)["RQ_wall_hi"] - req_10

    def lower_minus_req30(x: float) -> float:
        return compatibility_scan_row(x)["RQ_wall_lo"] - req_30

    def upper_minus_req30(x: float) -> float:
        return compatibility_scan_row(x)["RQ_wall_hi"] - req_30

    lam_lo_10 = bisect_root(lower_minus_req10, 0.4, 0.5)
    lam_hi_10 = bisect_root(upper_minus_req10, 0.5, 0.6)
    lam_lo_30 = bisect_root(lower_minus_req30, 0.75, 0.76)
    lam_hi_30 = bisect_root(upper_minus_req30, 0.89, 0.90)

    row_lo_10 = compatibility_scan_row(lam_lo_10)
    row_hi_10 = compatibility_scan_row(lam_hi_10)
    row_lo_30 = compatibility_scan_row(lam_lo_30)
    row_hi_30 = compatibility_scan_row(lam_hi_30)

    print("\n10%-loss threshold crossings:")
    print("lower wall:", row_lo_10)
    print("upper wall:", row_hi_10)
    print("\n30%-loss threshold crossings:")
    print("lower wall:", row_lo_30)
    print("upper wall:", row_hi_30)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    verify_exact_target_surface()
    verify_exact_primitive_specialization()
    verify_sample_compatibility_branch()
    verify_dynamic_window()

    banner("STAGE 006 VERDICT")
    print("1. The primitive same-charge branch can be forced onto the exact isotropic 5PN surface.")
    print("2. On that branch, simultaneous one-pole and normalization success reduce to one compatibility equation.")
    print("3. The dynamic corridor is not killed automatically by this 5PN calibration.")
    print("4. But it survives only inside a finite branch-compatible target window on the explicit primitive family.")


if __name__ == "__main__":
    main()
