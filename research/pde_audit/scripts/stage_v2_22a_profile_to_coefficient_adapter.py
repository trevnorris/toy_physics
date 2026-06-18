#!/usr/bin/env python3
"""
Stage V2-22A — Profile-to-coefficient adapter for the moving-throat PDE program.

This script turns radial/axial profile data on an open finite throat into the
V2-21 branch-extraction manifest format.  It is intended to sit between a PDE or
numerical profile solver and the V2-21 observable-packet fixture.

Main jobs:
  * compute overlap integrals I_eta,phi, I_eta,u, I_eta,w, I_u,w;
  * convert overlaps into BdG and Maxwell/mixed reduced couplings;
  * generate grouped real-P2 lane manifests for 20/21/22;
  * optionally apply a weak-axisymmetric lane signature (1, 1/2, -1);
  * run the same D/N/Z/N/P residual extraction used in V2-21;
  * write a V2-21-compatible branch manifest and observable packet.

Run:
  python stage_v2_22a_profile_to_coefficient_adapter.py \
      --write-profile-manifest stage_v2_22a_profile_input_manifest.json \
      --out-v21-manifest stage_v2_22a_generated_v21_manifest.json \
      --out-json stage_v2_22a_observable_packet.json

The built-in manifest contains three demonstration branches:
  1. analytic isotropic D/N-profile branch;
  2. sampled-grid copy of the same branch;
  3. weak-axisymmetric profile/parameter slope branch.

The built-in demo is not tuned to the GR normalization target; it is an adapter
and extraction fixture, not a branch-realization proof.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import sympy as sp

LANES = ("20", "21", "22")
AXISYM_LAMBDA = {"20": 1.0, "21": 0.5, "22": -1.0}
GROUP_WEIGHTS = {"20": 1.0, "21": 2.0, "22": 2.0}


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def near_zero(x: float, tol: float = 1e-9) -> bool:
    return abs(float(x)) <= tol


def grouped_decomposition(values: Mapping[str, float]) -> Dict[str, float]:
    x20 = float(values["20"])
    x21 = float(values["21"])
    x22 = float(values["22"])
    xbar = (x20 + 2.0 * x21 + 2.0 * x22) / 5.0
    a_x = (2.0 * x20 - x21 - x22) / 10.0
    b_x = (x21 - x22) / 2.0
    return {
        "bar": xbar,
        "a": a_x,
        "b": b_x,
        "anisotropy_norm_sq": 4.0 * a_x * a_x + (4.0 / 5.0) * b_x * b_x,
        "axisymmetric_b_minus_3a": b_x - 3.0 * a_x,
    }


def lane_values(packet_by_lane: Mapping[str, Mapping[str, float]], key: str) -> Dict[str, float]:
    return {lane: float(packet_by_lane[lane][key]) for lane in LANES}


# ---------------------------------------------------------------------------
# Profile representation and integration
# ---------------------------------------------------------------------------


def builtin_expr(name: str, s: sp.Symbol, L: sp.Symbol) -> sp.Expr:
    """Return a SymPy expression for a named finite-throat profile."""
    if name in {"chi_eta_sin", "wall_sin", "u_wall"}:
        return sp.sqrt(2 / L) * sp.sin(sp.pi * s / L)
    if name in {"phi_DN", "w_DN", "half_wave_DN"}:
        return sp.sqrt(2 / L) * sp.sin(sp.pi * s / (2 * L))
    if name in {"phi_ND", "cos_half"}:
        return sp.sqrt(2 / L) * sp.cos(sp.pi * s / (2 * L))
    if name in {"constant_unit", "one"}:
        return sp.Integer(1)
    raise ValueError(f"unknown built-in profile name: {name}")


def analytic_expr(profile: Mapping[str, Any], s: sp.Symbol, L: sp.Symbol) -> Optional[sp.Expr]:
    """Return an analytic expression for a profile spec, or None if sampled."""
    kind = profile.get("kind", "builtin")
    if kind == "builtin":
        return builtin_expr(str(profile["name"]), s, L)
    if kind == "expr":
        # Local scientific convenience.  This script is intended for trusted local manifests.
        namespace = {"s": s, "L": L, "pi": sp.pi, "sqrt": sp.sqrt, "sin": sp.sin, "cos": sp.cos, "exp": sp.exp}
        return sp.sympify(profile["expr"], locals=namespace)
    if kind == "sampled":
        return None
    raise ValueError(f"unknown profile kind: {kind}")


def make_numeric_evaluator(profile: Mapping[str, Any], L_value: float) -> Callable[[float], float]:
    """Return f(s) for analytic or sampled profile specs."""
    kind = profile.get("kind", "builtin")
    if kind in {"builtin", "expr"}:
        s_sym, L_sym = sp.symbols("s L", positive=True)
        expr = analytic_expr(profile, s_sym, L_sym)
        assert expr is not None
        f = sp.lambdify(s_sym, expr.subs(L_sym, L_value), "math")
        return lambda x: float(f(float(x)))
    if kind == "sampled":
        samples = sorted((float(p[0]), float(p[1])) for p in profile["samples"])
        if len(samples) < 2:
            raise ValueError("sampled profile requires at least two points")

        def interp(x: float) -> float:
            x = float(x)
            if x <= samples[0][0]:
                return samples[0][1]
            if x >= samples[-1][0]:
                return samples[-1][1]
            lo, hi = 0, len(samples) - 1
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if samples[mid][0] <= x:
                    lo = mid
                else:
                    hi = mid
            x0, y0 = samples[lo]
            x1, y1 = samples[hi]
            t = (x - x0) / (x1 - x0)
            return y0 * (1.0 - t) + y1 * t

        return interp
    raise ValueError(f"unknown profile kind: {kind}")


def integrate_product(
    profiles: Sequence[Mapping[str, Any]],
    weight: Mapping[str, Any],
    L_value: float,
    grid_points: int = 4001,
) -> Dict[str, Any]:
    """Integrate weight * product(profiles) over [0, L].

    Uses exact SymPy integration if all profiles and the weight are analytic;
    otherwise uses a composite trapezoidal rule.
    """
    s_sym, L_sym = sp.symbols("s L", positive=True)
    exprs: List[sp.Expr] = []
    analytic = True
    for spec in [weight, *profiles]:
        expr = analytic_expr(spec, s_sym, L_sym)
        if expr is None:
            analytic = False
            break
        exprs.append(expr)

    if analytic:
        integrand = sp.prod(exprs)
        exact = sp.integrate(integrand, (s_sym, 0, L_sym))
        exact_at_L = sp.simplify(exact.subs(L_sym, sp.nsimplify(L_value)))
        return {
            "method": "analytic_sympy",
            "exact": str(sp.simplify(exact)),
            "value": float(sp.N(exact_at_L, 30)),
        }

    evaluators = [make_numeric_evaluator(spec, L_value) for spec in [weight, *profiles]]
    n = int(grid_points)
    if n < 2:
        raise ValueError("grid_points must be >= 2")
    h = L_value / (n - 1)
    total = 0.0
    for i in range(n):
        x = i * h
        y = 1.0
        for ev in evaluators:
            y *= ev(x)
        total += (0.5 if i == 0 or i == n - 1 else 1.0) * y
    return {"method": "trapezoid", "grid_points": n, "value": total * h}


def sample_builtin_profile(name: str, L_value: float, points: int) -> Dict[str, Any]:
    spec = {"kind": "builtin", "name": name}
    f = make_numeric_evaluator(spec, L_value)
    samples = []
    for i in range(points):
        s = L_value * i / (points - 1)
        samples.append([s, f(s)])
    return {"kind": "sampled", "samples": samples}


# ---------------------------------------------------------------------------
# Symbolic audit
# ---------------------------------------------------------------------------


def run_symbolic_audit() -> Dict[str, Any]:
    checks: List[Tuple[str, bool, str]] = []
    s, L = sp.symbols("s L", positive=True)
    chi = builtin_expr("chi_eta_sin", s, L)
    phi = builtin_expr("phi_DN", s, L)
    one = sp.Integer(1)

    I_chi_chi = sp.simplify(sp.integrate(chi * chi, (s, 0, L)))
    I_chi_phi = sp.simplify(sp.integrate(chi * phi, (s, 0, L)))
    I_phi_phi = sp.simplify(sp.integrate(phi * phi, (s, 0, L)))
    checks.append(("I_chi_chi_is_1", I_chi_chi == 1, str(I_chi_chi)))
    checks.append(("I_phi_phi_is_1", I_phi_phi == 1, str(I_phi_phi)))
    checks.append(("I_chi_phi_is_8_over_3pi", sp.simplify(I_chi_phi - 8 / (3 * sp.pi)) == 0, str(I_chi_phi)))

    lam20, lam21, lam22 = sp.Integer(1), sp.Rational(1, 2), sp.Integer(-1)
    x0, eps, x1 = sp.symbols("x0 eps x1")
    y20, y21, y22 = x0 + eps * lam20 * x1, x0 + eps * lam21 * x1, x0 + eps * lam22 * x1
    abar = sp.simplify((y20 + 2 * y21 + 2 * y22) / 5)
    a = sp.simplify((2 * y20 - y21 - y22) / 10)
    b = sp.simplify((y21 - y22) / 2)
    checks.append(("axisym_trace_unchanged", sp.simplify(abar - x0) == 0, str(abar)))
    checks.append(("axisym_b_equals_3a", sp.simplify(b - 3 * a) == 0, f"a={a}, b={b}"))

    # Verify extraction formulas used downstream.
    D0, D2, D4, N0, N2, N4, w = sp.symbols("D0 D2 D4 N0 N2 N4 w", nonzero=True)
    D = D0 + D2 * w**2 + D4 * w**4
    N = N0 + N2 * w**2 + N4 * w**4
    Y = sp.series(D0 / D, w, 0, 6).removeO().expand()
    Pref = sp.series(D0 * N / D**2, w, 0, 6).removeO().expand()
    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2
    P0 = N0 / D0
    P2 = (D0 * N2 - 2 * D2 * N0) / D0**2
    P4 = (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3
    checks.append(("u2_formula", sp.simplify(Y.coeff(w, 2) - u2) == 0, str(Y.coeff(w, 2))))
    checks.append(("u4_formula", sp.simplify(Y.coeff(w, 4) - u4) == 0, str(Y.coeff(w, 4))))
    checks.append(("P0_formula", sp.simplify(Pref.coeff(w, 0) - P0) == 0, str(Pref.coeff(w, 0))))
    checks.append(("P2_formula", sp.simplify(Pref.coeff(w, 2) - P2) == 0, str(Pref.coeff(w, 2))))
    checks.append(("P4_formula", sp.simplify(Pref.coeff(w, 4) - P4) == 0, str(Pref.coeff(w, 4))))

    passed = sum(1 for _, ok, _ in checks if ok)
    return {
        "checks_total": len(checks),
        "checks_passed": passed,
        "checks_failed": len(checks) - passed,
        "details": [{"name": name, "pass": bool(ok), "note": note} for name, ok, note in checks],
    }


# ---------------------------------------------------------------------------
# Profile manifest -> V2-21 branch manifest
# ---------------------------------------------------------------------------


def default_profile_manifest() -> Dict[str, Any]:
    L_value = 1.0
    target = {"constants": {"G": 1.0, "c_s": 1.0, "c": 1.0, "a": 1.0, "mhat0": 1.0, "S_port": 1.0, "theta_tail": 1.0}}
    base_profiles = {
        "weight": {"kind": "builtin", "name": "one"},
        "wall": {"kind": "builtin", "name": "chi_eta_sin"},
        "bdg_phi": {"kind": "builtin", "name": "phi_DN"},
        "u": {"kind": "builtin", "name": "chi_eta_sin"},
        "w": {"kind": "builtin", "name": "phi_DN"},
    }
    sampled_profiles = {
        "weight": {"kind": "builtin", "name": "one"},
        "wall": sample_builtin_profile("chi_eta_sin", L_value, 801),
        "bdg_phi": sample_builtin_profile("phi_DN", L_value, 801),
        "u": sample_builtin_profile("chi_eta_sin", L_value, 801),
        "w": sample_builtin_profile("phi_DN", L_value, 801),
    }
    common_reduction = {
        "wall": {"K": 1.6, "M": 0.9},
        "bdg_modes": [{"name": "bdg_halfwave", "lambda_B": 0.42, "varpi": 2.7, "profile": "bdg_phi"}],
        "mixed_ports": [{"name": "one_mixed_port", "lambda_U": 0.28, "Omega_U": 3.0, "u_profile": "u", "lambda_W": 0.52, "Omega_W": 4.0, "w_profile": "w", "lambda_R": 0.65, "u_w_profiles": ["u", "w"]}],
    }
    weak_axisym = {
        "epsilon": 0.01,
        "slopes": {
            # These slopes are deliberately small and sector-mixed; the purpose is
            # to test b=3a transport, not to tune the target residuals.
            "K": 0.10,
            "M": -0.04,
            "lambda_B": 0.03,
            "lambda_U": -0.02,
            "lambda_W": 0.05,
            "lambda_R": 0.02,
            "Omega_U": 0.01,
            "Omega_W": -0.015,
        },
    }
    return {
        "schema": "stage_v2_22a_profile_adapter/v1",
        "description": "Built-in analytic/sampled profile-to-coefficient demonstration manifest.",
        "branches": [
            {
                "name": "analytic_isotropic_DN_profile_demo",
                "geometry": {"R_exit": 0.35, "boundary_class": "open_impedance", "Y_L_limit": 0.0, "L": L_value},
                "target": target,
                "profiles": base_profiles,
                "reduction": common_reduction,
            },
            {
                "name": "sampled_isotropic_DN_profile_demo",
                "geometry": {"R_exit": 0.35, "boundary_class": "open_impedance", "Y_L_limit": 0.0, "L": L_value},
                "target": target,
                "profiles": sampled_profiles,
                "integration": {"grid_points": 4001},
                "reduction": common_reduction,
            },
            {
                "name": "weak_axisymmetric_profile_slope_demo",
                "geometry": {"R_exit": 0.35, "boundary_class": "open_impedance", "Y_L_limit": 0.0, "L": L_value},
                "target": target,
                "profiles": base_profiles,
                "reduction": common_reduction,
                "weak_axisymmetric": weak_axisym,
            },
        ],
    }


def profile_lookup(branch: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    profiles = branch["profiles"]
    if name not in profiles:
        raise KeyError(f"profile {name!r} is not defined in branch {branch.get('name', 'unnamed')}")
    return profiles[name]


def lane_scaled(base: float, branch: Mapping[str, Any], lane: str, key: str) -> float:
    """Apply optional weak-axisymmetric scaling to a base primitive parameter."""
    weak = branch.get("weak_axisymmetric")
    if not weak:
        return float(base)
    eps = float(weak.get("epsilon", 0.0))
    slope = float(weak.get("slopes", {}).get(key, 0.0))
    return float(base) * (1.0 + eps * AXISYM_LAMBDA[lane] * slope)


def compute_branch_overlaps(branch: Mapping[str, Any]) -> Dict[str, Any]:
    L_value = float(branch.get("geometry", {}).get("L", 1.0))
    grid_points = int(branch.get("integration", {}).get("grid_points", 4001))
    weight = branch["profiles"].get("weight", {"kind": "builtin", "name": "one"})
    wall_profile = profile_lookup(branch, "wall")

    overlaps: Dict[str, Any] = {}
    for mode in branch["reduction"].get("bdg_modes", []):
        phi = profile_lookup(branch, mode["profile"])
        overlaps[f"I_eta_phi::{mode['name']}"] = integrate_product([wall_profile, phi], weight, L_value, grid_points)
    for port in branch["reduction"].get("mixed_ports", []):
        u = profile_lookup(branch, port["u_profile"])
        w = profile_lookup(branch, port["w_profile"])
        overlaps[f"I_eta_u::{port['name']}"] = integrate_product([wall_profile, u], weight, L_value, grid_points)
        overlaps[f"I_eta_w::{port['name']}"] = integrate_product([wall_profile, w], weight, L_value, grid_points)
        # R uses the u/w overlap by default.
        uw_names = port.get("u_w_profiles", [port["u_profile"], port["w_profile"]])
        uw_0 = profile_lookup(branch, uw_names[0])
        uw_1 = profile_lookup(branch, uw_names[1])
        overlaps[f"I_u_w::{port['name']}"] = integrate_product([uw_0, uw_1], weight, L_value, grid_points)
    return overlaps


def adapt_profile_branch_to_v21(branch: Mapping[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Convert one profile branch into a V2-21 branch manifest plus diagnostics."""
    reduction = branch["reduction"]
    if "direct_coefficients" in reduction:
        direct_coefficients = {k: float(v) for k, v in reduction["direct_coefficients"].items()}
        v21_branch = {
            "name": branch["name"],
            "geometry": {k: v for k, v in branch.get("geometry", {}).items() if k != "L"},
            "target": branch.get("target", {}),
            "lanes": {lane: {"direct_coefficients": copy.deepcopy(direct_coefficients)} for lane in LANES},
            "profile_adapter": {
                "source_branch_hash": stable_hash(branch),
                "coefficient_path": "direct_derived_coefficients",
                "derived_coefficients_provenance": reduction.get("derived_coefficients_provenance", {}),
                "weak_axisymmetric": branch.get("weak_axisymmetric", None),
            },
        }
        diagnostics = {
            "source_branch": branch["name"],
            "source_branch_hash": stable_hash(branch),
            "coefficient_path": "direct_derived_coefficients",
            "direct_coefficients": direct_coefficients,
            "v21_branch_hash": stable_hash(v21_branch),
        }
        return v21_branch, diagnostics

    overlaps = compute_branch_overlaps(branch)
    lanes: Dict[str, Any] = {}
    for lane in LANES:
        wall = reduction["wall"]
        lane_packet = {
            "K": lane_scaled(float(wall["K"]), branch, lane, "K"),
            "M": lane_scaled(float(wall["M"]), branch, lane, "M"),
            "bdg_modes": [],
            "mixed_ports": [],
        }
        for mode in reduction.get("bdg_modes", []):
            I = float(overlaps[f"I_eta_phi::{mode['name']}"]["value"])
            lam_B = lane_scaled(float(mode["lambda_B"]), branch, lane, "lambda_B")
            varpi = lane_scaled(float(mode["varpi"]), branch, lane, "varpi")
            lane_packet["bdg_modes"].append({
                "name": mode["name"],
                "coupling": lam_B * I,
                "varpi": varpi,
                "overlap_I_eta_phi": I,
                "lambda_B_eff": lam_B,
            })
        for port in reduction.get("mixed_ports", []):
            I_eta_u = float(overlaps[f"I_eta_u::{port['name']}"]["value"])
            I_eta_w = float(overlaps[f"I_eta_w::{port['name']}"]["value"])
            I_u_w = float(overlaps[f"I_u_w::{port['name']}"]["value"])
            lam_U = lane_scaled(float(port["lambda_U"]), branch, lane, "lambda_U")
            lam_W = lane_scaled(float(port["lambda_W"]), branch, lane, "lambda_W")
            lam_R = lane_scaled(float(port["lambda_R"]), branch, lane, "lambda_R")
            OU = lane_scaled(float(port["Omega_U"]), branch, lane, "Omega_U")
            OW = lane_scaled(float(port["Omega_W"]), branch, lane, "Omega_W")
            lane_packet["mixed_ports"].append({
                "name": port["name"],
                "Omega_U": OU,
                "Omega_W": OW,
                "R": lam_R * I_u_w,
                "g_U": lam_U * I_eta_u,
                "g_W": lam_W * I_eta_w,
                "overlap_I_eta_u": I_eta_u,
                "overlap_I_eta_w": I_eta_w,
                "overlap_I_u_w": I_u_w,
                "lambda_U_eff": lam_U,
                "lambda_W_eff": lam_W,
                "lambda_R_eff": lam_R,
            })
        lanes[lane] = lane_packet

    v21_branch = {
        "name": branch["name"],
        "geometry": {k: v for k, v in branch.get("geometry", {}).items() if k != "L"},
        "target": branch.get("target", {}),
        "lanes": lanes,
        "profile_adapter": {
            "source_branch_hash": stable_hash(branch),
            "overlap_summary": overlaps,
            "weak_axisymmetric": branch.get("weak_axisymmetric", None),
        },
    }
    diagnostics = {
        "source_branch": branch["name"],
        "source_branch_hash": stable_hash(branch),
        "overlaps": overlaps,
        "v21_branch_hash": stable_hash(v21_branch),
    }
    return v21_branch, diagnostics


# ---------------------------------------------------------------------------
# Extraction formulas copied intentionally from V2-21, kept local for robustness
# ---------------------------------------------------------------------------


def bdg_moments(modes: Iterable[Mapping[str, float]]) -> Dict[str, float]:
    B0 = B2 = B4 = 0.0
    for mode in modes:
        coupling = float(mode["coupling"])
        varpi = float(mode["varpi"])
        if varpi <= 0:
            raise ValueError(f"BdG mode has nonpositive varpi={varpi}")
        B0 += coupling**2 / varpi**2
        B2 += coupling**2 / varpi**4
        B4 += coupling**2 / varpi**6
    return {"B0": B0, "B2": B2, "B4": B4}


def mixed_port_moments(ports: Iterable[Mapping[str, float]]) -> Tuple[Dict[str, float], List[Dict[str, float]]]:
    Z0 = Z2 = Z4 = N0 = N2 = N4 = 0.0
    diagnostics: List[Dict[str, float]] = []
    for idx, port in enumerate(ports):
        OU = float(port["Omega_U"])
        OW = float(port["Omega_W"])
        R = float(port.get("R", 0.0))
        gU = float(port.get("g_U", 0.0))
        gW = float(port.get("g_W", 0.0))
        Delta = OU**2 * OW**2 - R**2
        if Delta == 0.0:
            raise ValueError(f"port {idx} has Delta=0")
        S = OU**2 + OW**2
        Q = gU**2 * OW**2 + 2.0 * gU * gW * R + gW**2 * OU**2
        H = gU**2 + gW**2
        P = OU**2 * gW + R * gU
        Z0_r = Q / Delta
        Z2_r = (Q * S - H * Delta) / Delta**2
        Z4_r = (Q * (S**2 - Delta) - S * H * Delta) / Delta**3
        N0_r = P**2 / Delta**2
        N2_r = 2.0 * P * (P * S - Delta * gW) / Delta**3
        N4_r = (Delta**2 * gW**2 - 2.0 * Delta * P**2 - 4.0 * Delta * P * S * gW + 3.0 * P**2 * S**2) / Delta**4
        Z0 += Z0_r; Z2 += Z2_r; Z4 += Z4_r
        N0 += N0_r; N2 += N2_r; N4 += N4_r
        diagnostics.append({
            "port_index": idx,
            "Delta": Delta,
            "S": S,
            "Q": Q,
            "H": H,
            "P": P,
            "Z0": Z0_r,
            "Z2": Z2_r,
            "Z4": Z4_r,
            "N0": N0_r,
            "N2": N2_r,
            "N4": N4_r,
        })
    return {"Z0": Z0, "Z2": Z2, "Z4": Z4, "N0": N0, "N2": N2, "N4": N4}, diagnostics


def lane_extract(lane: Mapping[str, Any]) -> Dict[str, Any]:
    if "direct_coefficients" in lane:
        coeff = {k: float(v) for k, v in lane["direct_coefficients"].items()}
        required = {"K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4"}
        missing = required - set(coeff)
        if missing:
            raise ValueError(f"direct_coefficients missing {sorted(missing)}")
        port_diag: List[Dict[str, float]] = []
    else:
        K = float(lane["K"])
        M = float(lane["M"])
        b = bdg_moments(lane.get("bdg_modes", []))
        z_n, port_diag = mixed_port_moments(lane.get("mixed_ports", []))
        coeff = {"K": K, "M": M, **b, **z_n}
    K = coeff["K"]
    M = coeff["M"]
    D0 = coeff["K"] - coeff["B0"] - coeff["Z0"]
    D2 = -(coeff["M"] + coeff["B2"] + coeff["Z2"])
    D4 = -(coeff["B4"] + coeff["Z4"])
    N0, N2, N4 = coeff["N0"], coeff["N2"], coeff["N4"]
    u2 = -D2 / D0 if D0 != 0 else math.nan
    u4 = (D2**2 - D0 * D4) / D0**2 if D0 != 0 else math.nan
    P0 = N0 / D0 if D0 != 0 else math.nan
    P2 = (D0 * N2 - 2.0 * D2 * N0) / D0**2 if D0 != 0 else math.nan
    P4 = (D0**2 * N4 - 2.0 * D0 * (D2 * N2 + D4 * N0) + 3.0 * D2**2 * N0) / D0**3 if D0 != 0 else math.nan
    A = coeff["M"] + coeff["B2"] + coeff["Z2"]
    C = coeff["B4"] + coeff["Z4"]
    R_pole = D0 * C - 3.0 * A**2
    stability = {
        "D0_positive": D0 > 0.0,
        "C_positive": C > 0.0,
        "N0_nonzero": N0 != 0.0,
        "wall_M_positive": M > 0.0,
        "all_Delta_positive": all(p["Delta"] > 0.0 for p in port_diag),
    }
    return {**coeff, "D0": D0, "D2": D2, "D4": D4, "u2": u2, "u4": u4, "P0": P0, "P2": P2, "P4": P4, "A_MplusB2plusZ2": A, "C_B4plusZ4": C, "R_pole": R_pole, "port_diagnostics": port_diag, "stability": stability}


def extract_branch(branch: Mapping[str, Any], tol: float = 1e-9) -> Dict[str, Any]:
    lanes_in = branch["lanes"]
    lane_packets = {lane: lane_extract(lanes_in[lane]) for lane in LANES}
    grouped: Dict[str, Dict[str, float]] = {}
    for key in ["K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4", "D0", "D2", "D4", "u2", "u4", "P0", "P2", "P4", "R_pole"]:
        grouped[key] = grouped_decomposition(lane_values(lane_packets, key))

    constants = branch.get("target", {}).get("constants", {})
    G = float(constants.get("G", 1.0)); c_s = float(constants.get("c_s", 1.0)); c = float(constants.get("c", 1.0)); a = float(constants.get("a", 1.0))
    mhat0 = float(constants.get("mhat0", 1.0)); S_port = float(constants.get("S_port", 1.0)); theta_tail = float(constants.get("theta_tail", 1.0))
    D0_bar = grouped["D0"]["bar"]
    A_bar = grouped["M"]["bar"] + grouped["B2"]["bar"] + grouped["Z2"]["bar"]
    C_bar = grouped["B4"]["bar"] + grouped["Z4"]["bar"]
    P0_bar = grouped["P0"]["bar"]
    P2_bar = grouped["P2"]["bar"]
    P4_bar = grouped["P4"]["bar"]
    P0_target = 54.0 * G * c_s**5 / (5.0 * a**5 * c**5)
    gamma_eff = mhat0**2 * S_port * P0_bar * a**5 / (27.0 * c_s**5)
    gamma_GR = 2.0 * G / (5.0 * c**5)
    R_pole = D0_bar * C_bar - 3.0 * A_bar**2
    R_norm = mhat0**2 * S_port * P0_bar - P0_target
    R_tail = theta_tail * (c / c_s) ** 3 - 1.0
    geometry = branch.get("geometry", {})
    open_gate = {
        "R_exit_positive": float(geometry.get("R_exit", 0.0)) > 0.0,
        "boundary_class_open_impedance": geometry.get("boundary_class") == "open_impedance",
        "hard_cap_forbidden": geometry.get("boundary_class") != "hard_cap" and float(geometry.get("R_exit", 0.0)) > 0.0,
    }
    lane_stability_all = {
        name: all(packet["stability"].get(name, False) for packet in lane_packets.values())
        for name in sorted({k for packet in lane_packets.values() for k in packet["stability"]})
    }
    residuals = {"R_pole": R_pole, "R_norm": R_norm, "R_P2": P2_bar, "R_P4": P4_bar, "R_tail": R_tail, "gamma_eff_minus_gamma_GR": gamma_eff - gamma_GR}
    flags = {
        "open_gate_pass": all(open_gate.values()),
        "stability_gate_pass": all(lane_stability_all.values()),
        "isotropic_D0_pass": near_zero(grouped["D0"]["anisotropy_norm_sq"], tol),
        "isotropic_N0_pass": near_zero(grouped["N0"]["anisotropy_norm_sq"], tol),
        "axisymmetric_P0_pattern_pass": near_zero(grouped["P0"]["axisymmetric_b_minus_3a"], tol),
        "one_pole_pass": near_zero(R_pole, tol),
        "normalization_pass": near_zero(R_norm, tol),
        "constant_prefactor_P2_pass": near_zero(P2_bar, tol),
        "constant_prefactor_P4_pass": near_zero(P4_bar, tol),
        "tail_transport_pass": near_zero(R_tail, tol),
    }
    flags["target_packet_pass"] = all(flags[k] for k in ["open_gate_pass", "stability_gate_pass", "isotropic_D0_pass", "isotropic_N0_pass", "one_pole_pass", "normalization_pass", "constant_prefactor_P2_pass", "constant_prefactor_P4_pass", "tail_transport_pass"])
    return {
        "name": branch["name"],
        "input_hash": stable_hash({"name": branch["name"], "geometry": branch.get("geometry"), "target": branch.get("target"), "lanes": branch.get("lanes")}),
        "open_gate": open_gate,
        "lane_packets": lane_packets,
        "grouped": grouped,
        "target_values": {"P0_target": P0_target, "gamma_eff": gamma_eff, "gamma_GR": gamma_GR},
        "residuals": residuals,
        "stability": lane_stability_all,
        "pass_flags": flags,
    }


def adapt_manifest(profile_manifest: Mapping[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    v21_branches: List[Dict[str, Any]] = []
    diagnostics: List[Dict[str, Any]] = []
    for branch in profile_manifest.get("branches", []):
        v21_branch, diag = adapt_profile_branch_to_v21(branch)
        v21_branches.append(v21_branch)
        diagnostics.append(diag)
    v21_manifest = {
        "schema": "stage_v2_21_branch_extraction_fixture/v1",
        "description": "Generated by stage_v2_22a_profile_to_coefficient_adapter.py from profile overlap data.",
        "source_profile_manifest_hash": stable_hash(profile_manifest),
        "branches": v21_branches,
    }
    adapter_diag = {
        "profile_manifest_hash": stable_hash(profile_manifest),
        "v21_manifest_hash": stable_hash(v21_manifest),
        "branch_diagnostics": diagnostics,
    }
    return v21_manifest, adapter_diag


def load_manifest(path: Optional[str]) -> Dict[str, Any]:
    if path is None:
        return default_profile_manifest()
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def validate_profile_manifest(profile_manifest: Mapping[str, Any], strict_profiles: bool = False) -> Dict[str, Any]:
    """Validate profile-manifest features that affect trusted-local parsing."""
    issues: List[Dict[str, str]] = []
    allowed = {"builtin", "sampled"} if strict_profiles else {"builtin", "sampled", "expr"}

    for bidx, branch in enumerate(profile_manifest.get("branches", [])):
        profiles = branch.get("profiles", {})
        if not isinstance(profiles, Mapping):
            issues.append({
                "severity": "error",
                "path": f"branches[{bidx}].profiles",
                "message": "profiles must be an object",
            })
            continue
        for pname, profile in profiles.items():
            if not isinstance(profile, Mapping):
                issues.append({
                    "severity": "error",
                    "path": f"branches[{bidx}].profiles.{pname}",
                    "message": "profile must be an object",
                })
                continue
            kind = str(profile.get("kind", "builtin"))
            if kind not in allowed:
                if strict_profiles and kind == "expr":
                    message = "expr profiles are disabled in strict profile mode"
                else:
                    message = f"unsupported profile kind {kind!r}"
                issues.append({
                    "severity": "error",
                    "path": f"branches[{bidx}].profiles.{pname}.kind",
                    "message": message,
                })

    error_count = sum(1 for item in issues if item["severity"] == "error")
    return {
        "schema": "stage_v2_22a_profile_manifest_validation/v1",
        "strict_profiles": strict_profiles,
        "validation_pass": error_count == 0,
        "error_count": error_count,
        "issues": issues,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage V2-22A profile-to-coefficient adapter")
    parser.add_argument("--profile-manifest", default=None, help="Optional profile input manifest JSON")
    parser.add_argument("--write-profile-manifest", default=None, help="Write the built-in profile manifest to this path")
    parser.add_argument("--out-v21-manifest", default=None, help="Write generated V2-21-compatible branch manifest")
    parser.add_argument("--out-json", default=None, help="Write adapter observable packet JSON")
    parser.add_argument("--tol", type=float, default=1e-9, help="Pass/fail tolerance")
    parser.add_argument("--strict-profiles", action="store_true", help="Reject trusted-local expr profiles; allow only builtin and sampled profiles")
    args = parser.parse_args()

    if args.write_profile_manifest:
        with open(args.write_profile_manifest, "w", encoding="utf-8") as f:
            json.dump(default_profile_manifest(), f, indent=2, sort_keys=True)

    profile_manifest = load_manifest(args.profile_manifest)
    manifest_validation = validate_profile_manifest(profile_manifest, strict_profiles=args.strict_profiles)
    if not manifest_validation["validation_pass"]:
        print("STAGE V2-22A PROFILE-TO-COEFFICIENT ADAPTER")
        print(f"profile_manifest_validation_pass: {manifest_validation['validation_pass']}")
        print(f"profile_manifest_error_count: {manifest_validation['error_count']}")
        for issue in manifest_validation["issues"]:
            print(f"  - {issue['severity']} {issue['path']}: {issue['message']}")
        return 2

    audit = run_symbolic_audit()
    v21_manifest, adapter_diag = adapt_manifest(profile_manifest)
    branch_packets = [extract_branch(branch, tol=args.tol) for branch in v21_manifest["branches"]]
    result = {
        "script": "stage_v2_22a_profile_to_coefficient_adapter.py",
        "profile_manifest_validation": manifest_validation,
        "symbolic_audit": audit,
        "adapter_diagnostics": adapter_diag,
        "v21_manifest_hash": stable_hash(v21_manifest),
        "branch_count": len(branch_packets),
        "branches": branch_packets,
    }

    if args.out_v21_manifest:
        with open(args.out_v21_manifest, "w", encoding="utf-8") as f:
            json.dump(v21_manifest, f, indent=2, sort_keys=True)
    if args.out_json:
        with open(args.out_json, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, sort_keys=True)

    print("STAGE V2-22A PROFILE-TO-COEFFICIENT ADAPTER")
    print(f"symbolic_checks: {audit['checks_passed']}/{audit['checks_total']} passed")
    print(f"profile_manifest_hash: {adapter_diag['profile_manifest_hash']}")
    print(f"generated_v21_manifest_hash: {stable_hash(v21_manifest)}")
    for diag in adapter_diag["branch_diagnostics"]:
        print("---")
        print(f"profile_branch: {diag['source_branch']}")
        for key, payload in sorted(diag["overlaps"].items()):
            val = payload["value"]
            exact = payload.get("exact")
            if exact:
                print(f"{key}: {val:.16g}  exact={exact}")
            else:
                print(f"{key}: {val:.16g}  method={payload['method']}")
    for packet in branch_packets:
        print("---")
        print(f"branch: {packet['name']}")
        print(f"open_gate_pass: {packet['pass_flags']['open_gate_pass']}")
        print(f"stability_gate_pass: {packet['pass_flags']['stability_gate_pass']}")
        print(f"target_packet_pass: {packet['pass_flags']['target_packet_pass']}")
        print(f"D0_bar: {packet['grouped']['D0']['bar']:.16g}")
        print(f"N0_bar: {packet['grouped']['N0']['bar']:.16g}")
        print(f"P0_bar: {packet['grouped']['P0']['bar']:.16g}")
        print(f"P0_target: {packet['target_values']['P0_target']:.16g}")
        print(f"R_pole: {packet['residuals']['R_pole']:.16g}")
        print(f"R_norm: {packet['residuals']['R_norm']:.16g}")
        print(f"R_P2: {packet['residuals']['R_P2']:.16g}")
        print(f"R_P4: {packet['residuals']['R_P4']:.16g}")
        print(f"P0_anisotropy_norm_sq: {packet['grouped']['P0']['anisotropy_norm_sq']:.16g}")
        print(f"P0_axisymmetric_b_minus_3a: {packet['grouped']['P0']['axisymmetric_b_minus_3a']:.16g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
