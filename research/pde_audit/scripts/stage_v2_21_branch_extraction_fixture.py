#!/usr/bin/env python3
"""
Stage V2-21 — Moving-throat branch extraction fixture.

This script is a reusable prototype for converting a frozen open-throat branch
packet into the observable grouped-P2 extraction packet used by Volume 2:

  * open-throat / no-hard-cap gate,
  * BdG support moments B0, B2, B4,
  * conservative Maxwell/mixed moments Z0, Z2, Z4,
  * outgoing-transfer moments N0, N2, N4,
  * wall operator moments D0, D2, D4,
  * normalized response moments u2, u4,
  * outgoing prefactor moments P0, P2, P4,
  * grouped trace/anomaly decompositions,
  * stability and target residuals.

The script supports two input styles per grouped lane:

  1. primitive mode data: wall K/M, BdG modes, Maxwell/mixed ports;
  2. direct coefficient data: K/M/Bn/Zn/Nn, useful for exact algebra tests.

Run:
  python stage_v2_21_branch_extraction_fixture.py \
      --out-json stage_v2_21_branch_extraction_packet.json

or with a manifest:
  python stage_v2_21_branch_extraction_fixture.py \
      --manifest stage_v2_21_sample_branch_manifest.json \
      --out-json packet.json
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

try:
    import sympy as sp
except Exception as exc:  # pragma: no cover
    raise RuntimeError("SymPy is required for the symbolic audit section.") from exc

LANES = ("20", "21", "22")
GROUP_WEIGHTS = {"20": 1.0, "21": 2.0, "22": 2.0}


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def stable_hash(obj: Any) -> str:
    """Return a stable SHA256 hash of a JSON-serializable object."""
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def near_zero(x: float, tol: float = 1e-10) -> bool:
    return abs(float(x)) <= tol


def finite_number(x: Any) -> bool:
    try:
        y = float(x)
        return math.isfinite(y)
    except Exception:
        return False


def lane_values(packet_by_lane: Mapping[str, Mapping[str, float]], key: str) -> Dict[str, float]:
    return {lane: float(packet_by_lane[lane][key]) for lane in LANES}


def grouped_decomposition(values: Mapping[str, float]) -> Dict[str, float]:
    """Weighted grouped-P2 decomposition for values keyed by 20/21/22."""
    x20 = float(values["20"])
    x21 = float(values["21"])
    x22 = float(values["22"])
    xbar = (x20 + 2.0 * x21 + 2.0 * x22) / 5.0
    a_x = (2.0 * x20 - x21 - x22) / 10.0
    b_x = (x21 - x22) / 2.0
    anisotropy_norm_sq = 4.0 * a_x * a_x + (4.0 / 5.0) * b_x * b_x
    return {
        "bar": xbar,
        "a": a_x,
        "b": b_x,
        "anisotropy_norm_sq": anisotropy_norm_sq,
        "axisymmetric_b_minus_3a": b_x - 3.0 * a_x,
    }


# ---------------------------------------------------------------------------
# Symbolic audit
# ---------------------------------------------------------------------------


def run_symbolic_audit() -> Dict[str, Any]:
    """Verify the algebraic formulas used by the fixture."""
    checks: List[Tuple[str, bool, str]] = []

    D0, D2, D4, N0, N2, N4, w = sp.symbols("D0 D2 D4 N0 N2 N4 w", nonzero=True)
    D = D0 + D2 * w**2 + D4 * w**4
    N = N0 + N2 * w**2 + N4 * w**4

    Y = sp.series(D0 / D, w, 0, 6).removeO().expand()
    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2
    checks.append(("response_u2", sp.simplify(Y.coeff(w, 2) - u2) == 0, str(Y.coeff(w, 2))))
    checks.append(("response_u4", sp.simplify(Y.coeff(w, 4) - u4) == 0, str(Y.coeff(w, 4))))

    Pref = sp.series(D0 * N / D**2, w, 0, 6).removeO().expand()
    P0 = N0 / D0
    P2 = (D0 * N2 - 2 * D2 * N0) / D0**2
    P4 = (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3
    checks.append(("prefactor_P0", sp.simplify(Pref.coeff(w, 0) - P0) == 0, str(Pref.coeff(w, 0))))
    checks.append(("prefactor_P2", sp.simplify(Pref.coeff(w, 2) - P2) == 0, str(Pref.coeff(w, 2))))
    checks.append(("prefactor_P4", sp.simplify(Pref.coeff(w, 4) - P4) == 0, str(Pref.coeff(w, 4))))

    N2_const = 2 * D2 * N0 / D0
    N4_const = N0 * (D2**2 + 2 * D0 * D4) / D0**2
    checks.append(("constant_prefactor_P2_zero", sp.simplify(P2.subs(N2, N2_const)) == 0, "P2=0"))
    checks.append(("constant_prefactor_P4_zero", sp.simplify(P4.subs({N2: N2_const, N4: N4_const})) == 0, "P4=0"))

    K, B0, B2, B4, Z0, Z2, Z4, M = sp.symbols("K B0 B2 B4 Z0 Z2 Z4 M", nonzero=True)
    D0_iso = K - B0 - Z0
    D2_iso = -(M + B2 + Z2)
    D4_iso = -(B4 + Z4)
    u2_iso = -D2_iso / D0_iso
    u4_iso = (D2_iso**2 - D0_iso * D4_iso) / D0_iso**2
    R_pole = D0_iso * (B4 + Z4) - 3 * (M + B2 + Z2) ** 2
    checks.append(("one_pole_surface_identity", sp.simplify((u4_iso - 4 * u2_iso**2) * D0_iso**2 - R_pole) == 0, "(u4-4u2^2)D0^2=R_pole"))

    x20, x21, x22 = sp.symbols("x20 x21 x22")
    xbar = (x20 + 2 * x21 + 2 * x22) / 5
    ax = (2 * x20 - x21 - x22) / 10
    bx = (x21 - x22) / 2
    inv20 = xbar + 4 * ax
    inv21 = xbar - ax + bx
    inv22 = xbar - ax - bx
    checks.append(("grouped_inverse_20", sp.simplify(inv20 - x20) == 0, "x20"))
    checks.append(("grouped_inverse_21", sp.simplify(inv21 - x21) == 0, "x21"))
    checks.append(("grouped_inverse_22", sp.simplify(inv22 - x22) == 0, "x22"))

    lam20, lam21, lam22, eps, x1 = 1, sp.Rational(1, 2), -1, sp.symbols("eps"), sp.symbols("x1")
    y20 = xbar + eps * lam20 * x1
    y21 = xbar + eps * lam21 * x1
    y22 = xbar + eps * lam22 * x1
    ay = sp.simplify((2 * y20 - y21 - y22) / 10)
    by = sp.simplify((y21 - y22) / 2)
    checks.append(("axisymmetric_b_equals_3a", sp.simplify(by - 3 * ay) == 0, "b=3a"))

    G, cs, a, c, mhat, Sport = sp.symbols("G cs a c mhat Sport", positive=True)
    Ptarget = 54 * G * cs**5 / (5 * a**5 * c**5)
    gamma_eff = mhat**2 * Sport * Ptarget * a**5 / (27 * cs**5)
    checks.append(("gamma_target_equivalence", sp.simplify(gamma_eff - Sport * mhat**2 * 2 * G / (5 * c**5)) == 0, "gamma_eff"))

    passed = sum(1 for _, ok, _ in checks if ok)
    return {
        "checks_total": len(checks),
        "checks_passed": passed,
        "checks_failed": len(checks) - passed,
        "details": [{"name": name, "pass": bool(ok), "note": note} for name, ok, note in checks],
    }


# ---------------------------------------------------------------------------
# Branch extraction formulas
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
    Z0 = Z2 = Z4 = 0.0
    N0 = N2 = N4 = 0.0
    diagnostics: List[Dict[str, float]] = []
    for idx, port in enumerate(ports):
        OU = float(port["Omega_U"])
        OW = float(port["Omega_W"])
        R = float(port.get("R", 0.0))
        gU = float(port.get("g_U", 0.0))
        gW = float(port.get("g_W", 0.0))
        Delta = OU**2 * OW**2 - R**2
        S = OU**2 + OW**2
        Q = gU**2 * OW**2 + 2.0 * gU * gW * R + gW**2 * OU**2
        H = gU**2 + gW**2
        P = OU**2 * gW + R * gU
        if Delta == 0:
            raise ValueError(f"mixed port {idx} has Delta=0")
        Z0_r = Q / Delta
        Z2_r = (Q * S - H * Delta) / Delta**2
        Z4_r = (Q * (S**2 - Delta) - S * H * Delta) / Delta**3
        N0_r = P**2 / Delta**2
        N2_r = 2.0 * P * (P * S - Delta * gW) / Delta**3
        N4_r = (Delta**2 * gW**2 - 2.0 * Delta * P**2 - 4.0 * Delta * P * S * gW + 3.0 * P**2 * S**2) / Delta**4
        Z0 += Z0_r
        Z2 += Z2_r
        Z4 += Z4_r
        N0 += N0_r
        N2 += N2_r
        N4 += N4_r
        diagnostics.append({
            "port_index": float(idx),
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
    """Extract all lane-level coefficients and diagnostics."""
    if "direct_coefficients" in lane:
        coeff = {k: float(v) for k, v in lane["direct_coefficients"].items()}
        required = {"K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4"}
        missing = required - set(coeff)
        if missing:
            raise ValueError(f"direct_coefficients missing {sorted(missing)}")
        port_diagnostics: List[Dict[str, float]] = []
    else:
        K = float(lane["K"])
        M = float(lane["M"])
        b = bdg_moments(lane.get("bdg_modes", []))
        z_n, port_diagnostics = mixed_port_moments(lane.get("mixed_ports", []))
        coeff = {"K": K, "M": M, **b, **z_n}

    K = coeff["K"]
    M = coeff["M"]
    B0 = coeff["B0"]
    B2 = coeff["B2"]
    B4 = coeff["B4"]
    Z0 = coeff["Z0"]
    Z2 = coeff["Z2"]
    Z4 = coeff["Z4"]
    N0 = coeff["N0"]
    N2 = coeff["N2"]
    N4 = coeff["N4"]

    D0 = K - B0 - Z0
    D2 = -(M + B2 + Z2)
    D4 = -(B4 + Z4)
    u2 = -D2 / D0 if D0 != 0 else math.nan
    u4 = (D2**2 - D0 * D4) / D0**2 if D0 != 0 else math.nan
    P0 = N0 / D0 if D0 != 0 else math.nan
    P2 = (D0 * N2 - 2.0 * D2 * N0) / D0**2 if D0 != 0 else math.nan
    P4 = (D0**2 * N4 - 2.0 * D0 * (D2 * N2 + D4 * N0) + 3.0 * D2**2 * N0) / D0**3 if D0 != 0 else math.nan

    A = M + B2 + Z2
    C = B4 + Z4
    R_pole = D0 * C - 3.0 * A**2
    constant_prefactor_N2_residual = N2 - 2.0 * D2 * N0 / D0 if D0 != 0 else math.nan
    constant_prefactor_N4_residual = N4 - N0 * (D2**2 + 2.0 * D0 * D4) / D0**2 if D0 != 0 else math.nan

    stability = {
        "D0_positive": D0 > 0.0,
        "C_positive": C > 0.0,
        "N0_nonzero": N0 != 0.0,
        "wall_M_positive": M > 0.0,
    }
    if port_diagnostics:
        stability["all_Delta_positive"] = all(p["Delta"] > 0.0 for p in port_diagnostics)
    else:
        stability["all_Delta_positive"] = True

    return {
        **coeff,
        "D0": D0,
        "D2": D2,
        "D4": D4,
        "u2": u2,
        "u4": u4,
        "P0": P0,
        "P2": P2,
        "P4": P4,
        "A_MplusB2plusZ2": A,
        "C_B4plusZ4": C,
        "R_pole": R_pole,
        "constant_prefactor_N2_residual": constant_prefactor_N2_residual,
        "constant_prefactor_N4_residual": constant_prefactor_N4_residual,
        "port_diagnostics": port_diagnostics,
        "stability": stability,
    }


def extract_branch(branch: Mapping[str, Any], tol: float = 1e-9) -> Dict[str, Any]:
    """Extract the full grouped packet for one branch."""
    geometry = branch.get("geometry", {})
    target = branch.get("target", {})
    lanes_in = branch["lanes"]
    if not all(lane in lanes_in for lane in LANES):
        raise ValueError(f"branch lanes must include {LANES}")

    lane_packets = {lane: lane_extract(lanes_in[lane]) for lane in LANES}

    grouped: Dict[str, Dict[str, float]] = {}
    for key in [
        "K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4",
        "D0", "D2", "D4", "u2", "u4", "P0", "P2", "P4", "R_pole",
    ]:
        grouped[key] = grouped_decomposition(lane_values(lane_packets, key))

    # Isotropic trace residuals from grouped trace quantities.
    D0_bar = grouped["D0"]["bar"]
    B4_bar = grouped["B4"]["bar"]
    Z4_bar = grouped["Z4"]["bar"]
    A_bar = grouped["M"]["bar"] + grouped["B2"]["bar"] + grouped["Z2"]["bar"]
    P0_bar = grouped["P0"]["bar"]
    P2_bar = grouped["P2"]["bar"]
    P4_bar = grouped["P4"]["bar"]

    constants = target.get("constants", {})
    G = float(constants.get("G", 1.0))
    c_s = float(constants.get("c_s", 1.0))
    c = float(constants.get("c", 1.0))
    a = float(constants.get("a", 1.0))
    mhat0 = float(constants.get("mhat0", 1.0))
    S_port = float(constants.get("S_port", 1.0))
    theta_tail = float(constants.get("theta_tail", 1.0))

    P0_target = 54.0 * G * c_s**5 / (5.0 * a**5 * c**5)
    gamma_eff = mhat0**2 * S_port * P0_bar * a**5 / (27.0 * c_s**5)
    gamma_GR = 2.0 * G / (5.0 * c**5)
    R_pole_iso = D0_bar * (B4_bar + Z4_bar) - 3.0 * A_bar**2
    R_norm = mhat0**2 * S_port * P0_bar - P0_target
    R_tail = theta_tail * (c / c_s) ** 3 - 1.0

    open_gate = {
        "R_exit_positive": float(geometry.get("R_exit", 0.0)) > 0.0,
        "boundary_class_open_impedance": geometry.get("boundary_class") == "open_impedance",
        "hard_cap_forbidden": geometry.get("boundary_class") != "hard_cap" and float(geometry.get("R_exit", 0.0)) > 0.0,
    }
    lane_stability_all = {
        name: all(packet["stability"].get(name, False) for packet in lane_packets.values())
        for name in sorted({k for packet in lane_packets.values() for k in packet["stability"]})
    }

    residuals = {
        "R_pole": R_pole_iso,
        "R_norm": R_norm,
        "R_P2": P2_bar,
        "R_P4": P4_bar,
        "R_tail": R_tail,
        "gamma_eff_minus_gamma_GR": gamma_eff - gamma_GR,
    }
    pass_flags = {
        "open_gate_pass": all(open_gate.values()),
        "stability_gate_pass": all(lane_stability_all.values()),
        "isotropic_D0_pass": near_zero(grouped["D0"]["anisotropy_norm_sq"], tol),
        "isotropic_N0_pass": near_zero(grouped["N0"]["anisotropy_norm_sq"], tol),
        "one_pole_pass": near_zero(R_pole_iso, tol),
        "normalization_pass": near_zero(R_norm, tol),
        "constant_prefactor_P2_pass": near_zero(P2_bar, tol),
        "constant_prefactor_P4_pass": near_zero(P4_bar, tol),
        "tail_transport_pass": near_zero(R_tail, tol),
    }
    pass_flags["target_packet_pass"] = all(pass_flags.values())

    frozen_packet = {
        "branch_name": branch.get("name", "unnamed"),
        "geometry": geometry,
        "target_constants": constants,
        "lanes_input": lanes_in,
    }

    return {
        "name": branch.get("name", "unnamed"),
        "input_hash": stable_hash(frozen_packet),
        "open_gate": open_gate,
        "lane_packets": lane_packets,
        "grouped": grouped,
        "target_values": {
            "P0_target": P0_target,
            "gamma_eff": gamma_eff,
            "gamma_GR": gamma_GR,
        },
        "residuals": residuals,
        "stability": lane_stability_all,
        "pass_flags": pass_flags,
    }


# ---------------------------------------------------------------------------
# Built-in fixtures
# ---------------------------------------------------------------------------


def calibrated_coefficient_branch() -> Dict[str, Any]:
    """A direct-coefficient branch that exactly satisfies the algebraic target surface."""
    # Normalized constants: G=c=c_s=a=mhat=S_port=theta_tail=1.
    # Choose D0=1, A=M+B2+Z2=1, C=B4+Z4=3, hence one-pole D0*C=3A^2.
    # Set B=Z=0 except M=1 and B4=3, K=1. Then D2=-1, D4=-3.
    # Target P0 = 54/5. Constant-prefactor branch gives N2 and N4 below.
    Ptarget = 54.0 / 5.0
    direct = {
        "K": 1.0,
        "M": 1.0,
        "B0": 0.0,
        "B2": 0.0,
        "B4": 3.0,
        "Z0": 0.0,
        "Z2": 0.0,
        "Z4": 0.0,
        "N0": Ptarget,
        "N2": -2.0 * Ptarget,
        "N4": -5.0 * Ptarget,
    }
    return {
        "name": "calibrated_direct_coefficient_fixture",
        "geometry": {"R_exit": 1.0, "boundary_class": "open_impedance", "Y_L_limit": 0.0},
        "target": {"constants": {"G": 1.0, "c_s": 1.0, "c": 1.0, "a": 1.0, "mhat0": 1.0, "S_port": 1.0, "theta_tail": 1.0}},
        "lanes": {lane: {"direct_coefficients": copy.deepcopy(direct)} for lane in LANES},
    }


def primitive_open_branch() -> Dict[str, Any]:
    """A stable open-throat primitive branch used as a falsifying extraction demo."""
    lane = {
        "K": 1.6,
        "M": 0.9,
        "bdg_modes": [
            {"coupling": 0.42, "varpi": 2.7},
            {"coupling": 0.18, "varpi": 4.1},
        ],
        "mixed_ports": [
            {"Omega_U": 3.0, "Omega_W": 4.0, "R": 0.65, "g_U": 0.28, "g_W": 0.52},
        ],
    }
    return {
        "name": "primitive_open_throat_demo_branch",
        "geometry": {"R_exit": 0.35, "boundary_class": "open_impedance", "Y_L_limit": 0.0},
        "target": {"constants": {"G": 1.0, "c_s": 1.0, "c": 1.0, "a": 1.0, "mhat0": 1.0, "S_port": 1.0, "theta_tail": 1.0}},
        "lanes": {lane_id: copy.deepcopy(lane) for lane_id in LANES},
    }


def default_manifest() -> Dict[str, Any]:
    return {
        "schema": "stage_v2_21_branch_extraction_fixture/v1",
        "description": "Built-in V2-21 manifest with one calibrated direct-coefficient test and one primitive open-throat demo branch.",
        "branches": [calibrated_coefficient_branch(), primitive_open_branch()],
    }


# ---------------------------------------------------------------------------
# Main program
# ---------------------------------------------------------------------------


def load_manifest(path: Optional[str]) -> Dict[str, Any]:
    if path is None:
        return default_manifest()
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage V2-21 branch extraction fixture")
    parser.add_argument("--manifest", default=None, help="Optional branch manifest JSON file")
    parser.add_argument("--out-json", default=None, help="Optional path for the extracted observable packet JSON")
    parser.add_argument("--write-sample-manifest", default=None, help="Optional path for writing the built-in sample manifest")
    parser.add_argument("--tol", type=float, default=1e-9, help="Residual tolerance for pass/fail flags")
    args = parser.parse_args()

    if args.write_sample_manifest:
        with open(args.write_sample_manifest, "w", encoding="utf-8") as f:
            json.dump(default_manifest(), f, indent=2, sort_keys=True)

    manifest = load_manifest(args.manifest)
    audit = run_symbolic_audit()
    branch_packets = [extract_branch(branch, tol=args.tol) for branch in manifest.get("branches", [])]
    result = {
        "script": "stage_v2_21_branch_extraction_fixture.py",
        "symbolic_audit": audit,
        "manifest_hash": stable_hash(manifest),
        "branch_count": len(branch_packets),
        "branches": branch_packets,
    }

    if args.out_json:
        with open(args.out_json, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, sort_keys=True)

    print("STAGE V2-21 BRANCH EXTRACTION FIXTURE")
    print(f"symbolic_checks: {audit['checks_passed']}/{audit['checks_total']} passed")
    print(f"manifest_hash: {result['manifest_hash']}")
    for packet in branch_packets:
        flags = packet["pass_flags"]
        residuals = packet["residuals"]
        target = packet["target_values"]
        print("---")
        print(f"branch: {packet['name']}")
        print(f"input_hash: {packet['input_hash']}")
        print(f"open_gate_pass: {flags['open_gate_pass']}")
        print(f"stability_gate_pass: {flags['stability_gate_pass']}")
        print(f"target_packet_pass: {flags['target_packet_pass']}")
        print(f"D0_bar: {packet['grouped']['D0']['bar']:.16g}")
        print(f"N0_bar: {packet['grouped']['N0']['bar']:.16g}")
        print(f"P0_bar: {packet['grouped']['P0']['bar']:.16g}")
        print(f"P0_target: {target['P0_target']:.16g}")
        print(f"R_pole: {residuals['R_pole']:.16g}")
        print(f"R_norm: {residuals['R_norm']:.16g}")
        print(f"R_P2: {residuals['R_P2']:.16g}")
        print(f"R_P4: {residuals['R_P4']:.16g}")
        print(f"R_tail: {residuals['R_tail']:.16g}")
        print(f"D0_anisotropy_norm_sq: {packet['grouped']['D0']['anisotropy_norm_sq']:.16g}")
        print(f"P0_anisotropy_norm_sq: {packet['grouped']['P0']['anisotropy_norm_sq']:.16g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
