#!/usr/bin/env python3
"""
Stage V2-22B — Numerical branch profile schema and solver handoff validator.

This script freezes the data format that a moving-throat PDE/continuation solver
must export before the V2-22A profile-to-coefficient adapter and V2-21 observable
extractor are allowed to evaluate target residuals.

The validator is intentionally conservative:
  * it enforces the open-throat patch R_exit > 0 and boundary_class=open_impedance;
  * it requires a pre-target freeze certificate;
  * it checks monotone finite-throat grids and sampled profile lengths;
  * it verifies the minimum profile set needed to populate BdG and Maxwell/mixed
    coefficients;
  * it builds a V2-22A-compatible profile manifest without inspecting any target
    residuals;
  * it rejects hard-cap, missing-profile, nonpositive-frequency, and nonpositive
    mixed-block Delta data before coefficient extraction.

Run examples:
  python stage_v2_22b_solver_handoff_validator.py \
      --write-schema stage_v2_22b_solver_output_schema.json \
      --write-valid stage_v2_22b_sample_solver_output_valid.json \
      --write-invalid stage_v2_22b_sample_solver_output_invalid_hardcap.json \
      --out-profile-manifest stage_v2_22b_generated_profile_manifest.json \
      --out-report stage_v2_22b_validation_report.json

The built-in valid sample is a solver-export version of the finite open-throat
D/N profile prototype.  It is not tuned to the GR normalization target; this
stage is a handoff and validation layer, not branch realization.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import sympy as sp

LANES = ("20", "21", "22")
AXISYM_LAMBDA = {"20": 1.0, "21": 0.5, "22": -1.0}
DIRECT_BDG_KEYS = ("B0", "B2", "B4")
DIRECT_TRANSFER_KEYS = ("Z0", "Z2", "Z4", "N0", "N2", "N4")
DIRECT_REQUIRED_CHECKS = (
    "current_frechet_matches_step8c",
    "outgoing_flux_positive",
    "open_not_hard_cap",
    "pure_gauge_zero_physical_transfer",
    "basis_invariance",
    "v2_09_regression",
    "green_residuals_small",
    "bdg_residuals_small",
    "N0_positive",
)


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def finite_number(x: Any) -> bool:
    try:
        y = float(x)
        return math.isfinite(y)
    except Exception:
        return False


def add_issue(issues: List[Dict[str, str]], severity: str, path: str, message: str) -> None:
    issues.append({"severity": severity, "path": path, "message": message})


def trapezoid(xs: Sequence[float], ys: Sequence[float]) -> float:
    if len(xs) != len(ys):
        raise ValueError("xs and ys lengths differ")
    if len(xs) < 2:
        raise ValueError("need at least two grid points")
    total = 0.0
    for i in range(len(xs) - 1):
        total += 0.5 * (float(ys[i]) + float(ys[i + 1])) * (float(xs[i + 1]) - float(xs[i]))
    return total


def product_integral(xs: Sequence[float], *ys_list: Sequence[float]) -> float:
    prod: List[float] = []
    for i in range(len(xs)):
        val = 1.0
        for ys in ys_list:
            val *= float(ys[i])
        prod.append(val)
    return trapezoid(xs, prod)


def sampled_profile(xs: Sequence[float], ys: Sequence[float]) -> Dict[str, Any]:
    return {"kind": "sampled", "samples": [[float(x), float(y)] for x, y in zip(xs, ys)]}


def direct_coefficients_from_packet(packet: Mapping[str, Any]) -> Dict[str, float]:
    """Build the V2-21 direct-coefficient lane from a validated direct packet."""
    bdg_wall = packet["derived_bdg_wall_coefficients"]["coefficients"]
    transfer = packet["derived_maxwell_transfer"]["coefficients"]
    return {
        "K": float(packet["wall"]["K"]),
        "M": float(packet["wall"]["M"]),
        "B0": float(bdg_wall["B0"]),
        "B2": float(bdg_wall["B2"]),
        "B4": float(bdg_wall["B4"]),
        "Z0": float(transfer["Z0"]),
        "Z2": float(transfer["Z2"]),
        "Z4": float(transfer["Z4"]),
        "N0": float(transfer["N0"]),
        "N2": float(transfer["N2"]),
        "N4": float(transfer["N4"]),
    }


# ---------------------------------------------------------------------------
# Symbolic audit
# ---------------------------------------------------------------------------


def run_symbolic_audit() -> Dict[str, Any]:
    """Verify the mathematical identities behind the handoff protocol."""
    checks: List[Tuple[str, bool, str]] = []
    s, L, k, YL, T, qL = sp.symbols("s L k YL T qL", positive=True)

    # D/N half-wave and endpoint conditions for the support field q(s)=sin(k s).
    q = sp.sin(k * s)
    mouth = sp.simplify(q.subs(s, 0))
    bottom_deriv = sp.simplify(sp.diff(q, s).subs({s: L, k: sp.pi / (2 * L)}))
    checks.append(("mouth_Dirichlet_q0", mouth == 0, str(mouth)))
    checks.append(("bottom_Neumann_halfwave", bottom_deriv == 0, str(bottom_deriv)))

    # Robin impedance endpoint reduces to Neumann as Y_L -> 0.
    robin = T * sp.Symbol("q_wL") + YL * qL
    robin_neumann = sp.simplify(sp.limit(robin, YL, 0))
    checks.append(("open_impedance_robin_to_neumann", robin_neumann == T * sp.Symbol("q_wL"), str(robin_neumann)))

    # Exact D/N overlap used by the profile handoff prototype.
    chi = sp.sqrt(2 / L) * sp.sin(sp.pi * s / L)
    phi = sp.sqrt(2 / L) * sp.sin(sp.pi * s / (2 * L))
    I_chi_chi = sp.simplify(sp.integrate(chi * chi, (s, 0, L)))
    I_phi_phi = sp.simplify(sp.integrate(phi * phi, (s, 0, L)))
    I_chi_phi = sp.simplify(sp.integrate(chi * phi, (s, 0, L)))
    checks.append(("chi_norm", I_chi_chi == 1, str(I_chi_chi)))
    checks.append(("phi_norm", I_phi_phi == 1, str(I_phi_phi)))
    checks.append(("DN_overlap", sp.simplify(I_chi_phi - 8 / (3 * sp.pi)) == 0, str(I_chi_phi)))

    # Weak-axisymmetric grouped signature survives the handoff schema.
    x0, eps, x1 = sp.symbols("x0 eps x1")
    y20 = x0 + eps * x1
    y21 = x0 + eps * sp.Rational(1, 2) * x1
    y22 = x0 - eps * x1
    abar = sp.simplify((y20 + 2 * y21 + 2 * y22) / 5)
    a = sp.simplify((2 * y20 - y21 - y22) / 10)
    b = sp.simplify((y21 - y22) / 2)
    checks.append(("axisym_trace_unchanged", sp.simplify(abar - x0) == 0, str(abar)))
    checks.append(("axisym_b_equals_3a", sp.simplify(b - 3 * a) == 0, f"a={a}, b={b}"))

    # Mixed block stability determinant and outgoing transfer are nonnegative if Delta>0.
    OU, OW, R, gU, gW = sp.symbols("OU OW R gU gW", positive=True)
    Delta = OU**2 * OW**2 - R**2
    N0 = (OU**2 * gW + R * gU) ** 2 / Delta**2
    checks.append(("N0_perfect_square_over_Delta2", sp.factor(sp.together(N0).as_numer_denom()[0]) == (OU**2 * gW + R * gU) ** 2, str(N0)))

    passed = sum(1 for _, ok, _ in checks if ok)
    return {
        "checks_total": len(checks),
        "checks_passed": passed,
        "checks_failed": len(checks) - passed,
        "details": [{"name": n, "pass": bool(ok), "note": note} for n, ok, note in checks],
    }


# ---------------------------------------------------------------------------
# Schema document
# ---------------------------------------------------------------------------


def solver_output_schema() -> Dict[str, Any]:
    """Return a human-readable JSON schema-like contract for solver exports."""
    return {
        "schema": "stage_v2_22b_solver_handoff_schema/v1",
        "required_top_level_keys": [
            "schema",
            "branch_id",
            "freeze",
            "geometry",
            "constants",
            "grid",
            "wall",
            "profiles",
            "bdg_modes",
            "solver_metadata",
        ],
        "coefficient_paths": {
            "mixed_ports": "legacy path; at least one Maxwell/A_w mixed port is required unless derived_maxwell_transfer is supplied",
            "derived_maxwell_transfer": "optional direct-derived path; when present without mixed_ports, requires derived_bdg_wall_coefficients and generates V2-21 direct_coefficients",
        },
        "freeze_required": {
            "pre_target_freeze": True,
            "target_blind": True,
            "parent_action": "description of GNLS + localized Maxwell + S_eta or effective wall closure",
            "gauge_convention": "localized/constrained Lorenz or other declared convention",
            "boundary_protocol": "open_impedance_AC_reflecting_DC_leaking",
        },
        "geometry_required": {
            "L": "positive throat length",
            "R_mouth": "positive mouth radius",
            "R_exit": "strictly positive open exit radius",
            "boundary_class": "must be open_impedance",
            "Y_L_limit": "finite load-admittance parameter; zero gives AC Neumann limit",
            "exit_model": "description of impedance mismatch / organ-pipe opening",
        },
        "grid_required": {
            "coordinate": "s",
            "points": "strictly increasing list from 0 to L",
        },
        "profiles_required": {
            "weight": "mu_s values on grid",
            "wall_chi_eta": "wall/support profile values on grid",
        },
        "bdg_modes_each_required": ["name", "lambda_B", "varpi", "profile_values"],
        "mixed_ports_each_required": [
            "name",
            "lambda_U",
            "Omega_U",
            "u_values",
            "lambda_W",
            "Omega_W",
            "w_values",
            "lambda_R",
        ],
        "derived_bdg_wall_coefficients_required_for_direct_path": {
            "status": "derived_bdg_wall_m1b or equivalent",
            "coefficients": ["B0", "B2", "B4"],
            "source_hashes": "hashes for upstream M1b packet/diagnostics",
        },
        "derived_maxwell_transfer_required_for_direct_path": {
            "status": "derived_green_function_transfer",
            "gauge_convention": "declared transfer gauge convention",
            "flux_normalization": "Gamma_port / sigma_Q^can convention",
            "coefficients": ["Z0", "Z2", "Z4", "N0", "N2", "N4"],
            "operator_gauge_residual_metrics": "finite residual metrics and boolean validation checks; do not use the forbidden key 'residuals'",
            "source_hashes": "hashes for upstream Spike-2 transfer diagnostics",
        },
        "optional": {
            "weak_axisymmetric": "epsilon and primitive slopes with grouped signature (20,21,22)=(1,1/2,-1)",
            "normalization_tolerances": "norm and overlap tolerances used by the validator",
        },
        "solver_metadata_required": {
            "exporter": "solver/exporter name",
            "mesh_points": "integer number of axial grid points exported",
            "nonlinear_residual_norm": "finite nonnegative frozen solver residual norm",
            "coefficient_family": "declared reduced-coefficient family",
            "source_commit": "solver source revision, fixture-only tag, or equivalent provenance id",
        },
        "forbidden": {
            "geometry.R_exit": "0 or negative",
            "geometry.boundary_class": "hard_cap",
            "freeze.pre_target_freeze": False,
            "any target residual": "must not be used to edit this solver output after validation",
        },
    }


# ---------------------------------------------------------------------------
# Built-in sample solver exports
# ---------------------------------------------------------------------------


def build_sample_solver_output(valid: bool = True, points: int = 801) -> Dict[str, Any]:
    L = 1.0
    xs = [L * i / (points - 1) for i in range(points)]
    wall = [math.sqrt(2.0 / L) * math.sin(math.pi * x / L) for x in xs]
    half = [math.sqrt(2.0 / L) * math.sin(math.pi * x / (2.0 * L)) for x in xs]
    weight = [1.0 for _ in xs]
    return {
        "schema": "stage_v2_22b_solver_handoff/v1",
        "branch_id": "sample_open_DN_solver_export" if valid else "invalid_hardcap_solver_export",
        "freeze": {
            "pre_target_freeze": True,
            "target_blind": True,
            "parent_action": "GNLS + localized Maxwell + effective S_eta^(2) wall closure",
            "gauge_convention": "localized_or_constrained_lorenz_declared",
            "boundary_protocol": "open_impedance_AC_reflecting_DC_leaking" if valid else "hard_cap_DC_blocked",
            "source_map_convention": "canonical real STF grouped P2 basis",
            "support_profile_family": "sampled finite-throat DN prototype",
        },
        "geometry": {
            "L": L,
            "R_mouth": 1.0,
            "R_exit": 0.35 if valid else 0.0,
            "boundary_class": "open_impedance" if valid else "hard_cap",
            "Y_L_limit": 0.0,
            "exit_model": "organ_pipe_low_impedance_expansion" if valid else "capped_tube",
        },
        "constants": {"G": 1.0, "c_s": 1.0, "c": 1.0, "a": 1.0, "mhat0": 1.0, "S_port": 1.0, "theta_tail": 1.0},
        "grid": {"coordinate": "s", "points": xs},
        "wall": {"K": 1.6, "M": 0.9},
        "profiles": {
            "weight": weight,
            "wall_chi_eta": wall,
        },
        "bdg_modes": [
            {"name": "bdg_halfwave", "lambda_B": 0.42, "varpi": 2.7, "profile_values": half}
        ],
        "mixed_ports": [
            {
                "name": "one_mixed_port",
                "lambda_U": 0.28,
                "Omega_U": 3.0,
                "u_values": wall,
                "lambda_W": 0.52,
                "Omega_W": 4.0,
                "w_values": half,
                "lambda_R": 0.65,
            }
        ],
        "weak_axisymmetric": None,
        "solver_metadata": {
            "exporter": "stage_v2_22b_builtin_sample",
            "mesh_points": points,
            "nonlinear_residual_norm": 0.0,
            "coefficient_family": "sampled_finite_throat_DN_validation_fixture_v1",
            "source_commit": "fixture-only",
            "export_timestamp": "2026-04-24T00:00:00Z",
            "notes": "This is a validation fixture, not a target-realized branch.",
        },
        "normalization_tolerances": {"profile_norm_tol": 5e-3, "Delta_tol": 1e-12},
    }


# ---------------------------------------------------------------------------
# Validation and conversion
# ---------------------------------------------------------------------------


def validate_solver_output(packet: Mapping[str, Any]) -> Dict[str, Any]:
    issues: List[Dict[str, str]] = []
    required = ["schema", "branch_id", "freeze", "geometry", "constants", "grid", "wall", "profiles", "bdg_modes", "solver_metadata"]
    for key in required:
        if key not in packet:
            add_issue(issues, "error", key, "missing required top-level key")

    has_direct_derived_transfer = isinstance(packet.get("derived_maxwell_transfer"), Mapping)
    if "derived_maxwell_transfer" in packet and not has_direct_derived_transfer:
        add_issue(issues, "error", "derived_maxwell_transfer", "derived transfer branch must be an object")

    forbidden_target_keys = {
        "P0_target",
        "R_pole",
        "R_norm",
        "R_P2",
        "R_P4",
        "gamma_eff",
        "gamma_GR",
        "pass_flags",
        "residuals",
        "target_packet_pass",
        "target_values",
    }

    def scan_for_target_leaks(obj: Any, path: str = "") -> None:
        if isinstance(obj, Mapping):
            for key, value in obj.items():
                key_str = str(key)
                child_path = f"{path}.{key_str}" if path else key_str
                if key_str in forbidden_target_keys:
                    add_issue(issues, "error", child_path, "target residual/output fields are forbidden in frozen solver exports")
                scan_for_target_leaks(value, child_path)
        elif isinstance(obj, list):
            for idx, value in enumerate(obj):
                scan_for_target_leaks(value, f"{path}[{idx}]")

    scan_for_target_leaks(packet)

    if packet.get("schema") != "stage_v2_22b_solver_handoff/v1":
        add_issue(issues, "error", "schema", "schema must be stage_v2_22b_solver_handoff/v1")

    freeze = packet.get("freeze", {}) if isinstance(packet.get("freeze"), Mapping) else {}
    if freeze.get("pre_target_freeze") is not True:
        add_issue(issues, "error", "freeze.pre_target_freeze", "branch packet must be frozen before target residual evaluation")
    if freeze.get("target_blind") is not True:
        add_issue(issues, "error", "freeze.target_blind", "solver export must be target-blind")
    if "gauge_convention" not in freeze:
        add_issue(issues, "error", "freeze.gauge_convention", "gauge convention must be declared")
    if freeze.get("boundary_protocol") != "open_impedance_AC_reflecting_DC_leaking":
        add_issue(issues, "error", "freeze.boundary_protocol", "must declare open impedance AC/DC split")

    geometry = packet.get("geometry", {}) if isinstance(packet.get("geometry"), Mapping) else {}
    L = geometry.get("L")
    if not finite_number(L) or float(L) <= 0.0:
        add_issue(issues, "error", "geometry.L", "L must be positive")
    if geometry.get("boundary_class") != "open_impedance":
        add_issue(issues, "error", "geometry.boundary_class", "boundary_class must be open_impedance")
    if not finite_number(geometry.get("R_exit")) or float(geometry.get("R_exit", 0.0)) <= 0.0:
        add_issue(issues, "error", "geometry.R_exit", "R_exit must be strictly positive; hard cap is forbidden")
    if geometry.get("boundary_class") == "hard_cap":
        add_issue(issues, "error", "geometry.boundary_class", "hard_cap geometry is forbidden by V2-04 patch")

    # Grid and profile shape checks.
    xs: List[float] = []
    grid = packet.get("grid", {}) if isinstance(packet.get("grid"), Mapping) else {}
    raw_points = grid.get("points", [])
    if not isinstance(raw_points, list) or len(raw_points) < 2:
        add_issue(issues, "error", "grid.points", "grid.points must contain at least two samples")
    else:
        try:
            xs = [float(x) for x in raw_points]
            if any(not math.isfinite(x) for x in xs):
                raise ValueError("non-finite grid value")
            if any(xs[i + 1] <= xs[i] for i in range(len(xs) - 1)):
                add_issue(issues, "error", "grid.points", "grid must be strictly increasing")
            if finite_number(L):
                if abs(xs[0]) > 1e-12:
                    add_issue(issues, "error", "grid.points[0]", "grid must begin at 0")
                if abs(xs[-1] - float(L)) > 1e-9:
                    add_issue(issues, "error", "grid.points[-1]", "grid must end at L")
        except Exception as exc:
            add_issue(issues, "error", "grid.points", f"invalid grid values: {exc}")

    def check_values(path: str, values: Any) -> Optional[List[float]]:
        if not isinstance(values, list):
            add_issue(issues, "error", path, "must be a list of numbers on the grid")
            return None
        if xs and len(values) != len(xs):
            add_issue(issues, "error", path, f"length {len(values)} does not match grid length {len(xs)}")
            return None
        out: List[float] = []
        for i, val in enumerate(values):
            if not finite_number(val):
                add_issue(issues, "error", f"{path}[{i}]", "non-finite profile value")
                return None
            out.append(float(val))
        return out

    profiles = packet.get("profiles", {}) if isinstance(packet.get("profiles"), Mapping) else {}
    weight = check_values("profiles.weight", profiles.get("weight"))
    wall = check_values("profiles.wall_chi_eta", profiles.get("wall_chi_eta"))

    diagnostics: Dict[str, Any] = {}
    tol = float(packet.get("normalization_tolerances", {}).get("profile_norm_tol", 5e-3)) if isinstance(packet.get("normalization_tolerances", {}), Mapping) else 5e-3
    if xs and weight and wall:
        wall_norm = product_integral(xs, weight, wall, wall)
        diagnostics["wall_norm"] = wall_norm
        if abs(wall_norm - 1.0) > tol:
            add_issue(issues, "warning", "profiles.wall_chi_eta", f"wall norm is {wall_norm:.8g}, expected near 1")

    wall_packet = packet.get("wall", {}) if isinstance(packet.get("wall"), Mapping) else {}
    for key in ["K", "M"]:
        if not finite_number(wall_packet.get(key)) or float(wall_packet.get(key, 0.0)) <= 0.0:
            add_issue(issues, "error", f"wall.{key}", f"wall.{key} must be positive")

    if has_direct_derived_transfer:
        bdg_wall_direct = packet.get("derived_bdg_wall_coefficients")
        if not isinstance(bdg_wall_direct, Mapping):
            add_issue(issues, "error", "derived_bdg_wall_coefficients", "direct-derived path requires derived BdG/wall coefficients")
        else:
            if not isinstance(bdg_wall_direct.get("status"), str) or not bdg_wall_direct.get("status", "").strip():
                add_issue(issues, "error", "derived_bdg_wall_coefficients.status", "derived BdG/wall status must be declared")
            bdg_coeff = bdg_wall_direct.get("coefficients")
            if not isinstance(bdg_coeff, Mapping):
                add_issue(issues, "error", "derived_bdg_wall_coefficients.coefficients", "must provide B0/B2/B4 coefficients")
            else:
                for key in DIRECT_BDG_KEYS:
                    if not finite_number(bdg_coeff.get(key)):
                        add_issue(issues, "error", f"derived_bdg_wall_coefficients.coefficients.{key}", f"{key} must be finite")
                    elif float(bdg_coeff[key]) < 0.0:
                        add_issue(issues, "error", f"derived_bdg_wall_coefficients.coefficients.{key}", f"{key} must be nonnegative")
                for key in ["K", "M"]:
                    if key in bdg_coeff:
                        if not finite_number(bdg_coeff.get(key)):
                            add_issue(issues, "error", f"derived_bdg_wall_coefficients.coefficients.{key}", f"{key} must be finite")
                        elif finite_number(wall_packet.get(key)) and abs(float(bdg_coeff[key]) - float(wall_packet[key])) > 1e-10:
                            add_issue(issues, "error", f"derived_bdg_wall_coefficients.coefficients.{key}", f"{key} must match wall.{key}")
            source_hashes = bdg_wall_direct.get("source_hashes")
            if not isinstance(source_hashes, Mapping) or not source_hashes:
                add_issue(issues, "error", "derived_bdg_wall_coefficients.source_hashes", "M1b source hashes must be recorded")

        transfer = packet["derived_maxwell_transfer"]
        if transfer.get("status") != "derived_green_function_transfer":
            add_issue(issues, "error", "derived_maxwell_transfer.status", "status must be derived_green_function_transfer")
        if not isinstance(transfer.get("gauge_convention"), str) or not transfer.get("gauge_convention", "").strip():
            add_issue(issues, "error", "derived_maxwell_transfer.gauge_convention", "transfer gauge convention must be declared")
        flux = transfer.get("flux_normalization")
        if not isinstance(flux, Mapping):
            add_issue(issues, "error", "derived_maxwell_transfer.flux_normalization", "Gamma_port flux normalization must be recorded")
        else:
            gamma_port = flux.get("Gamma_port", flux.get("value"))
            if not finite_number(gamma_port) or float(gamma_port) <= 0.0:
                add_issue(issues, "error", "derived_maxwell_transfer.flux_normalization.Gamma_port", "Gamma_port must be positive and finite")
        transfer_coeff = transfer.get("coefficients")
        if not isinstance(transfer_coeff, Mapping):
            add_issue(issues, "error", "derived_maxwell_transfer.coefficients", "must provide Z0/Z2/Z4 and N0/N2/N4")
        else:
            for key in DIRECT_TRANSFER_KEYS:
                if not finite_number(transfer_coeff.get(key)):
                    add_issue(issues, "error", f"derived_maxwell_transfer.coefficients.{key}", f"{key} must be finite")
            if finite_number(transfer_coeff.get("N0")) and float(transfer_coeff["N0"]) <= 0.0:
                add_issue(issues, "error", "derived_maxwell_transfer.coefficients.N0", "N0 must be positive for the outgoing transfer")
        metrics = transfer.get("operator_gauge_residual_metrics")
        if not isinstance(metrics, Mapping):
            add_issue(issues, "error", "derived_maxwell_transfer.operator_gauge_residual_metrics", "operator/gauge residual metrics must be recorded")
        else:
            for key in DIRECT_REQUIRED_CHECKS:
                if metrics.get(key) is not True:
                    add_issue(issues, "error", f"derived_maxwell_transfer.operator_gauge_residual_metrics.{key}", f"{key} must be true")
            for key, value in metrics.items():
                if isinstance(value, bool):
                    continue
                if isinstance(value, (int, float)) and (not math.isfinite(float(value)) or float(value) < 0.0):
                    add_issue(issues, "error", f"derived_maxwell_transfer.operator_gauge_residual_metrics.{key}", "numeric residual metric must be finite and nonnegative")
        source_hashes = transfer.get("source_hashes")
        if not isinstance(source_hashes, Mapping) or not source_hashes:
            add_issue(issues, "error", "derived_maxwell_transfer.source_hashes", "Spike-2 source hashes must be recorded")

        if (
            isinstance(bdg_wall_direct, Mapping)
            and isinstance(bdg_wall_direct.get("coefficients"), Mapping)
            and isinstance(transfer.get("coefficients"), Mapping)
        ):
            try:
                diagnostics["direct_derived_coefficients"] = direct_coefficients_from_packet(packet)
                diagnostics["coefficient_path"] = "direct_derived_coefficients"
            except Exception as exc:
                add_issue(issues, "error", "derived_maxwell_transfer", f"could not build direct coefficients: {exc}")

    bdg_modes = packet.get("bdg_modes", [])
    if not isinstance(bdg_modes, list) or len(bdg_modes) == 0:
        add_issue(issues, "error", "bdg_modes", "at least one BdG mode is required for this handoff")
    else:
        bdg_diag = []
        for idx, mode in enumerate(bdg_modes):
            path = f"bdg_modes[{idx}]"
            if not isinstance(mode, Mapping):
                add_issue(issues, "error", path, "mode must be an object")
                continue
            for key in ["name", "lambda_B", "varpi", "profile_values"]:
                if key not in mode:
                    add_issue(issues, "error", f"{path}.{key}", "missing required BdG field")
            if not finite_number(mode.get("lambda_B")):
                add_issue(issues, "error", f"{path}.lambda_B", "lambda_B must be finite")
            if not finite_number(mode.get("varpi")) or float(mode.get("varpi", 0.0)) <= 0.0:
                add_issue(issues, "error", f"{path}.varpi", "varpi must be positive")
            phi = check_values(f"{path}.profile_values", mode.get("profile_values"))
            if xs and weight and phi:
                norm = product_integral(xs, weight, phi, phi)
                bdg_diag.append({"name": mode.get("name", str(idx)), "profile_norm": norm})
                if abs(norm - 1.0) > tol:
                    add_issue(issues, "warning", f"{path}.profile_values", f"BdG profile norm is {norm:.8g}, expected near 1")
        diagnostics["bdg_profiles"] = bdg_diag

    mixed_ports = packet.get("mixed_ports", [])
    if not isinstance(mixed_ports, list) or len(mixed_ports) == 0:
        if not has_direct_derived_transfer:
            add_issue(issues, "error", "mixed_ports", "at least one mixed Maxwell/A_w port is required unless derived_maxwell_transfer is supplied")
    else:
        port_diag = []
        for idx, port in enumerate(mixed_ports):
            path = f"mixed_ports[{idx}]"
            if not isinstance(port, Mapping):
                add_issue(issues, "error", path, "port must be an object")
                continue
            for key in ["name", "lambda_U", "Omega_U", "u_values", "lambda_W", "Omega_W", "w_values", "lambda_R"]:
                if key not in port:
                    add_issue(issues, "error", f"{path}.{key}", "missing required mixed-port field")
            for key in ["lambda_U", "lambda_W", "lambda_R"]:
                if not finite_number(port.get(key)):
                    add_issue(issues, "error", f"{path}.{key}", f"{key} must be finite")
            for key in ["Omega_U", "Omega_W"]:
                if not finite_number(port.get(key)) or float(port.get(key, 0.0)) <= 0.0:
                    add_issue(issues, "error", f"{path}.{key}", f"{key} must be positive")
            u_values = check_values(f"{path}.u_values", port.get("u_values"))
            w_values = check_values(f"{path}.w_values", port.get("w_values"))
            if xs and weight and u_values and w_values:
                I_uw = product_integral(xs, weight, u_values, w_values)
                R_eff = float(port.get("lambda_R", 0.0)) * I_uw if finite_number(port.get("lambda_R")) else math.nan
                OU = float(port.get("Omega_U", math.nan)) if finite_number(port.get("Omega_U")) else math.nan
                OW = float(port.get("Omega_W", math.nan)) if finite_number(port.get("Omega_W")) else math.nan
                Delta = OU**2 * OW**2 - R_eff**2 if math.isfinite(OU) and math.isfinite(OW) and math.isfinite(R_eff) else math.nan
                port_diag.append({"name": port.get("name", str(idx)), "I_u_w": I_uw, "R_eff": R_eff, "Delta_eff": Delta})
                Delta_tol = float(packet.get("normalization_tolerances", {}).get("Delta_tol", 1e-12)) if isinstance(packet.get("normalization_tolerances", {}), Mapping) else 1e-12
                if not math.isfinite(Delta) or Delta <= Delta_tol:
                    add_issue(issues, "error", f"{path}.Delta_eff", "effective mixed-block Delta must be positive")
        diagnostics["mixed_ports"] = port_diag

    # Constants that are needed by V2-21/22A.
    constants = packet.get("constants", {}) if isinstance(packet.get("constants"), Mapping) else {}
    for key in ["G", "c_s", "c", "a", "mhat0", "S_port", "theta_tail"]:
        if not finite_number(constants.get(key)):
            add_issue(issues, "error", f"constants.{key}", "missing or non-finite target/extraction constant")

    metadata = packet.get("solver_metadata", {}) if isinstance(packet.get("solver_metadata"), Mapping) else {}
    if not isinstance(packet.get("solver_metadata", {}), Mapping):
        add_issue(issues, "error", "solver_metadata", "solver_metadata must be an object")
    else:
        for key in ["exporter", "coefficient_family", "source_commit"]:
            if not isinstance(metadata.get(key), str) or not str(metadata.get(key)).strip():
                add_issue(issues, "error", f"solver_metadata.{key}", f"{key} must be a nonempty provenance string")
        if not finite_number(metadata.get("mesh_points")) or int(float(metadata.get("mesh_points", 0))) <= 0:
            add_issue(issues, "error", "solver_metadata.mesh_points", "mesh_points must be a positive integer")
        elif xs and int(float(metadata["mesh_points"])) != len(xs):
            add_issue(issues, "error", "solver_metadata.mesh_points", "mesh_points must match grid length")
        for key in ["nonlinear_residual_norm", "linear_residual_norm", "continuation_residual_norm"]:
            if key in metadata and (not finite_number(metadata.get(key)) or float(metadata.get(key)) < 0.0):
                add_issue(issues, "error", f"solver_metadata.{key}", f"{key} must be finite and nonnegative")
        diagnostics["solver_metadata"] = {
            "exporter": metadata.get("exporter"),
            "coefficient_family": metadata.get("coefficient_family"),
            "source_commit": metadata.get("source_commit"),
            "mesh_points": metadata.get("mesh_points"),
        }

    error_count = sum(1 for item in issues if item["severity"] == "error")
    warning_count = sum(1 for item in issues if item["severity"] == "warning")
    report = {
        "schema": "stage_v2_22b_validation_report/v1",
        "branch_id": packet.get("branch_id"),
        "packet_hash": stable_hash(packet),
        "validation_pass": error_count == 0,
        "error_count": error_count,
        "warning_count": warning_count,
        "issues": issues,
        "diagnostics": diagnostics,
    }
    return report


def convert_solver_output_to_v22a_profile_manifest(packet: Mapping[str, Any]) -> Dict[str, Any]:
    """Convert a validated solver output into the V2-22A profile adapter schema."""
    report = validate_solver_output(packet)
    if not report["validation_pass"]:
        raise ValueError("solver output did not pass validation; refuse to convert")

    xs = [float(x) for x in packet["grid"]["points"]]
    profiles: Dict[str, Any] = {
        "weight": sampled_profile(xs, packet["profiles"]["weight"]),
        "wall": sampled_profile(xs, packet["profiles"]["wall_chi_eta"]),
    }
    bdg_modes = []
    for idx, mode in enumerate(packet["bdg_modes"]):
        pname = f"bdg_{mode.get('name', idx)}"
        profiles[pname] = sampled_profile(xs, mode["profile_values"])
        bdg_modes.append({"name": str(mode.get("name", idx)), "lambda_B": float(mode["lambda_B"]), "varpi": float(mode["varpi"]), "profile": pname})

    mixed_ports = []
    for idx, port in enumerate(packet.get("mixed_ports", [])):
        uname = f"u_{port.get('name', idx)}"
        wname = f"w_{port.get('name', idx)}"
        profiles[uname] = sampled_profile(xs, port["u_values"])
        profiles[wname] = sampled_profile(xs, port["w_values"])
        mixed_ports.append({
            "name": str(port.get("name", idx)),
            "lambda_U": float(port["lambda_U"]),
            "Omega_U": float(port["Omega_U"]),
            "u_profile": uname,
            "lambda_W": float(port["lambda_W"]),
            "Omega_W": float(port["Omega_W"]),
            "w_profile": wname,
            "lambda_R": float(port["lambda_R"]),
            "u_w_profiles": [uname, wname],
        })

    reduction: Dict[str, Any] = {
        "wall": {"K": float(packet["wall"]["K"]), "M": float(packet["wall"]["M"])},
        "bdg_modes": bdg_modes,
        "mixed_ports": mixed_ports,
    }
    if isinstance(packet.get("derived_maxwell_transfer"), Mapping) and not mixed_ports:
        reduction["direct_coefficients"] = direct_coefficients_from_packet(packet)
        reduction["derived_coefficients_provenance"] = {
            "bdg_wall": packet.get("derived_bdg_wall_coefficients", {}),
            "maxwell_transfer": packet.get("derived_maxwell_transfer", {}),
        }

    geometry = dict(packet["geometry"])
    geometry_for_adapter = {"R_exit": geometry["R_exit"], "boundary_class": geometry["boundary_class"], "Y_L_limit": geometry.get("Y_L_limit", 0.0), "L": geometry["L"]}
    branch: Dict[str, Any] = {
        "name": str(packet["branch_id"]),
        "geometry": geometry_for_adapter,
        "target": {"constants": dict(packet["constants"])},
        "profiles": profiles,
        "integration": {"grid_points": int(packet.get("solver_metadata", {}).get("mesh_points", len(xs)))},
        "reduction": reduction,
        "handoff_freeze": {
            "source_schema": packet["schema"],
            "source_packet_hash": stable_hash(packet),
            "validator_schema": "stage_v2_22b_solver_handoff/v1",
            "freeze": dict(packet["freeze"]),
        },
    }
    if packet.get("weak_axisymmetric"):
        branch["weak_axisymmetric"] = packet["weak_axisymmetric"]

    return {
        "schema": "stage_v2_22a_profile_adapter/v1",
        "description": "Generated by V2-22B solver-handoff validator from a frozen PDE solver export.",
        "source_solver_packet_hash": stable_hash(packet),
        "branches": [branch],
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def load_packet(path: Optional[str]) -> Dict[str, Any]:
    if path is None:
        return build_sample_solver_output(valid=True)
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage V2-22B solver handoff validator")
    parser.add_argument("--solver-output", default=None, help="Solver output JSON to validate; defaults to built-in valid sample")
    parser.add_argument("--write-schema", default=None, help="Write the frozen solver-output schema document")
    parser.add_argument("--write-valid", default=None, help="Write built-in valid sample solver output")
    parser.add_argument("--write-invalid", default=None, help="Write built-in invalid hard-cap sample solver output")
    parser.add_argument("--out-profile-manifest", default=None, help="Write V2-22A-compatible profile manifest if validation passes")
    parser.add_argument("--out-report", default=None, help="Write validation report JSON")
    args = parser.parse_args()

    if args.write_schema:
        with open(args.write_schema, "w", encoding="utf-8") as f:
            json.dump(solver_output_schema(), f, indent=2, sort_keys=True)
    if args.write_valid:
        with open(args.write_valid, "w", encoding="utf-8") as f:
            json.dump(build_sample_solver_output(valid=True), f, indent=2, sort_keys=True)
    if args.write_invalid:
        with open(args.write_invalid, "w", encoding="utf-8") as f:
            json.dump(build_sample_solver_output(valid=False), f, indent=2, sort_keys=True)

    packet = load_packet(args.solver_output)
    symbolic = run_symbolic_audit()
    report = validate_solver_output(packet)
    report["symbolic_audit"] = symbolic

    profile_manifest = None
    if report["validation_pass"]:
        profile_manifest = convert_solver_output_to_v22a_profile_manifest(packet)
        report["generated_profile_manifest_hash"] = stable_hash(profile_manifest)

    if args.out_report:
        with open(args.out_report, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, sort_keys=True)
    if args.out_profile_manifest and profile_manifest is not None:
        with open(args.out_profile_manifest, "w", encoding="utf-8") as f:
            json.dump(profile_manifest, f, indent=2, sort_keys=True)

    print("Stage V2-22B solver handoff validator")
    print(f"symbolic_checks: {symbolic['checks_passed']}/{symbolic['checks_total']} passed")
    print(f"branch_id: {report.get('branch_id')}")
    print(f"packet_hash: {report['packet_hash']}")
    print(f"validation_pass: {report['validation_pass']}")
    print(f"error_count: {report['error_count']}")
    print(f"warning_count: {report['warning_count']}")
    if profile_manifest is not None:
        print(f"generated_profile_manifest_hash: {stable_hash(profile_manifest)}")
    else:
        print("generated_profile_manifest_hash: <not generated due to validation failure>")
    if report["issues"]:
        print("issues:")
        for item in report["issues"]:
            print(f"  - {item['severity']} {item['path']}: {item['message']}")
    else:
        print("issues: none")
    return 0 if symbolic["checks_failed"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
