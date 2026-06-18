#!/usr/bin/env python3
"""Stage V2-24 — Negative controls for the PDE V2 audit pipeline.

This audit makes failure behavior executable.  Each control is expected to be
rejected by the relevant pre-target gate or manifest validator; a control that
unexpectedly passes is a hard failure of this audit.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping

import stage_v2_21_branch_extraction_fixture as v21
import stage_v2_22a_profile_to_coefficient_adapter as v22a
import stage_v2_22b_solver_handoff_validator as v22b


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_FIXTURE_DIR = SCRIPT_DIR / "fixtures" / "negative_controls"


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def issue_paths(report: Mapping[str, Any]) -> List[str]:
    return [str(item.get("path", "")) for item in report.get("issues", [])]


def has_path(report: Mapping[str, Any], expected_path: str) -> bool:
    return expected_path in issue_paths(report)


def build_direct_derived_packet() -> Dict[str, Any]:
    packet = v22b.build_sample_solver_output(valid=True)
    packet["branch_id"] = "direct_derived_positive_control"
    packet.pop("mixed_ports", None)
    packet["derived_bdg_wall_coefficients"] = {
        "status": "derived_bdg_wall_m1b_fixture_control",
        "coefficients": {
            "K": packet["wall"]["K"],
            "M": packet["wall"]["M"],
            "B0": 0.01,
            "B2": 0.02,
            "B4": 0.03,
        },
        "source_hashes": {"fixture_control": "direct-bdg-wall"},
    }
    packet["derived_maxwell_transfer"] = {
        "status": "derived_green_function_transfer",
        "gauge_convention": packet["freeze"]["gauge_convention"],
        "flux_normalization": {
            "Gamma_port": 4.0 / 27.0,
            "formula": "fixture_control_sigma_Q_can",
        },
        "coefficients": {
            "Z0": 0.001,
            "Z2": 0.002,
            "Z4": 0.003,
            "N0": 0.004,
            "N2": 0.005,
            "N4": 0.006,
        },
        "operator_gauge_residual_metrics": {
            "current_frechet_matches_step8c": True,
            "outgoing_flux_positive": True,
            "open_not_hard_cap": True,
            "pure_gauge_zero_physical_transfer": True,
            "basis_invariance": True,
            "v2_09_regression": True,
            "green_residuals_small": True,
            "bdg_residuals_small": True,
            "N0_positive": True,
            "max_green_residual": 1e-12,
            "max_gauge_physical_field_norm": 1e-12,
        },
        "source_hashes": {"fixture_control": "direct-maxwell-transfer"},
    }
    packet["solver_metadata"]["coefficient_family"] = "direct_derived_fixture_control_v1"
    packet["solver_metadata"]["source_commit"] = "fixture-only-direct-derived"
    return packet


def run_v22b_negative_controls(fixture_dir: Path) -> List[Dict[str, Any]]:
    cases = [
        ("nonmonotone_grid", "stage_v2_22b_negative_nonmonotone_grid.json", "grid.points"),
        ("profile_length_mismatch", "stage_v2_22b_negative_profile_length_mismatch.json", "profiles.wall_chi_eta"),
        ("pre_target_freeze_false", "stage_v2_22b_negative_pre_target_freeze_false.json", "freeze.pre_target_freeze"),
        ("target_blind_false", "stage_v2_22b_negative_target_blind_false.json", "freeze.target_blind"),
        ("nonpositive_delta", "stage_v2_22b_negative_nonpositive_delta.json", "mixed_ports[0].Delta_eff"),
        ("missing_gauge_convention", "stage_v2_22b_negative_missing_gauge_convention.json", "freeze.gauge_convention"),
        ("bad_boundary_protocol", "stage_v2_22b_negative_bad_boundary_protocol.json", "freeze.boundary_protocol"),
        ("nonfinite_solver_residual", "stage_v2_22b_negative_nonfinite_solver_residual.json", "solver_metadata.nonlinear_residual_norm"),
        ("target_leakage", "stage_v2_22b_negative_target_leakage.json", "target_values"),
    ]
    results: List[Dict[str, Any]] = []
    for name, filename, expected_path in cases:
        packet = load_json(fixture_dir / filename)
        report = v22b.validate_solver_output(packet)
        passed = (not bool(report.get("validation_pass"))) and has_path(report, expected_path)
        results.append({
            "name": f"v22b_{name}",
            "control_type": "negative",
            "expected_rejection_path": expected_path,
            "pass": bool(passed),
            "validation_pass": bool(report.get("validation_pass")),
            "error_count": int(report.get("error_count", 0)),
            "issue_paths": issue_paths(report),
            "packet_hash": report.get("packet_hash"),
        })
    return results


def run_v22b_direct_derived_controls() -> List[Dict[str, Any]]:
    packet = build_direct_derived_packet()
    report = v22b.validate_solver_output(packet)
    profile_manifest = v22b.convert_solver_output_to_v22a_profile_manifest(packet) if report["validation_pass"] else None
    v21_manifest, _ = v22a.adapt_manifest(profile_manifest) if profile_manifest is not None else (None, None)
    branch_packet = v21.extract_branch(v21_manifest["branches"][0]) if v21_manifest is not None else None
    direct_lane_pass = bool(
        branch_packet
        and branch_packet["pass_flags"]["open_gate_pass"]
        and branch_packet["pass_flags"]["stability_gate_pass"]
        and all("direct_coefficients" in v21_manifest["branches"][0]["lanes"][lane] for lane in v21.LANES)
    )

    leak_packet = copy.deepcopy(packet)
    leak_packet["derived_maxwell_transfer"]["R_norm"] = 0.0
    leak_report = v22b.validate_solver_output(leak_packet)
    leak_pass = (not leak_report["validation_pass"]) and has_path(leak_report, "derived_maxwell_transfer.R_norm")

    return [
        {
            "name": "v22b_direct_derived_path_positive_control",
            "control_type": "positive_guard",
            "pass": bool(report["validation_pass"] and direct_lane_pass),
            "validation_pass": bool(report["validation_pass"]),
            "error_count": int(report.get("error_count", 0)),
            "coefficient_path": "direct_derived_coefficients",
            "open_gate_pass": bool(branch_packet and branch_packet["pass_flags"]["open_gate_pass"]),
            "stability_gate_pass": bool(branch_packet and branch_packet["pass_flags"]["stability_gate_pass"]),
            "packet_hash": report.get("packet_hash"),
        },
        {
            "name": "v22b_direct_derived_target_leakage",
            "control_type": "negative",
            "expected_rejection_path": "derived_maxwell_transfer.R_norm",
            "pass": bool(leak_pass),
            "validation_pass": bool(leak_report.get("validation_pass")),
            "error_count": int(leak_report.get("error_count", 0)),
            "issue_paths": issue_paths(leak_report),
            "packet_hash": leak_report.get("packet_hash"),
        },
    ]


def run_v22a_strict_profile_controls() -> List[Dict[str, Any]]:
    expr_manifest = {
        "schema": "stage_v2_22a_profile_adapter/v1",
        "branches": [
            {
                "name": "expr_profile_control",
                "profiles": {
                    "weight": {"kind": "builtin", "name": "one"},
                    "wall": {"kind": "expr", "expr": "sqrt(2/L)*sin(pi*s/L)"},
                },
            }
        ],
    }

    strict_report = v22a.validate_profile_manifest(expr_manifest, strict_profiles=True)
    permissive_report = v22a.validate_profile_manifest(expr_manifest, strict_profiles=False)
    return [
        {
            "name": "v22a_strict_rejects_expr_profile",
            "control_type": "negative",
            "expected_rejection_path": "branches[0].profiles.wall.kind",
            "pass": (not strict_report["validation_pass"]) and has_path(strict_report, "branches[0].profiles.wall.kind"),
            "validation_pass": strict_report["validation_pass"],
            "error_count": strict_report["error_count"],
            "issue_paths": issue_paths(strict_report),
        },
        {
            "name": "v22a_permissive_keeps_trusted_local_expr_profile",
            "control_type": "positive_guard",
            "pass": bool(permissive_report["validation_pass"]),
            "validation_pass": permissive_report["validation_pass"],
            "error_count": permissive_report["error_count"],
            "issue_paths": issue_paths(permissive_report),
        },
    ]


def run_v21_unstable_direct_control() -> List[Dict[str, Any]]:
    branch = v21.calibrated_coefficient_branch()
    branch["name"] = "negative_unstable_direct_D0_branch"
    for lane in v21.LANES:
        branch["lanes"][lane] = copy.deepcopy(branch["lanes"][lane])
        branch["lanes"][lane]["direct_coefficients"]["K"] = -0.1

    packet = v21.extract_branch(branch)
    passed = (
        packet["pass_flags"]["open_gate_pass"]
        and not packet["pass_flags"]["stability_gate_pass"]
        and not packet["pass_flags"]["target_packet_pass"]
        and packet["grouped"]["D0"]["bar"] < 0.0
    )
    return [
        {
            "name": "v21_unstable_direct_D0_negative_control",
            "control_type": "negative",
            "pass": bool(passed),
            "D0_bar": packet["grouped"]["D0"]["bar"],
            "open_gate_pass": packet["pass_flags"]["open_gate_pass"],
            "stability_gate_pass": packet["pass_flags"]["stability_gate_pass"],
            "target_packet_pass": packet["pass_flags"]["target_packet_pass"],
            "stability": packet["stability"],
        }
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage V2-24 negative controls")
    parser.add_argument("--fixture-dir", default=str(DEFAULT_FIXTURE_DIR), help="Directory containing negative fixture JSON files")
    parser.add_argument("--out-json", default=None, help="Optional JSON report path")
    args = parser.parse_args()

    fixture_dir = Path(args.fixture_dir)
    results: List[Dict[str, Any]] = []
    results.extend(run_v22b_negative_controls(fixture_dir))
    results.extend(run_v22b_direct_derived_controls())
    results.extend(run_v22a_strict_profile_controls())
    results.extend(run_v21_unstable_direct_control())

    pass_count = sum(1 for item in results if item["pass"])
    total = len(results)
    report = {
        "schema": "stage_v2_24_negative_controls_report/v1",
        "fixture_dir": str(fixture_dir),
        "checks_total": total,
        "checks_passed": pass_count,
        "checks_failed": total - pass_count,
        "results": results,
    }

    if args.out_json:
        out_path = Path(args.out_json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    print("Stage V2-24: negative controls")
    print("=" * 48)
    for item in results:
        status = "PASS" if item["pass"] else "FAIL"
        print(f"{status}  {item['name']}")
        if item.get("expected_rejection_path"):
            print(f"  expected_rejection_path: {item['expected_rejection_path']}")
            print(f"  issue_paths: {', '.join(item.get('issue_paths', []))}")
    print("")
    print(f"negative_control_checks: {pass_count}/{total} passed")
    return 0 if pass_count == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
