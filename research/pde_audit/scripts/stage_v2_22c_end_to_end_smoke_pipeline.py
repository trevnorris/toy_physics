#!/usr/bin/env python3
"""
Stage V2-22C — End-to-end moving-throat branch-realization smoke pipeline.

This orchestrates the three executable handoff layers already built in Volume 2:

  V2-22B solver handoff validator
      frozen PDE/solver export -> validated open-throat branch packet

  V2-22A profile-to-coefficient adapter
      sampled/analytic profiles -> V2-21 branch manifest

  V2-21 branch extraction fixture
      branch manifest -> grouped-P2 observable packet and target residuals

The point of this script is not to claim that the sample branch realizes the GR
normalization target.  The built-in open-throat D/N sample is deliberately a
stable but uncalibrated smoke branch.  The point is to prove that the whole
branch-realization pipeline is now executable, target-blind, hashable, and able
to distinguish:

  * validation / handoff failures,
  * open-throat and stability failures,
  * and genuine target residual failures.

Run:
  python stage_v2_22c_end_to_end_smoke_pipeline.py \
    --out-report stage_v2_22c_pipeline_report.json \
    --out-profile-manifest stage_v2_22c_generated_profile_manifest.json \
    --out-v21-manifest stage_v2_22c_generated_v21_manifest.json \
    --out-observable-packet stage_v2_22c_observable_packet.json \
    --out-tolerance-budget stage_v2_22c_tolerance_budget.json
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import pathlib
import sys
import hashlib
from typing import Any, Dict, Mapping, Optional, Tuple

HERE = pathlib.Path(__file__).resolve().parent
LANES = ("20", "21", "22")


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def import_module_from_path(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_json(path: Optional[str]) -> Optional[Dict[str, Any]]:
    if path is None:
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Optional[str], payload: Any) -> None:
    if path is None:
        return
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def run_local_symbolic_audit() -> Dict[str, Any]:
    """Lightweight exact/numeric formula checks for the orchestration layer.

    This intentionally avoids rerunning the heavier component SymPy audits so
    the smoke pipeline remains fast and focused on handoff closure.
    """
    checks = []

    # Grouped inverse map on a representative exact triple.
    x20, x21, x22 = 7.0, -3.0, 11.0
    xbar = (x20 + 2*x21 + 2*x22)/5
    ax = (2*x20 - x21 - x22)/10
    bx = (x21 - x22)/2
    checks.append(("group_inverse_20", abs((xbar + 4*ax) - x20) < 1e-14))
    checks.append(("group_inverse_21", abs((xbar - ax + bx) - x21) < 1e-14))
    checks.append(("group_inverse_22", abs((xbar - ax - bx) - x22) < 1e-14))

    # Weak-axisymmetric signature.
    x0, eps, x1 = 2.0, 0.013, 5.0
    y20 = x0 + eps*x1
    y21 = x0 + eps*0.5*x1
    y22 = x0 - eps*x1
    trace = (y20 + 2*y21 + 2*y22)/5
    ay = (2*y20 - y21 - y22)/10
    by = (y21 - y22)/2
    checks.append(("axisym_trace_unchanged", abs(trace - x0) < 1e-14))
    checks.append(("axisym_b_equals_3a", abs(by - 3*ay) < 1e-14))

    # Constant-prefactor branch formulas on representative coefficients.
    D0, D2, D4, N0 = 1.7, -0.4, -0.9, 3.2
    N2_const = 2*D2*N0/D0
    N4_const = N0*(D2**2 + 2*D0*D4)/D0**2
    P2 = (D0*N2_const - 2*D2*N0)/D0**2
    P4 = (D0**2*N4_const - 2*D0*(D2*N2_const + D4*N0) + 3*D2**2*N0)/D0**3
    checks.append(("constant_prefactor_P2_zero", abs(P2) < 1e-14))
    checks.append(("constant_prefactor_P4_zero", abs(P4) < 1e-14))

    # Target equivalence in normalized constants.
    G = cs = a = c = mhat = Sport = 1.0
    Ptarget = 54*G*cs**5/(5*a**5*c**5)
    gamma_eff = mhat**2*Sport*Ptarget*a**5/(27*cs**5)
    gamma_gr = Sport*mhat**2*2*G/(5*c**5)
    checks.append(("target_gamma_equivalence", abs(gamma_eff - gamma_gr) < 1e-14))

    passed = sum(1 for _, ok in checks if ok)
    return {
        "checks_total": len(checks),
        "checks_passed": passed,
        "checks_failed": len(checks) - passed,
        "details": [{"name": name, "pass": bool(ok)} for name, ok in checks],
    }


def audit_summary(name: str, audit: Mapping[str, Any]) -> Dict[str, Any]:
    return {
        "name": name,
        "checks_total": int(audit.get("checks_total", 0)),
        "checks_passed": int(audit.get("checks_passed", 0)),
        "checks_failed": int(audit.get("checks_failed", 0)),
        "pass": int(audit.get("checks_failed", 0)) == 0,
    }


def residual_budget_for_branch(branch_packet: Mapping[str, Any], tol: float) -> Dict[str, Any]:
    residuals = branch_packet["residuals"]
    target = branch_packet["target_values"]
    grouped = branch_packet["grouped"]

    absolute = {k: float(v) for k, v in residuals.items()}
    scaled = {
        "R_pole_over_tol": abs(float(residuals["R_pole"])) / tol,
        "R_norm_over_P0_target": abs(float(residuals["R_norm"])) / max(abs(float(target["P0_target"])), tol),
        "R_norm_over_tol": abs(float(residuals["R_norm"])) / tol,
        "R_P2_over_tol": abs(float(residuals["R_P2"])) / tol,
        "R_P4_over_tol": abs(float(residuals["R_P4"])) / tol,
        "R_tail_over_tol": abs(float(residuals["R_tail"])) / tol,
        "P0_anisotropy_norm_sq_over_tol": abs(float(grouped["P0"]["anisotropy_norm_sq"])) / tol,
        "D0_anisotropy_norm_sq_over_tol": abs(float(grouped["D0"]["anisotropy_norm_sq"])) / tol,
    }
    return {
        "branch": branch_packet["name"],
        "tolerance": tol,
        "absolute_residuals": absolute,
        "scaled_residuals": scaled,
        "pass_flags": branch_packet["pass_flags"],
        "target_packet_pass": bool(branch_packet["pass_flags"].get("target_packet_pass", False)),
        "open_gate_pass": bool(branch_packet["pass_flags"].get("open_gate_pass", False)),
        "stability_gate_pass": bool(branch_packet["pass_flags"].get("stability_gate_pass", False)),
    }


def build_pipeline_report(
    v21,
    v22a,
    v22b,
    solver_packet: Mapping[str, Any],
    invalid_solver_packet: Mapping[str, Any],
    tol: float,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Run the full V2-22B -> V2-22A -> V2-21 chain."""

    # 1. Lightweight symbolic audit for the orchestration layer.
    # The component tools were audited in their own stages; rerunning those
    # full SymPy audits here is intentionally avoided so this smoke test remains
    # a fast handoff pipeline rather than another symbolic compiler.
    local_audit = run_local_symbolic_audit()
    formula_audits = [audit_summary("V2-22C orchestration identities", local_audit)]
    formula_pass = all(a["pass"] for a in formula_audits)
    component_hashes = {
        "V2-21_branch_extraction_script": stable_hash(open(HERE / "stage_v2_21_branch_extraction_fixture.py", "r", encoding="utf-8").read()),
        "V2-22A_profile_adapter_script": stable_hash(open(HERE / "stage_v2_22a_profile_to_coefficient_adapter.py", "r", encoding="utf-8").read()),
        "V2-22B_solver_handoff_script": stable_hash(open(HERE / "stage_v2_22b_solver_handoff_validator.py", "r", encoding="utf-8").read()),
    }

    # 2. Validate valid and invalid solver exports.
    validation_report = v22b.validate_solver_output(solver_packet)
    invalid_validation_report = v22b.validate_solver_output(invalid_solver_packet)
    hardcap_rejection_pass = (invalid_validation_report["validation_pass"] is False) and invalid_validation_report["error_count"] > 0

    # 3. Convert valid solver output to profile manifest and V2-21 manifest.
    if not validation_report["validation_pass"]:
        raise RuntimeError("Valid solver packet failed validation; cannot run smoke pipeline")
    profile_manifest = v22b.convert_solver_output_to_v22a_profile_manifest(solver_packet)
    v21_manifest, adapter_diag = v22a.adapt_manifest(profile_manifest)

    # 4. Extract observable packet with the V2-21 fixture.
    branch_packets = [v21.extract_branch(branch, tol=tol) for branch in v21_manifest.get("branches", [])]
    observable_packet = {
        "schema": "stage_v2_22c_observable_packet/v1",
        "source_v21_manifest_hash": stable_hash(v21_manifest),
        "branch_count": len(branch_packets),
        "branches": branch_packets,
    }

    # 5. Include a calibrated control to prove target pass/fail discrimination.
    calibrated_branch = v21.calibrated_coefficient_branch()
    calibrated_packet = v21.extract_branch(calibrated_branch, tol=tol)

    # 6. Gate accounting.
    smoke_branch = branch_packets[0] if branch_packets else None
    smoke_open_and_stable = bool(smoke_branch and smoke_branch["pass_flags"].get("open_gate_pass") and smoke_branch["pass_flags"].get("stability_gate_pass"))
    smoke_target_pass = bool(smoke_branch and smoke_branch["pass_flags"].get("target_packet_pass"))
    calibrated_target_pass = bool(calibrated_packet["pass_flags"].get("target_packet_pass"))

    mechanical_pipeline_pass = all([
        formula_pass,
        bool(validation_report["validation_pass"]),
        hardcap_rejection_pass,
        len(branch_packets) == 1,
        smoke_open_and_stable,
        calibrated_target_pass,
    ])

    tolerance_budget = {
        "schema": "stage_v2_22c_tolerance_budget/v1",
        "extraction_residual_tol": tol,
        "validator_profile_norm_tol": float(solver_packet.get("normalization_tolerances", {}).get("profile_norm_tol", 5e-3)) if isinstance(solver_packet.get("normalization_tolerances", {}), Mapping) else 5e-3,
        "validator_Delta_tol": float(solver_packet.get("normalization_tolerances", {}).get("Delta_tol", 1e-12)) if isinstance(solver_packet.get("normalization_tolerances", {}), Mapping) else 1e-12,
        "interpretation": {
            "mechanical_pipeline_pass": "Validation, conversion, extraction, open-throat gates, stability gates, calibrated-control target pass, and invalid hard-cap rejection all passed.",
            "branch_target_packet_pass": "Whether the actual tested solver branch hits the one-pole, normalization, constant-prefactor, and tail residuals within extraction_residual_tol.",
            "target_failure_is_not_pipeline_failure": "A stable open branch may legitimately fail target residuals; that is the point of the extraction fixture.",
        },
        "smoke_branch_budget": residual_budget_for_branch(smoke_branch, tol) if smoke_branch else None,
        "calibrated_control_budget": residual_budget_for_branch(calibrated_packet, tol),
    }

    report = {
        "schema": "stage_v2_22c_pipeline_report/v1",
        "stage": "V2-22C",
        "purpose": "End-to-end target-blind smoke pipeline from solver export to grouped-P2 observable packet.",
        "component_audits": formula_audits,
        "component_script_hashes": component_hashes,
        "formula_pass": formula_pass,
        "valid_solver_validation": validation_report,
        "invalid_hardcap_validation": invalid_validation_report,
        "hardcap_rejection_pass": hardcap_rejection_pass,
        "profile_manifest_hash": stable_hash(profile_manifest),
        "v21_manifest_hash": stable_hash(v21_manifest),
        "adapter_diagnostics": adapter_diag,
        "observable_packet_hash": stable_hash(observable_packet),
        "calibrated_control": calibrated_packet,
        "smoke_branch_summary": {
            "branch_name": smoke_branch["name"] if smoke_branch else None,
            "open_gate_pass": smoke_branch["pass_flags"]["open_gate_pass"] if smoke_branch else False,
            "stability_gate_pass": smoke_branch["pass_flags"]["stability_gate_pass"] if smoke_branch else False,
            "target_packet_pass": smoke_target_pass,
            "D0_bar": smoke_branch["grouped"]["D0"]["bar"] if smoke_branch else math.nan,
            "N0_bar": smoke_branch["grouped"]["N0"]["bar"] if smoke_branch else math.nan,
            "P0_bar": smoke_branch["grouped"]["P0"]["bar"] if smoke_branch else math.nan,
            "P0_target": smoke_branch["target_values"]["P0_target"] if smoke_branch else math.nan,
            "R_pole": smoke_branch["residuals"]["R_pole"] if smoke_branch else math.nan,
            "R_norm": smoke_branch["residuals"]["R_norm"] if smoke_branch else math.nan,
            "R_P2": smoke_branch["residuals"]["R_P2"] if smoke_branch else math.nan,
            "R_P4": smoke_branch["residuals"]["R_P4"] if smoke_branch else math.nan,
            "R_tail": smoke_branch["residuals"]["R_tail"] if smoke_branch else math.nan,
        },
        "mechanical_pipeline_pass": mechanical_pipeline_pass,
        "branch_target_realization_claimed": False,
        "notes": [
            "The smoke branch is expected to be open/stable but target-failing unless an actual calibrated PDE branch is supplied.",
            "The calibrated direct-coefficient control exists only to prove that the target gates can pass when the algebraic surface is satisfied.",
            "The hard-cap sample is intentionally rejected to enforce the V2-04 open-organ-pipe patch.",
        ],
        "hashes": {
            "solver_packet_hash": stable_hash(solver_packet),
            "invalid_solver_packet_hash": stable_hash(invalid_solver_packet),
            "profile_manifest_hash": stable_hash(profile_manifest),
            "v21_manifest_hash": stable_hash(v21_manifest),
            "observable_packet_hash": stable_hash(observable_packet),
            "tolerance_budget_hash": stable_hash(tolerance_budget),
        },
    }

    return report, profile_manifest, v21_manifest, observable_packet, tolerance_budget


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage V2-22C end-to-end smoke pipeline")
    parser.add_argument("--solver-output", default=None, help="Optional valid solver export JSON. Defaults to built-in V2-22B open-throat sample.")
    parser.add_argument("--invalid-solver-output", default=None, help="Optional invalid solver export JSON. Defaults to built-in V2-22B hard-cap sample.")
    parser.add_argument("--tol", type=float, default=1e-9, help="Extraction residual tolerance")
    parser.add_argument("--out-report", default=None, help="Write full pipeline report JSON")
    parser.add_argument("--out-profile-manifest", default=None, help="Write generated V2-22A profile manifest JSON")
    parser.add_argument("--out-v21-manifest", default=None, help="Write generated V2-21 branch manifest JSON")
    parser.add_argument("--out-observable-packet", default=None, help="Write grouped observable packet JSON")
    parser.add_argument("--out-tolerance-budget", default=None, help="Write tolerance budget JSON")
    parser.add_argument("--write-valid-solver", default=None, help="Write the valid solver packet actually used")
    parser.add_argument("--write-invalid-solver", default=None, help="Write the invalid hard-cap solver packet actually used")
    args = parser.parse_args()

    v21 = import_module_from_path("stage_v2_21", HERE / "stage_v2_21_branch_extraction_fixture.py")
    v22a = import_module_from_path("stage_v2_22a", HERE / "stage_v2_22a_profile_to_coefficient_adapter.py")
    v22b = import_module_from_path("stage_v2_22b", HERE / "stage_v2_22b_solver_handoff_validator.py")

    solver_packet = load_json(args.solver_output) or v22b.build_sample_solver_output(valid=True)
    invalid_solver_packet = load_json(args.invalid_solver_output) or v22b.build_sample_solver_output(valid=False)

    write_json(args.write_valid_solver, solver_packet)
    write_json(args.write_invalid_solver, invalid_solver_packet)

    report, profile_manifest, v21_manifest, observable_packet, tolerance_budget = build_pipeline_report(
        v21=v21,
        v22a=v22a,
        v22b=v22b,
        solver_packet=solver_packet,
        invalid_solver_packet=invalid_solver_packet,
        tol=args.tol,
    )

    write_json(args.out_report, report)
    write_json(args.out_profile_manifest, profile_manifest)
    write_json(args.out_v21_manifest, v21_manifest)
    write_json(args.out_observable_packet, observable_packet)
    write_json(args.out_tolerance_budget, tolerance_budget)

    print("STAGE V2-22C END-TO-END SMOKE PIPELINE")
    for audit in report["component_audits"]:
        print(f"{audit['name']}: {audit['checks_passed']}/{audit['checks_total']} formula checks passed")
    print(f"valid_solver_validation_pass: {report['valid_solver_validation']['validation_pass']}")
    print(f"invalid_hardcap_rejection_pass: {report['hardcap_rejection_pass']}")
    print(f"profile_manifest_hash: {report['profile_manifest_hash']}")
    print(f"v21_manifest_hash: {report['v21_manifest_hash']}")
    print(f"observable_packet_hash: {report['observable_packet_hash']}")
    s = report["smoke_branch_summary"]
    print("--- smoke branch ---")
    print(f"branch: {s['branch_name']}")
    print(f"open_gate_pass: {s['open_gate_pass']}")
    print(f"stability_gate_pass: {s['stability_gate_pass']}")
    print(f"target_packet_pass: {s['target_packet_pass']}")
    print(f"D0_bar: {s['D0_bar']:.16g}")
    print(f"N0_bar: {s['N0_bar']:.16g}")
    print(f"P0_bar: {s['P0_bar']:.16g}")
    print(f"P0_target: {s['P0_target']:.16g}")
    print(f"R_pole: {s['R_pole']:.16g}")
    print(f"R_norm: {s['R_norm']:.16g}")
    print(f"R_P2: {s['R_P2']:.16g}")
    print(f"R_P4: {s['R_P4']:.16g}")
    print(f"R_tail: {s['R_tail']:.16g}")
    c = report["calibrated_control"]
    print("--- calibrated direct-coefficient control ---")
    print(f"target_packet_pass: {c['pass_flags']['target_packet_pass']}")
    print(f"R_pole: {c['residuals']['R_pole']:.16g}")
    print(f"R_norm: {c['residuals']['R_norm']:.16g}")
    print(f"R_P2: {c['residuals']['R_P2']:.16g}")
    print(f"R_P4: {c['residuals']['R_P4']:.16g}")
    print(f"mechanical_pipeline_pass: {report['mechanical_pipeline_pass']}")
    print(f"branch_target_realization_claimed: {report['branch_target_realization_claimed']}")
    return 0 if report["mechanical_pipeline_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
