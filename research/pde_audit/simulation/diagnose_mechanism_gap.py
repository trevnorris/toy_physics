#!/usr/bin/env python3
"""Summarize the mechanism gap exposed by frozen simulation candidates.

This is a post-hoc diagnostic.  It reads already-generated evaluation,
deformation, and physical-model status reports.  It does not generate,
mutate, rank-for-generation, or refit solver packets.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Mapping


OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def dist_value(block: Mapping[str, Any], name: str, field: str) -> Any:
    return ((block.get("open_stable_candidates") or {}).get(name) or {}).get(field)


def lane_gap(name: str, deformation: Mapping[str, Any]) -> Dict[str, Any]:
    one_pole_max = dist_value(deformation, "one_pole_ratio_distribution", "max")
    c_min = dist_value(deformation, "required_C_multiplier_distribution", "min")
    c_median = dist_value(deformation, "required_C_multiplier_distribution", "median")
    p0_median = dist_value(deformation, "required_P0_multiplier_distribution", "median")
    local_continuation = (deformation.get("mechanism_assessment") or {}).get("local_continuation_looks_promising")
    dominant_counts = deformation.get("dominant_score_component_counts") or {}
    return {
        "name": name,
        "candidate_count": deformation.get("candidate_count"),
        "target_pass_count": deformation.get("target_pass_count"),
        "open_stable_count": deformation.get("open_stable_count"),
        "dominant_score_component_counts": dominant_counts,
        "primary_obstruction": (deformation.get("mechanism_assessment") or {}).get("primary_obstruction"),
        "open_stable_one_pole_ratio_max": one_pole_max,
        "open_stable_required_C_or_D0_multiplier_min": c_min,
        "open_stable_required_C_or_D0_multiplier_median": c_median,
        "open_stable_required_P0_multiplier_median": p0_median,
        "local_continuation_looks_promising": local_continuation,
        "gap_class": classify_gap(one_pole_max, c_min, local_continuation),
    }


def classify_gap(one_pole_max: Any, c_min: Any, local_continuation: Any) -> str:
    if local_continuation is True:
        return "near_surface_local_continuation_candidate"
    if isinstance(one_pole_max, (float, int)) and one_pole_max < 0.2:
        return "large_one_pole_support_deficit"
    if isinstance(c_min, (float, int)) and c_min > 2.0:
        return "moderate_to_large_one_pole_support_deficit"
    return "mixed_or_unclassified_gap"


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize physical mechanism gap from frozen diagnostics")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--report-prefix", default="mechanism_gap_report", help="Output basename without extension")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    paths = {
        "required_deformation": output_dir / "required_deformation_report.json",
        "nonlinear_required_deformation": output_dir / "nonlinear_required_deformation_report.json",
        "physical_model_status": output_dir / "physical_model_status.json",
        "evaluation_summary": output_dir / "evaluation_summary.json",
        "nonlinear_evaluation_summary": output_dir / "nonlinear_evaluation_summary.json",
    }
    missing = [name for name, path in paths.items() if not path.exists()]
    if missing:
        raise RuntimeError(f"missing required input report(s): {', '.join(missing)}")

    deformation = load_json(paths["required_deformation"])
    nonlinear_deformation = load_json(paths["nonlinear_required_deformation"])
    physical = load_json(paths["physical_model_status"])
    evaluation = load_json(paths["evaluation_summary"])
    nonlinear_evaluation = load_json(paths["nonlinear_evaluation_summary"])

    reduced_gap = lane_gap("reduced_operator_v1", deformation)
    nonlinear_gap = lane_gap("nonlinear_manufactured_readiness", nonlinear_deformation)
    blocker_ids = [
        item.get("id")
        for item in (physical.get("physical_model_status") or {}).get("blocking_reasons", [])
    ]
    physical_export_permitted = bool(physical.get("physical_export_permitted"))
    target_pass_count = int(evaluation.get("target_pass_count") or 0)
    nonlinear_target_pass_count = int(nonlinear_evaluation.get("target_pass_count") or 0)

    report = {
        "schema": "pde_audit_simulation_mechanism_gap_report/v1",
        "post_hoc_only": True,
        "target_residuals_used_to_generate_candidates": False,
        "candidate_generation_mutated": False,
        "input_reports": {
            name: {"path": str(path), "sha256": sha256_file(path)}
            for name, path in paths.items()
        },
        "overall": {
            "referee_harness_interpretation": "The algebraic target checks and handoff pipeline pass, but the frozen simulated families do not realize the target surface.",
            "reduced_target_pass_count": target_pass_count,
            "nonlinear_target_pass_count": nonlinear_target_pass_count,
            "any_target_pass": bool(target_pass_count or nonlinear_target_pass_count),
            "physical_export_permitted": physical_export_permitted,
            "physical_export_blocker_count": physical.get("blocker_count"),
            "physical_export_blockers": blocker_ids,
            "local_continuation_recommended": bool(
                reduced_gap["local_continuation_looks_promising"]
                or nonlinear_gap["local_continuation_looks_promising"]
            ),
        },
        "reduced_gap": reduced_gap,
        "nonlinear_manufactured_gap": nonlinear_gap,
        "mechanism_conclusion": {
            "primary_missing_channel": "one_pole_support",
            "coefficient_surface": "D0*C/(3*A^2)=1, with A=M+B2+Z2 and C=B4+Z4.",
            "why_algebra_can_pass_while_simulation_fails": (
                "The symbolic audit verifies the target surface and adapters. "
                "The simulation asks whether frozen branch families land on that surface. "
                "The current families under-supply the one-pole C/D0 channel."
            ),
            "why_not_retune_current_candidates": (
                "The reports are post-hoc classifications of frozen target-blind packets; "
                "using them to mutate the same packet family would violate the freeze/no-refit protocol."
            ),
            "next_physical_requirements": [
                "Promote and freeze parent S_eta/S_Sigma dynamics or explicitly declare an effective physical closure.",
                "Specify coupled nonlinear residual equations for wall geometry, stable BdG support, Maxwell/mixed modes, and outgoing port normalization.",
                "Demonstrate a target-blind mechanism that materially increases C or D0 relative to A before target residuals are evaluated.",
                "Keep the existing freeze hash, target-blind guard, evaluation chain, obstruction report, and deformation map unchanged for the next exporter.",
            ],
        },
    }
    report["pass"] = bool(
        report["post_hoc_only"]
        and not report["target_residuals_used_to_generate_candidates"]
        and not report["candidate_generation_mutated"]
        and not physical_export_permitted
    )
    report["report_hash"] = sha256_json(report)
    write_json(output_dir / f"{args.report_prefix}.json", report)

    lines = [
        "PDE audit simulation mechanism-gap report",
        "=" * 52,
        f"pass: {report['pass']}",
        f"reduced_target_pass_count: {target_pass_count}",
        f"nonlinear_target_pass_count: {nonlinear_target_pass_count}",
        f"physical_export_permitted: {physical_export_permitted}",
        f"reduced_gap_class: {reduced_gap['gap_class']}",
        f"reduced_C_or_D0_multiplier_min: {reduced_gap['open_stable_required_C_or_D0_multiplier_min']}",
        f"reduced_C_or_D0_multiplier_median: {reduced_gap['open_stable_required_C_or_D0_multiplier_median']}",
        f"nonlinear_gap_class: {nonlinear_gap['gap_class']}",
        f"nonlinear_C_or_D0_multiplier_min: {nonlinear_gap['open_stable_required_C_or_D0_multiplier_min']}",
        f"local_continuation_recommended: {report['overall']['local_continuation_recommended']}",
        f"primary_missing_channel: {report['mechanism_conclusion']['primary_missing_channel']}",
        f"report_hash: {report['report_hash']}",
    ]
    (output_dir / f"{args.report_prefix}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
