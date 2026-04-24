#!/usr/bin/env python3
"""Build a combined referee-facing PDE audit summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

from write_summary_json import parse_output


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def main() -> int:
    parser = argparse.ArgumentParser(description="Write combined PDE audit referee summary")
    parser.add_argument("--audit-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    audit_dir = Path(args.audit_dir)
    output_dir = Path(args.output_dir)
    top_outputs: List[Dict[str, Any]] = [
        parse_output(path)
        for path in sorted(output_dir.glob("*.txt"))
        if path.name != "_summary.txt"
    ]

    python_summary = load_json(audit_dir / "scripts" / "output" / "_summary.json")
    simulation_summary = load_json(audit_dir / "simulation" / "output" / "_summary.json")
    mathematica_summary = load_json(audit_dir / "mathematica" / "output" / "_summary.json")

    total_checks = len(top_outputs)
    top_passed = sum(1 for item in top_outputs if item["status"] == "PASS")
    top_failed = total_checks - top_passed

    root_json = sorted(path.name for path in audit_dir.glob("*.json"))
    referee_pass = (
        top_failed == 0
        and python_summary.get("failed", 1) == 0
        and simulation_summary.get("failed", 1) == 0
        and mathematica_summary.get("failed", 1) == 0
        and not root_json
    )

    combined = {
        "schema": "pde_audit_referee_summary/v1",
        "audit_dir": str(audit_dir),
        "referee_pass": referee_pass,
        "root_json_files": root_json,
        "top_level_checks": top_outputs,
        "python_summary": {
            "path": str(audit_dir / "scripts" / "output" / "_summary.json"),
            "total": python_summary.get("total"),
            "passed": python_summary.get("passed"),
            "failed": python_summary.get("failed"),
        },
        "python_environment": python_summary.get("environment"),
        "simulation_summary": {
            "path": str(audit_dir / "simulation" / "output" / "_summary.json"),
            "total": simulation_summary.get("total"),
            "passed": simulation_summary.get("passed"),
            "failed": simulation_summary.get("failed"),
            "candidate_count": simulation_summary.get("evaluation", {}).get("candidate_count"),
            "validation_pass_count": simulation_summary.get("evaluation", {}).get("validation_pass_count"),
            "open_stable_count": simulation_summary.get("evaluation", {}).get("open_stable_count"),
            "target_pass_count": simulation_summary.get("evaluation", {}).get("target_pass_count"),
            "best_post_hoc_candidate": simulation_summary.get("evaluation", {}).get("best_post_hoc_candidate"),
            "best_post_hoc_score": simulation_summary.get("evaluation", {}).get("best_post_hoc_score"),
            "score_distribution": simulation_summary.get("evaluation", {}).get("score_distribution"),
            "dominant_score_component_counts": simulation_summary.get("evaluation", {}).get("dominant_score_component_counts"),
            "pass_flag_counts": simulation_summary.get("evaluation", {}).get("pass_flag_counts"),
            "reduced_fem_verification": simulation_summary.get("reduced_fem_verification"),
            "nonlinear_readiness": simulation_summary.get("nonlinear_readiness"),
            "physical_model_status": simulation_summary.get("physical_model_status"),
            "nonlinear_export": simulation_summary.get("nonlinear_export"),
            "target_blind_guard": simulation_summary.get("target_blind_guard"),
            "obstruction": simulation_summary.get("obstruction"),
            "required_deformation": simulation_summary.get("required_deformation"),
            "mechanism_gap": simulation_summary.get("mechanism_gap"),
        },
        "mathematica_summary": {
            "path": str(audit_dir / "mathematica" / "output" / "_summary.json"),
            "total": mathematica_summary.get("total"),
            "passed": mathematica_summary.get("passed"),
            "failed": mathematica_summary.get("failed"),
        },
        "mathematica_environment": mathematica_summary.get("environment"),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "_summary.json").open("w", encoding="utf-8") as f:
        json.dump(combined, f, indent=2, sort_keys=True)
        f.write("\n")

    with (output_dir / "_summary.txt").open("w", encoding="utf-8") as f:
        for item in top_outputs:
            f.write(f"{item['status']}  {item['name']}  ({item.get('elapsed', 'unknown')})\n")
        f.write("\n")
        f.write(
            "PYTHON: "
            f"{python_summary.get('passed')}/{python_summary.get('total')} passed, "
            f"{python_summary.get('failed')} failed\n"
        )
        f.write(
            "SIMULATION: "
            f"{simulation_summary.get('passed')}/{simulation_summary.get('total')} passed, "
            f"{simulation_summary.get('failed')} failed; "
            f"{simulation_summary.get('evaluation', {}).get('target_pass_count')}/"
            f"{simulation_summary.get('evaluation', {}).get('candidate_count')} target-passing frozen candidates\n"
        )
        score_dist = simulation_summary.get("evaluation", {}).get("score_distribution") or {}
        if score_dist:
            f.write(
                "SIMULATION_SCORE: "
                f"min {score_dist.get('min')}, "
                f"median {score_dist.get('median')}, "
                f"max {score_dist.get('max')}\n"
            )
        reduced_fem = simulation_summary.get("reduced_fem_verification") or {}
        if reduced_fem:
            f.write(
                "SIMULATION_REDUCED_FEM: "
                f"pass {reduced_fem.get('pass')}, "
                f"failed {len(reduced_fem.get('failed_checks') or [])}\n"
            )
        nonlinear = simulation_summary.get("nonlinear_readiness") or {}
        if nonlinear:
            f.write(
                "SIMULATION_NONLINEAR_READINESS: "
                f"pass {nonlinear.get('pass')}, "
                f"failed {len(nonlinear.get('failed_checks') or [])}, "
                f"packets_emitted {nonlinear.get('packets_emitted')}\n"
            )
        physical_model = simulation_summary.get("physical_model_status") or {}
        if physical_model:
            f.write(
                "SIMULATION_PHYSICAL_MODEL: "
                f"pass {physical_model.get('pass')}, "
                f"export_permitted {physical_model.get('physical_export_permitted')}, "
                f"blockers {physical_model.get('blocker_count')}, "
                f"packets_emitted {physical_model.get('packets_emitted')}\n"
            )
        nonlinear_export = simulation_summary.get("nonlinear_export") or {}
        if nonlinear_export:
            f.write(
                "SIMULATION_NONLINEAR_EXPORT: "
                f"{nonlinear_export.get('target_pass_count')}/"
                f"{nonlinear_export.get('candidate_count')} target-passing frozen candidates, "
                f"guard {((nonlinear_export.get('target_blind_guard') or {}).get('pass'))}\n"
            )
        target_blind = simulation_summary.get("target_blind_guard") or {}
        if target_blind:
            f.write(
                "SIMULATION_TARGET_BLIND_GUARD: "
                f"pass {target_blind.get('pass')}, "
                f"issues {target_blind.get('issue_count')}, "
                f"packets {target_blind.get('packet_count')}\n"
            )
        obstruction = simulation_summary.get("obstruction") or {}
        ratio_dist = obstruction.get("open_stable_one_pole_ratio_distribution") or {}
        if ratio_dist:
            f.write(
                "SIMULATION_OPEN_STABLE_ONE_POLE_RATIO: "
                f"min {ratio_dist.get('min')}, "
                f"median {ratio_dist.get('median')}, "
                f"max {ratio_dist.get('max')}\n"
            )
        open_stable_near = obstruction.get("open_stable_near_one_pole_count_abs_ratio_minus_1_lt_0p1")
        if open_stable_near is not None:
            f.write(
                "SIMULATION_OPEN_STABLE_NEAR_ONE_POLE: "
                f"{open_stable_near} candidates within |ratio - 1| < 0.1\n"
            )
        deformation = simulation_summary.get("required_deformation") or {}
        c_mult_dist = deformation.get("open_stable_required_C_multiplier_distribution") or {}
        p0_mult_dist = deformation.get("open_stable_required_P0_multiplier_distribution") or {}
        if c_mult_dist:
            f.write(
                "SIMULATION_REQUIRED_DEFORMATION: "
                f"C_multiplier_min {c_mult_dist.get('min')}, "
                f"C_multiplier_median {c_mult_dist.get('median')}, "
                f"P0_multiplier_median {p0_mult_dist.get('median')}, "
                f"local_continuation {deformation.get('local_continuation_looks_promising')}\n"
            )
        nonlinear_req = (simulation_summary.get("nonlinear_export") or {}).get("required_deformation") or {}
        nonlinear_c_mult_dist = nonlinear_req.get("open_stable_required_C_multiplier_distribution") or {}
        if nonlinear_c_mult_dist:
            f.write(
                "SIMULATION_NONLINEAR_REQUIRED_DEFORMATION: "
                f"C_multiplier_min {nonlinear_c_mult_dist.get('min')}, "
                f"C_multiplier_median {nonlinear_c_mult_dist.get('median')}, "
                f"local_continuation {nonlinear_req.get('local_continuation_looks_promising')}\n"
            )
        mechanism_gap = simulation_summary.get("mechanism_gap") or {}
        if mechanism_gap:
            f.write(
                "SIMULATION_MECHANISM_GAP: "
                f"pass {mechanism_gap.get('pass')}, "
                f"primary {mechanism_gap.get('primary_missing_channel')}, "
                f"reduced_gap {mechanism_gap.get('reduced_gap_class')}, "
                f"nonlinear_gap {mechanism_gap.get('nonlinear_gap_class')}, "
                f"local_continuation {mechanism_gap.get('local_continuation_recommended')}\n"
            )
        f.write(
            "MATHEMATICA: "
            f"{mathematica_summary.get('passed')}/{mathematica_summary.get('total')} passed, "
            f"{mathematica_summary.get('failed')} failed\n"
        )
        f.write(f"ROOT_JSON_FILES: {len(root_json)}\n")
        f.write(f"REFEREE_PASS: {referee_pass}\n")

    print(f"Wrote {output_dir / '_summary.json'}")
    print(f"REFEREE_PASS: {referee_pass}")
    return 0 if referee_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
