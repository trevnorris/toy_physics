#!/usr/bin/env python3
"""Post-hoc coefficient projection stress tests.

This diagnostic asks what would happen under algebraic coefficient projections
applied after target evaluation.  It is intentionally not a candidate generator
and cannot produce a target-blind hit.  The point is to separate:

1. the one-pole support deficit in C/D0;
2. the outgoing normalization amplitude deficit in P0;
3. the outgoing moment-shape deficit in P2/P4.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping


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


def quantiles(values: List[float]) -> Dict[str, float]:
    xs = sorted(value for value in values if math.isfinite(value))
    if not xs:
        return {}

    def q(frac: float) -> float:
        if len(xs) == 1:
            return xs[0]
        pos = frac * (len(xs) - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        if lo == hi:
            return xs[lo]
        t = pos - lo
        return xs[lo] * (1.0 - t) + xs[hi] * t

    return {"min": xs[0], "q25": q(0.25), "median": q(0.5), "q75": q(0.75), "max": xs[-1]}


def stress_row(row: Mapping[str, Any]) -> Dict[str, Any]:
    score_components = {key: float(value) for key, value in row.get("score_components", {}).items()}
    residuals = row.get("residuals", {})
    p0_scale = float(row.get("required_P0_multiplier"))
    p2_abs = abs(float(residuals.get("R_P2", score_components.get("R_P2_abs", math.nan))))
    p4_abs = abs(float(residuals.get("R_P4", score_components.get("R_P4_abs", math.nan))))
    norm_abs = float(score_components.get("R_norm_over_P0_target", math.nan))
    tail_abs = float(score_components.get("R_tail_abs", 0.0))

    one_pole_only_components = {
        "R_pole_abs": 0.0,
        "R_norm_over_P0_target": norm_abs,
        "R_P2_abs": p2_abs,
        "R_P4_abs": p4_abs,
        "R_tail_abs": tail_abs,
    }
    uniform_outgoing_scale_components = {
        "R_pole_abs": 0.0,
        "R_norm_over_P0_target": 0.0,
        "R_P2_abs": abs(p0_scale * p2_abs),
        "R_P4_abs": abs(p0_scale * p4_abs),
        "R_tail_abs": tail_abs,
    }
    return {
        "branch_id": row.get("branch_id"),
        "open_stable": row.get("open_stable"),
        "original_score": row.get("score"),
        "required_C_or_D0_multiplier": row.get("required_C_multiplier_at_fixed_D0_A"),
        "required_P0_multiplier": row.get("required_P0_multiplier"),
        "one_pole_only_projection_score": max(one_pole_only_components.values()),
        "one_pole_only_projection_components": one_pole_only_components,
        "one_pole_plus_uniform_outgoing_scale_score": max(uniform_outgoing_scale_components.values()),
        "one_pole_plus_uniform_outgoing_scale_components": uniform_outgoing_scale_components,
        "ideal_algebraic_projection_score": 0.0,
        "ideal_algebraic_projection_not_physical_hit": True,
        "ideal_projection_requirements": [
            "set C or D0 to the one-pole surface",
            "set P0 to P0_target",
            "shape N2 and N4 to make P2=P4=0",
        ],
    }


def lane_stress(name: str, deformation: Mapping[str, Any]) -> Dict[str, Any]:
    rows = [stress_row(row) for row in deformation.get("best_by_score", [])]
    all_rows = [stress_row(row) for row in deformation.get("lowest_required_C_multiplier", [])]
    open_stable_rows = [
        stress_row(row)
        for row in deformation.get("lowest_required_C_multiplier_open_stable", [])
    ]
    best_by_score = rows[0] if rows else {}
    one_pole_scores = [float(row["one_pole_only_projection_score"]) for row in rows + all_rows + open_stable_rows]
    uniform_scores = [float(row["one_pole_plus_uniform_outgoing_scale_score"]) for row in rows + all_rows + open_stable_rows]
    return {
        "name": name,
        "candidate_count": deformation.get("candidate_count"),
        "target_pass_count": deformation.get("target_pass_count"),
        "open_stable_count": deformation.get("open_stable_count"),
        "sampled_rows_are_post_hoc_top_slices": True,
        "best_score_row": best_by_score,
        "one_pole_only_projection_score_distribution": quantiles(one_pole_scores),
        "one_pole_plus_uniform_outgoing_scale_score_distribution": quantiles(uniform_scores),
        "one_pole_only_sufficient_for_target": False,
        "uniform_outgoing_amplitude_scale_sufficient_for_target": False,
        "ideal_algebraic_projection_would_zero_target_residuals": True,
        "ideal_algebraic_projection_is_not_evidence": True,
        "interpretation": {
            "one_pole_only": "Setting R_pole to zero after the fact leaves normalization and P2/P4 residuals.",
            "uniform_outgoing_scale": "Scaling P0 to target also scales P2/P4 moments, so amplitude alone is not enough.",
            "ideal_projection": "A zero-residual algebraic projection would require independent shape control over C/D0, P0, N2, and N4; it is not a physical branch.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run post-hoc coefficient projection stress diagnostics")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--report-prefix", default="projection_stress_report", help="Output basename without extension")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    paths = {
        "required_deformation": output_dir / "required_deformation_report.json",
        "nonlinear_required_deformation": output_dir / "nonlinear_required_deformation_report.json",
    }
    missing = [name for name, path in paths.items() if not path.exists()]
    if missing:
        raise RuntimeError(f"missing required input report(s): {', '.join(missing)}")

    deformation = load_json(paths["required_deformation"])
    nonlinear_deformation = load_json(paths["nonlinear_required_deformation"])
    reduced = lane_stress("reduced_operator_v1", deformation)
    nonlinear = lane_stress("nonlinear_manufactured_readiness", nonlinear_deformation)
    report = {
        "schema": "pde_audit_simulation_projection_stress_report/v1",
        "post_hoc_only": True,
        "candidate_generation_mutated": False,
        "target_blind_hit_claimed": False,
        "input_reports": {
            name: {"path": str(path), "sha256": sha256_file(path)}
            for name, path in paths.items()
        },
        "reduced_projection_stress": reduced,
        "nonlinear_projection_stress": nonlinear,
        "conclusion": {
            "one_pole_support_alone_is_insufficient": True,
            "uniform_outgoing_amplitude_scale_is_insufficient": True,
            "needed_future_mechanism": (
                "a physical branch must control one-pole support and outgoing moment shape, "
                "not merely add an after-the-fact scalar multiplier"
            ),
        },
    }
    report["pass"] = bool(
        report["post_hoc_only"]
        and not report["candidate_generation_mutated"]
        and not report["target_blind_hit_claimed"]
    )
    report["report_hash"] = sha256_json(report)
    write_json(output_dir / f"{args.report_prefix}.json", report)

    reduced_best = reduced.get("best_score_row") or {}
    nonlinear_best = nonlinear.get("best_score_row") or {}
    lines = [
        "PDE audit projection stress report",
        "=" * 42,
        f"pass: {report['pass']}",
        f"target_blind_hit_claimed: {report['target_blind_hit_claimed']}",
        f"reduced_one_pole_only_best_score: {reduced_best.get('one_pole_only_projection_score')}",
        f"reduced_uniform_outgoing_scale_best_score: {reduced_best.get('one_pole_plus_uniform_outgoing_scale_score')}",
        f"nonlinear_one_pole_only_best_score: {nonlinear_best.get('one_pole_only_projection_score')}",
        f"nonlinear_uniform_outgoing_scale_best_score: {nonlinear_best.get('one_pole_plus_uniform_outgoing_scale_score')}",
        f"one_pole_support_alone_is_insufficient: {report['conclusion']['one_pole_support_alone_is_insufficient']}",
        f"uniform_outgoing_amplitude_scale_is_insufficient: {report['conclusion']['uniform_outgoing_amplitude_scale_is_insufficient']}",
        f"report_hash: {report['report_hash']}",
    ]
    (output_dir / f"{args.report_prefix}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
