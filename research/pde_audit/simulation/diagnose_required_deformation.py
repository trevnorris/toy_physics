#!/usr/bin/env python3
"""Post-hoc required-deformation map for frozen simulation candidates.

This script reads already-evaluated candidate reports.  It does not generate,
modify, rank for generation, or refit candidates.  Its purpose is to quantify
which physical coefficient changes would be required for the frozen candidates
to approach the target surface.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional


OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


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

    return {
        "min": xs[0],
        "q25": q(0.25),
        "median": q(0.5),
        "q75": q(0.75),
        "max": xs[-1],
    }


def candidate_name_from_branch_id(branch_id: str) -> str:
    for prefix in ("sim_reduced_open_throat_", "sim_nonlinear_open_throat_"):
        if branch_id.startswith(prefix):
            remainder = branch_id[len(prefix):]
            parts = remainder.split("_", 1)
            return parts[1] if len(parts) == 2 else remainder
    return branch_id


def parse_candidate_name(branch_id: str) -> Dict[str, str]:
    name = candidate_name_from_branch_id(branch_id)
    chunks = name.split("__")
    return {
        "candidate_name": name,
        "geometry": chunks[0] if len(chunks) > 0 else name,
        "couplings": chunks[1] if len(chunks) > 1 else "not_encoded",
        "operator": chunks[2] if len(chunks) > 2 else "not_encoded",
    }


def safe_div(num: float, den: float) -> float:
    if den == 0.0:
        return math.nan
    return num / den


def candidate_deformation(report: Mapping[str, Any]) -> Dict[str, Any]:
    obs = report["observable_packet"]
    grouped = obs["grouped"]
    result = report["result"]
    branch_id = str(result["branch_id"])
    parsed = parse_candidate_name(branch_id)

    D0 = float(grouped["D0"]["bar"])
    A = float(grouped["M"]["bar"] + grouped["B2"]["bar"] + grouped["Z2"]["bar"])
    C = float(grouped["B4"]["bar"] + grouped["Z4"]["bar"])
    one_pole_ratio = safe_div(D0 * C, 3.0 * A * A)
    required_C = safe_div(3.0 * A * A, D0)
    required_D0 = safe_div(3.0 * A * A, C)
    required_A = math.sqrt(max(D0 * C / 3.0, 0.0)) if D0 * C >= 0.0 else math.nan
    P0 = float(grouped["P0"]["bar"])
    P0_target = float(obs["target_values"]["P0_target"])
    score_components = {key: float(value) for key, value in result.get("score_components", {}).items()}
    max_component = max(score_components.values()) if score_components else math.nan

    return {
        "branch_id": branch_id,
        "candidate_name": parsed["candidate_name"],
        "geometry": parsed["geometry"],
        "couplings": parsed["couplings"],
        "operator": parsed["operator"],
        "score": float(result["score"]),
        "target_packet_pass": bool(result["target_packet_pass"]),
        "open_gate_pass": bool(result["open_gate_pass"]),
        "stability_gate_pass": bool(result["stability_gate_pass"]),
        "open_stable": bool(result["open_gate_pass"] and result["stability_gate_pass"]),
        "dominant_score_component": result.get("dominant_score_component"),
        "score_components": score_components,
        "pass_flags": result.get("pass_flags", {}),
        "D0": D0,
        "A_MplusB2plusZ2": A,
        "C_B4plusZ4": C,
        "one_pole_ratio_D0C_over_3A2": one_pole_ratio,
        "required_C_for_one_pole": required_C,
        "required_D0_for_one_pole": required_D0,
        "required_A_for_one_pole": required_A,
        "required_C_multiplier_at_fixed_D0_A": safe_div(required_C, C),
        "required_D0_multiplier_at_fixed_A_C": safe_div(required_D0, D0),
        "required_A_multiplier_at_fixed_D0_C": safe_div(required_A, A),
        "required_C_additive_shortfall": required_C - C if math.isfinite(required_C) else math.nan,
        "required_D0_additive_shortfall": required_D0 - D0 if math.isfinite(required_D0) else math.nan,
        "required_A_additive_change": required_A - A if math.isfinite(required_A) else math.nan,
        "P0": P0,
        "P0_target": P0_target,
        "P0_over_target": safe_div(P0, P0_target),
        "required_P0_multiplier": safe_div(P0_target, P0),
        "residuals": {
            "R_pole": float(obs["residuals"]["R_pole"]),
            "R_norm": float(obs["residuals"]["R_norm"]),
            "R_P2": float(obs["residuals"]["R_P2"]),
            "R_P4": float(obs["residuals"]["R_P4"]),
            "R_tail": float(obs["residuals"]["R_tail"]),
        },
        "max_score_component": max_component,
    }


def distribution_block(rows: List[Mapping[str, Any]]) -> Dict[str, Any]:
    one_pole_ratios = [float(row["one_pole_ratio_D0C_over_3A2"]) for row in rows]
    C_multipliers = [float(row["required_C_multiplier_at_fixed_D0_A"]) for row in rows]
    D0_multipliers = [float(row["required_D0_multiplier_at_fixed_A_C"]) for row in rows]
    A_multipliers = [float(row["required_A_multiplier_at_fixed_D0_C"]) for row in rows]
    P0_multipliers = [float(row["required_P0_multiplier"]) for row in rows]
    scores = [float(row["score"]) for row in rows]
    return {
        "candidate_count": len(rows),
        "target_pass_count": sum(1 for row in rows if row["target_packet_pass"]),
        "open_stable_count": sum(1 for row in rows if row["open_stable"]),
        "score_distribution": quantiles(scores),
        "one_pole_ratio_distribution": quantiles(one_pole_ratios),
        "required_C_multiplier_distribution": quantiles(C_multipliers),
        "required_D0_multiplier_distribution": quantiles(D0_multipliers),
        "required_A_multiplier_distribution": quantiles(A_multipliers),
        "required_P0_multiplier_distribution": quantiles(P0_multipliers),
        "near_one_pole_count_abs_ratio_minus_1_lt_0p1": sum(
            1 for value in one_pole_ratios if math.isfinite(value) and abs(value - 1.0) < 0.1
        ),
    }


def grouped_summary(rows: List[Mapping[str, Any]], key: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for group in sorted({str(row[key]) for row in rows}):
        subset = [row for row in rows if str(row[key]) == group]
        best = min(subset, key=lambda row: float(row["score"]))
        best_ratio = max(subset, key=lambda row: float(row["one_pole_ratio_D0C_over_3A2"]))
        out[group] = {
            "candidate_count": len(subset),
            "target_pass_count": sum(1 for row in subset if row["target_packet_pass"]),
            "open_stable_count": sum(1 for row in subset if row["open_stable"]),
            "best_score_branch_id": best["branch_id"],
            "best_score": best["score"],
            "best_score_required_C_multiplier": best["required_C_multiplier_at_fixed_D0_A"],
            "best_score_required_P0_multiplier": best["required_P0_multiplier"],
            "highest_one_pole_ratio_branch_id": best_ratio["branch_id"],
            "highest_one_pole_ratio": best_ratio["one_pole_ratio_D0C_over_3A2"],
            "lowest_required_C_multiplier": best_ratio["required_C_multiplier_at_fixed_D0_A"],
            "median_required_C_multiplier": (distribution_block(subset).get("required_C_multiplier_distribution") or {}).get("median"),
        }
    return out


def mechanism_assessment(open_stable: Mapping[str, Any], dominant_counts: Mapping[str, int]) -> Dict[str, Any]:
    c_dist = open_stable.get("required_C_multiplier_distribution") or {}
    p0_dist = open_stable.get("required_P0_multiplier_distribution") or {}
    ratio_dist = open_stable.get("one_pole_ratio_distribution") or {}
    min_c_mult = c_dist.get("min")
    median_p0_mult = p0_dist.get("median")
    max_ratio = ratio_dist.get("max")
    return {
        "primary_obstruction": "one_pole_gap" if dominant_counts.get("R_pole_abs", 0) else "mixed_residuals",
        "dominant_score_component_counts": dict(dominant_counts),
        "one_pole_gap_interpretation": (
            "At fixed A and D0, C must be multiplied by the reported required_C_multiplier. "
            "Equivalently, at fixed A and C, D0 must be multiplied by the same factor."
        ),
        "normalization_gap_interpretation": (
            "required_P0_multiplier measures the outgoing/normalization amplitude multiplier at fixed shape."
        ),
        "minimum_C_multiplier_open_stable": min_c_mult,
        "median_P0_multiplier_open_stable": median_p0_mult,
        "maximum_one_pole_ratio_open_stable": max_ratio,
        "local_continuation_looks_promising": bool(
            isinstance(max_ratio, float)
            and max_ratio > 0.9
            and open_stable.get("near_one_pole_count_abs_ratio_minus_1_lt_0p1", 0) > 0
        ),
        "suggested_missing_mechanism": (
            "new physical support that increases C or D0 by order required_C_multiplier, "
            "plus an outgoing/mixed-sector normalization mechanism if required_P0_multiplier remains large"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Map required post-hoc deformations for frozen candidates")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--candidate-dir", default="candidate_reports", help="Candidate report directory under output-dir")
    parser.add_argument("--report-prefix", default="required_deformation_report", help="Output basename without extension")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    detail_dir = output_dir / args.candidate_dir
    detail_paths = sorted(detail_dir.glob("*.json"))
    if not detail_paths:
        raise RuntimeError(f"no candidate reports found under {detail_dir}")

    rows = [candidate_deformation(load_json(path)) for path in detail_paths]
    open_stable_rows = [row for row in rows if row["open_stable"]]
    dominant_counts: Dict[str, int] = {}
    for row in rows:
        key = str(row["dominant_score_component"])
        dominant_counts[key] = dominant_counts.get(key, 0) + 1

    best_by_score = sorted(rows, key=lambda row: float(row["score"]))[:10]
    best_by_required_C_multiplier = sorted(
        rows,
        key=lambda row: float(row["required_C_multiplier_at_fixed_D0_A"])
        if math.isfinite(float(row["required_C_multiplier_at_fixed_D0_A"]))
        else math.inf,
    )[:10]
    best_open_stable_by_required_C_multiplier = [
        row for row in best_by_required_C_multiplier if row["open_stable"]
    ][:10]

    all_block = distribution_block(rows)
    open_stable_block = distribution_block(open_stable_rows)
    report = {
        "schema": "pde_audit_simulation_required_deformation_report/v1",
        "candidate_dir": str(detail_dir),
        "post_hoc_only": True,
        "target_residuals_used_to_generate_candidates": False,
        "candidate_count": len(rows),
        "target_pass_count": sum(1 for row in rows if row["target_packet_pass"]),
        "open_stable_count": len(open_stable_rows),
        "all_candidates": all_block,
        "open_stable_candidates": open_stable_block,
        "dominant_score_component_counts": dominant_counts,
        "mechanism_assessment": mechanism_assessment(open_stable_block, dominant_counts),
        "by_geometry": grouped_summary(rows, "geometry"),
        "by_couplings": grouped_summary(rows, "couplings"),
        "by_operator": grouped_summary(rows, "operator"),
        "best_by_score": best_by_score,
        "lowest_required_C_multiplier": best_by_required_C_multiplier,
        "lowest_required_C_multiplier_open_stable": best_open_stable_by_required_C_multiplier,
        "interpretation": {
            "A_definition": "A = M + B2 + Z2",
            "C_definition": "C = B4 + Z4",
            "one_pole_surface": "D0*C/(3*A^2) = 1",
            "required_C_multiplier_at_fixed_D0_A": "Multiplier on C required to satisfy the one-pole surface with D0 and A fixed.",
            "required_D0_multiplier_at_fixed_A_C": "Multiplier on D0 required to satisfy the one-pole surface with A and C fixed.",
            "required_A_multiplier_at_fixed_D0_C": "Multiplier on A required to satisfy the one-pole surface with D0 and C fixed.",
            "required_P0_multiplier": "Multiplier on P0 required to hit P0_target with the frozen candidate shape.",
        },
    }
    report["report_hash"] = sha256_json(report)
    write_json(output_dir / f"{args.report_prefix}.json", report)

    c_dist = open_stable_block.get("required_C_multiplier_distribution") or {}
    p0_dist = open_stable_block.get("required_P0_multiplier_distribution") or {}
    ratio_dist = open_stable_block.get("one_pole_ratio_distribution") or {}
    best = best_by_score[0]
    lines = [
        "PDE audit required deformation report",
        "=" * 45,
        f"candidate_count: {len(rows)}",
        f"target_pass_count: {report['target_pass_count']}",
        f"open_stable_count: {len(open_stable_rows)}",
        f"dominant_score_component_counts: {dominant_counts}",
        f"open_stable_one_pole_ratio_max: {ratio_dist.get('max')}",
        f"open_stable_required_C_multiplier_min: {c_dist.get('min')}",
        f"open_stable_required_C_multiplier_median: {c_dist.get('median')}",
        f"open_stable_required_P0_multiplier_median: {p0_dist.get('median')}",
        f"local_continuation_looks_promising: {report['mechanism_assessment']['local_continuation_looks_promising']}",
        f"best_score_branch_id: {best['branch_id']}",
        f"best_score: {best['score']:.12g}",
        f"best_score_required_C_multiplier: {best['required_C_multiplier_at_fixed_D0_A']:.12g}",
        f"best_score_required_P0_multiplier: {best['required_P0_multiplier']:.12g}",
        f"report_hash: {report['report_hash']}",
    ]
    (output_dir / f"{args.report_prefix}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
