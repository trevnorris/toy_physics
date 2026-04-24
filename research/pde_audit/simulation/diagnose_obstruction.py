#!/usr/bin/env python3
"""Post-hoc obstruction diagnostics for frozen simulation candidates.

This script reads already-evaluated candidate reports.  It does not generate or
modify candidates.
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
    if not values:
        return {}
    xs = sorted(values)

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


def parse_name(branch_id: str) -> Dict[str, str]:
    parts = branch_id.split("_", 5)
    tail = parts[-1] if len(parts) >= 6 else branch_id
    chunks = tail.split("__")
    return {
        "geometry": chunks[0] if len(chunks) > 0 else "",
        "couplings": chunks[1] if len(chunks) > 1 else "",
        "operator": chunks[2] if len(chunks) > 2 else "",
    }


def candidate_metrics(report: Mapping[str, Any]) -> Dict[str, Any]:
    obs = report["observable_packet"]
    grouped = obs["grouped"]
    result = report["result"]
    D0 = float(grouped["D0"]["bar"])
    A = float(grouped["M"]["bar"] + grouped["B2"]["bar"] + grouped["Z2"]["bar"])
    C = float(grouped["B4"]["bar"] + grouped["Z4"]["bar"])
    required_C = 3.0 * A * A / D0 if D0 != 0.0 else math.nan
    ratio = D0 * C / (3.0 * A * A) if A != 0.0 else math.nan
    C_shortfall = required_C - C if math.isfinite(required_C) else math.nan
    parsed = parse_name(str(result["branch_id"]))
    return {
        "branch_id": result["branch_id"],
        "geometry": parsed["geometry"],
        "couplings": parsed["couplings"],
        "operator": parsed["operator"],
        "score": float(result["score"]),
        "target_packet_pass": bool(result["target_packet_pass"]),
        "open_gate_pass": bool(result["open_gate_pass"]),
        "stability_gate_pass": bool(result["stability_gate_pass"]),
        "dominant_score_component": result["dominant_score_component"],
        "D0": D0,
        "A_MplusB2plusZ2": A,
        "C_B4plusZ4": C,
        "required_C_for_one_pole": required_C,
        "C_shortfall_for_one_pole": C_shortfall,
        "one_pole_ratio_D0C_over_3A2": ratio,
        "R_pole": float(obs["residuals"]["R_pole"]),
        "R_norm": float(obs["residuals"]["R_norm"]),
        "P0": float(obs["grouped"]["P0"]["bar"]),
        "P0_target": float(obs["target_values"]["P0_target"]),
        "P0_over_target": float(obs["grouped"]["P0"]["bar"]) / float(obs["target_values"]["P0_target"]),
        "R_P2": float(obs["residuals"]["R_P2"]),
        "R_P4": float(obs["residuals"]["R_P4"]),
    }


def grouped_minima(metrics: List[Mapping[str, Any]], key: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    groups = sorted({str(item[key]) for item in metrics})
    for group in groups:
        subset = [item for item in metrics if str(item[key]) == group]
        best = min(subset, key=lambda item: float(item["score"]))
        out[group] = {
            "count": len(subset),
            "best_branch_id": best["branch_id"],
            "best_score": best["score"],
            "best_one_pole_ratio": best["one_pole_ratio_D0C_over_3A2"],
            "best_C_shortfall": best["C_shortfall_for_one_pole"],
            "target_pass_count": sum(1 for item in subset if item["target_packet_pass"]),
            "open_stable_count": sum(1 for item in subset if item["open_gate_pass"] and item["stability_gate_pass"]),
        }
    return out


def distribution_block(metrics: List[Mapping[str, Any]]) -> Dict[str, Any]:
    ratios = [float(item["one_pole_ratio_D0C_over_3A2"]) for item in metrics if math.isfinite(float(item["one_pole_ratio_D0C_over_3A2"]))]
    shortfalls = [float(item["C_shortfall_for_one_pole"]) for item in metrics if math.isfinite(float(item["C_shortfall_for_one_pole"]))]
    p0_ratios = [float(item["P0_over_target"]) for item in metrics if math.isfinite(float(item["P0_over_target"]))]
    near_one_pole = [item for item in metrics if math.isfinite(float(item["one_pole_ratio_D0C_over_3A2"])) and abs(1.0 - float(item["one_pole_ratio_D0C_over_3A2"])) < 0.1]
    return {
        "candidate_count": len(metrics),
        "target_pass_count": sum(1 for item in metrics if item["target_packet_pass"]),
        "one_pole_ratio_distribution": quantiles(ratios),
        "C_shortfall_distribution": quantiles(shortfalls),
        "P0_over_target_distribution": quantiles(p0_ratios),
        "near_one_pole_count_abs_ratio_minus_1_lt_0p1": len(near_one_pole),
        "all_one_pole_ratios_below_one": all(r < 1.0 for r in ratios) if ratios else None,
        "all_C_shortfalls_positive": all(s > 0.0 for s in shortfalls) if shortfalls else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Diagnose one-pole obstruction in frozen simulation candidates")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--candidate-dir", default="candidate_reports", help="Candidate report directory under output-dir")
    parser.add_argument("--report-prefix", default="obstruction_report", help="Output report basename without extension")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    detail_dir = output_dir / args.candidate_dir
    detail_paths = sorted(detail_dir.glob("*.json"))
    metrics = [candidate_metrics(load_json(path)) for path in detail_paths]
    if not metrics:
        raise RuntimeError(f"no candidate reports found under {detail_dir}")

    open_stable_metrics = [item for item in metrics if item["open_gate_pass"] and item["stability_gate_pass"]]
    best_by_score = sorted(metrics, key=lambda item: float(item["score"]))[:10]
    best_by_ratio = sorted(metrics, key=lambda item: abs(1.0 - float(item["one_pole_ratio_D0C_over_3A2"])))[:10]

    all_candidates = distribution_block(metrics)
    open_stable = distribution_block(open_stable_metrics)
    obstruction_report = {
        "schema": "pde_audit_simulation_obstruction_report/v1",
        "candidate_count": len(metrics),
        "target_pass_count": sum(1 for item in metrics if item["target_packet_pass"]),
        "open_stable_count": sum(1 for item in metrics if item["open_gate_pass"] and item["stability_gate_pass"]),
        "all_candidates": all_candidates,
        "open_stable_candidates": open_stable,
        "one_pole_ratio_distribution": all_candidates["one_pole_ratio_distribution"],
        "C_shortfall_distribution": all_candidates["C_shortfall_distribution"],
        "P0_over_target_distribution": all_candidates["P0_over_target_distribution"],
        "near_one_pole_count_abs_ratio_minus_1_lt_0p1": all_candidates["near_one_pole_count_abs_ratio_minus_1_lt_0p1"],
        "all_one_pole_ratios_below_one": all_candidates["all_one_pole_ratios_below_one"],
        "all_C_shortfalls_positive": all_candidates["all_C_shortfalls_positive"],
        "best_by_score": best_by_score,
        "best_by_one_pole_ratio": best_by_ratio,
        "by_operator": grouped_minima(metrics, "operator"),
        "by_geometry": grouped_minima(metrics, "geometry"),
        "by_couplings": grouped_minima(metrics, "couplings"),
        "interpretation": {
            "one_pole_surface": "R_pole = D0*C - 3*A^2, where A=M+B2+Z2 and C=B4+Z4.",
            "ratio_meaning": "D0*C/(3*A^2)=1 is required for the one-pole condition.",
            "C_shortfall_meaning": "Positive C_shortfall means the branch needs a larger quartic reservoir C at fixed D0 and A.",
            "post_hoc_only": "This report classifies frozen candidates and does not generate or mutate candidates.",
        },
    }
    obstruction_report["report_hash"] = sha256_json(obstruction_report)
    write_json(output_dir / f"{args.report_prefix}.json", obstruction_report)

    lines = [
        "PDE audit simulation obstruction report",
        "=" * 52,
        f"candidate_count: {len(metrics)}",
        f"target_pass_count: {obstruction_report['target_pass_count']}",
        f"open_stable_count: {obstruction_report['open_stable_count']}",
        f"all_candidates_one_pole_ratios_below_one: {all_candidates['all_one_pole_ratios_below_one']}",
        f"open_stable_one_pole_ratios_below_one: {open_stable['all_one_pole_ratios_below_one']}",
        f"open_stable_C_shortfalls_positive: {open_stable['all_C_shortfalls_positive']}",
        f"open_stable_near_one_pole_count_abs_ratio_minus_1_lt_0p1: {open_stable['near_one_pole_count_abs_ratio_minus_1_lt_0p1']}",
    ]
    ratio_dist = open_stable["one_pole_ratio_distribution"]
    if ratio_dist:
        lines.extend([
            f"open_stable_one_pole_ratio_min: {ratio_dist['min']:.12g}",
            f"open_stable_one_pole_ratio_median: {ratio_dist['median']:.12g}",
            f"open_stable_one_pole_ratio_max: {ratio_dist['max']:.12g}",
        ])
    best = best_by_score[0]
    lines.extend([
        f"best_score_branch_id: {best['branch_id']}",
        f"best_score: {best['score']:.12g}",
        f"best_one_pole_ratio: {best['one_pole_ratio_D0C_over_3A2']:.12g}",
        f"best_C_shortfall: {best['C_shortfall_for_one_pole']:.12g}",
        f"report_hash: {obstruction_report['report_hash']}",
    ])
    (output_dir / f"{args.report_prefix}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
