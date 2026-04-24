#!/usr/bin/env python3
"""Evaluate frozen simulation packets through the PDE V2 audit chain."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

AUDIT_DIR = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = AUDIT_DIR / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import stage_v2_21_branch_extraction_fixture as v21
import stage_v2_22a_profile_to_coefficient_adapter as v22a
import stage_v2_22b_solver_handoff_validator as v22b


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


def residual_score_components(packet: Mapping[str, Any]) -> Dict[str, float]:
    residuals = packet["residuals"]
    target = packet["target_values"]
    return {
        "R_pole_abs": abs(float(residuals["R_pole"])),
        "R_norm_over_P0_target": abs(float(residuals["R_norm"])) / max(abs(float(target["P0_target"])), 1e-30),
        "R_P2_abs": abs(float(residuals["R_P2"])),
        "R_P4_abs": abs(float(residuals["R_P4"])),
        "R_tail_abs": abs(float(residuals["R_tail"])),
    }


def residual_score(packet: Mapping[str, Any]) -> float:
    return max(residual_score_components(packet).values())


def quantiles(values: List[float]) -> Dict[str, float]:
    if not values:
        return {}
    xs = sorted(values)
    mid = len(xs) // 2
    median = xs[mid] if len(xs) % 2 else 0.5 * (xs[mid - 1] + xs[mid])
    return {"min": xs[0], "median": median, "max": xs[-1]}


def evaluate_packet(packet_path: Path, output_dir: Path, tol: float, details_dir_name: str) -> Dict[str, Any]:
    solver_packet = load_json(packet_path)
    packet_content_hash = sha256_json(solver_packet)
    validation = v22b.validate_solver_output(solver_packet)
    result: Dict[str, Any] = {
        "branch_id": solver_packet.get("branch_id"),
        "packet_path": str(packet_path),
        "packet_sha256": sha256_file(packet_path),
        "packet_content_hash": packet_content_hash,
        "candidate_freeze_hash": solver_packet.get("freeze", {}).get("candidate_freeze_hash"),
        "validation_pass": bool(validation.get("validation_pass")),
        "validation_error_count": int(validation.get("error_count", 0)),
        "validation_warning_count": int(validation.get("warning_count", 0)),
        "validation_issue_paths": [item.get("path") for item in validation.get("issues", [])],
    }
    if not validation["validation_pass"]:
        return result

    profile_manifest = v22b.convert_solver_output_to_v22a_profile_manifest(solver_packet)
    manifest_validation = v22a.validate_profile_manifest(profile_manifest, strict_profiles=True)
    result["profile_manifest_validation_pass"] = bool(manifest_validation["validation_pass"])
    if not manifest_validation["validation_pass"]:
        result["profile_manifest_issues"] = manifest_validation["issues"]
        return result

    v21_manifest, adapter_diag = v22a.adapt_manifest(profile_manifest)
    branch_packets = [v21.extract_branch(branch, tol=tol) for branch in v21_manifest.get("branches", [])]
    if len(branch_packets) != 1:
        result["evaluation_error"] = f"expected one branch packet, found {len(branch_packets)}"
        return result

    branch_packet = branch_packets[0]
    result.update({
        "profile_manifest_hash": sha256_json(profile_manifest),
        "v21_manifest_hash": sha256_json(v21_manifest),
        "adapter_diagnostics_hash": sha256_json(adapter_diag),
        "open_gate_pass": bool(branch_packet["pass_flags"]["open_gate_pass"]),
        "stability_gate_pass": bool(branch_packet["pass_flags"]["stability_gate_pass"]),
        "target_packet_pass": bool(branch_packet["pass_flags"]["target_packet_pass"]),
        "pass_flags": branch_packet["pass_flags"],
        "residuals": branch_packet["residuals"],
        "target_values": branch_packet["target_values"],
        "score": residual_score(branch_packet),
        "score_components": residual_score_components(branch_packet),
        "dominant_score_component": max(residual_score_components(branch_packet).items(), key=lambda kv: kv[1])[0],
        "D0_bar": branch_packet["grouped"]["D0"]["bar"],
        "P0_bar": branch_packet["grouped"]["P0"]["bar"],
    })

    detail = {
        "schema": "pde_audit_simulation_candidate_evaluation/v1",
        "solver_packet_hash": result["packet_content_hash"],
        "validation": validation,
        "profile_manifest": profile_manifest,
        "v21_manifest": v21_manifest,
        "observable_packet": branch_packet,
        "result": result,
    }
    detail["detail_hash"] = sha256_json(detail)
    details_dir = output_dir / details_dir_name
    write_json(details_dir / f"{solver_packet['branch_id']}.json", detail)
    result["detail_report"] = str((details_dir / f"{solver_packet['branch_id']}.json").relative_to(output_dir))
    result["detail_hash"] = detail["detail_hash"]
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate frozen reduced-sweep simulation packets")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--manifest-name", default="manifest.json", help="Manifest filename under output-dir")
    parser.add_argument("--details-dir", default="candidate_reports", help="Per-candidate report directory under output-dir")
    parser.add_argument("--summary-prefix", default="evaluation_summary", help="Summary output basename without extension")
    parser.add_argument("--tol", type=float, default=1e-9, help="Target extraction tolerance")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    manifest_path = output_dir / args.manifest_name
    manifest = load_json(manifest_path)
    results: List[Dict[str, Any]] = []
    details_dir = output_dir / args.details_dir
    if details_dir.exists():
        for path in details_dir.glob("*.json"):
            path.unlink()

    print("PDE audit frozen simulation evaluator")
    print("=" * 64)
    print(f"manifest: {manifest_path}")
    print(f"candidate_count: {manifest['candidate_count']}")
    print("post_hoc_ranking_only: True")

    for entry in manifest["packets"]:
        packet_path = output_dir / entry["path"]
        result = evaluate_packet(packet_path, output_dir, args.tol, args.details_dir)
        results.append(result)
        status = "TARGET_PASS" if result.get("target_packet_pass") else "target_fail"
        gates = f"open={result.get('open_gate_pass')} stable={result.get('stability_gate_pass')}"
        score = result.get("score")
        score_text = f" score={score:.3e}" if isinstance(score, float) else ""
        print(f"{status}  {result.get('branch_id')}  {gates}{score_text}")

    validation_pass_count = sum(1 for item in results if item.get("validation_pass"))
    open_stable_count = sum(1 for item in results if item.get("open_gate_pass") and item.get("stability_gate_pass"))
    target_pass_count = sum(1 for item in results if item.get("target_packet_pass"))
    best: Optional[Dict[str, Any]] = None
    scored = [item for item in results if isinstance(item.get("score"), float)]
    if scored:
        best = min(scored, key=lambda item: item["score"])
    pass_flag_counts: Dict[str, int] = {}
    for item in scored:
        for key, value in item.get("pass_flags", {}).items():
            if bool(value):
                pass_flag_counts[key] = pass_flag_counts.get(key, 0) + 1
    dominant_component_counts: Dict[str, int] = {}
    for item in scored:
        key = str(item.get("dominant_score_component"))
        dominant_component_counts[key] = dominant_component_counts.get(key, 0) + 1

    summary = {
        "schema": "pde_audit_simulation_evaluation_summary/v1",
        "manifest_path": str(manifest_path),
        "manifest_hash": manifest.get("manifest_hash"),
        "evaluation_protocol": {
            "target_residuals_used_to_generate_candidates": False,
            "post_hoc_ranking_only": True,
            "no_post_residual_refit": True,
            "tolerance": args.tol,
        },
        "candidate_count": len(results),
        "validation_pass_count": validation_pass_count,
        "open_stable_count": open_stable_count,
        "target_pass_count": target_pass_count,
        "score_distribution": quantiles([float(item["score"]) for item in scored]),
        "dominant_score_component_counts": dominant_component_counts,
        "pass_flag_counts": pass_flag_counts,
        "best_post_hoc_candidate": best,
        "results": results,
    }
    summary["summary_hash"] = sha256_json(summary)
    write_json(output_dir / f"{args.summary_prefix}.json", summary)

    text_lines = [
        "PDE audit simulation evaluation summary",
        "=" * 48,
        f"candidate_count: {len(results)}",
        f"validation_pass_count: {validation_pass_count}",
        f"open_stable_count: {open_stable_count}",
        f"target_pass_count: {target_pass_count}",
    ]
    score_dist = summary["score_distribution"]
    if score_dist:
        text_lines.extend([
            f"score_min: {score_dist['min']:.12g}",
            f"score_median: {score_dist['median']:.12g}",
            f"score_max: {score_dist['max']:.12g}",
        ])
    if best is not None:
        text_lines.extend([
            f"best_post_hoc_branch_id: {best['branch_id']}",
            f"best_post_hoc_score: {best['score']:.12g}",
            f"best_post_hoc_target_pass: {best['target_packet_pass']}",
        ])
    (output_dir / f"{args.summary_prefix}.txt").write_text("\n".join(text_lines) + "\n", encoding="utf-8")
    print("")
    print(f"target_pass_count: {target_pass_count}")
    if best is not None:
        print(f"best_post_hoc_branch_id: {best['branch_id']}")
        print(f"best_post_hoc_score: {best['score']:.3e}")
    print(f"summary_hash: {summary['summary_hash']}")

    return 0 if validation_pass_count == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
