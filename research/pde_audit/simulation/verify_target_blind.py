#!/usr/bin/env python3
"""Executable target-blindness guard for frozen simulation packet generators."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Set


SIM_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SIM_DIR / "output"

BANNED_IMPORT_PREFIXES = ("stage_v2_",)
BANNED_PACKET_KEYS = {
    "dominant_score_component",
    "pass_flags",
    "R_norm",
    "R_P2",
    "R_P4",
    "R_pole",
    "R_tail",
    "residuals",
    "score",
    "score_components",
    "target_packet_pass",
    "target_values",
}
ALLOWED_IMPORTS_BY_MODE = {
    "reduced": {
        "generate_reduced_sweep.py": {
            "__future__",
            "argparse",
            "hashlib",
            "json",
            "math",
            "numpy",
            "pathlib",
            "reduced_fem",
            "typing",
        },
        "reduced_fem.py": {
            "__future__",
            "math",
            "numpy",
            "scipy",
            "typing",
        },
    },
    "nonlinear": {
        "nonlinear_protocol.py": {"__future__", "hashlib", "json", "typing"},
        "verify_nonlinear_solver.py": {
            "__future__",
            "argparse",
            "ast",
            "hashlib",
            "json",
            "math",
            "nonlinear_protocol",
            "numpy",
            "pathlib",
            "typing",
        },
        "generate_nonlinear_packets.py": {
            "__future__",
            "argparse",
            "hashlib",
            "json",
            "math",
            "nonlinear_protocol",
            "numpy",
            "pathlib",
            "typing",
            "verify_nonlinear_solver",
        },
    },
}


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def imported_roots(path: Path) -> Set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    roots: Set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                roots.add(alias.name.split(".", 1)[0])
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                roots.add(node.module.split(".", 1)[0] if node.level == 0 else "__future__")
    return roots


def find_forbidden_packet_keys(obj: Any, prefix: str = "") -> List[str]:
    paths: List[str] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            key_text = str(key)
            path = f"{prefix}.{key_text}" if prefix else key_text
            if key_text in BANNED_PACKET_KEYS:
                paths.append(path)
            paths.extend(find_forbidden_packet_keys(value, path))
    elif isinstance(obj, list):
        for idx, value in enumerate(obj):
            paths.extend(find_forbidden_packet_keys(value, f"{prefix}[{idx}]"))
    return paths


def verify_source_imports(mode: str) -> Dict[str, Any]:
    details: Dict[str, Any] = {}
    issues: List[str] = []
    for filename, allowed in ALLOWED_IMPORTS_BY_MODE[mode].items():
        path = SIM_DIR / filename
        imports = imported_roots(path)
        banned = sorted(root for root in imports if any(root.startswith(prefix) for prefix in BANNED_IMPORT_PREFIXES))
        unexpected = sorted(imports - allowed)
        details[filename] = {
            "imports": sorted(imports),
            "allowed_imports": sorted(allowed),
            "banned_imports": banned,
            "unexpected_imports": unexpected,
        }
        for root in banned:
            issues.append(f"{filename} imports banned target-evaluation module root {root!r}")
        for root in unexpected:
            issues.append(f"{filename} imports undeclared module root {root!r}")
    return {"pass": not issues, "issues": issues, "details": details}


def verify_manifest_and_packets(output_dir: Path, manifest_name: str) -> Dict[str, Any]:
    manifest_path = output_dir / manifest_name
    if not manifest_path.exists():
        return {"pass": False, "issues": [f"missing manifest {manifest_path}"], "details": {}}

    manifest = load_json(manifest_path)
    issues: List[str] = []
    stored_manifest_hash = manifest.get("manifest_hash")
    manifest_without_hash = dict(manifest)
    manifest_without_hash.pop("manifest_hash", None)
    computed_manifest_hash = sha256_json(manifest_without_hash)
    if stored_manifest_hash != computed_manifest_hash:
        issues.append("manifest_hash does not match manifest content")

    assertion = manifest.get("generator_target_blind_assertion", {})
    expected_false = {
        "target_evaluation_modules_imported",
        "target_residuals_computed",
        "packets_contain_target_outputs",
    }
    for key in sorted(expected_false):
        if assertion.get(key) is not False:
            issues.append(f"manifest generator_target_blind_assertion.{key} is not false")

    packets = manifest.get("packets", [])
    if manifest.get("candidate_count") != len(packets):
        issues.append("manifest candidate_count does not match packets length")

    branch_ids: Set[str] = set()
    packet_issue_count = 0
    packet_count = 0
    for entry in packets:
        packet_count += 1
        rel_path = entry.get("path")
        packet_path = output_dir / str(rel_path)
        if not packet_path.exists():
            issues.append(f"missing packet {rel_path}")
            packet_issue_count += 1
            continue
        packet = load_json(packet_path)
        branch_id = str(packet.get("branch_id"))
        if branch_id in branch_ids:
            issues.append(f"duplicate branch_id {branch_id!r}")
            packet_issue_count += 1
        branch_ids.add(branch_id)
        packet_hash = sha256_json(packet)
        if entry.get("packet_hash") != packet_hash:
            issues.append(f"packet_hash mismatch for {rel_path}")
            packet_issue_count += 1
        freeze = packet.get("freeze", {})
        for key in ("pre_target_freeze", "target_blind", "no_post_residual_refit"):
            if freeze.get(key) is not True:
                issues.append(f"{rel_path} freeze.{key} is not true")
                packet_issue_count += 1
        if entry.get("freeze_hash") != freeze.get("candidate_freeze_hash"):
            issues.append(f"freeze_hash mismatch for {rel_path}")
            packet_issue_count += 1
        forbidden_paths = find_forbidden_packet_keys(packet)
        if forbidden_paths:
            issues.append(f"{rel_path} contains target-output keys: {', '.join(forbidden_paths[:12])}")
            packet_issue_count += 1

    return {
        "pass": not issues,
        "issues": issues,
        "details": {
            "manifest_path": str(manifest_path),
            "stored_manifest_hash": stored_manifest_hash,
            "computed_manifest_hash": computed_manifest_hash,
            "packet_count": packet_count,
            "unique_branch_id_count": len(branch_ids),
            "packet_issue_count": packet_issue_count,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify target-blind simulation generation boundary")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    parser.add_argument("--mode", choices=sorted(ALLOWED_IMPORTS_BY_MODE), default="reduced", help="Generator source boundary to check")
    parser.add_argument("--manifest-name", default="manifest.json", help="Manifest filename under output-dir")
    parser.add_argument("--report-prefix", default="target_blind_report", help="Output report basename without extension")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    source_report = verify_source_imports(args.mode)
    packet_report = verify_manifest_and_packets(output_dir, args.manifest_name)
    report = {
        "schema": "pde_audit_simulation_target_blind_guard/v1",
        "mode": args.mode,
        "manifest_name": args.manifest_name,
        "source_imports": source_report,
        "manifest_and_packets": packet_report,
    }
    report["pass"] = bool(source_report["pass"] and packet_report["pass"])
    report["issue_count"] = len(source_report["issues"]) + len(packet_report["issues"])
    report["report_hash"] = sha256_json(report)
    write_json(output_dir / f"{args.report_prefix}.json", report)

    lines = [
        "PDE audit simulation target-blind guard",
        "=" * 48,
        f"mode: {args.mode}",
        f"manifest_name: {args.manifest_name}",
        f"pass: {report['pass']}",
        f"issue_count: {report['issue_count']}",
        f"source_imports_pass: {source_report['pass']}",
        f"manifest_and_packets_pass: {packet_report['pass']}",
        f"packet_count: {packet_report.get('details', {}).get('packet_count')}",
        f"report_hash: {report['report_hash']}",
    ]
    for issue in source_report["issues"] + packet_report["issues"]:
        lines.append(f"ISSUE: {issue}")
    (output_dir / f"{args.report_prefix}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
