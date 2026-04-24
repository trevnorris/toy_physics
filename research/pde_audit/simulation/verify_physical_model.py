#!/usr/bin/env python3
"""Verify the current physical nonlinear export boundary.

This check is expected to pass by confirming that a strict physical nonlinear
moving-throat exporter is not yet permitted by the frozen ledger.  It prevents
the manufactured nonlinear readiness lane from being silently reclassified as a
physical branch realization.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Set

import physical_nonlinear_model as physical_model


SIM_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SIM_DIR.parents[2]
OUTPUT_DIR = SIM_DIR / "output"
BANNED_IMPORT_PREFIXES = ("stage_v2_",)
ALLOWED_IMPORTS = {
    "physical_nonlinear_model.py": {"__future__", "hashlib", "json", "typing"},
    "verify_physical_model.py": {
        "__future__",
        "argparse",
        "ast",
        "hashlib",
        "json",
        "pathlib",
        "physical_nonlinear_model",
        "typing",
    },
}
EXPECTED_BLOCKERS = {
    "strict_parent_dynamic_wall_not_promoted",
    "full_stationary_branch_equations_not_frozen",
    "bdg_support_spectrum_branch_conditional",
    "maxwell_mixed_outgoing_not_unique_parent_theorem",
}
EXPECTED_EQUATIONS = {
    "promoted_confinement_variation",
    "effective_quadratic_wall_operator",
    "axisymmetric_collective_reduction",
    "stable_bdg_schur_kernel",
    "maxwell_mixed_self_energy",
}


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
        elif isinstance(node, ast.ImportFrom) and node.module:
            roots.add(node.module.split(".", 1)[0] if node.level == 0 else "__future__")
    return roots


def source_import_check() -> Dict[str, Any]:
    issues: List[str] = []
    details: Dict[str, Any] = {}
    for filename, allowed in ALLOWED_IMPORTS.items():
        imports = imported_roots(SIM_DIR / filename)
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
    return {"name": "source_import_boundary", "pass": not issues, "issues": issues, "details": details}


def status_hash_check(status: Mapping[str, Any]) -> Dict[str, Any]:
    status_without_hash = dict(status)
    status_without_hash.pop("status_hash", None)
    recomputed = physical_model.stable_hash(status_without_hash)
    stored = status.get("status_hash")
    return {
        "name": "status_hash",
        "pass": stored == recomputed,
        "stored_status_hash": stored,
        "recomputed_status_hash": recomputed,
    }


def export_boundary_check(status: Mapping[str, Any]) -> Dict[str, Any]:
    expected_values = {
        "strict_parent_dynamic_wall_pass": False,
        "effective_wall_closure_available": True,
        "physical_export_permitted": False,
        "packets_emitted": False,
        "target_residuals_computed": False,
    }
    mismatches = {
        key: {"expected": expected, "actual": status.get(key)}
        for key, expected in expected_values.items()
        if status.get(key) is not expected
    }
    return {
        "name": "physical_export_boundary",
        "pass": not mismatches,
        "expected_values": expected_values,
        "mismatches": mismatches,
    }


def inventory_check(status: Mapping[str, Any]) -> Dict[str, Any]:
    blockers = status.get("blocking_reasons") or []
    equations = status.get("equation_inventory") or []
    blocker_ids = {str(item.get("id")) for item in blockers}
    equation_names = {str(item.get("name")) for item in equations}
    missing_blockers = sorted(EXPECTED_BLOCKERS - blocker_ids)
    missing_equations = sorted(EXPECTED_EQUATIONS - equation_names)
    required_items = status.get("required_before_physical_export") or []
    issues: List[str] = []
    evidence_details: List[Dict[str, Any]] = []
    source_cache: Dict[str, Dict[str, str]] = {}
    if status.get("schema") != "pde_audit_physical_nonlinear_model_status/v1":
        issues.append("unexpected physical model status schema")
    for item in blockers + equations:
        source = item.get("source")
        item_name = item.get("id") or item.get("name")
        if not source:
            issues.append(f"{item_name} has no source citation")
            continue
        source_path = PROJECT_ROOT / str(source)
        if not source_path.exists():
            issues.append(f"{item_name} cites missing source {source}")
            continue
        if str(source_path) not in source_cache:
            text = source_path.read_text(encoding="utf-8")
            source_cache[str(source_path)] = {
                "text": text,
                "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
            }
        evidence_patterns = [str(pattern) for pattern in item.get("evidence_patterns", [])]
        missing_patterns = [
            pattern for pattern in evidence_patterns
            if pattern not in source_cache[str(source_path)]["text"]
        ]
        if not evidence_patterns:
            issues.append(f"{item_name} has no evidence patterns")
        for pattern in missing_patterns:
            issues.append(f"{item_name} evidence pattern not found in {source}: {pattern!r}")
        evidence_details.append({
            "item": item_name,
            "source": str(source),
            "source_sha256": source_cache[str(source_path)]["sha256"],
            "evidence_pattern_count": len(evidence_patterns),
            "missing_patterns": missing_patterns,
        })
    if missing_blockers:
        issues.append(f"missing blockers: {', '.join(missing_blockers)}")
    if missing_equations:
        issues.append(f"missing equations: {', '.join(missing_equations)}")
    if len(required_items) < 4:
        issues.append("physical export requirements list is too short")
    return {
        "name": "ledger_inventory_and_blockers",
        "pass": not issues,
        "issues": issues,
        "blocker_ids": sorted(blocker_ids),
        "equation_names": sorted(equation_names),
        "required_before_physical_export_count": len(required_items),
        "evidence_details": evidence_details,
    }


def no_stale_physical_packets_check(output_dir: Path) -> Dict[str, Any]:
    packet_dir = output_dir / "physical_nonlinear_packets"
    manifest_path = output_dir / "physical_nonlinear_manifest.json"
    packet_paths = sorted(packet_dir.glob("*.json")) if packet_dir.exists() else []
    issues: List[str] = []
    if packet_paths:
        issues.append(f"physical packet directory contains {len(packet_paths)} JSON packet(s)")
    if manifest_path.exists():
        issues.append("physical nonlinear manifest exists even though export is not permitted")
    return {
        "name": "no_physical_packets_emitted",
        "pass": not issues,
        "issues": issues,
        "packet_dir": str(packet_dir),
        "packet_count": len(packet_paths),
        "manifest_path": str(manifest_path),
        "manifest_exists": manifest_path.exists(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify that strict physical nonlinear export remains blocked")
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR), help="Simulation output directory")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    status = physical_model.build_status()
    checks = [
        source_import_check(),
        status_hash_check(status),
        export_boundary_check(status),
        inventory_check(status),
        no_stale_physical_packets_check(output_dir),
    ]
    report = {
        "schema": "pde_audit_physical_model_guard/v1",
        "physical_model_status": status,
        "strict_parent_dynamic_wall_pass": status["strict_parent_dynamic_wall_pass"],
        "effective_wall_closure_available": status["effective_wall_closure_available"],
        "physical_export_permitted": status["physical_export_permitted"],
        "packets_emitted": status["packets_emitted"],
        "target_residuals_computed": status["target_residuals_computed"],
        "blocker_count": len(status["blocking_reasons"]),
        "required_before_physical_export_count": len(status["required_before_physical_export"]),
        "checks": checks,
        "pass": all(bool(check["pass"]) for check in checks),
        "failed_checks": [check["name"] for check in checks if not check["pass"]],
    }
    report["report_hash"] = physical_model.stable_hash(report)
    write_json(output_dir / "physical_model_status.json", report)

    lines = [
        "PDE audit physical nonlinear model guard",
        "=" * 43,
        f"pass: {report['pass']}",
        f"checks_passed: {sum(1 for check in checks if check['pass'])}/{len(checks)}",
        f"strict_parent_dynamic_wall_pass: {report['strict_parent_dynamic_wall_pass']}",
        f"effective_wall_closure_available: {report['effective_wall_closure_available']}",
        f"physical_export_permitted: {report['physical_export_permitted']}",
        f"packets_emitted: {report['packets_emitted']}",
        f"target_residuals_computed: {report['target_residuals_computed']}",
        f"blocker_count: {report['blocker_count']}",
        f"status_hash: {status['status_hash']}",
        f"report_hash: {report['report_hash']}",
    ]
    for check in checks:
        lines.append(f"{'PASS' if check['pass'] else 'FAIL'}  {check['name']}")
    (output_dir / "physical_model_status.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
