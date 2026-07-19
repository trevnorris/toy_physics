#!/usr/bin/env python3
"""Preflight one orchestrator-owned U2 A9 external-verification leg.

This helper cannot mark a leg PASS.  It reasserts the ratified stage-0 bundle,
live environment, production code anchor, and exact every-object coverage map
immediately before control passes to a fresh external verifier.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


RATIFIED_DIGEST = "9eff1b0c49e89007aea1008cb6712b0ea495168d101ce43ddce1cffaf68749c4"
STAGE0_ANCHOR = "323b222846e2a9062330d2f25dd9cd28c57c7800"
LEGS = ("arbiter", "fidelity", "adversarial_recompute", "read_only_review")
EXPECTED_CATEGORIES = {
    "candidate_disposition": 144,
    "closure_adjudication": 112,
    "ensemble_level_1": 144,
    "ensemble_level_2": 144,
    "host_location": 144,
    "promotion": 16,
    "return_closure_ownership": 144,
    "same_route_fixture_control": 144,
    "topology_gate": 144,
}
EXPECTED_ZERO_SCOPES = {
    "EXCLUDED_decisive_incompatibility": 0,
    "positive_witness": 0,
    "selection_forcing_proof": 0,
    "closure_certificate": 0,
    "posed_BVP_template": 0,
}


class PreflightFailure(RuntimeError):
    pass


class AmendmentHalt(RuntimeError):
    pass


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=160),
        encoding="utf-8",
    )


def assert_anchored(repo: Path, anchor: str, path: Path) -> str:
    try:
        rel = path.resolve().relative_to(repo).as_posix()
    except ValueError as failure:
        raise PreflightFailure(f"task-authored A9 code outside repository: {path}") from failure
    result = subprocess.run(
        ["git", "cat-file", "blob", f"{anchor}:{rel}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise PreflightFailure(f"ASSERT_EVALUATED_CODE_CLOSURE: {rel} absent from production anchor")
    if hashlib.sha256(path.read_bytes()).digest() != hashlib.sha256(result.stdout).digest():
        raise PreflightFailure(f"ASSERT_EXECUTED_BYTES_MATCH_ANCHOR: {rel}")
    return rel


def reassert_stage0(repo: Path, source: Path, stage0: Path, digest: str) -> None:
    runtime = source / "_scratch/u2_production/a9_runtime"
    for path in (runtime, runtime / "wolfram_userbase/Licensing", runtime / "wolfram_base"):
        path.mkdir(parents=True, exist_ok=True)
    nsswitch = runtime / "nsswitch.conf"
    nsswitch.write_text("hosts: files\n", encoding="utf-8")
    verify_command = [
        "/usr/bin/python3", "-I", str(source / "u2_stage0_contract.py"),
        "--repo", str(repo), "--output-dir", str(stage0),
        "--startup-contract-commit", STAGE0_ANCHOR, "--verify",
        "--expected-digest", digest, "--recompute-environment",
    ]
    launched = [
        "/usr/bin/bwrap", "--ro-bind", "/", "/", "--dev-bind", "/dev", "/dev",
        "--proc", "/proc", "--unshare-net", "--bind", str(runtime), "/tmp",
        "--ro-bind", "/home/trevnorris/.Wolfram/Licensing", "/tmp/wolfram_userbase/Licensing",
        "--ro-bind", str(nsswitch), "/etc/nsswitch.conf", "--chdir", str(repo),
        "--setenv", "PYTHONNOUSERSITE", "1", "--setenv", "PYTHONDONTWRITEBYTECODE", "1",
        "--setenv", "TMPDIR", "/tmp", "--setenv", "WOLFRAM_USERBASE", "/tmp/wolfram_userbase",
        "--setenv", "WOLFRAM_BASE", "/tmp/wolfram_base", "--", *verify_command,
    ]
    verified = subprocess.run(
        launched, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=600,
    )
    combined = verified.stdout + verified.stderr
    if verified.returncode != 0 or "U2_STAGE0_CONTRACT_VERIFIED" not in combined:
        raise AmendmentHalt(f"stage-0/environment reassertion failed: {combined[-1600:]}")


def validate_coverage(results: dict[str, Any]) -> tuple[int, dict[str, int]]:
    coverage = results["A9_external_verification"]
    object_ids = coverage["object_ids"]
    coverage_map = coverage["coverage_map"]
    mapped_ids = [row["object_id"] for row in coverage_map]
    counted: dict[str, int] = {}
    for row in coverage_map:
        counted[row["category"]] = counted.get(row["category"], 0) + 1
        if not row["recompute_route"]:
            raise PreflightFailure(f"A9 object has no recomputation route: {row['object_id']}")
    if (
        coverage["object_count"] != len(object_ids)
        or len(object_ids) != len(set(object_ids))
        or sorted(mapped_ids) != sorted(object_ids)
        or len(mapped_ids) != len(set(mapped_ids))
        or counted != EXPECTED_CATEGORIES
        or coverage["coverage_category_counts"] != EXPECTED_CATEGORIES
        or sum(counted.values()) != len(object_ids)
        or not coverage["coverage_map_exact"]
        or not coverage["generated_category_union_exact"]
    ):
        raise PreflightFailure("A9 every-object coverage inventory is not exact")
    if (
        coverage["recomputation_granularity"] != "EVERY_OBJECT_INDIVIDUAL"
        or coverage["class_level_treatment_used"]
        or coverage["class_equivalence_proof_required"]
        or coverage["required_scope_zero_counts"] != EXPECTED_ZERO_SCOPES
        or not coverage["required_scope_zero_counts_measured"]
        or coverage["external_results_claimed"]
    ):
        raise PreflightFailure("A9 scope or external-ownership policy is not exact")
    return len(object_ids), counted


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True)
    parser.add_argument("--leg", choices=LEGS, required=True)
    parser.add_argument("--startup-contract-commit", required=True)
    parser.add_argument("--stage0-contract-digest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    repo = Path(args.repo).resolve()
    source = repo / "software/em_charge_attribute"
    stage0 = source / "reports/u2_boundary_adjudication_artifacts/stage_0_contract"
    production = source / "reports/u2_boundary_adjudication_artifacts/stage_1_production"
    anchor = args.startup_contract_commit
    if not re.fullmatch(r"[0-9a-f]{40}", anchor) or anchor == "HEAD":
        raise PreflightFailure("STARTUP_CONTRACT_COMMIT must be a full orchestrator-supplied commit")
    resolved = subprocess.run(
        ["git", "rev-parse", f"{anchor}^{{commit}}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if resolved.returncode != 0 or resolved.stdout.strip() != anchor:
        raise PreflightFailure("STARTUP_CONTRACT_COMMIT did not resolve identically")
    if args.stage0_contract_digest != RATIFIED_DIGEST:
        raise AmendmentHalt("STAGE0_CONTRACT_DIGEST differs from the ratified contract")
    anchored_paths = [
        assert_anchored(repo, anchor, path) for path in (
            Path(__file__),
            source / "run_u2_boundary_adjudication.sh",
            source / "u2_stage0_runner.py",
            source / "u2_stage0_contract.py",
            source / "u2_stage0_sympy.py",
            source / "u2_production_runner.py",
            source / "u2_production_sympy.py",
            source / "u2_production_dual.wl",
            source / "u2_production_compare.py",
            source / "u2_production_mutations.py",
            source / "u2_code_closure_guard.py",
            source / "u2_containment.py",
        )
    ]
    reassert_stage0(repo, source, stage0, args.stage0_contract_digest)
    results = load_yaml(production / "production_results.yaml")
    if results["stage0_contract_digest"] != RATIFIED_DIGEST:
        raise PreflightFailure("production result is not bound to the ratified stage-0 digest")
    object_count, categories = validate_coverage(results)
    record = {
        "schema_version": "U2_A9_EXTERNAL_LEG_PREFLIGHT_V1",
        "status": "READY_FOR_EXTERNAL_LEG_NOT_A_LEG_RESULT",
        "leg": args.leg,
        "startup_contract_commit": anchor,
        "stage0_contract_digest": args.stage0_contract_digest,
        "stage0_bundle_recomputed": True,
        "environment_identity_recomputed": True,
        "evaluated_code_anchor_paths": anchored_paths,
        "sanitized_closure_required_for_external_process": True,
        "python_launch": "python3 -I with PYTHONNOUSERSITE=1",
        "wolfram_launch": "WolframKernel -noinit with installation-only $Path",
        "network_policy": "unshare-net; external leg must measure process-tree network=0",
        "coverage_object_count": object_count,
        "coverage_category_counts": categories,
        "coverage_granularity": "EVERY_OBJECT_INDIVIDUAL",
        "zero_cardinality_required_scopes": EXPECTED_ZERO_SCOPES,
        "coverage_exact": True,
        "external_leg_pass_claimed": False,
        "external_leg_result_owner": "orchestrator/fresh external verifier",
    }
    dump_yaml(Path(args.output).resolve(), record)
    print(f"U2_A9_PREFLIGHT_READY leg={args.leg} objects={object_count}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AmendmentHalt as failure:
        print(f"U2_A9_AMENDMENT_HALT {failure}", file=sys.stderr)
        raise SystemExit(42)
    except (PreflightFailure, KeyError, FileNotFoundError, subprocess.TimeoutExpired) as failure:
        print(f"U2_A9_PREFLIGHT_BLOCKED {failure}", file=sys.stderr)
        raise SystemExit(1)
