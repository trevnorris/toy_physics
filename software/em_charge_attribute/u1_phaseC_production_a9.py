#!/usr/bin/env python3
"""Preflight an orchestrator-owned Phase-C A9 external-verification leg.

This helper never marks an external leg PASS.  It reasserts the ratified bundle,
live environment, production anchor, and exact recomputation coverage immediately
before handing control to the external arbiter/agent/reviewer.
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


RATIFIED_DIGEST = "e632a8d6729d0a1b3a4ade883c28f6b21f7a29fea566318cdd6fefec8c15d0da"
STAGE0_ANCHOR = "377eab17a4babc12847450956dc55fe3e16d33da"
LEGS = ("arbiter", "fidelity", "adversarial_recompute", "read_only_review")


class PreflightFailure(RuntimeError):
    pass


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=140),
        encoding="utf-8",
    )


def assert_anchored(repo: Path, anchor: str, path: Path) -> str:
    rel = path.resolve().relative_to(repo).as_posix()
    result = subprocess.run(
        ["git", "cat-file", "blob", f"{anchor}:{rel}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise PreflightFailure(f"ASSERT_EVALUATED_CODE_CLOSURE: {rel} absent from anchor")
    actual = hashlib.sha256(path.read_bytes()).digest()
    expected = hashlib.sha256(result.stdout).digest()
    if actual != expected:
        raise PreflightFailure(f"ASSERT_EXECUTED_BYTES_MATCH_ANCHOR: {rel}")
    return rel


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
    stage0 = source / "reports/u1_body_dynamics_artifacts/stage_c_0_tilt_coupling_contract"
    production = source / "reports/u1_body_dynamics_artifacts/stage_c_1_tilt_coupling_production"
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
        raise PreflightFailure("STAGE0_CONTRACT_DIGEST differs from the ratified contract")
    anchored_paths = [
        assert_anchored(repo, anchor, Path(__file__)),
        assert_anchored(repo, anchor, source / "u1_phaseC_stage0_runner.py"),
        assert_anchored(repo, anchor, source / "run_u1_body_phaseC.sh"),
        assert_anchored(repo, anchor, source / "u1_phaseC_stage0_contract.py"),
        assert_anchored(repo, anchor, source / "u1_phaseC_production_sympy.py"),
        assert_anchored(repo, anchor, source / "u1_phaseC_production_dual.wl"),
        assert_anchored(repo, anchor, source / "u1_phaseC_production_compare.py"),
    ]
    runtime = source / "_scratch/u1_phaseC_production/a9_runtime"
    (runtime / "wolfram_userbase/Licensing").mkdir(parents=True, exist_ok=True)
    (runtime / "wolfram_base").mkdir(parents=True, exist_ok=True)
    nsswitch = runtime / "nsswitch.conf"
    nsswitch.write_text("hosts: files\n", encoding="utf-8")
    verify_command = [
        "/usr/bin/python3", "-I", str(source / "u1_phaseC_stage0_contract.py"),
        "--repo", str(repo), "--output-dir", str(stage0), "--verify",
        "--expected-digest", args.stage0_contract_digest,
        "--startup-contract-commit", STAGE0_ANCHOR, "--recompute-environment",
    ]
    verify = subprocess.run(
        [
            "/usr/bin/bwrap", "--ro-bind", "/", "/", "--dev-bind", "/dev", "/dev",
            "--proc", "/proc", "--unshare-net", "--bind", str(runtime), "/tmp",
            "--ro-bind", "/home/trevnorris/.Wolfram/Licensing", "/tmp/wolfram_userbase/Licensing",
            "--ro-bind", str(nsswitch), "/etc/nsswitch.conf", "--chdir", str(repo),
            "--setenv", "PYTHONNOUSERSITE", "1", "--setenv", "PYTHONDONTWRITEBYTECODE", "1",
            "--setenv", "TMPDIR", "/tmp", "--setenv", "WOLFRAM_USERBASE", "/tmp/wolfram_userbase",
            "--setenv", "WOLFRAM_BASE", "/tmp/wolfram_base", "--", *verify_command,
        ],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=600,
    )
    if verify.returncode != 0:
        raise PreflightFailure(f"stage0/environment reassertion failed: {(verify.stdout + verify.stderr)[-1200:]}")
    results = load_yaml(production / "production_results.yaml")
    coverage = results["A9_external_verification"]
    object_ids = coverage["object_ids"]
    if coverage["object_count"] != len(object_ids) or len(object_ids) != len(set(object_ids)):
        raise PreflightFailure("A9 exact coverage inventory is not one-to-one")
    coverage_map = coverage["coverage_map"]
    mapped_ids = [row["object_id"] for row in coverage_map]
    if (
        sorted(mapped_ids) != sorted(object_ids)
        or len(mapped_ids) != len(set(mapped_ids))
        or sum(coverage["coverage_category_counts"].values()) != len(object_ids)
        or not coverage["generated_category_union_exact"]
        or not coverage["coverage_map_exact"]
        or not coverage["computed_class_equivalence_all"]
    ):
        raise PreflightFailure("A9 generated-category union or identity recomputation map is not exact")
    record = {
        "schema_version": "U1_PHASE_C_A9_EXTERNAL_LEG_PREFLIGHT_V1",
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
        "network_policy": "unshare-net; measured process-tree network=0",
        "coverage_object_count": len(object_ids),
        "coverage_exact": True,
        "external_leg_pass_claimed": False,
        "external_leg_result_owner": "orchestrator/fresh external verifier",
    }
    dump_yaml(Path(args.output).resolve(), record)
    print(f"PHASEC_A9_PREFLIGHT_READY leg={args.leg} objects={len(object_ids)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (PreflightFailure, KeyError, FileNotFoundError, subprocess.TimeoutExpired) as failure:
        print(f"PHASEC_A9_PREFLIGHT_BLOCKED {failure}", file=sys.stderr)
        raise SystemExit(1)
