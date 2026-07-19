#!/usr/bin/env python3
"""Out-of-process mutation campaign for the U2 production stack."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

import yaml


RATIFIED_DIGEST = "b01a1821e908589c3698512bbb9aff874b721af6dcbfa1c3b8b1f8d33119b32b"


class CampaignFailure(RuntimeError):
    pass


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=160), encoding="utf-8")


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command_digest(command: list[str]) -> str:
    return hashlib.sha256(json.dumps(command, separators=(",", ":")).encode()).hexdigest()


def run(command: list[str], timeout: int = 600) -> tuple[subprocess.CompletedProcess[str], float]:
    started = time.monotonic()
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
    return result, time.monotonic() - started


def require_death(command: list[str], assert_id: str, mutation_id: str) -> dict[str, Any]:
    result, elapsed = run(command)
    combined = result.stdout + result.stderr
    marker = f"ASSERTION_FAILED {assert_id}:"
    if result.returncode == 0 or marker not in combined or "ASSERTION_FAILED MUTATION_NOOP:" in combined:
        raise CampaignFailure(f"{mutation_id} did not die at {assert_id}: rc={result.returncode} {combined[-1000:]}")
    return {
        "mutation_id": mutation_id, "expected_assert_id": assert_id,
        "status": "DIED_AT_OWN_ASSERT", "return_code": result.returncode,
        "elapsed_seconds": round(elapsed, 6), "command_sha256": command_digest(command),
        "stdout_sha256": hashlib.sha256(result.stdout.encode()).hexdigest(),
        "stderr_sha256": hashlib.sha256(result.stderr.encode()).hexdigest(),
    }


def require_survival(command: list[str], control_id: str) -> dict[str, Any]:
    result, elapsed = run(command)
    combined = result.stdout + result.stderr
    if result.returncode != 0 or "ASSERTION_FAILED" in combined or "MUTATION_NOOP" in combined:
        raise CampaignFailure(f"control {control_id} failed: rc={result.returncode} {combined[-1000:]}")
    return {
        "control_id": control_id, "status": "DEFECT_ABSENT_REAL_CHECK_SURVIVED",
        "elapsed_seconds": round(elapsed, 6), "command_sha256": command_digest(command),
        "stdout_sha256": hashlib.sha256(result.stdout.encode()).hexdigest(),
        "stderr_sha256": hashlib.sha256(result.stderr.encode()).hexdigest(),
    }


def main() -> int:
    if "--closure-probe" in sys.argv:
        print("U2_PRODUCTION_MUTATION_CLOSURE_PROBE_PASS")
        return 0
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True); parser.add_argument("--sympy", required=True)
    parser.add_argument("--wolfram", required=True); parser.add_argument("--bundle-dir", required=True)
    parser.add_argument("--scratch", required=True); parser.add_argument("--output", required=True)
    parser.add_argument("--startup-contract-commit"); parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--closure-probe", action="store_true")
    args = parser.parse_args()
    repo = Path(args.repo).resolve(); source = repo / "software/em_charge_attribute"
    anchor = args.startup_contract_commit
    if not re.fullmatch(r"[0-9a-f]{40}", anchor or "") or anchor == "HEAD":
        raise CampaignFailure("production mutation campaign requires full orchestrator anchor, never HEAD")
    resolved = subprocess.run(
        ["git", "rev-parse", f"{anchor}^{{commit}}"], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if resolved.returncode != 0 or resolved.stdout.strip() != anchor:
        raise CampaignFailure("production mutation campaign anchor did not resolve identically")
    scratch = Path(args.scratch).resolve(); scratch.mkdir(parents=True, exist_ok=True)
    comparator = source / "u2_production_compare.py"
    probe = scratch / "comparator_mutation_probe.yaml"
    python = "/usr/bin/python3"
    python_prefix = [python, "-I", "-B"]
    probe_command = [
        *python_prefix, str(comparator), "--sympy", str(Path(args.sympy).resolve()),
        "--wolfram", str(Path(args.wolfram).resolve()), "--bundle-dir", str(Path(args.bundle_dir).resolve()),
        "--campaign-probe", "--campaign-output", str(probe),
    ]
    probe_control = require_survival(probe_command, "CONTROL_COMPARATOR_MUTATION_PROBE")
    campaign = load_yaml(probe)
    if campaign["status"] != "PASS" or campaign["vacuous_case_count"] or campaign["mutation_noop_count"]:
        raise CampaignFailure("comparator mutation probe is not green/nonvacuous")
    records = list(campaign["records"])
    controls = list(campaign["defect_absent_controls"])
    armor_records = []
    armor_controls = [probe_control]

    generation_output = scratch / "sympy_generation_control.yaml"
    generation_base = [
        *python_prefix, str(source / "u2_production_sympy.py"),
        "--repo", str(repo), "--bundle-dir", str(Path(args.bundle_dir).resolve()),
        "--stage0-contract-digest", RATIFIED_DIGEST, "--output", str(generation_output),
    ]
    armor_controls.append(require_survival(
        generation_base, "CONTROL_INDEPENDENT_GENERATION_CHECKS_CLEAN",
    ))
    generation_teeth = (
        ("TOOTH_EXCHANGE_EXPECTED_SIGNATURE_GENERATION", "ASSERT_EXCHANGE_SIGNATURE_GENERATION"),
        ("TOOTH_NO_PAIRING_CERTIFICATE_GENERATION", "ASSERT_NO_PAIRING_CERTIFICATE_GENERATION"),
        ("TOOTH_CLOSURE_CENSUS_GENERATION", "ASSERT_CLOSURE_CENSUS_GENERATION"),
        ("TOOTH_CLOSURE_TOTAL_GENERATION", "ASSERT_CLOSURE_TOTAL_GENERATION"),
        ("TOOTH_PHYSICS_RECORD_INVARIANCE", "ASSERT_PHYSICS_RECORD_INVARIANCE"),
    )
    for mutation_id, assert_id in generation_teeth:
        armor_records.append(require_death(
            [*generation_base, "--generation-mutation", mutation_id], assert_id, mutation_id,
        ))

    wolfram_output = scratch / "wolfram_generation_control.yaml"
    wolfram_base = [
        "/usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel",
        "-noinit", "-noprompt", "-script", str(source / "u2_production_dual.wl"),
        "--repo", str(repo), "--bundle-dir", str(Path(args.bundle_dir).resolve()),
        "--stage0-contract-digest", RATIFIED_DIGEST, "--output", str(wolfram_output),
    ]
    armor_controls.append(require_survival(
        wolfram_base, "CONTROL_WOLFRAM_INDEPENDENT_GENERATION_CHECKS_CLEAN",
    ))
    for mutation_id, assert_id in generation_teeth:
        armor_records.append(require_death(
            [*wolfram_base, "--generation-mutation", mutation_id],
            assert_id, f"WOLFRAM_{mutation_id}",
        ))

    empty_trace = scratch / "empty_containment.strace"
    empty_trace.write_text("", encoding="utf-8")
    containment_base = [
        *python_prefix, str(source / "u2_containment.py"), "--trace", str(empty_trace),
        "--cwd", str(repo), "--allow-write-root", str(scratch), "--output", str(scratch / "containment_control.yaml"),
    ]
    armor_records.append(require_death([*containment_base, "--probe-network"], "ASSERT_CONTAINMENT_NETWORK", "TOOTH_CONTAINMENT_NETWORK"))
    armor_records.append(require_death([*containment_base, "--probe-write"], "ASSERT_CONTAINMENT_WRITES", "TOOTH_CONTAINMENT_WRITE"))
    armor_controls.append(require_survival(containment_base, "CONTROL_CONTAINMENT_CLEAN"))

    verify_prefix: list[str] = []
    if args.dry_run:
        dry_runtime = scratch / "dry_runtime"
        for path in (dry_runtime, dry_runtime / "wolfram_userbase/Licensing", dry_runtime / "wolfram_base"):
            path.mkdir(parents=True, exist_ok=True)
        dry_nsswitch = dry_runtime / "nsswitch.conf"
        dry_nsswitch.write_text("hosts: files\n", encoding="utf-8")
        verify_prefix = [
            "/usr/bin/bwrap", "--ro-bind", "/", "/", "--dev-bind", "/dev", "/dev",
            "--proc", "/proc", "--unshare-net", "--bind", str(dry_runtime), "/tmp",
            "--ro-bind", "/home/trevnorris/.Wolfram/Licensing", "/tmp/wolfram_userbase/Licensing",
            "--ro-bind", str(dry_nsswitch), "/etc/nsswitch.conf", "--chdir", str(repo),
            "--setenv", "PYTHONNOUSERSITE", "1", "--setenv", "PYTHONDONTWRITEBYTECODE", "1",
            "--setenv", "TMPDIR", "/tmp", "--setenv", "WOLFRAM_USERBASE", "/tmp/wolfram_userbase",
            "--setenv", "WOLFRAM_BASE", "/tmp/wolfram_base", "--",
        ]
    verify_clean = [*verify_prefix, *python_prefix,
        str(source / "u2_stage0_contract.py"), "--repo", str(repo),
        "--output-dir", str(Path(args.bundle_dir).resolve()), "--startup-contract-commit", anchor,
        "--verify", "--expected-digest", RATIFIED_DIGEST, "--recompute-environment",
    ]
    armor_controls.append(require_survival(verify_clean, "CONTROL_PRODUCTION_STAGE0_AND_ENVIRONMENT"))
    armor_records.append(require_death(
        [*verify_clean[:-2], "0" * 64, "--recompute-environment"],
        "ASSERT_STAGE0_CONTRACT_DIGEST", "TOOTH_PRODUCTION_STAGE0_DIGEST",
    ))
    bundle = load_yaml(Path(args.bundle_dir) / "stage0_bundle.yaml")
    component_names = [row["name"] for row in bundle["components"]]
    if len(component_names) != 15 or len(component_names) != len(set(component_names)):
        raise CampaignFailure("stage-0 component inventory is not the ratified 15-component set")
    for name in component_names:
        armor_records.append(require_death(
            [*verify_clean, "--mutation", f"BUNDLE_COMPONENT:{name}"],
            "ASSERT_STAGE0_CONTRACT_DIGEST", f"TOOTH_PRODUCTION_BUNDLE_COMPONENT:{name}",
        ))
    armor_records.append(require_death(
        [*verify_clean, "--mutation", "TOOTH_ENVIRONMENT_IDENTITY"],
        "ASSERT_ENVIRONMENT_IDENTITY", "TOOTH_PRODUCTION_ENVIRONMENT_IDENTITY",
    ))

    if not args.dry_run:
        helper = scratch / "generated_helper.py"
        altered: Path | None = None
        try:
            helper.write_text("value = 1\n", encoding="utf-8")
            closure_base = [
                *python_prefix, str(source / "u2_code_closure_guard.py"), "--repo", str(repo),
                "--anchor", anchor, "--environment", str(Path(args.bundle_dir) / "environment_identity.yaml"),
            ]
            armor_records.append(require_death(
                [*closure_base, "--probe-code-path", str(helper)],
                "ASSERT_EVALUATED_CODE_CLOSURE", "TOOTH_HELPER_REPOSITORY_WRITABLE",
            ))
            with tempfile.TemporaryDirectory(prefix="u2_production_helper_", dir="/tmp") as outside_directory:
                outside_helper = Path(outside_directory) / "generated_helper.py"
                outside_helper.write_text("value = 1\n", encoding="utf-8")
                armor_records.append(require_death(
                    [*closure_base, "--probe-code-path", str(outside_helper)],
                    "ASSERT_EVALUATED_CODE_CLOSURE", "TOOTH_HELPER_OUTSIDE_REPOSITORY",
                ))
            anchored_rel = "software/em_charge_attribute/u2_production_sympy.py"
            altered = scratch / "altered_anchored_helper.py"
            altered.write_bytes((repo / anchored_rel).read_bytes() + b"\n# altered generated-helper tooth\n")
            armor_records.append(require_death(
                [
                    *closure_base, "--probe-code-path", str(repo / anchored_rel),
                    "--probe-content-path", str(altered), "--logical-repo-path", anchored_rel,
                ], "ASSERT_EXECUTED_BYTES_MATCH_ANCHOR", "TOOTH_HELPER_MASQUERADES_AS_ANCHORED"
            ))
            armor_controls.append(require_survival(
                [*closure_base, "--probe-code-path", str(repo / anchored_rel)], "CONTROL_ANCHORED_PRODUCTION_CODE"
            ))
        finally:
            helper.unlink(missing_ok=True)
            if altered is not None:
                altered.unlink(missing_ok=True)

    records.extend(armor_records)
    # Every armor death has a clean control; the comparator's 440 controls share
    # one full clean execution and containment/closure controls are explicit.
    controls.extend({
        "mutation_id": row["mutation_id"], "status": "DEFECT_ABSENT_REAL_CHECK_SURVIVED",
        "shared_control_execution": True,
    } for row in armor_records)
    result = {
        "schema_version": "U2_PRODUCTION_MUTATION_CAMPAIGN_V1", "status": "PASS",
        "execution_mode": "out_of_process_comparator_campaign_probe",
        "dry_run": args.dry_run, "tooth_count": len(records),
        "defect_absent_control_count": len(controls), "vacuous_case_count": 0,
        "unexpected_mutation_noop_count": 0, "control_failure_count": 0,
        "records": records, "defect_absent_controls": controls,
        "armor_control_executions": armor_controls,
        "comparator_probe_sha256": sha256_path(probe),
        "sealed_anchor_teeth_executed": not args.dry_run,
    }
    dump_yaml(Path(args.output).resolve(), result)
    print(f"U2_PRODUCTION_MUTATIONS_PASS teeth={len(records)} controls={len(controls)} dry_run={args.dry_run}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (CampaignFailure, subprocess.TimeoutExpired, KeyError, FileNotFoundError) as failure:
        print(f"U2_PRODUCTION_MUTATIONS_BLOCKED {failure}", file=sys.stderr)
        raise SystemExit(1)
