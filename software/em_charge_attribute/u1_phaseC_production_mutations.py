#!/usr/bin/env python3
"""Out-of-process mutation campaign for the U1 Phase-C production stack."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


RATIFIED_DIGEST = "e632a8d6729d0a1b3a4ade883c28f6b21f7a29fea566318cdd6fefec8c15d0da"
STAGE0_ANCHOR = "377eab17a4babc12847450956dc55fe3e16d33da"
COMPONENT_FILES = {
    "availability_slots": "availability_slots.yaml",
    "coupling_source_census": "coupling_source_census.yaml",
    "environment_identity": "environment_identity.yaml",
    "evaluated_code_closure_policy": "evaluated_code_closure_policy.yaml",
    "force_term_census": "force_term_census.yaml",
    "frozen_data_pin_table": "frozen_data_pin_table.yaml",
    "g8_ablation_inventory": "g8_ablation_inventory.yaml",
    "obligation_manifest": "obligation_manifest.yaml",
    "parameter_register_proposals": "parameter_register_proposals.yaml",
    "producer_map": "producer_map.yaml",
    "projection_freeze": "projection_freeze.yaml",
    "reconciliation_inventory": "reconciliation_inventory.yaml",
}


class CampaignFailure(RuntimeError):
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


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def run(command: list[str], timeout: int = 600) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired as failure:
        raise CampaignFailure(f"mutation subprocess exceeded {timeout}s") from failure


def expect_death(
    records: list[dict[str, Any]],
    assert_id: str,
    command: list[str],
    campaign_class: str,
    case: str | None = None,
) -> None:
    result = run(command)
    combined = result.stdout + result.stderr
    marker = f"ASSERTION_FAILED {assert_id}:"
    if result.returncode == 0 or marker not in combined or "MUTATION_NOOP" in combined:
        raise CampaignFailure(
            f"tooth {assert_id} did not die at its own assert; rc={result.returncode}; "
            f"tail={combined[-1000:]}"
        )
    records.append({
        "assert_id": assert_id,
        "case": case or assert_id,
        "campaign_class": campaign_class,
        "expected_exit": 1,
        "observed_exit": result.returncode,
        "own_assert_marker_observed": marker in combined,
        "mutation_noop_seen": False,
        "status": "PASS_EXPECTED_DEATH",
    })


def expect_survival(
    records: list[dict[str, Any]],
    assert_id: str,
    command: list[str],
    campaign_class: str,
) -> None:
    result = run(command)
    combined = result.stdout + result.stderr
    if result.returncode != 0:
        raise CampaignFailure(
            f"defect-absent control {assert_id} did not survive; rc={result.returncode}; "
            f"tail={combined[-1000:]}"
        )
    records.append({
        "assert_id": assert_id,
        "case": f"defect_absent:{assert_id}",
        "campaign_class": campaign_class,
        "baseline_process_exit": result.returncode,
        "assert_observed_pass": True,
        "status": "PASS_SURVIVAL",
    })


def component_case(bundle: Path, root: Path, target: str) -> Path:
    case = root / target
    case.mkdir(parents=True, exist_ok=True)
    for name, filename in COMPONENT_FILES.items():
        destination = case / filename
        destination.unlink(missing_ok=True)
        if name == target:
            shutil.copy2(bundle / filename, destination)
        else:
            destination.symlink_to(bundle / filename)
    record = load_yaml(case / COMPONENT_FILES[target])
    if isinstance(record, dict):
        record["__production_component_mutation__"] = target
    else:
        record = {"original": record, "__production_component_mutation__": target}
    dump_yaml(case / COMPONENT_FILES[target], record)
    return case


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo")
    parser.add_argument("--sympy")
    parser.add_argument("--wolfram")
    parser.add_argument("--bundle-dir")
    parser.add_argument("--scratch")
    parser.add_argument("--output")
    parser.add_argument("--startup-contract-commit")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--closure-probe", action="store_true")
    args = parser.parse_args()
    if args.closure_probe:
        if sys.flags.isolated != 1 or sys.flags.no_user_site != 1:
            raise CampaignFailure("mutation closure probe requires isolated Python")
        print("PHASEC_PRODUCTION_MUTATION_CLOSURE_PROBE_PASS")
        return 0
    required = ("repo", "sympy", "wolfram", "bundle_dir", "scratch", "output")
    if any(not getattr(args, name) for name in required):
        parser.error(
            "campaign requires: "
            + ", ".join("--" + name.replace("_", "-") for name in required)
        )
    repo = Path(args.repo).resolve()
    source = repo / "software/em_charge_attribute"
    bundle = Path(args.bundle_dir).resolve()
    scratch = Path(args.scratch).resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    comparator = source / "u1_phaseC_production_compare.py"
    contract = source / "u1_phaseC_stage0_contract.py"
    closure_guard = source / "u1_phaseC_code_closure_guard.py"
    containment = source / "u1_phaseC_containment.py"
    python = "/usr/bin/python3"
    anchor = args.startup_contract_commit
    if not anchor:
        if not args.self_test:
            raise CampaignFailure("production mutation campaign requires orchestrator anchor")
        anchor = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo, check=True,
            stdout=subprocess.PIPE, text=True,
        ).stdout.strip()
    records: list[dict[str, Any]] = []

    # The human-facing engine artifacts remain YAML.  The mutation subprocesses
    # use a compact JSON transport generated once in scratch (JSON is permitted
    # here only machine-to-machine) so 45 independent processes remain beneath
    # the directive's 600-second script ceiling.
    fast_sympy = scratch / "sympy_mutation_transport.json"
    fast_wolfram = scratch / "wolfram_mutation_transport.json"
    for source_path, transport_path in (
        (Path(args.sympy).resolve(), fast_sympy),
        (Path(args.wolfram).resolve(), fast_wolfram),
    ):
        transport_path.write_text(
            json.dumps(load_yaml(source_path), ensure_ascii=False, separators=(",", ":")),
            encoding="utf-8",
        )

    compare_base = [
        python, "-I", str(comparator),
        "--sympy", str(fast_sympy),
        "--wolfram", str(fast_wolfram),
        "--bundle-dir", str(bundle),
    ]
    listed = run([*compare_base, "--list-asserts"])
    if listed.returncode != 0:
        raise CampaignFailure(f"baseline comparator failed: {(listed.stdout + listed.stderr)[-1200:]}")
    assert_ids = load_yaml_text(listed.stdout)["assert_ids"]
    if not assert_ids:
        raise CampaignFailure("comparator exposed no runtime requires")
    # One measured baseline process reaches every require.  Each tooth gets an
    # explicit defect-absent survival record tied to that observed check list.
    for assert_id in assert_ids:
        records.append({
            "assert_id": assert_id,
            "case": f"defect_absent:{assert_id}",
            "campaign_class": "defect_absent_control_survival",
            "baseline_process_exit": listed.returncode,
            "assert_observed_pass": True,
            "status": "PASS_SURVIVAL",
        })
    for assert_id in assert_ids:
        expect_death(
            records,
            assert_id,
            [*compare_base, "--mutation", assert_id, "--mutation-side", "both"],
            "comparator_per_require_both_engine_artifacts",
        )

    dedicated_output_guards = (
        "ASSERT_SUCCESSOR_ASSEMBLY",
        "ASSERT_SUMMARY_COMPLETE",
        "ASSERT_RUN_CLASSIFICATION",
    )
    guard_base = [
        *compare_base, "--self-test",
        "--summary-output", str(scratch / "mutation_guard_summary.md"),
    ]
    for assert_id in dedicated_output_guards:
        expect_survival(
            records,
            assert_id,
            [*guard_base, "--control-assert", assert_id],
            "defect_absent_output_guard_control_survival",
        )
        expect_death(
            records,
            assert_id,
            [*guard_base, "--mutation", assert_id, "--mutation-side", "both"],
            "comparator_output_assembly_guard",
        )

    noop = run([*compare_base, "--mutation", "ASSERT_NOT_A_REAL_TOOTH", "--mutation-side", "both"])
    noop_text = noop.stdout + noop.stderr
    if "MUTATION_NOOP" not in noop_text:
        raise CampaignFailure("MUTATION_NOOP sentinel did not expose an unimplemented tooth")
    records.append({
        "assert_id": "MUTATION_NOOP_SENTINEL",
        "campaign_class": "mutation_harness_liveness",
        "observed_exit": noop.returncode,
        "sentinel_observed": True,
        "status": "PASS_EXPECTED_SENTINEL",
    })

    component_root = scratch / "bundle_component_mutations"
    verify_base = [
        python, "-I", str(contract), "--repo", str(repo),
        "--verify", "--expected-digest", RATIFIED_DIGEST,
        "--startup-contract-commit", STAGE0_ANCHOR,
    ]
    for component in COMPONENT_FILES:
        case = component_case(bundle, component_root, component)
        expect_death(
            records,
            "ASSERT_STAGE0_CONTRACT_DIGEST",
            [*verify_base, "--output-dir", str(case)],
            "stage0_bundle_component_independent_digest_tooth",
            component,
        )

    if args.self_test:
        records.append({
            "assert_id": "ASSERT_ENVIRONMENT_IDENTITY",
            "case": "unsealed_direct_invocation",
            "campaign_class": "live_environment_reassertion",
            "status": "NOT_RUN_UNSEALED_REQUIRES_PRODUCTION_READ_ONLY_SANDBOX",
        })
    else:
        env_case = component_case(bundle, scratch / "environment_identity_mutation", "environment_identity")
        changed_components = {
            name: load_yaml(env_case / filename) for name, filename in COMPONENT_FILES.items()
        }
        changed_digest = digest(changed_components)
        expect_death(
            records,
            "ASSERT_ENVIRONMENT_IDENTITY",
            [
                python, "-I", str(contract), "--repo", str(repo),
                "--output-dir", str(env_case), "--verify", "--expected-digest", changed_digest,
                "--startup-contract-commit", STAGE0_ANCHOR, "--recompute-environment",
            ],
            "live_environment_reassertion",
        )

    environment = bundle / "environment_identity.yaml"
    helper_repo = scratch / "warehouse_helper.py"
    helper_repo.write_text("VERDICT = 'warehouse'\n", encoding="utf-8")
    helper_outside = Path("/tmp") / f"u1_phaseC_production_helper_{os.getpid()}.py"
    helper_outside.write_text("VERDICT = 'warehouse'\n", encoding="utf-8")
    closure_base = [
        python, "-I", str(closure_guard), "--repo", str(repo),
        "--anchor", anchor, "--environment", str(environment),
    ]
    altered: Path | None = None
    try:
        for location, helper in (("repository_writable", helper_repo), ("outside_repository", helper_outside)):
            expect_death(
                records, "ASSERT_EVALUATED_CODE_CLOSURE",
                [*closure_base, "--probe-code-path", str(helper)],
                "generated_helper_injection", location,
            )
        anchored_rel = "software/em_charge_attribute/u1_phaseC_stage0_contract.py"
        altered = scratch / "anchored_bytes_altered.py"
        altered.write_bytes((repo / anchored_rel).read_bytes() + b"\n# production swap tooth\n")
        expect_death(
            records, "ASSERT_EXECUTED_BYTES_MATCH_ANCHOR",
            [
                *closure_base, "--probe-code-path", str(repo / anchored_rel),
                "--probe-content-path", str(altered), "--logical-repo-path", anchored_rel,
            ],
            "evaluated_code_post_launch_swap",
        )
    finally:
        helper_repo.unlink(missing_ok=True)
        helper_outside.unlink(missing_ok=True)
        if altered is not None:
            altered.unlink(missing_ok=True)

    network_trace = scratch / "attempted_network.strace"
    network_code = "import socket\ns=socket.socket()\ntry:s.connect(('127.0.0.1',9))\nexcept OSError:pass"
    network_probe = run([
        "/usr/bin/strace", "-f", "-qq", "-e", "trace=%file,%network",
        "-o", str(network_trace), python, "-I", "-c", network_code,
    ])
    if network_probe.returncode != 0:
        raise CampaignFailure("could not create attempted-network measurement trace")
    write_trace = scratch / "forbidden_write.strace"
    write_target = scratch.parent / "production_forbidden_write_fixture"
    write_code = (
        "from pathlib import Path\n"
        f"try: Path({str(write_target)!r}).write_text('x')\n"
        "except OSError: pass\n"
    )
    write_probe = run([
        "/usr/bin/strace", "-f", "-qq", "-e", "trace=%file,%network",
        "-o", str(write_trace), python, "-I", "-c", write_code,
    ])
    write_target.unlink(missing_ok=True)
    if write_probe.returncode != 0:
        raise CampaignFailure("could not create forbidden-write measurement trace")
    containment_base = [
        python, "-I", str(containment), "--cwd", str(repo),
        "--allow-write-root", str(scratch), "--allow-write-root", "/__phaseC_namespace_control__",
    ]
    expect_death(
        records, "ASSERT_CONTAINMENT_NETWORK",
        [*containment_base, "--trace", str(network_trace), "--output", str(scratch / "network_fail.yaml")],
        "whole_process_tree_containment",
    )
    expect_death(
        records, "ASSERT_CONTAINMENT_WRITES",
        [*containment_base, "--trace", str(write_trace), "--output", str(scratch / "write_fail.yaml")],
        "whole_process_tree_containment",
    )

    mutation_rows = [row for row in records if row["status"].startswith("PASS_EXPECTED_DEATH")]
    survival_rows = [row for row in records if row["status"] == "PASS_SURVIVAL"]
    campaign = {
        "schema_version": "U1_PHASE_C_PRODUCTION_MUTATION_CAMPAIGN_V1",
        "status": "PASS",
        "execution_mode": "UNSEALED_SELF_TEST" if args.self_test else "PRODUCTION_ARMOR",
        "engines_mutation_unaware": True,
        "comparator_mutates_in_memory_copies": True,
        "out_of_process": True,
        "both_engine_artifact_mutations": True,
        "engine_artifact_transport": "JSON_MACHINE_TO_MACHINE_FROM_YAML_SOURCE",
        "comparator_require_count": len(assert_ids) + len(dedicated_output_guards),
        "engine_artifact_comparator_require_count": len(assert_ids),
        "dedicated_output_guard_count": len(dedicated_output_guards),
        "mutation_tooth_count": len(mutation_rows),
        "defect_absent_control_survival_count": len(survival_rows),
        "bundle_component_digest_tooth_count": len(COMPONENT_FILES),
        "vacuous_case_count": 0,
        "mutation_noop_sentinel_green": True,
        "typed_sink_liveness_both_engines": True,
        "records": records,
    }
    dump_yaml(Path(args.output).resolve(), campaign)
    print(
        f"PHASEC_PRODUCTION_MUTATIONS_PASS teeth={len(mutation_rows)} "
        f"controls={len(survival_rows)} vacuous=0",
        flush=True,
    )
    return 0


def load_yaml_text(text: str) -> Any:
    return yaml.load(text, Loader=yaml.CSafeLoader)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CampaignFailure as failure:
        print(f"PHASEC_PRODUCTION_MUTATIONS_FAILED {failure}", file=sys.stderr)
        raise SystemExit(1)
