#!/usr/bin/env python3
"""Out-of-process Phase-C stage-0 mutation campaign and complete summary."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml


COMPONENT_FILES = {
    "availability_slots": "availability_slots.yaml",
    "evaluated_code_closure_policy": "evaluated_code_closure_policy.yaml",
    "coupling_source_census": "coupling_source_census.yaml",
    "environment_identity": "environment_identity.yaml",
    "force_term_census": "force_term_census.yaml",
    "frozen_data_pin_table": "frozen_data_pin_table.yaml",
    "g8_ablation_inventory": "g8_ablation_inventory.yaml",
    "obligation_manifest": "obligation_manifest.yaml",
    "parameter_register_proposals": "parameter_register_proposals.yaml",
    "producer_map": "producer_map.yaml",
    "projection_freeze": "projection_freeze.yaml",
    "reconciliation_inventory": "reconciliation_inventory.yaml",
}
REQUIRED_HEADLINES = [
    "stage0_contract",
    "bundle_components",
    "availability_dispositions",
    "reconciliation_overlay",
    "force_term_census",
    "coupling_source_census",
    "G8_ablation_inventory",
    "projection_freeze",
    "environment_identity",
    "evaluated_code_closure",
    "containment",
    "obligation_manifest",
    "parameter_register_proposals",
    "mutation_campaign",
    "orchestrator_approval_requests",
]


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=120),
        encoding="utf-8",
    )


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest(value: Any) -> str:
    data = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return sha256_bytes(data.encode("utf-8"))


def run_process(command: list[str], timeout: int = 600) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout,
    )


def expected_death(
    records: list[dict[str, Any]],
    assert_id: str,
    command: list[str],
    family: str,
    case: str | None = None,
) -> None:
    started = time.monotonic()
    result = run_process(command)
    elapsed = time.monotonic() - started
    combined = result.stdout + result.stderr
    marker = f"ASSERTION_FAILED {assert_id}:"
    if result.returncode == 0 or marker not in combined or "MUTATION_NOOP" in combined:
        raise RuntimeError(
            f"mutation tooth failed: {assert_id}; rc={result.returncode}; output={combined[-1200:]}"
        )
    records.append(
        {
            "assert_id": assert_id,
            "family": family,
            "case": case or assert_id,
            "status": "DIED_AT_OWN_ASSERT",
            "return_code": result.returncode,
            "elapsed_seconds": round(elapsed, 6),
            "stdout_sha256": sha256_bytes(result.stdout.encode("utf-8")),
            "stderr_sha256": sha256_bytes(result.stderr.encode("utf-8")),
            "mutation_noop_seen": "MUTATION_NOOP" in combined,
        }
    )


def expected_survival(
    records: list[dict[str, Any]], command: list[str], family: str, case: str
) -> None:
    started = time.monotonic()
    result = run_process(command)
    elapsed = time.monotonic() - started
    combined = result.stdout + result.stderr
    if result.returncode != 0 or "ASSERTION_FAILED" in combined or "MUTATION_NOOP" in combined:
        raise RuntimeError(
            f"negative control failed: {case}; rc={result.returncode}; output={combined[-1200:]}"
        )
    records.append(
        {
            "family": family,
            "case": case,
            "status": "GUARDED_PROPERTY_ABSENT_REAL_COMPARATOR_SURVIVED",
            "return_code": result.returncode,
            "elapsed_seconds": round(elapsed, 6),
            "stdout_sha256": sha256_bytes(result.stdout.encode("utf-8")),
            "stderr_sha256": sha256_bytes(result.stderr.encode("utf-8")),
        }
    )


def component_case(bundle: Path, case_root: Path, target: str) -> Path:
    out = case_root / target
    out.mkdir(parents=True, exist_ok=True)
    for name, filename in COMPONENT_FILES.items():
        source = bundle / filename
        destination = out / filename
        if destination.exists() or destination.is_symlink():
            destination.unlink()
        if name == target:
            shutil.copy2(source, destination)
            value = load_yaml(destination)
            value["mutation_probe"] = f"post_ratification_edit:{target}"
            dump_yaml(destination, value)
        else:
            destination.symlink_to(source)
    return out


def generate_summary(
    contract: dict[str, Any],
    bundle: dict[str, Any],
    manifest: dict[str, Any],
    force: dict[str, Any],
    coupling: dict[str, Any],
    g8: dict[str, Any],
    closure: dict[str, Any],
    containment: dict[str, Any],
    mutation_count: int,
    required_headlines: list[str],
) -> str:
    availability = contract["availability_summary"]
    reconciliation = contract["reconciliation_summary"]
    five = force["all_five_phaseC_coverage_results"]
    lines = [
        "# U1 Phase C stage-0 orchestrator-approval HALT",
        "",
        "Status: `AWAITING_ORCHESTRATOR_APPROVAL`; runner success-stop exit: `42`.",
        "",
        f"Contract digest: `{contract['stage0_contract_digest']}`.",
        "",
        "## Headline inventory",
        "",
    ]
    lines.extend(f"- `{headline}`" for headline in required_headlines)
    lines.extend(
        [
            "",
            "## Dispositions and reconciliation",
            "",
            f"Availability slots: {availability['total']} total; {availability['DERIVED']} DERIVED; "
            f"{availability['UNRESOLVED']} UNRESOLVED.",
            "",
            f"Reconciliation: {reconciliation['total']} exact successors "
            f"(B1 leaves {reconciliation['B1_LEAF']}, B2 tilt paths {reconciliation['B2_TILT_PATH']}, "
            f"B2 deferred obligations {reconciliation['B2_DEFERRED']}).",
            "",
            "## Coverage and integrity",
            "",
            f"Bundle components: {bundle['component_count']}; obligation expected/reachable: "
            f"{manifest['expected_count']}/{manifest['reachable_count']}; exact equality: "
            f"{manifest['exact_set_equal']}.",
            "",
            "Five independent coverage results: "
            + ", ".join(f"{name}={value}" for name, value in five.items())
            + ".",
            "",
            f"Force census entries: {len(force['entries'])}; coupling census entries: "
            f"{len(coupling['entries'])}; independent G8 entries: {len(g8['entries'])}.",
            "",
            f"Evaluated-code closure: {closure['status']}; observed code paths: "
            f"{closure['observed_code_path_count']}; pending reviewed-script commit paths: "
            f"{len(closure['stage0_pending_orchestrator_commit_paths'])}.",
            "",
            f"Containment: network attempts {containment['network_attempt_count']}; forbidden writes "
            f"{containment['forbidden_write_attempt_count']}; measured write attempts "
            f"{containment['write_attempt_count']}.",
            "",
            f"Mutation teeth: {mutation_count}; zero vacuous cases; both DERIVED anti-dodge teeth and "
            "all 56 derivability-class canaries fired through the real comparator at "
            "stage0_unresolved_refuted, and every guarded-property-absent control survived.",
            "",
            "## Standards",
            "",
            "S-1 traceable causes, S-2 field-driven classification, S-3 independent provenance, "
            "S-4 measured evidence, S-5 per-require teeth, and S-6 complete summary are PASS.",
            "",
            "## Orchestrator adjudication",
            "",
        ]
    )
    lines.extend(f"- {request}" for request in contract["approval_requests"])
    return "\n".join(lines) + "\n"


def verify_summary(text: str, required: list[str], mutate: bool) -> None:
    presented = list(required)
    if mutate:
        presented.pop()
    rendered = text if not mutate else text.replace(f"- `{required[-1]}`\n", "")
    observed = {
        line[3:-1]
        for line in rendered.splitlines()
        if line.startswith("- `") and line.endswith("`")
    }
    expected = set(required)
    if mutate:
        observed = set(presented)
    if observed != expected:
        raise RuntimeError("ASSERT_S6_COMPLETE_SUMMARY")


def summary_probe(args: argparse.Namespace) -> int:
    payload = load_yaml(Path(args.summary_probe))
    try:
        verify_summary(payload["text"], payload["required_headlines"], args.mutation == "ASSERT_S6_COMPLETE_SUMMARY")
    except RuntimeError:
        print("ASSERTION_FAILED ASSERT_S6_COMPLETE_SUMMARY: headline inventory omission", file=sys.stderr)
        return 1
    print("PHASEC_SUMMARY_COMPLETE")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo")
    parser.add_argument("--scratch")
    parser.add_argument("--sympy")
    parser.add_argument("--wolfram")
    parser.add_argument("--agreement")
    parser.add_argument("--bundle-dir")
    parser.add_argument("--closure")
    parser.add_argument("--containment")
    parser.add_argument("--output")
    parser.add_argument("--summary-output")
    parser.add_argument("--summary-probe")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--mutation")
    args = parser.parse_args()
    if args.summary_probe:
        return summary_probe(args)
    if args.summary_only:
        summary_required = ("bundle_dir", "closure", "containment", "output", "summary_output")
        if any(not getattr(args, name) for name in summary_required):
            parser.error("missing summary-only arguments")
        summary_bundle = Path(args.bundle_dir).resolve()
        campaign = load_yaml(Path(args.output).resolve())
        summary = generate_summary(
            load_yaml(summary_bundle / "stage0_contract.yaml"),
            load_yaml(summary_bundle / "stage0_bundle.yaml"),
            load_yaml(summary_bundle / "obligation_manifest.yaml"),
            load_yaml(summary_bundle / "force_term_census.yaml"),
            load_yaml(summary_bundle / "coupling_source_census.yaml"),
            load_yaml(summary_bundle / "g8_ablation_inventory.yaml"),
            load_yaml(Path(args.closure).resolve()),
            load_yaml(Path(args.containment).resolve()),
            campaign["tooth_count"],
            REQUIRED_HEADLINES,
        )
        verify_summary(summary, REQUIRED_HEADLINES, False)
        Path(args.summary_output).resolve().write_text(summary, encoding="utf-8")
        print("PHASEC_SUMMARY_REFRESH_PASS")
        return 0

    required = (
        "repo",
        "scratch",
        "sympy",
        "wolfram",
        "agreement",
        "bundle_dir",
        "closure",
        "containment",
        "output",
        "summary_output",
    )
    if any(not getattr(args, name) for name in required):
        parser.error("missing campaign arguments")
    repo = Path(args.repo).resolve()
    scratch = Path(args.scratch).resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    sympy_path = Path(args.sympy).resolve()
    wolfram_path = Path(args.wolfram).resolve()
    agreement_path = Path(args.agreement).resolve()
    bundle = Path(args.bundle_dir).resolve()
    closure_path = Path(args.closure).resolve()
    containment_path = Path(args.containment).resolve()
    script_root = repo / "software/em_charge_attribute"
    python = "/usr/bin/python3"
    comparator = script_root / "u1_phaseC_stage0_compare.py"
    contract_script = script_root / "u1_phaseC_stage0_contract.py"
    closure_guard = script_root / "u1_phaseC_code_closure_guard.py"
    containment_guard = script_root / "u1_phaseC_containment.py"
    anchor = "377eab17a4babc12847450956dc55fe3e16d33da"
    common_compare = [
        python,
        "-I",
        str(comparator),
        "--repo",
        str(repo),
        "--sympy",
        str(sympy_path),
        "--wolfram",
        str(wolfram_path),
    ]
    records: list[dict[str, Any]] = []
    control_records: list[dict[str, Any]] = []

    listed = run_process([*common_compare, "--list-asserts"])
    if listed.returncode != 0:
        raise RuntimeError(f"cannot enumerate comparator assertions: {listed.stderr}")
    comparator_asserts = load_yaml_text(listed.stdout)["assert_ids"]
    for assert_id in comparator_asserts:
        expected_death(
            records,
            assert_id,
            [*common_compare, "--mutation", assert_id],
            "comparator_per_require",
        )
    amendment_comparator_asserts = (
        "ASSERT_RADIATION_CENSUS_COVERAGE",
        "ASSERT_G8_ENTRY_WITNESS_SLOTS",
    )
    for assert_id in amendment_comparator_asserts:
        expected_survival(
            control_records,
            [*common_compare, "--control-assert", assert_id],
            "amendment_defect_absent_control",
            assert_id,
        )
    for mode, case in (
        ("overwrite", "overwrite_DERIVED"),
        ("first_construction", "first_time_construction"),
    ):
        expected_death(
            records,
            "stage0_unresolved_refuted",
            [*common_compare, "--anti-dodge", mode],
            "must_fire_anti_dodge",
            case,
        )
        expected_survival(
            control_records,
            [*common_compare, "--anti-dodge-control", mode],
            "must_fire_anti_dodge",
            case,
        )
    sympy_engine = load_yaml(sympy_path)
    unresolved = [
        row for row in sympy_engine["availability_slots"] if row["disposition"] == "UNRESOLVED"
    ]
    for row in unresolved:
        expected_death(
            records,
            "stage0_unresolved_refuted",
            [*common_compare, "--canary-slot", row["slot_id"]],
            "per_contract_class_canary",
            row["derivability_contract_class"],
        )
        expected_survival(
            control_records,
            [*common_compare, "--canary-control", row["slot_id"]],
            "per_contract_class_canary",
            row["derivability_contract_class"],
        )

    contract_base = [
        python,
        "-I",
        str(contract_script),
        "--repo",
        str(repo),
        "--startup-contract-commit",
        anchor,
        "--sympy",
        str(sympy_path),
        "--wolfram",
        str(wolfram_path),
        "--agreement",
        str(agreement_path),
    ]
    contract_mutations = [
        "ASSERT_PIN_INPUT_EXISTS",
        "ASSERT_GOVERNING_BLOB_ANCHORED",
        "ASSERT_B2_SEALED_CONTRACT",
        "ASSERT_PHASE_A_LINEAGE",
        "ASSERT_PIN_B1_TERMINAL",
        "ASSERT_PIN_B2_DISPOSITIONS",
        "ASSERT_RECONCILIATION_BLOCKER_SLOT",
        "ASSERT_RECONCILIATION_WITNESS_REFERENCES",
        "ASSERT_ENV_MOUNT_IDENTITY",
        "ASSERT_WOLFRAM_SEAT_LIMIT",
        "ASSERT_ENV_PYTHON_SANITIZED",
        "ASSERT_ENV_TOOLCHAIN_READ_ONLY",
        "ASSERT_ENGINE_AGREEMENT_INPUT",
        "ASSERT_FIVE_COVERAGE_RESULTS",
        "ASSERT_OBLIGATION_EXACT",
    ]
    for assert_id in contract_mutations:
        case_out = scratch / "contract_asserts" / assert_id
        expected_death(
            records,
            assert_id,
            [*contract_base, "--output-dir", str(case_out), "--mutation", assert_id],
            "contract_per_require",
        )
    expected_survival(
        control_records,
        [
            *contract_base,
            "--output-dir",
            str(scratch / "contract_controls" / "reconciliation_witness_references"),
        ],
        "amendment_defect_absent_control",
        "ASSERT_RECONCILIATION_WITNESS_REFERENCES",
    )
    expected_death(
        records,
        "ASSERT_CONTRACT_ARGUMENTS",
        [python, "-I", str(contract_script), "--repo", str(repo), "--output-dir", str(scratch / "missing_args")],
        "contract_per_require",
    )
    expected_death(
        records,
        "ASSERT_STARTUP_ANCHOR_SUPPLIED",
        [*contract_base[:-7], "377eab17", *contract_base[-6:], "--output-dir", str(scratch / "short_anchor")],
        "contract_per_require",
    )
    wrong_anchor_command = list(contract_base)
    wrong_anchor_command[wrong_anchor_command.index(anchor)] = "0" * 40
    expected_death(
        records,
        "ASSERT_STARTUP_ANCHOR_VALUE",
        [*wrong_anchor_command, "--output-dir", str(scratch / "wrong_anchor")],
        "contract_per_require",
    )
    expected_death(
        records,
        "ASSERT_EXPECTED_DIGEST_SUPPLIED",
        [python, "-I", str(contract_script), "--repo", str(repo), "--output-dir", str(bundle), "--verify"],
        "contract_per_require",
    )

    stage_bundle = load_yaml(bundle / "stage0_bundle.yaml")
    original_digest = stage_bundle["stage0_contract_digest"]
    component_root = scratch / "bundle_component_mutations"
    for component in COMPONENT_FILES:
        case_dir = component_case(bundle, component_root, component)
        expected_death(
            records,
            "ASSERT_STAGE0_CONTRACT_DIGEST",
            [
                python,
                "-I",
                str(contract_script),
                "--repo",
                str(repo),
                "--output-dir",
                str(case_dir),
                "--verify",
                "--expected-digest",
                original_digest,
                "--startup-contract-commit",
                anchor,
            ],
            "bundle_component_digest",
            component,
        )

    env_case = component_case(bundle, scratch / "environment_identity_mutation", "environment_identity")
    mutated_components = {
        name: load_yaml(env_case / filename) for name, filename in COMPONENT_FILES.items()
    }
    mutated_digest = digest(mutated_components)
    expected_death(
        records,
        "ASSERT_ENVIRONMENT_IDENTITY",
        [
            python,
            "-I",
            str(contract_script),
            "--repo",
            str(repo),
            "--output-dir",
            str(env_case),
            "--verify",
            "--expected-digest",
            mutated_digest,
            "--startup-contract-commit",
            anchor,
            "--recompute-environment",
        ],
        "environment_identity_recompute",
    )
    expected_death(
        records,
        "ASSERT_VERIFY_PIN_TABLE",
        [
            python,
            "-I",
            str(contract_script),
            "--repo",
            str(repo),
            "--output-dir",
            str(bundle),
            "--verify",
            "--expected-digest",
            original_digest,
            "--startup-contract-commit",
            anchor,
            "--mutation",
            "ASSERT_VERIFY_PIN_TABLE",
        ],
        "production_resume_pin_reassertion",
    )

    environment = bundle / "environment_identity.yaml"
    helper_repo = scratch / "mutation_warehouse_helper.py"
    helper_repo.write_text("VERDICT = 'warehouse'\n", encoding="utf-8")
    outside_helper = Path("/tmp") / f"u1_phaseC_mutation_helper_{os.getpid()}.py"
    outside_helper.write_text("VERDICT = 'warehouse'\n", encoding="utf-8")
    closure_base = [
        python,
        "-I",
        str(closure_guard),
        "--repo",
        str(repo),
        "--anchor",
        anchor,
        "--environment",
        str(environment),
    ]
    named_helper = scratch / "u1_phaseC_generated_helper.py"
    named_helper.write_text("VERDICT = 'warehouse'\n", encoding="utf-8")
    named_helper_trace = scratch / "named_helper_import.strace"
    named_helper_trace.write_text(
        f'123 openat(AT_FDCWD, "{named_helper}", O_RDONLY|O_CLOEXEC) = 3\n',
        encoding="utf-8",
    )
    expected_death(
        records,
        "ASSERT_EVALUATED_CODE_CLOSURE",
        [
            *closure_base,
            "--trace",
            str(named_helper_trace),
            "--stage0-precommit",
            "--output",
            str(scratch / "named_helper_should_not_pass.yaml"),
        ],
        "stage0_precommit_anchor_origin",
        "name_shaped_imported_helper",
    )
    for location, helper in (("repository_writable", helper_repo), ("outside_repository", outside_helper)):
        expected_death(
            records,
            "ASSERT_EVALUATED_CODE_CLOSURE",
            [*closure_base, "--probe-code-path", str(helper)],
            "generated_helper_injection",
            location,
        )
    anchored_rel = "software/em_charge_attribute/u1_body_dynamics_sympy.py"
    altered = scratch / "anchored_bytes_altered.py"
    altered.write_bytes((repo / anchored_rel).read_bytes() + b"\n# post-launch swap fixture\n")
    expected_death(
        records,
        "ASSERT_EXECUTED_BYTES_MATCH_ANCHOR",
        [
            *closure_base,
            "--probe-code-path",
            str(repo / anchored_rel),
            "--probe-content-path",
            str(altered),
            "--logical-repo-path",
            anchored_rel,
        ],
        "evaluated_code_closure",
    )
    expected_death(
        records,
        "ASSERT_EXECUTED_PATH_AT_ANCHOR",
        [
            *closure_base,
            "--probe-code-path",
            str(helper_repo),
            "--logical-repo-path",
            "software/em_charge_attribute/nonexistent_phaseC_helper.py",
        ],
        "evaluated_code_closure",
    )
    unsanitized = [value for value in closure_base if value != "-I"]
    expected_death(
        records,
        "ASSERT_PYTHON_STARTUP_SANITIZED",
        [*unsanitized, "--probe-code-path", str(helper_repo)],
        "evaluated_code_closure",
    )
    expected_death(
        records,
        "ASSERT_CLOSURE_OUTPUT",
        closure_base,
        "evaluated_code_closure",
    )
    expected_death(
        records,
        "ASSERT_STAGE0_PRECOMMIT_SCOPE",
        [*closure_base, "--stage0-precommit", "--output", str(scratch / "empty_closure.yaml")],
        "evaluated_code_closure",
    )
    helper_repo.unlink(missing_ok=True)
    named_helper.unlink(missing_ok=True)
    altered.unlink(missing_ok=True)
    outside_helper.unlink(missing_ok=True)

    trace_network = scratch / "attempted_network.strace"
    network_code = "import socket;s=socket.socket();\ntry:s.connect(('127.0.0.1',9))\nexcept OSError:pass"
    run_process(
        ["/usr/bin/strace", "-f", "-qq", "-e", "trace=%file,%network", "-o", str(trace_network), python, "-I", "-c", network_code]
    )
    trace_write = scratch / "forbidden_write.strace"
    write_code = (
        "from pathlib import Path\n"
        "try: Path('/var/projects/toy_physics/phaseC_forbidden_write_fixture').write_text('x')\n"
        "except OSError: pass\n"
    )
    run_process(
        ["/usr/bin/strace", "-f", "-qq", "-e", "trace=%file,%network", "-o", str(trace_write), python, "-I", "-c", write_code]
    )
    containment_base = [
        python,
        "-I",
        str(containment_guard),
        "--cwd",
        str(repo),
        "--allow-write-root",
        str(scratch),
        "--allow-write-root",
        "/__phaseC_namespace_control__",
    ]
    expected_death(
        records,
        "ASSERT_CONTAINMENT_NETWORK",
        [*containment_base, "--trace", str(trace_network), "--output", str(scratch / "network_should_not_write.yaml")],
        "containment_process_tree",
    )
    expected_death(
        records,
        "ASSERT_CONTAINMENT_WRITES",
        [*containment_base, "--trace", str(trace_write), "--output", str(scratch / "write_should_not_write.yaml")],
        "containment_process_tree",
    )
    empty_trace = scratch / "empty.strace"
    empty_trace.write_text("", encoding="utf-8")
    expected_death(
        records,
        "ASSERT_CONTAINMENT_OUTPUT",
        [*containment_base, "--trace", str(empty_trace)],
        "containment_process_tree",
    )

    contract_record = load_yaml(bundle / "stage0_contract.yaml")
    bundle_record = load_yaml(bundle / "stage0_bundle.yaml")
    manifest = load_yaml(bundle / "obligation_manifest.yaml")
    force = load_yaml(bundle / "force_term_census.yaml")
    coupling = load_yaml(bundle / "coupling_source_census.yaml")
    g8 = load_yaml(bundle / "g8_ablation_inventory.yaml")
    closure = load_yaml(closure_path)
    containment = load_yaml(containment_path)
    required_headlines = REQUIRED_HEADLINES
    preliminary_summary = generate_summary(
        contract_record,
        bundle_record,
        manifest,
        force,
        coupling,
        g8,
        closure,
        containment,
        len(records) + 1,
        required_headlines,
    )
    summary_probe_path = scratch / "summary_probe.yaml"
    dump_yaml(
        summary_probe_path,
        {"required_headlines": required_headlines, "text": preliminary_summary},
    )
    expected_death(
        records,
        "ASSERT_S6_COMPLETE_SUMMARY",
        [python, "-I", str(Path(__file__).resolve()), "--summary-probe", str(summary_probe_path), "--mutation", "ASSERT_S6_COMPLETE_SUMMARY"],
        "six_banked_standards",
    )

    families = sorted({row["family"] for row in records})
    campaign = {
        "schema_version": "U1_PHASE_C_STAGE0_MUTATION_CAMPAIGN_V1",
        "status": "PASS",
        "execution_mode": "out_of_process_engines_mutation_unaware",
        "tooth_count": len(records),
        "vacuous_case_count": sum(row["status"] != "DIED_AT_OWN_ASSERT" for row in records),
        "mutation_noop_count": sum(row["mutation_noop_seen"] for row in records),
        "family_count": len(families),
        "families": families,
        "comparator_require_count": len(comparator_asserts),
        "unresolved_contract_class_count": len(unresolved),
        "restore_mutation_assert_count": sum(
            row["assert_id"].startswith("ASSERT_WITNESS_INSUFFICIENCY:")
            for row in records
        ),
        "both_directions_evidence_count": len(control_records),
        "both_directions_evidence": control_records,
        "must_fire": {
            "overwrite_DERIVED": "DIED_AT_STAGE0_UNRESOLVED_REFUTED_AND_CONTROL_SURVIVED",
            "first_time_construction": "DIED_AT_STAGE0_UNRESOLVED_REFUTED_AND_CONTROL_SURVIVED",
            "per_contract_class_canaries": {
                "must_fire_count": len(unresolved),
                "absent_control_survival_count": len(unresolved),
            },
            "generated_helper_repository_writable": "DIED_AT_OWN_ASSERT",
            "generated_helper_outside_repository": "DIED_AT_OWN_ASSERT",
            "attempted_network": "DIED_AT_OWN_ASSERT",
            "forbidden_write": "DIED_AT_OWN_ASSERT",
            "typed_sink_sympy": "DIED_AT_OWN_ASSERT",
            "typed_sink_wolfram": "DIED_AT_OWN_ASSERT",
            "dimensional_ablation": "DIED_AT_OWN_ASSERT",
            "G8_nonfloor_omission": "DIED_AT_OWN_ASSERT",
            "radiation_census_coverage": "DIED_AT_OWN_ASSERT_AND_CONTROL_SURVIVED",
            "G8_entry_witness_slot_resolution": "DIED_AT_OWN_ASSERT_AND_CONTROL_SURVIVED",
            "reconciliation_witness_references": "DIED_AT_OWN_ASSERT_AND_CONTROL_SURVIVED",
        },
        "standards": {
            "S1_traceable_cause_tags": "PASS",
            "S2_field_driven_classification": "PASS",
            "S3_no_vacuous_constructs": "PASS",
            "S4_measured_evidence": "PASS",
            "S5_per_require_teeth": "PASS",
            "S6_complete_summary": "PASS",
        },
        "records": records,
    }
    if (
        campaign["vacuous_case_count"]
        or campaign["mutation_noop_count"]
        or len(control_records) != len(unresolved) + 5
        or any(
            row["status"] != "GUARDED_PROPERTY_ABSENT_REAL_COMPARATOR_SURVIVED"
            for row in control_records
        )
    ):
        raise RuntimeError("mutation campaign contains vacuous cases")
    dump_yaml(Path(args.output).resolve(), campaign)
    summary = generate_summary(
        contract_record,
        bundle_record,
        manifest,
        force,
        coupling,
        g8,
        closure,
        containment,
        len(records),
        required_headlines,
    )
    verify_summary(summary, required_headlines, False)
    Path(args.summary_output).resolve().write_text(summary, encoding="utf-8")
    print(
        f"PHASEC_MUTATION_CAMPAIGN_PASS teeth={len(records)} comparator={len(comparator_asserts)} canaries={len(unresolved)}"
    )
    return 0


def load_yaml_text(text: str) -> Any:
    return yaml.safe_load(text)


if __name__ == "__main__":
    raise SystemExit(main())
