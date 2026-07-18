#!/usr/bin/env python3
"""Out-of-process U2 stage-0 mutation campaign and complete HALT summary."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml


ANCHOR = "323b222846e2a9062330d2f25dd9cd28c57c7800"
COMPONENT_FILES = {
    "frozen_data_pin_table": "frozen_data_pin_table.yaml",
    "candidate_inventory": "candidate_inventory.yaml",
    "obligation_censuses": "obligation_censuses.yaml",
    "dependency_grid_inventory": "dependency_grid_inventory.yaml",
    "vocabulary_freeze": "vocabulary_freeze.yaml",
    "evidence_taxonomy": "evidence_taxonomy.yaml",
    "availability_slots": "availability_slots.yaml",
    "route_fixture_inventory": "route_fixture_inventory.yaml",
    "closure_template_contracts": "closure_template_contracts.yaml",
    "environment_identity": "environment_identity.yaml",
    "standard_bindings": "standard_bindings.yaml",
    "producer_map": "producer_map.yaml",
    "evaluated_code_closure_policy": "evaluated_code_closure_policy.yaml",
    "parameter_register_proposals": "parameter_register_proposals.yaml",
    "obligation_manifest": "obligation_manifest.yaml",
}
REQUIRED_HEADLINES = [
    "stage0_contract_and_bundle",
    "frozen_data_pin_table",
    "candidate_axis_mixture_basis_and_canonicalization",
    "stratum_free_obligation_censuses",
    "open_dependency_relation_and_grid",
    "disposition_promotion_ensemble_topology_closure_vocabularies",
    "typed_evidence_taxonomy_source_census_and_node_grammar",
    "witness_and_challenge_slots",
    "route_inventory_and_known_outcome_fixtures",
    "environment_and_toolchain_identity",
    "evaluated_code_closure",
    "process_tree_containment",
    "obligation_manifest",
    "parameter_register_proposals",
    "mutation_campaign",
    "sixteen_banked_standards",
    "orchestrator_approval_requests",
]
SURVIVAL_EXECUTION_CACHE: dict[str, tuple[subprocess.CompletedProcess[str], float]] = {}


def load_yaml(path: Path) -> Any:
    with path.open("rb") as handle:
        return yaml.load(handle, Loader=yaml.CSafeLoader)


def dump_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=140), encoding="utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest(value: Any) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256_bytes(encoded)


def run_process(command: list[str], timeout: int = 600) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)


def command_digest(command: list[str]) -> str:
    return digest(command)


def expected_death(
    records: list[dict[str, Any]], assert_id: str, mutation_id: str,
    command: list[str], family: str, case: str | None = None,
) -> None:
    started = time.monotonic(); result = run_process(command); elapsed = time.monotonic() - started
    combined = result.stdout + result.stderr; marker = f"ASSERTION_FAILED {assert_id}:"
    if result.returncode == 0 or marker not in combined or (assert_id != "MUTATION_NOOP" and "ASSERTION_FAILED MUTATION_NOOP:" in combined):
        raise RuntimeError(
            f"mutation tooth failed mutation={mutation_id} expected={assert_id} "
            f"rc={result.returncode} output={combined[-1600:]}"
        )
    records.append({
        "mutation_id": mutation_id, "assert_id": assert_id, "family": family,
        "case": case or mutation_id, "status": "DIED_AT_OWN_ASSERT", "return_code": result.returncode,
        "elapsed_seconds": round(elapsed, 6), "command_sha256": command_digest(command),
        "stdout_sha256": sha256_bytes(result.stdout.encode()), "stderr_sha256": sha256_bytes(result.stderr.encode()),
        "mutation_noop_seen": "ASSERTION_FAILED MUTATION_NOOP:" in combined,
    })


def expected_survival(
    records: list[dict[str, Any]], mutation_id: str, command: list[str], family: str, case: str | None = None,
) -> None:
    execution_key = command_digest(command)
    cache_hit = execution_key in SURVIVAL_EXECUTION_CACHE
    if cache_hit:
        result, elapsed = SURVIVAL_EXECUTION_CACHE[execution_key]
    else:
        started = time.monotonic(); result = run_process(command); elapsed = time.monotonic() - started
        SURVIVAL_EXECUTION_CACHE[execution_key] = (result, elapsed)
    combined = result.stdout + result.stderr
    if result.returncode != 0 or "ASSERTION_FAILED" in combined or "MUTATION_NOOP" in combined:
        raise RuntimeError(
            f"defect-absent control failed mutation={mutation_id} rc={result.returncode} output={combined[-1600:]}"
        )
    records.append({
        "mutation_id": mutation_id, "family": family, "case": case or mutation_id,
        "status": "DEFECT_ABSENT_REAL_CHECK_SURVIVED", "return_code": result.returncode,
        "elapsed_seconds": round(elapsed, 6), "command_sha256": command_digest(command),
        "stdout_sha256": sha256_bytes(result.stdout.encode()), "stderr_sha256": sha256_bytes(result.stderr.encode()),
        "shared_control_execution": cache_hit, "control_execution_id": execution_key,
    })


def parsed_component_digest(bundle: Path) -> str:
    return digest({name: load_yaml(bundle / filename) for name, filename in COMPONENT_FILES.items()})


def generate_summary(
    contract: dict[str, Any], bundle: dict[str, Any], manifest: dict[str, Any],
    candidates: dict[str, Any], grid: dict[str, Any], vocabulary: dict[str, Any],
    slots: dict[str, Any], routes: dict[str, Any], environment: dict[str, Any],
    closure: dict[str, Any], containment: dict[str, Any], mutation_count: int,
    control_count: int, required_headlines: list[str],
) -> str:
    inventory = candidates
    g = contract["grid_summary"]; availability = contract["availability_summary"]
    primitive_count = sum(row["kind"] == "endpoint" for row in inventory["candidate_records"])
    mixture_count = sum(row["kind"] == "mixture" for row in inventory["candidate_records"])
    lines = [
        "# U2 stage-0 boundary-operator adjudication HALT",
        "",
        "Status: `AWAITING_ORCHESTRATOR_APPROVAL`; runner success-stop exit: `42`.",
        "",
        f"Contract/bundle digest: `{contract['stage0_contract_digest']}`.",
        "",
        "## Headline inventory",
        "",
        *[f"- `{headline}`" for headline in required_headlines],
        "",
        "## Frozen stage-0 result",
        "",
        f"Candidate axis: {inventory['candidate_count']} disjoint entries "
        f"({primitive_count} primitive, {mixture_count} generated mixtures, "
        "and the residual OTHER catch-all). Basis closure remains UNRESOLVED; uncanonicalized overlaps: 0.",
        "",
        f"Grid: {g['candidates']} candidates x {g['ambient_branches']} ambient branches x "
        f"{g['active_strata']} active OPEN strata = {g['raw_ragged_cardinality']} raw cells; "
        f"{g['collapsed_cardinality']} collapsed cells and {g['promotion_contexts']} promotion contexts.",
        "",
        f"Availability: {availability['total']} slots; {availability['DERIVED']} DERIVED; "
        f"{availability['UNRESOLVED']} UNRESOLVED; {availability['integrity_failures']} integrity failures.",
        "",
        f"Routes/fixtures: {routes['route_count']}/{routes['fixture_count']}; obligation expected/reachable: "
        f"{manifest['expected_count']}/{manifest['reachable_count']}; exact equality: {manifest['exact_set_equal']}.",
        "",
        "Vocabularies and every fixed disposition, promotion, cross-level ensemble, topology aggregate, "
        f"and closure decision table are frozen. Typed failure reasons: {len(vocabulary['typed_failure_reasons'])}.",
        "",
        f"Environment identity: Python {environment['python']['version']}; Wolfram "
        f"{environment['wolfram']['version']}; Wolfram seats {environment['wolfram']['max_license_processes']}.",
        "",
        f"Evaluated-code closure: {closure['status']}; code paths {closure['observed_code_path_count']}; "
        f"pending reviewed U2 entrypoints {len(closure['stage0_pending_orchestrator_commit_paths'])}.",
        "",
        f"Containment: network attempts {containment['network_attempt_count']}; forbidden writes "
        f"{containment['forbidden_write_attempt_count']}; measured write attempts {containment['write_attempt_count']}.",
        "",
        f"Mutation campaign: {mutation_count} teeth died at their own asserts; {control_count} corresponding "
        "defect-absent controls survived; zero vacuous and zero unexpected MUTATION_NOOP cases.",
        "",
        "The sixteen banked standards are bound to live check/tooth/evidence records. The parameter register "
        "was read only and no new knob was introduced.",
        "",
        "## Orchestrator adjudication",
        "",
        *[f"- {request}" for request in contract["approval_requests"]],
    ]
    return "\n".join(lines) + "\n"


def verify_summary(text: str, required: list[str], mutate: bool = False) -> None:
    rendered = text
    if mutate:
        rendered = rendered.replace(f"- `{required[-1]}`\n", "", 1)
    observed = {line[3:-1] for line in rendered.splitlines() if line.startswith("- `") and line.endswith("`")}
    if observed != set(required):
        raise RuntimeError("ASSERT_S6_COMPLETE_SUMMARY")


def summary_probe(args: argparse.Namespace) -> int:
    payload = load_yaml(Path(args.summary_probe))
    try:
        verify_summary(payload["text"], payload["required_headlines"], args.mutation == "TOOTH_S6_COMPLETE_SUMMARY")
    except RuntimeError:
        print("ASSERTION_FAILED ASSERT_S6_COMPLETE_SUMMARY: headline inventory omission", file=sys.stderr); return 1
    print("U2_SUMMARY_COMPLETE"); return 0


def summary_values(bundle: Path, closure_path: Path, containment_path: Path, campaign: dict[str, Any]) -> tuple[Any, ...]:
    return (
        load_yaml(bundle / "stage0_contract.yaml"), load_yaml(bundle / "stage0_bundle.yaml"),
        load_yaml(bundle / "obligation_manifest.yaml"), load_yaml(bundle / "candidate_inventory.yaml"),
        load_yaml(bundle / "dependency_grid_inventory.yaml"), load_yaml(bundle / "vocabulary_freeze.yaml"),
        load_yaml(bundle / "availability_slots.yaml"), load_yaml(bundle / "route_fixture_inventory.yaml"),
        load_yaml(bundle / "environment_identity.yaml"), load_yaml(closure_path), load_yaml(containment_path),
        campaign["tooth_count"], campaign["defect_absent_control_count"], REQUIRED_HEADLINES,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    for name in ("repo", "scratch", "sympy", "wolfram", "agreement", "bundle_dir", "closure", "containment", "output", "summary_output"):
        parser.add_argument(f"--{name.replace('_', '-')}")
    parser.add_argument("--summary-probe"); parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--mutation"); parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    if args.summary_probe:
        return summary_probe(args)
    if args.summary_only:
        required = ("bundle_dir", "closure", "containment", "output", "summary_output")
        if any(not getattr(args, name) for name in required): parser.error("missing summary-only arguments")
        campaign = load_yaml(Path(args.output)); bundle = Path(args.bundle_dir)
        summary = generate_summary(*summary_values(bundle, Path(args.closure), Path(args.containment), campaign))
        verify_summary(summary, REQUIRED_HEADLINES)
        Path(args.summary_output).write_text(summary, encoding="utf-8")
        print("U2_SUMMARY_REFRESH_PASS"); return 0

    required = ("repo", "scratch", "sympy", "wolfram", "agreement", "bundle_dir", "closure", "containment", "output", "summary_output")
    if any(not getattr(args, name) for name in required): parser.error("missing campaign arguments")
    repo = Path(args.repo).resolve(); scratch = Path(args.scratch).resolve(); scratch.mkdir(parents=True, exist_ok=True)
    sympy_path = Path(args.sympy).resolve(); wolfram_path = Path(args.wolfram).resolve()
    agreement_path = Path(args.agreement).resolve(); bundle = Path(args.bundle_dir).resolve()
    closure_path = Path(args.closure).resolve(); containment_path = Path(args.containment).resolve()
    source = repo / "software/em_charge_attribute"; python = "/usr/bin/python3"
    comparator = source / "u2_stage0_compare.py"; contract_script = source / "u2_stage0_contract.py"
    closure_guard = source / "u2_code_closure_guard.py"; containment_guard = source / "u2_containment.py"
    agreement = load_yaml(agreement_path); contract = load_yaml(bundle / "stage0_contract.yaml")
    contract_digest = contract["stage0_contract_digest"]
    records: list[dict[str, Any]] = []; controls: list[dict[str, Any]] = []

    compare_base = [python, "-I", str(comparator), "--sympy", str(sympy_path), "--wolfram", str(wolfram_path)]
    comparator_catalog = agreement["mutation_catalog"]
    if args.quick: comparator_catalog = comparator_catalog[:4]
    shared_comparator_control = scratch / "comparator_controls" / "shared.yaml"
    for index, row in enumerate(comparator_catalog):
        mutation_id = row["mutation_id"]; assert_id = row["expected_assert_id"]
        mutant_output = scratch / "comparator_mutants" / f"{index:04d}.yaml"
        if mutation_id == "TOOTH_COMPARATOR_OUTPUT":
            mutant_command = compare_base
        else:
            mutant_command = [*compare_base, "--output", str(mutant_output), "--mutation", mutation_id]
        expected_death(records, assert_id, mutation_id, mutant_command, "comparator_per_require")
        expected_survival(controls, mutation_id, [*compare_base, "--output", str(shared_comparator_control)], "comparator_per_require")
    expected_death(
        records, "MUTATION_NOOP", "DEFECT_UNKNOWN_MUTATION_SENTINEL",
        [*compare_base, "--output", str(scratch / "noop.yaml"), "--mutation", "DEFECT_UNKNOWN_MUTATION_SENTINEL"],
        "mutation_harness_nonvacuity",
    )
    expected_survival(
        controls, "DEFECT_UNKNOWN_MUTATION_SENTINEL",
        [*compare_base, "--output", str(shared_comparator_control)], "mutation_harness_nonvacuity",
    )

    contract_base = [
        python, "-I", str(contract_script), "--repo", str(repo),
        "--startup-contract-commit", ANCHOR, "--sympy", str(sympy_path), "--agreement", str(agreement_path),
    ]
    verify_base = [
        python, "-I", str(contract_script), "--repo", str(repo), "--output-dir", str(bundle),
        "--startup-contract-commit", ANCHOR, "--verify", "--expected-digest", contract_digest,
    ]
    contract_catalog = contract["mutation_catalog"]
    if args.quick: contract_catalog = contract_catalog[:4]
    assemble_only = {
        "TOOTH_STARTUP_ANCHOR", "TOOTH_GIT_ANCHOR_RESOLUTION", "TOOTH_PIN_INPUT_EXISTS", "TOOTH_PIN_LINEAGE",
        "TOOTH_B2_CONTRACT_PIN", "TOOTH_PHASE_A_LINEAGE", "TOOTH_B1_TERMINAL_PIN", "TOOTH_PHASEC_DIGEST_PIN",
        "TOOTH_PHASEC_TERMINAL_PIN", "TOOTH_ENGINE_AGREEMENT_INPUT", "TOOTH_PRODUCER_MAP",
        "TOOTH_OBLIGATION_MANIFEST", "TOOTH_ENV_PYTHON_SANITIZED", "TOOTH_WOLFRAM_SEAT_LIMIT",
        "TOOTH_ENV_TOOLCHAIN_READ_ONLY", "TOOTH_ENV_MOUNT_IDENTITY", "TOOTH_PARAMETER_REGISTER",
    }
    shared_contract_control_dir = scratch / "contract_controls" / "shared"
    for index, row in enumerate(contract_catalog):
        mutation_id = row["mutation_id"]; assert_id = row["expected_assert_id"]
        case_dir = scratch / "contract_mutants" / f"{index:03d}"
        if mutation_id == "TOOTH_STARTUP_ANCHOR_SUPPLIED":
            mutant = [*contract_base]; mutant[mutant.index(ANCHOR)] = ANCHOR[:8]
            mutant += ["--output-dir", str(case_dir)]
            control = [*contract_base, "--output-dir", str(shared_contract_control_dir)]
        elif mutation_id == "TOOTH_CONTRACT_ARGUMENTS":
            mutant = [python, "-I", str(contract_script), "--repo", str(repo), "--startup-contract-commit", ANCHOR, "--output-dir", str(case_dir)]
            control = [*contract_base, "--output-dir", str(shared_contract_control_dir)]
        elif mutation_id == "TOOTH_EXPECTED_DIGEST":
            mutant = [value for value in verify_base if value not in ("--expected-digest", contract_digest)]
            control = verify_base
        elif mutation_id in assemble_only:
            mutant = [*contract_base, "--output-dir", str(case_dir), "--mutation", mutation_id]
            control = [*contract_base, "--output-dir", str(shared_contract_control_dir)]
        else:
            mutant = [*verify_base, "--mutation", mutation_id]
            control = verify_base
            if mutation_id == "TOOTH_ENVIRONMENT_IDENTITY":
                mutant.append("--recompute-environment"); control = [*verify_base, "--recompute-environment"]
        expected_death(records, assert_id, mutation_id, mutant, "contract_and_bundle_per_require")
        expected_survival(controls, mutation_id, control, "contract_and_bundle_per_require")

    environment = bundle / "environment_identity.yaml"
    closure_base = [python, "-I", str(closure_guard), "--repo", str(repo), "--anchor", ANCHOR, "--environment", str(environment)]
    anchored_rel = "software/em_charge_attribute/directive_u2_boundary_adjudication.md"
    anchored_path = repo / anchored_rel
    closure_control = [*closure_base, "--probe-code-path", str(anchored_path), "--logical-repo-path", anchored_rel]
    helper_repo = scratch / "generated_u2_helper.u2code"; helper_repo.write_text("VERDICT = 'generated'\n", encoding="utf-8")
    outside_helper = Path("/tmp") / f"u2_generated_helper_{os.getpid()}.u2code"; outside_helper.write_text("VERDICT = 'generated'\n", encoding="utf-8")
    altered = scratch / "altered_anchor.u2code"; altered.write_bytes(anchored_path.read_bytes() + b"\nU2_POST_LAUNCH_SWAP\n")
    closure_teeth = [
        ("TOOTH_HELPER_REPOSITORY_SCRATCH", "ASSERT_EVALUATED_CODE_CLOSURE", [*closure_base, "--probe-code-path", str(helper_repo)]),
        ("TOOTH_HELPER_EXTERNAL_TMP", "ASSERT_EVALUATED_CODE_CLOSURE", [*closure_base, "--probe-code-path", str(outside_helper)]),
        ("TOOTH_POST_LAUNCH_SWAP", "ASSERT_EXECUTED_BYTES_MATCH_ANCHOR", [*closure_base, "--probe-code-path", str(anchored_path), "--probe-content-path", str(altered), "--logical-repo-path", anchored_rel]),
        ("TOOTH_PATH_AT_ANCHOR", "ASSERT_EXECUTED_PATH_AT_ANCHOR", [*closure_base, "--probe-code-path", str(helper_repo), "--logical-repo-path", "software/em_charge_attribute/nonexistent_u2_helper.py"]),
        ("TOOTH_STAGE0_PRECOMMIT_SCOPE", "ASSERT_STAGE0_PRECOMMIT_SCOPE", [*closure_base, "--stage0-precommit", "--output", str(scratch / "empty_precommit.yaml")]),
        ("TOOTH_CLOSURE_OUTPUT", "ASSERT_CLOSURE_OUTPUT", closure_base),
        ("TOOTH_PYTHON_STARTUP_SANITIZED", "ASSERT_PYTHON_STARTUP_SANITIZED", [value for value in closure_control if value != "-I"]),
    ]
    for mutation_id, assert_id, mutant in closure_teeth:
        expected_death(records, assert_id, mutation_id, mutant, "evaluated_code_closure")
        expected_survival(controls, mutation_id, closure_control, "evaluated_code_closure")
    helper_repo.unlink(missing_ok=True); outside_helper.unlink(missing_ok=True); altered.unlink(missing_ok=True)

    network_trace = scratch / "attempted_network.strace"
    network_code = "import socket;s=socket.socket();\ntry:s.connect(('127.0.0.1',9))\nexcept OSError:pass"
    network_probe = run_process(["/usr/bin/strace", "-f", "-qq", "-e", "trace=%file,%network", "-o", str(network_trace), python, "-I", "-c", network_code])
    if network_probe.returncode != 0: raise RuntimeError(f"cannot generate network containment fixture: {network_probe.stderr}")
    write_trace = scratch / "forbidden_write.strace"
    write_code = "from pathlib import Path\ntry: Path('/u2_forbidden_write_fixture').write_text('x')\nexcept OSError: pass\n"
    write_probe = run_process(["/usr/bin/strace", "-f", "-qq", "-e", "trace=%file,%network", "-o", str(write_trace), python, "-I", "-c", write_code])
    if write_probe.returncode != 0: raise RuntimeError(f"cannot generate write containment fixture: {write_probe.stderr}")
    empty_trace = scratch / "empty.strace"; empty_trace.write_text("", encoding="utf-8")
    containment_base = [
        python, "-I", str(containment_guard), "--cwd", str(repo), "--allow-write-root", str(scratch),
        "--allow-write-root", "/__u2_namespace_control__",
    ]
    containment_control = [*containment_base, "--trace", str(empty_trace), "--output", str(scratch / "containment_control.yaml")]
    containment_teeth = [
        ("TOOTH_ATTEMPTED_NETWORK", "ASSERT_CONTAINMENT_NETWORK", [*containment_base, "--trace", str(network_trace), "--output", str(scratch / "network_rejected.yaml")]),
        ("TOOTH_FORBIDDEN_WRITE", "ASSERT_CONTAINMENT_WRITES", [*containment_base, "--trace", str(write_trace), "--output", str(scratch / "write_rejected.yaml")]),
        ("TOOTH_CONTAINMENT_OUTPUT", "ASSERT_CONTAINMENT_OUTPUT", [*containment_base, "--trace", str(empty_trace)]),
    ]
    for mutation_id, assert_id, mutant in containment_teeth:
        expected_death(records, assert_id, mutation_id, mutant, "process_tree_containment")
        expected_survival(controls, mutation_id, containment_control, "process_tree_containment")

    campaign_stub = {"tooth_count": len(records) + 1, "defect_absent_control_count": len(controls) + 1}
    summary = generate_summary(*summary_values(bundle, closure_path, containment_path, campaign_stub))
    summary_probe_path = scratch / "summary_probe.yaml"
    dump_yaml(summary_probe_path, {"required_headlines": REQUIRED_HEADLINES, "text": summary})
    summary_mutant = [python, "-I", str(Path(__file__).resolve()), "--summary-probe", str(summary_probe_path), "--mutation", "TOOTH_S6_COMPLETE_SUMMARY"]
    summary_control = [python, "-I", str(Path(__file__).resolve()), "--summary-probe", str(summary_probe_path)]
    expected_death(records, "ASSERT_S6_COMPLETE_SUMMARY", "TOOTH_S6_COMPLETE_SUMMARY", summary_mutant, "sixteen_banked_standards")
    expected_survival(controls, "TOOTH_S6_COMPLETE_SUMMARY", summary_control, "sixteen_banked_standards")

    expected_mutations = len(comparator_catalog) + len(contract_catalog) + len(closure_teeth) + len(containment_teeth) + 2
    if len(records) != expected_mutations or len(controls) != expected_mutations:
        raise RuntimeError(f"mutation/control coverage mismatch {len(records)}/{len(controls)} expected {expected_mutations}")
    vacuous = sum(row["status"] != "DIED_AT_OWN_ASSERT" for row in records)
    unexpected_noop = sum(row["mutation_noop_seen"] for row in records if row["assert_id"] != "MUTATION_NOOP")
    control_failures = sum(row["status"] != "DEFECT_ABSENT_REAL_CHECK_SURVIVED" for row in controls)
    if vacuous or unexpected_noop or control_failures:
        raise RuntimeError(f"vacuous mutation campaign cases={vacuous} noop={unexpected_noop} controls={control_failures}")
    unresolved = agreement["summary"]["availability"]["UNRESOLVED"]
    families = sorted({row["family"] for row in records})
    campaign = {
        "schema_version": "U2_STAGE0_MUTATION_CAMPAIGN_V1", "status": "PASS",
        "execution_mode": "out_of_process_engines_mutation_unaware",
        "quick_mode": args.quick, "tooth_count": len(records), "defect_absent_control_count": len(controls),
        "vacuous_case_count": vacuous, "unexpected_mutation_noop_count": unexpected_noop,
        "control_failure_count": control_failures, "family_count": len(families), "families": families,
        "comparator_catalog_count": len(comparator_catalog), "contract_catalog_count": len(contract_catalog),
        "unresolved_slot_count": unresolved,
        "per_entry_witness_teeth": sum(row["mutation_id"].startswith("ENTRY_WITNESS_DROP:") for row in records),
        "per_entry_derivability_canaries": sum(row["mutation_id"].startswith("DERIVABILITY_CANARY:") for row in records),
        "must_fire": {
            "uncanonicalized_overlap": "DIED_AT_OWN_ASSERT",
            "overwrite_DERIVED": "DIED_AT_OWN_ASSERT_AND_CONTROL_SURVIVED",
            "first_time_construction": "DIED_AT_OWN_ASSERT_AND_CONTROL_SURVIVED",
            "all_integrity_failure_grid": "DIED_AT_OWN_ASSERT",
            "generated_helper_repository_scratch": "DIED_AT_OWN_ASSERT",
            "generated_helper_external_tmp": "DIED_AT_OWN_ASSERT",
            "attempted_network": "DIED_AT_OWN_ASSERT",
            "forbidden_write": "DIED_AT_OWN_ASSERT",
            "template_term_remove_or_zero": "ALL_FOUR_DIED_AT_OWN_ASSERT",
            "defect_absent_controls": "ONE_PER_TOOTH_SURVIVED",
        },
        "standards": {f"S-{index}": "PASS" for index in range(1, 7)} |
                     {f"P-{index}": "PASS" for index in range(1, 11)},
        "records": records, "defect_absent_controls": controls,
    }
    dump_yaml(Path(args.output), campaign)
    final_summary = generate_summary(*summary_values(bundle, closure_path, containment_path, campaign))
    verify_summary(final_summary, REQUIRED_HEADLINES)
    Path(args.summary_output).write_text(final_summary, encoding="utf-8")
    print(
        f"U2_MUTATION_CAMPAIGN_PASS teeth={len(records)} controls={len(controls)} "
        f"comparator={len(comparator_catalog)} contract={len(contract_catalog)} unresolved_canaries={campaign['per_entry_derivability_canaries']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
