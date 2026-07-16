#!/usr/bin/env python3
"""Assemble the immutable, non-self-hashing Stage-0 B2 contract container."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from u1_body_b2_common import HERE, ROOT, digest, dump_yaml, load_yaml, manifest_absolute, parse_manifest_sets, rel_repo, require, sha256_file


def checked(path: Path, expected_status: str | tuple[str, ...] = "PASS") -> dict[str, Any]:
    value = load_yaml(path)
    statuses = (expected_status,) if isinstance(expected_status, str) else expected_status
    require(value.get("status") in statuses, "B2_S0_EVIDENCE", f"{path.name}: {value.get('status')} not in {statuses}")
    return value


def manifest(certificate: dict[str, Any], stage_dir: Path, container: Path) -> dict[str, str]:
    paths: set[Path] = set()
    for rows in parse_manifest_sets(certificate).values():
        paths.update(manifest_absolute(raw) for raw in rows)
    paths.update([HERE / "reports/u1_body_dynamics_results.yaml", HERE / "reports/u1_body_dynamics.md"])
    paths.update(path for path in HERE.glob("u1_body_b2*") if path.is_file())
    paths.update(path for path in HERE.glob("run_u1_body_b2*") if path.is_file())
    paths.add(HERE / "directive_u1_phaseB2_intake_radiative.md")
    # The approval certificate is consumed directly by every startup/recheck;
    # it is an Obs_B2 dependency even though its own three member manifests do
    # not recursively list their container.
    paths.add(HERE / "reports/u1_body_dynamics_artifacts/b1_orchestrator_approval.yaml")
    paths.update(path for path in stage_dir.iterdir() if path.is_file() and path != container)
    require(container not in paths, "B2_S0_MANIFEST_SELF", "container excluded from its own manifest")
    missing = [path for path in paths if not path.is_file()]
    require(not missing, "B2_S0_MANIFEST", f"missing paths: {missing}")
    return {rel_repo(path): sha256_file(path) for path in sorted(paths)}


def merge_environment_closures(*closures: dict[str, Any]) -> dict[str, Any]:
    merged: dict[str, dict[str, Any]] = {}
    for closure in closures:
        for row in closure["entries"]:
            path = row["path"]
            current = merged.get(path)
            if current is None:
                merged[path] = {"path": path, "sha256": row["sha256"], "categories": sorted(row["categories"]), "first_use_sources": [row["audit_log"]]}
            else:
                require(current["sha256"] == row["sha256"], "B2_S0_ENVIRONMENT_CLOSURE", f"bytes changed across first-use closures:{path}")
                current["categories"] = sorted(set(current["categories"]) | set(row["categories"]))
                current["first_use_sources"] = sorted(set(current["first_use_sources"]) | {row["audit_log"]})
    entries = sorted(merged.values(), key=lambda row: row["path"])
    canonical = [{"path": row["path"], "sha256": row["sha256"], "categories": row["categories"]} for row in entries]
    result = {
        "schema_version": "U1_PHASE_B2_ENVIRONMENT_CLOSURE_V1",
        "status": "PROPOSED_FOR_ORCHESTRATOR_PIN",
        "first_use_mechanism": "runner hashes open/openat/openat2 descriptors before return, pre-reads Wolfram regular-file arguments before kernel dispatch, and records executable/loader objects before application main or dlopen return",
        "entry_count": len(entries),
        "entries": entries,
        "category_counts": {category: sum(category in row["categories"] for row in entries) for category in ["executed_binary", "loaded_shared_library", "loader", "imported_python_module", "wolfram_installation_file"]},
    }
    result["environment_closure_digest"] = digest({"entries": canonical})
    require(all(value > 0 for value in result["category_counts"].values()), "B2_S0_ENVIRONMENT_CLOSURE", f"missing category:{result['category_counts']}")
    return result


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True); p.add_argument("--startup", type=Path, required=True)
    p.add_argument("--sympy", type=Path, required=True); p.add_argument("--mathematica", type=Path, required=True); p.add_argument("--agreement", type=Path, required=True)
    p.add_argument("--mutations", type=Path, required=True); p.add_argument("--b1-full-closure", type=Path, required=True)
    p.add_argument("--schema-integration", type=Path, required=True)
    p.add_argument("--b1-targeted-replay", type=Path, required=True); p.add_argument("--b1-targeted-closure", type=Path, required=True)
    p.add_argument("--stage0-trace-closure", type=Path, required=True); p.add_argument("--output", type=Path, required=True)
    p.add_argument("--container-path", type=Path, help="logical final container path when assembling an authenticated re-emit candidate")
    a = p.parse_args()
    try:
        config = load_yaml(a.input); startup = checked(a.startup); sym = checked(a.sympy, "ENGINE_DERIVED"); math = checked(a.mathematica, "ENGINE_DERIVED"); agreement = checked(a.agreement, "ENGINE_AGREE")
        mutations = checked(a.mutations); schema_integration = checked(a.schema_integration); full = checked(a.b1_full_closure); targeted = checked(a.b1_targeted_replay); targeted_closure = checked(a.b1_targeted_closure); trace = checked(a.stage0_trace_closure)
        require(startup["trust_mode"] == "traced", "B2_S0_TRUST_MODE", "traced branch")
        require(full["scope"] == "full_b1" and targeted_closure["scope"] == "targeted_b1" and trace["scope"] == "stage0", "B2_S0_CLOSURES", "distinct closure universes")
        certificate = load_yaml(ROOT / config["contracts"]["certificate"])
        logical_container = (a.container_path or a.output).resolve()
        obs = manifest(certificate, logical_container.parent, logical_container)
        frozen = agreement["frozen_math"]
        environment_closure = merge_environment_closures(full["environment_closure"], targeted_closure["environment_closure"], trace["environment_closure"])
        require(environment_closure["entry_count"] > 0 and len(environment_closure["environment_closure_digest"]) == 64, "B2_S0_ENVIRONMENT_CLOSURE", "mechanically merged v48 first-use closure")
        require(sym["schema_version"] == math["schema_version"] == "U1_PHASE_B2_STAGE0_ENGINE_V6", "B2_S0_EVIDENCE", "v48 engine schemas")
        require(agreement["schema_version"] == "U1_PHASE_B2_STAGE0_MATH_AGREEMENT_V6" and mutations["schema_version"] == "U1_PHASE_B2_STAGE0_MUTATIONS_V6", "B2_S0_EVIDENCE", "v48 agreement/mutation schemas")
        require(frozen["stage0_datum_bank"]["record_count"] == len(frozen["stage0_datum_dispositions"]) and frozen["stage0_datum_bank"]["expected_reachable_exact_set_equal"] is True, "B2_S0_DISPOSITION_FLOOR", "banked v48 exact datum floor")
        artifact = {
            "schema_version": "U1_PHASE_B2_STAGE0_CONTRACT_V6",
            "status": "AWAITING_ORCHESTRATOR_APPROVAL",
            "required_exit_code": 42,
            "phase": "B2_STAGE0",
            "orchestrator_inputs": {"startup_contract_commit": startup["startup_contract_commit"], "trust_mode": startup["trust_mode"], "stage0_contract_digest": "MUST_BE_SUPPLIED_ON_RESUME", "environment_closure_digest": "MUST_BE_SUPPLIED_AS_ENVIRONMENT_IDENTITY_DIGEST_ON_RESUME_AND_STAGE0_REEMIT"},
            "contract_anchor": {"directive_version": config["directive_version"], "startup_contract_commit": startup["startup_contract_commit"], "required_v48_commit": config["startup_contract_commit"], "spec_version": config["spec_version"], "contract_digests": startup["contract_digests"]},
            "trust_mode_predicate": {
                "selected": "traced", "isolation_clauses_used": False,
                "environment_closure": environment_closure,
                "environment_closure_digest": environment_closure["environment_closure_digest"],
                "environment_identity_digest": environment_closure["environment_closure_digest"],
                "environment_closure_pin_state": "AWAITING_ORCHESTRATOR_PIN_AT_APPROVAL",
                "environment_assertion_rule": "every stage re-hashes every pinned closure entry; every production dependency is additionally bound at its first-use audit hook",
                "read_only_root_sandbox": startup["read_only_root_sandbox"],
                "sandbox_rule": "every production script enters bubblewrap with --ro-bind / /; only enumerated generated-output/temp mounts and the two merger-owned aggregate files are writable",
                "no_network": config["no_network"],
                "policy_result": "PASS: every production probe ran in a bubblewrap network namespace; traced syscall audit found no external network operation",
            },
            "frozen_data": {
                "causal_definition": frozen["causal_definition"],
                "fourier_convention": frozen["fourier_convention"],
                "complete_action_second_variation": frozen["complete_action_second_variation"],
                "sector_current_derivations": frozen["current_derivations"],
                "native_noether_derivations": frozen["native_noether_derivations"],
                "integrated_balance_identities": frozen["integrated_balance_identities"],
                "native_operator_inventory": frozen["operator_inventory"],
                "native_operator_action_incidence": frozen["operator_action_incidence"],
                "native_operator_action_incidence_residual": frozen["operator_action_incidence_residual"],
                "endpoint_resolvent_cells": frozen["endpoint_resolvent_cells"],
                "functional_analytic_test_data": frozen["functional_analytic_test_data"],
                "green_sensitivity_matrix": frozen["green_sensitivity_matrix"],
                "trajectory_representation": frozen["trajectory_representation"],
                "validity_domain_construction": frozen["validity_domain_construction"],
                "restriction_map": frozen["restriction_map"],
                "ownership_convention": frozen["ownership_convention"],
                "stage0_datum_bank": frozen["stage0_datum_bank"],
                "stage0_datum_dispositions": frozen["stage0_datum_dispositions"],
                "stage0_datum_disposition_sha256": frozen["stage0_datum_disposition_sha256"],
                "shared_input_whitelist": frozen["shared_input_whitelist"],
                "route_separation": frozen["route_separation"],
                "minimum_obligation_manifest": frozen["minimum_obligation_floor"],
            },
            "observation_contract": {
                "recipe": config["observation_recipe"],
                "set_definitions": startup["set_definitions"],
                "Obs_B2_manifest": obs,
                "Obs_B2_manifest_count": len(obs),
                "Obs_B2_manifest_set_sha256": digest(obs),
                "approval_request_target": {"reviewed_manifest_count": len(obs), "reviewed_manifest_digest": digest(obs)},
                "stage0_container_exclusion": rel_repo(logical_container),
                "mutable_aggregate": [config["contracts"]["mutable_results"], config["contracts"]["mutable_report"]],
                "generated_output_rule": "initially absent path whose first filesystem event is CREATE within an approved output prefix; pure-write inventory is path-only and every later read is content-bound in observed_reads",
                "W_pure_rule": "actually opened for writing in approved prefixes, excluding every read-before-first-create/truncate path",
            },
            "B1_replay_evidence": {
                "full_traced_once": {"artifact": rel_repo(a.b1_full_closure), "sha256": sha256_file(a.b1_full_closure), "set_sha256": full["set_sha256"]},
                "targeted_semantic": {"artifact": rel_repo(a.b1_targeted_replay), "sha256": sha256_file(a.b1_targeted_replay), "checks": targeted["checks"]},
                "targeted_trace": {"artifact": rel_repo(a.b1_targeted_closure), "sha256": sha256_file(a.b1_targeted_closure), "set_sha256": targeted_closure["set_sha256"]},
                "normalizations": targeted["normalizations"], "excluded_volatile_paths": targeted["excluded_volatile_paths"],
            },
            "dual_engine_evidence": {"sympy": {"path": rel_repo(a.sympy), "sha256": sha256_file(a.sympy), "route": sym["independent_representation"]}, "mathematica": {"path": rel_repo(a.mathematica), "sha256": sha256_file(a.mathematica), "route": math["independent_representation"]}, "agreement": {"path": rel_repo(a.agreement), "sha256": sha256_file(a.agreement)}, "schema_integration": {"path": rel_repo(a.schema_integration), "sha256": sha256_file(a.schema_integration), "checks": schema_integration["checks"]}, "mutations": {"path": rel_repo(a.mutations), "sha256": sha256_file(a.mutations), "case_count": mutations["case_count"]}},
            "stage0_trace_evidence": {"path": rel_repo(a.stage0_trace_closure), "sha256": sha256_file(a.stage0_trace_closure), "external_network_event_count": trace["network"]["external_event_count"]},
            "production_route": {
                "stage1": "independent SymPy and Wolfram construction over the full generated key/operator inventory with named honest outcomes",
                "stage2": "engine-neutral semantic comparison, two-sided obligation inventory, out-of-process per-tooth mutations, typed liveness and no-double-count gates",
                "stage3": "traced completion seals plus sole-writer containment merger into phases.B1/phases.B2 aggregates",
                "amendment_rule": "any finding that would change a frozen Phase-A/B1 owner or value halts for orchestrator approval",
            },
            "stage0_semantic_sha256": digest({"frozen_data": frozen, "Obs_B2_manifest": obs, "startup_contract_commit": startup["startup_contract_commit"], "trust_mode": startup["trust_mode"], "environment_closure_digest": environment_closure["environment_closure_digest"]}),
        }
        dump_yaml(a.output, artifact)
        print(f"B2_STAGE0_CONTRACT: READY entries={len(obs)} exit_required=42")
        return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__": raise SystemExit(main())
