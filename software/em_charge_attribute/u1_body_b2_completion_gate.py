#!/usr/bin/env python3
"""Traced-mode Obs_B2 closure, environment, and producer-chain completion gate."""

from __future__ import annotations

import argparse
from pathlib import Path

from u1_body_b2_common import ROOT, dump_yaml, load_yaml, load_yaml_authenticated, require, sha256_file


def path_set(rows: list[dict]) -> set[str]:
    return {row["path"] for row in rows}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--stage0", type=Path, required=True); p.add_argument("--stage0-contract-digest", required=True)
    p.add_argument("--environment-identity-digest", required=True); p.add_argument("--trace-closure", type=Path, required=True)
    p.add_argument("--merger-provenance", type=Path, required=True); p.add_argument("--authentication-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True); a = p.parse_args()
    try:
        contract, stage_auth = load_yaml_authenticated(a.stage0, a.stage0_contract_digest, "completion_gate:stage0")
        require(contract["trust_mode_predicate"]["environment_identity_digest"] == a.environment_identity_digest, "B2_A1_ENVIRONMENT_PIN", "completion pin")
        trace, merger = load_yaml(a.trace_closure), load_yaml(a.merger_provenance)
        expected = set(contract["observation_contract"]["Obs_B2_manifest"])
        container = contract["observation_contract"]["stage0_container_exclusion"]
        generated = path_set(trace["generated_output"])
        observed = path_set(trace["observed_reads"])
        actual = observed - generated - {container}
        require(actual == expected, "B2_A1_OBS_B2_EXACT", f"missing={sorted(expected-actual)} extra={sorted(actual-expected)}")
        first_use_rows = {row["path"]: row for row in trace["observed_preexisting_reads"]}
        require(set(first_use_rows) - {container} == expected, "B2_A1_OBS_B2_FIRST_USE", "production preexisting first-use set")
        for raw, expected_sha in contract["observation_contract"]["Obs_B2_manifest"].items():
            require(first_use_rows[raw]["sha256"] == expected_sha and first_use_rows[raw]["first_use"]["binding_semantics"].startswith("regular-file bytes SHA-256 before"), "B2_A1_OBS_B2_FIRST_USE", raw)
        protected = expected - set(contract["observation_contract"]["mutable_aggregate"])
        writes = path_set(trace["observed_writes"])
        require(not (protected & writes), "B2_A1_PROTECTED_NO_WRITE", f"writes={sorted(protected&writes)}")

        pinned_environment = {row["path"]: row for row in contract["trust_mode_predicate"]["environment_closure"]["entries"]}
        production_environment = {row["path"]: row for row in trace["environment_closure"]["entries"]}
        require(set(production_environment) <= set(pinned_environment), "B2_A1_ENVIRONMENT_PIN", f"unapproved production environment paths={sorted(set(production_environment)-set(pinned_environment))}")
        for path, row in production_environment.items():
            require(row["sha256"] == pinned_environment[path]["sha256"], "B2_A1_ENVIRONMENT_PIN", f"production first-use digest:{path}")

        env_rows = [load_yaml(path) for path in sorted((a.authentication_root / "environment_assertions").glob("*.yaml"))]
        require(bool(env_rows) and all(row.get("status") == "PASS" and row.get("environment_identity_digest") == a.environment_identity_digest and row.get("live_equals_pinned_record") for row in env_rows), "B2_A1_ENVIRONMENT_PIN", "all stage assertions")
        env_consumers = {row["consumer"] for row in env_rows}
        first_rows = [load_yaml(path) for path in sorted((a.authentication_root / "first_use_authentication").glob("*.yaml"))]
        require(bool(first_rows) and all(row.get("status") == "PASS" and row.get("protected_path_rewrite_check") == "PASS" for row in first_rows), "B2_A1_PROTECTED_FIRST_USE", "authenticated stage executables")
        first_consumers = {row["consumer"] for row in first_rows}
        require(first_consumers <= env_consumers and {"resume_entry", "completion_gate", "production_trace_audit"} <= env_consumers, "B2_A1_ENVIRONMENT_PIN", f"missing assertions for {sorted(first_consumers-env_consumers)}")

        producer_records = []
        for path in sorted(a.authentication_root.glob("*producer*digest*.yaml")):
            record = load_yaml(path); require(record.get("status") == "PASS", "B2_A1_PRODUCER_DIGEST", path.name)
            for output in record["outputs"]:
                require(sha256_file(Path(output["path"])) == output["sha256"], "B2_A1_PRODUCER_DIGEST", output["path"])
            producer_records.append({"path": path.name, "producer": record["producer"], "output_count": len(record["outputs"])})
        require(len(producer_records) >= 3, "B2_A1_PRODUCER_DIGEST", "phase-A, targeted, and production producer records")

        mutation_path = Path(contract["dual_engine_evidence"]["mutations"]["path"])
        if not mutation_path.is_absolute():
            mutation_path = Path(__file__).resolve().parents[2] / mutation_path
        require(sha256_file(mutation_path) == contract["dual_engine_evidence"]["mutations"]["sha256"], "B2_A1_PROTECTED_REWRITE", "stage0 mutation evidence digest")
        mutations = load_yaml(mutation_path)
        rewrite = [row for row in mutations["cases"] if row["id"] == "trace_semantics::relative_dirfd_protected_rewrite"]
        require(len(rewrite) == 1 and rewrite[0]["killed_at_own_assert"] and rewrite[0]["expected_assert"] == "B2_A1_PROTECTED_REWRITE", "B2_A1_PROTECTED_REWRITE", "must-fire mid-run path swap")
        require(merger["status"] == "PASS" and merger["B1_results_semantic_containment"] and merger["B1_markdown_prefix_byte_identical"], "B2_A1_MERGER", "aggregate containment")
        dump_yaml(a.output, {
            "schema_version": "U1_PHASE_B2_COMPLETION_GATE_V3", "status": "PASS", "stage0_first_use_authentication": stage_auth,
            "environment_identity_digest": a.environment_identity_digest, "environment_assertion_consumers": sorted(env_consumers),
            "first_use_consumers": sorted(first_consumers), "producer_records": producer_records,
            "protected_rewrite_tooth": rewrite[0], "Obs_B2_exact": True, "Obs_B2_equation": "production Obs_B2 - production generated_output - stage0_container == frozen manifest", "historical_closure_union_used": False, "manifest_count": len(expected),
            "protected_write_intersection": [], "merger": merger,
        })
        print("B2_COMPLETION_GATE: PASS"); return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__":
    raise SystemExit(main())
