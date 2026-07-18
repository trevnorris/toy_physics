#!/usr/bin/env python3
"""Sole authorized B2 aggregate merger with B1 containment checks."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

import yaml

from u1_body_b2_common import dump_yaml, load_yaml, load_yaml_authenticated, rel_repo, require


def authenticated_bytes(path: Path, expected: str, consumer: str) -> tuple[bytes, dict]:
    fd = os.open(path.resolve(), os.O_RDONLY)
    try:
        before = os.fstat(fd); chunks = []
        while True:
            chunk = os.read(fd, 1 << 20)
            if not chunk: break
            chunks.append(chunk)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    raw = b"".join(chunks); observed = hashlib.sha256(raw).hexdigest()
    require(observed == expected, "B2_A1_PROTECTED_FIRST_USE", f"{consumer}:{path}")
    require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns), "B2_A1_PROTECTED_REWRITE", f"{consumer}:{path}")
    return raw, {"consumer": consumer, "path": rel_repo(path), "expected_sha256": expected, "consumed_sha256": observed, "held_descriptor": True}


def stage_and_write_in_place(target: Path, staging: Path, serialized: bytes) -> str:
    """Verify staged bytes, then overwrite an existing target inode in place."""
    staging.parent.mkdir(parents=True, exist_ok=True)
    with staging.open("wb") as stream:
        require(stream.write(serialized) == len(serialized), "B2_MERGE_STAGING", f"short staging write:{staging}")
        stream.flush()
        os.fsync(stream.fileno())
    staged_bytes = staging.read_bytes()
    require(staged_bytes == serialized, "B2_MERGE_STAGING", f"staged byte mismatch:{staging}")
    staged_sha256 = hashlib.sha256(staged_bytes).hexdigest()

    resolved_target = target.resolve(strict=True)
    require(resolved_target.is_file(), "B2_MERGE_TARGET", f"existing file required:{resolved_target}")
    with resolved_target.open("wb") as stream:
        require(stream.write(staged_bytes) == len(staged_bytes), "B2_MERGE_TARGET", f"short final write:{resolved_target}")
        stream.flush()
        os.fsync(stream.fileno())
    final_bytes = resolved_target.read_bytes()
    final_sha256 = hashlib.sha256(final_bytes).hexdigest()
    require(final_bytes == staged_bytes and final_sha256 == staged_sha256, "B2_MERGE_TARGET", f"final byte mismatch:{resolved_target}")
    return final_sha256


def main() -> int:
    p = argparse.ArgumentParser(); p.add_argument("--agreement", type=Path, required=True); p.add_argument("--stage0", type=Path, required=True); p.add_argument("--stage0-contract-digest", required=True); p.add_argument("--b1-results-snapshot", type=Path, required=True); p.add_argument("--b1-report-snapshot", type=Path, required=True); p.add_argument("--results", type=Path, required=True); p.add_argument("--report", type=Path, required=True); p.add_argument("--provenance", type=Path, required=True); a = p.parse_args()
    try:
        agreement = load_yaml(a.agreement)
        stage0, stage_auth = load_yaml_authenticated(a.stage0, a.stage0_contract_digest, "merger:stage0")
        manifest = stage0["observation_contract"]["Obs_B2_manifest"]
        mutable_aggregate = stage0["observation_contract"]["mutable_aggregate"]
        expected = lambda path: manifest[rel_repo(path.resolve())]

        def mutable_expected(path: Path) -> str:
            matches = [raw for raw in mutable_aggregate if Path(raw).name == path.name]
            require(len(matches) == 1, "B2_MERGE_TARGET", f"contract mutable aggregate:{path.name}")
            return manifest[matches[0]]

        snapshot, snapshot_auth = load_yaml_authenticated(a.b1_results_snapshot, expected(a.b1_results_snapshot), "merger:B1_results_snapshot")
        require(agreement["status"] == "PASS_WITH_HONEST_OUTCOMES", "B2_MERGE_AGREEMENT", "comparison status")
        current, current_auth = load_yaml_authenticated(a.results, mutable_expected(a.results), "merger:mutable_results_before_write")
        if "phases" in current:
            require(current["phases"]["B1"] == snapshot, "B2_MERGE_B1_CONTAINMENT", "resumed nested B1")
        else:
            require(current == snapshot, "B2_MERGE_B1_CONTAINMENT", "flat startup B1")
        frozen_md, frozen_md_auth = authenticated_bytes(a.b1_report_snapshot, expected(a.b1_report_snapshot), "merger:B1_report_snapshot")
        current_md, current_md_auth = authenticated_bytes(a.report, mutable_expected(a.report), "merger:mutable_report_before_write")
        require(current_md == frozen_md or current_md.startswith(frozen_md + b"\n"), "B2_MERGE_B1_MARKDOWN", "frozen prefix")
        merged = {"schema_version": "U1_BODY_DYNAMICS_AGGREGATE_V2", "phase": "B2", "current_headline": agreement["headline"], "phases": {"B1": snapshot, "B2": agreement}}
        result_bytes = yaml.safe_dump(merged, sort_keys=False, allow_unicode=True, width=220).encode("utf-8")
        result_sha256 = stage_and_write_in_place(a.results, a.provenance.with_name(f".{a.provenance.name}.results.staged"), result_bytes)
        section = "\n## Phase B2 — Intake response and radiative residues\n\n" + "- Status: `PASS_WITH_HONEST_OUTCOMES`\n" + "- Partition: `UNRESOLVED(return_closure)`\n- C_mdot: `UNRESOLVED(required_OPEN_leaves)`\n- Radiation: `UNRESOLVED(native_branch_inputs)`\n- Phase-C gates: `NOT_RUN(phase_C)`\n"
        report_sha256 = stage_and_write_in_place(a.report, a.provenance.with_name(f".{a.provenance.name}.report.staged"), frozen_md + section.encode())
        provenance = {"schema_version": "U1_PHASE_B2_MERGE_PROVENANCE_V2", "status": "PASS", "sole_authorized_writer": Path(__file__).name, "first_use_authentication": [stage_auth, snapshot_auth, current_auth, frozen_md_auth, current_md_auth], "B1_results_semantic_containment": True, "B1_markdown_prefix_byte_identical": True, "result_sha256": result_sha256, "report_sha256": report_sha256}; dump_yaml(a.provenance, provenance)
        print("B2_MERGE: PASS phases.B1/phases.B2 migration"); return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__": raise SystemExit(main())
