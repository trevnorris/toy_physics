#!/usr/bin/env python3
"""Sole authorized B2 aggregate merger with B1 containment checks."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

from u1_body_b2_common import dump_yaml, load_yaml, load_yaml_authenticated, rel_repo, require, sha256_file


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


def main() -> int:
    p = argparse.ArgumentParser(); p.add_argument("--agreement", type=Path, required=True); p.add_argument("--stage0", type=Path, required=True); p.add_argument("--stage0-contract-digest", required=True); p.add_argument("--b1-results-snapshot", type=Path, required=True); p.add_argument("--b1-report-snapshot", type=Path, required=True); p.add_argument("--results", type=Path, required=True); p.add_argument("--report", type=Path, required=True); p.add_argument("--provenance", type=Path, required=True); a = p.parse_args()
    try:
        agreement = load_yaml(a.agreement)
        stage0, stage_auth = load_yaml_authenticated(a.stage0, a.stage0_contract_digest, "merger:stage0")
        manifest = stage0["observation_contract"]["Obs_B2_manifest"]
        expected = lambda path: manifest[rel_repo(path.resolve())]
        snapshot, snapshot_auth = load_yaml_authenticated(a.b1_results_snapshot, expected(a.b1_results_snapshot), "merger:B1_results_snapshot")
        require(agreement["status"] == "PASS_WITH_HONEST_OUTCOMES", "B2_MERGE_AGREEMENT", "comparison status")
        current, current_auth = load_yaml_authenticated(a.results, expected(a.results), "merger:mutable_results_before_write")
        if "phases" in current:
            require(current["phases"]["B1"] == snapshot, "B2_MERGE_B1_CONTAINMENT", "resumed nested B1")
        else:
            require(current == snapshot, "B2_MERGE_B1_CONTAINMENT", "flat startup B1")
        frozen_md, frozen_md_auth = authenticated_bytes(a.b1_report_snapshot, expected(a.b1_report_snapshot), "merger:B1_report_snapshot")
        current_md, current_md_auth = authenticated_bytes(a.report, expected(a.report), "merger:mutable_report_before_write")
        require(current_md == frozen_md or current_md.startswith(frozen_md + b"\n"), "B2_MERGE_B1_MARKDOWN", "frozen prefix")
        merged = {"schema_version": "U1_BODY_DYNAMICS_AGGREGATE_V2", "phase": "B2", "current_headline": agreement["headline"], "phases": {"B1": snapshot, "B2": agreement}}
        temp_results = a.results.with_suffix(a.results.suffix + ".b2tmp"); dump_yaml(temp_results, merged); os.replace(temp_results, a.results)
        section = "\n## Phase B2 — Intake response and radiative residues\n\n" + "- Status: `PASS_WITH_HONEST_OUTCOMES`\n" + "- Partition: `UNRESOLVED(return_closure)`\n- C_mdot: `UNRESOLVED(required_OPEN_leaves)`\n- Radiation: `UNRESOLVED(native_branch_inputs)`\n- Phase-C gates: `NOT_RUN(phase_C)`\n"
        temp_report = a.report.with_suffix(a.report.suffix + ".b2tmp"); temp_report.write_bytes(frozen_md + section.encode()); os.replace(temp_report, a.report)
        provenance = {"schema_version": "U1_PHASE_B2_MERGE_PROVENANCE_V2", "status": "PASS", "sole_authorized_writer": Path(__file__).name, "first_use_authentication": [stage_auth, snapshot_auth, current_auth, frozen_md_auth, current_md_auth], "B1_results_semantic_containment": True, "B1_markdown_prefix_byte_identical": True, "result_sha256": sha256_file(a.results), "report_sha256": sha256_file(a.report)}; dump_yaml(a.provenance, provenance)
        print("B2_MERGE: PASS phases.B1/phases.B2 migration"); return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__": raise SystemExit(main())
