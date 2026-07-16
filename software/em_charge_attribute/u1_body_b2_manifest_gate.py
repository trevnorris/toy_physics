#!/usr/bin/env python3
"""First-use authentication of the Stage-0 container and every Obs_B2 entry."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

import yaml

from u1_body_b2_common import ROOT, UniqueKeyLoader, dump_yaml, require


def consume(path: Path) -> tuple[bytes, dict]:
    fd = os.open(path, os.O_RDONLY)
    try:
        chunks = []
        while True:
            block = os.read(fd, 1 << 20)
            if not block:
                break
            chunks.append(block)
        blob = b"".join(chunks)
        stat = os.fstat(fd)
        return blob, {"device": stat.st_dev, "inode": stat.st_ino, "consumed_fd_sha256": hashlib.sha256(blob).hexdigest()}
    finally:
        os.close(fd)


def main() -> int:
    p = argparse.ArgumentParser(); p.add_argument("--stage0", type=Path, required=True); p.add_argument("--stage0-contract-digest", required=True); p.add_argument("--output", type=Path, required=True); a = p.parse_args()
    try:
        stage_blob, stage_evidence = consume(a.stage0); actual = stage_evidence["consumed_fd_sha256"]; require(actual == a.stage0_contract_digest, "B2_A1_STAGE0_CONTAINER", "orchestrator-supplied digest")
        row = yaml.load(stage_blob.decode("utf-8"), Loader=UniqueKeyLoader); entries = row["observation_contract"]["Obs_B2_manifest"]; verified = []
        for raw, expected in entries.items():
            path = ROOT / raw; require(path.is_file(), "B2_A1_OBS_B2_FIRST_USE", f"missing:{raw}"); _, evidence = consume(path); computed = evidence["consumed_fd_sha256"]; require(computed == expected, "B2_A1_OBS_B2_FIRST_USE", f"digest:{raw}"); verified.append({"path": raw, "expected_sha256": expected, **evidence, "consumer": "manifest_gate_preflight"})
        dump_yaml(a.output, {"schema_version": "U1_PHASE_B2_MANIFEST_GATE_V2", "status": "PASS", "stage0_contract": stage_evidence, "stage0_contract_sha256": actual, "entry_count": len(verified), "entries": verified, "per_consumer_reauthentication_required": True})
        print(f"B2_MANIFEST_GATE: PASS entries={len(verified)}"); return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__": raise SystemExit(main())
