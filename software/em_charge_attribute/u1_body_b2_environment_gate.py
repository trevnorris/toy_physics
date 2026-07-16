#!/usr/bin/env python3
"""Assert the complete orchestrator-pinned executable environment closure."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from u1_body_b2_common import digest, dump_yaml, load_yaml_authenticated, require, sha256_file


def canonical(entries: list[dict]) -> list[dict]:
    return [{"path": row["path"], "sha256": row["sha256"], "categories": sorted(row["categories"])} for row in sorted(entries, key=lambda item: item["path"])]


def main() -> int:
    p = argparse.ArgumentParser(); p.add_argument("--stage0", type=Path, required=True); p.add_argument("--stage0-contract-digest", required=True); p.add_argument("--environment-identity-digest", required=True); p.add_argument("--no-network-shim", type=Path, required=True); p.add_argument("--consumer", required=True); p.add_argument("--output", type=Path, required=True); a = p.parse_args()
    try:
        stage0, auth = load_yaml_authenticated(a.stage0, a.stage0_contract_digest, f"environment_gate:{a.consumer}")
        closure = stage0["trust_mode_predicate"]["environment_closure"]
        pinned = stage0["trust_mode_predicate"]["environment_closure_digest"]
        require(pinned == a.environment_identity_digest == closure["environment_closure_digest"] == digest({"entries": canonical(closure["entries"])}), "B2_A1_ENVIRONMENT_PIN", f"{a.consumer}:orchestrator closure pin")
        require(closure["entry_count"] == len(closure["entries"]) > 0, "B2_A1_ENVIRONMENT_PIN", f"{a.consumer}:reviewed complete closure count")
        require(os.environ.get("B2_READ_ONLY_ROOT_SANDBOX") == "1", "B2_A1_READ_ONLY_ROOT", f"{a.consumer}:sandbox marker")
        root_mounts = [line.split() for line in Path("/proc/self/mountinfo").read_text(encoding="utf-8").splitlines() if len(line.split()) > 6 and line.split()[4] == "/"]
        require(root_mounts and all("ro" in row[5].split(",") for row in root_mounts), "B2_A1_READ_ONLY_ROOT", f"{a.consumer}:root mount")
        verified = []
        for row in closure["entries"]:
            path = Path(row["path"])
            require(path.is_file(), "B2_A1_ENVIRONMENT_PIN", f"{a.consumer}:missing:{path}")
            observed = sha256_file(path)
            require(observed == row["sha256"], "B2_A1_ENVIRONMENT_PIN", f"{a.consumer}:digest:{path}")
            verified.append({"path": row["path"], "sha256": observed, "categories": row["categories"]})
        shim_row = next((row for row in closure["entries"] if row["path"] == str(a.no_network_shim.resolve())), None)
        require(shim_row is not None and sha256_file(a.no_network_shim) == shim_row["sha256"] and "loaded_shared_library" in shim_row["categories"], "B2_A1_ENVIRONMENT_PIN", "network shim is a pinned loaded object")
        dump_yaml(a.output, {"schema_version": "U1_PHASE_B2_ENVIRONMENT_ASSERT_V4", "status": "PASS", "consumer": a.consumer, "stage0_contract_sha256": a.stage0_contract_digest, "stage0_first_use_authentication": auth, "environment_closure_digest": pinned, "environment_identity_digest": pinned, "verified_entry_count": len(verified), "live_equals_pinned_record": True, "read_only_root_sandbox": {"status": "PASS", "mechanism": "bubblewrap --ro-bind / /", "root_mount_options": [row[5] for row in root_mounts]}})
        print(f"B2_ENVIRONMENT_GATE: PASS consumer={a.consumer} entries={len(verified)}"); return 0
    except Exception as exc:
        print(str(exc)); return 1


if __name__ == "__main__":
    raise SystemExit(main())
