#!/usr/bin/env python3
"""External-copy-safe helper replacing the B1 runner's inline Python snippets."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import yaml


def main() -> int:
    p = argparse.ArgumentParser(); sub = p.add_subparsers(dest="command", required=True)
    agreement = sub.add_parser("agreement"); agreement.add_argument("source", type=Path); agreement.add_argument("target", type=Path)
    validate = sub.add_parser("validate"); validate.add_argument("source", type=Path)
    a = p.parse_args()
    row = yaml.safe_load(a.source.read_text())
    if a.command == "agreement":
        value = {"schema_version": "U1_PHASE_A_AMENDMENT_AGREEMENT_V3", "digest_gate": row["digest_gate"], "semantic_diff_gate": row["semantic_diff_gate"], "phase_a_acceptance_recheck": row["phase_a_acceptance_recheck"], "correction_finding": row["correction_finding"]}
        a.target.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True, width=220)); return 0
    expected = row.pop("resume_validation_sha256"); actual = hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
    if actual != expected or not row["digest_gate"]["agreement"] or row["final_b1_outputs_emitted"]: return 1
    return 0


if __name__ == "__main__": raise SystemExit(main())
