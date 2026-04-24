#!/usr/bin/env python3
"""Build a combined referee-facing PDE audit summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

from write_summary_json import parse_output


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def main() -> int:
    parser = argparse.ArgumentParser(description="Write combined PDE audit referee summary")
    parser.add_argument("--audit-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    audit_dir = Path(args.audit_dir)
    output_dir = Path(args.output_dir)
    top_outputs: List[Dict[str, Any]] = [
        parse_output(path)
        for path in sorted(output_dir.glob("*.txt"))
        if path.name != "_summary.txt"
    ]

    python_summary = load_json(audit_dir / "scripts" / "output" / "_summary.json")
    mathematica_summary = load_json(audit_dir / "mathematica" / "output" / "_summary.json")

    total_checks = len(top_outputs)
    top_passed = sum(1 for item in top_outputs if item["status"] == "PASS")
    top_failed = total_checks - top_passed

    root_json = sorted(path.name for path in audit_dir.glob("*.json"))
    referee_pass = (
        top_failed == 0
        and python_summary.get("failed", 1) == 0
        and mathematica_summary.get("failed", 1) == 0
        and not root_json
    )

    combined = {
        "schema": "pde_audit_referee_summary/v1",
        "audit_dir": str(audit_dir),
        "referee_pass": referee_pass,
        "root_json_files": root_json,
        "top_level_checks": top_outputs,
        "python_summary": {
            "path": str(audit_dir / "scripts" / "output" / "_summary.json"),
            "total": python_summary.get("total"),
            "passed": python_summary.get("passed"),
            "failed": python_summary.get("failed"),
        },
        "python_environment": python_summary.get("environment"),
        "mathematica_summary": {
            "path": str(audit_dir / "mathematica" / "output" / "_summary.json"),
            "total": mathematica_summary.get("total"),
            "passed": mathematica_summary.get("passed"),
            "failed": mathematica_summary.get("failed"),
        },
        "mathematica_environment": mathematica_summary.get("environment"),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "_summary.json").open("w", encoding="utf-8") as f:
        json.dump(combined, f, indent=2, sort_keys=True)
        f.write("\n")

    with (output_dir / "_summary.txt").open("w", encoding="utf-8") as f:
        for item in top_outputs:
            f.write(f"{item['status']}  {item['name']}  ({item.get('elapsed', 'unknown')})\n")
        f.write("\n")
        f.write(
            "PYTHON: "
            f"{python_summary.get('passed')}/{python_summary.get('total')} passed, "
            f"{python_summary.get('failed')} failed\n"
        )
        f.write(
            "MATHEMATICA: "
            f"{mathematica_summary.get('passed')}/{mathematica_summary.get('total')} passed, "
            f"{mathematica_summary.get('failed')} failed\n"
        )
        f.write(f"ROOT_JSON_FILES: {len(root_json)}\n")
        f.write(f"REFEREE_PASS: {referee_pass}\n")

    print(f"Wrote {output_dir / '_summary.json'}")
    print(f"REFEREE_PASS: {referee_pass}")
    return 0 if referee_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
