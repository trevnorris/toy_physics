#!/usr/bin/env python3
"""Synchronize generated graph formats.

The YAML graph is the maintained source. The JSON graph is generated from YAML
for tooling, browser embedding, and jq-friendly inspection.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import yaml


GRAPH_DIR = Path(__file__).resolve().parent
DEFAULT_YAML = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.yaml"
DEFAULT_JSON = GRAPH_DIR / "fluid_universe_derivation_atlas_graph.json"


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_yaml(path: Path) -> Any:
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        fail(f"{path} is not valid YAML: {exc}")


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        fail(f"{path} is not valid JSON: {exc}")


def write_json(path: Path, graph: Any) -> None:
    path.write_text(json.dumps(graph, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--yaml", type=Path, default=DEFAULT_YAML, help="Primary graph YAML path")
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON, help="Generated graph JSON path")
    parser.add_argument("--check", action="store_true", help="Fail if JSON is not in sync with YAML")
    args = parser.parse_args()

    graph = load_yaml(args.yaml)

    if args.check:
        if not args.json.exists():
            fail(f"{args.json} does not exist")
        generated = load_json(args.json)
        if generated != graph:
            fail(f"{args.json} is stale; run python3 graph/sync_graph_formats.py")
        print("OK: graph JSON is in sync with YAML")
        return

    write_json(args.json, graph)
    print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
