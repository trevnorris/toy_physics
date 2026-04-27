#!/usr/bin/env python3
"""Regression tests for graph-math normalization in the reader-site generator."""

from __future__ import annotations

import argparse
import difflib
import sys
from pathlib import Path
from typing import Any

import yaml

import generate_reader_site as reader


ATLAS_DIR = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURE = ATLAS_DIR / "fixtures" / "reader_math_cases.yaml"


def load_yaml(path: Path) -> Any:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def diff(expected: str, actual: str) -> str:
    return "\n".join(
        difflib.unified_diff(
            [expected + "\n"],
            [actual + "\n"],
            fromfile="expected",
            tofile="actual",
            lineterm="",
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    args = parser.parse_args()

    data = load_yaml(args.fixture)
    cases = data.get("cases", []) if isinstance(data, dict) else []
    if not cases:
        print(f"No reader math cases found in {args.fixture}", file=sys.stderr)
        sys.exit(1)

    failures = 0
    for case in cases:
        actual = reader.normalize_reader_math(case["input"])
        expected = case["expected"]
        if actual != expected:
            failures += 1
            print(f"FAIL: {case['name']}", file=sys.stderr)
            print(diff(expected, actual), file=sys.stderr)

    if failures:
        print(f"Reader math normalization failed: {failures}/{len(cases)} cases", file=sys.stderr)
        sys.exit(1)

    print(f"OK: reader math normalization fixtures passed ({len(cases)} cases)")


if __name__ == "__main__":
    main()
