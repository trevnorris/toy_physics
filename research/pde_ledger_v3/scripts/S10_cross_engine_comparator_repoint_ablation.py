#!/usr/bin/env python3
"""Re-point one shared PY name to a neighbouring PY object and show the row move."""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from S10_cross_engine_comparator import (
    ComparatorInputError,
    compare_records,
    load_pair,
    render_comparison,
    render_value,
)


TARGET_NAME = "S10_MAIN_D3_Q2_MATRIX_A"
NEIGHBOUR_NAME = "S10_MAIN_D3_Q2_MATRIX_B"


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Ablate one comparator binding without modifying either engine transcript."
    )
    parser.add_argument("out_files", nargs=2, type=Path, metavar="OUT")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        py, wl = load_pair(args.out_files)
        py_target = py.records[TARGET_NAME]
        py_neighbour = py.records[NEIGHBOUR_NAME]
        wl_target = wl.records[TARGET_NAME]
    except (ComparatorInputError, KeyError) as error:
        print(f"ABLATION_INPUT_ERROR: {error}")
        print("ABLATION_GUARD: FAIL")
        return 2

    baseline = compare_records(TARGET_NAME, py_target, wl_target)
    repointed_record = replace(py_target, raw=py_neighbour.raw)
    ablated = compare_records(TARGET_NAME, repointed_record, wl_target)

    print(f"ABLATION_TARGET_NAME: {TARGET_NAME}")
    print(f"ABLATION_REPOINT_SOURCE_NAME: {NEIGHBOUR_NAME}")
    print("CATEGORY: BASELINE_ROW")
    for line in render_comparison(baseline):
        print(line)
    print("CATEGORY: REPOINTED_ROW")
    for line in render_comparison(ablated):
        print(line)

    residual_moved = (
        baseline.reason is None
        and ablated.reason is None
        and render_value(baseline.residual) != render_value(ablated.residual)
    )
    guard_moved = baseline.passes and not ablated.passes
    print(f"RESIDUAL_MOVED: {residual_moved}")
    print(f"ROW_GUARD_MOVED: {guard_moved}")
    passed = residual_moved and guard_moved
    print(f"ABLATION_GUARD: {'PASS' if passed else 'FAIL'}")
    return int(not passed)


if __name__ == "__main__":
    raise SystemExit(main())
