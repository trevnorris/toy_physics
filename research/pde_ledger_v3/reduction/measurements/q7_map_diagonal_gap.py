#!/usr/bin/env python3
"""Does a declaration wrong on the DIAGONAL entries only survive both gates?

A review leg reported that corrupting g11/g22/g33 alone leaves both the unit
test and the ablation battery green, because the six off-diagonal renames are
the only ones either gate pins.  This script measures that claim directly.

Variants, applied only to the three diagonal gradient exceptions:

    AS_DECLARED       unchanged
    DIAGONAL_CYCLED   g{i}{i} <- WL g{i+1 mod 3}x{i+1 mod 3}

For each it prints the declared Q7 rows' status tally, the two predicates the
ablation battery asserts, and the status and applied-rename set of the one row
the unit test pins.  It draws no conclusion about whether the gap should be
closed.
"""

from __future__ import annotations

import copy
import re
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "reduction" / "checks_S10.yaml"
HARNESS = ROOT / "reduction" / "engine_output_checks.py"
WL_OUT = ROOT / "mathematica" / "out" / "S10_brane_mode_spectrum_mathematica_audit.out"
PY_OUT = ROOT / "scripts" / "out" / "S10_brane_mode_spectrum_sympy_audit.out"
WORK = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/q7_diagonal_gap")

GRADIENT = re.compile(r"^g([1-3])([1-3])$")
ROW = re.compile(r"^  ([a-z0-9_]+_q7_(?:stiffness|difference)): ([A-Z_]+) family=\S+ naming=\[([^\]]*)\]")

# The row the unit test pins, and the row whose payload carries all nine symbols.
PINNED_ROW = "main_d3_q7_stiffness"
ALL_NINE_ROW = "xform_fullgrad_d3_q7_stiffness"


def build(variant: str) -> Path:
    config = copy.deepcopy(yaml.safe_load(CONFIG.read_text()))
    for item in config["symbol_naming"]["exceptions"]:
        match = GRADIENT.fullmatch(str(item.get("canonical", "")))
        if not match:
            continue
        row, column = int(match.group(1)), int(match.group(2))
        if variant == "DIAGONAL_CYCLED" and row == column:
            shifted = row % 3 + 1
            item["spellings"]["wl"] = f"g{shifted}x{shifted}"
    WORK.mkdir(parents=True, exist_ok=True)
    path = WORK / f"checks_S10_{variant}.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True))
    return path


def measure(path: Path) -> tuple[dict[str, int], dict[str, tuple[str, str]]]:
    result = subprocess.run(
        [sys.executable, str(HARNESS), "--config", str(path),
         "--output", f"wl={WL_OUT}", "--output", f"py={PY_OUT}"],
        capture_output=True, text=True, timeout=600,
    )
    rows: dict[str, tuple[str, str]] = {}
    for line in result.stdout.splitlines():
        match = ROW.match(line)
        if match:
            rows[match.group(1)] = (match.group(2), match.group(3))
    tally: dict[str, int] = {}
    for status, _ in rows.values():
        tally[status] = tally.get(status, 0) + 1
    return dict(sorted(tally.items())), rows


def main() -> int:
    results = {}
    for variant in ("AS_DECLARED", "DIAGONAL_CYCLED"):
        path = build(variant)
        tally, rows = measure(path)
        results[variant] = tally
        print(f"VARIANT={variant}  config={path}")
        print(f"  q7_tally = {tally}")
        for name in (PINNED_ROW, ALL_NINE_ROW):
            status, naming = rows.get(name, ("<absent>", ""))
            print(f"  {name}: status={status}")
            print(f"    naming_applied = [{naming}]")
        print()

    print("THE TWO PREDICATES THE ABLATION BATTERY ASSERTS, evaluated with")
    print("DIAGONAL_CYCLED standing in for AS_DECLARED:")
    print(f"  DERANGED != AS_DECLARED : requires a separate DERANGED run; see the battery")
    print(f"  AS_DECLARED tally under DIAGONAL_CYCLED = {results['DIAGONAL_CYCLED']}")
    print(f"  AS_DECLARED tally as committed          = {results['AS_DECLARED']}")
    print(f"  tallies_equal = {results['DIAGONAL_CYCLED'] == results['AS_DECLARED']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
