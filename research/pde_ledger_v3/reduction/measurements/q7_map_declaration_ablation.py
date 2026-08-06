#!/usr/bin/env python3
"""Ablate the declared gradient-symbol map and print what the harness reports.

Variants, all applied only to the nine g{r}x{c} <-> g{r}{c} exceptions:

    AS_DECLARED  the committed declaration
    TRANSPOSE    g{r}x{c} -> g{c}{r}
    DERANGED     g{r}x{c} -> g{r}{(c mod 3) + 1}
    DROP_G12     the committed declaration with the g12 entry removed
    ABSENT       all nine entries removed

For each variant the script prints the cross-engine counters and the status of
every declared Q7 row.  It prints what each variant produces; it asserts
nothing about which variant is right.
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
WORK = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/q7_map_ablation")

GRADIENT = re.compile(r"^g[1-3][1-3]$")


def gradient_entries(exceptions: list[dict]) -> list[dict]:
    return [item for item in exceptions if GRADIENT.fullmatch(str(item.get("canonical", "")))]


def variant(base: dict, name: str) -> dict:
    config = copy.deepcopy(base)
    exceptions = config["symbol_naming"]["exceptions"]
    others = [item for item in exceptions if not GRADIENT.fullmatch(str(item.get("canonical", "")))]
    gradients = gradient_entries(exceptions)

    if name == "AS_DECLARED":
        pass
    elif name == "TRANSPOSE":
        for item in gradients:
            row, column = item["canonical"][1], item["canonical"][2]
            item["spellings"]["wl"] = f"g{column}x{row}"
    elif name == "DERANGED":
        for item in gradients:
            row, column = int(item["canonical"][1]), int(item["canonical"][2])
            item["spellings"]["wl"] = f"g{row}x{column % 3 + 1}"
    elif name == "DROP_G12":
        gradients = [item for item in gradients if item["canonical"] != "g12"]
        config["symbol_naming"]["exceptions"] = others + gradients
    elif name == "ABSENT":
        config["symbol_naming"]["exceptions"] = others
    else:
        raise ValueError(name)
    return config


def run(path: Path) -> str:
    result = subprocess.run(
        [
            sys.executable,
            str(HARNESS),
            "--config",
            str(path),
            "--output",
            f"wl={WL_OUT}",
            "--output",
            f"py={PY_OUT}",
        ],
        capture_output=True,
        text=True,
        timeout=600,
    )
    return result.stdout


def main() -> int:
    WORK.mkdir(parents=True, exist_ok=True)
    base = yaml.safe_load(CONFIG.read_text())
    print(f"CONFIG: {CONFIG}")
    print(f"GRADIENT_ENTRIES_IN_CONFIG: {len(gradient_entries(base['symbol_naming']['exceptions']))}")
    print()

    for name in ("AS_DECLARED", "TRANSPOSE", "DERANGED", "DROP_G12", "ABSENT"):
        config = variant(base, name)
        path = WORK / f"checks_S10_{name}.yaml"
        path.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True))
        report = run(path)
        counters = [line for line in report.splitlines() if line.startswith("CROSS_ENGINE:")]
        coverage = [line for line in report.splitlines() if line.startswith("CROSS_ENGINE_COVERAGE:")]
        q7_rows = [
            line.strip()
            for line in report.splitlines()
            if re.match(r"^\s+[a-z0-9_]+_q7_(stiffness|difference):", line)
        ]
        print(f"VARIANT={name}  config={path}")
        for line in counters + coverage:
            print(f"  {line}")
        statuses: dict[str, int] = {}
        for line in q7_rows:
            status = line.split(": ", 1)[1].split(" ", 1)[0]
            statuses[status] = statuses.get(status, 0) + 1
        print(f"  q7_rows={len(q7_rows)} statuses={dict(sorted(statuses.items()))}")
        for line in q7_rows:
            quantity, rest = line.split(": ", 1)
            print(f"    {quantity}: {rest.split(' family=')[0]}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
