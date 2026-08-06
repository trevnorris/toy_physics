#!/usr/bin/env python3
"""Measure how much of the cross-engine agreement each declaration carries.

For each declared naming exception and each declared symbol identity, this
script removes that one declaration, re-runs the comparator over the committed
outputs, and prints the resulting counters beside the undisturbed ones together
with the rows whose status changed.

It reports what each removal does.  It does not say whether any declaration is
justified -- that question is answered from the engine definitions, not from
these counts.

Usage:
    declaration_load_ablation.py S9|S10 [workdir]
"""

from __future__ import annotations

import copy
import re
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
HARNESS = ROOT / "reduction" / "engine_output_checks.py"
STEPS = {
    "S9": (
        ROOT / "reduction" / "checks_S9.yaml",
        ROOT / "mathematica" / "out" / "S9_light_requires_shear_mathematica_audit.out",
        ROOT / "scripts" / "out" / "S9_light_requires_shear_sympy_audit.out",
    ),
    "S10": (
        ROOT / "reduction" / "checks_S10.yaml",
        ROOT / "mathematica" / "out" / "S10_brane_mode_spectrum_mathematica_audit.out",
        ROOT / "scripts" / "out" / "S10_brane_mode_spectrum_sympy_audit.out",
    ),
}

ROW = re.compile(r"^  ([a-z0-9_]+): ([A-Z_]+) family=")


def run(config_path: Path, wl: Path, py: Path) -> tuple[str, dict[str, str]]:
    result = subprocess.run(
        [
            sys.executable,
            str(HARNESS),
            "--config",
            str(config_path),
            "--output",
            f"wl={wl}",
            "--output",
            f"py={py}",
        ],
        capture_output=True,
        text=True,
        timeout=600,
    )
    counters = ""
    rows: dict[str, str] = {}
    in_rows = False
    for line in result.stdout.splitlines():
        if line.startswith("CROSS_ENGINE:"):
            counters = line
        if line.startswith("CROSS_ENGINE_ROWS"):
            in_rows = True
            continue
        if in_rows:
            match = ROW.match(line)
            if match:
                rows[match.group(1)] = match.group(2)
            elif line and not line.startswith(" "):
                in_rows = False
    return counters, rows


def main() -> int:
    step = sys.argv[1] if len(sys.argv) > 1 else "S10"
    work = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("/tmp/declaration_ablation")
    config_path, wl, py = STEPS[step]
    work.mkdir(parents=True, exist_ok=True)

    base = yaml.safe_load(config_path.read_text())
    baseline_counters, baseline_rows = run(config_path, wl, py)
    print(f"STEP={step} CONFIG={config_path}")
    print(f"BASELINE {baseline_counters}")
    print()

    naming = base.get("symbol_naming", {}).get("exceptions", []) or []
    identities = base.get("symbol_identities", []) or []
    print(f"DECLARED_NAMING_EXCEPTIONS={len(naming)} DECLARED_SYMBOL_IDENTITIES={len(identities)}")
    print()

    for index, item in enumerate(naming):
        variant = copy.deepcopy(base)
        removed = variant["symbol_naming"]["exceptions"].pop(index)
        path = work / f"{step}_drop_naming_{removed['canonical']}.yaml"
        path.write_text(yaml.safe_dump(variant, sort_keys=False, allow_unicode=True))
        counters, rows = run(path, wl, py)
        changed = sorted(
            name for name, status in rows.items() if baseline_rows.get(name) != status
        )
        print(f"DROP naming[{removed['canonical']}]  {counters}")
        print(f"  changed_rows={len(changed)}")
        if changed:
            transitions: dict[str, int] = {}
            for name in changed:
                key = f"{baseline_rows.get(name)}->{rows[name]}"
                transitions[key] = transitions.get(key, 0) + 1
            print(f"  transitions={transitions}")
            print(f"  rows={', '.join(changed[:8])}{' ...' if len(changed) > 8 else ''}")

    for index, item in enumerate(identities):
        variant = copy.deepcopy(base)
        removed = variant["symbol_identities"].pop(index)
        path = work / f"{step}_drop_identity_{removed['symbol']}.yaml"
        path.write_text(yaml.safe_dump(variant, sort_keys=False, allow_unicode=True))
        counters, rows = run(path, wl, py)
        changed = sorted(
            name for name, status in rows.items() if baseline_rows.get(name) != status
        )
        print(f"DROP identity[{removed['engine']}:{removed['symbol']}]  {counters}")
        print(f"  changed_rows={len(changed)}")
        if changed:
            transitions = {}
            for name in changed:
                key = f"{baseline_rows.get(name)}->{rows[name]}"
                transitions[key] = transitions.get(key, 0) + 1
            print(f"  transitions={transitions}")
            print(f"  rows={', '.join(changed[:8])}{' ...' if len(changed) > 8 else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
