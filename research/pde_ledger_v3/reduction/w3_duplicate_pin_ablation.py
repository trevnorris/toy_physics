#!/usr/bin/env python3
"""Show that the committed engine operand catches a registry-only mutation."""

from __future__ import annotations

import copy
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import yaml


HERE = Path(__file__).resolve().parent


def _quantity(document: dict[str, Any], qid: str) -> dict[str, Any]:
    return next(row for row in document["quantities"] if row["qid"] == qid)


def _wrong_coefficient(document: dict[str, Any]) -> None:
    rho = _quantity(document, "Q.brane.rho_br")
    mu = _quantity(document, "Q.brane.mu_R")
    bulk = _quantity(document, "Q.brane.B_comp")
    rho["dimension"]["law"]["components"][0] = [
        "Sub",
        3,
        ["Mul", 2, ["Ref", "D"]],
    ]
    for row in (mu, bulk):
        row["dimension"]["law"]["components"][0] = [
            "Sub",
            5,
            ["Mul", 2, ["Ref", "D"]],
        ]


def _run(directory: Path) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PDE_LEDGER_ENGINE_OUTPUT_ROOT"] = str(HERE.parent)
    return subprocess.run(
        [sys.executable, str(directory / "dimension_law_check.py")],
        cwd=directory,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
        env=environment,
    )


def _print_run(label: str, completed: subprocess.CompletedProcess[str]) -> None:
    print(f"{label}_STDOUT_BEGIN")
    print(completed.stdout, end="")
    print(f"{label}_STDOUT_END")
    print(f"{label}_STDERR_BEGIN")
    print(completed.stderr, end="")
    print(f"{label}_STDERR_END")
    print(f"{label}_EXIT: {completed.returncode}")


def main() -> int:
    with (HERE / "quantities.yaml").open("r", encoding="utf-8") as handle:
        quantities = yaml.safe_load(handle)

    with tempfile.TemporaryDirectory(
        prefix=".w3_fix2_duplicate_", dir=HERE.parents[1]
    ) as temporary:
        baseline_dir = Path(temporary) / "baseline"
        wrong_dir = Path(temporary) / "wrong"
        shutil.copytree(HERE, baseline_dir)
        shutil.copytree(HERE, wrong_dir)
        wrong = copy.deepcopy(quantities)
        _wrong_coefficient(wrong)
        with (wrong_dir / "quantities.yaml").open("w", encoding="utf-8") as handle:
            yaml.safe_dump(wrong, handle, sort_keys=False)
        baseline = _run(baseline_dir)
        changed = _run(wrong_dir)

    _print_run("BASELINE", baseline)
    _print_run("WRONG_REGISTRY_ONLY", changed)
    print(f"DUPLICATE_PIN_EXIT_OPERAND baseline: {baseline.returncode}")
    print(f"DUPLICATE_PIN_EXIT_OPERAND changed: {changed.returncode}")
    print(
        "DUPLICATE_PIN_EXIT_RESIDUAL changed-minus-baseline: "
        f"{changed.returncode - baseline.returncode}"
    )
    passed = (
        baseline.returncode == 0
        and changed.returncode != 0
        and "D_COEFFICIENT_POLICED_IN_REDUCTION: YES" in baseline.stdout
        and "D_COEFFICIENT_POLICED_IN_REDUCTION: NO" in changed.stdout
    )
    print(f"W3_DUPLICATE_PIN_ABLATION: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
