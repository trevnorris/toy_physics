#!/usr/bin/env python3
"""Build four weaker pin implementations and require each focused test to fail."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from pytest import ExitCode


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


@dataclass(frozen=True)
class WeakerImplementation:
    name: str
    before: str
    after: str
    test: str


IMPLEMENTATIONS = (
    WeakerImplementation(
        "drop-covered-quantity-check",
        "            and self.covered_qids == REQUIRED_QIDS\n",
        "",
        "test_engine_dimension_pin_rejects_empty_configured_population",
    ),
    WeakerImplementation(
        "drop-population-count-check",
        "            and len(self.observations) == expected_count\n",
        "",
        "test_engine_dimension_pin_rejects_subset_population_with_required_qids",
    ),
    WeakerImplementation(
        "drop-errors-check",
        "            not self.errors\n            and len(self.observations)",
        "            len(self.observations)",
        "test_engine_dimension_pin_rejects_nonempty_errors",
    ),
    WeakerImplementation(
        "drop-unmapped-symbol-raise",
        '''    if unexpected:
        raise ValueError(f"unmapped engine dimension symbols: {sorted(map(str, unexpected))}")
''',
        "",
        "test_engine_dimension_pin_rejects_unmapped_symbol",
    ),
)


def _apply(reduction: Path, implementation: WeakerImplementation) -> Path:
    path = reduction / "engine_dimension_pin.py"
    source = path.read_text(encoding="utf-8")
    anchor_count = source.count(implementation.before)
    print(
        f"PIN_WEAKER_SOURCE_OPERAND anchor-count {implementation.name}: "
        f"{anchor_count}"
    )
    if anchor_count != 1:
        raise RuntimeError(f"mutation anchor drifted: {path}")
    path.write_text(
        source.replace(implementation.before, implementation.after, 1),
        encoding="utf-8",
    )
    return path


def main() -> int:
    results: list[bool] = []
    for implementation in IMPLEMENTATIONS:
        with tempfile.TemporaryDirectory(
            prefix=".w4_pin_completeness_", dir=REPO / "research"
        ) as temporary:
            reduction = Path(temporary) / "reduction"
            shutil.copytree(
                HERE,
                reduction,
                ignore=shutil.ignore_patterns("__pycache__", ".pytest_cache"),
            )
            print(f"PIN_WEAKER_IMPLEMENTATION: {implementation.name}")
            mutated_path = _apply(reduction, implementation)
            print(f"PIN_WEAKER_IMPLEMENTATION_PATH: {mutated_path}")
            test_target = (
                f"{reduction / 'test_dimension_laws.py'}::{implementation.test}"
            )
            completed = subprocess.run(
                [sys.executable, "-m", "pytest", "-q", test_target],
                cwd=reduction,
                check=False,
                capture_output=True,
                text=True,
                timeout=180,
                env={
                    **os.environ,
                    "PDE_LEDGER_ENGINE_OUTPUT_ROOT": str(HERE.parent),
                },
            )
            print(f"PIN_WEAKER_TEST_STDOUT_BEGIN {implementation.name}")
            print(completed.stdout, end="")
            print(f"PIN_WEAKER_TEST_STDOUT_END {implementation.name}")
            print(f"PIN_WEAKER_TEST_STDERR_BEGIN {implementation.name}")
            print(completed.stderr, end="")
            print(f"PIN_WEAKER_TEST_STDERR_END {implementation.name}")
            print(
                f"PIN_WEAKER_TEST_OPERAND observed-exit {implementation.name}: "
                f"{completed.returncode}"
            )
            failed_as_required = completed.returncode == int(ExitCode.TESTS_FAILED)
            print(
                f"PIN_WEAKER_TEST_GUARD {implementation.name}: "
                f"{'FAIL_AS_REQUIRED' if failed_as_required else 'SURVIVED'}"
            )
            results.append(failed_as_required)
    passed = all(results)
    print(f"W4_PIN_COMPLETENESS_WEAKER_IMPLEMENTATIONS: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
