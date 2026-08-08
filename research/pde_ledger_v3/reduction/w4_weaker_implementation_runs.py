#!/usr/bin/env python3
"""Apply each D4 weaker implementation and require its focused test to fail."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


@dataclass(frozen=True)
class WeakerImplementation:
    name: str
    path: str
    before: str
    after: str
    test: str


IMPLEMENTATIONS = (
    WeakerImplementation(
        "D4-1-acceptance-ignores-coefficient-subverdict",
        "w3_acceptance_ablations.py",
        "        case_passed = not any(residual.values())\n",
        '''        case_passed = not any(
            value
            for field, value in residual.items()
            if field != "coefficient_pin_passed"
        )
''',
        "test_dimension_laws.py::test_acceptance_subverdict_is_guarded",
    ),
    WeakerImplementation(
        "D4-2-demonstrate-escape-types-observation",
        "dimension_law_able_to_fail.py",
        '''def demonstrate_escape() -> int:
    weaker_case = "wrong-law"
    completed = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--case",
            weaker_case,
            "--weaken-case-exit",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    _echo_child(weaker_case, "weaker-stdout", completed.stdout)
    _echo_child(weaker_case, "weaker-stderr", completed.stderr)
    expected_exit = EXPECTED_FAILURE_EXIT
    observed_exit = completed.returncode
    residual = observed_exit - expected_exit
    status = classify_child(weaker_case, completed)
    print(f"ESCAPE_OPERAND expected-child-exit: {expected_exit}")
    print(f"ESCAPE_OPERAND observed-child-exit: {observed_exit}")
    print(f"ESCAPE_RESIDUAL observed-minus-expected: {residual}")
    print(f"DIMENSION_LAW_ABLE_TO_FAIL_{status}: {ESCAPE_CASE}")
    return observed_exit
''',
        '''def demonstrate_escape() -> int:
    print("ESCAPE_OPERAND expected-child-exit: 1")
    print("ESCAPE_OPERAND observed-child-exit: 0")
    print("ESCAPE_RESIDUAL observed-minus-expected: -1")
    print(f"DIMENSION_LAW_ABLE_TO_FAIL_ESCAPED: {ESCAPE_CASE}")
    return 0
''',
        "test_dimension_laws.py::test_demonstrate_escape_uses_the_observed_child_process",
    ),
    WeakerImplementation(
        "D4-3-status-print-is-literal-pass",
        "dimension_law_check.py",
        '''    print(f"DIMENSION_LAW_CHECK: {'PASS' if passed else 'FAIL'}")
''',
        '''    print("DIMENSION_LAW_CHECK: PASS")
''',
        "test_dimension_laws.py::test_printed_dimension_law_status_is_guarded",
    ),
)


def _apply(directory: Path, implementation: WeakerImplementation) -> None:
    path = directory / implementation.path
    source = path.read_text(encoding="utf-8")
    observed_count = source.count(implementation.before)
    print(f"WEAKER_SOURCE_OPERAND expected-anchor-count {implementation.name}: 1")
    print(
        f"WEAKER_SOURCE_OPERAND observed-anchor-count {implementation.name}: "
        f"{observed_count}"
    )
    print(
        f"WEAKER_SOURCE_RESIDUAL anchor-count {implementation.name}: "
        f"{observed_count - 1}"
    )
    if observed_count != 1:
        raise RuntimeError(f"mutation anchor drifted: {path}")
    path.write_text(
        source.replace(implementation.before, implementation.after, 1),
        encoding="utf-8",
    )


def main() -> int:
    results: list[tuple[str, int]] = []
    for implementation in IMPLEMENTATIONS:
        with tempfile.TemporaryDirectory(prefix=".w4_weaker_", dir=REPO / "research") as temporary:
            reduction = Path(temporary) / "reduction"
            shutil.copytree(
                HERE,
                reduction,
                ignore=shutil.ignore_patterns("__pycache__", ".pytest_cache"),
            )
            print(f"WEAKER_IMPLEMENTATION: {implementation.name}")
            _apply(reduction, implementation)
            completed = subprocess.run(
                [sys.executable, "-m", "pytest", "-q", implementation.test],
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
            print(f"WEAKER_TEST_STDOUT_BEGIN {implementation.name}")
            print(completed.stdout, end="")
            print(f"WEAKER_TEST_STDOUT_END {implementation.name}")
            print(f"WEAKER_TEST_STDERR_BEGIN {implementation.name}")
            print(completed.stderr, end="")
            print(f"WEAKER_TEST_STDERR_END {implementation.name}")
            print(f"WEAKER_TEST_OPERAND observed-exit {implementation.name}: {completed.returncode}")
            failed_as_required = completed.returncode == 1
            print(
                f"WEAKER_TEST_GUARD {implementation.name}: "
                f"{'FAIL_AS_REQUIRED' if failed_as_required else 'SURVIVED'}"
            )
            results.append((implementation.name, completed.returncode))
    passed = all(exit_code == 1 for _name, exit_code in results)
    print(f"W4_WEAKER_IMPLEMENTATIONS: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
