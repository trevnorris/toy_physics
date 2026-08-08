#!/usr/bin/env python3
"""Build weaker W3 implementations and show each pin failing."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


@dataclass(frozen=True)
class Mutation:
    name: str
    apply: Callable[[Path], None]
    command: tuple[str, ...]


def _replace_exact(path: Path, before: str, after: str) -> None:
    source = path.read_text(encoding="utf-8")
    observed_count = source.count(before)
    count_residual = observed_count - 1
    print(f"SOURCE_MUTATION_OPERAND expected-match-count {path.name}: 1")
    print(
        f"SOURCE_MUTATION_OPERAND observed-match-count {path.name}: "
        f"{observed_count}"
    )
    print(f"SOURCE_MUTATION_RESIDUAL match-count {path.name}: {count_residual}")
    if count_residual != 0:
        raise RuntimeError(f"mutation anchor drifted in {path}")
    path.write_text(source.replace(before, after, 1), encoding="utf-8")


def _weaken_symbolic_witness_guard(directory: Path) -> None:
    _replace_exact(
        directory / "dimension_law_check.py",
        "    passed = forms.passed and coverage.passed\n",
        "    passed = forms.passed\n",
    )


def _weaken_integer_resolution(directory: Path) -> None:
    _replace_exact(
        directory / "dimensional_homogeneity_gate.py",
        '''                if isinstance(value, bool) or not isinstance(value, int):
                    raise UndeterminedFinding(
                        f"{qid}: dimension-law binding {parameter!r} cannot resolve "
                        f"to an integer from {binding_qid}: observed={value!r}"
                    )
''',
        '''                if isinstance(value, bool) or not isinstance(value, int):
                    value = reference_values[parameter]
''',
    )


def _weaken_reference_binding(directory: Path) -> None:
    _replace_exact(
        directory / "registry_read.py",
        '''                reference = dict(law.reference_values)[parameter]
                declared = int(target.value)
                reference_residual = declared - reference
                if reference_residual != 0:
                    raise RegistryValidationError(
                        f"{quantity.qid}: dimension-law reference differs from bound "
                        f"quantity {binding_qid}: reference={reference} "
                        f"declared={declared} residual={reference_residual}"
                    )
''',
        "",
    )


def _weaken_reference_vector_binding(directory: Path) -> None:
    _replace_exact(
        directory / "dimension_laws.py",
        '''    if observed_reference != reference_vector:
        raise DimensionLawError(
            "dimension-law reference evaluation differs from exponents: "
            f"law={list(observed_reference)} exponents={list(reference_vector)}"
        )
''',
        "",
    )


def _weaken_structural_target(directory: Path) -> None:
    _replace_exact(
        directory / "dimensional_homogeneity_gate.py",
        '''                if (
                    target.get("kind") != "discrete-choice"
                    or target.get("counting_axis") != "discrete-structural"
                ):
                    raise UndeterminedFinding(
                        f"{qid}: dimension-law binding {parameter!r} targets "
                        f"non-structural quantity {binding_qid}"
                    )
''',
        "",
    )
    _replace_exact(
        directory / "registry_read.py",
        '''                if (
                    target.kind != "discrete-choice"
                    or target.counting_axis != "discrete-structural"
                ):
                    raise RegistryValidationError(
                        f"{quantity.qid}: dimension-law binding {parameter!r} targets "
                        f"non-structural quantity {binding_qid}"
                    )
''',
        "",
    )


def _weaken_runner_signal(directory: Path) -> None:
    target = directory / "test_runner_contract.py"
    present = target.is_file()
    print("RUNNER_SENTINEL_OPERAND expected-present: True")
    print(f"RUNNER_SENTINEL_OPERAND observed-present: {present}")
    print(f"RUNNER_SENTINEL_RESIDUAL equality: {present is True}")
    if not present:
        raise RuntimeError(f"runner sentinel missing before mutation: {target}")
    target.unlink()


def _weaken_aggregate_exit(directory: Path) -> None:
    _replace_exact(
        directory / "dimension_law_able_to_fail.py",
        '    return 0 if result == "PASS" else 1\n',
        '    return 0 if result in {"PASS", "FAIL"} else 1\n',
    )


MUTATIONS = (
    Mutation(
        "object-1-symbolic-witness-ignored",
        _weaken_symbolic_witness_guard,
        (sys.executable, "w3_acceptance_ablations.py"),
    ),
    Mutation(
        "object-2-reference-fallback",
        _weaken_integer_resolution,
        (
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "test_dimension_laws.py::test_gate_is_undetermined_when_a_binding_has_no_integer_value",
        ),
    ),
    Mutation(
        "object-3-reference-binding-removed",
        _weaken_reference_binding,
        (
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "test_dimension_laws.py::test_reference_values_must_equal_the_bound_quantities_declared_values",
        ),
    ),
    Mutation(
        "object-4-retained-triple-guard-removed",
        _weaken_reference_vector_binding,
        (
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "test_dimension_laws.py::test_retained_reference_vector_is_bound_to_the_law_evaluation",
        ),
    ),
    Mutation(
        "object-4-structural-target-guard-removed",
        _weaken_structural_target,
        (
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "test_dimension_laws.py::test_binding_target_must_be_structural",
        ),
    ),
    Mutation(
        "object-4-wrong-runner-signal-removed",
        _weaken_runner_signal,
        (sys.executable, "w3_runner_check.py"),
    ),
    Mutation(
        "object-5-escape-exits-zero",
        _weaken_aggregate_exit,
        (
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "test_dimension_laws.py::test_dimension_law_able_to_fail_escape_fails_the_aggregate",
        ),
    ),
)


def main() -> int:
    results: list[tuple[str, int]] = []
    for mutation in MUTATIONS:
        with tempfile.TemporaryDirectory(
            prefix=".w3_weaker_", dir=REPO / "research"
        ) as temporary:
            reduction = Path(temporary) / "reduction"
            shutil.copytree(
                HERE,
                reduction,
                ignore=shutil.ignore_patterns("__pycache__", ".pytest_cache"),
            )
            print(f"WEAKER_FIX: {mutation.name}")
            mutation.apply(reduction)
            completed = subprocess.run(
                mutation.command,
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
            print("WEAKER_TEST_STDOUT_BEGIN")
            print(completed.stdout, end="")
            print("WEAKER_TEST_STDOUT_END")
            print("WEAKER_TEST_STDERR_BEGIN")
            print(completed.stderr, end="")
            print("WEAKER_TEST_STDERR_END")
            print("WEAKER_TEST_OPERAND expected-exit: 1")
            print(f"WEAKER_TEST_OPERAND observed-exit: {completed.returncode}")
            exit_residual = completed.returncode - 1
            print(f"WEAKER_TEST_RESIDUAL exit: {exit_residual}")
            test_failed = completed.returncode == 1
            print(
                f"WEAKER_FIX_TEST {mutation.name}: "
                f"{'FAIL_AS_REQUIRED' if test_failed else 'SURVIVED'}"
            )
            results.append((mutation.name, completed.returncode))

    expected = tuple((mutation.name, 1) for mutation in MUTATIONS)
    observed = tuple(results)
    residual = tuple(
        (name, observed_exit - expected_exit)
        for (name, observed_exit), (_, expected_exit) in zip(observed, expected)
    )
    print(f"MUTATION_POPULATION_OPERAND expected: {expected}")
    print(f"MUTATION_POPULATION_OPERAND observed: {observed}")
    print(f"MUTATION_POPULATION_RESIDUAL exits: {residual}")
    passed = observed == expected
    print(f"W3_WEAKER_FIX_MUTATIONS: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
