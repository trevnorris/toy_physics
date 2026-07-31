#!/usr/bin/env python3
"""Mutation demonstrations for the required able-to-fail controls."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import sympy as sp

from registry_read import AdmissionError, Registry, load_raw_documents, load_registry


HERE = Path(__file__).resolve().parent
CASES = ("vacuous", "duplicate", "entailed", "independent", "provenance")
ERROR_CASE = "no-witness-crash"
SPOOF_CASE = "verdict-spoof"
EXPECTED_FAILURE_EXIT = 1
PARENT_VERDICT_TOKEN = "ABLE_TO_FAIL_HARNESS:"


def _baseline(registry: Registry) -> tuple[tuple, int]:
    constraints = tuple(registry.admitted_constraint_set)
    return constraints, registry.constraint_dimension(constraints)


def demonstrate_vacuous() -> int:
    registry = load_registry()
    baseline, before = _baseline(registry)
    r3 = registry.require_admitted("R3").residual
    assert r3 is not None
    output = registry.symbols[registry.resolve_qid("lambda_gamma")]
    mutated = tuple(expression for expression in baseline if expression != r3) + (output - output,)
    after = registry.constraint_dimension(mutated)
    if after == before:
        print(f"UNEXPECTED_PASS vacuous-flat: before={before} after={after}")
        return 0
    print(
        f"EXPECTED_FAILURE vacuous-flat: forbidden after={before}; "
        f"observed after={after} (increase={after - before})"
    )
    return 1


def demonstrate_duplicate() -> int:
    registry = load_registry()
    baseline, before = _baseline(registry)
    duplicate = registry.require_admitted("R1").residual
    assert duplicate is not None
    after = registry.constraint_dimension(baseline + (duplicate,))
    forbidden = before - 1
    if after == forbidden:
        print(f"UNEXPECTED_PASS duplicate-counted: forbidden={forbidden} observed={after}")
        return 0
    print(
        f"EXPECTED_FAILURE duplicate-counted: forbidden after={forbidden}; "
        f"observed after={after} (unchanged={after == before})"
    )
    return 1


def demonstrate_entailed() -> int:
    registry = load_registry()
    baseline, before = _baseline(registry)
    xi_h = registry.symbols[registry.resolve_qid("xi_h")]
    a_pin = registry.symbols[registry.resolve_qid("a")]
    entailed = xi_h**2 - 2 * a_pin**2
    assert entailed not in baseline
    after = registry.constraint_dimension(baseline + (entailed,))
    forbidden = before - 1
    if after == forbidden:
        print(f"UNEXPECTED_PASS semantic-entailment-counted: forbidden={forbidden} observed={after}")
        return 0
    print(
        f"EXPECTED_FAILURE semantic-entailment-counted: forbidden after={forbidden}; "
        f"observed after={after} (unchanged={after == before})"
    )
    return 1


def demonstrate_independent() -> int:
    registry = load_registry()
    baseline, before = _baseline(registry)
    big_k = registry.symbols[registry.resolve_qid("K")]
    rho0 = registry.symbols[registry.resolve_qid("rho0")]
    after = registry.constraint_dimension(baseline + (big_k - rho0,))
    if after == before:
        print(f"UNEXPECTED_PASS independent-ignored: before={before} after={after}")
        return 0
    print(
        f"EXPECTED_FAILURE independent-ignored: forbidden after={before}; "
        f"observed after={after} (drop={before - after})"
    )
    return 1


def demonstrate_provenance() -> int:
    quantities, relations, schema = load_raw_documents()
    row = next(row for row in relations["relations"] if row["relation_id"] == "R1")
    row["provenance_status"] = "CALIBRATED"
    mutated = Registry.from_documents(
        quantities,
        relations,
        schema,
        ledger_root=HERE.parent,
    )
    try:
        mutated.require_admitted("R1")
    except AdmissionError as exc:
        print(f"EXPECTED_FAILURE calibrated-promotion: {exc}")
        return 1
    print("UNEXPECTED_PASS calibrated-promotion: R1 entered the earned set")
    return 0


def demonstrate_no_witness_crash() -> int:
    """Exercise an unhandled live witness-certificate error in a child."""
    variable = sp.Symbol("no_positive_solution", real=True)
    registry = load_registry()
    registry.certify_positive_real_dimension(
        (variable + 1,),
        (variable,),
        dimension=0,
    )
    raise AssertionError("unreachable after missing positive witness")


def demonstrate_verdict_spoof() -> int:
    """Emit the reserved parent token to exercise child-output framing."""
    print(f"{PARENT_VERDICT_TOKEN} PASS")
    print("EXPECTED_FAILURE verdict-spoof: child output cannot forge the parent verdict")
    return EXPECTED_FAILURE_EXIT


def run_one(case_name: str) -> int:
    observed = {
        "vacuous": demonstrate_vacuous,
        "duplicate": demonstrate_duplicate,
        "entailed": demonstrate_entailed,
        "independent": demonstrate_independent,
        "provenance": demonstrate_provenance,
        ERROR_CASE: demonstrate_no_witness_crash,
        SPOOF_CASE: demonstrate_verdict_spoof,
    }[case_name]()
    marker = "CAUGHT" if observed == EXPECTED_FAILURE_EXIT else "ESCAPED"
    print(f"ABLE_TO_FAIL_{marker}: {case_name}")
    return observed


def classify_child(
    case_name: str,
    completed: subprocess.CompletedProcess[str],
) -> str:
    """Distinguish a deliberate caught failure from an escape or child error."""
    expected_marker = f"ABLE_TO_FAIL_CAUGHT: {case_name}"
    stdout_lines = completed.stdout.splitlines()
    if (
        completed.returncode == EXPECTED_FAILURE_EXIT
        and expected_marker in stdout_lines
        and not completed.stderr.strip()
    ):
        return "CAUGHT"
    if completed.returncode == 0:
        return "ESCAPED"
    return "ERROR"


def echo_child_output(case_name: str, stream_name: str, output: str) -> None:
    """Frame child lines and escape the token reserved for the parent verdict."""
    for line in output.splitlines():
        safe_line = line.replace(PARENT_VERDICT_TOKEN, "ABLE_TO_FAIL_CHILD_TEXT:")
        print(f"CHILD {case_name} {stream_name}: {safe_line}")


def run_harness(case_names: tuple[str, ...] = CASES) -> int:
    statuses: list[str] = []
    valid_case_count = len(case_names) == len(CASES)
    if not valid_case_count:
        print(
            f"HARNESS_CASE_COUNT: ERROR expected={len(CASES)} "
            f"observed={len(case_names)}"
        )
    for case_name in case_names:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--case", case_name],
            check=False,
            capture_output=True,
            text=True,
            timeout=600,
        )
        echo_child_output(case_name, "stdout", completed.stdout)
        echo_child_output(case_name, "stderr", completed.stderr)
        status = classify_child(case_name, completed)
        statuses.append(status)
        print(
            f"DEMONSTRATION {case_name}: status={status} "
            f"observed_exit={completed.returncode}"
        )
    if not valid_case_count:
        result = "ERROR"
    elif all(status == "CAUGHT" for status in statuses):
        result = "PASS"
    elif "ERROR" in statuses:
        result = "ERROR"
    else:
        result = "FAIL"
    print(f"{PARENT_VERDICT_TOKEN} {result}")
    return 0 if result == "PASS" else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=CASES + (ERROR_CASE, SPOOF_CASE))
    parser.add_argument(
        "--demonstrate-crash",
        action="store_true",
        help="run only the no-positive-witness child; success means harness exit 1/ERROR",
    )
    parser.add_argument(
        "--demonstrate-empty",
        action="store_true",
        help="run zero children; success means harness exit 1/ERROR",
    )
    parser.add_argument(
        "--demonstrate-spoof",
        action="store_true",
        help="run a child that emits the reserved parent verdict token",
    )
    arguments = parser.parse_args()
    if arguments.case:
        return run_one(arguments.case)
    if arguments.demonstrate_crash:
        return run_harness((ERROR_CASE,))
    if arguments.demonstrate_empty:
        return run_harness(())
    if arguments.demonstrate_spoof:
        return run_harness((SPOOF_CASE,))
    return run_harness()


if __name__ == "__main__":
    raise SystemExit(main())
