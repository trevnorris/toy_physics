#!/usr/bin/env python3
"""Mutation demonstrations for four required able-to-fail controls."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from registry_read import AdmissionError, Registry, load_raw_documents, load_registry


HERE = Path(__file__).resolve().parent
CASES = ("vacuous", "duplicate", "independent", "provenance")


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


def run_one(case_name: str) -> int:
    return {
        "vacuous": demonstrate_vacuous,
        "duplicate": demonstrate_duplicate,
        "independent": demonstrate_independent,
        "provenance": demonstrate_provenance,
    }[case_name]()


def run_harness() -> int:
    all_caught = True
    for case_name in CASES:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--case", case_name],
            check=False,
            capture_output=True,
            text=True,
            timeout=600,
        )
        output = completed.stdout.strip()
        if output:
            print(output)
        if completed.stderr.strip():
            print(completed.stderr.strip())
        print(f"DEMONSTRATION {case_name}: observed_exit={completed.returncode}")
        if completed.returncode != 1:
            all_caught = False
    print(f"ABLE_TO_FAIL_HARNESS: {'PASS' if all_caught else 'FAIL'}")
    return 0 if all_caught else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=CASES)
    arguments = parser.parse_args()
    return run_one(arguments.case) if arguments.case else run_harness()


if __name__ == "__main__":
    raise SystemExit(main())
