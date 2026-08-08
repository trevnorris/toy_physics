#!/usr/bin/env python3
"""Run the intact/absent/wrong-coefficient W3 acceptance surface."""

from __future__ import annotations

import copy
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import yaml


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CASES = ("laws-intact", "all-laws-absent", "wrong-D-coefficient")


@dataclass(frozen=True)
class AcceptanceObservation:
    case_name: str
    law_exit_zero: bool
    gate_exit_zero: bool
    declaration_form_passed: bool
    coefficient_pin_passed: bool
    combined_check_passed: bool


def _quantity(document: dict[str, Any], qid: str) -> dict[str, Any]:
    return next(row for row in document["quantities"] if row["qid"] == qid)


def _mutate(case_name: str, document: dict[str, Any]) -> None:
    rho = _quantity(document, "Q.brane.rho_br")
    mu = _quantity(document, "Q.brane.mu_R")
    bulk = _quantity(document, "Q.brane.B_comp")
    if case_name == "all-laws-absent":
        for row in (rho, mu, bulk):
            row["dimension"].pop("law")
    elif case_name == "wrong-D-coefficient":
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


def _run(command: str, directory: Path) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PDE_LEDGER_ENGINE_OUTPUT_ROOT"] = str(HERE.parent)
    completed = subprocess.run(
        [sys.executable, command],
        cwd=directory,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
        env=environment,
    )
    print(f"COMMAND {command} STDOUT BEGIN")
    print(completed.stdout, end="")
    print(f"COMMAND {command} STDOUT END")
    print(f"COMMAND {command} STDERR BEGIN")
    print(completed.stderr, end="")
    print(f"COMMAND {command} STDERR END")
    print(f"COMMAND {command} EXIT: {completed.returncode}")
    return completed


def _printed_boolean(stdout: str, label: str, true_word: str = "PASS") -> bool:
    prefix = f"{label}: "
    values = [line.removeprefix(prefix) for line in stdout.splitlines() if line.startswith(prefix)]
    if len(values) != 1:
        raise RuntimeError(f"expected one {label!r} verdict, observed {values!r}")
    return values[0] == true_word


def _observe(
    case_name: str,
    law: subprocess.CompletedProcess[str],
    gate: subprocess.CompletedProcess[str],
) -> AcceptanceObservation:
    return AcceptanceObservation(
        case_name,
        law.returncode == 0,
        gate.returncode == 0,
        _printed_boolean(law.stdout, "DECLARATION_FORM_CHECK"),
        _printed_boolean(law.stdout, "D_COEFFICIENT_POLICED_IN_REDUCTION", "YES"),
        _printed_boolean(law.stdout, "DIMENSION_LAW_CHECK"),
    )


def _expected(case_name: str) -> AcceptanceObservation:
    declaration_form_passed = case_name != "all-laws-absent"
    coefficient_pin_passed = case_name == "laws-intact"
    return AcceptanceObservation(
        case_name,
        coefficient_pin_passed,
        True,
        declaration_form_passed,
        coefficient_pin_passed,
        coefficient_pin_passed,
    )


def evaluate(observations: list[AcceptanceObservation]) -> bool:
    population_residual = tuple(observation.case_name for observation in observations) != CASES
    print(f"ACCEPTANCE_POPULATION_OPERAND configured: {list(CASES)}")
    print(
        "ACCEPTANCE_POPULATION_OPERAND observed: "
        f"{[observation.case_name for observation in observations]}"
    )
    print(f"ACCEPTANCE_POPULATION_RESIDUAL: {population_residual}")
    passed = not population_residual
    for observation in observations:
        expected = _expected(observation.case_name)
        residual = {
            field: getattr(observation, field) != getattr(expected, field)
            for field in (
                "law_exit_zero",
                "gate_exit_zero",
                "declaration_form_passed",
                "coefficient_pin_passed",
                "combined_check_passed",
            )
        }
        print(f"ACCEPTANCE_CASE_OPERAND expected {observation.case_name}: {expected}")
        print(f"ACCEPTANCE_CASE_OPERAND observed {observation.case_name}: {observation}")
        print(f"ACCEPTANCE_CASE_RESIDUAL {observation.case_name}: {residual}")
        case_passed = not any(residual.values())
        print(
            f"ACCEPTANCE_CASE_GUARD {observation.case_name}: "
            f"{'PASS' if case_passed else 'FAIL'}"
        )
        passed = passed and case_passed
    print(f"W3_ACCEPTANCE_ABLATIONS: {'PASS' if passed else 'FAIL'}")
    return passed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demonstrate-subverdict-escape", action="store_true")
    arguments = parser.parse_args()
    with (HERE / "quantities.yaml").open("r", encoding="utf-8") as handle:
        baseline = yaml.safe_load(handle)
    observed: list[AcceptanceObservation] = []
    for case_name in CASES:
        fixture = copy.deepcopy(baseline)
        _mutate(case_name, fixture)
        with tempfile.TemporaryDirectory(
            prefix=".w3_acceptance_", dir=REPO / "research"
        ) as temporary:
            reduction = Path(temporary) / "reduction"
            shutil.copytree(HERE, reduction)
            with (reduction / "quantities.yaml").open(
                "w", encoding="utf-8"
            ) as handle:
                yaml.safe_dump(fixture, handle, sort_keys=False)
            print(f"CASE: {case_name}")
            law = _run("dimension_law_check.py", reduction)
            gate = _run("dimensional_homogeneity_gate.py", reduction)
            observed.append(_observe(case_name, law, gate))

    if arguments.demonstrate_subverdict_escape:
        observed = [
            replace(
                observation,
                coefficient_pin_passed=not observation.coefficient_pin_passed,
            )
            if observation.case_name == "wrong-D-coefficient"
            else observation
            for observation in observed
        ]
        print("DEMONSTRATED_WEAKER_SUBVERDICT: toggled observed coefficient verdict")
    passed = evaluate(observed)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
