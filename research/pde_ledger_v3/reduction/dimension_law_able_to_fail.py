#!/usr/bin/env python3
"""Mutation population and aggregate harness for bound dimension laws."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping

import sympy as sp

from dimensional_homogeneity_gate import (
    HOMOGENEOUS,
    INHOMOGENEOUS,
    UNDETERMINED,
    DimensionAudit,
    UndeterminedFinding,
)
from dimension_law_check import (
    inspect_declaration_forms,
    inspect_symbolic_witness_coverage,
    print_dimension_law_check,
)
from dimension_laws import SymbolicDimension, dimension_residual
from registry_read import Registry, RegistryValidationError, load_raw_documents


CASES = (
    "absent-laws",
    "wrong-D-coefficient",
    "missing-binding-value",
    "noninteger-binding-value",
    "reference-value-mismatch",
    "reference-vector-mismatch",
    "nonstructural-binding",
    "wrong-law",
    "wrong-binding",
    "unbound",
    "unresolvable",
)
ESCAPE_CASE = "aggregate-escape"
EXPECTED_FAILURE_EXIT = 1
PARENT_VERDICT_TOKEN = "DIMENSION_LAW_ABLE_TO_FAIL_HARNESS:"
EXPECTED_CASE_STATUS = {case_name: "CAUGHT" for case_name in CASES}


def _quantity(document: dict[str, Any], qid: str) -> dict[str, Any]:
    return next(row for row in document["quantities"] if row["qid"] == qid)


def _render(dimension: SymbolicDimension | None, detail: str | None = None) -> str:
    if dimension is None:
        return f"<UNDETERMINED: {detail}>"
    return "[" + ",".join(sp.sstr(sp.simplify(value)) for value in dimension) + "]"


def _dimension_or_finding(
    audit: DimensionAudit, node: Any
) -> tuple[SymbolicDimension | None, str | None]:
    try:
        return audit.dimension_of(node), None
    except UndeterminedFinding as finding:
        return None, str(finding)


def _print_relation_probe(
    audit: DimensionAudit,
    relations: Mapping[str, Any],
    relation_id: str = "R4",
) -> tuple[dict[str, str], SymbolicDimension | None]:
    relation = next(
        row for row in relations["relations"] if row["relation_id"] == relation_id
    )
    left_node, right_node = relation["residual"][1:]
    left, left_detail = _dimension_or_finding(audit, left_node)
    right, right_detail = _dimension_or_finding(audit, right_node)
    if left is not None and right is not None:
        residual = dimension_residual(left, right)
        residual_detail = None
    else:
        residual = None
        residual_detail = left_detail or right_detail
    print(f"LEFT_OPERAND {relation_id}: {_render(left, left_detail)}")
    print(f"RIGHT_OPERAND {relation_id}: {_render(right, right_detail)}")
    print(
        f"DIMENSION_RESIDUAL {relation_id}: "
        f"{_render(residual, residual_detail)}"
    )
    results = audit.classify_relations()
    statuses = {result.relation_id: result.status for result in results}
    for selected in ("R4", "R5"):
        print(f"RELATION {selected}: {statuses[selected]}")
    return statuses, residual


def _loader_finding(
    quantities: Mapping[str, Any],
    relations: Mapping[str, Any],
    schema: Mapping[str, Any],
) -> str | None:
    try:
        Registry.from_documents(
            quantities,
            relations,
            schema,
            ledger_root=Path(__file__).resolve().parents[3],
        )
    except RegistryValidationError as finding:
        return str(finding)
    return None


def _mutate(case_name: str, quantities: dict[str, Any]) -> None:
    rho = _quantity(quantities, "Q.brane.rho_br")
    mu = _quantity(quantities, "Q.brane.mu_R")
    bulk = _quantity(quantities, "Q.brane.B_comp")
    structural = _quantity(quantities, "Q.brane.D_brane")
    if case_name == "absent-laws":
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
    elif case_name == "missing-binding-value":
        structural.pop("value")
    elif case_name == "noninteger-binding-value":
        structural["value"] = ["Rat", 3, 2]
    elif case_name == "reference-value-mismatch":
        structural["value"] = 4
    elif case_name == "reference-vector-mismatch":
        rho["dimension"]["exponents"] = [-9, 0, 1]
    elif case_name == "nonstructural-binding":
        rho["dimension"]["law"]["bindings"]["D"] = "Q.medium.hbar"
    elif case_name == "wrong-law":
        rho["dimension"]["law"]["components"][0] = [
            "Sub",
            3,
            ["Mul", 2, ["Ref", "D"]],
        ]
    elif case_name == "wrong-binding":
        rho["dimension"]["law"]["bindings"]["D"] = "Q.medium.n_eos"
    elif case_name == "unbound":
        rho["dimension"]["law"]["components"][0] = [
            "Neg",
            ["Ref", "missing_D"],
        ]
    elif case_name == "unresolvable":
        rho["dimension"]["law"]["bindings"]["D"] = "Q.missing.D"


def demonstrate(case_name: str) -> int:
    quantities, relations, schema = load_raw_documents()
    baseline_registry = Registry.from_documents(
        quantities,
        relations,
        schema,
        ledger_root=Path(__file__).resolve().parents[3],
    )
    baseline_rho = baseline_registry.dimension_law("rho_br")
    baseline_reference = (
        baseline_rho.evaluate_reference() if baseline_rho is not None else None
    )
    _mutate(case_name, quantities)
    print(f"CASE: {case_name}")

    if case_name in {"absent-laws", "wrong-D-coefficient"}:
        registry = Registry.from_documents(
            quantities,
            relations,
            schema,
            ledger_root=Path(__file__).resolve().parents[3],
        )
        forms = inspect_declaration_forms(registry)
        coverage = inspect_symbolic_witness_coverage(registry)
        check_passed = print_dimension_law_check(registry, forms, coverage)
        gate_statuses = {
            result.relation_id: result.status
            for result in DimensionAudit(quantities, relations).classify_relations()
        }
        print(
            "RELATION_GATE_OPERAND expected-blind-population: "
            f"{[HOMOGENEOUS] * len(gate_statuses)}"
        )
        print(
            "RELATION_GATE_OPERAND observed-population: "
            f"{list(gate_statuses.values())}"
        )
        gate_residual = sum(
            status != HOMOGENEOUS for status in gate_statuses.values()
        )
        print(f"RELATION_GATE_RESIDUAL nonhomogeneous-count: {gate_residual}")
        if case_name == "wrong-D-coefficient":
            caught = (
                forms.passed
                and not coverage.passed
                and not check_passed
                and gate_residual == 0
            )
            marker = "CAUGHT" if caught else "ESCAPED"
            print(f"DIMENSION_LAW_ABLE_TO_FAIL_{marker}: {case_name}")
            return EXPECTED_FAILURE_EXIT if caught else 0
        caught = not forms.passed and not check_passed and gate_residual == 0
    else:
        if case_name in {
            "missing-binding-value",
            "noninteger-binding-value",
            "reference-value-mismatch",
        }:
            structural = _quantity(quantities, "Q.brane.D_brane")
            declared = structural.get("value", "<missing>")
            if baseline_rho is None:
                raise RuntimeError("rho_br has no baseline dimension law")
            references_by_binding_qid = {
                binding_qid: dict(baseline_rho.reference_values)[parameter]
                for parameter, binding_qid in baseline_rho.bindings
            }
            reference = references_by_binding_qid["Q.brane.D_brane"]
            if isinstance(declared, int) and not isinstance(declared, bool):
                reference_residual: object = declared - reference
            else:
                reference_residual = "UNDETERMINED"
            print(f"REFERENCE_OPERAND law-reference: {reference}")
            print(f"REFERENCE_OPERAND bound-declared-value: {declared!r}")
            print(f"REFERENCE_RESIDUAL declared-minus-reference: {reference_residual}")
        elif case_name == "reference-vector-mismatch":
            declared_vector = tuple(
                _quantity(quantities, "Q.brane.rho_br")["dimension"]["exponents"]
            )
            vector_residual = tuple(
                left - right
                for left, right in zip(declared_vector, baseline_reference or ())
            )
            print(f"REFERENCE_VECTOR_OPERAND law-evaluation: {baseline_reference}")
            print(f"REFERENCE_VECTOR_OPERAND retained-triple: {declared_vector}")
            print(f"REFERENCE_VECTOR_RESIDUAL retained-minus-law: {vector_residual}")
        elif case_name == "nonstructural-binding":
            target = _quantity(quantities, "Q.medium.hbar")
            expected_target = ("discrete-choice", "discrete-structural")
            declared_target = (target["kind"], target["counting_axis"])
            target_residual = tuple(
                left == right
                for left, right in zip(declared_target, expected_target)
            )
            print(f"STRUCTURAL_TARGET_OPERAND expected: {expected_target}")
            print(f"STRUCTURAL_TARGET_OPERAND declared: {declared_target}")
            print(f"STRUCTURAL_TARGET_RESIDUAL field-equalities: {target_residual}")

        audit = DimensionAudit(quantities, relations)
        statuses, _ = _print_relation_probe(audit, relations)
        loader_finding = _loader_finding(quantities, relations, schema)
        print(
            "REGISTRY_LOADER: "
            + (
                f"CAUGHT {loader_finding}"
                if loader_finding is not None
                else "ESCAPED"
            )
        )
        expected_status = (
            INHOMOGENEOUS if case_name == "wrong-law" else UNDETERMINED
        )
        caught = (
            statuses.get("R4") == expected_status
            and statuses.get("R5") == expected_status
            and loader_finding is not None
        )

    marker = "CAUGHT" if caught else "ESCAPED"
    print(f"DIMENSION_LAW_ABLE_TO_FAIL_{marker}: {case_name}")
    return EXPECTED_FAILURE_EXIT if caught else 0


def demonstrate_escape() -> int:
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


def classify_child(
    case_name: str, completed: subprocess.CompletedProcess[str]
) -> str:
    if completed.returncode == 0:
        return "ESCAPED"
    if case_name not in EXPECTED_CASE_STATUS:
        return "ERROR"
    expected_status = EXPECTED_CASE_STATUS[case_name]
    expected_marker = f"DIMENSION_LAW_ABLE_TO_FAIL_{expected_status}: {case_name}"
    if (
        completed.returncode == EXPECTED_FAILURE_EXIT
        and expected_marker in completed.stdout.splitlines()
        and not completed.stderr.strip()
    ):
        return expected_status
    return "ERROR"


def _echo_child(case_name: str, stream_name: str, output: str) -> None:
    for line in output.splitlines():
        safe = line.replace(PARENT_VERDICT_TOKEN, "DIMENSION_LAW_CHILD_TEXT:")
        print(f"CHILD {case_name} {stream_name}: {safe}")


def run_harness(
    case_names: tuple[str, ...] = CASES, *, require_full_population: bool = True
) -> int:
    statuses: list[str] = []
    population_complete = not require_full_population or case_names == CASES
    print(f"HARNESS_POPULATION_OPERAND expected: {list(CASES)}")
    print(f"HARNESS_POPULATION_OPERAND observed: {list(case_names)}")
    population_residual = len(case_names) - len(CASES)
    print(f"HARNESS_POPULATION_RESIDUAL count: {population_residual}")
    for case_name in case_names:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--case", case_name],
            check=False,
            capture_output=True,
            text=True,
            timeout=60,
        )
        _echo_child(case_name, "stdout", completed.stdout)
        _echo_child(case_name, "stderr", completed.stderr)
        status = classify_child(case_name, completed)
        statuses.append(status)
        print(
            f"DEMONSTRATION {case_name}: status={status} "
            f"observed_exit={completed.returncode}"
        )
    if not population_complete or "ERROR" in statuses:
        result = "ERROR"
    elif all(
        status == EXPECTED_CASE_STATUS.get(case_name, "CAUGHT")
        for case_name, status in zip(case_names, statuses)
    ):
        result = "PASS"
    else:
        result = "FAIL"
    print(f"{PARENT_VERDICT_TOKEN} {result}")
    return 0 if result == "PASS" else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=CASES + (ESCAPE_CASE,))
    parser.add_argument("--demonstrate-escape", action="store_true")
    parser.add_argument("--weaken-case-exit", action="store_true", help=argparse.SUPPRESS)
    arguments = parser.parse_args()
    if arguments.case == ESCAPE_CASE:
        return demonstrate_escape()
    if arguments.case:
        observed = demonstrate(arguments.case)
        return 0 if arguments.weaken_case_exit else observed
    if arguments.demonstrate_escape:
        return run_harness((ESCAPE_CASE,), require_full_population=False)
    return run_harness()


if __name__ == "__main__":
    raise SystemExit(main())
