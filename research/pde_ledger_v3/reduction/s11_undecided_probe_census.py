#!/usr/bin/env python3
"""Kernel-free undecided-record probe census over one S11 record file.

Every raw ``UNDECIDED`` substring is reconciled.  Exact verdict tokens are
probed from their own emitted operands; strings mirrored into export metadata
remain visible as non-verdict text.  Exit status reports execution only.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
from pathlib import Path
import time
from typing import Any, Mapping

import sympy as sp

from s11_census_common import (
    EQUATIONS_SUFFIX,
    enumerate_probe_population,
    quote,
    read_record_map,
    undecided_occurrence_classes,
)
from s11_census_math import (
    OPERATION_MEMORY_BUDGET_BYTES,
    OPERATION_TIME_BUDGET_SECONDS,
    association,
    conjunction,
    enforce_memory_budget,
    equation_residuals,
    parse_payload,
    quoted,
    render_object,
    render_text,
    sample_counterexample,
    sample_existence,
    sequence,
    simplify_residual,
    split_delimited,
    split_top_relation,
    string_value,
)


def population_line(record: Path) -> str:
    population = enumerate_probe_population(record)
    verdict = (
        "POPULATION_RECONCILED"
        if population.unexplained_gap == 0
        else "UNRECONCILED_UNDECIDED"
    )
    operands = (
        f"record={quote(population.record)} dialect={population.dialect} "
        f"raw_undecided={population.raw_undecided} "
        f"raw_undecided_lines={population.raw_undecided_lines} "
        f"status_token_occurrences={population.status_token_occurrences} "
        f"status_token_lines={population.status_token_lines} "
        f"sign_token_occurrences={population.sign_token_occurrences} "
        f"sign_token_lines={population.sign_token_lines} "
        f"bare_exact_occurrences={population.bare_exact_occurrences} "
        f"embedded_text_occurrences={population.embedded_text_occurrences} "
        f"unexplained_gap={population.unexplained_gap}"
    )
    return f"PROBE_POPULATION {operands} verdict={verdict}"


def _wl_field_text(payload: str, key: str) -> str:
    for item in split_delimited(payload):
        left, right = split_top_relation(item, "->")
        if left == quoted(key):
            return right
    raise KeyError(key)


def _mapping_branch(value: object) -> dict[sp.Expr, sp.Expr]:
    current = value
    while isinstance(current, (tuple, list, sp.Tuple)) and len(current) == 1:
        only = current[0]
        if isinstance(only, sp.Basic) and only.func.__name__ == "Rule":
            break
        if isinstance(only, (tuple, list, sp.Tuple)):
            current = only
        else:
            break
    result: dict[sp.Expr, sp.Expr] = {}
    for entry in sequence(current):
        if isinstance(entry, sp.Basic) and entry.func.__name__ == "Rule":
            result[entry.args[0]] = entry.args[1]
        elif isinstance(entry, (tuple, list, sp.Tuple)) and len(entry) == 2:
            result[entry[0]] = entry[1]
        else:
            raise ValueError(f"unrecognized branch entry: {entry!r}")
    return result


def _residual_match(left: tuple[sp.Expr, ...], right: tuple[sp.Expr, ...]) -> bool:
    if len(left) != len(right):
        return False
    return all(simplify_residual(a - b)[1] in {"ZERO", "ZERO_SAMPLED"} for a, b in zip(left, right))


def _outcome(
    operand: object,
    recomputed: object,
    decision: object,
    verdict: str,
    **extra: object,
) -> dict[str, object]:
    return {
        "operand": render_object(operand),
        "recomputed": render_object(recomputed),
        "decision": render_object(decision),
        "verdict": verdict,
        **extra,
    }


def _identity_probe(dialect: str, payload: str, equation_payload: str | None) -> list[dict[str, object]]:
    parsed = parse_payload(dialect, payload)
    fields = association(parsed)
    operands = association(fields.get("OPERANDS"))
    emitted = operands.get("RESIDUALS")
    if emitted is None:
        return [_outcome(fields, (), (), "MISSING_OPERANDS")]
    emitted_residuals = tuple(sequence(emitted))
    if equation_payload is None:
        return [_outcome(emitted, (), (), "MISSING_OPERANDS")]
    equations = equation_residuals(parse_payload(dialect, equation_payload))
    match = _residual_match(tuple(emitted_residuals), equations)
    statuses = tuple(simplify_residual(item)[1] for item in emitted_residuals)
    if not match:
        verdict = "RESIDUAL_MISMATCH"
    elif all(status in {"ZERO", "ZERO_SAMPLED"} for status in statuses):
        verdict = "DECIDED_UNDECIDED_RECORD"
    elif any(status in {"NONZERO", "NONZERO_SAMPLED", "UNDEFINED"} for status in statuses):
        verdict = "DECIDED_UNDECIDED_RECORD"
    else:
        verdict = "STILL_UNDECIDED"
    return [
        _outcome(
            emitted,
            equations,
            statuses,
            verdict,
            residual_match=match,
        )
    ]


def _sign_probe(dialect: str, payload: str) -> list[dict[str, object]]:
    fields = association(parse_payload(dialect, payload))
    operand = fields.get("OPERAND")
    if operand is None:
        operands = association(fields.get("OPERANDS"))
        operand = operands.get("OPERAND")
    if not isinstance(operand, sp.Basic):
        return [_outcome(fields, (), (), "MISSING_OPERANDS")]
    simplified = sp.simplify(operand)
    signs = {
        "positive": simplified.is_positive,
        "negative": simplified.is_negative,
        "zero": simplified.is_zero,
    }
    verdict = (
        "DECIDED_UNDECIDED_RECORD"
        if any(value is True for value in signs.values())
        else "STILL_UNDECIDED"
    )
    return [_outcome(operand, simplified, signs, verdict)]


def _admissible_probe(dialect: str, payload: str) -> list[dict[str, object]]:
    parsed = parse_payload(dialect, payload)
    outcomes: list[dict[str, object]] = []
    for entry in sequence(parsed):
        fields = association(entry)
        if string_value(fields.get("STATUS_TOKEN")) != "UNDECIDED":
            continue
        operands_value = fields.get("OPERANDS")
        operands = association(operands_value)
        branch_value = operands.get("BRANCH", fields.get("BRANCH"))
        premises = operands.get("PREMISES")
        if dialect == "PY" and not operands:
            values = sequence(operands_value)
            if len(values) >= 2:
                branch_value, premises = values[0], values[1]
        if branch_value is None or premises is None:
            outcomes.append(_outcome(fields, (), (), "MISSING_OPERANDS"))
            continue
        branch = _mapping_branch(branch_value)
        branch_relations = tuple(
            sp.Eq(lhs, rhs, evaluate=False) for lhs, rhs in branch.items()
        )
        condition = conjunction((*sequence(premises), *branch_relations))
        decision = sample_existence(condition)
        verdict = (
            "DECIDED_UNDECIDED_RECORD"
            if decision["truth"] in {"TRUE", "FALSE"}
            else "STILL_UNDECIDED"
        )
        outcomes.append(_outcome((branch_value, premises), condition, decision, verdict))
    return outcomes


def _operand_equations_premises(dialect: str, payload: str) -> tuple[object, object, object]:
    parsed = parse_payload(dialect, payload)
    fields = association(parsed)
    if fields:
        return fields.get("EQUATIONS"), fields.get("SOLVE_VARIABLES"), fields.get("PREMISES")
    values = sequence(parsed)
    if len(values) < 3:
        return None, None, None
    return values[0], values[1], values[2]


def _real_status_probe(
    dialect: str,
    operand_payload: str | None,
    equation_payload: str | None,
) -> list[dict[str, object]]:
    if operand_payload is None or equation_payload is None:
        return [_outcome((), (), (), "MISSING_OPERANDS")]
    emitted_equations, _variables, premises = _operand_equations_premises(dialect, operand_payload)
    if emitted_equations is None or premises is None:
        return [_outcome((), (), (), "MISSING_OPERANDS")]
    operand_residuals = equation_residuals(emitted_equations)
    record_residuals = equation_residuals(parse_payload(dialect, equation_payload))
    match = _residual_match(operand_residuals, record_residuals)
    conditions = tuple(sp.Eq(item, 0, evaluate=False) for item in operand_residuals)
    condition = conjunction((*sequence(premises), *conditions))
    decision = sample_existence(condition)
    if not match:
        verdict = "RESIDUAL_MISMATCH"
    elif decision["truth"] in {"TRUE", "FALSE"}:
        verdict = "DECIDED_UNDECIDED_RECORD"
    else:
        verdict = "STILL_UNDECIDED"
    return [_outcome((emitted_equations, premises), condition, decision, verdict, residual_match=match)]


def _inconsistent_probe(dialect: str, payload: str) -> list[dict[str, object]]:
    fields = association(parse_payload(dialect, payload))
    operands_value = fields.get("OPERANDS")
    operands = association(operands_value)
    if operands:
        equations = operands.get("EQUATIONS")
        variables = operands.get("SOLVE_VARIABLES")
    else:
        values = sequence(operands_value)
        equations = values[0] if len(values) > 0 else None
        variables = values[1] if len(values) > 1 else None
    if equations is None or variables is None:
        return [_outcome(fields, (), (), "MISSING_OPERANDS")]
    residuals = equation_residuals(equations)
    try:
        solved = sp.nonlinsolve(residuals, tuple(sequence(variables)))
        if solved is sp.EmptySet:
            decision = "PROVED_INCONSISTENT"
            verdict = "DECIDED_UNDECIDED_RECORD"
        elif isinstance(solved, sp.FiniteSet) and len(solved) > 0:
            decision = "PROVED_CONSISTENT"
            verdict = "DECIDED_UNDECIDED_RECORD"
        else:
            decision = solved
            verdict = "STILL_UNDECIDED"
    except BaseException as exc:
        decision = f"{type(exc).__name__}:{exc}"
        verdict = "STILL_UNDECIDED"
    return [_outcome((equations, variables), residuals, decision, verdict)]


def _count_probe(dialect: str, payload: str) -> list[dict[str, object]]:
    if dialect == "WL":
        decision_text = _wl_field_text(payload, "DECISION_OPERANDS")
        test_text = _wl_field_text(decision_text, "TEST_OBJECT")
        premises_text = _wl_field_text(decision_text, "PREMISES")
        if test_text.startswith("Failure["):
            return [_outcome((test_text, premises_text), (), "EMITTED_FAILURE", "STILL_UNDECIDED")]
        test_object = parse_payload(dialect, test_text)
        premises = parse_payload(dialect, premises_text)
    else:
        fields = association(parse_payload(dialect, payload))
        decision_fields = association(fields.get("DECISION_OPERANDS"))
        test_object = decision_fields.get("TEST_OBJECT")
        premises = decision_fields.get("PREMISES")
        if test_object is None or premises is None:
            return [_outcome(fields, (), (), "MISSING_OPERANDS")]
    predicate = (
        test_object.args[-1]
        if isinstance(test_object, sp.Basic) and test_object.func.__name__ == "ForAll"
        else test_object
    )
    decision = sample_counterexample(predicate, premises)
    verdict = (
        "DECIDED_UNDECIDED_RECORD"
        if decision["truth"] in {"TRUE", "FALSE"}
        else "STILL_UNDECIDED"
    )
    return [_outcome((test_object, premises), predicate, decision, verdict)]


def _probe_worker(
    connection: Any,
    dialect: str,
    tag: str,
    payload: str,
    operand_payload: str | None,
    equation_payload: str | None,
) -> None:
    enforce_memory_budget()
    try:
        if tag.endswith("_IDENTICALLY_SATISFIED"):
            outcomes = _identity_probe(dialect, payload, equation_payload)
        elif tag.endswith("_SIGN"):
            outcomes = _sign_probe(dialect, payload)
        elif tag.endswith("_REAL_ADMISSIBLE"):
            outcomes = _admissible_probe(dialect, payload)
        elif tag.endswith("_REAL_STATUS"):
            outcomes = _real_status_probe(dialect, operand_payload, equation_payload)
        elif tag.endswith("_INCONSISTENT"):
            outcomes = _inconsistent_probe(dialect, payload)
        elif "_STRATUM" in tag:
            outcomes = _count_probe(dialect, payload)
        else:
            outcomes = [_outcome(payload, (), (), "MISSING_OPERANDS")]
        connection.send({"kind": "PROBE_RESULT", "outcomes": outcomes})
    except (KeyError, ValueError, TypeError, SyntaxError, sp.SympifyError) as exc:
        connection.send(
            {
                "kind": "OPERAND_PARSE_FAILURE",
                "detail": f"{type(exc).__name__}:{exc}",
                "failing_text": payload if operand_payload is None else operand_payload,
            }
        )
    except BaseException as exc:
        connection.send(
            {
                "kind": "PROBE_ERROR",
                "detail": f"{type(exc).__name__}:{exc}",
            }
        )
    finally:
        connection.close()


def _run_bounded(args: tuple[object, ...]) -> tuple[dict[str, object], float]:
    context = mp.get_context("spawn")
    parent, child = context.Pipe(duplex=False)
    process = context.Process(target=_probe_worker, args=(child, *args))
    started = time.monotonic()
    process.start()
    child.close()
    if parent.poll(OPERATION_TIME_BUDGET_SECONDS):
        try:
            result = parent.recv()
        except EOFError:
            result = {"kind": "WORKER_FAILURE"}
        process.join(5)
    else:
        process.kill()
        process.join()
        result = {"kind": "TIME_EXPIRED"}
    elapsed = time.monotonic() - started
    if process.is_alive():
        process.kill()
        process.join()
    parent.close()
    return result, elapsed


def _sibling_payloads(tag: str, data: Mapping[str, Any]) -> tuple[str | None, str | None]:
    operand = None
    equation = None
    if tag.endswith("_REAL_STATUS"):
        operand_item = data.get(tag + "_OPERANDS")
        stem = tag[: -len("_REAL_STATUS")]
        equation_item = data.get(stem + EQUATIONS_SUFFIX)
        operand = operand_item.payload if operand_item else None
        equation = equation_item.payload if equation_item else None
    elif tag.endswith("_IDENTICALLY_SATISFIED"):
        stem = tag[: -len("_IDENTICALLY_SATISFIED")]
        equation_item = data.get(stem + EQUATIONS_SUFFIX)
        equation = equation_item.payload if equation_item else None
    return operand, equation


def run_census(record: Path) -> None:
    population = enumerate_probe_population(record)
    print(population_line(record), flush=True)
    print(
        "PROBE_BUDGET "
        f"record={quote(record)} operation_seconds={OPERATION_TIME_BUDGET_SECONDS} "
        f"operation_memory_bytes={OPERATION_MEMORY_BUDGET_BYTES} "
        "solver_route=OPERAND_RECOMPUTE_NONLINSOLVE_EXACT_SAMPLE "
        "verdict=BUDGET_DECLARED",
        flush=True,
    )
    data = read_record_map(record)
    counts: dict[str, int] = {
        "occurrences": 0,
        "decided": 0,
        "still_undecided": 0,
        "missing_operands": 0,
        "unparseable_operands": 0,
        "residual_mismatches": 0,
        "non_verdict_text": 0,
        "resource_expiries": 0,
    }
    class_counts: dict[str, int] = {}
    expiry_counts: dict[str, int] = {}

    for item in sorted(data.values(), key=lambda value: value.line_no):
        if "UNDECIDED" not in item.payload:
            continue
        classes = undecided_occurrence_classes(population.dialect, item.payload)
        counts["occurrences"] += len(classes)
        for classification in classes:
            class_counts[classification] = class_counts.get(classification, 0) + 1
        if "_EXPORT_CANDIDATE_KEY_OPERANDS" in item.tag:
            for occurrence_index, classification in enumerate(classes, 1):
                counts["non_verdict_text"] += 1
                print(
                    "PROBE_RECORD "
                    f"record={quote(record)} tag={item.tag} line={item.line_no} "
                    f"occurrence={occurrence_index} class={classification} "
                    f"operand={render_text(item.payload)} recomputed=NOT_APPLICABLE "
                    "decision=NON_VERDICT_MIRROR verdict=NON_VERDICT_TEXT",
                    flush=True,
                )
            continue

        semantic_tag = item.tag
        semantic_payload = item.payload
        if item.tag.endswith("_STATUS") and not item.tag.endswith("_REAL_STATUS"):
            candidate_tag = item.tag[: -len("_STATUS")]
            candidate = data.get(candidate_tag)
            if candidate is not None:
                semantic_tag = candidate_tag
                semantic_payload = candidate.payload
        operand_payload, equation_payload = _sibling_payloads(semantic_tag, data)
        print(
            "PROBE_BEGIN "
            f"record={quote(record)} tag={item.tag} semantic_tag={semantic_tag} "
            f"line={item.line_no} occurrence_classes={list(classes)} "
            f"payload_operand={render_text(semantic_payload)} state=STARTED",
            flush=True,
        )
        result, elapsed = _run_bounded(
            (
                population.dialect,
                semantic_tag,
                semantic_payload,
                operand_payload,
                equation_payload,
            )
        )
        if result["kind"] in {"TIME_EXPIRED", "WORKER_FAILURE", "PROBE_ERROR"}:
            expiry_counts[semantic_tag.rsplit("_", 1)[-1]] = (
                expiry_counts.get(semantic_tag.rsplit("_", 1)[-1], 0) + len(classes)
            )
            counts["resource_expiries"] += len(classes)
            for occurrence_index, classification in enumerate(classes, 1):
                print(
                    "PROBE_RECORD "
                    f"record={quote(record)} tag={item.tag} line={item.line_no} "
                    f"occurrence={occurrence_index} class={classification} "
                    f"operand={render_text(semantic_payload)} recomputed=NOT_COMPUTED "
                    f"decision={result['kind']} elapsed_seconds={elapsed:.6f} "
                    "verdict=PROBE_RESOURCE_EXPIRED",
                    flush=True,
                )
            continue
        if result["kind"] == "OPERAND_PARSE_FAILURE":
            counts["unparseable_operands"] += len(classes)
            for occurrence_index, classification in enumerate(classes, 1):
                print(
                    "PROBE_RECORD "
                    f"record={quote(record)} tag={item.tag} line={item.line_no} "
                    f"occurrence={occurrence_index} class={classification} "
                    f"operand={quoted(result['failing_text'])} recomputed=NOT_COMPUTED "
                    f"decision={quoted(result['detail'])} verdict=UNPARSEABLE_OPERANDS",
                    flush=True,
                )
            continue

        outcomes = list(result["outcomes"])  # type: ignore[index]
        if len(outcomes) == 1 and len(classes) > 1:
            outcomes *= len(classes)
        if len(outcomes) != len(classes):
            outcomes = [
                _outcome(
                    semantic_payload,
                    (),
                    f"outcomes={len(outcomes)} occurrences={len(classes)}",
                    "MISSING_OPERANDS",
                )
                for _ in classes
            ]
        for occurrence_index, (classification, outcome) in enumerate(zip(classes, outcomes), 1):
            verdict = str(outcome["verdict"])
            if verdict == "DECIDED_UNDECIDED_RECORD":
                counts["decided"] += 1
            elif verdict == "STILL_UNDECIDED":
                counts["still_undecided"] += 1
            elif verdict == "MISSING_OPERANDS":
                counts["missing_operands"] += 1
            elif verdict == "RESIDUAL_MISMATCH":
                counts["residual_mismatches"] += 1
            print(
                "PROBE_RECORD "
                f"record={quote(record)} tag={item.tag} line={item.line_no} "
                f"occurrence={occurrence_index} class={classification} "
                f"operand={outcome['operand']} recomputed={outcome['recomputed']} "
                f"decision={outcome['decision']} verdict={verdict}",
                flush=True,
            )

    print(
        "PROBE_SUMMARY "
        f"record={quote(record)} "
        + " ".join(f"{key}={value}" for key, value in counts.items())
        + f" class_counts={class_counts} expiry_counts={expiry_counts} "
        "verdict=CENSUS_EXECUTED",
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("record", type=Path)
    parser.add_argument("--enumerate-only", action="store_true")
    args = parser.parse_args()
    record = args.record.resolve(strict=True)
    if args.enumerate_only:
        print(population_line(record), flush=True)
    else:
        run_census(record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
