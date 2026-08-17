#!/usr/bin/env python3
"""Independent containment/soundness census over one S11 record file.

Exit status reports execution of the instrument, never the physics findings.
The named acceptance reducer owns the round decision.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
from pathlib import Path
import time
from typing import Any

import sympy as sp

from s11_census_common import (
    EQUATIONS_SUFFIX,
    SOLUTION_SUFFIX,
    STATUS_SUFFIX,
    WITNESS_SUFFIX,
    enumerate_containment_population,
    proved_nonempty_payload,
    quote,
    read_record_map,
    sentinel_payload,
)
from s11_census_math import (
    OPERATION_MEMORY_BUDGET_BYTES,
    OPERATION_TIME_BUDGET_SECONDS,
    association,
    branch_memberships,
    completeness_candidates,
    enforce_memory_budget,
    equation_residuals,
    parse_payload,
    parse_solution,
    quoted,
    render_object,
    render_text,
    sequence,
    simplify_residual,
)


def population_line(record: Path) -> str:
    population = enumerate_containment_population(record)
    tags = "[" + ",".join(quote(tag) for tag in population.unreconciled_tags) + "]"
    verdict = (
        "POPULATION_RECONCILED"
        if not population.unreconciled_tags
        else "UNRECONCILED_LINE"
    )
    operands = (
        f"record={quote(population.record)} dialect={population.dialect} "
        f"raw_solution={population.raw_solution} "
        f"raw_equations={population.raw_equations} "
        f"raw_real_witness={population.raw_real_witness} "
        f"protocol_pairs={population.protocol_pairs} "
        f"semantic_pairs={population.semantic_pairs} "
        f"paired_not_applicable={population.paired_not_applicable} "
        f"witness_proved_nonempty={population.witness_proved_nonempty} "
        f"witness_not_applicable={population.witness_not_applicable} "
        f"excluded_dim_pairs={population.excluded_dim_pairs} "
        f"excluded_c1_equations={population.excluded_c1_equations} "
        f"excluded_stratum_defining_equations="
        f"{population.excluded_stratum_defining_equations} "
        f"unreconciled_count={len(population.unreconciled_tags)} "
        f"unreconciled_tags={tags}"
    )
    return f"CONTAINMENT_POPULATION {operands} verdict={verdict}"


def _branch_object(branch: dict[sp.Expr, sp.Expr]) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    return tuple(sorted(branch.items(), key=lambda item: sp.default_sort_key(item[0])))


def _py_variables(tag: str, parsed: Any, residuals: tuple[sp.Expr, ...]) -> tuple[sp.Symbol, ...]:
    lhs = list(parsed.variables)
    free = sorted(set().union(*(item.free_symbols for item in residuals)), key=sp.default_sort_key)
    marker = tag[: -len(SOLUTION_SUFFIX)]
    if marker.endswith("_K"):
        inferred = [symbol for symbol in free if str(symbol).startswith("k") and str(symbol)[1:].isdigit()]
    elif marker.endswith("_COEFF"):
        inferred = [symbol for symbol in free if not (str(symbol).startswith("k") and str(symbol)[1:].isdigit())]
    elif marker.endswith("_JOINT"):
        inferred = free
    else:
        inferred = lhs
    for symbol in inferred:
        if symbol not in lhs:
            lhs.append(symbol)
    return tuple(lhs)


def _pair_worker(
    connection: Any,
    dialect: str,
    tag: str,
    solution_payload: str,
    equation_payload: str,
) -> None:
    enforce_memory_budget()
    try:
        parsed = parse_solution(dialect, solution_payload)
    except BaseException as exc:
        connection.send(
            {
                "kind": "PARSE_FAILURE",
                "stage": "SOLUTION",
                "detail": f"{type(exc).__name__}:{exc}",
                "failing_text": solution_payload,
            }
        )
        connection.close()
        return
    try:
        equation_object = parse_payload(dialect, equation_payload)
        residuals = equation_residuals(equation_object)
    except BaseException as exc:
        connection.send(
            {
                "kind": "PARSE_FAILURE",
                "stage": "EQUATIONS",
                "detail": f"{type(exc).__name__}:{exc}",
                "failing_text": equation_payload,
            }
        )
        connection.close()
        return

    variables = parsed.variables
    if dialect == "PY":
        variables = _py_variables(tag, parsed, residuals)
    branch_results: list[dict[str, object]] = []
    for branch_index, branch in enumerate(parsed.branches, 1):
        rows, tested, total = branch_memberships(residuals, branch)
        row_verdicts = [str(row["verdict"]) for row in rows]
        if "SPURIOUS_BRANCH" in row_verdicts:
            verdict = "SPURIOUS_BRANCH"
        elif "SPURIOUS_BRANCH_SAMPLED" in row_verdicts:
            verdict = "SPURIOUS_BRANCH_SAMPLED"
        elif tested < total:
            verdict = "SHEET_INCOMPLETE"
        elif "BRANCH_MEMBERSHIP_UNDECIDED" in row_verdicts:
            verdict = "BRANCH_MEMBERSHIP_UNDECIDED"
        elif "BRANCH_CONTAINED_SAMPLED" in row_verdicts:
            verdict = "BRANCH_CONTAINED_SAMPLED"
        else:
            verdict = "BRANCH_CONTAINED"
        branch_results.append(
            {
                "index": branch_index,
                "branch": render_object(_branch_object(branch)),
                "rows": [
                    {
                        "sheet": row["sheet"],
                        "residuals": render_object(row["residuals"]),
                        "definedness": row["definedness"],
                        "verdict": row["verdict"],
                    }
                    for row in rows
                ],
                "tested": tested,
                "total": total,
                "verdict": verdict,
            }
        )

    try:
        completeness = completeness_candidates(residuals, variables, parsed.branches)
        completeness_result = {
            **completeness,
            "missing": [
                render_object(_branch_object(branch))
                for branch in completeness["missing"]  # type: ignore[index]
            ],
        }
    except BaseException as exc:
        completeness_result = {
            "factor_covers": 0,
            "cover_truncated": False,
            "candidate_count": 0,
            "unresolved_count": 1,
            "missing": [],
            "detail": f"{type(exc).__name__}:{exc}",
            "verdict": "COMPLETENESS_UNDECIDED",
        }
    connection.send(
        {
            "kind": "PAIR_RESULT",
            "solution_kind": parsed.kind,
            "variables": render_object(variables),
            "equation_count": len(residuals),
            "branches": branch_results,
            "completeness": completeness_result,
        }
    )
    connection.close()


def _premises_from_operand(dialect: str, payload: str) -> tuple[object, object]:
    parsed = parse_payload(dialect, payload)
    if dialect == "WL":
        fields = association(parsed)
        return fields.get("EQUATIONS"), fields.get("PREMISES")
    values = sequence(parsed)
    if len(values) < 3:
        return None, None
    return values[0], values[2]


def _residual_lists_equivalent(left: tuple[sp.Expr, ...], right: tuple[sp.Expr, ...]) -> bool:
    if len(left) != len(right):
        return False
    return all(simplify_residual(a - b)[1] == "ZERO" for a, b in zip(left, right))


def _witness_worker(
    connection: Any,
    dialect: str,
    witness_payload: str,
    equation_payload: str,
    operand_payload: str,
) -> None:
    enforce_memory_budget()
    stages = (
        ("WITNESS", witness_payload),
        ("EQUATIONS", equation_payload),
        ("REAL_STATUS_OPERANDS", operand_payload),
    )
    try:
        witness_object = parse_payload(dialect, witness_payload)
        if dialect == "WL":
            branch: dict[sp.Expr, sp.Expr] = {}
            for entry in sequence(witness_object):
                if not (isinstance(entry, sp.Basic) and entry.func.__name__ == "Rule"):
                    raise ValueError(f"unrecognized witness rule {entry!r}")
                branch[entry.args[0]] = entry.args[1]
        else:
            branch = {entry[0]: entry[1] for entry in sequence(witness_object)}  # type: ignore[index]
        equation_object = parse_payload(dialect, equation_payload)
        residuals = equation_residuals(equation_object)
        operand_equations, premises = _premises_from_operand(dialect, operand_payload)
        if operand_equations is None or premises is None:
            raise ValueError("missing EQUATIONS or PREMISES in REAL_STATUS_OPERANDS")
        operand_residuals = equation_residuals(operand_equations)
    except BaseException as exc:
        stage, failing = next(
            ((name, text) for name, text in stages if text in str(exc)),
            ("WITNESS_OR_OPERANDS", witness_payload + "\n" + equation_payload + "\n" + operand_payload),
        )
        connection.send(
            {
                "kind": "PARSE_FAILURE",
                "stage": stage,
                "detail": f"{type(exc).__name__}:{exc}",
                "failing_text": failing,
            }
        )
        connection.close()
        return

    residual_match = _residual_lists_equivalent(residuals, operand_residuals)
    rows, tested, total = branch_memberships(residuals, branch)
    substituted_premises = premises.subs(branch, simultaneous=True)  # type: ignore[union-attr]
    try:
        premise_truth = str(sp.simplify(substituted_premises))
    except BaseException:
        premise_truth = str(substituted_premises)
    membership_verdicts = [str(row["verdict"]) for row in rows]
    if not residual_match:
        verdict = "RESIDUAL_MISMATCH"
    elif any(item.startswith("SPURIOUS_BRANCH") for item in membership_verdicts):
        verdict = "WITNESS_FAILURE"
    elif tested < total:
        verdict = "WITNESS_SHEET_INCOMPLETE"
    elif "BRANCH_MEMBERSHIP_UNDECIDED" in membership_verdicts:
        verdict = "WITNESS_UNDECIDED"
    elif "BRANCH_CONTAINED_SAMPLED" in membership_verdicts:
        verdict = "WITNESS_VALIDATED_SAMPLED"
    elif substituted_premises is sp.false:
        verdict = "WITNESS_FAILURE"
    else:
        verdict = "WITNESS_VALIDATED"
    connection.send(
        {
            "kind": "WITNESS_RESULT",
            "witness": render_object(_branch_object(branch)),
            "residual_match": residual_match,
            "rows": [
                {
                    "sheet": row["sheet"],
                    "residuals": render_object(row["residuals"]),
                    "definedness": row["definedness"],
                    "verdict": row["verdict"],
                }
                for row in rows
            ],
            "tested": tested,
            "total": total,
            "substituted_premises": render_object(substituted_premises),
            "premise_truth": premise_truth,
            "verdict": verdict,
        }
    )
    connection.close()


def _run_bounded(worker: Any, args: tuple[object, ...]) -> tuple[dict[str, object], float]:
    context = mp.get_context("spawn")
    parent, child = context.Pipe(duplex=False)
    process = context.Process(target=worker, args=(child, *args))
    started = time.monotonic()
    process.start()
    child.close()
    if parent.poll(OPERATION_TIME_BUDGET_SECONDS):
        try:
            result = parent.recv()
        except EOFError:
            result = {"kind": "WORKER_FAILURE", "exitcode": process.exitcode}
        process.join(5)
    else:
        process.kill()
        process.join()
        result = {"kind": "TIME_EXPIRED", "exitcode": process.exitcode}
    elapsed = time.monotonic() - started
    if process.is_alive():
        process.kill()
        process.join()
    parent.close()
    return result, elapsed


def _print_parse_failure(record: Path, tag: str, result: dict[str, object]) -> None:
    print(
        "CONTAINMENT_PARSE "
        f"record={quote(record)} tag={tag} stage={result.get('stage')} "
        f"detail={quoted(result.get('detail'))} "
        f"failing_text={quoted(result.get('failing_text'))} "
        "membership=NOT_COMPUTED verdict=IN_POPULATION_PARSE_FAILURE",
        flush=True,
    )


def run_census(record: Path) -> None:
    population = enumerate_containment_population(record)
    print(population_line(record), flush=True)
    print(
        "CONTAINMENT_BUDGET "
        f"record={quote(record)} operation_seconds={OPERATION_TIME_BUDGET_SECONDS} "
        f"operation_memory_bytes={OPERATION_MEMORY_BUDGET_BYTES} "
        "solver_route=FACTOR_COVER_NONLINSOLVE verdict=BUDGET_DECLARED",
        flush=True,
    )
    data = read_record_map(record)
    sentinel = sentinel_payload(population.dialect)
    proved = proved_nonempty_payload(population.dialect)

    summary = {
        "semantic_pairs": 0,
        "paired_not_applicable": 0,
        "pair_parse_failures": 0,
        "branches": 0,
        "spurious_branches": 0,
        "omitted_records": 0,
        "sheet_incomplete_branches": 0,
        "pair_expiries": 0,
        "witnesses": 0,
        "witness_failures": 0,
        "witness_parse_failures": 0,
        "witness_expiries": 0,
        "residual_mismatches": 0,
    }

    solution_tags = sorted(tag for tag in data if tag.endswith(SOLUTION_SUFFIX))
    semantic_index = 0
    for tag in solution_tags:
        stem = tag[: -len(SOLUTION_SUFFIX)]
        equations_tag = stem + EQUATIONS_SUFFIX
        if stem.endswith("_DIM") or equations_tag not in data:
            continue
        solution = data[tag]
        equation = data[equations_tag]
        if solution.payload == sentinel and equation.payload == sentinel:
            summary["paired_not_applicable"] += 1
            print(
                "CONTAINMENT_SENTINEL "
                f"record={quote(record)} tag={tag} "
                f"solution_operand={quoted(solution.payload)} "
                f"equations_operand={quoted(equation.payload)} "
                "membership=EXACT_PAIRED_SENTINEL verdict=PAIRED_NOT_APPLICABLE",
                flush=True,
            )
            continue
        semantic_index += 1
        summary["semantic_pairs"] += 1
        print(
            "CONTAINMENT_BEGIN "
            f"record={quote(record)} index={semantic_index} tag={tag} "
            f"solution_line={solution.line_no} equations_line={equation.line_no} "
            f"solution_operand={render_text(solution.payload)} "
            f"equations_operand={render_text(equation.payload)} state=STARTED",
            flush=True,
        )
        result, elapsed = _run_bounded(
            _pair_worker,
            (population.dialect, tag, solution.payload, equation.payload),
        )
        if result["kind"] == "PARSE_FAILURE":
            summary["pair_parse_failures"] += 1
            _print_parse_failure(record, tag, result)
            continue
        if result["kind"] in {"TIME_EXPIRED", "WORKER_FAILURE"}:
            summary["pair_expiries"] += 1
            print(
                "CONTAINMENT_PAIR "
                f"record={quote(record)} tag={tag} elapsed_seconds={elapsed:.6f} "
                f"budget_seconds={OPERATION_TIME_BUDGET_SECONDS} "
                f"worker_state={result['kind']} membership=NOT_COMPUTED "
                "verdict=PAIR_RESOURCE_EXPIRED",
                flush=True,
            )
            continue
        for branch in result["branches"]:  # type: ignore[index]
            summary["branches"] += 1
            for row in branch["rows"]:  # type: ignore[index]
                print(
                    "CONTAINMENT_BRANCH_SHEET "
                    f"record={quote(record)} tag={tag} branch_index={branch['index']} "
                    f"branch_operand={branch['branch']} sheet_assignment={quoted(row['sheet'])} "
                    f"equation_operand_tag={equations_tag} residuals={row['residuals']} "
                    f"definedness={list(row['definedness'])} verdict={row['verdict']}",
                    flush=True,
                )
            if str(branch["verdict"]).startswith("SPURIOUS_BRANCH"):
                summary["spurious_branches"] += 1
            if branch["verdict"] == "SHEET_INCOMPLETE":
                summary["sheet_incomplete_branches"] += 1
            print(
                "CONTAINMENT_BRANCH "
                f"record={quote(record)} tag={tag} branch_index={branch['index']} "
                f"branch_operand={branch['branch']} tested_sheets={branch['tested']} "
                f"total_sheets={branch['total']} "
                f"sheet_memberships={[row['verdict'] for row in branch['rows']]} "
                f"verdict={branch['verdict']}",
                flush=True,
            )
        completeness = result["completeness"]  # type: ignore[index]
        if completeness["verdict"] == "OMITTED_BRANCH":
            summary["omitted_records"] += 1
        print(
            "CONTAINMENT_COMPLETENESS "
            f"record={quote(record)} tag={tag} variables={result['variables']} "
            f"equation_operand_tag={equations_tag} "
            f"factor_covers={completeness['factor_covers']} "
            f"cover_truncated={completeness['cover_truncated']} "
            f"candidate_count={completeness['candidate_count']} "
            f"unresolved_count={completeness['unresolved_count']} "
            f"missing_memberships={completeness['missing']} "
            f"verdict={completeness['verdict']}",
            flush=True,
        )

    witness_tags = sorted(tag for tag in data if tag.endswith(WITNESS_SUFFIX))
    witness_index = 0
    for tag in witness_tags:
        stem = tag[: -len(WITNESS_SUFFIX)]
        status_tag = stem + STATUS_SUFFIX
        equations_tag = stem + EQUATIONS_SUFFIX
        operand_tag = status_tag + "_OPERANDS"
        status = data.get(status_tag)
        if status is None or status.payload != proved:
            continue
        witness_index += 1
        summary["witnesses"] += 1
        witness = data[tag]
        equation = data.get(equations_tag)
        operand = data.get(operand_tag)
        if equation is None or operand is None:
            summary["witness_parse_failures"] += 1
            missing = [name for name, value in ((equations_tag, equation), (operand_tag, operand)) if value is None]
            print(
                "CONTAINMENT_WITNESS "
                f"record={quote(record)} tag={tag} witness_operand={render_text(witness.payload)} "
                f"missing_operands={missing} residuals=NOT_COMPUTED "
                "verdict=IN_POPULATION_PARSE_FAILURE",
                flush=True,
            )
            continue
        print(
            "CONTAINMENT_WITNESS_BEGIN "
            f"record={quote(record)} index={witness_index} tag={tag} "
            f"witness_operand={render_text(witness.payload)} "
            f"equations_operand={render_text(equation.payload)} "
            f"premise_operand_record={operand_tag} state=STARTED",
            flush=True,
        )
        result, elapsed = _run_bounded(
            _witness_worker,
            (population.dialect, witness.payload, equation.payload, operand.payload),
        )
        if result["kind"] == "PARSE_FAILURE":
            summary["witness_parse_failures"] += 1
            _print_parse_failure(record, tag, result)
            continue
        if result["kind"] in {"TIME_EXPIRED", "WORKER_FAILURE"}:
            summary["witness_expiries"] += 1
            print(
                "CONTAINMENT_WITNESS "
                f"record={quote(record)} tag={tag} witness_operand={render_text(witness.payload)} "
                f"elapsed_seconds={elapsed:.6f} worker_state={result['kind']} "
                "residuals=NOT_COMPUTED verdict=WITNESS_RESOURCE_EXPIRED",
                flush=True,
            )
            continue
        if result["verdict"] == "WITNESS_FAILURE":
            summary["witness_failures"] += 1
        if result["verdict"] == "RESIDUAL_MISMATCH":
            summary["residual_mismatches"] += 1
        for row in result["rows"]:  # type: ignore[index]
            print(
                "CONTAINMENT_WITNESS_SHEET "
                f"record={quote(record)} tag={tag} witness_operand={result['witness']} "
                f"sheet_assignment={quoted(row['sheet'])} "
                f"equation_operand_tag={equations_tag} residuals={row['residuals']} "
                f"definedness={list(row['definedness'])} verdict={row['verdict']}",
                flush=True,
            )
        print(
            "CONTAINMENT_WITNESS "
            f"record={quote(record)} tag={tag} witness_operand={result['witness']} "
            f"operand_residual_match={result['residual_match']} "
            f"substituted_premises={result['substituted_premises']} "
            f"premise_truth={quoted(result['premise_truth'])} "
            f"tested_sheets={result['tested']} total_sheets={result['total']} "
            f"verdict={result['verdict']}",
            flush=True,
        )

    expiry_counts = {
        "PAIR": summary["pair_expiries"],
        "WITNESS": summary["witness_expiries"],
    }
    print(
        "CONTAINMENT_SUMMARY "
        f"record={quote(record)} "
        + " ".join(f"{key}={value}" for key, value in summary.items())
        + f" expiry_counts={expiry_counts} verdict=CENSUS_EXECUTED",
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
