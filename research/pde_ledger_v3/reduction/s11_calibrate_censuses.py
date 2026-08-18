#!/usr/bin/env python3
"""Build byte-shaped planted records and run every production entrypoint first."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import shlex
import subprocess
import sys

from s11_census_common import read_record_map


REPO = Path("/var/projects/toy_physics")
REDUCTION = REPO / "research/pde_ledger_v3/reduction"
WL_RECORD = REPO / "research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out"
PY_RECORD = REPO / "research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out"
DEFAULT_SCRATCH = Path("/home/trevnorris/.s11_build/census_build2")
EMPTY_WITNESS_STEM = "WL_S11_MAIN_D2_STRATUM1_POINT_EVIDENCE_ROOT_COINCIDENCE_COEFF"


@dataclass(frozen=True)
class Assertion:
    name: str
    fragments: tuple[str, ...]
    present: bool = True


@dataclass(frozen=True)
class Case:
    name: str
    command: tuple[str, ...]
    assertions: tuple[Assertion, ...]


def _line(data: dict[str, object], tag: str, payload: str | None = None) -> str:
    item = data[tag]
    value = item.payload if payload is None else payload  # type: ignore[attr-defined]
    return f"{tag}: {value}\n"


def _literal(tag: str, payload: str) -> str:
    return f"{tag}: {payload}\n"


def _write_record(path: Path, lines: list[str]) -> None:
    path.write_text("".join(lines), encoding="utf-8")


def build_plants(scratch: Path) -> dict[str, Path]:
    wl = read_record_map(WL_RECORD)
    py = read_record_map(PY_RECORD)

    # Preserve the round-1 byte copies that exercise attempt provenance, a
    # planted omission, and the exact paired sentinel.
    omitted_stem = "WL_S11_MAIN_D2_ROOT1_RANK_DROP_COEFF"
    omitted_payload = wl[omitted_stem + "_SOLUTION"].payload  # type: ignore[attr-defined]
    old_fragment = '"SOLUTION_SET" -> {{muR -> bComp}}'
    if omitted_payload.count(old_fragment) != 1:
        raise RuntimeError("omitted-branch source shape changed")
    omitted_payload = omitted_payload.replace(old_fragment, '"SOLUTION_SET" -> {}', 1)
    provenance_stem = "WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS"
    na_stem = "WL_S11_MAIN_D2_STRATUM1_ROOT_COUNT_ALL_CHANGE_LOCUS"
    empty_witness_stem = EMPTY_WITNESS_STEM
    legacy_wl = scratch / "calibration_legacy_shapes_wl.out"
    _write_record(
        legacy_wl,
        [
            _line(wl, omitted_stem + "_EQUATIONS"),
            _line(wl, omitted_stem + "_SOLUTION", omitted_payload),
            _line(wl, provenance_stem + "_EQUATIONS"),
            _line(wl, provenance_stem + "_SOLUTION"),
            _line(wl, na_stem + "_EQUATIONS"),
            _line(wl, na_stem + "_SOLUTION"),
            _line(wl, empty_witness_stem + "_EQUATIONS"),
            _line(wl, empty_witness_stem + "_SOLUTION"),
            _line(wl, empty_witness_stem + "_REAL_STATUS"),
            _line(wl, empty_witness_stem + "_REAL_WITNESS"),
            _line(wl, empty_witness_stem + "_REAL_STATUS_OPERANDS"),
        ],
    )

    # Compact WL records use the same tag/payload grammar as the committed
    # lines while isolating the repaired radical, union-containment, and empty
    # branch semantics.
    semantic_wl = scratch / "calibration_semantics_wl.out"
    pairs = (
        (
            "WL_S11_CAL_RADICAL_ZERO",
            "{x == Sqrt[y]}",
            '{"SOLVE_VARIABLES" -> {x}, "SOLUTION_SET" -> {{x -> Sqrt[y]}}}',
        ),
        (
            "WL_S11_CAL_RADICAL_FAIL",
            "{x == Sqrt[y]}",
            '{"SOLVE_VARIABLES" -> {x}, "SOLUTION_SET" -> {{x -> 1 + Sqrt[y]}}}',
        ),
        (
            "WL_S11_CAL_CHART_COVERED",
            "{k1^2 + k2^2 == 0}",
            '{"SOLVE_VARIABLES" -> {k1, k2}, "SOLUTION_SET" -> {{k2 -> -I*k1}, {k2 -> I*k1}}}',
        ),
        (
            "WL_S11_CAL_CHART_OMITTED",
            "{k1^2 + k2^2 == 0}",
            '{"SOLVE_VARIABLES" -> {k1, k2}, "SOLUTION_SET" -> {{k2 -> -I*k1}}}',
        ),
        (
            "WL_S11_CAL_EMPTY_BRANCH",
            "{True}",
            '{"SOLVE_VARIABLES" -> {x}, "SOLUTION_SET" -> {{}}}',
        ),
    )
    semantic_lines: list[str] = []
    for stem, equations, solution in pairs:
        semantic_lines.extend(
            (
                _literal(stem + "_EQUATIONS", equations),
                _literal(stem + "_SOLUTION", solution),
            )
        )
    _write_record(semantic_wl, semantic_lines)

    probe_wl = scratch / "calibration_probe_parsers_wl.out"
    zero_sample = "f[x]*(x + 2)*(x + 1)*(2*x + 1)*(2*x - 1)"
    _write_record(
        probe_wl,
        [
            _literal(
                "WL_S11_CAL_UNEQUAL_REAL_ADMISSIBLE",
                '{{"STATUS_TOKEN" -> "UNDECIDED", "OPERANDS" -> {"BRANCH" -> {k1 -> 1, s -> 1/2}, "PREMISES" -> {Unequal[k1, 0], s != 1}}}}',
            ),
            _literal("WL_S11_CAL_ABS_EQUATIONS", "{Abs[0] == 0}"),
            _literal(
                "WL_S11_CAL_ABS_IDENTICALLY_SATISFIED",
                '{"STATUS_TOKEN" -> "UNDECIDED", "OPERANDS" -> {"RESIDUALS" -> {Abs[0]}}}',
            ),
            _literal(
                "WL_S11_CAL_ZERO_SAMPLE_EQUATIONS", f"{{{zero_sample} == 0}}"
            ),
            _literal(
                "WL_S11_CAL_ZERO_SAMPLE_IDENTICALLY_SATISFIED",
                '{"STATUS_TOKEN" -> "UNDECIDED", "OPERANDS" -> {"RESIDUALS" -> {'
                + zero_sample
                + "}}}",
            ),
        ],
    )

    emptyset_tag = "PY_S11_XFORM_CURLONLY_D2_ROOT2_KW_ZERO_LOCUS_INCONSISTENT"
    emptyset_py = scratch / "calibration_emptyset_py.out"
    _write_record(emptyset_py, [_line(py, emptyset_tag)])

    single_py = scratch / "calibration_single_assignment_py.out"
    single_payload = (
        "Tuple(Tuple(Tuple(Str('STATUS_TOKEN'), Str('UNDECIDED')), "
        "Tuple(Str('OPERANDS'), Tuple(Tuple(Tuple(Symbol('k1', real=True), "
        "Integer(1))), Tuple(Equality(Symbol('k1', real=True), Integer(1)))))))"
    )
    _write_record(
        single_py,
        [_literal("PY_S11_CAL_SINGLE_REAL_ADMISSIBLE", single_payload)],
    )

    premise_py = scratch / "calibration_premise_witness_py.out"
    premise_stem = "PY_S11_CAL_PREMISE"
    b_symbol = "Symbol('B_comp', positive=True)"
    equations = f"Tuple(Add({b_symbol}, Integer(1)))"
    _write_record(
        premise_py,
        [
            _literal(premise_stem + "_EQUATIONS", equations),
            _literal(
                premise_stem + "_SOLUTION",
                f"Tuple(Tuple(Tuple({b_symbol}, Integer(-1))))",
            ),
            _literal(premise_stem + "_REAL_STATUS", "Str('PROVED_NONEMPTY')"),
            _literal(
                premise_stem + "_REAL_WITNESS",
                f"Tuple(Tuple({b_symbol}, Integer(-1)))",
            ),
            _literal(
                premise_stem + "_REAL_STATUS_OPERANDS",
                f"Tuple({equations}, Tuple({b_symbol}), "
                f"Tuple(StrictGreaterThan({b_symbol}, Integer(0))))",
            ),
        ],
    )

    # Keep a real ConditionSet-bearing line whose live operands independently
    # produce a concrete consistency witness.
    condition_tag = "PY_S11_MAIN_D2_ROOT1_RANK_DROP_K_INCONSISTENT"
    condition_py = scratch / "calibration_condition_set_py.out"
    _write_record(
        condition_py,
        [_line(py, condition_tag)],
    )
    return {
        "legacy_wl": legacy_wl,
        "semantic_wl": semantic_wl,
        "probe_wl": probe_wl,
        "emptyset_py": emptyset_py,
        "single_py": single_py,
        "premise_py": premise_py,
        "condition_py": condition_py,
    }


def _run(command: tuple[str, ...], log: Path) -> tuple[int, str]:
    completed = subprocess.run(
        command,
        cwd=str(REPO),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    log.write_text(completed.stdout, encoding="utf-8")
    return completed.returncode, completed.stdout


def _matching_lines(stdout: str, fragments: tuple[str, ...]) -> list[str]:
    return [line for line in stdout.splitlines() if all(part in line for part in fragments)]


def _cases(plants: dict[str, Path]) -> tuple[Case, ...]:
    containment = str(REDUCTION / "s11_containment_census.py")
    probe = str(REDUCTION / "s11_undecided_probe_census.py")
    reducer = str(REDUCTION / "s11_acceptance_reducer.py")

    def census(script: str, plant: str) -> tuple[str, ...]:
        return (sys.executable, script, str(plants[plant]))

    return (
        Case(
            "legacy_byte_shapes_wl",
            census(containment, "legacy_wl"),
            (
                Assertion("planted_omission", ("CONTAINMENT_COMPLETENESS", "WL_S11_MAIN_D2_ROOT1_RANK_DROP_COEFF_SOLUTION", "verdict=OMITTED_BRANCH")),
                Assertion("attempt_provenance", ("CONTAINMENT_BEGIN", "WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_SOLUTION")),
                Assertion("paired_sentinel", ("CONTAINMENT_SENTINEL", "verdict=PAIRED_NOT_APPLICABLE")),
                Assertion("nested_empty_witness_operands", ("CONTAINMENT_WITNESS ", EMPTY_WITNESS_STEM + "_REAL_WITNESS", "verdict=WITNESS_UNDECIDED")),
            ),
        ),
        Case(
            "sheet_union_and_empty_branch_wl",
            census(containment, "semantic_wl"),
            (
                Assertion("as_written_radical_zero", ("CONTAINMENT_BRANCH ", "WL_S11_CAL_RADICAL_ZERO_SOLUTION", "verdict=BRANCH_CONTAINED")),
                Assertion("as_written_not_spurious", ("CONTAINMENT_BRANCH ", "WL_S11_CAL_RADICAL_ZERO_SOLUTION", "verdict=SPURIOUS_BRANCH"), False),
                Assertion("all_coherent_sheets_fail", ("CONTAINMENT_BRANCH ", "WL_S11_CAL_RADICAL_FAIL_SOLUTION", "verdict=SPURIOUS_BRANCH")),
                Assertion("chart_union_covered", ("CONTAINMENT_COMPLETENESS", "WL_S11_CAL_CHART_COVERED_SOLUTION", "verdict=COMPLETE_FACTOR_COVER")),
                Assertion("chart_not_omitted", ("CONTAINMENT_COMPLETENESS", "WL_S11_CAL_CHART_COVERED_SOLUTION", "verdict=OMITTED_BRANCH"), False),
                Assertion("genuine_component_omitted", ("CONTAINMENT_COMPLETENESS", "WL_S11_CAL_CHART_OMITTED_SOLUTION", "verdict=OMITTED_BRANCH")),
                Assertion("nested_empty_branch", ("CONTAINMENT_BRANCH ", "WL_S11_CAL_EMPTY_BRANCH_SOLUTION", "verdict=BRANCH_CONTAINED")),
            ),
        ),
        Case(
            "wl_boolean_abs_and_zero_sample",
            census(probe, "probe_wl"),
            (
                Assertion("unequal_boolean", ("PROBE_RECORD", "WL_S11_CAL_UNEQUAL_REAL_ADMISSIBLE", "Unequality(Symbol('s'), Integer(1))", "verdict=STILL_UNDECIDED")),
                Assertion("abs_zero", ("PROBE_RECORD", "WL_S11_CAL_ABS_IDENTICALLY_SATISFIED", "verdict=DECIDED_UNDECIDED_RECORD")),
                Assertion("zero_samples_not_proof", ("PROBE_RECORD", "WL_S11_CAL_ZERO_SAMPLE_IDENTICALLY_SATISFIED", "UNDECIDED_ZERO_SAMPLES", "verdict=STILL_UNDECIDED")),
            ),
        ),
        Case(
            "empty_solver_return_py",
            census(probe, "emptyset_py"),
            (
                Assertion("emptyset_is_opaque", ("PROBE_RECORD", "UNDECIDED_CONFIRMED", "verdict=STILL_UNDECIDED")),
                Assertion("emptyset_not_proof", ("PROBE_RECORD", "verdict=DECIDED_UNDECIDED_RECORD"), False),
            ),
        ),
        Case(
            "single_assignment_py",
            census(probe, "single_py"),
            (Assertion("single_assignment_branch", ("PROBE_RECORD", "PY_S11_CAL_SINGLE_REAL_ADMISSIBLE", "verdict=DECIDED_UNDECIDED_RECORD")),),
        ),
        Case(
            "premise_violating_witness_py",
            census(containment, "premise_py"),
            (Assertion("premise_failure", ("CONTAINMENT_WITNESS ", "PY_S11_CAL_PREMISE_REAL_WITNESS", 'premise_truth="FALSE"', "verdict=WITNESS_FAILURE")),),
        ),
        Case(
            "condition_set_and_decidable_py",
            census(probe, "condition_py"),
            (
                Assertion("condition_set_parsed", ("PROBE_BEGIN", "ConditionSet")),
                Assertion("decidable_undecided", ("PROBE_RECORD", "PY_S11_MAIN_D2_ROOT1_RANK_DROP_K_INCONSISTENT", "PROVED_CONSISTENT", "verdict=DECIDED_UNDECIDED_RECORD")),
            ),
        ),
        Case(
            "closed_reducer_taxonomy",
            (sys.executable, reducer, "--calibrate-taxonomy"),
            (
                Assertion("formerly_open_tokens_classified", ("TAXONOMY_CALIBRATION", "open_token_misses=0")),
                Assertion("unknown_is_failure", ("TAXONOMY_CALIBRATION", "unknown_bucket=FAILURE")),
            ),
        ),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH)
    args = parser.parse_args()
    scratch = args.scratch_dir.resolve(strict=True)
    if not os.access(scratch, os.W_OK):
        print(
            "CALIBRATION_ENVIRONMENT "
            f"scratch={shlex.quote(str(scratch))} writable=False "
            "verdict=CALIBRATION_EXECUTION_BLOCKED",
            flush=True,
        )
        return 73
    plants = build_plants(scratch)
    cases = _cases(plants)
    misses = 0
    assertions_run = 0
    for case in cases:
        print(
            "CALIBRATION_COMMAND "
            f"case={case.name} argv={shlex.join(case.command)} state=STARTED",
            flush=True,
        )
        returncode, stdout = _run(case.command, scratch / f"{case.name}.stdout")
        print(stdout, end="", flush=True)
        case_misses = int(returncode != 0)
        for assertion in case.assertions:
            assertions_run += 1
            matching = _matching_lines(stdout, assertion.fragments)
            detected = bool(matching) == assertion.present
            case_misses += int(not detected)
            print(
                "CALIBRATION_ASSERT "
                f"case={case.name} assertion={assertion.name} "
                f"required_presence={assertion.present} matching_lines={len(matching)} "
                f"verdict={'CALIBRATION_DETECTED' if detected else 'CALIBRATION_MISS'}",
                flush=True,
            )
        generic_failures = sum(
            stdout.count(token)
            for token in (
                "verdict=IN_POPULATION_PARSE_FAILURE",
                "verdict=UNPARSEABLE_OPERANDS",
                "verdict=MISSING_OPERANDS",
                "verdict=RESIDUAL_MISMATCH",
                "verdict=PROBE_RESOURCE_EXPIRED",
                "verdict=PAIR_RESOURCE_EXPIRED",
            )
        )
        case_misses += generic_failures
        detected = case_misses == 0
        print(
            "CALIBRATION_CASE "
            f"case={case.name} command={shlex.quote(shlex.join(case.command))} "
            f"returncode={returncode} assertions={len(case.assertions)} "
            f"generic_failures={generic_failures} "
            f"verdict={'CALIBRATION_DETECTED' if detected else 'CALIBRATION_MISS'}",
            flush=True,
        )
        if not detected:
            misses += case_misses
            break
    round_verdict = "CALIBRATION_PASS" if misses == 0 else "CALIBRATION_FAIL"
    print(
        f"CALIBRATION_SUMMARY cases={len(cases)} assertions={assertions_run} "
        f"misses={misses} verdict={round_verdict}",
        flush=True,
    )
    return 0 if misses == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
