#!/usr/bin/env python3
"""Shared record enumeration for the S11 locus-census instruments.

This module deliberately contains no physics expressions.  It reads the committed
record language, identifies the census population defined by the round brief, and
reports only values computed from record lines.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterator


SOLUTION_SUFFIX = "_SOLUTION"
EQUATIONS_SUFFIX = "_EQUATIONS"
WITNESS_SUFFIX = "_REAL_WITNESS"
STATUS_SUFFIX = "_REAL_STATUS"

_TAG_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")
_STRATUM_DEFINING_RE = re.compile(r"_STRATUM[0-9]+_DEFINING_EQUATIONS$")

_WL_STATUS_UNDECIDED_RE = re.compile(
    r'"STATUS_TOKEN"\s*->\s*"UNDECIDED"'
)
_WL_SIGN_UNDECIDED_RE = re.compile(
    r'"SIGN_TOKEN"\s*->\s*"UNDECIDED"'
)
_PY_STATUS_UNDECIDED_RE = re.compile(
    r"Tuple\(Str\((['\"])STATUS_TOKEN\1\),\s*Str\((['\"])UNDECIDED\2\)\)"
)
_PY_SIGN_UNDECIDED_RE = re.compile(
    r"Tuple\(Str\((['\"])SIGN_TOKEN\1\),\s*Str\((['\"])UNDECIDED\2\)\)"
)
_WL_EXACT_UNDECIDED_RE = re.compile(r'"UNDECIDED"')
_PY_EXACT_UNDECIDED_RE = re.compile(r"Str\((['\"])UNDECIDED\1\)")


@dataclass(frozen=True)
class RecordLine:
    line_no: int
    tag: str
    payload: str


@dataclass(frozen=True)
class ContainmentPopulation:
    record: Path
    dialect: str
    raw_solution: int
    raw_equations: int
    raw_real_witness: int
    protocol_pairs: int
    semantic_pairs: int
    paired_not_applicable: int
    witness_proved_nonempty: int
    witness_not_applicable: int
    excluded_dim_pairs: int
    excluded_c1_equations: int
    excluded_stratum_defining_equations: int
    unreconciled_tags: tuple[str, ...]


@dataclass(frozen=True)
class ProbePopulation:
    record: Path
    dialect: str
    raw_undecided: int
    raw_undecided_lines: int
    status_token_occurrences: int
    status_token_lines: int
    sign_token_occurrences: int
    sign_token_lines: int
    bare_exact_occurrences: int
    embedded_text_occurrences: int
    unexplained_gap: int


def quote(value: object) -> str:
    """Stable one-line string representation used by every instrument."""
    text = str(value)
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def iter_record_lines(record: Path) -> Iterator[RecordLine]:
    with record.open("r", encoding="utf-8", errors="strict") as source:
        for line_no, raw_line in enumerate(source, 1):
            line = raw_line.rstrip("\n")
            tag, separator, payload = line.partition(": ")
            if not separator or not _TAG_RE.fullmatch(tag):
                continue
            yield RecordLine(line_no=line_no, tag=tag, payload=payload)


def read_record_map(record: Path) -> dict[str, RecordLine]:
    result: dict[str, RecordLine] = {}
    for item in iter_record_lines(record):
        # The committed transcripts repeat a few run-level progress tags.  The
        # census populations themselves are duplicate-checked by their
        # enumerators; last-value semantics here is only for auxiliary operands.
        result[item.tag] = item
    return result


def record_dialect(first_tag: str) -> str:
    if first_tag.startswith("WL_"):
        return "WL"
    if first_tag.startswith("PY_"):
        return "PY"
    return "UNKNOWN"


def sentinel_payload(dialect: str) -> str:
    if dialect == "WL":
        return '"NOT_APPLICABLE"'
    if dialect == "PY":
        return "Str('NOT_APPLICABLE')"
    return ""


def proved_nonempty_payload(dialect: str) -> str:
    if dialect == "WL":
        return '"PROVED_NONEMPTY"'
    if dialect == "PY":
        return "Str('PROVED_NONEMPTY')"
    return ""


def enumerate_containment_population(record: Path) -> ContainmentPopulation:
    selected: dict[str, RecordLine] = {}
    first_tag = ""
    duplicate_tags: list[str] = []
    for item in iter_record_lines(record):
        if not first_tag:
            first_tag = item.tag
        if item.tag.endswith(
            (SOLUTION_SUFFIX, EQUATIONS_SUFFIX, WITNESS_SUFFIX, STATUS_SUFFIX)
        ):
            if item.tag in selected:
                duplicate_tags.append(item.tag)
            selected[item.tag] = item

    dialect = record_dialect(first_tag)
    sentinel = sentinel_payload(dialect)
    proved_nonempty = proved_nonempty_payload(dialect)
    solutions = {
        tag: item for tag, item in selected.items() if tag.endswith(SOLUTION_SUFFIX)
    }
    equations = {
        tag: item for tag, item in selected.items() if tag.endswith(EQUATIONS_SUFFIX)
    }
    witnesses = {
        tag: item for tag, item in selected.items() if tag.endswith(WITNESS_SUFFIX)
    }

    protocol_pairs = 0
    semantic_pairs = 0
    paired_na = 0
    excluded_dim = 0
    excluded_c1 = 0
    excluded_stratum = 0
    witness_proved = 0
    witness_na = 0
    reconciled: set[str] = set()
    unreconciled: set[str] = set(duplicate_tags)

    for solution_tag, solution in solutions.items():
        stem = solution_tag[: -len(SOLUTION_SUFFIX)]
        equations_tag = stem + EQUATIONS_SUFFIX
        equation = equations.get(equations_tag)
        if equation is None:
            unreconciled.add(solution_tag)
            continue
        if stem.endswith("_DIM"):
            excluded_dim += 1
            reconciled.update((solution_tag, equations_tag))
            continue

        protocol_pairs += 1
        solution_is_na = solution.payload == sentinel
        equations_is_na = equation.payload == sentinel
        if solution_is_na != equations_is_na:
            unreconciled.update((solution_tag, equations_tag))
        elif solution_is_na:
            paired_na += 1
            reconciled.update((solution_tag, equations_tag))
        else:
            semantic_pairs += 1
            reconciled.update((solution_tag, equations_tag))

    for equations_tag in equations:
        if equations_tag in reconciled or equations_tag in unreconciled:
            continue
        if equations_tag.endswith("_C1_EQUATIONS"):
            excluded_c1 += 1
            reconciled.add(equations_tag)
        elif _STRATUM_DEFINING_RE.search(equations_tag):
            excluded_stratum += 1
            reconciled.add(equations_tag)
        else:
            unreconciled.add(equations_tag)

    for witness_tag, witness in witnesses.items():
        stem = witness_tag[: -len(WITNESS_SUFFIX)]
        status_tag = stem + STATUS_SUFFIX
        status = selected.get(status_tag)
        if status is None:
            unreconciled.add(witness_tag)
        elif status.payload == proved_nonempty and witness.payload != sentinel:
            witness_proved += 1
            reconciled.add(witness_tag)
        elif status.payload != proved_nonempty and witness.payload == sentinel:
            witness_na += 1
            reconciled.add(witness_tag)
        else:
            unreconciled.add(witness_tag)

    for tag in (*solutions, *equations, *witnesses):
        if tag not in reconciled and tag not in unreconciled:
            unreconciled.add(tag)

    if dialect == "UNKNOWN":
        unreconciled.add("<UNKNOWN_DIALECT>")

    return ContainmentPopulation(
        record=record,
        dialect=dialect,
        raw_solution=len(solutions),
        raw_equations=len(equations),
        raw_real_witness=len(witnesses),
        protocol_pairs=protocol_pairs,
        semantic_pairs=semantic_pairs,
        paired_not_applicable=paired_na,
        witness_proved_nonempty=witness_proved,
        witness_not_applicable=witness_na,
        excluded_dim_pairs=excluded_dim,
        excluded_c1_equations=excluded_c1,
        excluded_stratum_defining_equations=excluded_stratum,
        unreconciled_tags=tuple(sorted(unreconciled)),
    )


def enumerate_probe_population(record: Path) -> ProbePopulation:
    dialect = "UNKNOWN"
    raw_undecided = 0
    raw_undecided_lines = 0
    status_token = 0
    status_token_lines = 0
    sign_token = 0
    sign_token_lines = 0
    exact_tokens = 0
    for item in iter_record_lines(record):
        if dialect == "UNKNOWN":
            dialect = record_dialect(item.tag)
        raw_here = item.payload.count("UNDECIDED")
        raw_undecided += raw_here
        raw_undecided_lines += int(raw_here > 0)
        if dialect == "WL":
            status_here = len(_WL_STATUS_UNDECIDED_RE.findall(item.payload))
            sign_here = len(_WL_SIGN_UNDECIDED_RE.findall(item.payload))
            exact_tokens += len(_WL_EXACT_UNDECIDED_RE.findall(item.payload))
        elif dialect == "PY":
            status_here = len(_PY_STATUS_UNDECIDED_RE.findall(item.payload))
            sign_here = len(_PY_SIGN_UNDECIDED_RE.findall(item.payload))
            exact_tokens += len(_PY_EXACT_UNDECIDED_RE.findall(item.payload))
        else:
            status_here = 0
            sign_here = 0
        status_token += status_here
        sign_token += sign_here
        status_token_lines += int(status_here > 0)
        sign_token_lines += int(sign_here > 0)

    classified = status_token + sign_token
    bare_exact = max(0, exact_tokens - classified)
    embedded_text = max(0, raw_undecided - exact_tokens)
    unexplained_gap = raw_undecided - classified - bare_exact - embedded_text
    if dialect == "UNKNOWN":
        unexplained_gap = raw_undecided
        bare_exact = 0
        embedded_text = 0
    return ProbePopulation(
        record=record,
        dialect=dialect,
        raw_undecided=raw_undecided,
        raw_undecided_lines=raw_undecided_lines,
        status_token_occurrences=status_token,
        status_token_lines=status_token_lines,
        sign_token_occurrences=sign_token,
        sign_token_lines=sign_token_lines,
        bare_exact_occurrences=bare_exact,
        embedded_text_occurrences=embedded_text,
        unexplained_gap=unexplained_gap,
    )


def undecided_occurrence_classes(dialect: str, payload: str) -> tuple[str, ...]:
    if dialect == "WL":
        status_spans = [match.span() for match in _WL_STATUS_UNDECIDED_RE.finditer(payload)]
        sign_spans = [match.span() for match in _WL_SIGN_UNDECIDED_RE.finditer(payload)]
        exact_spans = [match.span() for match in _WL_EXACT_UNDECIDED_RE.finditer(payload)]
    elif dialect == "PY":
        status_spans = [match.span() for match in _PY_STATUS_UNDECIDED_RE.finditer(payload)]
        sign_spans = [match.span() for match in _PY_SIGN_UNDECIDED_RE.finditer(payload)]
        exact_spans = [match.span() for match in _PY_EXACT_UNDECIDED_RE.finditer(payload)]
    else:
        return tuple("UNEXPLAINED" for _ in re.finditer("UNDECIDED", payload))

    result: list[str] = []
    for raw in re.finditer("UNDECIDED", payload):
        position = raw.start()
        if any(start <= position < end for start, end in status_spans):
            result.append("STATUS_TOKEN")
        elif any(start <= position < end for start, end in sign_spans):
            result.append("SIGN_TOKEN")
        elif any(start <= position < end for start, end in exact_spans):
            result.append("BARE")
        else:
            result.append("EMBEDDED_TEXT")
    return tuple(result)
