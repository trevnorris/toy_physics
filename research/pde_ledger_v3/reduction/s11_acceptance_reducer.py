#!/usr/bin/env python3
"""Reduce the complete S11 calibration and census transcripts to one verdict.

The instruments deliberately exit according to execution, not acceptance.  This
separate reducer reads their complete stdout files, verifies that every expected
terminal and population record is present, and applies obligation 6 without
turning real-locus findings into acceptance bits.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


VERDICT_RE = re.compile(r"(?:^| )verdict=([A-Z0-9_]+)(?:\r)?$")
TAG_RE = re.compile(r"(?:^| )tag=([^ ]+)")

FAILURE_VERDICTS = {
    "UNRECONCILED_LINE": "UNRECONCILED_POPULATION",
    "UNRECONCILED_UNDECIDED": "UNRECONCILED_POPULATION",
    "IN_POPULATION_PARSE_FAILURE": "IN_POPULATION_PARSE_FAILURE",
    "UNPARSEABLE_OPERANDS": "IN_POPULATION_PARSE_FAILURE",
    "MISSING_OPERANDS": "IN_POPULATION_PARSE_FAILURE",
    "DECIDED_UNDECIDED_RECORD": "DECIDED_UNDECIDED_RECORD",
    "RESIDUAL_MISMATCH": "RESIDUAL_MISMATCH",
}

LIMITATION_VERDICTS = {
    "PAIR_RESOURCE_EXPIRED",
    "WITNESS_RESOURCE_EXPIRED",
    "PROBE_RESOURCE_EXPIRED",
    "SHEET_INCOMPLETE",
}


@dataclass(frozen=True)
class Transcript:
    role: str
    path: Path
    raw: bytes
    lines: tuple[str, ...]

    @classmethod
    def read(cls, role: str, path: Path) -> "Transcript":
        resolved = path.resolve(strict=True)
        raw = resolved.read_bytes()
        return cls(role, resolved, raw, tuple(raw.decode("utf-8").splitlines()))

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.raw).hexdigest()


def verdict(line: str) -> str | None:
    match = VERDICT_RE.search(line)
    return match.group(1) if match else None


def tag(line: str) -> str:
    match = TAG_RE.search(line)
    return match.group(1) if match else "NONE"


def emit_item(prefix: str, transcript: Transcript, line_number: int, line: str, kind: str) -> None:
    print(
        f"{prefix} source={transcript.path} source_line={line_number} "
        f"kind={kind} tag={tag(line)} verdict={verdict(line) or 'NONE'}",
        flush=True,
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--calibration", required=True, type=Path)
    result.add_argument("--containment-wl", required=True, type=Path)
    result.add_argument("--containment-py", required=True, type=Path)
    result.add_argument("--probe-wl", required=True, type=Path)
    result.add_argument("--probe-py", required=True, type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    calibration = Transcript.read("calibration", args.calibration)
    censuses = (
        Transcript.read("containment_wl", args.containment_wl),
        Transcript.read("containment_py", args.containment_py),
        Transcript.read("probe_wl", args.probe_wl),
        Transcript.read("probe_py", args.probe_py),
    )

    for item in (calibration, *censuses):
        print(
            f"ACCEPTANCE_INPUT role={item.role} path={item.path} bytes={len(item.raw)} "
            f"lines={len(item.lines)} sha256={item.sha256} verdict=INPUT_READ",
            flush=True,
        )

    failure_counts: Counter[str] = Counter()
    finding_counts: Counter[str] = Counter()
    limitation_counts: Counter[str] = Counter()

    calibration_misses = [
        (number, line)
        for number, line in enumerate(calibration.lines, 1)
        if verdict(line) == "CALIBRATION_MISS"
    ]
    calibration_passes = [
        line
        for line in calibration.lines
        if line.startswith("CALIBRATION_SUMMARY ") and verdict(line) == "CALIBRATION_PASS"
    ]
    if len(calibration_passes) != 1:
        failure_counts["CALIBRATION_MISS"] += 1
        print(
            "ACCEPTANCE_FAILURE source="
            f"{calibration.path} source_line=0 kind=CALIBRATION_MISS tag=NONE "
            f"verdict=CALIBRATION_SUMMARY_COUNT_{len(calibration_passes)}",
            flush=True,
        )
    for number, line in calibration_misses:
        failure_counts["CALIBRATION_MISS"] += 1
        emit_item("ACCEPTANCE_FAILURE", calibration, number, line, "CALIBRATION_MISS")

    expected_summary = {
        "containment_wl": "CONTAINMENT_SUMMARY ",
        "containment_py": "CONTAINMENT_SUMMARY ",
        "probe_wl": "PROBE_SUMMARY ",
        "probe_py": "PROBE_SUMMARY ",
    }
    for item in censuses:
        population = [line for line in item.lines if line.endswith("verdict=POPULATION_RECONCILED")]
        summaries = [
            line
            for line in item.lines
            if line.startswith(expected_summary[item.role]) and verdict(line) == "CENSUS_EXECUTED"
        ]
        if len(population) != 1:
            failure_counts["UNRECONCILED_POPULATION"] += 1
            print(
                f"ACCEPTANCE_FAILURE source={item.path} source_line=0 "
                "kind=UNRECONCILED_POPULATION tag=NONE "
                f"verdict=POPULATION_RECONCILED_COUNT_{len(population)}",
                flush=True,
            )
        if len(summaries) != 1:
            failure_counts["INCOMPLETE_CENSUS"] += 1
            print(
                f"ACCEPTANCE_FAILURE source={item.path} source_line=0 "
                "kind=INCOMPLETE_CENSUS tag=NONE "
                f"verdict=CENSUS_EXECUTED_COUNT_{len(summaries)}",
                flush=True,
            )

        for number, line in enumerate(item.lines, 1):
            token = verdict(line)
            if token in FAILURE_VERDICTS:
                kind = FAILURE_VERDICTS[token]
                failure_counts[kind] += 1
                emit_item("ACCEPTANCE_FAILURE", item, number, line, kind)
            if token in LIMITATION_VERDICTS:
                limitation_counts[token] += 1
                emit_item("ACCEPTANCE_LIMITATION", item, number, line, token)

            finding_kind: str | None = None
            if line.startswith("CONTAINMENT_BRANCH ") and token and token.startswith("SPURIOUS_BRANCH"):
                finding_kind = "SPURIOUS_BRANCH"
            elif line.startswith("CONTAINMENT_COMPLETENESS ") and token == "OMITTED_BRANCH":
                finding_kind = "OMITTED_BRANCH"
            elif line.startswith("CONTAINMENT_WITNESS ") and token == "WITNESS_FAILURE":
                finding_kind = "WITNESS_FAILURE"
            if finding_kind is not None:
                finding_counts[finding_kind] += 1
                emit_item("ACCEPTANCE_FINDING", item, number, line, finding_kind)

    failure_total = sum(failure_counts.values())
    finding_total = sum(finding_counts.values())
    limitation_total = sum(limitation_counts.values())
    failure_fields = " ".join(
        f"{key.lower()}={failure_counts[key]}"
        for key in (
            "CALIBRATION_MISS",
            "UNRECONCILED_POPULATION",
            "IN_POPULATION_PARSE_FAILURE",
            "DECIDED_UNDECIDED_RECORD",
            "RESIDUAL_MISMATCH",
            "INCOMPLETE_CENSUS",
        )
    )
    finding_fields = " ".join(
        f"{key.lower()}={finding_counts[key]}"
        for key in ("SPURIOUS_BRANCH", "OMITTED_BRANCH", "WITNESS_FAILURE")
    )
    limitation_fields = " ".join(
        f"{key.lower()}={limitation_counts[key]}"
        for key in sorted(LIMITATION_VERDICTS)
    )
    print(
        f"ACCEPTANCE_COUNTS failures={failure_total} {failure_fields} "
        f"findings={finding_total} {finding_fields} "
        f"limitations={limitation_total} {limitation_fields} verdict=COUNTS_COMPUTED",
        flush=True,
    )
    round_verdict = "ROUND_PASS" if failure_total == 0 else "ROUND_FAIL"
    print(
        f"ACCEPTANCE_SUMMARY failures={failure_total} findings={finding_total} "
        f"limitations={limitation_total} verdict={round_verdict}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
