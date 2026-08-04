#!/usr/bin/env python3
"""Gate emitted engine records as symbol-derived or declared premises.

The public invocation is::

    python derived_or_declared.py path/to/engine.py

For ``engine.py``, supplied-premise tags are read from ``engine.premises``.
That sidecar is plain text: one exact tag per line, with blank lines and lines
beginning with ``#`` ignored.
"""

from __future__ import annotations

import argparse
import re
import runpy
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote, unquote


MAX_PERTURBATIONS = 6
RUN_TIMEOUT_SECONDS = 590
VERDICT_TOKEN = "DERIVED_OR_DECLARED:"
TAG_PATTERN = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)[ \t]*(?::|->)[ \t]*(.*)$")
VALID_TAG_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


@dataclass(frozen=True)
class SymbolRecord:
    name: str
    assumptions: str


@dataclass
class EngineRun:
    completed: subprocess.CompletedProcess[str] | None
    tags: dict[str, str] | None
    error: str | None


def _assumption_signature(symbol: object) -> str:
    assumptions = getattr(symbol, "assumptions0", {})
    return ",".join(
        f"{key}={value}" for key, value in sorted(assumptions.items())
    )


def _write_symbol_report(path: Path, records: list[SymbolRecord]) -> None:
    lines = [
        f"{quote(record.name, safe='')}\t{quote(record.assumptions, safe='')}"
        for record in records
    ]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def _run_engine_child(
    engine: Path,
    symbol_report: Path | None,
    collapse: tuple[str, str] | None,
) -> int:
    """Patch SymPy before runpy loads the engine; this is the child process."""
    import sympy as sp

    symbol_type = sp.Symbol
    original_new = symbol_type.__new__
    records: list[SymbolRecord] = []
    recorded: set[tuple[str, str]] = set()
    shared_symbol: object | None = None

    def patched_new(cls: type, name: str, **assumptions: object) -> object:
        nonlocal shared_symbol
        symbol = original_new(cls, name, **assumptions)
        if type(symbol) is symbol_type:
            record = SymbolRecord(str(symbol.name), _assumption_signature(symbol))
            record_key = (record.name, record.assumptions)
            if record_key not in recorded:
                recorded.add(record_key)
                records.append(record)
            if collapse is not None and record.name in collapse:
                if shared_symbol is None:
                    shared_symbol = symbol
                return shared_symbol
        return symbol

    symbol_type.__new__ = staticmethod(patched_new)
    old_argv = sys.argv
    old_path = list(sys.path)
    sys.argv = [str(engine)]
    sys.path.insert(0, str(engine.parent))
    try:
        runpy.run_path(str(engine), run_name="__main__")
    finally:
        symbol_type.__new__ = original_new
        sys.argv = old_argv
        sys.path[:] = old_path
        if symbol_report is not None:
            _write_symbol_report(symbol_report, records)
    return 0


def _parse_symbol_report(path: Path) -> list[SymbolRecord]:
    records: list[SymbolRecord] = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        fields = raw_line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"malformed symbol report line {line_number}")
        records.append(SymbolRecord(unquote(fields[0]), unquote(fields[1])))
    return records


def _parse_emissions(stdout: str) -> dict[str, str]:
    """Parse tagged records, retaining continuation lines in SymPy text."""
    tags: dict[str, str] = {}
    current_tag: str | None = None
    current_lines: list[str] = []

    def finish_record() -> None:
        nonlocal current_tag, current_lines
        if current_tag is None:
            return
        if current_tag in tags:
            raise ValueError(f"duplicate emitted tag {current_tag}")
        tags[current_tag] = "\n".join(current_lines)
        current_tag = None
        current_lines = []

    for line in stdout.splitlines():
        match = TAG_PATTERN.match(line)
        if match:
            finish_record()
            current_tag = match.group(1)
            current_lines = [match.group(2)]
        elif current_tag is not None:
            current_lines.append(line)
        elif line.strip():
            raise ValueError(f"untagged output before first record: {line!r}")
    finish_record()
    if not tags:
        raise ValueError("engine emitted no tagged records")
    return tags


def _compact_reason(output: str | bytes | None) -> str:
    if output is None:
        return "no diagnostic output"
    if isinstance(output, bytes):
        output = output.decode(errors="replace")
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    if not lines:
        return "no diagnostic output"
    last_line = lines[-1]
    if len(last_line) > 280:
        last_line = last_line[:277] + "..."
    return repr(last_line)


def _execute_engine(
    engine: Path,
    *,
    symbol_report: Path | None = None,
    collapse: tuple[str, str] | None = None,
) -> EngineRun:
    command = [sys.executable, str(Path(__file__).resolve()), str(engine), "--_child"]
    if symbol_report is not None:
        command.extend(("--_symbol-report", str(symbol_report)))
    if collapse is not None:
        command.extend(("--_collapse", collapse[0], collapse[1]))
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=RUN_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as exc:
        diagnostic = exc.stderr or exc.stdout
        return EngineRun(
            completed=None,
            tags=None,
            error=(
                f"exceeded {RUN_TIMEOUT_SECONDS}s; "
                f"last_output={_compact_reason(diagnostic)}"
            ),
        )
    if completed.returncode != 0:
        diagnostic = completed.stderr if completed.stderr.strip() else completed.stdout
        return EngineRun(
            completed=completed,
            tags=None,
            error=(
                f"exit={completed.returncode}; "
                f"last_output={_compact_reason(diagnostic)}"
            ),
        )
    try:
        tags = _parse_emissions(completed.stdout)
    except ValueError as exc:
        return EngineRun(completed=completed, tags=None, error=str(exc))
    return EngineRun(completed=completed, tags=tags, error=None)


def _select_perturbations(records: list[SymbolRecord]) -> list[tuple[str, str]]:
    """Pair the first distinct symbols having identical SymPy assumptions."""
    signatures_by_name: dict[str, set[str]] = {}
    for record in records:
        signatures_by_name.setdefault(record.name, set()).add(record.assumptions)

    unique_records: list[SymbolRecord] = []
    seen_names: set[str] = set()
    for record in records:
        if record.name in seen_names:
            continue
        seen_names.add(record.name)
        if len(signatures_by_name[record.name]) == 1:
            unique_records.append(record)

    waiting: dict[str, str] = {}
    pairs: list[tuple[str, str]] = []
    for record in unique_records:
        first_name = waiting.pop(record.assumptions, None)
        if first_name is None:
            waiting[record.assumptions] = record.name
            continue
        pairs.append((first_name, record.name))
        if len(pairs) == MAX_PERTURBATIONS:
            break
    return pairs


def _premise_path(engine: Path) -> Path:
    return engine.with_suffix(".premises")


def _read_premises(path: Path) -> tuple[set[str], str, str | None]:
    if not path.exists():
        return set(), "MISSING (treated as an empty declaration)", None
    premises: set[str] = set()
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if not VALID_TAG_PATTERN.fullmatch(line):
            return (
                set(),
                "INVALID",
                f"line {line_number} is not an exact emitted-tag name: {line!r}",
            )
        if line in premises:
            return set(), "INVALID", f"duplicate tag on line {line_number}: {line}"
        premises.add(line)
    return premises, "LOADED", None


def _print_constant_records(constants: list[tuple[str, str]]) -> None:
    if not constants:
        print("CONSTANTS_FOR_ADJUDICATION: (none)")
        return
    records = " || ".join(
        f"CONSTANT {tag}: {text!r}" for tag, text in constants
    )
    print(f"CONSTANTS_FOR_ADJUDICATION: {records}")


def run_gate(engine: Path) -> int:
    engine = engine.resolve()
    if not engine.is_file():
        print(f"BASELINE: ERROR engine is not a file: {engine}")
        print(f"{VERDICT_TOKEN} ERROR")
        return 1
    if engine.suffix != ".py":
        print(f"BASELINE: ERROR only .py engine scripts are supported: {engine}")
        print(f"{VERDICT_TOKEN} ERROR")
        return 1

    with tempfile.TemporaryDirectory(prefix="derived_or_declared_") as temp_dir:
        symbol_report = Path(temp_dir) / "symbols.txt"
        baseline = _execute_engine(engine, symbol_report=symbol_report)
        if baseline.error is not None or baseline.tags is None:
            print(f"BASELINE: ERROR {baseline.error}")
            print(f"{VERDICT_TOKEN} ERROR")
            return 1
        try:
            symbol_records = _parse_symbol_report(symbol_report)
        except (OSError, ValueError) as exc:
            print(f"BASELINE: ERROR cannot read runtime symbol inventory: {exc}")
            print(f"{VERDICT_TOKEN} ERROR")
            return 1

    print(
        f"BASELINE: RAN tags={len(baseline.tags)} "
        f"symbols={len(symbol_records)}"
    )
    perturbations = _select_perturbations(symbol_records)
    completed_perturbations: list[dict[str, str]] = []
    for first_name, second_name in perturbations:
        run = _execute_engine(engine, collapse=(first_name, second_name))
        label = f"collapse({first_name}={second_name})"
        if run.error is not None or run.tags is None:
            print(f"PERTURBATION {label}: SKIPPED reason={run.error}")
            continue
        completed_perturbations.append(run.tags)
        comparable = len(set(baseline.tags).intersection(run.tags))
        print(
            f"PERTURBATION {label}: RAN tags={len(run.tags)} "
            f"comparable={comparable}"
        )

    ran_count = len(completed_perturbations)
    print(f"PERTURBATIONS_RAN: {ran_count}/{len(perturbations)}")
    constants: list[tuple[str, str]] = []
    derived: list[str] = []
    for tag, baseline_text in baseline.tags.items():
        changed = any(
            tag in perturbed_tags and perturbed_tags[tag] != baseline_text
            for perturbed_tags in completed_perturbations
        )
        if changed:
            derived.append(tag)
        else:
            constants.append((tag, baseline_text))

    _print_constant_records(constants)
    print(f"CONSTANT_COUNT: {len(constants)}")
    print(
        "CONSTANT_ADJUDICATION: this list is for adjudication; a tag may "
        "legitimately not depend on any collapsed pair, and CONSTANT is not "
        "proof of literal source text"
    )
    print(f"DERIVED_TAGS: count={len(derived)} tags={','.join(derived) or '(none)'}")

    premise_path = _premise_path(engine)
    try:
        premises, premise_status, premise_error = _read_premises(premise_path)
    except OSError as exc:
        premises, premise_status, premise_error = set(), "ERROR", str(exc)
    print(
        f"SUPPLIED_PREMISES: path={premise_path} status={premise_status} "
        f"declared={len(premises)}"
    )
    undeclared = [tag for tag, _ in constants if tag not in premises]
    print(
        f"UNDECLARED_CONSTANTS: count={len(undeclared)} "
        f"tags={','.join(undeclared) or '(none)'}"
    )

    if premise_error is not None:
        print(f"PREMISE_FILE_ERROR: {premise_error}")
        verdict = "ERROR"
    elif ran_count == 0:
        verdict = "ERROR"
    elif undeclared:
        verdict = "FAIL"
    else:
        verdict = "PASS"
    print(f"{VERDICT_TOKEN} {verdict}")
    return 0 if verdict == "PASS" else 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Classify emitted engine tags by external symbol-collapse runs.",
        epilog=(
            "For ENGINE.py, declare supplied-premise tags one per line in "
            "ENGINE.premises."
        ),
    )
    parser.add_argument("engine", type=Path, help="Python engine script to audit")
    parser.add_argument("--_child", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--_symbol-report", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--_collapse", nargs=2, help=argparse.SUPPRESS)
    arguments = parser.parse_args()
    if arguments._child:
        collapse = tuple(arguments._collapse) if arguments._collapse else None
        return _run_engine_child(
            arguments.engine.resolve(), arguments._symbol_report, collapse
        )
    if arguments._symbol_report is not None or arguments._collapse is not None:
        parser.error("internal child options require --_child")
    return run_gate(arguments.engine)


if __name__ == "__main__":
    raise SystemExit(main())
