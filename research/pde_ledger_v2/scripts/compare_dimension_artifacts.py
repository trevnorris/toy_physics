#!/usr/bin/env python3
"""Compare a fresh Python dimension sidecar with committed Mathematica DIM records."""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Mapping

import check_ledger_dimensions_pin as module_pin


SCRIPT_DIR = Path(__file__).resolve().parent
LEDGER_DIR = SCRIPT_DIR.parent
LEDGER_DIMENSIONS_MODULE = SCRIPT_DIR / "ledger_dimensions.py"
STAGE_RE = re.compile(r"(?:stage)?(?P<number>\d{3})\Z", re.IGNORECASE)
DIM_LINE_RE = re.compile(
    r"^DIM\|axes=(?P<axes>[^|]+)"
    r"\|name=(?P<name>[^|]+)"
    r"\|exponents=\{(?P<exponents>.*)\}$"
)
HEADER_LINE_RE = re.compile(
    r"^DIMENSIONS\|(?:stage=(?P<stage>[^|]+)\|)?axes=(?P<axes>[^|]+)"
    r"(?:\|source_sha256=(?P<source_sha256>[0-9a-f]{64}))?"
    r"(?:\|ledger_dimensions_sha256="
    r"(?P<ledger_dimensions_sha256>[0-9a-f]{64}))?$"
)

# Per-stage, per-engine artifact-name waivers belong here.  Every waived name
# must be explicit.
ARTIFACT_NAME_WAIVERS: Mapping[str, Mapping[str, frozenset[str]]] = {}


@dataclass(frozen=True)
class LabelledDimension:
    axes: tuple[str, ...]
    exponents: Mapping[str, Fraction]


def one_matching_file(directory: Path, pattern: str, description: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"expected one {description} matching {pattern!r}, found "
            f"{len(matches)}: {[path.name for path in matches]}"
        )
    return matches[0]


def parse_stage(raw_stage: str) -> tuple[str, Path, Path]:
    match = STAGE_RE.fullmatch(raw_stage)
    if match is None:
        raise ValueError("stage must be a three-digit number or 'stage' plus one")
    number = match.group("number")
    sidecar = one_matching_file(
        SCRIPT_DIR,
        f"ledger_stage{number}_*_sympy_audit.dimensions.txt",
        "Python dimension sidecar",
    )
    mathematica_out = one_matching_file(
        LEDGER_DIR / "mathematica" / "out",
        f"ledger_stage{number}_*_mathematica_audit.out",
        "Mathematica output",
    )
    return f"stage{number}", sidecar, mathematica_out


def build_dimension(
    *,
    name: str,
    axes: list[str],
    exponent_values: list[str],
    source: str,
) -> LabelledDimension:
    if not axes or len(set(axes)) != len(axes):
        raise ValueError(f"{source} {name!r} has invalid axes: {axes}")
    if len(axes) != len(exponent_values):
        raise ValueError(
            f"{source} {name!r} has {len(axes)} axes but "
            f"{len(exponent_values)} exponents"
        )
    try:
        exponents = {
            axis: Fraction(value.strip())
            for axis, value in zip(axes, exponent_values, strict=True)
        }
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(
            f"{source} {name!r} has an invalid exact exponent: {exc}"
        ) from exc
    return LabelledDimension(tuple(axes), exponents)


def load_dimensions(path: Path) -> dict[str, LabelledDimension]:
    dimensions: dict[str, LabelledDimension] = {}
    declared_axes: list[str] | None = None
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if line.startswith("DIMENSIONS|"):
            if declared_axes is not None:
                raise ValueError(f"{path}:{line_number} repeats DIMENSIONS header")
            header_match = HEADER_LINE_RE.fullmatch(line)
            if header_match is None:
                raise ValueError(
                    f"{path}:{line_number} has a malformed DIMENSIONS header"
                )
            declared_axes = [
                axis.strip() for axis in header_match.group("axes").split(",")
            ]
            if not declared_axes or len(set(declared_axes)) != len(declared_axes):
                raise ValueError(
                    f"{path}:{line_number} has invalid declared axes: "
                    f"{declared_axes}"
                )
            continue
        if not line.startswith("DIM|"):
            continue
        match = DIM_LINE_RE.fullmatch(line)
        if match is None:
            raise ValueError(f"{path}:{line_number} has a malformed DIM record")
        name = match.group("name")
        if name in dimensions:
            raise ValueError(f"{path}:{line_number} repeats DIM name {name!r}")
        axes = [axis.strip() for axis in match.group("axes").split(",")]
        if declared_axes is None:
            raise ValueError(
                f"{path}:{line_number} has a DIM record before its header"
            )
        if axes != declared_axes:
            raise ValueError(
                f"{path}:{line_number} record axes {axes} do not match "
                f"declared axes {declared_axes}"
            )
        dimensions[name] = build_dimension(
            name=name,
            axes=axes,
            exponent_values=match.group("exponents").split(","),
            source=str(path),
        )
    if declared_axes is None:
        raise ValueError(f"{path} has no DIMENSIONS header")
    return dimensions


def python_source_for_sidecar(sidecar: Path) -> Path:
    suffix = ".dimensions.txt"
    if not sidecar.name.endswith(suffix):
        raise ValueError(f"not a Python dimension sidecar: {sidecar}")
    source = sidecar.with_name(sidecar.name.removesuffix(suffix) + ".py")
    if not source.is_file():
        raise ValueError(
            f"Python source for dimension sidecar does not exist: {source}"
        )
    return source


def asserted_source_sha256(sidecar: Path) -> str:
    for line_number, line in enumerate(
        sidecar.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not line.startswith("DIMENSIONS|"):
            continue
        header_match = HEADER_LINE_RE.fullmatch(line)
        if header_match is None:
            raise ValueError(
                f"{sidecar}:{line_number} has a malformed DIMENSIONS header"
            )
        digest = header_match.group("source_sha256")
        if digest is None:
            raise ValueError(
                f"{sidecar}:{line_number} has no source_sha256 assertion"
            )
        return digest
    raise ValueError(f"{sidecar} has no DIMENSIONS header")


def asserted_ledger_dimensions_sha256(sidecar: Path) -> str:
    for line_number, line in enumerate(
        sidecar.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not line.startswith("DIMENSIONS|"):
            continue
        header_match = HEADER_LINE_RE.fullmatch(line)
        if header_match is None:
            raise ValueError(
                f"{sidecar}:{line_number} has a malformed DIMENSIONS header"
            )
        digest = header_match.group("ledger_dimensions_sha256")
        if digest is None:
            raise ValueError(
                f"{sidecar}:{line_number} has no ledger_dimensions_sha256 "
                f"assertion for {LEDGER_DIMENSIONS_MODULE.name}"
            )
        return digest
    raise ValueError(f"{sidecar} has no DIMENSIONS header")


def require_fresh_stage_source(sidecar: Path) -> None:
    stage_source = python_source_for_sidecar(sidecar)
    asserted = asserted_source_sha256(sidecar)
    computed = hashlib.sha256(stage_source.read_bytes()).hexdigest()
    if asserted != computed:
        raise ValueError(
            f"stale Python dimension sidecar {sidecar.name}: "
            f"asserted source_sha256={asserted}, "
            f"computed source_sha256={computed} from {stage_source.name}"
        )


def require_fresh_ledger_dimensions_module(sidecar: Path) -> None:
    asserted = asserted_ledger_dimensions_sha256(sidecar)
    computed = hashlib.sha256(LEDGER_DIMENSIONS_MODULE.read_bytes()).hexdigest()
    if asserted != computed:
        raise ValueError(
            f"stale Python dimension sidecar {sidecar.name}: "
            f"asserted ledger_dimensions_sha256={asserted}, "
            f"computed ledger_dimensions_sha256={computed} "
            f"from {LEDGER_DIMENSIONS_MODULE.name}"
        )


def require_fresh_python_sidecar(sidecar: Path) -> None:
    require_fresh_stage_source(sidecar)
    require_fresh_ledger_dimensions_module(sidecar)


def require_accepted_ledger_dimensions(sidecar: Path) -> str:
    return module_pin.require_accepted_ledger_dimensions((sidecar,))


def format_names(names: list[str]) -> str:
    return "(none)" if not names else ", ".join(names)


def format_dimension(dimension: LabelledDimension) -> str:
    parts = (
        f"{axis}={dimension.exponents[axis]}"
        for axis in dimension.axes
    )
    return "{" + ", ".join(parts) + "}"


def compare(stage: str, sidecar: Path, mathematica_out: Path) -> int:
    try:
        require_accepted_ledger_dimensions(sidecar)
    except module_pin.ModulePinError as exc:
        print(f"CONTROL_FAILURE: {exc}")
        print("COMPARISON_SKIPPED: ledger-dimensions module pin failed")
        print(f"RESULT|stage={stage}|status=FAIL|mismatches=not_checked")
        return 1

    python_dimensions = load_dimensions(sidecar)
    mathematica_dimensions = load_dimensions(mathematica_out)
    freshness_failure: str | None = None
    try:
        require_fresh_python_sidecar(sidecar)
    except ValueError as exc:
        freshness_failure = str(exc)

    python_names = set(python_dimensions)
    mathematica_names = set(mathematica_dimensions)
    shared_names = sorted(python_names & mathematica_names)
    python_only = sorted(python_names - mathematica_names)
    mathematica_only = sorted(mathematica_names - python_names)
    stage_waivers = ARTIFACT_NAME_WAIVERS.get(stage, {})
    waived_python_only = stage_waivers.get("py_only", frozenset())
    waived_mathematica_only = stage_waivers.get("wl_only", frozenset())
    unwaived_python_only = sorted(set(python_only) - waived_python_only)
    unwaived_mathematica_only = sorted(
        set(mathematica_only) - waived_mathematica_only
    )

    print(
        f"ARTIFACT_NAME_SET|stage={stage}|py={len(python_names)}"
        f"|wl={len(mathematica_names)}|shared={len(shared_names)}"
        f"|py_only={len(python_only)}|wl_only={len(mathematica_only)}"
        "|source_coverage=not_checked"
    )
    print(f"PY_ONLY: {format_names(python_only)}")
    print(f"WL_ONLY: {format_names(mathematica_only)}")
    print(
        f"ARTIFACT_NAME_WAIVERS|stage={stage}"
        f"|py_only={format_names(sorted(waived_python_only))}"
        f"|wl_only={format_names(sorted(waived_mathematica_only))}"
    )

    mismatches: int | None = None
    if freshness_failure is None:
        mismatches = 0
        for name in shared_names:
            python_dimension = python_dimensions[name]
            mathematica_dimension = mathematica_dimensions[name]
            if python_dimension.exponents == mathematica_dimension.exponents:
                continue
            mismatches += 1
            print(
                f"MISMATCH {name}: "
                f"py={format_dimension(python_dimension)}; "
                f"wl={format_dimension(mathematica_dimension)}"
            )
    else:
        print("COMPARISON_SKIPPED: Python dimension sidecar freshness failed")

    failures: list[str] = []
    if freshness_failure is not None:
        failures.append(freshness_failure)
    if len(shared_names) == 0:
        failures.append("compared=0; no shared quantities were compared")
    if unwaived_python_only:
        failures.append(
            "unwaived py_only quantities: "
            + format_names(unwaived_python_only)
        )
    if unwaived_mathematica_only:
        failures.append(
            "unwaived wl_only quantities: "
            + format_names(unwaived_mathematica_only)
        )
    for failure in failures:
        print(f"FAIL: {failure}")

    result = "PASS" if mismatches == 0 and not failures else "FAIL"
    mismatch_summary = "not_checked" if mismatches is None else str(mismatches)
    print(f"RESULT|stage={stage}|status={result}|mismatches={mismatch_summary}")
    return 0 if result == "PASS" else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", help="stage number, for example 004 or stage004")
    args = parser.parse_args()
    try:
        stage, sidecar, mathematica_out = parse_stage(args.stage)
        return compare(stage, sidecar, mathematica_out)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
