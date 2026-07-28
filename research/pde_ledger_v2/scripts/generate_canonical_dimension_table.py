#!/usr/bin/env python3
# OUTPUT ONLY: the generated canonical table observes independent stage artifacts.
# No script may ever import the generated table; doing so would make agreement tautological.
# Never hand-edit the generated table; always regenerate it with this generator.
"""Generate the human-readable canonical dimension table from engine artifacts."""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import compare_dimension_artifacts as comparator


SCRIPT_DIR = Path(__file__).resolve().parent
LEDGER_DIR = SCRIPT_DIR.parent
OUTPUT_PATH = LEDGER_DIR / "CANONICAL_DIMENSIONS.md"
SIDECAR_RE = re.compile(r"^ledger_stage(?P<number>\d{3})_.*\.dimensions\.txt$")

# Source: notes/measure_rewrite_feasibility.md.  These are the 30 stages with
# non-zero Python dimension machinery in the committed 43-script survey.
DIMENSION_BEARING_STAGES = (
    "stage002",
    "stage003",
    "stage004",
    "stage005",
    "stage006",
    "stage007",
    "stage008",
    "stage009",
    "stage010",
    "stage011",
    "stage012",
    "stage013",
    "stage016",
    "stage018",
    "stage021",
    "stage023",
    "stage027",
    "stage030",
    "stage031",
    "stage032",
    "stage034",
    "stage035",
    "stage036",
    "stage037",
    "stage038",
    "stage039",
    "stage040",
    "stage041",
    "stage042",
    "stage044",
)

CANONICAL_AXIS_ORDER = ("M", "L", "T")
SUPERSCRIPT_TRANSLATION = str.maketrans("+-0123456789", "⁺⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


@dataclass(frozen=True)
class QuantityRow:
    stage: str
    scope: str
    name: str
    python: comparator.LabelledDimension | None
    wolfram: comparator.LabelledDimension | None
    engine_status: str

    @property
    def dimensions(self) -> tuple[comparator.LabelledDimension, ...]:
        values = tuple(
            dimension
            for dimension in (self.python, self.wolfram)
            if dimension is not None
        )
        if len(values) == 2 and dimension_key(values[0]) == dimension_key(values[1]):
            return (values[0],)
        return values


@dataclass(frozen=True)
class StageData:
    stage: str
    axes: tuple[str, ...]
    python_count: int
    wolfram_count: int
    rows: tuple[QuantityRow, ...]


@dataclass(frozen=True)
class CandidateGroup:
    key: str
    rows: tuple[QuantityRow, ...]
    status: str


def dimension_key(
    dimension: comparator.LabelledDimension,
) -> tuple[tuple[str, Fraction], ...]:
    """Return a position-independent key sorted by axis label."""
    return tuple(sorted(dimension.exponents.items()))


def split_scoped_name(emitted_name: str) -> tuple[str, str]:
    """Split the final dotted component from its optional scope."""
    scope, separator, name = emitted_name.rpartition(".")
    if not separator:
        return "", emitted_name
    return scope, name


def candidate_key(name: str) -> str:
    """Normalize separator boundaries in a scope-free name, preserving case."""
    # Case is semantic in this corpus: K_dim is the EOS constant K while
    # clean_walk.k_dim is the wavenumber k; documented c_E (L T^-1) and
    # C_E (M^-1 L^-4 T^2) are also distinct quantities.  Do not case-fold.
    segments = tuple(filter(None, re.split(r"[^A-Za-z0-9]+", name)))
    if not segments:
        return ""
    # Promoting the first character after a separator is structural
    # snake_case-to-CamelCase normalization.  Initial and other case is intact.
    return segments[0] + "".join(
        segment[0].upper() + segment[1:] for segment in segments[1:]
    )


def fraction_text(value: Fraction) -> str:
    return str(value)


def superscript(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator).translate(SUPERSCRIPT_TRANSLATION)
    numerator = str(value.numerator).translate(SUPERSCRIPT_TRANSLATION)
    denominator = str(value.denominator).translate(SUPERSCRIPT_TRANSLATION)
    return f"{numerator}⁄{denominator}"


def canonical_render(dimension: comparator.LabelledDimension) -> str:
    """Render by canonical axis label order, never by incoming tuple position."""
    known_axes = [
        axis for axis in CANONICAL_AXIS_ORDER if axis in dimension.exponents
    ]
    extra_axes = sorted(set(dimension.exponents) - set(CANONICAL_AXIS_ORDER))
    parts: list[str] = []
    for axis in known_axes + extra_axes:
        exponent = dimension.exponents[axis]
        if exponent == 0:
            continue
        parts.append(axis if exponent == 1 else f"{axis}{superscript(exponent)}")
    return "1" if not parts else " ".join(parts)


def axes_text(axes: Iterable[str]) -> str:
    return "(" + ",".join(axes) + ")"


def exponents_text(dimension: comparator.LabelledDimension) -> str:
    values = (fraction_text(dimension.exponents[axis]) for axis in dimension.axes)
    return "{" + ", ".join(values) + "}"


def row_axes_text(row: QuantityRow) -> str:
    assert row.python is not None or row.wolfram is not None
    if row.python is None:
        assert row.wolfram is not None
        return axes_text(row.wolfram.axes)
    if row.wolfram is None or row.python.axes == row.wolfram.axes:
        return axes_text(row.python.axes)
    return (
        f"py {axes_text(row.python.axes)}; "
        f"wl {axes_text(row.wolfram.axes)}"
    )


def row_exponents_text(row: QuantityRow) -> str:
    if row.python is None:
        assert row.wolfram is not None
        return exponents_text(row.wolfram)
    if row.wolfram is None or dimension_key(row.python) == dimension_key(row.wolfram):
        return exponents_text(row.python)
    return (
        f"py {exponents_text(row.python)}; "
        f"wl {exponents_text(row.wolfram)}"
    )


def row_render_text(row: QuantityRow) -> str:
    if row.python is None:
        assert row.wolfram is not None
        return canonical_render(row.wolfram)
    if row.wolfram is None or dimension_key(row.python) == dimension_key(row.wolfram):
        return canonical_render(row.python)
    return (
        f"py {canonical_render(row.python)}; "
        f"wl {canonical_render(row.wolfram)}"
    )


def classify_engine_status(
    *,
    stage: str,
    name: str,
    python: comparator.LabelledDimension | None,
    wolfram: comparator.LabelledDimension | None,
) -> str:
    """Use the comparator's dimensions and artifact-name waiver registry."""
    if python is not None and wolfram is not None:
        return (
            "AGREE"
            if python.exponents == wolfram.exponents
            else "DISAGREE"
        )

    stage_waivers = comparator.ARTIFACT_NAME_WAIVERS.get(stage, {})
    if python is not None:
        waived = name in stage_waivers.get("py_only", frozenset())
        return "ONE_SIDED_PY (WAIVED)" if waived else "ONE_SIDED_PY (UNWAIVED)"

    assert wolfram is not None
    waived = name in stage_waivers.get("wl_only", frozenset())
    return "ONE_SIDED_WL (WAIVED)" if waived else "ONE_SIDED_WL (UNWAIVED)"


def load_stage(stage: str) -> StageData:
    parsed_stage, sidecar, wolfram_out = comparator.parse_stage(stage)
    comparator.require_accepted_ledger_dimensions(sidecar)
    comparator.require_fresh_python_sidecar(sidecar)
    python_dimensions = comparator.load_dimensions(sidecar)
    wolfram_dimensions = comparator.load_dimensions(wolfram_out)
    emitted_names = sorted(set(python_dimensions) | set(wolfram_dimensions))

    rows: list[QuantityRow] = []
    for emitted_name in emitted_names:
        scope, name = split_scoped_name(emitted_name)
        rows.append(
            QuantityRow(
                stage=parsed_stage,
                scope=scope,
                name=name,
                python=python_dimensions.get(emitted_name),
                wolfram=wolfram_dimensions.get(emitted_name),
                engine_status=classify_engine_status(
                    stage=parsed_stage,
                    name=emitted_name,
                    python=python_dimensions.get(emitted_name),
                    wolfram=wolfram_dimensions.get(emitted_name),
                ),
            )
        )

    axis_orders = {
        dimension.axes
        for row in rows
        for dimension in (row.python, row.wolfram)
        if dimension is not None
    }
    if len(axis_orders) != 1:
        raise ValueError(
            f"{parsed_stage} has multiple declared axis orders: "
            f"{sorted(axis_orders)}"
        )

    return StageData(
        stage=parsed_stage,
        axes=next(iter(axis_orders)),
        python_count=len(python_dimensions),
        wolfram_count=len(wolfram_dimensions),
        rows=tuple(rows),
    )


def converted_stages() -> tuple[str, ...]:
    stages: list[str] = []
    for sidecar in sorted(SCRIPT_DIR.glob("*.dimensions.txt")):
        match = SIDECAR_RE.fullmatch(sidecar.name)
        if match is None:
            raise ValueError(f"unrecognized dimension sidecar name: {sidecar.name}")
        stages.append(f"stage{match.group('number')}")
    if len(stages) != len(set(stages)):
        raise ValueError(f"multiple dimension sidecars found for a stage: {stages}")
    unknown = sorted(set(stages) - set(DIMENSION_BEARING_STAGES))
    if unknown:
        raise ValueError(f"converted stages absent from the 30-stage corpus: {unknown}")
    return tuple(stages)


def build_candidate_groups(rows: tuple[QuantityRow, ...]) -> tuple[CandidateGroup, ...]:
    grouped: dict[str, list[QuantityRow]] = defaultdict(list)
    for row in rows:
        grouped[candidate_key(row.name)].append(row)

    candidates: list[CandidateGroup] = []
    for key in sorted(grouped):
        group_rows = tuple(
            sorted(grouped[key], key=lambda row: (row.stage, row.scope, row.name))
        )
        if len({row.stage for row in group_rows}) < 2:
            continue

        dimension_keys = {
            dimension_key(dimension)
            for row in group_rows
            for dimension in row.dimensions
        }
        if len(dimension_keys) > 1:
            status = "NEEDS_ADJUDICATION"
        else:
            one_sided = [
                row.engine_status
                for row in group_rows
                if row.engine_status != "AGREE"
            ]
            if not one_sided:
                status = "AGREE"
            elif any("UNWAIVED" in value for value in one_sided):
                status = "ONE_SIDED (UNWAIVED)"
            else:
                status = "ONE_SIDED (WAIVED)"
        candidates.append(CandidateGroup(key=key, rows=group_rows, status=status))
    return tuple(candidates)


def candidate_members_text(group: CandidateGroup) -> str:
    members = (
        f"{row.stage} scope `{row.scope or '(none)'}`, `{row.name}` "
        f"[{row_render_text(row)}; {row.engine_status}]"
        for row in group.rows
    )
    return "<br>".join(members)


def render_table(
    stages: tuple[StageData, ...],
    unrepresented: tuple[str, ...],
) -> str:
    rows = tuple(row for stage in stages for row in stage.rows)
    candidates = build_candidate_groups(rows)
    needs_adjudication = sum(
        group.status == "NEEDS_ADJUDICATION" for group in candidates
    )
    converted = tuple(stage.stage for stage in stages)

    lines = [
        "<!--",
        "OUTPUT ONLY: this table observes independent per-stage artifacts.",
        "No script may ever import this file; that would make cross-checks tautological.",
        "Never hand-edit this file; always regenerate it with",
        "`python3 scripts/generate_canonical_dimension_table.py`.",
        "-->",
        "",
        "# Canonical dimension table",
        "",
        "> **Coverage is incomplete:** this is the generated view of converted stages,",
        "> not a complete corpus view.",
        "",
        f"- Dimension-bearing stage corpus: **{len(DIMENSION_BEARING_STAGES)} stages**.",
        f"- Converted and represented: **{len(converted)} of "
        f"{len(DIMENSION_BEARING_STAGES)}** — {', '.join(converted)}.",
        f"- Not yet represented: {', '.join(unrepresented)}.",
        f"- Total quantity rows: **{len(rows)}** (one per exact emitted name and stage).",
        f"- Candidate-same-quantity groups: **{len(candidates)}**.",
        f"- `NEEDS_ADJUDICATION` groups: **{needs_adjudication}**.",
        "",
        "Values come only from committed `scripts/*.dimensions.txt` and",
        "`mathematica/out/*.out` artifacts. Exact emitted names are primary; candidate",
        "groups below are review flags, never automatic merges. Axis tuples are shown in",
        "each stage's declared order, while canonical renderings are always labelled and",
        "ordered `M`, `L`, `T`. Per-stage status uses the parser, labelled-axis comparison,",
        "and artifact-name waiver registry in `scripts/compare_dimension_artifacts.py`.",
        "A Python sidecar is rejected unless its source-digest assertions match",
        "independent hashes of the current stage source and `ledger_dimensions.py`;",
        "this is freshness, not source coverage.",
        "",
        "## Stage coverage",
        "",
        "| Stage | Axis order | Python quantities | Wolfram quantities | Quantity rows |",
        "|---|---:|---:|---:|---:|",
    ]
    lines.extend(
        f"| {stage.stage} | `{axes_text(stage.axes)}` | "
        f"{stage.python_count} | {stage.wolfram_count} | {len(stage.rows)} |"
        for stage in stages
    )

    lines.extend(
        [
            "",
            "## Quantities",
            "",
            "| Scope | Quantity (scope-stripped emitted name) | Candidate key | Stage | "
            "Axis order | Exponents | Canonical labelled dimension | "
            "Per-stage engine status |",
            "|---|---|---|---|---|---|---|---|",
        ]
    )
    lines.extend(
        f"| `{row.scope or '(none)'}` | `{row.name}` | "
        f"`{candidate_key(row.name)}` | {row.stage} | "
        f"`{row_axes_text(row)}` | `{row_exponents_text(row)}` | "
        f"{row_render_text(row)} | {row.engine_status} |"
        for row in rows
    )

    lines.extend(
        [
            "",
            "## Candidate-same-quantity groups",
            "",
            "A dotted prefix is parsed as scope and is not part of the quantity name.",
            "Candidate keys normalize separator boundaries in that scope-stripped name",
            "to CamelCase while preserving initial and all other case. Case is meaningful:",
            "`K_dim` is the EOS constant K while scope",
            "`clean_walk`, quantity `k_dim` is the wavenumber k; likewise, documented",
            "`c_E` (L T⁻¹) and `C_E` (M⁻¹ L⁻⁴ T²) are distinct quantities.",
            "A one-sided group is not counted as an agreement even when its visible",
            "dimensions are the same. A differing group is flagged `NEEDS_ADJUDICATION`;",
            "the generator never chooses a winner.",
            "",
            "| Case-sensitive candidate key | Members | Status |",
            "|---|---|---|",
        ]
    )
    lines.extend(
        f"| `{group.key}` | {candidate_members_text(group)} | {group.status} |"
        for group in candidates
    )
    lines.extend(
        [
            "",
            "## GROUPING LIMITATIONS",
            "",
            "The general rule actually finds `KDim` (stage011) and `K_dim`",
            "(stage012). It does not group the following known-same cross-stage",
            "quantities because their initial case differs after separator-boundary",
            "normalization:",
            "",
            "- `CsSquaredDim` (stage011) and scope `clean_walk`, quantity",
            "  `cs_squared_dim` (stage012);",
            "- `CorruptKDim` (stage011) and `corrupt_K_dim` (stage012);",
            "- `EnergyDim` (stage011) and `energy_dim` (stage012);",
            "- `FourVolumeDim` (stage011) and `four_volume_dim` (stage012);",
            "- `MassDim` (stage011) and `mass_dim` (stage012);",
            "- `OmegaDim` (stage011) and `omega_dim` (stage012);",
            "- `PressureDim` (stage011) and `pressure_dim` (stage012);",
            "- `RhoDim` (stage011) and `rho_dim` (stage012).",
            "",
            "These misses are reported only; they do not feed candidate grouping.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    try:
        comparator.module_pin.require_accepted_ledger_dimensions()
        converted = converted_stages()
        stages = tuple(load_stage(stage) for stage in converted)
        unrepresented = tuple(
            stage for stage in DIMENSION_BEARING_STAGES if stage not in converted
        )
        content = render_table(stages, unrepresented)
        OUTPUT_PATH.write_text(content, encoding="utf-8")
        rows = sum(len(stage.rows) for stage in stages)
        groups = build_candidate_groups(
            tuple(row for stage in stages for row in stage.rows)
        )
        needs = sum(group.status == "NEEDS_ADJUDICATION" for group in groups)
        print(
            f"WROTE|path={OUTPUT_PATH.relative_to(LEDGER_DIR)}"
            f"|quantities={rows}|candidate_groups={len(groups)}"
            f"|needs_adjudication={needs}"
        )
        return 0
    except (AssertionError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
