#!/usr/bin/env python3
"""Classify every residual relation by dimensional homogeneity.

The gate intentionally reads the transport documents without constructing a
Registry.  A damaged or incomplete dimension declaration must therefore land
in UNDETERMINED instead of being hidden behind an earlier schema exception.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml


HERE = Path(__file__).resolve().parent
QUANTITIES_PATH = HERE / "quantities.yaml"
RELATIONS_PATH = HERE / "relations.yaml"

HOMOGENEOUS = "HOMOGENEOUS"
INHOMOGENEOUS = "INHOMOGENEOUS"
UNDETERMINED = "UNDETERMINED"

Dimension = tuple[Fraction, ...]


class DimensionFinding(ValueError):
    """An expected, reportable dimensional finding."""


class InhomogeneousFinding(DimensionFinding):
    """Additive terms have unequal dimensions."""


class UndeterminedFinding(DimensionFinding):
    """The expression's dimension cannot be established."""


@dataclass(frozen=True)
class RelationResult:
    relation_id: str
    status: str
    dimension: Dimension | None
    detail: str
    row: Mapping[str, Any]


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        document = yaml.safe_load(handle)
    if not isinstance(document, dict):
        raise ValueError(f"{path}: top level must be a mapping")
    return document


def _format_locus(locus: Any) -> str:
    if not isinstance(locus, Mapping):
        return "<missing>"
    path = locus.get("path", "<missing-path>")
    start = locus.get("line_start", "?")
    end = locus.get("line_end", "?")
    return f"{path}:{start}-{end}"


def _format_dimension(dimension: Dimension | None) -> str:
    if dimension is None:
        return "<unknown>"

    def render(value: Fraction) -> str:
        return str(value.numerator) if value.denominator == 1 else str(value)

    return "[" + ",".join(render(value) for value in dimension) + "]"


class DimensionAudit:
    def __init__(
        self,
        quantities_document: Mapping[str, Any],
        relations_document: Mapping[str, Any],
    ) -> None:
        convention = quantities_document.get("dimension_convention")
        if not isinstance(convention, Mapping):
            raise ValueError("quantities.yaml: dimension_convention must be a mapping")
        self.convention_name = convention.get("name")
        bases = convention.get("ordered_bases")
        if not isinstance(self.convention_name, str) or not isinstance(bases, list) or not bases:
            raise ValueError("quantities.yaml: malformed declared dimension convention")
        self.bases = tuple(str(base) for base in bases)

        rows = quantities_document.get("quantities")
        if not isinstance(rows, list):
            raise ValueError("quantities.yaml: quantities must be a list")
        self.quantity_rows: dict[str, Mapping[str, Any]] = {}
        for index, row in enumerate(rows):
            if not isinstance(row, Mapping) or not isinstance(row.get("qid"), str):
                raise ValueError(f"quantities.yaml: malformed quantity row {index}")
            qid = str(row["qid"])
            if qid in self.quantity_rows:
                raise ValueError(f"quantities.yaml: duplicate QID {qid}")
            self.quantity_rows[qid] = row

        relation_rows = relations_document.get("relations")
        if not isinstance(relation_rows, list):
            raise ValueError("relations.yaml: relations must be a list")
        self.relation_rows = tuple(relation_rows)

    def _dimension_metadata(
        self, qid: str
    ) -> tuple[Dimension, str, bool, Mapping[str, Any]]:
        row = self.quantity_rows.get(qid)
        if row is None:
            raise UndeterminedFinding(f"{qid}: QID has no quantity record")
        declaration = row.get("dimension")
        if not isinstance(declaration, Mapping):
            raise UndeterminedFinding(f"{qid}: missing dimension declaration")
        if declaration.get("convention") != self.convention_name:
            raise UndeterminedFinding(
                f"{qid}: dimension convention {declaration.get('convention')!r} "
                f"does not match declared {self.convention_name!r}"
            )
        exponents = declaration.get("exponents")
        if (
            not isinstance(exponents, list)
            or len(exponents) != len(self.bases)
            or any(isinstance(value, bool) or not isinstance(value, int) for value in exponents)
        ):
            raise UndeterminedFinding(
                f"{qid}: dimension exponents must be {len(self.bases)} integers"
            )
        provenance = declaration.get("provenance")
        if not isinstance(provenance, Mapping):
            raise UndeterminedFinding(f"{qid}: missing dimension provenance")
        stage_id = provenance.get("stage_id")
        on_shared = provenance.get("stage_uses_shared_dimensions_module")
        locus = provenance.get("source_locus")
        if (
            not isinstance(stage_id, str)
            or not isinstance(on_shared, bool)
            or not isinstance(locus, Mapping)
        ):
            raise UndeterminedFinding(f"{qid}: incomplete dimension provenance")
        return (
            tuple(Fraction(value) for value in exponents),
            stage_id,
            on_shared,
            locus,
        )

    @staticmethod
    def _require_arity(operator: str, arguments: Sequence[Any], expected: int | str) -> None:
        valid = (
            len(arguments) >= int(expected[2:])
            if isinstance(expected, str)
            else len(arguments) == expected
        )
        if not valid:
            raise UndeterminedFinding(
                f"{operator}: expected arity {expected}, observed {len(arguments)}"
            )

    def dimension_of(self, node: Any) -> Dimension:
        zero = tuple(Fraction(0) for _ in self.bases)
        if isinstance(node, bool):
            raise UndeterminedFinding("boolean is not a numeric literal")
        if isinstance(node, int):
            return zero
        if not isinstance(node, list) or not node or not isinstance(node[0], str):
            raise UndeterminedFinding(f"malformed expression node: {node!r}")

        operator = node[0]
        arguments = node[1:]
        if operator == "Q":
            self._require_arity(operator, arguments, 1)
            if not isinstance(arguments[0], str):
                raise UndeterminedFinding(f"Q leaf is not a QID string: {node!r}")
            dimension, _, _, _ = self._dimension_metadata(arguments[0])
            return dimension
        if operator == "Rat":
            self._require_arity(operator, arguments, 2)
            if (
                any(isinstance(value, bool) or not isinstance(value, int) for value in arguments)
                or arguments[1] == 0
            ):
                raise UndeterminedFinding(f"Rat is not an exact numeric literal: {node!r}")
            return zero
        if operator in {"Add", "Mul"}:
            self._require_arity(operator, arguments, ">=2")
        elif operator in {"Sub", "Div", "Pow"}:
            self._require_arity(operator, arguments, 2)
        elif operator in {"Neg", "Sqrt"}:
            self._require_arity(operator, arguments, 1)
        else:
            raise UndeterminedFinding(f"unsupported operator {operator!r}")

        if operator == "Pow":
            exponent = arguments[1]
            if isinstance(exponent, bool) or not isinstance(exponent, int):
                raise UndeterminedFinding(
                    f"Pow exponent must be a bare integer, observed {exponent!r}"
                )
            base_dimension = self.dimension_of(arguments[0])
            return tuple(value * exponent for value in base_dimension)

        dimensions = tuple(self.dimension_of(argument) for argument in arguments)
        if operator in {"Add", "Sub"}:
            if any(value != dimensions[0] for value in dimensions[1:]):
                rendered = ", ".join(_format_dimension(value) for value in dimensions)
                raise InhomogeneousFinding(
                    f"{operator} additive term dimensions differ: {rendered}"
                )
            return dimensions[0]
        if operator == "Mul":
            return tuple(sum(values, Fraction(0)) for values in zip(*dimensions))
        if operator == "Div":
            return tuple(left - right for left, right in zip(*dimensions))
        if operator == "Neg":
            return dimensions[0]
        if operator == "Sqrt":
            return tuple(value / 2 for value in dimensions[0])
        raise AssertionError(f"unhandled supported operator {operator}")

    @staticmethod
    def _collect_qids(node: Any) -> tuple[str, ...]:
        if not isinstance(node, list) or not node:
            return ()
        if node[0] == "Q" and len(node) == 2 and isinstance(node[1], str):
            return (node[1],)
        result: list[str] = []
        for argument in node[1:]:
            result.extend(DimensionAudit._collect_qids(argument))
        return tuple(dict.fromkeys(result))

    def classify_relations(self) -> tuple[RelationResult, ...]:
        results: list[RelationResult] = []
        for index, row in enumerate(self.relation_rows):
            if not isinstance(row, Mapping):
                raise ValueError(f"relations.yaml: malformed relation row {index}")
            residual = row.get("residual")
            if residual is None:
                continue
            relation_id = str(row.get("relation_id", f"<row-{index}>"))
            try:
                dimension = self.dimension_of(residual)
            except InhomogeneousFinding as finding:
                results.append(
                    RelationResult(
                        relation_id, INHOMOGENEOUS, None, str(finding), row
                    )
                )
            except UndeterminedFinding as finding:
                results.append(
                    RelationResult(
                        relation_id, UNDETERMINED, None, str(finding), row
                    )
                )
            else:
                results.append(
                    RelationResult(
                        relation_id,
                        HOMOGENEOUS,
                        dimension,
                        "all additive terms agree",
                        row,
                    )
                )
        return tuple(results)

    def print_provenance(self) -> int:
        failures = 0
        for qid in self.quantity_rows:
            try:
                dimension, stage_id, on_shared, locus = self._dimension_metadata(qid)
            except UndeterminedFinding as finding:
                failures += 1
                print(f"DIMENSION_PROVENANCE {qid}: UNDETERMINED {finding}")
                continue
            print(
                f"DIMENSION_PROVENANCE {qid}: dimension={_format_dimension(dimension)} "
                f"stage={stage_id} shared_module={'yes' if on_shared else 'no'} "
                f"locus={_format_locus(locus)}"
            )
        return failures

    def print_candidate_culprits(self, result: RelationResult) -> None:
        row = result.row
        print(
            f"  CANDIDATE_CULPRIT relation={result.relation_id} "
            f"source={_format_locus(row.get('source_locus'))} "
            f"execution={_format_locus(row.get('execution_locus'))}"
        )
        qids = self._collect_qids(row.get("residual"))
        declarations: list[str] = []
        for qid in qids:
            quantity = self.quantity_rows.get(qid, {})
            declaration = quantity.get("dimension") if isinstance(quantity, Mapping) else None
            if not isinstance(declaration, Mapping):
                declarations.append(f"{qid}=<missing> locus=<missing>")
                continue
            raw_exponents = declaration.get("exponents", "<missing>")
            provenance = declaration.get("provenance")
            locus = provenance.get("source_locus") if isinstance(provenance, Mapping) else None
            declarations.append(
                f"{qid}={raw_exponents!r} locus={_format_locus(locus)}"
            )
        print("  CANDIDATE_CULPRIT dimensions: " + "; ".join(declarations))


def main() -> int:
    quantities = _read_yaml(QUANTITIES_PATH)
    relations = _read_yaml(RELATIONS_PATH)
    audit = DimensionAudit(quantities, relations)
    provenance_failures = audit.print_provenance()
    results = audit.classify_relations()

    populations = {
        status: tuple(result for result in results if result.status == status)
        for status in (HOMOGENEOUS, INHOMOGENEOUS, UNDETERMINED)
    }
    for result in results:
        if result.status == HOMOGENEOUS:
            print(
                f"RELATION {result.relation_id}: {result.status} "
                f"dimension={_format_dimension(result.dimension)}"
            )
        else:
            print(f"RELATION {result.relation_id}: {result.status} detail={result.detail}")
            audit.print_candidate_culprits(result)

    for status in (HOMOGENEOUS, INHOMOGENEOUS, UNDETERMINED):
        relation_ids = ",".join(result.relation_id for result in populations[status]) or "<none>"
        print(f"{status} ({len(populations[status])}): {relation_ids}")
    print(
        "POPULATION_COUNTS: "
        + " ".join(
            f"{status}={len(populations[status])}"
            for status in (HOMOGENEOUS, INHOMOGENEOUS, UNDETERMINED)
        )
    )

    failed = (
        bool(populations[INHOMOGENEOUS])
        or bool(populations[UNDETERMINED])
        or provenance_failures > 0
    )
    print(f"DIMENSIONAL_HOMOGENEITY_GATE: {'FAIL' if failed else 'PASS'}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
