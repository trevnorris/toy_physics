#!/usr/bin/env python3
"""Validated SymPy reader for the phase-1 QID + relation registry."""

from __future__ import annotations

import warnings
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import jsonschema
import sympy as sp
import yaml
from sympy.solvers.solvers import unrad


HERE = Path(__file__).resolve().parent
LEDGER_ROOT = HERE.parent
SCHEMA_PATH = HERE / "registry_schema.yaml"
QUANTITIES_PATH = HERE / "quantities.yaml"
RELATIONS_PATH = HERE / "relations.yaml"

EARNED_STATUS = "DERIVED-EXECUTED"
TERMINAL_SCALAR_KINDS = frozenset(
    {"parameter", "boundary-datum", "control", "discrete-choice"}
)


class RegistryError(ValueError):
    """Base class for registry/schema/admission/evaluation errors."""


class RegistryValidationError(RegistryError):
    """The transported data does not satisfy the versioned registry contract."""


class AdmissionError(RegistryError):
    """A caller tried to use a relation that the admission gate refused."""


class EvaluationError(RegistryError):
    """Forward evaluation could not be completed without freezing a derived value."""


class DeclaredValueDefaultWarning(UserWarning):
    """Forward propagation consumed an explicitly enabled declared-value default."""


@dataclass(frozen=True)
class Quantity:
    qid: str
    symbol_name: str
    kind: str
    scope: tuple[str, ...]
    regime: tuple[str, ...]
    state: str
    counting_axis: str
    dimension: tuple[int, int, int]
    value: sp.Rational | None
    aliases: tuple[str, ...]
    raw: Mapping[str, Any]


@dataclass(frozen=True)
class Relation:
    relation_id: str
    residual: sp.Expr | None
    designated_output: str | None
    input_qids: tuple[str, ...]
    provenance_status: str
    regime: tuple[str, ...]
    rhs: sp.Expr | None
    raw: Mapping[str, Any]


@dataclass(frozen=True)
class AdmissionDecision:
    relation_id: str
    admitted: bool
    reason: str


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        document = yaml.safe_load(handle)
    if not isinstance(document, dict):
        raise RegistryValidationError(f"{path}: top level must be a mapping")
    return document


def load_raw_documents(
    directory: Path | str = HERE,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Load detached raw documents, useful for mutation/able-to-fail checks."""
    base = Path(directory)
    return (
        deepcopy(_read_yaml(base / "quantities.yaml")),
        deepcopy(_read_yaml(base / "relations.yaml")),
        deepcopy(_read_yaml(base / "registry_schema.yaml")),
    )


def _schema_validate(document: Mapping[str, Any], schema: Mapping[str, Any], label: str) -> None:
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(document), key=lambda error: list(error.absolute_path))
    if errors:
        first = errors[0]
        location = ".".join(str(part) for part in first.absolute_path) or "<root>"
        raise RegistryValidationError(f"{label}:{location}: {first.message}")


class Registry:
    """Canonical quantities, parsed relations, gate decisions, rank, and evaluation."""

    def __init__(
        self,
        quantities_document: Mapping[str, Any],
        relations_document: Mapping[str, Any],
        schema: Mapping[str, Any],
        *,
        ledger_root: Path = LEDGER_ROOT,
        active_regime: str | None = None,
    ) -> None:
        _schema_validate(quantities_document, schema, "quantities.yaml")
        _schema_validate(relations_document, schema, "relations.yaml")
        self.quantities_document = deepcopy(dict(quantities_document))
        self.relations_document = deepcopy(dict(relations_document))
        self.ledger_root = Path(ledger_root)
        self.active_regime = active_regime or str(quantities_document["active_regime"])
        self._validate_transport_headers()

        self.quantities = self._build_quantities(quantities_document["quantities"])
        self.alias_to_qid = self._build_alias_map()
        self.symbols = {
            qid: sp.Symbol(quantity.symbol_name, real=True)
            for qid, quantity in self.quantities.items()
        }
        self.relations = self._build_relations(relations_document["relations"])
        self._validate_loci()
        self._validate_dimensions()
        self._validate_unique_outputs()
        self.admission_decisions, admitted_ids = self._apply_admission_gate()
        self.admitted_relations = tuple(self.relations[relation_id] for relation_id in admitted_ids)
        self.admitted_constraint_set = tuple(
            relation.residual
            for relation in self.admitted_relations
            if relation.residual is not None
        )
        self._admitted_by_output = {
            relation.designated_output: relation
            for relation in self.admitted_relations
            if relation.designated_output is not None
        }
        admitted_outputs = frozenset(self._admitted_by_output)
        self.residue = tuple(
            qid
            for qid, quantity in self.quantities.items()
            if quantity.state == "live"
            and self.active_regime in quantity.regime
            and quantity.counting_axis == "continuous-model"
            and qid not in admitted_outputs
        )

    @classmethod
    def from_documents(
        cls,
        quantities_document: Mapping[str, Any],
        relations_document: Mapping[str, Any],
        schema: Mapping[str, Any],
        *,
        ledger_root: Path = LEDGER_ROOT,
        active_regime: str | None = None,
    ) -> "Registry":
        return cls(
            quantities_document,
            relations_document,
            schema,
            ledger_root=ledger_root,
            active_regime=active_regime,
        )

    def _validate_transport_headers(self) -> None:
        if self.quantities_document["schema_version"] != self.relations_document["schema_version"]:
            raise RegistryValidationError("quantity/relation schema_version values differ")
        if self.relations_document["expression_language"] != "prefix-v1":
            raise RegistryValidationError("only expression_language=prefix-v1 is implemented")

    def _build_quantities(self, rows: Sequence[Mapping[str, Any]]) -> dict[str, Quantity]:
        result: dict[str, Quantity] = {}
        for row in rows:
            qid = str(row["qid"])
            if qid in result:
                raise RegistryValidationError(f"duplicate qid: {qid}")
            dimension = tuple(int(value) for value in row["dimension"]["exponents"])
            raw_value = row.get("value")
            if raw_value is None:
                value = None
            elif isinstance(raw_value, bool):
                raise RegistryValidationError(f"{qid}: boolean is not an exact value")
            elif isinstance(raw_value, int):
                value = sp.Rational(raw_value)
            else:
                numerator, denominator = raw_value[1:]
                value = sp.Rational(numerator, denominator)
            result[qid] = Quantity(
                qid=qid,
                symbol_name=str(row["symbol"]),
                kind=str(row["kind"]),
                scope=tuple(str(value) for value in row["scope"]),
                regime=tuple(str(value) for value in row["regime"]),
                state=str(row["state"]),
                counting_axis=str(row["counting_axis"]),
                dimension=dimension,  # type: ignore[arg-type]
                value=value,
                aliases=tuple(str(value) for value in row["aliases"]),
                raw=deepcopy(dict(row)),
            )
        return result

    def _build_alias_map(self) -> dict[str, str]:
        aliases: dict[str, str] = {}
        for qid, quantity in self.quantities.items():
            for name in (qid, quantity.symbol_name, *quantity.aliases):
                previous = aliases.get(name)
                if previous is not None and previous != qid:
                    raise RegistryValidationError(
                        f"alias collision: {name!r} resolves to both {previous} and {qid}"
                    )
                aliases[name] = qid
        return aliases

    def resolve_qid(self, name: str) -> str:
        """Resolve a canonical QID, symbol, or declared alias to one canonical QID."""
        try:
            return self.alias_to_qid[name]
        except KeyError as exc:
            raise RegistryValidationError(f"unknown QID/symbol/alias: {name!r}") from exc

    def _quantity_occurrences(self, node: Any) -> tuple[str, ...]:
        if isinstance(node, bool):
            raise RegistryValidationError("booleans are not expression integers")
        if isinstance(node, int):
            return ()
        if not isinstance(node, list) or not node:
            raise RegistryValidationError(f"expression node must be a nonempty prefix list: {node!r}")
        if node[0] == "Q":
            if len(node) != 2 or not isinstance(node[1], str):
                raise RegistryValidationError(f"Q node must be [Q, name]: {node!r}")
            return (self.resolve_qid(node[1]),)
        occurrences: list[str] = []
        for argument in node[1:]:
            occurrences.extend(self._quantity_occurrences(argument))
        return tuple(occurrences)

    def parse_expression(self, node: Any) -> sp.Expr:
        """Parse one prefix-v1 node; there is no eval/sympify engine syntax path."""
        if isinstance(node, bool):
            raise RegistryValidationError("booleans are not expression integers")
        if isinstance(node, int):
            return sp.Integer(node)
        if not isinstance(node, list) or not node or not isinstance(node[0], str):
            raise RegistryValidationError(f"malformed prefix-v1 expression: {node!r}")

        operator = node[0]
        arguments = node[1:]
        if operator == "Q":
            if len(arguments) != 1 or not isinstance(arguments[0], str):
                raise RegistryValidationError(f"Q expects one string argument: {node!r}")
            return self.symbols[self.resolve_qid(arguments[0])]
        if operator == "Rat":
            if (
                len(arguments) != 2
                or any(isinstance(value, bool) or not isinstance(value, int) for value in arguments)
                or arguments[1] == 0
            ):
                raise RegistryValidationError(f"Rat expects integer numerator/nonzero denominator: {node!r}")
            return sp.Rational(arguments[0], arguments[1])

        parsed = tuple(self.parse_expression(argument) for argument in arguments)
        if operator == "Add" and len(parsed) >= 2:
            return sp.Add(*parsed)
        if operator == "Mul" and len(parsed) >= 2:
            return sp.Mul(*parsed)
        if operator == "Sub" and len(parsed) == 2:
            return parsed[0] - parsed[1]
        if operator == "Div" and len(parsed) == 2:
            return parsed[0] / parsed[1]
        if operator == "Pow" and len(parsed) == 2:
            return sp.Pow(parsed[0], parsed[1])
        if operator == "Neg" and len(parsed) == 1:
            return -parsed[0]
        if operator == "Sqrt" and len(parsed) == 1:
            return sp.sqrt(parsed[0])
        raise RegistryValidationError(
            f"unknown operator or wrong arity for {operator}: {len(arguments)}"
        )

    def _canonical_qid_list(self, values: Iterable[str]) -> tuple[str, ...]:
        resolved = tuple(self.resolve_qid(value) for value in values)
        if len(set(resolved)) != len(resolved):
            raise RegistryValidationError(f"aliases duplicate a QID in list: {tuple(values)!r}")
        return resolved

    def _build_relations(self, rows: Sequence[Mapping[str, Any]]) -> dict[str, Relation]:
        result: dict[str, Relation] = {}
        for row in rows:
            relation_id = str(row["relation_id"])
            if relation_id in result:
                raise RegistryValidationError(f"duplicate relation_id: {relation_id}")
            raw_residual = row["residual"]
            residual = None if raw_residual is None else self.parse_expression(raw_residual)
            output = (
                None
                if row["designated_output"] is None
                else self.resolve_qid(str(row["designated_output"]))
            )
            inputs = self._canonical_qid_list(str(value) for value in row["input_qids"])

            occurrences = () if raw_residual is None else self._quantity_occurrences(raw_residual)
            rhs: sp.Expr | None = None
            if output is not None:
                expected_lhs = ["Q", output]
                if (
                    not isinstance(raw_residual, list)
                    or len(raw_residual) != 3
                    or raw_residual[0] != "Sub"
                    or raw_residual[1] != expected_lhs
                ):
                    raise RegistryValidationError(
                        f"{relation_id}: explicit residual must be [Sub,[Q,designated_output],rhs]"
                    )
                rhs_occurrences = self._quantity_occurrences(raw_residual[2])
                if output in rhs_occurrences or Counter(occurrences)[output] != 1:
                    raise RegistryValidationError(
                        f"{relation_id}: designated output recurs; implicit blocks are not phase-1 syntax"
                    )
                if set(rhs_occurrences) != set(inputs):
                    raise RegistryValidationError(
                        f"{relation_id}: input_qids do not exactly equal RHS leaves; "
                        f"declared={sorted(inputs)!r}, rhs={sorted(set(rhs_occurrences))!r}"
                    )
                rhs = self.parse_expression(raw_residual[2])
            elif set(occurrences) != set(inputs):
                raise RegistryValidationError(
                    f"{relation_id}: implicit input_qids do not exactly equal residual leaves"
                )

            for assumption in row["assumptions"]:
                self._canonical_qid_list(str(value) for value in assumption["qids"])
            for guard in row["denominator_guards"]:
                guard_qids = self._quantity_occurrences(guard)
                if not set(guard_qids) <= set(inputs):
                    raise RegistryValidationError(
                        f"{relation_id}: denominator guard uses a non-input QID"
                    )
                self.parse_expression(guard)
            self._validate_literal_consistency(relation_id, row)

            result[relation_id] = Relation(
                relation_id=relation_id,
                residual=residual,
                designated_output=output,
                input_qids=inputs,
                provenance_status=str(row["provenance_status"]),
                regime=tuple(str(value) for value in row["regime"]),
                rhs=rhs,
                raw=deepcopy(dict(row)),
            )
        return result

    @staticmethod
    def _node_at_path(node: Any, path: Sequence[int], label: str) -> Any:
        current = node
        for depth, index in enumerate(path):
            if (
                isinstance(current, bool)
                or not isinstance(current, list)
                or index < 0
                or index >= len(current)
            ):
                raise RegistryValidationError(
                    f"{label}: literal_path is invalid at depth {depth}: {tuple(path)!r}"
                )
            current = current[index]
        return current

    def _validate_literal_consistency(
        self,
        relation_id: str,
        row: Mapping[str, Any],
    ) -> None:
        """Check declared transport literals against valued registry quantities."""
        seen_paths: set[tuple[int, ...]] = set()
        for assertion in row["literal_consistency"]:
            path = tuple(int(index) for index in assertion["literal_path"])
            if path in seen_paths:
                raise RegistryValidationError(
                    f"{relation_id}: duplicate literal_consistency path: {path!r}"
                )
            seen_paths.add(path)
            observed = self._node_at_path(
                row["residual"],
                path,
                f"{relation_id}.literal_consistency",
            )
            if isinstance(observed, bool) or not isinstance(observed, int):
                raise RegistryValidationError(
                    f"{relation_id}: asserted node at {path!r} is not an integer literal: "
                    f"{observed!r}"
                )
            qid = self.resolve_qid(str(assertion["quantity"]))
            quantity_value = self.quantities[qid].value
            if quantity_value is None:
                raise RegistryValidationError(
                    f"{relation_id}: literal_consistency requires a declared value for {qid}"
                )
            expected = quantity_value + int(assertion["offset"])
            if sp.Integer(observed) != expected:
                raise RegistryValidationError(
                    f"{relation_id}: literal at {path!r} is {observed}, but "
                    f"{qid}({quantity_value}) + {int(assertion['offset'])} = {expected}"
                )

    def _all_loci(self) -> Iterable[tuple[str, Mapping[str, Any]]]:
        for quantity in self.quantities.values():
            for locus in quantity.raw["source_loci"]:
                yield quantity.qid, locus
            yield (
                f"{quantity.qid}.dimension.provenance",
                quantity.raw["dimension"]["provenance"]["source_locus"],
            )
        for relation in self.relations.values():
            yield f"{relation.relation_id}.source_locus", relation.raw["source_locus"]
            yield f"{relation.relation_id}.execution_locus", relation.raw["execution_locus"]

    def _validate_loci(self) -> None:
        line_counts: dict[Path, int] = {}
        for label, locus in self._all_loci():
            path = self.ledger_root / str(locus["path"])
            if not path.is_file():
                raise RegistryValidationError(f"{label}: stale locus; file does not exist: {path}")
            start, end = int(locus["line_start"]), int(locus["line_end"])
            if end < start:
                raise RegistryValidationError(f"{label}: locus line_end precedes line_start")
            if path not in line_counts:
                with path.open("r", encoding="utf-8") as handle:
                    line_counts[path] = sum(1 for _ in handle)
            if end > line_counts[path]:
                raise RegistryValidationError(
                    f"{label}: locus ends at {end}, but {path} has {line_counts[path]} lines"
                )

    def _dimension_of(self, node: Any) -> tuple[Fraction, Fraction, Fraction]:
        zero = (Fraction(0), Fraction(0), Fraction(0))
        if isinstance(node, bool):
            raise RegistryValidationError("booleans are not dimensionless integers")
        if isinstance(node, int):
            return zero
        if not isinstance(node, list) or not node:
            raise RegistryValidationError(f"malformed expression during dimension check: {node!r}")
        operator, arguments = node[0], node[1:]
        if operator == "Q":
            qid = self.resolve_qid(str(arguments[0]))
            return tuple(Fraction(value) for value in self.quantities[qid].dimension)  # type: ignore[return-value]
        if operator == "Rat":
            return zero
        dimensions = tuple(self._dimension_of(argument) for argument in arguments)
        if operator in {"Add", "Sub"}:
            if not dimensions or any(value != dimensions[0] for value in dimensions[1:]):
                raise RegistryValidationError(f"dimension mismatch within {operator}: {node!r}")
            return dimensions[0]
        if operator == "Mul":
            return tuple(sum(values, Fraction(0)) for values in zip(*dimensions))  # type: ignore[return-value]
        if operator == "Div":
            return tuple(left - right for left, right in zip(*dimensions))  # type: ignore[return-value]
        if operator == "Neg":
            return dimensions[0]
        if operator == "Sqrt":
            return tuple(value / 2 for value in dimensions[0])  # type: ignore[return-value]
        if operator == "Pow":
            exponent = self.parse_expression(arguments[1])
            if exponent.free_symbols or exponent.is_Rational is not True:
                raise RegistryValidationError(f"Pow exponent must be an exact rational: {node!r}")
            if dimensions[1] != zero:
                raise RegistryValidationError(f"Pow exponent must be dimensionless: {node!r}")
            factor = Fraction(int(exponent.p), int(exponent.q))
            return tuple(value * factor for value in dimensions[0])  # type: ignore[return-value]
        raise RegistryValidationError(f"dimension rule missing for operator {operator!r}")

    def _validate_dimensions(self) -> None:
        for relation in self.relations.values():
            raw_residual = relation.raw["residual"]
            if raw_residual is not None:
                self._dimension_of(raw_residual)
            for guard in relation.raw["denominator_guards"]:
                self._dimension_of(guard)

    def _validate_unique_outputs(self) -> None:
        seen: dict[str, str] = {}
        for relation in self.relations.values():
            output = relation.designated_output
            if output is None or self.active_regime not in relation.regime:
                continue
            previous = seen.get(output)
            if previous is not None:
                raise RegistryValidationError(
                    f"multiple active definitions for {output}: {previous}, {relation.relation_id}"
                )
            seen[output] = relation.relation_id

    def _base_rejection(self, relation: Relation) -> str | None:
        if relation.provenance_status != EARNED_STATUS:
            return (
                f"provenance_status={relation.provenance_status}; "
                f"only {EARNED_STATUS} enters the earned set"
            )
        if relation.residual is None:
            return "no parseable residual"
        if self.active_regime not in relation.regime:
            return f"regime incompatible with {self.active_regime}"
        if relation.raw["applied_functions"]:
            return "prefix-v1 has no admitted function definition for applied_functions"
        involved = set(relation.input_qids)
        if relation.designated_output is not None:
            involved.add(relation.designated_output)
        for qid in involved:
            quantity = self.quantities[qid]
            if quantity.state != "live":
                return f"references retired QID {qid}"
            if self.active_regime not in quantity.regime:
                return f"QID {qid} is outside active regime"
        return None

    def _apply_admission_gate(
        self,
    ) -> tuple[dict[str, AdmissionDecision], tuple[str, ...]]:
        decisions: dict[str, AdmissionDecision] = {}
        pending: dict[str, Relation] = {}
        for relation_id, relation in self.relations.items():
            rejection = self._base_rejection(relation)
            if rejection is None:
                pending[relation_id] = relation
            else:
                decisions[relation_id] = AdmissionDecision(relation_id, False, rejection)

        admitted: list[str] = []
        computable = {
            qid
            for qid, quantity in self.quantities.items()
            if quantity.state == "live"
            and self.active_regime in quantity.regime
            and quantity.kind in TERMINAL_SCALAR_KINDS
        }
        while pending:
            progress = False
            for relation_id, relation in tuple(pending.items()):
                missing = set(relation.input_qids) - computable
                if missing:
                    continue
                admitted.append(relation_id)
                decisions[relation_id] = AdmissionDecision(
                    relation_id, True, "admitted: derived, executed, parseable, and closure-complete"
                )
                if relation.designated_output is not None:
                    computable.add(relation.designated_output)
                del pending[relation_id]
                progress = True
            if not progress:
                break

        for relation_id, relation in pending.items():
            missing = sorted(set(relation.input_qids) - computable)
            decisions[relation_id] = AdmissionDecision(
                relation_id,
                False,
                f"transitive closure has unresolved inputs: {missing!r}",
            )
        ordered = tuple(
            relation_id
            for relation_id in self.relations
            if decisions[relation_id].admitted
        )
        return decisions, ordered

    @property
    def active_variables(self) -> tuple[sp.Symbol, ...]:
        """Live variables on the continuous-model axis, after quotienting conventions."""
        return tuple(
            self.symbols[qid]
            for qid, quantity in self.quantities.items()
            if quantity.state == "live"
            and self.active_regime in quantity.regime
            and quantity.counting_axis == "continuous-model"
        )

    @property
    def convention_variables(self) -> tuple[sp.Symbol, ...]:
        """Live coordinates reported on the separate convention-orbit axis."""
        return tuple(
            self.symbols[qid]
            for qid, quantity in self.quantities.items()
            if quantity.state == "live"
            and self.active_regime in quantity.regime
            and quantity.counting_axis == "convention-orbit"
        )

    @property
    def discrete_structural_variables(self) -> tuple[sp.Symbol, ...]:
        """Live choices reported on the separate discrete-structural axis."""
        return tuple(
            self.symbols[qid]
            for qid, quantity in self.quantities.items()
            if quantity.state == "live"
            and self.active_regime in quantity.regime
            and quantity.counting_axis == "discrete-structural"
        )

    def require_admitted(self, relation_id: str) -> Relation:
        try:
            decision = self.admission_decisions[relation_id]
        except KeyError as exc:
            raise AdmissionError(f"unknown relation_id: {relation_id}") from exc
        if not decision.admitted:
            raise AdmissionError(f"{relation_id} refused: {decision.reason}")
        return self.relations[relation_id]

    def constraint_dimension(
        self,
        constraints: Iterable[sp.Expr] | None = None,
        variables: Sequence[sp.Symbol] | None = None,
    ) -> int:
        """Return algebraic dimension from a grevlex initial-monomial ideal."""
        selected_variables = tuple(self.active_variables if variables is None else variables)
        selected_constraints = tuple(
            sp.simplify(expression)
            for expression in (
                self.admitted_constraint_set if constraints is None else tuple(constraints)
            )
            if sp.simplify(expression) != 0
        )
        if not selected_constraints:
            return len(selected_variables)
        for expression in selected_constraints:
            if not expression.free_symbols and expression != 0:
                return -1

        polynomials: list[sp.Expr] = []
        for expression in selected_constraints:
            if not expression.free_symbols <= set(selected_variables):
                outside = sorted(map(str, expression.free_symbols - set(selected_variables)))
                raise EvaluationError(
                    f"constraint uses symbols outside the dimension variables: {outside!r}"
                )
            equation = sp.together(expression)
            try:
                radical_result = unrad(equation, *selected_variables)
            except (NotImplementedError, TypeError, ValueError) as exc:
                raise EvaluationError(
                    f"could not polynomialize constraint: {expression}"
                ) from exc
            if radical_result is not None:
                equation, covariance = radical_result
                if covariance:
                    raise EvaluationError(
                        "constraint needs auxiliary radical-cover variables, which "
                        "the finite-scalar dimension helper does not support"
                    )
            numerator = sp.expand(sp.together(equation).as_numer_denom()[0])
            try:
                polynomial = sp.Poly(
                    numerator,
                    *selected_variables,
                    extension=True,
                ).as_expr()
            except (sp.PolynomialError, TypeError, ValueError) as exc:
                raise EvaluationError(
                    f"constraint is not polynomial after exact radical clearing: {expression}"
                ) from exc
            if polynomial != 0:
                polynomials.append(polynomial)

        if not polynomials:
            return len(selected_variables)
        try:
            basis = sp.groebner(
                tuple(polynomials),
                *selected_variables,
                order="grevlex",
                extension=True,
            )
        except (sp.PolynomialError, TypeError, ValueError) as exc:
            raise EvaluationError(
                "could not construct the grevlex Groebner basis"
            ) from exc

        leading_supports: list[frozenset[int]] = []
        for polynomial in basis.polys:
            leading_monomial = tuple(polynomial.LM(order=basis.order))
            support = frozenset(
                index
                for index, exponent in enumerate(leading_monomial)
                if exponent
            )
            if not support:
                return -1
            leading_supports.append(support)

        # A monomial quotient's dimension is the largest coordinate subset
        # containing no leading-monomial generator support in full.
        for size in range(len(selected_variables), -1, -1):
            for candidate in combinations(range(len(selected_variables)), size):
                candidate_set = frozenset(candidate)
                if all(
                    not support <= candidate_set
                    for support in leading_supports
                ):
                    return size
        raise AssertionError("initial-ideal dimension search failed")

    def certify_positive_real_dimension(
        self,
        constraints: Iterable[sp.Expr] | None = None,
        variables: Sequence[sp.Symbol] | None = None,
        dimension: int | None = None,
    ) -> None:
        """Certify an exact positive smooth point of the algebraic dimension."""
        selected_variables = tuple(self.active_variables if variables is None else variables)
        selected_constraints = tuple(
            sp.simplify(expression)
            for expression in (
                self.admitted_constraint_set if constraints is None else tuple(constraints)
            )
            if sp.simplify(expression) != 0
        )
        expected_dimension = (
            self.constraint_dimension(selected_constraints, selected_variables)
            if dimension is None
            else dimension
        )
        if expected_dimension < 0:
            raise EvaluationError("the algebraic constraint locus is empty")
        if not selected_constraints:
            return

        try:
            solution_branches = sp.solve(
                selected_constraints,
                selected_variables,
                dict=True,
                simplify=True,
                check=True,
            )
        except (NotImplementedError, TypeError, ValueError) as exc:
            raise EvaluationError(
                "could not construct an exact constraint-satisfying witness"
            ) from exc
        if isinstance(solution_branches, dict):
            solution_branches = [solution_branches]
        if not solution_branches:
            raise EvaluationError(
                "no exact constraint-satisfying witness was found; "
                "the locus may be empty or unsupported"
            )

        jacobian = sp.Matrix(selected_constraints).jacobian(selected_variables)
        witnessed_dimensions: set[int] = set()
        for branch in solution_branches:
            branch_symbols = set(selected_variables) | set().union(
                *(expression.free_symbols for expression in branch.values()),
                set(),
            )
            free_symbols = sorted(branch_symbols - set(branch), key=str)
            for offset in range(5):
                witness = {
                    symbol: sp.Integer(sp.prime(index + offset + 1))
                    for index, symbol in enumerate(free_symbols)
                }
                pending = dict(branch)
                while pending:
                    progress = False
                    for symbol, expression in tuple(pending.items()):
                        value = sp.simplify(expression.subs(witness))
                        if value.free_symbols:
                            continue
                        witness[symbol] = value
                        del pending[symbol]
                        progress = True
                    if not progress:
                        break
                if pending or not set(selected_variables) <= set(witness):
                    continue
                if any(
                    sp.simplify(witness[symbol]).is_positive is not True
                    for symbol in selected_variables
                ):
                    continue
                if any(
                    sp.simplify(expression.subs(witness)) != 0
                    for expression in selected_constraints
                ):
                    continue
                exact_jacobian = jacobian.subs(witness)
                if exact_jacobian.free_symbols:
                    continue
                rank = int(exact_jacobian.rank())
                local_dimension = len(selected_variables) - rank
                witnessed_dimensions.add(local_dimension)
                if local_dimension == expected_dimension:
                    return

        if not witnessed_dimensions:
            raise EvaluationError(
                "solver branches yielded no exact positive constraint-satisfying witness"
            )
        raise EvaluationError(
            "exact positive witnesses did not attain the algebraic dimension: "
            f"expected={expected_dimension}, witnessed={sorted(witnessed_dimensions)!r}"
        )

    @staticmethod
    def _predicate_holds(predicate: str, value: sp.Expr) -> bool:
        simplified = sp.simplify(value)
        if simplified.free_symbols or simplified.is_real is not True:
            return False
        if predicate == "real":
            return True
        if predicate == "positive":
            return simplified.is_positive is True
        if predicate == "nonnegative":
            return simplified.is_nonnegative is True
        if predicate == "nonzero":
            return simplified.is_zero is False
        raise EvaluationError(f"unknown assumption predicate: {predicate}")

    def evaluate_output(
        self,
        output: str,
        numeric_inputs: Mapping[str, Any],
        *,
        allow_declared_defaults: bool = False,
    ) -> sp.Expr:
        """Recursively compute an admitted designated output from primitive inputs.

        A caller may not supply any output that lies on the admitted dependency
        path: those values are recomputed, which prevents independently freezing
        the very quantity the relation claims to derive.
        """
        target = self.resolve_qid(output)
        if target not in self._admitted_by_output:
            relation_ids = [
                relation.relation_id
                for relation in self.relations.values()
                if relation.designated_output == target
            ]
            detail = (
                self.admission_decisions[relation_ids[0]].reason
                if relation_ids
                else "no relation designates this output"
            )
            raise AdmissionError(f"{target} is not an admitted output: {detail}")

        provided: dict[str, sp.Expr] = {}
        for name, raw_value in numeric_inputs.items():
            qid = self.resolve_qid(str(name))
            if qid in provided:
                raise EvaluationError(f"numeric inputs duplicate QID through aliases: {qid}")
            value = sp.sympify(raw_value)
            if value.free_symbols or value.is_number is not True:
                raise EvaluationError(f"{qid}: input is not numeric: {raw_value!r}")
            provided[qid] = value

        memo: dict[str, sp.Expr] = {}
        used_inputs: set[str] = set()
        defaulted_inputs: set[str] = set()
        visiting: list[str] = []

        def compute(qid: str) -> sp.Expr:
            if qid in memo:
                return memo[qid]
            relation = self._admitted_by_output.get(qid)
            if relation is None:
                if qid in provided:
                    used_inputs.add(qid)
                    memo[qid] = provided[qid]
                    return memo[qid]
                declared_value = self.quantities[qid].value
                if declared_value is None or not allow_declared_defaults:
                    opt_in_hint = (
                        "; declared value available with allow_declared_defaults=True"
                        if declared_value is not None
                        else ""
                    )
                    raise EvaluationError(
                        f"missing primitive/residue numeric input: {qid}{opt_in_hint}"
                    )
                used_inputs.add(qid)
                defaulted_inputs.add(qid)
                memo[qid] = declared_value
                return memo[qid]
            if qid in provided:
                raise EvaluationError(
                    f"{qid} is derived on this path and must be computed, not independently frozen"
                )
            if qid in visiting:
                raise EvaluationError(f"dataflow cycle requires implicit-block evaluation: {visiting + [qid]}")
            if relation.rhs is None:
                raise EvaluationError(f"{relation.relation_id}: no explicit RHS")

            visiting.append(qid)
            values = {input_qid: compute(input_qid) for input_qid in relation.input_qids}
            substitutions = {self.symbols[key]: value for key, value in values.items()}
            for guard_node in relation.raw["denominator_guards"]:
                guard = sp.simplify(self.parse_expression(guard_node).subs(substitutions))
                if guard.is_zero is not False:
                    raise EvaluationError(
                        f"{relation.relation_id}: denominator guard is not provably nonzero: {guard}"
                    )
            value = sp.simplify(relation.rhs.subs(substitutions))
            if value.free_symbols or value.is_number is not True:
                raise EvaluationError(f"{relation.relation_id}: RHS did not evaluate numerically: {value}")
            all_values = dict(values)
            all_values[qid] = value
            for assumption in relation.raw["assumptions"]:
                predicate = str(assumption["predicate"])
                for assumption_qid in assumption["qids"]:
                    canonical = self.resolve_qid(str(assumption_qid))
                    if canonical in all_values and not self._predicate_holds(
                        predicate, all_values[canonical]
                    ):
                        raise EvaluationError(
                            f"{relation.relation_id}: {predicate} assumption fails for "
                            f"{canonical}={all_values[canonical]}"
                        )
            residual = sp.simplify(
                relation.residual.subs(
                    {
                        **substitutions,
                        self.symbols[qid]: value,
                    }
                )
            )
            if residual != 0:
                raise EvaluationError(
                    f"{relation.relation_id}: computed output violates residual: {residual}"
                )
            visiting.pop()
            memo[qid] = value
            return value

        result = compute(target)
        unused = set(provided) - used_inputs
        if unused:
            raise EvaluationError(f"unused numeric inputs would mask dataflow mistakes: {sorted(unused)!r}")
        if defaulted_inputs:
            rendered_defaults = ", ".join(
                f"{qid}={memo[qid]}" for qid in sorted(defaulted_inputs)
            )
            warnings.warn(
                f"declared value defaults used: {rendered_defaults}",
                DeclaredValueDefaultWarning,
                stacklevel=2,
            )
        return result


def load_registry(
    directory: Path | str = HERE,
    *,
    active_regime: str | None = None,
) -> Registry:
    quantities, relations, schema = load_raw_documents(directory)
    base = Path(directory)
    ledger_root = base.parent
    return Registry.from_documents(
        quantities,
        relations,
        schema,
        ledger_root=ledger_root,
        active_regime=active_regime,
    )


def main() -> int:
    registry = load_registry()
    print("REGISTRY_SCHEMA: VALID")
    print("ACTIVE_REGIME:", registry.active_regime)
    print(
        "ADMITTED_CONSTRAINTS:",
        ",".join(relation.relation_id for relation in registry.admitted_relations),
    )
    for relation_id, decision in registry.admission_decisions.items():
        if not decision.admitted:
            print(f"REFUSED: {relation_id}: {decision.reason}")
    print("RESIDUE:", ",".join(registry.residue))
    print(
        f"CONVENTION_ORBIT: dimension={len(registry.convention_variables)} "
        f"coordinates={','.join(str(variable) for variable in registry.convention_variables)}"
    )
    print(
        "DISCRETE_STRUCTURAL: choices="
        + ",".join(str(variable) for variable in registry.discrete_structural_variables)
    )
    after = registry.constraint_dimension()
    registry.certify_positive_real_dimension(dimension=after)
    print(
        f"DIMENSION: ambient={len(registry.active_variables)} "
        f"after={after} "
        f"rank={len(registry.active_variables) - after}"
    )
    propagated = registry.evaluate_output(
        "lambda_gamma",
        {"K": 1, "rho0": 1, "mass": 5, "c_gamma": 1},
    )
    print("EVALUATION: lambda_gamma=", propagated, sep="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
