#!/usr/bin/env python3
"""Exact bound laws for registry dimension-vector declarations."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Mapping

import sympy as sp


DIMENSION_LAW_LANGUAGE = "dimension-prefix-v1"
_PARAMETER_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")

DimensionVector = tuple[int, ...]
SymbolicDimension = tuple[sp.Expr, ...]


class DimensionLawError(ValueError):
    """A dimension law is malformed, unbound, or cannot be evaluated exactly."""


def _require_arity(operator: str, arguments: list[Any], expected: int | str) -> None:
    valid = (
        len(arguments) >= int(expected[2:])
        if isinstance(expected, str)
        else len(arguments) == expected
    )
    if not valid:
        raise DimensionLawError(
            f"dimension law {operator}: expected arity {expected}, observed {len(arguments)}"
        )


def _parse_component(node: Any, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    if isinstance(node, bool):
        raise DimensionLawError("boolean is not a dimension-law integer")
    if isinstance(node, int):
        return sp.Integer(node)
    if not isinstance(node, list) or not node or not isinstance(node[0], str):
        raise DimensionLawError(f"malformed dimension-law expression: {node!r}")

    operator = node[0]
    arguments = node[1:]
    if operator == "Ref":
        _require_arity(operator, arguments, 1)
        name = arguments[0]
        if not isinstance(name, str) or _PARAMETER_NAME.fullmatch(name) is None:
            raise DimensionLawError(f"invalid dimension-law symbol: {name!r}")
        return symbols.setdefault(name, sp.Symbol(name, integer=True))
    if operator in {"Add", "Mul"}:
        _require_arity(operator, arguments, ">=2")
    elif operator == "Sub":
        _require_arity(operator, arguments, 2)
    elif operator == "Neg":
        _require_arity(operator, arguments, 1)
    else:
        raise DimensionLawError(f"unsupported dimension-law operator {operator!r}")

    parsed = tuple(_parse_component(argument, symbols) for argument in arguments)
    if operator == "Add":
        return sp.Add(*parsed)
    if operator == "Mul":
        return sp.Mul(*parsed)
    if operator == "Sub":
        return parsed[0] - parsed[1]
    return -parsed[0]


@dataclass(frozen=True)
class BoundDimensionLaw:
    """A symbolic dimension vector whose names are bound to structural QIDs."""

    bindings: tuple[tuple[str, str], ...]
    components: SymbolicDimension
    reference_values: tuple[tuple[str, int], ...]

    @property
    def binding_map(self) -> dict[str, str]:
        return dict(self.bindings)

    @property
    def binding_qids(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(qid for _, qid in self.bindings))

    @property
    def canonical_components(self) -> SymbolicDimension:
        substitutions = {
            sp.Symbol(name, integer=True): sp.Symbol(qid, integer=True)
            for name, qid in self.bindings
        }
        return tuple(sp.expand(component.subs(substitutions)) for component in self.components)

    def evaluate_parameters(self, values: Mapping[str, int]) -> DimensionVector:
        expected = {name for name, _ in self.bindings}
        observed = set(values)
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        if missing or extra:
            raise DimensionLawError(
                f"dimension-law parameter values differ: missing={missing!r} extra={extra!r}"
            )
        substitutions: dict[sp.Symbol, sp.Integer] = {}
        for name, _ in self.bindings:
            value = values[name]
            if isinstance(value, bool) or not isinstance(value, int):
                raise DimensionLawError(
                    f"dimension-law parameter {name!r} must be an integer, observed {value!r}"
                )
            substitutions[sp.Symbol(name, integer=True)] = sp.Integer(value)
        return self._evaluate(substitutions)

    def evaluate_bindings(self, values: Mapping[str, int]) -> DimensionVector:
        required = set(self.binding_qids)
        missing = sorted(required - set(values))
        if missing:
            raise DimensionLawError(f"missing dimension-law binding values: {missing!r}")
        parameters = {name: values[qid] for name, qid in self.bindings}
        return self.evaluate_parameters(parameters)

    def evaluate_reference(self) -> DimensionVector:
        return self.evaluate_parameters(dict(self.reference_values))

    def _evaluate(self, substitutions: Mapping[sp.Symbol, sp.Integer]) -> DimensionVector:
        evaluated: list[int] = []
        for component in self.components:
            value = sp.expand(component.subs(substitutions))
            if value.free_symbols or value.is_Integer is not True:
                raise DimensionLawError(
                    f"dimension-law component did not evaluate to an integer: {value}"
                )
            evaluated.append(int(value))
        return tuple(evaluated)


@dataclass(frozen=True)
class ParsedDimensionDeclaration:
    """The primary declaration plus its explicit numeric reference vector."""

    declaration: DimensionVector | BoundDimensionLaw
    reference_vector: DimensionVector


def parse_dimension_declaration(
    raw: Mapping[str, Any], component_count: int
) -> ParsedDimensionDeclaration:
    exponents = raw.get("exponents")
    if (
        not isinstance(exponents, list)
        or len(exponents) != component_count
        or any(isinstance(value, bool) or not isinstance(value, int) for value in exponents)
    ):
        raise DimensionLawError(
            f"dimension exponents must be {component_count} integers"
        )
    reference_vector = tuple(exponents)
    raw_law = raw.get("law")
    if raw_law is None:
        return ParsedDimensionDeclaration(reference_vector, reference_vector)
    if not isinstance(raw_law, Mapping):
        raise DimensionLawError("dimension law must be a mapping")
    if raw_law.get("expression_language") != DIMENSION_LAW_LANGUAGE:
        raise DimensionLawError(
            f"unsupported dimension-law language {raw_law.get('expression_language')!r}"
        )

    raw_bindings = raw_law.get("bindings")
    raw_components = raw_law.get("components")
    raw_reference = raw_law.get("reference_values")
    if not isinstance(raw_bindings, Mapping) or not raw_bindings:
        raise DimensionLawError("dimension law bindings must be a nonempty mapping")
    if not isinstance(raw_components, list) or len(raw_components) != component_count:
        raise DimensionLawError(
            f"dimension law components must contain {component_count} expressions"
        )
    if not isinstance(raw_reference, Mapping):
        raise DimensionLawError("dimension law reference_values must be a mapping")

    bindings: list[tuple[str, str]] = []
    for name, qid in raw_bindings.items():
        if not isinstance(name, str) or _PARAMETER_NAME.fullmatch(name) is None:
            raise DimensionLawError(f"invalid dimension-law binding name: {name!r}")
        if not isinstance(qid, str):
            raise DimensionLawError(f"dimension-law binding {name!r} has a non-string QID")
        bindings.append((name, qid))

    symbols: dict[str, sp.Symbol] = {}
    components = tuple(_parse_component(node, symbols) for node in raw_components)
    declared_names = set(raw_bindings)
    referenced_names = set(symbols)
    unbound = sorted(referenced_names - declared_names)
    unused = sorted(declared_names - referenced_names)
    if unbound:
        raise DimensionLawError(f"unbound dimension-law symbol(s): {unbound!r}")
    if unused:
        raise DimensionLawError(f"unused dimension-law binding(s): {unused!r}")

    reference_values: list[tuple[str, int]] = []
    missing_reference = sorted(declared_names - set(raw_reference))
    extra_reference = sorted(set(raw_reference) - declared_names)
    if missing_reference or extra_reference:
        raise DimensionLawError(
            "dimension-law reference values differ from bindings: "
            f"missing={missing_reference!r} extra={extra_reference!r}"
        )
    for name in raw_bindings:
        value = raw_reference[name]
        if isinstance(value, bool) or not isinstance(value, int):
            raise DimensionLawError(
                f"dimension-law reference value {name!r} must be an integer"
            )
        reference_values.append((str(name), value))

    law = BoundDimensionLaw(tuple(bindings), components, tuple(reference_values))
    observed_reference = law.evaluate_reference()
    if observed_reference != reference_vector:
        raise DimensionLawError(
            "dimension-law reference evaluation differs from exponents: "
            f"law={list(observed_reference)} exponents={list(reference_vector)}"
        )
    return ParsedDimensionDeclaration(law, reference_vector)


def as_symbolic_dimension(
    declaration: DimensionVector | BoundDimensionLaw,
) -> SymbolicDimension:
    if isinstance(declaration, BoundDimensionLaw):
        return declaration.canonical_components
    return tuple(sp.Integer(value) for value in declaration)


def dimension_residual(
    left: SymbolicDimension, right: SymbolicDimension
) -> SymbolicDimension:
    if len(left) != len(right):
        raise DimensionLawError(
            f"dimension lengths differ: left={len(left)} right={len(right)}"
        )
    return tuple(sp.simplify(a - b) for a, b in zip(left, right))


def dimensions_equal(left: SymbolicDimension, right: SymbolicDimension) -> bool:
    return all(value == 0 for value in dimension_residual(left, right))
