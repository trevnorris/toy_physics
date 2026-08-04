#!/usr/bin/env python3
"""Checks for tagged stdout emitted by symbolic engine audit scripts.

This module deliberately contains no expected roots, speeds, or other physics
results.  Comparison targets come from another emission, a control emission,
or the reduction registry at the time a check runs.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import sympy as sp
import yaml

try:
    from .registry_read import Registry, load_registry
except ImportError:  # Direct ``python engine_output_checks.py`` invocation.
    from registry_read import Registry, load_registry


TAG_PATTERN = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)[ \t]*(?::|->)[ \t]*(.*)$")
CONTROL_TAG_PATTERN = re.compile(r"^(?P<base>.+)_X(?P<control>[1-9][0-9]*)_(?P<suffix>.+)$")
TOKEN_PATTERN = re.compile(
    r"(?P<WS>\s+)|(?P<ARROW>->)|(?P<EQUAL>==)|(?P<INT>[0-9]+)|"
    r"(?P<IDENT>[A-Za-z][A-Za-z0-9_$]*)|(?P<PUNCT>[+\-*/^{}()[\],])"
)
ZERO_DIMENSION = (sp.Integer(0), sp.Integer(0), sp.Integer(0))


class HarnessError(RuntimeError):
    """An operational error that prevents a requested check from running."""


class CasParseError(ValueError):
    """Strict Mathematica InputForm parsing failed."""


class DimensionError(ValueError):
    """A parsed construct has no implemented dimensional rule."""


@dataclass(frozen=True)
class Unparsed:
    raw: str
    error: str


class OpaqueDerivative(sp.Symbol):
    """A Mathematica Derivative application kept as one atomic symbol."""

    def __new__(
        cls,
        name: str,
        orders: tuple[int, ...],
        function_name: str,
        variables: tuple[str, ...],
    ) -> "OpaqueDerivative":
        obj = sp.Symbol.__new__(cls, name)
        obj.orders = orders
        obj.function_name = function_name
        obj.variables = variables
        return obj


@dataclass(frozen=True)
class TermDimension:
    expression: str
    dimension: tuple[sp.Expr, sp.Expr, sp.Expr] | None


@dataclass(frozen=True)
class NonHomogeneous:
    tag: str
    expression: str
    summands: tuple[TermDimension, ...]


@dataclass(frozen=True)
class UnknownSymbol:
    tag: str
    symbol: str


@dataclass(frozen=True)
class DimensionReport:
    total_tags: int
    checked_tags: int
    homogeneous_tags: int
    non_homogeneous: tuple[NonHomogeneous, ...]
    unknown_symbols: tuple[UnknownSymbol, ...]
    errors: tuple[str, ...]
    unparsed: tuple[str, ...]


@dataclass(frozen=True)
class ControlRow:
    tag: str
    status: str
    control_tags: tuple[str, ...]
    differing_controls: tuple[str, ...]


@dataclass(frozen=True)
class ControlReport:
    controls: tuple[str, ...]
    responsive: tuple[ControlRow, ...]
    invariant: tuple[ControlRow, ...]
    unparsed: tuple[ControlRow, ...]
    unpaired: tuple[str, ...]


@dataclass(frozen=True)
class CrossEngineRow:
    quantity: str
    status: str
    tags: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class CrossEngineReport:
    rows: tuple[CrossEngineRow, ...]


@dataclass(frozen=True)
class RegistryResidualRow:
    relation_id: str
    status: str
    residual: sp.Expr | None
    detail: str = ""


@dataclass(frozen=True)
class RegistryResidualReport:
    rows: tuple[RegistryResidualRow, ...]


@dataclass(frozen=True)
class HarnessReport:
    controls: ControlReport
    dimensions: DimensionReport
    cross_engine: CrossEngineReport
    registry: RegistryResidualReport
    unparsed: tuple[tuple[str, str, Unparsed], ...]

    @property
    def operational_failures(self) -> tuple[str, ...]:
        failures: list[str] = []
        if self.unparsed:
            failures.append(f"{len(self.unparsed)} value(s) were UNPARSED")
        if self.dimensions.unknown_symbols:
            failures.append(
                f"{len(self.dimensions.unknown_symbols)} dimensional UNKNOWN_SYMBOL occurrence(s)"
            )
        failures.extend(self.dimensions.errors)
        failures.extend(
            f"cross-engine {row.quantity}: {row.status}"
            for row in self.cross_engine.rows
            if row.status in {"MISSING", "UNPARSED"}
        )
        failures.extend(
            f"registry {row.relation_id}: {row.status}"
            for row in self.registry.rows
            if row.status in {"MISSING", "UNPARSED"}
        )
        return tuple(failures)


def parse_tagged_output(stdout: str) -> dict[str, str]:
    """Parse the same tagged-record grammar as ``derived_or_declared.py``."""
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


@dataclass(frozen=True)
class _Token:
    kind: str
    text: str
    offset: int


def _tokenize(text: str) -> tuple[_Token, ...]:
    tokens: list[_Token] = []
    offset = 0
    while offset < len(text):
        match = TOKEN_PATTERN.match(text, offset)
        if match is None:
            excerpt = text[offset : offset + 24]
            raise CasParseError(f"unsupported token at offset {offset}: {excerpt!r}")
        kind = match.lastgroup
        assert kind is not None
        if kind != "WS":
            tokens.append(_Token(kind, match.group(), offset))
        offset = match.end()
    if not tokens:
        raise CasParseError("empty CAS value")
    tokens.append(_Token("EOF", "", len(text)))
    return tuple(tokens)


class _MathematicaParser:
    def __init__(self, text: str) -> None:
        self.text = text
        self.tokens = _tokenize(text)
        self.position = 0

    @property
    def token(self) -> _Token:
        return self.tokens[self.position]

    def accept(self, text: str) -> bool:
        if self.token.text == text:
            self.position += 1
            return True
        return False

    def expect(self, text: str) -> None:
        if not self.accept(text):
            raise CasParseError(
                f"expected {text!r} at offset {self.token.offset}, got {self.token.text!r}"
            )

    def parse(self) -> object:
        result = self.parse_relation()
        if self.token.kind != "EOF":
            raise CasParseError(
                f"trailing token at offset {self.token.offset}: {self.token.text!r}"
            )
        return result

    def parse_relation(self) -> object:
        left = self.parse_add()
        if self.accept("->"):
            return (left, self.parse_relation())
        if self.accept("=="):
            right = self.parse_relation()
            return sp.Eq(
                self._scalar(left, "equation left side"),
                self._scalar(right, "equation right side"),
                evaluate=False,
            )
        return left

    def parse_add(self) -> object:
        result = self.parse_mul()
        while self.token.text in {"+", "-"}:
            operator = self.token.text
            self.position += 1
            right = self.parse_mul()
            result = self._binary(operator, result, right)
        return result

    def parse_mul(self) -> object:
        result = self.parse_unary()
        while self.token.text in {"*", "/"}:
            operator = self.token.text
            self.position += 1
            right = self.parse_unary()
            result = self._binary(operator, result, right)
        return result

    def parse_unary(self) -> object:
        if self.accept("+"):
            return self.parse_unary()
        if self.accept("-"):
            operand = self.parse_unary()
            return -self._scalar(operand, "unary minus")
        return self.parse_power()

    def parse_power(self) -> object:
        base = self.parse_primary()
        if self.accept("^"):
            exponent = self.parse_unary()
            return self._scalar(base, "power base") ** self._scalar(exponent, "power exponent")
        return base

    def parse_primary(self) -> object:
        token = self.token
        if token.kind == "INT":
            self.position += 1
            return int(token.text)
        if token.kind == "IDENT":
            self.position += 1
            if token.text == "True":
                return True
            if token.text == "False":
                return False
            if token.text == "E":
                value: object = sp.E
            elif token.text == "I":
                value = sp.I
            elif token.text == "Derivative" and self.token.text == "[":
                return self.parse_derivative()
            else:
                value = sp.Symbol(token.text)
            while self.token.text == "[":
                arguments = self.parse_arguments()
                value = self._function_call(value, arguments)
            return value
        if self.accept("("):
            result = self.parse_relation()
            self.expect(")")
            return result
        if self.accept("{"):
            values: list[object] = []
            if not self.accept("}"):
                while True:
                    values.append(self.parse_relation())
                    if self.accept("}"):
                        break
                    self.expect(",")
            return values
        raise CasParseError(
            f"expected expression at offset {token.offset}, got {token.text!r}"
        )

    def parse_arguments(self) -> tuple[object, ...]:
        self.expect("[")
        arguments: list[object] = []
        if not self.accept("]"):
            while True:
                arguments.append(self.parse_relation())
                if self.accept("]"):
                    break
                self.expect(",")
        return tuple(arguments)

    def parse_derivative(self) -> OpaqueDerivative:
        orders_raw = self.parse_arguments()
        function_raw = self.parse_arguments()
        variables_raw = self.parse_arguments()
        if (
            len(orders_raw) != 4
            or any(isinstance(value, bool) or not isinstance(value, int) for value in orders_raw)
            or len(function_raw) != 1
            or not isinstance(function_raw[0], sp.Symbol)
            or len(variables_raw) != 4
            or any(not isinstance(value, sp.Symbol) for value in variables_raw)
        ):
            raise CasParseError(
                "Derivative must have Derivative[i,j,k,l][f][t,x,y,z] form"
            )
        orders = tuple(int(value) for value in orders_raw)
        function_name = str(function_raw[0])
        variables = tuple(str(value) for value in variables_raw)
        rendered = (
            f"Derivative[{', '.join(map(str, orders))}]"
            f"[{function_name}][{', '.join(variables)}]"
        )
        return OpaqueDerivative(rendered, orders, function_name, variables)

    @staticmethod
    def _scalar(value: object, context: str) -> sp.Expr:
        if isinstance(value, bool) or not isinstance(value, (int, sp.Basic)):
            raise CasParseError(f"{context} is not a scalar expression")
        return sp.Integer(value) if isinstance(value, int) else value

    def _binary(self, operator: str, left: object, right: object) -> sp.Expr:
        lhs = self._scalar(left, f"{operator} left operand")
        rhs = self._scalar(right, f"{operator} right operand")
        if operator == "+":
            return lhs + rhs
        if operator == "-":
            return lhs - rhs
        if operator == "*":
            return lhs * rhs
        if operator == "/":
            return lhs / rhs
        raise AssertionError(operator)

    def _function_call(self, head: object, arguments: tuple[object, ...]) -> sp.Expr:
        if not isinstance(head, sp.Symbol):
            raise CasParseError(f"unsupported function head {head!r}")
        scalars = tuple(self._scalar(value, f"argument of {head}") for value in arguments)
        name = str(head)
        if name == "Sqrt":
            if len(scalars) != 1:
                raise CasParseError("Sqrt expects one argument")
            return sp.sqrt(scalars[0])
        if name == "Rational":
            if (
                len(scalars) != 2
                or any(value.is_Integer is not True for value in scalars)
                or scalars[1] == 0
            ):
                raise CasParseError("Rational expects two integers and a nonzero denominator")
            return sp.Rational(int(scalars[0]), int(scalars[1]))
        if name == "Exp":
            if len(scalars) != 1:
                raise CasParseError("Exp expects one argument")
            return sp.exp(scalars[0])
        return sp.Function(name)(*scalars)


def normalize_mathematica(raw: str) -> object | Unparsed:
    """Normalize one InputForm value, preserving failures as data."""
    try:
        return _MathematicaParser(raw).parse()
    except (CasParseError, TypeError, ValueError) as exc:
        return Unparsed(raw=raw, error=str(exc))


def normalize_tags(tags: Mapping[str, str]) -> dict[str, object | Unparsed]:
    return {tag: normalize_mathematica(raw) for tag, raw in tags.items()}


def symbolic_equal(left: object, right: object) -> bool:
    """Compare normalized values symbolically; never compare their text."""
    if isinstance(left, Unparsed) or isinstance(right, Unparsed):
        raise ValueError("UNPARSED values cannot be compared")
    if isinstance(left, bool) or isinstance(right, bool):
        return isinstance(left, bool) and isinstance(right, bool) and left is right
    if isinstance(left, list) or isinstance(right, list):
        return (
            isinstance(left, list)
            and isinstance(right, list)
            and len(left) == len(right)
            and all(symbolic_equal(a, b) for a, b in zip(left, right))
        )
    if isinstance(left, tuple) or isinstance(right, tuple):
        return (
            isinstance(left, tuple)
            and isinstance(right, tuple)
            and len(left) == len(right)
            and all(symbolic_equal(a, b) for a, b in zip(left, right))
        )
    if isinstance(left, sp.Equality) or isinstance(right, sp.Equality):
        return (
            isinstance(left, sp.Equality)
            and isinstance(right, sp.Equality)
            and symbolic_equal(left.lhs, right.lhs)
            and symbolic_equal(left.rhs, right.rhs)
        )
    if isinstance(left, (int, sp.Basic)) and isinstance(right, (int, sp.Basic)):
        lhs = sp.Integer(left) if isinstance(left, int) else left
        rhs = sp.Integer(right) if isinstance(right, int) else right
        try:
            return sp.simplify(lhs - rhs) == 0
        except (TypeError, ValueError):
            return False
    return False


def symbolic_multiset_equal(left: object, right: object) -> bool:
    """Order-insensitive symbolic comparison, used only for root multisets."""
    if not isinstance(left, list) or not isinstance(right, list) or len(left) != len(right):
        return False
    unmatched = list(right)
    for item in left:
        for index, candidate in enumerate(unmatched):
            if symbolic_equal(item, candidate):
                del unmatched[index]
                break
        else:
            return False
    return not unmatched


def _is_root_multiset_tag(tag: str) -> bool:
    return tag.endswith("_ROOT_MULTISET") or tag.endswith("_FULL_ROOT_MULTISET")


def _compare_for_tag(tag: str, left: object, right: object) -> bool:
    if _is_root_multiset_tag(tag):
        return symbolic_multiset_equal(left, right)
    return symbolic_equal(left, right)


def check_control_response(values: Mapping[str, object | Unparsed]) -> ControlReport:
    """Partition fully-counterparted main tags into RESPONSIVE and INVARIANT."""
    controls_by_base: dict[str, set[str]] = {}
    for tag in values:
        match = CONTROL_TAG_PATTERN.match(tag)
        if match:
            controls_by_base.setdefault(match.group("base"), set()).add(match.group("control"))

    responsive: list[ControlRow] = []
    invariant: list[ControlRow] = []
    unparsed: list[ControlRow] = []
    unpaired: list[str] = []
    all_control_labels: list[str] = []

    for base, control_ids_set in sorted(controls_by_base.items()):
        control_ids = tuple(sorted(control_ids_set, key=int))
        all_control_labels.extend(f"{base}_X{control}" for control in control_ids)
        main_prefix = base + "_"
        for main_tag in sorted(values):
            if not main_tag.startswith(main_prefix) or CONTROL_TAG_PATTERN.match(main_tag):
                continue
            suffix = main_tag[len(main_prefix) :]
            control_tags = tuple(f"{base}_X{control}_{suffix}" for control in control_ids)
            if any(tag not in values for tag in control_tags):
                unpaired.append(main_tag)
                continue
            row = ControlRow(main_tag, "", control_tags, ())
            compared = (values[main_tag], *(values[tag] for tag in control_tags))
            if any(isinstance(value, Unparsed) for value in compared):
                unparsed.append(ControlRow(main_tag, "UNPARSED", control_tags, ()))
                continue
            differing = tuple(
                tag
                for tag in control_tags
                if not _compare_for_tag(main_tag, values[main_tag], values[tag])
            )
            if differing:
                responsive.append(ControlRow(main_tag, "RESPONSIVE", control_tags, differing))
            else:
                invariant.append(ControlRow(main_tag, "INVARIANT", control_tags, ()))

    return ControlReport(
        controls=tuple(all_control_labels),
        responsive=tuple(responsive),
        invariant=tuple(invariant),
        unparsed=tuple(unparsed),
        unpaired=tuple(sorted(set(unpaired))),
    )


def _as_dimension(value: object, label: str) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    if not isinstance(value, list) or len(value) != 3:
        raise HarnessError(f"{label} must carry a three-component dimension vector")
    result: list[sp.Expr] = []
    for component in value:
        if isinstance(component, bool) or not isinstance(component, (int, sp.Expr)):
            raise HarnessError(f"{label} has a nonscalar dimension component: {component!r}")
        result.append(sp.Integer(component) if isinstance(component, int) else component)
    return tuple(result)  # type: ignore[return-value]


def _dimension_equal(
    left: tuple[sp.Expr, sp.Expr, sp.Expr],
    right: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> bool:
    return all(sp.simplify(a - b) == 0 for a, b in zip(left, right))


def _dimension_add(
    left: tuple[sp.Expr, sp.Expr, sp.Expr],
    right: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return tuple(sp.simplify(a + b) for a, b in zip(left, right))  # type: ignore[return-value]


def _dimension_scale(
    dimension: tuple[sp.Expr, sp.Expr, sp.Expr], factor: sp.Expr
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return tuple(sp.simplify(factor * value) for value in dimension)  # type: ignore[return-value]


class _DimensionWalker:
    def __init__(
        self,
        dimensions: Mapping[str, tuple[sp.Expr, sp.Expr, sp.Expr]],
        tag: str,
    ) -> None:
        self.dimensions = dimensions
        self.tag = tag
        self.non_homogeneous: list[NonHomogeneous] = []
        self.unknown: set[str] = set()

    def walk_container(self, value: object) -> None:
        if isinstance(value, bool):
            return
        if isinstance(value, list):
            for item in value:
                self.walk_container(item)
            return
        if isinstance(value, tuple):
            if len(value) != 2:
                raise DimensionError("a normalized Rule must contain exactly two sides")
            sides = tuple(
                sp.Integer(side) if isinstance(side, int) and not isinstance(side, bool) else side
                for side in value
            )
            if any(not isinstance(side, sp.Basic) for side in sides):
                raise DimensionError("Rule sides must be scalar expressions")
            left = self.dimension(sides[0])
            right = self.dimension(sides[1])
            if left is not None and right is not None and not (
                sides[0] == 0 or sides[1] == 0 or _dimension_equal(left, right)
            ):
                self._record(f"{sides[0]} -> {sides[1]}", sides, (left, right))
            return
        if isinstance(value, sp.Equality):
            left = self.dimension(value.lhs)
            right = self.dimension(value.rhs)
            if left is not None and right is not None and not (
                value.lhs == 0 or value.rhs == 0 or _dimension_equal(left, right)
            ):
                self._record(value, (value.lhs, value.rhs), (left, right))
            return
        if isinstance(value, (int, sp.Basic)):
            self.dimension(sp.Integer(value) if isinstance(value, int) else value)
            return
        raise DimensionError(f"unsupported normalized value {type(value).__name__}")

    def dimension(self, expression: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr] | None:
        if isinstance(expression, OpaqueDerivative):
            function_dimension = self._named_dimension(expression.function_name)
            variable_dimensions = [self._named_dimension(name) for name in expression.variables]
            if function_dimension is None or any(value is None for value in variable_dimensions):
                return None
            result = function_dimension
            for order, variable_dimension in zip(expression.orders, variable_dimensions):
                assert variable_dimension is not None
                result = _dimension_add(result, _dimension_scale(variable_dimension, -sp.Integer(order)))
            return result
        if isinstance(expression, sp.Symbol):
            return self._named_dimension(str(expression))
        if expression.is_number is True:
            return ZERO_DIMENSION
        if isinstance(expression, sp.Add):
            terms = expression.args
            dimensions = tuple(self.dimension(term) for term in terms)
            known = tuple(
                (term, dimension)
                for term, dimension in zip(terms, dimensions)
                if dimension is not None and term != 0
            )
            if known and any(not _dimension_equal(known[0][1], value[1]) for value in known[1:]):
                self._record(expression, terms, dimensions)
            return known[0][1] if known else (ZERO_DIMENSION if not terms else None)
        if isinstance(expression, sp.Mul):
            result = ZERO_DIMENSION
            dimensions = tuple(self.dimension(factor) for factor in expression.args)
            if any(dimension is None for dimension in dimensions):
                return None
            for dimension in dimensions:
                assert dimension is not None
                result = _dimension_add(result, dimension)
            return result
        if isinstance(expression, sp.Pow):
            base, exponent = expression.args
            exponent_dimension = self.dimension(exponent)
            base_dimension = self.dimension(base)
            if exponent_dimension is not None and not _dimension_equal(exponent_dimension, ZERO_DIMENSION):
                self._record(expression, (exponent,), (exponent_dimension,))
            if exponent.is_Rational is not True:
                raise DimensionError(f"non-rational power exponent in {expression}")
            if base_dimension is None:
                return None
            return _dimension_scale(base_dimension, exponent)
        if expression.func == sp.exp:
            argument_dimension = self.dimension(expression.args[0])
            if argument_dimension is not None and not _dimension_equal(
                argument_dimension, ZERO_DIMENSION
            ):
                self._record(expression, (expression.args[0],), (argument_dimension,))
            return ZERO_DIMENSION
        if expression.is_Function:
            for argument in expression.args:
                self.dimension(argument)
            function_name = getattr(expression.func, "__name__", str(expression.func))
            return self._named_dimension(function_name)
        raise DimensionError(f"no dimension rule for {type(expression).__name__}: {expression}")

    def _named_dimension(self, name: str) -> tuple[sp.Expr, sp.Expr, sp.Expr] | None:
        dimension = self.dimensions.get(name)
        if dimension is None:
            self.unknown.add(name)
        return dimension

    def _record(
        self,
        expression: object,
        terms: Sequence[object],
        dimensions: Sequence[tuple[sp.Expr, sp.Expr, sp.Expr] | None],
    ) -> None:
        self.non_homogeneous.append(
            NonHomogeneous(
                tag=self.tag,
                expression=str(expression),
                summands=tuple(
                    TermDimension(str(term), dimension)
                    for term, dimension in zip(terms, dimensions)
                ),
            )
        )


def _dimension_table(
    config: Mapping[str, Any], values: Mapping[str, object | Unparsed]
) -> dict[str, tuple[sp.Expr, sp.Expr, sp.Expr]]:
    table: dict[str, tuple[sp.Expr, sp.Expr, sp.Expr]] = {}
    primitives = config.get("primitive_dimensions", {})
    if not isinstance(primitives, Mapping):
        raise HarnessError("primitive_dimensions must be a mapping")
    # Narrow exception: these are definitional conventions (coordinate, wave
    # number, frequency, displacement), never a derived physics result.
    for symbol, vector in primitives.items():
        if (
            not isinstance(symbol, str)
            or not isinstance(vector, list)
            or len(vector) != 3
            or any(isinstance(item, bool) or not isinstance(item, int) for item in vector)
        ):
            raise HarnessError(f"invalid definitional primitive dimension for {symbol!r}")
        table[symbol] = tuple(sp.Integer(item) for item in vector)  # type: ignore[assignment]

    dimensionless = config.get("dimensionless", [])
    if not isinstance(dimensionless, list) or any(not isinstance(name, str) for name in dimensionless):
        raise HarnessError("dimensionless must be a list of symbol names")
    for symbol in dimensionless:
        table[symbol] = ZERO_DIMENSION

    derived = config.get("derived_dimensions", {})
    if not isinstance(derived, Mapping):
        raise HarnessError("derived_dimensions must be a symbol-to-tag mapping")
    for symbol, tag in derived.items():
        if not isinstance(symbol, str) or not isinstance(tag, str):
            raise HarnessError("derived_dimensions entries must map names to tag names")
        if tag not in values:
            raise HarnessError(f"derived dimension tag is missing: {tag}")
        value = values[tag]
        if isinstance(value, Unparsed):
            raise HarnessError(f"derived dimension tag is UNPARSED: {tag}: {value.error}")
        table[symbol] = _as_dimension(value, tag)
    return table


def check_dimensions(
    values: Mapping[str, object | Unparsed], config: Mapping[str, Any]
) -> DimensionReport:
    """Check every parsed emission using dimensions generated by this run."""
    dimensions = _dimension_table(config, values)
    non_homogeneous: list[NonHomogeneous] = []
    unknown_symbols: list[UnknownSymbol] = []
    errors: list[str] = []
    unparsed: list[str] = []
    checked = 0
    homogeneous = 0
    for tag, value in values.items():
        if isinstance(value, Unparsed):
            unparsed.append(tag)
            continue
        walker = _DimensionWalker(dimensions, tag)
        try:
            walker.walk_container(value)
        except DimensionError as exc:
            errors.append(f"{tag}: {exc}")
            continue
        checked += 1
        non_homogeneous.extend(walker.non_homogeneous)
        unknown_symbols.extend(UnknownSymbol(tag, name) for name in sorted(walker.unknown))
        if not walker.non_homogeneous and not walker.unknown:
            homogeneous += 1
    return DimensionReport(
        total_tags=len(values),
        checked_tags=checked,
        homogeneous_tags=homogeneous,
        non_homogeneous=tuple(non_homogeneous),
        unknown_symbols=tuple(unknown_symbols),
        errors=tuple(errors),
        unparsed=tuple(unparsed),
    )


def check_cross_engine(
    rows: Sequence[Mapping[str, str]],
    outputs: Mapping[str, Mapping[str, object | Unparsed]],
) -> CrossEngineReport:
    results: list[CrossEngineRow] = []
    for row in rows:
        quantity = row.get("quantity")
        if not isinstance(quantity, str):
            raise HarnessError("each cross_engine row needs a quantity name")
        tag_map = tuple((engine, tag) for engine, tag in row.items() if engine != "quantity")
        if len(tag_map) < 2 or any(not isinstance(tag, str) for _, tag in tag_map):
            raise HarnessError(f"cross_engine {quantity} needs at least two engine tag names")
        if any(engine not in outputs or tag not in outputs[engine] for engine, tag in tag_map):
            status = "MISSING"
        else:
            compared = tuple(outputs[engine][tag] for engine, tag in tag_map)
            if any(isinstance(value, Unparsed) for value in compared):
                status = "UNPARSED"
            else:
                first_tag = tag_map[0][1]
                multiset = all(_is_root_multiset_tag(tag) for _, tag in tag_map)
                comparator = symbolic_multiset_equal if multiset else symbolic_equal
                status = "AGREE" if all(
                    comparator(compared[0], value) for value in compared[1:]
                ) else "DISAGREE"
        results.append(CrossEngineRow(quantity, status, tag_map))
    return CrossEngineReport(tuple(results))


def check_registry_residuals(
    rows: Sequence[Mapping[str, Any]],
    outputs: Mapping[str, Mapping[str, object | Unparsed]],
    registry: Registry,
    default_engine: str,
) -> RegistryResidualReport:
    results: list[RegistryResidualRow] = []
    for row in rows:
        relation_id = row.get("relation_id")
        qid_tags = row.get("qids", row.get("inputs"))
        engine = row.get("engine", default_engine)
        if not isinstance(relation_id, str) or relation_id not in registry.relations:
            raise HarnessError(f"unknown registry relation: {relation_id!r}")
        if not isinstance(qid_tags, Mapping) or any(
            not isinstance(qid, str) or not isinstance(tag, str) for qid, tag in qid_tags.items()
        ):
            raise HarnessError(f"registry {relation_id} qids must map QID names to tag names")
        if not isinstance(engine, str) or engine not in outputs:
            results.append(RegistryResidualRow(relation_id, "MISSING", None, "engine output"))
            continue
        relation = registry.relations[relation_id]
        if relation.residual is None:
            raise HarnessError(f"registry relation {relation_id} has no residual")
        canonical_tags: dict[str, str] = {}
        for qid, tag in qid_tags.items():
            canonical = registry.resolve_qid(qid)
            if canonical in canonical_tags:
                raise HarnessError(f"registry {relation_id} maps QID {canonical} more than once")
            canonical_tags[canonical] = tag
        involved = {
            qid
            for qid, symbol in registry.symbols.items()
            if symbol in relation.residual.free_symbols
        }
        missing_qids = sorted(involved - set(canonical_tags))
        missing_tags = sorted(tag for tag in canonical_tags.values() if tag not in outputs[engine])
        if missing_qids or missing_tags:
            detail = f"qids={missing_qids} tags={missing_tags}"
            results.append(RegistryResidualRow(relation_id, "MISSING", None, detail))
            continue
        values = {qid: outputs[engine][tag] for qid, tag in canonical_tags.items()}
        if any(isinstance(value, Unparsed) for value in values.values()):
            results.append(RegistryResidualRow(relation_id, "UNPARSED", None))
            continue
        if any(isinstance(value, (bool, list, tuple, sp.Equality)) for value in values.values()):
            raise HarnessError(f"registry {relation_id} substitutions must be scalar expressions")
        substitutions = {registry.symbols[qid]: value for qid, value in values.items()}
        residual = sp.simplify(relation.residual.subs(substitutions))
        status = "ZERO" if residual == 0 else "NONZERO"
        results.append(RegistryResidualRow(relation_id, status, residual))
    return RegistryResidualReport(tuple(results))


def load_config(path: Path | str) -> dict[str, Any]:
    config_path = Path(path)
    try:
        with config_path.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise HarnessError(f"cannot load config {config_path}: {exc}") from exc
    if not isinstance(config, dict):
        raise HarnessError("config top level must be a mapping")
    return config


def load_output(path: Path | str) -> tuple[dict[str, str], dict[str, object | Unparsed]]:
    output_path = Path(path)
    try:
        text = output_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise HarnessError(f"cannot read output {output_path}: {exc}") from exc
    try:
        raw = parse_tagged_output(text)
    except ValueError as exc:
        raise HarnessError(f"cannot parse tagged output {output_path}: {exc}") from exc
    return raw, normalize_tags(raw)


def run_checks(
    config: Mapping[str, Any],
    outputs: Mapping[str, Mapping[str, object | Unparsed]],
    *,
    registry_directory: Path | str | None = None,
) -> HarnessReport:
    default_engine = config.get("default_engine", "default")
    if not isinstance(default_engine, str) or default_engine not in outputs:
        raise HarnessError(f"default_engine {default_engine!r} has no loaded output")
    default_values = outputs[default_engine]
    cross_rows = config.get("cross_engine", [])
    registry_rows = config.get("registry_residual", [])
    if not isinstance(cross_rows, list) or not isinstance(registry_rows, list):
        raise HarnessError("cross_engine and registry_residual must be lists")
    registry = load_registry(registry_directory or Path(__file__).resolve().parent)
    all_unparsed = tuple(
        (engine, tag, value)
        for engine, values in outputs.items()
        for tag, value in values.items()
        if isinstance(value, Unparsed)
    )
    return HarnessReport(
        controls=check_control_response(default_values),
        dimensions=check_dimensions(default_values, config),
        cross_engine=check_cross_engine(cross_rows, outputs),
        registry=check_registry_residuals(
            registry_rows, outputs, registry, default_engine
        ),
        unparsed=all_unparsed,
    )


def _status_counts(rows: Iterable[object]) -> str:
    counts: dict[str, int] = {}
    for row in rows:
        status = str(getattr(row, "status"))
        counts[status] = counts.get(status, 0) + 1
    return " ".join(f"{name.lower()}={count}" for name, count in sorted(counts.items())) or "configured=0"


def _format_dimension(dimension: tuple[sp.Expr, sp.Expr, sp.Expr] | None) -> str:
    return "UNKNOWN" if dimension is None else "[" + ", ".join(map(str, dimension)) + "]"


def format_report(report: HarnessReport) -> str:
    control = report.controls
    dimension = report.dimensions
    lines = [
        (
            f"CONTROL_RESPONSE: compared={len(control.responsive) + len(control.invariant)} "
            f"responsive={len(control.responsive)} invariant={len(control.invariant)} "
            f"unparsed={len(control.unparsed)} unpaired={len(control.unpaired)}"
        ),
        (
            f"DIMENSIONS: total={dimension.total_tags} checked={dimension.checked_tags} "
            f"homogeneous={dimension.homogeneous_tags} "
            f"non_homogeneous={len(dimension.non_homogeneous)} "
            f"unknown_symbol={len(dimension.unknown_symbols)} unparsed={len(dimension.unparsed)}"
        ),
        f"CROSS_ENGINE: {_status_counts(report.cross_engine.rows)}",
        f"REGISTRY_RESIDUAL: {_status_counts(report.registry.rows)}",
        f"INVARIANT ({len(control.invariant)}):",
    ]
    lines.extend(f"  {row.tag}" for row in control.invariant)
    lines.append(f"NON_HOMOGENEOUS ({len(dimension.non_homogeneous)}):")
    for issue in dimension.non_homogeneous:
        rendered = "; ".join(
            f"{term.expression} => {_format_dimension(term.dimension)}"
            for term in issue.summands
        )
        lines.append(f"  {issue.tag}: {rendered}")
    lines.append(f"UNPARSED ({len(report.unparsed)}):")
    lines.extend(
        f"  {engine}:{tag}: {value.error}; raw={value.raw!r}"
        for engine, tag, value in report.unparsed
    )
    disagreements = [
        f"cross_engine:{row.quantity}"
        for row in report.cross_engine.rows
        if row.status == "DISAGREE"
    ] + [
        f"registry:{row.relation_id}: residual={row.residual}"
        for row in report.registry.rows
        if row.status == "NONZERO"
    ]
    lines.append(f"DISAGREE ({len(disagreements)}):")
    lines.extend(f"  {value}" for value in disagreements)
    lines.append(
        "LIMITATION: triage only; this run does not establish physical correctness, completeness, or derivation coverage."
    )
    return "\n".join(lines)


def _parse_output_arguments(
    specs: Sequence[str], default_engine: str
) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for spec in specs:
        if "=" in spec:
            engine, raw_path = spec.split("=", 1)
            if not engine or not raw_path:
                raise HarnessError(f"invalid --output specification: {spec!r}")
        elif len(specs) == 1:
            engine, raw_path = default_engine, spec
        else:
            raise HarnessError("multiple --output arguments require ENGINE=PATH form")
        if engine in result:
            raise HarnessError(f"duplicate output engine name: {engine}")
        result[engine] = Path(raw_path)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument(
        "--output",
        required=True,
        action="append",
        metavar="[ENGINE=]PATH",
        help="tagged stdout file; repeat with ENGINE=PATH for cross-engine checks",
    )
    parser.add_argument("--registry-dir", type=Path)
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config)
        default_engine = config.get("default_engine", "default")
        if not isinstance(default_engine, str):
            raise HarnessError("default_engine must be a string")
        paths = _parse_output_arguments(args.output, default_engine)
        outputs = {engine: load_output(path)[1] for engine, path in paths.items()}
        report = run_checks(
            config,
            outputs,
            registry_directory=args.registry_dir,
        )
        print(format_report(report))
        if report.operational_failures:
            print("OPERATIONAL_FAILURE: " + "; ".join(report.operational_failures), file=sys.stderr)
            return 2
        return 0
    except (HarnessError, OSError, ValueError) as exc:
        print(f"OPERATIONAL_FAILURE: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
