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
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import sympy as sp
import yaml

try:
    from .registry_read import Registry, load_registry
except ImportError:  # Direct ``python engine_output_checks.py`` invocation.
    from registry_read import Registry, load_registry


TAG_PATTERN = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)[ \t]*(?::|->)[ \t]*(.*)$")
TAG_NAME_PATTERN = re.compile(r"^(?:WL|PY)_S[1-9][0-9]*_[A-Z][A-Z0-9_]*$")
CONTROL_TAG_PATTERN = re.compile(r"^(?P<base>.+)_X(?P<control>[1-9][0-9]*)_(?P<suffix>.+)$")
TOKEN_PATTERN = re.compile(
    r"(?P<WS>\s+)|(?P<ASSOC_OPEN><\|)|(?P<ASSOC_CLOSE>\|>)|"
    r"(?P<AND>&&)|(?P<OR>\|\|)|(?P<ARROW>->)|(?P<GE>>=)|(?P<LE><=)|"
    r"(?P<EQUAL>==)|(?P<NE>!=)|(?P<GT>>)|(?P<LT><)|(?P<NOT>!)|"
    r"(?P<INT>[0-9]+)|(?P<IDENT>[A-Za-z][A-Za-z0-9_$]*)|"
    r"(?P<PUNCT>[+\-*/^{}()[\],])"
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


class ValueKind(str, Enum):
    EXPRESSION = "expression"
    MATRIX = "matrix"
    LIST = "list"
    MAPPING = "mapping"
    BOOLEAN = "boolean"
    RELATION = "relation"
    INTEGER = "integer"


DIMENSIONFUL_KINDS = frozenset(
    {
        ValueKind.EXPRESSION,
        ValueKind.MATRIX,
        ValueKind.LIST,
        ValueKind.RELATION,
        ValueKind.INTEGER,
    }
)


@dataclass(frozen=True)
class ParsedValue:
    kind: ValueKind
    value: object


@dataclass(frozen=True)
class CasRule:
    left: object
    right: object


def _is_boolean_value(value: object) -> bool:
    return (
        isinstance(value, bool)
        or getattr(value, "is_Boolean", False) is True
        or isinstance(value, sp.assumptions.assume.AppliedPredicate)
        or isinstance(value, (sp.core.relational.Relational, sp.Contains))
    )


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
class NonDimensional:
    tag: str
    kind: ValueKind


@dataclass(frozen=True)
class DimensionReport:
    total_tags: int
    checked_tags: int
    homogeneous_tags: int
    vacuous_tags: tuple[str, ...]
    non_homogeneous: tuple[NonHomogeneous, ...]
    unknown_symbols: tuple[UnknownSymbol, ...]
    errors: tuple[str, ...]
    unparsed: tuple[str, ...]
    non_dimensional: tuple[NonDimensional, ...]


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
    uncovered: tuple[str, ...] = ()
    no_ablation_dimensions: tuple[int, ...] = ()
    missing_declared_packages: tuple[str, ...] = ()
    coverage_required: bool = False


@dataclass(frozen=True)
class ParityRow:
    engine: str
    package: str
    dimension: int | None
    missing: tuple[str, ...]
    extra: tuple[str, ...]


@dataclass(frozen=True)
class ParityExclusionRow:
    engine: str
    tags: tuple[str, ...]
    by_pattern: tuple[tuple[str, int], ...]


@dataclass(frozen=True)
class ParityReport:
    rows: tuple[ParityRow, ...]
    exclusions: tuple[ParityExclusionRow, ...] = ()
    exclusion_warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class CrossEngineTagSetRow:
    left_engine: str
    right_engine: str
    left_only: tuple[str, ...]
    right_only: tuple[str, ...]


@dataclass(frozen=True)
class CrossEngineTagSetReport:
    rows: tuple[CrossEngineTagSetRow, ...]


@dataclass(frozen=True)
class CrossEngineRow:
    quantity: str
    status: str
    tags: tuple[tuple[str, str], ...]
    operands: tuple[tuple[str, str], ...] = ()
    detail: str = ""


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
    parity: ParityReport
    cross_engine_tags: CrossEngineTagSetReport
    dimensions: DimensionReport
    cross_engine: CrossEngineReport
    registry: RegistryResidualReport
    unparsed: tuple[tuple[str, str, Unparsed], ...]

    @property
    def operational_failures(self) -> tuple[str, ...]:
        failures: list[str] = []
        compared = len(self.controls.responsive) + len(self.controls.invariant)
        if self.controls.missing_declared_packages:
            failures.append(
                "declared packages matched no tags: "
                + ", ".join(self.controls.missing_declared_packages)
            )
        if self.controls.coverage_required and self.controls.uncovered:
            dimensions = ",".join(
                f"D{dimension}" for dimension in self.controls.no_ablation_dimensions
            ) or "none"
            failures.append(
                f"control coverage uncovered {len(self.controls.uncovered)} main tag(s); "
                f"dimensions with no ablation={dimensions}"
            )
        if compared >= 10 and not self.controls.responsive:
            failures.append(
                f"control response floor failed: compared={compared} responsive=0"
            )
        if self.unparsed:
            failures.append(f"{len(self.unparsed)} value(s) were UNPARSED")
        if self.dimensions.unknown_symbols:
            failures.extend(
                f"dimensional UNKNOWN_SYMBOL {issue.tag}:{issue.symbol}"
                for issue in self.dimensions.unknown_symbols
            )
        failures.extend(self.dimensions.errors)
        failures.extend(
            f"cross-engine {row.quantity}: {row.status}"
            for row in self.cross_engine.rows
            if row.status in {"MISSING", "UNPARSED", "SHAPE_MISMATCH"}
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
            if not TAG_NAME_PATTERN.fullmatch(match.group(1)):
                raise ValueError(f"rejected invalid tag line: {line!r}")
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
        result = self.parse_or()
        if self.token.kind != "EOF":
            raise CasParseError(
                f"trailing token at offset {self.token.offset}: {self.token.text!r}"
            )
        return result

    def parse_or(self) -> object:
        result = self.parse_and()
        while self.accept("||"):
            result = sp.Or(
                self._boolean(result, "|| left operand"),
                self._boolean(self.parse_and(), "|| right operand"),
                evaluate=False,
            )
        return result

    def parse_and(self) -> object:
        result = self.parse_not()
        while self.accept("&&"):
            result = sp.And(
                self._boolean(result, "&& left operand"),
                self._boolean(self.parse_not(), "&& right operand"),
                evaluate=False,
            )
        return result

    def parse_not(self) -> object:
        if self.accept("!"):
            return sp.Not(
                self._boolean(self.parse_not(), "! operand"), evaluate=False
            )
        return self.parse_relation()

    def parse_relation(self) -> object:
        left = self.parse_add()
        if self.accept("->"):
            return CasRule(left, self.parse_or())
        relations = {
            "==": sp.Eq,
            "!=": sp.Ne,
            ">": sp.Gt,
            ">=": sp.Ge,
            "<": sp.Lt,
            "<=": sp.Le,
        }
        relation = relations.get(self.token.text)
        if relation is not None:
            operator = self.token.text
            self.position += 1
            right = self.parse_add()
            return relation(
                self._scalar(left, f"{operator} left side"),
                self._scalar(right, f"{operator} right side"),
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
            result = self.parse_or()
            self.expect(")")
            return result
        if self.accept("{"):
            values: list[object] = []
            if not self.accept("}"):
                while True:
                    values.append(self.parse_or())
                    if self.accept("}"):
                        break
                    self.expect(",")
            return values
        if self.accept("<|"):
            values: dict[object, object] = {}
            if not self.accept("|>"):
                while True:
                    rule = self.parse_or()
                    if not isinstance(rule, CasRule):
                        raise CasParseError("association entries must have key -> value form")
                    try:
                        if rule.left in values:
                            raise CasParseError(
                                f"duplicate association key {rule.left!r}"
                            )
                        values[rule.left] = rule.right
                    except TypeError as exc:
                        raise CasParseError("association keys must be hashable") from exc
                    if self.accept("|>"):
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
                arguments.append(self.parse_or())
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

    @staticmethod
    def _boolean(value: object, context: str) -> object:
        if _is_boolean_value(value):
            return value
        raise CasParseError(f"{context} is not a boolean expression")

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
        name = str(head)
        if name == "Element":
            if (
                len(arguments) != 2
                or not isinstance(arguments[1], sp.Symbol)
                or str(arguments[1]) != "Reals"
            ):
                raise CasParseError("Element currently supports Element[x, Reals]")
            element = self._scalar(arguments[0], "Element first argument")
            return sp.Contains(element, sp.S.Reals, evaluate=False)
        scalars = tuple(self._scalar(value, f"argument of {head}") for value in arguments)
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


_SYMBOL_ALIASES: dict[str, sp.Expr] = {
    "mu_R": sp.Symbol("muR"),
    "rho_br": sp.Symbol("rhoBr"),
    "rho_z": sp.Symbol("rhoZ"),
    "mu_F": sp.Symbol("muF"),
    "mu_G": sp.Symbol("muG"),
    "lambda_rho": sp.Symbol("lambdaRho"),
    "lambda_mu": sp.Symbol("lambdaMu"),
    "lambda_scale": sp.Symbol("lambdaScale"),
    "omega2": sp.Symbol("omega") ** 2,
    "omegaSquared": sp.Symbol("omega") ** 2,
}


def _canonicalize(value: object) -> object:
    """Canonicalize spelling differences without supplying a comparison value."""
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(_canonicalize)
    if isinstance(value, Mapping):
        return {_canonicalize(key): _canonicalize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_canonicalize(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_canonicalize(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return [_canonicalize(item) for item in sorted(value, key=str)]
    if isinstance(value, CasRule):
        return CasRule(_canonicalize(value.left), _canonicalize(value.right))
    if isinstance(value, sp.Basic):
        replacements = {
            symbol: _SYMBOL_ALIASES[str(symbol)]
            for symbol in value.free_symbols
            if str(symbol) in _SYMBOL_ALIASES
        }
        return value.xreplace(replacements)
    return value


def value_kind(value: object) -> ValueKind:
    if isinstance(value, ParsedValue):
        return value.kind
    if isinstance(value, (CasRule, sp.core.relational.Relational, sp.Contains)):
        return ValueKind.RELATION
    if _is_boolean_value(value):
        return ValueKind.BOOLEAN
    if isinstance(value, Mapping):
        return ValueKind.MAPPING
    if isinstance(value, sp.MatrixBase):
        return ValueKind.MATRIX
    if isinstance(value, (list, tuple, set, frozenset)):
        return ValueKind.LIST
    if isinstance(value, int) or isinstance(value, sp.Integer):
        return ValueKind.INTEGER
    if isinstance(value, sp.Basic):
        return ValueKind.EXPRESSION
    raise CasParseError(f"unsupported normalized value {type(value).__name__}")


def _parsed(value: object) -> ParsedValue:
    canonical = _canonicalize(value)
    return ParsedValue(value_kind(canonical), canonical)


def normalize_mathematica(raw: str) -> ParsedValue | Unparsed:
    """Normalize one Mathematica InputForm value, preserving failures as data."""
    try:
        return _parsed(_MathematicaParser(raw).parse())
    except (CasParseError, TypeError, ValueError) as exc:
        return Unparsed(raw=raw, error=str(exc))


def normalize_sympy(raw: str) -> ParsedValue | Unparsed:
    """Normalize one SymPy ``str`` emission, including its container forms."""
    try:
        value = sp.sympify(
            raw,
            locals={
                "Matrix": sp.Matrix,
                "Eq": sp.Eq,
                "Derivative": sp.Derivative,
                "Q": sp.Q,
            },
            evaluate=False,
        )
        return _parsed(value)
    except Exception as exc:  # SymPy exposes several parser exception classes.
        return Unparsed(raw=raw, error=str(exc))


def _normalize_auto(raw: str) -> ParsedValue | Unparsed:
    mathematica = normalize_mathematica(raw)
    return normalize_sympy(raw) if isinstance(mathematica, Unparsed) else mathematica


def normalize_tags(
    tags: Mapping[str, str], syntax: str | None = None
) -> dict[str, ParsedValue | Unparsed]:
    if syntax not in {None, "wl", "py"}:
        raise ValueError(f"unknown CAS syntax {syntax!r}")
    normalized: dict[str, ParsedValue | Unparsed] = {}
    for tag, raw in tags.items():
        selected = syntax
        if selected is None:
            selected = "py" if tag.startswith("PY_") else "wl" if tag.startswith("WL_") else None
        normalizer = (
            normalize_sympy
            if selected == "py"
            else normalize_mathematica
            if selected == "wl"
            else _normalize_auto
        )
        normalized[tag] = normalizer(raw)
    return normalized


def _unwrap(value: object) -> object:
    return value.value if isinstance(value, ParsedValue) else value


def _sequence_form(value: object) -> list[object] | None:
    value = _unwrap(value)
    if isinstance(value, sp.MatrixBase):
        if value.rows == 1 or value.cols == 1:
            return list(value)
        return [[value[row, column] for column in range(value.cols)] for row in range(value.rows)]
    if isinstance(value, (list, tuple)):
        return list(value)
    return None


def symbolic_shapes_compatible(left: object, right: object) -> bool:
    left = _unwrap(left)
    right = _unwrap(right)
    left_sequence = _sequence_form(left)
    right_sequence = _sequence_form(right)
    if left_sequence is not None or right_sequence is not None:
        return (
            left_sequence is not None
            and right_sequence is not None
            and len(left_sequence) == len(right_sequence)
            and all(
                symbolic_shapes_compatible(a, b)
                for a, b in zip(left_sequence, right_sequence)
            )
        )
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        return isinstance(left, Mapping) and isinstance(right, Mapping) and len(left) == len(right)
    if isinstance(left, CasRule) or isinstance(right, CasRule):
        return isinstance(left, CasRule) and isinstance(right, CasRule)
    left_boolean = _is_boolean_value(left)
    right_boolean = _is_boolean_value(right)
    return left_boolean == right_boolean


def symbolic_equal(left: object, right: object) -> bool:
    """Compare normalized values symbolically; never compare their text."""
    if isinstance(left, Unparsed) or isinstance(right, Unparsed):
        raise ValueError("UNPARSED values cannot be compared")
    left = _unwrap(left)
    right = _unwrap(right)
    left_sequence = _sequence_form(left)
    right_sequence = _sequence_form(right)
    if left_sequence is not None or right_sequence is not None:
        return (
            left_sequence is not None
            and right_sequence is not None
            and len(left_sequence) == len(right_sequence)
            and all(symbolic_equal(a, b) for a, b in zip(left_sequence, right_sequence))
        )
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        if not isinstance(left, Mapping) or not isinstance(right, Mapping) or len(left) != len(right):
            return False
        unmatched = list(right.items())
        for left_key, left_value in left.items():
            for index, (right_key, right_value) in enumerate(unmatched):
                if symbolic_equal(left_key, right_key) and symbolic_equal(left_value, right_value):
                    del unmatched[index]
                    break
            else:
                return False
        return not unmatched
    if isinstance(left, CasRule) or isinstance(right, CasRule):
        return (
            isinstance(left, CasRule)
            and isinstance(right, CasRule)
            and symbolic_equal(left.left, right.left)
            and symbolic_equal(left.right, right.right)
        )
    if _is_boolean_value(left) or _is_boolean_value(right):
        try:
            return sp.simplify(sp.Equivalent(left, right)) is sp.true
        except (TypeError, ValueError):
            return left == right
    if isinstance(left, (int, sp.Basic)) and isinstance(right, (int, sp.Basic)):
        lhs = sp.Integer(left) if isinstance(left, int) else left
        rhs = sp.Integer(right) if isinstance(right, int) else right
        try:
            return sp.simplify(lhs - rhs) == 0
        except (TypeError, ValueError):
            return False
    return False


def symbolic_multiset_equal(left: object, right: object) -> bool:
    """Order-insensitive symbolic comparison, enabled explicitly by config."""
    left_sequence = _sequence_form(left)
    right_sequence = _sequence_form(right)
    if (
        left_sequence is None
        or right_sequence is None
        or len(left_sequence) != len(right_sequence)
    ):
        return False
    unmatched = list(right_sequence)
    for item in left_sequence:
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


@dataclass(frozen=True)
class _PackageLayout:
    base: str
    dimension: int | None
    main_package: str
    main: dict[str, str]
    controls: dict[str, dict[str, str]]
    declared: bool

    def package_label(self, package: str) -> str:
        return f"{self.base}_{package}" if self.declared else package


def _package_layout(
    values: Mapping[str, object],
    main_package: str | None = None,
    control_packages: Sequence[str] | None = None,
) -> tuple[_PackageLayout, ...]:
    """Return package layouts, split by dimension for declared S10 packages."""
    if main_package is not None:
        controls = tuple(control_packages or ())
        packages = tuple(dict.fromkeys((main_package, *controls)))
        alternatives = "|".join(
            re.escape(package)
            for package in sorted(packages, key=lambda value: (-len(value), value))
        )
        declared_pattern = re.compile(
            rf"^(?P<base>.+)_(?P<package>{alternatives})_"
            rf"D(?P<dimension>[1-9][0-9]*)_(?P<suffix>.+)$"
        )
        grouped: dict[
            tuple[str, int], dict[str, dict[str, str]]
        ] = {}
        for tag in values:
            if "_LOCAL_" in tag:
                continue
            match = declared_pattern.match(tag)
            if match is None:
                continue
            key = (match.group("base"), int(match.group("dimension")))
            grouped.setdefault(key, {}).setdefault(match.group("package"), {})[
                match.group("suffix")
            ] = tag
        return tuple(
            _PackageLayout(
                base=base,
                dimension=dimension,
                main_package=main_package,
                main=packages_at_dimension.get(main_package, {}),
                controls={
                    package: packages_at_dimension[package]
                    for package in controls
                    if package in packages_at_dimension
                },
                declared=True,
            )
            for (base, dimension), packages_at_dimension in sorted(grouped.items())
        )

    # S9 compatibility: infer the historical X<digits> controls when the
    # configuration does not declare package names.
    bases: set[str] = set()
    for tag in values:
        match = CONTROL_TAG_PATTERN.match(tag)
        if match:
            bases.add(match.group("base"))
    layouts: list[_PackageLayout] = []
    for base in sorted(bases):
        main_prefix = f"{base}_MAIN_" if any(
            tag.startswith(f"{base}_MAIN_") for tag in values
        ) else f"{base}_"
        main_package = main_prefix[:-1]
        main: dict[str, str] = {}
        controls: dict[str, dict[str, str]] = {}
        for tag in values:
            match = CONTROL_TAG_PATTERN.match(tag)
            if match and match.group("base") == base:
                package = f"{base}_X{match.group('control')}"
                controls.setdefault(package, {})[match.group("suffix")] = tag
            elif tag.startswith(main_prefix):
                main[tag[len(main_prefix) :]] = tag
        layouts.append(
            _PackageLayout(base, None, main_package, main, controls, False)
        )
    return tuple(layouts)


def check_control_response(
    values: Mapping[str, ParsedValue | Unparsed],
    main_package: str | None = None,
    control_packages: Sequence[str] | None = None,
) -> ControlReport:
    """Partition fully-counterparted main tags into RESPONSIVE and INVARIANT."""
    responsive: list[ControlRow] = []
    invariant: list[ControlRow] = []
    unparsed: list[ControlRow] = []
    unpaired: list[str] = []
    uncovered: list[str] = []
    all_control_labels: set[str] = set()
    no_ablation_dimensions: set[int] = set()
    observed_declared_packages: set[str] = set()
    declared_packages = (
        tuple(dict.fromkeys((main_package, *(control_packages or ()))))
        if main_package is not None
        else ()
    )

    layouts = _package_layout(values, main_package, control_packages)
    for layout in layouts:
        if layout.main:
            observed_declared_packages.add(layout.main_package)
        observed_declared_packages.update(layout.controls)
        all_control_labels.update(
            layout.package_label(package) for package in layout.controls
        )
        if layout.main and not layout.controls:
            uncovered.extend(layout.main.values())
            if layout.dimension is not None:
                no_ablation_dimensions.add(layout.dimension)
            continue
        for suffix, main_tag in sorted(layout.main.items()):
            control_tags = tuple(
                tags[suffix]
                for _package, tags in sorted(layout.controls.items())
                if suffix in tags
            )
            if len(control_tags) != len(layout.controls):
                unpaired.append(main_tag)
                uncovered.append(main_tag)
                continue
            compared = (values[main_tag], *(values[tag] for tag in control_tags))
            if any(isinstance(value, Unparsed) for value in compared):
                unparsed.append(ControlRow(main_tag, "UNPARSED", control_tags, ()))
                uncovered.append(main_tag)
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
        controls=tuple(sorted(all_control_labels)),
        responsive=tuple(responsive),
        invariant=tuple(invariant),
        unparsed=tuple(unparsed),
        unpaired=tuple(sorted(set(unpaired))),
        uncovered=tuple(sorted(set(uncovered))),
        no_ablation_dimensions=tuple(sorted(no_ablation_dimensions)),
        missing_declared_packages=tuple(
            package
            for package in declared_packages
            if package not in observed_declared_packages
        ),
        coverage_required=main_package is not None,
    )


def check_tag_parity(
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
    parity_exclude: Sequence[str] = (),
    main_package: str | None = None,
    control_packages: Sequence[str] | None = None,
) -> ParityReport:
    """Compare each control package's suffix set with its engine's main package."""
    if isinstance(parity_exclude, (str, bytes)):
        raise HarnessError("parity_exclude must be a sequence of substrings")
    configured_patterns = tuple(parity_exclude)
    if any(
        not isinstance(pattern, str) or not pattern
        for pattern in configured_patterns
    ):
        raise HarnessError("parity_exclude must contain only non-empty strings")
    patterns = tuple(dict.fromkeys(configured_patterns))
    rows: list[ParityRow] = []
    exclusions: list[ParityExclusionRow] = []
    exclusion_warnings: list[str] = []
    for engine, values in sorted(outputs.items()):
        layouts = _package_layout(values, main_package, control_packages)
        excluded_tags = {
            tag for tag in values if any(pattern in tag for pattern in patterns)
        }
        if patterns:
            pattern_counts = tuple(
                (pattern, sum(pattern in tag for tag in values))
                for pattern in patterns
            )
            exclusions.append(
                ParityExclusionRow(
                    engine=engine,
                    tags=tuple(sorted(excluded_tags)),
                    by_pattern=pattern_counts,
                )
            )
            exclusion_warnings.extend(
                f"{engine}: configured pattern {pattern!r} matched no tags"
                for pattern, count in pattern_counts
                if count == 0
            )
        for layout in layouts:
            main_suffixes = {
                suffix
                for suffix, tag in layout.main.items()
                if tag not in excluded_tags
            }
            for package, package_tags in sorted(layout.controls.items()):
                package_suffixes = {
                    suffix
                    for suffix, tag in package_tags.items()
                    if tag not in excluded_tags
                }
                rows.append(
                    ParityRow(
                        engine=engine,
                        package=layout.package_label(package),
                        dimension=layout.dimension,
                        missing=tuple(sorted(main_suffixes - package_suffixes)),
                        extra=tuple(sorted(package_suffixes - main_suffixes)),
                    )
                )
    return ParityReport(
        tuple(rows), tuple(exclusions), tuple(exclusion_warnings)
    )


def check_cross_engine_tag_sets(
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
    parity_exclude: Sequence[str] = (),
) -> CrossEngineTagSetReport:
    """Compare complete non-local tag-name sets after removing engine prefixes."""
    canonical: dict[str, dict[str, str]] = {}
    for engine, values in sorted(outputs.items()):
        canonical[engine] = {
            tag.split("_", 1)[1] if "_" in tag else tag: tag
            for tag in values
            if not any(pattern in tag for pattern in parity_exclude)
        }
    rows: list[CrossEngineTagSetRow] = []
    engines = sorted(canonical)
    for left_index, left_engine in enumerate(engines):
        for right_engine in engines[left_index + 1 :]:
            left = canonical[left_engine]
            right = canonical[right_engine]
            rows.append(
                CrossEngineTagSetRow(
                    left_engine,
                    right_engine,
                    tuple(left[name] for name in sorted(left.keys() - right.keys())),
                    tuple(right[name] for name in sorted(right.keys() - left.keys())),
                )
            )
    return CrossEngineTagSetReport(tuple(rows))


def _as_dimension(value: object, label: str) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    value = _unwrap(value)
    sequence = _sequence_form(value)
    if sequence is None or len(sequence) != 3:
        raise HarnessError(f"{label} must carry a three-component dimension vector")
    result: list[sp.Expr] = []
    for component in sequence:
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
        self.comparisons = 0

    def walk_container(self, value: object) -> None:
        value = _unwrap(value)
        if isinstance(value, CasRule):
            sides = (value.left, value.right)
            dimensions = tuple(self.dimension(self._expression(side)) for side in sides)
            if (
                all(dimension is not None for dimension in dimensions)
                and sides[0] != 0
                and sides[1] != 0
            ):
                self.comparisons += 1
                if not _dimension_equal(dimensions[0], dimensions[1]):  # type: ignore[arg-type]
                    self._record(f"{sides[0]} -> {sides[1]}", sides, dimensions)
            return
        if _is_boolean_value(value) and not isinstance(
            value, sp.core.relational.Relational
        ):
            return
        if isinstance(value, Mapping):
            return
        if isinstance(value, sp.MatrixBase):
            for item in value:
                self.walk_container(item)
            return
        if isinstance(value, (list, tuple)):
            for item in value:
                self.walk_container(item)
            return
        if isinstance(value, sp.Contains):
            return
        if isinstance(value, sp.core.relational.Relational):
            left = self.dimension(value.lhs)
            right = self.dimension(value.rhs)
            if (
                left is not None
                and right is not None
                and value.lhs != 0
                and value.rhs != 0
            ):
                self.comparisons += 1
                if not _dimension_equal(left, right):
                    self._record(value, (value.lhs, value.rhs), (left, right))
            return
        if isinstance(value, (int, sp.Basic)):
            self.dimension(sp.Integer(value) if isinstance(value, int) else value)
            return
        raise DimensionError(f"unsupported normalized value {type(value).__name__}")

    @staticmethod
    def _expression(value: object) -> sp.Expr:
        if isinstance(value, bool) or not isinstance(value, (int, sp.Basic)):
            raise DimensionError("Rule sides must be scalar expressions")
        return sp.Integer(value) if isinstance(value, int) else value

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
            if len(known) > 1:
                self.comparisons += len(known) - 1
                if any(
                    not _dimension_equal(known[0][1], value[1])
                    for value in known[1:]
                ):
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
            if exponent_dimension is not None:
                self.comparisons += 1
                if not _dimension_equal(exponent_dimension, ZERO_DIMENSION):
                    self._record(expression, (exponent,), (exponent_dimension,))
            if exponent.is_Rational is not True:
                raise DimensionError(f"non-rational power exponent in {expression}")
            if base_dimension is None:
                return None
            return _dimension_scale(base_dimension, exponent)
        if expression.func == sp.exp:
            argument_dimension = self.dimension(expression.args[0])
            if argument_dimension is not None:
                self.comparisons += 1
                if not _dimension_equal(argument_dimension, ZERO_DIMENSION):
                    self._record(
                        expression, (expression.args[0],), (argument_dimension,)
                    )
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
    config: Mapping[str, Any], values: Mapping[str, ParsedValue | Unparsed]
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
        # A coefficient may only exist in one control package.  Its absence
        # from the outputs supplied to this run is therefore not an error.
        if tag not in values:
            continue
        value = values[tag]
        if isinstance(value, Unparsed):
            raise HarnessError(f"derived dimension tag is UNPARSED: {tag}: {value.error}")
        raw_value = _unwrap(value)
        if isinstance(raw_value, (list, tuple)) and not raw_value:
            continue
        table[symbol] = _as_dimension(value, tag)
    return table


def check_dimensions(
    values: Mapping[str, ParsedValue | Unparsed], config: Mapping[str, Any]
) -> DimensionReport:
    """Check every parsed emission using dimensions generated by this run."""
    dimensions = _dimension_table(config, values)
    non_homogeneous: list[NonHomogeneous] = []
    unknown_symbols: list[UnknownSymbol] = []
    errors: list[str] = []
    unparsed: list[str] = []
    non_dimensional: list[NonDimensional] = []
    checked = 0
    homogeneous = 0
    vacuous: list[str] = []
    for tag, value in values.items():
        if isinstance(value, Unparsed):
            unparsed.append(tag)
            continue
        if value.kind not in DIMENSIONFUL_KINDS:
            non_dimensional.append(NonDimensional(tag, value.kind))
            continue
        walker = _DimensionWalker(dimensions, tag)
        try:
            walker.walk_container(value)
        except DimensionError as exc:
            errors.append(f"{tag}: {exc}")
            continue
        non_homogeneous.extend(walker.non_homogeneous)
        unknown_symbols.extend(UnknownSymbol(tag, name) for name in sorted(walker.unknown))
        if walker.comparisons:
            checked += 1
            if not walker.non_homogeneous and not walker.unknown:
                homogeneous += 1
        elif not walker.non_homogeneous and not walker.unknown:
            vacuous.append(tag)
    return DimensionReport(
        total_tags=len(values),
        checked_tags=checked,
        homogeneous_tags=homogeneous,
        vacuous_tags=tuple(vacuous),
        non_homogeneous=tuple(non_homogeneous),
        unknown_symbols=tuple(unknown_symbols),
        errors=tuple(errors),
        unparsed=tuple(unparsed),
        non_dimensional=tuple(non_dimensional),
    )


def check_cross_engine(
    rows: Sequence[Mapping[str, str]],
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
) -> CrossEngineReport:
    reserved = {"quantity", "mode"}

    def select_shape(value: ParsedValue, selector: str | None) -> object:
        selected = value.value
        if selector is None:
            return selected
        if selector == "scalar":
            if _sequence_form(selected) is not None or isinstance(selected, Mapping):
                raise ValueError("expected scalar")
            return selected
        if selector == "list":
            if not isinstance(selected, (list, tuple)) or len(selected) != 1:
                raise ValueError("expected singleton list")
            return selected[0]
        if selector == "matrix":
            if not isinstance(selected, sp.MatrixBase):
                raise ValueError("expected matrix")
            return selected
        if selector in {"list_of_pairs_second", "last_pair_second"}:
            if not isinstance(selected, (list, tuple)) or not selected:
                raise ValueError("expected nonempty list of pairs")
            pairs = [_sequence_form(item) for item in selected]
            if any(pair is None or len(pair) != 2 for pair in pairs):
                raise ValueError("expected list of pairs")
            seconds = [pair[1] for pair in pairs if pair is not None]
            if selector == "last_pair_second":
                return seconds[-1]
            return seconds[0] if len(seconds) == 1 else seconds
        if selector == "sequence_third":
            sequence = _sequence_form(selected)
            if sequence is None or len(sequence) != 3:
                raise ValueError("expected three-item sequence")
            return sequence[2]
        raise HarnessError(f"unknown cross-engine selector {selector!r}")

    results: list[CrossEngineRow] = []
    for row in rows:
        quantity = row.get("quantity")
        if not isinstance(quantity, str):
            raise HarnessError("each cross_engine row needs a quantity name")
        tag_map = tuple(
            (engine, tag)
            for engine, tag in row.items()
            if engine not in reserved and not engine.endswith("_select")
        )
        if len(tag_map) < 2 or any(not isinstance(tag, str) for _, tag in tag_map):
            raise HarnessError(f"cross_engine {quantity} needs at least two engine tag names")
        missing = tuple(
            f"{engine}:{tag}"
            for engine, tag in tag_map
            if engine not in outputs or tag not in outputs[engine]
        )
        operands: tuple[tuple[str, str], ...] = ()
        detail = ""
        if missing:
            status = "MISSING"
            detail = ", ".join(missing)
        else:
            compared = tuple(outputs[engine][tag] for engine, tag in tag_map)
            if any(isinstance(value, Unparsed) for value in compared):
                status = "UNPARSED"
                operands = tuple(
                    (
                        engine,
                        value.raw if isinstance(value, Unparsed) else str(value.value),
                    )
                    for (engine, _tag), value in zip(tag_map, compared)
                )
            else:
                parsed_compared = tuple(
                    value if isinstance(value, ParsedValue) else _parsed(value)
                    for value in compared
                )
                try:
                    selected = tuple(
                        select_shape(value, row.get(f"{engine}_select"))
                        for (engine, _tag), value in zip(tag_map, parsed_compared)
                    )
                except ValueError as exc:
                    status = "SHAPE_MISMATCH"
                    detail = str(exc)
                    selected = tuple(
                        value.value for value in parsed_compared
                    )
                operands = tuple(
                    (engine, str(value))
                    for (engine, _tag), value in zip(tag_map, selected)
                )
                if not detail:
                    mode = row.get("mode", "symbolic")
                    if mode not in {"symbolic", "multiset"}:
                        raise HarnessError(
                            f"cross_engine {quantity} has unknown mode {mode!r}"
                        )
                    shape_ok = all(
                        symbolic_shapes_compatible(selected[0], value)
                        for value in selected[1:]
                    )
                    if not shape_ok:
                        status = "SHAPE_MISMATCH"
                    else:
                        comparator = (
                            symbolic_multiset_equal if mode == "multiset" else symbolic_equal
                        )
                        status = (
                            "AGREE"
                            if all(comparator(selected[0], value) for value in selected[1:])
                            else "DISAGREE"
                        )
        results.append(CrossEngineRow(quantity, status, tag_map, operands, detail))
    return CrossEngineReport(tuple(results))


def check_registry_residuals(
    rows: Sequence[Mapping[str, Any]],
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
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
        raw_values = {qid: _unwrap(value) for qid, value in values.items()}
        if any(
            isinstance(value, (bool, list, tuple, Mapping, sp.MatrixBase, CasRule))
            or _is_boolean_value(value)
            for value in raw_values.values()
        ):
            raise HarnessError(f"registry {relation_id} substitutions must be scalar expressions")
        substitutions = {registry.symbols[qid]: value for qid, value in raw_values.items()}
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


def load_output(
    path: Path | str, syntax: str | None = None
) -> tuple[dict[str, str], dict[str, ParsedValue | Unparsed]]:
    output_path = Path(path)
    try:
        text = output_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise HarnessError(f"cannot read output {output_path}: {exc}") from exc
    try:
        raw = parse_tagged_output(text)
    except ValueError as exc:
        raise HarnessError(f"cannot parse tagged output {output_path}: {exc}") from exc
    return raw, normalize_tags(raw, syntax=syntax)


def run_checks(
    config: Mapping[str, Any],
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
    *,
    registry_directory: Path | str | None = None,
) -> HarnessReport:
    default_engine = config.get("default_engine", "default")
    if not isinstance(default_engine, str) or default_engine not in outputs:
        raise HarnessError(f"default_engine {default_engine!r} has no loaded output")
    default_values = outputs[default_engine]
    cross_rows = config.get("cross_engine", [])
    registry_rows = config.get("registry_residual", [])
    parity_exclude = config.get("parity_exclude", [])
    main_package = config.get("main_package")
    control_packages = config.get("control_packages")
    if not isinstance(cross_rows, list) or not isinstance(registry_rows, list):
        raise HarnessError("cross_engine and registry_residual must be lists")
    if not isinstance(parity_exclude, list):
        raise HarnessError("parity_exclude must be a list")
    if main_package is None and control_packages is not None:
        raise HarnessError("control_packages requires main_package")
    if main_package is not None and (
        not isinstance(main_package, str) or not main_package
    ):
        raise HarnessError("main_package must be a non-empty string")
    if control_packages is not None and (
        not isinstance(control_packages, list)
        or any(not isinstance(package, str) or not package for package in control_packages)
    ):
        raise HarnessError("control_packages must be a list of non-empty strings")
    declared_controls = control_packages if control_packages is not None else None
    registry = load_registry(registry_directory or Path(__file__).resolve().parent)
    all_unparsed = tuple(
        (engine, tag, value)
        for engine, values in outputs.items()
        for tag, value in values.items()
        if isinstance(value, Unparsed)
    )
    return HarnessReport(
        controls=check_control_response(
            default_values, main_package, declared_controls
        ),
        parity=check_tag_parity(
            outputs, parity_exclude, main_package, declared_controls
        ),
        cross_engine_tags=check_cross_engine_tag_sets(outputs, parity_exclude),
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
    parity_gaps = [row for row in report.parity.rows if row.missing or row.extra]
    compared_controls = len(control.responsive) + len(control.invariant)
    main_control_tags = compared_controls + len(control.uncovered)
    uncovered_fraction = (
        len(control.uncovered) / main_control_tags if main_control_tags else 0.0
    )
    no_ablation = ",".join(
        f"D{value}" for value in control.no_ablation_dimensions
    ) or "none"
    non_dimensional_counts: dict[str, int] = {}
    for row in dimension.non_dimensional:
        non_dimensional_counts[row.kind.value] = non_dimensional_counts.get(row.kind.value, 0) + 1
    non_dimensional_detail = ",".join(
        f"{name}={count}" for name, count in sorted(non_dimensional_counts.items())
    ) or "none"
    lines = [
        (
            f"CONTROL_RESPONSE: compared={compared_controls} "
            f"responsive={len(control.responsive)} invariant={len(control.invariant)} "
            f"unparsed={len(control.unparsed)} unpaired={len(control.unpaired)}"
        ),
        (
            f"CONTROL_COVERAGE: main={main_control_tags} compared={compared_controls} "
            f"uncovered={len(control.uncovered)} "
            f"uncovered_fraction={uncovered_fraction:.6f} "
            f"no_ablation_D=[{no_ablation}]"
        ),
        (
            "CONTROL_PACKAGES: "
            f"matched={len(control.controls)} "
            f"missing_declared=[{','.join(control.missing_declared_packages) or 'none'}]"
        ),
        (
            f"DIMENSIONS: total={dimension.total_tags} compared={dimension.checked_tags} "
            f"homogeneous={dimension.homogeneous_tags} "
            f"vacuous={len(dimension.vacuous_tags)} "
            f"non_homogeneous={len(dimension.non_homogeneous)} "
            f"unknown_symbol={len(dimension.unknown_symbols)} "
            f"unparsed={len(dimension.unparsed)}"
        ),
        (
            f"NON_DIMENSIONAL: total={len(dimension.non_dimensional)} "
            f"kinds=[{non_dimensional_detail}]"
        ),
        f"CROSS_ENGINE: {_status_counts(report.cross_engine.rows)}",
        f"REGISTRY_RESIDUAL: {_status_counts(report.registry.rows)}",
        (
            f"MAIN_CONTROL_TAG_PARITY: comparisons={len(report.parity.rows)} "
            f"gaps={len(parity_gaps)}"
        ),
    ]
    cross_tag_gap_count = sum(
        len(row.left_only) + len(row.right_only)
        for row in report.cross_engine_tags.rows
    )
    lines.append(
        f"CROSS_ENGINE_TAG_PARITY: comparisons={len(report.cross_engine_tags.rows)} "
        f"gaps={cross_tag_gap_count}"
    )
    for row in report.cross_engine_tags.rows:
        lines.append(
            f"  {row.left_engine}_only ({len(row.left_only)}): "
            + (", ".join(row.left_only) or "-")
        )
        lines.append(
            f"  {row.right_engine}_only ({len(row.right_only)}): "
            + (", ".join(row.right_only) or "-")
        )
    if report.parity.exclusions:
        excluded_total = sum(len(row.tags) for row in report.parity.exclusions)
        lines.append(f"PARITY_EXCLUDED ({excluded_total}):")
        for row in report.parity.exclusions:
            by_pattern = ", ".join(
                f"{pattern!r}:{count}" for pattern, count in row.by_pattern
            )
            lines.append(
                f"  {row.engine}: excluded={len(row.tags)} "
                f"by_pattern={{{by_pattern}}}"
            )
    if report.parity.exclusion_warnings:
        lines.append(
            f"PARITY_EXCLUSION_WARNING ({len(report.parity.exclusion_warnings)}):"
        )
        lines.extend(
            f"  {warning}" for warning in report.parity.exclusion_warnings
        )
    disagreements = [
        row for row in report.cross_engine.rows if row.status == "DISAGREE"
    ]
    registry_disagreements = [
        row for row in report.registry.rows if row.status == "NONZERO"
    ]
    lines.append(f"DISAGREE ({len(disagreements) + len(registry_disagreements)}):")
    for row in disagreements:
        operands = "; ".join(
            f"{engine}[{dict(row.tags).get(engine)}]={operand}"
            for engine, operand in row.operands
        )
        lines.append(f"  cross_engine:{row.quantity}: {operands}")
    lines.extend(
        f"  registry:{row.relation_id}: residual={row.residual}"
        for row in registry_disagreements
    )
    invariant_names = ", ".join(row.tag for row in control.invariant)
    lines.append(f"INVARIANT ({len(control.invariant)}): {invariant_names}")
    lines.append(
        f"CONTROL_UNCOVERED ({len(control.uncovered)}): "
        + ", ".join(control.uncovered)
    )
    lines.append(
        f"DIMENSION_VACUOUS ({len(dimension.vacuous_tags)}): "
        + ", ".join(dimension.vacuous_tags)
    )
    lines.append(f"NON_HOMOGENEOUS ({len(dimension.non_homogeneous)}):")
    for issue in dimension.non_homogeneous:
        rendered = "; ".join(
            f"{term.expression} => {_format_dimension(term.dimension)}"
            for term in issue.summands
        )
        lines.append(f"  {issue.tag}: {rendered}")
    lines.append(f"PARITY ({len(parity_gaps)}):")
    lines.append(
        "  why: present-and-identical is INVARIANT; absent is indistinguishable from never computed."
    )
    parity_groups: dict[
        tuple[str, int | None, tuple[str, ...], tuple[str, ...]], list[str]
    ] = {}
    for row in parity_gaps:
        parity_groups.setdefault(
            (row.engine, row.dimension, row.missing, row.extra), []
        ).append(row.package)
    for (
        engine,
        dimension_value,
        missing_values,
        extra_values,
    ), packages in parity_groups.items():
        missing = ",".join(missing_values) or "-"
        extra = ",".join(extra_values) or "-"
        dimension_label = (
            f"D{dimension_value}" if dimension_value is not None else "unswept"
        )
        lines.append(
            f"  {engine}:{dimension_label}:packages=[{','.join(packages)}]: "
            f"missing=[{missing}] extra=[{extra}]"
        )
    lines.append(f"UNPARSED ({len(report.unparsed)}):")
    lines.extend(
        f"  {engine}:{tag}: {value.error}; raw={value.raw!r}"
        for engine, tag, value in report.unparsed
    )
    lines.append(f"UNKNOWN_SYMBOL ({len(dimension.unknown_symbols)}):")
    lines.extend(
        f"  {issue.tag}: {issue.symbol}" for issue in dimension.unknown_symbols
    )
    lines.append(f"DIMENSION_ERROR ({len(dimension.errors)}):")
    lines.extend(f"  {error}" for error in dimension.errors)
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
