#!/usr/bin/env python3
"""Declaration-oriented checks for tagged symbolic-engine audit output.

The harness contains no expected physics result.  It distinguishes comparison
verdicts from failures to compare, and derives every coverage denominator from
the configuration before inspecting engine output.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, MutableMapping, Sequence

import sympy as sp
import yaml

try:
    from .registry_read import Registry, load_registry
except ImportError:  # Direct script invocation.
    from registry_read import Registry, load_registry


TAG_LINE = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)[ \t]*(?::|->)[ \t]*(.*)$")
TAG_NAME = re.compile(r"^(?:WL|PY)_S[1-9][0-9]*_[A-Z][A-Z0-9_]*$")
TOKEN = re.compile(
    r"(?P<WS>\s+)|(?P<ASSOC_OPEN><\|)|(?P<ASSOC_CLOSE>\|>)|"
    r"(?P<AND>&&)|(?P<OR>\|\|)|(?P<ARROW>->)|(?P<GE>>=)|(?P<LE><=)|"
    r"(?P<EQUAL>==)|(?P<NE>!=)|(?P<GT>>)|(?P<LT><)|(?P<NOT>!)|"
    r"(?P<STRING>\"(?:\\.|[^\"\\])*\")|"
    r"(?P<REAL>(?:[0-9]+\.[0-9]*|\.[0-9]+)(?:`[0-9.]*)?(?:\*\^[+-]?[0-9]+)?)|"
    r"(?P<INT>[0-9]+)|(?P<IDENT>[A-Za-z][A-Za-z0-9_$]*)|"
    r"(?P<PUNCT>[+\-*/^{}()[\],])"
)
ZERO_DIMENSION = (sp.Integer(0), sp.Integer(0), sp.Integer(0))
VERDICTS = frozenset({"AGREE", "DISAGREE"})
REGISTRY_VERDICTS = frozenset({"ZERO", "NONZERO"})
OPERATIONAL_CROSS_STATUSES = frozenset(
    {
        "MISSING",
        "UNPARSED",
        "SHAPE_MISMATCH",
        "EMPTY",
        "CARDINALITY_INVALID",
        "DECLARATION_ERROR",
        "NAMING_MISMATCH",
    }
)


class HarnessError(RuntimeError):
    """A configuration or runtime failure that prevents a requested check."""


class CasParseError(ValueError):
    """Strict Mathematica InputForm parsing failed."""


class DimensionError(ValueError):
    """A recognized dimensional path could not be assessed."""


class UnwalkedDimension(DimensionError):
    """A recognized container/path has no implemented walker."""


class TaggedRecords(dict[str, str]):
    """Parsed records plus non-record diagnostics retained as evidence."""

    def __init__(
        self,
        *args: object,
        ignored_lines: Sequence[str] = (),
        duplicate_tags: Sequence[str] = (),
        **kwargs: object,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.ignored_lines = tuple(ignored_lines)
        self.duplicate_tags = tuple(duplicate_tags)


class NormalizedRecords(dict[str, "ParsedValue | Unparsed"]):
    def __init__(
        self,
        *args: object,
        ignored_lines: Sequence[str] = (),
        duplicate_tags: Sequence[str] = (),
        **kwargs: object,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.ignored_lines = tuple(ignored_lines)
        self.duplicate_tags = tuple(duplicate_tags)


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
    REAL = "real"
    STRING = "string"
    CALL = "call"


@dataclass(frozen=True)
class ParsedValue:
    kind: ValueKind
    value: object


@dataclass(frozen=True)
class CasRule:
    left: object
    right: object


@dataclass(frozen=True)
class CasCall:
    head: str
    arguments: tuple[object, ...]


@dataclass(frozen=True)
class CasRelation:
    operator: str
    left: object
    right: object


@dataclass(frozen=True)
class CasInequality:
    operands: tuple[object, ...]
    operators: tuple[str, ...]


@dataclass(frozen=True)
class CasBoolean:
    operator: str
    operands: tuple[object, ...]


class OpaqueDerivative(sp.Symbol):
    """A WL derivative whose complete emitted identity remains attached."""

    def __new__(
        cls,
        rendered: str,
        orders: tuple[int, ...],
        function_name: str,
        variables: tuple[str, ...],
    ) -> "OpaqueDerivative":
        obj = sp.Symbol.__new__(cls, rendered)
        obj.orders = orders
        obj.function_name = function_name
        obj.variables = variables
        return obj


class CanonicalDerivative(sp.Symbol):
    """Engine-neutral derivative atom keyed only by mathematical identity."""

    def __new__(
        cls,
        function_name: str,
        function_arguments: tuple[str, ...],
        differentiated: tuple[tuple[str, int], ...],
    ) -> "CanonicalDerivative":
        arguments = ",".join(function_arguments)
        pairs = ",".join(f"{variable}^{order}" for variable, order in differentiated)
        obj = sp.Symbol.__new__(
            cls, f"DerivativeIdentity[{function_name}({arguments});{pairs}]"
        )
        obj.function_name = function_name
        obj.function_arguments = function_arguments
        obj.differentiated = differentiated
        return obj


def _is_boolean_value(value: object) -> bool:
    """Population whose truthiness exposure is observed, not repaired."""
    return (
        isinstance(value, bool)
        or getattr(value, "is_Boolean", False) is True
        or isinstance(value, sp.assumptions.assume.AppliedPredicate)
        or isinstance(value, (sp.core.relational.Relational, sp.Contains))
        or isinstance(value, (CasRelation, CasInequality, CasBoolean))
    )


def parse_tagged_output(stdout: str) -> TaggedRecords:
    """Parse one-record-per-line output; retain diagnostics and duplicates."""
    tags: dict[str, str] = {}
    ignored: list[str] = []
    duplicates: list[str] = []
    for line in stdout.splitlines():
        match = TAG_LINE.match(line)
        if match is None or TAG_NAME.fullmatch(match.group(1)) is None:
            if line.strip():
                ignored.append(line)
            continue
        tag, payload = match.groups()
        if tag in tags:
            duplicates.append(tag)
            continue
        tags[tag] = payload
    if not tags:
        raise ValueError("engine emitted no valid tagged records")
    return TaggedRecords(tags, ignored_lines=ignored, duplicate_tags=duplicates)


@dataclass(frozen=True)
class _Token:
    kind: str
    text: str
    offset: int


def _tokenize(text: str) -> tuple[_Token, ...]:
    tokens: list[_Token] = []
    offset = 0
    while offset < len(text):
        match = TOKEN.match(text, offset)
        if match is None:
            raise CasParseError(
                f"unsupported token at offset {offset}: {text[offset:offset + 24]!r}"
            )
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
            right = self.parse_and()
            self._boolean(result, "|| left")
            self._boolean(right, "|| right")
            result = self._combine_boolean("Or", result, right)
        return result

    def parse_and(self) -> object:
        result = self.parse_not()
        while self.accept("&&"):
            right = self.parse_not()
            self._boolean(result, "&& left")
            self._boolean(right, "&& right")
            result = self._combine_boolean("And", result, right)
        return result

    def parse_not(self) -> object:
        if self.accept("!"):
            value = self._boolean(self.parse_not(), "! operand")
            try:
                return sp.Not(value, evaluate=False)
            except (TypeError, sp.SympifyError):
                return CasBoolean("Not", (value,))
        return self.parse_relation()

    def parse_relation(self) -> object:
        left = self.parse_add()
        if self.accept("->"):
            return CasRule(left, self.parse_or())
        operators: list[str] = []
        operands: list[object] = [left]
        while self.token.text in {"==", "!=", ">", ">=", "<", "<="}:
            operators.append(self.token.text)
            self.position += 1
            operands.append(self.parse_add())
        if not operators:
            return left
        if len(operators) > 1:
            return CasInequality(tuple(operands), tuple(operators))
        operator = operators[0]
        lhs, rhs = operands
        if self._is_scalar(lhs) and self._is_scalar(rhs):
            relation = {
                "==": sp.Eq,
                "!=": sp.Ne,
                ">": sp.Gt,
                ">=": sp.Ge,
                "<": sp.Lt,
                "<=": sp.Le,
            }[operator]
            return relation(self._scalar(lhs, "relation lhs"), self._scalar(rhs, "relation rhs"), evaluate=False)
        return CasRelation(operator, lhs, rhs)

    def parse_add(self) -> object:
        result = self.parse_mul()
        while self.token.text in {"+", "-"}:
            operator = self.token.text
            self.position += 1
            result = self._binary(operator, result, self.parse_mul())
        return result

    def parse_mul(self) -> object:
        result = self.parse_unary()
        while self.token.text in {"*", "/"}:
            operator = self.token.text
            self.position += 1
            result = self._binary(operator, result, self.parse_unary())
        return result

    def parse_unary(self) -> object:
        if self.accept("+"):
            return self.parse_unary()
        if self.accept("-"):
            return -self._scalar(self.parse_unary(), "unary minus")
        return self.parse_power()

    def parse_power(self) -> object:
        base = self.parse_primary()
        if self.accept("^"):
            return self._scalar(base, "power base") ** self._scalar(
                self.parse_unary(), "power exponent"
            )
        return base

    def parse_primary(self) -> object:
        token = self.token
        if token.kind == "INT":
            self.position += 1
            return int(token.text)
        if token.kind == "REAL":
            self.position += 1
            rendered = token.text.split("`", 1)[0].replace("*^", "e")
            return sp.Float(rendered)
        if token.kind == "STRING":
            self.position += 1
            import ast

            return ast.literal_eval(token.text)
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
                value = self._function_call(value, self.parse_arguments())
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
                        raise CasParseError("association entries require key -> value")
                    try:
                        if rule.left in values:
                            raise CasParseError(f"duplicate association key {rule.left!r}")
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
            not orders_raw
            or any(isinstance(value, bool) or not isinstance(value, int) for value in orders_raw)
            or len(function_raw) != 1
            or not isinstance(function_raw[0], sp.Symbol)
            or len(variables_raw) != len(orders_raw)
            or any(not isinstance(value, sp.Symbol) for value in variables_raw)
        ):
            raise CasParseError(
                "Derivative requires equal nonzero order/argument arity and one symbolic function"
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
    def _is_scalar(value: object) -> bool:
        return not isinstance(value, bool) and isinstance(value, (int, sp.Basic))

    @classmethod
    def _scalar(cls, value: object, context: str) -> sp.Expr:
        if not cls._is_scalar(value):
            raise CasParseError(f"{context} is not a scalar expression")
        return sp.Integer(value) if isinstance(value, int) else value

    @staticmethod
    def _boolean(value: object, context: str) -> object:
        if _is_boolean_value(value):
            return value
        raise CasParseError(f"{context} is not boolean")

    @staticmethod
    def _combine_boolean(operator: str, left: object, right: object) -> object:
        try:
            constructor = sp.And if operator == "And" else sp.Or
            return constructor(left, right, evaluate=False)
        except (TypeError, sp.SympifyError):
            operands: list[object] = []
            for value in (left, right):
                if isinstance(value, CasBoolean) and value.operator == operator:
                    operands.extend(value.operands)
                else:
                    operands.append(value)
            return CasBoolean(operator, tuple(operands))

    def _binary(self, operator: str, left: object, right: object) -> sp.Expr:
        lhs = self._scalar(left, f"{operator} left")
        rhs = self._scalar(right, f"{operator} right")
        return {"+": lhs + rhs, "-": lhs - rhs, "*": lhs * rhs, "/": lhs / rhs}[operator]

    def _function_call(self, head: object, arguments: tuple[object, ...]) -> object:
        if not isinstance(head, sp.Symbol):
            raise CasParseError(f"unsupported function head {head!r}")
        name = str(head)
        if name == "Element":
            if len(arguments) != 2 or not isinstance(arguments[1], sp.Symbol):
                raise CasParseError("Element requires a member and named set")
            domain = {"Reals": sp.S.Reals, "Integers": sp.S.Integers}.get(str(arguments[1]))
            if domain is None:
                return CasCall(name, arguments)
            return sp.Contains(arguments[0], domain, evaluate=False)
        if name == "Inequality":
            if len(arguments) < 3 or len(arguments) % 2 == 0:
                raise CasParseError("Inequality requires alternating operands/operators")
            operators = arguments[1::2]
            if any(not isinstance(item, sp.Symbol) for item in operators):
                raise CasParseError("Inequality operators must be named symbols")
            return CasInequality(arguments[0::2], tuple(map(str, operators)))
        if name == "Piecewise":
            if len(arguments) not in {1, 2}:
                raise CasParseError("Piecewise requires branches and optional default")
            return CasCall(name, arguments)
        if name == "ConditionalExpression":
            if len(arguments) != 2:
                raise CasParseError("ConditionalExpression requires expression and condition")
            return CasCall(name, arguments)
        if name == "Sqrt" and len(arguments) == 1:
            return sp.sqrt(self._scalar(arguments[0], "Sqrt argument"))
        if name == "Rational" and len(arguments) == 2:
            left = self._scalar(arguments[0], "Rational numerator")
            right = self._scalar(arguments[1], "Rational denominator")
            if left.is_Integer is not True or right.is_Integer is not True or right == 0:
                raise CasParseError("Rational expects integers and nonzero denominator")
            return sp.Rational(int(left), int(right))
        if name == "Exp" and len(arguments) == 1:
            return sp.exp(self._scalar(arguments[0], "Exp argument"))
        if all(self._is_scalar(argument) for argument in arguments):
            return sp.Function(name)(*(self._scalar(arg, f"{name} argument") for arg in arguments))
        return CasCall(name, arguments)


def value_kind(value: object) -> ValueKind:
    if isinstance(value, ParsedValue):
        return value.kind
    if isinstance(value, str):
        return ValueKind.STRING
    if isinstance(value, CasCall):
        return ValueKind.CALL
    if isinstance(value, (CasRule, CasRelation, CasInequality, CasBoolean, sp.core.relational.Relational, sp.Contains)):
        return ValueKind.RELATION
    if _is_boolean_value(value):
        return ValueKind.BOOLEAN
    if isinstance(value, Mapping):
        return ValueKind.MAPPING
    if isinstance(value, sp.MatrixBase):
        return ValueKind.MATRIX
    if isinstance(value, (list, tuple, set, frozenset)):
        return ValueKind.LIST
    if isinstance(value, (int, sp.Integer)) and not isinstance(value, bool):
        return ValueKind.INTEGER
    if isinstance(value, sp.Float):
        return ValueKind.REAL
    if isinstance(value, sp.Basic):
        return ValueKind.EXPRESSION
    raise CasParseError(f"unsupported normalized value {type(value).__name__}")


def _parsed(value: object) -> ParsedValue:
    return ParsedValue(value_kind(value), value)


def normalize_mathematica(raw: str) -> ParsedValue | Unparsed:
    try:
        return _parsed(_MathematicaParser(raw).parse())
    except (CasParseError, TypeError, ValueError, SyntaxError) as exc:
        return Unparsed(raw, str(exc))


def normalize_sympy(raw: str) -> ParsedValue | Unparsed:
    try:
        value = sp.sympify(
            raw,
            locals={"Matrix": sp.Matrix, "Eq": sp.Eq, "Derivative": sp.Derivative, "Q": sp.Q},
            evaluate=False,
        )
        return _parsed(value)
    except Exception as exc:  # SymPy exposes several parser exception types.
        return Unparsed(raw, str(exc))


def _normalize_auto(raw: str) -> ParsedValue | Unparsed:
    result = normalize_mathematica(raw)
    return normalize_sympy(raw) if isinstance(result, Unparsed) else result


def normalize_tags(
    tags: Mapping[str, str], syntax: str | None = None
) -> NormalizedRecords:
    if syntax not in {None, "wl", "py"}:
        raise ValueError(f"unknown CAS syntax {syntax!r}")
    normalized: dict[str, ParsedValue | Unparsed] = {}
    for tag, raw in tags.items():
        selected = syntax or ("py" if tag.startswith("PY_") else "wl" if tag.startswith("WL_") else None)
        normalizer = normalize_sympy if selected == "py" else normalize_mathematica if selected == "wl" else _normalize_auto
        normalized[tag] = normalizer(raw)
    return NormalizedRecords(
        normalized,
        ignored_lines=getattr(tags, "ignored_lines", ()),
        duplicate_tags=getattr(tags, "duplicate_tags", ()),
    )


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
    left, right = _unwrap(left), _unwrap(right)
    left_sequence, right_sequence = _sequence_form(left), _sequence_form(right)
    if left_sequence is not None or right_sequence is not None:
        return (
            left_sequence is not None
            and right_sequence is not None
            and len(left_sequence) == len(right_sequence)
            and all(symbolic_shapes_compatible(a, b) for a, b in zip(left_sequence, right_sequence))
        )
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        return isinstance(left, Mapping) and isinstance(right, Mapping) and len(left) == len(right)
    structured = (CasRule, CasCall, CasRelation, CasInequality, CasBoolean)
    if isinstance(left, structured) or isinstance(right, structured):
        return type(left) is type(right)
    return _is_boolean_value(left) == _is_boolean_value(right)


def symbolic_equal(left: object, right: object) -> bool:
    """Existing equality semantics, intentionally unchanged by this rebuild."""
    if isinstance(left, Unparsed) or isinstance(right, Unparsed):
        raise ValueError("UNPARSED values cannot be compared")
    left, right = _unwrap(left), _unwrap(right)
    left_sequence, right_sequence = _sequence_form(left), _sequence_form(right)
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
        return isinstance(left, CasRule) and isinstance(right, CasRule) and symbolic_equal(left.left, right.left) and symbolic_equal(left.right, right.right)
    if isinstance(left, CasCall) or isinstance(right, CasCall):
        return isinstance(left, CasCall) and isinstance(right, CasCall) and left.head == right.head and symbolic_equal(left.arguments, right.arguments)
    if isinstance(left, CasRelation) or isinstance(right, CasRelation):
        return isinstance(left, CasRelation) and isinstance(right, CasRelation) and left.operator == right.operator and symbolic_equal(left.left, right.left) and symbolic_equal(left.right, right.right)
    if isinstance(left, CasInequality) or isinstance(right, CasInequality):
        return isinstance(left, CasInequality) and isinstance(right, CasInequality) and left.operators == right.operators and symbolic_equal(left.operands, right.operands)
    if isinstance(left, CasBoolean) or isinstance(right, CasBoolean):
        return isinstance(left, CasBoolean) and isinstance(right, CasBoolean) and left.operator == right.operator and symbolic_equal(left.operands, right.operands)
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
    return left == right


def symbolic_multiset_equal(left: object, right: object) -> bool:
    left_sequence, right_sequence = _sequence_form(left), _sequence_form(right)
    if left_sequence is None or right_sequence is None or len(left_sequence) != len(right_sequence):
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


# ---------------------------------------------------------------------------
# Declarations, naming, and cross-engine comparison
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Coverage:
    numerator: int
    denominator: int
    gaps: tuple[str, ...]
    formula: str


@dataclass(frozen=True)
class CrossEngineRow:
    quantity: str
    identity: tuple[tuple[str, str], ...]
    status: str
    family: str | None = None
    operands: tuple[tuple[str, str], ...] = ()
    detail: str = ""
    naming_applied: tuple[str, ...] = ()
    identities_applied: tuple[str, ...] = ()
    undeclared_spellings: tuple[str, ...] = ()
    selected_boolean: bool = False
    tree_boolean: bool = False
    legacy_algebraic_status: str | None = None


@dataclass(frozen=True)
class CrossEngineReport:
    rows: tuple[CrossEngineRow, ...]
    coverage: Coverage
    duplicates: tuple[str, ...]
    declaration_failure: str | None
    selected_boolean_rows: tuple[str, ...]
    tree_boolean_rows: tuple[str, ...]
    naming_before_agree: int
    naming_after_agree: int
    naming_changed_rows: tuple[str, ...]
    naming_exceptions: tuple[str, ...]


@dataclass(frozen=True)
class ActionFamilyReport:
    family: str
    coverage: Coverage
    statuses: tuple[tuple[str, str], ...]


def _snake_to_lower_camel(name: str) -> str:
    pieces = name.split("_")
    return pieces[0] + "".join(piece[:1].upper() + piece[1:] for piece in pieces[1:])


@dataclass(frozen=True)
class _NamingRule:
    engine_maps: Mapping[str, Mapping[str, str]]
    exceptions: tuple[str, ...]


@dataclass(frozen=True)
class _SymbolIdentities:
    engine_maps: Mapping[str, Mapping[sp.Symbol, sp.Expr]]


def _build_naming_rule(
    config: Mapping[str, Any], registry: Registry | None
) -> _NamingRule:
    declaration = config.get("symbol_naming")
    if not isinstance(declaration, Mapping):
        return _NamingRule({}, ())
    if declaration.get("rule") != "registry_snake_case_to_lower_camel":
        raise HarnessError("symbol_naming.rule must be registry_snake_case_to_lower_camel")
    styles = declaration.get("engine_styles")
    if not isinstance(styles, Mapping) or any(
        not isinstance(engine, str) or style not in {"canonical", "lower_camel"}
        for engine, style in styles.items()
    ):
        raise HarnessError("symbol_naming.engine_styles has an invalid declaration")
    canonical_names = (
        {str(symbol) for symbol in registry.symbols.values()} if registry is not None else set()
    )
    exceptions_raw = declaration.get("exceptions", [])
    if not isinstance(exceptions_raw, list):
        raise HarnessError("symbol_naming.exceptions must be a list")
    exceptions: list[str] = []
    exception_spellings: dict[str, dict[str, str]] = {}
    for item in exceptions_raw:
        if (
            not isinstance(item, Mapping)
            or not isinstance(item.get("canonical"), str)
            or not isinstance(item.get("reason"), str)
            or not item.get("reason")
            or not isinstance(item.get("spellings"), Mapping)
        ):
            raise HarnessError("every naming exception needs canonical, spellings, and reason")
        canonical = item["canonical"]
        spellings = item["spellings"]
        if any(not isinstance(engine, str) or not isinstance(name, str) for engine, name in spellings.items()):
            raise HarnessError(f"invalid spellings for naming exception {canonical}")
        exception_spellings[canonical] = dict(spellings)
        exceptions.append(f"{canonical}: {item['reason']}")
    engine_maps: dict[str, dict[str, str]] = {}
    for engine, style in styles.items():
        mapping: dict[str, str] = {}
        for canonical in canonical_names:
            emitted = canonical if style == "canonical" else _snake_to_lower_camel(canonical)
            if emitted != canonical:
                mapping[emitted] = canonical
        for canonical, spellings in exception_spellings.items():
            emitted = spellings.get(engine)
            if emitted and emitted != canonical:
                mapping[emitted] = canonical
        engine_maps[engine] = mapping
    return _NamingRule(engine_maps, tuple(exceptions))


def _build_symbol_identities(config: Mapping[str, Any]) -> _SymbolIdentities:
    declaration = config.get("symbol_identities")
    if declaration is None:
        return _SymbolIdentities({})
    if not isinstance(declaration, list):
        raise HarnessError("symbol_identities must be a list")
    engine_maps: dict[str, dict[sp.Symbol, sp.Expr]] = {}
    for item in declaration:
        if (
            not isinstance(item, Mapping)
            or not isinstance(item.get("engine"), str)
            or not isinstance(item.get("symbol"), str)
            or re.fullmatch(r"[A-Za-z][A-Za-z0-9_$]*", item.get("symbol", "")) is None
            or not isinstance(item.get("expression"), str)
            or not isinstance(item.get("reason"), str)
            or not item.get("reason")
        ):
            raise HarnessError(
                "every symbol identity needs engine, atomic symbol, expression, and reason"
            )
        engine = item["engine"]
        symbol = sp.Symbol(item["symbol"])
        expression = normalize_sympy(item["expression"])
        if (
            isinstance(expression, Unparsed)
            or expression.kind not in {ValueKind.EXPRESSION, ValueKind.INTEGER, ValueKind.REAL}
            or not isinstance(expression.value, sp.Expr)
        ):
            raise HarnessError(
                f"invalid expression for symbol identity {engine}:{item['symbol']}"
            )
        mapping = engine_maps.setdefault(engine, {})
        if symbol in mapping:
            raise HarnessError(f"duplicate symbol identity {engine}:{item['symbol']}")
        mapping[symbol] = expression.value
    return _SymbolIdentities(engine_maps)


def _map_tree(value: object, replacements: Mapping[str, str]) -> tuple[object, tuple[str, ...]]:
    """Apply a declared spelling map recursively and return mappings actually used."""
    applied: set[str] = set()

    def visit(item: object) -> object:
        item = _unwrap(item)
        if isinstance(item, sp.MatrixBase):
            return item.applyfunc(visit)
        if isinstance(item, Mapping):
            return {visit(key): visit(value) for key, value in item.items()}
        if isinstance(item, list):
            return [visit(value) for value in item]
        if isinstance(item, tuple):
            return tuple(visit(value) for value in item)
        if isinstance(item, CasRule):
            return CasRule(visit(item.left), visit(item.right))
        if isinstance(item, CasCall):
            return CasCall(item.head, tuple(visit(value) for value in item.arguments))
        if isinstance(item, CasRelation):
            return CasRelation(item.operator, visit(item.left), visit(item.right))
        if isinstance(item, CasInequality):
            return CasInequality(tuple(visit(value) for value in item.operands), item.operators)
        if isinstance(item, CasBoolean):
            return CasBoolean(item.operator, tuple(visit(value) for value in item.operands))
        if isinstance(item, OpaqueDerivative):
            function = replacements.get(item.function_name, item.function_name)
            variables = tuple(replacements.get(name, name) for name in item.variables)
            if function != item.function_name:
                applied.add(f"{item.function_name}->{function}")
            for old, new in zip(item.variables, variables):
                if old != new:
                    applied.add(f"{old}->{new}")
            return OpaqueDerivative(str(item), item.orders, function, variables)
        if isinstance(item, sp.Derivative):
            changes: dict[sp.Symbol, sp.Symbol] = {}
            for symbol in item.free_symbols:
                old = str(symbol)
                new = replacements.get(old)
                if new is not None and new != old:
                    changes[symbol] = sp.Symbol(new)
                    applied.add(f"{old}->{new}")
            mapped = item.xreplace(changes)
            base = mapped.expr
            if getattr(base, "is_Function", False):
                old = getattr(base.func, "__name__", str(base.func))
                new = replacements.get(old)
                if new is not None and new != old:
                    applied.add(f"{old}->{new}")
                    base = sp.Function(new)(*base.args)
                    mapped = sp.Derivative(base, *mapped.variable_count, evaluate=False)
            return mapped
        if isinstance(item, sp.Basic):
            changes: dict[sp.Symbol, sp.Symbol] = {}
            for symbol in item.free_symbols:
                old = str(symbol)
                new = replacements.get(old)
                if new is not None and new != old:
                    changes[symbol] = sp.Symbol(new)
                    applied.add(f"{old}->{new}")
            return item.xreplace(changes)
        return item

    result = visit(value)
    return result, tuple(sorted(applied))


def _canonicalize_derivatives(value: object) -> object:
    """Replace WL and SymPy derivatives by the same engine-neutral atom."""

    def canonical(item: OpaqueDerivative | sp.Derivative) -> CanonicalDerivative:
        if isinstance(item, OpaqueDerivative):
            function_name = item.function_name
            function_arguments = tuple(sorted(set(item.variables)))
            pairs = tuple(
                (variable, order)
                for variable, order in zip(item.variables, item.orders)
                if order
            )
        else:
            base = item.expr
            function_name = (
                getattr(base.func, "__name__", str(base.func))
                if getattr(base, "is_Function", False)
                else sp.srepr(base)
            )
            function_arguments = (
                tuple(sorted({str(argument) for argument in base.args}))
                if getattr(base, "is_Function", False)
                else ()
            )
            pairs = tuple((str(variable), int(order)) for variable, order in item.variable_count)
        return CanonicalDerivative(
            function_name, function_arguments, tuple(sorted(pairs))
        )

    def visit(item: object) -> object:
        item = _unwrap(item)
        if isinstance(item, sp.MatrixBase):
            return item.applyfunc(visit)
        if isinstance(item, Mapping):
            return {visit(key): visit(mapped) for key, mapped in item.items()}
        if isinstance(item, list):
            return [visit(mapped) for mapped in item]
        if isinstance(item, tuple):
            return tuple(visit(mapped) for mapped in item)
        if isinstance(item, CasRule):
            return CasRule(visit(item.left), visit(item.right))
        if isinstance(item, CasCall):
            return CasCall(item.head, tuple(visit(mapped) for mapped in item.arguments))
        if isinstance(item, CasRelation):
            return CasRelation(item.operator, visit(item.left), visit(item.right))
        if isinstance(item, CasInequality):
            return CasInequality(tuple(visit(mapped) for mapped in item.operands), item.operators)
        if isinstance(item, CasBoolean):
            return CasBoolean(item.operator, tuple(visit(mapped) for mapped in item.operands))
        if isinstance(item, (OpaqueDerivative, sp.Derivative)):
            return canonical(item)
        if isinstance(item, sp.Basic):
            derivatives = {
                node: canonical(node)
                for node in sp.preorder_traversal(item)
                if isinstance(node, (OpaqueDerivative, sp.Derivative))
            }
            return item.xreplace(derivatives)
        return item

    return visit(value)


def _apply_symbol_identities(
    value: object, replacements: Mapping[sp.Symbol, sp.Expr]
) -> tuple[object, tuple[str, ...]]:
    """Apply declared engine-local symbol identities and report those used."""
    applied: set[str] = set()

    def visit(item: object) -> object:
        item = _unwrap(item)
        if isinstance(item, sp.MatrixBase):
            return item.applyfunc(visit)
        if isinstance(item, Mapping):
            return {visit(key): visit(mapped) for key, mapped in item.items()}
        if isinstance(item, list):
            return [visit(mapped) for mapped in item]
        if isinstance(item, tuple):
            return tuple(visit(mapped) for mapped in item)
        if isinstance(item, CasRule):
            return CasRule(visit(item.left), visit(item.right))
        if isinstance(item, CasCall):
            return CasCall(item.head, tuple(visit(mapped) for mapped in item.arguments))
        if isinstance(item, CasRelation):
            return CasRelation(item.operator, visit(item.left), visit(item.right))
        if isinstance(item, CasInequality):
            return CasInequality(tuple(visit(mapped) for mapped in item.operands), item.operators)
        if isinstance(item, CasBoolean):
            return CasBoolean(item.operator, tuple(visit(mapped) for mapped in item.operands))
        if isinstance(item, sp.Basic):
            used = item.free_symbols.intersection(replacements)
            applied.update(f"{symbol}->{replacements[symbol]}" for symbol in used)
            return item.xreplace(replacements)
        return item

    return visit(value), tuple(sorted(applied))


def _tree_values(value: object) -> Iterable[object]:
    value = _unwrap(value)
    yield value
    if isinstance(value, sp.MatrixBase):
        for item in value:
            yield from _tree_values(item)
    elif isinstance(value, Mapping):
        for key, item in value.items():
            yield from _tree_values(key)
            yield from _tree_values(item)
    elif isinstance(value, (list, tuple, set, frozenset)):
        for item in value:
            yield from _tree_values(item)
    elif isinstance(value, CasRule):
        yield from _tree_values(value.left)
        yield from _tree_values(value.right)
    elif isinstance(value, CasCall):
        yield from _tree_values(value.arguments)
    elif isinstance(value, CasRelation):
        yield from _tree_values(value.left)
        yield from _tree_values(value.right)
    elif isinstance(value, CasInequality):
        yield from _tree_values(value.operands)
    elif isinstance(value, CasBoolean):
        yield from _tree_values(value.operands)
    elif isinstance(value, sp.Basic):
        for argument in value.args:
            yield from _tree_values(argument)


def _free_symbol_names(value: object) -> set[str]:
    names: set[str] = set()
    for item in _tree_values(value):
        if isinstance(item, OpaqueDerivative):
            names.add(item.function_name)
            names.update(item.variables)
        elif isinstance(item, sp.Symbol):
            names.add(str(item))
    return names


def _renamable_symbols(value: object) -> set[sp.Symbol]:
    symbols: set[sp.Symbol] = set()
    for item in _tree_values(value):
        if isinstance(item, sp.Basic):
            symbols.update(
                symbol
                for symbol in item.free_symbols
                if not isinstance(symbol, CanonicalDerivative)
            )
    return symbols


def _matchable_tree(value: object) -> sp.Basic:
    """Encode supported comparison containers as one positional SymPy tree."""
    value = _unwrap(value)
    if isinstance(value, sp.MatrixBase):
        return sp.Tuple(sp.Symbol("__matrix__"), value.rows, value.cols, *map(_matchable_tree, value))
    if isinstance(value, Mapping):
        pairs = [sp.Tuple(_matchable_tree(key), _matchable_tree(item)) for key, item in value.items()]
        return sp.Tuple(sp.Symbol("__mapping__"), *sorted(pairs, key=sp.default_sort_key))
    if isinstance(value, (list, tuple)):
        return sp.Tuple(sp.Symbol("__sequence__"), *map(_matchable_tree, value))
    if isinstance(value, CasRule):
        return sp.Tuple(sp.Symbol("__rule__"), _matchable_tree(value.left), _matchable_tree(value.right))
    if isinstance(value, CasCall):
        return sp.Tuple(sp.Symbol(f"__call__{value.head}"), *map(_matchable_tree, value.arguments))
    if isinstance(value, CasRelation):
        return sp.Tuple(sp.Symbol(f"__relation__{value.operator}"), _matchable_tree(value.left), _matchable_tree(value.right))
    if isinstance(value, CasInequality):
        return sp.Tuple(
            sp.Symbol("__inequality__"),
            sp.Tuple(*(sp.Symbol(operator) for operator in value.operators)),
            *map(_matchable_tree, value.operands),
        )
    if isinstance(value, CasBoolean):
        return sp.Tuple(sp.Symbol(f"__boolean__{value.operator}"), *map(_matchable_tree, value.operands))
    if isinstance(value, bool):
        return sp.true if value else sp.false
    if isinstance(value, str):
        return sp.Symbol(f"__string__{value!r}")
    if isinstance(value, int):
        return sp.Integer(value)
    if isinstance(value, sp.Basic):
        return value
    raise ValueError(f"cannot test symbol bijection for {type(value).__name__}")


def _symbol_bijection(
    left: object,
    right: object,
    comparator: Any,
) -> tuple[tuple[str, str], ...] | None:
    """Return an exact leftover-symbol bijection, without declaring it valid."""
    left_symbols = _renamable_symbols(left)
    right_symbols = _renamable_symbols(right)
    shared_names = {str(symbol) for symbol in left_symbols} & {
        str(symbol) for symbol in right_symbols
    }
    left_only = sorted((symbol for symbol in left_symbols if str(symbol) not in shared_names), key=str)
    right_only = sorted((symbol for symbol in right_symbols if str(symbol) not in shared_names), key=str)
    if not left_only or len(left_only) != len(right_only):
        return None
    try:
        encoded_left = _matchable_tree(left)
        encoded_right = _matchable_tree(right)
    except ValueError:
        return None
    wilds = {symbol: sp.Wild(f"__name_{index}") for index, symbol in enumerate(left_only)}
    matched = encoded_right.match(encoded_left.xreplace(wilds))
    if matched is None:
        return None
    replacements: dict[str, str] = {}
    for symbol in left_only:
        candidate = matched.get(wilds[symbol])
        if not isinstance(candidate, sp.Symbol) or candidate not in right_only:
            return None
        replacements[str(symbol)] = str(candidate)
    if len(set(replacements.values())) != len(right_only):
        return None
    mapped, _ = _map_tree(left, replacements)
    if not comparator(mapped, right):
        return None
    return tuple(sorted(replacements.items()))


def _select_shape(value: ParsedValue, selector: str | None) -> object:
    selected = value.value
    if selector is None:
        return selected
    if selector.startswith("free_symbol:"):
        name = selector.partition(":")[2]
        symbols = {symbol for item in _tree_values(selected) if isinstance(item, sp.Basic) for symbol in item.free_symbols}
        matches = [symbol for symbol in symbols if str(symbol) == name]
        if len(matches) != 1:
            raise ValueError(f"expected exactly one free symbol named {name}")
        return matches[0]
    if selector.startswith("pair_value:"):
        name = selector.partition(":")[2]
        sequence = _sequence_form(selected)
        if sequence is None:
            raise ValueError("expected sequence of symbol/value pairs")
        matches: list[object] = []
        for item in sequence:
            pair = _sequence_form(item)
            if pair is not None and len(pair) == 2 and str(_unwrap(pair[0])) == name:
                matches.append(pair[1])
        if len(matches) != 1:
            raise ValueError(f"expected exactly one pair for {name}")
        return matches[0]
    sequence = _sequence_form(selected)
    if selector == "scalar":
        if sequence is not None or isinstance(selected, Mapping):
            raise ValueError("expected scalar")
        return selected
    if selector in {"list", "list_first"}:
        if sequence is None or len(sequence) != 1:
            raise ValueError("expected singleton sequence")
        return sequence[0]
    if selector == "matrix":
        if not isinstance(selected, sp.MatrixBase):
            raise ValueError("expected matrix")
        return selected
    if selector in {"list_of_pairs_second", "last_pair_second"}:
        if sequence is None or not sequence:
            raise ValueError("expected nonempty list of pairs")
        pairs = [_sequence_form(item) for item in sequence]
        if any(pair is None or len(pair) != 2 for pair in pairs):
            raise ValueError("expected list of pairs")
        seconds = [pair[1] for pair in pairs if pair is not None]
        return seconds[-1] if selector == "last_pair_second" else seconds[0] if len(seconds) == 1 else seconds
    if selector == "sequence_third":
        if sequence is None or len(sequence) != 3:
            raise ValueError("expected three-item sequence")
        return sequence[2]
    if selector == "sequence_second":
        if sequence is None or len(sequence) < 2:
            raise ValueError("expected sequence with a second item")
        return sequence[1]
    if selector == "sequence_first":
        if sequence is None or not sequence:
            raise ValueError("expected nonempty sequence")
        return sequence[0]
    raise HarnessError(f"unknown cross-engine selector {selector!r}")


def _normalize_action(value: object, family: str | None) -> object:
    if family != "euler_lagrange":
        return value
    sequence = _sequence_form(value)
    if sequence is None:
        raise ValueError("Euler-Lagrange operand must be an ordered vector")
    residuals: list[object] = []
    for component in sequence:
        if isinstance(component, sp.Equality):
            residuals.append(component.lhs - component.rhs)
        elif isinstance(component, CasRelation) and component.operator == "==":
            if not isinstance(component.left, (int, sp.Basic)) or not isinstance(component.right, (int, sp.Basic)):
                raise ValueError("Euler-Lagrange equation sides must be scalar")
            residuals.append(component.left - component.right)
        elif isinstance(component, sp.Basic) and not _is_boolean_value(component):
            residuals.append(component)
        else:
            raise ValueError("Euler-Lagrange component is neither equation nor residual")
    return residuals


def _cardinality_error(value: object, declaration: object, mode: str) -> str | None:
    if not isinstance(declaration, Mapping) or not isinstance(declaration.get("kind"), str):
        return "missing structural cardinality declaration"
    kind = declaration["kind"]
    sequence = _sequence_form(value)
    if kind == "scalar":
        return None if sequence is None and not isinstance(value, Mapping) else "expected scalar cardinality one"
    if kind in {"sequence", "multiset"}:
        count = declaration.get("count", declaration.get("outer"))
        if not isinstance(count, int) or isinstance(count, bool) or count <= 0:
            return "sequence/multiset cardinality must be an exact positive integer"
        if kind == "multiset" and mode != "multiset":
            return "multiset cardinality requires mode: multiset"
        if sequence is None:
            return f"expected {kind}"
        if not sequence:
            return "empty operand"
        return None if len(sequence) == count else f"expected {count} outer elements, got {len(sequence)}"
    if kind == "mapping":
        count = declaration.get("entries")
        if not isinstance(count, int) or isinstance(count, bool) or count <= 0:
            return "mapping entries must be an exact positive integer"
        if not isinstance(value, Mapping):
            return "expected mapping"
        if not value:
            return "empty operand"
        return None if len(value) == count else f"expected {count} mapping entries, got {len(value)}"
    if kind == "matrix":
        shape = declaration.get("shape")
        if (
            not isinstance(shape, list)
            or len(shape) != 2
            or any(not isinstance(item, int) or isinstance(item, bool) or item <= 0 for item in shape)
        ):
            return "matrix shape must contain two positive integers"
        if not isinstance(value, sp.MatrixBase):
            return "expected matrix"
        return None if [value.rows, value.cols] == shape else f"expected matrix {shape}, got {[value.rows, value.cols]}"
    return f"unknown cardinality kind {kind!r}"


def _identity(row: Mapping[str, Any]) -> tuple[tuple[str, str], ...]:
    reserved = {"quantity", "mode", "cardinality", "family", "dimension", "package", "normalization"}
    return tuple(
        sorted(
            (engine, tag)
            for engine, tag in row.items()
            if engine not in reserved and not engine.endswith("_select") and isinstance(tag, str)
        )
    )


def _reject_unknown_keys(value: Mapping[str, Any], allowed: set[str], site: str) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise HarnessError(f"{site} has unrecognised key(s): {','.join(unknown)}")


def _validate_config_keys(config: Mapping[str, Any]) -> None:
    """Reject misspelled keys wherever a declaration has a fixed vocabulary."""
    _reject_unknown_keys(
        config,
        {
            "default_engine",
            "cells",
            "control",
            "symbol_naming",
            "symbol_identities",
            "dimension_sources",
            "primitive_dimensions",
            "dimensionless",
            "cross_engine",
            "registry_residual",
            "action_families",
            "parity_exclude",
        },
        "config",
    )
    cells = config.get("cells", [])
    if isinstance(cells, list):
        for index, item in enumerate(cells):
            if isinstance(item, Mapping):
                _reject_unknown_keys(item, {"package", "role", "dimension", "stiffness_form"}, f"cells[{index}]")
    control = config.get("control")
    engines: set[str] = set()
    if isinstance(control, Mapping):
        _reject_unknown_keys(control, {"required_suffixes", "tag_templates"}, "control")
        for key in ("required_suffixes", "tag_templates"):
            declared = control.get(key)
            if isinstance(declared, Mapping):
                engines.update(str(engine) for engine in declared)
        templates = control.get("tag_templates")
        if isinstance(templates, Mapping):
            for engine, item in templates.items():
                if isinstance(item, Mapping):
                    _reject_unknown_keys(item, {"main", "control"}, f"control.tag_templates.{engine}")
    naming = config.get("symbol_naming")
    if isinstance(naming, Mapping):
        _reject_unknown_keys(naming, {"rule", "engine_styles", "exceptions"}, "symbol_naming")
        exceptions = naming.get("exceptions", [])
        if isinstance(exceptions, list):
            for index, item in enumerate(exceptions):
                if isinstance(item, Mapping):
                    _reject_unknown_keys(item, {"canonical", "spellings", "reason"}, f"symbol_naming.exceptions[{index}]")
    symbol_identities = config.get("symbol_identities", [])
    if isinstance(symbol_identities, list):
        for index, item in enumerate(symbol_identities):
            if isinstance(item, Mapping):
                _reject_unknown_keys(item, {"engine", "symbol", "expression", "reason"}, f"symbol_identities[{index}]")
    dimension_sources = config.get("dimension_sources", [])
    if isinstance(dimension_sources, list):
        for index, item in enumerate(dimension_sources):
            if not isinstance(item, Mapping):
                continue
            _reject_unknown_keys(item, {"engine", "package", "symbols"}, f"dimension_sources[{index}]")
            symbols = item.get("symbols")
            if isinstance(symbols, Mapping):
                for symbol, source in symbols.items():
                    if isinstance(source, Mapping):
                        _reject_unknown_keys(
                            source,
                            {"tag", "tags", "tag_template", "dimensions", "shape", "select"},
                            f"dimension_sources[{index}].symbols.{symbol}",
                        )
    action_families = config.get("action_families")
    if isinstance(action_families, Mapping):
        _reject_unknown_keys(action_families, {"lagrangian", "euler_lagrange"}, "action_families")
        for family, item in action_families.items():
            if isinstance(item, Mapping) and engines:
                _reject_unknown_keys(item, engines, f"action_families.{family}")
    cross_rows = config.get("cross_engine", [])
    if isinstance(cross_rows, list):
        fixed = {"quantity", "mode", "cardinality", "family", "dimension", "package", "normalization"}
        for index, item in enumerate(cross_rows):
            if not isinstance(item, Mapping):
                continue
            allowed = fixed | engines | {f"{engine}_select" for engine in engines}
            _reject_unknown_keys(item, allowed, f"cross_engine[{index}]")
            cardinality = item.get("cardinality")
            if isinstance(cardinality, Mapping):
                _reject_unknown_keys(cardinality, {"kind", "count", "outer", "entries", "shape"}, f"cross_engine[{index}].cardinality")
    registry_rows = config.get("registry_residual", [])
    if isinstance(registry_rows, list):
        for index, item in enumerate(registry_rows):
            if not isinstance(item, Mapping):
                continue
            _reject_unknown_keys(item, {"relation_id", "engine", "qids"}, f"registry_residual[{index}]")
            qids = item.get("qids")
            if isinstance(qids, Mapping):
                for qid, source in qids.items():
                    if isinstance(source, Mapping):
                        _reject_unknown_keys(source, {"tag", "select"}, f"registry_residual[{index}].qids.{qid}")


def _action_family_declaration_error(
    config: Mapping[str, Any],
    row: Mapping[str, Any],
    identity: tuple[tuple[str, str], ...],
) -> str | None:
    family = row.get("family")
    if family is None:
        return None
    declarations = config.get("action_families")
    declared = declarations.get(family) if isinstance(declarations, Mapping) else None
    if not isinstance(declared, Mapping) or len(declared) < 2 or any(
        not isinstance(engine, str) or not isinstance(suffix, str) or not suffix
        for engine, suffix in declared.items()
    ):
        return f"family {family!r} has no valid declared engine tags"
    package, dimension = row.get("package"), row.get("dimension")
    cells = [
        cell
        for cell in _declared_cells(config)
        if cell.package == package and cell.dimension == dimension
    ]
    if len(cells) != 1:
        return f"family {family!r} row has no unique declared cell for {package}:D{dimension}"
    expected = tuple(
        sorted(
            (engine, _cell_tag(config, engine, cells[0], suffix))
            for engine, suffix in declared.items()
        )
    )
    if identity != expected:
        return f"family {family!r} tags do not match declared identity {expected!r}"
    return None


def _legacy_omega(value: object) -> object:
    replacements = {"omega2": "__legacy_omega_squared", "omegaSquared": "__legacy_omega_squared"}
    mapped, _ = _map_tree(value, replacements)
    if isinstance(mapped, sp.Basic):
        marker = sp.Symbol("__legacy_omega_squared")
        return mapped.xreplace({marker: sp.Symbol("omega") ** 2})
    return mapped


def check_cross_engine(
    rows: Sequence[Mapping[str, Any]],
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
    *,
    config: Mapping[str, Any] | None = None,
    registry: Registry | None = None,
) -> CrossEngineReport:
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
        raise HarnessError("cross_engine must be a list")
    _validate_config_keys(config or {})
    naming = _build_naming_rule(config or {}, registry)
    symbol_identities = _build_symbol_identities(config or {})
    identities = [_identity(row) for row in rows]
    counts = Counter(identities)
    duplicates = tuple(
        f"{identity}: declared {count} times"
        for identity, count in sorted(counts.items(), key=str)
        if identity and count > 1
    )
    seen: set[tuple[tuple[str, str], ...]] = set()
    results: list[CrossEngineRow] = []
    for configured in rows:
        quantity = configured.get("quantity")
        identity = _identity(configured)
        if not isinstance(quantity, str) or not quantity:
            quantity = "<unnamed>"
        if identity in seen:
            continue
        seen.add(identity)
        family = configured.get("family")
        if family not in {None, "lagrangian", "euler_lagrange"}:
            results.append(CrossEngineRow(quantity, identity, "DECLARATION_ERROR", detail=f"unknown family {family!r}"))
            continue
        family_error = _action_family_declaration_error(config or {}, configured, identity)
        if family_error is not None:
            results.append(CrossEngineRow(quantity, identity, "DECLARATION_ERROR", family=family, detail=family_error))
            continue
        if len(identity) < 2:
            results.append(CrossEngineRow(quantity, identity, "DECLARATION_ERROR", family=family, detail="needs at least two engine tags"))
            continue
        mode = configured.get("mode", "symbolic")
        if mode not in {"symbolic", "multiset"}:
            results.append(CrossEngineRow(quantity, identity, "DECLARATION_ERROR", family=family, detail=f"unknown mode {mode!r}"))
            continue
        if family is not None and mode == "multiset":
            results.append(CrossEngineRow(quantity, identity, "DECLARATION_ERROR", family=family, detail="action rows require positional comparison; mode: multiset is forbidden"))
            continue
        missing = tuple(f"{engine}:{tag}" for engine, tag in identity if engine not in outputs or tag not in outputs[engine])
        if missing:
            results.append(CrossEngineRow(quantity, identity, "MISSING", family=family, detail=", ".join(missing)))
            continue
        compared = tuple(outputs[engine][tag] for engine, tag in identity)
        if any(isinstance(value, Unparsed) for value in compared):
            operands = tuple((engine, value.raw if isinstance(value, Unparsed) else str(value.value)) for (engine, _), value in zip(identity, compared))
            results.append(CrossEngineRow(quantity, identity, "UNPARSED", family=family, operands=operands))
            continue
        raw_selected: list[object] = []
        selected: list[object] = []
        naming_applied: list[str] = []
        identities_applied: list[str] = []
        try:
            for (engine, _tag), raw_value in zip(identity, compared):
                assert isinstance(raw_value, ParsedValue)
                value = _select_shape(raw_value, configured.get(f"{engine}_select"))
                value = _normalize_action(value, family)
                raw_selected.append(_canonicalize_derivatives(value))
                value, applied = _map_tree(value, naming.engine_maps.get(engine, {}))
                naming_applied.extend(f"{engine}:{item}" for item in applied)
                value, applied = _apply_symbol_identities(
                    value, symbol_identities.engine_maps.get(engine, {})
                )
                identities_applied.extend(f"{engine}:{item}" for item in applied)
                selected.append(_canonicalize_derivatives(value))
        except ValueError as exc:
            operands = tuple((engine, str(_unwrap(value))) for (engine, _), value in zip(identity, compared))
            results.append(CrossEngineRow(quantity, identity, "SHAPE_MISMATCH", family=family, operands=operands, detail=str(exc)))
            continue
        cardinality_errors = [_cardinality_error(value, configured.get("cardinality"), mode) for value in selected]
        if any(error for error in cardinality_errors):
            detail = "; ".join(f"{engine}: {error}" for (engine, _), error in zip(identity, cardinality_errors) if error)
            status = "EMPTY" if "empty operand" in detail else "CARDINALITY_INVALID"
            results.append(CrossEngineRow(quantity, identity, status, family=family, operands=tuple((engine, str(value)) for (engine, _), value in zip(identity, selected)), detail=detail, naming_applied=tuple(sorted(set(naming_applied))), identities_applied=tuple(sorted(set(identities_applied)))))
            continue
        operands = tuple((engine, str(value)) for (engine, _), value in zip(identity, selected))
        selected_boolean = any(_is_boolean_value(value) for value in selected)
        tree_boolean = any(_is_boolean_value(item) for value in selected for item in _tree_values(value))
        undeclared: set[str] = set()
        shape_ok = all(symbolic_shapes_compatible(selected[0], value) for value in selected[1:])
        comparator = symbolic_multiset_equal if mode == "multiset" else symbolic_equal
        if not shape_ok:
            status = "SHAPE_MISMATCH"
        elif all(comparator(selected[0], value) for value in selected[1:]):
            status = "AGREE"
        else:
            reconciled = True
            for value in selected[1:]:
                if comparator(selected[0], value):
                    continue
                bijection = _symbol_bijection(selected[0], value, comparator)
                if bijection is None:
                    reconciled = False
                    break
                undeclared.update(f"{left}<->{right}" for left, right in bijection)
            status = "NAMING_MISMATCH" if reconciled and undeclared else "DISAGREE"
        legacy_status: str | None = None
        raw_shape_ok = all(
            symbolic_shapes_compatible(raw_selected[0], value)
            for value in raw_selected[1:]
        )
        if raw_shape_ok:
            legacy = [_legacy_omega(value) for value in raw_selected]
            legacy_status = "AGREE" if all(comparator(legacy[0], value) for value in legacy[1:]) else "DISAGREE"
        results.append(
            CrossEngineRow(
                quantity,
                identity,
                status,
                family=family,
                operands=operands,
                naming_applied=tuple(sorted(set(naming_applied))),
                identities_applied=tuple(sorted(set(identities_applied))),
                undeclared_spellings=tuple(sorted(undeclared)),
                selected_boolean=selected_boolean,
                tree_boolean=tree_boolean,
                legacy_algebraic_status=legacy_status,
            )
        )
    declaration_failure = None if rows else "cross_engine declaration is absent or empty"
    gaps = tuple(f"{row.quantity}:{row.status}" for row in results if row.status not in VERDICTS)
    denominator = len({identity for identity in identities if identity})
    numerator = sum(row.status in VERDICTS for row in results)
    selected_boolean_rows = tuple(row.quantity for row in results if row.selected_boolean)
    tree_boolean_rows = tuple(row.quantity for row in results if row.tree_boolean)
    before = sum(row.legacy_algebraic_status == "AGREE" for row in results)
    after = sum(row.status == "AGREE" for row in results)
    changed = tuple(
        row.quantity
        for row in results
        if row.legacy_algebraic_status is not None
        and row.legacy_algebraic_status != row.status
    )
    return CrossEngineReport(
        rows=tuple(results),
        coverage=Coverage(numerator, denominator, gaps, "distinct declared tag-pair identities with AGREE or DISAGREE / distinct declared tag-pair identities"),
        duplicates=duplicates,
        declaration_failure=declaration_failure,
        selected_boolean_rows=selected_boolean_rows,
        tree_boolean_rows=tree_boolean_rows,
        naming_before_agree=before,
        naming_after_agree=after,
        naming_changed_rows=changed,
        naming_exceptions=naming.exceptions,
    )


def action_family_reports(report: CrossEngineReport) -> tuple[ActionFamilyReport, ...]:
    result: list[ActionFamilyReport] = []
    for family in ("lagrangian", "euler_lagrange"):
        rows = [row for row in report.rows if row.family == family]
        gaps = tuple(f"{row.quantity}:{row.status}" for row in rows if row.status not in VERDICTS)
        result.append(
            ActionFamilyReport(
                family,
                Coverage(
                    sum(row.status in VERDICTS for row in rows),
                    len(rows),
                    gaps,
                    "declared action rows with AGREE or DISAGREE / declared action rows",
                ),
                tuple((row.quantity, row.status) for row in rows),
            )
        )
    return tuple(result)


# ---------------------------------------------------------------------------
# Declaration-derived control coverage and cell parity
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DeclaredCell:
    package: str
    role: str
    dimension: int | None
    stiffness_form: str | None

    @property
    def name(self) -> str:
        suffix = f":D{self.dimension}" if self.dimension is not None else ""
        return f"{self.package}{suffix}"


@dataclass(frozen=True)
class ControlRow:
    engine: str
    cell: str
    suffix: str
    main_tag: str
    control_tags: tuple[str, ...]
    status: str
    missing_tags: tuple[str, ...] = ()
    differing_controls: tuple[str, ...] = ()


@dataclass(frozen=True)
class ControlReport:
    engine: str
    rows: tuple[ControlRow, ...]
    coverage: Coverage
    missing_cells: tuple[str, ...]
    declaration_failure: str | None = None

    @property
    def responsive(self) -> tuple[ControlRow, ...]:
        return tuple(row for row in self.rows if row.status == "RESPONSIVE")

    @property
    def invariant(self) -> tuple[ControlRow, ...]:
        return tuple(row for row in self.rows if row.status == "INVARIANT")

    @property
    def partial(self) -> tuple[ControlRow, ...]:
        return tuple(row for row in self.rows if row.status == "PARTIALLY_PAIRED")


@dataclass(frozen=True)
class ParityRow:
    engine: str
    cell: str
    missing: tuple[str, ...]


@dataclass(frozen=True)
class ParityReport:
    engine: str
    rows: tuple[ParityRow, ...]
    declaration_failure: str | None = None


def _declared_cells(config: Mapping[str, Any]) -> tuple[DeclaredCell, ...]:
    raw = config.get("cells")
    if not isinstance(raw, list):
        return ()
    cells: list[DeclaredCell] = []
    for item in raw:
        if not isinstance(item, Mapping):
            raise HarnessError("each cell declaration must be a mapping")
        package, role = item.get("package"), item.get("role")
        dimension, stiffness = item.get("dimension"), item.get("stiffness_form")
        if not isinstance(package, str) or role not in {"main", "control"}:
            raise HarnessError("each cell needs package and main/control role")
        if dimension is not None and (not isinstance(dimension, int) or isinstance(dimension, bool) or dimension <= 0):
            raise HarnessError(f"cell {package} has invalid dimension")
        if stiffness is not None and not isinstance(stiffness, str):
            raise HarnessError(f"cell {package} has invalid stiffness_form")
        cells.append(DeclaredCell(package, role, dimension, stiffness))
    identities = [(cell.package, cell.dimension) for cell in cells]
    if len(set(identities)) != len(identities):
        raise HarnessError("duplicate package/dimension cell declaration")
    return tuple(cells)


def _control_declaration(config: Mapping[str, Any]) -> Mapping[str, Any]:
    declaration = config.get("control")
    if not isinstance(declaration, Mapping):
        return {}
    return declaration


def _required_suffixes(config: Mapping[str, Any], engine: str) -> tuple[str, ...]:
    raw = _control_declaration(config).get("required_suffixes")
    if not isinstance(raw, Mapping) or not isinstance(raw.get(engine), list):
        return ()
    values = raw[engine]
    if any(not isinstance(value, str) or not value for value in values):
        raise HarnessError(f"control.required_suffixes.{engine} must contain names")
    return tuple(dict.fromkeys(values))


def _cell_tag(config: Mapping[str, Any], engine: str, cell: DeclaredCell, suffix: str) -> str:
    raw = _control_declaration(config).get("tag_templates")
    if not isinstance(raw, Mapping) or not isinstance(raw.get(engine), Mapping):
        raise HarnessError(f"missing control tag template for engine {engine}")
    templates = raw[engine]
    template = templates.get(cell.role)
    if not isinstance(template, str):
        raise HarnessError(f"missing {cell.role} tag template for engine {engine}")
    try:
        return template.format(
            package=cell.package,
            dimension=cell.dimension if cell.dimension is not None else "",
            suffix=suffix,
        )
    except (KeyError, ValueError) as exc:
        raise HarnessError(f"invalid tag template for engine {engine}: {exc}") from exc


def _control_groups(cells: Sequence[DeclaredCell]) -> tuple[tuple[DeclaredCell, tuple[DeclaredCell, ...]], ...]:
    result: list[tuple[DeclaredCell, tuple[DeclaredCell, ...]]] = []
    for main in (cell for cell in cells if cell.role == "main"):
        controls = tuple(
            cell
            for cell in cells
            if cell.role == "control" and cell.dimension == main.dimension
        )
        if controls:
            result.append((main, controls))
    return tuple(result)


def check_control_response(
    values: Mapping[str, ParsedValue | Unparsed],
    config: Mapping[str, Any] | None = None,
    *,
    engine: str = "default",
) -> ControlReport:
    config = config or {}
    cells = _declared_cells(config)
    suffixes = _required_suffixes(config, engine)
    groups = _control_groups(cells)
    declaration_failure = None
    if not cells or not suffixes or not groups:
        declaration_failure = f"{engine} control declaration has no assessable declared cell"
    rows: list[ControlRow] = []
    for main, controls in groups:
        for suffix in suffixes:
            main_tag = _cell_tag(config, engine, main, suffix)
            control_tags = tuple(_cell_tag(config, engine, cell, suffix) for cell in controls)
            present_main = main_tag in values
            present_controls = tuple(tag in values for tag in control_tags)
            missing = tuple(
                tag
                for tag, present in ((main_tag, present_main), *zip(control_tags, present_controls))
                if not present
            )
            if not present_main and not any(present_controls):
                status = "MISSING_CELL"
            elif missing:
                status = "PARTIALLY_PAIRED"
            else:
                operands = (values[main_tag], *(values[tag] for tag in control_tags))
                if any(isinstance(value, Unparsed) for value in operands):
                    status = "UNPARSED"
                else:
                    differing = tuple(
                        tag
                        for tag in control_tags
                        if not symbolic_equal(values[main_tag], values[tag])
                    )
                    status = "RESPONSIVE" if differing else "INVARIANT"
                    rows.append(
                        ControlRow(engine, main.name, suffix, main_tag, control_tags, status, differing_controls=differing)
                    )
                    continue
            rows.append(ControlRow(engine, main.name, suffix, main_tag, control_tags, status, missing_tags=missing))
    gaps = tuple(f"{row.cell}:{row.suffix}:{row.status}" for row in rows if row.status not in {"RESPONSIVE", "INVARIANT"})
    coverage = Coverage(
        sum(row.status in {"RESPONSIVE", "INVARIANT"} for row in rows),
        len(rows),
        gaps,
        "fully paired declared control rows with RESPONSIVE or INVARIANT / declared control rows",
    )
    observed_cells = {
        cell.name
        for cell in cells
        if any(_cell_tag(config, engine, cell, suffix) in values for suffix in suffixes)
    } if suffixes else set()
    missing_cells = tuple(sorted(cell.name for cell in cells if cell.name not in observed_cells))
    return ControlReport(engine, tuple(rows), coverage, missing_cells, declaration_failure)


def check_tag_parity(
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
    config: Mapping[str, Any] | Sequence[str] | None = None,
) -> tuple[ParityReport, ...] | ParityReport:
    """Check required tag presence at every declared cell.

    The Sequence compatibility path is retained for small callers predating
    cell declarations; the rebuilt configs always take the mapping path.
    """
    if config is not None and not isinstance(config, Mapping):
        # Historical S9 helper: report inferred X-package suffix parity.
        rows: list[ParityRow] = []
        for engine, values in outputs.items():
            prefix = "WL_S9_" if engine == "wl" else "PY_S9_MAIN_"
            main = {tag[len(prefix):] for tag in values if tag.startswith(prefix)}
            control_pattern = re.compile(rf"^(?:WL|PY)_S9_(X[1-9][0-9]*)_(.+)$")
            grouped: dict[str, set[str]] = {}
            for tag in values:
                match = control_pattern.match(tag)
                if match:
                    grouped.setdefault(match.group(1), set()).add(match.group(2))
            rows.extend(ParityRow(engine, package, tuple(sorted(main - suffixes))) for package, suffixes in grouped.items())
        return ParityReport("all", tuple(rows))
    declaration = config or {}
    cells = _declared_cells(declaration)
    reports: list[ParityReport] = []
    for engine, values in sorted(outputs.items()):
        suffixes = _required_suffixes(declaration, engine)
        if not cells or not suffixes:
            reports.append(ParityReport(engine, (), f"{engine} parity declaration has no cells or suffixes"))
            continue
        rows = []
        for cell in cells:
            missing = tuple(
                suffix
                for suffix in suffixes
                if _cell_tag(declaration, engine, cell, suffix) not in values
            )
            rows.append(ParityRow(engine, cell.name, missing))
        reports.append(ParityReport(engine, tuple(rows)))
    return tuple(reports)


# ---------------------------------------------------------------------------
# Per-engine, per-package dimension tables and partition accounting
# ---------------------------------------------------------------------------


Dimension = tuple[sp.Expr, sp.Expr, sp.Expr]


@dataclass(frozen=True)
class TermDimension:
    expression: str
    dimension: Dimension | None


@dataclass(frozen=True)
class NonHomogeneous:
    tag: str
    expression: str
    site: str
    summands: tuple[TermDimension, ...]


@dataclass(frozen=True)
class DimensionTableReport:
    engine: str
    package: str
    assessable: bool
    vectors: tuple[tuple[str, Dimension], ...]
    sources: tuple[str, ...]
    failures: tuple[str, ...]


@dataclass(frozen=True)
class DimensionTagStatus:
    tag: str
    package: str | None
    status: str
    propositions: tuple[tuple[str, int], ...] = ()
    detail: str = ""


@dataclass(frozen=True)
class DimensionReport:
    engine: str
    total_tags: int
    compared: int
    no_comparison: int
    not_applicable: int
    unwalked: int
    unassessable: int
    unparsed: int
    homogeneous: int
    statuses: tuple[DimensionTagStatus, ...]
    non_homogeneous: tuple[NonHomogeneous, ...]
    tables: tuple[DimensionTableReport, ...]
    proposition_sites: tuple[tuple[str, int], ...]
    package_disagreements: tuple[str, ...]
    coverage: Coverage
    declaration_failure: str | None

    @property
    def checked_tags(self) -> int:  # compatibility name
        return self.compared

    @property
    def homogeneous_tags(self) -> int:
        return self.homogeneous

    @property
    def unknown_symbols(self) -> tuple[tuple[str, str], ...]:
        result: list[tuple[str, str]] = []
        for row in self.statuses:
            if row.status == "unassessable" and "UNKNOWN_SYMBOL" in row.detail:
                for name in row.detail.partition("UNKNOWN_SYMBOL=")[2].split(","):
                    if name:
                        result.append((row.tag, name))
        return tuple(result)


@lru_cache(maxsize=4096)
def _dimension_equal(left: Dimension, right: Dimension) -> bool:
    return all(sp.simplify(a - b) == 0 for a, b in zip(left, right))


@lru_cache(maxsize=4096)
def _dimension_add(left: Dimension, right: Dimension) -> Dimension:
    return tuple(sp.simplify(a + b) for a, b in zip(left, right))  # type: ignore[return-value]


@lru_cache(maxsize=4096)
def _dimension_scale(value: Dimension, factor: sp.Expr) -> Dimension:
    return tuple(sp.simplify(factor * item) for item in value)  # type: ignore[return-value]


def _as_vector(value: object, label: str) -> Dimension:
    sequence = _sequence_form(value)
    if sequence is None or len(sequence) != 3:
        raise ValueError(f"{label}: expected exactly three scalar components")
    result: list[sp.Expr] = []
    for component in sequence:
        component = _unwrap(component)
        if isinstance(component, bool) or not isinstance(component, (int, sp.Basic)):
            raise ValueError(f"{label}: nonscalar component {component!r}")
        result.append(sp.Integer(component) if isinstance(component, int) else component)
    return tuple(result)  # type: ignore[return-value]


def _dimension_source_value(value: object, shape: str, label: str) -> Dimension:
    value = _unwrap(value)
    if shape == "vector":
        return _as_vector(value, label)
    family = _sequence_form(value)
    if shape == "family":
        if family is None or not family:
            raise ValueError(f"{label}: family must be nonempty")
        vectors = [_as_vector(member, label) for member in family]
    elif shape == "nested_family":
        if family is None or not family:
            raise ValueError(f"{label}: nested_family must be nonempty")
        vectors = []
        for member in family:
            singleton = _sequence_form(member)
            if singleton is None or len(singleton) != 1:
                raise ValueError(f"{label}: nested_family member must be a singleton sequence")
            vectors.append(_as_vector(singleton[0], label))
    else:
        raise ValueError(f"{label}: unknown declared shape {shape!r}")
    if any(not _dimension_equal(vectors[0], vector) for vector in vectors[1:]):
        raise ValueError(f"{label}: family members disagree: {vectors}")
    return vectors[0]


def _expand_source_tags(source: Mapping[str, Any], package: str) -> tuple[str, ...]:
    if isinstance(source.get("tag"), str):
        return (source["tag"],)
    if isinstance(source.get("tags"), list) and all(isinstance(tag, str) for tag in source["tags"]):
        return tuple(source["tags"])
    template, dimensions = source.get("tag_template"), source.get("dimensions")
    if isinstance(template, str) and isinstance(dimensions, list) and all(isinstance(value, int) and not isinstance(value, bool) for value in dimensions):
        return tuple(template.format(package=package, dimension=dimension) for dimension in dimensions)
    raise HarnessError(f"dimension source for {package} needs tag, tags, or tag_template+dimensions")


def _primitive_table(config: Mapping[str, Any]) -> dict[str, Dimension]:
    table: dict[str, Dimension] = {}
    raw = config.get("primitive_dimensions", {})
    if not isinstance(raw, Mapping):
        raise HarnessError("primitive_dimensions must be a mapping")
    for symbol, vector in raw.items():
        if not isinstance(symbol, str) or not isinstance(vector, list):
            raise HarnessError(f"invalid primitive dimension {symbol!r}")
        table[symbol] = _as_vector(vector, f"primitive {symbol}")
    dimensionless = config.get("dimensionless", [])
    if not isinstance(dimensionless, list) or any(not isinstance(name, str) for name in dimensionless):
        raise HarnessError("dimensionless must be a list")
    table.update({name: ZERO_DIMENSION for name in dimensionless})
    return table


def _build_dimension_tables(
    engine: str,
    values: Mapping[str, ParsedValue | Unparsed],
    config: Mapping[str, Any],
) -> tuple[dict[str, dict[str, Dimension]], tuple[DimensionTableReport, ...], tuple[str, ...]]:
    raw = config.get("dimension_sources")
    if not isinstance(raw, list):
        return {}, (), (f"{engine}: dimension_sources declaration absent",)
    primitives = _primitive_table(config)
    tables: dict[str, dict[str, Dimension]] = {}
    reports: list[DimensionTableReport] = []
    for declaration in raw:
        if not isinstance(declaration, Mapping) or declaration.get("engine") != engine:
            continue
        package, symbols = declaration.get("package"), declaration.get("symbols")
        if not isinstance(package, str) or not isinstance(symbols, Mapping) or not symbols:
            raise HarnessError(f"{engine} dimension source needs package and symbols")
        table = dict(primitives)
        vectors: list[tuple[str, Dimension]] = []
        sources: list[str] = []
        failures: list[str] = []
        for symbol, source_raw in symbols.items():
            if not isinstance(symbol, str):
                raise HarnessError(f"{engine}:{package} has non-string dimension symbol")
            source: Mapping[str, Any]
            if isinstance(source_raw, str):
                source = {"tag": source_raw, "shape": "vector"}
            elif isinstance(source_raw, Mapping):
                source = source_raw
            else:
                failures.append(f"{symbol}: invalid source declaration")
                continue
            shape = source.get("shape", "vector")
            if shape not in {"vector", "family", "nested_family"}:
                failures.append(f"{symbol}: invalid shape {shape!r}")
                continue
            tags = _expand_source_tags(source, package)
            sources.extend(tags)
            observed: list[Dimension] = []
            for tag in tags:
                if tag not in values:
                    failures.append(f"{symbol}:{tag}: MISSING")
                    continue
                value = values[tag]
                if isinstance(value, Unparsed):
                    failures.append(f"{symbol}:{tag}: UNPARSED {value.error}")
                    continue
                try:
                    selected = _select_shape(value, source.get("select"))
                    observed.append(_dimension_source_value(selected, str(shape), tag))
                except ValueError as exc:
                    failures.append(f"{symbol}:{tag}: SHAPE {exc}")
            if observed and any(not _dimension_equal(observed[0], item) for item in observed[1:]):
                failures.append(f"{symbol}: source-tag family disagrees: {observed}")
            elif len(observed) == len(tags) and observed:
                table[symbol] = observed[0]
                vectors.append((symbol, observed[0]))
        assessable = not failures and len(vectors) == len(symbols)
        if assessable:
            tables[package] = table
        reports.append(DimensionTableReport(engine, package, assessable, tuple(vectors), tuple(sources), tuple(failures)))
    disagreements: list[str] = []
    assessable_reports = [report for report in reports if report.assessable]
    symbols = sorted({symbol for report in assessable_reports for symbol, _ in report.vectors})
    for symbol in symbols:
        by_vector: dict[str, list[str]] = {}
        for report in assessable_reports:
            mapping = dict(report.vectors)
            if symbol in mapping:
                by_vector.setdefault(str(mapping[symbol]), []).append(report.package)
        if len(by_vector) > 1:
            detail = "; ".join(f"{vector}: {','.join(packages)}" for vector, packages in sorted(by_vector.items()))
            disagreements.append(f"{engine}:{symbol}: {detail}")
    if not reports:
        disagreements.append(f"{engine}: no declared package dimension table")
    return tables, tuple(reports), tuple(disagreements)


def _tag_cell(config: Mapping[str, Any], engine: str, tag: str) -> DeclaredCell | None:
    candidates: list[tuple[int, DeclaredCell]] = []
    for cell in _declared_cells(config):
        try:
            marker = _cell_tag(config, engine, cell, "__SUFFIX__")
        except HarnessError:
            continue
        prefix, _, suffix = marker.partition("__SUFFIX__")
        if tag.startswith(prefix) and tag.endswith(suffix):
            candidates.append((len(prefix) + len(suffix), cell))
    return max(candidates, default=(0, None), key=lambda item: item[0])[1]  # type: ignore[return-value]


class _DimensionWalker:
    def __init__(self, dimensions: Mapping[str, Dimension], tag: str) -> None:
        self.dimensions = dimensions
        self.tag = tag
        self.unknown: set[str] = set()
        self.sites: Counter[str] = Counter()
        self.non_homogeneous: list[NonHomogeneous] = []

    def walk(self, value: object) -> None:
        value = _unwrap(value)
        if isinstance(value, str):
            return
        if isinstance(value, Mapping):
            for key, item in value.items():
                self.walk(key)
                self.walk(item)
            return
        if isinstance(value, sp.MatrixBase):
            for item in value:
                self.walk(item)
            return
        if isinstance(value, (list, tuple, set, frozenset)):
            for item in value:
                self.walk(item)
            return
        if isinstance(value, CasRule):
            self.walk(value.left)
            self.walk(value.right)
            return
        if isinstance(value, CasCall):
            if value.head == "ConditionalExpression":
                self.walk(value.arguments[0])
                raise DimensionError("ConditionalExpression condition requires a dimensional rule")
            raise UnwalkedDimension(f"unsupported marker/call path {value.head}")
        if isinstance(value, CasInequality):
            for left, right, operator in zip(value.operands, value.operands[1:], value.operators):
                self._compare_relation(left, right, f"inequality:{operator}")
            return
        if isinstance(value, CasBoolean):
            for operand in value.operands:
                self.walk(operand)
            return
        if isinstance(value, CasRelation):
            self._compare_relation(value.left, value.right, "relation")
            return
        if isinstance(value, sp.Contains):
            return  # Membership is semantic metadata, not homogeneity.
        if isinstance(value, sp.assumptions.assume.AppliedPredicate):
            return
        if isinstance(value, sp.logic.boolalg.BooleanFunction):
            for argument in value.args:
                self.walk(argument)
            return
        if isinstance(value, sp.core.relational.Relational):
            self._compare_relation(value.lhs, value.rhs, "relation")
            return
        if isinstance(value, bool):
            return
        if isinstance(value, (int, sp.Basic)):
            self.dimension(sp.Integer(value) if isinstance(value, int) else value)
            return
        raise UnwalkedDimension(f"unsupported normalized value {type(value).__name__}")

    def _compare_relation(self, left: object, right: object, site: str) -> None:
        if not isinstance(left, (int, sp.Basic)) or isinstance(left, bool) or not isinstance(right, (int, sp.Basic)) or isinstance(right, bool):
            raise UnwalkedDimension(f"{site} has non-scalar sides")
        lhs = sp.Integer(left) if isinstance(left, int) else left
        rhs = sp.Integer(right) if isinstance(right, int) else right
        left_dimension, right_dimension = self.dimension(lhs), self.dimension(rhs)
        if lhs == 0 or rhs == 0 or left_dimension is None or right_dimension is None:
            return
        self.sites[site] += 1
        if not _dimension_equal(left_dimension, right_dimension):
            self._record(f"{left} ? {right}", site, (left, right), (left_dimension, right_dimension))

    def dimension(self, expression: sp.Expr) -> Dimension | None:
        if isinstance(expression, OpaqueDerivative):
            result = self._named(expression.function_name)
            variables = [self._named(name) for name in expression.variables]
            if result is None or any(value is None for value in variables):
                return None
            for order, variable in zip(expression.orders, variables):
                assert variable is not None
                result = _dimension_add(result, _dimension_scale(variable, -sp.Integer(order)))
            return result
        if isinstance(expression, sp.Derivative):
            base = expression.expr
            if getattr(base, "is_Function", False):
                name = getattr(base.func, "__name__", str(base.func))
                result = self._named(name)
            else:
                result = self.dimension(base)
            variables = [(variable, count) for variable, count in expression.variable_count]
            variable_dimensions = [self.dimension(variable) for variable, _ in variables]
            if result is None or any(value is None for value in variable_dimensions):
                return None
            for (_, count), variable_dimension in zip(variables, variable_dimensions):
                assert variable_dimension is not None
                result = _dimension_add(result, _dimension_scale(variable_dimension, -sp.Integer(count)))
            return result
        if isinstance(expression, sp.Symbol):
            return self._named(str(expression))
        if expression.is_number is True:
            return ZERO_DIMENSION
        if isinstance(expression, sp.Add):
            terms = expression.args
            dimensions = tuple(self.dimension(term) for term in terms)
            if any(value is None for value in dimensions):
                return None
            known = tuple(value for value in dimensions if value is not None)
            if len(known) > 1:
                self.sites["add"] += 1
                if any(not _dimension_equal(known[0], value) for value in known[1:]):
                    self._record(expression, "add", terms, dimensions)
            return known[0] if known else ZERO_DIMENSION
        if isinstance(expression, sp.Mul):
            result = ZERO_DIMENSION
            dimensions = [self.dimension(factor) for factor in expression.args]
            if any(value is None for value in dimensions):
                return None
            for value in dimensions:
                assert value is not None
                result = _dimension_add(result, value)
            return result
        if isinstance(expression, sp.Pow):
            base, exponent = expression.args
            base_dimension = self.dimension(base)
            if exponent.is_number is not True:
                exponent_dimension = self.dimension(exponent)
                if exponent_dimension is not None and not _dimension_equal(exponent_dimension, ZERO_DIMENSION):
                    self.non_homogeneous.append(
                        NonHomogeneous(self.tag, str(expression), "power_exponent", (TermDimension(str(exponent), exponent_dimension),))
                    )
            if exponent.is_Rational is not True:
                raise DimensionError(f"non-rational exponent {exponent}")
            return None if base_dimension is None else _dimension_scale(base_dimension, exponent)
        if expression.func == sp.exp:
            argument = self.dimension(expression.args[0])
            if argument is not None and not _dimension_equal(argument, ZERO_DIMENSION):
                self.non_homogeneous.append(
                    NonHomogeneous(self.tag, str(expression), "function_argument", (TermDimension(str(expression.args[0]), argument),))
                )
            return ZERO_DIMENSION
        if expression.is_Function:
            for argument in expression.args:
                self.dimension(argument)
            return self._named(getattr(expression.func, "__name__", str(expression.func)))
        raise DimensionError(f"no dimension rule for {type(expression).__name__}: {expression}")

    def _named(self, name: str) -> Dimension | None:
        value = self.dimensions.get(name)
        if value is None:
            self.unknown.add(name)
        return value

    def _record(
        self,
        expression: object,
        site: str,
        terms: Sequence[object],
        dimensions: Sequence[Dimension | None],
    ) -> None:
        self.non_homogeneous.append(
            NonHomogeneous(
                self.tag,
                str(expression),
                site,
                tuple(TermDimension(str(term), dimension) for term, dimension in zip(terms, dimensions)),
            )
        )


def check_dimensions(
    values: Mapping[str, ParsedValue | Unparsed],
    config: Mapping[str, Any],
    *,
    engine: str | None = None,
) -> DimensionReport:
    engine = engine or str(config.get("default_engine", "default"))
    tables, table_reports, package_disagreements = _build_dimension_tables(engine, values, config)
    statuses: list[DimensionTagStatus] = []
    findings: list[NonHomogeneous] = []
    site_counts: Counter[str] = Counter()
    homogeneous = 0
    for tag, parsed in values.items():
        cell = _tag_cell(config, engine, tag)
        package = cell.package if cell is not None else None
        if isinstance(parsed, Unparsed):
            statuses.append(DimensionTagStatus(tag, package, "unparsed", detail=parsed.error))
            continue
        if package is None:
            statuses.append(DimensionTagStatus(tag, None, "not_applicable", detail="no declared package identity"))
            continue
        if parsed.kind == ValueKind.STRING or (
            parsed.kind == ValueKind.CALL
            and not (isinstance(parsed.value, CasCall) and parsed.value.head == "ConditionalExpression")
        ):
            statuses.append(DimensionTagStatus(tag, package, "not_applicable", detail=f"{parsed.kind.value} has no dimensional semantics"))
            continue
        table = tables.get(package)
        if table is None:
            statuses.append(DimensionTagStatus(tag, package, "unassessable", detail="required engine-package dimension table unavailable"))
            continue
        walker = _DimensionWalker(table, tag)
        try:
            walker.walk(parsed.value)
        except UnwalkedDimension as exc:
            statuses.append(DimensionTagStatus(tag, package, "unwalked", detail=str(exc)))
            continue
        except DimensionError as exc:
            statuses.append(DimensionTagStatus(tag, package, "unassessable", detail=str(exc)))
            continue
        if walker.unknown:
            statuses.append(DimensionTagStatus(tag, package, "unassessable", tuple(sorted(walker.sites.items())), "UNKNOWN_SYMBOL=" + ",".join(sorted(walker.unknown))))
            continue
        findings.extend(walker.non_homogeneous)
        if walker.sites:
            statuses.append(DimensionTagStatus(tag, package, "compared", tuple(sorted(walker.sites.items()))))
            site_counts.update(walker.sites)
            if not walker.non_homogeneous:
                homogeneous += 1
        else:
            statuses.append(DimensionTagStatus(tag, package, "no_comparison"))
    counts = Counter(row.status for row in statuses)
    source_declarations = config.get("dimension_sources")
    expected_packages = {
        item["package"]
        for item in source_declarations
        if isinstance(item, Mapping)
        and item.get("engine") == engine
        and isinstance(item.get("package"), str)
    } if isinstance(source_declarations, list) else set()
    compared_packages = {
        row.package for row in statuses if row.status == "compared" and row.package is not None
    }
    coverage_gaps = tuple(sorted(expected_packages - compared_packages))
    declaration_failure = (
        None
        if expected_packages
        else f"{engine}: dimension_sources declaration is absent or empty"
    )
    return DimensionReport(
        engine,
        len(values),
        counts["compared"],
        counts["no_comparison"],
        counts["not_applicable"],
        counts["unwalked"],
        counts["unassessable"],
        counts["unparsed"],
        homogeneous,
        tuple(statuses),
        tuple(findings),
        table_reports,
        tuple(sorted(site_counts.items())),
        package_disagreements,
        Coverage(
            len(expected_packages & compared_packages),
            len(expected_packages),
            coverage_gaps,
            "declared engine-package dimension tables with at least one comparison / declared engine-package dimension tables",
        ),
        declaration_failure,
    )


# ---------------------------------------------------------------------------
# Registry residual declaration coverage and top-level report
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RegistryResidualRow:
    identity: str
    relation_id: str
    status: str
    residual: sp.Expr | None = None
    detail: str = ""
    naming_applied: tuple[str, ...] = ()


@dataclass(frozen=True)
class RegistryResidualReport:
    rows: tuple[RegistryResidualRow, ...]
    coverage: Coverage
    duplicates: tuple[str, ...]
    declaration_failure: str | None


def _registry_identity(row: Mapping[str, Any], default_engine: str) -> str:
    engine = row.get("engine", default_engine)
    relation = row.get("relation_id")
    qids = row.get("qids", row.get("inputs"))
    rendered: list[tuple[str, str]] = []
    if isinstance(qids, Mapping):
        for qid, source in qids.items():
            tag = source if isinstance(source, str) else source.get("tag") if isinstance(source, Mapping) else None
            rendered.append((str(qid), str(tag)))
    return repr((engine, relation, tuple(sorted(rendered))))


def check_registry_residuals(
    rows: Sequence[Mapping[str, Any]],
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
    registry: Registry,
    default_engine: str,
    *,
    config: Mapping[str, Any] | None = None,
) -> RegistryResidualReport:
    identities = [_registry_identity(row, default_engine) for row in rows]
    counts = Counter(identities)
    duplicates = tuple(f"{identity}: declared {count} times" for identity, count in sorted(counts.items()) if count > 1)
    naming = _build_naming_rule(config or {}, registry)
    seen: set[str] = set()
    results: list[RegistryResidualRow] = []
    for row in rows:
        identity = _registry_identity(row, default_engine)
        if identity in seen:
            continue
        seen.add(identity)
        relation_id, engine = row.get("relation_id"), row.get("engine", default_engine)
        qids = row.get("qids", row.get("inputs"))
        if not isinstance(relation_id, str) or relation_id not in registry.relations:
            results.append(RegistryResidualRow(identity, str(relation_id), "DECLARATION_ERROR", detail="unknown relation"))
            continue
        if not isinstance(engine, str) or engine not in outputs:
            results.append(RegistryResidualRow(identity, relation_id, "MISSING", detail=f"engine {engine!r}"))
            continue
        if not isinstance(qids, Mapping) or not qids:
            results.append(RegistryResidualRow(identity, relation_id, "DECLARATION_ERROR", detail="qids mapping absent or empty"))
            continue
        relation = registry.relations[relation_id]
        if relation.residual is None:
            results.append(RegistryResidualRow(identity, relation_id, "DECLARATION_ERROR", detail="registry relation has no residual"))
            continue
        canonical_sources: dict[str, object] = {}
        declaration_error = ""
        for qid, source in qids.items():
            if not isinstance(qid, str):
                declaration_error = "QID names must be strings"
                break
            canonical = registry.resolve_qid(qid)
            if canonical in canonical_sources:
                declaration_error = f"QID {canonical} mapped more than once"
                break
            canonical_sources[canonical] = source
        if declaration_error:
            results.append(RegistryResidualRow(identity, relation_id, "DECLARATION_ERROR", detail=declaration_error))
            continue
        involved = {qid for qid, symbol in registry.symbols.items() if symbol in relation.residual.free_symbols}
        missing_qids = sorted(involved - set(canonical_sources))
        if missing_qids:
            results.append(RegistryResidualRow(identity, relation_id, "MISSING", detail=f"qids={missing_qids}"))
            continue
        substitutions: dict[sp.Symbol, object] = {}
        naming_applied: list[str] = []
        status = ""
        detail = ""
        for qid, source in canonical_sources.items():
            tag = source if isinstance(source, str) else source.get("tag") if isinstance(source, Mapping) else None
            selector = source.get("select") if isinstance(source, Mapping) else None
            if not isinstance(tag, str):
                status, detail = "DECLARATION_ERROR", f"{qid}: source needs tag"
                break
            if tag not in outputs[engine]:
                status, detail = "MISSING", f"{qid}:{tag}"
                break
            value = outputs[engine][tag]
            if isinstance(value, Unparsed):
                status, detail = "UNPARSED", f"{qid}:{tag}: {value.error}"
                break
            assert isinstance(value, ParsedValue)
            try:
                selected = _select_shape(value, selector)
            except ValueError as exc:
                status, detail = "INVALID_SHAPE", f"{qid}:{tag}: {exc}"
                break
            selected, applied = _map_tree(selected, naming.engine_maps.get(engine, {}))
            naming_applied.extend(f"{engine}:{item}" for item in applied)
            if isinstance(selected, bool) or not isinstance(selected, (int, sp.Basic)) or _is_boolean_value(selected):
                status, detail = "INVALID_SHAPE", f"{qid}:{tag}: substitution must be scalar"
                break
            substitutions[registry.symbols[qid]] = sp.Integer(selected) if isinstance(selected, int) else selected
        if status:
            results.append(
                RegistryResidualRow(
                    identity,
                    relation_id,
                    status,
                    detail=detail,
                    naming_applied=tuple(sorted(set(naming_applied))),
                )
            )
            continue
        try:
            residual = sp.simplify(relation.residual.subs(substitutions))
        except (TypeError, ValueError) as exc:
            results.append(RegistryResidualRow(identity, relation_id, "UNSUBSTITUTABLE", detail=str(exc)))
            continue
        results.append(
            RegistryResidualRow(
                identity,
                relation_id,
                "ZERO" if residual == 0 else "NONZERO",
                residual,
                naming_applied=tuple(sorted(set(naming_applied))),
            )
        )
    declaration_failure = None if rows else "registry_residual declaration is absent or empty"
    gaps = tuple(f"{row.relation_id}:{row.status}" for row in results if row.status not in REGISTRY_VERDICTS)
    return RegistryResidualReport(
        tuple(results),
        Coverage(
            sum(row.status in REGISTRY_VERDICTS for row in results),
            len(set(identities)),
            gaps,
            "distinct declared registry identities with ZERO or NONZERO / distinct declared registry identities",
        ),
        duplicates,
        declaration_failure,
    )


@dataclass(frozen=True)
class IgnoredOutputReport:
    engine: str
    ignored_lines: tuple[str, ...]
    duplicate_tags: tuple[str, ...]


@dataclass(frozen=True)
class HarnessReport:
    controls: tuple[ControlReport, ...]
    parity: tuple[ParityReport, ...]
    dimensions: tuple[DimensionReport, ...]
    cross_engine: CrossEngineReport
    registry: RegistryResidualReport
    actions: tuple[ActionFamilyReport, ...]
    output_diagnostics: tuple[IgnoredOutputReport, ...]

    @property
    def operational_failures(self) -> tuple[str, ...]:
        failures: list[str] = []
        for diagnostic in self.output_diagnostics:
            if diagnostic.ignored_lines:
                failures.append(f"{diagnostic.engine}: {len(diagnostic.ignored_lines)} ignored non-tag line(s)")
            if diagnostic.duplicate_tags:
                failures.append(f"{diagnostic.engine}: duplicate emitted tags {','.join(diagnostic.duplicate_tags)}")
        for control in self.controls:
            if control.declaration_failure:
                failures.append(control.declaration_failure)
            if control.coverage.gaps:
                failures.append(f"{control.engine}: control coverage gaps {','.join(control.coverage.gaps[:12])}")
            if control.missing_cells:
                failures.append(f"{control.engine}: missing declared cells {','.join(control.missing_cells)}")
        for parity in self.parity:
            if parity.declaration_failure:
                failures.append(parity.declaration_failure)
            missing = [f"{row.cell}:{','.join(row.missing)}" for row in parity.rows if row.missing]
            if missing:
                failures.append(f"{parity.engine}: parity missing {','.join(missing[:12])}")
        for dimension in self.dimensions:
            if dimension.declaration_failure:
                failures.append(dimension.declaration_failure)
            if dimension.coverage.gaps:
                failures.append(
                    f"{dimension.engine}: dimension coverage gaps {','.join(dimension.coverage.gaps)}"
                )
            if dimension.coverage.denominator and dimension.compared == 0:
                failures.append(
                    f"{dimension.engine}: zero dimension comparisons against non-empty declaration"
                )
            bad_tables = [f"{row.package}:{'|'.join(row.failures)}" for row in dimension.tables if not row.assessable]
            if bad_tables:
                failures.append(f"{dimension.engine}: invalid dimension tables {';'.join(bad_tables)}")
            bad_statuses = [f"{row.tag}:{row.status}" for row in dimension.statuses if row.status in {"unwalked", "unassessable", "unparsed"}]
            if bad_statuses:
                failures.append(f"{dimension.engine}: dimension operational rows {','.join(bad_statuses[:20])}")
        if self.cross_engine.declaration_failure:
            failures.append(self.cross_engine.declaration_failure)
        failures.extend(f"duplicate cross-engine declaration {item}" for item in self.cross_engine.duplicates)
        failures.extend(f"cross-engine {row.quantity}: {row.status}" for row in self.cross_engine.rows if row.status in OPERATIONAL_CROSS_STATUSES)
        if self.cross_engine.tree_boolean_rows:
            failures.append("tree boolean exposure in " + ",".join(self.cross_engine.tree_boolean_rows))
        action_required = any(action.coverage.denominator for action in self.actions)
        for action in self.actions:
            if action_required and (action.coverage.denominator == 0 or action.coverage.numerator == 0):
                failures.append(f"action {action.family} has no comparison verdict")
        if self.registry.declaration_failure:
            failures.append(self.registry.declaration_failure)
        failures.extend(f"duplicate registry declaration {item}" for item in self.registry.duplicates)
        failures.extend(f"registry {row.relation_id}: {row.status}" for row in self.registry.rows if row.status not in REGISTRY_VERDICTS)
        return tuple(failures)


def load_config(path: Path | str) -> dict[str, Any]:
    try:
        config = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise HarnessError(f"cannot load config {path}: {exc}") from exc
    if not isinstance(config, dict):
        raise HarnessError("config top level must be a mapping")
    _validate_config_keys(config)
    return config


def load_output(path: Path | str, syntax: str | None = None) -> tuple[TaggedRecords, NormalizedRecords]:
    try:
        raw = parse_tagged_output(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise HarnessError(f"cannot parse tagged output {path}: {exc}") from exc
    return raw, normalize_tags(raw, syntax)


def run_checks(
    config: Mapping[str, Any],
    outputs: Mapping[str, Mapping[str, ParsedValue | Unparsed]],
    *,
    registry_directory: Path | str | None = None,
) -> HarnessReport:
    _validate_config_keys(config)
    default_engine = config.get("default_engine", "default")
    if not isinstance(default_engine, str) or default_engine not in outputs:
        raise HarnessError(f"default_engine {default_engine!r} has no loaded output")
    cross_rows = config.get("cross_engine")
    registry_rows = config.get("registry_residual")
    if not isinstance(cross_rows, list):
        cross_rows = []
    if not isinstance(registry_rows, list):
        registry_rows = []
    registry = load_registry(registry_directory or Path(__file__).resolve().parent)
    controls = tuple(check_control_response(values, config, engine=engine) for engine, values in sorted(outputs.items()))
    parity_raw = check_tag_parity(outputs, config)
    assert isinstance(parity_raw, tuple)
    dimensions = tuple(check_dimensions(values, config, engine=engine) for engine, values in sorted(outputs.items()))
    cross = check_cross_engine(cross_rows, outputs, config=config, registry=registry)
    registry_report = check_registry_residuals(registry_rows, outputs, registry, default_engine, config=config)
    diagnostics = tuple(
        IgnoredOutputReport(
            engine,
            tuple(getattr(values, "ignored_lines", ())),
            tuple(getattr(values, "duplicate_tags", ())),
        )
        for engine, values in sorted(outputs.items())
    )
    return HarnessReport(controls, parity_raw, dimensions, cross, registry_report, action_family_reports(cross), diagnostics)


def _status_counts(rows: Iterable[object]) -> str:
    counts = Counter(str(getattr(row, "status")) for row in rows)
    return " ".join(f"{name.lower()}={count}" for name, count in sorted(counts.items())) or "configured=0"


def _format_vector(vector: Dimension) -> str:
    return "[" + ",".join(map(str, vector)) + "]"


def _coverage_line(label: str, coverage: Coverage) -> str:
    return (
        f"{label}: numerator={coverage.numerator} denominator={coverage.denominator} "
        f"formula=\"{coverage.formula}\""
    )


def format_report(report: HarnessReport) -> str:
    lines: list[str] = []
    for diagnostic in report.output_diagnostics:
        lines.append(
            f"IGNORED_LINES[{diagnostic.engine}]: count={len(diagnostic.ignored_lines)} "
            f"duplicate_tags={len(diagnostic.duplicate_tags)}"
        )
        lines.extend(f"  ignored[{index}]={line}" for index, line in enumerate(diagnostic.ignored_lines[:10], 1))
        if diagnostic.duplicate_tags:
            lines.append("  duplicates=" + ",".join(diagnostic.duplicate_tags))
    for control in report.controls:
        lines.append(
            f"CONTROL_RESPONSE[{control.engine}]: compared={control.coverage.numerator} "
            f"responsive={len(control.responsive)} invariant={len(control.invariant)} "
            f"partial={len(control.partial)} gaps={len(control.coverage.gaps)}"
        )
        lines.append(_coverage_line(f"CONTROL_COVERAGE[{control.engine}]", control.coverage))
        lines.append(
            f"CONTROL_MISSING_CELLS[{control.engine}] ({len(control.missing_cells)}): "
            + (", ".join(control.missing_cells) or "none")
        )
        lines.append(f"CONTROL_ROWS[{control.engine}] ({len(control.rows)}):")
        lines.extend(
            f"  {row.cell}:{row.suffix}: {row.status} main={row.main_tag} "
            f"controls=[{','.join(row.control_tags)}] missing=[{','.join(row.missing_tags)}]"
            for row in control.rows
        )
    for parity in report.parity:
        gaps = [row for row in parity.rows if row.missing]
        lines.append(f"TAG_PARITY[{parity.engine}]: cells={len(parity.rows)} gaps={len(gaps)}")
        lines.extend(f"  {row.cell}: missing=[{','.join(row.missing)}]" for row in gaps)
    for dimension in report.dimensions:
        arithmetic = (
            dimension.compared
            + dimension.no_comparison
            + dimension.not_applicable
            + dimension.unwalked
            + dimension.unassessable
            + dimension.unparsed
        )
        lines.append(
            f"DIMENSIONS[{dimension.engine}]: total_tags={dimension.total_tags} "
            f"compared={dimension.compared} no_comparison={dimension.no_comparison} "
            f"not_applicable={dimension.not_applicable} unwalked={dimension.unwalked} "
            f"unassessable={dimension.unassessable} unparsed={dimension.unparsed} "
            f"sum={arithmetic} homogeneous={dimension.homogeneous}"
        )
        lines.append(_coverage_line(f"DIMENSION_COVERAGE[{dimension.engine}]", dimension.coverage))
        lines.append(
            f"DIMENSION_GAPS[{dimension.engine}] ({len(dimension.coverage.gaps)}): "
            + (", ".join(dimension.coverage.gaps) or "none")
        )
        sites = ",".join(f"{site}={count}" for site, count in dimension.proposition_sites) or "none"
        lines.append(f"DIMENSION_PROPOSITIONS[{dimension.engine}]: {sites}")
        lines.append(f"DIMENSION_TABLES[{dimension.engine}] ({len(dimension.tables)}):")
        for table in dimension.tables:
            vectors = ",".join(f"{symbol}={_format_vector(vector)}" for symbol, vector in table.vectors) or "none"
            lines.append(
                f"  {table.package}: {'ASSESSABLE' if table.assessable else 'UNASSESSABLE'} "
                f"vectors=[{vectors}] sources=[{','.join(table.sources)}] failures=[{' | '.join(table.failures)}]"
            )
        lines.append(f"PACKAGE_DIMENSION_DISAGREEMENTS[{dimension.engine}] ({len(dimension.package_disagreements)}):")
        lines.extend(f"  {item}" for item in dimension.package_disagreements)
        bad = [row for row in dimension.statuses if row.status in {"unwalked", "unassessable", "unparsed"}]
        lines.append(f"DIMENSION_OPERATIONAL_ROWS[{dimension.engine}] ({len(bad)}):")
        lines.extend(f"  {row.tag}: {row.status} package={row.package} detail={row.detail}" for row in bad)
        lines.append(f"NON_HOMOGENEOUS[{dimension.engine}] ({len(dimension.non_homogeneous)}):")
        for issue in dimension.non_homogeneous:
            terms = "; ".join(
                f"{term.expression}=>{_format_vector(term.dimension) if term.dimension is not None else 'UNKNOWN'}"
                for term in issue.summands
            )
            lines.append(f"  {issue.tag}: site={issue.site} expression={issue.expression}; {terms}")
    cross = report.cross_engine
    lines.append(f"CROSS_ENGINE: {_status_counts(cross.rows)}")
    lines.append(_coverage_line("CROSS_ENGINE_COVERAGE", cross.coverage))
    lines.append(f"CROSS_ENGINE_GAPS ({len(cross.coverage.gaps)}): " + (", ".join(cross.coverage.gaps) or "none"))
    lines.append(f"DUPLICATE_CROSS_ENGINE_DECLARATIONS ({len(cross.duplicates)}): " + (", ".join(cross.duplicates) or "none"))
    lines.append(
        f"BOOLEAN_OBSERVER: selected_boolean_rows={len(cross.selected_boolean_rows)} "
        f"tree_boolean_rows={len(cross.tree_boolean_rows)}"
    )
    lines.append("  selected=" + (", ".join(cross.selected_boolean_rows) or "none"))
    lines.append("  tree=" + (", ".join(cross.tree_boolean_rows) or "none"))
    lines.append(
        f"NAMING_EFFECT: legacy_before_agree={cross.naming_before_agree} "
        f"declared_after_agree={cross.naming_after_agree} changed_rows={len(cross.naming_changed_rows)}"
    )
    lines.append("  changed=" + (", ".join(cross.naming_changed_rows) or "none"))
    lines.append(f"NAMING_EXCEPTIONS ({len(cross.naming_exceptions)}):")
    lines.extend(f"  {item}" for item in cross.naming_exceptions)
    lines.append(f"CROSS_ENGINE_ROWS ({len(cross.rows)}):")
    for row in cross.rows:
        naming = ",".join(row.naming_applied) or "none"
        identities = ",".join(row.identities_applied) or "none"
        undeclared = ",".join(row.undeclared_spellings) or "none"
        lines.append(
            f"  {row.quantity}: {row.status} family={row.family or 'other'} "
            f"naming=[{naming}] identities=[{identities}] "
            f"undeclared_spellings=[{undeclared}] detail={row.detail}"
        )
        if row.family:
            lines.extend(
                f"    normalized[{engine}]={operand}"
                for engine, operand in row.operands
            )
    for action in report.actions:
        label = action.family.upper() + "_ACTION_COVERAGE"
        prefix = action.family.upper() + "_ACTION"
        lines.append(_coverage_line(label, action.coverage))
        lines.append(f"{prefix}_ROWS ({len(action.statuses)}):")
        lines.extend(f"  {quantity}: {status}" for quantity, status in action.statuses)
        verdicts = [(quantity, status) for quantity, status in action.statuses if status in VERDICTS]
        gaps = [(quantity, status) for quantity, status in action.statuses if status not in VERDICTS]
        lines.append(f"{prefix}_VERDICTS ({len(verdicts)}):")
        lines.extend(f"  {quantity}: {status}" for quantity, status in verdicts)
        lines.append(f"{prefix}_GAPS ({len(gaps)}):")
        if gaps:
            for status in sorted({status for _, status in gaps}):
                names = ",".join(quantity for quantity, observed in gaps if observed == status)
                lines.append(f"  {status}=[{names}]")
        else:
            lines.append("  none")
        representative = next((row for row in cross.rows if row.family == action.family and row.operands), None)
        if representative is not None:
            lines.append(f"  representative={representative.quantity}")
            lines.extend(f"    {engine}={operand}" for engine, operand in representative.operands)
    registry = report.registry
    lines.append(f"REGISTRY_RESIDUAL: {_status_counts(registry.rows)}")
    lines.append(_coverage_line("REGISTRY_RESIDUAL_COVERAGE", registry.coverage))
    lines.append(f"REGISTRY_RESIDUAL_GAPS ({len(registry.coverage.gaps)}): " + (", ".join(registry.coverage.gaps) or "none"))
    lines.append(f"REGISTRY_ROWS ({len(registry.rows)}):")
    lines.extend(
        f"  {row.relation_id}: {row.status} residual={row.residual} "
        f"naming=[{','.join(row.naming_applied) or 'none'}] detail={row.detail}"
        for row in registry.rows
    )
    disagreements = [row for row in cross.rows if row.status == "DISAGREE"]
    registry_disagreements = [row for row in registry.rows if row.status == "NONZERO"]
    lines.append(f"PHYSICS_DISAGREEMENTS ({len(disagreements) + len(registry_disagreements)}):")
    for row in disagreements:
        operands = "; ".join(f"{engine}={operand}" for engine, operand in row.operands)
        lines.append(f"  cross_engine:{row.quantity}: {operands}")
    lines.extend(f"  registry:{row.relation_id}: residual={row.residual}" for row in registry_disagreements)
    lines.append(
        "LIMITATION: comparison and operational coverage only; a physics DISAGREE/NONZERO does not make the process fail."
    )
    return "\n".join(lines)


def _parse_output_arguments(specs: Sequence[str], default_engine: str) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for spec in specs:
        if "=" in spec:
            engine, raw_path = spec.split("=", 1)
        elif len(specs) == 1:
            engine, raw_path = default_engine, spec
        else:
            raise HarnessError("multiple --output arguments require ENGINE=PATH")
        if not engine or not raw_path or engine in result:
            raise HarnessError(f"invalid or duplicate output specification {spec!r}")
        result[engine] = Path(raw_path)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, action="append", metavar="[ENGINE=]PATH")
    parser.add_argument("--registry-dir", type=Path)
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config)
        default_engine = config.get("default_engine", "default")
        if not isinstance(default_engine, str):
            raise HarnessError("default_engine must be a string")
        paths = _parse_output_arguments(args.output, default_engine)
        outputs = {
            engine: load_output(path, syntax=engine if engine in {"wl", "py"} else None)[1]
            for engine, path in paths.items()
        }
        report = run_checks(config, outputs, registry_directory=args.registry_dir)
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
