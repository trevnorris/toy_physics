#!/usr/bin/env python3
"""Compare PY/WL S11b tag streams without interpreting their physics.

The comparator joins non-local tags by emitted object name, parses each
payload, applies only injective mechanical symbol transliteration, and prints
the operands, computed residual, and one of AGREE, DISAGREE, UNDECIDED, or
UNCOMPARED.  A disagreement is a reported result; only operational failures
make the process exit nonzero.
"""

from __future__ import annotations

import argparse
import ast
import io
import keyword
import re
import sys
import tokenize
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational
from sympy.core.symbol import Str
from sympy.parsing.mathematica import parse_mathematica


TAG_LINE = re.compile(r"^(PY|WL)_(S11B_[A-Z][A-Z0-9_]*): (.*)$")
WL_DERIVATIVE = re.compile(
    r"Derivative\[([0-9, ]+)\]\[([A-Za-z][A-Za-z0-9]*)\]"
    r"\[([A-Za-z][A-Za-z0-9]*(?:, [A-Za-z][A-Za-z0-9]*)*)\]"
)
WL_INEQUALITY = re.compile(r"Inequality\[([^\[\]]+)\]")
UNDECIDED_TOKEN = re.compile(r"(?:^|_)UNDECIDED(?:_|$)")


class ComparatorInputError(ValueError):
    """The two inputs are not an unambiguous PY/WL transcript pair."""


class AssociationParseError(ValueError):
    """A Wolfram Association did not have the required keyed shape."""


class DimensionVectorError(ValueError):
    """A DIM payload did not expose exactly one integer L,T,M vector."""


@dataclass(frozen=True)
class Record:
    engine: str
    name: str
    full_tag: str
    raw: str
    path: Path
    line_number: int


@dataclass
class Transcript:
    engine: str
    path: Path
    records: dict[str, Record]
    duplicates: dict[str, list[Record]]
    format_issues: list[str]


@dataclass(frozen=True)
class TextAtom:
    value: str


@dataclass(frozen=True)
class Association:
    entries: tuple[tuple[str, object], ...]

    def as_dict(self) -> dict[str, object]:
        return dict(self.entries)


@dataclass(frozen=True)
class ResidualAssociation:
    entries: tuple[tuple[str, object], ...]


@dataclass(frozen=True)
class Mismatch:
    kind: str
    py_value: object
    wl_value: object
    detail: str | None = None


@dataclass(frozen=True)
class BooleanNotResidualable:
    py_value: object
    wl_value: object


@dataclass(frozen=True)
class ResidualFailure:
    py_value: object
    wl_value: object
    reason: str


@dataclass(frozen=True)
class Comparison:
    name: str
    py_record: Record
    wl_record: Record
    py_object: object | None
    wl_object: object | None
    residual: object | None
    classification: str
    reason: str | None = None
    disagreement_kind: str | None = None
    dimension_difference: tuple[sp.Integer, sp.Integer, sp.Integer] | None = None


def mechanical_lower_camel(name: str) -> str:
    """Apply the S10 snake_case -> lowerCamel spelling rule mechanically."""
    pieces = name.split("_")
    return pieces[0] + "".join(piece[:1].upper() + piece[1:] for piece in pieces[1:])


def is_local_name(name: str) -> bool:
    """Recognize only the exact LOCAL_ infix immediately after S11B_."""
    return name.startswith("S11B_LOCAL_")


def read_transcript(path: Path) -> Transcript:
    if path.suffix != ".out":
        raise ComparatorInputError(f"input is not a .out file: {path}")
    if not path.is_file():
        raise ComparatorInputError(f"input does not exist: {path}")

    parsed: list[Record] = []
    format_issues: list[str] = []
    engines: set[str] = set()
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line:
            continue
        match = TAG_LINE.fullmatch(line)
        if match is None:
            format_issues.append(f"{path}:{line_number}: no S11b tagged-row grammar: {line}")
            continue
        engine, name, raw = match.groups()
        engines.add(engine)
        parsed.append(Record(engine, name, f"{engine}_{name}", raw, path, line_number))

    if len(engines) != 1:
        rendered = ", ".join(sorted(engines)) if engines else "none"
        raise ComparatorInputError(f"{path} has ambiguous engine prefixes: {rendered}")

    engine = next(iter(engines))
    records: dict[str, Record] = {}
    duplicates: dict[str, list[Record]] = {}
    for record in parsed:
        if record.name in records:
            duplicates.setdefault(record.name, [records[record.name]]).append(record)
        else:
            records[record.name] = record
    return Transcript(engine, path, records, duplicates, format_issues)


def load_pair(paths: Sequence[Path]) -> tuple[Transcript, Transcript]:
    if len(paths) != 2:
        raise ComparatorInputError("exactly two .out files are required")
    transcripts = [read_transcript(path) for path in paths]
    by_engine = {transcript.engine: transcript for transcript in transcripts}
    if set(by_engine) != {"PY", "WL"}:
        engines = ", ".join(transcript.engine for transcript in transcripts)
        raise ComparatorInputError(f"inputs must contain one PY and one WL transcript; got {engines}")
    return by_engine["PY"], by_engine["WL"]


def _wl_derivative_replacement(match: re.Match[str]) -> str:
    orders = [int(item.strip()) for item in match.group(1).split(",")]
    variables = [item.strip() for item in match.group(3).split(",")]
    if len(orders) != len(variables):
        return match.group(0)
    nonzero = sorted(
        ((variable, order) for variable, order in zip(variables, orders) if order),
        key=lambda item: item[0],
    )
    dependencies = ", ".join(sorted(set(variables)))
    pairs = ", ".join(f"{{{variable}, {order}}}" for variable, order in nonzero)
    return f"partialDerivative[{match.group(2)}, {{{dependencies}}}, {{{pairs}}}]"


def _wl_inequality_replacement(match: re.Match[str]) -> str:
    tokens = [item.strip() for item in match.group(1).split(",")]
    operators = {"Less": "<", "LessEqual": "<=", "Greater": ">", "GreaterEqual": ">="}
    if len(tokens) < 3 or len(tokens) % 2 == 0:
        return match.group(0)
    relations = []
    for index in range(1, len(tokens), 2):
        operator = operators.get(tokens[index])
        if operator is None:
            return match.group(0)
        relations.append(f"({tokens[index - 1]} {operator} {tokens[index + 1]})")
    return "(" + " && ".join(relations) + ")"


def _encloses(text: str, opening: str, closing: str) -> bool:
    """Return whether one outer enclosure spans the complete string."""
    if not text.startswith(opening) or not text.endswith(closing):
        return False
    stack: list[str] = []
    in_string = False
    escaped = False
    index = 0
    pairs = {"{": "}", "[": "]", "(": ")", "<|": "|>"}
    while index < len(text):
        character = text[index]
        if in_string:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == '"':
                in_string = False
            index += 1
            continue
        if character == '"':
            in_string = True
            index += 1
            continue
        token = text[index : index + 2]
        if token == "<|":
            stack.append("|>")
            index += 2
            continue
        if token == "|>":
            if not stack or stack.pop() != "|>":
                return False
            index += 2
            if not stack and index != len(text):
                return False
            continue
        if character in "{[(":
            stack.append(pairs[character])
        elif character in "}])":
            if not stack or stack.pop() != character:
                return False
            if not stack and index != len(text) - 1:
                return False
        index += 1
    return not stack and not in_string


def _split_top_level(text: str, delimiter: str) -> list[str]:
    """Split a Wolfram container while respecting nested syntax and strings."""
    parts: list[str] = []
    stack: list[str] = []
    in_string = False
    escaped = False
    start = 0
    index = 0
    pairs = {"{": "}", "[": "]", "(": ")"}
    while index < len(text):
        character = text[index]
        if in_string:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == '"':
                in_string = False
            index += 1
            continue
        if character == '"':
            in_string = True
            index += 1
            continue
        pair = text[index : index + 2]
        if pair == "<|":
            stack.append("|>")
            index += 2
            continue
        if pair == "|>":
            if not stack or stack.pop() != "|>":
                raise AssociationParseError("unbalanced |> in Wolfram payload")
            index += 2
            continue
        if character in pairs:
            stack.append(pairs[character])
            index += 1
            continue
        if character in "}])":
            if not stack or stack.pop() != character:
                raise AssociationParseError(f"unbalanced {character} in Wolfram payload")
            index += 1
            continue
        if not stack and text.startswith(delimiter, index):
            parts.append(text[start:index].strip())
            index += len(delimiter)
            start = index
            continue
        index += 1
    if stack or in_string:
        raise AssociationParseError("unbalanced Wolfram payload")
    parts.append(text[start:].strip())
    return parts


def _parse_association_key(raw: str) -> str:
    raw = raw.strip()
    if len(raw) >= 2 and raw[0] == raw[-1] == '"':
        value = ast.literal_eval(raw)
        if not isinstance(value, str):
            raise AssociationParseError(f"Association key is not a string: {raw}")
        return value
    if re.fullmatch(r"[A-Za-z$][A-Za-z0-9_$]*", raw):
        return raw
    raise AssociationParseError(f"Association key is not a string or symbol: {raw}")


def _parse_wolfram_value(raw: str) -> object:
    raw = raw.strip()
    if not raw:
        raise AssociationParseError("empty Wolfram value")
    if _encloses(raw, "<|", "|>"):
        body = raw[2:-2].strip()
        if not body:
            return Association(())
        entries: list[tuple[str, object]] = []
        seen: set[str] = set()
        for entry in _split_top_level(body, ","):
            rule = _split_top_level(entry, "->")
            if len(rule) != 2:
                raise AssociationParseError(f"Association entry is not one Rule: {entry}")
            key = _parse_association_key(rule[0])
            if key in seen:
                raise AssociationParseError(f"duplicate Association key: {key}")
            seen.add(key)
            entries.append((key, _parse_wolfram_value(rule[1])))
        return Association(tuple(entries))
    if _encloses(raw, "{", "}"):
        body = raw[1:-1].strip()
        if not body:
            return ()
        return tuple(_parse_wolfram_value(item) for item in _split_top_level(body, ","))
    if len(raw) >= 2 and raw[0] == raw[-1] == '"':
        value = ast.literal_eval(raw)
        if isinstance(value, str):
            return TextAtom(value)
    prepared = WL_DERIVATIVE.sub(_wl_derivative_replacement, raw)
    prepared = WL_INEQUALITY.sub(_wl_inequality_replacement, prepared)
    return parse_mathematica(prepared)


def parse_wolfram_payload(raw: str) -> object:
    """Parse InputForm, preprocessing every list/Association container."""
    return _parse_wolfram_value(raw)


def _sympy_payload_locals(raw: str) -> dict[str, object]:
    namespace: dict[str, object] = {
        "Str": Str,
        "true": sp.true,
        "false": sp.false,
        "S": sp.S,
    }
    tokens = [
        token
        for token in tokenize.generate_tokens(io.StringIO(raw).readline)
        if token.type
        not in {
            tokenize.ENDMARKER,
            tokenize.INDENT,
            tokenize.DEDENT,
            tokenize.NEWLINE,
            tokenize.NL,
        }
    ]
    for index, token in enumerate(tokens):
        if token.type != tokenize.NAME or keyword.iskeyword(token.string):
            continue
        if token.string in namespace:
            continue
        previous = tokens[index - 1].string if index else ""
        following = tokens[index + 1].string if index + 1 < len(tokens) else ""
        if previous == "." or following in {".", "("}:
            continue
        try:
            resolved = sp.sympify(token.string, evaluate=False)
        except Exception:
            continue
        is_existing_symbol = isinstance(resolved, sp.Symbol)
        is_mathematical_constant = isinstance(resolved, sp.Basic) and not resolved.free_symbols
        if not is_existing_symbol and not is_mathematical_constant:
            namespace[token.string] = sp.Symbol(token.string)
    return namespace


def parse_sympy_payload(raw: str) -> object:
    return sp.sympify(raw, locals=_sympy_payload_locals(raw), evaluate=False)


def _text_value(value: object) -> str | None:
    if isinstance(value, TextAtom):
        return value.value
    if isinstance(value, str):
        return value
    if isinstance(value, Str):
        return str(value)
    if isinstance(value, sp.Symbol):
        return value.name
    if isinstance(value, AppliedUndef) and len(value.args) == 1:
        if value.func.__name__ in {"Str", "_Str"}:
            argument = value.args[0]
            if isinstance(argument, sp.Symbol):
                return argument.name
    return None


def _association_key(value: object) -> str | None:
    text = _text_value(value)
    return text


def _convert_parsed_containers(value: object) -> object:
    if isinstance(value, Association):
        return Association(tuple((key, _convert_parsed_containers(item)) for key, item in value.entries))
    if isinstance(value, Mapping):
        entries: list[tuple[str, object]] = []
        seen: set[str] = set()
        for key, item in value.items():
            converted_key = _association_key(key)
            if converted_key is None:
                raise AssociationParseError(f"SymPy mapping key is not textual: {key!r}")
            if converted_key in seen:
                raise AssociationParseError(f"duplicate SymPy mapping key: {converted_key}")
            seen.add(converted_key)
            entries.append((converted_key, _convert_parsed_containers(item)))
        return Association(tuple(entries))
    if isinstance(value, sp.MatrixBase):
        if value.cols == 1:
            return tuple(_convert_parsed_containers(value[row, 0]) for row in range(value.rows))
        return tuple(
            tuple(_convert_parsed_containers(value[row, column]) for column in range(value.cols))
            for row in range(value.rows)
        )
    if isinstance(value, (tuple, list, sp.Tuple)):
        converted = tuple(_convert_parsed_containers(item) for item in value)
        if converted and all(
            isinstance(item, tuple) and len(item) == 2 and _association_key(item[0]) is not None
            for item in converted
        ):
            entries = []
            seen: set[str] = set()
            for key_object, item in converted:
                key = _association_key(key_object)
                assert key is not None
                if key in seen:
                    raise AssociationParseError(f"duplicate tuple-Association key: {key}")
                seen.add(key)
                entries.append((key, item))
            return Association(tuple(entries))
        return converted
    text = _text_value(value)
    if text is not None and not isinstance(value, sp.Symbol):
        return TextAtom(text)
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return sp.Integer(value)
    if isinstance(value, float):
        return sp.Float(value)
    return value


def _partial_derivative(expression: sp.Derivative) -> object:
    applied = expression.expr
    if not isinstance(applied, AppliedUndef):
        return expression
    orders: dict[sp.Symbol, int] = {}
    for variable, count in expression.variable_count:
        if not isinstance(variable, sp.Symbol):
            return expression
        orders[variable] = orders.get(variable, 0) + count
    pairs = sp.Tuple(
        *(sp.Tuple(variable, count) for variable, count in sorted(orders.items(), key=lambda item: item[0].name))
    )
    dependencies = sp.Tuple(*sorted(set(applied.args), key=sp.sstr))
    return sp.Function("partialDerivative")(sp.Symbol(applied.func.__name__), dependencies, pairs)


def _iter_basic_values(value: object) -> Iterable[sp.Basic]:
    if isinstance(value, Association):
        for _, item in value.entries:
            yield from _iter_basic_values(item)
    elif isinstance(value, tuple):
        for item in value:
            yield from _iter_basic_values(item)
    elif isinstance(value, sp.Basic):
        yield value


def transliteration_collisions(value: object) -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Return every non-injective target in one engine's source vocabulary."""
    by_target: dict[str, set[str]] = {}
    for basic in _iter_basic_values(value):
        for symbol in basic.atoms(sp.Symbol):
            by_target.setdefault(mechanical_lower_camel(symbol.name), set()).add(symbol.name)
    return tuple(
        (target, tuple(sorted(sources)))
        for target, sources in sorted(by_target.items())
        if len(sources) > 1
    )


def _transliterate_basic(value: sp.Basic) -> sp.Basic:
    value = value.replace(lambda item: isinstance(item, sp.Derivative), _partial_derivative)
    replacements = {
        symbol: sp.Symbol(
            mechanical_lower_camel(symbol.name),
            **getattr(symbol, "_assumptions_orig", {}),
        )
        for symbol in value.atoms(sp.Symbol)
        if "_" in symbol.name
    }
    if replacements:
        value = value.xreplace(replacements)
    return value.replace(
        lambda item: isinstance(item, AppliedUndef) and "_" in item.func.__name__,
        lambda item: sp.Function(mechanical_lower_camel(item.func.__name__))(*item.args),
    )


def transliterate(value: object) -> object:
    """Apply only container canonicalization and mechanical symbol spelling."""
    value = _convert_parsed_containers(value)
    if isinstance(value, Association):
        return Association(tuple((key, transliterate(item)) for key, item in value.entries))
    if isinstance(value, tuple):
        return tuple(transliterate(item) for item in value)
    if isinstance(value, sp.Basic):
        return _transliterate_basic(value)
    return value


def _is_native_boolean(value: object) -> bool:
    return isinstance(value, bool) or value is sp.true or value is sp.false


def _canonical_basic_same(left: sp.Basic, right: sp.Basic) -> bool:
    """Compare non-residualable SymPy structures without object equality."""
    return type(left) is type(right) and sp.srepr(left) == sp.srepr(right)


def residual(py_value: object, wl_value: object) -> object:
    """Compute a recursive exact residual without using operand equality."""
    if _is_native_boolean(py_value) or _is_native_boolean(wl_value):
        return BooleanNotResidualable(py_value, wl_value)

    if isinstance(py_value, Association) or isinstance(wl_value, Association):
        if not isinstance(py_value, Association) or not isinstance(wl_value, Association):
            return Mismatch(
                "STRUCTURE_DISAGREE",
                type(py_value).__name__,
                type(wl_value).__name__,
                "tuple/Association or scalar/Association",
            )
        py_items = py_value.as_dict()
        wl_items = wl_value.as_dict()
        py_keys = set(py_items)
        wl_keys = set(wl_items)
        if py_keys != wl_keys:
            missing_from_py = sorted(wl_keys - py_keys)
            missing_from_wl = sorted(py_keys - wl_keys)
            return Mismatch(
                "KEY_DISAGREE",
                tuple(sorted(py_keys)),
                tuple(sorted(wl_keys)),
                f"MISSING_FROM_PY={missing_from_py} MISSING_FROM_WL={missing_from_wl}",
            )
        return ResidualAssociation(
            tuple((key, residual(py_items[key], wl_items[key])) for key in sorted(py_keys))
        )

    if isinstance(py_value, tuple) or isinstance(wl_value, tuple):
        if not isinstance(py_value, tuple) or not isinstance(wl_value, tuple):
            return Mismatch(
                "STRUCTURE_DISAGREE",
                type(py_value).__name__,
                type(wl_value).__name__,
                "tuple/non-tuple",
            )
        if len(py_value) != len(wl_value):
            return Mismatch("STRUCTURE_DISAGREE", len(py_value), len(wl_value), "unequal tuple length")
        return tuple(residual(left, right) for left, right in zip(py_value, wl_value))

    if isinstance(py_value, Relational) or isinstance(wl_value, Relational):
        if not isinstance(py_value, Relational) or not isinstance(wl_value, Relational):
            return Mismatch(
                "STRUCTURE_DISAGREE",
                type(py_value).__name__,
                type(wl_value).__name__,
                "relational/non-relational",
            )
        if type(py_value) is not type(wl_value):
            return Mismatch(
                "STRUCTURE_DISAGREE",
                type(py_value).__name__,
                type(wl_value).__name__,
                "relational operator",
            )
        return (residual(py_value.lhs, wl_value.lhs), residual(py_value.rhs, wl_value.rhs))

    if isinstance(py_value, TextAtom) or isinstance(wl_value, TextAtom):
        if isinstance(py_value, TextAtom) and isinstance(wl_value, TextAtom):
            if py_value.value == wl_value.value:
                return sp.S.Zero
            return Mismatch("TOKEN_DISAGREE", py_value, wl_value)
        return Mismatch("STRUCTURE_DISAGREE", py_value, wl_value, "text/non-text")

    if isinstance(py_value, sp.Expr) and isinstance(wl_value, sp.Expr):
        try:
            return sp.factor(sp.cancel(sp.together(py_value - wl_value)))
        except Exception as error:
            return ResidualFailure(py_value, wl_value, f"{type(error).__name__}: {error}")

    if isinstance(py_value, sp.Basic) and isinstance(wl_value, sp.Basic):
        if _canonical_basic_same(py_value, wl_value):
            return sp.S.Zero
        return Mismatch("STRUCTURAL_ATOM_DISAGREE", py_value, wl_value)

    return ResidualFailure(
        py_value,
        wl_value,
        f"unsupported residual leaf types {type(py_value).__name__}/{type(wl_value).__name__}",
    )


def _iter_residual_leaves(value: object, path: str = "$") -> Iterable[tuple[str, object]]:
    if isinstance(value, ResidualAssociation):
        for key, item in value.entries:
            yield from _iter_residual_leaves(item, f"{path}.{key}")
    elif isinstance(value, tuple):
        for index, item in enumerate(value):
            yield from _iter_residual_leaves(item, f"{path}[{index}]")
    else:
        yield path, value


def _is_structural_zero(value: object) -> bool:
    return value is sp.S.Zero


def classify_residual(value: object) -> tuple[str, str | None]:
    """Use DISAGREE > UNCOMPARED > AGREE precedence across sibling leaves."""
    leaves = list(_iter_residual_leaves(value))
    disagreement_kinds: list[str] = []
    operational = False
    for _, leaf in leaves:
        if isinstance(leaf, BooleanNotResidualable):
            operational = True
        elif isinstance(leaf, ResidualFailure):
            operational = True
        elif isinstance(leaf, Mismatch):
            disagreement_kinds.append(leaf.kind)
        elif not _is_structural_zero(leaf):
            disagreement_kinds.append("NONZERO_RESIDUAL")
    if disagreement_kinds:
        if "KEY_DISAGREE" in disagreement_kinds:
            kind = "KEY"
        elif "STRUCTURE_DISAGREE" in disagreement_kinds:
            kind = "STRUCTURE"
        else:
            kind = "CONTENT"
        return "DISAGREE", kind
    if operational:
        return "UNCOMPARED", None
    return "AGREE", None


def _is_uncertain_token(token: str) -> bool:
    upper = token.upper()
    return bool(UNDECIDED_TOKEN.search(upper)) or upper.startswith("INCOMPLETE_") or upper in {
        "NOT_ESTABLISHED",
        "NOT_DEFINED_ON_COMPONENT",
    }


def authoritative_undecided_token(name: str, value: object) -> str | None:
    """Find only an explicit STATUS_TOKEN or coverage-status payload."""
    if isinstance(value, Association):
        items = value.as_dict()
        status = items.get("STATUS_TOKEN")
        status_text = _text_value(status)
        if status_text is not None and _is_uncertain_token(status_text):
            return f"STATUS_TOKEN={status_text}"
        for key, item in value.entries:
            if "COVERAGE" not in key.upper():
                continue
            coverage_text = _text_value(item)
            if coverage_text is not None and _is_uncertain_token(coverage_text):
                return f"{key}={coverage_text}"
        return None
    token = _text_value(value)
    if token is None or not _is_uncertain_token(token):
        return None
    if token.upper() == "UNDECIDED" or "COVERAGE" in name or "STATUS" in name:
        return token
    return None


def _strict_integer(value: object) -> sp.Integer | None:
    if _is_native_boolean(value):
        return None
    if isinstance(value, sp.Integer):
        return value
    if isinstance(value, int):
        return sp.Integer(value)
    return None


def _vector_candidate(value: object) -> tuple[sp.Integer, sp.Integer, sp.Integer] | None:
    if not isinstance(value, tuple):
        return None
    vector = value
    if len(vector) == 1 and isinstance(vector[0], tuple):
        vector = vector[0]
    elif len(vector) == 3 and all(isinstance(item, tuple) and len(item) == 1 for item in vector):
        vector = tuple(item[0] for item in vector)
    if len(vector) != 3:
        return None
    integers = tuple(_strict_integer(item) for item in vector)
    if any(item is None for item in integers):
        return None
    return integers  # type: ignore[return-value]


def extract_dimension_vector(value: object) -> tuple[sp.Integer, sp.Integer, sp.Integer]:
    direct = _vector_candidate(value)
    if direct is not None:
        return direct
    if isinstance(value, Association):
        candidates = [candidate for _, item in value.entries if (candidate := _vector_candidate(item)) is not None]
        if len(candidates) == 1:
            return candidates[0]
        raise DimensionVectorError(
            f"Association exposes {len(candidates)} direct integer-vector fields; expected exactly one"
        )
    raise DimensionVectorError(
        f"payload type {type(value).__name__} is neither an integer vector nor an Association with one"
    )


def _is_dimension_name(name: str) -> bool:
    return name.startswith("S11B_DIM_") and not name.startswith("S11B_DIM_ROUTE_KIND")


def _collision_reason(engine: str, collisions: tuple[tuple[str, tuple[str, ...]], ...]) -> str:
    rendered = "; ".join(f"TARGET={target} SOURCES={sources}" for target, sources in collisions)
    return f"TRANSLITERATION_COLLISION ENGINE={engine} {rendered}"


def compare_records(name: str, py_record: Record, wl_record: Record) -> Comparison:
    try:
        parsed_py = _convert_parsed_containers(parse_sympy_payload(py_record.raw))
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, None, None, None, "UNCOMPARED",
            f"PY_PARSE_{type(error).__name__}: {error}",
        )
    try:
        parsed_wl = _convert_parsed_containers(parse_wolfram_payload(wl_record.raw))
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, parsed_py, None, None, "UNCOMPARED",
            f"WL_PARSE_{type(error).__name__}: {error}",
        )

    py_collisions = transliteration_collisions(parsed_py)
    if py_collisions:
        return Comparison(
            name, py_record, wl_record, parsed_py, parsed_wl, None, "UNCOMPARED",
            _collision_reason("PY", py_collisions),
        )
    wl_collisions = transliteration_collisions(parsed_wl)
    if wl_collisions:
        return Comparison(
            name, py_record, wl_record, parsed_py, parsed_wl, None, "UNCOMPARED",
            _collision_reason("WL", wl_collisions),
        )

    try:
        py_object = transliterate(parsed_py)
        wl_object = transliterate(parsed_wl)
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, parsed_py, parsed_wl, None, "UNCOMPARED",
            f"TRANSLITERATE_{type(error).__name__}: {error}",
        )

    py_undecided = authoritative_undecided_token(name, py_object)
    wl_undecided = authoritative_undecided_token(name, wl_object)
    if py_undecided is not None or wl_undecided is not None:
        sources = []
        if py_undecided is not None:
            sources.append(f"PY={py_undecided}")
        if wl_undecided is not None:
            sources.append(f"WL={wl_undecided}")
        return Comparison(
            name, py_record, wl_record, py_object, wl_object, None, "UNDECIDED",
            "ENGINE_STATUS " + " ".join(sources),
        )

    if _is_dimension_name(name):
        try:
            py_vector = extract_dimension_vector(py_object)
            wl_vector = extract_dimension_vector(wl_object)
            difference = tuple(left - right for left, right in zip(py_vector, wl_vector))
            result = residual(py_vector, wl_vector)
            classification, _ = classify_residual(result)
            disagreement_kind = "DIM" if classification == "DISAGREE" else None
            return Comparison(
                name,
                py_record,
                wl_record,
                py_vector,
                wl_vector,
                result,
                classification,
                None,
                disagreement_kind,
                difference,
            )
        except Exception as error:
            return Comparison(
                name, py_record, wl_record, py_object, wl_object, None, "UNCOMPARED",
                f"DIM_VECTOR_{type(error).__name__}: {error}",
            )

    try:
        difference = residual(py_object, wl_object)
        classification, disagreement_kind = classify_residual(difference)
        reason = None
        if classification == "UNCOMPARED":
            reasons = []
            for path, leaf in _iter_residual_leaves(difference):
                if isinstance(leaf, BooleanNotResidualable):
                    reasons.append(f"{path}: BOOLEAN_NOT_RESIDUALABLE")
                elif isinstance(leaf, ResidualFailure):
                    reasons.append(f"{path}: RESIDUAL_FAILURE {leaf.reason}")
            reason = "; ".join(reasons)
        return Comparison(
            name,
            py_record,
            wl_record,
            py_object,
            wl_object,
            difference,
            classification,
            reason,
            disagreement_kind,
        )
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, py_object, wl_object, None, "UNCOMPARED",
            f"RESIDUAL_{type(error).__name__}: {error}",
        )


def render_value(value: object) -> str:
    if isinstance(value, TextAtom):
        return repr(value.value)
    if isinstance(value, Association):
        return "<|" + ", ".join(
            f"{key!r}: {render_value(item)}" for key, item in value.entries
        ) + "|>"
    if isinstance(value, ResidualAssociation):
        return "<|" + ", ".join(
            f"{key!r}: {render_value(item)}" for key, item in value.entries
        ) + "|>"
    if isinstance(value, Mismatch):
        detail = f", DETAIL={value.detail}" if value.detail is not None else ""
        return (
            f"{value.kind}(PY={render_value(value.py_value)}, "
            f"WL={render_value(value.wl_value)}{detail})"
        )
    if isinstance(value, BooleanNotResidualable):
        return (
            "BOOLEAN_NOT_RESIDUALABLE("
            f"PY={render_value(value.py_value)}, WL={render_value(value.wl_value)})"
        )
    if isinstance(value, ResidualFailure):
        return (
            "RESIDUAL_FAILURE("
            f"PY={render_value(value.py_value)}, WL={render_value(value.wl_value)}, "
            f"REASON={value.reason!r})"
        )
    if isinstance(value, tuple):
        contents = ", ".join(render_value(item) for item in value)
        if len(value) == 1:
            contents += ","
        return f"({contents})"
    if isinstance(value, sp.Basic):
        return sp.sstr(value)
    return repr(value)


def render_comparison(comparison: Comparison) -> list[str]:
    lines = [
        f"COMPARE_NAME: {comparison.name}",
        f"PY_OPERAND: {comparison.py_record.raw}",
        f"WL_OPERAND: {comparison.wl_record.raw}",
    ]
    if comparison.py_object is not None:
        lines.append(f"PY_PARSED_OBJECT: {render_value(comparison.py_object)}")
    if comparison.wl_object is not None:
        lines.append(f"WL_PARSED_OBJECT: {render_value(comparison.wl_object)}")
    if comparison.residual is None:
        lines.append("RESIDUAL: NOT_COMPUTED")
    else:
        lines.append(f"RESIDUAL: {render_value(comparison.residual)}")
        for path, leaf in _iter_residual_leaves(comparison.residual):
            if isinstance(leaf, BooleanNotResidualable):
                lines.append(
                    f"LEAF_OUTCOME: PATH={path} KIND=BOOLEAN_NOT_RESIDUALABLE "
                    f"PY={render_value(leaf.py_value)} WL={render_value(leaf.wl_value)}"
                )
            elif isinstance(leaf, ResidualFailure):
                lines.append(
                    f"LEAF_OUTCOME: PATH={path} KIND=RESIDUAL_FAILURE REASON={leaf.reason}"
                )
            elif isinstance(leaf, Mismatch):
                lines.append(f"LEAF_OUTCOME: PATH={path} KIND={leaf.kind}")
    if comparison.dimension_difference is not None:
        lines.append(f"DIM_DIFFERENCE_VECTOR: {render_value(comparison.dimension_difference)}")
    if comparison.reason is not None:
        label = "UNDECIDED_REASON" if comparison.classification == "UNDECIDED" else "UNCOMPARED_REASON"
        lines.append(f"{label}: {comparison.reason}")
    if comparison.disagreement_kind is not None:
        lines.append(f"DISAGREEMENT_KIND: {comparison.disagreement_kind}")
    lines.append(f"CLASSIFICATION: {comparison.classification}")
    return lines


def _render_duplicate_name(py: Transcript, wl: Transcript, name: str) -> list[str]:
    lines = [f"COMPARE_NAME: {name}"]
    for transcript in (py, wl):
        rows = transcript.duplicates.get(name)
        if rows is not None:
            for row in rows:
                lines.append(f"{transcript.engine}_OPERAND_AT_LINE_{row.line_number}: {row.raw}")
        elif name in transcript.records:
            row = transcript.records[name]
            lines.append(f"{transcript.engine}_OPERAND_AT_LINE_{row.line_number}: {row.raw}")
    lines.extend(
        (
            "RESIDUAL: NOT_COMPUTED",
            "UNCOMPARED_REASON: duplicate emitted non-local name",
            "CLASSIFICATION: UNCOMPARED",
        )
    )
    return lines


def _local_names(transcript: Transcript) -> list[str]:
    return sorted(name for name in transcript.records if is_local_name(name))


def run(py: Transcript, wl: Transcript) -> int:
    print(f"PY_INPUT: {py.path}")
    print(f"WL_INPUT: {wl.path}")
    print("JOIN_RULE: pair PY_S11B_<QUANTITY> with WL_S11B_<QUANTITY> by exact non-local object name")
    print("SYMBOL_RULE: injective snake_case to lowerCamel transliteration; collisions are UNCOMPARED")
    print("STRUCTURE_RULE: residual only equal-length tuples or Associations with equal key sets")
    print("BOOLEAN_RULE: a native boolean leaf is BOOLEAN_NOT_RESIDUALABLE and never scores AGREE")
    print("DIM_RULE: extract and subtract exact L,T,M integer vectors; no normalization")

    py_locals = _local_names(py)
    wl_locals = _local_names(wl)
    print("CATEGORY: PY_LOCAL_NAMES")
    print(f"PY_LOCAL_COUNT: {len(py_locals)}")
    for name in py_locals:
        print(f"PY_LOCAL_NAME: {name}")
    print("CATEGORY: WL_LOCAL_NAMES")
    print(f"WL_LOCAL_COUNT: {len(wl_locals)}")
    for name in wl_locals:
        print(f"WL_LOCAL_NAME: {name}")

    py_nonlocal = {name for name in py.records if not is_local_name(name)}
    wl_nonlocal = {name for name in wl.records if not is_local_name(name)}
    duplicate_names = {
        name
        for name in set(py.duplicates) | set(wl.duplicates)
        if not is_local_name(name)
    }
    shared_names = py_nonlocal & wl_nonlocal
    shared_duplicate_names = shared_names & duplicate_names
    comparable_names = sorted(shared_names - shared_duplicate_names)
    comparisons = [compare_records(name, py.records[name], wl.records[name]) for name in comparable_names]

    print("CATEGORY: SHARED_OBJECTS")
    for comparison in comparisons:
        for line in render_comparison(comparison):
            print(line)

    print("CATEGORY: DUPLICATE_NONLOCAL_NAMES")
    for name in sorted(duplicate_names):
        for line in _render_duplicate_name(py, wl, name):
            print(line)

    print("CATEGORY: UNPAIRED_NONLOCAL_COVERAGE")
    py_only = sorted(py_nonlocal - wl_nonlocal)
    wl_only = sorted(wl_nonlocal - py_nonlocal)
    print(f"PY_ONLY_NONLOCAL_COUNT: {len(py_only)}")
    for name in py_only:
        print(f"UNPAIRED_ENGINE: PY UNPAIRED_NAME: {name}")
    print(f"WL_ONLY_NONLOCAL_COUNT: {len(wl_only)}")
    for name in wl_only:
        print(f"UNPAIRED_ENGINE: WL UNPAIRED_NAME: {name}")

    print("CATEGORY: TRANSCRIPT_FORMAT_ISSUES")
    for issue in py.format_issues + wl.format_issues:
        print(f"FORMAT_ISSUE: {issue}")
        print("RESIDUAL: NOT_COMPUTED")
        print("UNCOMPARED_REASON: nonempty transcript line has no S11b emitted-name row grammar")
        print("CLASSIFICATION: UNCOMPARED")

    classifications = {label: 0 for label in ("AGREE", "DISAGREE", "UNDECIDED", "UNCOMPARED")}
    for comparison in comparisons:
        classifications[comparison.classification] += 1
    classifications["UNCOMPARED"] += len(shared_duplicate_names)
    format_issue_count = len(py.format_issues) + len(wl.format_issues)
    duplicate_row_count = sum(
        len(rows) - 1
        for name, rows in py.duplicates.items()
        if not is_local_name(name)
    ) + sum(
        len(rows) - 1
        for name, rows in wl.duplicates.items()
        if not is_local_name(name)
    )

    print("CATEGORY: SUMMARY")
    print(f"SHARED_NONLOCAL_NAME_COUNT: {len(shared_names)}")
    print(f"AGREE_COUNT: {classifications['AGREE']}")
    print(f"DISAGREE_COUNT: {classifications['DISAGREE']}")
    print(f"UNDECIDED_COUNT: {classifications['UNDECIDED']}")
    print(f"UNCOMPARED_SHARED_COUNT: {classifications['UNCOMPARED']}")
    print(f"UNPAIRED_NONLOCAL_COUNT: {len(py_only) + len(wl_only)}")
    print(f"LOCAL_NAME_COUNT: {len(py_locals) + len(wl_locals)}")
    print(f"DUPLICATE_NONLOCAL_ROW_COUNT: {duplicate_row_count}")
    print(f"FORMAT_ISSUE_COUNT: {format_issue_count}")

    operational_failure = bool(
        classifications["UNCOMPARED"] or duplicate_row_count or format_issue_count
    )
    print(f"FINAL_OPERATIONAL_STATUS: {'FAIL' if operational_failure else 'PASS'}")
    return int(operational_failure)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Join one synthetic-or-real PY and WL S11b .out transcript by emitted object name."
    )
    parser.add_argument("out_files", nargs=2, type=Path, metavar="OUT")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        py, wl = load_pair(args.out_files)
    except ComparatorInputError as error:
        print(f"INPUT_ERROR: {error}")
        print("FINAL_OPERATIONAL_STATUS: FAIL")
        return 2
    return run(py, wl)


if __name__ == "__main__":
    raise SystemExit(main())
