#!/usr/bin/env python3
"""Join two CAS transcripts on emitted object name and compare exact payloads."""

from __future__ import annotations

import argparse
import io
import keyword
import re
import sys
import tokenize
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational
from sympy.parsing.mathematica import parse_mathematica


TAG_LINE = re.compile(r"^(PY|WL)_([A-Z][A-Z0-9_]*): (.*)$")
WL_DERIVATIVE = re.compile(
    r"Derivative\[([0-9, ]+)\]\[([A-Za-z][A-Za-z0-9]*)\]"
    r"\[([A-Za-z][A-Za-z0-9]*(?:, [A-Za-z][A-Za-z0-9]*)*)\]"
)
WL_INEQUALITY = re.compile(r"Inequality\[([^\[\]]+)\]")
NULLSPACE_BASIS_SUFFIX = "N6_NULLSPACE_BASIS"


class ComparatorInputError(ValueError):
    """The two arguments are not an unambiguous PY/WL transcript pair."""


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
class Mismatch:
    kind: str
    py_value: object
    wl_value: object


@dataclass(frozen=True)
class Comparison:
    name: str
    py_record: Record
    wl_record: Record
    py_object: object | None
    wl_object: object | None
    residual: object | None
    passes: bool
    reason: str | None
    symbol_pairs: tuple[tuple[str, str, str], ...]
    disagreement_kind: str | None = None


def mechanical_lower_camel(name: str) -> str:
    """Apply snake_case -> lowerCamel without an exception table."""
    pieces = name.split("_")
    return pieces[0] + "".join(piece[:1].upper() + piece[1:] for piece in pieces[1:])


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
            format_issues.append(f"{path}:{line_number}: no tagged-row grammar: {line}")
            continue
        engine, name, raw = match.groups()
        engines.add(engine)
        parsed.append(
            Record(
                engine=engine,
                name=name,
                full_tag=f"{engine}_{name}",
                raw=raw,
                path=path,
                line_number=line_number,
            )
        )

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


def parse_wolfram_payload(raw: str) -> object:
    if raw == "{}":
        return ()
    prepared = raw.replace("<|", "Association[").replace("|>", "]")
    prepared = WL_DERIVATIVE.sub(_wl_derivative_replacement, prepared)
    prepared = WL_INEQUALITY.sub(_wl_inequality_replacement, prepared)
    return parse_mathematica(prepared)


def _sympy_payload_locals(raw: str) -> dict[str, sp.Symbol]:
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
    namespace: dict[str, sp.Symbol] = {}
    for index, token in enumerate(tokens):
        if token.type != tokenize.NAME or keyword.iskeyword(token.string):
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
        is_mathematical_constant = (
            isinstance(resolved, sp.Basic) and not resolved.free_symbols
        )
        if not is_existing_symbol and not is_mathematical_constant:
            namespace[token.string] = sp.Symbol(token.string)
    return namespace


def parse_sympy_payload(raw: str) -> object:
    return sp.sympify(raw, locals=_sympy_payload_locals(raw), evaluate=False)


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


def _renamed_symbol(symbol: sp.Symbol) -> sp.Symbol:
    name = mechanical_lower_camel(symbol.name)
    assumptions = getattr(symbol, "_assumptions_orig", {})
    return sp.Symbol(name, **assumptions)


def normalize(value: object) -> object:
    """Normalize only mechanical spelling, sequence, power, and parsed derivative syntax."""
    if isinstance(value, sp.MatrixBase):
        if value.cols == 1:
            return tuple(normalize(value[row, 0]) for row in range(value.rows))
        return tuple(
            tuple(normalize(value[row, column]) for column in range(value.cols))
            for row in range(value.rows)
        )
    if isinstance(value, (tuple, list, sp.Tuple)):
        return tuple(normalize(item) for item in value)
    if isinstance(value, sp.Basic):
        value = value.replace(lambda item: isinstance(item, sp.Derivative), _partial_derivative)
        replacements = {
            symbol: _renamed_symbol(symbol)
            for symbol in value.atoms(sp.Symbol)
            if "_" in symbol.name
        }
        if replacements:
            value = value.xreplace(replacements)
        value = value.replace(
            lambda item: isinstance(item, AppliedUndef) and "_" in item.func.__name__,
            lambda item: sp.Function(mechanical_lower_camel(item.func.__name__))(*item.args),
        )
    return value


def residual(py_value: object, wl_value: object) -> object:
    if isinstance(py_value, tuple) or isinstance(wl_value, tuple):
        if not isinstance(py_value, tuple) or not isinstance(wl_value, tuple):
            return Mismatch("TYPE", type(py_value).__name__, type(wl_value).__name__)
        if len(py_value) != len(wl_value):
            return Mismatch("LENGTH", len(py_value), len(wl_value))
        return tuple(residual(left, right) for left, right in zip(py_value, wl_value))

    if isinstance(py_value, Relational) or isinstance(wl_value, Relational):
        if not isinstance(py_value, Relational) or not isinstance(wl_value, Relational):
            return Mismatch("RELATIONAL_TYPE", type(py_value).__name__, type(wl_value).__name__)
        if type(py_value) is not type(wl_value):
            return Mismatch("RELATIONAL_OPERATOR", type(py_value).__name__, type(wl_value).__name__)
        return (
            residual(py_value.lhs, wl_value.lhs),
            residual(py_value.rhs, wl_value.rhs),
        )

    if py_value == wl_value:
        return sp.Integer(0)
    if isinstance(py_value, sp.Expr) and isinstance(wl_value, sp.Expr):
        try:
            return sp.factor(sp.cancel(sp.together(py_value - wl_value)))
        except Exception as error:
            return Mismatch(type(error).__name__, py_value, wl_value)
    return Mismatch("STRUCTURE", py_value, wl_value)


def residual_is_zero(value: object) -> bool:
    if isinstance(value, tuple):
        return all(residual_is_zero(item) for item in value)
    if isinstance(value, Mismatch):
        return False
    return value == 0


def _nested_free_symbols(value: object) -> set[sp.Symbol]:
    if isinstance(value, tuple):
        return set().union(*(_nested_free_symbols(item) for item in value)) if value else set()
    if isinstance(value, sp.Basic):
        return set(value.free_symbols)
    return set()


def _contains_mismatch(value: object) -> bool:
    if isinstance(value, tuple):
        return any(_contains_mismatch(item) for item in value)
    return isinstance(value, Mismatch)


def _symbol_spellings(value: object, *, mechanical: bool) -> dict[str, str]:
    if isinstance(value, sp.MatrixBase):
        values: Iterable[object] = tuple(value)
    elif isinstance(value, (tuple, list, sp.Tuple)):
        values = value
    elif isinstance(value, dict):
        values = tuple(value.keys()) + tuple(value.values())
    elif isinstance(value, sp.Basic):
        return {
            mechanical_lower_camel(symbol.name) if mechanical else symbol.name: symbol.name
            for symbol in value.free_symbols
        }
    else:
        return {}
    result: dict[str, str] = {}
    for item in values:
        result.update(_symbol_spellings(item, mechanical=mechanical))
    return result


def _paired_leaves(py_value: object, wl_value: object) -> Iterable[tuple[object, object]]:
    if isinstance(py_value, tuple) and isinstance(wl_value, tuple) and len(py_value) == len(wl_value):
        for left, right in zip(py_value, wl_value):
            yield from _paired_leaves(left, right)
    else:
        yield py_value, wl_value


def _nested_xreplace(value: object, replacements: Mapping[sp.Symbol, sp.Symbol]) -> object:
    if isinstance(value, tuple):
        return tuple(_nested_xreplace(item, replacements) for item in value)
    if isinstance(value, sp.Basic):
        return value.xreplace(replacements)
    return value


def infer_nonmechanical_symbol_pairs(
    py_value: object,
    wl_value: object,
    py_spellings: Mapping[str, str],
) -> tuple[tuple[str, str, str], ...]:
    py_symbols = _nested_free_symbols(py_value)
    wl_symbols = _nested_free_symbols(wl_value)
    py_only = py_symbols - wl_symbols
    wl_only = wl_symbols - py_symbols
    if not py_only or not wl_only:
        return ()

    wildcards = {symbol: sp.Wild(f"candidate_{index}") for index, symbol in enumerate(sorted(py_only, key=str))}
    wildcard_sources = {wildcard: symbol for symbol, wildcard in wildcards.items()}
    inferred: dict[sp.Symbol, sp.Symbol] = {}
    for py_leaf, wl_leaf in _paired_leaves(py_value, wl_value):
        if not isinstance(py_leaf, sp.Basic) or not isinstance(wl_leaf, sp.Basic):
            continue
        match = wl_leaf.match(py_leaf.xreplace(wildcards))
        if match is None:
            continue
        for wildcard, replacement in match.items():
            source = wildcard_sources.get(wildcard)
            if source is None or not isinstance(replacement, sp.Symbol) or replacement not in wl_only:
                continue
            prior = inferred.get(source)
            if prior is not None and prior != replacement:
                return ()
            inferred[source] = replacement

    if not inferred and len(py_only) == 1 and len(wl_only) == 1:
        source = next(iter(py_only))
        target = next(iter(wl_only))
        if residual_is_zero(residual(_nested_xreplace(py_value, {source: target}), wl_value)):
            inferred[source] = target

    if not inferred or len(set(inferred.values())) != len(inferred):
        return ()
    renamed = _nested_xreplace(py_value, inferred)
    if not residual_is_zero(residual(renamed, wl_value)):
        return ()
    return tuple(
        sorted(
            (
                py_spellings.get(source.name, source.name),
                source.name,
                target.name,
            )
            for source, target in inferred.items()
        )
    )


def _equation_system_as_expressions(value: object) -> object:
    if isinstance(value, tuple):
        return tuple(_equation_system_as_expressions(item) for item in value)
    if isinstance(value, sp.Equality):
        return value.lhs - value.rhs
    return value


def _has_equality(value: object) -> bool:
    if isinstance(value, tuple):
        return any(_has_equality(item) for item in value)
    return isinstance(value, sp.Equality)


def disagreement_kind(
    py_value: object,
    wl_value: object,
    symbol_pairs: tuple[tuple[str, str, str], ...],
    refuted_symbol_pairs: frozenset[tuple[str, str, str]] = frozenset(),
) -> str:
    refuted_normalized_pairs = {
        (py_mechanical, wl_name)
        for _, py_mechanical, wl_name in refuted_symbol_pairs
    }

    def contains_refuted_pair(pairs: tuple[tuple[str, str, str], ...]) -> bool:
        return any(
            (py_mechanical, wl_name) in refuted_normalized_pairs
            for _, py_mechanical, wl_name in pairs
        )

    if _has_equality(py_value) != _has_equality(wl_value):
        projected_py = _equation_system_as_expressions(py_value)
        projected_wl = _equation_system_as_expressions(wl_value)
        if residual_is_zero(residual(projected_py, projected_wl)):
            return "REPRESENTATIONAL"
        projected_pairs = infer_nonmechanical_symbol_pairs(projected_py, projected_wl, {})
        if projected_pairs and not contains_refuted_pair(projected_pairs):
            return "REPRESENTATIONAL"
    if symbol_pairs and not contains_refuted_pair(symbol_pairs):
        return "NAMING_ONLY"
    return "CONTENT"


def _matrix_from_nested(value: object, label: str) -> sp.Matrix:
    if not isinstance(value, tuple) or not value:
        raise ValueError(f"{label} is not a nonempty ordered collection of rows")
    if not all(isinstance(row, tuple) and row for row in value):
        raise ValueError(f"{label} does not contain nonempty ordered row vectors")
    widths = {len(row) for row in value}
    if len(widths) != 1:
        raise ValueError(f"{label} has ragged row vectors")
    return sp.Matrix(value)


def _residual_from_zero(value: object) -> object:
    if isinstance(value, tuple):
        return tuple(_residual_from_zero(item) for item in value)
    return residual(value, sp.Integer(0))


def compare_nullspace_basis(
    name: str,
    py_record: Record,
    wl_record: Record,
    py_matrix_record: Record,
    wl_matrix_record: Record,
) -> Comparison:
    try:
        parsed_py = parse_sympy_payload(py_record.raw)
        parsed_py_matrix = parse_sympy_payload(py_matrix_record.raw)
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, None, None, None, False,
            f"PY_PARSE_{type(error).__name__}: {error}", (),
        )
    try:
        parsed_wl = parse_wolfram_payload(wl_record.raw)
        parsed_wl_matrix = parse_wolfram_payload(wl_matrix_record.raw)
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, None, None, None, False,
            f"WL_PARSE_{type(error).__name__}: {error}", (),
        )
    try:
        py_basis = _matrix_from_nested(normalize(parsed_py), "PY nullspace basis")
        wl_basis = _matrix_from_nested(normalize(parsed_wl), "WL nullspace basis")
        py_matrix = _matrix_from_nested(normalize(parsed_py_matrix), "PY N1 matrix")
        wl_matrix = _matrix_from_nested(normalize(parsed_wl_matrix), "WL N1 matrix")
        py_canonical = py_basis.rref()[0]
        wl_canonical = wl_basis.rref()[0]
        if py_matrix.cols != py_canonical.cols:
            raise ValueError(
                f"PY N1 matrix width {py_matrix.cols} does not match basis width {py_canonical.cols}"
            )
        if wl_matrix.cols != wl_canonical.cols:
            raise ValueError(
                f"WL N1 matrix width {wl_matrix.cols} does not match basis width {wl_canonical.cols}"
            )

        py_canonical_object = normalize(py_canonical)
        wl_canonical_object = normalize(wl_canonical)
        subspace_difference = residual(py_canonical_object, wl_canonical_object)
        py_membership = _residual_from_zero(normalize(py_matrix * py_canonical.T))
        wl_membership = _residual_from_zero(normalize(wl_matrix * wl_canonical.T))
        difference = (subspace_difference, py_membership, wl_membership)
        passes = residual_is_zero(difference)
        return Comparison(
            name,
            py_record,
            wl_record,
            py_canonical_object,
            wl_canonical_object,
            difference,
            passes,
            None,
            (),
            None if passes else disagreement_kind(py_canonical_object, wl_canonical_object, ()),
        )
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, None, None, None, False,
            f"NORMALIZE_{type(error).__name__}: {error}", (),
        )


def compare_records(name: str, py_record: Record, wl_record: Record) -> Comparison:
    try:
        parsed_py = parse_sympy_payload(py_record.raw)
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, None, None, None, False,
            f"PY_PARSE_{type(error).__name__}: {error}", (),
        )
    try:
        parsed_wl = parse_wolfram_payload(wl_record.raw)
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, None, None, None, False,
            f"WL_PARSE_{type(error).__name__}: {error}", (),
        )
    try:
        py_object = normalize(parsed_py)
        wl_object = normalize(parsed_wl)
        difference = residual(py_object, wl_object)
        passes = residual_is_zero(difference)
        pairs = ()
        if not passes and not name.endswith(("_ROUTE", "_ROUTES")):
            pairs = infer_nonmechanical_symbol_pairs(
                py_object,
                wl_object,
                _symbol_spellings(parsed_py, mechanical=True),
            )
        return Comparison(
            name,
            py_record,
            wl_record,
            py_object,
            wl_object,
            difference,
            passes,
            None,
            pairs,
            None if passes else disagreement_kind(py_object, wl_object, pairs),
        )
    except Exception as error:
        return Comparison(
            name, py_record, wl_record, None, None, None, False,
            f"NORMALIZE_{type(error).__name__}: {error}", (),
        )


def render_value(value: object) -> str:
    if isinstance(value, Mismatch):
        return (
            f"{value.kind}(PY={render_value(value.py_value)}, "
            f"WL={render_value(value.wl_value)})"
        )
    if isinstance(value, tuple):
        contents = ", ".join(render_value(item) for item in value)
        if len(value) == 1:
            contents += ","
        return f"({contents})"
    if isinstance(value, sp.Basic):
        return sp.sstr(value)
    return repr(value)


def render_comparison(
    comparison: Comparison,
    refuted_symbol_pairs: frozenset[tuple[str, str, str]] | None = None,
) -> list[str]:
    lines = [
        f"COMPARE_NAME: {comparison.name}",
        f"PY_OPERAND: {comparison.py_record.raw}",
        f"WL_OPERAND: {comparison.wl_record.raw}",
    ]
    if comparison.reason is not None:
        lines.extend((f"UNCOMPARED_REASON: {comparison.reason}", "RESIDUAL: NOT_COMPUTED", "GUARD: FAIL"))
    else:
        if comparison.name.endswith(NULLSPACE_BASIS_SUFFIX):
            subspace_residual, py_membership, wl_membership = comparison.residual
            lines.extend(
                (
                    f"PY_CANONICAL_BASIS_SPAN: {render_value(comparison.py_object)}",
                    f"WL_CANONICAL_BASIS_SPAN: {render_value(comparison.wl_object)}",
                    f"SUBSPACE_RESIDUAL: {render_value(subspace_residual)}",
                    f"PY_N1_MEMBERSHIP_RESIDUAL: {render_value(py_membership)}",
                    f"WL_N1_MEMBERSHIP_RESIDUAL: {render_value(wl_membership)}",
                )
            )
        lines.append(f"RESIDUAL: {render_value(comparison.residual)}")
        for py_original, py_mechanical, wl_name in comparison.symbol_pairs:
            pair = (py_original, py_mechanical, wl_name)
            is_refuted = refuted_symbol_pairs is not None and pair in refuted_symbol_pairs
            label = (
                "D12_SYMBOL_NAME_REFUTED"
                if is_refuted
                else "D12_SYMBOL_NAME_CANDIDATE"
            )
            if refuted_symbol_pairs is None:
                status = "TRANSCRIPT_REFUTATION_NOT_ASSESSED"
            elif is_refuted:
                status = "TRANSCRIPT_REFUTED"
            else:
                status = "UNREFUTED"
            lines.append(
                f"{label}: "
                f"PY={py_original} PY_MECHANICAL={py_mechanical} WL={wl_name} "
                f"STATUS={status}"
            )
        if comparison.disagreement_kind is not None:
            lines.append(f"DISAGREEMENT_KIND: {comparison.disagreement_kind}")
        lines.append(f"GUARD: {'PASS' if comparison.passes else 'FAIL'}")
    return lines


def _render_duplicate_name(py: Transcript, wl: Transcript, name: str) -> list[str]:
    lines = [f"DUPLICATE_NAME: {name}"]
    for transcript in (py, wl):
        rows = transcript.duplicates.get(name)
        if rows is not None:
            for row in rows:
                lines.append(
                    f"{transcript.engine}_OPERAND_AT_LINE_{row.line_number}: {row.raw}"
                )
        elif name in transcript.records:
            row = transcript.records[name]
            lines.append(
                f"{transcript.engine}_COUNTERPART_OPERAND_AT_LINE_{row.line_number}: {row.raw}"
            )
    lines.extend(("UNCOMPARED_REASON: duplicate emitted name", "RESIDUAL: NOT_COMPUTED", "GUARD: FAIL"))
    return lines


def _symbol_pair_evidence(
    comparisons: Sequence[Comparison],
) -> tuple[
    dict[tuple[str, str, str], list[str]],
    dict[
        tuple[str, str, str],
        tuple[list[str], list[str], list[str], list[str]],
    ],
]:
    candidate_evidence: dict[tuple[str, str, str], list[str]] = {}
    for comparison in comparisons:
        candidate_pairs = list(comparison.symbol_pairs)
        if _has_equality(comparison.py_object) != _has_equality(comparison.wl_object):
            projected_py = _equation_system_as_expressions(comparison.py_object)
            projected_wl = _equation_system_as_expressions(comparison.wl_object)
            py_spellings = _symbol_spellings(
                parse_sympy_payload(comparison.py_record.raw),
                mechanical=True,
            )
            candidate_pairs.extend(
                infer_nonmechanical_symbol_pairs(
                    projected_py,
                    projected_wl,
                    py_spellings,
                )
            )
        for pair in dict.fromkeys(candidate_pairs):
            candidate_evidence.setdefault(pair, []).append(comparison.name)

    distinction_evidence: dict[
        tuple[str, str, str],
        tuple[list[str], list[str], list[str], list[str]],
    ] = {}
    for pair, candidate_names in candidate_evidence.items():
        _, py_mechanical, wl_name = pair
        py_source_names: list[str] = []
        py_target_names: list[str] = []
        wl_source_names: list[str] = []
        wl_target_names: list[str] = []
        for comparison in comparisons:
            if comparison.name in candidate_names or comparison.name.endswith(("_ROUTE", "_ROUTES")):
                continue
            py_symbols = {symbol.name for symbol in _nested_free_symbols(comparison.py_object)}
            wl_symbols = {symbol.name for symbol in _nested_free_symbols(comparison.wl_object)}
            if py_mechanical in py_symbols:
                py_source_names.append(comparison.name)
            if wl_name in py_symbols:
                py_target_names.append(comparison.name)
            if py_mechanical in wl_symbols:
                wl_source_names.append(comparison.name)
            if wl_name in wl_symbols:
                wl_target_names.append(comparison.name)
        distinction_evidence[pair] = (
            py_source_names,
            py_target_names,
            wl_source_names,
            wl_target_names,
        )
    return candidate_evidence, distinction_evidence


def run(py: Transcript, wl: Transcript) -> int:
    print(f"PY_INPUT: {py.path}")
    print(f"WL_INPUT: {wl.path}")
    print("JOIN_RULE: strip PY_/WL_ and match the remaining name exactly")
    print("SYMBOL_RULE: snake_case to lowerCamel mechanically, with no exceptions")
    print("SEQUENCE_RULE: ordered sequences share one recursive representation; a one-column Matrix is a sequence")
    print(
        "NULLSPACE_BASIS_RULE: compare canonical row spans and each span's action under its emitted N1 matrix"
    )

    duplicate_names = set(py.duplicates) | set(wl.duplicates)
    shared_names = set(py.records) & set(wl.records)
    shared_duplicate_names = shared_names & duplicate_names
    ordinary_shared = sorted(shared_names - shared_duplicate_names)
    nullspace_basis_names = [name for name in ordinary_shared if name.endswith(NULLSPACE_BASIS_SUFFIX)]
    ordinary_names = [name for name in ordinary_shared if name not in nullspace_basis_names]
    comparisons = [compare_records(name, py.records[name], wl.records[name]) for name in ordinary_names]
    for name in nullspace_basis_names:
        matrix_name = name.removesuffix(NULLSPACE_BASIS_SUFFIX) + "N1_MATRIX"
        missing = [
            engine
            for engine, transcript in (("PY", py), ("WL", wl))
            if matrix_name not in transcript.records or matrix_name in transcript.duplicates
        ]
        if missing:
            comparisons.append(
                Comparison(
                    name,
                    py.records[name],
                    wl.records[name],
                    None,
                    None,
                    None,
                    False,
                    f"N1_MATRIX_UNAVAILABLE_OR_DUPLICATE: {', '.join(missing)} {matrix_name}",
                    (),
                )
            )
        else:
            comparisons.append(
                compare_nullspace_basis(
                    name,
                    py.records[name],
                    wl.records[name],
                    py.records[matrix_name],
                    wl.records[matrix_name],
                )
            )
    comparisons.sort(key=lambda comparison: comparison.name)

    parsed_comparisons = [comparison for comparison in comparisons if comparison.reason is None]
    unparsed = [comparison for comparison in comparisons if comparison.reason is not None]
    symbol_evidence, distinction_evidence = _symbol_pair_evidence(parsed_comparisons)
    refuted_symbol_pairs = frozenset(
        pair
        for pair, evidence_names in distinction_evidence.items()
        if all(evidence_names)
    )
    parsed_comparisons = [
        replace(
            comparison,
            disagreement_kind=disagreement_kind(
                comparison.py_object,
                comparison.wl_object,
                comparison.symbol_pairs,
                refuted_symbol_pairs,
            ),
        )
        if not comparison.passes
        else comparison
        for comparison in parsed_comparisons
    ]

    print("CATEGORY: COMPARED_SHARED_NAMES")
    for comparison in parsed_comparisons:
        for line in render_comparison(comparison, refuted_symbol_pairs):
            print(line)

    print("CATEGORY: UNPARSED_SHARED_NAMES")
    for comparison in unparsed:
        for line in render_comparison(comparison):
            print(line)

    print("CATEGORY: DUPLICATE_NAMES")
    for name in sorted(duplicate_names):
        for line in _render_duplicate_name(py, wl, name):
            print(line)

    py_only = sorted(set(py.records) - set(wl.records))
    wl_only = sorted(set(wl.records) - set(py.records))
    print("CATEGORY: PY_ONLY_INVENTORY")
    print(f"PY_ONLY_COUNT: {len(py_only)}")
    for name in py_only:
        print(f"PY_ONLY_NAME: {name}")
    print("CATEGORY: WL_ONLY_INVENTORY")
    print(f"WL_ONLY_COUNT: {len(wl_only)}")
    for name in wl_only:
        print(f"WL_ONLY_NAME: {name}")

    print("CATEGORY: TRANSCRIPT_FORMAT_ISSUES")
    for issue in py.format_issues + wl.format_issues:
        print(f"FORMAT_ISSUE: {issue}")
        print("UNCOMPARED_REASON: nonempty transcript line has no emitted-name row grammar")
        print("RESIDUAL: NOT_COMPUTED")
        print("GUARD: FAIL")

    print("CATEGORY: D12_NONMECHANICAL_SYMBOL_WORKLIST")
    for pair, names in sorted(symbol_evidence.items()):
        if pair in refuted_symbol_pairs:
            continue
        py_original, py_mechanical, wl_name = pair
        print(
            "D12_SYMBOL_PAIR: "
            f"PY={py_original} PY_MECHANICAL={py_mechanical} WL={wl_name} "
            "STATUS=UNREFUTED "
            f"EVIDENCE_NAMES=({', '.join(names)})"
        )
    print("CATEGORY: D12_TRANSCRIPT_REFUTED_SYMBOL_PAIRS")
    for pair in sorted(refuted_symbol_pairs):
        py_original, py_mechanical, wl_name = pair
        py_source_names, py_target_names, wl_source_names, wl_target_names = distinction_evidence[pair]
        print(
            "D12_REFUTED_SYMBOL_PAIR: "
            f"PY={py_original} PY_MECHANICAL={py_mechanical} WL={wl_name} "
            f"CANDIDATE_EVIDENCE_NAMES=({', '.join(symbol_evidence[pair])}) "
            f"PY_SOURCE_EVIDENCE_NAMES=({', '.join(py_source_names)}) "
            f"PY_TARGET_EVIDENCE_NAMES=({', '.join(py_target_names)}) "
            f"WL_SOURCE_EVIDENCE_NAMES=({', '.join(wl_source_names)}) "
            f"WL_TARGET_EVIDENCE_NAMES=({', '.join(wl_target_names)})"
        )

    agreements = [comparison for comparison in parsed_comparisons if comparison.passes]
    disagreements = [comparison for comparison in parsed_comparisons if not comparison.passes]
    naming_only = [
        comparison for comparison in disagreements if comparison.disagreement_kind == "NAMING_ONLY"
    ]
    representational = [
        comparison
        for comparison in disagreements
        if comparison.disagreement_kind == "REPRESENTATIONAL"
    ]
    content_divergences = [
        comparison for comparison in disagreements if comparison.disagreement_kind == "CONTENT"
    ]
    content_shape_or_type = [
        comparison
        for comparison in content_divergences
        if _contains_mismatch(comparison.residual)
    ]
    content_route_tokens = [
        comparison
        for comparison in content_divergences
        if not _contains_mismatch(comparison.residual)
        and comparison.name.endswith(("_ROUTE", "_ROUTES"))
    ]
    content_numeric_or_algebraic = [
        comparison
        for comparison in content_divergences
        if not _contains_mismatch(comparison.residual)
        and not comparison.name.endswith(("_ROUTE", "_ROUTES"))
    ]
    bare_integer_agreements = [
        comparison for comparison in agreements if isinstance(comparison.py_object, sp.Integer)
    ]
    empty_container_agreements = [
        comparison
        for comparison in agreements
        if comparison.py_object == () and comparison.wl_object == ()
    ]
    symbolic_or_structured_agreements = [
        comparison
        for comparison in agreements
        if not isinstance(comparison.py_object, sp.Integer)
        and not (comparison.py_object == () and comparison.wl_object == ())
    ]
    format_issue_count = len(py.format_issues) + len(wl.format_issues)
    duplicate_count = sum(len(rows) - 1 for rows in py.duplicates.values()) + sum(
        len(rows) - 1 for rows in wl.duplicates.values()
    )
    shared_verdict_count = len(parsed_comparisons) + len(unparsed) + len(shared_duplicate_names)
    print("CATEGORY: SUMMARY")
    print(f"PY_NAME_COUNT: {len(py.records)}")
    print(f"WL_NAME_COUNT: {len(wl.records)}")
    print(f"SHARED_NAME_COUNT: {len(shared_names)}")
    print(f"COMPARED_NAME_COUNT: {len(parsed_comparisons)}")
    print(f"AGREEMENT_COUNT: {len(agreements)}")
    print(f"AGREEMENT_BARE_INTEGER_COUNT: {len(bare_integer_agreements)}")
    print(f"AGREEMENT_EMPTY_CONTAINER_COUNT: {len(empty_container_agreements)}")
    print(
        "AGREEMENT_SYMBOLIC_OR_STRUCTURED_COUNT: "
        f"{len(symbolic_or_structured_agreements)}"
    )
    print(f"DISAGREEMENT_COUNT: {len(disagreements)}")
    print(f"NAMING_ONLY_DISAGREEMENT_COUNT: {len(naming_only)}")
    print(f"REPRESENTATIONAL_DISAGREEMENT_COUNT: {len(representational)}")
    print(f"CONTENT_DIVERGENCE_COUNT: {len(content_divergences)}")
    print(f"CONTENT_SHAPE_OR_TYPE_MISMATCH_COUNT: {len(content_shape_or_type)}")
    print(f"CONTENT_ROUTE_TOKEN_COUNT: {len(content_route_tokens)}")
    print(
        "CONTENT_NUMERIC_OR_ALGEBRAIC_RESIDUAL_COUNT: "
        f"{len(content_numeric_or_algebraic)}"
    )
    print(f"UNPARSED_SHARED_COUNT: {len(unparsed)}")
    print(f"NULLSPACE_BASIS_COMPARED_COUNT: {len(nullspace_basis_names)}")
    print(f"D12_SYMBOL_PAIR_COUNT: {len(symbol_evidence) - len(refuted_symbol_pairs)}")
    print(f"D12_REFUTED_SYMBOL_PAIR_COUNT: {len(refuted_symbol_pairs)}")
    print(f"DUPLICATE_NAME_COUNT: {len(duplicate_names)}")
    print(f"SHARED_DUPLICATE_NAME_COUNT: {len(shared_duplicate_names)}")
    print(f"DUPLICATE_ROW_COUNT: {duplicate_count}")
    print(f"FORMAT_ISSUE_COUNT: {format_issue_count}")
    print(f"SHARED_VERDICT_COUNT: {shared_verdict_count}")

    failed = bool(
        disagreements
        or unparsed
        or duplicate_count
        or format_issue_count
    )
    print(f"FINAL_GUARD: {'FAIL' if failed else 'PASS'}")
    return int(failed)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Join exactly one PY and one WL .out transcript by emitted object name."
    )
    parser.add_argument("out_files", nargs=2, type=Path, metavar="OUT")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        py, wl = load_pair(args.out_files)
    except ComparatorInputError as error:
        print(f"INPUT_ERROR: {error}")
        print("FINAL_GUARD: FAIL")
        return 2
    return run(py, wl)


if __name__ == "__main__":
    raise SystemExit(main())
