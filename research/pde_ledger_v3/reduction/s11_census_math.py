#!/usr/bin/env python3
"""Kernel-free symbolic machinery shared by the S11 census instruments.

The producer routes are not reused.  In particular, completeness is built from
factor-cover decomposition followed by ``nonlinsolve``; the PY producer uses
``solve`` and the WL producer uses ``Solve``/a custom factor-lattice fallback.
No expression in this file is a physics expression: every symbolic operand comes
from a record payload supplied by the caller.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import itertools
import json
import math
import re
import resource
import signal
import time
from typing import Any, Iterable, Iterator, Mapping, Sequence

import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.parsing.mathematica import MathematicaParser


# Floors required by brief obligation 7.  These apply independently to each
# semantic locus, witness, and undecided-record operation.
OPERATION_TIME_BUDGET_SECONDS = 60
OPERATION_MEMORY_BUDGET_BYTES = 512 * 1024 * 1024
# Constant algebraic decisions are a mandatory sub-route of every operation.
# Bound an individual minimal-polynomial attempt so a hard algebraic number
# cannot consume the enclosing 60-second record budget without leaving a
# classifier-level expiry token in the printed residual statuses.
EXACT_ROUTE_TIME_BUDGET_SECONDS = 2.0
MAX_SHEET_ASSIGNMENTS = 256
MAX_FACTOR_COVERS = 256
MAX_SAMPLE_ASSIGNMENTS = 96
MAX_LITERAL_CHARS = 4096


@dataclass(frozen=True)
class ParsedSolution:
    variables: tuple[sp.Symbol, ...]
    branches: tuple[dict[sp.Expr, sp.Expr], ...]
    kind: str


def enforce_memory_budget() -> None:
    resource.setrlimit(
        resource.RLIMIT_AS,
        (OPERATION_MEMORY_BUDGET_BYTES, OPERATION_MEMORY_BUDGET_BYTES),
    )


def digest_text(text: str) -> str:
    return sha256(text.encode("utf-8")).hexdigest()


def quoted(text: object) -> str:
    return json.dumps(str(text), ensure_ascii=True, separators=(",", ":"))


def render_text(text: str) -> str:
    if len(text) <= MAX_LITERAL_CHARS:
        return quoted(text)
    preview = text[:240]
    return (
        "{bytes="
        + str(len(text.encode("utf-8")))
        + " sha256="
        + digest_text(text)
        + " preview="
        + quoted(preview)
        + "}"
    )


def render_object(value: object) -> str:
    try:
        text = sp.srepr(value)  # type: ignore[arg-type]
    except BaseException:
        text = repr(value)
    return render_text(text)


def _scan_top(text: str, separator: str) -> Iterator[int]:
    """Yield separator positions outside strings and all bracketed objects."""
    stack: list[str] = []
    in_string: str | None = None
    escaped = False
    pairs = {"(": ")", "[": "]", "{": "}"}
    index = 0
    while index < len(text):
        char = text[index]
        if in_string is not None:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == in_string:
                in_string = None
            index += 1
            continue
        if char in ("'", '"'):
            in_string = char
            index += 1
            continue
        if text.startswith("<|", index):
            stack.append("|>")
            index += 2
            continue
        if stack and text.startswith(stack[-1], index):
            index += len(stack.pop())
            continue
        if char in pairs:
            stack.append(pairs[char])
            index += 1
            continue
        if stack and char == stack[-1]:
            stack.pop()
            index += 1
            continue
        if not stack and text.startswith(separator, index):
            yield index
            index += len(separator)
            continue
        index += 1


def split_delimited(text: str, opening: str = "{", closing: str = "}") -> list[str]:
    source = text.strip()
    if not (source.startswith(opening) and source.endswith(closing)):
        raise ValueError(f"not a {opening}{closing} object: {source[:240]}")
    body = source[len(opening) : -len(closing)]
    positions = list(_scan_top(body, ","))
    if not positions:
        return [] if not body.strip() else [body.strip()]
    result: list[str] = []
    start = 0
    for position in positions:
        result.append(body[start:position].strip())
        start = position + 1
    result.append(body[start:].strip())
    return result


def split_top_relation(text: str, token: str) -> tuple[str, str]:
    positions = list(_scan_top(text, token))
    if len(positions) != 1:
        raise ValueError(
            f"expected one top-level {token!r}, found {len(positions)}: {text[:240]}"
        )
    position = positions[0]
    return text[:position].strip(), text[position + len(token) :].strip()


_WL_PROTECTED_SYMBOLS = {"beta": "s11BetaSymbol"}

# ``parse_mathematica`` deliberately leaves some Mathematica heads as undefined
# SymPy functions.  That is safe for opaque provenance heads, but not for heads
# whose semantics the census evaluates.  Keep this table explicit and small so
# every evaluated head is audited in one place.  Unequal is rewritten to a
# Boolean expression before parsing because an undefined function cannot be an
# And/Or operand.
_WL_EVALUATED_HEADS = {
    "Abs": sp.Abs,
    "Sign": sp.sign,
    "Re": sp.re,
    "Im": sp.im,
    "Conjugate": sp.conjugate,
    "Sin": sp.sin,
    "Cos": sp.cos,
    "Tan": sp.tan,
    "Exp": sp.exp,
    "Log": sp.log,
    "Min": sp.Min,
    "Max": sp.Max,
}


def _replace_wl_empty_lists(text: str) -> str:
    """Protect nested literal {} objects from the Mathematica parser."""
    result: list[str] = []
    index = 0
    in_string = False
    escaped = False
    while index < len(text):
        char = text[index]
        if in_string:
            result.append(char)
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                in_string = False
            index += 1
            continue
        if char == '"':
            in_string = True
            result.append(char)
            index += 1
            continue
        if text.startswith("{}", index):
            result.append("s11EmptyList[0]")
            index += 2
            continue
        result.append(char)
        index += 1
    return "".join(result)


def _wl_unequal(*operands: sp.Basic) -> sp.Basic:
    """Mathematica Unequal is pairwise inequality, including infix ``!=``."""
    if len(operands) < 2:
        raise ValueError("Unequal requires at least two operands")
    return sp.And(
        *(
            sp.Ne(operands[left], operands[right])
            for left in range(len(operands))
            for right in range(left + 1, len(operands))
        )
    )


def _parse_mathematica_evaluated(source: str) -> sp.Basic:
    """Parse with every Boolean-producing record head bound before parents."""
    parser = MathematicaParser()
    # Infix ``!=`` and explicit Unequal[...] both become the same FullForm
    # node.  It must be Boolean before an enclosing And/Or is constructed.
    parser._node_conversions["Unequal"] = _wl_unequal
    # Witness premises use these private heads when Element is deliberately
    # deferred.  Contains is a Boolean object, so it too must be constructed
    # before a surrounding And is built.  evaluate=False is essential here:
    # principal-value membership is decided only after witness substitution.
    parser._node_conversions["s11ElementReals"] = lambda value: sp.Contains(
        value, sp.S.Reals, evaluate=False
    )
    parser._node_conversions["s11ElementIntegers"] = lambda value: sp.Contains(
        value, sp.S.Integers, evaluate=False
    )
    return parser.parse(source)


def _bind_wl_evaluated_heads(value: sp.Basic) -> sp.Basic:
    value = value.replace(
        lambda item: isinstance(item, AppliedUndef)
        and item.func.__name__ == "s11EmptyList",
        lambda _item: sp.Tuple(),
    )

    def is_bound_head(item: sp.Basic) -> bool:
        return isinstance(item, AppliedUndef) and item.func.__name__ in _WL_EVALUATED_HEADS

    def bind(item: sp.Basic) -> sp.Basic:
        constructor = _WL_EVALUATED_HEADS[item.func.__name__]
        return constructor(*item.args)

    return value.replace(is_bound_head, bind)


def _replace_wl_elements(text: str, *, deferred: bool = False) -> str:
    source = text
    search_from = 0
    while True:
        start = source.find("Element[", search_from)
        if start < 0:
            return source
        depth = 0
        in_string = False
        escaped = False
        end = -1
        for index in range(start + len("Element"), len(source)):
            char = source[index]
            if in_string:
                if escaped:
                    escaped = False
                elif char == "\\":
                    escaped = True
                elif char == '"':
                    in_string = False
                continue
            if char == '"':
                in_string = True
            elif char == "[":
                depth += 1
            elif char == "]":
                depth -= 1
                if depth == 0:
                    end = index
                    break
        if end < 0:
            raise ValueError(f"unclosed Element head: {source[start:start + 240]}")
        body = source[start + len("Element[") : end]
        parts = list(_scan_top(body, ","))
        if len(parts) != 1:
            raise ValueError(f"unrecognized Element operands: {body[:240]}")
        comma = parts[0]
        expression = body[:comma].strip()
        domain = body[comma + 1 :].strip()
        if domain == "Reals":
            replacement = (
                f"s11ElementReals[{expression}]"
                if deferred
                else f"(Im[{expression}] == 0)"
            )
        elif domain == "Integers":
            replacement = (
                f"s11ElementIntegers[{expression}]"
                if deferred
                else f"(Im[{expression}] == 0 && Sin[Pi*({expression})] == 0)"
            )
        else:
            raise ValueError(f"unsupported Element domain: {domain}")
        source = source[:start] + replacement + source[end + 1 :]
        search_from = start + len(replacement)


def parse_wl(text: str, *, defer_elements: bool = False) -> sp.Basic:
    # SymPy's Mathematica translator represents unknown Element[...] heads as
    # non-Boolean functions, which then cannot be children of And.  Translate
    # the two record domains exactly to Boolean relations before parsing.
    stripped = text.strip()
    # The Mathematica parser rejects nested empty lists (notably the 62 real
    # ``{{}}`` solution payloads).  Parse list structure recursively; each
    # non-list leaf still goes through the same Mathematica expression parser.
    if stripped.startswith("{") and stripped.endswith("}"):
        return sp.Tuple(
            *(parse_wl(item, defer_elements=defer_elements) for item in split_delimited(stripped))
        )
    source = _replace_wl_empty_lists(
        _replace_wl_elements(stripped, deferred=defer_elements)
    )
    for original, protected in _WL_PROTECTED_SYMBOLS.items():
        source = re.sub(rf"\b{re.escape(original)}\b", protected, source)
    parsed = _parse_mathematica_evaluated(source)
    # Literal Mathematica True/False are returned as Python bools, unlike
    # evaluated relations which return SymPy booleans.
    if parsed is True:
        return sp.true
    if parsed is False:
        return sp.false
    replacements = {
        sp.Symbol(protected): sp.Symbol(original)
        for original, protected in _WL_PROTECTED_SYMBOLS.items()
    }
    return _bind_wl_evaluated_heads(parsed.xreplace(replacements))


_PY_ASSUMED_SYMBOL_RE = re.compile(
    r"Symbol\((?P<name>'(?:\\.|[^'])*'|\"(?:\\.|[^\"])*\")"
    r"(?:,\s*[A-Za-z_][A-Za-z0-9_]*=(?:True|False|None))+\)"
)


def parse_py(text: str, *, assumption_free: bool = False) -> sp.Basic:
    # srepr is executable constructor syntax, not user-authored source.  sympify
    # supplies the SymPy constructor namespace and evaluate=False preserves the
    # record's unevaluated relation forms where SymPy supports doing so.
    source = text
    if assumption_free:
        source = _PY_ASSUMED_SYMBOL_RE.sub(
            lambda match: f"Symbol({match.group('name')})", source
        )
    return sp.sympify(source, evaluate=False)


def parse_payload(
    dialect: str,
    text: str,
    *,
    assumption_free: bool = False,
    defer_wl_elements: bool = False,
) -> sp.Basic:
    if dialect == "WL":
        return parse_wl(text, defer_elements=defer_wl_elements)
    if dialect == "PY":
        return parse_py(text, assumption_free=assumption_free)
    raise ValueError(f"unknown dialect {dialect!r}")


def string_value(value: object) -> str | None:
    if isinstance(value, sp.Basic) and value.func.__name__ in {"Str", "_Str"}:
        if value.args:
            return str(value.args[0])
    return None


def association(value: object) -> dict[str, object]:
    result: dict[str, object] = {}
    if not isinstance(value, (tuple, list, sp.Tuple)):
        return result
    for entry in value:
        if isinstance(entry, sp.Basic) and entry.func.__name__ == "Rule" and len(entry.args) == 2:
            key = string_value(entry.args[0])
            if key is not None:
                result[key] = entry.args[1]
        elif isinstance(entry, (tuple, list, sp.Tuple)) and len(entry) == 2:
            key = string_value(entry[0])
            if key is not None:
                result[key] = entry[1]
    return result


def _wl_solution_fields(payload: str) -> tuple[str, str]:
    fields: dict[str, str] = {}
    for item in split_delimited(payload):
        left, right = split_top_relation(item, "->")
        if left.startswith('"') and left.endswith('"'):
            fields[left[1:-1]] = right
        if "SOLVE_VARIABLES" in fields and "SOLUTION_SET" in fields:
            # All later fields are attempt provenance.  Their delimiters were
            # validated by split_delimited, but they are intentionally not fed
            # into membership computation.
            continue
    return fields["SOLVE_VARIABLES"], fields["SOLUTION_SET"]


def _rule_mapping(branch: object) -> dict[sp.Expr, sp.Expr]:
    if not isinstance(branch, (tuple, list, sp.Tuple)):
        raise ValueError(f"branch is not a rule sequence: {repr(branch)[:240]}")
    result: dict[sp.Expr, sp.Expr] = {}
    for entry in branch:
        if isinstance(entry, sp.Basic) and entry.func.__name__ == "Rule" and len(entry.args) == 2:
            result[entry.args[0]] = entry.args[1]
        elif isinstance(entry, (tuple, list, sp.Tuple)) and len(entry) == 2:
            if not isinstance(entry[0], sp.Basic):
                raise ValueError(f"branch lhs is not symbolic: {repr(entry)[:240]}")
            result[entry[0]] = entry[1]
        else:
            raise ValueError(f"unrecognized branch entry: {repr(entry)[:240]}")
    return result


def parse_solution(dialect: str, payload: str) -> ParsedSolution:
    if dialect == "WL":
        variables_text, solutions_text = _wl_solution_fields(payload)
        variables_obj = parse_wl(variables_text)
        solution_obj = parse_wl(solutions_text)
        variables = tuple(variables_obj) if isinstance(variables_obj, sp.Tuple) else (variables_obj,)
        if not isinstance(solution_obj, sp.Tuple):
            # ConditionSet/FiniteSet/EmptySet are parsed objects rather than
            # opaque text.  They carry an explicit non-branch kind.
            return ParsedSolution(
                variables=tuple(v for v in variables if isinstance(v, sp.Symbol)),
                branches=(),
                kind=solution_obj.func.__name__,
            )
        branches = tuple(_rule_mapping(branch) for branch in solution_obj)
        return ParsedSolution(
            variables=tuple(v for v in variables if isinstance(v, sp.Symbol)),
            branches=branches,
            kind="BRANCH_LIST",
        )

    solution_obj = parse_py(payload)
    if not isinstance(solution_obj, sp.Tuple):
        return ParsedSolution(variables=(), branches=(), kind=solution_obj.func.__name__)
    branches = tuple(_rule_mapping(branch) for branch in solution_obj)
    variables: list[sp.Symbol] = []
    for branch in branches:
        for lhs in branch:
            if isinstance(lhs, sp.Symbol) and lhs not in variables:
                variables.append(lhs)
    return ParsedSolution(variables=tuple(variables), branches=branches, kind="BRANCH_LIST")


def relation_residual(value: object) -> sp.Expr:
    if value is True or value is sp.true:
        return sp.Integer(0)
    if value is False or value is sp.false:
        return sp.Integer(1)
    if isinstance(value, sp.Equality):
        return sp.Add(value.lhs, -value.rhs, evaluate=False)
    if isinstance(value, sp.Expr):
        return value
    raise ValueError(f"not an equation: {repr(value)[:240]}")


def equation_residuals(value: object) -> tuple[sp.Expr, ...]:
    if isinstance(value, (tuple, list, sp.Tuple)):
        return tuple(relation_residual(item) for item in value)
    return (relation_residual(value),)


def _radical_atoms(expression: sp.Basic) -> set[sp.Basic]:
    powers = {
        atom
        for atom in expression.atoms(sp.Pow)
        if atom.exp.is_Rational and atom.exp.q == 2 and atom.exp.p % 2 != 0
        and atom.base.is_number is not True
    }
    absolute_values = {
        atom for atom in expression.atoms(sp.Abs) if atom.args[0].is_number is not True
    }
    return powers | absolute_values


def sheet_residuals(
    residuals: Sequence[sp.Expr], branch: Mapping[sp.Expr, sp.Expr]
) -> tuple[list[tuple[str, tuple[sp.Expr, ...]]], int, int]:
    # Obligation 5 order is semantic: substitute the emitted branch exactly as
    # written, simplify that residual, and only then enumerate radicals that
    # remain.  Radicals originating in a branch RHS are fixed parts of the
    # emitted object and are excluded from the global sheet sweep.
    branch_radicals: set[sp.Basic] = set()
    for rhs in branch.values():
        branch_radicals.update(_radical_atoms(rhs))
    substituted = tuple(
        simplify_residual(residual.subs(branch, simultaneous=True))[0]
        for residual in residuals
    )
    sheet_atoms: set[sp.Basic] = set()
    for residual in substituted:
        sheet_atoms.update(_radical_atoms(residual))
    sheet_atoms.difference_update(branch_radicals)
    ordered_atoms = sorted(sheet_atoms, key=sp.default_sort_key)
    total = 2 ** len(ordered_atoms)
    tested = min(total, MAX_SHEET_ASSIGNMENTS)
    result: list[tuple[str, tuple[sp.Expr, ...]]] = []
    assignments = itertools.islice(itertools.product((-1, 1), repeat=len(ordered_atoms)), tested)
    for signs in assignments:
        replacements: dict[sp.Basic, sp.Basic] = {}
        labels: list[str] = []
        for atom, sign in zip(ordered_atoms, signs):
            if isinstance(atom, sp.Abs):
                replacements[atom] = sign * atom.args[0]
                labels.append(f"Abs({sp.sstr(atom.args[0])}):{sign:+d}")
            else:
                replacements[atom] = sign * atom
                labels.append(f"sqrt({sp.sstr(atom.base)}):{sign:+d}")
        values = tuple(item.xreplace(replacements) for item in substituted)
        result.append(("[" + ",".join(labels) + "]", values))
    return result, tested, total


def simplify_residual(expression: sp.Expr) -> tuple[sp.Expr, str]:
    try:
        value = sp.cancel(sp.together(expression))
        if value.has(sp.zoo, sp.nan, sp.oo, -sp.oo):
            return value, "UNDEFINED"
        if value == 0 or value.is_zero is True:
            return sp.Integer(0), "ZERO"
        simplified = value
        # Full simplify can explode on the radical minors this census was built
        # to inspect.  It is still useful for compact residuals; large ones go
        # directly to the explicitly named exact-sampling fallback.
        if len(str(value)) <= 2000:
            simplified = sp.simplify(value)
            if simplified == 0 or simplified.is_zero is True:
                return sp.Integer(0), "ZERO"
            if simplified.has(sp.zoo, sp.nan, sp.oo, -sp.oo):
                return simplified, "UNDEFINED"
            if not simplified.free_symbols:
                return simplified, _constant_zero_status(simplified)
            if simplified.is_zero is False:
                return simplified, "NONZERO"
        elif not simplified.free_symbols:
            return simplified, _constant_zero_status(simplified)
        # A single exact nonzero specialization is a proof that a symbolic
        # residual is not identically zero.  It is only a fallback, and the
        # status/verdict names that sampling explicitly.
        sampled_count = 0
        sampled_zero = 0
        for assignment in itertools.islice(
            exact_sample_assignments(tuple(simplified.free_symbols)), 4
        ):
            try:
                sampled = sp.simplify(simplified.subs(assignment, simultaneous=True))
            except BaseException:
                continue
            if sampled.has(sp.zoo, sp.nan, sp.oo, -sp.oo):
                continue
            sampled_count += 1
            if sampled.free_symbols:
                continue
            sampled_status = _constant_zero_status(sampled)
            if sampled_status == "ZERO":
                sampled_zero += 1
                continue
            if sampled_status == "NONZERO":
                return simplified, "NONZERO_SAMPLED"
        if sampled_count and sampled_zero == sampled_count:
            # Samples can refute an identity, never prove one.  Preserve the
            # explicit diagnostic while leaving the symbolic result undecided.
            return simplified, "UNDECIDED_ZERO_SAMPLES"
        return simplified, "UNDECIDED"
    except BaseException:
        if not expression.free_symbols:
            try:
                return expression, _constant_zero_status(expression)
            except BaseException:
                # The exact route itself was entered but could not return a
                # certificate or its own expiry token.
                return expression, "UNDECIDED_EXACT_ROUTE_FAILED"
        return expression, "UNDECIDED"


class _ExactRouteExpired(Exception):
    """Internal alarm used to make exact-constant resource expiry explicit."""


def _constant_zero_status(expression: sp.Expr) -> str:
    """Classify a constant only from an exact certificate.

    Rational equality and SymPy's three-valued exact zero predicate are exact
    equals-zero decisions.  When that predicate is undecided, the mandatory
    exact route continues through a guarded minimal-polynomial calculation. A
    successful minimal polynomial has zero constant term exactly for the
    algebraic number zero; a nonzero constant term certifies a nonzero
    algebraic number.  Floats never use the predicate as an exact certificate.
    No structural ``expression != 0`` or fixed-precision midpoint is a
    certificate.
    """
    if expression.free_symbols:
        raise ValueError("constant classifier received a symbolic expression")
    if expression.is_Rational:
        return "ZERO" if expression == 0 else "NONZERO"
    if not expression.has(sp.Float):
        exact_zero = expression.is_zero
        if exact_zero is True:
            return "ZERO"
        if exact_zero is False:
            return "NONZERO"

    previous_handler = signal.getsignal(signal.SIGALRM)
    previous_delay, previous_interval = signal.getitimer(signal.ITIMER_REAL)
    started = time.monotonic()

    def expire(_signum: int, _frame: object) -> None:
        raise _ExactRouteExpired

    guard_seconds = EXACT_ROUTE_TIME_BUDGET_SECONDS
    if previous_delay > 0:
        guard_seconds = min(guard_seconds, previous_delay)
    try:
        signal.signal(signal.SIGALRM, expire)
        signal.setitimer(signal.ITIMER_REAL, guard_seconds)
        generator = sp.Dummy("s11_zero_certificate")
        polynomial = sp.minpoly(expression, generator, polys=True)
    except _ExactRouteExpired:
        return "UNDECIDED_EXACT_ROUTE_EXPIRED"
    except BaseException:
        # The minimal-polynomial route was attempted after the exact predicate
        # declined to decide.  Do not turn either failure into NONZERO.
        if expression.has(sp.Float):
            return "UNDECIDED_EXACT_ROUTE_FAILED"
        if expression.is_zero is True:
            return "ZERO"
        if expression.is_zero is False:
            return "NONZERO"
        return "UNDECIDED_EXACT_ROUTE_FAILED"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
        if previous_delay > 0:
            elapsed = time.monotonic() - started
            signal.setitimer(
                signal.ITIMER_REAL,
                max(1e-6, previous_delay - elapsed),
                previous_interval,
            )

    constant_term = polynomial.TC()
    if polynomial.degree() == 1 and constant_term == 0:
        return "ZERO"
    if constant_term.is_zero is False:
        return "NONZERO"
    return "UNDECIDED_EXACT_ROUTE_FAILED"


def branch_memberships(
    residuals: Sequence[sp.Expr], branch: Mapping[sp.Expr, sp.Expr]
) -> tuple[list[dict[str, object]], int, int]:
    sheets, tested, total = sheet_residuals(residuals, branch)
    rows: list[dict[str, object]] = []
    for label, values in sheets:
        simplified = tuple(simplify_residual(value) for value in values)
        statuses = tuple(status for _value, status in simplified)
        if "NONZERO" in statuses or "UNDEFINED" in statuses:
            verdict = "SPURIOUS_BRANCH"
        elif "NONZERO_SAMPLED" in statuses:
            verdict = "SPURIOUS_BRANCH_SAMPLED"
        elif any(status.startswith("UNDECIDED") for status in statuses):
            verdict = "BRANCH_MEMBERSHIP_UNDECIDED"
        else:
            verdict = "BRANCH_CONTAINED"
        rows.append(
            {
                "sheet": label,
                "residuals": tuple(value for value, _status in simplified),
                "definedness": statuses,
                "verdict": verdict,
            }
        )
    return rows, tested, total


def _factor_choices(residual: sp.Expr) -> tuple[sp.Expr, ...] | None:
    simplified, status = simplify_residual(residual)
    if status == "ZERO":
        return None
    if status == "UNDEFINED":
        return ()
    numerator, _denominator = sp.fraction(sp.together(simplified))
    factored = sp.factor(numerator)
    factors = factored.args if isinstance(factored, sp.Mul) else (factored,)
    choices: list[sp.Expr] = []
    for factor in factors:
        base = factor.base if isinstance(factor, sp.Pow) else factor
        if base.is_number:
            continue
        normalized = sp.cancel(base)
        if normalized not in choices:
            choices.append(normalized)
    return tuple(choices)


def minimal_factor_covers(residuals: Sequence[sp.Expr]) -> tuple[list[frozenset[sp.Expr]], bool]:
    covers: list[frozenset[sp.Expr]] = [frozenset()]
    truncated = False
    for residual in residuals:
        choices = _factor_choices(residual)
        if choices is None:
            continue
        if not choices:
            return [], truncated
        expanded = {cover | {choice} for cover in covers for choice in choices}
        ordered = sorted(expanded, key=lambda item: (len(item), tuple(map(sp.default_sort_key, item))))
        minimal: list[frozenset[sp.Expr]] = []
        for candidate in ordered:
            if any(existing <= candidate for existing in minimal):
                continue
            minimal.append(candidate)
            if len(minimal) == MAX_FACTOR_COVERS:
                truncated = len(ordered) > len(minimal)
                break
        covers = minimal
    return covers, truncated


def _solution_tuple_to_branch(
    variables: Sequence[sp.Symbol], values: Sequence[sp.Expr]
) -> dict[sp.Expr, sp.Expr]:
    result: dict[sp.Expr, sp.Expr] = {}
    for variable, value in zip(variables, values):
        if value != variable:
            result[variable] = value
    return result


def _branch_constraints(branch: Mapping[sp.Expr, sp.Expr]) -> tuple[sp.Expr, ...]:
    return tuple(lhs - rhs for lhs, rhs in branch.items())


def _defined_union_constraints(
    candidate: Mapping[sp.Expr, sp.Expr],
    emitted: Sequence[Mapping[sp.Expr, sp.Expr]],
) -> tuple[tuple[sp.Expr, ...], ...]:
    """Return emitted branch equations that are defined on the candidate.

    Definedness is branch-local.  A sibling containing a pole after candidate
    substitution is non-covering and is removed; it must never contaminate the
    product ideal for the remaining union.  An empty branch is retained because
    it denotes the full ambient chart.
    """
    result: list[tuple[sp.Expr, ...]] = []
    for branch in emitted:
        substituted = tuple(
            constraint.subs(candidate, simultaneous=True)
            for constraint in _branch_constraints(branch)
        )
        statuses = tuple(simplify_residual(item)[1] for item in substituted)
        if "UNDEFINED" in statuses:
            continue
        result.append(substituted)
    return tuple(result)


def _sample_union_coverage(
    constraints: Sequence[Sequence[sp.Expr]],
) -> tuple[str, dict[sp.Symbol, sp.Rational]]:
    """Seek an exact refuting point for containment in a branch union.

    Each emitted branch is a conjunction.  A defined branch covers a point
    only when every constraint is certified zero.  A point is certified
    uncovered only when every defined branch has at least one certified
    nonzero constraint; undefined branches are absent from the union.  Any
    other status pattern is undecided.  Samples can therefore refute union
    coverage, but can never prove it.
    """
    symbols = set().union(
        *(constraint.free_symbols for branch in constraints for constraint in branch)
    )
    for assignment in itertools.islice(exact_sample_assignments(tuple(symbols)), 16):
        covered = False
        every_defined_branch_refuted = True
        for branch in constraints:
            statuses = tuple(
                simplify_residual(
                    constraint.subs(assignment, simultaneous=True)
                )[1]
                for constraint in branch
            )
            if any(status == "UNDEFINED" for status in statuses):
                continue
            if all(status == "ZERO" for status in statuses):
                covered = True
                break
            if not any(status in {"NONZERO", "NONZERO_SAMPLED"} for status in statuses):
                every_defined_branch_refuted = False
        if not covered and every_defined_branch_refuted:
            return "NOT_COVERED_SAMPLED", assignment
    return "COVERAGE_UNDECIDED", {}


def _branch_is_covered(
    candidate: Mapping[sp.Expr, sp.Expr], emitted: Sequence[Mapping[sp.Expr, sp.Expr]]
) -> tuple[str, dict[sp.Symbol, sp.Rational]]:
    """Decide containment in the union of emitted branch varieties.

    The union is represented by products choosing one equation from each
    emitted branch.  All such products vanish exactly on the branch union.
    Substitution into the candidate chart therefore recognizes algebraically
    equivalent radical/Abs charts that no single syntactic branch matches.
    """
    constraints = _defined_union_constraints(candidate, emitted)
    if not constraints:
        return "NOT_COVERED", {}
    if any(not items for items in constraints):
        return "COVERED_ALGEBRAIC", {}
    combination_count = math.prod(len(items) for items in constraints)
    if combination_count <= MAX_FACTOR_COVERS:
        statuses: list[str] = []
        for factors in itertools.product(*constraints):
            product = sp.Mul(*factors)
            _value, status = simplify_residual(product)
            statuses.append(status)
        if statuses and all(status == "ZERO" for status in statuses):
            return "COVERED_ALGEBRAIC", {}
        if any(status in {"NONZERO", "NONZERO_SAMPLED", "UNDEFINED"} for status in statuses):
            return _sample_union_coverage(constraints)
    return _sample_union_coverage(constraints)


def _is_point_candidate(values: Sequence[object]) -> bool:
    return all(isinstance(value, sp.Expr) and not isinstance(value, sp.Set) for value in values)


def completeness_candidates(
    residuals: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    emitted: Sequence[Mapping[sp.Expr, sp.Expr]],
) -> dict[str, object]:
    if not variables:
        variables = tuple(sorted(set().union(*(r.free_symbols for r in residuals)), key=sp.default_sort_key))
    covers, truncated = minimal_factor_covers(residuals)
    candidates: list[dict[sp.Expr, sp.Expr]] = []
    unresolved = 0
    artifact_count = 0
    for cover in covers:
        if not cover:
            candidates.append({})
            continue
        try:
            solved = sp.nonlinsolve(tuple(cover), tuple(variables))
        except BaseException:
            unresolved += 1
            continue
        if isinstance(solved, sp.FiniteSet):
            for item in solved:
                values = tuple(item) if isinstance(item, (tuple, sp.Tuple)) else (item,)
                if len(values) == len(variables):
                    if not _is_point_candidate(values):
                        artifact_count += 1
                        continue
                    candidate = _solution_tuple_to_branch(variables, values)
                    if candidate not in candidates:
                        candidates.append(candidate)
                else:
                    unresolved += 1
        elif solved is sp.EmptySet:
            continue
        else:
            # ConditionSet and other symbolic returns were parsed and retained,
            # but do not license a completeness conclusion.
            unresolved += 1
    surviving: list[dict[sp.Expr, sp.Expr]] = []
    excluded_artifacts: list[dict[str, object]] = []
    for candidate in candidates:
        try:
            substitutions = tuple(
                residual.subs(candidate, simultaneous=True) for residual in residuals
            )
            statuses = tuple(simplify_residual(item)[1] for item in substitutions)
        except BaseException:
            # Failure to compute is a limitation, not permission to call a
            # candidate an equation-domain artifact.
            surviving.append(candidate)
            continue
        if "UNDEFINED" in statuses:
            excluded_artifacts.append(
                {
                    "candidate": candidate,
                    "undefined_substitution": substitutions,
                    "definedness": statuses,
                }
            )
        else:
            surviving.append(candidate)
    coverage = [
        (candidate, *_branch_is_covered(candidate, emitted)) for candidate in surviving
    ]
    missing = [
        candidate
        for candidate, status, _witness in coverage
        if status in {"NOT_COVERED", "NOT_COVERED_SAMPLED"}
    ]
    coverage_undecided = sum(
        status == "COVERAGE_UNDECIDED"
        for _candidate, status, _witness in coverage
    )
    coverage_sampled = 0
    if missing:
        verdict = "OMITTED_BRANCH"
    elif truncated or unresolved or coverage_undecided:
        verdict = "COMPLETENESS_UNDECIDED"
    else:
        verdict = "COMPLETE_FACTOR_COVER"
    return {
        "factor_covers": len(covers),
        "cover_truncated": truncated,
        "generated_candidate_count": len(candidates),
        "candidate_count": len(surviving),
        "artifact_count": artifact_count,
        "excluded_artifact_count": len(excluded_artifacts),
        "excluded_artifacts": excluded_artifacts,
        "unresolved_count": unresolved,
        "coverage_sampled": coverage_sampled,
        "coverage_undecided": coverage_undecided,
        "coverage": coverage,
        "missing": missing,
        "verdict": verdict,
    }


def exact_sample_assignments(symbols: Sequence[sp.Symbol]) -> Iterator[dict[sp.Symbol, sp.Rational]]:
    general_seed = (
        sp.Integer(-2),
        sp.Integer(-1),
        sp.Rational(-1, 2),
        sp.Rational(1, 2),
        sp.Integer(1),
        sp.Integer(2),
        sp.Integer(3),
    )
    symbols = tuple(sorted(symbols, key=sp.default_sort_key))
    if not symbols:
        yield {}
        return
    positive_record_names = {
        "bComp",
        "B_comp",
        "muR",
        "mu_R",
        "muBr",
        "mu_br",
        "rhoBr",
        "rho_br",
        "sRho",
        "s_rho",
        "cs0",
        "c_s0",
    }
    pools: dict[str, list[sp.Symbol]] = {
        "positive": [],
        "nonnegative": [],
        "general": [],
    }
    for symbol in symbols:
        if symbol.is_positive or str(symbol) in positive_record_names:
            pools["positive"].append(symbol)
        elif symbol.is_nonnegative:
            pools["nonnegative"].append(symbol)
        else:
            pools["general"].append(symbol)

    def value_pool(kind: str, count: int) -> tuple[sp.Rational, ...]:
        if kind == "positive":
            result = [sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(3)]
            next_positive = 4
            while len(result) < count:
                result.append(sp.Integer(next_positive))
                next_positive += 1
            return tuple(result)
        if kind == "nonnegative":
            result = [
                sp.Integer(0),
                sp.Rational(1, 2),
                sp.Integer(1),
                sp.Integer(2),
                sp.Integer(3),
            ]
            next_nonnegative = 4
            while len(result) < count:
                result.append(sp.Integer(next_nonnegative))
                next_nonnegative += 1
            return tuple(result)
        result = list(general_seed)
        magnitude = 3
        while len(result) < count:
            result.extend((sp.Integer(-magnitude), sp.Integer(magnitude + 1)))
            magnitude += 1
        return tuple(result)

    for serial in range(MAX_SAMPLE_ASSIGNMENTS):
        assignment: dict[sp.Symbol, sp.Rational] = {}
        for kind, pool_symbols in pools.items():
            values = value_pool(kind, len(pool_symbols))
            for pool_index, symbol in enumerate(pool_symbols):
                assignment[symbol] = values[(serial + pool_index) % len(values)]
            # This invariant is the sampler's genericity contract: no
            # coincidence locus between two symbols in one assumption pool is
            # visited by any emitted assignment.
            assigned = [assignment[symbol] for symbol in pool_symbols]
            if len(assigned) != len(set(assigned)):
                raise AssertionError(f"non-generic {kind} sample {serial}: {assigned}")
        yield assignment


def exact_truth(value: object) -> str:
    if value is True or value is sp.true:
        return "TRUE"
    if value is False or value is sp.false:
        return "FALSE"
    try:
        simplified = sp.simplify(value)
    except BaseException:
        simplified = value
    if simplified is True or simplified is sp.true:
        return "TRUE"
    if simplified is False or simplified is sp.false:
        return "FALSE"
    return "UNDECIDED"


def premise_conjuncts(value: object) -> tuple[object, ...]:
    """Flatten only the premise-level And/list structure into atomic conjuncts."""
    if isinstance(value, sp.And):
        items: Iterable[object] = value.args
    elif isinstance(value, (tuple, list, sp.Tuple)):
        items = value
    else:
        return (value,)
    if isinstance(value, (tuple, list, sp.Tuple, sp.And)):
        return tuple(
            conjunct
            for item in items
            for conjunct in premise_conjuncts(item)
        )
    raise AssertionError("unreachable premise structure")


def premise_conjunct_truth(value: object) -> str:
    """Classify one already-substituted premise atom at its principal value."""
    direct = exact_truth(value)
    if direct != "UNDECIDED":
        return direct

    # Applied assumptions deliberately survive assumption-free srepr parsing.
    # refine, unlike simplify, decides Q.real(1/2), Q.positive(1), and their
    # false concrete counterparts without importing assumptions for free
    # symbols that remain after witness substitution.
    try:
        refined = sp.refine(value)  # type: ignore[arg-type]
    except BaseException:
        refined = value
    refined_truth = exact_truth(refined)
    if refined_truth != "UNDECIDED":
        return refined_truth

    # Deferred WL Element atoms are Contains objects.  For a concrete
    # principal value, doit/simplify decides realness or integrality.  With
    # unbound field symbols they remain intact and therefore contingent.
    try:
        evaluated = value.doit() if isinstance(value, sp.Basic) else value
    except BaseException:
        evaluated = value
    return exact_truth(evaluated)


def conjunction(items: Iterable[object]) -> object:
    values = tuple(items)
    if not values:
        return sp.true
    return sp.And(*values)


def sample_existence(condition: object) -> dict[str, object]:
    symbols = tuple(getattr(condition, "free_symbols", ()))
    direct = exact_truth(condition)
    if direct == "TRUE":
        return {"truth": "TRUE", "witness": {}}
    if direct == "FALSE":
        return {"truth": "FALSE", "witness": {}}
    for assignment in exact_sample_assignments(symbols):
        try:
            value = condition.subs(assignment, simultaneous=True)  # type: ignore[union-attr]
        except BaseException:
            continue
        if exact_truth(value) == "TRUE":
            return {"truth": "TRUE", "witness": assignment}
    return {"truth": "UNDECIDED", "witness": {}}


def sample_counterexample(predicate: object, premises: object = sp.true) -> dict[str, object]:
    symbols = tuple(
        set(getattr(predicate, "free_symbols", ()))
        | set(getattr(premises, "free_symbols", ()))
    )
    for assignment in exact_sample_assignments(symbols):
        try:
            premise_value = premises.subs(assignment, simultaneous=True)  # type: ignore[union-attr]
            predicate_value = predicate.subs(assignment, simultaneous=True)  # type: ignore[union-attr]
        except BaseException:
            continue
        if exact_truth(premise_value) == "TRUE" and exact_truth(predicate_value) == "FALSE":
            return {"truth": "FALSE", "counterexample": assignment}
    direct = exact_truth(predicate)
    return {"truth": direct, "counterexample": {}}


def sequence(value: object) -> tuple[object, ...]:
    if isinstance(value, (tuple, list, sp.Tuple)):
        return tuple(value)
    return (value,)
