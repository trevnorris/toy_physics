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
from typing import Any, Iterable, Iterator, Mapping, Sequence

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica


# Floors required by brief obligation 7.  These apply independently to each
# semantic locus, witness, and undecided-record operation.
OPERATION_TIME_BUDGET_SECONDS = 60
OPERATION_MEMORY_BUDGET_BYTES = 512 * 1024 * 1024
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


def _replace_wl_elements(text: str) -> str:
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
            replacement = f"(Im[{expression}] == 0)"
        elif domain == "Integers":
            replacement = f"(Im[{expression}] == 0 && Sin[Pi*({expression})] == 0)"
        else:
            raise ValueError(f"unsupported Element domain: {domain}")
        source = source[:start] + replacement + source[end + 1 :]
        search_from = start + len(replacement)


def parse_wl(text: str) -> sp.Basic:
    # SymPy's Mathematica translator represents unknown Element[...] heads as
    # non-Boolean functions, which then cannot be children of And.  Translate
    # the two record domains exactly to Boolean relations before parsing.
    if text.strip() == "{}":
        return sp.Tuple()
    source = _replace_wl_elements(text)
    for original, protected in _WL_PROTECTED_SYMBOLS.items():
        source = re.sub(rf"\b{re.escape(original)}\b", protected, source)
    parsed = parse_mathematica(source)
    replacements = {
        sp.Symbol(protected): sp.Symbol(original)
        for original, protected in _WL_PROTECTED_SYMBOLS.items()
    }
    return parsed.xreplace(replacements)


def parse_py(text: str) -> sp.Basic:
    # srepr is executable constructor syntax, not user-authored source.  sympify
    # supplies the SymPy constructor namespace and evaluate=False preserves the
    # record's unevaluated relation forms where SymPy supports doing so.
    return sp.sympify(text, evaluate=False)


def parse_payload(dialect: str, text: str) -> sp.Basic:
    if dialect == "WL":
        return parse_wl(text)
    if dialect == "PY":
        return parse_py(text)
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
    }
    return powers | set(expression.atoms(sp.Abs))


def _masked_branch(
    branch: Mapping[sp.Expr, sp.Expr],
) -> tuple[dict[sp.Expr, sp.Expr], dict[sp.Symbol, sp.Basic]]:
    masked: dict[sp.Expr, sp.Expr] = {}
    restore: dict[sp.Symbol, sp.Basic] = {}
    serial = 0
    for lhs, rhs in branch.items():
        replacements: dict[sp.Basic, sp.Symbol] = {}
        for atom in sorted(_radical_atoms(rhs), key=sp.default_sort_key):
            dummy = sp.Dummy(f"branch_radical_{serial}")
            serial += 1
            replacements[atom] = dummy
            restore[dummy] = atom
        masked[lhs] = rhs.xreplace(replacements)
    return masked, restore


def sheet_residuals(
    residuals: Sequence[sp.Expr], branch: Mapping[sp.Expr, sp.Expr]
) -> tuple[list[tuple[str, tuple[sp.Expr, ...]]], int, int]:
    masked_branch, restore = _masked_branch(branch)
    substituted = tuple(
        residual.subs(masked_branch, simultaneous=True) for residual in residuals
    )
    sheet_atoms: set[sp.Basic] = set()
    for residual in substituted:
        sheet_atoms.update(_radical_atoms(residual))
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
        values = tuple(item.xreplace(replacements).xreplace(restore) for item in substituted)
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
            if simplified.is_zero is False or (not simplified.free_symbols and simplified != 0):
                return simplified, "NONZERO"
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
            if sampled == 0 or sampled.is_zero is True:
                sampled_zero += 1
                continue
            if not sampled.free_symbols and sampled != 0 and sampled.is_zero is not True:
                return simplified, "NONZERO_SAMPLED"
        if sampled_count and sampled_zero == sampled_count:
            return simplified, "ZERO_SAMPLED"
        return simplified, "UNDECIDED"
    except BaseException:
        return expression, "UNDECIDED"


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
        elif "UNDECIDED" in statuses:
            verdict = "BRANCH_MEMBERSHIP_UNDECIDED"
        elif "ZERO_SAMPLED" in statuses:
            verdict = "BRANCH_CONTAINED_SAMPLED"
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


def _branch_is_covered(
    candidate: Mapping[sp.Expr, sp.Expr], emitted: Sequence[Mapping[sp.Expr, sp.Expr]]
) -> bool:
    for branch in emitted:
        statuses: list[str] = []
        for lhs, rhs in branch.items():
            value = (lhs - rhs).subs(candidate, simultaneous=True)
            _simplified, status = simplify_residual(value)
            statuses.append(status)
        if statuses and all(status == "ZERO" for status in statuses):
            return True
        if not statuses and not candidate:
            return True
    return False


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
    missing = [candidate for candidate in candidates if not _branch_is_covered(candidate, emitted)]
    if missing:
        verdict = "OMITTED_BRANCH"
    elif truncated or unresolved:
        verdict = "COMPLETENESS_UNDECIDED"
    else:
        verdict = "COMPLETE_FACTOR_COVER"
    return {
        "factor_covers": len(covers),
        "cover_truncated": truncated,
        "candidate_count": len(candidates),
        "unresolved_count": unresolved,
        "missing": missing,
        "verdict": verdict,
    }


def exact_sample_assignments(symbols: Sequence[sp.Symbol]) -> Iterator[dict[sp.Symbol, sp.Rational]]:
    values = (
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
    for serial in range(MAX_SAMPLE_ASSIGNMENTS):
        assignment: dict[sp.Symbol, sp.Rational] = {}
        for index, symbol in enumerate(symbols):
            if symbol.is_positive or str(symbol) in positive_record_names:
                positive = values[3:]
                assignment[symbol] = positive[(serial + index) % len(positive)]
            elif symbol.is_nonnegative:
                nonnegative = (sp.Integer(0),) + values[3:]
                assignment[symbol] = nonnegative[(serial + index) % len(nonnegative)]
            else:
                assignment[symbol] = values[(serial + 2 * index) % len(values)]
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
