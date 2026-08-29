#!/usr/bin/env python3
"""Sort-routed S11c-b cross-engine symbolic adjudication.

The committed comparator supplies parsing, case axes, lossless WL jet decoding,
and compact-support divergence reduction.  The committed hand reconciliation
supplies the spelling map and bound-integral linearity.  This layer adds the
single cited bRho normalization and one cited container-label restoration,
then prints the resulting symbolic objects and complete case accounting.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Iterator, Sequence

import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational
from sympy.logic.boolalg import Boolean

import S11c_b_cross_engine_comparator as C
import S11c_b_handcoded_comparison as H


# This is deliberately the same object, not a copied or extended table.
WL_TO_PY_RENAME = H.WL_TO_PY_RENAME

BRHO = sp.Symbol("bRho")
B_RHO_3 = sp.Symbol("B_rho_3")
W_0 = sp.Symbol("W_0")

# Bridge A, and no other algebraic normalization:
#   WL theta^2 coefficient: bRho*anchoredWidth/2, anchoredWidth=W_bg
#     mathematica/S11c_b_brane_operator_mathematica_audit.wl:472-479
#   WL homogeneous coefficient: bRho*WZero/2
#     mathematica/S11c_b_brane_operator_mathematica_audit.wl:1621-1630
#   PY theta^2 coefficient: B_rho_3*W_bg/(2*W0)
#     scripts/S11c_b_brane_operator_sympy_audit.py:1130-1140
#   B_rho^(3) == B_rho*W_0
#     directives/S11c_b_SHARED_PHYSICS.md:102
BRIDGE_A_RULE = (BRHO, B_RHO_3 / W_0)
BRIDGE_A_SOURCES = (
    "mathematica/S11c_b_brane_operator_mathematica_audit.wl:472-479",
    "mathematica/S11c_b_brane_operator_mathematica_audit.wl:1621-1630",
    "scripts/S11c_b_brane_operator_sympy_audit.py:1130-1140",
    "directives/S11c_b_SHARED_PHYSICS.md:102",
)

# The registry is intentionally closed.  PY emits KINETIC as
# (U_MOMENTUM_ROWS, THICKNESS_ROW), while WL emits an Association carrying
# exactly those labels.  Both sources construct these two and only these two
# semantic blocks at the cited sites.
#   scripts/S11c_b_brane_operator_sympy_audit.py:1573-1580
#   mathematica/S11c_b_brane_operator_mathematica_audit.wl:851-853
KINETIC_ADAPTER_NAME = "SLAB_OPERATOR_TERM_ORIGINS/KINETIC_LABEL_RESTORATION"
KINETIC_ADAPTER_SOURCES = (
    "scripts/S11c_b_brane_operator_sympy_audit.py:1573-1580",
    "mathematica/S11c_b_brane_operator_mathematica_audit.wl:851-853",
)


@dataclass
class JoinedCase:
    family: str
    key: C.Key
    operand_a: object | None
    operand_b: object | None

    @property
    def rendered_key(self) -> str:
        return C.serialise_key(self.key)

    @property
    def identity(self) -> tuple[str, C.Key]:
        return self.family, self.key

    @property
    def rendered_identity(self) -> str:
        return f"{self.family} {self.rendered_key}"


@dataclass(frozen=True)
class CollapseRule:
    source: str
    image: str
    reconciled_source: str
    reconciled_image: sp.Expr


@dataclass
class RunState:
    emitted_ids: Counter[tuple[str, C.Key]]
    classified_ids: Counter[tuple[str, C.Key]]
    comparator_counts: Counter[str]
    route_counts: Counter[str]
    jet_counts: Counter[str]
    captured_names: Counter[str]
    touched_cases: list[str]


def _map_value(value: object, basic_map: Callable[[sp.Basic], sp.Basic]) -> object:
    if isinstance(value, C.Association):
        return C.Association(
            tuple((label, _map_value(item, basic_map)) for label, item in value.entries)
        )
    if isinstance(value, tuple):
        return tuple(_map_value(item, basic_map) for item in value)
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(basic_map)
    if isinstance(value, sp.Basic):
        return basic_map(value)
    return value


def _bridge_a(value: object) -> object:
    return _map_value(value, lambda item: item.xreplace({BRIDGE_A_RULE[0]: BRIDGE_A_RULE[1]}))


def _normal_name(name: str, names: dict[str, str]) -> str:
    if "XJETX" in name:
        base, suffix = name.split("XJETX", 1)
        base = names.get(base, base)
        return C.s11ca.canon_jet_name(base + "_" + suffix.replace("X", "_"))
    return C.s11ca.canon_jet_name(names.get(name, name))


def _collapse(value: object, rule: CollapseRule) -> object:
    def replace(item: sp.Basic) -> sp.Basic:
        substitutions = {
            symbol: rule.reconciled_image
            for symbol in item.atoms(sp.Symbol)
            if symbol.name == rule.reconciled_source
        }
        return item.xreplace(substitutions) if substitutions else item

    return _map_value(value, replace)


def transform(
    value: object | None,
    names: dict[str, str],
    *,
    bridge_a: bool,
    collapse: CollapseRule | None,
) -> object | None:
    if value is None:
        return None
    output = H.reconcile(value, names)
    if bridge_a:
        output = _bridge_a(output)
    if collapse is not None:
        output = _collapse(output, collapse)
    return output


def _iter_basics(value: object | None) -> Iterator[sp.Basic]:
    if value is None:
        return
    if isinstance(value, C.Association):
        for _label, item in value.entries:
            yield from _iter_basics(item)
    elif isinstance(value, tuple):
        for item in value:
            yield from _iter_basics(item)
    elif isinstance(value, sp.MatrixBase):
        for item in value:
            yield item
    elif isinstance(value, sp.Basic):
        yield value


def _token_names(value: object | None) -> Counter[str]:
    output: Counter[str] = Counter()
    for basic in _iter_basics(value):
        for node in sp.preorder_traversal(basic):
            if isinstance(node, sp.Symbol):
                output[node.name] += 1
            elif isinstance(node, AppliedUndef):
                output[node.func.__name__] += 1
    return output


def _jet_id(name: str, names: dict[str, str]) -> tuple[str, tuple[object, ...]] | None:
    canonical = _normal_name(name, names)
    theta = re.fullmatch(r"grad_theta_([123])", canonical)
    if theta is not None:
        return "theta", (int(theta.group(1)),)

    parts = canonical.split("_")
    derivative_parts: list[str] = []
    while parts and (
        re.fullmatch(r"t+", parts[-1])
        or re.fullmatch(r"(?:d[123])+", parts[-1])
        or re.fullmatch(r"(?:dw)+", parts[-1])
    ):
        derivative_parts.append(parts.pop())
    if not derivative_parts or not parts:
        return None

    time_order = 0
    directions: list[int] = []
    for part in reversed(derivative_parts):
        if set(part) == {"t"}:
            time_order += len(part)
        elif part.startswith("dw"):
            directions.extend([99] * (len(part) // 2))
        else:
            directions.extend(int(item) for item in re.findall(r"d([123])", part))
    multiindex: tuple[object, ...] = tuple(["t"] * time_order + sorted(directions))
    return "_".join(parts), multiindex


def _jet_multiset(
    values: Iterable[object | None], names: dict[str, str]
) -> Counter[tuple[str, tuple[object, ...]]]:
    output: Counter[tuple[str, tuple[object, ...]]] = Counter()
    for value in values:
        for token, count in _token_names(value).items():
            decoded = _jet_id(token, names)
            if decoded is not None:
                output[decoded] += count
    return output


def _jet_object(items: Counter[tuple[str, tuple[object, ...]]]) -> C.Association:
    entries = []
    for (base, multiindex), count in sorted(items.items(), key=repr):
        rendered_index = ",".join(str(item) for item in multiindex)
        entries.append((f"{base}|({rendered_index})", sp.Integer(count)))
    return C.Association(tuple(entries))


def _leaf_sort(value: object) -> str:
    # Positive Expr recognition is load-bearing: Symbol is admitted here even
    # in SymPy versions where it also inherits from Boolean.
    if isinstance(value, sp.Expr):
        return "EXPR"
    if isinstance(value, Relational):
        return f"RELATIONAL:{type(value).__name__}"
    if isinstance(value, Boolean):
        return f"BOOLEAN:{type(value).__name__}"
    if isinstance(value, C.TextAtom):
        return "TEXT_ATOM"
    if isinstance(value, sp.Basic):
        return f"BASIC:{type(value).__name__}"
    return type(value).__name__.upper()


def _shape_leaves(value: object) -> Iterator[object]:
    if isinstance(value, C.Association):
        for _label, item in value.entries:
            yield from _shape_leaves(item)
    elif isinstance(value, tuple):
        for item in value:
            yield from _shape_leaves(item)
    elif isinstance(value, sp.MatrixBase):
        yield from value
    else:
        yield value


def _shape(value: object | None) -> C.Association:
    if value is None:
        return C.Association((("KIND", C.TextAtom("MISSING")),))
    labels: tuple[str, ...] = ()
    dimensions: tuple[int, ...] = ()
    if isinstance(value, C.Association):
        kind = "ASSOCIATION"
        labels = tuple(label for label, _item in value.entries)
        arity = len(value.entries)
    elif isinstance(value, tuple):
        kind = "TUPLE"
        arity = len(value)
    elif isinstance(value, sp.MatrixBase):
        kind = "MATRIX"
        dimensions = tuple(value.shape)
        arity = len(value)
    else:
        kind = _leaf_sort(value)
        arity = 1
    sorts = Counter(_leaf_sort(item) for item in _shape_leaves(value))
    entries: list[tuple[str, object]] = [
        ("KIND", C.TextAtom(kind)),
        ("ARITY", sp.Integer(arity)),
        (
            "LEAF_SORT_COUNTS",
            C.Association(tuple((name, sp.Integer(count)) for name, count in sorted(sorts.items()))),
        ),
    ]
    if labels:
        entries.append(("LABELS", tuple(C.TextAtom(item) for item in labels)))
    if dimensions:
        entries.append(("DIMENSIONS", tuple(sp.Integer(item) for item in dimensions)))
    return C.Association(tuple(entries))


def _shape_difference(left: object | None, right: object | None) -> C.Association:
    entries: list[tuple[str, object]] = [("A", _shape(left)), ("B", _shape(right))]
    if isinstance(left, C.Association) and isinstance(right, C.Association):
        left_labels = {label for label, _item in left.entries}
        right_labels = {label for label, _item in right.entries}
        entries.extend(
            (
                ("MISSING_FROM_A", tuple(C.TextAtom(item) for item in sorted(right_labels - left_labels))),
                ("MISSING_FROM_B", tuple(C.TextAtom(item) for item in sorted(left_labels - right_labels))),
            )
        )
    return C.Association(tuple(entries))


def _arithmetic_shape(value: object) -> object | None:
    if isinstance(value, sp.Expr):
        return "EXPR"
    if isinstance(value, sp.MatrixBase):
        if all(isinstance(item, sp.Expr) for item in value):
            return ("MATRIX", value.shape)
        return None
    if isinstance(value, tuple):
        children = tuple(_arithmetic_shape(item) for item in value)
        if all(item is not None for item in children):
            return ("TUPLE", children)
    return None


def _arithmetic_pair(left: object, right: object) -> bool:
    left_shape = _arithmetic_shape(left)
    return left_shape is not None and left_shape == _arithmetic_shape(right)


def _expr_residual(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return H._residual_value(left, right)  # Bridge C is applied inside this cited helper.


def _arithmetic_residual(left: object, right: object) -> object:
    if isinstance(left, sp.Expr) and isinstance(right, sp.Expr):
        return _expr_residual(left, right)
    if isinstance(left, sp.MatrixBase) and isinstance(right, sp.MatrixBase):
        return sp.ImmutableMatrix(
            left.rows,
            left.cols,
            [_expr_residual(a, b) for a, b in zip(left, right)],
        )
    if isinstance(left, tuple) and isinstance(right, tuple):
        return tuple(_arithmetic_residual(a, b) for a, b in zip(left, right))
    raise C.InputError("arithmetic residual received a non-arithmetic shape")


def _zero_object(value: object) -> bool:
    if isinstance(value, C.Association):
        return all(_zero_object(item) for _label, item in value.entries)
    if isinstance(value, tuple):
        return all(_zero_object(item) for item in value)
    if isinstance(value, sp.MatrixBase):
        return all(item == 0 for item in value)
    return isinstance(value, sp.Basic) and value == 0


def _kinetic_pairs(
    family: str, key: C.Key, left: object, right: object
) -> tuple[tuple[str, sp.Expr, sp.Expr], ...] | None:
    axes = C.key_dict(key)
    if family != "SLAB_OPERATOR_TERM_ORIGINS" or axes.get("OBJECT") != "KINETIC":
        return None
    if not isinstance(left, tuple) or len(left) != 2 or not isinstance(right, C.Association):
        return None
    right_labels = tuple(label for label, _item in right.entries)
    if len(right_labels) != 2 or len(set(right_labels)) != 2:
        return None
    right_items = right.as_dict()
    if set(right_items) != {"U_MOMENTUM_ROWS", "THICKNESS_ROW"}:
        return None
    left_u, left_w = left
    right_u, right_w = right_items["U_MOMENTUM_ROWS"], right_items["THICKNESS_ROW"]
    if not (
        isinstance(left_u, tuple)
        and isinstance(right_u, tuple)
        and len(left_u) == len(right_u) == 3
        and all(isinstance(item, sp.Expr) for item in (*left_u, *right_u, left_w, right_w))
    ):
        return None
    return tuple(
        (f"U_MOMENTUM_ROWS[{index}]", left_u[index], right_u[index])
        for index in range(3)
    ) + (("THICKNESS_ROW", left_w, right_w),)


def _kinetic_residual(pairs: tuple[tuple[str, sp.Expr, sp.Expr], ...]) -> C.Association:
    return C.Association(tuple((label, _expr_residual(left, right)) for label, left, right in pairs))


def _coverage(value: object | None) -> bool:
    return value is None or (
        isinstance(value, C.TextAtom) and value.value == "UNDEFINED_UNJOINED"
    )


def _comparison_payload(family: str, key: C.Key, left: object | None, right: object | None) -> object:
    if _coverage(left) or _coverage(right):
        return C.TextAtom("UNDEFINED_UNJOINED")
    if left is None or right is None:
        raise C.InputError("coverage routing did not retain its missing operand")
    if _arithmetic_pair(left, right):
        return _arithmetic_residual(left, right)
    pairs = _kinetic_pairs(family, key, left, right)
    if pairs is not None:
        return _kinetic_residual(pairs)
    return C.TextAtom("UNDEFINED_STRUCTURE_INCOMPLETE")


def _render_operand(value: object | None) -> str:
    return "<MISSING>" if value is None else C.serialise(value)


def _same_syntax(left: object | None, right: object | None) -> bool:
    return _render_operand(left) == _render_operand(right)


def _emit_case_objects(case: JoinedCase, left: object | None, right: object | None, residual: object) -> None:
    print(f"operand_A {case.family} {case.rendered_key} = {_render_operand(left)}", flush=True)
    print(f"operand_B {case.family} {case.rendered_key} = {_render_operand(right)}", flush=True)
    print(f"A_minus_B {case.family} {case.rendered_key} = {C.serialise(residual)}", flush=True)


def _emit_jet_line(
    case: JoinedCase,
    raw_values: tuple[object | None, object | None],
    transformed_values: tuple[object | None, object | None],
    names: dict[str, str],
    state: RunState,
) -> None:
    before = _jet_multiset(raw_values, names)
    after = _jet_multiset(transformed_values, names)
    label = "JET_CONSERVED" if before == after else "JET_LOST"
    state.jet_counts[label] += 1
    print(
        f"{label} {case.family} {case.rendered_key} "
        f"before={C.serialise(_jet_object(before))} after={C.serialise(_jet_object(after))}",
        flush=True,
    )


def _emit_ablation_case(
    case: JoinedCase,
    before_values: tuple[object | None, object | None],
    after_values: tuple[object | None, object | None],
    names: dict[str, str],
) -> None:
    before_residual = _comparison_payload(case.family, case.key, *before_values)
    after_residual = _comparison_payload(case.family, case.key, *after_values)
    before_jets = _jet_multiset(before_values, names)
    after_jets = _jet_multiset(after_values, names)
    print(f"ABLATION_CASE {case.family} {case.rendered_key}", flush=True)
    print(f"transformed_before_A = {_render_operand(before_values[0])}", flush=True)
    print(f"transformed_after_A = {_render_operand(after_values[0])}", flush=True)
    print(f"transformed_before_B = {_render_operand(before_values[1])}", flush=True)
    print(f"transformed_after_B = {_render_operand(after_values[1])}", flush=True)
    print(f"residual_before = {C.serialise(before_residual)}", flush=True)
    print(f"residual_after = {C.serialise(after_residual)}", flush=True)
    print(f"jet_multiset_before = {C.serialise(_jet_object(before_jets))}", flush=True)
    print(f"jet_multiset_after = {C.serialise(_jet_object(after_jets))}", flush=True)


def _classify_case(
    case: JoinedCase,
    left: object | None,
    right: object | None,
    state: RunState,
) -> None:
    if _coverage(left) or _coverage(right):
        residual: object = C.TextAtom("UNDEFINED_UNJOINED")
        state.route_counts["COVERAGE"] += 1
        print(f"COVERAGE {case.family} {case.rendered_key}", flush=True)
        _emit_case_objects(case, left, right, residual)
        return

    if left is None or right is None:
        raise C.InputError("coverage routing did not retain its missing operand")
    if _arithmetic_pair(left, right):
        residual = _arithmetic_residual(left, right)
        label = "MATCH" if _zero_object(residual) else "FLAG"
        state.route_counts[f"ALGEBRAIC_{label}"] += 1
        print(f"ALGEBRAIC {label} {case.family} {case.rendered_key}", flush=True)
        _emit_case_objects(case, left, right, residual)
        return

    pairs = _kinetic_pairs(case.family, case.key, left, right)
    if pairs is not None:
        residual = _kinetic_residual(pairs)
        label = "MATCH" if _zero_object(residual) else "FLAG"
        state.route_counts[f"CONTAINER_{label}"] += 1
        print(f"CONTAINER {label} {case.family} {case.rendered_key}", flush=True)
        print(f"container_adapter = {C.serialise(C.TextAtom(KINETIC_ADAPTER_NAME))}", flush=True)
        print(
            "container_adapter_sources = "
            + C.serialise(tuple(C.TextAtom(item) for item in KINETIC_ADAPTER_SOURCES)),
            flush=True,
        )
        _emit_case_objects(case, left, right, residual)
        for label_name, item in residual.entries:
            print(
                f"leaf_A_minus_B {case.family} {case.rendered_key} "
                f"label={label_name!r} = {C.serialise(item)}",
                flush=True,
            )
        return

    residual = C.TextAtom("UNDEFINED_STRUCTURE_INCOMPLETE")
    state.route_counts["STRUCTURE_INCOMPLETE"] += 1
    print(f"STRUCTURE_INCOMPLETE {case.family} {case.rendered_key}", flush=True)
    print(
        f"shape_diff {case.family} {case.rendered_key} = "
        f"{C.serialise(_shape_difference(left, right))}",
        flush=True,
    )
    _emit_case_objects(case, left, right, residual)


def _materialized(case: C.ParsedCase | None) -> object | None:
    if case is None:
        return None
    introduced_error = C.materialize(case)
    if introduced_error or case.error is not None:
        raise C.InputError(case.error or "case materialization error")
    return case.compare_value


def _family_cases(
    family: str, py_tags: dict[str, str], wl_tags: dict[str, str], state: RunState
) -> Iterator[JoinedCase]:
    extracted = C.extract_family(family, py_tags.get(family), wl_tags.get(family))
    py_by_key: dict[C.Key, list[C.ParsedCase]] = defaultdict(list)
    wl_by_key: dict[C.Key, list[C.ParsedCase]] = defaultdict(list)
    for item in extracted.py:
        py_by_key[item.key].append(item)
    for item in extracted.wl:
        wl_by_key[item.key].append(item)

    duplicate_keys = {
        key for key, values in py_by_key.items() if len(values) != 1
    } | {key for key, values in wl_by_key.items() if len(values) != 1}
    if duplicate_keys:
        rendered = tuple(C.serialise_key(key) for key in sorted(duplicate_keys, key=C.serialise_key))
        raise C.InputError(f"{family}: duplicate case axes {rendered}")

    py_keys, wl_keys = set(py_by_key), set(wl_by_key)
    common = py_keys & wl_keys
    state.comparator_counts["join"] += len(common)
    state.comparator_counts["py_only"] += len(py_keys - wl_keys)
    state.comparator_counts["wl_only"] += len(wl_keys - py_keys)
    for key in sorted(py_keys | wl_keys, key=C.serialise_key):
        py_item = py_by_key.get(key, [None])[0]
        wl_item = wl_by_key.get(key, [None])[0]
        left = _materialized(py_item)
        right = _materialized(wl_item)
        joined = JoinedCase(family, key, left, right)
        state.emitted_ids[joined.identity] += 1
        yield joined
        C.release_case(py_item)
        C.release_case(wl_item)


def _zero_order_base(name: str, names: dict[str, str]) -> str:
    canonical = _normal_name(name, names)
    decoded = _jet_id(canonical, {})
    return decoded[0] if decoded is not None else canonical


def _collapse_rule(raw: str) -> CollapseRule:
    if "=" not in raw:
        raise C.InputError("--collapse-jet requires <token>=<base>")
    source, image = raw.split("=", 1)
    if not source or not image:
        raise C.InputError("--collapse-jet requires nonempty names")
    source_id = _jet_id(source, dict(WL_TO_PY_RENAME))
    if source_id is None or not source_id[1]:
        raise C.InputError("--collapse-jet source is not a derivative jet token")
    reconciled_source = _normal_name(source, dict(WL_TO_PY_RENAME))
    if image == "0":
        reconciled_image: sp.Expr = sp.S.Zero
    else:
        image_id = _jet_id(image, dict(WL_TO_PY_RENAME))
        image_base = image_id[0] if image_id is not None else _zero_order_base(image, dict(WL_TO_PY_RENAME))
        image_order = len(image_id[1]) if image_id is not None else 0
        if image_base != source_id[0] or image_order >= len(source_id[1]):
            raise C.InputError("--collapse-jet image must have the same base and lower derivative order")
        reconciled_image = sp.Symbol(_normal_name(image, dict(WL_TO_PY_RENAME)))
    return CollapseRule(source, image, reconciled_source, reconciled_image)


def _configuration_objects(active_names: dict[str, str], bridge_a: bool, collapse: CollapseRule | None) -> None:
    print(f"RENAME_MAP_SIZE = {C.serialise(sp.Integer(len(active_names)))}", flush=True)
    bridge_object: object = (
        sp.Eq(BRIDGE_A_RULE[0], BRIDGE_A_RULE[1], evaluate=False)
        if bridge_a
        else C.TextAtom("DROPPED")
    )
    print(f"BRIDGE_A_SUBSTITUTION = {C.serialise(bridge_object)}", flush=True)
    print(
        "BRIDGE_A_SOURCES = "
        + C.serialise(tuple(C.TextAtom(item) for item in BRIDGE_A_SOURCES)),
        flush=True,
    )
    print(
        "CONTAINER_ADAPTER_REGISTRY = "
        + C.serialise(
            C.Association(
                (
                    (
                        KINETIC_ADAPTER_NAME,
                        tuple(C.TextAtom(item) for item in KINETIC_ADAPTER_SOURCES),
                    ),
                )
            )
        ),
        flush=True,
    )
    collapse_object: object = (
        C.TextAtom("NONE")
        if collapse is None
        else sp.Eq(sp.Symbol(collapse.reconciled_source), collapse.reconciled_image, evaluate=False)
    )
    print(f"COLLAPSE_JET_SUBSTITUTION = {C.serialise(collapse_object)}", flush=True)


def run(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--py", type=Path, default=C.DEFAULT_PY)
    parser.add_argument("--wl", type=Path, default=C.DEFAULT_WL)
    parser.add_argument(
        "--residual-leaf-budget",
        type=float,
        default=60.0,
        help="reported generous comparator leaf budget; routed arithmetic uses deterministic reconciliation",
    )
    hooks = parser.add_mutually_exclusive_group()
    hooks.add_argument("--collapse-jet", metavar="token=base")
    hooks.add_argument("--drop-bridge-a", action="store_true")
    hooks.add_argument("--drop-rename", metavar="WLname")
    arguments = parser.parse_args(argv)

    started = time.monotonic()
    try:
        if arguments.residual_leaf_budget <= 0:
            raise C.InputError("--residual-leaf-budget must be positive")
        if arguments.drop_rename is not None and arguments.drop_rename not in WL_TO_PY_RENAME:
            raise C.InputError(f"unknown spelling equivalence {arguments.drop_rename!r}")

        collapse = _collapse_rule(arguments.collapse_jet) if arguments.collapse_jet else None
        active_names = dict(WL_TO_PY_RENAME)
        if arguments.drop_rename is not None:
            active_names.pop(arguments.drop_rename)
            H._disable_imported_prepass_for_drop(arguments.drop_rename)

        py_tags = C.load_py(arguments.py)
        wl_tags = C.load_wl(arguments.wl)
        state = RunState(Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), [])
        bridge_a = not arguments.drop_bridge_a
        hook_active = collapse is not None or arguments.drop_bridge_a or arguments.drop_rename is not None

        _configuration_objects(active_names, bridge_a, collapse)
        print(
            f"RESIDUAL_LEAF_BUDGET_SECONDS = {C.serialise(sp.Float(arguments.residual_leaf_budget))}",
            flush=True,
        )

        for family in H.CORE_FAMILIES:
            for case in _family_cases(family, py_tags, wl_tags, state):
                raw_values = (case.operand_a, case.operand_b)
                state.captured_names.update(_token_names(case.operand_a))
                state.captured_names.update(_token_names(case.operand_b))

                transformed_values = (
                    transform(case.operand_a, active_names, bridge_a=bridge_a, collapse=collapse),
                    transform(case.operand_b, active_names, bridge_a=bridge_a, collapse=collapse),
                )
                if bridge_a and collapse is None and arguments.drop_rename is None:
                    for value in transformed_values:
                        if value is not None and any(BRHO in basic.atoms(sp.Symbol) for basic in _iter_basics(value)):
                            if _arithmetic_shape(value) is not None:
                                raise C.InputError(f"residual bRho in arithmetic operand {case.rendered_identity}")

                if hook_active:
                    baseline_values = (
                        transform(case.operand_a, dict(WL_TO_PY_RENAME), bridge_a=True, collapse=None),
                        transform(case.operand_b, dict(WL_TO_PY_RENAME), bridge_a=True, collapse=None),
                    )
                    if not (
                        _same_syntax(baseline_values[0], transformed_values[0])
                        and _same_syntax(baseline_values[1], transformed_values[1])
                    ):
                        state.touched_cases.append(case.rendered_identity)
                        _emit_ablation_case(case, baseline_values, transformed_values, active_names)

                _classify_case(case, *transformed_values, state)
                state.classified_ids[case.identity] += 1
                _emit_jet_line(case, raw_values, transformed_values, active_names, state)

        if collapse is not None and state.captured_names[collapse.source] == 0:
            raise C.InputError(f"non-occurring collapse token {collapse.source!r}")
        if arguments.drop_rename is not None and state.captured_names[arguments.drop_rename] == 0:
            raise C.InputError(f"non-occurring spelling name {arguments.drop_rename!r}")
        if hook_active and not state.touched_cases:
            raise C.InputError("ablation produced an empty touched-case set")

        if state.classified_ids != state.emitted_ids:
            missing = state.emitted_ids - state.classified_ids
            extra = state.classified_ids - state.emitted_ids
            raise C.InputError(f"case-ID multiset difference missing={missing} extra={extra}")

        for family in H.NAMESPACE_INCOMPLETE_FAMILIES:
            print(
                f"NAMESPACE_INCOMPLETE {family} "
                "(WL operand unparsed; cross-engine control comparison owed; "
                "each engine's internal control verified in the build legs)",
                flush=True,
            )

        if hook_active:
            print(
                "ABLATION_TOUCHED_CASES = "
                + C.serialise(tuple(C.TextAtom(item) for item in state.touched_cases)),
                flush=True,
            )
        print(
            "COMPARATOR_CORE_ACCOUNTING = "
            + C.serialise(
                C.Association(
                    tuple((name, sp.Integer(state.comparator_counts[name])) for name in ("join", "py_only", "wl_only"))
                )
            ),
            flush=True,
        )
        print(
            "ROUTE_ACCOUNTING = "
            + C.serialise(
                C.Association(
                    tuple(
                        (name, sp.Integer(state.route_counts[name]))
                        for name in (
                            "ALGEBRAIC_MATCH",
                            "ALGEBRAIC_FLAG",
                            "CONTAINER_MATCH",
                            "CONTAINER_FLAG",
                            "STRUCTURE_INCOMPLETE",
                            "COVERAGE",
                        )
                    )
                )
            ),
            flush=True,
        )
        emitted_total = sum(state.emitted_ids.values())
        classified_total = sum(state.classified_ids.values())
        print(f"EMITTED_CORE_CASES = {C.serialise(sp.Integer(emitted_total))}", flush=True)
        print(f"CLASSIFIED_CORE_CASES = {C.serialise(sp.Integer(classified_total))}", flush=True)
        print(f"CASE_ID_MULTISET_EQUAL = {C.serialise(sp.true)}", flush=True)
        print(
            "JET_ACCOUNTING = "
            + C.serialise(
                C.Association(
                    tuple((name, sp.Integer(state.jet_counts[name])) for name in ("JET_CONSERVED", "JET_LOST"))
                )
            ),
            flush=True,
        )
        print(f"RUNTIME_SECONDS = {C.serialise(sp.Float(time.monotonic() - started, 8))}", flush=True)
        return 0
    except Exception as error:
        print(f"OPERATIONAL_ERROR {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        return 2


if __name__ == "__main__":
    raise SystemExit(run())
