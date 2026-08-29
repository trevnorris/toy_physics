#!/usr/bin/env python3
"""Strong/weak-aware S11c-b cross-engine symbolic adjudication.

The committed comparator supplies parsing, case axes, and lossless WL jet
decoding.  The committed hand reconciliation supplies the spelling map and
bound-integral linearity.  This layer consumes raw case values, expands held
divergences in the anchored jet algebra, applies the engine's imported profile
grade map last, and constructs certificates only for source-proven weak scalar
densities.  It prints computed objects and residuals; interpretation is kept
outside this measurement instrument.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, Sequence

import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational
from sympy.logic.boolalg import Boolean

import S11c_b_cross_engine_comparator as C
import S11c_b_brane_operator_sympy_audit as P
import S11c_b_handcoded_comparison as H


# This is deliberately the same object, not a copied or extended table.
WL_TO_PY_RENAME = H.WL_TO_PY_RENAME

BRHO = sp.Symbol("bRho")
B_RHO_3 = sp.Symbol("B_rho_3")
W_0 = sp.Symbol("W_0")

# Bridge D is the engine object itself.  Do not copy, extend, or re-derive it.
PROFILE_GRADE_SUBS = P.PROFILE_GRADE_SUBS
PROFILE_GRADE_SYMBOLS_BY_NAME = {
    symbol.name: symbol
    for expression in (*PROFILE_GRADE_SUBS.keys(), *PROFILE_GRADE_SUBS.values())
    for symbol in expression.free_symbols
}
BRIDGE_D_SOURCES = (
    "scripts/S11c_b_brane_operator_sympy_audit.py:662-670",
    "scripts/S11c_b_brane_operator_sympy_audit.py:714",
    "scripts/S11c_b_brane_operator_sympy_audit.py:1850-1869",
)

PROTECTED_ATOM_NAMES = frozenset(
    {
        "gamma_s11cb_w_bg_07",
        "gamma_s11cb_w_bg_10",
        "gamma_s11cb_mu_r_bg_07",
        "gamma_s11cb_mu_r_bg_10",
        "gammaWidthDivGradTheta",
        "gammaWidthDivGradEw",
        "gammaModulusDivGradTheta",
        "gammaModulusDivGradEw",
    }
)

ENERGY_STRUCTURE_FAMILIES = frozenset(
    {
        "ENERGY_BASIS_VARIABLE",
        "ENERGY_BASIS_NEW_INVARIANTS",
        "ENERGY_BASIS_OMISSIONS",
    }
)
ENERGY_EXACT_FAMILY = "ENERGY_BASIS_COUNT"
WEAK_ORIGIN_FAMILY = "COUPLING_KERNEL_TERM_ORIGINS"
WEAK_DENSITY_SOURCES = (
    "mathematica/S11c_b_brane_operator_mathematica_audit.wl:1075-1124",
    "scripts/S11c_b_brane_operator_sympy_audit.py:1889-1996",
)
WEAK_ORIGIN_SOURCES = (
    "mathematica/S11c_b_brane_operator_mathematica_audit.wl:1193",
    "scripts/S11c_b_brane_operator_sympy_audit.py:1996",
)

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


@dataclass(frozen=True)
class FieldRegistry:
    """A closed formal spatial-jet algebra.

    ``derivative_maps`` uses directions 1,2,3.  ``variational_fields`` selects
    fields used by the Euler screen and homotopy witness; other registered
    fields remain differentiable coefficients.
    """

    field_bases: frozenset[str]
    variational_fields: frozenset[str]
    derivative_maps: Mapping[int, Mapping[sp.Symbol, sp.Expr]]


@dataclass(frozen=True)
class DivergenceClassification:
    route: str
    euler_signature: C.Association
    vector: tuple[sp.Expr, sp.Expr, sp.Expr] | None
    anchored_remainder: sp.Expr | None
    bridge_d_certificate: sp.Expr | None


PRODUCTION_FIELD_BASES = frozenset(
    {
        *C.DIFFERENTIABLE_BASES,
        "u_1",
        "u_2",
        "u_3",
        "theta",
        "e_W",
    }
)


def _production_derivative_maps() -> dict[int, dict[sp.Symbol, sp.Expr]]:
    """Reference the engine derivative map and add only generic profile jets.

    The engine map is keyed 0,1,2.  The comparator's declared profile fields
    are added for witness differentiation; the imported entries remain the
    authoritative mappings wherever present.
    """

    output = {direction + 1: dict(P.DERIVATIVE_MAP[direction]) for direction in range(3)}
    for base in C.DIFFERENTIABLE_BASES:
        for direction in range(1, 4):
            source = sp.Symbol(base)
            output[direction].setdefault(source, sp.Symbol(f"{base}_d{direction}"))
    return output


PRODUCTION_FIELD_REGISTRY = FieldRegistry(
    PRODUCTION_FIELD_BASES,
    frozenset(C.DEPENDENT_BASES),
    _production_derivative_maps(),
)


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


def _bridge_d(value: object) -> object:
    def apply(item: sp.Basic) -> sp.Basic:
        assumption_normalisation = {
            symbol: PROFILE_GRADE_SYMBOLS_BY_NAME[symbol.name]
            for symbol in item.free_symbols
            if symbol.name in PROFILE_GRADE_SYMBOLS_BY_NAME
            and symbol != PROFILE_GRADE_SYMBOLS_BY_NAME[symbol.name]
        }
        normalised = item.xreplace(assumption_normalisation)
        return normalised.subs(PROFILE_GRADE_SUBS, simultaneous=True)

    return _map_value(
        value,
        apply,
    )


def _decode_registry_jet(
    symbol: sp.Symbol, registry: FieldRegistry
) -> tuple[str, tuple[int, ...]] | None:
    name = symbol.name
    for base in sorted(registry.field_bases, key=len, reverse=True):
        if name == base:
            return base, ()
        prefix = base + "_"
        if not name.startswith(prefix):
            continue
        suffix = name[len(prefix) :]
        if re.fullmatch(r"d[123](?:d[123])*", suffix):
            return base, tuple(int(item) for item in re.findall(r"d([123])", suffix))
    return None


def _registry_jet_symbol(base: str, directions: Sequence[int]) -> sp.Symbol:
    suffix = "" if not directions else "_" + "".join(
        f"d{item}" for item in sorted(directions)
    )
    return sp.Symbol(base + suffix)


def formal_dx(expression: sp.Expr, direction: int, registry: FieldRegistry) -> sp.Expr:
    """Apply the engine-style total derivative in a declared jet algebra."""

    if direction not in (1, 2, 3):
        raise C.InputError(f"invalid formal derivative direction {direction}")
    result = sp.S.Zero
    exact_map = registry.derivative_maps[direction]
    for atom in expression.free_symbols:
        derivative = exact_map.get(atom)
        if derivative is None:
            decoded = _decode_registry_jet(atom, registry)
            if decoded is None:
                continue
            base, directions = decoded
            derivative = _registry_jet_symbol(base, (*directions, direction))
        partial = sp.diff(expression, atom)
        if partial != 0:
            result += partial * derivative
    return sp.expand(result)


def formal_divergence(
    vector: Sequence[sp.Expr], registry: FieldRegistry
) -> sp.Expr:
    if len(vector) != 3:
        raise C.InputError("formal divergence requires a three-component vector")
    return sp.expand(
        sp.Add(*(formal_dx(vector[index], index + 1, registry) for index in range(3)))
    )


def _held_divergence(node: AppliedUndef, registry: FieldRegistry) -> sp.Expr:
    if len(node.args) != 2:
        raise C.InputError("HeldDiv does not have vector and coordinate arguments")
    vector, coordinates = node.args
    if not isinstance(vector, (tuple, sp.Tuple)) or len(vector) != 3:
        raise C.InputError("HeldDiv vector is not three-dimensional")
    if not isinstance(coordinates, (tuple, sp.Tuple)) or len(coordinates) != 3:
        raise C.InputError("HeldDiv coordinates are not three-dimensional")
    if not all(isinstance(item, sp.Expr) for item in vector):
        raise C.InputError("HeldDiv contains a non-expression component")
    return formal_divergence(tuple(vector), registry)


def _expand_held_divergences(value: object, registry: FieldRegistry) -> object:
    def replace(item: sp.Basic) -> sp.Basic:
        output = item
        while True:
            nodes = tuple(
                node
                for node in sp.preorder_traversal(output)
                if isinstance(node, AppliedUndef) and node.func.__name__ == "HeldDiv"
            )
            if not nodes:
                return output
            innermost = tuple(
                node
                for node in nodes
                if not any(
                    isinstance(child, AppliedUndef)
                    and child.func.__name__ == "HeldDiv"
                    for argument in node.args
                    for child in sp.preorder_traversal(argument)
                )
            )
            if not innermost:
                raise C.InputError("HeldDiv nesting could not be ordered")
            substitutions = {
                node: _held_divergence(node, registry) for node in innermost
            }
            output = output.xreplace(substitutions)

    return _map_value(value, replace)


def _contains_held_divergence(value: object | None) -> bool:
    return any(
        isinstance(node, AppliedUndef) and node.func.__name__ == "HeldDiv"
        for basic in _iter_basics(value)
        for node in sp.preorder_traversal(basic)
    )


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
    bridge_d: bool,
    collapse: CollapseRule | None,
) -> object | None:
    if value is None:
        return None
    output = H.reconcile(value, names)
    if bridge_a:
        output = _bridge_a(output)
    if collapse is not None:
        output = _collapse(output, collapse)
    output = _expand_held_divergences(output, PRODUCTION_FIELD_REGISTRY)
    if bridge_d:
        output = _bridge_d(output)
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
    base = "_".join(parts)
    base = {"W_bg": "w1_profile", "mu_R_bg": "m1_profile"}.get(base, base)
    return base, multiindex


def _bridge_d_order_preserving_jets() -> frozenset[tuple[str, tuple[object, ...]]]:
    output: set[tuple[str, tuple[object, ...]]] = set()
    for source, image in PROFILE_GRADE_SUBS.items():
        source_id = _jet_id(source.name, {})
        if source_id is None:
            continue
        image_ids = {
            decoded
            for symbol in image.free_symbols
            if (decoded := _jet_id(symbol.name, {})) is not None
        }
        if source_id in image_ids:
            output.add(source_id)
    return frozenset(output)


BRIDGE_D_ORDER_PRESERVING_JETS = _bridge_d_order_preserving_jets()


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


def _normalise_exact(expression: sp.Expr) -> sp.Expr:
    expanded = sp.expand(expression)
    if expanded == 0:
        return sp.S.Zero
    # Large production residuals are already in the comparator's canonical
    # expanded rational form.  Global polynomial GCD cancellation can take
    # minutes while adding no cancellation for those expressions; reserve it
    # for the compact fixture/witness certificates that need the fallback.
    if sp.count_ops(expanded, visual=False) > 2000:
        return expanded
    try:
        cancelled = sp.cancel(expanded)
    except Exception:
        return expanded
    return sp.S.Zero if cancelled == 0 else cancelled


def _euler_signature(
    residual: sp.Expr, registry: FieldRegistry
) -> C.Association:
    entries: list[tuple[str, object]] = []
    decoded_symbols = {
        symbol: decoded
        for symbol in residual.free_symbols
        if (decoded := _decode_registry_jet(symbol, registry)) is not None
        and decoded[0] in registry.variational_fields
    }
    bases_present = {
        decoded[0] for decoded in decoded_symbols.values()
    }
    bases_by_cost = sorted(
        bases_present,
        key=lambda base: (
            sum(1 for decoded in decoded_symbols.values() if decoded[0] == base),
            base,
        ),
    )
    for base in bases_by_cost:
        euler = sp.S.Zero
        for symbol, (symbol_base, directions) in decoded_symbols.items():
            if symbol_base != base:
                continue
            term = sp.diff(residual, symbol)
            for direction in directions:
                term = -formal_dx(term, direction, registry)
            euler += term
        euler = _normalise_exact(euler)
        if euler != 0:
            entries.append((base, euler))
            # A single exact nonzero variational component is sufficient to
            # rule out total divergence.  Remaining components are not needed
            # for this classifier's route or printed certificate.
            break
    return C.Association(tuple(entries))


def _integrate_scaling_polynomial(expression: sp.Expr, scale: sp.Dummy) -> sp.Expr:
    result = sp.S.Zero
    for term in sp.Add.make_args(sp.expand(expression)):
        power = term.as_powers_dict().get(scale, sp.S.Zero)
        if not power.is_Integer or power < 0:
            raise C.InputError("homotopy scaling integrand is not polynomial")
        coefficient = term / scale**power
        if scale in coefficient.free_symbols:
            raise C.InputError("homotopy scaling coefficient retains the scale symbol")
        result += coefficient / (power + 1)
    return sp.expand(result)


def _homotopy_vector(
    residual: sp.Expr, registry: FieldRegistry
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Construct a divergence witness via the scaling/Lagrange identity."""

    scale = sp.Dummy("lambda_s11cb", positive=True)
    scaled_symbols: dict[sp.Symbol, sp.Expr] = {}
    decoded_symbols: dict[sp.Symbol, tuple[str, tuple[int, ...]]] = {}
    for symbol in residual.free_symbols:
        decoded = _decode_registry_jet(symbol, registry)
        if decoded is None or decoded[0] not in registry.variational_fields:
            continue
        scaled_symbols[symbol] = scale * symbol
        decoded_symbols[symbol] = decoded
    if not decoded_symbols:
        raise C.InputError("residual contains no registered variational field")

    scaled = residual.xreplace(scaled_symbols)
    lagrange_vector = [sp.S.Zero, sp.S.Zero, sp.S.Zero]
    for symbol, (base, directions) in sorted(
        decoded_symbols.items(), key=lambda item: item[0].name
    ):
        if not directions:
            continue
        coefficient = sp.diff(scaled, symbol)
        current_directions = list(directions)
        sign = sp.S.One
        while current_directions:
            direction = current_directions.pop()
            lower = _registry_jet_symbol(base, current_directions)
            lagrange_vector[direction - 1] += sign * lower * coefficient
            coefficient = formal_dx(coefficient, direction, registry)
            sign = -sign

    return tuple(
        _normalise_exact(_integrate_scaling_polynomial(component / scale, scale))
        for component in lagrange_vector
    )  # type: ignore[return-value]


def classify_total_divergence(
    residual_pre_bridge_d: sp.Expr,
    registry: FieldRegistry,
    *,
    apply_bridge_d: bool,
) -> DivergenceClassification:
    """Classify one raw scalar residual and, when possible, certify a vector."""

    residual_pre_bridge_d = _normalise_exact(residual_pre_bridge_d)
    signature = _euler_signature(residual_pre_bridge_d, registry)
    if signature.entries:
        return DivergenceClassification(
            "RESIDUAL_BULK", signature, None, None, None
        )
    try:
        vector = _homotopy_vector(residual_pre_bridge_d, registry)
        anchored_remainder = _normalise_exact(
            residual_pre_bridge_d - formal_divergence(vector, registry)
        )
        certificate = (
            _normalise_exact(_bridge_d(anchored_remainder))
            if apply_bridge_d
            else anchored_remainder
        )
    except Exception:
        return DivergenceClassification(
            "DIVERGENCE_INCOMPLETE", signature, None, None, None
        )
    if certificate != 0:
        return DivergenceClassification(
            "DIVERGENCE_INCOMPLETE",
            signature,
            vector,
            anchored_remainder,
            certificate,
        )
    return DivergenceClassification(
        "REPRESENTATIONAL_DIVERGENCE",
        signature,
        vector,
        anchored_remainder,
        certificate,
    )


def _free_symbol_names(value: object) -> frozenset[str]:
    return frozenset(
        symbol.name
        for basic in _iter_basics(value)
        for symbol in basic.free_symbols
    )


def _is_weak_scalar_density(case: JoinedCase, residual: object) -> bool:
    if case.family != "COUPLING_KERNEL" or not isinstance(residual, sp.Expr):
        return False
    return C.key_dict(case.key).get("OBJECT") != "ADJOINTNESS_RELATION"


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
    *,
    allow_profile_cancellation: bool,
) -> None:
    before = _jet_multiset(raw_values, names)
    after = _jet_multiset(transformed_values, names)
    # Bridge D can algebraically cancel duplicate occurrences while its
    # definitional images retain every derivative multi-index.  Conservation
    # therefore tracks ordered jet support; the full occurrence multisets are
    # still printed on both sides for audit.
    missing = set(before) - set(after)
    added = set(after) - set(before)
    profile_cancellation = (
        allow_profile_cancellation
        and not added
        and missing <= BRIDGE_D_ORDER_PRESERVING_JETS
    )
    label = (
        "JET_CONSERVED"
        if (not missing and not added) or profile_cancellation
        else "JET_LOST"
    )
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
    left_pre_bridge_d: object | None,
    right_pre_bridge_d: object | None,
    state: RunState,
    *,
    divergence: bool,
) -> None:
    if _coverage(left) or _coverage(right):
        residual: object = C.TextAtom("UNDEFINED_UNJOINED")
        state.route_counts["COVERAGE"] += 1
        print(f"COVERAGE {case.family} {case.rendered_key}", flush=True)
        _emit_case_objects(case, left, right, residual)
        return

    if left is None or right is None:
        raise C.InputError("coverage routing did not retain its missing operand")

    if case.family in ENERGY_STRUCTURE_FAMILIES or case.family == WEAK_ORIGIN_FAMILY:
        residual = C.TextAtom("UNDEFINED_STRUCTURE_INCOMPLETE")
        state.route_counts["STRUCTURE_INCOMPLETE"] += 1
        print(f"STRUCTURE_INCOMPLETE {case.family} {case.rendered_key}", flush=True)
        if case.family == WEAK_ORIGIN_FAMILY:
            print(
                "weak_origin_sources = "
                + C.serialise(tuple(C.TextAtom(item) for item in WEAK_ORIGIN_SOURCES)),
                flush=True,
            )
        print(
            f"shape_diff {case.family} {case.rendered_key} = "
            f"{C.serialise(_shape_difference(left, right))}",
            flush=True,
        )
        _emit_case_objects(case, left, right, residual)
        return

    if _arithmetic_pair(left, right):
        residual = _arithmetic_residual(left, right)
        if _zero_object(residual):
            state.route_counts["MATCH"] += 1
            print(f"MATCH {case.family} {case.rendered_key}", flush=True)
            _emit_case_objects(case, left, right, residual)
            return

        protected = sorted(_free_symbol_names(residual) & PROTECTED_ATOM_NAMES)
        if protected:
            state.route_counts["PROTECTED_UNREDUCED"] += 1
            print(f"PROTECTED_UNREDUCED {case.family} {case.rendered_key}", flush=True)
            print(
                "protected_atoms = "
                + C.serialise(tuple(C.TextAtom(item) for item in protected)),
                flush=True,
            )
            _emit_case_objects(case, left, right, residual)
            return

        if _is_weak_scalar_density(case, residual):
            print(
                "weak_density_sources = "
                + C.serialise(tuple(C.TextAtom(item) for item in WEAK_DENSITY_SOURCES)),
                flush=True,
            )
            if not divergence:
                state.route_counts["FLAG"] += 1
                print(f"FLAG {case.family} {case.rendered_key}", flush=True)
                print("divergence_classifier = TextAtom('DROPPED')", flush=True)
                _emit_case_objects(case, left, right, residual)
                return
            if not (
                isinstance(left_pre_bridge_d, sp.Expr)
                and isinstance(right_pre_bridge_d, sp.Expr)
            ):
                raise C.InputError("weak scalar density lost its anchored scalar operands")
            anchored_residual = _arithmetic_residual(
                left_pre_bridge_d, right_pre_bridge_d
            )
            if not isinstance(anchored_residual, sp.Expr):
                raise C.InputError("weak scalar density has a non-scalar anchored residual")
            result = classify_total_divergence(
                anchored_residual,
                PRODUCTION_FIELD_REGISTRY,
                apply_bridge_d=True,
            )
            state.route_counts[result.route] += 1
            print(f"{result.route} {case.family} {case.rendered_key}", flush=True)
            _emit_case_objects(case, left, right, residual)
            print(
                f"euler_signature {case.family} {case.rendered_key} = "
                f"{C.serialise(result.euler_signature)}",
                flush=True,
            )
            if result.vector is not None:
                print(
                    f"V_ANCHORED {case.family} {case.rendered_key} = "
                    f"{C.serialise(result.vector)}",
                    flush=True,
                )
            if result.anchored_remainder is not None:
                print(
                    f"anchored_R_minus_formal_div_V {case.family} {case.rendered_key} = "
                    f"{C.serialise(result.anchored_remainder)}",
                    flush=True,
                )
            if result.bridge_d_certificate is not None:
                print(
                    f"BridgeD_anchored_R_minus_formal_div_V {case.family} "
                    f"{case.rendered_key} = {C.serialise(result.bridge_d_certificate)}",
                    flush=True,
                )
            return

        state.route_counts["FLAG"] += 1
        print(f"FLAG {case.family} {case.rendered_key}", flush=True)
        _emit_case_objects(case, left, right, residual)
        return

    pairs = _kinetic_pairs(case.family, case.key, left, right)
    if pairs is not None:
        residual = _kinetic_residual(pairs)
        label = "MATCH" if _zero_object(residual) else "FLAG"
        protected = sorted(_free_symbol_names(residual) & PROTECTED_ATOM_NAMES)
        if protected:
            label = "PROTECTED_UNREDUCED"
        state.route_counts[label] += 1
        print(f"{label} {case.family} {case.rendered_key}", flush=True)
        if protected:
            print(
                "protected_atoms = "
                + C.serialise(tuple(C.TextAtom(item) for item in protected)),
                flush=True,
            )
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
    # v2 owns divergence handling.  Disable the comparator's legacy quotient
    # prepass before parsing so no HeldDiv can be discarded into a quotient.
    case.reduce_divergence = False
    introduced_error = C.materialize(case)
    if introduced_error or case.error is not None:
        raise C.InputError(case.error or "case materialization error")
    if any(label == "DIVERGENCE_REDUCED" for label, _value in case.context):
        raise C.InputError("DIVERGENCE_REDUCED context reached raw adjudication")
    return case.value


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


def _configuration_objects(
    active_names: dict[str, str],
    bridge_a: bool,
    bridge_d: bool,
    divergence: bool,
    collapse: CollapseRule | None,
) -> None:
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
    bridge_d_object: object = (
        tuple(
            sp.Eq(source, image, evaluate=False)
            for source, image in PROFILE_GRADE_SUBS.items()
        )
        if bridge_d
        else C.TextAtom("DROPPED")
    )
    print(f"BRIDGE_D_SUBSTITUTIONS = {C.serialise(bridge_d_object)}", flush=True)
    print(
        "BRIDGE_D_SOURCES = "
        + C.serialise(tuple(C.TextAtom(item) for item in BRIDGE_D_SOURCES)),
        flush=True,
    )
    print("BRIDGE_D_IMPORTED_OBJECT_IDENTITY = true", flush=True)
    print(
        "BRIDGE_D_ORDER_PRESERVING_JETS = "
        + C.serialise(
            _jet_object(Counter({item: 1 for item in BRIDGE_D_ORDER_PRESERVING_JETS}))
        ),
        flush=True,
    )
    print("LEGACY_DIVERGENCE_REDUCTION = TextAtom('DISABLED')", flush=True)
    print(
        f"DIVERGENCE_CLASSIFIER = {C.serialise(C.TextAtom('ENABLED' if divergence else 'DROPPED'))}",
        flush=True,
    )
    print("JET_CONSERVATION_BASIS = TextAtom('ORDERED_JET_SUPPORT')", flush=True)
    print(
        "PROTECTED_ATOMS = "
        + C.serialise(tuple(C.TextAtom(item) for item in sorted(PROTECTED_ATOM_NAMES))),
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


def _fixture_registry(*, include_profile: bool = False) -> FieldRegistry:
    bases = {"a", "φ", "f", "g", "h"}
    variational = set(bases)
    maps = {direction: {} for direction in range(1, 4)}
    if include_profile:
        bases.update({"W_bg", "w1_profile"})
        variational.discard("W_bg")
        variational.discard("w1_profile")
        maps = {
            direction: dict(PRODUCTION_FIELD_REGISTRY.derivative_maps[direction])
            for direction in range(1, 4)
        }
    return FieldRegistry(frozenset(bases), frozenset(variational), maps)


def _run_divergence_fixtures() -> None:
    a, phi = sp.symbols("a φ")
    a_d1, phi_d1 = sp.symbols("a_d1 φ_d1")
    registry = _fixture_registry()

    bulk_expression = a * phi_d1
    bulk = classify_total_divergence(
        bulk_expression, registry, apply_bridge_d=False
    )
    if bulk.route != "RESIDUAL_BULK" or not bulk.euler_signature.entries:
        raise C.InputError("variable-coefficient fixture route is inconsistent")
    print(f"FIXTURE_VARIABLE_COEFFICIENT = {C.serialise(C.TextAtom(bulk.route))}", flush=True)
    print(f"fixture_residual = {C.serialise(bulk_expression)}", flush=True)
    print(f"fixture_euler_signature = {C.serialise(bulk.euler_signature)}", flush=True)

    divergence_expression = a_d1 * phi + a * phi_d1
    exact = classify_total_divergence(
        divergence_expression, registry, apply_bridge_d=False
    )
    exact_vector = (a * phi, sp.S.Zero, sp.S.Zero)
    if (
        exact.route != "REPRESENTATIONAL_DIVERGENCE"
        or exact.vector is None
        or any(
            _normalise_exact(item - exact_item) != 0
            for item, exact_item in zip(exact.vector, exact_vector)
        )
        or exact.bridge_d_certificate != 0
    ):
        raise C.InputError("generic divergence fixture certificate is inconsistent")
    print(f"FIXTURE_GENERIC_DIVERGENCE = {C.serialise(C.TextAtom(exact.route))}", flush=True)
    print(f"fixture_residual = {C.serialise(divergence_expression)}", flush=True)
    print(f"fixture_V = {C.serialise(exact.vector)}", flush=True)
    print(
        f"fixture_R_minus_formal_div_V = {C.serialise(exact.bridge_d_certificate)}",
        flush=True,
    )

    strong_case = JoinedCase(
        "SLAB_OPERATOR",
        (("OBJECT", "DIVERGENCE_ROUTE_FIXTURE"),),
        divergence_expression,
        sp.S.Zero,
    )
    if _is_weak_scalar_density(strong_case, divergence_expression):
        raise C.InputError("strong route fixture became divergence-eligible")
    fixture_state = RunState(
        Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), []
    )
    _classify_case(
        strong_case,
        divergence_expression,
        sp.S.Zero,
        divergence_expression,
        sp.S.Zero,
        fixture_state,
        divergence=True,
    )
    strong_route = next(iter(fixture_state.route_counts), "")
    if strong_route != "FLAG":
        raise C.InputError("strong route fixture did not remain exact")
    print(f"FIXTURE_STRONG_ROUTE = {C.serialise(C.TextAtom(strong_route))}", flush=True)

    profile_registry = _fixture_registry(include_profile=True)
    W_bg = sp.Symbol("W_bg")
    anchored_vector = (W_bg * phi, sp.S.Zero, sp.S.Zero)
    anchored_expression = formal_divergence(anchored_vector, profile_registry)
    profile = classify_total_divergence(
        anchored_expression, profile_registry, apply_bridge_d=True
    )
    if profile.route != "REPRESENTATIONAL_DIVERGENCE" or profile.bridge_d_certificate != 0:
        raise C.InputError("anchored-profile fixture certificate is inconsistent")
    noncommutator = _normalise_exact(
        _bridge_d(anchored_expression)
        - formal_divergence(tuple(_bridge_d(item) for item in anchored_vector), profile_registry)
    )
    if noncommutator == 0 or sp.Symbol("W_bg_d1") not in anchored_expression.free_symbols:
        raise C.InputError("anchored-profile fixture lost non-commutation")
    print(f"FIXTURE_ANCHORED_PROFILE = {C.serialise(C.TextAtom(profile.route))}", flush=True)
    print(f"fixture_anchored_residual = {C.serialise(anchored_expression)}", flush=True)
    print(f"fixture_anchored_V = {C.serialise(profile.vector)}", flush=True)
    print(
        "fixture_BridgeD_R_minus_formal_div_V = "
        + C.serialise(profile.bridge_d_certificate),
        flush=True,
    )
    print(f"fixture_BridgeD_D_noncommutator = {C.serialise(noncommutator)}", flush=True)


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
    hooks.add_argument("--drop-bridge-d", action="store_true")
    hooks.add_argument("--drop-divergence", action="store_true")
    hooks.add_argument("--drop-rename", metavar="WLname")
    arguments = parser.parse_args(argv)

    started = time.monotonic()
    try:
        if arguments.residual_leaf_budget <= 0:
            raise C.InputError("--residual-leaf-budget must be positive")
        if arguments.drop_rename is not None and arguments.drop_rename not in WL_TO_PY_RENAME:
            raise C.InputError(f"unknown spelling equivalence {arguments.drop_rename!r}")
        if PROFILE_GRADE_SUBS is not P.PROFILE_GRADE_SUBS:
            raise C.InputError("Bridge D does not reference the imported profile-grade map")
        if P.e_W_bg in PROFILE_GRADE_SUBS:
            raise C.InputError("Bridge D contains the excluded local thickness fraction")
        if any(
            _bridge_d(source) != image
            for source, image in PROFILE_GRADE_SUBS.items()
        ):
            raise C.InputError("Bridge D application differs from its imported entries")

        collapse = _collapse_rule(arguments.collapse_jet) if arguments.collapse_jet else None
        active_names = dict(WL_TO_PY_RENAME)
        if arguments.drop_rename is not None:
            active_names.pop(arguments.drop_rename)
            H._disable_imported_prepass_for_drop(arguments.drop_rename)

        py_tags = C.load_py(arguments.py)
        wl_tags = C.load_wl(arguments.wl)
        state = RunState(Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), [])
        bridge_a = not arguments.drop_bridge_a
        bridge_d = not arguments.drop_bridge_d
        divergence = not arguments.drop_divergence
        hook_active = (
            collapse is not None
            or arguments.drop_bridge_a
            or arguments.drop_bridge_d
            or arguments.drop_divergence
            or arguments.drop_rename is not None
        )

        _configuration_objects(active_names, bridge_a, bridge_d, divergence, collapse)
        print(
            f"RESIDUAL_LEAF_BUDGET_SECONDS = {C.serialise(sp.Float(arguments.residual_leaf_budget))}",
            flush=True,
        )
        if divergence:
            _run_divergence_fixtures()
        else:
            print("DIVERGENCE_FIXTURES = TextAtom('DROPPED')", flush=True)

        for family in H.CORE_FAMILIES:
            for case in _family_cases(family, py_tags, wl_tags, state):
                raw_values = (case.operand_a, case.operand_b)
                state.captured_names.update(_token_names(case.operand_a))
                state.captured_names.update(_token_names(case.operand_b))

                energy_gated = family in ENERGY_STRUCTURE_FAMILIES or family == ENERGY_EXACT_FAMILY
                if energy_gated:
                    pre_bridge_d_values = raw_values
                    transformed_values = raw_values
                    baseline_pre_bridge_d_values = raw_values
                    baseline_values = raw_values
                else:
                    pre_bridge_d_values = (
                        transform(
                            case.operand_a,
                            active_names,
                            bridge_a=bridge_a,
                            bridge_d=False,
                            collapse=collapse,
                        ),
                        transform(
                            case.operand_b,
                            active_names,
                            bridge_a=bridge_a,
                            bridge_d=False,
                            collapse=collapse,
                        ),
                    )
                    transformed_values = (
                        _bridge_d(pre_bridge_d_values[0])
                        if bridge_d and pre_bridge_d_values[0] is not None
                        else pre_bridge_d_values[0],
                        _bridge_d(pre_bridge_d_values[1])
                        if bridge_d and pre_bridge_d_values[1] is not None
                        else pre_bridge_d_values[1],
                    )
                    if _contains_held_divergence(transformed_values[0]) or _contains_held_divergence(
                        transformed_values[1]
                    ):
                        raise C.InputError(
                            f"unexpanded HeldDiv in {case.rendered_identity}"
                        )
                    if hook_active and not arguments.drop_divergence:
                        baseline_pre_bridge_d_values = (
                            transform(
                                case.operand_a,
                                dict(WL_TO_PY_RENAME),
                                bridge_a=True,
                                bridge_d=False,
                                collapse=None,
                            ),
                            transform(
                                case.operand_b,
                                dict(WL_TO_PY_RENAME),
                                bridge_a=True,
                                bridge_d=False,
                                collapse=None,
                            ),
                        )
                        baseline_values = (
                            _bridge_d(baseline_pre_bridge_d_values[0])
                            if baseline_pre_bridge_d_values[0] is not None
                            else None,
                            _bridge_d(baseline_pre_bridge_d_values[1])
                            if baseline_pre_bridge_d_values[1] is not None
                            else None,
                        )
                    else:
                        baseline_pre_bridge_d_values = pre_bridge_d_values
                        baseline_values = transformed_values
                if bridge_a and collapse is None and arguments.drop_rename is None:
                    for value in transformed_values:
                        if value is not None and any(BRHO in basic.atoms(sp.Symbol) for basic in _iter_basics(value)):
                            if _arithmetic_shape(value) is not None:
                                raise C.InputError(f"residual bRho in arithmetic operand {case.rendered_identity}")

                if hook_active and not arguments.drop_divergence:
                    if not (
                        _same_syntax(baseline_values[0], transformed_values[0])
                        and _same_syntax(baseline_values[1], transformed_values[1])
                    ):
                        state.touched_cases.append(case.rendered_identity)
                        _emit_ablation_case(case, baseline_values, transformed_values, active_names)

                if arguments.drop_divergence and not (
                    _coverage(transformed_values[0]) or _coverage(transformed_values[1])
                ):
                    if (
                        transformed_values[0] is not None
                        and transformed_values[1] is not None
                        and _arithmetic_pair(*transformed_values)
                    ):
                        candidate = _arithmetic_residual(*transformed_values)
                        if (
                            not _zero_object(candidate)
                            and not (_free_symbol_names(candidate) & PROTECTED_ATOM_NAMES)
                            and _is_weak_scalar_density(case, candidate)
                        ):
                            state.touched_cases.append(case.rendered_identity)
                            _emit_ablation_case(
                                case, transformed_values, transformed_values, active_names
                            )

                _classify_case(
                    case,
                    *transformed_values,
                    *pre_bridge_d_values,
                    state,
                    divergence=divergence,
                )
                state.classified_ids[case.identity] += 1
                _emit_jet_line(
                    case,
                    baseline_pre_bridge_d_values,
                    transformed_values,
                    active_names,
                    state,
                    allow_profile_cancellation=(
                        collapse is None and arguments.drop_rename is None
                    ),
                )

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
        emitted_total = sum(state.emitted_ids.values())
        if sum(state.comparator_counts.values()) != emitted_total:
            raise C.InputError("join/one-sided accounting differs from emitted cases")
        if sum(state.route_counts.values()) != emitted_total:
            raise C.InputError("route accounting differs from emitted cases")

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
                            "MATCH",
                            "FLAG",
                            "REPRESENTATIONAL_DIVERGENCE",
                            "RESIDUAL_BULK",
                            "DIVERGENCE_INCOMPLETE",
                            "PROTECTED_UNREDUCED",
                            "STRUCTURE_INCOMPLETE",
                            "COVERAGE",
                        )
                    )
                )
            ),
            flush=True,
        )
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
