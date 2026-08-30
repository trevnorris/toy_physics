#!/usr/bin/env python3
"""Emit requested-truncation residuals for complete S11c-b semantic rows.

All physical operands are loaded from the two committed transcripts through
``S11c_b_adjudicated_comparison``.  This instrument aligns those operands,
Taylor-projects operands to the requested eta/sigma_W window appropriate to
each row, and dispatches exact, weak, or componentwise comparison according
to row type.  It emits measurements and instrument-integrity accounting; it
states no cross-engine verdict.
"""

from __future__ import annotations

import hashlib
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import sympy as sp
from sympy.core.relational import Relational


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import S11c_b_adjudicated_comparison as A  # noqa: E402


# Use the exact module objects imported by the reviewed layer.
C = A.C
P = A.P
eta_bg = P.eta_bg
sigma_W = P.sigma_W

STRONG_EXACT = "STRONG_EXACT"
WEAK_MODULO_DIVERGENCE = "WEAK_MODULO_EXACT_TOTAL_DIVERGENCE"
COMPONENTWISE = "COMPONENTWISE_ORDERED_PAIRING"
FACE_ATTRIBUTION = "FACE_ATTRIBUTION_ONLY"

SLAB_OBJECTS = frozenset(
    {
        "U_MOMENTUM_ROWS",
        "MU_THETA",
        "THICKNESS_ROW",
        "MASS_EVOLUTION_ROW",
    }
)
SLAB_ORIGIN_OBJECTS = frozenset({"FACE", "FLUX", "KINETIC"})


@dataclass(frozen=True)
class AlignedCase:
    case: A.JoinedCase
    py: object | None
    wl: object | None
    py_pre_bridge_d: object | None
    wl_pre_bridge_d: object | None


@dataclass(frozen=True)
class AccountingEntry:
    token: str
    engine: str
    role: str
    category: str


@dataclass(frozen=True)
class RowMeasurement:
    case: A.JoinedCase
    path: str
    equivalence: str
    epsilon_order: int
    py_trunc: sp.Expr
    wl_complete: sp.Expr
    face_attributed: sp.Expr | None
    residual: sp.Expr
    certificate: object | None
    accounting: tuple[AccountingEntry, ...]


@dataclass(frozen=True)
class CouplingMeasurement:
    case: A.JoinedCase
    path: str
    epsilon_order: int
    py_trunc: sp.Expr
    wl_trunc: sp.Expr
    residual: sp.Expr
    full_prebridge_route: str
    euler_signature_presentation: object
    euler_signature_sha256: str
    in_scope_weak_remainder: sp.Expr | None
    no_clean_quotient_error: str | None
    accounting: tuple[AccountingEntry, ...]


class AssemblyLedger:
    """Track discovered extraction leaves independently of row consumption."""

    def __init__(self) -> None:
        self.expected: Counter[str] = Counter()
        self.accounted: Counter[str] = Counter()
        self.categories: dict[str, str] = {}

    def discover(self, token: str, category: str) -> None:
        self.expected[token] += 1
        self.categories[token] = category

    def consume(self, entries: Iterable[AccountingEntry]) -> None:
        for entry in entries:
            self.accounted[entry.token] += 1

    def counts(self, category: str | None = None) -> tuple[int, int]:
        if category is None:
            tokens = set(self.expected) | set(self.accounted)
        else:
            tokens = {
                token
                for token in set(self.expected) | set(self.accounted)
                if self.categories.get(token) == category
            }
        return (
            sum(self.expected[token] for token in tokens),
            sum(self.accounted[token] for token in tokens),
        )

    def difference_object(self) -> C.Association:
        missing = self.expected - self.accounted
        extra = self.accounted - self.expected
        return C.Association(
            (
                (
                    "MISSING",
                    tuple(C.TextAtom(token) for token in sorted(missing.elements())),
                ),
                (
                    "EXTRA",
                    tuple(C.TextAtom(token) for token in sorted(extra.elements())),
                ),
            )
        )


def _emit(name: str, payload: object) -> None:
    print(f"{name} = {C.serialise(payload)}", flush=True)


def _map_arithmetic(value: object, scalar_function) -> object:
    if isinstance(value, C.Association):
        return C.Association(
            tuple(
                (label, _map_arithmetic(item, scalar_function))
                for label, item in value.entries
            )
        )
    if isinstance(value, tuple):
        return tuple(_map_arithmetic(item, scalar_function) for item in value)
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(scalar_function)
    if isinstance(value, Relational):
        return type(value)(
            _map_arithmetic(value.lhs, scalar_function),
            _map_arithmetic(value.rhs, scalar_function),
            evaluate=False,
        )
    if isinstance(value, sp.Expr):
        return scalar_function(value)
    return value


def _requested_truncation_scalar(expression: sp.Expr) -> sp.Expr:
    """Bivariate first-order Taylor projection by origin derivatives."""

    origin = {eta_bg: sp.S.Zero, sigma_W: sp.S.Zero}
    c00 = expression.subs(origin, simultaneous=True)
    c10 = sp.diff(expression, eta_bg).subs(origin, simultaneous=True)
    c01 = sp.diff(expression, sigma_W).subs(origin, simultaneous=True)
    c11 = sp.diff(expression, eta_bg, sigma_W).subs(
        origin, simultaneous=True
    )
    return sp.expand(
        c00 + eta_bg * c10 + sigma_W * c01 + eta_bg * sigma_W * c11
    )


def requested_truncation(value: object) -> object:
    """Project every arithmetic leaf to eta<=1 and sigma_W<=1."""

    return _map_arithmetic(value, _requested_truncation_scalar)


def independent_series_projection(expression: sp.Expr) -> sp.Expr:
    """Independent fixture route: sequential direct SymPy series calls."""

    eta_projected = sp.series(expression, eta_bg, 0, 2).removeO()
    sigma_projected = sp.series(eta_projected, sigma_W, 0, 2).removeO()
    return sp.expand(sigma_projected)


def _iter_expressions(value: object) -> Iterator[sp.Expr]:
    if isinstance(value, C.Association):
        for _label, item in value.entries:
            yield from _iter_expressions(item)
    elif isinstance(value, tuple):
        for item in value:
            yield from _iter_expressions(item)
    elif isinstance(value, sp.MatrixBase):
        for item in value:
            yield from _iter_expressions(item)
    elif isinstance(value, Relational):
        yield from _iter_expressions(value.lhs)
        yield from _iter_expressions(value.rhs)
    elif isinstance(value, sp.Expr):
        yield value


def _coefficient_at_grade(expression: sp.Expr, a: int, b: int) -> sp.Expr:
    coefficient = sp.diff(expression, eta_bg, a, sigma_W, b).subs(
        {eta_bg: sp.S.Zero, sigma_W: sp.S.Zero}, simultaneous=True
    )
    return sp.expand(coefficient / (sp.factorial(a) * sp.factorial(b)))


def multigrade_support(value: object) -> tuple[sp.Tuple, ...]:
    support: set[tuple[int, int]] = set()
    for expression in _iter_expressions(value):
        for a in range(2):
            for b in range(2):
                if _coefficient_at_grade(expression, a, b) != 0:
                    support.add((a, b))
    return tuple(
        sp.Tuple(sp.Integer(a), sp.Integer(b)) for a, b in sorted(support)
    )


def _cas_record(value: object, epsilon_order: int) -> C.Association:
    return C.Association(
        (
            ("VALUE", value),
            ("EPSILON_ORDER_FROM_FAMILY_METADATA", sp.Integer(epsilon_order)),
            ("ETA_SIGMAW_MULTIGRADE_SUPPORT", multigrade_support(value)),
        )
    )


def _row_case_object(
    case: A.JoinedCase, path: str, equivalence: str, epsilon_order: int
) -> C.Association:
    return C.Association(
        (
            ("FAMILY", C.TextAtom(case.family)),
            (
                "KEY",
                C.Association(
                    tuple(
                        (axis, C.TextAtom(value)) for axis, value in case.key
                    )
                ),
            ),
            ("LEAF_PATH", C.TextAtom(path)),
            ("EQUIVALENCE", C.TextAtom(equivalence)),
            ("EPSILON_ORDER_FROM_FAMILY_METADATA", sp.Integer(epsilon_order)),
        )
    )


def _accounting_object(
    entries: Sequence[AccountingEntry], epsilon_order: int
) -> C.Association:
    by_engine: dict[str, list[object]] = {"PY": [], "WL": []}
    for entry in entries:
        by_engine[entry.engine].append(
            C.Association(
                (
                    ("LAYER_EXTRACTED_OPERAND", C.TextAtom(entry.token)),
                    ("ASSEMBLY_ROLE", C.TextAtom(entry.role)),
                )
            )
        )
    return _cas_record(
        C.Association(
            (
                ("PY", tuple(by_engine["PY"])),
                ("WL", tuple(by_engine["WL"])),
            )
        ),
        epsilon_order,
    )


def _token(case: A.JoinedCase, engine: str, path: str) -> str:
    return (
        f"{engine}:{case.family}:{C.serialise_key(case.key)}:LEAF={path}"
    )


def _entry(
    case: A.JoinedCase,
    engine: str,
    path: str,
    role: str,
    category: str,
) -> AccountingEntry:
    return AccountingEntry(_token(case, engine, path), engine, role, category)


def _align_value(
    value: object | None,
) -> tuple[object | None, object | None]:
    if value is None:
        return None, None
    pre_bridge_d = A.transform(
        value,
        dict(A.WL_TO_PY_RENAME),
        bridge_a=True,
        bridge_d=False,
        collapse=None,
    )
    if pre_bridge_d is None:
        return None, None
    aligned = A._bridge_d(pre_bridge_d)
    if A._contains_held_divergence(aligned):
        raise C.InputError("aligned operand retains a held divergence")
    return pre_bridge_d, aligned


def _align_case(case: A.JoinedCase) -> AlignedCase:
    py_pre_bridge_d, py = _align_value(case.operand_a)
    wl_pre_bridge_d, wl = _align_value(case.operand_b)
    return AlignedCase(
        case,
        py,
        wl,
        py_pre_bridge_d,
        wl_pre_bridge_d,
    )


def _tuple_leaves(value: object, root: str = "ROOT") -> Mapping[str, sp.Expr]:
    output: dict[str, sp.Expr] = {}

    def visit(item: object, path: str) -> None:
        if isinstance(item, tuple):
            for index, child in enumerate(item):
                visit(child, f"{path}[{index}]")
            return
        if isinstance(item, sp.MatrixBase):
            for index, child in enumerate(item):
                visit(child, f"{path}[{index}]")
            return
        if not isinstance(item, sp.Expr):
            raise C.InputError(
                f"arithmetic leaf {path} has type {type(item).__name__}"
            )
        output[path] = item

    visit(value, root)
    return output


def _association_item(value: object, label: str) -> object:
    if not isinstance(value, C.Association):
        raise C.InputError(f"expected Association while selecting {label}")
    items = value.as_dict()
    if label not in items:
        raise C.InputError(f"Association lacks required label {label}")
    return items[label]


def _exact_scalar_residual(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    residual = A._arithmetic_residual(left, right)
    if not isinstance(residual, sp.Expr):
        raise C.InputError("scalar comparison produced a non-scalar residual")
    return residual


def _deterministic_witness_presentation(value: object) -> object:
    """Copy a witness into a hash-order-independent display form only."""

    if isinstance(value, C.Association):
        return C.Association(
            tuple(
                (
                    label,
                    _deterministic_witness_presentation(item),
                )
                for label, item in sorted(
                    value.entries, key=lambda entry: entry[0]
                )
            )
        )
    if isinstance(value, tuple):
        return tuple(
            _deterministic_witness_presentation(item) for item in value
        )
    if isinstance(value, sp.MatrixBase):
        return tuple(
            _deterministic_witness_presentation(item) for item in value
        )
    if isinstance(value, sp.Expr):
        canonical_copy = sp.signsimp(value)
        return C.RawCAS(
            "INSTRUMENT_WITNESS_PRESENTATION",
            sp.sstr(canonical_copy, order="lex"),
        )
    return value


def _witness_sha256(value: object) -> str:
    payload = C.serialise(value).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def equivalence_dispatch(
    equivalence: str,
    case: A.JoinedCase,
    py_trunc: sp.Expr,
    wl_complete: sp.Expr,
    face_attributed: sp.Expr | None,
) -> tuple[sp.Expr, object | None]:
    """Apply the directive's exact/weak/componentwise equivalence split."""

    if equivalence == STRONG_EXACT:
        wl_residual_operand = (
            _exact_scalar_residual(wl_complete, face_attributed)
            if face_attributed is not None
            else wl_complete
        )
        return _exact_scalar_residual(py_trunc, wl_residual_operand), None
    if equivalence == COMPONENTWISE:
        return _exact_scalar_residual(py_trunc, wl_complete), None
    raise C.InputError(f"unknown equivalence dispatch {equivalence}")


def _emit_row(measurement: RowMeasurement, ledger: AssemblyLedger) -> None:
    _emit(
        "CASE",
        _row_case_object(
            measurement.case,
            measurement.path,
            measurement.equivalence,
            measurement.epsilon_order,
        ),
    )
    _emit(
        "ROW_OPERAND_PY_TRUNC",
        _cas_record(measurement.py_trunc, measurement.epsilon_order),
    )
    _emit(
        "ROW_OPERAND_WL",
        _cas_record(measurement.wl_complete, measurement.epsilon_order),
    )
    if measurement.face_attributed is not None:
        _emit(
            "ROW_FACE_ATTRIBUTED",
            _cas_record(
                measurement.face_attributed, measurement.epsilon_order
            ),
        )
    _emit(
        "ROW_RESIDUAL",
        _cas_record(measurement.residual, measurement.epsilon_order),
    )
    if measurement.certificate is not None:
        _emit(
            "ROW_DIVERGENCE_CERTIFICATE",
            _cas_record(measurement.certificate, measurement.epsilon_order),
        )
    _emit(
        "ROW_ASSEMBLY_ACCOUNTING",
        _accounting_object(measurement.accounting, measurement.epsilon_order),
    )
    ledger.consume(measurement.accounting)


def _emit_coupling_row(
    measurement: CouplingMeasurement, ledger: AssemblyLedger
) -> None:
    _emit(
        "CASE",
        _row_case_object(
            measurement.case,
            measurement.path,
            WEAK_MODULO_DIVERGENCE,
            measurement.epsilon_order,
        ),
    )
    _emit(
        "FULL_PREBRIDGE_ROUTE",
        C.TextAtom(measurement.full_prebridge_route),
    )
    _emit(
        "EULER_SIGNATURE",
        measurement.euler_signature_presentation,
    )
    _emit(
        "EULER_SIGNATURE_SHA256",
        C.TextAtom(measurement.euler_signature_sha256),
    )
    _emit(
        "ROW_OPERAND_PY_TRUNC",
        _cas_record(measurement.py_trunc, measurement.epsilon_order),
    )
    _emit(
        "ROW_OPERAND_WL",
        _cas_record(measurement.wl_trunc, measurement.epsilon_order),
    )
    _emit(
        "ROW_RESIDUAL",
        _cas_record(measurement.residual, measurement.epsilon_order),
    )
    no_clean_quotient = measurement.in_scope_weak_remainder is None
    _emit("NO_CLEAN_QUOTIENT", sp.true if no_clean_quotient else sp.false)
    if measurement.in_scope_weak_remainder is not None:
        _emit(
            "IN_SCOPE_WEAK_REMAINDER",
            _cas_record(
                measurement.in_scope_weak_remainder,
                measurement.epsilon_order,
            ),
        )
    if measurement.no_clean_quotient_error is not None:
        _emit(
            "NO_CLEAN_QUOTIENT_ERROR",
            C.TextAtom(measurement.no_clean_quotient_error),
        )
    _emit(
        "ROW_ASSEMBLY_ACCOUNTING",
        _accounting_object(measurement.accounting, measurement.epsilon_order),
    )
    ledger.consume(measurement.accounting)


def _emit_face_only(
    case: A.JoinedCase,
    path: str,
    value: sp.Expr,
    accounting: tuple[AccountingEntry, ...],
    ledger: AssemblyLedger,
) -> None:
    _emit("CASE", _row_case_object(case, path, FACE_ATTRIBUTION, 1))
    _emit("ROW_FACE_ATTRIBUTED", _cas_record(value, 1))
    _emit("ROW_ASSEMBLY_ACCOUNTING", _accounting_object(accounting, 1))
    ledger.consume(accounting)


def _load_aligned_cases(
    py_tags: dict[str, str],
    wl_tags: dict[str, str],
    state: A.RunState,
) -> tuple[
    dict[tuple[str, str, str], AlignedCase],
    dict[tuple[str, str, str], AlignedCase],
    list[AlignedCase],
    list[AlignedCase],
]:
    slab: dict[tuple[str, str, str], AlignedCase] = {}
    origins: dict[tuple[str, str, str], AlignedCase] = {}
    coupling: list[AlignedCase] = []
    admissibility: list[AlignedCase] = []

    for case in A._family_cases("SLAB_OPERATOR", py_tags, wl_tags, state):
        axes = C.key_dict(case.key)
        object_name = axes.get("OBJECT")
        if object_name not in SLAB_OBJECTS:
            continue
        key = (axes["BRANCH"], axes["DENSITY"], object_name)
        slab[key] = _align_case(case)

    for case in A._family_cases(
        "SLAB_OPERATOR_TERM_ORIGINS", py_tags, wl_tags, state
    ):
        axes = C.key_dict(case.key)
        object_name = axes.get("OBJECT")
        if object_name not in SLAB_ORIGIN_OBJECTS:
            continue
        key = (axes["BRANCH"], axes["DENSITY"], object_name)
        origins[key] = _align_case(case)

    for case in A._family_cases("COUPLING_KERNEL", py_tags, wl_tags, state):
        axes = C.key_dict(case.key)
        if axes.get("OBJECT") == "ADJOINTNESS_RELATION":
            continue
        aligned = _align_case(case)
        if aligned.py is None or aligned.wl is None:
            raise C.InputError(
                f"coupling row is not aligned: {case.rendered_identity}"
            )
        coupling.append(aligned)

    for case in A._family_cases(
        "ADMISSIBILITY_OPERATOR_OPERAND", py_tags, wl_tags, state
    ):
        aligned = _align_case(case)
        if aligned.py is None or aligned.wl is None:
            raise C.InputError(
                f"admissibility component is not aligned: {case.rendered_identity}"
            )
        admissibility.append(aligned)

    return slab, origins, coupling, admissibility


def _discover_slab(
    slab: Mapping[tuple[str, str, str], AlignedCase],
    origins: Mapping[tuple[str, str, str], AlignedCase],
    ledger: AssemblyLedger,
) -> dict[tuple[str, str], tuple[tuple[str, sp.Expr, sp.Expr], ...]]:
    kinetic_pairs: dict[
        tuple[str, str], tuple[tuple[str, sp.Expr, sp.Expr], ...]
    ] = {}
    for (branch, density, object_name), aligned in slab.items():
        if aligned.py is None or aligned.wl is None:
            raise C.InputError(
                f"strong row is not aligned: {aligned.case.rendered_identity}"
            )
        py_leaves = _tuple_leaves(aligned.py)
        wl_leaves = _tuple_leaves(aligned.wl)
        if set(py_leaves) != set(wl_leaves):
            raise C.InputError(
                f"strong row leaf mismatch: {aligned.case.rendered_identity}"
            )
        category = (
            "MASS_ROW"
            if object_name == "MASS_EVOLUTION_ROW"
            else "THICKNESS_KINETIC_ROW"
            if object_name == "THICKNESS_ROW"
            else "SLAB_OTHER"
        )
        for path in py_leaves:
            ledger.discover(_token(aligned.case, "PY", path), category)
            ledger.discover(_token(aligned.case, "WL", path), category)

    for (branch, density, object_name), aligned in origins.items():
        if object_name == "KINETIC":
            if aligned.py is None or aligned.wl is None:
                raise C.InputError("kinetic origin lacks one engine operand")
            pairs = A._kinetic_pairs(
                aligned.case.family,
                aligned.case.key,
                aligned.py,
                aligned.wl,
            )
            if pairs is None:
                raise C.InputError(
                    f"kinetic adapter rejected {aligned.case.rendered_identity}"
                )
            kinetic_pairs[(branch, density)] = pairs
            for label, _py, _wl in pairs:
                category = (
                    "THICKNESS_KINETIC_ROW"
                    if label == "THICKNESS_ROW"
                    else "SLAB_OTHER"
                )
                ledger.discover(_token(aligned.case, "PY", label), category)
                ledger.discover(_token(aligned.case, "WL", label), category)
        elif object_name == "FACE":
            if aligned.wl is None:
                raise C.InputError("WL FACE origin is missing")
            for label in (
                "U_MOMENTUM_ROWS",
                "THICKNESS_ROW",
                "CENTER_FACE_GENERALIZED_ROW",
            ):
                leaves = _tuple_leaves(_association_item(aligned.wl, label), label)
                category = (
                    "THICKNESS_KINETIC_ROW"
                    if label == "THICKNESS_ROW"
                    else "SLAB_OTHER"
                )
                for path in leaves:
                    ledger.discover(_token(aligned.case, "WL", path), category)
        elif object_name == "FLUX":
            if aligned.wl is None:
                raise C.InputError("WL FLUX origin is missing")
            for path in _tuple_leaves(aligned.wl, "MASS_EVOLUTION_ROW"):
                ledger.discover(_token(aligned.case, "WL", path), "MASS_ROW")
    return kinetic_pairs


def _discover_scalar_aligned(
    aligned_cases: Iterable[AlignedCase],
    ledger: AssemblyLedger,
    category: str,
) -> None:
    for aligned in aligned_cases:
        if aligned.py is None or aligned.wl is None:
            raise C.InputError(
                f"discovery received an unaligned case: {aligned.case.rendered_identity}"
            )
        py_leaves = _tuple_leaves(aligned.py)
        wl_leaves = _tuple_leaves(aligned.wl)
        if set(py_leaves) != set(wl_leaves):
            raise C.InputError(
                f"component leaf mismatch: {aligned.case.rendered_identity}"
            )
        for path in py_leaves:
            ledger.discover(_token(aligned.case, "PY", path), category)
            ledger.discover(_token(aligned.case, "WL", path), category)


def _strong_rows(
    slab: Mapping[tuple[str, str, str], AlignedCase],
    origins: Mapping[tuple[str, str, str], AlignedCase],
    kinetic_pairs: Mapping[
        tuple[str, str], tuple[tuple[str, sp.Expr, sp.Expr], ...]
    ],
) -> Iterator[RowMeasurement]:
    for key in sorted(slab):
        branch, density, object_name = key
        aligned = slab[key]
        if aligned.py is None or aligned.wl is None:
            raise C.InputError("strong row lost an aligned operand")
        py_leaves = _tuple_leaves(aligned.py)
        wl_leaves = _tuple_leaves(aligned.wl)

        face_case: AlignedCase | None = None
        face_leaves: Mapping[str, sp.Expr] = {}
        face_label: str | None = None
        if object_name in {"U_MOMENTUM_ROWS", "THICKNESS_ROW"}:
            face_case = origins[(branch, density, "FACE")]
            face_label = object_name
            face_value = _association_item(face_case.wl, face_label)
            face_leaves = _tuple_leaves(face_value, face_label)
        elif object_name == "MASS_EVOLUTION_ROW":
            face_case = origins[(branch, density, "FLUX")]
            face_label = "MASS_EVOLUTION_ROW"
            face_leaves = _tuple_leaves(face_case.wl, face_label)

        kinetic_by_label = {
            label: (py, wl)
            for label, py, wl in kinetic_pairs[(branch, density)]
        }
        for path in sorted(py_leaves):
            py_trunc = requested_truncation(py_leaves[path])
            if not isinstance(py_trunc, sp.Expr):
                raise C.InputError("requested truncation lost a scalar leaf")
            face_path = (
                path.replace("ROOT", face_label, 1)
                if face_label is not None
                else None
            )
            face_attributed = (
                face_leaves[face_path]
                if face_path is not None and face_path in face_leaves
                else None
            )
            residual, certificate = equivalence_dispatch(
                STRONG_EXACT,
                aligned.case,
                py_trunc,
                wl_leaves[path],
                face_attributed,
            )
            category = (
                "MASS_ROW"
                if object_name == "MASS_EVOLUTION_ROW"
                else "THICKNESS_KINETIC_ROW"
                if object_name == "THICKNESS_ROW"
                else "SLAB_OTHER"
            )
            accounting = [
                _entry(
                    aligned.case,
                    "PY",
                    path,
                    "COMPLETE_ROW_INPUT",
                    category,
                ),
                _entry(
                    aligned.case,
                    "WL",
                    path,
                    "COMPLETE_ROW_INPUT",
                    category,
                ),
            ]
            if face_case is not None and face_path is not None:
                accounting.append(
                    _entry(
                        face_case.case,
                        "WL",
                        face_path,
                        "FACE_ATTRIBUTED_AND_EXCLUDED_FROM_RESIDUAL",
                        category,
                    )
                )
            kinetic_label = (
                "THICKNESS_ROW"
                if object_name == "THICKNESS_ROW"
                else path.replace("ROOT", "U_MOMENTUM_ROWS", 1)
                if object_name == "U_MOMENTUM_ROWS"
                else None
            )
            if kinetic_label is not None:
                if kinetic_label not in kinetic_by_label:
                    raise C.InputError(
                        f"missing kinetic adapter leaf {kinetic_label}"
                    )
                kinetic_case = origins[(branch, density, "KINETIC")].case
                accounting.extend(
                    (
                        _entry(
                            kinetic_case,
                            "PY",
                            kinetic_label,
                            "PROVENANCE_ONLY_ALREADY_IN_COMPLETE_ROW",
                            category,
                        ),
                        _entry(
                            kinetic_case,
                            "WL",
                            kinetic_label,
                            "PROVENANCE_ONLY_ALREADY_IN_COMPLETE_ROW",
                            category,
                        ),
                    )
                )
            yield RowMeasurement(
                aligned.case,
                path,
                STRONG_EXACT,
                1,
                py_trunc,
                wl_leaves[path],
                face_attributed,
                residual,
                certificate,
                tuple(accounting),
            )


def _coupling_rows(
    coupling: Iterable[AlignedCase],
) -> Iterator[CouplingMeasurement]:
    for aligned in coupling:
        if not isinstance(aligned.py_pre_bridge_d, sp.Expr) or not isinstance(
            aligned.wl_pre_bridge_d, sp.Expr
        ):
            raise C.InputError(
                "weak row lacks scalar pre-bridge operands: "
                f"{aligned.case.rendered_identity}"
            )

        residual_pre_bridge_d = A._arithmetic_residual(
            aligned.py_pre_bridge_d,
            aligned.wl_pre_bridge_d,
        )
        if not isinstance(residual_pre_bridge_d, sp.Expr):
            raise C.InputError("weak pre-bridge residual is not scalar")

        full = A.classify_total_divergence(
            residual_pre_bridge_d,
            A.PRODUCTION_FIELD_REGISTRY,
            apply_bridge_d=True,
        )

        py_trunc = requested_truncation(
            A._bridge_d(aligned.py_pre_bridge_d)
        )
        wl_trunc = requested_truncation(
            A._bridge_d(aligned.wl_pre_bridge_d)
        )
        residual = requested_truncation(
            A._bridge_d(residual_pre_bridge_d)
        )
        if not all(
            isinstance(value, sp.Expr)
            for value in (py_trunc, wl_trunc, residual)
        ):
            raise C.InputError("weak bridged truncation lost a scalar operand")

        in_scope_weak_remainder: sp.Expr | None
        no_clean_quotient_error: str | None = None
        try:
            vector = A._homotopy_vector(
                residual_pre_bridge_d,
                A.PRODUCTION_FIELD_REGISTRY,
            )
            remainder_pre_bridge_d = A._normalise_exact(
                residual_pre_bridge_d
                - A.formal_divergence(
                    vector,
                    A.PRODUCTION_FIELD_REGISTRY,
                )
            )
            in_scope_weak_remainder = requested_truncation(
                A._bridge_d(remainder_pre_bridge_d)
            )
            if not isinstance(in_scope_weak_remainder, sp.Expr):
                raise C.InputError(
                    "weak remainder truncation did not return a scalar"
                )
        except Exception as error:
            in_scope_weak_remainder = None
            no_clean_quotient_error = type(error).__name__

        euler_signature_presentation = (
            _deterministic_witness_presentation(full.euler_signature)
        )
        accounting = (
            _entry(
                aligned.case,
                "PY",
                "ROOT",
                "WEAK_DENSITY_INPUT",
                "COUPLING",
            ),
            _entry(
                aligned.case,
                "WL",
                "ROOT",
                "WEAK_DENSITY_INPUT",
                "COUPLING",
            ),
        )
        yield CouplingMeasurement(
            aligned.case,
            "ROOT",
            1,
            py_trunc,
            wl_trunc,
            residual,
            full.route,
            euler_signature_presentation,
            _witness_sha256(euler_signature_presentation),
            in_scope_weak_remainder,
            no_clean_quotient_error,
            accounting,
        )


def _admissibility_rows(
    admissibility: Iterable[AlignedCase],
) -> Iterator[RowMeasurement]:
    for aligned in admissibility:
        if aligned.py is None or aligned.wl is None:
            raise C.InputError("admissibility row lost an aligned operand")
        py_leaves = _tuple_leaves(aligned.py)
        wl_leaves = _tuple_leaves(aligned.wl)
        for path in sorted(py_leaves):
            py_trunc = requested_truncation(py_leaves[path])
            if not isinstance(py_trunc, sp.Expr):
                raise C.InputError("component truncation lost a scalar operand")
            residual, certificate = equivalence_dispatch(
                COMPONENTWISE,
                aligned.case,
                py_trunc,
                wl_leaves[path],
                None,
            )
            accounting = (
                _entry(
                    aligned.case,
                    "PY",
                    path,
                    "ORDERED_COMPONENT_INPUT",
                    "ADMISSIBILITY",
                ),
                _entry(
                    aligned.case,
                    "WL",
                    path,
                    "ORDERED_COMPONENT_INPUT",
                    "ADMISSIBILITY",
                ),
            )
            yield RowMeasurement(
                aligned.case,
                path,
                COMPONENTWISE,
                0,
                py_trunc,
                wl_leaves[path],
                None,
                residual,
                certificate,
                accounting,
            )


def _emit_center_face_rows(
    origins: Mapping[tuple[str, str, str], AlignedCase],
    ledger: AssemblyLedger,
) -> None:
    for branch, density, object_name in sorted(origins):
        if object_name != "FACE":
            continue
        aligned = origins[(branch, density, object_name)]
        center = _association_item(aligned.wl, "CENTER_FACE_GENERALIZED_ROW")
        leaves = _tuple_leaves(center, "CENTER_FACE_GENERALIZED_ROW")
        for path, value in sorted(leaves.items()):
            accounting = (
                _entry(
                    aligned.case,
                    "WL",
                    path,
                    "FACE_ATTRIBUTION_ONLY",
                    "SLAB_OTHER",
                ),
            )
            _emit_face_only(
                aligned.case, path, value, accounting, ledger
            )


def _fixture_objects() -> tuple[sp.Expr, sp.Expr, sp.Expr, bool]:
    f0, f1, f2, f3, f4 = sp.symbols(
        "fixture_f0 fixture_f1 fixture_f2 fixture_f3 fixture_f4"
    )
    fixture = (
        (sp.S.One + eta_bg * f1) ** 2
        / (sp.S.One + eta_bg * f2)
        + sigma_W * f3
        + eta_bg * sigma_W * f4
        + eta_bg**2 * f0
        + sigma_W**2 * f1
    )
    primary = requested_truncation(fixture)
    independent = independent_series_projection(fixture)
    if not isinstance(primary, sp.Expr):
        raise C.InputError("fixture truncation did not return an expression")
    agrees = sp.expand(primary - independent) == 0
    return fixture, primary, independent, agrees


def _corruption_objects() -> tuple[sp.Expr, sp.Expr]:
    kept, eta_dropped, sigma_dropped = sp.symbols(
        "fixture_kept fixture_eta_dropped fixture_sigma_dropped"
    )
    before = (
        eta_bg * sigma_W * kept
        + eta_bg**2 * eta_dropped
        + sigma_W**2 * sigma_dropped
    )
    after = requested_truncation(before)
    if not isinstance(after, sp.Expr):
        raise C.InputError("corruption truncation did not return an expression")
    return before, after


def _guard_count_object(
    ledger: AssemblyLedger, category: str
) -> C.Association:
    expected, accounted = ledger.counts(category)
    return C.Association(
        (
            ("EXPECTED_LAYER_OPERANDS", sp.Integer(expected)),
            ("ACCOUNTED_LAYER_OPERANDS", sp.Integer(accounted)),
            ("COUNT_DIFFERENCE", sp.Integer(expected - accounted)),
        )
    )


def _coupling_accounting_object(
    coupling: Iterable[AlignedCase],
) -> C.Association:
    counts: Counter[tuple[str, str]] = Counter()
    for aligned in coupling:
        axes = C.key_dict(aligned.case.key)
        sector = axes.get("SECTOR", "NO_SECTOR")
        object_name = axes.get("OBJECT", "BLOCK")
        counts[(object_name, sector)] += 1
    return C.Association(
        tuple(
            (
                f"{object_name}:{sector}",
                sp.Integer(counts[(object_name, sector)]),
            )
            for object_name, sector in sorted(counts)
        )
    )


def run() -> int:
    try:
        fixture, primary, independent, truncation_agrees = _fixture_objects()
        corruption_before, corruption_after = _corruption_objects()

        _emit("A_MODULE_PATH", C.TextAtom(str(Path(A.__file__).resolve())))
        _emit("C_MODULE_PATH", C.TextAtom(str(Path(C.__file__).resolve())))
        _emit("PY_TRANSCRIPT_PATH", C.TextAtom(str(C.DEFAULT_PY.resolve())))
        _emit("WL_TRANSCRIPT_PATH", C.TextAtom(str(C.DEFAULT_WL.resolve())))
        _emit("REQUESTED_TRUNCATION_FIXTURE", fixture)
        _emit("REQUESTED_TRUNCATION_PRIMARY", primary)
        _emit("REQUESTED_TRUNCATION_INDEPENDENT_SERIES", independent)
        _emit("REQUESTED_TRUNCATION_CROSSCHECK", sp.true if truncation_agrees else sp.false)
        _emit("ONE_SIDED_CORRUPTION_BEFORE", corruption_before)
        _emit("ONE_SIDED_CORRUPTION_AFTER", corruption_after)

        py_tags = C.load_py(C.DEFAULT_PY)
        wl_tags = C.load_wl(C.DEFAULT_WL)
        state = A.RunState(
            Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), []
        )
        slab, origins, coupling, admissibility = _load_aligned_cases(
            py_tags, wl_tags, state
        )

        ledger = AssemblyLedger()
        kinetic_pairs = _discover_slab(slab, origins, ledger)
        _discover_scalar_aligned(coupling, ledger, "COUPLING")
        _discover_scalar_aligned(admissibility, ledger, "ADMISSIBILITY")

        emitted: Counter[str] = Counter()
        for measurement in _strong_rows(slab, origins, kinetic_pairs):
            _emit_row(measurement, ledger)
            emitted[STRONG_EXACT] += 1
        _emit_center_face_rows(origins, ledger)
        for measurement in _coupling_rows(coupling):
            _emit_coupling_row(measurement, ledger)
            emitted[WEAK_MODULO_DIVERGENCE] += 1
        for measurement in _admissibility_rows(admissibility):
            _emit_row(measurement, ledger)
            emitted[COMPONENTWISE] += 1

        _emit(
            "COUPLING_BLOCK_ACCOUNTING",
            _coupling_accounting_object(coupling),
        )
        _emit(
            "ASSEMBLY_ACCOUNTING_GUARD_MASS_ROW",
            _guard_count_object(ledger, "MASS_ROW"),
        )
        _emit(
            "ASSEMBLY_ACCOUNTING_GUARD_THICKNESS_KINETIC_ROW",
            _guard_count_object(ledger, "THICKNESS_KINETIC_ROW"),
        )
        expected_total, accounted_total = ledger.counts()
        _emit(
            "ASSEMBLY_ACCOUNTING_GUARD",
            C.Association(
                (
                    ("EXPECTED_LAYER_OPERANDS", sp.Integer(expected_total)),
                    ("ACCOUNTED_LAYER_OPERANDS", sp.Integer(accounted_total)),
                    (
                        "COUNT_DIFFERENCE",
                        sp.Integer(expected_total - accounted_total),
                    ),
                    ("MULTISET_DIFFERENCE", ledger.difference_object()),
                )
            ),
        )
        _emit(
            "EMITTED_ROW_ACCOUNTING",
            C.Association(
                tuple(
                    (name, sp.Integer(emitted[name]))
                    for name in (
                        STRONG_EXACT,
                        WEAK_MODULO_DIVERGENCE,
                        COMPONENTWISE,
                    )
                )
            ),
        )

        # Arithmetic and assembly guards intentionally follow every measured
        # residual emission.  No residual value is tested here.
        if not truncation_agrees:
            raise C.InputError(
                "requested_truncation differs from independent series projection"
            )
        if ledger.expected != ledger.accounted:
            raise C.InputError("assembly-accounting multiset mismatch")

        coupling_counts = _coupling_accounting_object(coupling).as_dict()
        required_coupling_counts = {
            "BLOCK:TRANSVERSE_TO_THICKNESS": 4,
            "BLOCK:THICKNESS_TO_TRANSVERSE": 4,
            "ADJOINTNESS_OPERAND_FORWARD:TRANSVERSE_TO_THICKNESS": 4,
            "ADJOINTNESS_OPERAND_REVERSE:THICKNESS_TO_TRANSVERSE": 4,
        }
        if any(
            coupling_counts.get(label) != sp.Integer(count)
            for label, count in required_coupling_counts.items()
        ):
            raise C.InputError("coupling block/adjoint assembly is incomplete")
        return 0
    except Exception as error:
        print(
            f"OPERATIONAL_ERROR {type(error).__name__}: {error}",
            file=sys.stderr,
            flush=True,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(run())
