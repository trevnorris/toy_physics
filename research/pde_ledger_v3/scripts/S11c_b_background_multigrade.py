#!/usr/bin/env python3
"""Measure background-bookkeeper multigrades for the twenty S11c-b cases.

The aligned operands and their residuals come exclusively from the committed
adjudication layer.  This instrument extracts exact Taylor coefficients in
``eta_bg`` and ``sigma_W`` on a fixed rectangular window and retains the exact
uncaptured remainder.  It emits decomposition objects and arithmetic guards;
it performs no physical classification.
"""

from __future__ import annotations

import sys
from collections import Counter
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Iterable, Mapping

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import S11c_b_adjudicated_comparison as A  # noqa: E402


# Use the comparator and engine module objects imported by the reviewed layer.
C = A.C
P = A.P

N = 4
eta_bg = P.eta_bg
sigma_W = P.sigma_W

# Symbol assumptions are part of SymPy symbol identity.  These assertions are
# intentionally object-identity checks rather than same-spelling checks.
assert eta_bg is P.eta_bg
assert sigma_W is P.sigma_W


def _key(**axes: str) -> C.Key:
    unknown = set(axes) - set(C.AXIS_ORDER)
    if unknown:
        raise C.InputError(f"unknown target axes {sorted(unknown)}")
    return tuple((axis, axes[axis]) for axis in C.AXIS_ORDER if axis in axes)


_BRANCHES = ("LAB_HELD", "MATERIAL_ADVECTED")
_DENSITIES = ("RHO4_CONSTANT", "RHOBR_CONSTANT")

# This is the closed case registry supplied by the build directive.  It is not
# inferred from routes or residual values.
TARGET_IDENTITIES: frozenset[tuple[str, C.Key]] = frozenset(
    [
        (
            "SLAB_OPERATOR_TERM_ORIGINS",
            _key(OBJECT=object_name, BRANCH=branch, DENSITY=density),
        )
        for object_name, branch, density in product(
            ("ADVECTIVE", "KINETIC"), _BRANCHES, _DENSITIES
        )
    ]
    + [
        (
            "ADMISSIBILITY_OPERATOR_OPERAND",
            _key(
                OBJECT="BODY_FORCE",
                BRANCH=branch,
                DENSITY=density,
                DOF="THETA",
            ),
        )
        for branch, density in product(_BRANCHES, _DENSITIES)
    ]
    + [
        (
            "COUPLING_KERNEL",
            _key(
                BRANCH=branch,
                DENSITY=density,
                SECTOR="TRANSVERSE_TO_THICKNESS",
            ),
        )
        for branch, density in product(_BRANCHES, _DENSITIES)
    ]
    + [
        (
            "COUPLING_KERNEL",
            _key(
                OBJECT="ADJOINTNESS_OPERAND_FORWARD",
                BRANCH=branch,
                DENSITY=density,
                SECTOR="TRANSVERSE_TO_THICKNESS",
            ),
        )
        for branch, density in product(_BRANCHES, _DENSITIES)
    ]
)

TARGET_FAMILIES = (
    "SLAB_OPERATOR_TERM_ORIGINS",
    "ADMISSIBILITY_OPERATOR_OPERAND",
    "COUPLING_KERNEL",
)

if len(TARGET_IDENTITIES) != 20:
    raise C.InputError("the closed target registry does not contain twenty identities")


@dataclass(frozen=True)
class LeafMultigrade:
    coefficients: Mapping[tuple[int, int], sp.Expr]
    remainder: sp.Expr
    remainder_leading_grades: sp.Set
    window_clean: Mapping[tuple[int, int], sp.Expr]


def _normalise(expression: sp.Expr) -> sp.Expr:
    """Return one exact rational normal form, including the exact zero."""

    # ``factor_terms`` is an exact structural compaction.  It keeps SymPy from
    # promoting the many non-bookkeeper atoms in an expanded routed residual
    # to a prohibitively large dense multivariate polynomial during cancel.
    return sp.cancel(sp.together(sp.factor_terms(expression)))


def _at_origin(expression: sp.Expr, a: int, b: int) -> sp.Expr:
    return expression.diff(eta_bg, a).diff(sigma_W, b).subs(
        {eta_bg: sp.S.Zero, sigma_W: sp.S.Zero}
    )


def _coefficient(expression: sp.Expr, a: int, b: int) -> sp.Expr:
    coefficient = _at_origin(expression, a, b) / (
        sp.factorial(a) * sp.factorial(b)
    )
    if not coefficient.free_symbols.isdisjoint({eta_bg, sigma_W}):
        raise C.InputError(
            f"coefficient at grade {a},{b} retains a background bookkeeper"
        )
    return coefficient


def _window_polynomial(coefficients: Mapping[tuple[int, int], sp.Expr]) -> sp.Expr:
    return sp.Add(
        *(
            coefficients[(a, b)] * eta_bg**a * sigma_W**b
            for a in range(N + 1)
            for b in range(N + 1)
        )
    )


def _remainder_leading_grades(remainder: sp.Expr) -> sp.Set:
    if remainder == 0:
        return sp.S.EmptySet

    # Traverse total order increasingly.  Once a nonempty order is found it is
    # exactly the requested lowest-total-order set inside [0, 2N]^2.
    for total_order in range(4 * N + 1):
        grades: list[sp.Tuple] = []
        for a in range(2 * N + 1):
            b = total_order - a
            if not 0 <= b <= 2 * N:
                continue
            derivative = _normalise(_at_origin(remainder, a, b))
            if derivative != 0:
                grades.append(sp.Tuple(sp.Integer(a), sp.Integer(b)))
        if grades:
            return sp.FiniteSet(*grades)
    return sp.S.EmptySet


def _grade_leaf(expression: sp.Expr) -> LeafMultigrade:
    if not isinstance(expression, sp.Expr):
        raise C.InputError(
            f"multigrade leaf is not a scalar SymPy expression: {type(expression).__name__}"
        )

    coefficients = {
        (a, b): _coefficient(expression, a, b)
        for a in range(N + 1)
        for b in range(N + 1)
    }
    polynomial = _window_polynomial(coefficients)
    difference = expression - polynomial
    if expression.is_polynomial(eta_bg, sigma_W):
        difference = sp.expand(difference)
    remainder = _normalise(difference)
    # No RECONSTRUCTION guard: R := _normalise(leaf - window), so
    # (leaf - window - R) is identically zero and tests nothing physical (both
    # build legs flagged it tautological).  WINDOW_CLEAN below is the genuine
    # check — it proves the extracted window coefficients captured all in-window
    # content, and it goes nonzero under every coefficient corruption.
    window_clean = {
        (a, b): _normalise(_at_origin(remainder, a, b))
        for a in range(N + 1)
        for b in range(N + 1)
    }
    return LeafMultigrade(
        coefficients=coefficients,
        remainder=remainder,
        remainder_leading_grades=_remainder_leading_grades(remainder),
        window_clean=window_clean,
    )


def _text_key(value: str) -> C.TextAtom:
    return C.TextAtom(value)


def _grade_association(values: Mapping[tuple[int, int], sp.Expr]) -> C.Association:
    return C.Association(
        tuple(
            (_text_key(f"{a},{b}"), values[(a, b)])
            for a in range(N + 1)
            for b in range(N + 1)
        )
    )


def _multigrade_association(
    paths: Iterable[str], multigrades: Mapping[str, LeafMultigrade]
) -> C.Association:
    entries: list[tuple[C.TextAtom, object]] = []
    for path in paths:
        grade = multigrades[path]
        leaf_entries: list[tuple[C.TextAtom, object]] = list(
            _grade_association(grade.coefficients).entries
        )
        leaf_entries.extend(
            (
                (_text_key("REMAINDER"), grade.remainder),
                (
                    _text_key("REMAINDER_LEADING_GRADES"),
                    grade.remainder_leading_grades,
                ),
            )
        )
        entries.append(
            (_text_key(path), C.Association(tuple(leaf_entries)))
        )
    return C.Association(tuple(entries))


def _window_clean_association(
    paths: Iterable[str],
    multigrades_a: Mapping[str, LeafMultigrade],
    multigrades_b: Mapping[str, LeafMultigrade],
    multigrades_residual: Mapping[str, LeafMultigrade],
) -> C.Association:
    return C.Association(
        tuple(
            (
                _text_key(path),
                C.Association(
                    (
                        (
                            _text_key("A"),
                            _grade_association(multigrades_a[path].window_clean),
                        ),
                        (
                            _text_key("B"),
                            _grade_association(multigrades_b[path].window_clean),
                        ),
                        (
                            _text_key("RESIDUAL"),
                            _grade_association(
                                multigrades_residual[path].window_clean
                            ),
                        ),
                    )
                ),
            )
            for path in paths
        )
    )


def _grade_difference_association(
    paths: Iterable[str],
    multigrades_a: Mapping[str, LeafMultigrade],
    multigrades_b: Mapping[str, LeafMultigrade],
    multigrades_residual: Mapping[str, LeafMultigrade],
) -> C.Association:
    entries: list[tuple[C.TextAtom, object]] = []
    for path in paths:
        differences = {
            (a, b): _normalise(
                multigrades_a[path].coefficients[(a, b)]
                - multigrades_b[path].coefficients[(a, b)]
                - multigrades_residual[path].coefficients[(a, b)]
            )
            for a in range(N + 1)
            for b in range(N + 1)
        }
        entries.append((_text_key(path), _grade_association(differences)))
    return C.Association(tuple(entries))


def _remainder_difference_association(
    paths: Iterable[str],
    multigrades_a: Mapping[str, LeafMultigrade],
    multigrades_b: Mapping[str, LeafMultigrade],
    multigrades_residual: Mapping[str, LeafMultigrade],
) -> C.Association:
    return C.Association(
        tuple(
            (
                _text_key(path),
                _normalise(
                    multigrades_a[path].remainder
                    - multigrades_b[path].remainder
                    - multigrades_residual[path].remainder
                ),
            )
            for path in paths
        )
    )


def _case_object(case: A.JoinedCase) -> C.Association:
    return C.Association(
        (
            ("FAMILY", C.TextAtom(case.family)),
            (
                "KEY",
                C.Association(
                    tuple((axis, C.TextAtom(value)) for axis, value in case.key)
                ),
            ),
        )
    )


def _emit(name: str, payload: object) -> None:
    print(f"{name} = {C.serialise(payload)}", flush=True)


def _aligned_leaves(
    case: A.JoinedCase,
) -> tuple[
    tuple[str, ...],
    Mapping[str, sp.Expr],
    Mapping[str, sp.Expr],
    Mapping[str, sp.Expr],
]:
    active_names = dict(A.WL_TO_PY_RENAME)
    pre_a = A.transform(
        case.operand_a,
        active_names,
        bridge_a=True,
        bridge_d=False,
        collapse=None,
    )
    pre_b = A.transform(
        case.operand_b,
        active_names,
        bridge_a=True,
        bridge_d=False,
        collapse=None,
    )
    if pre_a is None or pre_b is None:
        raise C.InputError(f"missing aligned operand for {case.rendered_identity}")
    left = A._bridge_d(pre_a)
    right = A._bridge_d(pre_b)

    if C.key_dict(case.key).get("OBJECT") == "KINETIC":
        pairs = A._kinetic_pairs(case.family, case.key, left, right)
        if pairs is None:
            raise C.InputError(
                f"kinetic adapter rejected {case.rendered_identity}"
            )
        residual = A._kinetic_residual(pairs)
        residual_items = residual.as_dict()
        paths = tuple(label for label, _left, _right in pairs)
        if tuple(label for label, _item in residual.entries) != paths:
            raise C.InputError(
                f"kinetic residual labels differ for {case.rendered_identity}"
            )
        leaves_a = {label: item_a for label, item_a, _item_b in pairs}
        leaves_b = {label: item_b for label, _item_a, item_b in pairs}
        leaves_residual = {label: residual_items[label] for label in paths}
    else:
        if not isinstance(left, sp.Expr) or not isinstance(right, sp.Expr):
            raise C.InputError(
                f"scalar target has a non-scalar operand: {case.rendered_identity}"
            )
        residual = A._arithmetic_residual(left, right)
        if not isinstance(residual, sp.Expr):
            raise C.InputError(
                f"scalar target has a non-scalar residual: {case.rendered_identity}"
            )
        paths = ("ROOT",)
        leaves_a = {"ROOT": left}
        leaves_b = {"ROOT": right}
        leaves_residual = {"ROOT": residual}

    if not (
        set(paths) == set(leaves_a) == set(leaves_b) == set(leaves_residual)
    ):
        raise C.InputError(f"leaf-path mismatch for {case.rendered_identity}")
    return paths, leaves_a, leaves_b, leaves_residual


def _emit_case(case: A.JoinedCase) -> None:
    paths, leaves_a, leaves_b, leaves_residual = _aligned_leaves(case)
    multigrades_a = {path: _grade_leaf(leaves_a[path]) for path in paths}
    multigrades_b = {path: _grade_leaf(leaves_b[path]) for path in paths}
    multigrades_residual = {
        path: _grade_leaf(leaves_residual[path]) for path in paths
    }

    _emit("CASE", _case_object(case))
    _emit("MULTIGRADE_A", _multigrade_association(paths, multigrades_a))
    _emit("MULTIGRADE_B", _multigrade_association(paths, multigrades_b))
    _emit(
        "MULTIGRADE_RESIDUAL",
        _multigrade_association(paths, multigrades_residual),
    )
    _emit(
        "WINDOW_CLEAN",
        _window_clean_association(
            paths, multigrades_a, multigrades_b, multigrades_residual
        ),
    )
    _emit(
        "GRADE_DIFFERENCE",
        _grade_difference_association(
            paths, multigrades_a, multigrades_b, multigrades_residual
        ),
    )
    _emit(
        "REMAINDER_DIFFERENCE",
        _remainder_difference_association(
            paths, multigrades_a, multigrades_b, multigrades_residual
        ),
    )


def run() -> int:
    try:
        _emit("A_MODULE_PATH", C.TextAtom(str(Path(A.__file__).resolve())))
        _emit("C_MODULE_PATH", C.TextAtom(str(Path(C.__file__).resolve())))
        _emit("P_MODULE_PATH", C.TextAtom(str(Path(P.__file__).resolve())))
        _emit("PY_TRANSCRIPT_PATH", C.TextAtom(str(C.DEFAULT_PY.resolve())))
        _emit("WL_TRANSCRIPT_PATH", C.TextAtom(str(C.DEFAULT_WL.resolve())))
        _emit("GRADE_WINDOW_N", sp.Integer(N))

        py_tags = C.load_py(C.DEFAULT_PY)
        wl_tags = C.load_wl(C.DEFAULT_WL)
        state = A.RunState(
            Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), []
        )
        seen: Counter[tuple[str, C.Key]] = Counter()

        for family in TARGET_FAMILIES:
            for case in A._family_cases(family, py_tags, wl_tags, state):
                if case.identity not in TARGET_IDENTITIES:
                    continue
                seen[case.identity] += 1
                _emit_case(case)

        expected = Counter({identity: 1 for identity in TARGET_IDENTITIES})
        if seen != expected:
            missing = expected - seen
            extra = seen - expected
            raise C.InputError(
                f"target case multiset mismatch missing={missing} extra={extra}"
            )
        _emit("EMITTED_CASES", sp.Integer(sum(seen.values())))
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
