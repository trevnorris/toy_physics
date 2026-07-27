"""Shared exact, axis-labelled dimensions for PDE-ledger audit scripts.

This v0 intentionally contains only the operations required by stages 004
and 011.  A stage binds one ``DimensionBasis`` and constructs every
``Dimension`` through that binding.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import sympy as sp


Exponent = int | float | str | sp.Rational


def _exact(value: Exponent) -> sp.Rational:
    """Return an exact SymPy exponent, including for inputs such as 0.5."""

    return sp.Rational(value)


class DimensionBasis:
    """A stage's one positional-construction and rendering convention."""

    __slots__ = ("_axes", "_render")

    def __init__(self, *axes: str, render: str) -> None:
        if not axes or any(not axis for axis in axes):
            raise ValueError("a dimension basis requires named axes")
        if len(set(axes)) != len(axes):
            raise ValueError(f"dimension basis axes must be unique: {axes}")
        if render not in {"symbolic", "tuple"}:
            raise ValueError(f"unsupported dimension rendering: {render}")
        self._axes = tuple(axes)
        self._render = render

    @property
    def axes(self) -> tuple[str, ...]:
        return self._axes

    def __call__(self, *exponents: Exponent) -> "Dimension":
        if not exponents:
            exponents = (0,) * len(self._axes)
        if len(exponents) != len(self._axes):
            raise ValueError(
                f"basis {self._axes} needs {len(self._axes)} exponents, "
                f"got {len(exponents)}"
            )
        return Dimension(
            self,
            dict(zip(self._axes, (_exact(exponent) for exponent in exponents))),
        )

    def from_mapping(self, exponents: Mapping[str, Exponent]) -> "Dimension":
        missing = set(self._axes) - set(exponents)
        extra = set(exponents) - set(self._axes)
        if missing or extra:
            raise ValueError(
                f"dimension axes do not match basis {self._axes}: "
                f"missing={sorted(missing)}, extra={sorted(extra)}"
            )
        return Dimension(
            self,
            {axis: _exact(exponents[axis]) for axis in self._axes},
        )

    def render(self, dimension: "Dimension") -> str:
        self._require_axes(dimension)
        if self._render == "tuple":
            return "(" + ",".join(sp.sstr(dimension.exponents[axis]) for axis in self._axes) + ")"
        pieces: list[str] = []
        for axis in self._axes:
            exponent = dimension.exponents[axis]
            if exponent == 0:
                continue
            pieces.append(axis if exponent == 1 else f"{axis}^{exponent}")
        return "1" if not pieces else " ".join(pieces)

    def _require_axes(self, dimension: "Dimension") -> None:
        if set(dimension.exponents) != set(self._axes):
            raise ValueError(
                f"dimension axes {tuple(dimension.exponents)} "
                f"do not match basis {self._axes}"
            )


class Dimension:
    """An immutable mapping from axis label to exact exponent."""

    __slots__ = ("_basis", "_exponents")

    def __init__(
        self,
        basis: DimensionBasis,
        exponents: Mapping[str, sp.Rational],
    ) -> None:
        if set(exponents) != set(basis.axes):
            raise ValueError("dimension exponent mapping does not match its basis")
        normalized = {axis: _exact(exponents[axis]) for axis in basis.axes}
        self._basis = basis
        self._exponents = MappingProxyType(normalized)

    @property
    def basis(self) -> DimensionBasis:
        return self._basis

    @property
    def exponents(self) -> Mapping[str, sp.Rational]:
        return self._exponents

    def components(self) -> tuple[sp.Rational, ...]:
        return tuple(self._exponents[axis] for axis in self._basis.axes)

    def without(self, axis: str) -> "Dimension":
        if axis not in self._exponents:
            raise ValueError(f"dimension has no {axis!r} axis")
        values = dict(self._exponents)
        values[axis] = sp.Rational(0)
        return self._basis.from_mapping(values)

    def _combine(self, other: "Dimension", sign: int) -> "Dimension":
        if not isinstance(other, Dimension):
            return NotImplemented
        if set(self._exponents) != set(other._exponents):
            raise ValueError("cannot combine dimensions with different axis sets")
        return self._basis.from_mapping(
            {
                axis: self._exponents[axis] + sign * other._exponents[axis]
                for axis in self._basis.axes
            }
        )

    def __mul__(self, other: "Dimension") -> "Dimension":
        return self._combine(other, 1)

    def __truediv__(self, other: "Dimension") -> "Dimension":
        return self._combine(other, -1)

    def __pow__(self, power: Exponent) -> "Dimension":
        exact_power = _exact(power)
        return self._basis.from_mapping(
            {
                axis: exponent * exact_power
                for axis, exponent in self._exponents.items()
            }
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Dimension):
            return NotImplemented
        return dict(self._exponents) == dict(other._exponents)

    def __hash__(self) -> int:
        return hash(frozenset(self._exponents.items()))

    def __str__(self) -> str:
        return self._basis.render(self)

    def __repr__(self) -> str:
        labelled = ", ".join(
            f"{axis}={exponent}" for axis, exponent in self._exponents.items()
        )
        return f"Dimension({labelled})"


def dim_residual(actual: Dimension, expected: Dimension) -> sp.Expr:
    """Squared residual paired by axis label, never by tuple position."""

    if set(actual.exponents) != set(expected.exponents):
        raise ValueError("cannot compare dimensions with different axis sets")
    return sp.simplify(
        sum(
            (actual.exponents[axis] - expected.exponents[axis]) ** 2
            for axis in actual.exponents
        )
    )


def emit_dimension_sidecar(
    stage_file: str,
    dimensions: Mapping[str, Dimension],
) -> Path:
    """Write deterministic DIM records and source-digest assertions beside a stage."""

    lines: list[str] = []
    declared_axes: tuple[str, ...] | None = None
    for name, dimension in dimensions.items():
        if not name:
            raise ValueError("dimension sidecar quantity names must be nonempty")
        if not isinstance(dimension, Dimension):
            raise TypeError(f"dimension sidecar value for {name!r} is not a Dimension")
        axes = dimension.basis.axes
        if declared_axes is None:
            declared_axes = axes
        elif axes != declared_axes:
            raise ValueError("dimension sidecar records must use one declared axis order")
        axis_text = ",".join(axes)
        exponent_text = ", ".join(
            sp.sstr(dimension.exponents[axis]) for axis in axes
        )
        lines.append(
            f"DIM|axes={axis_text}|name={name}|exponents={{{exponent_text}}}"
        )

    if declared_axes is None:
        raise ValueError("dimension sidecar needs at least one record")
    source = Path(stage_file)
    destination = source.with_suffix(".dimensions.txt")
    source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    ledger_dimensions_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    header = (
        f"DIMENSIONS|stage={source.stem}|axes={','.join(declared_axes)}"
        f"|source_sha256={source_sha256}"
        f"|ledger_dimensions_sha256={ledger_dimensions_sha256}"
    )
    destination.write_text("\n".join([header, *lines]) + "\n", encoding="utf-8")
    return destination
