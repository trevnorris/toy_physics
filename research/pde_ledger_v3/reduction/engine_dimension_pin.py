#!/usr/bin/env python3
"""Compare committed engine-derived dimension laws with registry declarations."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp

from dimension_laws import (
    BoundDimensionLaw,
    SymbolicDimension,
    dimension_residual,
    dimensions_equal,
)
from engine_output_checks import (
    ParsedValue,
    Unparsed,
    normalize_tags,
    parse_tagged_output,
)
from registry_read import Registry, load_registry


HERE = Path(__file__).resolve().parent
LEDGER = Path(os.environ.get("PDE_LEDGER_ENGINE_OUTPUT_ROOT", HERE.parent)).resolve()
REQUIRED_QIDS = frozenset(
    {"Q.brane.rho_br", "Q.brane.mu_R", "Q.brane.B_comp"}
)


@dataclass(frozen=True)
class QuantityRecord:
    qid: str
    tag: str
    extraction: str = "direct"
    coefficient: str | None = None


@dataclass(frozen=True)
class EngineRecord:
    engine: str
    output: Path
    syntax: str
    dimension_symbol: str
    quantities: tuple[QuantityRecord, ...]


ENGINE_RECORDS = (
    EngineRecord(
        "S9-py",
        LEDGER / "scripts/out/S9_light_requires_shear_sympy_audit.out",
        "py",
        "D",
        (
            QuantityRecord("Q.brane.rho_br", "PY_S9_MAIN_DIM_PRIMARY_INERTIA"),
            QuantityRecord("Q.brane.mu_R", "PY_S9_MAIN_DIM_PRIMARY_STIFFNESS"),
        ),
    ),
    EngineRecord(
        "S9-wl",
        LEDGER / "mathematica/out/S9_light_requires_shear_mathematica_audit.out",
        "wl",
        "D",
        (
            QuantityRecord("Q.brane.rho_br", "WL_S9_RHO_DIMENSION"),
            QuantityRecord("Q.brane.mu_R", "WL_S9_MU_DIMENSION"),
        ),
    ),
    EngineRecord(
        "S10-py",
        LEDGER / "scripts/out/S10_brane_mode_spectrum_sympy_audit.out",
        "py",
        "D",
        (
            QuantityRecord(
                "Q.brane.rho_br",
                "PY_S10_LOCAL_REGISTRY_RHO_BR_DERIVED_DIMENSION_SYMBOLIC",
            ),
            QuantityRecord(
                "Q.brane.mu_R",
                "PY_S10_LOCAL_REGISTRY_MU_R_DERIVED_DIMENSION_SYMBOLIC",
            ),
        ),
    ),
    EngineRecord(
        "S10-wl",
        LEDGER / "mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out",
        "wl",
        "braneDimension",
        (
            QuantityRecord(
                "Q.brane.rho_br",
                "WL_S10_MAIN_D2_Q6_INERTIAL_COEFFICIENT_DIMENSIONS",
                "first",
            ),
            QuantityRecord(
                "Q.brane.mu_R",
                "WL_S10_MAIN_D2_Q6_STIFFNESS_COEFFICIENT_DIMENSIONS",
                "first",
            ),
        ),
    ),
    EngineRecord(
        "S11-py",
        LEDGER / "scripts/out/S11_stray_longitudinal_sympy_audit.out",
        "py",
        "D",
        tuple(
            QuantityRecord(
                qid,
                "PY_S11_MAIN_D3_DIM_COEFFICIENTS",
                "coefficient",
                symbol,
            )
            for qid, symbol in (
                ("Q.brane.B_comp", "B_comp"),
                ("Q.brane.mu_R", "mu_R"),
                ("Q.brane.rho_br", "rho_br"),
            )
        ),
    ),
)


@dataclass(frozen=True)
class EngineDimensionObservation:
    engine: str
    qid: str
    tag: str
    engine_operand: SymbolicDimension
    registry_operand: SymbolicDimension
    residual: SymbolicDimension

    @property
    def passed(self) -> bool:
        return dimensions_equal(self.engine_operand, self.registry_operand)


@dataclass(frozen=True)
class EngineDimensionPin:
    observations: tuple[EngineDimensionObservation, ...]
    errors: tuple[str, ...]

    @property
    def covered_qids(self) -> frozenset[str]:
        return frozenset(observation.qid for observation in self.observations)

    @property
    def passed(self) -> bool:
        expected_count = sum(len(record.quantities) for record in ENGINE_RECORDS)
        return (
            not self.errors
            and len(self.observations) == expected_count
            and self.covered_qids == REQUIRED_QIDS
            and all(observation.passed for observation in self.observations)
        )


def _sequence(value: object) -> list[object]:
    if isinstance(value, sp.MatrixBase):
        if value.rows != 1 and value.cols != 1:
            raise ValueError(f"expected a vector, observed matrix shape {value.shape}")
        return list(value)
    if isinstance(value, (tuple, list)):
        return list(value)
    raise ValueError(f"expected a sequence, observed {type(value).__name__}")


def _vector(value: object) -> SymbolicDimension:
    values = _sequence(value)
    if len(values) != 3 or any(isinstance(item, (tuple, list, dict)) for item in values):
        raise ValueError(f"expected three scalar dimension components, observed {value!r}")
    try:
        return tuple(sp.sympify(item) for item in values)
    except sp.SympifyError as exc:
        raise ValueError(f"dimension component is not scalar: {value!r}") from exc


def _extract(value: object, record: QuantityRecord) -> SymbolicDimension:
    if record.extraction == "direct":
        return _vector(value)
    if record.extraction == "first":
        values = _sequence(value)
        if not values:
            raise ValueError("expected a nonempty vector population")
        return _vector(values[0])
    if record.extraction == "coefficient":
        for row in _sequence(value):
            fields = _sequence(row)
            if len(fields) == 2 and str(fields[0]) == record.coefficient:
                return _vector(fields[1])
        raise ValueError(f"coefficient {record.coefficient!r} is absent")
    raise ValueError(f"unknown extraction route {record.extraction!r}")


def _normalise_engine_operand(
    operand: SymbolicDimension,
    engine_dimension_symbol: str,
    law: BoundDimensionLaw,
) -> SymbolicDimension:
    if len(law.bindings) != 1:
        raise ValueError(
            f"pin adapter requires one binding, observed {len(law.bindings)}"
        )
    _parameter, binding_qid = law.bindings[0]
    engine_symbol = sp.Symbol(engine_dimension_symbol)
    canonical_symbol = sp.Symbol(binding_qid, integer=True)
    normalised = tuple(sp.expand(component.subs(engine_symbol, canonical_symbol)) for component in operand)
    unexpected = set().union(*(component.free_symbols for component in normalised)) - {
        canonical_symbol
    }
    if unexpected:
        raise ValueError(f"unmapped engine dimension symbols: {sorted(map(str, unexpected))}")
    return normalised


def inspect_engine_dimension_pin(
    registry: Registry, engine_records: Iterable[EngineRecord] = ENGINE_RECORDS
) -> EngineDimensionPin:
    """Read committed outputs, compute both operands, and retain every residual."""
    observations: list[EngineDimensionObservation] = []
    errors: list[str] = []
    for engine_record in engine_records:
        try:
            tagged = parse_tagged_output(engine_record.output.read_text(encoding="utf-8"))
            selected_tags = {
                record.tag: tagged[record.tag]
                for record in engine_record.quantities
                if record.tag in tagged
            }
            normalised = normalize_tags(selected_tags, engine_record.syntax)
        except (OSError, ValueError) as exc:
            errors.append(f"{engine_record.engine}: {exc}")
            continue
        if tagged.duplicate_tags:
            errors.append(
                f"{engine_record.engine}: duplicate tags {sorted(tagged.duplicate_tags)}"
            )
        for quantity_record in engine_record.quantities:
            try:
                parsed = normalised[quantity_record.tag]
                if isinstance(parsed, Unparsed):
                    raise ValueError(f"unparsed payload: {parsed.error}")
                if not isinstance(parsed, ParsedValue):
                    raise ValueError("normalizer returned no parsed value")
                declaration = registry.dimension_declaration(quantity_record.qid)
                if not isinstance(declaration, BoundDimensionLaw):
                    raise ValueError("registry declaration is not a bound law")
                engine_operand = _normalise_engine_operand(
                    _extract(parsed.value, quantity_record),
                    engine_record.dimension_symbol,
                    declaration,
                )
                registry_operand = declaration.canonical_components
                residual = dimension_residual(engine_operand, registry_operand)
                observations.append(
                    EngineDimensionObservation(
                        engine_record.engine,
                        quantity_record.qid,
                        quantity_record.tag,
                        engine_operand,
                        registry_operand,
                        residual,
                    )
                )
            except (KeyError, ValueError) as exc:
                errors.append(
                    f"{engine_record.engine} {quantity_record.qid} "
                    f"tag={quantity_record.tag}: {exc}"
                )
    return EngineDimensionPin(tuple(observations), tuple(errors))


def _render(values: SymbolicDimension) -> str:
    return "[" + ",".join(sp.sstr(sp.simplify(value)) for value in values) + "]"


def _declared_binding_values(
    registry: Registry, law: BoundDimensionLaw
) -> dict[str, int]:
    values: dict[str, int] = {}
    for _parameter, qid in law.bindings:
        declared = registry.quantities[qid].value
        if declared is None or declared.is_Integer is not True:
            raise ValueError(f"declared dimension binding is not an integer: {qid}")
        values[qid] = int(declared)
    return values


def print_engine_dimension_pin(pin: EngineDimensionPin, registry: Registry) -> bool:
    """Print operands and residuals before their computed guards."""
    for observation in pin.observations:
        label = f"{observation.engine} {observation.qid}"
        print(
            f"ENGINE_DIMENSION_COMPARISON {label} tag={observation.tag}: "
            f"engine-operand={_render(observation.engine_operand)}; "
            f"registry-operand={_render(observation.registry_operand)}; "
            f"residual={_render(observation.residual)}"
        )
        if observation.engine == "S11-py":
            declaration = registry.dimension_declaration(observation.qid)
            if not isinstance(declaration, BoundDimensionLaw):
                raise ValueError("registry declaration is not a bound law")
            numeric_view = tuple(
                sp.Integer(value)
                for value in registry.quantities[observation.qid].dimension
            )
            numeric_reference_residual = dimension_residual(
                observation.engine_operand, numeric_view
            )
            declared_bindings = _declared_binding_values(registry, declaration)
            substitutions = {
                sp.Symbol(qid, integer=True): sp.Integer(value)
                for qid, value in declared_bindings.items()
            }
            residual_at_declared_bindings = tuple(
                sp.expand(component.subs(substitutions))
                for component in numeric_reference_residual
            )
            rendered_bindings = ",".join(
                f"{qid}={value}" for qid, value in sorted(declared_bindings.items())
            )
            print(
                f"S11_NUMERIC_REFERENCE_COMPARISON {label}: "
                f"engine-operand={_render(observation.engine_operand)}; "
                f"declared-reference-operand={_render(numeric_view)}; "
                f"residual={_render(numeric_reference_residual)}; "
                f"declared-bindings={rendered_bindings}; "
                f"residual-at-declared-bindings="
                f"{_render(residual_at_declared_bindings)}"
            )
        print(
            f"ENGINE_DIMENSION_GUARD {label}: "
            f"{'PASS' if observation.passed else 'NONZERO'}"
        )
    for error in pin.errors:
        print(f"ENGINE_DIMENSION_PIN_ERROR: {error}")
    print(f"ENGINE_DIMENSION_PIN: {'PASS' if pin.passed else 'FAIL'}")
    return pin.passed


def main() -> int:
    registry = load_registry()
    pin = inspect_engine_dimension_pin(registry)
    passed = print_engine_dimension_pin(pin, registry)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
