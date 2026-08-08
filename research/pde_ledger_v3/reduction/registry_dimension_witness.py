#!/usr/bin/env python3
"""Build guarded witnesses from engine dimension emissions to registry declarations.

The manifest classifies every selected source as DERIVED or ECHOED, declares
the exponent multiplier carried by that source, and binds symbolic
specialisation to that source.  ECHOED operands remain visible but can never
produce AGREEMENT.  Axis order is measured from labelled companion emissions
where available; an unlabelled engine vector is reported as UNDETERMINED.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import sympy as sp
import yaml

try:
    from . import engine_output_checks as engine_checks
except ImportError:  # Direct script invocation.
    import engine_output_checks as engine_checks


HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = HERE / "registry_dimension_witness.yaml"
DEFAULT_QUANTITIES = HERE / "quantities.yaml"
STATUS_ORDER = (
    "AGREEMENT",
    "DISAGREEMENT",
    "ECHOED",
    "ECHOED_MISMATCH",
    "BRANCH_DIMENSION_MISMATCH",
    "CONVENTION_MISMATCH",
    "UNDETERMINED",
    "NO_TAGS",
    "NO_ROWS",
    "NOT_EMITTED",
)
GUARDED_STATUSES = frozenset(
    {
        "DISAGREEMENT",
        "ECHOED_MISMATCH",
        "BRANCH_DIMENSION_MISMATCH",
        "CONVENTION_MISMATCH",
        "UNDETERMINED",
        "NO_TAGS",
        "NO_ROWS",
    }
)
SOURCE_CLASSES = frozenset({"DERIVED", "ECHOED"})


class WitnessError(RuntimeError):
    """The witness declaration or an input artifact is not operational."""


Dimension = tuple[sp.Expr, sp.Expr, sp.Expr]


@dataclass(frozen=True)
class RegistryQuantity:
    qid: str
    symbol: str
    convention: str | None
    exponents_raw: object
    value_raw: object


@dataclass(frozen=True)
class Emission:
    tag: str
    raw: Dimension | None
    source_class: str = "DERIVED"
    multiplier: sp.Rational = sp.Rational(1)
    branch_dimension: int | None = None
    specialisations: Mapping[str, object] | None = None
    axis_raw: Dimension | None = None
    axis_tag: str | None = None
    axis_error: str | None = None
    error: str | None = None


@dataclass(frozen=True)
class SpecialisedEmission:
    tag: str
    raw: Dimension | None
    value: Dimension | None
    expected: Dimension | None
    difference: Dimension | None
    applied: tuple[str, ...]
    source_class: str
    multiplier: sp.Rational
    branch_dimension: int | None
    registry_branch_dimension: sp.Rational | None
    axis_raw: Dimension | None
    axis_value: Dimension | None
    axis_tag: str | None
    axis_evidence: str
    error: str | None = None


@dataclass(frozen=True)
class WitnessRow:
    artifact: str
    quantity: RegistryQuantity
    status: str
    emitted_convention_input: str | None
    emissions: tuple[SpecialisedEmission, ...]
    note: str


@dataclass(frozen=True)
class ArtifactReport:
    artifact: str
    rows: tuple[WitnessRow, ...]
    not_emitted: tuple[str, ...]
    selected_quantities: tuple[str, ...]
    status: str


def _read_yaml(path: Path) -> dict[str, Any]:
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise WitnessError(f"cannot load YAML {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise WitnessError(f"{path}: top level must be a mapping")
    return document


def _exact_number(raw: object, label: str) -> sp.Rational:
    if isinstance(raw, bool):
        raise ValueError(f"{label}: boolean is not an exact number")
    if isinstance(raw, int):
        return sp.Rational(raw)
    if (
        isinstance(raw, list)
        and len(raw) == 3
        and raw[0] == "Rat"
        and isinstance(raw[1], int)
        and not isinstance(raw[1], bool)
        and isinstance(raw[2], int)
        and not isinstance(raw[2], bool)
    ):
        return sp.Rational(raw[1], raw[2])
    raise ValueError(f"{label}: expected an integer or [Rat,p,q]")


def _declared_vector(quantity: RegistryQuantity) -> Dimension:
    raw = quantity.exponents_raw
    if not isinstance(raw, list) or len(raw) != 3:
        raise ValueError(f"{quantity.qid}: declared exponents are not a triple")
    values = tuple(
        _exact_number(component, f"{quantity.qid}.dimension.exponents[{index}]")
        for index, component in enumerate(raw)
    )
    return values  # type: ignore[return-value]


def _load_registry_inputs(
    path: Path,
) -> tuple[tuple[RegistryQuantity, ...], str, tuple[str, str, str]]:
    document = _read_yaml(path)
    convention = document.get("dimension_convention")
    if not isinstance(convention, Mapping):
        raise WitnessError(f"{path}: dimension_convention must be a mapping")
    convention_name, raw_bases = convention.get("name"), convention.get("ordered_bases")
    if not isinstance(convention_name, str) or not convention_name:
        raise WitnessError(f"{path}: dimension_convention.name must be a string")
    if (
        not isinstance(raw_bases, list)
        or len(raw_bases) != 3
        or not all(isinstance(base, str) and base for base in raw_bases)
        or len(set(raw_bases)) != 3
    ):
        raise WitnessError(f"{path}: dimension_convention.ordered_bases must be three unique strings")

    rows = document.get("quantities")
    if not isinstance(rows, list) or not rows:
        raise WitnessError(f"{path}: quantities must be a nonempty list")
    result: list[RegistryQuantity] = []
    seen: set[str] = set()
    for index, raw in enumerate(rows):
        if not isinstance(raw, Mapping):
            raise WitnessError(f"{path}: quantities[{index}] must be a mapping")
        qid, symbol = raw.get("qid"), raw.get("symbol")
        if not isinstance(qid, str) or not isinstance(symbol, str):
            raise WitnessError(f"{path}: quantities[{index}] needs string qid and symbol")
        if qid in seen:
            raise WitnessError(f"{path}: duplicate qid {qid}")
        seen.add(qid)
        dimension = raw.get("dimension")
        row_convention: str | None = None
        exponents: object = None
        if isinstance(dimension, Mapping):
            candidate = dimension.get("convention")
            row_convention = candidate if isinstance(candidate, str) and candidate else None
            exponents = dimension.get("exponents")
        result.append(
            RegistryQuantity(qid, symbol, row_convention, exponents, raw.get("value"))
        )
    return (
        tuple(result),
        convention_name,
        (str(raw_bases[0]), str(raw_bases[1]), str(raw_bases[2])),
    )


def load_quantities(path: Path) -> tuple[RegistryQuantity, ...]:
    return _load_registry_inputs(path)[0]


def _snake_to_lower_camel(name: str) -> str:
    pieces = name.split("_")
    return pieces[0] + "".join(piece[:1].upper() + piece[1:] for piece in pieces[1:])


def _emitted_symbol_to_qid(
    quantities: Sequence[RegistryQuantity],
    config: Mapping[str, Any],
    engine: str,
) -> dict[str, str]:
    naming = config.get("symbol_naming")
    style = "canonical"
    exceptions: list[Mapping[str, Any]] = []
    if isinstance(naming, Mapping):
        styles = naming.get("engine_styles")
        if isinstance(styles, Mapping) and isinstance(styles.get(engine), str):
            style = str(styles[engine])
        raw_exceptions = naming.get("exceptions")
        if isinstance(raw_exceptions, list):
            exceptions = [item for item in raw_exceptions if isinstance(item, Mapping)]
    if style not in {"canonical", "lower_camel"}:
        raise WitnessError(f"{engine}: unsupported symbol naming style {style!r}")

    by_canonical = {quantity.symbol: quantity.qid for quantity in quantities}
    result: dict[str, str] = {}
    for canonical, qid in by_canonical.items():
        emitted = canonical if style == "canonical" else _snake_to_lower_camel(canonical)
        if emitted in result and result[emitted] != qid:
            raise WitnessError(f"{engine}: emitted symbol collision {emitted}")
        result[emitted] = qid
    for exception in exceptions:
        canonical, spellings = exception.get("canonical"), exception.get("spellings")
        if not isinstance(canonical, str) or canonical not in by_canonical:
            continue
        if isinstance(spellings, Mapping) and isinstance(spellings.get(engine), str):
            result[str(spellings[engine])] = by_canonical[canonical]
    return result


def _source_identity(source: Mapping[str, Any]) -> str:
    if isinstance(source.get("tag"), str):
        return str(source["tag"])
    if isinstance(source.get("tag_template"), str):
        return str(source["tag_template"])
    tags = source.get("tags")
    if isinstance(tags, list) and all(isinstance(tag, str) for tag in tags):
        return ",".join(tags)
    return "<invalid-source>"


def _source_with_metadata(
    artifact: Mapping[str, Any], source: Mapping[str, Any]
) -> dict[str, Any]:
    result = dict(source)
    metadata = artifact.get("source_metadata", {})
    if not isinstance(metadata, Mapping):
        raise WitnessError(f"{artifact.get('id')}: source_metadata must be a mapping")
    selected = metadata.get(_source_identity(source), {})
    if selected and not isinstance(selected, Mapping):
        raise WitnessError(
            f"{artifact.get('id')}: metadata for {_source_identity(source)} must be a mapping"
        )
    if isinstance(selected, Mapping):
        result.update(selected)
    return result


def _selected_source_tags(
    source: Mapping[str, Any],
    package: str,
    selected_dimensions: frozenset[int] | None,
) -> tuple[tuple[str, int | None], ...]:
    if isinstance(source.get("tag"), str):
        return ((str(source["tag"]), None),)
    if isinstance(source.get("tags"), list):
        tags = source["tags"]
        if not all(isinstance(tag, str) for tag in tags):
            raise WitnessError(f"dimension source for {package}: tags must be strings")
        return tuple((str(tag), None) for tag in tags)
    template, raw_dimensions = source.get("tag_template"), source.get("dimensions")
    if not isinstance(template, str) or not isinstance(raw_dimensions, list):
        raise WitnessError(
            f"dimension source for {package} needs tag, tags, or tag_template+dimensions"
        )
    if not all(
        isinstance(value, int) and not isinstance(value, bool)
        for value in raw_dimensions
    ):
        raise WitnessError(f"dimension source for {package}: dimensions must be integers")
    dimensions = [
        value
        for value in raw_dimensions
        if selected_dimensions is None or value in selected_dimensions
    ]
    return tuple(
        (template.format(package=package, dimension=dimension), dimension)
        for dimension in dimensions
    )


def _rules_mapping(value: object, label: str) -> dict[str, sp.Expr]:
    if isinstance(value, engine_checks.ParsedValue):
        value = value.value
    while isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    if isinstance(value, Mapping):
        pairs = list(value.items())
    elif isinstance(value, (list, tuple)) and all(
        isinstance(item, engine_checks.CasRule) for item in value
    ):
        pairs = [(item.left, item.right) for item in value]
    else:
        raise ValueError(f"{label}: expected one emitted mapping or rule list")
    result: dict[str, sp.Expr] = {}
    for key, raw in pairs:
        if isinstance(raw, bool) or not isinstance(raw, (int, sp.Basic)):
            raise ValueError(f"{label}: component {key} is not a scalar expression")
        result[str(key)] = sp.Integer(raw) if isinstance(raw, int) else raw
    return result


def _axis_measurement(
    source: Mapping[str, Any],
    output: Mapping[str, engine_checks.ParsedValue | engine_checks.Unparsed],
    package: str,
    branch_dimension: int | None,
    ordered_bases: Sequence[str],
) -> tuple[Dimension | None, str | None, str | None]:
    raw_measurement = source.get("axis_measurement")
    if raw_measurement is None:
        return None, None, "NO_AXIS_LABELLED_COMPONENTS"
    if not isinstance(raw_measurement, Mapping):
        return None, None, "axis_measurement must be a mapping"
    raw_tag = raw_measurement.get("tag")
    if raw_tag is None:
        raw_tag = raw_measurement.get("tag_template")
    if not isinstance(raw_tag, str):
        return None, None, "axis_measurement needs tag or tag_template"
    try:
        axis_tag = raw_tag.format(package=package, dimension=branch_dimension)
    except (KeyError, ValueError) as exc:
        return None, raw_tag, f"axis tag expansion failed: {exc}"
    emitted = output.get(axis_tag)
    if emitted is None:
        return None, axis_tag, "MISSING_AXIS_TAG"
    if isinstance(emitted, engine_checks.Unparsed):
        return None, axis_tag, f"UNPARSED_AXIS_TAG {emitted.error}"
    components = raw_measurement.get("components")
    if not isinstance(components, Mapping):
        return None, axis_tag, "axis_measurement.components must be a mapping"
    missing_axes = [base for base in ordered_bases if not isinstance(components.get(base), str)]
    if missing_axes:
        return None, axis_tag, "missing axis component expressions " + ",".join(missing_axes)
    try:
        rules = _rules_mapping(emitted, axis_tag)
        values: list[sp.Expr] = []
        symbolic_names = {name: sp.Symbol(name) for name in rules}
        for base in ordered_bases:
            expression_text = str(components[base])
            symbolic = sp.sympify(expression_text, locals=symbolic_names)
            unknown = sorted(str(symbol) for symbol in symbolic.free_symbols if str(symbol) not in rules)
            if unknown:
                raise ValueError(
                    f"component {base} references absent labels {','.join(unknown)}"
                )
            values.append(sp.simplify(sp.sympify(expression_text, locals=rules)))
    except (TypeError, ValueError, sp.SympifyError) as exc:
        return None, axis_tag, str(exc)
    return tuple(values), axis_tag, None  # type: ignore[return-value]


def _source_emission(
    *,
    artifact: Mapping[str, Any],
    source: Mapping[str, Any],
    tag: str,
    branch_dimension: int | None,
    package: str,
    output: Mapping[str, engine_checks.ParsedValue | engine_checks.Unparsed],
    ordered_bases: Sequence[str],
) -> Emission:
    source_class = source.get("source_class")
    if source_class not in SOURCE_CLASSES:
        raise WitnessError(
            f"{artifact.get('id')}:{_source_identity(source)} needs source_class DERIVED or ECHOED"
        )
    try:
        multiplier = _exact_number(source.get("multiplier", 1), f"{tag}.multiplier")
    except ValueError as exc:
        raise WitnessError(str(exc)) from exc
    specialisations = source.get("specialisations", {})
    if not isinstance(specialisations, Mapping):
        raise WitnessError(f"{tag}: specialisations must be a mapping")
    axis_raw, axis_tag, axis_error = _axis_measurement(
        source, output, package, branch_dimension, ordered_bases
    )
    value = output.get(tag)
    common = {
        "source_class": str(source_class),
        "multiplier": multiplier,
        "branch_dimension": branch_dimension,
        "specialisations": specialisations,
        "axis_raw": axis_raw,
        "axis_tag": axis_tag,
        "axis_error": axis_error,
    }
    if value is None:
        return Emission(tag, None, error="MISSING", **common)
    if isinstance(value, engine_checks.Unparsed):
        return Emission(tag, None, error=f"UNPARSED {value.error}", **common)
    try:
        selector = source.get("select")
        if isinstance(selector, str) and selector.startswith("mapping_value:"):
            name = selector.partition(":")[2]
            if not isinstance(value.value, Mapping):
                raise ValueError("expected an emitted mapping")
            matches = [item for key, item in value.value.items() if str(key) == name]
            if len(matches) != 1:
                raise ValueError(f"expected exactly one mapping value for {name}")
            selected = matches[0]
        else:
            selected = engine_checks._select_shape(value, selector)
        vector = engine_checks._dimension_source_value(
            selected, str(source.get("shape", "vector")), tag
        )
    except ValueError as exc:
        return Emission(tag, None, error=f"SHAPE {exc}", **common)
    return Emission(tag, vector, **common)


def _source_emissions(
    *,
    artifact: Mapping[str, Any],
    config: Mapping[str, Any],
    output: Mapping[str, engine_checks.ParsedValue | engine_checks.Unparsed],
    quantities: Sequence[RegistryQuantity],
    ordered_bases: Sequence[str],
) -> dict[str, list[Emission]]:
    engine = artifact.get("engine")
    if not isinstance(engine, str):
        raise WitnessError("artifact engine must be a string")
    mapping = _emitted_symbol_to_qid(quantities, config, engine)
    selection = artifact.get("selection", {})
    if not isinstance(selection, Mapping):
        raise WitnessError(f"{artifact.get('id')}: selection must be a mapping")
    raw_packages, raw_dimensions = selection.get("packages"), selection.get("dimensions")
    packages = (
        frozenset(raw_packages)
        if isinstance(raw_packages, list) and all(isinstance(value, str) for value in raw_packages)
        else None
    )
    dimensions = (
        frozenset(raw_dimensions)
        if isinstance(raw_dimensions, list)
        and all(isinstance(value, int) and not isinstance(value, bool) for value in raw_dimensions)
        else None
    )
    if raw_packages is not None and packages is None:
        raise WitnessError(f"{artifact.get('id')}: selection.packages must be strings")
    if raw_dimensions is not None and dimensions is None:
        raise WitnessError(f"{artifact.get('id')}: selection.dimensions must be integers")

    result: dict[str, list[Emission]] = {}
    declarations = config.get("dimension_sources", [])
    if not isinstance(declarations, list):
        raise WitnessError(f"{artifact.get('id')}: dimension_sources is not a list")
    for declaration in declarations:
        if not isinstance(declaration, Mapping) or declaration.get("engine") != engine:
            continue
        package, symbols = declaration.get("package"), declaration.get("symbols")
        if not isinstance(package, str) or not isinstance(symbols, Mapping):
            raise WitnessError(f"{artifact.get('id')}: malformed dimension source")
        if packages is not None and package not in packages:
            continue
        for emitted_symbol, raw_source in symbols.items():
            if not isinstance(emitted_symbol, str) or emitted_symbol not in mapping:
                continue
            qid = mapping[emitted_symbol]
            if isinstance(raw_source, str):
                base_source: Mapping[str, Any] = {"tag": raw_source, "shape": "vector"}
            elif isinstance(raw_source, Mapping):
                base_source = raw_source
            else:
                raise WitnessError(
                    f"{artifact.get('id')}:{package}:{emitted_symbol} invalid source declaration"
                )
            source = _source_with_metadata(artifact, base_source)
            result.setdefault(qid, [])
            tags = _selected_source_tags(source, package, dimensions)
            if not tags:
                source_class = source.get("source_class")
                if source_class not in SOURCE_CLASSES:
                    raise WitnessError(
                        f"{artifact.get('id')}:{_source_identity(source)} needs source_class DERIVED or ECHOED"
                    )
                multiplier = _exact_number(
                    source.get("multiplier", 1), f"{_source_identity(source)}.multiplier"
                )
                result[qid].append(
                    Emission(
                        _source_identity(source),
                        None,
                        str(source_class),
                        multiplier,
                        specialisations=source.get("specialisations", {}),
                        axis_error="NO_TAGS",
                        error="NO_TAGS",
                    )
                )
                continue
            for tag, branch_dimension in tags:
                result[qid].append(
                    _source_emission(
                        artifact=artifact,
                        source=source,
                        tag=tag,
                        branch_dimension=branch_dimension,
                        package=package,
                        output=output,
                        ordered_bases=ordered_bases,
                    )
                )

    extras = artifact.get("extra_dimension_sources", {})
    if not isinstance(extras, Mapping):
        raise WitnessError(f"{artifact.get('id')}: extra_dimension_sources must be a mapping")
    known_qids = {quantity.qid for quantity in quantities}
    for qid, raw_sources in extras.items():
        if qid not in known_qids or not isinstance(raw_sources, list):
            raise WitnessError(f"{artifact.get('id')}: invalid extra source for {qid}")
        result.setdefault(str(qid), [])
        for raw_source in raw_sources:
            if not isinstance(raw_source, Mapping):
                raise WitnessError(f"{artifact.get('id')}: malformed extra source for {qid}")
            source = _source_with_metadata(artifact, raw_source)
            tags = _selected_source_tags(source, "EXTRA", dimensions)
            if not tags:
                source_class = source.get("source_class")
                if source_class not in SOURCE_CLASSES:
                    raise WitnessError(
                        f"{artifact.get('id')}:{_source_identity(source)} needs source_class DERIVED or ECHOED"
                    )
                multiplier = _exact_number(
                    source.get("multiplier", 1), f"{_source_identity(source)}.multiplier"
                )
                result[str(qid)].append(
                    Emission(
                        _source_identity(source),
                        None,
                        str(source_class),
                        multiplier,
                        specialisations=source.get("specialisations", {}),
                        axis_error="NO_TAGS",
                        error="NO_TAGS",
                    )
                )
                continue
            for tag, branch_dimension in tags:
                result[str(qid)].append(
                    _source_emission(
                        artifact=artifact,
                        source=source,
                        tag=tag,
                        branch_dimension=branch_dimension,
                        package="EXTRA",
                        output=output,
                        ordered_bases=ordered_bases,
                    )
                )
    return result


def _resolve_specialisations(
    emission: Emission,
    quantities_by_qid: Mapping[str, RegistryQuantity],
) -> tuple[dict[sp.Symbol, sp.Rational], tuple[str, ...], sp.Rational | None, str | None]:
    raw_bindings = emission.specialisations or {}
    free_symbols = sorted(
        {
            str(symbol)
            for vector in (emission.raw, emission.axis_raw)
            if vector is not None
            for component in vector
            for symbol in component.free_symbols
        }
    )
    substitutions: dict[sp.Symbol, sp.Rational] = {}
    applied: list[str] = []
    registry_branch_dimension: sp.Rational | None = None
    unresolved: list[str] = []
    for symbol_name in free_symbols:
        binding = raw_bindings.get(symbol_name)
        if not isinstance(binding, Mapping):
            unresolved.append(symbol_name)
            continue
        registry_qid = binding.get("registry_value")
        branch_qid = binding.get("branch_dimension_from")
        if isinstance(registry_qid, str) and branch_qid is None:
            quantity = quantities_by_qid.get(registry_qid)
            if quantity is None:
                return {}, (), None, f"unknown specialisation quantity {registry_qid}"
            try:
                value = _exact_number(
                    quantity.value_raw, f"specialisation {symbol_name} from {registry_qid}"
                )
            except ValueError as exc:
                return {}, (), None, str(exc)
            applied.append(f"{symbol_name}={value}<-registry({registry_qid})")
        elif isinstance(branch_qid, str) and registry_qid is None:
            quantity = quantities_by_qid.get(branch_qid)
            if quantity is None:
                return {}, (), None, f"unknown branch quantity {branch_qid}"
            try:
                registry_branch_dimension = _exact_number(
                    quantity.value_raw, f"branch registry value {branch_qid}"
                )
            except ValueError as exc:
                return {}, (), None, str(exc)
            if emission.branch_dimension is None:
                return {}, (), registry_branch_dimension, "source has no selected branch dimension"
            value = sp.Rational(emission.branch_dimension)
            applied.append(f"{symbol_name}={value}<-branch({branch_qid})")
        else:
            return {}, (), None, f"invalid source binding for {symbol_name}"
        substitutions[sp.Symbol(symbol_name)] = value
    if unresolved:
        return (
            substitutions,
            tuple(applied),
            registry_branch_dimension,
            "unresolved symbolic exponents " + ",".join(unresolved),
        )
    return substitutions, tuple(applied), registry_branch_dimension, None


def _substitute_vector(vector: Dimension, substitutions: Mapping[sp.Symbol, sp.Expr]) -> Dimension:
    return tuple(sp.simplify(component.subs(substitutions)) for component in vector)  # type: ignore[return-value]


def _specialise_emission(
    emission: Emission,
    declared: Dimension | None,
    quantities_by_qid: Mapping[str, RegistryQuantity],
    ordered_bases: Sequence[str],
) -> SpecialisedEmission:
    if emission.error is not None or emission.raw is None:
        return SpecialisedEmission(
            emission.tag,
            emission.raw,
            None,
            None,
            None,
            (),
            emission.source_class,
            emission.multiplier,
            emission.branch_dimension,
            None,
            emission.axis_raw,
            None,
            emission.axis_tag,
            "UNDETERMINED(" + (emission.axis_error or "missing-vector") + ")",
            emission.error or "missing vector",
        )
    substitutions, applied, registry_branch, binding_error = _resolve_specialisations(
        emission, quantities_by_qid
    )
    if binding_error:
        return SpecialisedEmission(
            emission.tag,
            emission.raw,
            None,
            None,
            None,
            applied,
            emission.source_class,
            emission.multiplier,
            emission.branch_dimension,
            registry_branch,
            emission.axis_raw,
            None,
            emission.axis_tag,
            "UNDETERMINED(" + (emission.axis_error or "specialisation-failed") + ")",
            binding_error,
        )
    value = _substitute_vector(emission.raw, substitutions)
    if any(component.free_symbols for component in value):
        return SpecialisedEmission(
            emission.tag,
            emission.raw,
            None,
            None,
            None,
            applied,
            emission.source_class,
            emission.multiplier,
            emission.branch_dimension,
            registry_branch,
            emission.axis_raw,
            None,
            emission.axis_tag,
            "UNDETERMINED(symbolic-exponent-remains)",
            "symbolic exponent remains after specialisation",
        )

    axis_value: Dimension | None = None
    if emission.axis_error is not None:
        axis_evidence = "UNDETERMINED(" + emission.axis_error + ")"
    elif emission.axis_raw is None:
        axis_evidence = "UNDETERMINED(NO_AXIS_LABELLED_COMPONENTS)"
    else:
        axis_value = _substitute_vector(emission.axis_raw, substitutions)
        if any(component.free_symbols for component in axis_value):
            axis_evidence = "UNDETERMINED(symbolic-axis-component-remains)"
        elif all(
            sp.simplify(left - right) == 0 for left, right in zip(value, axis_value)
        ):
            axis_evidence = "MEASURED[" + ",".join(ordered_bases) + "]"
        else:
            axis_evidence = "MISMATCH[" + ",".join(ordered_bases) + "]"

    expected = None
    difference = None
    if declared is not None:
        expected = tuple(sp.simplify(emission.multiplier * item) for item in declared)
        difference = tuple(
            sp.simplify(emitted - target) for emitted, target in zip(value, expected)
        )
    return SpecialisedEmission(
        emission.tag,
        emission.raw,
        value,
        expected,  # type: ignore[arg-type]
        difference,  # type: ignore[arg-type]
        applied,
        emission.source_class,
        emission.multiplier,
        emission.branch_dimension,
        registry_branch,
        emission.axis_raw,
        axis_value,
        emission.axis_tag,
        axis_evidence,
    )


def _has_nonzero_residual(emissions: Sequence[SpecialisedEmission]) -> bool:
    return any(
        item.difference is not None
        and any(sp.simplify(component) != 0 for component in item.difference)
        for item in emissions
    )


def _build_row(
    artifact_id: str,
    artifact: Mapping[str, Any],
    quantity: RegistryQuantity,
    raw_emissions: Sequence[Emission],
    quantities_by_qid: Mapping[str, RegistryQuantity],
    registry_convention: str = "LTM-exponent-vector-v1",
    ordered_bases: Sequence[str] = ("L", "T", "M"),
) -> WitnessRow:
    convention_input_raw = artifact.get("declared_emitted_convention")
    convention_input = (
        convention_input_raw
        if isinstance(convention_input_raw, str) and convention_input_raw
        else None
    )
    try:
        declared = _declared_vector(quantity)
        declared_error = ""
    except ValueError as exc:
        declared = None
        declared_error = str(exc)
    emissions = tuple(
        _specialise_emission(item, declared, quantities_by_qid, ordered_bases)
        for item in raw_emissions
    )
    errors = [item.error for item in emissions if item.error]
    if declared_error:
        errors.insert(0, declared_error)
    if any(item.error == "NO_TAGS" for item in emissions):
        return WitnessRow(
            artifact_id, quantity, "NO_TAGS", convention_input, emissions, "tag_template expanded to zero tags"
        )
    if errors:
        return WitnessRow(
            artifact_id,
            quantity,
            "UNDETERMINED",
            convention_input,
            emissions,
            " | ".join(str(error) for error in errors),
        )

    mismatched_branches = [
        item
        for item in emissions
        if item.branch_dimension is not None
        and item.registry_branch_dimension is not None
        and sp.Rational(item.branch_dimension) != item.registry_branch_dimension
    ]
    if mismatched_branches:
        return WitnessRow(
            artifact_id,
            quantity,
            "BRANCH_DIMENSION_MISMATCH",
            convention_input,
            emissions,
            "selected branch dimension differs from registry value",
        )

    source_classes = {item.source_class for item in emissions}
    if source_classes == {"ECHOED"}:
        status = "ECHOED_MISMATCH" if _has_nonzero_residual(emissions) else "ECHOED"
        return WitnessRow(
            artifact_id,
            quantity,
            status,
            convention_input,
            emissions,
            "registry/declaration echo; excluded from agreement",
        )
    if source_classes != {"DERIVED"}:
        return WitnessRow(
            artifact_id,
            quantity,
            "UNDETERMINED",
            convention_input,
            emissions,
            "mixed source classes in one artifact/quantity row",
        )

    if quantity.convention is None or convention_input is None:
        missing = []
        if quantity.convention is None:
            missing.append("registry")
        if convention_input is None:
            missing.append("emitted-convention-input")
        return WitnessRow(
            artifact_id,
            quantity,
            "UNDETERMINED",
            convention_input,
            emissions,
            "missing convention input=" + ",".join(missing),
        )
    if quantity.convention != registry_convention:
        return WitnessRow(
            artifact_id,
            quantity,
            "CONVENTION_MISMATCH",
            convention_input,
            emissions,
            "quantity convention differs from registry ordered-bases convention",
        )
    if quantity.convention != convention_input:
        return WitnessRow(
            artifact_id,
            quantity,
            "CONVENTION_MISMATCH",
            convention_input,
            emissions,
            "emitted convention input differs from registry convention input",
        )
    if any(item.axis_evidence.startswith("MISMATCH") for item in emissions):
        return WitnessRow(
            artifact_id,
            quantity,
            "CONVENTION_MISMATCH",
            convention_input,
            emissions,
            "axis-labelled components do not reproduce the bare emitted vector",
        )
    unavailable = [
        item.axis_evidence
        for item in emissions
        if not item.axis_evidence.startswith("MEASURED")
    ]
    if unavailable:
        return WitnessRow(
            artifact_id,
            quantity,
            "UNDETERMINED",
            convention_input,
            emissions,
            "axis order not measured: " + ",".join(sorted(set(unavailable))),
        )
    return WitnessRow(
        artifact_id,
        quantity,
        "DISAGREEMENT" if _has_nonzero_residual(emissions) else "AGREEMENT",
        convention_input,
        emissions,
        "",
    )


def _resolve_path(base: Path, raw: object, label: str) -> Path:
    if not isinstance(raw, str) or not raw:
        raise WitnessError(f"{label} must be a nonempty path string")
    path = Path(raw)
    return path if path.is_absolute() else base / path


def build_report(
    manifest_path: Path,
    quantities_path: Path,
    output_overrides: Mapping[str, Path] | None = None,
) -> tuple[ArtifactReport, ...]:
    manifest = _read_yaml(manifest_path)
    if manifest.get("schema_version") != "registry-dimension-witness-v2":
        raise WitnessError("manifest schema_version must be registry-dimension-witness-v2")
    raw_artifacts = manifest.get("artifacts")
    if not isinstance(raw_artifacts, list) or not raw_artifacts:
        raise WitnessError("manifest artifacts must be a nonempty list")
    quantities, registry_convention, ordered_bases = _load_registry_inputs(quantities_path)
    quantities_by_qid = {quantity.qid: quantity for quantity in quantities}
    base = manifest_path.parent
    reports: list[ArtifactReport] = []
    seen_artifacts: set[str] = set()
    overrides = output_overrides or {}
    for raw_artifact in raw_artifacts:
        if not isinstance(raw_artifact, Mapping):
            raise WitnessError("every artifact must be a mapping")
        artifact_id, engine = raw_artifact.get("id"), raw_artifact.get("engine")
        if not isinstance(artifact_id, str) or not artifact_id or artifact_id in seen_artifacts:
            raise WitnessError(f"invalid or duplicate artifact id {artifact_id!r}")
        if not isinstance(engine, str) or not engine:
            raise WitnessError(f"{artifact_id}: engine must be a nonempty string")
        seen_artifacts.add(artifact_id)
        raw_config = raw_artifact.get("config")
        if raw_config is None:
            config: Mapping[str, Any] = {}
        else:
            config_path = _resolve_path(base, raw_config, f"{artifact_id}.config")
            config = engine_checks.load_config(config_path)
        output_path = overrides.get(artifact_id) or _resolve_path(
            base, raw_artifact.get("output"), f"{artifact_id}.output"
        )
        output = engine_checks.load_output(
            output_path, syntax=engine if engine in {"wl", "py"} else None
        )[1]
        emissions_by_qid = _source_emissions(
            artifact=raw_artifact,
            config=config,
            output=output,
            quantities=quantities,
            ordered_bases=ordered_bases,
        )
        rows: list[WitnessRow] = []
        not_emitted: list[str] = []
        for quantity in quantities:
            if quantity.qid not in emissions_by_qid:
                not_emitted.append(quantity.qid)
                continue
            rows.append(
                _build_row(
                    artifact_id,
                    raw_artifact,
                    quantity,
                    emissions_by_qid[quantity.qid],
                    quantities_by_qid,
                    registry_convention,
                    ordered_bases,
                )
            )
        selected_quantities = tuple(emissions_by_qid)
        reports.append(
            ArtifactReport(
                artifact_id,
                tuple(rows),
                tuple(not_emitted),
                selected_quantities,
                "OK" if rows else "NO_ROWS",
            )
        )
    unknown_overrides = sorted(set(overrides) - seen_artifacts)
    if unknown_overrides:
        raise WitnessError("output override names unknown artifact(s): " + ",".join(unknown_overrides))
    return tuple(reports)


def _expr_text(value: sp.Expr) -> str:
    return sp.sstr(value).replace(" ", "")


def _vector_text(vector: Dimension | None) -> str:
    if vector is None:
        return "NOT_COMPUTED"
    return "[" + ",".join(_expr_text(component) for component in vector) + "]"


def _group_values(
    emissions: Sequence[SpecialisedEmission],
    render: Any,
) -> str:
    groups: dict[str, list[str]] = {}
    for emission in emissions:
        value = str(render(emission)).replace("|", "/")
        groups.setdefault(value, []).append(emission.tag)
    return "{" + ";".join(
        f"{value}@{','.join(tags)}" for value, tags in groups.items()
    ) + "}"


def _row_line(row: WitnessRow) -> str:
    try:
        declared = _declared_vector(row.quantity)
    except ValueError:
        declared = None
    applied = sorted({item for emission in row.emissions for item in emission.applied})
    note = row.note.replace("|", "/") or "none"
    return (
        f"WITNESS|artifact={row.artifact}|quantity={row.quantity.qid}|status={row.status}"
        f"|source_class={_group_values(row.emissions, lambda item: item.source_class)}"
        f"|declared_convention={row.quantity.convention or 'MISSING'}"
        f"|emitted_convention_input={row.emitted_convention_input or 'MISSING'}"
        f"|axis_evidence={_group_values(row.emissions, lambda item: item.axis_evidence)}"
        f"|axis_operand={_group_values(row.emissions, lambda item: _vector_text(item.axis_value))}"
        f"|branch_dimension={_group_values(row.emissions, lambda item: item.branch_dimension if item.branch_dimension is not None else 'none')}"
        f"|registry_branch_dimension={_group_values(row.emissions, lambda item: _expr_text(item.registry_branch_dimension) if item.registry_branch_dimension is not None else 'none')}"
        f"|specialisation={','.join(applied) or 'none-required'}"
        f"|declared={_vector_text(declared)}"
        f"|multiplier={_group_values(row.emissions, lambda item: _expr_text(item.multiplier))}"
        f"|expected(multiplier*declared)={_group_values(row.emissions, lambda item: _vector_text(item.expected))}"
        f"|emitted_raw={_group_values(row.emissions, lambda item: _vector_text(item.raw))}"
        f"|emitted={_group_values(row.emissions, lambda item: _vector_text(item.value))}"
        f"|residual(emitted-multiplier*declared)={_group_values(row.emissions, lambda item: _vector_text(item.difference))}"
        f"|note={note}"
    )


def format_report(reports: Sequence[ArtifactReport]) -> str:
    artifact_scope = "[" + ",".join(report.artifact for report in reports) + "]"
    lines = [
        "REGISTRY_DIMENSION_WITNESS|schema=registry-dimension-witness-v2"
        "|residual=emitted-multiplier*declared"
        f"|scope_artifacts={artifact_scope}"
    ]
    counts: Counter[str] = Counter()
    for report in reports:
        if report.status == "NO_ROWS":
            counts["NO_ROWS"] += 1
        lines.append(
            f"ARTIFACT_COVERAGE|artifact={report.artifact}|status={report.status}"
            f"|selected_quantities=[{','.join(report.selected_quantities)}]"
            f"|row_quantities=[{','.join(row.quantity.qid for row in report.rows)}]"
            f"|not_emitted_count={len(report.not_emitted)}"
            f"|not_emitted_quantities=[{','.join(report.not_emitted)}]"
        )
        for row in report.rows:
            counts[row.status] += 1
            lines.append(_row_line(row))
        counts["NOT_EMITTED"] += len(report.not_emitted)
    guarded = sum(counts[status] for status in GUARDED_STATUSES)
    lines.append(
        f"STATUS_COUNTS|scope_artifacts={artifact_scope}|"
        + "|".join(f"{status}={counts[status]}" for status in STATUS_ORDER)
        + "|guard_statuses=["
        + ",".join(status for status in STATUS_ORDER if status in GUARDED_STATUSES)
        + f"]|guard_count={guarded}"
    )
    return "\n".join(lines)


def _output_overrides(specs: Sequence[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for spec in specs:
        artifact, separator, raw_path = spec.partition("=")
        if not separator or not artifact or not raw_path or artifact in result:
            raise WitnessError(f"invalid or duplicate output override {spec!r}")
        result[artifact] = Path(raw_path)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--quantities", type=Path, default=DEFAULT_QUANTITIES)
    parser.add_argument("--output", action="append", default=[], metavar="ARTIFACT=PATH")
    args = parser.parse_args(argv)
    try:
        reports = build_report(
            args.manifest,
            args.quantities,
            _output_overrides(args.output),
        )
        rendered = format_report(reports)
        print(rendered)
        guarded = sum(
            row.status in GUARDED_STATUSES
            for report in reports
            for row in report.rows
        ) + sum(report.status in GUARDED_STATUSES for report in reports)
        return 1 if guarded else 0
    except (WitnessError, engine_checks.HarnessError, OSError, ValueError) as exc:
        print(f"OPERATIONAL_FAILURE: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
