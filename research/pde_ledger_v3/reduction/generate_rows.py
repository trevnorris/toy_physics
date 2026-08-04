#!/usr/bin/env python3
"""Generate cross-engine-gated quantity-row proposals.

This is a post-audit text consumer.  It never imports either audit program, and
the audit programs must never import it.  Tag correspondence and registry QIDs
are literal data below: no tag, symbol, or QID is inferred from spelling.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import jsonschema
import yaml


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]

DEFAULT_WL_OUTPUT = Path("/home/trevnorris/s11bB_build/wl_run5.txt")
DEFAULT_PY_OUTPUT = Path("/home/trevnorris/s11bB_build/py_run3.txt")
DEFAULT_SCHEMA = HERE / "registry_schema.yaml"
DEFAULT_QUANTITIES = HERE / "quantities.yaml"
DEFAULT_RELATIONS = HERE / "relations.yaml"
DEFAULT_OUTPUT_DIR = HERE / "_generated"
DEFAULT_WL_SOURCE = (
    PROJECT_ROOT
    / "research/pde_ledger_v3/mathematica/"
    / "S11bB_interface_assembly_mathematica_audit.wl"
)
DEFAULT_PY_SOURCE = (
    PROJECT_ROOT
    / "research/pde_ledger_v3/scripts/"
    / "S11bB_interface_assembly_sympy_audit.py"
)

# The existing brane row supplies only the registry-wide row context (scope,
# regime, state, counting axis, and stage id).  QIDs, symbols, kinds, aliases,
# dimensions, and source loci are not copied from it.
ROW_CONTEXT_QID = "Q.brane.B_comp"
STAGE_USES_SHARED_DIMENSIONS_MODULE = False

TAG_START = re.compile(r"^(WL_[A-Za-z0-9_]+|S11BB_[A-Za-z0-9_]+):(?:\s?(.*))?$")
DIMENSION_TRIPLE = re.compile(
    r"^\s*[({]\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*[)}]\s*$"
)
ROUTE_KIND = re.compile(r"^\s*(independent|definitional)(?::|\s)")
LABEL = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class SourceLines:
    dimension: int
    route: int
    token: str


@dataclass(frozen=True)
class QuantitySpec:
    name: str
    qid: str
    symbol: str
    kind: str
    aliases: tuple[str, ...]
    wl_dimension_tag: str | None
    wl_route_tag: str | None
    py_dimension_tag: str | None
    py_route_tag: str | None
    wl_source: SourceLines | None
    py_source: SourceLines | None


@dataclass(frozen=True)
class TagRecord:
    value: str
    line_start: int
    line_end: int


@dataclass(frozen=True)
class BlockedQuantity:
    name: str
    reason: str
    details: tuple[str, ...]


def qspec(
    name: str,
    qid: str,
    symbol: str,
    kind: str,
    aliases: Sequence[str],
    wl_name: str | None,
    py_name: str | None,
    wl_lines: tuple[int, int] | None,
    py_lines: tuple[int, int] | None,
) -> QuantitySpec:
    """Build one literal correspondence entry without deriving any name."""

    return QuantitySpec(
        name=name,
        qid=qid,
        symbol=symbol,
        kind=kind,
        aliases=tuple(aliases),
        wl_dimension_tag=None if wl_name is None else f"WL_DIM_{wl_name}",
        wl_route_tag=None if wl_name is None else f"WL_DIM_ROUTE_KIND_{wl_name}",
        py_dimension_tag=None if py_name is None else f"S11BB_DIM_{py_name}",
        py_route_tag=None if py_name is None else f"S11BB_DIM_ROUTE_KIND_{py_name}",
        wl_source=(
            None
            if wl_name is None or wl_lines is None
            else SourceLines(wl_lines[0], wl_lines[1], f'"{wl_name}"')
        ),
        py_source=(
            None
            if py_name is None or py_lines is None
            else SourceLines(py_lines[0], py_lines[1], f'"{py_name}"')
        ),
    )


# EXPLICIT QUANTITY CORRESPONDENCE TABLE.  Every differing spelling is stated.
# The helper only adds the two fixed, engine-specific tag wrappers; it never
# transforms or compares the quantity names.
QUANTITY_SPECS: tuple[QuantitySpec, ...] = (
    qspec("B_rho", "Q.brane.B_rho", "B_rho", "parameter", ("B_RHO",),
          "B_RHO", "B_rho", (291, 304), (265, 501)),
    qspec("B_rho3", "Q.brane.B_rho3", "B_rho3", "parameter", ("b3", "B_RHO3"),
          "B_RHO3", "B_rho3", (291, 305), (265, 502)),
    qspec("mu_W", "Q.brane.mu_W", "mu_W", "parameter", ("muW", "MU_W"),
          "MU_W", "mu_W", (291, 306), (266, 503)),
    qspec("k_W", "Q.brane.k_W", "k_W", "parameter", ("kW", "K_W"),
          "K_W", "k_W", (292, 307), (266, 504)),
    qspec("kappa_W", "Q.brane.kappa_W", "kappa_W", "parameter", ("kapW", "KAPPA_W"),
          "KAPPA_W", "kappa_W", (292, 308), (267, 505)),
    qspec("C", "Q.brane.C", "C", "parameter", ("cc",),
          "C", "C", (292, 309), (267, 506)),
    qspec(
        "B3_response_deltaW_over_force",
        "Q.brane.B3_response_deltaW_over_force",
        "B3_response_deltaW_over_force",
        "intermediate",
        ("B3_response_deltaW_over_pressure",),
        "B3_RESPONSE_DELTAW_OVER_PRESSURE",
        "B3_response_deltaW_over_force",
        (293, 310),
        (274, 507),
    ),
    qspec(
        "B4_response_stress_over_divu",
        "Q.brane.B4_response_stress_over_divu",
        "B4_response_stress_over_divu",
        "intermediate",
        ("B4_response_Pi_over_minus_divu",),
        "B4_RESPONSE_PI_OVER_MINUS_DIVU",
        "B4_response_stress_over_divu",
        (293, 311),
        (275, 508),
    ),
    qspec("B6_coefficient", "Q.brane.B6_coefficient", "B6_coefficient", "intermediate", (),
          "B6_COEFFICIENT", "B6_coefficient", (294, 312), (533, 534)),
    qspec("K_L", "Q.brane.K_L", "K_L", "parameter", ("kD", "K_D"),
          "K_D", "K_L", (295, 313), (276, 520)),
    qspec("D_theta", "Q.brane.D_theta", "D_theta", "parameter", ("aTheta", "A_THETA"),
          "A_THETA", "D_theta", (295, 313), (276, 521)),
    qspec("D_e", "Q.brane.D_e", "D_e", "parameter", ("aE", "A_E"),
          "A_E", "D_e", (295, 314), (277, 522)),
    qspec("A_theta", "Q.brane.A_theta", "A_theta", "parameter", ("gTheta", "G_THETA"),
          "G_THETA", "A_theta", (296, 314), (277, 523)),
    qspec("A_theta_e", "Q.brane.A_theta_e", "A_theta_e", "parameter", ("gThetaE", "G_THETAE"),
          "G_THETAE", "A_theta_e", (296, 315), (278, 524)),
    qspec("Lambda_A0", "Q.interface.Lambda_A0", "Lambda_A0", "parameter", ("la0",),
          "LAMBDA_A0", "Lambda_A0", (297, 316), (269, 509)),
    qspec("Lambda_V0", "Q.interface.Lambda_V0", "Lambda_V0", "parameter", ("lv0",),
          "LAMBDA_V0", "Lambda_V0", (297, 316), (270, 510)),
    qspec("Lambda_X0", "Q.interface.Lambda_X0", "Lambda_X0", "parameter", ("lx0",),
          "LAMBDA_X0", "Lambda_X0", (297, 317), (270, 511)),
    qspec("tau_A", "Q.interface.tau_A", "tau_A", "parameter", ("tauA",),
          "TAU_A", "tau_A", (298, 318), (271, 512)),
    qspec("tau_V", "Q.interface.tau_V", "tau_V", "parameter", ("tauV",),
          "TAU_V", "tau_V", (298, 318), (271, 513)),
    qspec("tau_X", "Q.interface.tau_X", "tau_X", "parameter", ("tauX",),
          "TAU_X", "tau_X", (298, 318), (271, 514)),
    qspec("affinity", "Q.interface.affinity", "affinity", "intermediate", ("A", "AFFINITY"),
          "AFFINITY", "affinity", (299, 319), (269, 515)),
    qspec("mu_theta", "Q.brane.mu_theta", "mu_theta", "intermediate", ("muTheta", "MU_THETA"),
          "MU_THETA", "mu_theta", (299, 319), (268, 516)),
    qspec("mu_s", "Q.interface.mu_s", "mu_s", "intermediate", ("mus", "MU_S"),
          "MU_S", "mu_s", (299, 320), (268, 517)),
    qspec("face_V_coeff", "Q.interface.face_V_coeff", "face_V_coeff", "intermediate", ("pV", "P_V"),
          "FACE_P_OVER_V", "face_V_coeff", (300, 320), (272, 518)),
    qspec(
        "face_mu_theta_coeff",
        "Q.interface.face_mu_theta_coeff",
        "face_mu_theta_coeff",
        "intermediate",
        ("pMuTheta", "P_muTheta"),
        "FACE_P_OVER_MUTHETA",
        "face_mu_theta_coeff",
        (300, 321),
        (273, 519),
    ),
    # SymPy emits a composite summary in addition to the two scalar response
    # coefficients.  Mathematica has no single-triple counterpart, so it must
    # remain PY_ONLY; it is explicit here and is never decomposed heuristically.
    qspec(
        "face_response",
        "Q.interface.face_response",
        "face_response",
        "intermediate",
        (),
        None,
        "face_response",
        None,
        (535, 537),
    ),
)


# EXPLICIT NON-QUANTITY/RELATION-CANDIDATE CORRESPONDENCE TABLE.  Values in
# these records are deliberately not compared; the report puts them side by
# side for review.  Exact pairs are listed even where the suffixes happen to
# look alike.
RELATION_TAG_PAIRS: tuple[tuple[str, str], ...] = (
    ("WL_ENERGY_BASIS", "S11BB_ENERGY_BASIS"),
    ("WL_ENERGY_BASIS_INDEPENDENT_TERMS", "S11BB_ENERGY_BASIS_INDEPENDENT_TERMS"),
    ("WL_ENERGY_BASIS_OMISSIONS", "S11BB_ENERGY_BASIS_OMISSIONS"),
    ("WL_BASIS_REDUNDANCY_UNDER_CONSTRAINT", "S11BB_BASIS_REDUNDANCY_UNDER_CONSTRAINT"),
    ("WL_FACE_RESPONSE", "S11BB_FACE_RESPONSE"),
    ("WL_FACE_RESPONSE_MU_COEFF", "S11BB_FACE_RESPONSE_MU_COEFF"),
    ("WL_ZPERM_REDUCTION_CHECK", "S11BB_ZPERM_REDUCTION_CHECK"),
    ("WL_TWO_PORT_POWER_IDENTITY", "S11BB_TWO_PORT_POWER_IDENTITY"),
    ("WL_PORT_DISSIPATIVITY", "S11BB_PORT_DISSIPATIVITY"),
    ("WL_PORT_CONDITION_KIND", "S11BB_PORT_CONDITION_KIND"),
    ("WL_COEFFICIENT_ADMISSIBILITY", "S11BB_COEFFICIENT_ADMISSIBILITY"),
    ("WL_ONSAGER_CONDITION", "S11BB_ONSAGER_CONDITION"),
    ("WL_ONSAGER_RECIPROCITY", "S11BB_ONSAGER_RECIPROCITY"),
    ("WL_ONSAGER_DETERMINABLE", "S11BB_ONSAGER_DETERMINABLE"),
    ("WL_RELAXATION_TIME_RELATIONS", "S11BB_RELAXATION_TIME_RELATIONS"),
    ("WL_GROWTH_INSIDE_ADMISSIBLE", "S11BB_GROWTH_INSIDE_ADMISSIBLE"),
    ("WL_DECAY_INSIDE_ADMISSIBLE", "S11BB_DECAY_INSIDE_ADMISSIBLE"),
    ("WL_CONSTRAINT", "S11BB_CONSTRAINT"),
    ("WL_CONSTRAINT_TERM_ORIGINS", "S11BB_CONSTRAINT_TERM_ORIGINS"),
    ("WL_INTERNAL_DOF_COUNT", "S11BB_INTERNAL_DOF_COUNT"),
    ("WL_DOF_COUNTING_CONVENTION", "S11BB_DOF_COUNTING_CONVENTION"),
    ("WL_INPLANE_EOM", "S11BB_INPLANE_EOM"),
    ("WL_THICKNESS_EOM", "S11BB_THICKNESS_EOM"),
    ("WL_BULK_FORCE_ON_THICKNESS", "S11BB_BULK_FORCE_ON_THICKNESS"),
    ("WL_RECIPROCAL_TRACTION_THICKNESS_EFFECT", "S11BB_RECIPROCAL_TRACTION_THICKNESS_EFFECT"),
    ("WL_THICKNESS_RESPONSE", "S11BB_THICKNESS_RESPONSE"),
    ("WL_RESPONSE_NORMALIZATION", "S11BB_RESPONSE_NORMALIZATION"),
    ("WL_BULK_OPERATOR_BY_REGIME", "S11BB_BULK_OPERATOR_BY_REGIME"),
    ("WL_MASS_INTERPRETATION_VALID_WHERE", "S11BB_MASS_INTERPRETATION_VALID_WHERE"),
    ("WL_COMPRESSIONAL_RESPONSE", "S11BB_COMPRESSIONAL_RESPONSE"),
    ("WL_LIMITS_AND_PATH", "S11BB_LIMITS_AND_PATH"),
    ("WL_FROZEN_THICKNESS_IDENTIFICATION", "S11BB_FROZEN_THICKNESS_IDENTIFICATION"),
    ("WL_LONGITUDINAL_DISPERSION", "S11BB_LONGITUDINAL_DISPERSION"),
    ("WL_ROOTS", "S11BB_ROOTS"),
    ("WL_ROOT_STABILITY_CLASS", "S11BB_ROOT_STABILITY_CLASS"),
    ("WL_STABILITY_CONDITION", "S11BB_STABILITY_CONDITION"),
    ("WL_IMAGINARY_PART", "S11BB_IMAGINARY_PART"),
    ("WL_DISSIPATION_ORIGIN", "S11BB_DISSIPATION_ORIGIN"),
    ("WL_RECIPROCAL_TRACTION_ROOT_EFFECT", "S11BB_RECIPROCAL_TRACTION_ROOT_EFFECT"),
    ("WL_BRANCH_REALAXIS_CHECK", "S11BB_BRANCH_REALAXIS_CHECK"),
    ("WL_BRANCH_DEGENERATE_POINT", "S11BB_BRANCH_DEGENERATE_POINT"),
    ("WL_BRANCH_SENSITIVITY", "S11BB_BRANCH_SENSITIVITY"),
    ("WL_SHEET_OF_EACH_ROOT", "S11BB_SHEET_OF_EACH_ROOT"),
    ("WL_KERNEL_ORIENTATION_IDENTITIES", "S11BB_KERNEL_ORIENTATION_IDENTITIES"),
    ("WL_KERNEL_PROPAGATION_RESIDUALS", "S11BB_KERNEL_PROPAGATION_RESIDUALS"),
    ("WL_KERNEL_POLE_LOCATIONS", "S11BB_KERNEL_POLE_LOCATIONS"),
    ("WL_CAUSALITY_CHECK", "S11BB_CAUSALITY_CHECK"),
    ("WL_GROWTH_ARTIFACT_DIAGNOSTICS", "S11BB_GROWTH_ARTIFACT_DIAGNOSTICS"),
    ("WL_DECAY_ARTIFACT_DIAGNOSTICS", "S11BB_DECAY_ARTIFACT_DIAGNOSTICS"),
    ("WL_CONVENTION_CHECK_INPLANE", "S11BB_CONVENTION_CHECK_INPLANE"),
    ("WL_CONVENTION_CHECK_CONSERVATIVE", "S11BB_CONVENTION_CHECK_CONSERVATIVE"),
    ("WL_CONSERVATIVE_POSITIVITY_INEQUALITY", "S11BB_CONSERVATIVE_POSITIVITY_INEQUALITY"),
    ("WL_ENERGY_SINKS", "S11BB_ENERGY_SINKS"),
    ("WL_ENERGY_SOURCES", "S11BB_ENERGY_SOURCES"),
    ("WL_UNATTRIBUTED_SINK_TERMS", "S11BB_UNATTRIBUTED_SINK_TERMS"),
    ("WL_UNATTRIBUTED_EXCHANGE_TERMS", "S11BB_UNATTRIBUTED_EXCHANGE_TERMS"),
    ("WL_PRESSURE_WORK_SIGN_CHECK", "S11BB_PRESSURE_WORK_SIGN_CHECK"),
    ("WL_FULL_TWO_PORT_BALANCE_CHECK", "S11BB_FULL_TWO_PORT_BALANCE_CHECK"),
    ("WL_TRANSVERSE_COUPLING", "S11BB_TRANSVERSE_COUPLING"),
    ("WL_TRANSVERSE_DISPERSION", "S11BB_TRANSVERSE_DISPERSION"),
    ("WL_TRANSVERSE_DISSIPATION", "S11BB_TRANSVERSE_DISSIPATION"),
    ("WL_HOMOGENEITY_INPLANE_EQUATION", "S11BB_HOMOGENEITY_INPLANE_EOM"),
    ("WL_HOMOGENEITY_THICKNESS_EQUATION", "S11BB_HOMOGENEITY_THICKNESS_EOM"),
    ("WL_HOMOGENEITY_MASS_BALANCE", "S11BB_HOMOGENEITY_MASS_BALANCE"),
    ("WL_HOMOGENEITY_AFFINITY", "S11BB_HOMOGENEITY_AFFINITY"),
    ("WL_HOMOGENEITY_CLOSURE", "S11BB_HOMOGENEITY_CLOSURE"),
    ("WL_HOMOGENEITY_FACE_RESPONSE", "S11BB_HOMOGENEITY_FACE_RESPONSE"),
    ("WL_HOMOGENEITY_TWO_PORT_POWER_IDENTITY", "S11BB_HOMOGENEITY_TWO_PORT_POWER_IDENTITY"),
    ("WL_HOMOGENEITY_DISPERSION_DETERMINANT", "S11BB_HOMOGENEITY_DISPERSION_DETERMINANT"),
    ("WL_HOMOGENEITY_ABLATION_DEMO", "S11BB_HOMOGENEITY_ABLATION_DEMO"),
    ("WL_CONTROL_NO_THICKNESS", "S11BB_CONTROL_NO_THICKNESS"),
    ("WL_CONTROL_A_ATTRIBUTION", "S11BB_CONTROL_A_ATTRIBUTION"),
    ("WL_CONTROL_NO_GRADIENT_STIFFNESS", "S11BB_CONTROL_NO_GRADIENT_STIFFNESS"),
    ("WL_CONTROL_IMPERMEABLE", "S11BB_CONTROL_IMPERMEABLE"),
    ("WL_CONTROL_NO_CROSS_TERM", "S11BB_CONTROL_NO_CROSS_TERM"),
    ("WL_CONTROL_NO_MU_COUPLING", "S11BB_CONTROL_NO_MU_COUPLING"),
    ("WL_CONTROL_NO_RECIPROCAL_TRACTION", "S11BB_CONTROL_NO_RECIPROCAL_TRACTION"),
    ("WL_CONTROLS_ON_TRANSVERSE", "S11BB_CONTROLS_ON_TRANSVERSE"),
    ("WL_VALIDITY_CONDITIONS", "S11BB_VALIDITY_CONDITIONS"),
    ("WL_VALIDITY_FAILURE_REGION", "S11BB_VALIDITY_FAILURE_REGION"),
)


class GeneratorError(ValueError):
    """Input, mapping, or output-safety failure."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wl-output", type=Path, default=DEFAULT_WL_OUTPUT)
    parser.add_argument("--py-output", type=Path, default=DEFAULT_PY_OUTPUT)
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument("--quantities", type=Path, default=DEFAULT_QUANTITIES)
    parser.add_argument("--relations", type=Path, default=DEFAULT_RELATIONS)
    parser.add_argument("--wl-source", type=Path, default=DEFAULT_WL_SOURCE)
    parser.add_argument("--py-source", type=Path, default=DEFAULT_PY_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--label", default="S11bB")
    args = parser.parse_args(argv)
    if not LABEL.fullmatch(args.label):
        parser.error("--label must contain only letters, digits, dot, underscore, or hyphen")
    return args


def read_yaml(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise GeneratorError(f"cannot read YAML {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise GeneratorError(f"{path}: top level must be a mapping")
    return value


def schema_errors(document: Mapping[str, Any], schema: Mapping[str, Any]) -> list[str]:
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(document), key=lambda item: list(item.absolute_path))
    result: list[str] = []
    for error in errors:
        location = ".".join(str(part) for part in error.absolute_path) or "<root>"
        result.append(f"{location}: {error.message}")
    return result


def require_schema_valid(document: Mapping[str, Any], schema: Mapping[str, Any], label: str) -> None:
    errors = schema_errors(document, schema)
    if errors:
        raise GeneratorError(f"{label} violates registry schema: {errors[0]}")


def parse_tag_file(path: Path, required_prefix: str) -> dict[str, TagRecord]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise GeneratorError(f"cannot read tag output {path}: {exc}") from exc

    records: dict[str, TagRecord] = {}
    active_tag: str | None = None
    active_start = 0
    active_parts: list[str] = []

    def finish(line_end: int) -> None:
        nonlocal active_tag, active_start, active_parts
        if active_tag is None:
            return
        if active_tag in records:
            first = records[active_tag].line_start
            raise GeneratorError(
                f"{path}: duplicate tag {active_tag} at lines {first} and {active_start}"
            )
        records[active_tag] = TagRecord("\n".join(active_parts), active_start, line_end)
        active_tag = None
        active_start = 0
        active_parts = []

    for number, line in enumerate(lines, start=1):
        match = TAG_START.match(line)
        if match:
            tag = match.group(1)
            if tag.startswith(required_prefix):
                finish(number - 1)
                active_tag = tag
                active_start = number
                active_parts = [match.group(2) or ""]
                continue
        if active_tag is not None:
            # A verdict terminates the last record; other non-tag lines are a
            # verbatim continuation (needed for pretty-printed matrices).
            if line.startswith("VERDICT:"):
                finish(number - 1)
            else:
                active_parts.append(line)
    finish(len(lines))
    return records


def parse_dimension(value: str) -> tuple[int, int, int] | None:
    match = DIMENSION_TRIPLE.fullmatch(value)
    if not match:
        return None
    return tuple(int(match.group(index)) for index in (1, 2, 3))  # type: ignore[return-value]


def parse_route_kind(value: str) -> str | None:
    match = ROUTE_KIND.match(value)
    return None if match is None else match.group(1)


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            while chunk := handle.read(1024 * 1024):
                digest.update(chunk)
    except OSError as exc:
        raise GeneratorError(f"cannot hash protected input {path}: {exc}") from exc
    return digest.hexdigest()


def unique_paths(paths: Iterable[Path]) -> tuple[Path, ...]:
    result: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            result.append(resolved)
    return tuple(result)


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(resolved)


def locus(path: Path, line: int) -> dict[str, Any]:
    return {"path": display_path(path), "line_start": line, "line_end": line}


def validate_source_lines(path: Path, source: SourceLines) -> str | None:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        return f"cannot read provenance source {path}: {exc}"
    for label, number in (("dimension", source.dimension), ("route", source.route)):
        if number < 1 or number > len(lines):
            return f"{label} provenance line {number} is outside {path} (1..{len(lines)})"
        if source.token not in lines[number - 1]:
            return (
                f"{label} provenance line {number} in {path} does not contain "
                f"the explicit token {source.token}"
            )
    return None


def explicit_coverage() -> tuple[set[str], set[str]]:
    wl: set[str] = set()
    py: set[str] = set()

    def add(target: set[str], tag: str | None) -> None:
        if tag is None:
            return
        if tag in target:
            raise GeneratorError(f"correspondence table contains duplicate tag {tag}")
        target.add(tag)

    for spec in QUANTITY_SPECS:
        add(wl, spec.wl_dimension_tag)
        add(wl, spec.wl_route_tag)
        add(py, spec.py_dimension_tag)
        add(py, spec.py_route_tag)
    for wl_tag, py_tag in RELATION_TAG_PAIRS:
        add(wl, wl_tag)
        add(py, py_tag)
    return wl, py


def find_context(live_quantities: Mapping[str, Any]) -> Mapping[str, Any]:
    rows = live_quantities.get("quantities")
    if not isinstance(rows, list):
        raise GeneratorError("live quantities document has no quantity list")
    for row in rows:
        if isinstance(row, dict) and row.get("qid") == ROW_CONTEXT_QID:
            return row
    raise GeneratorError(f"live registry lacks row-context quantity {ROW_CONTEXT_QID}")


def make_document_header(live_quantities: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": live_quantities["schema_version"],
        "document_kind": "quantities",
        "dimension_convention": live_quantities["dimension_convention"],
        "active_regime": live_quantities["active_regime"],
    }


def build_row(
    spec: QuantitySpec,
    dimension: tuple[int, int, int],
    context: Mapping[str, Any],
    wl_source_path: Path,
    py_source_path: Path,
) -> dict[str, Any]:
    assert spec.wl_source is not None and spec.py_source is not None
    stage_id = context["dimension"]["provenance"]["stage_id"]
    source_loci = [
        locus(wl_source_path, spec.wl_source.dimension),
        locus(wl_source_path, spec.wl_source.route),
        locus(py_source_path, spec.py_source.dimension),
        locus(py_source_path, spec.py_source.route),
    ]
    return {
        "qid": spec.qid,
        "symbol": spec.symbol,
        "kind": spec.kind,
        "scope": list(context["scope"]),
        "regime": list(context["regime"]),
        "state": context["state"],
        "counting_axis": context["counting_axis"],
        "dimension": {
            "convention": "LTM-exponent-vector-v1",
            "exponents": list(dimension),
            "provenance": {
                "stage_id": stage_id,
                "stage_uses_shared_dimensions_module": STAGE_USES_SHARED_DIMENSIONS_MODULE,
                # The schema has one primary locus.  All four exact dimension
                # and route loci from both engines remain in source_loci.
                "source_locus": dict(source_loci[0]),
            },
        },
        "aliases": list(spec.aliases),
        "source_loci": source_loci,
    }


def missing_side_reason(
    wl_dimension: TagRecord | None,
    wl_route: TagRecord | None,
    py_dimension: TagRecord | None,
    py_route: TagRecord | None,
) -> str:
    wl_complete = wl_dimension is not None and wl_route is not None
    py_complete = py_dimension is not None and py_route is not None
    if wl_complete and not py_complete:
        return "WL_ONLY"
    if py_complete and not wl_complete:
        return "PY_ONLY"
    return "UNMAPPED"


def get_record(records: Mapping[str, TagRecord], tag: str | None) -> TagRecord | None:
    return None if tag is None else records.get(tag)


def compare_quantities(
    wl_records: Mapping[str, TagRecord],
    py_records: Mapping[str, TagRecord],
    schema: Mapping[str, Any],
    live_quantities: Mapping[str, Any],
    wl_source_path: Path,
    py_source_path: Path,
) -> tuple[list[dict[str, Any]], list[BlockedQuantity]]:
    context = find_context(live_quantities)
    existing_qids = {
        str(row.get("qid"))
        for row in live_quantities["quantities"]
        if isinstance(row, dict) and "qid" in row
    }
    seen_candidate_qids: set[str] = set()
    header = make_document_header(live_quantities)
    rows: list[dict[str, Any]] = []
    blocked: list[BlockedQuantity] = []

    for spec in QUANTITY_SPECS:
        wl_dimension = get_record(wl_records, spec.wl_dimension_tag)
        wl_route = get_record(wl_records, spec.wl_route_tag)
        py_dimension = get_record(py_records, spec.py_dimension_tag)
        py_route = get_record(py_records, spec.py_route_tag)

        records = (wl_dimension, wl_route, py_dimension, py_route)
        if any(record is None for record in records):
            missing = []
            for tag, record in (
                (spec.wl_dimension_tag, wl_dimension),
                (spec.wl_route_tag, wl_route),
                (spec.py_dimension_tag, py_dimension),
                (spec.py_route_tag, py_route),
            ):
                if tag is None:
                    missing.append("no explicit counterpart tag")
                elif record is None:
                    missing.append(f"missing field/tag {tag}")
            blocked.append(
                BlockedQuantity(
                    spec.name,
                    missing_side_reason(wl_dimension, wl_route, py_dimension, py_route),
                    tuple(dict.fromkeys(missing)),
                )
            )
            continue

        assert wl_dimension is not None and wl_route is not None
        assert py_dimension is not None and py_route is not None
        wl_triple = parse_dimension(wl_dimension.value)
        py_triple = parse_dimension(py_dimension.value)
        wl_kind = parse_route_kind(wl_route.value)
        py_kind = parse_route_kind(py_route.value)

        disagreement: list[str] = []
        if wl_triple != py_triple:
            disagreement.append(
                f"dimension: WL={wl_dimension.value!r}; PY={py_dimension.value!r}"
            )
        if wl_kind != py_kind:
            disagreement.append(
                f"route kind: WL={wl_kind or wl_route.value!r}; PY={py_kind or py_route.value!r}"
            )
        if disagreement:
            blocked.append(BlockedQuantity(spec.name, "DISAGREEMENT", tuple(disagreement)))
            continue
        if wl_triple is None:
            blocked.append(
                BlockedQuantity(
                    spec.name,
                    "UNMAPPED",
                    (
                        "missing field: neither engine supplied a dimension exponent triple",
                        f"WL={wl_dimension.value!r}",
                        f"PY={py_dimension.value!r}",
                    ),
                )
            )
            continue
        if wl_kind is None:
            blocked.append(
                BlockedQuantity(
                    spec.name,
                    "UNMAPPED",
                    (
                        "missing field: route kind is neither independent nor definitional",
                        f"WL={wl_route.value!r}",
                        f"PY={py_route.value!r}",
                    ),
                )
            )
            continue
        if spec.qid in existing_qids or spec.qid in seen_candidate_qids:
            blocked.append(
                BlockedQuantity(
                    spec.name,
                    "UNMAPPED",
                    (f"qid is not usable because it duplicates {spec.qid}",),
                )
            )
            continue
        if spec.wl_source is None or spec.py_source is None:
            blocked.append(
                BlockedQuantity(spec.name, "UNMAPPED", ("missing provenance source mapping",))
            )
            continue
        provenance_errors = tuple(
            error
            for error in (
                validate_source_lines(wl_source_path, spec.wl_source),
                validate_source_lines(py_source_path, spec.py_source),
            )
            if error is not None
        )
        if provenance_errors:
            blocked.append(BlockedQuantity(spec.name, "UNMAPPED", provenance_errors))
            continue

        row = build_row(spec, wl_triple, context, wl_source_path, py_source_path)
        candidate_document = {**header, "quantities": [row]}
        errors = schema_errors(candidate_document, schema)
        if errors:
            blocked.append(
                BlockedQuantity(
                    spec.name,
                    "UNMAPPED",
                    tuple(f"schema constraint: {error}" for error in errors),
                )
            )
            continue
        rows.append(row)
        seen_candidate_qids.add(spec.qid)

    return rows, blocked


def indented(value: str) -> str:
    if not value:
        return "    [EMPTY]\n"
    return "".join(f"    {line}\n" for line in value.splitlines())


def relation_candidates_markdown(
    wl_records: Mapping[str, TagRecord],
    py_records: Mapping[str, TagRecord],
    wl_output: Path,
    py_output: Path,
) -> tuple[str, list[str]]:
    lines = [
        "# Relation candidates — review required\n",
        "**REQUIRES REVIEW BEFORE INSERTION.** These prose/algebra records were not auto-compared, and no relation YAML rows were generated.\n",
        f"- Mathematica transcript: `{wl_output}`\n",
        f"- SymPy transcript: `{py_output}`\n",
    ]
    missing: list[str] = []
    for index, (wl_tag, py_tag) in enumerate(RELATION_TAG_PAIRS, start=1):
        wl_record = wl_records.get(wl_tag)
        py_record = py_records.get(py_tag)
        if wl_record is None or py_record is None:
            if wl_record is None and py_record is None:
                status = "UNMAPPED"
            elif wl_record is None:
                status = "PY_ONLY"
            else:
                status = "WL_ONLY"
            missing.append(f"{status}: {wl_tag} <-> {py_tag}")
        lines.append(f"\n## Candidate {index}: `{wl_tag}` ↔ `{py_tag}`\n")
        lines.append("**Status: REQUIRES REVIEW BEFORE INSERTION**\n")
        lines.append(f"\n### Mathematica — `{wl_tag}`\n\n")
        lines.append(indented("[MISSING]" if wl_record is None else wl_record.value))
        lines.append(f"\n### SymPy — `{py_tag}`\n\n")
        lines.append(indented("[MISSING]" if py_record is None else py_record.value))
    return "".join(lines), missing


def blocked_markdown(
    emitted_count: int,
    blocked: Sequence[BlockedQuantity],
    unmapped_tags: Sequence[str],
    relation_mapping_failures: Sequence[str],
) -> str:
    counts = Counter(item.reason for item in blocked)
    lines = [
        "# Blocked quantity proposals\n\n",
        f"Emitted: **{emitted_count}**\n\n",
        "Quantity block counts: "
        + ", ".join(
            f"{reason}=**{counts.get(reason, 0)}**"
            for reason in ("DISAGREEMENT", "WL_ONLY", "PY_ONLY", "UNMAPPED")
        )
        + "\n",
    ]
    if blocked:
        lines.append("\n## Quantities not emitted\n")
        for item in blocked:
            lines.append(f"\n### `{item.name}` — {item.reason}\n")
            for detail in item.details:
                lines.append(f"\n- {detail}\n")
    else:
        lines.append("\nNo quantity was blocked.\n")

    lines.append("\n## Output tags absent from the explicit correspondence tables\n")
    if unmapped_tags:
        for tag in unmapped_tags:
            lines.append(f"\n- `UNMAPPED`: `{tag}`\n")
    else:
        lines.append("\nNone.\n")

    lines.append("\n## Explicit relation pairs missing an engine side\n")
    if relation_mapping_failures:
        for failure in relation_mapping_failures:
            lines.append(f"\n- {failure}\n")
    else:
        lines.append("\nNone.\n")
    return "".join(lines)


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        if temporary_name is not None:
            temporary = Path(temporary_name)
            if temporary.exists():
                temporary.unlink()


def refuse_protected_outputs(outputs: Sequence[Path], protected: Sequence[Path]) -> None:
    protected_set = {path.resolve() for path in protected}
    for output in outputs:
        if output.resolve() in protected_set:
            raise GeneratorError(f"refusing to overwrite protected input {output}")


def run(args: argparse.Namespace) -> int:
    protected_paths = unique_paths(
        (
            args.wl_output,
            args.py_output,
            args.schema,
            args.quantities,
            args.relations,
            args.wl_source,
            args.py_source,
        )
    )
    before = {path: file_digest(path) for path in protected_paths}

    schema = read_yaml(args.schema)
    live_quantities = read_yaml(args.quantities)
    live_relations = read_yaml(args.relations)
    require_schema_valid(live_quantities, schema, str(args.quantities))
    require_schema_valid(live_relations, schema, str(args.relations))

    wl_records = parse_tag_file(args.wl_output, "WL_")
    py_records = parse_tag_file(args.py_output, "S11BB_")
    covered_wl, covered_py = explicit_coverage()
    unmapped_tags = sorted(
        (set(wl_records) - covered_wl) | (set(py_records) - covered_py)
    )

    rows, blocked = compare_quantities(
        wl_records,
        py_records,
        schema,
        live_quantities,
        args.wl_source,
        args.py_source,
    )
    relation_text, relation_mapping_failures = relation_candidates_markdown(
        wl_records, py_records, args.wl_output, args.py_output
    )

    quantity_document = {**make_document_header(live_quantities), "quantities": rows}
    if rows:
        errors = schema_errors(quantity_document, schema)
        if errors:
            raise GeneratorError(
                "internal error: combined proposed quantity document violates schema: " + errors[0]
            )
    else:
        # The schema requires at least one quantity.  Never fake a row merely
        # to make an all-blocked run conform.
        raise GeneratorError("all quantities are blocked; schema forbids an empty proposal document")

    quantities_output = args.output_dir / f"quantities_{args.label}.yaml"
    blocked_output = args.output_dir / "BLOCKED.md"
    relations_output = args.output_dir / "RELATION_CANDIDATES.md"
    output_paths = (quantities_output, blocked_output, relations_output)
    refuse_protected_outputs(output_paths, protected_paths)

    quantity_text = yaml.safe_dump(
        quantity_document,
        sort_keys=False,
        allow_unicode=True,
        width=1000,
    )
    blocked_text = blocked_markdown(
        len(rows), blocked, unmapped_tags, relation_mapping_failures
    )
    atomic_write(quantities_output, quantity_text)
    atomic_write(blocked_output, blocked_text)
    atomic_write(relations_output, relation_text)

    after = {path: file_digest(path) for path in protected_paths}
    changed_inputs = [str(path) for path in protected_paths if before[path] != after[path]]
    if changed_inputs:
        raise GeneratorError("protected input changed during run: " + ", ".join(changed_inputs))

    counts = Counter(item.reason for item in blocked)
    print(f"EMITTED: {len(rows)}")
    for reason in ("DISAGREEMENT", "WL_ONLY", "PY_ONLY", "UNMAPPED"):
        print(f"BLOCKED_{reason}: {counts.get(reason, 0)}")
    for item in blocked:
        print(f"BLOCKED: {item.name} — {item.reason}")
    print(f"UNMAPPED_OUTPUT_TAGS: {len(unmapped_tags)}")
    for tag in unmapped_tags:
        print(f"UNMAPPED_TAG: {tag}")
    print(f"QUANTITIES_ARTIFACT: {quantities_output}")
    print(f"BLOCKED_ARTIFACT: {blocked_output}")
    print(f"RELATION_CANDIDATES_ARTIFACT: {relations_output}")
    print("PROTECTED_INPUTS_UNCHANGED: true")

    anything_blocked = bool(blocked or unmapped_tags or relation_mapping_failures)
    return 1 if anything_blocked else 0


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except GeneratorError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
