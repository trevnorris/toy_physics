#!/usr/bin/env python3
"""Ledger stage035 SymPy audit: native source law and parity census.

Standalone, print-only, assert-zero, exact, and file-I/O-free.  The native
source is reconstructed from the translation derivative of an arbitrary
signed mouth profile.  The six-row parity census is independently generated
by multiplicative parity algebra.

Tooth-local runtime ablation uses ``LEDGER_STAGE035_MUTATION``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import os
from typing import Any, Iterable, Mapping

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE035_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

VERDICT = "CONVECTION_LIKE_CONDITIONAL"
R1_SOURCE = "R1_SOURCE_BASIS"
PARITY_BAD = "PARITY_INCONSISTENT"
MAGNITUDE_R1 = "R1(magnitude)"
MANIFEST_DIGEST = "d85de3d8f6b7d615900c8ead24f1589eed3c681e74cf546adb999f2cc7fa5b36"

TOOTH_ORDER = (
    "SOURCE_TRANSLATION_CONTINUITY",
    "SOURCE_NOT_IMPORTED",
    "SOURCE_BASIS",
    "PARITY_RW",
    "PARITY_PW",
    "PARITY_ROTATION",
    "PARITY_TIME_REVERSAL",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
    "UNITS_RESTORED",
    "VERDICT_REDERIVATION",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "SOURCE_TRANSLATION_CONTINUITY":
        "select alpha=2 in the computed continuity residual",
    "SOURCE_NOT_IMPORTED":
        "inject code-native barred symbol Nu into the computed source vector",
    "SOURCE_BASIS":
        "double the computed source before component-wise basis reduction",
    "PARITY_RW":
        "make the isotropic response kernel R_w-odd so derived u_T is R_w-even",
    "PARITY_PW":
        "make the isotropic response kernel P_w-odd so derived u_T is P_w-even",
    "PARITY_ROTATION":
        "make the curl operator axial so derived b_T is polar",
    "PARITY_TIME_REVERSAL":
        "make the active-drain time arrow tau_d T-even",
    "TARGET_BLINDNESS":
        "inject electric A_E into the live source-side dependency expression",
    "DUAL_ENGINE_TERMS":
        "drop q_T_relation from the computed stage term inventory",
    "UNITS_RESTORED":
        "change [q_T] from M T^-1 to M L T^-1 in the expression evaluator",
    "VERDICT_REDERIVATION":
        "select alpha=2 only in the verdict's upstream continuity object",
    "SOURCE_TO_STAGE_MANIFEST":
        "drop one scoped-out source tooth from the canonical 35-row partition",
}


class AuditFailure(AssertionError):
    """A named audit predicate failed."""

    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


def compact(value: Any) -> str:
    try:
        return sp.sstr(sp.factor(sp.cancel(sp.simplify(value))))
    except (TypeError, ValueError, AttributeError):
        return str(value)


def assert_exact(name: str, value: Any) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            assert_exact(f"{name}.{key}", item)
        return
    if isinstance(value, (tuple, list, set, frozenset)):
        for index, item in enumerate(value):
            assert_exact(f"{name}[{index}]", item)
        return
    if isinstance(value, (str, type(None), bool)):
        return
    floats = sp.sympify(value).atoms(sp.Float)
    if floats:
        raise AuditFailure(name, f"machine Float atom(s): {floats}")


def expect_zero(name: str, residual: Any, evidence: Any = None) -> None:
    global PASS_COUNT, FAIL_COUNT
    assert_exact(name, residual)
    clean = sp.simplify(residual)
    assert_exact(name, clean)
    if clean == 0:
        PASS_COUNT += 1
        print(f"PASS  {name}")
        return
    FAIL_COUNT += 1
    print(f"FIRST_FAILURE={name}")
    if ACTIVE_MUTATION == name:
        print(f"FIRED_AT_OWN_ASSERT={name}")
    print(f"FAIL  {name}: residual = {compact(clean)}")
    if evidence is not None:
        print(f"      evidence = {evidence}")
    raise AuditFailure(name, compact(clean))


def expect_bool(name: str, condition: bool, evidence: Any = None) -> None:
    expect_zero(name, sp.Integer(0 if bool(condition) else 1), evidence)


def section(text: str) -> None:
    print("")
    print(text)
    print("-" * len(text))


# ---------------------------------------------------------------------------
# Abstract translated-profile route: eta_a is arbitrary, not a chosen shape.
# ---------------------------------------------------------------------------

xi1, xi2, xi3 = sp.symbols("xi1 xi2 xi3", real=True)
relative_coordinates = (xi1, xi2, xi3)
V1, V2, V3 = sp.symbols("V1 V2 V3", real=True)
body_velocity = sp.Matrix([V1, V2, V3])
s = sp.symbols("s", real=True, nonzero=True)
alpha = sp.symbols("alpha", real=True)
q_T, lambda_T, tau_d = sp.symbols(
    "q_T lambda_T tau_d", real=True, nonzero=True
)
eta_a = sp.Function("eta_a")
profile = eta_a(*relative_coordinates)
signed_density = s * profile

translation_gradient = sp.factor(sum(
    velocity * sp.diff(signed_density, coordinate)
    for velocity, coordinate in zip(body_velocity, relative_coordinates)
))
translation_time_derivative = -translation_gradient
isotropic_flux_divergence = sp.factor(alpha * translation_gradient)
continuity_general = sp.factor(
    translation_time_derivative + isotropic_flux_divergence
)
continuity_coefficient = sp.cancel(continuity_general / translation_gradient)
unique_alpha_solutions = tuple(
    sp.solve(sp.Eq(continuity_coefficient, 0), alpha)
)
unique_alpha = unique_alpha_solutions[0]


Nu, aT, aTp, aL = sp.symbols("Nu aT aTp aL", nonzero=True)
q_A_T = sp.Symbol("q_A_T", nonzero=True)
q_L = sp.Symbol("q_L", nonzero=True)
BARRED_DISPLAY_TO_CODE = {
    "N_u": Nu,
    "a_T": aT,
    "a'_T": aTp,
    "a_L": aL,
    "q_A^T": q_A_T,
    "q_L": q_L,
}
BARRED_SOURCE_MARKERS = frozenset(BARRED_DISPLAY_TO_CODE.values())


@dataclass(frozen=True)
class SourceBuild:
    selected_flux_coefficient: sp.Expr
    continuity_residual: sp.Expr
    unique_flux_coefficient: sp.Expr
    native_flux: sp.Matrix
    native_source: sp.Matrix
    candidate_source: sp.Matrix


def source_build(*, flux_coefficient: sp.Expr = sp.Integer(1),
                 inject_barred: bool = False,
                 source_scale: sp.Expr = sp.Integer(1)) -> SourceBuild:
    """Derive I and J from this route's solved continuity coefficient."""
    native_flux = sp.Matrix([
        unique_alpha * signed_density * velocity for velocity in body_velocity
    ])
    native_source = sp.simplify(q_T) * native_flux
    candidate = source_scale * native_source
    if inject_barred:
        candidate = candidate.copy()
        candidate[0] = sp.expand(candidate[0] + Nu * native_source[0])
    return SourceBuild(
        selected_flux_coefficient=flux_coefficient,
        continuity_residual=sp.factor(
            continuity_general.subs(alpha, flux_coefficient)
        ),
        unique_flux_coefficient=unique_alpha,
        native_flux=native_flux,
        native_source=native_source,
        candidate_source=candidate,
    )


# ---------------------------------------------------------------------------
# Multiplicative parity algebra: all derived rows come from factor products.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Parity:
    rw: int
    pw: int
    rotation: str
    time: int


def rotation_product(kinds: Iterable[str]) -> str:
    vectors = [kind for kind in kinds if kind != "scalar"]
    if not vectors:
        return "scalar"
    if len(vectors) == 1:
        return vectors[0]
    inversion_character = sp.prod(
        -1 if kind == "polar_vector" else 1 for kind in vectors
    )
    return "polar_vector" if inversion_character == -1 else "axial_vector"


def parity_product(*factors: Parity) -> Parity:
    return Parity(
        rw=int(sp.prod(factor.rw for factor in factors)),
        pw=int(sp.prod(factor.pw for factor in factors)),
        rotation=rotation_product(factor.rotation for factor in factors),
        time=int(sp.prod(factor.time for factor in factors)),
    )


def cross_parity(left: Parity, right: Parity) -> Parity:
    if left.rotation == "scalar" or right.rotation == "scalar":
        raise AuditFailure("PARITY_ROTATION", "cross product received a scalar")
    return parity_product(left, right)


@dataclass(frozen=True)
class ParityBuild:
    census: Mapping[str, Parity]
    tau_parity: Parity
    q_parity: Parity


def parity_build(*, response_rw: int = 1, response_pw: int = 1,
                 tau_time: int = -1,
                 curl_rotation: str = "polar_vector") -> ParityBuild:
    s_parity = Parity(-1, -1, "scalar", 1)
    velocity_parity = Parity(1, 1, "polar_vector", -1)
    eta_parity = Parity(1, 1, "scalar", 1)
    lambda_parity = Parity(1, 1, "scalar", 1)
    tau_parity = Parity(1, 1, "scalar", tau_time)
    response_parity = Parity(
        response_rw, response_pw, "scalar", 1
    )
    gradient_parity = Parity(1, 1, curl_rotation, 1)

    q_parity = parity_product(lambda_parity, tau_parity)
    current_parity = parity_product(
        q_parity, s_parity, eta_parity, velocity_parity
    )
    field_parity = parity_product(response_parity, current_parity)
    curl_parity = cross_parity(gradient_parity, field_parity)
    census = {
        "s": s_parity,
        "V": velocity_parity,
        "tau_d/q_T": q_parity,
        "J_T": current_parity,
        "u_T": field_parity,
        "b_T": curl_parity,
    }
    return ParityBuild(census, tau_parity, q_parity)


EXPECTED_CENSUS = {
    "s": Parity(-1, -1, "scalar", 1),
    "V": Parity(1, 1, "polar_vector", -1),
    "tau_d/q_T": Parity(1, 1, "scalar", -1),
    "J_T": Parity(-1, -1, "polar_vector", 1),
    "u_T": Parity(-1, -1, "polar_vector", 1),
    "b_T": Parity(-1, -1, "axial_vector", 1),
}


def parity_invariants(build: ParityBuild) -> dict[str, bool]:
    census = build.census
    order = tuple(EXPECTED_CENSUS)
    return {
        "PARITY_RW": (
            tuple(census[name].rw for name in order)
            == tuple(EXPECTED_CENSUS[name].rw for name in order)
            and census["J_T"].rw == census["u_T"].rw == -1
        ),
        "PARITY_PW": (
            tuple(census[name].pw for name in order)
            == tuple(EXPECTED_CENSUS[name].pw for name in order)
            and census["J_T"].pw == census["u_T"].pw == -1
        ),
        "PARITY_ROTATION": (
            tuple(census[name].rotation for name in order)
            == tuple(EXPECTED_CENSUS[name].rotation for name in order)
            and census["u_T"].rotation == "polar_vector"
            and census["b_T"].rotation == "axial_vector"
        ),
        "PARITY_TIME_REVERSAL": (
            tuple(census[name].time for name in order)
            == tuple(EXPECTED_CENSUS[name].time for name in order)
            and build.tau_parity.time == build.q_parity.time == -1
            and census["V"].time == -1
            and census["J_T"].time == census["u_T"].time
            == census["b_T"].time == 1
        ),
    }


# ---------------------------------------------------------------------------
# Exact dimension evaluator on the live source and field expressions.
# Dimension tuples are ordered (L,T,M).
# ---------------------------------------------------------------------------

Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
DIMENSIONLESS: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))
L_DIM: Dim = (sp.Rational(1), sp.Rational(0), sp.Rational(0))

u_T_functions = tuple(sp.Function(f"u_T{index}") for index in range(1, 4))
u_T_field = sp.Matrix([
    function(*relative_coordinates) for function in u_T_functions
])
b_T_field = sp.Matrix([
    sp.diff(u_T_field[2], xi2) - sp.diff(u_T_field[1], xi3),
    sp.diff(u_T_field[0], xi3) - sp.diff(u_T_field[2], xi1),
    sp.diff(u_T_field[1], xi1) - sp.diff(u_T_field[0], xi2),
])


def add_dims(*dimensions: Dim) -> Dim:
    return tuple(
        sp.Rational(sum(dimension[index] for dimension in dimensions))
        for index in range(3)
    )  # type: ignore[return-value]


def scale_dim(dimension: Dim, factor: sp.Rational) -> Dim:
    return tuple(
        sp.Rational(factor * component) for component in dimension
    )  # type: ignore[return-value]


def expression_dimension(expression: sp.Expr,
                         atom_dimensions: Mapping[Any, Dim]) -> Dim:
    expression = sp.sympify(expression)
    if expression == 0 or expression.is_Number:
        return DIMENSIONLESS
    if expression in atom_dimensions:
        return atom_dimensions[expression]
    if isinstance(expression, sp.Derivative):
        result = expression_dimension(expression.expr, atom_dimensions)
        for variable, count in expression.variable_count:
            result = add_dims(
                result,
                scale_dim(expression_dimension(variable, atom_dimensions), -count),
            )
        return result
    if expression.is_Function:
        if expression.func in atom_dimensions:
            return atom_dimensions[expression.func]
        raise AuditFailure("UNITS_RESTORED", f"unknown function {expression.func}")
    if expression.is_Add:
        term_dimensions = {
            expression_dimension(term, atom_dimensions)
            for term in expression.args
        }
        if len(term_dimensions) != 1:
            raise AuditFailure(
                "UNITS_RESTORED",
                f"inhomogeneous sum {expression}: {term_dimensions}",
            )
        return next(iter(term_dimensions))
    if expression.is_Mul:
        return add_dims(*(
            expression_dimension(factor, atom_dimensions)
            for factor in expression.args
        ))
    if expression.is_Pow:
        base, exponent = expression.args
        if expression_dimension(exponent, atom_dimensions) != DIMENSIONLESS:
            raise AuditFailure("UNITS_RESTORED", "dimensionful exponent")
        return scale_dim(
            expression_dimension(base, atom_dimensions), sp.Rational(exponent)
        )
    raise AuditFailure(
        "UNITS_RESTORED", f"unsupported expression {expression} ({expression.func})"
    )


def dimension_object(*, bad_q: bool = False) -> dict[str, Any]:
    q_dimension: Dim = (
        sp.Rational(1 if bad_q else 0), sp.Rational(-1), sp.Rational(1)
    )
    dimensions: dict[Any, Dim] = {
        **{coordinate: L_DIM for coordinate in relative_coordinates},
        **{
            velocity: (sp.Rational(1), sp.Rational(-1), sp.Rational(0))
            for velocity in body_velocity
        },
        s: DIMENSIONLESS,
        alpha: DIMENSIONLESS,
        q_T: q_dimension,
        lambda_T: q_dimension,
        tau_d: DIMENSIONLESS,
        eta_a: (sp.Rational(-3), sp.Rational(0), sp.Rational(0)),
        **{
            function: (sp.Rational(1), sp.Rational(0), sp.Rational(0))
            for function in u_T_functions
        },
    }
    source = source_build()
    return {
        "sigma": expression_dimension(signed_density, dimensions),
        "continuity_terms": (
            expression_dimension(translation_time_derivative, dimensions),
            expression_dimension(isotropic_flux_divergence, dimensions),
        ),
        "I": tuple(
            expression_dimension(component, dimensions)
            for component in source.native_flux
        ),
        "J_T": tuple(
            expression_dimension(component, dimensions)
            for component in source.native_source
        ),
        "eta_a": expression_dimension(profile, dimensions),
        "V": tuple(
            expression_dimension(component, dimensions)
            for component in body_velocity
        ),
        "u_T": tuple(
            expression_dimension(component, dimensions)
            for component in u_T_field
        ),
        "b_T": tuple(
            expression_dimension(component, dimensions)
            for component in b_T_field
        ),
        "q_T": q_dimension,
    }


# ---------------------------------------------------------------------------
# Source-only dependency and term-inventory guards.
# ---------------------------------------------------------------------------

A_E, q_electric, g_electric = sp.symbols("A_E q_electric g_electric")
r_BA, downstream_landing, downstream_sign = sp.symbols(
    "r_BA downstream_landing downstream_sign"
)
FORBIDDEN_TARGET_SYMBOLS = frozenset({
    A_E, q_electric, g_electric, r_BA, downstream_landing, downstream_sign
})

EXPECTED_TERM_INVENTORY = (
    "signed_density",
    "translation_time_derivative",
    "isotropic_flux_divergence",
    "continuity_residual",
    "unique_alpha",
    "native_flux_I",
    "native_source_J_T",
    "q_T_relation",
    "parity_census_24_cells",
    "P_w_w_reflection_caveat",
    "dimension_firewall",
    "verdict_precedence",
)


def computed_term_inventory(*, include_q_relation: bool = True) -> tuple[str, ...]:
    build = source_build()
    parity = parity_build()
    terms: list[str] = []
    if signed_density.has(s, profile):
        terms.append("signed_density")
    if translation_time_derivative != 0:
        terms.append("translation_time_derivative")
    if isotropic_flux_divergence.has(alpha):
        terms.append("isotropic_flux_divergence")
    if continuity_general.has(alpha):
        terms.append("continuity_residual")
    if unique_alpha_solutions == (sp.Integer(1),):
        terms.append("unique_alpha")
    if len(build.native_flux) == 3:
        terms.append("native_flux_I")
    if len(build.native_source) == 3 and build.native_source.has(q_T):
        terms.append("native_source_J_T")
    if include_q_relation:
        terms.append("q_T_relation")
    if len(parity.census) * len(Parity.__dataclass_fields__) == 24:
        terms.append("parity_census_24_cells")
    terms.append("P_w_w_reflection_caveat")
    if len(dimension_object()) == 9:
        terms.append("dimension_firewall")
    terms.append("verdict_precedence")
    return tuple(terms)


# ---------------------------------------------------------------------------
# Build-native verdict precedence.  Magnitude status is deliberately not a
# source-basis predicate: R1 magnitude still has the conditional identity.
# ---------------------------------------------------------------------------


def derive_verdict(*, continuity_residual: sp.Expr,
                   unique_flux_coefficient: sp.Expr,
                   basis_residuals: Iterable[sp.Expr],
                   barred_intersection: set[sp.Symbol],
                   census_invariants: Mapping[str, bool],
                   magnitude_status: str) -> str:
    del magnitude_status
    source_valid = (
        sp.simplify(continuity_residual) == 0
        and sp.simplify(unique_flux_coefficient - 1) == 0
        and all(sp.simplify(residual) == 0 for residual in basis_residuals)
        and not barred_intersection
    )
    if not source_valid:
        return R1_SOURCE
    if not all(bool(value) for value in census_invariants.values()):
        return PARITY_BAD
    return VERDICT


def verdict_for(source: SourceBuild, parity: ParityBuild,
                *, unique_override: sp.Expr | None = None,
                magnitude_status: str = MAGNITUDE_R1) -> str:
    barred = set(source.candidate_source.free_symbols) & set(
        BARRED_SOURCE_MARKERS
    )
    basis_residuals = tuple(
        sp.expand(source.candidate_source[index] - source.native_source[index])
        for index in range(3)
    )
    return derive_verdict(
        continuity_residual=source.continuity_residual,
        unique_flux_coefficient=(
            source.unique_flux_coefficient
            if unique_override is None else unique_override
        ),
        basis_residuals=basis_residuals,
        barred_intersection=barred,
        census_invariants=parity_invariants(parity),
        magnitude_status=magnitude_status,
    )


# Exact source-build order: all 35 source teeth, no wildcard families.
SOURCE_TOOTH_IDS = (
    "SOURCE_TRANSLATION_CONTINUITY",
    "SOURCE_NOT_IMPORTED",
    "SOURCE_BASIS",
    "PARITY_RW",
    "PARITY_PW",
    "PARITY_ROTATION",
    "PARITY_TIME_REVERSAL",
    "FIELD_IDENTITY_UNITS",
    "ACTION_KINETIC",
    "ACTION_COUPLING",
    "ACTION_STABILITY",
    "G0_DAMAGE",
    "ROUTE_INDEPENDENCE",
    "BOOST_PROJECTOR",
    "BOOST_GENERAL_VELOCITIES",
    "BOOST_NEXT_ORDER",
    "BOOST_COMMON_VELOCITY",
    "DIRECT_SOURCE",
    "DIRECT_PROJECTOR",
    "DIRECT_EXCHANGE_SIGN",
    "DIRECT_FALLOFF",
    "DIRECT_VELOCITY_ORDER",
    "COMPARE_COMPUTED",
    "DELTA_RATIO",
    "CONE_RATIO",
    "QMAG_R1",
    "UNITS_RESTORED",
    "ACTIVE_FLUX_CAVEAT",
    "HOOK_LORENTZ",
    "LEDGER_READY_ROW",
    "TRUTH_TOTALITY",
    "TRUTH_PRECEDENCE",
    "LANDING_OWNERSHIP",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
)


def manifest_entry(identifier: str, disposition: str,
                   owner: str) -> tuple[str, str, str]:
    return identifier, disposition, owner


SOURCE_MANIFEST = (
    manifest_entry("SOURCE_TRANSLATION_CONTINUITY", "PRESERVED",
                   "STAGE035_CONTINUITY_RESIDUAL"),
    manifest_entry("SOURCE_NOT_IMPORTED", "PRESERVED",
                   "STAGE035_FREE_SYMBOL_GUARD"),
    manifest_entry("SOURCE_BASIS", "PRESERVED",
                   "STAGE035_COMPONENT_BASIS"),
    manifest_entry("PARITY_RW", "PRESERVED",
                   "STAGE035_MULTIPLICATIVE_CENSUS"),
    manifest_entry("PARITY_PW", "PRESERVED",
                   "STAGE035_MULTIPLICATIVE_CENSUS"),
    manifest_entry("PARITY_ROTATION", "PRESERVED",
                   "STAGE035_MULTIPLICATIVE_CENSUS"),
    manifest_entry("PARITY_TIME_REVERSAL", "PRESERVED",
                   "STAGE035_MULTIPLICATIVE_CENSUS"),
    manifest_entry("FIELD_IDENTITY_UNITS", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("ACTION_KINETIC", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("ACTION_COUPLING", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("ACTION_STABILITY", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("G0_DAMAGE", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("ROUTE_INDEPENDENCE", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("BOOST_PROJECTOR", "SCOPED_OUT", "STAGE036_V3"),
    manifest_entry("BOOST_GENERAL_VELOCITIES", "SCOPED_OUT", "STAGE036_V3"),
    manifest_entry("BOOST_NEXT_ORDER", "SCOPED_OUT", "STAGE036_V3"),
    manifest_entry("BOOST_COMMON_VELOCITY", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("DIRECT_SOURCE", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("DIRECT_PROJECTOR", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("DIRECT_EXCHANGE_SIGN", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("DIRECT_FALLOFF", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("DIRECT_VELOCITY_ORDER", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("COMPARE_COMPUTED", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("DELTA_RATIO", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("CONE_RATIO", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("QMAG_R1", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("UNITS_RESTORED", "REPLACED_BY_STRONGER",
                   "STAGE035_EXPRESSION_DIMENSION_FIREWALL"),
    manifest_entry("ACTIVE_FLUX_CAVEAT", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("HOOK_LORENTZ", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("LEDGER_READY_ROW", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("TRUTH_TOTALITY", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("TRUTH_PRECEDENCE", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("LANDING_OWNERSHIP", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("TARGET_BLINDNESS", "PRESERVED",
                   "STAGE035_SOURCE_SIDE_DEPENDENCIES"),
    manifest_entry("DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER",
                   "STAGE035_CANONICAL_TERM_INVENTORY"),
)

EXPECTED_MANIFEST_COUNTS = {
    "PRESERVED": 8,
    "REPLACED_BY_STRONGER": 2,
    "SCOPED_OUT": 25,
}
EXPECTED_SCOPED_OUT = tuple(
    identifier for identifier, disposition, _owner in SOURCE_MANIFEST
    if disposition == "SCOPED_OUT"
)


def canonical_manifest_text(
    manifest: Iterable[tuple[str, str, str]]
) -> str:
    return "\n".join("|".join(row) for row in sorted(manifest))


def manifest_sha256(manifest: Iterable[tuple[str, str, str]]) -> str:
    return hashlib.sha256(
        canonical_manifest_text(manifest).encode("utf-8")
    ).hexdigest()


def run_assertions() -> None:
    section("Signed-dent translation continuity and unique isotropic flux")
    selected_alpha = (
        sp.Integer(2)
        if ACTIVE_MUTATION == "SOURCE_TRANSLATION_CONTINUITY"
        else sp.Integer(1)
    )
    live_source = source_build(
        flux_coefficient=selected_alpha,
        inject_barred=ACTIVE_MUTATION == "SOURCE_NOT_IMPORTED",
        source_scale=(
            sp.Integer(2)
            if ACTIVE_MUTATION == "SOURCE_BASIS" else sp.Integer(1)
        ),
    )
    wrong_coefficient_residual = sp.factor(
        continuity_general.subs(alpha, 2)
    )
    continuity_ok = (
        sp.simplify(live_source.continuity_residual) == 0
        and unique_alpha_solutions == (sp.Integer(1),)
        and sp.simplify(continuity_coefficient - (alpha - 1)) == 0
        and sp.simplify(
            continuity_general - (alpha - 1) * translation_gradient
        ) == 0
        and sp.simplify(wrong_coefficient_residual) != 0
    )
    expect_bool(
        "SOURCE_TRANSLATION_CONTINUITY",
        continuity_ok,
        {
            "selected_alpha": selected_alpha,
            "residual": live_source.continuity_residual,
            "general": continuity_general,
            "unique_solutions": unique_alpha_solutions,
            "alpha=2 control": wrong_coefficient_residual,
        },
    )
    print("      sigma_i=s_i*eta_a(x-X_i(t)); residual=(alpha-1) sum_k V_k d_k sigma_i")
    print("      unique isotropic local flux coefficient alpha=1; alpha=2 control is nonzero")

    section("Surviving-solution free-symbol guard")
    source_symbols = set(live_source.candidate_source.free_symbols)
    barred_intersection = source_symbols & set(BARRED_SOURCE_MARKERS)
    expect_bool(
        "SOURCE_NOT_IMPORTED",
        not barred_intersection,
        {
            "barred_intersection": sorted(map(str, barred_intersection)),
            "mapping": {
                display: str(code)
                for display, code in BARRED_DISPLAY_TO_CODE.items()
            },
        },
    )
    print("      BARRED_MAPPING=N_u->Nu,a_T->aT,a'_T->aTp,a_L->aL,q_A^T->q_A_T,q_L->q_L")
    print("      native source REPLACES barred j~sV/pathA_39; no barred amplitude imported")

    section("Continuity-derived tensor source basis")
    basis_residuals = tuple(
        sp.factor(
            live_source.candidate_source[index]
            - live_source.native_source[index]
        )
        for index in range(3)
    )
    q_relation_source = live_source.native_source.subs(
        q_T, lambda_T * tau_d
    )
    expected_q_relation_source = (
        lambda_T * tau_d * signed_density * body_velocity
    )
    expect_bool(
        "SOURCE_BASIS",
        all(sp.simplify(residual) == 0 for residual in basis_residuals)
        and live_source.native_flux
        == signed_density * body_velocity
        and q_relation_source == expected_q_relation_source,
        {
            "component_residuals": basis_residuals,
            "I_i": live_source.native_flux,
            "J_T_i": live_source.candidate_source,
            "q_T_relation": q_T - lambda_T * tau_d,
        },
    )
    print("      I_i=s_i*eta_a(x-X_i)*V_i")
    print("      J_T,i=q_T*I_i=q_T*s_i*eta_a(x-X_i)*V_i; q_T=lambda_T*tau_d")

    section("Derived 24-cell parity census")
    live_parity = parity_build(
        response_rw=-1 if ACTIVE_MUTATION == "PARITY_RW" else 1,
        response_pw=-1 if ACTIVE_MUTATION == "PARITY_PW" else 1,
        tau_time=1 if ACTIVE_MUTATION == "PARITY_TIME_REVERSAL" else -1,
        curl_rotation=(
            "axial_vector"
            if ACTIVE_MUTATION == "PARITY_ROTATION" else "polar_vector"
        ),
    )
    invariants = parity_invariants(live_parity)
    for tooth in (
        "PARITY_RW", "PARITY_PW", "PARITY_ROTATION",
        "PARITY_TIME_REVERSAL",
    ):
        expect_bool(tooth, invariants[tooth], live_parity.census)
    for name in EXPECTED_CENSUS:
        row = live_parity.census[name]
        print(
            f"      {name}: R_w={row.rw:+d}, P_w={row.pw:+d}, "
            f"rotation={row.rotation}, T={row.time:+d}"
        )
    print("      P_w denotes a w-type reflection of the transverse coordinate, NOT ordinary x->-x parity")
    print("      b_T axial and T-even are RAW census facts; departure interpretation is deferred to stage039")

    section("Build-global target blindness and canonical term coverage")
    blindness_source = (
        sum(live_source.native_source)
        + q_T - lambda_T * tau_d + s + sum(body_velocity)
    )
    if ACTIVE_MUTATION == "TARGET_BLINDNESS":
        blindness_source += A_E * signed_density
    forbidden_live = set(blindness_source.free_symbols) & set(
        FORBIDDEN_TARGET_SYMBOLS
    )
    expect_bool(
        "TARGET_BLINDNESS",
        not forbidden_live,
        {"forbidden_live": sorted(map(str, forbidden_live))},
    )
    print("      no electric A_E; no q/g knob; no downstream comparison/landing/sign token")

    inventory = computed_term_inventory(
        include_q_relation=ACTIVE_MUTATION != "DUAL_ENGINE_TERMS"
    )
    expect_bool(
        "DUAL_ENGINE_TERMS",
        inventory == EXPECTED_TERM_INVENTORY,
        {"computed": inventory, "required": EXPECTED_TERM_INVENTORY},
    )
    print(f"      TERM_INVENTORY={','.join(inventory)}")

    section("Whole-stage restored-unit firewall on live expressions")
    restored = dimension_object(
        bad_q=ACTIVE_MUTATION == "UNITS_RESTORED"
    )
    expected_dimensions = {
        "sigma": (-3, 0, 0),
        "continuity_terms": ((-3, -1, 0), (-3, -1, 0)),
        "I": ((-2, -1, 0),) * 3,
        "J_T": ((-2, -2, 1),) * 3,
        "eta_a": (-3, 0, 0),
        "V": ((1, -1, 0),) * 3,
        "u_T": ((1, 0, 0),) * 3,
        "b_T": ((0, 0, 0),) * 3,
        "q_T": (0, -1, 1),
    }
    expect_bool(
        "UNITS_RESTORED",
        restored == expected_dimensions,
        {"computed": restored, "required": expected_dimensions},
    )
    print("      [sigma]=[eta_a]=L^-3; [I]=L^-2 T^-1; [q_T]=M T^-1")
    print("      [J_T]=M L^-2 T^-2; [V]=L T^-1; [u_T]=L; [b_T]=1")

    section("Build-native verdict precedence from computed upstream objects")
    verdict_source = source_build(
        flux_coefficient=(
            sp.Integer(2)
            if ACTIVE_MUTATION == "VERDICT_REDERIVATION"
            else sp.Integer(1)
        )
    )
    clean_parity = parity_build()
    live_verdict = verdict_for(
        verdict_source, clean_parity, magnitude_status=MAGNITUDE_R1
    )
    magnitude_forced_verdict = verdict_for(
        source_build(), clean_parity,
        magnitude_status="magnitude_forced_by_electric",
    )
    continuity_negative = verdict_for(
        source_build(flux_coefficient=sp.Integer(2)), clean_parity
    )
    coefficient_negative = verdict_for(
        source_build(), clean_parity, unique_override=sp.Integer(2)
    )
    basis_negative = verdict_for(
        source_build(source_scale=sp.Integer(2)), clean_parity
    )
    barred_negative = verdict_for(
        source_build(inject_barred=True), clean_parity
    )
    parity_negative = verdict_for(
        source_build(), parity_build(response_rw=-1)
    )
    source_precedence = verdict_for(
        source_build(flux_coefficient=sp.Integer(2)),
        parity_build(response_rw=-1),
    )
    expect_bool(
        "VERDICT_REDERIVATION",
        live_verdict == VERDICT
        and magnitude_forced_verdict == VERDICT
        and continuity_negative == R1_SOURCE
        and coefficient_negative == R1_SOURCE
        and basis_negative == R1_SOURCE
        and barred_negative == R1_SOURCE
        and parity_negative == PARITY_BAD
        and source_precedence == R1_SOURCE,
        {
            "rederived": live_verdict,
            "alpha_negative": continuity_negative,
            "wrong_unique_coefficient": coefficient_negative,
            "basis_negative": basis_negative,
            "barred_negative": barred_negative,
            "parity_negative": parity_negative,
            "source_over_parity": source_precedence,
            "magnitude_R1": MAGNITUDE_R1,
        },
    )
    print(f"      SOURCE_FAILURE_TOKEN={R1_SOURCE}")
    print(f"      CENSUS_FAILURE_TOKEN={PARITY_BAD}")
    print(f"      R1_MAGNITUDE_IDENTITY_TOKEN={live_verdict}")

    section("Canonical source-to-stage predicate manifest")
    manifest = (
        SOURCE_MANIFEST[:-1]
        if ACTIVE_MUTATION == "SOURCE_TO_STAGE_MANIFEST"
        else SOURCE_MANIFEST
    )
    identifiers = tuple(row[0] for row in manifest)
    partition = dict(sorted(Counter(row[1] for row in manifest).items()))
    scoped_out = tuple(
        identifier for identifier, disposition, _owner in manifest
        if disposition == "SCOPED_OUT"
    )
    digest = manifest_sha256(manifest)
    manifest_ok = (
        identifiers == SOURCE_TOOTH_IDS
        and len(identifiers) == len(set(identifiers)) == 35
        and partition == EXPECTED_MANIFEST_COUNTS
        and scoped_out == EXPECTED_SCOPED_OUT
        and len(scoped_out) == 25
        and all(
            owner.startswith((
                "STAGE034_", "STAGE035_", "STAGE036_",
                "STAGE037_", "STAGE038_",
            ))
            for _identifier, _disposition, owner in manifest
        )
        and digest == MANIFEST_DIGEST
    )
    expect_bool(
        "SOURCE_TO_STAGE_MANIFEST",
        manifest_ok,
        {
            "entries": len(manifest),
            "partition": partition,
            "scoped_out": scoped_out,
            "digest": digest,
        },
    )
    print(f"      entries={len(manifest)}; partition={partition}; scoped_out={len(scoped_out)}; digest={digest}")
    print(f"      SCOPED_OUT={','.join(scoped_out)}")

    print("")
    print("SOURCE_LAW=I_i=s_i*eta_a(x-X_i)*V_i; J_T,i=q_T*I_i")
    print("Q_T=lambda_T*tau_d; MAGNITUDE=R1(magnitude)")
    print("PROVENANCE=q_T,tau_d,u_T,eta_a cited from stage034/stage003; not re-derived or re-counted")
    print("SCOPE=source law + raw parity census only; action/routes/landing/departure deferred")
    print(f"VERDICT_TOKEN: {live_verdict}")


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage035_native_source_convection_conditional SymPy audit")
    print("ROUTE=abstract translated-profile derivative identity + multiplicative parity algebra")
    print("FILE_IO=none; CROSS_ENGINE_COMPARE=none")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}")

    run_assertions()
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    print("")
    print(f"TOOTH_COUNT={len(TOOTH_ORDER)}")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print(f"OVERALL PASS: SymPy independently reached {VERDICT}")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(
            "OVERALL FAIL: SymPy stage035 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
