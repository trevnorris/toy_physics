#!/usr/bin/env python3
"""Ledger stage036 SymPy audit: Route-A Maxwell--Darwin reference.

Standalone, print-only, assert-zero, exact, and file-I/O-free.  This engine
derives the k^-4 radial seed coefficient from its Poisson equation, takes its
Cartesian Hessian, reconstructs the transverse projector component by
component, and differentiates the Lorentz-completed coordinate potential.

Tooth-local runtime ablation uses ``LEDGER_STAGE036_MUTATION``.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import os
from typing import Any, Iterable, Mapping

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE036_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

VERDICT = "MAXWELL_DARWIN_REFERENCE"
UNCERTIFIED = "MAXWELL_DARWIN_REFERENCE_UNCERTIFIED"
TIER = "tier_A_conditional"
ELECTRIC_R1 = "R1_REQUIRED(bc_selection)"

TOOTH_ORDER = (
    "BOOST_PROJECTOR",
    "BOOST_GENERAL_VELOCITIES",
    "BOOST_NEXT_ORDER",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
    "UNITS_RESTORED",
    "VERDICT_REDERIVATION",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "BOOST_PROJECTOR":
        "double only the reconstructed kernel checked by BOOST_PROJECTOR",
    "BOOST_GENERAL_VELOCITIES":
        "drop A_V only from the live Lorentz-anchor interaction",
    "BOOST_NEXT_ORDER":
        "raise only the claimed computed velocity order from 2 to 4",
    "TARGET_BLINDNESS":
        "inject a Route-B q_T^2/mu_R term only into the dependency object",
    "DUAL_ENGINE_TERMS":
        "drop radial_force_A2 only from the canonical computed-term inventory",
    "UNITS_RESTORED":
        "make R dimensionless only inside the restored-unit evaluator",
    "VERDICT_REDERIVATION":
        "double a computed projector input consumed only by verdict derivation",
    "SOURCE_TO_STAGE_MANIFEST":
        "drop one scoped-out row only from the canonical 35-tooth manifest",
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


def reduce_unit_sphere(expression: sp.Expr) -> sp.Expr:
    """Use the source build's n_z^2 elimination convention exactly."""
    return sp.factor(
        sp.simplify(expression).subs(nz**2, 1 - nx**2 - ny**2)
    )


# ---------------------------------------------------------------------------
# Exact symbols and velocity invariants.
# ---------------------------------------------------------------------------

R = sp.symbols("R", positive=True)
nx, ny, nz = sp.symbols("n_x n_y n_z", real=True)
nvec = sp.Matrix([nx, ny, nz])
v1x, v1y, v1z, v2x, v2y, v2z = sp.symbols(
    "V1_x V1_y V1_z V2_x V2_y V2_z", real=True
)
velocity_symbols = (v1x, v1y, v1z, v2x, v2y, v2z)
V1 = sp.Matrix([v1x, v1y, v1z])
V2 = sp.Matrix([v2x, v2y, v2z])
s1, s2 = sp.symbols("s_1 s_2", real=True, nonzero=True)
A_E = sp.symbols("A_E", real=True, nonzero=True)
c_gamma = sp.symbols("c_gamma", positive=True)

D_V = sp.expand(V1.dot(V2))
V1n = sp.expand(V1.dot(nvec))
V2n = sp.expand(V2.dot(nvec))
A_V = sp.expand(V1n * V2n)


# ---------------------------------------------------------------------------
# SymPy route: derive the radial seed coefficient, then take its Hessian.
# ---------------------------------------------------------------------------

rx, ry, rz = sp.symbols("r_x r_y r_z", real=True)
rvec = sp.Matrix([rx, ry, rz])
radius = sp.sqrt(rvec.dot(rvec))
seed_coefficient = sp.symbols("seed_coefficient", real=True)
radial_trial = seed_coefficient * radius
radial_laplacian = sp.simplify(sum(
    sp.diff(radial_trial, coordinate, 2) for coordinate in rvec
))
seed_solutions = tuple(sp.solve(
    sp.Eq(radial_laplacian, -1 / (4 * sp.pi * radius)),
    seed_coefficient,
))
SEED_K4 = sp.factor(seed_solutions[0] * radius)
KK_OVER_K4 = sp.Matrix(
    3, 3,
    lambda i, j: -sp.diff(SEED_K4, rvec[i], rvec[j]),
)
TRANSVERSE_CARTESIAN = sp.simplify(
    sp.eye(3) / (4 * sp.pi * radius) - KK_OVER_K4
)
EXPECTED_KERNEL = (
    sp.eye(3) + nvec * nvec.T
) / (8 * sp.pi * R)

KERNEL_AT_R = sp.Matrix(3, 3, lambda i, j: reduce_unit_sphere(
    TRANSVERSE_CARTESIAN[i, j].subs({
        rx: R * nx,
        ry: R * ny,
        rz: R * nz,
    })
))
KERNEL_COMPONENT_RESIDUALS = tuple(
    reduce_unit_sphere(KERNEL_AT_R[i, j] - EXPECTED_KERNEL[i, j])
    for i in range(3) for j in range(3)
)


# ---------------------------------------------------------------------------
# Lorentz-completed electric anchor and coordinate-force differentiation.
# ---------------------------------------------------------------------------

ELECTRIC_U0 = s1 * s2 * A_E / (4 * sp.pi * R)
EXPECTED_A2 = (
    -s1 * s2 * A_E * (D_V + A_V)
    / (8 * sp.pi * c_gamma**2 * R)
)
ROUTE_A_U2 = sp.factor(
    -s1 * s2 * A_E / c_gamma**2
    * (V1.T * EXPECTED_KERNEL * V2)[0]
)
FULL_ANCHOR = sp.factor(ELECTRIC_U0 + ROUTE_A_U2)
EXPECTED_FULL_ANCHOR = sp.factor(
    ELECTRIC_U0 * (1 - (D_V + A_V) / (2 * c_gamma**2))
)

force_coefficient = -s1 * s2 * A_E / (8 * sp.pi * c_gamma**2)
V1r = V1.dot(rvec) / radius
V2r = V2.dot(rvec) / radius
coordinate_u2 = force_coefficient * (D_V + V1r * V2r) / radius
FORCE_CARTESIAN = sp.Matrix([
    -sp.diff(coordinate_u2, coordinate) for coordinate in rvec
])
FORCE_EXPECTED_CARTESIAN = (
    -force_coefficient / radius**2 * (
        V2r * V1 + V1r * V2
        - (D_V + 3 * V1r * V2r) * rvec / radius
    )
)
FORCE_DERIVATION_RESIDUALS = tuple(
    sp.factor(FORCE_CARTESIAN[i] - FORCE_EXPECTED_CARTESIAN[i])
    for i in range(3)
)
FORCE_A2 = sp.Matrix([
    sp.factor(
        s1 * s2 * A_E / (8 * sp.pi * c_gamma**2 * R**2)
        * (
            V2n * V1[i] + V1n * V2[i]
            - (D_V + 3 * A_V) * nvec[i]
        )
    )
    for i in range(3)
])
FORCE_A2_FROM_DERIVATIVE = sp.Matrix([
    reduce_unit_sphere(FORCE_CARTESIAN[i].subs({
        rx: R * nx,
        ry: R * ny,
        rz: R * nz,
    }))
    for i in range(3)
])
FORCE_COMPONENT_RESIDUALS = tuple(
    reduce_unit_sphere(FORCE_A2_FROM_DERIVATIVE[i] - FORCE_A2[i])
    for i in range(3)
)
RADIAL_FORCE_A2 = sp.factor(
    -s1 * s2 * A_E * (D_V + A_V)
    / (8 * sp.pi * c_gamma**2 * R**2)
)
RADIAL_RESIDUAL = reduce_unit_sphere(
    sp.factor(nvec.dot(FORCE_A2) - RADIAL_FORCE_A2)
)


# ---------------------------------------------------------------------------
# Explicit velocity-order runtime object.
# ---------------------------------------------------------------------------

epsilon_v = sp.symbols("epsilon_v", real=True)
velocity_scaled_anchor = sp.expand(
    FULL_ANCHOR.subs({
        symbol: epsilon_v * symbol for symbol in velocity_symbols
    })
)
VELOCITY_POLYNOMIAL = sp.Poly(velocity_scaled_anchor, epsilon_v)
COMPUTED_VELOCITY_ORDERS = tuple(sorted(
    monomial[0] for monomial, coefficient
    in VELOCITY_POLYNOMIAL.terms() if coefficient != 0
))
COMPUTED_VELOCITY_ORDER = max(COMPUTED_VELOCITY_ORDERS)
NEXT_UNCOMPUTED_ORDER = 4


# ---------------------------------------------------------------------------
# Exact dimension evaluator on the live kernel/anchor/force expressions.
# Tuples are ordered (L,T,M).
# ---------------------------------------------------------------------------

Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
DIMENSIONLESS: Dim = (sp.Rational(0),) * 3
L_DIM: Dim = (sp.Rational(1), sp.Rational(0), sp.Rational(0))
V_DIM: Dim = (sp.Rational(1), sp.Rational(-1), sp.Rational(0))
AE_DIM: Dim = (sp.Rational(3), sp.Rational(-2), sp.Rational(1))
E_DIM: Dim = (sp.Rational(2), sp.Rational(-2), sp.Rational(1))
F_DIM: Dim = (sp.Rational(1), sp.Rational(-2), sp.Rational(1))


def add_dims(*dimensions: Dim) -> Dim:
    return tuple(
        sp.Rational(sum(dimension[index] for dimension in dimensions))
        for index in range(3)
    )  # type: ignore[return-value]


def scale_dim(dimension: Dim, factor: sp.Rational) -> Dim:
    return tuple(
        sp.Rational(factor * component) for component in dimension
    )  # type: ignore[return-value]


def expression_dimension(
    expression: sp.Expr,
    atom_dimensions: Mapping[Any, Dim],
) -> Dim:
    expression = sp.sympify(expression)
    if expression == 0 or expression.is_number:
        return DIMENSIONLESS
    if expression in atom_dimensions:
        return atom_dimensions[expression]
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
            expression_dimension(base, atom_dimensions),
            sp.Rational(exponent),
        )
    raise AuditFailure(
        "UNITS_RESTORED",
        f"unsupported expression {expression} ({expression.func})",
    )


def dimension_object(*, bad_radius: bool = False) -> dict[str, Any]:
    radius_dimension = DIMENSIONLESS if bad_radius else L_DIM
    dimensions: dict[Any, Dim] = {
        R: radius_dimension,
        **{component: DIMENSIONLESS for component in nvec},
        **{velocity: V_DIM for velocity in velocity_symbols},
        s1: DIMENSIONLESS,
        s2: DIMENSIONLESS,
        A_E: AE_DIM,
        c_gamma: V_DIM,
        epsilon_v: DIMENSIONLESS,
    }
    return {
        "I_ij": tuple(
            expression_dimension(KERNEL_AT_R[i, j], dimensions)
            for i in range(3) for j in range(3)
        ),
        "A_E": expression_dimension(A_E, dimensions),
        "c_gamma_squared": expression_dimension(c_gamma**2, dimensions),
        "V": tuple(
            expression_dimension(symbol, dimensions)
            for symbol in velocity_symbols
        ),
        "D_V": expression_dimension(D_V, dimensions),
        "A_V": expression_dimension(A_V, dimensions),
        "electric_U0": expression_dimension(ELECTRIC_U0, dimensions),
        "route_A_U2": expression_dimension(ROUTE_A_U2, dimensions),
        "full_anchor": expression_dimension(FULL_ANCHOR, dimensions),
        "force_A2": tuple(
            expression_dimension(component, dimensions)
            for component in FORCE_A2
        ),
        "radial_force_A2": expression_dimension(
            RADIAL_FORCE_A2, dimensions
        ),
        "s": (
            expression_dimension(s1, dimensions),
            expression_dimension(s2, dimensions),
        ),
        "velocity_order": expression_dimension(epsilon_v**2, dimensions),
    }


EXPECTED_DIMENSIONS = {
    "I_ij": ((-1, 0, 0),) * 9,
    "A_E": (3, -2, 1),
    "c_gamma_squared": (2, -2, 0),
    "V": ((1, -1, 0),) * 6,
    "D_V": (2, -2, 0),
    "A_V": (2, -2, 0),
    "electric_U0": (2, -2, 1),
    "route_A_U2": (2, -2, 1),
    "full_anchor": (2, -2, 1),
    "force_A2": ((1, -2, 1),) * 3,
    "radial_force_A2": (1, -2, 1),
    "s": ((0, 0, 0),) * 2,
    "velocity_order": (0, 0, 0),
}


# ---------------------------------------------------------------------------
# Runtime-only target blindness, computed-term coverage, and verdict.
# ---------------------------------------------------------------------------

q_T, mu_R, r_BA, r_cone, delta_u = sp.symbols(
    "q_T mu_R r_BA r_cone Delta_U"
)
sealed_landing = sp.symbols("sealed_landing")
FORBIDDEN_TARGET_SYMBOLS = frozenset({
    q_T, mu_R, r_BA, r_cone, delta_u, sealed_landing,
})

EXPECTED_TERM_INVENTORY = (
    "darwin_kernel_9_components",
    "electric_U0",
    "route_A_U2",
    "full_anchor_UA",
    "force_A2_3_components",
    "radial_force_A2",
    "velocity_orders_0_2_next_4",
    "dimension_firewall",
    "verdict_precedence",
)


def computed_term_inventory(
    *, include_radial_force: bool = True
) -> tuple[str, ...]:
    terms: list[str] = []
    if KERNEL_AT_R.shape == (3, 3):
        terms.append("darwin_kernel_9_components")
    if ELECTRIC_U0.has(A_E, R):
        terms.append("electric_U0")
    if ROUTE_A_U2.has(A_E, c_gamma, R):
        terms.append("route_A_U2")
    if (
        FULL_ANCHOR.has(A_E, R, c_gamma)
        and sp.simplify(FULL_ANCHOR - EXPECTED_FULL_ANCHOR) == 0
    ):
        terms.append("full_anchor_UA")
    if FORCE_A2.shape == (3, 1):
        terms.append("force_A2_3_components")
    if include_radial_force and RADIAL_FORCE_A2.has(R):
        terms.append("radial_force_A2")
    if COMPUTED_VELOCITY_ORDERS == (0, 2):
        terms.append("velocity_orders_0_2_next_4")
    terms.append("dimension_firewall")
    terms.append("verdict_precedence")
    return tuple(terms)


def residuals_are_zero(residuals: Iterable[sp.Expr]) -> bool:
    return all(sp.simplify(residual) == 0 for residual in residuals)


def derive_verdict(
    projector_residuals: Iterable[sp.Expr],
    anchor_residual: sp.Expr,
    velocity_order_is_two: bool,
) -> str:
    certified = (
        residuals_are_zero(projector_residuals)
        and sp.simplify(anchor_residual) == 0
        and bool(velocity_order_is_two)
    )
    return VERDICT if certified else UNCERTIFIED


# Exact source-build order: all 35 teeth, no wildcard families.
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


def manifest_entry(
    identifier: str,
    disposition: str,
    owner: str,
) -> tuple[str, str, str]:
    return identifier, disposition, owner


SOURCE_MANIFEST = (
    manifest_entry("SOURCE_TRANSLATION_CONTINUITY", "SCOPED_OUT",
                   "STAGE035_V2_DONE"),
    manifest_entry("SOURCE_NOT_IMPORTED", "SCOPED_OUT",
                   "STAGE035_V2_DONE"),
    manifest_entry("SOURCE_BASIS", "SCOPED_OUT", "STAGE035_V2_DONE"),
    manifest_entry("PARITY_RW", "SCOPED_OUT", "STAGE035_V2_DONE"),
    manifest_entry("PARITY_PW", "SCOPED_OUT", "STAGE035_V2_DONE"),
    manifest_entry("PARITY_ROTATION", "SCOPED_OUT", "STAGE035_V2_DONE"),
    manifest_entry("PARITY_TIME_REVERSAL", "SCOPED_OUT",
                   "STAGE035_V2_DONE"),
    manifest_entry("FIELD_IDENTITY_UNITS", "SCOPED_OUT",
                   "STAGE034_V1_DONE"),
    manifest_entry("ACTION_KINETIC", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("ACTION_COUPLING", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("ACTION_STABILITY", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("G0_DAMAGE", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("ROUTE_INDEPENDENCE", "SCOPED_OUT", "STAGE037_V4"),
    manifest_entry("BOOST_PROJECTOR", "REPLACED_BY_STRONGER",
                   "STAGE036_INDEPENDENT_PROJECTOR_RECONSTRUCTION"),
    manifest_entry("BOOST_GENERAL_VELOCITIES", "REPLACED_BY_STRONGER",
                   "STAGE036_ANCHOR_FORCE_RADIAL_RECONSTRUCTION"),
    manifest_entry("BOOST_NEXT_ORDER", "PRESERVED",
                   "STAGE036_RUNTIME_VELOCITY_ORDER"),
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
                   "STAGE036_EXPRESSION_DIMENSION_FIREWALL"),
    manifest_entry("ACTIVE_FLUX_CAVEAT", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("HOOK_LORENTZ", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("LEDGER_READY_ROW", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("TRUTH_TOTALITY", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("TRUTH_PRECEDENCE", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("LANDING_OWNERSHIP", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("TARGET_BLINDNESS", "PRESERVED",
                   "STAGE036_REFERENCE_DEPENDENCY_OBJECT"),
    manifest_entry("DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER",
                   "STAGE036_CANONICAL_TERM_INVENTORY"),
)

EXPECTED_MANIFEST_COUNTS = {
    "PRESERVED": 2,
    "REPLACED_BY_STRONGER": 4,
    "SCOPED_OUT": 29,
}
EXPECTED_IN_SCOPE = (
    "BOOST_PROJECTOR",
    "BOOST_GENERAL_VELOCITIES",
    "BOOST_NEXT_ORDER",
    "UNITS_RESTORED",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
)
EXPECTED_SCOPED_OUT = (
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
    "ACTIVE_FLUX_CAVEAT",
    "HOOK_LORENTZ",
    "LEDGER_READY_ROW",
    "TRUTH_TOTALITY",
    "TRUTH_PRECEDENCE",
    "LANDING_OWNERSHIP",
)


def canonical_manifest_text(
    manifest: Iterable[tuple[str, str, str]]
) -> str:
    return "\n".join("|".join(row) for row in sorted(manifest))


def manifest_sha256(
    manifest: Iterable[tuple[str, str, str]]
) -> str:
    return hashlib.sha256(
        canonical_manifest_text(manifest).encode("utf-8")
    ).hexdigest()


MANIFEST_DIGEST = (
    "f8b1569834c1c5cfd404fee4fba7bc49d3d7263f332e551270ad499fd3585cac"
)


def run_assertions() -> None:
    section("Inverse-FT transverse projector from the derived k^-4 seed")
    live_kernel = (
        2 * KERNEL_AT_R
        if ACTIVE_MUTATION == "BOOST_PROJECTOR"
        else KERNEL_AT_R
    )
    live_projector_residuals = tuple(
        reduce_unit_sphere(live_kernel[i, j] - EXPECTED_KERNEL[i, j])
        for i in range(3) for j in range(3)
    )
    seed_ok = (
        seed_solutions == (-1 / (8 * sp.pi),)
        and sp.simplify(SEED_K4 + radius / (8 * sp.pi)) == 0
        and sp.simplify(
            sum(sp.diff(SEED_K4, coordinate, 2) for coordinate in rvec)
            + 1 / (4 * sp.pi * radius)
        ) == 0
    )
    expect_bool(
        "BOOST_PROJECTOR",
        seed_ok and residuals_are_zero(live_projector_residuals),
        {
            "seed": SEED_K4,
            "seed_solutions": seed_solutions,
            "component_residuals": live_projector_residuals,
        },
    )
    print("      F^-1[k^-4]=-R/(8*pi), coefficient derived from Laplacian(seed)=-1/(4*pi*R)")
    print("      I_ij=(delta_ij+n_i*n_j)/(8*pi*R), nine component residuals zero")

    section("Independent-velocity Lorentz anchor and force")
    live_u2 = (
        -s1 * s2 * A_E * D_V
        / (8 * sp.pi * c_gamma**2 * R)
        if ACTIVE_MUTATION == "BOOST_GENERAL_VELOCITIES"
        else ROUTE_A_U2
    )
    live_full_anchor = sp.factor(ELECTRIC_U0 + live_u2)
    anchor_residual = sp.factor(live_u2 - EXPECTED_A2)
    full_anchor_residual = sp.factor(
        live_full_anchor - EXPECTED_FULL_ANCHOR
    )
    expect_bool(
        "BOOST_GENERAL_VELOCITIES",
        sp.simplify(anchor_residual) == 0
        and sp.simplify(full_anchor_residual) == 0
        and residuals_are_zero(FORCE_DERIVATION_RESIDUALS)
        and residuals_are_zero(FORCE_COMPONENT_RESIDUALS)
        and sp.simplify(RADIAL_RESIDUAL) == 0,
        {
            "U_A2_residual": anchor_residual,
            "full_anchor_residual": full_anchor_residual,
            "force_derivation": FORCE_DERIVATION_RESIDUALS,
            "force_components": FORCE_COMPONENT_RESIDUALS,
            "radial_residual": RADIAL_RESIDUAL,
        },
    )
    print("      D_V=V1.V2; A_V=(V1.n)(V2.n)")
    print("      U_A=(s1*s2*A_E/(4*pi*R))*(1-(D_V+A_V)/(2*c_gamma^2))")
    print("      F_A2=(s1*s2*A_E/(8*pi*c_gamma^2*R^2))*[(V2.n)V1+(V1.n)V2-(D_V+3A_V)n]")
    print("      F_A2,r=-s1*s2*A_E*(D_V+A_V)/(8*pi*c_gamma^2*R^2)")

    section("Explicit computed orders and named uncomputed remainder")
    claimed_order = (
        4 if ACTIVE_MUTATION == "BOOST_NEXT_ORDER"
        else COMPUTED_VELOCITY_ORDER
    )
    order_ok = (
        COMPUTED_VELOCITY_ORDERS == (0, 2)
        and claimed_order == 2
        and VELOCITY_POLYNOMIAL.coeff_monomial(epsilon_v**4) == 0
    )
    expect_bool(
        "BOOST_NEXT_ORDER",
        order_ok,
        {
            "computed_orders": COMPUTED_VELOCITY_ORDERS,
            "claimed_order": claimed_order,
            "next_uncomputed": NEXT_UNCOMPUTED_ORDER,
            "epsilon^4_coefficient":
                VELOCITY_POLYNOMIAL.coeff_monomial(epsilon_v**4),
        },
    )
    print("      explicit=O(1)+O(v^2/c_gamma^2); next_uncomputed=O(v^4/c_gamma^4)")

    section("Build-global reference-only dependency and term coverage")
    dependency_expression = (
        sum(KERNEL_AT_R)
        + ELECTRIC_U0 + ROUTE_A_U2 + FULL_ANCHOR
        + sum(FORCE_A2) + RADIAL_FORCE_A2
    )
    if ACTIVE_MUTATION == "TARGET_BLINDNESS":
        dependency_expression += q_T**2 * (D_V + A_V) / (mu_R * R)
    forbidden_live = (
        set(dependency_expression.free_symbols) & set(FORBIDDEN_TARGET_SYMBOLS)
    )
    expect_bool(
        "TARGET_BLINDNESS",
        not forbidden_live,
        {"forbidden_live": sorted(map(str, forbidden_live))},
    )
    print("      reference depends on cited A_E,c_gamma only; no Route-B/comparison/ratio/landing input")

    inventory = computed_term_inventory(
        include_radial_force=ACTIVE_MUTATION != "DUAL_ENGINE_TERMS"
    )
    expect_bool(
        "DUAL_ENGINE_TERMS",
        inventory == EXPECTED_TERM_INVENTORY,
        {"computed": inventory, "required": EXPECTED_TERM_INVENTORY},
    )
    print(f"      TERM_INVENTORY={','.join(inventory)}")

    section("Whole-stage restored-unit firewall on live expressions")
    restored = dimension_object(
        bad_radius=ACTIVE_MUTATION == "UNITS_RESTORED"
    )
    expect_bool(
        "UNITS_RESTORED",
        restored == EXPECTED_DIMENSIONS,
        {"computed": restored, "required": EXPECTED_DIMENSIONS},
    )
    print("      [I_ij]=L^-1; [A_E]=M L^3 T^-2; [c_gamma^2]=L^2 T^-2")
    print("      [D_V]=[A_V]=L^2 T^-2; [U_A]=M L^2 T^-2; [F_A2]=M L T^-2")

    section("Verdict re-derived from projector, anchor, and order objects")
    verdict_kernel = (
        2 * KERNEL_AT_R
        if ACTIVE_MUTATION == "VERDICT_REDERIVATION"
        else KERNEL_AT_R
    )
    verdict_projector_residuals = tuple(
        reduce_unit_sphere(
            verdict_kernel[i, j] - EXPECTED_KERNEL[i, j]
        )
        for i in range(3) for j in range(3)
    )
    live_verdict = derive_verdict(
        verdict_projector_residuals,
        sp.factor(ROUTE_A_U2 - EXPECTED_A2),
        COMPUTED_VELOCITY_ORDER == 2,
    )
    projector_negative = derive_verdict(
        tuple(
            reduce_unit_sphere(
                2 * KERNEL_AT_R[i, j] - EXPECTED_KERNEL[i, j]
            )
            for i in range(3) for j in range(3)
        ),
        sp.Integer(0),
        True,
    )
    anchor_negative = derive_verdict(
        KERNEL_COMPONENT_RESIDUALS,
        sp.factor(
            -s1 * s2 * A_E * D_V
            / (8 * sp.pi * c_gamma**2 * R)
            - EXPECTED_A2
        ),
        True,
    )
    order_negative = derive_verdict(
        KERNEL_COMPONENT_RESIDUALS,
        sp.Integer(0),
        False,
    )
    print(f"      REDERIVED_TOKEN={live_verdict}")
    expect_bool(
        "VERDICT_REDERIVATION",
        live_verdict == VERDICT
        and projector_negative == UNCERTIFIED
        and anchor_negative == UNCERTIFIED
        and order_negative == UNCERTIFIED,
        {
            "live": live_verdict,
            "projector_negative": projector_negative,
            "anchor_negative": anchor_negative,
            "order_negative": order_negative,
        },
    )
    print(f"      FAILURE_TOKEN={UNCERTIFIED}")
    print(f"      OPEN_ELECTRIC_R1={ELECTRIC_R1} is a scope tag, not a certification failure")

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
    in_scope = tuple(
        identifier for identifier, disposition, _owner in manifest
        if disposition != "SCOPED_OUT"
    )
    digest = manifest_sha256(manifest)
    manifest_ok = (
        identifiers == SOURCE_TOOTH_IDS
        and len(identifiers) == len(set(identifiers)) == 35
        and partition == EXPECTED_MANIFEST_COUNTS
        and scoped_out == EXPECTED_SCOPED_OUT
        and len(scoped_out) == 29
        and in_scope == EXPECTED_IN_SCOPE
        and digest == MANIFEST_DIGEST
        and all(
            owner.startswith((
                "STAGE034_", "STAGE035_", "STAGE036_",
                "STAGE037_", "STAGE038_",
            ))
            for _identifier, _disposition, owner in manifest
        )
    )
    expect_bool(
        "SOURCE_TO_STAGE_MANIFEST",
        manifest_ok,
        {
            "entries": len(manifest),
            "partition": partition,
            "in_scope": in_scope,
            "scoped_out": scoped_out,
            "digest": digest,
        },
    )
    print(f"      entries={len(manifest)}; partition={partition}; scoped_out={len(scoped_out)}; digest={digest}")
    print(f"      IN_SCOPE={','.join(in_scope)}")
    print(f"      SCOPED_OUT={','.join(scoped_out)}")

    print("")
    print("REFERENCE_ONLY=Route-A electric boost; Route-B/comparison/ratios deferred to stage037")
    print("PROVENANCE=A_E carries R1_REQUIRED(bc_selection); c_gamma cited from stage003")
    print("SCOPE=EARNED reference kernel at tier_A_conditional; electric sector not re-derived")
    print("REMAINDER=O(v^4/c_gamma^4) named but not computed")
    print(f"VERDICT_TOKEN: {live_verdict}")


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage036_route_a_maxwell_darwin_reference SymPy audit")
    print("ROUTE=Poisson-derived k^-4 seed + Cartesian Hessian + coordinate-force differentiation")
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
            "OVERALL FAIL: SymPy stage036 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
