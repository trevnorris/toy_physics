#!/usr/bin/env python3
"""Ledger stage037 SymPy audit: blind Route B and structural comparison.

Standalone, print-only, assert-zero, exact, and file-I/O-free.  Route B is
constructed before the cited Route-A record and with ``foreign_payload=None``.
This engine solves the direct isotropic tensor constraints, differentiates its
own coordinate potential, performs the prefactor-stripped structural
comparison, and derives the R1-valued coefficient ratios.

Tooth-local runtime ablation uses ``LEDGER_STAGE037_MUTATION``.
"""

from __future__ import annotations

from collections import Counter, OrderedDict
from dataclasses import dataclass
from enum import Enum
import hashlib
import os
from typing import Any, Iterable, Mapping

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE037_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

VERDICT = "BOOST_STRUCTURAL_RELATION_HOLDS"
UNCERTIFIED = "BOOST_STRUCTURAL_RELATION_UNCERTIFIED"
COBLOCKER = "R1_REQUIRED(direct_moving_throat)"
COMPARISON_FACT = "route_B_R1"
RELATIVE_SIGN_FACT = "relative_sign_anchor_conditional"

TOOTH_ORDER = (
    "DIRECT_SOURCE",
    "DIRECT_PROJECTOR",
    "DIRECT_EXCHANGE_SIGN",
    "DIRECT_FALLOFF",
    "DIRECT_VELOCITY_ORDER",
    "ROUTE_INDEPENDENCE",
    "BOOST_COMMON_VELOCITY",
    "COMPARE_COMPUTED",
    "DELTA_RATIO",
    "CONE_RATIO",
    "QMAG_R1",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
    "UNITS_RESTORED",
    "VERDICT_REDERIVATION",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "DIRECT_SOURCE":
        "double the live direct-source coupling in U_B, F_B, and F_B,r",
    "DIRECT_PROJECTOR":
        "add an excess n_i*n_j/(8*pi*mu_R*R) to the live direct kernel",
    "DIRECT_EXCHANGE_SIGN":
        "reverse the live interaction used for the computed parallel/antiparallel signature",
    "DIRECT_FALLOFF":
        "divide the live derived radial force by one extra R",
    "DIRECT_VELOCITY_ORDER":
        "replace the second velocity by a fixed direction in the live interaction",
    "ROUTE_INDEPENDENCE":
        "derive a mutation-only Route-B copy from the already-instantiated Route A",
    "BOOST_COMMON_VELOCITY":
        "add D_V/c_gamma^2 to the Route-A/electric correction ratio",
    "COMPARE_COMPUTED":
        "add a D_V tensor term only to the live prefactor-stripped Route-B comparison",
    "DELTA_RATIO":
        "double the computed r_BA used by the ratio/delta/Delta_U checks",
    "CONE_RATIO":
        "replace c_E^2/c_gamma^2 by c_E/sqrt(c_gamma^2)",
    "QMAG_R1":
        "replace R1(magnitude) by magnitude_forced_by_electric",
    "TARGET_BLINDNESS":
        "inject barred and sealed-decision symbols and remove A_E from the comparison lane",
    "DUAL_ENGINE_TERMS":
        "drop Delta_U and corrupt every remaining canonical symbolic term",
    "UNITS_RESTORED":
        "corrupt every live base dimension, including [q_T]: M*T^-1 -> L*T^-1",
    "VERDICT_REDERIVATION":
        "corrupt verdict-local kernel, comparison, falloff, order, and ancestry objects",
    "SOURCE_TO_STAGE_MANIFEST":
        "remove a scoped-out row and mis-scope DIRECT_SOURCE in the live manifest",
}


class AuditFailure(AssertionError):
    """A named audit predicate failed."""

    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


class DimensionFailure(ValueError):
    """A live expression is dimensionally inhomogeneous or unsupported."""


class MagnitudeFact(Enum):
    R1 = "R1(magnitude)"
    FORCED = "magnitude_forced_by_electric"


@dataclass(frozen=True)
class RouteRecord:
    name: str
    kernel: sp.Matrix
    interaction: sp.Expr
    force: sp.Matrix
    radial_force: sp.Expr
    dependencies: tuple[str, ...]
    ancestry_symbols: frozenset[sp.Symbol]
    foreign_payload_used: bool
    build_ordinal: int


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
    if isinstance(value, (str, type(None), bool, Enum)):
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


def sphere_reduce(expression: sp.Expr) -> sp.Expr:
    return sp.factor(
        sp.simplify(sp.expand(expression).subs(nz**2, 1 - nx**2 - ny**2))
    )


def residuals_are_zero(residuals: Iterable[sp.Expr]) -> bool:
    return all(sp.simplify(residual) == 0 for residual in residuals)


# ---------------------------------------------------------------------------
# Exact symbols, invariants, and cited inputs.
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
q_T = sp.symbols("q_T", real=True, nonzero=True)
mu_R, rho_br = sp.symbols("mu_R rho_br", positive=True)
A_E = sp.symbols("A_E", real=True, nonzero=True)
c_gamma2, c_E = sp.symbols("c_gamma_squared c_E", positive=True)

D_V = sp.expand(V1.dot(V2))
V1n = sp.expand(V1.dot(nvec))
V2n = sp.expand(V2.dot(nvec))
A_V = sp.expand(V1n * V2n)
SHAPE = sp.expand(D_V + A_V)
C_GAMMA_RELATION = {c_gamma2: mu_R / rho_br}

EXPECTED_KERNEL_B = (
    sp.eye(3) + nvec * nvec.T
) / (8 * sp.pi * mu_R * R)
EXPECTED_KERNEL_A = (
    sp.eye(3) + nvec * nvec.T
) / (8 * sp.pi * R)
EXPECTED_U_B = -s1 * s2 * q_T**2 * SHAPE / (
    8 * sp.pi * mu_R * R
)
EXPECTED_FORCE_B = sp.Matrix([
    sp.factor(
        s1 * s2 * q_T**2 / (8 * sp.pi * mu_R * R**2)
        * (
            V2n * V1[index] + V1n * V2[index]
            - (D_V + 3 * A_V) * nvec[index]
        )
    )
    for index in range(3)
])
EXPECTED_RADIAL_B = -s1 * s2 * q_T**2 * SHAPE / (
    8 * sp.pi * mu_R * R**2
)
EXPECTED_U_A2 = -s1 * s2 * A_E * SHAPE / (
    8 * sp.pi * c_gamma2 * R
)
ELECTRIC_RADIAL_FORCE = s1 * s2 * A_E / (4 * sp.pi * R**2)


# ---------------------------------------------------------------------------
# SymPy Route B: solve the two ansatz constraints, then differentiate U_B.
# ---------------------------------------------------------------------------

BUILD_LOG: list[str] = []
Q_CURRENT_SOURCE_TAG, PATHA36_EOM_TAG, DIRECT_ANSATZ_TAG = sp.symbols(
    "Q_CURRENT_source_tag pathA36_transverse_EOM_tag "
    "direct_transverse_tensor_ansatz_tag"
)
ILLICIT_ROUTE_A_TAG = sp.symbols("ILLICIT_ROUTE_A_READ_tag")


def build_stamp(label: str) -> int:
    ordinal = len(BUILD_LOG)
    BUILD_LOG.append(label)
    return ordinal


def force_from_coordinate_potential(
    coupling: sp.Expr,
) -> tuple[sp.Matrix, sp.Expr]:
    """Compute -grad(U) with n=r/|r| and return the unit-vector form."""
    rx, ry, rz = sp.symbols("r_x r_y r_z", real=True)
    rvec = sp.Matrix([rx, ry, rz])
    radius = sp.sqrt(rvec.dot(rvec))
    v1r = V1.dot(rvec) / radius
    v2r = V2.dot(rvec) / radius
    potential = coupling * (D_V + v1r * v2r) / radius
    force_cartesian = sp.Matrix([
        -sp.diff(potential, coordinate) for coordinate in rvec
    ])
    force_at_r = sp.Matrix([
        sphere_reduce(force_cartesian[index].subs({
            rx: R * nx,
            ry: R * ny,
            rz: R * nz,
        }))
        for index in range(3)
    ])
    radial = sphere_reduce(nvec.dot(force_at_r))
    return force_at_r, radial


def derive_route_b(
    foreign_payload: RouteRecord | None,
    *,
    stamp_label: str = "ROUTE_B",
) -> RouteRecord:
    coefficient_a, coefficient_b = sp.symbols(
        "route_B_a route_B_b", real=True
    )
    solution = sp.solve(
        (
            sp.Eq(coefficient_b - coefficient_a, 0),
            sp.Eq(3 * coefficient_a + coefficient_b, 1 / (2 * sp.pi)),
        ),
        (coefficient_a, coefficient_b),
        dict=True,
    )[0]
    kernel = sp.simplify(
        (
            solution[coefficient_a] * sp.eye(3)
            + solution[coefficient_b] * nvec * nvec.T
        ) / (mu_R * R)
    )
    direct_interaction = sp.factor(
        -s1 * s2 * q_T**2 * (V1.T * kernel * V2)[0]
    )
    dependencies = (
        "Q_CURRENT.source",
        "pathA36.transverse_EOM",
        "direct_transverse_tensor_ansatz",
    )
    ancestry = frozenset({
        Q_CURRENT_SOURCE_TAG, PATHA36_EOM_TAG, DIRECT_ANSATZ_TAG,
    })
    interaction = direct_interaction
    name = "B_DIRECT"
    if foreign_payload is not None:
        interaction = sp.factor(
            foreign_payload.interaction * q_T**2 * c_gamma2 / A_E
        )
        dependencies += ("ILLICIT_ROUTE_A_READ",)
        ancestry = ancestry | frozenset({ILLICIT_ROUTE_A_TAG})
        name = "B_DIRECT_FROM_ROUTE_A"
    force, radial_force = force_from_coordinate_potential(
        -s1 * s2 * q_T**2 / (8 * sp.pi * mu_R)
    )
    return RouteRecord(
        name=name,
        kernel=kernel,
        interaction=interaction,
        force=force,
        radial_force=radial_force,
        dependencies=dependencies,
        ancestry_symbols=ancestry,
        foreign_payload_used=foreign_payload is not None,
        build_ordinal=build_stamp(stamp_label),
    )


def cite_route_a() -> RouteRecord:
    """Instantiate the stage036 kernel/anchor as cited data, without re-deriving it."""
    radial_force = -s1 * s2 * A_E * SHAPE / (
        8 * sp.pi * c_gamma2 * R**2
    )
    force = sp.Matrix([
        sp.factor(
            s1 * s2 * A_E / (8 * sp.pi * c_gamma2 * R**2)
            * (
                V2n * V1[index] + V1n * V2[index]
                - (D_V + 3 * A_V) * nvec[index]
            )
        )
        for index in range(3)
    ])
    return RouteRecord(
        name="A_CITED_STAGE036",
        kernel=EXPECTED_KERNEL_A,
        interaction=EXPECTED_U_A2,
        force=force,
        radial_force=radial_force,
        dependencies=("stage036.R70.kernel", "stage036.R70.U_A2"),
        ancestry_symbols=frozenset({sp.Symbol("STAGE036_R70_CITATION")}),
        foreign_payload_used=False,
        build_ordinal=build_stamp("ROUTE_A"),
    )


# Production construction order is load-bearing: Route B exists first and sees
# no Route-A payload.  The cited Route-A record is instantiated only afterward.
ROUTE_B = derive_route_b(None)
ROUTE_A = cite_route_a()


# ---------------------------------------------------------------------------
# Derived direct-route objects, comparison, and coefficient ratios.
# ---------------------------------------------------------------------------

KERNEL_B_RESIDUALS = tuple(
    sphere_reduce(ROUTE_B.kernel[i, j] - EXPECTED_KERNEL_B[i, j])
    for i in range(3) for j in range(3)
)
FORCE_B_RESIDUALS = tuple(
    sphere_reduce(ROUTE_B.force[index] - EXPECTED_FORCE_B[index])
    for index in range(3)
)
RADIAL_B_RESIDUAL = sphere_reduce(
    ROUTE_B.radial_force - EXPECTED_RADIAL_B
)

TENSOR_A = sp.factor(
    ROUTE_A.interaction / (-s1 * s2 * A_E / c_gamma2)
)
TENSOR_B = sp.factor(
    ROUTE_B.interaction / (-s1 * s2 * q_T**2 / mu_R)
)
TENSOR_COMPARISON_RESIDUAL = sp.factor(TENSOR_B - TENSOR_A)

RATIO_BA = sp.factor(ROUTE_B.interaction / ROUTE_A.interaction)
DELTA_BA = sp.factor(RATIO_BA - 1)
CONE_RATIO = sp.factor(c_E**2 / c_gamma2)
DELTA_U = sp.factor(ROUTE_B.interaction - ROUTE_A.interaction)

PARALLEL_RATIO_A = sp.factor(
    (ROUTE_A.radial_force / ELECTRIC_RADIAL_FORCE).subs({
        nx: 1,
        ny: 0,
        nz: 0,
        v1x: 0,
        v2x: 0,
        v1z: 0,
        v2z: 0,
    })
)
PARALLEL_D = v1y * v2y


def radial_power(expression: sp.Expr) -> int:
    exponent = sp.factor(expression).as_powers_dict().get(R, sp.Integer(0))
    if not exponent.is_Integer:
        raise AuditFailure("DIRECT_FALLOFF", f"noninteger R exponent {exponent}")
    return int(-exponent)


epsilon_v = sp.symbols("epsilon_v", real=True)


def velocity_degree_set(expression: sp.Expr) -> tuple[int, ...]:
    scaled = sp.expand(expression.subs({
        symbol: epsilon_v * symbol for symbol in velocity_symbols
    }))
    polynomial = sp.Poly(scaled, epsilon_v)
    return tuple(sorted({
        int(monomial[0])
        for monomial, coefficient in polynomial.terms()
        if coefficient != 0
    }))


COMPUTED_FORCE_POWER = radial_power(ROUTE_B.radial_force)
COMPUTED_POTENTIAL_POWER = radial_power(ROUTE_B.interaction)
COMPUTED_VELOCITY_DEGREES = velocity_degree_set(ROUTE_B.interaction)


def exchange_signature(interaction: sp.Expr) -> tuple[int, int]:
    common = {
        nx: 1, ny: 0, nz: 0,
        v1x: 0, v1y: 1, v1z: 0,
        v2x: 0, v2z: 0,
        s1: 1, s2: 1, q_T: 1, mu_R: 1, R: 1,
    }
    parallel = sp.simplify(interaction.subs({**common, v2y: 1}))
    antiparallel = sp.simplify(interaction.subs({**common, v2y: -1}))
    return int(sp.sign(parallel)), int(sp.sign(antiparallel))


EXCHANGE_SIGNATURE = exchange_signature(ROUTE_B.interaction)


# ---------------------------------------------------------------------------
# Computed dependency/scope objects.
# ---------------------------------------------------------------------------

ROUTE_A_ANCESTRY_SYMBOLS = frozenset({ILLICIT_ROUTE_A_TAG})


def independence_violations(
    route_b: RouteRecord,
    route_a: RouteRecord,
) -> tuple[str, ...]:
    violations: list[str] = []
    if any("ROUTE_A" in dependency for dependency in route_b.dependencies):
        violations.append("ROUTE_A_DEPENDENCY_TAG")
    if route_b.ancestry_symbols & ROUTE_A_ANCESTRY_SYMBOLS:
        violations.append("ROUTE_A_ANCESTRY_SYMBOL")
    if route_b.foreign_payload_used:
        violations.append("FOREIGN_PAYLOAD_USED")
    if route_b.build_ordinal >= route_a.build_ordinal:
        violations.append("ROUTE_B_NOT_FIRST")
    if route_b.name != "B_DIRECT":
        violations.append("NOT_B_DIRECT")
    return tuple(violations)


N_u, a_T, a_T_prime, a_L, q_A_T, q_L = sp.symbols(
    "N_u a_T a_T_prime a_L q_A_T q_L"
)
SEALED_MAGNETISM_LORENTZ, SEALED_AMENDMENT_EXCLUDED = sp.symbols(
    "MAGNETISM_LORENTZ_CONSISTENT AMENDMENT_EXCLUDED"
)
SEALED_MAGNETISM_DEPARTURE, SEALED_NO_GO = sp.symbols(
    "MAGNETISM_DEPARTURE_CHARACTERIZED NO_GO_sector"
)
BARRED_PATHA39_SYMBOLS = frozenset({
    N_u, a_T, a_T_prime, a_L, q_A_T, q_L,
})
SEALED_SECTION4_SYMBOLS = frozenset({
    SEALED_MAGNETISM_LORENTZ,
    SEALED_AMENDMENT_EXCLUDED,
    SEALED_MAGNETISM_DEPARTURE,
    SEALED_NO_GO,
})


def target_blindness_violations(*, mutate: bool) -> tuple[str, ...]:
    direct_expression = (
        sum(ROUTE_B.kernel)
        + ROUTE_B.interaction
        + sum(ROUTE_B.force)
        + ROUTE_B.radial_force
    )
    comparison_expression = (
        ROUTE_A.interaction + ROUTE_B.interaction
        + RATIO_BA + DELTA_BA + CONE_RATIO + DELTA_U
    )
    decision_expression = sp.Integer(0)
    if mutate:
        direct_expression += N_u + SEALED_MAGNETISM_LORENTZ + A_E
        comparison_expression = comparison_expression.subs(A_E, 1)
        decision_expression += SEALED_AMENDMENT_EXCLUDED
    direct_symbols = direct_expression.free_symbols
    comparison_symbols = comparison_expression.free_symbols
    decision_symbols = decision_expression.free_symbols
    violations: list[str] = []
    if direct_symbols & BARRED_PATHA39_SYMBOLS:
        violations.append("BARRED_PATHA39_DIRECT")
    if direct_symbols & SEALED_SECTION4_SYMBOLS:
        violations.append("SEALED_TOKEN_DIRECT")
    if decision_symbols & SEALED_SECTION4_SYMBOLS:
        violations.append("SEALED_DECISION_CHANNEL")
    if A_E in direct_symbols:
        violations.append("A_E_ENTERED_DIRECT_ROUTE")
    if A_E not in comparison_symbols:
        violations.append("A_E_MISSING_FROM_CITED_COMPARISON")
    return tuple(violations)


# ---------------------------------------------------------------------------
# Canonical dual-engine symbolic term contract.
# ---------------------------------------------------------------------------

EXPECTED_TERM_KEYS = (
    "routeB_kernel00",
    "routeB_kernel01",
    "routeB_U2",
    "routeB_Fr",
    "routeA_U2_cited",
    "ratio_BA",
    "delta_BA",
    "cone_ratio",
    "Delta_U",
)


def computed_terms() -> OrderedDict[str, sp.Expr]:
    return OrderedDict((
        ("routeB_kernel00", ROUTE_B.kernel[0, 0]),
        ("routeB_kernel01", ROUTE_B.kernel[0, 1]),
        ("routeB_U2", ROUTE_B.interaction),
        ("routeB_Fr", ROUTE_B.radial_force),
        ("routeA_U2_cited", ROUTE_A.interaction),
        ("ratio_BA", RATIO_BA),
        ("delta_BA", DELTA_BA),
        ("cone_ratio", CONE_RATIO),
        ("Delta_U", DELTA_U),
    ))


EXPECTED_TERMS = OrderedDict((
    ("routeB_kernel00", (1 + nx**2) / (8 * sp.pi * mu_R * R)),
    ("routeB_kernel01", nx * ny / (8 * sp.pi * mu_R * R)),
    ("routeB_U2", EXPECTED_U_B),
    ("routeB_Fr", sphere_reduce(EXPECTED_RADIAL_B)),
    ("routeA_U2_cited", EXPECTED_U_A2),
    ("ratio_BA", q_T**2 * c_gamma2 / (mu_R * A_E)),
    ("delta_BA", q_T**2 * c_gamma2 / (mu_R * A_E) - 1),
    ("cone_ratio", c_E**2 / c_gamma2),
    (
        "Delta_U",
        -s1 * s2 * A_E / (8 * sp.pi * c_gamma2 * R)
        * (q_T**2 * c_gamma2 / (mu_R * A_E) - 1) * SHAPE,
    ),
))


def term_violations(*, mutate: bool) -> tuple[str, ...]:
    live = computed_terms()
    if mutate:
        live = OrderedDict(
            (key, 2 * value) for key, value in live.items()
            if key != "Delta_U"
        )
    violations: list[str] = []
    if tuple(live) != EXPECTED_TERM_KEYS:
        violations.append("TERM_INVENTORY")
    for key in EXPECTED_TERM_KEYS:
        if key not in live:
            violations.append(f"MISSING:{key}")
        elif sp.simplify(live[key] - EXPECTED_TERMS[key]) != 0:
            violations.append(f"SYMBOLIC_MISMATCH:{key}")
    return tuple(violations)


# ---------------------------------------------------------------------------
# Whole-stage exact dimensional firewall on the live expressions.
# Dimensions are ordered (L,T,M).
# ---------------------------------------------------------------------------

Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
DIMENSIONLESS: Dim = (sp.Rational(0),) * 3
L_DIM: Dim = (sp.Rational(1), sp.Rational(0), sp.Rational(0))
V_DIM: Dim = (sp.Rational(1), sp.Rational(-1), sp.Rational(0))
Q_DIM: Dim = (sp.Rational(0), sp.Rational(-1), sp.Rational(1))
MU_DIM: Dim = (sp.Rational(-1), sp.Rational(-2), sp.Rational(1))
RHO_DIM: Dim = (sp.Rational(-3), sp.Rational(0), sp.Rational(1))
AE_DIM: Dim = (sp.Rational(3), sp.Rational(-2), sp.Rational(1))
E_DIM: Dim = (sp.Rational(2), sp.Rational(-2), sp.Rational(1))
F_DIM: Dim = (sp.Rational(1), sp.Rational(-2), sp.Rational(1))
KERNEL_B_DIM: Dim = (
    sp.Rational(0), sp.Rational(2), sp.Rational(-1)
)


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
        dimensions = {
            expression_dimension(term, atom_dimensions)
            for term in expression.args
        }
        if len(dimensions) != 1:
            raise DimensionFailure(
                f"inhomogeneous sum {expression}: {dimensions}"
            )
        return next(iter(dimensions))
    if expression.is_Mul:
        return add_dims(*(
            expression_dimension(factor, atom_dimensions)
            for factor in expression.args
        ))
    if expression.is_Pow:
        base, exponent = expression.args
        if expression_dimension(exponent, atom_dimensions) != DIMENSIONLESS:
            raise DimensionFailure("dimensionful exponent")
        return scale_dim(
            expression_dimension(base, atom_dimensions),
            sp.Rational(exponent),
        )
    raise DimensionFailure(
        f"unsupported expression {expression} ({expression.func})"
    )


UNIT_EXPECTATIONS: OrderedDict[str, Dim] = OrderedDict((
    ("q_T", Q_DIM),
    ("mu_R", MU_DIM),
    ("rho_br", RHO_DIM),
    ("A_E", AE_DIM),
    ("c_gamma_squared", (2, -2, 0)),
    ("c_E", V_DIM),
    ("D_V", (2, -2, 0)),
    ("A_V", (2, -2, 0)),
    ("kernel_B00", KERNEL_B_DIM),
    ("U_B", E_DIM),
    ("F_B0", F_DIM),
    ("F_B1", F_DIM),
    ("F_B2", F_DIM),
    ("F_Br", F_DIM),
    ("U_A2", E_DIM),
    ("r_BA", DIMENSIONLESS),
    ("delta_BA", DIMENSIONLESS),
    ("r_cone", DIMENSIONLESS),
    ("Delta_U", E_DIM),
    ("s_1", DIMENSIONLESS),
    ("s_2", DIMENSIONLESS),
))


def unit_violations(*, bad_units: bool) -> tuple[str, ...]:
    atom_dimensions: dict[Any, Dim] = {
        R: DIMENSIONLESS if bad_units else L_DIM,
        **{
            component: L_DIM if bad_units else DIMENSIONLESS
            for component in nvec
        },
        **{
            velocity: DIMENSIONLESS if bad_units else V_DIM
            for velocity in velocity_symbols
        },
        s1: L_DIM if bad_units else DIMENSIONLESS,
        s2: (0, 1, 0) if bad_units else DIMENSIONLESS,
        q_T: V_DIM if bad_units else Q_DIM,
        mu_R: DIMENSIONLESS if bad_units else MU_DIM,
        rho_br: L_DIM if bad_units else RHO_DIM,
        A_E: L_DIM if bad_units else AE_DIM,
        c_gamma2: DIMENSIONLESS if bad_units else (2, -2, 0),
        c_E: (0, 0, 1) if bad_units else V_DIM,
    }
    expressions: OrderedDict[str, sp.Expr] = OrderedDict((
        ("q_T", q_T),
        ("mu_R", mu_R),
        ("rho_br", rho_br),
        ("A_E", A_E),
        ("c_gamma_squared", c_gamma2),
        ("c_E", c_E),
        ("D_V", D_V),
        ("A_V", A_V),
        ("kernel_B00", ROUTE_B.kernel[0, 0]),
        ("U_B", ROUTE_B.interaction),
        ("F_B0", ROUTE_B.force[0]),
        ("F_B1", ROUTE_B.force[1]),
        ("F_B2", ROUTE_B.force[2]),
        ("F_Br", ROUTE_B.radial_force),
        ("U_A2", ROUTE_A.interaction),
        ("r_BA", RATIO_BA),
        ("delta_BA", DELTA_BA),
        ("r_cone", CONE_RATIO),
        ("Delta_U", DELTA_U),
        ("s_1", s1),
        ("s_2", s2),
    ))
    violations: list[str] = []
    for label, expression in expressions.items():
        try:
            computed = expression_dimension(expression, atom_dimensions)
        except DimensionFailure as exc:
            violations.append(f"{label}:INHOMOGENEOUS:{exc}")
            continue
        if computed != UNIT_EXPECTATIONS[label]:
            violations.append(
                f"{label}:{computed}!={UNIT_EXPECTATIONS[label]}"
            )
    return tuple(violations)


# ---------------------------------------------------------------------------
# Verdict precedence from exactly the computed certification objects.
# ---------------------------------------------------------------------------

def derive_verdict(
    projector_residuals: Iterable[sp.Expr],
    comparison_residual: sp.Expr,
    force_power: int,
    velocity_degrees: tuple[int, ...],
    independence_state: tuple[str, ...],
) -> str:
    certified = (
        residuals_are_zero(projector_residuals)
        and sp.simplify(comparison_residual) == 0
        and force_power == 2
        and velocity_degrees == (2,)
        and independence_state == ()
    )
    return VERDICT if certified else UNCERTIFIED


# ---------------------------------------------------------------------------
# Exact 35-item source-to-stage predicate manifest.
# ---------------------------------------------------------------------------

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
    manifest_entry("ACTION_STABILITY", "SCOPED_OUT",
                   "STAGE034_V1_DONE"),
    manifest_entry("G0_DAMAGE", "SCOPED_OUT", "STAGE034_V1_DONE"),
    manifest_entry("ROUTE_INDEPENDENCE", "REPLACED_BY_STRONGER",
                   "STAGE037_TAG_SYMBOL_ORDER_GUARD"),
    manifest_entry("BOOST_PROJECTOR", "SCOPED_OUT",
                   "STAGE036_V3_CITED"),
    manifest_entry("BOOST_GENERAL_VELOCITIES", "SCOPED_OUT",
                   "STAGE036_V3_CITED"),
    manifest_entry("BOOST_NEXT_ORDER", "SCOPED_OUT",
                   "STAGE036_V3_CITED"),
    manifest_entry("BOOST_COMMON_VELOCITY", "PRESERVED",
                   "STAGE037_ROUTE_A_ELECTRIC_CROSSCHECK"),
    manifest_entry("DIRECT_SOURCE", "REPLACED_BY_STRONGER",
                   "STAGE037_DIRECT_U_FORCE_RECONSTRUCTION"),
    manifest_entry("DIRECT_PROJECTOR", "REPLACED_BY_STRONGER",
                   "STAGE037_CONSTRAINT_SOLVE_COMPONENTS"),
    manifest_entry("DIRECT_EXCHANGE_SIGN", "REPLACED_BY_STRONGER",
                   "STAGE037_COMPUTED_SIGN_SIGNATURE"),
    manifest_entry("DIRECT_FALLOFF", "REPLACED_BY_STRONGER",
                   "STAGE037_EXPRESSION_R_POWER"),
    manifest_entry("DIRECT_VELOCITY_ORDER", "REPLACED_BY_STRONGER",
                   "STAGE037_EXPRESSION_VELOCITY_DEGREE"),
    manifest_entry("COMPARE_COMPUTED", "REPLACED_BY_STRONGER",
                   "STAGE037_TENSOR_ONLY_COMPARISON"),
    manifest_entry("DELTA_RATIO", "REPLACED_BY_STRONGER",
                   "STAGE037_RATIO_DELTA_POTENTIAL_DIFFERENCE"),
    manifest_entry("CONE_RATIO", "REPLACED_BY_STRONGER",
                   "STAGE037_CONE_DENSITY_EQUIVALENCE"),
    manifest_entry("QMAG_R1", "PRESERVED",
                   "STAGE037_R1_MAGNITUDE_ENUM"),
    manifest_entry("UNITS_RESTORED", "REPLACED_BY_STRONGER",
                   "STAGE037_EXPRESSION_DIMENSION_FIREWALL"),
    manifest_entry("ACTIVE_FLUX_CAVEAT", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("HOOK_LORENTZ", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("LEDGER_READY_ROW", "SCOPED_OUT",
                   "STAGE034_V1_DONE"),
    manifest_entry("TRUTH_TOTALITY", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("TRUTH_PRECEDENCE", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("LANDING_OWNERSHIP", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("TARGET_BLINDNESS", "REPLACED_BY_STRONGER",
                   "STAGE037_RUNTIME_SCOPE_CHANNELS"),
    manifest_entry("DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER",
                   "STAGE037_SYMBOLIC_TERM_CONTRACT"),
)

EXPECTED_MANIFEST_COUNTS = {
    "PRESERVED": 2,
    "REPLACED_BY_STRONGER": 12,
    "SCOPED_OUT": 21,
}
EXPECTED_IN_SCOPE = (
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
    "BOOST_PROJECTOR",
    "BOOST_GENERAL_VELOCITIES",
    "BOOST_NEXT_ORDER",
    "ACTIVE_FLUX_CAVEAT",
    "HOOK_LORENTZ",
    "LEDGER_READY_ROW",
    "TRUTH_TOTALITY",
    "TRUTH_PRECEDENCE",
    "LANDING_OWNERSHIP",
)
MANIFEST_DIGEST = (
    "3c88849c5f4f5b7fe05c0d06c59004acccf8c6e85f6823b0e839c41174c8adf4"
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


def manifest_state(
    manifest: tuple[tuple[str, str, str], ...],
) -> tuple[Any, ...]:
    identifiers = tuple(row[0] for row in manifest)
    partition = dict(sorted(Counter(row[1] for row in manifest).items()))
    in_scope = tuple(
        identifier for identifier, disposition, _owner in manifest
        if disposition != "SCOPED_OUT"
    )
    scoped_out = tuple(
        identifier for identifier, disposition, _owner in manifest
        if disposition == "SCOPED_OUT"
    )
    return (
        identifiers,
        partition,
        in_scope,
        scoped_out,
        manifest_sha256(manifest),
    )


EXPECTED_MANIFEST_STATE = (
    SOURCE_TOOTH_IDS,
    EXPECTED_MANIFEST_COUNTS,
    EXPECTED_IN_SCOPE,
    EXPECTED_SCOPED_OUT,
    MANIFEST_DIGEST,
)


# ---------------------------------------------------------------------------
# Assertions and report.
# ---------------------------------------------------------------------------

def run_assertions() -> str:
    section("Blind direct tensor solve and Route-B source/force reconstruction")
    source_scale = (
        sp.Integer(2) if ACTIVE_MUTATION == "DIRECT_SOURCE"
        else sp.Integer(1)
    )
    live_source_u = sp.factor(source_scale * ROUTE_B.interaction)
    live_source_force = sp.Matrix([
        sp.factor(source_scale * component) for component in ROUTE_B.force
    ])
    live_source_radial = sp.factor(
        source_scale * ROUTE_B.radial_force
    )
    source_residuals = (
        sp.factor(live_source_u - EXPECTED_U_B),
        *(
            sphere_reduce(
                live_source_force[index] - EXPECTED_FORCE_B[index]
            )
            for index in range(3)
        ),
        sphere_reduce(live_source_radial - EXPECTED_RADIAL_B),
    )
    expect_bool(
        "DIRECT_SOURCE",
        residuals_are_zero(source_residuals),
        {"residuals": source_residuals},
    )
    print("      U_B=-s1*s2*q_T^2*(D_V+A_V)/(8*pi*mu_R*R)")
    print("      F_B=(s1*s2*q_T^2/(8*pi*mu_R*R^2))*[(V2.n)V1+(V1.n)V2-(D_V+3*A_V)n]")
    print("      F_B,r=-s1*s2*q_T^2*(D_V+A_V)/(8*pi*mu_R*R^2)")

    live_kernel = ROUTE_B.kernel
    if ACTIVE_MUTATION == "DIRECT_PROJECTOR":
        live_kernel = (
            live_kernel
            + nvec * nvec.T / (8 * sp.pi * mu_R * R)
        )
    projector_residuals = tuple(
        sphere_reduce(live_kernel[i, j] - EXPECTED_KERNEL_B[i, j])
        for i in range(3) for j in range(3)
    )
    expect_bool(
        "DIRECT_PROJECTOR",
        residuals_are_zero(projector_residuals),
        {"component_residuals": projector_residuals},
    )
    print("      transversality b=a; trace 3*a+b=1/(2*pi); a=b=1/(8*pi)")
    print("      G_B,ij=(delta_ij+n_i*n_j)/(8*pi*mu_R*R)")

    section("Computed sign, falloff, and velocity-order objects")
    sign_interaction = (
        -ROUTE_B.interaction
        if ACTIVE_MUTATION == "DIRECT_EXCHANGE_SIGN"
        else ROUTE_B.interaction
    )
    live_signature = exchange_signature(sign_interaction)
    expect_bool(
        "DIRECT_EXCHANGE_SIGN",
        live_signature == (-1, 1),
        {"parallel_antiparallel_signature": live_signature},
    )
    print(f"      computed parallel/antiparallel interaction signs={live_signature}")

    falloff_force = (
        ROUTE_B.radial_force / R
        if ACTIVE_MUTATION == "DIRECT_FALLOFF"
        else ROUTE_B.radial_force
    )
    live_force_power = radial_power(falloff_force)
    expect_bool(
        "DIRECT_FALLOFF",
        live_force_power == 2,
        {"computed_force_power": live_force_power},
    )
    print(
        f"      computed U power=R^-{COMPUTED_POTENTIAL_POWER}; "
        f"computed F power=R^-{live_force_power}"
    )

    order_interaction = ROUTE_B.interaction
    if ACTIVE_MUTATION == "DIRECT_VELOCITY_ORDER":
        order_interaction = order_interaction.subs({
            v2x: 1, v2y: 0, v2z: 0,
        })
    live_velocity_degrees = velocity_degree_set(order_interaction)
    expect_bool(
        "DIRECT_VELOCITY_ORDER",
        live_velocity_degrees == (2,),
        {"computed_velocity_degrees": live_velocity_degrees},
    )
    print(f"      computed velocity degrees={live_velocity_degrees}; O(V1*V2)")

    section("Dependency-tag, ancestry-symbol, and execution-order blindness")
    independence_route = ROUTE_B
    if ACTIVE_MUTATION == "ROUTE_INDEPENDENCE":
        independence_route = derive_route_b(
            ROUTE_A,
            stamp_label="ROUTE_B_ILLICIT_COPY",
        )
    independence_state = independence_violations(
        independence_route, ROUTE_A
    )
    expect_bool(
        "ROUTE_INDEPENDENCE",
        independence_state == (),
        {
            "violations": independence_state,
            "dependencies": independence_route.dependencies,
            "build_log": tuple(BUILD_LOG),
        },
    )
    print("      foreign_payload=None; production build order=ROUTE_B then ROUTE_A")
    print(f"      dependencies={ROUTE_B.dependencies}")

    section("Route-A/electric equal-velocity cross-check")
    common_ratio = PARALLEL_RATIO_A
    if ACTIVE_MUTATION == "BOOST_COMMON_VELOCITY":
        common_ratio = sp.factor(common_ratio + PARALLEL_D / c_gamma2)
    common_residual = sp.factor(
        common_ratio + PARALLEL_D / (2 * c_gamma2)
    )
    expect_zero(
        "BOOST_COMMON_VELOCITY",
        common_residual,
        {"parallel_correction_ratio": common_ratio},
    )
    print("      F_A2,r/F_E,r=-(V1.V2)/(2*c_gamma^2), using Route A + electric anchor only")
    print("      equal v: F_r/F_E,r=1-v^2/(2*c_gamma^2)+O(v^4)")

    section("Prefactor-stripped tensor structural comparison")
    live_tensor_b = TENSOR_B
    if ACTIVE_MUTATION == "COMPARE_COMPUTED":
        live_tensor_b = sp.factor(
            live_tensor_b + D_V / (8 * sp.pi * R)
        )
    comparison_residual = sp.factor(live_tensor_b - TENSOR_A)
    expect_zero(
        "COMPARE_COMPUTED",
        comparison_residual,
        {
            "route_A_stripped": TENSOR_A,
            "route_B_stripped": live_tensor_b,
        },
    )
    print("      U_A/(-s1*s2*A_E/c_gamma^2)=U_B/(-s1*s2*q_T^2/mu_R)")
    print("      tensor structure agrees; computed Route-B falloff/order match cited Route A")

    section("Computed coefficient, delta, cone, and potential-difference ratios")
    live_ratio = (
        2 * RATIO_BA
        if ACTIVE_MUTATION == "DELTA_RATIO"
        else RATIO_BA
    )
    live_delta = sp.factor(live_ratio - 1)
    ratio_expected = q_T**2 * c_gamma2 / (mu_R * A_E)
    ratio_residuals = (
        sp.factor(live_ratio - ratio_expected),
        sp.factor(
            live_ratio.subs(C_GAMMA_RELATION)
            - q_T**2 / (rho_br * A_E)
        ),
        sp.factor(live_delta - (ratio_expected - 1)),
        sp.factor(
            DELTA_U
            + s1 * s2 * A_E / (8 * sp.pi * c_gamma2 * R)
            * live_delta * SHAPE
        ),
    )
    expect_bool(
        "DELTA_RATIO",
        residuals_are_zero(ratio_residuals),
        {
            "r_BA": live_ratio,
            "delta_BA": live_delta,
            "Delta_U": DELTA_U,
            "residuals": ratio_residuals,
        },
    )
    print("      r_BA=q_T^2*c_gamma^2/(mu_R*A_E)=q_T^2/(rho_br*A_E)")
    print("      delta_BA=r_BA-1; Delta_U=-(s1*s2*A_E/(8*pi*c_gamma^2*R))*delta_BA*(D_V+A_V)")

    live_cone = (
        c_E / sp.sqrt(c_gamma2)
        if ACTIVE_MUTATION == "CONE_RATIO"
        else CONE_RATIO
    )
    cone_residuals = (
        sp.factor(live_cone - c_E**2 / c_gamma2),
        sp.factor(
            live_cone.subs(C_GAMMA_RELATION)
            - c_E**2 * rho_br / mu_R
        ),
    )
    expect_bool(
        "CONE_RATIO",
        residuals_are_zero(cone_residuals),
        {"r_cone": live_cone, "residuals": cone_residuals},
    )
    print("      r_cone=c_E^2/c_gamma^2=c_E^2*rho_br/mu_R (OPEN)")

    qmag = (
        MagnitudeFact.FORCED
        if ACTIVE_MUTATION == "QMAG_R1"
        else MagnitudeFact.R1
    )
    expect_bool(
        "QMAG_R1",
        qmag is MagnitudeFact.R1,
        {"magnitude_fact": qmag.value},
    )
    print(f"      magnitude={qmag.value}; comparison={COMPARISON_FACT}; relative_sign={RELATIVE_SIGN_FACT}")

    section("Build-global runtime scope, dual-term, and unit firewalls")
    blindness_state = target_blindness_violations(
        mutate=ACTIVE_MUTATION == "TARGET_BLINDNESS"
    )
    expect_bool(
        "TARGET_BLINDNESS",
        blindness_state == (),
        {"violations": blindness_state},
    )
    print("      no SEALED section-4 landing channel; pathA_39 markers absent from Route B")
    print("      A_E enters only the cited comparison/ratio lane")

    term_state = term_violations(
        mutate=ACTIVE_MUTATION == "DUAL_ENGINE_TERMS"
    )
    expect_bool(
        "DUAL_ENGINE_TERMS",
        term_state == (),
        {"violations": term_state},
    )
    print(f"      TERM_INVENTORY={','.join(EXPECTED_TERM_KEYS)}")
    print("      every listed term is checked by exact symbolic residual, never by text")

    units_state = unit_violations(
        bad_units=ACTIVE_MUTATION == "UNITS_RESTORED"
    )
    expect_bool(
        "UNITS_RESTORED",
        units_state == (),
        {"violations": units_state},
    )
    print("      [q_T]=M*T^-1; [mu_R]=M*L^-1*T^-2; [rho_br]=M*L^-3")
    print("      [G_B]=M^-1*T^2; [U_B]=E; [F_B]=E/L; [r_BA]=[delta_BA]=[r_cone]=1")

    section("Verdict re-derived from computed certification objects")
    verdict_kernel = ROUTE_B.kernel
    verdict_tensor_b = TENSOR_B
    verdict_force = ROUTE_B.radial_force
    verdict_interaction = ROUTE_B.interaction
    verdict_independence = independence_violations(ROUTE_B, ROUTE_A)
    if ACTIVE_MUTATION == "VERDICT_REDERIVATION":
        verdict_kernel = (
            verdict_kernel
            + nvec * nvec.T / (8 * sp.pi * mu_R * R)
        )
        verdict_tensor_b = sp.factor(
            verdict_tensor_b + D_V / (8 * sp.pi * R)
        )
        verdict_force = verdict_force / R
        verdict_interaction = verdict_interaction.subs({
            v2x: 1, v2y: 0, v2z: 0,
        })
        verdict_independence = ("ILLICIT_ROUTE_A_READ",)
    verdict_projector_residuals = tuple(
        sphere_reduce(
            verdict_kernel[i, j] - EXPECTED_KERNEL_B[i, j]
        )
        for i in range(3) for j in range(3)
    )
    verdict_compare_residual = sp.factor(verdict_tensor_b - TENSOR_A)
    live_verdict = derive_verdict(
        verdict_projector_residuals,
        verdict_compare_residual,
        radial_power(verdict_force),
        velocity_degree_set(verdict_interaction),
        verdict_independence,
    )
    print(f"      REDERIVED_TOKEN={live_verdict}")
    expect_bool(
        "VERDICT_REDERIVATION",
        live_verdict == VERDICT,
        {
            "live": live_verdict,
            "failure_token": UNCERTIFIED,
            "projector_residuals": verdict_projector_residuals,
            "comparison_residual": verdict_compare_residual,
            "force_power": radial_power(verdict_force),
            "velocity_degrees": velocity_degree_set(verdict_interaction),
            "independence": verdict_independence,
        },
    )
    print(f"      FAILURE_TOKEN={UNCERTIFIED}")
    print("      unresolved q_T/A_E magnitude does not enter structural certification")

    section("Canonical source-to-stage predicate manifest")
    live_manifest = SOURCE_MANIFEST
    if ACTIVE_MUTATION == "SOURCE_TO_STAGE_MANIFEST":
        mutated_rows = list(SOURCE_MANIFEST[1:])
        mutated_rows = [
            (
                row[0],
                "SCOPED_OUT" if row[0] == "DIRECT_SOURCE" else row[1],
                row[2],
            )
            for row in mutated_rows
        ]
        live_manifest = tuple(mutated_rows)
    live_manifest_state = manifest_state(live_manifest)
    expect_bool(
        "SOURCE_TO_STAGE_MANIFEST",
        live_manifest_state == EXPECTED_MANIFEST_STATE,
        {
            "partition": live_manifest_state[1],
            "in_scope": live_manifest_state[2],
            "scoped_out": live_manifest_state[3],
            "digest": live_manifest_state[4],
        },
    )
    print(
        f"      entries=35; partition={EXPECTED_MANIFEST_COUNTS}; "
        f"scoped_out=21; digest={MANIFEST_DIGEST}"
    )
    print(f"      IN_SCOPE={','.join(EXPECTED_IN_SCOPE)}")
    print(f"      SCOPED_OUT={','.join(EXPECTED_SCOPED_OUT)}")

    print("")
    print("BUILD_ORDER=Route B first with foreign_payload=None; Route A cited afterward")
    print("STRUCTURAL_AGREEMENT=tensor + R^-2 falloff + O(V1*V2) velocity order")
    print("RATIO_STATUS=DECIDED expressions; R1-valued through q_T and A_E")
    print(f"COMPARISON_FACT={COMPARISON_FACT}")
    print(f"RELATIVE_SIGN_FACT={RELATIVE_SIGN_FACT}")
    print(f"COBLOCKER={COBLOCKER}")
    print("STAGE038_BOUNDARY=no SEALED section-4 landing is adjudicated here")
    print(f"VERDICT_TOKEN: {live_verdict}")
    return live_verdict


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage037_route_b_boost_structural_relation SymPy audit")
    print("ROUTE=linear ansatz solve + coordinate-potential force differentiation + direct symbolic quotients")
    print("BLIND_BUILD=Route B instantiated before cited Route A")
    print("FILE_IO=none; CROSS_ENGINE_COMPARE=none")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}")

    live_verdict = run_assertions()
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    print("")
    print(f"TOOTH_COUNT={len(TOOTH_ORDER)}")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print(
        "OVERALL PASS: SymPy independently reached "
        f"{live_verdict}"
    )


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(
            "OVERALL FAIL: SymPy stage037 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
