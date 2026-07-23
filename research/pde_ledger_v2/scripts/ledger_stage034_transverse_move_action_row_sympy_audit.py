#!/usr/bin/env python3
"""Ledger stage034 SymPy audit: the moving-throat transverse action row.

Standalone, print-only, assert-zero, exact, and file-I/O-free.  The imported
stage003/pathA_36 transverse row is instantiated but not re-earned.  Stage034
earns only the finite-profile moving coupling.  Tooth-local runtime ablation
uses ``LEDGER_STAGE034_MUTATION``.

This engine's physics route is an explicit two-polarization kinetic/gradient
Hessian construction followed by generalized-eigenvalue dispersion.
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
MUTATION_ENV = "LEDGER_STAGE034_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

VERDICT = "TRANSVERSE_MOVE_ACTION_ROW"
NO_INCONSISTENCY = "none"
MANIFEST_DIGEST = "4343bd60cd974f653a0a8ac2eeced6c7aca15b1831c81be20bd358f449c454af"

TOOTH_ORDER = (
    "ACTION_KINETIC",
    "ACTION_COUPLING",
    "ACTION_STABILITY",
    "G0_DAMAGE",
    "LEDGER_READY_ROW",
    "FIELD_IDENTITY_UNITS",
    "GUARD_IMPORT_VS_EARN",
    "TARGET_BLINDNESS",
    "DUAL_ENGINE_TERMS",
    "UNITS_RESTORED",
    "VERDICT_REDERIVATION",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "ACTION_KINETIC": "rho_br/2 -> rho_br/3 in the differentiated action density",
    "ACTION_COUPLING": "q_T=lambda_T*tau_d -> lambda_T in the differentiated source",
    "ACTION_STABILITY": "rho_br -> -rho_br before recomputing transverse Hessian positivity",
    "G0_DAMAGE": "absorb the moving row into the parsed active F_flux ledger row",
    "LEDGER_READY_ROW": "make the local moving coupling quadratic in u_T",
    "FIELD_IDENTITY_UNITS": "[b_T]=1 -> L^-1 in the field-identity dimension object",
    "GUARD_IMPORT_VS_EARN": "misclassify imported u_T_kinetic as earned by stage034",
    "TARGET_BLINDNESS": "inject electric A_E into the live action dependency graph",
    "DUAL_ENGINE_TERMS": "drop q_T_relation from the computed term inventory",
    "UNITS_RESTORED": "[q_T]=M T^-1 -> M T^-2 in the whole-density firewall",
    "VERDICT_REDERIVATION": "make the verdict pipeline's kinetic Hessian non-PD",
    "SOURCE_TO_STAGE_MANIFEST": "drop one deferred source tooth from the canonical partition",
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


rho_br, mu_R = sp.symbols("rho_br mu_R", positive=True)
q_T, lambda_T, tau_d = sp.symbols("q_T lambda_T tau_d", nonzero=True)
s, eta_a = sp.symbols("s eta_a", nonzero=True)
k = sp.symbols("k", positive=True)
omega_sq = sp.symbols("omega_sq", real=True)
c_gamma_sq = mu_R / rho_br

u = sp.Matrix(sp.symbols("u1:4", real=True))
u_dot = sp.Matrix(sp.symbols("ud1:4", real=True))
curl_u = sp.Matrix(sp.symbols("b1:4", real=True))
velocity = sp.Matrix(sp.symbols("V1:4", real=True))

A_E, q_electric, g_electric = sp.symbols("A_E q_electric g_electric")


@dataclass(frozen=True)
class ActionBuild:
    density: sp.Expr
    kinetic: sp.Expr
    gradient: sp.Expr
    coupling: sp.Expr
    q_definition: sp.Expr


def action_build(*, kinetic_denominator: int = 2, rho_sign: int = 1,
                 mu_sign: int = 1, nonlinear_coupling: bool = False,
                 q_definition: sp.Expr | None = None,
                 electric_leak: bool = False) -> ActionBuild:
    """Construct this engine's own local action density."""
    q_def = lambda_T * tau_d if q_definition is None else q_definition
    kinetic = sp.Rational(1, kinetic_denominator) * rho_sign * rho_br * u_dot.dot(u_dot)
    gradient = -sp.Rational(1, 2) * mu_sign * mu_R * curl_u.dot(curl_u)
    linear_form = velocity.dot(u)
    coupling = q_T * s * eta_a * linear_form
    if nonlinear_coupling:
        coupling *= u[0]
    density = kinetic + gradient + coupling
    if electric_leak:
        density += A_E * u[0]
    return ActionBuild(
        density=sp.expand(density),
        kinetic=kinetic,
        gradient=gradient,
        coupling=coupling,
        q_definition=q_def,
    )


def transverse_hessian_data(build: ActionBuild) -> dict[str, Any]:
    """Differentiate the action and reduce to k along z, u_T in span(e_x,e_y)."""
    kinetic_hessian = sp.hessian(build.kinetic, tuple(u_dot))
    curl_stiffness = -sp.hessian(build.gradient, tuple(curl_u))
    polarization_basis = sp.Matrix([[1, 0], [0, 1], [0, 0]])
    wave_vector = sp.Matrix([0, 0, k])
    transverse_checks = tuple(
        sp.simplify(wave_vector.dot(polarization_basis[:, index]))
        for index in range(polarization_basis.cols)
    )
    kinetic_T = sp.simplify(polarization_basis.T * kinetic_hessian * polarization_basis)
    gradient_T = sp.simplify(
        k**2 * polarization_basis.T * curl_stiffness * polarization_basis
    )
    wave_operator = gradient_T - omega_sq * kinetic_T
    roots = tuple(
        sp.solve(sp.Eq(wave_operator[index, index], 0), omega_sq)[0]
        for index in range(2)
    )
    kinetic_eigenvalues = tuple(
        eigenvalue
        for eigenvalue, multiplicity in kinetic_T.eigenvals().items()
        for _ in range(multiplicity)
    )
    gradient_eigenvalues = tuple(
        eigenvalue
        for eigenvalue, multiplicity in gradient_T.eigenvals().items()
        for _ in range(multiplicity)
    )
    return {
        "kinetic_hessian": kinetic_hessian,
        "curl_stiffness": curl_stiffness,
        "kinetic_T": kinetic_T,
        "gradient_T": gradient_T,
        "wave_operator": wave_operator,
        "dispersion_roots": roots,
        "kinetic_eigenvalues": kinetic_eigenvalues,
        "gradient_eigenvalues": gradient_eigenvalues,
        "transverse_checks": transverse_checks,
        "polarization_count": polarization_basis.cols,
    }


def stability_object(rho_sign: int = 1, mu_sign: int = 1) -> dict[str, Any]:
    data = transverse_hessian_data(
        action_build(rho_sign=rho_sign, mu_sign=mu_sign)
    )
    kinetic_pd = all(
        True if value.is_positive is True else False
        for value in data["kinetic_eigenvalues"]
    )
    gradient_pd = all(
        True if value.is_positive is True else False
        for value in data["gradient_eigenvalues"]
    )
    dispersion_nonnegative = all(
        True if value.is_positive is True else False
        for value in data["dispersion_roots"]
    )
    return {
        **data,
        "kinetic_pd": kinetic_pd,
        "gradient_pd": gradient_pd,
        "dispersion_nonnegative": dispersion_nonnegative,
        "stable": kinetic_pd and gradient_pd and dispersion_nonnegative,
    }


G0_ROW_TRANSCRIPT = """
bulk_scalar|rho,theta retained scalar action
electric_scalar|localized H/h scalar row
longitudinal|u_L retained scalar-longitudinal row
drain|active drain throughput Gamma_0
return_F_flux|remote return momentum ledger ACTIVE_DEFERRED
wall_r_B|holonomic r_B wall and passive reaction
geon|held-out geon rest constant
zero_bulk_scalar_couplings|r_BH,r_B^2H^2,Hrho,Hdelta_rho,Hdt_theta,Hgrad_theta=0
zero_source_modulation|delta_J_m and neighbor response=0
zero_scalar_transverse_mixing|r_Bu_T,Hu_T,u_Lu_T,two-gradient mixing=0
zero_cross_kinetic|dt_u_L_dt_h and Berry rows=0
zero_scalar_masses|u_L^2,h^2,higher gradients and nonlinearities=0
zero_bulk_modulus|independent B(divu)^2=0
zero_brane_phase|theta_B and brane-phase drain=0
zero_wall_dynamics|bending,anchoring,surface storage,dissipation=0
zero_drain_response|dynamic drain and return responses=0
zero_direct_drain_sources|direct h,u_L drain sources=0
zero_geon_derivatives|field-dependent geon derivatives=0
zero_viscosity|viscosity,drag,no-slip,permeability,phase-jump=0
zero_other_branches|E4,E5,E1 and mixture terms=0
prohibited_ancestry|Maxwell/gauge fields,point sources,native current law,Coulomb prior absent
""".strip()


def parse_g0_rows(transcript: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    for raw_line in transcript.splitlines():
        key, payload = raw_line.split("|", 1)
        if key in rows:
            raise AuditFailure("G0_DAMAGE", f"duplicate parsed G0 row: {key}")
        rows[key] = payload
    return rows


def g0_damage_object(*, absorb_flux: bool = False) -> dict[str, Any]:
    before = parse_g0_rows(G0_ROW_TRANSCRIPT)
    amendment = {
        "delta_transverse_action": "imported pathA_36 u_T row",
        "delta_moving_coupling": "q_T sum_i s_i eta_a V_i.u_T",
    }
    if absorb_flux:
        amendment["return_F_flux"] = "SILENTLY_ABSORBED_INTO_CONSERVATIVE_ROW"
    after = dict(before)
    after.update(amendment)
    preexisting = set(before)
    changed = tuple(sorted(
        key for key in preexisting if after.get(key) != before[key]
    ))
    overlap = tuple(sorted(preexisting.intersection(amendment)))
    flux_untouched = after["return_F_flux"] == before["return_F_flux"]
    internal = NO_INCONSISTENCY if not changed and not overlap and flux_untouched else "g0-damage"
    return {
        "parsed_row_count": len(before),
        "changed_preexisting": changed,
        "amendment_overlap": overlap,
        "flux_untouched": flux_untouched,
        "new_rows": tuple(sorted(set(after) - preexisting)),
        "internal_inconsistency": internal,
    }


def ledger_ready_object(*, nonlinear: bool = False) -> dict[str, Any]:
    build = action_build(nonlinear_coupling=nonlinear)
    source = sp.Matrix([sp.diff(build.coupling, component) for component in u])
    polynomial_degree = sp.Poly(build.coupling, *tuple(u)).total_degree()
    local = not build.density.has(sp.Integral, sp.Derivative)
    source_linear = all(not component.has(*tuple(u)) for component in source)
    variational_responses = {
        "momentum": sp.Matrix([sp.diff(build.density, component) for component in u_dot]),
        "curl_response": sp.Matrix([sp.diff(build.density, component) for component in curl_u]),
        "field_source": source,
    }
    variational = (
        variational_responses["momentum"] == rho_br * u_dot
        and variational_responses["curl_response"] == -mu_R * curl_u
        and any(component != 0 for component in source)
    )
    return {
        "local": local,
        "coupling_degree": polynomial_degree,
        "source_linear": source_linear,
        "variational": variational,
        "responses": variational_responses,
        "well_formed": local and polynomial_degree == 1 and source_linear and variational,
    }


Dimension = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dimension = (sp.Rational(0), sp.Rational(0), sp.Rational(0))
L_DIM: Dimension = (sp.Rational(1), sp.Rational(0), sp.Rational(0))
T_DIM: Dimension = (sp.Rational(0), sp.Rational(1), sp.Rational(0))
M_DIM: Dimension = (sp.Rational(0), sp.Rational(0), sp.Rational(1))
ENERGY_DIM: Dimension = (sp.Rational(2), sp.Rational(-2), sp.Rational(1))
DENSITY_DIM: Dimension = (sp.Rational(-1), sp.Rational(-2), sp.Rational(1))
ACTION_DIM: Dimension = (sp.Rational(2), sp.Rational(-1), sp.Rational(1))


def dim_add(*dimensions: Dimension) -> Dimension:
    return tuple(sum(items, sp.Rational(0)) for items in zip(*dimensions))  # type: ignore[return-value]


def dim_scale(dimension: Dimension, power: sp.Rational) -> Dimension:
    return tuple(power * item for item in dimension)  # type: ignore[return-value]


class InhomogeneousDimension(ValueError):
    pass


def dimension_of(expression: sp.Expr, dimensions: Mapping[sp.Symbol, Dimension]) -> Dimension:
    expression = sp.sympify(expression)
    if expression.is_Number:
        return ZERO_DIM
    if expression.is_Symbol:
        if expression not in dimensions:
            raise InhomogeneousDimension(f"unregistered symbol {expression}")
        return dimensions[expression]
    if expression.is_Add:
        term_dimensions = tuple(dimension_of(term, dimensions) for term in expression.args)
        if len(set(term_dimensions)) != 1:
            raise InhomogeneousDimension(str(term_dimensions))
        return term_dimensions[0]
    if expression.is_Mul:
        return dim_add(*(dimension_of(factor, dimensions) for factor in expression.args))
    if expression.is_Pow and expression.exp.is_Rational:
        return dim_scale(dimension_of(expression.base, dimensions), expression.exp)
    raise InhomogeneousDimension(f"unsupported expression {expression}")


def primitive_dimensions(*, bad_b: bool = False, bad_q: bool = False) -> dict[sp.Symbol, Dimension]:
    dims: dict[sp.Symbol, Dimension] = {
        rho_br: (-3, 0, 1),
        mu_R: (-1, -2, 1),
        q_T: (0, -2 if bad_q else -1, 1),
        lambda_T: (0, -1, 1),
        tau_d: ZERO_DIM,
        s: ZERO_DIM,
        eta_a: (-3, 0, 0),
    }
    dims.update({component: L_DIM for component in u})
    dims.update({component: (1, -1, 0) for component in u_dot})
    dims.update({component: ((-1, 0, 0) if bad_b else ZERO_DIM) for component in curl_u})
    dims.update({component: (1, -1, 0) for component in velocity})
    return dims


def dimension_object(*, bad_b: bool = False, bad_q: bool = False) -> dict[str, Any]:
    dims = primitive_dimensions(bad_b=bad_b, bad_q=bad_q)
    build = action_build()
    term_dimensions: dict[str, Any] = {}
    for name, expression in (
        ("kinetic", build.kinetic),
        ("gradient", build.gradient),
        ("coupling", build.coupling),
        ("density", build.density),
    ):
        try:
            term_dimensions[name] = dimension_of(expression, dims)
        except InhomogeneousDimension as exc:
            term_dimensions[name] = ("INHOMOGENEOUS", str(exc))
    density_dimension = term_dimensions["density"]
    action_dimension = (
        dim_add(density_dimension, (3, 1, 0))
        if density_dimension == DENSITY_DIM else ("INHOMOGENEOUS",)
    )
    field_identity = {
        "u_T": dims[u[0]],
        "u_dot_T": dims[u_dot[0]],
        "curl_u_T": dims[curl_u[0]],
        "rho_br": dims[rho_br],
        "mu_R": dims[mu_R],
        "q_T": dims[q_T],
        "eta_a": dims[eta_a],
        "b_T": dims[curl_u[0]],
    }
    return {
        "field_identity": field_identity,
        "term_dimensions": term_dimensions,
        "action_dimension": action_dimension,
    }


def computed_term_inventory(*, include_q_relation: bool = True) -> tuple[str, ...]:
    build = action_build()
    groups: set[str] = set()
    for term in sp.Add.make_args(sp.expand(build.density)):
        if any(term.has(component) for component in u_dot):
            groups.add("u_T_kinetic")
        elif any(term.has(component) for component in curl_u):
            groups.add("u_T_gradient")
        elif any(term.has(component) for component in u):
            groups.add("moving_coupling")
        else:
            groups.add("unclassified")
    groups.update({"c_gamma_relation", "transverse_constraint"})
    if include_q_relation:
        groups.add("q_T_relation")
    return tuple(sorted(groups))


EXPECTED_TERM_INVENTORY = (
    "c_gamma_relation",
    "moving_coupling",
    "q_T_relation",
    "transverse_constraint",
    "u_T_gradient",
    "u_T_kinetic",
)


def accounting_object(*, overclaim: bool = False) -> dict[str, Any]:
    rows = {
        "u_T_kinetic": "IMPORTED_STAGE003_PATHA36",
        "u_T_gradient": "IMPORTED_STAGE003_PATHA36",
        "c_gamma_relation": "IMPORTED_STAGE003_PATHA36",
        "moving_coupling": "EARNED_STAGE034",
    }
    if overclaim:
        rows["u_T_kinetic"] = "EARNED_STAGE034"
    earned = tuple(sorted(key for key, tag in rows.items() if tag == "EARNED_STAGE034"))
    imported = tuple(sorted(
        key for key, tag in rows.items() if tag == "IMPORTED_STAGE003_PATHA36"
    ))
    return {"earned_new_terms": earned, "imported_terms": imported, "tags": rows}


def dependency_object(*, electric_leak: bool = False) -> dict[str, Any]:
    build = action_build(electric_leak=electric_leak)
    forbidden = {A_E, q_electric, g_electric}
    live_symbols = set(build.density.free_symbols)
    forbidden_live = tuple(sorted((live_symbols & forbidden), key=str))
    allowed_source_symbols = {
        rho_br, mu_R, q_T, s, eta_a, *tuple(u), *tuple(u_dot), *tuple(curl_u), *tuple(velocity)
    }
    unknown = tuple(sorted(live_symbols - allowed_source_symbols, key=str))
    return {
        "forbidden_live": forbidden_live,
        "unknown_live": unknown,
        "source_side_only": not forbidden_live and not unknown,
    }


def derive_verdict(stability: Mapping[str, Any], g0: Mapping[str, Any],
                   row: Mapping[str, Any]) -> tuple[str, str]:
    sectors: list[str] = []
    if not bool(stability["stable"]):
        sectors.append("stability-fail")
    if g0["internal_inconsistency"] != NO_INCONSISTENCY:
        sectors.append("g0-damage")
    if not bool(row["well_formed"]):
        sectors.append("row-malformed")
    internal = NO_INCONSISTENCY if not sectors else "+".join(sectors)
    if "stability-fail" in sectors:
        token = "ROW_UNSTABLE"
    elif "g0-damage" in sectors:
        token = "G0_DAMAGED"
    elif sectors:
        token = "AMENDMENT_INCONSISTENT"
    else:
        token = VERDICT
    return token, internal


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


def manifest_entry(identifier: str, disposition: str, owner: str) -> tuple[str, str, str]:
    return identifier, disposition, owner


SOURCE_MANIFEST = (
    manifest_entry("SOURCE_TRANSLATION_CONTINUITY", "SCOPED_OUT", "STAGE035_V2"),
    manifest_entry("SOURCE_NOT_IMPORTED", "SCOPED_OUT", "STAGE035_V2"),
    manifest_entry("SOURCE_BASIS", "SCOPED_OUT", "STAGE035_V2"),
    manifest_entry("PARITY_RW", "SCOPED_OUT", "STAGE035_V2"),
    manifest_entry("PARITY_PW", "SCOPED_OUT", "STAGE035_V2"),
    manifest_entry("PARITY_ROTATION", "SCOPED_OUT", "STAGE035_V2"),
    manifest_entry("PARITY_TIME_REVERSAL", "SCOPED_OUT", "STAGE035_V2"),
    manifest_entry("FIELD_IDENTITY_UNITS", "REPLACED_BY_STRONGER", "STAGE034_DIMENSION_OBJECT"),
    manifest_entry("ACTION_KINETIC", "REPLACED_BY_STRONGER", "STAGE034_HESSIAN_DISPERSION"),
    manifest_entry("ACTION_COUPLING", "REPLACED_BY_STRONGER", "STAGE034_DIFFERENTIATED_SOURCE"),
    manifest_entry("ACTION_STABILITY", "REPLACED_BY_STRONGER", "STAGE034_PD_HESSIANS"),
    manifest_entry("G0_DAMAGE", "REPLACED_BY_STRONGER", "STAGE034_PARSED_G0_DIFF"),
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
    manifest_entry("UNITS_RESTORED", "REPLACED_BY_STRONGER", "STAGE034_WHOLE_DENSITY_FIREWALL"),
    manifest_entry("ACTIVE_FLUX_CAVEAT", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("HOOK_LORENTZ", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("LEDGER_READY_ROW", "REPLACED_BY_STRONGER", "STAGE034_LOCAL_VARIATIONAL_ROW"),
    manifest_entry("TRUTH_TOTALITY", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("TRUTH_PRECEDENCE", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("LANDING_OWNERSHIP", "SCOPED_OUT", "STAGE038_V5"),
    manifest_entry("TARGET_BLINDNESS", "PRESERVED", "STAGE034_SOURCE_SIDE_DEPENDENCIES"),
    manifest_entry("DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER", "STAGE034_CANONICAL_TERM_INVENTORY"),
)

EXPECTED_MANIFEST_COUNTS = {
    "PRESERVED": 1,
    "REPLACED_BY_STRONGER": 8,
    "SCOPED_OUT": 26,
}


def canonical_manifest_text(manifest: Iterable[tuple[str, str, str]]) -> str:
    return "\n".join("|".join(row) for row in sorted(manifest))


def manifest_sha256(manifest: Iterable[tuple[str, str, str]]) -> str:
    return hashlib.sha256(canonical_manifest_text(manifest).encode("utf-8")).hexdigest()


def run_assertions() -> None:
    section("Imported pathA_36 transverse row: coefficients, dispersion, and two polarizations")
    kinetic_denominator = 3 if ACTIVE_MUTATION == "ACTION_KINETIC" else 2
    kinetic_build = action_build(kinetic_denominator=kinetic_denominator)
    kinetic_data = transverse_hessian_data(kinetic_build)
    kinetic_coefficient = sp.diff(kinetic_build.kinetic, u_dot[0], 2) / 2
    gradient_coefficient = sp.diff(kinetic_build.gradient, curl_u[0], 2) / 2
    kinetic_ok = (
        sp.simplify(kinetic_coefficient - rho_br / 2) == 0
        and sp.simplify(gradient_coefficient + mu_R / 2) == 0
        and kinetic_data["kinetic_T"] == rho_br * sp.eye(2)
        and kinetic_data["gradient_T"] == mu_R * k**2 * sp.eye(2)
        and all(sp.simplify(root - c_gamma_sq * k**2) == 0
                for root in kinetic_data["dispersion_roots"])
        and kinetic_data["polarization_count"] == 2
        and kinetic_data["transverse_checks"] == (0, 0)
    )
    expect_bool(
        "ACTION_KINETIC",
        kinetic_ok,
        {
            "kinetic_coefficient": kinetic_coefficient,
            "gradient_coefficient": gradient_coefficient,
            "dispersion": kinetic_data["dispersion_roots"],
            "c_gamma^2": c_gamma_sq,
        },
    )
    print("      imported provenance=stage003/pathA_36; omega^2=(mu_R/rho_br) k^2; polarizations=2")

    section("Magnetism-new finite-profile moving coupling")
    q_definition = lambda_T if ACTIVE_MUTATION == "ACTION_COUPLING" else lambda_T * tau_d
    coupling_build = action_build(q_definition=q_definition)
    source_vector = sp.Matrix([
        sp.diff(coupling_build.coupling.subs(q_T, coupling_build.q_definition), component)
        for component in u
    ])
    expected_source = lambda_T * tau_d * s * eta_a * velocity
    radial_coordinate, mouth_radius = sp.symbols("r mouth_radius", positive=True)
    normalized_profile = sp.exp(-radial_coordinate**2 / mouth_radius**2) / (
        sp.pi ** sp.Rational(3, 2) * mouth_radius**3
    )
    profile_integral = sp.integrate(
        4 * sp.pi * radial_coordinate**2 * normalized_profile,
        (radial_coordinate, 0, sp.oo),
    )
    coupling_ok = (
        source_vector == expected_source
        and sp.simplify(coupling_build.q_definition - lambda_T * tau_d) == 0
        and sp.simplify(profile_integral - 1) == 0
        and coupling_build.coupling.has(q_T, s, eta_a, *tuple(velocity))
    )
    expect_bool(
        "ACTION_COUPLING",
        coupling_ok,
        {
            "dL/du_T": source_vector,
            "q_definition": coupling_build.q_definition,
            "profile_integral": profile_integral,
        },
    )
    print("      earned_new=moving_coupling; q_T=lambda_T*tau_d; integral eta_a d^3x=1")

    section("Computed transverse Hessian positivity")
    production_stability = stability_object(
        rho_sign=-1 if ACTIVE_MUTATION == "ACTION_STABILITY" else 1
    )
    ghost_control = stability_object(rho_sign=-1)
    tachyon_control = stability_object(mu_sign=-1)
    stability_ok = (
        production_stability["stable"]
        and not ghost_control["stable"]
        and not tachyon_control["stable"]
        and production_stability["kinetic_eigenvalues"] == (rho_br, rho_br)
        and production_stability["gradient_eigenvalues"] == (k**2 * mu_R, k**2 * mu_R)
    )
    expect_bool(
        "ACTION_STABILITY",
        stability_ok,
        {
            "kinetic_eigenvalues": production_stability["kinetic_eigenvalues"],
            "gradient_eigenvalues": production_stability["gradient_eigenvalues"],
            "ghost_control": ghost_control["stable"],
            "tachyon_control": tachyon_control["stable"],
        },
    )
    print("      kinetic_Hessian=PD; gradient_Hessian=PD; no_ghost; no_tachyon")

    section("Parsed G0 row diff and active-flux preservation")
    damage = g0_damage_object(absorb_flux=ACTIVE_MUTATION == "G0_DAMAGE")
    expect_bool(
        "G0_DAMAGE",
        damage["parsed_row_count"] == 21
        and damage["changed_preexisting"] == ()
        and damage["amendment_overlap"] == ()
        and damage["flux_untouched"]
        and damage["new_rows"] == ("delta_moving_coupling", "delta_transverse_action")
        and damage["internal_inconsistency"] == NO_INCONSISTENCY,
        damage,
    )
    print("      scalar/drain/return_F_flux/wall_r_B/geon/declared-zero rows unchanged")
    print("      active F_flux=untouched,deferred(V-5/Part-VII)")

    section("Local variational G0+delta row")
    ledger = ledger_ready_object(nonlinear=ACTIVE_MUTATION == "LEDGER_READY_ROW")
    expect_bool(
        "LEDGER_READY_ROW",
        ledger["well_formed"]
        and ledger["coupling_degree"] == 1
        and ledger["source_linear"]
        and ledger["local"]
        and ledger["variational"],
        ledger,
    )
    print("      one local density; linear current.field source; imported+new pieces variational")

    section("Field identity and restored units")
    field_dimensions = dimension_object(bad_b=ACTIVE_MUTATION == "FIELD_IDENTITY_UNITS")
    expected_field_identity = {
        "u_T": (1, 0, 0),
        "u_dot_T": (1, -1, 0),
        "curl_u_T": (0, 0, 0),
        "rho_br": (-3, 0, 1),
        "mu_R": (-1, -2, 1),
        "q_T": (0, -1, 1),
        "eta_a": (-3, 0, 0),
        "b_T": (0, 0, 0),
    }
    expect_bool(
        "FIELD_IDENTITY_UNITS",
        field_dimensions["field_identity"] == expected_field_identity
        and all(value == DENSITY_DIM
                for value in field_dimensions["term_dimensions"].values())
        and field_dimensions["action_dimension"] == ACTION_DIM,
        field_dimensions,
    )
    print("      [u_T]=L; [u_dot_T]=L/T; [curl u_T]=[b_T]=1; [S]=E*T")

    section("Provenance accounting and build-global guards")
    accounting = accounting_object(overclaim=ACTIVE_MUTATION == "GUARD_IMPORT_VS_EARN")
    expect_bool(
        "GUARD_IMPORT_VS_EARN",
        accounting["earned_new_terms"] == ("moving_coupling",)
        and accounting["imported_terms"] == (
            "c_gamma_relation", "u_T_gradient", "u_T_kinetic"
        ),
        accounting,
    )
    print("      earned_new_terms={moving_coupling}; imported=stage003/pathA_36")

    dependencies = dependency_object(electric_leak=ACTIVE_MUTATION == "TARGET_BLINDNESS")
    expect_bool(
        "TARGET_BLINDNESS",
        dependencies["source_side_only"]
        and dependencies["forbidden_live"] == ()
        and dependencies["unknown_live"] == (),
        dependencies,
    )
    print("      no A_E; no electric q/g knob; no downstream sign/landing token")

    inventory = computed_term_inventory(
        include_q_relation=ACTIVE_MUTATION != "DUAL_ENGINE_TERMS"
    )
    expect_bool(
        "DUAL_ENGINE_TERMS",
        inventory == EXPECTED_TERM_INVENTORY and "unclassified" not in inventory,
        {"computed": inventory, "required": EXPECTED_TERM_INVENTORY},
    )
    print(f"      TERM_INVENTORY={','.join(inventory)}")

    restored = dimension_object(bad_q=ACTIVE_MUTATION == "UNITS_RESTORED")
    expect_bool(
        "UNITS_RESTORED",
        restored["term_dimensions"] == {
            "kinetic": DENSITY_DIM,
            "gradient": DENSITY_DIM,
            "coupling": DENSITY_DIM,
            "density": DENSITY_DIM,
        }
        and restored["action_dimension"] == ACTION_DIM
        and restored["field_identity"]["q_T"] == (0, -1, 1),
        restored,
    )
    print("      whole action density=[E L^-3]; dt d^3x density=[E T]")

    section("Computed verdict re-derivation")
    verdict_stability = stability_object(
        rho_sign=-1 if ACTIVE_MUTATION == "VERDICT_REDERIVATION" else 1
    )
    live_verdict, live_internal = derive_verdict(
        verdict_stability, g0_damage_object(), ledger_ready_object()
    )
    named_alternative, alternative_internal = derive_verdict(
        stability_object(rho_sign=-1), g0_damage_object(), ledger_ready_object()
    )
    expect_bool(
        "VERDICT_REDERIVATION",
        live_verdict == VERDICT
        and live_internal == NO_INCONSISTENCY
        and named_alternative == "ROW_UNSTABLE"
        and alternative_internal == "stability-fail",
        {
            "rederived": live_verdict,
            "internal_inconsistency": live_internal,
            "named_negative_control": named_alternative,
            "negative_internal": alternative_internal,
        },
    )

    section("Canonical source-to-stage predicate manifest")
    manifest = (
        SOURCE_MANIFEST[:-1]
        if ACTIVE_MUTATION == "SOURCE_TO_STAGE_MANIFEST"
        else SOURCE_MANIFEST
    )
    identifiers = tuple(row[0] for row in manifest)
    counts = dict(sorted(Counter(row[1] for row in manifest).items()))
    digest = manifest_sha256(manifest)
    deferred = tuple(row[0] for row in manifest if row[1] == "SCOPED_OUT")
    expected_deferred = tuple(
        identifier for identifier, disposition, _owner in SOURCE_MANIFEST
        if disposition == "SCOPED_OUT"
    )
    manifest_ok = (
        identifiers == SOURCE_TOOTH_IDS
        and len(identifiers) == len(set(identifiers)) == 35
        and counts == EXPECTED_MANIFEST_COUNTS
        and deferred == expected_deferred
        and len(deferred) == 26
        and all(
            owner.startswith(("STAGE034_", "STAGE035_", "STAGE036_", "STAGE037_", "STAGE038_"))
            for _identifier, _disposition, owner in manifest
        )
        and digest == MANIFEST_DIGEST
    )
    expect_bool(
        "SOURCE_TO_STAGE_MANIFEST",
        manifest_ok,
        {"entries": len(manifest), "partition": counts, "deferred": deferred, "digest": digest},
    )
    print(f"      entries={len(manifest)}; partition={counts}; deferred={len(deferred)}; digest={digest}")

    print("")
    print("ACTION_ROW=S_T+move:int dt d^3x [rho_br/2 |u_dot_T|^2 - mu_R/2 |curl u_T|^2 + q_T sum_i s_i eta_a V_i.u_T]")
    print("TRANSVERSE_CONSTRAINT=div u_T=0; POLARIZATIONS=2")
    print("C_GAMMA_SQUARED=mu_R/rho_br")
    print("Q_T=lambda_T*tau_d")
    print("PROVENANCE=stage003/pathA_36 imported-not-earned")
    print(f"internal_inconsistency={live_internal}")
    print(f"VERDICT_TOKEN: {live_verdict}")


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage034_transverse_move_action_row SymPy audit")
    print("ROUTE=explicit transverse kinetic/gradient Hessians + generalized-eigenvalue dispersion")
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
    print("OVERALL PASS: SymPy independently reached TRANSVERSE_MOVE_ACTION_ROW")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage034 audit did not close ({exc.predicate})")
        raise SystemExit(1)
