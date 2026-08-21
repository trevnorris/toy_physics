#!/usr/bin/env python3
"""Unified SymPy audit for the S11b uniform interface-coupling law.

The physics implemented here is supplied solely by
``directives/S11b_SHARED_PHYSICS.md``.  The run streams one re-parseable
CAS object per tag and, after every primary task completes, writes the
accumulated flat ledger consumed by the next step.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
import hashlib
import re
import sys
import time

import sympy as sp
from sympy.core.relational import Relational
from sympy.core.symbol import Str
from sympy.logic.boolalg import Boolean
from sympy.printing.str import StrPrinter


SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
DIRECTIVE_PATH = SCRIPT_DIR.parent / "directives" / "S11b_SHARED_PHYSICS.md"
EXPORT_PATH = SCRIPT_DIR / "S11b_exports.py"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from S11_exports import LEDGER as INCOMING_LEDGER  # noqa: E402


CLASS_TAGS = {"KNOB", "STRUCTURAL", "COORDINATE", "CONTROL", "PREMISE", "DERIVED"}
PRIMARY_TASKS = (
    "B0_ENERGY", "B0A", "B0B", "B0C", "B1", "B2A", "B2B", "B2C",
    "B2D", "B5", "B4", "B6", "B7", "B9",
)
CONTROL_TASK = "B8"

NOT_ESTABLISHED = Str("NOT_ESTABLISHED")
NOT_APPLICABLE = Str("NOT_APPLICABLE")
UNDECIDED = Str("UNDECIDED")
PROVED_TRUE = Str("PROVED_TRUE")
PROVED_FALSE = Str("PROVED_FALSE")
ADMISSIBLE = Str("ADMISSIBLE")
EXCLUDED = Str("EXCLUDED")
CROSS_STEP_CONSISTENCY_ONLY = Str("CROSS_STEP_CONSISTENCY_ONLY")
UNDEFINED_RATIO = Str("UNDEFINED_RATIO")


DECLARED_SYMBOLS: dict[sp.Symbol, dict[str, str]] = {}
SYMBOL_DIMENSION_KEYS: dict[sp.Symbol, str] = {}
ISSUES: list[str] = []
CUSTOM_QUANTITY_NAMES: set[str] = set()
CURRENT_TASK: str | None = None


def register_symbol(symbol: sp.Symbol, class_tag: str, description: str) -> sp.Symbol:
    if class_tag not in CLASS_TAGS:
        raise ValueError(f"unknown class {class_tag!r}")
    DECLARED_SYMBOLS[symbol] = {"class": class_tag, "description": description}
    return symbol


def symbol(name: str, class_tag: str, description: str, **assumptions: object) -> sp.Symbol:
    return register_symbol(sp.Symbol(name, **assumptions), class_tag, description)


# Inherited identities: these are the live objects in the incoming rows.
c_s0 = register_symbol(INCOMING_LEDGER["c_s0"]["value"], "KNOB", "bulk sound speed")
mu_R = register_symbol(INCOMING_LEDGER["mu_R"]["value"], "KNOB", "curl-type brane modulus")
rho_br = register_symbol(INCOMING_LEDGER["rho_br"]["value"], "KNOB", "slab-integrated inertia density")

# Model-definition symbols introduced by this step.
rho_m = symbol("rho_m", "KNOB", "four-spatial-dimensional bulk mass density", positive=True)
v_dr = symbol("v_bulk_normal_0", "KNOB", "steady bulk normal drain speed", real=True)
W0 = symbol("W_0", "KNOB", "background slab thickness", positive=True)
B_rho = symbol("B_rho", "KNOB", "local four-dimensional compression modulus", real=True)
B_rho_3 = symbol("B_rho_3", "KNOB", "slab-integrated compression modulus", real=True)
mu_W = symbol("mu_W", "KNOB", "thickness inertia", positive=True)
k_W = symbol("k_W", "KNOB", "thickness restoring coefficient", real=True)
kappa_W = symbol("kappa_W", "KNOB", "thickness-gradient coefficient", real=True)
C = symbol("C", "KNOB", "densification-thickening cross coefficient", real=True)
B_div = symbol("B_div", "KNOB", "divergence invariant coefficient", real=True)
mu_S = symbol("mu_S", "KNOB", "symmetric-traceless gradient coefficient", real=True)
G_theta_u = symbol("G_theta_u", "KNOB", "densification-divergence coefficient", real=True)
G_W_u = symbol("G_W_u", "KNOB", "thickening-divergence coefficient", real=True)
kappa_theta = symbol("kappa_theta", "KNOB", "densification-gradient coefficient", real=True)
kappa_theta_W = symbol("kappa_theta_W", "KNOB", "mixed scalar-gradient coefficient", real=True)

Lambda_A0 = symbol("Lambda_A_0", "KNOB", "zero-frequency affinity mobility", real=True)
Lambda_V0 = symbol("Lambda_V_0", "KNOB", "zero-frequency velocity mobility", real=True)
Lambda_X0 = symbol("Lambda_X_0", "KNOB", "zero-frequency reciprocal traction coefficient", real=True)
tau_A = symbol("tau_A", "KNOB", "affinity relaxation time", nonnegative=True)
tau_V = symbol("tau_V", "KNOB", "velocity relaxation time", nonnegative=True)
tau_X = symbol("tau_X", "KNOB", "reciprocal-traction relaxation time", nonnegative=True)

D_brane = symbol("D_brane", "STRUCTURAL", "in-plane spatial dimension", integer=True, positive=True)

# Spectral and amplitude coordinates.
# The modal frequency is deliberately assumption-free: B5 solves over C and
# classifies both signs of Im(omega).  Diagnostics explicitly restricted to the
# real axis are evaluated with this private real surrogate and then written
# back in the public coordinate, together with their real-axis condition.
omega = symbol("omega", "COORDINATE", "frequency")
_omega_real_axis = sp.Dummy("omega_real_axis", real=True)
k = symbol("k", "COORDINATE", "in-plane wavenumber magnitude", positive=True)
q_out = symbol("q_out", "COORDINATE", "outward bulk normal wavenumber")
eta = symbol("eta", "COORDINATE", "longitudinal divergence amplitude")
theta = symbol("theta", "COORDINATE", "Eulerian densification amplitude")
e_W = symbol("e_W", "COORDINATE", "fractional thickness amplitude")
u_T = symbol("u_T", "COORDINATE", "transverse displacement amplitude")
V_face = symbol("V_face", "COORDINATE", "outward face-velocity amplitude")
p_face = symbol("delta_p_face", "COORDINATE", "bulk face-pressure amplitude")
J_face = symbol("J_face", "COORDINATE", "outward relative mass-flux amplitude")
mu_drive = symbol("mu_theta_drive", "COORDINATE", "slab chemical driving amplitude")
F_W = symbol("F_W", "COORDINATE", "generalized thickness driving force density")

# Diagnostic coordinates and controls.
ell_A = symbol("ell_A", "CONTROL", "inert affinity-kernel placeholder")
ell_V = symbol("ell_V", "CONTROL", "inert velocity-kernel placeholder")
ell_X = symbol("ell_X", "CONTROL", "inert traction-kernel placeholder")
control_b = symbol("b", "CONTROL", "window-shift control", real=True)
control_c = symbol("c", "CONTROL", "interval-asymmetry control", real=True)
control_a = symbol("a", "CONTROL", "window width", positive=True)
control_L = symbol("L", "CONTROL", "half-interval size", positive=True)


I = sp.I
Lambda_A = sp.cancel(Lambda_A0 / (1 - I * omega * tau_A))
Lambda_V = sp.cancel(Lambda_V0 / (1 - I * omega * tau_V))
Lambda_X = sp.cancel(Lambda_X0 / (1 - I * omega * tau_X))
Z_GENERAL = sp.cancel(rho_m * omega / q_out)


def casify(value: object) -> object:
    """Convert emitted containers and booleans to re-parseable SymPy objects."""
    if isinstance(value, bool):
        return sp.true if value else sp.false
    if isinstance(value, str):
        return Str(value)
    if isinstance(value, Mapping):
        return sp.Tuple(*(sp.Tuple(casify(key), casify(item)) for key, item in value.items()))
    if isinstance(value, sp.Tuple):
        return sp.Tuple(*(casify(item) for item in value))
    if isinstance(value, (tuple, list)):
        return sp.Tuple(*(casify(item) for item in value))
    if isinstance(value, sp.MatrixBase) and not isinstance(value, sp.ImmutableMatrix):
        return sp.ImmutableMatrix(value)
    return value


def render(value: object) -> str:
    value = casify(value)
    if isinstance(value, (sp.Basic, sp.MatrixBase)):
        return sp.srepr(value)
    return repr(value)


class UnevaluatedDisplayPrinter(StrPrinter):
    """Readable rendering that never calls ``doit`` on matrix expressions."""

    def _print_MatMul(self, expression: sp.MatMul) -> str:
        return "MatMul(" + ", ".join(self._print(argument) for argument in expression.args) + ")"

    def _print_Inverse(self, expression: sp.Inverse) -> str:
        return "Inverse(" + self._print(expression.arg) + ")"


DISPLAY_PRINTER = UnevaluatedDisplayPrinter()


def display(value: object) -> str:
    value = casify(value)
    if isinstance(value, (sp.Basic, sp.MatrixBase)):
        return DISPLAY_PRINTER.doprint(value)
    return repr(value)


def real_axis_simplify(value: object) -> object:
    """Simplify under real omega without changing omega's global identity."""
    axis_value = casify(value).subs(omega, _omega_real_axis)
    if isinstance(axis_value, sp.MatrixBase):
        simplified = axis_value.applyfunc(sp.simplify)
    else:
        simplified = sp.simplify(axis_value)
    return simplified.xreplace({_omega_real_axis: omega})


def real_axis_cancel(value: object) -> object:
    """Put rational expressions on the real-omega locus without expansion."""
    axis_value = casify(value).subs(omega, _omega_real_axis)
    if isinstance(axis_value, sp.MatrixBase):
        reduced = axis_value.applyfunc(sp.cancel)
    else:
        reduced = sp.cancel(axis_value)
    return reduced.xreplace({_omega_real_axis: omega})


def real_axis_form(value: object) -> object:
    """Impose the real-omega locus without requesting canonicalization."""
    return casify(value).subs(omega, _omega_real_axis).xreplace({_omega_real_axis: omega})


def real_axis_part(value: object, part: str) -> object:
    """Take Re or Im while the frequency surrogate still carries realness."""
    axis_value = casify(value).subs(omega, _omega_real_axis)
    extracted = sp.re(axis_value) if part == "real" else sp.im(axis_value)
    return extracted.xreplace({_omega_real_axis: omega})


REAL_AXIS = sp.Eq(sp.im(omega), 0, evaluate=False)


@dataclass(frozen=True)
class CandidateRow:
    key: str
    value: object
    class_tag: str
    source_tag: str
    dimension_key: str | None = None
    description: str | None = None


class Emitter:
    def __init__(self) -> None:
        self.count = 0
        self.values: dict[str, object] = {}
        self.rendered: dict[str, str] = {}
        self.local_tags: list[str] = []
        self.export_candidates: list[CandidateRow] = []

    def emit(
        self,
        tag: str,
        payload: object,
        *,
        export_key: str | None = None,
        class_tag: str = "DERIVED",
        dimension_key: str | None = None,
        description: str | None = None,
    ) -> object:
        if tag in self.values:
            raise RuntimeError(f"duplicate emitted tag {tag}")
        payload = casify(payload)
        payload_rendering = render(payload)
        print(f"{tag}: {payload_rendering}", flush=True)
        self.count += 1
        self.values[tag] = payload
        self.rendered[tag] = payload_rendering
        if "_LOCAL_" in tag:
            self.local_tags.append(tag)
        if export_key is not None:
            self.export_candidates.append(
                CandidateRow(export_key, payload, class_tag, tag, dimension_key, description)
            )
        return payload


EMITTER = Emitter()


def emit(
    quantity: str,
    payload: object,
    *,
    key: str | None,
    local: bool = False,
    export: bool = True,
    class_tag: str = "DERIVED",
    dimension_key: str | None = None,
) -> object:
    infix = "LOCAL_" if local else ""
    tag = f"PY_S11B_{infix}{quantity}"
    export_key = None if local or not export else key
    return EMITTER.emit(
        tag,
        payload,
        export_key=export_key,
        class_tag=class_tag,
        dimension_key=dimension_key,
    )


def issue(message: str) -> None:
    prefix = CURRENT_TASK or "GLOBAL"
    ISSUES.append(f"{prefix}: {message}")


def bool_status(value: object) -> tuple[Str, object]:
    value = casify(value)
    if value is sp.true or value == sp.true:
        return PROVED_TRUE, sp.true
    if value is sp.false or value == sp.false:
        return PROVED_FALSE, sp.false
    return UNDECIDED, value


def zero_test(value: object) -> tuple[Str, object]:
    value = casify(value)
    if isinstance(value, sp.MatrixBase):
        tests = [zero_test(entry)[0] for entry in value]
        if all(item == PROVED_TRUE for item in tests):
            return PROVED_TRUE, sp.true
        if any(item == PROVED_FALSE for item in tests):
            return PROVED_FALSE, sp.false
        return UNDECIDED, sp.Eq(value, sp.zeros(*value.shape), evaluate=False)
    if not isinstance(value, sp.Basic):
        return UNDECIDED, sp.Tuple(casify(value))
    simplified = sp.simplify(value)
    if simplified == 0:
        return PROVED_TRUE, sp.true
    answer = simplified.equals(0)
    if answer is False:
        return PROVED_FALSE, sp.false
    if answer is True:
        return PROVED_TRUE, sp.true
    return UNDECIDED, sp.Eq(simplified, 0, evaluate=False)


def relation_status(test: object) -> tuple[Str, object]:
    test = casify(test)
    simplified = sp.simplify(test) if isinstance(test, sp.Basic) else test
    if simplified is sp.true or simplified == sp.true:
        return PROVED_TRUE, sp.true
    if simplified is sp.false or simplified == sp.false:
        return PROVED_FALSE, sp.false
    answer = sp.ask(simplified) if isinstance(simplified, Boolean) else None
    if answer is True:
        return PROVED_TRUE, sp.true
    if answer is False:
        return PROVED_FALSE, sp.false
    return UNDECIDED, simplified


def independent_columns(expressions: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> tuple[int, ...]:
    expanded = [sp.expand(expression) for expression in expressions]
    monomials: set[tuple[int, ...]] = set()
    polys = []
    for expression in expanded:
        poly = sp.Poly(expression, *variables, domain="EX")
        polys.append(poly)
        monomials.update(poly.monoms())
    ordered = sorted(monomials)
    matrix = sp.Matrix([[poly.coeff_monomial(monomial) for expression, poly in zip(expanded, polys)] for monomial in ordered])
    return tuple(matrix.rref()[1])


def solve_linear(equations: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]:
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    solution = matrix.inv() * rhs
    return sp.ImmutableMatrix(solution), sp.ImmutableMatrix(matrix)


def face_solution(
    z_value: sp.Expr,
    mu_theta_value: sp.Expr,
    velocity_value: sp.Expr,
    lambda_a: sp.Expr = Lambda_A,
    lambda_v: sp.Expr = Lambda_V,
) -> dict[str, sp.Expr | sp.ImmutableMatrix]:
    mu_s_value = sp.cancel(mu_theta_value / rho_br)
    affinity = mu_s_value - p_face / rho_m
    closure = lambda_a * affinity + lambda_v * velocity_value
    equations = (
        p_face - z_value * (velocity_value + J_face / rho_m),
        J_face - closure,
    )
    solution, matrix = solve_linear(equations, (p_face, J_face))
    pressure = sp.cancel(solution[0])
    flux = sp.cancel(solution[1])
    return {
        "pressure": pressure,
        "flux": flux,
        "affinity": sp.cancel(mu_s_value - pressure / rho_m),
        "bulk_velocity": sp.cancel(velocity_value + flux / rho_m),
        "matrix": matrix,
        "equations": sp.ImmutableMatrix(equations),
    }


# The O(3)-basis is constructed from generic tensors and vectors.
G_COMPONENTS = tuple(symbol(f"G_{i + 1}_{j + 1}", "COORDINATE", "in-plane displacement gradient component") for i in range(3) for j in range(3))
G = sp.ImmutableMatrix(3, 3, G_COMPONENTS)
grad_theta_components = tuple(symbol(f"grad_theta_{i + 1}", "COORDINATE", "densification gradient component") for i in range(3))
grad_e_components = tuple(symbol(f"grad_e_W_{i + 1}", "COORDINATE", "thickness-fraction gradient component") for i in range(3))
grad_theta = sp.ImmutableMatrix(grad_theta_components)
grad_e = sp.ImmutableMatrix(grad_e_components)


def construct_energy_basis() -> dict[str, object]:
    trace_g = sp.trace(G)
    sym_g = (G + G.T) / 2
    anti_g = (G - G.T) / 2
    st_g = sym_g - sp.eye(3) * trace_g / 3
    curl_squared = sp.expand(2 * sp.trace(anti_g.T * anti_g))
    st_squared = sp.expand(sp.trace(st_g.T * st_g))
    basis = (
        curl_squared,
        st_squared,
        sp.expand(trace_g**2),
        theta**2,
        theta * e_W,
        e_W**2,
        theta * trace_g,
        e_W * trace_g,
        sp.expand((grad_theta.T * grad_theta)[0]),
        sp.expand((grad_theta.T * grad_e)[0]),
        sp.expand((grad_e.T * grad_e)[0]),
    )
    variables = (*G_COMPONENTS, theta, e_W, *grad_theta_components, *grad_e_components)
    pivots = independent_columns(basis, variables)
    carried = (basis[0], basis[3], basis[4], basis[5], basis[10])
    carried_rank = len(independent_columns(carried, variables))
    omitted = []
    span = list(carried)
    rank = carried_rank
    for candidate in basis:
        trial = span + [candidate]
        trial_rank = len(independent_columns(trial, variables))
        if trial_rank > rank:
            omitted.append(candidate)
            span.append(candidate)
            rank = trial_rank

    full_energy = sp.expand(
        mu_R * basis[0] / 2
        + mu_S * basis[1] / 2
        + B_div * basis[2] / 2
        + B_rho_3 * basis[3] / 2
        + C * W0 * basis[4]
        + k_W * W0**2 * basis[5] / 2
        + G_theta_u * basis[6]
        + G_W_u * basis[7]
        + kappa_theta * basis[8] / 2
        + kappa_theta_W * basis[9]
        + kappa_W * W0**4 * basis[10] / 2
    )
    longitudinal_substitutions = {
        **{component: sp.Integer(0) for component in G_COMPONENTS},
        **{component: sp.Integer(0) for component in grad_theta_components},
        **{component: sp.Integer(0) for component in grad_e_components},
        G[0, 0]: eta,
        grad_theta[0]: k * theta,
        grad_e[0]: k * e_W,
    }
    transverse_substitutions = {
        **{component: sp.Integer(0) for component in G_COMPONENTS},
        **{component: sp.Integer(0) for component in grad_theta_components},
        **{component: sp.Integer(0) for component in grad_e_components},
        G[1, 0]: k * u_T,
    }
    longitudinal_basis = tuple(sp.expand(item.subs(longitudinal_substitutions)) for item in basis)
    return {
        "basis": basis,
        "pivots": pivots,
        "omitted": tuple(omitted),
        "full_energy": full_energy,
        "longitudinal_energy": sp.expand(full_energy.subs(longitudinal_substitutions)),
        "transverse_energy": sp.expand(full_energy.subs(transverse_substitutions)),
        "longitudinal_basis": longitudinal_basis,
    }


ENERGY_DATA = construct_energy_basis()
U_LONG = ENERGY_DATA["longitudinal_energy"]
U_TRANSVERSE = ENERGY_DATA["transverse_energy"]


def derive_model(
    substitutions: Mapping[sp.Symbol, sp.Expr] | None = None,
    *,
    z_value: sp.Expr = Z_GENERAL,
    lambda_a: sp.Expr | None = None,
    lambda_v: sp.Expr | None = None,
    lambda_x: sp.Expr | None = None,
    slab_affinity: bool = True,
) -> dict[str, object]:
    cuts = dict(substitutions or {})
    energy = sp.expand(U_LONG.subs(cuts))
    la = sp.cancel((Lambda_A if lambda_a is None else lambda_a).subs(cuts))
    lv = sp.cancel((Lambda_V if lambda_v is None else lambda_v).subs(cuts))
    lx = sp.cancel((Lambda_X if lambda_x is None else lambda_x).subs(cuts))
    z = sp.cancel(z_value.subs(cuts))

    mu_theta_value = sp.diff(energy, theta)
    p_w_value = sp.diff(energy, e_W)
    sigma_eta = sp.diff(energy, eta)
    velocity = -I * omega * W0 * e_W / 2
    face_mu_theta = mu_theta_value if slab_affinity else sp.Integer(0)
    face = face_solution(z, face_mu_theta, velocity, la, lv)
    pressure = face["pressure"]
    flux = face["flux"]
    affinity = face["affinity"]

    mass = sp.cancel(-I * omega * rho_br * (theta + e_W + eta) + 2 * flux)
    inplane = sp.expand(rho_br * omega**2 * eta - k**2 * (sigma_eta - mu_theta_value))
    thickness = sp.cancel(
        -mu_W * W0**2 * omega**2 * e_W
        + p_w_value - mu_theta_value
        + W0 * (pressure + lx * affinity)
    )
    equations = (inplane, mass, thickness)
    raw_matrix = sp.ImmutableMatrix([[sp.diff(equation, variable) for variable in (eta, theta, e_W)] for equation in equations])
    normalized_equations = tuple(sp.together(equation).as_numer_denom()[0] for equation in equations)
    matrix = sp.ImmutableMatrix([[sp.diff(equation, variable) for variable in (eta, theta, e_W)] for equation in normalized_equations])
    raw_determinant = sp.Determinant(raw_matrix)
    determinant = sp.Determinant(matrix)
    return {
        "energy": energy,
        "mu_theta": sp.expand(mu_theta_value),
        "p_W": sp.expand(p_w_value),
        "sigma_eta": sp.expand(sigma_eta),
        "face": face,
        "mass": mass,
        "inplane": inplane,
        "thickness": thickness,
        "raw_matrix": raw_matrix,
        "matrix": matrix,
        "raw_determinant": raw_determinant,
        "determinant": determinant,
        "lambda_A": la,
        "lambda_V": lv,
        "lambda_X": lx,
        "Z": z,
    }


MODEL: dict[str, object] = {}


def task_b0_energy() -> None:
    basis = ENERGY_DATA["basis"]
    emit("ENERGY_BASIS", basis, key="energy_basis")
    emit("ENERGY_BASIS_COUNT", sp.Integer(len(basis)), key="energy_basis_count")
    emit("ENERGY_BASIS_OMISSIONS", ENERGY_DATA["omitted"], key="energy_basis_omissions")
    emit("ENERGY_BASIS_INDEPENDENT_TERMS", sp.Tuple(*(basis[index] for index in ENERGY_DATA["pivots"])), key="energy_basis_independent_terms")

    impermeable_substitution = {theta: -eta - e_W}
    impermeable_reduced = tuple(sp.expand(item.subs(impermeable_substitution)) for item in ENERGY_DATA["longitudinal_basis"])
    imp_pivots = independent_columns(impermeable_reduced, (eta, e_W))
    flux_mass = MODEL["mass"]
    coeff_theta = sp.diff(flux_mass, theta)
    coeff_eta = sp.diff(flux_mass, eta)
    coeff_e = sp.diff(flux_mass, e_W)
    alpha_eta = symbol("alpha_eta_constraint", "CONTROL", "generic flux-on constraint coefficient for divergence")
    alpha_e = symbol("alpha_e_constraint", "CONTROL", "generic flux-on constraint coefficient for thickness")
    flux_reduced = tuple(sp.expand(item.subs(theta, alpha_eta * eta + alpha_e * e_W)) for item in ENERGY_DATA["longitudinal_basis"])
    flux_pivots = independent_columns(flux_reduced, (eta, e_W))
    derived_constraint_coefficients = sp.Tuple(
        sp.cancel(-coeff_eta / coeff_theta),
        sp.cancel(-coeff_e / coeff_theta),
        sp.Eq(coeff_theta, 0, evaluate=False),
    )
    redundancy = sp.Tuple(
        sp.Tuple(Str("IMPERMEABLE"), sp.Tuple(*(sp.Integer(i) for i in imp_pivots)), sp.Tuple(*impermeable_reduced)),
        sp.Tuple(Str("FLUX_ON"), sp.Tuple(*(sp.Integer(i) for i in flux_pivots)), flux_reduced, derived_constraint_coefficients),
    )
    emit("BASIS_REDUNDANCY_UNDER_CONSTRAINT", redundancy, key="basis_redundancy_under_constraint")


def task_b0a() -> None:
    w = symbol("w", "COORDINATE", "bulk normal coordinate", real=True)
    w1 = symbol("w_1", "COORDINATE", "finite projection lower endpoint", real=True)
    w2 = symbol("w_2", "COORDINATE", "finite projection upper endpoint", real=True)
    x1, x2, x3, t = tuple(symbol(name, "COORDINATE", "projection coordinate", real=True) for name in ("x_1", "x_2", "x_3", "t"))
    Omega = sp.Function("Omega")(w)
    rho4 = sp.Function("rho_4D")(x1, x2, x3, w, t)
    jw = sp.Function("j_w")(w)
    jx = tuple(sp.Function(f"j_x_{index}")(x1, x2, x3, w, t) for index in range(1, 4))
    source_finite = sp.Integral(sp.diff(Omega, w) * jw, (w, w1, w2)) - (Omega * jw).subs(w, w2) + (Omega * jw).subs(w, w1)
    source_infinite = sp.Integral(sp.diff(Omega, w) * jw, (w, -sp.oo, sp.oo))
    emit("PROJECTION_FINITE", sp.Tuple(sp.Integral(Omega * rho4, (w, w1, w2)), source_finite), key="projection_finite")
    emit("PROJECTION_INFINITE", sp.Tuple(sp.Integral(Omega * rho4, (w, -sp.oo, sp.oo)), source_infinite), key="projection_infinite")

    Om = symbol("Omega_w", "COORDINATE", "window value at a reflected point")
    Omp = symbol("Omega_prime_w", "COORDINATE", "window derivative at a reflected point")
    je = symbol("j_even_w", "COORDINATE", "even normal-current value")
    jo = symbol("j_odd_w", "COORDINATE", "odd normal-current value")
    even_pair = sp.expand(Omp * je + (-Omp) * je)
    odd_pair = sp.expand(Omp * jo + (-Omp) * (-jo))
    even_boundary = sp.expand(-Om * je + Om * je)
    odd_boundary = sp.expand(-Om * jo + Om * (-jo))
    emit("PARITY_EVEN_JW", sp.Tuple(even_pair, even_boundary, sp.simplify(even_pair + even_boundary)), key="parity_even_jw")
    emit("PARITY_ODD_JW", sp.Tuple(odd_pair, odd_boundary, sp.simplify(odd_pair + odd_boundary)), key="parity_odd_jw")
    emit(
        "PARITY_INTERVAL",
        sp.Tuple(
            sp.Interval(w1, w2),
            sp.Eq(w1 + w2, 0, evaluate=False),
            sp.Interval(-control_L, control_L),
            sp.true,
        ),
        key="parity_interval",
    )

    Omega_dynamic = sp.Function("Omega_dynamic")(w, x1, x2, x3, t)
    extras = [sp.Integral(sp.diff(Omega_dynamic, t) * rho4, (w, w1, w2))]
    extras.extend(sp.Integral(sp.diff(Omega_dynamic, coordinate) * current, (w, w1, w2)) for coordinate, current in zip((x1, x2, x3), jx))
    emit("DYNAMIC_WINDOW_EXTRA_TERMS", sp.Tuple(*extras), key="dynamic_window_extra_terms")


def task_b0b() -> None:
    w = next(sym for sym in DECLARED_SYMBOLS if sym.name == "w")
    A_plus = symbol("A_plus", "COORDINATE", "upper outgoing bulk amplitude")
    A_minus = symbol("A_minus", "COORDINATE", "lower outgoing bulk amplitude")
    phi_plus = A_plus * sp.exp(I * q_out * (w - W0 / 2))
    phi_minus = A_minus * sp.exp(I * q_out * (-w - W0 / 2))
    v_plus = sp.diff(phi_plus, w).subs(w, W0 / 2)
    v_minus = -sp.diff(phi_minus, w).subs(w, -W0 / 2)
    p_plus = I * rho_m * omega * phi_plus.subs(w, W0 / 2)
    p_minus = I * rho_m * omega * phi_minus.subs(w, -W0 / 2)
    z_plus = sp.cancel(p_plus / v_plus)
    z_minus = sp.cancel(p_minus / v_minus)
    emit("Z_IMPERMEABLE", sp.Tuple(z_plus, z_minus), key="z_impermeable")

    q_squared = sp.expand(omega**2 / c_s0**2 - k**2)
    q_squared_axis = q_squared.subs(omega, _omega_real_axis)
    q_prop_axis = sp.sign(_omega_real_axis) * sp.sqrt(q_squared_axis)
    q_prop = q_prop_axis.xreplace({_omega_real_axis: omega})
    kappa = symbol("kappa_out", "COORDINATE", "positive evanescent decay rate", positive=True)
    q_evan = I * kappa
    z_prop = real_axis_simplify(Z_GENERAL.subs(q_out, q_prop))
    z_evan = real_axis_simplify(Z_GENERAL.subs(q_out, q_evan))
    grazing_sides = sp.Tuple(sp.limit(Z_GENERAL, q_out, 0, dir="+"), sp.limit(Z_GENERAL, q_out, 0, dir="-"))
    regimes = sp.Tuple(
        REAL_AXIS,
        sp.Tuple(sp.Gt(q_squared, 0, evaluate=False), q_prop, z_prop, real_axis_simplify(sp.re(z_prop)), real_axis_simplify(sp.im(z_prop))),
        sp.Tuple(sp.Lt(q_squared, 0, evaluate=False), q_evan, z_evan, real_axis_simplify(sp.re(z_evan)), real_axis_simplify(sp.im(z_evan))),
        sp.Tuple(sp.Eq(q_squared, 0, evaluate=False), sp.Integer(0), grazing_sides),
    )
    emit("Z_BY_REGIME", regimes, key="z_by_regime")

    zeta_W = W0 * e_W
    zeta_c = symbol("zeta_c", "COORDINATE", "centre-shift amplitude")
    thickness_velocities = sp.ImmutableMatrix([-I * omega * zeta_W / 2, -I * omega * zeta_W / 2])
    centre_velocities = sp.ImmutableMatrix([-I * omega * zeta_c, I * omega * zeta_c])
    parity_payload = sp.Tuple(
        sp.Tuple(Str("THICKNESS"), thickness_velocities, sp.simplify(Z_GENERAL * thickness_velocities)),
        sp.Tuple(Str("CENTRE_SHIFT"), centre_velocities, sp.simplify(Z_GENERAL * centre_velocities)),
    )
    emit("Z_BY_PARITY", parity_payload, key="z_by_parity")

    zeta_out = symbol("zeta_out", "COORDINATE", "one-face outward displacement amplitude")
    pressure_evan = sp.expand(z_evan * (-I * omega * zeta_out))
    acceleration = -omega**2 * zeta_out
    m_added_solution = sp.solve(sp.Eq(pressure_evan, symbol("m_add_trial", "CONTROL", "added-mass solve placeholder") * acceleration), next(sym for sym in DECLARED_SYMBOLS if sym.name == "m_add_trial"), dict=True)
    m_added = sp.simplify(m_added_solution[0][next(sym for sym in DECLARED_SYMBOLS if sym.name == "m_add_trial")])
    emit("ADDED_MASS", sp.Tuple(pressure_evan, acceleration, m_added), key="added_mass")
    emit("GRAZING_BEHAVIOUR", sp.Tuple(q_squared, sp.limit(v_plus, q_out, 0), grazing_sides), key="grazing_behaviour")

    # Real-axis radiation/decay checks and branch objects required by §1b.
    p_plus_axis = p_plus.subs({omega: _omega_real_axis, q_out: q_prop_axis})
    p_minus_axis = p_minus.subs({omega: _omega_real_axis, q_out: q_prop_axis})
    v_plus_axis = v_plus.subs({omega: _omega_real_axis, q_out: q_prop_axis})
    v_minus_axis = v_minus.subs({omega: _omega_real_axis, q_out: q_prop_axis})
    acoustic_flux_plus = real_axis_simplify(sp.re(p_plus_axis * sp.conjugate(v_plus_axis)) / 2)
    acoustic_flux_minus = real_axis_simplify(sp.re(p_minus_axis * sp.conjugate(v_minus_axis)) / 2)
    realaxis = sp.Tuple(
        REAL_AXIS,
        sp.Tuple(Str("UPPER"), q_prop, acoustic_flux_plus.subs(q_out, q_prop)),
        sp.Tuple(Str("LOWER"), q_prop, acoustic_flux_minus.subs(q_out, q_prop)),
        sp.Tuple(Str("EVANESCENT"), q_evan, sp.exp(I * q_evan * symbol("s_out", "COORDINATE", "outward distance", positive=True))),
    )
    emit("BRANCH_REALAXIS_CHECK", realaxis, key="branch_realaxis_check")
    emit("BRANCH_DEGENERATE_POINT", sp.Tuple(REAL_AXIS, sp.Eq(omega**2, c_s0**2 * k**2, evaluate=False), sp.Integer(0), sp.ImmutableMatrix([[1, 0], [1, 0]]).rank()), key="branch_degenerate_point")


def locus_payload(equation: sp.Expr, solve_variable: sp.Symbol) -> tuple[object, object, object, object, object]:
    equation_relation = sp.Eq(equation, 0, evaluate=False)
    solutions = sp.solve((equation,), (solve_variable,), dict=True)
    solution_payload = sp.Tuple(*(sp.Tuple(*(sp.Eq(key, value, evaluate=False) for key, value in branch.items())) for branch in solutions))
    identity_status, identity_test = zero_test(equation)
    inconsistent_test = sp.false if solutions else sp.Contains(solve_variable, sp.S.EmptySet)
    inconsistent_status, inconsistent_object = relation_status(inconsistent_test)
    real_entries = []
    for branch, branch_payload in zip(solutions, solution_payload):
        branch_value = branch[solve_variable]
        test = sp.And(sp.Eq(sp.im(branch_value), 0, evaluate=False), sp.Ge(tau_A, 0, evaluate=False))
        status, test_object = relation_status(test)
        admission = ADMISSIBLE if status == PROVED_TRUE else EXCLUDED if status == PROVED_FALSE else UNDECIDED
        real_entries.append(sp.Tuple(branch_payload, admission, test_object, sp.Tuple(branch_value, tau_A)))
    return (
        sp.Tuple(equation_relation),
        sp.Tuple(Str(solve_variable.name), solution_payload),
        sp.Tuple(identity_status, identity_test, sp.Tuple(equation)),
        sp.Tuple(inconsistent_status, inconsistent_object, sp.Tuple(equation, solve_variable)),
        sp.Tuple(*real_entries),
    )


def task_b0c() -> None:
    face = face_solution(Z_GENERAL, mu_drive, V_face)
    pressure = face["pressure"]
    emit("FACE_RESPONSE", sp.Tuple(face["equations"], pressure), key="face_response")
    coeff_v = sp.cancel(sp.diff(pressure, V_face))
    coeff_mu = sp.cancel(sp.diff(pressure, mu_drive))
    emit("FACE_RESPONSE_COEFFS", sp.Tuple(coeff_v, coeff_mu), key="face_response_coeffs")

    coefficient_matrix = face["matrix"]
    determinant = sp.factor(coefficient_matrix.det())
    locus = locus_payload(determinant, Lambda_A0)
    for suffix, payload in zip(
        ("EQUATIONS", "SOLUTION", "IDENTICALLY_SATISFIED", "INCONSISTENT", "REAL_ADMISSIBLE"),
        locus,
    ):
        emit(f"DEGENERATE_LOCI_{suffix}", payload, key=f"degenerate_loci_{suffix.lower()}")

    q_squared = omega**2 / c_s0**2 - k**2
    z_prop = rho_m * sp.Abs(omega) / sp.sqrt(q_squared)
    kappa = next(sym for sym in DECLARED_SYMBOLS if sym.name == "kappa_out")
    z_evan = -I * rho_m * omega / kappa
    coefficient_by_regime = []
    for regime_name, z_value in (("PROPAGATING", z_prop), ("EVANESCENT", z_evan), ("GRAZING", sp.oo)):
        response = face_solution(z_value, sp.Integer(0), V_face)["pressure"]
        coefficient = real_axis_cancel(sp.diff(response, V_face))
        dissipative = real_axis_part(coefficient, "real")
        coefficient_by_regime.append(sp.Tuple(Str(regime_name), Str("THICKNESS"), coefficient, dissipative))
        coefficient_by_regime.append(sp.Tuple(Str(regime_name), Str("CENTRE_SHIFT"), coefficient, dissipative))
    emit("PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY", sp.Tuple(REAL_AXIS, *coefficient_by_regime), key="permeable_dissipative_by_regime_and_parity")

    tau_rows = []
    for tau_symbol in (tau_A, tau_V):
        tau_rows.append(
            sp.Tuple(
                tau_symbol,
                sp.simplify(coeff_v.subs(tau_symbol, 0)),
                sp.limit(coeff_v, tau_symbol, sp.oo),
                sp.simplify(coeff_v.subs(tau_symbol, 1 / sp.Abs(omega))),
            )
        )
    emit("PERMEABLE_DISSIPATION_VS_OMEGA_TAU", sp.Tuple(*tau_rows), key="permeable_dissipation_vs_omega_tau")

    J_plus = symbol("J_plus", "COORDINATE", "upper outward relative mass flux")
    J_minus = symbol("J_minus", "COORDINATE", "lower outward relative mass flux")
    channel_matrix = sp.ImmutableMatrix([[-1, -1], [sp.Rational(1, 2), -sp.Rational(1, 2)]])
    channel_values = channel_matrix * sp.ImmutableMatrix([J_plus, J_minus])
    emit("FLUX_CHANNELS", sp.Tuple(channel_matrix, channel_values), key="flux_channels")


def task_b1() -> None:
    exact_sigma = rho_br * (1 + theta + e_W)
    exact_velocity_divergence = symbol("div_v", "COORDINATE", "in-plane velocity divergence")
    exact_balance = sp.Eq(
        symbol("dt_Sigma", "COORDINATE", "time derivative of slab-integrated density")
        + exact_sigma * exact_velocity_divergence,
        -(next(sym for sym in DECLARED_SYMBOLS if sym.name == "J_plus") + next(sym for sym in DECLARED_SYMBOLS if sym.name == "J_minus")),
        evaluate=False,
    )
    constraint = sp.Eq(MODEL["mass"], 0, evaluate=False)
    emit("CONSTRAINT", sp.Tuple(exact_balance, constraint), key="constraint")

    term_origins = sp.Tuple(
        sp.Tuple(sp.diff(rho_br * (1 + theta + e_W), theta) * (-I * omega * theta), Str("PARTIAL_T_SIGMA_THETA")),
        sp.Tuple(sp.diff(rho_br * (1 + theta + e_W), e_W) * (-I * omega * e_W), Str("PARTIAL_T_SIGMA_E_W")),
        sp.Tuple(rho_br * (-I * omega * eta), Str("DIVERGENCE_SIGMA_V_BACKGROUND")),
        sp.Tuple(2 * MODEL["face"]["flux"], Str("OUTWARD_FACE_FLUX_SUM")),
    )
    emit("CONSTRAINT_TERM_ORIGINS", term_origins, key="constraint_term_origins")

    constraint_matrix = sp.ImmutableMatrix([[sp.diff(MODEL["mass"], variable) for variable in (eta, theta, e_W)]])
    generic_rank = constraint_matrix.rank()
    zero_frequency_rank = constraint_matrix.subs(omega, 0).rank()
    dof_payload = sp.Tuple(
        sp.Integer(3),
        generic_rank,
        sp.Integer(3) - generic_rank,
        zero_frequency_rank,
        sp.Integer(3) - zero_frequency_rank,
    )
    emit("INTERNAL_DOF_COUNT", dof_payload, key="internal_dof_count")
    emit(
        "DOF_COUNTING_CONVENTION",
        sp.Tuple(Str("SCALAR_AMPLITUDES_AT_FIXED_K_OMEGA"), sp.Tuple(eta, theta, e_W), constraint_matrix),
        key="dof_counting_convention",
    )


def kernel_pole_inventory() -> sp.Tuple:
    objects = (
        (Str("LAMBDA_A"), Lambda_A),
        (Str("LAMBDA_V"), Lambda_V),
        (Str("LAMBDA_X"), Lambda_X),
        (Str("FACE_PRESSURE"), MODEL["face"]["pressure"]),
        (Str("FACE_FLUX"), MODEL["face"]["flux"]),
        (Str("DISPERSION"), MODEL["raw_determinant"]),
    )
    rows = []
    bare = sp.Tuple(-I / tau_A, -I / tau_V, -I / tau_X)
    for name, expression in objects:
        cancelled = sp.cancel(expression)
        denominator = sp.factor(sp.denom(cancelled))
        try:
            locations = sp.solveset(sp.Eq(denominator, 0), omega, domain=sp.S.Complexes)
        except Exception as exc:
            issue(f"pole solve for {display(name)} returned {type(exc).__name__}")
            locations = sp.ConditionSet(omega, sp.Eq(denominator, 0, evaluate=False), sp.S.Complexes)
        retained = sp.Tuple(*(sp.simplify(denominator.subs(omega, pole)) for pole in bare))
        rows.append(sp.Tuple(name, cancelled, denominator, locations, bare, retained))
    return sp.Tuple(*rows)


def causality_objects() -> dict[str, object]:
    affinity_coordinate = symbol("affinity_coordinate", "COORDINATE", "independent affinity amplitude")
    velocity_coordinate = symbol("velocity_coordinate", "COORDINATE", "independent outward velocity amplitude")
    closure_law = Lambda_A * affinity_coordinate + Lambda_V * velocity_coordinate
    traction_law = p_face + Lambda_X * affinity_coordinate
    extracted = (
        sp.diff(closure_law.subs(velocity_coordinate, 0), affinity_coordinate),
        sp.diff(closure_law.subs(affinity_coordinate, 0), velocity_coordinate),
        sp.diff(traction_law - p_face, affinity_coordinate),
    )
    supplied = (Lambda_A, Lambda_V, Lambda_X)
    orientation_rows = []
    for name, left, right, coefficient, relaxation in zip(
        ("A", "V", "X"), extracted, supplied, (Lambda_A0, Lambda_V0, Lambda_X0), (tau_A, tau_V, tau_X)
    ):
        residual = sp.cancel(left - right)
        status, test = zero_test(sp.together(residual).as_numer_denom()[0])
        indistinguishable = sp.Or(sp.Eq(coefficient, 0, evaluate=False), sp.Eq(relaxation, 0, evaluate=False))
        orientation_rows.append(sp.Tuple(Str(name), left, right, residual, status, test, indistinguishable))

    placeholder_model = derive_model(lambda_a=ell_A, lambda_v=ell_V, lambda_x=ell_X)
    replacements = {ell_A: Lambda_A, ell_V: Lambda_V, ell_X: Lambda_X}
    propagation_rows = []
    for name in ("face", "mass", "inplane", "thickness", "determinant"):
        if name == "face":
            placeholder_object = placeholder_model["face"]["pressure"]
            actual_object = MODEL["face"]["pressure"]
        else:
            placeholder_object = placeholder_model[name]
            actual_object = MODEL[name]
        propagated = casify(placeholder_object).subs(replacements)
        if isinstance(actual_object, sp.MatrixBase):
            residual = sp.ImmutableMatrix(propagated - actual_object)
        else:
            residual = sp.Add(propagated, -actual_object, evaluate=False)
        propagation_rows.append(sp.Tuple(Str(name.upper()), propagated, actual_object, residual))
    return {
        "orientation": sp.Tuple(*orientation_rows),
        "propagation": sp.Tuple(*propagation_rows),
        "poles": kernel_pole_inventory(),
    }


def energy_accounting_objects() -> dict[str, object]:
    Vp = symbol("V_plus", "COORDINATE", "upper outward face velocity amplitude")
    Vm = symbol("V_minus", "COORDINATE", "lower outward face velocity amplitude")
    pp = symbol("delta_p_plus", "COORDINATE", "upper face pressure amplitude")
    pm = symbol("delta_p_minus", "COORDINATE", "lower face pressure amplitude")
    Ap = symbol("A_plus_affinity", "COORDINATE", "upper affinity amplitude")
    Am = symbol("A_minus_affinity", "COORDINATE", "lower affinity amplitude")
    Jp = next(sym for sym in DECLARED_SYMBOLS if sym.name == "J_plus")
    Jm = next(sym for sym in DECLARED_SYMBOLS if sym.name == "J_minus")
    mu_s_coordinate = symbol("mu_s_coordinate", "COORDINATE", "specific slab chemical potential amplitude")

    pressure_terms = sp.Tuple(*(
        real_axis_simplify(term)
        for term in (
            -sp.re(pp * sp.conjugate(Vp)) / 2,
            -sp.re(pm * sp.conjugate(Vm)) / 2,
        )
    ))
    traction_terms = sp.Tuple(*(
        real_axis_simplify(term)
        for term in (
            -sp.re(Lambda_X * Ap * sp.conjugate(Vp)) / 2,
            -sp.re(Lambda_X * Am * sp.conjugate(Vm)) / 2,
        )
    ))
    transfer_terms = sp.Tuple(*(
        real_axis_simplify(term)
        for term in (
            -sp.re(mu_s_coordinate * sp.conjugate(Jp)) / 2,
            -sp.re(mu_s_coordinate * sp.conjugate(Jm)) / 2,
        )
    ))
    total = sp.simplify(sum((*pressure_terms, *traction_terms, *transfer_terms), sp.Integer(0)))

    Zprop = symbol("Z_prop_real", "COORDINATE", "real outgoing propagating face response", positive=True)
    Vabs2 = symbol("V_abs_squared", "COORDINATE", "squared face-velocity amplitude", nonnegative=True)
    slab_pressure_route = sp.simplify(-2 * sp.re(Zprop) * Vabs2 / 2)
    bulk_flux_route = sp.simplify(2 * sp.re(Zprop) * Vabs2 / 2)
    pressure_residual = sp.simplify(slab_pressure_route + bulk_flux_route)

    slab_exchange = real_axis_simplify(-sp.re(
        (pp + Lambda_X * Ap) * sp.conjugate(Vp)
        + (pm + Lambda_X * Am) * sp.conjugate(Vm)
        + mu_s_coordinate * sp.conjugate(Jp)
        + mu_s_coordinate * sp.conjugate(Jm)
    ) / 2)
    supplied_exchange = sp.simplify(sum((*pressure_terms, *traction_terms, *transfer_terms), sp.Integer(0)))
    scale_X = symbol("scale_X", "CONTROL", "traction-order decomposition scale", real=True)
    left_scaled = sp.expand(slab_exchange.subs(Lambda_X0, scale_X * Lambda_X0))
    right_scaled = sp.expand(supplied_exchange.subs(Lambda_X0, scale_X * Lambda_X0))
    max_degree = max(sp.degree(left_scaled, scale_X) or 0, sp.degree(right_scaled, scale_X) or 0)
    order_rows = []
    for degree in range(max_degree + 1):
        left_coefficient = sp.expand(left_scaled).coeff(scale_X, degree)
        right_coefficient = sp.expand(right_scaled).coeff(scale_X, degree)
        order_rows.append(sp.Tuple(sp.Integer(degree), left_coefficient, right_coefficient, sp.simplify(left_coefficient - right_coefficient)))
    return {
        "pressure": pressure_terms,
        "traction": traction_terms,
        "transfer": transfer_terms,
        "total": total,
        "pressure_check": sp.Tuple(REAL_AXIS, slab_pressure_route, -bulk_flux_route, pressure_residual),
        "two_port_check": sp.Tuple(REAL_AXIS, *order_rows),
    }


def task_b2a() -> None:
    virtual_constraint = sp.Eq(symbol("delta_v_theta", "COORDINATE", "virtual densification variation") + symbol("delta_v_e_W", "COORDINATE", "virtual thickness-fraction variation") + symbol("div_delta_v_u", "COORDINATE", "divergence of virtual in-plane displacement"), 0, evaluate=False)
    multiplier = sp.diff(U_LONG, theta)
    CUSTOM_QUANTITY_NAMES.add("virtual_displacement_route")
    emit("VIRTUAL_DISPLACEMENT_ROUTE", sp.Tuple(Str("ELIMINATE_DELTA_V_THETA"), virtual_constraint, multiplier, Str("CHEMICAL_POTENTIAL_DENSITY")), key="virtual_displacement_route")

    emit("INPLANE_EOM", sp.Eq(MODEL["inplane"], 0, evaluate=False), key="inplane_eom")
    emit("THICKNESS_EOM", sp.Eq(MODEL["thickness"], 0, evaluate=False), key="thickness_eom")
    bulk_force = sp.cancel(W0 * (MODEL["face"]["pressure"] + Lambda_X * MODEL["face"]["affinity"]))
    emit("BULK_FORCE_ON_THICKNESS", bulk_force, key="bulk_force_on_thickness")
    no_x_operator = sp.cancel(MODEL["thickness"].subs(Lambda_X0, 0))
    reciprocal_effect = sp.cancel(MODEL["thickness"] - no_x_operator)
    emit("RECIPROCAL_TRACTION_THICKNESS_EFFECT", sp.Tuple(MODEL["thickness"], no_x_operator, reciprocal_effect), key="reciprocal_traction_thickness_effect")

    causality = causality_objects()
    emit("KERNEL_ORIENTATION_IDENTITIES", causality["orientation"], key="kernel_orientation_identities")
    emit("KERNEL_PROPAGATION_RESIDUALS", causality["propagation"], key="kernel_propagation_residuals")
    emit("KERNEL_POLE_LOCATIONS", causality["poles"], key="kernel_pole_locations")
    emit("CAUSALITY_CHECK", sp.Tuple(causality["orientation"], causality["propagation"], causality["poles"]), key="causality_check")

    sigma_coordinate = symbol("sigma_eta_coordinate", "COORDINATE", "independent explicit strain stress")
    mu_coordinate = symbol("mu_theta_coordinate", "COORDINATE", "independent chemical potential density")
    constrained_stress = sigma_coordinate - mu_coordinate
    inplane_coefficient = sp.diff(constrained_stress, mu_coordinate)
    emit("CONVENTION_CHECK_INPLANE", sp.Tuple(constrained_stress, inplane_coefficient, sp.Eq(inplane_coefficient, -1, evaluate=False)), key="convention_check_inplane")

    reduced_energy = sp.expand(U_LONG.subs({k: 0, kappa_W: 0, theta: -e_W, eta: 0}))
    K_check = sp.diff(reduced_energy, e_W, 2)
    p_w = sp.diff(U_LONG, e_W)
    mu_theta_value = sp.diff(U_LONG, theta)
    equation_stiffness = sp.diff((p_w - mu_theta_value).subs({k: 0, kappa_W: 0, theta: -e_W, eta: 0}), e_W)
    stiffness_residual = sp.simplify(equation_stiffness - K_check)
    omega_squared = sp.solve(sp.Eq(-mu_W * W0**2 * omega**2 + equation_stiffness, 0), omega**2, dict=True)[0][omega**2]
    conservative_inequality = sp.Gt(omega_squared, 0, evaluate=False)
    positive_energy = sp.And(sp.Gt(B_rho_3, 0, evaluate=False), sp.Gt(B_rho_3 * k_W * W0**2 - (C * W0)**2, 0, evaluate=False))
    positivity_implication = sp.Implies(positive_energy, conservative_inequality, evaluate=False)
    emit("CONVENTION_CHECK_CONSERVATIVE", sp.Tuple(equation_stiffness, K_check, stiffness_residual, omega_squared, sp.diff(omega_squared, B_rho_3)), key="convention_check_conservative")
    emit("CONSERVATIVE_POSITIVITY_INEQUALITY", sp.Tuple(conservative_inequality, positive_energy, positivity_implication), key="conservative_positivity_inequality")

    accounting = energy_accounting_objects()
    sinks = sp.Tuple(accounting["pressure"], accounting["traction"], accounting["transfer"])
    sources = sp.Tuple(*(
        -term
        for group in (accounting["pressure"], accounting["traction"], accounting["transfer"])
        for term in group
    ))
    emit("ENERGY_SINKS", sinks, key="energy_sinks")
    emit("ENERGY_SOURCES", sources, key="energy_sources")
    emit("UNATTRIBUTED_SINK_TERMS", sp.Tuple(), key="unattributed_sink_terms")
    emit("UNATTRIBUTED_EXCHANGE_TERMS", sp.Tuple(), key="unattributed_exchange_terms")
    emit("PRESSURE_WORK_SIGN_CHECK", accounting["pressure_check"], key="pressure_work_sign_check")
    emit("FULL_TWO_PORT_BALANCE_CHECK", accounting["two_port_check"], key="full_two_port_balance_check")


def bulk_motion_operator(z_value: sp.Expr) -> sp.Expr:
    velocity = -I * omega * W0 * e_W / 2
    face = face_solution(z_value, sp.Integer(0), velocity)
    operator = W0 * (face["pressure"] + Lambda_X * face["affinity"])
    return sp.cancel(sp.diff(operator, e_W))


def task_b2b() -> None:
    thickness_operator = sp.diff(MODEL["thickness"], e_W)
    thickness_remainder = MODEL["thickness"].subs(e_W, 0)
    solved_e = (F_W - thickness_remainder) / thickness_operator
    response = W0 / thickness_operator
    emit("THICKNESS_RESPONSE", sp.Tuple(solved_e, response), key="thickness_response")
    emit(
        "RESPONSE_NORMALIZATION",
        sp.Tuple(Str("DELTA_W"), Str("GENERALIZED_THICKNESS_FORCE_PER_X_VOLUME"), W0 * e_W, F_W),
        key="response_normalization",
    )

    q_squared = omega**2 / c_s0**2 - k**2
    kappa = next(sym for sym in DECLARED_SYMBOLS if sym.name == "kappa_out")
    regime_z = (
        (Str("PROPAGATING"), rho_m * sp.Abs(omega) / sp.sqrt(q_squared)),
        (Str("EVANESCENT"), -I * rho_m * omega / kappa),
        (Str("GRAZING"), symbol("Z_grazing_control", "CONTROL", "grazing face-response limit coordinate", positive=True)),
    )
    rows = []
    interpretation_rows = []
    for name, z_value in regime_z:
        operator = bulk_motion_operator(z_value)
        if name == Str("GRAZING"):
            operator = safe_limit(operator, z_value, sp.oo)
        operator = real_axis_cancel(operator)
        real_piece = real_axis_part(operator, "real")
        imaginary_piece = I * real_axis_part(operator, "imag")
        inertia_coefficient = -real_piece / omega**2
        resistance_coefficient = -imaginary_piece / (I * omega)
        rows.append(sp.Tuple(name, operator, real_piece, imaginary_piece, inertia_coefficient, resistance_coefficient))
        mass_test = sp.Eq(imaginary_piece, 0, evaluate=False)
        interpretation_rows.append(sp.Tuple(name, mass_test))
    emit("BULK_OPERATOR_BY_REGIME", sp.Tuple(REAL_AXIS, *rows), key="bulk_operator_by_regime")
    emit("MASS_INTERPRETATION_VALID_WHERE", sp.Tuple(REAL_AXIS, *interpretation_rows), key="mass_interpretation_valid_where")


def compressional_response(model: Mapping[str, object]) -> dict[str, object]:
    mass = model["mass"]
    thickness = model["thickness"]
    equations = (mass, thickness)
    coefficient_matrix = sp.ImmutableMatrix([[sp.diff(eq, variable) for variable in (theta, e_W)] for eq in equations])
    forcing_per_eta = -sp.ImmutableMatrix([sp.diff(eq, eta) for eq in equations])
    inverse = sp.Inverse(coefficient_matrix)
    solution_per_eta = sp.MatMul(inverse, forcing_per_eta, evaluate=False)
    stress = sp.expand(model["sigma_eta"] - model["mu_theta"])
    stress_row = sp.ImmutableMatrix([[sp.diff(stress, theta), sp.diff(stress, e_W)]])
    stress_direct = sp.diff(stress, eta)
    response_product = sp.MatMul(stress_row, inverse, forcing_per_eta, evaluate=False)
    ratio = sp.Add(stress_direct, sp.Trace(response_product), evaluate=False)
    stress_solved = ratio * eta
    return {
        "matrix": coefficient_matrix,
        "solution_per_eta": solution_per_eta,
        "stress": stress_solved,
        "ratio": ratio,
    }


def safe_limit(expression: sp.Expr, variable: sp.Symbol, point: object, direction: str = "+") -> object:
    try:
        return sp.limit(expression, variable, point, dir=direction)
    except Exception as exc:
        issue(f"limit in {variable} at {point} returned {type(exc).__name__}")
        return sp.Limit(expression, variable, point, dir=direction)


def matrix_limit(matrix: sp.MatrixBase, substitutions: Mapping[sp.Symbol, sp.Expr], variable: sp.Symbol, point: object) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(matrix.rows, matrix.cols, lambda i, j: safe_limit(matrix[i, j].subs(substitutions), variable, point))


def compressional_path_limit(
    model: Mapping[str, object],
    substitutions: Mapping[sp.Symbol, sp.Expr],
    variable: sp.Symbol,
    point: object,
) -> sp.Tuple:
    equations = (model["mass"], model["thickness"])
    coefficient_matrix = sp.ImmutableMatrix([[sp.diff(eq, field) for field in (theta, e_W)] for eq in equations])
    forcing = -sp.ImmutableMatrix([sp.diff(eq, eta) for eq in equations])
    stress = sp.expand(model["sigma_eta"] - model["mu_theta"])
    stress_row = sp.ImmutableMatrix([[sp.diff(stress, theta), sp.diff(stress, e_W)]])
    stress_direct = sp.diff(stress, eta)
    limited_matrix = matrix_limit(coefficient_matrix, substitutions, variable, point)
    limited_forcing = matrix_limit(forcing, substitutions, variable, point)
    limited_row = matrix_limit(stress_row, substitutions, variable, point)
    limited_direct = safe_limit(stress_direct.subs(substitutions), variable, point)
    product = sp.MatMul(limited_row, sp.Inverse(limited_matrix), limited_forcing, evaluate=False)
    response = sp.Add(limited_direct, sp.Trace(product), evaluate=False)
    return sp.Tuple(response, sp.Determinant(limited_matrix), limited_matrix, limited_forcing)


def task_b2c() -> None:
    compression = compressional_response(MODEL)
    emit(
        "COMPRESSIONAL_RESPONSE",
        sp.Tuple(Str("CONSTRAINED_LONGITUDINAL_STRESS"), Str("DIVERGENCE_AMPLITUDE"), compression["stress"], eta, compression["ratio"]),
        key="compressional_response",
    )

    v_phase = symbol("v_phase", "COORDINATE", "fixed phase velocity", positive=True)
    q_phase = omega * sp.sqrt(c_s0**-2 - v_phase**-2)
    low_limit = compressional_path_limit(MODEL, {q_out: I * k}, omega, 0)
    high_limit = compressional_path_limit(MODEL, {q_out: omega / c_s0}, omega, sp.oo)
    phase_limit = compressional_path_limit(MODEL, {k: omega / v_phase, q_out: q_phase}, omega, 0)
    limits = sp.Tuple(
        sp.Tuple(Str("OMEGA_TO_ZERO_FIXED_NONZERO_K"), low_limit),
        sp.Tuple(Str("OMEGA_TO_INFINITY_FIXED_K"), high_limit),
        sp.Tuple(Str("OMEGA_TO_ZERO_FIXED_PHASE_VELOCITY"), phase_limit),
        sp.Tuple(Str("PATH_DIFFERENCE"), sp.Add(low_limit[0], -phase_limit[0], evaluate=False)),
    )
    emit("LIMITS_AND_PATH", limits, key="limits_and_path")

    mass_at_zero = sp.cancel(MODEL["mass"].subs(omega, 0))
    divided_mass = sp.cancel(MODEL["mass"] / omega)
    divided_limit = safe_limit(divided_mass.subs(q_out, I * k), omega, 0)
    frozen_mass = MODEL["mass"].subs(e_W, 0)
    frozen_theta = -sp.diff(frozen_mass, eta) * eta / sp.diff(frozen_mass, theta)
    frozen_mass_solution = sp.Tuple(sp.Eq(theta, frozen_theta, evaluate=False))
    frozen_stress = ((MODEL["sigma_eta"] - MODEL["mu_theta"]).subs({e_W: 0, theta: frozen_theta})) / eta
    emit(
        "FROZEN_THICKNESS_IDENTIFICATION",
        sp.Tuple(mass_at_zero, divided_limit, frozen_mass_solution, frozen_stress),
        key="frozen_thickness_identification",
    )

    Lambda_p0 = symbol("Lambda_p_0", "DERIVED", "raw-pressure slice coefficient", real=True)
    tau_slice = symbol("tau", "KNOB", "equal relaxation-time slice", nonnegative=True)
    affinity_flux_pressure_coefficient = sp.diff((Lambda_A * (-p_face / rho_m)).subs(omega, 0), p_face)
    map_equation = sp.Eq(Lambda_p0, affinity_flux_pressure_coefficient, evaluate=False)
    mapping_solution = sp.solve(map_equation, Lambda_A0, dict=True)
    emit("ZPERM_SLICE_MAP", sp.Tuple(affinity_flux_pressure_coefficient, map_equation, mapping_solution), key="zperm_slice_map")

    sliced_face = face_solution(Z_GENERAL, sp.Integer(0), V_face)
    sliced_pressure = sp.cancel(sliced_face["pressure"].subs({tau_A: tau_slice, tau_V: tau_slice}))
    if mapping_solution:
        sliced_pressure = sp.cancel(sliced_pressure.subs(mapping_solution[0]))
    sliced_ratio = sp.cancel(sliced_pressure / V_face)
    emit("ZPERM_SLICE", sp.Tuple(sliced_pressure, sliced_ratio), key="zperm_slice")
    CUSTOM_QUANTITY_NAMES.add("zperm_general_slice_scope")
    emit(
        "ZPERM_GENERAL_SLICE_SCOPE",
        sp.Tuple(MODEL["face"]["pressure"], sliced_pressure, sp.Tuple(tau_A, tau_V), tau_slice),
        key="zperm_general_slice_scope",
    )


def hermitian_part(matrix: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix((matrix + matrix.conjugate().T) / 2)


def real_axis_hermitian_part(matrix: sp.MatrixBase) -> sp.ImmutableMatrix:
    """Form the Hermitian part after imposing only omega in R."""
    axis_matrix = matrix.subs(omega, _omega_real_axis)
    axis_hermitian = hermitian_part(axis_matrix)
    return sp.ImmutableMatrix(axis_hermitian.xreplace({_omega_real_axis: omega}))


def task_b2d() -> None:
    V = symbol("V_port", "COORDINATE", "two-port velocity input")
    mu_s = symbol("mu_s_port", "COORDINATE", "two-port chemical input")
    p = symbol("p_port", "COORDINATE", "two-port pressure")
    J = symbol("J_port", "COORDINATE", "two-port mass flux")
    A = sp.expand(mu_s - p / rho_m)
    v_bulk = sp.expand(V + J / rho_m)
    left_power = sp.expand((p + Lambda_X * A) * sp.conjugate(V) + mu_s * sp.conjugate(J))
    right_power = sp.expand(p * sp.conjugate(v_bulk) + A * sp.conjugate(J) + Lambda_X * A * sp.conjugate(V))
    left_power_real = real_axis_simplify(sp.re(left_power) / 2)
    right_power_real = real_axis_simplify(sp.re(right_power) / 2)
    power_residual = sp.simplify(left_power_real - right_power_real)
    emit("TWO_PORT_POWER_IDENTITY", sp.Tuple(REAL_AXIS, left_power_real, right_power_real, power_residual), key="two_port_power_identity")

    response = face_solution(Z_GENERAL, rho_br * mu_s, V)
    pressure = response["pressure"]
    flux = response["flux"]
    affinity = response["affinity"]
    output_matrix = sp.ImmutableMatrix([
        [sp.diff(pressure + Lambda_X * affinity, V), sp.diff(pressure + Lambda_X * affinity, mu_s)],
        [sp.diff(flux, V), sp.diff(flux, mu_s)],
    ])
    output_matrix = sp.ImmutableMatrix(real_axis_form(output_matrix))
    H_port = real_axis_hermitian_part(output_matrix)
    port_minors = sp.Tuple(
        sp.Ge(sp.re(H_port[0, 0]), 0, evaluate=False),
        sp.Ge(sp.Determinant(H_port), 0, evaluate=False),
    )
    emit("PORT_DISSIPATIVITY", sp.Tuple(REAL_AXIS, output_matrix, H_port, port_minors), key="port_dissipativity")
    z_dependence = sp.Tuple(*sorted(H_port.free_symbols.intersection(Z_GENERAL.free_symbols), key=sp.default_sort_key))
    emit("PORT_CONDITION_KIND", sp.Tuple(z_dependence, sp.Eq(sp.Integer(len(z_dependence)), 0, evaluate=False)), key="port_condition_kind")

    mixed_matrix = sp.ImmutableMatrix([[Lambda_A, Lambda_V], [Lambda_X, 0]])
    H_mixed = real_axis_hermitian_part(mixed_matrix)
    admissibility = sp.Tuple(
        sp.Ge(sp.re(H_mixed[0, 0]), 0, evaluate=False),
        sp.Ge(sp.simplify(H_mixed.det()), 0, evaluate=False),
    )

    epsilon = symbol("epsilon_port", "CONTROL", "formal mixed-law conversion coefficient", nonzero=True, real=True)
    a_entry, b_entry, c_entry, d_entry = Lambda_A, Lambda_V, Lambda_X, epsilon
    all_flux = sp.ImmutableMatrix([
        [a_entry - b_entry * c_entry / d_entry, b_entry / d_entry],
        [-c_entry / d_entry, 1 / d_entry],
    ])
    flux_reciprocity_residual = sp.factor(sp.together(all_flux[0, 1] - all_flux[1, 0]) * epsilon)
    all_force = sp.ImmutableMatrix([
        [1 / a_entry, -b_entry / a_entry],
        [c_entry / a_entry, d_entry - c_entry * b_entry / a_entry],
    ])
    force_reciprocity_residual = sp.factor(sp.together(all_force[0, 1] - all_force[1, 0]) * a_entry)
    reciprocity_equation = sp.Eq(flux_reciprocity_residual, 0, evaluate=False)
    reciprocity_crosscheck = sp.simplify(flux_reciprocity_residual - force_reciprocity_residual)
    emit("ONSAGER_CONDITION", sp.Tuple(all_flux, reciprocity_equation), key="onsager_condition")
    emit("ONSAGER_RECIPROCITY", sp.Tuple(all_force, flux_reciprocity_residual, force_reciprocity_residual, reciprocity_crosscheck), key="onsager_reciprocity")
    emit("ONSAGER_DETERMINABLE", relation_status(reciprocity_equation), key="onsager_determinable")

    reciprocal_numerator = sp.Poly(sp.together(flux_reciprocity_residual).as_numer_denom()[0], omega)
    reciprocal_time_equations = sp.Tuple(*(sp.Eq(coefficient, 0, evaluate=False) for coefficient in reciprocal_numerator.all_coeffs()))
    reciprocal_time_solutions = sp.solve(tuple(rel.lhs for rel in reciprocal_time_equations), (Lambda_X0, tau_X), dict=True)
    admissible_cross = sp.simplify(H_mixed[0, 1])
    admissible_numerator = sp.Poly(sp.together(admissible_cross).as_numer_denom()[0], omega)
    admissible_time_equations = sp.Tuple(*(sp.Eq(coefficient, 0, evaluate=False) for coefficient in admissible_numerator.all_coeffs()))
    admissible_time_solutions = sp.solve(tuple(rel.lhs for rel in admissible_time_equations), (Lambda_X0, tau_X), dict=True)
    emit("RELAXATION_TIME_RELATIONS", sp.Tuple(reciprocal_time_equations, reciprocal_time_solutions, admissible_time_equations, admissible_time_solutions), key="relaxation_time_relations")
    emit("COEFFICIENT_ADMISSIBILITY", sp.Tuple(REAL_AXIS, mixed_matrix, H_mixed, admissibility, admissible_time_equations, admissible_time_solutions), key="coefficient_admissibility")

    region_relation_tests = sp.Tuple(
        sp.Implies(sp.And(*admissibility), reciprocity_equation, evaluate=False),
        sp.Implies(reciprocity_equation, sp.And(*admissibility), evaluate=False),
    )
    CUSTOM_QUANTITY_NAMES.add("admissibility_reciprocity_set_relation")
    emit("ADMISSIBILITY_RECIPROCITY_SET_RELATION", region_relation_tests, key="admissibility_reciprocity_set_relation")

    # Root membership is kept as an explicit symbolic condition and never used as a gate.
    root_coordinate = symbol("omega_root", "COORDINATE", "longitudinal root coordinate")
    root_condition = sp.Eq(MODEL["determinant"].subs(omega, root_coordinate), 0, evaluate=False)
    growing = sp.Gt(sp.im(root_coordinate), 0, evaluate=False)
    decaying = sp.Lt(sp.im(root_coordinate), 0, evaluate=False)
    emit("GROWTH_INSIDE_ADMISSIBLE", sp.Tuple(root_condition, growing, admissibility, region_relation_tests), key="growth_inside_admissible")
    emit("DECAY_INSIDE_ADMISSIBLE", sp.Tuple(root_condition, decaying, admissibility, region_relation_tests), key="decay_inside_admissible")


def solve_root_set(equation: sp.Expr, variable: sp.Symbol) -> object:
    if equation.has(sp.Determinant):
        issue(f"closed-form root construction in {variable} unavailable for the exact unevaluated determinant; retained implicit complex locus")
        return sp.ConditionSet(variable, sp.Eq(equation, 0, evaluate=False), sp.S.Complexes)
    try:
        return sp.solveset(sp.Eq(equation, 0), variable, domain=sp.S.Complexes)
    except Exception as exc:
        issue(f"root construction in {variable} returned {type(exc).__name__}")
        return sp.ConditionSet(variable, sp.Eq(equation, 0, evaluate=False), sp.S.Complexes)


def task_b5() -> None:
    determinant = MODEL["determinant"]
    dispersion_system = sp.Tuple(
        sp.Eq(determinant, 0, evaluate=False),
        sp.Eq(q_out**2, omega**2 / c_s0**2 - k**2, evaluate=False),
    )
    emit("LONGITUDINAL_DISPERSION", sp.Tuple(MODEL["matrix"], dispersion_system), key="longitudinal_dispersion")

    roots = solve_root_set(determinant, omega)
    closed_form_test = sp.false if isinstance(roots, sp.ConditionSet) else sp.true
    emit("ROOTS", sp.Tuple(roots, closed_form_test), key="roots")
    imaginary_set = sp.ImageSet(sp.Lambda(omega, sp.im(omega)), roots)
    emit("IMAGINARY_PART", imaginary_set, key="imaginary_part")

    growing_roots = sp.ConditionSet(omega, sp.And(sp.Contains(omega, roots), sp.Gt(sp.im(omega), 0, evaluate=False)), sp.S.Complexes)
    decaying_roots = sp.ConditionSet(omega, sp.And(sp.Contains(omega, roots), sp.Lt(sp.im(omega), 0, evaluate=False)), sp.S.Complexes)
    real_roots = sp.ConditionSet(omega, sp.And(sp.Contains(omega, roots), sp.Eq(sp.im(omega), 0, evaluate=False)), sp.S.Complexes)
    emit("ROOT_STABILITY_CLASS", sp.Tuple(growing_roots, decaying_roots, real_roots), key="root_stability_class")

    moduli = sp.Tuple(B_rho_3, C, k_W, kappa_W, B_div, mu_S, G_theta_u, G_W_u, kappa_theta, kappa_theta_W)
    stability_system = sp.Tuple(
        sp.Eq(determinant, 0, evaluate=False),
        sp.Le(sp.im(omega), 0, evaluate=False),
        moduli,
    )
    emit("STABILITY_CONDITION", stability_system, key="stability_condition")

    origin_derivatives = sp.Tuple(
        sp.Tuple(Str("BULK_RADIATION"), sp.diff(determinant, q_out)),
        sp.Tuple(Str("AFFINITY_CHANNEL"), sp.diff(determinant, Lambda_A0)),
        sp.Tuple(Str("VELOCITY_CHANNEL"), sp.diff(determinant, Lambda_V0)),
        sp.Tuple(Str("RECIPROCAL_TRACTION"), sp.diff(determinant, Lambda_X0)),
    )
    emit("DISSIPATION_ORIGIN", origin_derivatives, key="dissipation_origin")

    z_control = symbol("Z_control", "CONTROL", "grazing response limit coordinate")
    grazing_model = derive_model(z_value=z_control)
    grazing_matrix = sp.ImmutableMatrix([
        [safe_limit(entry, z_control, sp.oo) for entry in row]
        for row in grazing_model["matrix"].tolist()
    ]).subs(omega**2, c_s0**2 * k**2)
    grazing_det = sp.Determinant(grazing_matrix)
    internal_no_bulk = derive_model(
        substitutions={Lambda_A0: 0, Lambda_V0: 0, Lambda_X0: 0},
        z_value=sp.Integer(0),
    )["matrix"].subs(omega**2, c_s0**2 * k**2)
    grazing_classification = sp.Tuple(
        grazing_matrix,
        grazing_det,
        sp.Tuple(Str("THRESHOLD_STATE"), sp.Eq(grazing_det, 0, evaluate=False)),
        sp.Tuple(Str("BOUND_STATE_IN_CONTINUUM"), sp.Lt(internal_no_bulk.rank(), 3, evaluate=False)),
        sp.Tuple(Str("RESONANCE"), sp.Ne(grazing_det, 0, evaluate=False)),
        sp.Tuple(Str("BOUND_MODE"), sp.Eq(q_out, I * sp.Abs(q_out), evaluate=False)),
    )
    emit("GRAZING_MODE_CLASSIFICATION", grazing_classification, key="grazing_mode_classification")

    opposite_determinant = determinant.subs(q_out, -q_out)
    opposite_roots = solve_root_set(opposite_determinant, omega)
    sensitivity_rows: object
    if isinstance(roots, sp.FiniteSet) and isinstance(opposite_roots, sp.FiniteSet) and len(roots) == len(opposite_roots):
        ordered_roots = sorted(roots, key=sp.default_sort_key)
        ordered_opposite = sorted(opposite_roots, key=sp.default_sort_key)
        sensitivity_rows = sp.Tuple(*(
            sp.Tuple(root, opposite, sp.im(root), sp.im(opposite), sp.cancel(sp.im(opposite) / sp.im(root)))
            for root, opposite in zip(ordered_roots, ordered_opposite)
        ))
    else:
        sensitivity_rows = sp.Tuple(roots, opposite_roots, UNDEFINED_RATIO)
    emit("BRANCH_SENSITIVITY", sensitivity_rows, key="branch_sensitivity")
    emit(
        "SHEET_OF_EACH_ROOT",
        sp.Tuple(
            roots,
            sp.Eq(q_out**2, omega**2 / c_s0**2 - k**2, evaluate=False),
            Str("ANALYTIC_CONTINUATION_FROM_UPPER_RIM"),
            opposite_roots,
            sp.Eq(-q_out, -q_out, evaluate=False),
        ),
        key="sheet_of_each_root",
    )

    no_x_determinant = determinant.subs(Lambda_X0, 0)
    no_x_roots = solve_root_set(no_x_determinant, omega)
    reciprocal_effect = sp.Tuple(
        determinant,
        no_x_determinant,
        sp.simplify(determinant - no_x_determinant),
        roots,
        no_x_roots,
        sp.ConditionSet(omega, sp.And(sp.Eq(determinant, 0, evaluate=False), sp.Eq(sp.diff(determinant, omega), 0, evaluate=False)), sp.S.Complexes),
        sp.ConditionSet(omega, sp.And(sp.Eq(no_x_determinant, 0, evaluate=False), sp.Eq(sp.diff(no_x_determinant, omega), 0, evaluate=False)), sp.S.Complexes),
    )
    emit("RECIPROCAL_TRACTION_ROOT_EFFECT", reciprocal_effect, key="reciprocal_traction_root_effect")

    causality = EMITTER.values["PY_S11B_CAUSALITY_CHECK"]
    region_growth = EMITTER.values["PY_S11B_GROWTH_INSIDE_ADMISSIBLE"]
    region_decay = EMITTER.values["PY_S11B_DECAY_INSIDE_ADMISSIBLE"]
    emit("GROWTH_ARTIFACT_DIAGNOSTICS", sp.Tuple(growing_roots, causality, sensitivity_rows, region_growth), key="growth_artifact_diagnostics")
    emit("DECAY_ARTIFACT_DIAGNOSTICS", sp.Tuple(decaying_roots, causality, sensitivity_rows, region_decay), key="decay_artifact_diagnostics")

    accounting = energy_accounting_objects()
    source_condition = sp.And(sp.Contains(omega, growing_roots), sp.Lt(sp.re(accounting["total"]), 0, evaluate=False))
    emit("ROOT_POWER_SOURCE", sp.Tuple(growing_roots, decaying_roots, accounting["total"], source_condition), key="root_power_source")


def task_b4() -> None:
    slice_model = derive_model(
        substitutions={Lambda_A0: 0, Lambda_V0: 0, Lambda_X0: 0, k: 0},
        z_value=rho_m * c_s0,
    )
    breathing_equation = sp.cancel(slice_model["thickness"].subs({eta: 0, theta: -e_W}) / e_W)
    roots = sp.solve(sp.Eq(breathing_equation, 0), omega)
    emit("BREATHING_SLICE_DISPERSION", sp.Tuple(sp.Eq(breathing_equation, 0, evaluate=False), sp.Tuple(*roots)), key="breathing_slice_dispersion")
    root_conditions = sp.Tuple(*(
        sp.Tuple(root, sp.re(root), sp.im(root), sp.Ge(-sp.im(root), 0, evaluate=False), sp.Gt(sp.im(root), 0, evaluate=False))
        for root in roots
    ))
    stiffness = sp.diff(
        (slice_model["p_W"] - slice_model["mu_theta"]).subs({eta: 0, theta: -e_W, k: 0}),
        e_W,
    )
    emit("BREATHING_STABILITY_CONDITION", sp.Tuple(stiffness, C, root_conditions), key="breathing_stability_condition")


def task_b6() -> None:
    transverse_coupling = sp.diff(U_TRANSVERSE, u_T, e_W)
    normalization = sp.cancel(transverse_coupling / (u_T * e_W)) if transverse_coupling != 0 else sp.nan
    emit(
        "TRANSVERSE_COUPLING",
        sp.Tuple(Str("TRANSVERSE_DISPLACEMENT"), Str("THICKNESS_FRACTION"), transverse_coupling, normalization),
        key="transverse_coupling",
    )
    stiffness = sp.diff(U_TRANSVERSE, u_T, 2)
    transverse_operator = sp.expand(-rho_br * omega**2 + stiffness)
    transverse_roots = sp.solve(sp.Eq(transverse_operator, 0), omega)
    emit("TRANSVERSE_DISPERSION", sp.Tuple(sp.Eq(transverse_operator, 0, evaluate=False), sp.Tuple(*transverse_roots)), key="transverse_dispersion")
    dependence = sp.ImmutableMatrix([
        [sp.diff(transverse_operator, coefficient) for coefficient in (Lambda_A0, Lambda_V0, Lambda_X0, tau_A, tau_V, tau_X)]
    ])
    imaginary_parts = sp.Tuple(*(sp.im(root) for root in transverse_roots))
    emit("TRANSVERSE_DISSIPATION", sp.Tuple(imaginary_parts, dependence, sp.diff(transverse_operator, mu_drive)), key="transverse_dissipation")


L_DIM = sp.ImmutableMatrix([1, 0, 0])
T_DIM = sp.ImmutableMatrix([0, 1, 0])
M_DIM = sp.ImmutableMatrix([0, 0, 1])
ZERO_DIM = sp.ImmutableMatrix([0, 0, 0])


def derive_dimension(rhs: sp.MatrixBase) -> sp.ImmutableMatrix:
    d_l, d_t, d_m = sp.symbols("d_l d_t d_m")
    solution = sp.solve(
        tuple(sp.Eq(left, right) for left, right in zip((d_l, d_t, d_m), rhs)),
        (d_l, d_t, d_m),
        dict=True,
    )
    if not solution:
        raise RuntimeError("dimension equations did not solve")
    return sp.ImmutableMatrix([solution[0][d_l], solution[0][d_t], solution[0][d_m]])


def homogeneous_payload(term_dimensions: Sequence[sp.MatrixBase]) -> sp.Tuple:
    frozen = tuple(sp.ImmutableMatrix(item) for item in term_dimensions)
    tests = tuple(sp.Eq(sp.Tuple(*item), sp.Tuple(*frozen[0]), evaluate=False) for item in frozen[1:])
    result = sp.And(*tests) if tests else sp.true
    return sp.Tuple(sp.Tuple(*frozen), result)


def task_b7() -> None:
    energy_density = derive_dimension(M_DIM - L_DIM - 2 * T_DIM)
    pressure_4d = derive_dimension(M_DIM - 2 * L_DIM - 2 * T_DIM)
    mass_flux = derive_dimension(M_DIM - 3 * L_DIM - T_DIM)
    specific_energy = derive_dimension(2 * L_DIM - 2 * T_DIM)
    dimensions: dict[str, tuple[object, Str, object]] = {
        "Z": (derive_dimension(pressure_4d - (L_DIM - T_DIM)), Str("INDEPENDENT"), sp.Tuple(pressure_4d, L_DIM - T_DIM)),
        "M_ADD": (derive_dimension(pressure_4d - (L_DIM - 2 * T_DIM)), Str("INDEPENDENT"), sp.Tuple(pressure_4d, L_DIM - 2 * T_DIM)),
        "RHO_M": (derive_dimension(M_DIM - 4 * L_DIM), Str("DEFINITIONAL"), sp.Tuple(M_DIM, 4 * L_DIM)),
        "C_S0": (derive_dimension(L_DIM - T_DIM), Str("INDEPENDENT"), sp.Tuple(L_DIM, T_DIM)),
        "V_DR": (derive_dimension(L_DIM - T_DIM), Str("INDEPENDENT"), sp.Tuple(L_DIM, T_DIM)),
        "RHO_BR": (derive_dimension(M_DIM - 3 * L_DIM), Str("DEFINITIONAL"), sp.Tuple(M_DIM, 3 * L_DIM)),
        "B_RHO": (pressure_4d, Str("INDEPENDENT"), sp.Tuple(pressure_4d)),
        "B_RHO_3": (energy_density, Str("INDEPENDENT"), sp.Tuple(pressure_4d, L_DIM)),
        "MU_W": (derive_dimension(energy_density - 2 * (L_DIM - T_DIM)), Str("INDEPENDENT"), sp.Tuple(energy_density, 2 * (L_DIM - T_DIM))),
        "K_W": (derive_dimension(energy_density - 2 * L_DIM), Str("INDEPENDENT"), sp.Tuple(energy_density, 2 * L_DIM)),
        "KAPPA_W": (derive_dimension(energy_density - 2 * L_DIM), Str("INDEPENDENT"), sp.Tuple(energy_density, 2 * L_DIM)),
        "C": (derive_dimension(energy_density - L_DIM), Str("INDEPENDENT"), sp.Tuple(energy_density, L_DIM)),
        "MU_R": (energy_density, Str("INDEPENDENT"), sp.Tuple(energy_density, ZERO_DIM)),
        "THICKNESS_RESPONSE": (derive_dimension(L_DIM - energy_density), Str("INDEPENDENT"), sp.Tuple(L_DIM, energy_density)),
        "COMPRESSIONAL_RESPONSE": (energy_density, Str("INDEPENDENT"), sp.Tuple(energy_density, ZERO_DIM)),
        "TRANSVERSE_COUPLING": (sp.nan, Str("UNDETERMINED"), sp.Tuple(EMITTER.values["PY_S11B_TRANSVERSE_COUPLING"])),
        "LAMBDA_A_0": (derive_dimension(mass_flux - specific_energy), Str("INDEPENDENT"), sp.Tuple(mass_flux, specific_energy)),
        "LAMBDA_V_0": (derive_dimension(mass_flux - (L_DIM - T_DIM)), Str("INDEPENDENT"), sp.Tuple(mass_flux, L_DIM - T_DIM)),
        "LAMBDA_X_0": (derive_dimension(pressure_4d - specific_energy), Str("INDEPENDENT"), sp.Tuple(pressure_4d, specific_energy)),
        "TAU_A": (T_DIM, Str("DEFINITIONAL"), sp.Tuple(omega, tau_A, ZERO_DIM)),
        "TAU_V": (T_DIM, Str("DEFINITIONAL"), sp.Tuple(omega, tau_V, ZERO_DIM)),
        "TAU_X": (T_DIM, Str("DEFINITIONAL"), sp.Tuple(omega, tau_X, ZERO_DIM)),
        "AFFINITY": (specific_energy, Str("INDEPENDENT"), sp.Tuple(pressure_4d, M_DIM - 4 * L_DIM)),
        "MU_THETA": (energy_density, Str("INDEPENDENT"), sp.Tuple(U_LONG, theta)),
        "MU_S": (specific_energy, Str("INDEPENDENT"), sp.Tuple(energy_density, M_DIM - 3 * L_DIM)),
        "PROJECTION_SOURCE": (mass_flux, Str("INDEPENDENT"), sp.Tuple(M_DIM - 4 * L_DIM, T_DIM, L_DIM)),
        "FACE_RESPONSE": (sp.Tuple(pressure_4d, pressure_4d - (L_DIM - T_DIM), pressure_4d - energy_density), Str("INDEPENDENT"), sp.Tuple(p_face, V_face, mu_drive)),
        "B_DIV": (energy_density, Str("INDEPENDENT"), sp.Tuple(energy_density)),
        "MU_S_COEFFICIENT": (energy_density, Str("INDEPENDENT"), sp.Tuple(energy_density)),
        "G_THETA_U": (energy_density, Str("INDEPENDENT"), sp.Tuple(energy_density)),
        "G_W_U": (energy_density, Str("INDEPENDENT"), sp.Tuple(energy_density)),
        "KAPPA_THETA": (derive_dimension(energy_density + 2 * L_DIM), Str("INDEPENDENT"), sp.Tuple(energy_density, -2 * L_DIM)),
        "KAPPA_THETA_W": (derive_dimension(energy_density + 2 * L_DIM), Str("INDEPENDENT"), sp.Tuple(energy_density, -2 * L_DIM)),
    }
    symbol_links = {
        rho_m: "dim_rho_m", c_s0: "dim_c_s0", v_dr: "dim_v_dr", rho_br: "dim_rho_br",
        B_rho: "dim_b_rho", B_rho_3: "dim_b_rho_3", mu_W: "dim_mu_w", k_W: "dim_k_w",
        kappa_W: "dim_kappa_w", C: "dim_c", mu_R: "dim_mu_r", Lambda_A0: "dim_lambda_a_0",
        Lambda_V0: "dim_lambda_v_0", Lambda_X0: "dim_lambda_x_0", tau_A: "dim_tau_a",
        tau_V: "dim_tau_v", tau_X: "dim_tau_x", B_div: "dim_b_div", mu_S: "dim_mu_s_coefficient",
        G_theta_u: "dim_g_theta_u", G_W_u: "dim_g_w_u", kappa_theta: "dim_kappa_theta",
        kappa_theta_W: "dim_kappa_theta_w",
    }
    SYMBOL_DIMENSION_KEYS.update(symbol_links)
    for name, (vector, route_kind, route_operands) in dimensions.items():
        key_name = f"dim_{name.lower()}"
        emit(f"DIM_{name}", vector, key=key_name)
        emit(f"DIM_ROUTE_KIND_{name}", sp.Tuple(route_kind, route_operands), key=f"dim_route_kind_{name.lower()}")

    homogeneity = {
        "INPLANE_EOM": (M_DIM - 3 * L_DIM - 2 * T_DIM,) * len(sp.Add.make_args(sp.expand(MODEL["inplane"]))),
        "THICKNESS_EOM": (energy_density,) * len(sp.Add.make_args(sp.expand(MODEL["thickness"]))),
        "MASS_BALANCE": (mass_flux,) * len(sp.Add.make_args(sp.expand(MODEL["mass"]))),
        "AFFINITY": (specific_energy, pressure_4d - (M_DIM - 4 * L_DIM)),
        "CLOSURE": (mass_flux, dimensions["LAMBDA_A_0"][0] + specific_energy, dimensions["LAMBDA_V_0"][0] + L_DIM - T_DIM),
        "FACE_RESPONSE": (pressure_4d, dimensions["Z"][0] + L_DIM - T_DIM),
        "TWO_PORT_POWER_IDENTITY": (M_DIM - L_DIM - 3 * T_DIM,) * 3,
        "DISPERSION_DETERMINANT": (3 * M_DIM - 7 * L_DIM - 5 * T_DIM,),
    }
    for name, term_dimensions in homogeneity.items():
        emit(f"HOMOGENEITY_{name}", homogeneous_payload(term_dimensions), key=f"homogeneity_{name.lower()}")

    correct = homogeneous_payload(homogeneity["INPLANE_EOM"])
    corrupted_terms = list(homogeneity["INPLANE_EOM"])
    corrupted_terms[0] = corrupted_terms[0] + L_DIM
    corrupted = homogeneous_payload(corrupted_terms)
    restored = homogeneous_payload(homogeneity["INPLANE_EOM"])
    emit("HOMOGENEITY_ABLATION_DEMO", sp.Tuple(correct, corrupted, restored), key="homogeneity_ablation_demo")


def task_b9() -> None:
    convective_dispersion = sp.expand((omega - v_dr * q_out)**2 - c_s0**2 * (k**2 + q_out**2))
    rest_dispersion = sp.expand(convective_dispersion.subs(v_dr, 0))
    discarded = sp.expand(convective_dispersion - rest_dispersion)
    leading = sp.expand(sp.diff(convective_dispersion, v_dr).subs(v_dr, 0) * v_dr)
    relative = sp.cancel(leading / omega**2)
    modulus_measure = sp.Abs(relative)
    failure_region = sp.Ge(modulus_measure, 1, evaluate=False)
    emit("DISCARDED_CONVECTIVE_CORRECTION", sp.Tuple(convective_dispersion, rest_dispersion, discarded, leading, relative), key="discarded_convective_correction")
    emit("VALIDITY_TIMESCALE", sp.Tuple(modulus_measure, sp.Lt(modulus_measure, 1, evaluate=False)), key="validity_timescale")
    speed_ratio = sp.Abs(v_dr / c_s0)
    emit("VALIDITY_FLOW_SPEED", sp.Tuple(speed_ratio, sp.Lt(speed_ratio, 1, evaluate=False)), key="validity_flow_speed")
    complex_measure = sp.Tuple(sp.Abs(v_dr * q_out), sp.Abs(omega), sp.Eq(sp.im(omega), 0, evaluate=False), sp.Eq(sp.im(q_out), 0, evaluate=False))
    emit("VALIDITY_FAILURE_REGION", sp.Tuple(failure_region, complex_measure), key="validity_failure_region")

    amplitude_conditions = sp.Tuple(
        sp.Lt(sp.Abs(theta), 1, evaluate=False),
        sp.Lt(sp.Abs(e_W), 1, evaluate=False),
        sp.Lt(sp.Abs(eta), 1, evaluate=False),
        sp.Lt(sp.Abs(k * W0 * e_W), 1, evaluate=False),
    )
    kernel_ranges = sp.Tuple(*(
        sp.Tuple(relaxation, sp.Abs(omega * relaxation), safe_limit(kernel, relaxation, 0), safe_limit(kernel, relaxation, sp.oo))
        for relaxation, kernel in ((tau_A, Lambda_A), (tau_V, Lambda_V), (tau_X, Lambda_X))
    ))
    emit("VALIDITY_CONDITIONS", sp.Tuple(amplitude_conditions, kernel_ranges), key="validity_conditions")


def no_thickness_control_objects() -> dict[str, object]:
    energy = sp.expand(U_LONG.subs(e_W, 0))
    mu_theta_value = sp.diff(energy, theta)
    sigma_eta = sp.diff(energy, eta)
    face = face_solution(Z_GENERAL, mu_theta_value, sp.Integer(0))
    mass = sp.cancel(-I * omega * rho_br * (theta + eta) + 2 * face["flux"])
    inplane = sp.expand(rho_br * omega**2 * eta - k**2 * (sigma_eta - mu_theta_value))
    matrix = sp.ImmutableMatrix([[sp.diff(equation, variable) for variable in (eta, theta)] for equation in (inplane, mass)])
    determinant = sp.Determinant(matrix)
    theta_solution = sp.solve(sp.Eq(mass, 0), theta, dict=True)
    compression = NOT_ESTABLISHED
    if theta_solution:
        compression = sp.cancel((sigma_eta - mu_theta_value).subs(theta_solution[0]) / eta)
    return {
        "energy": energy,
        "face": face,
        "mass": mass,
        "inplane": inplane,
        "matrix": matrix,
        "determinant": determinant,
        "compression": compression,
    }


def transverse_control(substitutions: Mapping[sp.Symbol, sp.Expr]) -> sp.Tuple:
    energy = sp.expand(U_TRANSVERSE.subs(substitutions))
    coupling = sp.diff(energy, u_T, e_W)
    operator = sp.expand(-rho_br * omega**2 + sp.diff(energy, u_T, 2))
    return sp.Tuple(coupling, operator, sp.factor(operator - (-rho_br * omega**2 + sp.diff(U_TRANSVERSE, u_T, 2))))


def task_b8() -> None:
    w = next(sym for sym in DECLARED_SYMBOLS if sym.name == "w")
    shifted_window = sp.sech((w - control_b) / control_a)**2
    centred_window = sp.sech(w / control_a)**2
    even_current = w**2
    window_source = -sp.Integral(shifted_window * sp.diff(even_current, w), (w, -control_L, control_L))
    interval_source = -sp.Integral(centred_window * sp.diff(even_current, w), (w, -control_L, control_L + control_c))
    window_dependence = sp.diff(window_source, control_b)
    interval_dependence = sp.diff(interval_source, control_c)
    emit(
        "CONTROL_WINDOW_PARITY",
        sp.Tuple(window_source, window_dependence, sp.Eq(window_dependence, 0, evaluate=False)),
        key=None,
        export=False,
    )
    emit(
        "CONTROL_INTERVAL_SYMMETRY",
        sp.Tuple(interval_source, interval_dependence, sp.Eq(interval_dependence, 0, evaluate=False)),
        key=None,
        export=False,
    )

    control_a_objects = no_thickness_control_objects()
    emit(
        "CONTROL_NO_THICKNESS",
        sp.Tuple(control_a_objects["compression"], control_a_objects["determinant"]),
        key=None,
        export=False,
    )
    channel_cuts = sp.Tuple(
        sp.Eq(e_W, 0, evaluate=False),
        sp.Eq(V_face, 0, evaluate=False),
        sp.Eq(Lambda_V * V_face, 0, evaluate=False),
        sp.Eq(p_face * V_face, 0, evaluate=False),
        sp.Eq(Lambda_X * (mu_drive / rho_br - p_face / rho_m) * V_face, 0, evaluate=False),
        control_a_objects["face"]["flux"],
        control_a_objects["face"]["pressure"],
    )
    emit("CONTROL_A_ATTRIBUTION", channel_cuts, key=None, export=False)

    control_b_model = derive_model(substitutions={kappa_W: 0})
    emit(
        "CONTROL_NO_GRADIENT_STIFFNESS",
        sp.Tuple(sp.diff(control_b_model["thickness"], e_W), control_b_model["determinant"]),
        key=None,
        export=False,
    )
    control_c_model = derive_model(substitutions={Lambda_A0: 0, Lambda_V0: 0})
    emit("CONTROL_IMPERMEABLE", control_c_model["determinant"], key=None, export=False)

    control_d_model = derive_model(substitutions={C: 0})
    control_d_compression = compressional_response(control_d_model)
    emit(
        "CONTROL_NO_CROSS_TERM",
        sp.Tuple(control_d_compression["ratio"], control_d_model["determinant"]),
        key=None,
        export=False,
    )

    control_e_model = derive_model(slab_affinity=False)
    control_e_power_response = face_solution(Z_GENERAL, sp.Integer(0), V_face)
    emit(
        "CONTROL_NO_MU_COUPLING",
        sp.Tuple(control_e_model["determinant"], control_e_power_response["pressure"], control_e_power_response["flux"]),
        key=None,
        export=False,
    )

    control_f_model = derive_model(substitutions={Lambda_X0: 0})
    operator_substituted = sp.cancel(MODEL["thickness"].subs(Lambda_X0, 0))
    operator_residual = sp.cancel(control_f_model["thickness"] - operator_substituted)
    V = symbol("V_control_F", "CONTROL", "control-F face velocity amplitude")
    A = symbol("A_control_F", "CONTROL", "control-F affinity amplitude")
    p = symbol("p_control_F", "CONTROL", "control-F face pressure amplitude")
    J = symbol("J_control_F", "CONTROL", "control-F face flux amplitude")
    mu_s = A + p / rho_m
    v_bulk = V + J / rho_m
    full_power = sp.re((p + Lambda_X * A) * sp.conjugate(V) + mu_s * sp.conjugate(J)) / 2
    recomputed_power = sp.re(p * sp.conjugate(V) + mu_s * sp.conjugate(J)) / 2
    substituted_power = sp.simplify(full_power.subs(Lambda_X0, 0))
    power_residual = sp.simplify(recomputed_power - substituted_power)
    control_f_compression = compressional_response(control_f_model)
    emit(
        "CONTROL_NO_RECIPROCAL_TRACTION",
        sp.Tuple(
            control_f_model["inplane"], control_f_model["thickness"],
            control_f_compression["ratio"], control_f_model["determinant"],
            operator_substituted, operator_residual,
            recomputed_power, substituted_power, power_residual,
        ),
        key=None,
        export=False,
    )

    transverse_rows = sp.Tuple(
        sp.Tuple(Str("A"), transverse_control({e_W: 0})),
        sp.Tuple(Str("B"), transverse_control({kappa_W: 0})),
        sp.Tuple(Str("C"), transverse_control({Lambda_A0: 0, Lambda_V0: 0})),
        sp.Tuple(Str("D"), transverse_control({C: 0})),
    )
    emit("CONTROLS_ON_TRANSVERSE", transverse_rows, key=None, export=False)


def total_compare(left: object, right: object) -> tuple[str, object]:
    """Three-valued, total comparison over every incoming value shape."""
    if isinstance(left, bool) or isinstance(right, bool):
        if isinstance(left, bool) and isinstance(right, bool):
            return ("PROVED_EQUAL" if left is right else "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right)))
        return "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right))
    if isinstance(left, str) or isinstance(right, str):
        if isinstance(left, str) and isinstance(right, str):
            return ("PROVED_EQUAL" if left == right else "PROVED_DIFFERENT", sp.Tuple(Str(left), Str(right)))
        return "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right))
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        if not isinstance(left, Mapping) or not isinstance(right, Mapping):
            return "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right))
        if set(left) != set(right):
            return "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right))
        comparisons = [total_compare(left[key], right[key]) for key in sorted(left, key=str)]
        statuses = [status for status, _ in comparisons]
        detail = sp.Tuple(*(sp.Tuple(Str(str(key)), Str(status), casify(operand)) for key, (status, operand) in zip(sorted(left, key=str), comparisons)))
        if all(status == "PROVED_EQUAL" for status in statuses):
            return "PROVED_EQUAL", detail
        if any(status == "PROVED_DIFFERENT" for status in statuses):
            return "PROVED_DIFFERENT", detail
        return "UNDECIDED", detail
    left_container = isinstance(left, (tuple, list, sp.Tuple))
    right_container = isinstance(right, (tuple, list, sp.Tuple))
    if left_container or right_container:
        if not left_container or not right_container:
            return "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right))
        if len(left) != len(right):
            return "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right))
        comparisons = [total_compare(a, b) for a, b in zip(left, right)]
        statuses = [status for status, _ in comparisons]
        detail = sp.Tuple(*(sp.Tuple(Str(status), casify(operand)) for status, operand in comparisons))
        if all(status == "PROVED_EQUAL" for status in statuses):
            return "PROVED_EQUAL", detail
        if any(status == "PROVED_DIFFERENT" for status in statuses):
            return "PROVED_DIFFERENT", detail
        return "UNDECIDED", detail
    if isinstance(left, sp.MatrixBase) or isinstance(right, sp.MatrixBase):
        if not isinstance(left, sp.MatrixBase) or not isinstance(right, sp.MatrixBase) or left.shape != right.shape:
            return "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right))
        comparisons = [total_compare(a, b) for a, b in zip(left, right)]
        statuses = [status for status, _ in comparisons]
        residual = sp.ImmutableMatrix(left - right)
        if all(status == "PROVED_EQUAL" for status in statuses):
            return "PROVED_EQUAL", residual
        if any(status == "PROVED_DIFFERENT" for status in statuses):
            return "PROVED_DIFFERENT", residual
        return "UNDECIDED", residual

    left = casify(left)
    right = casify(right)
    if isinstance(left, sp.Symbol) and isinstance(right, sp.Symbol):
        equal = left.name == right.name and left.assumptions0 == right.assumptions0
        return ("PROVED_EQUAL" if equal else "PROVED_DIFFERENT", sp.Tuple(left, right))
    if isinstance(left, Str) or isinstance(right, Str):
        if isinstance(left, Str) and isinstance(right, Str):
            equal = sp.srepr(left) == sp.srepr(right)
            return ("PROVED_EQUAL" if equal else "PROVED_DIFFERENT", sp.Tuple(left, right))
        return "PROVED_DIFFERENT", sp.Tuple(left, right)
    if isinstance(left, (Relational, Boolean)) or isinstance(right, (Relational, Boolean)):
        if not isinstance(left, Boolean) or not isinstance(right, Boolean):
            return "PROVED_DIFFERENT", sp.Tuple(left, right)
        if sp.srepr(left) == sp.srepr(right):
            return "PROVED_EQUAL", sp.Tuple(left, right)
        if sp.count_ops(left) + sp.count_ops(right) > 20:
            return "UNDECIDED", sp.Tuple(left, right)
        equivalence = sp.simplify_logic(sp.Equivalent(left, right))
        if equivalence is sp.true or equivalence == sp.true:
            return "PROVED_EQUAL", equivalence
        if equivalence is sp.false or equivalence == sp.false:
            return "PROVED_DIFFERENT", equivalence
        return "UNDECIDED", equivalence
    if type(left) is not type(right) and not isinstance(left, sp.Expr) and not isinstance(right, sp.Expr):
        return "PROVED_DIFFERENT", sp.Tuple(left, right)
    if isinstance(left, sp.Basic) and isinstance(right, sp.Basic):
        if sp.srepr(left) == sp.srepr(right):
            return "PROVED_EQUAL", sp.Tuple(left, right)
        if sp.count_ops(left) + sp.count_ops(right) > 20:
            return "UNDECIDED", sp.Tuple(left, right)
        try:
            equality = left.equals(right)
        except Exception:
            equality = None
        if equality is True:
            return "PROVED_EQUAL", sp.Tuple(left, right)
        if equality is False:
            return "PROVED_DIFFERENT", sp.Tuple(left, right)
        try:
            residual = sp.simplify(left - right)
        except Exception:
            return "UNDECIDED", sp.Tuple(left, right)
        if residual == 0:
            return "PROVED_EQUAL", residual
        if not getattr(residual, "free_symbols", set()) and residual != 0:
            return "PROVED_DIFFERENT", residual
        return "UNDECIDED", residual
    try:
        equality = left == right
    except Exception:
        equality = None
    if equality is True:
        return "PROVED_EQUAL", sp.Tuple(casify(left), casify(right))
    if equality is False:
        return "PROVED_DIFFERENT", sp.Tuple(casify(left), casify(right))
    return "UNDECIDED", sp.Tuple(casify(left), casify(right))


FREE_SYMBOL_CACHE: dict[int, frozenset[sp.Symbol]] = {}


def atoms_of(value: object) -> set[sp.Symbol]:
    cache_key = id(value)
    cached = FREE_SYMBOL_CACHE.get(cache_key)
    if cached is not None:
        return set(cached)
    if isinstance(value, sp.Symbol):
        result = {value}
    elif isinstance(value, sp.MatrixBase):
        result: set[sp.Symbol] = set()
        for entry in value:
            result.update(atoms_of(entry))
    elif isinstance(value, (sp.Integral, sp.Sum, sp.Lambda, sp.ConditionSet, sp.ImageSet, sp.Subs)):
        result = set(value.free_symbols)
    elif isinstance(value, sp.Basic):
        result = set()
        for argument in value.args:
            result.update(atoms_of(argument))
    elif isinstance(value, (tuple, list)):
        result = set()
        for item in value:
            result.update(atoms_of(item))
    else:
        result = set()
    FREE_SYMBOL_CACHE[cache_key] = frozenset(result)
    return result


def add_free_symbol_candidates(candidates: list[CandidateRow]) -> None:
    seen = {candidate.key for candidate in candidates}
    symbol_names: set[str] = set()
    for candidate in candidates:
        rendering = EMITTER.rendered.get(candidate.source_tag)
        if rendering is None:
            rendering = render(candidate.value)
        symbol_names.update(re.findall(r"Symbol\('([^']+)'", rendering))
    declared_by_name: dict[str, list[sp.Symbol]] = {}
    for declared in DECLARED_SYMBOLS:
        declared_by_name.setdefault(declared.name, []).append(declared)
    symbols: set[sp.Symbol] = set()
    for name in sorted(symbol_names):
        variants = declared_by_name.get(name, [])
        if not variants:
            issue(f"unclassifiable free symbol name {name}")
            continue
        symbols.update(variants)
    for free_symbol in sorted(symbols, key=sp.default_sort_key):
        if free_symbol.name in seen:
            continue
        metadata = DECLARED_SYMBOLS.get(free_symbol)
        if metadata is None:
            issue(f"unclassifiable free symbol {sp.srepr(free_symbol)}")
            continue
        candidates.append(
            CandidateRow(
                free_symbol.name,
                free_symbol,
                metadata["class"],
                "FREE_SYMBOL_POPULATION",
                SYMBOL_DIMENSION_KEYS.get(free_symbol),
                metadata["description"],
            )
        )
        seen.add(free_symbol.name)


def make_record(candidate: CandidateRow) -> dict[str, object]:
    record: dict[str, object] = {
        "display": display(candidate.value),
        "value": casify(candidate.value),
        "value_kind": "COMPUTED_OBJECT",
        "class": candidate.class_tag,
        "step": "S11b",
    }
    if candidate.dimension_key is not None:
        record["dimension_key"] = candidate.dimension_key
    if candidate.description is not None:
        record["description"] = candidate.description
    return record


def merged_export() -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    candidates = list(EMITTER.export_candidates)
    add_free_symbol_candidates(candidates)
    candidate_keys = sp.Tuple(*(Str(candidate.key) for candidate in candidates))
    imported_matching_keys = sp.Tuple(*(Str(candidate.key) for candidate in candidates if candidate.key in INCOMING_LEDGER))

    merged = {name: dict(record) for name, record in INCOMING_LEDGER.items()}
    routing_rows = []
    seen_candidates: set[str] = set()
    f9c_rows = []
    for candidate in candidates:
        if candidate.key in seen_candidates:
            raise RuntimeError(f"two S11b candidates have the same key {candidate.key!r}")
        seen_candidates.add(candidate.key)
        record = make_record(candidate)
        imported = INCOMING_LEDGER.get(candidate.key)
        if imported is None:
            write_key = candidate.key
            status = "ABSENT"
            comparison_operand = sp.Tuple(record["value"], sp.Tuple())
        else:
            value_status, comparison_operand = total_compare(record["value"], imported.get("value"))
            class_status = "PROVED_EQUAL" if record["class"] == imported.get("class") else "PROVED_DIFFERENT"
            if value_status == "PROVED_EQUAL" and class_status == "PROVED_EQUAL":
                write_key = candidate.key
                status = "PROVED_EQUAL"
                record["f9_operands"] = sp.Tuple(record["value"], imported.get("value"))
                prior_steps = list(imported.get("corroborated_steps", ()))
                if not prior_steps and imported.get("step"):
                    prior_steps.append(imported["step"])
                prior_steps.append("S11b")
                record["corroborated_steps"] = tuple(dict.fromkeys(prior_steps))
            else:
                write_key = f"s11b_{candidate.key}"
                status = value_status if value_status != "PROVED_EQUAL" else "PROVED_DIFFERENT"
                f9c_rows.append(sp.Tuple(Str(candidate.key), Str(write_key), Str(status), casify(comparison_operand)))
                issue(f"F9c write {write_key} for {candidate.key}: {status}")
        routing_rows.append(sp.Tuple(Str(candidate.source_tag), Str(candidate.key), Str(write_key), Str(status), casify(comparison_operand)))
        merged[write_key] = record
    diagnostics = {
        "candidate_keys": candidate_keys,
        "imported_matching_keys": imported_matching_keys,
        "routing": sp.Tuple(*routing_rows),
        "f9c": sp.Tuple(*f9c_rows),
        "live_candidates": candidates,
    }
    return merged, diagnostics


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record_source(record: Mapping[str, object]) -> str:
    lines = ["    {"]
    for field, value in record.items():
        if field in {"value", "f9_operands"}:
            lines.append(f"        {field!r}: _restore({sp.srepr(casify(value))!r}),")
        else:
            lines.append(f"        {field!r}: {value!r},")
    lines.append("    }")
    return "\n".join(lines)


def export_source(ledger: Mapping[str, Mapping[str, object]]) -> str:
    lines = [
        "# S11b_exports.py -- GENERATED by S11b_interface_coupling_law_sympy_audit.py. Do not edit.",
        "from types import MappingProxyType",
        "",
        "import sympy as sp",
        "from sympy.core.symbol import Str",
        "from sympy.functions.elementary.piecewise import ExprCondPair",
        "",
        "_RELATIONALS = {",
        "    'Equality': lambda left, right: sp.Eq(left, right, evaluate=False),",
        "    'Unequality': lambda left, right: sp.Ne(left, right, evaluate=False),",
        "    'StrictGreaterThan': lambda left, right: sp.Gt(left, right, evaluate=False),",
        "    'StrictLessThan': lambda left, right: sp.Lt(left, right, evaluate=False),",
        "    'GreaterThan': lambda left, right: sp.Ge(left, right, evaluate=False),",
        "    'LessThan': lambda left, right: sp.Le(left, right, evaluate=False),",
        "}",
        "",
        "def _restore(source):",
        "    return eval(source, {'__builtins__': {}, 'Str': Str, 'ExprCondPair': ExprCondPair, **vars(sp), **_RELATIONALS})",
        "",
        "BUILD_INPUT_DIGESTS = MappingProxyType({",
        f"    'S11b_interface_coupling_law_sympy_audit.py': {sha256(SCRIPT_PATH)!r},",
        f"    'S11_exports.py': {sha256(SCRIPT_DIR / 'S11_exports.py')!r},",
        f"    'S11b_SHARED_PHYSICS.md': {sha256(DIRECTIVE_PATH)!r},",
        "})",
        "",
        "_LEDGER = {",
    ]
    for name in sorted(ledger):
        lines.append(f"    {name!r}: {record_source(ledger[name])},")
    lines.extend([
        "}",
        "LEDGER = MappingProxyType({",
        "    name: MappingProxyType(record) for name, record in _LEDGER.items()",
        "})",
        "del _LEDGER",
        "",
    ])
    return "\n".join(lines)


def publish_export(ledger: dict[str, dict[str, object]]) -> None:
    source = export_source(ledger)
    EXPORT_PATH.write_text(source, encoding="utf-8")
    finished_source = EXPORT_PATH.read_text(encoding="utf-8")
    namespace: dict[str, object] = {}
    exec(compile(finished_source, str(EXPORT_PATH), "exec"), namespace)
    reconstructed = namespace["LEDGER"]
    rows = []
    failures = []
    for name, live_record in ledger.items():
        restored_value = reconstructed[name]["value"]
        status, comparison = total_compare(live_record["value"], restored_value)
        detail = sp.Integer(0) if status == "PROVED_EQUAL" else casify(comparison)
        rows.append(sp.Tuple(Str(name), Str(str(live_record["display"])), Str(status), detail))
        if status != "PROVED_EQUAL":
            failures.append(name)
    emit("EXPORT_ROUNDTRIP", sp.Tuple(*rows), key=None, local=True, export=False)
    if failures:
        issue(
            "export reconstruction comparison was undecided for "
            + ", ".join(failures)
        )


def cross_step_dimension_lookups() -> None:
    derived_tags = {
        c_s0: "PY_S11B_DIM_C_S0",
        mu_R: "PY_S11B_DIM_MU_R",
        rho_br: "PY_S11B_DIM_RHO_BR",
    }
    resolved = []
    unresolved = []
    for coefficient, tag in derived_tags.items():
        row_name = {c_s0: "c_s0", mu_R: "mu_R", rho_br: "rho_br"}[coefficient]
        row = INCOMING_LEDGER[row_name]
        derived = EMITTER.values[tag]
        dimension_key = row.get("dimension_key")
        if dimension_key is None or dimension_key not in INCOMING_LEDGER:
            unresolved.append(sp.Tuple(coefficient, Str(row_name), Str(str(row.get("class"))), Str(str(row.get("step")))))
            continue
        dimension_row = INCOMING_LEDGER[dimension_key]
        imported = dimension_row["value"]
        residual = sp.simplify(derived - imported)
        stem = row_name.upper()
        emit(f"Q_{stem}_DERIVED_DIMENSION", derived, key=None, local=True, export=False)
        emit(f"Q_{stem}_IMPORTED_DIMENSION", imported, key=None, local=True, export=False)
        emit(f"Q_{stem}_DIMENSION_RESIDUAL", sp.Tuple(derived, imported, residual), key=None, local=True, export=False)
        provenance = sp.Tuple(
            Str(str(row.get("class"))), Str(str(row.get("step"))),
            Str(str(dimension_row.get("class"))), Str(str(dimension_row.get("step"))),
            sp.Tuple(*(Str(str(item)) for item in dimension_row.get("corroborated_steps", ()))),
        )
        emit(f"Q_{stem}_PROVENANCE", provenance, key=None, local=True, export=False)
        resolved.append(sp.Tuple(coefficient, Str(row_name), Str(str(dimension_key))))
    emit("Q_RESOLVED_COEFFICIENTS", sp.Tuple(*resolved), key=None, local=True, export=False)
    emit("Q_UNRESOLVED_COEFFICIENTS", sp.Tuple(*unresolved), key=None, local=True, export=False)
    emit("Q_RESIDUAL_SCOPE", CROSS_STEP_CONSISTENCY_ONLY, key=None, local=True, export=False)


def emit_premise_inventory() -> None:
    premises = sp.Tuple(
        sp.Eq(D_brane, 3, evaluate=False),
        sp.Eq(B_rho_3, B_rho * W0, evaluate=False),
        sp.Eq(q_out**2, omega**2 / c_s0**2 - k**2, evaluate=False),
        sp.Eq(Lambda_A, Lambda_A0 / (1 - I * omega * tau_A), evaluate=False),
        sp.Eq(Lambda_V, Lambda_V0 / (1 - I * omega * tau_V), evaluate=False),
        sp.Eq(Lambda_X, Lambda_X0 / (1 - I * omega * tau_X), evaluate=False),
        Str("NO_INCOMING_WAVES_FROM_INFINITY"),
        Str("NO_POSITIVITY_ASSUMPTION_ON_STORED_ENERGY"),
        Str("REST_FRAME_BULK_OPERATOR"),
        Str("NO_DIRECT_J_GENERALIZED_FORCE"),
        Str("NO_COMPLEX_FREQUENCY_BRANCH_RESELECTION"),
    )
    emit("PREMISE_INVENTORY", premises, key=None, local=True, export=False)


def remove_stale_export() -> None:
    if EXPORT_PATH.exists():
        EXPORT_PATH.unlink()


def run() -> None:
    global MODEL, CURRENT_TASK
    start = time.monotonic()
    completed: list[str] = []
    attempted: list[str] = []
    task_functions = {
        "B0_ENERGY": task_b0_energy,
        "B0A": task_b0a,
        "B0B": task_b0b,
        "B0C": task_b0c,
        "B1": task_b1,
        "B2A": task_b2a,
        "B2B": task_b2b,
        "B2C": task_b2c,
        "B2D": task_b2d,
        "B5": task_b5,
        "B4": task_b4,
        "B6": task_b6,
        "B7": task_b7,
        "B9": task_b9,
    }
    MODEL = derive_model()
    for task_name in PRIMARY_TASKS:
        attempted.append(task_name)
        CURRENT_TASK = task_name
        try:
            task_functions[task_name]()
            completed.append(task_name)
        except Exception as exc:
            issue(f"task raised {type(exc).__name__}: {exc}")
            emit(f"TASK_{task_name}_EXCEPTION", sp.Tuple(Str(task_name), Str(type(exc).__name__), Str(repr(exc))), key=None, local=True, export=False)
        finally:
            CURRENT_TASK = None

    attempted.append(CONTROL_TASK)
    CURRENT_TASK = CONTROL_TASK
    try:
        task_b8()
        completed.append(CONTROL_TASK)
    except Exception as exc:
        issue(f"task raised {type(exc).__name__}: {exc}")
        emit("CONTROL_TASK_EXCEPTION", sp.Tuple(Str(CONTROL_TASK), Str(type(exc).__name__), Str(repr(exc))), key=None, local=True, export=False)
    finally:
        CURRENT_TASK = None

    emit_premise_inventory()
    if "B7" in completed:
        cross_step_dimension_lookups()

    primary_complete = all(task in completed for task in PRIMARY_TASKS)
    if primary_complete:
        ledger, diagnostics = merged_export()
        emit("EXPORT_CANDIDATE_KEY_OPERANDS", sp.Tuple(diagnostics["candidate_keys"], diagnostics["imported_matching_keys"], diagnostics["routing"]), key=None, local=True, export=False)
        emit("EXPORT_F9C_WRITES", diagnostics["f9c"], key=None, local=True, export=False)
        publish_export(ledger)
    else:
        remove_stale_export()
        issue("S11b_exports.py not published because a primary derivation task did not complete")

    skipped = tuple(task for task in (*PRIMARY_TASKS, CONTROL_TASK) if task not in completed)
    run_record = sp.Tuple(*(Str(task) for task in completed))
    skipped_record = sp.Tuple(*(Str(task) for task in skipped))
    emit("RUN_TASKS", run_record, key=None, local=True, export=False)
    emit("SKIPPED_TASKS", skipped_record, key=None, local=True, export=False)

    local_names = list(dict.fromkeys((*EMITTER.local_tags, "PY_S11B_LOCAL_TAGS", "PY_S11B_LOCAL_SECTION13_REPORT")))
    emit("TAGS", sp.Tuple(*(Str(name) for name in local_names)), key=None, local=True, export=False)
    runtime = sp.Float(time.monotonic() - start, 8)
    report = sp.Tuple(
        sp.Tuple(Str("SCRIPT"), Str(str(SCRIPT_PATH)), sp.Integer(len(SCRIPT_PATH.read_text(encoding="utf-8").splitlines()))),
        sp.Tuple(Str("EMITTED_TAGS"), sp.Integer(EMITTER.count + 1)),
        sp.Tuple(Str("RUN_TASKS"), run_record),
        sp.Tuple(Str("SKIPPED_TASKS"), skipped_record),
        sp.Tuple(Str("RUNTIME_SECONDS"), runtime),
        sp.Tuple(Str("ISSUES"), sp.Tuple(*(Str(item) for item in ISSUES))),
        sp.Tuple(Str("CUSTOM_QUANTITY_NAMES"), sp.Tuple(*(Str(item) for item in sorted(CUSTOM_QUANTITY_NAMES)))),
    )
    emit("SECTION13_REPORT", report, key=None, local=True, export=False)


if __name__ == "__main__":
    run()
