#!/usr/bin/env python3
"""SymPy audit for S11c-a non-uniform interface geometry.

The only physics input for this engine is
``directives/S11c_a_SHARED_PHYSICS.md``.  The script streams CAS payloads and,
when every primary T-0...T-i task completes, publishes the accumulated flat
ledger for the next sub-step.
"""

from __future__ import annotations

from collections.abc import Mapping
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
DIRECTIVE_PATH = SCRIPT_DIR.parent / "directives" / "S11c_a_SHARED_PHYSICS.md"
EXPORT_PATH = SCRIPT_DIR / "S11c_a_exports.py"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from S11b_exports import LEDGER as INCOMING_LEDGER  # noqa: E402


CLASS_TAGS = {"KNOB", "STRUCTURAL", "COORDINATE", "CONTROL", "PREMISE", "DERIVED"}
PRIMARY_TASKS = (
    "T0", "TA", "TA_PRIME", "TA_DOUBLE", "TB", "TC", "TC_PRIME",
    "TD", "TE", "TF", "TG", "TH", "TI",
)
CONTROL_TASKS = ("REP_INVARIANCE", "INDEPENDENCE", "FORM", "UNIFORM", "HOMOGENEITY")
BRANCHES = ("LAB_HELD", "MATERIAL_ADVECTED")
FACES = (1, -1)
DOFS = ("DELTA_W", "ZETA_C")
DENSITY_REPS = ("RHO4_CONSTANT", "RHOBR_CONSTANT")

PROVED_EQUAL = "PROVED_EQUAL"
PROVED_DIFFERENT = "PROVED_DIFFERENT"
UNDECIDED = "UNDECIDED"


DECLARED_SYMBOLS: dict[sp.Symbol, dict[str, str]] = {}
ISSUES: list[str] = []
OPERATIONAL_FAILURES: list[str] = []
CURRENT_TASK: str | None = None
RAW_PRIMARY: dict[str, object] = {}
DIMENSION_PRIMARY: dict[str, object] = {}


def register_symbol(value: sp.Symbol, class_tag: str, description: str) -> sp.Symbol:
    if class_tag not in CLASS_TAGS:
        raise ValueError(f"unknown class {class_tag!r}")
    if not isinstance(value, sp.Symbol):
        raise TypeError(f"declared object is not a Symbol: {value!r}")
    DECLARED_SYMBOLS[value] = {"class": class_tag, "description": description}
    return value


def symbol(name: str, class_tag: str, description: str, **assumptions: object) -> sp.Symbol:
    return register_symbol(sp.Symbol(name, **assumptions), class_tag, description)


def inherited(name: str, class_tag: str, description: str) -> sp.Symbol:
    """Bind a supplied quantity to the live object carried by S11b."""
    return register_symbol(INCOMING_LEDGER[name]["value"], class_tag, description)


# Imported identities named by the build directive and supplied S11b laws.
c_s0 = inherited("c_s0", "KNOB", "bulk sound speed")
mu_R = inherited("mu_R", "KNOB", "uniform curl-type brane modulus")
rho_br = inherited("rho_br", "KNOB", "uniform slab-integrated inertia density")
W0 = inherited("W_0", "KNOB", "uniform reference thickness")
e_W = inherited("e_W", "COORDINATE", "uniform-reference fractional thickness field")
rho_m = inherited("rho_m", "KNOB", "bulk mass density")
v_dr = inherited("v_bulk_normal_0", "KNOB", "steady bulk-normal drain scope-limit speed")
omega = inherited("omega", "COORDINATE", "frequency operand of inherited response kernels")
Lambda_A0 = inherited("Lambda_A_0", "KNOB", "zero-frequency affinity mobility")
Lambda_V0 = inherited("Lambda_V_0", "KNOB", "zero-frequency face-velocity mobility")
Lambda_X0 = inherited("Lambda_X_0", "KNOB", "zero-frequency reciprocal traction coefficient")
tau_A = inherited("tau_A", "KNOB", "affinity relaxation time")
tau_V = inherited("tau_V", "KNOB", "velocity relaxation time")
tau_X = inherited("tau_X", "KNOB", "traction relaxation time")
theta = inherited("theta", "COORDINATE", "Eulerian fractional densification")
zeta_c = inherited("zeta_c", "COORDINATE", "face-centre displacement")
delta_v_theta = inherited("delta_v_theta", "COORDINATE", "virtual densification variation")
delta_v_e_W = inherited("delta_v_e_W", "COORDINATE", "virtual fractional-thickness variation")
delta_p_plus = inherited("delta_p_plus", "COORDINATE", "upper perturbation face pressure")
delta_p_minus = inherited("delta_p_minus", "COORDINATE", "lower perturbation face pressure")
x1 = inherited("x_1", "COORDINATE", "first Eulerian in-plane coordinate")
x2 = inherited("x_2", "COORDINATE", "second Eulerian in-plane coordinate")
x3 = inherited("x_3", "COORDINATE", "third Eulerian in-plane coordinate")
t = inherited("t", "COORDINATE", "time coordinate")
w = inherited("w", "COORDINATE", "bulk normal coordinate")
x = (x1, x2, x3)


# Fresh S11c-a model-definition objects and profile jets.
D_brane = symbol("D_brane", "STRUCTURAL", "in-plane spatial dimension", integer=True, positive=True)
epsilon = symbol("epsilon_shape", "CONTROL", "wave-amplitude shape bookkeeper")
eta_bg = symbol("eta_bg", "KNOB", "zero-jet background contrast bookkeeper", real=True)
sigma_W = symbol("sigma_W", "DERIVED", "first-jet background bookkeeper", real=True)
L_W = symbol("L_W", "KNOB", "independent background profile length", positive=True)
w1_profile = symbol("w1_profile", "KNOB", "dimensionless thickness profile value", real=True)
m1_profile = symbol("m1_profile", "KNOB", "dimensionless modulus profile value", real=True)
w1_grad = tuple(symbol(f"w1_profile_d{i}", "KNOB", f"dimensionless thickness-profile first jet {i}", real=True) for i in range(1, 4))
m1_grad = tuple(symbol(f"m1_profile_d{i}", "KNOB", f"dimensionless modulus-profile first jet {i}", real=True) for i in range(1, 4))
w1_hess = tuple(tuple(symbol(f"w1_profile_d{i}d{j}", "KNOB", f"dimensionless thickness-profile second jet {i},{j}", real=True) for j in range(1, 4)) for i in range(1, 4))

W_bg = symbol("W_bg", "DERIVED", "spatially varying background thickness")
mu_R_bg = symbol("mu_R_bg", "DERIVED", "spatially varying twist modulus")
rho4_rho4 = symbol("rho_4D_bg_rho4_constant", "DERIVED", "RHO4-CONSTANT background four-density")
rhobr_rho4 = symbol("rho_br_bg_rho4_constant", "DERIVED", "RHO4-CONSTANT integrated background density")
rho4_rhobr = symbol("rho_4D_bg_rhobr_constant", "DERIVED", "RHOBR-CONSTANT background four-density")
rhobr_rhobr = symbol("rho_br_bg_rhobr_constant", "DERIVED", "RHOBR-CONSTANT integrated background density")
e_W_bg = symbol("e_W_bg", "DERIVED", "local-background fractional thickness field")

grad_W = tuple(symbol(f"W_bg_d{i}", "DERIVED", f"background thickness first jet {i}") for i in range(1, 4))
hess_W = tuple(tuple(symbol(f"W_bg_d{i}d{j}", "DERIVED", f"background thickness second jet {i},{j}") for j in range(1, 4)) for i in range(1, 4))
grad_mu_R_bg = tuple(symbol(f"mu_R_bg_d{i}", "DERIVED", f"background modulus first jet {i}") for i in range(1, 4))

# Wave and virtual jets.  They remain distinct operands; in particular, no
# physical u is substituted for delta_v_u in T-g.
u = tuple(symbol(f"u_{i}", "COORDINATE", f"in-plane displacement component {i}") for i in range(1, 4))
u_t = tuple(symbol(f"u_{i}_t", "COORDINATE", f"time derivative of displacement component {i}") for i in range(1, 4))
grad_u = tuple(tuple(symbol(f"u_{i}_d{j}", "COORDINATE", f"displacement jet {i},{j}") for j in range(1, 4)) for i in range(1, 4))
grad_u_t = tuple(tuple(symbol(f"u_{i}_t_d{j}", "COORDINATE", f"velocity jet {i},{j}") for j in range(1, 4)) for i in range(1, 4))
theta_t = symbol("theta_t", "COORDINATE", "time derivative of densification")
grad_theta = tuple(symbol(f"theta_d{i}", "COORDINATE", f"densification first jet {i}") for i in range(1, 4))
e_W_t = symbol("e_W_t", "COORDINATE", "time derivative of uniform-reference fractional thickness")
grad_e_W = tuple(symbol(f"e_W_d{i}", "COORDINATE", f"fractional-thickness first jet {i}") for i in range(1, 4))
zeta_c_t = symbol("zeta_c_t", "COORDINATE", "time derivative of centre displacement")
grad_zeta_c = tuple(symbol(f"zeta_c_d{i}", "COORDINATE", f"centre-displacement first jet {i}") for i in range(1, 4))

delta_v_u = tuple(symbol(f"delta_v_u_{i}", "COORDINATE", f"virtual displacement component {i}") for i in range(1, 4))
grad_delta_v_u = tuple(tuple(symbol(f"delta_v_u_{i}_d{j}", "COORDINATE", f"virtual-displacement jet {i},{j}") for j in range(1, 4)) for i in range(1, 4))
delta_v_zeta_c = symbol("delta_v_zeta_c", "COORDINATE", "virtual centre-face displacement")

# Traced bulk perturbations and background normal derivatives used by the
# supplied shifted-trace law.
delta_v_bulk = {
    s: tuple(symbol(f"delta_v_bulk_{'plus' if s == 1 else 'minus'}_{i}", "COORDINATE", f"bulk velocity perturbation at face {s}, component {i}") for i in range(1, 5))
    for s in FACES
}
dw_v_bulk_0 = {
    s: tuple(symbol(f"d_w_v_bulk_0_{'plus' if s == 1 else 'minus'}_{i}", "PREMISE", f"background bulk-velocity normal jet at face {s}, component {i}") for i in range(1, 5))
    for s in FACES
}
dw_delta_p_0 = {
    1: symbol("d_w_delta_p_0_plus", "PREMISE", "upper background pressure normal jet"),
    -1: symbol("d_w_delta_p_0_minus", "PREMISE", "lower background pressure normal jet"),
}
rho4_bulk_1 = symbol("delta_rho_4D_bulk", "COORDINATE", "bulk density perturbation in the projection")
rho4_bulk_1_t = symbol("delta_rho_4D_bulk_t", "COORDINATE", "time derivative of bulk density perturbation")
j_bulk = tuple(symbol(f"delta_j_bulk_{i}", "COORDINATE", f"bulk current perturbation component {i}") for i in range(1, 5))
grad_j_bulk = tuple(tuple(symbol(f"delta_j_bulk_{i}_d{j}", "COORDINATE", f"bulk-current in-plane jet {i},{j}") for j in range(1, 4)) for i in range(1, 5))

mu_theta_branch = {
    "LAB_HELD": symbol("mu_theta_L", "PREMISE", "reserved S11c-b LAB_HELD variable-coefficient operand"),
    "MATERIAL_ADVECTED": symbol("mu_theta_M", "PREMISE", "reserved S11c-b MATERIAL_ADVECTED branch-anchored operand"),
}

support_body = symbol("f_hold_0", "PREMISE", "declared background body support")
support_face_plus = symbol("t_hold_plus_0", "PREMISE", "declared upper background traction support")
support_face_minus = symbol("t_hold_minus_0", "PREMISE", "declared lower background traction support")
theta_0 = symbol("theta_0", "PREMISE", "declared background densification")


I = sp.I
Lambda_A = sp.cancel(Lambda_A0 / (1 - I * omega * tau_A))
Lambda_V = sp.cancel(Lambda_V0 / (1 - I * omega * tau_V))
Lambda_X = sp.cancel(Lambda_X0 / (1 - I * omega * tau_X))


# Dimension vectors are ordered [L,T,M].  They are operands of a small CAS
# algebra rather than annotations of results.
DIM_ZERO = sp.ImmutableMatrix([0, 0, 0])
DIM_L = sp.ImmutableMatrix([1, 0, 0])
DIM_T = sp.ImmutableMatrix([0, 1, 0])
DIM_M = sp.ImmutableMatrix([0, 0, 1])


def dim_mul(*dimensions: sp.MatrixBase) -> sp.ImmutableMatrix:
    total = sp.zeros(3, 1)
    for dimension in dimensions:
        total += dimension
    return sp.ImmutableMatrix(total)


def dim_div(left: sp.MatrixBase, right: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(left - right)


DIM_VELOCITY = dim_div(DIM_L, DIM_T)
DIM_RHO4 = dim_mul(DIM_M, -4 * DIM_L)
DIM_RHOBR = dim_mul(DIM_M, -3 * DIM_L)
DIM_PRESSURE = dim_mul(DIM_RHO4, 2 * DIM_VELOCITY)
DIM_AFFINITY = dim_div(DIM_PRESSURE, DIM_RHO4)
DIM_FLUX = dim_mul(DIM_RHO4, DIM_VELOCITY)
DIM_WORK = dim_mul(DIM_PRESSURE, DIM_L)
DIM_MASS_RATE = dim_div(DIM_RHOBR, DIM_T)

SYMBOL_DIMENSIONS: dict[sp.Symbol, sp.ImmutableMatrix] = {
    W0: DIM_L,
    L_W: DIM_L,
    rho_m: DIM_RHO4,
    rho_br: DIM_RHOBR,
    delta_p_plus: DIM_PRESSURE,
    delta_p_minus: DIM_PRESSURE,
    mu_theta_branch["LAB_HELD"]: dim_mul(DIM_RHOBR, DIM_AFFINITY),
    mu_theta_branch["MATERIAL_ADVECTED"]: dim_mul(DIM_RHOBR, DIM_AFFINITY),
    Lambda_A0: dim_div(DIM_FLUX, DIM_AFFINITY),
    Lambda_V0: dim_div(DIM_FLUX, DIM_VELOCITY),
    Lambda_X0: dim_div(DIM_PRESSURE, DIM_AFFINITY),
    tau_A: DIM_T,
    rho4_bulk_1_t: dim_div(DIM_RHO4, DIM_T),
    delta_v_zeta_c: DIM_L,
    theta: DIM_ZERO,
    theta_t: -DIM_T,
    e_W: DIM_ZERO,
    e_W_t: -DIM_T,
    delta_v_theta: DIM_ZERO,
    delta_v_e_W: DIM_ZERO,
}
SYMBOL_DIMENSIONS.update({item: DIM_VELOCITY for item in u_t})
SYMBOL_DIMENSIONS.update({item: DIM_VELOCITY for values in delta_v_bulk.values() for item in values})
SYMBOL_DIMENSIONS.update({item: DIM_FLUX for item in j_bulk})
SYMBOL_DIMENSIONS.update({item: dim_div(DIM_FLUX, DIM_L) for row in grad_j_bulk for item in row})
SYMBOL_DIMENSIONS.update({item: -DIM_T for row in grad_u_t for item in row})


def dimension_of(expression: sp.Expr) -> sp.ImmutableMatrix:
    if expression.is_Number:
        return DIM_ZERO
    if isinstance(expression, sp.Symbol):
        if expression not in SYMBOL_DIMENSIONS:
            raise KeyError(f"dimension unavailable for {expression}")
        return SYMBOL_DIMENSIONS[expression]
    if isinstance(expression, sp.Mul):
        return sp.ImmutableMatrix(sum((dimension_of(argument) for argument in expression.args), sp.zeros(3, 1)))
    if isinstance(expression, sp.Pow) and expression.exp.is_Number:
        return sp.ImmutableMatrix(expression.exp * dimension_of(expression.base))
    if isinstance(expression, sp.Add):
        dimensions = tuple(dimension_of(argument) for argument in expression.args)
        first = dimensions[0]
        if any(item != first for item in dimensions[1:]):
            return sp.ImmutableMatrix([sp.nan, sp.nan, sp.nan])
        return first
    raise TypeError(f"dimension trace has no rule for {sp.srepr(expression)}")


def casify(value: object) -> object:
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


@dataclass(frozen=True)
class CandidateRow:
    key: str
    value: object
    class_tag: str
    source_tag: str
    description: str | None = None


class Emitter:
    def __init__(self) -> None:
        self.count = 0
        self.values: dict[str, object] = {}
        self.rendered: dict[str, str] = {}
        self.local_tags: list[str] = []
        self.primary_tags: list[str] = []
        self.export_candidates: list[CandidateRow] = []

    def emit(
        self,
        tag: str,
        payload: object,
        *,
        export_key: str | None = None,
        class_tag: str = "DERIVED",
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
            self.export_candidates.append(CandidateRow(export_key, payload, class_tag, tag, description))
        return payload


EMITTER = Emitter()


def emit(
    quantity: str,
    payload: object,
    *,
    key: str | None = None,
    local: bool = False,
    export: bool = False,
    class_tag: str = "DERIVED",
) -> object:
    infix = "LOCAL_" if local else ""
    tag = f"PY_S11CA_{infix}{quantity}"
    export_key = key if export and not local else None
    emitted = EMITTER.emit(tag, payload, export_key=export_key, class_tag=class_tag)
    if export_key is not None:
        EMITTER.primary_tags.append(tag)
    return emitted


def issue(message: str) -> None:
    ISSUES.append(f"{CURRENT_TASK or 'GLOBAL'}: {message}")


def dot(left: tuple[sp.Expr, ...], right: tuple[sp.Expr, ...]) -> sp.Expr:
    return sp.Add(*(a * b for a, b in zip(left, right)))


def trace_grad(matrix: tuple[tuple[sp.Expr, ...], ...]) -> sp.Expr:
    return sp.Add(*(matrix[i][i] for i in range(3)))


def determinant_jacobian(jet: tuple[tuple[sp.Expr, ...], ...], parameter: sp.Symbol) -> sp.Expr:
    matrix = sp.eye(3) + parameter * sp.Matrix(jet)
    return sp.det(matrix)


def profile_substitutions() -> dict[sp.Symbol, sp.Expr]:
    substitutions: dict[sp.Symbol, sp.Expr] = {
        W_bg: W0 * (1 + eta_bg * w1_profile),
        mu_R_bg: mu_R * (1 + eta_bg * m1_profile),
    }
    substitutions.update({grad_W[i]: sigma_W * w1_grad[i] for i in range(3)})
    substitutions.update({hess_W[i][j]: sigma_W * w1_hess[i][j] / L_W for i in range(3) for j in range(3)})
    return substitutions


PROFILE_SUBS = profile_substitutions()


def map_object(value: object, scalar_function) -> object:
    if isinstance(value, Mapping):
        return {key: map_object(item, scalar_function) for key, item in value.items()}
    if isinstance(value, sp.MatrixBase):
        return sp.ImmutableMatrix(value.rows, value.cols, [scalar_function(item) for item in value])
    if isinstance(value, Relational):
        return type(value)(
            map_object(value.lhs, scalar_function),
            map_object(value.rhs, scalar_function),
            evaluate=False,
        )
    if isinstance(value, sp.Tuple):
        return sp.Tuple(*(map_object(item, scalar_function) for item in value))
    if isinstance(value, (tuple, list)):
        return tuple(map_object(item, scalar_function) for item in value)
    if isinstance(value, sp.Basic):
        return scalar_function(value)
    return value


_FINAL_CACHE: dict[sp.Basic, sp.Basic] = {}


def final_scalar(expression: sp.Basic) -> sp.Basic:
    if isinstance(expression, (Str, Boolean)):
        return expression
    cached = _FINAL_CACHE.get(expression)
    if cached is not None:
        return cached
    substituted = expression.xreplace(PROFILE_SUBS)
    result = sp.Integer(0)
    for eta_degree in (0, 1):
        for sigma_degree in (0, 1):
            coefficient = sp.diff(substituted, eta_bg, eta_degree, sigma_W, sigma_degree)
            coefficient = coefficient.subs({eta_bg: 0, sigma_W: 0})
            result += coefficient * eta_bg**eta_degree * sigma_W**sigma_degree
    result = sp.factor_terms(result)
    _FINAL_CACHE[expression] = result
    return result


def finalize(value: object) -> object:
    return map_object(value, final_scalar)


def multigrade(value: object) -> sp.Tuple:
    value = casify(value)
    expressions: list[sp.Expr] = []
    if isinstance(value, sp.MatrixBase):
        expressions.extend(value)
    elif isinstance(value, sp.Tuple):
        pass
    elif isinstance(value, Relational):
        expressions.extend((value.lhs, value.rhs))
    elif isinstance(value, sp.Basic):
        expressions.append(value)
    grades: set[tuple[int, int, int]] = set()
    for expression in expressions:
        if isinstance(expression, (Str, Boolean)):
            continue
        for term in sp.Add.make_args(sp.expand(expression)):
            powers = term.as_powers_dict()
            degrees = []
            admissible = True
            for bookkeeper in (epsilon, eta_bg, sigma_W):
                degree = powers.get(bookkeeper, sp.Integer(0))
                if not degree.is_Integer:
                    admissible = False
                    break
                degrees.append(int(degree))
            if admissible:
                grades.add(tuple(degrees))
    if isinstance(value, sp.Tuple):
        for item in value:
            for grade in multigrade(item):
                grades.add(tuple(int(entry) for entry in grade))
    if not grades:
        grades.add((0, 0, 0))
    return sp.Tuple(*(sp.Tuple(*grade) for grade in sorted(grades)))


def case_payload(raw_cases: Mapping[object, object], dimensions: object) -> object:
    rows = {}
    for case, raw in raw_cases.items():
        dimension = dimensions.get(case, dimensions) if isinstance(dimensions, Mapping) else dimensions
        rows[case] = sp.Tuple(
            sp.Tuple(Str("VALUE"), casify(raw)),
            sp.Tuple(Str("MULTIGRADE"), multigrade(raw)),
            sp.Tuple(Str("DIMENSION_L_T_M"), casify(dimension)),
        )
    return rows


def source_jet_scales(direction: int | None = None, reverse_upper_x1: bool = False, face: int | None = None) -> tuple[int, int, int]:
    scales = [1, 1, 1]
    if direction is not None:
        scales[direction - 1] = 0
    if reverse_upper_x1 and face == 1:
        scales[0] = -1
    return tuple(scales)


def dx(expression: sp.Expr, direction: int, scales: tuple[int, int, int] = (1, 1, 1)) -> sp.Expr:
    """Total in-plane derivative on the finite jet algebra."""
    i = direction
    derivative_map: dict[sp.Symbol, sp.Expr] = {
        W_bg: scales[i] * grad_W[i],
        theta: grad_theta[i],
        e_W: grad_e_W[i],
        zeta_c: grad_zeta_c[i],
    }
    for j in range(3):
        derivative_map[grad_W[j]] = scales[j] * hess_W[j][i]
        derivative_map[u[j]] = grad_u[j][i]
        derivative_map[u_t[j]] = grad_u_t[j][i]
    for j in range(4):
        derivative_map[j_bulk[j]] = grad_j_bulk[j][i]
    return sp.Add(*(sp.diff(expression, atom) * derivative for atom, derivative in derivative_map.items()))


def dt(expression: sp.Expr) -> sp.Expr:
    derivative_map: dict[sp.Symbol, sp.Expr] = {
        theta: theta_t,
        e_W: e_W_t,
        zeta_c: zeta_c_t,
        rho4_bulk_1: rho4_bulk_1_t,
    }
    derivative_map.update({u[i]: u_t[i] for i in range(3)})
    derivative_map.update({grad_u[i][j]: grad_u_t[i][j] for i in range(3) for j in range(3)})
    return sp.Add(*(sp.diff(expression, atom) * derivative for atom, derivative in derivative_map.items()))


def shape(expression: object, parameter: sp.Symbol) -> object:
    return map_object(expression, lambda item: sp.diff(item, parameter).subs(parameter, 0))


def density_pair(representative: str, W_value: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    rho4_ref = rho_br / W0
    if representative == "RHO4_CONSTANT":
        return rho4_ref, rho4_ref * W_value
    return rho_br / W_value, rho_br


def branch_profile(expression: sp.Expr, branch: str, parameter: sp.Symbol, scales: tuple[int, int, int]) -> sp.Expr:
    if branch == "LAB_HELD":
        return expression
    displacement_derivative = sp.Add(*(u[i] * dx(expression, i, scales) for i in range(3)))
    return expression - parameter * displacement_derivative


def dof_fields(dof: str, face: int) -> tuple[sp.Expr, sp.Expr, tuple[sp.Expr, ...], sp.Expr, sp.Expr]:
    if dof == "DELTA_W":
        return (
            face * W0 * e_W / 2,
            face * W0 * e_W_t / 2,
            tuple(face * W0 * grad_e_W[i] / 2 for i in range(3)),
            W0 * e_W,
            W0 * delta_v_e_W,
        )
    return zeta_c, zeta_c_t, grad_zeta_c, sp.Integer(0), sp.Integer(0)


@dataclass(frozen=True)
class FaceSource:
    branch: str
    face: int
    dof: str
    representative: str
    route: str
    parameter: sp.Symbol
    h0: sp.Expr
    dh: sp.Expr
    normal_exact: tuple[sp.Expr, ...]
    measure_exact: sp.Expr
    virtual_displacement: tuple[sp.Expr, ...]
    face_velocity_exact: tuple[sp.Expr, ...]
    pressure_trace_first: sp.Expr
    bulk_velocity_trace_first: tuple[sp.Expr, ...]
    rho4_bg_exact: sp.Expr
    rhobr_bg_exact: sp.Expr


_FACE_CACHE: dict[tuple[object, ...], FaceSource] = {}


def build_face_source(
    branch: str,
    face: int,
    dof: str,
    representative: str = "RHO4_CONSTANT",
    *,
    route: str = "EULERIAN",
    ablate_direction: int | None = None,
    reverse_upper_x1: bool = False,
) -> FaceSource:
    cache_key = (branch, face, dof, representative, route, ablate_direction, reverse_upper_x1)
    cached = _FACE_CACHE.get(cache_key)
    if cached is not None:
        return cached
    parameter = sp.Dummy("shape_parameter", real=True)
    scales = source_jet_scales(ablate_direction, reverse_upper_x1, face)
    zeta, zeta_t, _, _, virtual_delta_W = dof_fields(dof, face)
    advective_height = sp.Add(*(u[i] * scales[i] * grad_W[i] for i in range(3)))
    h0 = face * W_bg / 2
    dh = zeta if branch == "LAB_HELD" else zeta - face * advective_height / 2
    h_exact = h0 + parameter * dh

    if route == "EULERIAN":
        grad_h = tuple(dx(h_exact, i, scales) for i in range(3))
        denominator = sp.sqrt(1 + dot(grad_h, grad_h))
        normal_exact = tuple([-face * component / denominator for component in grad_h] + [face / denominator])
        measure_exact = denominator
    else:
        # Tangents of the fixed-X parametric map give an independent normal
        # construction.  Their cofactor is mapped back to Eulerian jets.
        tangent_vertical = tuple(dx(h_exact, i, scales) for i in range(3))
        cofactor = tuple([-face * component for component in tangent_vertical] + [face])
        norm = sp.sqrt(dot(cofactor, cofactor))
        normal_exact = tuple(component / norm for component in cofactor)
        measure_exact = norm

    virtual_zeta = delta_v_zeta_c if dof == "ZETA_C" else face * virtual_delta_W / 2
    if branch == "LAB_HELD":
        virtual_vertical = virtual_zeta + face * dot(delta_v_u, tuple(scales[i] * grad_W[i] for i in range(3))) / 2
        velocity_vertical = zeta_t + face * dot(u_t, tuple(scales[i] * grad_W[i] for i in range(3))) / 2
    else:
        virtual_vertical = virtual_zeta
        velocity_vertical = zeta_t
    virtual_displacement = tuple(delta_v_u) + (virtual_vertical,)
    face_velocity_exact = tuple(parameter * item for item in tuple(u_t) + (velocity_vertical,))

    pressure_face = delta_p_plus if face == 1 else delta_p_minus
    pressure_trace_first = pressure_face + dh * dw_delta_p_0[face]
    bulk_velocity_trace_first = tuple(delta_v_bulk[face][i] + dh * dw_v_bulk_0[face][i] for i in range(4))

    W_alpha = branch_profile(W_bg, branch, parameter, scales)
    rho4_exact, rhobr_exact = density_pair(representative, W_alpha)
    source = FaceSource(
        branch, face, dof, representative, route, parameter, h0, dh,
        normal_exact, measure_exact, virtual_displacement, face_velocity_exact,
        pressure_trace_first, bulk_velocity_trace_first, rho4_exact, rhobr_exact,
    )
    _FACE_CACHE[cache_key] = source
    return source


def face_normal_raw(source: FaceSource) -> sp.Tuple:
    background = tuple(component.subs(source.parameter, 0) for component in source.normal_exact)
    derivative = shape(source.normal_exact, source.parameter)
    return sp.Tuple(sp.ImmutableMatrix(background), epsilon * sp.ImmutableMatrix(derivative), sp.ImmutableMatrix(background) + epsilon * sp.ImmutableMatrix(derivative))


def face_measure_raw(source: FaceSource) -> sp.Tuple:
    background = source.measure_exact.subs(source.parameter, 0)
    derivative = shape(source.measure_exact, source.parameter)
    return sp.Tuple(background, epsilon * derivative, background + epsilon * derivative)


def face_velocity_raw(source: FaceSource) -> sp.Expr:
    return epsilon * shape(dot(source.normal_exact, source.face_velocity_exact), source.parameter)


def relative_flux_raw(source: FaceSource) -> sp.Expr:
    bulk_exact = tuple(source.parameter * item for item in source.bulk_velocity_trace_first)
    relative = tuple(bulk_exact[i] - source.face_velocity_exact[i] for i in range(4))
    exact = rho_m * dot(relative, source.normal_exact)
    return epsilon * shape(exact, source.parameter)


def true_area_flux_raw(source: FaceSource) -> sp.Expr:
    bulk_exact = tuple(source.parameter * item for item in source.bulk_velocity_trace_first)
    relative = tuple(bulk_exact[i] - source.face_velocity_exact[i] for i in range(4))
    exact_flux = rho_m * dot(relative, source.normal_exact)
    return epsilon * shape(source.measure_exact * exact_flux, source.parameter)


def pressure_trace_raw(source: FaceSource) -> sp.Expr:
    return epsilon * source.pressure_trace_first


def mu_specific_raw(source: FaceSource) -> sp.Expr:
    exact = source.parameter * mu_theta_branch[source.branch] / source.rhobr_bg_exact
    return epsilon * shape(exact, source.parameter)


def affinity_raw(source: FaceSource) -> sp.Expr:
    return mu_specific_raw(source) - pressure_trace_raw(source) / rho_m


def traction_raw(source: FaceSource) -> sp.ImmutableMatrix:
    exact_pressure = source.parameter * source.pressure_trace_first
    exact_mu = source.parameter * mu_theta_branch[source.branch] / source.rhobr_bg_exact
    exact_affinity = exact_mu - exact_pressure / rho_m
    exact = -(exact_pressure + Lambda_X * exact_affinity) * sp.ImmutableMatrix(source.normal_exact)
    return epsilon * sp.ImmutableMatrix(shape(exact, source.parameter))


def closure_raw(source: FaceSource) -> sp.Expr:
    return sp.factor_terms(relative_flux_raw(source) - Lambda_A * affinity_raw(source) - Lambda_V * face_velocity_raw(source))


def kinematic_raw(source: FaceSource) -> Relational:
    lhs = epsilon * shape(
        dot(source.normal_exact, tuple(source.parameter * item for item in source.bulk_velocity_trace_first)),
        source.parameter,
    )
    residual = sp.factor_terms(lhs - face_velocity_raw(source) - relative_flux_raw(source) / rho_m)
    return sp.Eq(residual, 0, evaluate=False)


def conormal_raw(source: FaceSource) -> sp.Tuple:
    grad_f = tuple(symbol(f"trace_grad_f_{i}", "COORDINATE", f"generic traced-field derivative component {i}") for i in range(1, 5))
    dw_grad_f = tuple(symbol(f"d_w_trace_grad_f_{i}", "COORDINATE", f"normal jet of generic traced-field derivative {i}") for i in range(1, 5))
    traced_gradient = tuple(grad_f[i] + source.parameter * source.dh * dw_grad_f[i] for i in range(4))
    exact = dot(source.normal_exact, traced_gradient)
    background = exact.subs(source.parameter, 0)
    derivative = shape(exact, source.parameter)
    return sp.Tuple(background, epsilon * derivative, background + epsilon * derivative)


def virtual_work_cases(
    branch: str,
    dof: str,
    representative: str,
    *,
    route: str = "EULERIAN",
    ablate_direction: int | None = None,
    reverse_upper_x1: bool = False,
) -> sp.Tuple:
    face_rows = []
    total = sp.Integer(0)
    for face in FACES:
        source = build_face_source(
            branch, face, dof, representative, route=route,
            ablate_direction=ablate_direction, reverse_upper_x1=reverse_upper_x1,
        )
        exact_pressure = source.parameter * source.pressure_trace_first
        exact_mu = source.parameter * mu_theta_branch[branch] / source.rhobr_bg_exact
        exact_affinity = exact_mu - exact_pressure / rho_m
        exact_traction = -(exact_pressure + Lambda_X * exact_affinity) * sp.ImmutableMatrix(source.normal_exact)
        exact_work = source.measure_exact * dot(tuple(exact_traction), source.virtual_displacement)
        contribution = epsilon * shape(exact_work, source.parameter)
        face_rows.append(sp.Tuple(Str("UPPER" if face == 1 else "LOWER"), contribution))
        total += contribution
    return sp.Tuple(sp.Tuple(*face_rows), sp.factor_terms(total))


def profile_context() -> sp.Tuple:
    rho4_ref = rho_br / W0
    return sp.Tuple(
        sp.Eq(W_bg, W0 * (1 + eta_bg * w1_profile), evaluate=False),
        sp.Eq(mu_R_bg, mu_R * (1 + eta_bg * m1_profile), evaluate=False),
        sp.Eq(sigma_W, eta_bg * W0 / L_W, evaluate=False),
        sp.Tuple(*(sp.Eq(grad_W[i], sigma_W * w1_grad[i], evaluate=False) for i in range(3))),
        sp.Tuple(*(sp.Eq(grad_mu_R_bg[i], mu_R * sigma_W * m1_grad[i] / W0, evaluate=False) for i in range(3))),
        sp.Eq(rho4_rho4, rho4_ref, evaluate=False),
        sp.Eq(rhobr_rho4, rho4_ref * W_bg, evaluate=False),
        sp.Eq(rhobr_rhobr, rho_br, evaluate=False),
        sp.Eq(rho4_rhobr, rho_br / W_bg, evaluate=False),
        sp.Eq(e_W_bg, W0 * e_W / W_bg, evaluate=False),
    )


def build_background_density_raw() -> dict[object, object]:
    cases: dict[object, object] = {}
    for branch in BRANCHES:
        for representative in DENSITY_REPS:
            rho4_value, _ = density_pair(representative, W_bg)
            sigma_e = sp.factor_terms(rho4_value * W_bg)
            gradient = sp.ImmutableMatrix([dx(sigma_e, i) for i in range(3)])
            cases[(branch, representative)] = sp.Tuple(profile_context(), finalize(sigma_e), finalize(gradient))
    return cases


def build_geometry_quantity(
    quantity: str,
    *,
    route: str = "EULERIAN",
    ablate_direction: int | None = None,
    reverse_upper_x1: bool = False,
) -> dict[object, object]:
    cases: dict[object, object] = {}
    representatives = DENSITY_REPS if quantity in {"TRACTION", "VIRTUAL_WORK_SHAPE_DERIV", "CLOSURE_SHAPE_DERIV"} else ("RHO4_CONSTANT",)
    for branch in BRANCHES:
        for dof in DOFS:
            if quantity == "VIRTUAL_WORK_SHAPE_DERIV":
                for representative in representatives:
                    cases[(branch, dof, representative)] = virtual_work_cases(
                        branch, dof, representative, route=route,
                        ablate_direction=ablate_direction, reverse_upper_x1=reverse_upper_x1,
                    )
                continue
            for face in FACES:
                for representative in representatives:
                    source = build_face_source(
                        branch, face, dof, representative, route=route,
                        ablate_direction=ablate_direction, reverse_upper_x1=reverse_upper_x1,
                    )
                    if quantity == "FACE_NORMAL":
                        raw = face_normal_raw(source)
                    elif quantity == "CONORMAL_DERIV":
                        raw = conormal_raw(source)
                    elif quantity == "FACE_MEASURE_SHAPE_DERIV":
                        raw = face_measure_raw(source)
                    elif quantity == "FACE_VELOCITY":
                        raw = face_velocity_raw(source)
                    elif quantity == "RELATIVE_FLUX":
                        raw = relative_flux_raw(source)
                    elif quantity == "KINEMATIC_BALANCE":
                        raw = kinematic_raw(source)
                    elif quantity == "TRACTION":
                        raw = traction_raw(source)
                    elif quantity == "CLOSURE_SHAPE_DERIV":
                        raw = closure_raw(source)
                    else:
                        raise KeyError(quantity)
                    case = (branch, face, dof) if len(representatives) == 1 else (branch, face, dof, representative)
                    cases[case] = raw
    return finalize(cases)


def build_face_shift_raw(
    *, ablate_direction: int | None = None, reverse_upper_x1: bool = False
) -> dict[object, object]:
    cases: dict[object, object] = {}
    for branch in BRANCHES:
        for face in FACES:
            for dof in DOFS:
                source = build_face_source(
                    branch, face, dof, ablate_direction=ablate_direction,
                    reverse_upper_x1=reverse_upper_x1,
                )
                pressure = pressure_trace_raw(source)
                velocity = epsilon * sp.ImmutableMatrix(source.bulk_velocity_trace_first)
                rho_trace_1 = symbol(f"delta_rho_4D_face_{'plus' if face == 1 else 'minus'}", "COORDINATE", f"density perturbation at face {face}")
                dw_rho_0 = symbol(f"d_w_rho_4D_0_{'plus' if face == 1 else 'minus'}", "PREMISE", f"background density normal jet at face {face}")
                current_trace = sp.ImmutableMatrix([
                    epsilon * (j_bulk[i] + source.dh * symbol(f"d_w_j_0_{'plus' if face == 1 else 'minus'}_{i + 1}", "PREMISE", f"background-current normal jet at face {face}, component {i + 1}"))
                    for i in range(4)
                ])
                density = epsilon * (rho_trace_1 + source.dh * dw_rho_0)
                cases[(branch, face, dof)] = sp.Tuple(pressure, velocity, density, current_trace)
    return finalize(cases)


def flat_window_jets() -> tuple[sp.Expr, ...]:
    first = sp.Dummy("window_argument_plus")
    second = sp.Dummy("window_argument_minus")
    G_plus = w - W0 / 2
    G_minus = -w - W0 / 2
    O_source = sp.Function("O_window")(first, second)

    def evaluated(derivative: sp.Expr) -> sp.Expr:
        return sp.Subs(derivative, (first, second), (G_plus, G_minus))

    return (
        sp.Function("O_window")(G_plus, G_minus),
        evaluated(sp.diff(O_source, first)),
        evaluated(sp.diff(O_source, second)),
        evaluated(sp.diff(O_source, first, 2)),
        evaluated(sp.diff(O_source, first, second)),
        evaluated(sp.diff(O_source, second, 2)),
    )


WINDOW_FLAT, WINDOW_PLUS, WINDOW_MINUS, WINDOW_PP, WINDOW_PM, WINDOW_MM = flat_window_jets()


def projection_terms(
    branch: str,
    dof: str,
    representative: str,
    *,
    dynamic: bool,
    ablate_direction: int | None = None,
) -> tuple[sp.Expr, dict[str, sp.Expr]]:
    scales = source_jet_scales(ablate_direction)
    bounds = (w, -sp.oo, sp.oo)
    zero_jet_shift = -eta_bg * W0 * w1_profile / 2
    plus_bg = WINDOW_PLUS + zero_jet_shift * (WINDOW_PP + WINDOW_PM)
    minus_bg = WINDOW_MINUS + zero_jet_shift * (WINDOW_PM + WINDOW_MM)
    window_bg = WINDOW_FLAT + zero_jet_shift * (WINDOW_PLUS + WINDOW_MINUS)

    zeta_plus, _, _, _, _ = dof_fields(dof, 1)
    zeta_minus, _, _, _, _ = dof_fields(dof, -1)
    advective = sigma_W * dot(u, tuple(scales[i] * w1_grad[i] for i in range(3)))
    dh_plus = zeta_plus if branch == "LAB_HELD" else zeta_plus - advective / 2
    dh_minus = zeta_minus if branch == "LAB_HELD" else zeta_minus + advective / 2
    delta_window = -plus_bg * dh_plus + minus_bg * dh_minus
    delta_window_t = dt(delta_window)

    window_gradient = tuple(
        -scales[i] * sigma_W * w1_grad[i] * (
            WINDOW_PLUS + WINDOW_MINUS
            + zero_jet_shift * (WINDOW_PP + 2 * WINDOW_PM + WINDOW_MM)
        ) / 2
        for i in range(3)
    )
    window_normal = plus_bg - minus_bg

    parameter = sp.Dummy("projection_density_shape", real=True)
    W_alpha = branch_profile(W_bg, branch, parameter, scales)
    rho4_bg_alpha, _ = density_pair(representative, W_alpha)
    rho_shape = finalize(shape(rho4_bg_alpha + parameter * rho4_bulk_1, parameter))
    rho_shape_t = dt(rho_shape)
    rho0 = finalize(density_pair(representative, W_bg)[0])
    current_divergence = sp.Add(*(grad_j_bulk[i][i] for i in range(3)))

    if dynamic:
        origins = {
            "PROJECTED_MASS_TIME": epsilon * sp.Integral(window_bg * rho_shape_t + rho0 * delta_window_t, bounds),
            "PROJECTED_INPLANE_CURRENT": epsilon * sp.Integral(window_bg * current_divergence + sp.Add(*(j_bulk[i] * window_gradient[i] for i in range(3))), bounds),
            "DYNAMIC_WINDOW_TIME": -epsilon * sp.Integral(rho0 * delta_window_t, bounds),
            "DYNAMIC_WINDOW_INPLANE": -epsilon * sp.Integral(sp.Add(*(j_bulk[i] * window_gradient[i] for i in range(3))), bounds),
            "WINDOW_NORMAL_CURRENT": -epsilon * sp.Integral(j_bulk[3] * window_normal, bounds),
        }
    else:
        origins = {
            "PROJECTED_MASS_TIME": epsilon * sp.Integral(WINDOW_FLAT * rho4_bulk_1_t, bounds),
            "PROJECTED_INPLANE_CURRENT": epsilon * sp.Integral(WINDOW_FLAT * current_divergence, bounds),
            "DYNAMIC_WINDOW_TIME": sp.Integer(0),
            "DYNAMIC_WINDOW_INPLANE": sp.Integer(0),
            "WINDOW_NORMAL_CURRENT": -epsilon * sp.Integral(j_bulk[3] * (WINDOW_PLUS - WINDOW_MINUS), bounds),
        }
    operand = sp.Add(*origins.values())
    return sp.factor_terms(operand), origins


_PROJECTION_CACHE: dict[int | None, dict[str, dict[object, object]]] = {}


def build_projection_raw(
    quantity: str, *, ablate_direction: int | None = None
) -> dict[object, object]:
    cached = _PROJECTION_CACHE.get(ablate_direction)
    if cached is None:
        cached = {
            "PROJECTION_SHAPE_DERIV": {},
            "PROJECTION_STATIC_OPERAND": {},
            "PROJECTION_DYNAMIC_OPERAND": {},
            "PROJECTION_RESIDUAL": {},
            "PROJECTION_TERM_ORIGINS": {},
        }
        for branch in BRANCHES:
            for dof in DOFS:
                for representative in DENSITY_REPS:
                    dynamic, origins = projection_terms(branch, dof, representative, dynamic=True, ablate_direction=ablate_direction)
                    static, static_origins = projection_terms(branch, dof, representative, dynamic=False, ablate_direction=ablate_direction)
                    case = (branch, dof, representative)
                    cached["PROJECTION_SHAPE_DERIV"][case] = dynamic
                    cached["PROJECTION_DYNAMIC_OPERAND"][case] = dynamic
                    cached["PROJECTION_STATIC_OPERAND"][case] = static
                    cached["PROJECTION_RESIDUAL"][case] = sp.factor_terms(dynamic - static)
                    cached["PROJECTION_TERM_ORIGINS"][case] = sp.Tuple(casify(origins), casify(static_origins))
        _PROJECTION_CACHE[ablate_direction] = cached
    return cached[quantity]


def virtual_constraint_route(
    branch: str,
    dof: str,
    representative: str,
    *,
    route: str,
    corrupt_direct: bool = False,
    ablate_direction: int | None = None,
) -> sp.Expr:
    parameter = sp.Dummy("virtual_parameter", real=True)
    scales = source_jet_scales(ablate_direction)
    virtual_thickness = W0 * delta_v_e_W if dof == "DELTA_W" else sp.Integer(0)
    transported = dot(delta_v_u, tuple(scales[i] * grad_W[i] for i in range(3)))
    if branch == "LAB_HELD" and not (route == "EULERIAN" and corrupt_direct):
        W_anchor = W_bg + parameter * transported
    else:
        W_anchor = W_bg
    rho4_anchor, _ = density_pair(representative, W_anchor)
    thickness = W_anchor + parameter * virtual_thickness
    density = rho4_anchor * (1 + parameter * delta_v_theta)
    jacobian = determinant_jacobian(grad_delta_v_u, parameter)
    if route == "MATERIAL":
        # The face-flattening Jacobian supplies exactly the local thickness.
        flattened_mass = density * thickness * jacobian
        background_mass = density_pair(representative, W_bg)[0] * W_bg
        return epsilon * shape(flattened_mass / background_mass, parameter)
    sigma_e = density * thickness
    material_mass = sigma_e * jacobian
    background_mass = density_pair(representative, W_bg)[0] * W_bg
    return epsilon * shape(material_mass / background_mass, parameter)


def build_virtual_constraint_raw(
    *, route: str = "EULERIAN", corrupt_direct: bool = False, ablate_direction: int | None = None
) -> dict[object, object]:
    cases = {}
    for branch in BRANCHES:
        for dof in DOFS:
            for representative in DENSITY_REPS:
                value = virtual_constraint_route(
                    branch, dof, representative, route=route,
                    corrupt_direct=corrupt_direct, ablate_direction=ablate_direction,
                )
                cases[(branch, dof, representative)] = value
    return finalize(cases)


def evolution_route(
    branch: str,
    dof: str,
    representative: str,
    *,
    route: str,
    corrupt_direct: bool = False,
    ablate_direction: int | None = None,
) -> tuple[sp.Expr, dict[str, sp.Expr]]:
    parameter = sp.Dummy("evolution_shape", real=True)
    scales = source_jet_scales(ablate_direction)
    _, _, _, physical_delta_W, _ = dof_fields(dof, 1)
    W_alpha = branch_profile(W_bg, branch, parameter, scales)
    rho4_alpha, _ = density_pair(representative, W_alpha)
    sigma_exact = rho4_alpha * (1 + parameter * theta) * (W_alpha + parameter * physical_delta_W)
    sigma0 = density_pair(representative, W_bg)[0] * W_bg
    if route == "EULERIAN":
        density_time = epsilon * dt(shape(sigma_exact, parameter))
        dilatation = epsilon * sigma0 * trace_grad(grad_u_t)
        advection = epsilon * dot(u_t, tuple(dx(sigma0, i, scales) for i in range(3)))
        if corrupt_direct:
            advection = sp.Integer(0)
    else:
        evaluation_shift = parameter * dot(u, tuple(dx(sigma0, i, scales) for i in range(3)))
        material_mass = (sigma_exact + evaluation_shift) * determinant_jacobian(grad_u, parameter)
        density_time = epsilon * dt(shape(material_mass, parameter))
        dilatation = sp.Integer(0)
        advection = sp.Integer(0)
    flux_sum = sp.Integer(0)
    for face in FACES:
        source = build_face_source(branch, face, dof, representative, route=route, ablate_direction=ablate_direction)
        flux_sum += true_area_flux_raw(source)
    origins = {
        "DENSITY_TIME": density_time,
        "VELOCITY_DILATATION": dilatation,
        "BACKGROUND_ADVECTION": advection,
        "TRUE_AREA_FACE_FLUX": flux_sum,
    }
    residual = sp.Add(*origins.values())
    return finalize(residual), finalize(origins)


def build_evolution_raw(
    quantity: str,
    *, route: str = "EULERIAN", corrupt_direct: bool = False, ablate_direction: int | None = None,
) -> dict[object, object]:
    cases = {}
    for branch in BRANCHES:
        for dof in DOFS:
            for representative in DENSITY_REPS:
                residual, origins = evolution_route(
                    branch, dof, representative, route=route,
                    corrupt_direct=corrupt_direct, ablate_direction=ablate_direction,
                )
                cases[(branch, dof, representative)] = residual if quantity == "EVOLUTION_MASS_BALANCE" else origins
    return cases


def emit_primary(quantity: str, raw: Mapping[object, object], dimensions: object, key: str) -> None:
    RAW_PRIMARY[quantity] = raw
    DIMENSION_PRIMARY[quantity] = dimensions
    emit(quantity, case_payload(raw, dimensions), key=key, export=True)


def task_t0() -> None:
    raw = build_background_density_raw()
    dimensions = {case: sp.Tuple(DIM_ZERO, DIM_RHOBR, dim_div(DIM_RHOBR, DIM_L)) for case in raw}
    emit_primary("BACKGROUND_DENSITY_MAP", raw, dimensions, "background_density_map")


def task_ta() -> None:
    raw = build_geometry_quantity("FACE_NORMAL")
    emit_primary("FACE_NORMAL", raw, DIM_ZERO, "face_normal")


def task_ta_prime() -> None:
    raw = build_geometry_quantity("CONORMAL_DERIV")
    emit_primary("CONORMAL_DERIV", raw, -DIM_L, "conormal_deriv")


def task_ta_double() -> None:
    raw = build_geometry_quantity("FACE_MEASURE_SHAPE_DERIV")
    emit_primary("FACE_MEASURE_SHAPE_DERIV", raw, DIM_ZERO, "face_measure_shape_deriv")


def task_tb() -> None:
    raw = build_geometry_quantity("FACE_VELOCITY")
    emit_primary("FACE_VELOCITY", raw, DIM_VELOCITY, "face_velocity")


def task_tc() -> None:
    raw = build_geometry_quantity("RELATIVE_FLUX")
    emit_primary("RELATIVE_FLUX", raw, DIM_FLUX, "relative_flux")


def task_tc_prime() -> None:
    raw = build_geometry_quantity("KINEMATIC_BALANCE")
    emit_primary("KINEMATIC_BALANCE", raw, DIM_VELOCITY, "kinematic_balance")


def task_td() -> None:
    traction = build_geometry_quantity("TRACTION")
    work = build_geometry_quantity("VIRTUAL_WORK_SHAPE_DERIV")
    emit_primary("TRACTION", traction, DIM_PRESSURE, "traction")
    emit_primary("VIRTUAL_WORK_SHAPE_DERIV", work, DIM_WORK, "virtual_work_shape_deriv")


def task_te() -> None:
    raw = build_face_shift_raw()
    dimensions = sp.Tuple(DIM_PRESSURE, DIM_VELOCITY, DIM_RHO4, DIM_FLUX)
    emit_primary("FACE_SHIFT", raw, dimensions, "face_shift")


def task_tf() -> None:
    quantities = (
        "PROJECTION_SHAPE_DERIV", "PROJECTION_STATIC_OPERAND", "PROJECTION_DYNAMIC_OPERAND",
        "PROJECTION_RESIDUAL", "PROJECTION_TERM_ORIGINS",
    )
    for quantity in quantities:
        raw = build_projection_raw(quantity)
        emit_primary(quantity, raw, DIM_MASS_RATE, quantity.lower())


def task_tg() -> None:
    raw = build_virtual_constraint_raw()
    emit_primary("VIRTUAL_CONSTRAINT", raw, DIM_ZERO, "virtual_constraint")


def task_th() -> None:
    balance = build_evolution_raw("EVOLUTION_MASS_BALANCE")
    origins = build_evolution_raw("EVOLUTION_TERM_ORIGINS")
    emit_primary("EVOLUTION_MASS_BALANCE", balance, DIM_MASS_RATE, "evolution_mass_balance")
    emit_primary("EVOLUTION_TERM_ORIGINS", origins, DIM_MASS_RATE, "evolution_term_origins")


def task_ti() -> None:
    raw = build_geometry_quantity("CLOSURE_SHAPE_DERIV")
    emit_primary("CLOSURE_SHAPE_DERIV", raw, DIM_FLUX, "closure_shape_deriv")


def supplied_face_maps() -> tuple[object, object]:
    X = tuple(symbol(f"X_{i}", "COORDINATE", f"material coordinate {i}") for i in range(1, 4))
    x_of_X = tuple(X[i] + epsilon * u[i] for i in range(3))
    W_at_xX = sp.Function("W_bg")(*x_of_X)
    W_at_X = sp.Function("W_bg")(*X)
    lab = {}
    material = {}
    for face in FACES:
        zeta = zeta_c + face * W0 * e_W / 2
        lab[face] = sp.Tuple(sp.Tuple(*x_of_X), face * W_at_xX / 2 + epsilon * zeta)
        material[face] = sp.Tuple(sp.Tuple(*x_of_X), face * W_at_X / 2 + epsilon * zeta)
    return lab, material


def emit_supplied_objects() -> None:
    state = {
        (branch, representative): sp.Tuple(
            profile_context(), sp.Eq(theta_0, 0, evaluate=False),
            sp.Eq(symbol(f"V_0_{branch}_{representative}", "PREMISE", "supplied zero background face velocity"), 0, evaluate=False),
            sp.Eq(symbol(f"J_0_{branch}_{representative}", "PREMISE", "supplied zero background face flux"), 0, evaluate=False),
            sp.Eq(symbol(f"A_0_{branch}_{representative}", "PREMISE", "supplied zero background affinity"), 0, evaluate=False),
        )
        for branch in BRANCHES for representative in DENSITY_REPS
    }
    premise = sp.Tuple(
        Str("SUPPORT_STABILISED_BACKGROUND"),
        sp.Tuple(support_body, support_face_plus, support_face_minus),
        Str("STATIONARITY_NOT_TESTED_IN_S11C_A"),
        Str("S11CB_ADMISSIBILITY_PACKAGE_RESERVED"),
    )
    lab, material = supplied_face_maps()
    emit("BACKGROUND_STATE", state)
    emit("ADMISSIBILITY_PREMISE", premise)
    emit("FACE_MAP_LAB_HELD", lab)
    emit("FACE_MAP_MATERIAL_ADVECTED", material)


def object_difference(left: object, right: object) -> object:
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        if not isinstance(left, Mapping) or not isinstance(right, Mapping):
            return sp.Tuple(casify(left), casify(right))
        keys = tuple(dict.fromkeys((*left.keys(), *right.keys())))
        return {key: object_difference(left.get(key, sp.Integer(0)), right.get(key, sp.Integer(0))) for key in keys}
    if isinstance(left, Relational) or isinstance(right, Relational):
        if isinstance(left, Relational) and isinstance(right, Relational):
            return sp.factor_terms((left.lhs - left.rhs) - (right.lhs - right.rhs))
        return sp.Tuple(casify(left), casify(right))
    if isinstance(left, sp.MatrixBase) and isinstance(right, sp.MatrixBase) and left.shape == right.shape:
        return sp.ImmutableMatrix(left - right)
    if isinstance(left, (tuple, list, sp.Tuple)) and isinstance(right, (tuple, list, sp.Tuple)) and len(left) == len(right):
        return sp.Tuple(*(casify(object_difference(a, b)) for a, b in zip(left, right)))
    left_cas = casify(left)
    right_cas = casify(right)
    try:
        return sp.factor_terms(left_cas - right_cas)
    except Exception:
        status, detail = total_compare(left_cas, right_cas)
        return sp.Tuple(Str(status), casify(detail))


def task_rep_invariance() -> None:
    eulerian: dict[object, object] = {}
    material: dict[object, object] = {}
    for quantity in ("RELATIVE_FLUX", "TRACTION", "CLOSURE_SHAPE_DERIV"):
        eulerian[quantity] = RAW_PRIMARY[quantity]
        material[quantity] = build_geometry_quantity(quantity, route="MATERIAL")
    eulerian["VIRTUAL_WORK_SHAPE_DERIV"] = RAW_PRIMARY["VIRTUAL_WORK_SHAPE_DERIV"]
    material["VIRTUAL_WORK_SHAPE_DERIV"] = build_geometry_quantity("VIRTUAL_WORK_SHAPE_DERIV", route="MATERIAL")
    eulerian["VIRTUAL_CONSTRAINT"] = RAW_PRIMARY["VIRTUAL_CONSTRAINT"]
    material["VIRTUAL_CONSTRAINT"] = build_virtual_constraint_raw(route="MATERIAL")
    eulerian["EVOLUTION_MASS_BALANCE"] = RAW_PRIMARY["EVOLUTION_MASS_BALANCE"]
    material["EVOLUTION_MASS_BALANCE"] = build_evolution_raw("EVOLUTION_MASS_BALANCE", route="MATERIAL")
    residual = object_difference(eulerian, material)
    emit("REP_INVARIANCE_EULERIAN_OPERAND", eulerian)
    emit("REP_INVARIANCE_MATERIAL_OPERAND", material)
    emit("REP_INVARIANCE_RESIDUAL", residual)


def task_independence() -> None:
    base: dict[object, object] = {}
    corrupted: dict[object, object] = {}
    for quantity in ("RELATIVE_FLUX", "TRACTION", "CLOSURE_SHAPE_DERIV", "VIRTUAL_WORK_SHAPE_DERIV"):
        base[quantity] = build_geometry_quantity(quantity)
        corrupted[quantity] = build_geometry_quantity(quantity, reverse_upper_x1=True)
    base["VIRTUAL_CONSTRAINT"] = build_virtual_constraint_raw()
    corrupted["VIRTUAL_CONSTRAINT"] = build_virtual_constraint_raw(corrupt_direct=True)
    base["EVOLUTION_MASS_BALANCE"] = build_evolution_raw("EVOLUTION_MASS_BALANCE")
    corrupted["EVOLUTION_MASS_BALANCE"] = build_evolution_raw("EVOLUTION_MASS_BALANCE", corrupt_direct=True)
    emit("CONTROL_INDEPENDENCE_BASE_OPERAND", base)
    emit("CONTROL_INDEPENDENCE_CORRUPTED_OPERAND", corrupted)
    emit("CONTROL_INDEPENDENCE_RESIDUAL", object_difference(base, corrupted))


FORM_BUILDERS = {
    "FACE_NORMAL": lambda direction: build_geometry_quantity("FACE_NORMAL", ablate_direction=direction),
    "CONORMAL_DERIV": lambda direction: build_geometry_quantity("CONORMAL_DERIV", ablate_direction=direction),
    "FACE_MEASURE_SHAPE_DERIV": lambda direction: build_geometry_quantity("FACE_MEASURE_SHAPE_DERIV", ablate_direction=direction),
    "FACE_VELOCITY": lambda direction: build_geometry_quantity("FACE_VELOCITY", ablate_direction=direction),
    "RELATIVE_FLUX": lambda direction: build_geometry_quantity("RELATIVE_FLUX", ablate_direction=direction),
    "KINEMATIC_BALANCE": lambda direction: build_geometry_quantity("KINEMATIC_BALANCE", ablate_direction=direction),
    "TRACTION": lambda direction: build_geometry_quantity("TRACTION", ablate_direction=direction),
    "VIRTUAL_WORK_SHAPE_DERIV": lambda direction: build_geometry_quantity("VIRTUAL_WORK_SHAPE_DERIV", ablate_direction=direction),
    "FACE_SHIFT": lambda direction: build_face_shift_raw(ablate_direction=direction),
    "PROJECTION_SHAPE_DERIV": lambda direction: build_projection_raw("PROJECTION_SHAPE_DERIV", ablate_direction=direction),
    "PROJECTION_STATIC_OPERAND": lambda direction: build_projection_raw("PROJECTION_STATIC_OPERAND", ablate_direction=direction),
    "PROJECTION_DYNAMIC_OPERAND": lambda direction: build_projection_raw("PROJECTION_DYNAMIC_OPERAND", ablate_direction=direction),
    "PROJECTION_RESIDUAL": lambda direction: build_projection_raw("PROJECTION_RESIDUAL", ablate_direction=direction),
    "PROJECTION_TERM_ORIGINS": lambda direction: build_projection_raw("PROJECTION_TERM_ORIGINS", ablate_direction=direction),
    "VIRTUAL_CONSTRAINT": lambda direction: build_virtual_constraint_raw(ablate_direction=direction),
    "EVOLUTION_MASS_BALANCE": lambda direction: build_evolution_raw("EVOLUTION_MASS_BALANCE", ablate_direction=direction),
    "EVOLUTION_TERM_ORIGINS": lambda direction: build_evolution_raw("EVOLUTION_TERM_ORIGINS", ablate_direction=direction),
    "CLOSURE_SHAPE_DERIV": lambda direction: build_geometry_quantity("CLOSURE_SHAPE_DERIV", ablate_direction=direction),
}


def task_form_control() -> None:
    base: dict[object, object] = {}
    ablated: dict[object, object] = {}
    residual: dict[object, object] = {}
    for direction in (1, 2, 3):
        for quantity, builder in FORM_BUILDERS.items():
            baseline = RAW_PRIMARY[quantity]
            changed = builder(direction)
            for case, baseline_value in baseline.items():
                changed_value = changed[case]
                case_tuple = case if isinstance(case, tuple) else (case,)
                if len(case_tuple) >= 2 and case_tuple[1] in FACES:
                    keyed_case = case_tuple
                else:
                    keyed_case = (case_tuple[0], "BOTH_FACES", *case_tuple[1:])
                key = (quantity, *keyed_case, direction)
                base[key] = baseline_value
                ablated[key] = changed_value
                residual[key] = object_difference(baseline_value, changed_value)
    emit("CONTROL_FORM_BASE_OPERAND", base)
    emit("CONTROL_FORM_ABLATED_OPERAND", ablated)
    emit("CONTROL_FORM_RESIDUAL", residual)


def uniform_substitution(value: object) -> object:
    substitutions = {
        eta_bg: 0,
        sigma_W: 0,
        W_bg: W0,
        **{grad_W[i]: 0 for i in range(3)},
        **{hess_W[i][j]: 0 for i in range(3) for j in range(3)},
    }
    return map_object(value, lambda item: sp.factor_terms(item.subs(substitutions)))


def uniform_reference_geometry(quantity: str) -> dict[object, object]:
    """Independent flat-face construction from the supplied S11b laws."""
    cases: dict[object, object] = {}
    reps = DENSITY_REPS if quantity in {"TRACTION", "VIRTUAL_WORK_SHAPE_DERIV", "CLOSURE_SHAPE_DERIV"} else ("RHO4_CONSTANT",)
    for branch in BRANCHES:
        for dof in DOFS:
            if quantity == "VIRTUAL_WORK_SHAPE_DERIV":
                for representative in reps:
                    contributions = []
                    for face in FACES:
                        zeta, zeta_t, grad_zeta, _, virtual_delta_W = dof_fields(dof, face)
                        normal = sp.ImmutableMatrix([0, 0, 0, face])
                        p = epsilon * ((delta_p_plus if face == 1 else delta_p_minus) + zeta * dw_delta_p_0[face])
                        mu_s = epsilon * mu_theta_branch[branch] / rho_br
                        tr = -(p + Lambda_X * (mu_s - p / rho_m)) * normal
                        virtual_zeta = delta_v_zeta_c if dof == "ZETA_C" else face * virtual_delta_W / 2
                        virtual = sp.ImmutableMatrix([*delta_v_u, virtual_zeta])
                        contributions.append(sp.Tuple(Str("UPPER" if face == 1 else "LOWER"), dot(tuple(tr), tuple(virtual))))
                    cases[(branch, dof, representative)] = sp.Tuple(sp.Tuple(*contributions), sp.Add(*(row[1] for row in contributions)))
                continue
            for face in FACES:
                zeta, zeta_t, grad_zeta, _, _ = dof_fields(dof, face)
                normal = sp.ImmutableMatrix([0, 0, 0, face])
                dn = sp.ImmutableMatrix([-face * component for component in grad_zeta] + [0])
                V = epsilon * face * zeta_t
                vb = delta_v_bulk[face]
                flux = epsilon * rho_m * face * (vb[3] + zeta * dw_v_bulk_0[face][3] - zeta_t)
                p = epsilon * ((delta_p_plus if face == 1 else delta_p_minus) + zeta * dw_delta_p_0[face])
                mu_s = epsilon * mu_theta_branch[branch] / rho_br
                affinity = mu_s - p / rho_m
                if quantity == "FACE_NORMAL":
                    raw = sp.Tuple(normal, epsilon * dn, normal + epsilon * dn)
                elif quantity == "CONORMAL_DERIV":
                    grad_f = tuple(symbol(f"trace_grad_f_{i}", "COORDINATE", f"generic traced-field derivative component {i}") for i in range(1, 5))
                    dw_grad_f = tuple(symbol(f"d_w_trace_grad_f_{i}", "COORDINATE", f"normal jet of generic traced-field derivative {i}") for i in range(1, 5))
                    background = face * grad_f[3]
                    derivative = -face * dot(tuple(grad_zeta), tuple(grad_f[:3])) + face * zeta * dw_grad_f[3]
                    raw = sp.Tuple(background, epsilon * derivative, background + epsilon * derivative)
                elif quantity == "FACE_MEASURE_SHAPE_DERIV":
                    raw = sp.Tuple(1, 0, 1)
                elif quantity == "FACE_VELOCITY":
                    raw = V
                elif quantity == "RELATIVE_FLUX":
                    raw = flux
                elif quantity == "KINEMATIC_BALANCE":
                    raw = sp.Eq(epsilon * face * (vb[3] + zeta * dw_v_bulk_0[face][3]) - V - flux / rho_m, 0, evaluate=False)
                elif quantity == "TRACTION":
                    raw = -(p + Lambda_X * affinity) * normal
                elif quantity == "CLOSURE_SHAPE_DERIV":
                    raw = flux - Lambda_A * affinity - Lambda_V * V
                else:
                    raise KeyError(quantity)
                for representative in reps:
                    key = (branch, face, dof) if len(reps) == 1 else (branch, face, dof, representative)
                    cases[key] = sp.factor_terms(raw) if isinstance(raw, sp.Expr) else raw
    return cases


def uniform_face_shift_reference() -> dict[object, object]:
    cases: dict[object, object] = {}
    for branch in BRANCHES:
        for face in FACES:
            for dof in DOFS:
                zeta, _, _, _, _ = dof_fields(dof, face)
                pressure_face = delta_p_plus if face == 1 else delta_p_minus
                pressure = epsilon * (pressure_face + zeta * dw_delta_p_0[face])
                velocity = epsilon * sp.ImmutableMatrix([
                    delta_v_bulk[face][i] + zeta * dw_v_bulk_0[face][i] for i in range(4)
                ])
                rho_trace_1 = symbol(f"delta_rho_4D_face_{'plus' if face == 1 else 'minus'}", "COORDINATE", f"density perturbation at face {face}")
                dw_rho_0 = symbol(f"d_w_rho_4D_0_{'plus' if face == 1 else 'minus'}", "PREMISE", f"background density normal jet at face {face}")
                density = epsilon * (rho_trace_1 + zeta * dw_rho_0)
                current = epsilon * sp.ImmutableMatrix([
                    j_bulk[i] + zeta * symbol(f"d_w_j_0_{'plus' if face == 1 else 'minus'}_{i + 1}", "PREMISE", f"background-current normal jet at face {face}, component {i + 1}")
                    for i in range(4)
                ])
                cases[(branch, face, dof)] = sp.Tuple(pressure, velocity, density, current)
    return cases


def uniform_projection_reference(quantity: str) -> dict[object, object]:
    cases: dict[object, object] = {}
    bounds = (w, -sp.oo, sp.oo)
    current_divergence = sp.Add(*(grad_j_bulk[i][i] for i in range(3)))
    rho0 = rho_br / W0
    for branch in BRANCHES:
        for dof in DOFS:
            zeta_plus, _, _, _, _ = dof_fields(dof, 1)
            zeta_minus, _, _, _, _ = dof_fields(dof, -1)
            delta_window = -WINDOW_PLUS * zeta_plus + WINDOW_MINUS * zeta_minus
            delta_window_t = dt(delta_window)
            dynamic_origins = {
                "PROJECTED_MASS_TIME": epsilon * sp.Integral(WINDOW_FLAT * rho4_bulk_1_t + rho0 * delta_window_t, bounds),
                "PROJECTED_INPLANE_CURRENT": epsilon * sp.Integral(WINDOW_FLAT * current_divergence, bounds),
                "DYNAMIC_WINDOW_TIME": -epsilon * sp.Integral(rho0 * delta_window_t, bounds),
                "DYNAMIC_WINDOW_INPLANE": sp.Integer(0),
                "WINDOW_NORMAL_CURRENT": -epsilon * sp.Integral(j_bulk[3] * (WINDOW_PLUS - WINDOW_MINUS), bounds),
            }
            static_origins = {
                "PROJECTED_MASS_TIME": epsilon * sp.Integral(WINDOW_FLAT * rho4_bulk_1_t, bounds),
                "PROJECTED_INPLANE_CURRENT": epsilon * sp.Integral(WINDOW_FLAT * current_divergence, bounds),
                "DYNAMIC_WINDOW_TIME": sp.Integer(0),
                "DYNAMIC_WINDOW_INPLANE": sp.Integer(0),
                "WINDOW_NORMAL_CURRENT": -epsilon * sp.Integral(j_bulk[3] * (WINDOW_PLUS - WINDOW_MINUS), bounds),
            }
            dynamic = sp.Add(*dynamic_origins.values())
            static = sp.Add(*static_origins.values())
            for representative in DENSITY_REPS:
                case = (branch, dof, representative)
                if quantity in {"PROJECTION_SHAPE_DERIV", "PROJECTION_DYNAMIC_OPERAND"}:
                    cases[case] = dynamic
                elif quantity == "PROJECTION_STATIC_OPERAND":
                    cases[case] = static
                elif quantity == "PROJECTION_RESIDUAL":
                    cases[case] = sp.factor_terms(dynamic - static)
                elif quantity == "PROJECTION_TERM_ORIGINS":
                    cases[case] = sp.Tuple(casify(dynamic_origins), casify(static_origins))
                else:
                    raise KeyError(quantity)
    return cases


def uniform_evolution_reference(quantity: str) -> dict[object, object]:
    cases: dict[object, object] = {}
    flat_flux = uniform_reference_geometry("RELATIVE_FLUX")
    for branch in BRANCHES:
        for dof in DOFS:
            thickness_rate = e_W_t if dof == "DELTA_W" else sp.Integer(0)
            flux_sum = flat_flux[(branch, 1, dof)] + flat_flux[(branch, -1, dof)]
            origins = {
                "DENSITY_TIME": epsilon * rho_br * (theta_t + thickness_rate),
                "VELOCITY_DILATATION": epsilon * rho_br * trace_grad(grad_u_t),
                "BACKGROUND_ADVECTION": sp.Integer(0),
                "TRUE_AREA_FACE_FLUX": flux_sum,
            }
            for representative in DENSITY_REPS:
                case = (branch, dof, representative)
                cases[case] = sp.Add(*origins.values()) if quantity == "EVOLUTION_MASS_BALANCE" else origins
    return cases


def independent_uniform_reference(quantity: str) -> object:
    if quantity == "BACKGROUND_DENSITY_MAP":
        cases = {}
        for branch in BRANCHES:
            for representative in DENSITY_REPS:
                cases[(branch, representative)] = sp.Tuple(profile_context(), rho_br, sp.zeros(3, 1))
        return uniform_substitution(cases)
    if quantity in {"FACE_NORMAL", "CONORMAL_DERIV", "FACE_MEASURE_SHAPE_DERIV", "FACE_VELOCITY", "RELATIVE_FLUX", "KINEMATIC_BALANCE", "TRACTION", "VIRTUAL_WORK_SHAPE_DERIV", "CLOSURE_SHAPE_DERIV"}:
        return uniform_reference_geometry(quantity)
    if quantity == "FACE_SHIFT":
        return uniform_face_shift_reference()
    if quantity.startswith("PROJECTION_"):
        return uniform_projection_reference(quantity)
    if quantity == "VIRTUAL_CONSTRAINT":
        cases = {}
        for branch in BRANCHES:
            for dof in DOFS:
                for representative in DENSITY_REPS:
                    thickness = delta_v_e_W if dof == "DELTA_W" else 0
                    cases[(branch, dof, representative)] = epsilon * (delta_v_theta + thickness + trace_grad(grad_delta_v_u))
        return cases
    if quantity.startswith("EVOLUTION_"):
        return uniform_evolution_reference(quantity)
    raise KeyError(quantity)


def task_uniform_control() -> None:
    s11ca: dict[object, object] = {}
    s11b: dict[object, object] = {}
    for quantity, raw in RAW_PRIMARY.items():
        s11ca[quantity] = uniform_substitution(raw)
        s11b[quantity] = independent_uniform_reference(quantity)
    emit("UNIFORM_LIMIT_S11CA_OPERAND", s11ca)
    emit("UNIFORM_LIMIT_S11B_OPERAND", s11b)
    emit("UNIFORM_LIMIT_RESIDUAL", object_difference(s11ca, s11b))


def task_homogeneity() -> None:
    dimensions: dict[object, object] = {
        "BACKGROUND_STATE": sp.Tuple(DIM_L, dim_mul(DIM_M, -DIM_L, -2 * DIM_T), DIM_RHO4, DIM_RHOBR, DIM_ZERO),
        "ADMISSIBILITY_PREMISE": DIM_ZERO,
        "FACE_MAP_LAB_HELD": DIM_L,
        "FACE_MAP_MATERIAL_ADVECTED": DIM_L,
    }
    dimensions.update({quantity: DIMENSION_PRIMARY[quantity] for quantity in RAW_PRIMARY})
    representation_dimensions = {
        quantity: DIMENSION_PRIMARY[quantity]
        for quantity in (
            "RELATIVE_FLUX", "TRACTION", "VIRTUAL_WORK_SHAPE_DERIV",
            "VIRTUAL_CONSTRAINT", "EVOLUTION_MASS_BALANCE", "CLOSURE_SHAPE_DERIV",
        )
    }
    dimensions.update({
        "REP_INVARIANCE_EULERIAN_OPERAND": representation_dimensions,
        "REP_INVARIANCE_MATERIAL_OPERAND": representation_dimensions,
        "REP_INVARIANCE_RESIDUAL": representation_dimensions,
        "CONTROL_INDEPENDENCE_BASE_OPERAND": representation_dimensions,
        "CONTROL_INDEPENDENCE_CORRUPTED_OPERAND": representation_dimensions,
        "CONTROL_INDEPENDENCE_RESIDUAL": representation_dimensions,
        "CONTROL_FORM_BASE_OPERAND": {quantity: DIMENSION_PRIMARY[quantity] for quantity in FORM_BUILDERS},
        "CONTROL_FORM_ABLATED_OPERAND": {quantity: DIMENSION_PRIMARY[quantity] for quantity in FORM_BUILDERS},
        "CONTROL_FORM_RESIDUAL": {quantity: DIMENSION_PRIMARY[quantity] for quantity in FORM_BUILDERS},
        "UNIFORM_LIMIT_S11CA_OPERAND": {quantity: DIMENSION_PRIMARY[quantity] for quantity in RAW_PRIMARY},
        "UNIFORM_LIMIT_S11B_OPERAND": {quantity: DIMENSION_PRIMARY[quantity] for quantity in RAW_PRIMARY},
        "UNIFORM_LIMIT_RESIDUAL": {quantity: DIMENSION_PRIMARY[quantity] for quantity in RAW_PRIMARY},
        "HOMOGENEITY_BASE_OPERAND": DIM_ZERO,
        "HOMOGENEITY_CONTROL_OPERAND": DIM_ZERO,
        "HOMOGENEITY_RESIDUAL": DIM_ZERO,
    })
    velocity_operand = delta_v_bulk[1][3]
    affinity_operand = mu_theta_branch["LAB_HELD"] / rho_br - delta_p_plus / rho_m
    source_terms = {
        "RELATIVE_FLUX_LAW": (rho_m * velocity_operand, rho_m * u_t[0]),
        "KINEMATIC_BALANCE": (velocity_operand, u_t[0], rho_m * velocity_operand / rho_m),
        "TRACTION": (delta_p_plus, Lambda_X0 * affinity_operand),
        "VIRTUAL_WORK": (delta_p_plus * delta_v_zeta_c, Lambda_X0 * affinity_operand * delta_v_zeta_c),
        "PROJECTION": (rho4_bulk_1_t, grad_j_bulk[0][0], grad_j_bulk[3][0]),
        "VIRTUAL_CONSTRAINT": (delta_v_theta, delta_v_e_W, grad_u_t[0][0] * tau_A),
        "EVOLUTION_BALANCE": (rho_br * theta_t, rho_br * grad_u_t[0][0], j_bulk[3]),
        "CLOSURE": (rho_m * velocity_operand, Lambda_A0 * affinity_operand, Lambda_V0 * u_t[0]),
    }
    corrupt_terms = dict(source_terms)
    corrupt_terms["CLOSURE"] = (
        source_terms["CLOSURE"][0],
        source_terms["CLOSURE"][1],
        source_terms["CLOSURE"][2] * W0,
    )

    def traced(terms: Mapping[str, tuple[sp.Expr, ...]]) -> dict[str, sp.Tuple]:
        return {
            name: sp.Tuple(*(sp.Tuple(term, dimension_of(term)) for term in row))
            for name, row in terms.items()
        }

    base = traced(source_terms)
    control = traced(corrupt_terms)
    emit("DIMENSIONS", dimensions)
    emit("HOMOGENEITY_BASE_OPERAND", base)
    emit("HOMOGENEITY_CONTROL_OPERAND", control)
    emit("HOMOGENEITY_RESIDUAL", object_difference(base, control))


def total_compare(left: object, right: object) -> tuple[str, object]:
    """F9's three-valued comparison, total over imported row shapes."""
    if isinstance(left, bool) or isinstance(right, bool):
        if isinstance(left, bool) and isinstance(right, bool):
            return (PROVED_EQUAL if left is right else PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right)))
        return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
    if isinstance(left, str) or isinstance(right, str):
        if isinstance(left, str) and isinstance(right, str):
            return (PROVED_EQUAL if left == right else PROVED_DIFFERENT, sp.Tuple(Str(left), Str(right)))
        return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        if not isinstance(left, Mapping) or not isinstance(right, Mapping) or set(left) != set(right):
            return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
        ordered = sorted(left, key=str)
        comparisons = [total_compare(left[key], right[key]) for key in ordered]
        statuses = [status for status, _ in comparisons]
        detail = sp.Tuple(*(sp.Tuple(Str(str(key)), Str(status), casify(operand)) for key, (status, operand) in zip(ordered, comparisons)))
        if all(status == PROVED_EQUAL for status in statuses):
            return PROVED_EQUAL, detail
        if any(status == PROVED_DIFFERENT for status in statuses):
            return PROVED_DIFFERENT, detail
        return UNDECIDED, detail
    left_container = isinstance(left, (tuple, list, sp.Tuple))
    right_container = isinstance(right, (tuple, list, sp.Tuple))
    if left_container or right_container:
        if not left_container or not right_container or len(left) != len(right):
            return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
        comparisons = [total_compare(a, b) for a, b in zip(left, right)]
        statuses = [status for status, _ in comparisons]
        detail = sp.Tuple(*(sp.Tuple(Str(status), casify(operand)) for status, operand in comparisons))
        if all(status == PROVED_EQUAL for status in statuses):
            return PROVED_EQUAL, detail
        if any(status == PROVED_DIFFERENT for status in statuses):
            return PROVED_DIFFERENT, detail
        return UNDECIDED, detail
    if isinstance(left, sp.MatrixBase) or isinstance(right, sp.MatrixBase):
        if not isinstance(left, sp.MatrixBase) or not isinstance(right, sp.MatrixBase) or left.shape != right.shape:
            return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
        comparisons = [total_compare(a, b) for a, b in zip(left, right)]
        statuses = [status for status, _ in comparisons]
        difference = sp.ImmutableMatrix(left - right)
        if all(status == PROVED_EQUAL for status in statuses):
            return PROVED_EQUAL, difference
        if any(status == PROVED_DIFFERENT for status in statuses):
            return PROVED_DIFFERENT, difference
        return UNDECIDED, difference
    left = casify(left)
    right = casify(right)
    if isinstance(left, sp.Symbol) and isinstance(right, sp.Symbol):
        equal = left.name == right.name and left.assumptions0 == right.assumptions0
        return (PROVED_EQUAL if equal else PROVED_DIFFERENT, sp.Tuple(left, right))
    if isinstance(left, Str) or isinstance(right, Str):
        equal = isinstance(left, Str) and isinstance(right, Str) and sp.srepr(left) == sp.srepr(right)
        return (PROVED_EQUAL if equal else PROVED_DIFFERENT, sp.Tuple(left, right))
    if isinstance(left, (Relational, Boolean)) or isinstance(right, (Relational, Boolean)):
        if not isinstance(left, Boolean) or not isinstance(right, Boolean):
            return PROVED_DIFFERENT, sp.Tuple(left, right)
        if sp.srepr(left) == sp.srepr(right):
            return PROVED_EQUAL, sp.Tuple(left, right)
        if sp.count_ops(left) + sp.count_ops(right) > 20:
            return UNDECIDED, sp.Tuple(left, right)
        equivalent = sp.simplify_logic(sp.Equivalent(left, right))
        if equivalent is sp.true:
            return PROVED_EQUAL, equivalent
        if equivalent is sp.false:
            return PROVED_DIFFERENT, equivalent
        return UNDECIDED, equivalent
    if isinstance(left, sp.Basic) and isinstance(right, sp.Basic):
        if sp.srepr(left) == sp.srepr(right):
            return PROVED_EQUAL, sp.Tuple(left, right)
        try:
            equality = left.equals(right)
        except Exception:
            equality = None
        if equality is True:
            return PROVED_EQUAL, sp.Tuple(left, right)
        if equality is False:
            return PROVED_DIFFERENT, sp.Tuple(left, right)
        if sp.count_ops(left) + sp.count_ops(right) > 40:
            return UNDECIDED, sp.Tuple(left, right)
        try:
            residual = sp.simplify(left - right)
        except Exception:
            return UNDECIDED, sp.Tuple(left, right)
        if residual == 0:
            return PROVED_EQUAL, residual
        if not residual.free_symbols:
            return PROVED_DIFFERENT, residual
        return UNDECIDED, residual
    try:
        equality = left == right
    except Exception:
        equality = None
    if equality is True:
        return PROVED_EQUAL, sp.Tuple(casify(left), casify(right))
    if equality is False:
        return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
    return UNDECIDED, sp.Tuple(casify(left), casify(right))


def atoms_of(value: object) -> set[sp.Symbol]:
    if isinstance(value, sp.Symbol):
        return {value}
    if isinstance(value, Mapping):
        result: set[sp.Symbol] = set()
        for key, item in value.items():
            result.update(atoms_of(key))
            result.update(atoms_of(item))
        return result
    if isinstance(value, sp.MatrixBase):
        result = set()
        for item in value:
            result.update(atoms_of(item))
        return result
    if isinstance(value, sp.Basic):
        return set(value.free_symbols)
    if isinstance(value, (tuple, list)):
        result = set()
        for item in value:
            result.update(atoms_of(item))
        return result
    return set()


def add_free_symbol_candidates(candidates: list[CandidateRow]) -> None:
    seen = {candidate.key for candidate in candidates}
    free_symbols: set[sp.Symbol] = set()
    for candidate in candidates:
        free_symbols.update(atoms_of(candidate.value))
    for free_symbol in sorted(free_symbols, key=sp.default_sort_key):
        if free_symbol.name in seen:
            continue
        metadata = DECLARED_SYMBOLS.get(free_symbol)
        if metadata is None:
            issue(f"unclassifiable free symbol {sp.srepr(free_symbol)}")
            continue
        candidates.append(CandidateRow(free_symbol.name, free_symbol, metadata["class"], "FREE_SYMBOL_POPULATION", metadata["description"]))
        seen.add(free_symbol.name)


def make_record(candidate: CandidateRow) -> dict[str, object]:
    record: dict[str, object] = {
        "display": display(candidate.value),
        "value": casify(candidate.value),
        "value_kind": "COMPUTED_OBJECT",
        "class": candidate.class_tag,
        "step": "S11c-a",
    }
    if candidate.description is not None:
        record["description"] = candidate.description
    return record


def merged_export() -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    candidates = list(EMITTER.export_candidates)
    add_free_symbol_candidates(candidates)
    candidate_keys = sp.Tuple(*(Str(candidate.key) for candidate in candidates))
    imported_matching_keys = sp.Tuple(*(Str(candidate.key) for candidate in candidates if candidate.key in INCOMING_LEDGER))
    merged = {name: dict(record) for name, record in INCOMING_LEDGER.items()}
    routing = []
    f9c_rows = []
    seen: set[str] = set()
    for candidate in candidates:
        if candidate.key in seen:
            raise RuntimeError(f"two S11c-a candidates have the same key {candidate.key!r}")
        seen.add(candidate.key)
        record = make_record(candidate)
        imported_record = INCOMING_LEDGER.get(candidate.key)
        if imported_record is None:
            write_key = candidate.key
            status = "ABSENT"
            comparison_operand = sp.Tuple(record["value"], sp.Tuple())
        else:
            value_status, comparison_operand = total_compare(record["value"], imported_record.get("value"))
            class_status = PROVED_EQUAL if record["class"] == imported_record.get("class") else PROVED_DIFFERENT
            if value_status == PROVED_EQUAL and class_status == PROVED_EQUAL:
                write_key = candidate.key
                status = PROVED_EQUAL
                record["f9_operands"] = sp.Tuple(record["value"], imported_record.get("value"))
                prior_steps = list(imported_record.get("corroborated_steps", ()))
                if not prior_steps and imported_record.get("step"):
                    prior_steps.append(imported_record["step"])
                prior_steps.append("S11c-a")
                record["corroborated_steps"] = tuple(dict.fromkeys(prior_steps))
            else:
                write_key = f"s11c_a_{candidate.key}"
                status = value_status if value_status != PROVED_EQUAL else PROVED_DIFFERENT
                f9c_rows.append(sp.Tuple(Str(candidate.key), Str(write_key), Str(status), casify(comparison_operand)))
                issue(f"F9c write {write_key} for {candidate.key}: {status}")
        routing.append(sp.Tuple(Str(candidate.source_tag), Str(candidate.key), Str(write_key), Str(status), casify(comparison_operand)))
        merged[write_key] = record
    diagnostics = {
        "candidate_keys": candidate_keys,
        "imported_matching_keys": imported_matching_keys,
        "routing": sp.Tuple(*routing),
        "f9c": sp.Tuple(*f9c_rows),
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
        "# S11c_a_exports.py -- GENERATED by S11c_a_interface_geometry_sympy_audit.py. Do not edit.",
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
        f"    'S11c_a_interface_geometry_sympy_audit.py': {sha256(SCRIPT_PATH)!r},",
        f"    'S11b_exports.py': {sha256(SCRIPT_DIR / 'S11b_exports.py')!r},",
        f"    'S11c_a_SHARED_PHYSICS.md': {sha256(DIRECTIVE_PATH)!r},",
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
        restored = reconstructed[name]["value"]
        status, operands = total_compare(live_record["value"], restored)
        residual = sp.Integer(0) if status == PROVED_EQUAL else object_difference(live_record["value"], restored)
        rows.append(sp.Tuple(Str(name), Str(str(live_record["display"])), Str(status), casify(residual), casify(operands)))
        if status != PROVED_EQUAL:
            failures.append(name)
    emit("EXPORT_ROUNDTRIP", sp.Tuple(*rows), local=True)
    if failures:
        raise RuntimeError("export reconstruction was not proved equal for " + ", ".join(failures))


def remove_stale_export() -> None:
    if EXPORT_PATH.exists():
        EXPORT_PATH.unlink()


def run_primary(name: str, function) -> list[str]:
    global CURRENT_TASK
    CURRENT_TASK = name
    try:
        function()
        return [name]
    except Exception as exc:
        message = f"{name}: {type(exc).__name__}: {exc}"
        OPERATIONAL_FAILURES.append(message)
        return []
    finally:
        CURRENT_TASK = None


def run_control(name: str, function) -> list[str]:
    global CURRENT_TASK
    CURRENT_TASK = name
    try:
        function()
        return [name]
    except Exception as exc:
        message = f"{name}: {type(exc).__name__}: {exc}"
        OPERATIONAL_FAILURES.append(message)
        return []
    finally:
        CURRENT_TASK = None


def run() -> None:
    start = time.monotonic()
    completed: list[str] = []
    attempted: list[str] = []
    emit_supplied_objects()

    attempted.append("T0"); completed += run_primary("T0", task_t0)
    attempted.append("TA"); completed += run_primary("TA", task_ta)
    attempted.append("TA_PRIME"); completed += run_primary("TA_PRIME", task_ta_prime)
    attempted.append("TA_DOUBLE"); completed += run_primary("TA_DOUBLE", task_ta_double)
    attempted.append("TB"); completed += run_primary("TB", task_tb)
    attempted.append("TC"); completed += run_primary("TC", task_tc)
    attempted.append("TC_PRIME"); completed += run_primary("TC_PRIME", task_tc_prime)
    attempted.append("TD"); completed += run_primary("TD", task_td)
    attempted.append("TE"); completed += run_primary("TE", task_te)
    attempted.append("TF"); completed += run_primary("TF", task_tf)
    attempted.append("TG"); completed += run_primary("TG", task_tg)
    attempted.append("TH"); completed += run_primary("TH", task_th)
    attempted.append("TI"); completed += run_primary("TI", task_ti)

    attempted.append("REP_INVARIANCE"); completed += run_control("REP_INVARIANCE", task_rep_invariance)
    attempted.append("INDEPENDENCE"); completed += run_control("INDEPENDENCE", task_independence)
    attempted.append("FORM"); completed += run_control("FORM", task_form_control)
    attempted.append("UNIFORM"); completed += run_control("UNIFORM", task_uniform_control)
    attempted.append("HOMOGENEITY"); completed += run_control("HOMOGENEITY", task_homogeneity)
    primary_complete = all(task in completed for task in PRIMARY_TASKS)
    f9c_report = sp.Tuple()
    if primary_complete:
        try:
            ledger, diagnostics = merged_export()
            emit("EXPORT_CANDIDATE_KEY_OPERANDS", sp.Tuple(diagnostics["candidate_keys"], diagnostics["imported_matching_keys"], diagnostics["routing"]), local=True)
            f9c_report = diagnostics["f9c"]
            emit("EXPORT_F9C_WRITES", f9c_report, local=True)
            publish_export(ledger)
        except Exception as exc:
            OPERATIONAL_FAILURES.append(f"EXPORT: {type(exc).__name__}: {exc}")
            issue(f"export operational failure: {type(exc).__name__}: {exc}")
            remove_stale_export()
    else:
        remove_stale_export()
        issue("S11c_a_exports.py not published because a primary task did not complete")

    emit("OPERATIONAL_EXCEPTIONS", sp.Tuple(*(Str(item) for item in OPERATIONAL_FAILURES)), local=True)

    skipped = tuple(task for task in (*PRIMARY_TASKS, *CONTROL_TASKS) if task not in completed)
    run_record = sp.Tuple(*(Str(task) for task in completed))
    skipped_record = sp.Tuple(*(Str(task) for task in skipped))
    emit("RUN_TASKS", run_record, local=True)
    emit("SKIPPED_TASKS", skipped_record, local=True)

    local_names = list(dict.fromkeys((*EMITTER.local_tags, "PY_S11CA_LOCAL_TAGS", "PY_S11CA_LOCAL_SECTION8_REPORT")))
    emit("TAGS", sp.Tuple(*(Str(name) for name in local_names)), local=True)
    runtime = sp.Float(time.monotonic() - start, 8)
    export_lines = len(EXPORT_PATH.read_text(encoding="utf-8").splitlines()) if EXPORT_PATH.exists() else 0
    report = sp.Tuple(
        sp.Tuple(Str("FILES_WRITTEN"), sp.Tuple(Str(str(SCRIPT_PATH)), Str(str(EXPORT_PATH)) if EXPORT_PATH.exists() else Str("EXPORT_NOT_PUBLISHED"))),
        sp.Tuple(Str("SCRIPT_LINES"), sp.Integer(len(SCRIPT_PATH.read_text(encoding="utf-8").splitlines()))),
        sp.Tuple(Str("EXPORT_LINES"), sp.Integer(export_lines)),
        sp.Tuple(Str("EMITTED_TAGS"), sp.Integer(EMITTER.count + 1)),
        sp.Tuple(Str("RUN_TASKS"), run_record),
        sp.Tuple(Str("SKIPPED_TASKS"), skipped_record),
        sp.Tuple(Str("RUNTIME_SECONDS"), runtime),
        sp.Tuple(Str("TAG_NAMES"), sp.Tuple(*(Str(name) for name in (*EMITTER.values.keys(), "PY_S11CA_LOCAL_SECTION8_REPORT")))),
        sp.Tuple(Str("ISSUES_OR_NONCOMPUTABLE"), sp.Tuple(*(Str(item) for item in (*ISSUES, *OPERATIONAL_FAILURES)))),
        sp.Tuple(Str("F9C_WRITES"), f9c_report),
        sp.Tuple(Str("SUPPLIED_UNFALSIFIABLE"), Str("SECTIONS_1_TO_3_AND_ADMISSIBILITY_PREMISE")),
    )
    emit("SECTION8_REPORT", report, local=True)
    if OPERATIONAL_FAILURES:
        raise RuntimeError("; ".join(OPERATIONAL_FAILURES))


if __name__ == "__main__":
    run()
