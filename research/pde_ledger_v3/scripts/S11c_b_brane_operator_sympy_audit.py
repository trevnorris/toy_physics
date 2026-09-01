#!/usr/bin/env python3
"""SymPy audit for the S11c-b variable-coefficient brane operator.

The sole physics input is ``directives/S11c_b_SHARED_PHYSICS.md``.  The
incoming accumulated ledger supplies object identity and the S11c-a
shape-derivative substrate.  This program streams CAS objects and, after all
four primary tasks complete, writes the accumulated S11c-b export module.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from itertools import combinations_with_replacement, product
from pathlib import Path
import atexit
import hashlib
import multiprocessing as mp
import os
import sys
import time

import sympy as sp
from sympy.core.relational import Relational
from sympy.core.symbol import Str
from sympy.logic.boolalg import Boolean
from sympy.printing.str import StrPrinter


SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
DIRECTIVE_PATH = SCRIPT_DIR.parent / "directives" / "S11c_b_SHARED_PHYSICS.md"
EXPORT_PATH = SCRIPT_DIR / "S11c_b_exports.py"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from S11c_a_exports import LEDGER as INCOMING_LEDGER  # noqa: E402


CLASS_TAGS = {"KNOB", "STRUCTURAL", "COORDINATE", "CONTROL", "PREMISE", "DERIVED"}
PRIMARY_TASKS = ("ENERGY_BASIS", "SLAB_OPERATOR", "COUPLING_KERNEL", "ADMISSIBILITY")
CONTROL_TASKS = (
    "PROJECTION_EQUIVALENCE",
    "REP_INVARIANCE",
    "INDEPENDENCE",
    "FORM",
    "HESSIAN_FREEZE",
    "COEFFICIENT",
    "JET_ZERO",
    "TOWER_DEPTH",
    "UNIFORM",
    "HOMOGENEITY",
)
BRANCHES = ("LAB_HELD", "MATERIAL_ADVECTED")
DENSITY_REPS = ("RHO4_CONSTANT", "RHOBR_CONSTANT")
FACES = (1, -1)
DIRECTIONS = range(3)
HESSIAN_FREEZE_DEPTH = 1
STRONG_ROW_JET_DEPTH = 2
COUPLING_JET_DEPTH = 3
DEPTH_CONTROL_JET_DEPTH = 4

PROVED_EQUAL = "PROVED_EQUAL"
PROVED_DIFFERENT = "PROVED_DIFFERENT"
UNDECIDED = "UNDECIDED"

DECLARED_SYMBOLS: dict[sp.Symbol, dict[str, str]] = {}
SYMBOL_DIMENSIONS: dict[sp.Symbol, sp.ImmutableMatrix] = {}
ISSUES: list[str] = []
OPERATIONAL_FAILURES: list[str] = []
CURRENT_TASK: str | None = None
RAW_PRIMARY: dict[str, object] = {}
DIMENSION_PRIMARY: dict[str, object] = {}
ENERGY_CACHE: dict[tuple[object, ...], "EnergyBuild"] = {}
OPERATOR_CACHE: dict[tuple[object, ...], tuple[object, object, object]] = {}
KERNEL_CACHE: dict[tuple[object, ...], tuple[object, object]] = {}
KERNEL_BLOCK_CACHE: dict[tuple[object, ...], tuple[object, object]] = {}
KERNEL_ORIGIN_CACHE: dict[tuple[object, ...], object] = {}
FIRST_SHAPE_SERIES_CACHE: dict[sp.Basic, sp.Expr] = {}
PROJECTION_POOL: object | None = None
TASK_TIMINGS: dict[str, float] = {}


DIM_ZERO = sp.ImmutableMatrix([0, 0, 0])
DIM_L = sp.ImmutableMatrix([1, 0, 0])
DIM_T = sp.ImmutableMatrix([0, 1, 0])
DIM_M = sp.ImmutableMatrix([0, 0, 1])


def dim_mul(*dimensions: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(sum((sp.Matrix(d) for d in dimensions), sp.zeros(3, 1)))


def dim_div(left: sp.MatrixBase, right: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(sp.Matrix(left) - sp.Matrix(right))


DIM_ENERGY = dim_mul(DIM_M, -DIM_L, -2 * DIM_T)
DIM_BODY_U = dim_div(DIM_ENERGY, DIM_L)
DIM_PRESSURE = dim_mul(DIM_M, -2 * DIM_L, -2 * DIM_T)
DIM_RHO4 = dim_mul(DIM_M, -4 * DIM_L)
DIM_RHOBR = dim_mul(DIM_M, -3 * DIM_L)
DIM_VELOCITY = dim_div(DIM_L, DIM_T)


def register_symbol(
    value: sp.Symbol,
    class_tag: str,
    description: str,
    dimension: sp.MatrixBase | None = None,
) -> sp.Symbol:
    if class_tag not in CLASS_TAGS:
        raise ValueError(f"unknown class {class_tag!r}")
    if not isinstance(value, sp.Symbol):
        raise TypeError(f"declared object is not a Symbol: {value!r}")
    DECLARED_SYMBOLS[value] = {"class": class_tag, "description": description}
    if dimension is not None:
        SYMBOL_DIMENSIONS[value] = sp.ImmutableMatrix(dimension)
    return value


def symbol(
    name: str,
    class_tag: str,
    description: str,
    dimension: sp.MatrixBase | None = None,
    **assumptions: object,
) -> sp.Symbol:
    return register_symbol(sp.Symbol(name, **assumptions), class_tag, description, dimension)


def inherited_symbol(
    name: str,
    class_tag: str,
    description: str,
    dimension: sp.MatrixBase | None = None,
) -> sp.Symbol:
    value = INCOMING_LEDGER[name]["value"]
    return register_symbol(value, class_tag, description, dimension)


# Exact inherited identities required by the chain directive.
c_s0 = inherited_symbol("c_s0", "KNOB", "bulk sound speed", DIM_VELOCITY)
mu_R = inherited_symbol("mu_R", "KNOB", "uniform curl modulus", DIM_ENERGY)
rho_br = inherited_symbol("rho_br", "KNOB", "uniform integrated brane density", DIM_RHOBR)
W0 = inherited_symbol("W_0", "KNOB", "uniform reference thickness", DIM_L)
e_W = inherited_symbol("e_W", "COORDINATE", "uniform-reference thickness fraction", DIM_ZERO)
rho_m = inherited_symbol("rho_m", "KNOB", "bulk mass density", DIM_RHO4)
v_dr = inherited_symbol(
    "v_bulk_normal_0", "KNOB", "inert bulk-normal drain scope-limit speed", DIM_VELOCITY
)

# Other supplied constants needed by the independently constructed basis and balance laws.
B_rho_3 = inherited_symbol("B_rho_3", "KNOB", "uniform integrated compression modulus", DIM_ENERGY)
C = inherited_symbol("C", "KNOB", "densification-thickness coefficient", dim_div(DIM_ENERGY, DIM_L))
k_W = inherited_symbol("k_W", "KNOB", "thickness restoring coefficient", dim_div(DIM_ENERGY, 2 * DIM_L))
kappa_W = inherited_symbol(
    "kappa_W", "KNOB", "thickness-gradient coefficient", dim_div(DIM_ENERGY, 2 * DIM_L)
)
mu_W = inherited_symbol("mu_W", "KNOB", "thickness inertia", dim_mul(DIM_M, -3 * DIM_L))
mu_S = inherited_symbol("mu_S", "KNOB", "symmetric-gradient modulus", DIM_ENERGY)
G_theta_u = inherited_symbol("G_theta_u", "KNOB", "densification-divergence coefficient", DIM_ENERGY)
G_W_u = inherited_symbol("G_W_u", "KNOB", "thickness-divergence coefficient", DIM_ENERGY)
kappa_theta = inherited_symbol(
    "kappa_theta", "KNOB", "densification-gradient coefficient", dim_mul(DIM_ENERGY, 2 * DIM_L)
)
kappa_theta_W = inherited_symbol(
    "kappa_theta_W", "KNOB", "mixed scalar-gradient coefficient", dim_mul(DIM_ENERGY, 2 * DIM_L)
)

# The S11c-a profile rows are bound, never re-originated.
W_bg = inherited_symbol("W_bg", "DERIVED", "spatially varying background thickness", DIM_L)
mu_R_bg = inherited_symbol("mu_R_bg", "DERIVED", "spatially varying curl modulus", DIM_ENERGY)
w1_profile = inherited_symbol("w1_profile", "KNOB", "dimensionless thickness profile", DIM_ZERO)
m1_profile = inherited_symbol("m1_profile", "KNOB", "dimensionless modulus profile", DIM_ZERO)
L_W = inherited_symbol("L_W", "KNOB", "independent profile length", DIM_L)
sigma_W = inherited_symbol("sigma_W", "DERIVED", "first-background-jet bookkeeper", DIM_ZERO)
eta_bg = inherited_symbol("eta_bg", "KNOB", "zero-jet contrast bookkeeper", DIM_ZERO)
e_W_bg = inherited_symbol("e_W_bg", "DERIVED", "local-background thickness fraction", DIM_ZERO)

grad_W = tuple(
    inherited_symbol(f"W_bg_d{i}", "DERIVED", f"background thickness first jet {i}", DIM_ZERO)
    for i in range(1, 4)
)
grad_mu = tuple(
    inherited_symbol(
        f"mu_R_bg_d{i}", "DERIVED", f"background modulus first jet {i}", dim_div(DIM_ENERGY, DIM_L)
    )
    for i in range(1, 4)
)
w1_grad = tuple(
    inherited_symbol(f"w1_profile_d{i}", "KNOB", f"thickness-profile first jet {i}", DIM_ZERO)
    for i in range(1, 4)
)
m1_grad = tuple(
    inherited_symbol(f"m1_profile_d{i}", "KNOB", f"modulus-profile first jet {i}", DIM_ZERO)
    for i in range(1, 4)
)

rho4_rho4 = inherited_symbol(
    "rho_4D_bg_rho4_constant", "DERIVED", "RHO4-CONSTANT background four-density", DIM_RHO4
)
rhobr_rho4 = inherited_symbol(
    "rho_br_bg_rho4_constant", "DERIVED", "RHO4-CONSTANT integrated density", DIM_RHOBR
)
rho4_rhobr = inherited_symbol(
    "rho_4D_bg_rhobr_constant", "DERIVED", "RHOBR-CONSTANT background four-density", DIM_RHO4
)
rhobr_rhobr = inherited_symbol(
    "rho_br_bg_rhobr_constant", "DERIVED", "RHOBR-CONSTANT integrated density", DIM_RHOBR
)

# Wave coordinates use the S11c-a identities where available.
theta = inherited_symbol("theta", "COORDINATE", "Eulerian densification", DIM_ZERO)
u = tuple(
    inherited_symbol(f"u_{a}", "COORDINATE", f"displacement component {a}", DIM_L)
    for a in range(1, 4)
)
grad_u = tuple(
    tuple(
        inherited_symbol(f"u_{a}_d{i}", "COORDINATE", f"displacement first jet {a},{i}", DIM_ZERO)
        for i in range(1, 4)
    )
    for a in range(1, 4)
)
grad_theta = tuple(
    inherited_symbol(f"grad_theta_{i}", "COORDINATE", f"densification first jet {i}", -DIM_L)
    for i in range(1, 4)
)
grad_e = tuple(
    inherited_symbol(f"e_W_d{i}", "COORDINATE", f"thickness-fraction first jet {i}", -DIM_L)
    for i in range(1, 4)
)
u_t = tuple(
    inherited_symbol(f"u_{a}_t", "COORDINATE", f"displacement velocity {a}", DIM_VELOCITY)
    for a in range(1, 4)
)
e_t = inherited_symbol("e_W_t", "COORDINATE", "thickness-fraction time derivative", -DIM_T)
epsilon = inherited_symbol("epsilon_shape", "CONTROL", "wave-amplitude bookkeeper", DIM_ZERO)


def symmetric_jet(prefix: str, rank: int, dimension: sp.MatrixBase) -> dict[tuple[int, ...], sp.Symbol]:
    jets: dict[tuple[int, ...], sp.Symbol] = {}
    if rank not in (2, 3, 4):
        raise ValueError(rank)
    indices = combinations_with_replacement(DIRECTIONS, rank)
    for index in indices:
        suffix = "d" + "d".join(str(i + 1) for i in index)
        jets[index] = symbol(
            f"{prefix}_{suffix}", "COORDINATE", f"spatial jet {prefix} {index}", dimension
        )
    return jets


u_second = tuple(symmetric_jet(f"u_{a + 1}", 2, -DIM_L) for a in DIRECTIONS)
u_third = tuple(symmetric_jet(f"u_{a + 1}", 3, -2 * DIM_L) for a in DIRECTIONS)
theta_second = symmetric_jet("theta", 2, -2 * DIM_L)
theta_third = symmetric_jet("theta", 3, -3 * DIM_L)
e_second = symmetric_jet("e_W", 2, -2 * DIM_L)
e_third = symmetric_jet("e_W", 3, -3 * DIM_L)
u_tt = tuple(
    symbol(f"u_{a + 1}_tt", "COORDINATE", f"displacement acceleration {a + 1}", dim_div(DIM_L, 2 * DIM_T))
    for a in DIRECTIONS
)
e_tt = symbol("e_W_tt", "COORDINATE", "thickness-fraction second time derivative", -2 * DIM_T)

# Local sector probes.  These are differential labels, not global projections.
uT = tuple(symbol(f"u_T_{a + 1}", "COORDINATE", f"local curl-sector probe {a + 1}", DIM_L) for a in DIRECTIONS)
GT = tuple(
    tuple(symbol(f"u_T_{a + 1}_d{i + 1}", "COORDINATE", f"curl-sector first jet {a + 1},{i + 1}", DIM_ZERO) for i in DIRECTIONS)
    for a in DIRECTIONS
)
uL = tuple(symbol(f"u_L_{a + 1}", "COORDINATE", f"local divergence-sector probe {a + 1}", DIM_L) for a in DIRECTIONS)
GL = tuple(
    tuple(symbol(f"u_L_{a + 1}_d{i + 1}", "COORDINATE", f"divergence-sector first jet {a + 1},{i + 1}", DIM_ZERO) for i in DIRECTIONS)
    for a in DIRECTIONS
)
theta_probe = symbol("theta_probe", "COORDINATE", "densification-sector probe", DIM_ZERO)
theta_probe_grad = tuple(symbol(f"theta_probe_d{i + 1}", "COORDINATE", f"densification probe jet {i + 1}", -DIM_L) for i in DIRECTIONS)
e_probe = symbol("e_W_probe", "COORDINATE", "thickness-sector probe", DIM_ZERO)
e_probe_grad = tuple(symbol(f"e_W_probe_d{i + 1}", "COORDINATE", f"thickness probe jet {i + 1}", -DIM_L) for i in DIRECTIONS)

uT_t = tuple(symbol(f"u_T_{a + 1}_t", "COORDINATE", f"curl-sector velocity probe {a + 1}", DIM_VELOCITY) for a in DIRECTIONS)
uL_t = tuple(symbol(f"u_L_{a + 1}_t", "COORDINATE", f"divergence-sector velocity probe {a + 1}", DIM_VELOCITY) for a in DIRECTIONS)
e_probe_t = symbol("e_W_probe_t", "COORDINATE", "thickness probe time derivative", -DIM_T)

# The support bundle is supplied; its components are operands of the computed package.
support_body_u = tuple(symbol(f"f_hold_u_{a + 1}_0", "PREMISE", f"held background u-force {a + 1}", DIM_BODY_U) for a in DIRECTIONS)
support_body_theta = symbol("f_hold_theta_0", "PREMISE", "held background theta-force", DIM_ENERGY)
support_body_e = symbol("f_hold_e_W_0", "PREMISE", "held background thickness-force", DIM_ENERGY)
support_face = {
    face: tuple(
        symbol(
            f"t_hold_{'plus' if face == 1 else 'minus'}_0_{a + 1}",
            "PREMISE",
            f"held face traction {face},{a + 1}",
            DIM_PRESSURE,
        )
        for a in range(4)
    )
    for face in FACES
}

# Each potential spurion invariant has a fresh coefficient; only coefficients
# selected by the quotient computation enter the primary objects and exports.
NEW_COEFFICIENTS = {
    source: tuple(
        symbol(
            f"gamma_s11cb_{source.lower()}_{index:02d}",
            "KNOB",
            f"coefficient for a constructed {source} first-jet invariant {index}",
        )
        for index in range(1, 32)
    )
    for source in ("W_BG", "MU_R_BG")
}


SHAPE_SUBSTRATE_KEYS = (
    "background_density_map",
    "face_normal",
    "conormal_deriv",
    "face_measure_shape_deriv",
    "face_velocity",
    "relative_flux",
    "kinematic_balance",
    "traction",
    "virtual_work_shape_deriv",
    "face_shift",
    "projection_shape_deriv",
    "projection_term_origins",
    "projection_static_operand",
    "projection_dynamic_operand",
    "projection_residual",
    "virtual_constraint",
    "evolution_mass_balance",
    "evolution_term_origins",
    "closure_shape_deriv",
)
S11CA_SUBSTRATE = {name: INCOMING_LEDGER[name]["value"] for name in SHAPE_SUBSTRATE_KEYS}


def bind_additional_inherited(
    name: str,
    class_tag: str,
    description: str,
    dimension: sp.MatrixBase,
) -> None:
    value = INCOMING_LEDGER[name]["value"]
    if value not in DECLARED_SYMBOLS:
        register_symbol(value, class_tag, description, dimension)


DIM_AFFINITY = dim_mul(2 * DIM_L, -2 * DIM_T)
DIM_FLUX = dim_mul(DIM_RHO4, DIM_VELOCITY)
for name, dimension in (
    ("Lambda_A_0", dim_div(DIM_FLUX, DIM_AFFINITY)),
    ("Lambda_V_0", dim_div(DIM_FLUX, DIM_VELOCITY)),
    ("Lambda_X_0", dim_div(DIM_PRESSURE, DIM_AFFINITY)),
    ("tau_A", DIM_T),
    ("tau_V", DIM_T),
    ("tau_X", DIM_T),
):
    bind_additional_inherited(name, "KNOB", f"inherited face-law knob {name}", dimension)

for i in range(1, 4):
    for j in range(i, 4):
        bind_additional_inherited(
            f"w1_profile_d{i}d{j}",
            "KNOB",
            f"dimensionless thickness-profile second jet {i},{j}",
            DIM_ZERO,
        )


def symmetric_profile_jet(
    prefix: str,
    rank: int,
    description: str,
    *,
    inherited: bool = False,
) -> dict[tuple[int, ...], sp.Symbol]:
    table = {}
    for index in combinations_with_replacement(DIRECTIONS, rank):
        suffix = "d".join(str(i + 1) for i in index)
        name = f"{prefix}_d{suffix}"
        if inherited:
            value = INCOMING_LEDGER[name]["value"]
            if value not in DECLARED_SYMBOLS:
                register_symbol(value, "KNOB", f"{description} {index}", DIM_ZERO)
        else:
            value = symbol(name, "KNOB", f"{description} {index}", DIM_ZERO)
        table[index] = value
    return table


W_PROFILE_JETS: dict[int, Mapping[tuple[int, ...], sp.Symbol]] = {
    1: {(i,): w1_grad[i] for i in DIRECTIONS},
    2: symmetric_profile_jet(
        "w1_profile", 2, "dimensionless thickness-profile second jet", inherited=True
    ),
    3: symmetric_profile_jet(
        "w1_profile", 3, "dimensionless thickness-profile third jet"
    ),
    4: symmetric_profile_jet(
        "w1_profile", 4, "dimensionless thickness-profile fourth jet depth control"
    ),
}
M_PROFILE_JETS: dict[int, Mapping[tuple[int, ...], sp.Symbol]] = {
    1: {(i,): m1_grad[i] for i in DIRECTIONS},
    2: symmetric_profile_jet(
        "m1_profile", 2, "dimensionless modulus-profile second jet"
    ),
    3: symmetric_profile_jet(
        "m1_profile", 3, "dimensionless modulus-profile third jet"
    ),
    4: symmetric_profile_jet(
        "m1_profile", 4, "dimensionless modulus-profile fourth jet depth control"
    ),
}

for name in ("mu_theta_L", "mu_theta_M"):
    bind_additional_inherited(
        name, "PREMISE", f"S11c-a reserved variable-coefficient operand {name}", DIM_ENERGY
    )

for face_name in ("plus", "minus"):
    bind_additional_inherited(
        f"delta_p_{face_name}", "COORDINATE", f"{face_name} face pressure perturbation", DIM_PRESSURE
    )
    bind_additional_inherited(
        f"d_w_delta_p_{face_name}",
        "COORDINATE",
        f"{face_name} face pressure normal jet",
        dim_div(DIM_PRESSURE, DIM_L),
    )
    bind_additional_inherited(
        f"delta_rho_4D_face_{face_name}",
        "COORDINATE",
        f"{face_name} face density perturbation",
        DIM_RHO4,
    )
    bind_additional_inherited(
        f"d_w_delta_rho_4D_face_{face_name}",
        "COORDINATE",
        f"{face_name} face density normal jet",
        dim_div(DIM_RHO4, DIM_L),
    )
    for component in range(1, 5):
        bind_additional_inherited(
            f"delta_v_bulk_{face_name}_{component}",
            "COORDINATE",
            f"{face_name} bulk velocity perturbation {component}",
            DIM_VELOCITY,
        )
        bind_additional_inherited(
            f"d_w_delta_v_bulk_{face_name}_{component}",
            "COORDINATE",
            f"{face_name} bulk velocity normal jet {component}",
            dim_div(DIM_VELOCITY, DIM_L),
        )
        bind_additional_inherited(
            f"d_w_delta_j_bulk_{face_name}_{component}",
            "COORDINATE",
            f"{face_name} bulk-current normal jet {component}",
            dim_div(DIM_FLUX, DIM_L),
        )

for component in range(1, 5):
    bind_additional_inherited(
        f"delta_j_bulk_{component}",
        "COORDINATE",
        f"bulk-current perturbation {component}",
        DIM_FLUX,
    )
    bind_additional_inherited(
        f"trace_grad_f_{component}",
        "COORDINATE",
        f"generic traced-field derivative {component}",
        -DIM_L,
    )
    bind_additional_inherited(
        f"d_w_trace_grad_f_{component}",
        "COORDINATE",
        f"generic traced-field normal second derivative {component}",
        -2 * DIM_L,
    )

bind_additional_inherited(
    "delta_rho_4D_bulk_t",
    "COORDINATE",
    "bulk-density perturbation time derivative",
    dim_div(DIM_RHO4, DIM_T),
)
for name in ("delta_v_theta", "delta_v_e_W"):
    bind_additional_inherited(name, "COORDINATE", f"virtual scalar variation {name}", DIM_ZERO)
for component in range(1, 4):
    bind_additional_inherited(
        f"delta_v_u_{component}",
        "COORDINATE",
        f"virtual displacement component {component}",
        DIM_L,
    )
    bind_additional_inherited(
        f"delta_v_u_{component}_d{component}",
        "COORDINATE",
        f"virtual displacement diagonal jet {component}",
        DIM_ZERO,
    )
    bind_additional_inherited(
        f"u_{component}_t_d{component}",
        "COORDINATE",
        f"velocity diagonal jet {component}",
        -DIM_T,
    )

bind_additional_inherited("theta_t", "COORDINATE", "densification time derivative", -DIM_T)
bind_additional_inherited("zeta_c", "COORDINATE", "face-centre displacement", DIM_L)
bind_additional_inherited("zeta_c_t", "COORDINATE", "face-centre velocity", DIM_VELOCITY)
bind_additional_inherited(
    "delta_v_zeta_c", "COORDINATE", "virtual face-centre displacement", DIM_L
)
for i in range(1, 4):
    bind_additional_inherited(
        f"zeta_c_d{i}", "COORDINATE", f"face-centre displacement jet {i}", DIM_ZERO
    )


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
        self.local_tags: list[str] = []
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
        print(f"{tag}: {render(payload)}", flush=True)
        self.count += 1
        self.values[tag] = payload
        if "_LOCAL_" in tag:
            self.local_tags.append(tag)
        if export_key is not None:
            self.export_candidates.append(
                CandidateRow(export_key, payload, class_tag, tag, description)
            )
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
    tag = f"PY_S11CB_{infix}{quantity}"
    return EMITTER.emit(
        tag,
        payload,
        export_key=key if export and not local else None,
        class_tag=class_tag,
    )


def issue(message: str) -> None:
    ISSUES.append(f"{CURRENT_TASK or 'GLOBAL'}: {message}")


def dot(left: tuple[sp.Expr, ...], right: tuple[sp.Expr, ...]) -> sp.Expr:
    return sp.Add(*(a * b for a, b in zip(left, right)))


def trace(matrix: tuple[tuple[sp.Expr, ...], ...]) -> sp.Expr:
    return sp.Add(*(matrix[i][i] for i in DIRECTIONS))


def sorted_index(*indices: int) -> tuple[int, ...]:
    return tuple(sorted(indices))


DERIVATIVE_MAP: dict[int, dict[sp.Symbol, sp.Expr]] = {i: {} for i in DIRECTIONS}


def add_field_derivatives(
    field: sp.Symbol,
    first: tuple[sp.Symbol, ...],
    second: Mapping[tuple[int, ...], sp.Symbol],
    third: Mapping[tuple[int, ...], sp.Symbol],
) -> None:
    for i in DIRECTIONS:
        DERIVATIVE_MAP[i][field] = first[i]
        for j in DIRECTIONS:
            DERIVATIVE_MAP[i][first[j]] = second[sorted_index(i, j)]
        for j in DIRECTIONS:
            for k in range(j, 3):
                second_symbol = second[sorted_index(j, k)]
                DERIVATIVE_MAP[i][second_symbol] = third[sorted_index(i, j, k)]


for a in DIRECTIONS:
    add_field_derivatives(u[a], grad_u[a], u_second[a], u_third[a])
add_field_derivatives(theta, grad_theta, theta_second, theta_third)
add_field_derivatives(e_W, grad_e, e_second, e_third)
for i in DIRECTIONS:
    DERIVATIVE_MAP[i][W_bg] = grad_W[i]
    DERIVATIVE_MAP[i][mu_R_bg] = grad_mu[i]


BACKGROUND_FIRST_JETS = {"W_BG": grad_W, "MU_R_BG": grad_mu}
BACKGROUND_PROFILE_JETS = {"W_BG": W_PROFILE_JETS, "MU_R_BG": M_PROFILE_JETS}
BACKGROUND_ZERO_JETS = {"W_BG": W_bg, "MU_R_BG": mu_R_bg}


def background_jet_expression(
    source: str,
    order: int,
    index: tuple[int, ...],
) -> sp.Expr:
    if order == 1:
        return BACKGROUND_FIRST_JETS[source][index[0]]
    profile_atom = BACKGROUND_PROFILE_JETS[source][order][sorted_index(*index)]
    if source == "W_BG":
        return sigma_W * profile_atom / L_W ** (order - 1)
    if source == "MU_R_BG":
        return mu_R * sigma_W * profile_atom / (W0 * L_W ** (order - 1))
    raise ValueError(source)


def total_derivative(
    expression: sp.Expr,
    direction: int,
    *,
    background_depth: int = COUPLING_JET_DEPTH,
    background_first: Mapping[str, tuple[sp.Expr, ...]] | None = None,
    zero_sources: frozenset[str] = frozenset(),
) -> sp.Expr:
    if background_depth < HESSIAN_FREEZE_DEPTH or background_depth > DEPTH_CONTROL_JET_DEPTH:
        raise ValueError(background_depth)
    derivative_map = dict(DERIVATIVE_MAP[direction])
    first_jets = BACKGROUND_FIRST_JETS if background_first is None else background_first
    for source in BACKGROUND_ZERO_JETS:
        zero_jet = BACKGROUND_ZERO_JETS[source]
        first = BACKGROUND_FIRST_JETS[source]
        if source in zero_sources:
            derivative_map.pop(zero_jet, None)
            for atom in first:
                derivative_map.pop(atom, None)
            continue
        derivative_map[zero_jet] = first_jets[source][direction]
        if background_depth >= 2:
            for jet_direction in DIRECTIONS:
                derivative_map[first[jet_direction]] = background_jet_expression(
                    source, 2, sorted_index(direction, jet_direction)
                )
        else:
            for atom in first:
                derivative_map.pop(atom, None)
        for order in range(2, background_depth):
            for index, atom in BACKGROUND_PROFILE_JETS[source][order].items():
                derivative_map[atom] = (
                    BACKGROUND_PROFILE_JETS[source][order + 1][
                        sorted_index(direction, *index)
                    ]
                    / L_W
                )
    result = sp.Integer(0)
    for atom, derivative in derivative_map.items():
        if atom in expression.free_symbols:
            result += sp.diff(expression, atom) * derivative
    return sp.expand(result)


def frozen_derivative(expression: sp.Expr, direction: int) -> sp.Expr:
    result = sp.Integer(0)
    for atom, derivative in DERIVATIVE_MAP[direction].items():
        if atom in expression.free_symbols:
            result += sp.diff(expression, atom) * derivative
    return sp.expand(result)


def dx(expression: sp.Expr, direction: int) -> sp.Expr:
    return total_derivative(expression, direction)


def euler_derivative(
    density: sp.Expr,
    field: sp.Symbol,
    first_jets: tuple[sp.Symbol, ...],
    background_depth: int = COUPLING_JET_DEPTH,
) -> sp.Expr:
    return sp.expand(
        sp.diff(density, field)
        - sp.Add(
            *(
                total_derivative(
                    sp.diff(density, first_jets[i]),
                    i,
                    background_depth=background_depth,
                )
                for i in DIRECTIONS
            )
        )
    )


def curl(
    vector: tuple[sp.Expr, ...],
    background_depth: int = COUPLING_JET_DEPTH,
) -> tuple[sp.Expr, ...]:
    return (
        sp.expand(
            total_derivative(vector[2], 1, background_depth=background_depth)
            - total_derivative(vector[1], 2, background_depth=background_depth)
        ),
        sp.expand(
            total_derivative(vector[0], 2, background_depth=background_depth)
            - total_derivative(vector[2], 0, background_depth=background_depth)
        ),
        sp.expand(
            total_derivative(vector[1], 0, background_depth=background_depth)
            - total_derivative(vector[0], 1, background_depth=background_depth)
        ),
    )


def divergence(
    vector: tuple[sp.Expr, ...],
    background_depth: int = COUPLING_JET_DEPTH,
) -> sp.Expr:
    return sp.expand(
        sp.Add(
            *(
                total_derivative(vector[i], i, background_depth=background_depth)
                for i in DIRECTIONS
            )
        )
    )


def profile_definitions() -> sp.Tuple:
    return sp.Tuple(
        sp.Eq(W_bg, W0 * (1 + eta_bg * w1_profile), evaluate=False),
        sp.Eq(mu_R_bg, mu_R * (1 + eta_bg * m1_profile), evaluate=False),
        sp.Eq(e_W_bg, W0 * e_W / W_bg, evaluate=False),
        sp.Tuple(*(sp.Eq(grad_W[i], sigma_W * w1_grad[i], evaluate=False) for i in DIRECTIONS)),
        sp.Tuple(
            *(
                sp.Eq(grad_mu[i], mu_R * sigma_W * m1_grad[i] / W0, evaluate=False)
                for i in DIRECTIONS
            )
        ),
    )


PROFILE_GRADE_SUBS = {
    W_bg: W0 * (1 + eta_bg * w1_profile),
    mu_R_bg: mu_R * (1 + eta_bg * m1_profile),
    **{grad_W[i]: sigma_W * w1_grad[i] for i in DIRECTIONS},
    **{grad_mu[i]: mu_R * sigma_W * m1_grad[i] / W0 for i in DIRECTIONS},
    rho4_rho4: rho_br / W0,
    rhobr_rho4: rho_br * (1 + eta_bg * w1_profile),
    rho4_rhobr: rho_br / (W0 * (1 + eta_bg * w1_profile)),
    rhobr_rhobr: rho_br,
}


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
    if isinstance(value, Str):
        return value
    if isinstance(value, sp.Basic):
        return scalar_function(value)
    return value


WAVE_SYMBOLS = {
    theta,
    e_W,
    *u,
    *grad_theta,
    *grad_e,
    *(item for row in grad_u for item in row),
    *(item for row in u_second for item in row.values()),
    *theta_second.values(),
    *e_second.values(),
    *u_t,
    e_t,
    *u_tt,
    e_tt,
}


def first_shape_series_reference(expression: sp.Expr) -> sp.Expr:
    expression = sp.sympify(expression).subs(PROFILE_GRADE_SUBS, simultaneous=True)
    try:
        expression = sp.series(expression, eta_bg, 0, 2).removeO()
    except (NotImplementedError, ValueError, TypeError):
        pass
    expression = sp.expand(expression)
    kept = []
    for term in sp.Add.make_args(expression):
        powers = term.as_powers_dict()
        if powers.get(eta_bg, 0) <= 1 and powers.get(sigma_W, 0) <= 1:
            kept.append(term)
    return sp.Add(*kept)


def first_shape_series_fast(expression: sp.Expr) -> sp.Expr:
    # Lever C.  The retained grade keeps only the eta_bg^0 and eta_bg^1 parts, i.e.
    # series(...,eta_bg,0,2).removeO() == f(0) + eta_bg*f'(0).  For these
    # analytic-at-0 profile expressions the direct first-order Taylor equals the
    # sp.series result exactly, while avoiding the heavyweight series machinery
    # (measured 5-8x+ on large scalars, byte-identical).  Any error falls back to
    # the sp.series reference so the failure/swallow behaviour is never altered.
    substituted = sp.sympify(expression).subs(PROFILE_GRADE_SUBS, simultaneous=True)
    # Integral/Derivative-bearing scalars: the direct-Taylor equivalence is not
    # numerically checkable and differentiation-under-integral is an edge case, so
    # use the exact sp.series reference on them (rare; keeps them byte-identical).
    if substituted.has(sp.Integral) or substituted.has(sp.Derivative):
        return first_shape_series_reference(expression)
    try:
        f0 = substituted.subs(eta_bg, 0)
        f1 = substituted.diff(eta_bg).subs(eta_bg, 0)
        truncated = f0 + eta_bg * f1
    except (NotImplementedError, ValueError, TypeError, RecursionError, AttributeError):
        return first_shape_series_reference(expression)
    truncated = sp.expand(truncated)
    kept = []
    for term in sp.Add.make_args(truncated):
        powers = term.as_powers_dict()
        if powers.get(eta_bg, 0) <= 1 and powers.get(sigma_W, 0) <= 1:
            kept.append(term)
    return sp.Add(*kept)


def first_shape_series(expression: sp.Expr) -> sp.Expr:
    try:
        return FIRST_SHAPE_SERIES_CACHE[expression]
    except KeyError:
        projected = first_shape_series_fast(expression)
        FIRST_SHAPE_SERIES_CACHE[expression] = projected
        return projected


def projection_worker_count() -> int:
    configured = os.environ.get("S11CB_PROJECTION_WORKERS")
    if configured is not None:
        workers = int(configured)
        if workers < 1:
            raise ValueError("S11CB_PROJECTION_WORKERS must be at least 1")
        return workers
    logical_cores = os.cpu_count() or 1
    return max(1, logical_cores // 2)


def object_scalars(value: object) -> tuple[sp.Basic, ...]:
    scalars: list[sp.Basic] = []

    def collect(expression: sp.Basic) -> sp.Basic:
        scalars.append(expression)
        return expression

    map_object(value, collect)
    return tuple(scalars)


def project_scalar_worker(expression: sp.Expr) -> sp.Expr:
    return first_shape_series(expression)


def projection_pool() -> object:
    global PROJECTION_POOL
    if PROJECTION_POOL is None:
        PROJECTION_POOL = mp.get_context("fork").Pool(projection_worker_count())
    return PROJECTION_POOL


def close_projection_pool() -> None:
    global PROJECTION_POOL
    if PROJECTION_POOL is not None:
        PROJECTION_POOL.close()
        PROJECTION_POOL.join()
        PROJECTION_POOL = None


atexit.register(close_projection_pool)


def project_scalars_ordered(
    scalars: tuple[sp.Basic, ...], *, allow_pool: bool
) -> tuple[sp.Expr, ...]:
    missing = tuple(
        dict.fromkeys(
            expression
            for expression in scalars
            if expression not in FIRST_SHAPE_SERIES_CACHE
        )
    )
    can_fork = (
        allow_pool
        and projection_worker_count() > 1
        and mp.current_process().name == "MainProcess"
    )
    if missing and can_fork:
        projected = projection_pool().map(project_scalar_worker, missing, chunksize=1)
        FIRST_SHAPE_SERIES_CACHE.update(zip(missing, projected))
    else:
        for expression in missing:
            first_shape_series(expression)
    return tuple(FIRST_SHAPE_SERIES_CACHE[expression] for expression in scalars)


def reassemble_projected_object(
    value: object, projected: tuple[sp.Expr, ...]
) -> object:
    iterator = iter(projected)
    result = map_object(value, lambda _expression: next(iterator))
    try:
        next(iterator)
    except StopIteration:
        return result
    raise RuntimeError("projection reassembly left unused scalar results")


def retained_grade_serial(value: object) -> object:
    return map_object(value, first_shape_series)


def retained_grade_parallel(value: object) -> object:
    scalars = object_scalars(value)
    projected = project_scalars_ordered(scalars, allow_pool=True)
    return reassemble_projected_object(value, projected)


def retained_grade(value: object) -> object:
    if (
        projection_worker_count() == 1
        or mp.current_process().name != "MainProcess"
    ):
        return retained_grade_serial(value)
    return retained_grade_parallel(value)


def retained_energy_basis(value: sp.Tuple) -> sp.Tuple:
    return sp.Tuple(value[0], *(retained_grade(row) for row in value[1:]))


def multigrade(value: object) -> sp.Tuple:
    expressions: list[sp.Expr] = []

    def collect(item: object) -> None:
        if isinstance(item, Mapping):
            for key, member in item.items():
                collect(key)
                collect(member)
        elif isinstance(item, sp.MatrixBase):
            for member in item:
                collect(member)
        elif isinstance(item, Relational):
            collect(item.lhs)
            collect(item.rhs)
        elif isinstance(item, sp.Tuple) or isinstance(item, (tuple, list)):
            for member in item:
                collect(member)
        elif isinstance(item, sp.Expr):
            expressions.append(item)

    collect(value)
    def combine(
        left: set[tuple[int, int, int]], right: set[tuple[int, int, int]]
    ) -> set[tuple[int, int, int]]:
        return {
            (a[0] + b[0], a[1] + b[1], a[2] + b[2])
            for a in left
            for b in right
            if a[1] + b[1] <= 1 and a[2] + b[2] <= 1
        }

    def expression_grades(
        expression: sp.Basic, count_wave_symbols: bool
    ) -> set[tuple[int, int, int]]:
        if expression.is_Number or isinstance(expression, Str):
            return {(0, 0, 0)}
        if isinstance(expression, sp.Symbol):
            if expression == epsilon:
                return {(1, 0, 0)}
            if count_wave_symbols and expression in WAVE_SYMBOLS:
                return {(1, 0, 0)}
            if expression in {W_bg, mu_R_bg, rhobr_rho4, rho4_rhobr}:
                return {(0, 0, 0), (0, 1, 0)}
            if expression in {*grad_W, *grad_mu}:
                return {(0, 0, 1)}
            return {(0, 0, 0)}
        if isinstance(expression, sp.Add):
            return set().union(
                *(expression_grades(argument, count_wave_symbols) for argument in expression.args)
            )
        if isinstance(expression, sp.Mul):
            result = {(0, 0, 0)}
            for argument in expression.args:
                result = combine(result, expression_grades(argument, count_wave_symbols))
            return result
        if isinstance(expression, sp.Pow):
            base_grades = expression_grades(expression.base, count_wave_symbols)
            if expression.exp.is_Integer and expression.exp >= 0:
                result = {(0, 0, 0)}
                for _ in range(int(expression.exp)):
                    result = combine(result, base_grades)
                return result
            # The first-shape-order series of an inverse or fractional power
            # has the base's zero-jet and first-jet supports.
            return base_grades
        if isinstance(expression, Relational):
            return expression_grades(expression.lhs, count_wave_symbols) | expression_grades(
                expression.rhs, count_wave_symbols
            )
        if expression.args:
            return set().union(
                *(expression_grades(argument, count_wave_symbols) for argument in expression.args)
            )
        return {(0, 0, 0)}

    grades: set[tuple[int, int, int]] = set()
    for expression in expressions:
        grades.update(expression_grades(expression, not expression.has(epsilon)))
    return sp.Tuple(*(sp.Tuple(*grade) for grade in sorted(grades)))


def dimension_of(expression: sp.Expr, dimensions: Mapping[sp.Symbol, sp.ImmutableMatrix] | None = None) -> sp.ImmutableMatrix:
    dimension_map = SYMBOL_DIMENSIONS if dimensions is None else dimensions
    expression = sp.sympify(expression)
    if expression.is_Number or isinstance(expression, Str):
        return DIM_ZERO
    if isinstance(expression, sp.Symbol):
        if expression not in dimension_map:
            raise KeyError(f"dimension unavailable for {sp.srepr(expression)}")
        return sp.ImmutableMatrix(dimension_map[expression])
    if isinstance(expression, sp.Add):
        items = tuple(dimension_of(argument, dimension_map) for argument in expression.args)
        return items[0] if all(item == items[0] for item in items[1:]) else sp.ImmutableMatrix([sp.nan] * 3)
    if isinstance(expression, sp.Mul):
        return sp.ImmutableMatrix(
            sum((sp.Matrix(dimension_of(argument, dimension_map)) for argument in expression.args), sp.zeros(3, 1))
        )
    if isinstance(expression, sp.Pow) and expression.exp.is_Number:
        return sp.ImmutableMatrix(expression.exp * sp.Matrix(dimension_of(expression.base, dimension_map)))
    if isinstance(expression, Relational):
        left = dimension_of(expression.lhs, dimension_map)
        right = dimension_of(expression.rhs, dimension_map)
        return left if left == right else sp.ImmutableMatrix([sp.nan] * 3)
    if isinstance(expression, Boolean):
        return DIM_ZERO
    raise TypeError(f"dimension trace has no rule for {sp.srepr(expression)}")


def dimension_object(value: object, dimensions: Mapping[sp.Symbol, sp.ImmutableMatrix] | None = None) -> object:
    if isinstance(value, Mapping):
        return {key: dimension_object(item, dimensions) for key, item in value.items()}
    if isinstance(value, sp.MatrixBase):
        return sp.Tuple(*(dimension_object(item, dimensions) for item in value))
    if isinstance(value, sp.Tuple) or isinstance(value, (tuple, list)):
        return sp.Tuple(*(dimension_object(item, dimensions) for item in value))
    if isinstance(value, (sp.Expr, Relational, Boolean)):
        try:
            return dimension_of(value, dimensions)
        except (KeyError, TypeError):
            return sp.ImmutableMatrix([sp.nan] * 3)
    return DIM_ZERO


def case_payload(raw_cases: Mapping[object, object], dimensions: Mapping[object, object] | None = None) -> sp.Tuple:
    rows = []
    for key, value in raw_cases.items():
        dimension = dimension_object(value) if dimensions is None else dimensions[key]
        rows.append(
            sp.Tuple(
                casify(key),
                sp.Tuple(
                    sp.Tuple(Str("VALUE"), casify(value)),
                    sp.Tuple(Str("MULTIGRADE"), multigrade(value)),
                    sp.Tuple(Str("DIMENSION_L_T_M"), casify(dimension)),
                ),
            )
        )
    return sp.Tuple(*rows)


def object_difference(left: object, right: object) -> object:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        keys = tuple(dict.fromkeys((*left.keys(), *right.keys())))
        return {
            key: object_difference(left.get(key, sp.Tuple()), right.get(key, sp.Tuple()))
            for key in keys
        }
    if isinstance(left, sp.MatrixBase) and isinstance(right, sp.MatrixBase) and left.shape == right.shape:
        return sp.ImmutableMatrix(left - right)
    if isinstance(left, Relational) and isinstance(right, Relational):
        return sp.expand((left.lhs - left.rhs) - (right.lhs - right.rhs))
    if isinstance(left, Boolean) or isinstance(right, Boolean):
        return sp.Equivalent(casify(left), casify(right), evaluate=False)
    left_container = isinstance(left, (tuple, list, sp.Tuple))
    right_container = isinstance(right, (tuple, list, sp.Tuple))
    if left_container and right_container and len(left) == len(right):
        return sp.Tuple(*(object_difference(a, b) for a, b in zip(left, right)))
    if isinstance(left, sp.Basic) and isinstance(right, sp.Basic):
        try:
            return sp.expand(left - right)
        except TypeError:
            return sp.Tuple(casify(left), casify(right))
    return sp.Tuple(casify(left), casify(right))


def structural_residual(left: object, right: object) -> object:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        keys = tuple(dict.fromkeys((*left.keys(), *right.keys())))
        return {
            key: structural_residual(
                left.get(key, sp.Tuple()), right.get(key, sp.Tuple())
            )
            for key in keys
        }
    left_container = isinstance(left, (tuple, list, sp.Tuple))
    right_container = isinstance(right, (tuple, list, sp.Tuple))
    if left_container and right_container and len(left) == len(right):
        return sp.Tuple(*(structural_residual(a, b) for a, b in zip(left, right)))
    if isinstance(left, sp.Basic) and isinstance(right, sp.Basic):
        if sp.srepr(left) == sp.srepr(right):
            return sp.Integer(0)
    return object_difference(left, right)


def perfect_matchings(items: tuple[tuple[int, int], ...]) -> tuple[tuple[tuple[tuple[int, int], tuple[int, int]], ...], ...]:
    if not items:
        return ((),)
    first = items[0]
    result = []
    for index in range(1, len(items)):
        pair = (first, items[index])
        rest = items[1:index] + items[index + 1 :]
        for matching in perfect_matchings(rest):
            result.append((pair,) + matching)
    return tuple(result)


def tensor_component(data: object, indices: tuple[int, ...]) -> sp.Expr:
    if not indices:
        return sp.sympify(data)
    if len(indices) == 1:
        return data[indices[0]]  # type: ignore[index]
    return data[indices[0]][indices[1]]  # type: ignore[index]


def delta_contractions(factors: tuple[tuple[str, int, object], ...]) -> tuple[sp.Expr, ...]:
    slots = tuple((factor, index) for factor, (_, rank, _) in enumerate(factors) for index in range(rank))
    if len(slots) % 2:
        return ()
    expressions = []
    for matching in perfect_matchings(slots):
        total = sp.Integer(0)
        for assignment in product(range(3), repeat=len(matching)):
            index_map = {
                slot: assignment[pair_index]
                for pair_index, pair in enumerate(matching)
                for slot in pair
            }
            term = sp.Integer(1)
            for factor_index, (_, rank, data) in enumerate(factors):
                indices = tuple(index_map[(factor_index, slot)] for slot in range(rank))
                term *= tensor_component(data, indices)
            total += term
        expressions.append(sp.expand(total))
    return tuple(dict.fromkeys(expressions))


def basis_euler_signatures(
    candidates: tuple[sp.Expr, ...],
    fields: tuple[tuple[sp.Symbol, tuple[sp.Symbol, ...], Mapping[tuple[int, ...], sp.Symbol]], ...],
    *,
    background_first_jets: tuple[sp.Symbol, ...] = (),
    background_second_jets: Mapping[tuple[int, ...], sp.Symbol] | None = None,
    background_depth: int = STRONG_ROW_JET_DEPTH,
    background_jet_scale: sp.Expr = sp.Integer(1),
) -> tuple[tuple[sp.Expr, ...], ...]:
    derivative_maps = {i: {} for i in DIRECTIONS}
    for field, first, second in fields:
        for i in DIRECTIONS:
            derivative_maps[i][field] = first[i]
            for j in DIRECTIONS:
                derivative_maps[i][first[j]] = second[sorted_index(i, j)]
    if background_first_jets and background_depth >= 2:
        if background_second_jets is None:
            raise ValueError("background second-jet table is required for a live spurion")
        for i in DIRECTIONS:
            for j in DIRECTIONS:
                derivative_maps[i][background_first_jets[j]] = (
                    background_jet_scale
                    * background_second_jets[sorted_index(i, j)]
                )

    def basis_dx(expression: sp.Expr, direction: int) -> sp.Expr:
        return sp.expand(
            sp.Add(
                *(
                    sp.diff(expression, atom) * derivative
                    for atom, derivative in derivative_maps[direction].items()
                    if atom in expression.free_symbols
                )
            )
        )

    signatures = []
    for candidate in candidates:
        signature = []
        for field, first, _ in fields:
            signature.append(
                sp.expand(
                    sp.diff(candidate, field)
                    - sp.Add(*(basis_dx(sp.diff(candidate, first[i]), i) for i in DIRECTIONS))
                )
            )
        signatures.append(tuple(signature))
    return tuple(signatures)


def quotient_independent_indices(
    candidates: tuple[sp.Expr, ...],
    signatures: tuple[tuple[sp.Expr, ...], ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if not candidates:
        return (), ()
    variables = sorted(
        {symbol for signature in signatures for expression in signature for symbol in expression.free_symbols},
        key=sp.default_sort_key,
    )
    monomials = sorted(
        {
            monomial
            for signature in signatures
            for expression in signature
            for monomial in sp.Poly(expression, *variables).monoms()
        }
    )
    matrix = sp.zeros(len(signatures[0]) * len(monomials), 0)
    rank = 0
    selected = []
    for index, signature in enumerate(signatures):
        column = []
        for expression in signature:
            polynomial = sp.Poly(expression, *variables)
            column.extend(polynomial.coeff_monomial(monomial) for monomial in monomials)
        trial = matrix.row_join(sp.Matrix(column))
        trial_rank = trial.rank()
        if trial_rank > rank:
            selected.append(index)
            matrix = trial
            rank = trial_rank
    omitted = tuple(index for index in range(len(candidates)) if index not in selected)
    return tuple(selected), omitted


# Abstract tensor data used only to enumerate and rank invariants.  Every
# selected expression is subsequently mapped to the live imported fields.
bu = sp.symbols("bu_1:4")
bG = tuple(tuple(sp.symbols(f"bG_{a + 1}_{i + 1}") for i in DIRECTIONS) for a in DIRECTIONS)
btheta = sp.Symbol("btheta")
bq = sp.symbols("bq_1:4")
be = sp.Symbol("be_local")
br = sp.symbols("br_local_1:4")
bg = sp.symbols("bg_1:4")


def basis_second(prefix: str, component: int | None = None) -> dict[tuple[int, ...], sp.Symbol]:
    return {
        (i, j): sp.Symbol(f"{prefix}_{'' if component is None else component + 1}_{i + 1}_{j + 1}")
        for i in DIRECTIONS
        for j in range(i, 3)
    }


basis_fields = tuple(
    [(bu[a], bG[a], basis_second("bU2", a)) for a in DIRECTIONS]
    + [(btheta, bq, basis_second("bTheta2")), (be, br, basis_second("bE2"))]
)
basis_background_second = basis_second("bBackground2")


def unique_expressions(expressions: list[sp.Expr]) -> tuple[sp.Expr, ...]:
    result = []
    for expression in expressions:
        expression = sp.expand(expression)
        if expression != 0 and expression not in result:
            result.append(expression)
    return tuple(result)


def enumerate_uniform_candidates() -> tuple[tuple[str, sp.Expr], ...]:
    data = (
        ("GRAD_U", 2, bG),
        ("THETA", 0, btheta),
        ("GRAD_THETA", 1, bq),
        ("E_LOCAL", 0, be),
        ("GRAD_E_LOCAL", 1, br),
    )
    raw = []
    for left_index, left in enumerate(data):
        for right in data[left_index:]:
            raw.extend(delta_contractions((left, right)))
    antisymmetric = sp.expand(
        sp.Add(*((bG[i][j] - bG[j][i]) ** 2 for i in DIRECTIONS for j in range(i + 1, 3)))
    )
    symmetric_tracefree = sp.expand(
        sp.Rational(1, 2)
        * sp.Add(*((bG[i][j] + bG[j][i]) ** 2 for i in DIRECTIONS for j in DIRECTIONS))
        - sp.Rational(2, 3) * trace(bG) ** 2
    )
    ordered = unique_expressions([antisymmetric, symmetric_tracefree, *raw])
    labels = []
    for index, expression in enumerate(ordered):
        if expression == antisymmetric:
            label = "CURL_SQUARED"
        elif expression == symmetric_tracefree:
            label = "SYMMETRIC_TRACEFREE_SQUARED"
        else:
            label = f"UNIFORM_CONTRACTION_{index + 1:02d}"
        labels.append((label, expression))
    return tuple(labels)


def enumerate_new_candidates(g_vector: tuple[sp.Expr, ...]) -> tuple[tuple[str, sp.Expr], ...]:
    data = (
        ("U", 1, bu),
        ("GRAD_U", 2, bG),
        ("THETA", 0, btheta),
        ("GRAD_THETA", 1, bq),
        ("E_LOCAL", 0, be),
        ("GRAD_E_LOCAL", 1, br),
    )
    spurion = ("BACKGROUND_FIRST_JET", 1, g_vector)
    raw = []
    for left_index, left in enumerate(data):
        for right in data[left_index:]:
            raw.extend(delta_contractions((spurion, left, right)))
    return tuple(
        (f"FIRST_JET_CONTRACTION_{index + 1:02d}", expression)
        for index, expression in enumerate(unique_expressions(raw))
    )


UNIFORM_CANDIDATES = enumerate_uniform_candidates()
UNIFORM_SIGNATURES = basis_euler_signatures(
    tuple(expression for _, expression in UNIFORM_CANDIDATES), basis_fields
)
UNIFORM_SELECTED, UNIFORM_OMITTED = quotient_independent_indices(
    tuple(expression for _, expression in UNIFORM_CANDIDATES), UNIFORM_SIGNATURES
)


@dataclass(frozen=True)
class EnergyBuild:
    density: sp.Expr
    basis: sp.Tuple
    new_invariants: sp.Tuple
    omissions: sp.Tuple
    count: sp.Integer
    new_counts: sp.Tuple
    terms: tuple[tuple[str, sp.Expr], ...]


# Probe jets are installed only after the common differential machinery is
# available.  Their local curl/divergence constraints are algebraic first-jet
# identities and introduce no projector.
uT_second = tuple(symmetric_jet(f"u_T_{a + 1}", 2, -DIM_L) for a in DIRECTIONS)
uT_third = tuple(symmetric_jet(f"u_T_{a + 1}", 3, -2 * DIM_L) for a in DIRECTIONS)
uL_second = tuple(symmetric_jet(f"u_L_{a + 1}", 2, -DIM_L) for a in DIRECTIONS)
uL_third = tuple(symmetric_jet(f"u_L_{a + 1}", 3, -2 * DIM_L) for a in DIRECTIONS)
theta_probe_second = symmetric_jet("theta_probe", 2, -2 * DIM_L)
theta_probe_third = symmetric_jet("theta_probe", 3, -3 * DIM_L)
e_probe_second = symmetric_jet("e_W_probe", 2, -2 * DIM_L)
e_probe_third = symmetric_jet("e_W_probe", 3, -3 * DIM_L)
for a in DIRECTIONS:
    add_field_derivatives(uT[a], GT[a], uT_second[a], uT_third[a])
    add_field_derivatives(uL[a], GL[a], uL_second[a], uL_third[a])
add_field_derivatives(theta_probe, theta_probe_grad, theta_probe_second, theta_probe_third)
add_field_derivatives(e_probe, e_probe_grad, e_probe_second, e_probe_third)


def uniform_coefficient(label: str, invariant: sp.Expr) -> sp.Expr:
    div_bu = trace(bG)
    catalogue = (
        (btheta * div_bu, G_theta_u),
        (be * div_bu, G_W_u),
        (btheta**2, B_rho_3 * W_bg / (2 * W0)),
        (btheta * be, C * W_bg),
        (be**2, k_W * W_bg**2 / 2),
        (dot(bq, bq), kappa_theta / 2),
        (dot(bq, br), kappa_theta_W),
        (dot(br, br), kappa_W * W_bg**4 / 2),
    )
    if label == "CURL_SQUARED":
        return mu_R_bg / 2
    if label == "SYMMETRIC_TRACEFREE_SQUARED":
        return mu_S / 2
    for candidate, coefficient in catalogue:
        if sp.expand(invariant - candidate) == 0:
            return coefficient
    # This branch records an independently admitted uniform invariant without
    # assigning it the identity of an upstream coefficient.
    index = next(
        i for i, (candidate_label, _) in enumerate(UNIFORM_CANDIDATES, 1) if candidate_label == label
    )
    abstract_dimensions = {
        btheta: DIM_ZERO,
        be: DIM_ZERO,
        **{bu[i]: DIM_L for i in DIRECTIONS},
        **{bq[i]: -DIM_L for i in DIRECTIONS},
        **{br[i]: -DIM_L for i in DIRECTIONS},
        **{bG[a][i]: DIM_ZERO for a in DIRECTIONS for i in DIRECTIONS},
    }
    invariant_dimension = dimension_of(invariant, abstract_dimensions)
    coefficient = symbol(
        f"a_s11cb_uniform_{index:02d}",
        "KNOB",
        f"coefficient for constructed uniform invariant {index}",
        dim_div(DIM_ENERGY, invariant_dimension),
    )
    return coefficient


def local_thickness_map(
    background_depth: int = STRONG_ROW_JET_DEPTH,
    background_first: Mapping[str, tuple[sp.Expr, ...]] | None = None,
    zero_sources: frozenset[str] = frozenset(),
) -> tuple[sp.Expr, tuple[sp.Expr, ...]]:
    local_value = W0 * e_W / W_bg
    local_gradient = tuple(
        total_derivative(
            local_value,
            i,
            background_depth=background_depth,
            background_first=background_first,
            zero_sources=zero_sources,
        )
        for i in DIRECTIONS
    )
    return local_value, local_gradient


def live_basis_substitution(
    g_vector: tuple[sp.Expr, ...],
    *,
    source: str | None = None,
    background_depth: int = STRONG_ROW_JET_DEPTH,
    background_first: Mapping[str, tuple[sp.Expr, ...]] | None = None,
    zero_sources: frozenset[str] = frozenset(),
) -> dict[sp.Symbol, sp.Expr]:
    local_value, local_gradient = local_thickness_map(
        background_depth, background_first, zero_sources
    )
    substitution: dict[sp.Symbol, sp.Expr] = {
        btheta: theta,
        be: local_value,
        **{bu[i]: u[i] for i in DIRECTIONS},
        **{bq[i]: grad_theta[i] for i in DIRECTIONS},
        **{br[i]: local_gradient[i] for i in DIRECTIONS},
        **{bg[i]: g_vector[i] for i in DIRECTIONS},
    }
    substitution.update(
        {bG[a][i]: grad_u[a][i] for a in DIRECTIONS for i in DIRECTIONS}
    )
    for a in DIRECTIONS:
        abstract_second = basis_fields[a][2]
        substitution.update(
            {
                abstract_second[index]: u_second[a][index]
                for index in abstract_second
            }
        )
    abstract_theta_second = basis_fields[3][2]
    substitution.update(
        {
            abstract_theta_second[index]: theta_second[index]
            for index in abstract_theta_second
        }
    )
    abstract_e_second = basis_fields[4][2]
    substitution.update(
        {
            abstract_e_second[index]: total_derivative(
                local_gradient[index[0]],
                index[1],
                background_depth=background_depth,
                background_first=background_first,
                zero_sources=zero_sources,
            )
            for index in abstract_e_second
        }
    )
    if source is not None:
        substitution.update(
            {
                basis_background_second[index]: (
                    sp.Integer(0)
                    if source in zero_sources or background_depth < 2
                    else background_jet_expression(source, 2, index)
                )
                for index in basis_background_second
            }
        )
    return substitution


def ablated_jets(
    source: str | None,
    direction: int | None,
    full_zero_source: str | None = None,
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]]:
    g_w: list[sp.Expr] = list(grad_W)
    g_mu: list[sp.Expr] = list(grad_mu)
    if full_zero_source == "W_BG":
        g_w = [sp.Integer(0)] * 3
    if full_zero_source == "MU_R_BG":
        g_mu = [sp.Integer(0)] * 3
    if source == "W_BG" and direction is not None:
        g_w[direction] = sp.Integer(0)
    if source == "MU_R_BG" and direction is not None:
        g_mu[direction] = sp.Integer(0)
    return tuple(g_w), tuple(g_mu)


def construct_energy(
    branch: str,
    ablation_source: str | None = None,
    ablation_direction: int | None = None,
    *,
    background_depth: int = STRONG_ROW_JET_DEPTH,
    full_zero_source: str | None = None,
) -> EnergyBuild:
    cache_key = (
        branch,
        ablation_source,
        ablation_direction,
        background_depth,
        full_zero_source,
    )
    if cache_key in ENERGY_CACHE:
        return ENERGY_CACHE[cache_key]

    g_w, g_mu = ablated_jets(
        ablation_source, ablation_direction, full_zero_source
    )
    background_first = {"W_BG": g_w, "MU_R_BG": g_mu}
    zero_sources = (
        frozenset((full_zero_source,))
        if full_zero_source is not None
        else frozenset()
    )
    uniform_substitution = live_basis_substitution(
        g_w,
        background_depth=background_depth,
        background_first=background_first,
        zero_sources=zero_sources,
    )
    selected_terms: list[tuple[str, sp.Expr]] = []
    basis_rows = []

    for index in UNIFORM_SELECTED:
        label, abstract_invariant = UNIFORM_CANDIDATES[index]
        invariant = sp.expand(abstract_invariant.subs(uniform_substitution, simultaneous=True))
        coefficient = uniform_coefficient(label, abstract_invariant)
        if ablation_source == "W_BG" and ablation_direction is not None:
            coefficient = coefficient.subs(grad_W[ablation_direction], 0)
        if ablation_source == "MU_R_BG" and ablation_direction is not None:
            coefficient = coefficient.subs(grad_mu[ablation_direction], 0)
        term = sp.expand(coefficient * invariant)
        selected_terms.append((label, term))
        coefficient_jet = sp.Tuple(
            *(
                total_derivative(
                    coefficient,
                    i,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
                for i in DIRECTIONS
            )
        )
        basis_rows.append(
            sp.Tuple(
                Str(label),
                invariant,
                coefficient,
                epsilon**2 * term,
                coefficient_jet,
                dimension_of(epsilon**2 * term),
            )
        )

    new_rows = []
    new_omissions = []
    for source, g_vector, actual_vector in (
        ("W_BG", tuple(bg[i] if g_w[i] != 0 else sp.Integer(0) for i in DIRECTIONS), g_w),
        ("MU_R_BG", tuple(bg[i] if g_mu[i] != 0 else sp.Integer(0) for i in DIRECTIONS), g_mu),
    ):
        candidates = enumerate_new_candidates(g_vector)
        expressions = tuple(expression for _, expression in candidates)
        signatures = basis_euler_signatures(
            expressions,
            basis_fields,
            background_first_jets=bg,
            background_second_jets=basis_background_second,
            background_depth=background_depth,
        )
        selected, omitted = quotient_independent_indices(expressions, signatures)
        substitution = live_basis_substitution(
            actual_vector,
            source=source,
            background_depth=background_depth,
            background_first=background_first,
            zero_sources=zero_sources,
        )
        for index in selected:
            label, abstract_invariant = candidates[index]
            invariant = sp.expand(abstract_invariant.subs(substitution, simultaneous=True))
            coefficient = NEW_COEFFICIENTS[source][index]
            invariant_dimension = dimension_of(invariant)
            SYMBOL_DIMENSIONS[coefficient] = dim_div(DIM_ENERGY, invariant_dimension)
            term = sp.expand(coefficient * invariant)
            full_label = f"{source}_{label}"
            selected_terms.append((full_label, term))
            new_rows.append(
                sp.Tuple(
                    Str(source),
                    Str(label),
                    invariant,
                    coefficient,
                    epsilon**2 * term,
                    dimension_of(epsilon**2 * term),
                )
            )
        for index in omitted:
            new_omissions.append(
                sp.Tuple(
                    Str(source),
                    Str(candidates[index][0]),
                    expressions[index].subs(substitution, simultaneous=True),
                    sp.Tuple(*signatures[index]).subs(substitution, simultaneous=True),
                )
            )

    uniform_omissions = tuple(
        sp.Tuple(
            Str(UNIFORM_CANDIDATES[index][0]),
            UNIFORM_CANDIDATES[index][1].subs(uniform_substitution, simultaneous=True),
            sp.Tuple(*UNIFORM_SIGNATURES[index]).subs(uniform_substitution, simultaneous=True),
        )
        for index in UNIFORM_OMITTED
    )
    density = sp.expand(sp.Add(*(term for _, term in selected_terms)))
    build = EnergyBuild(
        density=density,
        basis=sp.Tuple(profile_definitions(), *basis_rows, *new_rows),
        new_invariants=sp.Tuple(*new_rows),
        omissions=sp.Tuple(*uniform_omissions, *new_omissions),
        count=sp.Integer(len(basis_rows) + len(new_rows)),
        new_counts=sp.Tuple(
            *(
                sp.Tuple(
                    Str(source),
                    sp.Integer(sum(1 for row in new_rows if str(row[0]) == source)),
                )
                for source in ("W_BG", "MU_R_BG")
            )
        ),
        terms=tuple(selected_terms),
    )
    ENERGY_CACHE[cache_key] = build
    return build


def density_pair(
    representative: str,
    *,
    background_depth: int = STRONG_ROW_JET_DEPTH,
    background_first: Mapping[str, tuple[sp.Expr, ...]] | None = None,
    zero_sources: frozenset[str] = frozenset(),
) -> tuple[sp.Expr, sp.Expr, tuple[sp.Expr, ...]]:
    if representative == "RHO4_CONSTANT":
        rho4 = rho_br / W0
        rhobr = rho_br * W_bg / W0
    elif representative == "RHOBR_CONSTANT":
        rho4 = rho_br / W_bg
        rhobr = rho_br
    else:
        raise ValueError(representative)
    return rho4, rhobr, tuple(
        total_derivative(
            rhobr,
            i,
            background_depth=background_depth,
            background_first=background_first,
            zero_sources=zero_sources,
        )
        for i in DIRECTIONS
    )


def wave_truncate(expression: sp.Expr, degree: int) -> sp.Expr:
    scale = sp.Dummy("wave_scale")
    substitution = {item: scale * item for item in WAVE_SYMBOLS if item in expression.free_symbols}
    scaled = sp.expand(expression.subs(substitution, simultaneous=True))
    try:
        scaled = sp.series(scaled, scale, 0, degree + 1).removeO()
    except (NotImplementedError, ValueError, TypeError):
        pass
    return sp.expand(scaled.subs(scale, 1))


def material_pullback(
    density: sp.Expr,
    branch: str,
    representative: str,
    *,
    background_depth: int = STRONG_ROW_JET_DEPTH,
    background_first: Mapping[str, tuple[sp.Expr, ...]] | None = None,
    zero_sources: frozenset[str] = frozenset(),
) -> sp.Expr:
    rho4, _, _ = density_pair(
        representative,
        background_depth=background_depth,
        background_first=background_first,
        zero_sources=zero_sources,
    )
    rho4_gradient = tuple(
        total_derivative(
            rho4,
            i,
            background_depth=background_depth,
            background_first=background_first,
            zero_sources=zero_sources,
        )
        for i in DIRECTIONS
    )
    theta_material = theta + dot(u, rho4_gradient) / rho4
    if branch == "LAB_HELD":
        e_material = e_W + dot(u, grad_W) / W_bg
    else:
        e_material = e_W
    substitutions: dict[sp.Symbol, sp.Expr] = {
        theta: theta_material,
        e_W: e_material,
        **{
            grad_theta[i]: total_derivative(
                theta_material,
                i,
                background_depth=background_depth,
                background_first=background_first,
                zero_sources=zero_sources,
            )
            for i in DIRECTIONS
        },
        **{
            grad_e[i]: total_derivative(
                e_material,
                i,
                background_depth=background_depth,
                background_first=background_first,
                zero_sources=zero_sources,
            )
            for i in DIRECTIONS
        },
    }
    pulled = density.subs(substitutions, simultaneous=True)
    jacobian = 1 + trace(grad_u)
    scale = sp.Dummy("material_wave_scale")
    scaled = (jacobian * pulled).subs(
        {
            item: scale * item
            for item in WAVE_SYMBOLS
            if item in (jacobian * pulled).free_symbols
        },
        simultaneous=True,
    )
    return sp.expand(sp.diff(scaled, scale, 2).subs(scale, 0) / 2)


def case_axes(row: sp.Tuple) -> tuple[object, ...]:
    return tuple(int(item) if isinstance(item, sp.Integer) else str(item) for item in row[0])


def metadata_value(row: sp.Tuple) -> object:
    for key, value in row[1]:
        if str(key) == "VALUE":
            return value
    raise KeyError("VALUE")


def substrate_substitutions(
    branch: str,
    mu_theta_amplitude: sp.Expr,
    ablation_source: str | None,
    ablation_direction: int | None,
    full_zero_source: str | None,
) -> dict[sp.Symbol, sp.Expr]:
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    if ablation_source == "W_BG" and ablation_direction is not None:
        substitutions[w1_grad[ablation_direction]] = sp.Integer(0)
        substitutions[grad_W[ablation_direction]] = sp.Integer(0)
    if ablation_source == "MU_R_BG" and ablation_direction is not None:
        substitutions[m1_grad[ablation_direction]] = sp.Integer(0)
        substitutions[grad_mu[ablation_direction]] = sp.Integer(0)
    if full_zero_source is not None:
        substitutions.update(
            {
                atom: sp.Integer(0)
                for atom in (
                    *BACKGROUND_FIRST_JETS[full_zero_source],
                    *(
                        profile_atom
                        for table in BACKGROUND_PROFILE_JETS[full_zero_source].values()
                        for profile_atom in table.values()
                    ),
                )
            }
        )
    return substitutions


def filtered_substrate(
    name: str,
    branch: str,
    representative: str,
    mu_theta_amplitude: sp.Expr,
    ablation_source: str | None = None,
    ablation_direction: int | None = None,
    full_zero_source: str | None = None,
) -> sp.Tuple:
    substitutions = substrate_substitutions(
        branch,
        mu_theta_amplitude,
        ablation_source,
        ablation_direction,
        full_zero_source,
    )
    rows = []
    for row in S11CA_SUBSTRATE[name]:
        axes = case_axes(row)
        if any(axis in BRANCHES for axis in axes) and branch not in axes:
            continue
        if any(axis in DENSITY_REPS for axis in axes) and representative not in axes:
            continue
        value = metadata_value(row)
        value = map_object(
            value,
            lambda expression: expression.subs(substitutions, simultaneous=True),
        )
        rows.append(sp.Tuple(row[0], casify(value)))
    return sp.Tuple(*rows)


def substrate_bundle(
    branch: str,
    representative: str,
    mu_theta_amplitude: sp.Expr,
    ablation_source: str | None = None,
    ablation_direction: int | None = None,
    full_zero_source: str | None = None,
) -> sp.Tuple:
    return sp.Tuple(
        *(
            sp.Tuple(
                Str(name),
                filtered_substrate(
                    name,
                    branch,
                    representative,
                    mu_theta_amplitude,
                    ablation_source,
                    ablation_direction,
                    full_zero_source,
                ),
            )
            for name in SHAPE_SUBSTRATE_KEYS
        )
    )


def operator_from_density(
    density: sp.Expr,
    representative: str,
    include_kinetic: bool = True,
    *,
    background_depth: int = STRONG_ROW_JET_DEPTH,
    background_first: Mapping[str, tuple[sp.Expr, ...]] | None = None,
    zero_sources: frozenset[str] = frozenset(),
) -> tuple[dict[str, object], sp.Expr]:
    _, rhobr, rhobr_gradient = density_pair(
        representative,
        background_depth=background_depth,
        background_first=background_first,
        zero_sources=zero_sources,
    )
    u_local = tuple(epsilon * sp.diff(density, u[a]) for a in DIRECTIONS)
    u_flux = tuple(
        tuple(epsilon * sp.diff(density, grad_u[a][i]) for i in DIRECTIONS)
        for a in DIRECTIONS
    )
    u_expanded = tuple(
        sp.expand(
            u_local[a]
            - sp.Add(
                *(
                    total_derivative(
                        u_flux[a][i],
                        i,
                        background_depth=background_depth,
                        background_first=background_first,
                        zero_sources=zero_sources,
                    )
                    for i in DIRECTIONS
                )
            )
            - (epsilon * rhobr * u_tt[a] if include_kinetic else 0)
        )
        for a in DIRECTIONS
    )
    theta_local = epsilon * sp.diff(density, theta)
    theta_flux = tuple(epsilon * sp.diff(density, grad_theta[i]) for i in DIRECTIONS)
    mu_theta_amplitude = sp.expand(
        sp.diff(density, theta)
        - sp.Add(
            *(
                total_derivative(
                    sp.diff(density, grad_theta[i]),
                    i,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
                for i in DIRECTIONS
            )
        )
    )
    theta_expanded = epsilon * mu_theta_amplitude
    e_local = epsilon * sp.diff(density, e_W)
    e_flux = tuple(epsilon * sp.diff(density, grad_e[i]) for i in DIRECTIONS)
    e_expanded = sp.expand(
        e_local
        - sp.Add(
            *(
                total_derivative(
                    e_flux[i],
                    i,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
                for i in DIRECTIONS
            )
        )
        - (epsilon * mu_W * W_bg**2 * e_tt if include_kinetic else 0)
    )
    advective_constraint = epsilon * dot(u_t, rhobr_gradient)
    operator = {
        "U_BODY_BALANCE": sp.Tuple(
            sp.Tuple(Str("LOCAL"), sp.Tuple(*u_local)),
            sp.Tuple(Str("DIVERGENCE_FLUX"), sp.Tuple(*(sp.Tuple(*row) for row in u_flux))),
            sp.Tuple(Str("EXPANDED"), sp.Tuple(*u_expanded)),
        ),
        "THETA_BALANCE": sp.Tuple(
            sp.Tuple(Str("LOCAL"), theta_local),
            sp.Tuple(Str("DIVERGENCE_FLUX"), sp.Tuple(*theta_flux)),
            sp.Tuple(Str("EXPANDED"), theta_expanded),
        ),
        "E_W_BALANCE": sp.Tuple(
            sp.Tuple(Str("LOCAL"), e_local),
            sp.Tuple(Str("DIVERGENCE_FLUX"), sp.Tuple(*e_flux)),
            sp.Tuple(Str("EXPANDED"), e_expanded),
        ),
        "ADVECTIVE_MASS_OPERAND": advective_constraint,
    }
    return operator, mu_theta_amplitude


def named_tuple_row(rows: sp.Tuple, name: str) -> object:
    for row_name, row_value in rows:
        if str(row_name) == name:
            return row_value
    raise KeyError(name)


def strong_rows(operator: object) -> object:
    rows = casify(operator)
    return sp.Tuple(
        *(
            sp.Tuple(
                Str(name),
                named_tuple_row(named_tuple_row(rows, balance), "EXPANDED"),
            )
            for name, balance in (
                ("U", "U_BODY_BALANCE"),
                ("THETA", "THETA_BALANCE"),
                ("E_W", "E_W_BALANCE"),
            )
        )
    )


def committed_material_pullback(
    density: sp.Expr,
    branch: str,
    representative: str,
) -> sp.Expr:
    rho4, _, _ = density_pair(
        representative, background_depth=HESSIAN_FREEZE_DEPTH
    )
    rho4_gradient = tuple(frozen_derivative(rho4, i) for i in DIRECTIONS)
    theta_material = theta + dot(u, rho4_gradient) / rho4
    e_material = e_W + dot(u, grad_W) / W_bg if branch == "LAB_HELD" else e_W
    substitutions: dict[sp.Symbol, sp.Expr] = {
        theta: theta_material,
        e_W: e_material,
        **{
            grad_theta[i]: frozen_derivative(theta_material, i)
            for i in DIRECTIONS
        },
        **{grad_e[i]: frozen_derivative(e_material, i) for i in DIRECTIONS},
    }
    pulled = density.subs(substitutions, simultaneous=True)
    jacobian = 1 + trace(grad_u)
    scale = sp.Dummy("committed_material_wave_scale")
    scaled = (jacobian * pulled).subs(
        {
            item: scale * item
            for item in WAVE_SYMBOLS
            if item in (jacobian * pulled).free_symbols
        },
        simultaneous=True,
    )
    return sp.expand(sp.diff(scaled, scale, 2).subs(scale, 0) / 2)


def committed_strong_rows(
    density: sp.Expr,
    branch: str,
    representative: str,
    route: str,
) -> object:
    if route == "MATERIAL":
        density = committed_material_pullback(density, branch, representative)
    _, rhobr, _ = density_pair(
        representative, background_depth=HESSIAN_FREEZE_DEPTH
    )
    u_rows = tuple(
        sp.expand(
            epsilon * sp.diff(density, u[a])
            - sp.Add(
                *(
                    frozen_derivative(
                        epsilon * sp.diff(density, grad_u[a][i]), i
                    )
                    for i in DIRECTIONS
                )
            )
            - epsilon * rhobr * u_tt[a]
        )
        for a in DIRECTIONS
    )
    theta_row = sp.expand(
        epsilon
        * (
            sp.diff(density, theta)
            - sp.Add(
                *(
                    frozen_derivative(sp.diff(density, grad_theta[i]), i)
                    for i in DIRECTIONS
                )
            )
        )
    )
    e_row = sp.expand(
        epsilon * sp.diff(density, e_W)
        - sp.Add(
            *(
                frozen_derivative(epsilon * sp.diff(density, grad_e[i]), i)
                for i in DIRECTIONS
            )
        )
        - epsilon * mu_W * W_bg**2 * e_tt
    )
    return sp.Tuple(
        sp.Tuple(Str("U"), sp.Tuple(*u_rows)),
        sp.Tuple(Str("THETA"), theta_row),
        sp.Tuple(Str("E_W"), e_row),
    )


def live_strong_rows(
    density: sp.Expr,
    branch: str,
    representative: str,
    route: str,
    background_depth: int,
) -> object:
    if route == "MATERIAL":
        density = material_pullback(
            density,
            branch,
            representative,
            background_depth=background_depth,
        )
    operator, _ = operator_from_density(
        density,
        representative,
        background_depth=background_depth,
    )
    return strong_rows(casify(operator))


def coupling_outer_rows(
    rows: object,
    background_depth: int,
) -> object:
    u_row = tuple(named_tuple_row(casify(rows), "U"))
    return sp.Tuple(
        sp.Tuple(Str("DIVERGENCE"), divergence(u_row, background_depth)),
        sp.Tuple(Str("CURL"), sp.Tuple(*curl(u_row, background_depth))),
    )


def build_operator(
    branch: str,
    representative: str,
    route: str = "EULERIAN",
    ablation_source: str | None = None,
    ablation_direction: int | None = None,
    corrupt_material_constraint: bool = False,
    *,
    background_depth: int = STRONG_ROW_JET_DEPTH,
    full_zero_source: str | None = None,
) -> tuple[object, object, object]:
    cache_key = (
        branch,
        representative,
        route,
        ablation_source,
        ablation_direction,
        corrupt_material_constraint,
        background_depth,
        full_zero_source,
    )
    if cache_key in OPERATOR_CACHE:
        return OPERATOR_CACHE[cache_key]
    energy = construct_energy(
        branch,
        ablation_source,
        ablation_direction,
        background_depth=background_depth,
        full_zero_source=full_zero_source,
    )
    g_w, g_mu = ablated_jets(
        ablation_source, ablation_direction, full_zero_source
    )
    background_first = {"W_BG": g_w, "MU_R_BG": g_mu}
    zero_sources = (
        frozenset((full_zero_source,))
        if full_zero_source is not None
        else frozenset()
    )
    density = energy.density
    if route == "MATERIAL":
        density = material_pullback(
            density,
            branch,
            representative,
            background_depth=background_depth,
            background_first=background_first,
            zero_sources=zero_sources,
        )
    operator, mu_theta_amplitude = operator_from_density(
        density,
        representative,
        background_depth=background_depth,
        background_first=background_first,
        zero_sources=zero_sources,
    )
    if corrupt_material_constraint:
        operator["ADVECTIVE_MASS_OPERAND"] = sp.Integer(0)
    if branch == "MATERIAL_ADVECTED":
        mu_theta_value = sp.Tuple(Str("mu_theta^alpha"), epsilon * mu_theta_amplitude)
    else:
        mu_theta_value = sp.Tuple(Str("mu_theta"), epsilon * mu_theta_amplitude)
    faces = substrate_bundle(
        branch,
        representative,
        mu_theta_amplitude,
        ablation_source,
        ablation_direction,
        full_zero_source,
    )
    operator["FACE_FLUX_BOUNDARY_OPERANDS"] = faces
    reserved_mu_theta = INCOMING_LEDGER[
        "mu_theta_M" if branch == "MATERIAL_ADVECTED" else "mu_theta_L"
    ]["value"]
    operator["MU_THETA_FACE_BINDING"] = sp.Tuple(
        reserved_mu_theta, mu_theta_amplitude
    )

    energy_origins = []
    for label, term in energy.terms:
        term_operator, _ = operator_from_density(
            term,
            representative,
            include_kinetic=False,
            background_depth=background_depth,
            background_first=background_first,
            zero_sources=zero_sources,
        )
        energy_origins.append(
            sp.Tuple(
                Str(label),
                casify(
                    {
                        "U": term_operator["U_BODY_BALANCE"],
                        "THETA": term_operator["THETA_BALANCE"],
                        "E_W": term_operator["E_W_BALANCE"],
                    }
                ),
            )
        )
    origins = sp.Tuple(
        sp.Tuple(Str("BULK_ENERGY"), sp.Tuple(*energy_origins)),
        sp.Tuple(
            Str("KINETIC"),
            sp.Tuple(
                sp.Tuple(
                    *(
                        epsilon
                        * density_pair(
                            representative,
                            background_depth=background_depth,
                            background_first=background_first,
                            zero_sources=zero_sources,
                        )[1]
                        * u_tt[a]
                        for a in DIRECTIONS
                    )
                ),
                epsilon * mu_W * W_bg**2 * e_tt,
            ),
        ),
        sp.Tuple(Str("FACE_FLUX"), faces),
        sp.Tuple(Str("ADVECTIVE"), operator["ADVECTIVE_MASS_OPERAND"]),
    )
    result = casify(operator), origins, mu_theta_value
    OPERATOR_CACHE[cache_key] = result
    return result


# Curl-free local probes are represented as derivatives of a local scalar
# potential.  This is finite-order jet algebra; it does not solve a Poisson
# equation or construct a global projection.
phi_L_hessian = {
    (i, j): symbol(
        f"phi_L_d{i + 1}d{j + 1}",
        "COORDINATE",
        f"local longitudinal potential Hessian {i + 1},{j + 1}",
        DIM_ZERO,
    )
    for i in DIRECTIONS
    for j in range(i, 3)
}
phi_L_third = {
    (i, j, k): symbol(
        f"phi_L_d{i + 1}d{j + 1}d{k + 1}",
        "COORDINATE",
        f"local longitudinal potential third jet {i + 1},{j + 1},{k + 1}",
        -DIM_L,
    )
    for i in DIRECTIONS
    for j in range(i, 3)
    for k in range(j, 3)
}
phi_L_fourth = {
    (i, j, k, ell): symbol(
        f"phi_L_d{i + 1}d{j + 1}d{k + 1}d{ell + 1}",
        "COORDINATE",
        f"local longitudinal potential fourth jet {i + 1},{j + 1},{k + 1},{ell + 1}",
        -2 * DIM_L,
    )
    for i in DIRECTIONS
    for j in range(i, 3)
    for k in range(j, 3)
    for ell in range(k, 3)
}


def all_actual_wave_symbols() -> set[sp.Symbol]:
    result = set(WAVE_SYMBOLS)
    result.update(item for row in u_third for item in row.values())
    result.update(theta_third.values())
    result.update(e_third.values())
    return result


ACTUAL_WAVE_SYMBOLS = all_actual_wave_symbols()


def transverse_constraint_substitutions() -> dict[sp.Symbol, sp.Expr]:
    substitutions: dict[sp.Symbol, sp.Expr] = {
        GT[2][2]: -GT[0][0] - GT[1][1]
    }
    for i in DIRECTIONS:
        substitutions[uT_second[2][sorted_index(2, i)]] = -uT_second[0][
            sorted_index(0, i)
        ] - uT_second[1][sorted_index(1, i)]
    for i in DIRECTIONS:
        for j in range(i, 3):
            substitutions[uT_third[2][sorted_index(2, i, j)]] = -uT_third[0][
                sorted_index(0, i, j)
            ] - uT_third[1][sorted_index(1, i, j)]
    return substitutions


TRANSVERSE_CONSTRAINT_SUBS = transverse_constraint_substitutions()


def sector_substitution(sector: str) -> dict[sp.Symbol, sp.Expr]:
    substitutions: dict[sp.Symbol, sp.Expr] = {item: sp.Integer(0) for item in ACTUAL_WAVE_SYMBOLS}
    if sector == "TRANSVERSE":
        substitutions.update({u[a]: uT[a] for a in DIRECTIONS})
        substitutions.update({u_t[a]: uT_t[a] for a in DIRECTIONS})
        substitutions.update(
            {grad_u[a][i]: GT[a][i] for a in DIRECTIONS for i in DIRECTIONS}
        )
        substitutions.update(
            {
                u_second[a][index]: uT_second[a][index]
                for a in DIRECTIONS
                for index in u_second[a]
            }
        )
        substitutions.update(
            {
                u_third[a][index]: uT_third[a][index]
                for a in DIRECTIONS
                for index in u_third[a]
            }
        )
    elif sector == "LONGITUDINAL":
        substitutions.update({u[a]: uL[a] for a in DIRECTIONS})
        substitutions.update({u_t[a]: uL_t[a] for a in DIRECTIONS})
        substitutions.update(
            {
                grad_u[a][i]: phi_L_hessian[sorted_index(a, i)]
                for a in DIRECTIONS
                for i in DIRECTIONS
            }
        )
        substitutions.update(
            {
                u_second[a][index]: phi_L_third[sorted_index(a, *index)]
                for a in DIRECTIONS
                for index in u_second[a]
            }
        )
        substitutions.update(
            {
                u_third[a][index]: phi_L_fourth[sorted_index(a, *index)]
                for a in DIRECTIONS
                for index in u_third[a]
            }
        )
    elif sector == "THETA":
        substitutions[theta] = theta_probe
        substitutions.update({grad_theta[i]: theta_probe_grad[i] for i in DIRECTIONS})
        substitutions.update(
            {theta_second[index]: theta_probe_second[index] for index in theta_second}
        )
        substitutions.update(
            {theta_third[index]: theta_probe_third[index] for index in theta_third}
        )
    elif sector == "E_W":
        substitutions[e_W] = e_probe
        substitutions[e_t] = e_probe_t
        substitutions.update({grad_e[i]: e_probe_grad[i] for i in DIRECTIONS})
        substitutions.update({e_second[index]: e_probe_second[index] for index in e_second})
        substitutions.update({e_third[index]: e_probe_third[index] for index in e_third})
    else:
        raise ValueError(sector)
    return substitutions


def restrict_expression(expression: sp.Expr, sector: str) -> sp.Expr:
    restricted = sp.expand(expression.subs(sector_substitution(sector), simultaneous=True))
    if sector == "TRANSVERSE":
        restricted = sp.expand(restricted.subs(TRANSVERSE_CONSTRAINT_SUBS, simultaneous=True))
    return restricted


def restrict_object(value: object, sector: str) -> object:
    return map_object(value, lambda expression: restrict_expression(expression, sector))


def bulk_kernel_from_density(density: sp.Expr) -> dict[str, object]:
    u_el = tuple(euler_derivative(density, u[a], grad_u[a]) for a in DIRECTIONS)
    theta_el = euler_derivative(density, theta, grad_theta)
    e_el = euler_derivative(density, e_W, grad_e)
    u_curl = curl(u_el)
    u_divergence = divergence(u_el)
    return {
        "TRANSVERSE_TO_THICKNESS": {
            "THETA": epsilon * restrict_expression(theta_el, "TRANSVERSE"),
            "E_W": epsilon * restrict_expression(e_el, "TRANSVERSE"),
            "DIV_U": epsilon * restrict_expression(u_divergence, "TRANSVERSE"),
        },
        "THICKNESS_TO_TRANSVERSE": {
            "THETA": sp.Tuple(*(epsilon * restrict_expression(item, "THETA") for item in u_curl)),
            "E_W": sp.Tuple(*(epsilon * restrict_expression(item, "E_W") for item in u_curl)),
            "DIV_U": sp.Tuple(
                *(epsilon * restrict_expression(item, "LONGITUDINAL") for item in u_curl)
            ),
        },
    }


def combined_sector_substitution(
    target: str,
    transverse_scale: sp.Symbol,
    target_scale: sp.Symbol,
) -> dict[sp.Symbol, sp.Expr]:
    transverse = sector_substitution("TRANSVERSE")
    target_map = sector_substitution(target)
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    for atom in ACTUAL_WAVE_SYMBOLS:
        transverse_value = transverse.get(atom, 0)
        target_value = target_map.get(atom, 0)
        substitutions[atom] = transverse_scale * transverse_value + target_scale * target_value
    return substitutions


def mixed_variation(density: sp.Expr, target: str, reverse: bool = False) -> sp.Expr:
    first = sp.Dummy("pair_first")
    second = sp.Dummy("pair_second")
    if reverse:
        substitution = combined_sector_substitution(target, second, first)
        varied = sp.diff(sp.diff(density.subs(substitution, simultaneous=True), second), first)
    else:
        substitution = combined_sector_substitution(target, first, second)
        varied = sp.diff(sp.diff(density.subs(substitution, simultaneous=True), first), second)
    varied = sp.expand(varied.subs({first: 0, second: 0}))
    return sp.expand(varied.subs(TRANSVERSE_CONSTRAINT_SUBS, simultaneous=True))


def paired_kernel_from_density(density: sp.Expr) -> tuple[dict[str, object], sp.Tuple]:
    forward = {}
    reciprocal = {}
    adjoint_rows = []
    for target in ("THETA", "E_W", "LONGITUDINAL"):
        forward_value = epsilon * mixed_variation(density, target)
        reciprocal_value = epsilon * mixed_variation(density, target, reverse=True)
        target_label = "DIV_U" if target == "LONGITUDINAL" else target
        forward[target_label] = forward_value
        reciprocal[target_label] = reciprocal_value
        adjoint_rows.append(
            sp.Tuple(
                Str(target_label),
                sp.Tuple(Str("FORWARD_PAIRING_OPERAND"), forward_value),
                sp.Tuple(Str("RECIPROCAL_PAIRING_OPERAND"), reciprocal_value),
                sp.Tuple(
                    Str("PAIRING_RESIDUAL"),
                    sp.expand(forward_value - reciprocal_value),
                ),
            )
        )
    return {
        "TRANSVERSE_TO_THICKNESS": forward,
        "THICKNESS_TO_TRANSVERSE": reciprocal,
    }, sp.Tuple(*adjoint_rows)


def build_kernel(
    branch: str,
    representative: str,
    route: str = "EULERIAN",
    ablation_source: str | None = None,
    ablation_direction: int | None = None,
    include_term_origins: bool = True,
    corrupt_material_constraint: bool = False,
    *,
    background_depth: int = COUPLING_JET_DEPTH,
    full_zero_source: str | None = None,
) -> tuple[object, object]:
    cache_key = (
        branch,
        representative,
        route,
        ablation_source,
        ablation_direction,
        include_term_origins,
        corrupt_material_constraint,
        background_depth,
        full_zero_source,
    )
    if cache_key in KERNEL_CACHE:
        return KERNEL_CACHE[cache_key]
    operator, operator_origins, mu_theta_value = build_operator(
        branch,
        representative,
        route,
        ablation_source,
        ablation_direction,
        corrupt_material_constraint,
        background_depth=min(background_depth, STRONG_ROW_JET_DEPTH),
        full_zero_source=full_zero_source,
    )
    g_w, g_mu = ablated_jets(
        ablation_source, ablation_direction, full_zero_source
    )
    background_first = {"W_BG": g_w, "MU_R_BG": g_mu}
    zero_sources = (
        frozenset((full_zero_source,))
        if full_zero_source is not None
        else frozenset()
    )

    def named_row(rows: sp.Tuple, name: str) -> object:
        for row_name, row_value in rows:
            if str(row_name) == name:
                return row_value
        raise KeyError(name)

    def expanded_balance(balance: sp.Tuple) -> object:
        return named_row(balance, "EXPANDED")

    def operator_divergence(vector: tuple[sp.Expr, ...]) -> sp.Expr:
        return sp.expand(
            sp.Add(
                *(
                    total_derivative(
                        vector[i],
                        i,
                        background_depth=background_depth,
                        background_first=background_first,
                        zero_sources=zero_sources,
                    )
                    for i in DIRECTIONS
                )
            )
        )

    def operator_curl(vector: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
        return (
            sp.expand(
                total_derivative(
                    vector[2],
                    1,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
                - total_derivative(
                    vector[1],
                    2,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
            ),
            sp.expand(
                total_derivative(
                    vector[0],
                    2,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
                - total_derivative(
                    vector[2],
                    0,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
            ),
            sp.expand(
                total_derivative(
                    vector[1],
                    0,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
                - total_derivative(
                    vector[0],
                    1,
                    background_depth=background_depth,
                    background_first=background_first,
                    zero_sources=zero_sources,
                )
            ),
        )

    theta_test = symbol(
        "v_theta_s11cb",
        "COORDINATE",
        "independent compact-support densification test field",
        DIM_ZERO,
    )
    e_test = symbol(
        "v_e_W_s11cb",
        "COORDINATE",
        "independent compact-support thickness test field",
        DIM_ZERO,
    )
    longitudinal_test_potential = symbol(
        "psi_L_s11cb",
        "COORDINATE",
        "compact-support longitudinal test potential",
        2 * DIM_L,
    )
    transverse_test_potential = tuple(
        symbol(
            f"A_T_s11cb_{a + 1}",
            "COORDINATE",
            f"compact-support solenoidal test potential {a + 1}",
            2 * DIM_L,
        )
        for a in DIRECTIONS
    )

    def weak_operator_blocks(
        u_balance: sp.Tuple,
        theta_balance: sp.Tuple,
        e_balance: sp.Tuple,
    ) -> dict[str, object]:
        u_residual = tuple(expanded_balance(u_balance))
        theta_residual = sp.sympify(expanded_balance(theta_balance))
        e_residual = sp.sympify(expanded_balance(e_balance))
        u_transverse = tuple(
            restrict_expression(item, "TRANSVERSE") for item in u_residual
        )
        transverse_divergence = restrict_expression(
            operator_divergence(u_transverse), "TRANSVERSE"
        )
        target_curls = {
            target: operator_curl(
                tuple(restrict_expression(item, target) for item in u_residual)
            )
            for target in ("THETA", "E_W", "LONGITUDINAL")
        }
        return {
            "TRANSVERSE_TO_THICKNESS": {
                "THETA": sp.expand(
                    theta_test * restrict_expression(theta_residual, "TRANSVERSE")
                ),
                "E_W": sp.expand(
                    e_test * restrict_expression(e_residual, "TRANSVERSE")
                ),
                "DIV_U": sp.expand(
                    -longitudinal_test_potential
                    * transverse_divergence
                ),
            },
            "THICKNESS_TO_TRANSVERSE": {
                "THETA": dot(
                    transverse_test_potential,
                    target_curls["THETA"],
                ),
                "E_W": dot(
                    transverse_test_potential,
                    target_curls["E_W"],
                ),
                "DIV_U": dot(
                    transverse_test_potential,
                    target_curls["LONGITUDINAL"],
                ),
            },
        }

    operator_u = named_row(operator, "U_BODY_BALANCE")
    operator_theta = named_row(operator, "THETA_BALANCE")
    operator_e = named_row(operator, "E_W_BALANCE")
    core_key = (
        branch,
        route,
        representative if route == "MATERIAL" else None,
        ablation_source,
        ablation_direction,
        corrupt_material_constraint,
        background_depth,
        full_zero_source,
    )
    if core_key in KERNEL_BLOCK_CACHE:
        bulk, adjointness = KERNEL_BLOCK_CACHE[core_key]
    else:
        bulk = weak_operator_blocks(operator_u, operator_theta, operator_e)
        adjointness = sp.Tuple(
            Str("NO_INDEPENDENT_SECOND_ROUTE"),
            sp.Tuple(
                Str("TRANSVERSE_TO_THICKNESS_OPERATOR_BLOCK"),
                casify(bulk["TRANSVERSE_TO_THICKNESS"]),
            ),
            sp.Tuple(
                Str("THICKNESS_TO_TRANSVERSE_OPERATOR_BLOCK"),
                casify(bulk["THICKNESS_TO_TRANSVERSE"]),
            ),
        )
        KERNEL_BLOCK_CACHE[core_key] = bulk, adjointness
    if include_term_origins and core_key in KERNEL_ORIGIN_CACHE:
        origin_rows_value = KERNEL_ORIGIN_CACHE[core_key]
    elif include_term_origins:
        origin_rows = []
        bulk_energy_origins = named_row(operator_origins, "BULK_ENERGY")
        for label, term_operator in bulk_energy_origins:
            term_blocks = weak_operator_blocks(
                named_row(term_operator, "U"),
                named_row(term_operator, "THETA"),
                named_row(term_operator, "E_W"),
            )
            origin_rows.append(sp.Tuple(label, casify(term_blocks)))
        origin_rows_value = sp.Tuple(*origin_rows)
        KERNEL_ORIGIN_CACHE[core_key] = origin_rows_value
    else:
        origin_rows_value = sp.Tuple()
    mu_theta_amplitude = mu_theta_value[1] / epsilon
    face_bundle = named_row(operator, "FACE_FLUX_BOUNDARY_OPERANDS")
    boundary_source = sp.Tuple(
        *(
            sp.Tuple(
                Str(name),
                named_row(face_bundle, name),
            )
            for name in ("traction", "virtual_work_shape_deriv", "evolution_mass_balance")
        ),
        sp.Tuple(
            Str("MU_THETA_FACE_BINDING"),
            named_row(operator, "MU_THETA_FACE_BINDING"),
        ),
    )
    mass_source = named_row(face_bundle, "evolution_mass_balance")
    advective_source = named_row(operator, "ADVECTIVE_MASS_OPERAND")
    boundary = sp.Tuple(
        sp.Tuple(Str("SOURCE_SUBSTRATE"), boundary_source),
        sp.Tuple(
            Str("TRANSVERSE_TO_FACE_FLUX"),
            casify(restrict_object(mass_source, "TRANSVERSE")),
        ),
        sp.Tuple(
            Str("MU_THETA_TO_TRANSVERSE"),
            casify(
                {
                    "THETA": restrict_expression(mu_theta_amplitude, "THETA"),
                    "E_W": restrict_expression(mu_theta_amplitude, "E_W"),
                    "DIV_U": restrict_expression(mu_theta_amplitude, "LONGITUDINAL"),
                    "TRANSVERSE": restrict_expression(mu_theta_amplitude, "TRANSVERSE"),
                }
            ),
        ),
        sp.Tuple(
            Str("ADVECTIVE_MASS_BLOCK"),
            casify(
                {
                    "TRANSVERSE_TRIAL": restrict_expression(
                        advective_source, "TRANSVERSE"
                    ),
                    "LONGITUDINAL_TRIAL": restrict_expression(
                        advective_source, "LONGITUDINAL"
                    ),
                }
            ),
        ),
    )
    kernel = sp.Tuple(
        sp.Tuple(Str("SECTOR_LABELS"), sp.Tuple(Str("CURL_LOCAL"), Str("DIVERGENCE_LOCAL"))),
        sp.Tuple(Str("BULK_BLOCKS"), casify(bulk)),
        sp.Tuple(Str("FACE_FLUX_BLOCKS"), boundary),
        sp.Tuple(Str("VARIATIONAL_ADJOINTNESS"), adjointness),
    )
    origins = sp.Tuple(
        sp.Tuple(Str("BULK_ENERGY"), origin_rows_value),
        sp.Tuple(Str("FACE_FLUX"), boundary),
        sp.Tuple(
            Str("ADVECTIVE"),
            advective_source,
        ),
    )
    result = kernel, origins
    KERNEL_CACHE[cache_key] = result
    return result


def zero_wave(value: object) -> object:
    substitutions = {item: sp.Integer(0) for item in ACTUAL_WAVE_SYMBOLS}
    substitutions[epsilon] = sp.Integer(0)
    return map_object(
        value,
        lambda expression: sp.expand(expression.subs(substitutions, simultaneous=True)),
    )


def exact_substrate_value(name: str, axes: tuple[object, ...]) -> object:
    for row in S11CA_SUBSTRATE[name]:
        if case_axes(row) == axes:
            return metadata_value(row)
    raise KeyError((name, axes))


def admissibility_operator(
    branch: str,
    representative: str,
    route: str = "EULERIAN",
    ablation_source: str | None = None,
    ablation_direction: int | None = None,
    corrupt_material_constraint: bool = False,
    *,
    background_depth: int = STRONG_ROW_JET_DEPTH,
    full_zero_source: str | None = None,
) -> object:
    energy = construct_energy(
        branch,
        ablation_source,
        ablation_direction,
        background_depth=background_depth,
        full_zero_source=full_zero_source,
    )
    selected_labels = {label for label, _ in energy.terms}
    g_w, g_mu = ablated_jets(
        ablation_source, ablation_direction, full_zero_source
    )
    background_first = {"W_BG": g_w, "MU_R_BG": g_mu}
    zero_sources = (
        frozenset((full_zero_source,))
        if full_zero_source is not None
        else frozenset()
    )
    rho4, rhobr, rhobr_gradient = density_pair(
        representative,
        background_depth=background_depth,
        background_first=background_first,
        zero_sources=zero_sources,
    )

    def background_dx(expression: sp.Expr, direction: int) -> sp.Expr:
        return total_derivative(
            expression,
            direction,
            background_depth=background_depth,
            background_first=background_first,
            zero_sources=zero_sources,
        )

    rho4_gradient = tuple(background_dx(rho4, i) for i in DIRECTIONS)
    constraint_gradient = (
        (sp.Integer(0),) * 3 if corrupt_material_constraint else rhobr_gradient
    )
    theta_configuration = theta
    e_configuration = e_W
    if route == "EULERIAN":
        theta_configuration += dot(u, constraint_gradient) / rhobr
    elif route == "MATERIAL":
        theta_configuration += dot(u, rho4_gradient) / rho4
        if branch == "LAB_HELD":
            e_configuration += dot(u, g_w) / W_bg
    else:
        raise ValueError(route)

    theta_perturbation = theta_configuration
    theta_perturbation_gradient = tuple(
        background_dx(theta_configuration, i) for i in DIRECTIONS
    )
    thickness_perturbation = W0 * e_configuration / W_bg
    thickness_full_gradient = tuple(
        sp.expand((g_w[i] + W0 * background_dx(e_configuration, i)) / W_bg)
        for i in DIRECTIONS
    )

    def full_field_substitution(g_vector: tuple[sp.Expr, ...]) -> dict[sp.Symbol, sp.Expr]:
        substitutions: dict[sp.Symbol, sp.Expr] = {
            btheta: theta_perturbation,
            be: thickness_perturbation,
            **{bu[i]: u[i] for i in DIRECTIONS},
            **{bq[i]: theta_perturbation_gradient[i] for i in DIRECTIONS},
            **{br[i]: thickness_full_gradient[i] for i in DIRECTIONS},
            **{bg[i]: g_vector[i] for i in DIRECTIONS},
        }
        substitutions.update(
            {bG[a][i]: grad_u[a][i] for a in DIRECTIONS for i in DIRECTIONS}
        )
        return substitutions

    full_terms = []
    uniform_substitution = full_field_substitution(g_w)
    for index in UNIFORM_SELECTED:
        label, abstract_invariant = UNIFORM_CANDIDATES[index]
        if label not in selected_labels:
            continue
        full_terms.append(
            uniform_coefficient(label, abstract_invariant)
            * abstract_invariant.subs(uniform_substitution, simultaneous=True)
        )
    for source, abstract_g, actual_g in (
        ("W_BG", tuple(bg[i] if g_w[i] != 0 else sp.Integer(0) for i in DIRECTIONS), g_w),
        ("MU_R_BG", tuple(bg[i] if g_mu[i] != 0 else sp.Integer(0) for i in DIRECTIONS), g_mu),
    ):
        substitution = full_field_substitution(actual_g)
        for index, (label, abstract_invariant) in enumerate(
            enumerate_new_candidates(abstract_g)
        ):
            if f"{source}_{label}" not in selected_labels:
                continue
            full_terms.append(
                NEW_COEFFICIENTS[source][index]
                * abstract_invariant.subs(substitution, simultaneous=True)
            )
    full_density = sp.expand(sp.Add(*full_terms))
    if route == "MATERIAL":
        full_density = sp.expand((1 + trace(grad_u)) * full_density)

    def background_euler(
        field: sp.Symbol,
        first_jets: tuple[sp.Symbol, ...],
    ) -> sp.Expr:
        return sp.expand(
            sp.diff(full_density, field)
            - sp.Add(
                *(
                    background_dx(sp.diff(full_density, first_jets[i]), i)
                    for i in DIRECTIONS
                )
            )
        )

    background_substitutions = {
        item: sp.Integer(0) for item in ACTUAL_WAVE_SYMBOLS
    }

    def at_background(expression: sp.Expr) -> sp.Expr:
        value = sp.expand(
            expression.subs(background_substitutions, simultaneous=True)
        )
        return first_shape_series(value)

    u_force = tuple(
        at_background(background_euler(u[a], grad_u[a])) for a in DIRECTIONS
    )
    theta_force = at_background(background_euler(theta, grad_theta))
    e_force = at_background(background_euler(e_W, grad_e))
    body = sp.Tuple(
        sp.Tuple(Str("U"), sp.Tuple(*u_force)),
        sp.Tuple(Str("THETA"), theta_force),
        sp.Tuple(Str("E_W"), e_force),
    )
    face_rows = []
    for face in FACES:
        face_rows.append(
            sp.Tuple(
                sp.Integer(face),
                sp.Tuple(
                    sp.Integer(0),
                    sp.Integer(0),
                    sp.Integer(0),
                    sp.expand(face * e_force / W0),
                ),
            )
        )
    return sp.Tuple(
        sp.Tuple(Str("BODY_FORCE"), body),
        sp.Tuple(Str("PER_FACE_TRACTION"), sp.Tuple(*face_rows)),
    )


def admissibility_support() -> object:
    return sp.Tuple(
        sp.Tuple(
            Str("BODY_FORCE"),
            sp.Tuple(
                sp.Tuple(Str("U"), sp.Tuple(*support_body_u)),
                sp.Tuple(Str("THETA"), support_body_theta),
                sp.Tuple(Str("E_W"), support_body_e),
            ),
        ),
        sp.Tuple(
            Str("PER_FACE_TRACTION"),
            sp.Tuple(
                sp.Tuple(sp.Integer(1), sp.Tuple(*support_face[1])),
                sp.Tuple(sp.Integer(-1), sp.Tuple(*support_face[-1])),
            ),
        ),
    )


ENERGY_PRIMARY_CASES: dict[object, EnergyBuild] = {}
OPERATOR_PRIMARY_CASES: dict[object, object] = {}
ORIGINS_PRIMARY_CASES: dict[object, object] = {}
MU_PRIMARY_CASES: dict[object, object] = {}
KERNEL_PRIMARY_CASES: dict[object, object] = {}
KERNEL_ORIGINS_PRIMARY_CASES: dict[object, object] = {}
ADMISSIBILITY_OPERATOR_CASES: dict[object, object] = {}
ADMISSIBILITY_SUPPORT_CASES: dict[object, object] = {}
ADMISSIBILITY_RESIDUAL_CASES: dict[object, object] = {}


def emit_primary(quantity: str, raw: Mapping[object, object], key: str) -> None:
    dimensions = {case: dimension_object(value) for case, value in raw.items()}
    RAW_PRIMARY[quantity] = raw
    DIMENSION_PRIMARY[quantity] = dimensions
    emit(quantity, case_payload(raw, dimensions), key=key, export=True)


def task_energy_basis() -> None:
    variable = {}
    counts = {}
    new_counts = {}
    new = {}
    omissions = {}
    for branch in BRANCHES:
        build = construct_energy(branch)
        ENERGY_PRIMARY_CASES[branch] = build
        variable[branch] = retained_energy_basis(build.basis)
        counts[branch] = build.count
        new_counts[branch] = build.new_counts
        new[branch] = retained_grade(build.new_invariants)
        omissions[branch] = retained_grade(build.omissions)
    emit_primary("ENERGY_BASIS_VARIABLE", variable, "energy_basis_variable")
    emit_primary("ENERGY_BASIS_COUNT", counts, "energy_basis_count")
    emit_primary(
        "ENERGY_BASIS_NEW_INVARIANT_COUNT",
        new_counts,
        "energy_basis_new_invariant_count",
    )
    emit_primary("ENERGY_BASIS_NEW_INVARIANTS", new, "energy_basis_new_invariants")
    emit_primary("ENERGY_BASIS_OMISSIONS", omissions, "energy_basis_omissions")


def task_slab_operator() -> None:
    for branch in BRANCHES:
        for representative in DENSITY_REPS:
            case = (branch, representative)
            operator, origins, mu_theta_value = build_operator(branch, representative)
            OPERATOR_PRIMARY_CASES[case] = retained_grade(operator)
            ORIGINS_PRIMARY_CASES[case] = retained_grade(origins)
            MU_PRIMARY_CASES[case] = retained_grade(mu_theta_value)
    emit_primary("SLAB_OPERATOR", OPERATOR_PRIMARY_CASES, "slab_operator")
    emit_primary("SLAB_OPERATOR_TERM_ORIGINS", ORIGINS_PRIMARY_CASES, "slab_operator_term_origins")
    emit_primary("MU_THETA_OPERATOR", MU_PRIMARY_CASES, "mu_theta_operator")


def task_coupling_kernel() -> None:
    for branch in BRANCHES:
        for representative in DENSITY_REPS:
            case = (branch, representative)
            kernel, origins = build_kernel(branch, representative)
            KERNEL_PRIMARY_CASES[case] = retained_grade(kernel)
            KERNEL_ORIGINS_PRIMARY_CASES[case] = retained_grade(origins)
    emit_primary("COUPLING_KERNEL", KERNEL_PRIMARY_CASES, "coupling_kernel")
    emit_primary(
        "COUPLING_KERNEL_TERM_ORIGINS",
        KERNEL_ORIGINS_PRIMARY_CASES,
        "coupling_kernel_term_origins",
    )


def task_admissibility() -> None:
    support = admissibility_support()
    for branch in BRANCHES:
        for representative in DENSITY_REPS:
            case = (branch, representative)
            operator = retained_grade(admissibility_operator(branch, representative))
            ADMISSIBILITY_OPERATOR_CASES[case] = operator
            ADMISSIBILITY_SUPPORT_CASES[case] = retained_grade(support)
            ADMISSIBILITY_RESIDUAL_CASES[case] = object_difference(
                operator, retained_grade(support)
            )
    emit_primary(
        "ADMISSIBILITY_OPERATOR_OPERAND",
        ADMISSIBILITY_OPERATOR_CASES,
        "admissibility_operator_operand",
    )
    emit_primary(
        "ADMISSIBILITY_SUPPORT_OPERAND",
        ADMISSIBILITY_SUPPORT_CASES,
        "admissibility_support_operand",
    )
    emit_primary(
        "ADMISSIBILITY_RESIDUAL",
        ADMISSIBILITY_RESIDUAL_CASES,
        "admissibility_residual",
    )


def eta_denominator_signature(
    expression: sp.Expr,
) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    substituted = expression.subs(PROFILE_GRADE_SUBS, simultaneous=True)
    factors = {
        (node.base, node.exp)
        for node in sp.preorder_traversal(substituted)
        if isinstance(node, sp.Pow)
        and node.base.has(eta_bg)
        and node.exp.is_negative is True
    }
    return tuple(sorted(factors, key=sp.default_sort_key))


def projection_form_samples(value: object) -> tuple[tuple[object, sp.Expr], ...]:
    groups: dict[tuple[tuple[sp.Expr, sp.Expr], ...], list[sp.Expr]] = {}
    for scalar in object_scalars(value):
        if isinstance(scalar, sp.Expr):
            groups.setdefault(eta_denominator_signature(scalar), []).append(scalar)
    samples = []
    for signature, expressions in groups.items():
        eligible = [
            expression
            for expression in expressions
            if not expression.is_Atom
            and expression.has(
                eta_bg,
                sigma_W,
                W_bg,
                mu_R_bg,
                rho4_rho4,
                rhobr_rho4,
                rho4_rhobr,
                rhobr_rhobr,
            )
        ]
        candidates = eligible or [
            expression for expression in expressions if not expression.is_Atom
        ]
        candidates = candidates or expressions
        sample = min(candidates, key=lambda item: (sp.count_ops(item), sp.default_sort_key(item)))
        displayed_signature = (
            sp.Tuple(*(sp.Tuple(base, exponent) for base, exponent in signature))
            if signature
            else Str("NO_ETA_DENOMINATOR")
        )
        samples.append((displayed_signature, sample))
    return tuple(samples)


def exact_structural_residual(left: object, right: object) -> object:
    if sp.srepr(casify(left)) == sp.srepr(casify(right)):
        return sp.Integer(0)
    return structural_residual(left, right)


def task_projection_equivalence() -> None:
    samples = []
    for branch in BRANCHES:
        for representative in DENSITY_REPS:
            for route in ("EULERIAN", "MATERIAL"):
                objects = (
                    (
                        "SLAB_OPERATOR",
                        build_operator(branch, representative, route)[0],
                    ),
                    (
                        "COUPLING_KERNEL",
                        build_kernel(
                            branch,
                            representative,
                            route,
                            include_term_origins=False,
                        )[0],
                    ),
                    (
                        "ADMISSIBILITY_OPERATOR",
                        admissibility_operator(branch, representative, route),
                    ),
                )
                for object_name, value in objects:
                    for signature, expression in projection_form_samples(value):
                        samples.append(
                            (
                                object_name,
                                branch,
                                representative,
                                route,
                                signature,
                                expression,
                            )
                        )

    sample_expressions = tuple(row[-1] for row in samples)
    parallel_values = project_scalars_ordered(sample_expressions, allow_pool=True)
    reference_values = tuple(
        first_shape_series_reference(expression) for expression in sample_expressions
    )
    residual_rows = []
    operand_rows = []
    for sample, reference_value, parallel_value in zip(
        samples, reference_values, parallel_values
    ):
        key = sample[:-1]
        residual = exact_structural_residual(reference_value, parallel_value)
        operand_rows.append(
            sp.Tuple(
                *(Str(str(item)) if isinstance(item, str) else item for item in key),
                sample[-1],
                reference_value,
                parallel_value,
            )
        )
        residual_rows.append(
            sp.Tuple(
                *(Str(str(item)) if isinstance(item, str) else item for item in key),
                residual,
            )
        )

    failure_input = sp.Order(eta_bg)
    failure_exception = Str("NO_FAILURE")
    try:
        sp.series(
            failure_input.subs(PROFILE_GRADE_SUBS, simultaneous=True),
            eta_bg,
            0,
            2,
        ).removeO()
    except (NotImplementedError, ValueError, TypeError) as exc:
        failure_exception = Str(type(exc).__name__)
    failure_parallel = project_scalars_ordered((failure_input,), allow_pool=True)[0]
    failure_reference = first_shape_series_reference(failure_input)
    failure_residual = exact_structural_residual(
        failure_reference, failure_parallel
    )

    operator_raw = build_operator("LAB_HELD", "RHO4_CONSTANT", "EULERIAN")[0]
    kernel_raw = build_kernel("LAB_HELD", "RHO4_CONSTANT", "EULERIAN")[0]
    operator_serial = retained_grade_serial(operator_raw)
    operator_parallel = retained_grade_parallel(operator_raw)
    kernel_serial = retained_grade_serial(kernel_raw)
    kernel_parallel = retained_grade_parallel(kernel_raw)
    operator_residual = exact_structural_residual(
        operator_serial, operator_parallel
    )
    kernel_residual = exact_structural_residual(kernel_serial, kernel_parallel)

    emit("PROJECTION_EQUIVALENCE_OPERAND", sp.Tuple(*operand_rows))
    emit("PROJECTION_EQUIVALENCE_RESIDUAL", sp.Tuple(*residual_rows))
    emit(
        "PROJECTION_FAILURE_PATH_OPERAND",
        sp.Tuple(
            failure_input,
            failure_exception,
            failure_reference,
            failure_parallel,
        ),
    )
    emit("PROJECTION_FAILURE_PATH_RESIDUAL", failure_residual)
    emit("PROJECTION_OPERATOR_ORDERING_RESIDUAL", operator_residual)
    emit("PROJECTION_KERNEL_ORDERING_RESIDUAL", kernel_residual)

    if any(row[-1] != 0 for row in residual_rows):
        raise AssertionError("scalar projection equivalence changed")
    if failure_exception == Str("NO_FAILURE") or failure_residual != 0:
        raise AssertionError("series failure fallback equivalence changed")
    if operator_residual != 0 or kernel_residual != 0:
        raise AssertionError("parallel projection reassembly changed")


def task_rep_invariance() -> None:
    eulerian = {}
    material = {}
    for branch in BRANCHES:
        for representative in DENSITY_REPS:
            for object_name in ("SLAB_OPERATOR", "COUPLING_KERNEL", "ADMISSIBILITY_OPERATOR"):
                case = (object_name, branch, representative)
                if object_name == "SLAB_OPERATOR":
                    eulerian[case] = retained_grade(
                        build_operator(branch, representative, "EULERIAN")[0]
                    )
                    material[case] = retained_grade(
                        build_operator(branch, representative, "MATERIAL")[0]
                    )
                elif object_name == "COUPLING_KERNEL":
                    eulerian[case] = retained_grade(
                        build_kernel(
                            branch,
                            representative,
                            "EULERIAN",
                            include_term_origins=False,
                        )[0]
                    )
                    material[case] = retained_grade(
                        build_kernel(
                            branch,
                            representative,
                            "MATERIAL",
                            include_term_origins=False,
                        )[0]
                    )
                else:
                    eulerian[case] = retained_grade(
                        admissibility_operator(branch, representative, "EULERIAN")
                    )
                    material[case] = retained_grade(
                        admissibility_operator(branch, representative, "MATERIAL")
                    )
    residual = object_difference(eulerian, material)
    emit("REP_INVARIANCE_EULERIAN_OPERAND", case_payload(eulerian))
    emit("REP_INVARIANCE_MATERIAL_OPERAND", case_payload(material))
    emit("REP_INVARIANCE_RESIDUAL", case_payload(residual))


def task_independence() -> None:
    baseline = {}
    corrupted = {}

    def record_if_source_present(
        object_name: str,
        branch: str,
        representative: str,
        base_value: object,
        corrupted_value: object,
    ) -> None:
        if sp.srepr(casify(base_value)) == sp.srepr(casify(corrupted_value)):
            return
        key = (object_name, branch, representative)
        baseline[key] = base_value
        corrupted[key] = corrupted_value

    for branch in BRANCHES:
        for representative in DENSITY_REPS:
            operator_base = retained_grade(
                build_operator(branch, representative, "EULERIAN")[0]
            )
            operator_corrupt = retained_grade(
                build_operator(
                    branch,
                    representative,
                    "EULERIAN",
                    corrupt_material_constraint=True,
                )[0]
            )
            record_if_source_present(
                "SLAB_OPERATOR",
                branch,
                representative,
                operator_base,
                operator_corrupt,
            )
            kernel_base = retained_grade(
                build_kernel(
                    branch,
                    representative,
                    "EULERIAN",
                    include_term_origins=False,
                )[0]
            )
            kernel_corrupt = retained_grade(
                build_kernel(
                    branch,
                    representative,
                    "EULERIAN",
                    include_term_origins=False,
                    corrupt_material_constraint=True,
                )[0]
            )
            record_if_source_present(
                "COUPLING_KERNEL",
                branch,
                representative,
                kernel_base,
                kernel_corrupt,
            )
            admissibility_base = retained_grade(
                admissibility_operator(branch, representative, "EULERIAN")
            )
            admissibility_corrupt = retained_grade(
                admissibility_operator(
                    branch,
                    representative,
                    "EULERIAN",
                    corrupt_material_constraint=True,
                )
            )
            record_if_source_present(
                "ADMISSIBILITY_OPERATOR",
                branch,
                representative,
                admissibility_base,
                admissibility_corrupt,
            )
    residual = object_difference(baseline, corrupted)
    emit("CONTROL_INDEPENDENCE_BASE_OPERAND", case_payload(baseline))
    emit("CONTROL_INDEPENDENCE_CORRUPTED_OPERAND", case_payload(corrupted))
    emit("CONTROL_INDEPENDENCE_RESIDUAL", case_payload(residual))


def task_form_control() -> None:
    baseline = {}
    ablated = {}
    for source in ("W_BG", "MU_R_BG"):
        for direction in DIRECTIONS:
            for branch in BRANCHES:
                for representative in DENSITY_REPS:
                    cases = (
                        (
                            "ENERGY_BASIS",
                            retained_energy_basis(construct_energy(branch).basis),
                            retained_energy_basis(
                                construct_energy(branch, source, direction).basis
                            ),
                        ),
                        (
                            "SLAB_OPERATOR",
                            retained_grade(build_operator(branch, representative)[0]),
                            retained_grade(
                                build_operator(
                                    branch,
                                    representative,
                                    "EULERIAN",
                                    source,
                                    direction,
                                )[0]
                            ),
                        ),
                        (
                            "COUPLING_KERNEL",
                            retained_grade(build_kernel(branch, representative)[0]),
                            retained_grade(
                                build_kernel(
                                    branch,
                                    representative,
                                    "EULERIAN",
                                    source,
                                    direction,
                                    include_term_origins=False,
                                )[0]
                            ),
                        ),
                    )
                    for object_name, base_value, ablated_value in cases:
                        key = (object_name, branch, representative, source, direction + 1)
                        baseline[key] = base_value
                        ablated[key] = ablated_value
    residual = object_difference(baseline, ablated)
    emit("CONTROL_FORM_BASE_OPERAND", case_payload(baseline))
    emit("CONTROL_FORM_ABLATED_OPERAND", case_payload(ablated))
    emit("CONTROL_FORM_RESIDUAL", case_payload(residual))


def source_new_rows(build: EnergyBuild, source: str) -> sp.Tuple:
    return sp.Tuple(*(row for row in build.new_invariants if str(row[0]) == source))


def source_new_count(build: EnergyBuild, source: str) -> sp.Integer:
    return sp.Integer(len(source_new_rows(build, source)))


def source_new_label_set(build: EnergyBuild, source: str) -> sp.FiniteSet:
    return sp.FiniteSet(
        *(row[1] for row in source_new_rows(build, source))
    )


def task_hessian_freeze_control() -> None:
    target = sp.Integer(26)
    counts = {}
    count_targets = {}
    count_residuals = {}
    live_rows = {}
    frozen_rows = {}
    committed_rows = {}
    row_residuals = {}
    form_residuals = {}
    for branch in BRANCHES:
        live_energy = construct_energy(
            branch, background_depth=STRONG_ROW_JET_DEPTH
        )
        frozen_energy = construct_energy(
            branch, background_depth=HESSIAN_FREEZE_DEPTH
        )
        counts[branch] = frozen_energy.count
        count_targets[branch] = target
        count_residuals[branch] = sp.expand(frozen_energy.count - target)
        for representative in DENSITY_REPS:
            for route in ("EULERIAN", "MATERIAL"):
                case = (branch, representative, route)
                live_value = retained_grade(
                    live_strong_rows(
                        live_energy.density,
                        branch,
                        representative,
                        route,
                        STRONG_ROW_JET_DEPTH,
                    )
                )
                frozen_value = retained_grade(
                    live_strong_rows(
                        frozen_energy.density,
                        branch,
                        representative,
                        route,
                        HESSIAN_FREEZE_DEPTH,
                    )
                )
                committed_value = retained_grade(
                    committed_strong_rows(
                        frozen_energy.density,
                        branch,
                        representative,
                        route,
                    )
                )
                live_rows[case] = live_value
                frozen_rows[case] = frozen_value
                committed_rows[case] = committed_value
                row_residuals[case] = structural_residual(
                    frozen_value, committed_value
                )
                form_residuals[case] = structural_residual(
                    live_value, frozen_value
                )
    emit("HESSIAN_FREEZE_COUNT_OPERAND", case_payload(counts))
    emit("HESSIAN_FREEZE_COUNT_REFERENCE_OPERAND", case_payload(count_targets))
    emit("HESSIAN_FREEZE_COUNT_RESIDUAL", case_payload(count_residuals))
    emit("HESSIAN_FREEZE_STRONG_ROW_OPERAND", case_payload(frozen_rows))
    emit(
        "HESSIAN_FREEZE_COMMITTED_STRONG_ROW_OPERAND",
        case_payload(committed_rows),
    )
    emit("HESSIAN_FREEZE_STRONG_ROW_RESIDUAL", case_payload(row_residuals))
    emit("HESSIAN_FORM_LIVE_STRONG_ROW_OPERAND", case_payload(live_rows))
    emit("HESSIAN_FORM_STRONG_ROW_RESIDUAL", case_payload(form_residuals))
    if any(residual != 0 for residual in count_residuals.values()):
        raise AssertionError("Hessian-freeze count regression changed")


def task_coefficient_control() -> None:
    baseline = {}
    scaled = {}
    residual = {}
    for source in ("W_BG", "MU_R_BG"):
        candidates = enumerate_new_candidates(tuple(bg))
        expressions = tuple(expression for _, expression in candidates)
        scaled_expressions = tuple(
            sp.expand(
                expression.subs(
                    {bg[i]: 2 * bg[i] for i in DIRECTIONS}, simultaneous=True
                )
            )
            for expression in expressions
        )
        baseline_signatures = basis_euler_signatures(
            expressions,
            basis_fields,
            background_first_jets=bg,
            background_second_jets=basis_background_second,
            background_depth=STRONG_ROW_JET_DEPTH,
        )
        scaled_signatures = basis_euler_signatures(
            scaled_expressions,
            basis_fields,
            background_first_jets=bg,
            background_second_jets=basis_background_second,
            background_depth=STRONG_ROW_JET_DEPTH,
        )
        baseline_selected, _ = quotient_independent_indices(
            expressions, baseline_signatures
        )
        scaled_selected, _ = quotient_independent_indices(
            scaled_expressions, scaled_signatures
        )
        baseline[source] = sp.Tuple(
            sp.Integer(len(baseline_selected)),
            sp.Tuple(*(sp.Integer(index + 1) for index in baseline_selected)),
            sp.Tuple(*(sp.Tuple(*row) for row in baseline_signatures)),
        )
        scaled[source] = sp.Tuple(
            sp.Integer(len(scaled_selected)),
            sp.Tuple(*(sp.Integer(index + 1) for index in scaled_selected)),
            sp.Tuple(*(sp.Tuple(*row) for row in scaled_signatures)),
        )
        residual[source] = sp.Tuple(
            sp.Integer(len(scaled_selected) - len(baseline_selected)),
            object_difference(
                sp.Tuple(*(sp.Integer(index + 1) for index in scaled_selected)),
                sp.Tuple(*(sp.Integer(index + 1) for index in baseline_selected)),
            ),
            object_difference(scaled_signatures, baseline_signatures),
        )
    emit("COEFFICIENT_CONTROL_BASE_OPERAND", case_payload(baseline))
    emit("COEFFICIENT_CONTROL_SCALED_OPERAND", case_payload(scaled))
    emit("COEFFICIENT_CONTROL_RESIDUAL", case_payload(residual))


def task_full_source_jet_zero_control() -> None:
    live = {}
    frozen = {}
    zeroed = {}
    inclusion_residual = {}
    zero_substitutions = {}
    for branch in BRANCHES:
        for source in ("W_BG", "MU_R_BG"):
            case = (branch, source)
            live_build = construct_energy(
                branch, background_depth=STRONG_ROW_JET_DEPTH
            )
            frozen_build = construct_energy(
                branch, background_depth=HESSIAN_FREEZE_DEPTH
            )
            zero_build = construct_energy(
                branch,
                background_depth=STRONG_ROW_JET_DEPTH,
                full_zero_source=source,
            )
            live_set = source_new_label_set(live_build, source)
            frozen_set = source_new_label_set(frozen_build, source)
            zero_set = source_new_label_set(zero_build, source)
            live[case] = sp.Tuple(
                source_new_count(live_build, source), live_set
            )
            frozen[case] = sp.Tuple(
                source_new_count(frozen_build, source), frozen_set
            )
            zeroed[case] = sp.Tuple(
                source_new_count(zero_build, source),
                zero_set,
                retained_grade(source_new_rows(zero_build, source)),
            )
            inclusion_residual[case] = sp.Tuple(
                sp.Complement(zero_set, frozen_set),
                sp.Complement(frozen_set, live_set),
                sp.Integer(source_new_count(zero_build, source)),
            )
            zero_substitutions[case] = sp.Tuple(
                *(
                    sp.Tuple(atom, sp.Integer(0))
                    for atom in dict.fromkeys(
                        (
                            *BACKGROUND_FIRST_JETS[source],
                            *(
                                profile_atom
                                for table in BACKGROUND_PROFILE_JETS[source].values()
                                for profile_atom in table.values()
                            ),
                        )
                    )
                )
            )
    emit("FULL_SOURCE_JET_ZERO_LIVE_OPERAND", case_payload(live))
    emit("FULL_SOURCE_JET_ZERO_FROZEN_OPERAND", case_payload(frozen))
    emit("FULL_SOURCE_JET_ZERO_OPERAND", case_payload(zeroed))
    emit("FULL_SOURCE_JET_ZERO_SUBSTITUTION", case_payload(zero_substitutions))
    emit("FULL_SOURCE_JET_ZERO_INCLUSION_RESIDUAL", case_payload(inclusion_residual))


def task_tower_depth_control() -> None:
    energy_chosen = {}
    energy_deeper = {}
    energy_shallower = {}
    pullback_chosen = {}
    pullback_deeper = {}
    pullback_shallower = {}
    strong_chosen = {}
    strong_deeper = {}
    strong_shallower = {}
    coupling_chosen = {}
    coupling_deeper = {}
    coupling_shallower = {}
    for branch in BRANCHES:
        chosen_build = construct_energy(
            branch, background_depth=STRONG_ROW_JET_DEPTH
        )
        deeper_build = construct_energy(
            branch, background_depth=COUPLING_JET_DEPTH
        )
        shallower_build = construct_energy(
            branch, background_depth=HESSIAN_FREEZE_DEPTH
        )
        energy_chosen[branch] = retained_energy_basis(chosen_build.basis)
        energy_deeper[branch] = retained_energy_basis(deeper_build.basis)
        energy_shallower[branch] = retained_energy_basis(shallower_build.basis)
        for representative in DENSITY_REPS:
            pullback_case = (branch, representative)
            pullback_chosen[pullback_case] = retained_grade(
                material_pullback(
                    chosen_build.density,
                    branch,
                    representative,
                    background_depth=STRONG_ROW_JET_DEPTH,
                )
            )
            pullback_deeper[pullback_case] = retained_grade(
                material_pullback(
                    deeper_build.density,
                    branch,
                    representative,
                    background_depth=COUPLING_JET_DEPTH,
                )
            )
            pullback_shallower[pullback_case] = retained_grade(
                material_pullback(
                    shallower_build.density,
                    branch,
                    representative,
                    background_depth=HESSIAN_FREEZE_DEPTH,
                )
            )
            case = (branch, representative)
            chosen_rows = retained_grade(
                live_strong_rows(
                    chosen_build.density,
                    branch,
                    representative,
                    "EULERIAN",
                    STRONG_ROW_JET_DEPTH,
                )
            )
            deeper_rows = retained_grade(
                live_strong_rows(
                    deeper_build.density,
                    branch,
                    representative,
                    "EULERIAN",
                    COUPLING_JET_DEPTH,
                )
            )
            shallower_rows = retained_grade(
                live_strong_rows(
                    shallower_build.density,
                    branch,
                    representative,
                    "EULERIAN",
                    HESSIAN_FREEZE_DEPTH,
                )
            )
            strong_chosen[case] = chosen_rows
            strong_deeper[case] = deeper_rows
            strong_shallower[case] = shallower_rows
            coupling_chosen[case] = retained_grade(
                coupling_outer_rows(chosen_rows, COUPLING_JET_DEPTH)
            )
            coupling_deeper[case] = retained_grade(
                coupling_outer_rows(chosen_rows, DEPTH_CONTROL_JET_DEPTH)
            )
            coupling_shallower[case] = retained_grade(
                coupling_outer_rows(chosen_rows, STRONG_ROW_JET_DEPTH)
            )
    depth_operands = sp.Tuple(
        sp.Tuple(Str("ENERGY_BASIS"), sp.Integer(STRONG_ROW_JET_DEPTH)),
        sp.Tuple(Str("MATERIAL_PULLBACK"), sp.Integer(STRONG_ROW_JET_DEPTH)),
        sp.Tuple(Str("STRONG_ROWS"), sp.Integer(STRONG_ROW_JET_DEPTH)),
        sp.Tuple(Str("COUPLING_CASCADE"), sp.Integer(COUPLING_JET_DEPTH)),
        sp.Tuple(Str("DEEPER_CONTROL"), sp.Integer(DEPTH_CONTROL_JET_DEPTH)),
    )
    deeper_residual = {
        "ENERGY_BASIS": structural_residual(energy_deeper, energy_chosen),
        "MATERIAL_PULLBACK": structural_residual(pullback_deeper, pullback_chosen),
        "STRONG_ROWS": structural_residual(strong_deeper, strong_chosen),
        "COUPLING_CASCADE": structural_residual(coupling_deeper, coupling_chosen),
    }
    shallower_residual = {
        "ENERGY_BASIS": structural_residual(energy_chosen, energy_shallower),
        "MATERIAL_PULLBACK": structural_residual(pullback_chosen, pullback_shallower),
        "STRONG_ROWS": structural_residual(strong_chosen, strong_shallower),
        "COUPLING_CASCADE": structural_residual(coupling_chosen, coupling_shallower),
    }
    emit("TOWER_DEPTH_OPERAND", depth_operands)
    emit("TOWER_DEPTH_DEEPER_RESIDUAL", case_payload(deeper_residual))
    emit("TOWER_DEPTH_SHALLOWER_RESIDUAL", case_payload(shallower_residual))


def uniformize(value: object) -> object:
    substitutions = {
        W_bg: W0,
        mu_R_bg: mu_R,
        eta_bg: sp.Integer(0),
        sigma_W: sp.Integer(0),
        **{grad_W[i]: sp.Integer(0) for i in DIRECTIONS},
        **{grad_mu[i]: sp.Integer(0) for i in DIRECTIONS},
        **{w1_grad[i]: sp.Integer(0) for i in DIRECTIONS},
        **{m1_grad[i]: sp.Integer(0) for i in DIRECTIONS},
        rhobr_rho4: rho_br,
        rhobr_rhobr: rho_br,
        rho4_rho4: rho_br / W0,
        rho4_rhobr: rho_br / W0,
    }
    return map_object(
        value,
        lambda expression: sp.expand(expression.subs(substitutions, simultaneous=True)),
    )


k_uniform = inherited_symbol("k", "COORDINATE", "uniform-regression wave number", -DIM_L)
omega_uniform = inherited_symbol("omega", "COORDINATE", "uniform-regression frequency", -DIM_T)


def uniform_transverse_dispersion() -> object:
    density = uniformize(construct_energy("LAB_HELD").density)
    amplitude = sp.Dummy("transverse_amplitude")
    substitutions = {item: sp.Integer(0) for item in ACTUAL_WAVE_SYMBOLS}
    substitutions[grad_u[1][0]] = sp.I * k_uniform * amplitude
    substitutions[u_tt[1]] = -omega_uniform**2 * amplitude
    operator, _ = operator_from_density(density, "RHOBR_CONSTANT")
    expanded = operator["U_BODY_BALANCE"][2][1][1]
    operand = sp.simplify(expanded.subs(substitutions, simultaneous=True) / (epsilon * amplitude))
    return sp.Eq(operand, 0, evaluate=False)


def task_uniform_control() -> None:
    s11cb = {}
    s11b = {}
    for branch in BRANCHES:
        for representative in DENSITY_REPS:
            for object_name in ("SLAB_OPERATOR", "COUPLING_KERNEL", "TRANSVERSE_DISPERSION"):
                key = (object_name, branch, representative)
                if object_name == "SLAB_OPERATOR":
                    s11cb[key] = uniformize(build_operator(branch, representative)[0])
                    s11b[key] = sp.Tuple(
                        INCOMING_LEDGER["inplane_eom"]["value"],
                        INCOMING_LEDGER["thickness_eom"]["value"],
                        INCOMING_LEDGER["constraint"]["value"],
                    )
                elif object_name == "COUPLING_KERNEL":
                    s11cb[key] = uniformize(build_kernel(branch, representative)[0])
                    s11b[key] = INCOMING_LEDGER["transverse_coupling"]["value"]
                else:
                    s11cb[key] = uniform_transverse_dispersion()
                    s11b[key] = INCOMING_LEDGER["transverse_dispersion"]["value"]
    residual = object_difference(s11cb, s11b)
    emit("UNIFORM_LIMIT_S11CB_OPERAND", case_payload(s11cb))
    emit("UNIFORM_LIMIT_S11B_OPERAND", case_payload(s11b))
    emit("UNIFORM_LIMIT_RESIDUAL", case_payload(residual))


def task_homogeneity() -> None:
    dimensions = {
        quantity: DIMENSION_PRIMARY[quantity]
        for quantity in (
            "ENERGY_BASIS_VARIABLE",
            "ENERGY_BASIS_COUNT",
            "ENERGY_BASIS_NEW_INVARIANT_COUNT",
            "ENERGY_BASIS_NEW_INVARIANTS",
            "ENERGY_BASIS_OMISSIONS",
            "SLAB_OPERATOR",
            "SLAB_OPERATOR_TERM_ORIGINS",
            "MU_THETA_OPERATOR",
            "COUPLING_KERNEL",
            "COUPLING_KERNEL_TERM_ORIGINS",
            "ADMISSIBILITY_OPERATOR_OPERAND",
            "ADMISSIBILITY_SUPPORT_OPERAND",
            "ADMISSIBILITY_RESIDUAL",
        )
    }
    used_new_coefficients = sorted(
        (
            symbol
            for value in RAW_PRIMARY.values()
            for symbol in atoms_of(value)
            if symbol.name.startswith("gamma_s11cb_")
        ),
        key=sp.default_sort_key,
    )
    corrupted_dimensions = dict(SYMBOL_DIMENSIONS)
    if not used_new_coefficients:
        raise RuntimeError("no constructed first-jet coefficient reached a primary object")
    corrupted_symbol = used_new_coefficients[0]
    corrupted_dimensions[corrupted_symbol] = dim_mul(
        SYMBOL_DIMENSIONS[corrupted_symbol], DIM_L
    )
    base = {
        quantity: {case: dimension_object(value) for case, value in raw.items()}
        for quantity, raw in RAW_PRIMARY.items()
    }
    control = {
        quantity: {
            case: dimension_object(value, corrupted_dimensions) for case, value in raw.items()
        }
        for quantity, raw in RAW_PRIMARY.items()
    }
    residual = object_difference(base, control)
    emit("DIMENSIONS", casify(dimensions))
    emit("HOMOGENEITY_BASE_OPERAND", casify(base))
    emit("HOMOGENEITY_CONTROL_OPERAND", casify(control))
    emit("HOMOGENEITY_RESIDUAL", casify(residual))


def total_compare(left: object, right: object) -> tuple[str, object]:
    """F9's three-valued comparison, total over imported row shapes."""
    if isinstance(left, bool) or isinstance(right, bool):
        if isinstance(left, bool) and isinstance(right, bool):
            status = PROVED_EQUAL if left is right else PROVED_DIFFERENT
            return status, sp.Tuple(casify(left), casify(right))
        return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
    if isinstance(left, str) or isinstance(right, str):
        if isinstance(left, str) and isinstance(right, str):
            status = PROVED_EQUAL if left == right else PROVED_DIFFERENT
            return status, sp.Tuple(Str(left), Str(right))
        return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        if not isinstance(left, Mapping) or not isinstance(right, Mapping) or set(left) != set(right):
            return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
        ordered = sorted(left, key=str)
        comparisons = [total_compare(left[key], right[key]) for key in ordered]
        statuses = [status for status, _ in comparisons]
        detail = sp.Tuple(
            *(
                sp.Tuple(Str(str(key)), Str(status), casify(operand))
                for key, (status, operand) in zip(ordered, comparisons)
            )
        )
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
        detail = sp.Tuple(
            *(sp.Tuple(Str(status), casify(operand)) for status, operand in comparisons)
        )
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
        return (PROVED_EQUAL if equal else PROVED_DIFFERENT), sp.Tuple(left, right)
    if isinstance(left, Str) or isinstance(right, Str):
        equal = isinstance(left, Str) and isinstance(right, Str) and sp.srepr(left) == sp.srepr(right)
        return (PROVED_EQUAL if equal else PROVED_DIFFERENT), sp.Tuple(left, right)
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
        candidates.append(
            CandidateRow(
                free_symbol.name,
                free_symbol,
                metadata["class"],
                "FREE_SYMBOL_POPULATION",
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
        "step": "S11c-b",
    }
    if candidate.description is not None:
        record["description"] = candidate.description
    return record


def merged_export() -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    candidates = list(EMITTER.export_candidates)
    add_free_symbol_candidates(candidates)
    candidate_keys = sp.Tuple(*(Str(candidate.key) for candidate in candidates))
    imported_matching_keys = sp.Tuple(
        *(Str(candidate.key) for candidate in candidates if candidate.key in INCOMING_LEDGER)
    )
    merged = {name: dict(record) for name, record in INCOMING_LEDGER.items()}
    routing = []
    f9c_rows = []
    seen: set[str] = set()
    for candidate in candidates:
        if candidate.key in seen:
            raise RuntimeError(f"two S11c-b candidates have the same key {candidate.key!r}")
        seen.add(candidate.key)
        record = make_record(candidate)
        imported_record = INCOMING_LEDGER.get(candidate.key)
        if imported_record is None:
            write_key = candidate.key
            status = "ABSENT"
            comparison_operand = sp.Tuple(record["value"], sp.Tuple())
        else:
            value_status, comparison_operand = total_compare(
                record["value"], imported_record.get("value")
            )
            class_status = (
                PROVED_EQUAL
                if record["class"] == imported_record.get("class")
                else PROVED_DIFFERENT
            )
            if value_status == PROVED_EQUAL and class_status == PROVED_EQUAL:
                write_key = candidate.key
                status = PROVED_EQUAL
                if "dimension_key" in imported_record:
                    record["dimension_key"] = imported_record["dimension_key"]
                record["f9_operands"] = sp.Tuple(
                    record["value"], imported_record.get("value")
                )
                prior_steps = list(imported_record.get("corroborated_steps", ()))
                if not prior_steps and imported_record.get("step"):
                    prior_steps.append(imported_record["step"])
                prior_steps.append("S11c-b")
                record["corroborated_steps"] = tuple(dict.fromkeys(prior_steps))
            else:
                write_key = f"s11c_b_{candidate.key}"
                status = value_status if value_status != PROVED_EQUAL else PROVED_DIFFERENT
                f9c_rows.append(
                    sp.Tuple(
                        Str(candidate.key),
                        Str(write_key),
                        Str(status),
                        casify(comparison_operand),
                    )
                )
                issue(f"F9c write {write_key} for {candidate.key}: {status}")
        routing.append(
            sp.Tuple(
                Str(candidate.source_tag),
                Str(candidate.key),
                Str(write_key),
                Str(status),
                casify(comparison_operand),
            )
        )
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
        "# S11c_b_exports.py -- GENERATED by S11c_b_brane_operator_sympy_audit.py. Do not edit.",
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
        f"    'S11c_b_brane_operator_sympy_audit.py': {sha256(SCRIPT_PATH)!r},",
        f"    'S11c_a_exports.py': {sha256(SCRIPT_DIR / 'S11c_a_exports.py')!r},",
        f"    'S11c_b_SHARED_PHYSICS.md': {sha256(DIRECTIVE_PATH)!r},",
        "})",
        "",
        "_LEDGER = {",
    ]
    for name in sorted(ledger):
        lines.append(f"    {name!r}: {record_source(ledger[name])},")
    lines.extend(
        [
            "}",
            "LEDGER = MappingProxyType({",
            "    name: MappingProxyType(record) for name, record in _LEDGER.items()",
            "})",
            "del _LEDGER",
            "",
        ]
    )
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
        residual = (
            sp.Integer(0)
            if status == PROVED_EQUAL
            else object_difference(live_record["value"], restored)
        )
        rows.append(
            sp.Tuple(
                Str(name),
                Str(str(live_record["display"])),
                Str(status),
                casify(residual),
                casify(operands),
            )
        )
        if status != PROVED_EQUAL:
            failures.append(name)
    emit("EXPORT_ROUNDTRIP", sp.Tuple(*rows), local=True)
    if failures:
        raise RuntimeError(
            "export reconstruction was not proved equal for " + ", ".join(failures)
        )


def remove_stale_export() -> None:
    if EXPORT_PATH.exists():
        EXPORT_PATH.unlink()


def run_task(name: str, function) -> list[str]:
    global CURRENT_TASK
    task_start = time.monotonic()
    CURRENT_TASK = name
    try:
        function()
        return [name]
    except Exception as exc:
        OPERATIONAL_FAILURES.append(f"{name}: {type(exc).__name__}: {exc}")
        return []
    finally:
        TASK_TIMINGS[name] = time.monotonic() - task_start
        CURRENT_TASK = None


def run() -> None:
    start = time.monotonic()
    TASK_TIMINGS.clear()
    completed: list[str] = []

    completed += run_task("ENERGY_BASIS", task_energy_basis)
    completed += run_task("SLAB_OPERATOR", task_slab_operator)
    completed += run_task("COUPLING_KERNEL", task_coupling_kernel)
    completed += run_task("ADMISSIBILITY", task_admissibility)

    # S11CB_PRIMARIES_ONLY skips the build-heavy tasks (PROJECTION_EQUIVALENCE and the
    # control tasks both rebuild/sample the expensive MATERIAL kernels — the deferred
    # cost); the primaries are emitted and S11c_b_exports.py is still published (gated on
    # PRIMARY_TASKS below). lever-C equivalence is verified out-of-band; see _measurements.
    if not os.environ.get("S11CB_PRIMARIES_ONLY"):
        completed += run_task("PROJECTION_EQUIVALENCE", task_projection_equivalence)
        completed += run_task("REP_INVARIANCE", task_rep_invariance)
        completed += run_task("INDEPENDENCE", task_independence)
        completed += run_task("FORM", task_form_control)
        completed += run_task("HESSIAN_FREEZE", task_hessian_freeze_control)
        completed += run_task("COEFFICIENT", task_coefficient_control)
        completed += run_task("JET_ZERO", task_full_source_jet_zero_control)
        completed += run_task("TOWER_DEPTH", task_tower_depth_control)
        completed += run_task("UNIFORM", task_uniform_control)
        completed += run_task("HOMOGENEITY", task_homogeneity)
    close_projection_pool()

    primary_complete = all(task in completed for task in PRIMARY_TASKS)
    f9c_report = sp.Tuple()
    export_start = time.monotonic()
    if primary_complete:
        try:
            ledger, diagnostics = merged_export()
            emit(
                "EXPORT_CANDIDATE_KEY_OPERANDS",
                sp.Tuple(
                    diagnostics["candidate_keys"],
                    diagnostics["imported_matching_keys"],
                    diagnostics["routing"],
                ),
                local=True,
            )
            f9c_report = diagnostics["f9c"]
            emit("EXPORT_F9C_WRITES", f9c_report, local=True)
            publish_export(ledger)
        except Exception as exc:
            OPERATIONAL_FAILURES.append(f"EXPORT: {type(exc).__name__}: {exc}")
            issue(f"export operational failure: {type(exc).__name__}: {exc}")
            remove_stale_export()
    else:
        remove_stale_export()
        issue("S11c_b_exports.py not published because a primary task did not complete")
    TASK_TIMINGS["EXPORT"] = time.monotonic() - export_start

    emit(
        "OPERATIONAL_EXCEPTIONS",
        sp.Tuple(*(Str(item) for item in OPERATIONAL_FAILURES)),
        local=True,
    )
    skipped = tuple(task for task in (*PRIMARY_TASKS, *CONTROL_TASKS) if task not in completed)
    run_record = sp.Tuple(*(Str(task) for task in completed))
    skipped_record = sp.Tuple(*(Str(task) for task in skipped))
    emit("RUN_TASKS", run_record, local=True)
    emit("SKIPPED_TASKS", skipped_record, local=True)
    timing_record = sp.Tuple(
        *(
            sp.Tuple(Str(task), sp.Float(TASK_TIMINGS[task], 8))
            for task in (*PRIMARY_TASKS, *CONTROL_TASKS, "EXPORT")
            if task in TASK_TIMINGS
        ),
        sp.Tuple(Str("TOTAL"), sp.Float(time.monotonic() - start, 8)),
    )
    emit("TASK_TIMING_SECONDS", timing_record, local=True)

    local_names = list(
        dict.fromkeys(
            (
                *EMITTER.local_tags,
                "PY_S11CB_LOCAL_TAGS",
                "PY_S11CB_LOCAL_SECTION8_REPORT",
            )
        )
    )
    emit("TAGS", sp.Tuple(*(Str(name) for name in local_names)), local=True)
    runtime = sp.Float(time.monotonic() - start, 8)
    export_lines = (
        len(EXPORT_PATH.read_text(encoding="utf-8").splitlines())
        if EXPORT_PATH.exists()
        else 0
    )
    report = sp.Tuple(
        sp.Tuple(
            Str("FILES_WRITTEN"),
            sp.Tuple(
                Str(str(SCRIPT_PATH)),
                Str(str(EXPORT_PATH)) if EXPORT_PATH.exists() else Str("EXPORT_NOT_PUBLISHED"),
            ),
        ),
        sp.Tuple(
            Str("SCRIPT_LINES"),
            sp.Integer(len(SCRIPT_PATH.read_text(encoding="utf-8").splitlines())),
        ),
        sp.Tuple(Str("EXPORT_LINES"), sp.Integer(export_lines)),
        sp.Tuple(Str("EMITTED_TAGS"), sp.Integer(EMITTER.count + 1)),
        sp.Tuple(Str("RUN_TASKS"), run_record),
        sp.Tuple(Str("SKIPPED_TASKS"), skipped_record),
        sp.Tuple(Str("RUNTIME_SECONDS"), runtime),
        sp.Tuple(
            Str("TAG_NAMES"),
            sp.Tuple(
                *(Str(name) for name in (*EMITTER.values.keys(), "PY_S11CB_LOCAL_SECTION8_REPORT"))
            ),
        ),
        sp.Tuple(
            Str("ISSUES_OR_NONCOMPUTABLE"),
            sp.Tuple(*(Str(item) for item in (*ISSUES, *OPERATIONAL_FAILURES))),
        ),
        sp.Tuple(Str("F9C_WRITES"), f9c_report),
        sp.Tuple(
            Str("SUPPLIED_UNFALSIFIABLE"),
            Str("SECTIONS_1_TO_2_AND_S11CA_SHAPE_DERIVATIVE_SUBSTRATE"),
        ),
    )
    emit("SECTION8_REPORT", report, local=True)
    if OPERATIONAL_FAILURES:
        raise RuntimeError("; ".join(OPERATIONAL_FAILURES))


if __name__ == "__main__":
    run()
