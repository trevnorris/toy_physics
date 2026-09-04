#!/usr/bin/env python3
"""SymPy audit for S11c-c1 curved-bulk closure.

The sole physics input is ``directives/S11c_c1_SHARED_PHYSICS.md``.  The
inherited model is read through ``ledger_fold.load_model`` over the atomic
S11c-b base.  This program streams every requested CAS object and writes only
the S11c-c1 bind-closure delta ``S11c_c1_exports.py``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
import hashlib
import sys
import time

import sympy as sp
from sympy.core.relational import Relational
from sympy.core.symbol import Str
from sympy.logic.boolalg import Boolean
from sympy.printing.str import StrPrinter


SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
DIRECTIVE_PATH = SCRIPT_DIR.parent / "directives" / "S11c_c1_SHARED_PHYSICS.md"
BASE_PATH = SCRIPT_DIR / "S11c_b_exports.py"
EXPORT_PATH = SCRIPT_DIR / "S11c_c1_exports.py"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from ledger_fold import (  # noqa: E402
    assert_delta_is_minimal,
    assert_lookups_equal_manifest,
    check_consumer,
    load_model,
)


IMPORT_KEYS = (
    "mu_theta_operator",
    "c_s0",
    "rho_m",
    "rho_br",
    "W_0",
    "e_W",
    "v_bulk_normal_0",
    "q_out",
    "omega",
    "epsilon_shape",
    "Lambda_A_0",
    "Lambda_V_0",
    "Lambda_X_0",
    "tau_A",
    "tau_V",
    "tau_X",
    "face_normal",
    "conormal_deriv",
    "face_measure_shape_deriv",
    "face_velocity",
    "relative_flux",
    "kinematic_balance",
    "traction",
    "face_shift",
    "closure_shape_deriv",
    "W_bg",
    "w1_profile",
    "L_W",
    "sigma_W",
    "eta_bg",
    "rho_br_bg_rho4_constant",
    "z_impermeable",
    "z_by_regime",
    "z_by_parity",
    "added_mass",
    "grazing_behaviour",
    "face_response",
    "face_response_coeffs",
    "permeable_dissipative_by_regime_and_parity",
    "degenerate_loci_equations",
    "degenerate_loci_solution",
    "degenerate_loci_identically_satisfied",
    "degenerate_loci_inconsistent",
    "degenerate_loci_real_admissible",
)

CLASS_TAGS = {"KNOB", "STRUCTURAL", "COORDINATE", "CONTROL", "PREMISE", "DERIVED"}
PRIMARY_TASKS = ("DTN", "FACE_RESPONSE")
CONTROL_TASKS = (
    "REP_INVARIANCE",
    "INDEPENDENCE",
    "FORM",
    "UNIFORM",
    "ZERO_JET",
    "BRANCH",
    "HOMOGENEITY",
)
INTENDED_EXPORT_CANDIDATE_ROOTS = frozenset(
    {
        "dtn_flat_symbol",
        "dtn_operator",
        "dtn_kernel",
        "face_response",
        "face_response_coeffs",
    }
)
INTENDED_EXPORT_WRITE_ROOTS = frozenset(
    {
        "dtn_flat_symbol",
        "dtn_operator",
        "dtn_kernel",
        "s11c_c1_face_response",
        "s11c_c1_face_response_coeffs",
    }
)
ANCHORINGS = ("LAB_HELD", "MATERIAL_ADVECTED")
FACES = (1, -1)
DENSITIES = ("RHO4_CONSTANT", "RHOBR_CONSTANT")
DIRECTIONS = (1, 2, 3)
PARITIES = ("THICKNESS", "CENTRE_SHIFT")
REGIMES = ("PROPAGATING", "EVANESCENT", "GRAZING")

PROVED_EQUAL = "PROVED_EQUAL"
PROVED_DIFFERENT = "PROVED_DIFFERENT"
UNDECIDED = "UNDECIDED"

DIM_ZERO = sp.ImmutableMatrix([0, 0, 0])
DIM_L = sp.ImmutableMatrix([1, 0, 0])
DIM_T = sp.ImmutableMatrix([0, 1, 0])
DIM_M = sp.ImmutableMatrix([0, 0, 1])


def dim_add(*dimensions: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(sum((sp.Matrix(d) for d in dimensions), sp.zeros(3, 1)))


def dim_sub(left: sp.MatrixBase, right: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(sp.Matrix(left) - sp.Matrix(right))


DIM_VELOCITY = dim_sub(DIM_L, DIM_T)
DIM_RHO4 = dim_add(DIM_M, -4 * DIM_L)
DIM_RHOBR = dim_add(DIM_M, -3 * DIM_L)
DIM_PRESSURE = dim_add(DIM_M, -2 * DIM_L, -2 * DIM_T)
DIM_AFFINITY = dim_add(2 * DIM_L, -2 * DIM_T)
DIM_FLUX = dim_add(DIM_M, -3 * DIM_L, -DIM_T)
DIM_IMPEDANCE = dim_sub(DIM_PRESSURE, DIM_VELOCITY)
DIM_KERNEL = dim_add(DIM_IMPEDANCE, 3 * DIM_L)
DIM_MU_THETA = dim_add(DIM_M, -DIM_L, -2 * DIM_T)
DIM_LAMBDA_A = dim_sub(DIM_FLUX, DIM_AFFINITY)
DIM_LAMBDA_V = dim_sub(DIM_FLUX, DIM_VELOCITY)
DIM_LAMBDA_X = dim_sub(DIM_PRESSURE, DIM_AFFINITY)


NEW_SYMBOLS: dict[sp.Basic, dict[str, object]] = {}
SYMBOL_DIMENSIONS: dict[sp.Basic, sp.ImmutableMatrix] = {}


def new_symbol(
    name: str,
    class_tag: str,
    description: str,
    dimension: sp.MatrixBase,
    **assumptions: object,
) -> sp.Symbol:
    if class_tag not in CLASS_TAGS:
        raise ValueError(class_tag)
    value = sp.Symbol(name, **assumptions)
    NEW_SYMBOLS[value] = {"class": class_tag, "description": description}
    SYMBOL_DIMENSIONS[value] = sp.ImmutableMatrix(dimension)
    return value


# Fresh kernel legs and Fourier-profile representatives.  The imported q_out
# below remains the flat S11b branch object; these are the two live legs of the
# curved operator and are tied to it in the emitted branch package.
k_out = tuple(
    new_symbol(f"s11cc1_k_output_{i}", "COORDINATE", f"DtN output momentum leg {i}", -DIM_L, real=True)
    for i in DIRECTIONS
)
k_in = tuple(
    new_symbol(f"s11cc1_k_input_{i}", "COORDINATE", f"DtN input momentum leg {i}", -DIM_L, real=True)
    for i in DIRECTIONS
)
q_out_k = new_symbol(
    "s11cc1_q_out_output", "COORDINATE", "radiation-selected output-leg normal momentum", -DIM_L
)
q_out_kp = new_symbol(
    "s11cc1_q_out_input", "COORDINATE", "radiation-selected input-leg normal momentum", -DIM_L
)
w1_hat = new_symbol(
    "s11cc1_w1_profile_hat_transfer",
    "DERIVED",
    "Fourier transform of the dimensionless thickness profile at output minus input momentum",
    3 * DIM_L,
)
w1_jet_hat = tuple(
    new_symbol(
        f"s11cc1_w1_profile_jet_hat_{i}",
        "DERIVED",
        f"Fourier transform of the dimensionless face-tilt profile in direction {i}",
        3 * DIM_L,
    )
    for i in DIRECTIONS
)

delta_k = sp.Mul(*(sp.DiracDelta(a - b) for a, b in zip(k_out, k_in)))


def operator_symbol(name: str, description: str) -> sp.Symbol:
    return new_symbol(name, "DERIVED", description, DIM_ZERO, commutative=False)


identity_operator = new_symbol(
    "s11cc1_identity_operator", "STRUCTURAL", "identity on one face-field Hilbert space", DIM_ZERO,
    commutative=False,
)
flat_normal_operator = operator_symbol(
    "s11cc1_flat_normal_dtn_operator", "flat outgoing potential-to-conormal operator"
)
flat_normal_inverse = operator_symbol(
    "s11cc1_flat_normal_dtn_inverse", "inverse flat outgoing conormal operator"
)
height_multiplication_operator = operator_symbol(
    "s11cc1_height_multiplication_operator", "multiplication by outward background face height"
)
div_height_grad_operator = operator_symbol(
    "s11cc1_div_height_grad_operator", "divergence of height times in-plane gradient"
)
flat_impedance_operator = operator_symbol(
    "s11cc1_flat_impedance_operator", "flat outgoing pressure-to-normal-velocity impedance operator"
)
first_shape_impedance_operator = operator_symbol(
    "s11cc1_first_shape_impedance_operator", "first-shape-order impedance operator"
)

z_case_operators: dict[tuple[str, int], sp.Symbol] = {}
for anchoring in ANCHORINGS:
    for face in FACES:
        label = "plus" if face == 1 else "minus"
        z_case_operators[(anchoring, face)] = operator_symbol(
            f"s11cc1_dtn_operator_{anchoring.lower()}_{label}",
            f"curved DtN operator for {anchoring} face {label}",
        )

response_resolvents: dict[tuple[str, int, str], sp.Symbol] = {}
velocity_inputs: dict[tuple[str, int], sp.Symbol] = {}
mu_theta_inputs: dict[tuple[str, int], sp.Symbol] = {}
for anchoring in ANCHORINGS:
    for face in FACES:
        label = "plus" if face == 1 else "minus"
        velocity_inputs[(anchoring, face)] = new_symbol(
            f"s11cc1_V_{anchoring.lower()}_{label}",
            "COORDINATE",
            f"face-velocity input amplitude for {anchoring} face {label}",
            DIM_VELOCITY,
        )
        mu_theta_inputs[(anchoring, face)] = new_symbol(
            f"s11cc1_mu_theta_{anchoring.lower()}_{label}",
            "COORDINATE",
            f"opaque mu_theta input amplitude for {anchoring} face {label}",
            DIM_MU_THETA,
        )
        for density in DENSITIES:
            response_resolvents[(anchoring, face, density)] = operator_symbol(
                f"s11cc1_response_resolvent_{anchoring.lower()}_{label}_{density.lower()}",
                f"closed-face response resolvent for {anchoring}, face {label}, {density}",
            )


class UnevaluatedDisplayPrinter(StrPrinter):
    def _print_MatMul(self, expression: sp.MatMul) -> str:
        return "MatMul(" + ", ".join(self._print(argument) for argument in expression.args) + ")"

    def _print_Inverse(self, expression: sp.Inverse) -> str:
        return "Inverse(" + self._print(expression.arg) + ")"


DISPLAY_PRINTER = UnevaluatedDisplayPrinter()


def casify(value: object) -> object:
    if isinstance(value, bool):
        return sp.true if value else sp.false
    if isinstance(value, str):
        return Str(value)
    if isinstance(value, Mapping):
        return sp.Tuple(*(sp.Tuple(casify(k), casify(v)) for k, v in value.items()))
    if isinstance(value, sp.Tuple):
        return sp.Tuple(*(casify(item) for item in value))
    if isinstance(value, (tuple, list, set, frozenset)):
        return sp.Tuple(*(casify(item) for item in value))
    if isinstance(value, sp.MatrixBase) and not isinstance(value, sp.ImmutableMatrix):
        return sp.ImmutableMatrix(value)
    return value


def render(value: object) -> str:
    value = casify(value)
    if isinstance(value, (sp.Basic, sp.MatrixBase)):
        return sp.srepr(value)
    return repr(value)


def display(value: object) -> str:
    value = casify(value)
    if isinstance(value, (sp.Basic, sp.MatrixBase)):
        return DISPLAY_PRINTER.doprint(value)
    return repr(value)


def object_difference(left: object, right: object) -> object:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        keys = tuple(dict.fromkeys((*left.keys(), *right.keys())))
        return {key: object_difference(left.get(key, sp.Tuple()), right.get(key, sp.Tuple())) for key in keys}
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
            return sp.Add(left, -right, evaluate=False)
        except TypeError:
            return sp.Function("S11CC1ObjectDifference")(casify(left), casify(right))
    return sp.Function("S11CC1ObjectDifference")(casify(left), casify(right))


def simplified_object_difference(left: object, right: object) -> object:
    """Reduce a CAS comparison while preserving its container topology."""
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        keys = tuple(dict.fromkeys((*left.keys(), *right.keys())))
        return {
            key: simplified_object_difference(left.get(key, sp.Tuple()), right.get(key, sp.Tuple()))
            for key in keys
        }
    if isinstance(left, sp.MatrixBase) and isinstance(right, sp.MatrixBase) and left.shape == right.shape:
        return sp.ImmutableMatrix(left - right).applyfunc(sp.simplify)
    left_container = isinstance(left, (tuple, list, sp.Tuple))
    right_container = isinstance(right, (tuple, list, sp.Tuple))
    if left_container and right_container and len(left) == len(right):
        return sp.Tuple(*(simplified_object_difference(a, b) for a, b in zip(left, right)))
    if isinstance(left, sp.Basic) and isinstance(right, sp.Basic):
        try:
            return sp.simplify(left - right)
        except TypeError:
            return object_difference(left, right)
    return object_difference(left, right)


def total_compare(left: object, right: object) -> tuple[str, object]:
    """F9's three-valued comparison, total over imported row shapes."""
    if isinstance(left, bool) or isinstance(right, bool):
        if isinstance(left, bool) and isinstance(right, bool):
            return (PROVED_EQUAL if left is right else PROVED_DIFFERENT), sp.Tuple(casify(left), casify(right))
        return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
    if isinstance(left, str) or isinstance(right, str):
        if isinstance(left, str) and isinstance(right, str):
            return (PROVED_EQUAL if left == right else PROVED_DIFFERENT), sp.Tuple(Str(left), Str(right))
        return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        if not isinstance(left, Mapping) or not isinstance(right, Mapping) or set(left) != set(right):
            return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
        comparisons = [(key, total_compare(left[key], right[key])) for key in sorted(left, key=str)]
        statuses = [item[1][0] for item in comparisons]
        detail = sp.Tuple(
            *(sp.Tuple(Str(str(key)), Str(status), casify(operands)) for key, (status, operands) in comparisons)
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
        detail = sp.Tuple(*(sp.Tuple(Str(status), casify(operands)) for status, operands in comparisons))
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
        if all(status == PROVED_EQUAL for status in statuses):
            return PROVED_EQUAL, sp.ImmutableMatrix(left - right)
        if any(status == PROVED_DIFFERENT for status in statuses):
            return PROVED_DIFFERENT, sp.ImmutableMatrix(left - right)
        return UNDECIDED, sp.ImmutableMatrix(left - right)
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
        try:
            equivalent = sp.simplify_logic(sp.Equivalent(left, right))
        except Exception:
            return UNDECIDED, sp.Tuple(left, right)
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
        if sp.count_ops(left) + sp.count_ops(right) <= 40:
            try:
                residual = sp.simplify(left - right)
            except Exception:
                residual = None
            if residual == 0:
                return PROVED_EQUAL, residual
            if isinstance(residual, sp.Basic) and not residual.free_symbols:
                return PROVED_DIFFERENT, residual
        return UNDECIDED, sp.Tuple(left, right)
    try:
        equality = left == right
    except Exception:
        equality = None
    if equality is True:
        return PROVED_EQUAL, sp.Tuple(casify(left), casify(right))
    if equality is False:
        return PROVED_DIFFERENT, sp.Tuple(casify(left), casify(right))
    return UNDECIDED, sp.Tuple(casify(left), casify(right))


def package(value: object, grades: tuple[tuple[int, int, int], ...], dimension: object, **fields: object) -> sp.Tuple:
    rows = [
        sp.Tuple(Str("VALUE"), casify(value)),
        sp.Tuple(Str("MULTIGRADE"), casify(grades)),
        sp.Tuple(Str("DIMENSION_L_T_M"), casify(dimension)),
    ]
    rows.extend(sp.Tuple(Str(name), casify(field)) for name, field in fields.items())
    return sp.Tuple(*rows)


def case_map(cases: Mapping[object, object]) -> sp.Tuple:
    return sp.Tuple(*(sp.Tuple(casify(key), casify(value)) for key, value in cases.items()))


def tuple_record(value: sp.Tuple) -> dict[str, object]:
    return {str(key): item for key, item in value}


def atom_named(value: object, name: str) -> sp.Symbol:
    atoms = [atom for atom in casify(value).atoms(sp.Symbol) if atom.name == name]
    if len(atoms) != 1:
        raise RuntimeError(f"expected one {name!r} atom, found {len(atoms)}")
    return atoms[0]


@dataclass(frozen=True)
class Inputs:
    rows: Mapping[str, object]
    c_s0: sp.Symbol
    rho_m: sp.Symbol
    rho_br: sp.Symbol
    W0: sp.Symbol
    e_W: sp.Symbol
    q_out: sp.Symbol
    omega: sp.Symbol
    epsilon: sp.Symbol
    lambda_A0: sp.Symbol
    lambda_V0: sp.Symbol
    lambda_X0: sp.Symbol
    tau_A: sp.Symbol
    tau_V: sp.Symbol
    tau_X: sp.Symbol
    W_bg: sp.Symbol
    w1_profile: sp.Symbol
    L_W: sp.Symbol
    sigma_W: sp.Symbol
    eta: sp.Symbol
    rho_br_bg_rho4: sp.Symbol


@dataclass(frozen=True)
class AuditModel:
    inputs: Inputs
    objects: Mapping[str, object]


def ledger_value(ledger: Mapping[str, Mapping[str, object]], key: str) -> object:
    return ledger[key]["value"]


def bind_inputs(ledger: Mapping[str, Mapping[str, object]]) -> Inputs:
    # Every access here is a direct construction/control bind from §1–§5.
    rows = {
        "mu_theta_operator": ledger_value(ledger, "mu_theta_operator"),
        "face_normal": ledger_value(ledger, "face_normal"),
        "conormal_deriv": ledger_value(ledger, "conormal_deriv"),
        "face_measure_shape_deriv": ledger_value(ledger, "face_measure_shape_deriv"),
        "face_velocity": ledger_value(ledger, "face_velocity"),
        "relative_flux": ledger_value(ledger, "relative_flux"),
        "kinematic_balance": ledger_value(ledger, "kinematic_balance"),
        "traction": ledger_value(ledger, "traction"),
        "face_shift": ledger_value(ledger, "face_shift"),
        "closure_shape_deriv": ledger_value(ledger, "closure_shape_deriv"),
        "z_impermeable": ledger_value(ledger, "z_impermeable"),
        "z_by_regime": ledger_value(ledger, "z_by_regime"),
        "z_by_parity": ledger_value(ledger, "z_by_parity"),
        "added_mass": ledger_value(ledger, "added_mass"),
        "grazing_behaviour": ledger_value(ledger, "grazing_behaviour"),
        "face_response": ledger_value(ledger, "face_response"),
        "face_response_coeffs": ledger_value(ledger, "face_response_coeffs"),
        "permeable_dissipative_by_regime_and_parity": ledger_value(
            ledger, "permeable_dissipative_by_regime_and_parity"
        ),
        "degenerate_loci_equations": ledger_value(ledger, "degenerate_loci_equations"),
        "degenerate_loci_solution": ledger_value(ledger, "degenerate_loci_solution"),
        "degenerate_loci_identically_satisfied": ledger_value(
            ledger, "degenerate_loci_identically_satisfied"
        ),
        "degenerate_loci_inconsistent": ledger_value(ledger, "degenerate_loci_inconsistent"),
        "degenerate_loci_real_admissible": ledger_value(ledger, "degenerate_loci_real_admissible"),
        "v_bulk_normal_0": ledger_value(ledger, "v_bulk_normal_0"),
    }
    return Inputs(
        MappingProxyType(rows),
        ledger_value(ledger, "c_s0"),
        ledger_value(ledger, "rho_m"),
        ledger_value(ledger, "rho_br"),
        ledger_value(ledger, "W_0"),
        ledger_value(ledger, "e_W"),
        ledger_value(ledger, "q_out"),
        ledger_value(ledger, "omega"),
        ledger_value(ledger, "epsilon_shape"),
        ledger_value(ledger, "Lambda_A_0"),
        ledger_value(ledger, "Lambda_V_0"),
        ledger_value(ledger, "Lambda_X_0"),
        ledger_value(ledger, "tau_A"),
        ledger_value(ledger, "tau_V"),
        ledger_value(ledger, "tau_X"),
        ledger_value(ledger, "W_bg"),
        ledger_value(ledger, "w1_profile"),
        ledger_value(ledger, "L_W"),
        ledger_value(ledger, "sigma_W"),
        ledger_value(ledger, "eta_bg"),
        ledger_value(ledger, "rho_br_bg_rho4_constant"),
    )


def lambda_kernels(inputs: Inputs, overrides: Mapping[str, sp.Expr] | None = None) -> tuple[sp.Expr, ...]:
    overrides = {} if overrides is None else overrides
    values = {
        "A": inputs.lambda_A0 / (1 - sp.I * inputs.omega * inputs.tau_A),
        "V": inputs.lambda_V0 / (1 - sp.I * inputs.omega * inputs.tau_V),
        "X": inputs.lambda_X0 / (1 - sp.I * inputs.omega * inputs.tau_X),
    }
    values.update(overrides)
    return values["A"], values["V"], values["X"]


def density_value(inputs: Inputs, density: str) -> sp.Symbol:
    return inputs.rho_br_bg_rho4 if density == "RHO4_CONSTANT" else inputs.rho_br


def shape_source(inputs: Inputs, *, flip_upper_x1: bool = False, face: int = 1, ablate: int | None = None) -> dict[str, object]:
    height = inputs.eta * inputs.W0 * w1_hat / 2
    tilt = []
    for direction, jet in zip(DIRECTIONS, w1_jet_hat):
        factor = sp.Integer(0) if ablate == direction else sp.Integer(1)
        if flip_upper_x1 and face == 1 and direction == 1:
            factor = -factor
        tilt.append(factor * inputs.sigma_W * jet / 2)
    return {"height_hat": height, "tilt_hat": tuple(tilt)}


def dtn_flat_symbol(inputs: Inputs, q_leg: sp.Expr) -> sp.Expr:
    return inputs.rho_m * inputs.omega / q_leg


def dtn_first_kernel(
    inputs: Inputs,
    source: Mapping[str, object],
    q_output: sp.Expr,
    q_input: sp.Expr,
    *,
    route: str = "EULERIAN",
) -> sp.Expr:
    height = source["height_hat"]
    tilt = source["tilt_hat"]
    k_input_sq = sp.Add(*(component**2 for component in k_in))
    tilt_dot_input = sp.Add(*(component * momentum for component, momentum in zip(tilt, k_in)))
    kappa_sq = inputs.omega**2 / inputs.c_s0**2
    remainder = height * (k_input_sq - kappa_sq) - sp.I * tilt_dot_input
    if route == "EULERIAN":
        return (
            sp.I
            * inputs.rho_m
            * inputs.omega
            * (height * q_output * q_input + remainder)
            / (q_output * q_input)
        )
    if route == "HANZAWA":
        return (
            sp.I
            * inputs.rho_m
            * inputs.omega
            * (height + remainder / (q_output * q_input))
        )
    raise ValueError(route)


def hanzawa_first_kernel(
    inputs: Inputs,
    source: Mapping[str, object],
    q_output: sp.Expr,
    q_input: sp.Expr,
) -> tuple[sp.Expr, sp.Tuple]:
    """Construct the second route with the outgoing Neumann layer potential."""
    normal_coordinate = sp.Symbol("s11cc1_layer_normal_coordinate", nonnegative=True, real=True)
    phi_flat_amplitude = sp.Symbol("s11cc1_layer_phi_flat_amplitude")
    height = source["height_hat"]
    tilt = source["tilt_hat"]
    k_input_sq = sp.Add(*(component**2 for component in k_in))
    tilt_dot_input = sp.Add(*(component * momentum for component, momentum in zip(tilt, k_in)))
    kappa_sq = inputs.omega**2 / inputs.c_s0**2

    conormal_source = height * (k_input_sq - kappa_sq) - sp.I * tilt_dot_input
    flat_boundary_equation = sp.Eq(
        sp.I * q_input * phi_flat_amplitude, 1, evaluate=False
    )
    flat_amplitudes = sp.solve(
        (flat_boundary_equation,),
        (phi_flat_amplitude,),
        dict=True,
    )
    if len(flat_amplitudes) != 1:
        raise RuntimeError(
            f"layer-potential flat boundary solve returned {len(flat_amplitudes)} branches"
        )
    phi_flat_trace = flat_amplitudes[0][phi_flat_amplitude]
    layer_density = -conormal_source * phi_flat_trace
    outgoing_neumann_poisson_kernel = sp.exp(
        sp.I * q_output * normal_coordinate
    ) / (sp.I * q_output)
    outgoing_scattered_field = layer_density * outgoing_neumann_poisson_kernel
    phi_scattered_trace = outgoing_scattered_field.subs(normal_coordinate, 0)
    shifted_pressure_trace = sp.I * inputs.rho_m * inputs.omega * (
        sp.I * q_input * height * phi_flat_trace + phi_scattered_trace
    )
    kernel = sp.factor(shifted_pressure_trace)
    construction = sp.Tuple(
        sp.Tuple(
            Str("SECOND_CONSTRUCTION"),
            Str("OUTGOING_NEUMANN_LAYER_POTENTIAL"),
        ),
        sp.Tuple(
            Str("PULLED_BACK_BULK_DISPERSION"),
            sp.Eq(q_input**2 + k_input_sq, kappa_sq, evaluate=False),
            sp.Eq(
                q_output**2 + sp.Add(*(component**2 for component in k_out)),
                kappa_sq,
                evaluate=False,
            ),
        ),
        sp.Tuple(Str("FLAT_BOUNDARY_EQUATION"), flat_boundary_equation),
        sp.Tuple(Str("CURVED_CONORMAL_LAYER_DENSITY"), layer_density),
        sp.Tuple(
            Str("OUTGOING_NEUMANN_POISSON_KERNEL"),
            outgoing_neumann_poisson_kernel,
        ),
        sp.Tuple(
            Str("OUTGOING_SCATTERED_FIELD"),
            outgoing_scattered_field,
        ),
        sp.Tuple(
            Str("LAYER_POTENTIAL_RADIATION_CONDITION"),
            Str("OUTGOING_OR_DECAYING_ALREADY_SELECTED_BRANCH"),
            q_output,
        ),
    )
    return kernel, construction


def reduce_on_shell(
    expression: sp.Expr,
    inputs: Inputs,
    q_leg: sp.Symbol,
    momentum_leg: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Apply the bulk dispersion through the normal-momentum square."""
    tangential_sq = sp.Add(*(component**2 for component in momentum_leg))
    dispersion = inputs.omega**2 / inputs.c_s0**2 - tangential_sq
    collected = sp.together(expression)
    return sp.factor(sp.cancel(collected.subs(q_leg**2, dispersion)))


def closed_coefficients(
    inputs: Inputs,
    z0_output: sp.Expr,
    z0_input: sp.Expr,
    z1: sp.Expr,
    rho_background: sp.Expr,
    lambdas: tuple[sp.Expr, sp.Expr, sp.Expr] | None = None,
) -> dict[str, object]:
    lambda_A, lambda_V, lambda_X = lambda_kernels(inputs) if lambdas is None else lambdas
    response_output = 1 / (1 + lambda_A * z0_output / inputs.rho_m**2)
    response_input = 1 / (1 + lambda_A * z0_input / inputs.rho_m**2)
    pressure_V0 = response_output * z0_output * (1 + lambda_V / inputs.rho_m)
    pressure_mu0 = response_output * z0_output * lambda_A / (inputs.rho_m * rho_background)
    pressure_V1 = response_output * z1 * response_input * (1 + lambda_V / inputs.rho_m)
    pressure_mu1 = response_output * z1 * response_input * lambda_A / (inputs.rho_m * rho_background)
    flux_V0 = lambda_V - lambda_A * pressure_V0 / inputs.rho_m
    flux_mu0 = lambda_A / rho_background - lambda_A * pressure_mu0 / inputs.rho_m
    flux_V1 = -lambda_A * pressure_V1 / inputs.rho_m
    flux_mu1 = -lambda_A * pressure_mu1 / inputs.rho_m
    traction_V0 = (1 - lambda_X / inputs.rho_m) * pressure_V0
    traction_mu0 = (1 - lambda_X / inputs.rho_m) * pressure_mu0 + lambda_X / rho_background
    traction_V1 = (1 - lambda_X / inputs.rho_m) * pressure_V1
    traction_mu1 = (1 - lambda_X / inputs.rho_m) * pressure_mu1
    return {
        "PRESSURE": ((pressure_V0, pressure_mu0), (pressure_V1, pressure_mu1)),
        "RELATIVE_FLUX": ((flux_V0, flux_mu0), (flux_V1, flux_mu1)),
        "TRACTION_NORMAL_AMPLITUDE": (
            (traction_V0, traction_mu0),
            (traction_V1, traction_mu1),
        ),
        "RESOLVENT_SYMBOL_OUTPUT": response_output,
        "RESOLVENT_SYMBOL_INPUT": response_input,
    }


def first_shape_true_area_weight(source: Mapping[str, object]) -> sp.Expr:
    area_order = sp.Symbol("s11cc1_true_area_expansion_order", real=True)
    slope_sq = sp.Add(*(component**2 for component in source["tilt_hat"]))
    area = sp.sqrt(1 + area_order**2 * slope_sq)
    return sp.series(area, area_order, 0, 2).removeO().subs(area_order, 1)


def outgoing_farfield_poynting(
    inputs: Inputs,
    source: Mapping[str, object],
    q_output: sp.Expr,
    q_input: sp.Expr,
    velocity: sp.Expr,
) -> tuple[sp.Expr, sp.Tuple]:
    """Solve the outgoing bulk amplitudes and evaluate their remote flux."""
    control_radius = sp.Symbol("s11cc1_farfield_control_radius", positive=True, real=True)
    phi_flat_amplitude, phi_scattered_amplitude = sp.symbols(
        "s11cc1_farfield_phi_flat_amplitude s11cc1_farfield_phi_scattered_amplitude"
    )
    height = source["height_hat"]
    tilt = source["tilt_hat"]
    k_input_sq = sp.Add(*(component**2 for component in k_in))
    kappa_sq = inputs.omega**2 / inputs.c_s0**2
    tilt_dot_input = sp.Add(*(component * momentum for component, momentum in zip(tilt, k_in)))
    conormal_source = height * (k_input_sq - kappa_sq) - sp.I * tilt_dot_input
    boundary_system = (
        sp.Eq(sp.I * q_input * phi_flat_amplitude, velocity, evaluate=False),
        sp.Eq(
            sp.I * q_output * phi_scattered_amplitude
            + conormal_source * phi_flat_amplitude,
            0,
            evaluate=False,
        ),
    )
    amplitudes = sp.solve(
        boundary_system,
        (phi_flat_amplitude, phi_scattered_amplitude),
        dict=True,
    )
    if len(amplitudes) != 1:
        raise RuntimeError(f"far-field amplitude solve returned {len(amplitudes)} branches")
    phi_flat = amplitudes[0][phi_flat_amplitude] * sp.exp(
        sp.I * q_input * control_radius
    )
    phi_scattered = amplitudes[0][phi_scattered_amplitude] * sp.exp(
        sp.I * q_output * control_radius
    )
    pressure_flat = sp.I * inputs.rho_m * inputs.omega * phi_flat
    pressure_scattered = sp.I * inputs.rho_m * inputs.omega * phi_scattered
    normal_velocity_flat = sp.diff(phi_flat, control_radius)
    normal_velocity_scattered = sp.diff(phi_scattered, control_radius)
    poynting_through_first_shape = sp.re(
        pressure_flat * sp.conjugate(normal_velocity_flat)
        + pressure_flat * sp.conjugate(normal_velocity_scattered)
        + pressure_scattered * sp.conjugate(normal_velocity_flat),
        evaluate=False,
    ) / 2
    control_surface_flux = sp.simplify(poynting_through_first_shape)
    farfield_flux = sp.limit(control_surface_flux, control_radius, sp.oo)
    construction = sp.Tuple(
        sp.Tuple(Str("BULK_OPERATOR"), Str("REST_FRAME_HELMHOLTZ"), kappa_sq),
        sp.Tuple(Str("DRIVEN_BOUNDARY_SYSTEM"), casify(boundary_system)),
        sp.Tuple(Str("OUTGOING_PHI_FLAT"), phi_flat),
        sp.Tuple(Str("OUTGOING_PHI_SCATTERED"), phi_scattered),
        sp.Tuple(
            Str("RADIATION_CONDITION"),
            Str("OUTGOING_ALREADY_SELECTED_BRANCH"),
            q_output,
            q_input,
        ),
        sp.Tuple(Str("CONTROL_SURFACE_FLUX"), control_surface_flux),
        sp.Tuple(Str("FARFIELD_LIMIT"), farfield_flux),
    )
    return farfield_flux, construction


def response_operator_case(inputs: Inputs, anchoring: str, face: int, density: str) -> sp.Tuple:
    lambda_A, lambda_V, lambda_X = lambda_kernels(inputs)
    rho_background = density_value(inputs, density)
    z_operator = z_case_operators[(anchoring, face)]
    resolvent = response_resolvents[(anchoring, face, density)]
    velocity = inputs.epsilon * velocity_inputs[(anchoring, face)]
    mu_theta = inputs.epsilon * mu_theta_inputs[(anchoring, face)]
    mu_s = mu_theta / rho_background
    source = (1 + lambda_V / inputs.rho_m) * velocity + lambda_A * mu_s / inputs.rho_m
    pressure = resolvent * z_operator * source
    flux = lambda_A * (mu_s - pressure / inputs.rho_m) + lambda_V * velocity
    traction_amplitude = pressure + lambda_X * (mu_s - pressure / inputs.rho_m)
    return sp.Tuple(
        sp.Tuple(Str("RESOLVENT"), resolvent),
        sp.Tuple(
            Str("RESOLVENT_DEFINITION"),
            sp.Tuple(resolvent, identity_operator + lambda_A * z_operator / inputs.rho_m**2, Str("OPERATOR_INVERSE")),
        ),
        sp.Tuple(Str("MU_S"), mu_s),
        sp.Tuple(Str("DELTA_P"), pressure),
        sp.Tuple(Str("J"), flux),
        sp.Tuple(
            Str("T"),
            sp.Tuple(-traction_amplitude, Str("TIMES_CURVED_FACE_NORMAL"), Str("face_normal")),
        ),
    )


def regime_replacements(inputs: Inputs) -> dict[str, tuple[sp.Expr, sp.Expr]]:
    q_prop_out = sp.Symbol("s11cc1_q_prop_output", positive=True, real=True) * sp.sign(inputs.omega)
    q_prop_in = sp.Symbol("s11cc1_q_prop_input", positive=True, real=True) * sp.sign(inputs.omega)
    kappa_out = sp.Symbol("s11cc1_kappa_evan_output", positive=True, real=True)
    kappa_in = sp.Symbol("s11cc1_kappa_evan_input", positive=True, real=True)
    return {
        "PROPAGATING": (q_prop_out, q_prop_in),
        "EVANESCENT": (sp.I * kappa_out, sp.I * kappa_in),
        "GRAZING": (sp.Integer(0), sp.Integer(0)),
    }


def substitute_regime(
    expression: sp.Expr,
    output_regime: str,
    input_regime: str,
    inputs: Inputs,
    *,
    evaluate_grazing: bool = True,
) -> sp.Expr:
    replacements = regime_replacements(inputs)
    output_value = replacements[output_regime][0]
    input_value = replacements[input_regime][1]
    if not evaluate_grazing and output_regime == "GRAZING" and input_regime == "GRAZING":
        return sp.Limit(sp.Limit(expression, q_out_k, 0, dir="+"), q_out_kp, 0, dir="+")
    if not evaluate_grazing and output_regime == "GRAZING":
        return sp.Limit(expression.xreplace({q_out_kp: input_value}), q_out_k, 0, dir="+")
    if not evaluate_grazing and input_regime == "GRAZING":
        return sp.Limit(expression.xreplace({q_out_k: output_value}), q_out_kp, 0, dir="+")
    if output_regime == "GRAZING" and input_regime == "GRAZING":
        return sp.Tuple(
            sp.limit(sp.limit(expression, q_out_k, 0, dir="+"), q_out_kp, 0, dir="+"),
            sp.limit(sp.limit(expression, q_out_kp, 0, dir="+"), q_out_k, 0, dir="+"),
        )
    if output_regime == "GRAZING":
        return sp.limit(expression.xreplace({q_out_kp: input_value}), q_out_k, 0, dir="+")
    if input_regime == "GRAZING":
        return sp.limit(expression.xreplace({q_out_k: output_value}), q_out_kp, 0, dir="+")
    return expression.xreplace({q_out_k: output_value, q_out_kp: input_value})


def hermitian_kernel(expression: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    swap = {q_out_k: q_out_kp, q_out_kp: q_out_k}
    swap.update({a: b for a, b in zip(k_out, k_in)})
    swap.update({b: a for a, b in zip(k_out, k_in)})
    adjoint = sp.conjugate(expression.xreplace(swap))
    return adjoint, (expression + adjoint) / 2, (expression - adjoint) / (2 * sp.I)


def status_of(test_object: object) -> str:
    try:
        simplified = sp.simplify(test_object)
    except Exception:
        simplified = test_object
    if simplified is sp.true or simplified is True:
        return PROVED_EQUAL
    if simplified is sp.false or simplified is False:
        return PROVED_DIFFERENT
    return UNDECIDED


def locus_package(inputs: Inputs) -> dict[str, object]:
    lambda_A, _, _ = lambda_kernels(inputs)
    denominator = 1 + lambda_A * dtn_flat_symbol(inputs, q_out_k) / inputs.rho_m**2
    equation = sp.Eq(denominator, 0, evaluate=False)
    solutions = sp.solve((equation,), (inputs.lambda_A0,), dict=True)
    coefficient = sp.diff(denominator, inputs.lambda_A0)
    constant = denominator.xreplace({inputs.lambda_A0: sp.Integer(0)})
    identical_test = sp.And(
        sp.Eq(coefficient, 0, evaluate=False), sp.Eq(constant, 0, evaluate=False), evaluate=False
    )
    inconsistent_test = sp.And(
        sp.Eq(coefficient, 0, evaluate=False), sp.Ne(constant, 0, evaluate=False), evaluate=False
    )
    identical_status = status_of(identical_test)
    inconsistent_status = status_of(inconsistent_test)
    real_rows = []
    for branch in solutions:
        branch_value = branch[inputs.lambda_A0]
        test = sp.Eq(sp.im(branch_value), 0, evaluate=False)
        branch_status = status_of(test)
        token = "ADMISSIBLE" if branch_status == PROVED_EQUAL else (
            "EXCLUDED" if branch_status == PROVED_DIFFERENT else "UNDECIDED"
        )
        real_rows.append(
            sp.Tuple(
                casify(branch), Str(token), test, sp.Tuple(branch_value, inputs.tau_A, inputs.omega, q_out_k)
            )
        )
    return {
        "EQUATIONS": sp.Tuple(equation),
        "SOLUTION": sp.Tuple(
            inputs.lambda_A0,
            sp.Tuple(*(sp.Tuple(*(sp.Eq(key, value, evaluate=False) for key, value in branch.items())) for branch in solutions)),
        ),
        "IDENTICALLY_SATISFIED": sp.Tuple(
            Str(identical_status), identical_test, sp.Tuple(coefficient, constant)
        ),
        "INCONSISTENT": sp.Tuple(
            Str(inconsistent_status), inconsistent_test, sp.Tuple(coefficient, constant, inputs.lambda_A0)
        ),
        "REAL_ADMISSIBLE": sp.Tuple(*real_rows),
        "DIAGONAL_DENOMINATOR": denominator,
    }


def port_matrix(
    inputs: Inputs,
    z0_output: sp.Expr,
    z0_input: sp.Expr,
    bare_dtn: sp.Expr,
    lambdas: tuple[sp.Expr, sp.Expr, sp.Expr] | None = None,
) -> sp.ImmutableMatrix:
    lambda_A, lambda_V, lambda_X = lambda_kernels(inputs) if lambdas is None else lambdas
    response_output = 1 / (1 + lambda_A * z0_output / inputs.rho_m**2)
    response_input = 1 / (1 + lambda_A * z0_input / inputs.rho_m**2)
    z0_kernel = z0_output * delta_k
    z1_kernel = sp.expand(bare_dtn - z0_kernel)
    closed_dtn = response_output * z0_kernel + response_output * z1_kernel * response_input
    p_v = closed_dtn * (1 + lambda_V / inputs.rho_m)
    p_mu = closed_dtn * lambda_A / inputs.rho_m
    j_v = lambda_V * delta_k - lambda_A * p_v / inputs.rho_m
    j_mu = lambda_A * delta_k - lambda_A * p_mu / inputs.rho_m
    return sp.ImmutableMatrix(
        [
            [
                (1 - lambda_X / inputs.rho_m) * p_v,
                (1 - lambda_X / inputs.rho_m) * p_mu + lambda_X * delta_k,
            ],
            [j_v, j_mu],
        ]
    )


def port_hermitian(
    inputs: Inputs,
    bare_dtn: sp.Expr,
    lambdas: tuple[sp.Expr, ...] | None = None,
) -> sp.ImmutableMatrix:
    z0_output = dtn_flat_symbol(inputs, q_out_k)
    z0_input = dtn_flat_symbol(inputs, q_out_kp)
    matrix = port_matrix(inputs, z0_output, z0_input, bare_dtn, lambdas)
    swap = {q_out_k: q_out_kp, q_out_kp: q_out_k}
    swap.update({a: b for a, b in zip(k_out, k_in)})
    swap.update({b: a for a, b in zip(k_out, k_in)})
    reverse = sp.ImmutableMatrix(matrix.xreplace(swap))
    return sp.ImmutableMatrix(
        2,
        2,
        lambda row, column: sp.Mul(
            sp.Rational(1, 2),
            sp.Add(
                matrix[row, column],
                sp.conjugate(reverse[column, row], evaluate=False),
                evaluate=False,
            ),
            evaluate=False,
        ),
    )


def rederive_flat_face_response(inputs: Inputs) -> tuple[object, object]:
    target = inputs.rows["face_response"]
    p = atom_named(target, "delta_p_face")
    flux = atom_named(target, "J_face")
    velocity = atom_named(target, "V_face")
    mu_drive = atom_named(target, "mu_theta_drive")
    lambda_A, lambda_V, _ = lambda_kernels(inputs)
    zflat = dtn_flat_symbol(inputs, inputs.q_out)
    equations = sp.ImmutableMatrix(
        [
            p - zflat * (velocity + flux / inputs.rho_m),
            flux - lambda_A * (mu_drive / inputs.rho_br - p / inputs.rho_m) - lambda_V * velocity,
        ]
    )
    solved = sp.solve(tuple(equations), (p, flux), dict=True)
    if len(solved) != 1:
        raise RuntimeError(f"flat response solve returned {len(solved)} branches")
    pressure = sp.factor(solved[0][p])
    response = sp.Tuple(equations, pressure)
    coeffs = sp.Tuple(sp.factor(sp.diff(pressure, velocity)), sp.factor(sp.diff(pressure, mu_drive)))
    return response, coeffs


def build_model(ledger: Mapping[str, Mapping[str, object]]) -> AuditModel:
    inputs = bind_inputs(ledger)
    objects: dict[str, object] = {}
    kappa_sq = inputs.omega**2 / inputs.c_s0**2
    z0_out = dtn_flat_symbol(inputs, q_out_k)
    z0_in = dtn_flat_symbol(inputs, q_out_kp)

    flat_cases = {}
    kernel_cases = {}
    operator_cases = {}
    direct_kernels: dict[tuple[str, int], sp.Expr] = {}
    hanzawa_kernels: dict[tuple[str, int], sp.Expr] = {}
    hanzawa_constructions: dict[tuple[str, int], sp.Tuple] = {}
    g1_operator = (
        -flat_normal_operator * height_multiplication_operator * flat_normal_operator
        - div_height_grad_operator
        - kappa_sq * height_multiplication_operator
    )
    z1_composition = sp.I * inputs.rho_m * inputs.omega * (
        -flat_normal_inverse * g1_operator * flat_normal_inverse
    )
    z0_composition = sp.I * inputs.rho_m * inputs.omega * flat_normal_inverse
    for anchoring in ANCHORINGS:
        for face in FACES:
            key = (anchoring, face)
            source = shape_source(inputs, face=face)
            direct = dtn_first_kernel(inputs, source, q_out_k, q_out_kp, route="EULERIAN")
            hanzawa, hanzawa_construction = hanzawa_first_kernel(
                inputs, source, q_out_k, q_out_kp
            )
            direct_kernels[key] = direct
            hanzawa_kernels[key] = hanzawa
            hanzawa_constructions[key] = hanzawa_construction
            flat_cases[key] = package(
                z0_out,
                ((0, 0, 0),),
                DIM_IMPEDANCE,
                RATIO=sp.Tuple(Str("BULK_PRESSURE"), Str("BULK_OUTWARD_NORMAL_VELOCITY")),
                BRANCH_BINDING=sp.Tuple(inputs.q_out, q_out_k, Str("OUTPUT_LEG")),
            )
            kernel_cases[key] = package(
                sp.Tuple(
                    sp.Tuple(Str("FLAT_DIAGONAL"), z0_out * delta_k),
                    sp.Tuple(Str("FIRST_SHAPE"), direct),
                    sp.Tuple(Str("ASSEMBLED"), z0_out * delta_k + direct),
                ),
                ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
                sp.Tuple(DIM_KERNEL, DIM_KERNEL, DIM_KERNEL),
                LEG_CONVENTION=sp.Tuple(Str("OUTPUT"), casify(k_out), Str("INPUT"), casify(k_in)),
                PROFILE_TRANSFER=sp.Tuple(w1_hat, *w1_jet_hat),
                BRANCH_BINDING=sp.Tuple(inputs.q_out, q_out_k, q_out_kp),
            )
            z_operator = z_case_operators[key]
            operator_cases[key] = package(
                sp.Tuple(
                    sp.Tuple(flat_impedance_operator, z0_composition),
                    sp.Tuple(first_shape_impedance_operator, z1_composition),
                    sp.Eq(
                        z_operator,
                        flat_impedance_operator + first_shape_impedance_operator,
                        evaluate=False,
                    ),
                ),
                ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
                DIM_IMPEDANCE,
                HEIGHT_OPERATOR_BINDING=sp.Tuple(
                    height_multiplication_operator, (inputs.W_bg - inputs.W0) / 2, inputs.eta, inputs.sigma_W
                ),
                TWO_SIDED_COMPOSITION=g1_operator,
                RADIATION_CONDITION=sp.Tuple(Str("OUTGOING_OR_DECAYING"), inputs.q_out),
            )

    objects["DTN_FLAT_SYMBOL"] = case_map(flat_cases)
    objects["DTN_KERNEL"] = case_map(kernel_cases)
    objects["DTN_OPERATOR"] = case_map(operator_cases)

    rigid_height = sp.Symbol("s11cc1_rigid_translation_height")
    rigid_source = {"height_hat": rigid_height, "tilt_hat": (sp.Integer(0),) * 3}
    rigid_pre = dtn_first_kernel(inputs, rigid_source, q_out_k, q_out_kp)
    rigid_diagonal = rigid_pre.xreplace(
        {q_out_kp: q_out_k, **{b: a for a, b in zip(k_out, k_in)}}
    )
    rigid_after_dispersion = reduce_on_shell(rigid_diagonal, inputs, q_out_k, k_out)
    objects["DTN_RIGID_SHIFT_OPERAND"] = case_map(
        {key: package(sp.Tuple(rigid_diagonal, sp.Integer(0)), ((0, 1, 0),), DIM_KERNEL) for key in direct_kernels}
    )
    objects["DTN_RIGID_SHIFT_RESIDUAL"] = case_map(
        {key: package(rigid_after_dispersion, ((0, 1, 0),), DIM_KERNEL) for key in direct_kernels}
    )

    regime_cases = {}
    grazing_cases = {}
    for key, kernel in direct_kernels.items():
        for output_regime in REGIMES:
            for input_regime in REGIMES:
                regime_cases[(*key, output_regime, input_regime)] = package(
                    substitute_regime(kernel, output_regime, input_regime, inputs),
                    ((0, 1, 0), (0, 0, 1)),
                    DIM_KERNEL,
                    CONDITIONS=sp.Tuple(Str(output_regime), Str(input_regime)),
                )
        grazing_cases[key] = package(
            sp.Tuple(
                sp.Tuple(Str("OUTPUT_LEG"), sp.limit(kernel, q_out_k, 0, dir="+")),
                sp.Tuple(Str("INPUT_LEG"), sp.limit(kernel, q_out_kp, 0, dir="+")),
                sp.Tuple(
                    Str("BOTH_OUTPUT_THEN_INPUT"),
                    sp.limit(sp.limit(kernel, q_out_k, 0, dir="+"), q_out_kp, 0, dir="+"),
                ),
                sp.Tuple(
                    Str("BOTH_INPUT_THEN_OUTPUT"),
                    sp.limit(sp.limit(kernel, q_out_kp, 0, dir="+"), q_out_k, 0, dir="+"),
                ),
            ),
            ((0, 1, 0), (0, 0, 1)),
            DIM_KERNEL,
        )
    objects["DTN_BY_REGIME_PAIR"] = case_map(regime_cases)
    objects["DTN_GRAZING_BEHAVIOUR"] = case_map(grazing_cases)

    parity_cases = {}
    hermitian_cases = {}
    reactive_cases = {}
    for anchoring in ANCHORINGS:
        plus = direct_kernels[(anchoring, 1)]
        minus = direct_kernels[(anchoring, -1)]
        face_matrix = sp.diag(plus, minus)
        parity_transform = sp.ImmutableMatrix([[1, 1], [1, -1]]) / sp.sqrt(2)
        parity_matrix = sp.ImmutableMatrix(parity_transform.T * face_matrix * parity_transform)
        parity_cases[anchoring] = package(
            sp.Tuple(
                sp.Tuple(Str("FACE_BASIS"), face_matrix),
                sp.Tuple(Str("THICKNESS_CENTRE_BASIS"), parity_matrix),
                sp.Tuple(Str("OFF_DIAGONAL_BLOCKS"), parity_matrix[0, 1], parity_matrix[1, 0]),
            ),
            ((0, 1, 0), (0, 0, 1)),
            DIM_KERNEL,
        )
        for face in FACES:
            key = (anchoring, face)
            bare_dtn = z0_out * delta_k + direct_kernels[key]
            adjoint, hermitian, reactive = hermitian_kernel(bare_dtn)
            hermitian_cases[key] = package(
                sp.Tuple(
                    sp.Tuple(Str("FULL_BARE_DTN"), bare_dtn),
                    sp.Tuple(Str("TRUE_AREA_ADJOINT"), adjoint),
                    sp.Tuple(Str("HERMITIAN_PART"), hermitian),
                    sp.Tuple(Str("PAIRING_MEASURE_SOURCE"), Str("face_measure_shape_deriv")),
                ),
                ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
                DIM_KERNEL,
            )
            reactive_cases[key] = package(
                reactive,
                ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
                DIM_KERNEL,
            )
    objects["DTN_BY_PARITY"] = case_map(parity_cases)
    objects["DTN_HERMITIAN_PART"] = case_map(hermitian_cases)
    objects["DTN_REACTIVE_PART"] = case_map(reactive_cases)

    kappa_evan = sp.Symbol("s11cc1_kappa_evanescent", positive=True, real=True)
    m_add = sp.factor(inputs.rho_m / kappa_evan)
    objects["DTN_INERTIAL_LOADING"] = case_map(
        {
            (anchoring, face, "EVANESCENT_PURELY_REACTIVE"): package(
                sp.Tuple(
                    sp.Tuple(Str("PRESSURE"), -inputs.rho_m * inputs.omega**2 / kappa_evan),
                    sp.Tuple(Str("OUTWARD_ACCELERATION"), -inputs.omega**2),
                    sp.Tuple(Str("M_ADD"), m_add),
                ),
                ((0, 0, 0),),
                sp.Tuple(DIM_PRESSURE, dim_sub(DIM_VELOCITY, DIM_T), dim_sub(DIM_IMPEDANCE, -DIM_T)),
            )
            for anchoring in ANCHORINGS
            for face in FACES
        }
    )
    objects["DTN_TERM_ORIGINS"] = case_map(
        {
            "FLAT_OUTGOING_SOLUTION": package(
                sp.Tuple(inputs.rows["conormal_deriv"], inputs.q_out), ((0, 0, 0),), DIM_IMPEDANCE
            ),
            "SHIFTED_TRACE": package(inputs.rows["face_shift"], ((1, 1, 0),), DIM_ZERO),
            "CURVED_NORMAL": package(inputs.rows["face_normal"], ((0, 0, 1),), DIM_ZERO),
            "TRUE_AREA": package(inputs.rows["face_measure_shape_deriv"], ((0, 0, 0),), DIM_ZERO),
            "TWO_SIDED_COMPOSITION": package(g1_operator, ((0, 1, 0), (0, 0, 1)), -DIM_L),
            "PROFILE_BINDING": package(
                sp.Tuple(inputs.W_bg, inputs.w1_profile, inputs.L_W, inputs.eta, inputs.sigma_W),
                ((0, 1, 0), (0, 0, 1)),
                sp.Tuple(DIM_L, DIM_ZERO, DIM_L, DIM_ZERO, DIM_ZERO),
            ),
        }
    )

    response_cases = {}
    response_coeff_cases = {}
    for anchoring in ANCHORINGS:
        for face in FACES:
            kernel = direct_kernels[(anchoring, face)]
            for density in DENSITIES:
                key = (anchoring, face, density)
                rho_background = density_value(inputs, density)
                response_cases[key] = package(
                    response_operator_case(inputs, anchoring, face, density),
                    ((1, 0, 0), (1, 1, 0), (1, 0, 1)),
                    sp.Tuple(DIM_PRESSURE, DIM_FLUX, DIM_PRESSURE),
                    CURVED_FACE_SOURCES=sp.Tuple(
                        Str("face_velocity"), Str("kinematic_balance"), Str("relative_flux"),
                        Str("closure_shape_deriv"), Str("traction"), Str("face_normal")
                    ),
                )
                coefficients = closed_coefficients(inputs, z0_out, z0_in, kernel, rho_background)
                response_coeff_cases[key] = package(
                    coefficients,
                    ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
                    sp.Tuple(DIM_IMPEDANCE, DIM_RHO4, DIM_RHO4, DIM_LAMBDA_A),
                    INPUTS=sp.Tuple(Str("V_S"), Str("MU_THETA")),
                )
    objects["FACE_RESPONSE"] = sp.Tuple(
        sp.Tuple(Str("OPAQUE_MU_THETA_OPERATOR"), inputs.rows["mu_theta_operator"]),
        sp.Tuple(Str("CASES"), case_map(response_cases)),
    )
    objects["FACE_RESPONSE_COEFFS"] = case_map(response_coeff_cases)

    loci = locus_package(inputs)
    noninvert_cases = {}
    lambda_A, _, _ = lambda_kernels(inputs)
    for anchoring in ANCHORINGS:
        for face in FACES:
            z_operator = z_case_operators[(anchoring, face)]
            response_operator = identity_operator + lambda_A * z_operator / inputs.rho_m**2
            zero_in_spectrum = sp.Eq(
                sp.Function("S11CC1ZeroInSpectrum")(response_operator), sp.Integer(1), evaluate=False
            )
            noninvert_cases[(anchoring, face)] = package(
                sp.Tuple(
                    sp.Tuple(Str("OPERATOR"), response_operator),
                    sp.Tuple(Str("NONINVERTIBILITY"), zero_in_spectrum),
                    sp.Tuple(Str("FLAT_DIAGONAL_SYMBOL"), loci["DIAGONAL_DENOMINATOR"]),
                    sp.Tuple(
                        Str("PROFILE_CONDITIONED_RESOLVENT"), Str("RESERVED_FOR_S11C_D")
                    ),
                ),
                ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
                DIM_ZERO,
            )
    objects["NONINVERTIBILITY_CONDITION"] = case_map(noninvert_cases)
    for name in (
        "EQUATIONS", "SOLUTION", "IDENTICALLY_SATISFIED", "INCONSISTENT", "REAL_ADMISSIBLE"
    ):
        objects[f"DEGENERATE_LOCI_{name}"] = loci[name]

    port_cases = {}
    dissipation_limits = {}
    for anchoring in ANCHORINGS:
        bare_dtn = z0_out * delta_k + direct_kernels[(anchoring, 1)]
        hermitian = port_hermitian(inputs, bare_dtn)
        for parity in PARITIES:
            for output_regime in REGIMES:
                for input_regime in REGIMES:
                    key = (anchoring, parity, output_regime, input_regime)
                    regime_h = hermitian.applyfunc(
                        lambda entry: substitute_regime(
                            entry,
                            output_regime,
                            input_regime,
                            inputs,
                            evaluate_grazing=False,
                        )
                    )
                    determinant = sp.Add(
                        regime_h[0, 0] * regime_h[1, 1],
                        -regime_h[0, 1] * regime_h[1, 0],
                        evaluate=False,
                    )
                    sign_tests = sp.Tuple(
                        sp.Ge(sp.re(regime_h[0, 0], evaluate=False), 0, evaluate=False),
                        sp.Ge(sp.re(determinant, evaluate=False), 0, evaluate=False),
                    )
                    nondegenerate = output_regime == input_regime == "PROPAGATING"
                    token = Str("SIGN_TEST_ON_NONDEGENERATE_FLAT_SUBSPACE") if nondegenerate else Str(
                        "NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER"
                    )
                    form_fields = [
                        sp.Tuple(Str("FULL_BARE_DTN"), bare_dtn),
                        sp.Tuple(Str("BLOCK_HERMITIAN_FORM"), regime_h),
                        sp.Tuple(Str("STATUS_TOKEN"), token),
                        sp.Tuple(Str("PAIRING_MEASURE_SOURCE"), Str("face_measure_shape_deriv")),
                    ]
                    if nondegenerate:
                        form_fields.insert(2, sp.Tuple(Str("SIGN_TEST_OBJECTS"), sign_tests))
                    port_cases[key] = package(
                        sp.Tuple(*form_fields),
                        ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
                        sp.Tuple(DIM_PRESSURE, DIM_FLUX),
                    )
        normal_lambdas = dict(zip(("A", "V", "X"), lambda_kernels(inputs)))
        for channel, zero_value in (
            ("A", inputs.lambda_A0), ("V", inputs.lambda_V0), ("X", inputs.lambda_X0)
        ):
            for limit_name, replacement in (("OMEGA_TAU_TO_ZERO", zero_value), ("OMEGA_TAU_TO_INFINITY", sp.Integer(0))):
                changed = dict(normal_lambdas)
                changed[channel] = replacement
                changed_tuple = (changed["A"], changed["V"], changed["X"])
                rebuilt = port_hermitian(inputs, bare_dtn, changed_tuple)
                for parity in PARITIES:
                    for output_regime in REGIMES:
                        for input_regime in REGIMES:
                            dissipation_limits[(anchoring, parity, output_regime, input_regime, channel, limit_name)] = package(
                                rebuilt.applyfunc(
                                    lambda entry: substitute_regime(
                                        entry,
                                        output_regime,
                                        input_regime,
                                        inputs,
                                        evaluate_grazing=False,
                                    )
                                ),
                                ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
                                sp.Tuple(DIM_PRESSURE, DIM_FLUX),
                            )
    objects["PERMEABLE_PORT_HERMITIAN"] = case_map(port_cases)
    objects["PERMEABLE_DISSIPATION_VS_OMEGA_TAU"] = case_map(dissipation_limits)

    energy_face = {}
    energy_bulk = {}
    energy_residual = {}
    qprop = regime_replacements(inputs)["PROPAGATING"][0]
    for anchoring in ANCHORINGS:
        for face in FACES:
            velocity = velocity_inputs[(anchoring, face)]
            source = shape_source(inputs, face=face)
            z0_propagating = dtn_flat_symbol(inputs, qprop)
            z1_propagating = direct_kernels[(anchoring, face)].xreplace(
                {q_out_k: qprop, q_out_kp: qprop}
            )
            impermeable_traction = closed_coefficients(
                inputs,
                z0_propagating,
                z0_propagating,
                z1_propagating,
                inputs.rho_br,
                (sp.Integer(0), sp.Integer(0), sp.Integer(0)),
            )["TRACTION_NORMAL_AMPLITUDE"]
            traction_v0 = impermeable_traction[0][0]
            traction_v1 = impermeable_traction[1][0]
            true_area_weight = first_shape_true_area_weight(source)
            traction_pairing = (
                -true_area_weight
                * (traction_v0 + traction_v1)
                * velocity
                * sp.conjugate(velocity)
            )
            face_power = sp.re(traction_pairing, evaluate=False) / 2
            outgoing_flux, farfield_construction = outgoing_farfield_poynting(
                inputs, source, qprop, qprop, velocity
            )
            branch_reversed_flux, branch_reversed_construction = outgoing_farfield_poynting(
                inputs, source, -qprop, -qprop, velocity
            )
            bulk_comparison = -outgoing_flux
            branch_reversed_comparison = -branch_reversed_flux
            key = (anchoring, face, "REAL_OMEGA_PROPAGATING_IMPERMEABLE_LAMBDA_X0_ZERO")
            energy_face[key] = package(
                sp.Tuple(
                    sp.Tuple(Str("BASELINE"), face_power),
                    sp.Tuple(Str("TRACTION_SIGN_REVERSED"), -face_power),
                    sp.Tuple(Str("FARFIELD_BRANCH_SIGN_REVERSED"), face_power),
                ),
                ((2, 0, 0), (2, 1, 0), (2, 0, 1)),
                dim_add(DIM_PRESSURE, DIM_VELOCITY),
                TRACTION_FROM_CLOSED_RESPONSE=casify(impermeable_traction),
                TRUE_AREA_WEIGHT=true_area_weight,
                MEASURE_SOURCE=Str("face_measure_shape_deriv"),
            )
            energy_bulk[key] = package(
                sp.Tuple(
                    sp.Tuple(
                        Str("BASELINE"),
                        sp.Tuple(
                            Str("OUTGOING_POYNTING_FLUX"),
                            outgoing_flux,
                            Str("COMPARISON_ORIENTATION"),
                            bulk_comparison,
                        ),
                    ),
                    sp.Tuple(
                        Str("TRACTION_SIGN_REVERSED"),
                        sp.Tuple(
                            Str("OUTGOING_POYNTING_FLUX"),
                            outgoing_flux,
                            Str("COMPARISON_ORIENTATION"),
                            bulk_comparison,
                        ),
                    ),
                    sp.Tuple(
                        Str("FARFIELD_BRANCH_SIGN_REVERSED"),
                        sp.Tuple(
                            Str("OUTGOING_POYNTING_FLUX"),
                            branch_reversed_flux,
                            Str("COMPARISON_ORIENTATION"),
                            branch_reversed_comparison,
                        ),
                    ),
                ),
                ((2, 0, 0), (2, 1, 0), (2, 0, 1)),
                dim_add(DIM_PRESSURE, DIM_VELOCITY),
                OUTGOING_PHI_CONSTRUCTION=farfield_construction,
                BRANCH_REVERSED_PHI_CONSTRUCTION=branch_reversed_construction,
            )
            energy_residual[key] = package(
                sp.Tuple(
                    sp.Tuple(Str("BASELINE"), sp.simplify(face_power - bulk_comparison)),
                    sp.Tuple(
                        Str("TRACTION_SIGN_REVERSED"),
                        sp.simplify(-face_power - bulk_comparison),
                    ),
                    sp.Tuple(
                        Str("FARFIELD_BRANCH_SIGN_REVERSED"),
                        sp.simplify(face_power - branch_reversed_comparison),
                    ),
                ),
                ((2, 0, 0), (2, 1, 0), (2, 0, 1)),
                dim_add(DIM_PRESSURE, DIM_VELOCITY),
            )
    objects["ENERGY_FACE_TRACTION_OPERAND"] = case_map(energy_face)
    objects["ENERGY_BULK_FARFIELD_FLUX_OPERAND"] = case_map(energy_bulk)
    objects["ENERGY_RESIDUAL"] = case_map(energy_residual)

    rep_eulerian = {}
    rep_hanzawa = {}
    rep_residual = {}
    for anchoring in ANCHORINGS:
        direct_map = {face: direct_kernels[(anchoring, face)] for face in FACES}
        hanzawa_map = {face: hanzawa_kernels[(anchoring, face)] for face in FACES}
        eulerian_operand = {
            face: sp.Tuple(
                sp.Tuple(Str("VALUE"), direct_map[face]),
                sp.Tuple(
                    Str("CONSTRUCTION"),
                    Str("DIRECT_EULERIAN_LEVEL_SET_BOUNDARY_PERTURBATION"),
                    Str("face_normal"),
                    Str("conormal_deriv"),
                    Str("face_shift"),
                ),
            )
            for face in FACES
        }
        hanzawa_operand = {
            face: sp.Tuple(
                sp.Tuple(Str("VALUE"), hanzawa_map[face]),
                sp.Tuple(Str("CONSTRUCTION"), hanzawa_constructions[(anchoring, face)]),
            )
            for face in FACES
        }
        rep_eulerian[("DTN", anchoring)] = case_map(eulerian_operand)
        rep_hanzawa[("DTN", anchoring)] = case_map(hanzawa_operand)
        rep_residual[("DTN", anchoring)] = simplified_object_difference(
            direct_map, hanzawa_map
        )

        one_sided_eulerian = {
            face: dtn_first_kernel(
                inputs,
                shape_source(inputs, face=face, flip_upper_x1=True),
                q_out_k,
                q_out_kp,
            )
            for face in FACES
        }
        rep_eulerian[("DTN_EULERIAN_ONE_SIDED_CORRUPTION", anchoring)] = case_map(
            one_sided_eulerian
        )
        rep_hanzawa[("DTN_EULERIAN_ONE_SIDED_CORRUPTION", anchoring)] = case_map(
            hanzawa_map
        )
        rep_residual[("DTN_EULERIAN_ONE_SIDED_CORRUPTION", anchoring)] = (
            simplified_object_difference(one_sided_eulerian, hanzawa_map)
        )
        for density in DENSITIES:
            rho_background = density_value(inputs, density)
            direct_response = {
                face: closed_coefficients(inputs, z0_out, z0_in, direct_map[face], rho_background)
                for face in FACES
            }
            hanzawa_response = {
                face: closed_coefficients(inputs, z0_out, z0_in, hanzawa_map[face], rho_background)
                for face in FACES
            }
            rep_eulerian[("FACE_RESPONSE", anchoring, density)] = case_map(direct_response)
            rep_hanzawa[("FACE_RESPONSE", anchoring, density)] = case_map(hanzawa_response)
            rep_residual[("FACE_RESPONSE", anchoring, density)] = simplified_object_difference(
                direct_response, hanzawa_response
            )
            one_sided_response = {
                face: closed_coefficients(
                    inputs,
                    z0_out,
                    z0_in,
                    one_sided_eulerian[face],
                    rho_background,
                )
                for face in FACES
            }
            corruption_key = (
                "FACE_RESPONSE_EULERIAN_ONE_SIDED_CORRUPTION",
                anchoring,
                density,
            )
            rep_eulerian[corruption_key] = case_map(one_sided_response)
            rep_hanzawa[corruption_key] = case_map(hanzawa_response)
            rep_residual[corruption_key] = simplified_object_difference(
                one_sided_response, hanzawa_response
            )
    objects["REP_INVARIANCE_EULERIAN_OPERAND"] = case_map(rep_eulerian)
    objects["REP_INVARIANCE_HANZAWA_OPERAND"] = case_map(rep_hanzawa)
    objects["REP_INVARIANCE_RESIDUAL"] = case_map(rep_residual)

    independence_base = {}
    independence_corrupt = {}
    independence_residual = {}
    for anchoring in ANCHORINGS:
        base_map = {}
        corrupt_map = {}
        for face in FACES:
            base_source = shape_source(inputs, face=face)
            corrupt_source = shape_source(inputs, face=face, flip_upper_x1=True)
            base_map[face] = dtn_first_kernel(inputs, base_source, q_out_k, q_out_kp)
            corrupt_map[face] = dtn_first_kernel(inputs, corrupt_source, q_out_k, q_out_kp)
        independence_base[("DTN", anchoring)] = case_map(base_map)
        independence_corrupt[("DTN", anchoring)] = case_map(corrupt_map)
        independence_residual[("DTN", anchoring)] = object_difference(base_map, corrupt_map)
        for density in DENSITIES:
            rho_background = density_value(inputs, density)
            base_response = {
                face: closed_coefficients(inputs, z0_out, z0_in, base_map[face], rho_background)
                for face in FACES
            }
            corrupt_response = {
                face: closed_coefficients(inputs, z0_out, z0_in, corrupt_map[face], rho_background)
                for face in FACES
            }
            independence_base[("FACE_RESPONSE", anchoring, density)] = case_map(base_response)
            independence_corrupt[("FACE_RESPONSE", anchoring, density)] = case_map(corrupt_response)
            independence_residual[("FACE_RESPONSE", anchoring, density)] = object_difference(
                base_response, corrupt_response
            )
    objects["CONTROL_INDEPENDENCE_BASE_OPERAND"] = case_map(independence_base)
    objects["CONTROL_INDEPENDENCE_CORRUPTED_OPERAND"] = case_map(independence_corrupt)
    objects["CONTROL_INDEPENDENCE_RESIDUAL"] = case_map(independence_residual)

    form_base = {}
    form_ablated = {}
    form_residual = {}
    for anchoring in ANCHORINGS:
        for density in DENSITIES:
            rho_background = density_value(inputs, density)
            for direction in DIRECTIONS:
                base_map = {}
                ablated_map = {}
                base_response = {}
                ablated_response = {}
                for face in FACES:
                    base_kernel = dtn_first_kernel(
                        inputs, shape_source(inputs, face=face), q_out_k, q_out_kp
                    )
                    ablated_kernel = dtn_first_kernel(
                        inputs, shape_source(inputs, face=face, ablate=direction), q_out_k, q_out_kp
                    )
                    base_map[face] = base_kernel
                    ablated_map[face] = ablated_kernel
                    base_response[face] = closed_coefficients(
                        inputs, z0_out, z0_in, base_kernel, rho_background
                    )
                    ablated_response[face] = closed_coefficients(
                        inputs, z0_out, z0_in, ablated_kernel, rho_background
                    )
                for object_name, base_object, ablated_object in (
                    ("DTN", base_map, ablated_map),
                    ("FACE_RESPONSE", base_response, ablated_response),
                ):
                    key = (object_name, anchoring, density, direction)
                    form_base[key] = case_map(base_object)
                    form_ablated[key] = case_map(ablated_object)
                    form_residual[key] = object_difference(base_object, ablated_object)
    objects["CONTROL_FORM_BASE_OPERAND"] = case_map(form_base)
    objects["CONTROL_FORM_ABLATED_OPERAND"] = case_map(form_ablated)
    objects["CONTROL_FORM_RESIDUAL"] = case_map(form_residual)

    c1_flat_response, c1_flat_coeffs = rederive_flat_face_response(inputs)
    uniform_c1_package = {
        "Z_IMPERMEABLE": sp.Tuple(
            dtn_flat_symbol(inputs, inputs.q_out), dtn_flat_symbol(inputs, inputs.q_out)
        ),
        "Z_BY_REGIME": sp.Tuple(
            sp.Eq(sp.im(inputs.omega), 0, evaluate=False),
            sp.Tuple(Str("PROPAGATING"), dtn_flat_symbol(inputs, inputs.q_out)),
            sp.Tuple(Str("EVANESCENT"), dtn_flat_symbol(inputs, inputs.q_out)),
            sp.Tuple(Str("GRAZING"), sp.limit(dtn_flat_symbol(inputs, q_out_k), q_out_k, 0, dir="+")),
        ),
        "Z_BY_PARITY": sp.diag(
            dtn_flat_symbol(inputs, inputs.q_out), dtn_flat_symbol(inputs, inputs.q_out)
        ),
        "ADDED_MASS": m_add,
        "GRAZING_BEHAVIOUR": sp.limit(dtn_flat_symbol(inputs, q_out_k), q_out_k, 0, dir="+"),
        "FACE_RESPONSE": c1_flat_response,
        "FACE_RESPONSE_COEFFS": c1_flat_coeffs,
        "PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY": objects["PERMEABLE_PORT_HERMITIAN"],
        "DEGENERATE_LOCI_EQUATIONS": loci["EQUATIONS"],
        "DEGENERATE_LOCI_SOLUTION": loci["SOLUTION"],
        "DEGENERATE_LOCI_IDENTICALLY_SATISFIED": loci["IDENTICALLY_SATISFIED"],
        "DEGENERATE_LOCI_INCONSISTENT": loci["INCONSISTENT"],
        "DEGENERATE_LOCI_REAL_ADMISSIBLE": loci["REAL_ADMISSIBLE"],
    }
    uniform_s11b_package = {
        "Z_IMPERMEABLE": inputs.rows["z_impermeable"],
        "Z_BY_REGIME": inputs.rows["z_by_regime"],
        "Z_BY_PARITY": inputs.rows["z_by_parity"],
        "ADDED_MASS": inputs.rows["added_mass"],
        "GRAZING_BEHAVIOUR": inputs.rows["grazing_behaviour"],
        "FACE_RESPONSE": inputs.rows["face_response"],
        "FACE_RESPONSE_COEFFS": inputs.rows["face_response_coeffs"],
        "PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY": inputs.rows[
            "permeable_dissipative_by_regime_and_parity"
        ],
        "DEGENERATE_LOCI_EQUATIONS": inputs.rows["degenerate_loci_equations"],
        "DEGENERATE_LOCI_SOLUTION": inputs.rows["degenerate_loci_solution"],
        "DEGENERATE_LOCI_IDENTICALLY_SATISFIED": inputs.rows[
            "degenerate_loci_identically_satisfied"
        ],
        "DEGENERATE_LOCI_INCONSISTENT": inputs.rows["degenerate_loci_inconsistent"],
        "DEGENERATE_LOCI_REAL_ADMISSIBLE": inputs.rows["degenerate_loci_real_admissible"],
    }
    uniform_c1 = {}
    uniform_s11b = {}
    uniform_residual = {}
    for object_name in uniform_c1_package:
        for output_regime in REGIMES:
            for input_regime in REGIMES:
                for parity in PARITIES:
                    key = (object_name, output_regime, input_regime, parity)
                    uniform_c1[key] = uniform_c1_package[object_name]
                    uniform_s11b[key] = uniform_s11b_package[object_name]
                    uniform_residual[key] = object_difference(
                        uniform_c1_package[object_name], uniform_s11b_package[object_name]
                    )
    objects["UNIFORM_LIMIT_S11CC1_OPERAND"] = case_map(uniform_c1)
    objects["UNIFORM_LIMIT_S11B_OPERAND"] = case_map(uniform_s11b)
    objects["UNIFORM_LIMIT_RESIDUAL"] = case_map(uniform_residual)

    constant_height = sp.Symbol("s11cc1_zero_jet_constant_height")
    zero_jet_source = {"height_hat": constant_height, "tilt_hat": (sp.Integer(0),) * 3}
    zero_jet_first = dtn_first_kernel(inputs, zero_jet_source, q_out_k, q_out_kp)
    zero_jet_diagonal = zero_jet_first.xreplace(
        {q_out_kp: q_out_k, **{b: a for a, b in zip(k_out, k_in)}}
    )
    zero_jet_diagonal = reduce_on_shell(
        zero_jet_diagonal, inputs, q_out_k, k_out
    )
    zero_jet_value = sp.Tuple(
        z0_out + sp.factor(zero_jet_diagonal), z0_out + sp.factor(zero_jet_diagonal)
    ).xreplace({q_out_k: inputs.q_out})
    zero_jet_target = inputs.rows["z_impermeable"]
    objects["ZERO_JET_S11CC1_OPERAND"] = package(
        zero_jet_value, ((0, 0, 0), (0, 1, 0)), DIM_IMPEDANCE,
        CONSTANT_PROFILE=sp.Tuple(inputs.eta, inputs.w1_profile),
    )
    objects["ZERO_JET_S11B_OPERAND"] = package(
        zero_jet_target, ((0, 0, 0),), DIM_IMPEDANCE
    )
    objects["ZERO_JET_RESIDUAL"] = package(
        object_difference(zero_jet_value, zero_jet_target),
        ((0, 0, 0), (0, 1, 0)),
        DIM_IMPEDANCE,
    )

    branch_base = {}
    branch_corruptions = {name: {} for name in (
        "SIGNFLIP_INPUT", "SIGNFLIP_OUTPUT", "MOMENTUMFREEZE_OUTPUT", "MOMENTUMFREEZE_INPUT"
    )}
    branch_residual = {}
    q_mutations = {
        "SIGNFLIP_INPUT": (q_out_k, -q_out_kp),
        "SIGNFLIP_OUTPUT": (-q_out_k, q_out_kp),
        "MOMENTUMFREEZE_OUTPUT": (q_out_kp, q_out_kp),
        "MOMENTUMFREEZE_INPUT": (q_out_k, q_out_k),
    }
    for anchoring in ANCHORINGS:
        base_kernel = dtn_first_kernel(
            inputs, shape_source(inputs, face=1), q_out_k, q_out_kp
        )
        branch_base[("DTN", anchoring)] = base_kernel
        for mutation, (q_output, q_input) in q_mutations.items():
            corrupted_kernel = dtn_first_kernel(
                inputs, shape_source(inputs, face=1), q_output, q_input
            )
            branch_corruptions[mutation][("DTN", anchoring)] = corrupted_kernel
            branch_residual[(mutation, "DTN", anchoring)] = object_difference(
                base_kernel, corrupted_kernel
            )
        for density in DENSITIES:
            rho_background = density_value(inputs, density)
            base_response = closed_coefficients(
                inputs, z0_out, z0_in, base_kernel, rho_background
            )
            branch_base[("FACE_RESPONSE", anchoring, density)] = base_response
            for mutation, (q_output, q_input) in q_mutations.items():
                mutated_z0_output = dtn_flat_symbol(inputs, q_output)
                mutated_z0_input = dtn_flat_symbol(inputs, q_input)
                corrupted_kernel = dtn_first_kernel(
                    inputs, shape_source(inputs, face=1), q_output, q_input
                )
                corrupted_response = closed_coefficients(
                    inputs, mutated_z0_output, mutated_z0_input, corrupted_kernel, rho_background
                )
                key = ("FACE_RESPONSE", anchoring, density)
                branch_corruptions[mutation][key] = corrupted_response
                branch_residual[(mutation, *key)] = object_difference(
                    base_response, corrupted_response
                )
    branch_base[("RADIATION_REAL_AXIS", "PROPAGATING")] = sp.re(
        inputs.omega * regime_replacements(inputs)["PROPAGATING"][0]
    )
    branch_base[("RADIATION_REAL_AXIS", "EVANESCENT")] = sp.im(
        regime_replacements(inputs)["EVANESCENT"][0]
    )
    objects["CONTROL_BRANCH_BASE_OPERAND"] = case_map(branch_base)
    objects["CONTROL_BRANCH_SIGNFLIP_INPUT_OPERAND"] = case_map(branch_corruptions["SIGNFLIP_INPUT"])
    objects["CONTROL_BRANCH_SIGNFLIP_OUTPUT_OPERAND"] = case_map(branch_corruptions["SIGNFLIP_OUTPUT"])
    objects["CONTROL_BRANCH_MOMENTUMFREEZE_OUTPUT_OPERAND"] = case_map(
        branch_corruptions["MOMENTUMFREEZE_OUTPUT"]
    )
    objects["CONTROL_BRANCH_MOMENTUMFREEZE_INPUT_OPERAND"] = case_map(
        branch_corruptions["MOMENTUMFREEZE_INPUT"]
    )
    objects["CONTROL_BRANCH_RESIDUAL"] = case_map(branch_residual)

    dimensions = {
        "DTN_FLAT_SYMBOL": DIM_IMPEDANCE,
        "DTN_OPERATOR": DIM_IMPEDANCE,
        "DTN_KERNEL": DIM_KERNEL,
        "PRESSURE_OVER_VELOCITY": DIM_IMPEDANCE,
        "PRESSURE_OVER_MU_THETA": dim_sub(DIM_PRESSURE, DIM_MU_THETA),
        "FLUX_OVER_VELOCITY": DIM_RHO4,
        "FLUX_OVER_MU_THETA": dim_sub(DIM_FLUX, DIM_MU_THETA),
        "FACE_RESPONSE": sp.Tuple(DIM_PRESSURE, DIM_FLUX, DIM_PRESSURE),
        "PORT_POWER_DENSITY": dim_add(DIM_PRESSURE, DIM_VELOCITY),
        "LAMBDA_A_0": DIM_LAMBDA_A,
        "LAMBDA_V_0": DIM_LAMBDA_V,
        "LAMBDA_X_0": DIM_LAMBDA_X,
    }
    objects["DIMENSIONS"] = case_map(dimensions)
    height_term = inputs.eta * inputs.W0 * w1_hat * q_out_k * q_out_kp / 2
    laplace_term = inputs.eta * inputs.W0 * w1_hat * (
        sp.Add(*(component**2 for component in k_in)) - kappa_sq
    ) / 2
    tilt_terms = tuple(-sp.I * inputs.sigma_W * jet * momentum / 2 for jet, momentum in zip(w1_jet_hat, k_in))
    base_dimensions = {
        "HEIGHT_NORMAL_NORMAL": dim_add(DIM_L, 3 * DIM_L, -DIM_L, -DIM_L),
        "HEIGHT_TANGENTIAL_HELMHOLTZ": dim_add(DIM_L, 3 * DIM_L, -2 * DIM_L),
        **{
            f"TILT_DIRECTION_{direction}": dim_add(3 * DIM_L, -DIM_L)
            for direction in DIRECTIONS
        },
    }
    control_dimensions = dict(base_dimensions)
    control_dimensions["TILT_DIRECTION_1"] = dim_add(4 * DIM_L, -DIM_L)
    homogeneity_base = {
        name: sp.Tuple(expression, base_dimensions[name])
        for name, expression in {
            "HEIGHT_NORMAL_NORMAL": height_term,
            "HEIGHT_TANGENTIAL_HELMHOLTZ": laplace_term,
            **{f"TILT_DIRECTION_{i}": value for i, value in zip(DIRECTIONS, tilt_terms)},
        }.items()
    }
    homogeneity_control = {
        name: sp.Tuple(homogeneity_base[name][0], control_dimensions[name]) for name in homogeneity_base
    }
    objects["HOMOGENEITY_BASE_OPERAND"] = case_map(homogeneity_base)
    objects["HOMOGENEITY_CONTROL_OPERAND"] = case_map(homogeneity_control)
    objects["HOMOGENEITY_RESIDUAL"] = case_map(
        {
            name: sp.ImmutableMatrix(base_dimensions[name] - control_dimensions[name])
            for name in homogeneity_base
        }
    )

    # These inherited laws are construction sources, not re-derived rows.  They
    # remain live in the model and in the emit-only provenance package.
    objects["SUPPLIED_FACE_LAW_OPERANDS"] = case_map(
        {
            "FACE_VELOCITY": inputs.rows["face_velocity"],
            "RELATIVE_FLUX": inputs.rows["relative_flux"],
            "KINEMATIC_BALANCE": inputs.rows["kinematic_balance"],
            "TRACTION": inputs.rows["traction"],
            "CLOSURE_SHAPE_DERIV": inputs.rows["closure_shape_deriv"],
            "E_W_DOF": inputs.e_W,
        }
    )
    return AuditModel(inputs, MappingProxyType(objects))


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
ISSUES: list[str] = []
OPERATIONAL_FAILURES: list[str] = []
TASK_TIMINGS: dict[str, float] = {}
CURRENT_TASK: str | None = None


def emit(
    quantity: str,
    payload: object,
    *,
    export_key: str | None = None,
    local: bool = False,
    class_tag: str = "DERIVED",
    description: str | None = None,
) -> object:
    infix = "LOCAL_" if local else ""
    return EMITTER.emit(
        f"PY_S11CC1_{infix}{quantity}",
        payload,
        export_key=export_key if not local else None,
        class_tag=class_tag,
        description=description,
    )


def issue(message: str) -> None:
    ISSUES.append(f"{CURRENT_TASK or 'GLOBAL'}: {message}")


def task_dtn(model: AuditModel) -> None:
    for name in (
        "DTN_FLAT_SYMBOL",
        "DTN_OPERATOR",
        "DTN_KERNEL",
        "DTN_RIGID_SHIFT_OPERAND",
        "DTN_RIGID_SHIFT_RESIDUAL",
        "DTN_BY_REGIME_PAIR",
        "DTN_BY_PARITY",
        "DTN_HERMITIAN_PART",
        "DTN_REACTIVE_PART",
        "DTN_INERTIAL_LOADING",
        "DTN_GRAZING_BEHAVIOUR",
        "DTN_TERM_ORIGINS",
    ):
        export_key = {
            "DTN_FLAT_SYMBOL": "dtn_flat_symbol",
            "DTN_OPERATOR": "dtn_operator",
            "DTN_KERNEL": "dtn_kernel",
        }.get(name)
        emit(name, model.objects[name], export_key=export_key)


def task_face_response(model: AuditModel) -> None:
    for name in (
        "FACE_RESPONSE",
        "FACE_RESPONSE_COEFFS",
        "NONINVERTIBILITY_CONDITION",
        "DEGENERATE_LOCI_EQUATIONS",
        "DEGENERATE_LOCI_SOLUTION",
        "DEGENERATE_LOCI_IDENTICALLY_SATISFIED",
        "DEGENERATE_LOCI_INCONSISTENT",
        "DEGENERATE_LOCI_REAL_ADMISSIBLE",
        "PERMEABLE_PORT_HERMITIAN",
        "PERMEABLE_DISSIPATION_VS_OMEGA_TAU",
        "ENERGY_FACE_TRACTION_OPERAND",
        "ENERGY_BULK_FARFIELD_FLUX_OPERAND",
        "ENERGY_RESIDUAL",
    ):
        export_key = {
            "FACE_RESPONSE": "face_response",
            "FACE_RESPONSE_COEFFS": "face_response_coeffs",
        }.get(name)
        emit(name, model.objects[name], export_key=export_key)


def task_rep_invariance(model: AuditModel) -> None:
    for name in (
        "REP_INVARIANCE_EULERIAN_OPERAND",
        "REP_INVARIANCE_HANZAWA_OPERAND",
        "REP_INVARIANCE_RESIDUAL",
    ):
        emit(name, model.objects[name])


def task_independence(model: AuditModel) -> None:
    for name in (
        "CONTROL_INDEPENDENCE_BASE_OPERAND",
        "CONTROL_INDEPENDENCE_CORRUPTED_OPERAND",
        "CONTROL_INDEPENDENCE_RESIDUAL",
    ):
        emit(name, model.objects[name])


def task_form(model: AuditModel) -> None:
    for name in (
        "CONTROL_FORM_BASE_OPERAND",
        "CONTROL_FORM_ABLATED_OPERAND",
        "CONTROL_FORM_RESIDUAL",
    ):
        emit(name, model.objects[name])


def task_uniform(model: AuditModel) -> None:
    for name in (
        "UNIFORM_LIMIT_S11CC1_OPERAND",
        "UNIFORM_LIMIT_S11B_OPERAND",
        "UNIFORM_LIMIT_RESIDUAL",
    ):
        emit(name, model.objects[name])


def task_zero_jet(model: AuditModel) -> None:
    for name in ("ZERO_JET_S11CC1_OPERAND", "ZERO_JET_S11B_OPERAND", "ZERO_JET_RESIDUAL"):
        emit(name, model.objects[name])


def task_branch(model: AuditModel) -> None:
    for name in (
        "CONTROL_BRANCH_BASE_OPERAND",
        "CONTROL_BRANCH_SIGNFLIP_INPUT_OPERAND",
        "CONTROL_BRANCH_SIGNFLIP_OUTPUT_OPERAND",
        "CONTROL_BRANCH_MOMENTUMFREEZE_OUTPUT_OPERAND",
        "CONTROL_BRANCH_MOMENTUMFREEZE_INPUT_OPERAND",
        "CONTROL_BRANCH_RESIDUAL",
    ):
        emit(name, model.objects[name])


def task_homogeneity(model: AuditModel) -> None:
    for name in (
        "DIMENSIONS",
        "HOMOGENEITY_BASE_OPERAND",
        "HOMOGENEITY_CONTROL_OPERAND",
        "HOMOGENEITY_RESIDUAL",
    ):
        emit(name, model.objects[name])


TASK_FUNCTIONS = {
    "DTN": task_dtn,
    "FACE_RESPONSE": task_face_response,
    "REP_INVARIANCE": task_rep_invariance,
    "INDEPENDENCE": task_independence,
    "FORM": task_form,
    "UNIFORM": task_uniform,
    "ZERO_JET": task_zero_jet,
    "BRANCH": task_branch,
    "HOMOGENEITY": task_homogeneity,
}


def run_task(name: str, model: AuditModel) -> list[str]:
    global CURRENT_TASK
    started = time.monotonic()
    CURRENT_TASK = name
    try:
        TASK_FUNCTIONS[name](model)
        return [name]
    except Exception as exc:
        OPERATIONAL_FAILURES.append(f"{name}: {type(exc).__name__}: {exc}")
        return []
    finally:
        TASK_TIMINGS[name] = time.monotonic() - started
        CURRENT_TASK = None


def atoms_of(value: object) -> set[sp.Basic]:
    value = casify(value)
    atoms: set[sp.Basic] = set(value.atoms(sp.Symbol))
    atoms.update(value.atoms(sp.MatrixSymbol))
    return atoms


def make_record(candidate: CandidateRow) -> dict[str, object]:
    record: dict[str, object] = {
        "display": display(candidate.value),
        "value": casify(candidate.value),
        "value_kind": "COMPUTED_OBJECT",
        "class": candidate.class_tag,
        "step": "S11c-c1",
    }
    if candidate.description is not None:
        record["description"] = candidate.description
    return record


def intended_own_bind_closure(candidates: list[CandidateRow]) -> frozenset[str]:
    candidate_roots = {candidate.key for candidate in candidates}
    if candidate_roots != INTENDED_EXPORT_CANDIDATE_ROOTS:
        raise RuntimeError(
            "export candidate roots differ from intended roots: "
            + repr(
                sorted(
                    candidate_roots.symmetric_difference(INTENDED_EXPORT_CANDIDATE_ROOTS)
                )
            )
        )
    root_values = tuple(candidate.value for candidate in candidates)
    used_atoms = set().union(*(atoms_of(value) for value in root_values))
    introduced_root_atoms = {
        atom
        for atom in used_atoms.intersection(NEW_SYMBOLS)
        if atom.name.startswith("s11cc1_")
    }
    return frozenset(
        INTENDED_EXPORT_WRITE_ROOTS
        | {atom.name for atom in introduced_root_atoms}
    )


def build_delta(
    fold: Mapping[str, Mapping[str, object]],
) -> tuple[dict[str, dict[str, object]], dict[str, object], frozenset[str]]:
    candidates = list(EMITTER.export_candidates)
    own_closure = intended_own_bind_closure(candidates)
    root_values = tuple(candidate.value for candidate in candidates)
    used_atoms = set().union(*(atoms_of(value) for value in root_values))
    existing_candidate_names = {candidate.key for candidate in candidates}
    for atom in sorted(used_atoms.intersection(NEW_SYMBOLS), key=sp.default_sort_key):
        if atom.name in existing_candidate_names:
            continue
        metadata = NEW_SYMBOLS[atom]
        candidates.append(
            CandidateRow(
                atom.name,
                atom,
                str(metadata["class"]),
                "NEW_SYMBOL_BIND_CLOSURE",
                str(metadata["description"]),
            )
        )
        existing_candidate_names.add(atom.name)

    rows: dict[str, dict[str, object]] = {}
    routing = []
    f9c = []
    key_operands = []
    resolved_keys: set[str] = set()
    for candidate in candidates:
        imported_spellings = sp.Tuple(Str(candidate.key)) if candidate.key in fold else sp.Tuple()
        key_operands.append(
            sp.Tuple(Str(candidate.source_tag), Str(candidate.key), imported_spellings)
        )
        record = make_record(candidate)
        if candidate.key not in fold:
            write_key = candidate.key
            status = "ABSENT"
            operands = sp.Tuple(record["value"], sp.Tuple())
            record["route"] = "F9A_ABSENT"
        else:
            imported_record = fold[candidate.key]
            value_status, operands = total_compare(record["value"], imported_record.get("value"))
            class_status = (
                PROVED_EQUAL if record["class"] == imported_record.get("class") else PROVED_DIFFERENT
            )
            if value_status == PROVED_EQUAL and class_status == PROVED_EQUAL:
                write_key = candidate.key
                status = PROVED_EQUAL
                record["f9_operands"] = sp.Tuple(record["value"], imported_record.get("value"))
                record["route"] = "F9B_PROVED_EQUAL"
                if "dimension_key" in imported_record:
                    record["dimension_key"] = imported_record["dimension_key"]
                prior_steps = list(imported_record.get("corroborated_steps", ()))
                if not prior_steps and imported_record.get("step"):
                    prior_steps.append(imported_record["step"])
                prior_steps.append("S11c-c1")
                record["corroborated_steps"] = tuple(dict.fromkeys(prior_steps))
            else:
                write_key = f"s11c_c1_{candidate.key}"
                status = value_status if value_status != PROVED_EQUAL else PROVED_DIFFERENT
                record["f9_operands"] = sp.Tuple(record["value"], imported_record.get("value"))
                record["route"] = f"F9C_{status}"
                f9c.append(sp.Tuple(Str(candidate.key), Str(write_key), Str(status), casify(operands)))
                issue(f"F9c write {write_key} for {candidate.key}: {status}")
                if write_key in fold:
                    raise RuntimeError(f"F9c write key collides with imported key {write_key!r}")
        if write_key in rows:
            raise RuntimeError(f"two candidates route to {write_key!r}")
        rows[write_key] = record
        resolved_keys.add(write_key)
        routing.append(
            sp.Tuple(
                Str(candidate.source_tag), Str(candidate.key), Str(write_key), Str(status), casify(operands)
            )
        )

    missing_roots = INTENDED_EXPORT_WRITE_ROOTS.difference(resolved_keys)
    if missing_roots:
        raise RuntimeError("missing resolved bind roots: " + ", ".join(sorted(missing_roots)))
    diagnostics = {
        "KEY_OPERANDS": sp.Tuple(*key_operands),
        "ROUTING": sp.Tuple(*routing),
        "F9C": sp.Tuple(*f9c),
    }
    return rows, diagnostics, own_closure


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
    digest_paths = (
        SCRIPT_PATH,
        SCRIPT_DIR / "S11b_exports.py",
        SCRIPT_DIR / "S11c_a_exports.py",
        SCRIPT_DIR / "S11c_b_exports.py",
        DIRECTIVE_PATH,
        SCRIPT_DIR / "ledger_fold.py",
    )
    lines = [
        "# S11c_c1_exports.py -- GENERATED by S11c_c1_bulk_closure_sympy_audit.py. Do not edit.",
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
        "IMPORT_KEYS = (",
    ]
    lines.extend(f"    {key!r}," for key in IMPORT_KEYS)
    lines.extend([")", "", "BUILD_INPUT_DIGESTS = MappingProxyType({"])
    lines.extend(f"    {path.name!r}: {sha256(path)!r}," for path in digest_paths)
    lines.extend(["})", "", "_LEDGER = {"])
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


def publish_export(ledger: dict[str, dict[str, object]]) -> sp.Tuple:
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
        rows.append(
            sp.Tuple(
                Str(name), Str(str(live_record["display"])), Str(status), casify(residual), casify(operands)
            )
        )
        if status != PROVED_EQUAL:
            failures.append(name)
    if failures:
        raise RuntimeError("export reconstruction not proved equal for " + ", ".join(failures))
    return sp.Tuple(*rows)


def remove_stale_export() -> None:
    if EXPORT_PATH.exists():
        EXPORT_PATH.unlink()


def assertion_package(operand_a: object, operand_b: object, residual: object) -> sp.Tuple:
    return sp.Tuple(
        sp.Tuple(Str("OPERAND_A"), casify(operand_a)),
        sp.Tuple(Str("OPERAND_B"), casify(operand_b)),
        sp.Tuple(Str("RESIDUAL"), casify(residual)),
    )


def run() -> None:
    start = time.monotonic()
    remove_stale_export()
    fold, fold_audit = load_model(BASE_PATH)
    consumer_audit = check_consumer(fold, IMPORT_KEYS)
    emit(
        "FOLD_AUDIT",
        sp.Tuple(
            sp.Tuple(Str("SOURCE_ROW_COUNTS"), casify(fold_audit["source_row_counts"])),
            sp.Tuple(Str("OVERWRITES"), casify(fold_audit["overwrites"])),
        ),
        local=True,
    )
    emit(
        "CONSUMER_CLOSURE",
        sp.Tuple(
            sp.Tuple(Str("CLOSURE_COUNT"), len(consumer_audit["closure"])),
            sp.Tuple(Str("SYMBOL_EDGE_COUNT"), len(consumer_audit["symbol_edges"])),
            sp.Tuple(Str("DIMENSION_EDGE_COUNT"), len(consumer_audit["dimension_edges"])),
        ),
        local=True,
    )

    lookup_audit = assert_lookups_equal_manifest(build_model, fold, IMPORT_KEYS)
    model = lookup_audit["result"]
    observed = tuple(sorted(lookup_audit["lookups"]))
    manifest = tuple(sorted(IMPORT_KEYS))
    lookup_residual = sp.Tuple(
        sp.Tuple(Str("OBSERVED_MINUS_MANIFEST"), casify(sorted(set(observed).difference(manifest)))),
        sp.Tuple(Str("MANIFEST_MINUS_OBSERVED"), casify(sorted(set(manifest).difference(observed)))),
    )
    emit(
        "LOOKUP_MANIFEST_ASSERTION",
        assertion_package(casify(observed), casify(manifest), lookup_residual),
        local=True,
    )

    completed: list[str] = []
    for task in (*PRIMARY_TASKS, *CONTROL_TASKS):
        completed += run_task(task, model)

    f9c_report = sp.Tuple()
    export_started = time.monotonic()
    if all(task in completed for task in PRIMARY_TASKS):
        try:
            delta, diagnostics, own_closure = build_delta(fold)
            emit(
                "EXPORT_CANDIDATE_KEY_OPERANDS",
                sp.Tuple(diagnostics["KEY_OPERANDS"], diagnostics["ROUTING"]),
                local=True,
            )
            f9c_report = diagnostics["F9C"]
            emit("EXPORT_F9C_WRITES", f9c_report, local=True)
            minimal = assert_delta_is_minimal(delta, own_closure, ())
            minimal_residual = sp.Tuple(
                sp.Tuple(
                    Str("ACTUAL_MINUS_REQUIRED"),
                    casify(sorted(minimal["exported_keys"].difference(minimal["required_keys"]))),
                ),
                sp.Tuple(
                    Str("REQUIRED_MINUS_ACTUAL"),
                    casify(sorted(minimal["required_keys"].difference(minimal["exported_keys"]))),
                ),
            )
            emit(
                "DELTA_MINIMAL_ASSERTION",
                assertion_package(
                    casify(sorted(minimal["exported_keys"])),
                    casify(sorted(minimal["required_keys"] | minimal["allowed_infra_keys"])),
                    minimal_residual,
                ),
                local=True,
            )
            roundtrip = publish_export(delta)
            emit("EXPORT_ROUNDTRIP", roundtrip, local=True)
        except Exception as exc:
            OPERATIONAL_FAILURES.append(f"EXPORT: {type(exc).__name__}: {exc}")
            issue(f"export operational failure: {type(exc).__name__}: {exc}")
            remove_stale_export()
    else:
        issue("S11c_c1_exports.py not published because a primary task did not complete")
        remove_stale_export()
    TASK_TIMINGS["EXPORT"] = time.monotonic() - export_started

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
                "PY_S11CC1_LOCAL_TAGS",
                "PY_S11CC1_LOCAL_SECTION8_REPORT",
            )
        )
    )
    emit("TAGS", sp.Tuple(*(Str(name) for name in local_names)), local=True)
    export_lines = len(EXPORT_PATH.read_text(encoding="utf-8").splitlines()) if EXPORT_PATH.exists() else 0
    report = sp.Tuple(
        sp.Tuple(
            Str("FILES_WRITTEN"),
            sp.Tuple(
                Str(str(SCRIPT_PATH)),
                Str(str(EXPORT_PATH)) if EXPORT_PATH.exists() else Str("EXPORT_NOT_PUBLISHED"),
            ),
        ),
        sp.Tuple(Str("SCRIPT_LINES"), len(SCRIPT_PATH.read_text(encoding="utf-8").splitlines())),
        sp.Tuple(Str("EXPORT_LINES"), export_lines),
        sp.Tuple(Str("EMITTED_TAGS"), EMITTER.count + 1),
        sp.Tuple(Str("RUN_TASKS"), run_record),
        sp.Tuple(Str("SKIPPED_TASKS"), skipped_record),
        sp.Tuple(Str("RUNTIME_SECONDS"), sp.Float(time.monotonic() - start, 8)),
        sp.Tuple(
            Str("TAG_NAMES"),
            sp.Tuple(*(Str(name) for name in (*EMITTER.values.keys(), "PY_S11CC1_LOCAL_SECTION8_REPORT"))),
        ),
        sp.Tuple(
            Str("ISSUES_OR_NONCOMPUTABLE"),
            sp.Tuple(*(Str(item) for item in (*ISSUES, *OPERATIONAL_FAILURES))),
        ),
        sp.Tuple(Str("F9C_WRITES"), f9c_report),
        sp.Tuple(
            Str("SUPPLIED_UNFALSIFIABLE"),
            Str("SECTIONS_1_TO_2_AND_S11CA_S11CB_SUBSTRATE"),
        ),
        sp.Tuple(
            Str("REST_FRAME_SCOPE"),
            Str(
                "NONUNIFORM_VALIDITY_DOMAIN_CONDITIONS_EVERY_RESULT; GRAZING_IS_STRICT_V_BULK_NORMAL_0_ZERO"
            ),
        ),
    )
    emit("SECTION8_REPORT", report, local=True)
    if OPERATIONAL_FAILURES:
        raise RuntimeError("; ".join(OPERATIONAL_FAILURES))


if __name__ == "__main__":
    run()
