#!/usr/bin/env python3
"""S11 independent SymPy audit, derived only from S11_SHARED_PHYSICS.md."""

from __future__ import annotations

import itertools
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import sympy as sp
from sympy.core.symbol import Str
from sympy.polys.matrices import DomainMatrix


SCRIPT_PATH = Path(__file__).resolve()
REDUCTION_DIR = SCRIPT_PATH.parent.parent / "reduction"
sys.path.insert(0, str(REDUCTION_DIR))
from registry_read import RegistryValidationError, load_registry  # noqa: E402


D_symbol = sp.Symbol("D", integer=True, positive=True)
rho_br = sp.Symbol("rho_br", real=True, positive=True)
mu_R = sp.Symbol("mu_R", real=True, positive=True)
B_comp = sp.Symbol("B_comp", real=True, positive=True)
mu_br = sp.Symbol("mu_br", real=True, positive=True)
beta = sp.Symbol("beta", real=True)
s_scale = sp.Symbol("s", real=True, positive=True)
c_s0 = sp.Symbol("c_s0", real=True, positive=True)
omegaSquared = sp.Symbol("omegaSquared")
lambdaScale = sp.Symbol("lambdaScale", real=True, positive=True)

ENERGY_DENSITY_DIM = sp.Tuple(2 - D_symbol, -2, 1)
U_DIMENSION = sp.Tuple(1, 0, 0)

PACKAGE_SWEEPS: dict[str, tuple[int, ...]] = {
    "MAIN": (2, 3, 4, 5),
    "XFORM_CURLONLY": (2, 3, 4, 5),
    "XFORM_EXTRA": (2, 3, 4, 5),
    "XFORM_DIVONLY": (3, 4),
    "XFORM_TRACELESS": (3, 4),
    "XCOEF_BSCALE": (3,),
    "XCOEF_BSIGN": (3,),
}

emitted_names: set[str] = set()
local_names: list[str] = []


def emit(name: str, payload: Any) -> None:
    if name in emitted_names:
        raise RuntimeError(f"duplicate tag name: {name}")
    emitted_names.add(name)
    if "PY_S11_LOCAL_" in name:
        local_names.append(name)
    if isinstance(payload, str):
        rendered = sp.srepr(Str(payload))
    elif ((isinstance(payload, sp.Basic) and payload.atoms(Str))
          or beta in getattr(payload, "free_symbols", set())):
        rendered = sp.srepr(payload)
    else:
        rendered = str(payload)
        if "\n" in rendered or "\r" in rendered:
            rendered = sp.srepr(payload)
    print(f"{name}: {rendered}", flush=True)


def tag(package: str, dimension: int, quantity: str, *, stratum: int | None = None,
        root: int | None = None, local: bool = False) -> str:
    head = "PY_S11_LOCAL" if local else "PY_S11"
    parts = [head, package, f"D{dimension}"]
    if stratum is not None:
        parts.append(f"STRATUM{stratum}")
    if root is not None:
        parts.append(f"ROOT{root}")
    parts.append(quantity)
    return "_".join(parts)


def canonical_solver_object(value: Any) -> Any:
    if isinstance(value, dict):
        return sp.Dict({canonical_solver_object(key): canonical_solver_object(item)
                        for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return sp.Tuple(*(canonical_solver_object(item) for item in value))
    if value is True:
        return sp.true
    if value is False:
        return sp.false
    if isinstance(value, str):
        return Str(value)
    return value


@dataclass(frozen=True)
class Q9Result:
    payloads: tuple[tuple[str, Any], ...]
    g_symbols: tuple[sp.Symbol, ...]
    monomials: tuple[sp.Expr, ...]
    v6_matrix: sp.Matrix
    pd_form: sp.Expr


@dataclass(frozen=True)
class ActionData:
    coordinates: tuple[sp.Symbol, ...]
    time: sp.Symbol
    fields: tuple[sp.Expr, ...]
    gradient: sp.Matrix
    densities: Mapping[str, sp.Expr]
    stiffness_records: tuple[tuple[str, sp.Expr, sp.Expr], ...]
    stiffness_terms: tuple[sp.Expr, ...]
    stiffness: sp.Expr
    kinetic: sp.Expr
    lagrangian: sp.Expr
    action_terms: tuple[sp.Expr, ...]
    coefficient_factors: tuple[sp.Expr, ...]
    eom: tuple[sp.Expr, ...]
    pd_term: sp.Expr
    configuration: Mapping[str, Any]


@dataclass(frozen=True)
class AssumptionData:
    joint: sp.Expr
    common_parts: tuple[sp.Expr, ...]
    package_parts: tuple[sp.Expr, ...]
    dimension_declarations: tuple["DimensionPremise", ...]


@dataclass(frozen=True)
class DimensionPremise:
    symbol: sp.Symbol
    vector: sp.Tuple


@dataclass(frozen=True)
class PhaseAverage:
    variable: sp.Symbol
    lower: sp.Expr
    upper: sp.Expr
    normalization: sp.Expr

    def apply(self, expression: sp.Expr) -> sp.Expr:
        return self.normalization * sp.integrate(
            expression, (self.variable, self.lower, self.upper))

    def declaration(self) -> sp.Expr:
        integrand = sp.Function("F")(self.variable)
        return self.normalization * sp.Integral(
            integrand, (self.variable, self.lower, self.upper))


@dataclass(frozen=True)
class PlaneWaveAnsatz:
    equations: sp.Tuple
    phase: sp.Expr
    amplitudes: tuple[sp.Symbol, ...]
    wavevector: tuple[sp.Symbol, ...]
    average: PhaseAverage

    def substitutions(self) -> dict[sp.Expr, sp.Expr]:
        return {equation.lhs: equation.rhs for equation in self.equations}


@dataclass(frozen=True)
class BulkModeData:
    ansatz: sp.Equality
    amplitude: sp.Symbol
    amplitude_reality: sp.Expr
    sound_speed_domain: sp.Expr
    dispersion: sp.Equality
    normal_wavevector: sp.Symbol
    interface_equations: sp.Tuple


@dataclass(frozen=True)
class MatrixRoute:
    name: sp.Symbol
    matrix: sp.Matrix


@dataclass(frozen=True)
class RegistryState:
    value: Any | None
    error: str | None

    def unavailable_payload(self) -> sp.Tuple:
        detail: Any = self.error if self.error is not None else sp.Symbol("NO_ERROR_REPORTED")
        return sp.Tuple(sp.Symbol("REGISTRY_UNAVAILABLE"), canonical_solver_object(detail))


@dataclass
class LocusResult:
    equations: tuple[sp.Equality, ...]
    variables: tuple[sp.Symbol, ...]
    raw_branches: list[dict[sp.Symbol, sp.Expr]]
    solution_payload: sp.Tuple
    identically_satisfied: sp.logic.boolalg.Boolean
    inconsistent: sp.logic.boolalg.Boolean
    real_admissible: sp.Tuple


def canonical_rref_rows(vectors: Sequence[sp.MatrixBase], width: int) -> sp.Matrix:
    if not vectors:
        return sp.zeros(0, width)
    rows = sp.Matrix.hstack(*(sp.Matrix(vector) for vector in vectors)).T
    reduced = rows.rref()[0]
    nonzero = [list(reduced.row(index)) for index in range(reduced.rows)
               if any(entry != 0 for entry in reduced.row(index))]
    return sp.Matrix(nonzero) if nonzero else sp.zeros(0, width)


def build_q9(dimension: int) -> Q9Result:
    """Compute invariant subspaces without inserting their known closed forms."""
    n_entries = dimension * dimension
    g_symbols = sp.symbols(f"g1:{n_entries + 1}", real=True)
    G = sp.Matrix(dimension, dimension, g_symbols)
    pairs = tuple((left, right) for left in range(n_entries)
                  for right in range(left, n_entries))
    monomials = tuple(g_symbols[left] * g_symbols[right] for left, right in pairs)
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    width = len(monomials)

    generator_operators: list[sp.SparseMatrix] = []
    for first in range(dimension):
        for second in range(first + 1, dimension):
            generator = sp.zeros(dimension)
            generator[first, second] = 1
            generator[second, first] = -1
            delta = generator * G - G * generator
            delta_terms: list[dict[int, sp.Expr]] = []
            for flat in range(n_entries):
                row, column = divmod(flat, dimension)
                expression = sp.expand(delta[row, column])
                delta_terms.append({
                    index: expression.coeff(symbol)
                    for index, symbol in enumerate(g_symbols)
                    if expression.coeff(symbol) != 0
                })
            entries: dict[tuple[int, int], sp.Expr] = {}
            for source, (left, right) in enumerate(pairs):
                for replacement, value in delta_terms[left].items():
                    destination = pair_index[tuple(sorted((replacement, right)))]
                    entries[(destination, source)] = entries.get((destination, source), 0) + value
                for replacement, value in delta_terms[right].items():
                    destination = pair_index[tuple(sorted((left, replacement)))]
                    entries[(destination, source)] = entries.get((destination, source), 0) + value
            generator_operators.append(sp.SparseMatrix(width, width, entries))

    so_constraints = sp.SparseMatrix.vstack(*generator_operators)
    v1_raw = so_constraints.nullspace()
    v1_matrix = canonical_rref_rows(v1_raw, width)

    signs = tuple([-1] + [1] * (dimension - 1))
    reflection_diagonal = []
    for left, right in pairs:
        li, lj = divmod(left, dimension)
        ri, rj = divmod(right, dimension)
        reflection_diagonal.append(signs[li] * signs[lj] * signs[ri] * signs[rj])
    reflection_operator_ambient = sp.diag(*reflection_diagonal)
    reflection_constraints = reflection_operator_ambient - sp.eye(width)
    o_constraints = sp.SparseMatrix.vstack(so_constraints, sp.SparseMatrix(reflection_constraints))
    v2_raw = o_constraints.nullspace()
    v2_matrix = canonical_rref_rows(v2_raw, width)

    v1_forms = tuple(sp.expand(sum(v1_matrix[row, column] * monomials[column]
                                   for column in range(width)))
                     for row in range(v1_matrix.rows))
    reflected_forms: list[sp.Expr] = []
    reflection_coordinates: list[sp.Matrix] = []
    for row, form in enumerate(v1_forms):
        reflected_row = sp.Matrix([reflection_diagonal[column] * v1_matrix[row, column]
                                   for column in range(width)]).T
        reflected_form = sp.expand(sum(reflected_row[column] * monomials[column]
                                       for column in range(width)))
        reflected_forms.append(reflected_form)
        reflection_coordinates.append(v1_matrix.T.gauss_jordan_solve(reflected_row.T)[0])
    v6_operator = (sp.Matrix.hstack(*reflection_coordinates)
                   if reflection_coordinates else sp.zeros(0, 0))
    minus_eigenvectors = (v6_operator + sp.eye(v6_operator.rows)).nullspace()
    v6_ambient_vectors = [v1_matrix.T * vector for vector in minus_eigenvectors]
    v6_matrix = canonical_rref_rows(v6_ambient_vectors, width)
    v6_forms = tuple(sp.expand(sum(v6_matrix[row, column] * monomials[column]
                                   for column in range(width)))
                     for row in range(v6_matrix.rows))
    pd_form = sp.Add(*v6_forms)

    coordinates = sp.symbols(f"x1:{dimension + 1}", real=True)
    fields = tuple(sp.Function(f"u{index + 1}")(*coordinates) for index in range(dimension))
    gradient_substitution = {
        g_symbols[row * dimension + column]: sp.Derivative(fields[column], coordinates[row])
        for row in range(dimension) for column in range(dimension)
    }
    euler_lagrange: list[sp.Tuple] = []
    for form in v1_forms:
        components = []
        for field_column in range(dimension):
            component = 0
            for derivative_row in range(dimension):
                gij = g_symbols[derivative_row * dimension + field_column]
                first_variation = sp.diff(form, gij).xreplace(gradient_substitution)
                component += sp.diff(first_variation, coordinates[derivative_row])
            components.append(sp.expand(component.doit()))
        euler_lagrange.append(sp.Tuple(*components))

    R0 = sp.diag(*signs)
    payloads: tuple[tuple[str, Any], ...] = (
        ("MONOMIAL_ORDERING", sp.Tuple(*monomials)),
        ("V1_BASIS", v1_matrix),
        ("V1_DIM", sp.Integer(v1_matrix.rows)),
        ("V2_BASIS", v2_matrix),
        ("V2_DIM", sp.Integer(v2_matrix.rows)),
        ("V3_DIFFERENCE", sp.Integer(v1_matrix.rows - v2_matrix.rows)),
        ("R0_MATRIX", R0),
        ("R0_ORTHOGONALITY_RESIDUAL", R0.T * R0 - sp.eye(dimension)),
        ("R0_DETERMINANT", R0.det()),
        ("V4_REFLECTED", sp.Tuple(*reflected_forms)),
        ("V4_RESIDUAL", sp.Tuple(*(sp.expand(reflected - original)
                                    for reflected, original in zip(reflected_forms, v1_forms)))),
        ("V4_SUM", sp.Tuple(*(sp.expand(reflected + original)
                               for reflected, original in zip(reflected_forms, v1_forms)))),
        ("V5_EULER_LAGRANGE", sp.Tuple(*euler_lagrange)),
        ("V6_OPERATOR", v6_operator),
        ("V6_BASIS", v6_matrix),
        ("V6_DIM", sp.Integer(v6_matrix.rows)),
        ("V7_RESIDUAL", sp.Integer(v6_matrix.rows - (v1_matrix.rows - v2_matrix.rows))),
    )
    return Q9Result(payloads, tuple(g_symbols), monomials, v6_matrix, pd_form)


def density_forms(gradient: sp.Matrix, dimension: int) -> dict[str, sp.Expr]:
    curl = sp.Rational(1, 2) * sp.Add(*(
        (gradient[row, column] - gradient[column, row]) ** 2
        for row in range(dimension) for column in range(dimension)
    ))
    divergence = sp.Add(*(gradient[index, index] for index in range(dimension))) ** 2
    symmetric_traceless = ((gradient + gradient.T) / 2
                            - sp.trace(gradient) * sp.eye(dimension) / dimension)
    symtl = sp.Add(*(symmetric_traceless[row, column] ** 2
                     for row in range(dimension) for column in range(dimension)))
    return {"curl": curl, "div": divergence, "symtl": symtl}


def package_stiffness_records(package: str, densities: Mapping[str, sp.Expr],
                              pd_density: sp.Expr) -> tuple[tuple[str, sp.Expr, sp.Expr], ...]:
    if package == "MAIN":
        return (("curl", mu_R / 2, densities["curl"]),
                ("div", B_comp / 2, densities["div"]))
    if package == "XFORM_CURLONLY":
        return (("curl", mu_R / 2, densities["curl"]),)
    if package == "XFORM_DIVONLY":
        return (("div", B_comp / 2, densities["div"]),)
    if package == "XFORM_TRACELESS":
        return (("curl", mu_R / 2, densities["curl"]),
                ("symtl", mu_br, densities["symtl"]))
    if package == "XFORM_EXTRA":
        return (("curl", mu_R / 2, densities["curl"]),
                ("div", B_comp / 2, densities["div"]),
                ("pd", beta / 2, pd_density))
    if package == "XCOEF_BSCALE":
        return (("curl", mu_R / 2, densities["curl"]),
                ("div", s_scale * B_comp / 2, densities["div"]))
    if package == "XCOEF_BSIGN":
        return (("curl", mu_R / 2, densities["curl"]),
                ("div", -B_comp / 2, densities["div"]))
    raise ValueError(f"unknown package: {package}")


def build_action(package: str, dimension: int, q9: Q9Result) -> ActionData:
    coordinates = sp.symbols(f"x1:{dimension + 1}", real=True)
    time_coordinate = sp.Symbol("t", real=True)
    field_arguments = (*coordinates, time_coordinate)
    fields = tuple(sp.Function(f"u{index + 1}")(*field_arguments) for index in range(dimension))
    gradient = sp.Matrix(dimension, dimension,
                         lambda row, column: sp.Derivative(fields[column], coordinates[row]))
    densities = density_forms(gradient, dimension)
    pd_substitution = {
        q9.g_symbols[row * dimension + column]: gradient[row, column]
        for row in range(dimension) for column in range(dimension)
    }
    pd_term = sp.expand(q9.pd_form.xreplace(pd_substitution))
    records = package_stiffness_records(package, densities, pd_term)
    stiffness_terms = tuple(factor * density for _, factor, density in records)
    stiffness = sp.Add(*stiffness_terms)
    kinetic = rho_br / 2 * sp.Add(*(sp.Derivative(field, time_coordinate) ** 2
                                    for field in fields))
    configuration = {
        "background_velocity": sp.zeros(dimension, 1),
        "dissipative_terms": sp.Tuple(),
        "nonlinear_terms": sp.Tuple(),
        "wall_width_fields": sp.Tuple(),
    }
    convective = sp.Add(*tuple(configuration["background_velocity"]))
    dissipative = sp.Add(*tuple(configuration["dissipative_terms"]))
    nonlinear = sp.Add(*tuple(configuration["nonlinear_terms"]))
    width_terms = sp.Add(*tuple(configuration["wall_width_fields"]))
    lagrangian = sp.expand(kinetic - stiffness + convective + dissipative + nonlinear + width_terms)
    action_terms = (kinetic,) + tuple(-term for term in stiffness_terms)
    coefficient_factors = (rho_br / 2,) + tuple(-factor for _, factor, _ in records)
    eom: list[sp.Expr] = []
    for field in fields:
        expression = sp.diff(sp.diff(lagrangian, sp.Derivative(field, time_coordinate)), time_coordinate)
        for coordinate in coordinates:
            expression += sp.diff(sp.diff(lagrangian, sp.Derivative(field, coordinate)), coordinate)
        expression -= sp.diff(lagrangian, field)
        eom.append(sp.expand(expression.doit()))
    return ActionData(tuple(coordinates), time_coordinate, fields, gradient, densities, records,
                      stiffness_terms, stiffness, kinetic, lagrangian, action_terms,
                      coefficient_factors, tuple(eom), pd_term, configuration)


def build_assumptions(package: str, dimension: int, k_symbols: Sequence[sp.Symbol],
                      amplitudes: Sequence[sp.Symbol]) -> AssumptionData:
    common = (
        sp.Q.positive(rho_br), sp.Q.positive(mu_R), sp.Q.positive(B_comp),
        sp.Q.positive(sp.Add(*(component ** 2 for component in k_symbols))),
        *(sp.Q.real(component) for component in k_symbols),
        *(sp.Q.real(component) for component in amplitudes),
        sp.Q.positive(sp.Integer(dimension)), sp.Q.integer(sp.Integer(dimension)),
        sp.Q.positive(c_s0),
    )
    package_parts: tuple[sp.Expr, ...] = ()
    dimension_declarations: tuple[DimensionPremise, ...] = ()
    if package == "XFORM_TRACELESS":
        package_parts = (sp.Q.positive(mu_br),)
    elif package == "XFORM_EXTRA":
        package_parts = (sp.Q.real(beta),)
    elif package == "XCOEF_BSCALE":
        package_parts = (sp.Q.positive(s_scale), sp.Q.nonzero(s_scale - 1))
        dimension_declarations = (DimensionPremise(s_scale, sp.Tuple(0, 0, 0)),)
    return AssumptionData(sp.And(*(common + package_parts)), tuple(common), package_parts,
                          dimension_declarations)


def build_plane_wave_ansatz(action: ActionData, amplitudes: Sequence[sp.Symbol],
                            k_symbols: Sequence[sp.Symbol]) -> PlaneWaveAnsatz:
    phase = (sp.Add(*(wave * coordinate for wave, coordinate in
                      zip(k_symbols, action.coordinates)))
             - sp.sqrt(omegaSquared) * action.time)
    equations = sp.Tuple(*(sp.Eq(field, amplitude * sp.cos(phase), evaluate=False)
                           for field, amplitude in zip(action.fields, amplitudes)))
    average_variable = sp.Symbol("theta", real=True)
    average = PhaseAverage(average_variable, sp.Integer(0), 2 * sp.pi,
                           1 / (2 * sp.pi))
    return PlaneWaveAnsatz(equations, phase, tuple(amplitudes), tuple(k_symbols), average)


def build_bulk_mode(action: ActionData, k_symbols: Sequence[sp.Symbol],
                    bulk_amplitude: sp.Symbol) -> BulkModeData:
    normal_coordinate = sp.Symbol("w", real=True)
    normal_wavevector = sp.Symbol("k_w", real=True)
    scalar_field = sp.Function("phi")(
        *action.coordinates, normal_coordinate, action.time)
    phase = (sp.Add(*(wave * coordinate for wave, coordinate in
                      zip(k_symbols, action.coordinates)))
             + normal_wavevector * normal_coordinate
             - sp.sqrt(omegaSquared) * action.time)
    scalar_ansatz = sp.Eq(
        scalar_field, bulk_amplitude * sp.cos(phase), evaluate=False)
    dispersion = sp.Eq(
        omegaSquared,
        c_s0 ** 2 * (sp.Add(*(wave ** 2 for wave in k_symbols))
                     + normal_wavevector ** 2),
        evaluate=False)
    return BulkModeData(
        scalar_ansatz, bulk_amplitude, sp.Q.real(bulk_amplitude),
        sp.Q.positive(c_s0), dispersion, normal_wavevector, sp.Tuple())


def emit_premises(package: str, dimension: int, action: ActionData,
                  assumptions: AssumptionData, ansatz: PlaneWaveAnsatz,
                  bulk_mode: BulkModeData) -> None:
    amplitudes = ansatz.amplitudes
    k_symbols = ansatz.wavevector
    emit(tag(package, dimension, "PREMISE_U_DIMENSION"), U_DIMENSION)
    emit(tag(package, dimension, "PREMISE_IN_PLANE_FIELD"), sp.Tuple(*amplitudes))
    emit(tag(package, dimension, "PREMISE_REAL_ANSATZ"), ansatz.equations)
    emit(tag(package, dimension, "PREMISE_PHASE_AVERAGE"),
         ansatz.average.declaration())
    emit(tag(package, dimension, "PREMISE_CURL_DENSITY_FORM"),
         action.densities["curl"])
    emit(tag(package, dimension, "PREMISE_COMPRESSION_DENSITY_FORM"),
         action.densities["div"])
    emit(tag(package, dimension, "PREMISE_S_DIMENSION"),
         sp.Tuple(*(sp.Tuple(premise.symbol, premise.vector)
                    for premise in assumptions.dimension_declarations)))
    emit(tag(package, dimension, "PREMISE_JOINT_ASSUMPTIONS"), assumptions.joint)
    emit(tag(package, dimension, "PREMISE_RHO_BR_DOMAIN"), assumptions.common_parts[0])
    emit(tag(package, dimension, "PREMISE_MU_R_DOMAIN"), assumptions.common_parts[1])
    emit(tag(package, dimension, "PREMISE_B_COMP_DOMAIN"), assumptions.common_parts[2])
    emit(tag(package, dimension, "PREMISE_K_NORM_DOMAIN"), assumptions.common_parts[3])
    emit(tag(package, dimension, "PREMISE_K_REAL"),
         sp.Tuple(*assumptions.common_parts[4:4 + dimension]))
    emit(tag(package, dimension, "PREMISE_AMPLITUDES_REAL"),
         sp.Tuple(*assumptions.common_parts[4 + dimension:4 + 2 * dimension]))
    emit(tag(package, dimension, "PREMISE_D_DOMAIN"),
         sp.Tuple(*assumptions.common_parts[4 + 2 * dimension:6 + 2 * dimension]))
    emit(tag(package, dimension, "PREMISE_C_S0_DOMAIN"),
         assumptions.common_parts[-1])
    emit(tag(package, dimension, "PREMISE_PACKAGE_DOMAIN"), sp.Tuple(*assumptions.package_parts))
    emit(tag(package, dimension, "PREMISE_BACKGROUND_REST"),
         action.configuration["background_velocity"])
    emit(tag(package, dimension, "PREMISE_NO_DISSIPATION"),
         action.configuration["dissipative_terms"])
    emit(tag(package, dimension, "PREMISE_LINEAR_RESPONSE"),
         action.configuration["nonlinear_terms"])
    emit(tag(package, dimension, "PREMISE_FROZEN_WALL_WIDTH"),
         action.configuration["wall_width_fields"])
    emit(tag(package, dimension, "PREMISE_BULK_SCALAR_MODE"),
         sp.Tuple(bulk_mode.ansatz, bulk_mode.dispersion,
                  bulk_mode.amplitude_reality, bulk_mode.sound_speed_domain))
    emit(tag(package, dimension, "PREMISE_INTERFACE_EQUATIONS"),
         bulk_mode.interface_equations)


def coefficient_ordering(action: ActionData) -> tuple[sp.Symbol, ...]:
    inventory: set[sp.Symbol] = set()
    for factor in action.coefficient_factors:
        inventory.update(factor.free_symbols)
    return tuple(sorted(inventory, key=str))


def substitute_plane_wave(expression: sp.Expr, ansatz: PlaneWaveAnsatz) -> sp.Expr:
    return sp.expand(expression.xreplace(ansatz.substitutions()).doit())


def dynamical_matrices(action: ActionData,
                       ansatz: PlaneWaveAnsatz) -> tuple[sp.Matrix, sp.Matrix, sp.Expr, sp.Expr]:
    plane_eom = tuple(substitute_plane_wave(expression, ansatz)
                      for expression in action.eom)
    trigonometric_factor = sp.cos(ansatz.phase)
    if any(sp.simplify(expression.xreplace({trigonometric_factor: 0})) != 0
           for expression in plane_eom):
        raise RuntimeError("plane-wave EOM lacks the common trigonometric factor")
    stripped_factor = trigonometric_factor
    stripped = [sp.cancel(expression / stripped_factor) for expression in plane_eom]
    matrix_a = sp.Matrix(stripped).jacobian(ansatz.amplitudes)

    plane_lagrangian = substitute_plane_wave(action.lagrangian, ansatz)
    phase_lagrangian = plane_lagrangian.xreplace({ansatz.phase: ansatz.average.variable})
    averaged = ansatz.average.apply(sp.expand(phase_lagrangian))
    matrix_b = sp.hessian(sp.simplify(averaged), ansatz.amplitudes)
    ratio = sp.cancel(matrix_a[0, 0] / matrix_b[0, 0])
    return sp.simplify(matrix_a), sp.simplify(matrix_b), stripped_factor, ratio


def sign_test(expression: sp.Expr, assumptions: sp.Expr) -> sp.Tuple:
    operand = sp.refine(sp.simplify(expression), assumptions)
    zero = sp.ask(sp.Q.zero(operand), assumptions)
    positive = sp.ask(sp.Q.positive(operand), assumptions)
    negative = sp.ask(sp.Q.negative(operand), assumptions)
    if zero is True:
        token = sp.Symbol("ZERO")
    elif positive is True:
        token = sp.Symbol("POSITIVE")
    elif negative is True:
        token = sp.Symbol("NEGATIVE")
    else:
        token = sp.Symbol("UNDECIDED")
    return sp.Tuple(operand, token)


def exact_iszero(expression: sp.Expr) -> bool:
    """Exact zero oracle for Matrix rank/nullspace elimination."""
    return bool(sp.cancel(expression) == 0)


def syntactic_iszero(expression: sp.Expr) -> bool:
    return bool(expression == 0)


def row_primitive_matrix(matrix: sp.Matrix) -> sp.Matrix:
    """Exact nonzero row rescaling used only to precondition nullspace()."""
    rows: list[list[sp.Expr]] = []
    for row_index in range(matrix.rows):
        values = list(matrix.row(row_index))
        denominator = sp.lcm([sp.fraction(sp.cancel(value))[1]
                              for value in values])
        cleared = [sp.factor(sp.cancel(value * denominator))
                   for value in values]
        content = sp.gcd_list(cleared)
        if content != 0:
            cleared = [sp.factor(sp.cancel(value / content))
                       for value in cleared]
        rows.append(cleared)
    return sp.Matrix(rows)


def exact_nullspace_basis(matrix: sp.Matrix) -> list[sp.Matrix]:
    """Row-only exact preprocessing followed by Matrix.nullspace()."""
    primitive = row_primitive_matrix(matrix)
    reduced_domain, _pivots = DomainMatrix.from_Matrix(primitive).rref()
    reduced = reduced_domain.to_Matrix()
    return reduced.nullspace(
        simplify=False, iszerofunc=syntactic_iszero)


def equations_identically_satisfied(equations: Sequence[sp.Equality],
                                     variables: Sequence[sp.Symbol]) -> sp.logic.boolalg.Boolean:
    del variables
    residuals = [sp.simplify(equation.lhs - equation.rhs) for equation in equations]
    return sp.true if all(residual == 0 for residual in residuals) else sp.false


def equations_inconsistent(equations: Sequence[sp.Equality], variables: Sequence[sp.Symbol],
                           branches: Sequence[Mapping[sp.Symbol, sp.Expr]]
                           ) -> sp.logic.boolalg.Boolean:
    residuals = [sp.simplify(equation.lhs - equation.rhs) for equation in equations]
    if not residuals or all(residual == 0 for residual in residuals) or branches:
        return sp.false
    if any(not residual.free_symbols and residual != 0 for residual in residuals):
        return sp.true
    external = set().union(*(residual.free_symbols for residual in residuals), set()) - set(variables)
    if external or not variables:
        return sp.false
    try:
        numerators = [sp.together(residual).as_numer_denom()[0]
                      for residual in residuals if residual != 0]
        basis = sp.groebner(numerators, *variables)
        return sp.true if any(poly.as_expr() == 1 for poly in basis.polys) else sp.false
    except (sp.PolynomialError, NotImplementedError, ValueError):
        return sp.false


def branch_real_admissibility(branch: Mapping[sp.Symbol, sp.Expr],
                              assumptions: sp.Expr) -> sp.Tuple:
    branch_relations = sp.Tuple(*(sp.Eq(variable, value, evaluate=False)
                                  for variable, value in sorted(branch.items(),
                                                                key=lambda pair: str(pair[0]))))
    try:
        substituted = sp.refine(
            assumptions.subs(branch, simultaneous=True), assumptions)
        predicate = sp.And(
            substituted, *(sp.Eq(variable, value, evaluate=False)
                           for variable, value in branch.items()))
        try:
            result = sp.satisfiable(predicate, use_lra_theory=True)
        except Exception:
            result = sp.satisfiable(predicate)
        return sp.Tuple(branch_relations, predicate,
                        canonical_solver_object(result))
    except Exception:
        return sp.Tuple(branch_relations, assumptions, sp.Symbol("UNDECIDED"))


def locus_protocol(equations: Sequence[sp.Equality], variables: Sequence[sp.Symbol],
                   assumptions: sp.Expr) -> LocusResult:
    equation_tuple = tuple(equations)
    variable_tuple = tuple(variables)
    residuals = [sp.simplify(equation.lhs - equation.rhs) for equation in equation_tuple]
    solve_residuals = [residual for residual in residuals if residual != 0]
    unrestricted = tuple(sp.Symbol(str(variable)) for variable in variable_tuple)
    to_unrestricted = dict(zip(variable_tuple, unrestricted))
    from_unrestricted = dict(zip(unrestricted, variable_tuple))
    unrestricted_residuals = [residual.xreplace(to_unrestricted)
                              for residual in solve_residuals]
    if not unrestricted_residuals:
        raw: list[dict[sp.Symbol, sp.Expr]] = []
    else:
        returned = sp.solve(unrestricted_residuals, unrestricted, dict=True,
                            simplify=True, check=False)
        unrestricted_raw = [returned] if isinstance(returned, dict) else list(returned)
        raw = [{from_unrestricted.get(key, key): value.xreplace(from_unrestricted)
                for key, value in branch.items()}
               for branch in unrestricted_raw]
    solution_payload = sp.Tuple(sp.Tuple(*variable_tuple), canonical_solver_object(raw))
    identity = equations_identically_satisfied(equation_tuple, variable_tuple)
    inconsistent = equations_inconsistent(equation_tuple, variable_tuple, raw)
    admissible = sp.Tuple(*(branch_real_admissibility(branch, assumptions) for branch in raw))
    return LocusResult(equation_tuple, variable_tuple, raw, solution_payload,
                       identity, inconsistent, admissible)


def emit_locus(package: str, dimension: int, base: str, result: LocusResult,
               *, root: int | None = None, stratum: int | None = None) -> None:
    emit(tag(package, dimension, f"{base}_EQUATIONS", root=root, stratum=stratum),
         sp.Tuple(*result.equations))
    emit(tag(package, dimension, f"{base}_SOLUTION", root=root, stratum=stratum),
         result.solution_payload)
    emit(tag(package, dimension, f"{base}_IDENTICALLY_SATISFIED", root=root, stratum=stratum),
         result.identically_satisfied)
    emit(tag(package, dimension, f"{base}_INCONSISTENT", root=root, stratum=stratum),
         result.inconsistent)
    emit(tag(package, dimension, f"{base}_REAL_ADMISSIBLE", root=root, stratum=stratum),
         result.real_admissible)


def roots_with_multiplicity(determinant: sp.Expr) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]]:
    polynomial = sp.Poly(determinant, omegaSquared)
    root_dictionary = sp.roots(polynomial, omegaSquared)
    if sum(int(multiplicity) for multiplicity in root_dictionary.values()) != polynomial.degree():
        raise RuntimeError("SymPy did not return a multiplicity-complete root set")
    ordered_items = sorted(root_dictionary.items(), key=lambda pair: sp.default_sort_key(pair[0]))
    all_roots = tuple(root for root, multiplicity in ordered_items
                      for _ in range(int(multiplicity)))
    distinct_roots = tuple(root for root, _ in ordered_items)
    return all_roots, distinct_roots


def deduplicate_roots(roots: Sequence[sp.Expr], assumptions: sp.Expr) -> tuple[sp.Expr, ...]:
    result: list[sp.Expr] = []
    for candidate in roots:
        if not any(sp.ask(sp.Q.zero(sp.refine(candidate - existing, assumptions)),
                          assumptions) is True for existing in result):
            result.append(candidate)
    return tuple(result)


def coincidence_equations(roots: Sequence[sp.Expr]) -> tuple[sp.Equality, ...]:
    return tuple(sp.Eq(sp.simplify(left - right), 0, evaluate=False)
                 for left, right in itertools.combinations(roots, 2))


@dataclass(frozen=True)
class SpectrumData:
    determinant: sp.Expr
    all_roots: tuple[sp.Expr, ...]
    roots: tuple[sp.Expr, ...]
    root_matrices: tuple[sp.Matrix, ...]
    ranks: tuple[int, ...]
    nullities: tuple[int, ...]


def compute_emit_spectrum_and_modes(package: str, dimension: int, matrix: sp.Matrix,
                                    coefficient_symbols: Sequence[sp.Symbol],
                                    k_symbols: Sequence[sp.Symbol], assumptions: sp.Expr,
                                    *, stratum: int | None = None) -> SpectrumData:
    determinant = sp.factor(matrix.det(method="berkowitz"))
    emit(tag(package, dimension, "DET_M", stratum=stratum), determinant)
    all_roots, solver_distinct = roots_with_multiplicity(determinant)
    roots = deduplicate_roots(solver_distinct, assumptions)
    emit(tag(package, dimension, "ROOT_SOLUTION_SET", stratum=stratum), sp.Tuple(*all_roots))
    emit(tag(package, dimension, "ROOT_DISTINCT", stratum=stratum), sp.Tuple(*roots))
    emit(tag(package, dimension, "ROOT_COUNT_ALL", stratum=stratum), sp.Integer(len(all_roots)))
    emit(tag(package, dimension, "ROOT_COUNT_DISTINCT", stratum=stratum), sp.Integer(len(roots)))
    emit(tag(package, dimension, "ROOT_ORDERING", stratum=stratum), sp.Tuple(*roots))
    for index, root_value in enumerate(roots, 1):
        emit(tag(package, dimension, "VALUE", root=index, stratum=stratum), root_value)
        emit(tag(package, dimension, "SIGN", root=index, stratum=stratum),
             sign_test(root_value, assumptions))

    coincidence = coincidence_equations(roots)
    emit_locus(package, dimension, "ROOT_COINCIDENCE_K",
               locus_protocol(coincidence, k_symbols, assumptions), stratum=stratum)
    emit_locus(package, dimension, "ROOT_COINCIDENCE_COEFF",
               locus_protocol(coincidence, coefficient_symbols, assumptions), stratum=stratum)

    root_matrices: list[sp.Matrix] = []
    ranks: list[int] = []
    nullities: list[int] = []
    k_column = sp.Matrix(k_symbols)
    k_norm_squared = sp.Add(*(component ** 2 for component in k_symbols))
    for index, root_value in enumerate(roots, 1):
        root_matrix = matrix.subs(omegaSquared, root_value).applyfunc(
            lambda entry: sp.factor(sp.cancel(entry)))
        emit(tag(package, dimension, "N1", root=index, stratum=stratum), root_matrix)
        rank_input = row_primitive_matrix(root_matrix)
        rank = int(rank_input.rank(iszerofunc=exact_iszero, simplify=False))
        emit(tag(package, dimension, "N2_RANK", root=index, stratum=stratum), sp.Integer(rank))
        nullity = dimension - rank
        emit(tag(package, dimension, "N2_NULLITY", root=index, stratum=stratum), sp.Integer(nullity))
        stacked = root_matrix.col_join(sp.Matrix([list(k_symbols)]))
        stacked_rank_input = row_primitive_matrix(stacked)
        stacked_rank = int(stacked_rank_input.rank(
            iszerofunc=exact_iszero, simplify=False))
        emit(tag(package, dimension, "N3_STACKED_RANK", root=index, stratum=stratum),
             sp.Integer(stacked_rank))
        transverse_nullity = dimension - stacked_rank
        emit(tag(package, dimension, "N3_TRANSVERSE_NULLITY", root=index, stratum=stratum),
             sp.Integer(transverse_nullity))
        emit(tag(package, dimension, "N4_NULLITY_DIFFERENCE", root=index, stratum=stratum),
             sp.Integer(nullity - transverse_nullity))
        m_dot_k = root_matrix * k_column
        emit(tag(package, dimension, "N5_M_DOT_K", root=index, stratum=stratum),
             m_dot_k.applyfunc(sp.simplify))
        basis = exact_nullspace_basis(root_matrix)
        emit(tag(package, dimension, "N6_BASIS", root=index, stratum=stratum), sp.Tuple(*basis))
        basis_verification = all(
            all(exact_iszero(entry) for entry in root_matrix * vector)
            for vector in basis)
        if not basis_verification:
            raise RuntimeError("row-preconditioned nullspace basis does not annihilate M_r")
        dot_products = tuple((vector.T * k_column)[0] for vector in basis)
        emit(tag(package, dimension, "N6_DOT_K", root=index, stratum=stratum),
             sp.Tuple(*dot_products))
        residuals = tuple(sp.simplify(k_norm_squared * vector - dot * k_column)
                          for vector, dot in zip(basis, dot_products))
        emit(tag(package, dimension, "N6_RESIDUAL", root=index, stratum=stratum),
             sp.Tuple(*residuals))
        emit(tag(package, dimension, "N7_BASIS_COUNT", root=index, stratum=stratum),
             sp.Integer(len(basis)))
        emit(tag(package, dimension, "N7_RESIDUAL", root=index, stratum=stratum),
             sp.Integer(len(basis) - nullity))
        root_matrices.append(root_matrix)
        ranks.append(rank)
        nullities.append(nullity)
    return SpectrumData(determinant, all_roots, roots, tuple(root_matrices),
                        tuple(ranks), tuple(nullities))


def emit_scaling(package: str, dimension: int, roots: Sequence[sp.Expr],
                 k_symbols: Sequence[sp.Symbol], assumptions: sp.Expr) -> None:
    scaling = {component: lambdaScale * component for component in k_symbols}
    for index, root in enumerate(roots, 1):
        scaled = sp.simplify(root.subs(scaling, simultaneous=True))
        emit(tag(package, dimension, "SCALED", root=index), scaled)
        emit(tag(package, dimension, "UNSCALED", root=index), root)
        root_is_zero = sp.ask(sp.Q.zero(sp.refine(root, assumptions)), assumptions) is True
        if root_is_zero:
            ratio: sp.Expr = sp.Symbol("UNDEFINED_RATIO")
            exponent: sp.Expr = sp.Symbol("UNDEFINED_RATIO")
        else:
            ratio = sp.factor(sp.cancel(scaled / root))
            numerator, denominator = sp.fraction(ratio)
            try:
                exponent_candidate = (sp.degree(numerator, lambdaScale)
                                      - sp.degree(denominator, lambdaScale))
                exponent = (sp.Integer(exponent_candidate)
                            if sp.simplify(ratio / lambdaScale ** exponent_candidate) == 1
                            else sp.Symbol("NOT_A_PURE_POWER"))
            except (sp.PolynomialError, TypeError):
                exponent = sp.Symbol("NOT_A_PURE_POWER")
        emit(tag(package, dimension, "SCALE_RATIO", root=index), ratio)
        emit(tag(package, dimension, "SCALE_EXPONENT", root=index), exponent)


def vector_add(left: Sequence[sp.Expr], right: Sequence[sp.Expr]) -> tuple[sp.Expr, ...]:
    return tuple(sp.simplify(a + b) for a, b in zip(left, right))


def vector_scale(vector: Sequence[sp.Expr], scalar: sp.Expr) -> tuple[sp.Expr, ...]:
    return tuple(sp.simplify(scalar * entry) for entry in vector)


class DimensionWalker:
    def __init__(self, action: ActionData, amplitudes: Sequence[sp.Symbol],
                 k_symbols: Sequence[sp.Symbol],
                 coefficient_dimensions: Mapping[sp.Symbol, Sequence[sp.Expr]]):
        self.action = action
        self.symbol_dimensions: dict[sp.Symbol, tuple[sp.Expr, ...]] = {
            **{coordinate: (sp.Integer(1), 0, 0) for coordinate in action.coordinates},
            action.time: (0, 1, 0),
            **{amplitude: tuple(U_DIMENSION) for amplitude in amplitudes},
            **{wave: (-1, 0, 0) for wave in k_symbols},
            omegaSquared: (0, -2, 0),
            lambdaScale: (0, 0, 0),
            c_s0: (1, -1, 0),
            D_symbol: (0, 0, 0),
        }
        for coefficient, vector in coefficient_dimensions.items():
            self.symbol_dimensions[coefficient] = tuple(vector)

    def dimensions(self, expression: sp.Expr) -> list[tuple[sp.Expr, ...]]:
        expression = sp.sympify(expression)
        if expression == 0:
            return []
        if isinstance(expression, sp.Add):
            result: list[tuple[sp.Expr, ...]] = []
            for term in expression.args:
                result.extend(self.dimensions(term))
            return result
        return [self._single(expression)]

    def _single(self, expression: sp.Expr) -> tuple[sp.Expr, ...]:
        if expression.is_Number:
            return (0, 0, 0)
        if expression in self.action.fields:
            return tuple(U_DIMENSION)
        if isinstance(expression, sp.Derivative):
            result = (tuple(U_DIMENSION) if expression.expr in self.action.fields
                      else self._single(expression.expr))
            for variable, count in expression.variable_count:
                derivative_dimension = ((sp.Integer(-1), 0, 0)
                                        if variable in self.action.coordinates else (0, -1, 0))
                result = vector_add(result, vector_scale(derivative_dimension, sp.Integer(count)))
            return result
        if isinstance(expression, sp.Symbol):
            if expression in self.symbol_dimensions:
                return self.symbol_dimensions[expression]
            if expression.name.startswith("g"):
                return (0, 0, 0)
            raise ValueError(f"no dimension assigned to symbol {expression}")
        if isinstance(expression, sp.Add):
            dimensions = self.dimensions(expression)
            if not dimensions:
                return (0, 0, 0)
            reference = dimensions[0]
            if any(any(sp.simplify(left - right) != 0
                       for left, right in zip(vector, reference))
                   for vector in dimensions[1:]):
                raise ValueError(f"inhomogeneous additive subexpression: {expression}")
            return reference
        if isinstance(expression, sp.Mul):
            result: tuple[sp.Expr, ...] = (0, 0, 0)
            for factor in expression.args:
                result = vector_add(result, self._single(factor))
            return result
        if isinstance(expression, sp.Pow):
            if not expression.exp.is_number:
                raise ValueError(f"dimensionful nonnumeric exponent: {expression}")
            return vector_scale(self._single(expression.base), expression.exp)
        if expression.func in (sp.sin, sp.cos):
            argument_dimensions = self.dimensions(expression.args[0])
            if any(any(sp.simplify(entry) != 0 for entry in vector)
                   for vector in argument_dimensions):
                raise ValueError(f"dimensionful trigonometric argument: {expression}")
            return (0, 0, 0)
        if isinstance(expression, sp.Function):
            if expression in self.action.fields:
                return tuple(U_DIMENSION)
            return (0, 0, 0)
        raise ValueError(f"no dimension rule for {type(expression).__name__}: {expression}")


def dimension_solve(action: ActionData, coefficient_symbols: Sequence[sp.Symbol],
                    dimension_declarations: Sequence[DimensionPremise]) -> tuple[
        tuple[sp.Equality, ...], sp.FiniteSet, dict[sp.Symbol, tuple[sp.Expr, ...]],
        int, int, sp.Symbol]:
    declared_dimensions = {premise.symbol: tuple(premise.vector)
                           for premise in dimension_declarations}
    unknown_coefficients = tuple(symbol for symbol in coefficient_symbols
                                 if symbol not in declared_dimensions)
    dimension_unknowns: dict[sp.Symbol, tuple[sp.Symbol, ...]] = {
        coefficient: tuple(sp.Symbol(f"dim_{coefficient}_{slot}")
                           for slot in ("L", "T", "M"))
        for coefficient in unknown_coefficients
    }
    provisional_dimensions: dict[sp.Symbol, tuple[sp.Expr, ...]] = {
        **dimension_unknowns,
        **{coefficient: declared_dimensions[coefficient]
           for coefficient in declared_dimensions if coefficient in coefficient_symbols},
    }
    dummy_amplitudes = tuple(sp.Symbol(f"a{index + 1}", real=True)
                             for index in range(len(action.fields)))
    dummy_k = tuple(sp.Symbol(f"k{index + 1}", real=True)
                    for index in range(len(action.fields)))
    walker = DimensionWalker(action, dummy_amplitudes, dummy_k, provisional_dimensions)
    equations: list[sp.Equality] = []
    incidence_rows: list[list[sp.Expr]] = []
    for term, factor in zip(action.action_terms, action.coefficient_factors):
        term_dimensions = walker.dimensions(term)
        if not term_dimensions:
            continue
        representative = term_dimensions[0]
        equations.extend(sp.Eq(representative[slot], ENERGY_DENSITY_DIM[slot], evaluate=False)
                         for slot in range(3))
        incidence_rows.append([factor.as_powers_dict().get(symbol, 0)
                               for symbol in unknown_coefficients])
    flat_unknowns = tuple(value for coefficient in unknown_coefficients
                          for value in dimension_unknowns[coefficient])
    solution = sp.linsolve(equations, flat_unknowns)
    if solution is sp.EmptySet:
        raise RuntimeError("dimension system is unsolvable")
    solution_tuple = next(iter(solution)) if flat_unknowns else tuple()
    solved_dimensions: dict[sp.Symbol, tuple[sp.Expr, ...]] = {}
    cursor = 0
    for coefficient in unknown_coefficients:
        solved_dimensions[coefficient] = tuple(solution_tuple[cursor:cursor + 3])
        cursor += 3
    for coefficient in coefficient_symbols:
        if coefficient in declared_dimensions:
            solved_dimensions[coefficient] = declared_dimensions[coefficient]
    incidence = (sp.Matrix(incidence_rows)
                 if incidence_rows else sp.zeros(0, len(unknown_coefficients)))
    equation_count = int(incidence.rank())
    unknown_count = len(unknown_coefficients)
    if equation_count < unknown_count:
        determination = sp.Symbol("UNDER_DETERMINED")
    elif equation_count == unknown_count:
        determination = sp.Symbol("EXACTLY_DETERMINED")
    else:
        determination = sp.Symbol("OVER_DETERMINED")
    return (tuple(equations), solution, solved_dimensions, equation_count,
            unknown_count, determination)


def homogeneity_record(expression: Any, walker: DimensionWalker) -> sp.Tuple:
    if isinstance(expression, sp.MatrixBase):
        components = list(expression)
    elif isinstance(expression, (tuple, list, sp.Tuple)):
        components = []
        for item in expression:
            components.extend(list(item) if isinstance(item, sp.MatrixBase) else [sp.sympify(item)])
    else:
        components = [sp.sympify(expression)]
    all_dimensions: list[sp.Tuple] = []
    homogeneous = True
    reference: tuple[sp.Expr, ...] | None = None
    for component in components:
        term_dimensions = walker.dimensions(component)
        all_dimensions.append(sp.Tuple(*(sp.Tuple(*vector) for vector in term_dimensions)))
        if term_dimensions:
            local_reference = term_dimensions[0]
            if any(any(sp.simplify(a - b) != 0 for a, b in zip(vector, local_reference))
                   for vector in term_dimensions[1:]):
                homogeneous = False
            if reference is None:
                reference = local_reference
            elif any(sp.simplify(a - b) != 0 for a, b in zip(local_reference, reference)):
                homogeneous = False
    return sp.Tuple(expression, sp.Tuple(*all_dimensions),
                    sp.true if homogeneous else sp.false)


def emit_dimensions(package: str, dimension: int, action: ActionData,
                    coefficient_symbols: Sequence[sp.Symbol], spectrum: SpectrumData,
                    k_symbols: Sequence[sp.Symbol], amplitudes: Sequence[sp.Symbol],
                    q7_residual: sp.Expr | None, kw_squared_values: Sequence[sp.Expr],
                    assumptions: AssumptionData, registry_state: RegistryState) -> None:
    equations, solution, solved, equation_count, unknown_count, determination = dimension_solve(
        action, coefficient_symbols, assumptions.dimension_declarations)
    ordered_dimensions = sp.Tuple(*(sp.Tuple(coefficient, sp.Tuple(*solved[coefficient]))
                                    for coefficient in coefficient_symbols))
    emit(tag(package, dimension, "DIM_COEFFICIENTS"), ordered_dimensions)
    emit(tag(package, dimension, "DIM_EQUATIONS"), sp.Tuple(*equations))
    emit(tag(package, dimension, "DIM_SOLUTION"), solution)
    emit(tag(package, dimension, "DIM_EQUATION_COUNT"), sp.Integer(equation_count))
    emit(tag(package, dimension, "DIM_UNKNOWN_COUNT"), sp.Integer(unknown_count))
    emit(tag(package, dimension, "DIM_COUNT_DIFFERENCE"),
         sp.Integer(equation_count - unknown_count))
    emit(tag(package, dimension, "DIM_DETERMINACY"), determination)
    walker = DimensionWalker(action, amplitudes, k_symbols, solved)
    k_squared = sp.Add(*(component ** 2 for component in k_symbols))
    for index, root in enumerate(spectrum.roots, 1):
        dimensions = walker.dimensions(sp.cancel(root / k_squared))
        payload = sp.nan if not dimensions else sp.Tuple(*dimensions[0])
        emit(tag(package, dimension, "DIM_OVER_KSQ", root=index), payload)

    emit(tag(package, dimension, "DIM_HOMOGENEITY_ACTION"),
         sp.Tuple(*(homogeneity_record(term, walker) for term in action.action_terms)))
    derived: list[sp.Tuple] = [homogeneity_record(spectrum.determinant, walker)]
    for root, root_matrix, kw_value in zip(spectrum.roots,
                                           spectrum.root_matrices, kw_squared_values):
        derived.append(homogeneity_record(root, walker))
        derived.append(homogeneity_record(root_matrix * sp.Matrix(k_symbols), walker))
        basis = exact_nullspace_basis(root_matrix)
        k_column = sp.Matrix(k_symbols)
        k_norm = (k_column.T * k_column)[0]
        n6 = tuple(sp.simplify(k_norm * vector - (vector.T * k_column)[0] * k_column)
                   for vector in basis)
        derived.append(homogeneity_record(n6, walker))
        derived.append(homogeneity_record(kw_value, walker))
    if q7_residual is not None:
        derived.append(homogeneity_record(q7_residual, walker))
    emit(tag(package, dimension, "DIM_HOMOGENEITY_DERIVED"), sp.Tuple(*derived))

    if registry_state.value is not None:
        by_symbol = {quantity.symbol_name: quantity
                     for quantity in registry_state.value.quantities.values()}
        derived_rows = []
        declared_rows = []
        residual_rows = []
        locus_rows = []
        for coefficient in coefficient_symbols:
            quantity = by_symbol.get(str(coefficient))
            if quantity is None:
                continue
            derived_vector = sp.Tuple(*solved[coefficient])
            declared_vector = sp.Tuple(*(sp.Integer(value) for value in quantity.dimension))
            residual_vector = sp.Tuple(*(sp.simplify(left - right)
                                         for left, right in zip(derived_vector, declared_vector)))
            opaque_locus = quantity.raw["dimension"]["provenance"]["source_locus"]
            derived_rows.append(sp.Tuple(coefficient, derived_vector))
            declared_rows.append(sp.Tuple(coefficient, declared_vector))
            residual_rows.append(sp.Tuple(coefficient, residual_vector))
            locus_rows.append(sp.Tuple(coefficient, canonical_solver_object(opaque_locus)))
        registry_payloads = (
            sp.Tuple(*derived_rows), sp.Tuple(*declared_rows),
            sp.Tuple(*residual_rows), sp.Tuple(*locus_rows))
    else:
        unavailable = registry_state.unavailable_payload()
        registry_payloads = (unavailable,) * 4
    for quantity, payload in zip(
            ("DIM_REGISTRY_DERIVED", "DIM_REGISTRY_DECLARED",
             "DIM_REGISTRY_RESIDUAL", "DIM_REGISTRY_SOURCE_LOCUS"),
            registry_payloads):
        emit(tag(package, dimension, quantity, local=True), payload)


def rank_drop_minors(matrix: sp.Matrix, rank: int) -> tuple[sp.Expr, ...]:
    if rank == 0:
        return tuple()
    return tuple(sp.factor(matrix.extract(rows, columns).det(method="berkowitz"))
                 for rows in itertools.combinations(range(matrix.rows), rank)
                 for columns in itertools.combinations(range(matrix.cols), rank))


def numeric_real(expr: sp.Expr) -> bool:
    simplified = sp.simplify(expr)
    return not simplified.free_symbols and simplified.is_real is True


def witness_for_branch(branch: Mapping[sp.Symbol, sp.Expr],
                       equations: Sequence[sp.Equality], package: str,
                       k_symbols: Sequence[sp.Symbol], amplitudes: Sequence[sp.Symbol],
                       coefficient_symbols: Sequence[sp.Symbol]
                       ) -> dict[sp.Symbol, sp.Expr] | None:
    positive_symbols = {rho_br, mu_R, B_comp, c_s0}
    if package == "XFORM_TRACELESS":
        positive_symbols.add(mu_br)
    if package == "XCOEF_BSCALE":
        positive_symbols.add(s_scale)
    real_symbols = (set(k_symbols) | set(amplitudes) | set(coefficient_symbols)
                    | {beta, mu_br, s_scale, c_s0})
    all_symbols = real_symbols | set().union(*(equation.free_symbols
                                               for equation in equations), set())
    all_symbols |= set().union(*(value.free_symbols for value in branch.values()), set())
    all_symbols |= set(branch)
    for trial in range(1, 32):
        values: dict[sp.Symbol, sp.Expr] = {}
        for offset, symbol in enumerate(sorted(all_symbols - set(branch), key=str)):
            candidate = sp.Integer(sp.prime(trial + offset + 1))
            if symbol == beta:
                candidate = sp.Integer(trial - 1)
            if symbol == s_scale and candidate == 1:
                candidate = sp.Integer(2)
            values[symbol] = candidate
        pending = dict(branch)
        for _ in range(len(pending) + 2):
            progressed = False
            for symbol, expression in tuple(pending.items()):
                value = sp.simplify(expression.subs(values))
                if value.free_symbols:
                    continue
                values[symbol] = value
                del pending[symbol]
                progressed = True
            if not pending or not progressed:
                break
        if pending:
            continue
        if any(not numeric_real(values.get(symbol, symbol)) for symbol in all_symbols):
            continue
        if any(sp.simplify(values[symbol]) <= 0
               for symbol in positive_symbols if symbol in values):
            continue
        if s_scale in values and sp.simplify(values[s_scale] - 1) == 0:
            continue
        if sp.simplify(sum(values[symbol] ** 2 for symbol in k_symbols)) <= 0:
            continue
        if any(sp.simplify((equation.lhs - equation.rhs).subs(values)) != 0
               for equation in equations):
            continue
        return values
    return None


@dataclass(frozen=True)
class Stratum:
    equations: tuple[sp.Equality, ...]
    point: Mapping[sp.Symbol, sp.Expr]
    source_root: int


def compute_rank_drop_and_strata(package: str, dimension: int,
                                 spectrum: SpectrumData,
                                 coefficient_symbols: Sequence[sp.Symbol],
                                 k_symbols: Sequence[sp.Symbol],
                                 amplitudes: Sequence[sp.Symbol],
                                 assumptions: sp.Expr) -> tuple[Stratum, ...]:
    candidates: list[Stratum] = []
    seen: set[str] = set()
    for root_index, (root_matrix, rank) in enumerate(
            zip(spectrum.root_matrices, spectrum.ranks), 1):
        minors = rank_drop_minors(root_matrix, rank)
        emit(tag(package, dimension, "RANK_DROP_MINORS", root=root_index),
             sp.Tuple(*minors))
        equations = tuple(sp.Eq(minor, 0, evaluate=False) for minor in minors)
        protocols = {
            "RANK_DROP_K": locus_protocol(equations, k_symbols, assumptions),
            "RANK_DROP_COEFF": locus_protocol(equations, coefficient_symbols, assumptions),
            "RANK_DROP_JOINT": locus_protocol(
                equations, (*k_symbols, *coefficient_symbols), assumptions),
        }
        for base, protocol in protocols.items():
            emit_locus(package, dimension, base, protocol, root=root_index)
        for branch in protocols["RANK_DROP_JOINT"].raw_branches:
            defining = tuple(sp.Eq(variable, value, evaluate=False)
                             for variable, value in sorted(branch.items(),
                                                           key=lambda pair: str(pair[0])))
            point = witness_for_branch(branch, equations, package, k_symbols,
                                       amplitudes, coefficient_symbols)
            if point is None or int(root_matrix.subs(point).rank()) >= rank:
                continue
            key = sp.srepr(sp.Tuple(*defining))
            if key not in seen:
                seen.add(key)
                candidates.append(Stratum(defining, point, root_index))
    emit(tag(package, dimension, "STRATUM_ORDERING"),
         sp.Tuple(*(sp.Tuple(*stratum.equations) for stratum in candidates)))
    return tuple(candidates)


def emit_strata(package: str, dimension: int, strata: Sequence[Stratum],
                matrix: sp.Matrix, coefficient_symbols: Sequence[sp.Symbol],
                k_symbols: Sequence[sp.Symbol], assumptions: sp.Expr,
                generic_root_jacobian: sp.Matrix) -> None:
    for stratum_index, stratum in enumerate(strata, 1):
        point_payload = sp.Dict({symbol: value for symbol, value in sorted(
            stratum.point.items(), key=lambda pair: str(pair[0]))})
        residual = sp.Tuple(*(sp.simplify((equation.lhs - equation.rhs)
                                          .subs(stratum.point))
                              for equation in stratum.equations))
        emit(tag(package, dimension, "DEFINING_EQUATIONS", stratum=stratum_index),
             sp.Tuple(*stratum.equations))
        emit(tag(package, dimension, "POINT", stratum=stratum_index), point_payload)
        emit(tag(package, dimension, "POINT_RESIDUAL", stratum=stratum_index), residual)
        specialized_matrix = sp.simplify(matrix.subs(stratum.point))
        specialized = compute_emit_spectrum_and_modes(
            package, dimension, specialized_matrix, coefficient_symbols, k_symbols,
            assumptions.subs(stratum.point), stratum=stratum_index)
        restricted = generic_root_jacobian.subs(stratum.point)
        recomputed = sp.Matrix([[sp.diff(root, coefficient)
                                 for coefficient in coefficient_symbols]
                                for root in specialized.roots])
        emit(tag(package, dimension, "ROOT_COEFFICIENT_JACOBIAN_RESTRICTED",
                 stratum=stratum_index), restricted)
        emit(tag(package, dimension, "ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED",
                 stratum=stratum_index), recomputed)


def q7_objects(dimension: int, action: ActionData
               ) -> tuple[sp.Expr | None, tuple[tuple[str, Any], ...]]:
    if dimension != 3:
        return None, tuple()
    g = sp.symbols("g11 g12 g13 g21 g22 g23 g31 g32 g33", real=True)
    G = sp.Matrix(3, 3, g)
    substitution = {action.gradient[row, column]: G[row, column]
                    for row in range(3) for column in range(3)}
    w_full = action.stiffness.xreplace(substitution)
    selected_curl_terms = tuple(term.xreplace(substitution)
                                for (kind, _, _), term in zip(
                                    action.stiffness_records, action.stiffness_terms)
                                if kind == "curl")
    curl_term = sp.Add(*selected_curl_terms)
    curl_density = density_forms(G, 3)["curl"]
    c_vector = sp.Matrix([
        sp.Add(*(sp.LeviCivita(index, row, column) * G[row, column]
                 for row in range(3) for column in range(3)))
        for index in range(3)
    ])
    curl_reference = sp.expand((c_vector.T * c_vector)[0])
    residual = sp.expand(curl_density - curl_reference)
    return residual, (
        ("Q7_W_FULL", w_full),
        ("Q7_CURL_TERM", curl_term),
        ("Q7_CURL_DENSITY", curl_density),
        ("Q7_CURL_REFERENCE", curl_reference),
        ("Q7_RESIDUAL", residual),
    )


def q11_objects(package: str, dimension: int, spectrum: SpectrumData,
                coefficient_symbols: Sequence[sp.Symbol],
                ansatz: PlaneWaveAnsatz, assumptions: sp.Expr,
                bulk_mode: BulkModeData) -> tuple[sp.Expr, ...]:
    kw_values: list[sp.Expr] = []
    for root_index, root in enumerate(spectrum.roots, 1):
        equation = bulk_mode.dispersion.xreplace({omegaSquared: root})
        solution = sp.solve(equation, bulk_mode.normal_wavevector ** 2,
                            dict=False, simplify=True)
        if not solution:
            raise RuntimeError("bulk normal-wavevector equation is unsolved")
        kw_value = sp.simplify(solution[0])
        kw_values.append(kw_value)
        emit(tag(package, dimension, "KW_EQUATION", root=root_index), equation)
        emit(tag(package, dimension, "KW_SQUARED", root=root_index), kw_value)
        emit(tag(package, dimension, "KW_SIGN", root=root_index),
             sign_test(kw_value, assumptions))
        zero_locus = locus_protocol(
            (sp.Eq(kw_value, 0, evaluate=False),),
            (*coefficient_symbols, c_s0), assumptions)
        emit_locus(package, dimension, "KW_ZERO_LOCUS", zero_locus,
                   root=root_index)

    c1 = tuple(bulk_mode.interface_equations)
    unknowns = tuple(ansatz.amplitudes) + (bulk_mode.amplitude,)
    coefficient_matrix = (sp.linear_eq_to_matrix(
        [equation.lhs - equation.rhs for equation in c1], unknowns)[0]
        if c1 else sp.zeros(0, len(unknowns)))
    c2_count = len(unknowns)
    c3_rank = int(coefficient_matrix.rank())
    emit(tag(package, dimension, "C1_EQUATIONS"), sp.Tuple(*c1))
    emit(tag(package, dimension, "C2_UNKNOWNS"), sp.Tuple(*unknowns))
    emit(tag(package, dimension, "C2_COUNT"), sp.Integer(c2_count))
    emit(tag(package, dimension, "C3_RANK"), sp.Integer(c3_rank))
    emit(tag(package, dimension, "C4_DIFFERENCE"),
         sp.Integer(c2_count - c3_rank))
    return tuple(kw_values)


def emit_q9(package: str, dimension: int, q9: Q9Result) -> None:
    for quantity, payload in q9.payloads:
        emit(tag(package, dimension, quantity), payload)


def run_cell(package: str, dimension: int, q9: Q9Result,
             registry_state: RegistryState) -> None:
    emit_q9(package, dimension, q9)
    action = build_action(package, dimension, q9)
    k_symbols = tuple(sp.Symbol(f"k{index + 1}", real=True)
                      for index in range(dimension))
    amplitudes = tuple(sp.Symbol(f"a{index + 1}", real=True)
                       for index in range(dimension))
    bulk_amplitude = sp.Symbol("A", real=True)
    assumptions = build_assumptions(package, dimension, k_symbols, amplitudes)
    ansatz = build_plane_wave_ansatz(action, amplitudes, k_symbols)
    bulk_mode = build_bulk_mode(action, k_symbols, bulk_amplitude)
    emit_premises(package, dimension, action, assumptions, ansatz, bulk_mode)
    emit(tag(package, dimension, "PD_TERM"), action.pd_term)
    emit(tag(package, dimension, "L"), action.lagrangian)
    emit(tag(package, dimension, "STIFFNESS_TERMS"),
         sp.Tuple(*action.stiffness_terms))
    emit(tag(package, dimension, "EULER_LAGRANGE"),
         sp.Tuple(*(sp.Eq(expression, 0, evaluate=False)
                    for expression in action.eom)))
    coefficients = coefficient_ordering(action)
    emit(tag(package, dimension, "COEFFICIENT_ORDERING"),
         sp.Tuple(*coefficients))

    matrix_a, matrix_b, stripped_factor, ratio = dynamical_matrices(action, ansatz)
    emit(tag(package, dimension, "M_ROUTE_A_STRIPPED_FACTOR"), stripped_factor)
    emit(tag(package, dimension, "M_A"), matrix_a)
    emit(tag(package, dimension, "M_B"), matrix_b)
    emit(tag(package, dimension, "M_A_MINUS_M_B"),
         sp.simplify(matrix_a - matrix_b))
    emit(tag(package, dimension, "M_A11_OVER_M_B11"), ratio)
    emit(tag(package, dimension, "M_RESIDUAL_TEST_SCOPE"),
         sp.Tuple(sp.Symbol("SUPPLIED_CLASSIFICATION"),
                  sp.Symbol("CODING_CONSISTENCY")))
    route_a = MatrixRoute(sp.Symbol("M_A"), matrix_a)
    route_b = MatrixRoute(sp.Symbol("M_B"), matrix_b)
    downstream_route = route_b
    matrix = downstream_route.matrix
    emit(tag(package, dimension, "M_ROUTE_USED"), downstream_route.name)
    coefficient_jacobian = sp.Tuple(*(
        matrix.applyfunc(lambda entry: sp.diff(entry, coefficient))
        for coefficient in coefficients))
    emit(tag(package, dimension, "M_COEFFICIENT_JACOBIAN"),
         coefficient_jacobian)

    spectrum = compute_emit_spectrum_and_modes(
        package, dimension, matrix, coefficients, k_symbols, assumptions.joint)
    emit_scaling(package, dimension, spectrum.roots,
                 k_symbols, assumptions.joint)
    root_jacobian = sp.Matrix([
        [sp.diff(root, coefficient) for coefficient in coefficients]
        for root in spectrum.roots
    ])
    emit(tag(package, dimension, "ROOT_COEFFICIENT_JACOBIAN"),
         root_jacobian)

    q7_residual, q7_payloads = q7_objects(dimension, action)
    for quantity, payload in q7_payloads:
        emit(tag(package, dimension, quantity), payload)

    strata = compute_rank_drop_and_strata(
        package, dimension, spectrum, coefficients, k_symbols,
        amplitudes, assumptions.joint)
    emit_strata(package, dimension, strata, matrix, coefficients,
                k_symbols, assumptions.joint, root_jacobian)

    kw_values = q11_objects(
        package, dimension, spectrum, coefficients, ansatz,
        assumptions.joint, bulk_mode)
    emit_dimensions(package, dimension, action, coefficients, spectrum,
                    k_symbols, amplitudes, q7_residual, kw_values,
                    assumptions, registry_state)


def main() -> int:
    start = time.monotonic()
    registry_error: str | None = None
    try:
        registry = load_registry(REDUCTION_DIR)
    except RegistryValidationError as error:
        registry = None
        registry_error = str(error)
    registry_state = RegistryState(registry, registry_error)

    declared = tuple((package, dimension)
                     for package, dimensions in PACKAGE_SWEEPS.items()
                     for dimension in dimensions)
    completed: list[tuple[str, int]] = []
    q9_cache: dict[int, Q9Result] = {}
    for dimension in sorted({dimension for _, dimension in declared}):
        q9_cache[dimension] = build_q9(dimension)
        for package, dimensions in PACKAGE_SWEEPS.items():
            if dimension not in dimensions:
                continue
            run_cell(package, dimension, q9_cache[dimension], registry_state)
            completed.append((package, dimension))

    skipped = tuple(pair for pair in declared if pair not in completed)
    emit("PY_S11_RUN_PAIRS",
         sp.Tuple(*(sp.Tuple(sp.Symbol(package), dimension)
                    for package, dimension in completed)))
    emit("PY_S11_SKIPPED_PAIRS",
         sp.Tuple(*(sp.Tuple(sp.Symbol(package), dimension)
                    for package, dimension in skipped)))
    emit("PY_S11_LOCAL_REGISTRY_VALIDATION_ERROR",
         registry_error if registry_error is not None else sp.Symbol("NO_ERROR"))
    local_listing_name = "PY_S11_LOCAL_TAG_NAMES"
    local_names.append(local_listing_name)
    emit(local_listing_name, sp.Tuple(*(sp.Symbol(name) for name in local_names)))
    _runtime = time.monotonic() - start
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
