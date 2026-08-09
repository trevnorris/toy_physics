#!/usr/bin/env python3
"""S10 engine 2: derive the brane mode audit from the shared action with SymPy."""

from __future__ import annotations

from pathlib import Path


EXPORT_PATH = Path(__file__).with_name("S10_exports.py")  # CONTROL · generated-module publication path
if __name__ == "__main__":
    EXPORT_PATH.unlink(missing_ok=True)


import itertools
import os
import re
import signal
from dataclasses import dataclass
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence

import sympy as sp
from sympy.core.symbol import Str

from extract_knob_inventory import (
    CLASS_TAGS,
    assumption_channel_operands_and_residual,
    declaration_classes,
    file_sha256,
)
from S9_exports import LEDGER as S9_LEDGER


D = S9_LEDGER["D"]["value"]  # STRUCTURAL · imported brane spatial-dimension symbol
rho_br = S9_LEDGER["rho_br"]["value"]  # KNOB · imported brane inertia-density symbol
mu_R = S9_LEDGER["mu_R"]["value"]  # KNOB · imported brane shear-modulus symbol
length_dimension = S9_LEDGER["length_dimension"]["value"]  # STRUCTURAL · imported length unit-basis marker
time_dimension = S9_LEDGER["time_dimension"]["value"]  # STRUCTURAL · imported time unit-basis marker
field_dimension = S9_LEDGER["field_dimension"]["value"]  # PREMISE · imported displacement-field unit premise
energy_density_dimension = S9_LEDGER["dim_energy_density"]["value"]  # PREMISE · imported action-density unit premise
wavevector_norm_dimension = S9_LEDGER["wavevector_norm_dimension"]["value"]  # PREMISE · imported wave-norm unit premise
squared_velocity_dimension = S9_LEDGER["dim_squared_velocity"]["value"]  # PREMISE · imported squared-speed unit premise
s_rho = sp.Symbol("s_rho", real=True)  # CONTROL · anisotropic-inertia ablation coefficient
s = sp.Symbol("s", real=True)  # CONTROL · stiffness-coefficient ablation scale
omegaSquared = S9_LEDGER["omegaSquared"]["value"]  # COORDINATE · imported squared-frequency spectral variable
lambdaScale = S9_LEDGER["lambdaScale"]["value"]  # CONTROL · imported positive wavevector homogeneity scale

ZERO_DIM = tuple(length_dimension - length_dimension)  # DERIVED · neutral unit vector built from the imported basis
LENGTH_DIM = tuple(length_dimension)  # DERIVED · sequence view of the imported length marker
TIME_DIM = tuple(time_dimension)  # DERIVED · sequence view of the imported time marker
FIELD_DIM = tuple(field_dimension)  # DERIVED · sequence view of the imported field-unit premise
ENERGY_DENSITY_DIM = tuple(energy_density_dimension)  # DERIVED · sequence view of the imported action-density premise
WAVENUMBER_DIM = tuple(wavevector_norm_dimension / 2)  # DERIVED · component wave-number units built from the imported norm premise
OMEGA_SQUARED_DIM = tuple(squared_velocity_dimension + wavevector_norm_dimension)  # DERIVED · spectral units built from imported unit objects
OPERATION_TIMEOUT_SECONDS = 60  # CONTROL · symbolic-operation timeout
declaration_class_map = declaration_classes(Path(__file__))  # CONTROL · source-derived declaration class lookup
declared_symbol_classes = {}  # CONTROL · source-derived classes for symbols reaching the export boundary


def directly_declared_symbols(value: object) -> tuple[sp.Symbol, ...]:
    if isinstance(value, sp.Symbol):
        return (value,)
    if isinstance(value, (list, tuple)):
        nested = [directly_declared_symbols(item) for item in value]
        if all(items for items in nested):
            return tuple(symbol for items in nested for symbol in items)
    return ()


def register_declared_symbols(namespace: Mapping[str, object]) -> None:
    for variable_name, class_tag in declaration_class_map.items():
        if variable_name not in namespace:
            continue
        for symbol in directly_declared_symbols(namespace[variable_name]):
            declared_symbol_classes.setdefault(symbol.name, class_tag)


@dataclass(frozen=True)
class Package:
    name: str
    dimensions: tuple[int, ...]
    stiffness: str
    stiffness_sign: int = -1
    anisotropic_kinetic: bool = False
    coefficient_scale: bool = False


PACKAGES = (  # CONTROL · action-package and construction-sweep specifications
    Package("MAIN", (2, 3, 4, 5), "curl"),
    Package("XFORM_FULLGRAD", (3, 4), "fullgrad"),
    Package("XFORM_DIVONLY", (3, 4), "divonly"),
    Package("XFORM_SIGNFLIP", (3, 4), "curl", stiffness_sign=1),
    Package("XFORM_ANISO", (3, 4), "curl", anisotropic_kinetic=True),
    Package("XCOEF_SCALE", (3,), "curl", coefficient_scale=True),
)

INDEXED_EXPORT_TOKEN = re.compile(r"(?P<label>ROOT|STRATUM|COEFFICIENT)(?P<index>[0-9]+)")  # CONTROL · indexed-emission grouping grammar
SYMBOLIC_D_EXPORT_SUFFIXES = (  # CONTROL · ordered symbolic dimension-path export names
    "PREMISE_U_DIMENSION",
    "Q6_ENERGY_DENSITY_DIMENSION",
    "Q6_DIMENSION_SOLUTION",
)
POSITED_OUTPUT_CLASSES = {  # CONTROL · posited output class assignments
    "PREMISE_FIELD_SECTOR": "PREMISE",
    "PREMISE_U_DIMENSION": "PREMISE",
    "PREMISE_ANSATZ": "PREMISE",
    "PREMISE_PERIOD_AVERAGE": "PREMISE",
    "PREMISE_BACKGROUND_VELOCITY": "PREMISE",
    "PREMISE_TIME_ODD_KERNEL": "PREMISE",
    "PREMISE_RESPONSE_DEGREE": "PREMISE",
    "PREMISE_STIFFNESS_INPUT": "PREMISE",
    "PREMISE_RHO_DOMAIN": "PREMISE",
    "PREMISE_MU_DOMAIN": "PREMISE",
    "PREMISE_WAVEVECTOR_NORM_DOMAIN": "PREMISE",
    "PREMISE_WAVEVECTOR_REAL_DOMAIN": "PREMISE",
    "PREMISE_AMPLITUDE_REAL_DOMAIN": "PREMISE",
    "PREMISE_BRANE_DIMENSION_DOMAIN": "PREMISE",
    "ASSUMPTION_JOINT_PREDICATE": "PREMISE",
    "Q6_CONTROL_DIMENSION_PREMISES": "PREMISE",
    "Q6_ENERGY_DENSITY_DIMENSION": "PREMISE",
    "Q2_RESIDUAL_TEST_SCOPE": "CONTROL",
    "Q2_DOWNSTREAM_ROUTE": "CONTROL",
}


class OperationTimeout(RuntimeError):
    pass


def _timeout_handler(_signum: int, _frame: object) -> None:
    raise OperationTimeout


class Emitter:
    def __init__(self) -> None:
        self.names: set[str] = set()
        self.local_names: list[str] = []
        self.values: dict[str, object] = {}
        self.dimensions: dict[str, object] = {}

    def emit(self, name: str, payload: object, *, dimension: object | None = None) -> None:
        if name in self.names:
            raise RuntimeError(f"duplicate emitted tag: {name}")
        self.names.add(name)
        self.values[name] = payload
        if dimension is not None:
            self.dimensions[name] = dimension
        if name.startswith("PY_S10_LOCAL_"):
            self.local_names.append(name)
        rendered = str(payload).replace("\r", " ").replace("\n", " ")
        print(f"{name}: {rendered}")


def assumed_simplify(expr: object, assumptions: sp.logic.boolalg.Boolean) -> object:
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(lambda item: sp.refine(sp.simplify(item), assumptions))
    if isinstance(expr, sp.Basic):
        return sp.refine(sp.simplify(expr), assumptions)
    return expr


def algebraic_normalize(expr: sp.Expr, assumptions: sp.logic.boolalg.Boolean) -> sp.Expr:
    return sp.refine(sp.factor(sp.cancel(sp.together(expr))), assumptions)


def polynomial_iszero(expr: sp.Expr) -> bool:
    numerator = sp.cancel(sp.together(expr)).as_numer_denom()[0]
    return numerator == 0


def exact_row_equivalent(matrix: sp.Matrix) -> sp.Matrix:
    domain_matrix = sp.polys.matrices.DomainMatrix.from_Matrix(matrix)
    reduced, _denominator, _pivots = domain_matrix.rref_den()
    return reduced.to_Matrix()


def tuple_expr(values: Iterable[object]) -> sp.Tuple:
    return sp.Tuple(*(value if isinstance(value, sp.Basic) else sp.sympify(value) for value in values))


def dim_expr(value: Sequence[sp.Expr]) -> sp.Tuple:
    return sp.Tuple(*value)


def add_dims(left: Sequence[sp.Expr], right: Sequence[sp.Expr]) -> tuple[sp.Expr, ...]:
    return tuple(sp.Add(a, b) for a, b in zip(left, right))


def scale_dim(value: Sequence[sp.Expr], exponent: sp.Expr) -> tuple[sp.Expr, ...]:
    return tuple(sp.Mul(exponent, component) for component in value)


class DimensionWalker:
    """Dimension expressions by walking all factors, including bare fields."""

    def __init__(
        self,
        symbol_dimensions: Mapping[sp.Symbol, Sequence[sp.Expr]],
        fields: Sequence[sp.Expr],
        coordinate_dimensions: Mapping[sp.Symbol, Sequence[sp.Expr]],
        assumptions: sp.logic.boolalg.Boolean,
    ) -> None:
        self.symbol_dimensions = dict(symbol_dimensions)
        self.fields = set(fields)
        self.coordinate_dimensions = dict(coordinate_dimensions)
        self.assumptions = assumptions

    @staticmethod
    def indeterminate_dimension(expr: sp.Expr) -> tuple[sp.Expr, ...]:
        failed_head = Str(getattr(expr.func, "__name__", type(expr).__name__))  # DERIVED · unsupported expression-head marker
        marker = sp.Function("indeterminate_dimension")(failed_head)  # DERIVED · unresolved unit-walk record
        return (marker, marker, marker)

    def dimension(self, expr: object) -> tuple[sp.Expr, ...]:
        if isinstance(expr, Str):
            return ZERO_DIM
        if expr is sp.zoo or expr is sp.nan or expr in (sp.oo, -sp.oo):
            return ZERO_DIM
        if expr == 0 or expr.is_Number:
            return ZERO_DIM
        if expr in self.fields:
            return FIELD_DIM
        if isinstance(expr, sp.Derivative):
            result = self.dimension(expr.expr)
            for variable, count in expr.variable_count:
                coordinate_dim = self.coordinate_dimensions[variable]
                result = add_dims(result, scale_dim(coordinate_dim, -sp.Integer(count)))
            return result
        if isinstance(expr, sp.Symbol):
            return tuple(self.symbol_dimensions.get(expr, ZERO_DIM))
        if isinstance(expr, sp.Pow):
            if expr.exp.free_symbols:
                raise RuntimeError(f"dimensionful expression has symbolic exponent: {expr}")
            return scale_dim(self.dimension(expr.base), expr.exp)
        if isinstance(expr, sp.Mul):
            result = ZERO_DIM
            for factor in expr.args:
                result = add_dims(result, self.dimension(factor))
            return result
        if isinstance(expr, sp.Add):
            dimensions = [self.dimension(term) for term in expr.args if term != 0]
            if not dimensions:
                return ZERO_DIM
            return dimensions[0]
        if expr.func is sp.Abs:
            return self.dimension(expr.args[0])
        if expr.func in (sp.sin, sp.cos, sp.tan, sp.exp, sp.log, sp.sign):
            return ZERO_DIM
        if expr.is_Function:
            return self.indeterminate_dimension(expr)
        return self.indeterminate_dimension(expr)

    def component_term_dimensions(self, expr: object) -> tuple[tuple[sp.Expr, ...], ...]:
        if isinstance(expr, Str):
            return (ZERO_DIM,)
        expanded = sp.expand(expr)
        if expanded == 0:
            return ()
        return tuple(self.dimension(term) for term in sp.Add.make_args(expanded))

    def nested_add_dimensions(self, expr: object) -> tuple[sp.Tuple, ...]:
        reports: list[sp.Tuple] = []
        seen: set[sp.Add] = set()
        for node in sorted(expr.atoms(sp.Add), key=sp.default_sort_key):
            if node in seen:
                continue
            seen.add(node)
            term_dimensions = tuple(self.dimension(term) for term in node.args if term != 0)
            reports.append(
                sp.Tuple(node, sp.Tuple(*(dim_expr(value) for value in term_dimensions)))
            )
        return tuple(reports)

    def report(self, obj: object) -> tuple[sp.Tuple, sp.Tuple, sp.Tuple, sp.logic.boolalg.Boolean]:
        if isinstance(obj, sp.MatrixBase):
            components = list(obj)
        elif isinstance(obj, (list, tuple, sp.Tuple)):
            components = []
            for item in obj:
                if isinstance(item, sp.MatrixBase):
                    components.extend(list(item))
                elif isinstance(item, sp.Expr):
                    components.append(item)
        elif isinstance(obj, (sp.Expr, Str)):
            components = [obj]
        else:
            components = []

        all_term_dimensions: list[sp.Tuple] = []
        all_component_dimensions: list[sp.Tuple] = []
        all_nested_add_dimensions: list[sp.Tuple] = []
        homogeneous = True
        for component in components:
            term_dimensions = self.component_term_dimensions(component)
            rendered_terms = sp.Tuple(*(dim_expr(value) for value in term_dimensions))
            all_term_dimensions.append(rendered_terms)
            unique: list[tuple[sp.Expr, ...]] = []
            for value in term_dimensions:
                simplified = tuple(
                    sp.refine(sp.simplify(entry), self.assumptions) for entry in value
                )
                if simplified not in unique:
                    unique.append(simplified)
            all_component_dimensions.append(sp.Tuple(*(dim_expr(value) for value in unique)))
            if len(unique) > 1:
                homogeneous = False
            nested_reports = self.nested_add_dimensions(component)
            all_nested_add_dimensions.append(sp.Tuple(*nested_reports))
            for report in nested_reports:
                node_dimensions = report[1]
                if len(set(node_dimensions)) > 1:
                    homogeneous = False
        return (
            sp.Tuple(*all_term_dimensions),
            sp.Tuple(*all_component_dimensions),
            sp.Tuple(*all_nested_add_dimensions),
            sp.true if homogeneous else sp.false,
        )


def emit_physical(
    emitter: Emitter,
    walker: DimensionWalker,
    name: str,
    payload: object,
    homogeneity_class: str = "UNSOLVED_EXPRESSION",
    *,
    dimension_payload: object | None = None,
) -> None:
    dimension_source = payload if dimension_payload is None else dimension_payload
    term_dimensions, dimensions, nested_add_dimensions, homogeneous = walker.report(
        dimension_source
    )
    emitter.emit(name, payload, dimension=dimensions)
    emitter.emit(name + "_Q6_TERM_DIMENSIONS", term_dimensions)
    emitter.emit(name + "_Q6_DIMENSIONS", dimensions)
    emitter.emit(name + "_Q6_NESTED_ADD_DIMENSIONS", nested_add_dimensions)
    emitter.emit(name + f"_Q6_{homogeneity_class}_HOMOGENEITY", homogeneous)


def stiffness_density(kind: str, gradient: Sequence[Sequence[sp.Expr]]) -> sp.Expr:
    n = len(gradient)
    if kind == "curl":
        return sp.Rational(1, 2) * sp.Add(
            *((gradient[i][j] - gradient[j][i]) ** 2 for i in range(n) for j in range(n))
        )
    if kind == "fullgrad":
        return sp.Add(*(gradient[i][j] ** 2 for i in range(n) for j in range(n)))
    if kind == "divonly":
        return sp.Add(*(gradient[i][i] for i in range(n))) ** 2
    raise RuntimeError(f"unknown stiffness package: {kind}")


def build_joint_assumptions(
    package: Package,
    kvec: Sequence[sp.Symbol],
    avec: Sequence[sp.Symbol],
) -> sp.logic.boolalg.Boolean:
    k_squared = sp.Add(*(component**2 for component in kvec))
    predicates: list[sp.logic.boolalg.Boolean] = [
        sp.Q.positive(rho_br),
        sp.Q.positive(mu_R),
        sp.Q.positive(k_squared),
        *(sp.Q.real(component) for component in kvec),
        *(sp.Q.real(component) for component in avec),
        sp.Q.integer(D),
        sp.Q.positive(D),
        sp.Q.positive(lambdaScale),
    ]
    if package.anisotropic_kinetic:
        predicates.extend((sp.Q.positive(s_rho), sp.Ne(s_rho, 1)))
    if package.coefficient_scale:
        predicates.append(sp.Q.positive(s))
    return sp.And(*predicates)


def build_action(
    package: Package,
    n: int,
    t: sp.Symbol,
    xvec: Sequence[sp.Symbol],
) -> tuple[list[sp.Expr], sp.Expr, list[sp.Expr], list[sp.Expr], sp.Expr]:
    fields = [sp.Function(f"u{index + 1}")(t, *xvec) for index in range(n)]  # PREMISE · brane displacement field heads
    time_derivatives = [sp.diff(field, t) for field in fields]
    gradient = [[sp.diff(fields[j], xvec[i]) for j in range(n)] for i in range(n)]
    inertial_coefficients = [rho_br for _ in range(n)]
    if package.anisotropic_kinetic:
        inertial_coefficients[0] = s_rho * rho_br
    stiffness_coefficient = s * mu_R if package.coefficient_scale else mu_R
    kinetic = sp.Add(
        *(
            sp.Rational(1, 2) * inertial_coefficients[j] * time_derivatives[j] ** 2
            for j in range(n)
        )
    )
    stiffness = stiffness_density(package.stiffness, gradient)
    lagrangian = kinetic + sp.Integer(package.stiffness_sign) * sp.Rational(1, 2) * stiffness_coefficient * stiffness
    return fields, sp.expand(lagrangian), inertial_coefficients, [stiffness_coefficient], stiffness


def euler_lagrange_system(
    lagrangian: sp.Expr,
    fields: Sequence[sp.Expr],
    t: sp.Symbol,
    xvec: Sequence[sp.Symbol],
) -> sp.Matrix:
    equations = []
    for field in fields:
        equation = sp.diff(sp.diff(lagrangian, sp.diff(field, t)), t)
        equation += sp.Add(
            *(sp.diff(sp.diff(lagrangian, sp.diff(field, coordinate)), coordinate) for coordinate in xvec)
        )
        equation -= sp.diff(lagrangian, field)
        equations.append(sp.expand(equation))
    return sp.Matrix(equations)


def substitute_ansatz(
    expression: object,
    fields: Sequence[sp.Expr],
    ansatz: Sequence[sp.Expr],
) -> object:
    substitutions = dict(zip(fields, ansatz))
    if isinstance(expression, sp.MatrixBase):
        return expression.applyfunc(lambda item: item.subs(substitutions).doit())
    if isinstance(expression, sp.Expr):
        return expression.subs(substitutions).doit()
    raise RuntimeError("ansatz substitution received an unsupported object")


def derive_dimension_solution(
    package: Package,
    lagrangian: sp.Expr,
    fields: Sequence[sp.Expr],
    t: sp.Symbol,
    xvec: Sequence[sp.Symbol],
    assumptions: sp.logic.boolalg.Boolean,
) -> tuple[
    list[sp.Equality],
    tuple[sp.Symbol, ...],
    sp.Tuple,
    dict[sp.Symbol, tuple[sp.Expr, ...]],
    int,
    int,
    sp.Symbol,
    sp.logic.boolalg.Boolean,
]:
    dimensionful_coefficient_symbols = tuple(
        sorted(
            lagrangian.free_symbols.difference({t, *xvec, D, s_rho, s}),
            key=sp.default_sort_key,
        )
    )
    coefficient_unknowns = {
        coefficient: sp.symbols(  # DERIVED · coefficient-axis solver unknowns
            " ".join(f"dim_{coefficient}_{axis}" for axis in ("L", "T", "M"))
        )
        for coefficient in dimensionful_coefficient_symbols
    }
    unknowns = tuple(  # DERIVED · coefficient-axis solver unknown collection
        component
        for coefficient in dimensionful_coefficient_symbols
        for component in coefficient_unknowns[coefficient]
    )
    register_declared_symbols(locals())
    unknown_map = {
        **coefficient_unknowns,
        s_rho: ZERO_DIM,
        s: ZERO_DIM,
    }
    coordinates = {t: TIME_DIM, **{coordinate: LENGTH_DIM for coordinate in xvec}}
    walker = DimensionWalker(unknown_map, fields, coordinates, assumptions)
    target = ENERGY_DENSITY_DIM
    equations: list[sp.Equality] = []
    for term in sp.Add.make_args(sp.expand(lagrangian)):
        term_dimension = walker.dimension(term)
        equations.extend(
            sp.Eq(component, required, evaluate=False)
            for component, required in zip(term_dimension, target)
        )
    equations = list(dict.fromkeys(equations))
    solution = sp.Tuple(*(
        sp.Dict(item)
        for item in sp.solve(equations, unknowns, dict=True)
    ))
    coefficient_matrix, right_hand_side = sp.linear_eq_to_matrix(equations, unknowns)
    coefficient_rank = int(coefficient_matrix.rank())
    augmented_rank = int(coefficient_matrix.row_join(right_hand_side).rank())
    independent_equation_count = augmented_rank
    difference = independent_equation_count - len(unknowns)
    determination = Str(  # DERIVED · constraint-system classification marker
        "over_determined" if difference > 0 else "exactly_determined" if difference == 0 else "under_determined"
    )
    consistent = sp.true if coefficient_rank == augmented_rank else sp.false
    if solution:
        solved_map = {
            coefficient: tuple(
                solution[0].get(symbol, symbol)
                for symbol in coefficient_unknowns[coefficient]
            )
            for coefficient in dimensionful_coefficient_symbols
        }
        solved_map.update({
            s_rho: ZERO_DIM,
            s: ZERO_DIM,
        })
    else:
        solved_map = {
            **coefficient_unknowns,
            s_rho: ZERO_DIM,
            s: ZERO_DIM,
        }
    return (
        equations,
        unknowns,
        solution,
        solved_map,
        independent_equation_count,
        difference,
        determination,
        consistent,
    )


def factor_with_timeout(
    expr: sp.Expr, seconds: int = OPERATION_TIMEOUT_SECONDS
) -> tuple[sp.Expr, Str]:
    previous_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(seconds)
    try:
        result = sp.factor(expr)
        route = Str("factor_returned")  # DERIVED · determinant-factorization route record
    except OperationTimeout:
        result = expr
        route = Str("factor_timeout_unfactored")  # DERIVED · determinant-factorization route record
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_handler)
    return result, route


def derived_sign(expr: sp.Expr, assumptions: sp.logic.boolalg.Boolean) -> object:
    refined = sp.refine(sp.factor(expr), assumptions)
    positive = sp.ask(sp.Q.positive(refined), assumptions)
    zero = sp.ask(sp.Q.zero(refined), assumptions)
    negative = sp.ask(sp.Q.negative(refined), assumptions)
    if positive is True:
        return sp.Integer(1)
    if zero is True:
        return sp.Integer(0)
    if negative is True:
        return sp.Integer(-1)
    return Str("undecided_under_joint_assumptions")  # DERIVED · sign-query status marker


def deduplicate_roots(
    raw_solutions: Sequence[dict[sp.Symbol, sp.Expr]],
    assumptions: sp.logic.boolalg.Boolean,
) -> list[sp.Expr]:
    roots: list[sp.Expr] = []
    for solution in raw_solutions:
        if omegaSquared not in solution:
            continue
        root = assumed_simplify(solution[omegaSquared], assumptions)
        if not isinstance(root, sp.Expr):
            raise RuntimeError("spectrum solver returned a non-expression root")
        duplicate = False
        for existing in roots:
            difference = assumed_simplify(root - existing, assumptions)
            duplicate = difference == 0 or sp.ask(sp.Q.zero(root - existing), assumptions) is True
            if duplicate:
                break
        if not duplicate:
            roots.append(root)
    return roots


def solve_spectrum(
    matrix: sp.Matrix,
    assumptions: sp.logic.boolalg.Boolean,
) -> tuple[sp.Expr, sp.Expr, Str, list[dict[sp.Symbol, sp.Expr]]]:
    determinant = sp.cancel(matrix.det(method="berkowitz"))
    factored, factor_route = factor_with_timeout(determinant)
    raw_solutions = sp.solve(factored, omegaSquared, dict=True)
    return determinant, factored, factor_route, raw_solutions


@dataclass
class RealLocus:
    equations: tuple[sp.Expr, ...]
    raw_solver_output: object
    branches: tuple[dict[sp.Symbol, sp.Expr], ...]
    branch_conditions: tuple[tuple[sp.Expr, ...], ...]
    reality_filter: sp.Tuple
    solver_route: Str
    fallback_expression: object = sp.EmptySet

    def branch_expression(self, index: int) -> sp.Tuple:
        branch = self.branches[index]
        assignments = tuple(
            sp.Eq(symbol, value, evaluate=False)
            for symbol, value in sorted(branch.items(), key=lambda item: str(item[0]))
        )
        return sp.Tuple(*assignments, *self.branch_conditions[index])

    def expression(self) -> object:
        if self.branches:
            return sp.Tuple(*(self.branch_expression(index) for index in range(len(self.branches))))
        return self.fallback_expression


def _zero_component_deductions(
    expression: sp.Expr,
    variables: Sequence[sp.Symbol],
    assumptions: sp.logic.boolalg.Boolean,
) -> dict[sp.Symbol, sp.Expr]:
    reduced = algebraic_normalize(expression, assumptions)
    if isinstance(reduced, sp.Mul):
        dependent_factors = [factor for factor in reduced.args if factor.has(*variables)]
        independent_factor = sp.Mul(
            *(factor for factor in reduced.args if not factor.has(*variables))
        )
        if (
            len(dependent_factors) == 1
            and sp.ask(sp.Q.nonzero(independent_factor), assumptions) is True
        ):
            reduced = dependent_factors[0]
    if isinstance(reduced, sp.Pow) and reduced.exp.is_positive:
        reduced = reduced.base
    try:
        polynomial = sp.Poly(reduced, *variables)
    except sp.PolynomialError:
        return {}
    terms = polynomial.terms()
    if not terms:
        return {}
    if len(terms) == 1:
        monomial, coefficient = terms[0]
        active_indices = [index for index, exponent in enumerate(monomial) if exponent != 0]
        if (
            len(active_indices) == 1
            and sp.ask(sp.Q.nonzero(coefficient), assumptions) is True
        ):
            return {variables[active_indices[0]]: sp.Integer(0)}
        return {}
    active_variables: list[sp.Symbol] = []
    coefficient_signs: list[int] = []
    for monomial, coefficient in terms:
        active_indices = [index for index, exponent in enumerate(monomial) if exponent != 0]
        if len(active_indices) != 1:
            return {}
        active_index = active_indices[0]
        if monomial[active_index] % 2 != 0:
            return {}
        if sp.ask(sp.Q.positive(coefficient), assumptions) is True:
            coefficient_signs.append(1)
        elif sp.ask(sp.Q.negative(coefficient), assumptions) is True:
            coefficient_signs.append(-1)
        else:
            return {}
        active_variables.append(variables[active_index])
    if len(set(coefficient_signs)) != 1:
        return {}
    return {symbol: sp.Integer(0) for symbol in active_variables}


def solve_real_locus(
    equations: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    assumptions: sp.logic.boolalg.Boolean,
    impossible: bool = False,
) -> RealLocus:
    normalized: list[sp.Expr] = []
    for equation in equations:
        numerator = sp.together(equation).as_numer_denom()[0]
        simplified = algebraic_normalize(numerator, assumptions)
        if isinstance(simplified, sp.Expr) and simplified != 0 and simplified not in normalized:
            normalized.append(simplified)
    if impossible:
        return RealLocus(
            tuple(normalized),
            sp.EmptySet,
            (),
            (),
            sp.Tuple(),
            Str("rank_floor_empty_locus"),  # DERIVED · real-locus solver route marker
        )

    if not normalized:
        real_domain = sp.ProductSet(*(sp.S.Reals for _ in variables))
        return RealLocus(
            (),
            real_domain,
            ({},),
            ((),),
            sp.Tuple(sp.Tuple(sp.Tuple(), sp.Tuple(), Str("retained_real_branch"))),  # DERIVED · real-branch filter record
            Str("universal_real_locus"),  # DERIVED · real-locus solver route marker
            real_domain,
        )

    previous_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(OPERATION_TIMEOUT_SECONDS)
    try:
        raw = sp.solve(normalized, list(variables), dict=True)
        route = Str("solve_then_explicit_real_filter")  # DERIVED · real-locus solver route record
        fallback_expression: object = sp.EmptySet
    except OperationTimeout:
        raw = sp.ConditionSet(
            sp.Tuple(*variables),
            sp.And(*(sp.Eq(eq, 0) for eq in normalized)),
            sp.ProductSet(*(sp.S.Reals for _ in variables)),
        )
        route = Str("solve_timeout_explicit_real_conditionset")  # DERIVED · real-locus solver route record
        fallback_expression = raw
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_handler)

    raw_branches: list[dict[sp.Symbol, sp.Expr]] = []
    if isinstance(raw, dict):
        raw_branches = [raw]
    elif isinstance(raw, list):
        raw_branches = [branch for branch in raw if isinstance(branch, dict)]

    branches: list[dict[sp.Symbol, sp.Expr]] = []
    branch_conditions: list[tuple[sp.Expr, ...]] = []
    filter_records: list[sp.Tuple] = []
    seen: set[tuple[tuple[str, str], tuple[str, ...]]] = set()
    for raw_branch in raw_branches:
        raw_normalized_branch = {
            symbol: algebraic_normalize(value, assumptions)
            for symbol, value in raw_branch.items()
        }
        conditions: list[sp.Expr] = []
        for value in raw_normalized_branch.values():
            imaginary_part = sp.simplify(sp.im(value))
            if imaginary_part != 0:
                conditions.append(sp.Eq(imaginary_part, 0, evaluate=False))
            denominator = sp.together(value).as_numer_denom()[1]
            if denominator != 1:
                conditions.append(sp.Ne(denominator, 0, evaluate=False))
        for equation in normalized:
            residual = algebraic_normalize(
                equation.subs(raw_normalized_branch, simultaneous=True), assumptions
            )
            if residual != 0:
                conditions.append(sp.Eq(residual, 0, evaluate=False))
        conditions = list(dict.fromkeys(conditions))
        raw_conditions = tuple(conditions)
        raw_assignments = sp.Tuple(
            *(
                sp.Eq(symbol, value, evaluate=False)
                for symbol, value in sorted(
                    raw_normalized_branch.items(), key=lambda item: str(item[0])
                )
            )
        )
        if any(condition is sp.false for condition in conditions):
            filter_records.append(
                sp.Tuple(
                    raw_assignments,
                    sp.Tuple(*conditions),
                    Str("discarded_nonreal_branch"),  # DERIVED · real-branch filter status
                )
            )
            continue

        deductions: dict[sp.Symbol, sp.Expr] = {}
        for condition in conditions:
            if isinstance(condition, sp.Equality):
                deductions.update(
                    _zero_component_deductions(
                        condition.lhs - condition.rhs, variables, assumptions
                    )
                )
        branch = dict(raw_normalized_branch)
        branch.update(deductions)
        for _ in range(len(variables) + 1):
            branch = {
                symbol: algebraic_normalize(value.subs(branch, simultaneous=True), assumptions)
                for symbol, value in branch.items()
            }
            branch.update(deductions)
        projected_conditions: list[sp.Expr] = []
        for condition in conditions:
            if isinstance(condition, (sp.Equality, sp.Unequality)):
                residual = algebraic_normalize(
                    (condition.lhs - condition.rhs).subs(branch), assumptions
                )
                if isinstance(condition, sp.Equality) and residual == 0:
                    continue
                if isinstance(condition, sp.Unequality) and residual != 0:
                    continue
            projected_conditions.append(condition.subs(branch))
        conditions = list(dict.fromkeys(projected_conditions))
        assignments = sp.Tuple(
            *(
                sp.Eq(symbol, value, evaluate=False)
                for symbol, value in sorted(branch.items(), key=lambda item: str(item[0]))
            )
        )
        key = (
            tuple(sorted((str(symbol), str(value)) for symbol, value in branch.items())),
            tuple(str(condition) for condition in conditions),
        )
        filter_records.append(
            sp.Tuple(
                raw_assignments,
                sp.Tuple(*raw_conditions),
                Str("retained_real_branch"),  # DERIVED · real-branch filter status
                assignments,
            )
        )
        if key not in seen:
            seen.add(key)
            branches.append(branch)
            branch_conditions.append(tuple(conditions))
    return RealLocus(
        tuple(normalized),
        raw,
        tuple(branches),
        tuple(branch_conditions),
        sp.Tuple(*filter_records),
        route,
        fallback_expression,
    )


@dataclass
class LocusAllowedData:
    operands: sp.Tuple
    allowed: object
    witnesses: sp.Tuple
    branch_operands: sp.Tuple
    branch_tests: sp.Tuple
    branch_witnesses: sp.Tuple
    branch_points: tuple[dict[sp.Symbol, sp.Expr] | None, ...]
    branch_reasons: sp.Tuple


def _materialize_branch_point(
    branch: Mapping[sp.Symbol, sp.Expr],
    variables: Sequence[sp.Symbol],
    free_values: Mapping[sp.Symbol, sp.Expr],
    assumptions: sp.logic.boolalg.Boolean,
) -> dict[sp.Symbol, sp.Expr] | None:
    point = dict(free_values)
    for _ in range(len(variables) + 1):
        for symbol, value in branch.items():
            resolved = assumed_simplify(value.subs(point, simultaneous=True), assumptions)
            if not isinstance(resolved, sp.Expr):
                return None
            point[symbol] = resolved
    if any(component not in point for component in variables):
        return None
    if any(point[component].has(*variables) for component in variables):
        return None
    return point


def _condition_at_point(
    condition: sp.Expr,
    point: Mapping[sp.Symbol, sp.Expr],
    assumptions: sp.logic.boolalg.Boolean,
) -> bool | None:
    if isinstance(condition, sp.Equality):
        residual = algebraic_normalize((condition.lhs - condition.rhs).subs(point), assumptions)
        if residual == 0:
            return True
        if not residual.free_symbols:
            return False
        answer = sp.ask(sp.Q.zero(residual), assumptions)
        return answer
    if isinstance(condition, sp.Unequality):
        residual = algebraic_normalize((condition.lhs - condition.rhs).subs(point), assumptions)
        if residual == 0:
            return False
        if not residual.free_symbols:
            return True
        answer = sp.ask(sp.Q.nonzero(residual), assumptions)
        return answer
    refined = sp.refine(condition.subs(point), assumptions)
    if refined is sp.true:
        return True
    if refined is sp.false:
        return False
    return None


def locus_allowed_data(
    locus: RealLocus,
    kvec: Sequence[sp.Symbol],
    k_squared: sp.Expr,
    assumptions: sp.logic.boolalg.Boolean,
) -> LocusAllowedData:
    operands = sp.Tuple(locus.expression(), sp.Gt(k_squared, 0, evaluate=False), assumptions)
    branch_operands: list[sp.Tuple] = []
    branch_tests: list[sp.Expr] = []
    branch_witnesses: list[sp.Tuple] = []
    branch_points: list[dict[sp.Symbol, sp.Expr] | None] = []
    branch_reasons: list[Str] = []
    witnesses: list[sp.Tuple] = []

    norm_conflict = False
    for equation in locus.equations:
        zero_deductions = _zero_component_deductions(equation, kvec, assumptions)
        if set(zero_deductions) == set(kvec):
            norm_conflict = True
            break
        quotient = sp.cancel(equation / k_squared)
        if quotient != 0 and not quotient.has(*kvec):
            nonzero = sp.ask(sp.Q.nonzero(quotient), assumptions)
            if nonzero is True:
                norm_conflict = True
                break

    for branch_index, branch in enumerate(locus.branches):
        branch_expression = locus.branch_expression(branch_index)
        branch_operands.append(
            sp.Tuple(branch_expression, sp.Gt(k_squared, 0, evaluate=False), assumptions)
        )
        if norm_conflict:
            branch_tests.append(sp.false)
            branch_witnesses.append(sp.Tuple())
            branch_points.append(None)
            branch_reasons.append(Str("locus_conflicts_with_positive_wavevector_norm"))  # DERIVED · locus-admissibility status
            continue
        free = [component for component in kvec if component not in branch]
        trial_values: list[dict[sp.Symbol, sp.Expr]] = []
        zero_values = {component: sp.Integer(0) for component in free}
        trial_values.append(zero_values)
        for selected in free:
            trial_values.append(
                {
                    component: sp.Integer(1) if component == selected else sp.Integer(0)
                    for component in free
                }
            )
        if free:
            trial_values.append({component: sp.Integer(1) for component in free})
            trial_values.append({component: sp.Integer(2) for component in free})

        witness_point: dict[sp.Symbol, sp.Expr] | None = None
        for free_values in trial_values:
            point = _materialize_branch_point(branch, kvec, free_values, assumptions)
            if point is None:
                continue
            equation_values = [
                algebraic_normalize(equation.subs(point), assumptions)
                for equation in locus.equations
            ]
            condition_values = [
                _condition_at_point(condition, point, assumptions)
                for condition in locus.branch_conditions[branch_index]
            ]
            k_value = assumed_simplify(k_squared.subs(point), assumptions)
            if (
                all(value == 0 for value in equation_values)
                and all(value is True for value in condition_values)
                and sp.ask(sp.Q.positive(k_value), assumptions) is True
            ):
                witness_point = point
                break

        if witness_point is not None:
            witness = sp.Tuple(
                *(
                    sp.Eq(component, witness_point[component], evaluate=False)
                    for component in kvec
                )
            )
            witnesses.append(witness)
            branch_tests.append(sp.true)
            branch_witnesses.append(witness)
            branch_points.append(witness_point)
            branch_reasons.append(Str("allowed_witness_found"))  # DERIVED · locus-admissibility status
        elif not free:
            branch_tests.append(Str("undecided_no_explicit_real_witness"))  # DERIVED · locus-admissibility status
            branch_witnesses.append(sp.Tuple())
            branch_points.append(None)
            branch_reasons.append(Str("fully_fixed_branch_not_materialized"))  # DERIVED · locus-admissibility status
        else:
            branch_tests.append(Str("undecided_no_explicit_real_witness"))  # DERIVED · locus-admissibility status
            branch_witnesses.append(sp.Tuple())
            branch_points.append(None)
            branch_reasons.append(Str("no_explicit_real_witness"))  # DERIVED · locus-admissibility status

    if any(test is sp.true for test in branch_tests):
        allowed: object = sp.true
    elif branch_tests and all(test is sp.false for test in branch_tests):
        allowed = sp.false
    elif locus.fallback_expression is not sp.EmptySet and not locus.branches:
        allowed = Str("undecided_solver_conditionset")  # DERIVED · locus-admissibility status
    else:
        allowed = Str("undecided_no_explicit_real_witness") if branch_tests else sp.false  # DERIVED · locus-admissibility status
    return LocusAllowedData(
        operands,
        allowed,
        sp.Tuple(*witnesses),
        sp.Tuple(*branch_operands),
        sp.Tuple(*branch_tests),
        sp.Tuple(*branch_witnesses),
        tuple(branch_points),
        sp.Tuple(*branch_reasons),
    )


def rank_drop_minors(matrix: sp.Matrix, generic_rank: int) -> list[sp.Expr]:
    if generic_rank == 0:
        return []
    minors: list[sp.Expr] = []
    seen: set[sp.Expr] = set()
    for rows in itertools.combinations(range(matrix.rows), generic_rank):
        for columns in itertools.combinations(range(matrix.cols), generic_rank):
            minor = sp.factor(matrix.extract(rows, columns).det(method="berkowitz"))
            if minor not in seen:
                seen.add(minor)
                minors.append(minor)
    return minors


def q4_objects(
    matrix: sp.Matrix,
    root: sp.Expr,
    k_column: sp.Matrix,
    k_squared: sp.Expr,
    assumptions: sp.logic.boolalg.Boolean,
) -> dict[str, object]:
    root_matrix = sp.Matrix(matrix.subs(omegaSquared, root)).applyfunc(
        lambda item: algebraic_normalize(item, assumptions)
    )
    rank_input = exact_row_equivalent(root_matrix)
    rank = int(rank_input.rank(iszerofunc=polynomial_iszero, simplify=False))
    nullity = sp.Integer(root_matrix.cols - rank)
    stacked = root_matrix.col_join(sp.Matrix([[*list(k_column)]]))
    stacked_rank_input = exact_row_equivalent(stacked)
    stacked_rank = int(stacked_rank_input.rank(iszerofunc=polynomial_iszero, simplify=False))
    transverse_nullity = sp.Integer(root_matrix.cols - stacked_rank)
    matrix_times_k = assumed_simplify(root_matrix * k_column, assumptions)
    basis = rank_input.nullspace(simplify=False, iszerofunc=polynomial_iszero)
    direct_nullspace = sp.polys.matrices.DomainMatrix.from_Matrix(root_matrix).nullspace()
    dots: list[sp.Expr] = []
    residuals: list[sp.Matrix] = []
    for vector in basis:
        dot = assumed_simplify((vector.T * k_column)[0], assumptions)
        residual = assumed_simplify(k_squared * vector - dot * k_column, assumptions)
        if not isinstance(dot, sp.Expr) or not isinstance(residual, sp.MatrixBase):
            raise RuntimeError("null-space display computation returned an unsupported object")
        dots.append(dot)
        residuals.append(sp.Matrix(residual))
    basis_count = sp.Integer(direct_nullspace.shape[0])
    return {
        "N1_MATRIX": root_matrix,
        "N2_RANK": sp.Integer(rank),
        "N2_NULLITY": nullity,
        "N3_STACKED_MATRIX": stacked,
        "N3_STACKED_RANK": sp.Integer(stacked_rank),
        "N3_TRANSVERSE_NULLITY": transverse_nullity,
        "N4_NULLITY_DIFFERENCE": nullity - transverse_nullity,
        "N5_MATRIX_TIMES_K": matrix_times_k,
        "N6_NULLSPACE_BASIS": sp.Tuple(*basis),
        "N6_BASIS_DOT_K": sp.Tuple(*dots),
        "N6_BASIS_VECTOR_RESIDUALS": sp.Tuple(*residuals),
        "N7_BASIS_COUNT": basis_count,
        "N7_BASIS_COUNT_RESIDUAL": basis_count - nullity,
    }


def emit_q4(
    emitter: Emitter,
    walker: DimensionWalker,
    root_prefix: str,
    objects: Mapping[str, object],
    dimension_objects: Mapping[str, object] | None = None,
) -> None:
    physical = {
        "N1_MATRIX",
        "N3_STACKED_MATRIX",
        "N5_MATRIX_TIMES_K",
        "N6_NULLSPACE_BASIS",
        "N6_BASIS_DOT_K",
        "N6_BASIS_VECTOR_RESIDUALS",
    }
    for suffix, payload in objects.items():
        name = root_prefix + suffix
        if suffix in physical:
            emit_physical(
                emitter,
                walker,
                name,
                payload,
                dimension_payload=(
                    dimension_objects[suffix]
                    if dimension_objects is not None
                    else None
                ),
            )
        else:
            emitter.emit(name, payload)


def root_scale_objects(
    root: sp.Expr,
    kvec: Sequence[sp.Symbol],
    assumptions: sp.logic.boolalg.Boolean,
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.logic.boolalg.Boolean, sp.Expr]:
    substitution = {component: lambdaScale * component for component in kvec}
    scaled = assumed_simplify(root.subs(substitution, simultaneous=True), assumptions)
    original = assumed_simplify(root, assumptions)
    if not isinstance(scaled, sp.Expr) or not isinstance(original, sp.Expr):
        raise RuntimeError("root scaling returned a non-expression")
    is_zero = sp.ask(sp.Q.zero(original), assumptions)
    if is_zero is True:
        ratio = Str("undefined_zero_root_ratio")  # DERIVED · root-scaling domain status
        ratio_defined = sp.false
        exponent = Str("undefined_zero_root_ratio")  # DERIVED · root-scaling classification status
    else:
        ratio_object = sp.Mul(scaled, sp.Pow(original, -1, evaluate=False), evaluate=False)
        ratio = assumed_simplify(ratio_object, assumptions)
        ratio_defined = sp.true
        try:
            polynomial = sp.Poly(ratio, lambdaScale)
            exponent = (
                sp.Integer(polynomial.degree())
                if len(polynomial.terms()) == 1
                else Str("not_a_pure_lambdaScale_power")  # DERIVED · root-scaling classification status
            )
        except (sp.PolynomialError, TypeError, ValueError):
            exponent = Str("not_a_pure_lambdaScale_power")  # DERIVED · root-scaling classification status
    if ratio_defined is sp.true and not isinstance(ratio, sp.Expr):
        raise RuntimeError("root scaling ratio returned a non-expression")
    return scaled, original, ratio, ratio_defined, exponent


def emit_premises(
    emitter: Emitter,
    prefix: str,
    ansatz: Sequence[sp.Expr],
    kvec: Sequence[sp.Symbol],
    avec: Sequence[sp.Symbol],
    assumptions: sp.logic.boolalg.Boolean,
    package: Package,
    n: int,
) -> None:
    emitter.emit(prefix + "PREMISE_FIELD_SECTOR", sp.Eq(sp.Symbol("u_component_count"), n, evaluate=False))  # PREMISE · supplied field-sector declaration
    emitter.emit(prefix + "PREMISE_U_DIMENSION", dim_expr(FIELD_DIM))
    emitter.emit(prefix + "PREMISE_ANSATZ", sp.Tuple(*ansatz))
    emitter.emit(
        prefix + "PREMISE_PERIOD_AVERAGE",
        sp.Tuple(sp.Symbol("phase"), sp.Integer(0), 2 * sp.pi, sp.Rational(1, 2) / sp.pi),  # PREMISE · supplied phase-average prescription
    )
    emitter.emit(prefix + "PREMISE_BACKGROUND_VELOCITY", sp.Eq(sp.Symbol("v_0"), 0, evaluate=False))  # PREMISE · supplied background-velocity condition
    emitter.emit(prefix + "PREMISE_TIME_ODD_KERNEL", sp.Eq(sp.Symbol("time_odd_kernel"), 0, evaluate=False))  # PREMISE · supplied time-parity condition
    emitter.emit(prefix + "PREMISE_RESPONSE_DEGREE", sp.Eq(sp.Symbol("response_degree"), 2, evaluate=False))  # PREMISE · supplied response-order condition
    emitter.emit(prefix + "PREMISE_STIFFNESS_INPUT", Str("S_curl_supplied"))  # PREMISE · supplied stiffness-form label
    emitter.emit(prefix + "PREMISE_RHO_DOMAIN", sp.Q.positive(rho_br))
    emitter.emit(prefix + "PREMISE_MU_DOMAIN", sp.Q.positive(mu_R))
    emitter.emit(prefix + "PREMISE_WAVEVECTOR_NORM_DOMAIN", sp.Q.positive(sp.Add(*(item**2 for item in kvec))))
    emitter.emit(prefix + "PREMISE_WAVEVECTOR_REAL_DOMAIN", sp.And(*(sp.Q.real(item) for item in kvec)))
    emitter.emit(prefix + "PREMISE_AMPLITUDE_REAL_DOMAIN", sp.And(*(sp.Q.real(item) for item in avec)))
    emitter.emit(prefix + "PREMISE_BRANE_DIMENSION_DOMAIN", sp.And(sp.Q.integer(D), sp.Q.positive(D)))
    emitter.emit(prefix + "ASSUMPTION_JOINT_PREDICATE", assumptions)
    controls: list[sp.Expr] = []
    if package.anisotropic_kinetic:
        controls.append(sp.Eq(sp.Symbol("dimension_s_rho"), sp.Symbol("dimensionless"), evaluate=False))  # PREMISE · supplied anisotropy-control unit condition
    if package.coefficient_scale:
        controls.append(sp.Eq(sp.Symbol("dimension_s"), sp.Symbol("dimensionless"), evaluate=False))  # PREMISE · supplied scaling-control unit condition
    emitter.emit(prefix + "Q6_CONTROL_DIMENSION_PREMISES", sp.Tuple(*controls))


def emit_q3(
    emitter: Emitter,
    walker: DimensionWalker,
    prefix: str,
    matrix: sp.Matrix,
    roots_data: tuple[sp.Expr, sp.Expr, Str, list[dict[sp.Symbol, sp.Expr]]],
    assumptions: sp.logic.boolalg.Boolean,
    dimension_roots_data: tuple[
        sp.Expr, sp.Expr, sp.Symbol, list[dict[sp.Symbol, sp.Expr]]
    ]
    | None = None,
) -> list[sp.Expr]:
    determinant, factored, factor_route, raw_solutions = roots_data
    dimension_determinant = (
        dimension_roots_data[0] if dimension_roots_data is not None else None
    )
    dimension_factored = (
        dimension_roots_data[1] if dimension_roots_data is not None else None
    )
    emit_physical(
        emitter,
        walker,
        prefix + "Q3_DETERMINANT",
        determinant,
        homogeneity_class="SOLVED",
        dimension_payload=dimension_determinant,
    )
    emit_physical(
        emitter,
        walker,
        prefix + "Q3_DETERMINANT_FACTORED",
        factored,
        homogeneity_class="SOLVED",
        dimension_payload=dimension_factored,
    )
    emitter.emit("PY_S10_LOCAL_" + prefix.removeprefix("PY_S10_") + "Q3_DETERMINANT_FACTOR_ROUTE", factor_route)
    emitter.emit(prefix + "Q3_ROOT_SOLUTIONS_RAW", raw_solutions)
    roots = deduplicate_roots(raw_solutions, assumptions)
    emitter.emit(prefix + "Q3_ROOTS_DISTINCT", sp.Tuple(*roots))
    emitter.emit(prefix + "Q3_ROOT_COUNT", sp.Integer(len(roots)))
    emitter.emit(prefix + "ROOT_ORDERING", sp.Tuple(*roots))
    if not roots:
        emitter.emit(prefix + "Q3_SPECTRUM_SOLVE_CONDITION", Str("no_roots_returned"))  # DERIVED · spectrum-solver status
        emitter.emit(
            prefix + "Q3_SPECTRUM_SOLVE_CONDITION_OPERANDS",
            sp.Tuple(factored, omegaSquared),
        )
    else:
        emitter.emit(prefix + "Q3_SPECTRUM_SOLVE_CONDITION", Str("roots_returned"))  # DERIVED · spectrum-solver status
        emitter.emit(
            prefix + "Q3_SPECTRUM_SOLVE_CONDITION_OPERANDS",
            sp.Tuple(factored, omegaSquared),
        )
    return roots


def q3_coincidence_data(
    roots: Sequence[sp.Expr],
    kvec: Sequence[sp.Symbol],
    k_squared: sp.Expr,
    assumptions: sp.logic.boolalg.Boolean,
) -> list[tuple[tuple[int, int], sp.Expr, RealLocus, LocusAllowedData]]:
    result = []
    for left, right in itertools.combinations(range(len(roots)), 2):
        equation = assumed_simplify(roots[left] - roots[right], assumptions)
        if not isinstance(equation, sp.Expr):
            raise RuntimeError("root-coincidence equation is not an expression")
        locus = solve_real_locus([equation], kvec, assumptions)
        allowed_data = locus_allowed_data(locus, kvec, k_squared, assumptions)
        result.append(((left + 1, right + 1), equation, locus, allowed_data))
    return result


def emit_coincidences(
    emitter: Emitter,
    walker: DimensionWalker,
    prefix: str,
    data: Sequence[tuple[tuple[int, int], sp.Expr, RealLocus, LocusAllowedData]],
    dimension_equations: Sequence[sp.Expr] | None = None,
) -> None:
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_PAIR_INDICES", sp.Tuple(*(sp.Tuple(*pair) for pair, *_ in data)))
    equations = sp.Tuple(*(equation for _, equation, *_ in data))
    emit_physical(
        emitter,
        walker,
        prefix + "Q3_ROOT_COINCIDENCE_EQUATIONS",
        equations,
        dimension_payload=(
            sp.Tuple(*dimension_equations)
            if dimension_equations is not None
            else None
        ),
    )
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_LOCI", sp.Tuple(*(locus.expression() for _, _, locus, *_ in data)))
    emitter.emit(
        prefix + "Q3_ROOT_COINCIDENCE_REALITY_FILTERS",
        sp.Tuple(*(locus.reality_filter for _, _, locus, _ in data)),
    )
    emitter.emit(
        prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_OPERANDS",
        sp.Tuple(*(allowed_data.operands for *_, allowed_data in data)),
    )
    emitter.emit(
        prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_TESTS",
        sp.Tuple(*(allowed_data.allowed for *_, allowed_data in data)),
    )
    emitter.emit(
        prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_WITNESSES",
        sp.Tuple(*(allowed_data.witnesses for *_, allowed_data in data)),
    )
    emitter.emit(
        prefix + "Q3_ROOT_COINCIDENCE_BRANCH_ALLOWED_OPERANDS",
        sp.Tuple(*(allowed_data.branch_operands for *_, allowed_data in data)),
    )
    emitter.emit(
        prefix + "Q3_ROOT_COINCIDENCE_BRANCH_ALLOWED_TESTS",
        sp.Tuple(*(allowed_data.branch_tests for *_, allowed_data in data)),
    )
    emitter.emit(
        prefix + "Q3_ROOT_COINCIDENCE_BRANCH_ALLOWED_WITNESSES",
        sp.Tuple(*(allowed_data.branch_witnesses for *_, allowed_data in data)),
    )
    local_prefix = "PY_S10_LOCAL_" + prefix.removeprefix("PY_S10_")
    emitter.emit(local_prefix + "Q3_ROOT_COINCIDENCE_SOLVER_OUTPUTS", sp.Tuple(*(sp.sympify(locus.raw_solver_output) for _, _, locus, *_ in data)))
    emitter.emit(local_prefix + "Q3_ROOT_COINCIDENCE_SOLVER_ROUTES", sp.Tuple(*(locus.solver_route for _, _, locus, *_ in data)))


def emit_stratum_q3_q4(
    emitter: Emitter,
    walker: DimensionWalker,
    prefix: str,
    matrix: sp.Matrix,
    point: Mapping[sp.Symbol, sp.Expr],
    kvec: Sequence[sp.Symbol],
    assumptions: sp.logic.boolalg.Boolean,
) -> None:
    point_matrix = sp.Matrix(matrix.subs(point))
    point_k = sp.Matrix([component.subs(point) for component in kvec])
    point_k_squared = assumed_simplify((point_k.T * point_k)[0], assumptions)
    point_assumptions = sp.refine(assumptions.subs(point), assumptions)
    wave_number_scale = sp.Symbol(  # COORDINATE · stratum wave-number scaling coordinate
        prefix.removeprefix("PY_S10_") + "WAVE_NUMBER_SCALE",
        positive=True,
        real=True,
    )
    dimension_point = {
        component: point[component] * wave_number_scale for component in kvec
    }
    dimension_matrix = sp.Matrix(matrix.subs(dimension_point))
    dimension_k = sp.Matrix(
        [component.subs(dimension_point) for component in kvec]
    )
    dimension_k_squared = assumed_simplify(
        (dimension_k.T * dimension_k)[0], assumptions
    )
    dimension_walker = DimensionWalker(
        {**walker.symbol_dimensions, wave_number_scale: WAVENUMBER_DIM},
        tuple(walker.fields),
        walker.coordinate_dimensions,
        assumptions,
    )
    if not isinstance(point_k_squared, sp.Expr):
        raise RuntimeError("stratum wavevector norm is not an expression")
    if not isinstance(dimension_k_squared, sp.Expr):
        raise RuntimeError("dimension-preserving stratum wavevector norm is not an expression")
    spectrum = solve_spectrum(point_matrix, point_assumptions)
    dimension_spectrum = solve_spectrum(dimension_matrix, assumptions)
    roots = emit_q3(
        emitter,
        dimension_walker,
        prefix,
        point_matrix,
        spectrum,
        point_assumptions,
        dimension_roots_data=dimension_spectrum,
    )
    dimension_roots = deduplicate_roots(dimension_spectrum[3], assumptions)
    for index, root in enumerate(roots, 1):
        root_prefix = prefix + f"ROOT{index}_"
        dimension_root = dimension_roots[index - 1]
        emit_physical(
            emitter,
            dimension_walker,
            root_prefix + "Q3_ROOT",
            root,
            dimension_payload=dimension_root,
        )
        emitter.emit(root_prefix + "Q3_SIGN", derived_sign(root, point_assumptions))
        objects = q4_objects(point_matrix, root, point_k, point_k_squared, point_assumptions)
        dimension_objects = q4_objects(
            dimension_matrix,
            dimension_root,
            dimension_k,
            dimension_k_squared,
            assumptions,
        )
        emit_q4(
            emitter,
            dimension_walker,
            root_prefix,
            objects,
            dimension_objects,
        )
    coincidence_data = q3_coincidence_data(roots, kvec, point_k_squared, point_assumptions)
    dimension_coincidence_equations = [
        assumed_simplify(left - right, assumptions)
        for left, right in itertools.combinations(dimension_roots, 2)
    ]
    emit_coincidences(
        emitter,
        dimension_walker,
        prefix,
        coincidence_data,
        dimension_coincidence_equations,
    )


def run_package_dimension(
    emitter: Emitter,
    package: Package,
    n: int,
) -> dict[sp.Symbol, tuple[sp.Expr, ...]]:
    prefix = f"PY_S10_{package.name}_D{n}_"
    t = sp.Symbol("t", real=True)  # COORDINATE · time coordinate
    xvec = sp.symbols(f"x1:{n + 1}", real=True)  # COORDINATE · spatial coordinate collection
    kvec = sp.symbols(f"k1:{n + 1}", real=True)  # COORDINATE · wavevector component collection
    avec = sp.symbols(f"a1:{n + 1}", real=True)  # PREMISE · plane-wave amplitude collection
    phase = sp.Symbol("phase", real=True)  # COORDINATE · period-average phase coordinate
    theta = sp.Add(*(kvec[index] * xvec[index] for index in range(n))) - sp.sqrt(omegaSquared) * t
    ansatz = [avec[index] * sp.cos(theta) for index in range(n)]
    k_column = sp.Matrix(kvec)
    k_squared = sp.Add(*(component**2 for component in kvec))
    assumptions = build_joint_assumptions(package, kvec, avec)

    fields, lagrangian, inertial_coefficients, stiffness_coefficients, stiffness = build_action(package, n, t, xvec)
    (
        dimension_equations,
        dimension_unknowns,
        dimension_solution,
        solved_coefficient_dims,
        independent_dimension_equation_count,
        dimension_count_difference,
        dimension_determination,
        dimension_consistent,
    ) = derive_dimension_solution(package, lagrangian, fields, t, xvec, assumptions)
    symbol_dimensions = {
        **solved_coefficient_dims,
        omegaSquared: OMEGA_SQUARED_DIM,
        lambdaScale: ZERO_DIM,
        D: ZERO_DIM,
        **{component: WAVENUMBER_DIM for component in kvec},
        **{component: LENGTH_DIM for component in avec},
        **{coordinate: LENGTH_DIM for coordinate in xvec},
        t: TIME_DIM,
    }
    walker = DimensionWalker(
        symbol_dimensions,
        fields,
        {t: TIME_DIM, **{coordinate: LENGTH_DIM for coordinate in xvec}},
        assumptions,
    )

    emit_premises(emitter, prefix, ansatz, kvec, avec, assumptions, package, n)
    emit_physical(
        emitter,
        walker,
        prefix + "Q1_LAGRANGIAN_EXPANDED",
        lagrangian,
        homogeneity_class="SOLVED",
    )
    equations = euler_lagrange_system(lagrangian, fields, t, xvec)
    emit_physical(
        emitter,
        walker,
        prefix + "Q1_EULER_LAGRANGE_SYSTEM",
        equations,
        homogeneity_class="SOLVED",
    )

    ansatz_equations = substitute_ansatz(equations, fields, ansatz)
    if not isinstance(ansatz_equations, sp.MatrixBase):
        raise RuntimeError("route-A ansatz equations are not a matrix")
    expanded_ansatz_equations = [sp.expand(item) for item in ansatz_equations]
    amplitude_equations = sp.Matrix(
        [assumed_simplify(item.coeff(sp.cos(theta)), assumptions) for item in expanded_ansatz_equations]
    )
    route_a_discarded_remainder = sp.Matrix(
        [
            assumed_simplify(item - coefficient * sp.cos(theta), assumptions)
            for item, coefficient in zip(expanded_ansatz_equations, amplitude_equations)
        ]
    )
    matrix_a = sp.Matrix(amplitude_equations).jacobian(avec)

    ansatz_lagrangian = substitute_ansatz(lagrangian, fields, ansatz)
    if not isinstance(ansatz_lagrangian, sp.Expr):
        raise RuntimeError("route-B ansatz Lagrangian is not an expression")
    phase_density = sp.expand(ansatz_lagrangian.xreplace({theta: phase}))
    averaged_lagrangian = sp.integrate(phase_density, (phase, 0, 2 * sp.pi)) / (2 * sp.pi)
    averaged_lagrangian = assumed_simplify(averaged_lagrangian, assumptions)
    if not isinstance(averaged_lagrangian, sp.Expr):
        raise RuntimeError("period average is not an expression")
    matrix_b = sp.hessian(averaged_lagrangian, avec)
    matrix_residual = matrix_a - matrix_b
    matrix_ratio_numerator = assumed_simplify(matrix_a[0, 0], assumptions)
    matrix_ratio_denominator = assumed_simplify(matrix_b[0, 0], assumptions)
    denominator_is_zero = matrix_ratio_denominator == 0 or sp.ask(
        sp.Q.zero(matrix_ratio_denominator), assumptions
    ) is True
    if denominator_is_zero:
        matrix_ratio = Str("undefined_zero_denominator")  # DERIVED · route-ratio domain status
    else:
        computed_ratio = assumed_simplify(
            matrix_ratio_numerator / matrix_ratio_denominator, assumptions
        )
        if isinstance(computed_ratio, sp.Expr) and computed_ratio.has(sp.zoo, sp.nan):
            matrix_ratio = Str("undefined_nonfinite_ratio")  # DERIVED · route-ratio domain status
        else:
            matrix_ratio = computed_ratio

    emit_physical(
        emitter,
        walker,
        prefix + "Q2_PERIOD_AVERAGED_LAGRANGIAN",
        averaged_lagrangian,
        homogeneity_class="SOLVED",
    )
    emit_physical(
        emitter,
        walker,
        prefix + "Q2_MATRIX_A",
        matrix_a,
        homogeneity_class="SOLVED",
    )
    emit_physical(
        emitter,
        walker,
        prefix + "Q2_ROUTE_A_DISCARDED_REMAINDER",
        route_a_discarded_remainder,
        homogeneity_class="SOLVED",
    )
    emit_physical(
        emitter,
        walker,
        prefix + "Q2_MATRIX_B",
        matrix_b,
        homogeneity_class="SOLVED",
    )
    emit_physical(
        emitter,
        walker,
        prefix + "Q2_MATRIX_RESIDUAL",
        matrix_residual,
        homogeneity_class="SOLVED",
    )
    emit_physical(
        emitter,
        walker,
        prefix + "Q2_MATRIX_ENTRY_RATIO",
        matrix_ratio,
        homogeneity_class="SOLVED",
    )
    emitter.emit(
        prefix + "Q2_MATRIX_ENTRY_RATIO_OPERANDS",
        sp.Tuple(matrix_ratio_numerator, matrix_ratio_denominator),
    )
    emitter.emit(prefix + "Q2_RESIDUAL_TEST_SCOPE", Str("same_action_variational_identity"))  # CONTROL · route-comparison scope label
    emitter.emit(prefix + "Q2_DOWNSTREAM_ROUTE", Str("M_B"))  # CONTROL · downstream matrix-route selector

    spectrum = solve_spectrum(matrix_b, assumptions)
    roots = emit_q3(emitter, walker, prefix, matrix_b, spectrum, assumptions)
    coincidences = q3_coincidence_data(roots, kvec, k_squared, assumptions)
    emit_coincidences(emitter, walker, prefix, coincidences)

    q4_by_root: list[dict[str, object]] = []
    for index, root in enumerate(roots, 1):
        root_prefix = prefix + f"ROOT{index}_"
        emit_physical(emitter, walker, root_prefix + "Q3_ROOT", root)
        emitter.emit(root_prefix + "Q3_SIGN", derived_sign(root, assumptions))
        objects = q4_objects(matrix_b, root, k_column, k_squared, assumptions)
        q4_by_root.append(objects)
        emit_q4(emitter, walker, root_prefix, objects)

        scaled, original, ratio, ratio_defined, exponent = root_scale_objects(root, kvec, assumptions)
        emit_physical(emitter, walker, root_prefix + "Q5_ROOT_SCALED", scaled)
        emit_physical(emitter, walker, root_prefix + "Q5_ROOT_ORIGINAL", original)
        emit_physical(emitter, walker, root_prefix + "Q5_SCALE_RATIO", ratio)
        emitter.emit(root_prefix + "Q5_SCALE_RATIO_DOMAIN_TEST", ratio_defined)
        emitter.emit(root_prefix + "Q5_SCALE_EXPONENT", exponent)

        quotient = assumed_simplify(root / k_squared, assumptions)
        emit_physical(emitter, walker, root_prefix + "Q6_ROOT_OVER_WAVENUMBER_NORM", quotient)

    emitter.emit(prefix + "Q6_ENERGY_DENSITY_DIMENSION", dim_expr(ENERGY_DENSITY_DIM))
    emitter.emit(prefix + "Q6_DIMENSION_EQUATIONS", sp.Tuple(*dimension_equations))
    emitter.emit(prefix + "Q6_DIMENSION_SOLUTION", dimension_solution)
    emitter.emit(
        prefix + "Q6_INDEPENDENT_DIMENSION_EQUATION_COUNT",
        sp.Integer(independent_dimension_equation_count),
    )
    emitter.emit(
        prefix + "Q6_UNKNOWN_COEFFICIENT_DIMENSION_COUNT",
        sp.Integer(len(dimension_unknowns)),
    )
    emitter.emit(prefix + "Q6_DIMENSION_COUNT_DIFFERENCE", sp.Integer(dimension_count_difference))
    emitter.emit(prefix + "Q6_DIMENSION_DETERMINATION", dimension_determination)
    emitter.emit(
        prefix + "Q6_SOLVED_HOMOGENEITY_VACUOUS",
        sp.true if dimension_count_difference <= 0 else sp.false,
    )
    emitter.emit(
        prefix + "Q6_SOLVED_HOMOGENEITY_VACUOUS_OPERANDS",
        sp.Tuple(
            sp.Integer(independent_dimension_equation_count),
            sp.Integer(len(dimension_unknowns)),
            sp.Integer(dimension_count_difference),
        ),
    )
    emitter.emit(
        prefix + "Q6_DIMENSION_SOLVE_CONDITION",
        Str("solution_returned") if dimension_solution else Str("no_solution_returned"),  # DERIVED · dimension-solver status
    )
    emitter.emit(
        prefix + "Q6_DIMENSION_SOLVE_CONDITION_OPERANDS",
        sp.Tuple(sp.Tuple(*dimension_equations), sp.Tuple(*dimension_unknowns), dimension_consistent),
    )
    emitter.emit(prefix + "Q6_INERTIAL_COEFFICIENTS", sp.Tuple(*inertial_coefficients))
    for index, coefficient in enumerate(inertial_coefficients, 1):
        emit_physical(
            emitter,
            walker,
            prefix + f"Q6_INERTIAL_COEFFICIENT{index}",
            coefficient,
            homogeneity_class="SOLVED",
        )
    emitter.emit(prefix + "Q6_STIFFNESS_COEFFICIENTS", sp.Tuple(*stiffness_coefficients))
    for index, coefficient in enumerate(stiffness_coefficients, 1):
        emit_physical(
            emitter,
            walker,
            prefix + f"Q6_STIFFNESS_COEFFICIENT{index}",
            coefficient,
            homogeneity_class="SOLVED",
        )

    gradient_symbols = [[sp.Symbol(f"g{i + 1}{j + 1}", real=True) for j in range(3)] for i in range(3)]  # COORDINATE · auxiliary gradient-component coordinates
    if n == 3:
        gradient_substitution = {
            sp.diff(fields[j], xvec[i]): gradient_symbols[i][j]
            for i in range(3)
            for j in range(3)
        }
        package_stiffness = sp.expand(stiffness.xreplace(gradient_substitution))
        curl_vector = sp.Matrix(
            [
                gradient_symbols[1][2] - gradient_symbols[2][1],
                gradient_symbols[2][0] - gradient_symbols[0][2],
                gradient_symbols[0][1] - gradient_symbols[1][0],
            ]
        )
        curl_dot = sp.expand((curl_vector.T * curl_vector)[0])
        curl_difference = sp.expand(package_stiffness - curl_dot)
        gradient_dim_map = {symbol: ZERO_DIM for row in gradient_symbols for symbol in row}
        q7_walker = DimensionWalker(gradient_dim_map, (), {}, assumptions)
        emit_physical(emitter, q7_walker, prefix + "Q7_STIFFNESS", package_stiffness)
        emit_physical(emitter, q7_walker, prefix + "Q7_CURL_DOT", curl_dot)
        emit_physical(emitter, q7_walker, prefix + "Q7_DIFFERENCE", curl_difference)
    else:
        emitter.emit(prefix + "Q7_OBJECTS", sp.Tuple())

    stratum_candidates: dict[str, tuple[sp.Tuple, object, dict[sp.Symbol, sp.Expr] | None, Str]] = {}

    def enroll_locus(locus: RealLocus, allowed_data: LocusAllowedData) -> None:
        for branch_index in range(len(locus.branches)):
            branch_expression = locus.branch_expression(branch_index)
            candidate = (
                branch_expression,
                allowed_data.branch_tests[branch_index],
                allowed_data.branch_points[branch_index],
                allowed_data.branch_reasons[branch_index],
            )
            key = str(branch_expression)
            existing = stratum_candidates.get(key)
            if existing is None or (
                candidate[1] is sp.true and existing[1] is not sp.true
            ):
                stratum_candidates[key] = candidate

    coincidence_loci = sp.Tuple(*(item[2].expression() for item in coincidences))
    coincidence_operands = sp.Tuple(*(item[3].operands for item in coincidences))
    coincidence_tests = sp.Tuple(*(item[3].allowed for item in coincidences))
    for item in coincidences:
        enroll_locus(item[2], item[3])

    for index, (root, objects) in enumerate(zip(roots, q4_by_root), 1):
        root_prefix = prefix + f"ROOT{index}_"
        root_matrix = objects["N1_MATRIX"]
        generic_rank = int(objects["N2_RANK"])
        if not isinstance(root_matrix, sp.MatrixBase):
            raise RuntimeError("Q8 root matrix is not a matrix")
        minors = rank_drop_minors(sp.Matrix(root_matrix), generic_rank)
        emit_physical(emitter, walker, root_prefix + "Q8_RANK_DROP_MINORS", sp.Tuple(*minors))
        rank_locus = solve_real_locus(minors, kvec, assumptions, impossible=(generic_rank == 0))
        rank_allowed = locus_allowed_data(rank_locus, kvec, k_squared, assumptions)
        emitter.emit(root_prefix + "Q8_RANK_DROP_LOCUS", rank_locus.expression())
        emitter.emit(root_prefix + "Q8_RANK_DROP_REALITY_FILTER", rank_locus.reality_filter)
        emitter.emit(root_prefix + "Q8_RANK_DROP_ALLOWED_OPERANDS", rank_allowed.operands)
        emitter.emit(root_prefix + "Q8_RANK_DROP_ALLOWED_TEST", rank_allowed.allowed)
        emitter.emit(root_prefix + "Q8_RANK_DROP_ALLOWED_WITNESSES", rank_allowed.witnesses)
        emitter.emit(
            root_prefix + "Q8_RANK_DROP_BRANCH_ALLOWED_OPERANDS",
            rank_allowed.branch_operands,
        )
        emitter.emit(
            root_prefix + "Q8_RANK_DROP_BRANCH_ALLOWED_TESTS",
            rank_allowed.branch_tests,
        )
        emitter.emit(
            root_prefix + "Q8_RANK_DROP_BRANCH_ALLOWED_WITNESSES",
            rank_allowed.branch_witnesses,
        )
        emitter.emit(root_prefix + "Q8_ROOT_COINCIDENCE_LOCI", coincidence_loci)
        emitter.emit(root_prefix + "Q8_ROOT_COINCIDENCE_ALLOWED_OPERANDS", coincidence_operands)
        emitter.emit(root_prefix + "Q8_ROOT_COINCIDENCE_ALLOWED_TESTS", coincidence_tests)
        local_prefix = "PY_S10_LOCAL_" + root_prefix.removeprefix("PY_S10_")
        emitter.emit(local_prefix + "Q8_RANK_DROP_SOLVER_OUTPUT", rank_locus.raw_solver_output)
        emitter.emit(local_prefix + "Q8_RANK_DROP_SOLVER_ROUTE", rank_locus.solver_route)
        enroll_locus(rank_locus, rank_allowed)

    allowed_strata = [
        candidate
        for candidate in stratum_candidates.values()
        if candidate[1] is sp.true and candidate[2] is not None
    ]
    skipped_strata = [
        candidate
        for candidate in stratum_candidates.values()
        if candidate[1] is not sp.true or candidate[2] is None
    ]
    emitter.emit(
        prefix + "Q8_ALLOWED_STRATA",
        sp.Tuple(*(candidate[0] for candidate in allowed_strata)),
    )
    emitter.emit(
        prefix + "Q8_SKIPPED_STRATA",
        sp.Tuple(
            *(
                sp.Tuple(candidate[0], candidate[1], candidate[3])
                for candidate in skipped_strata
            )
        ),
    )
    for skipped_index, candidate in enumerate(skipped_strata, 1):
        skipped_prefix = prefix + f"Q8_SKIPPED_STRATUM{skipped_index}_"
        emitter.emit(skipped_prefix + "BRANCH", candidate[0])
        emitter.emit(skipped_prefix + "ALLOWED_TEST", candidate[1])
        emitter.emit(skipped_prefix + "REASON", candidate[3])

    for stratum_index, candidate in enumerate(allowed_strata, 1):
        point = candidate[2]
        if point is None:
            raise RuntimeError("allowed stratum has no explicit point")
        stratum_prefix = prefix + f"Q8_STRATUM{stratum_index}_"
        emitter.emit(stratum_prefix + "SKIP_STATUS", Str("not_skipped_allowed_branch"))  # DERIVED · stratum-processing status
        emitter.emit(
            stratum_prefix + "POINT",
            sp.Tuple(*(sp.Eq(component, point[component], evaluate=False) for component in kvec)),
        )
        emit_stratum_q3_q4(emitter, walker, stratum_prefix, matrix_b, point, kvec, assumptions)

    register_declared_symbols(locals())
    return solved_coefficient_dims


@dataclass(frozen=True)
class ExportRecord:
    name: str
    value: object
    class_tag: str
    computation_dimension: sp.Expr
    dimension: object | None = None
    overwrites_upstream: bool = False


def output_class(suffix: str) -> str:
    return POSITED_OUTPUT_CLASSES.get(suffix, "DERIVED")


def stable_indexed_suffix(suffix: str) -> tuple[str, tuple[sp.Integer, ...]]:
    indices: list[sp.Integer] = []

    def replace(match: re.Match[str]) -> str:
        indices.append(sp.Integer(match.group("index")))
        return match.group("label") + "S"

    return INDEXED_EXPORT_TOKEN.sub(replace, suffix), tuple(indices)


def collect_main_export_records(
    emitter: Emitter,
    component_counts: Sequence[int],
    derived_dimensions: Mapping[sp.Symbol, Sequence[sp.Expr]],
) -> list[ExportRecord]:
    records = [
        ExportRecord(
            "inertia_coefficient_dimension",
            sp.Matrix(derived_dimensions[rho_br]),
            "DERIVED",
            D,
            overwrites_upstream=True,
        ),
        ExportRecord(
            "stiffness_coefficient_dimension",
            sp.Matrix(derived_dimensions[mu_R]),
            "DERIVED",
            D,
            overwrites_upstream=True,
        ),
        ExportRecord(
            "coefficient_dimension_difference",
            sp.simplify(
                sp.Matrix(derived_dimensions[mu_R])
                - sp.Matrix(derived_dimensions[rho_br])
            ),
            "DERIVED",
            D,
            overwrites_upstream=True,
        ),
    ]
    symbolic_slots: dict[str, list[sp.Tuple]] = {
        suffix: [] for suffix in SYMBOLIC_D_EXPORT_SUFFIXES
    }
    for component_count in component_counts:
        prefix = f"PY_S10_MAIN_D{component_count}_"
        indexed_rows: list[sp.Tuple] = []
        for tag, payload in emitter.values.items():
            if not tag.startswith(prefix):
                continue
            suffix = tag.removeprefix(prefix)
            stable_suffix, indices = stable_indexed_suffix(suffix)
            dimension = emitter.dimensions.get(tag)
            if indices:
                row_parts: list[object] = [
                    Str(stable_suffix.lower()),  # DERIVED · authored indexed-object field name
                    sp.Tuple(*indices),
                    payload,
                ]
                if dimension is not None:
                    row_parts.append(dimension)
                indexed_rows.append(sp.Tuple(*row_parts))
                continue
            if suffix in SYMBOLIC_D_EXPORT_SUFFIXES:
                slot_parts: list[object] = [sp.Integer(component_count), payload]
                if dimension is not None:
                    slot_parts.append(dimension)
                symbolic_slots[suffix].append(sp.Tuple(*slot_parts))
                continue
            records.append(
                ExportRecord(
                    suffix.lower(),
                    payload,
                    output_class(suffix),
                    sp.Integer(component_count),
                    dimension,
                )
            )
        records.append(
            ExportRecord(
                "indexed_derivations",
                sp.Tuple(*indexed_rows),
                "DERIVED",
                sp.Integer(component_count),
            )
        )
    records.extend(
        ExportRecord(
            suffix.lower(),
            sp.Tuple(*slots),
            output_class(suffix),
            D,
        )
        for suffix, slots in symbolic_slots.items()
    )
    return records


def reconstruction_expression(value: object) -> str:
    return sp.srepr(value)


def validate_export_tree(value: object) -> None:
    if isinstance(value, sp.MutableMatrix):
        raise TypeError("export tree contains a mutable matrix")
    if isinstance(value, sp.MatrixBase):
        for item in value:
            validate_export_tree(item)
        return
    if isinstance(value, sp.Basic):
        for argument in value.args:
            if not isinstance(argument, (sp.Basic, sp.MatrixBase)):
                raise TypeError("SymPy export object contains a raw Python argument")
            validate_export_tree(argument)
        return
    raise TypeError("export tree contains a non-SymPy object")


def traversable_export_value(value: object) -> object:
    if isinstance(value, sp.MatrixBase):
        normalized = sp.ImmutableMatrix(
            value.rows,
            value.cols,
            [traversable_export_value(item) for item in value],
        )
        validate_export_tree(normalized)
        return normalized
    if isinstance(value, sp.Basic):
        validate_export_tree(value)
        return value
    if isinstance(value, dict):
        normalized = sp.Dict({
            traversable_export_value(key): traversable_export_value(item)
            for key, item in value.items()
        })
    elif isinstance(value, (list, tuple)):
        normalized = sp.Tuple(*(traversable_export_value(item) for item in value))
    elif isinstance(value, (set, frozenset)):
        normalized = sp.FiniteSet(*(traversable_export_value(item) for item in value))
    elif isinstance(value, str):
        normalized = Str(value)
    else:
        normalized = sp.sympify(value)
    validate_export_tree(normalized)
    return normalized


def exact_reconstruction_match(live_value: object, reconstructed_value: object) -> bool:
    return (
        type(live_value) is type(reconstructed_value)
        and sp.srepr(live_value) == sp.srepr(reconstructed_value)
    )


def exact_value_residual(left: object, right: object) -> sp.Integer:
    return sp.Integer(not exact_reconstruction_match(left, right))


def export_value_kind(value: object) -> str:
    return "AUTHORED_WORD" if isinstance(value, Str) else "COMPUTED_OBJECT"


def generated_record_lines(
    name: str,
    value: object,
    class_tag: str,
    step: str,
    value_kind: str | None = None,
    dimension_key: str | None = None,
    corroborated_steps: tuple[str, ...] | None = None,
    display: str | None = None,
) -> list[str]:
    lines = [
        f"    {name!r}: {{",
        f"        'display': {(sp.sstr(value) if display is None else display)!r},",
        f"        'value': _restore({reconstruction_expression(value)!r}),",
        f"        'value_kind': {(export_value_kind(value) if value_kind is None else value_kind)!r},",
    ]
    if dimension_key is not None:
        lines.append(f"        'dimension_key': {dimension_key!r},")
    if corroborated_steps is not None:
        lines.append(f"        'corroborated_steps': {corroborated_steps!r},")
    lines.extend([
        f"        'class': {class_tag!r},",
        f"        'step': {step!r},",
        "    },",
    ])
    return lines


def generated_ledger_key(name: str, computation_dimension: sp.Expr) -> str:
    base_name = re.sub(r"_d[0-9]+$", "", name)
    if computation_dimension == D:
        return base_name
    return f"{base_name}_d{computation_dimension}"


def d_partition(records: Sequence[tuple[str, sp.Expr]]) -> sp.Tuple:
    grouped: dict[sp.Expr, list[Str]] = {}
    for name, computation_dimension in records:
        grouped.setdefault(computation_dimension, []).append(Str(name))  # DERIVED · production-partition member name
    return sp.Tuple(*(
        sp.Tuple(
            computation_dimension,
            sp.Tuple(*sorted(names, key=sp.default_sort_key)),
        )
        for computation_dimension, names in sorted(
            grouped.items(), key=lambda item: sp.default_sort_key(item[0])
        )
    ))


def dimension_record_key(
    owner_name: str,
    dimension: object,
    records: Mapping[str, Mapping[str, object]],
) -> str:
    normalized_dimension = traversable_export_value(dimension)
    candidates = [
        name
        for name, record in records.items()
        if exact_reconstruction_match(normalized_dimension, record["value"])
    ]
    if not candidates:
        raise RuntimeError("export dimension has no ledger record")
    dimension_named = [
        name for name in candidates
        if name.startswith("dim_") or "dimension" in name
    ]
    if dimension_named:
        candidates = dimension_named

    def rank(name: str) -> tuple[int, int, str]:
        common_prefix_length = len(os.path.commonprefix((owner_name, name)))
        return (
            -common_prefix_length,
            0 if name.endswith("_coefficient_dimension") else
            1 if name.startswith("dim_") else
            2 if "_dimension" in name else
            3,
            name,
        )

    return min(candidates, key=rank)


def add_symbol_records(merged: dict[str, dict[str, object]]) -> None:
    occurrence_classes: dict[str, set[str]] = {}
    occurrence_objects: dict[str, dict[str, sp.Symbol]] = {}
    for record in merged.values():
        value = record["value"]
        if not isinstance(value, (sp.Basic, sp.MatrixBase)):
            raise TypeError("export value is not a traversable SymPy object")
        for symbol in value.atoms(sp.Symbol):
            occurrence_classes.setdefault(symbol.name, set()).add(str(record["class"]))
            occurrence_objects.setdefault(symbol.name, {})[sp.srepr(symbol)] = symbol
    for symbol_name in sorted(occurrence_objects):
        if symbol_name in merged:
            continue
        classes = occurrence_classes[symbol_name]
        class_tag = declared_symbol_classes.get(
            symbol_name,
            next(iter(classes)) if len(classes) == 1 else "DERIVED",
        )
        symbol = occurrence_objects[symbol_name][sorted(occurrence_objects[symbol_name])[0]]
        merged[symbol_name] = {
            "display": sp.sstr(symbol),
            "value": symbol,
            "value_kind": export_value_kind(symbol),
            "class": class_tag,
            "step": "S10",
        }


def symbol_binding_operands_and_residual(
    ledger: Mapping[str, Mapping[str, object]],
) -> tuple[sp.Tuple, sp.Integer]:
    symbols_by_name: dict[str, dict[str, sp.Symbol]] = {}
    for record in ledger.values():
        value = record["value"]
        if not isinstance(value, (sp.Basic, sp.MatrixBase)):
            raise TypeError("written export value is not a traversable SymPy object")
        for symbol in value.atoms(sp.Symbol):
            symbols_by_name.setdefault(symbol.name, {})[sp.srepr(symbol)] = symbol
    rows: list[sp.Tuple] = []
    residual = sp.Integer(0)
    for symbol_name in sorted(symbols_by_name):
        variants = symbols_by_name[symbol_name]
        binding = ledger.get(symbol_name, {}).get(
            "value", Str("missing_ledger_binding")
        )
        rows.append(
            sp.Tuple(
                Str(symbol_name),
                binding,
                sp.Tuple(*(variants[source] for source in sorted(variants))),
            )
        )
        residual += sp.Integer(max(len(variants) - 1, 0))
        residual += sp.Integer(
            not isinstance(binding, sp.Symbol)
            or binding.name != symbol_name
            or sp.srepr(binding) not in variants
        )
    return sp.Tuple(*rows), residual


def write_exports(emitter: Emitter, own_records: Sequence[ExportRecord]) -> None:
    register_declared_symbols(globals())
    keyed_own_records: list[tuple[str, ExportRecord]] = []
    seen_own_keys: set[str] = set()
    for record in own_records:
        name = generated_ledger_key(record.name, record.computation_dimension)
        if name in seen_own_keys:
            raise RuntimeError(f"duplicate S10 export key: {name}")
        seen_own_keys.add(name)
        keyed_own_records.append((name, record))

    overwrite_rows: list[sp.Tuple] = []
    overwrite_residuals: list[sp.Integer] = []
    for name, record in keyed_own_records:
        if name not in S9_LEDGER:
            continue
        upstream = S9_LEDGER[name]
        downstream_value = traversable_export_value(record.value)
        value_residual = exact_value_residual(upstream["value"], downstream_value)
        class_residual = sp.Integer(upstream["class"] != record.class_tag)
        provenance_residual = sp.Integer(
            not record.overwrites_upstream or upstream["step"] != "S9"
        )
        overwrite_rows.append(
            sp.Tuple(
                Str(name),  # DERIVED · authored overwritten-object name
                upstream["value"],
                downstream_value,
                value_residual,
                Str(str(upstream["class"])),
                Str(record.class_tag),
                class_residual,
                Str(str(upstream["step"])),
                Str("S10"),
                sp.true if record.overwrites_upstream else sp.false,
                provenance_residual,
            )
        )
        overwrite_residuals.extend(
            (value_residual, class_residual, provenance_residual)
        )
    emitter.emit(
        "PY_S10_EXPORT_OVERWRITE_OPERANDS_AND_RESIDUALS",
        sp.Tuple(*overwrite_rows),
    )
    assert all(residual == 0 for residual in overwrite_residuals)

    merged: dict[str, dict[str, object]] = {
        name: dict(record) for name, record in S9_LEDGER.items()
    }
    for name, record in keyed_own_records:
        stored_value = traversable_export_value(record.value)
        generated = {
            "display": sp.sstr(stored_value),
            "value": stored_value,
            "value_kind": export_value_kind(stored_value),
            "class": record.class_tag,
            "step": "S10",
        }
        if name in S9_LEDGER and record.overwrites_upstream:
            generated["corroborated_steps"] = ("S9", "S10")
        merged[name] = generated
    add_symbol_records(merged)
    for name, record in keyed_own_records:
        if record.dimension is not None:
            merged[name]["dimension_key"] = dimension_record_key(
                name,
                record.dimension,
                merged,
            )

    build_input_digests = {
        Path(__file__).name: file_sha256(Path(__file__)),
        "extract_knob_inventory.py": file_sha256(
            Path(__file__).with_name("extract_knob_inventory.py")
        ),
        "S9_exports.py": file_sha256(Path(__file__).with_name("S9_exports.py")),
    }
    source_lines = [
        "# S10_exports.py — GENERATED by S10_brane_mode_spectrum_sympy_audit.py. Do not edit.",
        "from types import MappingProxyType",
        "",
        "import sympy as sp",
        "from sympy.core.symbol import Str",
        "",
        f"BUILD_INPUT_DIGESTS = MappingProxyType({build_input_digests!r})",
        "",
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
        "",
        "def _restore(source):",
        "    return eval(source, {'__builtins__': {}, 'Str': Str, **vars(sp), **_RELATIONALS})",
        "",
        "",
        "_LEDGER = {",
    ]
    for name, record in merged.items():
        source_lines.extend(
            generated_record_lines(
                name,
                record["value"],
                str(record["class"]),
                str(record["step"]),
                str(record["value_kind"]),
                record.get("dimension_key"),
                record.get("corroborated_steps"),
                str(record["display"]),
            )
        )
    source_lines.extend((
        "}",
        "LEDGER = MappingProxyType({",
        "    name: MappingProxyType(record) for name, record in _LEDGER.items()",
        "})",
        "del _LEDGER",
        "",
    ))
    export_source = "\n".join(source_lines)
    export_path = EXPORT_PATH

    reconstructed_namespace: dict[str, object] = {}
    exec(
        compile(export_source, str(export_path), "exec"),
        reconstructed_namespace,
    )
    reconstructed_ledger = reconstructed_namespace["LEDGER"]
    reconstructed_build_inputs = reconstructed_namespace["BUILD_INPUT_DIGESTS"]
    if not isinstance(reconstructed_ledger, Mapping):
        raise RuntimeError("generated S10 ledger did not reconstruct as a mapping")
    live_count = len(merged)
    residuals: list[sp.Integer] = []
    for name, record in merged.items():
        reconstructed_record = reconstructed_ledger[name]
        residuals.append(
            exact_value_residual(record["value"], reconstructed_record["value"])
        )
        residuals.append(sp.Integer(
            record.get("dimension_key") != reconstructed_record.get("dimension_key")
        ))
        residuals.append(sp.Integer(
            record.get("corroborated_steps")
            != reconstructed_record.get("corroborated_steps")
        ))
        residuals.append(sp.Integer(record["class"] != reconstructed_record["class"]))
        residuals.append(sp.Integer(record["step"] != reconstructed_record["step"]))
    roundtrip_residual = sum(residuals, sp.Integer(0))
    roundtrip_count_residual = sp.Integer(live_count - len(reconstructed_ledger))
    class_tally = sp.Tuple(*(
        sp.Tuple(
            Str(class_tag),  # STRUCTURAL · exported class label
            sp.Integer(
                sum(
                    record["class"] == class_tag
                    for record in reconstructed_ledger.values()
                )
            ),
        )
        for class_tag in CLASS_TAGS
    ))
    class_tally_residual = sp.Integer(len(reconstructed_ledger)) - sum(
        (count for _, count in class_tally), sp.Integer(0)
    )
    production_d_partition = d_partition([
        (name, record.computation_dimension) for name, record in keyed_own_records
    ])
    symbol_binding_operands, symbol_binding_residual = (
        symbol_binding_operands_and_residual(reconstructed_ledger)
    )
    assumption_channel_operands, assumption_channel_residual = (
        assumption_channel_operands_and_residual(
            reconstructed_ledger,
            [
                value
                for name, value in emitter.values.items()
                if name.startswith("PY_S10_MAIN_")
                and name.endswith("_ASSUMPTION_JOINT_PREDICATE")
            ],
        )
    )
    authored_word_records = sp.Tuple(*(
        Str(name)
        for name, record in reconstructed_ledger.items()
        if isinstance(record["value"], Str)
    ))
    dimension_links = sp.Tuple(*(
        sp.Tuple(Str(name), Str(str(record["dimension_key"])))
        for name, record in reconstructed_ledger.items()
        if "dimension_key" in record
    ))
    mapping_immutability_operands = sp.Tuple(
        sp.true if isinstance(reconstructed_ledger, MappingProxyType) else sp.false,
        sp.Tuple(*(
            sp.true if isinstance(record, MappingProxyType) else sp.false
            for record in reconstructed_ledger.values()
        )),
    )
    mapping_immutability_residual = sp.Integer(
        not isinstance(reconstructed_ledger, MappingProxyType)
    ) + sum(
        (
            sp.Integer(not isinstance(record, MappingProxyType))
            for record in reconstructed_ledger.values()
        ),
        sp.Integer(0),
    )
    build_input_operands = sp.Tuple(*(
        sp.Tuple(Str(name), Str(digest))
        for name, digest in reconstructed_build_inputs.items()
    ))
    emitter.emit("PY_S10_EXPORT_ROUNDTRIP_LIVE_COUNT", sp.Integer(live_count))
    emitter.emit(
        "PY_S10_EXPORT_ROUNDTRIP_RECONSTRUCTED_COUNT",
        sp.Integer(len(reconstructed_ledger)),
    )
    emitter.emit(
        "PY_S10_EXPORT_ROUNDTRIP_COUNT_RESIDUAL", roundtrip_count_residual
    )
    emitter.emit("PY_S10_EXPORT_ROUNDTRIP_RESIDUAL", roundtrip_residual)
    emitter.emit("PY_S10_EXPORT_CLASS_TALLY", class_tally)
    emitter.emit("PY_S10_EXPORT_CLASS_TALLY_RESIDUAL", class_tally_residual)
    emitter.emit("PY_S10_EXPORT_D_PARTITION", production_d_partition)
    emitter.emit(
        "PY_S10_EXPORT_SYMBOL_BINDING_OPERANDS", symbol_binding_operands
    )
    emitter.emit("PY_S10_EXPORT_SYMBOL_BINDING_RESIDUAL", symbol_binding_residual)
    emitter.emit(
        "PY_S10_EXPORT_ASSUMPTION_CHANNEL_OPERANDS", assumption_channel_operands
    )
    emitter.emit(
        "PY_S10_EXPORT_ASSUMPTION_CHANNEL_RESIDUAL", assumption_channel_residual
    )
    emitter.emit("PY_S10_EXPORT_AUTHORED_WORD_RECORDS", authored_word_records)
    emitter.emit(
        "PY_S10_EXPORT_AUTHORED_WORD_RECORD_COUNT",
        sp.Integer(len(authored_word_records)),
    )
    emitter.emit("PY_S10_EXPORT_DIMENSION_LINKS", dimension_links)
    emitter.emit(
        "PY_S10_EXPORT_DIMENSION_LINK_COUNT", sp.Integer(len(dimension_links))
    )
    emitter.emit(
        "PY_S10_EXPORT_MAPPING_IMMUTABILITY_OPERANDS",
        mapping_immutability_operands,
    )
    emitter.emit(
        "PY_S10_EXPORT_MAPPING_IMMUTABILITY_RESIDUAL",
        mapping_immutability_residual,
    )
    emitter.emit("PY_S10_EXPORT_BUILD_INPUT_DIGESTS", build_input_operands)
    assert roundtrip_count_residual == 0
    assert roundtrip_residual == 0
    assert class_tally_residual == 0
    assert symbol_binding_residual == 0
    assert assumption_channel_residual == 0
    assert mapping_immutability_residual == 0
    export_path.write_text(export_source, encoding="utf-8")


def main() -> int:
    emitter = Emitter()
    declared_pairs = [
        (package.name, n) for package in PACKAGES for n in package.dimensions
    ]
    completed_pairs: list[tuple[str, int]] = []
    derived_dimensions: dict[sp.Symbol, tuple[sp.Expr, ...]] | None = None
    for package in PACKAGES:
        for n in package.dimensions:
            current = run_package_dimension(emitter, package, n)
            completed_pairs.append((package.name, n))
            if derived_dimensions is None:
                derived_dimensions = current
    if derived_dimensions is None:
        raise RuntimeError("no package/dimension pairs were run")
    completed_set = set(completed_pairs)
    skipped_pairs = [pair for pair in declared_pairs if pair not in completed_set]
    emitter.emit(
        "PY_S10_RUN_PAIRS",
        sp.Tuple(*(Str(f"{name}_D{n}") for name, n in completed_pairs)),  # DERIVED · completed construction-pair collection
    )
    emitter.emit(
        "PY_S10_SKIPPED_PAIRS",
        sp.Tuple(*(Str(f"{name}_D{n}") for name, n in skipped_pairs)),  # DERIVED · skipped construction-pair collection
    )
    local_list_tag = "PY_S10_LOCAL_TAG_NAMES"
    local_names = [*emitter.local_names, local_list_tag]
    emitter.emit(local_list_tag, sp.Tuple(*(Str(name) for name in local_names)))  # DERIVED · local emission-name inventory
    main_component_counts = next(
        package.dimensions for package in PACKAGES if package.name == "MAIN"
    )
    own_records = collect_main_export_records(
        emitter, main_component_counts, derived_dimensions
    )
    write_exports(emitter, own_records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
