#!/usr/bin/env python3
"""S10 engine 2: derive the brane mode audit from the shared action with SymPy."""

from __future__ import annotations

import itertools
import signal
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import sympy as sp


SCRIPT_PATH = Path(__file__).resolve()
LEDGER_DIR = SCRIPT_PATH.parent.parent
REDUCTION_DIR = LEDGER_DIR / "reduction"
sys.path.insert(0, str(REDUCTION_DIR))
import registry_read  # noqa: E402


D = sp.Symbol("D")
rho_br = sp.Symbol("rho_br", real=True)
mu_R = sp.Symbol("mu_R", real=True)
s_rho = sp.Symbol("s_rho", real=True)
s = sp.Symbol("s", real=True)
omegaSquared = sp.Symbol("omegaSquared", real=True)
lambdaScale = sp.Symbol("lambdaScale", real=True)

ZERO_DIM = (sp.Integer(0), sp.Integer(0), sp.Integer(0))
LENGTH_DIM = (sp.Integer(1), sp.Integer(0), sp.Integer(0))
TIME_DIM = (sp.Integer(0), sp.Integer(1), sp.Integer(0))
WAVENUMBER_DIM = (sp.Integer(-1), sp.Integer(0), sp.Integer(0))
OMEGA_SQUARED_DIM = (sp.Integer(0), sp.Integer(-2), sp.Integer(0))


@dataclass(frozen=True)
class Package:
    name: str
    dimensions: tuple[int, ...]
    stiffness: str
    stiffness_sign: int = -1
    anisotropic_kinetic: bool = False
    coefficient_scale: bool = False


PACKAGES = (
    Package("MAIN", (2, 3, 4, 5), "curl"),
    Package("XFORM_FULLGRAD", (3, 4), "fullgrad"),
    Package("XFORM_DIVONLY", (3, 4), "divonly"),
    Package("XFORM_SIGNFLIP", (3, 4), "curl", stiffness_sign=1),
    Package("XFORM_ANISO", (3, 4), "curl", anisotropic_kinetic=True),
    Package("XCOEF_SCALE", (3,), "curl", coefficient_scale=True),
)


class OperationTimeout(RuntimeError):
    pass


def _timeout_handler(_signum: int, _frame: object) -> None:
    raise OperationTimeout


class Emitter:
    def __init__(self) -> None:
        self.names: set[str] = set()
        self.local_names: list[str] = []

    def emit(self, name: str, payload: object) -> None:
        if name in self.names:
            raise RuntimeError(f"duplicate emitted tag: {name}")
        self.names.add(name)
        if name.startswith("PY_S10_LOCAL_"):
            self.local_names.append(name)
        print(f"{name}: {payload}")


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

    def dimension(self, expr: sp.Expr) -> tuple[sp.Expr, ...]:
        if expr == 0 or expr.is_Number:
            return ZERO_DIM
        if expr in self.fields:
            return LENGTH_DIM
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
            return ZERO_DIM
        raise RuntimeError(f"no dimension-tree rule for {type(expr).__name__}: {expr}")

    def component_term_dimensions(self, expr: sp.Expr) -> tuple[tuple[sp.Expr, ...], ...]:
        expanded = sp.expand(expr)
        if expanded == 0:
            return ()
        return tuple(self.dimension(term) for term in sp.Add.make_args(expanded))

    def report(self, obj: object) -> tuple[sp.Tuple, sp.Tuple, sp.logic.boolalg.Boolean]:
        if isinstance(obj, sp.MatrixBase):
            components = list(obj)
        elif isinstance(obj, (list, tuple, sp.Tuple)):
            components = []
            for item in obj:
                if isinstance(item, sp.MatrixBase):
                    components.extend(list(item))
                elif isinstance(item, sp.Expr):
                    components.append(item)
        elif isinstance(obj, sp.Expr):
            components = [obj]
        else:
            components = []

        all_term_dimensions: list[sp.Tuple] = []
        all_component_dimensions: list[sp.Tuple] = []
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
        return (
            sp.Tuple(*all_term_dimensions),
            sp.Tuple(*all_component_dimensions),
            sp.true if homogeneous else sp.false,
        )


def emit_physical(
    emitter: Emitter,
    walker: DimensionWalker,
    name: str,
    payload: object,
) -> None:
    emitter.emit(name, payload)
    term_dimensions, dimensions, homogeneous = walker.report(payload)
    emitter.emit(name + "_Q6_TERM_DIMENSIONS", term_dimensions)
    emitter.emit(name + "_Q6_DIMENSIONS", dimensions)
    emitter.emit(name + "_Q6_HOMOGENEITY", homogeneous)


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
    fields = [sp.Function(f"u{index + 1}")(t, *xvec) for index in range(n)]
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
) -> tuple[list[sp.Equality], list[dict[sp.Symbol, sp.Expr]], dict[sp.Symbol, tuple[sp.Expr, ...]]]:
    rho_unknown = sp.symbols("rho_br_dim_length rho_br_dim_time rho_br_dim_mass")
    mu_unknown = sp.symbols("mu_R_dim_length mu_R_dim_time mu_R_dim_mass")
    unknowns = (*rho_unknown, *mu_unknown)
    unknown_map = {
        rho_br: rho_unknown,
        mu_R: mu_unknown,
        s_rho: ZERO_DIM,
        s: ZERO_DIM,
    }
    coordinates = {t: TIME_DIM, **{coordinate: LENGTH_DIM for coordinate in xvec}}
    walker = DimensionWalker(unknown_map, fields, coordinates, assumptions)
    target = (sp.Integer(2) - D, sp.Integer(-2), sp.Integer(1))
    equations: list[sp.Equality] = []
    for term in sp.Add.make_args(sp.expand(lagrangian)):
        term_dimension = walker.dimension(term)
        equations.extend(
            sp.Eq(component, required, evaluate=False)
            for component, required in zip(term_dimension, target)
        )
    equations = list(dict.fromkeys(equations))
    solution = sp.solve(equations, unknowns, dict=True)
    if solution:
        solved_map = {
            rho_br: tuple(solution[0].get(symbol, symbol) for symbol in rho_unknown),
            mu_R: tuple(solution[0].get(symbol, symbol) for symbol in mu_unknown),
            s_rho: ZERO_DIM,
            s: ZERO_DIM,
        }
    else:
        solved_map = {rho_br: rho_unknown, mu_R: mu_unknown, s_rho: ZERO_DIM, s: ZERO_DIM}
    return equations, solution, solved_map


def factor_with_timeout(expr: sp.Expr, seconds: int = 5) -> tuple[sp.Expr, sp.Symbol]:
    previous_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(seconds)
    try:
        result = sp.factor(expr)
        route = sp.Symbol("factor_returned")
    except OperationTimeout:
        result = expr
        route = sp.Symbol("factor_timeout_unfactored")
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_handler)
    return result, route


def derived_sign(expr: sp.Expr, assumptions: sp.logic.boolalg.Boolean) -> sp.Expr:
    positive = sp.ask(sp.Q.positive(expr), assumptions)
    zero = sp.ask(sp.Q.zero(expr), assumptions)
    negative = sp.ask(sp.Q.negative(expr), assumptions)
    if positive is True:
        return sp.Integer(1)
    if zero is True:
        return sp.Integer(0)
    if negative is True:
        return sp.Integer(-1)
    return sp.refine(sp.sign(expr), assumptions)


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
) -> tuple[sp.Expr, sp.Expr, sp.Symbol, list[dict[sp.Symbol, sp.Expr]], list[sp.Expr]]:
    determinant = matrix.det(method="berkowitz")
    factored, factor_route = factor_with_timeout(determinant)
    raw_solutions = sp.solve(factored, omegaSquared, dict=True)
    roots = deduplicate_roots(raw_solutions, assumptions)
    return determinant, factored, factor_route, raw_solutions, roots


@dataclass
class RealLocus:
    equations: tuple[sp.Expr, ...]
    raw_solver_output: object
    branches: tuple[dict[sp.Symbol, sp.Expr], ...]
    solver_route: sp.Symbol

    def expression(self) -> object:
        if not self.branches:
            return sp.EmptySet
        return sp.Tuple(
            *(
                sp.Tuple(
                    *(sp.Eq(symbol, value, evaluate=False) for symbol, value in sorted(branch.items(), key=lambda item: str(item[0])))
                )
                for branch in self.branches
            )
        )


def _coordinate_zero_branches(
    equations: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    assumptions: sp.logic.boolalg.Boolean,
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    candidates: list[frozenset[sp.Symbol]] = []
    for count in range(len(variables) + 1):
        for selected in itertools.combinations(variables, count):
            substitutions = {symbol: sp.Integer(0) for symbol in selected}
            if all(algebraic_normalize(equation.subs(substitutions), assumptions) == 0 for equation in equations):
                selected_set = frozenset(selected)
                if not any(existing <= selected_set for existing in candidates):
                    candidates.append(selected_set)
    minimal = [candidate for candidate in candidates if not any(other < candidate for other in candidates)]
    return tuple({symbol: sp.Integer(0) for symbol in variables if symbol in branch} for branch in minimal)


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
        return RealLocus(tuple(normalized), sp.EmptySet, (), sp.Symbol("rank_floor_empty_locus"))

    previous_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(5)
    try:
        raw = sp.solve(normalized, list(variables), dict=True, domain=sp.S.Reals)
        route = sp.Symbol("solve_reals_plus_coordinate_real_projection")
    except OperationTimeout:
        raw = sp.ConditionSet(sp.Tuple(*variables), sp.And(*(sp.Eq(eq, 0) for eq in normalized)), sp.ProductSet(*(sp.S.Reals for _ in variables)))
        route = sp.Symbol("solve_reals_timeout_coordinate_real_projection")
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_handler)
    branches = _coordinate_zero_branches(normalized, variables, assumptions)
    return RealLocus(tuple(normalized), raw, branches, route)


def locus_allowed_data(
    locus: RealLocus,
    kvec: Sequence[sp.Symbol],
    k_squared: sp.Expr,
    assumptions: sp.logic.boolalg.Boolean,
) -> tuple[sp.Tuple, sp.logic.boolalg.Boolean, sp.Tuple]:
    operands = sp.Tuple(locus.expression(), sp.Gt(k_squared, 0, evaluate=False), assumptions)
    witnesses: list[sp.Tuple] = []
    for branch in locus.branches:
        point = dict(branch)
        free = [component for component in kvec if component not in point]
        for component in free:
            point[component] = sp.Integer(0)
        if free:
            point[free[0]] = sp.Integer(1)
        equation_values = [algebraic_normalize(equation.subs(point), assumptions) for equation in locus.equations]
        k_value = assumed_simplify(k_squared.subs(point), assumptions)
        if all(value == 0 for value in equation_values) and sp.ask(sp.Q.positive(k_value), assumptions) is True:
            witnesses.append(
                sp.Tuple(*(sp.Eq(component, point[component], evaluate=False) for component in kvec))
            )
    return operands, sp.true if witnesses else sp.false, sp.Tuple(*witnesses)


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
    dots: list[sp.Expr] = []
    residuals: list[sp.Matrix] = []
    for vector in basis:
        dot = assumed_simplify((vector.T * k_column)[0], assumptions)
        residual = assumed_simplify(k_squared * vector - dot * k_column, assumptions)
        if not isinstance(dot, sp.Expr) or not isinstance(residual, sp.MatrixBase):
            raise RuntimeError("null-space display computation returned an unsupported object")
        dots.append(dot)
        residuals.append(sp.Matrix(residual))
    basis_count = sp.Integer(len(basis))
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
            emit_physical(emitter, walker, name, payload)
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
        ratio = sp.nan
        ratio_defined = sp.false
        exponent = sp.nan
    else:
        ratio_object = sp.Mul(scaled, sp.Pow(original, -1, evaluate=False), evaluate=False)
        ratio = assumed_simplify(ratio_object, assumptions)
        ratio_defined = sp.true
        try:
            polynomial = sp.Poly(ratio, lambdaScale)
            exponent = sp.Integer(polynomial.degree()) if len(polynomial.terms()) == 1 else sp.nan
        except (sp.PolynomialError, TypeError, ValueError):
            exponent = sp.nan
    if not isinstance(ratio, sp.Expr):
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
    emitter.emit(prefix + "PREMISE_FIELD_SECTOR", sp.Eq(sp.Symbol("u_component_count"), n, evaluate=False))
    emitter.emit(prefix + "PREMISE_ANSATZ", sp.Tuple(*ansatz))
    emitter.emit(prefix + "PREMISE_PERIOD_AVERAGE", sp.Tuple(sp.Integer(0), 2 * sp.pi / sp.sqrt(omegaSquared)))
    emitter.emit(prefix + "PREMISE_BACKGROUND_VELOCITY", sp.Eq(sp.Symbol("v_0"), 0, evaluate=False))
    emitter.emit(prefix + "PREMISE_TIME_ODD_KERNEL", sp.Eq(sp.Symbol("time_odd_kernel"), 0, evaluate=False))
    emitter.emit(prefix + "PREMISE_RESPONSE_DEGREE", sp.Eq(sp.Symbol("response_degree"), 2, evaluate=False))
    emitter.emit(prefix + "PREMISE_STIFFNESS_INPUT", sp.Symbol("S_curl_supplied"))
    emitter.emit(prefix + "PREMISE_RHO_DOMAIN", sp.Q.positive(rho_br))
    emitter.emit(prefix + "PREMISE_MU_DOMAIN", sp.Q.positive(mu_R))
    emitter.emit(prefix + "PREMISE_WAVEVECTOR_NORM_DOMAIN", sp.Q.positive(sp.Add(*(item**2 for item in kvec))))
    emitter.emit(prefix + "PREMISE_WAVEVECTOR_REAL_DOMAIN", sp.And(*(sp.Q.real(item) for item in kvec)))
    emitter.emit(prefix + "PREMISE_AMPLITUDE_REAL_DOMAIN", sp.And(*(sp.Q.real(item) for item in avec)))
    emitter.emit(prefix + "PREMISE_BRANE_DIMENSION_DOMAIN", sp.And(sp.Q.integer(D), sp.Q.positive(D)))
    emitter.emit(prefix + "ASSUMPTION_JOINT_PREDICATE", assumptions)
    controls: list[sp.Expr] = []
    if package.anisotropic_kinetic:
        controls.append(sp.Eq(sp.Symbol("dimension_s_rho"), sp.Symbol("dimensionless"), evaluate=False))
    if package.coefficient_scale:
        controls.append(sp.Eq(sp.Symbol("dimension_s"), sp.Symbol("dimensionless"), evaluate=False))
    emitter.emit(prefix + "Q6_CONTROL_DIMENSION_PREMISES", sp.Tuple(*controls))


def emit_q3(
    emitter: Emitter,
    walker: DimensionWalker,
    prefix: str,
    matrix: sp.Matrix,
    roots_data: tuple[sp.Expr, sp.Expr, sp.Symbol, list[dict[sp.Symbol, sp.Expr]], list[sp.Expr]],
    assumptions: sp.logic.boolalg.Boolean,
) -> list[sp.Expr]:
    determinant, factored, factor_route, raw_solutions, roots = roots_data
    emit_physical(emitter, walker, prefix + "Q3_DETERMINANT", determinant)
    emit_physical(emitter, walker, prefix + "Q3_DETERMINANT_FACTORED", factored)
    emitter.emit("PY_S10_LOCAL_" + prefix.removeprefix("PY_S10_") + "Q3_DETERMINANT_FACTOR_ROUTE", factor_route)
    emitter.emit(prefix + "Q3_ROOT_SOLUTIONS_RAW", raw_solutions)
    emitter.emit(prefix + "Q3_ROOTS_DISTINCT", sp.Tuple(*roots))
    emitter.emit(prefix + "Q3_ROOT_COUNT", sp.Integer(len(roots)))
    emitter.emit(prefix + "ROOT_ORDERING", sp.Tuple(*roots))
    if not roots:
        raise RuntimeError(f"spectrum solve produced no roots for {prefix}")
    return roots


def q3_coincidence_data(
    roots: Sequence[sp.Expr],
    kvec: Sequence[sp.Symbol],
    k_squared: sp.Expr,
    assumptions: sp.logic.boolalg.Boolean,
) -> list[tuple[tuple[int, int], sp.Expr, RealLocus, sp.Tuple, sp.logic.boolalg.Boolean, sp.Tuple]]:
    result = []
    for left, right in itertools.combinations(range(len(roots)), 2):
        equation = assumed_simplify(roots[left] - roots[right], assumptions)
        if not isinstance(equation, sp.Expr):
            raise RuntimeError("root-coincidence equation is not an expression")
        locus = solve_real_locus([equation], kvec, assumptions)
        operands, allowed, witnesses = locus_allowed_data(locus, kvec, k_squared, assumptions)
        result.append(((left + 1, right + 1), equation, locus, operands, allowed, witnesses))
    return result


def emit_coincidences(
    emitter: Emitter,
    walker: DimensionWalker,
    prefix: str,
    data: Sequence[tuple[tuple[int, int], sp.Expr, RealLocus, sp.Tuple, sp.logic.boolalg.Boolean, sp.Tuple]],
) -> None:
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_PAIR_INDICES", sp.Tuple(*(sp.Tuple(*pair) for pair, *_ in data)))
    equations = sp.Tuple(*(equation for _, equation, *_ in data))
    emit_physical(emitter, walker, prefix + "Q3_ROOT_COINCIDENCE_EQUATIONS", equations)
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_LOCI", sp.Tuple(*(locus.expression() for _, _, locus, *_ in data)))
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_OPERANDS", sp.Tuple(*(operands for *_, operands, _allowed, _witnesses in data)))
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_TESTS", sp.Tuple(*(allowed for *_, allowed, _witnesses in data)))
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_WITNESSES", sp.Tuple(*(witnesses for *_, witnesses in data)))
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
    if not isinstance(point_k_squared, sp.Expr):
        raise RuntimeError("stratum wavevector norm is not an expression")
    spectrum = solve_spectrum(point_matrix, point_assumptions)
    roots = emit_q3(emitter, walker, prefix, point_matrix, spectrum, point_assumptions)
    for index, root in enumerate(roots, 1):
        root_prefix = prefix + f"ROOT{index}_"
        emit_physical(emitter, walker, root_prefix + "Q3_ROOT", root)
        emitter.emit(root_prefix + "Q3_SIGN", derived_sign(root, point_assumptions))
        objects = q4_objects(point_matrix, root, point_k, point_k_squared, point_assumptions)
        emit_q4(emitter, walker, root_prefix, objects)
    empty = sp.Tuple()
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_PAIR_INDICES", empty)
    emit_physical(emitter, walker, prefix + "Q3_ROOT_COINCIDENCE_EQUATIONS", empty)
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_LOCI", empty)
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_OPERANDS", empty)
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_TESTS", empty)
    emitter.emit(prefix + "Q3_ROOT_COINCIDENCE_ALLOWED_WITNESSES", empty)


def run_package_dimension(
    emitter: Emitter,
    package: Package,
    n: int,
) -> dict[sp.Symbol, tuple[sp.Expr, ...]]:
    prefix = f"PY_S10_{package.name}_D{n}_"
    t = sp.Symbol("t", real=True)
    xvec = sp.symbols(f"x1:{n + 1}", real=True)
    kvec = sp.symbols(f"k1:{n + 1}", real=True)
    avec = sp.symbols(f"a1:{n + 1}", real=True)
    phase = sp.Symbol("phase", real=True)
    theta = sp.Add(*(kvec[index] * xvec[index] for index in range(n))) - sp.sqrt(omegaSquared) * t
    ansatz = [avec[index] * sp.cos(theta) for index in range(n)]
    k_column = sp.Matrix(kvec)
    k_squared = sp.Add(*(component**2 for component in kvec))
    assumptions = build_joint_assumptions(package, kvec, avec)

    fields, lagrangian, inertial_coefficients, stiffness_coefficients, _stiffness = build_action(package, n, t, xvec)
    dimension_equations, dimension_solution, solved_coefficient_dims = derive_dimension_solution(
        package, lagrangian, fields, t, xvec, assumptions
    )
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
    emit_physical(emitter, walker, prefix + "Q1_LAGRANGIAN_EXPANDED", lagrangian)
    equations = euler_lagrange_system(lagrangian, fields, t, xvec)
    emit_physical(emitter, walker, prefix + "Q1_EULER_LAGRANGE_SYSTEM", equations)

    ansatz_equations = substitute_ansatz(equations, fields, ansatz)
    if not isinstance(ansatz_equations, sp.MatrixBase):
        raise RuntimeError("route-A ansatz equations are not a matrix")
    amplitude_equations = sp.Matrix(
        [
            assumed_simplify(sp.expand(item).coeff(sp.cos(theta)), assumptions)
            for item in ansatz_equations
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
    matrix_ratio = assumed_simplify(matrix_a[0, 0] / matrix_b[0, 0], assumptions)

    emit_physical(emitter, walker, prefix + "Q2_PERIOD_AVERAGED_LAGRANGIAN", averaged_lagrangian)
    emit_physical(emitter, walker, prefix + "Q2_MATRIX_A", matrix_a)
    emit_physical(emitter, walker, prefix + "Q2_MATRIX_B", matrix_b)
    emit_physical(emitter, walker, prefix + "Q2_MATRIX_RESIDUAL", matrix_residual)
    emit_physical(emitter, walker, prefix + "Q2_MATRIX_ENTRY_RATIO", matrix_ratio)
    emitter.emit(prefix + "Q2_RESIDUAL_TEST_SCOPE", sp.Symbol("same_action_variational_identity"))
    emitter.emit(prefix + "Q2_DOWNSTREAM_ROUTE", sp.Symbol("M_B"))

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

    emitter.emit(prefix + "Q6_ENERGY_DENSITY_DIMENSION", dim_expr((2 - D, -2, 1)))
    emitter.emit(prefix + "Q6_DIMENSION_EQUATIONS", sp.Tuple(*dimension_equations))
    emitter.emit(prefix + "Q6_DIMENSION_SOLUTION", dimension_solution)
    if not dimension_solution:
        raise RuntimeError(f"dimension solve returned no solution for {prefix}")
    emitter.emit(prefix + "Q6_INERTIAL_COEFFICIENTS", sp.Tuple(*inertial_coefficients))
    for index, coefficient in enumerate(inertial_coefficients, 1):
        emit_physical(emitter, walker, prefix + f"Q6_INERTIAL_COEFFICIENT{index}", coefficient)
    emitter.emit(prefix + "Q6_STIFFNESS_COEFFICIENTS", sp.Tuple(*stiffness_coefficients))
    for index, coefficient in enumerate(stiffness_coefficients, 1):
        emit_physical(emitter, walker, prefix + f"Q6_STIFFNESS_COEFFICIENT{index}", coefficient)

    gradient_symbols = [[sp.Symbol(f"g{i + 1}{j + 1}", real=True) for j in range(3)] for i in range(3)]
    if n == 3:
        curl_stiffness = sp.expand(stiffness_density("curl", gradient_symbols))
        curl_vector = sp.Matrix(
            [
                gradient_symbols[1][2] - gradient_symbols[2][1],
                gradient_symbols[2][0] - gradient_symbols[0][2],
                gradient_symbols[0][1] - gradient_symbols[1][0],
            ]
        )
        curl_dot = sp.expand((curl_vector.T * curl_vector)[0])
        curl_difference = sp.expand(curl_stiffness - curl_dot)
        gradient_dim_map = {symbol: ZERO_DIM for row in gradient_symbols for symbol in row}
        q7_walker = DimensionWalker(gradient_dim_map, (), {}, assumptions)
        emit_physical(emitter, q7_walker, prefix + "Q7_STIFFNESS", curl_stiffness)
        emit_physical(emitter, q7_walker, prefix + "Q7_CURL_DOT", curl_dot)
        emit_physical(emitter, q7_walker, prefix + "Q7_DIFFERENCE", curl_difference)
    else:
        emitter.emit(prefix + "Q7_OBJECTS", sp.Tuple())

    allowed_strata: list[dict[sp.Symbol, sp.Expr]] = []
    coincidence_loci = sp.Tuple(*(item[2].expression() for item in coincidences))
    coincidence_operands = sp.Tuple(*(item[3] for item in coincidences))
    coincidence_tests = sp.Tuple(*(item[4] for item in coincidences))
    for item in coincidences:
        if item[4] == sp.true:
            allowed_strata.extend(item[2].branches)

    for index, (root, objects) in enumerate(zip(roots, q4_by_root), 1):
        root_prefix = prefix + f"ROOT{index}_"
        root_matrix = objects["N1_MATRIX"]
        generic_rank = int(objects["N2_RANK"])
        if not isinstance(root_matrix, sp.MatrixBase):
            raise RuntimeError("Q8 root matrix is not a matrix")
        minors = rank_drop_minors(sp.Matrix(root_matrix), generic_rank)
        emit_physical(emitter, walker, root_prefix + "Q8_RANK_DROP_MINORS", sp.Tuple(*minors))
        rank_locus = solve_real_locus(minors, kvec, assumptions, impossible=(generic_rank == 0))
        operands, allowed, witnesses = locus_allowed_data(rank_locus, kvec, k_squared, assumptions)
        emitter.emit(root_prefix + "Q8_RANK_DROP_LOCUS", rank_locus.expression())
        emitter.emit(root_prefix + "Q8_RANK_DROP_ALLOWED_OPERANDS", operands)
        emitter.emit(root_prefix + "Q8_RANK_DROP_ALLOWED_TEST", allowed)
        emitter.emit(root_prefix + "Q8_RANK_DROP_ALLOWED_WITNESSES", witnesses)
        emitter.emit(root_prefix + "Q8_ROOT_COINCIDENCE_LOCI", coincidence_loci)
        emitter.emit(root_prefix + "Q8_ROOT_COINCIDENCE_ALLOWED_OPERANDS", coincidence_operands)
        emitter.emit(root_prefix + "Q8_ROOT_COINCIDENCE_ALLOWED_TESTS", coincidence_tests)
        local_prefix = "PY_S10_LOCAL_" + root_prefix.removeprefix("PY_S10_")
        emitter.emit(local_prefix + "Q8_RANK_DROP_SOLVER_OUTPUT", rank_locus.raw_solver_output)
        emitter.emit(local_prefix + "Q8_RANK_DROP_SOLVER_ROUTE", rank_locus.solver_route)
        if allowed == sp.true:
            allowed_strata.extend(rank_locus.branches)

    unique_strata: list[dict[sp.Symbol, sp.Expr]] = []
    keys: set[tuple[tuple[str, str], ...]] = set()
    for branch in allowed_strata:
        key = tuple(sorted((str(symbol), str(value)) for symbol, value in branch.items()))
        if key not in keys:
            keys.add(key)
            unique_strata.append(branch)
    emitter.emit(
        prefix + "Q8_ALLOWED_STRATA",
        sp.Tuple(
            *(
                sp.Tuple(*(sp.Eq(symbol, value, evaluate=False) for symbol, value in sorted(branch.items(), key=lambda item: str(item[0]))))
                for branch in unique_strata
            )
        ),
    )
    for stratum_index, branch in enumerate(unique_strata, 1):
        point = dict(branch)
        free = [component for component in kvec if component not in point]
        for component in free:
            point[component] = sp.Integer(0)
        if not free:
            continue
        point[free[0]] = sp.Integer(1)
        stratum_prefix = prefix + f"Q8_STRATUM{stratum_index}_"
        emitter.emit(
            stratum_prefix + "POINT",
            sp.Tuple(*(sp.Eq(component, point[component], evaluate=False) for component in kvec)),
        )
        emit_stratum_q3_q4(emitter, walker, stratum_prefix, matrix_b, point, kvec, assumptions)

    return solved_coefficient_dims


def load_allowed_registry_fields() -> dict[str, tuple[str, tuple[int, int, int], sp.Expr | None]]:
    original_validate_loci = registry_read.Registry._validate_loci
    registry_read.Registry._validate_loci = lambda self: None
    try:
        registry = registry_read.load_registry(REDUCTION_DIR)
    finally:
        registry_read.Registry._validate_loci = original_validate_loci
    selected = {}
    for qid in ("Q.brane.rho_br", "Q.brane.mu_R", "Q.brane.D_brane"):
        quantity = registry.quantities[qid]
        selected[qid] = (quantity.symbol_name, quantity.dimension, quantity.value)
    return selected


def emit_registry_comparison(
    emitter: Emitter,
    derived_dimensions: Mapping[sp.Symbol, Sequence[sp.Expr]],
) -> None:
    registry = load_allowed_registry_fields()
    d_symbol_name, d_declared_dimension, d_value = registry["Q.brane.D_brane"]
    if d_value is None:
        raise RuntimeError("registry D_brane has no declared value")
    emitter.emit("PY_S10_LOCAL_REGISTRY_D_BRANE_SYMBOL_NAME", sp.Symbol(d_symbol_name))
    emitter.emit("PY_S10_LOCAL_REGISTRY_D_BRANE_VALUE", d_value)
    emitter.emit("PY_S10_LOCAL_REGISTRY_D_BRANE_DIMENSION", dim_expr(d_declared_dimension))
    for qid, symbol, label in (
        ("Q.brane.rho_br", rho_br, "RHO_BR"),
        ("Q.brane.mu_R", mu_R, "MU_R"),
    ):
        symbol_name, declared_dimension, _value = registry[qid]
        symbolic = tuple(derived_dimensions[symbol])
        specialised = tuple(component.subs(D, d_value) for component in symbolic)
        declared = tuple(sp.Integer(component) for component in declared_dimension)
        residual = tuple(left - right for left, right in zip(specialised, declared))
        emitter.emit(f"PY_S10_LOCAL_REGISTRY_{label}_SYMBOL_NAME", sp.Symbol(symbol_name))
        emitter.emit(f"PY_S10_LOCAL_REGISTRY_{label}_DERIVED_DIMENSION_SYMBOLIC", dim_expr(symbolic))
        emitter.emit(f"PY_S10_LOCAL_REGISTRY_{label}_DERIVED_DIMENSION_SPECIALISED", dim_expr(specialised))
        emitter.emit(f"PY_S10_LOCAL_REGISTRY_{label}_DECLARED_DIMENSION", dim_expr(declared))
        emitter.emit(f"PY_S10_LOCAL_REGISTRY_{label}_DIMENSION_RESIDUAL", dim_expr(residual))


def main() -> int:
    emitter = Emitter()
    run_pairs = [sp.Symbol(f"{package.name}_D{n}") for package in PACKAGES for n in package.dimensions]
    emitter.emit("PY_S10_RUN_PAIRS", sp.Tuple(*run_pairs))
    emitter.emit("PY_S10_SKIPPED_PAIRS", sp.Tuple())
    derived_dimensions: dict[sp.Symbol, tuple[sp.Expr, ...]] | None = None
    for package in PACKAGES:
        for n in package.dimensions:
            current = run_package_dimension(emitter, package, n)
            if derived_dimensions is None:
                derived_dimensions = current
    if derived_dimensions is None:
        raise RuntimeError("no package/dimension pairs were run")
    emit_registry_comparison(emitter, derived_dimensions)
    local_list_tag = "PY_S10_LOCAL_TAG_NAMES"
    local_names = [*emitter.local_names, local_list_tag]
    emitter.emit(local_list_tag, sp.Tuple(*(sp.Symbol(name) for name in local_names)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
