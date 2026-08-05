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
OPERATION_TIMEOUT_SECONDS = 60


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
        failed_head = sp.Symbol(getattr(expr.func, "__name__", type(expr).__name__))
        marker = sp.Function("indeterminate_dimension")(failed_head)
        return (marker, marker, marker)

    def dimension(self, expr: sp.Expr) -> tuple[sp.Expr, ...]:
        if expr is sp.zoo or expr is sp.nan or expr in (sp.oo, -sp.oo):
            return ZERO_DIM
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
            return self.indeterminate_dimension(expr)
        return self.indeterminate_dimension(expr)

    def component_term_dimensions(self, expr: sp.Expr) -> tuple[tuple[sp.Expr, ...], ...]:
        expanded = sp.expand(expr)
        if expanded == 0:
            return ()
        return tuple(self.dimension(term) for term in sp.Add.make_args(expanded))

    def nested_add_dimensions(self, expr: sp.Expr) -> tuple[sp.Tuple, ...]:
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
        elif isinstance(obj, sp.Expr):
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
    emitter.emit(name, payload)
    dimension_source = payload if dimension_payload is None else dimension_payload
    term_dimensions, dimensions, nested_add_dimensions, homogeneous = walker.report(
        dimension_source
    )
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
) -> tuple[list[sp.Expr], sp.Expr, list[sp.Expr], list[sp.Expr]]:
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
    return fields, sp.expand(lagrangian), inertial_coefficients, [stiffness_coefficient]


def action_stiffness_density(
    lagrangian: sp.Expr,
    fields: Sequence[sp.Expr],
    t: sp.Symbol,
    xvec: Sequence[sp.Symbol],
    gradient_symbols: Sequence[Sequence[sp.Symbol]],
    stiffness_coefficient: sp.Expr,
    stiffness_sign: int,
) -> sp.Expr:
    substitutions = {
        **{sp.diff(field, t): sp.Integer(0) for field in fields},
        **{
            sp.diff(fields[j], xvec[i]): gradient_symbols[i][j]
            for i in range(len(xvec))
            for j in range(len(fields))
        },
    }
    stiffness_action_term = sp.expand(lagrangian.xreplace(substitutions))
    action_prefactor = sp.Integer(stiffness_sign) * sp.Rational(1, 2) * stiffness_coefficient
    return sp.expand(sp.cancel(stiffness_action_term / action_prefactor))


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
    list[dict[sp.Symbol, sp.Expr]],
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
        coefficient: sp.symbols(
            f"{coefficient}_dim_length {coefficient}_dim_time {coefficient}_dim_mass"
        )
        for coefficient in dimensionful_coefficient_symbols
    }
    unknowns = tuple(
        component
        for coefficient in dimensionful_coefficient_symbols
        for component in coefficient_unknowns[coefficient]
    )
    unknown_map = {
        **coefficient_unknowns,
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
    coefficient_matrix, right_hand_side = sp.linear_eq_to_matrix(equations, unknowns)
    coefficient_rank = int(coefficient_matrix.rank())
    augmented_rank = int(coefficient_matrix.row_join(right_hand_side).rank())
    independent_equation_count = augmented_rank
    difference = independent_equation_count - len(unknowns)
    determination = sp.Symbol(
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
) -> tuple[sp.Expr, sp.Symbol]:
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
    return sp.Symbol("undecided_under_joint_assumptions")


def undecided_sign_subexpression(
    expr: sp.Expr,
    assumptions: sp.logic.boolalg.Boolean,
    sign_result: sp.Expr,
) -> sp.Expr:
    if sign_result != sp.Symbol("undecided_under_joint_assumptions"):
        return sp.Symbol("none_sign_decided")
    refined = sp.refine(sp.factor(expr), assumptions)
    unsettled = []
    for factor in sp.Mul.make_args(refined):
        answers = (
            sp.ask(sp.Q.positive(factor), assumptions),
            sp.ask(sp.Q.zero(factor), assumptions),
            sp.ask(sp.Q.negative(factor), assumptions),
        )
        if all(answer is None for answer in answers):
            unsettled.append(factor)
    if len(unsettled) == 1:
        return unsettled[0]
    if unsettled:
        return sp.Tuple(*unsettled)
    return refined


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
) -> tuple[sp.Expr, sp.Expr, sp.Symbol, list[dict[sp.Symbol, sp.Expr]]]:
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
    solver_route: sp.Symbol
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
            sp.Symbol("rank_floor_empty_locus"),
        )

    if not normalized:
        real_domain = sp.ProductSet(*(sp.S.Reals for _ in variables))
        return RealLocus(
            (),
            real_domain,
            ({},),
            ((),),
            sp.Tuple(sp.Tuple(sp.Tuple(), sp.Tuple(), sp.Symbol("retained_real_branch"))),
            sp.Symbol("universal_real_locus"),
            real_domain,
        )

    previous_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(OPERATION_TIMEOUT_SECONDS)
    try:
        raw = sp.solve(normalized, list(variables), dict=True)
        route = sp.Symbol("solve_then_explicit_real_filter")
        fallback_expression: object = sp.EmptySet
    except OperationTimeout:
        raw = sp.ConditionSet(
            sp.Tuple(*variables),
            sp.And(*(sp.Eq(eq, 0) for eq in normalized)),
            sp.ProductSet(*(sp.S.Reals for _ in variables)),
        )
        route = sp.Symbol("solve_timeout_explicit_real_conditionset")
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
                    sp.Symbol("discarded_nonreal_branch"),
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
                sp.Symbol("retained_real_branch"),
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
    allowed: sp.Expr
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
    branch_reasons: list[sp.Symbol] = []
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
            branch_reasons.append(sp.Symbol("locus_conflicts_with_positive_wavevector_norm"))
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
            branch_reasons.append(sp.Symbol("allowed_witness_found"))
        elif not free:
            branch_tests.append(sp.Symbol("undecided_no_explicit_real_witness"))
            branch_witnesses.append(sp.Tuple())
            branch_points.append(None)
            branch_reasons.append(sp.Symbol("fully_fixed_branch_not_materialized"))
        else:
            branch_tests.append(sp.Symbol("undecided_no_explicit_real_witness"))
            branch_witnesses.append(sp.Tuple())
            branch_points.append(None)
            branch_reasons.append(sp.Symbol("no_explicit_real_witness"))

    if any(test is sp.true for test in branch_tests):
        allowed: sp.Expr = sp.true
    elif branch_tests and all(test is sp.false for test in branch_tests):
        allowed = sp.false
    elif locus.fallback_expression is not sp.EmptySet and not locus.branches:
        allowed = sp.Symbol("undecided_solver_conditionset")
    else:
        allowed = sp.Symbol("undecided_no_explicit_real_witness") if branch_tests else sp.false
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
        ratio = sp.Symbol("undefined_zero_root_ratio")
        ratio_defined = sp.false
        exponent = sp.Symbol("undefined_zero_root_ratio")
    else:
        ratio_object = sp.Mul(scaled, sp.Pow(original, -1, evaluate=False), evaluate=False)
        ratio = assumed_simplify(ratio_object, assumptions)
        ratio_defined = sp.true
        try:
            polynomial = sp.Poly(ratio, lambdaScale)
            exponent = (
                sp.Integer(polynomial.degree())
                if len(polynomial.terms()) == 1
                else sp.Symbol("not_a_pure_lambdaScale_power")
            )
        except (sp.PolynomialError, TypeError, ValueError):
            exponent = sp.Symbol("not_a_pure_lambdaScale_power")
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
    emitter.emit(prefix + "PREMISE_U_DIMENSION", dim_expr(LENGTH_DIM))
    emitter.emit(prefix + "PREMISE_ANSATZ", sp.Tuple(*ansatz))
    emitter.emit(
        prefix + "PREMISE_PERIOD_AVERAGE",
        sp.Tuple(sp.Symbol("phase"), sp.Integer(0), 2 * sp.pi, sp.Rational(1, 2) / sp.pi),
    )
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
    roots_data: tuple[sp.Expr, sp.Expr, sp.Symbol, list[dict[sp.Symbol, sp.Expr]]],
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
        emitter.emit(prefix + "Q3_SPECTRUM_SOLVE_CONDITION", sp.Symbol("no_roots_returned"))
        emitter.emit(
            prefix + "Q3_SPECTRUM_SOLVE_CONDITION_OPERANDS",
            sp.Tuple(factored, sp.Symbol("omegaSquared")),
        )
    else:
        emitter.emit(prefix + "Q3_SPECTRUM_SOLVE_CONDITION", sp.Symbol("roots_returned"))
        emitter.emit(
            prefix + "Q3_SPECTRUM_SOLVE_CONDITION_OPERANDS",
            sp.Tuple(factored, sp.Symbol("omegaSquared")),
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
    wave_number_scale = sp.Symbol(
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

    fields, lagrangian, inertial_coefficients, stiffness_coefficients = build_action(package, n, t, xvec)
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
        matrix_ratio = sp.Symbol("undefined_zero_denominator")
    else:
        computed_ratio = assumed_simplify(
            matrix_ratio_numerator / matrix_ratio_denominator, assumptions
        )
        if isinstance(computed_ratio, sp.Expr) and computed_ratio.has(sp.zoo, sp.nan):
            matrix_ratio = sp.Symbol("undefined_nonfinite_ratio")
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
        sign_result = derived_sign(root, assumptions)
        emitter.emit(root_prefix + "Q3_SIGN", sign_result)
        local_root_prefix = root_prefix.replace("PY_S10_", "PY_S10_LOCAL_", 1)
        emitter.emit(local_root_prefix + "Q3_SIGN_JOINT_ASSUMPTIONS", assumptions)
        emitter.emit(
            local_root_prefix + "Q3_SIGN_UNDECIDED_SUBEXPRESSION",
            undecided_sign_subexpression(root, assumptions, sign_result),
        )
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
        sp.Symbol("solution_returned") if dimension_solution else sp.Symbol("no_solution_returned"),
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

    if n == 3:
        gradient_symbols = [
            [sp.Symbol(f"g{i + 1}{j + 1}", real=True) for j in range(n)]
            for i in range(n)
        ]
        package_stiffness = action_stiffness_density(
            lagrangian,
            fields,
            t,
            xvec,
            gradient_symbols,
            stiffness_coefficients[0],
            package.stiffness_sign,
        )
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

    stratum_candidates: dict[str, tuple[sp.Tuple, sp.Expr, dict[sp.Symbol, sp.Expr] | None, sp.Symbol]] = {}

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
        emitter.emit(stratum_prefix + "SKIP_STATUS", sp.Symbol("not_skipped_allowed_branch"))
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
        sp.Tuple(*(sp.Symbol(f"{name}_D{n}") for name, n in completed_pairs)),
    )
    emitter.emit(
        "PY_S10_SKIPPED_PAIRS",
        sp.Tuple(*(sp.Symbol(f"{name}_D{n}") for name, n in skipped_pairs)),
    )
    emit_registry_comparison(emitter, derived_dimensions)
    local_list_tag = "PY_S10_LOCAL_TAG_NAMES"
    local_names = [*emitter.local_names, local_list_tag]
    emitter.emit(local_list_tag, sp.Tuple(*(sp.Symbol(name) for name in local_names)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
