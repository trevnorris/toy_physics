#!/usr/bin/env python3
"""Independent SymPy audit for ledger step S9."""

import re
from pathlib import Path

import sympy as sp

from extract_knob_inventory import CLASS_TAGS, declaration_classes


# Field, coordinate, amplitude, and material symbols.
t, x, y, z = sp.symbols("t x y z", real=True)  # COORDINATE · spacetime coordinates
coordinates = (t, x, y, z)  # COORDINATE · spacetime coordinate collection
spatial_coordinates = (x, y, z)  # COORDINATE · spatial coordinate collection
u_functions = tuple(sp.Function(name)(t, x, y, z) for name in ("u1", "u2", "u3"))  # PREMISE · brane displacement field heads
a_symbols = sp.symbols("a1 a2 a3", real=True)  # PREMISE · plane-wave amplitudes
b_symbols = sp.symbols("b1 b2 b3", real=True)  # PREMISE · paired plane-wave amplitudes
omega = sp.symbols("omega", real=True)  # COORDINATE · spectral-frequency variable
omega2 = sp.symbols("omega2", real=True)  # COORDINATE · squared-frequency polynomial variable
k_input = sp.symbols("kx ky kz", real=True)  # COORDINATE · wavevector components
rho_br = sp.Symbol("rho_br", positive=True)  # KNOB · brane inertia density not derived here
mu_R = sp.Symbol("mu_R", positive=True)  # KNOB · brane shear modulus not derived here
rho_z = sp.Symbol("rho_z", positive=True)  # CONTROL · anisotropic inertia coefficient
mu_F = sp.Symbol("mu_F", positive=True)  # CONTROL · flexural ablation coefficient
mu_G = sp.Symbol("mu_G", positive=True)  # CONTROL · bare-field ablation coefficient
lambda_rho, lambda_mu = sp.symbols("lambda_rho lambda_mu", positive=True)  # CONTROL · coefficient ablation scales
D = sp.symbols("D", integer=True, positive=True)  # STRUCTURAL · brane spatial dimension
q = sp.symbols("q", positive=True)  # COORDINATE · wave-norm placeholder
lambda_scale = sp.symbols("lambda_scale", positive=True)  # CONTROL · wavevector homogeneity scale


def curl_of(field):
    """Construct a three-dimensional curl from a field triple."""
    return sp.Matrix(
        [
            sp.diff(field[2], y) - sp.diff(field[1], z),
            sp.diff(field[0], z) - sp.diff(field[2], x),
            sp.diff(field[1], x) - sp.diff(field[0], y),
        ]
    )


def construct_curl_action(inertia, stiffness, stiffness_sign=-1):
    velocity = sp.Matrix([sp.diff(component, t) for component in u_functions])
    curl = curl_of(u_functions)
    return sp.expand(
        sp.Rational(1, 2) * (velocity.T * inertia * velocity)[0]
        + stiffness_sign * sp.Rational(1, 2) * stiffness * (curl.T * curl)[0]
    )


def construct_divergence_action():
    velocity = sp.Matrix([sp.diff(component, t) for component in u_functions])
    divergence = sum(sp.diff(u_functions[index], spatial_coordinates[index]) for index in range(3))
    return sp.expand(
        sp.Rational(1, 2) * rho_br * (velocity.T * velocity)[0]
        - sp.Rational(1, 2) * mu_R * divergence**2
    )


def construct_gradient_action():
    velocity = sp.Matrix([sp.diff(component, t) for component in u_functions])
    gradient_square = sum(
        sp.diff(component, coordinate) ** 2
        for component in u_functions
        for coordinate in spatial_coordinates
    )
    return sp.expand(
        sp.Rational(1, 2) * rho_br * (velocity.T * velocity)[0]
        - sp.Rational(1, 2) * mu_R * gradient_square
    )


def construct_flexural_action():
    velocity = sp.Matrix([sp.diff(component, t) for component in u_functions])
    laplacian = sp.Matrix(
        [sum(sp.diff(component, coordinate, 2) for coordinate in spatial_coordinates) for component in u_functions]
    )
    return sp.expand(
        sp.Rational(1, 2) * rho_br * (velocity.T * velocity)[0]
        - sp.Rational(1, 2) * mu_F * (laplacian.T * laplacian)[0]
    )


def construct_bare_field_action():
    field = sp.Matrix(u_functions)
    return sp.expand(
        construct_curl_action(rho_br * sp.eye(3), mu_R)
        - sp.Rational(1, 2) * mu_G * (field.T * field)[0]
    )


# Every physical modification below enters through an action constructor.
identity3 = sp.eye(len(spatial_coordinates))  # DERIVED · spatial identity tensor
main_action = construct_curl_action(rho_br * identity3, mu_R)  # PREMISE · curl-only MacCullagh action
actions = {  # CONTROL · action package control specifications
    "MAIN": (main_action, (rho_br, mu_R), rho_br, mu_R, ()),
    "X1": (
        construct_curl_action(lambda_rho * rho_br * identity3, mu_R),
        (rho_br, mu_R),
        rho_br,
        mu_R,
        (lambda_rho,),
    ),
    "X2": (
        construct_curl_action(rho_br * identity3, lambda_mu * mu_R),
        (rho_br, mu_R),
        rho_br,
        mu_R,
        (lambda_mu,),
    ),
    "X3": (construct_divergence_action(), (rho_br, mu_R), rho_br, mu_R, ()),
    "X4": (construct_gradient_action(), (rho_br, mu_R), rho_br, mu_R, ()),
    "X5": (construct_curl_action(rho_br * identity3, mu_R, stiffness_sign=1), (rho_br, mu_R), rho_br, mu_R, ()),
    "X6": (
        construct_curl_action(sp.diag(rho_br, rho_br, rho_z), mu_R),
        (rho_br, rho_z, mu_R),
        rho_br,
        mu_R,
        (),
    ),
    "X7": (construct_flexural_action(), (rho_br, mu_F), rho_br, mu_F, ()),
    "X8": (construct_bare_field_action(), (rho_br, mu_R, mu_G), rho_br, mu_R, ()),
}


# The prescribed plane wave is the second and final hand-constructed physical expression.
position = sp.Matrix(spatial_coordinates)  # COORDINATE · spatial position vector
input_wavevector = sp.Matrix(k_input)  # COORDINATE · wavevector component collection
phase_exponent = sp.I * ((input_wavevector.T * position)[0] - omega * t)  # PREMISE · plane-wave phase prescription
phase = sp.exp(phase_exponent)  # PREMISE · plane-wave phase factor
plane_wave = sp.Matrix([amplitude * phase for amplitude in a_symbols])  # PREMISE · prescribed plane-wave ansatz

# Recover wave data by differentiating the constructed phase rather than retyping its combinations.
derived_wavevector = sp.Matrix([-sp.I * sp.diff(phase, coordinate) / phase for coordinate in spatial_coordinates])  # DERIVED · phase-gradient wavevector
wave_norm = sp.factor((derived_wavevector.T * derived_wavevector)[0])  # DERIVED · wavevector norm form
transverse_generator = sp.simplify(wave_norm * identity3 - derived_wavevector * derived_wavevector.T)  # DERIVED · transverse generator
longitudinal_generator = sp.simplify(derived_wavevector * derived_wavevector.T)  # DERIVED · longitudinal generator


def euler_lagrange_residual(action):
    entries = []
    for field in u_functions:
        entry = sp.diff(action, field)
        field_derivatives = sorted(
            (derivative for derivative in action.atoms(sp.Derivative) if derivative.expr == field),
            key=sp.default_sort_key,
        )
        for derivative in field_derivatives:
            varied_term = sp.diff(action, derivative)
            for variable, count in derivative.variable_count:
                varied_term = sp.diff(varied_term, variable, count)
            entry += (-1) ** derivative.derivative_count * varied_term
        entries.append(sp.simplify(entry))
    return sp.Matrix(entries)


def route_a(action):
    residual = euler_lagrange_residual(action)
    substitution = dict(zip(u_functions, plane_wave))
    substituted = residual.subs(substitution).doit()
    reduced = sp.Matrix([sp.simplify(entry / phase).subs(omega**2, omega2) for entry in substituted])
    matrix = sp.simplify(reduced.jacobian(sp.Matrix(a_symbols)))
    return residual, reduced, matrix


def route_b(action):
    opposite_phase = sp.exp(-phase_exponent)
    paired_wave = tuple(
        a_symbols[index] * phase + b_symbols[index] * opposite_phase for index in range(3)
    )
    substituted_action = action.subs(dict(zip(u_functions, paired_wave))).doit()
    mixed_hessian = sp.Matrix(
        3,
        3,
        lambda row, column: sp.simplify(
            sp.diff(substituted_action, b_symbols[row], a_symbols[column]).subs(omega**2, omega2)
        ),
    )
    return mixed_hessian


def roots_from_factorization(factored_determinant):
    factor_coefficient, factors = sp.factor_list(factored_determinant, omega2)
    del factor_coefficient
    multiplicities = {}
    for factor, factor_multiplicity in factors:
        factor_roots = sp.roots(factor, omega2)
        for root, root_multiplicity in factor_roots.items():
            multiplicities[sp.factor(root)] = factor_multiplicity * root_multiplicity
    ordered = sorted(multiplicities, key=sp.default_sort_key)
    multiset = [root for root in ordered for _ in range(multiplicities[root])]
    return ordered, multiplicities, multiset


def symbolic_zero(value):
    if isinstance(value, sp.MatrixBase):
        return sp.S.true if all(sp.simplify(entry) == 0 for entry in value) else sp.S.false
    return sp.S.true if sp.simplify(value) == 0 else sp.S.false


def root_tests(matrix, roots, specialised_wavevector):
    row = specialised_wavevector.T
    specialised_norm = sp.factor((row * specialised_wavevector)[0])
    specialised_transverse = sp.simplify(specialised_norm * identity3 - specialised_wavevector * row)
    records = []
    for root in roots:
        root_matrix = matrix.subs(omega2, root).applyfunc(sp.factor)
        stacked = root_matrix.col_join(row)
        stacked_rank = stacked.rank()
        e1 = sp.S.true if stacked_rank < 3 else sp.S.false
        e2 = sp.Integer(3 - stacked_rank)
        mk = (root_matrix * specialised_wavevector).applyfunc(sp.factor)
        mt = (root_matrix * specialised_transverse).applyfunc(sp.factor)
        e3 = symbolic_zero(mk)
        e4 = symbolic_zero(mt)
        nullity = sp.Integer(3 - root_matrix.rank())
        records.append((root, stacked, mk, mt, e1, e2, e3, e4, nullity))
    return records


def evaluated_root_test_block(matrix, specialised_wavevector):
    determinant = sp.factor(matrix.det())
    roots, _, _ = roots_from_factorization(determinant)
    records = root_tests(matrix, roots, specialised_wavevector)
    row = specialised_wavevector.T
    specialised_norm = sp.factor((row * specialised_wavevector)[0])
    specialised_transverse = sp.simplify(specialised_norm * identity3 - specialised_wavevector * row)
    return {
        "E1": [sp.Tuple(record[0], record[1], record[4]) for record in records],
        "E2": [sp.Tuple(record[0], record[1].rank(), record[5]) for record in records],
        "E3": [sp.Tuple(record[0], specialised_wavevector, record[2], record[6]) for record in records],
        "E4": [sp.Tuple(record[0], specialised_transverse, record[3], record[7]) for record in records],
    }


def q_form(expression):
    return sp.factor(sp.factor(expression).subs(wave_norm, q))


dkey_computation_dimensions = "_DKEY_COMPUTATION_DIMENSIONS"  # CONTROL · per-output production-dimension metadata key


def produced_outputs(entries, computation_dimension):
    output = dict(entries)
    output[dkey_computation_dimensions] = {
        name: computation_dimension
        for name in entries
        if not name.startswith("_")
    }
    return output


def extend_produced_outputs(target, source):
    target.update({
        name: value
        for name, value in source.items()
        if name != dkey_computation_dimensions
    })
    target[dkey_computation_dimensions].update(source[dkey_computation_dimensions])


def add_produced_outputs(target, entries, computation_dimension):
    extend_produced_outputs(target, produced_outputs(entries, computation_dimension))


length_dimension = sp.Matrix([1, 0, 0])  # STRUCTURAL · length unit-basis marker
time_dimension = sp.Matrix([0, 1, 0])  # STRUCTURAL · time unit-basis marker
mass_dimension = sp.Matrix([0, 0, 1])  # STRUCTURAL · mass unit-basis marker
acceleration_dimension = length_dimension - 2 * time_dimension  # DERIVED · computed on the way to a supplied reference
force_dimension = mass_dimension + acceleration_dimension  # DERIVED · computed on the way to a supplied reference
energy_dimension = force_dimension + length_dimension  # DERIVED · computed on the way to a supplied reference
energy_density_dimension = energy_dimension - D * length_dimension  # PREMISE · supplied reference used as the dimension-solve target
field_dimension = length_dimension  # PREMISE · supplied reference used by the dimension walk
q_dimension = -2 * length_dimension  # PREMISE · supplied reference used by the dimension walk
squared_velocity_dimension = 2 * (length_dimension - time_dimension)  # PREMISE · supplied reference used for the speed-dimension residual


class DimensionWalkError(Exception):
    def __init__(self, kind, payload):
        super().__init__(kind)
        self.kind = kind
        self.payload = payload


def derivative_order(derivative):
    counts = {coordinate: 0 for coordinate in coordinates}
    for variable, count in derivative.variable_count:
        counts[variable] += count
    return sp.Tuple(*(sp.Integer(counts[coordinate]) for coordinate in coordinates))


def field_factor_multiorders(expression):
    if expression in u_functions:
        return [sp.Tuple(0, 0, 0, 0)]
    if isinstance(expression, sp.Derivative) and expression.expr in u_functions:
        return [derivative_order(expression)]
    if expression.is_Pow:
        base, exponent = expression.args
        if exponent.is_integer and exponent.is_nonnegative:
            return field_factor_multiorders(base) * int(exponent)
        return []
    orders = []
    for argument in expression.args:
        orders.extend(field_factor_multiorders(argument))
    return orders


def expression_dimension(expression, dimension_lookup, dimensionless_symbols):
    zero_dimension = sp.zeros(3, 1)
    symbol_dimensions = {
        t: time_dimension,
        x: length_dimension,
        y: length_dimension,
        z: length_dimension,
        omega: -time_dimension,
        omega2: -2 * time_dimension,
        q: q_dimension,
        D: zero_dimension,
    }
    symbol_dimensions.update({component: -length_dimension for component in k_input})
    symbol_dimensions.update({component: field_dimension for component in (*a_symbols, *b_symbols)})
    symbol_dimensions.update({field: field_dimension for field in u_functions})
    symbol_dimensions.update({symbol: zero_dimension for symbol in dimensionless_symbols})
    symbol_dimensions.update(dimension_lookup)

    def walk(node):
        if node in symbol_dimensions:
            return symbol_dimensions[node]
        if node.is_number:
            return zero_dimension
        if isinstance(node, sp.Symbol):
            raise DimensionWalkError("unknown_symbol", node)
        if isinstance(node, sp.Derivative):
            dimension = walk(node.expr)
            for variable, count in node.variable_count:
                dimension -= count * walk(variable)
            return sp.simplify(dimension)
        if node.is_Add:
            summand_records = [(summand, walk(summand)) for summand in node.args]
            reference = summand_records[0][1]
            if any(symbolic_zero(dimension - reference) != sp.S.true for _, dimension in summand_records[1:]):
                raise DimensionWalkError("sum", summand_records)
            return sp.simplify(reference)
        if node.is_Mul:
            return sp.simplify(sum((walk(factor) for factor in node.args), zero_dimension))
        if node.is_Pow:
            base, exponent = node.args
            if not exponent.is_Rational:
                raise DimensionWalkError("power", sp.Tuple(base, exponent))
            return sp.simplify(exponent * walk(base))
        if node.func == sp.exp:
            argument_dimension = walk(node.args[0])
            if symbolic_zero(argument_dimension) != sp.S.true:
                raise DimensionWalkError("function_argument", sp.Tuple(node.args[0], argument_dimension))
            return zero_dimension
        raise DimensionWalkError("unknown_node", node)

    return walk(sp.sympify(expression))


def dimension_block(action, dimensional_coefficients, primary_inertia, primary_stiffness, dimensionless_symbols):
    dimension_variables = {
        coefficient: sp.Matrix(
            [sp.Symbol(f"dim_{coefficient}_{axis}") for axis in ("L", "T", "M")]  # DERIVED · coefficient-axis solver unknowns
        )
        for coefficient in dimensional_coefficients
    }
    terms = sorted(sp.Add.make_args(sp.expand(action)), key=sp.default_sort_key)
    orders_by_term = []
    term_dimensions = []
    unsolved_dimension_lookup = dict(dimension_variables)
    for term in terms:
        orders_by_term.append(sp.Tuple(*field_factor_multiorders(term)))
        term_dimensions.append(expression_dimension(term, unsolved_dimension_lookup, dimensionless_symbols))
    equations = [
        sp.Eq(term_dimension[index], energy_density_dimension[index], evaluate=False)
        for term_dimension in term_dimensions
        for index in range(3)
    ]
    unknowns = [entry for coefficient in dimensional_coefficients for entry in dimension_variables[coefficient]]
    solutions = sp.solve(equations, unknowns, dict=True)
    solve_inputs = produced_outputs({
        "DIM_TERMS": terms,
        "DIM_FIELD_MULTIORDERS": orders_by_term,
        "DIM_TERM_EXPRESSIONS": term_dimensions,
        "DIM_LINEAR_SYSTEM": equations,
    }, sp.Integer(len(spatial_coordinates)))
    solve_outputs = produced_outputs({
        "DIM_ENERGY_DENSITY": energy_density_dimension,
        "FIELD_DIMENSION": field_dimension,
        "WAVEVECTOR_NORM_DIMENSION": q_dimension,
        "DIM_SQUARED_VELOCITY": squared_velocity_dimension,
    }, D)
    add_produced_outputs(solve_outputs, {
        "DIM_SOLUTION": solutions,
    }, D)
    extend_produced_outputs(solve_inputs, solve_outputs)
    if not solutions:
        solve_inputs["_DIMENSION_SOLVE_FAILED"] = sp.S.true
        return solve_inputs
    solution = solutions[0]
    coefficient_dimensions = [
        sp.Tuple(coefficient, sp.Matrix([solution[entry] for entry in dimension_variables[coefficient]]))
        for coefficient in dimensional_coefficients
    ]
    inertia_dimension = sp.Matrix([solution[entry] for entry in dimension_variables[primary_inertia]])
    stiffness_dimension = sp.Matrix([solution[entry] for entry in dimension_variables[primary_stiffness]])
    add_produced_outputs(solve_inputs, {
        "DIM_COEFFICIENTS": coefficient_dimensions,
        "DIM_PRIMARY_INERTIA": inertia_dimension,
        "DIM_PRIMARY_STIFFNESS": stiffness_dimension,
        "DIM_STIFFNESS_MINUS_INERTIA": sp.simplify(stiffness_dimension - inertia_dimension),
    }, D)
    fixed_component_count = sp.Integer(len(spatial_coordinates))
    add_produced_outputs(solve_inputs, {
        "DIM_SOLUTION_D3": [
            {key: value.subs(D, fixed_component_count) for key, value in item.items()}
            for item in solutions
        ],
        "DIM_PRIMARY_INERTIA_D3": inertia_dimension.subs(D, fixed_component_count),
        "DIM_PRIMARY_STIFFNESS_D3": stiffness_dimension.subs(D, fixed_component_count),
        "DIM_STIFFNESS_MINUS_INERTIA_D3": (
            stiffness_dimension - inertia_dimension
        ).subs(D, fixed_component_count),
    }, fixed_component_count)
    solve_inputs["_DIMENSION_LOOKUP"] = {
        coefficient: sp.Matrix([solution[entry] for entry in dimension_variables[coefficient]])
        for coefficient in dimensional_coefficients
    }
    return solve_inputs


def assumptions_for(dimensional_coefficients, dimensionless_symbols):
    predicates = [
        sp.Q.real(symbol) for symbol in (*a_symbols, *k_input, omega2)
    ]
    predicates.extend(sp.Q.positive(symbol) for symbol in dimensional_coefficients)
    predicates.extend(sp.Q.positive(symbol) for symbol in dimensionless_symbols)
    predicates.extend((sp.Q.integer(D), sp.Q.positive(D), sp.Q.nonzero(wave_norm), sp.Q.positive(q)))
    return sp.And(*predicates)


def direction_block(matrix):
    substitutions = {
        "DIR_GENERIC": {},
        "DIR_XY": {k_input[2]: 0},
        "DIR_Z": {k_input[0]: 0, k_input[1]: 0},
        "DIR_X": {k_input[1]: 0, k_input[2]: 0},
        "DIR_DIAGONAL": {k_input[1]: k_input[0], k_input[2]: k_input[0]},
    }
    output = {}
    for name, substitution in substitutions.items():
        specialised_matrix = matrix.subs(substitution).applyfunc(sp.factor)
        specialised_wavevector = derived_wavevector.subs(substitution)
        determinant = sp.factor(specialised_matrix.det())
        roots, multiplicities, multiset = roots_from_factorization(determinant)
        records = root_tests(specialised_matrix, roots, specialised_wavevector)
        output[f"{name}_ROOTS"] = roots
        output[f"{name}_ROOT_MULTISET"] = multiset
        output[f"{name}_ROOT_MULTIPLICITIES"] = [sp.Tuple(root, multiplicities[root]) for root in roots]
        output[f"{name}_STACKED"] = [sp.Tuple(record[0], record[1]) for record in records]
        output[f"{name}_MK"] = [sp.Tuple(record[0], record[2]) for record in records]
        output[f"{name}_MT"] = [sp.Tuple(record[0], record[3]) for record in records]
        for offset, test_name in enumerate(("E1", "E2", "E3", "E4"), start=4):
            output[f"{name}_{test_name}"] = [sp.Tuple(record[0], record[offset]) for record in records]
        for offset, test_name in ((4, "E1"), (6, "E3"), (7, "E4")):
            output[f"{name}_ROOTS_PASSING_{test_name}"] = [record[0] for record in records if record[offset] == sp.S.true]
        output[f"{name}_ROOTS_PASSING_E2"] = [record[0] for record in records if record[5] > 0]
        output[f"{name}_ROOT_NULLITIES"] = [sp.Tuple(record[0], record[8]) for record in records]
    return produced_outputs(output, sp.Integer(len(spatial_coordinates)))


def dimension_walk_failure_output(error):
    if error.kind == "sum":
        return {"DIM_WALK_SUMMAND_DIMENSIONS": [sp.Tuple(summand, dimension) for summand, dimension in error.payload]}
    if error.kind == "unknown_symbol":
        return {"DIM_WALK_UNKNOWN_SYMBOL": error.payload}
    if error.kind == "power":
        return {"DIM_WALK_POWER": error.payload}
    if error.kind == "function_argument":
        return {"DIM_WALK_FUNCTION_ARGUMENT": error.payload}
    return {"DIM_WALK_NODE": error.payload}


def derive(
    action,
    dimensional_coefficients,
    primary_inertia,
    primary_stiffness,
    dimensionless_symbols,
    assumption_object,
):
    el_residual, plane_wave_residual, matrix_a = route_a(action)
    matrix_b = route_b(action)
    matrix_residual = sp.simplify(matrix_a - matrix_b)
    determinant = sp.factor(matrix_a.det())
    solved_roots = sorted((sp.factor(root) for root in sp.solve(determinant, omega2)), key=sp.default_sort_key)
    roots, multiplicities, root_multiset = roots_from_factorization(determinant)
    records = root_tests(matrix_a, roots, derived_wavevector)
    passing_e1 = [record[0] for record in records if record[4] == sp.S.true]
    passing_e2 = [record[0] for record in records if record[5] > 0]
    passing_e3 = [record[0] for record in records if record[6] == sp.S.true]
    passing_e4 = [record[0] for record in records if record[7] == sp.S.true]
    transverse_q_roots = [q_form(root) for root in passing_e1]
    transverse_q_wavevector_occurrences = [
        sp.Tuple(
            *(sp.Tuple(symbol, sp.S.true if root.has(symbol) else sp.S.false) for symbol in k_input)
        )
        for root in transverse_q_roots
    ]
    scaled_roots = [
        root.subs(
            {component: lambda_scale * component for component in k_input},
            simultaneous=True,
        )
        for root in roots
    ]
    quadratic_scaled_roots = [lambda_scale**2 * root for root in roots]
    scaling_residuals = [
        sp.simplify(scaled - quadratic_scaled)
        for scaled, quadratic_scaled in zip(scaled_roots, quadratic_scaled_roots)
    ]
    root_sign_assumptions = sp.And(assumption_object, sp.Q.positive(wave_norm))
    root_signs = [sp.refine(sp.sign(root), root_sign_assumptions) for root in roots]
    speed_candidates = [sp.factor(root / q) for root in transverse_q_roots]
    homogeneity_defects = [sp.factor(q * sp.diff(root, q) - root) for root in transverse_q_roots]
    speed_derivatives = [sp.factor(sp.diff(speed, q)) for speed in speed_candidates]
    degree_pairs = [
        (
            sp.Poly(sp.together(root).as_numer_denom()[0], q).degree(),
            sp.Poly(sp.together(root).as_numer_denom()[1], q).degree(),
        )
        for root in transverse_q_roots
    ]
    try:
        dimensions = dimension_block(
            action,
            dimensional_coefficients,
            primary_inertia,
            primary_stiffness,
            dimensionless_symbols,
        )
    except DimensionWalkError as error:
        return {}, {}, dimension_walk_failure_output(error)
    if dimensions.get("_DIMENSION_SOLVE_FAILED") == sp.S.true:
        return {}, {}, {"DIM_SOLUTION": dimensions["DIM_SOLUTION"]}
    try:
        speed_dimensions = [
            expression_dimension(speed, dimensions["_DIMENSION_LOOKUP"], dimensionless_symbols)
            for speed in speed_candidates
        ]
    except DimensionWalkError as error:
        return {}, {}, dimension_walk_failure_output(error)
    output = produced_outputs({
        "LAGRANGIAN": action,
        "EL_RESIDUAL": el_residual,
        "EQUATION_OF_MOTION": sp.Tuple(*(sp.Eq(entry, 0) for entry in el_residual)),
        "PLANE_WAVE_ANSATZ": plane_wave,
        "PLANE_WAVE_RESIDUAL": plane_wave_residual,
        "ROUTE_OPERANDS_AND_RESIDUAL": sp.Tuple(matrix_a, matrix_b, matrix_residual),
        "DET_M_FACTORED": determinant,
        "OMEGA2_SOLUTIONS": sp.FiniteSet(*solved_roots),
        "T_GENERATOR": transverse_generator,
        "LAMBDA_GENERATOR": longitudinal_generator,
        "ROOT_STACKED": [sp.Tuple(record[0], record[1]) for record in records],
        "ROOT_MK": [sp.Tuple(record[0], record[2]) for record in records],
        "ROOT_MT": [sp.Tuple(record[0], record[3]) for record in records],
        "ROOT_E1": [sp.Tuple(record[0], record[4]) for record in records],
        "ROOT_E2": [sp.Tuple(record[0], record[5]) for record in records],
        "ROOT_E3": [sp.Tuple(record[0], record[6]) for record in records],
        "ROOT_E4": [sp.Tuple(record[0], record[7]) for record in records],
        "ROOTS_PASSING_E1": passing_e1,
        "ROOTS_PASSING_E2": passing_e2,
        "ROOTS_PASSING_E3": passing_e3,
        "ROOTS_PASSING_E4": passing_e4,
        "ROOT_MULTISET": root_multiset,
        "ROOT_MULTIPLICITIES": [sp.Tuple(root, multiplicities[root]) for root in roots],
        "ROOT_NULLITIES": [sp.Tuple(record[0], record[8]) for record in records],
        "ROOT_SCALING_SCALED": scaled_roots,
        "ROOT_SCALING_QUADRATIC": quadratic_scaled_roots,
        "ROOT_SCALING_RESIDUAL": scaling_residuals,
        "ROOT_OMEGA2_SIGNS": [sp.Tuple(root, sign) for root, sign in zip(roots, root_signs)],
        "TRANSVERSE_ROOTS_Q": transverse_q_roots,
        "TRANSVERSE_ROOTS_Q_WAVEVECTOR_OCCURRENCES": transverse_q_wavevector_occurrences,
        "SPEED_SQUARED_CANDIDATES": speed_candidates,
        "HOMOGENEITY_DEFECTS": homogeneity_defects,
        "SPEED_Q_DERIVATIVES": speed_derivatives,
        "ROOT_Q_NUMERATOR_DEGREES": [pair[0] for pair in degree_pairs],
        "ROOT_Q_DENOMINATOR_DEGREES": [pair[1] for pair in degree_pairs],
        "DIM_SPEED_FROM_EXPRESSION": [
            sp.Tuple(speed, dimension) for speed, dimension in zip(speed_candidates, speed_dimensions)
        ],
        "DIM_SPEED_DIFFERENCE": [
            sp.Tuple(speed, sp.simplify(dimension - squared_velocity_dimension))
            for speed, dimension in zip(speed_candidates, speed_dimensions)
        ],
        "ASSUMPTIONS": assumption_object,
        "_M_A": matrix_a,
    }, sp.Integer(len(spatial_coordinates)))
    extend_produced_outputs(output, dimensions)
    guards = {
        "matrix_residual": matrix_residual,
        "root_degree": (len(root_multiset), sp.Poly(determinant, omega2).degree()),
        "root_sets": (sp.FiniteSet(*roots), sp.FiniteSet(*solved_roots)),
        "dimension_solution": dimensions["DIM_SOLUTION"],
    }
    return output, guards, None


all_outputs = {}  # DERIVED · end derivation package collection
all_guards = {}  # DERIVED · derivation guard collection
fatal_outputs = {}  # DERIVED · failed derivation output collection
main_assumptions = assumptions_for(actions["MAIN"][1], actions["MAIN"][4])  # PREMISE · stated assumption set
for control_name, action_specification in actions.items():
    assumption_object = (  # PREMISE · package assumption set
        main_assumptions
        if control_name == "MAIN"
        else assumptions_for(action_specification[1], action_specification[4])
    )
    with sp.assuming(assumption_object):
        result, guards, failure = derive(*action_specification, assumption_object)  # DERIVED · package derivation products
    if failure is not None:
        fatal_outputs[control_name] = failure  # DERIVED · failed package output
        continue
    all_outputs[control_name] = result  # DERIVED · completed package output
    all_guards[control_name] = guards  # DERIVED · completed package guards

for control_name in ("MAIN", "X6"):
    if control_name in all_outputs:
        assumption_object = all_outputs[control_name]["ASSUMPTIONS"]  # PREMISE · package assumption set
        with sp.assuming(assumption_object):
            extend_produced_outputs(
                all_outputs[control_name],
                direction_block(all_outputs[control_name]["_M_A"]),
            )

if "MAIN" in all_outputs:
    zero_wavevector_substitution = {component: 0 for component in k_input}  # CONTROL · selected wavevector locus
    zero_matrix = all_outputs["MAIN"]["_M_A"].subs(zero_wavevector_substitution).applyfunc(sp.factor)  # DERIVED · specialised dynamical matrix
    zero_wavevector = derived_wavevector.subs(zero_wavevector_substitution)  # DERIVED · specialised wavevector
    zero_test_block = evaluated_root_test_block(zero_matrix, zero_wavevector)  # DERIVED · specialised root tests
    add_produced_outputs(
        all_outputs["MAIN"],
        {f"K_ZERO_{key}": value for key, value in zero_test_block.items()},
        sp.Integer(len(spatial_coordinates)),
    )

if "X6" in all_outputs:
    anisotropic_assumption = all_outputs["X6"]["ASSUMPTIONS"]  # PREMISE · anisotropic-package assumption set
    with sp.assuming(anisotropic_assumption):
        isotropic_parameter_matrix = all_outputs["X6"]["_M_A"].subs(rho_z, rho_br).applyfunc(sp.factor)  # CONTROL · isotropic-parameter specialisation
        isotropic_parameter_test_block = evaluated_root_test_block(isotropic_parameter_matrix, derived_wavevector)  # DERIVED · specialised root tests
    add_produced_outputs(
        all_outputs["X6"],
        {f"RHO_Z_TO_RHO_BR_{key}": value for key, value in isotropic_parameter_test_block.items()},
        sp.Integer(len(spatial_coordinates)),
    )


emitted_tags = set()  # CONTROL · emitted-name uniqueness state


def contains_authored_text(value):
    if isinstance(value, str):
        return True
    if isinstance(value, dict):
        return any(contains_authored_text(item) for pair in value.items() for item in pair)
    if isinstance(value, (list, tuple, set)):
        return any(contains_authored_text(item) for item in value)
    return False


def emit(tag, value):
    full_tag = f"PY_S9_{tag}"
    assert re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", full_tag)
    assert full_tag not in emitted_tags
    assert not contains_authored_text(value)
    rendered = sp.sstr(value).replace("\n", "")
    emitted_tags.add(full_tag)
    print(f"{full_tag}: {rendered}")


def unique_item(tag, items):
    materialised = list(items)
    emit(f"{tag}_CANDIDATES_FOUND", materialised)
    assert len(materialised) == 1
    return materialised[0]


main_name_changes = {  # CONTROL · standard emission-name substitutions
    "DET_M_FACTORED": "FACTORED_DETERMINANT",
    "ROOT_MULTISET": "FULL_ROOT_MULTISET",
    "DIM_PRIMARY_INERTIA": "INERTIA_COEFFICIENT_DIMENSION",
    "DIM_PRIMARY_STIFFNESS": "STIFFNESS_COEFFICIENT_DIMENSION",
    "DIM_STIFFNESS_MINUS_INERTIA": "COEFFICIENT_DIMENSION_DIFFERENCE",
}
standard_emission_names = frozenset({  # CONTROL · standard emission-name collection
    "FACTORED_DETERMINANT",
    "FULL_ROOT_MULTISET",
    "TRANSVERSE_MULTIPLICITY",
    "TRANSVERSE_SPEED_SQUARED",
    "DISPERSION_SCALING_RESIDUAL_FLEXURAL",
    "INERTIA_COEFFICIENT_DIMENSION",
    "STIFFNESS_COEFFICIENT_DIMENSION",
    "COEFFICIENT_DIMENSION_DIFFERENCE",
    "IMPLIED_SPEED_DIMENSION",
    "SPEED_DIMENSION_DIFFERENCE",
    "DYNAMICAL_MATRIX_ROUTE_RESIDUAL",
    "BARE_FIELD_COEFFICIENT_DIMENSION",
})
posited_output_classes = {  # CONTROL · posited output class assignments
    "LAGRANGIAN": "PREMISE",
    "PLANE_WAVE_ANSATZ": "PREMISE",
    "ASSUMPTIONS": "PREMISE",
    "DIM_ENERGY_DENSITY": "PREMISE",
    "FIELD_DIMENSION": "PREMISE",
    "WAVEVECTOR_NORM_DIMENSION": "PREMISE",
    "DIM_SQUARED_VELOCITY": "PREMISE",
}
if "MAIN" in all_outputs:
    main_outputs = all_outputs["MAIN"]  # DERIVED · main-package output collection
    main_output_computation_dimensions = main_outputs[dkey_computation_dimensions]  # DERIVED · main-package production-dimension collection
    main_outputs = {  # DERIVED · standard-named main-package output collection
        main_name_changes.get(name, name): value
        for name, value in main_outputs.items()
        if name != dkey_computation_dimensions
    }
    main_outputs[dkey_computation_dimensions] = {  # DERIVED · standard-named production-dimension collection
        main_name_changes.get(name, name): computation_dimension
        for name, computation_dimension in main_output_computation_dimensions.items()
    }
    all_outputs["MAIN"] = main_outputs  # DERIVED · standard-named main-package output collection

    transverse_roots = set(main_outputs["ROOTS_PASSING_E1"])  # DERIVED · transverse-root collection
    transverse_multiplicity = unique_item(  # DERIVED · transverse-root multiplicity
        "TRANSVERSE_MULTIPLICITY",
        [
            multiplicity
            for root, multiplicity in main_outputs["ROOT_E2"]
            if root in transverse_roots
        ],
    )
    transverse_speed_squared = unique_item(  # DERIVED · transverse speed candidate
        "TRANSVERSE_SPEED_SQUARED",
        main_outputs["SPEED_SQUARED_CANDIDATES"],
    )
    route_operand_a, route_operand_b, route_residual = main_outputs["ROUTE_OPERANDS_AND_RESIDUAL"]  # DERIVED · dynamical-matrix route comparison
    del route_operand_a, route_operand_b
    implied_speed_dimensions = [dimension for _, dimension in main_outputs["DIM_SPEED_FROM_EXPRESSION"]]  # DERIVED · per-candidate unit results
    speed_dimension_differences = [dimension for _, dimension in main_outputs["DIM_SPEED_DIFFERENCE"]]  # DERIVED · per-candidate unit residuals
    add_produced_outputs(main_outputs, {
        "TRANSVERSE_MULTIPLICITY": transverse_multiplicity,
        "TRANSVERSE_SPEED_SQUARED": transverse_speed_squared,
        "DYNAMICAL_MATRIX_ROUTE_RESIDUAL": route_residual,
        "IMPLIED_SPEED_DIMENSION": implied_speed_dimensions,
        "SPEED_DIMENSION_DIFFERENCE": speed_dimension_differences,
    }, sp.Integer(len(spatial_coordinates)))

if "X7" in all_outputs:
    all_outputs["X7"]["DISPERSION_SCALING_RESIDUAL_FLEXURAL"] = unique_item(  # DERIVED · flexural scaling residual
        "DISPERSION_SCALING_RESIDUAL_FLEXURAL",
        all_outputs["X7"]["ROOT_SCALING_RESIDUAL"]
    )

if "X8" in all_outputs:
    all_outputs["X8"]["BARE_FIELD_COEFFICIENT_DIMENSION"] = unique_item(  # DERIVED · bare-field coefficient unit result
        "BARE_FIELD_COEFFICIENT_DIMENSION",
        [
            dimension
            for coefficient, dimension in all_outputs["X8"]["DIM_COEFFICIENTS"]
            if coefficient == mu_G
        ],
    )


for control_name in actions:
    if control_name in all_outputs:
        for suffix, value in all_outputs[control_name].items():
            if not suffix.startswith("_"):
                emitted_name = suffix if suffix in standard_emission_names else f"{control_name}_{suffix}"  # CONTROL · emitted object name
                emit(emitted_name, value)
    if control_name in fatal_outputs:
        for suffix, value in fatal_outputs[control_name].items():
            emit(f"{control_name}_{suffix}", value)


# Guards deliberately follow emission of both operands and their residuals.
for control_name, guards in all_guards.items():
    assert guards["matrix_residual"] == sp.zeros(3)
    assert guards["root_degree"][0] == guards["root_degree"][1]
    assert guards["root_sets"][0] == guards["root_sets"][1]
    assert guards["dimension_solution"]

if fatal_outputs:
    raise SystemExit(1)


def reconstruction_expression(value):
    return sp.srepr(value)


def exact_reconstruction_match(live_value, reconstructed_value):
    return type(live_value) is type(reconstructed_value) and sp.srepr(live_value) == sp.srepr(reconstructed_value)


def generated_record_lines(name, value, class_tag, dimension=None):
    lines = [
        f"    {name!r}: {{",
        f"        'display': {sp.sstr(value)!r},",
        f"        'value': _restore({reconstruction_expression(value)!r}),",
    ]
    if dimension is not None:
        lines.append(f"        'dim': _restore({reconstruction_expression(dimension)!r}),")
    lines.extend([
        f"        'class': {class_tag!r},",
        "        'step': 'S9',",
        "    },",
    ])
    return lines


def generated_ledger_key(name, computation_dimension):
    base_name = re.sub(r"_d[0-9]+$", "", name)
    if computation_dimension == D:
        return base_name
    return f"{base_name}_d{computation_dimension}"


def d_partition(records):
    grouped = {}
    for name, computation_dimension in records:
        grouped.setdefault(computation_dimension, []).append(sp.Symbol(name))  # DERIVED · D-partition member collection
    return sp.Tuple(*(
        sp.Tuple(
            computation_dimension,
            sp.Tuple(*sorted(names, key=sp.default_sort_key)),
        )
        for computation_dimension, names in sorted(
            grouped.items(),
            key=lambda item: sp.default_sort_key(item[0]),
        )
    ))


def write_exports():
    main_outputs = all_outputs["MAIN"]
    main_output_computation_dimensions = main_outputs[dkey_computation_dimensions]
    assert set(posited_output_classes) <= set(main_outputs)
    classes = declaration_classes(Path(__file__))
    exported_inputs = [
        (name, globals()[name], class_tag)
        for name, class_tag in classes.items()
        if class_tag in ("KNOB", "STRUCTURAL")
    ]
    computed_dimensions = [
        *main_outputs["DIM_COEFFICIENTS"],
        *main_outputs["DIM_SPEED_FROM_EXPRESSION"],
    ]

    def computed_dimension_for(value):
        for dimensioned_object, dimension in computed_dimensions:
            if exact_reconstruction_match(value, dimensioned_object):
                return dimension
        return None

    source_lines = [
        "# S9_exports.py — GENERATED by S9_light_requires_shear_sympy_audit.py. Do not edit.",
        "import sympy as sp",
        "",
        "",
        "def _restore(source):",
        "    return eval(source, {'__builtins__': {}, **vars(sp)})",
        "",
        "",
        "LEDGER = {",
    ]
    live_records = []
    main_production_records = []
    for name, value, class_tag in exported_inputs:
        dimension = computed_dimension_for(value)
        live_records.append((name, value, class_tag, dimension))
        source_lines.extend(generated_record_lines(name, value, class_tag, dimension))
    for suffix, value in main_outputs.items():
        if suffix.startswith("_"):
            continue
        base_name = suffix.lower()
        if suffix not in main_output_computation_dimensions:
            raise RuntimeError(f"unclassified export production path: {base_name}")
        computation_dimension = main_output_computation_dimensions[suffix]
        name = generated_ledger_key(base_name, computation_dimension)
        class_tag = posited_output_classes.get(suffix, "DERIVED")
        dimension = computed_dimension_for(value)
        live_records.append((name, value, class_tag, dimension))
        main_production_records.append((name, computation_dimension))
        source_lines.extend(generated_record_lines(name, value, class_tag, dimension))
    source_lines.append("}")
    source_lines.append("")
    export_source = "\n".join(source_lines)
    export_path = Path(__file__).with_name("S9_exports.py")
    export_path.write_text(export_source, encoding="utf-8")

    reconstructed_namespace = {}
    exec(compile(export_path.read_text(encoding="utf-8"), str(export_path), "exec"), reconstructed_namespace)
    reconstructed_ledger = reconstructed_namespace["LEDGER"]
    live_count = len(live_records)
    residuals = []
    for name, value, class_tag, dimension in live_records:
        record = reconstructed_ledger[name]
        residuals.append(sp.Integer(not exact_reconstruction_match(value, record["value"])))
        if dimension is not None:
            residuals.append(sp.Integer(not exact_reconstruction_match(dimension, record["dim"])))
        residuals.append(sp.Integer(class_tag != record["class"]))
        residuals.append(sp.Integer("S9" != record["step"]))
    roundtrip_residual = sum(residuals, sp.Integer(0))
    roundtrip_count_residual = sp.Integer(live_count - len(reconstructed_ledger))
    computed_dimension_count = sum("dim" in record for record in reconstructed_ledger.values())
    absent_dimension_count = len(reconstructed_ledger) - computed_dimension_count
    class_tally = sp.Tuple(*(
        sp.Tuple(
            sp.Symbol(class_tag),  # STRUCTURAL · exported class label
            sp.Integer(sum(record["class"] == class_tag for record in reconstructed_ledger.values())),
        )
        for class_tag in CLASS_TAGS
    ))
    class_tally_residual = sp.Integer(len(reconstructed_ledger)) - sum(
        (count for _, count in class_tally),
        sp.Integer(0),
    )
    production_d_partition = d_partition(main_production_records)
    emit("EXPORT_ROUNDTRIP_LIVE_COUNT", sp.Integer(live_count))
    emit("EXPORT_ROUNDTRIP_RECONSTRUCTED_COUNT", sp.Integer(len(reconstructed_ledger)))
    emit("EXPORT_ROUNDTRIP_COUNT_RESIDUAL", roundtrip_count_residual)
    emit("EXPORT_ROUNDTRIP_RESIDUAL", roundtrip_residual)
    emit("EXPORT_COMPUTED_DIMENSION_COUNT", sp.Integer(computed_dimension_count))
    emit("EXPORT_ABSENT_DIMENSION_COUNT", sp.Integer(absent_dimension_count))
    emit("EXPORT_CLASS_TALLY", class_tally)
    emit("EXPORT_CLASS_TALLY_RESIDUAL", class_tally_residual)
    emit("EXPORT_D_PARTITION", production_d_partition)
    assert roundtrip_count_residual == 0
    assert roundtrip_residual == 0
    assert class_tally_residual == 0


write_exports()
