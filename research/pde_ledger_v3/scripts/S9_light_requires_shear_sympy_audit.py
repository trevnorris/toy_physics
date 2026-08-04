#!/usr/bin/env python3
"""Independent SymPy audit for ledger step S9."""

import re

import sympy as sp


# Field, coordinate, amplitude, and material symbols.
t, x, y, z = sp.symbols("t x y z", real=True)
coordinates = (t, x, y, z)
spatial_coordinates = (x, y, z)
u_functions = tuple(sp.Function(name)(t, x, y, z) for name in ("u1", "u2", "u3"))
a_symbols = sp.symbols("a1 a2 a3", real=True)
b_symbols = sp.symbols("b1 b2 b3", real=True)
omega = sp.symbols("omega", real=True)
omega2 = sp.symbols("omega2", real=True)
k_input = sp.symbols("kx ky kz", real=True)
rho_br, mu_R, rho_z, mu_F, mu_G = sp.symbols("rho_br mu_R rho_z mu_F mu_G", positive=True)
lambda_rho, lambda_mu = sp.symbols("lambda_rho lambda_mu", positive=True)
D = sp.symbols("D", integer=True, positive=True)
q = sp.symbols("q", positive=True)


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
identity3 = sp.eye(3)
actions = {
    "MAIN": (construct_curl_action(rho_br * identity3, mu_R), (rho_br, mu_R), rho_br, mu_R, ()),
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
position = sp.Matrix(spatial_coordinates)
input_wavevector = sp.Matrix(k_input)
phase_exponent = sp.I * ((input_wavevector.T * position)[0] - omega * t)
phase = sp.exp(phase_exponent)
plane_wave = sp.Matrix([amplitude * phase for amplitude in a_symbols])

# Recover wave data by differentiating the constructed phase rather than retyping its combinations.
derived_wavevector = sp.Matrix([-sp.I * sp.diff(phase, coordinate) / phase for coordinate in spatial_coordinates])
wave_norm = sp.factor((derived_wavevector.T * derived_wavevector)[0])
transverse_generator = sp.simplify(wave_norm * identity3 - derived_wavevector * derived_wavevector.T)
longitudinal_generator = sp.simplify(derived_wavevector * derived_wavevector.T)


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


length_dimension = sp.Matrix([1, 0, 0])
time_dimension = sp.Matrix([0, 1, 0])
mass_dimension = sp.Matrix([0, 0, 1])
acceleration_dimension = length_dimension - 2 * time_dimension
force_dimension = mass_dimension + acceleration_dimension
energy_dimension = force_dimension + length_dimension
energy_density_dimension = energy_dimension - D * length_dimension
field_dimension = length_dimension
q_dimension = -2 * length_dimension
squared_velocity_dimension = 2 * (length_dimension - time_dimension)


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
            [sp.Symbol(f"dim_{coefficient}_{axis}") for axis in ("L", "T", "M")]
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
    common_output = {
        "DIM_ENERGY_DENSITY": energy_density_dimension,
        "DIM_TERMS": terms,
        "DIM_FIELD_MULTIORDERS": orders_by_term,
        "DIM_TERM_EXPRESSIONS": term_dimensions,
        "DIM_LINEAR_SYSTEM": equations,
        "DIM_SOLUTION": solutions,
    }
    if not solutions:
        common_output["_DIMENSION_SOLVE_FAILED"] = sp.S.true
        return common_output
    solution = solutions[0]
    coefficient_dimensions = [
        sp.Tuple(coefficient, sp.Matrix([solution[entry] for entry in dimension_variables[coefficient]]))
        for coefficient in dimensional_coefficients
    ]
    inertia_dimension = sp.Matrix([solution[entry] for entry in dimension_variables[primary_inertia]])
    stiffness_dimension = sp.Matrix([solution[entry] for entry in dimension_variables[primary_stiffness]])
    common_output.update({
        "DIM_COEFFICIENTS": coefficient_dimensions,
        "DIM_PRIMARY_INERTIA": inertia_dimension,
        "DIM_PRIMARY_STIFFNESS": stiffness_dimension,
        "DIM_STIFFNESS_MINUS_INERTIA": sp.simplify(stiffness_dimension - inertia_dimension),
        "DIM_SOLUTION_D3": [{key: value.subs(D, 3) for key, value in item.items()} for item in solutions],
        "DIM_PRIMARY_INERTIA_D3": inertia_dimension.subs(D, 3),
        "DIM_PRIMARY_STIFFNESS_D3": stiffness_dimension.subs(D, 3),
        "DIM_STIFFNESS_MINUS_INERTIA_D3": (stiffness_dimension - inertia_dimension).subs(D, 3),
        "_DIMENSION_LOOKUP": {
            coefficient: sp.Matrix([solution[entry] for entry in dimension_variables[coefficient]])
            for coefficient in dimensional_coefficients
        },
    })
    return common_output


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
    return output


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
    output = {
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
        "TRANSVERSE_ROOTS_Q": transverse_q_roots,
        "SPEED_SQUARED_CANDIDATES": speed_candidates,
        "HOMOGENEITY_DEFECTS": homogeneity_defects,
        "SPEED_Q_DERIVATIVES": speed_derivatives,
        "ROOT_Q_NUMERATOR_DEGREES": [pair[0] for pair in degree_pairs],
        "ROOT_Q_DENOMINATOR_DEGREES": [pair[1] for pair in degree_pairs],
        "DIM_SPEED_FROM_EXPRESSION": [
            sp.Tuple(speed, dimension) for speed, dimension in zip(speed_candidates, speed_dimensions)
        ],
        "DIM_SQUARED_VELOCITY": squared_velocity_dimension,
        "DIM_SPEED_DIFFERENCE": [
            sp.Tuple(speed, sp.simplify(dimension - squared_velocity_dimension))
            for speed, dimension in zip(speed_candidates, speed_dimensions)
        ],
        "ASSUMPTIONS": assumption_object,
        "_M_A": matrix_a,
    }
    output.update({key: value for key, value in dimensions.items() if not key.startswith("_")})
    guards = {
        "matrix_residual": matrix_residual,
        "root_degree": (len(root_multiset), sp.Poly(determinant, omega2).degree()),
        "root_sets": (sp.FiniteSet(*roots), sp.FiniteSet(*solved_roots)),
        "dimension_solution": dimensions["DIM_SOLUTION"],
    }
    return output, guards, None


all_outputs = {}
all_guards = {}
fatal_output = None
fatal_control = None
for control_name, action_specification in actions.items():
    assumption_object = assumptions_for(action_specification[1], action_specification[4])
    with sp.assuming(assumption_object):
        result, guards, failure = derive(*action_specification, assumption_object)
    if failure is not None:
        fatal_control = control_name
        fatal_output = failure
        break
    all_outputs[control_name] = result
    all_guards[control_name] = guards

if fatal_output is None:
    for control_name in ("MAIN", "X6"):
        assumption_object = all_outputs[control_name]["ASSUMPTIONS"]
        with sp.assuming(assumption_object):
            all_outputs[control_name].update(direction_block(all_outputs[control_name]["_M_A"]))

    zero_wavevector_substitution = {component: 0 for component in k_input}
    zero_matrix = all_outputs["MAIN"]["_M_A"].subs(zero_wavevector_substitution).applyfunc(sp.factor)
    zero_wavevector = derived_wavevector.subs(zero_wavevector_substitution)
    zero_test_block = evaluated_root_test_block(zero_matrix, zero_wavevector)
    all_outputs["MAIN"].update({f"K_ZERO_{key}": value for key, value in zero_test_block.items()})

    anisotropic_assumption = all_outputs["X6"]["ASSUMPTIONS"]
    with sp.assuming(anisotropic_assumption):
        isotropic_parameter_matrix = all_outputs["X6"]["_M_A"].subs(rho_z, rho_br).applyfunc(sp.factor)
        isotropic_parameter_test_block = evaluated_root_test_block(isotropic_parameter_matrix, derived_wavevector)
    all_outputs["X6"].update(
        {f"RHO_Z_TO_RHO_BR_{key}": value for key, value in isotropic_parameter_test_block.items()}
    )


emitted_tags = set()


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


if fatal_output is not None:
    for suffix, value in fatal_output.items():
        emit(f"{fatal_control}_{suffix}", value)
    raise SystemExit(1)


for control_name in actions:
    for suffix, value in all_outputs[control_name].items():
        if not suffix.startswith("_"):
            emit(f"{control_name}_{suffix}", value)


# Guards deliberately follow emission of both operands and their residuals.
for control_name, guards in all_guards.items():
    assert guards["matrix_residual"] == sp.zeros(3)
    assert guards["root_degree"][0] == guards["root_degree"][1]
    assert guards["root_sets"][0] == guards["root_sets"][1]
    assert guards["dimension_solution"]
