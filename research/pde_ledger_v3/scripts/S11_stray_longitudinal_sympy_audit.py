#!/usr/bin/env python3
"""S11: independent SymPy audit of the curl-plus-compression brane action."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
REDUCTION_DIR = HERE.parent / "reduction"
sys.path.insert(0, str(REDUCTION_DIR))

# The shared reader is part of the audit: A5 and A7 compare against transported data.
from registry_read import Registry, load_registry  # noqa: E402


D_CASES = (2, 3, 4, 5)


def emit(tag: str, text: object) -> None:
    print(f"S11_{tag}: {text}")


def clean(expression: sp.Expr) -> sp.Expr:
    return sp.factor(sp.simplify(expression))


def rendered(expression: object) -> str:
    if isinstance(expression, sp.MatrixBase):
        return "Matrix(" + sp.sstr(expression.applyfunc(clean).tolist()) + ")"
    if isinstance(expression, sp.Basic):
        return sp.sstr(clean(expression))
    return str(expression)


def vector_rendered(vector: sp.Matrix) -> str:
    return "(" + ",".join(rendered(value) for value in vector) + ")"


def is_zero_matrix(matrix: sp.MatrixBase) -> bool:
    return all(clean(value) == 0 for value in matrix)


@dataclass(frozen=True)
class AssertionResult:
    name: str
    passed: bool
    detail: str


class Checks:
    def __init__(self) -> None:
        self.results: list[AssertionResult] = []

    def record(self, name: str, condition: bool, detail: str) -> bool:
        passed = bool(condition)
        self.results.append(AssertionResult(name, passed, detail))
        return passed


@dataclass(frozen=True)
class InvariantCensus:
    dimension: int
    n_so: int
    n_o: int
    standard_span: int
    constructed_span: int
    standard_so_invariant: bool
    standard_reflection_even: bool
    extra_so_residual: bool
    extra_reflection_odd: bool


def representation_generator(dimension: int, left: int, right: int) -> sp.Matrix:
    """Infinitesimal action on G -> R G R.T, in row-major coordinates."""
    generator = sp.zeros(dimension)
    generator[left, right] = 1
    generator[right, left] = -1
    size = dimension**2
    representation = sp.zeros(size)
    for i in range(dimension):
        for j in range(dimension):
            output = i * dimension + j
            for source in range(dimension):
                representation[output, source * dimension + j] += generator[i, source]
                representation[output, i * dimension + source] += generator[j, source]
    return representation


def quadratic_matrix(polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Matrix:
    """Return C for the homogeneous quadratic polynomial g.T*C*g."""
    return sp.hessian(polynomial, variables) / 2


def upper_vector(matrix: sp.Matrix, pairs: tuple[tuple[int, int], ...]) -> sp.Matrix:
    return sp.Matrix([matrix[row, column] for row, column in pairs])


def invariant_census(dimension: int) -> InvariantCensus:
    """Construct every invariant symmetric bilinear form by Lie-algebra nullity."""
    size = dimension**2
    pairs = tuple((row, column) for row in range(size) for column in range(row, size))
    coefficients = sp.symbols(f"c0:{len(pairs)}_D{dimension}")
    candidate = sp.MutableSparseMatrix(size, size, {})
    for coefficient, (row, column) in zip(coefficients, pairs):
        candidate[row, column] = coefficient
        candidate[column, row] = coefficient

    equations: list[sp.Expr] = []
    generators: list[sp.Matrix] = []
    for left in range(dimension):
        for right in range(left + 1, dimension):
            representation = representation_generator(dimension, left, right)
            generators.append(representation)
            residual = representation.T * candidate + candidate * representation
            equations.extend(
                residual[row, column]
                for row, column in pairs
                if residual[row, column] != 0
            )
    linear_system, _ = sp.linear_eq_to_matrix(equations, coefficients)
    so_basis = tuple(linear_system.nullspace())
    n_so = len(so_basis)

    # SO(D) plus one reflection generates O(D).  Restrict that reflection to
    # the just-computed SO(D)-invariant nullspace rather than solving twice.
    coordinate_signs = tuple(-1 if index == 0 else 1 for index in range(dimension))
    representation_signs = tuple(
        coordinate_signs[i] * coordinate_signs[j]
        for i in range(dimension)
        for j in range(dimension)
    )
    reflection_odd_rows = tuple(
        index
        for index, (row, column) in enumerate(pairs)
        if representation_signs[row] * representation_signs[column] == -1
    )
    reflection_restriction = sp.Matrix(
        [
            [basis[index] for basis in so_basis]
            for index in reflection_odd_rows
        ]
    )
    n_o = n_so - reflection_restriction.rank()

    entries = sp.symbols(f"g0:{size}_D{dimension}")
    gradient = sp.Matrix(dimension, dimension, entries)
    trace = sp.trace(gradient)
    standard_polynomials = (
        trace**2,
        sum(gradient[i, j] ** 2 for i in range(dimension) for j in range(dimension)),
        sum(
            gradient[i, j] * gradient[j, i]
            for i in range(dimension)
            for j in range(dimension)
        ),
    )
    standard_matrices = tuple(
        quadratic_matrix(polynomial, entries) for polynomial in standard_polynomials
    )
    standard_span = sp.Matrix.hstack(
        *(upper_vector(matrix, pairs) for matrix in standard_matrices)
    ).rank()
    reflection = sp.diag(*representation_signs)
    standard_so_invariant = all(
        is_zero_matrix(generator.T * matrix + matrix * generator)
        for matrix in standard_matrices
        for generator in generators
    )
    standard_reflection_even = all(
        is_zero_matrix(reflection.T * matrix * reflection - matrix)
        for matrix in standard_matrices
    )

    extra_matrix: sp.Matrix | None = None
    if dimension == 2:
        pseudoscalar = sum(
            sp.LeviCivita(i, j) * gradient[i, j]
            for i in range(dimension)
            for j in range(dimension)
        )
        extra_matrix = quadratic_matrix(trace * pseudoscalar, entries)
    elif dimension == 4:
        epsilon_contraction = sum(
            sp.LeviCivita(i, j, k, ell) * gradient[i, j] * gradient[k, ell]
            for i in range(dimension)
            for j in range(dimension)
            for k in range(dimension)
            for ell in range(dimension)
        )
        extra_matrix = quadratic_matrix(epsilon_contraction, entries)

    constructed_matrices = standard_matrices + (() if extra_matrix is None else (extra_matrix,))
    constructed_span = sp.Matrix.hstack(
        *(upper_vector(matrix, pairs) for matrix in constructed_matrices)
    ).rank()
    if extra_matrix is None:
        extra_so_residual = True
        extra_reflection_odd = True
    else:
        extra_so_residual = all(
            is_zero_matrix(generator.T * extra_matrix + extra_matrix * generator)
            for generator in generators
        )
        reflected_extra = reflection.T * extra_matrix * reflection
        extra_reflection_odd = is_zero_matrix(reflected_extra + extra_matrix)

    return InvariantCensus(
        dimension,
        n_so,
        n_o,
        standard_span,
        constructed_span,
        standard_so_invariant,
        standard_reflection_even,
        extra_so_residual,
        extra_reflection_odd,
    )


@dataclass(frozen=True)
class FieldDerivation:
    dimension: int
    lagrangian: sp.Expr
    eom: tuple[sp.Expr, ...]
    eom_formula_residuals: tuple[sp.Expr, ...] | None
    matrix: sp.Matrix
    wavevector: sp.Matrix


def field_derivation(
    dimension: int,
    rho_br: sp.Symbol,
    mu_r: sp.Symbol,
    b_comp: sp.Symbol,
    omega2: sp.Symbol,
    *,
    form_control: bool = False,
    mu_br: sp.Symbol | None = None,
) -> FieldDerivation:
    """Build L, vary it componentwise, and extract the plane-wave matrix."""
    t = sp.Symbol(f"t_D{dimension}", real=True)
    coordinates = sp.symbols(f"x0:{dimension}_D{dimension}", real=True)
    fields = tuple(
        sp.Function(f"u{component}_D{dimension}")(t, *coordinates)
        for component in range(dimension)
    )
    gradient = tuple(
        tuple(sp.diff(fields[component], coordinates[direction]) for component in range(dimension))
        for direction in range(dimension)
    )
    divergence = sum(gradient[index][index] for index in range(dimension))
    curl2 = sp.Rational(1, 2) * sum(
        (gradient[i][j] - gradient[j][i]) ** 2
        for i in range(dimension)
        for j in range(dimension)
    )
    kinetic = sp.Rational(1, 2) * rho_br * sum(
        sp.diff(field, t) ** 2 for field in fields
    )
    if form_control:
        if mu_br is None:
            raise ValueError("form control requires mu_br")
        symmetric_traceless = tuple(
            tuple(
                (gradient[i][j] + gradient[j][i]) / 2
                - (divergence / dimension if i == j else 0)
                for j in range(dimension)
            )
            for i in range(dimension)
        )
        replacement = mu_br * sum(
            symmetric_traceless[i][j] ** 2
            for i in range(dimension)
            for j in range(dimension)
        )
    else:
        replacement = sp.Rational(1, 2) * b_comp * divergence**2
    lagrangian = kinetic - sp.Rational(1, 2) * mu_r * curl2 - replacement

    eom = tuple(
        clean(
            sp.diff(sp.diff(lagrangian, sp.diff(field, t)), t)
            + sum(
                sp.diff(
                    sp.diff(lagrangian, sp.diff(field, coordinates[direction])),
                    coordinates[direction],
                )
                for direction in range(dimension)
            )
            - sp.diff(lagrangian, field)
        )
        for field in fields
    )
    if form_control:
        formula_residuals = None
    else:
        expected_eom = tuple(
            rho_br * sp.diff(fields[component], t, 2)
            - mu_r
            * sum(
                sp.diff(fields[component], coordinate, 2)
                for coordinate in coordinates
            )
            + (mu_r - b_comp) * sp.diff(divergence, coordinates[component])
            for component in range(dimension)
        )
        formula_residuals = tuple(
            clean(actual - expected)
            for actual, expected in zip(eom, expected_eom)
        )

    amplitudes = sp.symbols(f"a0:{dimension}_D{dimension}")
    wavevector_symbols = sp.symbols(f"k0:{dimension}_D{dimension}", real=True)
    wavevector = sp.Matrix(wavevector_symbols)
    substitutions: dict[sp.Derivative, sp.Expr] = {}
    field_indices = {field: index for index, field in enumerate(fields)}
    coordinate_indices = {coordinate: index for index, coordinate in enumerate(coordinates)}
    for equation in eom:
        for derivative in equation.atoms(sp.Derivative):
            factor: sp.Expr = amplitudes[field_indices[derivative.expr]]
            for variable, count in derivative.variable_count:
                if variable == t:
                    if count != 2:
                        raise ValueError(f"unexpected time derivative: {derivative}")
                    factor *= -omega2
                else:
                    factor *= (sp.I * wavevector_symbols[coordinate_indices[variable]]) ** count
            substitutions[derivative] = factor
    plane_eom = sp.Matrix(
        [sp.expand(equation.xreplace(substitutions)) for equation in eom]
    )
    amplitude_vector = sp.Matrix(amplitudes)
    stiffness_action = plane_eom + rho_br * omega2 * amplitude_vector
    matrix = stiffness_action.jacobian(amplitude_vector).applyfunc(clean)
    return FieldDerivation(
        dimension,
        lagrangian,
        eom,
        formula_residuals,
        matrix,
        wavevector,
    )


def representative_direction(dimension: int) -> sp.Matrix:
    directions = {
        2: (sp.Rational(3, 5), sp.Rational(4, 5)),
        3: (sp.Rational(2, 3), sp.Rational(1, 3), sp.Rational(2, 3)),
        4: (sp.Rational(1, 2),) * 4,
        5: (sp.Rational(2, 5),) * 4 + (sp.Rational(3, 5),),
    }
    direction = sp.Matrix(directions[dimension])
    if clean(direction.dot(direction) - 1) != 0:
        raise ValueError("representative direction is not a unit vector")
    return direction


@dataclass(frozen=True)
class Mode:
    root: sp.Expr
    nullity: int
    basis: tuple[sp.Matrix, ...]
    orientations: tuple[str, ...]
    dot_residuals: tuple[sp.Expr, ...]
    parallel_residuals: tuple[sp.Matrix, ...]


def spectrum_by_nullity(
    matrix: sp.Matrix,
    wavevector: sp.Matrix,
    dimension: int,
    rho_br: sp.Symbol,
    omega2: sp.Symbol,
    kappa: sp.Symbol,
) -> tuple[Mode, ...]:
    direction = representative_direction(dimension)
    representative = matrix.subs(
        dict(zip(tuple(wavevector), tuple(kappa * direction)))
    ).applyfunc(clean)
    characteristic = clean((representative - rho_br * omega2 * sp.eye(dimension)).det())
    roots = tuple(dict.fromkeys(sp.solve(characteristic, omega2)))
    modes: list[Mode] = []
    for root in roots:
        kernel_matrix = (representative - rho_br * root * sp.eye(dimension)).applyfunc(clean)
        basis = tuple(kernel_matrix.nullspace())
        orientations: list[str] = []
        dots: list[sp.Expr] = []
        parallel_residuals: list[sp.Matrix] = []
        for vector in basis:
            dot = clean(direction.dot(vector))
            parallel_residual = (vector - direction * direction.dot(vector)).applyfunc(clean)
            dots.append(dot)
            parallel_residuals.append(parallel_residual)
            if dot == 0:
                orientations.append("perpendicular")
            elif is_zero_matrix(parallel_residual):
                orientations.append("parallel")
            else:
                orientations.append("mixed")
        modes.append(
            Mode(
                clean(root),
                len(basis),
                basis,
                tuple(orientations),
                tuple(dots),
                tuple(parallel_residuals),
            )
        )
    return tuple(modes)


def mode_by_orientation(modes: tuple[Mode, ...], orientation: str) -> Mode:
    matches = tuple(
        mode for mode in modes if mode.orientations and set(mode.orientations) == {orientation}
    )
    if len(matches) != 1:
        raise ValueError(f"expected one {orientation} mode, observed {len(matches)}")
    return matches[0]


def dimension_text(vector: sp.Matrix) -> str:
    return "[" + ",".join(rendered(value) for value in vector) + "]"


def audit() -> int:
    checks = Checks()
    rho_br, mu_r, b_comp, mu_br = sp.symbols(
        "rho_br mu_R B_comp mu_br", positive=True
    )
    omega2 = sp.Symbol("omega2", real=True)
    kappa = sp.Symbol("kappa", positive=True)

    # A0: solve the full invariant-tensor problem, then display explicit bases.
    censuses: dict[int, InvariantCensus] = {}
    for dimension in D_CASES:
        census = invariant_census(dimension)
        censuses[dimension] = census
        emit(
            f"A0_D{dimension}",
            f"N_SO={census.n_so} N_O={census.n_o} "
            f"standard_span={census.standard_span} constructed_span={census.constructed_span} "
            f"standard_SO_invariant={census.standard_so_invariant} "
            f"standard_reflection_even={census.standard_reflection_even}",
        )
        expected_constructed = census.n_so
        checks.record(
            f"A0_D{dimension}_complete_construction",
            census.constructed_span == expected_constructed
            and census.standard_span == census.n_o
            and census.standard_so_invariant
            and census.standard_reflection_even,
            f"N_SO={census.n_so},N_O={census.n_o}",
        )
    emit(
        "A0_EXTRA_D2",
        "invariant=(tr G)*Sum(epsilon_ij*G_ij) "
        f"SO_residual_zero={censuses[2].extra_so_residual} "
        f"reflection_sign=-1:{censuses[2].extra_reflection_odd}",
    )
    emit(
        "A0_EXTRA_D4",
        "invariant=Sum(epsilon_ijkl*G_ij*G_kl) "
        f"SO_residual_zero={censuses[4].extra_so_residual} "
        f"reflection_sign=-1:{censuses[4].extra_reflection_odd}",
    )
    checks.record(
        "A0_pseudoscalars",
        all(
            censuses[d].extra_so_residual and censuses[d].extra_reflection_odd
            for d in (2, 4)
        ),
        "D=2 and D=4 extras are SO-even and reflection-odd",
    )

    # A1: form the field Lagrangian, vary it, and only then extract M.
    derivations = {
        dimension: field_derivation(
            dimension, rho_br, mu_r, b_comp, omega2
        )
        for dimension in D_CASES
    }
    eom_residuals = {
        dimension: tuple(derivation.eom_formula_residuals or ())
        for dimension, derivation in derivations.items()
    }
    emit(
        "A1_EULER_LAGRANGE",
        "rho_br*u_j,tt-mu_R*Delta(u_j)+(mu_R-B_comp)*partial_j(div u)=0 "
        f"computed_residuals={{{','.join(f'D{d}:{rendered(res)}' for d, res in eom_residuals.items())}}}",
    )
    checks.record(
        "A1_euler_lagrange",
        all(all(value == 0 for value in residuals) for residuals in eom_residuals.values()),
        "componentwise field variation agrees with displayed equation",
    )
    matrix_residuals: dict[int, sp.Matrix] = {}
    for dimension, derivation in derivations.items():
        k_squared = derivation.wavevector.dot(derivation.wavevector)
        expected_matrix = (
            mu_r * k_squared * sp.eye(dimension)
            + (b_comp - mu_r) * derivation.wavevector * derivation.wavevector.T
        )
        matrix_residuals[dimension] = (derivation.matrix - expected_matrix).applyfunc(clean)
    emit(
        "A1_DYNAMICAL_MATRIX_FORM",
        "M_ij=mu_R*(k.k)*delta_ij+(B_comp-mu_R)*k_i*k_j "
        f"residual_zero_by_D={{{','.join(f'D{d}:{is_zero_matrix(r)}' for d, r in matrix_residuals.items())}}}",
    )
    emit("A1_M_D3", rendered(derivations[3].matrix))
    checks.record(
        "A1_dynamical_matrix",
        all(is_zero_matrix(residual) for residual in matrix_residuals.values()),
        "all D=2..5 matrices equal the displayed tensor form",
    )

    spectra: dict[int, tuple[Mode, ...]] = {}
    for dimension, derivation in derivations.items():
        modes = spectrum_by_nullity(
            derivation.matrix,
            derivation.wavevector,
            dimension,
            rho_br,
            omega2,
            kappa,
        )
        spectra[dimension] = modes
        emit(f"A1_D{dimension}_DISTINCT_ROOTS", len(modes))
        for mode in modes:
            orientation_label = (
                next(iter(set(mode.orientations)))
                if len(set(mode.orientations)) == 1
                else "mixed"
            )
            basis_text = ";".join(vector_rendered(vector) for vector in mode.basis)
            emit(
                f"A1_D{dimension}_{orientation_label.upper()}",
                f"omega2={rendered(mode.root)} nullity={mode.nullity} "
                f"kernel_orientations={','.join(mode.orientations)} basis={basis_text}",
            )
        perpendicular = mode_by_orientation(modes, "perpendicular")
        parallel = mode_by_orientation(modes, "parallel")
        checks.record(
            f"A1_D{dimension}_spectrum_by_nullity",
            len(modes) == 2
            and perpendicular.nullity == dimension - 1
            and parallel.nullity == 1
            and all(dot == 0 for dot in perpendicular.dot_residuals)
            and all(is_zero_matrix(res) for res in parallel.parallel_residuals),
            f"nullities={perpendicular.nullity}+{parallel.nullity}",
        )

    perpendicular_root = mode_by_orientation(spectra[3], "perpendicular").root
    parallel_root = mode_by_orientation(spectra[3], "parallel").root
    perpendicular_speed = clean(sp.sqrt(perpendicular_root / kappa**2))
    longitudinal_speed = clean(sp.sqrt(parallel_root / kappa**2))
    speed_residuals = (
        clean(perpendicular_speed**2 - perpendicular_root / kappa**2),
        clean(longitudinal_speed**2 - parallel_root / kappa**2),
    )
    emit(
        "A1_PHASE_SPEEDS",
        f"c_perpendicular={rendered(perpendicular_speed)} "
        f"c_L={rendered(longitudinal_speed)} residuals={rendered(speed_residuals)}",
    )
    checks.record(
        "A1_phase_speeds",
        speed_residuals == (0, 0),
        "both displayed speeds square to their computed omega2/kappa2 roots",
    )

    registry: Registry = load_registry()
    registry_speed_residuals: dict[str, sp.Expr] = {}
    registry_speed_derivations = (
        (
            "R4",
            "c_gamma",
            perpendicular_speed,
            {
                "Q.brane.c_gamma": perpendicular_speed,
                "Q.brane.mu_R": mu_r,
                "Q.brane.rho_br": rho_br,
            },
        ),
        (
            "R5",
            "c_L",
            longitudinal_speed,
            {
                "Q.brane.c_L": longitudinal_speed,
                "Q.brane.B_comp": b_comp,
                "Q.brane.rho_br": rho_br,
            },
        ),
    )
    registry_speed_details: list[str] = []
    for relation_id, speed_name, derived_speed, substitutions in registry_speed_derivations:
        registry_residual = registry.require_admitted(relation_id).residual
        if registry_residual is None:
            raise ValueError(f"{relation_id} has no registry residual")
        residual = clean(
            registry_residual.subs(
                {
                    registry.symbols[qid]: value
                    for qid, value in substitutions.items()
                }
            )
        )
        registry_speed_residuals[relation_id] = residual
        registry_speed_details.append(
            f"{relation_id}_derived_{speed_name}={rendered(derived_speed)} "
            f"{relation_id}_registry_residual={rendered(residual)}"
        )
    emit("A1_REGISTRY_PHASE_SPEEDS", " ".join(registry_speed_details))

    # A2: differentiate the roots with respect to the cross-sector coefficient.
    perpendicular_b_residual = clean(sp.diff(perpendicular_root, b_comp))
    parallel_mu_residual = clean(sp.diff(parallel_root, mu_r))
    emit(
        "A2_PERP_DEPENDS_ON_B_COMP",
        f"{perpendicular_b_residual != 0} derivative_residual={rendered(perpendicular_b_residual)}",
    )
    emit(
        "A2_PARALLEL_DEPENDS_ON_MU_R",
        f"{parallel_mu_residual != 0} derivative_residual={rendered(parallel_mu_residual)}",
    )
    checks.record(
        "A2_cross_sector_residuals",
        perpendicular_b_residual == 0 and parallel_mu_residual == 0,
        "both cross derivatives vanish",
    )

    # A3: solve equality of the two computed roots and recompute the kernel there.
    degeneracy_solution = sp.solve(
        sp.Eq(perpendicular_root, parallel_root), b_comp, dict=True
    )
    emit("A3_DEGENERACY_LOCUS", rendered(degeneracy_solution))
    degeneracy_ok = degeneracy_solution == [{b_comp: mu_r}]
    for dimension, derivation in derivations.items():
        direction = representative_direction(dimension)
        representative = derivation.matrix.subs(
            dict(zip(tuple(derivation.wavevector), tuple(kappa * direction)))
        )
        locus_matrix = (
            representative.subs(b_comp, mu_r)
            - rho_br
            * perpendicular_root.subs(b_comp, mu_r)
            * sp.eye(dimension)
        ).applyfunc(clean)
        full_kernel = tuple(locus_matrix.nullspace())
        perpendicular_subspace = tuple(direction.T.nullspace())
        parallel_residual = locus_matrix * direction
        perpendicular_residuals = tuple(locus_matrix * vector for vector in perpendicular_subspace)
        emit(
            f"A3_D{dimension}_ON_LOCUS",
            f"nullity={len(full_kernel)} parallel_dim=1 "
            f"perpendicular_dim={len(perpendicular_subspace)}",
        )
        degeneracy_ok = (
            degeneracy_ok
            and len(full_kernel) == dimension
            and is_zero_matrix(parallel_residual)
            and all(is_zero_matrix(value) for value in perpendicular_residuals)
        )
    checks.record(
        "A3_degeneracy",
        degeneracy_ok,
        "B_comp=mu_R gives a D-dimensional kernel with 1 parallel and D-1 perpendicular",
    )

    # A4: dimensions follow from energy/brane-volume and the displacement premise.
    d_symbol = sp.Symbol("D", integer=True, positive=True)
    length = sp.Matrix([1, 0, 0])
    time = sp.Matrix([0, 1, 0])
    energy = sp.Matrix([2, -2, 1])
    energy_density = energy - d_symbol * length
    displacement = length
    velocity = displacement - time
    gradient_dimension = displacement - length
    rho_dimension = (energy_density - 2 * velocity).applyfunc(clean)
    modulus_dimension = (energy_density - 2 * gradient_dimension).applyfunc(clean)
    phase_speed_dimension = ((modulus_dimension - rho_dimension) / 2).applyfunc(clean)
    emit(
        "A4_U_USE",
        f"[u]={dimension_text(displacement)} used in [partial_t u]={dimension_text(velocity)} "
        f"and [partial_i u]={dimension_text(gradient_dimension)}",
    )
    emit("A4_RHO_BR_D", dimension_text(rho_dimension))
    emit("A4_MU_R_D", dimension_text(modulus_dimension))
    emit("A4_B_COMP_D", dimension_text(modulus_dimension))
    emit("A4_PERP_PHASE_SPEED_D", dimension_text(phase_speed_dimension))
    emit("A4_PARALLEL_PHASE_SPEED_D", dimension_text(phase_speed_dimension))
    derived_d3 = {
        "rho_br": tuple(int(value.subs(d_symbol, 3)) for value in rho_dimension),
        "mu_R": tuple(int(value.subs(d_symbol, 3)) for value in modulus_dimension),
        "B_comp": tuple(int(value.subs(d_symbol, 3)) for value in modulus_dimension),
        "c_perp": tuple(int(value) for value in phase_speed_dimension),
        "c_L": tuple(int(value) for value in phase_speed_dimension),
    }
    emit(
        "A4_D3",
        " ".join(f"{name}={list(value)}" for name, value in derived_d3.items()),
    )
    checks.record(
        "A4_dimensions",
        rho_dimension == sp.Matrix([-d_symbol, 0, 1])
        and modulus_dimension == sp.Matrix([2 - d_symbol, -2, 1])
        and phase_speed_dimension == sp.Matrix([1, -1, 0]),
        "closed-D dimensions and both phase-speed dimensions derived",
    )

    # A5: external registry comparison (no adjustment of the derivation).
    registry_dimension_results: list[bool] = []
    for name, qid in (
        ("rho_br", "Q.brane.rho_br"),
        ("mu_R", "Q.brane.mu_R"),
    ):
        declared = registry.quantities[qid].dimension
        derived = derived_d3[name]
        passed = declared == derived
        registry_dimension_results.append(passed)
        emit(
            f"A5_{name.upper()}",
            f"{'PASS' if passed else 'FAIL'} derived={list(derived)} declared={list(declared)}",
        )
    checks.record(
        "A5_registry_dimensions",
        all(registry_dimension_results),
        "rho_br and mu_R D=3 vectors match independently derived values",
    )

    # A6: solve only the bulk dispersion relation and state the exact scope.
    c = sp.Symbol("c", positive=True)
    c_s0 = sp.Symbol("c_s0", positive=True)
    k_in = sp.Symbol("k_in", positive=True)
    k_w2 = sp.Symbol("k_w2", real=True)
    matching_solution = sp.solve(
        sp.Eq(c**2 * k_in**2, c_s0**2 * (k_in**2 + k_w2)),
        k_w2,
    )[0]
    real_condition = sp.solve_univariate_inequality(matching_solution > 0, c)
    imaginary_condition = sp.solve_univariate_inequality(matching_solution < 0, c)
    threshold_condition = sp.solve_univariate_inequality(
        sp.Eq(matching_solution, 0), c
    )
    emit("A6_KW2", f"k_w2={rendered(matching_solution)}")
    emit(
        "A6_REAL",
        f"condition={real_condition}; real k_w, propagating scalar bulk channel exists, "
        "so spatial binding is not kinematically protected",
    )
    emit(
        "A6_IMAGINARY",
        f"condition={imaginary_condition}; imaginary k_w, no propagating scalar channel; "
        "normal evanescence is compatible with binding but does not prove a bound mode",
    )
    emit(
        "A6_THRESHOLD",
        f"condition={threshold_condition}; k_w=0, grazing threshold with no normal decay, "
        "hence not spatially localized by this calculation",
    )
    emit(
        "A6_SCOPE",
        "kinematic phase matching only; no brane-bulk coupling law, interface boundary "
        "condition, or energy-flux calculation is present",
    )
    emit(
        "A6_NECESSARY",
        "the dispersion algebra determines whether a propagating scalar channel exists "
        "at the shared omega and k_in",
    )
    emit(
        "A6_REQUIRES_COUPLING_LAW",
        "whether an existing channel is used and at what rate (and whether an evanescent "
        "solution forms a bound eigenmode) is not established",
    )
    perpendicular_matching = matching_solution.subs(
        c**2, perpendicular_root / kappa**2
    ).subs(k_in, kappa)
    emit(
        "A6_PERPENDICULAR_ATTEMPT",
        f"formal_scalar_dispersion_gives_k_w2={rendered(perpendicular_matching)}; "
        "the scalar bulk mode has no transverse vector polarization to match, and without "
        "a coupling law this algebra neither establishes nor forbids coupling",
    )
    matching_residual = clean(
        c**2 * k_in**2
        - c_s0**2 * (k_in**2 + matching_solution)
    )
    checks.record(
        "A6_phase_matching",
        matching_residual == 0
        and str(real_condition) == "c_s0 < c"
        and str(imaginary_condition) == "c < c_s0"
        and str(threshold_condition) == "Eq(c, c_s0)",
        "solved residual is zero and the three positive-speed cases are exhaustive",
    )

    # A7: compare the inequality's two sides using registry dimensions only.
    c_l_registry_dimension = registry.quantities["Q.brane.c_L"].dimension
    c_s0_registry_dimension = registry.quantities["Q.medium.c_s0"].dimension
    commensurable = c_l_registry_dimension == c_s0_registry_dimension
    emit(
        "A7_COMMENSURABILITY",
        f"{'PASS' if commensurable else 'FAIL'} "
        f"[c_L]={list(c_l_registry_dimension)} [c_s0]={list(c_s0_registry_dimension)}",
    )
    checks.record(
        "A7_registry_commensurability",
        commensurable,
        "registry-imported c_L and c_s0 dimensions agree",
    )

    # A8: replace the compression invariant by the symmetric-traceless form.
    form_derivation = field_derivation(
        3,
        rho_br,
        mu_r,
        b_comp,
        omega2,
        form_control=True,
        mu_br=mu_br,
    )
    form_modes = spectrum_by_nullity(
        form_derivation.matrix,
        form_derivation.wavevector,
        3,
        rho_br,
        omega2,
        kappa,
    )
    form_perpendicular = mode_by_orientation(form_modes, "perpendicular")
    form_parallel = mode_by_orientation(form_modes, "parallel")
    form_movements = {
        "perpendicular": clean(form_perpendicular.root - perpendicular_root),
        "parallel": clean(form_parallel.root - parallel_root),
    }
    form_cross_residuals = {
        "perpendicular_wrt_mu_br": clean(sp.diff(form_perpendicular.root, mu_br)),
        "parallel_wrt_mu_R": clean(sp.diff(form_parallel.root, mu_r)),
    }
    for label, mode in (
        ("PERPENDICULAR", form_perpendicular),
        ("PARALLEL", form_parallel),
    ):
        emit(
            f"A8_{label}",
            f"omega2={rendered(mode.root)} nullity={mode.nullity} "
            f"kernel_orientations={','.join(mode.orientations)}",
        )
    emit(
        "A8_MOVEMENT",
        " ".join(
            f"{name}=MOVED(residual={rendered(residual)})"
            if residual != 0
            else f"{name}=UNCHANGED(residual=0)"
            for name, residual in form_movements.items()
        ),
    )
    emit(
        "A8_CROSS_DEPENDENCE",
        f"perpendicular_depends_on_mu_br={form_cross_residuals['perpendicular_wrt_mu_br'] != 0} "
        f"derivative={rendered(form_cross_residuals['perpendicular_wrt_mu_br'])} "
        f"parallel_depends_on_mu_R={form_cross_residuals['parallel_wrt_mu_R'] != 0} "
        f"derivative={rendered(form_cross_residuals['parallel_wrt_mu_R'])}",
    )
    checks.record(
        "A8_form_control",
        form_perpendicular.nullity == 2
        and form_parallel.nullity == 1
        and all(residual != 0 for residual in form_movements.values())
        and form_cross_residuals["perpendicular_wrt_mu_br"] != 0
        and form_cross_residuals["parallel_wrt_mu_R"] == 0,
        "both generic roots move; perpendicular gains mu_br dependence, parallel stays mu_R-independent",
    )

    # A9: change only the coefficient multiplying the original compression form.
    coefficient_matrix = derivations[3].matrix.subs(b_comp, 2 * b_comp).applyfunc(clean)
    coefficient_modes = spectrum_by_nullity(
        coefficient_matrix,
        derivations[3].wavevector,
        3,
        rho_br,
        omega2,
        kappa,
    )
    coefficient_perpendicular = mode_by_orientation(coefficient_modes, "perpendicular")
    coefficient_parallel = mode_by_orientation(coefficient_modes, "parallel")
    coefficient_movements = {
        "perpendicular": clean(coefficient_perpendicular.root - perpendicular_root),
        "parallel": clean(coefficient_parallel.root - parallel_root),
    }
    for label, mode in (
        ("PERPENDICULAR", coefficient_perpendicular),
        ("PARALLEL", coefficient_parallel),
    ):
        emit(
            f"A9_{label}",
            f"omega2={rendered(mode.root)} nullity={mode.nullity} "
            f"kernel_orientations={','.join(mode.orientations)}",
        )
    emit(
        "A9_MOVEMENT",
        " ".join(
            f"{name}=MOVED(residual={rendered(residual)})"
            if residual != 0
            else f"{name}=UNCHANGED(residual=0)"
            for name, residual in coefficient_movements.items()
        ),
    )
    checks.record(
        "A9_coefficient_control",
        coefficient_perpendicular.nullity == 2
        and coefficient_parallel.nullity == 1
        and coefficient_movements["perpendicular"] == 0
        and coefficient_movements["parallel"] != 0,
        "only the parallel root moves under B_comp -> 2 B_comp",
    )

    checks.record(
        "A1_registry_phase_speed_relations",
        all(residual == 0 for residual in registry_speed_residuals.values()),
        ",".join(
            f"{relation_id}_registry_residual={rendered(residual)}"
            for relation_id, residual in registry_speed_residuals.items()
        ),
    )

    # A10: the verdict is exactly the conjunction printed here, not a physics verdict.
    emit(
        "A10_ASSERTION_LIST",
        ",".join(f"{index:02d}={result.name}" for index, result in enumerate(checks.results, 1)),
    )
    for index, result in enumerate(checks.results, 1):
        emit(
            f"ASSERTION_{index:02d}",
            f"{'PASS' if result.passed else 'FAIL'} {result.name} -- {result.detail}",
        )
    verdict = all(result.passed for result in checks.results)
    emit(
        "A10_SCOPE",
        "the verdict certifies only that these SymPy checks are mutually consistent; "
        "it is not a verdict on the physics",
    )
    emit("VERDICT", "PASS" if verdict else "FAIL")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(audit())
