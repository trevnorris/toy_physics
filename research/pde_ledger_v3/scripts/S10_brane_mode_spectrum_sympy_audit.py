#!/usr/bin/env python3
"""S10: derive and count the modes of MacCullagh's curl-only brane action."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import sympy as sp


HERE = Path(__file__).resolve().parent
REDUCTION_DIR = HERE.parent / "reduction"
sys.path.insert(0, str(REDUCTION_DIR))

# This is deliberately the ledger's shared reader, not a private YAML parser.
from registry_read import Registry, load_registry  # noqa: E402


REQUIRED_QIDS = (
    "Q.brane.rho_br",
    "Q.brane.mu_R",
    "Q.brane.D_brane",
)
D_CASES = (2, 3, 4, 5)


@dataclass(frozen=True)
class Mode:
    root: sp.Expr
    rank: int
    nullity: int
    basis: tuple[sp.Matrix, ...]
    dot_products: tuple[sp.Expr, ...]
    orientations: tuple[str, ...]


@dataclass(frozen=True)
class Spectrum:
    label: str
    dimension: int
    lagrangian: sp.Expr
    equations: tuple[sp.Expr, ...]
    matrix_generic: sp.Matrix
    wavevector: sp.Matrix
    determinant: sp.Expr
    modes: tuple[Mode, ...]


class Checks:
    def __init__(self) -> None:
        self.ok = True

    def record(self, name: str, condition: bool, detail: str = "") -> bool:
        passed = bool(condition)
        self.ok = self.ok and passed
        suffix = f" {detail}" if detail else ""
        print(f"S10_CHECK_{name}: {'PASS' if passed else 'FAIL'}{suffix}")
        return passed


def rendered(expression: object) -> str:
    if isinstance(expression, sp.MatrixBase):
        return sp.sstr(expression.applyfunc(sp.factor))
    if isinstance(expression, sp.Basic):
        return sp.sstr(sp.factor(expression))
    return str(expression)


def vector_rendered(vector: sp.Matrix) -> str:
    return "(" + ",".join(rendered(value) for value in vector) + ")"


def representative_direction(dimension: int) -> sp.Matrix:
    """Explicit non-axis unit directions; no modal-count pattern is encoded here."""
    directions = {
        2: (sp.Rational(3, 5), sp.Rational(4, 5)),
        3: (sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(2, 3)),
        4: (sp.Rational(1, 2),) * 4,
        5: (
            sp.Rational(1, 3),
            sp.Rational(1, 3),
            sp.Rational(1, 3),
            sp.Rational(1, 3),
            sp.sqrt(5) / 3,
        ),
    }
    return sp.Matrix(directions[dimension])


def derive_lagrangian_and_eom(
    dimension: int,
    rho_br: sp.Symbol,
    stiffness_coefficient: sp.Expr,
    stiffness_form: str,
) -> tuple[sp.Expr, tuple[sp.Expr, ...], tuple[sp.Expr, ...], tuple[sp.Symbol, ...]]:
    """Construct the density and take its Euler--Lagrange variation componentwise."""
    t = sp.Symbol(f"t_D{dimension}", real=True)
    coordinates = sp.symbols(f"x0:{dimension}_D{dimension}", real=True)
    fields = tuple(
        sp.Function(f"u{component}_D{dimension}")(t, *coordinates)
        for component in range(dimension)
    )
    gradients = tuple(
        tuple(sp.diff(fields[component], coordinates[direction]) for component in range(dimension))
        for direction in range(dimension)
    )
    kinetic_square = sum(sp.diff(field, t) ** 2 for field in fields)

    if stiffness_form == "antisymmetric":
        stiffness_square = sp.Rational(1, 2) * sum(
            (gradients[i][j] - gradients[j][i]) ** 2
            for i in range(dimension)
            for j in range(dimension)
        )
    elif stiffness_form == "full-gradient":
        stiffness_square = sum(
            gradients[i][j] ** 2
            for i in range(dimension)
            for j in range(dimension)
        )
    else:
        raise ValueError(f"unknown stiffness form: {stiffness_form}")

    lagrangian = sp.expand(
        rho_br * kinetic_square / 2 - stiffness_coefficient * stiffness_square / 2
    )
    equations = []
    for field in fields:
        varied = sp.diff(sp.diff(lagrangian, sp.diff(field, t)), t)
        varied += sum(
            sp.diff(sp.diff(lagrangian, sp.diff(field, coordinate)), coordinate)
            for coordinate in coordinates
        )
        equations.append(sp.simplify(varied))
    return lagrangian, tuple(equations), fields, coordinates


def plane_wave_matrix(
    equations: tuple[sp.Expr, ...],
    fields: tuple[sp.Expr, ...],
    coordinates: tuple[sp.Symbol, ...],
    dimension: int,
) -> tuple[sp.Matrix, tuple[sp.Symbol, ...], tuple[sp.Symbol, ...], sp.Symbol]:
    """Substitute u=a exp(i(k.x-omega t)); the common exponential cancels."""
    t = next(iter(fields[0].args))
    amplitudes = sp.symbols(f"a0:{dimension}_D{dimension}", real=True)
    wave_numbers = sp.symbols(f"k0:{dimension}_D{dimension}", real=True)
    omega_squared = sp.Symbol(f"omega_squared_D{dimension}", real=True)

    substitutions: dict[sp.Expr, sp.Expr] = {}
    for component, field in enumerate(fields):
        substitutions[sp.diff(field, t, 2)] = -omega_squared * amplitudes[component]
        for left, coordinate_left in enumerate(coordinates):
            for right, coordinate_right in enumerate(coordinates):
                substitutions[sp.diff(field, coordinate_left, coordinate_right)] = (
                    -wave_numbers[left] * wave_numbers[right] * amplitudes[component]
                )

    plane_equations = tuple(
        sp.expand(equation.subs(substitutions, simultaneous=True)) for equation in equations
    )
    matrix = sp.Matrix(plane_equations).jacobian(amplitudes).applyfunc(sp.simplify)
    return matrix, amplitudes, wave_numbers, omega_squared


def classify_orientation(vector: sp.Matrix, wavevector: sp.Matrix) -> tuple[sp.Expr, str]:
    """Classify with dot products (the Gram equality detects parallel vectors)."""
    dot_product = sp.simplify(wavevector.dot(vector))
    gram_residual = sp.simplify(
        dot_product**2 - wavevector.dot(wavevector) * vector.dot(vector)
    )
    if dot_product == 0:
        orientation = "perpendicular"
    elif gram_residual == 0:
        orientation = "parallel"
    else:
        orientation = "neither"
    return dot_product, orientation


def roots_are_complete(
    determinant: sp.Expr, omega_squared: sp.Symbol, roots: Iterable[sp.Expr]
) -> bool:
    """Compare the square-free determinant with the product over distinct roots."""
    polynomial = sp.Poly(determinant, omega_squared)
    square_free = polynomial.sqf_part().monic().as_expr()
    root_product = sp.Poly(
        sp.prod(omega_squared - root for root in roots), omega_squared
    ).monic().as_expr()
    return sp.simplify(square_free - root_product) == 0


def compute_spectrum(
    label: str,
    dimension: int,
    rho_br: sp.Symbol,
    stiffness_coefficient: sp.Expr,
    stiffness_form: str,
    kappa: sp.Symbol,
    checks: Checks,
) -> Spectrum:
    lagrangian, equations, fields, coordinates = derive_lagrangian_and_eom(
        dimension, rho_br, stiffness_coefficient, stiffness_form
    )
    matrix_generic, amplitudes, wave_numbers, omega_squared = plane_wave_matrix(
        equations, fields, coordinates, dimension
    )
    no_unsubstituted_fields = all(
        not entry.has(sp.Derivative) and not any(entry.has(field) for field in fields)
        for entry in matrix_generic
    )
    checks.record(f"{label}_PLANE_WAVE_SUBSTITUTION", no_unsubstituted_fields)
    checks.record(
        f"{label}_LINEAR_AMPLITUDE_SYSTEM",
        all(sp.diff(entry, amplitude) == 0 for entry in matrix_generic for amplitude in amplitudes),
    )

    direction = representative_direction(dimension)
    checks.record(
        f"{label}_UNIT_DIRECTION",
        sp.simplify(direction.dot(direction) - 1) == 0,
        f"direction={vector_rendered(direction)}",
    )
    wavevector = kappa * direction
    specialized = matrix_generic.subs(dict(zip(wave_numbers, wavevector))).applyfunc(sp.simplify)
    determinant = sp.factor(specialized.det())
    roots = tuple(
        sorted(
            set(sp.solve(sp.Eq(determinant, 0), omega_squared)),
            key=sp.default_sort_key,
        )
    )
    checks.record(
        f"{label}_DISTINCT_ROOT_COVERAGE",
        bool(roots) and roots_are_complete(determinant, omega_squared, roots),
        f"root_count={len(roots)}",
    )

    modes = []
    for root_index, root in enumerate(roots):
        root_matrix = specialized.subs(omega_squared, root).applyfunc(sp.simplify)
        rank = int(root_matrix.rank())
        nullity = dimension - rank
        basis = tuple(vector.applyfunc(sp.simplify) for vector in root_matrix.nullspace())
        checks.record(
            f"{label}_ROOT{root_index}_RANK_NULLSPACE",
            len(basis) == nullity
            and all(
                all(sp.simplify(value) == 0 for value in root_matrix * vector)
                for vector in basis
            ),
            f"rank={rank} nullity={nullity} basis_size={len(basis)}",
        )
        classified = tuple(classify_orientation(vector, wavevector) for vector in basis)
        dot_products = tuple(item[0] for item in classified)
        orientations = tuple(item[1] for item in classified)
        modes.append(Mode(root, rank, nullity, basis, dot_products, orientations))

    checks.record(
        f"{label}_NULLITY_SUM",
        sum(mode.nullity for mode in modes) == dimension,
        f"observed={sum(mode.nullity for mode in modes)} D={dimension}",
    )
    return Spectrum(
        label,
        dimension,
        lagrangian,
        equations,
        matrix_generic,
        wavevector,
        determinant,
        tuple(modes),
    )


def print_spectrum(spectrum: Spectrum, prefix: str) -> None:
    print(
        f"S10_{prefix}_ROOTS: "
        + ", ".join(rendered(mode.root) for mode in spectrum.modes)
    )
    for root_index, mode in enumerate(spectrum.modes):
        print(
            f"S10_{prefix}_ROOT_{root_index}: omega_squared={rendered(mode.root)} "
            f"rank={mode.rank} nullity={mode.nullity}"
        )
        for basis_index, (vector, dot_product, orientation) in enumerate(
            zip(mode.basis, mode.dot_products, mode.orientations)
        ):
            print(
                f"S10_{prefix}_ROOT_{root_index}_ORIENTATION_{basis_index}: "
                f"basis={vector_rendered(vector)} dot_k={rendered(dot_product)} "
                f"class={orientation}"
            )
    print(
        f"S10_{prefix}_NULLITY_SUM: "
        f"{sum(mode.nullity for mode in spectrum.modes)}"
    )


def dimension_subtract(
    left: tuple[sp.Expr, sp.Expr, sp.Expr],
    right: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return tuple(sp.simplify(a - b) for a, b in zip(left, right))  # type: ignore[return-value]


def dimension_scale(
    coefficient: int, dimension: tuple[sp.Expr, sp.Expr, sp.Expr]
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return tuple(sp.simplify(coefficient * value) for value in dimension)  # type: ignore[return-value]


def derive_coefficient_dimensions(
    D_brane: sp.Symbol,
) -> tuple[
    tuple[sp.Expr, sp.Expr, sp.Expr],
    tuple[sp.Expr, sp.Expr, sp.Expr],
    tuple[sp.Expr, sp.Expr, sp.Expr],
]:
    energy = (sp.Integer(2), sp.Integer(-2), sp.Integer(1))
    brane_volume = (D_brane, sp.Integer(0), sp.Integer(0))
    lagrangian_density = dimension_subtract(energy, brane_volume)
    displacement = (sp.Integer(1), sp.Integer(0), sp.Integer(0))
    inverse_time = (sp.Integer(0), sp.Integer(-1), sp.Integer(0))
    inverse_length = (sp.Integer(-1), sp.Integer(0), sp.Integer(0))
    time_derivative = tuple(
        sp.simplify(a + b) for a, b in zip(displacement, inverse_time)
    )
    spatial_derivative = tuple(
        sp.simplify(a + b) for a, b in zip(displacement, inverse_length)
    )
    rho_dimension = dimension_subtract(
        lagrangian_density, dimension_scale(2, time_derivative)  # type: ignore[arg-type]
    )
    mu_dimension = dimension_subtract(
        lagrangian_density, dimension_scale(2, spatial_derivative)  # type: ignore[arg-type]
    )
    return lagrangian_density, rho_dimension, mu_dimension


def registry_inputs(
    registry: Registry, checks: Checks
) -> tuple[sp.Symbol, sp.Symbol, sp.Symbol] | None:
    missing = tuple(qid for qid in REQUIRED_QIDS if qid not in registry.quantities)
    if not checks.record(
        "REGISTRY_REQUIRED_QUANTITIES",
        not missing,
        f"missing={','.join(missing) if missing else 'none'}",
    ):
        return None
    rho_br, mu_R, D_brane = (registry.symbols[qid] for qid in REQUIRED_QIDS)
    checks.record(
        "REGISTRY_SYMBOL_IDENTITY",
        all(
            str(symbol) == registry.quantities[qid].symbol_name
            for qid, symbol in zip(REQUIRED_QIDS, (rho_br, mu_R, D_brane))
        ),
    )
    return rho_br, mu_R, D_brane


def spectrum_signature(spectrum: Spectrum) -> tuple[tuple[sp.Expr, int, tuple[str, ...]], ...]:
    return tuple(
        (mode.root, mode.nullity, mode.orientations) for mode in spectrum.modes
    )


def main() -> int:
    checks = Checks()
    registry = load_registry(REDUCTION_DIR)
    inputs = registry_inputs(registry, checks)
    if inputs is None:
        print("S10_VERDICT: FAIL")
        return 1
    rho_br, mu_R, D_brane = inputs

    d_quantity = registry.quantities["Q.brane.D_brane"]
    checks.record(
        "D_BRANE_DISCRETE_AXIS",
        d_quantity.kind == "discrete-choice"
        and d_quantity.counting_axis == "discrete-structural",
        f"kind={d_quantity.kind} axis={d_quantity.counting_axis}",
    )
    checks.record(
        "D_BRANE_DIMENSIONLESS",
        d_quantity.dimension == (0, 0, 0),
        f"declared={d_quantity.dimension}",
    )
    if not checks.record(
        "D_BRANE_INTEGER_VALUE",
        d_quantity.value is not None and bool(d_quantity.value.is_Integer),
        f"declared={d_quantity.value}",
    ):
        print("S10_VERDICT: FAIL")
        return 1
    registry_dimension = int(d_quantity.value)

    lagrangian_dimension, rho_derived, mu_derived = derive_coefficient_dimensions(D_brane)
    rho_at_registry_d = tuple(sp.simplify(value.subs(D_brane, registry_dimension)) for value in rho_derived)
    mu_at_registry_d = tuple(sp.simplify(value.subs(D_brane, registry_dimension)) for value in mu_derived)
    rho_declared = registry.quantities["Q.brane.rho_br"].dimension
    mu_declared = registry.quantities["Q.brane.mu_R"].dimension
    print(
        "S10_DIMENSION_DERIVATION: "
        f"energy_density_D={rendered(sp.Tuple(*lagrangian_dimension))} "
        f"rho_br={rendered(sp.Tuple(*rho_derived))} "
        f"mu_R={rendered(sp.Tuple(*mu_derived))}"
    )
    checks.record(
        "RHO_BR_REGISTRY_DIMENSION",
        rho_at_registry_d == rho_declared,
        f"derived_at_D={rho_at_registry_d} declared={rho_declared}",
    )
    checks.record(
        "MU_R_REGISTRY_DIMENSION",
        mu_at_registry_d == mu_declared,
        f"derived_at_D={mu_at_registry_d} declared={mu_declared}",
    )

    gradient = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"G{i}{j}", real=True))
    antisymmetric_square = sp.Rational(1, 2) * sum(
        (gradient[i, j] - gradient[j, i]) ** 2 for i in range(3) for j in range(3)
    )
    ordinary_curl = sp.Matrix(
        (
            gradient[1, 2] - gradient[2, 1],
            gradient[2, 0] - gradient[0, 2],
            gradient[0, 1] - gradient[1, 0],
        )
    )
    ordinary_curl_square = sp.expand(ordinary_curl.dot(ordinary_curl))
    curl_residual = sp.simplify(antisymmetric_square - ordinary_curl_square)
    print(
        "S10_CURL_REDUCTION_D3: "
        f"antisymmetric={rendered(antisymmetric_square)} "
        f"ordinary_curl_norm={rendered(ordinary_curl_square)} "
        f"residual={rendered(curl_residual)}"
    )
    checks.record("CURL_REDUCTION_D3", curl_residual == 0)

    kappa = sp.Symbol("kappa", positive=True, real=True)
    main_spectra: dict[int, Spectrum] = {}
    for dimension in D_CASES:
        spectrum = compute_spectrum(
            f"MAIN_D{dimension}",
            dimension,
            rho_br,
            mu_R,
            "antisymmetric",
            kappa,
            checks,
        )
        main_spectra[dimension] = spectrum
        print_spectrum(spectrum, f"D{dimension}")
        if dimension == registry_dimension:
            print(
                f"S10_EOM_D{dimension}: "
                + " ; ".join(rendered(equation) + " = 0" for equation in spectrum.equations)
            )
            print(
                f"S10_DYNAMICAL_MATRIX_D{dimension}: "
                f"{rendered(spectrum.matrix_generic)}"
            )

    print("S10_D_TABLE:")
    for dimension in D_CASES:
        entries = []
        for mode in main_spectra[dimension].modes:
            classes = "/".join(sorted(set(mode.orientations)))
            entries.append(
                f"root={rendered(mode.root)};nullity={mode.nullity};orientation={classes}"
            )
        print(f"S10_D_TABLE_ROW: D={dimension} " + " | ".join(entries))

    baseline_d3 = main_spectra[3]
    form_control = compute_spectrum(
        "FORM_CONTROL_D3",
        3,
        rho_br,
        mu_R,
        "full-gradient",
        kappa,
        checks,
    )
    print("S10_FORM_CONTROL: stiffness=full-gradient-squared D=3")
    print_spectrum(form_control, "FORM_CONTROL_D3")
    print(
        "S10_FORM_CONTROL_EOM_D3: "
        + " ; ".join(rendered(equation) + " = 0" for equation in form_control.equations)
    )
    checks.record(
        "FORM_CONTROL_SENSITIVITY",
        spectrum_signature(form_control) != spectrum_signature(baseline_d3),
        "control_signature_compared_to_main_D3",
    )

    coefficient_control = compute_spectrum(
        "COEFFICIENT_CONTROL_D3",
        3,
        rho_br,
        2 * mu_R,
        "antisymmetric",
        kappa,
        checks,
    )
    print("S10_COEFFICIENT_CONTROL: stiffness_coefficient=2*mu_R D=3")
    print(
        "S10_COEFFICIENT_CONTROL_ROOTS: "
        + ", ".join(rendered(mode.root) for mode in coefficient_control.modes)
    )
    checks.record(
        "COEFFICIENT_CONTROL_ARITHMETIC_SENSITIVITY",
        tuple(mode.root for mode in coefficient_control.modes)
        != tuple(mode.root for mode in baseline_d3.modes),
        "control_roots_compared_to_main_D3",
    )

    print(f"S10_VERDICT: {'PASS' if checks.ok else 'FAIL'}")
    return 0 if checks.ok else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # Every unexpected failure is still a named audit check.
        print(
            "S10_CHECK_RUNTIME: FAIL "
            f"exception={type(exc).__name__} diagnostic={str(exc)}"
        )
        print("S10_VERDICT: FAIL")
        raise SystemExit(1)
