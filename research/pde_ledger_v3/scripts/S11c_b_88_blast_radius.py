#!/usr/bin/env python3
"""Compute the S11c-b #88 retained-grade strong-row disturbance objects.

This is a compute-and-print instrument.  It imports the committed S11c-b
SymPy engine for its candidate enumeration, emitted LAB_HELD density, and
differential-jet conventions, then constructs only the specified basis
completion and Hessian-retaining spatial derivative.
"""

from __future__ import annotations

from collections.abc import Callable
from itertools import product
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.matrices import DomainMatrix


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import S11c_b_brane_operator_sympy_audit as M  # noqa: E402


DX = Callable[[sp.Expr, int], sp.Expr]
SOURCES = (
    ("W_BG", M.grad_W),
    ("MU_R_BG", M.grad_mu),
)


def emit(tag: str, payload: object) -> object:
    """Print one tagged CAS object in the engine's serialization."""
    payload = M.casify(payload)
    print(f"{tag}: {M.render(payload)}", flush=True)
    return payload


def termcount(expression: sp.Expr) -> int:
    """Count expanded additive terms, assigning zero a count of zero."""
    expression = sp.expand(expression)
    return 0 if expression == 0 else len(sp.Add.make_args(expression))


def el_amplitude(
    density: sp.Expr,
    field: sp.Symbol,
    first_jets: tuple[sp.Symbol, ...],
    derivative: DX,
) -> sp.Expr:
    """Stored-energy Euler--Lagrange amplitude for one strong row."""
    algebraic = sp.diff(density, field)
    divergence = sp.Add(
        *(
            derivative(sp.diff(density, first_jets[direction]), direction)
            for direction in M.DIRECTIONS
        )
    )
    return sp.expand(algebraic - divergence)


def named_form(balance: sp.Tuple, name: str) -> object:
    for form_name, value in balance:
        if str(form_name) == name:
            return value
    raise KeyError(name)


def grade_support(expression: sp.Expr) -> sp.Tuple:
    polynomial = sp.Poly(sp.expand(expression), M.eta_bg, M.sigma_W)
    support = []
    for eta_power, sigma_power in product((0, 1), repeat=2):
        coefficient = polynomial.coeff_monomial(
            M.eta_bg**eta_power * M.sigma_W**sigma_power
        )
        if coefficient != 0:
            support.append(sp.Tuple(eta_power, sigma_power))
    return sp.Tuple(*support)


def monomial_expression(
    generators: tuple[sp.Symbol, ...], monomial: tuple[int, ...]
) -> sp.Expr:
    return sp.Mul(
        *(generator**power for generator, power in zip(generators, monomial))
    )


def constant_span_data(
    f_templates: tuple[sp.Expr, ...],
    r_templates: tuple[sp.Expr, ...],
    forbidden_coefficients: set[sp.Symbol],
) -> dict[str, object]:
    """Build and rank F and F+R over background-independent constants.

    DOF coordinates and all profile values/jets/bookkeepers are polynomial
    generators.  Every other symbol is admitted only as a scalar-field
    generator after checking that no free energy coefficient remains there.
    """
    templates = tuple(sp.expand(item) for item in f_templates + r_templates)
    appearing = set().union(*(item.free_symbols for item in templates)) if templates else set()
    dof_atoms = {
        atom
        for atom in appearing
        if M.DECLARED_SYMBOLS.get(atom, {}).get("class") == "COORDINATE"
    }
    background_atoms = {
        atom
        for atom in appearing
        if atom in {M.eta_bg, M.sigma_W}
        or atom.name.startswith("w1_profile")
        or atom.name.startswith("m1_profile")
    }
    joint_generators = tuple(
        sorted(dof_atoms | background_atoms, key=sp.default_sort_key)
    )
    scalar_generators = tuple(
        sorted(appearing - set(joint_generators), key=sp.default_sort_key)
    )
    leaked_coefficients = set(scalar_generators) & forbidden_coefficients
    if leaked_coefficients:
        leaked = ", ".join(str(item) for item in sorted(leaked_coefficients, key=sp.default_sort_key))
        raise RuntimeError(f"free coefficient in constant-span matrix entry: {leaked}")

    domain = sp.QQ.frac_field(*scalar_generators) if scalar_generators else sp.QQ
    polynomials = tuple(
        sp.Poly(item, *joint_generators, domain=domain) for item in templates
    )
    monomials = tuple(
        sorted({monomial for polynomial in polynomials for monomial in polynomial.monoms()})
    )
    rows = [
        [polynomial.coeff_monomial(monomial) for polynomial in polynomials]
        for monomial in monomials
    ]
    if not rows:
        dim_f = 0
        dim_f_plus_r = 0
    else:
        f_width = len(f_templates)
        f_rows = [row[:f_width] for row in rows]
        dim_f = DomainMatrix.from_list(f_rows, domain).rank() if f_width else 0
        dim_f_plus_r = DomainMatrix.from_list(rows, domain).rank()

    return {
        "JOINT_GENERATORS": joint_generators,
        "JOINT_MONOMIALS": tuple(
            monomial_expression(joint_generators, monomial) for monomial in monomials
        ),
        "SCALAR_FIELD_GENERATORS": scalar_generators,
        "F_TEMPLATES": f_templates,
        "R_TEMPLATES": r_templates,
        "DIM_F": sp.Integer(dim_f),
        "DIM_F_PLUS_R": sp.Integer(dim_f_plus_r),
        "RANK_GAIN": sp.Integer(dim_f_plus_r - dim_f),
    }


def frozen_free_coefficients(
    selected: tuple[int, ...],
) -> tuple[sp.Symbol, ...]:
    """Recover the independent constant amplitudes carried by the density."""
    coefficients: list[sp.Symbol] = []
    excluded_structural_knobs = {
        M.W0,
        M.L_W,
        M.eta_bg,
        M.sigma_W,
        M.w1_profile,
        M.m1_profile,
        *M.w1_grad,
        *M.m1_grad,
    }
    for index in M.UNIFORM_SELECTED:
        label, invariant = M.UNIFORM_CANDIDATES[index]
        coefficient = M.uniform_coefficient(label, invariant)
        for atom in sorted(coefficient.free_symbols, key=sp.default_sort_key):
            if (
                M.DECLARED_SYMBOLS.get(atom, {}).get("class") == "KNOB"
                and atom not in excluded_structural_knobs
                and atom not in coefficients
            ):
                coefficients.append(atom)
    for source, _ in SOURCES:
        coefficients.extend(M.NEW_COEFFICIENTS[source][index] for index in selected)
    return tuple(coefficients)


def main() -> None:
    # Enumerate once on the abstract spurion, then retain the engine's quotient
    # split only as provenance for selected and completion coefficients.
    candidates = M.enumerate_new_candidates(M.bg)
    expressions = tuple(expression for _, expression in candidates)
    signatures = M.basis_euler_signatures(expressions, M.basis_fields)
    selected, omitted = M.quotient_independent_indices(expressions, signatures)

    energy = M.construct_energy("LAB_HELD")
    density_frozen = energy.density
    density_uniform = sp.expand(
        sp.Add(
            *(
                term
                for label, term in energy.terms
                if not label.startswith(("W_BG_", "MU_R_BG_"))
            )
        )
    )
    density_selected_spurion = sp.expand(
        sp.Add(
            *(
                term
                for label, term in energy.terms
                if label.startswith(("W_BG_", "MU_R_BG_"))
            )
        )
    )

    omitted_coefficients: dict[tuple[str, int], sp.Symbol] = {}
    omitted_terms: list[sp.Expr] = []
    for source, actual_vector in SOURCES:
        substitution = M.live_basis_substitution(actual_vector)
        for index in omitted:
            coefficient = sp.Symbol(
                f"c_s11cb_88_{source.lower()}_{index + 1:02d}"
            )
            template = sp.expand(
                candidates[index][1].subs(substitution, simultaneous=True)
            )
            omitted_coefficients[(source, index)] = coefficient
            omitted_terms.append(sp.expand(coefficient * template))
    density_omitted = sp.expand(sp.Add(*omitted_terms))
    density_correct = sp.expand(density_frozen + density_omitted)

    # Extend the committed global derivative map only by the committed
    # operator_dx/background_dx second-background-jet assignments.
    w_hessian: dict[tuple[int, int], sp.Symbol] = {}
    m_hessian: dict[tuple[int, int], sp.Symbol] = {}
    for i in M.DIRECTIONS:
        for j in range(i, 3):
            index = M.sorted_index(i, j)
            w_hessian[index] = M.INCOMING_LEDGER[
                f"w1_profile_d{i + 1}d{j + 1}"
            ]["value"]
            m_hessian[index] = M.symbol(
                f"m1_profile_d{i + 1}d{j + 1}",
                "KNOB",
                f"dimensionless modulus-profile second jet {i + 1},{j + 1}",
                M.DIM_ZERO,
            )

    def frozen_dx(expression: sp.Expr, direction: int) -> sp.Expr:
        return M.dx(expression, direction)

    def hessian_dx(expression: sp.Expr, direction: int) -> sp.Expr:
        derivative_map = dict(M.DERIVATIVE_MAP[direction])
        for jet_direction in M.DIRECTIONS:
            jet_index = M.sorted_index(direction, jet_direction)
            derivative_map[M.grad_W[jet_direction]] = (
                M.sigma_W * w_hessian[jet_index] / M.L_W
            )
            derivative_map[M.grad_mu[jet_direction]] = (
                M.mu_R * M.sigma_W * m_hessian[jet_index] / (M.W0 * M.L_W)
            )
        result = sp.Integer(0)
        for atom, derivative in derivative_map.items():
            if atom in expression.free_symbols:
                result += sp.diff(expression, atom) * derivative
        return sp.expand(result)

    # CONTROL_SENTINEL: print raw and retained-grade chain-rule operands.
    sentinel_rows = []
    for source, expression in (
        ("W_BG", M.grad_W[0] * M.u[0]),
        ("MU_R_BG", M.grad_mu[0] * M.u[0]),
    ):
        frozen_raw = frozen_dx(expression, 0)
        hessian_raw = hessian_dx(expression, 0)
        sentinel_rows.append(
            {
                "SOURCE": source,
                "OPERAND": expression,
                "FROZEN_DX": frozen_raw,
                "HESSIAN_DX": hessian_raw,
                "FROZEN_DX_RETAINED": M.first_shape_series(frozen_raw),
                "HESSIAN_DX_RETAINED": M.first_shape_series(hessian_raw),
            }
        )
    emit("S11CB_88_CONTROL_SENTINEL", tuple(sentinel_rows))

    # CONTROL_JACOBIAN: the completion templates come from live-substituting
    # only after the abstract omitted indices have been selected.
    jacobian_rows = []
    jacobian_templates = tuple(
        sp.expand(sp.diff(density_omitted, omitted_coefficients[key]))
        for key in omitted_coefficients
    )
    jacobian_termcount = sum(termcount(item) for item in jacobian_templates)
    distinct_templates = len(set(jacobian_templates)) == len(jacobian_templates)
    for key, template in zip(omitted_coefficients, jacobian_templates):
        source, index = key
        correct_vector = M.grad_W if source == "W_BG" else M.grad_mu
        wrong_vector = M.grad_mu if source == "W_BG" else M.grad_W
        jacobian_rows.append(
            {
                "SOURCE": source,
                "INDEX": sp.Integer(index),
                "LABEL": candidates[index][0],
                "COEFFICIENT": omitted_coefficients[key],
                "TEMPLATE": template,
                "NONZERO": template != 0,
                "CARRIES_SOURCE_JET": bool(set(correct_vector) & template.free_symbols),
                "CARRIES_OTHER_SOURCE_JET": bool(set(wrong_vector) & template.free_symbols),
            }
        )
    emit(
        "S11CB_88_CONTROL_JACOBIAN",
        {
            "SELECTED_ZERO_BASED": selected,
            "OMITTED_ZERO_BASED": omitted,
            "JACOBIAN_TERMCOUNT": sp.Integer(jacobian_termcount),
            "DISTINCT_TEMPLATES": distinct_templates,
            "ROWS": tuple(jacobian_rows),
        },
    )
    assert jacobian_termcount > 0

    row_specs = (
        ("U_MOMENTUM_1", M.u[0], M.grad_u[0]),
        ("U_MOMENTUM_2", M.u[1], M.grad_u[1]),
        ("U_MOMENTUM_3", M.u[2], M.grad_u[2]),
        ("MU_THETA", M.theta, M.grad_theta),
        ("THICKNESS_EW", M.e_W, M.grad_e),
    )

    engine_operator, _ = M.operator_from_density(
        density_frozen, "RHO4_CONSTANT", include_kinetic=False
    )
    engine_u = named_form(engine_operator["U_BODY_BALANCE"], "EXPANDED")
    engine_theta = named_form(engine_operator["THETA_BALANCE"], "EXPANDED")
    engine_e = named_form(engine_operator["E_W_BALANCE"], "EXPANDED")
    engine_rows = {
        "U_MOMENTUM_1": engine_u[0],
        "U_MOMENTUM_2": engine_u[1],
        "U_MOMENTUM_3": engine_u[2],
        "MU_THETA": engine_theta,
        "THICKNESS_EW": engine_e,
    }

    frozen_coefficients = frozen_free_coefficients(selected)
    all_free_coefficients = set(frozen_coefficients) | set(omitted_coefficients.values())
    second_jet_atoms = set(w_hessian.values()) | set(m_hessian.values())
    engine_controls = []
    recon_controls = []

    for row_name, field, first_jets in row_specs:
        frozen_raw = el_amplitude(density_frozen, field, first_jets, frozen_dx)
        correct_raw = el_amplitude(density_correct, field, first_jets, hessian_dx)
        pre_truncation_residual = sp.expand(correct_raw - frozen_raw)
        el_frozen = M.first_shape_series(frozen_raw)
        el_correct = M.first_shape_series(correct_raw)
        residual = sp.expand(el_correct - el_frozen)

        uniform_hessian = M.first_shape_series(
            sp.expand(
                el_amplitude(density_uniform, field, first_jets, hessian_dx)
                - el_amplitude(density_uniform, field, first_jets, frozen_dx)
            )
        )
        selected_spurion_hessian = M.first_shape_series(
            sp.expand(
                el_amplitude(
                    density_selected_spurion, field, first_jets, hessian_dx
                )
                - el_amplitude(
                    density_selected_spurion, field, first_jets, frozen_dx
                )
            )
        )
        omitted_correct = M.first_shape_series(
            el_amplitude(density_omitted, field, first_jets, hessian_dx)
        )
        recon_difference = sp.expand(
            residual
            - uniform_hessian
            - selected_spurion_hessian
            - omitted_correct
        )
        new_second_jets = tuple(
            sorted(
                (second_jet_atoms & residual.free_symbols) - el_frozen.free_symbols,
                key=sp.default_sort_key,
            )
        )

        emit(
            f"S11CB_88_ROW_{row_name}",
            {
                "EL_FROZEN": el_frozen,
                "EL_CORRECT": el_correct,
                "RESIDUAL": residual,
                "RESIDUAL_TERMCOUNT": sp.Integer(termcount(residual)),
                "UNIFORM_HESSIAN": uniform_hessian,
                "UNIFORM_HESSIAN_TERMCOUNT": sp.Integer(termcount(uniform_hessian)),
                "SELECTED_SPURION_HESSIAN": selected_spurion_hessian,
                "SELECTED_SPURION_HESSIAN_TERMCOUNT": sp.Integer(
                    termcount(selected_spurion_hessian)
                ),
                "OMITTED_CORRECT": omitted_correct,
                "OMITTED_CORRECT_TERMCOUNT": sp.Integer(termcount(omitted_correct)),
                "NEW_SECOND_JET_ATOMS": new_second_jets,
                "GRADE_SUPPORT": grade_support(residual),
                "PRE_TRUNCATION_TERMCOUNT": sp.Integer(
                    termcount(pre_truncation_residual)
                ),
            },
        )

        recon_controls.append(
            {"ROW": row_name, "RECON_DIFFERENCE": recon_difference}
        )

        engine_amplitude = M.first_shape_series(
            sp.expand(engine_rows[row_name] / M.epsilon)
        )
        engine_difference = sp.expand(engine_amplitude - el_frozen)
        engine_controls.append(
            {
                "ROW": row_name,
                "ENGINE_AMPLITUDE": engine_amplitude,
                "INSTRUMENT_EL_FROZEN": el_frozen,
                "DIFFERENCE": engine_difference,
            }
        )

        f_templates = tuple(
            sp.expand(sp.diff(el_frozen, coefficient))
            for coefficient in frozen_coefficients
        )
        r_templates = tuple(
            sp.expand(sp.diff(residual, coefficient))
            for coefficient in omitted_coefficients.values()
        )
        span_data = constant_span_data(
            f_templates, r_templates, all_free_coefficients
        )
        emit(
            f"S11CB_88_ABSORB_{row_name}",
            {
                "F_COEFFICIENTS": frozen_coefficients,
                "R_COEFFICIENTS": tuple(omitted_coefficients.values()),
                **span_data,
            },
        )

    emit("S11CB_88_CONTROL_ENGINE", tuple(engine_controls))
    for control in engine_controls:
        assert control["DIFFERENCE"] == 0

    emit("S11CB_88_CONTROL_RECON", tuple(recon_controls))
    for control in recon_controls:
        assert control["RECON_DIFFERENCE"] == 0

    emit(
        "S11CB_88_SCOPE",
        {
            "BRANCH": "LAB_HELD",
            "OBJECT": "STORED_ENERGY_EL_ROWS",
            "ENGINE_SIDE": "PY",
            "MEASUREMENT": "PY_SIDE_DISTURBANCE_WITNESS",
            "NONZERO_ROW_SCOPE": "EMITTED_PY_ROW_VS_SECTION_1D_OPERATOR",
            "ZERO_ROW_SCOPE": (
                "NO_WL_FROZEN_OR_MATERIAL_ADVECTED_CLEARANCE;_RESERVED_FOR_89"
            ),
            "NOT_MEASURED": ("WL_FROZEN_ROWS", "MATERIAL_ADVECTED_ROWS"),
        },
    )


if __name__ == "__main__":
    main()
