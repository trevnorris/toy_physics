#!/usr/bin/env python3
"""pathA_36 C5 phase-potential test, SymPy engine.

This script starts from the primitive first-order theta/Josephson Lagrangian
required by the directive.  The Maxwell square is only reconstructed as a
diagnostic after the canonical constraints have been derived.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"
SYM_OUT = SCRATCH / "pathA_36_c5_sympy.json"
WL_OUT = SCRATCH / "pathA_36_c5_mathematica.json"
YAML_OUT = REPORTS / "pathA_36_c5_phase_potential_results.yaml"

SCHEMA = "pathA_36_c5_phase_potential/v1"


Dim = tuple[int, int, int]


def dadd(*dims: Dim) -> Dim:
    return tuple(sum(dim[i] for dim in dims) for i in range(3))  # type: ignore[return-value]


def dmul(n: int, dim: Dim) -> Dim:
    return (n * dim[0], n * dim[1], n * dim[2])


def dsub(left: Dim, right: Dim) -> Dim:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def dim_str(dim: Dim) -> str:
    labels = ("M", "L", "T")
    parts: list[str] = []
    for label, power in zip(labels, dim):
        if power == 0:
            continue
        parts.append(label if power == 1 else f"{label}^{power}")
    return " ".join(parts) if parts else "1"


def dim_record(dim: Dim) -> dict[str, Any]:
    return {"triple_MLT": list(dim), "string": dim_str(dim)}


class DimChecker:
    def __init__(self) -> None:
        self.records: list[dict[str, Any]] = []
        self.ablations: list[dict[str, Any]] = []

    def check(self, category: str, name: str, actual: Dim, expected: Dim, expression: str) -> None:
        if actual != expected:
            raise AssertionError(
                f"{category}:{name}: expected {expected} ({dim_str(expected)}), "
                f"got {actual} ({dim_str(actual)})"
            )
        self.records.append(
            {
                "category": category,
                "name": name,
                "expression": expression,
                "dimension": dim_record(actual),
                "expected": dim_record(expected),
                "status": "PASS",
            }
        )

    def expect_fail(self, category: str, name: str, actual: Dim, expected: Dim, expression: str) -> None:
        if actual == expected:
            raise AssertionError(f"dimension ablation did not fire: {category}:{name}")
        self.ablations.append(
            {
                "category": category,
                "name": name,
                "expression": expression,
                "actual": dim_record(actual),
                "expected": dim_record(expected),
                "status": "FIRED",
            }
        )


def s(expr: Any) -> str:
    return sp.sstr(sp.factor(sp.simplify(expr)))


def is_zero(expr: sp.Expr) -> bool:
    return sp.simplify(expr) == 0


rho_br, mu_R, B, J, rho_B0, K_theta, chi_c = sp.symbols(
    "rho_br mu_R B J rho_B0 K_theta chi_c", positive=True
)
k, omega = sp.symbols("k omega", positive=True)
kappa_phase, beta_B, m_theta2 = sp.symbols("kappa_phase beta_B m_theta2", positive=True)

q, theta, vq, vt, p_q, pi_theta, delta_rho = sp.symbols(
    "u_L theta dot_u_L dot_theta p_u_L pi_theta delta_rho_B"
)

C_J = -J * rho_B0
C2 = sp.simplify(C_J**2)


def build_dimensions() -> dict[str, Any]:
    check = DimChecker()
    brane_lag: Dim = (1, -1, -2)
    Z: Dim = (0, 0, 0)

    d_u: Dim = (0, 1, 0)
    d_theta: Dim = Z
    d_grad: Dim = (0, -1, 0)
    d_dt: Dim = (0, 0, -1)
    d_k: Dim = d_grad
    d_omega: Dim = d_dt
    d_rho_br: Dim = (1, -3, 0)
    d_mu_R: Dim = brane_lag
    d_B: Dim = brane_lag
    d_rho_B0: Dim = d_rho_br
    d_J: Dim = (0, 2, -1)
    d_CJ: Dim = (1, -1, -1)
    d_K_theta: Dim = (1, 1, -2)
    d_chi_c: Dim = (1, -5, 2)
    d_mtheta2: Dim = brane_lag

    d_div_u = dadd(d_grad, d_u)
    d_curl_u = d_div_u
    d_grad_theta = dadd(d_grad, d_theta)
    d_dt_u = dadd(d_dt, d_u)
    d_dt_theta = dadd(d_dt, d_theta)
    d_delta_rho = d_rho_B0

    check.check("primitive", "brane_inertia", dadd(d_rho_br, dmul(2, d_dt_u)), brane_lag, "rho_br (partial_t u)^2")
    check.check("primitive", "MacCullagh_curl", dadd(d_mu_R, dmul(2, d_curl_u)), brane_lag, "mu_R (nabla x u)^2")
    check.check("primitive", "Cauchy_bulk", dadd(d_B, dmul(2, d_div_u)), brane_lag, "B (nabla dot u)^2")
    check.check("primitive", "Josephson_density", dadd(d_J, d_dt_theta, d_delta_rho), brane_lag, "J (partial_t theta) delta_rho_B")
    check.check("primitive", "slaved_delta_rho", dadd(d_rho_B0, d_div_u), d_delta_rho, "rho_B0 nabla dot u")
    check.check("primitive", "phase_gradient_signed", dadd(d_K_theta, dmul(2, d_grad_theta)), brane_lag, "K_theta (nabla theta)^2")
    check.check("primitive", "density_compressibility", dsub(dmul(2, d_delta_rho), d_chi_c), brane_lag, "delta_rho_B^2/chi_c")
    check.check("primitive", "theta_mass_locking", dadd(d_mtheta2, dmul(2, d_theta)), brane_lag, "m_theta^2 theta^2")

    check.check("ibp", "C_J_definition", dadd(d_J, d_rho_B0), d_CJ, "C_J = -J rho_B0")
    check.check("ibp", "Josephson_cross", dadd(d_CJ, d_dt_u, d_grad_theta), brane_lag, "C_J partial_t u dot nabla theta")
    check.check("ibp", "electric_square_velocity_piece", dadd(d_rho_br, dmul(2, d_dt_u)), brane_lag, "rho_br (partial_t u)^2")
    check.check("ibp", "electric_square_mixed_piece", dadd(d_CJ, d_dt_u, d_grad_theta), brane_lag, "C_J partial_t u dot nabla theta")
    check.check("ibp", "electric_square_gradient_piece", dadd(d_K_theta, dmul(2, d_grad_theta)), brane_lag, "K_theta (nabla theta)^2")
    check.check("derived", "Maxwell_locus_CJ2", dmul(2, d_CJ), dadd(d_rho_br, d_K_theta), "C_J^2 = rho_br K_theta")
    check.check("derived", "c_gamma_squared", dsub(d_mu_R, d_rho_br), (0, 2, -2), "c_gamma^2 = mu_R/rho_br")
    check.check("branch_a", "induced_theta_kinetic_without_continuity", dadd(d_chi_c, dmul(2, d_J), dmul(2, d_dt_theta)), brane_lag, "chi_c J^2 (partial_t theta)^2")
    check.check("branch_b", "slaved_compressibility_B_eff", dsub(dmul(2, d_rho_B0), d_chi_c), d_B, "rho_B0^2/chi_c")
    check.check("fourier", "omega2_transverse", dadd(dsub(d_mu_R, d_rho_br), dmul(2, d_k)), dmul(2, d_omega), "(mu_R/rho_br) k^2")

    check.expect_fail("ablation", "drop_rho_B0_from_Josephson_cross", dadd(d_J, d_dt_u, d_grad_theta), brane_lag, "J partial_t u dot nabla theta")
    check.expect_fail("ablation", "phase_stiffness_without_gradient", dadd(d_K_theta, dmul(2, d_theta)), brane_lag, "K_theta theta^2")
    check.expect_fail("ablation", "compressibility_multiplied_not_divided", dadd(d_chi_c, dmul(2, d_delta_rho)), brane_lag, "chi_c delta_rho_B^2")

    return {
        "pass": True,
        "checked_expression_count": len(check.records),
        "records": check.records,
        "ablations": check.ablations,
        "constants": {
            "rho_br": dim_record(d_rho_br),
            "mu_R": dim_record(d_mu_R),
            "B": dim_record(d_B),
            "J": dim_record(d_J),
            "rho_B0": dim_record(d_rho_B0),
            "C_J": dim_record(d_CJ),
            "K_theta": dim_record(d_K_theta),
            "chi_c": dim_record(d_chi_c),
            "m_theta2_control": dim_record(d_mtheta2),
        },
    }


def poisson(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.simplify(
        sp.diff(f, q) * sp.diff(g, p_q)
        - sp.diff(f, p_q) * sp.diff(g, q)
        + sp.diff(f, theta) * sp.diff(g, pi_theta)
        - sp.diff(f, pi_theta) * sp.diff(g, theta)
    )


def sign_label(expr: sp.Expr) -> str:
    expr = sp.simplify(expr)
    if expr == 0:
        return "zero"
    if expr.is_positive:
        return "positive"
    if expr.is_negative:
        return "negative"
    return "symbolic"


def first_order_analysis(name: str, K_expr: sp.Expr, B_eff_expr: sp.Expr, *, theta_mass_expr: sp.Expr = sp.Integer(0)) -> dict[str, Any]:
    """Analyze L = 1/2 rho vq^2 - C k q vt + 1/2 S theta^2 - 1/2 B k^2 q^2."""

    S = sp.simplify(K_expr * k**2 - theta_mass_expr)
    primitive_L = sp.Rational(1, 2) * rho_br * vq**2 - C_J * k * q * vt + sp.Rational(1, 2) * S * theta**2 - sp.Rational(1, 2) * B_eff_expr * k**2 * q**2
    p_u = sp.diff(primitive_L, vq)
    p_th = sp.diff(primitive_L, vt)
    primary = sp.simplify(pi_theta - p_th)
    H = sp.simplify(p_q**2 / (2 * rho_br) + sp.Rational(1, 2) * B_eff_expr * k**2 * q**2 - sp.Rational(1, 2) * S * theta**2)
    secondary = sp.simplify(poisson(primary, H))
    bracket = sp.simplify(poisson(primary, secondary))
    secondary_preservation_no_multiplier = sp.simplify(poisson(secondary, H))
    determinant = sp.factor((rho_br * S - C2 * k**2) * omega**2 - B_eff_expr * k**2 * S)
    A_eff = sp.simplify(rho_br - C2 * k**2 / S) if not is_zero(S) else sp.oo
    if is_zero(B_eff_expr):
        omega2_text = "0"
    elif A_eff != sp.oo and not is_zero(A_eff):
        omega2_text = s(sp.simplify(B_eff_expr * k**2 / A_eff))
    else:
        omega2_text = "singular_second_class_no_propagating_pole"

    bracket_zero = is_zero(bracket)
    b_zero = is_zero(B_eff_expr)
    mass_zero = is_zero(theta_mass_expr)
    all_k_gauge = bracket_zero and b_zero and mass_zero

    if all_k_gauge:
        first_class = 2
        second_class = 0
        physical_dof = 0
        initial_data_functions = 0
        constraint_structure = "FIRST_CLASS_MAXWELL_CHAIN"
        gauge_closure = "OFF_SHELL_ALL_K_LOCAL"
        verdict = "C5_RESOLVED_MAXWELL_BY_TUNING"
        pole_count = 0
        residue_sign = "no_physical_longitudinal_pole"
        bounded = True
        hamiltonian_status = "bounded_on_Gauss_constraint"
    elif bracket_zero and not b_zero:
        first_class = 0
        second_class = 4
        physical_dof = 0
        initial_data_functions = 0
        constraint_structure = "TERTIARY_SECOND_CLASS_CHAIN"
        gauge_closure = "BROKEN_BY_CAUCHY_TERM"
        verdict = "FAIL_SECOND_CLASS_NOT_MAXWELL"
        pole_count = 0
        residue_sign = "no_gauge_pole_but_not_first_class"
        bounded = True
        hamiltonian_status = "bounded_but_second_class"
    else:
        first_class = 0
        second_class = 2
        physical_dof = 1
        initial_data_functions = 2
        constraint_structure = "SECOND_CLASS_PAIR"
        gauge_closure = "BROKEN_OFF_LOCUS"
        pole_count = 1
        residue_sign = sign_label(1 / A_eff) if A_eff not in (sp.oo, -sp.oo) else "undefined"
        bounded = A_eff.is_positive is True and (b_zero or B_eff_expr.is_positive is True)
        hamiltonian_status = "bounded_reduced_Hamiltonian" if bounded else "unbounded_or_negative_residue"
        if theta_mass_expr != 0:
            verdict = "FAIL_SECOND_CLASS_NOT_MAXWELL"
        elif bounded is False:
            verdict = "FAIL_GHOST_OR_NEGATIVE_NORM"
        elif b_zero:
            verdict = "FAIL_C5_LONGITUDINAL_ZERO_MODE"
        else:
            verdict = "FAIL_CAUCHY_STRAY_LONGITUDINAL"

    square_residual = sp.simplify(K_expr - C2 / rho_br)
    local_gauge_generator = None
    if all_k_gauge:
        local_gauge_generator = {
            "generator": "G[chi]=(rho_br/C_J)*(chi*Phi_2-dot(chi)*Phi_1)",
            "delta_u_L": "k chi",
            "delta_theta": "-(rho_br/C_J) dot(chi)",
            "local": True,
            "inverse_k_or_inverse_omega": False,
        }

    return {
        "name": name,
        "primitive_lagrangian": "1/2 rho_br dot_u_L^2 - C_J k u_L dot_theta + 1/2 (K_theta k^2 - m_theta^2) theta^2 - 1/2 B_eff k^2 u_L^2",
        "coefficients": {
            "C_J": "-J*rho_B0",
            "K_theta": s(K_expr),
            "B_eff": s(B_eff_expr),
            "theta_mass2": s(theta_mass_expr),
            "S_theta": s(S),
        },
        "momenta": {
            "p_u_L": s(p_u),
            "pi_theta": s(p_th),
        },
        "constraints": {
            "primary": s(primary),
            "secondary": s(secondary),
            "constraint_bracket": s(bracket),
            "secondary_preservation_no_multiplier": s(secondary_preservation_no_multiplier),
            "first_class_count": first_class,
            "second_class_count": second_class,
            "classification": constraint_structure,
        },
        "gauge": {
            "closure": gauge_closure,
            "square_residual_K_minus_CJ2_over_rho": s(square_residual),
            "B_eff_required_zero": b_zero,
            "theta_mass_required_zero": mass_zero,
            "generator": local_gauge_generator,
        },
        "hamiltonian": {
            "canonical_H": "p_u_L^2/(2 rho_br) + 1/2 B_eff k^2 u_L^2 - 1/2 (K_theta k^2 - m_theta^2) theta^2",
            "reduced_kinetic_A": s(A_eff),
            "bounded": bounded,
            "status": hamiltonian_status,
        },
        "dispersion": {
            "determinant": s(determinant),
            "omega2": omega2_text,
            "pole_count": pole_count,
            "pole_residue_sign": residue_sign,
        },
        "mode_count": {
            "physical_dof_per_finite_k": physical_dof,
            "independent_initial_data_functions_per_k": initial_data_functions,
        },
        "verdict": verdict,
    }


def elastic_longitudinal_control(name: str, B_expr: sp.Expr) -> dict[str, Any]:
    omega2_expr = sp.simplify(B_expr * k**2 / rho_br)
    propagating = not is_zero(B_expr)
    return {
        "name": name,
        "primitive_lagrangian": "1/2 rho_br dot_u_L^2 - 1/2 B_eff k^2 u_L^2",
        "constraints": {
            "primary": [],
            "secondary": [],
            "first_class_count": 0,
            "second_class_count": 0,
            "classification": "UNCONSTRAINED_ELASTIC_COORDINATE",
        },
        "dispersion": {
            "omega2": s(omega2_expr),
            "pole_count": 1,
            "pole_residue_sign": "positive",
        },
        "mode_count": {
            "physical_dof_per_finite_k": 1,
            "independent_initial_data_functions_per_k": 2,
        },
        "verdict": "FAIL_CAUCHY_STRAY_LONGITUDINAL" if propagating else "FAIL_C5_LONGITUDINAL_ZERO_MODE",
    }


def decoupled_second_order_theta_control() -> dict[str, Any]:
    M_theta = sp.simplify(chi_c * J**2)
    theta_omega2 = sp.simplify(kappa_phase * k**2 / M_theta)
    return {
        "name": "decoupled_theta_independent_density_no_continuity",
        "primitive_lagrangian": "1/2 rho_br dot_u_L^2 + 1/2 chi_c J^2 dot_theta^2 - 1/2 kappa_phase k^2 theta^2",
        "hessian_rank": 2,
        "constraints": {
            "first_class_count": 0,
            "second_class_count": 0,
            "classification": "TWO_UNCONSTRAINED_SECOND_ORDER_FIELDS",
        },
        "dispersion": {
            "u_L": "omega^2 = 0",
            "theta": f"omega^2 = {s(theta_omega2)}",
            "pole_count": 2,
            "pole_residue_sign": "positive",
        },
        "mode_count": {
            "physical_dof_per_finite_k": 2,
            "independent_initial_data_functions_per_k": 4,
        },
        "verdict": "FAIL_EXTRA_SCALAR_DOF",
    }


def decoupled_slaved_theta_control() -> dict[str, Any]:
    """C_J=0 with no density/amplitude theta kinetic: theta is algebraic."""

    primitive_L = sp.Rational(1, 2) * rho_br * vq**2 + sp.Rational(1, 2) * (-kappa_phase) * k**2 * theta**2
    p_u = sp.diff(primitive_L, vq)
    p_th = sp.diff(primitive_L, vt)
    primary = sp.simplify(pi_theta - p_th)
    H = sp.simplify(p_q**2 / (2 * rho_br) + sp.Rational(1, 2) * kappa_phase * k**2 * theta**2)
    secondary = sp.simplify(poisson(primary, H))
    bracket = sp.simplify(poisson(primary, secondary))
    return {
        "name": "decoupled_theta_slaved_CJ_zero",
        "primitive_lagrangian": "1/2 rho_br dot_u_L^2 - 1/2 kappa_phase k^2 theta^2",
        "momenta": {"p_u_L": s(p_u), "pi_theta": s(p_th)},
        "constraints": {
            "primary": s(primary),
            "secondary": s(secondary),
            "constraint_bracket": s(bracket),
            "first_class_count": 0,
            "second_class_count": 2,
            "classification": "THETA_ALGEBRAIC_SECOND_CLASS_PAIR_PLUS_U_ZERO_MODE",
        },
        "dispersion": {
            "u_L": "omega^2 = 0",
            "theta": "algebraic",
            "pole_count": 1,
            "pole_residue_sign": "positive",
        },
        "mode_count": {
            "physical_dof_per_finite_k": 1,
            "independent_initial_data_functions_per_k": 2,
        },
        "decoupled_theta_status": "NON_DYNAMICAL_ALGEBRAIC_PHASE",
        "verdict": "FAIL_C5_LONGITUDINAL_ZERO_MODE",
    }


def epsilon_mismatch_control() -> dict[str, Any]:
    """A closed longitudinal square with epsilon != rho_br must not pass."""

    rho_eps = 2 * rho_br
    K_eps = sp.simplify(C2 / rho_eps)
    S = sp.simplify(K_eps * k**2)
    primitive_L = sp.Rational(1, 2) * rho_eps * vq**2 - C_J * k * q * vt + sp.Rational(1, 2) * S * theta**2
    p_u = sp.diff(primitive_L, vq)
    p_th = sp.diff(primitive_L, vt)
    primary = sp.simplify(pi_theta - p_th)
    H = sp.simplify(p_q**2 / (2 * rho_eps) - sp.Rational(1, 2) * S * theta**2)
    secondary = sp.simplify(poisson(primary, H))
    bracket = sp.simplify(poisson(primary, secondary))
    shifted_transverse = sp.simplify(mu_R * k**2 / rho_eps)
    frozen_transverse = sp.simplify(mu_R * k**2 / rho_br)
    return {
        "name": "epsilon_mismatch_square_closes_wrong_inertia",
        "primitive_lagrangian": "1/2 (2 rho_br) dot_u_L^2 - C_J k u_L dot_theta + 1/2 (C_J^2/(2 rho_br)) k^2 theta^2",
        "constraints": {
            "primary": s(primary),
            "secondary": s(secondary),
            "constraint_bracket": s(bracket),
            "first_class_count": 2 if is_zero(bracket) else 0,
            "second_class_count": 0 if is_zero(bracket) else 2,
            "classification": "FIRST_CLASS_LONGITUDINAL_BUT_WRONG_EPSILON" if is_zero(bracket) else "NOT_CLOSED",
        },
        "longitudinal_square": {
            "epsilon": "2*rho_br",
            "closes": is_zero(bracket),
            "physical_dof_per_finite_k": 0 if is_zero(bracket) else 1,
        },
        "transverse_check": {
            "frozen_omega2": s(frozen_transverse),
            "mismatched_omega2": s(shifted_transverse),
            "speed_shift": s(sp.simplify(shifted_transverse - frozen_transverse)),
        },
        "mode_count": {
            "physical_dof_per_finite_k": 0 if is_zero(bracket) else 1,
            "independent_initial_data_functions_per_k": 0 if is_zero(bracket) else 2,
        },
        "verdict": "FAIL_TRANSVERSE_DISTURBED",
    }


def branch_a_integration_proof() -> dict[str, Any]:
    r = delta_rho
    L_density = sp.simplify(J * vt * r - r**2 / (2 * chi_c))
    r_gaussian = sp.solve(sp.Eq(sp.diff(L_density, r), 0), r)[0]
    L_no_continuity = sp.simplify(L_density.subs(r, r_gaussian))
    r_continuity = -rho_B0 * k * q
    L_with_continuity = sp.simplify(L_density.subs(r, r_continuity))
    B_increment = sp.simplify(rho_B0**2 / chi_c)

    return {
        "independent_field": "delta_rho_B",
        "continuity_constraint_fourier": "omega*(delta_rho_B + rho_B0*k*u_L)=0",
        "finite_frequency_solution_fixed_number_sector": "delta_rho_B = -rho_B0*k*u_L",
        "gaussian_solution_if_continuity_removed": f"delta_rho_B = {s(r_gaussian)}",
        "theta_kinetic_if_continuity_removed": s(L_no_continuity),
        "effective_density_L_with_continuity": s(L_with_continuity),
        "effective_C_J": "-J*rho_B0",
        "effective_B_increment": s(B_increment),
        "proof_status": "CONTINUITY_FORCES_SAME_SLAVED_SECTOR",
    }


def transverse_sector() -> dict[str, Any]:
    return {
        "basis": ["u_T1", "u_T2"],
        "primitive_lagrangian_per_polarization": "1/2 rho_br dot_u_T^2 - 1/2 mu_R k^2 u_T^2",
        "hessian_rank_per_polarization": 1,
        "physical_dof": 2,
        "omega2": "mu_R*k^2/rho_br",
        "c_gamma_squared": "mu_R/rho_br",
        "massless": True,
        "theta_couplings": 0,
        "verdict": "PASS_TRANSVERSE_UNDISTURBED",
    }


def locus_evaluation(maxwell_case: dict[str, Any]) -> dict[str, Any]:
    return {
        "conditions": {
            "i_square_closes": "C_J^2 = rho_br*K_theta",
            "i_requires_K_theta_positive": True,
            "ii_epsilon_equals_frozen_rho_br": True,
            "iii_B_eff_zero": "B + rho_B0^2/chi_c (if finite stiffness included) = 0",
            "iv_bounded_Hamiltonian": maxwell_case["hamiltonian"]["bounded"],
            "derived_C_J": "C_J = -J*rho_B0",
        },
        "free_locus_result": {
            "verdict": maxwell_case["verdict"],
            "physical_dof_per_finite_k": maxwell_case["mode_count"]["physical_dof_per_finite_k"],
            "constraint_classification": maxwell_case["constraints"]["classification"],
            "gauge_closure": maxwell_case["gauge"]["closure"],
        },
        "provenance_status": "BY_TUNING",
        "with_provenance": False,
        "not_forced_by_frozen_definitions": [
            "no frozen definition forces K_theta = J^2*rho_B0^2/rho_br",
            "the order-parameter phase has conventional signed gradient K_theta=-kappa_phase<0, while the Maxwell square requires K_theta>0",
            "finite conjugate-density stiffness contributes rho_B0^2/chi_c to B_eff and is not forced to vanish",
        ],
        "exact_tuning_conditions": [
            "K_theta -> J^2*rho_B0^2/rho_br",
            "B_eff -> 0, i.e. B + rho_B0^2/chi_c -> 0 if finite density stiffness is present",
            "m_theta2 -> 0",
        ],
    }


def control_status(name: str, analysis: dict[str, Any], expected: list[str]) -> dict[str, Any]:
    verdict = analysis["verdict"]
    if verdict not in expected:
        raise AssertionError(f"{name}: expected one of {expected}, got {verdict}")
    return {
        "name": name,
        "status": "FIRED",
        "verdict": verdict,
        "physical_dof_per_finite_k": analysis["mode_count"]["physical_dof_per_finite_k"],
        "initial_data_functions_per_k": analysis["mode_count"]["independent_initial_data_functions_per_k"],
    }


def control_agreement_label(item: dict[str, Any]) -> str:
    if "verdict" in item:
        return str(item["verdict"])
    if "free_locus_verdict" in item:
        return f"{item['fixed_coefficients_verdict']}->{item['free_locus_verdict']}:{item['provenance_status']}"
    if "absent_verdict" in item:
        return f"{item['absent_verdict']}->{item['included_verdict']}"
    raise KeyError(f"control lacks an agreement verdict: {item['name']}")


def build_payload() -> dict[str, Any]:
    dimensions = build_dimensions()

    K_conventional = -kappa_phase
    B_none = sp.Integer(0)
    B_finite_compressibility = sp.simplify(rho_B0**2 / chi_c)
    K_locus = sp.simplify(C2 / rho_br)

    branch_b_curl_only = first_order_analysis("branch_b_slaved_curl_only_conventional_K", K_conventional, B_none)
    branch_b_finite = first_order_analysis("branch_b_slaved_finite_compressibility_conventional_K", K_conventional, B_finite_compressibility)
    branch_b_tuned = first_order_analysis("branch_b_slaved_tuned_Maxwell_locus", K_locus, B_none)

    branch_a_proof = branch_a_integration_proof()
    branch_a_finite = first_order_analysis("branch_a_independent_with_continuity_integrated_out", K_conventional, B_finite_compressibility)
    branch_a_no_continuity = decoupled_second_order_theta_control()

    no_theta = elastic_longitudinal_control("no_theta_curl_only", sp.Integer(0))
    cauchy = elastic_longitudinal_control("cauchy_bulk_no_theta", beta_B)
    detuned_K = first_order_analysis("mismatched_positive_K_no_Cauchy", sp.simplify(2 * C2 / rho_br), sp.Integer(0))
    wrong_sign_K = first_order_analysis("mismatched_K_theta_le_0", K_conventional, sp.Integer(0))
    detuned_with_B = first_order_analysis("mismatched_positive_K_with_Cauchy", sp.simplify(2 * C2 / rho_br), beta_B)
    ghost_K = first_order_analysis("mismatched_positive_K_negative_residue", sp.simplify(C2 / (2 * rho_br)), beta_B)
    B_on_locus = first_order_analysis("B_nonzero_on_square_locus", K_locus, beta_B)
    epsilon_mismatch = epsilon_mismatch_control()
    decoupled_slaved = decoupled_slaved_theta_control()
    theta_mass = first_order_analysis("theta_mass_breaks_gauge", K_locus, sp.Integer(0), theta_mass_expr=m_theta2)
    transverse = transverse_sector()

    controls = [
        control_status("1_no_theta", no_theta, ["FAIL_C5_LONGITUDINAL_ZERO_MODE"]),
        control_status("2_cauchy_bulk", cauchy, ["FAIL_CAUCHY_STRAY_LONGITUDINAL"]),
        control_status("3_mismatched_positive_K_no_B", detuned_K, ["FAIL_C5_LONGITUDINAL_ZERO_MODE"]),
        control_status("3_mismatched_K_theta_le_0", wrong_sign_K, ["FAIL_C5_LONGITUDINAL_ZERO_MODE"]),
        control_status("3_mismatched_positive_K_with_B", detuned_with_B, ["FAIL_CAUCHY_STRAY_LONGITUDINAL"]),
        control_status("3_mismatched_positive_K_ghost", ghost_K, ["FAIL_GHOST_OR_NEGATIVE_NORM"]),
        control_status("3_B_nonzero_on_locus", B_on_locus, ["FAIL_SECOND_CLASS_NOT_MAXWELL"]),
        {
            "name": "3_epsilon_mismatch",
            "status": "FIRED",
            "verdict": epsilon_mismatch["verdict"],
            "longitudinal_dof_per_finite_k": epsilon_mismatch["mode_count"]["physical_dof_per_finite_k"],
            "frozen_transverse_omega2": epsilon_mismatch["transverse_check"]["frozen_omega2"],
            "mismatched_transverse_omega2": epsilon_mismatch["transverse_check"]["mismatched_omega2"],
        },
        control_status("4_decoupled_theta_slaved", decoupled_slaved, ["FAIL_C5_LONGITUDINAL_ZERO_MODE"]),
        control_status("4_decoupled_theta_independent", branch_a_no_continuity, ["FAIL_EXTRA_SCALAR_DOF"]),
        {
            "name": "5_transverse_undisturbed",
            "status": "FIRED",
            "verdict": transverse["verdict"],
            "physical_dof": transverse["physical_dof"],
            "omega2": transverse["omega2"],
        },
        {
            "name": "6_provenance_ablation",
            "status": "FIRED",
            "fixed_coefficients_verdict": branch_b_finite["verdict"],
            "free_locus_verdict": branch_b_tuned["verdict"],
            "provenance_status": "BY_TUNING",
        },
        {
            "name": "7_compressibility_absent_vs_included",
            "status": "FIRED",
            "absent_verdict": branch_b_curl_only["verdict"],
            "included_verdict": branch_b_finite["verdict"],
            "B_increment": "rho_B0^2/chi_c",
        },
        control_status("8_theta_mass", theta_mass, ["FAIL_SECOND_CLASS_NOT_MAXWELL", "FAIL_C5_LONGITUDINAL_ZERO_MODE"]),
    ]

    locus = locus_evaluation(branch_b_tuned)
    drift = {
        "branch_b_slaved": {
            "new_fields": ["theta"],
            "new_constants": ["B", "J", "rho_B0", "K_theta", "chi_c"],
            "count": 6,
            "tag": "DRIFT(6)",
        },
        "branch_a_independent": {
            "new_fields": ["theta", "delta_rho_B"],
            "new_constants": ["B", "J", "rho_B0", "K_theta", "chi_c"],
            "count": 7,
            "tag": "DRIFT(7)",
        },
    }

    headline = {
        "verdict": branch_b_finite["verdict"],
        "main_branch": "branch_b_slaved_finite_compressibility_conventional_K",
        "physical_reason": "finite density stiffness creates B_eff=rho_B0^2/chi_c and conventional K_theta=-kappa_phase<0 misses the electric Maxwell sign",
        "curl_only_conventional_subcase_verdict": branch_b_curl_only["verdict"],
        "tuned_locus_verdict": branch_b_tuned["verdict"],
        "provenance_status": locus["provenance_status"],
        "longitudinal_dof_main": branch_b_finite["mode_count"]["physical_dof_per_finite_k"],
        "longitudinal_initial_data_main": branch_b_finite["mode_count"]["independent_initial_data_functions_per_k"],
        "transverse_dof": transverse["physical_dof"],
        "engine_agreement": "PENDING",
    }

    agreement_payload = {
        "headline_verdict": headline["verdict"],
        "main_constraint_classification": branch_b_finite["constraints"]["classification"],
        "main_longitudinal_dof": branch_b_finite["mode_count"]["physical_dof_per_finite_k"],
        "main_initial_data_functions": branch_b_finite["mode_count"]["independent_initial_data_functions_per_k"],
        "main_pole_count": branch_b_finite["dispersion"]["pole_count"],
        "curl_only_verdict": branch_b_curl_only["verdict"],
        "curl_only_longitudinal_dof": branch_b_curl_only["mode_count"]["physical_dof_per_finite_k"],
        "tuned_locus_verdict": branch_b_tuned["verdict"],
        "tuned_constraint_classification": branch_b_tuned["constraints"]["classification"],
        "tuned_longitudinal_dof": branch_b_tuned["mode_count"]["physical_dof_per_finite_k"],
        "tuned_first_class_count": branch_b_tuned["constraints"]["first_class_count"],
        "locus_provenance": locus["provenance_status"],
        "branch_a_integrated_verdict": branch_a_finite["verdict"],
        "transverse_dof": transverse["physical_dof"],
        "transverse_omega2": transverse["omega2"],
        "control_verdicts": {item["name"]: control_agreement_label(item) for item in controls},
        "dimensional_ablations_fired": len(dimensions["ablations"]),
    }

    return {
        "schema": SCHEMA,
        "engine": "sympy",
        "directive": "software/stage1_solver/directives/pathA_36_c5_phase_potential.md",
        "frozen_baseline": "T0_SHEAR_FROZEN(d9520d3819c3)",
        "primitive_start": {
            "lagrangian": "1/2 rho_br dot(u)^2 - 1/2 mu_R (curl u)^2 - 1/2 B (div u)^2 + J dot(theta) delta_rho_B + 1/2 K_theta (grad theta)^2",
            "square_completion_used_as_input": False,
            "derived_C_J": "C_J = -J*rho_B0",
        },
        "headline": headline,
        "branches": {
            "branch_b_slaved_curl_only_conventional_K": branch_b_curl_only,
            "branch_b_slaved_finite_compressibility_conventional_K": branch_b_finite,
            "branch_b_slaved_tuned_Maxwell_locus": branch_b_tuned,
            "branch_a_integration_proof": branch_a_proof,
            "branch_a_independent_with_continuity_integrated_out": branch_a_finite,
            "branch_a_no_continuity_second_medium_ablation": branch_a_no_continuity,
            "epsilon_mismatch_control": epsilon_mismatch,
        },
        "locus": locus,
        "controls": controls,
        "transverse": transverse,
        "dimensional_firewall": dimensions,
        "drift": drift,
        "absence_tags": {
            "AXIS_RE_ADMITTED": False,
            "U_W_COLLISION": False,
        },
        "agreement_payload": agreement_payload,
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def compare_and_write_yaml(payload: dict[str, Any]) -> dict[str, Any]:
    if not WL_OUT.exists():
        raise FileNotFoundError(f"Mathematica payload missing: {WL_OUT}")
    other = json.loads(WL_OUT.read_text(encoding="utf-8"))
    if payload["agreement_payload"] != other["agreement_payload"]:
        diff = {
            "sympy": payload["agreement_payload"],
            "mathematica": other["agreement_payload"],
        }
        raise AssertionError("ENGINE_DISAGREE\n" + json.dumps(diff, indent=2, sort_keys=True))

    payload["headline"]["engine_agreement"] = "ENGINE_AGREE"
    results = {
        "schema": SCHEMA,
        "verdict": payload["headline"]["verdict"],
        "headline": payload["headline"],
        "engine_agreement": {
            "status": "ENGINE_AGREE",
            "compared_payload": payload["agreement_payload"],
            "sympy_payload": str(SYM_OUT),
            "mathematica_payload": str(WL_OUT),
        },
        "per_branch_sub_results": payload["branches"],
        "locus": payload["locus"],
        "controls": payload["controls"],
        "DRIFT": payload["drift"],
        "dimensional_firewall": {
            "pass": payload["dimensional_firewall"]["pass"],
            "checked_expression_count": payload["dimensional_firewall"]["checked_expression_count"],
            "ablations": payload["dimensional_firewall"]["ablations"],
        },
        "transverse": payload["transverse"],
    }
    YAML_OUT.write_text(yaml.safe_dump(results, sort_keys=False, width=120), encoding="utf-8")
    write_json(SYM_OUT, payload)
    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true", help="compare with Mathematica payload and write final YAML")
    args = parser.parse_args()

    payload = build_payload()
    write_json(SYM_OUT, payload)
    if args.compare:
        results = compare_and_write_yaml(payload)
        print(json.dumps({"engine": "sympy", "status": "ENGINE_AGREE", "verdict": results["verdict"]}, sort_keys=True))
    else:
        print(json.dumps({"engine": "sympy", "status": "OK", "verdict": payload["headline"]["verdict"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
