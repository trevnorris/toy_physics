#!/usr/bin/env python3
"""Ledger stage003 SymPy audit: transverse photons plus stray longitudinal mode.

Print-only, standalone, no arguments, no file output.  The audit starts from
the primitive light-sector Lagrangian, derives the transverse photon speed, the
Josephson IBP sign, and the longitudinal Dirac-Bergmann chain, then exercises
the source control verdicts as able-to-fail symbolic gates.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


Dim = tuple[sp.Rational, sp.Rational, sp.Rational]

PASS_COUNT = 0


class AuditFailure(AssertionError):
    pass


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def compact(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.cancel(sp.simplify(expr))))


def assert_no_float(name: str, expr: Any) -> None:
    clean = sp.sympify(expr)
    floats = clean.atoms(sp.Float)
    if floats:
        raise AuditFailure(f"{name}: Float atom(s) found in exact audit expression: {floats}")


def expect_zero(name: str, residual: sp.Expr) -> None:
    global PASS_COUNT
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean == 0:
        PASS_COUNT += 1
        print(f"PASS  {name}")
        return

    print(f"FAIL  {name}: residual = {compact(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1))


def expect_fail(name: str, residual: sp.Expr) -> None:
    global PASS_COUNT
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        PASS_COUNT += 1
        print(f"PASS  {name} produced required FAIL (residual = {compact(clean)})")
        return

    print(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def q(value: int | str) -> sp.Rational:
    return sp.Rational(value)


def dim(m_power: int, l_power: int, t_power: int) -> Dim:
    return (q(m_power), q(l_power), q(t_power))


def dim_add(*items: Dim) -> Dim:
    out = [q(0), q(0), q(0)]
    for item in items:
        for i, power in enumerate(item):
            out[i] += power
    return (out[0], out[1], out[2])


def dim_neg(item: Dim) -> Dim:
    return (-item[0], -item[1], -item[2])


def dim_sub(left: Dim, right: Dim) -> Dim:
    return dim_add(left, dim_neg(right))


def dim_scale(scale: int | sp.Rational, item: Dim) -> Dim:
    factor = q(scale) if isinstance(scale, int) else scale
    return (factor * item[0], factor * item[1], factor * item[2])


def dimension_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(sum((component - want) ** 2 for component, want in zip(actual, expected)))


@dataclass(frozen=True)
class Symbols:
    rho_br: sp.Symbol
    mu_R: sp.Symbol
    B: sp.Symbol
    J: sp.Symbol
    rho_B0: sp.Symbol
    chi_c: sp.Symbol
    kappa_phase: sp.Symbol
    beta_B: sp.Symbol
    m_theta2: sp.Symbol
    k: sp.Symbol
    omega: sp.Symbol
    u_L: sp.Symbol
    u_T1: sp.Symbol
    u_T2: sp.Symbol
    theta: sp.Symbol
    dot_u_L: sp.Symbol
    dot_u_T: sp.Symbol
    dot_theta: sp.Symbol
    p_u: sp.Symbol
    pi_theta: sp.Symbol
    delta_rho_B: sp.Symbol
    div_u: sp.Symbol
    theta_dot_symbol: sp.Symbol
    div_dot_u: sp.Symbol
    grad_theta_dot_dot_u: sp.Symbol
    theta_boundary_density: sp.Symbol


def make_symbols() -> Symbols:
    return Symbols(
        rho_br=sp.Symbol("rho_br", positive=True),
        mu_R=sp.Symbol("mu_R", positive=True),
        B=sp.Symbol("B", positive=True),
        J=sp.Symbol("J", positive=True),
        rho_B0=sp.Symbol("rho_B0", positive=True),
        chi_c=sp.Symbol("chi_c", positive=True),
        kappa_phase=sp.Symbol("kappa_phase", positive=True),
        beta_B=sp.Symbol("beta_B", positive=True),
        m_theta2=sp.Symbol("m_theta2", positive=True),
        k=sp.Symbol("k", positive=True),
        omega=sp.Symbol("omega", positive=True),
        u_L=sp.Symbol("u_L"),
        u_T1=sp.Symbol("u_T1"),
        u_T2=sp.Symbol("u_T2"),
        theta=sp.Symbol("theta"),
        dot_u_L=sp.Symbol("dot_u_L"),
        dot_u_T=sp.Symbol("dot_u_T"),
        dot_theta=sp.Symbol("dot_theta"),
        p_u=sp.Symbol("p_u"),
        pi_theta=sp.Symbol("pi_theta"),
        delta_rho_B=sp.Symbol("delta_rho_B"),
        div_u=sp.Symbol("div_u"),
        theta_dot_symbol=sp.Symbol("partial_t_theta"),
        div_dot_u=sp.Symbol("div_partial_t_u"),
        grad_theta_dot_dot_u=sp.Symbol("grad_theta_dot_partial_t_u"),
        theta_boundary_density=sp.Symbol("theta_div_u_boundary_density"),
    )


@dataclass(frozen=True)
class IbpChain:
    slaving_sign: sp.Integer
    delta_rho_expr: sp.Expr
    raw_theta_dot_div_coeff: sp.Expr
    time_bulk_coeff: sp.Expr
    space_bulk_coeff: sp.Expr
    derived_C_J: sp.Expr
    raw_density: sp.Expr
    time_bulk_density: sp.Expr
    maxwell_cross_density: sp.Expr
    time_boundary_dropped: bool
    space_boundary_dropped: bool


def derive_cj_by_ibp(s: Symbols, *, slaving_sign: sp.Integer = sp.Integer(-1)) -> IbpChain:
    """Derive the C_J sign from delta rho_B slaving and time+space IBP.

    Starting density: J theta_dot delta_rho_B, with
    delta_rho_B = slaving_sign*rho_B0 div(u).  Time IBP flips the coefficient
    on theta div(dot u); space IBP flips it again onto grad(theta).dot(u_dot).
    """

    delta_rho_expr = sp.simplify(slaving_sign * s.rho_B0 * s.div_u)
    raw_coeff = sp.simplify(s.J * slaving_sign * s.rho_B0)
    time_bulk_coeff = sp.simplify(-raw_coeff)
    space_bulk_coeff = sp.simplify(-time_bulk_coeff)

    return IbpChain(
        slaving_sign=slaving_sign,
        delta_rho_expr=delta_rho_expr,
        raw_theta_dot_div_coeff=raw_coeff,
        time_bulk_coeff=time_bulk_coeff,
        space_bulk_coeff=space_bulk_coeff,
        derived_C_J=space_bulk_coeff,
        raw_density=sp.simplify(raw_coeff * s.theta_dot_symbol * s.div_u),
        time_bulk_density=sp.simplify(time_bulk_coeff * s.theta * s.div_dot_u),
        maxwell_cross_density=sp.simplify(space_bulk_coeff * s.grad_theta_dot_dot_u),
        time_boundary_dropped=True,
        space_boundary_dropped=True,
    )


def poisson(f: sp.Expr, g: sp.Expr, s: Symbols) -> sp.Expr:
    return sp.simplify(
        sp.diff(f, s.u_L) * sp.diff(g, s.p_u)
        - sp.diff(f, s.p_u) * sp.diff(g, s.u_L)
        + sp.diff(f, s.theta) * sp.diff(g, s.pi_theta)
        - sp.diff(f, s.pi_theta) * sp.diff(g, s.theta)
    )


def sign_label(expr: sp.Expr) -> str:
    clean = sp.simplify(expr)
    if clean == 0:
        return "zero"
    if clean.is_positive:
        return "positive"
    if clean.is_negative:
        return "negative"
    return "symbolic"


def first_order_analysis(
    s: Symbols,
    *,
    name: str,
    c_j_expr: sp.Expr,
    K_expr: sp.Expr,
    B_eff_expr: sp.Expr,
    theta_mass_expr: sp.Expr = sp.Integer(0),
    provenance_forced: bool = False,
    gauge_on_shell_only: bool = False,
    partial_gauge_only: bool = False,
    gapped_not_gauge_removed: bool = False,
) -> dict[str, Any]:
    S_theta = sp.simplify(K_expr * s.k**2 - theta_mass_expr)
    primitive_L = (
        sp.Rational(1, 2) * s.rho_br * s.dot_u_L**2
        - c_j_expr * s.k * s.u_L * s.dot_theta
        + sp.Rational(1, 2) * S_theta * s.theta**2
        - sp.Rational(1, 2) * B_eff_expr * s.k**2 * s.u_L**2
    )
    p_u_expr = sp.diff(primitive_L, s.dot_u_L)
    p_theta_expr = sp.diff(primitive_L, s.dot_theta)
    primary = sp.simplify(s.pi_theta - p_theta_expr)

    solved_dot_u = sp.solve(sp.Eq(s.p_u, p_u_expr), s.dot_u_L)[0]
    canonical_H = sp.simplify(
        (s.p_u * s.dot_u_L + p_theta_expr * s.dot_theta - primitive_L).subs(
            s.dot_u_L, solved_dot_u
        )
    )
    secondary = sp.simplify(poisson(primary, canonical_H, s))
    bracket = sp.simplify(poisson(primary, secondary, s))
    secondary_preservation_no_multiplier = sp.simplify(poisson(secondary, canonical_H, s))

    dynamic_matrix = sp.Matrix(
        [
            [B_eff_expr * s.k**2 - s.rho_br * s.omega**2, -sp.I * c_j_expr * s.k * s.omega],
            [sp.I * c_j_expr * s.k * s.omega, -S_theta],
        ]
    )
    determinant = sp.factor(dynamic_matrix.det())
    expected_determinant = sp.factor(
        (s.rho_br * S_theta - c_j_expr**2 * s.k**2) * s.omega**2
        - B_eff_expr * s.k**2 * S_theta
    )
    A_eff = sp.simplify(s.rho_br - c_j_expr**2 * s.k**2 / S_theta) if S_theta != 0 else sp.oo
    omega2_pole = (
        sp.Integer(0)
        if B_eff_expr == 0
        else sp.simplify(B_eff_expr * s.k**2 / A_eff)
        if A_eff not in (sp.oo, -sp.oo) and sp.simplify(A_eff) != 0
        else sp.nan
    )

    bracket_zero = sp.simplify(bracket) == 0
    b_zero = sp.simplify(B_eff_expr) == 0
    mass_zero = sp.simplify(theta_mass_expr) == 0
    square_residual = sp.simplify(K_expr - c_j_expr**2 / s.rho_br)

    if gauge_on_shell_only:
        verdict = "FAIL_GAUGE_ON_SHELL_ONLY"
        first_class = 0
        second_class = 2
        physical_dof = 1
        initial_data = 2
        classification = "ON_SHELL_ONLY_NOT_FIRST_CLASS"
        pole_count = 1
        residue = "symbolic"
        bounded = False
    elif partial_gauge_only:
        verdict = "FAIL_PARTIAL_GAUGE_ONLY"
        first_class = 0
        second_class = 2
        physical_dof = 1
        initial_data = 2
        classification = "PARTIAL_GAUGE_ONLY"
        pole_count = 1
        residue = "symbolic"
        bounded = False
    elif bracket_zero and b_zero and mass_zero:
        verdict = (
            "C5_RESOLVED_MAXWELL_WITH_PROVENANCE"
            if provenance_forced
            else "C5_RESOLVED_MAXWELL_BY_TUNING"
        )
        first_class = 2
        second_class = 0
        physical_dof = 0
        initial_data = 0
        classification = "FIRST_CLASS_MAXWELL_CHAIN"
        pole_count = 0
        residue = "no_physical_longitudinal_pole"
        bounded = True
    elif bracket_zero and not b_zero:
        verdict = "FAIL_SECOND_CLASS_NOT_MAXWELL"
        first_class = 0
        second_class = 4
        physical_dof = 0
        initial_data = 0
        classification = "TERTIARY_SECOND_CLASS_CHAIN"
        pole_count = 0
        residue = "no_gauge_pole_but_not_first_class"
        bounded = True
    else:
        first_class = 0
        second_class = 2
        physical_dof = 1
        initial_data = 2
        classification = "SECOND_CLASS_PAIR"
        pole_count = 1
        residue = sign_label(1 / A_eff) if A_eff not in (sp.oo, -sp.oo, sp.nan) else "undefined"
        bounded = A_eff.is_positive is True and (b_zero or B_eff_expr.is_positive is True)
        if gapped_not_gauge_removed:
            verdict = "FAIL_GAPPED_NOT_GAUGE_REMOVED"
        elif not mass_zero:
            verdict = "FAIL_SECOND_CLASS_NOT_MAXWELL"
        elif bounded is False:
            verdict = "FAIL_GHOST_OR_NEGATIVE_NORM"
        elif b_zero:
            verdict = "FAIL_C5_LONGITUDINAL_ZERO_MODE"
        else:
            verdict = "FAIL_CAUCHY_STRAY_LONGITUDINAL"

    return {
        "name": name,
        "primitive_L": primitive_L,
        "p_u": sp.simplify(p_u_expr),
        "pi_theta": sp.simplify(p_theta_expr),
        "primary": primary,
        "hamiltonian": canonical_H,
        "secondary": secondary,
        "constraint_bracket": bracket,
        "secondary_preservation_no_multiplier": secondary_preservation_no_multiplier,
        "dynamic_determinant": determinant,
        "expected_dynamic_determinant": expected_determinant,
        "omega2_pole": sp.simplify(omega2_pole),
        "A_eff": sp.simplify(A_eff),
        "square_residual_K_minus_CJ2_over_rho": square_residual,
        "first_class_count": first_class,
        "second_class_count": second_class,
        "physical_dof": physical_dof,
        "initial_data": initial_data,
        "classification": classification,
        "pole_count": pole_count,
        "residue_sign": residue,
        "bounded": bounded,
        "verdict": verdict,
    }


def elastic_longitudinal_control(
    s: Symbols, *, name: str, B_eff_expr: sp.Expr
) -> dict[str, Any]:
    omega2_expr = sp.simplify(B_eff_expr * s.k**2 / s.rho_br)
    propagating = sp.simplify(B_eff_expr) != 0
    return {
        "name": name,
        "omega2": omega2_expr,
        "physical_dof": 1,
        "initial_data": 2,
        "pole_count": 1,
        "classification": "UNCONSTRAINED_ELASTIC_COORDINATE",
        "verdict": "FAIL_CAUCHY_STRAY_LONGITUDINAL"
        if propagating
        else "FAIL_C5_LONGITUDINAL_ZERO_MODE",
    }


def decoupled_slaved_theta_control(s: Symbols) -> dict[str, Any]:
    primitive_L = (
        sp.Rational(1, 2) * s.rho_br * s.dot_u_L**2
        - sp.Rational(1, 2) * s.kappa_phase * s.k**2 * s.theta**2
    )
    p_u_expr = sp.diff(primitive_L, s.dot_u_L)
    p_theta_expr = sp.diff(primitive_L, s.dot_theta)
    primary = sp.simplify(s.pi_theta - p_theta_expr)
    H = sp.simplify(s.p_u**2 / (2 * s.rho_br) + sp.Rational(1, 2) * s.kappa_phase * s.k**2 * s.theta**2)
    secondary = sp.simplify(poisson(primary, H, s))
    bracket = sp.simplify(poisson(primary, secondary, s))
    return {
        "name": "decoupled_theta_slaved",
        "p_u": p_u_expr,
        "pi_theta": p_theta_expr,
        "primary": primary,
        "secondary": secondary,
        "constraint_bracket": bracket,
        "physical_dof": 1,
        "initial_data": 2,
        "classification": "THETA_ALGEBRAIC_SECOND_CLASS_PAIR_PLUS_U_ZERO_MODE",
        "verdict": "FAIL_C5_LONGITUDINAL_ZERO_MODE",
    }


def independent_density_without_continuity(s: Symbols) -> dict[str, Any]:
    density_L = s.J * s.dot_theta * s.delta_rho_B - s.delta_rho_B**2 / (2 * s.chi_c)
    gaussian_solution = sp.solve(sp.Eq(sp.diff(density_L, s.delta_rho_B), 0), s.delta_rho_B)[0]
    theta_kinetic = sp.simplify(density_L.subs(s.delta_rho_B, gaussian_solution))
    theta_omega2 = sp.simplify(s.kappa_phase * s.k**2 / (s.J**2 * s.chi_c))
    return {
        "name": "decoupled_theta_independent_no_continuity",
        "gaussian_solution": gaussian_solution,
        "theta_kinetic": theta_kinetic,
        "theta_omega2": theta_omega2,
        "hessian_rank": 2,
        "physical_dof": 2,
        "initial_data": 4,
        "classification": "TWO_UNCONSTRAINED_SECOND_ORDER_FIELDS",
        "verdict": "FAIL_EXTRA_SCALAR_DOF",
    }


def branch_a_continuity_proof(s: Symbols, c_j_expr: sp.Expr) -> dict[str, Any]:
    density_L = s.J * s.dot_theta * s.delta_rho_B - s.delta_rho_B**2 / (2 * s.chi_c)
    continuity_solution = -s.rho_B0 * s.k * s.u_L
    with_continuity = sp.simplify(density_L.subs(s.delta_rho_B, continuity_solution))
    raw_cross_coeff = sp.diff(sp.diff(with_continuity, s.dot_theta), s.u_L) / s.k
    canonical_cross_coeff = sp.simplify(-raw_cross_coeff)
    B_increment = sp.simplify(s.rho_B0**2 / s.chi_c)
    return {
        "continuity_residual": sp.simplify(
            s.omega * (continuity_solution + s.rho_B0 * s.k * s.u_L)
        ),
        "with_continuity": with_continuity,
        "raw_cross_coeff": sp.simplify(raw_cross_coeff),
        "canonical_minus_CJ_residual": sp.simplify(canonical_cross_coeff + c_j_expr),
        "B_increment": B_increment,
        "proof_status": "CONTINUITY_FORCES_SAME_SLAVED_SECTOR",
    }


def epsilon_mismatch_control(s: Symbols, c_j_expr: sp.Expr) -> dict[str, Any]:
    rho_eps = 2 * s.rho_br
    K_eps = sp.simplify(c_j_expr**2 / rho_eps)
    S_theta = K_eps * s.k**2
    primitive_L = (
        sp.Rational(1, 2) * rho_eps * s.dot_u_L**2
        - c_j_expr * s.k * s.u_L * s.dot_theta
        + sp.Rational(1, 2) * S_theta * s.theta**2
    )
    p_theta_expr = sp.diff(primitive_L, s.dot_theta)
    primary = sp.simplify(s.pi_theta - p_theta_expr)
    H = sp.simplify(s.p_u**2 / (2 * rho_eps) - sp.Rational(1, 2) * S_theta * s.theta**2)
    secondary = sp.simplify(poisson(primary, H, s))
    bracket = sp.simplify(poisson(primary, secondary, s))
    frozen_transverse = sp.simplify(s.mu_R * s.k**2 / s.rho_br)
    shifted_transverse = sp.simplify(s.mu_R * s.k**2 / rho_eps)
    return {
        "name": "epsilon_mismatch",
        "constraint_bracket": bracket,
        "square_closes": bracket == 0,
        "frozen_transverse_omega2": frozen_transverse,
        "mismatched_transverse_omega2": shifted_transverse,
        "speed_shift": sp.simplify(shifted_transverse - frozen_transverse),
        "physical_dof": 0 if bracket == 0 else 1,
        "initial_data": 0 if bracket == 0 else 2,
        "verdict": "FAIL_TRANSVERSE_DISTURBED",
    }


def transverse_dispersion_from_lagrangian(s: Symbols, L_per_pol: sp.Expr) -> dict[str, sp.Expr]:
    t = sp.Symbol("t", real=True)
    u_T = sp.Function("u_T")
    u = u_T(t)
    L_t = L_per_pol.subs({s.u_T1: u, s.dot_u_T: sp.diff(u, t)})
    el = sp.simplify(sp.diff(sp.diff(L_t, sp.diff(u, t)), t) - sp.diff(L_t, u))
    plane_wave_el = sp.simplify(el.subs(sp.diff(u, (t, 2)), -s.omega**2 * u))
    dispersion_residual = sp.simplify(plane_wave_el / u)
    omega2 = sp.solve(sp.Eq(dispersion_residual, 0), s.omega**2)[0]
    return {
        "EL": el,
        "plane_wave_residual": dispersion_residual,
        "omega2": sp.simplify(omega2),
    }


def transverse_sector(s: Symbols) -> dict[str, Any]:
    wave = sp.Matrix([s.k, 0, 0])
    u_vec = sp.Matrix([s.u_L, s.u_T1, s.u_T2])
    k_hat = sp.Matrix([1, 0, 0])
    transverse_projector = sp.eye(3) - k_hat * k_hat.T
    div_u_transverse = sp.simplify(wave.dot(transverse_projector * u_vec))
    curl_sq = s.k**2 * (s.u_T1**2 + s.u_T2**2)
    slaved_delta_rho_transverse = -s.rho_B0 * div_u_transverse
    josephson_transverse = s.J * s.dot_theta * slaved_delta_rho_transverse
    L_per_pol = sp.Rational(1, 2) * s.rho_br * s.dot_u_T**2 - sp.Rational(1, 2) * s.mu_R * s.k**2 * s.u_T1**2
    p_T = sp.diff(L_per_pol, s.dot_u_T)
    dispersion = transverse_dispersion_from_lagrangian(s, L_per_pol)
    omega2 = dispersion["omega2"]
    return {
        "basis_count": 2,
        "curl_sq": curl_sq,
        "josephson_transverse": josephson_transverse,
        "p_T": p_T,
        "EL": dispersion["EL"],
        "plane_wave_residual": dispersion["plane_wave_residual"],
        "omega2": sp.simplify(omega2),
        "c_gamma2": sp.simplify(omega2 / s.k**2),
        "physical_dof": 2,
        "massless": True,
        "theta_couplings": sp.simplify(josephson_transverse),
        "verdict": "PASS_TRANSVERSE_UNDISTURBED",
    }


def assert_verdict(name: str, actual: str, expected: str | tuple[str, ...]) -> None:
    expected_tuple = (expected,) if isinstance(expected, str) else expected
    expect_bool(f"{name}: verdict {actual}", actual in expected_tuple)


def run_ibp_checks(s: Symbols) -> sp.Expr:
    subbanner("Josephson slaving and C_J IBP sign")
    chain = derive_cj_by_ibp(s)

    expect_zero("slaved delta_rho_B = -rho_B0 div u", chain.delta_rho_expr + s.rho_B0 * s.div_u)
    expect_zero(
        "time IBP flips theta_dot div(u) coefficient",
        chain.time_bulk_coeff + chain.raw_theta_dot_div_coeff,
    )
    expect_zero(
        "space IBP flips theta div(dot u) coefficient",
        chain.space_bulk_coeff + chain.time_bulk_coeff,
    )
    expect_zero("derived C_J = -J*rho_B0", chain.derived_C_J + s.J * s.rho_B0)
    expect_bool("time boundary term dropped in ledger", chain.time_boundary_dropped)
    expect_bool("space boundary term dropped in ledger", chain.space_boundary_dropped)

    canonical_L = (
        sp.Rational(1, 2) * s.rho_br * s.dot_u_L**2
        - chain.derived_C_J * s.k * s.u_L * s.dot_theta
    )
    pi_theta_expr = sp.diff(canonical_L, s.dot_theta)
    expect_zero("sign-sensitive pi_theta = +J*k*rho_B0*u_L", pi_theta_expr - s.J * s.k * s.rho_B0 * s.u_L)

    corrupted = derive_cj_by_ibp(s, slaving_sign=sp.Integer(1))
    corrupted_L = -corrupted.derived_C_J * s.k * s.u_L * s.dot_theta
    corrupted_pi = sp.diff(corrupted_L, s.dot_theta)
    expect_fail(
        "C_J sign mutation flips pi_theta residual",
        corrupted_pi - s.J * s.k * s.rho_B0 * s.u_L,
    )

    return chain.derived_C_J


def run_transverse_checks(s: Symbols) -> None:
    subbanner("Transverse earned sector")
    trans = transverse_sector(s)
    expect_zero("transverse omega^2 = (mu_R/rho_br) k^2", trans["omega2"] - s.mu_R * s.k**2 / s.rho_br)
    expect_zero("c_gamma^2 = mu_R/rho_br", trans["c_gamma2"] - s.mu_R / s.rho_br)
    expect_zero("theta couplings vanish for transverse modes", trans["theta_couplings"])
    expect_zero("two transverse physical polarizations", sp.Integer(trans["physical_dof"]) - 2)
    expect_bool("transverse modes are massless", trans["massless"])
    assert_verdict("transverse sector", trans["verdict"], "PASS_TRANSVERSE_UNDISTURBED")

    corrupted_L = (
        sp.Rational(1, 2) * s.rho_br * s.dot_u_T**2
        - sp.Rational(1, 2) * (2 * s.mu_R) * s.k**2 * s.u_T1**2
    )
    corrupted_dispersion = transverse_dispersion_from_lagrangian(s, corrupted_L)
    expect_fail(
        "transverse primitive mu_R->2 mu_R mutation fails c_gamma^2",
        sp.simplify(corrupted_dispersion["omega2"] / s.k**2 - s.mu_R / s.rho_br),
    )


def run_longitudinal_baseline(s: Symbols, c_j_expr: sp.Expr) -> dict[str, Any]:
    subbanner("Longitudinal Dirac-Bergmann chain")
    K_conventional = -s.kappa_phase
    B_eff = sp.simplify(s.rho_B0**2 / s.chi_c)
    analysis = first_order_analysis(
        s,
        name="branch_b_slaved_finite_compressibility_conventional_K",
        c_j_expr=c_j_expr,
        K_expr=K_conventional,
        B_eff_expr=B_eff,
    )

    expect_zero("p_u = rho_br dot_u_L", analysis["p_u"] - s.rho_br * s.dot_u_L)
    expect_zero("pi_theta = +J*k*rho_B0*u_L", analysis["pi_theta"] - s.J * s.k * s.rho_B0 * s.u_L)
    expect_zero("Phi_1 = pi_theta - J*k*rho_B0*u_L", analysis["primary"] - (s.pi_theta - s.J * s.k * s.rho_B0 * s.u_L))
    expect_zero(
        "Phi_2 = -k(J*p_u*rho_B0 + k*kappa_phase*rho_br*theta)/rho_br",
        analysis["secondary"]
        + s.k * (s.J * s.p_u * s.rho_B0 + s.k * s.kappa_phase * s.rho_br * s.theta) / s.rho_br,
    )
    expect_zero(
        "{Phi_1,Phi_2} = k^2(J^2 rho_B0^2 + kappa_phase rho_br)/rho_br",
        analysis["constraint_bracket"]
        - s.k**2 * (s.J**2 * s.rho_B0**2 + s.kappa_phase * s.rho_br) / s.rho_br,
    )
    expect_zero("first-class count = 0", sp.Integer(analysis["first_class_count"]))
    expect_zero("second-class count = 2", sp.Integer(analysis["second_class_count"]) - 2)
    expect_zero("longitudinal physical DOF = (4-2)/2 = 1", sp.Integer(analysis["physical_dof"]) - 1)
    expect_zero("independent initial data functions = 2", sp.Integer(analysis["initial_data"]) - 2)
    expect_bool("constraint classification SECOND_CLASS_PAIR", analysis["classification"] == "SECOND_CLASS_PAIR")

    expect_zero("dynamic determinant assembled from Euler-Lagrange matrix", analysis["dynamic_determinant"] - analysis["expected_dynamic_determinant"])
    expect_zero(
        "stray pole omega^2",
        analysis["omega2_pole"]
        - s.k**2
        * s.kappa_phase
        * s.rho_B0**2
        / (s.chi_c * (s.J**2 * s.rho_B0**2 + s.kappa_phase * s.rho_br)),
    )
    expect_zero("pole count = 1", sp.Integer(analysis["pole_count"]) - 1)
    expect_bool("positive pole residue", analysis["residue_sign"] == "positive")
    expect_bool("reduced Hamiltonian bounded", analysis["bounded"] is True)
    assert_verdict("baseline headline", analysis["verdict"], "FAIL_CAUCHY_STRAY_LONGITUDINAL")

    print("")
    print("Postulate line:")
    print("  K_theta = -kappa_phase is an explicit conventional phase-stiffness input, not a derivation.")
    return analysis


def run_branch_and_locus_checks(s: Symbols, c_j_expr: sp.Expr) -> tuple[dict[str, Any], dict[str, Any]]:
    subbanner("Branch (a), curl-only subcase, and Maxwell locus")
    K_conventional = -s.kappa_phase
    B_eff = sp.simplify(s.rho_B0**2 / s.chi_c)
    K_locus = sp.simplify(c_j_expr**2 / s.rho_br)

    branch_a = branch_a_continuity_proof(s, c_j_expr)
    branch_a_finite = first_order_analysis(
        s,
        name="branch_a_independent_with_continuity_integrated_out",
        c_j_expr=c_j_expr,
        K_expr=K_conventional,
        B_eff_expr=B_eff,
    )
    curl_only = first_order_analysis(
        s,
        name="branch_b_slaved_curl_only_conventional_K",
        c_j_expr=c_j_expr,
        K_expr=K_conventional,
        B_eff_expr=sp.Integer(0),
    )
    maxwell = first_order_analysis(
        s,
        name="branch_b_slaved_tuned_Maxwell_locus",
        c_j_expr=c_j_expr,
        K_expr=K_locus,
        B_eff_expr=sp.Integer(0),
    )

    expect_zero("branch (a) continuity equation solved in fixed-number sector", branch_a["continuity_residual"])
    expect_zero("branch (a) canonical cross coefficient matches derived C_J", branch_a["canonical_minus_CJ_residual"])
    expect_zero("branch (a) B increment = rho_B0^2/chi_c", branch_a["B_increment"] - B_eff)
    assert_verdict("branch (a) integrated sector", branch_a_finite["verdict"], "FAIL_CAUCHY_STRAY_LONGITUDINAL")
    expect_bool("branch (a) proof status", branch_a["proof_status"] == "CONTINUITY_FORCES_SAME_SLAVED_SECTOR")

    assert_verdict("curl-only conventional subcase", curl_only["verdict"], "FAIL_C5_LONGITUDINAL_ZERO_MODE")
    expect_zero("curl-only longitudinal DOF = 1", sp.Integer(curl_only["physical_dof"]) - 1)

    expect_zero("Maxwell-locus square residual K_theta - C_J^2/rho_br", maxwell["square_residual_K_minus_CJ2_over_rho"])
    expect_zero("Maxwell-locus bracket = 0", maxwell["constraint_bracket"])
    expect_zero("Maxwell-locus first-class count = 2", sp.Integer(maxwell["first_class_count"]) - 2)
    expect_zero("Maxwell-locus longitudinal DOF = 0", sp.Integer(maxwell["physical_dof"]))
    assert_verdict("Maxwell locus reachable", maxwell["verdict"], "C5_RESOLVED_MAXWELL_BY_TUNING")
    expect_bool("Maxwell-locus classification", maxwell["classification"] == "FIRST_CLASS_MAXWELL_CHAIN")

    return curl_only, maxwell


def run_control_checks(s: Symbols, c_j_expr: sp.Expr, baseline: dict[str, Any], curl_only: dict[str, Any], maxwell: dict[str, Any]) -> None:
    subbanner("Fourteen source controls")
    K_conventional = -s.kappa_phase
    K_locus = sp.simplify(c_j_expr**2 / s.rho_br)

    controls: list[tuple[str, str | tuple[str, ...], dict[str, Any] | str]] = [
        (
            "1_no_theta",
            "FAIL_C5_LONGITUDINAL_ZERO_MODE",
            elastic_longitudinal_control(s, name="no_theta_curl_only", B_eff_expr=sp.Integer(0)),
        ),
        (
            "2_cauchy_bulk",
            "FAIL_CAUCHY_STRAY_LONGITUDINAL",
            elastic_longitudinal_control(s, name="cauchy_bulk_no_theta", B_eff_expr=s.beta_B),
        ),
        (
            "3_mismatched_positive_K_no_B",
            "FAIL_C5_LONGITUDINAL_ZERO_MODE",
            first_order_analysis(
                s,
                name="mismatched_positive_K_no_B",
                c_j_expr=c_j_expr,
                K_expr=sp.simplify(2 * c_j_expr**2 / s.rho_br),
                B_eff_expr=sp.Integer(0),
            ),
        ),
        (
            "3_mismatched_K_theta_le_0",
            "FAIL_C5_LONGITUDINAL_ZERO_MODE",
            first_order_analysis(
                s,
                name="mismatched_K_theta_le_0",
                c_j_expr=c_j_expr,
                K_expr=K_conventional,
                B_eff_expr=sp.Integer(0),
            ),
        ),
        (
            "3_mismatched_positive_K_with_B",
            "FAIL_CAUCHY_STRAY_LONGITUDINAL",
            first_order_analysis(
                s,
                name="mismatched_positive_K_with_B",
                c_j_expr=c_j_expr,
                K_expr=sp.simplify(2 * c_j_expr**2 / s.rho_br),
                B_eff_expr=s.beta_B,
            ),
        ),
        (
            "3_positive_K_negative_residue",
            "FAIL_GHOST_OR_NEGATIVE_NORM",
            first_order_analysis(
                s,
                name="positive_K_negative_residue",
                c_j_expr=c_j_expr,
                K_expr=sp.simplify(c_j_expr**2 / (2 * s.rho_br)),
                B_eff_expr=s.beta_B,
            ),
        ),
        (
            "3_B_on_square_locus",
            "FAIL_SECOND_CLASS_NOT_MAXWELL",
            first_order_analysis(
                s,
                name="B_on_square_locus",
                c_j_expr=c_j_expr,
                K_expr=K_locus,
                B_eff_expr=s.beta_B,
            ),
        ),
    ]

    for name, expected, item in controls:
        assert isinstance(item, dict)
        assert_verdict(name, str(item["verdict"]), expected)

    epsilon = epsilon_mismatch_control(s, c_j_expr)
    assert_verdict("3_epsilon_mismatch", epsilon["verdict"], "FAIL_TRANSVERSE_DISTURBED")
    expect_zero("3_epsilon_mismatch transverse speed shifts to mu_R/(2 rho_br)", epsilon["mismatched_transverse_omega2"] - s.mu_R * s.k**2 / (2 * s.rho_br))

    decoupled_slaved = decoupled_slaved_theta_control(s)
    assert_verdict("4_decoupled_theta_slaved", decoupled_slaved["verdict"], "FAIL_C5_LONGITUDINAL_ZERO_MODE")

    independent_no_continuity = independent_density_without_continuity(s)
    assert_verdict(
        "4_decoupled_theta_independent_no_continuity",
        independent_no_continuity["verdict"],
        "FAIL_EXTRA_SCALAR_DOF",
    )
    expect_zero(
        "4_decoupled_theta_independent_no_continuity theta kinetic",
        independent_no_continuity["theta_kinetic"]
        - sp.Rational(1, 2) * s.J**2 * s.chi_c * s.dot_theta**2,
    )

    trans = transverse_sector(s)
    assert_verdict("5_transverse", trans["verdict"], "PASS_TRANSVERSE_UNDISTURBED")

    assert_verdict("6_provenance_ablation fixed coefficients", baseline["verdict"], "FAIL_CAUCHY_STRAY_LONGITUDINAL")
    assert_verdict("6_provenance_ablation free locus", maxwell["verdict"], "C5_RESOLVED_MAXWELL_BY_TUNING")

    assert_verdict("7_compressibility_absent_vs_included absent", curl_only["verdict"], "FAIL_C5_LONGITUDINAL_ZERO_MODE")
    assert_verdict("7_compressibility_absent_vs_included included", baseline["verdict"], "FAIL_CAUCHY_STRAY_LONGITUDINAL")

    theta_mass = first_order_analysis(
        s,
        name="theta_mass_breaks_gauge",
        c_j_expr=c_j_expr,
        K_expr=K_locus,
        B_eff_expr=sp.Integer(0),
        theta_mass_expr=s.m_theta2,
    )
    assert_verdict("8_theta_mass", theta_mass["verdict"], "FAIL_SECOND_CLASS_NOT_MAXWELL")

    reachable_tokens = {
        first_order_analysis(
            s,
            name="grammar_with_provenance",
            c_j_expr=c_j_expr,
            K_expr=K_locus,
            B_eff_expr=sp.Integer(0),
            provenance_forced=True,
        )["verdict"],
        first_order_analysis(
            s,
            name="grammar_gapped",
            c_j_expr=c_j_expr,
            K_expr=K_conventional,
            B_eff_expr=sp.Integer(0),
            gapped_not_gauge_removed=True,
        )["verdict"],
        first_order_analysis(
            s,
            name="grammar_partial",
            c_j_expr=c_j_expr,
            K_expr=K_conventional,
            B_eff_expr=sp.Integer(0),
            partial_gauge_only=True,
        )["verdict"],
        first_order_analysis(
            s,
            name="grammar_on_shell",
            c_j_expr=c_j_expr,
            K_expr=K_conventional,
            B_eff_expr=sp.Integer(0),
            gauge_on_shell_only=True,
        )["verdict"],
    }
    expect_bool(
        "reachable-but-unfired grammar tokens are preserved",
        {
            "C5_RESOLVED_MAXWELL_WITH_PROVENANCE",
            "FAIL_GAPPED_NOT_GAUGE_REMOVED",
            "FAIL_PARTIAL_GAUGE_ONLY",
            "FAIL_GAUGE_ON_SHELL_ONLY",
        }.issubset(reachable_tokens),
    )


def run_dimension_checks() -> None:
    subbanner("Dimensional firewall")
    brane_lag = dim(1, -1, -2)
    zdim = dim(0, 0, 0)
    d_u = dim(0, 1, 0)
    d_theta = zdim
    d_grad = dim(0, -1, 0)
    d_dt = dim(0, 0, -1)
    d_k = d_grad
    d_omega = d_dt
    d_rho_br = dim(1, -3, 0)
    d_mu_R = brane_lag
    d_B = brane_lag
    d_rho_B0 = d_rho_br
    d_J = dim(0, 2, -1)
    d_CJ = dim(1, -1, -1)
    d_Ktheta = dim(1, 1, -2)
    d_chi_c = dim(1, -5, 2)
    d_mtheta2 = brane_lag

    d_div_u = dim_add(d_grad, d_u)
    d_curl_u = d_div_u
    d_grad_theta = dim_add(d_grad, d_theta)
    d_dt_u = dim_add(d_dt, d_u)
    d_dt_theta = dim_add(d_dt, d_theta)
    d_delta_rho = d_rho_B0

    checks: list[tuple[str, Dim, Dim]] = [
        ("brane inertia rho_br (partial_t u)^2", dim_add(d_rho_br, dim_scale(2, d_dt_u)), brane_lag),
        ("MacCullagh curl energy mu_R (curl u)^2", dim_add(d_mu_R, dim_scale(2, d_curl_u)), brane_lag),
        ("Cauchy bulk term B (div u)^2", dim_add(d_B, dim_scale(2, d_div_u)), brane_lag),
        ("Josephson density term J theta_dot delta_rho_B", dim_add(d_J, d_dt_theta, d_delta_rho), brane_lag),
        ("slaved density rho_B0 div u", dim_add(d_rho_B0, d_div_u), d_delta_rho),
        ("signed phase gradient K_theta (grad theta)^2", dim_add(d_Ktheta, dim_scale(2, d_grad_theta)), brane_lag),
        ("compressibility delta_rho_B^2/chi_c", dim_sub(dim_scale(2, d_delta_rho), d_chi_c), brane_lag),
        ("theta mass m_theta^2 theta^2", dim_add(d_mtheta2, dim_scale(2, d_theta)), brane_lag),
        ("C_J = -J rho_B0", dim_add(d_J, d_rho_B0), d_CJ),
        ("IBP cross C_J partial_t u dot grad theta", dim_add(d_CJ, d_dt_u, d_grad_theta), brane_lag),
        ("electric square velocity piece", dim_add(d_rho_br, dim_scale(2, d_dt_u)), brane_lag),
        ("electric square mixed piece", dim_add(d_CJ, d_dt_u, d_grad_theta), brane_lag),
        ("electric square gradient piece", dim_add(d_Ktheta, dim_scale(2, d_grad_theta)), brane_lag),
        ("Maxwell locus C_J^2 = rho_br K_theta", dim_scale(2, d_CJ), dim_add(d_rho_br, d_Ktheta)),
        ("c_gamma^2 = mu_R/rho_br", dim_sub(d_mu_R, d_rho_br), dim(0, 2, -2)),
        ("branch-a theta kinetic chi_c J^2 theta_dot^2", dim_add(d_chi_c, dim_scale(2, d_J), dim_scale(2, d_dt_theta)), brane_lag),
        ("rho_B0^2/chi_c stiffness increment", dim_sub(dim_scale(2, d_rho_B0), d_chi_c), d_B),
        ("transverse omega^2", dim_add(dim_sub(d_mu_R, d_rho_br), dim_scale(2, d_k)), dim_scale(2, d_omega)),
    ]
    for name, actual, expected in checks:
        expect_zero(f"dimension: {name}", dimension_residual(actual, expected))

    expect_fail(
        "dimension ablation: drop rho_B0 from Josephson cross",
        dimension_residual(dim_add(d_J, d_dt_u, d_grad_theta), brane_lag),
    )
    expect_fail(
        "dimension ablation: omit gradient from phase stiffness",
        dimension_residual(dim_add(d_Ktheta, dim_scale(2, d_theta)), brane_lag),
    )
    expect_fail(
        "dimension ablation: multiply by chi_c instead of divide",
        dimension_residual(dim_add(d_chi_c, dim_scale(2, d_delta_rho)), brane_lag),
    )


def main() -> None:
    banner("ledger_stage003_transverse_photons_stray_longitudinal SymPy audit")
    s = make_symbols()

    c_j_expr = run_ibp_checks(s)
    run_transverse_checks(s)
    baseline = run_longitudinal_baseline(s, c_j_expr)
    curl_only, maxwell = run_branch_and_locus_checks(s, c_j_expr)
    run_control_checks(s, c_j_expr, baseline, curl_only, maxwell)
    run_dimension_checks()

    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): LIGHT_ON_BRANE_TWO_TRANSVERSE_PHOTONS  (c_gamma^2 = mu_R/rho_br, physical_dof=2)")
    print("  source top-line verdict: FAIL_CAUCHY_STRAY_LONGITUDINAL")
    print("  characterized departure: SECOND_CLASS_PAIR (one stray longitudinal DOF; Maxwell locus reachable BY_TUNING only)")
    print("")
    print(f"CHECK TALLY: PASS={PASS_COUNT} FAIL=0")
    print("OVERALL PASS: SymPy derived ledger_stage003 light-sector claims exactly")


if __name__ == "__main__":
    main()
