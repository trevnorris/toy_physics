#!/usr/bin/env python3
"""Stage 044 dual-engine leg: conservative parent-action reconciliation.

This audit is deliberately self-contained.  Ratified source facts are represented
as structured manifests, while physics identities are re-derived with SymPy.
The production verdict is selected by documented first-match over the ordered
tooth results and is written only on a clean, unmutated run.

Mutation protocol
-----------------
``LEDGER_STAGE044_MUTATION`` selects one primitive corruption.  Every registered
corruption must terminate non-zero at its mapped tooth and print
``FIRED_AT_OWN_ASSERT=<tooth>``.  ``--mutation-map`` exposes the complete
mutation-to-tooth map to the ablation runner.
"""

from __future__ import annotations

from collections import OrderedDict
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Callable, Iterable

import sympy as sp


MUTATION_ENV = "LEDGER_STAGE044_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()
PASS_COUNT = 0
FAIL_COUNT = 0
CHECK_RESULTS: OrderedDict[str, bool] = OrderedDict()
EVIDENCE: dict[str, Any] = {}

TOOTH_ORDER = (
    "A_WALL_STATIC_ENERGY_DEDUP",
    "B_WALL_HESSIAN_ZERO_MODE",
    "C1_CHARGE_FACTORIZATION_ZERO_MODE",
    "C2_CHARGE_ZERO_MODE_NORM",
    "D1_BULK_SOUND_SPEED",
    "D2_SCALAR_WAVE_SPEEDS",
    "D3_LIGHT_WAVE_SPEED",
    "E1_SCALAR_SYLVESTER",
    "E2_WALL_MINIMA",
    "E3_TRANSVERSE_GAP",
    "F_DIMENSIONAL_HOMOGENEITY",
    "G_SOURCE_FREE_CONTINUITY",
    "H_DRAIN_PROVENANCE",
    "I_P_RETIREMENT",
    "J_FIELD_UNION_SCHUR",
    "K1_DEFERRED_GATE_COMPLETENESS",
    "K2_DRAFT_OPEN_ACTION",
    "L_ACTION_INCIDENCE",
    "M_PROVENANCE_STATUS_BINDING",
    "REPRODUCTION",
)

MUTATION_TO_TOOTH: OrderedDict[str, str] = OrderedDict([
    ("A_MAP_LAMBDA_2A", "A_WALL_STATIC_ENERGY_DEDUP"),
    ("A_FALSE_FULL_FUNCTIONAL_IDENTITY", "A_WALL_STATIC_ENERGY_DEDUP"),
    ("B_WRONG_A_CHI_SIGN", "B_WALL_HESSIAN_ZERO_MODE"),
    ("C_PERTURB_F0", "C1_CHARGE_FACTORIZATION_ZERO_MODE"),
    ("C_NORM_TARGET_7_OVER_3ELL", "C2_CHARGE_ZERO_MODE_NORM"),
    ("D_BULK_PERTURB_RHO0", "D1_BULK_SOUND_SPEED"),
    ("D_SCALAR_PERTURB_K", "D2_SCALAR_WAVE_SPEEDS"),
    ("D_LIGHT_PERTURB_MU_R", "D3_LIGHT_WAVE_SPEED"),
    ("E_BEFF_PERTURB", "E1_SCALAR_SYLVESTER"),
    ("E_CHU_3_OVER_2", "E1_SCALAR_SYLVESTER"),
    ("E_FLIP_WALL_POTENTIAL", "E2_WALL_MINIMA"),
    ("E_NEGATIVE_OMEGA_W_SQ", "E3_TRANSVERSE_GAP"),
    ("F_CORRUPT_QT_DIMENSION", "F_DIMENSIONAL_HOMOGENEITY"),
    ("G_INJECT_DRAIN", "G_SOURCE_FREE_CONTINUITY"),
    ("H_GAMMA_B_VARIATIONAL", "H_DRAIN_PROVENANCE"),
    ("H_G0_DRAIN_VARIATIONAL", "H_DRAIN_PROVENANCE"),
    ("H_SUPPLY_EQUIVALENCE_MAP", "H_DRAIN_PROVENANCE"),
    ("H_SELECT_DRAIN_IN_044", "H_DRAIN_PROVENANCE"),
    ("I_RETIRED_IN_OPERATIVE", "I_P_RETIREMENT"),
    ("I_DOF_8_TO_5", "I_P_RETIREMENT"),
    ("I_STAGE007_DRIFT_11_TO_8", "I_P_RETIREMENT"),
    ("I_STAGE006_DRIFT_6_TO_6", "I_P_RETIREMENT"),
    ("J_RETAIN_THETA_B", "J_FIELD_UNION_SCHUR"),
    ("J_DROP_U_W", "J_FIELD_UNION_SCHUR"),
    ("J_CONFLATE_PHASE_INCIDENCE", "J_FIELD_UNION_SCHUR"),
    ("J_COUNT_AUX_AS_PHYSICAL_DOF", "J_FIELD_UNION_SCHUR"),
    ("J_DROP_TRANSVERSE_CONSTRAINT", "J_FIELD_UNION_SCHUR"),
    ("K_DROP_DEFERRED_GATE", "K1_DEFERRED_GATE_COMPLETENESS"),
    ("K_RESOLVE_ALL_GATES", "K2_DRAFT_OPEN_ACTION"),
    ("K_CORRUPT_CARD_HEADER", "K2_DRAFT_OPEN_ACTION"),
    ("L_DELETE_S_MOVE", "L_ACTION_INCIDENCE"),
    ("L_SUM_BOTH_SCALAR_BRANCHES", "L_ACTION_INCIDENCE"),
    ("M_G_ELL_WRONG_EDGE", "M_PROVENANCE_STATUS_BINDING"),
    ("M_LIGHT_GATING_WRONG_EDGE", "M_PROVENANCE_STATUS_BINDING"),
    ("M_Z_CHI_DROP_POSTULATE", "M_PROVENANCE_STATUS_BINDING"),
    ("M_Z_CHI_INVENT_ROUTE", "M_PROVENANCE_STATUS_BINDING"),
    ("M_Z_CHI_WRONG_ENDPOINT", "M_PROVENANCE_STATUS_BINDING"),
    ("M_GM_DROP_PN_CITATION", "M_PROVENANCE_STATUS_BINDING"),
    ("M_THROAT_WRONG_STATUS", "M_PROVENANCE_STATUS_BINDING"),
    ("M_EDGE_SEQUENCE_GAP", "M_PROVENANCE_STATUS_BINDING"),
    ("M_Z_MASTER_WRONG_DIM", "M_PROVENANCE_STATUS_BINDING"),
    ("REPRODUCTION_FIRST_MATCH", "REPRODUCTION"),
    ("REPRODUCTION_FIELD_BINDING", "REPRODUCTION"),
])

if len(TOOTH_ORDER) != 20:
    raise RuntimeError("stage044 detailed tooth declaration is not exactly 20")


class AuditFailure(AssertionError):
    """A named stage-044 predicate failed."""

    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


def canonical_json(value: Any) -> str:
    """RFC-8259 JSON with recursively sorted keys and compact separators."""

    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ) + "\n"


def section(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def expect_equal(name: str, actual: Any, expected: Any) -> None:
    global PASS_COUNT, FAIL_COUNT
    passed = actual == expected
    CHECK_RESULTS[name] = passed
    if passed:
        PASS_COUNT += 1
        print(f"PASS  {name}")
        return
    FAIL_COUNT += 1
    print(f"FIRST_FAILURE={name}")
    if MUTATION_TO_TOOTH.get(ACTIVE_MUTATION) == name:
        print(f"FIRED_AT_OWN_ASSERT={name}")
    print(f"FAIL  {name}")
    print("      actual   = " + repr(actual))
    print("      expected = " + repr(expected))
    raise AuditFailure(name)


def z(expr: sp.Expr) -> bool:
    return sp.simplify(sp.trigsimp(expr)) == 0


def exact_float(expr: sp.Expr) -> float:
    return float(sp.N(expr, 16))


def add_dim(*dims: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(sum(parts) for parts in zip(*dims))  # type: ignore[return-value]


def scale_dim(factor: int, dim: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(factor * value for value in dim)  # type: ignore[return-value]


def first_match(results: OrderedDict[str, bool]) -> str:
    """First failed tooth wins; only an all-clean sequence earns the floor token."""

    for name in TOOTH_ORDER[:-1]:
        if not results.get(name, False):
            return "ISSUES_FOUND_FIRST_" + name
    return "_".join(("PARENT_ACTION", "ASSEMBLED", "AT", "COMPLETENESS", "FLOOR"))


def wall_static_energy_tooth() -> dict[str, Any]:
    r, q = sp.symbols("r q", real=True)
    a_b, kappa_b, kappa_chi = sp.symbols(
        "a_B kappa_B kappa_chi", positive=True
    )
    lambda_chi = (
        2 * a_b if ACTIVE_MUTATION == "A_MAP_LAMBDA_2A" else 4 * a_b
    )
    ell = sp.sqrt(2 * kappa_chi / lambda_chi)
    delta = sp.sqrt(kappa_b / (2 * a_b))
    sigma_chi = kappa_chi / (6 * ell)
    sigma_wall = sp.sqrt(2 * a_b * kappa_b) / 6
    e_chi = kappa_chi * q**2 / 2 + lambda_chi * r**2 * (1 - r)**2 / 4
    e_part_i = kappa_b * q**2 / 2 + a_b * r**2 * (1 - r)**2
    mapped = {kappa_chi: kappa_b}
    v = a_b * r**2 * (1 - r)**2
    stationary = sp.solve(sp.diff(v, r), r)
    minima = [
        point for point in stationary
        if sp.simplify(sp.diff(v, r, 2).subs(r, point)) > 0
    ]
    g0_dynamics = {
        "type": "INERTIAL",
        "kinetic": "Z_chi*(dt r_B)^2/2",
    }
    part_i_dynamics = {
        "type": "DISSIPATIVE_ADJUNCT",
        "kinetic": None,
    }
    if ACTIVE_MUTATION == "A_FALSE_FULL_FUNCTIONAL_IDENTITY":
        part_i_dynamics = dict(g0_dynamics)
    full_functionals_equal = g0_dynamics == part_i_dynamics
    actual = {
        "ell_minus_delta": sp.simplify((ell - delta).subs(mapped)),
        "sigma_minus_wall": sp.simplify((sigma_chi - sigma_wall).subs(mapped)),
        "density_difference": sp.expand((e_chi - e_part_i).subs(mapped)),
        "minima": minima,
        "full_functionals_equal": full_functionals_equal,
    }
    expected = {
        "ell_minus_delta": sp.Integer(0),
        "sigma_minus_wall": sp.Integer(0),
        "density_difference": sp.Integer(0),
        "minima": [sp.Integer(0), sp.Integer(1)],
        "full_functionals_equal": False,
    }
    expect_equal("A_WALL_STATIC_ENERGY_DEDUP", actual, expected)
    EVIDENCE["A"] = actual
    return {
        "wall_dedup":
            "STATIC_WALL_ENERGY_IDENTITY(r_B==chi_B; "
            "kappa_chi==kappa_B; lambda_chi==4*a_B => "
            "ell==delta, sigma matches); "
            "FULL_FUNCTIONALS_DIFFER(inertial Z_chi vs dissipative M_chi)"
    }


def wall_hessian_tooth() -> dict[str, Any]:
    x, ell = sp.symbols("x ell", real=True, positive=True)
    sign = -1 if ACTIVE_MUTATION == "B_WRONG_A_CHI_SIGN" else 1
    W = sign * sp.tanh(x / (2 * ell)) / ell
    factor_potential = sp.simplify(W**2 - sp.diff(W, x))
    target_potential = (
        1 - sp.Rational(3, 2) * sp.sech(x / (2 * ell))**2
    ) / ell**2
    r0 = (1 + sp.tanh(x / (2 * ell))) / 2
    slope = sp.diff(r0, x)
    A_slope = sp.simplify(sp.diff(slope, x) + W * slope)
    actual = {
        "factorization_residual": sp.simplify(
            sp.trigsimp(factor_potential - target_potential)
        ),
        "A_chi_r0prime": sp.simplify(sp.trigsimp(A_slope)),
    }
    expected = {
        "factorization_residual": sp.Integer(0),
        "A_chi_r0prime": sp.Integer(0),
    }
    expect_equal("B_WALL_HESSIAN_ZERO_MODE", actual, expected)
    EVIDENCE["B"] = actual
    return {
        "wall_factorization":
            "L_chi_2ndVar/kappa_chi == A_chi_dag A_chi "
            "(annihilates r0')"
    }


def charge_factorization_teeth() -> dict[str, Any]:
    w, ell = sp.symbols("w ell", real=True, positive=True)
    W = 2 * sp.tanh(w / ell) / ell
    factor_potential = sp.simplify(W**2 - sp.diff(W, w))
    target_potential = (4 - 6 * sp.sech(w / ell)**2) / ell**2
    if ACTIVE_MUTATION == "C_PERTURB_F0":
        f0 = sp.sech(w / ell) / ell
    else:
        f0 = sp.sech(w / ell)**2 / ell
    o_f0 = sp.simplify(
        sp.trigsimp(-sp.diff(f0, w, 2) + target_potential * f0)
    )
    c1_actual = {
        "factorization_residual": sp.simplify(
            sp.trigsimp(factor_potential - target_potential)
        ),
        "O_perp_f0": o_f0,
    }
    c1_expected = {
        "factorization_residual": sp.Integer(0),
        "O_perp_f0": sp.Integer(0),
    }
    expect_equal("C1_CHARGE_FACTORIZATION_ZERO_MODE", c1_actual, c1_expected)

    # t=tanh(w/ell) turns sech^4(w/ell) dw into
    # ell*(1-t^2) dt, avoiding an unevaluated SymPy improper integral.
    t_sub = sp.symbols("t_sub", real=True)
    norm = sp.simplify(
        2 / ell * sp.integrate(1 - t_sub**2, (t_sub, -1, 1))
    )
    target_norm = (
        sp.Rational(7, 3) / ell
        if ACTIVE_MUTATION == "C_NORM_TARGET_7_OVER_3ELL"
        else sp.Rational(8, 3) / ell
    )
    c2_actual = sp.simplify(norm - target_norm)
    expect_equal("C2_CHARGE_ZERO_MODE_NORM", c2_actual, sp.Integer(0))
    EVIDENCE["C"] = {
        "factorization": c1_actual,
        "N0": norm,
    }
    return {"N0": norm}


def wave_speed_teeth() -> dict[str, Any]:
    rho, rho0, K_bulk, mass = sp.symbols(
        "rho rho0 K m", positive=True
    )
    U = K_bulk * rho**5 / 4
    background = 2 * rho0 if ACTIVE_MUTATION == "D_BULK_PERTURB_RHO0" else rho0
    c_s_sq = sp.simplify(background * sp.diff(U, rho, 2).subs(rho, background) / mass)
    bulk_actual = c_s_sq
    bulk_expected = 5 * K_bulk * rho0**4 / mass
    expect_equal("D1_BULK_SOUND_SPEED", bulk_actual, bulk_expected)

    b_entry = 3 if ACTIVE_MUTATION == "D_SCALAR_PERTURB_K" else 2
    K_matrix = sp.Matrix([[b_entry, sp.Rational(1, 2)],
                          [sp.Rational(1, 2), 1]])
    M_matrix = sp.eye(2)
    roots = sorted(
        sp.solve(sp.det(K_matrix - sp.Symbol("z") * M_matrix), sp.Symbol("z")),
        key=lambda value: float(sp.N(value)),
    )
    scalar_actual = roots
    scalar_expected = [
        (3 - sp.sqrt(2)) / 2,
        (3 + sp.sqrt(2)) / 2,
    ]
    expect_equal("D2_SCALAR_WAVE_SPEEDS", scalar_actual, scalar_expected)

    rho_br, mu_r, c2 = sp.symbols("rho_br mu_R c2", positive=True)
    stiffness = 2 * mu_r if ACTIVE_MUTATION == "D_LIGHT_PERTURB_MU_R" else mu_r
    light_dispersion = rho_br * c2 - stiffness
    light_speed = sp.solve(light_dispersion, c2)[0]
    expect_equal("D3_LIGHT_WAVE_SPEED", light_speed, mu_r / rho_br)
    EVIDENCE["D"] = {
        "bulk": c_s_sq,
        "scalar": roots,
        "light": light_speed,
    }
    return {
        "wave_speeds": {
            "c_s0_sq": "5*K*rho0^4/m",
            "c_pm_sq": [
                round(exact_float(scalar_expected[0]), 10),
                round(exact_float(scalar_expected[1]), 10),
            ],
            "c_gamma_sq": "mu_R/rho_br",
            "D_star": 1.75,
        }
    }


def stability_teeth() -> dict[str, Any]:
    b_eff = sp.Integer(-1) if ACTIVE_MUTATION == "E_BEFF_PERTURB" else sp.Integer(2)
    c_hu = (
        sp.Rational(3, 2)
        if ACTIVE_MUTATION == "E_CHU_3_OVER_2"
        else sp.Rational(1, 2)
    )
    scalar_k = sp.Matrix([[b_eff, c_hu], [c_hu, 1]])
    eigenvalues = list(scalar_k.eigenvals().keys())
    e1_actual = {
        "B_eff": scalar_k[0, 0],
        "D_star": sp.det(scalar_k),
        "positive_eigenvalues": all(value > 0 for value in eigenvalues),
    }
    e1_expected = {
        "B_eff": sp.Integer(2),
        "D_star": sp.Rational(7, 4),
        "positive_eigenvalues": True,
    }
    expect_equal("E1_SCALAR_SYLVESTER", e1_actual, e1_expected)

    r, a_b = sp.symbols("r a_B", real=True, positive=True)
    potential_sign = -1 if ACTIVE_MUTATION == "E_FLIP_WALL_POTENTIAL" else 1
    potential = potential_sign * a_b * r**2 * (1 - r)**2
    curvatures = [
        sp.simplify(sp.diff(potential, r, 2).subs(r, point))
        for point in (0, 1)
    ]
    e2_actual = {
        "curvatures": curvatures,
        "both_positive": all(curvature > 0 for curvature in curvatures),
    }
    e2_expected = {
        "curvatures": [2 * a_b, 2 * a_b],
        "both_positive": True,
    }
    expect_equal("E2_WALL_MINIMA", e2_actual, e2_expected)

    omega_sq = -sp.Integer(1) if ACTIVE_MUTATION == "E_NEGATIVE_OMEGA_W_SQ" else sp.Symbol(
        "Omega_w_sq", nonnegative=True
    )
    e3_actual = bool(omega_sq.is_nonnegative)
    expect_equal("E3_TRANSVERSE_GAP", e3_actual, True)
    EVIDENCE["E"] = {
        "scalar": e1_actual,
        "wall": e2_actual,
        "transverse_gap_nonnegative": e3_actual,
    }
    return {
        "stability":
            "POSITIVE_DEFINITE(B_eff>0 AND D_star=7/4>0; "
            "wall minima {0,1}; Omega_w^2>=0)"
    }


def dimensional_tooth() -> None:
    # Dimension order is [L,T,M].
    L = (1, 0, 0)
    inv_l = (-1, 0, 0)
    inv_t = (0, -1, 0)
    action = (2, -1, 1)
    energy = (2, -2, 1)
    dt_d4 = (4, 1, 0)
    dt_d3 = (3, 1, 0)
    dt_surface = (3, 1, 0)
    dt_only = (0, 1, 0)
    dims = {
        "hbar": (2, -1, 1),
        "rho": (-4, 0, 0),
        "m": (0, 0, 1),
        "U": (-2, -2, 1),
        "Z_chi": (-2, 0, 1),
        "kappa_chi": (0, -2, 1),
        "lambda_chi": (-2, -2, 1),
        "A_eff": (-3, 0, 1),
        "M_h": (-1, 0, 1),
        "B_eff": (-1, -2, 1),
        "K_h": (1, -2, 1),
        "C_hu": (0, -2, 1),
        "M4": (0, 0, 1),
        "K4": (2, -2, 1),
        "H": (-1, 0, 0),
        "K_m": (4, -2, 1),
        "J_m": (3, -2, 1),
        "eta": (-3, 0, 0),
        "lambda_Sigma": (-1, -2, 1),
        "E_g": energy,
        "g_ell": (-1, 0, 0),
        "rho_br": (-3, 0, 1),
        "mu_R": (-1, -2, 1),
        "Omega_w": (0, -1, 0),
        "u": L,
        "q_T": (0, -1, 1),
        "V": (1, -1, 0),
    }
    if ACTIVE_MUTATION == "F_CORRUPT_QT_DIMENSION":
        dims["q_T"] = (1, -1, 1)
    d_theta_t = inv_t
    d_grad_theta = inv_l
    d_grad_rho = add_dim(dims["rho"], inv_l)
    d_r_t = inv_t
    d_grad_r = inv_l
    d_u_t = add_dim(dims["u"], inv_t)
    d_grad_u = add_dim(dims["u"], inv_l)
    d_h_t = inv_t
    d_grad_h = inv_l
    d_H_t = add_dim(dims["H"], inv_t)
    d_grad_H = add_dim(dims["H"], inv_l)
    d_O_H = add_dim(dims["H"], scale_dim(2, inv_l))
    d_omega_u = d_grad_u
    integrand_dims = {
        "S_bulk.phase_time": add_dim(dims["hbar"], dims["rho"], d_theta_t),
        "S_bulk.phase_gradient": add_dim(
            scale_dim(2, dims["hbar"]), dims["rho"],
            scale_dim(-1, dims["m"]), scale_dim(2, d_grad_theta)
        ),
        "S_bulk.U": dims["U"],
        "S_bulk.quantum": add_dim(
            scale_dim(2, dims["hbar"]), scale_dim(2, d_grad_rho),
            scale_dim(-1, dims["m"]), scale_dim(-1, dims["rho"])
        ),
        "S_chi.kinetic": add_dim(dims["Z_chi"], scale_dim(2, d_r_t)),
        "S_chi.gradient": add_dim(dims["kappa_chi"], scale_dim(2, d_grad_r)),
        "S_chi.potential": dims["lambda_chi"],
        "S_scalar_reduced.u_kinetic": add_dim(dims["A_eff"], scale_dim(2, d_u_t)),
        "S_scalar_reduced.h_kinetic": add_dim(dims["M_h"], scale_dim(2, d_h_t)),
        "S_scalar_reduced.u_gradient": add_dim(dims["B_eff"], scale_dim(2, d_grad_u)),
        "S_scalar_reduced.h_gradient": add_dim(dims["K_h"], scale_dim(2, d_grad_h)),
        "S_scalar_reduced.mix": add_dim(dims["C_hu"], d_grad_u, d_grad_h),
        "S_scalar_parent.H_kinetic": add_dim(dims["M4"], scale_dim(2, d_H_t)),
        "S_scalar_parent.H_gradient": add_dim(dims["K4"], scale_dim(2, d_grad_H)),
        "S_scalar_parent.H_operator": add_dim(dims["K4"], dims["H"], d_O_H),
        "S_mouth.robin": add_dim(dims["eta"], dims["K_m"], scale_dim(2, dims["H"])),
        "S_mouth.source": add_dim(dims["eta"], dims["J_m"], dims["H"]),
        "S_hold": dims["lambda_Sigma"],
        "S_geon_const": dims["E_g"],
        "S_brane.Mac_kinetic": add_dim(dims["g_ell"], dims["rho_br"], scale_dim(2, d_u_t)),
        "S_brane.Mac_shear": add_dim(dims["g_ell"], dims["mu_R"], scale_dim(2, d_omega_u)),
        "S_brane.uw_kinetic": add_dim(dims["g_ell"], dims["rho_br"], scale_dim(2, d_u_t)),
        "S_brane.uw_gap": add_dim(
            dims["g_ell"], dims["rho_br"], scale_dim(2, dims["Omega_w"]),
            scale_dim(2, dims["u"])
        ),
        "S_move": add_dim(dims["q_T"], dims["eta"], dims["V"], dims["u"]),
    }
    measures = {
        key: (
            dt_d4 if key.startswith(("S_bulk", "S_chi", "S_scalar_parent", "S_brane"))
            else dt_surface if key == "S_hold"
            else dt_only if key == "S_geon_const"
            else dt_d3
        )
        for key in integrand_dims
    }
    action_dims = {
        key: add_dim(value, measures[key])
        for key, value in integrand_dims.items()
    }
    bad = {key: value for key, value in action_dims.items() if value != action}
    # Free-carrier-independence is structural: dimensions are integer [L, T, M] tuples.
    actual = {
        "bad_terms": bad,
        "term_count": len(action_dims),
    }
    expected = {
        "bad_terms": {},
        "term_count": 24,
    }
    expect_equal("F_DIMENSIONAL_HOMOGENEITY", actual, expected)
    EVIDENCE["F"] = actual


def continuity_tooth() -> dict[str, Any]:
    x, t = sp.symbols("x t", real=True)
    hbar, mass, drain = sp.symbols("hbar m S_drain", nonzero=True)
    rho = sp.Function("rho")(x, t)
    theta = sp.Function("theta")(x, t)
    source_coupling = (
        hbar * theta * drain if ACTIVE_MUTATION == "G_INJECT_DRAIN" else 0
    )
    lagrangian = (
        -hbar * rho * sp.diff(theta, t)
        -hbar**2 * rho * sp.diff(theta, x)**2 / (2 * mass)
        + source_coupling
    )
    el_theta = sp.simplify(
        sp.diff(lagrangian, theta)
        - sp.diff(sp.diff(lagrangian, sp.diff(theta, t)), t)
        - sp.diff(sp.diff(lagrangian, sp.diff(theta, x)), x)
    )
    v = hbar * sp.diff(theta, x) / mass
    conservative_lhs = hbar * (sp.diff(rho, t) + sp.diff(rho * v, x))
    source_residual = sp.simplify(el_theta - conservative_lhs)
    actual = {
        "EL_minus_source_free_lhs": source_residual,
        "drain_in_action": source_residual != 0,
    }
    expected = {
        "EL_minus_source_free_lhs": sp.Integer(0),
        "drain_in_action": False,
    }
    expect_equal("G_SOURCE_FREE_CONTINUITY", actual, expected)
    EVIDENCE["G"] = actual
    return {"conservative_continuity": "SOURCE_FREE"}


def drain_provenance_tooth() -> dict[str, Any]:
    source_facts = [
        {
            "interface": "part_I_Gamma_B",
            "citation": "stage006 balance manifest",
            "placement": "ADDED_RHS_OUTSIDE_F",
            "variational": False,
            "balance_type": "internal_order_conversion_total_n_conserved",
            "total_carrier_conserved": True,
        },
        {
            "interface": "g0_S_drain",
            "citation": "G0-card sections 2 and 6",
            "placement": "DELIBERATELY_NONVARIATIONAL",
            "variational": False,
            "balance_type": "total_rho_sink_plus_remote_return",
            "total_carrier_conserved": False,
        },
    ]
    if ACTIVE_MUTATION == "H_GAMMA_B_VARIATIONAL":
        source_facts[0]["variational"] = True
    if ACTIVE_MUTATION == "H_G0_DRAIN_VARIATIONAL":
        source_facts[1]["placement"] = "VARIATIONAL_ACTION_TERM"
    expected_facts = [
        {
            "interface": "part_I_Gamma_B",
            "citation": "stage006 balance manifest",
            "placement": "ADDED_RHS_OUTSIDE_F",
            "variational": False,
            "balance_type": "internal_order_conversion_total_n_conserved",
            "total_carrier_conserved": True,
        },
        {
            "interface": "g0_S_drain",
            "citation": "G0-card sections 2 and 6",
            "placement": "DELIBERATELY_NONVARIATIONAL",
            "variational": False,
            "balance_type": "total_rho_sink_plus_remote_return",
            "total_carrier_conserved": False,
        },
    ]
    types = {
        fact["interface"]: fact["balance_type"]
        for fact in source_facts
    }
    both_nonvariational = all(
        not fact["variational"]
        and fact["placement"] in {
            "ADDED_RHS_OUTSIDE_F", "DELIBERATELY_NONVARIATIONAL"
        }
        for fact in source_facts
    )
    mapping_supplied = ACTIVE_MUTATION == "H_SUPPLY_EQUIVALENCE_MAP"
    equivalence = (
        "EQUIVALENT"
        if mapping_supplied
        else "UNRESOLVED"
        if len(set(types.values())) == 2
        else "NOT_ADJUDICATED"
    )
    handoff = {
        "nonvariational_stage": "045",
        "standing_user_gate": True,
    }
    if ACTIVE_MUTATION == "H_SELECT_DRAIN_IN_044":
        handoff["nonvariational_stage"] = "044"
    selection = (
        "DEFERRED_045_AND_USER_GATE"
        if handoff == {
            "nonvariational_stage": "045",
            "standing_user_gate": True,
        }
        else "SELECTED_IN_044"
    )
    actual = {
        "facts": source_facts,
        "interface_token":
            "TWO_NONVARIATIONAL_NAMED" if both_nonvariational else "PROVENANCE_MISMATCH",
        "types": types,
        "mapping_supplied": mapping_supplied,
        "equivalence": equivalence,
        "selection": selection,
    }
    expected = {
        "facts": expected_facts,
        "interface_token": "TWO_NONVARIATIONAL_NAMED",
        "types": {
            "part_I_Gamma_B": "internal_order_conversion_total_n_conserved",
            "g0_S_drain": "total_rho_sink_plus_remote_return",
        },
        "mapping_supplied": False,
        "equivalence": "UNRESOLVED",
        "selection": "DEFERRED_045_AND_USER_GATE",
    }
    expect_equal("H_DRAIN_PROVENANCE", actual, expected)
    EVIDENCE["H"] = actual
    return {
        "drain_interfaces": actual["interface_token"],
        "drain_types": actual["types"],
        "drain_equivalence": actual["equivalence"],
        "drain_selection": actual["selection"],
    }


def p_retirement_tooth() -> dict[str, Any]:
    operative = {"S_GNLS", "gL_Mac", "gL_uw"}
    retired = {"L_pol", "gL_Pu"}
    if ACTIVE_MUTATION == "I_RETIRED_IN_OPERATIVE":
        operative.add("L_pol")
    dof_before, dof_after = 8, 4
    if ACTIVE_MUTATION == "I_DOF_8_TO_5":
        dof_after = 5
    s007_before, s007_after = 11, 7
    if ACTIVE_MUTATION == "I_STAGE007_DRIFT_11_TO_8":
        s007_after = 8
    s006_before, s006_after = 6, 5
    if ACTIVE_MUTATION == "I_STAGE006_DRIFT_6_TO_6":
        s006_after = 6
    stage007_drift = s007_after - s007_before
    stage006_drift = s006_after - s006_before
    actual = {
        "stage007_action_operative": sorted(operative),
        "stage007_action_retired": sorted(retired),
        "intersection": sorted(operative & retired),
        "stage007_DOF": f"{dof_before}->{dof_after}",
        "DOF_removed": dof_before - dof_after,
        "stage007_drift": f"{s007_before}->{s007_after}",
        "stage007_drift_delta": stage007_drift,
        "stage006_drift": f"{s006_before}->{s006_after}",
        "stage006_drift_delta": stage006_drift,
        "net_routeless_delta": stage007_drift + stage006_drift,
    }
    expected = {
        "stage007_action_operative": ["S_GNLS", "gL_Mac", "gL_uw"],
        "stage007_action_retired": ["L_pol", "gL_Pu"],
        "intersection": [],
        "stage007_DOF": "8->4",
        "DOF_removed": 4,
        "stage007_drift": "11->7",
        "stage007_drift_delta": -4,
        "stage006_drift": "6->5",
        "stage006_drift_delta": -1,
        "net_routeless_delta": -5,
    }
    expect_equal("I_P_RETIREMENT", actual, expected)
    EVIDENCE["I"] = actual
    return {
        "P_retirement": {
            "stage007_action_operative": actual["stage007_action_operative"],
            "stage007_action_retired": actual["stage007_action_retired"],
            "stage007_DOF": actual["stage007_DOF"],
            "stage007_drift": actual["stage007_drift"],
            "stage006_drift": actual["stage006_drift"],
            "net_routeless": "-5 == -4(s007 drift) -1(s006 drift)",
        }
    }


def field_union_schur_tooth() -> dict[str, Any]:
    rho_br, c_j, kappa = sp.symbols("rho_br C_J kappa_phase", positive=True)
    pre_schur = sp.Matrix([[rho_br, c_j], [c_j, -kappa]])
    a_eff = sp.simplify(
        pre_schur[0, 0]
        - pre_schur[0, 1] * pre_schur[1, 1]**-1 * pre_schur[1, 0]
    )
    incidence = {
        "theta": {
            "sectors": ["S_bulk"],
            "role": "Madelung_flow_v=grad(theta)",
            "pre_schur_only": False,
        },
        "theta_B": {
            "sectors": ["pre_Schur_brane_phase_block"],
            "role": "brane_phase",
            "pre_schur_only": True,
        },
    }
    if ACTIVE_MUTATION == "J_CONFLATE_PHASE_INCIDENCE":
        incidence["theta_B"] = dict(incidence["theta"])
    theta_distinct = incidence["theta"] != incidence["theta_B"]
    reduced = [
        "rho", "theta", "r_B==chi_B", "u_L", "h",
        "u_T(divfree)", "u_w", "lambda_Sigma(aux)",
    ]
    parent = [
        "rho", "theta", "r_B==chi_B", "H", "u_L",
        "u_T(divfree)", "u_w", "lambda_Sigma(aux)",
    ]
    theta_b_status = "ELIMINATED_SCHUR_into_A_eff"
    if ACTIVE_MUTATION == "J_RETAIN_THETA_B":
        reduced.insert(2, "theta_B")
        parent.insert(2, "theta_B")
        theta_b_status = "INDEPENDENT"
    if ACTIVE_MUTATION == "J_DROP_U_W":
        reduced.remove("u_w")
        parent.remove("u_w")
    physical_fields = {
        "rho", "theta", "r_B==chi_B", "u_L", "h",
        "H", "u_T(divfree)", "u_w",
    }
    if ACTIVE_MUTATION == "J_COUNT_AUX_AS_PHYSICAL_DOF":
        physical_fields.add("lambda_Sigma(aux)")
    lambda_sigma_physical = "lambda_Sigma(aux)" in physical_fields
    transverse_constraint = (
        "none"
        if ACTIVE_MUTATION == "J_DROP_TRANSVERSE_CONSTRAINT"
        else "div(u_T)=0"
    )
    actual = {
        "incidence": incidence,
        "theta_distinct_derived": theta_distinct,
        "A_eff": a_eff,
        "theta_B": theta_b_status,
        "field_set_reduced": reduced,
        "field_set_parent": parent,
        "lambda_Sigma_physical_DOF": lambda_sigma_physical,
        "u_T_constraint": transverse_constraint,
    }
    expected = {
        "incidence": {
            "theta": {
                "sectors": ["S_bulk"],
                "role": "Madelung_flow_v=grad(theta)",
                "pre_schur_only": False,
            },
            "theta_B": {
                "sectors": ["pre_Schur_brane_phase_block"],
                "role": "brane_phase",
                "pre_schur_only": True,
            },
        },
        "theta_distinct_derived": True,
        "A_eff": rho_br + c_j**2 / kappa,
        "theta_B": "ELIMINATED_SCHUR_into_A_eff",
        "field_set_reduced": [
            "rho", "theta", "r_B==chi_B", "u_L", "h",
            "u_T(divfree)", "u_w", "lambda_Sigma(aux)",
        ],
        "field_set_parent": [
            "rho", "theta", "r_B==chi_B", "H", "u_L",
            "u_T(divfree)", "u_w", "lambda_Sigma(aux)",
        ],
        "lambda_Sigma_physical_DOF": False,
        "u_T_constraint": "div(u_T)=0",
    }
    expect_equal("J_FIELD_UNION_SCHUR", actual, expected)
    EVIDENCE["J"] = actual
    return {
        "field_set_reduced": reduced,
        "field_set_parent": parent,
        "theta_B": theta_b_status,
    }


DEFERRED_GATE_NAMES = (
    "CURVED_HELD_WALL_LOWEST_EIGENVALUE",
    "SOURCED_BULK_HESSIAN_LOWEST_EIGENVALUE",
    "DRESSED_H_MONOPOLE",
    "HADAMARD_10",
    "HADAMARD_01",
    "ONE_BODY_FINITE_VOLUME_CLOSURE",
    "DRAIN_ABLATION",
    "PAIR_MASS_MOMENTUM_ENERGY_CLOSURE",
    "OUTER_RETURN_ABLATION",
    "PAIR_FORCE_INTEGRABILITY",
    "PHYSICAL_TOTAL_CHANNEL_READOUTS",
)


def status_teeth() -> dict[str, Any]:
    card_gate_rows = [
        {"gate": gate, "class": 2 if index in {0, 1, 5, 6} else 3,
         "status": "DEFERRED"}
        for index, gate in enumerate(DEFERRED_GATE_NAMES)
    ]
    if ACTIVE_MUTATION == "K_DROP_DEFERRED_GATE":
        card_gate_rows.pop(0)
    completeness_manifest = {
        row["gate"]: row["status"] for row in card_gate_rows
    }
    expect_equal(
        "K1_DEFERRED_GATE_COMPLETENESS",
        completeness_manifest,
        {gate: "DEFERRED" for gate in DEFERRED_GATE_NAMES},
    )

    closed_rows = [dict(row) for row in card_gate_rows]
    if ACTIVE_MUTATION == "K_RESOLVE_ALL_GATES":
        for row in closed_rows:
            row["status"] = "RESOLVED"
    closed_action = not any(row["status"] == "DEFERRED" for row in closed_rows)
    header_facts = {
        "document": "G0 shared minimal closure card v0",
        "status_flag": "DRAFT v0",
        "final_model": False,
    }
    if ACTIVE_MUTATION == "K_CORRUPT_CARD_HEADER":
        header_facts["status_flag"] = "FINAL"
    g0_status = (
        "DRAFT_V0"
        if header_facts["status_flag"] == "DRAFT v0"
        and not header_facts["final_model"]
        else "FINAL"
    )
    k2_actual = {
        "closed_action": closed_action,
        "g0_card": g0_status,
        "remaining_deferred_count": sum(
            row["status"] == "DEFERRED" for row in closed_rows
        ),
    }
    k2_expected = {
        "closed_action": False,
        "g0_card": "DRAFT_V0",
        "remaining_deferred_count": len(DEFERRED_GATE_NAMES),
    }
    expect_equal("K2_DRAFT_OPEN_ACTION", k2_actual, k2_expected)
    EVIDENCE["K"] = {
        "manifest": completeness_manifest,
        "derived": k2_actual,
    }
    return {
        "g0_card": g0_status,
        "closed_action": closed_action,
    }


EXPECTED_INCIDENCE = OrderedDict([
    ("S_bulk", ["rho", "theta(v=grad_theta)"]),
    ("S_chi", ["r_B==chi_B"]),
    ("S_scalar", ["exclusive_scalar_branch", "u_L", "H_or_h"]),
    ("S_mouth", ["Q_chi_fixed_source", "H_at_w0_or_h"]),
    ("S_hold", ["r_B", "lambda_Sigma_aux", "Gamma_Sigma_amb"]),
    ("S_geon_const", ["E_g_i(no_fields)"]),
    ("S_brane", ["g_ell", "u_T(divfree)", "u_w"]),
    ("S_move", ["q_T*s_i*eta_a*V_i_dot_u_T"]),
])


def action_incidence_tooth() -> dict[str, Any]:
    incidence = OrderedDict(
        (key, list(value)) for key, value in EXPECTED_INCIDENCE.items()
    )
    if ACTIVE_MUTATION == "L_DELETE_S_MOVE":
        incidence.pop("S_move")
    selected_branches = {"REDUCED_S_Lh"}
    if ACTIVE_MUTATION == "L_SUM_BOTH_SCALAR_BRANCHES":
        selected_branches.add("PARENT_S_H_plus_S_u_plus_S_mix")
    scalar_xor = len(selected_branches) == 1
    cons = [
        name for name in incidence
        if name in {
            "S_bulk", "S_chi", "S_scalar", "S_mouth",
            "S_hold", "S_geon_const",
        }
    ]
    extra = [
        {
            "S_brane": "S_brane_light",
            "S_move": "S_move_magnetism_moving_coupling",
        }[name]
        for name in incidence if name in {"S_brane", "S_move"}
    ]
    actual = {
        "incidence": incidence,
        "summand_keys": list(incidence),
        "scalar_selected": sorted(selected_branches),
        "scalar_xor": scalar_xor,
        "S_cons_G0_summands": cons,
        "S_assembled_extra": extra,
    }
    expected = {
        "incidence": EXPECTED_INCIDENCE,
        "summand_keys": list(EXPECTED_INCIDENCE),
        "scalar_selected": ["REDUCED_S_Lh"],
        "scalar_xor": True,
        "S_cons_G0_summands": [
            "S_bulk", "S_chi", "S_scalar", "S_mouth",
            "S_hold", "S_geon_const",
        ],
        "S_assembled_extra": [
            "S_brane_light", "S_move_magnetism_moving_coupling",
        ],
    }
    expect_equal("L_ACTION_INCIDENCE", actual, expected)
    EVIDENCE["L"] = actual
    return {
        "S_cons_G0_summands": cons,
        "S_assembled_extra": extra,
        "scalar_branch":
            "REDUCED_S_Lh | PARENT_S_H_plus_S_u_plus_S_mix (exclusive)",
    }


def provenance_status_tooth() -> dict[str, Any]:
    provenance = {
        "g_ell": {
            "edge": "R21",
            "relation": "SCOPE_SPLIT_NOT_REDUCTION",
            "persists": "LIGHT_BRANE_LOCALIZATION_ENVELOPE",
        },
        "light_gating": {
            "edge": "R17",
            "status": "PENDING",
            "relation": "chi_B_projection_distinct_from_g_ell",
        },
        "Z_chi": {
            "source": "G0-card section 2.2",
            "flag": "[POSTULATE]",
            "register_route": None,
            "action_class": "ACTION",
        },
        "gravitomagnetism": {
            "citation": "research/4d_*pn*",
            "coverage": "GR_MATCHED_1PN_TO_4PN_INCLUDING_FRAME_DRAGGING",
            "em_asymmetry": "EM_DEPARTS_EXACT_MAXWELL",
            "stage044_register_row": False,
        },
        "throat_solve": {
            "citation": "ratified-plan/register central shared debt",
            "status": "SIM_DEFERRED_CENTRAL_DEBT",
        },
    }
    base_range = [40, 49]
    shift = [1, 1]
    shifted_range = [left + right for left, right in zip(base_range, shift)]
    count_fact = {
        "base": base_range,
        "shift": shift,
        "sensitivity": shifted_range,
        "spread_before": base_range[1] - base_range[0],
        "spread_after": shifted_range[1] - shifted_range[0],
        "committed_recount": "DEFERRED_046",
        "possible_substitution": "dissipative_M_chi",
    }
    edge_tail = 92
    edge_relations = [
        "wall_static_energy_dedup",
        "field_set_union",
        "drain_named_both",
        "P_retirement_four_manifests",
        "Z_chi_count_reconciliation",
    ]
    edge_sequence = [f"R{edge_tail + offset}" for offset in range(1, 1 + len(edge_relations))]
    master_row = {
        "Param": "Z_chi",
        "[L,T,M]": [-2, 0, 1],
        "Enters": "S_chi kinetic",
        "Class": "ACTION",
        "Depends on / relation": "inertial wall normalization",
        "Reduction route + status":
            "G0 DRAFT-v0 [POSTULATE]; no reduction route",
    }
    if ACTIVE_MUTATION == "M_G_ELL_WRONG_EDGE":
        provenance["g_ell"]["edge"] = "R17"
    if ACTIVE_MUTATION == "M_LIGHT_GATING_WRONG_EDGE":
        provenance["light_gating"]["edge"] = "R21"
    if ACTIVE_MUTATION == "M_Z_CHI_DROP_POSTULATE":
        provenance["Z_chi"]["flag"] = "[DERIVED]"
    if ACTIVE_MUTATION == "M_Z_CHI_INVENT_ROUTE":
        provenance["Z_chi"]["register_route"] = "field_rescaling"
    if ACTIVE_MUTATION == "M_Z_CHI_WRONG_ENDPOINT":
        count_fact["sensitivity"][1] = 49
        count_fact["spread_after"] = (
            count_fact["sensitivity"][1] - count_fact["sensitivity"][0]
        )
    if ACTIVE_MUTATION == "M_GM_DROP_PN_CITATION":
        provenance["gravitomagnetism"]["citation"] = "uncited"
    if ACTIVE_MUTATION == "M_THROAT_WRONG_STATUS":
        provenance["throat_solve"]["status"] = "RESOLVED"
    if ACTIVE_MUTATION == "M_EDGE_SEQUENCE_GAP":
        edge_sequence[2] = "R96"
    if ACTIVE_MUTATION == "M_Z_MASTER_WRONG_DIM":
        master_row["[L,T,M]"] = [0, 0, 0]

    actual = {
        "provenance": provenance,
        "Z_chi_count": count_fact,
        "prospective_edge_sequence": edge_sequence,
        "prospective_edge_tail": edge_sequence[-1],
        "Z_chi_master_parameter_row": master_row,
    }
    expected = {
        "provenance": {
            "g_ell": {
                "edge": "R21",
                "relation": "SCOPE_SPLIT_NOT_REDUCTION",
                "persists": "LIGHT_BRANE_LOCALIZATION_ENVELOPE",
            },
            "light_gating": {
                "edge": "R17",
                "status": "PENDING",
                "relation": "chi_B_projection_distinct_from_g_ell",
            },
            "Z_chi": {
                "source": "G0-card section 2.2",
                "flag": "[POSTULATE]",
                "register_route": None,
                "action_class": "ACTION",
            },
            "gravitomagnetism": {
                "citation": "research/4d_*pn*",
                "coverage": "GR_MATCHED_1PN_TO_4PN_INCLUDING_FRAME_DRAGGING",
                "em_asymmetry": "EM_DEPARTS_EXACT_MAXWELL",
                "stage044_register_row": False,
            },
            "throat_solve": {
                "citation": "ratified-plan/register central shared debt",
                "status": "SIM_DEFERRED_CENTRAL_DEBT",
            },
        },
        "Z_chi_count": {
            "base": [40, 49],
            "shift": [1, 1],
            "sensitivity": [41, 50],
            "spread_before": 9,
            "spread_after": 9,
            "committed_recount": "DEFERRED_046",
            "possible_substitution": "dissipative_M_chi",
        },
        "prospective_edge_sequence": ["R93", "R94", "R95", "R96", "R97"],
        "prospective_edge_tail": "R97",
        "Z_chi_master_parameter_row": {
            "Param": "Z_chi",
            "[L,T,M]": [-2, 0, 1],
            "Enters": "S_chi kinetic",
            "Class": "ACTION",
            "Depends on / relation": "inertial wall normalization",
            "Reduction route + status":
                "G0 DRAFT-v0 [POSTULATE]; no reduction route",
        },
    }
    expect_equal("M_PROVENANCE_STATUS_BINDING", actual, expected)
    EVIDENCE["M"] = actual
    return {
        "g_ell": "SCOPE_SPLIT_R21_NOT_MERGED",
        "light_gating":
            "g_ell_LOCALIZED; chi_B_projection R17 PENDING "
            "(distinct, not double-counted)",
        "gravitomagnetism":
            "BOOST_OF_GRAVITOELECTRIC_FLOW(PN-ladder-covered; honest "
            "asymmetry: gravity matches GR incl GM, EM departs Maxwell); "
            "prose-only in 044",
        "new_draft_knob": {
            "Z_chi":
                "POSTULATE_inertial_normalization_no_reduction_route "
                "(G0 DRAFT-v0)"
        },
        "Z_chi_count_consequence":
            "SHIFTS_CONTINUOUS_[40,49]->[41,50] (043-rule; both endpoints "
            "+1; spread 9 unchanged); committed re-count DEFERRED_046 "
            "(may substitute dissipative M_chi)",
        "throat_solve": provenance["throat_solve"]["status"],
    }


def expected_verdict_object() -> dict[str, Any]:
    return {
        "verdict": "PARENT_ACTION_ASSEMBLED_AT_COMPLETENESS_FLOOR",
        "S_cons_G0_summands": [
            "S_bulk", "S_chi", "S_scalar", "S_mouth",
            "S_hold", "S_geon_const",
        ],
        "S_assembled_extra": [
            "S_brane_light", "S_move_magnetism_moving_coupling",
        ],
        "scalar_branch":
            "REDUCED_S_Lh | PARENT_S_H_plus_S_u_plus_S_mix (exclusive)",
        "field_set_reduced": [
            "rho", "theta", "r_B==chi_B", "u_L", "h",
            "u_T(divfree)", "u_w", "lambda_Sigma(aux)",
        ],
        "field_set_parent": [
            "rho", "theta", "r_B==chi_B", "H", "u_L",
            "u_T(divfree)", "u_w", "lambda_Sigma(aux)",
        ],
        "theta_B": "ELIMINATED_SCHUR_into_A_eff",
        "wall_dedup":
            "STATIC_WALL_ENERGY_IDENTITY(r_B==chi_B; "
            "kappa_chi==kappa_B; lambda_chi==4*a_B => "
            "ell==delta, sigma matches); "
            "FULL_FUNCTIONALS_DIFFER(inertial Z_chi vs dissipative M_chi)",
        "wall_factorization":
            "L_chi_2ndVar/kappa_chi == A_chi_dag A_chi "
            "(annihilates r0')",
        "g_ell": "SCOPE_SPLIT_R21_NOT_MERGED",
        "light_gating":
            "g_ell_LOCALIZED; chi_B_projection R17 PENDING "
            "(distinct, not double-counted)",
        "drain_interfaces": "TWO_NONVARIATIONAL_NAMED",
        "drain_types": {
            "part_I_Gamma_B": "internal_order_conversion_total_n_conserved",
            "g0_S_drain": "total_rho_sink_plus_remote_return",
        },
        "drain_equivalence": "UNRESOLVED",
        "drain_selection": "DEFERRED_045_AND_USER_GATE",
        "conservative_continuity": "SOURCE_FREE",
        "wave_speeds": {
            "c_s0_sq": "5*K*rho0^4/m",
            "c_pm_sq": [0.7928932188, 2.2071067812],
            "c_gamma_sq": "mu_R/rho_br",
            "D_star": 1.75,
        },
        "stability":
            "POSITIVE_DEFINITE(B_eff>0 AND D_star=7/4>0; "
            "wall minima {0,1}; Omega_w^2>=0)",
        "P_retirement": {
            "stage007_action_operative": ["S_GNLS", "gL_Mac", "gL_uw"],
            "stage007_action_retired": ["L_pol", "gL_Pu"],
            "stage007_DOF": "8->4",
            "stage007_drift": "11->7",
            "stage006_drift": "6->5",
            "net_routeless": "-5 == -4(s007 drift) -1(s006 drift)",
        },
        "gravitomagnetism":
            "BOOST_OF_GRAVITOELECTRIC_FLOW(PN-ladder-covered; honest "
            "asymmetry: gravity matches GR incl GM, EM departs Maxwell); "
            "prose-only in 044",
        "new_draft_knob": {
            "Z_chi":
                "POSTULATE_inertial_normalization_no_reduction_route "
                "(G0 DRAFT-v0)"
        },
        "Z_chi_count_consequence":
            "SHIFTS_CONTINUOUS_[40,49]->[41,50] (043-rule; both endpoints "
            "+1; spread 9 unchanged); committed re-count DEFERRED_046 "
            "(may substitute dissipative M_chi)",
        "g0_card": "DRAFT_V0",
        "closed_action": False,
        "throat_solve": "SIM_DEFERRED_CENTRAL_DEBT",
    }


def reproduction_tooth(
    pieces: Iterable[dict[str, Any]],
    field_producers: dict[str, str],
) -> dict[str, Any]:
    verdict = {}
    for piece in pieces:
        verdict.update(piece)
    first_match_results = OrderedDict(CHECK_RESULTS)
    if ACTIVE_MUTATION == "REPRODUCTION_FIRST_MATCH":
        first_match_results["A_WALL_STATIC_ENERGY_DEDUP"] = False
    verdict["verdict"] = first_match(first_match_results)

    # Reorder only for human readability; canonical serialization sorts recursively.
    display_order = [
        "verdict", "S_cons_G0_summands", "S_assembled_extra", "scalar_branch",
        "field_set_reduced", "field_set_parent", "theta_B", "wall_dedup",
        "wall_factorization", "g_ell", "light_gating", "drain_interfaces",
        "drain_types", "drain_equivalence", "drain_selection",
        "conservative_continuity", "wave_speeds", "stability", "P_retirement",
        "gravitomagnetism", "new_draft_knob", "Z_chi_count_consequence",
        "g0_card", "closed_action", "throat_solve",
    ]
    verdict = {key: verdict[key] for key in display_order}
    binding = dict(field_producers)
    if ACTIVE_MUTATION == "REPRODUCTION_FIELD_BINDING":
        binding.pop("throat_solve")
    actual = {
        "verdict_object": verdict,
        "producer_fields": sorted(binding),
        "object_fields": sorted(verdict),
        "all_prior_teeth_pass": all(first_match_results.values()),
    }
    expected_object = expected_verdict_object()
    expected = {
        "verdict_object": expected_object,
        "producer_fields": sorted(expected_object),
        "object_fields": sorted(expected_object),
        "all_prior_teeth_pass": True,
    }
    expect_equal("REPRODUCTION", actual, expected)
    EVIDENCE["REPRODUCTION"] = {
        "first_match_order": list(TOOTH_ORDER[:-1]),
        "field_producers": binding,
    }
    return verdict


def run_audit() -> dict[str, Any]:
    pieces: list[dict[str, Any]] = []
    producers: dict[str, str] = {}

    def take(piece: dict[str, Any], tooth: str) -> None:
        pieces.append(piece)
        for key in piece:
            if key in producers:
                raise RuntimeError(f"duplicate verdict producer for {key}")
            producers[key] = tooth

    section("A-B. Static wall reconciliation and wall Hessian")
    take(wall_static_energy_tooth(), "A")
    take(wall_hessian_tooth(), "B")

    section("C-E. Localized scalar, wave speeds, and stability")
    charge_factorization_teeth()
    take(wave_speed_teeth(), "D")
    take(stability_teeth(), "E")

    section("F-H. Units, conservative continuity, and drain provenance")
    dimensional_tooth()
    take(continuity_tooth(), "G")
    take(drain_provenance_tooth(), "H")

    section("I-L. Retirement, field union, open status, and action incidence")
    take(p_retirement_tooth(), "I")
    take(field_union_schur_tooth(), "J")
    take(status_teeth(), "K")
    take(action_incidence_tooth(), "L")

    section("M. Provenance/status binding and prospective register handoff")
    take(provenance_status_tooth(), "M")

    # C's N0 is a checked intermediate, not a verdict field.  The directive maps
    # wave_speeds/stability to C/D/E collectively; D/E produce the JSON fields.
    expected_producer = {
        "verdict": "REPRODUCTION",
        "S_cons_G0_summands": "L",
        "S_assembled_extra": "L",
        "scalar_branch": "L",
        "field_set_reduced": "J",
        "field_set_parent": "J",
        "theta_B": "J",
        "wall_dedup": "A",
        "wall_factorization": "B",
        "g_ell": "M",
        "light_gating": "M",
        "drain_interfaces": "H",
        "drain_types": "H",
        "drain_equivalence": "H",
        "drain_selection": "H",
        "conservative_continuity": "G",
        "wave_speeds": "D",
        "stability": "E",
        "P_retirement": "I",
        "gravitomagnetism": "M",
        "new_draft_knob": "M",
        "Z_chi_count_consequence": "M",
        "g0_card": "K",
        "closed_action": "K",
        "throat_solve": "M",
    }
    producers["verdict"] = "REPRODUCTION"
    if producers != expected_producer:
        raise AuditFailure("PRODUCER_MAP_INTERNAL", repr(producers))

    section("Reproduction. First-match verdict and canonical object")
    return reproduction_tooth(pieces, producers)


def main() -> None:
    if len(sys.argv) == 2 and sys.argv[1] == "--mutation-map":
        print(canonical_json(MUTATION_TO_TOOTH), end="")
        return
    if len(sys.argv) != 1:
        raise AuditFailure("CLI", "only --mutation-map is accepted")
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in MUTATION_TO_TOOTH:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage044_parent_action_reconciliation SymPy audit")
    print(
        "ROUTE=SymPy exact differential algebra + symbolic EL variation + "
        "dimension-vector assembly + structured manifests"
    )
    print("VERDICT_DERIVATION=DOCUMENTED_FIRST_MATCH")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"EXPECTED_TOOTH={MUTATION_TO_TOOTH[ACTIVE_MUTATION]}")

    verdict = run_audit()
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    root = Path(__file__).resolve().parents[3]
    verdict_path = root / "research/pde_ledger_v2/_scratch/stage044/verdict_py.json"
    verdict_path.parent.mkdir(parents=True, exist_ok=True)
    verdict_path.write_text(canonical_json(verdict), encoding="utf-8")
    print("")
    print("CANONICAL_JSON=" + canonical_json(verdict).rstrip("\n"))
    print(f"VERDICT_FILE={verdict_path}")
    print(f"TOOTH_COUNT={len(TOOTH_ORDER)}")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy stage044 parent-action reconciliation")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(
            "OVERALL FAIL: SymPy stage044 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
