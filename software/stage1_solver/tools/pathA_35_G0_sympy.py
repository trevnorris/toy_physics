#!/usr/bin/env python3
"""SymPy G0 freeze checker for pathA_35.

Scope: G0 only.  This script checks freeze-byte fidelity, restored-unit
homogeneity for every expression it constructs, able-to-fail ablations, and the
flat-brane linear DOF count.  It does not compute any Gate L verdict.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
SCRATCH = STAGE1_ROOT / "_scratch"
T0_REPORT = STAGE1_ROOT / "reports" / "pathA_24_T0_freeze.md"
G0_REPORT = STAGE1_ROOT / "reports" / "pathA_35_G0_freeze.md"

EXPECTED_T0_HASH = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064"
EXPECTED_G0_HASH = "d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3"
EXPECTED_G0_SHORT = "d9520d3819c3"

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
        if power == 1:
            parts.append(label)
        else:
            parts.append(f"{label}^{power}")
    return " ".join(parts) if parts else "1"


def dim_record(dim: Dim) -> dict[str, Any]:
    return {"triple_MLT": list(dim), "string": dim_str(dim)}


def extract_fence_bytes(path: Path, label: str) -> bytes:
    start = f"```{label}\n".encode("utf-8")
    end = b"```\n"
    lines = path.read_bytes().splitlines(keepends=True)
    blocks: list[bytes] = []
    in_block = False
    current: list[bytes] = []
    for line in lines:
        if not in_block and line == start:
            in_block = True
            current = []
            continue
        if in_block and line == end:
            blocks.append(b"".join(current))
            in_block = False
            current = []
            continue
        if in_block:
            current.append(line)
    if in_block:
        raise AssertionError(f"unterminated {label!r} fence in {path}")
    if len(blocks) != 1:
        raise AssertionError(f"expected exactly one {label!r} fence in {path}, found {len(blocks)}")
    return blocks[0]


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def byte_range(path: Path, label: str) -> dict[str, int | list[int]]:
    marker = f"```{label}\n".encode("utf-8")
    data = path.read_bytes()
    marker_at = data.index(marker)
    start = marker_at + len(marker)
    end = data.index(b"\n```", start)
    return {
        "start_0_based_inclusive": start,
        "end_0_based_exclusive": end,
        "range_0_based": [start, end],
        "start_1_based_inclusive": start + 1,
        "end_1_based_inclusive": end,
        "length_bytes": end - start,
    }


def check_freeze_fidelity() -> dict[str, Any]:
    t0_block = extract_fence_bytes(T0_REPORT, "freeze-action")
    g0_block = extract_fence_bytes(G0_REPORT, "freeze-action")
    t0_hash = sha256_hex(t0_block)
    g0_hash = sha256_hex(g0_block)
    if t0_hash != EXPECTED_T0_HASH:
        raise AssertionError(f"T0 hash mismatch: expected {EXPECTED_T0_HASH}, got {t0_hash}")
    if g0_hash != EXPECTED_G0_HASH:
        raise AssertionError(f"G0 hash mismatch: expected {EXPECTED_G0_HASH}, got {g0_hash}")
    if t0_block not in g0_block:
        raise AssertionError("byte-identical T0 freeze-action block is not embedded in G0 block")
    return {
        "t0_hash": t0_hash,
        "g0_hash": g0_hash,
        "g0_short_hash": EXPECTED_G0_SHORT,
        "t0_bytes_embedded_in_g0": True,
        "g0_byte_range": byte_range(G0_REPORT, "freeze-action"),
    }


class DimChecker:
    def __init__(self) -> None:
        self.records: list[dict[str, Any]] = []
        self.ablations: list[dict[str, Any]] = []

    def check(self, category: str, name: str, actual: Dim, expected: Dim, expression: str) -> Dim:
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
        return actual

    def check_same(
        self, category: str, name: str, dims: list[Dim], expected: Dim, expression: str
    ) -> Dim:
        for index, dim in enumerate(dims):
            self.check(category, f"{name}.part_{index}", dim, expected, expression)
        self.records.append(
            {
                "category": category,
                "name": name,
                "expression": expression,
                "dimension": dim_record(expected),
                "expected": dim_record(expected),
                "status": "PASS",
                "homogeneous_parts": len(dims),
            }
        )
        return expected

    def expect_fail(self, category: str, name: str, actual: Dim, expected: Dim, expression: str) -> None:
        fired = actual != expected
        if not fired:
            raise AssertionError(f"ablation did not fire: {category}:{name}")
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


def build_dimension_payload() -> dict[str, Any]:
    check = DimChecker()

    M: Dim = (1, 0, 0)
    L: Dim = (0, 1, 0)
    T: Dim = (0, 0, 1)
    Z: Dim = (0, 0, 0)

    bulk_lag = (1, -2, -2)
    brane_lag = (1, -1, -2)
    action_dim = (1, 2, -1)
    stress = bulk_lag
    eom_u_op = (1, -3, -2)

    d_m = M
    d_hbar = (1, 2, -1)
    d_rho = (0, -4, 0)
    d_K = (1, 18, -2)
    d_a = L
    d_u = L
    d_P = Z
    d_w = L
    d_ell_g = L
    d_g = (0, -1, 0)
    d_grad = (0, -1, 0)
    d_dt = (0, 0, -1)
    d_dt_measure = (0, 0, 1)
    d_d4x = (0, 4, 0)
    d_v = (0, 1, -1)
    d_k = (0, -1, 0)
    d_omega = (0, 0, -1)
    d_rho_br = (1, -3, 0)
    d_mu_R = brane_lag
    d_lambda_Pu = brane_lag
    d_Omega_w = d_omega

    d_cs2 = dsub(dadd(d_K, dmul(4, d_rho)), d_m)
    check.check("kept_gnls", "c_s^2(rho)", d_cs2, (0, 2, -2), "5 K rho^4 / m")
    check.check("kept_gnls", "U(rho)", dadd(d_K, dmul(5, d_rho)), bulk_lag, "(K/4) rho^5")
    check.check(
        "kept_gnls",
        "quantum_pressure",
        dadd(dmul(2, d_hbar), dmul(-1, d_m), dmul(-1, d_rho), dmul(2, dadd(d_grad, d_rho))),
        bulk_lag,
        "(hbar^2/(8 m rho)) (partial_i rho)^2",
    )
    check.check(
        "kept_gnls",
        "bulk_kinetic_stress_scale",
        dadd(d_m, d_rho, dmul(2, d_v)),
        stress,
        "m rho v_i v_j",
    )

    d_DtP = d_dt
    d_gradP = dadd(d_grad, d_P)
    check.check(
        "kept_T0",
        "P_kinetic",
        dadd(d_m, d_rho, dmul(2, d_a), dmul(2, d_DtP)),
        bulk_lag,
        "m rho a^2 (D_t^v P)^2",
    )
    check.check(
        "kept_T0",
        "P_Frank",
        dadd(d_m, d_rho, d_cs2, dmul(2, d_a), dmul(2, d_gradP)),
        bulk_lag,
        "m rho c_s^2 a^2 (partial_j P^i)^2",
    )
    check.check(
        "kept_T0",
        "P_radial_potential",
        dadd(d_m, d_rho, d_cs2),
        bulk_lag,
        "m rho c_s^2 (P^i P^i - 1)^2",
    )
    check.check("kept_T0", "bulk_couple_inertia", dadd(d_m, d_rho, dmul(2, d_a)), (1, -2, 0), "m rho a^2")
    check.check(
        "kept_T0",
        "bulk_couple_stiffness",
        dadd(d_m, d_rho, d_cs2, dmul(2, d_a)),
        (1, 0, -2),
        "m rho c_s^2 a^2",
    )
    check.check("kept_T0", "bulk_radial_scale", dadd(d_m, d_rho, d_cs2), bulk_lag, "m rho c_s^2")

    check.check("confinement", "w_over_ell_g", dsub(d_w, d_ell_g), Z, "w/ell_g")
    check.check("confinement", "g_ell", d_g, (0, -1, 0), "exp(-(w/ell_g)^2)/(sqrt(pi) ell_g)")
    check.check("confinement", "dw_g_ell", dadd(L, d_g), Z, "dw g_ell(w)")

    d_dtu = dadd(d_dt, d_u)
    d_curlu = dadd(d_grad, d_u)
    d_uw = d_u
    d_dtuw = dadd(d_dt, d_uw)
    d_varpi = d_P
    check.check("brane_surface", "rho_br", d_rho_br, (1, -3, 0), "rho_br")
    mac_kin = check.check(
        "brane_surface",
        "u_kinetic",
        dadd(d_rho_br, dmul(2, d_dtu)),
        brane_lag,
        "rho_br (partial_t u)^2",
    )
    mac_curl = check.check(
        "brane_surface",
        "MacCullagh_curl",
        dadd(d_mu_R, dmul(2, d_curlu)),
        brane_lag,
        "mu_R (nabla_parallel x u)^2",
    )
    check.check_same("brane_surface", "L_Mac_homogeneous", [mac_kin, mac_curl], brane_lag, "L_Mac")
    pu_term = check.check(
        "brane_surface",
        "P_u_coupling",
        dadd(d_lambda_Pu, d_varpi, d_curlu),
        brane_lag,
        "lambda_Pu varpi_a Omega_u^a",
    )
    uw_kin = check.check(
        "brane_surface",
        "u_w_kinetic",
        dadd(d_rho_br, dmul(2, d_dtuw)),
        brane_lag,
        "rho_br (partial_t u_w)^2",
    )
    uw_gap = check.check(
        "brane_surface",
        "u_w_gap",
        dadd(d_rho_br, dmul(2, d_Omega_w), dmul(2, d_uw)),
        brane_lag,
        "rho_br Omega_w^2 u_w^2",
    )
    check.check_same("brane_surface", "L_uw_homogeneous", [uw_kin, uw_gap], brane_lag, "L_uw")
    for name, dim in {
        "g_L_Mac_kinetic": mac_kin,
        "g_L_Mac_curl": mac_curl,
        "g_L_Pu": pu_term,
        "g_L_uw_kinetic": uw_kin,
        "g_L_uw_gap": uw_gap,
    }.items():
        check.check("brane_bulk_representation", name, dadd(d_g, dim), bulk_lag, f"g_ell(w) {name}")

    check.check("action_measure", "bulk_action_integral", dadd(d_dt_measure, d_d4x, bulk_lag), action_dim, "int dt d^4X bulk_lag")
    check.check(
        "action_measure",
        "brane_bulk_action_integral",
        dadd(d_dt_measure, d_d4x, d_g, brane_lag),
        action_dim,
        "int dt d^4X g_ell(w) L_brane",
    )

    t_wa = check.check("traction", "T_wa", dadd(d_m, d_rho, d_v, d_v), stress, "m rho v_w v_a")
    slope = check.check("traction", "partial_b_u_w", dadd(d_grad, d_uw), Z, "partial_b u_w")
    stress_slope = check.check(
        "traction",
        "slope_mixing",
        dadd(stress, slope),
        stress,
        "(T_ww delta_ab - T_ab) partial_b u_w",
    )
    check.check_same(
        "traction",
        "T_na_projected",
        [t_wa, stress_slope],
        stress,
        "T_wa + (T_ww delta_ab - T_ab) partial_b u_w",
    )

    check.expect_fail(
        "ablation",
        "drop_m_from_T_wa",
        dadd(d_rho, d_v, d_v),
        stress,
        "rho v_w v_a",
    )
    check.expect_fail(
        "ablation",
        "MacCullagh_without_curl",
        dadd(d_mu_R, dmul(2, d_u)),
        brane_lag,
        "mu_R u^2",
    )

    op_rho = check.check("linearization", "rho_br_omega2", dadd(d_rho_br, dmul(2, d_omega)), eom_u_op, "rho_br omega^2")
    op_mu = check.check("linearization", "mu_R_k2", dadd(d_mu_R, dmul(2, d_k)), eom_u_op, "mu_R k^2")
    check.check_same("linearization", "O_u_homogeneous", [op_rho, op_mu], eom_u_op, "rho_br omega^2 I - mu_R(k^2 I-k k^T)")
    check.check("linearization", "c_gamma_squared", dsub(d_mu_R, d_rho_br), (0, 2, -2), "mu_R/rho_br")
    check.check("linearization", "omega_T_squared", dadd(dsub(d_mu_R, d_rho_br), dmul(2, d_k)), (0, 0, -2), "c_gamma^2 k^2")
    check.check("linearization", "P_spin_omega_squared", dadd(d_cs2, dmul(2, d_k)), (0, 0, -2), "c_s^2 k^2")
    check.check("linearization", "P_radial_gap_term", dsub(d_cs2, dmul(2, d_a)), (0, 0, -2), "2 c_s^2/a^2")
    check.check_same(
        "linearization",
        "P_radial_omega_squared",
        [dadd(d_cs2, dmul(2, d_k)), dsub(d_cs2, dmul(2, d_a))],
        (0, 0, -2),
        "c_s^2 k^2 + 2 c_s^2/a^2",
    )
    check.check("linearization", "u_w_bare_omega_squared", dmul(2, d_Omega_w), (0, 0, -2), "Omega_w^2")
    check.check(
        "linearization",
        "Fourier_P_u_monomial",
        dadd(d_lambda_Pu, d_P, d_k, d_u),
        brane_lag,
        "lambda_Pu P k u",
    )

    parity = {
        "direct_P_dot_curlu": {
            "operator": "P_parallel . Omega_u",
            "parity": "odd_pseudoscalar",
            "status": "excluded",
        },
        "repaired_w_dot_P_cross_curlu": {
            "operator": "w_hat . (P_parallel x Omega_u)",
            "parity": "even_scalar_using_imposed_normal",
            "status": "frozen",
        },
    }

    constants = {
        "rho": dim_record(d_rho),
        "m": dim_record(d_m),
        "K": dim_record(d_K),
        "hbar": dim_record(d_hbar),
        "a": dim_record(d_a),
        "rho_br": dim_record(d_rho_br),
        "mu_R": dim_record(d_mu_R),
        "lambda_Pu": dim_record(d_lambda_Pu),
        "Omega_w": dim_record(d_Omega_w),
        "ell_g": dim_record(d_ell_g),
        "g_ell": dim_record(d_g),
    }

    return {
        "base_dimensions": {"M": [1, 0, 0], "L": [0, 1, 0], "T": [0, 0, 1]},
        "targets": {
            "bulk_action_density": dim_record(bulk_lag),
            "brane_action_density": dim_record(brane_lag),
            "action": dim_record(action_dim),
            "stress": dim_record(stress),
            "u_eom_operator": dim_record(eom_u_op),
        },
        "constants": constants,
        "checks": check.records,
        "ablations": check.ablations,
        "parity": parity,
    }


def matrix_rank(matrix: sp.Matrix) -> int:
    return int(matrix.rank())


def compute_flat_brane_dof(
    *,
    u_active: tuple[int, int, int] = (1, 1, 1),
    p_kinetic_present: bool = True,
    p_frank_present: bool = True,
    p_radial_present: bool = True,
    u_w_kinetic_present: bool = True,
    u_w_gap_present: bool = True,
    phi_present: bool = False,
) -> dict[str, Any]:
    """Count G0 flat-brane field-space DOF from frozen quadratic terms only."""
    eye_u = sp.eye(3)
    active_u = sp.diag(*u_active)
    curl_projector = sp.diag(0, 1, 1)
    u_kinetic_form = active_u * eye_u * active_u
    u_curl_form = active_u * curl_projector * active_u
    u_kinetic_rank = matrix_rank(u_kinetic_form)
    u_curl_rank = matrix_rank(u_curl_form)
    u_curl_nullity = u_kinetic_rank - u_curl_rank
    if u_curl_nullity < 0:
        raise AssertionError("u curl rank exceeds active u kinetic rank")

    eye_p = sp.eye(4)
    zero_p = sp.zeros(4, 4)
    tangent_p = sp.diag(1, 1, 1, 0)
    radial_p = sp.diag(0, 0, 0, 1)
    p_kinetic_form = eye_p if p_kinetic_present else zero_p
    p_frank_form = eye_p if p_frank_present else zero_p
    p_radial_hessian = radial_p if p_radial_present else zero_p
    p_tangent_kinetic_rank = matrix_rank(tangent_p * p_kinetic_form * tangent_p)
    p_tangent_frank_rank = matrix_rank(tangent_p * p_frank_form * tangent_p)
    p_radial_kinetic_rank = matrix_rank(radial_p * p_kinetic_form * radial_p)
    p_radial_hessian_rank = matrix_rank(radial_p * p_radial_hessian * radial_p)

    u_w_kinetic_rank = matrix_rank(sp.Matrix([[1 if u_w_kinetic_present else 0]]))
    u_w_gap_rank = matrix_rank(sp.Matrix([[1 if u_w_gap_present else 0]]))
    phi_kinetic_rank = matrix_rank(sp.Matrix([[1 if phi_present else 0]]))

    computed_counts = {
        "u_transverse": u_curl_rank,
        "u_longitudinal": u_curl_nullity,
        "P_spin_wave": min(p_tangent_kinetic_rank, p_tangent_frank_rank),
        "P_soft_spin_radial": min(p_radial_kinetic_rank, p_radial_hessian_rank),
        "u_w": min(u_w_kinetic_rank, u_w_gap_rank),
        "phi": phi_kinetic_rank,
    }
    computed_total = sum(computed_counts.values())
    return {
        "counts": computed_counts,
        "total": computed_total,
        "rank_bookkeeping": {
            "u_kinetic_rank": u_kinetic_rank,
            "u_curl_rank": u_curl_rank,
            "u_curl_nullity_within_active_kinetic_space": u_curl_nullity,
            "P_tangent_kinetic_rank": p_tangent_kinetic_rank,
            "P_tangent_Frank_rank": p_tangent_frank_rank,
            "P_radial_kinetic_rank": p_radial_kinetic_rank,
            "P_radial_soft_spin_hessian_rank": p_radial_hessian_rank,
            "u_w_kinetic_rank": u_w_kinetic_rank,
            "u_w_gap_rank": u_w_gap_rank,
            "phi_kinetic_rank": phi_kinetic_rank,
        },
    }


def dof_count_ablation(name: str, mutation: str, baseline: dict[str, Any], mutated: dict[str, Any]) -> dict[str, Any]:
    if mutated["total"] == baseline["total"]:
        raise AssertionError(f"DOF-count ablation did not fire: {name}")
    return {
        "name": name,
        "mutation": mutation,
        "baseline_total": baseline["total"],
        "mutated_total": mutated["total"],
        "mutated_counts": mutated["counts"],
        "status": "FIRED",
    }


def flat_brane_modes() -> dict[str, Any]:
    k = sp.symbols("k", nonzero=True)
    curl_stiffness = sp.Matrix([[0, 0, 0], [0, k**2, 0], [0, 0, k**2]])
    rank = int(curl_stiffness.rank())
    nullity = 3 - rank
    if rank != 2 or nullity != 1:
        raise AssertionError(f"unexpected curl rank/nullity: rank={rank}, nullity={nullity}")

    computed = compute_flat_brane_dof()
    reported_counts = dict(computed["counts"])
    reported_total = sum(reported_counts.values())
    if reported_counts != computed["counts"] or reported_total != computed["total"]:
        raise AssertionError("reported flat-brane DOF count is not wired to computed count")
    if reported_counts["u_transverse"] != rank or reported_counts["u_longitudinal"] != nullity:
        raise AssertionError("reported u-sector count is not wired to curl rank/nullity")
    if reported_total != 8:
        raise AssertionError(f"unexpected flat-brane DOF count: {reported_total}")

    dof_ablations = [
        dof_count_ablation(
            "drop_u_w_gap_term",
            "set the frozen u_w gap matrix Omega_w^2 u_w^2 to rank 0",
            computed,
            compute_flat_brane_dof(u_w_gap_present=False),
        ),
        dof_count_ablation(
            "drop_P_soft_spin_radial_term",
            "set the T0 (P^i P^i - 1)^2 radial Hessian to rank 0",
            computed,
            compute_flat_brane_dof(p_radial_present=False),
        ),
        dof_count_ablation(
            "zero_u_longitudinal_component",
            "remove the k-parallel u component from the kinetic field space",
            computed,
            compute_flat_brane_dof(u_active=(0, 1, 1)),
        ),
    ]

    return {
        "computed_from_quadratic_terms": True,
        "reported_counts_wired_to_computation": True,
        "basis_order": {
            "u": ["u_x_parallel_to_k", "u_y_transverse", "u_z_transverse"],
            "P": ["P_tangent_1", "P_tangent_2", "P_tangent_3", "P_radial"],
            "u_w": ["u_w"],
            "phi": [],
        },
        "quadratic_terms_used": {
            "u": "rho_br (partial_t u^a)^2 and mu_R (nabla_parallel x u)^2",
            "P": "m rho a^2 (D_t P^i)^2, m rho c_s^2 a^2 (partial_j P^i)^2, and m rho c_s^2 (P^i P^i - 1)^2",
            "u_w": "rho_br (partial_t u_w)^2 and rho_br Omega_w^2 u_w^2",
            "phi": "absent in active baseline",
        },
        "curl_stiffness_matrix_for_k_along_x": [[sp.sstr(x) for x in row] for row in curl_stiffness.tolist()],
        "curl_rank": rank,
        "curl_nullity": nullity,
        "rank_bookkeeping": computed["rank_bookkeeping"],
        "counts": reported_counts,
        "total": reported_total,
        "dof_count_ablations": dof_ablations,
        "bulk_shear_status": "kept_GNLS_bulk_shear_free",
        "classification_guard": "counts_only_no_gate_verdict_no_boundedness_no_gauge_no_leak_claim",
    }


def ledger() -> dict[str, Any]:
    return {
        "field_dof_subcount_new_at_G0": 4,
        "constant_subcount": 4,
        "function_subcount": 1,
        "structural_subcount": 6,
        "independent_new_input_count_n": 11,
        "drift_verdict": "SECOND_MEDIUM_DRIFT_AT_FREEZE(11)",
        "flat_brane_total_listed_dof": 8,
    }


def build_report() -> dict[str, Any]:
    freeze = check_freeze_fidelity()
    dimensions = build_dimension_payload()
    modes = flat_brane_modes()
    machine_ledger = ledger()
    agreement_payload = {
        "freeze": freeze,
        "dimension_checks": dimensions,
        "flat_brane_modes": modes,
        "ledger": machine_ledger,
    }
    return {
        "schema": "pathA_35_G0_sympy/v1",
        "engine": "sympy",
        "pass": True,
        "scope": "G0 freeze only; no Gate L verdict computed",
        "agreement_payload": agreement_payload,
        "verdicts": [
            f"T0_SHEAR_FROZEN({EXPECTED_G0_SHORT})",
            "SECOND_MEDIUM_DRIFT_AT_FREEZE(11)",
        ],
    }


def compare_payloads() -> dict[str, Any]:
    sympy_path = SCRATCH / "pathA_35_G0_sympy.json"
    math_path = SCRATCH / "pathA_35_G0_mathematica.json"
    if not sympy_path.exists():
        raise AssertionError(f"missing SymPy output: {sympy_path}")
    if not math_path.exists():
        raise AssertionError(f"missing Mathematica output: {math_path}")
    sympy_payload = json.loads(sympy_path.read_text(encoding="utf-8"))["agreement_payload"]
    math_payload = json.loads(math_path.read_text(encoding="utf-8"))["agreement_payload"]
    if sympy_payload != math_payload:
        raise AssertionError("ENGINE_DISAGREE: SymPy and Mathematica agreement_payload differ")
    result = {
        "schema": "pathA_35_G0_engine_agreement/v1",
        "pass": True,
        "verdict": "ENGINE_AGREE",
        "compared_files": [str(sympy_path), str(math_path)],
    }
    out = SCRATCH / "pathA_35_G0_engine_agreement.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true", help="compare SymPy and Mathematica scratch payloads")
    args = parser.parse_args()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    if args.compare:
        result = compare_payloads()
        print(f"wrote {SCRATCH / 'pathA_35_G0_engine_agreement.json'}")
        print(result["verdict"])
        return
    report = build_report()
    out = SCRATCH / "pathA_35_G0_sympy.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    print("pathA_35 G0 SymPy: PASS")


if __name__ == "__main__":
    main()
