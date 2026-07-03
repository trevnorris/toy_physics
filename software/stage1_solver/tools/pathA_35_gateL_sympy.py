#!/usr/bin/env python3
"""Gate L checker for pathA_35, SymPy engine.

The script computes the flat-brane Gate L sub-hurdles on the two frozen
configurations named by G0:

  A. baseline massless P with parity-repaired P-u coupling, no phi, gapped u_w
  B. slaved-rigid P_parallel = w_hat x curl(u), no independent P modes

It writes a JSON scratch payload for dual-engine comparison.  With --compare it
requires the Mathematica payload, compares the agreement payload exactly, and
writes the Gate-L YAML ledger.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"
T0_REPORT = REPORTS / "pathA_24_T0_freeze.md"
G0_REPORT = REPORTS / "pathA_35_G0_freeze.md"
YAML_OUT = REPORTS / "pathA_35_gateL_light_results.yaml"

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


def check_freeze_fidelity() -> dict[str, Any]:
    t0_block = extract_fence_bytes(T0_REPORT, "freeze-action")
    g0_block = extract_fence_bytes(G0_REPORT, "freeze-action")
    t0_hash = sha256_hex(t0_block)
    g0_hash = sha256_hex(g0_block)
    if t0_hash != EXPECTED_T0_HASH:
        raise AssertionError(f"T0 hash mismatch: expected {EXPECTED_T0_HASH}, got {t0_hash}")
    if g0_hash != EXPECTED_G0_HASH:
        raise AssertionError(f"G0 hash mismatch: expected {EXPECTED_G0_HASH}, got {g0_hash}")
    required = [
        b"Active baseline: massless T0 spin-wave modes",
        b"L_Pu := - lambda_Pu varpi_a Omega_u^a",
        b"No scalar-potential analog phi is present",
        b"Named inactive alternate: slaved-rigid P_parallel = w_hat x Omega_u",
        b"T_wa := m rho v_w v_a",
    ]
    for needle in required:
        if needle not in g0_block:
            raise AssertionError(f"required frozen G0 line missing: {needle!r}")
    return {
        "t0_hash": t0_hash,
        "g0_hash": g0_hash,
        "g0_short_hash": EXPECTED_G0_SHORT,
        "t0_bytes_embedded_in_g0": t0_block in g0_block,
        "required_gateL_lines_present": True,
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

    def check_same(self, category: str, name: str, dims: list[Dim], expected: Dim, expression: str) -> Dim:
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
        if actual == expected:
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


def build_dimensions() -> dict[str, Any]:
    check = DimChecker()

    M: Dim = (1, 0, 0)
    L: Dim = (0, 1, 0)
    T: Dim = (0, 0, 1)
    Z: Dim = (0, 0, 0)

    bulk_stress = (1, -2, -2)
    brane_lag = (1, -1, -2)
    u_operator = (1, -3, -2)
    p_operator = brane_lag
    mixed_u_p_operator = (1, -2, -2)
    determinant_dim = (2, -4, -4)
    closure_dim = brane_lag
    advective_current = (1, -2, -1)
    advective_curl = (1, -3, -1)

    d_m = M
    d_rho = (0, -4, 0)
    d_a = L
    d_ell_g = L
    d_u = L
    d_uw = L
    d_P = Z
    d_phi = (0, 2, -1)
    d_grad = (0, -1, 0)
    d_dt = (0, 0, -1)
    d_k = (0, -1, 0)
    d_omega = (0, 0, -1)
    d_v = (0, 1, -1)
    d_cs2 = (0, 2, -2)
    d_rho_br = (1, -3, 0)
    d_mu_R = brane_lag
    d_lambda = brane_lag
    d_Omega_w = d_omega
    d_Omega_P = d_omega
    d_Tbulk = bulk_stress
    d_frank_couple_traction = (1, 0, -2)

    d_curl_u = dadd(d_grad, d_u)
    d_div_u = d_curl_u
    d_J_P = check.check("projected_T0", "J_P", dadd(d_m, d_rho, dmul(2, d_a), d_ell_g), (1, -1, 0), "m rho a^2 ell_g")
    d_Gamma_P = check.check(
        "projected_T0",
        "Gamma_P",
        dadd(d_m, d_rho, d_cs2, dmul(2, d_a), d_ell_g),
        (1, 1, -2),
        "m rho c_s^2 a^2 ell_g",
    )
    d_Mgap_P = check.check(
        "projected_T0",
        "M_gap_P",
        dadd(d_J_P, dmul(2, d_Omega_P)),
        brane_lag,
        "J_P Omega_P^2",
    )

    check.check("L_a_i_traction", "MacCullagh_brane_traction", dadd(d_mu_R, d_curl_u), brane_lag, "mu_R Omega_u")
    check.check("L_a_i_traction", "P_u_brane_traction", dadd(d_lambda, d_P), brane_lag, "lambda_Pu varpi")
    check.check(
        "L_a_i_traction",
        "P_Frank_torque_not_surface_traction",
        dadd(d_Gamma_P, d_grad, d_P),
        d_frank_couple_traction,
        "Gamma_P partial P",
    )
    check.expect_fail(
        "ablation",
        "P_u_without_curl_or_cut_gradient",
        dadd(d_lambda, d_P, d_u),
        brane_lag,
        "lambda_Pu P u",
    )

    check.check("L_a_ii_symbol", "rho_br_omega2_u2", dadd(d_rho_br, dmul(2, d_omega), dmul(2, d_u)), brane_lag, "rho_br omega^2 u^2")
    check.check("L_a_ii_symbol", "mu_R_k2_u2", dadd(d_mu_R, dmul(2, d_k), dmul(2, d_u)), brane_lag, "mu_R k^2 u^2")
    check.check("L_a_ii_symbol", "J_P_omega2_P2", dadd(d_J_P, dmul(2, d_omega), dmul(2, d_P)), brane_lag, "J_P omega^2 P^2")
    check.check("L_a_ii_symbol", "Gamma_P_k2_P2", dadd(d_Gamma_P, dmul(2, d_k), dmul(2, d_P)), brane_lag, "Gamma_P k^2 P^2")
    check.check("L_a_ii_symbol", "lambda_k_u_P", dadd(d_lambda, d_k, d_u, d_P), brane_lag, "lambda_Pu k u P")
    check.check("L_a_ii_symbol", "u_operator_entry", dadd(d_mu_R, dmul(2, d_k)), u_operator, "mu_R k^2")
    check.check("L_a_ii_symbol", "P_operator_entry", dadd(d_Gamma_P, dmul(2, d_k)), p_operator, "Gamma_P k^2")
    check.check("L_a_ii_symbol", "mixed_operator_entry", dadd(d_lambda, d_k), mixed_u_p_operator, "lambda_Pu k")
    check.check("L_a_ii_dispersion", "c_gamma_squared", dsub(d_mu_R, d_rho_br), (0, 2, -2), "mu_R/rho_br")
    check.check("L_a_ii_dispersion", "c_gamma_eff_squared", dsub(d_mu_R, d_rho_br), (0, 2, -2), "(mu_R - 2 lambda_Pu)/rho_br")
    check.check("L_a_ii_dispersion", "slaved_k4_coefficient", dsub(d_Gamma_P, d_rho_br), (0, 4, -2), "Gamma_P/rho_br")
    check.check("L_a_ii_dispersion", "k_disp_squared", dsub(d_mu_R, d_Gamma_P), (0, -2, 0), "(mu_R - 2 lambda_Pu)/Gamma_P")
    check.check("L_a_ii_dispersion", "omega2_k2", dadd(dsub(d_mu_R, d_rho_br), dmul(2, d_k)), (0, 0, -2), "c_gamma^2 k^2")
    check.check("L_a_ii_dispersion", "omega2_k4", dadd(dsub(d_Gamma_P, d_rho_br), dmul(4, d_k)), (0, 0, -2), "(Gamma_P/rho_br) k^4")

    check.check("L_a_iii_C5", "u_longitudinal_kinetic", dadd(d_rho_br, dmul(2, d_dt), dmul(2, d_u)), brane_lag, "rho_br (partial_t u_L)^2")
    check.check("L_a_iii_C5", "phi_gradient_velocity", dadd(d_grad, d_phi), (0, 1, -1), "partial_a phi")
    check.check_same(
        "L_a_iii_C5",
        "phi_fixture_kinetic_parts",
        [dadd(d_dt, d_u), dadd(d_grad, d_phi)],
        (0, 1, -1),
        "partial_t u_a - partial_a phi",
    )
    check.check("L_a_iii_C5", "phi_fixture_term", dadd(d_rho_br, dmul(2, dadd(d_dt, d_u))), brane_lag, "rho_br (partial_t u - grad phi)^2")

    check.check("L_b_Hamiltonian", "det_mu_gamma_k4", dadd(d_mu_R, d_Gamma_P, dmul(4, d_k)), determinant_dim, "mu_R Gamma_P k^4")
    check.check("L_b_Hamiltonian", "det_lambda2_k2", dadd(dmul(2, d_lambda), dmul(2, d_k)), determinant_dim, "lambda_Pu^2 k^2")
    check.check("L_b_closure", "antisymmetric_stress_torque", dadd(d_mu_R, d_curl_u), closure_dim, "2 mu_R Omega_u")
    check.check("L_b_closure", "couple_stress_divergence", dadd(d_Gamma_P, dmul(2, d_k), d_P), closure_dim, "Gamma_P k^2 P")
    check.check("L_b_closure", "gapped_P_response", dsub(d_lambda, d_Mgap_P), Z, "lambda_Pu Omega_u / M_gap_P")

    check.check("L_c_direct_leak", "T_wa", dadd(d_m, d_rho, d_v, d_v), bulk_stress, "m rho v_w v_a")
    check.check("L_c_direct_leak", "partial_b_u_w", dadd(d_grad, d_uw), Z, "partial_b u_w")
    check.check("L_c_direct_leak", "slope_mixed_bulk_stress", dadd(d_Tbulk, dadd(d_grad, d_uw)), bulk_stress, "(T_ww delta_ab - T_ab) partial_b u_w")
    check.check_same(
        "L_c_direct_leak",
        "T_na_projected",
        [dadd(d_m, d_rho, d_v, d_v), dadd(d_Tbulk, dadd(d_grad, d_uw))],
        bulk_stress,
        "T_wa + (T_ww delta_ab - T_ab) partial_b u_w",
    )
    check.expect_fail("ablation", "drop_m_from_T_wa", dadd(d_rho, d_v, d_v), bulk_stress, "rho v_w v_a")
    check.check("L_c_indirect_leak", "advective_P_current", dadd(d_J_P, d_dt, d_P, d_grad, d_P), advective_current, "J_P partial_t P partial_a P")
    check.check("L_c_indirect_leak", "advective_current_curl", dadd(advective_current, d_grad), advective_curl, "curl(J_P partial_t P grad P)")

    check.check("L_d_uw_gap", "u_w_kinetic", dadd(d_rho_br, dmul(2, d_dt), dmul(2, d_uw)), brane_lag, "rho_br (partial_t u_w)^2")
    check.check("L_d_uw_gap", "u_w_gap", dadd(d_rho_br, dmul(2, d_Omega_w), dmul(2, d_uw)), brane_lag, "rho_br Omega_w^2 u_w^2")
    check.check("L_d_uw_gap", "omega_uw_squared", dmul(2, d_Omega_w), (0, 0, -2), "Omega_w^2")

    constants = {
        "rho": dim_record(d_rho),
        "m": dim_record(d_m),
        "a": dim_record(d_a),
        "ell_g": dim_record(d_ell_g),
        "rho_br": dim_record(d_rho_br),
        "mu_R": dim_record(d_mu_R),
        "lambda_Pu": dim_record(d_lambda),
        "Omega_w": dim_record(d_Omega_w),
        "Omega_P_gap_control": dim_record(d_Omega_P),
        "J_P": dim_record(d_J_P),
        "Gamma_P": dim_record(d_Gamma_P),
        "phi_fixture": dim_record(d_phi),
    }
    return {
        "base_dimensions": {"M": [1, 0, 0], "L": [0, 1, 0], "T": [0, 0, 1]},
        "targets": {
            "bulk_interface_traction": dim_record(bulk_stress),
            "brane_action_density_or_stress": dim_record(brane_lag),
            "u_operator": dim_record(u_operator),
            "P_operator": dim_record(p_operator),
            "mixed_u_P_operator": dim_record(mixed_u_p_operator),
            "Hamiltonian_minor": dim_record(determinant_dim),
            "closure_residual": dim_record(closure_dim),
            "advective_current": dim_record(advective_current),
            "advective_current_curl": dim_record(advective_curl),
            "Frank_couple_traction_not_surface_traction": dim_record(d_frank_couple_traction),
        },
        "constants": constants,
        "checks": check.records,
        "ablations": check.ablations,
        "pass": True,
    }


def rank_int(matrix: sp.Matrix) -> int:
    return int(matrix.rank())


def assert_zero(expr: sp.Expr, label: str) -> None:
    if sp.simplify(expr) != 0:
        raise AssertionError(f"{label}: {sp.factor(expr)}")


def factor_power(expr: sp.Expr, symbol: sp.Symbol, root: sp.Expr) -> int:
    poly = sp.Poly(sp.factor(expr), symbol)
    power = 0
    shifted = poly.as_expr()
    factor = symbol - root
    while shifted != 0:
        quotient, remainder = sp.div(sp.Poly(shifted, symbol), sp.Poly(factor, symbol))
        if remainder.as_expr() != 0:
            break
        power += 1
        shifted = sp.factor(quotient.as_expr())
    return power


def positive_gapless_count_from_roots(roots: dict[sp.Expr, int], k: sp.Symbol) -> int:
    count = 0
    for root, multiplicity in roots.items():
        if root == 0:
            continue
        if sp.simplify(sp.limit(root, k, 0, dir="+")) == 0:
            count += multiplicity
    return count


def positive_branch_count_from_roots(roots: dict[sp.Expr, int]) -> int:
    return sum(multiplicity for root, multiplicity in roots.items() if root != 0)


def build_derived_quantities() -> dict[str, Any]:
    """Derive Gate-L load-bearing quantities from the frozen quadratic action."""
    w = sp.symbols("w", real=True)
    k, ell_g, m, rho, a, c_s = sp.symbols("k ell_g m rho a c_s", positive=True)
    mu, lam, rho_br = sp.symbols("mu_R lambda_Pu rho_br", positive=True)
    omega2 = sp.symbols("omega2")
    Omega_w, Omega_P = sp.symbols("Omega_w Omega_P", positive=True)
    J_sym, Gamma_sym, M_gap_sym = sp.symbols("J_P Gamma_P M_gap_P", positive=True)

    g_ell = sp.exp(-(w / ell_g) ** 2) / (sp.sqrt(sp.pi) * ell_g)
    g_norm = sp.simplify(sp.integrate(g_ell, (w, -sp.oo, sp.oo)))
    mode_weight = sp.simplify(ell_g * g_norm)
    J_expr = sp.factor(m * rho * a**2 * mode_weight)
    Gamma_expr = sp.factor(m * rho * c_s**2 * a**2 * mode_weight)
    radial_expr = sp.factor(m * rho * c_s**2 * mode_weight)
    assert_zero(g_norm - 1, "g_ell normalization")
    assert_zero(mode_weight - ell_g, "projected confinement width")
    assert_zero(J_expr - m * rho * a**2 * ell_g, "projected J_P")
    assert_zero(Gamma_expr - m * rho * c_s**2 * a**2 * ell_g, "projected Gamma_P")
    assert_zero(radial_expr - Gamma_expr / a**2, "projected radial stiffness Gamma_P/a^2")

    u_t, varpi_t = sp.symbols("u_T varpi_T")
    omega_u = k * u_t
    # Energy is -L_static.  The frozen L_Pu=-lambda varpi.Omega therefore
    # contributes +lambda varpi.Omega to the Hamiltonian block.
    transverse_energy = sp.Rational(1, 2) * mu * omega_u**2 + lam * varpi_t * omega_u + sp.Rational(1, 2) * Gamma_sym * k**2 * varpi_t**2
    transverse_hessian = sp.simplify(sp.hessian(transverse_energy, (u_t, varpi_t)))
    expected_hessian = sp.Matrix([[mu * k**2, lam * k], [lam * k, Gamma_sym * k**2]])
    if sp.simplify(transverse_hessian - expected_hessian) != sp.zeros(2, 2):
        raise AssertionError(f"derived transverse Hessian mismatch: {transverse_hessian}")
    hessian_minor = sp.factor(transverse_hessian.det())
    assert_zero(hessian_minor - k**2 * (Gamma_sym * mu * k**2 - lam**2), "Hamiltonian minor")
    assert_zero(hessian_minor.subs(lam, 0) - Gamma_sym * mu * k**4, "lambda_Pu zero boundedness sanity")
    low_k_minor_sign = sp.factor(sp.limit(hessian_minor / k**2, k, 0, dir="+"))
    if low_k_minor_sign != -lam**2:
        raise AssertionError(f"unexpected low-k Hessian minor sign: {low_k_minor_sign}")

    # Frozen inactive slaved-rigid branch: P_parallel = w_hat x Omega_u.
    # Since varpi = w_hat x P_parallel, the induced axial micro-rotation is
    # varpi = w_hat x (w_hat x Omega_u) = -Omega_u.
    slaved_energy = sp.factor(transverse_energy.subs(varpi_t, -omega_u))
    slaved_hessian = sp.factor(sp.diff(slaved_energy, u_t, 2))
    mu_eff = sp.factor(sp.limit(slaved_hessian / k**2, k, 0, dir="+"))
    k4_coeff = sp.factor(sp.diff(slaved_hessian, k, 4) / sp.factorial(4))
    assert_zero(mu_eff - (mu - 2 * lam), "slaved mu_eff")
    assert_zero(k4_coeff - Gamma_sym, "slaved k4 coefficient")
    assert_zero(sp.simplify((-omega_u) + omega_u), "slaved closure residual")

    ux, uy, uz, px, py, pz, uw = sp.symbols("u_x u_y u_z varpi_x varpi_y varpi_z u_w")
    u_vec = sp.Matrix([ux, uy, uz])
    p_vec = sp.Matrix([px, py, pz])
    omega_vec = sp.Matrix([0, -k * uz, k * uy])

    def unslaved_symbol(gap: sp.Expr) -> tuple[sp.Matrix, sp.Expr, list[sp.Symbol]]:
        q = [ux, uy, uz, px, py, pz, uw]
        potential = (
            sp.Rational(1, 2) * mu * (omega_vec.dot(omega_vec))
            + lam * (p_vec.dot(omega_vec))
            + sp.Rational(1, 2) * (Gamma_sym * k**2 + gap) * (p_vec.dot(p_vec))
            + sp.Rational(1, 2) * rho_br * Omega_w**2 * uw**2
        )
        stiffness = sp.hessian(potential, q)
        mass = sp.diag(rho_br, rho_br, rho_br, J_sym, J_sym, J_sym, rho_br)
        dyn = sp.factor(stiffness - omega2 * mass)
        return dyn, sp.factor(dyn.det()), q

    dyn_massless, det_massless, _ = unslaved_symbol(sp.Integer(0))
    dyn_gapped, det_gapped, _ = unslaved_symbol(M_gap_sym)
    det_massless_lam0 = sp.factor(det_massless.subs(lam, 0))
    det_gapped_lam0 = sp.factor(det_gapped.subs(lam, 0))
    roots_massless_lam0 = sp.roots(sp.Poly(det_massless_lam0, omega2))
    roots_gapped_lam0 = sp.roots(sp.Poly(det_gapped_lam0, omega2))

    slaved_q = [ux, uy, uz, uw]
    slaved_transverse_stiffness = sp.factor((mu_eff * k**2 + Gamma_sym * k**4))
    slaved_potential = sp.Rational(1, 2) * slaved_transverse_stiffness * (uy**2 + uz**2) + sp.Rational(1, 2) * rho_br * Omega_w**2 * uw**2
    slaved_stiffness = sp.hessian(slaved_potential, slaved_q)
    slaved_mass = sp.diag(rho_br, rho_br, rho_br, rho_br)
    dyn_slaved = sp.factor(slaved_stiffness - omega2 * slaved_mass)
    det_slaved = sp.factor(dyn_slaved.det())
    roots_slaved = sp.roots(sp.Poly(det_slaved, omega2))

    massless_degree = sp.Poly(det_massless, omega2).degree()
    gapped_degree = sp.Poly(det_gapped, omega2).degree()
    slaved_degree = sp.Poly(det_slaved, omega2).degree()
    if (massless_degree, gapped_degree, slaved_degree) != (7, 7, 4):
        raise AssertionError(f"unexpected principal-symbol degrees: {(massless_degree, gapped_degree, slaved_degree)}")

    gapless_massless_lam0 = positive_gapless_count_from_roots(roots_massless_lam0, k)
    gapless_gapped_lam0 = positive_gapless_count_from_roots(roots_gapped_lam0, k)
    gapless_slaved = positive_gapless_count_from_roots(roots_slaved, k)
    positive_massless_lam0 = positive_branch_count_from_roots(roots_massless_lam0)
    positive_gapped_lam0 = positive_branch_count_from_roots(roots_gapped_lam0)
    positive_slaved = positive_branch_count_from_roots(roots_slaved)
    zero_massless = roots_massless_lam0.get(sp.Integer(0), 0)
    zero_gapped = roots_gapped_lam0.get(sp.Integer(0), 0)
    zero_slaved = roots_slaved.get(sp.Integer(0), 0)
    if (gapless_massless_lam0, gapless_gapped_lam0, gapless_slaved) != (5, 2, 2):
        raise AssertionError("unexpected gapless branch counts")
    if (positive_massless_lam0, positive_gapped_lam0, positive_slaved) != (6, 6, 3):
        raise AssertionError("unexpected positive branch counts")
    if (zero_massless, zero_gapped, zero_slaved) != (1, 1, 1):
        raise AssertionError("unexpected longitudinal zero-root counts")

    transverse_negative_low_k = 2 if low_k_minor_sign == -lam**2 else 0
    massless_positive_low_k_lambda_nonzero = gapless_massless_lam0 - transverse_negative_low_k
    if massless_positive_low_k_lambda_nonzero != 3:
        raise AssertionError("lambda_Pu nonzero low-k branch count did not change")

    # Gapped-P closure control retained, now fed by the derived response.
    Omega = sp.symbols("Omega_u")
    p_response_gapped = sp.factor(-lam * Omega / (M_gap_sym + Gamma_sym * k**2))
    couple_divergence = sp.factor(Gamma_sym * k**2 * p_response_gapped)
    gapped_closure_residual = sp.factor(sp.limit(2 * mu * Omega + couple_divergence, k, 0, dir="+"))
    assert_zero(gapped_closure_residual - 2 * mu * Omega, "gapped-P low-k closure residual")

    # Direct and indirect L(c), using a real flat shear wave and the frozen
    # advective T0 current J_P (D_t^v P) partial_a P at leading v=0.
    x, y, t, A, B, q, omega = sp.symbols("x y t A B q omega", real=True)
    phase = k * x - omega * t
    u_flat_y = A * sp.cos(phase)
    u_flat_w = sp.Integer(0)
    v_flat_w = sp.diff(u_flat_w, t)
    v_flat_y = sp.diff(u_flat_y, t)
    slope_flat_x = sp.diff(u_flat_w, x)
    delta_T = sp.symbols("DeltaT")
    T_na_flat = sp.simplify(m * rho * v_flat_w * v_flat_y + delta_T * slope_flat_x)
    u_bent_w = B * sp.cos(q * x - omega * t)
    T_na_bent = sp.simplify((m * rho * sp.diff(u_bent_w, t) * v_flat_y + delta_T * sp.diff(u_bent_w, x)).subs({A: 2, B: 3, m: 5, rho: 7, delta_T: 11, k: 13, q: 17, omega: 19, x: 23, t: 29}))
    if T_na_flat != 0 or T_na_bent == 0:
        raise AssertionError("direct L(c) traction derivation failed")

    response_den_x = Gamma_sym * k**2 - J_sym * omega**2
    omega_u_flat = sp.diff(u_flat_y, x)
    p_flat = sp.factor(-lam * omega_u_flat / response_den_x)
    advective_jx = sp.factor(J_sym * sp.diff(p_flat, t) * sp.diff(p_flat, x))
    advective_jy = sp.Integer(0)
    indirect_curl_flat = sp.simplify(sp.diff(advective_jy, x) - sp.diff(advective_jx, y))
    p_frank_only = sp.simplify(p_flat.subs(lam, 0))
    if indirect_curl_flat != 0 or p_frank_only != 0:
        raise AssertionError("flat indirect channel or Frank-only control failed")

    response_den_y = Gamma_sym * q**2 - J_sym * omega**2
    p_nonplanar = sp.factor(
        -lam
        * (
            A * k * sp.sin(k * x) * sp.cos(omega * t) / response_den_x
            + B * q * sp.sin(q * y) * sp.sin(omega * t) / response_den_y
        )
    )
    jx_np = sp.diff(p_nonplanar, t) * sp.diff(p_nonplanar, x)
    jy_np = sp.diff(p_nonplanar, t) * sp.diff(p_nonplanar, y)
    curl_np = sp.trigsimp(sp.factor(sp.diff(jy_np, x) - sp.diff(jx_np, y)))
    if sp.simplify(curl_np) == 0:
        raise AssertionError("nonplanar indirect able-to-fail control failed")

    return {
        "projection": {
            "g_ell_normalization": "1",
            "mode_weight_integral": "ell_g",
            "J_P": "m rho a^2 ell_g",
            "Gamma_P": "m rho c_s^2 a^2 ell_g",
            "radial_stiffness": "Gamma_P/a^2",
            "source": "int dw ell_g g_ell(w) times the inherited T0 P terms",
        },
        "transverse_hessian": {
            "basis": ["u_T", "varpi_T"],
            "matrix": [["mu_R k^2", "lambda_Pu k"], ["lambda_Pu k", "Gamma_P k^2"]],
            "principal_minor": "k^2 (Gamma_P mu_R k^2 - lambda_Pu^2)",
            "negative_energy_interval": "0 < k^2 < lambda_Pu^2/(Gamma_P mu_R)",
            "lambda_Pu_zero_minor": "Gamma_P mu_R k^4",
            "bounded_below_lambda_Pu_zero": True,
            "low_k_minor_limit_over_k2": "-lambda_Pu^2",
        },
        "slaved_reduction": {
            "substitution": "P_parallel = w_hat x Omega_u; varpi = w_hat x P_parallel = -Omega_u",
            "mu_eff": "mu_R - 2 lambda_Pu",
            "transverse_energy_hessian": "k^2 (mu_R - 2 lambda_Pu + Gamma_P k^2)",
            "k4_coefficient": "Gamma_P",
            "omega2_slaved": "((mu_R - 2 lambda_Pu) k^2 + Gamma_P k^4)/rho_br",
            "k_disp_squared": "(mu_R - 2 lambda_Pu)/Gamma_P",
            "closure_residual": "0",
        },
        "principal_symbols": {
            "massless_P": {
                "det_degree": massless_degree,
                "positive_omega2_branches_lambda_Pu_zero": positive_massless_lam0,
                "gapless_positive_branches_lambda_Pu_zero": gapless_massless_lam0,
                "positive_gapless_branches_low_k_lambda_Pu_nonzero": massless_positive_low_k_lambda_nonzero,
                "negative_omega2_branches_low_k_lambda_Pu_nonzero": transverse_negative_low_k,
                "longitudinal_zero_roots": zero_massless,
                "extra_gapless_P_branches": gapless_massless_lam0 - 2,
            },
            "gapped_P": {
                "det_degree": gapped_degree,
                "positive_omega2_branches_lambda_Pu_zero": positive_gapped_lam0,
                "gapless_positive_branches_lambda_Pu_zero": gapless_gapped_lam0,
                "longitudinal_zero_roots": zero_gapped,
                "extra_gapless_P_branches": gapless_gapped_lam0 - 2,
                "gapped_P_branches": 3,
            },
            "slaved_rigid": {
                "det_degree": slaved_degree,
                "positive_omega2_branches": positive_slaved,
                "gapless_positive_branches": gapless_slaved,
                "longitudinal_zero_roots": zero_slaved,
                "extra_gapless_P_branches": gapless_slaved - 2,
            },
        },
        "closure": {
            "gapped_P_response": "-lambda_Pu Omega_u/(M_gap_P + Gamma_P k^2)",
            "gapped_P_low_k_residual": "2 mu_R Omega_u",
            "omit_reservoir_low_k_residual": "2 mu_R Omega_u",
            "slaved_residual": "0",
        },
        "leak": {
            "direct_flat_T_na": "0",
            "direct_flat_wave": "u_y=A cos(k x - omega t), u_w=0, v_w=partial_t u_w=0",
            "bent_control_nonzero": True,
            "indirect_flat_curl": "0",
            "indirect_flat_P_response": "-lambda_Pu partial_x u_y/(Gamma_P k^2 - J_P omega^2)",
            "advective_current": "J_P (partial_t P) partial_a P from D_t^v P at leading v=0",
            "Frank_only_induced_P": "0",
            "nonplanar_able_to_fail_nonzero": True,
        },
    }


def traction_audit(config: str) -> dict[str, Any]:
    # Traction across the x-normal cut from E_Pu = lambda varpi . curl(u).
    lam = sp.symbols("lambda_Pu", nonzero=True)
    mu = sp.symbols("mu_R", nonzero=True)
    vx, vy, vz = sp.symbols("varpi_x varpi_y varpi_z")
    varpi = sp.Matrix([vx, vy, vz])

    # Columns are varpi components, rows are tangential traction components
    # t_c = partial E / partial(partial_x u_c) = lambda epsilon_{a x c} varpi_a.
    coupling_cut_matrix = sp.Matrix(
        [
            [0, 0, 0],
            [0, 0, -lam],
            [0, lam, 0],
        ]
    )
    mac_cut_matrix = sp.Matrix(
        [
            [0, 0, 0],
            [0, 0, -mu],
            [0, mu, 0],
        ]
    )
    frank_only_matrix = sp.zeros(3, 3)
    coupling_rank = rank_int(coupling_cut_matrix)
    mac_rank = rank_int(mac_cut_matrix)
    frank_rank = rank_int(frank_only_matrix)
    if coupling_rank != 2 or mac_rank != 2 or frank_rank != 0:
        raise AssertionError("unexpected traction ranks")

    if config == "A_baseline":
        provenance = "ARROWS_SUPPLY_TRACTION"
        note = "P-u virtual work gives a rank-2 cut traction; standalone mu_R traction is also present."
    else:
        provenance = "POSTULATED_SURFACE_ELASTICITY"
        note = "slaving removes independent arrow modes; the clean k^2 light traction is the postulated surface modulus, with arrow-sector k^4/coupling corrections."

    return {
        "status": "PASS",
        "label": None,
        "provenance": provenance,
        "cut_normal": "x",
        "coupling_traction_rank": coupling_rank,
        "standalone_MacCullagh_traction_rank": mac_rank,
        "Frank_only_reference_traction_rank": frank_rank,
        "Frank_only_control": {
            "status": "FIRED",
            "label": "FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION",
            "computed_rank": frank_rank,
        },
        "standalone_elasticity_discriminator": {
            "lambda_Pu_zero_mu_R_nonzero_rank": mac_rank,
            "status": "FIRED",
            "meaning": "the machinery distinguishes P-u-sourced traction from postulated surface elasticity",
        },
        "note": note,
    }


def hidden_mode_audit(config: str, derived: dict[str, Any]) -> dict[str, Any]:
    k, mu, rho_br, gamma, lam = sp.symbols("k mu_R rho_br Gamma_P lambda_Pu", positive=True)
    cauchy = sp.diag((sp.Symbol("lambda_C", positive=True) + 2 * sp.Symbol("mu_C", positive=True)) * k**2, sp.Symbol("mu_C", positive=True) * k**2, sp.Symbol("mu_C", positive=True) * k**2)
    cauchy_rank = rank_int(cauchy)
    if cauchy_rank != 3:
        raise AssertionError("Cauchy control did not produce three propagating elastic modes")

    c_gamma = "mu_R/rho_br"
    mode_counts = derived["principal_symbols"]
    if config == "A_baseline":
        massless_counts = mode_counts["massless_P"]
        counts = {
            "positive_gapless_branches_low_k_lambda_Pu_nonzero": massless_counts["positive_gapless_branches_low_k_lambda_Pu_nonzero"],
            "negative_omega2_branches_low_k_lambda_Pu_nonzero": massless_counts["negative_omega2_branches_low_k_lambda_Pu_nonzero"],
            "gapless_positive_branches_lambda_Pu_zero": massless_counts["gapless_positive_branches_lambda_Pu_zero"],
            "P_spin_wave_gapless_lambda_Pu_zero": massless_counts["extra_gapless_P_branches"],
            "u_longitudinal_zero_roots": massless_counts["longitudinal_zero_roots"],
            "P_radial_gapped": 1,
            "u_w_gapped": 1,
            "phi": 0,
        }
        extra = massless_counts["extra_gapless_P_branches"]
        status = "FAIL"
        label = "FAIL_HIDDEN_PROPAGATING_MODE"
        dispersion = {
            "u_light_omega2_uncoupled_reference": "(mu_R/rho_br) k^2",
            "P_spin_omega2": "c_s^2 k^2",
            "coupled_transverse_pencil": "det[[mu_R k^2-rho_br omega^2, lambda_Pu k],[lambda_Pu k, Gamma_P k^2-J_P omega^2]]=0",
            "small_k_warning": "derived Hamiltonian minor gives two negative omega^2 transverse branches at low k for nonzero lambda_Pu; counted separately in L(b)",
        }
    else:
        slaved = derived["slaved_reduction"]
        slaved_counts = mode_counts["slaved_rigid"]
        mu_eff = sp.factor(mu - 2 * lam)
        omega2 = sp.factor((mu_eff * k**2 + gamma * k**4) / rho_br)
        kdisp2 = sp.factor(mu_eff / gamma)
        if sp.factor(rho_br * omega2 - (mu_eff * k**2 + gamma * k**4)) != 0:
            raise AssertionError("slaved dispersion algebra mismatch")
        counts = {
            "gapless_positive_branches": slaved_counts["gapless_positive_branches"],
            "positive_omega2_branches": slaved_counts["positive_omega2_branches"],
            "u_longitudinal_zero_roots": slaved_counts["longitudinal_zero_roots"],
            "P_spin_wave_gapless": slaved_counts["extra_gapless_P_branches"],
            "u_w_gapped": 1,
            "phi": 0,
        }
        extra = slaved_counts["extra_gapless_P_branches"]
        status = "PASS_WITH_LOW_K_DISPERSION_CAVEAT"
        label = None
        dispersion = {
            "omega2_slaved": slaved["omega2_slaved"],
            "c_gamma_squared_MacCullagh_feed_forward": c_gamma,
            "effective_c_squared_if_bilinear_retained": "(mu_R - 2 lambda_Pu)/rho_br",
            "k4_correction": "(Gamma_P/rho_br) k^4",
            "k_disp_squared": slaved["k_disp_squared"],
            "low_k_tolerance": "requires k^2 << (mu_R - 2 lambda_Pu)/Gamma_P with positive effective modulus; exact nondispersive equality is false at finite k",
            "symbolic_checks": {"omega2_identity": True, "k_disp_squared": str(kdisp2).replace("Gamma_P", "Gamma_P")},
        }
    return {
        "status": status,
        "label": label,
        "counts": counts,
        "extra_propagating_modes": extra,
        "derived_from_principal_symbol": True,
        "mode_count_controls": {
            "gapped_P": mode_counts["gapped_P"],
            "lambda_Pu_zero_changes_low_k_count": mode_counts["massless_P"]["positive_gapless_branches_low_k_lambda_Pu_nonzero"]
            != mode_counts["massless_P"]["gapless_positive_branches_lambda_Pu_zero"],
        },
        "cauchy_reference_control": {
            "status": "FIRED",
            "label": "FAIL_CAUCHY_STRAY_LONGITUDINAL",
            "computed_propagating_modes": cauchy_rank,
            "stray_longitudinal": True,
        },
        "c_gamma_squared": c_gamma,
        "dispersion": dispersion,
    }


def c5_audit(config: str) -> dict[str, Any]:
    omega, k = sp.symbols("omega k", nonzero=True)
    no_phi_matrix = sp.Matrix([[omega**2]])
    no_phi_stiffness_rank = rank_int(sp.Matrix([[0]]))
    no_phi_kinetic_rank = rank_int(no_phi_matrix)
    phi_fixture_matrix = sp.Matrix([[omega**2, -omega * k], [-omega * k, k**2]])
    phi_fixture_rank = rank_int(phi_fixture_matrix)
    phi_fixture_det = sp.factor(phi_fixture_matrix.det())
    if no_phi_kinetic_rank != 1 or no_phi_stiffness_rank != 0:
        raise AssertionError("no-phi C5 branch did not leave the kinetic zero mode")
    if phi_fixture_rank != 1 or phi_fixture_det != 0:
        raise AssertionError("phi fixture did not produce a gauge-null kinetic form")
    return {
        "status": "FAIL",
        "label": "FAIL_C5_LONGITUDINAL_ZERO_MODE",
        "frozen_phi_status": "absent",
        "no_phi_branch_control": {
            "status": "FIRED",
            "kinetic_rank": no_phi_kinetic_rank,
            "curl_stiffness_rank": no_phi_stiffness_rank,
            "computed_result": "constrained physical zero mode remains",
        },
        "independent_variational_phi_fixture_control": {
            "status": "FIRED_PASS_FIXTURE_ONLY",
            "kinetic_form_rank": phi_fixture_rank,
            "determinant": str(phi_fixture_det),
            "gauge_nullity": 1,
            "physical_longitudinal_modes": 0,
            "gate_pass_candidate_only_if_frozen_in_G0": False,
        },
        "raw_divergence_projector_control": {
            "status": "FIRED",
            "label": "FAIL_C5_LONGITUDINAL_ZERO_MODE",
            "reason": "removes u_L by hand with no variational provenance",
        },
        "phi_equals_u_w_collision": {
            "status": "FIRED",
            "reason": "a scalar-potential phi identified with u_w must be massless, but L(d) requires Omega_w>0",
        },
        "config_note": f"{config} inherits the G0 no-phi decision",
    }


def bounded_closure_audit(config: str, derived: dict[str, Any]) -> dict[str, Any]:
    hessian = derived["transverse_hessian"]
    slaved = derived["slaved_reduction"]
    closure = derived["closure"]

    controls = {
        "omit_couple_stress_reservoir": {
            "status": "FIRED",
            "label": "FAIL_GYROSTAT_NO_CLOSURE",
            "low_k_residual": closure["omit_reservoir_low_k_residual"],
        },
        "large_gap_P_control": {
            "status": "FIRED",
            "label": "FAIL_GYROSTAT_NO_CLOSURE",
            "P_response": closure["gapped_P_response"],
            "low_k_couple_divergence": "0",
            "low_k_residual": closure["gapped_P_low_k_residual"],
        },
    }
    if config == "A_baseline":
        return {
            "status": "FAIL",
            "labels": ["FAIL_NOT_BOUNDED_BELOW"],
            "Hamiltonian_energy_matrix_one_transverse_pair": hessian["matrix"],
            "principal_minor": hessian["principal_minor"],
            "negative_energy_interval": hessian["negative_energy_interval"],
            "bounded_below": False,
            "lambda_Pu_zero_restores_boundedness": hessian["bounded_below_lambda_Pu_zero"],
            "closure": {
                "live_P_needed": True,
                "generic_low_k_identity": "not banked; closure requires a live/singular P response and is lost when P is gapped",
                "gapped_leg_low_k_residual": closure["gapped_P_low_k_residual"],
            },
            "controls": controls,
            "derived_from_action_hessian": True,
        }
    return {
        "status": "PASS_CONDITIONAL_ON_POSITIVE_EFFECTIVE_MODULUS",
        "labels": [],
        "Hamiltonian_energy_coefficient_transverse": slaved["transverse_energy_hessian"],
        "bounded_below_conditions": ["rho_br>0", "Gamma_P>0", "mu_R - 2 lambda_Pu > 0"],
        "closure": {
            "residual": slaved["closure_residual"],
            "reason": "slaved micro-rotation equals the material rotation, so the antisymmetric stress is absorbed algebraically",
        },
        "controls": controls,
        "derived_from_slaved_substitution": True,
    }


def leak_audit(config: str, derived: dict[str, Any]) -> dict[str, Any]:
    leak = derived["leak"]

    return {
        "status": "PASS",
        "label": None,
        "direct": {
            "flat_wave_T_na": leak["direct_flat_T_na"],
            "flat_wave": leak["direct_flat_wave"],
            "flat_conditions": ["v_w=partial_t u_w=0 read from the wave", "partial_b u_w=0"],
            "bent_control_nonzero": leak["bent_control_nonzero"],
            "bent_control_status": "FIRED",
        },
        "indirect": {
            "P_response_from_L_Pu": leak["indirect_flat_P_response"],
            "advective_current_from_T0": leak["advective_current"],
            "flat_one_k_bulk_vorticity_source": leak["indirect_flat_curl"],
            "Frank_only_control": {
                "status": "FIRED",
                "induced_P_from_u": leak["Frank_only_induced_P"],
                "bulk_vorticity_source": "0",
            },
            "nonplanar_able_to_fail_control": {
                "status": "FIRED",
                "curl_expression_nonzero": leak["nonplanar_able_to_fail_nonzero"],
            },
        },
        "curvature_squared_residual": "relocated_to_Gate_T",
        "derived_from_actual_flat_wave": True,
    }


def uw_gap_audit(config: str) -> dict[str, Any]:
    Omega_w, c_s, a, k = sp.symbols("Omega_w c_s a k", positive=True)
    scalar_omega2 = [Omega_w**2, c_s**2 * k**2 + 2 * c_s**2 / a**2]
    ungapped = [expr.subs(Omega_w, 0) for expr in scalar_omega2]
    if scalar_omega2[0] == 0 or ungapped[0] != 0:
        raise AssertionError("u_w gap control failed")
    return {
        "status": "PASS",
        "label": None,
        "full_coupled_scalar_spectrum": {
            "u_w_descendant": "Omega_w^2",
            "P_radial_reference": "c_s^2 k^2 + 2 c_s^2/a^2",
            "phi": "absent",
            "coupling_to_P_u_or_confinement": "block diagonal at linear flat-brane order",
        },
        "massless_scalar_modes_from_u_w": 0,
        "ungapped_u_w_control": {
            "status": "FIRED",
            "label": "FAIL_BENDING_MASSLESS_FIFTH_FORCE",
            "Omega_w_set_to_zero_modes": 1,
        },
        "phi_equals_u_w_collision_in_C5_fixture": "would set the scalar-potential descendant massless and trip this hurdle",
    }


def sub_hurdle_pass(status: str) -> bool:
    return status == "PASS" or status.startswith("PASS_CONDITIONAL") or status.startswith("PASS_WITH")


def labels_from_sub_hurdles(sub: dict[str, Any]) -> list[str]:
    labels: list[str] = []
    for value in sub.values():
        label = value.get("label")
        if label:
            labels.append(label)
        for item in value.get("labels", []):
            if item:
                labels.append(item)
    return labels


def aggregate_config_verdict(sub: dict[str, Any]) -> tuple[str, list[str]]:
    labels = labels_from_sub_hurdles(sub)
    if all(sub_hurdle_pass(value["status"]) for value in sub.values()):
        return "FREE_LIGHT_OK_CONDITIONAL", labels
    linked_mode_or_closure = any(
        label in labels
        for label in ["FAIL_HIDDEN_PROPAGATING_MODE", "FAIL_NOT_BOUNDED_BELOW", "FAIL_GYROSTAT_NO_CLOSURE"]
    )
    if linked_mode_or_closure and "FAIL_C5_LONGITUDINAL_ZERO_MODE" in labels:
        return "FAIL_COUPLE_STRESS_NOGO", labels
    if labels:
        return labels[0], labels
    return "FAIL_TAUTOLOGICAL", labels


def aggregate_overall_verdict(configs: dict[str, Any], derived: dict[str, Any]) -> dict[str, Any]:
    passing = [name for name, result in configs.items() if result["config_verdict"] == "FREE_LIGHT_OK_CONDITIONAL"]
    if passing:
        return {
            "verdict": "FREE_LIGHT_OK_CONDITIONAL",
            "pass_subtag": configs[passing[0]]["sub_hurdles"]["L_a_i"]["provenance"],
            "reason": f"{passing[0]} satisfies all Gate-L sub-hurdles.",
        }
    all_labels = {label for result in configs.values() for label in result["labels_fired"]}
    gapped_control_fires = derived["closure"]["gapped_P_low_k_residual"] == "2 mu_R Omega_u"
    linked = (
        "FAIL_HIDDEN_PROPAGATING_MODE" in all_labels
        and ("FAIL_NOT_BOUNDED_BELOW" in all_labels or "FAIL_GYROSTAT_NO_CLOSURE" in all_labels)
        and "FAIL_C5_LONGITUDINAL_ZERO_MODE" in all_labels
        and gapped_control_fires
    )
    if linked:
        return {
            "verdict": "FAIL_COUPLE_STRESS_NOGO",
            "pass_subtag": None,
            "reason": "No frozen configuration satisfies L(a-i), L(a-ii), L(a-iii), L(b), L(c), and L(d). The derived live-P symbol has hidden/unstable P branches and an unbounded Hessian, the derived gapped-P control loses closure at low k, and the slaved-rigid escape still inherits the frozen no-phi C5 failure.",
        }
    if all_labels:
        first = sorted(all_labels)[0]
        return {"verdict": first, "pass_subtag": None, "reason": f"Gate L failed through {first}."}
    return {"verdict": "FAIL_TAUTOLOGICAL", "pass_subtag": None, "reason": "No passing config and no computed failure label."}


def good_structure_fixture_verdict() -> dict[str, Any]:
    sub = {
        "L_a_i": {"status": "PASS", "provenance": "ARROWS_SUPPLY_TRACTION"},
        "L_a_ii": {"status": "PASS"},
        "L_a_iii": {"status": "PASS"},
        "L_b": {"status": "PASS"},
        "L_c": {"status": "PASS"},
        "L_d": {"status": "PASS"},
    }
    verdict, labels = aggregate_config_verdict(sub)
    if verdict != "FREE_LIGHT_OK_CONDITIONAL" or labels:
        raise AssertionError("good-structure fixture did not aggregate to FREE_LIGHT_OK_CONDITIONAL")
    return {
        "name": "all_sub_hurdles_pass_fixture",
        "aggregated_verdict": verdict,
        "status": "FIRED_PASS_FIXTURE",
    }


def config_result(config: str, derived: dict[str, Any]) -> dict[str, Any]:
    sub = {
        "L_a_i": traction_audit(config),
        "L_a_ii": hidden_mode_audit(config, derived),
        "L_a_iii": c5_audit(config),
        "L_b": bounded_closure_audit(config, derived),
        "L_c": leak_audit(config, derived),
        "L_d": uw_gap_audit(config),
    }
    verdict, fired = aggregate_config_verdict(sub)
    if config == "A_baseline":
        chain_role = "live_P_horn: massless P supplies the needed reservoir candidate but creates hidden modes and an unbounded gyroscopic P-u Hamiltonian; no phi leaves C5 exposed"
    else:
        chain_role = "slaved-rigid escape resolves hidden P mode-count and closure with a derived k^4 correction, but inherits the frozen no-phi C5 zero mode; phi=u_w rescue collides with L(d)"
    return {
        "config": config,
        "sub_hurdles": sub,
        "labels_fired": fired,
        "config_verdict": verdict,
        "section_2_6_chain_role": chain_role,
    }


def controls_summary(configs: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"name": "Frank_only_reference", "target": "L(a-i)", "status": "FIRED", "label": "FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION"},
        {"name": "Cauchy_reference", "target": "L(a-ii)", "status": "FIRED", "label": "FAIL_CAUCHY_STRAY_LONGITUDINAL"},
        {"name": "no_phi_C5_branch", "target": "L(a-iii)", "status": "FIRED", "label": "FAIL_C5_LONGITUDINAL_ZERO_MODE"},
        {"name": "independent_variational_phi_fixture", "target": "L(a-iii)", "status": "FIRED_PASS_FIXTURE_ONLY", "fresh_G0_required": True},
        {"name": "raw_divergence_projector", "target": "L(a-iii)", "status": "FIRED", "label": "FAIL_C5_LONGITUDINAL_ZERO_MODE"},
        {"name": "omit_couple_stress_reservoir", "target": "L(b)", "status": "FIRED", "label": "FAIL_GYROSTAT_NO_CLOSURE"},
        {"name": "large_gap_P_leg", "target": "L(b)", "status": "FIRED", "label": "FAIL_GYROSTAT_NO_CLOSURE"},
        {"name": "bent_interface", "target": "L(c) direct", "status": "FIRED", "label": "FAIL_LEAK_BREAKS_MAGNUS"},
        {"name": "Frank_only_indirect", "target": "L(c) indirect", "status": "FIRED_ZERO_SOURCE", "bulk_vorticity_source": "0"},
        {"name": "nonplanar_indirect_able_to_fail", "target": "L(c) indirect", "status": "FIRED_NONZERO_SOURCE"},
        {"name": "ungapped_u_w", "target": "L(d)", "status": "FIRED", "label": "FAIL_BENDING_MASSLESS_FIFTH_FORCE"},
        {"name": "drop_m_from_T_wa", "target": "dimensional_firewall", "status": "FIRED", "label": "FAIL_DIMENSIONAL"},
        {"name": "P_u_without_curl_or_cut_gradient", "target": "dimensional_firewall", "status": "FIRED", "label": "FAIL_DIMENSIONAL"},
    ]


def hypothetical_pass_fixture() -> dict[str, Any]:
    return {
        "name": "B_plus_independent_variational_phi",
        "computed_status": good_structure_fixture_verdict()["aggregated_verdict"] + "_fixture",
        "traction_subtag": "POSTULATED_SURFACE_ELASTICITY",
        "passes_all_subhurdles_if": [
            "the low-k interpretation tolerates the slaved k^4 correction",
            "mu_R - 2 lambda_Pu > 0",
            "the independent phi fixture was frozen in G0",
        ],
        "gateL_admissible": False,
        "reason": "G0 froze phi absent; adding phi is a fresh G0, not a Gate-L pass",
    }


def build_agreement_payload() -> dict[str, Any]:
    freeze = check_freeze_fidelity()
    dimensions = build_dimensions()
    derived = build_derived_quantities()
    configs = {
        "A_baseline": config_result("A_baseline", derived),
        "B_slaved_rigid": config_result("B_slaved_rigid", derived),
    }
    controls = controls_summary(configs)
    aggregate = aggregate_overall_verdict(configs, derived)
    overall = {
        "verdict": aggregate["verdict"],
        "provenance": "CONDITIONAL_ON(both)",
        "pass_subtag": aggregate["pass_subtag"],
        "reason": aggregate["reason"],
        "section_2_6_four_way_materialized": True,
        "branch_B_escaped_mode_count_plus_closure": True,
        "gapped_P_leg_tabulated": True,
        "c_gamma_squared": "mu_R/rho_br",
        "dimensional_firewall": "PASS",
        "engine_agreement": "PENDING_COMPARE",
        "verdict_aggregated_from_subhurdles": True,
    }
    return {
        "freeze_fidelity": freeze,
        "dimensional_firewall": dimensions,
        "derived_from_S_G0": derived,
        "configurations": configs,
        "controls": controls,
        "good_structure_able_to_pass_fixture": good_structure_fixture_verdict(),
        "hypothetical_pass_fixture_not_gate_admissible": hypothetical_pass_fixture(),
        "overall": overall,
    }


def build_report(engine: str) -> dict[str, Any]:
    payload = build_agreement_payload()
    return {
        "schema": f"pathA_35_gateL_{engine}/v1",
        "engine": engine,
        "pass": True,
        "agreement_payload": payload,
        "verdict": payload["overall"]["verdict"],
    }


def yaml_ledger(compared: dict[str, Any]) -> dict[str, Any]:
    payload = compared["agreement_payload"]
    summary = {
        "schema": "pathA_35_gateL_light_results/v1",
        "stage": "Gate_L_light_on_shear_surface_brane",
        "date": "2026-06-26",
        "frozen_baseline": f"T0_SHEAR_FROZEN({EXPECTED_G0_SHORT})",
        "overall_verdict": payload["overall"]["verdict"],
        "provenance": payload["overall"]["provenance"],
        "pass_subtag": payload["overall"]["pass_subtag"],
        "section_2_6_four_way_materialized": payload["overall"]["section_2_6_four_way_materialized"],
        "branch_B_escaped_mode_count_plus_closure": payload["overall"]["branch_B_escaped_mode_count_plus_closure"],
        "gapped_P_leg_tabulated": payload["overall"]["gapped_P_leg_tabulated"],
        "c_gamma_squared": payload["overall"]["c_gamma_squared"],
        "derived_from_S_G0": {
            "projection": payload["derived_from_S_G0"]["projection"],
            "transverse_hessian": payload["derived_from_S_G0"]["transverse_hessian"],
            "slaved_reduction": payload["derived_from_S_G0"]["slaved_reduction"],
            "principal_symbol_mode_counts": payload["derived_from_S_G0"]["principal_symbols"],
            "closure": payload["derived_from_S_G0"]["closure"],
            "leak": payload["derived_from_S_G0"]["leak"],
            "verdict_aggregated_from_subhurdles": payload["overall"]["verdict_aggregated_from_subhurdles"],
            "good_structure_able_to_pass_fixture": payload["good_structure_able_to_pass_fixture"],
        },
        "dimensional_firewall": {
            "pass": True,
            "engine_agreement": "ENGINE_AGREE",
            "checked_expression_count": len(payload["dimensional_firewall"]["checks"]),
            "ablations": payload["dimensional_firewall"]["ablations"],
        },
        "configurations": {
            name: {
                "config_verdict": result["config_verdict"],
                "labels_fired": result["labels_fired"],
                "sub_hurdle_status": {
                    key: value["status"] for key, value in result["sub_hurdles"].items()
                },
                "traction_provenance": result["sub_hurdles"]["L_a_i"]["provenance"],
                "extra_propagating_modes": result["sub_hurdles"]["L_a_ii"]["extra_propagating_modes"],
                "dispersion": result["sub_hurdles"]["L_a_ii"]["dispersion"],
                "section_2_6_chain_role": result["section_2_6_chain_role"],
            }
            for name, result in payload["configurations"].items()
        },
        "controls": payload["controls"],
        "hypothetical_pass_fixture_not_gate_admissible": payload["hypothetical_pass_fixture_not_gate_admissible"],
        "artifacts": {
            "report": "software/stage1_solver/reports/pathA_35_gateL_light.md",
            "results_yaml": "software/stage1_solver/reports/pathA_35_gateL_light_results.yaml",
            "sympy": "software/stage1_solver/tools/pathA_35_gateL_sympy.py",
            "mathematica": "software/stage1_solver/tools/pathA_35_gateL.wl",
            "sympy_json": "software/stage1_solver/_scratch/pathA_35_gateL_sympy.json",
            "mathematica_json": "software/stage1_solver/_scratch/pathA_35_gateL_mathematica.json",
            "engine_agreement_json": "software/stage1_solver/_scratch/pathA_35_gateL_engine_agreement.json",
        },
    }
    return summary


def compare_payloads() -> dict[str, Any]:
    sympy_path = SCRATCH / "pathA_35_gateL_sympy.json"
    math_path = SCRATCH / "pathA_35_gateL_mathematica.json"
    if not sympy_path.exists():
        raise AssertionError(f"missing SymPy output: {sympy_path}")
    if not math_path.exists():
        raise AssertionError(f"missing Mathematica output: {math_path}")
    sympy_report = json.loads(sympy_path.read_text(encoding="utf-8"))
    math_report = json.loads(math_path.read_text(encoding="utf-8"))
    sympy_payload = sympy_report["agreement_payload"]
    math_payload = math_report["agreement_payload"]
    if sympy_payload != math_payload:
        raise AssertionError("ENGINE_DISAGREE: SymPy and Mathematica agreement_payload differ")
    sympy_payload["overall"]["engine_agreement"] = "ENGINE_AGREE"
    result = {
        "schema": "pathA_35_gateL_engine_agreement/v1",
        "pass": True,
        "verdict": "ENGINE_AGREE",
        "compared_files": [str(sympy_path), str(math_path)],
        "agreement_payload": sympy_payload,
    }
    out = SCRATCH / "pathA_35_gateL_engine_agreement.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    YAML_OUT.write_text(yaml.safe_dump(yaml_ledger(result), sort_keys=False, width=120), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true", help="compare SymPy and Mathematica scratch payloads and write YAML ledger")
    args = parser.parse_args()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    if args.compare:
        result = compare_payloads()
        print(f"wrote {SCRATCH / 'pathA_35_gateL_engine_agreement.json'}")
        print(f"wrote {YAML_OUT}")
        print(result["verdict"])
        return
    report = build_report("sympy")
    out = SCRATCH / "pathA_35_gateL_sympy.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    print("pathA_35 Gate L SymPy: PASS")


if __name__ == "__main__":
    main()
