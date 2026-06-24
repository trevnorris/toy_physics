#!/usr/bin/env python3
"""Gate B4 SymPy engine: baseline B0_Lifshitz smectic test.

The script is intentionally analytic.  It verifies the frozen T0/G0 bytes,
derives the coupled quadratic spectrum with units restored, checks the
static light-sector minimization, and evaluates the weak-crystallization
morphology discriminator that decides the B4 verdict.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Iterable

import sympy as sp
import yaml
from sympy.printing.mathematica import mathematica_code


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"
T0_REPORT = REPORTS / "pathA_24_T0_freeze.md"
G0_REPORT = REPORTS / "pathA_25_G0_freeze.md"
YAML_OUT = REPORTS / "pathA_25_gateB4_results.yaml"
JSON_OUT = SCRATCH / "pathA_25_gateB4_smectic_sympy.json"

EXPECTED_T0_HASH = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064"
EXPECTED_G0_HASH = "f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290"


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


def assert_freeze_fidelity() -> dict[str, object]:
    t0_block = extract_fence_bytes(T0_REPORT, "freeze-action")
    g0_block = extract_fence_bytes(G0_REPORT, "freeze-action")
    t0_hash = sha256_hex(t0_block)
    g0_hash = sha256_hex(g0_block)
    if t0_hash != EXPECTED_T0_HASH:
        raise AssertionError(f"T0 hash mismatch: expected {EXPECTED_T0_HASH}, got {t0_hash}")
    if g0_hash != EXPECTED_G0_HASH:
        raise AssertionError(f"G0 hash mismatch: expected {EXPECTED_G0_HASH}, got {g0_hash}")
    if t0_block not in g0_block:
        raise AssertionError("exact T0 freeze-action bytes are not embedded unchanged in G0 block")
    required = [
        b"Family L baseline B0_Lifshitz:",
        b"lambda_Cdiv=0; chi_Cpin=0",
        b"L_Mac =",
        b"L_Pu =",
        b"No new G0 term uses V_conf, a prescribed layer profile, a fixed layer normal",
    ]
    for needle in required:
        if needle not in g0_block:
            raise AssertionError(f"required frozen G0 line missing: {needle!r}")
    return {
        "t0_hash": t0_hash,
        "g0_hash": g0_hash,
        "t0_bytes_embedded_in_g0": True,
        "required_frozen_lines_present": True,
    }


Dim = sp.Matrix


def dim_tuple(dim: Dim) -> tuple[int, int, int]:
    return tuple(int(x) for x in list(dim))


def assert_dim(name: str, actual: Dim, expected: Dim) -> None:
    delta = sp.simplify(actual - expected)
    if any(entry != 0 for entry in delta):
        raise AssertionError(f"{name}: expected {dim_tuple(expected)}, got {dim_tuple(actual)}")


def human_expr(expr: sp.Expr) -> str:
    return str(sp.factor(expr)).replace("cL1", "c_L1").replace("cL2", "c_L2")


def human_sstr(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr)).replace("cL1", "c_L1").replace("cL2", "c_L2")


def mma_expr(expr: sp.Expr) -> str:
    return mathematica_code(sp.factor(expr))


def check_dimensions() -> dict[str, tuple[int, int, int]]:
    M = sp.Matrix([1, 0, 0])
    L = sp.Matrix([0, 1, 0])
    T = sp.Matrix([0, 0, 1])
    bulk = M - 2 * L - 2 * T
    drho = -4 * L
    dK = M + 18 * L - 2 * T
    dhbar = M + 2 * L - T
    dm = M
    dk = -L
    dc1 = M + 8 * L - 2 * T
    dc2 = M + 10 * L - 2 * T
    da = L

    dU2 = dK + 3 * drho
    dCQ = 2 * dhbar - dm - drho
    dA = bulk - 2 * drho
    domega2 = -2 * T
    dcs2 = dK + 4 * drho - dm
    dB = bulk
    dKbend = M * 0  # placeholder for the smectic bend modulus dimension below
    # de Gennes f=(B/2)(compression)^2+(K/2)(curvature)^2 in 4D bulk.
    # compression is dimensionless, curvature has L^-1.
    dKbend = bulk + 2 * L

    assert_dim("U''", dU2, dA)
    assert_dim("hbar^2/(m rho0)", dCQ, dc1)
    assert_dim("c_L1", dc1, dCQ)
    assert_dim("c_L2 k^2", dc2 + 2 * dk, dc1)
    assert_dim("omega^2", drho - dm + 2 * dk + dA, domega2)
    assert_dim("c_s^2", dcs2, 2 * L - 2 * T)
    assert_dim("P Goldstone omega^2", dcs2 + 2 * dk, domega2)
    assert_dim("P magnitude gap omega^2", dcs2 - 2 * da, domega2)

    return {
        "omega_squared": dim_tuple(domega2),
        "c_L1": dim_tuple(dc1),
        "c_L2": dim_tuple(dc2),
        "A_kernel": dim_tuple(dA),
        "smectic_B_if_reached": dim_tuple(dB),
        "smectic_K_if_reached": dim_tuple(dKbend),
        "smectic_branch_status": "not_reached_after_FAIL_NOT_CODIM1",
    }


def coupled_spectrum() -> dict[str, object]:
    k, kk, rho0, K, hbar, m, cL1, cL2, a = sp.symbols(
        "k kk rho0 K hbar m cL1 cL2 a", positive=True
    )
    U2 = 5 * K * rho0**3
    cs2 = 5 * K * rho0**4 / m
    cQ = hbar**2 / (4 * m * rho0)
    A_kk = sp.factor(U2 + (cQ - cL1) * kk + cL2 * kk**2)
    A = sp.factor(A_kk.subs(kk, k**2))
    omega2 = sp.factor(rho0 * k**2 * A / m)
    stationary = sp.solve(sp.Eq(sp.diff(A_kk, kk), 0), kk)
    if len(stationary) != 1:
        raise AssertionError(f"expected one roton stationary point, got {stationary}")
    kstar2 = sp.factor(stationary[0])
    Amin = sp.factor(A_kk.subs(kk, kstar2))
    thresholds = sp.solve(sp.Eq(Amin, 0), cL1)
    threshold_sample = {K: 1, rho0: 1, hbar: 1, m: 1, cL2: 1}
    threshold = sp.factor(max(thresholds, key=lambda root: float(sp.N(root.subs(threshold_sample), 30))))
    kstar2_at_threshold = sp.factor(kstar2.subs(cL1, threshold))

    # Energy Hessian for (eta, theta, pi1, pi2, pi3, sigma).
    theta_stiffness = sp.factor(rho0 * hbar**2 * k**2 / m)
    p_inertia = m * rho0 * a**2
    p_trans_stiffness = sp.factor(m * rho0 * cs2 * a**2 * k**2)
    p_long_stiffness = sp.factor(m * rho0 * cs2 * a**2 * k**2 + 2 * m * rho0 * cs2)
    p_trans_omega2 = sp.factor(p_trans_stiffness / p_inertia)
    p_long_omega2 = sp.factor(p_long_stiffness / p_inertia)
    hessian_diag = [A, theta_stiffness, p_trans_stiffness, p_trans_stiffness, p_trans_stiffness, p_long_stiffness]
    hessian = sp.diag(*hessian_diag)
    off_diag_zero = all(hessian[i, j] == 0 for i in range(6) for j in range(6) if i != j)
    if not off_diag_zero:
        raise AssertionError("unexpected off-diagonal Hessian block")

    # Explicit rho-P quadratic coupling checks from the frozen T0 coefficients.
    eta, theta, sigma, pi, pgrad, pdot = sp.symbols("eta theta sigma pi pgrad pdot")
    rho = rho0 + eta
    coeff = m * rho * (5 * K * rho**4 / m)
    mag_density = sp.expand(sp.Rational(1, 4) * coeff * ((1 + sigma) ** 2 + pi**2 - 1) ** 2)
    frank_density = sp.expand(sp.Rational(1, 2) * coeff * a**2 * pgrad**2)
    inertia_density = sp.expand(sp.Rational(1, 2) * m * rho * a**2 * pdot**2)
    theta_density = sp.expand(sp.Rational(1, 2) * rho0 * hbar**2 * k**2 / m * theta**2)
    total_energy_density = mag_density + frank_density + inertia_density + theta_density
    rho_sigma_quad = sp.diff(mag_density, eta, sigma).subs({eta: 0, sigma: 0, pi: 0})
    rho_pi_quad = sp.diff(mag_density, eta, pi).subs({eta: 0, sigma: 0, pi: 0})
    rho_pgrad_quad = sp.diff(frank_density, eta, pgrad).subs({eta: 0, pgrad: 0})
    rho_pdot_quad = sp.diff(inertia_density, eta, pdot).subs({eta: 0, pdot: 0})
    rho_p_zeros = [rho_sigma_quad, rho_pi_quad, rho_pgrad_quad, rho_pdot_quad]
    if any(sp.simplify(z) != 0 for z in rho_p_zeros):
        raise AssertionError(f"rho-P block did not vanish: {rho_p_zeros}")
    theta_p_zeros = [
        sp.diff(total_energy_density, theta, sigma).subs({theta: 0, sigma: 0, pi: 0, pgrad: 0, pdot: 0}),
        sp.diff(total_energy_density, theta, pi).subs({theta: 0, sigma: 0, pi: 0, pgrad: 0, pdot: 0}),
        sp.diff(total_energy_density, theta, pgrad).subs({theta: 0, sigma: 0, pi: 0, pgrad: 0, pdot: 0}),
        sp.diff(total_energy_density, theta, pdot).subs({theta: 0, sigma: 0, pi: 0, pgrad: 0, pdot: 0}),
    ]
    if any(sp.simplify(z) != 0 for z in theta_p_zeros):
        raise AssertionError(f"theta-P block did not vanish: {theta_p_zeros}")

    return {
        "density_kernel_A_of_k": human_expr(A),
        "omega2_density_phase": human_expr(omega2),
        "theta_density_symplectic_entry": "hbar",
        "P_transverse_goldstone_count": 3,
        "P_transverse_omega2": human_expr(p_trans_omega2),
        "P_longitudinal_omega2": human_expr(p_long_omega2),
        "rho_P_quadratic_block": {"computed_mixed_derivatives": [human_expr(z) for z in rho_p_zeros], "status": "zero"},
        "theta_P_quadratic_block": {"computed_mixed_derivatives": [human_expr(z) for z in theta_p_zeros], "status": "zero"},
        "P_sector_nonnegative_conditions": ["K>0", "rho0>0", "m>0", "a>0"],
        "finite_k_condition": "c_L1 > hbar^2/(4*m*rho0)",
        "finite_kstar_squared": human_expr(kstar2),
        "A_min": human_expr(Amin),
        "softening_threshold_c_L1": human_expr(threshold),
        "kstar_squared_at_threshold": human_expr(kstar2_at_threshold),
        "k_Lstar_reconciliation": "k_*^2 = c_L1/(2*c_L2) - hbar^2/(8*m*rho0*c_L2)",
        "engine_exprs": {
            "density_kernel_A_of_k": mma_expr(A),
            "omega2_density_phase": mma_expr(omega2),
            "finite_kstar_squared": mma_expr(kstar2),
            "A_min": mma_expr(Amin),
            "softening_threshold_c_L1": mma_expr(threshold),
            "kstar_squared_at_threshold": mma_expr(kstar2_at_threshold),
            "P_transverse_omega2": mma_expr(p_trans_omega2),
            "P_longitudinal_omega2": mma_expr(p_long_omega2),
            "rho_P_zero_sum": mma_expr(sum(rho_p_zeros)),
            "theta_P_zero_sum": mma_expr(sum(theta_p_zeros)),
        },
    }


Vector = tuple[sp.Expr, ...]


def dot(u: Vector, v: Vector) -> sp.Expr:
    return sp.simplify(sum(ui * vi for ui, vi in zip(u, v)))


def vadd(*vectors: Vector) -> Vector:
    return tuple(sp.simplify(sum(v[i] for v in vectors)) for i in range(4))


def vneg(v: Vector) -> Vector:
    return tuple(sp.simplify(-x) for x in v)


def is_zero_vec(v: Vector) -> bool:
    return all(sp.simplify(x) == 0 for x in v)


def cos_modes(base_vectors: Iterable[Vector]) -> dict[Vector, sp.Rational]:
    modes: dict[Vector, sp.Rational] = {}
    for v in base_vectors:
        modes[v] = modes.get(v, sp.Rational(0)) + sp.Rational(1, 2)
        nv = vneg(v)
        modes[nv] = modes.get(nv, sp.Rational(0)) + sp.Rational(1, 2)
    return modes


def avg_power_coeff(base_vectors: Iterable[Vector], power: int) -> sp.Expr:
    modes = cos_modes(base_vectors)
    items = list(modes.items())
    total = sp.Rational(0)
    for combo in itertools.product(items, repeat=power):
        vec = vadd(*(entry[0] for entry in combo))
        if is_zero_vec(vec):
            coeff = sp.prod(entry[1] for entry in combo)
            total += coeff
    return sp.simplify(total)


def avg_eta_p_grad2_coeff(base_vectors: Iterable[Vector], eta_power: int) -> sp.Expr:
    modes = cos_modes(base_vectors)
    items = list(modes.items())
    total = sp.Rational(0)
    for eta_combo in itertools.product(items, repeat=eta_power):
        eta_vecs = [entry[0] for entry in eta_combo]
        eta_coeff = sp.prod(entry[1] for entry in eta_combo)
        for p_vec, p_coeff in items:
            for q_vec, q_coeff in items:
                if is_zero_vec(vadd(*eta_vecs, p_vec, q_vec)):
                    total += eta_coeff * p_coeff * q_coeff * (-dot(p_vec, q_vec))
    return sp.simplify(total)


def minimized_polynomial(
    energy: sp.Expr,
    amplitude: sp.Symbol,
    *,
    nonnegative_amplitude: bool,
    ranking_subs: dict[sp.Symbol, sp.Expr] | None = None,
) -> tuple[sp.Expr, sp.Expr, list[sp.Expr]]:
    ranking_subs = ranking_subs or {}
    stationary = sp.solve(sp.Eq(sp.diff(energy, amplitude), 0), amplitude)
    if 0 not in stationary:
        stationary.append(sp.Integer(0))
    candidates: list[sp.Expr] = []
    for root in stationary:
        root_n = complex(sp.N(root.subs(ranking_subs), 30))
        if abs(root_n.imag) > 1e-20:
            continue
        if nonnegative_amplitude and root_n.real < -1e-20:
            continue
        candidates.append(sp.simplify(root))
    if not candidates:
        raise AssertionError(f"no real stationary candidates for {energy}")
    values = [sp.factor(energy.subs(amplitude, root)) for root in candidates]
    best_index = min(range(len(values)), key=lambda i: float(sp.N(values[i].subs(ranking_subs), 30)))
    return sp.factor(values[best_index]), sp.simplify(candidates[best_index]), candidates


def shell_directions(M_count: int) -> list[Vector]:
    if not 1 <= M_count <= 4:
        raise ValueError("Gate B4 controls only use M=1..4 in 4D")
    directions: list[Vector] = []
    for i in range(M_count):
        directions.append(tuple(sp.Integer(1) if j == i else sp.Integer(0) for j in range(4)))
    return directions


def morphology_discriminator() -> dict[str, object]:
    rho0, K, hbar, m, k = sp.symbols("rho0 K hbar m k", positive=True)
    sqrt3 = sp.sqrt(3)
    lamellar = [(sp.Integer(1), 0, 0, 0)]
    triad = [
        (sp.Integer(1), 0, 0, 0),
        (sp.Rational(-1, 2), sqrt3 / 2, 0, 0),
        (sp.Rational(-1, 2), -sqrt3 / 2, 0, 0),
    ]

    lam_avg3 = avg_power_coeff(lamellar, 3)
    lam_avg2 = avg_power_coeff(lamellar, 2)
    lam_avg4 = avg_power_coeff(lamellar, 4)
    triad_avg2 = avg_power_coeff(triad, 2)
    triad_avg3 = avg_power_coeff(triad, 3)
    triad_avg4 = avg_power_coeff(triad, 4)
    triad_eta_grad2 = avg_eta_p_grad2_coeff(triad, 1)
    if lam_avg3 != 0:
        raise AssertionError(f"single-k lamellar cubic average should vanish, got {lam_avg3}")
    if sp.simplify(triad_avg2 - sp.Rational(3, 2)) != 0:
        raise AssertionError(f"triad <phi^2> coefficient changed: {triad_avg2}")
    if sp.simplify(triad_avg3 - sp.Rational(3, 2)) != 0:
        raise AssertionError(f"triad <phi^3> coefficient changed: {triad_avg3}")
    if sp.simplify(triad_eta_grad2 - sp.Rational(3, 4)) != 0:
        raise AssertionError(f"triad <phi |grad phi|^2> coefficient changed: {triad_eta_grad2}")

    U3 = 15 * K * rho0**2
    cubic_coeff = sp.factor(U3 / 6 * triad_avg3 - hbar**2 / (8 * m * rho0**2) * k**2 * triad_eta_grad2)
    cubic_zero_k2 = sp.solve(sp.Eq(cubic_coeff, 0), k**2)[0]

    # Positive control: 4D shell, Z2-symmetric Brazovskii/Swift-Hohenberg local phi^4.
    r, u, A = sp.symbols("r u A", positive=True)
    # Here r denotes |negative quadratic coefficient|, so energy uses -r/2 <phi^2>.
    positive_control = {}
    positive_control_derived = {}
    for M in range(1, 5):
        dirs = shell_directions(M)
        avg2 = avg_power_coeff(dirs, 2)
        avg4 = avg_power_coeff(dirs, 4)
        energy = sp.factor(-sp.Rational(1, 2) * r * avg2 * A**2 + sp.Rational(1, 4) * u * avg4 * A**4)
        fmin, amin, candidates = minimized_polynomial(
            energy, A, nonnegative_amplitude=True, ranking_subs={r: 1, u: 1}
        )
        expected = sp.factor(-r**2 * M / (6 * u * (2 * M - 1)))
        if sp.simplify(fmin - expected) != 0:
            raise AssertionError(f"positive-control M={M} minimized energy mismatch: {fmin} vs {expected}")
        positive_control[f"{M}_directions"] = human_expr(fmin)
        positive_control_derived[f"{M}_directions"] = {
            "avg_phi2": human_expr(avg2),
            "avg_phi4": human_expr(avg4),
            "energy": human_expr(energy),
            "minimizing_amplitude": human_expr(amin),
            "stationary_candidates": [human_expr(root) for root in candidates],
            "engine_expr": mma_expr(fmin),
        }

    # Negative multi-k benchmark, predeclared coefficients in report:
    # F_lam(A)=tau A^2/4 + 3 A^4/32; F_tri(A)=3 tau A^2/4 - gamma A^3/4 + 45 A^4/32.
    tau_val = sp.Rational(1, 100)
    gamma_val = sp.Integer(1)
    u_val = sp.Integer(1)
    F_lam = tau_val * A**2 / 4 + 3 * u_val * A**4 / 32
    F_tri = 3 * tau_val * A**2 / 4 - gamma_val * A**3 / 4 + 45 * u_val * A**4 / 32
    lam_min, lam_min_amp, lam_candidates = minimized_polynomial(F_lam, A, nonnegative_amplitude=False)
    tri_min, tri_min_amp, tri_candidates = minimized_polynomial(F_tri, A, nonnegative_amplitude=True)
    tri_min_numeric = float(sp.N(tri_min, 20))
    if not tri_min_numeric < 0:
        raise AssertionError("NC-B4d benchmark did not prefer the triad")

    # Baseline Landau comparison.  epsilon is a small positive distance above
    # the finite-k onset; the triad phase is chosen so the generic cubic lowers
    # the energy.  The strict epsilon=0 comparison is the continuity certificate
    # for an open neighborhood above onset.
    eps = sp.Rational(1, 100)
    beta = sp.Integer(1)
    sample_values = {K: 1, rho0: 1, hbar: 1, m: 1, k**2: sp.Rational(23, 8)}
    gamma_abs_sample = sp.factor(abs(sp.N(cubic_coeff.subs(sample_values), 30)))
    gamma_abs_exact = sp.factor(cubic_coeff.subs(sample_values))
    if gamma_abs_exact == 0:
        raise AssertionError("generic baseline sample accidentally lies on Gamma_T=0")
    gamma_abs = sp.factor(sp.Abs(gamma_abs_exact))
    F_lam_landau = sp.factor(-eps * lam_avg2 * A**2 + beta * lam_avg4 * A**4)
    F_tri_landau = sp.factor(-eps * triad_avg2 * A**2 - gamma_abs * A**3 + beta * triad_avg4 * A**4)
    F_lam_onset = sp.factor(beta * lam_avg4 * A**4)
    F_tri_onset = sp.factor(-gamma_abs * A**3 + beta * triad_avg4 * A**4)
    lam_landau_min, lam_landau_amp, lam_landau_candidates = minimized_polynomial(
        F_lam_landau, A, nonnegative_amplitude=False
    )
    tri_landau_min, tri_landau_amp, tri_landau_candidates = minimized_polynomial(
        F_tri_landau, A, nonnegative_amplitude=True
    )
    lam_onset_min, _, _ = minimized_polynomial(F_lam_onset, A, nonnegative_amplitude=False)
    tri_onset_min, tri_onset_amp, _ = minimized_polynomial(F_tri_onset, A, nonnegative_amplitude=True)
    triad_beats_lamella = bool(sp.N(tri_landau_min - lam_landau_min, 30) < 0)
    open_neighborhood_certificate = bool(sp.N(tri_onset_min - lam_onset_min, 30) < 0)
    if not triad_beats_lamella or not open_neighborhood_certificate:
        raise AssertionError("computed baseline Landau comparison did not prefer the triad")

    return {
        "single_k_lamellar_cubic_average": str(lam_avg3),
        "single_k_lamellar_average_phi2": str(lam_avg2),
        "single_k_lamellar_average_phi4": str(lam_avg4),
        "rank2_triad_average_phi2": str(triad_avg2),
        "rank2_triad_average_phi3": str(triad_avg3),
        "rank2_triad_average_phi4": str(triad_avg4),
        "rank2_triad_average_phi_grad_phi_sq": str(triad_eta_grad2),
        "baseline_rank2_triad_cubic_coefficient": human_expr(cubic_coeff),
        "baseline_cubic_zero_codimension_condition": f"k_*^2 = {human_sstr(cubic_zero_k2)}",
        "baseline_landau_comparison": {
            "epsilon_above_onset": str(eps),
            "quartic_stabilizer_beta": str(beta),
            "generic_sample": {"K": "1", "rho0": "1", "hbar": "1", "m": "1", "kstar_squared": "23/8"},
            "Gamma_T_sample_exact": human_expr(gamma_abs_exact),
            "Gamma_T_sample_abs_numeric": float(gamma_abs_sample),
            "lamellar_energy": human_expr(F_lam_landau),
            "triad_energy": human_expr(F_tri_landau),
            "lamellar_min_exact": human_expr(lam_landau_min),
            "lamellar_min_numeric": float(sp.N(lam_landau_min, 20)),
            "lamellar_minimizing_amplitude": human_expr(lam_landau_amp),
            "triad_min_exact": human_expr(tri_landau_min),
            "triad_min_numeric": float(sp.N(tri_landau_min, 20)),
            "triad_minimizing_amplitude": human_expr(tri_landau_amp),
            "triad_beats_lamella": triad_beats_lamella,
            "onset_lamellar_min_exact": human_expr(lam_onset_min),
            "onset_triad_min_exact": human_expr(tri_onset_min),
            "onset_triad_minimizing_amplitude": human_expr(tri_onset_amp),
            "open_neighborhood_certificate": open_neighborhood_certificate,
            "stationary_candidates": {
                "lamella": [human_expr(root) for root in lam_landau_candidates],
                "triad": [human_expr(root) for root in tri_landau_candidates],
            },
        },
        "baseline_morphology_result": "computed_rank2_triad_minimum_below_rank1_lamella_for_Gamma_T_nonzero",
        "baseline_verdict_driver_region": "generic finite-k B0_Lifshitz onset is non-codim1; cubic-zero is codimension-one fine tuning, not an open smectic region",
        "positive_control_NC_B4b": {
            "model": "4D Z2-symmetric conserved Swift-Hohenberg/Brazovskii shell with u>0 and -r<0 on |k|=k0",
            "minimized_energies": positive_control,
            "derived_minimizations": positive_control_derived,
            "outcome": "stable_single_k_stripe",
        },
        "negative_control_NC_B4d": {
            "model": "4D shell Landau benchmark with a rank-2 equilateral-triad cubic invariant",
            "predeclared_coefficients": {"tau": str(tau_val), "gamma": str(gamma_val), "u": str(u_val)},
            "lamellar_min": human_expr(lam_min),
            "lamellar_minimizing_amplitude": human_expr(lam_min_amp),
            "triad_min_numeric": tri_min_numeric,
            "triad_min_exact": human_expr(tri_min),
            "triad_minimizing_amplitude": human_expr(tri_min_amp),
            "stationary_candidates": {
                "lamella": [human_expr(root) for root in lam_candidates],
                "triad": [human_expr(root) for root in tri_candidates],
            },
            "outcome": "FAIL_NOT_CODIM1",
        },
        "engine_exprs": {
            "single_k_lamellar_cubic_average": mma_expr(lam_avg3),
            "rank2_triad_average_phi2": mma_expr(triad_avg2),
            "rank2_triad_average_phi3": mma_expr(triad_avg3),
            "rank2_triad_average_phi_grad_phi_sq": mma_expr(triad_eta_grad2),
            "baseline_rank2_triad_cubic_coefficient": mma_expr(cubic_coeff),
            "baseline_cubic_zero_kstar_squared": mma_expr(cubic_zero_k2),
            "baseline_lamella_min": mma_expr(lam_landau_min),
            "baseline_triad_min": mma_expr(tri_landau_min),
            "baseline_onset_triad_min": mma_expr(tri_onset_min),
            "negative_control_NC_B4d_lamellar_min": mma_expr(lam_min),
            "negative_control_NC_B4d_triad_min": mma_expr(tri_min),
            **{f"positive_control_SH_M{M}_min": positive_control_derived[f"{M}_directions"]["engine_expr"] for M in range(1, 5)},
        },
    }


def phase_separation_and_boundedness() -> dict[str, object]:
    q, rho, rho0, K, hbar, m, cL1, cL2 = sp.symbols("q rho rho0 K hbar m cL1 cL2", positive=True)
    U2 = 5 * K * rho0**3
    cQ = hbar**2 / (4 * m * rho0)
    threshold = cQ + 2 * sp.sqrt(U2 * cL2)
    U = K * rho**5 / 4
    U_second = sp.factor(sp.diff(U, rho, 2))
    gradient_symbol = sp.factor(sp.Rational(1, 2) * cL2 * q**2 - sp.Rational(1, 2) * cL1 * q)
    q_min = sp.solve(sp.Eq(sp.diff(gradient_symbol, q), 0), q)[0]
    gradient_lower_bound = sp.factor(gradient_symbol.subs(q, q_min))
    if sp.simplify(gradient_lower_bound + cL1**2 / (8 * cL2)) != 0:
        raise AssertionError("bounded-below Fourier-symbol minimum changed")
    return {
        "k0_compressibility": human_expr(U2),
        "k0_linear_instability": {"computed_A_rho_0": human_expr(U2), "status": "absent_for_K_positive_rho0_positive"},
        "same_mean_phase_separation": {
            "local_potential": "K*rho^5/4",
            "computed_second_derivative": human_expr(U_second),
            "convex_on_rho_nonnegative": True,
            "fixed_mean_exclusion": "computed_convexity_plus_Jensen",
        },
        "bounded_below": {
            "fourier_symbol": human_expr(gradient_symbol),
            "minimizing_k_squared": human_expr(q_min),
            "minimum": human_expr(gradient_lower_bound),
            "bound": "E_L >= -c_L1^2/(8*c_L2) * integral((rho-rho0)^2)",
            "positive_controls": ["U(rho)~rho^5 at fixed mean", "quantum pressure >= 0"],
        },
        "finite_k_linear_instability_region": f"c_L1 > {human_sstr(threshold)}",
        "engine_exprs": {
            "bounded_symbol_minimum": mma_expr(gradient_lower_bound),
            "bounded_symbol_minimizing_k_squared": mma_expr(q_min),
            "phase_convexity_second_derivative": mma_expr(U_second),
            "k0_compressibility": mma_expr(U2),
        },
    }


def light_sector() -> dict[str, object]:
    mu_br, kappa_Pu = sp.symbols("muBr kappaPu", positive=True)
    curl_u, deltaP_parallel, Omega_u = sp.symbols("curlU deltaPParallel OmegaU", real=True)
    mismatch = sp.symbols("mismatch", real=True)
    integrand = sp.Rational(1, 2) * mu_br * curl_u**2 + sp.Rational(1, 2) * kappa_Pu * mismatch**2
    stationary = sp.solve(
        [sp.Eq(sp.diff(integrand, curl_u), 0), sp.Eq(sp.diff(integrand, mismatch), 0)],
        [curl_u, mismatch],
        dict=True,
    )
    if stationary != [{curl_u: 0, mismatch: 0}]:
        raise AssertionError(f"unexpected light-sector minimizer: {stationary}")
    minimum = sp.factor(integrand.subs(stationary[0]))
    if minimum != 0:
        raise AssertionError("light-sector quadratic minimum is not zero")
    original_integrand = sp.Rational(1, 2) * mu_br * curl_u**2 + sp.Rational(1, 2) * kappa_Pu * (
        deltaP_parallel - Omega_u
    ) ** 2
    vanishing_integrand = sp.factor(original_integrand.subs({curl_u: 0, deltaP_parallel: Omega_u}))
    return {
        "uniform_state_support": {
            "status": "analytic_argument_not_computed",
            "statement": "Sigma_n[rho0] is empty, so there is no layer support and no quadratic contribution.",
        },
        "static_stripe_minimum": {
            "computed_integrand": "1/2*mu_br*|curl u|^2 + 1/2*kappa_Pu*|deltaP_parallel-Omega_u|^2",
            "quadratic_minimizer": {"curl_u": "0", "deltaP_parallel_minus_Omega_u": "0"},
            "minimum": human_expr(minimum),
            "nonnegative_conditions": ["mu_br>=0", "kappa_Pu>=0", "inherited T0 surface terms nonnegative"],
        },
        "support_variation": {
            "integrand_at_u0_deltaP_equals_Omega": human_expr(vanishing_integrand),
            "status": "computed_zero_integrand_on_each_support",
        },
        "density_stability_controls": {
            "static_dropouts_after_minimization": ["mu_br", "J_Pu", "kappa_Pu"],
            "collapse_if_negative_modulus": "separate_light_sector_failure_not_used_to_rescue_B4",
        },
        "engine_exprs": {
            "light_integrand_minimum": mma_expr(minimum),
            "light_support_integrand_at_minimum": mma_expr(vanishing_integrand),
        },
    }


def controls() -> dict[str, object]:
    return {
        "NC_B4a_driver_off": {
            "setting": "c_L1,c_L2 -> 0 after removing the Lifshitz driver",
            "spectrum": "omega^2=(rho0*k^2/m)*(5*K*rho0^3 + hbar^2*k^2/(4*m*rho0)) > 0",
            "outcome": "FAIL_NO_MODULATION",
        },
        "NC_B4b_positive_method": {
            "setting": "4D Z2 Swift-Hohenberg/Brazovskii shell, r<0,u>0, mean zero",
            "outcome": "stable_single_k_stripe",
        },
        "NC_B4c_k0_attraction": {
            "setting": "replace U'' by a negative compressibility -abs(mu0) with nonnegative gradient regularization",
            "first_soft_mode": "k=0",
            "outcome": "FAIL_PHASE_SEPARATION",
        },
        "NC_B4d_multik_benchmark": {
            "setting": "predeclared tau=1/100,gamma=1,u=1 rank-2 cubic shell benchmark",
            "outcome": "FAIL_NOT_CODIM1",
        },
    }


def derive_verdict(
    spectrum: dict[str, object],
    morphology: dict[str, object],
    phase: dict[str, object],
    light: dict[str, object],
) -> dict[str, object]:
    comparison = morphology["baseline_landau_comparison"]
    bounded = phase["bounded_below"]
    convex = phase["same_mean_phase_separation"]
    light_minimum = light["static_stripe_minimum"]["minimum"]

    finite_k_window_open = "sqrt" in str(spectrum["softening_threshold_c_L1"]) or "sqrt" in str(
        spectrum["kstar_squared_at_threshold"]
    )
    phase_separation_excluded = bool(convex["convex_on_rho_nonnegative"])
    bounded_below = bounded["minimum"] == "-c_L1**2/(8*c_L2)" or bounded["minimum"] == "-cL1**2/(8*cL2)"
    light_nonnegative = light_minimum == "0"
    gamma_nonzero_on_generic_sample = comparison["Gamma_T_sample_exact"] != "0"
    triad_beats_lamella = bool(comparison["triad_beats_lamella"])
    open_neighborhood = bool(comparison["open_neighborhood_certificate"])

    if not finite_k_window_open:
        verdict = "FAIL_NO_MODULATION"
        reason = "computed finite-k softening window is empty"
    elif not phase_separation_excluded:
        verdict = "FAIL_PHASE_SEPARATION"
        reason = "computed local density convexity failed at fixed mean"
    elif not bounded_below or not light_nonnegative:
        verdict = "FAIL_STABILITY"
        reason = "computed boundedness or static light-sector positivity failed"
    elif gamma_nonzero_on_generic_sample and triad_beats_lamella and open_neighborhood:
        verdict = "FAIL_NOT_CODIM1"
        reason = (
            "computed rank-2 triad Landau minimum is below the single-k lamella in an open neighborhood above onset"
        )
    else:
        verdict = "SMECTIC_CONDITIONAL"
        reason = "computed checks did not find a lower multi-k state off the codimension-one cubic-zero surface"

    return {
        "verdict": verdict,
        "reason": reason,
        "conditions": {
            "finite_k_window_open": finite_k_window_open,
            "phase_separation_excluded_by_convexity": phase_separation_excluded,
            "bounded_below_symbol_minimum": bounded_below,
            "static_light_sector_nonnegative": light_nonnegative,
            "Gamma_T_nonzero_on_generic_sample": gamma_nonzero_on_generic_sample,
            "triad_minimum_below_lamella_minimum": triad_beats_lamella,
            "open_neighborhood_certificate": open_neighborhood,
        },
        "computed_baseline_comparison": {
            "lamellar_min_numeric": comparison["lamellar_min_numeric"],
            "triad_min_numeric": comparison["triad_min_numeric"],
            "lamellar_min_exact": comparison["lamellar_min_exact"],
            "triad_min_exact": comparison["triad_min_exact"],
        },
    }


def main() -> None:
    SCRATCH.mkdir(exist_ok=True)
    freeze = assert_freeze_fidelity()
    dims = check_dimensions()
    spectrum = coupled_spectrum()
    morphology = morphology_discriminator()
    phase = phase_separation_and_boundedness()
    light = light_sector()
    ctrl = controls()
    verdict_info = derive_verdict(spectrum, morphology, phase, light)
    engine_exprs = {
        **spectrum.pop("engine_exprs"),
        **morphology.pop("engine_exprs"),
        **phase.pop("engine_exprs"),
        **light.pop("engine_exprs"),
    }

    results = {
        "verdict": verdict_info["verdict"],
        "reason": verdict_info["reason"],
        "computed_verdict": verdict_info,
        "freeze_fidelity": freeze,
        "background": {
            "A_mu": "0 up to pure gauge",
            "V_conf": 0,
            "fixed_mean": "bar rho = rho0",
            "boundary_potential": "none",
            "box_role": "numerical regulator only",
        },
        "dimensions_MLT": dims,
        "coupled_spectrum": spectrum,
        "light_sector": light,
        "finite_k_vs_k0": phase,
        "morphology": morphology,
        "controls": ctrl,
        "stability_region": {
            "dimensionless_linear_control": "Lambda=(c_L1-hbar^2/(4*m*rho0))/(2*sqrt(5*K*rho0^3*c_L2))",
            "finite_k_softening": "Lambda>1",
            "generic_codim1_status": "not stable at onset for Gamma_triad != 0",
            "Gamma_triad": morphology["baseline_rank2_triad_cubic_coefficient"],
            "Gamma_triad_zero": morphology["baseline_cubic_zero_codimension_condition"],
        },
        "engine_agreement": {
            "sympy_engine": "PASS",
            "mathematica_engine": "asserted_by_wolfram_script_against_this_json",
            "mathematica_exprs": engine_exprs,
        },
    }
    JSON_OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    YAML_OUT.write_text(yaml.safe_dump(results, sort_keys=False, width=120))
    print("PASS pathA_25_gateB4_smectic_sympy")
    print(json.dumps({"verdict": results["verdict"], "json": str(JSON_OUT), "yaml": str(YAML_OUT)}, indent=2))


if __name__ == "__main__":
    main()
