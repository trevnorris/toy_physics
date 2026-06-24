#!/usr/bin/env python3
"""pathA_24 T1 SymPy/numpy cross-check for the frozen polar-arrow wall.

This script is intentionally target-blind: it transcribes only the frozen T0
polar action, checks the freeze-action hash, derives the static wall
functional with units restored, and tests whether the O(4)-isotropic baseline
has a localized stable wall.  Negative-control spectra are included so the
stability gate can return both stable and unstable outcomes.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import sympy as sp


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
PROJECT_ROOT = STAGE1_ROOT.parents[1]
SCRATCH = STAGE1_ROOT / "_scratch"
FREEZE_REPORT = STAGE1_ROOT / "reports" / "pathA_24_T0_freeze.md"
PDE_TEX = PROJECT_ROOT / "research" / "pde" / "paper" / "pde.tex"
CONCEPTUAL = PROJECT_ROOT / "docs" / "conceptual_foundation.md"
PATHA23_STAGE0 = STAGE1_ROOT / "reports" / "pathA_23_stage0_action_and_contracts.md"

EXPECTED_FREEZE_SHA256 = (
    "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064"
)
EXPECTED_L_POL_LINES = [
    "  L_pol =",
    "      (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)",
    "    - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)",
    "    - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2.",
]

SPECTRUM_TOL = 1.0e-3
TAU_THRESHOLD = "1000 * tau_Hubble"


def extract_freeze_block(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"^```freeze-action\n(?P<block>.*?)\n```", text, re.M | re.S)
    if not match:
        raise AssertionError(f"freeze-action block not found in {path}")
    return match.group("block")


def assert_freeze_fidelity() -> dict:
    block = extract_freeze_block(FREEZE_REPORT)
    digest = hashlib.sha256((block + "\n").encode("utf-8")).hexdigest()
    if digest != EXPECTED_FREEZE_SHA256:
        raise AssertionError(
            f"freeze-action SHA-256 mismatch: expected {EXPECTED_FREEZE_SHA256}, got {digest}"
        )
    line_positions: list[int] = []
    lines = block.splitlines()
    for expected in EXPECTED_L_POL_LINES:
        if expected not in lines:
            raise AssertionError(f"frozen L_pol line missing or changed: {expected!r}")
        line_positions.append(lines.index(expected))
    if line_positions != sorted(line_positions):
        raise AssertionError(f"frozen L_pol lines are out of order: {line_positions}")

    forbidden_action_fragments = [
        "gamma_w",
        "w_hat",
        "P dot w",
        "delta(w)",
        "Z(w)|nabla P|",
        "A_i P^i",
    ]
    forbidden_present = [frag for frag in forbidden_action_fragments if frag in block]
    if forbidden_present:
        raise AssertionError(f"forbidden frozen-action fragments present: {forbidden_present}")
    if "no explicit easy axis" not in block:
        raise AssertionError("baseline no-explicit-easy-axis line missing from freeze-action")

    return {
        "path": str(FREEZE_REPORT),
        "sha256": digest,
        "checked_lines": EXPECTED_L_POL_LINES,
        "forbidden_fragments_absent": forbidden_action_fragments,
    }


Dim = sp.Matrix


def dim_tuple(dim: Dim) -> tuple[int, int, int]:
    return tuple(int(x) for x in list(dim))


def dim_str(dim: Dim) -> str:
    labels = ("M", "L", "T")
    parts: list[str] = []
    for label, power in zip(labels, dim_tuple(dim)):
        if power == 0:
            continue
        parts.append(label if power == 1 else f"{label}^{power}")
    return " ".join(parts) if parts else "1"


def assert_dim(name: str, actual: Dim, expected: Dim) -> None:
    delta = sp.simplify(actual - expected)
    if any(entry != 0 for entry in delta):
        raise AssertionError(
            f"{name}: expected {dim_tuple(expected)} ({dim_str(expected)}), "
            f"got {dim_tuple(actual)} ({dim_str(actual)})"
        )


def check_dimensions() -> dict:
    M = sp.Matrix([1, 0, 0])
    L = sp.Matrix([0, 1, 0])
    T = sp.Matrix([0, 0, 1])
    d_m = M
    d_rho = -4 * L
    d_K = M + 18 * L - 2 * T
    d_hbar = M + 2 * L - T
    d_a = L
    d_w = L
    d_rhop = d_rho - L
    d_gradP = -L
    dimless = sp.Matrix([0, 0, 0])

    d_cs2 = d_K + 4 * d_rho - d_m
    assert_dim("c_s^2(rho)", d_cs2, 2 * L - 2 * T)

    qp = 2 * d_hbar - d_m - d_rho + 2 * d_rhop
    u = d_K + 5 * d_rho
    op_grad = d_m + d_rho + d_cs2 + 2 * d_a + 2 * d_gradP
    op_pot = d_m + d_rho + d_cs2 + dimless
    lag_density = M - 2 * L - 2 * T
    sigma_dim = lag_density + d_w
    energy_per_3area = M - L - 2 * T
    k_p = d_m + d_rho + d_cs2 + 2 * d_a
    b_depth = d_m + d_rho + d_cs2

    for name, actual in {
        "GNLS quantum pressure": qp,
        "GNLS U(rho)": u,
        "OP gradient": op_grad,
        "OP potential": op_pot,
    }.items():
        assert_dim(name, actual, lag_density)
    assert_dim("surface tension", sigma_dim, energy_per_3area)
    assert_dim("spread sigma K_P/L", k_p - L, energy_per_3area)
    assert_dim("radial saddle sigma B*a", b_depth + d_a, energy_per_3area)
    assert_dim("omega^2", d_cs2 - 2 * d_a, -2 * T)

    return {
        "energy_density": {"triple_MLT": dim_tuple(lag_density), "string": dim_str(lag_density)},
        "sigma_brane": {"triple_MLT": dim_tuple(sigma_dim), "string": dim_str(sigma_dim)},
        "terms": {
            "gnls_quantum_pressure": {"triple_MLT": dim_tuple(qp), "string": dim_str(qp)},
            "gnls_U": {"triple_MLT": dim_tuple(u), "string": dim_str(u)},
            "op_gradient": {"triple_MLT": dim_tuple(op_grad), "string": dim_str(op_grad)},
            "op_potential": {"triple_MLT": dim_tuple(op_pot), "string": dim_str(op_pot)},
        },
    }


def finite_difference_spectrum(
    potential, x_min: float, x_max: float, n: int, count: int
) -> list[float]:
    h = (x_max - x_min) / (n + 1)
    xs = x_min + h * np.arange(1, n + 1)
    diag = 2.0 / h**2 + potential(xs)
    off = -np.ones(n - 1) / h**2
    mat = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    vals = np.linalg.eigvalsh(mat)
    return [float(v) for v in vals[:count]]


def count_negative(vals: list[float], tol: float = SPECTRUM_TOL) -> int:
    return sum(1 for val in vals if val < -tol)


def scan_imposed_terms() -> dict:
    sources = [PDE_TEX, CONCEPTUAL, PATHA23_STAGE0]
    terms = {
        "V_conf": ["V_conf", "V_{\\mathrm{conf}}"],
        "Z(w)": ["Z(w)"],
        "W(w)": ["W(w)"],
        "B_ell": ["B_\\ell", "B_ℓ"],
        "k_w": ["k_w"],
    }
    found: dict[str, list[str]] = {}
    for term, needles in terms.items():
        hits: list[str] = []
        for source in sources:
            if not source.exists():
                continue
            lines = source.read_text(encoding="utf-8").splitlines()
            for lineno, line in enumerate(lines, start=1):
                if any(needle in line for needle in needles):
                    hits.append(f"{source.relative_to(PROJECT_ROOT)}:{lineno}")
                    break
        found[term] = hits
    return found


def symbolic_derivations() -> dict:
    w, rho = sp.symbols("w rho", nonnegative=True, real=True)
    rho0, K, m, a, hbar, Lhalf = sp.symbols(
        "rho0 K m a hbar Lhalf", positive=True
    )
    pnorm2, pprime2, rhop = sp.symbols("P2 Pprime2 rhop", nonnegative=True, real=True)
    cs2 = sp.Rational(5) * K * rho**4 / m
    cs20 = sp.Rational(5) * K * rho0**4 / m
    U = K * rho**5 / 4
    static_density = (
        hbar**2 * rhop**2 / (8 * m * rho)
        + U
        + sp.Rational(1, 2) * m * rho * cs2 * a**2 * pprime2
        + sp.Rational(1, 4) * m * rho * cs2 * (pnorm2 - 1) ** 2
    )
    static_density_substituted = sp.factor(static_density)

    dU = sp.diff(U, rho)
    d2U = sp.diff(U, rho, 2)
    stationary_points = sp.solve(sp.Eq(dU, 0), rho)
    if stationary_points != [0]:
        raise AssertionError(f"unexpected U stationary points: {stationary_points}")
    if sp.simplify(d2U - 5 * K * rho**3) != 0:
        raise AssertionError(f"unexpected U'' = {d2U}")

    Kp = sp.factor(m * rho0 * cs20 * a**2)
    B = sp.factor(m * rho0 * cs20)
    sigma_spread = sp.factor(Kp * sp.pi**2 / (4 * Lhalf))
    if sp.limit(sigma_spread, Lhalf, sp.oo) != 0:
        raise AssertionError("spread minimizing sequence does not have zero tension limit")
    dsigma_dL = sp.diff(sigma_spread, Lhalf)
    if sp.simplify(dsigma_dL + Kp * sp.pi**2 / (4 * Lhalf**2)) != 0:
        raise AssertionError("unexpected d sigma_spread / dL")

    sigma_radial_saddle = sp.factor(sp.Rational(2, 3) * sp.sqrt(2) * B * a)
    radial_negative_omega2 = sp.factor(-cs20 / (2 * a**2))
    confinement_gap_spread = sp.limit(cs20 * sp.pi**2 / (4 * Lhalf**2), Lhalf, sp.oo)
    if confinement_gap_spread != 0:
        raise AssertionError("delocalized confinement gap should vanish")

    return {
        "static_energy_density": sp.sstr(static_density_substituted),
        "U": {
            "expr": sp.sstr(U),
            "dU": sp.sstr(dU),
            "d2U": sp.sstr(d2U),
            "stationary_points_nonnegative": [sp.sstr(x) for x in stationary_points],
            "single_nonnegative_minimum": "rho=0; no second scalar vacuum",
        },
        "coefficients": {
            "c_s0_squared": sp.sstr(cs20),
            "K_P0": sp.sstr(Kp),
            "B0": sp.sstr(B),
        },
        "fixed_background_minimizing_sequence": {
            "domain": "w in [-Lhalf, Lhalf]",
            "profile": "P_L(w)=(sin(pi*(w+Lhalf)/(2*Lhalf)),0,0,-cos(pi*(w+Lhalf)/(2*Lhalf)))",
            "rho": "rho0",
            "core_at_w0": "P_L(0)=(1,0,0,0), flat/in-plane for this finite-box representative",
            "sigma_L": sp.sstr(sigma_spread),
            "d_sigma_d_Lhalf": sp.sstr(sp.factor(dsigma_dL)),
            "limit_Lhalf_to_infinity": "0",
            "width_status": "arbitrary_box_size_dependent_no_natural_localized_width",
        },
        "radial_amplitude_kink_saddle": {
            "profile": "P=(0,0,0,tanh(w/(sqrt(2)*a)))",
            "sigma_saddle": sp.sstr(sigma_radial_saddle),
            "core_texture": "amplitude_collapse_P=0_not_flat_orientation",
            "transverse_negative_omega_squared": sp.sstr(radial_negative_omega2),
            "not_minimizer_reason": "finite fixed sigma but beaten by sigma_L -> 0 spread sequence",
        },
        "topology": {
            "vacuum_manifold": "{P in R^4 | P.P = 1} = S^3",
            "sphere_dimension": 3,
            "pi0": 0,
            "connected": True,
            "topologically_protected": False,
        },
        "unwinding": {
            "deltaE_unwind": "0",
            "tau_unwind": "not_computed_no_metastable_minimum",
            "tau_threshold_preregistered_before_scan": TAU_THRESHOLD,
            "local_minimum_against_delocalization": False,
        },
        "confinement": {
            "gap_spread_limit_omega_squared": sp.sstr(confinement_gap_spread),
            "bound_zero_modes_from_wall": False,
        },
    }


def spectrum_checks() -> dict:
    n = 401
    x_max = 14.0
    phi4_vals = finite_difference_spectrum(
        lambda x: 2.0 - 3.0 / np.cosh(x / math.sqrt(2.0)) ** 2,
        -x_max,
        x_max,
        n,
        4,
    )
    unstable_vals = finite_difference_spectrum(
        lambda x: -np.ones_like(x),
        -x_max,
        x_max,
        n,
        4,
    )
    Lgeo = 10.0
    q = math.pi / (2.0 * Lgeo)
    connected_vals = finite_difference_spectrum(
        lambda x: -q**2 * np.ones_like(x),
        -Lgeo,
        Lgeo,
        n,
        4,
    )
    radial_vals = finite_difference_spectrum(
        lambda y: -2.0 / np.cosh(y) ** 2,
        -x_max,
        x_max,
        n,
        4,
    )

    controls = {
        "phi4_kink_topologically_stable": {
            "operator": "-d^2/dx^2 + 2 - 3 sech^2(x/sqrt(2))",
            "pi0": "Z2",
            "lowest_eigenvalues": phi4_vals,
            "negative_count": count_negative(phi4_vals),
            "expected": "no negative modes; translation near-zero; protected",
            "pass": count_negative(phi4_vals) == 0 and phi4_vals[0] > -SPECTRUM_TOL,
        },
        "known_unstable_saddle_phi_at_zero": {
            "operator": "-d^2/dx^2 - 1",
            "lowest_eigenvalues": unstable_vals,
            "negative_count": count_negative(unstable_vals),
            "expected": "negative mode found",
            "pass": count_negative(unstable_vals) >= 1,
        },
        "connected_vacuum_clamped_antipodal_geodesic": {
            "operator": "-d^2/dw^2 - (pi/(2L))^2, Dirichlet endpoints",
            "L": Lgeo,
            "pi0": 0,
            "lowest_eigenvalues": connected_vals,
            "negative_count_below_tol": count_negative(connected_vals),
            "expected": "clean clamped spectrum but topological diagnosis flags unwinding",
            "pass": count_negative(connected_vals) == 0 and abs(connected_vals[0]) < 5.0e-4,
        },
    }
    baseline = {
        "finite_box_geodesic_sensitivity": controls["connected_vacuum_clamped_antipodal_geodesic"],
        "radial_soft_spin_saddle_transverse": {
            "operator_dimensionless": "-d^2/dy^2 - 2 sech^2(y)",
            "analytic_lowest_eigenvalue_dimensionless": "-1",
            "lowest_eigenvalues_dimensionless": radial_vals,
            "negative_count": count_negative(radial_vals),
            "expected": "negative transverse mode for the amplitude-collapse saddle",
            "pass": count_negative(radial_vals) >= 1 and abs(radial_vals[0] + 1.0) < 5.0e-3,
        },
    }
    all_pass = all(item["pass"] for item in controls.values()) and baseline[
        "radial_soft_spin_saddle_transverse"
    ]["pass"]
    if not all_pass:
        raise AssertionError(f"spectrum fixture failure: {json.dumps(controls, indent=2)}")
    return {
        "tolerance_negative_threshold": SPECTRUM_TOL,
        "controls": controls,
        "baseline": baseline,
    }


def main() -> int:
    freeze_guard = assert_freeze_fidelity()
    dimensions = check_dimensions()
    derivations = symbolic_derivations()
    spectra = spectrum_checks()
    imposed_terms = scan_imposed_terms()

    labels = {
        "T1a": "BRANE_IMPOSED_NOT_DERIVED",
        "T1b": "FAIL_NO_WALL_PROFILE",
        "T1c": "FAIL_WALL_UNWINDS_SPHERE_VACUA",
        "T1d": "W_EMERGENT_BUT_UNSTABLE",
        "T1e": "FAIL_NO_BOUND_ZERO_MODES",
        "rollup": "T1_FAIL_NO_STABLE_WALL",
    }
    report = {
        "schema": "stage1_pathA_24_T1_wall_sympy/v1",
        "engine": "sympy+numpy",
        "pass": True,
        "freeze_fidelity_guard": freeze_guard,
        "labels": labels,
        "dimensions": dimensions,
        "derivations": derivations,
        "spectra": spectra,
        "imposed_existing_brane_terms_source_scan": imposed_terms,
    }
    SCRATCH.mkdir(parents=True, exist_ok=True)
    out_path = SCRATCH / "pathA_24_T1_wall_sympy.json"
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("PATHA_24_T1_WALL_SYMPY: PASS")
    print(f"freeze_sha256: {freeze_guard['sha256']}")
    print(f"T1_rollup: {labels['rollup']}")
    print(f"T1a: {labels['T1a']}")
    print(f"T1b: {labels['T1b']}")
    print(f"T1c: {labels['T1c']}")
    print(f"T1d: {labels['T1d']}")
    print(f"T1e: {labels['T1e']}")
    print(f"sigma_L: {derivations['fixed_background_minimizing_sequence']['sigma_L']}")
    print(f"sigma_L_limit: {derivations['fixed_background_minimizing_sequence']['limit_Lhalf_to_infinity']}")
    print(f"radial_saddle_sigma: {derivations['radial_amplitude_kink_saddle']['sigma_saddle']}")
    print(f"pi0_vacuum: {derivations['topology']['pi0']}")
    print(f"deltaE_unwind: {derivations['unwinding']['deltaE_unwind']}")
    print(f"tau_threshold: {TAU_THRESHOLD}")
    print(f"confinement_gap_limit_omega2: {derivations['confinement']['gap_spread_limit_omega_squared']}")
    print(f"controls_pass: {all(item['pass'] for item in spectra['controls'].values())}")
    print(f"wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
