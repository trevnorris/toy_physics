from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

from stage1_solver.patha22b_gate3 import (
    DEFAULT_DOMAIN_RADIAL_SCALES,
    DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE,
    DEFAULT_GRID_SWEEP,
    DEFAULT_NW_TAIL_BUDGET_MIN_NW,
    DEFAULT_NW_TAIL_CENTRAL_MIN_NW,
    OUTCOME_NOT_EXTRACTABLE,
    OUTCOME_NW_CONVERGES_PREFIX,
    OUTCOME_PREFIX_DELTA,
    assess_domain_plateau,
    assess_grid_convergence,
    assemble_gate3_sparse_vsh_maxwell_system,
    canonical_fingerprint_coefficient,
    canonical_sigma_q,
    characterize_even_nw_axis,
    chi_q_from_omega5_coefficient,
    compare_nw_trends_across_nr,
    apply_gate3_w_lane_cleaner,
    centered_w_derivative_null_diagnostic,
    dimensional_checks,
    domain_truncation_sweep,
    exact_outgoing_yhat_l2,
    fit_nw_free_power_law,
    fit_nw_power_law,
    extract_iomega5_coefficient,
    frozen_branch_inventory,
    gate3_report_payload,
    default_nw_characterization_path,
    jackknife_nw_limit,
    load_frozen_solve_inputs,
    make_branch_geometry_reference_background,
    make_uniform_vacuum_background,
    outgoing_hankel_lambda_series,
    radius_invariance_sweep,
    richardson_geometric_extrapolation,
    solve_outgoing_dtn_sweep,
    summarize_nw_characterization,
    target_blindness_guard,
    tiny_defect_linearity_scan,
    vacuum_sanity_check,
)
from stage1_solver import patha_b2b_maxwell as b2b


def test_gate3_hankel_fingerprint_and_sigma_match_source_formula() -> None:
    a, c_s = sp.symbols("a c_s", positive=True)
    hankel = outgoing_hankel_lambda_series()

    assert "I*z**5/9" in hankel["lambda_series"]
    assert "I*z**5/27" in hankel["yhat_series"]
    assert canonical_fingerprint_coefficient(a, c_s) == a**5 / (27 * c_s**5)
    assert canonical_sigma_q(a, c_s) == 4 * a**5 / (27 * c_s**5)


def test_gate3_omega5_extraction_is_discriminating() -> None:
    omega, a, c_s, q = sp.symbols("omega a c_s q", positive=True)
    fingerprint = canonical_fingerprint_coefficient(a, c_s)
    response = 1 + a**2 * omega**2 / (9 * c_s**2) + sp.I * q * fingerprint * omega**5

    coeff = extract_iomega5_coefficient(response, omega)
    chi_q = chi_q_from_omega5_coefficient(coeff, a=a, c_s=c_s)

    assert sp.simplify(coeff - q * fingerprint) == 0
    assert sp.simplify(chi_q - q) == 0


def test_gate3_exact_outgoing_yhat_numeric_moves_with_radius_argument() -> None:
    small = exact_outgoing_yhat_l2(0.02)
    larger = exact_outgoing_yhat_l2(0.04)

    assert small.imag > 0.0
    assert larger.imag > small.imag


def test_gate3_centered_w_derivative_has_odd_checkerboard_null() -> None:
    odd = centered_w_derivative_null_diagnostic(13, 1.0)
    even = centered_w_derivative_null_diagnostic(14, 1.0)

    assert odd["is_exact_null"] is True
    assert odd["residual_norm"] < 1.0e-12
    assert odd["even_index_norm"] > 0.99
    assert odd["odd_index_norm"] < 1.0e-12
    assert even["is_exact_null"] is False
    assert even["min_singular"] > 0.0


def test_gate3_w_null_lift_is_gate3_local_and_exact_null_scoped() -> None:
    inputs = load_frozen_solve_inputs()
    odd_sys = b2b.assemble_vsh_maxwell_system(
        inputs["maxwell_background"],
        nr=5,
        nw=5,
        omega=0.034,
        radial_scale=1.0,
    )
    even_sys = b2b.assemble_vsh_maxwell_system(
        inputs["maxwell_background"],
        nr=5,
        nw=6,
        omega=0.034,
        radial_scale=1.0,
    )

    odd_clean = apply_gate3_w_lane_cleaner(odd_sys, cleaner="exact_null_phi_lift", strength=1.0e-2)
    even_clean = apply_gate3_w_lane_cleaner(even_sys, cleaner="exact_null_phi_lift", strength=1.0e-2)

    assert odd_clean["Gate3WNullLift"]["active"] is True
    assert even_clean["Gate3WNullLift"]["active"] is False
    assert odd_clean["KCons"].shape == odd_sys["KCons"].shape
    assert even_clean["KCons"].shape == even_sys["KCons"].shape
    assert np.linalg.norm(odd_clean["KCons"] - odd_sys["KCons"]) > 0.0
    assert np.linalg.norm(even_clean["KCons"] - even_sys["KCons"]) == 0.0


def test_gate3_inventory_records_legacy_seed_but_not_as_solve_input() -> None:
    inventory = frozen_branch_inventory()

    assert inventory["background_schema"] == "stage1_m1c_physical_wp1_background/v1"
    assert inventory["background_residual_linf"] < 5.0e-8
    assert inventory["closure_radius"] == 1.5
    assert inventory["mouth_radius"] == 1.0
    assert inventory["legacy_flux_normalization_present_but_not_used"]["Gamma_port"] is not None


def test_gate3_outgoing_dtn_solve_reports_raw_coefficient_without_fixed_a_chi_q() -> None:
    inputs = load_frozen_solve_inputs()
    solve = solve_outgoing_dtn_sweep(
        maxwell_background=inputs["maxwell_background"],
        source_packet=inputs["source_packet"],
        closure_radius=1.5,
        nr=6,
        nw=6,
    )

    assert solve["closure"]["type"] == "spherical_hankel_l2_normalized_dtn"
    assert "chi_Q" not in solve["fit"]
    assert solve["fit"]["effective_omega5_coefficient"] > 0.0
    assert solve["fit"]["uncalibrated_fixed_a_ratio"] > 0.0
    assert solve["diagnostics"]["max_linear_residual"] < 1.0e-12
    assert solve["fit"]["fit_rms"] < 1.0e-8


def test_gate3_sparse_assembly_matches_shared_dense_operator_on_small_grid() -> None:
    inputs = load_frozen_solve_inputs()
    dense = b2b.assemble_vsh_maxwell_system(
        inputs["maxwell_background"],
        nr=4,
        nw=4,
        omega=0.034,
        radial_scale=1.0,
        discretization="primary_second_order",
    )
    sparse = assemble_gate3_sparse_vsh_maxwell_system(
        inputs["maxwell_background"],
        nr=4,
        nw=4,
        omega=0.034,
        radial_scale=1.0,
    )

    assert sparse["matrix_assembly"] == "gate3_sparse"
    for key in ("KCons", "ACons", "BoundaryB", "WLaneInv"):
        assert np.max(np.abs(dense[key] - sparse[key].toarray())) < 1.0e-12


def test_gate3_sparse_direct_solve_uses_sparse_assembly_and_keeps_residual_small() -> None:
    inputs = load_frozen_solve_inputs()
    solve = solve_outgoing_dtn_sweep(
        maxwell_background=inputs["maxwell_background"],
        source_packet=inputs["source_packet"],
        closure_radius=1.5,
        nr=4,
        nw=4,
        linear_solver="sparse_direct",
        matrix_assembly="gate3_sparse",
    )

    assert solve["w_lane_control"]["linear_solver"] == "sparse_direct"
    assert solve["w_lane_control"]["matrix_assembly"] == "gate3_sparse"
    assert solve["diagnostics"]["max_linear_residual"] < 1.0e-12
    assert solve["fit"]["effective_omega5_coefficient"] > 0.0


def test_gate3_radius_consistent_chi_q_is_invariant_to_closure_radius() -> None:
    inputs = load_frozen_solve_inputs()
    vacuum = make_uniform_vacuum_background(inputs["maxwell_background"])
    rows = radius_invariance_sweep(
        defect_background=inputs["maxwell_background"],
        vacuum_background=vacuum,
        source_packet=inputs["source_packet"],
        radii=(1.0, 1.5, 2.0),
        grid=(6, 6),
    )
    chi_values = [row["chi_Q"] for row in rows]
    legacy_values = [row["defect_uncalibrated_fixed_a_ratio"] for row in rows]

    assert max(chi_values) - min(chi_values) < 1.0e-8
    assert legacy_values[2] > legacy_values[1] > legacy_values[0]


def test_gate3_branch_geometry_reference_keeps_zw_and_r0() -> None:
    inputs = load_frozen_solve_inputs()
    reference = make_branch_geometry_reference_background(inputs["maxwell_background"])

    assert reference["derived"]["Z_w"] == inputs["maxwell_background"]["derived"]["Z_w"]
    assert reference["fields"]["R0_w"] == inputs["maxwell_background"]["fields"]["R0_w"]
    assert reference["fields"]["A_00"][0][0] == 0.0
    assert reference["fields"]["psi_R0"][0][0] == 1.0


def test_gate3_domain_truncation_sweep_is_not_closure_invariance() -> None:
    inputs = load_frozen_solve_inputs()
    reference = make_branch_geometry_reference_background(inputs["maxwell_background"])
    rows = domain_truncation_sweep(
        defect_background=inputs["maxwell_background"],
        vacuum_background=reference,
        source_packet=inputs["source_packet"],
        radial_scales=(1.0, 1.25, 1.5),
        base_grid=(6, 6),
    )
    assessment = assess_domain_plateau(rows)

    assert [row["grid"]["r_max"] for row in rows] == [2.0, 2.5, 3.0]
    assert rows[-1]["grid"]["nr"] > rows[0]["grid"]["nr"]
    assert rows[-1]["chi_Q"] > rows[0]["chi_Q"]
    assert assessment["plateau"] is False


def test_gate3_extended_rmax_sweep_reaches_plateau() -> None:
    inputs = load_frozen_solve_inputs()
    reference = make_branch_geometry_reference_background(inputs["maxwell_background"])
    rows = domain_truncation_sweep(
        defect_background=inputs["maxwell_background"],
        vacuum_background=reference,
        source_packet=inputs["source_packet"],
        radial_scales=DEFAULT_DOMAIN_RADIAL_SCALES,
        base_grid=(6, 6),
    )
    assessment = assess_domain_plateau(rows)

    assert [row["grid"]["r_max"] for row in rows] == [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    assert rows[-1]["grid"]["nr"] > rows[0]["grid"]["nr"]
    assert assessment["plateau"] is True
    assert assessment["plateau_onset_r_max"] == 6.0
    assert assessment["strict_1e_minus_3_onset_r_max"] == 7.0
    assert assessment["abs_deltas"][-1] < 1.0e-3


def test_gate3_richardson_extrapolates_below_descending_finest_value() -> None:
    rows = [
        {"chi_Q": 0.1807112082040},
        {"chi_Q": 0.1692434814186},
        {"chi_Q": 0.1670570231357},
        {"chi_Q": 0.1666270056256},
    ]

    richardson = richardson_geometric_extrapolation(rows)

    assert richardson["limit"] < richardson["finest_value"]
    assert richardson["limit"] < richardson["two_grid_mean"]
    assert 0.0 < richardson["observed_delta_ratio"] < 1.0


def test_gate3_post_plateau_grid_sequence_extrapolates_large_rmax_limit() -> None:
    rows = [
        {"chi_Q": 0.7411559775808990},
        {"chi_Q": 0.7220740387947787},
        {"chi_Q": 0.7159045107266220},
    ]

    richardson = richardson_geometric_extrapolation(rows)
    assessment = assess_grid_convergence(rows, richardson)

    assert DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE == 3.0
    assert DEFAULT_GRID_SWEEP[-1][1] > 8
    assert assessment["converging"] is True
    assert richardson["limit"] < richardson["finest_value"]
    assert abs(richardson["limit"] - 0.7129567206187051) < 1.0e-12
    assert assessment["rel_tail_correction"] < 5.0e-3


def _synthetic_nw_rows(
    *,
    nr: int,
    chi_inf: float = 0.72,
    c: float = 1.5,
    p: float = 1.25,
    offset: float = 0.0,
) -> list[dict[str, object]]:
    return [
        {
            "grid": {"nr": nr, "nw": nw},
            "chi_Q": chi_inf + offset + c / float(nw) ** p,
            "max_linear_residual": 1.0e-14,
            "max_fit_rms": 1.0e-10,
        }
        for nw in (8, 10, 12, 14, 16, 20, 24, 28, 32, 40)
    ]


def test_gate3_nw_power_fit_jackknife_and_nr_confirmation_on_synthetic_limit() -> None:
    primary = _synthetic_nw_rows(nr=48)
    secondary = _synthetic_nw_rows(nr=61, offset=1.0e-4)

    fixed_1 = fit_nw_power_law(primary, p=1.0)
    free = fit_nw_free_power_law(primary)
    jackknife = jackknife_nw_limit(primary, free_p=True)
    comparison = compare_nw_trends_across_nr(primary, secondary, free, fit_nw_free_power_law(secondary))
    characterization = characterize_even_nw_axis(primary_rows=primary, secondary_rows=secondary)

    assert fixed_1["rms"] > free["rms"]
    assert abs(free["chi_inf"] - 0.72) < 1.0e-8
    assert abs(free["p"] - 1.25) < 1.0e-6
    assert jackknife["stable"] is True
    assert comparison["nr_independent"] is True
    assert characterization["outcome_class"].startswith(OUTCOME_NW_CONVERGES_PREFIX)


def test_gate3_stored_even_nw_characterization_sets_production_central_and_budget() -> None:
    path = default_nw_characterization_path(Path(__file__).resolve().parents[3])
    summary = summarize_nw_characterization(json.loads(path.read_text(encoding="utf-8")), source_path=path)

    assert summary["converged"] is True
    assert summary["outcome_class"].startswith(OUTCOME_NW_CONVERGES_PREFIX)
    assert summary["central_tail"]["min_nw"] == DEFAULT_NW_TAIL_CENTRAL_MIN_NW
    assert summary["budget_tail"]["min_nw"] == DEFAULT_NW_TAIL_BUDGET_MIN_NW
    assert summary["chi_Q_reported"] == 0.712
    assert summary["numerical_tail_supported"] == 0.0008
    assert abs(summary["full_sweep_fit_uncertainty"] - 0.0042) < 5.0e-5


def test_gate3_vacuum_sanity_returns_one_and_tiny_defect_tends_to_one() -> None:
    inputs = load_frozen_solve_inputs()
    sanity = vacuum_sanity_check(
        physical_background=inputs["maxwell_background"],
        source_packet=inputs["source_packet"],
        radii=(1.0, 1.5, 2.0),
        grid=(6, 6),
    )

    assert sanity["passed"] is True
    assert sanity["max_uniform_abs_error"] < 1.0e-10
    assert sanity["tiny_defect_abs_error"] < 1.0e-3


def test_gate3_tiny_defect_linearity_scan_has_stable_small_fraction_slope() -> None:
    inputs = load_frozen_solve_inputs()
    reference = make_branch_geometry_reference_background(inputs["maxwell_background"])
    scan = tiny_defect_linearity_scan(
        physical_background=inputs["maxwell_background"],
        reference_background=reference,
        source_packet=inputs["source_packet"],
        fractions=(1.0e-3, 3.0e-3, 1.0e-2),
        grid=(6, 6),
    )

    assert scan["linear_toward_zero"] is True
    assert len(scan["rows"]) == 3
    assert max(abs(row["chi_Q_minus_1"]) for row in scan["rows"]) < 2.0e-2


def test_gate3_payload_reports_domain_nonplateau_as_not_extractable() -> None:
    payload = gate3_report_payload(
        radius_sweep=(1.0, 1.5),
        grid_sweep=((4, 4), (5, 5), (6, 6)),
        vacuum_sanity_radii=(1.5,),
        domain_base_grid=(4, 4),
        largest_domain_grid_bases=((4, 4), (5, 5), (6, 6)),
        tiny_defect_grid=(4, 4),
    )

    assert not payload["gate3_outcome"].startswith(OUTCOME_PREFIX_DELTA)
    assert payload["gate3_outcome"].startswith(OUTCOME_NOT_EXTRACTABLE)
    assert payload["blocker"] != ""
    assert payload["chi_Q_value"] is None
    assert payload["result"]["can_report_chi_Q_number"] is False
    assert payload["result"]["convergence"]["radius_pass"] is False
    assert payload["result"]["tiny_defect_linearity"]["linear_toward_zero"] is True
    assert "calibration knob" in payload["result"]["reason"]


def test_gate3_dimensional_checks_confirm_chi_q_dimensionless() -> None:
    checks = {check["name"]: check for check in dimensional_checks()["checks"]}

    assert checks["defect omega^5 coefficient"]["status"] == "CONSISTENT"
    assert checks["defect omega^5 coefficient"]["actual"] == "T^5"
    assert checks["vacuum omega^5 coefficient"]["status"] == "CONSISTENT"
    assert checks["vacuum omega^5 coefficient"]["actual"] == "T^5"
    assert checks["closure placement factor (R_exit/c_s)^5"]["status"] == "CONSISTENT"
    assert checks["closure placement factor (R_exit/c_s)^5"]["actual"] == "T^5"
    assert checks["radius-consistent chi_Q ratio"]["status"] == "CONSISTENT"
    assert checks["radius-consistent chi_Q ratio"]["actual"] == "1"


def test_gate3_target_blindness_guard_scans_new_sources() -> None:
    root = Path(__file__).resolve().parents[1]
    paths = [
        root / "src/stage1_solver/patha22b_gate3.py",
        root / "tools/pathA_22b_gate3_crosscheck.wl",
        root / "tests/test_patha22b_gate3.py",
        root.parents[1] / "_scratch/pathA_22b_gate3_nw_characterization.py",
    ]

    result = target_blindness_guard(paths)

    assert result["status"] == "TARGET_BLIND_PASS"
    assert result["hits"] == []
