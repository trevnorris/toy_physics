from __future__ import annotations

from pathlib import Path

import pytest

from stage1_solver.patha22b_gate1 import (
    exported_z_quadrature,
    gate1_report_payload,
    measure_determination,
    resolution_diagnostics,
    target_blindness_guard,
)


def test_gate1_uses_pde_flat_dw_measure_and_flags_sqrtg_discrepancy() -> None:
    quad = exported_z_quadrature()
    measure = measure_determination(quad)

    assert measure["pde_z_int_measure"] == "flat_dw_no_sqrt_g_w"
    assert measure["sqrt_g_w_exported"] is False
    assert measure["sqrt_g_w_variant_status"] == "NOT_INDEPENDENTLY_COMPUTABLE_FROM_EXPORTED_Z_W"
    assert "pde.tex defines Z_int with flat dw" in measure["discrepancy_flag"]
    assert measure["sqrt_g_w_variant_if_flat_or_densitized"] == pytest.approx(
        quad["z_int_finite_domain"]
    )


def test_gate1_quadrature_reports_finite_exported_grid_value_and_error_bar() -> None:
    payload = gate1_report_payload()

    assert payload["z_int_value"] == pytest.approx(2.031114372358842)
    assert payload["z_int_error_bar_grid"] == pytest.approx(0.001118740362587922)
    assert payload["z_int_units"] == "L"
    assert payload["gate1_outcome"] == "BLOCKED_NEEDS_DECAYING_Z_PROFILE_OR_FLOOR_PROVENANCE"


def test_gate1_tail_diagnostic_blocks_small_ideal_tail_claim() -> None:
    quad = exported_z_quadrature()

    assert quad["ideal_infinite_integral_status"] == "BLOCKED_NEEDS_DECAYING_OR_ZERO_FLOOR_EXTENSION"
    assert quad["z_edge_max_abs"] == pytest.approx(0.940421257709774)
    assert quad["edge_to_peak_abs_ratio"] > 0.7
    assert quad["edge_cell_integral_fraction_abs"] > 0.2


def test_gate1_reports_floor_dominated_finite_box_decomposition() -> None:
    quad = exported_z_quadrature()

    assert quad["localization_floor"] == pytest.approx(0.8)
    assert quad["grid"]["domain_width"] == pytest.approx(1.85)
    assert quad["floor_integral_finite_domain"] == pytest.approx(1.48)
    assert quad["gaussian_integral_finite_domain"] == pytest.approx(
        quad["z_int_finite_domain"] - 1.48
    )
    assert quad["floor_fraction_finite_domain"] == pytest.approx(1.48 / 2.031114372358842)
    assert quad["floor_fraction_finite_domain"] > 0.7
    assert quad["gaussian_fraction_finite_domain"] < 0.3


def test_gate1_resolution_diagnostics_are_from_existing_background_exports() -> None:
    diagnostics = resolution_diagnostics()
    rows = diagnostics["rows"]

    assert [row["nw"] for row in rows] == [4, 6, 8]
    assert rows[-1]["z_int_finite_domain"] == pytest.approx(2.031114372358842)
    assert diagnostics["nearest_grid_delta"] == pytest.approx(0.001118740362587922)
    assert diagnostics["ladder_spread"] == pytest.approx(0.004299348740518966)


def test_gate1_dimensional_checks_confirm_z_dimensionless_and_zint_length() -> None:
    payload = gate1_report_payload()
    checks = {check["name"]: check for check in payload["dimensional_checks"]["checks"]}

    assert checks["exported Z_w localization weight"]["status"] == "CONSISTENT"
    assert checks["exported Z_w localization weight"]["actual"] == "1"
    assert checks["Z_int=int Z(w) dw"]["status"] == "CONSISTENT"
    assert checks["Z_int=int Z(w) dw"]["actual"] == "L"


def test_gate1_scope_is_coupling_artifact_not_lambda_gamma() -> None:
    payload = gate1_report_payload()

    assert "mu0_eff=mu0/Z_int" in payload["scope_statement"]
    assert "q_eff=q_star/sqrt(Z_int)" in payload["scope_statement"]
    assert "P0*chi_Q*g_mhat^2*lambda_gamma^5/g_G" in payload["scope_statement"]
    assert "does not gate the xi verdict" in payload["scope_statement"]
    assert "never as the numeric finite-box value" in payload["scope_statement"]
    assert "does not promote Z_int to lambda_gamma" in payload["scope_statement"]
    assert any("lambda_gamma" in residual for residual in payload["residuals"])
    assert any("localization_floor=0.8" in residual for residual in payload["residuals"])


def test_gate1_target_blindness_guard_scans_new_sources() -> None:
    root = Path(__file__).resolve().parents[1]
    paths = [
        root / "src/stage1_solver/patha22b_gate1.py",
        root / "tools/pathA_22b_gate1_crosscheck.wl",
    ]

    result = target_blindness_guard(paths)

    assert result["status"] == "TARGET_BLIND_PASS"
    assert result["hits"] == []
