from __future__ import annotations

from pathlib import Path

from stage1_solver.patha22b_gate2 import (
    OUTCOME_CLASS,
    gate2_report_payload,
    localized_maxwell_reduction,
    negative_control_cgamma_equals_cs,
    speed_normalization_result,
    target_blindness_guard,
)


def test_gate2_reduces_localized_maxwell_coefficients_and_cancels_zint() -> None:
    reduction = localized_maxwell_reduction()

    assert reduction["measure"] == "flat_dw_no_sqrt_g_w"
    assert reduction["metric_finding"] == "FLAT_UNWARPED_BULK_METRIC"
    assert reduction["C_E"] == "Z_int/mu0"
    assert reduction["C_B"] == "Z_int/mu0"
    assert reduction["C_B_over_C_E"] == "1"
    assert reduction["F_squared_zero_mode"] == "-2*F_tx**2 - 2*F_ty**2 - 2*F_tz**2 + 2*F_xy**2 + 2*F_xz**2 + 2*F_yz**2"
    assert reduction["F_squared_zero_mode_check"] is True
    assert reduction["z_cancel_control"]["pass"] is True
    assert reduction["z_cancel_control"]["actual"] == "1"


def test_gate2_carries_speed_map_and_conditional_lambda_gamma() -> None:
    speed = speed_normalization_result()

    assert speed["speed_map_status"] == "UNPINNED_BY_SOURCES"
    assert speed["open_knob_count"] == 1
    assert speed["closed_result_available"] is False
    assert speed["outcome"] == OUTCOME_CLASS
    assert speed["c_bulk_squared_reduced"] == "1"
    assert speed["c_gamma_squared"] == "beta_bulk_to_brane**2"
    assert (
        speed["lambda_gamma_carried"]
        == "beta_bulk_to_brane*sqrt(m_GNLS/(5*K*rho0**4))"
    )


def test_gate2_negative_control_rejects_forced_cgamma_equals_cs() -> None:
    control = negative_control_cgamma_equals_cs()

    assert control["verdict"] == "GUARD_AGAINST_FORCING_LAMBDAGAMMA_TO_ONE"
    assert control["guard_not_discriminating_test"] is True
    assert control["residual_is_identity_zero"] is False
    assert control["forced_valid"] is False
    assert control["respects_pathA_20b_negative_control"] is True
    assert "5*K*rho0**4" in control["symbolic_residual"]


def test_gate2_dimensional_checks_confirm_lambda_gamma_dimensionless() -> None:
    payload = gate2_report_payload()
    checks = {check["name"]: check for check in payload["dimensional_checks"]["checks"]}

    assert checks["C_B/C_E flat bulk cone ratio"]["status"] == "CONSISTENT"
    assert checks["C_B/C_E flat bulk cone ratio"]["actual"] == "1"
    assert checks["speed-normalization beta_bulk_to_brane"]["status"] == "CONSISTENT"
    assert checks["speed-normalization beta_bulk_to_brane"]["actual"] == "L T^-1"
    assert checks["c_gamma=beta*sqrt(C_B/C_E)"]["status"] == "CONSISTENT"
    assert checks["c_gamma=beta*sqrt(C_B/C_E)"]["actual"] == "L T^-1"
    assert checks["c_s=sqrt(5*K*rho0^4/m_GNLS)"]["status"] == "CONSISTENT"
    assert checks["c_s=sqrt(5*K*rho0^4/m_GNLS)"]["actual"] == "L T^-1"
    assert checks["lambda_gamma=c_gamma/c_s"]["status"] == "CONSISTENT"
    assert checks["lambda_gamma=c_gamma/c_s"]["actual"] == "1"


def test_gate2_payload_reports_conditional_outcome_and_patha20b_group() -> None:
    payload = gate2_report_payload()

    assert payload["gate2_outcome"] == OUTCOME_CLASS
    assert "pathA_20b_cgamma_cs_linearization.md:24-50" in payload["existing_pathA_20b_group"]
    assert payload["negative_control"]["verdict"] == "GUARD_AGAINST_FORCING_LAMBDAGAMMA_TO_ONE"
    assert payload["speed_normalization"]["open_knob_count"] == 1
    assert any("Z_int cancels" in residual for residual in payload["residuals"])
    assert any("beta_bulk_to_brane" in residual for residual in payload["residuals"])
    assert any("FAIL_ABLE_PENDING_LAMBDAGAMMA_SPEED_MAP" in residual for residual in payload["residuals"])


def test_gate2_target_blindness_guard_scans_new_sources() -> None:
    root = Path(__file__).resolve().parents[1]
    paths = [
        root / "src/stage1_solver/patha22b_gate2.py",
        root / "tools/pathA_22b_gate2_crosscheck.wl",
        root / "tests/test_patha22b_gate2.py",
    ]

    result = target_blindness_guard(paths)

    assert result["status"] == "TARGET_BLIND_PASS"
    assert result["hits"] == []
