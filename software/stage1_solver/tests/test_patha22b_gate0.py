from __future__ import annotations

from pathlib import Path

import sympy as sp

from stage1_solver.patha22b_gate0 import (
    DOES_NOT_CANCEL_SOURCE_LABEL,
    classify_kernel_pair,
    gate0_report_payload,
    mhat_reconciliation,
    target_blindness_guard,
    z_kernel_cancellation_source_assessment,
)


def test_gate0a_confirms_dimensionful_mhat_and_rejects_direct_dimensionless_reading() -> None:
    report = mhat_reconciliation()
    checks = {check["name"]: check for check in report["checks"]}

    assert report["outcome"] == "MHAT_DIMENSIONFUL_CONFIRMED"
    assert checks["required mhat from mhat^2*Gamma5=G/c^5"]["status"] == "CONSISTENT"
    assert (
        checks["direct dimensionless mhat reading in odd-coefficient law"]["status"]
        == "INCONSISTENT"
    )
    assert checks["natural source-map correction a^2/r^2"]["status"] == "CONSISTENT"
    assert "mhat0^2 remains the scale carrier" in report["pathA_22a_effect"]


def test_gate0b_source_assessment_does_not_cancel() -> None:
    report = z_kernel_cancellation_source_assessment()

    assert report["outcome"] == "DOES_NOT_CANCEL"
    assert report["outcome_label"] == DOES_NOT_CANCEL_SOURCE_LABEL
    assert report["sources_establish_condition"] is False
    assert report["source_establishes_shared_scalar_route"] is False
    assert report["source_establishes_pointwise_proportional_route"] is False
    assert "W_eff/full transverse profile" in report["implication_for_gate4"]
    assert "both cancellation routes" in report["implication_for_gate4"]


def test_gate0b_negative_control_is_reachable() -> None:
    w = sp.Symbol("w")

    control = classify_kernel_pair(1 + w, 1 + w**2, shared_scalar_factored=False, w=w)

    assert control["outcome"] == "DOES_NOT_CANCEL"
    assert control["z_independent_under_available_route"] is False
    assert control["kernels_proportional_for_arbitrary_Z"] is False
    assert control["proportionality_condition"] != "0"


def test_gate0b_proportional_kernels_are_valid_route_b() -> None:
    w = sp.Symbol("w")

    control = classify_kernel_pair(2 * (1 + w**2), 1 + w**2, shared_scalar_factored=False, w=w)

    assert control["outcome"] == "CANCELS"
    assert control["route_a_shared_factored_scalar"] is False
    assert control["route_b_pointwise_proportional_kernels"] is True
    assert control["z_independent_under_available_route"] is True
    assert control["kernels_proportional_for_arbitrary_Z"] is True
    assert control["proportionality_condition"] == "0"


def test_gate0b_positive_common_scalar_control_is_reachable() -> None:
    report = z_kernel_cancellation_source_assessment()

    assert report["factored_control"]["outcome"] == "CANCELS"
    assert report["factored_control"]["ratio_after_cancel"] == "K_stress/K_source"


def test_gate0_payload_carries_required_outcomes_and_residuals() -> None:
    payload = gate0_report_payload()

    assert payload["gate0_outcomes"]["0a"] == "MHAT_DIMENSIONFUL_CONFIRMED"
    assert payload["gate0_outcomes"]["0b"] == "DOES_NOT_CANCEL"
    assert payload["gate0_outcomes"]["0b_label"] == DOES_NOT_CANCEL_SOURCE_LABEL
    assert "z_kernel_cancellation_source_assessment" in payload
    assert payload["target_blindness"]["status"] == "TARGET_BLIND_PASS"
    assert payload["residuals"]


def test_gate0_target_blindness_guard_scans_new_sources() -> None:
    root = Path(__file__).resolve().parents[1]
    paths = [
        root / "src/stage1_solver/patha22b_gate0.py",
        root / "tools/pathA_22b_gate0_crosscheck.wl",
    ]

    result = target_blindness_guard(paths)

    assert result["status"] == "TARGET_BLIND_PASS"
    assert result["hits"] == []
