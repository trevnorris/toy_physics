from __future__ import annotations

from stage1_solver.dimensional_check import (
    D,
    patha22a_classify_headline,
    patha22a_dimensional_skeleton_report,
    patha22a_homogeneity_checks,
    patha22a_knob_ledger,
)


def test_patha22a_planted_dimensional_mismatch_is_caught() -> None:
    target = D["G_3_spatial"] * (D["c_s"] ** 5) / ((D["a"] ** 5) * (D["c"] ** 5))
    planted_missing_a5 = target * (D["a"] ** 5)

    report = patha22a_homogeneity_checks(planted_rhs=planted_missing_a5)

    assert report["homogeneity_status"] == "HOMOGENEITY_FAILURE"
    r_norm = [
        check
        for check in report["checks"]
        if check["name"].startswith("R_norm=")
    ][0]
    assert r_norm["status"] == "INCONSISTENT"
    assert r_norm["factor_needed_to_reach_expected"] != "1"


def test_patha22a_unresolved_dimensionless_factor_is_preserved() -> None:
    ledger = patha22a_knob_ledger(include_control_free_factor=True)

    headline = patha22a_classify_headline("HOMOGENEITY_PASS", ledger)

    assert headline == "TUNABILITY_CHANNEL_PRESENT"
    assert any(item.name == "g_control_unresolved" for item in ledger)


def test_patha22a_all_headlines_are_reachable() -> None:
    non_free_ledger = [
        item
        for item in patha22a_knob_ledger()
        if not item.classification.startswith("(d)")
    ]

    assert patha22a_classify_headline("HOMOGENEITY_FAILURE", []) == "HOMOGENEITY_FAILURE"
    assert (
        patha22a_classify_headline("HOMOGENEITY_PASS", patha22a_knob_ledger())
        == "TUNABILITY_CHANNEL_PRESENT"
    )
    assert (
        patha22a_classify_headline("HOMOGENEITY_PASS", non_free_ledger)
        == "INDETERMINATE_NEEDS_FORMS"
    )


def test_patha22a_report_carries_expected_headline_and_p0_finding() -> None:
    report = patha22a_dimensional_skeleton_report()

    assert report["homogeneity"]["homogeneity_status"] == "HOMOGENEITY_PASS"
    assert (
        report["homogeneity"]["p0_dimensionless_finding"]["status"]
        == "P0_DIMENSIONLESS_AFTER_EXPLICIT_FREQUENCY_NORMALIZATION"
    )
    assert report["headline"] == "TUNABILITY_CHANNEL_PRESENT"
    assert "chi_Q / S_port" in report["class_d_free_knobs"]
    assert set(report["tunability_channels"]) >= {"chi_Q / S_port", "lambda_gamma=c/c_s"}
    assert report["tunability_channel_count_lower_bound"] >= 2


def test_patha22a_n_tower_is_derived_from_code_formula_not_asserted() -> None:
    report = patha22a_homogeneity_checks()

    checks = {check["name"]: check for check in report["checks"]}

    assert checks["N0=P^2/Delta^2 derived from code formula"]["actual"] == "L^-1 M"
    assert checks["N0=P^2/Delta^2 derived from code formula"]["expected"] == "L^-1 M"
    assert checks["N2 derived from code formula"]["actual"] == "L^-1 T^2 M"
    assert checks["N4 derived from code formula"]["actual"] == "L^-1 T^4 M"
    assert checks["faithful P0=N0/D0 before frequency normalization"]["actual"] == "T^2"
    assert checks["normalized P0=(c_s/a)^2*N0/D0 dimension"]["actual"] == "1"

    wrong_assertion = report["formula_negative_controls"]["wrong_n0_reduced_stiffness_assertion"]
    assert report["formula_negative_controls"]["status"] == "CAUGHT_WRONG_N0_ASSERTION"
    assert wrong_assertion["status"] == "INCONSISTENT"
    assert wrong_assertion["expected"] == "L^-1 T^-2 M"
    assert wrong_assertion["actual"] == "L^-1 M"
