from __future__ import annotations

from pathlib import Path

import sympy as sp

from stage1_solver.patha22b_gate4 import (
    SOURCE_BLOCKER,
    cancellation_assessment,
    classify_kernel_pair,
    comparator_controls,
    gate4_report_payload,
    gravity_single_defect_mini_lemma,
    source_kernel_derivation,
    stress_kernel_derivation,
    target_blindness_guard,
)


def test_gate4_stress_kernel_is_derived_from_reduction_density() -> None:
    stress = stress_kernel_derivation()
    checks = {check["name"]: check for check in stress["checks"]}

    assert stress["status"] == "DERIVED_TARGET_BLIND"
    assert stress["kernel"] == "chi_N(w)*rho_inf4(w)"
    assert stress["integral_label"] == "N_infty,3"
    assert checks["K_stress(w)=chi_N(w)*rho_infty,4(w)"]["status"] == "CONSISTENT"
    assert checks["integral K_stress(w) dw = N_infty,3"]["actual"] == "L^-3"
    assert "alpha_J1" in stress["alpha_lane_power_in_g_G"]


def test_gate4_source_kernel_is_blocked_before_real_cancellation() -> None:
    source = source_kernel_derivation()

    assert source["status"] == SOURCE_BLOCKER
    assert source["kernel"] is None
    assert source["kernel_ast"] is None
    assert source["source_map_kernel_exists_target_blind"] is False
    assert "no independently derived alpha_J" in source["alpha_lane_power_in_g_mhat_squared"]


def test_gate4_pair_to_scalar_lemma_retains_pair_factors_without_roots() -> None:
    lemma = gravity_single_defect_mini_lemma()
    algebra = {check["name"]: check for check in lemma["checks"]}
    dims = {check["name"]: check for check in lemma["dimensional_checks"]}

    assert lemma["status"] == "PAIR_TO_SCALAR_G_CARRIED_WITH_PRODUCTS"
    assert algebra["G_cond pair expression"]["status"] == "CONSISTENT"
    assert "sqrt" not in lemma["mini_lemma"].lower()
    assert lemma["power_ledger"]["alpha_J1"] == -1
    assert lemma["power_ledger"]["alpha_J2"] == -1
    assert dims["conditional G from stress lane"]["status"] == "CONSISTENT"
    assert dims["g_G=G*m_GNLS/(a*c_s^2)"]["actual"] == "1"


def test_gate4_negative_control_is_reachable() -> None:
    w = sp.Symbol("w")

    control = classify_kernel_pair(1 + w, 1 + w**2, shared_scalar_factored=False, w=w)

    assert control["outcome"] == "DOES_NOT_CANCEL"
    assert control["proportionality_residual"] != "0"
    assert control["route_b_pointwise_proportional_kernels"] is False


def test_gate4_mutated_real_stress_kernel_control_is_discriminating() -> None:
    controls = comparator_controls()
    mutated = controls["mutated_kernel_control"]

    assert mutated["status"] == "PASS"
    assert mutated["outcome"] == "DOES_NOT_CANCEL"
    assert mutated["proportionality_residual"] != "0"
    assert "epsilon" in mutated["proportionality_residual"]


def test_gate4_real_routes_stop_on_undefined_source_kernel() -> None:
    assessment = cancellation_assessment()

    assert assessment["overall_verdict"] == SOURCE_BLOCKER
    assert assessment["route_a"]["outcome"] == "NOT_RUN_UNDEFINED_K_SOURCE"
    assert assessment["route_b"]["outcome"] == "NOT_RUN_UNDEFINED_K_SOURCE"
    assert assessment["h_k_source"]["finding"] == "CONFIRMED"
    assert assessment["h_alpha_J"]["finding"] == "CONFIRMED_NOT_CANCELLED"
    assert assessment["step_4iii_not_run"] is True


def test_gate4_independent_provenance_objects_are_not_reused_ast() -> None:
    assessment = cancellation_assessment()
    independence = assessment["independence_check"]
    stress = assessment["stress_kernel"]
    source = assessment["source_kernel"]

    assert independence["status"] == "PASS_SEPARATE_PROVENANCE"
    assert independence["same_symbolic_object_reused"] is False
    assert stress["provenance_id"] != source["provenance_id"]
    assert stress["kernel_ast"] is not None
    assert source["kernel_ast"] is None


def test_gate4_target_blindness_guard_scans_new_sources() -> None:
    root = Path(__file__).resolve().parents[1]
    paths = [
        root / "src/stage1_solver/patha22b_gate4.py",
        root / "tools/pathA_22b_gate4_crosscheck.wl",
        root / "tests/test_patha22b_gate4.py",
    ]

    result = target_blindness_guard(paths)

    assert result["status"] == "TARGET_BLIND_PASS"
    assert result["hits"] == []


def test_gate4_payload_carries_stop_fork_verdict_and_controls() -> None:
    payload = gate4_report_payload()
    assessment = payload["cancellation_assessment"]

    assert payload["gate4_outcome"] == SOURCE_BLOCKER
    assert payload["stop_at_decisive_fork"] is True
    assert payload["target_blindness"]["status"] == "TARGET_BLIND_PASS"
    assert assessment["controls"]["negative_control"]["status"] == "PASS"
    assert assessment["controls"]["mutated_kernel_control"]["status"] == "PASS"
