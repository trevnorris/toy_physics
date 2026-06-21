from __future__ import annotations

import math

import numpy as np

from stage1_solver import patha_c0f_globalization_probe as c0f


def test_c0f_zero_epsilon_summary_does_not_conflate_linf_pass_with_solver_convergence() -> None:
    row = {
        "epsilon_attempts": [
            {
                "core_epsilon": 0.0,
                "k1_radius_epsilon": 0.0,
                "converged": False,
                "iterations": 20,
                "line_search_alphas": [1.0, 0.25],
            }
        ],
        "stage_converged": False,
    }

    zero = c0f._final_zero_epsilon_attempt(row)
    alphas = c0f._accepted_alphas(zero)

    assert zero is not None
    assert not zero["converged"]
    assert min(alphas) == 0.25
    assert c0f._line_search_backtrack_count(alphas, shrink=0.5) == 2


def test_c0f_fold_growth_requires_monotone_threshold_growth() -> None:
    no_growth = c0f._fold_growth_from_values([2.0, 3.0, 4.0], threshold=10.0)
    growth = c0f._fold_growth_from_values([1.0, 4.0, 12.0], threshold=10.0)
    nonmonotone = c0f._fold_growth_from_values([1.0, 12.0, 11.0], threshold=10.0)

    assert no_growth["status"] == "MEASURED"
    assert not no_growth["growth_condition"]
    assert math.isclose(no_growth["growth_factor"], 2.0)
    assert growth["growth_condition"]
    assert not nonmonotone["growth_condition"]


def test_c0f_verdict_wall_was_config_requires_target_and_bounded_fold() -> None:
    config = c0f.C0fConfig()
    per_tau_rows = [
        {"tau": 0.03, "accepted_default_success": True, "b2c_linf_pass": True},
        {"tau": 0.029, "accepted_default_success": True, "b2c_linf_pass": True},
    ]
    bounded_fold = {
        "status": "MEASURED",
        "call": "FOLD_RISK",
        "growth_condition": False,
        "growth_factor": 1.5,
    }

    verdict, support, _next = c0f._determine_c0f_verdict(
        per_tau_rows=per_tau_rows,
        fold=bounded_fold,
        merit={"status": "SKIPPED"},
        config=config,
    )

    assert verdict == "DIAGNOSTIC_INCOMPLETE"
    assert support["accepted_default_crawl_through_target"] is False

    per_tau_rows.append(
        {"tau": 0.028, "accepted_default_success": True, "b2c_linf_pass": True}
    )
    verdict, support, _next = c0f._determine_c0f_verdict(
        per_tau_rows=per_tau_rows,
        fold=bounded_fold,
        merit={"status": "SKIPPED"},
        config=config,
    )

    assert verdict == "WALL_WAS_CONFIG"
    assert support["accepted_default_crawl_through_tau"] == 0.028


def test_c0f_mixed_control_math_reports_epsilon_squared_fraction() -> None:
    basis = np.eye(3, 1)
    gradient = np.asarray([1.0, 0.0, 0.0])
    transverse = np.asarray([0.0, 1.0, 0.0])
    eps = 0.1
    mixed = c0f._unit(gradient + eps * transverse)
    p_g = float(np.linalg.norm(basis.T @ mixed) ** 2)

    assert np.isclose(1.0 - p_g, eps * eps / (1.0 + eps * eps))
