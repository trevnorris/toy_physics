from __future__ import annotations

from dataclasses import replace
import math
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.sparse import csc_matrix
import torch

from stage1_solver import patha_c0_conditioning_spike as c0
from stage1_solver.backend import configure_backend
from stage1_solver.config import BackendConfig, PreconditionerConfig
from stage1_solver.coupled_branch import (
    _create_branch_grid,
    branch_boundary_conditions,
    initial_closed_branch_state,
    pack_closed_coupled_fields,
    patha_closed_branch_residual,
    unpack_closed_coupled_fields,
)
from stage1_solver.patha_b2a_bdg import FROZEN_A, FROZEN_L, frozen_s_sigma_spec
from stage1_solver.patha_static_balance import resolve_s_sigma


def _small_branch(tau: float = 0.03):
    branch = c0.c0_frozen_branch(tau=tau, grid=(4, 4))
    return replace(
        branch,
        r_max=1.4,
        newton=replace(
            branch.newton,
            preconditioner=PreconditionerConfig(
                type="colored_sparse_jacobian_lu",
                side="left",
                rebuild_policy="every_newton_step",
                stencil_radius=1,
                color_separation=3,
                factorization="splu",
                diagonal_shift=0.0,
            ),
        ),
    )


def test_c0_zero_epsilon_preconditioner_residual_matches_physical_residual() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    boundaries = branch_boundary_conditions(branch)
    provider = resolve_s_sigma(frozen_s_sigma_spec(0.03))
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    eos_K = float(branch.continuation_K_values[-1])

    physical = patha_closed_branch_residual(
        state,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=provider,
    )
    conditioned = c0.c0_preconditioner_residual(
        state,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=provider,
        aid=c0.C0AidParameters(core_epsilon=0.0, k1_radius_epsilon=0.0),
    )

    assert torch.max(torch.abs(physical - conditioned)).item() <= 1.0e-12


def test_c0_jacobi_scaling_solves_equivalent_linear_system() -> None:
    matrix = csc_matrix(
        np.asarray(
            [
                [1.0e-3, 2.0, 0.0],
                [0.0, -3.0e2, 4.0],
                [5.0e-2, 0.0, 6.0e1],
            ],
            dtype=np.float64,
        )
    )
    rhs = np.asarray([1.0, -2.0, 0.5], dtype=np.float64)
    aid = c0.C0AidParameters(use_jacobi_scaling=True)
    row_scale, col_scale = c0.jacobi_row_col_scales(matrix, aid=aid)

    direct = np.linalg.solve(matrix.toarray(), rhs)
    scaled_matrix = (np.diag(row_scale) @ matrix.toarray() @ np.diag(col_scale))
    scaled_solution = np.linalg.solve(scaled_matrix, row_scale * rhs)
    recovered = col_scale * scaled_solution

    assert np.all(row_scale > 0.0)
    assert np.all(col_scale > 0.0)
    assert np.allclose(recovered, direct, rtol=1.0e-11, atol=1.0e-11)


def test_c0_tau_depth_branch_preserves_frozen_physics() -> None:
    shallow = c0.c0_frozen_branch(tau=0.03, grid=(4, 4))
    deep = c0.c0_frozen_branch(tau=0.028, grid=(4, 4))

    for branch in (shallow, deep):
        assert branch.w_max == FROZEN_L
        assert branch.r_mouth == FROZEN_A
        assert branch.r_exit == FROZEN_A
        assert branch.geometry_profile == "cubic_smoothstep"
        assert branch.matter_mouth_boundary == shallow.matter_mouth_boundary
        assert branch.a0_mouth_boundary == shallow.a0_mouth_boundary

    shallow_spec = frozen_s_sigma_spec(0.03).to_dict()
    deep_spec = frozen_s_sigma_spec(0.028).to_dict()
    assert shallow_spec["family"] == deep_spec["family"] == "homogeneous_isotropic_hooke_v1"
    assert shallow_spec["parameters"]["a"] == deep_spec["parameters"]["a"] == FROZEN_A
    assert shallow_spec["parameters"]["w_max"] == deep_spec["parameters"]["w_max"] == FROZEN_L
    assert shallow_spec["parameters"]["tau"] != deep_spec["parameters"]["tau"]


def test_c0_below_floor_disables_existing_b2c_predictor() -> None:
    config = c0.C0Config(prefer_existing_b2c_background_predictor=True)

    assert c0._may_use_existing_b2c(c0.PRIOR_TAU_FLOOR + 1.0e-4, config)
    assert not c0._may_use_existing_b2c(c0.PRIOR_TAU_FLOOR - 1.0e-4, config)


def test_c0_true_sigma_triplet_uses_dense_svd_and_closed_layout() -> None:
    branch = _small_branch()
    dtype = configure_backend(BackendConfig())
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    dims = c0._closed_layout_dimensions(grid)
    values = np.linspace(2.0, 4.0, dims["state_size"], dtype=np.float64)
    values[dims["state_size"] - 1] = 1.0e-6
    matrix = csc_matrix(np.diag(values))

    diagnostic = c0._true_singular_triplet_diagnostics(
        matrix,
        grid=grid,
        dense_max_size=10_000,
        triplet_count=3,
        triplet_residual_rtol=1.0e-10,
    )

    assert diagnostic["status"] == "MEASURED"
    assert diagnostic["sigma_method"] == "dense_svd"
    assert diagnostic["sigma_min"] == np.min(values)
    assert diagnostic["triplet_residual_rel"] <= 1.0e-10
    assert diagnostic["v_min_energy_fractions"]["mu"] > 0.999
    assert diagnostic["u_min_energy_fractions"]["mass_constraint_row"] > 0.999
    assert diagnostic["cond_ratio"] < 0.1


def test_c0g_complement_svd_does_not_report_planted_gauge_zero(tmp_path) -> None:
    grid = SimpleNamespace(
        spec=SimpleNamespace(nr=1, nw=1),
        r_centers=torch.as_tensor([0.5], dtype=torch.float64),
        w_centers=torch.as_tensor([0.25], dtype=torch.float64),
    )
    matrix = csc_matrix(np.diag([0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]))
    gauge = np.zeros((7, 1), dtype=np.float64)
    gauge[0, 0] = 1.0

    result = c0._c0g_dense_complement_svd(
        scaled_matrix=matrix,
        gauge_matrix=gauge,
        col_scale=np.ones(7, dtype=np.float64),
        grid=grid,
        config=c0.C0gConfig(near_null_mode_count=2),
        mode_artifact_path=tmp_path / "modes.npz",
    )

    assert result["status"] == "MEASURED"
    assert result["forbidden_construction_not_used"] == "SVD((I-P_G)J)"
    assert np.isclose(result["sigma_min"], 2.0)
    assert result["sigma_min"] > 0.0


def test_c0g_gauge_fixed_step_solves_only_complement_coordinates() -> None:
    matrix = csc_matrix(np.diag([0.0, 2.0, 3.0]))
    gauge = np.asarray([[1.0], [0.0], [0.0]], dtype=np.float64)
    rhs = np.asarray([1.0, 4.0, 9.0], dtype=np.float64)

    step, info = c0._solve_scaled_gauge_complement_step(
        scaled_matrix=matrix,
        rhs_scaled=rhs,
        gauge_matrix=gauge,
        col_scale=np.ones(3, dtype=np.float64),
        gradient_rank_rtol=1.0e-12,
        lstsq_rcond=1.0e-12,
    )

    assert info["gauge_rank"] == 1
    assert info["method"] == "iterative_lsmr_scaled_J_times_Q_perp"
    assert abs(step[0]) <= 1.0e-14
    assert np.allclose(step[1:], [2.0, 3.0], rtol=1.0e-12, atol=1.0e-12)


def test_c0g_premise_gate_thresholds_are_executable() -> None:
    config = c0.C0gConfig()

    assert c0._c0g_premise_gate_call(9.9e-3, config) == "PREMISE_FAILED"
    assert c0._c0g_premise_gate_call(5.0e-2, config) == "PREMISE_GRAY"
    assert c0._c0g_premise_gate_call(1.1e-1, config) == "PREMISE_HOLDS"


def test_c0g_provenance_prefers_real_crawl_config_field() -> None:
    value, source = c0._c0g_prefer_existing_flag_from_c0f2_payload(
        {
            "config": {"prefer_existing_b2c_background_predictor": True},
            "crawl_config": {"prefer_existing_b2c_background_predictor": False},
            "scope_guard": {"prefer_existing_b2c_background_predictor": True},
        }
    )

    assert value is False
    assert source == "crawl_config.prefer_existing_b2c_background_predictor"


def test_c0g_provenance_fails_if_real_prefer_existing_field_absent() -> None:
    with pytest.raises(KeyError):
        c0._c0g_prefer_existing_flag_from_c0f2_payload({"config": {}})


def test_c0g_scaled_lane_overlap_keeps_branch_direction() -> None:
    grid = SimpleNamespace(spec=SimpleNamespace(nr=1, nw=1))
    origin = np.zeros(7, dtype=np.float64)
    reference = np.ones(7, dtype=np.float64)
    candidate = -reference

    overlap = c0._c0g_scaled_lane_overlaps(
        candidate=candidate,
        reference=reference,
        origin=origin,
        col_scale=np.ones(7, dtype=np.float64),
        grid=grid,
    )

    assert overlap["aggregate_overlap_call"] == "POST_FOLD_BRANCH"
    assert all(row["scaled_overlap"] < 0.5 for row in overlap["rows"])


def test_c0g_step8_does_not_force_clean_verdict_on_gray_overlap() -> None:
    steps0_3 = {
        "step0": {"premise_gate": {"call": "PREMISE_HOLDS"}},
        "state_results": [
            {"ftau": {"status": "MEASURED", "call": "FOLD_SUPPORT"}},
        ],
        "sigma_min_squared_fit": {
            "status": "MEASURED",
            "fit_reliability": "RELIABLE",
            "call": "LINEAR_MONOTONE_FOLD_SUPPORT",
        },
    }
    steps4_6_7 = {
        "step4_mach": {"sonic_context_label": "NO_SONIC"},
        "step6_bordered_conditioning": {"call": "FOLD_SUPPORT"},
    }
    step5 = {"branch_overlap": {"aggregate_overlap_call": "INCONCLUSIVE"}}

    verdict = c0._c0g_determine_step8_verdict(
        steps0_3=steps0_3,
        steps4_6_7=steps4_6_7,
        step5=step5,
    )

    assert verdict["verdict"] == "MIXED/INCONCLUSIVE"
    assert "Step-5 branch overlap" in verdict["primary_evidence"]["gray_items"]


def test_c0g_step8_allows_fold_when_primary_clean_and_step6_confirms() -> None:
    steps0_3 = {
        "step0": {"premise_gate": {"call": "PREMISE_HOLDS"}},
        "state_results": [
            {"ftau": {"status": "MEASURED", "call": "FOLD_SUPPORT"}},
        ],
        "sigma_min_squared_fit": {
            "status": "MEASURED",
            "fit_reliability": "RELIABLE",
            "call": "LINEAR_MONOTONE_FOLD_SUPPORT",
        },
    }
    steps4_6_7 = {
        "step4_mach": {"sonic_context_label": "NO_SONIC"},
        "step6_bordered_conditioning": {"call": "FOLD_SUPPORT"},
    }
    step5 = {"branch_overlap": {"aggregate_overlap_call": "NOT_COMPUTED_NO_PROGRESS"}}

    verdict = c0._c0g_determine_step8_verdict(
        steps0_3=steps0_3,
        steps4_6_7=steps4_6_7,
        step5=step5,
    )

    assert verdict["verdict"] == "FOLD_CONFIRMED"
    assert verdict["sub_label"] == "NON_SONIC_FOLD"


def test_c0g_b2_gate_confirms_only_relative_trends() -> None:
    crawl = {
        "config": {"max_tau_backtracks": 5, "max_newton_iters_override": None},
        "tau_attempts": [
            {"target_tau": 0.03, "final_physical_converged": True},
            {"target_tau": 0.0295, "final_physical_converged": True},
            {"target_tau": 0.02925, "final_physical_converged": True},
            {"target_tau": 0.029125, "final_physical_converged": True},
            {"target_tau": 0.029124, "final_physical_converged": True},
            {"target_tau": 0.0291225, "final_physical_converged": False, "backtrack_index": 0},
            {"target_tau": 0.02912325, "final_physical_converged": False, "backtrack_index": 1},
        ]
    }
    b2 = {
        "state_results": [
            {
                "tau": tau,
                "crawl_final_physical_converged": True,
                "ftau": {
                    "status": "MEASURED",
                    "stepsize_stability_call": "STABLE",
                    "representative_cos_theta": 0.5,
                    "call": "FOLD_SUPPORT",
                },
            }
            for tau in (0.03, 0.0295, 0.02925, 0.029125)
        ],
        "sigma_min_squared_fit": {
            "status": "MEASURED",
            "call": "LINEAR_MONOTONE_FOLD_SUPPORT",
            "sigma_min_monotone_decreasing_toward_stall": True,
            "linear": {
                "r2": 0.99,
                "slope": 1.0,
                "tau_fold_zero_crossing": 0.0291233,
            },
            "rows": [],
        },
        "step6_bordered_conditioning": {
            "rows": [
                {"status": "MEASURED", "tau": 0.03, "cond_JQ_perp": 1.0e3, "cond_Jb": 5.0e2},
                {"status": "MEASURED", "tau": 0.0295, "cond_JQ_perp": 4.0e3, "cond_Jb": 5.0e2},
                {"status": "MEASURED", "tau": 0.02925, "cond_JQ_perp": 1.0e4, "cond_Jb": 5.0e2},
                {"status": "MEASURED", "tau": 0.029125, "cond_JQ_perp": 3.0e4, "cond_Jb": 5.0e2},
            ]
        },
    }

    verdict = c0._c0g_b2_gate_verdict(
        crawl=crawl,
        b2=b2,
        config=c0.C0gBuildB1B2Config(),
    )

    assert verdict["verdict"] == "FOLD_CONFIRMED"
    assert verdict["absolute_cond_bar_used"] is False
    assert verdict["full_budget_dissolve_contest"]["status"] == "PASS"
    assert verdict["fold_confirmed_components"]["bordered_conditioning_trend"]["ratio_growth_factor"] == 30.0


def test_c0g_b2_gate_stays_inconclusive_without_full_budget_dissolve_contest() -> None:
    crawl = {
        "config": {"max_tau_backtracks": 0, "max_newton_iters_override": 1},
        "tau_attempts": [
            {"target_tau": 0.03, "final_physical_converged": True},
            {"target_tau": 0.0295, "final_physical_converged": True},
            {"target_tau": 0.02925, "final_physical_converged": True},
            {"target_tau": 0.029125, "final_physical_converged": True},
            {"target_tau": 0.029124, "final_physical_converged": False, "backtrack_index": 0},
            {"target_tau": 0.0291225, "final_physical_converged": False, "backtrack_index": 0},
        ],
    }
    b2 = {
        "state_results": [
            {
                "tau": tau,
                "crawl_final_physical_converged": True,
                "ftau": {
                    "status": "MEASURED",
                    "stepsize_stability_call": "STABLE",
                    "representative_cos_theta": 0.5,
                    "call": "FOLD_SUPPORT",
                },
            }
            for tau in (0.03, 0.0295, 0.02925, 0.029125)
        ],
        "sigma_min_squared_fit": {
            "status": "MEASURED",
            "call": "LINEAR_MONOTONE_FOLD_SUPPORT",
            "sigma_min_monotone_decreasing_toward_stall": True,
            "linear": {
                "r2": 0.99,
                "slope": 1.0,
                "tau_fold_zero_crossing": 0.0291233,
            },
        },
        "step6_bordered_conditioning": {
            "rows": [
                {"status": "MEASURED", "tau": 0.03, "cond_JQ_perp": 1.0e3, "cond_Jb": 5.0e2},
                {"status": "MEASURED", "tau": 0.0295, "cond_JQ_perp": 4.0e3, "cond_Jb": 5.0e2},
                {"status": "MEASURED", "tau": 0.02925, "cond_JQ_perp": 1.0e4, "cond_Jb": 5.0e2},
                {"status": "MEASURED", "tau": 0.029125, "cond_JQ_perp": 3.0e4, "cond_Jb": 5.0e2},
            ]
        },
    }

    verdict = c0._c0g_b2_gate_verdict(
        crawl=crawl,
        b2=b2,
        config=c0.C0gBuildB1B2Config(),
    )

    assert verdict["verdict"] == "STILL_INCONCLUSIVE"
    assert verdict["full_budget_dissolve_contest"]["status"] == "FAIL"


def test_c0g_b2_dissolved_accepts_gauge_fixed_previous_source() -> None:
    crawl = {
        "config": {"max_tau_backtracks": 5, "max_newton_iters_override": None},
        "tau_attempts": [
            {
                "target_tau": 0.0291225,
                "final_physical_converged": True,
                "final_original_residual_linf": 1.0e-8,
                "used_existing_b2c": False,
                "init": {"source": "previous_c0g_gauge_fixed_converged_state"},
            },
            {
                "target_tau": 0.029122,
                "final_physical_converged": True,
                "final_original_residual_linf": 2.0e-8,
                "used_existing_b2c": False,
                "init": {"source": "previous_c0g_gauge_fixed_converged_state"},
            },
        ],
    }
    b2 = {
        "state_results": [],
        "sigma_min_squared_fit": {
            "status": "MEASURED",
            "call": "NO_LINEAR_MONOTONE_FOLD_SUPPORT",
            "sigma_min_monotone_decreasing_toward_stall": False,
            "linear": {"r2": 0.2},
        },
        "step6_bordered_conditioning": {"rows": []},
    }

    verdict = c0._c0g_b2_gate_verdict(
        crawl=crawl,
        b2=b2,
        config=c0.C0gBuildB1B2Config(),
    )

    assert verdict["verdict"] == "FOLD_DISSOLVED"
    assert verdict["fold_dissolved"]["branch_continuous_count"] == 2


def test_c0g_deep_verdict_locks_no_fold_with_admissible_finite_sigma() -> None:
    rows = []
    ladder = []
    state_results = []
    admissibility = []
    for index, tau in enumerate((0.0291225, 0.02912, 0.0290, 0.0288, 0.0285)):
        rows.append(
            {
                "target_tau": tau,
                "deepcrawl_role": "deep_crawl_attempt" if index >= 1 else "below_fold_seed",
                "final_physical_converged": True,
                "final_original_residual_linf": 1.0e-9,
                "measurement_status": "MEASURED",
            }
        )
        ladder.append(
            {
                "tau": tau,
                "converged": True,
                "min_R0": 0.79 + 0.01 * index,
                "min_rho": 7.0e-6 + 1.0e-7 * index,
                "mu": 2.04 - 0.01 * index,
            }
        )
        state_results.append({"tau": tau, "svd": {"status": "MEASURED", "sigma_min": 3.0e-4}})
        admissibility.append({"tau": tau, "status": "PASS"})

    verdict = c0._c0g_deep_verdict(
        rows=rows,
        ladder=ladder,
        state_results=state_results,
        admissibility_rows=admissibility,
        sigma_fit={"call": "FIT_UNRELIABLE_NON_MONOTONE_SIGMA_MIN", "linear": {"r2": 0.2}},
        tangent={"aggregate_call": "NO_TANGENT_REVERSAL"},
        config=c0.C0gDeepCrawlConfig(no_fold_required_new_converged=4),
    )

    assert verdict["verdict"] == "NO_FOLD_LOCKED"
    assert verdict["sigma_min_finite"] is True
    assert verdict["all_converged_below_admissible"] is True


def test_c0g_deep_verdict_does_not_lock_if_admissibility_fails() -> None:
    rows = [
        {
            "target_tau": 0.0290 - index * 1.0e-4,
            "deepcrawl_role": "deep_crawl_attempt",
            "final_physical_converged": True,
            "final_original_residual_linf": 1.0e-9,
            "measurement_status": "MEASURED",
        }
        for index in range(5)
    ]
    ladder = [
        {
            "tau": row["target_tau"],
            "converged": True,
            "min_R0": 0.8 + 0.01 * index,
            "min_rho": 7.0e-6 + 1.0e-7 * index,
            "mu": 2.0 - 0.01 * index,
        }
        for index, row in enumerate(rows)
    ]
    state_results = [
        {"tau": row["target_tau"], "svd": {"status": "MEASURED", "sigma_min": 2.0e-4}}
        for row in rows
    ]
    admissibility = [
        {"tau": row["target_tau"], "status": "FAIL" if index == 2 else "PASS"}
        for index, row in enumerate(rows)
    ]

    verdict = c0._c0g_deep_verdict(
        rows=rows,
        ladder=ladder,
        state_results=state_results,
        admissibility_rows=admissibility,
        sigma_fit={"call": "FIT_UNRELIABLE_NON_MONOTONE_SIGMA_MIN", "linear": {"r2": 0.2}},
        tangent={"aggregate_call": "NO_TANGENT_REVERSAL"},
        config=c0.C0gDeepCrawlConfig(no_fold_required_new_converged=5),
    )

    assert verdict["verdict"] == "STILL_GOING"
    assert verdict["all_converged_below_admissible"] is False


def test_c0_verdict_is_incomplete_when_sigma_or_fold_not_measured() -> None:
    result = {
        "tau_attempts": [
            {
                "target_tau": c0.PRIOR_TAU_FLOOR - 1.0e-4,
                "final_physical_converged": False,
            }
        ],
        "sigma_diagnostics": {"status": "NOT_MEASURED"},
        "fold_diagnostic": {"status": "NOT_MEASURED"},
    }

    verdict, support, _next_step = c0._determine_verdict(result)

    assert verdict == "DIAGNOSTIC_INCOMPLETE"
    assert support["reason"] == "required_sigma_or_fold_evidence_not_measured"


def test_c0_tau_attempt_schema_catches_below_floor_short_circuit() -> None:
    rows = [
        {
            "target_tau": c0.PRIOR_TAU_FLOOR - 1.0e-4,
            "delta_tau": -1.0e-4,
            "backtrack_index": 0,
            "init": {"source": "existing_b2c_background_target_predictor"},
            "used_existing_b2c": True,
            "start_from_tau": c0.PRIOR_TAU_FLOOR,
            "epsilon_attempts": [
                {"epsilon": 8.0e-2},
                {"epsilon": 2.0e-2},
                {"epsilon": 0.0},
            ],
        }
    ]

    errors = c0._validate_tau_attempts_schema(rows)

    assert any("used existing B2c" in error for error in errors)
    assert any("did not warm-start" in error for error in errors)


def test_c0_production_verdict_requires_persistent_backtracked_full_budget_failure() -> None:
    result = {
        "config": {"max_tau_backtracks": 0, "max_newton_iters_override": 2},
        "tau_attempts": [
            {
                "target_tau": c0.PRIOR_TAU_FLOOR - 1.0e-5,
                "backtrack_index": 0,
                "final_physical_converged": False,
            }
        ],
        "sigma_diagnostics": {
            "status": "MEASURED",
            "deepest": {
                "v_min_energy_fractions": {"field": 1.0, "R0": 0.0, "mu": 0.0},
                "u_min_energy_fractions": {
                    "pde_rows": 1.0,
                    "wall_rows": 0.0,
                    "mass_constraint_row": 0.0,
                },
                "cond_ratio": 1.0,
                "sigma_min": 1.0e-12,
                "sigma_method": "dense_svd",
            },
        },
        "fold_diagnostic": {
            "status": "MEASURED",
            "call": "NO_TRACKED_COMPLEMENT_CROSSING",
        },
    }

    verdict, support, _next_step = c0._determine_verdict(result)

    assert verdict == "DIAGNOSTIC_INCOMPLETE"
    assert not support["crawl_persistent_failure_evidence"]
    assert not support["attempted_backtracking"]
    assert not support["full_newton_budget"]


def test_c0c_phase_generator_uses_confirmed_closed_layout() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    assert mu is not None
    imag = torch.linspace(
        0.1,
        0.4,
        fields.psi_imag.numel(),
        dtype=state.dtype,
        device=state.device,
    ).reshape_as(fields.psi_imag)
    state = pack_closed_coupled_fields(
        type(fields)(
            psi_real=fields.psi_real,
            psi_imag=imag,
            a0=fields.a0,
            ar=fields.ar,
            aw=fields.aw,
            r0=fields.r0,
        ),
        mu,
    )

    generators = c0._c0c_generators_for_state(state, grid)
    phase = next(generator for generator in generators if generator.name == "phase")
    n = grid.spec.nr * grid.spec.nw

    assert phase.symmetry_status == "EXACT_SYMMETRY"
    assert np.allclose(phase.vector[:n], -imag.detach().cpu().numpy().reshape(-1))
    assert np.allclose(phase.vector[n : 2 * n], fields.psi_real.detach().cpu().numpy().reshape(-1))
    assert np.allclose(phase.vector[2 * n :], 0.0)


def test_c0c_planted_generator_overlap_reports_one() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    generators = c0._c0c_generators_for_state(state, grid)
    phase = next(generator for generator in generators if generator.name == "phase")
    planted, _norm = c0._unit_vector(phase.vector)

    modes, span_rank = c0._c0c_overlap_diagnostics(
        right_vectors=planted.reshape(1, -1),
        singular_values=[1.0e-12],
        generators=generators,
        grid=grid,
    )

    assert span_rank >= 1
    assert modes[0]["overlaps"]["phase"] >= 1.0 - 1.0e-12
    assert modes[0]["unexplained_residual_fraction"] <= 1.0e-12


def test_c0c_whole_cluster_verdict_is_mixed_when_only_one_mode_is_phase() -> None:
    config = c0.C0cConfig(cluster_mode_count=2)
    modes = [
        {
            "mode_index": 0,
            "sigma": 1.0e-14,
            "overlaps": {"phase": 0.95, "translation_r": 0.0},
            "unexplained_residual_fraction": 0.05,
            "lane_energy_fractions": {"psi_real": 0.5, "psi_imag": 0.5},
        },
        {
            "mode_index": 1,
            "sigma": 2.0e-10,
            "overlaps": {"phase": 0.1, "translation_r": 0.2},
            "unexplained_residual_fraction": 0.96,
            "lane_energy_fractions": {"a0": 0.9, "ar": 0.1},
        },
    ]
    generators = [
        {
            "name": "phase",
            "classification": "GAUGE_PHASE",
            "symmetry_status": "EXACT_SYMMETRY",
        },
        {
            "name": "translation_r",
            "classification": "TRANSLATION",
            "symmetry_status": "PROBE",
        },
    ]
    annihilation = [
        {"generator": "phase", "status": "MEASURED", "null_gate_pass": True},
        {"generator": "translation_r", "status": "MEASURED", "null_gate_pass": False},
    ]

    classified, support, verdict, recommended = c0._classify_c0c_modes(
        modes=modes,
        generators=generators,
        annihilation_rows=annihilation,
        config=config,
    )

    assert verdict == "MIXED"
    assert classified[0]["classification"] == "GAUGE_PHASE"
    assert classified[1]["classification"] == "UNEXPLAINED_STIFFNESS"
    assert support["explained_mode_count"] == 1
    assert "pin or deflate" in recommended


def test_c0d_gauge_subspace_controls_are_anti_hardcode() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    config = c0.C0dConfig()

    subspace = c0._c0d_build_gauge_subspace(grid, branch, config)
    controls = c0._c0d_control_diagnostics(subspace, grid, config)

    assert subspace["dim_g"] > 1
    assert subspace["dim_g"] <= grid.spec.nr * grid.spec.nw
    assert controls["status"] == "PASS"
    assert controls["positive_control"]["p_g_fraction"] >= 1.0 - 1.0e-10
    assert controls["positive_control"]["non_a_remainder"] <= 1.0e-12
    assert controls["negative_control"]["p_g_fraction"] <= 1.0e-10
    assert controls["negative_control"]["non_a_remainder"] <= 1.0e-12


def test_c0d_weighted_gauge_residual_matches_operator_definition() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    r_scaled = grid.r_centers[:, None] / float(grid.spec.r_max)
    w_scaled = (grid.w_centers[None, :] - float(grid.spec.w_min)) / (
        float(grid.spec.w_max) - float(grid.spec.w_min)
    )
    lambda_field = torch.sin(math.pi * r_scaled) * (1.0 + 0.25 * w_scaled)
    vector = c0._c0d_scalar_gradient_vector(lambda_field, grid)
    ar, aw = c0._c0d_extract_spatial_a_fields(vector, grid)

    divergence = c0.axisymmetric_vector_divergence(ar, aw, grid)
    z = c0.localization_weight_torch(grid.w_centers, branch)[None, :]
    weighted = z * divergence
    manual_weighted_gradient = torch.cat(
        [
            c0.tensor_center_gradient_r(weighted, grid).reshape(-1),
            c0.tensor_center_gradient_w(weighted, grid).reshape(-1),
        ]
    )
    a_norm = c0._c0d_spatial_a_norm(vector, grid)
    metrics = c0._c0d_spatial_a_metrics(vector, grid, branch)

    assert np.isclose(
        metrics["raw_divergence"],
        float(torch.linalg.vector_norm(divergence).item()) / a_norm,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    assert np.isclose(
        metrics["weighted_gauge_residual"],
        float(torch.linalg.vector_norm(manual_weighted_gradient).item()) / a_norm,
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_c0e_coupled_generator_uses_real_lane_charge_formula() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    assert mu is not None
    imag = torch.linspace(
        0.2,
        0.5,
        fields.psi_imag.numel(),
        dtype=state.dtype,
        device=state.device,
    ).reshape_as(fields.psi_imag)
    state = pack_closed_coupled_fields(
        type(fields)(
            psi_real=fields.psi_real,
            psi_imag=imag,
            a0=fields.a0,
            ar=fields.ar,
            aw=fields.aw,
            r0=fields.r0,
        ),
        mu,
    )
    fields, _mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    r_scaled = grid.r_centers[:, None] / float(grid.spec.r_max)
    w_scaled = (grid.w_centers[None, :] - float(grid.spec.w_min)) / (
        float(grid.spec.w_max) - float(grid.spec.w_min)
    )
    lambda_field = r_scaled * (1.0 - w_scaled)
    vector = c0._c0e_coupled_gauge_vector(lambda_field, state, grid, branch)
    lanes = c0._closed_lane_slices(grid)
    alpha = branch.gauge_charge / branch.hbar

    assert np.allclose(
        vector[lanes["psi_real"][0] : lanes["psi_real"][1]],
        (-alpha * lambda_field * fields.psi_imag).detach().cpu().numpy().reshape(-1),
    )
    assert np.allclose(
        vector[lanes["psi_imag"][0] : lanes["psi_imag"][1]],
        (alpha * lambda_field * fields.psi_real).detach().cpu().numpy().reshape(-1),
    )
    assert np.allclose(vector[lanes["a0"][0] : lanes["a0"][1]], 0.0)
    assert np.allclose(
        vector[lanes["ar"][0] : lanes["ar"][1]],
        c0.tensor_center_gradient_r(lambda_field, grid).detach().cpu().numpy().reshape(-1),
    )
    assert np.allclose(
        vector[lanes["aw"][0] : lanes["aw"][1]],
        c0.tensor_center_gradient_w(lambda_field, grid).detach().cpu().numpy().reshape(-1),
    )


def test_c0e_curl_controls_are_non_circular_and_separated() -> None:
    dtype = configure_backend(BackendConfig())
    branch = replace(c0.c0_frozen_branch(tau=0.03, grid=(12, 12)), r_max=1.4)
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    config = c0.C0eConfig()
    a_subspace = c0._c0d_build_gauge_subspace(
        grid,
        branch,
        c0.C0dConfig(
            gradient_rank_rtol=config.gradient_rank_rtol,
            harmonic_weighted_divergence_rtol=config.harmonic_weighted_divergence_rtol,
        ),
    )
    coupled = c0._c0e_build_coupled_subspace(state, grid, branch, config)
    controls = c0._c0e_control_diagnostics(
        a_subspace=a_subspace,
        coupled_subspace=coupled,
        state=state,
        grid=grid,
        branch=branch,
    )

    assert controls["status"] == "PASS"
    assert controls["bands"]["unweighted"]["status"] == "SEPARATED"
    assert controls["bands"]["cell_volume_weighted"]["status"] == "SEPARATED"
    assert controls["bands"]["unweighted"]["separation_log10"] > 10.0
    assert controls["bands"]["cell_volume_weighted"]["separation_log10"] > 10.0
    positives = [row for row in controls["controls"] if row["family"] == "positive_gradient"]
    negatives = [
        row for row in controls["controls"] if row["family"] == "negative_stream_function_transverse"
    ]
    assert len(positives) >= 3
    assert len(negatives) >= 3
    assert all(row["unweighted"]["curl_fraction"] < 1.0e-12 for row in positives)
    assert all(row["cell_volume_weighted"]["curl_fraction"] < 1.0e-12 for row in positives)
    assert all(row["a_only_p_g_fraction"] >= 1.0 - 1.0e-10 for row in positives)
    assert all(row["coupled_capture_fraction"] >= 1.0 - 1.0e-10 for row in positives)
    assert all(row["non_a_remainder"] <= 1.0e-12 for row in positives)
    assert all(row["unweighted"]["curl_fraction"] > 1.0 for row in negatives)
    assert all(row["cell_volume_weighted"]["curl_fraction"] > 1.0 for row in negatives)
    assert all(row["a_only_p_g_fraction"] < 0.5 for row in negatives)
    assert all(row["construction"].startswith("independent_stream_function_A") for row in negatives)
    assert all("not_random_minus_own_gradient_projection" in row["construction"] for row in negatives)
