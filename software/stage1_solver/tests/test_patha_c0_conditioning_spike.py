from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy.sparse import csc_matrix
import torch

from stage1_solver import patha_c0_conditioning_spike as c0
from stage1_solver.backend import configure_backend
from stage1_solver.config import BackendConfig, PreconditionerConfig
from stage1_solver.coupled_branch import (
    _create_branch_grid,
    branch_boundary_conditions,
    initial_closed_branch_state,
    patha_closed_branch_residual,
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
