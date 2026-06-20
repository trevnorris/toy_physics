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
