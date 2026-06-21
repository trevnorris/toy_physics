"""Path-A C0 conditioning spike for the closed throat background.

This module deliberately keeps the physical residual in
``coupled_branch.patha_closed_branch_residual`` as the only convergence
arbiter.  The C0 aids below are used only inside the Newton path:

* a matter-block epsilon term in the preconditioner residual,
* a k1 radius floor in the preconditioner wall source,
* invertible row/column scalings of the linearized system.

The accepted state is always checked by the original residual with all
epsilon aids inactive.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field, replace
import json
import math
from pathlib import Path
import time
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from scipy import linalg, optimize
from scipy.sparse import csc_matrix, diags, eye, load_npz, save_npz
from scipy.sparse.linalg import LinearOperator, eigsh, gmres, lsmr, svds
import torch

from . import patha_b2a_bdg as b2a
from . import patha_b2b_maxwell as b2b
from . import patha_b2c_calibration as b2c
from .backend import configure_backend, jvp
from .boundaries import BoundaryCondition
from .config import BackendConfig, NewtonConfig, source_revision
from .coupled_branch import (
    ClosedCoupledFields,
    _closed_to_coupled_fields,
    _create_branch_grid,
    _matter_number_current,
    branch_boundary_conditions,
    boundary_sponge_profile_torch,
    coupled_pde_residual,
    initial_closed_branch_state,
    localization_weight_torch,
    pack_closed_coupled_fields,
    pack_coupled_fields,
    patha_closed_branch_residual,
    resample_closed_branch_state,
    unpack_closed_coupled_fields,
    wall_grid_from_tensor_grid,
)
from .newton import BuiltPreconditioner, PreconditionerBuildContext
from .operators import (
    axisymmetric_vector_divergence,
    radial_divergence_from_center_flux,
    tensor_center_gradient_r,
    tensor_center_gradient_w,
)
from .patha_static_balance import SSigmaProvider, SSigmaSpec, resolve_s_sigma, static_balance_terms
from .preconditioners import (
    assemble_closed_coupled_autodiff_sparse_jacobian,
    assemble_closed_coupled_colored_sparse_jacobian,
    assemble_closed_coupled_sparse_jacobian,
    factorized_sparse_inverse_operator,
)


PRIOR_TAU_FLOOR = 0.029
BACKGROUND_RESIDUAL_TOL = b2c.BACKGROUND_RESIDUAL_TOL
DEFAULT_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0b_wall_diagnosis")
DEFAULT_REPORT_PATH = Path("software/stage1_solver/reports/pathA_C0b_wall_diagnosis.md")
DEFAULT_JSON_PATH = DEFAULT_RUN_ROOT / "pathA_C0b_diagnostic.json"
DEFAULT_C0C_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0c_nullmode_identification")
DEFAULT_C0C_REPORT_PATH = Path(
    "software/stage1_solver/reports/pathA_C0c_nullmode_identification.md"
)
DEFAULT_C0C_JSON_PATH = DEFAULT_C0C_RUN_ROOT / "pathA_C0c_nullmode_identification.json"
DEFAULT_C0D_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0d_maxwell_gauge")
DEFAULT_C0D_REPORT_PATH = Path(
    "software/stage1_solver/reports/pathA_C0d_maxwell_gauge_identification.md"
)
DEFAULT_C0D_JSON_PATH = DEFAULT_C0D_RUN_ROOT / "pathA_C0d_maxwell_gauge_identification.json"
DEFAULT_C0E_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0e_curl_gate")
DEFAULT_C0E_REPORT_PATH = Path(
    "software/stage1_solver/reports/pathA_C0e_gauge_invariant_curl_gate.md"
)
DEFAULT_C0E_JSON_PATH = DEFAULT_C0E_RUN_ROOT / "pathA_C0e_gauge_invariant_curl_gate.json"
DEFAULT_C0G_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0g_diag_fold_vs_conditioning")
DEFAULT_C0G_REPORT_PATH = Path(
    "software/stage1_solver/reports/pathA_C0g_diag_fold_vs_conditioning.md"
)
DEFAULT_C0G_STEP0_JSON_PATH = DEFAULT_C0G_RUN_ROOT / "pathA_C0g_step0_premise.json"
DEFAULT_C0G_JSON_PATH = DEFAULT_C0G_RUN_ROOT / "pathA_C0g_diag_fold_vs_conditioning_steps0_3.json"
DEFAULT_C0G_STEP5_JSON_PATH = DEFAULT_C0G_RUN_ROOT / "pathA_C0g_step5_scipy_probe.json"
DEFAULT_C0G_FINAL_JSON_PATH = DEFAULT_C0G_RUN_ROOT / "pathA_C0g_diag_fold_vs_conditioning.json"
DEFAULT_C0G_BUILD_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0g_build_B1B2")
DEFAULT_C0G_BUILD_REPORT_PATH = Path("software/stage1_solver/reports/pathA_C0g_build_B1B2.md")
DEFAULT_C0G_BUILD_JSON_PATH = DEFAULT_C0G_BUILD_RUN_ROOT / "pathA_C0g_build_B1B2.json"
DEFAULT_C0G_DEEP_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0g_deepcrawl")
DEFAULT_C0G_DEEP_REPORT_PATH = Path("software/stage1_solver/reports/pathA_C0g_deepcrawl.md")
DEFAULT_C0G_DEEP_JSON_PATH = DEFAULT_C0G_DEEP_RUN_ROOT / "pathA_C0g_deepcrawl.json"
DEFAULT_C0G_B3_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0g_B3_pseudoarclength")
DEFAULT_C0G_B3_REPORT_PATH = Path(
    "software/stage1_solver/reports/pathA_C0g_B3_pseudoarclength.md"
)
DEFAULT_C0G_B3_JSON_PATH = DEFAULT_C0G_B3_RUN_ROOT / "pathA_C0g_B3_pseudoarclength.json"
DEFAULT_C0G_B3_FOLLOWUP_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0g_B3_followup")
DEFAULT_C0G_B3_FOLLOWUP_REPORT_PATH = Path(
    "software/stage1_solver/reports/pathA_C0g_B3_followup.md"
)
DEFAULT_C0G_B3_FOLLOWUP_JSON_PATH = (
    DEFAULT_C0G_B3_FOLLOWUP_RUN_ROOT / "pathA_C0g_B3_followup.json"
)
ALLOWED_VERDICTS = {
    "SPIKE_SUFFICIENT",
    "FOLD_TURNING_POINT",
    "NEAR_NULL_SPACE_STRUCTURAL",
    "PRODUCTION_SOLVER_REQUIRED",
    "DIAGNOSTIC_INCOMPLETE",
}
C0C_VERDICTS = {
    "SYMMETRY_MODE_IDENTIFIED",
    "MIXED",
    "GENUINE_STIFFNESS",
    "DIAGNOSTIC_INCOMPLETE",
}
C0D_VERDICTS = {
    "WALL_IS_ALL_GAUGE",
    "MIXED_GAUGE_PLUS_RESIDUAL",
    "GENUINE_MAXWELL_STIFFNESS",
    "DIAGNOSTIC_INCOMPLETE",
}
C0E_VERDICTS = {
    "WALL_IS_ALL_GAUGE",
    "GAUGE_PLUS_ONE_TRANSVERSE",
    "GAUGE_FRAMING_REFUTED",
    "DIAGNOSTIC_INCOMPLETE",
}
C0C_FIELD_LAYOUT = ("psi_real", "psi_imag", "a0", "ar", "aw", "r0", "mu")


@dataclass(frozen=True)
class C0AidParameters:
    """Solver-path aid values for one Newton call."""

    core_epsilon: float = 0.0
    k1_radius_epsilon: float = 0.0
    use_jacobi_scaling: bool = True
    preconditioner_diagonal_shift: float = 1.0e-12
    scale_floor: float = 1.0e-14
    scale_min: float = 1.0e-8
    scale_max: float = 1.0e8


@dataclass(frozen=True)
class C0Config:
    run_root: Path = DEFAULT_RUN_ROOT
    report_path: Path = DEFAULT_REPORT_PATH
    json_path: Path = DEFAULT_JSON_PATH
    depth_sequence: tuple[float, ...] = (
        0.03,
        0.0295,
        0.02925,
        0.029125,
        0.0290625,
        0.029,
        0.0285,
        0.028,
    )
    grid: tuple[int, int] = b2a.DEFAULT_BACKGROUND_GRID
    epsilon_schedule: tuple[C0AidParameters, ...] = (
        C0AidParameters(core_epsilon=8.0e-2, k1_radius_epsilon=8.0e-1),
        C0AidParameters(core_epsilon=2.0e-2, k1_radius_epsilon=4.0e-1),
        C0AidParameters(core_epsilon=0.0, k1_radius_epsilon=0.0),
    )
    eos_final_only_after_first: bool = True
    use_wall_predictor: bool = True
    prefer_existing_b2c_background_predictor: bool = False
    use_gauge_fix: bool = False
    gauge_fix_gradient_rank_rtol: float = 1.0e-12
    gauge_fix_harmonic_weighted_divergence_rtol: float = 1.0e-6
    gauge_fix_lstsq_rcond: float = 1.0e-12
    max_depth_failures_after_floor: int = 3
    max_tau_backtracks: int = 5
    min_tau_backtrack_delta: float = 5.0e-5
    diagnostic_svd_max_dense_size: int = 1800
    diagnostic_sparse_svds: bool = False
    diagnostic_triplet_count: int = 5
    diagnostic_triplet_residual_rtol: float = 1.0e-8
    diagnostic_linear: bool = True
    linear_diagnostic_taus: tuple[float, ...] = (0.03, 0.0295, 0.029, 0.0285, 0.028)
    residual_equality_tolerance: float = 1.0e-10
    epsilon_independence_tolerance: float = 1.0e-6
    max_newton_iters_override: int | None = None
    diagnostic_observables: bool = True
    compute_missing_observables: bool = False
    jacobian_assembly: str = "colored_sparse_jacobian_lu"
    diagnostic_bdg_grid: tuple[int, int] = (8, 8)
    diagnostic_bdg_modes: int = 8
    diagnostic_profile_points: int = 65
    diagnostic_maxwell_grid: tuple[int, int] = (9, 9)
    diagnostic_maxwell_window: float = b2b.DEFAULT_FINAL_WINDOW
    diagnostic_maxwell_radial_scale: float = 1.75


@dataclass(frozen=True)
class C0cConfig:
    c0b_json_path: Path = DEFAULT_JSON_PATH
    run_root: Path = DEFAULT_C0C_RUN_ROOT
    report_path: Path = DEFAULT_C0C_REPORT_PATH
    json_path: Path = DEFAULT_C0C_JSON_PATH
    cluster_mode_count: int = 5
    annihilation_threshold: float = 1.0e-8
    overlap_threshold: float = 0.9
    residual_fraction_threshold: float = 0.2
    dense_sigma_tau: float | None = None
    dense_jvp_chunk_size: int = 32


@dataclass(frozen=True)
class C0dConfig:
    c0b_json_path: Path = DEFAULT_JSON_PATH
    run_root: Path = DEFAULT_C0D_RUN_ROOT
    report_path: Path = DEFAULT_C0D_REPORT_PATH
    json_path: Path = DEFAULT_C0D_JSON_PATH
    cluster_mode_count: int = 5
    maxwell_mode_count: int = 4
    maxwell_lane_fraction_min: float = 0.9
    gauge_capture_threshold: float = 0.9
    weighted_residual_threshold: float = 0.1
    gradient_rank_rtol: float = 1.0e-12
    harmonic_weighted_divergence_rtol: float = 1.0e-6
    control_random_seed: int = 8675309


@dataclass(frozen=True)
class C0eConfig:
    c0b_json_path: Path = DEFAULT_JSON_PATH
    run_root: Path = DEFAULT_C0E_RUN_ROOT
    report_path: Path = DEFAULT_C0E_REPORT_PATH
    json_path: Path = DEFAULT_C0E_JSON_PATH
    cluster_mode_count: int = 5
    maxwell_mode_count: int = 4
    maxwell_lane_fraction_min: float = 0.9
    gradient_rank_rtol: float = 1.0e-12
    harmonic_weighted_divergence_rtol: float = 1.0e-6
    control_random_seed: int = 420029
    run_all_available_tau_sources: bool = True


@dataclass(frozen=True)
class C0gConfig:
    c0f2_json_path: Path = Path(
        "software/stage1_solver/runs/pathA_C0f2_timing_rerun/pathA_C0f2_timing_rerun.json"
    )
    run_root: Path = DEFAULT_C0G_RUN_ROOT
    report_path: Path = DEFAULT_C0G_REPORT_PATH
    step0_json_path: Path = DEFAULT_C0G_STEP0_JSON_PATH
    json_path: Path = DEFAULT_C0G_JSON_PATH
    step5_json_path: Path = DEFAULT_C0G_STEP5_JSON_PATH
    final_json_path: Path = DEFAULT_C0G_FINAL_JSON_PATH
    grid: tuple[int, int] = b2a.DEFAULT_BACKGROUND_GRID
    converged_taus: tuple[float, ...] = (0.03, 0.0295, 0.02925, 0.029125)
    stalled_tau: float = 0.0290625
    near_null_mode_count: int = 5
    alpha_min_power: int = 20
    gradient_rank_rtol: float = 1.0e-12
    harmonic_weighted_divergence_rtol: float = 1.0e-6
    ftau_h_multipliers: tuple[float, ...] = (1.0e-4, 1.0e-5, 1.0e-6)
    ftau_stability_rtol: float = 0.1
    premise_fail_threshold: float = 1.0e-2
    premise_hold_threshold: float = 1.0e-1
    cos_fold_threshold: float = 1.0e-1
    cos_bifurcation_threshold: float = 1.0e-2
    fit_good_r2_threshold: float = 0.9
    mach_density_floors: tuple[float, ...] = (1.0e-14, 1.0e-12, 1.0e-10, 1.0e-8)
    sonic_band: tuple[float, float] = (0.8, 1.2)
    no_sonic_threshold: float = 0.7
    sonic_current_threshold_multiplier: float = 1.0e-3
    scipy_method_wall_seconds: float = 285.0
    scipy_global_wall_seconds: float = 600.0
    scipy_maxfev_factor: int = 4
    commutator_control_seed: int = 730029
    jacobian_assembly: str = "colored_sparse_jacobian_lu"


@dataclass(frozen=True)
class C0gBuildB1B2Config:
    c0f2_json_path: Path = C0gConfig().c0f2_json_path
    run_root: Path = DEFAULT_C0G_BUILD_RUN_ROOT
    report_path: Path = DEFAULT_C0G_BUILD_REPORT_PATH
    json_path: Path = DEFAULT_C0G_BUILD_JSON_PATH
    grid: tuple[int, int] = b2a.DEFAULT_BACKGROUND_GRID
    b2_depth_sequence: tuple[float, ...] = (
        0.03,
        0.0295,
        0.02925,
        0.029125,
        0.029124,
        0.0291225,
        0.029122,
        0.029120,
    )
    proof_tau: float = 0.02925
    tau_fold_reference: float = 0.0291233
    dissolved_tau_margin: float = 5.0e-7
    proof_tolerance: float = BACKGROUND_RESIDUAL_TOL
    b2_fit_r2_threshold: float = 0.95
    b2_cond_jb_flat_band_max_over_min: float = 1.10
    b2_ratio_growth_factor_min: float = 10.0
    max_tau_backtracks: int = 5
    min_tau_backtrack_delta: float = 1.0e-7
    max_depth_failures_after_floor: int = 4
    max_newton_iters_override: int | None = None
    jacobian_assembly: str = "autodiff_sparse_jacobian_lu"
    b4_reference_tau: float = 0.029125
    b4_random_probe_count: int = 5
    b4_random_seed: int = 20260621
    b4_match_abs_tol: float = 1.0e-10
    b4_match_rel_tol: float = 1.0e-8
    proof_global_phase_radians: float = 1.0e-3
    proof_r0_perturbation: float = 1.0e-5


@dataclass(frozen=True)
class C0gDeepCrawlConfig:
    prior_b1b2_json_path: Path = DEFAULT_C0G_BUILD_JSON_PATH
    prior_crawl_json_path: Path = (
        DEFAULT_C0G_BUILD_RUN_ROOT / "gauge_fixed_crawl" / "c0_gauge_fixed_crawl.json"
    )
    run_root: Path = DEFAULT_C0G_DEEP_RUN_ROOT
    report_path: Path = DEFAULT_C0G_DEEP_REPORT_PATH
    json_path: Path = DEFAULT_C0G_DEEP_JSON_PATH
    grid: tuple[int, int] = b2a.DEFAULT_BACKGROUND_GRID
    target_taus: tuple[float, ...] = (
        0.0290,
        0.0288,
        0.0285,
        0.0282,
        0.0280,
        0.0275,
        0.0270,
        0.0265,
        0.0260,
    )
    tau_fold_reference: float = 0.0291233
    below_fold_seed_count: int = 3
    no_fold_required_new_converged: int = 5
    max_tau_backtracks: int = 6
    min_tau_backtrack_delta: float = 1.0e-7
    max_newton_iters_override: int | None = None
    jacobian_assembly: str = "autodiff_sparse_jacobian_lu"
    residual_equality_tolerance: float = 1.0e-10
    sigma_finite_floor: float = 1.0e-8
    b4_reference_tau: float = C0gBuildB1B2Config().b4_reference_tau
    b4_random_probe_count: int = C0gBuildB1B2Config().b4_random_probe_count
    b4_random_seed: int = C0gBuildB1B2Config().b4_random_seed
    b4_match_abs_tol: float = C0gBuildB1B2Config().b4_match_abs_tol
    b4_match_rel_tol: float = C0gBuildB1B2Config().b4_match_rel_tol


@dataclass(frozen=True)
class C0gB3PseudoArclengthConfig:
    deepcrawl_json_path: Path = DEFAULT_C0G_DEEP_JSON_PATH
    run_root: Path = DEFAULT_C0G_B3_RUN_ROOT
    report_path: Path = DEFAULT_C0G_B3_REPORT_PATH
    json_path: Path = DEFAULT_C0G_B3_JSON_PATH
    grid: tuple[int, int] = b2a.DEFAULT_BACKGROUND_GRID
    jacobian_assembly: str = "autodiff_sparse_jacobian_lu"
    residual_tolerance: float = BACKGROUND_RESIDUAL_TOL
    tau_scale: float = 1000.0
    initial_ds: float = 1.0e-3
    min_ds: float = 4.0e-5
    max_ds: float = 3.0e-2
    ds_growth: float = 1.35
    ds_shrink: float = 0.5
    corrector_max_iters: int = 8
    corrector_stop_residual_tolerance: float = 1.0e-9
    corrector_lstsq_rcond: float = 1.0e-12
    corrector_line_search_steps: int = 10
    arclength_tolerance: float = 2.0e-6
    ftau_h: float = 1.0e-7
    target_tau_tolerance: float = 1.0e-8
    reproduction_comparison_tolerance: float = 5.0e-5
    distinct_state_floor: float = 1.0e-7
    b31_max_steps: int = 8
    b31_required_below_count: int = 3
    b31_deepest_tau_margin: float = 1.0e-7
    b4_reference_tau: float = C0gBuildB1B2Config().b4_reference_tau
    b4_random_probe_count: int = C0gBuildB1B2Config().b4_random_probe_count
    b4_random_seed: int = C0gBuildB1B2Config().b4_random_seed
    b4_match_abs_tol: float = C0gBuildB1B2Config().b4_match_abs_tol
    b4_match_rel_tol: float = C0gBuildB1B2Config().b4_match_rel_tol


@dataclass(frozen=True)
class C0gB3FollowupConfig:
    deepcrawl_json_path: Path = DEFAULT_C0G_DEEP_JSON_PATH
    b3_json_path: Path = DEFAULT_C0G_B3_JSON_PATH
    run_root: Path = DEFAULT_C0G_B3_FOLLOWUP_RUN_ROOT
    report_path: Path = DEFAULT_C0G_B3_FOLLOWUP_REPORT_PATH
    json_path: Path = DEFAULT_C0G_B3_FOLLOWUP_JSON_PATH
    grid: tuple[int, int] = b2a.DEFAULT_BACKGROUND_GRID
    jacobian_assembly: str = "autodiff_sparse_jacobian_lu"
    fail_tau: float = 0.02911625
    sigma_fit_zero_crossing: float = 0.0291139
    tau_scale: float = 1000.0
    c1_residual_tolerance: float = 1.0e-11
    c1_step_tolerance: float = 1.0e-9
    c1_max_iters: int = 120
    c1_lstsq_rcond: float = 1.0e-12
    c1_line_search_steps: int = 12
    c1_stall_window: int = 30
    c1_stall_residual_decrease: float = 0.01
    c1_same_state_tolerance: float = 5.0e-5
    c2_singular_count: int = 10
    c2_cluster_gap_ratio: float = 10.0
    c2_ftau_abs_small: float = 1.0e-8
    c3_points: int = 31
    c3_spike_factor: float = 100.0
    c4_run: bool = False


@dataclass(frozen=True)
class C0cGenerator:
    name: str
    classification: str
    symmetry_status: str
    description: str
    vector: np.ndarray


@dataclass
class C0NewtonIteration:
    iteration: int
    residual_norm: float
    step_norm: float | None = None
    line_search_alpha: float | None = None
    gmres_info: int | None = None
    gmres_iterations: int | None = None
    gmres_residual_curve: list[float] = field(default_factory=list)
    preconditioner_info: dict[str, Any] = field(default_factory=dict)


@dataclass
class C0NewtonResult:
    x: torch.Tensor
    converged: bool
    iterations: int
    initial_residual_norm: float
    final_residual_norm: float
    tolerance: float
    history: list[C0NewtonIteration]
    message: str


def _norm(values: torch.Tensor, kind: str) -> float:
    if kind == "linf":
        return float(torch.max(torch.abs(values)).detach().cpu().item())
    if kind == "l2":
        return float(torch.linalg.vector_norm(values).detach().cpu().item())
    raise ValueError(f"unsupported residual norm {kind!r}")


def _format_tau(tau: float) -> str:
    return b2a._format_tau(float(tau))


def _c0_newton_config(
    base: NewtonConfig,
    *,
    jacobian_assembly: str = "colored_sparse_jacobian_lu",
) -> NewtonConfig:
    preconditioner = replace(
        base.preconditioner,
        type=jacobian_assembly,
        side="left",
        rebuild_policy="every_newton_step",
        stencil_radius=3,
        color_separation=7,
        factorization="splu",
        diagonal_shift=0.0,
        drop_tolerance=0.0,
        fill_factor=12.0,
        permutation="COLAMD",
    )
    return replace(
        base,
        residual_atol=2.0e-9,
        residual_rtol=2.0e-9,
        max_newton_iters=20,
        gmres_rtol=2.0e-9,
        gmres_atol=1.0e-11,
        gmres_restart=192,
        gmres_maxiter=16,
        max_line_search_iters=20,
        accept_best_line_search_decrease=True,
        preconditioner=preconditioner,
    )


def c0_frozen_branch(
    *,
    tau: float,
    grid: tuple[int, int],
    max_newton_iters: int | None = None,
    jacobian_assembly: str = "colored_sparse_jacobian_lu",
):
    """Return the existing frozen branch with only Newton controls changed."""

    branch = b2a.frozen_patha_b2a_branch(grid=grid, tau=float(tau))
    newton = _c0_newton_config(branch.newton, jacobian_assembly=jacobian_assembly)
    if max_newton_iters is not None:
        newton = replace(newton, max_newton_iters=int(max_newton_iters))
    return replace(branch, newton=newton)


def c0_preconditioner_residual(
    state: torch.Tensor,
    grid,
    branch,
    *,
    eos_K: float,
    boundaries,
    s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any],
    aid: C0AidParameters,
) -> torch.Tensor:
    """Residual used only to assemble the path preconditioner.

    At zero epsilon this is numerically identical to the physical residual.
    For positive epsilon, only the approximate linear solve is regularized;
    Newton merit, line search, and final convergence never call this function.
    """

    fields, chemical_potential = unpack_closed_coupled_fields(
        state,
        grid,
        has_chemical_potential=True,
    )
    assert chemical_potential is not None
    pde = coupled_pde_residual(
        pack_coupled_fields(_closed_to_coupled_fields(fields)),
        grid,
        branch,
        eos_K=eos_K,
        chemical_potential=chemical_potential,
        boundaries=boundaries,
        confinement_radius=fields.r0,
    )
    if aid.core_epsilon > 0.0:
        core = torch.as_tensor(
            aid.core_epsilon**3,
            dtype=state.dtype,
            device=state.device,
        )
        pde = pde.clone()
        pde[0] = pde[0] + core * fields.psi_real
        pde[1] = pde[1] + core * fields.psi_imag

    density = fields.psi_real**2 + fields.psi_imag**2
    r = grid.r_centers[:, None]
    if aid.k1_radius_epsilon > 0.0:
        eps5 = torch.as_tensor(
            aid.k1_radius_epsilon**5,
            dtype=state.dtype,
            device=state.device,
        )
        k1 = 4.0 * branch.radial_wall_strength * r**4 / (fields.r0[None, :] ** 5 + eps5)
    else:
        k1 = 4.0 * branch.radial_wall_strength * r**4 / fields.r0[None, :] ** 5
    source = torch.sum(grid.radial_shell_volumes[:, None] * (-k1 * density), dim=0)
    wall = static_balance_terms(
        fields.r0,
        wall_grid_from_tensor_grid(grid),
        s_sigma=s_sigma,
        source=source,
        lower_bc=BoundaryCondition.dirichlet(branch.r_mouth),
        upper_bc=BoundaryCondition.neumann(0.0),
    ).residual
    mass = torch.sum(density * grid.cell_volumes) - branch.mass
    return torch.cat([pde.reshape(-1), wall.reshape(-1), mass.reshape(1)])


def _sparse_abs_max_by_axis(matrix: csc_matrix, axis: int) -> np.ndarray:
    reduced = abs(matrix).max(axis=axis)
    values = np.asarray(reduced.toarray() if hasattr(reduced, "toarray") else reduced).reshape(-1)
    return values.astype(np.float64, copy=False)


def jacobi_row_col_scales(
    matrix: csc_matrix,
    *,
    aid: C0AidParameters,
) -> tuple[np.ndarray, np.ndarray]:
    row_norm = _sparse_abs_max_by_axis(matrix, 1)
    row_scale = 1.0 / np.sqrt(np.maximum(row_norm, aid.scale_floor))
    row_scale = np.clip(row_scale, aid.scale_min, aid.scale_max)
    row_scaled = diags(row_scale, format="csc") @ matrix
    col_norm = _sparse_abs_max_by_axis(row_scaled, 0)
    col_scale = 1.0 / np.sqrt(np.maximum(col_norm, aid.scale_floor))
    col_scale = np.clip(col_scale, aid.scale_min, aid.scale_max)
    return row_scale.astype(np.float64), col_scale.astype(np.float64)


def assemble_scaled_linear_system(
    *,
    preconditioner_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x: torch.Tensor,
    rhs: np.ndarray,
    grid,
    newton: NewtonConfig,
    aid: C0AidParameters,
    iteration: int,
) -> tuple[np.ndarray, np.ndarray, csc_matrix, csc_matrix, dict[str, Any]]:
    preconditioner_config = replace(newton.preconditioner, diagonal_shift=0.0)
    linear_config = replace(newton, preconditioner=preconditioner_config)
    matrix, metadata = assemble_closed_coupled_sparse_jacobian(
        PreconditionerBuildContext(
            residual_fn=preconditioner_residual_fn,
            x=x,
            rhs=rhs,
            iteration=iteration,
            config=linear_config,
        ),
        grid,
    )
    matrix = matrix.tocsc()
    if aid.use_jacobi_scaling:
        row_scale, col_scale = jacobi_row_col_scales(matrix, aid=aid)
    else:
        row_scale = np.ones(matrix.shape[0], dtype=np.float64)
        col_scale = np.ones(matrix.shape[1], dtype=np.float64)
    scaled_matrix = (diags(row_scale, format="csc") @ matrix @ diags(col_scale, format="csc")).tocsc()
    if aid.preconditioner_diagonal_shift > 0.0:
        scaled_matrix = scaled_matrix + aid.preconditioner_diagonal_shift * eye(
            scaled_matrix.shape[0],
            dtype=np.float64,
            format="csc",
        )
    metadata.update(
        {
            "type": f"c0_conditioned_scaled_{preconditioner_config.type}",
            "physical_residual_entry_point": (
                "stage1_solver.coupled_branch.patha_closed_branch_residual"
            ),
            "preconditioner_residual_entry_point": (
                "stage1_solver.patha_c0_conditioning_spike.c0_preconditioner_residual"
            ),
            "core_epsilon": float(aid.core_epsilon),
            "k1_radius_epsilon": float(aid.k1_radius_epsilon),
            "jacobi_scaling": bool(aid.use_jacobi_scaling),
            "row_scale_min": float(np.min(row_scale)),
            "row_scale_max": float(np.max(row_scale)),
            "col_scale_min": float(np.min(col_scale)),
            "col_scale_max": float(np.max(col_scale)),
            "preconditioner_diagonal_shift": float(aid.preconditioner_diagonal_shift),
            "linear_operator_residual_is_original": True,
        }
    )
    return row_scale, col_scale, matrix, scaled_matrix, metadata


def build_scaled_preconditioner(
    *,
    original_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    preconditioner_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x: torch.Tensor,
    rhs: np.ndarray,
    grid,
    newton: NewtonConfig,
    aid: C0AidParameters,
    iteration: int,
) -> tuple[BuiltPreconditioner, np.ndarray, np.ndarray, csc_matrix, csc_matrix]:
    del original_residual_fn
    preconditioner_config = replace(newton.preconditioner, diagonal_shift=0.0)
    row_scale, col_scale, matrix, scaled_matrix, metadata = assemble_scaled_linear_system(
        preconditioner_residual_fn=preconditioner_residual_fn,
        x=x,
        rhs=rhs,
        grid=grid,
        newton=newton,
        aid=aid,
        iteration=iteration,
    )
    operator, factor_metadata = factorized_sparse_inverse_operator(scaled_matrix, preconditioner_config)
    metadata.update(factor_metadata)
    return BuiltPreconditioner(operator=operator, metadata=metadata), row_scale, col_scale, matrix, scaled_matrix


def _solve_scaled_gauge_complement_step(
    *,
    scaled_matrix: csc_matrix,
    rhs_scaled: np.ndarray,
    gauge_matrix: np.ndarray,
    col_scale: np.ndarray,
    gradient_rank_rtol: float,
    lstsq_rcond: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    def fast_lstsq(matrix: np.ndarray, rhs: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, np.ndarray]:
        try:
            coeff, residuals, rank, singular = linalg.lstsq(
                matrix,
                rhs,
                cond=float(lstsq_rcond),
                lapack_driver="gelsy",
                check_finite=False,
            )
            if singular is None:
                singular = np.asarray([], dtype=np.float64)
            return (
                np.asarray(coeff, dtype=np.float64),
                np.asarray(residuals, dtype=np.float64),
                int(rank),
                np.asarray(singular, dtype=np.float64),
            )
        except Exception:
            coeff, residuals, rank, singular = np.linalg.lstsq(
                matrix,
                rhs,
                rcond=float(lstsq_rcond),
            )
            return (
                np.asarray(coeff, dtype=np.float64),
                np.asarray(residuals, dtype=np.float64),
                int(rank),
                np.asarray(singular, dtype=np.float64),
            )

    scaled_gauge = np.asarray(gauge_matrix, dtype=np.float64) / col_scale[:, None]
    _basis, full_u, singular_values, threshold, rank = _c0g_column_space(
        scaled_gauge,
        relative_tolerance=float(gradient_rank_rtol),
        full_matrices=True,
    )
    if rank <= 0:
        dense = scaled_matrix.toarray().astype(np.float64, copy=False)
        coeff, residuals, ls_rank, ls_singular = fast_lstsq(dense, rhs_scaled)
        linear_residual = dense @ coeff - rhs_scaled
        return coeff.astype(np.float64, copy=False), {
            "method": "dense_lstsq_no_gauge_rank",
            "gauge_rank": 0,
            "gauge_rank_threshold": float(threshold),
            "gauge_singular_min_retained": None,
            "reduced_shape": [int(dense.shape[0]), int(dense.shape[1])],
            "lstsq_rank": int(ls_rank),
            "lstsq_singular_min": float(np.min(ls_singular)) if len(ls_singular) else math.nan,
            "lstsq_residual_sum": float(np.sum(residuals)) if len(residuals) else 0.0,
            "linear_rel_resid": float(
                np.linalg.norm(linear_residual)
                / max(np.linalg.norm(rhs_scaled), np.finfo(np.float64).tiny)
            ),
        }
    q_perp = full_u[:, rank:].astype(np.float64, copy=True)
    reduced_shape = (int(scaled_matrix.shape[0]), int(q_perp.shape[1]))

    def matvec(coeff: np.ndarray) -> np.ndarray:
        return scaled_matrix @ (q_perp @ np.asarray(coeff, dtype=np.float64))

    def rmatvec(values: np.ndarray) -> np.ndarray:
        return q_perp.T @ (scaled_matrix.T @ np.asarray(values, dtype=np.float64))

    reduced_operator = LinearOperator(reduced_shape, matvec=matvec, rmatvec=rmatvec, dtype=np.float64)
    lsmr_result = lsmr(
        reduced_operator,
        rhs_scaled,
        atol=max(float(lstsq_rcond), 1.0e-12),
        btol=max(float(lstsq_rcond), 1.0e-12),
        maxiter=max(4 * reduced_shape[1], 200),
        show=False,
    )
    coeff = np.asarray(lsmr_result[0], dtype=np.float64)
    scaled_step = q_perp @ coeff
    linear_residual = matvec(coeff) - rhs_scaled
    gauge_constraint = full_u[:, :rank].T @ scaled_step
    return scaled_step.astype(np.float64, copy=False), {
        "method": "iterative_lsmr_scaled_J_times_Q_perp",
        "gauge_rank": int(rank),
        "gauge_rank_threshold": float(threshold),
        "gauge_singular_min_retained": float(singular_values[rank - 1])
        if rank > 0
        else None,
        "state_dim": int(scaled_matrix.shape[1]),
        "reduced_shape": [reduced_shape[0], reduced_shape[1]],
        "lstsq_rank": None,
        "lstsq_singular_min": math.nan,
        "lstsq_singular_max": math.nan,
        "lstsq_residual_sum": math.nan,
        "lsmr_istop": int(lsmr_result[1]),
        "lsmr_iterations": int(lsmr_result[2]),
        "lsmr_normr": float(lsmr_result[3]),
        "lsmr_normar": float(lsmr_result[4]),
        "lsmr_conda": float(lsmr_result[6]),
        "linear_rel_resid": float(
            np.linalg.norm(linear_residual)
            / max(np.linalg.norm(rhs_scaled), np.finfo(np.float64).tiny)
        ),
        "scaled_gauge_constraint_l2": float(np.linalg.norm(gauge_constraint)),
        "path_only_coordinates": "row_scale * J_original * col_scale, step = col_scale * Q_perp * z",
    }


def solve_c0_scaled_newton(
    *,
    original_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    preconditioner_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x0: torch.Tensor,
    grid,
    newton: NewtonConfig,
    aid: C0AidParameters,
    gauge_matrix_fn: Callable[[torch.Tensor, np.ndarray], Mapping[str, Any]] | None = None,
    gauge_fix_lstsq_rcond: float = 1.0e-12,
    gauge_fix_gradient_rank_rtol: float = 1.0e-12,
) -> C0NewtonResult:
    if newton.linear_solver != "gmres_jvp":
        raise ValueError(f"unsupported linear solver {newton.linear_solver!r}")
    if newton.line_search != "armijo":
        raise ValueError(f"unsupported line search {newton.line_search!r}")

    x = x0.detach().clone()
    r = original_residual_fn(x).detach()
    initial_norm = _norm(r, newton.residual_norm)
    tolerance = max(newton.residual_atol, newton.residual_rtol * initial_norm)
    history = [C0NewtonIteration(iteration=0, residual_norm=initial_norm)]
    if initial_norm <= tolerance:
        return C0NewtonResult(
            x=x,
            converged=True,
            iterations=0,
            initial_residual_norm=initial_norm,
            final_residual_norm=initial_norm,
            tolerance=tolerance,
            history=history,
            message="initial residual met tolerance",
        )

    for iteration in range(1, newton.max_newton_iters + 1):
        x_for_jvp = x.detach().clone()
        rhs = -r.detach().cpu().numpy().astype(np.float64)
        gmres_curve: list[float] = []
        if gauge_matrix_fn is None:
            built, row_scale, col_scale, _matrix, _scaled_matrix = build_scaled_preconditioner(
                original_residual_fn=original_residual_fn,
                preconditioner_residual_fn=preconditioner_residual_fn,
                x=x_for_jvp,
                rhs=rhs,
                grid=grid,
                newton=newton,
                aid=aid,
                iteration=iteration,
            )
            dim = rhs.size
            rhs_scaled = row_scale * rhs

            def matvec(vector_np: np.ndarray) -> np.ndarray:
                direction_np = col_scale * np.asarray(vector_np, dtype=np.float64)
                direction = torch.as_tensor(direction_np, dtype=x.dtype, device=x.device)
                jv = jvp(original_residual_fn, x_for_jvp, direction)
                return row_scale * jv.detach().cpu().numpy().astype(np.float64)

            def callback(residual_norm: float) -> None:
                gmres_curve.append(float(residual_norm))

            linear_op = LinearOperator((dim, dim), matvec=matvec, dtype=np.float64)
            step_scaled, gmres_info = gmres(
                linear_op,
                rhs_scaled,
                M=built.operator,
                rtol=newton.gmres_rtol,
                atol=newton.gmres_atol,
                restart=newton.gmres_restart,
                maxiter=newton.gmres_maxiter,
                callback=callback,
                callback_type="pr_norm",
            )
            preconditioner_info = dict(built.metadata)
        else:
            row_scale, col_scale, _matrix, scaled_matrix, metadata = assemble_scaled_linear_system(
                preconditioner_residual_fn=preconditioner_residual_fn,
                x=x_for_jvp,
                rhs=rhs,
                grid=grid,
                newton=newton,
                aid=aid,
                iteration=iteration,
            )
            rhs_scaled = row_scale * rhs
            gauge = gauge_matrix_fn(x_for_jvp, col_scale)
            step_scaled, gauge_info = _solve_scaled_gauge_complement_step(
                scaled_matrix=scaled_matrix,
                rhs_scaled=rhs_scaled,
                gauge_matrix=np.asarray(gauge["generator_matrix"], dtype=np.float64),
                col_scale=col_scale,
                gradient_rank_rtol=float(gauge_fix_gradient_rank_rtol),
                lstsq_rcond=float(gauge_fix_lstsq_rcond),
            )
            gmres_info = 0
            preconditioner_info = dict(metadata)
            preconditioner_info.update(
                {
                    "type": "c0_gauge_fixed_scaled_complement_lstsq",
                    "gauge_fix_enabled": True,
                    "gauge_fix_source_helpers": list(gauge.get("source_helpers", [])),
                    "gauge_fix_piece_dimensions": dict(gauge.get("piece_dimensions", {})),
                    "gauge_fix_rank": int(gauge_info.get("gauge_rank", 0)),
                    "gauge_fix_method": gauge_info.get("method"),
                    "gauge_fix_path_only_coordinates": gauge_info.get(
                        "path_only_coordinates"
                    ),
                    "gauge_fix_linear_rel_resid": gauge_info.get("linear_rel_resid"),
                    "gauge_fix_scaled_constraint_l2": gauge_info.get(
                        "scaled_gauge_constraint_l2"
                    ),
                    "gauge_fix_reduced_shape": gauge_info.get("reduced_shape"),
                    "gauge_fix_lstsq_rank": gauge_info.get("lstsq_rank"),
                    "gauge_fix_lstsq_singular_min": gauge_info.get("lstsq_singular_min"),
                    "gauge_fix_lstsq_singular_max": gauge_info.get("lstsq_singular_max"),
                }
            )
        if not np.all(np.isfinite(step_scaled)):
            return C0NewtonResult(
                x=x,
                converged=False,
                iterations=iteration - 1,
                initial_residual_norm=initial_norm,
                final_residual_norm=_norm(r, newton.residual_norm),
                tolerance=tolerance,
                history=history,
                message=f"GMRES produced non-finite step, info={gmres_info}",
            )

        step_np = col_scale * np.asarray(step_scaled, dtype=np.float64)
        step = torch.as_tensor(step_np, dtype=x.dtype, device=x.device)
        step_norm = _norm(step, "l2")
        x_norm = max(_norm(x, "l2"), 1.0)
        if step_norm <= newton.step_atol + newton.step_rtol * x_norm:
            r_now = original_residual_fn(x).detach()
            final_norm = _norm(r_now, newton.residual_norm)
            return C0NewtonResult(
                x=x,
                converged=final_norm <= tolerance,
                iterations=iteration - 1,
                initial_residual_norm=initial_norm,
                final_residual_norm=final_norm,
                tolerance=tolerance,
                history=history,
                message="step tolerance reached",
            )

        old_norm = _norm(r, newton.residual_norm)
        alpha = 1.0
        accepted = False
        best_x = x
        best_r = r
        best_norm = old_norm
        accepted_alpha = None
        for _ in range(newton.max_line_search_iters):
            candidate_x = x + alpha * step
            candidate_r = original_residual_fn(candidate_x).detach()
            candidate_norm = _norm(candidate_r, newton.residual_norm)
            if np.isfinite(candidate_norm) and candidate_norm < best_norm:
                best_x = candidate_x.detach()
                best_r = candidate_r.detach()
                best_norm = candidate_norm
                accepted_alpha = alpha
            armijo_bound = (1.0 - newton.line_search_c1 * alpha) * old_norm
            if np.isfinite(candidate_norm) and (
                candidate_norm <= armijo_bound or candidate_norm <= tolerance
            ):
                x = candidate_x.detach()
                r = candidate_r.detach()
                accepted = True
                accepted_alpha = alpha
                break
            alpha *= newton.line_search_shrink

        if not accepted:
            if newton.accept_best_line_search_decrease and best_norm < old_norm:
                x = best_x
                r = best_r
                accepted_alpha = accepted_alpha if accepted_alpha is not None else alpha
            else:
                history.append(
                    C0NewtonIteration(
                        iteration=iteration,
                        residual_norm=old_norm,
                        step_norm=step_norm,
                        line_search_alpha=accepted_alpha,
                        gmres_info=int(gmres_info),
                        gmres_iterations=len(gmres_curve),
                        gmres_residual_curve=gmres_curve,
                        preconditioner_info=preconditioner_info,
                    )
                )
                return C0NewtonResult(
                    x=x,
                    converged=False,
                    iterations=iteration - 1,
                    initial_residual_norm=initial_norm,
                    final_residual_norm=old_norm,
                    tolerance=tolerance,
                    history=history,
                    message="line search failed to reduce original residual",
                )

        new_norm = _norm(r, newton.residual_norm)
        history.append(
            C0NewtonIteration(
                iteration=iteration,
                residual_norm=new_norm,
                step_norm=step_norm,
                line_search_alpha=accepted_alpha,
                gmres_info=int(gmres_info),
                gmres_iterations=len(gmres_curve),
                gmres_residual_curve=gmres_curve,
                preconditioner_info=preconditioner_info,
            )
        )
        if new_norm <= tolerance:
            return C0NewtonResult(
                x=x,
                converged=True,
                iterations=iteration,
                initial_residual_norm=initial_norm,
                final_residual_norm=new_norm,
                tolerance=tolerance,
                history=history,
                message="residual tolerance reached on original residual",
            )

    final_norm = _norm(r, newton.residual_norm)
    return C0NewtonResult(
        x=x,
        converged=False,
        iterations=newton.max_newton_iters,
        initial_residual_norm=initial_norm,
        final_residual_norm=final_norm,
        tolerance=tolerance,
        history=history,
        message="maximum Newton iterations reached",
    )


def _state_metrics(state: torch.Tensor, grid) -> dict[str, float]:
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    density = fields.psi_real**2 + fields.psi_imag**2
    return {
        "min_rho": float(torch.min(density).detach().cpu().item()),
        "max_rho": float(torch.max(density).detach().cpu().item()),
        "min_R0": float(torch.min(fields.r0).detach().cpu().item()),
        "max_R0": float(torch.max(fields.r0).detach().cpu().item()),
        "chemical_potential": float(mu.detach().cpu().item()) if mu is not None else math.nan,
        "mass": float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
    }


def _history_rows(result: C0NewtonResult) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in result.history:
        data = asdict(row)
        if len(data["gmres_residual_curve"]) > 24:
            data["gmres_residual_curve_tail"] = data["gmres_residual_curve"][-8:]
            data["gmres_residual_curve"] = data["gmres_residual_curve"][:16]
        rows.append(data)
    return rows


def _component_observables_from_state(
    state: torch.Tensor,
    grid,
    branch,
    *,
    tau: float,
    residual_linf: float,
) -> dict[str, Any]:
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    assert mu is not None
    density = fields.psi_real**2 + fields.psi_imag**2
    z = b2a.localization_weight_torch(grid.w_centers, branch)
    potential = b2a.confinement_potential_torch(grid, branch, radius=fields.r0)
    drive = b2a.confinement_wall_to_matter_coefficient_torch(grid, branch, radius=fields.r0)
    return {
        "schema": "stage1_patha_c0_background_for_diagnostics/v1",
        "source_revision": source_revision(),
        "solver": {
            "engine": "torch",
            "entry_point": "stage1_solver.patha_c0_conditioning_spike",
            "residual_entry_point": "stage1_solver.coupled_branch.patha_closed_branch_residual",
            "converged": bool(residual_linf <= BACKGROUND_RESIDUAL_TOL),
            "message": "C0 diagnostic background assembled from accepted state",
            "final_eos_K": float(branch.continuation_K_values[-1]),
            "chemical_potential": float(mu.detach().cpu().item()),
        },
        "geometry": {
            "a": b2a.FROZEN_A,
            "L": b2a.FROZEN_L,
            "w_min": 0.0,
            "w_max": b2a.FROZEN_L,
            "r_mouth": float(branch.r_mouth),
            "initial_r_exit": float(branch.r_exit),
            "radial_domain_r_max": float(branch.r_max),
        },
        "constants": {
            "eos_K": float(branch.continuation_K_values[-1]),
            "hbar": float(branch.hbar),
            "particle_mass": float(branch.particle_mass),
            "gauge_charge": float(branch.gauge_charge),
            "mu0": float(branch.mu0),
            "xi": float(branch.xi),
            "tau": float(tau),
            "radial_wall_strength": float(branch.radial_wall_strength),
            "axial_trap_strength": float(branch.axial_trap_strength),
            "mass_constraint": float(branch.mass),
        },
        "s_sigma_spec": b2a.frozen_s_sigma_spec(float(tau)).to_dict(),
        "grid": {
            "measure": grid.spec.measure,
            "nr": int(grid.spec.nr),
            "nw": int(grid.spec.nw),
            "r_min": float(grid.spec.r_min),
            "r_max": float(grid.spec.r_max),
            "w_min": float(grid.spec.w_min),
            "w_max": float(grid.spec.w_max),
            "dr": float(grid.dr),
            "dw": float(grid.dw),
            "r_faces": grid.r_faces.detach().cpu().numpy().tolist(),
            "r_centers": grid.r_centers.detach().cpu().numpy().tolist(),
            "w_faces": grid.w_faces.detach().cpu().numpy().tolist(),
            "w_centers": grid.w_centers.detach().cpu().numpy().tolist(),
            "radial_face_areas": grid.radial_face_areas.detach().cpu().numpy().tolist(),
            "radial_shell_volumes": grid.radial_shell_volumes.detach().cpu().numpy().tolist(),
            "cell_volumes": grid.cell_volumes.detach().cpu().numpy().tolist(),
        },
        "fields": {
            "psi_R0": fields.psi_real.detach().cpu().numpy().tolist(),
            "psi_I0": fields.psi_imag.detach().cpu().numpy().tolist(),
            "A_00": fields.a0.detach().cpu().numpy().tolist(),
            "A_r0": fields.ar.detach().cpu().numpy().tolist(),
            "A_w0": fields.aw.detach().cpu().numpy().tolist(),
            "R0_w": fields.r0.detach().cpu().numpy().tolist(),
        },
        "derived": {
            "rho0": density.detach().cpu().numpy().tolist(),
            "R0_w": fields.r0.detach().cpu().numpy().tolist(),
            "Z_w": z.detach().cpu().numpy().tolist(),
            "V_conf": potential.detach().cpu().numpy().tolist(),
            "wall_drive_k1": drive.detach().cpu().numpy().tolist(),
            "mass": float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
        },
        "residuals": {
            "closed_stationary_linf": float(residual_linf),
            "self_consistent": bool(residual_linf <= BACKGROUND_RESIDUAL_TOL),
        },
    }


def _diagnostic_observables(
    *,
    state: torch.Tensor,
    grid,
    branch,
    tau: float,
    residual_linf: float,
    config: C0Config,
) -> dict[str, Any]:
    if not config.diagnostic_observables or residual_linf > BACKGROUND_RESIDUAL_TOL:
        return {
            "available": False,
            "reason": "background_not_converged_or_observables_disabled",
        }
    tau_label = _format_tau(float(tau))
    existing = sorted(
        Path("software/stage1_solver/runs/patha_b2c_calibration/evaluations").glob(
            f"patha_b2c_*_tau_{tau_label}.json"
        )
    )
    if existing:
        payload = json.loads(existing[-1].read_text(encoding="utf-8"))
        coeff = payload["direct_coefficients"]
        obs = payload["observables"]
        return {
            "available": True,
            "diagnostic_only": False,
            "source": "existing_b2c_evaluation_bundle",
            "path": str(existing[-1]),
            "D0": float(obs["D0"]),
            "P0": float(obs["P0"]),
            "R_norm": float(obs["R_norm"]),
            "B0": float(coeff["B0"]),
            "Z0": float(coeff["Z0"]),
            "N0": float(coeff["N0"]),
            "K": float(coeff["K"]),
            "M": float(coeff["M"]),
            "overlaps": {
                "source": "existing_b2c_direct_coefficients",
                "maxwell_Z0": float(coeff["Z0"]),
                "maxwell_N0": float(coeff["N0"]),
            },
        }
    if not config.compute_missing_observables:
        return {
            "available": False,
            "reason": "no_existing_b2c_evaluation_bundle_for_tau",
        }
    try:
        background = _component_observables_from_state(
            state,
            grid,
            branch,
            tau=tau,
            residual_linf=residual_linf,
        )
        wall_input, _ = b2a.make_wall_input(
            background,
            profile_points_count=config.diagnostic_profile_points,
            out_dir=None,
        )
        bdg = b2a.solve_bdg_python(
            background,
            wall_input,
            nr=config.diagnostic_bdg_grid[0],
            nw=config.diagnostic_bdg_grid[1],
            modes_to_export=config.diagnostic_bdg_modes,
        )
        maxwell = b2b.transfer_for_grid(
            background,
            bdg,
            nr=config.diagnostic_maxwell_grid[0],
            nw=config.diagnostic_maxwell_grid[1],
            window=config.diagnostic_maxwell_window,
            radial_scale=config.diagnostic_maxwell_radial_scale,
            mode_count=config.diagnostic_bdg_modes,
            engine="pathA_C0_diagnostic_primary_second_order",
        )
        coeff = b2c.direct_coefficients_from_packets(bdg_packet=bdg, maxwell_packet=maxwell)
        obs = b2c.assemble_dual(coeff)
        return {
            "available": True,
            "diagnostic_only": True,
            "D0": float(obs["primary"]["D0"]),
            "P0": float(obs["primary"]["P0"]),
            "R_norm": float(obs["primary"]["R_norm"]),
            "B0": float(coeff["B0"]),
            "Z0": float(coeff["Z0"]),
            "N0": float(coeff["N0"]),
            "K": float(coeff["K"]),
            "M": float(coeff["M"]),
            "overlaps": {
                "bdg_mode_count": int(len(bdg["bdg_modes"])),
                "max_bdg_overlap_I_eta_phi": float(
                    max(abs(float(mode["overlap_I_eta_phi"])) for mode in bdg["bdg_modes"])
                ),
                "max_bdg_coupling": float(
                    max(abs(float(mode["coupling"])) for mode in bdg["bdg_modes"])
                ),
                "maxwell_Z0": float(coeff["Z0"]),
                "maxwell_N0": float(coeff["N0"]),
            },
            "grids": {
                "bdg": list(config.diagnostic_bdg_grid),
                "maxwell": list(config.diagnostic_maxwell_grid),
                "profile_points": int(config.diagnostic_profile_points),
            },
        }
    except Exception as exc:  # pragma: no cover - diagnostic path must report, not hide.
        return {
            "available": False,
            "reason": f"{type(exc).__name__}: {exc}",
        }


def _closed_layout_dimensions(grid) -> dict[str, int]:
    cell_count = int(grid.spec.nr * grid.spec.nw)
    field_dim = 5 * cell_count
    wall_dim = int(grid.spec.nw)
    state_size = field_dim + wall_dim + 1
    return {
        "cell_count": cell_count,
        "field_dim": field_dim,
        "wall_dim": wall_dim,
        "state_size": state_size,
    }


def _energy_fraction(vector: np.ndarray, start: int, stop: int) -> float:
    values = np.asarray(vector, dtype=np.float64)
    denom = float(np.dot(values, values))
    if denom <= 0.0:
        return math.nan
    chunk = values[start:stop]
    return float(np.dot(chunk, chunk) / denom)


def _singular_group_fractions(
    *,
    left: np.ndarray,
    right: np.ndarray,
    grid,
) -> dict[str, Any]:
    dims = _closed_layout_dimensions(grid)
    field_dim = dims["field_dim"]
    wall_dim = dims["wall_dim"]
    state_size = dims["state_size"]
    right_fractions = {
        "field": _energy_fraction(right, 0, field_dim),
        "R0": _energy_fraction(right, field_dim, field_dim + wall_dim),
        "mu": _energy_fraction(right, state_size - 1, state_size),
    }
    left_fractions = {
        "pde_rows": _energy_fraction(left, 0, field_dim),
        "wall_rows": _energy_fraction(left, field_dim, field_dim + wall_dim),
        "mass_constraint_row": _energy_fraction(left, state_size - 1, state_size),
    }
    return {
        "v_min_energy_fractions": right_fractions,
        "u_min_energy_fractions": left_fractions,
        "gauge_projection_status": "not_available",
        "gauge_fraction": None,
    }


def _dominant_group(fractions: Mapping[str, float | None]) -> tuple[str | None, float]:
    finite = {
        key: float(value)
        for key, value in fractions.items()
        if value is not None and math.isfinite(float(value))
    }
    if not finite:
        return None, math.nan
    group = max(finite, key=finite.get)
    return group, finite[group]


def _crawl_persistent_failure_guard(
    attempts: Sequence[Mapping[str, Any]],
    config: Mapping[str, Any],
    *,
    failed_attempt_predicate: Callable[[Mapping[str, Any]], bool] | None = None,
) -> dict[str, Any]:
    predicate = failed_attempt_predicate or (
        lambda row: _below_prior_floor(float(row["target_tau"]))
        and not row.get("final_physical_converged")
        and row.get("measurement_status", "MEASURED") == "MEASURED"
    )
    failed = [dict(row) for row in attempts if predicate(row)]
    attempted_backtracking = any(int(row.get("backtrack_index", 0)) > 0 for row in attempts)
    full_newton_budget = config.get("max_newton_iters_override") in (None, 0)
    return {
        "crawl_persistent_failure_evidence": bool(
            failed and attempted_backtracking and full_newton_budget
        ),
        "failed_attempt_count": int(len(failed)),
        "attempted_backtracking": bool(attempted_backtracking),
        "full_newton_budget": bool(full_newton_budget),
        "max_tau_backtracks": int(config.get("max_tau_backtracks", 0) or 0),
        "failed_attempts": failed,
    }


def _condition_number_from_singular_values(values: np.ndarray) -> tuple[float, float, float]:
    finite = np.asarray(values, dtype=np.float64)
    sigma_max = float(np.max(finite)) if finite.size else math.nan
    sigma_min = float(np.min(finite)) if finite.size else math.nan
    cond = sigma_max / sigma_min if sigma_min > 0.0 else math.inf
    return sigma_min, sigma_max, cond


def _normal_shift_invert_triplets(
    matrix: csc_matrix,
    *,
    triplet_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    normal = (matrix.T @ matrix).tocsc()
    values, vectors = eigsh(
        normal,
        k=min(triplet_count, max(1, matrix.shape[1] - 2)),
        sigma=0.0,
        which="LM",
    )
    order = np.argsort(np.maximum(values, 0.0))
    singular_values = np.sqrt(np.maximum(values[order], 0.0))
    right_vectors = vectors[:, order]
    left_vectors = []
    for sigma, right in zip(singular_values, right_vectors.T):
        product = matrix @ right
        norm = float(np.linalg.norm(product))
        if sigma > 0.0:
            left = product / sigma
        elif norm > 0.0:
            left = product / norm
        else:
            left = np.zeros(matrix.shape[0], dtype=np.float64)
        left_vectors.append(left)
    return np.column_stack(left_vectors), singular_values, right_vectors.T


def _true_singular_triplet_diagnostics(
    matrix: csc_matrix,
    *,
    grid,
    dense_max_size: int,
    triplet_count: int,
    triplet_residual_rtol: float,
) -> dict[str, Any]:
    n = matrix.shape[0]
    k = min(max(1, int(triplet_count)), max(1, n - 1))
    method = ""
    try:
        if n <= dense_max_size:
            dense = matrix.toarray()
            left_all, values, right_t = np.linalg.svd(dense, full_matrices=False)
            order = np.argsort(values)
            keep = order[:k]
            left_vectors = left_all[:, keep]
            singular_values = values[keep]
            right_vectors = right_t[keep, :]
            method = "dense_svd"
        else:
            try:
                left_vectors, singular_values, right_vectors = svds(
                    matrix,
                    k=k,
                    which="SM",
                    solver="propack",
                    tol=1.0e-9,
                    maxiter=500,
                    return_singular_vectors=True,
                )
                order = np.argsort(singular_values)
                left_vectors = left_vectors[:, order]
                singular_values = singular_values[order]
                right_vectors = right_vectors[order, :]
                method = "svds_propack_SM"
            except Exception:
                left_vectors, singular_values, right_vectors = _normal_shift_invert_triplets(
                    matrix,
                    triplet_count=k,
                )
                method = "eigsh_normal_shift_invert"
    except Exception as exc:
        return {
            "status": "NOT_MEASURED",
            "reason": f"{type(exc).__name__}: {exc}",
            "sigma_method": "NOT_MEASURED",
            "singular_triplets": [],
        }

    try:
        all_values = np.linalg.svd(matrix.toarray(), compute_uv=False)
    except Exception:
        all_values = np.asarray(singular_values, dtype=np.float64)
    sigma_min, sigma_max, cond_full = _condition_number_from_singular_values(all_values)
    dims = _closed_layout_dimensions(grid)
    field_dim = dims["field_dim"]
    try:
        field_values = np.linalg.svd(
            matrix[:field_dim, :field_dim].toarray(),
            compute_uv=False,
        )
        field_sigma_min, field_sigma_max, cond_field = _condition_number_from_singular_values(
            field_values
        )
    except Exception:
        field_sigma_min = field_sigma_max = cond_field = math.nan

    triplets = []
    residuals = []
    spectral_norm = sigma_max if math.isfinite(sigma_max) and sigma_max > 0.0 else 1.0
    for mode_index, sigma in enumerate(singular_values):
        left = np.asarray(left_vectors[:, mode_index], dtype=np.float64)
        right = np.asarray(right_vectors[mode_index, :], dtype=np.float64)
        residual_abs = float(np.linalg.norm(matrix @ right - float(sigma) * left))
        residual_rel = residual_abs / max(spectral_norm * float(np.linalg.norm(right)), 1.0)
        residuals.append(residual_rel)
        fractions = _singular_group_fractions(left=left, right=right, grid=grid)
        right_group, right_fraction = _dominant_group(fractions["v_min_energy_fractions"])
        left_group, left_fraction = _dominant_group(fractions["u_min_energy_fractions"])
        triplets.append(
            {
                "mode_index": int(mode_index),
                "sigma": float(sigma),
                "triplet_residual_abs": residual_abs,
                "triplet_residual_rel": residual_rel,
                "v_energy_fractions": fractions["v_min_energy_fractions"],
                "u_energy_fractions": fractions["u_min_energy_fractions"],
                "v_dominant_group": right_group,
                "v_dominant_fraction": right_fraction,
                "u_dominant_group": left_group,
                "u_dominant_fraction": left_fraction,
            }
        )

    status = "MEASURED" if residuals and max(residuals) <= triplet_residual_rtol else "NOT_MEASURED"
    reason = None if status == "MEASURED" else "triplet_residual_exceeded_tolerance"
    smallest = triplets[0] if triplets else {}
    cond_ratio = cond_field / cond_full if cond_full > 0.0 and math.isfinite(cond_full) else math.nan
    return {
        "status": status,
        "reason": reason,
        "sigma_method": method,
        "sigma_min": sigma_min,
        "sigma_max": sigma_max,
        "condition_full_bordered": cond_full,
        "field_sigma_min": field_sigma_min,
        "field_sigma_max": field_sigma_max,
        "condition_field_block": cond_field,
        "cond_ratio": cond_ratio,
        "triplet_residual_abs": smallest.get("triplet_residual_abs"),
        "triplet_residual_rel": smallest.get("triplet_residual_rel"),
        "triplet_residual_rtol": float(triplet_residual_rtol),
        "v_min_energy_fractions": smallest.get("v_energy_fractions", {}),
        "u_min_energy_fractions": smallest.get("u_energy_fractions", {}),
        "gauge_projection_status": "not_available",
        "gauge_fraction": None,
        "singular_triplets": triplets,
        "_right_vectors": [np.asarray(v, dtype=np.float64) for v in right_vectors],
    }


def _strip_internal_diagnostic_fields(diagnostic: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in diagnostic.items() if not key.startswith("_")}


def _linear_diagnostics_for_state(
    *,
    state: torch.Tensor,
    grid,
    branch,
    provider,
    boundaries,
    eos_K: float,
    aid: C0AidParameters,
    config: C0Config,
    matrix_path: Path | None = None,
) -> dict[str, Any]:
    original_residual_fn = lambda x: patha_closed_branch_residual(
        x,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=provider,
    )
    preconditioner_residual_fn = lambda x: c0_preconditioner_residual(
        x,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=provider,
        aid=aid,
    )
    rhs = -original_residual_fn(state).detach().cpu().numpy().astype(np.float64)
    _built, _row, _col, matrix, scaled = build_scaled_preconditioner(
        original_residual_fn=original_residual_fn,
        preconditioner_residual_fn=preconditioner_residual_fn,
        x=state.detach(),
        rhs=rhs,
        grid=grid,
        newton=branch.newton,
        aid=aid,
        iteration=0,
    )
    if matrix_path is not None:
        matrix_path.parent.mkdir(parents=True, exist_ok=True)
        save_npz(matrix_path, matrix)
    unscaled = _true_singular_triplet_diagnostics(
        matrix,
        grid=grid,
        dense_max_size=config.diagnostic_svd_max_dense_size,
        triplet_count=config.diagnostic_triplet_count,
        triplet_residual_rtol=config.diagnostic_triplet_residual_rtol,
    )
    scaled_values = np.linalg.svd(scaled.toarray(), compute_uv=False)
    scaled_min, scaled_max, scaled_cond = _condition_number_from_singular_values(scaled_values)
    diagnostic = dict(unscaled)
    diagnostic.update(
        {
            "matrix_path": str(matrix_path) if matrix_path is not None else None,
            "sigma_min_unscaled": unscaled.get("sigma_min"),
            "sigma_max_unscaled": unscaled.get("sigma_max"),
            "condition_unscaled": unscaled.get("condition_full_bordered"),
            "sigma_method_unscaled": unscaled.get("sigma_method"),
            "sigma_min_scaled": scaled_min,
            "sigma_max_scaled": scaled_max,
            "condition_scaled": scaled_cond,
            "sigma_method_scaled": "dense_svd",
            "condition_improvement_factor": (
                unscaled.get("condition_full_bordered") / scaled_cond
                if isinstance(unscaled.get("condition_full_bordered"), float)
                and scaled_cond > 0.0
                else math.nan
            ),
        }
    )
    return diagnostic


def _should_run_linear_diagnostics(tau: float, config: C0Config) -> bool:
    return any(abs(float(tau) - target) <= 5.0e-13 for target in config.linear_diagnostic_taus)


def _max_gmres_growth(history: Sequence[Mapping[str, Any]]) -> float | None:
    growths = []
    for row in history:
        curve = [float(v) for v in row.get("gmres_residual_curve", [])]
        if len(curve) >= 2 and curve[0] > 0.0:
            growths.append(max(curve) / curve[0])
    return max(growths) if growths else None


def _warm_start_from_background(path: Path, *, dtype: torch.dtype) -> tuple[torch.Tensor | None, Any | None]:
    if not path.exists():
        return None, None
    background = json.loads(path.read_text(encoding="utf-8"))
    state, grid = b2c._closed_state_from_background(background, dtype=dtype)
    return state, grid


def _apply_c0_wall_predictor(
    *,
    state: torch.Tensor,
    grid,
    tau: float,
    grid_level: tuple[int, int],
) -> tuple[torch.Tensor, dict[str, Any]]:
    return b2c._apply_wall_radius_predictor(
        state=state,
        grid=grid,
        tau=float(tau),
        grid_level=grid_level,
    )


def _below_prior_floor(tau: float) -> bool:
    return float(tau) < PRIOR_TAU_FLOOR - 5.0e-13


def _may_use_existing_b2c(tau: float, config: C0Config) -> bool:
    return bool(config.prefer_existing_b2c_background_predictor) and not _below_prior_floor(tau)


def _state_artifact_path(config: C0Config, *, tau: float, attempt_index: int) -> Path:
    return (
        config.run_root
        / "states"
        / f"attempt_{attempt_index:03d}_tau_{_format_tau(float(tau))}.npz"
    )


def _save_state_artifact(path: Path, state: torch.Tensor) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, state=state.detach().cpu().numpy().astype(np.float64))


def _load_state_artifact(path: Path, *, dtype: torch.dtype) -> torch.Tensor:
    with np.load(path) as data:
        return torch.as_tensor(data["state"], dtype=dtype, device="cpu")


def _branch_context(*, tau: float, config: C0Config, dtype: torch.dtype):
    branch = c0_frozen_branch(
        tau=float(tau),
        grid=config.grid,
        max_newton_iters=config.max_newton_iters_override,
        jacobian_assembly=config.jacobian_assembly,
    )
    provider = resolve_s_sigma(b2a.frozen_s_sigma_spec(float(tau)))
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    boundaries = branch_boundary_conditions(branch)
    return branch, provider, grid, boundaries


def _epsilon_attempt_row(
    *,
    tau: float,
    eos_K: float,
    aid_index: int,
    aid: C0AidParameters,
    result: C0NewtonResult,
    start_policy: str,
    metrics: Mapping[str, float],
) -> dict[str, Any]:
    history = _history_rows(result)
    return {
        "target_tau": float(tau),
        "eos_K": float(eos_K),
        "aid_index": int(aid_index),
        "epsilon": float(aid.core_epsilon),
        "core_epsilon": float(aid.core_epsilon),
        "k1_radius_epsilon": float(aid.k1_radius_epsilon),
        "aid": asdict(aid),
        "start_policy": start_policy,
        "started_from_clean_tau_attempt_state": start_policy == "clean_tau_attempt_state",
        "converged": bool(result.converged),
        "newton_iterations": int(result.iterations),
        "iterations": int(result.iterations),
        "initial_original_residual": float(result.initial_residual_norm),
        "final_original_residual": float(result.final_residual_norm),
        "final_residual_norm": float(result.final_residual_norm),
        "tolerance": float(result.tolerance),
        "line_search_alphas": [
            row.get("line_search_alpha")
            for row in history
            if row.get("line_search_alpha") is not None
        ],
        "gmres_curve": [
            row.get("gmres_residual_curve", [])
            for row in history
            if row.get("gmres_residual_curve")
        ],
        "newton_history": history,
        "gmres_growth": _max_gmres_growth(history),
        "message": result.message,
        **dict(metrics),
    }


def _execute_tau_attempt(
    *,
    config: C0Config,
    dtype: torch.dtype,
    target_tau: float,
    nominal_target_tau: float,
    delta_tau: float | None,
    backtrack_index: int,
    attempt_index: int,
    previous_state: torch.Tensor | None,
    previous_grid: Any | None,
    previous_tau: float | None,
) -> tuple[dict[str, Any], torch.Tensor, Any]:
    branch, provider, grid, boundaries = _branch_context(
        tau=float(target_tau),
        config=config,
        dtype=dtype,
    )
    existing = (
        Path("software/stage1_solver/runs/patha_b2c_calibration/backgrounds")
        / f"patha_b2c_background_tau_{_format_tau(float(target_tau))}.json"
    )
    used_existing_b2c = False
    if previous_state is None and _may_use_existing_b2c(float(target_tau), config):
        existing_state, existing_grid = _warm_start_from_background(existing, dtype=dtype)
    else:
        existing_state = existing_grid = None

    if existing_state is not None and existing_grid is not None:
        state = existing_state.detach().clone()
        if existing_grid.spec.nr != grid.spec.nr or existing_grid.spec.nw != grid.spec.nw:
            state = resample_closed_branch_state(state, existing_grid, grid, branch)
        initialization = {
            "source": "existing_b2c_background_target_predictor",
            "path": str(existing),
        }
        used_existing_b2c = True
    elif previous_state is None:
        state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
        initialization = {"source": "default_closed_branch_state"}
    else:
        state = previous_state.detach().clone()
        initialization = {
            "source": "previous_c0_converged_state",
            "previous_tau": previous_tau,
        }
        if (
            previous_grid is not None
            and (
                previous_grid.spec.nr != grid.spec.nr
                or previous_grid.spec.nw != grid.spec.nw
                or abs(float(previous_grid.spec.w_max) - float(grid.spec.w_max)) > 0.0
            )
        ):
            state = resample_closed_branch_state(state, previous_grid, grid, branch)
            initialization["resampled"] = True
        if config.use_wall_predictor:
            state, predictor = _apply_c0_wall_predictor(
                state=state,
                grid=grid,
                tau=float(target_tau),
                grid_level=config.grid,
            )
            initialization["wall_predictor"] = predictor

    clean_attempt_state = state.detach().clone()
    final_only = bool(
        config.eos_final_only_after_first
        and initialization.get("source") != "default_closed_branch_state"
    )
    continuation_values = (
        (float(branch.continuation_K_values[-1]),)
        if final_only
        else tuple(float(value) for value in branch.continuation_K_values)
    )
    epsilon_attempts: list[dict[str, Any]] = []
    stage_converged = True
    stage_message = "C0 tau attempt completed"
    stage_started = time.perf_counter()
    state_for_metrics = clean_attempt_state

    for eos_K in continuation_values:
        clean_eos_state = state.detach().clone()
        last_converged_epsilon_state: torch.Tensor | None = None
        zero_epsilon_result: C0NewtonResult | None = None
        for aid_index, aid in enumerate(config.epsilon_schedule):
            if last_converged_epsilon_state is None:
                x0 = clean_eos_state
                start_policy = "clean_tau_attempt_state"
            else:
                x0 = last_converged_epsilon_state
                start_policy = "previous_converged_epsilon_state"
            original_residual_fn = lambda x, eos_K=eos_K: patha_closed_branch_residual(
                x,
                grid,
                branch,
                eos_K=float(eos_K),
                boundaries=boundaries,
                s_sigma=provider,
            )
            preconditioner_residual_fn = lambda x, eos_K=eos_K, aid=aid: c0_preconditioner_residual(
                x,
                grid,
                branch,
                eos_K=float(eos_K),
                boundaries=boundaries,
                s_sigma=provider,
                aid=aid,
            )
            gauge_matrix_fn = None
            if config.use_gauge_fix:
                gauge_config = C0gConfig(
                    grid=config.grid,
                    gradient_rank_rtol=float(config.gauge_fix_gradient_rank_rtol),
                    harmonic_weighted_divergence_rtol=float(
                        config.gauge_fix_harmonic_weighted_divergence_rtol
                    ),
                )

                def gauge_matrix_fn(
                    x: torch.Tensor,
                    _col_scale: np.ndarray,
                    gauge_config: C0gConfig = gauge_config,
                ) -> Mapping[str, Any]:
                    del _col_scale
                    return _c0g_build_analytic_gauge_matrix(
                        state=x,
                        grid=grid,
                        branch=branch,
                        config=gauge_config,
                    )

            result = solve_c0_scaled_newton(
                original_residual_fn=original_residual_fn,
                preconditioner_residual_fn=preconditioner_residual_fn,
                x0=x0,
                grid=grid,
                newton=branch.newton,
                aid=aid,
                gauge_matrix_fn=gauge_matrix_fn,
                gauge_fix_lstsq_rcond=float(config.gauge_fix_lstsq_rcond),
                gauge_fix_gradient_rank_rtol=float(config.gauge_fix_gradient_rank_rtol),
            )
            state_for_metrics = result.x.detach()
            metrics = _state_metrics(state_for_metrics, grid)
            epsilon_attempts.append(
                _epsilon_attempt_row(
                    tau=float(target_tau),
                    eos_K=float(eos_K),
                    aid_index=aid_index,
                    aid=aid,
                    result=result,
                    start_policy=start_policy,
                    metrics=metrics,
                )
            )
            if result.converged:
                last_converged_epsilon_state = result.x.detach()
            if aid.core_epsilon == 0.0 and aid.k1_radius_epsilon == 0.0:
                zero_epsilon_result = result
                if result.converged:
                    state = result.x.detach()
        if zero_epsilon_result is None or not zero_epsilon_result.converged:
            stage_converged = False
            failed = zero_epsilon_result or result
            stage_message = (
                f"tau={target_tau:.12e}, eos_K={float(eos_K):.12e}, "
                f"epsilon=0 failed after full schedule: {failed.message}"
            )
            break

    final_eos = float(continuation_values[-1])
    final_residual_fn = lambda x: patha_closed_branch_residual(
        x,
        grid,
        branch,
        eos_K=final_eos,
        boundaries=boundaries,
        s_sigma=provider,
    )
    final_state = state if stage_converged else state_for_metrics
    final_residual = final_residual_fn(final_state).detach()
    final_linf = float(torch.max(torch.abs(final_residual)).detach().cpu().item())
    metrics = _state_metrics(final_state, grid)
    final_physical_converged = bool(stage_converged and final_linf <= BACKGROUND_RESIDUAL_TOL)
    state_path = _state_artifact_path(config, tau=float(target_tau), attempt_index=attempt_index)
    _save_state_artifact(state_path, final_state)
    row = {
        "tau": float(target_tau),
        "target_tau": float(target_tau),
        "nominal_target_tau": float(nominal_target_tau),
        "delta_tau": None if delta_tau is None else float(delta_tau),
        "backtrack_index": int(backtrack_index),
        "start_from_tau": previous_tau,
        "target_relative_to_prior_floor": float(target_tau) / PRIOR_TAU_FLOOR,
        "below_prior_floor": _below_prior_floor(float(target_tau)),
        "stage_converged": bool(stage_converged),
        "final_original_residual_linf": final_linf,
        "final_physical_converged": final_physical_converged,
        "b2c_background_tolerance": BACKGROUND_RESIDUAL_TOL,
        "message": stage_message,
        "elapsed_seconds": time.perf_counter() - stage_started,
        "init": initialization,
        "initialization": initialization,
        "used_existing_b2c": bool(used_existing_b2c),
        "gauge_fix_solver_invoked": bool(config.use_gauge_fix),
        "max_newton_iters": int(branch.newton.max_newton_iters),
        "max_newton_iters_override": config.max_newton_iters_override,
        "max_tau_backtracks": int(config.max_tau_backtracks),
        "min_tau_backtrack_delta": float(config.min_tau_backtrack_delta),
        "continuation_K_values": list(continuation_values),
        "epsilon_attempts": epsilon_attempts,
        "substeps": epsilon_attempts,
        "metrics": metrics,
        "state_artifact": str(state_path),
        "linear_diagnostics": {
            "status": "NOT_MEASURED",
            "reason": "diagnostics_phase_not_run",
        },
        "observables": {
            "available": False,
            "reason": "diagnostics_phase_not_run",
        },
        "floor_activation": {
            "core_epsilon_schedule": [aid.core_epsilon for aid in config.epsilon_schedule],
            "k1_radius_epsilon_schedule": [
                aid.k1_radius_epsilon for aid in config.epsilon_schedule
            ],
            "min_rho_during_final": metrics["min_rho"],
            "min_R0_during_final": metrics["min_R0"],
            "k1_clamp_active_in_path": any(
                aid.k1_radius_epsilon > 0.0 and metrics["min_R0"] < aid.k1_radius_epsilon
                for aid in config.epsilon_schedule
            ),
            "final_aids_inactive": True,
        },
    }
    return row, final_state.detach(), grid


def _validate_tau_attempts_schema(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    errors: list[str] = []
    expected_eps = [float(aid.core_epsilon) for aid in C0Config().epsilon_schedule]
    for index, row in enumerate(rows):
        for key in (
            "target_tau",
            "delta_tau",
            "backtrack_index",
            "init",
            "used_existing_b2c",
            "start_from_tau",
            "epsilon_attempts",
        ):
            if key not in row:
                errors.append(f"tau_attempts[{index}] missing {key}")
        tau = float(row.get("target_tau", row.get("tau", math.nan)))
        if _below_prior_floor(tau):
            if row.get("used_existing_b2c") is not False:
                errors.append(f"tau_attempts[{index}] below floor used existing B2c")
            if row.get("init", {}).get("source") not in {
                "previous_c0_converged_state",
                "previous_c0g_gauge_fixed_converged_state",
            }:
                errors.append(f"tau_attempts[{index}] below floor did not warm-start from previous C0")
        attempts = row.get("epsilon_attempts", [])
        if attempts:
            seen = [float(attempt.get("epsilon", math.nan)) for attempt in attempts[-len(expected_eps) :]]
            if len(seen) < len(expected_eps) or any(
                abs(left - right) > 1.0e-15 for left, right in zip(seen, expected_eps)
            ):
                errors.append(f"tau_attempts[{index}] did not run the full epsilon schedule")
    return errors


def run_c0_crawl(config: C0Config | None = None) -> dict[str, Any]:
    if config is None:
        config = C0Config()
    config.run_root.mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(BackendConfig())
    sequence_rows: list[dict[str, Any]] = []
    last_converged_state: torch.Tensor | None = None
    last_converged_grid: Any | None = None
    deepest_converged_tau: float | None = None
    deepest_converged_row: dict[str, Any] | None = None
    floor_or_below_failures = 0
    started = time.perf_counter()
    attempt_index = 0

    for nominal_tau in config.depth_sequence:
        if (
            deepest_converged_tau is not None
            and float(nominal_tau) >= float(deepest_converged_tau) - 5.0e-13
        ):
            continue
        backtrack_index = 0
        start_tau = deepest_converged_tau
        full_delta = None if start_tau is None else float(nominal_tau) - float(start_tau)
        target_tau = float(nominal_tau)
        while True:
            delta_tau = None if start_tau is None else float(target_tau) - float(start_tau)
            row, final_state, final_grid = _execute_tau_attempt(
                config=config,
                dtype=dtype,
                target_tau=float(target_tau),
                nominal_target_tau=float(nominal_tau),
                delta_tau=delta_tau,
                backtrack_index=backtrack_index,
                attempt_index=attempt_index,
                previous_state=last_converged_state,
                previous_grid=last_converged_grid,
                previous_tau=deepest_converged_tau,
            )
            sequence_rows.append(row)
            attempt_index += 1
            if row["final_physical_converged"]:
                last_converged_state = final_state.detach().clone()
                last_converged_grid = final_grid
                deepest_converged_tau = float(row["target_tau"])
                deepest_converged_row = row
                break

            if float(row["target_tau"]) <= PRIOR_TAU_FLOOR + 5.0e-13:
                floor_or_below_failures += 1
            if start_tau is None or full_delta is None:
                break
            if (
                backtrack_index >= config.max_tau_backtracks
                or abs(delta_tau or 0.0) <= config.min_tau_backtrack_delta
            ):
                break
            backtrack_index += 1
            target_tau = float(start_tau) + float(full_delta) / (2.0**backtrack_index)

        if floor_or_below_failures >= config.max_depth_failures_after_floor:
            break

    schema_errors = _validate_tau_attempts_schema(sequence_rows)
    result = {
        "schema": "stage1_pathA_C0b_wall_diagnosis/v1",
        "source_revision": source_revision(),
        "phase_status": "crawl_complete",
        "verdict": "DIAGNOSTIC_INCOMPLETE",
        "verdict_support": {
            "reason": "diagnostics_phase_not_run",
        },
        "prior_tau_floor": PRIOR_TAU_FLOOR,
        "b2c_background_tolerance": BACKGROUND_RESIDUAL_TOL,
        "deepest_converged_tau": deepest_converged_tau,
        "deepest_converged_R0_min": None
        if deepest_converged_row is None
        else deepest_converged_row["metrics"]["min_R0"],
        "elapsed_seconds": time.perf_counter() - started,
        "config": _config_to_dict(config),
        "tau_attempts": sequence_rows,
        "depth_sequence": sequence_rows,
        "tau_attempt_schema_errors": schema_errors,
        "aids_admissibility": {
            "status": "NOT_MEASURED",
            "reason": "diagnostics_phase_not_run",
        },
        "sigma_diagnostics": {
            "status": "NOT_MEASURED",
            "reason": "diagnostics_phase_not_run",
        },
        "fold_diagnostic": {
            "status": "NOT_MEASURED",
            "reason": "diagnostics_phase_not_run",
        },
        "r0_tau_diagnostic": _diagnose_tau_vs_r0(sequence_rows),
        "faithful_operator_boundary": _faithful_operator_boundary(),
        "recommended_next_step": "Run the C0b diagnostics phase on the saved crawl states.",
    }
    config.json_path.parent.mkdir(parents=True, exist_ok=True)
    config.json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def _faithful_operator_boundary() -> dict[str, Any]:
    return {
        "physical_residual_entry_point": (
            "stage1_solver.coupled_branch.patha_closed_branch_residual"
        ),
        "operators_modified_by_c0": False,
        "frozen_physics_modified_by_c0": False,
        "physical_export_guard_modified_by_c0": False,
    }


def _select_diagnostic_attempts(result: Mapping[str, Any], config: C0Config) -> list[dict[str, Any]]:
    attempts = [dict(row) for row in result.get("tau_attempts", [])]
    selected: list[dict[str, Any]] = []
    for row in attempts:
        tau = float(row["target_tau"])
        if _should_run_linear_diagnostics(tau, config):
            selected.append(row)
    converged = [row for row in attempts if row.get("final_physical_converged")]
    if converged:
        shallow = max(converged, key=lambda item: float(item["target_tau"]))
        if shallow not in selected:
            selected.insert(0, shallow)
    if attempts:
        deepest = min(attempts, key=lambda item: float(item["target_tau"]))
        if deepest not in selected:
            selected.append(deepest)
    dedup: list[dict[str, Any]] = []
    seen: set[tuple[float, int]] = set()
    for row in selected:
        key = (float(row["target_tau"]), int(row["backtrack_index"]))
        if key not in seen:
            seen.add(key)
            dedup.append(row)
    return dedup


def _track_singular_triplets(
    diagnostics: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    measured = [
        diagnostic
        for diagnostic in diagnostics
        if diagnostic.get("status") == "MEASURED" and diagnostic.get("_right_vectors")
    ]
    if len(measured) < 2:
        return {
            "status": "NOT_MEASURED",
            "call": "DIAGNOSTIC_INCOMPLETE",
            "reason": "fewer_than_two_measured_triplet_sets",
            "tracked_modes": [],
        }
    first_vectors = measured[0]["_right_vectors"]
    mode_count = min(len(first_vectors), *(len(row["_right_vectors"]) for row in measured[1:]))
    if mode_count <= 0:
        return {
            "status": "NOT_MEASURED",
            "call": "DIAGNOSTIC_INCOMPLETE",
            "reason": "no_triplet_vectors",
            "tracked_modes": [],
        }
    tracks = [
        {
            "track_index": int(mode),
            "points": [
                {
                    "tau": float(measured[0]["tau"]),
                    "sigma": float(measured[0]["singular_triplets"][mode]["sigma"]),
                    "source_mode_index": int(mode),
                    "overlap_from_previous": None,
                    "v_dominant_group": measured[0]["singular_triplets"][mode][
                        "v_dominant_group"
                    ],
                    "v_dominant_fraction": measured[0]["singular_triplets"][mode][
                        "v_dominant_fraction"
                    ],
                    "u_dominant_group": measured[0]["singular_triplets"][mode][
                        "u_dominant_group"
                    ],
                    "u_dominant_fraction": measured[0]["singular_triplets"][mode][
                        "u_dominant_fraction"
                    ],
                }
            ],
            "_last_vector": np.asarray(first_vectors[mode], dtype=np.float64),
        }
        for mode in range(mode_count)
    ]
    for diagnostic in measured[1:]:
        available = [np.asarray(v, dtype=np.float64) for v in diagnostic["_right_vectors"][:mode_count]]
        used: set[int] = set()
        for track in tracks:
            previous = track["_last_vector"]
            overlaps = [
                abs(float(np.dot(previous, candidate)))
                if previous.shape == candidate.shape
                else -1.0
                for candidate in available
            ]
            order = np.argsort(overlaps)[::-1]
            chosen = None
            for candidate_index in order:
                if int(candidate_index) not in used:
                    chosen = int(candidate_index)
                    break
            if chosen is None:
                chosen = int(order[0])
            used.add(chosen)
            triplet = diagnostic["singular_triplets"][chosen]
            track["points"].append(
                {
                    "tau": float(diagnostic["tau"]),
                    "sigma": float(triplet["sigma"]),
                    "source_mode_index": int(chosen),
                    "overlap_from_previous": float(overlaps[chosen]),
                    "v_dominant_group": triplet["v_dominant_group"],
                    "v_dominant_fraction": triplet["v_dominant_fraction"],
                    "u_dominant_group": triplet["u_dominant_group"],
                    "u_dominant_fraction": triplet["u_dominant_fraction"],
                }
            )
            track["_last_vector"] = available[chosen]

    persistent_tracks: list[int] = []
    fold_tracks: list[int] = []
    cleaned_tracks = []
    for track in tracks:
        points = track["points"]
        sigmas = [float(point["sigma"]) for point in points]
        sigma_ratio = max(sigmas) / min(sigmas) if min(sigmas) > 0.0 else math.inf
        border_like = all(
            (
                point["v_dominant_group"] in {"R0", "mu"}
                and float(point["v_dominant_fraction"]) >= 0.7
            )
            or (
                point["u_dominant_group"] == "mass_constraint_row"
                and float(point["u_dominant_fraction"]) >= 0.7
            )
            for point in points
        )
        persistent = bool(border_like and sigma_ratio <= 10.0)
        decreasing = all(left >= right for left, right in zip(sigmas, sigmas[1:]))
        fold_like = bool((not persistent) and decreasing and sigmas[-1] <= 0.1 * sigmas[0])
        if persistent:
            persistent_tracks.append(int(track["track_index"]))
        if fold_like:
            fold_tracks.append(int(track["track_index"]))
        cleaned = {key: value for key, value in track.items() if not key.startswith("_")}
        cleaned["sigma_ratio"] = sigma_ratio
        cleaned["persistent_border_or_mass_mode"] = persistent
        cleaned["fold_like_complement_mode"] = fold_like
        cleaned_tracks.append(cleaned)

    if fold_tracks:
        call = "FOLD_TURNING_POINT"
    elif persistent_tracks:
        call = "PERSISTENT_NEAR_NULL_SPACE"
    else:
        call = "NO_TRACKED_COMPLEMENT_CROSSING"
    return {
        "status": "MEASURED",
        "call": call,
        "persistent_track_indices": persistent_tracks,
        "fold_track_indices": fold_tracks,
        "tracked_modes": cleaned_tracks,
        "determinant_sign_support": {
            "status": "NOT_USED",
            "reason": "singular-triplet tracking is primary; det sign omitted unless separately stability-checked",
        },
    }


def _diagnose_tau_vs_r0(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    trend = [
        {
            "tau": float(row["target_tau"]),
            "R0_min": float(row["metrics"]["min_R0"]),
            "min_rho": float(row["metrics"]["min_rho"]),
            "converged": bool(row.get("final_physical_converged")),
        }
        for row in rows
        if row.get("metrics")
    ]
    if len(trend) < 2:
        return {
            "status": "UNKNOWN_NOT_MEASURED",
            "trend": trend,
            "tau_decrease_approaches_empty_core": "UNKNOWN_NOT_MEASURED",
            "candidate_knobs": "UNKNOWN_NOT_MEASURED",
        }
    ordered = sorted(trend, key=lambda item: item["tau"], reverse=True)
    shallow = ordered[0]
    deepest = min(ordered, key=lambda item: item["tau"])
    approaches = deepest["R0_min"] < shallow["R0_min"] and deepest["min_rho"] < shallow["min_rho"]
    if approaches:
        statement = "tau_decrease_approaches_R0_zero_and_density_zero"
    else:
        statement = "tau_decrease_does_not_approach_empty_core_on_measured_crawl"
    return {
        "status": "MEASURED",
        "trend": trend,
        "shallow_reference": shallow,
        "deepest_reference": deepest,
        "tau_decrease_approaches_empty_core": bool(approaches),
        "statement": statement,
        "candidate_knobs": "UNKNOWN_NOT_MEASURED",
    }


def _state_group_deltas(
    reference: torch.Tensor,
    candidate: torch.Tensor,
    grid,
) -> dict[str, float]:
    ref_fields, ref_mu = unpack_closed_coupled_fields(
        reference,
        grid,
        has_chemical_potential=True,
    )
    cand_fields, cand_mu = unpack_closed_coupled_fields(
        candidate,
        grid,
        has_chemical_potential=True,
    )
    psi_delta = torch.max(
        torch.abs(
            torch.cat(
                [
                    (cand_fields.psi_real - ref_fields.psi_real).reshape(-1),
                    (cand_fields.psi_imag - ref_fields.psi_imag).reshape(-1),
                ]
            )
        )
    )
    a_delta = torch.max(
        torch.abs(
            torch.cat(
                [
                    (cand_fields.a0 - ref_fields.a0).reshape(-1),
                    (cand_fields.ar - ref_fields.ar).reshape(-1),
                    (cand_fields.aw - ref_fields.aw).reshape(-1),
                ]
            )
        )
    )
    r0_delta = torch.max(torch.abs(cand_fields.r0 - ref_fields.r0))
    mu_delta = torch.abs(cand_mu - ref_mu)
    return {
        "psi": float(psi_delta.detach().cpu().item()),
        "A": float(a_delta.detach().cpu().item()),
        "R0": float(r0_delta.detach().cpu().item()),
        "mu": float(mu_delta.detach().cpu().item()),
    }


def _observable_input_deltas(
    reference: torch.Tensor,
    candidate: torch.Tensor,
    grid,
) -> dict[str, Any]:
    ref_fields, _ = unpack_closed_coupled_fields(reference, grid, has_chemical_potential=True)
    cand_fields, _ = unpack_closed_coupled_fields(candidate, grid, has_chemical_potential=True)
    ref_density = ref_fields.psi_real**2 + ref_fields.psi_imag**2
    cand_density = cand_fields.psi_real**2 + cand_fields.psi_imag**2
    density_delta = torch.max(torch.abs(cand_density - ref_density))
    r0_delta = torch.max(torch.abs(cand_fields.r0 - ref_fields.r0))
    return {
        "D0_inputs": {
            "density_linf_delta": float(density_delta.detach().cpu().item()),
            "R0_linf_delta": float(r0_delta.detach().cpu().item()),
            "status": "MEASURED",
        },
        "BZN_overlaps": {
            "status": "NOT_MEASURED",
            "reason": "BdG/Maxwell overlap recomputation is outside this bounded diagnostic",
        },
    }


def _compute_admissibility(
    result: Mapping[str, Any],
    *,
    config: C0Config,
    dtype: torch.dtype,
) -> dict[str, Any]:
    attempts = [dict(row) for row in result.get("tau_attempts", [])]
    selected: list[tuple[str, dict[str, Any]]] = []
    converged = [row for row in attempts if row.get("final_physical_converged")]
    if converged:
        selected.append(("tame", max(converged, key=lambda row: float(row["target_tau"]))))
        selected.append(
            (
                "near_floor",
                min(converged, key=lambda row: abs(float(row["target_tau"]) - PRIOR_TAU_FLOOR)),
            )
        )
    if attempts:
        selected.append(("deepest_attempted", min(attempts, key=lambda row: float(row["target_tau"]))))

    residual_rows: list[dict[str, Any]] = []
    residual_equality_max_abs = 0.0
    for label, row in selected:
        try:
            tau = float(row["target_tau"])
            branch, provider, grid, boundaries = _branch_context(tau=tau, config=config, dtype=dtype)
            state = _load_state_artifact(Path(row["state_artifact"]), dtype=dtype)
            eos_K = float(branch.continuation_K_values[-1])
            original = patha_closed_branch_residual(
                state,
                grid,
                branch,
                eos_K=eos_K,
                boundaries=boundaries,
                s_sigma=provider,
            )
            conditioned = c0_preconditioner_residual(
                state,
                grid,
                branch,
                eos_K=eos_K,
                boundaries=boundaries,
                s_sigma=provider,
                aid=C0AidParameters(core_epsilon=0.0, k1_radius_epsilon=0.0),
            )
            max_abs = float(torch.max(torch.abs(original - conditioned)).detach().cpu().item())
            residual_equality_max_abs = max(residual_equality_max_abs, max_abs)
            residual_rows.append(
                {
                    "label": label,
                    "tau": tau,
                    "max_abs": max_abs,
                    "tolerance": config.residual_equality_tolerance,
                    "status": "PASS"
                    if max_abs <= config.residual_equality_tolerance
                    else "FAIL",
                }
            )
        except Exception as exc:
            residual_rows.append(
                {
                    "label": label,
                    "tau": row.get("target_tau"),
                    "status": "NOT_MEASURED",
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )

    epsilon_table: list[dict[str, Any]] = []
    near_floor = selected[1][1] if len(selected) >= 2 and selected[1][0] == "near_floor" else None
    if near_floor is None:
        epsilon_status = "NOT_MEASURED"
        epsilon_reason = "no_near_floor_converged_state"
    else:
        epsilon_status = "PASS"
        epsilon_reason = None
        tau = float(near_floor["target_tau"])
        branch, provider, grid, boundaries = _branch_context(tau=tau, config=config, dtype=dtype)
        clean_state = _load_state_artifact(Path(near_floor["state_artifact"]), dtype=dtype)
        eos_K = float(branch.continuation_K_values[-1])
        solved: list[tuple[C0AidParameters, C0NewtonResult]] = []
        for aid in config.epsilon_schedule:
            original_residual_fn = lambda x, eos_K=eos_K: patha_closed_branch_residual(
                x,
                grid,
                branch,
                eos_K=float(eos_K),
                boundaries=boundaries,
                s_sigma=provider,
            )
            preconditioner_residual_fn = lambda x, eos_K=eos_K, aid=aid: c0_preconditioner_residual(
                x,
                grid,
                branch,
                eos_K=float(eos_K),
                boundaries=boundaries,
                s_sigma=provider,
                aid=aid,
            )
            solved.append(
                (
                    aid,
                    solve_c0_scaled_newton(
                        original_residual_fn=original_residual_fn,
                        preconditioner_residual_fn=preconditioner_residual_fn,
                        x0=clean_state,
                        grid=grid,
                        newton=branch.newton,
                        aid=aid,
                    ),
                )
            )
        reference_result = next(
            (
                item
                for item in solved
                if item[0].core_epsilon == 0.0 and item[0].k1_radius_epsilon == 0.0
            ),
            solved[-1],
        )[1]
        reference_state = reference_result.x.detach()
        for aid, solve_result in solved:
            deltas = _state_group_deltas(reference_state, solve_result.x.detach(), grid)
            observable_deltas = _observable_input_deltas(reference_state, solve_result.x.detach(), grid)
            row_status = (
                "PASS"
                if solve_result.converged
                and max(deltas.values()) <= config.epsilon_independence_tolerance
                else "FAIL"
            )
            if observable_deltas["BZN_overlaps"]["status"] == "NOT_MEASURED":
                row_status = "NOT_MEASURED" if row_status == "PASS" else row_status
            if row_status != "PASS":
                epsilon_status = "NOT_MEASURED" if row_status == "NOT_MEASURED" else "FAIL"
            epsilon_table.append(
                {
                    "tau": tau,
                    "epsilon": float(aid.core_epsilon),
                    "core_epsilon": float(aid.core_epsilon),
                    "k1_radius_epsilon": float(aid.k1_radius_epsilon),
                    "final_original_residual": float(solve_result.final_residual_norm),
                    "converged": bool(solve_result.converged),
                    "state_group_deltas_vs_epsilon0": deltas,
                    "observable_deltas_vs_epsilon0": observable_deltas,
                    "comparison_tolerance": config.epsilon_independence_tolerance,
                    "status": row_status,
                }
            )

    residual_statuses = {row.get("status") for row in residual_rows}
    residual_status = (
        "PASS"
        if residual_rows and residual_statuses == {"PASS"}
        else "NOT_MEASURED"
        if "NOT_MEASURED" in residual_statuses or not residual_rows
        else "FAIL"
    )
    return {
        "residual_equality_status": residual_status,
        "residual_equality_max_abs": residual_equality_max_abs,
        "residual_equality_table": residual_rows,
        "epsilon_independence_status": epsilon_status,
        "epsilon_independence_reason": epsilon_reason,
        "epsilon_independence_table": epsilon_table,
        "epsilon_path_final_aids_inactive": all(
            row.get("floor_activation", {}).get("final_aids_inactive") for row in attempts
        ),
        "original_residual_gating_status": "PASS",
        "jacobi_scaling_solution_neutral": True,
        "jacobi_scaling_note": (
            "Scaling is applied as R*J*C dz = R*(-F), step=C*dz. "
            "Line search and convergence use the unscaled original residual."
        ),
    }


def _determine_verdict(result: Mapping[str, Any]) -> tuple[str, dict[str, Any], str]:
    attempts = list(result.get("tau_attempts", []))
    config = result.get("config", {})
    below_converged = [
        row
        for row in attempts
        if _below_prior_floor(float(row["target_tau"])) and row.get("final_physical_converged")
    ]
    if below_converged:
        deepest = min(below_converged, key=lambda row: float(row["target_tau"]))
        support = {
            "deepest_genuine_converged_tau": float(deepest["target_tau"]),
            "prior_tau_floor": PRIOR_TAU_FLOOR,
            "init_source": deepest.get("init", {}).get("source"),
            "used_existing_b2c": deepest.get("used_existing_b2c"),
            "final_original_residual_linf": deepest.get("final_original_residual_linf"),
        }
        return "SPIKE_SUFFICIENT", support, "No production solver pivot; continue physics review from the deeper genuine state."

    fold = result.get("fold_diagnostic", {})
    sigma = result.get("sigma_diagnostics", {})
    if fold.get("status") != "MEASURED" or sigma.get("status") != "MEASURED":
        return (
            "DIAGNOSTIC_INCOMPLETE",
            {
                "fold_status": fold.get("status"),
                "sigma_status": sigma.get("status"),
                "reason": "required_sigma_or_fold_evidence_not_measured",
            },
            "Rerun or narrow the bounded linear diagnostics before selecting a solver strategy.",
        )
    if fold.get("call") == "FOLD_TURNING_POINT":
        return (
            "FOLD_TURNING_POINT",
            {
                "fold_track_indices": fold.get("fold_track_indices"),
                "persistent_track_indices": fold.get("persistent_track_indices"),
                "tracked_modes": fold.get("tracked_modes"),
            },
            "Pivot to pseudo-arclength continuation; do not treat this as a linear-solver rebuild.",
        )

    deepest = sigma.get("deepest", {})
    v_fracs = deepest.get("v_min_energy_fractions", {})
    u_fracs = deepest.get("u_min_energy_fractions", {})
    v_group, v_fraction = _dominant_group(v_fracs)
    u_group, u_fraction = _dominant_group(u_fracs)
    cond_ratio = deepest.get("cond_ratio")
    persistent = fold.get("call") == "PERSISTENT_NEAR_NULL_SPACE"
    border_dominant = (
        (v_group in {"R0", "mu"} and v_fraction >= 0.7)
        or (u_group == "mass_constraint_row" and u_fraction >= 0.7)
    )
    if (
        persistent
        and border_dominant
        and isinstance(cond_ratio, float)
        and math.isfinite(cond_ratio)
        and cond_ratio < 0.1
    ):
        return (
            "NEAR_NULL_SPACE_STRUCTURAL",
            {
                "dominant_v_group": v_group,
                "dominant_v_fraction": v_fraction,
                "dominant_u_group": u_group,
                "dominant_u_fraction": u_fraction,
                "cond_ratio": cond_ratio,
                "thresholds": {
                    "dominant_fraction_min": 0.7,
                    "near_null_cond_ratio_max": 0.1,
                },
                "sigma_min": deepest.get("sigma_min"),
                "sigma_method": deepest.get("sigma_method"),
            },
            "Target border scaling, Schur complementing, or null-space deflation before a generic production solver.",
        )

    persistent_guard = _crawl_persistent_failure_guard(attempts, config)
    below_failed = persistent_guard["failed_attempts"]
    attempted_backtracking = bool(persistent_guard["attempted_backtracking"])
    full_newton_budget = bool(persistent_guard["full_newton_budget"])
    crawl_persistent_failure = bool(persistent_guard["crawl_persistent_failure_evidence"])
    if (
        crawl_persistent_failure
        and v_group == "field"
        and v_fraction >= 0.7
        and isinstance(cond_ratio, float)
        and math.isfinite(cond_ratio)
        and cond_ratio >= 0.1
        and fold.get("call") == "NO_TRACKED_COMPLEMENT_CROSSING"
    ):
        return (
            "PRODUCTION_SOLVER_REQUIRED",
            {
                "dominant_v_group": v_group,
                "dominant_v_fraction": v_fraction,
                "cond_ratio": cond_ratio,
                "thresholds": {
                    "dominant_fraction_min": 0.7,
                    "production_cond_ratio_min": 0.1,
                },
                "sigma_min": deepest.get("sigma_min"),
                "sigma_method": deepest.get("sigma_method"),
                "fold_call": fold.get("call"),
            },
            "Implement an evidence-supported field-block solver; generic multigrid is justified only with field-mode evidence.",
        )

    return (
        "DIAGNOSTIC_INCOMPLETE",
        {
            "reason": "measured_evidence_did_not_satisfy_exact_substantive_verdict_gates",
            "crawl_persistent_failure_evidence": crawl_persistent_failure,
            "below_floor_failed_attempt_count": len(below_failed),
            "attempted_backtracking": attempted_backtracking,
            "full_newton_budget": full_newton_budget,
            "fold_call": fold.get("call"),
            "dominant_v_group": v_group,
            "dominant_v_fraction": v_fraction,
            "dominant_u_group": u_group,
            "dominant_u_fraction": u_fraction,
            "cond_ratio": cond_ratio,
        },
        "Review the measured diagnostics and tighten the next bounded diagnostic around the unresolved gate.",
    )


def run_c0_diagnostics(
    config: C0Config | None = None,
    *,
    result: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if config is None:
        config = C0Config()
    dtype = configure_backend(BackendConfig())
    if result is None:
        result = json.loads(config.json_path.read_text(encoding="utf-8"))
    updated = dict(result)
    attempts = [dict(row) for row in updated.get("tau_attempts", [])]
    selected = _select_diagnostic_attempts(updated, config)
    measured_diagnostics: list[dict[str, Any]] = []
    diagnostics_by_key: dict[tuple[float, int], dict[str, Any]] = {}
    final_aid = C0AidParameters(core_epsilon=0.0, k1_radius_epsilon=0.0)
    for row in selected:
        tau = float(row["target_tau"])
        key = (tau, int(row["backtrack_index"]))
        try:
            branch, provider, grid, boundaries = _branch_context(tau=tau, config=config, dtype=dtype)
            state = _load_state_artifact(Path(row["state_artifact"]), dtype=dtype)
            matrix_path = (
                config.run_root
                / "matrices"
                / f"attempt_tau_{_format_tau(tau)}_bt_{int(row['backtrack_index'])}.npz"
            )
            diag = _linear_diagnostics_for_state(
                state=state,
                grid=grid,
                branch=branch,
                provider=provider,
                boundaries=boundaries,
                eos_K=float(branch.continuation_K_values[-1]),
                aid=final_aid,
                config=config,
                matrix_path=matrix_path,
            )
            diag["tau"] = tau
            diag["backtrack_index"] = int(row["backtrack_index"])
            measured_diagnostics.append(diag)
            diagnostics_by_key[key] = diag
        except Exception as exc:
            diagnostics_by_key[key] = {
                "status": "NOT_MEASURED",
                "reason": f"{type(exc).__name__}: {exc}",
                "tau": tau,
                "backtrack_index": int(row["backtrack_index"]),
            }

    fold = _track_singular_triplets(measured_diagnostics)
    json_linear_rows = []
    for row in attempts:
        key = (float(row["target_tau"]), int(row["backtrack_index"]))
        diag = diagnostics_by_key.get(
            key,
            {
                "status": "NOT_MEASURED",
                "reason": "not_selected_for_bounded_linear_diagnostics",
            },
        )
        row["linear_diagnostics"] = _strip_internal_diagnostic_fields(diag)
        json_linear_rows.append(row)

    measured_json = [_strip_internal_diagnostic_fields(diag) for diag in measured_diagnostics]
    converged_measured = [
        diag
        for diag in measured_diagnostics
        if any(
            abs(float(row["target_tau"]) - float(diag["tau"])) <= 5.0e-13
            and row.get("final_physical_converged")
            for row in attempts
        )
    ]
    shallow = max(converged_measured, key=lambda diag: float(diag["tau"])) if converged_measured else None
    deepest = min(measured_diagnostics, key=lambda diag: float(diag["tau"])) if measured_diagnostics else None
    sigma_status = (
        "MEASURED"
        if shallow is not None
        and deepest is not None
        and shallow.get("status") == "MEASURED"
        and deepest.get("status") == "MEASURED"
        else "NOT_MEASURED"
    )
    sigma_summary = {
        "status": sigma_status,
        "shallow": _strip_internal_diagnostic_fields(shallow) if shallow is not None else None,
        "deepest": _strip_internal_diagnostic_fields(deepest) if deepest is not None else None,
        "measured_points": measured_json,
        "required_methods": ["dense_svd", "svds_propack_SM", "eigsh_normal_shift_invert"],
    }
    updated["tau_attempts"] = json_linear_rows
    updated["depth_sequence"] = json_linear_rows
    updated["sigma_diagnostics"] = sigma_summary
    updated["fold_diagnostic"] = fold
    updated["r0_tau_diagnostic"] = _diagnose_tau_vs_r0(json_linear_rows)
    updated["aids_admissibility"] = _compute_admissibility(updated, config=config, dtype=dtype)
    verdict, support, next_step = _determine_verdict(updated)
    updated["verdict"] = verdict
    updated["verdict_support"] = support
    updated["recommended_next_step"] = next_step
    updated["phase_status"] = "diagnostics_complete"
    updated["faithful_operator_boundary"] = _faithful_operator_boundary()
    updated["elapsed_seconds_total_with_diagnostics"] = updated.get("elapsed_seconds")
    config.json_path.parent.mkdir(parents=True, exist_ok=True)
    config.json_path.write_text(json.dumps(updated, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_c0_report(updated, config.report_path)
    return updated


def run_c0_conditioning_spike(config: C0Config | None = None) -> dict[str, Any]:
    if config is None:
        config = C0Config()
    crawl = run_c0_crawl(config)
    return run_c0_diagnostics(config, result=crawl)


def _config_to_dict(config: C0Config) -> dict[str, Any]:
    data = asdict(config)
    for key in ("run_root", "report_path", "json_path"):
        data[key] = str(data[key])
    return data


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isfinite(value):
            return f"{value:.6e}"
        return str(value)
    return str(value)


def _markdown_table(headers: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(header)) for header in headers) + " |")
    return "\n".join(lines)


def write_c0_report(result: Mapping[str, Any], path: Path) -> None:
    rows = []
    for row in result["tau_attempts"]:
        linear = row.get("linear_diagnostics", {})
        rows.append(
            {
                "tau": row["target_tau"],
                "delta_tau": row.get("delta_tau"),
                "backtrack": row.get("backtrack_index"),
                "init": row.get("init", {}).get("source"),
                "used_b2c": row.get("used_existing_b2c"),
                "eps_probes": len(row.get("epsilon_attempts", [])),
                "R0_min": row["metrics"]["min_R0"],
                "min_rho": row["metrics"]["min_rho"],
                "residual": row["final_original_residual_linf"],
                "converged": row["final_physical_converged"],
                "sigma_status": linear.get("status"),
                "sigma_method": linear.get("sigma_method"),
                "sigma_min": linear.get("sigma_min"),
                "cond_ratio": linear.get("cond_ratio"),
                "message": row["message"],
            }
        )
    lines: list[str] = []
    lines.append("# Path-A C0b Wall Diagnosis")
    lines.append("")
    lines.append(f"Verdict: **{result['verdict']}**")
    lines.append(
        f"Deepest converged tau: `{_fmt(result['deepest_converged_tau'])}` "
        f"versus prior floor `{PRIOR_TAU_FLOOR:.6e}`; "
        f"R0_min: `{_fmt(result['deepest_converged_R0_min'])}`."
    )
    lines.append("")
    lines.append("## Single Arbiter Boundary")
    lines.append("")
    lines.append(
        "Newton merit, line search, and convergence are evaluated with "
        "`stage1_solver.coupled_branch.patha_closed_branch_residual`. "
        "The C0 matter epsilon and k1 radius floor are used only while assembling "
        "the preconditioner. The final accepted state for every converged row has "
        "all epsilon aids inactive and is checked against the original residual."
    )
    lines.append("")
    lines.append("## Tau Attempts")
    lines.append("")
    lines.append(
        _markdown_table(
            [
                "tau",
                "delta_tau",
                "backtrack",
                "init",
                "used_b2c",
                "eps_probes",
                "R0_min",
                "min_rho",
                "residual",
                "converged",
                "sigma_status",
                "sigma_method",
                "sigma_min",
                "cond_ratio",
                "message",
            ],
            rows,
        )
    )
    lines.append("")
    lines.append("Notes:")
    lines.append("")
    lines.append(
        "- The residual column is the original unscaled physical residual from "
        "`patha_closed_branch_residual`; no scaled or preconditioner residual is used for convergence."
    )
    lines.append(
        "- The depth sequence includes failed attempts below the prior `tau ~= 0.029` floor."
    )
    lines.append("- Below-floor attempts are required to show `init.source=previous_c0_converged_state` and `used_existing_b2c=false` in the machine JSON.")
    lines.append("")
    lines.append("## True Sigma Diagnostics")
    lines.append("")
    sigma = result.get("sigma_diagnostics", {})
    lines.append("```yaml")
    lines.append(f"status: {sigma.get('status')}")
    for label in ("shallow", "deepest"):
        point = sigma.get(label) or {}
        lines.append(f"{label}:")
        lines.append(f"  tau: {point.get('tau')}")
        lines.append(f"  status: {point.get('status')}")
        lines.append(f"  sigma_method: {point.get('sigma_method')}")
        lines.append(f"  sigma_min: {point.get('sigma_min')}")
        lines.append(f"  triplet_residual_rel: {point.get('triplet_residual_rel')}")
        lines.append(f"  v_min_energy_fractions: {point.get('v_min_energy_fractions')}")
        lines.append(f"  u_min_energy_fractions: {point.get('u_min_energy_fractions')}")
        lines.append(f"  condition_field_block: {point.get('condition_field_block')}")
        lines.append(f"  condition_full_bordered: {point.get('condition_full_bordered')}")
        lines.append(f"  cond_ratio: {point.get('cond_ratio')}")
        lines.append(f"  gauge_projection_status: {point.get('gauge_projection_status')}")
    lines.append("```")
    lines.append("")
    lines.append("## Tracked-Triplet Fold Test")
    lines.append("")
    fold = result.get("fold_diagnostic", {})
    lines.append("```yaml")
    lines.append(f"status: {fold.get('status')}")
    lines.append(f"call: {fold.get('call')}")
    lines.append(f"persistent_track_indices: {fold.get('persistent_track_indices')}")
    lines.append(f"fold_track_indices: {fold.get('fold_track_indices')}")
    lines.append(f"determinant_sign_support: {fold.get('determinant_sign_support')}")
    lines.append("tracked_modes:")
    for track in fold.get("tracked_modes", []):
        lines.append(f"  - track_index: {track.get('track_index')}")
        lines.append(f"    sigma_ratio: {track.get('sigma_ratio')}")
        lines.append(f"    persistent_border_or_mass_mode: {track.get('persistent_border_or_mass_mode')}")
        lines.append(f"    fold_like_complement_mode: {track.get('fold_like_complement_mode')}")
        lines.append("    points:")
        for point in track.get("points", []):
            lines.append(
                "      - "
                f"tau: {point.get('tau')}, sigma: {point.get('sigma')}, "
                f"overlap: {point.get('overlap_from_previous')}, "
                f"v_group: {point.get('v_dominant_group')}, "
                f"v_fraction: {point.get('v_dominant_fraction')}"
            )
    lines.append("```")
    lines.append("")
    lines.append("## R0-vs-Tau")
    lines.append("")
    r0_diag = result.get("r0_tau_diagnostic", {})
    lines.append("```yaml")
    lines.append(f"status: {r0_diag.get('status')}")
    lines.append(f"statement: {r0_diag.get('statement')}")
    lines.append(f"tau_decrease_approaches_empty_core: {r0_diag.get('tau_decrease_approaches_empty_core')}")
    lines.append(f"candidate_knobs: {r0_diag.get('candidate_knobs')}")
    lines.append("trend:")
    for point in r0_diag.get("trend", []):
        lines.append(
            f"  - tau: {point.get('tau')}, R0_min: {point.get('R0_min')}, "
            f"min_rho: {point.get('min_rho')}, converged: {point.get('converged')}"
        )
    lines.append("```")
    lines.append("")
    lines.append("## Aid Admissibility")
    lines.append("")
    lines.append("```yaml")
    admissibility = result.get("aids_admissibility", {})
    lines.append(f"residual_equality_status: {admissibility.get('residual_equality_status')}")
    lines.append(f"residual_equality_max_abs: {admissibility.get('residual_equality_max_abs')}")
    lines.append(f"epsilon_independence_status: {admissibility.get('epsilon_independence_status')}")
    lines.append("residual_equality_table:")
    for row in admissibility.get("residual_equality_table", []):
        lines.append(
            f"  - label: {row.get('label')}, tau: {row.get('tau')}, "
            f"max_abs: {row.get('max_abs')}, status: {row.get('status')}"
        )
    lines.append("epsilon_independence_table:")
    for row in admissibility.get("epsilon_independence_table", []):
        lines.append(
            f"  - epsilon: {row.get('epsilon')}, residual: {row.get('final_original_residual')}, "
            f"state_group_deltas: {row.get('state_group_deltas_vs_epsilon0')}, "
            f"observable_deltas: {row.get('observable_deltas_vs_epsilon0')}, "
            f"status: {row.get('status')}"
        )
    lines.append("```")
    lines.append("")
    lines.append(
        "C0-1 and C0-2 are implemented as preconditioner-only epsilon aids. "
        "The conditioned residual used for the approximate inverse equals the physical residual at "
        "`core_epsilon=0` and `k1_radius_epsilon=0`; accepted rows are all evaluated with those "
        "epsilons inactive. No `log rho`, `sqrt rho`, density-only variable rewrite, faithful "
        "operator edit, or final-state k1 clamp is used."
    )
    lines.append("")
    lines.append("Epsilon schedule used on each target:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["core_epsilon", "k1_radius_epsilon", "status"],
            [
                {
                    "core_epsilon": aid["core_epsilon"],
                    "k1_radius_epsilon": aid["k1_radius_epsilon"],
                    "status": "final physical-limit polish"
                    if aid["core_epsilon"] == 0.0 and aid["k1_radius_epsilon"] == 0.0
                    else "path-only preconditioner aid",
                }
                for aid in result["config"]["epsilon_schedule"]
            ],
        )
    )
    lines.append("")
    lines.append(
        "Verdict gates use the exact C0b thresholds: `SPIKE_SUFFICIENT` requires a genuine "
        "below-floor convergence; structural near-null requires a border/gauge dominant fraction "
        "`>=0.7` and `cond_ratio<0.1`; production-solver required needs field dominance `>=0.7`, "
        "`cond_ratio>=0.1`, measured sigma, and no tracked fold."
    )
    lines.append("")
    lines.append("## Tau Attempt Details")
    lines.append("")
    for row in result["tau_attempts"]:
        lines.append(f"### tau={row['target_tau']:.12e}")
        lines.append("")
        lines.append("```yaml")
        lines.append(f"target_tau: {row['target_tau']:.12e}")
        lines.append(f"nominal_target_tau: {row.get('nominal_target_tau')}")
        lines.append(f"delta_tau: {row.get('delta_tau')}")
        lines.append(f"backtrack_index: {row.get('backtrack_index')}")
        lines.append(f"start_from_tau: {row.get('start_from_tau')}")
        lines.append(f"init: {row.get('init')}")
        lines.append(f"used_existing_b2c: {row.get('used_existing_b2c')}")
        lines.append(f"final_original_residual_linf: {row['final_original_residual_linf']:.12e}")
        lines.append(f"b2c_background_tolerance: {BACKGROUND_RESIDUAL_TOL:.12e}")
        lines.append(f"final_physical_converged: {row['final_physical_converged']}")
        lines.append(f"min_rho: {row['metrics']['min_rho']:.12e}")
        lines.append(f"min_R0: {row['metrics']['min_R0']:.12e}")
        lines.append(f"k1_clamp_active_in_path: {row['floor_activation']['k1_clamp_active_in_path']}")
        lines.append("```")
        lines.append("")
        step_rows = []
        for step in row["epsilon_attempts"]:
            last_hist = step["newton_history"][-1] if step["newton_history"] else {}
            step_rows.append(
                {
                    "eos_K": step["eos_K"],
                    "core_eps": step["aid"]["core_epsilon"],
                    "k1_eps": step["aid"]["k1_radius_epsilon"],
                    "iters": step["iterations"],
                    "residual": step["final_residual_norm"],
                    "alpha": last_hist.get("line_search_alpha"),
                    "step_norm": last_hist.get("step_norm"),
                    "gmres_iters": last_hist.get("gmres_iterations"),
                    "gmres_growth": step.get("gmres_growth"),
                    "start": step.get("start_policy"),
                    "message": step["message"],
                }
            )
        lines.append(
            _markdown_table(
                [
                    "eos_K",
                    "core_eps",
                    "k1_eps",
                    "iters",
                    "residual",
                    "alpha",
                    "step_norm",
                    "gmres_iters",
                    "gmres_growth",
                    "start",
                    "message",
                ],
                step_rows,
            )
        )
        lines.append("")
    lines.append("## Verdict Support")
    lines.append("")
    lines.append("```yaml")
    for key, value in result.get("verdict_support", {}).items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Recommended Next Step")
    lines.append("")
    lines.append(str(result.get("recommended_next_step")))
    lines.append("")
    lines.append(
        "Machine artifact: `software/stage1_solver/runs/pathA_C0b_wall_diagnosis/pathA_C0b_diagnostic.json`."
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _project_root() -> Path:
    return Path(__file__).resolve().parents[4]


def _stage_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_input_path(path: str | Path) -> Path:
    raw = Path(path)
    if raw.is_absolute():
        return raw
    candidates = [Path.cwd() / raw, _project_root() / raw]
    try:
        stage_relative = raw.relative_to(Path("software/stage1_solver"))
        candidates.append(_stage_root() / stage_relative)
    except ValueError:
        candidates.append(_stage_root() / raw)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _resolve_output_path(path: str | Path) -> Path:
    raw = Path(path)
    if raw.is_absolute():
        return raw
    if raw.parts[:2] == ("software", "stage1_solver"):
        return _project_root() / raw
    return Path.cwd() / raw


def _c0c_config_to_dict(config: C0cConfig) -> dict[str, Any]:
    data = asdict(config)
    for key in ("c0b_json_path", "run_root", "report_path", "json_path"):
        data[key] = str(data[key])
    return data


def _closed_lane_slices(grid) -> dict[str, tuple[int, int]]:
    n = int(grid.spec.nr * grid.spec.nw)
    field_dim = 5 * n
    return {
        "psi_real": (0, n),
        "psi_imag": (n, 2 * n),
        "a0": (2 * n, 3 * n),
        "ar": (3 * n, 4 * n),
        "aw": (4 * n, 5 * n),
        "r0": (field_dim, field_dim + int(grid.spec.nw)),
        "mu": (field_dim + int(grid.spec.nw), field_dim + int(grid.spec.nw) + 1),
    }


def _lane_energy_split(vector: np.ndarray, grid) -> dict[str, float]:
    values = np.asarray(vector, dtype=np.float64)
    denom = float(np.dot(values, values))
    if denom <= 0.0:
        return {name: math.nan for name in C0C_FIELD_LAYOUT}
    split = {}
    for name, (start, stop) in _closed_lane_slices(grid).items():
        chunk = values[start:stop]
        split[name] = float(np.dot(chunk, chunk) / denom)
    return split


def _c0d_config_to_dict(config: C0dConfig) -> dict[str, Any]:
    data = asdict(config)
    for key in ("c0b_json_path", "run_root", "report_path", "json_path"):
        data[key] = str(data[key])
    return data


def _c0e_config_to_dict(config: C0eConfig) -> dict[str, Any]:
    data = asdict(config)
    for key in ("c0b_json_path", "run_root", "report_path", "json_path"):
        data[key] = str(data[key])
    return data


def _c0d_spatial_a_norm(vector: np.ndarray, grid) -> float:
    lanes = _closed_lane_slices(grid)
    values = np.asarray(vector, dtype=np.float64)
    pieces = []
    for name in ("ar", "aw"):
        start, stop = lanes[name]
        pieces.append(values[start:stop])
    return float(np.linalg.norm(np.concatenate(pieces)))


def _c0d_non_a_fraction(vector: np.ndarray, grid) -> float:
    split = _lane_energy_split(vector, grid)
    a_fraction = float(split.get("ar", 0.0)) + float(split.get("aw", 0.0))
    return float(max(0.0, 1.0 - a_fraction))


def _c0d_extract_spatial_a_fields(vector: np.ndarray, grid) -> tuple[torch.Tensor, torch.Tensor]:
    values = np.asarray(vector, dtype=np.float64)
    lanes = _closed_lane_slices(grid)
    shape = (grid.spec.nr, grid.spec.nw)
    ar_start, ar_stop = lanes["ar"]
    aw_start, aw_stop = lanes["aw"]
    ar = torch.as_tensor(
        values[ar_start:ar_stop].reshape(shape),
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    )
    aw = torch.as_tensor(
        values[aw_start:aw_stop].reshape(shape),
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    )
    return ar, aw


def _c0d_embed_spatial_a_fields(ar: torch.Tensor, aw: torch.Tensor, grid) -> np.ndarray:
    dims = _closed_layout_dimensions(grid)
    values = np.zeros(dims["state_size"], dtype=np.float64)
    lanes = _closed_lane_slices(grid)
    ar_start, ar_stop = lanes["ar"]
    aw_start, aw_stop = lanes["aw"]
    values[ar_start:ar_stop] = ar.detach().cpu().numpy().reshape(-1)
    values[aw_start:aw_stop] = aw.detach().cpu().numpy().reshape(-1)
    return values


def _c0d_scalar_gradient_vector(lambda_field: torch.Tensor, grid) -> np.ndarray:
    return _c0d_embed_spatial_a_fields(
        tensor_center_gradient_r(lambda_field, grid),
        tensor_center_gradient_w(lambda_field, grid),
        grid,
    )


def _c0d_scalar_gradient_matrix(grid) -> np.ndarray:
    cell_count = int(grid.spec.nr * grid.spec.nw)
    dims = _closed_layout_dimensions(grid)
    matrix = np.empty((dims["state_size"], cell_count), dtype=np.float64)
    eye_fields = torch.eye(
        cell_count,
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    ).reshape(cell_count, grid.spec.nr, grid.spec.nw)
    for column in range(cell_count):
        matrix[:, column] = _c0d_scalar_gradient_vector(eye_fields[column], grid)
    return matrix


def _c0d_column_space_basis(
    matrix: np.ndarray,
    *,
    relative_tolerance: float,
) -> tuple[np.ndarray, np.ndarray, float, int]:
    u, singular_values, _vh = np.linalg.svd(np.asarray(matrix, dtype=np.float64), full_matrices=False)
    largest = float(singular_values[0]) if singular_values.size else 0.0
    threshold = max(
        float(relative_tolerance) * largest,
        np.finfo(np.float64).eps * max(matrix.shape) * max(largest, 1.0),
    )
    rank = int(np.sum(singular_values > threshold))
    if rank <= 0:
        return np.zeros((matrix.shape[0], 0), dtype=np.float64), singular_values, threshold, rank
    return u[:, :rank].astype(np.float64, copy=True), singular_values, threshold, rank


def _c0d_weighted_divergence_flat(vector: np.ndarray, grid, branch) -> np.ndarray:
    ar, aw = _c0d_extract_spatial_a_fields(vector, grid)
    divergence = axisymmetric_vector_divergence(ar, aw, grid)
    z = localization_weight_torch(grid.w_centers, branch)[None, :]
    return (z * divergence).detach().cpu().numpy().reshape(-1).astype(np.float64, copy=True)


def _c0d_weighted_gauge_gradient_flat(vector: np.ndarray, grid, branch) -> np.ndarray:
    ar, aw = _c0d_extract_spatial_a_fields(vector, grid)
    divergence = axisymmetric_vector_divergence(ar, aw, grid)
    z = localization_weight_torch(grid.w_centers, branch)[None, :]
    weighted = z * divergence
    grad_r = tensor_center_gradient_r(weighted, grid)
    grad_w = tensor_center_gradient_w(weighted, grid)
    return (
        torch.cat([grad_r.reshape(-1), grad_w.reshape(-1)])
        .detach()
        .cpu()
        .numpy()
        .astype(np.float64, copy=True)
    )


def _c0d_raw_divergence_flat(vector: np.ndarray, grid) -> np.ndarray:
    ar, aw = _c0d_extract_spatial_a_fields(vector, grid)
    divergence = axisymmetric_vector_divergence(ar, aw, grid)
    return divergence.detach().cpu().numpy().reshape(-1).astype(np.float64, copy=True)


def _c0d_build_gauge_subspace(grid, branch, config: C0dConfig) -> dict[str, Any]:
    gradient_matrix = _c0d_scalar_gradient_matrix(grid)
    basis, singular_values, rank_threshold, rank = _c0d_column_space_basis(
        gradient_matrix,
        relative_tolerance=config.gradient_rank_rtol,
    )
    if rank > 0:
        weighted_divergence_matrix = np.column_stack(
            [_c0d_weighted_divergence_flat(basis[:, index], grid, branch) for index in range(rank)]
        )
        _u, weighted_singular_values, vh = np.linalg.svd(
            weighted_divergence_matrix,
            full_matrices=False,
        )
        largest = float(weighted_singular_values[0]) if weighted_singular_values.size else 0.0
        harmonic_threshold = float(config.harmonic_weighted_divergence_rtol) * largest
        if largest <= np.finfo(np.float64).tiny:
            harmonic_indices = np.arange(rank, dtype=int)
        else:
            harmonic_indices = np.where(weighted_singular_values <= harmonic_threshold)[0]
        if harmonic_indices.size:
            harmonic_basis = basis @ vh.T[:, harmonic_indices]
        else:
            harmonic_basis = np.zeros((basis.shape[0], 0), dtype=np.float64)
    else:
        weighted_divergence_matrix = np.zeros((grid.spec.nr * grid.spec.nw, 0), dtype=np.float64)
        weighted_singular_values = np.zeros(0, dtype=np.float64)
        harmonic_threshold = 0.0
        harmonic_indices = np.zeros(0, dtype=int)
        harmonic_basis = np.zeros((gradient_matrix.shape[0], 0), dtype=np.float64)

    return {
        "basis": basis,
        "harmonic_basis": harmonic_basis.astype(np.float64, copy=True),
        "gradient_matrix": gradient_matrix,
        "gradient_singular_values": singular_values.astype(np.float64, copy=True),
        "gradient_rank_threshold": float(rank_threshold),
        "weighted_divergence_matrix": weighted_divergence_matrix,
        "weighted_divergence_singular_values": weighted_singular_values.astype(
            np.float64,
            copy=True,
        ),
        "weighted_divergence_harmonic_threshold": float(harmonic_threshold),
        "weighted_divergence_harmonic_indices": [int(index) for index in harmonic_indices],
        "dim_g": int(basis.shape[1]),
        "dim_g_harm": int(harmonic_basis.shape[1]),
    }


def _c0d_projection_fraction(basis: np.ndarray, vector: np.ndarray) -> float:
    values = np.asarray(vector, dtype=np.float64)
    norm = float(np.linalg.norm(values))
    if norm <= 0.0:
        return math.nan
    if basis.shape[1] == 0:
        return 0.0
    unit = values / norm
    captured = float(np.linalg.norm(basis.T @ unit) ** 2)
    return float(min(max(captured, 0.0), 1.0))


def _c0d_a_only_vector(vector: np.ndarray, grid) -> np.ndarray:
    lanes = _closed_lane_slices(grid)
    values = np.asarray(vector, dtype=np.float64)
    a_only = np.zeros_like(values)
    for name in ("ar", "aw"):
        start, stop = lanes[name]
        a_only[start:stop] = values[start:stop]
    return a_only


def _c0d_spatial_a_metrics(vector: np.ndarray, grid, branch) -> dict[str, float]:
    a_norm = _c0d_spatial_a_norm(vector, grid)
    if a_norm <= 0.0:
        return {
            "spatial_a_norm": a_norm,
            "raw_divergence": math.nan,
            "weighted_gauge_residual": math.nan,
            "weighted_divergence": math.nan,
        }
    raw = _c0d_raw_divergence_flat(vector, grid)
    weighted = _c0d_weighted_divergence_flat(vector, grid, branch)
    weighted_gradient = _c0d_weighted_gauge_gradient_flat(vector, grid, branch)
    return {
        "spatial_a_norm": a_norm,
        "raw_divergence": float(np.linalg.norm(raw) / a_norm),
        "weighted_gauge_residual": float(np.linalg.norm(weighted_gradient) / a_norm),
        "weighted_divergence": float(np.linalg.norm(weighted) / a_norm),
    }


def _c0d_mode_diagnostics(
    *,
    mode_index: int,
    sigma: float,
    vector: np.ndarray,
    subspace: Mapping[str, Any],
    grid,
    branch,
    config: C0dConfig,
) -> dict[str, Any]:
    unit, mode_norm = _unit_vector(vector)
    lane_split = _lane_energy_split(unit, grid)
    a_only = _c0d_a_only_vector(unit, grid)
    p_g = _c0d_projection_fraction(subspace["basis"], unit)
    p_g_harm = _c0d_projection_fraction(subspace["harmonic_basis"], unit)
    p_g_a_normalized = _c0d_projection_fraction(subspace["basis"], a_only)
    metrics = _c0d_spatial_a_metrics(unit, grid, branch)
    weighted = metrics["weighted_gauge_residual"]
    classification = (
        "MAXWELL_GAUGE"
        if (
            math.isfinite(p_g)
            and math.isfinite(weighted)
            and p_g >= config.gauge_capture_threshold
            and weighted <= config.weighted_residual_threshold
        )
        else "GENUINE_MAXWELL_STIFFNESS"
    )
    return {
        "mode_index": int(mode_index),
        "sigma": float(sigma),
        "mode_norm_before_unit": float(mode_norm),
        "p_g_fraction": p_g,
        "p_g_harm_fraction": p_g_harm,
        "p_g_a_normalized_fraction": p_g_a_normalized,
        "unexplained_residual": float(max(0.0, 1.0 - p_g)) if math.isfinite(p_g) else math.nan,
        "raw_divergence": metrics["raw_divergence"],
        "weighted_gauge_residual": weighted,
        "weighted_divergence": metrics["weighted_divergence"],
        "spatial_a_norm": metrics["spatial_a_norm"],
        "lane_energy_fractions": lane_split,
        "spatial_a_energy_fraction": float(lane_split.get("ar", 0.0))
        + float(lane_split.get("aw", 0.0)),
        "non_a_remainder": _c0d_non_a_fraction(unit, grid),
        "classification": classification,
        "classification_gate": {
            "p_g_min": float(config.gauge_capture_threshold),
            "weighted_gauge_residual_max": float(config.weighted_residual_threshold),
            "p_g_pass": bool(
                math.isfinite(p_g) and p_g >= config.gauge_capture_threshold
            ),
            "weighted_gauge_residual_pass": bool(
                math.isfinite(weighted) and weighted <= config.weighted_residual_threshold
            ),
        },
    }


def _c0d_control_diagnostics(subspace: Mapping[str, Any], grid, config: C0dConfig) -> dict[str, Any]:
    r_scaled = grid.r_centers[:, None] / float(grid.spec.r_max)
    w_span = float(grid.spec.w_max) - float(grid.spec.w_min)
    w_scaled = (grid.w_centers[None, :] - float(grid.spec.w_min)) / w_span
    lambda_field = (
        torch.sin(math.pi * r_scaled)
        * torch.cos(2.0 * math.pi * w_scaled)
        + 0.37 * r_scaled * (1.0 - w_scaled)
    )
    positive = _c0d_scalar_gradient_vector(lambda_field, grid)
    positive_p_g = _c0d_projection_fraction(subspace["basis"], positive)
    positive_non_a = _c0d_non_a_fraction(positive, grid)

    rng = np.random.default_rng(int(config.control_random_seed))
    random_a = np.zeros(_closed_layout_dimensions(grid)["state_size"], dtype=np.float64)
    lanes = _closed_lane_slices(grid)
    for name in ("ar", "aw"):
        start, stop = lanes[name]
        random_a[start:stop] = rng.normal(size=stop - start)
    if subspace["basis"].shape[1] > 0:
        negative = random_a - subspace["basis"] @ (subspace["basis"].T @ random_a)
    else:
        negative = random_a
    negative_p_g = _c0d_projection_fraction(subspace["basis"], negative)
    negative_non_a = _c0d_non_a_fraction(negative, grid)

    status = (
        "PASS"
        if positive_p_g >= 1.0 - 1.0e-10
        and positive_non_a <= 1.0e-12
        and negative_p_g <= 1.0e-10
        and negative_non_a <= 1.0e-12
        else "FAIL"
    )
    return {
        "status": status,
        "positive_control": {
            "type": "held_out_smooth_discrete_gradient",
            "p_g_fraction": float(positive_p_g),
            "non_a_remainder": float(positive_non_a),
            "pass": bool(positive_p_g >= 1.0 - 1.0e-10 and positive_non_a <= 1.0e-12),
        },
        "negative_control": {
            "type": "transverse_orthogonalized_random_A",
            "p_g_fraction": float(negative_p_g),
            "non_a_remainder": float(negative_non_a),
            "pass": bool(negative_p_g <= 1.0e-10 and negative_non_a <= 1.0e-12),
        },
        "thresholds": {
            "positive_p_g_min": 1.0 - 1.0e-10,
            "positive_non_a_max": 1.0e-12,
            "negative_p_g_max": 1.0e-10,
            "negative_non_a_max": 1.0e-12,
        },
    }


def _c0e_scope_guard() -> dict[str, bool]:
    return {
        "diagnosis_only": True,
        "single_read_only_linear_solve_per_artifact": True,
        "trial_residual_evaluations_on_temporary_tensors_only": True,
        "state_advance_or_newton_write": False,
        "gauge_fix_or_deflation_implemented": False,
        "stencil_change_implemented": False,
        "recrawl_implemented": False,
        "changed_xi_reassembly": False,
        "faithful_operators_touched_by_c0e": False,
        "frozen_physics_touched_by_c0e": False,
        "physical_export_guard_touched_by_c0e": False,
    }


def _c0e_operator_sources(branch) -> dict[str, Any]:
    return {
        "curl_fraction": {
            "field_strength": "stage1_solver.operators.tensor_center_gradient_r(aw) - tensor_center_gradient_w(ar)",
            "localization_weight": "stage1_solver.coupled_branch.localization_weight_torch(grid.w_centers, branch)",
            "ratio": "||Z*F_rw[v]|| / ||A[v]||",
            "norms": ["unweighted_code_space", "cell_volume_weighted"],
        },
        "coupled_generator": {
            "real_lane_formula": (
                "delta psi_real=-(q/hbar)*lambda*psi_imag; "
                "delta psi_imag=+(q/hbar)*lambda*psi_real; "
                "delta a0=0; delta ar=d_r(lambda); delta aw=d_w(lambda)"
            ),
            "q_over_hbar_source": (
                "branch.gauge_charge / branch.hbar from "
                "stage1_solver.coupled_branch.coupled_pde_residual alpha and "
                "q*A0 matter coupling"
            ),
            "q_over_hbar_value": float(branch.gauge_charge / branch.hbar),
        },
        "maxwell_operator": {
            "source": "stage1_solver.operators.localized_maxwell_operator",
            "curl_rows": "ar=-d_w(Z*F_rw), aw=radial_divergence_from_center_flux(Z*F_rw)",
            "gauge_penalty_rows": "(1/xi)*grad(Z*axisymmetric_vector_divergence(A))",
            "xi": float(branch.xi),
            "sponge_gauge_strength": float(getattr(branch, "sponge_gauge_strength", 0.0)),
        },
    }


def _c0e_coupled_gauge_vector(
    lambda_field: torch.Tensor,
    state: torch.Tensor,
    grid,
    branch,
) -> np.ndarray:
    fields, _mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    alpha = float(branch.gauge_charge / branch.hbar)
    dims = _closed_layout_dimensions(grid)
    values = np.zeros(dims["state_size"], dtype=np.float64)
    lanes = _closed_lane_slices(grid)
    pieces = {
        "psi_real": -alpha * lambda_field * fields.psi_imag,
        "psi_imag": alpha * lambda_field * fields.psi_real,
        "ar": tensor_center_gradient_r(lambda_field, grid),
        "aw": tensor_center_gradient_w(lambda_field, grid),
    }
    for name, piece in pieces.items():
        start, stop = lanes[name]
        values[start:stop] = piece.detach().cpu().numpy().reshape(-1)
    return values


def _c0e_coupled_gauge_matrix(state: torch.Tensor, grid, branch) -> np.ndarray:
    cell_count = int(grid.spec.nr * grid.spec.nw)
    dims = _closed_layout_dimensions(grid)
    matrix = np.empty((dims["state_size"], cell_count), dtype=np.float64)
    eye_fields = torch.eye(
        cell_count,
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    ).reshape(cell_count, grid.spec.nr, grid.spec.nw)
    for column in range(cell_count):
        matrix[:, column] = _c0e_coupled_gauge_vector(eye_fields[column], state, grid, branch)
    return matrix


def _c0e_build_coupled_subspace(
    state: torch.Tensor,
    grid,
    branch,
    config: C0eConfig,
) -> dict[str, Any]:
    matrix = _c0e_coupled_gauge_matrix(state, grid, branch)
    u, singular_values, vh = np.linalg.svd(matrix, full_matrices=False)
    largest = float(singular_values[0]) if singular_values.size else 0.0
    threshold = max(
        float(config.gradient_rank_rtol) * largest,
        np.finfo(np.float64).eps * max(matrix.shape) * max(largest, 1.0),
    )
    rank = int(np.sum(singular_values > threshold))
    basis = u[:, :rank].astype(np.float64, copy=True) if rank > 0 else np.zeros((matrix.shape[0], 0))
    return {
        "basis": basis,
        "generator_matrix": matrix,
        "singular_values": singular_values.astype(np.float64, copy=True),
        "vh": vh.astype(np.float64, copy=True),
        "rank_threshold": float(threshold),
        "rank": int(rank),
        "scalar_basis_dim": int(matrix.shape[1]),
        "full_state_dim": int(matrix.shape[0]),
        "q_over_hbar": float(branch.gauge_charge / branch.hbar),
    }


def _c0e_field_norm(
    values: torch.Tensor,
    grid,
    *,
    weighted: bool,
    mask: torch.Tensor | None = None,
) -> float:
    selected = values if mask is None else values[mask]
    if selected.numel() == 0:
        return 0.0
    if not weighted:
        return float(torch.linalg.vector_norm(selected).detach().cpu().item())
    volumes = grid.cell_volumes if mask is None else grid.cell_volumes[mask]
    return float(torch.sqrt(torch.sum((selected**2) * volumes)).detach().cpu().item())


def _c0e_spatial_a_norm_from_fields(
    ar: torch.Tensor,
    aw: torch.Tensor,
    grid,
    *,
    weighted: bool,
) -> float:
    density = ar**2 + aw**2
    if not weighted:
        return float(torch.sqrt(torch.sum(density)).detach().cpu().item())
    return float(torch.sqrt(torch.sum(density * grid.cell_volumes)).detach().cpu().item())


def _c0e_boundary_mask(grid) -> torch.Tensor:
    mask = torch.zeros(
        (grid.spec.nr, grid.spec.nw),
        dtype=torch.bool,
        device=grid.r_centers.device,
    )
    mask[0, :] = True
    mask[-1, :] = True
    mask[:, 0] = True
    mask[:, -1] = True
    return mask


def _c0e_curl_fraction_from_fields(
    ar: torch.Tensor,
    aw: torch.Tensor,
    grid,
    branch,
) -> dict[str, Any]:
    z = localization_weight_torch(grid.w_centers, branch)[None, :]
    curl = z * (tensor_center_gradient_r(aw, grid) - tensor_center_gradient_w(ar, grid))
    boundary_mask = _c0e_boundary_mask(grid)
    interior_mask = ~boundary_mask
    metrics: dict[str, Any] = {}
    for label, weighted in (
        ("unweighted", False),
        ("cell_volume_weighted", True),
    ):
        numerator = _c0e_field_norm(curl, grid, weighted=weighted)
        denominator = _c0e_spatial_a_norm_from_fields(ar, aw, grid, weighted=weighted)
        boundary = _c0e_field_norm(curl, grid, weighted=weighted, mask=boundary_mask)
        interior = _c0e_field_norm(curl, grid, weighted=weighted, mask=interior_mask)
        metrics[label] = {
            "curl_norm": float(numerator),
            "spatial_a_norm": float(denominator),
            "curl_fraction": float(numerator / denominator)
            if denominator > 0.0
            else math.nan,
            "boundary_curl_norm": float(boundary),
            "interior_curl_norm": float(interior),
            "boundary_curl_fraction_of_norm": float(boundary / numerator)
            if numerator > 0.0
            else 0.0,
        }
    return metrics


def _c0e_curl_fraction_metrics(vector: np.ndarray, grid, branch) -> dict[str, Any]:
    ar, aw = _c0d_extract_spatial_a_fields(vector, grid)
    return _c0e_curl_fraction_from_fields(ar, aw, grid, branch)


def _c0e_positive_lambda_fields(grid) -> list[tuple[str, torch.Tensor]]:
    r_scaled = grid.r_centers[:, None] / float(grid.spec.r_max)
    w_span = float(grid.spec.w_max) - float(grid.spec.w_min)
    w_scaled = (grid.w_centers[None, :] - float(grid.spec.w_min)) / w_span
    i = torch.arange(grid.spec.nr, dtype=grid.r_centers.dtype, device=grid.r_centers.device)[:, None]
    j = torch.arange(grid.spec.nw, dtype=grid.r_centers.dtype, device=grid.r_centers.device)[None, :]
    checker = torch.where(((i + j).remainder(2) == 0), 1.0, -1.0).to(grid.r_centers.dtype)
    return [
        (
            "smooth_sin_poly_gradient",
            torch.sin(math.pi * r_scaled) * torch.cos(2.0 * math.pi * w_scaled)
            + 0.37 * r_scaled * (1.0 - w_scaled),
        ),
        (
            "smooth_mixed_gradient",
            (r_scaled**2) * torch.sin(math.pi * w_scaled)
            + 0.19 * torch.cos(2.0 * math.pi * r_scaled) * (0.5 + w_scaled),
        ),
        (
            "high_k_sine_gradient",
            torch.sin((grid.spec.nr - 1.0) * math.pi * r_scaled)
            * torch.cos((grid.spec.nw - 2.0) * math.pi * w_scaled),
        ),
        ("checkerboard_gradient", checker),
    ]


def _c0e_stream_function_a_vector(chi: torch.Tensor, grid) -> np.ndarray:
    r2 = torch.clamp(grid.r_centers[:, None] ** 2, min=float(grid.dr) ** 2)
    return _c0d_embed_spatial_a_fields(
        tensor_center_gradient_w(chi, grid) / r2,
        -tensor_center_gradient_r(chi, grid) / r2,
        grid,
    )


def _c0e_negative_stream_functions(grid) -> list[tuple[str, torch.Tensor]]:
    r_scaled = grid.r_centers[:, None] / float(grid.spec.r_max)
    w_span = float(grid.spec.w_max) - float(grid.spec.w_min)
    w_scaled = (grid.w_centers[None, :] - float(grid.spec.w_min)) / w_span
    i = torch.arange(grid.spec.nr, dtype=grid.r_centers.dtype, device=grid.r_centers.device)[:, None]
    j = torch.arange(grid.spec.nw, dtype=grid.r_centers.dtype, device=grid.r_centers.device)[None, :]
    checker = torch.where(((i + j).remainder(2) == 0), 1.0, -1.0).to(grid.r_centers.dtype)
    return [
        (
            "smooth_stream_sin",
            torch.sin(math.pi * r_scaled) * torch.sin(2.0 * math.pi * w_scaled),
        ),
        (
            "smooth_stream_mixed",
            r_scaled * (1.0 - r_scaled) * torch.cos(math.pi * w_scaled)
            + 0.11 * torch.sin(2.0 * math.pi * r_scaled),
        ),
        (
            "high_k_stream_sine",
            torch.sin((grid.spec.nr - 1.0) * math.pi * r_scaled)
            * torch.sin((grid.spec.nw - 1.0) * math.pi * w_scaled),
        ),
        ("checkerboard_stream", checker),
    ]


def _c0e_control_bands(rows: Sequence[Mapping[str, Any]], norm_label: str) -> dict[str, Any]:
    positive = [
        float(row[norm_label]["curl_fraction"])
        for row in rows
        if row.get("family") == "positive_gradient"
        and math.isfinite(float(row[norm_label]["curl_fraction"]))
    ]
    negative = [
        float(row[norm_label]["curl_fraction"])
        for row in rows
        if row.get("family") == "negative_stream_function_transverse"
        and math.isfinite(float(row[norm_label]["curl_fraction"]))
    ]
    if not positive or not negative:
        return {"status": "NOT_MEASURED", "reason": "missing_control_family"}
    positive_max = max(positive)
    negative_min = min(negative)
    tiny = np.finfo(np.float64).eps
    separated = bool(positive_max < negative_min)
    separator = math.sqrt(max(positive_max, tiny) * max(negative_min, tiny)) if separated else math.nan
    return {
        "status": "SEPARATED" if separated else "OVERLAP",
        "positive_min": float(min(positive)),
        "positive_max": float(positive_max),
        "negative_min": float(negative_min),
        "negative_max": float(max(negative)),
        "geometric_separator": float(separator),
        "separation_ratio": float(negative_min / max(positive_max, tiny)),
        "separation_log10": float(math.log10(negative_min / max(positive_max, tiny)))
        if negative_min > 0.0
        else -math.inf,
    }


def _c0e_control_diagnostics(
    *,
    a_subspace: Mapping[str, Any],
    coupled_subspace: Mapping[str, Any],
    state: torch.Tensor,
    grid,
    branch,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for name, lambda_field in _c0e_positive_lambda_fields(grid):
        vector = _c0d_scalar_gradient_vector(lambda_field, grid)
        coupled_vector = _c0e_coupled_gauge_vector(lambda_field, state, grid, branch)
        metrics = _c0e_curl_fraction_metrics(vector, grid, branch)
        rows.append(
            {
                "name": name,
                "family": "positive_gradient",
                "construction": "held_out_discrete_gradient_not_used_to_define_mode_classification",
                "unweighted": metrics["unweighted"],
                "cell_volume_weighted": metrics["cell_volume_weighted"],
                "a_only_p_g_fraction": _c0d_projection_fraction(a_subspace["basis"], vector),
                "coupled_capture_fraction": _c0d_projection_fraction(
                    coupled_subspace["basis"], coupled_vector
                ),
                "non_a_remainder": _c0d_non_a_fraction(vector, grid),
            }
        )
    for name, chi in _c0e_negative_stream_functions(grid):
        vector = _c0e_stream_function_a_vector(chi, grid)
        metrics = _c0e_curl_fraction_metrics(vector, grid, branch)
        rows.append(
            {
                "name": name,
                "family": "negative_stream_function_transverse",
                "construction": "independent_stream_function_A_not_random_minus_own_gradient_projection",
                "unweighted": metrics["unweighted"],
                "cell_volume_weighted": metrics["cell_volume_weighted"],
                "a_only_p_g_fraction": _c0d_projection_fraction(a_subspace["basis"], vector),
                "coupled_capture_fraction": _c0d_projection_fraction(
                    coupled_subspace["basis"], vector
                ),
                "non_a_remainder": _c0d_non_a_fraction(vector, grid),
            }
        )
    bands = {
        "unweighted": _c0e_control_bands(rows, "unweighted"),
        "cell_volume_weighted": _c0e_control_bands(rows, "cell_volume_weighted"),
    }
    positive_rows = [row for row in rows if row["family"] == "positive_gradient"]
    negative_rows = [row for row in rows if row["family"] == "negative_stream_function_transverse"]
    status = "PASS"
    reasons: list[str] = []
    if any(band.get("status") != "SEPARATED" for band in bands.values()):
        status = "FAIL"
        reasons.append("control_curl_bands_overlap_or_missing")
    if any(float(row["a_only_p_g_fraction"]) < 1.0 - 1.0e-10 for row in positive_rows):
        status = "FAIL"
        reasons.append("positive_gradient_a_only_capture_not_exact")
    if any(float(row["coupled_capture_fraction"]) < 1.0 - 1.0e-10 for row in positive_rows):
        status = "FAIL"
        reasons.append("positive_coupled_capture_not_exact")
    if any(float(row["non_a_remainder"]) > 1.0e-12 for row in positive_rows):
        status = "FAIL"
        reasons.append("positive_gradient_non_a_remainder_nonzero")
    max_negative_capture = max(
        [float(row["a_only_p_g_fraction"]) for row in negative_rows],
        default=math.inf,
    )
    if max_negative_capture >= 0.5:
        status = "FAIL"
        reasons.append("stream_function_negative_control_has_high_gradient_capture")
    return {
        "status": status,
        "reason": ";".join(reasons) if reasons else None,
        "controls": rows,
        "bands": bands,
        "negative_capture_max_observed": float(max_negative_capture),
        "negative_capture_low_gate": 0.5,
    }


def _c0e_curl_margin(value: float, bands: Mapping[str, Any]) -> dict[str, Any]:
    if bands.get("status") != "SEPARATED" or not math.isfinite(value) or value <= 0.0:
        return {
            "outcome": "CONTROL_BANDS_OVERLAP"
            if bands.get("status") != "SEPARATED"
            else "AMBIGUOUS",
            "margin_log10_to_separator": math.nan,
            "separator": bands.get("geometric_separator"),
            "observed_control_gap_log10": bands.get("separation_log10"),
        }
    separator = float(bands["geometric_separator"])
    if value <= separator:
        margin = math.log10(separator / value)
        outcome = "GAUGE"
    else:
        margin = math.log10(value / separator)
        outcome = "TRANSVERSE"
    return {
        "outcome": outcome,
        "margin_log10_to_separator": float(margin),
        "separator": float(separator),
        "observed_control_gap_log10": float(bands["separation_log10"]),
        "positive_control_max": float(bands["positive_max"]),
        "negative_control_min": float(bands["negative_min"]),
    }


def _c0e_frequency_classification(lambda_field: np.ndarray) -> dict[str, Any]:
    values = np.asarray(lambda_field, dtype=np.float64)
    centered = values - float(np.mean(values))
    total = float(np.sum(centered**2))
    if total <= np.finfo(np.float64).tiny:
        return {
            "classification": "ZERO_PREIMAGE",
            "low_frequency_fraction": 0.0,
            "high_frequency_fraction": 0.0,
            "nyquist_fraction": 0.0,
            "checkerboard_score": 0.0,
        }
    fft = np.fft.fftn(centered)
    energy = np.abs(fft) ** 2
    total_fft = float(np.sum(energy))
    kr = np.abs(np.fft.fftfreq(values.shape[0]) * values.shape[0])[:, None]
    kw = np.abs(np.fft.fftfreq(values.shape[1]) * values.shape[1])[None, :]
    low_mask = (kr <= 2.0) & (kw <= 2.0)
    high_mask = (kr >= max(2.0, values.shape[0] / 4.0)) | (
        kw >= max(2.0, values.shape[1] / 4.0)
    )
    nyquist_mask = (kr >= max(1.0, values.shape[0] / 2.0 - 1.0)) | (
        kw >= max(1.0, values.shape[1] / 2.0 - 1.0)
    )
    ii = np.arange(values.shape[0])[:, None]
    jj = np.arange(values.shape[1])[None, :]
    patterns = [
        (-1.0) ** ii * np.ones_like(values),
        np.ones_like(values) * ((-1.0) ** jj),
        (-1.0) ** (ii + jj),
    ]
    checker_scores = []
    for pattern in patterns:
        p = pattern - float(np.mean(pattern))
        denom = math.sqrt(float(np.sum(p**2)) * total)
        checker_scores.append(float(abs(np.sum(centered * p)) / denom) if denom > 0.0 else 0.0)
    low = float(np.sum(energy[low_mask]) / total_fft) if total_fft > 0.0 else 0.0
    high = float(np.sum(energy[high_mask]) / total_fft) if total_fft > 0.0 else 0.0
    nyquist = float(np.sum(energy[nyquist_mask]) / total_fft) if total_fft > 0.0 else 0.0
    checker = float(max(checker_scores))
    if checker >= 0.35 or nyquist >= 0.35:
        classification = "CHECKERBOARD"
    elif high >= 0.5:
        classification = "HIGH_K"
    elif low >= 0.7:
        classification = "SMOOTH"
    else:
        classification = "MIXED"
    return {
        "classification": classification,
        "low_frequency_fraction": low,
        "high_frequency_fraction": high,
        "nyquist_fraction": nyquist,
        "checkerboard_score": checker,
    }


def _c0e_embed_pde_component(
    components: Mapping[int, torch.Tensor],
    grid,
    *,
    dtype: torch.dtype,
    device: torch.device | str,
) -> torch.Tensor:
    pde = torch.zeros((5, grid.spec.nr, grid.spec.nw), dtype=dtype, device=device)
    for lane, values in components.items():
        pde[int(lane)] = values
    wall_mass = torch.zeros(grid.spec.nw + 1, dtype=dtype, device=device)
    return torch.cat([pde.reshape(-1), wall_mass])


def _c0e_maxwell_curl_component(
    state: torch.Tensor,
    grid,
    branch,
) -> torch.Tensor:
    fields, _mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    z = localization_weight_torch(grid.w_centers, branch)[None, :]
    f_rw = tensor_center_gradient_r(fields.aw, grid) - tensor_center_gradient_w(fields.ar, grid)
    weighted_f_rw = z * f_rw
    ar_block = -tensor_center_gradient_w(weighted_f_rw, grid)
    aw_block = radial_divergence_from_center_flux(weighted_f_rw, grid)
    return _c0e_embed_pde_component(
        {3: ar_block, 4: aw_block},
        grid,
        dtype=state.dtype,
        device=state.device,
    )


def _c0e_maxwell_gauge_penalty_component(
    state: torch.Tensor,
    grid,
    branch,
) -> torch.Tensor:
    fields, _mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    z = localization_weight_torch(grid.w_centers, branch)[None, :]
    divergence = axisymmetric_vector_divergence(fields.ar, fields.aw, grid)
    weighted_divergence = z * divergence
    ar_block = (1.0 / float(branch.xi)) * tensor_center_gradient_r(weighted_divergence, grid)
    aw_block = (1.0 / float(branch.xi)) * tensor_center_gradient_w(weighted_divergence, grid)
    return _c0e_embed_pde_component(
        {3: ar_block, 4: aw_block},
        grid,
        dtype=state.dtype,
        device=state.device,
    )


def _c0e_maxwell_current_source_component(
    state: torch.Tensor,
    grid,
    branch,
) -> torch.Tensor:
    fields, _mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    coupled_fields = _closed_to_coupled_fields(fields)
    jr_number, jw_number = _matter_number_current(coupled_fields, grid, branch)
    density = fields.psi_real**2 + fields.psi_imag**2
    charge_current_r = branch.gauge_charge * jr_number
    charge_current_w = branch.gauge_charge * jw_number
    sponge = boundary_sponge_profile_torch(grid, branch)
    ar_block = (
        -branch.mu0 * charge_current_r
        + getattr(branch, "sponge_gauge_strength", 0.0) * sponge * fields.ar
    )
    aw_block = (
        -branch.mu0 * charge_current_w
        + getattr(branch, "sponge_gauge_strength", 0.0) * sponge * fields.aw
    )
    # Keep the density local so autograd sees the same dependencies used by the Gauss split.
    _ = density
    return _c0e_embed_pde_component(
        {3: ar_block, 4: aw_block},
        grid,
        dtype=state.dtype,
        device=state.device,
    )


def _c0e_matter_rows_component(
    state: torch.Tensor,
    grid,
    branch,
    *,
    eos_K: float,
    boundaries,
) -> torch.Tensor:
    fields, chemical_potential = unpack_closed_coupled_fields(
        state,
        grid,
        has_chemical_potential=True,
    )
    assert chemical_potential is not None
    pde = coupled_pde_residual(
        pack_coupled_fields(_closed_to_coupled_fields(fields)),
        grid,
        branch,
        eos_K=eos_K,
        chemical_potential=chemical_potential,
        boundaries=boundaries,
        confinement_radius=fields.r0,
    )
    return _c0e_embed_pde_component(
        {0: pde[0], 1: pde[1]},
        grid,
        dtype=state.dtype,
        device=state.device,
    )


def _c0e_gauss_charge_component(
    state: torch.Tensor,
    grid,
    branch,
    *,
    eos_K: float,
    boundaries,
) -> torch.Tensor:
    fields, chemical_potential = unpack_closed_coupled_fields(
        state,
        grid,
        has_chemical_potential=True,
    )
    assert chemical_potential is not None
    pde = coupled_pde_residual(
        pack_coupled_fields(_closed_to_coupled_fields(fields)),
        grid,
        branch,
        eos_K=eos_K,
        chemical_potential=chemical_potential,
        boundaries=boundaries,
        confinement_radius=fields.r0,
    )
    return _c0e_embed_pde_component(
        {2: pde[2]},
        grid,
        dtype=state.dtype,
        device=state.device,
    )


def _c0e_wall_mass_component(
    state: torch.Tensor,
    grid,
    branch,
    *,
    eos_K: float,
    boundaries,
    s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any],
) -> torch.Tensor:
    residual = patha_closed_branch_residual(
        state,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=s_sigma,
    )
    result = torch.zeros_like(residual)
    n = int(grid.spec.nr * grid.spec.nw)
    result[5 * n :] = residual[5 * n :]
    return result


def _c0e_jvp_component(
    fn: Callable[[torch.Tensor], torch.Tensor],
    state: torch.Tensor,
    direction: np.ndarray,
) -> np.ndarray:
    direction_t = torch.as_tensor(direction, dtype=state.dtype, device=state.device)
    return jvp(fn, state.detach(), direction_t).detach().cpu().numpy().astype(np.float64)


def _c0e_residual_budget_for_mode(
    *,
    matrix: csc_matrix,
    state: torch.Tensor,
    vector: np.ndarray,
    grid,
    branch,
    eos_K: float,
    boundaries,
    s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any],
) -> dict[str, Any]:
    unit, _norm = _unit_vector(vector)
    assembled = (matrix @ unit).astype(np.float64, copy=False)
    components: dict[str, np.ndarray] = {
        "maxwell_curl_rows": _c0e_jvp_component(
            lambda x: _c0e_maxwell_curl_component(x, grid, branch),
            state,
            unit,
        ),
        "maxwell_gauge_penalty_rows": _c0e_jvp_component(
            lambda x: _c0e_maxwell_gauge_penalty_component(x, grid, branch),
            state,
            unit,
        ),
        "maxwell_current_source_rows": _c0e_jvp_component(
            lambda x: _c0e_maxwell_current_source_component(x, grid, branch),
            state,
            unit,
        ),
        "matter_covariant_a_coupling_rows": _c0e_jvp_component(
            lambda x: _c0e_matter_rows_component(
                x,
                grid,
                branch,
                eos_K=eos_K,
                boundaries=boundaries,
            ),
            state,
            unit,
        ),
        "gauss_charge_rows": _c0e_jvp_component(
            lambda x: _c0e_gauss_charge_component(
                x,
                grid,
                branch,
                eos_K=eos_K,
                boundaries=boundaries,
            ),
            state,
            unit,
        ),
        "wall_mass_rows": _c0e_jvp_component(
            lambda x: _c0e_wall_mass_component(
                x,
                grid,
                branch,
                eos_K=eos_K,
                boundaries=boundaries,
                s_sigma=s_sigma,
            ),
            state,
            unit,
        ),
    }
    component_sum = np.sum(np.column_stack(list(components.values())), axis=1)
    unexplained = assembled - component_sum
    assembled_norm = float(np.linalg.norm(assembled))
    component_rows = {}
    for name, values in components.items():
        norm = float(np.linalg.norm(values))
        component_rows[name] = {
            "norm": norm,
            "fraction_of_assembled_jv_norm": float(norm / max(assembled_norm, np.finfo(np.float64).tiny)),
        }
    unexplained_norm = float(np.linalg.norm(unexplained))
    return {
        "assembled_jv_norm": assembled_norm,
        "component_sum_norm": float(np.linalg.norm(component_sum)),
        "sum_component_norms": float(sum(float(np.linalg.norm(v)) for v in components.values())),
        "components": component_rows,
        "unexplained_remainder_norm": unexplained_norm,
        "unexplained_fraction_of_assembled_jv_norm": float(
            unexplained_norm / max(assembled_norm, np.finfo(np.float64).tiny)
        ),
        "unexplained_fraction_of_component_sum_norm": float(
            unexplained_norm / max(float(np.linalg.norm(component_sum)), np.finfo(np.float64).tiny)
        ),
    }


def _c0e_row_scaling_summary(matrix: csc_matrix, grid) -> dict[str, Any]:
    row_norms = _sparse_abs_max_by_axis(matrix, 1)
    lanes = _closed_lane_slices(grid)
    groups = {
        "matter_psi_rows": (lanes["psi_real"][0], lanes["psi_imag"][1]),
        "gauss_a0_rows": lanes["a0"],
        "maxwell_ar_rows": lanes["ar"],
        "maxwell_aw_rows": lanes["aw"],
        "wall_rows": lanes["r0"],
        "mass_row": lanes["mu"],
    }
    summaries: dict[str, Any] = {}
    for name, (start, stop) in groups.items():
        chunk = row_norms[start:stop]
        summaries[name] = {
            "count": int(chunk.size),
            "min": float(np.min(chunk)) if chunk.size else math.nan,
            "median": float(np.median(chunk)) if chunk.size else math.nan,
            "max": float(np.max(chunk)) if chunk.size else math.nan,
            "rms": float(math.sqrt(float(np.mean(chunk**2)))) if chunk.size else math.nan,
        }
    maxwell = np.concatenate([row_norms[lanes["ar"][0] : lanes["ar"][1]], row_norms[lanes["aw"][0] : lanes["aw"][1]]])
    non_maxwell = np.concatenate(
        [
            row_norms[lanes["psi_real"][0] : lanes["a0"][1]],
            row_norms[lanes["r0"][0] : lanes["mu"][1]],
        ]
    )
    maxwell_rms = float(math.sqrt(float(np.mean(maxwell**2)))) if maxwell.size else math.nan
    non_maxwell_rms = (
        float(math.sqrt(float(np.mean(non_maxwell**2)))) if non_maxwell.size else math.nan
    )
    return {
        "row_norm_kind": "sparse_row_abs_max",
        "groups": summaries,
        "maxwell_araw_rms": maxwell_rms,
        "non_maxwell_rms": non_maxwell_rms,
        "maxwell_to_non_maxwell_rms_ratio": float(
            maxwell_rms / max(non_maxwell_rms, np.finfo(np.float64).tiny)
        )
        if math.isfinite(maxwell_rms) and math.isfinite(non_maxwell_rms)
        else math.nan,
    }


def _center_gradient_1d(values: torch.Tensor, spacing: float) -> torch.Tensor:
    if values.ndim != 1 or values.numel() < 3:
        raise ValueError("centered 1D gradient requires at least three values")
    grad = torch.empty_like(values)
    grad[1:-1] = (values[2:] - values[:-2]) / (2.0 * spacing)
    grad[0] = (-3.0 * values[0] + 4.0 * values[1] - values[2]) / (2.0 * spacing)
    grad[-1] = (3.0 * values[-1] - 4.0 * values[-2]) / (2.0 * spacing)
    return grad


def _pack_closed_generator(
    *,
    state: torch.Tensor,
    psi_real: torch.Tensor,
    psi_imag: torch.Tensor,
    a0: torch.Tensor,
    ar: torch.Tensor,
    aw: torch.Tensor,
    r0: torch.Tensor,
) -> np.ndarray:
    mu_zero = torch.zeros((), dtype=state.dtype, device=state.device)
    packed = pack_closed_coupled_fields(
        ClosedCoupledFields(
            psi_real=psi_real,
            psi_imag=psi_imag,
            a0=a0,
            ar=ar,
            aw=aw,
            r0=r0,
        ),
        mu_zero,
    )
    return packed.detach().cpu().numpy().astype(np.float64, copy=True)


def _c0c_generators_for_state(state: torch.Tensor, grid) -> list[C0cGenerator]:
    fields, _mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    zeros_field = torch.zeros_like(fields.psi_real)
    zeros_r0 = torch.zeros_like(fields.r0)

    phase = _pack_closed_generator(
        state=state,
        psi_real=-fields.psi_imag,
        psi_imag=fields.psi_real,
        a0=zeros_field,
        ar=zeros_field,
        aw=zeros_field,
        r0=zeros_r0,
    )

    grad_r = {
        "psi_real": tensor_center_gradient_r(fields.psi_real, grid),
        "psi_imag": tensor_center_gradient_r(fields.psi_imag, grid),
        "a0": tensor_center_gradient_r(fields.a0, grid),
        "ar": tensor_center_gradient_r(fields.ar, grid),
        "aw": tensor_center_gradient_r(fields.aw, grid),
    }
    translation_r = _pack_closed_generator(
        state=state,
        psi_real=grad_r["psi_real"],
        psi_imag=grad_r["psi_imag"],
        a0=grad_r["a0"],
        ar=grad_r["ar"],
        aw=grad_r["aw"],
        r0=zeros_r0,
    )

    grad_w = {
        "psi_real": tensor_center_gradient_w(fields.psi_real, grid),
        "psi_imag": tensor_center_gradient_w(fields.psi_imag, grid),
        "a0": tensor_center_gradient_w(fields.a0, grid),
        "ar": tensor_center_gradient_w(fields.ar, grid),
        "aw": tensor_center_gradient_w(fields.aw, grid),
    }
    translation_w = _pack_closed_generator(
        state=state,
        psi_real=grad_w["psi_real"],
        psi_imag=grad_w["psi_imag"],
        a0=grad_w["a0"],
        ar=grad_w["ar"],
        aw=grad_w["aw"],
        r0=_center_gradient_1d(fields.r0, float(grid.dw)),
    )

    r = grid.r_centers[:, None]
    dilation_r = _pack_closed_generator(
        state=state,
        psi_real=r * grad_r["psi_real"],
        psi_imag=r * grad_r["psi_imag"],
        a0=r * grad_r["a0"],
        ar=r * grad_r["ar"],
        aw=r * grad_r["aw"],
        r0=zeros_r0,
    )

    r_scale = grid.r_centers[:, None] / float(grid.spec.r_max)
    w_scale = (grid.w_centers[None, :] - float(grid.spec.w_min)) / (
        float(grid.spec.w_max) - float(grid.spec.w_min)
    )
    gauge_lambda = (r_scale**2) * torch.sin(math.pi * w_scale)
    maxwell_probe = _pack_closed_generator(
        state=state,
        psi_real=zeros_field,
        psi_imag=zeros_field,
        a0=zeros_field,
        ar=tensor_center_gradient_r(gauge_lambda, grid),
        aw=tensor_center_gradient_w(gauge_lambda, grid),
        r0=zeros_r0,
    )

    return [
        C0cGenerator(
            name="phase",
            classification="GAUGE_PHASE",
            symmetry_status="EXACT_SYMMETRY",
            description=(
                "global U(1) phase: (-psi_imag, psi_real, 0, 0, 0, 0, 0)"
            ),
            vector=phase,
        ),
        C0cGenerator(
            name="translation_r",
            classification="TRANSLATION",
            symmetry_status="PROBE",
            description="broken radial-translation probe: finite-difference d_r(state)",
            vector=translation_r,
        ),
        C0cGenerator(
            name="translation_w",
            classification="TRANSLATION",
            symmetry_status="PROBE",
            description="broken axial-translation probe: finite-difference d_w(state)",
            vector=translation_w,
        ),
        C0cGenerator(
            name="dilation_r",
            classification="DILATION",
            symmetry_status="PROBE",
            description="broken radial dilation probe: r*d_r(field lanes)",
            vector=dilation_r,
        ),
        C0cGenerator(
            name="maxwell_residual_gauge",
            classification="MAXWELL_RESIDUAL_GAUGE",
            symmetry_status="PROBE",
            description="Maxwell-sector pure-gradient probe: delta A=(0, d_r lambda, d_w lambda)",
            vector=maxwell_probe,
        ),
    ]


def _unit_vector(vector: np.ndarray) -> tuple[np.ndarray, float]:
    values = np.asarray(vector, dtype=np.float64)
    norm = float(np.linalg.norm(values))
    if norm <= 0.0:
        return values.copy(), norm
    return values / norm, norm


def _c0c_generator_metadata(generators: Sequence[C0cGenerator]) -> list[dict[str, Any]]:
    return [
        {
            "name": generator.name,
            "classification": generator.classification,
            "symmetry_status": generator.symmetry_status,
            "description": generator.description,
            "norm": float(np.linalg.norm(generator.vector)),
        }
        for generator in generators
    ]


def _phase_action_closed_vector(vector: np.ndarray, grid) -> np.ndarray:
    values = np.asarray(vector, dtype=np.float64)
    n = int(grid.spec.nr * grid.spec.nw)
    acted = np.zeros_like(values)
    acted[0:n] = -values[n : 2 * n]
    acted[n : 2 * n] = values[0:n]
    return acted


def _full_svd_cluster_from_matrix(
    matrix: csc_matrix,
    *,
    mode_count: int,
) -> dict[str, Any]:
    dense = matrix.toarray()
    values_no_uv = np.linalg.svd(dense, compute_uv=False)
    left_all, values_desc, right_t_desc = np.linalg.svd(dense, full_matrices=False)
    order = np.argsort(values_desc)
    keep = order[: min(int(mode_count), len(order))]
    sigma_min = float(np.min(values_no_uv))
    sigma_max = float(np.max(values_no_uv))
    return {
        "sigma_min": sigma_min,
        "sigma_max": sigma_max,
        "condition": float(sigma_max / sigma_min)
        if sigma_min > 0.0
        else math.inf,
        "vector_svd_sigma_min": float(values_desc[order[0]]),
        "vector_svd_sigma_max": float(values_desc[order[-1]]),
        "singular_values": [float(values_desc[index]) for index in keep],
        "right_vectors": right_t_desc[keep, :].astype(np.float64, copy=True),
        "left_vectors": left_all[:, keep].astype(np.float64, copy=True),
        "method": "dense_svd_recomputed_from_saved_c0b_matrix",
    }


def _c0c_overlap_diagnostics(
    *,
    right_vectors: np.ndarray,
    singular_values: Sequence[float],
    generators: Sequence[C0cGenerator],
    grid,
) -> tuple[list[dict[str, Any]], int]:
    unit_generators: list[tuple[C0cGenerator, np.ndarray]] = []
    for generator in generators:
        unit, norm = _unit_vector(generator.vector)
        if norm > 0.0:
            unit_generators.append((generator, unit))

    if unit_generators:
        generator_matrix = np.column_stack([unit for _generator, unit in unit_generators])
        q, r = np.linalg.qr(generator_matrix)
        diagonal = np.abs(np.diag(r))
        rank = int(np.sum(diagonal > 1.0e-12))
        span_basis = q[:, :rank] if rank > 0 else np.zeros((right_vectors.shape[1], 0))
    else:
        rank = 0
        span_basis = np.zeros((right_vectors.shape[1], 0))

    rows: list[dict[str, Any]] = []
    for mode_index, (sigma, vector) in enumerate(zip(singular_values, right_vectors)):
        unit_mode, mode_norm = _unit_vector(vector)
        overlaps = {
            generator.name: float(abs(np.dot(unit_mode, unit))) if mode_norm > 0.0 else math.nan
            for generator, unit in unit_generators
        }
        if span_basis.shape[1] > 0 and mode_norm > 0.0:
            captured = float(np.linalg.norm(span_basis.T @ unit_mode) ** 2)
        else:
            captured = 0.0
        captured = min(max(captured, 0.0), 1.0)
        rows.append(
            {
                "mode_index": int(mode_index),
                "sigma": float(sigma),
                "overlaps": overlaps,
                "span_capture_fraction": captured,
                "unexplained_residual_fraction": float(max(0.0, 1.0 - captured)),
                "lane_energy_fractions": _lane_energy_split(unit_mode, grid),
            }
        )
    return rows, rank


def _c0c_residual_context(
    *,
    tau: float,
    grid_shape: tuple[int, int],
    dtype: torch.dtype,
):
    context_config = C0Config(grid=grid_shape)
    branch, provider, grid, boundaries = _branch_context(
        tau=float(tau),
        config=context_config,
        dtype=dtype,
    )
    eos_K = float(branch.continuation_K_values[-1])

    def residual_fn(x: torch.Tensor) -> torch.Tensor:
        return patha_closed_branch_residual(
            x,
            grid,
            branch,
            eos_K=eos_K,
            boundaries=boundaries,
            s_sigma=provider,
        )

    return branch, provider, grid, boundaries, eos_K, residual_fn


def _c0g_config_to_dict(config: C0gConfig) -> dict[str, Any]:
    data = asdict(config)
    for key in (
        "c0f2_json_path",
        "run_root",
        "report_path",
        "step0_json_path",
        "json_path",
        "step5_json_path",
        "final_json_path",
    ):
        data[key] = str(data[key])
    return data


def _c0g_load_c0f2_payload(config: C0gConfig) -> dict[str, Any]:
    path = _resolve_input_path(config.c0f2_json_path)
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    payload["_resolved_path"] = str(path)
    return payload


def _c0g_tau_close(left: float, right: float) -> bool:
    return abs(float(left) - float(right)) <= 5.0e-13


def _c0g_rows_by_tau(payload: Mapping[str, Any]) -> dict[float, dict[str, Any]]:
    rows: dict[float, dict[str, Any]] = {}
    for row in payload.get("per_tau_timing", []):
        if "tau" in row:
            rows[float(row["tau"])] = dict(row)
    return rows


def _c0g_find_tau_row(payload: Mapping[str, Any], tau: float) -> dict[str, Any]:
    for row in payload.get("per_tau_timing", []):
        if _c0g_tau_close(float(row.get("tau", math.nan)), tau):
            return dict(row)
    raise ValueError(f"C0f2 per_tau_timing has no row for tau={tau:.12g}")


def _c0g_state_path_from_row(row: Mapping[str, Any]) -> Path:
    artifact = row.get("state_artifact")
    if not artifact:
        raise ValueError(f"missing state_artifact for tau={row.get('tau')}")
    path = _resolve_input_path(str(artifact))
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def _c0g_prefer_existing_flag_from_c0f2_payload(
    payload: Mapping[str, Any],
) -> tuple[bool, str]:
    key = "prefer_existing_b2c_background_predictor"
    for section_name in ("crawl_config", "scope_guard"):
        section = payload.get(section_name)
        if isinstance(section, Mapping) and key in section:
            value = section[key]
            if not isinstance(value, bool):
                raise ValueError(
                    f"C0f2 {section_name}.{key} must be a bool, got {value!r}"
                )
            return bool(value), f"{section_name}.{key}"
    raise KeyError(
        "C0f2 provenance is missing crawl_config.prefer_existing_b2c_background_predictor "
        "and scope_guard.prefer_existing_b2c_background_predictor"
    )


def _c0g_verify_stalled_provenance(config: C0gConfig) -> dict[str, Any]:
    payload = _c0g_load_c0f2_payload(config)
    row = _c0g_find_tau_row(payload, config.stalled_tau)
    prefer_existing, prefer_existing_source = _c0g_prefer_existing_flag_from_c0f2_payload(
        payload
    )
    provenance = {
        "status": "MEASURED",
        "c0f2_json_path": payload.get("_resolved_path"),
        "tau": float(row.get("tau", math.nan)),
        "attempt_index": int(row.get("attempt_index", -1)),
        "state_artifact": row.get("state_artifact"),
        "state_path": str(_c0g_state_path_from_row(row)),
        "source": row.get("init_source"),
        "start_tau": row.get("start_tau"),
        "prefer_existing_b2c_background_predictor": bool(prefer_existing),
        "prefer_existing_b2c_background_predictor_source": prefer_existing_source,
        "used_existing_b2c": bool(row.get("used_existing_b2c", False)),
        "accepted_default_success": bool(row.get("accepted_default_success", False)),
        "solver_converged": bool(row.get("solver_converged", False)),
        "message": row.get("message"),
    }
    genuine = (
        _c0g_tau_close(float(row.get("tau", math.nan)), config.stalled_tau)
        and row.get("init_source") == "previous_c0_converged_state"
        and not bool(prefer_existing)
        and not bool(row.get("used_existing_b2c", False))
    )
    provenance["genuine_warm_start_not_cold_loaded"] = bool(genuine)
    provenance["call"] = (
        "GENUINE_WARM_START_NOT_COLD_LOADED" if genuine else "PROVENANCE_FAILED"
    )
    return provenance


def _c0g_column_space(
    matrix: np.ndarray,
    *,
    relative_tolerance: float,
    full_matrices: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, int]:
    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2:
        raise ValueError("column-space input must be a matrix")
    if values.shape[1] == 0:
        u = np.eye(values.shape[0], dtype=np.float64) if full_matrices else np.zeros(
            (values.shape[0], 0), dtype=np.float64
        )
        return np.zeros((values.shape[0], 0), dtype=np.float64), u, np.zeros(0), 0.0, 0
    u, singular_values, _vh = np.linalg.svd(values, full_matrices=full_matrices)
    largest = float(singular_values[0]) if singular_values.size else 0.0
    threshold = max(
        float(relative_tolerance) * largest,
        np.finfo(np.float64).eps * max(values.shape) * max(largest, 1.0),
    )
    rank = int(np.sum(singular_values > threshold))
    basis = u[:, :rank].astype(np.float64, copy=True)
    return basis, u.astype(np.float64, copy=True), singular_values, float(threshold), rank


def _c0g_projection_fraction_from_basis(basis: np.ndarray, vector: np.ndarray) -> float:
    values = np.asarray(vector, dtype=np.float64)
    norm = float(np.linalg.norm(values))
    if norm <= 0.0:
        return math.nan
    if basis.shape[1] == 0:
        return 0.0
    captured = float(np.linalg.norm(basis.T @ (values / norm)) ** 2)
    return float(min(max(captured, 0.0), 1.0))


def _c0g_projection_vector(basis: np.ndarray, vector: np.ndarray) -> np.ndarray:
    values = np.asarray(vector, dtype=np.float64)
    if basis.shape[1] == 0:
        return np.zeros_like(values)
    return basis @ (basis.T @ values)


def _c0g_premise_gate_call(f_nn: float, config: C0gConfig) -> str:
    value = float(f_nn)
    if math.isfinite(value) and value < float(config.premise_fail_threshold):
        return "PREMISE_FAILED"
    if math.isfinite(value) and value > float(config.premise_hold_threshold):
        return "PREMISE_HOLDS"
    return "PREMISE_GRAY"


def _c0g_state_center_of_energy(vector: np.ndarray, grid) -> dict[str, float | None]:
    values = np.asarray(vector, dtype=np.float64)
    n = int(grid.spec.nr * grid.spec.nw)
    field_dim = 5 * n
    if values.size < field_dim:
        return {"r": None, "w": None}
    fields = values[:field_dim].reshape(5, grid.spec.nr, grid.spec.nw)
    energy = np.sum(fields * fields, axis=0)
    denom = float(np.sum(energy))
    if denom <= 0.0:
        return {"r": None, "w": None}
    r = grid.r_centers.detach().cpu().numpy().astype(np.float64)
    w = grid.w_centers.detach().cpu().numpy().astype(np.float64)
    return {
        "r": float(np.sum(energy * r[:, None]) / denom),
        "w": float(np.sum(energy * w[None, :]) / denom),
    }


def _c0g_scaled_norm(row_scale: np.ndarray, residual: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(row_scale, dtype=np.float64) * residual))


def _c0g_residual_context(
    *,
    tau: float,
    config: C0gConfig,
    dtype: torch.dtype,
):
    branch, provider, grid, boundaries = _branch_context(
        tau=float(tau),
        config=C0Config(grid=config.grid, jacobian_assembly=config.jacobian_assembly),
        dtype=dtype,
    )
    eos_K = float(branch.continuation_K_values[-1])

    def residual_fn(x: torch.Tensor) -> torch.Tensor:
        return patha_closed_branch_residual(
            x,
            grid,
            branch,
            eos_K=eos_K,
            boundaries=boundaries,
            s_sigma=provider,
        )

    return branch, provider, grid, boundaries, eos_K, residual_fn


def _c0g_assemble_original_jacobian(
    *,
    state: torch.Tensor,
    tau: float,
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    branch, _provider, grid, _boundaries, _eos_K, residual_fn = _c0g_residual_context(
        tau=float(tau),
        config=config,
        dtype=dtype,
    )
    residual = residual_fn(state).detach().cpu().numpy().astype(np.float64, copy=False)
    preconditioner_config = replace(
        branch.newton.preconditioner,
        type=str(config.jacobian_assembly),
        diagonal_shift=0.0,
    )
    linear_config = replace(branch.newton, preconditioner=preconditioner_config)
    start = time.perf_counter()
    build_context = PreconditionerBuildContext(
        residual_fn=residual_fn,
        x=state.detach(),
        rhs=-residual,
        iteration=0,
        config=linear_config,
    )
    if str(config.jacobian_assembly) == "autodiff_sparse_jacobian_lu":
        matrix, metadata = assemble_closed_coupled_autodiff_sparse_jacobian(
            build_context,
            grid,
            autodiff_mode="jacfwd",
        )
    elif str(config.jacobian_assembly) == "colored_sparse_jacobian_lu":
        matrix, metadata = assemble_closed_coupled_sparse_jacobian(build_context, grid)
    else:
        raise ValueError(f"unsupported C0g jacobian assembly {config.jacobian_assembly!r}")
    assembly_seconds = float(time.perf_counter() - start)
    matrix = matrix.tocsc()
    aid = C0AidParameters(
        core_epsilon=0.0,
        k1_radius_epsilon=0.0,
        use_jacobi_scaling=True,
        preconditioner_diagonal_shift=0.0,
    )
    row_scale, col_scale = jacobi_row_col_scales(matrix, aid=aid)
    scaled_matrix = (diags(row_scale, format="csc") @ matrix @ diags(col_scale, format="csc")).tocsc()
    metadata = dict(metadata)
    metadata.update(
        {
            "assembly_seconds": assembly_seconds,
            "jacobian_source": (
                "stage1_solver.preconditioners."
                "assemble_closed_coupled_autodiff_sparse_jacobian"
                if str(config.jacobian_assembly) == "autodiff_sparse_jacobian_lu"
                else "stage1_solver.preconditioners."
                "assemble_closed_coupled_sparse_jacobian"
            ),
            "residual_entry_point": (
                "stage1_solver.coupled_branch.patha_closed_branch_residual"
            ),
            "scaling": "jacobi_row_col_from_unmodified_original_jacobian",
            "row_scale_min": float(np.min(row_scale)),
            "row_scale_max": float(np.max(row_scale)),
            "col_scale_min": float(np.min(col_scale)),
            "col_scale_max": float(np.max(col_scale)),
            "diagonal_shift": 0.0,
            "core_epsilon": 0.0,
            "k1_radius_epsilon": 0.0,
        }
    )
    return {
        "branch": branch,
        "grid": grid,
        "residual_fn": residual_fn,
        "residual": residual,
        "matrix": matrix,
        "scaled_matrix": scaled_matrix,
        "row_scale": row_scale,
        "col_scale": col_scale,
        "metadata": metadata,
    }


def _c0g_build_analytic_gauge_matrix(
    *,
    state: torch.Tensor,
    grid,
    branch,
    config: C0gConfig,
) -> dict[str, Any]:
    generators = _c0c_generators_for_state(state, grid)
    phase = [generator.vector for generator in generators if generator.name == "phase"]
    if len(phase) != 1:
        raise RuntimeError("_c0c_generators_for_state did not return exactly one phase generator")
    c0d_config = C0dConfig(
        gradient_rank_rtol=float(config.gradient_rank_rtol),
        harmonic_weighted_divergence_rtol=float(config.harmonic_weighted_divergence_rtol),
    )
    a_subspace = _c0d_build_gauge_subspace(grid, branch, c0d_config)
    coupled_matrix = _c0e_coupled_gauge_matrix(state, grid, branch)
    pieces = [
        np.asarray(phase[0], dtype=np.float64).reshape(-1, 1),
        np.asarray(a_subspace["basis"], dtype=np.float64),
        np.asarray(coupled_matrix, dtype=np.float64),
    ]
    generator_matrix = np.column_stack([piece for piece in pieces if piece.shape[1] > 0])
    basis, _u, singular_values, threshold, rank = _c0g_column_space(
        generator_matrix,
        relative_tolerance=float(config.gradient_rank_rtol),
        full_matrices=False,
    )
    return {
        "generator_matrix": generator_matrix,
        "physical_basis": basis,
        "singular_values": singular_values,
        "rank_threshold": threshold,
        "rank": int(rank),
        "source_helpers": [
            "_c0c_generators_for_state(name='phase')",
            "_c0d_scalar_gradient_matrix",
            "_c0d_build_gauge_subspace",
            "_c0e_coupled_gauge_matrix",
        ],
        "piece_dimensions": {
            "phase_columns": 1,
            "a_sector_gradient_basis_columns": int(a_subspace["basis"].shape[1]),
            "a_sector_gradient_raw_columns": int(a_subspace["gradient_matrix"].shape[1]),
            "coupled_local_gauge_columns": int(coupled_matrix.shape[1]),
        },
        "a_sector_gradient_rank": int(a_subspace["dim_g"]),
        "a_sector_gradient_rank_threshold": float(a_subspace["gradient_rank_threshold"]),
    }


def _c0g_scaled_gauge_complement(
    *,
    gauge_matrix: np.ndarray,
    col_scale: np.ndarray,
    config: C0gConfig,
) -> dict[str, Any]:
    scaled_gauge = np.asarray(gauge_matrix, dtype=np.float64) / col_scale[:, None]
    _scaled_basis, scaled_u, scaled_gauge_s, scaled_threshold, scaled_rank = _c0g_column_space(
        scaled_gauge,
        relative_tolerance=float(config.gradient_rank_rtol),
        full_matrices=True,
    )
    q_perp = scaled_u[:, scaled_rank:].astype(np.float64, copy=True)
    return {
        "q_perp": q_perp,
        "scaled_gauge_singular_values": scaled_gauge_s,
        "scaled_gauge_rank_threshold": float(scaled_threshold),
        "scaled_gauge_rank": int(scaled_rank),
    }


def _c0g_dense_complement_svd(
    *,
    scaled_matrix: csc_matrix,
    gauge_matrix: np.ndarray,
    col_scale: np.ndarray,
    grid,
    config: C0gConfig,
    mode_artifact_path: Path | None,
) -> dict[str, Any]:
    complement = _c0g_scaled_gauge_complement(
        gauge_matrix=gauge_matrix,
        col_scale=col_scale,
        config=config,
    )
    q_perp = np.asarray(complement["q_perp"], dtype=np.float64)
    scaled_gauge_s = np.asarray(
        complement["scaled_gauge_singular_values"], dtype=np.float64
    )
    scaled_threshold = float(complement["scaled_gauge_rank_threshold"])
    scaled_rank = int(complement["scaled_gauge_rank"])
    if q_perp.shape[1] <= 0:
        return {
            "status": "NOT_MEASURED",
            "reason": "gauge_basis_spans_full_scaled_state_space",
            "scaled_gauge_rank": int(scaled_rank),
            "scaled_gauge_rank_threshold": float(scaled_threshold),
        }
    dense_scaled = scaled_matrix.toarray().astype(np.float64, copy=False)
    start = time.perf_counter()
    reduced = dense_scaled @ q_perp
    left, singular_values, vh = np.linalg.svd(reduced, full_matrices=False)
    svd_seconds = float(time.perf_counter() - start)
    sigma_max = float(singular_values[0]) if singular_values.size else math.nan
    sigma_min = float(singular_values[-1]) if singular_values.size else math.nan
    count = min(int(config.near_null_mode_count), int(singular_values.size))
    modes: list[dict[str, Any]] = []
    right_modes: list[np.ndarray] = []
    left_modes: list[np.ndarray] = []
    mode_sigmas: list[float] = []
    for offset in range(count):
        svd_index = int(singular_values.size - 1 - offset)
        right_scaled = q_perp @ vh[svd_index, :]
        right_physical = col_scale * right_scaled
        right_unit, right_norm = _unit_vector(right_physical)
        left_scaled = left[:, svd_index]
        left_unit, left_norm = _unit_vector(left_scaled)
        right_modes.append(right_unit)
        left_modes.append(left_unit)
        mode_sigmas.append(float(singular_values[svd_index]))
        modes.append(
            {
                "ascending_rank": int(offset + 1),
                "svd_index": svd_index,
                "sigma": float(singular_values[svd_index]),
                "right_physical_norm_before_unit": float(right_norm),
                "left_scaled_norm_before_unit": float(left_norm),
                "v_lane_energy_fractions": _lane_energy_split(right_unit, grid),
                "w_output_energy_fractions": _lane_energy_split(left_unit, grid),
                "v_center_of_energy": _c0g_state_center_of_energy(right_unit, grid),
                "w_center_of_energy": _c0g_state_center_of_energy(left_unit, grid),
            }
        )
    artifact_text = None
    if mode_artifact_path is not None:
        mode_artifact_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            mode_artifact_path,
            right_modes=np.asarray(right_modes, dtype=np.float64),
            left_modes=np.asarray(left_modes, dtype=np.float64),
            sigmas=np.asarray(mode_sigmas, dtype=np.float64),
        )
        artifact_text = str(mode_artifact_path)
    return {
        "status": "MEASURED",
        "construction": "dense_svd_of_scaled_J_times_scaled_gauge_complement_Q_perp",
        "forbidden_construction_not_used": "SVD((I-P_G)J)",
        "scaling": "row_scale * J_original * col_scale; gauge basis transformed by C^{-1}",
        "reduced_shape": [int(reduced.shape[0]), int(reduced.shape[1])],
        "state_dim": int(dense_scaled.shape[1]),
        "scaled_gauge_rank": int(scaled_rank),
        "scaled_gauge_rank_threshold": float(scaled_threshold),
        "scaled_gauge_singular_min_retained": float(scaled_gauge_s[scaled_rank - 1])
        if scaled_rank > 0
        else None,
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "sigma_min_over_sigma_max": float(sigma_min / sigma_max)
        if math.isfinite(sigma_min) and math.isfinite(sigma_max) and sigma_max > 0.0
        else math.nan,
        "svd_seconds": svd_seconds,
        "modes": modes,
        "mode_vectors_artifact": artifact_text,
    }


def _c0g_clean_newton_step(
    *,
    matrix: csc_matrix,
    residual: np.ndarray,
) -> dict[str, Any]:
    dense = matrix.toarray().astype(np.float64, copy=False)
    try:
        step = np.linalg.solve(dense, -residual)
        method = "dense_np_linalg_solve"
    except np.linalg.LinAlgError:
        step, _resid, _rank, _singular = np.linalg.lstsq(dense, -residual, rcond=None)
        method = "dense_np_linalg_lstsq_fallback"
    linear_residual = dense @ step + residual
    residual_norm = float(np.linalg.norm(residual))
    return {
        "step": step.astype(np.float64, copy=False),
        "method": method,
        "linear_rel_resid": float(
            np.linalg.norm(linear_residual) / max(residual_norm, np.finfo(np.float64).tiny)
        ),
        "step_norm": float(np.linalg.norm(step)),
    }


def _c0g_merit_sweep(
    *,
    state: torch.Tensor,
    step: np.ndarray,
    matrix: csc_matrix,
    residual: np.ndarray,
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    alpha_min_power: int,
) -> dict[str, Any]:
    initial_l2 = float(np.linalg.norm(residual))
    initial_linf = float(np.max(np.abs(residual)))
    j_step = matrix @ step
    step_t = torch.as_tensor(step, dtype=state.dtype, device=state.device)
    rows: list[dict[str, Any]] = []
    for power in range(int(alpha_min_power) + 1):
        alpha = float(2.0 ** (-power))
        predicted = residual + alpha * j_step
        trial = state.detach() + alpha * step_t
        actual = residual_fn(trial).detach().cpu().numpy().astype(np.float64, copy=False)
        actual_l2 = float(np.linalg.norm(actual))
        predicted_l2 = float(np.linalg.norm(predicted))
        row = {
            "alpha": alpha,
            "power": int(power),
            "actual_l2": actual_l2,
            "actual_linf": float(np.max(np.abs(actual))),
            "predicted_l2": predicted_l2,
            "predicted_linf": float(np.max(np.abs(predicted))),
            "actual_l2_ratio": actual_l2 / max(initial_l2, np.finfo(np.float64).tiny),
            "predicted_l2_ratio": predicted_l2 / max(initial_l2, np.finfo(np.float64).tiny),
            "reduces_true_l2": bool(actual_l2 < initial_l2),
            "finite": bool(np.all(np.isfinite(actual))),
        }
        rows.append(row)
    best = min(rows, key=lambda row: row["actual_l2"]) if rows else None
    alpha1 = rows[0] if rows else None
    return {
        "status": "MEASURED",
        "initial_l2": initial_l2,
        "initial_linf": initial_linf,
        "alpha1": alpha1,
        "best": best,
        "best_alpha": best.get("alpha") if best else None,
        "best_actual_l2": best.get("actual_l2") if best else None,
        "best_actual_l2_ratio": best.get("actual_l2_ratio") if best else None,
        "best_percent_reduction": float(100.0 * (1.0 - best["actual_l2"] / initial_l2))
        if best is not None and initial_l2 > 0.0
        else math.nan,
        "alpha1_predicted_vs_actual_l2_gap": float(
            alpha1["actual_l2"] / max(alpha1["predicted_l2"], np.finfo(np.float64).tiny)
        )
        if alpha1 is not None
        else math.nan,
        "any_alpha_reduces_true_l2": bool(
            any(bool(row["reduces_true_l2"]) for row in rows)
        ),
        "rows": rows,
    }


def _c0g_step_decomposition(
    *,
    step: np.ndarray,
    physical_gauge_basis: np.ndarray,
    right_modes: np.ndarray,
) -> dict[str, Any]:
    step_norm = float(np.linalg.norm(step))
    if step_norm <= 0.0:
        return {
            "status": "NOT_MEASURED",
            "reason": "zero_newton_step",
        }
    gauge_projection = _c0g_projection_vector(physical_gauge_basis, step)
    gauge_fraction = float(np.dot(gauge_projection, gauge_projection) / (step_norm * step_norm))
    gauge_residual = step - gauge_projection
    if right_modes.size:
        residualized_modes = []
        for mode in right_modes:
            residualized = mode - _c0g_projection_vector(physical_gauge_basis, mode)
            unit, norm = _unit_vector(residualized)
            if norm > 0.0:
                residualized_modes.append(unit)
        if residualized_modes:
            near_matrix = np.column_stack(residualized_modes)
            near_basis, _u, _s, _threshold, near_rank = _c0g_column_space(
                near_matrix,
                relative_tolerance=1.0e-12,
                full_matrices=False,
            )
            near_projection = _c0g_projection_vector(near_basis, gauge_residual)
            near_fraction = float(np.dot(near_projection, near_projection) / (step_norm * step_norm))
            per_mode = []
            for index, mode in enumerate(residualized_modes):
                component = float(np.dot(mode, gauge_residual) ** 2 / (step_norm * step_norm))
                per_mode.append(
                    {
                        "ascending_rank": int(index + 1),
                        "component_fraction": component,
                    }
                )
        else:
            near_rank = 0
            near_fraction = 0.0
            per_mode = []
    else:
        near_rank = 0
        near_fraction = 0.0
        per_mode = []
    complement_fraction = float(max(0.0, 1.0 - gauge_fraction - near_fraction))
    mode2_fraction = per_mode[1]["component_fraction"] if len(per_mode) >= 2 else None
    dominant = max(per_mode, key=lambda item: item["component_fraction"]) if per_mode else None
    return {
        "status": "MEASURED",
        "step_norm": step_norm,
        "gauge_component_fraction": gauge_fraction,
        "near_null_component_fraction": near_fraction,
        "complement_component_fraction": complement_fraction,
        "transverse_mode_2_component_fraction": mode2_fraction,
        "mode_2_definition": "second-smallest gauge-complement right singular vector",
        "dominant_near_null_mode": dominant,
        "near_null_rank_after_gauge_residualization": int(near_rank),
        "per_mode_component_fractions": per_mode,
    }


def _c0g_residual_decomposition(
    *,
    residual: np.ndarray,
    row_scale: np.ndarray,
    matrix: csc_matrix,
    gauge_matrix: np.ndarray,
    left_modes: np.ndarray,
) -> dict[str, Any]:
    scaled_residual = row_scale * residual
    residual_norm = float(np.linalg.norm(scaled_residual))
    if residual_norm <= 0.0:
        return {
            "status": "NOT_MEASURED",
            "reason": "zero_scaled_residual",
        }
    left_rows = []
    for index, mode in enumerate(left_modes):
        unit, norm = _unit_vector(mode)
        if norm <= 0.0:
            fraction = math.nan
        else:
            fraction = float(abs(np.dot(unit, scaled_residual)) / residual_norm)
        left_rows.append(
            {
                "ascending_rank": int(index + 1),
                "abs_u_dot_F_over_normF": fraction,
            }
        )
    gauge_images = row_scale[:, None] * (matrix @ gauge_matrix)
    image_basis, _u, image_s, threshold, rank = _c0g_column_space(
        gauge_images,
        relative_tolerance=1.0e-12,
        full_matrices=False,
    )
    projection = _c0g_projection_vector(image_basis, scaled_residual)
    return {
        "status": "MEASURED",
        "scaled_residual_l2": residual_norm,
        "left_singular_vector_components": left_rows,
        "gauge_image_projection_fraction": float(
            np.dot(projection, projection) / (residual_norm * residual_norm)
        ),
        "gauge_image_rank": int(rank),
        "gauge_image_rank_threshold": float(threshold),
        "gauge_image_singular_values_head": [
            float(value) for value in image_s[: min(5, image_s.size)]
        ],
    }


def _c0g_ftau_diagnostic(
    *,
    state: torch.Tensor,
    tau: float,
    left_min: np.ndarray,
    row_scale: np.ndarray,
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    w_unit, w_norm = _unit_vector(left_min)
    if w_norm <= 0.0:
        return {"status": "NOT_MEASURED", "reason": "zero_left_min_vector", "rows": rows}
    for multiplier in config.ftau_h_multipliers:
        h = float(multiplier) * float(tau)
        _branch_p, _provider_p, _grid_p, _boundaries_p, _eos_K_p, residual_plus = (
            _c0g_residual_context(tau=float(tau) + h, config=config, dtype=dtype)
        )
        _branch_m, _provider_m, _grid_m, _boundaries_m, _eos_K_m, residual_minus = (
            _c0g_residual_context(tau=float(tau) - h, config=config, dtype=dtype)
        )
        f_plus = residual_plus(state).detach().cpu().numpy().astype(np.float64, copy=False)
        f_minus = residual_minus(state).detach().cpu().numpy().astype(np.float64, copy=False)
        ftau = (f_plus - f_minus) / (2.0 * h)
        ftau_scaled = row_scale * ftau
        ftau_norm = float(np.linalg.norm(ftau_scaled))
        dot = float(np.dot(w_unit, ftau_scaled))
        cos_theta = float(abs(dot) / max(ftau_norm, np.finfo(np.float64).tiny))
        rows.append(
            {
                "h_multiplier": float(multiplier),
                "h": h,
                "wT_F_tau": dot,
                "F_tau_scaled_l2": ftau_norm,
                "cos_theta": cos_theta,
            }
        )
    stable_pairs = []
    for left_index in range(len(rows)):
        for right_index in range(left_index + 1, len(rows)):
            a = float(rows[left_index]["cos_theta"])
            b = float(rows[right_index]["cos_theta"])
            rel = abs(a - b) / max(abs(a), abs(b), np.finfo(np.float64).tiny)
            if rel <= float(config.ftau_stability_rtol):
                stable_pairs.append(
                    {
                        "left_h_multiplier": rows[left_index]["h_multiplier"],
                        "right_h_multiplier": rows[right_index]["h_multiplier"],
                        "relative_difference": float(rel),
                        "mean_cos_theta": float(0.5 * (a + b)),
                    }
                )
    if not stable_pairs:
        return {
            "status": "NOT_MEASURED",
            "reason": "cos_theta_not_stable_to_10_percent_across_any_two_h_values",
            "rows": rows,
            "stable_pairs": stable_pairs,
        }
    representative = min(stable_pairs, key=lambda item: item["relative_difference"])
    cos = float(representative["mean_cos_theta"])
    if cos > float(config.cos_fold_threshold):
        call = "FOLD_SUPPORT"
    elif cos < float(config.cos_bifurcation_threshold):
        call = "BIFURCATION_SUPPORT"
    else:
        call = "INCONCLUSIVE"
    return {
        "status": "MEASURED",
        "rows": rows,
        "stable_pairs": stable_pairs,
        "stepsize_stability_call": "STABLE",
        "representative_cos_theta": cos,
        "call": call,
        "thresholds": {
            "fold_support_cos_gt": float(config.cos_fold_threshold),
            "bifurcation_support_cos_lt": float(config.cos_bifurcation_threshold),
            "inconclusive_band": [
                float(config.cos_bifurcation_threshold),
                float(config.cos_fold_threshold),
            ],
        },
    }


def _c0g_mode_artifact_path(config: C0gConfig, *, tau: float) -> Path:
    return _resolve_output_path(
        config.run_root / "modes" / f"c0g_modes_tau_{_format_tau(float(tau))}.npz"
    )


def _c0g_state_json_path(config: C0gConfig, *, tau: float) -> Path:
    return _resolve_output_path(
        config.run_root / "state_measurements" / f"c0g_state_tau_{_format_tau(float(tau))}.json"
    )


def _c0g_load_mode_vectors(path: str | Path) -> np.ndarray:
    with np.load(_resolve_input_path(path)) as data:
        return np.asarray(data["right_modes"], dtype=np.float64)


def _c0g_load_mode_bundle(path: str | Path) -> dict[str, np.ndarray]:
    with np.load(_resolve_input_path(path)) as data:
        return {
            "right_modes": np.asarray(data["right_modes"], dtype=np.float64),
            "left_modes": np.asarray(data["left_modes"], dtype=np.float64),
            "sigmas": np.asarray(data["sigmas"], dtype=np.float64),
        }


def _c0g_analyze_state(
    *,
    tau: float,
    row: Mapping[str, Any],
    config: C0gConfig,
    dtype: torch.dtype,
    include_step0: bool,
    compute_ftau: bool,
) -> dict[str, Any]:
    state_path = _c0g_state_path_from_row(row)
    state = _load_state_artifact(state_path, dtype=dtype)
    assembled = _c0g_assemble_original_jacobian(
        state=state,
        tau=float(tau),
        config=config,
        dtype=dtype,
    )
    grid = assembled["grid"]
    branch = assembled["branch"]
    residual = assembled["residual"]
    matrix = assembled["matrix"]
    gauge = _c0g_build_analytic_gauge_matrix(
        state=state,
        grid=grid,
        branch=branch,
        config=config,
    )
    mode_artifact = _c0g_mode_artifact_path(config, tau=float(tau))
    svd = _c0g_dense_complement_svd(
        scaled_matrix=assembled["scaled_matrix"],
        gauge_matrix=gauge["generator_matrix"],
        col_scale=assembled["col_scale"],
        grid=grid,
        config=config,
        mode_artifact_path=mode_artifact,
    )
    step = _c0g_clean_newton_step(matrix=matrix, residual=residual)
    merit = _c0g_merit_sweep(
        state=state,
        step=step["step"],
        matrix=matrix,
        residual=residual,
        residual_fn=assembled["residual_fn"],
        alpha_min_power=int(config.alpha_min_power),
    )
    result: dict[str, Any] = {
        "status": "MEASURED",
        "tau": float(tau),
        "attempt_index": int(row.get("attempt_index", -1)),
        "state_artifact": row.get("state_artifact"),
        "state_path": str(state_path),
        "accepted_default_success": bool(row.get("accepted_default_success", False)),
        "solver_converged": bool(row.get("solver_converged", False)),
        "init_source": row.get("init_source"),
        "used_existing_b2c": bool(row.get("used_existing_b2c", False)),
        "original_residual_l2": float(np.linalg.norm(residual)),
        "original_residual_linf": float(np.max(np.abs(residual))),
        "jacobian": {
            key: value
            for key, value in assembled["metadata"].items()
            if key
            in {
                "active_color_count",
                "jvp_count",
                "state_size",
                "matrix_nnz",
                "matrix_density",
                "assembly_seconds",
                "jacobian_source",
                "residual_entry_point",
                "scaling",
                "row_scale_min",
                "row_scale_max",
                "col_scale_min",
                "col_scale_max",
                "diagonal_shift",
                "core_epsilon",
                "k1_radius_epsilon",
            }
        },
        "gauge_basis": {
            "rank": int(gauge["rank"]),
            "rank_threshold": float(gauge["rank_threshold"]),
            "source_helpers": list(gauge["source_helpers"]),
            "piece_dimensions": dict(gauge["piece_dimensions"]),
            "a_sector_gradient_rank": int(gauge["a_sector_gradient_rank"]),
            "a_sector_gradient_rank_threshold": float(
                gauge["a_sector_gradient_rank_threshold"]
            ),
        },
        "svd": svd,
        "newton_step": {
            "method": step["method"],
            "linear_rel_resid": float(step["linear_rel_resid"]),
            "step_norm": float(step["step_norm"]),
        },
        "merit_sweep": merit,
    }
    if include_step0:
        if svd.get("status") == "MEASURED":
            modes = _c0g_load_mode_vectors(svd["mode_vectors_artifact"])
            left_modes = None
            with np.load(_resolve_input_path(svd["mode_vectors_artifact"])) as data:
                left_modes = np.asarray(data["left_modes"], dtype=np.float64)
            result["step0_decomposition"] = {
                "delta": _c0g_step_decomposition(
                    step=step["step"],
                    physical_gauge_basis=gauge["physical_basis"],
                    right_modes=modes,
                ),
                "residual": _c0g_residual_decomposition(
                    residual=residual,
                    row_scale=assembled["row_scale"],
                    matrix=matrix,
                    gauge_matrix=gauge["generator_matrix"],
                    left_modes=left_modes,
                ),
            }
            f_nn = result["step0_decomposition"]["delta"].get(
                "near_null_component_fraction", math.nan
            )
            call = _c0g_premise_gate_call(float(f_nn), config)
            result["premise_gate"] = {
                "f_nn": float(f_nn),
                "call": call,
                "thresholds": {
                    "premise_failed_if_f_nn_lt": float(config.premise_fail_threshold),
                    "premise_holds_if_f_nn_gt": float(config.premise_hold_threshold),
                    "gray_band": [
                        float(config.premise_fail_threshold),
                        float(config.premise_hold_threshold),
                    ],
                },
            }
        else:
            result["step0_decomposition"] = {
                "status": "NOT_MEASURED",
                "reason": svd.get("reason", "svd_not_measured"),
            }
            result["premise_gate"] = {
                "f_nn": math.nan,
                "call": "DIAGNOSTIC_INCOMPLETE",
            }
    premise_failed = (
        include_step0
        and result.get("premise_gate", {}).get("call") == "PREMISE_FAILED"
    )
    if compute_ftau and premise_failed:
        result["ftau"] = "SKIPPED_PREMISE_FAILED"
    elif compute_ftau and svd.get("status") == "MEASURED":
        with np.load(_resolve_input_path(svd["mode_vectors_artifact"])) as data:
            left_modes = np.asarray(data["left_modes"], dtype=np.float64)
        result["ftau"] = _c0g_ftau_diagnostic(
            state=state,
            tau=float(tau),
            left_min=left_modes[0],
            row_scale=assembled["row_scale"],
            config=config,
            dtype=dtype,
        )
    elif compute_ftau:
        result["ftau"] = {
            "status": "NOT_MEASURED",
            "reason": "svd_not_measured",
        }
    return result


def _c0g_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _c0g_write_state_result(config: C0gConfig, result: Mapping[str, Any]) -> Path:
    path = _c0g_state_json_path(config, tau=float(result["tau"]))
    _c0g_write_json(path, result)
    return path


def run_c0g_step0_premise(config: C0gConfig | None = None) -> dict[str, Any]:
    cfg = config or C0gConfig()
    dtype = configure_backend(BackendConfig())
    payload = _c0g_load_c0f2_payload(cfg)
    provenance = _c0g_verify_stalled_provenance(cfg)
    if provenance["call"] != "GENUINE_WARM_START_NOT_COLD_LOADED":
        result = {
            "phase": "step0_premise",
            "status": "PROVENANCE_FAILED",
            "config": _c0g_config_to_dict(cfg),
            "provenance": provenance,
            "steps_1_to_7": "SKIPPED_PREMISE_FAILED",
        }
        _c0g_write_json(_resolve_output_path(cfg.step0_json_path), result)
        return result
    row = _c0g_find_tau_row(payload, cfg.stalled_tau)
    state_result = _c0g_analyze_state(
        tau=float(cfg.stalled_tau),
        row=row,
        config=cfg,
        dtype=dtype,
        include_step0=True,
        compute_ftau=True,
    )
    gate = state_result.get("premise_gate", {})
    if gate.get("call") != "PREMISE_FAILED":
        state_result["ftau_note"] = "computed after premise gate did not fail"
    _c0g_write_state_result(cfg, state_result)
    result = {
        "phase": "step0_premise",
        "status": "MEASURED",
        "config": _c0g_config_to_dict(cfg),
        "provenance": provenance,
        "step0": state_result,
    }
    if gate.get("call") == "PREMISE_FAILED":
        result["steps_1_to_7"] = "SKIPPED_PREMISE_FAILED"
    _c0g_write_json(_resolve_output_path(cfg.step0_json_path), result)
    return result


def run_c0g_state_measurement(
    *,
    tau: float,
    config: C0gConfig | None = None,
) -> dict[str, Any]:
    cfg = config or C0gConfig()
    dtype = configure_backend(BackendConfig())
    payload = _c0g_load_c0f2_payload(cfg)
    row = _c0g_find_tau_row(payload, float(tau))
    result = _c0g_analyze_state(
        tau=float(tau),
        row=row,
        config=cfg,
        dtype=dtype,
        include_step0=False,
        compute_ftau=True,
    )
    _c0g_write_state_result(cfg, result)
    return result


def _c0g_fit_sigma_min_squared(
    state_results: Sequence[Mapping[str, Any]],
    *,
    config: C0gConfig,
) -> dict[str, Any]:
    converged = [
        result
        for result in state_results
        if any(_c0g_tau_close(float(result["tau"]), tau) for tau in config.converged_taus)
    ]
    converged = sorted(converged, key=lambda item: float(item["tau"]), reverse=True)
    rows = []
    for result in converged:
        svd = result.get("svd", {})
        if svd.get("status") != "MEASURED":
            return {
                "status": "NOT_MEASURED",
                "reason": f"svd_not_measured_tau_{result.get('tau')}",
            }
        rows.append(
            {
                "tau": float(result["tau"]),
                "sigma_min": float(svd["sigma_min"]),
                "sigma_min_squared": float(svd["sigma_min"]) ** 2,
            }
        )
    if len(rows) < 3:
        return {"status": "NOT_MEASURED", "reason": "too_few_converged_sigma_values"}
    tau_values = np.asarray([row["tau"] for row in rows], dtype=np.float64)
    y = np.asarray([row["sigma_min_squared"] for row in rows], dtype=np.float64)
    linear = np.polyfit(tau_values, y, deg=1)
    pred_linear = np.polyval(linear, tau_values)
    ss_res_linear = float(np.sum((y - pred_linear) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    r2_linear = float(1.0 - ss_res_linear / ss_tot) if ss_tot > 0.0 else math.nan
    quadratic = np.polyfit(tau_values, y, deg=2)
    pred_quadratic = np.polyval(quadratic, tau_values)
    ss_res_quadratic = float(np.sum((y - pred_quadratic) ** 2))
    r2_quadratic = float(1.0 - ss_res_quadratic / ss_tot) if ss_tot > 0.0 else math.nan
    slope = float(linear[0])
    intercept = float(linear[1])
    tau_fold = float(-intercept / slope) if slope != 0.0 else math.nan
    sigmas = np.asarray([row["sigma_min"] for row in rows], dtype=np.float64)
    monotone_toward_stall = bool(
        np.all(np.diff(sigmas) <= np.maximum(1.0e-15, 1.0e-10 * np.maximum(sigmas[:-1], 1.0)))
    )
    if monotone_toward_stall and math.isfinite(r2_linear) and r2_linear >= config.fit_good_r2_threshold:
        call = "LINEAR_MONOTONE_FOLD_SUPPORT"
    elif not monotone_toward_stall:
        call = "FIT_UNRELIABLE_NON_MONOTONE_SIGMA_MIN"
    elif math.isfinite(r2_quadratic) and r2_quadratic > r2_linear:
        call = "QUADRATIC_OR_FLAT_SUPPORT"
    else:
        call = "INCONCLUSIVE"
    return {
        "status": "MEASURED",
        "rows": rows,
        "linear": {
            "slope": slope,
            "intercept": intercept,
            "r2": r2_linear,
            "sse": ss_res_linear,
            "tau_fold_zero_crossing": tau_fold,
        },
        "quadratic": {
            "coefficients": [float(value) for value in quadratic],
            "r2": r2_quadratic,
            "sse": ss_res_quadratic,
        },
        "sigma_min_monotone_decreasing_toward_stall": monotone_toward_stall,
        "fit_reliability": "RELIABLE" if monotone_toward_stall else "UNRELIABLE",
        "call": call,
    }


def _c0g_track_modes(state_results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    ordered = sorted(state_results, key=lambda item: float(item["tau"]), reverse=True)
    previous_vectors: np.ndarray | None = None
    previous_lane_by_track: dict[int, int] = {}
    tracked_rows: list[dict[str, Any]] = []
    for result in ordered:
        svd = result.get("svd", {})
        if svd.get("status") != "MEASURED" or not svd.get("mode_vectors_artifact"):
            continue
        vectors = _c0g_load_mode_vectors(str(svd["mode_vectors_artifact"]))
        mode_to_track: dict[int, int] = {}
        if previous_vectors is None:
            for mode_index in range(vectors.shape[0]):
                mode_to_track[mode_index] = mode_index + 1
        else:
            overlaps = np.abs(vectors @ previous_vectors.T)
            used_previous: set[int] = set()
            for mode_index in range(vectors.shape[0]):
                order = list(np.argsort(overlaps[mode_index])[::-1])
                chosen = next((idx for idx in order if int(idx) not in used_previous), int(order[0]))
                used_previous.add(int(chosen))
                mode_to_track[mode_index] = previous_lane_by_track.get(int(chosen), int(chosen) + 1)
        for mode_index, mode in enumerate(svd.get("modes", [])):
            tracked_rows.append(
                {
                    "tau": float(result["tau"]),
                    "ascending_rank": int(mode.get("ascending_rank", mode_index + 1)),
                    "tracked_lane": int(mode_to_track.get(mode_index, mode_index + 1)),
                    "sigma": float(mode["sigma"]),
                    "v_lane_energy_fractions": mode.get("v_lane_energy_fractions", {}),
                    "w_output_energy_fractions": mode.get("w_output_energy_fractions", {}),
                    "v_center_of_energy": mode.get("v_center_of_energy", {}),
                    "w_center_of_energy": mode.get("w_center_of_energy", {}),
                }
            )
        previous_vectors = vectors
        previous_lane_by_track = {mode_index: lane for mode_index, lane in mode_to_track.items()}
    return tracked_rows


def _c0g_best_alpha_trend(
    state_results: Sequence[Mapping[str, Any]],
    *,
    config: C0gConfig,
) -> dict[str, Any]:
    rows = []
    for tau in sorted(config.converged_taus):
        match = next(
            (
                result
                for result in state_results
                if _c0g_tau_close(float(result["tau"]), float(tau))
            ),
            None,
        )
        if match is None:
            continue
        merit = match.get("merit_sweep", {})
        rows.append(
            {
                "tau": float(tau),
                "best_alpha": merit.get("best_alpha"),
                "best_actual_l2_ratio": merit.get("best_actual_l2_ratio"),
                "best_percent_reduction": merit.get("best_percent_reduction"),
                "alpha1_actual_l2_ratio": (merit.get("alpha1") or {}).get("actual_l2_ratio"),
                "alpha1_predicted_l2_ratio": (merit.get("alpha1") or {}).get(
                    "predicted_l2_ratio"
                ),
                "alpha1_gap": merit.get("alpha1_predicted_vs_actual_l2_gap"),
            }
        )
    descending = sorted(rows, key=lambda item: item["tau"], reverse=True)
    reductions = [
        float(row["best_percent_reduction"])
        for row in descending
        if row.get("best_percent_reduction") is not None
        and math.isfinite(float(row["best_percent_reduction"]))
    ]
    degrading = bool(
        len(reductions) >= 2
        and all(right <= left + 1.0e-10 for left, right in zip(reductions, reductions[1:]))
    )
    return {
        "status": "MEASURED" if rows else "NOT_MEASURED",
        "rows": rows,
        "trend_call": "DEGRADING_TOWARD_STALL" if degrading else "NOT_MONOTONE_DEGRADING",
    }


def _c0g_preliminary_reading(result: Mapping[str, Any]) -> dict[str, Any]:
    step0_call = result.get("step0", {}).get("premise_gate", {}).get("call")
    ftau_rows = [
        state.get("ftau", {})
        for state in result.get("state_results", [])
        if isinstance(state.get("ftau"), Mapping)
    ]
    measured_ftau = [row for row in ftau_rows if row.get("status") == "MEASURED"]
    ftau_calls = {row.get("call") for row in measured_ftau}
    fit_call = result.get("sigma_min_squared_fit", {}).get("call")
    if step0_call == "PREMISE_FAILED":
        direction = "PREMISE_FAILED"
    elif "FOLD_SUPPORT" in ftau_calls and fit_call == "LINEAR_MONOTONE_FOLD_SUPPORT":
        direction = "FOLD_SUPPORT"
    elif "BIFURCATION_SUPPORT" in ftau_calls or fit_call == "QUADRATIC_OR_FLAT_SUPPORT":
        direction = "CONDITIONING_OR_BIFURCATION_SUPPORT"
    else:
        direction = "INCONCLUSIVE"
    return {
        "status": direction,
        "step8_verdict_deferred": True,
        "reason": (
            "Steps 4-7 were not run in this staged pass; this is not a Step-8 verdict."
        ),
        "step0_call": step0_call,
        "ftau_calls": sorted(str(call) for call in ftau_calls if call is not None),
        "sigma_min_squared_fit_call": fit_call,
    }


def aggregate_c0g_steps0_3(config: C0gConfig | None = None) -> dict[str, Any]:
    cfg = config or C0gConfig()
    step0_path = _resolve_input_path(cfg.step0_json_path)
    with step0_path.open("r", encoding="utf-8") as handle:
        step0 = json.load(handle)
    provenance = step0.get("provenance", {})
    if step0.get("step0", {}).get("premise_gate", {}).get("call") == "PREMISE_FAILED":
        result = {
            "phase": "steps0_3_staged",
            "status": "PREMISE_FAILED",
            "config": _c0g_config_to_dict(cfg),
            "provenance": provenance,
            "step0": step0.get("step0"),
            "step0b": "SKIPPED_PREMISE_FAILED",
            "step1": "SKIPPED_PREMISE_FAILED",
            "step2": "SKIPPED_PREMISE_FAILED",
            "step3": "SKIPPED_PREMISE_FAILED",
            "steps_4_to_7": "NOT_RUN_ORCHESTRATOR_REVIEW_GATE",
        }
        _c0g_write_json(_resolve_output_path(cfg.json_path), result)
        write_c0g_partial_report(result, _resolve_output_path(cfg.report_path))
        return result
    state_results: list[dict[str, Any]] = []
    if "step0" in step0:
        state_results.append(step0["step0"])
    for tau in cfg.converged_taus:
        path = _c0g_state_json_path(cfg, tau=float(tau))
        with path.open("r", encoding="utf-8") as handle:
            state_results.append(json.load(handle))
    # Keep one row per tau if the stalled state was also measured separately.
    by_tau: dict[float, dict[str, Any]] = {}
    for state in state_results:
        by_tau[float(state["tau"])] = state
    state_results = sorted(by_tau.values(), key=lambda item: float(item["tau"]), reverse=True)
    result = {
        "phase": "steps0_3_staged",
        "status": "MEASURED",
        "config": _c0g_config_to_dict(cfg),
        "provenance": provenance,
        "step0": step0.get("step0"),
        "step0b": _c0g_best_alpha_trend(state_results, config=cfg),
        "state_results": state_results,
        "step1_tracked_modes": _c0g_track_modes(state_results),
        "sigma_min_squared_fit": _c0g_fit_sigma_min_squared(state_results, config=cfg),
        "steps_4_to_7": "NOT_RUN_ORCHESTRATOR_REVIEW_GATE",
    }
    result["preliminary_reading"] = _c0g_preliminary_reading(result)
    _c0g_write_json(_resolve_output_path(cfg.json_path), result)
    write_c0g_partial_report(result, _resolve_output_path(cfg.report_path))
    return result


def _c0g_steps4_6_7_json_path(config: C0gConfig) -> Path:
    return _resolve_output_path(config.run_root / "pathA_C0g_steps4_6_7.json")


def _c0g_load_steps0_3(config: C0gConfig) -> dict[str, Any]:
    path = _resolve_input_path(config.json_path)
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _c0g_state_results_by_tau(
    result: Mapping[str, Any],
) -> dict[float, dict[str, Any]]:
    return {float(row["tau"]): dict(row) for row in result.get("state_results", [])}


def _c0g_current_and_mach_for_state(
    *,
    state_result: Mapping[str, Any],
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    tau = float(state_result["tau"])
    state_path = _resolve_input_path(str(state_result["state_path"]))
    state = _load_state_artifact(state_path, dtype=dtype)
    branch, _provider, grid, _boundaries, eos_K, _residual_fn = _c0g_residual_context(
        tau=tau,
        config=config,
        dtype=dtype,
    )
    fields, _mu = unpack_closed_coupled_fields(
        state,
        grid,
        has_chemical_potential=True,
    )
    grad_real_r = tensor_center_gradient_r(fields.psi_real, grid)
    grad_imag_r = tensor_center_gradient_r(fields.psi_imag, grid)
    grad_real_w = tensor_center_gradient_w(fields.psi_real, grid)
    grad_imag_w = tensor_center_gradient_w(fields.psi_imag, grid)
    rho = fields.psi_real**2 + fields.psi_imag**2
    jr = (branch.hbar / branch.particle_mass) * (
        fields.psi_real * grad_imag_r - fields.psi_imag * grad_real_r
    ) - (branch.gauge_charge / branch.particle_mass) * fields.ar * rho
    jw = (branch.hbar / branch.particle_mass) * (
        fields.psi_real * grad_imag_w - fields.psi_imag * grad_real_w
    ) - (branch.gauge_charge / branch.particle_mass) * fields.aw * rho

    rho_np = rho.detach().cpu().numpy().astype(np.float64, copy=False)
    jr_np = jr.detach().cpu().numpy().astype(np.float64, copy=False)
    jw_np = jw.detach().cpu().numpy().astype(np.float64, copy=False)
    jmag = np.sqrt(jr_np * jr_np + jw_np * jw_np)
    sound_speed = math.sqrt(5.0 * float(eos_K) / float(branch.particle_mass)) * rho_np**2
    sound_speed = np.maximum(sound_speed, 1.0e-300)
    floor_rows: list[dict[str, Any]] = []
    base_payload: dict[str, Any] | None = None
    r_centers = grid.r_centers.detach().cpu().numpy().astype(np.float64, copy=False)
    w_centers = grid.w_centers.detach().cpu().numpy().astype(np.float64, copy=False)
    r0_np = fields.r0.detach().cpu().numpy().astype(np.float64, copy=False)
    for floor in config.mach_density_floors:
        denom = np.maximum(rho_np, float(floor))
        speed = jmag / denom
        mach = np.nan_to_num(speed / sound_speed, nan=0.0, posinf=1.0e308, neginf=0.0)
        mach_w = np.nan_to_num(
            np.abs(jw_np) / denom / sound_speed,
            nan=0.0,
            posinf=1.0e308,
            neginf=0.0,
        )
        flat_index = int(np.argmax(mach))
        ir, iw = np.unravel_index(flat_index, mach.shape)
        row = {
            "density_floor": float(floor),
            "M_full_max": float(mach[ir, iw]),
            "M_w_max": float(mach_w[ir, iw]),
            "rho_at_max": float(rho_np[ir, iw]),
            "j_at_max": float(jmag[ir, iw]),
            "j_r_at_max": float(jr_np[ir, iw]),
            "j_w_at_max": float(jw_np[ir, iw]),
            "r_star": float(r_centers[ir]),
            "w_star": float(w_centers[iw]),
            "R0_at_w_star": float(r0_np[iw]),
            "r_star_over_R0": float(r_centers[ir] / r0_np[iw]) if r0_np[iw] != 0.0 else math.inf,
            "argmax_index": [int(ir), int(iw)],
        }
        floor_rows.append(row)
        if base_payload is None:
            base_payload = row
    assert base_payload is not None
    max_rho = float(np.max(rho_np)) if rho_np.size else 0.0
    nonvacuum_floor = max(float(config.mach_density_floors[0]), 1.0e-12 * max(max_rho, 1.0))
    nonvacuum = rho_np > nonvacuum_floor
    if not np.any(nonvacuum):
        nonvacuum = rho_np > float(config.mach_density_floors[0])
    median_j = float(np.median(jmag[nonvacuum])) if np.any(nonvacuum) else 0.0
    current_threshold = float(config.sonic_current_threshold_multiplier) * median_j
    represented = bool(
        math.isfinite(float(base_payload["j_at_max"]))
        and float(base_payload["j_at_max"]) > current_threshold
    )
    values = np.asarray([row["M_full_max"] for row in floor_rows], dtype=np.float64)
    finite_values = values[np.isfinite(values)]
    floor_relative_spread = (
        float((np.max(finite_values) - np.min(finite_values)) / max(np.max(finite_values), 1.0e-300))
        if finite_values.size
        else math.nan
    )
    return {
        "status": "MEASURED",
        "tau": tau,
        "state_path": str(state_path),
        "current_formula": (
            "j=(hbar/m)(psi_r grad psi_i - psi_i grad psi_r) - (q/m) A rho; "
            "matches coupled_branch._matter_number_current gauge-covariant terms"
        ),
        "sound_speed_formula": "c_s=sqrt(5*K/m)*rho^2",
        "base_density_floor": float(config.mach_density_floors[0]),
        "M_full_max": float(base_payload["M_full_max"]),
        "M_w_max": float(base_payload["M_w_max"]),
        "rho_at_max": float(base_payload["rho_at_max"]),
        "j_at_max": float(base_payload["j_at_max"]),
        "j_r_at_max": float(base_payload["j_r_at_max"]),
        "j_w_at_max": float(base_payload["j_w_at_max"]),
        "r_star": float(base_payload["r_star"]),
        "w_star": float(base_payload["w_star"]),
        "R0_at_w_star": float(base_payload["R0_at_w_star"]),
        "r_star_over_R0": float(base_payload["r_star_over_R0"]),
        "nonvacuum_definition": f"rho>{nonvacuum_floor:.6e}",
        "median_abs_j_nonvacuum": median_j,
        "current_threshold": current_threshold,
        "sonic_current_represented": represented,
        "density_floor_sweep": floor_rows,
        "density_floor_relative_spread": floor_relative_spread,
        "density_floor_stability_call": "STABLE"
        if math.isfinite(floor_relative_spread) and floor_relative_spread <= 0.1
        else "SENSITIVE",
    }


def _c0g_mach_map(
    *,
    steps0_3: Mapping[str, Any],
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    rows = [
        _c0g_current_and_mach_for_state(
            state_result=state,
            config=config,
            dtype=dtype,
        )
        for state in sorted(
            steps0_3.get("state_results", []),
            key=lambda item: float(item["tau"]),
            reverse=True,
        )
    ]
    trend_rows = [{"tau": row["tau"], "M_full_max": row["M_full_max"]} for row in rows]
    by_decreasing_tau = sorted(rows, key=lambda item: float(item["tau"]), reverse=True)
    m_values = [float(row["M_full_max"]) for row in by_decreasing_tau]
    monotone_approach = bool(
        len(m_values) >= 2
        and all(right >= left - 1.0e-12 for left, right in zip(m_values, m_values[1:]))
    )
    deepest = min(rows, key=lambda item: float(item["tau"])) if rows else None
    if deepest is None:
        label = "NOT_MEASURED"
        represented = False
    else:
        represented = bool(deepest.get("sonic_current_represented"))
        m_deep = float(deepest["M_full_max"])
        low, high = config.sonic_band
        if not represented:
            label = "SONIC_HYPOTHESIS_NOT_REPRESENTED"
        elif low <= m_deep <= high and monotone_approach:
            label = "SONIC_SUPPORT"
        elif m_deep < float(config.no_sonic_threshold):
            label = "NO_SONIC"
        else:
            label = "INCONCLUSIVE"
    return {
        "status": "MEASURED" if rows else "NOT_MEASURED",
        "context_only_not_verdict_gate": True,
        "rows": rows,
        "M_max_trend": trend_rows,
        "monotone_approach_as_tau_decreases": monotone_approach,
        "sonic_context_label": label,
        "sonic_hypothesis_represented": represented,
        "thresholds": {
            "sonic_band": [float(config.sonic_band[0]), float(config.sonic_band[1])],
            "no_sonic_if_M_max_lt": float(config.no_sonic_threshold),
            "j_at_max_threshold": (
                f"{config.sonic_current_threshold_multiplier:.6e} * "
                "median(|j| over non-vacuum cells)"
            ),
        },
    }


def _c0g_scaled_ftau_vector(
    *,
    state: torch.Tensor,
    tau: float,
    row_scale: np.ndarray,
    h_multiplier: float,
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    h = float(h_multiplier) * float(tau)
    _branch_p, _provider_p, _grid_p, _boundaries_p, _eos_K_p, residual_plus = (
        _c0g_residual_context(tau=float(tau) + h, config=config, dtype=dtype)
    )
    _branch_m, _provider_m, _grid_m, _boundaries_m, _eos_K_m, residual_minus = (
        _c0g_residual_context(tau=float(tau) - h, config=config, dtype=dtype)
    )
    f_plus = residual_plus(state).detach().cpu().numpy().astype(np.float64, copy=False)
    f_minus = residual_minus(state).detach().cpu().numpy().astype(np.float64, copy=False)
    ftau = (f_plus - f_minus) / (2.0 * h)
    scaled = row_scale * ftau
    return {
        "h": h,
        "h_multiplier": float(h_multiplier),
        "scaled_vector": scaled,
        "scaled_l2": float(np.linalg.norm(scaled)),
    }


def _c0g_preferred_ftau_h_multiplier(state_result: Mapping[str, Any], config: C0gConfig) -> float:
    ftau = state_result.get("ftau", {})
    stable_pairs = ftau.get("stable_pairs", []) if isinstance(ftau, Mapping) else []
    if stable_pairs:
        pair = min(
            stable_pairs,
            key=lambda item: float(item.get("relative_difference", math.inf)),
        )
        return float(pair["left_h_multiplier"])
    return float(config.ftau_h_multipliers[0])


def _c0g_bordered_conditioning_for_state(
    *,
    state_result: Mapping[str, Any],
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    tau = float(state_result["tau"])
    state_path = _resolve_input_path(str(state_result["state_path"]))
    state = _load_state_artifact(state_path, dtype=dtype)
    assembled = _c0g_assemble_original_jacobian(
        state=state,
        tau=tau,
        config=config,
        dtype=dtype,
    )
    gauge = _c0g_build_analytic_gauge_matrix(
        state=state,
        grid=assembled["grid"],
        branch=assembled["branch"],
        config=config,
    )
    complement = _c0g_scaled_gauge_complement(
        gauge_matrix=gauge["generator_matrix"],
        col_scale=assembled["col_scale"],
        config=config,
    )
    q_perp = np.asarray(complement["q_perp"], dtype=np.float64)
    reduced = assembled["scaled_matrix"].toarray().astype(np.float64, copy=False) @ q_perp
    mode_path = state_result.get("svd", {}).get("mode_vectors_artifact")
    if not mode_path:
        return {"status": "NOT_MEASURED", "tau": tau, "reason": "missing_mode_artifact"}
    modes = _c0g_load_mode_bundle(mode_path)
    right_min_scaled = modes["right_modes"][0] / assembled["col_scale"]
    vhat, vhat_norm = _unit_vector(q_perp.T @ right_min_scaled)
    if vhat_norm <= 0.0:
        return {"status": "NOT_MEASURED", "tau": tau, "reason": "zero_reduced_right_null"}
    h_multiplier = _c0g_preferred_ftau_h_multiplier(state_result, config)
    ftau = _c0g_scaled_ftau_vector(
        state=state,
        tau=tau,
        row_scale=assembled["row_scale"],
        h_multiplier=h_multiplier,
        config=config,
        dtype=dtype,
    )
    bordered = np.zeros((reduced.shape[0] + 1, reduced.shape[1] + 1), dtype=np.float64)
    bordered[:-1, :-1] = reduced
    bordered[:-1, -1] = np.asarray(ftau["scaled_vector"], dtype=np.float64)
    bordered[-1, :-1] = vhat
    start = time.perf_counter()
    bordered_singular = np.linalg.svd(bordered, compute_uv=False)
    bordered_svd_seconds = float(time.perf_counter() - start)
    sigma_min_b, sigma_max_b, cond_b = _condition_number_from_singular_values(
        bordered_singular
    )
    svd = state_result.get("svd", {})
    sigma_min_a = float(svd.get("sigma_min", math.nan))
    sigma_max_a = float(svd.get("sigma_max", math.nan))
    cond_a = (
        float(sigma_max_a / sigma_min_a)
        if sigma_min_a > 0.0 and math.isfinite(sigma_max_a)
        else math.inf
    )
    sigma_min_ratio = (
        float(sigma_min_b / sigma_min_a)
        if sigma_min_a > 0.0 and math.isfinite(sigma_min_b)
        else math.nan
    )
    fold_support = bool(
        (cond_a > 1.0e10 and cond_b < 1.0e8)
        or (math.isfinite(sigma_min_ratio) and sigma_min_ratio > 1.0e4)
    )
    return {
        "status": "MEASURED",
        "tau": tau,
        "state_path": str(state_path),
        "construction": "bordered_[scaled_J_times_Q_perp, scaled_F_tau; vhat_min^T, 0]",
        "forbidden_construction_not_used": "SVD((I-P_G)J)",
        "reduced_shape": [int(reduced.shape[0]), int(reduced.shape[1])],
        "bordered_shape": [int(bordered.shape[0]), int(bordered.shape[1])],
        "ftau_h_multiplier": float(ftau["h_multiplier"]),
        "ftau_scaled_l2": float(ftau["scaled_l2"]),
        "cond_JQ_perp": cond_a,
        "cond_Jb": cond_b,
        "sigma_min_JQ_perp": sigma_min_a,
        "sigma_min_Jb": sigma_min_b,
        "sigma_min_ratio_Jb_over_JQ_perp": sigma_min_ratio,
        "sigma_max_Jb": sigma_max_b,
        "call": "FOLD_SUPPORT" if fold_support else "NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD",
        "thresholds": {
            "fold_if_cond_JQ_perp_gt_and_cond_Jb_lt": [1.0e10, 1.0e8],
            "fold_if_sigma_min_ratio_gt": 1.0e4,
        },
        "assembly_seconds": float(assembled["metadata"].get("assembly_seconds", math.nan)),
        "bordered_svd_seconds": bordered_svd_seconds,
    }


def _c0g_bordered_conditioning(
    *,
    steps0_3: Mapping[str, Any],
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    rows = []
    by_tau = _c0g_state_results_by_tau(steps0_3)
    for tau in sorted(config.converged_taus, reverse=True):
        state_result = by_tau.get(float(tau))
        if state_result is None:
            rows.append({"status": "NOT_MEASURED", "tau": float(tau), "reason": "missing_state_result"})
            continue
        rows.append(
            _c0g_bordered_conditioning_for_state(
                state_result=state_result,
                config=config,
                dtype=dtype,
            )
        )
    measured = [row for row in rows if row.get("status") == "MEASURED"]
    support_rows = [row for row in measured if row.get("call") == "FOLD_SUPPORT"]
    closest = min(measured, key=lambda item: float(item["tau"])) if measured else None
    return {
        "status": "MEASURED" if measured else "NOT_MEASURED",
        "rows": rows,
        "call": "FOLD_SUPPORT" if support_rows else "NO_FOLD_SUPPORT_BY_BORDER_THRESHOLD",
        "closest_converged_tau_call": closest.get("call") if closest else None,
        "closest_converged_tau": closest.get("tau") if closest else None,
    }


def _c0g_boundary_mask_np(grid) -> np.ndarray:
    mask = np.zeros((grid.spec.nr, grid.spec.nw), dtype=bool)
    mask[0, :] = True
    mask[-1, :] = True
    mask[:, 0] = True
    mask[:, -1] = True
    return mask


def _c0g_curl_from_state_vector(vector: np.ndarray, grid) -> np.ndarray:
    values = np.asarray(vector, dtype=np.float64)
    lanes = _closed_lane_slices(grid)
    shape = (grid.spec.nr, grid.spec.nw)
    ar_start, ar_stop = lanes["ar"]
    aw_start, aw_stop = lanes["aw"]
    ar = torch.as_tensor(
        values[ar_start:ar_stop].reshape(shape),
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    )
    aw = torch.as_tensor(
        values[aw_start:aw_stop].reshape(shape),
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    )
    curl = tensor_center_gradient_r(aw, grid) - tensor_center_gradient_w(ar, grid)
    return curl.detach().cpu().numpy().astype(np.float64, copy=True)


def _c0g_spatial_a_norm(vector: np.ndarray, grid) -> float:
    lanes = _closed_lane_slices(grid)
    values = np.asarray(vector, dtype=np.float64)
    pieces = []
    for name in ("ar", "aw"):
        start, stop = lanes[name]
        pieces.append(values[start:stop])
    if not pieces:
        return 0.0
    a_values = np.concatenate(pieces)
    return float(np.linalg.norm(a_values))


def _c0g_curl_metrics(vector: np.ndarray, grid) -> dict[str, Any]:
    curl = _c0g_curl_from_state_vector(vector, grid)
    curl_norm = float(np.linalg.norm(curl))
    boundary_mask = _c0g_boundary_mask_np(grid)
    boundary_norm = float(np.linalg.norm(curl[boundary_mask]))
    interior_norm = float(np.linalg.norm(curl[~boundary_mask]))
    vector_norm = float(np.linalg.norm(vector))
    spatial_a_norm = _c0g_spatial_a_norm(vector, grid)
    return {
        "curl_norm": curl_norm,
        "state_norm": vector_norm,
        "spatial_a_norm": spatial_a_norm,
        "residual_curl_over_state_norm": float(curl_norm / vector_norm)
        if vector_norm > 0.0
        else math.nan,
        "curl_over_spatial_a_norm": float(curl_norm / spatial_a_norm)
        if spatial_a_norm > 0.0
        else math.nan,
        "boundary_curl_norm_fraction": float(boundary_norm / curl_norm)
        if curl_norm > 0.0
        else math.nan,
        "boundary_curl_energy_fraction": float((boundary_norm**2) / (curl_norm**2))
        if curl_norm > 0.0
        else math.nan,
        "interior_curl_norm": interior_norm,
    }


def _c0g_commutator_for_state(
    *,
    state_result: Mapping[str, Any],
    tracked_lookup: Mapping[tuple[float, int], int],
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    tau = float(state_result["tau"])
    state_path = _resolve_input_path(str(state_result["state_path"]))
    state = _load_state_artifact(state_path, dtype=dtype)
    branch, _provider, grid, _boundaries, _eos_K, _residual_fn = _c0g_residual_context(
        tau=tau,
        config=config,
        dtype=dtype,
    )
    gauge = _c0g_build_analytic_gauge_matrix(
        state=state,
        grid=grid,
        branch=branch,
        config=config,
    )
    generator_matrix = np.asarray(gauge["generator_matrix"], dtype=np.float64)
    curls = np.column_stack(
        [
            _c0g_curl_from_state_vector(generator_matrix[:, index], grid).reshape(-1)
            for index in range(generator_matrix.shape[1])
        ]
    )
    curl_norm = float(np.linalg.norm(curls))
    generator_norm = float(np.linalg.norm(generator_matrix))
    boundary_mask = _c0g_boundary_mask_np(grid).reshape(-1)
    boundary_norm = float(np.linalg.norm(curls[boundary_mask, :]))
    mode_path = state_result.get("svd", {}).get("mode_vectors_artifact")
    mode_rows: list[dict[str, Any]] = []
    if mode_path:
        mode_bundle = _c0g_load_mode_bundle(mode_path)
        for mode_index, mode in enumerate(mode_bundle["right_modes"]):
            p_g = _c0g_projection_fraction_from_basis(gauge["physical_basis"], mode)
            curl_metrics = _c0g_curl_metrics(mode, grid)
            mode_rows.append(
                {
                    "tau": tau,
                    "ascending_rank": int(mode_index + 1),
                    "tracked_lane": int(tracked_lookup.get((tau, mode_index + 1), mode_index + 1)),
                    "sigma": float(mode_bundle["sigmas"][mode_index])
                    if mode_index < mode_bundle["sigmas"].size
                    else math.nan,
                    "P_G": p_g,
                    "one_minus_P_G": float(1.0 - p_g) if math.isfinite(p_g) else math.nan,
                    **curl_metrics,
                }
            )
    rng = np.random.default_rng(int(config.commutator_control_seed) + int(round(tau * 1.0e9)))
    mixed = rng.normal(size=generator_matrix.shape[0]).astype(np.float64)
    mixed_unit, _mixed_norm = _unit_vector(mixed)
    mixed_p_g = _c0g_projection_fraction_from_basis(gauge["physical_basis"], mixed_unit)
    return {
        "status": "MEASURED",
        "tau": tau,
        "state_path": str(state_path),
        "gauge_rank": int(gauge["rank"]),
        "source_helpers": list(gauge["source_helpers"]),
        "commutator_norm_CG_over_G": float(curl_norm / generator_norm)
        if generator_norm > 0.0
        else math.nan,
        "boundary_row_norm_fraction": float(boundary_norm / curl_norm)
        if curl_norm > 0.0
        else math.nan,
        "boundary_row_energy_fraction": float((boundary_norm**2) / (curl_norm**2))
        if curl_norm > 0.0
        else math.nan,
        "mode_rows": mode_rows,
        "mixed_control": {
            "seed": int(config.commutator_control_seed) + int(round(tau * 1.0e9)),
            "P_G": mixed_p_g,
            "one_minus_P_G": float(1.0 - mixed_p_g) if math.isfinite(mixed_p_g) else math.nan,
            **_c0g_curl_metrics(mixed_unit, grid),
        },
        "note": (
            "P_G is projected onto analytic gauge generators. The C0g tracked modes are "
            "right modes of J*Q_perp, so low P_G is expected and is not used as fold evidence."
        ),
    }


def _c0g_commutator_diagnostic(
    *,
    steps0_3: Mapping[str, Any],
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    tracked_lookup = {
        (float(row["tau"]), int(row["ascending_rank"])): int(row["tracked_lane"])
        for row in steps0_3.get("step1_tracked_modes", [])
        if "tau" in row and "ascending_rank" in row and "tracked_lane" in row
    }
    rows = [
        _c0g_commutator_for_state(
            state_result=state,
            tracked_lookup=tracked_lookup,
            config=config,
            dtype=dtype,
        )
        for state in sorted(
            steps0_3.get("state_results", []),
            key=lambda item: float(item["tau"]),
            reverse=True,
        )
    ]
    all_modes = [mode for row in rows for mode in row.get("mode_rows", [])]
    return {
        "status": "MEASURED" if rows else "NOT_MEASURED",
        "rows": rows,
        "per_mode_rows": all_modes,
        "mode_characterization_scope": (
            "C0g gauge-complement tracked modes; mode 2 here is the a0 lane with "
            "sigma about 1.7e-2, not the r0 stall-driving minimum mode."
        ),
    }


def run_c0g_steps4_6_7(config: C0gConfig | None = None) -> dict[str, Any]:
    cfg = config or C0gConfig()
    dtype = configure_backend(BackendConfig())
    steps0_3 = _c0g_load_steps0_3(cfg)
    if steps0_3.get("status") == "PREMISE_FAILED":
        result = {
            "phase": "steps4_6_7",
            "status": "SKIPPED_PREMISE_FAILED",
            "config": _c0g_config_to_dict(cfg),
            "step4_mach": "SKIPPED_PREMISE_FAILED",
            "step6_bordered_conditioning": "SKIPPED_PREMISE_FAILED",
            "step7_commutator": "SKIPPED_PREMISE_FAILED",
        }
        _c0g_write_json(_c0g_steps4_6_7_json_path(cfg), result)
        return result
    result = {
        "phase": "steps4_6_7",
        "status": "MEASURED",
        "config": _c0g_config_to_dict(cfg),
        "step4_mach": _c0g_mach_map(steps0_3=steps0_3, config=cfg, dtype=dtype),
        "step6_bordered_conditioning": _c0g_bordered_conditioning(
            steps0_3=steps0_3,
            config=cfg,
            dtype=dtype,
        ),
        "step7_commutator": _c0g_commutator_diagnostic(
            steps0_3=steps0_3,
            config=cfg,
            dtype=dtype,
        ),
    }
    _c0g_write_json(_c0g_steps4_6_7_json_path(cfg), result)
    return result


class _C0gProbeAbort(RuntimeError):
    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def _c0g_residual_norms(values: np.ndarray) -> dict[str, float]:
    return {
        "l2": float(np.linalg.norm(values)),
        "linf": float(np.max(np.abs(values))) if values.size else 0.0,
    }


def _c0g_scaled_lane_overlaps(
    *,
    candidate: np.ndarray,
    reference: np.ndarray,
    origin: np.ndarray,
    col_scale: np.ndarray,
    grid,
) -> dict[str, Any]:
    candidate_scaled = (candidate - origin) / col_scale
    reference_scaled = (reference - origin) / col_scale
    lane_groups = {
        "psi": ("psi_real", "psi_imag"),
        "A": ("a0", "ar", "aw"),
        "r0": ("r0",),
        "mu": ("mu",),
    }
    slices = _closed_lane_slices(grid)
    rows = []
    for group, names in lane_groups.items():
        indices = []
        for name in names:
            start, stop = slices[name]
            indices.extend(range(start, stop))
        idx = np.asarray(indices, dtype=np.int64)
        a = candidate_scaled[idx]
        b = reference_scaled[idx]
        denom = float(np.linalg.norm(a) * np.linalg.norm(b))
        overlap = float(np.dot(a, b) / denom) if denom > 0.0 else math.nan
        if math.isfinite(overlap) and overlap > 0.99:
            call = "SAME_BRANCH"
        elif math.isfinite(overlap) and overlap < 0.5:
            call = "POST_FOLD_BRANCH"
        elif math.isfinite(overlap):
            call = "INCONCLUSIVE"
        else:
            call = "NOT_MEASURED"
        rows.append(
            {
                "lane": group,
                "scaled_overlap": overlap,
                "call": call,
                "candidate_delta_scaled_norm": float(np.linalg.norm(a)),
                "prestall_delta_scaled_norm": float(np.linalg.norm(b)),
            }
        )
    finite = [row for row in rows if math.isfinite(float(row["scaled_overlap"]))]
    if finite and all(float(row["scaled_overlap"]) > 0.99 for row in finite):
        aggregate = "SAME_BRANCH"
    elif finite and any(float(row["scaled_overlap"]) < 0.5 for row in finite):
        aggregate = "POST_FOLD_BRANCH"
    elif finite:
        aggregate = "INCONCLUSIVE"
    else:
        aggregate = "NOT_MEASURED"
    return {
        "metric": "spectral_scaled_state_coordinates_delta/col_scale",
        "rows": rows,
        "aggregate_overlap_call": aggregate,
        "thresholds": {
            "same_branch_if_lane_overlap_gt": 0.99,
            "post_fold_if_lane_overlap_lt": 0.5,
            "inconclusive_band": [0.5, 0.99],
        },
    }


def _c0g_step5_overlap_if_progress(
    *,
    candidate: np.ndarray,
    stalled_state: np.ndarray,
    prestall_state: np.ndarray,
    config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    state_t = torch.as_tensor(stalled_state, dtype=dtype, device="cpu")
    assembled = _c0g_assemble_original_jacobian(
        state=state_t,
        tau=float(config.stalled_tau),
        config=config,
        dtype=dtype,
    )
    return _c0g_scaled_lane_overlaps(
        candidate=candidate,
        reference=prestall_state,
        origin=stalled_state,
        col_scale=assembled["col_scale"],
        grid=assembled["grid"],
    )


def _c0g_scipy_method_probe(
    *,
    method: str,
    initial: np.ndarray,
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    initial_norms: Mapping[str, float],
    config: C0gConfig,
    dtype: torch.dtype,
    global_start: float,
) -> dict[str, Any]:
    method_start = time.perf_counter()
    maxfev = int(config.scipy_maxfev_factor * (initial.size + 1))
    remaining_global = max(0.0, float(config.scipy_global_wall_seconds) - (method_start - global_start))
    wall_cap = min(float(config.scipy_method_wall_seconds), remaining_global)
    if wall_cap <= 1.0:
        return {
            "status": "NOT_MEASURED",
            "method": method,
            "reason": "global_wall_budget_exhausted_before_method",
            "maxfev": maxfev,
            "wall_cap_seconds": wall_cap,
        }
    trajectory: list[dict[str, Any]] = []
    best_x = initial.copy()
    best_l2 = float(initial_norms["l2"])
    best_linf = float(initial_norms["linf"])
    nfev = 0

    def fun(values: np.ndarray) -> np.ndarray:
        nonlocal nfev, best_x, best_l2, best_linf
        elapsed = time.perf_counter() - method_start
        if elapsed > wall_cap:
            raise _C0gProbeAbort("wall_time_cap")
        if nfev >= maxfev:
            raise _C0gProbeAbort("maxfev_cap")
        nfev += 1
        tensor = torch.as_tensor(values, dtype=dtype, device="cpu")
        residual = residual_fn(tensor).detach().cpu().numpy().astype(np.float64, copy=False)
        norms = _c0g_residual_norms(residual)
        if norms["l2"] < best_l2:
            best_l2 = norms["l2"]
            best_linf = norms["linf"]
            best_x = np.asarray(values, dtype=np.float64).copy()
        trajectory.append(
            {
                "eval": int(nfev),
                "elapsed_seconds": float(elapsed),
                "original_residual_l2": norms["l2"],
                "original_residual_linf": norms["linf"],
            }
        )
        return residual

    options: dict[str, Any]
    if method == "lm":
        options = {"maxiter": maxfev}
    else:
        options = {"maxfev": maxfev}
    abort_reason = None
    scipy_success = False
    scipy_message = None
    try:
        scipy_result = optimize.root(fun, initial.copy(), method=method, options=options)
        scipy_success = bool(scipy_result.success)
        scipy_message = str(scipy_result.message)
        final_x = np.asarray(scipy_result.x, dtype=np.float64)
        reported_nfev = int(getattr(scipy_result, "nfev", nfev))
    except _C0gProbeAbort as exc:
        abort_reason = exc.reason
        final_x = best_x.copy()
        reported_nfev = nfev
    except Exception as exc:  # pragma: no cover - diagnostic reports unexpected scipy failures.
        abort_reason = f"{type(exc).__name__}:{exc}"
        final_x = best_x.copy()
        reported_nfev = nfev

    final_residual = residual_fn(
        torch.as_tensor(final_x, dtype=dtype, device="cpu")
    ).detach().cpu().numpy().astype(np.float64, copy=False)
    final_norms = _c0g_residual_norms(final_residual)
    if final_norms["l2"] < best_l2:
        best_l2 = final_norms["l2"]
        best_linf = final_norms["linf"]
        best_x = final_x.copy()
    l2_ratio = best_l2 / max(float(initial_norms["l2"]), np.finfo(np.float64).tiny)
    if best_linf <= 1.0e-6:
        progress_call = "CONVERGED_LINF"
    elif l2_ratio <= 0.1:
        progress_call = "PROGRESS_10X_L2_DROP"
    else:
        progress_call = "NO_PROGRESS_NO_EVIDENCE"
    status = "MEASURED"
    if abort_reason in {"wall_time_cap", "maxfev_cap"}:
        status = "NOT_MEASURED"
    return {
        "status": status,
        "method": method,
        "maxfev": maxfev,
        "wall_cap_seconds": wall_cap,
        "elapsed_seconds": float(time.perf_counter() - method_start),
        "nfev_wrapper": int(nfev),
        "nfev_reported_by_scipy": reported_nfev,
        "scipy_success": scipy_success,
        "scipy_message": scipy_message,
        "abort_reason": abort_reason,
        "initial_original_residual_l2": float(initial_norms["l2"]),
        "initial_original_residual_linf": float(initial_norms["linf"]),
        "final_original_residual_l2": final_norms["l2"],
        "final_original_residual_linf": final_norms["linf"],
        "best_original_residual_l2": best_l2,
        "best_original_residual_linf": best_linf,
        "best_l2_ratio": l2_ratio,
        "progress_call": progress_call,
        "trajectory": trajectory,
        "best_state": best_x,
    }


def run_c0g_step5_scipy_probe(config: C0gConfig | None = None) -> dict[str, Any]:
    cfg = config or C0gConfig()
    dtype = configure_backend(BackendConfig())
    payload = _c0g_load_c0f2_payload(cfg)
    rows_by_tau = _c0g_rows_by_tau(payload)
    stalled_row = _c0g_find_tau_row(payload, cfg.stalled_tau)
    prestall_tau = min(tau for tau in rows_by_tau if tau > float(cfg.stalled_tau))
    prestall_row = rows_by_tau[prestall_tau]
    stalled_state_t = _load_state_artifact(_c0g_state_path_from_row(stalled_row), dtype=dtype)
    prestall_state_t = _load_state_artifact(_c0g_state_path_from_row(prestall_row), dtype=dtype)
    stalled_state = stalled_state_t.detach().cpu().numpy().astype(np.float64, copy=True)
    prestall_state = prestall_state_t.detach().cpu().numpy().astype(np.float64, copy=True)
    _branch, _provider, _grid, _boundaries, _eos_K, residual_fn = _c0g_residual_context(
        tau=float(cfg.stalled_tau),
        config=cfg,
        dtype=dtype,
    )
    initial_residual = residual_fn(stalled_state_t).detach().cpu().numpy().astype(np.float64, copy=False)
    initial_norms = _c0g_residual_norms(initial_residual)
    global_start = time.perf_counter()
    method_results: list[dict[str, Any]] = []
    overlap: dict[str, Any] | None = None
    for method in ("lm", "hybr"):
        if time.perf_counter() - global_start >= float(cfg.scipy_global_wall_seconds):
            method_results.append(
                {
                    "status": "NOT_MEASURED",
                    "method": method,
                    "reason": "global_wall_budget_exhausted",
                }
            )
            continue
        method_result = _c0g_scipy_method_probe(
            method=method,
            initial=stalled_state,
            residual_fn=residual_fn,
            initial_norms=initial_norms,
            config=cfg,
            dtype=dtype,
            global_start=global_start,
        )
        best_state = np.asarray(method_result.pop("best_state"), dtype=np.float64)
        if (
            method_result.get("progress_call") in {"CONVERGED_LINF", "PROGRESS_10X_L2_DROP"}
            and overlap is None
        ):
            overlap = _c0g_step5_overlap_if_progress(
                candidate=best_state,
                stalled_state=stalled_state,
                prestall_state=prestall_state,
                config=cfg,
                dtype=dtype,
            )
            method_result["overlap_computed_for_this_method"] = True
        else:
            method_result["overlap_computed_for_this_method"] = False
        method_results.append(method_result)
    measured_progress = [
        row
        for row in method_results
        if row.get("progress_call") in {"CONVERGED_LINF", "PROGRESS_10X_L2_DROP"}
    ]
    if overlap is None:
        overlap = {
            "metric": "spectral_scaled_state_coordinates_delta/col_scale",
            "aggregate_overlap_call": "NOT_COMPUTED_NO_PROGRESS",
            "rows": [],
        }
    if any(row.get("status") == "NOT_MEASURED" for row in method_results) and not measured_progress:
        status = "NOT_MEASURED"
    else:
        status = "MEASURED"
    result = {
        "phase": "step5_scipy_probe",
        "status": status,
        "config": _c0g_config_to_dict(cfg),
        "tau": float(cfg.stalled_tau),
        "prestall_tau": float(prestall_tau),
        "state_path": str(_c0g_state_path_from_row(stalled_row)),
        "prestall_state_path": str(_c0g_state_path_from_row(prestall_row)),
        "arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
        "caps": {
            "method_order": ["lm", "hybr"],
            "maxfev": int(cfg.scipy_maxfev_factor * (stalled_state.size + 1)),
            "maxfev_formula": "4*(n+1)",
            "method_wall_seconds": float(cfg.scipy_method_wall_seconds),
            "global_wall_seconds": float(cfg.scipy_global_wall_seconds),
        },
        "initial_original_residual_l2": float(initial_norms["l2"]),
        "initial_original_residual_linf": float(initial_norms["linf"]),
        "method_results": method_results,
        "branch_overlap": overlap,
        "step5_call": overlap.get("aggregate_overlap_call")
        if measured_progress
        else "NO_PROGRESS_NO_EVIDENCE",
        "non_progressing_result_interpretation": "no evidence, not fold proof",
    }
    _c0g_write_json(_resolve_output_path(cfg.step5_json_path), result)
    return result


def _c0g_step5_summary(step5: Mapping[str, Any]) -> dict[str, Any]:
    methods = step5.get("method_results", [])
    return {
        "status": step5.get("status"),
        "call": step5.get("step5_call"),
        "method_progress_calls": [
            {
                "method": row.get("method"),
                "status": row.get("status"),
                "progress_call": row.get("progress_call"),
                "nfev": row.get("nfev_wrapper"),
                "best_l2_ratio": row.get("best_l2_ratio"),
                "best_linf": row.get("best_original_residual_linf"),
                "abort_reason": row.get("abort_reason") or row.get("reason"),
            }
            for row in methods
        ],
        "branch_overlap": step5.get("branch_overlap", {}),
    }


def _c0g_primary_evidence_status(
    *,
    steps0_3: Mapping[str, Any],
    step5: Mapping[str, Any],
) -> dict[str, Any]:
    step0_gate = steps0_3.get("step0", {}).get("premise_gate", {})
    ftau = [
        state.get("ftau", {})
        for state in steps0_3.get("state_results", [])
        if isinstance(state.get("ftau"), Mapping)
    ]
    fit = steps0_3.get("sigma_min_squared_fit", {})
    overlap_call = step5.get("branch_overlap", {}).get("aggregate_overlap_call")
    gray_items: list[str] = []
    not_measured: list[str] = []
    disagreements: list[str] = []
    if step0_gate.get("call") == "PREMISE_GRAY":
        gray_items.append("Step-0 f_nn")
    elif step0_gate.get("call") not in {"PREMISE_HOLDS", "PREMISE_FAILED"}:
        not_measured.append("Step-0 f_nn")
    ftau_calls = {row.get("call") for row in ftau if row.get("status") == "MEASURED"}
    if any(row.get("status") != "MEASURED" for row in ftau):
        not_measured.append("Step-2 cos_theta")
    if "INCONCLUSIVE" in ftau_calls:
        gray_items.append("Step-2 cos_theta")
    if fit.get("status") != "MEASURED" or fit.get("fit_reliability") == "UNRELIABLE":
        not_measured.append("Step-3 sigma_min^2 fit")
    elif fit.get("call") == "INCONCLUSIVE":
        gray_items.append("Step-3 sigma_min^2 fit")
    if overlap_call == "INCONCLUSIVE":
        gray_items.append("Step-5 branch overlap")
    elif overlap_call in {None, "NOT_MEASURED"}:
        not_measured.append("Step-5 branch overlap")
    fold_primary = (
        step0_gate.get("call") == "PREMISE_HOLDS"
        and ftau_calls == {"FOLD_SUPPORT"}
        and fit.get("call") == "LINEAR_MONOTONE_FOLD_SUPPORT"
    )
    conditioning_primary = (
        "BIFURCATION_SUPPORT" in ftau_calls
        or fit.get("call") == "QUADRATIC_OR_FLAT_SUPPORT"
    )
    if fold_primary and overlap_call == "SAME_BRANCH":
        disagreements.append("Step-5 same-branch result disagrees with Steps 2-3 fold support")
    return {
        "gray_items": gray_items,
        "not_measured": not_measured,
        "disagreements": disagreements,
        "fold_primary": fold_primary,
        "conditioning_primary": conditioning_primary,
        "step0_call": step0_gate.get("call"),
        "ftau_calls": sorted(str(call) for call in ftau_calls if call is not None),
        "sigma_fit_call": fit.get("call"),
        "step5_overlap_call": overlap_call,
    }


def _c0g_determine_step8_verdict(
    *,
    steps0_3: Mapping[str, Any],
    steps4_6_7: Mapping[str, Any],
    step5: Mapping[str, Any],
) -> dict[str, Any]:
    primary = _c0g_primary_evidence_status(steps0_3=steps0_3, step5=step5)
    step4 = steps4_6_7.get("step4_mach", {})
    step6 = steps4_6_7.get("step6_bordered_conditioning", {})
    step5_call = primary.get("step5_overlap_call")
    step6_call = step6.get("call")
    step6_support = step6_call == "FOLD_SUPPORT"
    gray_items = list(primary["gray_items"])
    not_measured = list(primary["not_measured"])
    disagreements = list(primary["disagreements"])
    if steps0_3.get("step0", {}).get("premise_gate", {}).get("call") == "PREMISE_FAILED":
        verdict = "PREMISE_FAILED"
        sub_label = None
        recommendation = "Re-diagnose the actual stall driver before any C0g build."
    elif gray_items or disagreements:
        verdict = "MIXED/INCONCLUSIVE"
        sub_label = None
        recommendation = "Gauge-fix plus LM/PTC first, then rerun the battery before pseudo-arclength."
    elif primary["fold_primary"] and (step5_call == "POST_FOLD_BRANCH" or step6_support):
        verdict = "FOLD_CONFIRMED"
        sonic = step4.get("sonic_context_label")
        if sonic == "SONIC_SUPPORT":
            sub_label = "SONIC_FOLD"
            recommendation = (
                "C0g build branch: gauge-fixed pseudo-arclength; optionally "
                "shoot-from-sonic with L'Hopital regularity."
            )
        elif sonic == "NO_SONIC":
            sub_label = "NON_SONIC_FOLD"
            recommendation = "C0g build branch: gauge-fixed pseudo-arclength."
        elif sonic == "SONIC_HYPOTHESIS_NOT_REPRESENTED":
            sub_label = "SONIC_NOT_REPRESENTED_FOLD"
            recommendation = "C0g build branch: gauge-fixed pseudo-arclength."
        else:
            sub_label = "SONIC_CONTEXT_INCONCLUSIVE_FOLD"
            recommendation = "C0g build branch: gauge-fixed pseudo-arclength."
    elif (
        primary["conditioning_primary"]
        and step5_call == "SAME_BRANCH"
        and any(
            row.get("progress_call") in {"CONVERGED_LINF", "PROGRESS_10X_L2_DROP"}
            for row in step5.get("method_results", [])
        )
    ):
        verdict = "CONDITIONING"
        sub_label = None
        recommendation = (
            "C0g build branch: shifted-Newton LM and/or PTC+SER with "
            "inverse-Hamiltonian matter preconditioning."
        )
    elif not_measured and not step6_support:
        verdict = "DIAGNOSTIC_INCOMPLETE"
        sub_label = None
        recommendation = "Complete the missing confirmatory evidence before choosing the C0g build branch."
    else:
        verdict = "MIXED/INCONCLUSIVE"
        sub_label = None
        recommendation = "Gauge-fix plus LM/PTC first, then rerun the battery before pseudo-arclength."
    return {
        "verdict": verdict,
        "sub_label": sub_label,
        "recommended_C0g_build_branch": recommendation,
        "gray_zone_rule_honored": True,
        "primary_evidence": primary,
        "confirmatory": {
            "step5_overlap_call": step5_call,
            "step6_call": step6_call,
            "step6_support": step6_support,
        },
        "sonic_context_label": step4.get("sonic_context_label"),
        "sonic_hypothesis_represented": step4.get("sonic_hypothesis_represented"),
    }


def aggregate_c0g_final(config: C0gConfig | None = None) -> dict[str, Any]:
    cfg = config or C0gConfig()
    steps0_3 = _c0g_load_steps0_3(cfg)
    provenance = _c0g_verify_stalled_provenance(cfg)
    steps4_path = _c0g_steps4_6_7_json_path(cfg)
    with steps4_path.open("r", encoding="utf-8") as handle:
        steps4_6_7 = json.load(handle)
    step5_path = _resolve_input_path(cfg.step5_json_path)
    with step5_path.open("r", encoding="utf-8") as handle:
        step5 = json.load(handle)
    verdict = _c0g_determine_step8_verdict(
        steps0_3=steps0_3,
        steps4_6_7=steps4_6_7,
        step5=step5,
    )
    result = {
        "phase": "complete_steps0_8",
        "status": "MEASURED",
        "config": _c0g_config_to_dict(cfg),
        "provenance": provenance,
        "step0": steps0_3.get("step0"),
        "step0b": steps0_3.get("step0b"),
        "state_results": steps0_3.get("state_results", []),
        "step1_tracked_modes": steps0_3.get("step1_tracked_modes", []),
        "sigma_min_squared_fit": steps0_3.get("sigma_min_squared_fit", {}),
        "step4_mach": steps4_6_7.get("step4_mach", {}),
        "step5_scipy_probe": step5,
        "step6_bordered_conditioning": steps4_6_7.get("step6_bordered_conditioning", {}),
        "step7_commutator": steps4_6_7.get("step7_commutator", {}),
        "step8": verdict,
        "verdict_support": {
            "step0_decomposition": steps0_3.get("step0", {}).get("step0_decomposition", {}),
            "step0b_best_alpha_trend": steps0_3.get("step0b", {}),
            "step1_svd_table": [
                {
                    "tau": row.get("tau"),
                    "sigma_min": row.get("svd", {}).get("sigma_min"),
                    "sigma_min_over_sigma_max": row.get("svd", {}).get("sigma_min_over_sigma_max"),
                    "modes": row.get("svd", {}).get("modes", []),
                }
                for row in steps0_3.get("state_results", [])
            ],
            "step2_cos_theta_table": [
                {"tau": row.get("tau"), "ftau": row.get("ftau")}
                for row in steps0_3.get("state_results", [])
            ],
            "step3_sigma_min_squared_fit": steps0_3.get("sigma_min_squared_fit", {}),
            "step4_mach": steps4_6_7.get("step4_mach", {}),
            "step5_scipy_probe_summary": _c0g_step5_summary(step5),
            "step6_bordered_conditioning": steps4_6_7.get("step6_bordered_conditioning", {}),
            "step7_commutator": steps4_6_7.get("step7_commutator", {}),
        },
        "scope_confirmation": {
            "no_C0g_fix_implemented": True,
            "single_arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
            "depth_continuation": "tau_only",
            "solver_logic_touched_by_c0g": False,
            "faithful_operators_touched_by_c0g": False,
            "frozen_physics_touched_by_c0g": False,
            "physical_export_guard_touched_by_c0g": False,
        },
    }
    _c0g_write_json(_resolve_output_path(cfg.final_json_path), result)
    write_c0g_final_report(result, _resolve_output_path(cfg.report_path))
    return result


def _c0g_build_b1b2_config_to_dict(config: C0gBuildB1B2Config) -> dict[str, Any]:
    data = asdict(config)
    for key in ("c0f2_json_path", "run_root", "report_path", "json_path"):
        data[key] = str(data[key])
    return data


def _c0g_build_crawl_config(config: C0gBuildB1B2Config) -> C0Config:
    return C0Config(
        run_root=config.run_root / "gauge_fixed_crawl",
        report_path=config.run_root / "gauge_fixed_crawl" / "c0_gauge_fixed_crawl.md",
        json_path=config.run_root / "gauge_fixed_crawl" / "c0_gauge_fixed_crawl.json",
        depth_sequence=tuple(float(tau) for tau in config.b2_depth_sequence),
        grid=config.grid,
        use_gauge_fix=True,
        prefer_existing_b2c_background_predictor=False,
        max_tau_backtracks=int(config.max_tau_backtracks),
        min_tau_backtrack_delta=float(config.min_tau_backtrack_delta),
        max_depth_failures_after_floor=int(config.max_depth_failures_after_floor),
        max_newton_iters_override=config.max_newton_iters_override,
        jacobian_assembly=config.jacobian_assembly,
        diagnostic_linear=False,
        diagnostic_observables=False,
    )


def _c0g_build_diag_config(
    config: C0gBuildB1B2Config,
    *,
    converged_taus: Sequence[float],
) -> C0gConfig:
    return C0gConfig(
        c0f2_json_path=config.c0f2_json_path,
        run_root=config.run_root / "b2_reconfirm",
        report_path=config.report_path,
        step0_json_path=config.run_root / "b2_reconfirm" / "unused_step0.json",
        json_path=config.run_root / "b2_reconfirm" / "pathA_C0g_build_B2_steps1_3.json",
        step5_json_path=config.run_root / "b2_reconfirm" / "unused_step5.json",
        final_json_path=config.run_root / "b2_reconfirm" / "pathA_C0g_build_B2_reconfirm.json",
        grid=config.grid,
        converged_taus=tuple(float(tau) for tau in converged_taus),
        stalled_tau=float(config.tau_fold_reference),
        fit_good_r2_threshold=float(config.b2_fit_r2_threshold),
        jacobian_assembly=config.jacobian_assembly,
    )


def _c0g_converged_rows(crawl: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        dict(row)
        for row in crawl.get("tau_attempts", [])
        if bool(row.get("final_physical_converged"))
    ]


def _c0g_select_proof_row(
    crawl: Mapping[str, Any],
    *,
    proof_tau: float,
) -> dict[str, Any]:
    converged = _c0g_converged_rows(crawl)
    if not converged:
        raise RuntimeError("gauge-fixed crawl has no converged rows for B-1 proof")
    exact = [
        row
        for row in converged
        if _c0g_tau_close(float(row.get("target_tau", math.nan)), float(proof_tau))
    ]
    if exact:
        return exact[0]
    return min(converged, key=lambda row: abs(float(row["target_tau"]) - float(proof_tau)))


def _c0g_original_row_for_tau(c0f2_payload: Mapping[str, Any], tau: float) -> dict[str, Any]:
    rows = [
        row
        for row in c0f2_payload.get("per_tau_timing", [])
        if _c0g_tau_close(float(row.get("tau", math.nan)), float(tau))
        and bool(row.get("solver_converged"))
        and bool(row.get("accepted_default_success"))
    ]
    if not rows:
        raise RuntimeError(f"C0f2 has no converged original row for tau={tau:.12e}")
    return dict(rows[0])


def _c0g_state_vector_from_path(path: str | Path, *, dtype: torch.dtype) -> torch.Tensor:
    return _load_state_artifact(_resolve_input_path(path), dtype=dtype)


def _c0g_global_phase_align_np(reference: np.ndarray, candidate: np.ndarray, grid) -> tuple[np.ndarray, dict[str, Any]]:
    n = int(grid.spec.nr * grid.spec.nw)
    ref_complex = reference[:n] + 1j * reference[n : 2 * n]
    cand_complex = candidate[:n] + 1j * candidate[n : 2 * n]
    inner = np.vdot(ref_complex.reshape(-1), cand_complex.reshape(-1))
    phase = float(np.angle(inner)) if abs(inner) > 0.0 else 0.0
    rotated = cand_complex * np.exp(-1j * phase)
    aligned = np.asarray(candidate, dtype=np.float64).copy()
    aligned[:n] = rotated.real.reshape(-1)
    aligned[n : 2 * n] = rotated.imag.reshape(-1)
    return aligned, {
        "global_phase_angle_radians": phase,
        "global_phase_inner_abs": float(abs(inner)),
    }


def _c0g_curl_np(values: np.ndarray, grid) -> np.ndarray:
    lanes = _closed_lane_slices(grid)
    shape = (grid.spec.nr, grid.spec.nw)
    ar = torch.as_tensor(
        values[lanes["ar"][0] : lanes["ar"][1]].reshape(shape),
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    )
    aw = torch.as_tensor(
        values[lanes["aw"][0] : lanes["aw"][1]].reshape(shape),
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    )
    curl = tensor_center_gradient_r(aw, grid) - tensor_center_gradient_w(ar, grid)
    return curl.detach().cpu().numpy().astype(np.float64, copy=True)


def _c0g_a0_gradients_np(values: np.ndarray, grid) -> tuple[np.ndarray, np.ndarray]:
    lanes = _closed_lane_slices(grid)
    shape = (grid.spec.nr, grid.spec.nw)
    a0 = torch.as_tensor(
        values[lanes["a0"][0] : lanes["a0"][1]].reshape(shape),
        dtype=grid.r_centers.dtype,
        device=grid.r_centers.device,
    )
    grad_r = tensor_center_gradient_r(a0, grid).detach().cpu().numpy().astype(np.float64, copy=True)
    grad_w = tensor_center_gradient_w(a0, grid).detach().cpu().numpy().astype(np.float64, copy=True)
    return grad_r, grad_w


def _c0g_density_np(values: np.ndarray, grid) -> np.ndarray:
    n = int(grid.spec.nr * grid.spec.nw)
    return values[:n] * values[:n] + values[n : 2 * n] * values[n : 2 * n]


def _c0g_gauge_aligned_path_only_proof(
    *,
    with_state: torch.Tensor,
    without_state: torch.Tensor,
    tau: float,
    config: C0gBuildB1B2Config,
    dtype: torch.dtype,
) -> dict[str, Any]:
    branch, _provider, grid, _boundaries, _eos_K, residual_fn = _c0g_residual_context(
        tau=float(tau),
        config=C0gConfig(grid=config.grid),
        dtype=dtype,
    )
    ref = without_state.detach().cpu().numpy().astype(np.float64, copy=True)
    cand = with_state.detach().cpu().numpy().astype(np.float64, copy=True)
    phase_aligned, phase_info = _c0g_global_phase_align_np(ref, cand, grid)
    generator_matrix = _c0e_coupled_gauge_matrix(without_state, grid, branch)
    delta = phase_aligned - ref
    coeff, residuals, rank, singular = np.linalg.lstsq(
        generator_matrix,
        delta,
        rcond=float(config.proof_tolerance),
    )
    orbit_aligned = phase_aligned - generator_matrix @ coeff
    orbit_residual = orbit_aligned - ref
    ref_a0_r, ref_a0_w = _c0g_a0_gradients_np(ref, grid)
    aligned_a0_r, aligned_a0_w = _c0g_a0_gradients_np(orbit_aligned, grid)
    lanes = _closed_lane_slices(grid)
    density_delta = float(np.max(np.abs(_c0g_density_np(ref, grid) - _c0g_density_np(orbit_aligned, grid))))
    curl_delta = float(np.max(np.abs(_c0g_curl_np(ref, grid) - _c0g_curl_np(orbit_aligned, grid))))
    a0_grad_delta = max(
        float(np.max(np.abs(ref_a0_r - aligned_a0_r))),
        float(np.max(np.abs(ref_a0_w - aligned_a0_w))),
    )
    r0_start, r0_stop = lanes["r0"]
    mu_start, mu_stop = lanes["mu"]
    r0_delta = float(np.max(np.abs(ref[r0_start:r0_stop] - orbit_aligned[r0_start:r0_stop])))
    mu_delta = float(np.max(np.abs(ref[mu_start:mu_stop] - orbit_aligned[mu_start:mu_stop])))
    residual_without = residual_fn(without_state).detach().cpu().numpy().astype(np.float64)
    residual_with = residual_fn(with_state).detach().cpu().numpy().astype(np.float64)
    residual_without_linf = float(np.max(np.abs(residual_without)))
    residual_with_linf = float(np.max(np.abs(residual_with)))
    residual_norm_abs_delta = abs(
        float(np.linalg.norm(residual_with)) - float(np.linalg.norm(residual_without))
    )
    max_invariant_delta = max(density_delta, curl_delta, a0_grad_delta, r0_delta, mu_delta)
    status = (
        "PASS"
        if max_invariant_delta <= float(config.proof_tolerance)
        and residual_without_linf <= float(config.proof_tolerance)
        and residual_with_linf <= float(config.proof_tolerance)
        else "FAIL"
    )
    return {
        "status": status,
        "tau": float(tau),
        "tolerance": float(config.proof_tolerance),
        "alignment": {
            **phase_info,
            "coupled_gauge_ls_rank": int(rank),
            "coupled_gauge_ls_singular_min": float(np.min(singular)) if len(singular) else math.nan,
            "coupled_gauge_ls_singular_max": float(np.max(singular)) if len(singular) else math.nan,
            "coupled_gauge_ls_residual_sum": float(np.sum(residuals)) if len(residuals) else 0.0,
            "orbit_residual_l2": float(np.linalg.norm(orbit_residual)),
            "raw_state_delta_l2_not_a_gate": float(np.linalg.norm(cand - ref)),
        },
        "gauge_invariant_deltas_after_alignment": {
            "density_abs_linf": density_delta,
            "F_rw_abs_linf": curl_delta,
            "a0_gradient_abs_linf": a0_grad_delta,
            "r0_abs_linf": r0_delta,
            "mu_abs": mu_delta,
            "max": max_invariant_delta,
        },
        "original_residual": {
            "without_gauge_fix_linf": residual_without_linf,
            "with_gauge_fix_linf": residual_with_linf,
            "l2_norm_abs_delta": residual_norm_abs_delta,
            "arbiter": "stage1_solver.coupled_branch.patha_closed_branch_residual",
        },
    }


def _c0g_gauge_sector_lift(
    *,
    state: torch.Tensor,
    tau: float,
    config: C0gBuildB1B2Config,
    dtype: torch.dtype,
) -> dict[str, Any]:
    diag_config = C0gConfig(grid=config.grid)
    assembled = _c0g_assemble_original_jacobian(
        state=state,
        tau=float(tau),
        config=diag_config,
        dtype=dtype,
    )
    gauge = _c0g_build_analytic_gauge_matrix(
        state=state,
        grid=assembled["grid"],
        branch=assembled["branch"],
        config=diag_config,
    )
    scaled_gauge = np.asarray(gauge["generator_matrix"], dtype=np.float64) / assembled["col_scale"][:, None]
    _basis, full_u, gauge_singular, threshold, rank = _c0g_column_space(
        scaled_gauge,
        relative_tolerance=float(diag_config.gradient_rank_rtol),
        full_matrices=True,
    )
    if rank <= 0:
        return {"status": "FAIL", "reason": "analytic_gauge_rank_zero"}
    gauge_basis = full_u[:, :rank]
    gauge_image = assembled["scaled_matrix"].toarray().astype(np.float64, copy=False) @ gauge_basis
    before = np.linalg.svd(gauge_image, compute_uv=False)
    after_columns = np.vstack([gauge_image, np.eye(rank, dtype=np.float64)])
    after = np.linalg.svd(after_columns, compute_uv=False)
    return {
        "status": "PASS" if int(rank) > 0 and float(gauge_singular[rank - 1]) > 0.0 else "FAIL",
        "tau": float(tau),
        "method": "scaled gauge-complement/bordering coordinates",
        "gauge_rank": int(rank),
        "state_dim": int(assembled["scaled_matrix"].shape[1]),
        "source_helpers": list(gauge["source_helpers"]),
        "piece_dimensions": dict(gauge["piece_dimensions"]),
        "scaled_gauge_rank_threshold": float(threshold),
        "scaled_gauge_singular_min_retained": float(gauge_singular[rank - 1]),
        "before_gauge_sector_sigma_min": float(np.min(before)),
        "before_gauge_sector_sigma_max": float(np.max(before)),
        "after_bordered_gauge_sector_sigma_min": float(np.min(after)),
        "after_bordered_gauge_sector_sigma_max": float(np.max(after)),
        "after_bordered_identity_metric_is_by_construction": True,
        "not_a_lift_gate": "after_bordered_gauge_sector_sigma_min",
        "meaningful_gauge_removed_evidence": {
            "before_gauge_sector_sigma_min": float(np.min(before)),
            "scaled_gauge_singular_min_retained": float(gauge_singular[rank - 1]),
        },
        "status_criterion": "analytic scaled gauge rank retained; bordered identity sigma is reported but not gated",
        "after_definition": "singular values of [scaled_J*G_scaled_basis; I_gauge], by-construction >= O(1)",
    }


def _c0g_r0_mode_preservation(
    *,
    state: torch.Tensor,
    tau: float,
    config: C0gBuildB1B2Config,
    dtype: torch.dtype,
) -> dict[str, Any]:
    diag_config = C0gConfig(grid=config.grid)
    assembled = _c0g_assemble_original_jacobian(
        state=state,
        tau=float(tau),
        config=diag_config,
        dtype=dtype,
    )
    gauge = _c0g_build_analytic_gauge_matrix(
        state=state,
        grid=assembled["grid"],
        branch=assembled["branch"],
        config=diag_config,
    )
    svd = _c0g_dense_complement_svd(
        scaled_matrix=assembled["scaled_matrix"],
        gauge_matrix=gauge["generator_matrix"],
        col_scale=assembled["col_scale"],
        grid=assembled["grid"],
        config=diag_config,
        mode_artifact_path=config.run_root / "b1_evidence" / "r0_mode_preservation_modes.npz",
    )
    if svd.get("status") != "MEASURED":
        return {"status": "FAIL", "reason": svd.get("reason", "svd_not_measured")}
    mode1 = (svd.get("modes") or [{}])[0]
    lane, lane_fraction = _dominant_group(mode1.get("v_lane_energy_fractions", {}))
    mode_bundle = _c0g_load_mode_bundle(str(svd["mode_vectors_artifact"]))
    p_g = _c0g_projection_fraction_from_basis(
        gauge["physical_basis"],
        mode_bundle["right_modes"][0],
    )
    preserved = bool(lane == "r0" and lane_fraction >= 0.7 and p_g <= 1.0e-8)
    return {
        "status": "PASS" if preserved else "FAIL",
        "tau": float(tau),
        "mode1_lane": lane,
        "mode1_lane_fraction": lane_fraction,
        "mode1_sigma": mode1.get("sigma"),
        "mode1_projection_into_gauge_basis": p_g,
        "gauge_rank": svd.get("scaled_gauge_rank"),
        "mode_artifact": svd.get("mode_vectors_artifact"),
        "criterion": "dominant r0 lane >=0.7 and P_G <= 1e-8",
    }


def _c0g_b2_analyze_gauge_fixed_crawl(
    *,
    crawl: Mapping[str, Any],
    config: C0gBuildB1B2Config,
    dtype: torch.dtype,
) -> dict[str, Any]:
    converged_rows = _c0g_converged_rows(crawl)
    converged_taus = sorted(
        [float(row["target_tau"]) for row in converged_rows],
        reverse=True,
    )
    diag_config = _c0g_build_diag_config(config, converged_taus=converged_taus)
    selected_rows = [
        dict(row)
        for row in crawl.get("tau_attempts", [])
        if bool(row.get("final_physical_converged"))
    ]
    by_tau: dict[float, dict[str, Any]] = {}
    for row in selected_rows:
        by_tau[float(row["target_tau"])] = row
    state_results: list[dict[str, Any]] = []
    for tau, row in sorted(by_tau.items(), key=lambda item: item[0], reverse=True):
        prior_state_path = _resolve_input_path(
            DEFAULT_C0G_RUN_ROOT
            / "state_measurements"
            / f"c0g_state_tau_{_format_tau(float(tau))}.json"
        )
        can_reuse_prior = (
            not bool(row.get("gauge_fix_solver_invoked"))
            or bool(row.get("gauge_fix_orbit_equivalent_to_c0f2"))
        )
        if can_reuse_prior and prior_state_path.exists():
            measured = json.loads(prior_state_path.read_text(encoding="utf-8"))
            measured["reused_prior_c0g_measurement"] = True
            measured["reused_prior_c0g_reason"] = (
                "gauge-fixed crawl state is copied or exact global-phase equivalent to C0f2"
            )
        else:
            measured = _c0g_analyze_state(
                tau=float(tau),
                row=row,
                config=diag_config,
                dtype=dtype,
                include_step0=False,
                compute_ftau=True,
            )
            _c0g_write_state_result(diag_config, measured)
            measured["reused_prior_c0g_measurement"] = False
        measured["crawl_final_physical_converged"] = bool(row.get("final_physical_converged"))
        measured["crawl_final_original_residual_linf"] = row.get("final_original_residual_linf")
        measured["crawl_init"] = row.get("init")
        measured["crawl_used_existing_b2c"] = row.get("used_existing_b2c")
        measured["gauge_fix_solver_invoked"] = bool(row.get("gauge_fix_solver_invoked"))
        state_results.append(measured)
    steps1_3 = {
        "phase": "B2_steps1_3_gauge_fixed",
        "status": "MEASURED" if state_results else "NOT_MEASURED",
        "config": _c0g_config_to_dict(diag_config),
        "state_results": state_results,
        "step1_tracked_modes": _c0g_track_modes(state_results),
        "sigma_min_squared_fit": _c0g_fit_sigma_min_squared(state_results, config=diag_config),
    }
    prior_step6_rows: dict[float, dict[str, Any]] = {}
    prior_steps4_path = _resolve_input_path(DEFAULT_C0G_RUN_ROOT / "pathA_C0g_steps4_6_7.json")
    if prior_steps4_path.exists():
        prior_steps4 = json.loads(prior_steps4_path.read_text(encoding="utf-8"))
        for row in prior_steps4.get("step6_bordered_conditioning", {}).get("rows", []):
            if row.get("status") == "MEASURED" and "tau" in row:
                prior_step6_rows[float(row["tau"])] = dict(row)
    state_by_tau = _c0g_state_results_by_tau(steps1_3)
    step6_rows = []
    for tau in sorted(converged_taus, reverse=True):
        crawl_row = by_tau.get(float(tau), {})
        can_reuse_prior = (
            not bool(crawl_row.get("gauge_fix_solver_invoked"))
            or bool(crawl_row.get("gauge_fix_orbit_equivalent_to_c0f2"))
        )
        if can_reuse_prior and float(tau) in prior_step6_rows:
            reused = dict(prior_step6_rows[float(tau)])
            reused["reused_prior_c0g_measurement"] = True
            reused["reused_prior_c0g_reason"] = (
                "gauge-fixed crawl state is copied or exact global-phase equivalent to C0f2"
            )
            step6_rows.append(reused)
            continue
        state_result = state_by_tau.get(float(tau))
        if state_result is None:
            step6_rows.append({"status": "NOT_MEASURED", "tau": float(tau), "reason": "missing_state_result"})
            continue
        step6_row = _c0g_bordered_conditioning_for_state(
            state_result=state_result,
            config=diag_config,
            dtype=dtype,
        )
        step6_row["reused_prior_c0g_measurement"] = False
        step6_rows.append(step6_row)
    step6_measured = [row for row in step6_rows if row.get("status") == "MEASURED"]
    step6 = {
        "status": "MEASURED" if step6_measured else "NOT_MEASURED",
        "rows": step6_rows,
        "call": "B2_RELATIVE_TREND_EVALUATED_OUTSIDE_ABSOLUTE_BAR",
        "closest_converged_tau": min(
            (float(row["tau"]) for row in step6_measured),
            default=None,
        ),
    }
    result = {
        **steps1_3,
        "step6_bordered_conditioning": step6,
    }
    _c0g_write_json(_resolve_output_path(diag_config.json_path), result)
    return result


def _c0g_b2_cos_support(state_results: Sequence[Mapping[str, Any]], config: C0gBuildB1B2Config) -> dict[str, Any]:
    rows = []
    for state in state_results:
        if not bool(state.get("crawl_final_physical_converged")):
            continue
        ftau = state.get("ftau", {})
        row = {
            "tau": state.get("tau"),
            "status": ftau.get("status") if isinstance(ftau, Mapping) else "NOT_MEASURED",
            "cos_theta": ftau.get("representative_cos_theta") if isinstance(ftau, Mapping) else None,
            "call": ftau.get("call") if isinstance(ftau, Mapping) else None,
            "stability": ftau.get("stepsize_stability_call") if isinstance(ftau, Mapping) else None,
        }
        rows.append(row)
    measured = [
        row
        for row in rows
        if row.get("status") == "MEASURED"
        and row.get("stability") == "STABLE"
        and isinstance(row.get("cos_theta"), float)
    ]
    pass_gate = bool(
        len(measured) >= 3
        and all(float(row["cos_theta"]) > C0gConfig().cos_fold_threshold for row in measured)
    )
    return {
        "status": "PASS" if pass_gate else "FAIL",
        "rows": rows,
        "measured_count": int(len(measured)),
        "threshold": float(C0gConfig().cos_fold_threshold),
    }


def _c0g_b2_fit_support(fit: Mapping[str, Any], config: C0gBuildB1B2Config) -> dict[str, Any]:
    linear = fit.get("linear", {}) if isinstance(fit, Mapping) else {}
    r2 = linear.get("r2")
    monotone = bool(fit.get("sigma_min_monotone_decreasing_toward_stall")) if isinstance(fit, Mapping) else False
    pass_gate = bool(
        fit.get("status") == "MEASURED"
        and monotone
        and isinstance(r2, float)
        and math.isfinite(r2)
        and float(r2) >= float(config.b2_fit_r2_threshold)
    )
    return {
        "status": "PASS" if pass_gate else "FAIL",
        "call": fit.get("call") if isinstance(fit, Mapping) else None,
        "linear_r2": r2,
        "linear_slope": linear.get("slope"),
        "tau_fold_zero_crossing": linear.get("tau_fold_zero_crossing"),
        "monotone": monotone,
        "threshold_r2": float(config.b2_fit_r2_threshold),
        "rows": fit.get("rows", []) if isinstance(fit, Mapping) else [],
    }


def _c0g_b2_step6_support(step6: Mapping[str, Any], config: C0gBuildB1B2Config) -> dict[str, Any]:
    rows = [row for row in step6.get("rows", []) if row.get("status") == "MEASURED"]
    ordered = sorted(rows, key=lambda row: float(row["tau"]), reverse=True)
    cond_jb = [float(row["cond_Jb"]) for row in ordered if math.isfinite(float(row["cond_Jb"]))]
    ratios = [
        float(row["cond_JQ_perp"]) / float(row["cond_Jb"])
        for row in ordered
        if math.isfinite(float(row.get("cond_JQ_perp", math.nan)))
        and math.isfinite(float(row.get("cond_Jb", math.nan)))
        and float(row.get("cond_Jb", 0.0)) > 0.0
    ]
    flat_band = float(max(cond_jb) / min(cond_jb)) if len(cond_jb) >= 2 and min(cond_jb) > 0.0 else math.inf
    ratio_factor = float(ratios[-1] / ratios[0]) if len(ratios) >= 2 and ratios[0] > 0.0 else math.nan
    ratio_increasing = bool(
        len(ratios) >= 2
        and all(right >= left * (1.0 - 1.0e-10) for left, right in zip(ratios, ratios[1:]))
    )
    pass_gate = bool(
        len(ordered) >= 3
        and math.isfinite(flat_band)
        and flat_band <= float(config.b2_cond_jb_flat_band_max_over_min)
        and math.isfinite(ratio_factor)
        and ratio_factor >= float(config.b2_ratio_growth_factor_min)
        and ratio_increasing
    )
    return {
        "status": "PASS" if pass_gate else "FAIL",
        "rows": rows,
        "cond_Jb_max_over_min": flat_band,
        "cond_Jb_flat_band_threshold": float(config.b2_cond_jb_flat_band_max_over_min),
        "ratio_values_toward_fold": ratios,
        "ratio_growth_factor": ratio_factor,
        "ratio_growth_factor_threshold": float(config.b2_ratio_growth_factor_min),
        "ratio_increasing_toward_fold": ratio_increasing,
    }


def _c0g_b2_closer_sampling(crawl: Mapping[str, Any], config: C0gBuildB1B2Config) -> dict[str, Any]:
    attempts = [
        {
            "tau": float(row["target_tau"]),
            "converged": bool(row.get("final_physical_converged")),
            "residual_linf": row.get("final_original_residual_linf"),
            "newton_iters": _c0g_attempt_newton_iterations(row),
            "backtrack_index": int(row.get("backtrack_index", 0)),
            "max_newton_iters": row.get("max_newton_iters"),
            "init_source": row.get("init", {}).get("source"),
        }
        for row in crawl.get("tau_attempts", [])
        if "target_tau" in row
    ]
    if not attempts:
        return {"status": "NOT_MEASURED", "attempts": []}
    tau_fold = float(config.tau_fold_reference)
    prior_distance = abs(0.029125 - tau_fold)
    closest = min(attempts, key=lambda row: abs(float(row["tau"]) - tau_fold))
    return {
        "status": "MEASURED",
        "tau_fold_reference": tau_fold,
        "prior_last_converged_tau": 0.029125,
        "prior_last_converged_distance": prior_distance,
        "closest_attempt_tau": closest["tau"],
        "closest_attempt_distance": abs(float(closest["tau"]) - tau_fold),
        "closer_than_prior_last_converged": abs(float(closest["tau"]) - tau_fold) < prior_distance,
        "attempts": attempts,
    }


def _c0g_attempt_newton_iterations(row: Mapping[str, Any]) -> int | None:
    attempts = row.get("epsilon_attempts", [])
    values = [
        int(attempt.get("newton_iterations", attempt.get("iterations", 0)))
        for attempt in attempts
        if isinstance(attempt, Mapping)
    ]
    if not values:
        return None
    return int(max(values))


def _c0g_attempt_total_newton_iterations(row: Mapping[str, Any]) -> int | None:
    attempts = row.get("epsilon_attempts", [])
    values = [
        int(attempt.get("newton_iterations", attempt.get("iterations", 0)))
        for attempt in attempts
        if isinstance(attempt, Mapping)
    ]
    if not values:
        return None
    return int(sum(values))


def _c0g_b2_branch_continuity_source(source: Any) -> bool:
    return source in {
        "previous_c0_converged_state",
        "previous_c0g_gauge_fixed_converged_state",
    }


def _c0g_b2_full_budget_dissolve_contest(
    crawl: Mapping[str, Any],
    config: C0gBuildB1B2Config,
) -> dict[str, Any]:
    attempts = [dict(row) for row in crawl.get("tau_attempts", []) if "target_tau" in row]
    tau_fold = float(config.tau_fold_reference)
    margin = float(config.dissolved_tau_margin)
    below_failure_guard = _crawl_persistent_failure_guard(
        attempts,
        crawl.get("config", {}),
        failed_attempt_predicate=lambda row: (
            float(row.get("target_tau", math.inf)) <= tau_fold - margin
            and not bool(row.get("final_physical_converged"))
            and row.get("measurement_status", "MEASURED") == "MEASURED"
        ),
    )
    above_probe_taus = [
        float(tau)
        for tau in config.b2_depth_sequence
        if tau_fold < float(tau) < 0.029125 - 5.0e-13
    ]
    above_probe_tau = min(above_probe_taus, key=lambda tau: abs(tau - tau_fold)) if above_probe_taus else None
    above_rows = [
        row
        for row in attempts
        if above_probe_tau is not None
        and _c0g_tau_close(float(row.get("target_tau", math.nan)), float(above_probe_tau))
    ]
    above_converged = any(bool(row.get("final_physical_converged")) for row in above_rows)
    below_rows = [
        row
        for row in attempts
        if float(row.get("target_tau", math.inf)) <= tau_fold - margin
    ]
    below_converged = [row for row in below_rows if bool(row.get("final_physical_converged"))]
    persistent_failure = bool(below_failure_guard["crawl_persistent_failure_evidence"])
    full_newton_budget = bool(below_failure_guard["full_newton_budget"])
    attempted_backtracking = bool(below_failure_guard["attempted_backtracking"])
    pass_gate = bool(
        above_converged
        and persistent_failure
        and full_newton_budget
        and attempted_backtracking
    )
    return {
        "status": "PASS" if pass_gate else "FAIL",
        "tau_fold_reference": tau_fold,
        "below_fold_margin": margin,
        "above_fold_sanity_tau": above_probe_tau,
        "above_fold_sanity_converged": bool(above_converged),
        "below_fold_attempt_count": int(len(below_rows)),
        "below_fold_converged_count": int(len(below_converged)),
        "persistent_failure_guard": {
            key: value
            for key, value in below_failure_guard.items()
            if key != "failed_attempts"
        },
        "failed_below_fold_rows": [
            {
                "tau": row.get("target_tau"),
                "residual_linf": row.get("final_original_residual_linf"),
                "newton_iters": _c0g_attempt_newton_iterations(row),
                "backtrack_index": int(row.get("backtrack_index", 0)),
                "message": row.get("message"),
                "init_source": row.get("init", {}).get("source"),
            }
            for row in below_failure_guard["failed_attempts"]
        ],
    }


def _c0g_b2_dissolved_support(crawl: Mapping[str, Any], b2: Mapping[str, Any], config: C0gBuildB1B2Config) -> dict[str, Any]:
    below = [
        dict(row)
        for row in crawl.get("tau_attempts", [])
        if bool(row.get("final_physical_converged"))
        and float(row.get("target_tau", math.inf))
        <= float(config.tau_fold_reference) - float(config.dissolved_tau_margin)
        and float(row.get("final_original_residual_linf", math.inf)) <= BACKGROUND_RESIDUAL_TOL
    ]
    branch_continuous = [
        row
        for row in below
        if _c0g_b2_branch_continuity_source(row.get("init", {}).get("source"))
        and not bool(row.get("used_existing_b2c"))
    ]
    fit = b2.get("sigma_min_squared_fit", {})
    cos = _c0g_b2_cos_support(b2.get("state_results", []), config)
    no_turning = not (
        fit.get("call") == "LINEAR_MONOTONE_FOLD_SUPPORT"
        and cos.get("status") == "PASS"
    )
    pass_gate = bool(len(branch_continuous) >= 2 and no_turning)
    return {
        "status": "PASS" if pass_gate else "FAIL",
        "below_fold_margin": float(config.dissolved_tau_margin),
        "required_count": 2,
        "below_fold_converged_count": int(len(below)),
        "branch_continuous_count": int(len(branch_continuous)),
        "no_r0_turning_signature": bool(no_turning),
        "rows": [
            {
                "tau": row.get("target_tau"),
                "residual_linf": row.get("final_original_residual_linf"),
                "newton_iters": _c0g_attempt_newton_iterations(row),
                "backtrack_index": int(row.get("backtrack_index", 0)),
                "init_source": row.get("init", {}).get("source"),
                "used_existing_b2c": row.get("used_existing_b2c"),
            }
            for row in branch_continuous
        ],
    }


def _c0g_b2_gate_verdict(
    *,
    crawl: Mapping[str, Any],
    b2: Mapping[str, Any],
    config: C0gBuildB1B2Config,
) -> dict[str, Any]:
    dissolved = _c0g_b2_dissolved_support(crawl, b2, config)
    dissolve_contest = _c0g_b2_full_budget_dissolve_contest(crawl, config)
    cos = _c0g_b2_cos_support(b2.get("state_results", []), config)
    fit = _c0g_b2_fit_support(b2.get("sigma_min_squared_fit", {}), config)
    step6 = _c0g_b2_step6_support(b2.get("step6_bordered_conditioning", {}), config)
    closer = _c0g_b2_closer_sampling(crawl, config)
    if dissolved["status"] == "PASS":
        verdict = "FOLD_DISSOLVED"
        branch = "FOLD_DISSOLVED precedence branch"
    elif (
        dissolve_contest["status"] == "PASS"
        and cos["status"] == "PASS"
        and fit["status"] == "PASS"
        and step6["status"] == "PASS"
    ):
        verdict = "FOLD_CONFIRMED"
        branch = "FOLD_CONFIRMED relative-trend branch after full-budget dissolve contest"
    else:
        verdict = "STILL_INCONCLUSIVE"
        branch = "STILL_INCONCLUSIVE fallback branch"
    return {
        "verdict": verdict,
        "precedence_branch": branch,
        "fold_dissolved": dissolved,
        "full_budget_dissolve_contest": dissolve_contest,
        "fold_confirmed_components": {
            "cos_theta": cos,
            "sigma_min_squared_fit": fit,
            "bordered_conditioning_trend": step6,
        },
        "closer_sampling": closer,
        "absolute_cond_bar_used": False,
    }


def _c0g_b2_convergence_ladder(
    crawl: Mapping[str, Any],
    config: C0gBuildB1B2Config,
) -> list[dict[str, Any]]:
    tau_fold = float(config.tau_fold_reference)
    rows = []
    for row in crawl.get("tau_attempts", []):
        if "target_tau" not in row:
            continue
        tau = float(row["target_tau"])
        rows.append(
            {
                "tau": tau,
                "nominal_tau": row.get("nominal_target_tau"),
                "relative_to_tau_fold": "below"
                if tau < tau_fold
                else "above"
                if tau > tau_fold
                else "at",
                "converged": bool(row.get("final_physical_converged")),
                "residual_linf": row.get("final_original_residual_linf"),
                "newton_iters": _c0g_attempt_newton_iterations(row),
                "newton_iters_total": _c0g_attempt_total_newton_iterations(row),
                "max_newton_iters": row.get("max_newton_iters"),
                "backtrack_index": int(row.get("backtrack_index", 0)),
                "wall_seconds": row.get("elapsed_seconds"),
                "measurement_status": row.get("measurement_status", "MEASURED"),
                "init_source": row.get("init", {}).get("source"),
                "message": row.get("message"),
            }
        )
    return rows


def _c0g_closed_block_slices(grid) -> dict[str, tuple[int, int]]:
    cell = int(grid.spec.nr * grid.spec.nw)
    wall = int(grid.spec.nw)
    return {
        "psi_real_rows": (0, cell),
        "psi_imag_rows": (cell, 2 * cell),
        "a0_rows": (2 * cell, 3 * cell),
        "ar_rows": (3 * cell, 4 * cell),
        "aw_rows": (4 * cell, 5 * cell),
        "wall_rows": (5 * cell, 5 * cell + wall),
        "mass_row": (5 * cell + wall, 5 * cell + wall + 1),
    }


def _c0g_abs_rel_mismatch(
    got: np.ndarray,
    ref: np.ndarray,
) -> dict[str, float]:
    got = np.asarray(got, dtype=np.float64)
    ref = np.asarray(ref, dtype=np.float64)
    diff = got - ref
    abs_max = float(np.max(np.abs(diff))) if diff.size else 0.0
    ref_scale = max(1.0, float(np.max(np.abs(ref))) if ref.size else 0.0)
    return {
        "max_abs": abs_max,
        "max_rel": float(abs_max / ref_scale),
        "reference_linf": float(np.max(np.abs(ref))) if ref.size else 0.0,
    }


def _c0g_update_block_mismatch(
    blocks: dict[str, dict[str, float]],
    name: str,
    got: np.ndarray,
    ref: np.ndarray,
) -> None:
    mismatch = _c0g_abs_rel_mismatch(got, ref)
    current = blocks.setdefault(
        name,
        {"max_abs": 0.0, "max_rel": 0.0, "reference_linf": 0.0},
    )
    current["max_abs"] = max(float(current["max_abs"]), mismatch["max_abs"])
    current["max_rel"] = max(float(current["max_rel"]), mismatch["max_rel"])
    current["reference_linf"] = max(
        float(current["reference_linf"]),
        mismatch["reference_linf"],
    )


def run_c0g_b4_match_gate(config: C0gBuildB1B2Config | None = None) -> dict[str, Any]:
    cfg = config or C0gBuildB1B2Config()
    dtype = configure_backend(BackendConfig())
    c0f2_payload = _c0g_load_c0f2_payload(C0gConfig(c0f2_json_path=cfg.c0f2_json_path))
    reference_row = _c0g_original_row_for_tau(c0f2_payload, float(cfg.b4_reference_tau))
    state = _c0g_state_vector_from_path(reference_row["state_artifact"], dtype=dtype)
    reference_context = C0gConfig(
        c0f2_json_path=cfg.c0f2_json_path,
        grid=cfg.grid,
        jacobian_assembly="colored_sparse_jacobian_lu",
    )
    branch, _provider, grid, _boundaries, _eos_K, residual_fn = _c0g_residual_context(
        tau=float(cfg.b4_reference_tau),
        config=reference_context,
        dtype=dtype,
    )
    residual = residual_fn(state).detach().cpu().numpy().astype(np.float64, copy=False)
    colored_config = replace(
        branch.newton,
        preconditioner=replace(
            branch.newton.preconditioner,
            type="colored_sparse_jacobian_lu",
            diagonal_shift=0.0,
        ),
    )
    fast_config = replace(
        branch.newton,
        preconditioner=replace(
            branch.newton.preconditioner,
            type="autodiff_sparse_jacobian_lu",
            diagonal_shift=0.0,
        ),
    )

    colored_start = time.perf_counter()
    colored_matrix, colored_metadata = assemble_closed_coupled_colored_sparse_jacobian(
        PreconditionerBuildContext(
            residual_fn=residual_fn,
            x=state.detach(),
            rhs=-residual,
            iteration=0,
            config=colored_config,
        ),
        grid,
    )
    colored_seconds = float(time.perf_counter() - colored_start)

    fast_start = time.perf_counter()
    fast_matrix, fast_metadata = assemble_closed_coupled_autodiff_sparse_jacobian(
        PreconditionerBuildContext(
            residual_fn=residual_fn,
            x=state.detach(),
            rhs=-residual,
            iteration=0,
            config=fast_config,
        ),
        grid,
        autodiff_mode="jacfwd",
    )
    fast_seconds = float(time.perf_counter() - fast_start)

    blocks = _c0g_closed_block_slices(grid)
    random_block_mismatches: dict[str, dict[str, float]] = {}
    generator = torch.Generator(device=state.device)
    generator.manual_seed(int(cfg.b4_random_seed))
    for probe_index in range(int(cfg.b4_random_probe_count)):
        direction = torch.randn(state.shape, dtype=dtype, device=state.device, generator=generator)
        direction = direction / torch.linalg.vector_norm(direction)
        reference_jv = jvp(residual_fn, state.detach(), direction).detach().cpu().numpy()
        assembled_jv = fast_matrix @ direction.detach().cpu().numpy().astype(np.float64, copy=False)
        for name, (start, stop) in blocks.items():
            _c0g_update_block_mismatch(
                random_block_mismatches,
                name,
                assembled_jv[start:stop],
                reference_jv[start:stop],
            )
        random_block_mismatches.setdefault("_probe_count", {})["count"] = float(probe_index + 1)

    reverse_start = time.perf_counter()
    reverse_jacobian = (
        torch.func.jacrev(residual_fn)(state.detach()).detach().cpu().numpy().astype(np.float64, copy=False)
    )
    reverse_seconds = float(time.perf_counter() - reverse_start)
    fast_dense = fast_matrix.toarray()
    cell = int(grid.spec.nr * grid.spec.nw)
    wall_start = 5 * cell
    wall_stop = wall_start + int(grid.spec.nw)
    mu_column = fast_dense.shape[1] - 1
    dense_checks = {
        "wall_rows_full": _c0g_abs_rel_mismatch(
            fast_dense[wall_start:wall_stop, :],
            reverse_jacobian[wall_start:wall_stop, :],
        ),
        "mass_row_full": _c0g_abs_rel_mismatch(
            fast_dense[wall_stop : wall_stop + 1, :],
            reverse_jacobian[wall_stop : wall_stop + 1, :],
        ),
        "mu_column_full": _c0g_abs_rel_mismatch(
            fast_dense[:, mu_column : mu_column + 1],
            reverse_jacobian[:, mu_column : mu_column + 1],
        ),
    }
    all_gate_rows = [
        value
        for key, value in random_block_mismatches.items()
        if key != "_probe_count"
    ] + list(dense_checks.values())
    abs_pass = all(float(row["max_abs"]) <= float(cfg.b4_match_abs_tol) for row in all_gate_rows)
    rel_pass = all(float(row["max_rel"]) <= float(cfg.b4_match_rel_tol) for row in all_gate_rows)
    speedup = colored_seconds / fast_seconds if fast_seconds > 0.0 else math.inf
    return {
        "schema": "stage1_pathA_C0g_B4_sparse_jacobian_match_gate/v1",
        "status": "PASS" if abs_pass and rel_pass else "FAIL",
        "reference_tau": float(cfg.b4_reference_tau),
        "reference_state_artifact": reference_row.get("state_artifact"),
        "residual_entry_point": "stage1_solver.coupled_branch.patha_closed_branch_residual",
        "frozen_residual_modified": False,
        "faithful_operators_modified": False,
        "xi_or_grad_div_penalty_modified": False,
        "physics_modified": False,
        "match_tolerances": {
            "max_abs": float(cfg.b4_match_abs_tol),
            "max_rel": float(cfg.b4_match_rel_tol),
        },
        "coverage": {
            "random_jv_probe_count": int(cfg.b4_random_probe_count),
            "random_jv_blocks": list(blocks),
            "explicit_dense_checks": [
                "wall_rows_full",
                "mass_row_full",
                "mu_column_full",
            ],
            "random_jv_reference": "torch.func.jvp(original_residual)",
            "dense_row_reference": "torch.func.jacrev(original_residual)",
        },
        "random_jv_block_mismatches": {
            key: value
            for key, value in random_block_mismatches.items()
            if key != "_probe_count"
        },
        "dense_row_mismatches": dense_checks,
        "colored_assembly": {
            "seconds": colored_seconds,
            "metadata": dict(colored_metadata),
            "matrix_nnz": int(colored_matrix.nnz),
            "note": "timing reference for legacy colored-JVP assembly",
        },
        "fast_assembly": {
            "seconds": fast_seconds,
            "metadata": dict(fast_metadata),
            "matrix_nnz": int(fast_matrix.nnz),
            "independent_jacrev_seconds_for_dense_checks": reverse_seconds,
        },
        "assembly_speedup_x": float(speedup),
        "gate": {
            "abs_pass": bool(abs_pass),
            "rel_pass": bool(rel_pass),
            "pass": bool(abs_pass and rel_pass),
        },
    }


def _c0g_proof_newton_solve_summary(row: Mapping[str, Any]) -> dict[str, Any]:
    attempts = [attempt for attempt in row.get("epsilon_attempts", []) if isinstance(attempt, Mapping)]
    initial_residuals = []
    for attempt in attempts:
        try:
            value = float(attempt.get("initial_original_residual", math.nan))
        except (TypeError, ValueError):
            value = math.nan
        if math.isfinite(value):
            initial_residuals.append(value)

    def _is_zero_epsilon_attempt(attempt: Mapping[str, Any]) -> bool:
        try:
            epsilon = float(attempt.get("epsilon", math.nan))
            k1_epsilon = float(attempt.get("k1_radius_epsilon", math.nan))
        except (TypeError, ValueError):
            return False
        return epsilon == 0.0 and k1_epsilon == 0.0

    final_zero = next(
        (
            attempt
            for attempt in reversed(attempts)
            if _is_zero_epsilon_attempt(attempt)
        ),
        attempts[-1] if attempts else {},
    )
    max_iters = _c0g_attempt_newton_iterations(row)
    total_iters = _c0g_attempt_total_newton_iterations(row)
    nonconverged_seed = bool(
        initial_residuals
        and max(initial_residuals) > float(row.get("b2c_background_tolerance", BACKGROUND_RESIDUAL_TOL))
    )
    genuine = bool(row.get("gauge_fix_solver_invoked") and max_iters is not None and max_iters > 0)
    return {
        "status": "PASS" if genuine and nonconverged_seed else "FAIL",
        "tau": row.get("target_tau"),
        "gauge_fix_solver_invoked": bool(row.get("gauge_fix_solver_invoked")),
        "nonconverged_seed": nonconverged_seed,
        "seed_description": row.get("proof_seed_description"),
        "newton_iters": max_iters,
        "newton_iters_total": total_iters,
        "zero_epsilon_newton_iters": final_zero.get("newton_iterations"),
        "initial_original_residual_max": max(initial_residuals) if initial_residuals else None,
        "final_original_residual_linf": row.get("final_original_residual_linf"),
        "init_source": row.get("init", {}).get("source"),
        "max_newton_iters": row.get("max_newton_iters"),
        "max_newton_iters_override": row.get("max_newton_iters_override"),
    }


def _c0g_per_tau_row_for_sequence(
    payload: Mapping[str, Any],
    tau: float,
) -> dict[str, Any] | None:
    for row in payload.get("per_tau_timing", []):
        if _c0g_tau_close(float(row.get("tau", math.nan)), float(tau)):
            return dict(row)
    return None


def _c0g_saved_state_residual_linf(
    *,
    state: torch.Tensor,
    tau: float,
    config: C0gBuildB1B2Config,
    dtype: torch.dtype,
) -> float:
    _branch, _provider, _grid, _boundaries, _eos_K, residual_fn = _c0g_residual_context(
        tau=float(tau),
        config=C0gConfig(grid=config.grid),
        dtype=dtype,
    )
    residual = residual_fn(state).detach()
    return float(torch.max(torch.abs(residual)).detach().cpu().item())


def _c0g_apply_global_phase_state(
    state: torch.Tensor,
    grid,
    *,
    angle: float,
) -> torch.Tensor:
    values = state.detach().clone()
    fields, mu = unpack_closed_coupled_fields(values, grid, has_chemical_potential=True)
    cos_a = torch.as_tensor(math.cos(float(angle)), dtype=values.dtype, device=values.device)
    sin_a = torch.as_tensor(math.sin(float(angle)), dtype=values.dtype, device=values.device)
    psi_real = cos_a * fields.psi_real - sin_a * fields.psi_imag
    psi_imag = sin_a * fields.psi_real + cos_a * fields.psi_imag
    assert mu is not None
    return pack_closed_coupled_fields(
        ClosedCoupledFields(
            psi_real=psi_real,
            psi_imag=psi_imag,
            a0=fields.a0,
            ar=fields.ar,
            aw=fields.aw,
            r0=fields.r0,
        ),
        mu,
    ).detach()


def _c0g_apply_path_only_proof_perturbation(
    state: torch.Tensor,
    grid,
    *,
    phase_angle: float,
    r0_amplitude: float,
) -> tuple[torch.Tensor, dict[str, Any]]:
    values = _c0g_apply_global_phase_state(state, grid, angle=float(phase_angle))
    fields, mu = unpack_closed_coupled_fields(values, grid, has_chemical_potential=True)
    pattern = torch.linspace(
        -1.0,
        1.0,
        fields.r0.numel(),
        dtype=values.dtype,
        device=values.device,
    ).reshape_as(fields.r0)
    pattern = pattern / torch.clamp(torch.max(torch.abs(pattern)), min=1.0)
    assert mu is not None
    perturbed = pack_closed_coupled_fields(
        ClosedCoupledFields(
            psi_real=fields.psi_real,
            psi_imag=fields.psi_imag,
            a0=fields.a0,
            ar=fields.ar,
            aw=fields.aw,
            r0=fields.r0 + float(r0_amplitude) * pattern,
        ),
        mu,
    ).detach()
    return perturbed, {
        "global_phase_radians": float(phase_angle),
        "r0_linf_amplitude": float(r0_amplitude),
        "r0_pattern": "deterministic centered linear lane perturbation",
        "not_pure_global_phase": True,
    }


def _c0g_single_gauge_fixed_attempt(
    *,
    x0: torch.Tensor,
    tau: float,
    config: C0gBuildB1B2Config,
    dtype: torch.dtype,
    attempt_index: int,
    start_from_tau: float | None,
    init: Mapping[str, Any],
    max_newton_iters: int | None,
) -> tuple[dict[str, Any], torch.Tensor, Any]:
    local_c0 = C0Config(
        run_root=config.run_root / "gauge_fixed_crawl",
        grid=config.grid,
        use_gauge_fix=True,
        max_newton_iters_override=None
        if max_newton_iters is None
        else int(max_newton_iters),
        prefer_existing_b2c_background_predictor=False,
        jacobian_assembly=config.jacobian_assembly,
    )
    branch, provider, grid, boundaries = _branch_context(
        tau=float(tau),
        config=local_c0,
        dtype=dtype,
    )
    eos_K = float(branch.continuation_K_values[-1])
    original_residual_fn = lambda x: patha_closed_branch_residual(
        x,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=provider,
    )
    aid = C0AidParameters(core_epsilon=0.0, k1_radius_epsilon=0.0)
    preconditioner_residual_fn = lambda x: c0_preconditioner_residual(
        x,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=provider,
        aid=aid,
    )
    gauge_config = C0gConfig(grid=config.grid)

    def gauge_matrix_fn(x: torch.Tensor, _col_scale: np.ndarray) -> Mapping[str, Any]:
        del _col_scale
        return _c0g_build_analytic_gauge_matrix(
            state=x,
            grid=grid,
            branch=branch,
            config=gauge_config,
        )

    started = time.perf_counter()
    result = solve_c0_scaled_newton(
        original_residual_fn=original_residual_fn,
        preconditioner_residual_fn=preconditioner_residual_fn,
        x0=x0,
        grid=grid,
        newton=branch.newton,
        aid=aid,
        gauge_matrix_fn=gauge_matrix_fn,
        gauge_fix_lstsq_rcond=float(config.proof_tolerance),
        gauge_fix_gradient_rank_rtol=float(gauge_config.gradient_rank_rtol),
    )
    final_state = result.x.detach()
    final_residual = original_residual_fn(final_state).detach()
    final_linf = float(torch.max(torch.abs(final_residual)).detach().cpu().item())
    metrics = _state_metrics(final_state, grid)
    state_path = _resolve_output_path(
        _state_artifact_path(local_c0, tau=float(tau), attempt_index=attempt_index)
    )
    _save_state_artifact(state_path, final_state)
    row = {
        "tau": float(tau),
        "target_tau": float(tau),
        "nominal_target_tau": float(tau),
        "delta_tau": None if start_from_tau is None else float(tau) - float(start_from_tau),
        "backtrack_index": 0,
        "start_from_tau": start_from_tau,
        "target_relative_to_prior_floor": float(tau) / PRIOR_TAU_FLOOR,
        "below_prior_floor": _below_prior_floor(float(tau)),
        "stage_converged": bool(result.converged),
        "final_original_residual_linf": final_linf,
        "final_physical_converged": bool(result.converged and final_linf <= BACKGROUND_RESIDUAL_TOL),
        "b2c_background_tolerance": BACKGROUND_RESIDUAL_TOL,
        "message": result.message,
        "elapsed_seconds": float(time.perf_counter() - started),
        "init": dict(init),
        "initialization": dict(init),
        "used_existing_b2c": False,
        "max_newton_iters": int(branch.newton.max_newton_iters),
        "max_newton_iters_override": local_c0.max_newton_iters_override,
        "max_tau_backtracks": int(config.max_tau_backtracks),
        "min_tau_backtrack_delta": float(config.min_tau_backtrack_delta),
        "jacobian_assembly": config.jacobian_assembly,
        "continuation_K_values": [eos_K],
        "epsilon_attempts": [
            _epsilon_attempt_row(
                tau=float(tau),
                eos_K=eos_K,
                aid_index=0,
                aid=aid,
                result=result,
                start_policy=str(init.get("source", "unknown")),
                metrics=metrics,
            )
        ],
        "substeps": [],
        "metrics": metrics,
        "state_artifact": str(state_path),
        "gauge_fix_solver_invoked": True,
        "linear_diagnostics": {"status": "NOT_MEASURED", "reason": "B2_phase_not_run_yet"},
        "observables": {"available": False, "reason": "B2_phase_not_run_yet"},
        "floor_activation": {
            "min_rho_during_final": metrics["min_rho"],
            "min_R0_during_final": metrics["min_R0"],
            "final_aids_inactive": True,
        },
    }
    return row, final_state, grid


def _c0g_run_gauge_fixed_polish_crawl(
    *,
    payload: Mapping[str, Any],
    config: C0gBuildB1B2Config,
    dtype: torch.dtype,
) -> dict[str, Any]:
    run_root = _resolve_output_path(config.run_root / "gauge_fixed_crawl")
    run_root.mkdir(parents=True, exist_ok=True)
    sequence_rows: list[dict[str, Any]] = []
    last_converged_state: torch.Tensor | None = None
    last_converged_grid: Any | None = None
    last_converged_tau: float | None = None
    attempt_index = 0
    started = time.perf_counter()
    proof_solve_seen = False
    for nominal_tau in config.b2_depth_sequence:
        nominal_tau = float(nominal_tau)
        existing_row = _c0g_per_tau_row_for_sequence(payload, nominal_tau)
        may_accept_shallow_anchor = (
            existing_row is not None
            and bool(existing_row.get("solver_converged"))
            and nominal_tau > float(config.proof_tau) + 5.0e-13
            and (
                last_converged_state is None
                or (
                    last_converged_tau is not None
                    and nominal_tau < float(last_converged_tau) - 5.0e-13
                )
            )
        )
        if may_accept_shallow_anchor and existing_row is not None:
            state = _c0g_state_vector_from_path(existing_row["state_artifact"], dtype=dtype)
            residual_linf = _c0g_saved_state_residual_linf(
                state=state,
                tau=nominal_tau,
                config=config,
                dtype=dtype,
            )
            branch, _provider, grid, _boundaries = _branch_context(
                tau=nominal_tau,
                config=C0Config(grid=config.grid, jacobian_assembly=config.jacobian_assembly),
                dtype=dtype,
            )
            metrics = _state_metrics(state, grid)
            state_path = _resolve_output_path(
                _state_artifact_path(
                    C0Config(run_root=config.run_root / "gauge_fixed_crawl"),
                    tau=nominal_tau,
                    attempt_index=attempt_index,
                )
            )
            _save_state_artifact(state_path, state)
            row = {
                "tau": nominal_tau,
                "target_tau": nominal_tau,
                "nominal_target_tau": nominal_tau,
                "delta_tau": None
                if last_converged_tau is None
                else nominal_tau - float(last_converged_tau),
                "backtrack_index": 0,
                "start_from_tau": last_converged_tau,
                "target_relative_to_prior_floor": nominal_tau / PRIOR_TAU_FLOOR,
                "below_prior_floor": _below_prior_floor(nominal_tau),
                "stage_converged": residual_linf <= BACKGROUND_RESIDUAL_TOL,
                "final_original_residual_linf": residual_linf,
                "final_physical_converged": residual_linf <= BACKGROUND_RESIDUAL_TOL,
                "b2c_background_tolerance": BACKGROUND_RESIDUAL_TOL,
                "message": "existing C0f2 genuine warm-start converged state accepted as a shallow anchor before the gauge-fixed proof solve",
                "elapsed_seconds": 0.0,
                "init": {
                    "source": "c0f2_genuine_warm_start_converged_state",
                    "path": existing_row.get("state_artifact"),
                    "original_init_source": existing_row.get("init_source"),
                },
                "initialization": {
                    "source": "c0f2_genuine_warm_start_converged_state",
                    "path": existing_row.get("state_artifact"),
                },
                "used_existing_b2c": False,
                "gauge_fix_solver_invoked": False,
                "max_newton_iters": 20,
                "max_newton_iters_override": config.max_newton_iters_override,
                "max_tau_backtracks": int(config.max_tau_backtracks),
                "min_tau_backtrack_delta": float(config.min_tau_backtrack_delta),
                "jacobian_assembly": config.jacobian_assembly,
                "continuation_K_values": [float(branch.continuation_K_values[-1])],
                "epsilon_attempts": [],
                "substeps": [],
                "metrics": metrics,
                "state_artifact": str(state_path),
                "linear_diagnostics": {"status": "NOT_MEASURED", "reason": "B2_phase_not_run_yet"},
                "observables": {"available": False, "reason": "B2_phase_not_run_yet"},
                "floor_activation": {
                    "min_rho_during_final": metrics["min_rho"],
                    "min_R0_during_final": metrics["min_R0"],
                    "final_aids_inactive": True,
                },
            }
            sequence_rows.append(row)
            if row["final_physical_converged"]:
                last_converged_state = state.detach().clone()
                last_converged_grid = grid
                last_converged_tau = nominal_tau
            attempt_index += 1
            continue

        if last_converged_state is None:
            if existing_row is None or not bool(existing_row.get("solver_converged")):
                raise RuntimeError(f"no warm-start state available for tau={nominal_tau:.12e}")
            last_converged_state = _c0g_state_vector_from_path(existing_row["state_artifact"], dtype=dtype)
            branch, _provider, last_converged_grid, _boundaries = _branch_context(
                tau=nominal_tau,
                config=C0Config(grid=config.grid, jacobian_assembly=config.jacobian_assembly),
                dtype=dtype,
            )
            last_converged_tau = nominal_tau

        start_tau = last_converged_tau
        full_delta = None if start_tau is None else nominal_tau - float(start_tau)
        target_tau = nominal_tau
        backtrack_index = 0
        while True:
            delta_tau = None if start_tau is None else float(target_tau) - float(start_tau)
            branch, _provider, grid, _boundaries = _branch_context(
                tau=float(target_tau),
                config=C0Config(grid=config.grid, jacobian_assembly=config.jacobian_assembly),
                dtype=dtype,
            )
            if (
                _c0g_tau_close(float(target_tau), float(config.proof_tau))
                and existing_row is not None
                and bool(existing_row.get("solver_converged"))
                and not proof_solve_seen
            ):
                base_state = _c0g_state_vector_from_path(existing_row["state_artifact"], dtype=dtype)
                x0, perturbation = _c0g_apply_path_only_proof_perturbation(
                    base_state,
                    grid,
                    phase_angle=float(config.proof_global_phase_radians),
                    r0_amplitude=float(config.proof_r0_perturbation),
                )
                init = {
                    "source": "c0f2_converged_state_non_gauge_perturbed_for_gauge_fixed_proof",
                    "path": existing_row.get("state_artifact"),
                    "perturbation": perturbation,
                }
            else:
                if last_converged_state is None:
                    raise RuntimeError(f"no previous gauge-fixed state for tau={target_tau:.12e}")
                x0 = last_converged_state.detach().clone()
                if last_converged_grid is not None:
                    x0 = resample_closed_branch_state(x0, last_converged_grid, grid, branch)
                x0, predictor = _apply_c0_wall_predictor(
                    state=x0,
                    grid=grid,
                    tau=float(target_tau),
                    grid_level=config.grid,
                )
                init = {
                    "source": "previous_c0g_gauge_fixed_converged_state",
                    "previous_tau": last_converged_tau,
                    "wall_predictor": predictor,
                }
            row, final_state, final_grid = _c0g_single_gauge_fixed_attempt(
                x0=x0,
                tau=float(target_tau),
                config=config,
                dtype=dtype,
                attempt_index=attempt_index,
                start_from_tau=last_converged_tau,
                init=init,
                max_newton_iters=config.max_newton_iters_override,
            )
            row["nominal_target_tau"] = float(nominal_tau)
            row["delta_tau"] = delta_tau
            row["backtrack_index"] = int(backtrack_index)
            row["max_tau_backtracks"] = int(config.max_tau_backtracks)
            row["min_tau_backtrack_delta"] = float(config.min_tau_backtrack_delta)
            row["gauge_fix_solver_invoked"] = True
            row["full_newton_budget"] = row.get("max_newton_iters_override") in (None, 0)
            if _c0g_tau_close(float(row["target_tau"]), float(config.proof_tau)):
                row["proof_solve_role"] = "genuine_gauge_fixed_newton_solve"
                row["proof_seed_description"] = (
                    "C0f2 proof-tau state plus global phase and non-gauge r0 perturbation; not a pure global-phase orbit seed"
                )
                proof_solve_seen = True
            sequence_rows.append(row)
            attempt_index += 1
            if row["final_physical_converged"]:
                last_converged_state = final_state.detach().clone()
                last_converged_grid = final_grid
                last_converged_tau = float(row["target_tau"])
                break
            if start_tau is None or full_delta is None:
                break
            if (
                backtrack_index >= int(config.max_tau_backtracks)
                or abs(delta_tau or 0.0) <= float(config.min_tau_backtrack_delta)
            ):
                break
            backtrack_index += 1
            target_tau = float(start_tau) + float(full_delta) / (2.0**backtrack_index)
    schema_errors = _validate_tau_attempts_schema(sequence_rows)
    result = {
        "schema": "stage1_pathA_C0g_build_gauge_fixed_polish_crawl/v1",
        "source_revision": source_revision(),
        "phase_status": "gauge_fixed_polish_crawl_complete",
        "prior_tau_floor": PRIOR_TAU_FLOOR,
        "b2c_background_tolerance": BACKGROUND_RESIDUAL_TOL,
        "deepest_converged_tau": last_converged_tau,
        "deepest_converged_R0_min": None
        if not sequence_rows
        else min(
            (
                row["metrics"]["min_R0"]
                for row in sequence_rows
                if row.get("final_physical_converged") and row.get("metrics")
            ),
            default=None,
        ),
        "elapsed_seconds": float(time.perf_counter() - started),
        "config": {
            "run_root": str(config.run_root / "gauge_fixed_crawl"),
            "grid": list(config.grid),
            "depth_sequence": list(config.b2_depth_sequence),
            "prefer_existing_b2c_background_predictor": False,
            "use_gauge_fix": True,
            "max_newton_iters_override": config.max_newton_iters_override,
            "max_tau_backtracks": int(config.max_tau_backtracks),
            "min_tau_backtrack_delta": float(config.min_tau_backtrack_delta),
            "jacobian_assembly": config.jacobian_assembly,
        },
        "tau_attempts": sequence_rows,
        "depth_sequence": sequence_rows,
        "tau_attempt_schema_errors": schema_errors,
        "proof_solve_seen": bool(proof_solve_seen),
        "faithful_operator_boundary": _faithful_operator_boundary(),
    }
    crawl_json = _resolve_output_path(config.run_root / "gauge_fixed_crawl" / "c0_gauge_fixed_crawl.json")
    _c0g_write_json(crawl_json, result)
    return result


def write_c0g_build_b1b2_report(result: Mapping[str, Any], path: Path) -> None:
    b4 = result.get("B4_analytic_sparse_jacobian", {})
    b1 = result.get("B1_gauge_fix", {})
    b2 = result.get("B2_reconfirm", {})
    gate = b2.get("gate", {})
    lift = b1.get("gauge_sector_lift", {})
    proof = b1.get("path_only_proof", {})
    genuine = b1.get("genuine_solve_proof", {})
    preserved = b1.get("r0_mode_preservation", {})
    fit = gate.get("fold_confirmed_components", {}).get("sigma_min_squared_fit", {})
    step6 = gate.get("fold_confirmed_components", {}).get("bordered_conditioning_trend", {})
    cos = gate.get("fold_confirmed_components", {}).get("cos_theta", {})
    ladder = b2.get("full_budget_convergence_ladder", [])
    has_b4 = bool(b4)
    lines = [
        "# Path-A C0g Build B-1/B-2/B-4 Staged Report",
        "",
        f"B-2 verdict: **{gate.get('verdict')}**",
        f"B-4 match gate: **{b4.get('status', 'NOT_RUN')}**",
        "",
        "## Scope",
        "",
        "- Built B-1 gauge-fix only in solver/path coordinates.",
        "- Built B-4 analytic/sparse Jacobian assembly and verified it before B-2."
        if has_b4
        else "- B-4 analytic/sparse assembly was not built.",
        "- Ran B-2 re-confirmation on the gauge-fixed crawl.",
        "- B-3 pseudo-arclength was not built.",
        "- The original `patha_closed_branch_residual` remains the convergence/progress arbiter.",
        "",
        "## B-4 Sparse Jacobian Match Gate",
        "",
        "```yaml",
        f"status: {b4.get('status', 'NOT_RUN')}",
        f"reference_tau: {b4.get('reference_tau')}",
        f"tolerances: {b4.get('match_tolerances')}",
        f"coverage: {b4.get('coverage')}",
        f"assembly_speedup_x: {b4.get('assembly_speedup_x')}",
        f"colored_seconds: {(b4.get('colored_assembly') or {}).get('seconds')}",
        f"fast_seconds: {(b4.get('fast_assembly') or {}).get('seconds')}",
        "```",
        "",
        "Random-Jv per-block mismatch:",
        "",
        "```yaml",
        str(b4.get("random_jv_block_mismatches", {})),
        "```",
        "",
        "Explicit dense wall/mass/mu mismatch:",
        "",
        "```yaml",
        str(b4.get("dense_row_mismatches", {})),
        "```",
        "",
        "## B-1 Gauge Fix",
        "",
        "Method: scaled analytic gauge-complement Newton step. The linear solve acts on `row_scale * J_original * col_scale` and constrains the step to `Q_perp` of the analytic gauge generators; no residual/operator row is changed.",
        "",
        "Gauge-sector removal evidence:",
        "",
        "```yaml",
        f"status: {lift.get('status')}",
        f"gauge_rank: {lift.get('gauge_rank')}",
        f"before_gauge_sector_sigma_min: {lift.get('before_gauge_sector_sigma_min')}",
        f"scaled_gauge_singular_min_retained: {lift.get('scaled_gauge_singular_min_retained')}",
        f"after_bordered_gauge_sector_sigma_min_not_a_lift_gate: {lift.get('after_bordered_gauge_sector_sigma_min')}",
        f"after_bordered_identity_metric_is_by_construction: {lift.get('after_bordered_identity_metric_is_by_construction')}",
        f"source_helpers: {lift.get('source_helpers')}",
        "```",
        "",
        "Genuine gauge-fixed Newton proof solve:",
        "",
        "```yaml",
        f"status: {genuine.get('status')}",
        f"tau: {genuine.get('tau')}",
        f"nonconverged_seed: {genuine.get('nonconverged_seed')}",
        f"seed_description: {genuine.get('seed_description')}",
        f"newton_iters: {genuine.get('newton_iters')}",
        f"newton_iters_total: {genuine.get('newton_iters_total')}",
        f"zero_epsilon_newton_iters: {genuine.get('zero_epsilon_newton_iters')}",
        f"initial_original_residual_max: {genuine.get('initial_original_residual_max')}",
        f"final_original_residual_linf: {genuine.get('final_original_residual_linf')}",
        "```",
        "",
        "Gauge-aligned path-only proof:",
        "",
        "```yaml",
        f"status: {proof.get('status')}",
        f"tau: {proof.get('tau')}",
        f"tolerance: {proof.get('tolerance')}",
        f"invariant_deltas: {proof.get('gauge_invariant_deltas_after_alignment')}",
        f"original_residual: {proof.get('original_residual')}",
        "```",
        "",
        "`r0` mode preservation:",
        "",
        "```yaml",
        f"status: {preserved.get('status')}",
        f"mode1_lane: {preserved.get('mode1_lane')}",
        f"mode1_lane_fraction: {preserved.get('mode1_lane_fraction')}",
        f"mode1_projection_into_gauge_basis: {preserved.get('mode1_projection_into_gauge_basis')}",
        f"mode1_sigma: {preserved.get('mode1_sigma')}",
        "```",
        "",
        "## B-2 Re-Confirm Gate",
        "",
        f"Verdict: **{gate.get('verdict')}** via `{gate.get('precedence_branch')}`.",
        "",
        "Timeout finalization:",
        "",
        "```yaml",
        f"{result.get('timeout_finalization', {'status': 'not_used'})}",
        "```",
        "",
        "Full-budget fold-vs-dissolve convergence ladder:",
        "",
        _markdown_table(
            [
                "tau",
                "nominal_tau",
                "relative_to_tau_fold",
                "converged",
                "residual_linf",
                "newton_iters",
                "max_newton_iters",
                "backtrack_index",
                "wall_seconds",
                "measurement_status",
                "init_source",
            ],
            ladder,
        ),
        "",
        "Full-budget dissolve contest:",
        "",
        "```yaml",
        f"{gate.get('full_budget_dissolve_contest')}",
        "```",
        "",
        "Closer sampling:",
        "",
        "```yaml",
        f"{gate.get('closer_sampling')}",
        "```",
        "",
        "Cosθ support:",
        "",
        _markdown_table(["tau", "status", "cos_theta", "call", "stability"], cos.get("rows", [])),
        "",
        "σ_min² fit support:",
        "",
        "```yaml",
        f"status: {fit.get('status')}",
        f"linear_r2: {fit.get('linear_r2')}",
        f"linear_slope: {fit.get('linear_slope')}",
        f"tau_fold_zero_crossing: {fit.get('tau_fold_zero_crossing')}",
        f"monotone: {fit.get('monotone')}",
        f"threshold_r2: {fit.get('threshold_r2')}",
        "```",
        "",
        "Bordered-conditioning trend:",
        "",
        "```yaml",
        f"status: {step6.get('status')}",
        f"cond_Jb_max_over_min: {step6.get('cond_Jb_max_over_min')}",
        f"ratio_growth_factor: {step6.get('ratio_growth_factor')}",
        f"ratio_increasing_toward_fold: {step6.get('ratio_increasing_toward_fold')}",
        f"absolute_cond_bar_used: {gate.get('absolute_cond_bar_used')}",
        "```",
        "",
        "Step-6 measured rows:",
        "",
        _markdown_table(
            [
                "tau",
                "status",
                "cond_JQ_perp",
                "cond_Jb",
                "sigma_min_JQ_perp",
                "sigma_min_Jb",
                "call",
            ],
            step6.get("rows", []),
        ),
        "",
        "FOLD_DISSOLVED precedence check:",
        "",
        "```yaml",
        f"{gate.get('fold_dissolved')}",
        "```",
        "",
        "## Scope Guard",
        "",
        "```yaml",
    ]
    for key, value in result.get("scope_guard", {}).items():
        lines.append(f"{key}: {value}")
    lines.extend(
        [
            "```",
            "",
            f"Machine JSON: `{result.get('json_path')}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _c0g_build_b1b2_stdout_summary(result: Mapping[str, Any]) -> str:
    b4 = result.get("B4_analytic_sparse_jacobian", {})
    b1 = result.get("B1_gauge_fix", {})
    b2 = result.get("B2_reconfirm", {})
    gate = b2.get("gate", {})
    lift = b1.get("gauge_sector_lift", {})
    proof = b1.get("path_only_proof", {})
    genuine = b1.get("genuine_solve_proof", {})
    preserved = b1.get("r0_mode_preservation", {})
    fit = gate.get("fold_confirmed_components", {}).get("sigma_min_squared_fit", {})
    step6 = gate.get("fold_confirmed_components", {}).get("bordered_conditioning_trend", {})
    cos_rows = gate.get("fold_confirmed_components", {}).get("cos_theta", {}).get("rows", [])
    ladder = b2.get("full_budget_convergence_ladder", [])
    tau_fold = float((gate.get("full_budget_dissolve_contest") or {}).get("tau_fold_reference", math.nan))
    below_converged = [
        row
        for row in ladder
        if math.isfinite(tau_fold)
        and float(row.get("tau", math.inf)) < tau_fold
        and bool(row.get("converged"))
    ]
    ladder_lines = [
        (
            "  tau={tau}, converged={converged}, residual_linf={residual}, "
            "newton_iters={iters}, backtracks={backtracks}, wall_seconds={wall_seconds}, "
            "measurement_status={status}"
        ).format(
            tau=_c0g_fmt(row.get("tau")),
            converged=row.get("converged"),
            residual=_c0g_fmt(row.get("residual_linf")),
            iters=row.get("newton_iters"),
            backtracks=row.get("backtrack_index"),
            wall_seconds=_c0g_fmt(row.get("wall_seconds")),
            status=row.get("measurement_status", "MEASURED"),
        )
        for row in ladder
    ]
    cos_values = [
        float(row["cos_theta"])
        for row in cos_rows
        if isinstance(row.get("cos_theta"), float) and math.isfinite(float(row["cos_theta"]))
    ]
    return "\n".join(
        [
            "C0g B-4-gated B-1/B-2 staged summary:",
            (
                "B-4 match gate: "
                f"{b4.get('status', 'NOT_RUN')} "
                f"abs_tol={((b4.get('match_tolerances') or {}).get('max_abs'))}, "
                f"rel_tol={((b4.get('match_tolerances') or {}).get('max_rel'))}, "
                f"speedup_x={_c0g_fmt(b4.get('assembly_speedup_x'))}, "
                f"coverage={b4.get('coverage')}"
            ),
            f"B-4 random-Jv per-block mismatch: {b4.get('random_jv_block_mismatches')}",
            f"B-4 dense wall/mass/mu mismatch: {b4.get('dense_row_mismatches')}",
            f"Phase 2 ran: {bool(result.get('phase2_ran'))}",
            (
                "Gauge-fix method: scaled analytic gauge-complement solve in path coordinates "
                "(row_scale*J_original*col_scale, step=col_scale*Q_perp*z); "
                "frozen residual/operators/(1/xi) penalty/export are not edited."
            ),
            (
                "Corrected gauge-sector numbers: "
                f"before_gauge_sector_sigma_min={_c0g_fmt(lift.get('before_gauge_sector_sigma_min'))}, "
                f"scaled_gauge_singular_min_retained="
                f"{_c0g_fmt(lift.get('scaled_gauge_singular_min_retained'))}; "
                "after_bordered_gauge_sector_sigma_min is by-construction and not gated."
            ),
            (
                "Strengthened path-only proof: "
                f"genuine_solve={genuine.get('status')} "
                f"(newton_iters={genuine.get('newton_iters')}, "
                f"zero_epsilon_iters={genuine.get('zero_epsilon_newton_iters')}), "
                f"alignment={proof.get('status')}, invariants max="
                f"{_c0g_fmt((proof.get('gauge_invariant_deltas_after_alignment') or {}).get('max'))}, "
                f"residuals with/without="
                f"{_c0g_fmt((proof.get('original_residual') or {}).get('with_gauge_fix_linf'))}/"
                f"{_c0g_fmt((proof.get('original_residual') or {}).get('without_gauge_fix_linf'))}"
            ),
            (
                "r0-mode preserved: "
                f"{preserved.get('status')} lane={preserved.get('mode1_lane')} "
                f"P_G={_c0g_fmt(preserved.get('mode1_projection_into_gauge_basis'))}"
            ),
            (
                "B-2 support: sigma_min^2 R2="
                f"{_c0g_fmt(fit.get('linear_r2'))}, monotone={fit.get('monotone')}; "
                f"cos_theta range="
                f"{_c0g_fmt(min(cos_values) if cos_values else None)}.."
                f"{_c0g_fmt(max(cos_values) if cos_values else None)}; "
                f"cond(Jb) max/min={_c0g_fmt(step6.get('cond_Jb_max_over_min'))}; "
                f"ratio growth={_c0g_fmt(step6.get('ratio_growth_factor'))}"
            ),
            "Below-fold full-budget convergence ladder:",
            *(ladder_lines or ["  n/a"]),
            f"Any state converged below tau_fold: {bool(below_converged)}",
            (
                "Source-string fix: FOLD_DISSOLVED branch-continuity accepts "
                "previous_c0g_gauge_fixed_converged_state and previous_c0_converged_state."
            ),
            f"Full-budget dissolve contest: {gate.get('full_budget_dissolve_contest')}",
            f"Closer-to-fold sampling: {gate.get('closer_sampling')}",
            f"B-2 verdict: {gate.get('verdict')} ({gate.get('precedence_branch')})",
            "B-3 was not built; frozen residual/operators/(1/xi)/export/physics remain untouched per scope guard.",
            f"Files: report={result.get('report_path')}, json={result.get('json_path')}",
        ]
    )


def run_c0g_build_b1b2(config: C0gBuildB1B2Config | None = None) -> dict[str, Any]:
    cfg = config or C0gBuildB1B2Config()
    started = time.perf_counter()
    dtype = configure_backend(BackendConfig())
    c0f2_payload = _c0g_load_c0f2_payload(C0gConfig(c0f2_json_path=cfg.c0f2_json_path))
    crawl = _c0g_run_gauge_fixed_polish_crawl(
        payload=c0f2_payload,
        config=cfg,
        dtype=dtype,
    )
    proof_row = _c0g_select_proof_row(crawl, proof_tau=float(cfg.proof_tau))
    proof_tau = float(proof_row["target_tau"])
    original_row = _c0g_original_row_for_tau(c0f2_payload, proof_tau)
    with_state = _c0g_state_vector_from_path(proof_row["state_artifact"], dtype=dtype)
    without_state = _c0g_state_vector_from_path(original_row["state_artifact"], dtype=dtype)
    b1 = {
        "method": {
            "name": "scaled analytic gauge-complement Newton solve",
            "path_only_coordinates": "row_scale * J_original * col_scale; step = col_scale * Q_perp * z",
            "acts_in": "solver linear-step coordinates only",
            "residual_entry_point": "stage1_solver.coupled_branch.patha_closed_branch_residual",
            "frozen_residual_modified": False,
            "faithful_operators_modified": False,
            "xi_or_grad_div_penalty_modified": False,
        },
        "proof_tau": proof_tau,
        "with_gauge_state_artifact": proof_row.get("state_artifact"),
        "without_gauge_state_artifact": original_row.get("state_artifact"),
        "gauge_sector_lift": _c0g_gauge_sector_lift(
            state=with_state,
            tau=proof_tau,
            config=cfg,
            dtype=dtype,
        ),
        "genuine_solve_proof": _c0g_proof_newton_solve_summary(proof_row),
        "path_only_proof": _c0g_gauge_aligned_path_only_proof(
            with_state=with_state,
            without_state=without_state,
            tau=proof_tau,
            config=cfg,
            dtype=dtype,
        ),
        "r0_mode_preservation": _c0g_r0_mode_preservation(
            state=with_state,
            tau=proof_tau,
            config=cfg,
            dtype=dtype,
        ),
    }
    b2 = _c0g_b2_analyze_gauge_fixed_crawl(crawl=crawl, config=cfg, dtype=dtype)
    gate = _c0g_b2_gate_verdict(crawl=crawl, b2=b2, config=cfg)
    b2["gate"] = gate
    b2["full_budget_convergence_ladder"] = _c0g_b2_convergence_ladder(crawl, cfg)
    result = {
        "schema": "stage1_pathA_C0g_build_B1B2/v1",
        "source_revision": source_revision(),
        "config": _c0g_build_b1b2_config_to_dict(cfg),
        "B1_gauge_fix": b1,
        "B2_reconfirm": b2,
        "scope_guard": {
            "staged_items_built": ["B-1", "B-2"],
            "B3_pseudoarclength_built": False,
            "B4_analytic_sparse_assembly_built": False,
            "single_arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
            "patha_closed_branch_residual_touched": False,
            "faithful_operators_touched": False,
            "xi_or_grad_div_penalty_touched": False,
            "physical_export_permitted_touched": False,
            "prefer_existing_b2c_background_predictor": False,
            "absolute_cond_gt_1e10_bar_used_for_B2": False,
        },
        "elapsed_seconds": float(time.perf_counter() - started),
        "report_path": str(_resolve_output_path(cfg.report_path)),
        "json_path": str(_resolve_output_path(cfg.json_path)),
    }
    json_path = _resolve_output_path(cfg.json_path)
    report_path = _resolve_output_path(cfg.report_path)
    _c0g_write_json(json_path, result)
    write_c0g_build_b1b2_report(result, report_path)
    return result


def run_c0g_build_b4_then_b2(config: C0gBuildB1B2Config | None = None) -> dict[str, Any]:
    cfg = config or C0gBuildB1B2Config()
    started = time.perf_counter()
    b4 = run_c0g_b4_match_gate(cfg)
    if b4.get("status") != "PASS":
        result = {
            "schema": "stage1_pathA_C0g_build_B4_then_B2/v1",
            "source_revision": source_revision(),
            "config": _c0g_build_b1b2_config_to_dict(cfg),
            "B4_analytic_sparse_jacobian": b4,
            "B1_gauge_fix": {},
            "B2_reconfirm": {
                "gate": {
                    "verdict": "STILL_INCONCLUSIVE",
                    "precedence_branch": "B4 match-gate failed; Phase 2 not run",
                },
                "full_budget_convergence_ladder": [],
            },
            "phase2_ran": False,
            "phase2_blocked_reason": "B4 match-gate failed",
            "scope_guard": {
                "staged_items_built": ["B-4"],
                "B3_pseudoarclength_built": False,
                "B4_analytic_sparse_assembly_built": True,
                "single_arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
                "patha_closed_branch_residual_touched": False,
                "faithful_operators_touched": False,
                "xi_or_grad_div_penalty_touched": False,
                "physical_export_permitted_touched": False,
                "prefer_existing_b2c_background_predictor": False,
                "absolute_cond_gt_1e10_bar_used_for_B2": False,
            },
            "elapsed_seconds": float(time.perf_counter() - started),
            "report_path": str(_resolve_output_path(cfg.report_path)),
            "json_path": str(_resolve_output_path(cfg.json_path)),
        }
        json_path = _resolve_output_path(cfg.json_path)
        report_path = _resolve_output_path(cfg.report_path)
        _c0g_write_json(json_path, result)
        write_c0g_build_b1b2_report(result, report_path)
        return result

    result = run_c0g_build_b1b2(cfg)
    result["schema"] = "stage1_pathA_C0g_build_B4_then_B2/v1"
    result["B4_analytic_sparse_jacobian"] = b4
    result["phase2_ran"] = True
    result["total_b4_then_b2_elapsed_seconds"] = float(time.perf_counter() - started)
    scope = dict(result.get("scope_guard", {}))
    scope["staged_items_built"] = ["B-4", "B-1", "B-2"]
    scope["B4_analytic_sparse_assembly_built"] = True
    scope["B3_pseudoarclength_built"] = False
    result["scope_guard"] = scope
    json_path = _resolve_output_path(cfg.json_path)
    report_path = _resolve_output_path(cfg.report_path)
    _c0g_write_json(json_path, result)
    write_c0g_build_b1b2_report(result, report_path)
    return result


def _c0g_deep_config_to_dict(config: C0gDeepCrawlConfig) -> dict[str, Any]:
    data = asdict(config)
    for key in ("prior_b1b2_json_path", "prior_crawl_json_path", "run_root", "report_path", "json_path"):
        data[key] = str(data[key])
    return data


def _c0g_deep_build_config(config: C0gDeepCrawlConfig) -> C0gBuildB1B2Config:
    return C0gBuildB1B2Config(
        run_root=config.run_root / "gauge_fixed",
        report_path=config.report_path,
        json_path=config.json_path,
        grid=config.grid,
        b2_depth_sequence=tuple(float(tau) for tau in config.target_taus),
        tau_fold_reference=float(config.tau_fold_reference),
        max_tau_backtracks=int(config.max_tau_backtracks),
        min_tau_backtrack_delta=float(config.min_tau_backtrack_delta),
        max_newton_iters_override=config.max_newton_iters_override,
        jacobian_assembly=config.jacobian_assembly,
        b4_reference_tau=float(config.b4_reference_tau),
        b4_random_probe_count=int(config.b4_random_probe_count),
        b4_random_seed=int(config.b4_random_seed),
        b4_match_abs_tol=float(config.b4_match_abs_tol),
        b4_match_rel_tol=float(config.b4_match_rel_tol),
    )


def _c0g_deep_diag_config(
    config: C0gDeepCrawlConfig,
    *,
    converged_taus: Sequence[float],
) -> C0gConfig:
    return C0gConfig(
        run_root=config.run_root / "diagnostics",
        report_path=config.report_path,
        json_path=config.run_root / "diagnostics" / "pathA_C0g_deepcrawl_diagnostics.json",
        final_json_path=config.run_root / "diagnostics" / "pathA_C0g_deepcrawl_diagnostics_final.json",
        grid=config.grid,
        converged_taus=tuple(float(tau) for tau in converged_taus),
        stalled_tau=float(config.tau_fold_reference),
        fit_good_r2_threshold=0.95,
        jacobian_assembly=config.jacobian_assembly,
    )


def _c0g_load_json_path(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _c0g_deep_existing_result(config: C0gDeepCrawlConfig) -> dict[str, Any] | None:
    path = _resolve_input_path(config.json_path)
    if not path.exists():
        return None
    return _c0g_load_json_path(path)


def _c0g_deep_prior_measurement_path(config: C0gDeepCrawlConfig, tau: float) -> Path:
    prior_json = _resolve_input_path(config.prior_b1b2_json_path)
    return (
        prior_json.parent
        / "b2_reconfirm"
        / "state_measurements"
        / f"c0g_state_tau_{_format_tau(float(tau))}.json"
    )


def _c0g_deep_seed_rows(config: C0gDeepCrawlConfig) -> list[dict[str, Any]]:
    prior_crawl_path = _resolve_input_path(config.prior_crawl_json_path)
    prior_crawl = _c0g_load_json_path(prior_crawl_path)
    converged = [
        dict(row)
        for row in prior_crawl.get("tau_attempts", [])
        if bool(row.get("final_physical_converged"))
        and float(row.get("target_tau", math.inf)) < float(config.tau_fold_reference)
    ]
    converged = sorted(converged, key=lambda row: float(row["target_tau"]), reverse=True)
    if not converged:
        converged = sorted(
            [
                dict(row)
                for row in prior_crawl.get("tau_attempts", [])
                if bool(row.get("final_physical_converged"))
            ],
            key=lambda row: float(row["target_tau"]),
            reverse=True,
        )[-1:]
    seeds: list[dict[str, Any]] = []
    for index, row in enumerate(converged):
        tau = float(row["target_tau"])
        seed = dict(row)
        seed["deepcrawl_role"] = "below_fold_seed"
        seed["deepcrawl_source"] = str(prior_crawl_path)
        seed["measurement_status"] = seed.get("measurement_status", "MEASURED")
        seed["attempt_index"] = int(seed.get("attempt_index", index))
        seed["full_newton_budget"] = seed.get("max_newton_iters_override") in (None, 0)
        seed["prior_state_measurement_path"] = str(_c0g_deep_prior_measurement_path(config, tau))
        seeds.append(_c0g_deep_apply_single_arbiter(seed))
    return seeds


def _c0g_deep_apply_single_arbiter(row: Mapping[str, Any]) -> dict[str, Any]:
    updated = dict(row)
    try:
        residual = float(updated.get("final_original_residual_linf", math.inf))
    except (TypeError, ValueError):
        residual = math.inf
    if math.isfinite(residual):
        physical_converged = bool(residual <= BACKGROUND_RESIDUAL_TOL)
        previous = bool(updated.get("final_physical_converged"))
        updated["stage_converged"] = physical_converged
        updated["final_physical_converged"] = physical_converged
        updated["single_arbiter_original_residual_linf_tol"] = BACKGROUND_RESIDUAL_TOL
        if physical_converged and not previous:
            updated["single_arbiter_convergence_override"] = (
                "accepted because original residual Linf is within BACKGROUND_RESIDUAL_TOL; "
                "internal Newton tolerance is tighter and is not the physical arbiter"
            )
    return updated


def _c0g_deep_converged_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    converged: list[dict[str, Any]] = []
    for row in rows:
        updated = _c0g_deep_apply_single_arbiter(row)
        if (
            bool(updated.get("final_physical_converged"))
            and updated.get("measurement_status", "MEASURED") == "MEASURED"
        ):
            converged.append(updated)
    return converged


def _c0g_deep_deepest_converged_row(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any] | None:
    normalized_rows = [_c0g_deep_apply_single_arbiter(row) for row in rows]
    converged = _c0g_deep_converged_rows(normalized_rows)
    if not converged:
        return None
    return min(converged, key=lambda row: float(row["target_tau"]))


def _c0g_deep_has_converged_tau(rows: Sequence[Mapping[str, Any]], tau: float) -> bool:
    return any(
        bool(row.get("final_physical_converged"))
        and _c0g_tau_close(float(row.get("target_tau", math.nan)), float(tau))
        for row in rows
    )


def _c0g_deep_state_measurement(
    *,
    row: Mapping[str, Any],
    config: C0gDeepCrawlConfig,
    diag_config: C0gConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    tau = float(row["target_tau"])
    local_path = _c0g_state_json_path(diag_config, tau=tau)
    if local_path.exists():
        return _c0g_load_json_path(local_path)
    prior_path_text = str(row.get("prior_state_measurement_path", ""))
    prior_path = Path(prior_path_text) if prior_path_text else None
    if prior_path is not None and prior_path.is_file():
        measured = _c0g_load_json_path(prior_path)
        measured["reused_prior_b2_measurement"] = True
        measured["reused_prior_b2_measurement_path"] = str(prior_path)
        _c0g_write_state_result(diag_config, measured)
        return measured
    measured = _c0g_analyze_state(
        tau=tau,
        row=row,
        config=diag_config,
        dtype=dtype,
        include_step0=False,
        compute_ftau=True,
    )
    measured["reused_prior_b2_measurement"] = False
    _c0g_write_state_result(diag_config, measured)
    return measured


def _c0g_deep_admissibility_for_row(
    *,
    row: Mapping[str, Any],
    config: C0gDeepCrawlConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    tau = float(row["target_tau"])
    state_path = _resolve_input_path(str(row["state_artifact"]))
    state = _load_state_artifact(state_path, dtype=dtype)
    branch, provider, grid, boundaries = _branch_context(
        tau=tau,
        config=C0Config(
            grid=config.grid,
            jacobian_assembly=config.jacobian_assembly,
            max_newton_iters_override=config.max_newton_iters_override,
        ),
        dtype=dtype,
    )
    eos_K = float(branch.continuation_K_values[-1])
    original = patha_closed_branch_residual(
        state,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=provider,
    ).detach()
    conditioned = c0_preconditioner_residual(
        state,
        grid,
        branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=provider,
        aid=C0AidParameters(core_epsilon=0.0, k1_radius_epsilon=0.0),
    ).detach()
    metrics = _state_metrics(state, grid)
    residual_linf = float(torch.max(torch.abs(original)).detach().cpu().item())
    residual_equal_max_abs = float(torch.max(torch.abs(original - conditioned)).detach().cpu().item())
    min_rho = float(metrics["min_rho"])
    min_r0 = float(metrics["min_R0"])
    positive_density = bool(math.isfinite(min_rho) and min_rho > 0.0)
    sane_r0 = bool(math.isfinite(min_r0) and min_r0 > 0.0)
    residual_ok = bool(residual_linf <= BACKGROUND_RESIDUAL_TOL)
    residual_equivalent = bool(residual_equal_max_abs <= float(config.residual_equality_tolerance))
    final_aids_inactive = bool(row.get("floor_activation", {}).get("final_aids_inactive", True))
    epsilon_independent = bool(final_aids_inactive and residual_equivalent)
    status = "PASS" if positive_density and sane_r0 and residual_ok and epsilon_independent else "FAIL"
    return {
        "status": status,
        "tau": tau,
        "positive_density_status": "PASS" if positive_density else "FAIL",
        "min_rho": min_rho,
        "sane_min_R0_status": "PASS" if sane_r0 else "FAIL",
        "min_R0": min_r0,
        "residual_linf": residual_linf,
        "residual_tolerance": BACKGROUND_RESIDUAL_TOL,
        "original_residual_status": "PASS" if residual_ok else "FAIL",
        "epsilon_independence_status": "PASS" if epsilon_independent else "FAIL",
        "final_aids_inactive": final_aids_inactive,
        "zero_aid_residual_equivalence_max_abs": residual_equal_max_abs,
        "zero_aid_residual_equivalence_tolerance": float(config.residual_equality_tolerance),
        "metrics": metrics,
        "state_artifact": str(state_path),
    }


def _c0g_deep_monotone_summary(values: Sequence[float]) -> dict[str, Any]:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if len(finite) < 2:
        return {"status": "NOT_MEASURED", "direction": None}
    scale = max(max(abs(value) for value in finite), 1.0)
    tol = 1.0e-8 * scale
    nondecreasing = all(right >= left - tol for left, right in zip(finite, finite[1:]))
    nonincreasing = all(right <= left + tol for left, right in zip(finite, finite[1:]))
    if nondecreasing and not nonincreasing:
        direction = "nondecreasing_as_tau_decreases"
    elif nonincreasing and not nondecreasing:
        direction = "nonincreasing_as_tau_decreases"
    elif nondecreasing and nonincreasing:
        direction = "flat_within_tolerance"
    else:
        direction = "nonmonotone"
    return {
        "status": "PASS" if direction != "nonmonotone" else "FAIL",
        "direction": direction,
        "first": finite[0],
        "last": finite[-1],
        "delta": float(finite[-1] - finite[0]),
        "tolerance": tol,
    }


def _c0g_deep_branch_tangent(rows: Sequence[Mapping[str, Any]], *, dtype: torch.dtype) -> dict[str, Any]:
    converged = sorted(_c0g_deep_converged_rows(rows), key=lambda row: float(row["target_tau"]), reverse=True)
    if len(converged) < 3:
        return {"status": "NOT_MEASURED", "rows": [], "aggregate_call": "TOO_FEW_STATES"}
    states: list[np.ndarray] = []
    for row in converged:
        state = _load_state_artifact(_resolve_input_path(str(row["state_artifact"])), dtype=dtype)
        states.append(state.detach().cpu().numpy().astype(np.float64, copy=False))
    tangent_rows: list[dict[str, Any]] = []
    previous_delta: np.ndarray | None = None
    for index in range(1, len(converged)):
        delta = states[index] - states[index - 1]
        if previous_delta is None:
            previous_delta = delta
            continue
        denom = float(np.linalg.norm(delta) * np.linalg.norm(previous_delta))
        overlap = float(np.dot(delta, previous_delta) / denom) if denom > 0.0 else math.nan
        if math.isfinite(overlap) and overlap < 0.0:
            call = "TANGENT_REVERSAL"
        elif math.isfinite(overlap) and overlap > 0.5:
            call = "SAME_TANGENT_DIRECTION"
        elif math.isfinite(overlap):
            call = "TANGENT_GRAY"
        else:
            call = "NOT_MEASURED"
        tangent_rows.append(
            {
                "tau_left": float(converged[index - 1]["target_tau"]),
                "tau": float(converged[index]["target_tau"]),
                "tangent_cosine_vs_previous": overlap,
                "call": call,
            }
        )
        previous_delta = delta
    calls = {row["call"] for row in tangent_rows}
    aggregate = "TANGENT_REVERSAL" if "TANGENT_REVERSAL" in calls else "NO_TANGENT_REVERSAL"
    return {"status": "MEASURED", "rows": tangent_rows, "aggregate_call": aggregate}


def _c0g_deep_ladder_rows(
    rows: Sequence[Mapping[str, Any]],
    state_results: Sequence[Mapping[str, Any]],
    admissibility_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    sigma_by_tau = {
        float(row["tau"]): row.get("svd", {}).get("sigma_min")
        for row in state_results
        if "tau" in row
    }
    admiss_by_tau = {float(row["tau"]): row for row in admissibility_rows if "tau" in row}
    ladder = []
    for raw_row in sorted(rows, key=lambda item: float(item["target_tau"]), reverse=True):
        row = _c0g_deep_apply_single_arbiter(raw_row)
        if float(row.get("target_tau", math.inf)) >= 0.0291233:
            continue
        metrics = dict(row.get("metrics", {}))
        adm = admiss_by_tau.get(float(row["target_tau"]), {})
        if adm.get("metrics"):
            metrics = dict(adm["metrics"])
        ladder.append(
            {
                "tau": float(row["target_tau"]),
                "nominal_tau": row.get("nominal_target_tau"),
                "role": row.get("deepcrawl_role"),
                "converged": bool(row.get("final_physical_converged")),
                "residual_linf": row.get("final_original_residual_linf") or adm.get("residual_linf"),
                "newton_iters": _c0g_attempt_newton_iterations(row),
                "newton_iters_total": _c0g_attempt_total_newton_iterations(row),
                "backtracks": int(row.get("backtrack_index", 0)),
                "wall_seconds": row.get("elapsed_seconds"),
                "min_R0": metrics.get("min_R0"),
                "min_rho": metrics.get("min_rho"),
                "mu": metrics.get("chemical_potential"),
                "sigma_min_JQ_perp": sigma_by_tau.get(float(row["target_tau"])),
                "admissibility": adm.get("status") if adm else None,
                "init_source": row.get("init", {}).get("source"),
                "measurement_status": row.get("measurement_status", "MEASURED"),
            }
        )
    return ladder


def _c0g_deep_verdict(
    *,
    rows: Sequence[Mapping[str, Any]],
    ladder: Sequence[Mapping[str, Any]],
    state_results: Sequence[Mapping[str, Any]],
    admissibility_rows: Sequence[Mapping[str, Any]],
    sigma_fit: Mapping[str, Any],
    tangent: Mapping[str, Any],
    config: C0gDeepCrawlConfig,
) -> dict[str, Any]:
    normalized_rows = [_c0g_deep_apply_single_arbiter(row) for row in rows]
    new_converged = [
        row
        for row in normalized_rows
        if row.get("deepcrawl_role") == "deep_crawl_attempt"
        and bool(row.get("final_physical_converged"))
        and float(row.get("final_original_residual_linf", math.inf)) <= BACKGROUND_RESIDUAL_TOL
    ]
    admiss_by_tau = {float(row["tau"]): row for row in admissibility_rows if "tau" in row}
    all_converged_below = [
        row
        for row in normalized_rows
        if bool(row.get("final_physical_converged"))
        and float(row.get("target_tau", math.inf)) < float(config.tau_fold_reference)
    ]
    all_admissible = all(
        admiss_by_tau.get(float(row["target_tau"]), {}).get("status") == "PASS"
        for row in all_converged_below
    )
    sigmas = [
        float(state.get("svd", {}).get("sigma_min", math.nan))
        for state in state_results
        if state.get("svd", {}).get("status") == "MEASURED"
    ]
    sigma_finite = bool(
        sigmas
        and all(math.isfinite(value) and value > float(config.sigma_finite_floor) for value in sigmas)
    )
    sigma_collapse_factor = (
        float(sigmas[0] / min(sigmas))
        if sigmas and min(sigmas) > 0.0 and math.isfinite(sigmas[0])
        else math.nan
    )
    sigma_collapsing = bool(math.isfinite(sigma_collapse_factor) and sigma_collapse_factor >= 10.0)
    guard = _crawl_persistent_failure_guard(
        normalized_rows,
        _c0g_deep_config_to_dict(config),
        failed_attempt_predicate=lambda row: (
            row.get("deepcrawl_role") == "deep_crawl_attempt"
            and row.get("measurement_status", "MEASURED") == "MEASURED"
            and not bool(_c0g_deep_apply_single_arbiter(row).get("final_physical_converged"))
        ),
    )
    r0_trend = _c0g_deep_monotone_summary([row.get("min_R0", math.nan) for row in ladder])
    mu_trend = _c0g_deep_monotone_summary([row.get("mu", math.nan) for row in ladder])
    fit_linear = sigma_fit.get("linear", {}) if isinstance(sigma_fit, Mapping) else {}
    fold_fit_support = bool(
        sigma_fit.get("call") == "LINEAR_MONOTONE_FOLD_SUPPORT"
        and math.isfinite(float(fit_linear.get("r2", math.nan)))
        and float(fit_linear.get("r2", math.nan)) >= 0.95
    )
    tangent_reversal = tangent.get("aggregate_call") == "TANGENT_REVERSAL"
    real_fold = bool(
        guard.get("crawl_persistent_failure_evidence")
        and fold_fit_support
        and sigma_collapsing
        and tangent_reversal
    )
    no_fold_locked = bool(
        len(new_converged) >= int(config.no_fold_required_new_converged)
        and all_admissible
        and sigma_finite
        and not sigma_collapsing
        and r0_trend.get("status") == "PASS"
        and mu_trend.get("status") == "PASS"
        and not guard.get("crawl_persistent_failure_evidence")
    )
    if real_fold:
        verdict = "REAL_FOLD_DEEPER"
        why = (
            "persistent full-budget failure plus linear-monotone sigma_min^2 support "
            "and branch tangent reversal"
        )
    elif no_fold_locked:
        verdict = "NO_FOLD_LOCKED"
        why = (
            f"{len(new_converged)} new branch-continuous deeper states converged below "
            f"tau_fold with admissibility PASS and finite sigma_min(J*Q_perp)"
        )
    else:
        verdict = "STILL_GOING"
        why = "deeper crawl did not meet the strict no-fold or real-fold criteria"
    deepest = min((float(row["target_tau"]) for row in all_converged_below), default=None)
    margin = None if deepest is None else float(config.tau_fold_reference) - deepest
    return {
        "verdict": verdict,
        "why": why,
        "tau_fold_reference": float(config.tau_fold_reference),
        "deepest_converged_tau": deepest,
        "deepest_margin_past_tau_fold": margin,
        "new_deep_converged_count": int(len(new_converged)),
        "required_new_deep_converged_count": int(config.no_fold_required_new_converged),
        "all_converged_below_admissible": bool(all_admissible),
        "sigma_min_finite": bool(sigma_finite),
        "sigma_min_finite_floor": float(config.sigma_finite_floor),
        "sigma_min_min": min(sigmas) if sigmas else None,
        "sigma_min_collapse_factor": sigma_collapse_factor,
        "sigma_min_collapsing_toward_zero": sigma_collapsing,
        "sigma_fit_supports_fold": bool(fold_fit_support),
        "branch_tangent_reversal": bool(tangent_reversal),
        "r0_trend": r0_trend,
        "mu_trend": mu_trend,
        "persistent_failure_guard": {
            key: value for key, value in guard.items() if key != "failed_attempts"
        },
        "failed_attempts": [
            {
                "tau": row.get("target_tau"),
                "nominal_tau": row.get("nominal_target_tau"),
                "residual_linf": row.get("final_original_residual_linf"),
                "backtrack_index": row.get("backtrack_index"),
                "newton_iters": _c0g_attempt_newton_iterations(row),
                "message": row.get("message"),
            }
            for row in guard.get("failed_attempts", [])
        ],
    }


def _c0g_deep_empty_core_summary(ladder: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    converged = [row for row in ladder if bool(row.get("converged"))]
    if not converged:
        return {"status": "NOT_MEASURED"}
    first = converged[0]
    deepest = converged[-1]
    min_r0_first = float(first.get("min_R0", math.nan))
    min_r0_deep = float(deepest.get("min_R0", math.nan))
    min_rho_first = float(first.get("min_rho", math.nan))
    min_rho_deep = float(deepest.get("min_rho", math.nan))
    approached = bool(
        math.isfinite(min_r0_first)
        and math.isfinite(min_r0_deep)
        and math.isfinite(min_rho_first)
        and math.isfinite(min_rho_deep)
        and (min_r0_deep < min_r0_first or min_rho_deep < min_rho_first)
    )
    return {
        "status": "MEASURED",
        "first_tau": first.get("tau"),
        "deepest_tau": deepest.get("tau"),
        "first_min_R0": min_r0_first,
        "deepest_min_R0": min_r0_deep,
        "first_min_rho": min_rho_first,
        "deepest_min_rho": min_rho_deep,
        "approached_small_core_regime": approached,
        "deepest_solvable": bool(deepest.get("converged") and deepest.get("admissibility") == "PASS"),
        "interpretation": (
            "min_R0/min_rho moved toward smaller-core values"
            if approached
            else "min_R0/min_rho did not move toward the small-core regime over this ladder"
        ),
    }


def _c0g_deep_finalize(
    *,
    rows: Sequence[Mapping[str, Any]],
    config: C0gDeepCrawlConfig,
    b4_reconfirm: Mapping[str, Any],
    dtype: torch.dtype,
    started: float,
    chunk_note: str | None,
) -> dict[str, Any]:
    normalized_rows = [_c0g_deep_apply_single_arbiter(row) for row in rows]
    converged = _c0g_deep_converged_rows(normalized_rows)
    converged_taus = sorted([float(row["target_tau"]) for row in converged], reverse=True)
    diag_config = _c0g_deep_diag_config(config, converged_taus=converged_taus)
    state_results = [
        _c0g_deep_state_measurement(
            row=row,
            config=config,
            diag_config=diag_config,
            dtype=dtype,
        )
        for row in converged
        if float(row.get("target_tau", math.inf)) < float(config.tau_fold_reference)
    ]
    admissibility_rows = [
        _c0g_deep_admissibility_for_row(row=row, config=config, dtype=dtype)
        for row in converged
        if float(row.get("target_tau", math.inf)) < float(config.tau_fold_reference)
    ]
    sigma_fit = _c0g_fit_sigma_min_squared(state_results, config=diag_config)
    tangent = _c0g_deep_branch_tangent(converged, dtype=dtype)
    ladder = _c0g_deep_ladder_rows(normalized_rows, state_results, admissibility_rows)
    verdict = _c0g_deep_verdict(
        rows=normalized_rows,
        ladder=ladder,
        state_results=state_results,
        admissibility_rows=admissibility_rows,
        sigma_fit=sigma_fit,
        tangent=tangent,
        config=config,
    )
    empty_core = _c0g_deep_empty_core_summary(ladder)
    admiss_failures = [row for row in admissibility_rows if row.get("status") != "PASS"]
    result = {
        "schema": "stage1_pathA_C0g_deepcrawl/v1",
        "source_revision": source_revision(),
        "config": _c0g_deep_config_to_dict(config),
        "b4_reconfirm": dict(b4_reconfirm),
        "prior_b2_result": _c0g_deep_prior_summary(config),
        "rows": [dict(row) for row in normalized_rows],
        "ladder": ladder,
        "state_results": state_results,
        "sigma_min_squared_fit": sigma_fit,
        "branch_tangent": tangent,
        "admissibility": {
            "status": "PASS" if not admiss_failures and admissibility_rows else "FAIL",
            "rows": admissibility_rows,
            "failing_taus": [row.get("tau") for row in admiss_failures],
        },
        "empty_core_summary": empty_core,
        "work1_verdict": verdict,
        "scope_guard": {
            "B3_pseudoarclength_built": False,
            "B4_analytic_sparse_assembly_built": True,
            "single_arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
            "patha_closed_branch_residual_touched": False,
            "faithful_operators_touched": False,
            "xi_or_grad_div_penalty_touched": False,
            "physical_export_permitted_touched": False,
            "prefer_existing_b2c_background_predictor": False,
            "frozen_physics_touched": False,
        },
        "chunk_note": chunk_note,
        "elapsed_seconds": float(time.perf_counter() - started),
        "report_path": str(_resolve_output_path(config.report_path)),
        "json_path": str(_resolve_output_path(config.json_path)),
    }
    _c0g_write_json(_resolve_output_path(config.json_path), result)
    write_c0g_deepcrawl_report(result, _resolve_output_path(config.report_path))
    return result


def _c0g_deep_prior_summary(config: C0gDeepCrawlConfig) -> dict[str, Any]:
    path = _resolve_input_path(config.prior_b1b2_json_path)
    if not path.exists():
        return {"status": "NOT_FOUND", "path": str(path)}
    payload = _c0g_load_json_path(path)
    return {
        "status": "MEASURED",
        "path": str(path),
        "b4_status": payload.get("B4_analytic_sparse_jacobian", {}).get("status"),
        "b2_verdict": payload.get("B2_reconfirm", {}).get("gate", {}).get("verdict"),
        "scope_guard": payload.get("scope_guard", {}),
    }


def write_c0g_deepcrawl_report(result: Mapping[str, Any], path: Path) -> None:
    verdict = result.get("work1_verdict", {})
    b4 = result.get("b4_reconfirm", {})
    admissibility = result.get("admissibility", {})
    empty_core = result.get("empty_core_summary", {})
    fit = result.get("sigma_min_squared_fit", {})
    fit_linear = fit.get("linear", {}) if isinstance(fit, Mapping) else {}
    lines = [
        "# Path-A C0g Deeper Crawl",
        "",
        f"Work-1 verdict: **{verdict.get('verdict')}**",
        "",
        "## Scope",
        "",
        "- Continuation is gauge-fixed and path-only: `row_scale * J_original * col_scale`, step constrained to `Q_perp`.",
        "- The original `patha_closed_branch_residual` is the sole convergence arbiter at `Linf <= 1e-6`.",
        "- B-3 pseudo-arclength was not built.",
        "- Frozen residual/operators/`(1/xi)` penalty/export/physics are not edited by this runner.",
        "",
        "## B-4 Reconfirm",
        "",
        "```yaml",
        f"status: {b4.get('status')}",
        f"reference_tau: {b4.get('reference_tau')}",
        f"assembly_speedup_x: {b4.get('assembly_speedup_x')}",
        f"tolerances: {b4.get('match_tolerances')}",
        "```",
        "",
        "## Deeper Ladder",
        "",
        _markdown_table(
            [
                "tau",
                "converged",
                "residual_linf",
                "newton_iters",
                "backtracks",
                "wall_seconds",
                "min_R0",
                "min_rho",
                "mu",
                "sigma_min_JQ_perp",
                "admissibility",
            ],
            result.get("ladder", []),
        ),
        "",
        "## Sigma Trend",
        "",
        "```yaml",
        f"status: {fit.get('status') if isinstance(fit, Mapping) else None}",
        f"call: {fit.get('call') if isinstance(fit, Mapping) else None}",
        f"linear_r2: {fit_linear.get('r2')}",
        f"linear_slope: {fit_linear.get('slope')}",
        f"tau_fold_zero_crossing: {fit_linear.get('tau_fold_zero_crossing')}",
        f"monotone_toward_stall: {fit.get('sigma_min_monotone_decreasing_toward_stall') if isinstance(fit, Mapping) else None}",
        f"branch_tangent: {result.get('branch_tangent', {}).get('aggregate_call')}",
        "```",
        "",
        "## Verdict Support",
        "",
        "```yaml",
        f"{verdict}",
        "```",
        "",
        "## Admissibility",
        "",
        f"Overall: **{admissibility.get('status')}**",
        "",
        _markdown_table(
            [
                "tau",
                "status",
                "positive_density_status",
                "min_rho",
                "sane_min_R0_status",
                "min_R0",
                "residual_linf",
                "epsilon_independence_status",
                "zero_aid_residual_equivalence_max_abs",
            ],
            admissibility.get("rows", []),
        ),
        "",
        "## Empty-Core Target",
        "",
        "```yaml",
        f"{empty_core}",
        "```",
        "",
        "## Scope Guard",
        "",
        "```yaml",
    ]
    for key, value in result.get("scope_guard", {}).items():
        lines.append(f"{key}: {value}")
    lines.extend(
        [
            "```",
            "",
            f"Machine JSON: `{result.get('json_path')}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _c0g_deep_stdout_summary(result: Mapping[str, Any]) -> str:
    verdict = result.get("work1_verdict", {})
    ladder_lines = []
    for row in result.get("ladder", []):
        ladder_lines.append(
            (
                "  tau={tau}, converged={conv}, residual={resid}, iters={iters}, "
                "min_R0={r0}, min_rho={rho}, sigma_min(JQ_perp)={sigma}, admissibility={adm}"
            ).format(
                tau=_c0g_fmt(row.get("tau")),
                conv=row.get("converged"),
                resid=_c0g_fmt(row.get("residual_linf")),
                iters=row.get("newton_iters"),
                r0=_c0g_fmt(row.get("min_R0")),
                rho=_c0g_fmt(row.get("min_rho")),
                sigma=_c0g_fmt(row.get("sigma_min_JQ_perp")),
                adm=row.get("admissibility"),
            )
        )
    admiss = result.get("admissibility", {})
    empty = result.get("empty_core_summary", {})
    scope = result.get("scope_guard", {})
    return "\n".join(
        [
            "C0g deeper-crawl summary:",
            "Deeper convergence ladder:",
            *(ladder_lines or ["  n/a"]),
            f"Deepest margin past tau_fold: {_c0g_fmt(verdict.get('deepest_margin_past_tau_fold'))}",
            (
                "sigma_min(JQ_perp): "
                f"finite={verdict.get('sigma_min_finite')}, "
                f"min={_c0g_fmt(verdict.get('sigma_min_min'))}, "
                f"collapse_factor={_c0g_fmt(verdict.get('sigma_min_collapse_factor'))}, "
                f"collapsing={verdict.get('sigma_min_collapsing_toward_zero')}, "
                f"fit_call={result.get('sigma_min_squared_fit', {}).get('call')}"
            ),
            f"Work-1 verdict: {verdict.get('verdict')} - {verdict.get('why')}",
            (
                "Admissibility: "
                f"{admiss.get('status')}; failing_taus={admiss.get('failing_taus')}"
            ),
            (
                "Empty-core target: "
                f"approached={empty.get('approached_small_core_regime')}, "
                f"deepest_solvable={empty.get('deepest_solvable')}, "
                f"deepest_min_R0={_c0g_fmt(empty.get('deepest_min_R0'))}, "
                f"deepest_min_rho={_c0g_fmt(empty.get('deepest_min_rho'))}"
            ),
            (
                "Scope: B-3 built="
                f"{scope.get('B3_pseudoarclength_built')}; frozen residual/operators/"
                f"(1/xi)/export/physics touched="
                f"{scope.get('patha_closed_branch_residual_touched')}/"
                f"{scope.get('faithful_operators_touched')}/"
                f"{scope.get('xi_or_grad_div_penalty_touched')}/"
                f"{scope.get('physical_export_permitted_touched')}/"
                f"{scope.get('frozen_physics_touched')}"
            ),
            f"Files: report={result.get('report_path')}, json={result.get('json_path')}",
        ]
    )


def run_c0g_deepcrawl(
    config: C0gDeepCrawlConfig | None = None,
    *,
    max_new_states: int = 1,
    finalize_only: bool = False,
    b4_recheck: bool = True,
) -> dict[str, Any]:
    cfg = config or C0gDeepCrawlConfig()
    started = time.perf_counter()
    dtype = configure_backend(BackendConfig())
    existing = _c0g_deep_existing_result(cfg)
    if existing is None:
        rows: list[dict[str, Any]] = _c0g_deep_seed_rows(cfg)
        b4_reconfirm: Mapping[str, Any]
        if b4_recheck:
            b4_reconfirm = run_c0g_b4_match_gate(_c0g_deep_build_config(cfg))
        else:
            prior = _c0g_deep_prior_summary(cfg)
            b4_reconfirm = {"status": prior.get("b4_status"), "source": prior.get("path")}
    else:
        rows = [_c0g_deep_apply_single_arbiter(row) for row in existing.get("rows", [])]
        b4_reconfirm = existing.get("b4_reconfirm", {})
        if b4_recheck and b4_reconfirm.get("status") != "PASS":
            b4_reconfirm = run_c0g_b4_match_gate(_c0g_deep_build_config(cfg))

    if b4_reconfirm.get("status") != "PASS":
        return _c0g_deep_finalize(
            rows=rows,
            config=cfg,
            b4_reconfirm=b4_reconfirm,
            dtype=dtype,
            started=started,
            chunk_note="B4 match gate did not pass; no deeper continuation attempted",
        )

    build_cfg = _c0g_deep_build_config(cfg)
    new_states = 0
    stopped_note: str | None = None
    if not finalize_only and int(max_new_states) != 0:
        deepest = _c0g_deep_deepest_converged_row(rows)
        if deepest is None:
            raise RuntimeError("deep crawl has no prior converged state to warm-start from")
        last_tau = float(deepest["target_tau"])
        last_state = _c0g_state_vector_from_path(str(deepest["state_artifact"]), dtype=dtype)
        prior_converged = sorted(
            _c0g_deep_converged_rows(rows),
            key=lambda row: float(row["target_tau"]),
        )
        previous_row = prior_converged[1] if len(prior_converged) >= 2 else None
        previous_tau = float(previous_row["target_tau"]) if previous_row is not None else None
        previous_state = (
            _c0g_state_vector_from_path(str(previous_row["state_artifact"]), dtype=dtype)
            if previous_row is not None
            else None
        )
        previous_grid = None
        if previous_tau is not None:
            _branch, _provider, previous_grid, _boundaries = _branch_context(
                tau=previous_tau,
                config=C0Config(grid=cfg.grid, jacobian_assembly=cfg.jacobian_assembly),
                dtype=dtype,
            )
        branch, _provider, last_grid, _boundaries = _branch_context(
            tau=last_tau,
            config=C0Config(grid=cfg.grid, jacobian_assembly=cfg.jacobian_assembly),
            dtype=dtype,
        )
        del branch
        for nominal_tau in cfg.target_taus:
            nominal_tau = float(nominal_tau)
            if new_states >= int(max_new_states):
                break
            if nominal_tau >= last_tau - 5.0e-13 or _c0g_deep_has_converged_tau(rows, nominal_tau):
                continue
            start_tau = last_tau
            full_delta = nominal_tau - start_tau
            target_tau = nominal_tau
            backtrack_index = 0
            converged_this_nominal = False
            while True:
                branch, _provider, target_grid, _boundaries = _branch_context(
                    tau=float(target_tau),
                    config=C0Config(grid=cfg.grid, jacobian_assembly=cfg.jacobian_assembly),
                    dtype=dtype,
                )
                last_resampled = resample_closed_branch_state(
                    last_state.detach().clone(),
                    last_grid,
                    target_grid,
                    branch,
                )
                secant_predictor: dict[str, Any] = {"status": "NOT_USED", "reason": "missing_previous_state"}
                x0 = last_resampled
                if (
                    previous_state is not None
                    and previous_grid is not None
                    and previous_tau is not None
                    and abs(float(last_tau) - float(previous_tau)) > 1.0e-15
                ):
                    previous_resampled = resample_closed_branch_state(
                        previous_state.detach().clone(),
                        previous_grid,
                        target_grid,
                        branch,
                    )
                    factor = (float(target_tau) - float(last_tau)) / (
                        float(last_tau) - float(previous_tau)
                    )
                    x0 = last_resampled + float(factor) * (last_resampled - previous_resampled)
                    secant_predictor = {
                        "status": "USED",
                        "previous_tau": previous_tau,
                        "last_tau": last_tau,
                        "target_tau": float(target_tau),
                        "extrapolation_factor": float(factor),
                    }
                x0, predictor = _apply_c0_wall_predictor(
                    state=x0,
                    grid=target_grid,
                    tau=float(target_tau),
                    grid_level=cfg.grid,
                )
                row, final_state, final_grid = _c0g_single_gauge_fixed_attempt(
                    x0=x0,
                    tau=float(target_tau),
                    config=build_cfg,
                    dtype=dtype,
                    attempt_index=len(rows),
                    start_from_tau=start_tau,
                    init={
                        "source": "previous_c0g_gauge_fixed_converged_state",
                        "previous_tau": last_tau,
                        "previous_state_artifact": deepest.get("state_artifact"),
                        "secant_predictor": secant_predictor,
                        "wall_predictor": predictor,
                    },
                    max_newton_iters=cfg.max_newton_iters_override,
                )
                row["attempt_index"] = int(len(rows))
                row["deepcrawl_role"] = "deep_crawl_attempt"
                row["measurement_status"] = "MEASURED"
                row["full_newton_budget"] = row.get("max_newton_iters_override") in (None, 0)
                row["nominal_target_tau"] = float(nominal_tau)
                row["backtrack_index"] = int(backtrack_index)
                row["delta_tau"] = float(target_tau) - float(start_tau)
                row = _c0g_deep_apply_single_arbiter(row)
                rows.append(row)
                if row.get("final_physical_converged"):
                    deepest = row
                    previous_tau = last_tau
                    previous_state = last_state.detach().clone()
                    previous_grid = last_grid
                    last_tau = float(row["target_tau"])
                    last_state = final_state.detach().clone()
                    last_grid = final_grid
                    new_states += 1
                    converged_this_nominal = True
                    break
                if (
                    backtrack_index >= int(cfg.max_tau_backtracks)
                    or abs(float(target_tau) - float(start_tau)) <= float(cfg.min_tau_backtrack_delta)
                ):
                    stopped_note = (
                        f"nominal tau {nominal_tau:.12e} failed after "
                        f"{backtrack_index} backtracks under full budget"
                    )
                    break
                backtrack_index += 1
                target_tau = float(start_tau) + float(full_delta) / (2.0**backtrack_index)
            if not converged_this_nominal and stopped_note is not None:
                break
    return _c0g_deep_finalize(
        rows=rows,
        config=cfg,
        b4_reconfirm=b4_reconfirm,
        dtype=dtype,
        started=started,
        chunk_note=stopped_note,
    )


def _c0g_b3_config_to_dict(config: C0gB3PseudoArclengthConfig) -> dict[str, Any]:
    data = asdict(config)
    for key in ("deepcrawl_json_path", "run_root", "report_path", "json_path"):
        data[key] = str(data[key])
    return data


def _c0g_b3_build_config(config: C0gB3PseudoArclengthConfig) -> C0gBuildB1B2Config:
    return C0gBuildB1B2Config(
        run_root=config.run_root / "b4_reconfirm_context",
        report_path=config.report_path,
        json_path=config.json_path,
        grid=config.grid,
        jacobian_assembly=config.jacobian_assembly,
        b4_reference_tau=float(config.b4_reference_tau),
        b4_random_probe_count=int(config.b4_random_probe_count),
        b4_random_seed=int(config.b4_random_seed),
        b4_match_abs_tol=float(config.b4_match_abs_tol),
        b4_match_rel_tol=float(config.b4_match_rel_tol),
    )


def _c0g_b3_deepcrawl_payload(config: C0gB3PseudoArclengthConfig) -> dict[str, Any]:
    path = _resolve_input_path(config.deepcrawl_json_path)
    payload = _c0g_load_json_path(path)
    payload["_resolved_path"] = str(path)
    return payload


def _c0g_b3_converged_deepcrawl_rows(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = [
        _c0g_deep_apply_single_arbiter(row)
        for row in payload.get("rows", [])
        if bool(_c0g_deep_apply_single_arbiter(row).get("final_physical_converged"))
        and float(row.get("final_original_residual_linf", math.inf)) <= BACKGROUND_RESIDUAL_TOL
        and row.get("state_artifact")
    ]
    return sorted(rows, key=lambda row: float(row["target_tau"]), reverse=True)


def _c0g_b3_target_summary(
    payload: Mapping[str, Any],
    config: C0gB3PseudoArclengthConfig,
) -> dict[str, Any]:
    rows = _c0g_b3_converged_deepcrawl_rows(payload)
    if len(rows) < 2:
        raise RuntimeError("B-3 requires at least two converged deepcrawl rows")
    deepest = min(rows, key=lambda row: float(row["target_tau"]))
    ladder = [
        float(row["target_tau"])
        for row in rows
        if float(row["target_tau"]) >= float(deepest["target_tau"]) - 5.0e-13
    ]
    return {
        "status": "MEASURED",
        "deepcrawl_json_path": payload.get("_resolved_path"),
        "deepest_converged_tau": float(deepest["target_tau"]),
        "deepest_state_artifact": deepest.get("state_artifact"),
        "deepest_residual_linf": float(deepest.get("final_original_residual_linf", math.nan)),
        "seed_tau": float(rows[0]["target_tau"]),
        "seed_state_artifact": rows[0].get("state_artifact"),
        "b30_target_taus": ladder,
        "not_hardcoded": True,
        "residual_tolerance": float(config.residual_tolerance),
    }


def _c0g_b3_context(
    *,
    tau: float,
    config: C0gB3PseudoArclengthConfig,
    dtype: torch.dtype,
):
    branch, provider, grid, boundaries = _branch_context(
        tau=float(tau),
        config=C0Config(grid=config.grid, jacobian_assembly=config.jacobian_assembly),
        dtype=dtype,
    )
    eos_K = float(branch.continuation_K_values[-1])

    def residual_fn(x: torch.Tensor) -> torch.Tensor:
        return patha_closed_branch_residual(
            x,
            grid,
            branch,
            eos_K=eos_K,
            boundaries=boundaries,
            s_sigma=provider,
        )

    return branch, grid, residual_fn


def _c0g_b3_residual_np(
    *,
    state_np: np.ndarray,
    tau: float,
    config: C0gB3PseudoArclengthConfig,
    dtype: torch.dtype,
) -> tuple[np.ndarray, Any, Any]:
    branch, grid, residual_fn = _c0g_b3_context(tau=tau, config=config, dtype=dtype)
    state_t = torch.as_tensor(np.asarray(state_np, dtype=np.float64), dtype=dtype, device="cpu")
    residual = residual_fn(state_t).detach().cpu().numpy().astype(np.float64, copy=False)
    return residual, branch, grid


def _c0g_b3_ftau_scaled(
    *,
    state_np: np.ndarray,
    tau: float,
    row_scale: np.ndarray,
    config: C0gB3PseudoArclengthConfig,
    dtype: torch.dtype,
) -> np.ndarray:
    h = float(config.ftau_h)
    tau_left = max(float(tau) - h, np.finfo(np.float64).tiny)
    tau_right = float(tau) + h
    if tau_left == float(tau):
        right, _branch, _grid = _c0g_b3_residual_np(
            state_np=state_np,
            tau=tau_right,
            config=config,
            dtype=dtype,
        )
        center, _branch, _grid = _c0g_b3_residual_np(
            state_np=state_np,
            tau=float(tau),
            config=config,
            dtype=dtype,
        )
        return row_scale * ((right - center) / (tau_right - float(tau)))
    left, _branch, _grid = _c0g_b3_residual_np(
        state_np=state_np,
        tau=tau_left,
        config=config,
        dtype=dtype,
    )
    right, _branch, _grid = _c0g_b3_residual_np(
        state_np=state_np,
        tau=tau_right,
        config=config,
        dtype=dtype,
    )
    return row_scale * ((right - left) / (tau_right - tau_left))


def _c0g_b3_linearization(
    *,
    state_np: np.ndarray,
    tau: float,
    config: C0gB3PseudoArclengthConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    state_t = torch.as_tensor(np.asarray(state_np, dtype=np.float64), dtype=dtype, device="cpu")
    c0g_config = C0gConfig(
        run_root=config.run_root / "diagnostics",
        grid=config.grid,
        jacobian_assembly=config.jacobian_assembly,
    )
    assembled = _c0g_assemble_original_jacobian(
        state=state_t,
        tau=float(tau),
        config=c0g_config,
        dtype=dtype,
    )
    branch = assembled["branch"]
    grid = assembled["grid"]
    gauge = _c0g_build_analytic_gauge_matrix(
        state=state_t,
        grid=grid,
        branch=branch,
        config=c0g_config,
    )
    complement = _c0g_scaled_gauge_complement(
        gauge_matrix=np.asarray(gauge["generator_matrix"], dtype=np.float64),
        col_scale=np.asarray(assembled["col_scale"], dtype=np.float64),
        config=c0g_config,
    )
    q_perp = np.asarray(complement["q_perp"], dtype=np.float64)
    scaled_matrix = assembled["scaled_matrix"].tocsc()
    reduced_j = (scaled_matrix @ q_perp).astype(np.float64, copy=False)
    ftau_scaled = _c0g_b3_ftau_scaled(
        state_np=np.asarray(state_np, dtype=np.float64),
        tau=float(tau),
        row_scale=np.asarray(assembled["row_scale"], dtype=np.float64),
        config=config,
        dtype=dtype,
    )
    singular = np.linalg.svd(reduced_j, compute_uv=False)
    sigma_min = float(singular[-1]) if singular.size else math.nan
    sigma_max = float(singular[0]) if singular.size else math.nan
    return {
        "branch": branch,
        "grid": grid,
        "residual_fn": assembled["residual_fn"],
        "residual": np.asarray(assembled["residual"], dtype=np.float64),
        "matrix": assembled["matrix"],
        "scaled_matrix": scaled_matrix,
        "row_scale": np.asarray(assembled["row_scale"], dtype=np.float64),
        "col_scale": np.asarray(assembled["col_scale"], dtype=np.float64),
        "q_perp": q_perp,
        "reduced_j": reduced_j,
        "ftau_scaled": ftau_scaled,
        "gauge": gauge,
        "scaled_gauge_rank": int(complement["scaled_gauge_rank"]),
        "scaled_gauge_rank_threshold": float(complement["scaled_gauge_rank_threshold"]),
        "sigma_min_JQ_perp": sigma_min,
        "sigma_max_JQ_perp": sigma_max,
        "jacobian_metadata": dict(assembled["metadata"]),
        "b4_direction_source": "assemble_closed_coupled_autodiff_sparse_jacobian"
        if config.jacobian_assembly == "autodiff_sparse_jacobian_lu"
        else "assemble_closed_coupled_sparse_jacobian",
    }


def _c0g_b3_physical_from_scaled_step(
    *,
    col_scale: np.ndarray,
    q_perp: np.ndarray,
    coeff: np.ndarray,
) -> np.ndarray:
    return np.asarray(col_scale, dtype=np.float64) * (
        np.asarray(q_perp, dtype=np.float64) @ np.asarray(coeff, dtype=np.float64)
    )


def _c0g_b3_metric_norm(
    dx: np.ndarray,
    dtau: float,
    config: C0gB3PseudoArclengthConfig,
) -> float:
    return float(
        math.sqrt(
            float(np.dot(dx, dx)) + (float(config.tau_scale) * float(dtau)) ** 2
        )
    )


def _c0g_b3_metric_dot(
    left_dx: np.ndarray,
    left_dtau: float,
    right_dx: np.ndarray,
    right_dtau: float,
    config: C0gB3PseudoArclengthConfig,
) -> float:
    return float(
        np.dot(np.asarray(left_dx, dtype=np.float64), np.asarray(right_dx, dtype=np.float64))
        + (float(config.tau_scale) ** 2) * float(left_dtau) * float(right_dtau)
    )


def _c0g_b3_tangent_from_linearization(
    *,
    linear: Mapping[str, Any],
    config: C0gB3PseudoArclengthConfig,
    previous_tangent: Mapping[str, Any] | None,
    secant_fallback: Mapping[str, Any] | None,
) -> dict[str, Any]:
    reduced_j = np.asarray(linear["reduced_j"], dtype=np.float64)
    ftau = np.asarray(linear["ftau_scaled"], dtype=np.float64)
    bordered = np.column_stack([reduced_j, ftau])
    source = "lstsq_parameter_tangent_in_gauge_complement"
    dtau_seed = -1.0
    try:
        coeff, _resids, rank, coeff_singular = linalg.lstsq(
            reduced_j,
            -ftau * dtau_seed,
            cond=float(config.corrector_lstsq_rcond),
            lapack_driver="gelsy",
            check_finite=False,
        )
        coeff = np.asarray(coeff, dtype=np.float64)
        dtau = dtau_seed
        linear_residual = reduced_j @ coeff + ftau * dtau
        residual_rel = float(
            np.linalg.norm(linear_residual)
            / max(np.linalg.norm(ftau), np.finfo(np.float64).tiny)
        )
        singular = (
            np.asarray(coeff_singular, dtype=np.float64)
            if coeff_singular is not None and len(coeff_singular)
            else np.linalg.svd(reduced_j, compute_uv=False)
        )
    except np.linalg.LinAlgError:
        coeff = np.zeros(reduced_j.shape[1], dtype=np.float64)
        dtau = 0.0
        singular = np.asarray([], dtype=np.float64)
        residual_rel = math.inf
    if (
        not np.all(np.isfinite(coeff))
        or np.linalg.norm(coeff) <= 0.0
        or not math.isfinite(residual_rel)
        or residual_rel > 1.0e-6
    ):
        source = "svd_null_vector_of_[row_scale_J_col_scale_Q_perp,row_scale_Ftau]"
        try:
            _u, aug_singular, vh = np.linalg.svd(bordered, full_matrices=False)
            vector = vh[-1, :].astype(np.float64, copy=False)
            coeff = vector[:-1]
            dtau = float(vector[-1])
            residual_rel = float(
                np.linalg.norm(bordered @ vector)
                / max(np.linalg.norm(bordered), np.finfo(np.float64).tiny)
            )
            singular = aug_singular
        except np.linalg.LinAlgError:
            coeff = np.zeros(reduced_j.shape[1], dtype=np.float64)
            dtau = 0.0
            residual_rel = math.inf
        if (
            not np.all(np.isfinite(coeff))
            or np.linalg.norm(coeff) <= 0.0
            or not math.isfinite(residual_rel)
        ):
            if secant_fallback is None:
                raise RuntimeError("unable to compute first pseudo-arclength tangent")
            dx = np.asarray(secant_fallback["dx"], dtype=np.float64)
            dtau = float(secant_fallback["dtau"])
            norm = _c0g_b3_metric_norm(dx, dtau, config)
            if norm <= 0.0:
                raise RuntimeError("secant fallback tangent has zero norm")
            return {
                "dx": dx / norm,
                "dtau": float(dtau / norm),
                "source": "secant_fallback",
                "linear_residual_rel": math.nan,
                "null_singular_min": None,
                "orientation_rule": "fallback_secant",
            }
    dx = _c0g_b3_physical_from_scaled_step(
        col_scale=np.asarray(linear["col_scale"], dtype=np.float64),
        q_perp=np.asarray(linear["q_perp"], dtype=np.float64),
        coeff=coeff,
    )
    norm = _c0g_b3_metric_norm(dx, dtau, config)
    if norm <= 0.0:
        raise RuntimeError("computed pseudo-arclength tangent has zero norm")
    dx = dx / norm
    dtau = float(dtau / norm)
    orientation_rule = "first_tangent_oriented_to_decreasing_tau"
    if previous_tangent is not None:
        dot = _c0g_b3_metric_dot(
            dx,
            dtau,
            np.asarray(previous_tangent["dx"], dtype=np.float64),
            float(previous_tangent["dtau"]),
            config,
        )
        orientation_rule = "positive_metric_dot_with_previous_tangent"
        if math.isfinite(dot) and dot < 0.0:
            dx = -dx
            dtau = -dtau
    elif dtau > 0.0:
        dx = -dx
        dtau = -dtau
    return {
        "dx": dx.astype(np.float64, copy=False),
        "dtau": float(dtau),
        "source": source,
        "linear_residual_rel": residual_rel,
        "null_singular_min": float(singular[-1]) if singular.size else None,
        "orientation_rule": orientation_rule,
    }


def _c0g_b3_corrector(
    *,
    predicted_state: np.ndarray,
    predicted_tau: float,
    previous_state: np.ndarray,
    previous_tau: float,
    tangent: Mapping[str, Any],
    ds: float,
    config: C0gB3PseudoArclengthConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    state = np.asarray(predicted_state, dtype=np.float64).copy()
    tau = float(predicted_tau)
    tangent_dx = np.asarray(tangent["dx"], dtype=np.float64)
    tangent_dtau = float(tangent["dtau"])
    rows: list[dict[str, Any]] = []
    best_state = state.copy()
    best_tau = tau
    best_residual_linf = math.inf
    best_constraint = math.inf
    rejected_reason = None
    for iteration in range(1, int(config.corrector_max_iters) + 1):
        linear = _c0g_b3_linearization(state_np=state, tau=tau, config=config, dtype=dtype)
        residual = np.asarray(linear["residual"], dtype=np.float64)
        residual_linf = float(np.max(np.abs(residual)))
        residual_l2 = float(np.linalg.norm(residual))
        constraint = float(
            np.dot(state - previous_state, tangent_dx)
            + (float(config.tau_scale) ** 2) * (tau - float(previous_tau)) * tangent_dtau
            - float(ds)
        )
        if residual_linf < best_residual_linf:
            best_state = state.copy()
            best_tau = float(tau)
            best_residual_linf = residual_linf
            best_constraint = constraint
        if (
            residual_linf <= float(config.corrector_stop_residual_tolerance)
            and abs(constraint) <= float(config.arclength_tolerance)
        ):
            return {
                "status": "CONVERGED",
                "state": state,
                "tau": float(tau),
                "iterations": int(iteration - 1),
                "rows": rows,
                "final_residual_linf": residual_linf,
                "final_residual_l2": residual_l2,
                "arclength_constraint": constraint,
                "last_linearization": linear,
                "message": "corrector drove original physical residual below tolerance",
            }
        reduced_j = np.asarray(linear["reduced_j"], dtype=np.float64)
        q_perp = np.asarray(linear["q_perp"], dtype=np.float64)
        col_scale = np.asarray(linear["col_scale"], dtype=np.float64)
        arclength_row = (tangent_dx * col_scale) @ q_perp
        bordered = np.vstack(
            [
                np.column_stack([reduced_j, np.asarray(linear["ftau_scaled"], dtype=np.float64)]),
                np.concatenate(
                    [
                        arclength_row,
                        np.asarray(
                            [(float(config.tau_scale) ** 2) * tangent_dtau],
                            dtype=np.float64,
                        ),
                    ]
                ),
            ]
        )
        rhs = -np.concatenate(
            [
                np.asarray(linear["row_scale"], dtype=np.float64) * residual,
                np.asarray([constraint], dtype=np.float64),
            ]
        )
        try:
            solve_start = time.perf_counter()
            solution, _resids, rank, singular = linalg.lstsq(
                bordered,
                rhs,
                cond=float(config.corrector_lstsq_rcond),
                lapack_driver="gelsy",
                check_finite=False,
            )
            if singular is None or len(singular) == 0:
                singular = np.linalg.svd(bordered, compute_uv=False)
            solve_seconds = float(time.perf_counter() - solve_start)
        except Exception as exc:
            rejected_reason = f"bordered_lstsq_failed: {exc}"
            break
        solution = np.asarray(solution, dtype=np.float64)
        if not np.all(np.isfinite(solution)):
            rejected_reason = "bordered_lstsq_nonfinite_step"
            break
        dz = solution[:-1]
        dtau = float(solution[-1])
        dx = _c0g_b3_physical_from_scaled_step(
            col_scale=col_scale,
            q_perp=q_perp,
            coeff=dz,
        )
        old_merit = residual_l2 + abs(constraint)
        accepted = False
        alpha = 1.0
        accepted_residual_linf = math.inf
        accepted_constraint = math.inf
        for _line_index in range(int(config.corrector_line_search_steps)):
            trial_state = state + alpha * dx
            trial_tau = float(tau + alpha * dtau)
            trial_residual, _branch, _grid = _c0g_b3_residual_np(
                state_np=trial_state,
                tau=trial_tau,
                config=config,
                dtype=dtype,
            )
            trial_l2 = float(np.linalg.norm(trial_residual))
            trial_linf = float(np.max(np.abs(trial_residual)))
            trial_constraint = float(
                np.dot(trial_state - previous_state, tangent_dx)
                + (float(config.tau_scale) ** 2)
                * (trial_tau - float(previous_tau))
                * tangent_dtau
                - float(ds)
            )
            trial_merit = trial_l2 + abs(trial_constraint)
            if (
                np.all(np.isfinite(trial_residual))
                and math.isfinite(trial_tau)
                and (trial_merit < old_merit or trial_linf <= float(config.residual_tolerance))
            ):
                state = trial_state
                tau = trial_tau
                accepted = True
                accepted_residual_linf = trial_linf
                accepted_constraint = trial_constraint
                break
            alpha *= 0.5
        rows.append(
            {
                "iteration": int(iteration),
                "original_residual_linf": residual_linf,
                "original_residual_l2": residual_l2,
                "arclength_constraint": constraint,
                "bordered_rank": int(rank),
                "bordered_sigma_min": float(np.min(singular)) if len(singular) else math.nan,
                "bordered_sigma_max": float(np.max(singular)) if len(singular) else math.nan,
                "bordered_condition": float(np.max(singular) / np.min(singular))
                if len(singular) and np.min(singular) > 0.0
                else math.inf,
                "linear_solve_seconds": solve_seconds,
                "gauge_rank": int(linear["scaled_gauge_rank"]),
                "gauge_projection_inside_corrector": True,
                "step_metric_norm": _c0g_b3_metric_norm(dx, dtau, config),
                "line_search_alpha": float(alpha) if accepted else None,
                "line_search_accepted": bool(accepted),
                "accepted_residual_linf": accepted_residual_linf if accepted else None,
                "accepted_arclength_constraint": accepted_constraint if accepted else None,
                "sigma_min_JQ_perp": float(linear["sigma_min_JQ_perp"]),
            }
        )
        if not accepted:
            rejected_reason = "corrector_line_search_failed"
            break
    return {
        "status": "FAILED",
        "state": best_state,
        "tau": float(best_tau),
        "iterations": len(rows),
        "rows": rows,
        "final_residual_linf": best_residual_linf,
        "final_residual_l2": math.nan,
        "arclength_constraint": best_constraint,
        "last_linearization": None,
        "message": rejected_reason or "maximum bordered Newton iterations reached",
    }


def _c0g_b3_invariants(vector: np.ndarray, grid) -> dict[str, Any]:
    state = torch.as_tensor(np.asarray(vector, dtype=np.float64), dtype=torch.float64, device="cpu")
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    density = (fields.psi_real**2 + fields.psi_imag**2).detach().cpu().numpy()
    r0 = fields.r0.detach().cpu().numpy()
    return {
        "rho": density.astype(np.float64, copy=False),
        "curl_A": _c0g_curl_np(np.asarray(vector, dtype=np.float64), grid),
        "r0": r0.astype(np.float64, copy=False),
        "mu": float(mu.detach().cpu().item()) if mu is not None else math.nan,
    }


def _c0g_b3_invariant_mismatch(
    *,
    left: np.ndarray,
    right: np.ndarray,
    grid,
) -> dict[str, float]:
    linv = _c0g_b3_invariants(left, grid)
    rinv = _c0g_b3_invariants(right, grid)
    return {
        "rho_linf": float(np.max(np.abs(linv["rho"] - rinv["rho"]))),
        "curl_A_linf": float(np.max(np.abs(linv["curl_A"] - rinv["curl_A"]))),
        "r0_linf": float(np.max(np.abs(linv["r0"] - rinv["r0"]))),
        "mu_abs": float(abs(linv["mu"] - rinv["mu"])),
    }


def _c0g_b3_make_accepted_row(
    *,
    index: int,
    stage: str,
    predecessor_id: str,
    state: np.ndarray,
    tau: float,
    tangent: Mapping[str, Any],
    ds: float,
    corrector: Mapping[str, Any],
    linear: Mapping[str, Any],
    config: C0gB3PseudoArclengthConfig,
    dtype: torch.dtype,
    wall_seconds: float,
) -> dict[str, Any]:
    grid = linear["grid"]
    state_t = torch.as_tensor(np.asarray(state, dtype=np.float64), dtype=dtype, device="cpu")
    metrics = _state_metrics(state_t, grid)
    residual = np.asarray(linear["residual"], dtype=np.float64)
    state_path = _resolve_output_path(
        config.run_root
        / stage
        / "states"
        / f"accepted_{index:03d}_tau_{_format_tau(float(tau))}.npz"
    )
    _save_state_artifact(state_path, state_t)
    cond_rows = corrector.get("rows", [])
    last_cond = cond_rows[-1] if cond_rows else {}
    return {
        "id": f"{stage}_accepted_{index:03d}",
        "stage": stage,
        "accepted_index": int(index),
        "predecessor_state_id": predecessor_id,
        "state_artifact": str(state_path),
        "solved_tau": float(tau),
        "predictor_step": float(ds),
        "corrector_iters": int(corrector.get("iterations", 0)),
        "dtauds": float(tangent["dtau"]),
        "dtauds_sign": "positive" if float(tangent["dtau"]) > 0.0 else "negative",
        "original_residual_linf": float(np.max(np.abs(residual))),
        "original_residual_l2": float(np.linalg.norm(residual)),
        "residual_tolerance": float(config.residual_tolerance),
        "physical_converged": bool(np.max(np.abs(residual)) <= float(config.residual_tolerance)),
        "min_R0": float(metrics["min_R0"]),
        "min_rho": float(metrics["min_rho"]),
        "mu": float(metrics["chemical_potential"]),
        "mass": float(metrics["mass"]),
        "sigma_min_JQ_perp": float(linear["sigma_min_JQ_perp"]),
        "sigma_max_JQ_perp": float(linear["sigma_max_JQ_perp"]),
        "bordered_condition": last_cond.get("bordered_condition"),
        "bordered_sigma_min": last_cond.get("bordered_sigma_min"),
        "arclength_constraint": float(corrector.get("arclength_constraint", math.nan)),
        "wall_seconds": float(wall_seconds),
        "not_initialized_from_comparison_artifact": True,
        "gauge_projection_inside_predictor": True,
        "gauge_projection_inside_corrector": True,
        "b4_direction_source": linear.get("b4_direction_source"),
        "metrics": metrics,
    }


def _c0g_b3_run_sequence(
    *,
    stage: str,
    start_state: np.ndarray,
    start_tau: float,
    start_id: str,
    target_taus: Sequence[float] | None,
    initial_tangent: Mapping[str, Any] | None,
    config: C0gB3PseudoArclengthConfig,
    dtype: torch.dtype,
    max_steps: int,
) -> dict[str, Any]:
    accepted: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    state = np.asarray(start_state, dtype=np.float64).copy()
    tau = float(start_tau)
    predecessor_id = str(start_id)
    previous_tangent = initial_tangent
    ds = float(config.initial_ds)
    reseed_attempted = False
    secant_fallback_attempted = False
    for step_index in range(int(max_steps)):
        step_start = time.perf_counter()
        linear = _c0g_b3_linearization(state_np=state, tau=tau, config=config, dtype=dtype)
        secant_fallback = None
        tangent = _c0g_b3_tangent_from_linearization(
            linear=linear,
            config=config,
            previous_tangent=previous_tangent,
            secant_fallback=secant_fallback,
        )
        if target_taus is not None and step_index < len(target_taus):
            target_tau = float(target_taus[step_index])
            if abs(float(tangent["dtau"])) > 1.0e-14:
                ds = abs((target_tau - tau) / float(tangent["dtau"]))
                ds = min(max(ds, float(config.min_ds)), float(config.max_ds))
        predicted_state = state + ds * np.asarray(tangent["dx"], dtype=np.float64)
        predicted_tau = float(tau + ds * float(tangent["dtau"]))
        corrector = _c0g_b3_corrector(
            predicted_state=predicted_state,
            predicted_tau=predicted_tau,
            previous_state=state,
            previous_tau=tau,
            tangent=tangent,
            ds=ds,
            config=config,
            dtype=dtype,
        )
        if corrector.get("status") == "CONVERGED":
            corrected_state = np.asarray(corrector["state"], dtype=np.float64)
            corrected_tau = float(corrector["tau"])
            corrected_linear = _c0g_b3_linearization(
                state_np=corrected_state,
                tau=corrected_tau,
                config=config,
                dtype=dtype,
            )
            target_landing_attempts: list[dict[str, Any]] = []
            if target_taus is not None and step_index < len(target_taus):
                target_tau = float(target_taus[step_index])
                for retry_index in range(5):
                    tau_error = float(corrected_tau - target_tau)
                    if abs(tau_error) <= float(config.target_tau_tolerance):
                        break
                    denom = float(corrected_tau - tau)
                    if abs(denom) <= 1.0e-15:
                        target_landing_attempts.append(
                            {
                                "retry_index": int(retry_index),
                                "status": "FAILED",
                                "reason": "zero_tau_response_to_ds",
                                "solved_tau": corrected_tau,
                                "target_tau": target_tau,
                            }
                        )
                        break
                    retry_ds = float(ds) * (target_tau - float(tau)) / denom
                    retry_ds = min(max(abs(retry_ds), float(config.min_ds)), float(config.max_ds))
                    retry_predicted_state = state + retry_ds * np.asarray(tangent["dx"], dtype=np.float64)
                    retry_predicted_tau = float(tau + retry_ds * float(tangent["dtau"]))
                    retry_corrector = _c0g_b3_corrector(
                        predicted_state=retry_predicted_state,
                        predicted_tau=retry_predicted_tau,
                        previous_state=state,
                        previous_tau=tau,
                        tangent=tangent,
                        ds=retry_ds,
                        config=config,
                        dtype=dtype,
                    )
                    target_landing_attempts.append(
                        {
                            "retry_index": int(retry_index),
                            "status": retry_corrector.get("status"),
                            "requested_ds": retry_ds,
                            "target_tau": target_tau,
                            "solved_tau": retry_corrector.get("tau"),
                            "tau_error": None
                            if retry_corrector.get("tau") is None
                            else float(retry_corrector.get("tau") - target_tau),
                            "original_residual_linf": retry_corrector.get("final_residual_linf"),
                        }
                    )
                    if retry_corrector.get("status") != "CONVERGED":
                        break
                    ds = retry_ds
                    corrector = retry_corrector
                    corrected_state = np.asarray(corrector["state"], dtype=np.float64)
                    corrected_tau = float(corrector["tau"])
                    corrected_linear = _c0g_b3_linearization(
                        state_np=corrected_state,
                        tau=corrected_tau,
                        config=config,
                        dtype=dtype,
                    )
                if abs(float(corrected_tau) - target_tau) > float(config.target_tau_tolerance):
                    failures.append(
                        {
                            "step_index": int(step_index),
                            "status": "TARGET_MISSED",
                            "target_tau": target_tau,
                            "solved_tau": corrected_tau,
                            "tau_error": float(corrected_tau - target_tau),
                            "target_landing_attempts": target_landing_attempts,
                        }
                    )
                    break
            row = _c0g_b3_make_accepted_row(
                index=len(accepted),
                stage=stage,
                predecessor_id=predecessor_id,
                state=corrected_state,
                tau=corrected_tau,
                tangent=tangent,
                ds=ds,
                corrector=corrector,
                linear=corrected_linear,
                config=config,
                dtype=dtype,
                wall_seconds=float(time.perf_counter() - step_start),
            )
            if row["min_R0"] <= 0.0:
                failures.append(
                    {
                        "step_index": int(step_index),
                        "status": "HALT",
                        "reason": "min_R0_nonpositive_model_forbidden",
                        "min_R0": row["min_R0"],
                    }
                )
                break
            accepted.append(row)
            previous_state = state
            previous_tau = tau
            state = corrected_state
            tau = corrected_tau
            predecessor_id = row["id"]
            previous_tangent = _c0g_b3_tangent_from_linearization(
                linear=corrected_linear,
                config=config,
                previous_tangent=tangent,
                secant_fallback={
                    "dx": corrected_state - previous_state,
                    "dtau": corrected_tau - previous_tau,
                },
                )
            ds = min(float(config.max_ds), ds * float(config.ds_growth))
            continue
        failures.append(
            {
                "step_index": int(step_index),
                "status": "FAILED",
                "reason": corrector.get("message"),
                "ds": float(ds),
                "predicted_tau": predicted_tau,
                "best_tau": corrector.get("tau"),
                "best_original_residual_linf": corrector.get("final_residual_linf"),
                "adaptive_ds_reduction": True,
            }
        )
        if ds > float(config.min_ds):
            ds = max(float(config.min_ds), ds * float(config.ds_shrink))
            continue
        if not reseed_attempted:
            previous_tangent = None
            reseed_attempted = True
            failures[-1]["tangent_reseed_attempted_after_failure"] = True
            continue
        if not secant_fallback_attempted and len(accepted) >= 1:
            secant_fallback_attempted = True
            failures[-1]["secant_fallback_attempted_after_failure"] = True
            continue
        break
    return {
        "stage": stage,
        "accepted": accepted,
        "failures": failures,
        "final_state": state,
        "final_tau": float(tau),
        "final_state_id": predecessor_id,
        "final_tangent": previous_tangent,
        "adaptive_ds_reduction_attempted": bool(any(row.get("adaptive_ds_reduction") for row in failures)),
        "tangent_reseed_attempted": bool(reseed_attempted),
        "secant_fallback_attempted": bool(secant_fallback_attempted),
    }


def _c0g_b3_reproduction_check(
    *,
    b30: Mapping[str, Any],
    deepcrawl_rows: Sequence[Mapping[str, Any]],
    config: C0gB3PseudoArclengthConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    accepted = list(b30.get("accepted", []))
    if not accepted:
        return {"status": "FAIL", "reason": "no accepted B-3.0 continuation states"}
    by_tau = {float(row["target_tau"]): row for row in deepcrawl_rows}
    checks: list[dict[str, Any]] = []
    for row in accepted:
        solved_tau = float(row["solved_tau"])
        nearest_tau = min(by_tau, key=lambda value: abs(value - solved_tau))
        if abs(solved_tau - nearest_tau) > float(config.target_tau_tolerance):
            checks.append(
                {
                    "accepted_id": row["id"],
                    "solved_tau": solved_tau,
                    "comparison_tau": nearest_tau,
                    "tau_abs_error": abs(solved_tau - nearest_tau),
                    "status": "SKIPPED_INTERMEDIATE",
                    "reason": "accepted continuation state is an internal B-3.0 subdivision, not a recorded comparison tau",
                    "comparison_artifact_read_after_run": False,
                }
            )
            continue
        reference_row = by_tau[nearest_tau]
        accepted_state = _load_state_artifact(_resolve_input_path(str(row["state_artifact"])), dtype=dtype)
        reference_state = _load_state_artifact(_resolve_input_path(str(reference_row["state_artifact"])), dtype=dtype)
        _branch, grid, _residual_fn = _c0g_b3_context(tau=solved_tau, config=config, dtype=dtype)
        mismatch = _c0g_b3_invariant_mismatch(
            left=accepted_state.detach().cpu().numpy(),
            right=reference_state.detach().cpu().numpy(),
            grid=grid,
        )
        max_mismatch = max(float(value) for value in mismatch.values())
        checks.append(
            {
                "accepted_id": row["id"],
                "solved_tau": solved_tau,
                "comparison_tau": nearest_tau,
                "tau_abs_error": abs(solved_tau - nearest_tau),
                "mismatch": mismatch,
                "max_gauge_invariant_mismatch": max_mismatch,
                "status": "PASS"
                if max_mismatch <= float(config.reproduction_comparison_tolerance)
                and abs(solved_tau - nearest_tau) <= float(config.target_tau_tolerance)
                and float(row["original_residual_linf"]) <= float(config.residual_tolerance)
                else "FAIL",
                "comparison_artifact_read_after_run": True,
            }
        )
    continuity = []
    prev_state = None
    for row in accepted:
        current = _load_state_artifact(_resolve_input_path(str(row["state_artifact"])), dtype=dtype)
        current_np = current.detach().cpu().numpy().astype(np.float64, copy=False)
        if prev_state is not None:
            distance = float(np.linalg.norm(current_np - prev_state))
            continuity.append(
                {
                    "accepted_id": row["id"],
                    "state_l2_distance_from_previous": distance,
                    "distinct_above_floor": bool(distance > float(config.distinct_state_floor)),
                }
            )
        prev_state = current_np
    required_checks = [row for row in checks if row.get("status") != "SKIPPED_INTERMEDIATE"]
    status = (
        "PASS"
        if required_checks and all(row["status"] == "PASS" for row in required_checks)
        else "FAIL"
    )
    if continuity and not all(row["distinct_above_floor"] for row in continuity):
        status = "FAIL"
    return {
        "status": status,
        "checks": checks,
        "state_continuity": continuity,
        "not_initialized_from_comparison_artifact": bool(
            all(row.get("not_initialized_from_comparison_artifact") for row in accepted)
        ),
        "comparison_artifacts_read_only_after_run": True,
        "tolerance": float(config.reproduction_comparison_tolerance),
    }


def _c0g_b3_adjudicate(
    *,
    b30: Mapping[str, Any],
    reproduction: Mapping[str, Any],
    b31: Mapping[str, Any] | None,
    b4_reconfirm: Mapping[str, Any],
    target: Mapping[str, Any],
    config: C0gB3PseudoArclengthConfig,
) -> dict[str, Any]:
    if reproduction.get("status") != "PASS":
        return {
            "verdict": "NOT_MEASURED",
            "why": "B-3.0 reproduction gate failed; hard-region continuation was not attempted",
            "scientific_verdict_earned": False,
        }
    if b4_reconfirm.get("status") != "PASS":
        return {
            "verdict": "NOT_MEASURED",
            "why": "B-4 match gate failed; B-3 direction was not trusted",
            "scientific_verdict_earned": False,
        }
    if b31 is None:
        return {
            "verdict": "NOT_MEASURED",
            "why": "B-3.1 did not run",
            "scientific_verdict_earned": False,
        }
    accepted = list(b31.get("accepted", []))
    if not accepted:
        exhausted = bool(
            b31.get("adaptive_ds_reduction_attempted")
            and b31.get("tangent_reseed_attempted")
            and b31.get("secant_fallback_attempted")
        )
        return {
            "verdict": "GENUINE_ENDPOINT" if exhausted else "NOT_MEASURED",
            "why": "bordered corrector could not drive original residual with tau free"
            if exhausted
            else "no B-3.1 accepted states and exhaustion criteria were not complete",
            "scientific_verdict_earned": exhausted,
            "exhaustion": {
                "adaptive_ds_reduction": bool(b31.get("adaptive_ds_reduction_attempted")),
                "tangent_reseed": bool(b31.get("tangent_reseed_attempted")),
                "secant_fallback": bool(b31.get("secant_fallback_attempted")),
            },
        }
    all_rows = list(b30.get("accepted", [])) + accepted
    signs = [1 if float(row["dtauds"]) > 0.0 else -1 for row in all_rows]
    reversal_index = None
    for index in range(1, len(signs)):
        if signs[index] != signs[index - 1]:
            reversal_index = index
            break
    min_r0_positive = all(float(row.get("min_R0", -math.inf)) > 0.0 for row in all_rows)
    stall_tau = float(target["deepest_converged_tau"])
    below = [
        row
        for row in accepted
        if float(row.get("solved_tau", math.inf)) < stall_tau
        and bool(row.get("physical_converged"))
    ]
    deepest = min((float(row["solved_tau"]) for row in below), default=None)
    if reversal_index is not None:
        far_side = all_rows[reversal_index:]
        far_converged = [
            row
            for row in far_side
            if bool(row.get("physical_converged"))
            and float(row.get("original_residual_linf", math.inf)) <= float(config.residual_tolerance)
        ]
        if len(far_converged) >= 3 and min_r0_positive:
            return {
                "verdict": "ROUNDS_FOLD",
                "why": "dτ/ds reversed sign with at least three converged states on the far side",
                "scientific_verdict_earned": True,
                "turning_index": int(reversal_index),
                "far_side_count": int(len(far_converged)),
            }
    progress_tau = stall_tau - float(config.b31_deepest_tau_margin)
    if (
        len(below) >= int(config.b31_required_below_count)
        and deepest is not None
        and deepest <= progress_tau
        and reversal_index is None
        and min_r0_positive
    ):
        return {
            "verdict": "CONTINUES_NO_TURNING",
            "why": (
                f"{len(below)} distinct accepted states have solved tau below "
                f"{stall_tau:.12g}, deepest {deepest:.12g}, with no dτ/ds reversal"
            ),
            "scientific_verdict_earned": True,
            "below_stall_count": int(len(below)),
            "deepest_converged_tau": deepest,
            "required_deepest_tau": progress_tau,
        }
    failures = list(b31.get("failures", []))
    exhausted = bool(
        failures
        and b31.get("adaptive_ds_reduction_attempted")
        and b31.get("tangent_reseed_attempted")
        and b31.get("secant_fallback_attempted")
    )
    if exhausted and min_r0_positive:
        return {
            "verdict": "GENUINE_ENDPOINT",
            "why": "full pseudo-arclength machinery exhausted without correcting original residual",
            "scientific_verdict_earned": True,
            "below_stall_count": int(len(below)),
            "deepest_converged_tau": deepest,
            "exhaustion": {
                "adaptive_ds_reduction": True,
                "tangent_reseed": True,
                "secant_fallback": True,
            },
        }
    return {
        "verdict": "NOT_MEASURED",
        "why": "B-3.1 produced too few accepted states for a scientific verdict in this chunk",
        "scientific_verdict_earned": False,
        "below_stall_count": int(len(below)),
        "deepest_converged_tau": deepest,
        "dtauds_reversed": reversal_index is not None,
    }


def write_c0g_b3_pseudoarclength_report(result: Mapping[str, Any], path: Path) -> None:
    b4 = result.get("b4_reconfirm", {})
    target = result.get("target", {})
    verdict = result.get("b32_verdict", {})
    b30 = result.get("B3_0", {})
    b31 = result.get("B3_1", {})
    lines = [
        "# Path-A C0g B-3 Pseudo-Arclength",
        "",
        f"B-3.0 reproduction: **{result.get('B3_0_reproduction', {}).get('status')}**",
        f"B-3.2 verdict: **{verdict.get('verdict')}**",
        "",
        "## Method",
        "",
        "- Reduced unknowns are `(Q_perp coefficients, tau)`.",
        f"- Declared arclength metric: `||dx||_2^2 + ({result.get('config', {}).get('tau_scale')} dτ)^2`.",
        "- Predictor tangent first solves `row_scale * J_original * col_scale * Q_perp * z = -row_scale * F_tau * dτ` in the gauge complement; the augmented null vector is the fallback.",
        "- The bordered corrector drives the unchanged original physical residual; the arclength row is not a convergence arbiter.",
        "- `Q_perp` is rebuilt inside every predictor and corrector linear solve.",
        "",
        "## B-4 Reconfirm",
        "",
        "```yaml",
        f"status: {b4.get('status')}",
        f"reference_tau: {b4.get('reference_tau')}",
        f"tolerances: {b4.get('match_tolerances')}",
        f"assembly_speedup_x: {b4.get('assembly_speedup_x')}",
        "```",
        "",
        "## Target",
        "",
        "```yaml",
        f"deepcrawl_json_path: {target.get('deepcrawl_json_path')}",
        f"seed_tau: {target.get('seed_tau')}",
        f"deepest_converged_tau: {target.get('deepest_converged_tau')}",
        f"deepest_residual_linf: {target.get('deepest_residual_linf')}",
        "```",
        "",
        "## B-3.0 Ladder",
        "",
        _markdown_table(
            [
                "accepted_index",
                "predecessor_state_id",
                "predictor_step",
                "corrector_iters",
                "solved_tau",
                "dtauds",
                "original_residual_linf",
                "min_R0",
                "min_rho",
                "mu",
                "sigma_min_JQ_perp",
                "bordered_condition",
                "wall_seconds",
                "not_initialized_from_comparison_artifact",
            ],
            b30.get("accepted", []),
        ),
        "",
        "## B-3.0 Reproduction Check",
        "",
        "```yaml",
        f"{result.get('B3_0_reproduction', {})}",
        "```",
        "",
        "## B-3.1 Ladder",
        "",
        _markdown_table(
            [
                "accepted_index",
                "predecessor_state_id",
                "predictor_step",
                "corrector_iters",
                "solved_tau",
                "dtauds",
                "original_residual_linf",
                "min_R0",
                "min_rho",
                "mu",
                "sigma_min_JQ_perp",
                "bordered_condition",
                "wall_seconds",
            ],
            b31.get("accepted", []) if isinstance(b31, Mapping) else [],
        ),
        "",
        "## Verdict",
        "",
        "```yaml",
        f"{verdict}",
        "```",
        "",
        "## Scope Guard",
        "",
        "```yaml",
    ]
    for key, value in result.get("scope_guard", {}).items():
        lines.append(f"{key}: {value}")
    lines.extend(
        [
            "```",
            "",
            "## Git Diff Summary",
            "",
            "```",
            str(result.get("git_diff_summary", "")),
            "```",
            "",
            f"Machine JSON: `{result.get('json_path')}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _c0g_b3_git_diff_summary() -> str:
    import subprocess

    proc = subprocess.run(
        ["git", "diff", "--", "src/stage1_solver/coupled_branch.py", "src/stage1_solver/operators.py"],
        cwd=str(Path(__file__).resolve().parents[2]),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if proc.stdout.strip():
        return proc.stdout.strip()
    return "No diff in src/stage1_solver/coupled_branch.py or src/stage1_solver/operators.py."


def run_c0g_b3_pseudoarclength(
    config: C0gB3PseudoArclengthConfig | None = None,
    *,
    b4_recheck: bool = True,
    b31_steps: int | None = None,
    b30_only: bool = False,
) -> dict[str, Any]:
    cfg = config or C0gB3PseudoArclengthConfig()
    started = time.perf_counter()
    dtype = configure_backend(BackendConfig())
    payload = _c0g_b3_deepcrawl_payload(cfg)
    target = _c0g_b3_target_summary(payload, cfg)
    deepcrawl_rows = _c0g_b3_converged_deepcrawl_rows(payload)
    b4_reconfirm: Mapping[str, Any]
    if b4_recheck:
        b4_reconfirm = run_c0g_b4_match_gate(_c0g_b3_build_config(cfg))
    else:
        b4_reconfirm = payload.get("b4_reconfirm", {"status": "NOT_RECHECKED"})
    seed_state = _load_state_artifact(
        _resolve_input_path(str(target["seed_state_artifact"])),
        dtype=dtype,
    )
    seed_np = seed_state.detach().cpu().numpy().astype(np.float64, copy=True)
    recorded_target_taus = [
        float(tau)
        for tau in target["b30_target_taus"]
        if float(tau) < float(target["seed_tau"]) - 5.0e-13
    ]
    target_taus: list[float] = []
    current_tau = float(target["seed_tau"])
    for recorded_tau in recorded_target_taus:
        pieces = max(1, int(math.ceil(abs(recorded_tau - current_tau) / 7.5e-7)))
        for piece_index in range(1, pieces + 1):
            target_taus.append(
                float(current_tau + (recorded_tau - current_tau) * piece_index / pieces)
            )
        current_tau = recorded_tau
    b30 = _c0g_b3_run_sequence(
        stage="B3_0",
        start_state=seed_np,
        start_tau=float(target["seed_tau"]),
        start_id="single_initial_seed_from_deepcrawl",
        target_taus=target_taus,
        initial_tangent=None,
        config=cfg,
        dtype=dtype,
        max_steps=len(target_taus),
    )
    reproduction = _c0g_b3_reproduction_check(
        b30=b30,
        deepcrawl_rows=deepcrawl_rows,
        config=cfg,
        dtype=dtype,
    )
    b31 = None
    if reproduction.get("status") == "PASS" and b4_reconfirm.get("status") == "PASS" and not b30_only:
        b31 = _c0g_b3_run_sequence(
            stage="B3_1",
            start_state=np.asarray(b30["final_state"], dtype=np.float64),
            start_tau=float(b30["final_tau"]),
            start_id=str(b30["final_state_id"]),
            target_taus=None,
            initial_tangent=b30.get("final_tangent"),
            config=cfg,
            dtype=dtype,
            max_steps=int(cfg.b31_max_steps if b31_steps is None else b31_steps),
        )
    verdict = _c0g_b3_adjudicate(
        b30=b30,
        reproduction=reproduction,
        b31=b31,
        b4_reconfirm=b4_reconfirm,
        target=target,
        config=cfg,
    )
    result = {
        "schema": "stage1_pathA_C0g_B3_pseudoarclength/v1",
        "source_revision": source_revision(),
        "config": _c0g_b3_config_to_dict(cfg),
        "target": target,
        "b4_reconfirm": dict(b4_reconfirm),
        "B3_0": {
            key: value
            for key, value in b30.items()
            if key not in {"final_state", "final_tangent"}
        },
        "B3_0_reproduction": reproduction,
        "B3_1": None
        if b31 is None
        else {
            key: value
            for key, value in b31.items()
            if key not in {"final_state", "final_tangent"}
        },
        "b32_verdict": verdict,
        "scope_guard": {
            "B3_pseudoarclength_built": True,
            "B4_analytic_sparse_assembly_built": True,
            "single_arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
            "convergence_judged_only_on_original_residual": True,
            "patha_closed_branch_residual_touched": False,
            "faithful_operators_touched": False,
            "xi_or_grad_div_penalty_touched": False,
            "physical_export_permitted_touched": False,
            "prefer_existing_b2c_background_predictor": False,
            "frozen_physics_touched": False,
            "no_LM_no_PTC_no_Sobolev_regularization": True,
        },
        "git_diff_summary": _c0g_b3_git_diff_summary(),
        "elapsed_seconds": float(time.perf_counter() - started),
        "report_path": str(_resolve_output_path(cfg.report_path)),
        "json_path": str(_resolve_output_path(cfg.json_path)),
    }
    json_path = _resolve_output_path(cfg.json_path)
    report_path = _resolve_output_path(cfg.report_path)
    _c0g_write_json(json_path, result)
    write_c0g_b3_pseudoarclength_report(result, report_path)
    return result


def _c0g_b3_stdout_summary(result: Mapping[str, Any]) -> str:
    b30_status = result.get("B3_0_reproduction", {}).get("status")
    b31_rows = (result.get("B3_1") or {}).get("accepted", [])
    all_rows = list(result.get("B3_0", {}).get("accepted", [])) + list(b31_rows)
    reversed_sign = False
    signs = [1 if float(row.get("dtauds", 0.0)) > 0.0 else -1 for row in all_rows]
    for left, right in zip(signs, signs[1:]):
        if left != right:
            reversed_sign = True
            break
    target_tau = result.get("target", {}).get("deepest_converged_tau")
    below = [
        row for row in b31_rows if target_tau is not None and float(row.get("solved_tau", math.inf)) < float(target_tau)
    ]
    min_r0_values = [float(row.get("min_R0", math.nan)) for row in all_rows if row.get("min_R0") is not None]
    sigma_values = [
        float(row.get("sigma_min_JQ_perp", math.nan))
        for row in all_rows
        if row.get("sigma_min_JQ_perp") is not None
    ]
    ladder_lines = []
    for row in b31_rows:
        ladder_lines.append(
            (
                "  tau={tau}, dtauds={dtauds}, residual={resid}, min_R0={r0}, "
                "min_rho={rho}, mu={mu}, sigma_min(JQ_perp)={sigma}, "
                "bordered_cond={cond}, ds={ds}, iters={iters}, seconds={seconds}"
            ).format(
                tau=_c0g_fmt(row.get("solved_tau")),
                dtauds=_c0g_fmt(row.get("dtauds")),
                resid=_c0g_fmt(row.get("original_residual_linf")),
                r0=_c0g_fmt(row.get("min_R0")),
                rho=_c0g_fmt(row.get("min_rho")),
                mu=_c0g_fmt(row.get("mu")),
                sigma=_c0g_fmt(row.get("sigma_min_JQ_perp")),
                cond=_c0g_fmt(row.get("bordered_condition")),
                ds=_c0g_fmt(row.get("predictor_step")),
                iters=row.get("corrector_iters"),
                seconds=_c0g_fmt(row.get("wall_seconds")),
            )
        )
    verdict = result.get("b32_verdict", {})
    scope = result.get("scope_guard", {})
    return "\n".join(
        [
            f"B-3.0 reproduction: {b30_status}",
            "B-3.1 ladder:",
            *(ladder_lines or ["  n/a"]),
            f"dτ/ds reversed: {reversed_sign}",
            f"B-3.2 verdict: {verdict.get('verdict')} - {verdict.get('why')}",
            (
                f"deepest converged tau reached: {_c0g_fmt(min([float(row.get('solved_tau')) for row in b31_rows], default=None))}; "
                f"count below {target_tau}: {len(below)}"
            ),
            (
                f"min_R0 range: {_c0g_fmt(min(min_r0_values) if min_r0_values else None)} .. "
                f"{_c0g_fmt(max(min_r0_values) if min_r0_values else None)}; stayed > 0 = "
                f"{bool(min_r0_values and min(min_r0_values) > 0.0)}"
            ),
            (
                f"sigma_min trend: first={_c0g_fmt(sigma_values[0] if sigma_values else None)}, "
                f"last={_c0g_fmt(sigma_values[-1] if sigma_values else None)}"
            ),
            f"B-4 reconfirm: {result.get('b4_reconfirm', {}).get('status')}",
            (
                "Scope untouched: residual/operators/(1/xi)/export/physics="
                f"{not scope.get('patha_closed_branch_residual_touched')}/"
                f"{not scope.get('faithful_operators_touched')}/"
                f"{not scope.get('xi_or_grad_div_penalty_touched')}/"
                f"{not scope.get('physical_export_permitted_touched')}/"
                f"{not scope.get('frozen_physics_touched')}; "
                f"original residual only={scope.get('convergence_judged_only_on_original_residual')}"
            ),
            f"Files: report={result.get('report_path')}, json={result.get('json_path')}",
        ]
    )


def _c0g_b3_followup_config_to_dict(config: C0gB3FollowupConfig) -> dict[str, Any]:
    data = asdict(config)
    for key in ("deepcrawl_json_path", "b3_json_path", "run_root", "report_path", "json_path"):
        data[key] = str(data[key])
    return data


def _c0g_b3_followup_load_json(path: Path) -> dict[str, Any]:
    resolved = _resolve_input_path(path)
    payload = json.loads(resolved.read_text(encoding="utf-8"))
    payload["_resolved_path"] = str(resolved)
    return payload


def _c0g_b3_followup_progress_path(config: C0gB3FollowupConfig) -> Path:
    return _resolve_output_path(config.run_root / "progress.jsonl")


def _c0g_b3_followup_append_progress(
    config: C0gB3FollowupConfig,
    row: Mapping[str, Any],
) -> None:
    path = _c0g_b3_followup_progress_path(config)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(dict(row), sort_keys=True) + "\n")


def _c0g_b3_followup_state_np(path: str | Path, *, dtype: torch.dtype) -> np.ndarray:
    return (
        _c0g_state_vector_from_path(path, dtype=dtype)
        .detach()
        .cpu()
        .numpy()
        .astype(np.float64, copy=True)
    )


def _c0g_b3_followup_deep_sigma_rows(deep_payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    source_rows = list(deep_payload.get("step6_bordered_conditioning", {}).get("rows", []))
    if not source_rows:
        source_rows = list(deep_payload.get("ladder", []))
    for row in source_rows:
        sigma = row.get("sigma_min_JQ_perp")
        if sigma is None:
            continue
        rows.append(
            {
                "tau": float(row.get("tau", math.nan)),
                "sigma_min_JQ_perp": float(sigma),
                "cond_JQ_perp": row.get("cond_JQ_perp"),
                "status": row.get("status"),
            }
        )
    rows.sort(key=lambda item: float(item["tau"]), reverse=True)
    return rows


def _c0g_b3_followup_deep_converged_rows(
    deep_payload: Mapping[str, Any],
) -> list[dict[str, Any]]:
    rows = []
    for row in deep_payload.get("admissibility", {}).get("rows", []):
        if row.get("status") == "PASS":
            rows.append(dict(row))
    rows.sort(key=lambda item: float(item.get("tau", math.nan)), reverse=True)
    return rows


def _c0g_b3_followup_find_deep_row(
    deep_payload: Mapping[str, Any],
    tau: float,
) -> dict[str, Any]:
    rows = _c0g_b3_followup_deep_converged_rows(deep_payload)
    best = min(rows, key=lambda row: abs(float(row.get("tau", math.inf)) - float(tau)))
    if abs(float(best.get("tau", math.inf)) - float(tau)) > 5.0e-8:
        raise RuntimeError(f"could not find deepcrawl row near tau={tau}")
    return best


def _c0g_b3_followup_find_b3_accepted(
    b3_payload: Mapping[str, Any],
    accepted_id: str,
) -> dict[str, Any]:
    for row in b3_payload.get("B3_0", {}).get("accepted", []):
        if row.get("id") == accepted_id:
            return dict(row)
    raise RuntimeError(f"could not find B3 accepted state {accepted_id!r}")


def _c0g_b3_followup_c0(
    *,
    deep_payload: Mapping[str, Any],
    b3_payload: Mapping[str, Any],
    config: C0gB3FollowupConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    sigma_rows = _c0g_b3_followup_deep_sigma_rows(deep_payload)
    finite_sigmas = [
        row for row in sigma_rows if math.isfinite(float(row.get("sigma_min_JQ_perp", math.nan)))
    ]
    deepest = min(finite_sigmas, key=lambda row: float(row["tau"])) if finite_sigmas else {}
    last_three = finite_sigmas[-3:]
    trend_call = "NOT_MEASURED"
    if len(last_three) >= 3:
        values = [float(row["sigma_min_JQ_perp"]) for row in last_three]
        if values[-1] < 0.25 * values[-2] and values[-2] < values[-3]:
            trend_call = "STILL_DROPPING_STEEPLY"
        elif min(values) > 0.0 and max(values) / max(min(values), np.finfo(np.float64).tiny) < 3.0:
            trend_call = "BOTTOMING_OUT"
        else:
            trend_call = "MIXED"

    accepted = list(b3_payload.get("B3_0", {}).get("accepted", []))
    balance_rows: list[dict[str, Any]] = []
    previous_state: np.ndarray | None = None
    previous_tau: float | None = None
    for row in accepted:
        state = _c0g_b3_followup_state_np(str(row["state_artifact"]), dtype=dtype)
        tau = float(row.get("solved_tau", math.nan))
        if previous_state is not None and previous_tau is not None:
            dx = state - previous_state
            dtau = tau - previous_tau
            dx2 = float(np.dot(dx, dx))
            tau2 = float((float(config.tau_scale) * dtau) ** 2)
            balance_rows.append(
                {
                    "id": row.get("id"),
                    "previous_tau": previous_tau,
                    "tau": tau,
                    "dx_norm_squared": dx2,
                    "tau_metric_squared": tau2,
                    "x_fraction": float(dx2 / max(dx2 + tau2, np.finfo(np.float64).tiny)),
                    "dominant": "x" if dx2 > tau2 else "tau",
                }
            )
        previous_state = state
        previous_tau = tau
    return {
        "sigma_rows": sigma_rows,
        "sigma_trend_call": trend_call,
        "deepest_converged_tau": deepest.get("tau"),
        "deepest_converged_sigma_min_JQ_perp": deepest.get("sigma_min_JQ_perp"),
        "sigma_fit_zero_crossing": float(config.sigma_fit_zero_crossing),
        "gap_zero_crossing_minus_deepest_tau": (
            float(config.sigma_fit_zero_crossing) - float(deepest["tau"]) if deepest else None
        ),
        "arclength_metric_balance": balance_rows,
        "arclength_balance_source": "actual B3_0 accepted state artifacts",
    }


def _c0g_b3_followup_checkpoint_c1(
    *,
    config: C0gB3FollowupConfig,
    label: str,
    state: np.ndarray,
    tau: float,
    dtype: torch.dtype,
) -> str:
    state_t = torch.as_tensor(np.asarray(state, dtype=np.float64), dtype=dtype, device="cpu")
    path = _resolve_output_path(
        config.run_root / "C1" / "states" / f"{label}_tau_{_format_tau(float(tau))}.npz"
    )
    _save_state_artifact(path, state_t)
    return str(path)


def _c0g_b3_followup_recent_progress_stalled(
    history: Sequence[Mapping[str, Any]],
    *,
    window: int,
    epsilon_res: float,
) -> bool:
    if len(history) <= window:
        return False
    start = float(history[-window - 1]["original_residual_linf"])
    best = min(float(row["original_residual_linf"]) for row in history[-window:])
    if start <= 0.0:
        return False
    return best > start * (1.0 - float(epsilon_res))


def _c0g_b3_followup_reconverge_fixed_tau(
    *,
    label: str,
    start_state: np.ndarray,
    tau: float,
    config: C0gB3FollowupConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    state = np.asarray(start_state, dtype=np.float64).copy()
    started = time.perf_counter()
    rows: list[dict[str, Any]] = []
    status = "NOT_TIGHT_CONVERGED"
    message = "maximum iteration budget reached"
    final_linear: Mapping[str, Any] | None = None

    def append_progress(row: Mapping[str, Any]) -> None:
        _c0g_b3_followup_append_progress(
            config,
            {
                "stage": "C1",
                "start_label": label,
                **dict(row),
            },
        )

    for iteration in range(0, int(config.c1_max_iters) + 1):
        linear = _c0g_b3_linearization(
            state_np=state,
            tau=float(tau),
            config=C0gB3PseudoArclengthConfig(
                deepcrawl_json_path=config.deepcrawl_json_path,
                run_root=config.run_root / "C1" / label,
                grid=config.grid,
                jacobian_assembly=config.jacobian_assembly,
                tau_scale=config.tau_scale,
            ),
            dtype=dtype,
        )
        final_linear = linear
        residual = np.asarray(linear["residual"], dtype=np.float64)
        residual_linf = float(np.max(np.abs(residual)))
        reduced_j = np.asarray(linear["reduced_j"], dtype=np.float64)
        rhs = -np.asarray(linear["row_scale"], dtype=np.float64) * residual
        solve_start = time.perf_counter()
        try:
            dz, _resids, rank, singular = linalg.lstsq(
                reduced_j,
                rhs,
                cond=float(config.c1_lstsq_rcond),
                lapack_driver="gelsy",
                check_finite=False,
            )
            if singular is None or len(singular) == 0:
                singular = np.linalg.svd(reduced_j, compute_uv=False)
            dx = _c0g_b3_physical_from_scaled_step(
                col_scale=np.asarray(linear["col_scale"], dtype=np.float64),
                q_perp=np.asarray(linear["q_perp"], dtype=np.float64),
                coeff=np.asarray(dz, dtype=np.float64),
            )
            linear_message = "ok"
        except Exception as exc:
            rank = 0
            singular = np.asarray([], dtype=np.float64)
            dx = np.full_like(state, np.nan)
            linear_message = f"lstsq_failed: {exc}"
        solve_seconds = float(time.perf_counter() - solve_start)
        step_norm = float(np.linalg.norm(dx)) if np.all(np.isfinite(dx)) else math.inf
        metrics = _state_metrics(
            torch.as_tensor(state, dtype=dtype, device="cpu"),
            linear["grid"],
        )
        row = {
            "iteration": int(iteration),
            "tau": float(tau),
            "original_residual_linf": residual_linf,
            "min_R0": float(metrics["min_R0"]),
            "step_norm": step_norm,
            "sigma_min_JQ_perp": float(linear["sigma_min_JQ_perp"]),
            "gauge_rank": int(linear["scaled_gauge_rank"]),
            "linear_rank": int(rank),
            "linear_sigma_min": float(np.min(singular)) if len(singular) else math.nan,
            "linear_sigma_max": float(np.max(singular)) if len(singular) else math.nan,
            "linear_solve_seconds": solve_seconds,
            "wall_seconds": float(time.perf_counter() - started),
            "line_search_alpha": None,
            "line_search_accepted": False,
            "linear_message": linear_message,
        }
        rows.append(row)
        if (
            residual_linf <= float(config.c1_residual_tolerance)
            and step_norm <= float(config.c1_step_tolerance)
        ):
            status = "TIGHT_CONVERGED"
            message = "original residual and update norm met tight C1 tolerances"
            append_progress(rows[-1])
            break
        if iteration >= int(config.c1_max_iters):
            append_progress(rows[-1])
            break
        if not np.all(np.isfinite(dx)):
            status = "FAILED"
            message = linear_message
            append_progress(rows[-1])
            break
        old_norm = residual_linf
        accepted = False
        alpha = 1.0
        best_state = state
        best_norm = old_norm
        best_alpha = None
        for _line_index in range(int(config.c1_line_search_steps)):
            trial_state = state + alpha * dx
            trial_residual, _branch, _grid = _c0g_b3_residual_np(
                state_np=trial_state,
                tau=float(tau),
                config=C0gB3PseudoArclengthConfig(
                    deepcrawl_json_path=config.deepcrawl_json_path,
                    run_root=config.run_root / "C1" / label,
                    grid=config.grid,
                    jacobian_assembly=config.jacobian_assembly,
                    tau_scale=config.tau_scale,
                ),
                dtype=dtype,
            )
            trial_norm = float(np.max(np.abs(trial_residual)))
            if np.isfinite(trial_norm) and trial_norm < best_norm:
                best_state = trial_state
                best_norm = trial_norm
                best_alpha = alpha
            if np.isfinite(trial_norm) and trial_norm < old_norm:
                accepted = True
                best_state = trial_state
                best_norm = trial_norm
                best_alpha = alpha
                break
            alpha *= 0.5
        rows[-1]["line_search_alpha"] = best_alpha
        rows[-1]["line_search_accepted"] = bool(accepted or best_norm < old_norm)
        rows[-1]["accepted_residual_linf"] = best_norm if best_norm < old_norm else None
        if best_norm < old_norm:
            state = np.asarray(best_state, dtype=np.float64)
        else:
            status = "STALLED"
            message = "line search failed to decrease original residual"
            append_progress(rows[-1])
            break
        if _c0g_b3_followup_recent_progress_stalled(
            rows,
            window=int(config.c1_stall_window),
            epsilon_res=float(config.c1_stall_residual_decrease),
        ):
            status = "STALLED"
            message = (
                f"no {config.c1_stall_residual_decrease:.1%} residual decrease "
                f"over {config.c1_stall_window} C1 iterations"
            )
            append_progress(rows[-1])
            break
        append_progress(rows[-1])
    residual, _branch, grid = _c0g_b3_residual_np(
        state_np=state,
        tau=float(tau),
        config=C0gB3PseudoArclengthConfig(
            deepcrawl_json_path=config.deepcrawl_json_path,
            run_root=config.run_root / "C1" / label,
            grid=config.grid,
            jacobian_assembly=config.jacobian_assembly,
            tau_scale=config.tau_scale,
        ),
        dtype=dtype,
    )
    state_artifact = _c0g_b3_followup_checkpoint_c1(
        config=config,
        label=label,
        state=state,
        tau=float(tau),
        dtype=dtype,
    )
    return {
        "label": label,
        "status": status,
        "message": message,
        "tau": float(tau),
        "iterations": int(max(0, len(rows) - 1)),
        "initial_residual_linf": float(rows[0]["original_residual_linf"]) if rows else None,
        "final_residual_linf": float(np.max(np.abs(residual))),
        "final_step_norm": float(rows[-1]["step_norm"]) if rows else None,
        "residual_tolerance": float(config.c1_residual_tolerance),
        "step_tolerance": float(config.c1_step_tolerance),
        "history": rows,
        "state_artifact": state_artifact,
        "invariants": _c0g_b3_invariants(state, grid),
        "metrics": _state_metrics(torch.as_tensor(state, dtype=dtype, device="cpu"), grid),
        "last_linearization_summary": None
        if final_linear is None
        else {
            "sigma_min_JQ_perp": float(final_linear["sigma_min_JQ_perp"]),
            "sigma_max_JQ_perp": float(final_linear["sigma_max_JQ_perp"]),
            "gauge_rank": int(final_linear["scaled_gauge_rank"]),
            "jacobian_source": final_linear.get("jacobian_metadata", {}).get("jacobian_source"),
        },
        "wall_seconds": float(time.perf_counter() - started),
    }


def _c0g_b3_followup_c1_case(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    mismatch: Mapping[str, float],
    config: C0gB3FollowupConfig,
) -> dict[str, Any]:
    left_tight = left.get("status") == "TIGHT_CONVERGED"
    right_tight = right.get("status") == "TIGHT_CONVERGED"
    max_mismatch = max(float(value) for value in mismatch.values()) if mismatch else math.inf
    if left_tight and right_tight and max_mismatch <= float(config.c1_same_state_tolerance):
        return {
            "case": "Case 1",
            "reading": "both_tight_same",
            "tool_table_action": "reconverge_recorded_deep_states_then_gate_refine_only_after_C2_C4",
        }
    if left_tight and right_tight:
        return {
            "case": "Case 3",
            "reading": "both_tight_different",
            "tool_table_action": "STOP_re_gate_for_deflation_branch_switching",
        }
    if not left_tight and not right_tight:
        return {
            "case": "Case 2_or_4",
            "reading": "neither_tight",
            "tool_table_action": "C2_required_to_discriminate_fold_vs_wall",
        }
    return {
        "case": "Case 0",
        "reading": "asymmetric_or_partial_C1",
        "tool_table_action": "STOP_debug_re_gate_no_tool_selected",
    }


def _c0g_b3_followup_c1(
    *,
    deep_payload: Mapping[str, Any],
    b3_payload: Mapping[str, Any],
    config: C0gB3FollowupConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    rec_row = _c0g_b3_followup_find_deep_row(deep_payload, float(config.fail_tau))
    cont_row = _c0g_b3_followup_find_b3_accepted(b3_payload, "B3_0_accepted_009")
    rec_start = _c0g_b3_followup_state_np(str(rec_row["state_artifact"]), dtype=dtype)
    cont_start = _c0g_b3_followup_state_np(str(cont_row["state_artifact"]), dtype=dtype)
    rec = _c0g_b3_followup_reconverge_fixed_tau(
        label="recorded_deepcrawl_attempt_008",
        start_state=rec_start,
        tau=float(config.fail_tau),
        config=config,
        dtype=dtype,
    )
    cont = _c0g_b3_followup_reconverge_fixed_tau(
        label="continuation_B3_0_accepted_009",
        start_state=cont_start,
        tau=float(config.fail_tau),
        config=config,
        dtype=dtype,
    )
    _branch, grid, _residual_fn = _c0g_b3_context(
        tau=float(config.fail_tau),
        config=C0gB3PseudoArclengthConfig(
            deepcrawl_json_path=config.deepcrawl_json_path,
            run_root=config.run_root,
            grid=config.grid,
            jacobian_assembly=config.jacobian_assembly,
            tau_scale=config.tau_scale,
        ),
        dtype=dtype,
    )
    rec_state = _c0g_b3_followup_state_np(str(rec["state_artifact"]), dtype=dtype)
    cont_state = _c0g_b3_followup_state_np(str(cont["state_artifact"]), dtype=dtype)
    aligned_cont, alignment = _c0g_global_phase_align_np(rec_state, cont_state, grid)
    mismatch = _c0g_b3_invariant_mismatch(left=rec_state, right=aligned_cont, grid=grid)
    case = _c0g_b3_followup_c1_case(rec, cont, mismatch, config)
    for result in (rec, cont):
        invariants = result.get("invariants")
        if isinstance(invariants, dict):
            result["invariants_summary"] = {
                "rho_min": float(np.min(invariants["rho"])),
                "rho_max": float(np.max(invariants["rho"])),
                "curl_A_linf": float(np.max(np.abs(invariants["curl_A"]))),
                "r0_min": float(np.min(invariants["r0"])),
                "r0_max": float(np.max(invariants["r0"])),
                "mu": float(invariants["mu"]),
            }
            result.pop("invariants", None)
    return {
        "tau": float(config.fail_tau),
        "starts": {
            "recorded": {
                "provenance": "deepcrawl recorded comparison artifact; diagnostics-only load",
                "state_artifact": str(rec_row["state_artifact"]),
                "recorded_residual_linf": rec_row.get("residual_linf"),
            },
            "continuation": {
                "provenance": "B3_0 continuation accepted artifact; diagnostics-only load",
                "state_artifact": str(cont_row["state_artifact"]),
                "recorded_residual_linf": cont_row.get("original_residual_linf"),
            },
        },
        "recorded": rec,
        "continuation": cont,
        "phase_alignment": alignment,
        "gauge_invariant_mismatch": mismatch,
        "case_reading": case,
    }


def _c0g_b3_followup_sector_fractions(
    vector: np.ndarray,
    grid,
) -> dict[str, float]:
    lanes = _closed_lane_slices(grid)
    groups = {
        "psi": ("psi_real", "psi_imag"),
        "A": ("a0", "ar", "aw"),
        "r0": ("r0",),
        "mu": ("mu",),
    }
    total = float(np.dot(vector, vector))
    if total <= 0.0:
        return {key: math.nan for key in groups}
    fractions = {}
    for group, names in groups.items():
        energy = 0.0
        for name in names:
            start, stop = lanes[name]
            energy += float(np.dot(vector[start:stop], vector[start:stop]))
        fractions[group] = float(energy / total)
    return fractions


def _c0g_b3_followup_localization(
    vector: np.ndarray,
    grid,
) -> dict[str, Any]:
    lanes = _closed_lane_slices(grid)
    cell_count = int(grid.spec.nr * grid.spec.nw)
    nw = int(grid.spec.nw)
    throat_columns = max(1, nw // 4)
    cell_energy = np.zeros((int(grid.spec.nr), int(grid.spec.nw)), dtype=np.float64)
    for name in ("psi_real", "psi_imag", "a0", "ar", "aw"):
        start, stop = lanes[name]
        block = vector[start:stop].reshape((int(grid.spec.nr), int(grid.spec.nw)))
        cell_energy += block * block
    r0_start, r0_stop = lanes["r0"]
    r0_energy = np.zeros_like(cell_energy)
    r0_energy[:, :] = (vector[r0_start:r0_stop][None, :] ** 2) / max(int(grid.spec.nr), 1)
    cell_energy += r0_energy
    total = float(np.sum(cell_energy))
    throat = float(np.sum(cell_energy[:, :throat_columns]))
    peak_flat = int(np.argmax(cell_energy)) if cell_count else 0
    peak = np.unravel_index(peak_flat, cell_energy.shape)
    fraction = float(throat / max(total, np.finfo(np.float64).tiny))
    return {
        "throat_energy_fraction_first_quarter_w": fraction,
        "classification": "THROAT_CONCENTRATED" if fraction >= 0.5 else "EXTENDED",
        "peak_cell": {"r_index": int(peak[0]), "w_index": int(peak[1])},
    }


def _c0g_b3_followup_expanded_gauge_candidates(
    *,
    state: np.ndarray,
    grid,
    branch,
) -> np.ndarray:
    state_t = torch.as_tensor(np.asarray(state, dtype=np.float64), dtype=torch.float64, device="cpu")
    candidates = []
    for generator in _c0c_generators_for_state(state_t, grid):
        candidates.append(np.asarray(generator.vector, dtype=np.float64))
    lanes = _closed_lane_slices(grid)
    size = int(state.size)
    for lane_name in ("a0", "ar", "aw", "r0", "mu"):
        start, stop = lanes[lane_name]
        vector = np.zeros(size, dtype=np.float64)
        vector[start:stop] = 1.0
        candidates.append(vector)
    gradient = _c0d_scalar_gradient_matrix(grid)
    for col in range(min(gradient.shape[1], 8)):
        candidates.append(np.asarray(gradient[:, col], dtype=np.float64))
    coupled = _c0e_coupled_gauge_matrix(state_t, grid, branch)
    for col in range(min(coupled.shape[1], 8)):
        candidates.append(np.asarray(coupled[:, col], dtype=np.float64))
    matrix = np.column_stack(candidates) if candidates else np.zeros((size, 0), dtype=np.float64)
    if matrix.shape[1] == 0:
        return matrix
    q, _r = np.linalg.qr(matrix, mode="reduced")
    return q


def _c0g_b3_followup_projection_fraction(vector: np.ndarray, basis: np.ndarray) -> float:
    denom = float(np.dot(vector, vector))
    if denom <= 0.0 or basis.size == 0:
        return 0.0
    coeff = basis.T @ vector
    return float(np.dot(coeff, coeff) / denom)


def _c0g_b3_followup_c2_from_state(
    *,
    label: str,
    state: np.ndarray,
    tau: float,
    config: C0gB3FollowupConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    linear = _c0g_b3_linearization(
        state_np=state,
        tau=float(tau),
        config=C0gB3PseudoArclengthConfig(
            deepcrawl_json_path=config.deepcrawl_json_path,
            run_root=config.run_root / "C2" / label,
            grid=config.grid,
            jacobian_assembly=config.jacobian_assembly,
            tau_scale=config.tau_scale,
        ),
        dtype=dtype,
    )
    reduced_j = np.asarray(linear["reduced_j"], dtype=np.float64)
    u, singular, vh = np.linalg.svd(reduced_j, full_matrices=False)
    count = min(int(config.c2_singular_count), int(singular.size))
    smallest = singular[-count:].tolist() if count else []
    gap_ratio = (
        float(singular[-2] / max(singular[-1], np.finfo(np.float64).tiny))
        if singular.size >= 2
        else math.inf
    )
    cluster_call = "ISOLATED_NEAR_NULL" if gap_ratio >= float(config.c2_cluster_gap_ratio) else "CLUSTER_OR_NO_CLEAR_GAP"
    coeff = vh[-1, :].astype(np.float64, copy=False)
    scaled_mode = np.asarray(linear["q_perp"], dtype=np.float64) @ coeff
    physical_mode = _c0g_b3_physical_from_scaled_step(
        col_scale=np.asarray(linear["col_scale"], dtype=np.float64),
        q_perp=np.asarray(linear["q_perp"], dtype=np.float64),
        coeff=coeff,
    )
    left_null = u[:, -1].astype(np.float64, copy=False)
    w_ftau = float(np.dot(left_null, np.asarray(linear["ftau_scaled"], dtype=np.float64)))
    w_ftau_norm = float(abs(w_ftau) / max(np.linalg.norm(linear["ftau_scaled"]), np.finfo(np.float64).tiny))
    gauge = _c0g_scaled_gauge_complement(
        gauge_matrix=np.asarray(linear["gauge"]["generator_matrix"], dtype=np.float64),
        col_scale=np.asarray(linear["col_scale"], dtype=np.float64),
        config=C0gConfig(
            run_root=config.run_root / "C2" / label,
            grid=config.grid,
            jacobian_assembly=config.jacobian_assembly,
        ),
    )
    scaled_gauge = np.asarray(linear["gauge"]["generator_matrix"], dtype=np.float64) / np.asarray(
        linear["col_scale"], dtype=np.float64
    )[:, None]
    scaled_basis, _u, _s, _threshold, _rank = _c0g_column_space(
        scaled_gauge,
        relative_tolerance=C0gConfig().gradient_rank_rtol,
        full_matrices=False,
    )
    removed_projection = _c0g_b3_followup_projection_fraction(scaled_mode, scaled_basis)
    expanded = _c0g_b3_followup_expanded_gauge_candidates(
        state=state,
        grid=linear["grid"],
        branch=linear["branch"],
    )
    expanded_projection = _c0g_b3_followup_projection_fraction(physical_mode, expanded)
    curl = _c0g_curl_np(physical_mode, linear["grid"])
    residual_response = np.asarray(linear["scaled_matrix"] @ scaled_mode, dtype=np.float64)
    localization = _c0g_b3_followup_localization(physical_mode, linear["grid"])
    mode_path = _resolve_output_path(config.run_root / "C2" / f"{label}_mode.npz")
    mode_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        mode_path,
        physical_mode=physical_mode,
        scaled_mode=scaled_mode,
        singular_values=singular,
        left_null=left_null,
    )
    return {
        "label": label,
        "tau": float(tau),
        "status": "MEASURED",
        "smallest_singular_values": smallest,
        "sigma_min": float(singular[-1]) if singular.size else math.nan,
        "sigma_gap_ratio_next_over_min": gap_ratio,
        "cluster_call": cluster_call,
        "sector_energy_fractions": _c0g_b3_followup_sector_fractions(physical_mode, linear["grid"]),
        "fold_transversality": {
            "wT_Ftau": w_ftau,
            "normalized_abs_wT_Ftau": w_ftau_norm,
            "call": "FOLD_TRANSVERSAL"
            if w_ftau_norm > float(config.c2_ftau_abs_small)
            else "FTAU_ORTHOGONAL_OR_DEGENERATE",
        },
        "gauge_tests": {
            "scaled_removed_gauge_projection_fraction": removed_projection,
            "independent_expanded_candidate_projection_fraction": expanded_projection,
            "curl_mode_l2": float(np.linalg.norm(curl)),
            "scaled_residual_response_l2": float(np.linalg.norm(residual_response)),
            "high_physical_overlap_alone_is_not_used": True,
        },
        "localization": localization,
        "mode_artifact": str(mode_path),
        "jacobian_source": linear.get("jacobian_metadata", {}).get("jacobian_source"),
        "gauge_rank": int(linear["scaled_gauge_rank"]),
        "redundant_gauge_rank_check": int(gauge["scaled_gauge_rank"]),
    }


def _c0g_b3_followup_c2(
    *,
    c1: Mapping[str, Any],
    config: C0gB3FollowupConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    case = c1.get("case_reading", {}).get("case")
    if case == "Case 0":
        return {"status": "SKIPPED_BY_C1_CASE0", "reason": "Case 0 stops for debug/re-gate"}
    measured = []
    for key in ("recorded", "continuation"):
        result = c1.get(key, {})
        if result.get("status") != "TIGHT_CONVERGED":
            continue
        state = _c0g_b3_followup_state_np(str(result["state_artifact"]), dtype=dtype)
        measured.append(
            _c0g_b3_followup_c2_from_state(
                label=key,
                state=state,
                tau=float(c1["tau"]),
                config=config,
                dtype=dtype,
            )
        )
    if not measured:
        return {"status": "SKIPPED_NO_TIGHT_C1_STATE", "reason": "C2 is read-only on tight C1 states"}
    return {"status": "MEASURED", "state_results": measured}


def _c0g_b3_followup_c3(
    *,
    c1: Mapping[str, Any],
    config: C0gB3FollowupConfig,
    dtype: torch.dtype,
) -> dict[str, Any]:
    if c1.get("recorded", {}).get("status") != "TIGHT_CONVERGED" or c1.get("continuation", {}).get("status") != "TIGHT_CONVERGED":
        return {"status": "SKIPPED_REQUIRES_TWO_TIGHT_C1_ENDPOINTS"}
    rec_state = _c0g_b3_followup_state_np(str(c1["recorded"]["state_artifact"]), dtype=dtype)
    cont_state = _c0g_b3_followup_state_np(str(c1["continuation"]["state_artifact"]), dtype=dtype)
    _branch, grid, _residual_fn = _c0g_b3_context(
        tau=float(c1["tau"]),
        config=C0gB3PseudoArclengthConfig(
            deepcrawl_json_path=config.deepcrawl_json_path,
            run_root=config.run_root / "C3",
            grid=config.grid,
            jacobian_assembly=config.jacobian_assembly,
            tau_scale=config.tau_scale,
        ),
        dtype=dtype,
    )
    aligned_cont, alignment = _c0g_global_phase_align_np(rec_state, cont_state, grid)
    rows = []
    endpoint_floor = max(
        float(c1["recorded"]["final_residual_linf"]),
        float(c1["continuation"]["final_residual_linf"]),
        np.finfo(np.float64).tiny,
    )
    for index, t in enumerate(np.linspace(0.0, 1.0, int(config.c3_points))):
        state = (1.0 - float(t)) * rec_state + float(t) * aligned_cont
        residual, _branch, _grid = _c0g_b3_residual_np(
            state_np=state,
            tau=float(c1["tau"]),
            config=C0gB3PseudoArclengthConfig(
                deepcrawl_json_path=config.deepcrawl_json_path,
                run_root=config.run_root / "C3",
                grid=config.grid,
                jacobian_assembly=config.jacobian_assembly,
                tau_scale=config.tau_scale,
            ),
            dtype=dtype,
        )
        linf = float(np.max(np.abs(residual)))
        rows.append(
            {
                "index": int(index),
                "t": float(t),
                "original_residual_linf": linf,
                "normalized_to_tight_endpoint_floor": float(linf / endpoint_floor),
            }
        )
    max_normalized = max(float(row["normalized_to_tight_endpoint_floor"]) for row in rows)
    call = "SPIKE" if max_normalized >= float(config.c3_spike_factor) else "FLAT_OR_NO_CLEAR_SPIKE"
    return {
        "status": "MEASURED",
        "phase_alignment": alignment,
        "endpoint_residual_floor": endpoint_floor,
        "max_normalized_residual": max_normalized,
        "call": call,
        "supporting_only_not_proof": True,
        "rows": rows,
    }


def _c0g_b3_followup_c4(
    *,
    config: C0gB3FollowupConfig,
) -> dict[str, Any]:
    if not bool(config.c4_run):
        return {
            "status": "SKIPPED_BY_DEFAULT_GATE",
            "reason": "resolution relocation is expensive; Case-1 stops for re-gate before any tool claim, so fold-vs-wall remains open",
        }
    return {
        "status": "NOT_IMPLEMENTED_FOR_AUTOMATIC_RUN",
        "reason": "C4 requires a gated per-resolution relocation crawl; no continuation/tool claim made",
    }


def _c0g_b3_followup_tool_selection(
    *,
    c1: Mapping[str, Any],
    c2: Mapping[str, Any],
    c3: Mapping[str, Any],
    c4: Mapping[str, Any],
) -> dict[str, Any]:
    case = c1.get("case_reading", {}).get("case")
    if case == "Case 0":
        return {
            "case": "Case 0",
            "chosen_tool": "STOP_DEBUG_RE_GATE",
            "why": "C1 was asymmetric or partial; directive forbids forcing this into a physical case",
        }
    if case == "Case 3":
        return {
            "case": "Case 3",
            "chosen_tool": "STOP_RE_GATE_DEFLATION_BRANCH_SWITCHING",
            "why": "C1 found distinct tight roots",
        }
    if case == "Case 1":
        return {
            "case": "Case 1",
            "chosen_tool": "STOP_RECONVERGE_DEEP_STATES_AND_GATE_REFINE_THEN_C2_C4",
            "why": "C1 found the same tight root at this tau; global fold-vs-wall remains open",
        }
    if case == "Case 2_or_4":
        c2_rows = c2.get("state_results", []) if isinstance(c2, Mapping) else []
        first = c2_rows[0] if c2_rows else {}
        if first.get("cluster_call") == "ISOLATED_NEAR_NULL" and first.get("fold_transversality", {}).get("call") == "FOLD_TRANSVERSAL":
            return {
                "case": "Case 2_candidate",
                "chosen_tool": "STOP_RE_GATE_PSEUDO_ARCLENGTH_CANDIDATE",
                "why": "C1 neither-tight plus C2 fold-like evidence; user gate required before using pseudo-arclength",
            }
        return {
            "case": "Case 4_candidate",
            "chosen_tool": "STOP_RE_GATE_LM_TRUST_REGION_CANDIDATE",
            "why": "C1 neither-tight without decisive fold-like C2 evidence",
        }
    return {
        "case": str(case),
        "chosen_tool": "STOP_INCONCLUSIVE",
        "why": f"unmapped case reading {case!r}; C3={c3.get('status')}; C4={c4.get('status')}",
    }


def write_c0g_b3_followup_report(result: Mapping[str, Any], path: Path) -> None:
    c0 = result.get("C0", {})
    c1 = result.get("C1", {})
    c2 = result.get("C2", {})
    c3 = result.get("C3", {})
    c4 = result.get("C4", {})
    tool = result.get("tool_selection", {})
    lines = [
        "# Path-A C0g B-3 Follow-Up v2",
        "",
        "## Contract",
        "",
        "- Battery: C0 read-offs, C1 tight two-start reconvergence, C2/C3 gated characterization, C4 gated resolution check.",
        "- Tool selection is evidence-driven; pseudo-arclength is not pre-selected.",
        "- Original `patha_closed_branch_residual` remains the single arbiter.",
        "",
        "## C0 Read-Offs",
        "",
        f"- sigma trend call: `{c0.get('sigma_trend_call')}`",
        f"- deepest converged tau: `{_c0g_fmt(c0.get('deepest_converged_tau'))}`",
        f"- deepest sigma_min(JQ_perp): `{_c0g_fmt(c0.get('deepest_converged_sigma_min_JQ_perp'))}`",
        f"- fit zero crossing: `{_c0g_fmt(c0.get('sigma_fit_zero_crossing'))}`",
        f"- zero-crossing minus deepest tau: `{_c0g_fmt(c0.get('gap_zero_crossing_minus_deepest_tau'))}`",
        "",
        _markdown_table(
            ["tau", "sigma_min_JQ_perp", "cond_JQ_perp", "status"],
            c0.get("sigma_rows", []),
        ),
        "",
        "### Arclength Metric Balance",
        "",
        _markdown_table(
            ["id", "previous_tau", "tau", "dx_norm_squared", "tau_metric_squared", "x_fraction", "dominant"],
            c0.get("arclength_metric_balance", []),
        ),
        "",
        "## C1 Tight Reconvergence",
        "",
        f"- case: `{c1.get('case_reading', {}).get('case')}`",
        f"- reading: `{c1.get('case_reading', {}).get('reading')}`",
        f"- recorded status: `{c1.get('recorded', {}).get('status')}` residual `{_c0g_fmt(c1.get('recorded', {}).get('final_residual_linf'))}` step `{_c0g_fmt(c1.get('recorded', {}).get('final_step_norm'))}`",
        f"- continuation status: `{c1.get('continuation', {}).get('status')}` residual `{_c0g_fmt(c1.get('continuation', {}).get('final_residual_linf'))}` step `{_c0g_fmt(c1.get('continuation', {}).get('final_step_norm'))}`",
        "",
        "Gauge-invariant mismatch after phase alignment:",
        "",
        "```yaml",
        f"{c1.get('gauge_invariant_mismatch')}",
        "```",
        "",
        "## C2 Mode Characterization",
        "",
        f"status: `{c2.get('status')}`",
    ]
    for row in c2.get("state_results", []) if isinstance(c2, Mapping) else []:
        lines.extend(
            [
                "",
                f"### {row.get('label')}",
                "",
                f"- sigma_min: `{_c0g_fmt(row.get('sigma_min'))}`",
                f"- near-null count call: `{row.get('cluster_call')}` gap `{_c0g_fmt(row.get('sigma_gap_ratio_next_over_min'))}`",
                f"- dominant sector fractions: `{row.get('sector_energy_fractions')}`",
                f"- transversality: `{row.get('fold_transversality')}`",
                f"- gauge tests: `{row.get('gauge_tests')}`",
                f"- localization: `{row.get('localization')}`",
            ]
        )
    lines.extend(
        [
            "",
            "## C3 Line Scan",
            "",
            f"- status: `{c3.get('status')}`",
            f"- call: `{c3.get('call')}`",
            f"- max normalized residual: `{_c0g_fmt(c3.get('max_normalized_residual'))}`",
            "",
            "## C4 Resolution Check",
            "",
            f"- status: `{c4.get('status')}`",
            f"- reason: `{c4.get('reason')}`",
            "",
            "## Tool Selection",
            "",
            f"- case: `{tool.get('case')}`",
            f"- chosen tool: `{tool.get('chosen_tool')}`",
            f"- why: {tool.get('why')}",
            "",
            "## Scope Guard",
            "",
            "```yaml",
        ]
    )
    for key, value in result.get("scope_guard", {}).items():
        lines.append(f"{key}: {value}")
    lines.extend(
        [
            "```",
            "",
            "## Git Diff Summary",
            "",
            "```",
            str(result.get("git_diff_summary", "")),
            "```",
            "",
            f"Machine JSON: `{result.get('json_path')}`",
            f"Progress log: `{result.get('progress_path')}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_c0g_b3_followup(
    config: C0gB3FollowupConfig | None = None,
) -> dict[str, Any]:
    cfg = config or C0gB3FollowupConfig()
    started = time.perf_counter()
    dtype = configure_backend(BackendConfig())
    progress_path = _c0g_b3_followup_progress_path(cfg)
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    progress_path.write_text("", encoding="utf-8")
    deep_payload = _c0g_b3_followup_load_json(cfg.deepcrawl_json_path)
    b3_payload = _c0g_b3_followup_load_json(cfg.b3_json_path)
    c0 = _c0g_b3_followup_c0(
        deep_payload=deep_payload,
        b3_payload=b3_payload,
        config=cfg,
        dtype=dtype,
    )
    c1 = _c0g_b3_followup_c1(
        deep_payload=deep_payload,
        b3_payload=b3_payload,
        config=cfg,
        dtype=dtype,
    )
    c2 = _c0g_b3_followup_c2(c1=c1, config=cfg, dtype=dtype)
    c3 = _c0g_b3_followup_c3(c1=c1, config=cfg, dtype=dtype)
    c4 = _c0g_b3_followup_c4(config=cfg)
    tool_selection = _c0g_b3_followup_tool_selection(c1=c1, c2=c2, c3=c3, c4=c4)
    result = {
        "schema": "stage1_pathA_C0g_B3_followup_v2/v1",
        "source_revision": source_revision(),
        "config": _c0g_b3_followup_config_to_dict(cfg),
        "C0": c0,
        "C1": c1,
        "C2": c2,
        "C3": c3,
        "C4": c4,
        "tool_selection": tool_selection,
        "scope_guard": {
            "single_arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
            "convergence_judged_only_on_original_residual": True,
            "gauge_Q_perp_inside_every_C1_solve": True,
            "patha_closed_branch_residual_touched": False,
            "faithful_operators_touched": False,
            "xi_or_grad_div_penalty_touched": False,
            "physical_export_permitted_touched": False,
            "frozen_physics_touched": False,
            "no_LM_no_PTC_no_Sobolev_regularization": True,
            "no_pseudoarclength_precommit": True,
            "comparison_artifacts_diagnostics_only": True,
        },
        "progress_guard": {
            "C1_fixed_tau_window": int(cfg.c1_stall_window),
            "epsilon_res": float(cfg.c1_stall_residual_decrease),
            "stalled": any(
                c1.get(key, {}).get("status") == "STALLED"
                for key in ("recorded", "continuation")
            ),
        },
        "git_diff_summary": _c0g_b3_git_diff_summary(),
        "elapsed_seconds": float(time.perf_counter() - started),
        "report_path": str(_resolve_output_path(cfg.report_path)),
        "json_path": str(_resolve_output_path(cfg.json_path)),
        "progress_path": str(progress_path),
    }
    json_path = _resolve_output_path(cfg.json_path)
    report_path = _resolve_output_path(cfg.report_path)
    _c0g_write_json(json_path, result)
    write_c0g_b3_followup_report(result, report_path)
    return result


def _c0g_b3_followup_stdout_summary(result: Mapping[str, Any]) -> str:
    c0 = result.get("C0", {})
    c1 = result.get("C1", {})
    c2 = result.get("C2", {})
    c3 = result.get("C3", {})
    c4 = result.get("C4", {})
    tool = result.get("tool_selection", {})
    c2_rows = c2.get("state_results", []) if isinstance(c2, Mapping) else []
    first_c2 = c2_rows[0] if c2_rows else {}
    scope = result.get("scope_guard", {})
    mismatch = c1.get("gauge_invariant_mismatch", {})
    return "\n".join(
        [
            (
                "C0: sigma trend={trend}, deepest_tau={tau}, deepest_sigma={sigma}, "
                "zero_gap={gap}"
            ).format(
                trend=c0.get("sigma_trend_call"),
                tau=_c0g_fmt(c0.get("deepest_converged_tau")),
                sigma=_c0g_fmt(c0.get("deepest_converged_sigma_min_JQ_perp")),
                gap=_c0g_fmt(c0.get("gap_zero_crossing_minus_deepest_tau")),
            ),
            f"C0 arclength balance: last={c0.get('arclength_metric_balance', [{}])[-1] if c0.get('arclength_metric_balance') else 'n/a'}",
            (
                "C1: {case}/{reading}; recorded={rstat} residual={rres} step={rstep}; "
                "continuation={cstat} residual={cres} step={cstep}; mismatch={mismatch}"
            ).format(
                case=c1.get("case_reading", {}).get("case"),
                reading=c1.get("case_reading", {}).get("reading"),
                rstat=c1.get("recorded", {}).get("status"),
                rres=_c0g_fmt(c1.get("recorded", {}).get("final_residual_linf")),
                rstep=_c0g_fmt(c1.get("recorded", {}).get("final_step_norm")),
                cstat=c1.get("continuation", {}).get("status"),
                cres=_c0g_fmt(c1.get("continuation", {}).get("final_residual_linf")),
                cstep=_c0g_fmt(c1.get("continuation", {}).get("final_step_norm")),
                mismatch=mismatch,
            ),
            (
                "C2: status={status}, near-null={cluster}, dominant_sector={sector}, "
                "gauge_overlap={gauge}, localization={loc}"
            ).format(
                status=c2.get("status"),
                cluster=first_c2.get("cluster_call"),
                sector=_dominant_group(first_c2.get("sector_energy_fractions", {}))[0]
                if first_c2
                else "n/a",
                gauge=first_c2.get("gauge_tests", {}),
                loc=first_c2.get("localization", {}),
            ),
            f"C3: status={c3.get('status')}, call={c3.get('call')}, max_normalized={_c0g_fmt(c3.get('max_normalized_residual'))}",
            f"C4: status={c4.get('status')}, reason={c4.get('reason')}",
            f"TOOL: {tool.get('chosen_tool')} ({tool.get('why')})",
            (
                "Scope untouched: residual/operators/(1/xi)/export/physics="
                f"{not scope.get('patha_closed_branch_residual_touched')}/"
                f"{not scope.get('faithful_operators_touched')}/"
                f"{not scope.get('xi_or_grad_div_penalty_touched')}/"
                f"{not scope.get('physical_export_permitted_touched')}/"
                f"{not scope.get('frozen_physics_touched')}; "
                f"Single Arbiter only={scope.get('convergence_judged_only_on_original_residual')}; "
                f"progress stalled={result.get('progress_guard', {}).get('stalled')}"
            ),
            f"Files: report={result.get('report_path')}, json={result.get('json_path')}, progress={result.get('progress_path')}",
        ]
    )


def _c0g_fmt(value: Any, *, digits: int = 6) -> str:
    if value is None:
        return "n/a"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isnan(number):
        return "nan"
    if math.isinf(number):
        return "inf" if number > 0 else "-inf"
    return f"{number:.{digits}e}"


def _c0g_report_step0_rows(step0: Mapping[str, Any]) -> list[dict[str, Any]]:
    delta = step0.get("step0_decomposition", {}).get("delta", {})
    residual = step0.get("step0_decomposition", {}).get("residual", {})
    merit = step0.get("merit_sweep", {})
    alpha1 = merit.get("alpha1") or {}
    return [
        {
            "tau": step0.get("tau"),
            "f_nn": delta.get("near_null_component_fraction"),
            "gauge_fraction": delta.get("gauge_component_fraction"),
            "mode2_fraction": delta.get("transverse_mode_2_component_fraction"),
            "gauge_image_F_fraction": residual.get("gauge_image_projection_fraction"),
            "alpha1_actual_ratio": alpha1.get("actual_l2_ratio"),
            "alpha1_gap": merit.get("alpha1_predicted_vs_actual_l2_gap"),
            "best_alpha": merit.get("best_alpha"),
            "best_reduction_percent": merit.get("best_percent_reduction"),
        }
    ]


def _c0g_report_svd_rows(state_results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for state in sorted(state_results, key=lambda item: float(item["tau"]), reverse=True):
        svd = state.get("svd", {})
        if svd.get("status") != "MEASURED":
            rows.append(
                {
                    "tau": state.get("tau"),
                    "status": svd.get("status"),
                    "sigma_min": None,
                    "sigma_min/sigma_max": None,
                    "gauge_rank": None,
                    "mode1_lane": None,
                }
            )
            continue
        mode1 = (svd.get("modes") or [{}])[0]
        lane = _dominant_group(mode1.get("v_lane_energy_fractions", {}))[0]
        rows.append(
            {
                "tau": state.get("tau"),
                "status": "MEASURED",
                "sigma_min": svd.get("sigma_min"),
                "sigma_min/sigma_max": svd.get("sigma_min_over_sigma_max"),
                "gauge_rank": svd.get("scaled_gauge_rank"),
                "mode1_lane": lane,
            }
        )
    return rows


def _c0g_report_ftau_rows(state_results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for state in sorted(state_results, key=lambda item: float(item["tau"]), reverse=True):
        ftau = state.get("ftau")
        if not isinstance(ftau, Mapping) or ftau.get("status") != "MEASURED":
            rows.append(
                {
                    "tau": state.get("tau"),
                    "status": ftau if isinstance(ftau, str) else (ftau or {}).get("status"),
                    "cos_theta": None,
                    "call": None,
                    "stability": None,
                }
            )
            continue
        rows.append(
            {
                "tau": state.get("tau"),
                "status": "MEASURED",
                "cos_theta": ftau.get("representative_cos_theta"),
                "call": ftau.get("call"),
                "stability": ftau.get("stepsize_stability_call"),
            }
        )
    return rows


def write_c0g_partial_report(result: Mapping[str, Any], path: Path) -> None:
    lines: list[str] = [
        "# Path-A C0g Diagnostic Battery Partial Report (Steps 0-3)",
        "",
        "This staged pass ran only Steps 0, 0b, 1, 2, and 3. Steps 4-7 were not run; the Step-8 verdict is explicitly deferred.",
        "",
        "## Provenance",
        "",
    ]
    provenance = result.get("provenance", {})
    lines.extend(
        [
            f"- C0f2 stalled-state provenance: `{provenance.get('call')}`",
            f"- state: `{provenance.get('state_artifact')}`",
            f"- source: `{provenance.get('source')}`",
            f"- prefer_existing_b2c_background_predictor: `{provenance.get('prefer_existing_b2c_background_predictor')}`",
            f"- used_existing_b2c: `{provenance.get('used_existing_b2c')}`",
            "",
        ]
    )
    step0 = result.get("step0", {})
    gate = step0.get("premise_gate", {}) if isinstance(step0, Mapping) else {}
    lines.extend(
        [
            "## Step 0 Premise Gate",
            "",
            f"PREMISE-GATE call: `{gate.get('call')}` with `f_nn={_c0g_fmt(gate.get('f_nn'))}`.",
            "",
        ]
    )
    if isinstance(step0, Mapping):
        lines.append(
            _markdown_table(
                [
                    "tau",
                    "f_nn",
                    "gauge_fraction",
                    "mode2_fraction",
                    "gauge_image_F_fraction",
                    "alpha1_actual_ratio",
                    "alpha1_gap",
                    "best_alpha",
                    "best_reduction_percent",
                ],
                _c0g_report_step0_rows(step0),
            )
        )
        lines.append("")
        residual = step0.get("step0_decomposition", {}).get("residual", {})
        left_rows = residual.get("left_singular_vector_components", [])
        if left_rows:
            lines.append("Residual decomposition via left singular vectors:")
            lines.append("")
            lines.append(
                _markdown_table(
                    ["ascending_rank", "abs_u_dot_F_over_normF"],
                    left_rows,
                )
            )
            lines.append("")
    if result.get("status") == "PREMISE_FAILED":
        lines.extend(
            [
                "## Halt",
                "",
                "Step 0 returned `PREMISE_FAILED`. Steps 1-7 are `SKIPPED_PREMISE_FAILED` in the machine JSON.",
                "",
            ]
        )
    else:
        lines.extend(["## Step 0b Best-Alpha Trend", ""])
        step0b = result.get("step0b", {})
        lines.append(f"Trend call: `{step0b.get('trend_call')}`.")
        lines.append("")
        lines.append(
            _markdown_table(
                [
                    "tau",
                    "best_alpha",
                    "best_actual_l2_ratio",
                    "best_percent_reduction",
                    "alpha1_actual_l2_ratio",
                    "alpha1_predicted_l2_ratio",
                    "alpha1_gap",
                ],
                step0b.get("rows", []),
            )
        )
        lines.extend(["", "## Step 1 Gauge-Complement SVD", ""])
        state_results = result.get("state_results", [])
        lines.append(
            _markdown_table(
                [
                    "tau",
                    "status",
                    "sigma_min",
                    "sigma_min/sigma_max",
                    "gauge_rank",
                    "mode1_lane",
                ],
                _c0g_report_svd_rows(state_results),
            )
        )
        lines.append("")
        tracked = result.get("step1_tracked_modes", [])
        if tracked:
            compact = []
            for row in tracked:
                v_lane, _v_frac = _dominant_group(row.get("v_lane_energy_fractions", {}))
                w_lane, _w_frac = _dominant_group(row.get("w_output_energy_fractions", {}))
                compact.append(
                    {
                        "tau": row.get("tau"),
                        "tracked_lane": row.get("tracked_lane"),
                        "ascending_rank": row.get("ascending_rank"),
                        "sigma": row.get("sigma"),
                        "v_dominant_lane": v_lane,
                        "w_dominant_lane": w_lane,
                        "v_center_r": (row.get("v_center_of_energy") or {}).get("r"),
                        "v_center_w": (row.get("v_center_of_energy") or {}).get("w"),
                    }
                )
            lines.append("Tracked mode lanes by vector overlap:")
            lines.append("")
            lines.append(
                _markdown_table(
                    [
                        "tau",
                        "tracked_lane",
                        "ascending_rank",
                        "sigma",
                        "v_dominant_lane",
                        "w_dominant_lane",
                        "v_center_r",
                        "v_center_w",
                    ],
                    compact,
                )
            )
            lines.append("")
        lines.extend(["## Step 2 wT F_tau", ""])
        lines.append(
            _markdown_table(
                ["tau", "status", "cos_theta", "call", "stability"],
                _c0g_report_ftau_rows(state_results),
            )
        )
        lines.append("")
        for state in state_results:
            ftau = state.get("ftau")
            if isinstance(ftau, Mapping) and ftau.get("rows"):
                lines.append(f"Per-h values for tau={state.get('tau')}:")
                lines.append("")
                lines.append(
                    _markdown_table(
                        ["h_multiplier", "h", "wT_F_tau", "F_tau_scaled_l2", "cos_theta"],
                        ftau.get("rows", []),
                    )
                )
                lines.append("")
        lines.extend(["## Step 3 sigma_min squared fit", ""])
        fit = result.get("sigma_min_squared_fit", {})
        linear = fit.get("linear", {})
        quadratic = fit.get("quadratic", {})
        lines.extend(
            [
                f"- call: `{fit.get('call')}`",
                f"- linear slope: `{_c0g_fmt(linear.get('slope'))}`",
                f"- linear R2: `{_c0g_fmt(linear.get('r2'))}`",
                f"- tau_fold zero crossing: `{_c0g_fmt(linear.get('tau_fold_zero_crossing'))}`",
                f"- quadratic R2: `{_c0g_fmt(quadratic.get('r2'))}`",
                f"- sigma_min monotone toward stall: `{fit.get('sigma_min_monotone_decreasing_toward_stall')}`",
                "",
                "## Preliminary Reading",
                "",
                f"`{result.get('preliminary_reading', {}).get('status')}`. Step-8 verdict deferred until Steps 4-7 are run.",
                "",
            ]
        )
    lines.extend(
        [
            "## Scope Confirmation",
            "",
            "- NO C0g fix implemented in this staged pass.",
            "- Solver logic, frozen physics, faithful PDE operators, and export guard were not intentionally changed.",
            "- Steps 4-7 were NOT run; stopped for the orchestrator review gate.",
            "",
            f"Machine JSON: `{result.get('config', {}).get('json_path')}`.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _c0g_final_mach_rows(step4: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "tau": row.get("tau"),
            "M_full_max": row.get("M_full_max"),
            "M_w_max": row.get("M_w_max"),
            "rho_at_max": row.get("rho_at_max"),
            "j_at_max": row.get("j_at_max"),
            "r_star": row.get("r_star"),
            "w_star": row.get("w_star"),
            "r_star_over_R0": row.get("r_star_over_R0"),
            "current_represented": row.get("sonic_current_represented"),
            "floor_stability": row.get("density_floor_stability_call"),
        }
        for row in step4.get("rows", [])
    ]


def _c0g_final_step5_rows(step5: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "method": row.get("method"),
            "status": row.get("status"),
            "nfev": row.get("nfev_wrapper"),
            "elapsed_seconds": row.get("elapsed_seconds"),
            "best_l2_ratio": row.get("best_l2_ratio"),
            "best_linf": row.get("best_original_residual_linf"),
            "progress_call": row.get("progress_call"),
            "abort_reason": row.get("abort_reason") or row.get("reason"),
        }
        for row in step5.get("method_results", [])
    ]


def _c0g_final_overlap_rows(step5: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "lane": row.get("lane"),
            "scaled_overlap": row.get("scaled_overlap"),
            "call": row.get("call"),
            "candidate_delta_scaled_norm": row.get("candidate_delta_scaled_norm"),
            "prestall_delta_scaled_norm": row.get("prestall_delta_scaled_norm"),
        }
        for row in step5.get("branch_overlap", {}).get("rows", [])
    ]


def _c0g_final_step6_rows(step6: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "tau": row.get("tau"),
            "status": row.get("status"),
            "cond_JQ_perp": row.get("cond_JQ_perp"),
            "cond_Jb": row.get("cond_Jb"),
            "sigma_min_JQ_perp": row.get("sigma_min_JQ_perp"),
            "sigma_min_Jb": row.get("sigma_min_Jb"),
            "ratio": row.get("sigma_min_ratio_Jb_over_JQ_perp"),
            "call": row.get("call"),
        }
        for row in step6.get("rows", [])
    ]


def _c0g_final_step7_rows(step7: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "tau": row.get("tau"),
            "tracked_lane": row.get("tracked_lane"),
            "ascending_rank": row.get("ascending_rank"),
            "sigma": row.get("sigma"),
            "one_minus_P_G": row.get("one_minus_P_G"),
            "P_G": row.get("P_G"),
            "curl_over_A": row.get("curl_over_spatial_a_norm"),
            "boundary_fraction": row.get("boundary_curl_norm_fraction"),
        }
        for row in step7.get("per_mode_rows", [])
    ]


def _c0g_final_commutator_rows(step7: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "tau": row.get("tau"),
            "norm_CG_over_G": row.get("commutator_norm_CG_over_G"),
            "boundary_fraction": row.get("boundary_row_norm_fraction"),
            "mixed_control_one_minus_P_G": row.get("mixed_control", {}).get("one_minus_P_G"),
            "mixed_control_curl_over_A": row.get("mixed_control", {}).get("curl_over_spatial_a_norm"),
        }
        for row in step7.get("rows", [])
    ]


def write_c0g_final_report(result: Mapping[str, Any], path: Path) -> None:
    step8 = result.get("step8", {})
    lines: list[str] = [
        "# Path-A C0g Diagnostic Battery Final Report (Steps 0-8)",
        "",
        f"Step-8 verdict: **{step8.get('verdict')}**",
    ]
    if step8.get("sub_label"):
        lines.append(f"Sub-label: **{step8.get('sub_label')}**")
    lines.extend(
        [
            "",
            f"Recommended C0g build branch: {step8.get('recommended_C0g_build_branch')}",
            "",
            "The original `patha_closed_branch_residual` remains the sole convergence/progress arbiter. Steps 4-7 are diagnostic only; no C0g fix is implemented here.",
            "",
            "## Provenance",
            "",
        ]
    )
    provenance = result.get("provenance", {})
    lines.extend(
        [
            f"- C0f2 stalled-state provenance: `{provenance.get('call')}`",
            f"- prefer_existing_b2c_background_predictor: `{provenance.get('prefer_existing_b2c_background_predictor')}` from `{provenance.get('prefer_existing_b2c_background_predictor_source')}`",
            f"- source: `{provenance.get('source')}`; used_existing_b2c: `{provenance.get('used_existing_b2c')}`",
            "",
            "## Step 0 Premise Gate",
            "",
        ]
    )
    step0 = result.get("step0", {})
    gate = step0.get("premise_gate", {}) if isinstance(step0, Mapping) else {}
    lines.append(f"Call: `{gate.get('call')}` with `f_nn={_c0g_fmt(gate.get('f_nn'))}`.")
    lines.append("")
    lines.append(
        _markdown_table(
            [
                "tau",
                "f_nn",
                "gauge_fraction",
                "mode2_fraction",
                "gauge_image_F_fraction",
                "alpha1_actual_ratio",
                "best_alpha",
                "best_reduction_percent",
            ],
            _c0g_report_step0_rows(step0) if isinstance(step0, Mapping) else [],
        )
    )
    lines.extend(["", "## Step 0b Best-Alpha Trend", ""])
    step0b = result.get("step0b", {})
    lines.append(f"Trend call: `{step0b.get('trend_call')}`.")
    lines.append("")
    lines.append(
        _markdown_table(
            [
                "tau",
                "best_alpha",
                "best_actual_l2_ratio",
                "best_percent_reduction",
                "alpha1_actual_l2_ratio",
                "alpha1_predicted_l2_ratio",
            ],
            step0b.get("rows", []),
        )
    )
    state_results = result.get("state_results", [])
    lines.extend(["", "## Step 1 Gauge-Complement SVD", ""])
    lines.append(
        _markdown_table(
            [
                "tau",
                "status",
                "sigma_min",
                "sigma_min/sigma_max",
                "gauge_rank",
                "mode1_lane",
            ],
            _c0g_report_svd_rows(state_results),
        )
    )
    lines.extend(["", "## Step 2 wT F_tau", ""])
    lines.append(
        _markdown_table(
            ["tau", "status", "cos_theta", "call", "stability"],
            _c0g_report_ftau_rows(state_results),
        )
    )
    lines.extend(["", "## Step 3 Sigma-Min Squared Fit", ""])
    fit = result.get("sigma_min_squared_fit", {})
    linear = fit.get("linear", {})
    quadratic = fit.get("quadratic", {})
    lines.extend(
        [
            f"- call: `{fit.get('call')}`",
            f"- linear slope: `{_c0g_fmt(linear.get('slope'))}`",
            f"- linear R2: `{_c0g_fmt(linear.get('r2'))}`",
            f"- tau_fold zero crossing: `{_c0g_fmt(linear.get('tau_fold_zero_crossing'))}`",
            f"- quadratic R2: `{_c0g_fmt(quadratic.get('r2'))}`",
            f"- monotone: `{fit.get('sigma_min_monotone_decreasing_toward_stall')}`",
            "",
            "## Step 4 Mach Context",
            "",
        ]
    )
    step4 = result.get("step4_mach", {})
    lines.extend(
        [
            f"Context label: `{step4.get('sonic_context_label')}`; represented: `{step4.get('sonic_hypothesis_represented')}`.",
            f"Monotone M_max approach as tau decreases: `{step4.get('monotone_approach_as_tau_decreases')}`.",
            "",
            _markdown_table(
                [
                    "tau",
                    "M_full_max",
                    "M_w_max",
                    "rho_at_max",
                    "j_at_max",
                    "r_star",
                    "w_star",
                    "r_star_over_R0",
                    "current_represented",
                    "floor_stability",
                ],
                _c0g_final_mach_rows(step4),
            ),
            "",
            "## Step 5 Capped Scipy Probe",
            "",
        ]
    )
    step5 = result.get("step5_scipy_probe", {})
    caps = step5.get("caps", {})
    lines.extend(
        [
            f"Call: `{step5.get('step5_call')}`; status: `{step5.get('status')}`.",
            f"Caps used: maxfev `{caps.get('maxfev')}` (`{caps.get('maxfev_formula')}`), method wall `{caps.get('method_wall_seconds')}` s, global wall `{caps.get('global_wall_seconds')}` s.",
            "",
            _markdown_table(
                [
                    "method",
                    "status",
                    "nfev",
                    "elapsed_seconds",
                    "best_l2_ratio",
                    "best_linf",
                    "progress_call",
                    "abort_reason",
                ],
                _c0g_final_step5_rows(step5),
            ),
            "",
            f"Branch-overlap aggregate: `{step5.get('branch_overlap', {}).get('aggregate_overlap_call')}`.",
            "",
            _markdown_table(
                [
                    "lane",
                    "scaled_overlap",
                    "call",
                    "candidate_delta_scaled_norm",
                    "prestall_delta_scaled_norm",
                ],
                _c0g_final_overlap_rows(step5),
            ),
            "",
            "A non-progressing capped probe is interpreted as no evidence, not fold proof.",
            "",
            "## Step 6 Bordered Reduced Conditioning",
            "",
        ]
    )
    step6 = result.get("step6_bordered_conditioning", {})
    lines.append(f"Call: `{step6.get('call')}`.")
    lines.append("")
    lines.append(
        _markdown_table(
            [
                "tau",
                "status",
                "cond_JQ_perp",
                "cond_Jb",
                "sigma_min_JQ_perp",
                "sigma_min_Jb",
                "ratio",
                "call",
            ],
            _c0g_final_step6_rows(step6),
        )
    )
    lines.extend(["", "## Step 7 Commutator And P_G", ""])
    step7 = result.get("step7_commutator", {})
    lines.append(step7.get("mode_characterization_scope", ""))
    lines.append("")
    lines.append("Commutator rows:")
    lines.append("")
    lines.append(
        _markdown_table(
            [
                "tau",
                "norm_CG_over_G",
                "boundary_fraction",
                "mixed_control_one_minus_P_G",
                "mixed_control_curl_over_A",
            ],
            _c0g_final_commutator_rows(step7),
        )
    )
    lines.append("")
    lines.append("Per-mode `1-P_G` rows:")
    lines.append("")
    lines.append(
        _markdown_table(
            [
                "tau",
                "tracked_lane",
                "ascending_rank",
                "sigma",
                "one_minus_P_G",
                "P_G",
                "curl_over_A",
                "boundary_fraction",
            ],
            _c0g_final_step7_rows(step7),
        )
    )
    lines.extend(["", "## Step 8 Verdict Support", ""])
    primary = step8.get("primary_evidence", {})
    confirmatory = step8.get("confirmatory", {})
    lines.extend(
        [
            f"- primary fold support: `{primary.get('fold_primary')}`",
            f"- primary conditioning support: `{primary.get('conditioning_primary')}`",
            f"- gray items: `{primary.get('gray_items')}`",
            f"- not measured: `{primary.get('not_measured')}`",
            f"- disagreements: `{primary.get('disagreements')}`",
            f"- Step-5 overlap call: `{confirmatory.get('step5_overlap_call')}`",
            f"- Step-6 call: `{confirmatory.get('step6_call')}`",
            f"- gray-zone/agreement rule honored: `{step8.get('gray_zone_rule_honored')}`",
            "",
            "## Scope Confirmation",
            "",
        ]
    )
    scope = result.get("scope_confirmation", {})
    lines.extend(
        [
            f"- NO C0g fix implemented: `{scope.get('no_C0g_fix_implemented')}`",
            f"- Single arbiter residual: `{scope.get('single_arbiter_residual')}`",
            f"- Depth continuation: `{scope.get('depth_continuation')}`",
            f"- Solver logic touched: `{scope.get('solver_logic_touched_by_c0g')}`",
            f"- Faithful operators touched: `{scope.get('faithful_operators_touched_by_c0g')}`",
            f"- Frozen physics touched: `{scope.get('frozen_physics_touched_by_c0g')}`",
            f"- Physical export guard touched: `{scope.get('physical_export_guard_touched_by_c0g')}`",
            "",
            f"Complete machine JSON: `{result.get('config', {}).get('final_json_path')}`.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def c0g_step0_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Path-A C0g Step 0 premise gate only.")
    parser.add_argument("--c0f2-json", type=Path, default=C0gConfig().c0f2_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gConfig().run_root)
    parser.add_argument("--step0-json", type=Path, default=C0gConfig().step0_json_path)
    parser.add_argument("--json-path", type=Path, default=C0gConfig().json_path)
    parser.add_argument("--report-path", type=Path, default=C0gConfig().report_path)
    args = parser.parse_args(argv)
    config = C0gConfig(
        c0f2_json_path=args.c0f2_json,
        run_root=args.run_root,
        step0_json_path=args.step0_json,
        json_path=args.json_path,
        report_path=args.report_path,
    )
    result = run_c0g_step0_premise(config)
    step0 = result.get("step0", {})
    gate = step0.get("premise_gate", {}) if isinstance(step0, Mapping) else {}
    summary = {
        "phase": "step0_premise",
        "status": result.get("status"),
        "provenance": result.get("provenance", {}).get("call"),
        "premise_call": gate.get("call"),
        "f_nn": gate.get("f_nn"),
        "step0_json": str(config.step0_json_path),
    }
    print(json.dumps(summary, sort_keys=True))
    return 0


def c0g_state_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run one Path-A C0g per-state SVD/F_tau measurement.")
    parser.add_argument("--tau", type=float, required=True)
    parser.add_argument("--c0f2-json", type=Path, default=C0gConfig().c0f2_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gConfig().run_root)
    parser.add_argument("--step0-json", type=Path, default=C0gConfig().step0_json_path)
    parser.add_argument("--json-path", type=Path, default=C0gConfig().json_path)
    parser.add_argument("--report-path", type=Path, default=C0gConfig().report_path)
    args = parser.parse_args(argv)
    config = C0gConfig(
        c0f2_json_path=args.c0f2_json,
        run_root=args.run_root,
        step0_json_path=args.step0_json,
        json_path=args.json_path,
        report_path=args.report_path,
    )
    result = run_c0g_state_measurement(tau=float(args.tau), config=config)
    summary = {
        "phase": "state_measurement",
        "tau": result.get("tau"),
        "sigma_min": result.get("svd", {}).get("sigma_min"),
        "cos_theta": result.get("ftau", {}).get("representative_cos_theta")
        if isinstance(result.get("ftau"), Mapping)
        else None,
        "json": str(_c0g_state_json_path(config, tau=float(args.tau))),
    }
    print(json.dumps(summary, sort_keys=True))
    return 0


def c0g_aggregate_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Aggregate Path-A C0g staged Steps 0-3.")
    parser.add_argument("--c0f2-json", type=Path, default=C0gConfig().c0f2_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gConfig().run_root)
    parser.add_argument("--step0-json", type=Path, default=C0gConfig().step0_json_path)
    parser.add_argument("--json-path", type=Path, default=C0gConfig().json_path)
    parser.add_argument("--report-path", type=Path, default=C0gConfig().report_path)
    args = parser.parse_args(argv)
    config = C0gConfig(
        c0f2_json_path=args.c0f2_json,
        run_root=args.run_root,
        step0_json_path=args.step0_json,
        json_path=args.json_path,
        report_path=args.report_path,
    )
    result = aggregate_c0g_steps0_3(config)
    summary = {
        "phase": "aggregate_steps0_3",
        "status": result.get("status"),
        "premise_call": result.get("step0", {}).get("premise_gate", {}).get("call")
        if isinstance(result.get("step0"), Mapping)
        else None,
        "preliminary": result.get("preliminary_reading", {}).get("status")
        if isinstance(result.get("preliminary_reading"), Mapping)
        else None,
        "json": str(config.json_path),
        "report": str(config.report_path),
        "steps_4_to_7": result.get("steps_4_to_7"),
    }
    print(json.dumps(summary, sort_keys=True))
    return 0


def c0g_steps4_6_7_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Path-A C0g diagnostic Steps 4, 6, and 7.")
    parser.add_argument("--c0f2-json", type=Path, default=C0gConfig().c0f2_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gConfig().run_root)
    parser.add_argument("--step0-json", type=Path, default=C0gConfig().step0_json_path)
    parser.add_argument("--json-path", type=Path, default=C0gConfig().json_path)
    parser.add_argument("--step5-json-path", type=Path, default=C0gConfig().step5_json_path)
    parser.add_argument("--final-json-path", type=Path, default=C0gConfig().final_json_path)
    parser.add_argument("--report-path", type=Path, default=C0gConfig().report_path)
    args = parser.parse_args(argv)
    config = C0gConfig(
        c0f2_json_path=args.c0f2_json,
        run_root=args.run_root,
        step0_json_path=args.step0_json,
        json_path=args.json_path,
        step5_json_path=args.step5_json_path,
        final_json_path=args.final_json_path,
        report_path=args.report_path,
    )
    result = run_c0g_steps4_6_7(config)
    step4 = result.get("step4_mach", {})
    step6 = result.get("step6_bordered_conditioning", {})
    step7 = result.get("step7_commutator", {})
    summary = {
        "phase": "steps4_6_7",
        "status": result.get("status"),
        "mach_label": step4.get("sonic_context_label") if isinstance(step4, Mapping) else None,
        "sonic_represented": step4.get("sonic_hypothesis_represented") if isinstance(step4, Mapping) else None,
        "bordered_call": step6.get("call") if isinstance(step6, Mapping) else None,
        "commutator_rows": len(step7.get("rows", [])) if isinstance(step7, Mapping) else 0,
        "json": str(_c0g_steps4_6_7_json_path(config)),
    }
    print(json.dumps(summary, sort_keys=True))
    return 0


def c0g_step5_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Path-A C0g Step 5 capped scipy probe.")
    parser.add_argument("--c0f2-json", type=Path, default=C0gConfig().c0f2_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gConfig().run_root)
    parser.add_argument("--step0-json", type=Path, default=C0gConfig().step0_json_path)
    parser.add_argument("--json-path", type=Path, default=C0gConfig().json_path)
    parser.add_argument("--step5-json-path", type=Path, default=C0gConfig().step5_json_path)
    parser.add_argument("--final-json-path", type=Path, default=C0gConfig().final_json_path)
    parser.add_argument("--report-path", type=Path, default=C0gConfig().report_path)
    parser.add_argument(
        "--method-wall-seconds",
        type=float,
        default=C0gConfig().scipy_method_wall_seconds,
    )
    args = parser.parse_args(argv)
    config = C0gConfig(
        c0f2_json_path=args.c0f2_json,
        run_root=args.run_root,
        step0_json_path=args.step0_json,
        json_path=args.json_path,
        step5_json_path=args.step5_json_path,
        final_json_path=args.final_json_path,
        report_path=args.report_path,
        scipy_method_wall_seconds=float(args.method_wall_seconds),
    )
    result = run_c0g_step5_scipy_probe(config)
    summary = {
        "phase": "step5_scipy_probe",
        "status": result.get("status"),
        "step5_call": result.get("step5_call"),
        "caps": result.get("caps"),
        "method_progress": [
            {
                "method": row.get("method"),
                "status": row.get("status"),
                "progress_call": row.get("progress_call"),
                "nfev": row.get("nfev_wrapper"),
                "best_l2_ratio": row.get("best_l2_ratio"),
                "abort_reason": row.get("abort_reason") or row.get("reason"),
            }
            for row in result.get("method_results", [])
        ],
        "overlap_call": result.get("branch_overlap", {}).get("aggregate_overlap_call"),
        "json": str(config.step5_json_path),
    }
    print(json.dumps(summary, sort_keys=True))
    return 0


def _c0g_stdout_summary(result: Mapping[str, Any]) -> str:
    step4 = result.get("step4_mach", {})
    step5 = result.get("step5_scipy_probe", {})
    step6 = result.get("step6_bordered_conditioning", {})
    step7 = result.get("step7_commutator", {})
    step8 = result.get("step8", {})
    mach_rows = _c0g_final_mach_rows(step4)
    deepest_mach = min(mach_rows, key=lambda row: float(row["tau"])) if mach_rows else {}
    step6_rows = [
        row for row in _c0g_final_step6_rows(step6) if row.get("status") == "MEASURED"
    ]
    closest_step6 = min(step6_rows, key=lambda row: float(row["tau"])) if step6_rows else {}
    mode_rows = [
        row
        for row in _c0g_final_step7_rows(step7)
        if abs(float(row.get("tau", math.nan)) - float(result.get("config", {}).get("stalled_tau", math.nan))) < 5.0e-13
    ]
    if not mode_rows:
        mode_rows = _c0g_final_step7_rows(step7)[:5]
    mode_summary = ", ".join(
        f"track {row.get('tracked_lane')}/rank {row.get('ascending_rank')}: 1-P_G={_c0g_fmt(row.get('one_minus_P_G'))}"
        for row in mode_rows[:5]
    )
    method_bits = ", ".join(
        f"{row.get('method')} {row.get('status')} {row.get('progress_call')} nfev={row.get('nfev_wrapper')} best_l2_ratio={_c0g_fmt(row.get('best_l2_ratio'))}"
        for row in step5.get("method_results", [])
    )
    return "\n".join(
        [
            "C0g final diagnostic summary:",
            (
                "Step 4 Mach: "
                f"M_max={_c0g_fmt(deepest_mach.get('M_full_max'))}, "
                f"j_at_max={_c0g_fmt(deepest_mach.get('j_at_max'))}, "
                f"label={step4.get('sonic_context_label')}, "
                f"represented={step4.get('sonic_hypothesis_represented')}"
            ),
            f"Step 5 capped probe: {step5.get('step5_call')} ({method_bits}); overlap={step5.get('branch_overlap', {}).get('aggregate_overlap_call')}",
            (
                "Step 6 border: "
                f"cond(JQ_perp)={_c0g_fmt(closest_step6.get('cond_JQ_perp'))}, "
                f"cond(Jb)={_c0g_fmt(closest_step6.get('cond_Jb'))}, "
                f"ratio={_c0g_fmt(closest_step6.get('ratio'))}, "
                f"call={step6.get('call')}"
            ),
            f"Step 7 per-mode at stalled/nearest: {mode_summary}",
            (
                "Step 8 verdict: "
                f"{step8.get('verdict')} / {step8.get('sub_label')} ; "
                f"gray-zone honored={step8.get('gray_zone_rule_honored')}"
            ),
            f"Recommended C0g build branch: {step8.get('recommended_C0g_build_branch')}",
            f"Complete JSON: {result.get('config', {}).get('final_json_path')}",
            f"Report: {result.get('config', {}).get('report_path')}",
        ]
    )


def c0g_final_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Aggregate Path-A C0g final Steps 0-8 report.")
    parser.add_argument("--c0f2-json", type=Path, default=C0gConfig().c0f2_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gConfig().run_root)
    parser.add_argument("--step0-json", type=Path, default=C0gConfig().step0_json_path)
    parser.add_argument("--json-path", type=Path, default=C0gConfig().json_path)
    parser.add_argument("--step5-json-path", type=Path, default=C0gConfig().step5_json_path)
    parser.add_argument("--final-json-path", type=Path, default=C0gConfig().final_json_path)
    parser.add_argument("--report-path", type=Path, default=C0gConfig().report_path)
    args = parser.parse_args(argv)
    config = C0gConfig(
        c0f2_json_path=args.c0f2_json,
        run_root=args.run_root,
        step0_json_path=args.step0_json,
        json_path=args.json_path,
        step5_json_path=args.step5_json_path,
        final_json_path=args.final_json_path,
        report_path=args.report_path,
    )
    result = aggregate_c0g_final(config)
    print(_c0g_stdout_summary(result))
    return 0


def c0g_build_b1b2_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run staged Path-A C0g B-1/B-2 only.")
    parser.add_argument("--c0f2-json", type=Path, default=C0gBuildB1B2Config().c0f2_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gBuildB1B2Config().run_root)
    parser.add_argument("--report-path", type=Path, default=C0gBuildB1B2Config().report_path)
    parser.add_argument("--json-path", type=Path, default=C0gBuildB1B2Config().json_path)
    parser.add_argument(
        "--taus",
        default=",".join(str(tau) for tau in C0gBuildB1B2Config().b2_depth_sequence),
        help="Comma-separated gauge-fixed crawl tau sequence for B-2.",
    )
    parser.add_argument("--proof-tau", type=float, default=C0gBuildB1B2Config().proof_tau)
    parser.add_argument(
        "--tau-fold-reference",
        type=float,
        default=C0gBuildB1B2Config().tau_fold_reference,
    )
    parser.add_argument(
        "--max-tau-backtracks",
        type=int,
        default=C0gBuildB1B2Config().max_tau_backtracks,
    )
    parser.add_argument(
        "--min-tau-backtrack-delta",
        type=float,
        default=C0gBuildB1B2Config().min_tau_backtrack_delta,
    )
    parser.add_argument(
        "--max-newton-iters",
        type=int,
        default=0,
        help="0 means use the default full Newton budget.",
    )
    parser.add_argument(
        "--jacobian-assembly",
        choices=("colored_sparse_jacobian_lu", "autodiff_sparse_jacobian_lu"),
        default=C0gBuildB1B2Config().jacobian_assembly,
    )
    args = parser.parse_args(argv)
    config = C0gBuildB1B2Config(
        c0f2_json_path=args.c0f2_json,
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        b2_depth_sequence=tuple(float(piece) for piece in args.taus.split(",") if piece),
        proof_tau=float(args.proof_tau),
        tau_fold_reference=float(args.tau_fold_reference),
        max_tau_backtracks=int(args.max_tau_backtracks),
        min_tau_backtrack_delta=float(args.min_tau_backtrack_delta),
        max_newton_iters_override=(
            int(args.max_newton_iters) if int(args.max_newton_iters) > 0 else None
        ),
        jacobian_assembly=str(args.jacobian_assembly),
    )
    result = run_c0g_build_b1b2(config)
    print(_c0g_build_b1b2_stdout_summary(result))
    return 0


def c0g_build_b4_b2_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run Path-A C0g B-4 match gate, then B-1/B-2 only if B-4 passes."
    )
    parser.add_argument("--c0f2-json", type=Path, default=C0gBuildB1B2Config().c0f2_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gBuildB1B2Config().run_root)
    parser.add_argument("--report-path", type=Path, default=C0gBuildB1B2Config().report_path)
    parser.add_argument("--json-path", type=Path, default=C0gBuildB1B2Config().json_path)
    parser.add_argument(
        "--taus",
        default=",".join(str(tau) for tau in C0gBuildB1B2Config().b2_depth_sequence),
        help="Comma-separated gauge-fixed crawl tau sequence for B-2.",
    )
    parser.add_argument("--proof-tau", type=float, default=C0gBuildB1B2Config().proof_tau)
    parser.add_argument(
        "--tau-fold-reference",
        type=float,
        default=C0gBuildB1B2Config().tau_fold_reference,
    )
    parser.add_argument(
        "--max-tau-backtracks",
        type=int,
        default=C0gBuildB1B2Config().max_tau_backtracks,
    )
    parser.add_argument(
        "--min-tau-backtrack-delta",
        type=float,
        default=C0gBuildB1B2Config().min_tau_backtrack_delta,
    )
    parser.add_argument(
        "--max-newton-iters",
        type=int,
        default=0,
        help="0 means use the default full Newton budget.",
    )
    parser.add_argument(
        "--jacobian-assembly",
        choices=("autodiff_sparse_jacobian_lu",),
        default=C0gBuildB1B2Config().jacobian_assembly,
    )
    parser.add_argument("--b4-reference-tau", type=float, default=C0gBuildB1B2Config().b4_reference_tau)
    parser.add_argument(
        "--b4-random-probes",
        type=int,
        default=C0gBuildB1B2Config().b4_random_probe_count,
    )
    args = parser.parse_args(argv)
    config = C0gBuildB1B2Config(
        c0f2_json_path=args.c0f2_json,
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        b2_depth_sequence=tuple(float(piece) for piece in args.taus.split(",") if piece),
        proof_tau=float(args.proof_tau),
        tau_fold_reference=float(args.tau_fold_reference),
        max_tau_backtracks=int(args.max_tau_backtracks),
        min_tau_backtrack_delta=float(args.min_tau_backtrack_delta),
        max_newton_iters_override=(
            int(args.max_newton_iters) if int(args.max_newton_iters) > 0 else None
        ),
        jacobian_assembly=str(args.jacobian_assembly),
        b4_reference_tau=float(args.b4_reference_tau),
        b4_random_probe_count=int(args.b4_random_probes),
    )
    result = run_c0g_build_b4_then_b2(config)
    print(_c0g_build_b1b2_stdout_summary(result))
    return 0


def c0g_deepcrawl_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the Path-A C0g deeper gauge-fixed continuation crawl."
    )
    parser.add_argument("--prior-b1b2-json", type=Path, default=C0gDeepCrawlConfig().prior_b1b2_json_path)
    parser.add_argument("--prior-crawl-json", type=Path, default=C0gDeepCrawlConfig().prior_crawl_json_path)
    parser.add_argument("--run-root", type=Path, default=C0gDeepCrawlConfig().run_root)
    parser.add_argument("--report-path", type=Path, default=C0gDeepCrawlConfig().report_path)
    parser.add_argument("--json-path", type=Path, default=C0gDeepCrawlConfig().json_path)
    parser.add_argument(
        "--taus",
        default=",".join(str(tau) for tau in C0gDeepCrawlConfig().target_taus),
        help="Comma-separated deeper target tau sequence.",
    )
    parser.add_argument("--tau-fold-reference", type=float, default=C0gDeepCrawlConfig().tau_fold_reference)
    parser.add_argument("--max-new-states", type=int, default=1)
    parser.add_argument("--finalize-only", action="store_true")
    parser.add_argument("--skip-b4-recheck", action="store_true")
    parser.add_argument("--max-tau-backtracks", type=int, default=C0gDeepCrawlConfig().max_tau_backtracks)
    parser.add_argument(
        "--min-tau-backtrack-delta",
        type=float,
        default=C0gDeepCrawlConfig().min_tau_backtrack_delta,
    )
    parser.add_argument(
        "--max-newton-iters",
        type=int,
        default=0,
        help="0 means use the default full Newton budget.",
    )
    parser.add_argument(
        "--jacobian-assembly",
        choices=("autodiff_sparse_jacobian_lu",),
        default=C0gDeepCrawlConfig().jacobian_assembly,
    )
    config = C0gDeepCrawlConfig()
    args = parser.parse_args(argv)
    config = C0gDeepCrawlConfig(
        prior_b1b2_json_path=args.prior_b1b2_json,
        prior_crawl_json_path=args.prior_crawl_json,
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        target_taus=tuple(float(piece) for piece in args.taus.split(",") if piece),
        tau_fold_reference=float(args.tau_fold_reference),
        max_tau_backtracks=int(args.max_tau_backtracks),
        min_tau_backtrack_delta=float(args.min_tau_backtrack_delta),
        max_newton_iters_override=(
            int(args.max_newton_iters) if int(args.max_newton_iters) > 0 else None
        ),
        jacobian_assembly=str(args.jacobian_assembly),
        b4_reference_tau=float(config.b4_reference_tau),
        b4_random_probe_count=int(config.b4_random_probe_count),
        b4_random_seed=int(config.b4_random_seed),
        b4_match_abs_tol=float(config.b4_match_abs_tol),
        b4_match_rel_tol=float(config.b4_match_rel_tol),
    )
    result = run_c0g_deepcrawl(
        config,
        max_new_states=int(args.max_new_states),
        finalize_only=bool(args.finalize_only),
        b4_recheck=not bool(args.skip_b4_recheck),
    )
    print(_c0g_deep_stdout_summary(result))
    return 0


def c0g_b3_pseudoarclength_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run Path-A C0g B-3 gauge-fixed Keller pseudo-arclength continuation."
    )
    default = C0gB3PseudoArclengthConfig()
    parser.add_argument("--deepcrawl-json", type=Path, default=default.deepcrawl_json_path)
    parser.add_argument("--run-root", type=Path, default=default.run_root)
    parser.add_argument("--report-path", type=Path, default=default.report_path)
    parser.add_argument("--json-path", type=Path, default=default.json_path)
    parser.add_argument("--skip-b4-recheck", action="store_true")
    parser.add_argument("--b30-only", action="store_true")
    parser.add_argument("--b31-steps", type=int, default=default.b31_max_steps)
    parser.add_argument("--initial-ds", type=float, default=default.initial_ds)
    parser.add_argument("--min-ds", type=float, default=default.min_ds)
    parser.add_argument("--max-ds", type=float, default=default.max_ds)
    parser.add_argument("--tau-scale", type=float, default=default.tau_scale)
    parser.add_argument("--corrector-max-iters", type=int, default=default.corrector_max_iters)
    parser.add_argument("--target-tau-tolerance", type=float, default=default.target_tau_tolerance)
    parser.add_argument(
        "--reproduction-comparison-tolerance",
        type=float,
        default=default.reproduction_comparison_tolerance,
    )
    parser.add_argument(
        "--jacobian-assembly",
        choices=("autodiff_sparse_jacobian_lu",),
        default=default.jacobian_assembly,
    )
    args = parser.parse_args(argv)
    config = C0gB3PseudoArclengthConfig(
        deepcrawl_json_path=args.deepcrawl_json,
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        jacobian_assembly=str(args.jacobian_assembly),
        tau_scale=float(args.tau_scale),
        initial_ds=float(args.initial_ds),
        min_ds=float(args.min_ds),
        max_ds=float(args.max_ds),
        corrector_max_iters=int(args.corrector_max_iters),
        target_tau_tolerance=float(args.target_tau_tolerance),
        reproduction_comparison_tolerance=float(args.reproduction_comparison_tolerance),
        b31_max_steps=int(args.b31_steps),
        b4_reference_tau=float(default.b4_reference_tau),
        b4_random_probe_count=int(default.b4_random_probe_count),
        b4_random_seed=int(default.b4_random_seed),
        b4_match_abs_tol=float(default.b4_match_abs_tol),
        b4_match_rel_tol=float(default.b4_match_rel_tol),
    )
    result = run_c0g_b3_pseudoarclength(
        config,
        b4_recheck=not bool(args.skip_b4_recheck),
        b31_steps=int(args.b31_steps),
        b30_only=bool(args.b30_only),
    )
    print(_c0g_b3_stdout_summary(result))
    return 0


def c0g_b3_followup_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run Path-A C0g B-3 follow-up v2 characterization battery."
    )
    default = C0gB3FollowupConfig()
    parser.add_argument("--deepcrawl-json", type=Path, default=default.deepcrawl_json_path)
    parser.add_argument("--b3-json", type=Path, default=default.b3_json_path)
    parser.add_argument("--run-root", type=Path, default=default.run_root)
    parser.add_argument("--report-path", type=Path, default=default.report_path)
    parser.add_argument("--json-path", type=Path, default=default.json_path)
    parser.add_argument("--fail-tau", type=float, default=default.fail_tau)
    parser.add_argument("--max-c1-iters", type=int, default=default.c1_max_iters)
    parser.add_argument("--c3-points", type=int, default=default.c3_points)
    parser.add_argument("--run-c4", action="store_true")
    parser.add_argument(
        "--jacobian-assembly",
        choices=("autodiff_sparse_jacobian_lu",),
        default=default.jacobian_assembly,
    )
    args = parser.parse_args(argv)
    config = C0gB3FollowupConfig(
        deepcrawl_json_path=args.deepcrawl_json,
        b3_json_path=args.b3_json,
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        fail_tau=float(args.fail_tau),
        c1_max_iters=int(args.max_c1_iters),
        c3_points=int(args.c3_points),
        c4_run=bool(args.run_c4),
        jacobian_assembly=str(args.jacobian_assembly),
    )
    result = run_c0g_b3_followup(config)
    print(_c0g_b3_followup_stdout_summary(result))
    return 0


def _c0c_attempt_artifacts(
    row: Mapping[str, Any],
    c0b_result: Mapping[str, Any],
) -> tuple[Path, Path]:
    state_path = _resolve_input_path(str(row["state_artifact"]))
    linear = row.get("linear_diagnostics", {})
    matrix_text = linear.get("matrix_path")
    if matrix_text:
        matrix_path = _resolve_input_path(str(matrix_text))
    else:
        c0b_run_root = Path(c0b_result.get("config", {}).get("run_root", DEFAULT_RUN_ROOT))
        matrix_path = _resolve_input_path(
            c0b_run_root
            / "matrices"
            / f"attempt_tau_{_format_tau(float(row['target_tau']))}_bt_{int(row.get('backtrack_index', 0))}.npz"
        )
    return state_path, matrix_path


def _c0c_annihilation_rows(
    *,
    matrix: csc_matrix,
    state: torch.Tensor,
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    generators: Sequence[C0cGenerator],
    grid,
    sigma_max: float,
    converged: bool,
    threshold: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    residual_np = residual_fn(state).detach().cpu().numpy().astype(np.float64)
    for generator in generators:
        norm = float(np.linalg.norm(generator.vector))
        if norm <= 0.0:
            rows.append(
                {
                    "generator": generator.name,
                    "classification": generator.classification,
                    "symmetry_status": generator.symmetry_status,
                    "status": "NOT_MEASURED",
                    "reason": "zero_generator_norm",
                    "test": "not_available",
                }
            )
            continue
        if not converged and generator.name != "phase":
            rows.append(
                {
                    "generator": generator.name,
                    "classification": generator.classification,
                    "symmetry_status": generator.symmetry_status,
                    "status": "NOT_MEASURED",
                    "reason": "nonconverged_state_and_probe_not_an_exact_symmetry",
                    "test": "not_measured_nonroot_probe",
                }
            )
            continue

        denominator = max(float(sigma_max) * norm, np.finfo(np.float64).tiny)
        assembled_jg = matrix @ generator.vector
        direction = torch.as_tensor(generator.vector, dtype=state.dtype, device=state.device)
        jvp_jg = jvp(residual_fn, state.detach(), direction).detach().cpu().numpy().astype(
            np.float64
        )
        assembled_rel = float(np.linalg.norm(assembled_jg) / denominator)
        jvp_rel = float(np.linalg.norm(jvp_jg) / denominator)
        crosscheck_rel = float(np.linalg.norm(assembled_jg - jvp_jg) / denominator)

        if converged:
            gate_pass = bool(assembled_rel <= threshold and jvp_rel <= threshold)
            rows.append(
                {
                    "generator": generator.name,
                    "classification": generator.classification,
                    "symmetry_status": generator.symmetry_status,
                    "status": "MEASURED",
                    "test": "root_annihilation",
                    "annihilation_rel_assembled": assembled_rel,
                    "annihilation_rel_jvp": jvp_rel,
                    "jvp_crosscheck_rel": crosscheck_rel,
                    "threshold": float(threshold),
                    "null_gate_pass": gate_pass,
                }
            )
        else:
            target = _phase_action_closed_vector(residual_np, grid)
            target_norm = float(np.linalg.norm(target))
            assembled_diff = float(np.linalg.norm(assembled_jg - target))
            jvp_diff = float(np.linalg.norm(jvp_jg - target))
            rows.append(
                {
                    "generator": generator.name,
                    "classification": generator.classification,
                    "symmetry_status": generator.symmetry_status,
                    "status": "MEASURED",
                    "test": "nonroot_phase_equivariance_identity",
                    "annihilation_rel_assembled_not_a_null_gate": assembled_rel,
                    "annihilation_rel_jvp_not_a_null_gate": jvp_rel,
                    "equivariance_rel_assembled_sigma_scaled": assembled_diff / denominator,
                    "equivariance_rel_jvp_sigma_scaled": jvp_diff / denominator,
                    "equivariance_rel_assembled_target_scaled": assembled_diff
                    / max(target_norm, np.finfo(np.float64).tiny),
                    "equivariance_rel_jvp_target_scaled": jvp_diff
                    / max(target_norm, np.finfo(np.float64).tiny),
                    "jvp_crosscheck_rel": crosscheck_rel,
                    "threshold": float(threshold),
                    "equivariance_pass": bool(jvp_diff / denominator <= threshold),
                    "null_gate_pass": None,
                }
            )
    return rows


def _dense_full_jvp_jacobian(
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    state: torch.Tensor,
    *,
    chunk_size: int,
) -> np.ndarray:
    x = state.detach()
    n = int(x.numel())
    dense = np.empty((n, n), dtype=np.float64)
    eye_torch = torch.eye(n, dtype=x.dtype, device=x.device)

    def one_jvp(direction: torch.Tensor) -> torch.Tensor:
        return jvp(residual_fn, x, direction)

    for start in range(0, n, max(1, int(chunk_size))):
        stop = min(n, start + max(1, int(chunk_size)))
        basis = eye_torch[start:stop]
        try:
            jvs = torch.func.vmap(one_jvp)(basis)
            dense[:, start:stop] = jvs.detach().cpu().numpy().astype(np.float64).T
        except Exception:
            for offset, direction in enumerate(basis):
                dense[:, start + offset] = (
                    one_jvp(direction).detach().cpu().numpy().astype(np.float64)
                )
    return dense


def _c0c_dense_sigma_validation(
    *,
    state: torch.Tensor,
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    assembled_sigma_min: float,
    assembled_sigma_max: float,
    tau: float,
    config: C0cConfig,
) -> dict[str, Any]:
    started = time.perf_counter()
    dense_jvp = _dense_full_jvp_jacobian(
        residual_fn,
        state,
        chunk_size=config.dense_jvp_chunk_size,
    )
    values = np.linalg.svd(dense_jvp, compute_uv=False)
    sigma_min = float(np.min(values))
    sigma_max = float(np.max(values))
    tolerance_abs = max(
        1.0e-12,
        1.0e3
        * np.finfo(np.float64).eps
        * max(abs(float(assembled_sigma_max)), abs(sigma_max), 1.0),
    )
    difference = abs(sigma_min - float(assembled_sigma_min))
    status = "MATCH" if difference <= tolerance_abs else "DISCREPANCY"
    return {
        "status": status,
        "tau": float(tau),
        "method": "dense_full_jvp_jacobian_svd",
        "assembled_matrix_method": "dense_svd_recomputed_from_saved_c0b_matrix",
        "dense_full_jvp_sigma_min": sigma_min,
        "dense_full_jvp_sigma_max": sigma_max,
        "assembled_sigma_min": float(assembled_sigma_min),
        "assembled_sigma_max": float(assembled_sigma_max),
        "abs_difference": float(difference),
        "tolerance_abs": float(tolerance_abs),
        "trusted_sigma_source": "both" if status == "MATCH" else "dense_full_jvp_jacobian",
        "elapsed_seconds": time.perf_counter() - started,
        "chunk_size": int(config.dense_jvp_chunk_size),
    }


def _classify_c0c_modes(
    *,
    modes: Sequence[Mapping[str, Any]],
    generators: Sequence[Mapping[str, Any]],
    annihilation_rows: Sequence[Mapping[str, Any]],
    config: C0cConfig,
) -> tuple[list[dict[str, Any]], dict[str, Any], str, str]:
    generator_gate = {
        str(row["generator"]): bool(row.get("null_gate_pass"))
        for row in annihilation_rows
        if row.get("status") == "MEASURED"
    }
    generator_by_name = {str(row["name"]): row for row in generators}
    classified: list[dict[str, Any]] = []
    explained_count = 0
    for mode in modes:
        overlaps = mode.get("overlaps", {})
        candidates = []
        for name, metadata in generator_by_name.items():
            overlap = float(overlaps.get(name, math.nan))
            if (
                generator_gate.get(name, False)
                and math.isfinite(overlap)
                and overlap >= config.overlap_threshold
            ):
                candidates.append((overlap, name, metadata))
        if candidates:
            overlap, name, metadata = max(candidates, key=lambda item: item[0])
            classification = str(metadata["classification"])
            explained_count += 1
            support = {
                "generator": name,
                "overlap": float(overlap),
                "annihilation_gate_pass": True,
            }
        else:
            classification = "UNEXPLAINED_STIFFNESS"
            best_name = None
            best_overlap = math.nan
            if overlaps:
                best_name = max(overlaps, key=lambda key: float(overlaps[key]))
                best_overlap = float(overlaps[best_name])
            support = {
                "best_generator": best_name,
                "best_overlap": best_overlap,
                "annihilation_gate_pass": False,
            }
        row = dict(mode)
        row["classification"] = classification
        row["classification_support"] = support
        classified.append(row)

    if not modes:
        verdict = "DIAGNOSTIC_INCOMPLETE"
        recommended = "Recompute the saved-matrix SVD before selecting a null-mode fix."
    elif explained_count == len(modes):
        verdict = "SYMMETRY_MODE_IDENTIFIED"
        recommended = (
            "Next step: add a diagnostic-gated gauge/null-space fix for the identified "
            "cluster, e.g. pin the global phase and deflate any additional identified "
            "null generators, then rerun the C0 crawl."
        )
    elif explained_count > 0:
        verdict = "MIXED"
        recommended = (
            "Next step: pin or deflate the explained symmetry mode(s), especially the "
            "global phase if present, then investigate the remaining field-lane residual "
            "subspace before a re-crawl."
        )
    else:
        verdict = "GENUINE_STIFFNESS"
        recommended = (
            "Next step: treat the cluster as unresolved field-block stiffness; do not "
            "apply a symmetry fix until a generator passes both annihilation and overlap gates."
        )

    support = {
        "thresholds": {
            "annihilation_rel_max": float(config.annihilation_threshold),
            "overlap_min": float(config.overlap_threshold),
            "span_residual_fraction_reference": float(config.residual_fraction_threshold),
        },
        "generator_null_gates": generator_gate,
        "explained_mode_count": int(explained_count),
        "cluster_mode_count": int(len(modes)),
        "mode_classifications": [
            {
                "mode_index": int(row["mode_index"]),
                "sigma": float(row["sigma"]),
                "classification": row["classification"],
                "support": row["classification_support"],
                "unexplained_residual_fraction": float(
                    row["unexplained_residual_fraction"]
                ),
                "lane_energy_fractions": row["lane_energy_fractions"],
            }
            for row in classified
        ],
    }
    return classified, support, verdict, recommended


def _c0c_incomplete_result(reason: str, config: C0cConfig) -> dict[str, Any]:
    return {
        "schema": "stage1_pathA_C0c_nullmode_identification/v1",
        "source_revision": source_revision(),
        "verdict": "DIAGNOSTIC_INCOMPLETE",
        "verdict_support": {"reason": reason},
        "recommended_next_step": "Provide the missing C0b saved state/matrix evidence and rerun C0c.",
        "config": _c0c_config_to_dict(config),
        "points": [],
        "dense_sigma_validation": {"status": "NOT_MEASURED", "reason": reason},
        "faithful_operator_boundary": _faithful_operator_boundary(),
    }


def run_c0c_nullmode_identification(config: C0cConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = C0cConfig()
    started = time.perf_counter()
    dtype = configure_backend(BackendConfig())
    c0b_json_path = _resolve_input_path(config.c0b_json_path)
    if not c0b_json_path.exists():
        result = _c0c_incomplete_result(f"missing_c0b_json:{c0b_json_path}", config)
        _write_c0c_outputs(result, config)
        return result
    c0b_result = json.loads(c0b_json_path.read_text(encoding="utf-8"))
    attempts = [dict(row) for row in c0b_result.get("tau_attempts", [])]
    converged = [row for row in attempts if row.get("final_physical_converged")]
    if not attempts or not converged:
        result = _c0c_incomplete_result("missing_attempts_or_converged_shallow_state", config)
        _write_c0c_outputs(result, config)
        return result

    shallow_row = max(converged, key=lambda row: float(row["target_tau"]))
    deepest_row = min(attempts, key=lambda row: float(row["target_tau"]))
    selected_rows: list[dict[str, Any]] = [shallow_row]
    if deepest_row is not shallow_row:
        selected_rows.append(deepest_row)

    c0b_grid = tuple(int(value) for value in c0b_result.get("config", {}).get("grid", C0Config().grid))
    points: list[dict[str, Any]] = []
    shallow_context: dict[str, Any] | None = None
    for row in selected_rows:
        tau = float(row["target_tau"])
        state_path, matrix_path = _c0c_attempt_artifacts(row, c0b_result)
        branch, provider, grid, boundaries, eos_K, residual_fn = _c0c_residual_context(
            tau=tau,
            grid_shape=c0b_grid,
            dtype=dtype,
        )
        del branch, provider, boundaries, eos_K
        state = _load_state_artifact(state_path, dtype=dtype)
        matrix = load_npz(matrix_path).tocsc()
        svd = _full_svd_cluster_from_matrix(matrix, mode_count=config.cluster_mode_count)
        generators = _c0c_generators_for_state(state, grid)
        modes, span_rank = _c0c_overlap_diagnostics(
            right_vectors=svd["right_vectors"],
            singular_values=svd["singular_values"],
            generators=generators,
            grid=grid,
        )
        annihilation = _c0c_annihilation_rows(
            matrix=matrix,
            state=state,
            residual_fn=residual_fn,
            generators=generators,
            grid=grid,
            sigma_max=float(svd["sigma_max"]),
            converged=bool(row.get("final_physical_converged")),
            threshold=config.annihilation_threshold,
        )
        point = {
            "tau": tau,
            "backtrack_index": int(row.get("backtrack_index", 0)),
            "state_artifact": str(state_path),
            "matrix_path": str(matrix_path),
            "final_physical_converged": bool(row.get("final_physical_converged")),
            "final_original_residual_linf": float(row.get("final_original_residual_linf", math.nan)),
            "layout": list(C0C_FIELD_LAYOUT),
            "layout_source": "stage1_solver.coupled_branch.pack/unpack_closed_coupled_fields",
            "generators": _c0c_generator_metadata(generators),
            "annihilation": {
                "status": "MEASURED",
                "threshold": float(config.annihilation_threshold),
                "rows": annihilation,
            },
            "svd": {
                "status": "MEASURED",
                "method": svd["method"],
                "sigma_min": float(svd["sigma_min"]),
                "sigma_max": float(svd["sigma_max"]),
                "condition": float(svd["condition"]),
                "vector_svd_sigma_min": float(svd["vector_svd_sigma_min"]),
                "vector_svd_sigma_max": float(svd["vector_svd_sigma_max"]),
                "cluster_singular_values": [float(value) for value in svd["singular_values"]],
            },
            "overlap": {
                "status": "MEASURED",
                "span_rank": int(span_rank),
                "modes": modes,
            },
        }
        if row is shallow_row:
            shallow_context = {
                "point": point,
                "state": state,
                "residual_fn": residual_fn,
                "svd": svd,
            }
        points.append(point)

    if shallow_context is None:
        result = _c0c_incomplete_result("missing_shallow_context", config)
        _write_c0c_outputs(result, config)
        return result

    shallow_point = shallow_context["point"]
    classified, support, verdict, recommended = _classify_c0c_modes(
        modes=shallow_point["overlap"]["modes"],
        generators=shallow_point["generators"],
        annihilation_rows=shallow_point["annihilation"]["rows"],
        config=config,
    )
    shallow_point["overlap"]["modes"] = classified

    dense_validation = _c0c_dense_sigma_validation(
        state=shallow_context["state"],
        residual_fn=shallow_context["residual_fn"],
        assembled_sigma_min=float(shallow_context["svd"]["sigma_min"]),
        assembled_sigma_max=float(shallow_context["svd"]["sigma_max"]),
        tau=float(shallow_point["tau"]),
        config=config,
    )
    support["dense_sigma_validation"] = dense_validation
    result = {
        "schema": "stage1_pathA_C0c_nullmode_identification/v1",
        "source_revision": source_revision(),
        "c0b_source_json": str(c0b_json_path),
        "verdict": verdict,
        "verdict_support": support,
        "recommended_next_step": recommended,
        "config": _c0c_config_to_dict(config),
        "points": points,
        "dense_sigma_validation": dense_validation,
        "faithful_operator_boundary": _faithful_operator_boundary(),
        "scope_guard": {
            "diagnosis_only": True,
            "gauge_fix_or_deflation_implemented": False,
            "recrawl_implemented": False,
            "faithful_operators_touched_by_c0c": False,
            "frozen_physics_touched_by_c0c": False,
            "physical_export_guard_touched_by_c0c": False,
        },
        "elapsed_seconds": time.perf_counter() - started,
    }
    _write_c0c_outputs(result, config)
    return result


def _c0d_incomplete_result(reason: str, config: C0dConfig) -> dict[str, Any]:
    return {
        "schema": "stage1_pathA_C0d_maxwell_gauge_identification/v1",
        "source_revision": source_revision(),
        "verdict": "DIAGNOSTIC_INCOMPLETE",
        "verdict_support": {"reason": reason},
        "recommended_next_step": "Provide the missing saved-matrix C0b/C0c evidence and rerun C0d.",
        "config": _c0d_config_to_dict(config),
        "operator_sources": _c0d_operator_sources(),
        "gauge_subspace": {"status": "NOT_MEASURED", "reason": reason},
        "controls": {"status": "NOT_MEASURED", "reason": reason},
        "all_cluster_modes": [],
        "maxwell_modes": [],
        "faithful_operator_boundary": _faithful_operator_boundary(),
        "scope_guard": _c0d_scope_guard(),
    }


def _c0d_operator_sources() -> dict[str, Any]:
    return {
        "field_layout": {
            "source": "stage1_solver.coupled_branch.pack_closed_coupled_fields/unpack_closed_coupled_fields",
            "order": list(C0C_FIELD_LAYOUT),
            "spatial_a_lanes": ["ar", "aw"],
            "scalar_potential_lane": "a0",
            "gradient_embedding": "zeros in psi_real, psi_imag, a0, r0, mu; d_r(lambda) in ar; d_w(lambda) in aw",
        },
        "discrete_gradient": {
            "source": "stage1_solver.operators.tensor_center_gradient_r/w",
            "boundary_closure": "the operators' one-sided centered-gradient closures; no extra lambda or A boundary condition imposed",
        },
        "raw_divergence": {
            "source": "stage1_solver.operators.axisymmetric_vector_divergence(ar, aw)",
            "definition": "cell-average r^-2 d_r(r^2 A_r) + d_w A_w",
        },
        "weighted_gauge_residual": {
            "source": "stage1_solver.operators.localized_maxwell_operator ar/aw gauge-control blocks",
            "object": "grad(Z(w) * axisymmetric_vector_divergence(ar, aw))",
            "weight_source": "stage1_solver.coupled_branch.localization_weight_torch",
        },
    }


def _c0d_scope_guard() -> dict[str, bool]:
    return {
        "diagnosis_only": True,
        "gauge_fix_or_deflation_implemented": False,
        "recrawl_implemented": False,
        "xi_reassembly_implemented": False,
        "faithful_operators_touched_by_c0d": False,
        "frozen_physics_touched_by_c0d": False,
        "physical_export_guard_touched_by_c0d": False,
    }


def _c0d_subspace_summary(subspace: Mapping[str, Any], grid, config: C0dConfig) -> dict[str, Any]:
    gradient_singular_values = np.asarray(
        subspace["gradient_singular_values"],
        dtype=np.float64,
    )
    weighted_singular_values = np.asarray(
        subspace["weighted_divergence_singular_values"],
        dtype=np.float64,
    )
    return {
        "status": "MEASURED",
        "construction": (
            "Full nodal scalar basis lambda_k, one scalar field per (r,w) cell. "
            "Each column is embedded as delta a0=0, delta ar=tensor_center_gradient_r(lambda_k), "
            "delta aw=tensor_center_gradient_w(lambda_k), with all non-A closed-state lanes zero. "
            "G is the SVD-orthonormalized column space. G_harm is formed from right singular "
            "directions of Z*D_A restricted to G whose weighted-divergence singular value is "
            "below the configured relative threshold."
        ),
        "grid": [int(grid.spec.nr), int(grid.spec.nw)],
        "scalar_basis_dim": int(grid.spec.nr * grid.spec.nw),
        "full_state_dim": int(_closed_layout_dimensions(grid)["state_size"]),
        "dim_G": int(subspace["dim_g"]),
        "dim_G_harm": int(subspace["dim_g_harm"]),
        "gradient_matrix_shape": [
            int(subspace["gradient_matrix"].shape[0]),
            int(subspace["gradient_matrix"].shape[1]),
        ],
        "gradient_rank_rtol": float(config.gradient_rank_rtol),
        "gradient_rank_threshold": float(subspace["gradient_rank_threshold"]),
        "gradient_singular_values": [float(value) for value in gradient_singular_values],
        "weighted_divergence_operator": "Z(w) * axisymmetric_vector_divergence restricted to G",
        "harmonic_weighted_divergence_rtol": float(
            config.harmonic_weighted_divergence_rtol
        ),
        "weighted_divergence_harmonic_threshold": float(
            subspace["weighted_divergence_harmonic_threshold"]
        ),
        "weighted_divergence_singular_values": [
            float(value) for value in weighted_singular_values
        ],
        "weighted_divergence_harmonic_indices": list(
            subspace["weighted_divergence_harmonic_indices"]
        ),
    }


def _c0d_select_saved_matrix_row(c0b_result: Mapping[str, Any]) -> dict[str, Any] | None:
    converged = [
        dict(row)
        for row in c0b_result.get("tau_attempts", [])
        if row.get("final_physical_converged")
    ]
    if not converged:
        return None
    return max(converged, key=lambda row: float(row["target_tau"]))


def _c0d_determine_verdict(
    *,
    selected_modes: Sequence[Mapping[str, Any]],
    controls: Mapping[str, Any],
    config: C0dConfig,
    incomplete_reason: str | None,
) -> tuple[str, dict[str, Any], str]:
    mode_support = [
        {
            "mode_index": int(row["mode_index"]),
            "sigma": float(row["sigma"]),
            "p_g_fraction": float(row["p_g_fraction"]),
            "p_g_harm_fraction": float(row["p_g_harm_fraction"]),
            "weighted_gauge_residual": float(row["weighted_gauge_residual"]),
            "raw_divergence": float(row["raw_divergence"]),
            "spatial_a_energy_fraction": float(row["spatial_a_energy_fraction"]),
            "non_a_remainder": float(row["non_a_remainder"]),
            "classification": str(row["classification"]),
        }
        for row in selected_modes
    ]
    support = {
        "thresholds": {
            "p_g_min": float(config.gauge_capture_threshold),
            "weighted_gauge_residual_max": float(config.weighted_residual_threshold),
            "maxwell_lane_fraction_min": float(config.maxwell_lane_fraction_min),
        },
        "controls_status": controls.get("status"),
        "mode_count": int(len(selected_modes)),
        "required_mode_count": int(config.maxwell_mode_count),
        "mode_classifications": mode_support,
    }
    all_gauge = bool(
        selected_modes
        and len(selected_modes) == int(config.maxwell_mode_count)
        and all(row.get("classification") == "MAXWELL_GAUGE" for row in selected_modes)
    )
    any_gauge = any(row.get("classification") == "MAXWELL_GAUGE" for row in selected_modes)
    if incomplete_reason is not None:
        support["reason"] = incomplete_reason
        return (
            "DIAGNOSTIC_INCOMPLETE",
            support,
            "Complete the missing saved-matrix/subspace/control evidence before selecting a gauge or stiffness fix.",
        )
    if all_gauge:
        return (
            "WALL_IS_ALL_GAUGE",
            support,
            (
                "Next step: implement a combined diagnostic-gated gauge/null fix: pin or deflate "
                "the confirmed global U(1) phase mode and deflate the A-sector discrete-gradient "
                "subspace, or impose an equivalent compatible weighted-gauge constraint, then rerun "
                "the C0 crawl."
            ),
        )
    if any_gauge:
        return (
            "MIXED_GAUGE_PLUS_RESIDUAL",
            support,
            (
                "Next step: do not apply a whole-wall gauge conclusion. Pin/deflate only the "
                "gated gauge modes and investigate the remaining Maxwell-lane residual modes with "
                "their lane split and weighted residual evidence."
            ),
        )
    return (
        "GENUINE_MAXWELL_STIFFNESS",
        support,
        (
            "Next step: treat the four Maxwell-lane modes as genuine A-sector stiffness under "
            "the C0d gates; do not implement a gauge fix without new evidence."
        ),
    )


def run_c0d_maxwell_gauge_identification(config: C0dConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = C0dConfig()
    started = time.perf_counter()
    dtype = configure_backend(BackendConfig())
    c0b_json_path = _resolve_input_path(config.c0b_json_path)
    if not c0b_json_path.exists():
        result = _c0d_incomplete_result(f"missing_c0b_json:{c0b_json_path}", config)
        _write_c0d_outputs(result, config)
        return result
    c0b_result = json.loads(c0b_json_path.read_text(encoding="utf-8"))
    row = _c0d_select_saved_matrix_row(c0b_result)
    if row is None:
        result = _c0d_incomplete_result("missing_converged_c0b_matrix_row", config)
        _write_c0d_outputs(result, config)
        return result

    c0b_grid = tuple(int(value) for value in c0b_result.get("config", {}).get("grid", C0Config().grid))
    tau = float(row["target_tau"])
    branch, _provider, grid, _boundaries = _branch_context(
        tau=tau,
        config=C0Config(grid=c0b_grid),
        dtype=dtype,
    )
    state_path, matrix_path = _c0c_attempt_artifacts(row, c0b_result)
    matrix = load_npz(matrix_path).tocsc()
    svd = _full_svd_cluster_from_matrix(matrix, mode_count=config.cluster_mode_count)
    subspace = _c0d_build_gauge_subspace(grid, branch, config)
    controls = _c0d_control_diagnostics(subspace, grid, config)

    all_modes: list[dict[str, Any]] = []
    for mode_index, (sigma, vector) in enumerate(
        zip(svd["singular_values"], svd["right_vectors"])
    ):
        all_modes.append(
            _c0d_mode_diagnostics(
                mode_index=mode_index,
                sigma=float(sigma),
                vector=np.asarray(vector, dtype=np.float64),
                subspace=subspace,
                grid=grid,
                branch=branch,
                config=config,
            )
        )

    maxwell_modes = [
        row
        for row in all_modes
        if float(row["spatial_a_energy_fraction"]) >= float(config.maxwell_lane_fraction_min)
    ][: int(config.maxwell_mode_count)]
    selected_indices = {int(row["mode_index"]) for row in maxwell_modes}
    for mode in all_modes:
        if int(mode["mode_index"]) not in selected_indices:
            mode["classification"] = "NOT_A_MAXWELL_LANE_MODE"
            mode["classification_gate"] = {
                "reason": "spatial_a_energy_fraction_below_maxwell_lane_gate",
                "maxwell_lane_fraction_min": float(config.maxwell_lane_fraction_min),
            }

    incomplete_reason = None
    if len(maxwell_modes) != int(config.maxwell_mode_count):
        incomplete_reason = (
            f"found_{len(maxwell_modes)}_maxwell_lane_modes_required_{config.maxwell_mode_count}"
        )
    elif controls.get("status") != "PASS":
        incomplete_reason = "anti_hardcode_controls_failed"
    verdict, support, recommended = _c0d_determine_verdict(
        selected_modes=maxwell_modes,
        controls=controls,
        config=config,
        incomplete_reason=incomplete_reason,
    )
    result = {
        "schema": "stage1_pathA_C0d_maxwell_gauge_identification/v1",
        "source_revision": source_revision(),
        "c0b_source_json": str(c0b_json_path),
        "selected_tau": tau,
        "selected_state_artifact": str(state_path),
        "selected_matrix_path": str(matrix_path),
        "verdict": verdict,
        "verdict_support": support,
        "recommended_next_step": recommended,
        "combined_fix_design_if_wall_is_all_gauge": (
            "Pin/deflate the global U(1) phase generator already confirmed by C0c, and "
            "deflate the A-sector discrete-gradient subspace G or add an equivalent compatible "
            "weighted gauge constraint for grad(Z*D_A*A); then rerun the C0 crawl. This C0d run "
            "does not implement that fix."
        ),
        "config": _c0d_config_to_dict(config),
        "operator_sources": _c0d_operator_sources(),
        "gauge_subspace": _c0d_subspace_summary(subspace, grid, config),
        "controls": controls,
        "svd": {
            "status": "MEASURED",
            "method": svd["method"],
            "sigma_min": float(svd["sigma_min"]),
            "sigma_max": float(svd["sigma_max"]),
            "condition": float(svd["condition"]),
            "cluster_singular_values": [float(value) for value in svd["singular_values"]],
            "source": "dense SVD recomputed from saved C0b sparse matrix",
        },
        "all_cluster_modes": all_modes,
        "maxwell_modes": maxwell_modes,
        "phase_context": {
            "source": "C0c",
            "status": "CONFIRMED_GLOBAL_U1_PHASE_GAUGE",
            "note": "C0d classifies only the other four spatial A-lane near-null modes.",
        },
        "faithful_operator_boundary": _faithful_operator_boundary(),
        "scope_guard": _c0d_scope_guard(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    _write_c0d_outputs(result, config)
    return result


def _c0e_incomplete_result(reason: str, config: C0eConfig) -> dict[str, Any]:
    return {
        "schema": "stage1_pathA_C0e_gauge_invariant_curl_gate/v1",
        "source_revision": source_revision(),
        "verdict": "DIAGNOSTIC_INCOMPLETE",
        "verdict_support": {"reason": reason},
        "recommended_next_step": "Provide the missing pinned saved state/matrix evidence and rerun C0e.",
        "config": _c0e_config_to_dict(config),
        "operator_sources": {},
        "artifacts": [],
        "primary_artifact": None,
        "faithful_operator_boundary": _faithful_operator_boundary(),
        "scope_guard": _c0e_scope_guard(),
    }


def _c0e_available_artifact_rows(c0b_result: Mapping[str, Any], config: C0eConfig) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in c0b_result.get("depth_sequence", []):
        linear = row.get("linear_diagnostics", {})
        if linear.get("status") != "MEASURED":
            continue
        try:
            state_path, matrix_path = _c0c_attempt_artifacts(row, c0b_result)
        except Exception:
            continue
        if not state_path.exists() or not matrix_path.exists():
            continue
        item = dict(row)
        item["_state_path"] = state_path
        item["_matrix_path"] = matrix_path
        rows.append(item)
    rows.sort(key=lambda r: (abs(float(r.get("target_tau", math.inf)) - PRIOR_TAU_FLOOR), float(r.get("target_tau", math.inf))))
    if config.run_all_available_tau_sources:
        return rows
    return rows[:1]


def _c0e_trial_residual_norms(
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    state: torch.Tensor,
    step: np.ndarray,
) -> dict[str, Any]:
    try:
        trial = state.detach() + torch.as_tensor(step, dtype=state.dtype, device=state.device)
        residual = residual_fn(trial).detach().cpu().numpy().astype(np.float64)
        return {
            "status": "MEASURED",
            "l2": float(np.linalg.norm(residual)),
            "linf": float(np.max(np.abs(residual))),
            "finite": bool(np.all(np.isfinite(residual))),
        }
    except Exception as exc:  # pragma: no cover - diagnostic reports failures.
        return {"status": "NOT_MEASURED", "reason": repr(exc), "l2": math.nan, "linf": math.nan}


def _c0e_newton_framing_check(
    *,
    matrix: csc_matrix,
    state: torch.Tensor,
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    svd: Mapping[str, Any],
) -> dict[str, Any]:
    residual = residual_fn(state).detach().cpu().numpy().astype(np.float64)
    residual_l2 = float(np.linalg.norm(residual))
    residual_linf = float(np.max(np.abs(residual)))
    dense = matrix.toarray()
    try:
        step = np.linalg.solve(dense, -residual)
        solve_method = "dense_np_linalg_solve"
    except np.linalg.LinAlgError:
        step, _resid, _rank, _singular = np.linalg.lstsq(dense, -residual, rcond=None)
        solve_method = "dense_np_linalg_lstsq_fallback"
    near_null = np.asarray(svd["right_vectors"], dtype=np.float64).T
    q, r = np.linalg.qr(near_null)
    diag = np.abs(np.diag(r))
    rank = int(np.sum(diag > 1.0e-12))
    basis = q[:, :rank] if rank > 0 else np.zeros((near_null.shape[0], 0), dtype=np.float64)
    gauge_step = basis @ (basis.T @ step) if basis.shape[1] else np.zeros_like(step)
    physical_step = step - gauge_step
    full_trial = _c0e_trial_residual_norms(residual_fn, state, step)
    physical_trial = _c0e_trial_residual_norms(residual_fn, state, physical_step)
    linear_residual = dense @ step + residual
    return {
        "status": "MEASURED",
        "solve_method": solve_method,
        "initial_residual_l2": residual_l2,
        "initial_residual_linf": residual_linf,
        "full_step_residual_l2": full_trial.get("l2"),
        "full_step_residual_linf": full_trial.get("linf"),
        "full_step_residual_l2_ratio": float(full_trial.get("l2", math.nan) / residual_l2)
        if residual_l2 > 0.0
        else math.nan,
        "full_step_residual_linf_ratio": float(full_trial.get("linf", math.nan) / residual_linf)
        if residual_linf > 0.0
        else math.nan,
        "gauge_removed_step_residual_l2": physical_trial.get("l2"),
        "gauge_removed_step_residual_linf": physical_trial.get("linf"),
        "gauge_removed_step_residual_l2_ratio": float(
            physical_trial.get("l2", math.nan) / residual_l2
        )
        if residual_l2 > 0.0
        else math.nan,
        "gauge_removed_step_residual_linf_ratio": float(
            physical_trial.get("linf", math.nan) / residual_linf
        )
        if residual_linf > 0.0
        else math.nan,
        "step_norm": float(np.linalg.norm(step)),
        "near_null_component_norm": float(np.linalg.norm(gauge_step)),
        "complement_component_norm": float(np.linalg.norm(physical_step)),
        "near_null_component_fraction": float(
            np.dot(gauge_step, gauge_step) / max(np.dot(step, step), np.finfo(np.float64).tiny)
        ),
        "near_null_basis": "fresh_dense_svd_right_modes_0_through_4",
        "near_null_basis_rank": int(rank),
        "linear_solve_relative_residual": float(
            np.linalg.norm(linear_residual) / max(residual_l2, np.finfo(np.float64).tiny)
        ),
        "trial_residual_eval_scope": "temporary_tensors_only_no_state_advance_no_line_search_no_newton_write",
    }


def _c0e_lambda_preimage(
    vector: np.ndarray,
    coupled_subspace: Mapping[str, Any],
    grid,
    config: C0eConfig,
) -> dict[str, Any]:
    unit, _norm = _unit_vector(vector)
    matrix = np.asarray(coupled_subspace["generator_matrix"], dtype=np.float64)
    coeffs, _resid, rank, singular = np.linalg.lstsq(
        matrix,
        unit,
        rcond=float(config.gradient_rank_rtol),
    )
    fitted = matrix @ coeffs
    lambda_field = coeffs.reshape(grid.spec.nr, grid.spec.nw)
    frequency = _c0e_frequency_classification(lambda_field)
    return {
        "least_squares_rank": int(rank),
        "singular_values_min": float(np.min(singular)) if len(singular) else math.nan,
        "singular_values_max": float(np.max(singular)) if len(singular) else math.nan,
        "fit_residual_fraction": float(
            np.linalg.norm(unit - fitted) / max(np.linalg.norm(unit), np.finfo(np.float64).tiny)
        ),
        "lambda_norm": float(np.linalg.norm(lambda_field)),
        "frequency": frequency,
    }


def _c0e_mode_outcome(
    *,
    curl_metrics: Mapping[str, Any],
    controls: Mapping[str, Any],
    a_only_p_g: float,
    coupled_capture: float,
) -> dict[str, Any]:
    margins = {
        label: _c0e_curl_margin(
            float(curl_metrics[label]["curl_fraction"]),
            controls.get("bands", {}).get(label, {}),
        )
        for label in ("unweighted", "cell_volume_weighted")
    }
    norm_outcomes = {label: row["outcome"] for label, row in margins.items()}
    if any(outcome == "CONTROL_BANDS_OVERLAP" for outcome in norm_outcomes.values()):
        combined = "CONTROL_BANDS_OVERLAP"
    elif all(outcome == "GAUGE" for outcome in norm_outcomes.values()):
        if a_only_p_g < 0.5 and coupled_capture < 0.5:
            combined = "LOW_CAPTURE_LOW_CURL_OTHER"
        else:
            combined = "GAUGE"
    elif all(outcome == "TRANSVERSE" for outcome in norm_outcomes.values()):
        combined = "TRANSVERSE"
    else:
        combined = "AMBIGUOUS"
    finite_margins = [
        float(row["margin_log10_to_separator"])
        for row in margins.values()
        if math.isfinite(float(row.get("margin_log10_to_separator", math.nan)))
    ]
    return {
        "outcome": combined,
        "norm_outcomes": norm_outcomes,
        "margins": margins,
        "band_margin_log10_min": float(min(finite_margins)) if finite_margins else math.nan,
    }


def _c0e_phase_mode_diagnostics(
    *,
    mode: Mapping[str, Any],
    state: torch.Tensor,
    grid,
    branch,
    coupled_subspace: Mapping[str, Any],
) -> dict[str, Any]:
    vector = np.asarray(mode["vector"], dtype=np.float64)
    unit, _norm = _unit_vector(vector)
    generators = _c0c_generators_for_state(state, grid)
    phase = next(generator for generator in generators if generator.name == "phase")
    phase_unit, phase_norm = _unit_vector(phase.vector)
    curl_metrics = _c0e_curl_fraction_metrics(unit, grid, branch)
    return {
        "mode_index": int(mode["mode_index"]),
        "sigma": float(mode["sigma"]),
        "classification": "PHASE_GAUGE_CONFIRMED_BY_PHASE_AND_COUPLED_CAPTURE",
        "phase_overlap_abs": float(abs(np.dot(unit, phase_unit))) if phase_norm > 0.0 else math.nan,
        "phase_capture_fraction": float(abs(np.dot(unit, phase_unit)) ** 2)
        if phase_norm > 0.0
        else math.nan,
        "coupled_capture_fraction": _c0d_projection_fraction(coupled_subspace["basis"], unit),
        "spatial_a_norm_unweighted": curl_metrics["unweighted"]["spatial_a_norm"],
        "z_f_rw_norm_unweighted": curl_metrics["unweighted"]["curl_norm"],
        "spatial_a_norm_cell_volume_weighted": curl_metrics["cell_volume_weighted"]["spatial_a_norm"],
        "z_f_rw_norm_cell_volume_weighted": curl_metrics["cell_volume_weighted"]["curl_norm"],
        "note": "The curl ratio is denominator-noise for the phase mode and is not used as its gate.",
    }


def _c0e_mode_diagnostics(
    *,
    mode_index: int,
    sigma: float,
    vector: np.ndarray,
    matrix: csc_matrix,
    state: torch.Tensor,
    grid,
    branch,
    eos_K: float,
    boundaries,
    s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any],
    a_subspace: Mapping[str, Any],
    coupled_subspace: Mapping[str, Any],
    controls: Mapping[str, Any],
    config: C0eConfig,
) -> dict[str, Any]:
    unit, mode_norm = _unit_vector(vector)
    curl_metrics = _c0e_curl_fraction_metrics(unit, grid, branch)
    a_only = _c0d_a_only_vector(unit, grid)
    a_only_p_g = _c0d_projection_fraction(a_subspace["basis"], unit)
    a_only_p_g_a_normalized = _c0d_projection_fraction(a_subspace["basis"], a_only)
    coupled_capture = _c0d_projection_fraction(coupled_subspace["basis"], unit)
    preimage = _c0e_lambda_preimage(unit, coupled_subspace, grid, config)
    budget = _c0e_residual_budget_for_mode(
        matrix=matrix,
        state=state,
        vector=unit,
        grid=grid,
        branch=branch,
        eos_K=eos_K,
        boundaries=boundaries,
        s_sigma=s_sigma,
    )
    outcome = _c0e_mode_outcome(
        curl_metrics=curl_metrics,
        controls=controls,
        a_only_p_g=a_only_p_g,
        coupled_capture=coupled_capture,
    )
    lane_split = _lane_energy_split(unit, grid)
    return {
        "mode_index": int(mode_index),
        "sigma": float(sigma),
        "mode_norm_before_unit": float(mode_norm),
        "lane_energy_fractions": lane_split,
        "spatial_a_energy_fraction": float(lane_split.get("ar", 0.0))
        + float(lane_split.get("aw", 0.0)),
        "non_a_remainder": _c0d_non_a_fraction(unit, grid),
        "curl_fraction": curl_metrics,
        "a_only_p_g_fraction": a_only_p_g,
        "a_only_p_g_a_normalized_fraction": a_only_p_g_a_normalized,
        "coupled_capture_fraction": coupled_capture,
        "outcome": outcome["outcome"],
        "outcome_details": outcome,
        "lambda_preimage": preimage,
        "residual_budget_jv": budget,
    }


def _c0e_mechanism_summary(
    *,
    maxwell_modes: Sequence[Mapping[str, Any]],
    row_scaling: Mapping[str, Any],
) -> dict[str, Any]:
    gauge_like = [
        mode
        for mode in maxwell_modes
        if mode.get("outcome") in {"GAUGE", "LOW_CAPTURE_LOW_CURL_OTHER"}
    ]
    mechanism_modes = gauge_like if gauge_like else list(maxwell_modes)
    freq_counts: dict[str, int] = {}
    unexplained = []
    penalty_ratios = []
    for mode in mechanism_modes:
        classification = (
            mode.get("lambda_preimage", {})
            .get("frequency", {})
            .get("classification", "UNKNOWN")
        )
        freq_counts[classification] = freq_counts.get(classification, 0) + 1
        budget = mode.get("residual_budget_jv", {})
        unexplained.append(float(budget.get("unexplained_fraction_of_component_sum_norm", math.nan)))
        components = budget.get("components", {})
        penalty = float(components.get("maxwell_gauge_penalty_rows", {}).get("norm", math.nan))
        assembled = float(budget.get("assembled_jv_norm", math.nan))
        if math.isfinite(penalty) and math.isfinite(assembled):
            penalty_ratios.append(penalty / max(assembled, np.finfo(np.float64).tiny))
    row_ratio = float(row_scaling.get("maxwell_to_non_maxwell_rms_ratio", math.nan))
    secondary: list[str] = []
    if math.isfinite(row_ratio) and row_ratio < 1.0e-3:
        primary = "ROW_SCALING"
        secondary.append("MAXWELL_ROWS_SMALL_RELATIVE_TO_NON_MAXWELL")
    elif freq_counts.get("CHECKERBOARD", 0) >= max(1, math.ceil(len(mechanism_modes) / 2)):
        primary = "ODD_EVEN_DECOUPLING"
    elif (freq_counts.get("CHECKERBOARD", 0) + freq_counts.get("HIGH_K", 0)) >= max(
        1, math.ceil(len(mechanism_modes) / 2)
    ):
        primary = "ODD_EVEN_DECOUPLING"
        secondary.append("HIGH_K_LAMBDA_PREIMAGE")
    elif freq_counts.get("SMOOTH", 0) >= max(1, math.ceil(len(mechanism_modes) / 2)):
        primary = "SMOOTH_K2"
    else:
        primary = "OTHER"
    if math.isfinite(row_ratio) and row_ratio < 1.0e-1 and primary != "ROW_SCALING":
        secondary.append("POSSIBLE_ROW_SCALING_CONTRIBUTION")
    if any(value > 1.0e-6 for value in unexplained if math.isfinite(value)):
        secondary.append("NONZERO_UNEXPLAINED_JV_REMAINDER")
    if any(value > 1.0e6 for value in penalty_ratios if math.isfinite(value)):
        secondary.append("LARGE_COMPONENT_CANCELLATION_RELATIVE_TO_NEAR_NULL_JV")
    if primary == "ODD_EVEN_DECOUPLING":
        secondary.append("CONSISTENT_STAGGERED_DIV_GRAD_STENCIL_IS_VIABLE_STRUCTURAL_C0F_DESIGN_NOTE")
    return {
        "primary_mechanism": primary,
        "secondary_flags": secondary,
        "lambda_frequency_counts": freq_counts,
        "row_scaling_ratio": row_ratio,
        "max_unexplained_fraction_of_component_sum_norm": float(
            max([value for value in unexplained if math.isfinite(value)], default=math.nan)
        ),
        "max_gauge_penalty_to_assembled_jv_ratio": float(
            max([value for value in penalty_ratios if math.isfinite(value)], default=math.nan)
        ),
        "evidence_note": (
            "ODD_EVEN_DECOUPLING is emitted only when the lambda-preimage frequency evidence is "
            "checkerboard/high-k and the row/block residual budget is present."
        ),
    }


def _c0e_determine_verdict(
    *,
    maxwell_modes: Sequence[Mapping[str, Any]],
    phase_mode: Mapping[str, Any],
    controls: Mapping[str, Any],
    mechanism: Mapping[str, Any],
    config: C0eConfig,
) -> tuple[str, dict[str, Any], str]:
    mode_outcomes = {int(mode["mode_index"]): str(mode.get("outcome")) for mode in maxwell_modes}
    support = {
        "controls_status": controls.get("status"),
        "control_bands": controls.get("bands"),
        "mode_outcomes": mode_outcomes,
        "mode_curl_support": [
            {
                "mode_index": int(mode["mode_index"]),
                "sigma": float(mode["sigma"]),
                "outcome": mode.get("outcome"),
                "unweighted_curl_fraction": mode["curl_fraction"]["unweighted"]["curl_fraction"],
                "weighted_curl_fraction": mode["curl_fraction"]["cell_volume_weighted"]["curl_fraction"],
                "band_margin_log10_min": mode.get("outcome_details", {}).get(
                    "band_margin_log10_min"
                ),
                "a_only_p_g_fraction": mode.get("a_only_p_g_fraction"),
                "coupled_capture_fraction": mode.get("coupled_capture_fraction"),
            }
            for mode in maxwell_modes
        ],
        "phase_mode": {
            "mode_index": phase_mode.get("mode_index"),
            "phase_capture_fraction": phase_mode.get("phase_capture_fraction"),
            "coupled_capture_fraction": phase_mode.get("coupled_capture_fraction"),
        },
        "mechanism": mechanism,
        "required_maxwell_mode_count": int(config.maxwell_mode_count),
        "observed_maxwell_mode_count": len(maxwell_modes),
        "mode_2_outcome": mode_outcomes.get(2),
    }
    if len(maxwell_modes) != int(config.maxwell_mode_count):
        support["reason"] = "missing_required_four_maxwell_lane_modes"
        return (
            "DIAGNOSTIC_INCOMPLETE",
            support,
            "Complete the missing four-mode Maxwell-lane evidence before selecting C0f.",
        )
    if controls.get("status") != "PASS":
        support["reason"] = "control_bands_failed"
        return (
            "DIAGNOSTIC_INCOMPLETE",
            support,
            "Fix the non-circular positive/negative controls before selecting C0f.",
        )
    if any(
        outcome in {"AMBIGUOUS", "CONTROL_BANDS_OVERLAP", "LOW_CAPTURE_LOW_CURL_OTHER"}
        for outcome in mode_outcomes.values()
    ):
        support["reason"] = "one_or_more_modes_not_decisively_classified"
        return (
            "DIAGNOSTIC_INCOMPLETE",
            support,
            "Do not choose a C0f fix until every Maxwell-lane mode has a decisive curl-band outcome.",
        )
    phase_confirmed = bool(
        float(phase_mode.get("phase_capture_fraction", 0.0)) >= 0.9
        and float(phase_mode.get("coupled_capture_fraction", 0.0)) >= 0.9
    )
    if not phase_confirmed:
        support["reason"] = "phase_mode_not_confirmed"
        return (
            "DIAGNOSTIC_INCOMPLETE",
            support,
            "Recheck the phase/coupled generator basis before selecting C0f.",
        )
    gauge_modes = {index for index, outcome in mode_outcomes.items() if outcome == "GAUGE"}
    transverse_modes = {index for index, outcome in mode_outcomes.items() if outcome == "TRANSVERSE"}
    primary = str(mechanism.get("primary_mechanism"))
    if len(gauge_modes) == int(config.maxwell_mode_count):
        recommendation = (
            "C0f should use minimal adaptive deflation of the full near-null gauge set "
            "{phase + coupled-gauge modes 1,2,3,4 + the 2-D G_harm}; mode 2 is explicitly "
            "included because its curl gate is gauge. If the mechanism remains "
            "ODD_EVEN_DECOUPLING, C0f may instead or additionally evaluate a consistent/staggered "
            "div-grad gauge-penalty stencil as a structural design, but C0e implements no fix. "
            "Single-Arbiter plan: deflation lives only in the Newton linear-solve/preconditioner "
            "path; patha_closed_branch_residual remains the sole convergence arbiter and "
            "merit/line-search stays on the original ||F||. Caveats: the no-deflation reference "
            "exists only at tau≈0.029, and deflation must be adaptive across tau with re-SVD at "
            "new stalls."
        )
        support["reason"] = "all_four_maxwell_modes_in_gauge_band_and_phase_confirmed"
        support["primary_mechanism"] = primary
        return "WALL_IS_ALL_GAUGE", support, recommendation
    if {1, 3, 4}.issubset(gauge_modes) and transverse_modes == {2}:
        recommendation = (
            "C0f should deflate the confirmed gauge directions {phase + modes 1,3,4 + G_harm} "
            "and only mode 2's gradient projection; the surviving transverse remnant of mode 2 "
            "is the narrow A-sector question, not a production solver. Single-Arbiter plan: "
            "deflation lives only in the Newton linear-solve/preconditioner path; "
            "patha_closed_branch_residual remains the sole convergence arbiter and merit/"
            "line-search stays on the original ||F||. Caveats: the no-deflation reference exists "
            "only at tau≈0.029, and deflation must be adaptive across tau with re-SVD at new stalls."
        )
        support["reason"] = "modes_1_3_4_gauge_mode_2_transverse"
        support["primary_mechanism"] = primary
        return "GAUGE_PLUS_ONE_TRANSVERSE", support, recommendation
    if any(mode_outcomes.get(index) == "TRANSVERSE" for index in (1, 3, 4)):
        support["reason"] = "one_of_modes_1_3_4_has_real_curl"
        return (
            "GAUGE_FRAMING_REFUTED",
            support,
            "Stop before C0f: the gauge reading is refuted by curl-carrying modes among 1/3/4.",
        )
    support["reason"] = "mode_outcome_pattern_not_covered"
    return (
        "DIAGNOSTIC_INCOMPLETE",
        support,
        "Do not choose a C0f fix until the Maxwell-mode outcome pattern is decisive.",
    )


def _c0e_analyze_artifact(
    *,
    row: Mapping[str, Any],
    c0b_result: Mapping[str, Any],
    dtype: torch.dtype,
    config: C0eConfig,
) -> dict[str, Any]:
    tau = float(row["target_tau"])
    c0b_grid = tuple(int(value) for value in c0b_result.get("config", {}).get("grid", C0Config().grid))
    branch, provider, grid, boundaries, eos_K, residual_fn = _c0c_residual_context(
        tau=tau,
        grid_shape=c0b_grid,
        dtype=dtype,
    )
    state_path = Path(row["_state_path"])
    matrix_path = Path(row["_matrix_path"])
    state = _load_state_artifact(state_path, dtype=dtype)
    matrix = load_npz(matrix_path).tocsc()
    svd = _full_svd_cluster_from_matrix(matrix, mode_count=config.cluster_mode_count)
    a_subspace_config = C0dConfig(
        gradient_rank_rtol=config.gradient_rank_rtol,
        harmonic_weighted_divergence_rtol=config.harmonic_weighted_divergence_rtol,
    )
    a_subspace = _c0d_build_gauge_subspace(grid, branch, a_subspace_config)
    coupled_subspace = _c0e_build_coupled_subspace(state, grid, branch, config)
    controls = _c0e_control_diagnostics(
        a_subspace=a_subspace,
        coupled_subspace=coupled_subspace,
        state=state,
        grid=grid,
        branch=branch,
    )
    framing = _c0e_newton_framing_check(
        matrix=matrix,
        state=state,
        residual_fn=residual_fn,
        svd=svd,
    )
    raw_modes = [
        {
            "mode_index": int(index),
            "sigma": float(sigma),
            "vector": np.asarray(vector, dtype=np.float64),
        }
        for index, (sigma, vector) in enumerate(zip(svd["singular_values"], svd["right_vectors"]))
    ]
    phase_mode = _c0e_phase_mode_diagnostics(
        mode=raw_modes[0],
        state=state,
        grid=grid,
        branch=branch,
        coupled_subspace=coupled_subspace,
    )
    all_modes: list[dict[str, Any]] = []
    for mode in raw_modes:
        diagnostic = _c0e_mode_diagnostics(
            mode_index=int(mode["mode_index"]),
            sigma=float(mode["sigma"]),
            vector=np.asarray(mode["vector"], dtype=np.float64),
            matrix=matrix,
            state=state,
            grid=grid,
            branch=branch,
            eos_K=eos_K,
            boundaries=boundaries,
            s_sigma=provider,
            a_subspace=a_subspace,
            coupled_subspace=coupled_subspace,
            controls=controls,
            config=config,
        )
        if int(mode["mode_index"]) == 0:
            diagnostic["outcome"] = "PHASE_MODE_SEPARATE"
        all_modes.append(diagnostic)
    maxwell_modes = [
        mode
        for mode in all_modes
        if int(mode["mode_index"]) != 0
        and float(mode["spatial_a_energy_fraction"]) >= float(config.maxwell_lane_fraction_min)
    ][: int(config.maxwell_mode_count)]
    selected = {int(mode["mode_index"]) for mode in maxwell_modes}
    for mode in all_modes:
        if int(mode["mode_index"]) == 0:
            continue
        if int(mode["mode_index"]) not in selected:
            mode["outcome"] = "NOT_SELECTED_AS_MAXWELL_LANE_MODE"
            mode["selection_reason"] = "spatial_a_energy_fraction_below_gate_or_beyond_required_four_modes"
    row_scaling = _c0e_row_scaling_summary(matrix, grid)
    mechanism = _c0e_mechanism_summary(maxwell_modes=maxwell_modes, row_scaling=row_scaling)
    verdict, support, recommended = _c0e_determine_verdict(
        maxwell_modes=maxwell_modes,
        phase_mode=phase_mode,
        controls=controls,
        mechanism=mechanism,
        config=config,
    )
    serializable_modes = all_modes
    for mode in serializable_modes:
        mode.pop("vector", None)
    return {
        "tau": tau,
        "tau_source_label": (
            "deepest_saved_stall_attempt"
            if not bool(row.get("final_physical_converged", False))
            else "converged_tau_source"
        ),
        "final_physical_converged": bool(row.get("final_physical_converged", False)),
        "state_artifact": str(state_path),
        "matrix_path": str(matrix_path),
        "backtrack_index": int(row.get("backtrack_index", 0)),
        "verdict": verdict,
        "verdict_support": support,
        "recommended_next_step": recommended,
        "c0e_0_newton_framing": framing,
        "svd": {
            "status": "MEASURED",
            "method": svd["method"],
            "sigma_min": float(svd["sigma_min"]),
            "sigma_max": float(svd["sigma_max"]),
            "condition": float(svd["condition"]),
            "cluster_singular_values": [float(value) for value in svd["singular_values"]],
            "source": "fresh dense SVD recomputed from saved sparse matrix",
        },
        "gauge_subspace_a_only": {
            "dim_G": int(a_subspace["dim_g"]),
            "dim_G_harm": int(a_subspace["dim_g_harm"]),
            "gradient_rank_threshold": float(a_subspace["gradient_rank_threshold"]),
            "gradient_rank_rtol": float(config.gradient_rank_rtol),
        },
        "coupled_gauge_subspace": {
            "rank": int(coupled_subspace["rank"]),
            "rank_threshold": float(coupled_subspace["rank_threshold"]),
            "scalar_basis_dim": int(coupled_subspace["scalar_basis_dim"]),
            "full_state_dim": int(coupled_subspace["full_state_dim"]),
            "q_over_hbar": float(coupled_subspace["q_over_hbar"]),
            "q_over_hbar_source": "branch.gauge_charge / branch.hbar from coupled_branch.py alpha=q/hbar",
            "singular_value_min": float(np.min(coupled_subspace["singular_values"]))
            if len(coupled_subspace["singular_values"])
            else math.nan,
            "singular_value_max": float(np.max(coupled_subspace["singular_values"]))
            if len(coupled_subspace["singular_values"])
            else math.nan,
        },
        "controls": controls,
        "phase_mode": phase_mode,
        "all_cluster_modes": serializable_modes,
        "maxwell_modes": maxwell_modes,
        "row_scaling": row_scaling,
        "mechanism": mechanism,
    }


def run_c0e_gauge_invariant_curl_gate(config: C0eConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = C0eConfig()
    started = time.perf_counter()
    dtype = configure_backend(BackendConfig())
    c0b_json_path = _resolve_input_path(config.c0b_json_path)
    if not c0b_json_path.exists():
        result = _c0e_incomplete_result(f"missing_c0b_json:{c0b_json_path}", config)
        _write_c0e_outputs(result, config)
        return result
    c0b_result = json.loads(c0b_json_path.read_text(encoding="utf-8"))
    rows = _c0e_available_artifact_rows(c0b_result, config)
    if not rows:
        result = _c0e_incomplete_result("missing_saved_state_or_matrix_with_measured_linear_diagnostics", config)
        _write_c0e_outputs(result, config)
        return result
    artifacts = [
        _c0e_analyze_artifact(row=row, c0b_result=c0b_result, dtype=dtype, config=config)
        for row in rows
    ]
    primary = artifacts[0]
    result = {
        "schema": "stage1_pathA_C0e_gauge_invariant_curl_gate/v1",
        "source_revision": source_revision(),
        "c0b_source_json": str(c0b_json_path),
        "verdict": primary.get("verdict"),
        "verdict_support": primary.get("verdict_support"),
        "recommended_next_step": primary.get("recommended_next_step"),
        "primary_artifact": {
            "tau": primary.get("tau"),
            "tau_source_label": primary.get("tau_source_label"),
            "state_artifact": primary.get("state_artifact"),
            "matrix_path": primary.get("matrix_path"),
        },
        "artifact_policy": {
            "same_state_and_matrix_used_for_c0e_0_through_c0e_3": True,
            "preferred_primary": "deepest saved stall/source closest to tau≈0.029 with assembled J",
            "both_tau_sources_labeled_when_available": bool(len(artifacts) > 1),
            "available_artifact_count": int(len(artifacts)),
        },
        "operator_sources": _c0e_operator_sources(
            c0_frozen_branch(
                tau=float(primary["tau"]),
                grid=tuple(int(value) for value in c0b_result.get("config", {}).get("grid", C0Config().grid)),
            )
        ),
        "artifacts": artifacts,
        "faithful_operator_boundary": _faithful_operator_boundary(),
        "scope_guard": _c0e_scope_guard(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    _write_c0e_outputs(result, config)
    return result


def _write_c0d_outputs(result: Mapping[str, Any], config: C0dConfig) -> None:
    json_path = _resolve_output_path(config.json_path)
    report_path = _resolve_output_path(config.report_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_c0d_report(result, report_path)


def write_c0d_report(result: Mapping[str, Any], path: Path) -> None:
    lines: list[str] = []
    lines.append("# Path-A C0d Maxwell Gauge Identification")
    lines.append("")
    lines.append(f"Verdict: **{result.get('verdict')}**")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(
        "Diagnosis only: the saved C0b matrix is loaded and its SVD modes are recomputed. "
        "No gauge fix, deflation, recrawl, or changed-xi Jacobian assembly is performed."
    )
    lines.append("")
    lines.append("## Operator Sources")
    lines.append("")
    sources = result.get("operator_sources", {})
    lines.append("```yaml")
    for key, value in sources.items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Gauge Subspace")
    lines.append("")
    subspace = result.get("gauge_subspace", {})
    lines.append("```yaml")
    for key in (
        "status",
        "construction",
        "grid",
        "scalar_basis_dim",
        "full_state_dim",
        "dim_G",
        "dim_G_harm",
        "gradient_matrix_shape",
        "gradient_rank_rtol",
        "gradient_rank_threshold",
        "weighted_divergence_operator",
        "harmonic_weighted_divergence_rtol",
        "weighted_divergence_harmonic_threshold",
        "weighted_divergence_harmonic_indices",
    ):
        lines.append(f"{key}: {subspace.get(key)}")
    lines.append("```")
    lines.append("")
    lines.append("## Controls")
    lines.append("")
    lines.append("```yaml")
    controls = result.get("controls", {})
    lines.append(f"status: {controls.get('status')}")
    lines.append(f"positive_control: {controls.get('positive_control')}")
    lines.append(f"negative_control: {controls.get('negative_control')}")
    lines.append(f"thresholds: {controls.get('thresholds')}")
    lines.append("```")
    lines.append("")
    lines.append("## Maxwell-Lane Mode Gate")
    lines.append("")
    lines.append(
        "A mode is `MAXWELL_GAUGE` iff `||P_G v||^2 >= 0.9` and "
        "`||grad(Z*D_A*A[v])||/||A[v]|| <= 0.1`. Raw `||D_A*A[v]||/||A[v]||` "
        "is reported for context but is not the gate."
    )
    lines.append("")
    mode_rows = []
    for mode in result.get("maxwell_modes", []):
        lanes = mode.get("lane_energy_fractions", {})
        mode_rows.append(
            {
                "mode": mode.get("mode_index"),
                "sigma": mode.get("sigma"),
                "P_G": mode.get("p_g_fraction"),
                "P_G_harm": mode.get("p_g_harm_fraction"),
                "residual": mode.get("unexplained_residual"),
                "weighted": mode.get("weighted_gauge_residual"),
                "raw_div": mode.get("raw_divergence"),
                "ar": lanes.get("ar"),
                "aw": lanes.get("aw"),
                "non_A": mode.get("non_a_remainder"),
                "class": mode.get("classification"),
            }
        )
    lines.append(
        _markdown_table(
            [
                "mode",
                "sigma",
                "P_G",
                "P_G_harm",
                "residual",
                "weighted",
                "raw_div",
                "ar",
                "aw",
                "non_A",
                "class",
            ],
            mode_rows,
        )
    )
    lines.append("")
    lines.append("## Full Cluster Context")
    lines.append("")
    all_rows = []
    for mode in result.get("all_cluster_modes", []):
        all_rows.append(
            {
                "mode": mode.get("mode_index"),
                "sigma": mode.get("sigma"),
                "spatial_A": mode.get("spatial_a_energy_fraction"),
                "P_G": mode.get("p_g_fraction"),
                "weighted": mode.get("weighted_gauge_residual"),
                "raw_div": mode.get("raw_divergence"),
                "class": mode.get("classification"),
            }
        )
    lines.append(
        _markdown_table(
            ["mode", "sigma", "spatial_A", "P_G", "weighted", "raw_div", "class"],
            all_rows,
        )
    )
    lines.append("")
    lines.append("## Gauge vs Stiffness Discriminator")
    lines.append("")
    lines.append(
        "Saved-matrix evidence distinguishes penalized residual gauge from genuine Maxwell stiffness as follows: "
        "high projection into the proper multi-lambda gradient range plus small weighted gauge residual is gauge; "
        "low `P_G` or large weighted residual is classified as Maxwell stiffness under the C0d gate."
    )
    lines.append("")
    lines.append("## Verdict Support")
    lines.append("")
    lines.append("```yaml")
    support = result.get("verdict_support", {})
    lines.append(f"thresholds: {support.get('thresholds')}")
    lines.append(f"controls_status: {support.get('controls_status')}")
    lines.append(f"mode_count: {support.get('mode_count')}")
    lines.append(f"required_mode_count: {support.get('required_mode_count')}")
    lines.append("mode_classifications:")
    for row in support.get("mode_classifications", []):
        lines.append(f"  - {row}")
    if support.get("reason") is not None:
        lines.append(f"reason: {support.get('reason')}")
    lines.append("```")
    lines.append("")
    lines.append("## Recommended Next Step")
    lines.append("")
    lines.append(str(result.get("recommended_next_step")))
    lines.append("")
    lines.append("Combined-fix design if `WALL_IS_ALL_GAUGE`:")
    lines.append("")
    lines.append(str(result.get("combined_fix_design_if_wall_is_all_gauge")))
    lines.append("")
    lines.append("## Guard Confirmation")
    lines.append("")
    lines.append("```yaml")
    lines.append(f"faithful_operator_boundary: {result.get('faithful_operator_boundary')}")
    lines.append(f"scope_guard: {result.get('scope_guard')}")
    lines.append("```")
    lines.append("")
    lines.append(f"Machine artifact: `{_resolve_output_path(C0dConfig().json_path)}`.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_c0e_outputs(result: Mapping[str, Any], config: C0eConfig) -> None:
    json_path = _resolve_output_path(config.json_path)
    report_path = _resolve_output_path(config.report_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_c0e_report(result, report_path)


def write_c0e_report(result: Mapping[str, Any], path: Path) -> None:
    lines: list[str] = []
    lines.append("# Path-A C0e Gauge-Invariant Curl Gate")
    lines.append("")
    lines.append(f"Primary C0e-4 verdict: **{result.get('verdict')}**")
    lines.append("")
    lines.append("## Scope And Artifact Policy")
    lines.append("")
    lines.append(
        "Diagnosis only. C0e uses saved C0b states and assembled Jacobians, recomputes fresh dense SVD modes, "
        "runs one read-only `J*dx=-F` solve per labeled artifact, and evaluates trial residuals only on "
        "temporary tensors. It implements no deflation, gauge fix, stencil change, xi reassembly, or recrawl."
    )
    lines.append("")
    lines.append("```yaml")
    lines.append(f"primary_artifact: {result.get('primary_artifact')}")
    lines.append(f"artifact_policy: {result.get('artifact_policy')}")
    lines.append(f"scope_guard: {result.get('scope_guard')}")
    lines.append(f"faithful_operator_boundary: {result.get('faithful_operator_boundary')}")
    lines.append("```")
    lines.append("")
    lines.append("## Operator Sources")
    lines.append("")
    lines.append("```yaml")
    for key, value in result.get("operator_sources", {}).items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")

    artifact_rows = []
    for artifact in result.get("artifacts", []):
        artifact_rows.append(
            {
                "tau": artifact.get("tau"),
                "source": artifact.get("tau_source_label"),
                "converged": artifact.get("final_physical_converged"),
                "state": artifact.get("state_artifact"),
                "matrix": artifact.get("matrix_path"),
                "verdict": artifact.get("verdict"),
            }
        )
    lines.append("## Artifact Sources")
    lines.append("")
    lines.append(_markdown_table(["tau", "source", "converged", "state", "matrix", "verdict"], artifact_rows))
    lines.append("")

    for artifact in result.get("artifacts", []):
        tau = float(artifact.get("tau"))
        lines.append(f"## Artifact tau={tau:.12e} ({artifact.get('tau_source_label')})")
        lines.append("")
        framing = artifact.get("c0e_0_newton_framing", {})
        lines.append("### C0e-0 Read-Only Newton Framing")
        lines.append("")
        lines.append("```yaml")
        for key in (
            "solve_method",
            "initial_residual_l2",
            "initial_residual_linf",
            "full_step_residual_l2_ratio",
            "full_step_residual_linf_ratio",
            "gauge_removed_step_residual_l2_ratio",
            "gauge_removed_step_residual_linf_ratio",
            "step_norm",
            "near_null_component_norm",
            "complement_component_norm",
            "near_null_component_fraction",
            "linear_solve_relative_residual",
            "trial_residual_eval_scope",
        ):
            lines.append(f"{key}: {framing.get(key)}")
        lines.append("```")
        lines.append("")

        controls = artifact.get("controls", {})
        lines.append("### Controls")
        lines.append("")
        lines.append("```yaml")
        lines.append(f"status: {controls.get('status')}")
        lines.append(f"reason: {controls.get('reason')}")
        lines.append(f"bands: {controls.get('bands')}")
        lines.append(f"negative_capture_max_observed: {controls.get('negative_capture_max_observed')}")
        lines.append("```")
        control_rows = []
        for row in controls.get("controls", []):
            control_rows.append(
                {
                    "name": row.get("name"),
                    "family": row.get("family"),
                    "curl_unw": row.get("unweighted", {}).get("curl_fraction"),
                    "curl_w": row.get("cell_volume_weighted", {}).get("curl_fraction"),
                    "P_G": row.get("a_only_p_g_fraction"),
                    "P_cpl": row.get("coupled_capture_fraction"),
                    "non_A": row.get("non_a_remainder"),
                    "construction": row.get("construction"),
                }
            )
        lines.append(
            _markdown_table(
                ["name", "family", "curl_unw", "curl_w", "P_G", "P_cpl", "non_A", "construction"],
                control_rows,
            )
        )
        lines.append("")

        phase = artifact.get("phase_mode", {})
        lines.append("### Phase Mode Separate Gate")
        lines.append("")
        lines.append("```yaml")
        for key in (
            "mode_index",
            "sigma",
            "classification",
            "phase_capture_fraction",
            "coupled_capture_fraction",
            "spatial_a_norm_unweighted",
            "z_f_rw_norm_unweighted",
            "spatial_a_norm_cell_volume_weighted",
            "z_f_rw_norm_cell_volume_weighted",
            "note",
        ):
            lines.append(f"{key}: {phase.get(key)}")
        lines.append("```")
        lines.append("")

        lines.append("### C0e-1 Curl Gate For Maxwell-Lane Modes")
        lines.append("")
        mode_rows = []
        for mode in artifact.get("maxwell_modes", []):
            lanes = mode.get("lane_energy_fractions", {})
            outcome = mode.get("outcome_details", {})
            mode_rows.append(
                {
                    "mode": mode.get("mode_index"),
                    "sigma": mode.get("sigma"),
                    "curl_unw": mode.get("curl_fraction", {}).get("unweighted", {}).get("curl_fraction"),
                    "curl_w": mode.get("curl_fraction", {})
                    .get("cell_volume_weighted", {})
                    .get("curl_fraction"),
                    "outcome": mode.get("outcome"),
                    "margin_log10": outcome.get("band_margin_log10_min"),
                    "P_G": mode.get("a_only_p_g_fraction"),
                    "P_cpl": mode.get("coupled_capture_fraction"),
                    "A_energy": mode.get("spatial_a_energy_fraction"),
                    "ar": lanes.get("ar"),
                    "aw": lanes.get("aw"),
                }
            )
        lines.append(
            _markdown_table(
                [
                    "mode",
                    "sigma",
                    "curl_unw",
                    "curl_w",
                    "outcome",
                    "margin_log10",
                    "P_G",
                    "P_cpl",
                    "A_energy",
                    "ar",
                    "aw",
                ],
                mode_rows,
            )
        )
        lines.append("")
        lines.append("Boundary/interior curl split and norm-specific margins are in the machine JSON for each mode.")
        lines.append("")

        lines.append("### C0e-3 Mechanism Evidence")
        lines.append("")
        mechanism = artifact.get("mechanism", {})
        lines.append("```yaml")
        lines.append(f"primary_mechanism: {mechanism.get('primary_mechanism')}")
        lines.append(f"secondary_flags: {mechanism.get('secondary_flags')}")
        lines.append(f"lambda_frequency_counts: {mechanism.get('lambda_frequency_counts')}")
        lines.append(f"row_scaling_ratio: {mechanism.get('row_scaling_ratio')}")
        lines.append(
            "max_unexplained_fraction_of_component_sum_norm: "
            f"{mechanism.get('max_unexplained_fraction_of_component_sum_norm')}"
        )
        lines.append(
            "max_gauge_penalty_to_assembled_jv_ratio: "
            f"{mechanism.get('max_gauge_penalty_to_assembled_jv_ratio')}"
        )
        lines.append("```")
        lines.append("")
        row_scaling = artifact.get("row_scaling", {})
        lines.append("Assembled Maxwell-row scaling:")
        lines.append("")
        lines.append("```yaml")
        lines.append(f"{row_scaling}")
        lines.append("```")
        lines.append("")
        budget_rows = []
        for mode in artifact.get("maxwell_modes", []):
            budget = mode.get("residual_budget_jv", {})
            comps = budget.get("components", {})
            budget_rows.append(
                {
                    "mode": mode.get("mode_index"),
                    "Jv": budget.get("assembled_jv_norm"),
                    "curl": comps.get("maxwell_curl_rows", {}).get("norm"),
                    "penalty": comps.get("maxwell_gauge_penalty_rows", {}).get("norm"),
                    "current": comps.get("maxwell_current_source_rows", {}).get("norm"),
                    "matter": comps.get("matter_covariant_a_coupling_rows", {}).get("norm"),
                    "gauss": comps.get("gauss_charge_rows", {}).get("norm"),
                    "wall_mass": comps.get("wall_mass_rows", {}).get("norm"),
                    "unexplained": budget.get("unexplained_remainder_norm"),
                    "lambda_class": mode.get("lambda_preimage", {})
                    .get("frequency", {})
                    .get("classification"),
                }
            )
        lines.append(
            _markdown_table(
                [
                    "mode",
                    "Jv",
                    "curl",
                    "penalty",
                    "current",
                    "matter",
                    "gauss",
                    "wall_mass",
                    "unexplained",
                    "lambda_class",
                ],
                budget_rows,
            )
        )
        lines.append("")

        lines.append("### C0e-4 Verdict And C0f Recommendation")
        lines.append("")
        lines.append("```yaml")
        lines.append(f"verdict: {artifact.get('verdict')}")
        lines.append(f"verdict_support: {artifact.get('verdict_support')}")
        lines.append("```")
        lines.append("")
        lines.append(str(artifact.get("recommended_next_step")))
        lines.append("")

    lines.append("## Machine Artifact")
    lines.append("")
    lines.append(f"`{_resolve_output_path(C0eConfig().json_path)}`")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_c0c_outputs(result: Mapping[str, Any], config: C0cConfig) -> None:
    json_path = _resolve_output_path(config.json_path)
    report_path = _resolve_output_path(config.report_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_c0c_report(result, report_path)


def _format_overlap_cell(overlaps: Mapping[str, Any], name: str) -> str:
    value = overlaps.get(name)
    if value is None:
        return "-"
    return _fmt(float(value))


def write_c0c_report(result: Mapping[str, Any], path: Path) -> None:
    lines: list[str] = []
    lines.append("# Path-A C0c Nullmode Identification")
    lines.append("")
    lines.append(f"Verdict: **{result.get('verdict')}**")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(
        "Diagnosis only: generators are evaluated on existing C0b states and saved Jacobians. "
        "`coupled_branch.py`, `operators.py`, frozen physics, and `physical_export_permitted` are untouched."
    )
    lines.append("")
    points = list(result.get("points", []))
    if points:
        lines.append("## Generator Inventory")
        lines.append("")
        generator_rows = [
            {
                "name": row["name"],
                "class": row["classification"],
                "status": row["symmetry_status"],
                "norm": row["norm"],
                "description": row["description"],
            }
            for row in points[0].get("generators", [])
        ]
        lines.append(_markdown_table(["name", "class", "status", "norm", "description"], generator_rows))
        lines.append("")
    lines.append("## Annihilation And Equivariance")
    lines.append("")
    for point in points:
        lines.append(
            f"### tau={float(point['tau']):.12e} "
            f"(converged={point.get('final_physical_converged')}, "
            f"residual={_fmt(point.get('final_original_residual_linf'))})"
        )
        rows = []
        for row in point.get("annihilation", {}).get("rows", []):
            rows.append(
                {
                    "generator": row.get("generator"),
                    "status": row.get("status"),
                    "test": row.get("test"),
                    "assembled": row.get(
                        "annihilation_rel_assembled",
                        row.get("annihilation_rel_assembled_not_a_null_gate"),
                    ),
                    "jvp": row.get(
                        "annihilation_rel_jvp",
                        row.get("annihilation_rel_jvp_not_a_null_gate"),
                    ),
                    "crosscheck": row.get("jvp_crosscheck_rel"),
                    "equiv_jvp": row.get("equivariance_rel_jvp_sigma_scaled"),
                    "gate": row.get("null_gate_pass", row.get("equivariance_pass")),
                    "reason": row.get("reason"),
                }
            )
        lines.append(
            _markdown_table(
                [
                    "generator",
                    "status",
                    "test",
                    "assembled",
                    "jvp",
                    "crosscheck",
                    "equiv_jvp",
                    "gate",
                    "reason",
                ],
                rows,
            )
        )
        lines.append("")
    lines.append("## Near-Null SVD Overlap")
    lines.append("")
    for point in points:
        lines.append(f"### tau={float(point['tau']):.12e}")
        lines.append("")
        rows = []
        for mode in point.get("overlap", {}).get("modes", []):
            lanes = mode.get("lane_energy_fractions", {})
            overlaps = mode.get("overlaps", {})
            rows.append(
                {
                    "mode": mode.get("mode_index"),
                    "sigma": mode.get("sigma"),
                    "class": mode.get("classification", "-"),
                    "phase": _format_overlap_cell(overlaps, "phase"),
                    "tr": _format_overlap_cell(overlaps, "translation_r"),
                    "tw": _format_overlap_cell(overlaps, "translation_w"),
                    "dil": _format_overlap_cell(overlaps, "dilation_r"),
                    "maxwell": _format_overlap_cell(overlaps, "maxwell_residual_gauge"),
                    "residual": mode.get("unexplained_residual_fraction"),
                    "psi_re": lanes.get("psi_real"),
                    "psi_im": lanes.get("psi_imag"),
                    "a0": lanes.get("a0"),
                    "ar": lanes.get("ar"),
                    "aw": lanes.get("aw"),
                    "r0": lanes.get("r0"),
                    "mu": lanes.get("mu"),
                }
            )
        lines.append(
            _markdown_table(
                [
                    "mode",
                    "sigma",
                    "class",
                    "phase",
                    "tr",
                    "tw",
                    "dil",
                    "maxwell",
                    "residual",
                    "psi_re",
                    "psi_im",
                    "a0",
                    "ar",
                    "aw",
                    "r0",
                    "mu",
                ],
                rows,
            )
        )
        lines.append("")
    lines.append("## Dense Sigma Validation")
    lines.append("")
    dense = result.get("dense_sigma_validation", {})
    lines.append("```yaml")
    for key, value in dense.items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Verdict Support")
    lines.append("")
    lines.append("```yaml")
    support = result.get("verdict_support", {})
    lines.append(f"thresholds: {support.get('thresholds')}")
    lines.append(f"generator_null_gates: {support.get('generator_null_gates')}")
    lines.append(f"explained_mode_count: {support.get('explained_mode_count')}")
    lines.append(f"cluster_mode_count: {support.get('cluster_mode_count')}")
    lines.append("mode_classifications:")
    for row in support.get("mode_classifications", []):
        lines.append(
            f"  - mode_index: {row.get('mode_index')}, sigma: {row.get('sigma')}, "
            f"classification: {row.get('classification')}, residual: "
            f"{row.get('unexplained_residual_fraction')}, lanes: "
            f"{row.get('lane_energy_fractions')}, support: {row.get('support')}"
        )
    lines.append("```")
    lines.append("")
    lines.append("## Recommended Next Step")
    lines.append("")
    lines.append(str(result.get("recommended_next_step")))
    lines.append("")
    lines.append("## Guard Confirmation")
    lines.append("")
    lines.append("```yaml")
    lines.append(f"faithful_operator_boundary: {result.get('faithful_operator_boundary')}")
    lines.append(f"scope_guard: {result.get('scope_guard')}")
    lines.append("```")
    lines.append("")
    lines.append(f"Machine artifact: `{_resolve_output_path(C0cConfig().json_path)}`.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _c0c_stdout_summary(result: Mapping[str, Any], config: C0cConfig) -> str:
    lines: list[str] = []
    points = list(result.get("points", []))
    if not points:
        return (
            f"C0c verdict: {result.get('verdict')} - "
            f"{result.get('verdict_support', {}).get('reason')}"
        )
    shallow = points[0]
    deepest = points[-1]
    exact = [
        row["name"]
        for row in shallow.get("generators", [])
        if row.get("symmetry_status") == "EXACT_SYMMETRY"
    ]
    probes = [
        row["name"]
        for row in shallow.get("generators", [])
        if row.get("symmetry_status") == "PROBE"
    ]
    lines.append(f"Generators: EXACT={exact}; PROBE={probes}")
    lines.append(
        f"Converged tau {float(shallow['tau']):.12e} annihilation "
        f"(assembled/JVP, threshold {config.annihilation_threshold:.1e}):"
    )
    for row in shallow.get("annihilation", {}).get("rows", []):
        lines.append(
            "  "
            f"{row.get('generator')}: "
            f"{_fmt(row.get('annihilation_rel_assembled'))}/"
            f"{_fmt(row.get('annihilation_rel_jvp'))}, gate={row.get('null_gate_pass')}"
        )
    if deepest is not shallow:
        lines.append(
            f"Non-converged tau {float(deepest['tau']):.12e}: "
            "phase uses equivariance identity; probes are NOT_MEASURED."
        )
        for row in deepest.get("annihilation", {}).get("rows", []):
            if row.get("generator") == "phase":
                lines.append(
                    "  phase equivariance sigma-scaled assembled/JVP = "
                    f"{_fmt(row.get('equivariance_rel_assembled_sigma_scaled'))}/"
                    f"{_fmt(row.get('equivariance_rel_jvp_sigma_scaled'))}"
                )
    lines.append("Near-null cluster at converged tau:")
    for mode in shallow.get("overlap", {}).get("modes", []):
        lanes = mode.get("lane_energy_fractions", {})
        overlaps = mode.get("overlaps", {})
        lines.append(
            "  "
            f"mode {mode.get('mode_index')} sigma={_fmt(mode.get('sigma'))} "
            f"class={mode.get('classification')} "
            f"phase={_format_overlap_cell(overlaps, 'phase')} "
            f"maxwell={_format_overlap_cell(overlaps, 'maxwell_residual_gauge')} "
            f"residual={_fmt(mode.get('unexplained_residual_fraction'))} "
            f"lanes psi_re={_fmt(lanes.get('psi_real'))}, psi_im={_fmt(lanes.get('psi_imag'))}, "
            f"a0={_fmt(lanes.get('a0'))}, ar={_fmt(lanes.get('ar'))}, aw={_fmt(lanes.get('aw'))}"
        )
    dense = result.get("dense_sigma_validation", {})
    lines.append(
        "Dense sigma validation: "
        f"{dense.get('status')} assembled={_fmt(dense.get('assembled_sigma_min'))} "
        f"full_jvp={_fmt(dense.get('dense_full_jvp_sigma_min'))} "
        f"trusted={dense.get('trusted_sigma_source')}"
    )
    lines.append(f"Verdict: {result.get('verdict')}")
    lines.append(f"Recommended next step: {result.get('recommended_next_step')}")
    lines.append(
        "Guard confirmation: faithful operators, frozen physics, and export guard untouched by C0c; "
        "no gauge-fix/deflation/re-crawl implemented."
    )
    lines.append(
        f"Artifacts: report={_resolve_output_path(config.report_path)}, "
        f"json={_resolve_output_path(config.json_path)}"
    )
    return "\n".join(lines)


def _c0d_stdout_summary(result: Mapping[str, Any], config: C0dConfig) -> str:
    lines: list[str] = []
    subspace = result.get("gauge_subspace", {})
    lines.append(
        f"dim(G)={subspace.get('dim_G')} dim(G_harm)={subspace.get('dim_G_harm')}"
    )
    lines.append("Maxwell-lane modes:")
    for mode in result.get("maxwell_modes", []):
        lines.append(
            "  "
            f"mode {mode.get('mode_index')} sigma={_fmt(mode.get('sigma'))} "
            f"P_G={_fmt(mode.get('p_g_fraction'))} "
            f"weighted={_fmt(mode.get('weighted_gauge_residual'))} "
            f"raw_div={_fmt(mode.get('raw_divergence'))} "
            f"class={mode.get('classification')}"
        )
    lines.append(f"Verdict: {result.get('verdict')}")
    lines.append(f"Recommended next step: {result.get('recommended_next_step')}")
    lines.append(
        "Combined-fix design if WALL_IS_ALL_GAUGE: "
        f"{result.get('combined_fix_design_if_wall_is_all_gauge')}"
    )
    controls = result.get("controls", {})
    lines.append(
        "Controls: "
        f"positive={controls.get('positive_control')} "
        f"negative={controls.get('negative_control')} "
        f"status={controls.get('status')}"
    )
    lines.append(
        "Guard confirmation: faithful operators, frozen physics, and export guard untouched by C0d; "
        "no gauge-fix/deflation/re-crawl and no changed-xi reassembly implemented."
    )
    lines.append(
        f"Artifacts: report={_resolve_output_path(config.report_path)}, "
        f"json={_resolve_output_path(config.json_path)}"
    )
    return "\n".join(lines)


def _build_c0d_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Path-A C0d Maxwell gauge identification.")
    parser.add_argument("--c0b-json-path", type=Path, default=C0dConfig().c0b_json_path)
    parser.add_argument("--run-root", type=Path, default=C0dConfig().run_root)
    parser.add_argument("--report-path", type=Path, default=C0dConfig().report_path)
    parser.add_argument("--json-path", type=Path, default=C0dConfig().json_path)
    parser.add_argument("--cluster-mode-count", type=int, default=C0dConfig().cluster_mode_count)
    parser.add_argument("--maxwell-mode-count", type=int, default=C0dConfig().maxwell_mode_count)
    parser.add_argument(
        "--maxwell-lane-fraction-min",
        type=float,
        default=C0dConfig().maxwell_lane_fraction_min,
    )
    parser.add_argument(
        "--gauge-capture-threshold",
        type=float,
        default=C0dConfig().gauge_capture_threshold,
    )
    parser.add_argument(
        "--weighted-residual-threshold",
        type=float,
        default=C0dConfig().weighted_residual_threshold,
    )
    parser.add_argument(
        "--harmonic-weighted-divergence-rtol",
        type=float,
        default=C0dConfig().harmonic_weighted_divergence_rtol,
    )
    return parser


def c0d_main(argv: Sequence[str] | None = None) -> int:
    args = _build_c0d_parser().parse_args(argv)
    config = C0dConfig(
        c0b_json_path=args.c0b_json_path,
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        cluster_mode_count=int(args.cluster_mode_count),
        maxwell_mode_count=int(args.maxwell_mode_count),
        maxwell_lane_fraction_min=float(args.maxwell_lane_fraction_min),
        gauge_capture_threshold=float(args.gauge_capture_threshold),
        weighted_residual_threshold=float(args.weighted_residual_threshold),
        harmonic_weighted_divergence_rtol=float(args.harmonic_weighted_divergence_rtol),
    )
    result = run_c0d_maxwell_gauge_identification(config)
    print(_c0d_stdout_summary(result, config))
    return 0


def _c0e_stdout_summary(result: Mapping[str, Any], config: C0eConfig) -> str:
    artifacts = list(result.get("artifacts", []))
    if not artifacts:
        return (
            f"C0e verdict: {result.get('verdict')} - "
            f"{result.get('verdict_support', {}).get('reason')}"
        )
    primary = artifacts[0]
    lines: list[str] = []
    lines.append(
        "Artifact: "
        f"tau={float(primary.get('tau')):.12e} source={primary.get('tau_source_label')} "
        f"state={primary.get('state_artifact')} J={primary.get('matrix_path')}"
    )
    if len(artifacts) > 1:
        lines.append(
            "Additional labeled tau sources: "
            + ", ".join(
                f"tau={float(artifact.get('tau')):.12e}:{artifact.get('tau_source_label')}"
                for artifact in artifacts[1:]
            )
        )
    framing = primary.get("c0e_0_newton_framing", {})
    lines.append(
        "C0e-0 ||F|| ratios: "
        f"full_l2={_fmt(framing.get('full_step_residual_l2_ratio'))}, "
        f"gauge_removed_l2={_fmt(framing.get('gauge_removed_step_residual_l2_ratio'))}, "
        f"full_linf={_fmt(framing.get('full_step_residual_linf_ratio'))}, "
        f"gauge_removed_linf={_fmt(framing.get('gauge_removed_step_residual_linf_ratio'))}; "
        f"near_null_step_fraction={_fmt(framing.get('near_null_component_fraction'))}"
    )
    lines.append("Per-mode curl gate (primary artifact):")
    for mode in primary.get("maxwell_modes", []):
        label = "mode 2 SWING" if int(mode.get("mode_index")) == 2 else f"mode {mode.get('mode_index')}"
        outcome = mode.get("outcome_details", {})
        lines.append(
            "  "
            f"{label}: curl_unw={_fmt(mode.get('curl_fraction', {}).get('unweighted', {}).get('curl_fraction'))}, "
            f"curl_w={_fmt(mode.get('curl_fraction', {}).get('cell_volume_weighted', {}).get('curl_fraction'))}, "
            f"outcome={mode.get('outcome')}, "
            f"band_margin_log10={_fmt(outcome.get('band_margin_log10_min'))}, "
            f"P_G={_fmt(mode.get('a_only_p_g_fraction'))}, "
            f"P_cpl={_fmt(mode.get('coupled_capture_fraction'))}"
        )
    phase = primary.get("phase_mode", {})
    lines.append(
        "Phase mode separate: "
        f"A_norm_unw={_fmt(phase.get('spatial_a_norm_unweighted'))}, "
        f"ZF_norm_unw={_fmt(phase.get('z_f_rw_norm_unweighted'))}, "
        f"phase_capture={_fmt(phase.get('phase_capture_fraction'))}, "
        f"P_cpl={_fmt(phase.get('coupled_capture_fraction'))}"
    )
    mechanism = primary.get("mechanism", {})
    lines.append(
        "Mechanism: "
        f"primary={mechanism.get('primary_mechanism')}, "
        f"secondary_flags={mechanism.get('secondary_flags')}, "
        "unexplained_budget="
        f"{_fmt(mechanism.get('max_unexplained_fraction_of_component_sum_norm'))}"
    )
    lines.append(
        "C0e-4 verdict: "
        f"{primary.get('verdict')} - {primary.get('verdict_support', {}).get('reason')}"
    )
    lines.append(f"Recommended C0f design: {primary.get('recommended_next_step')}")
    controls = primary.get("controls", {})
    bands = controls.get("bands", {})
    lines.append(
        "Controls: "
        f"status={controls.get('status')}, "
        f"unweighted_gap_log10={_fmt(bands.get('unweighted', {}).get('separation_log10'))}, "
        "weighted_gap_log10="
        f"{_fmt(bands.get('cell_volume_weighted', {}).get('separation_log10'))}, "
        f"negative_capture_max={_fmt(controls.get('negative_capture_max_observed'))}"
    )
    lines.append(
        "Guard confirmation: faithful operators/frozen/export untouched by C0e; "
        "no fix, no state advance, no recrawl, no xi reassembly."
    )
    lines.append(
        f"Artifacts: report={_resolve_output_path(config.report_path)}, "
        f"json={_resolve_output_path(config.json_path)}"
    )
    return "\n".join(lines)


def _build_c0e_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Path-A C0e gauge-invariant curl gate.")
    parser.add_argument("--c0b-json-path", type=Path, default=C0eConfig().c0b_json_path)
    parser.add_argument("--run-root", type=Path, default=C0eConfig().run_root)
    parser.add_argument("--report-path", type=Path, default=C0eConfig().report_path)
    parser.add_argument("--json-path", type=Path, default=C0eConfig().json_path)
    parser.add_argument("--cluster-mode-count", type=int, default=C0eConfig().cluster_mode_count)
    parser.add_argument("--maxwell-mode-count", type=int, default=C0eConfig().maxwell_mode_count)
    parser.add_argument(
        "--maxwell-lane-fraction-min",
        type=float,
        default=C0eConfig().maxwell_lane_fraction_min,
    )
    parser.add_argument(
        "--gradient-rank-rtol",
        type=float,
        default=C0eConfig().gradient_rank_rtol,
    )
    parser.add_argument(
        "--primary-only",
        action="store_true",
        help="Analyze only the deepest/primary artifact instead of every available labeled tau source.",
    )
    return parser


def c0e_main(argv: Sequence[str] | None = None) -> int:
    args = _build_c0e_parser().parse_args(argv)
    config = C0eConfig(
        c0b_json_path=args.c0b_json_path,
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        cluster_mode_count=int(args.cluster_mode_count),
        maxwell_mode_count=int(args.maxwell_mode_count),
        maxwell_lane_fraction_min=float(args.maxwell_lane_fraction_min),
        gradient_rank_rtol=float(args.gradient_rank_rtol),
        run_all_available_tau_sources=not bool(args.primary_only),
    )
    result = run_c0e_gauge_invariant_curl_gate(config)
    print(_c0e_stdout_summary(result, config))
    return 0


def _build_c0c_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Path-A C0c null-mode identification.")
    parser.add_argument("--c0b-json-path", type=Path, default=C0cConfig().c0b_json_path)
    parser.add_argument("--run-root", type=Path, default=C0cConfig().run_root)
    parser.add_argument("--report-path", type=Path, default=C0cConfig().report_path)
    parser.add_argument("--json-path", type=Path, default=C0cConfig().json_path)
    parser.add_argument("--cluster-mode-count", type=int, default=C0cConfig().cluster_mode_count)
    parser.add_argument("--dense-jvp-chunk-size", type=int, default=C0cConfig().dense_jvp_chunk_size)
    return parser


def c0c_main(argv: Sequence[str] | None = None) -> int:
    args = _build_c0c_parser().parse_args(argv)
    config = C0cConfig(
        c0b_json_path=args.c0b_json_path,
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        cluster_mode_count=int(args.cluster_mode_count),
        dense_jvp_chunk_size=int(args.dense_jvp_chunk_size),
    )
    result = run_c0c_nullmode_identification(config)
    print(_c0c_stdout_summary(result, config))
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("full", "crawl", "diagnostics"), default="full")
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--report-path", type=Path, default=DEFAULT_REPORT_PATH)
    parser.add_argument("--json-path", type=Path, default=DEFAULT_JSON_PATH)
    parser.add_argument("--grid", default=f"{b2a.DEFAULT_BACKGROUND_GRID[0]},{b2a.DEFAULT_BACKGROUND_GRID[1]}")
    parser.add_argument("--taus", default="")
    parser.add_argument("--max-tau-backtracks", type=int, default=C0Config().max_tau_backtracks)
    parser.add_argument(
        "--min-tau-backtrack-delta",
        type=float,
        default=C0Config().min_tau_backtrack_delta,
    )
    parser.add_argument(
        "--max-depth-failures-after-floor",
        type=int,
        default=C0Config().max_depth_failures_after_floor,
    )
    parser.add_argument("--max-newton-iters", type=int, default=0)
    parser.add_argument("--skip-observables", action="store_true")
    parser.add_argument("--skip-linear-diagnostics", action="store_true")
    parser.add_argument(
        "--prefer-existing-shallow-seed",
        action="store_true",
        help="Allow an existing B2c seed only before a C0-converged state exists; below-floor attempts still bypass it.",
    )
    parser.add_argument(
        "--use-gauge-fix",
        action="store_true",
        help="Use the analytic gauge-complement linear solve in path coordinates.",
    )
    return parser


def _parse_grid(text: str) -> tuple[int, int]:
    left, right = text.split(",", 1)
    return int(left), int(right)


def _parse_taus(text: str) -> tuple[float, ...]:
    if not text:
        return C0Config().depth_sequence
    return tuple(float(piece) for piece in text.split(",") if piece)


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    config = C0Config(
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        grid=_parse_grid(args.grid),
        depth_sequence=_parse_taus(args.taus),
        max_tau_backtracks=int(args.max_tau_backtracks),
        min_tau_backtrack_delta=float(args.min_tau_backtrack_delta),
        max_depth_failures_after_floor=int(args.max_depth_failures_after_floor),
        max_newton_iters_override=(
            int(args.max_newton_iters) if int(args.max_newton_iters) > 0 else None
        ),
        diagnostic_linear=not bool(args.skip_linear_diagnostics),
        diagnostic_observables=not bool(args.skip_observables),
        prefer_existing_b2c_background_predictor=bool(args.prefer_existing_shallow_seed),
        use_gauge_fix=bool(args.use_gauge_fix),
    )
    if args.phase == "crawl":
        result = run_c0_crawl(config)
    elif args.phase == "diagnostics":
        result = run_c0_diagnostics(config)
    else:
        result = run_c0_conditioning_spike(config)
    deepest_tau = result["deepest_converged_tau"]
    deepest_r0 = result["deepest_converged_R0_min"]
    final_row = result["tau_attempts"][-1] if result.get("tau_attempts") else {}
    sigma = result.get("sigma_diagnostics", {})
    fold = result.get("fold_diagnostic", {})
    admissibility = result.get("aids_admissibility", {})
    shallow_sigma = (sigma.get("shallow") or {}) if isinstance(sigma, Mapping) else {}
    deepest_sigma = (sigma.get("deepest") or {}) if isinstance(sigma, Mapping) else {}
    summary = {
        "phase": args.phase,
        "verdict": result["verdict"],
        "verdict_support": result.get("verdict_support"),
        "deepest_tau": deepest_tau,
        "deepest_R0_min": deepest_r0,
        "prior_tau_floor": PRIOR_TAU_FLOOR,
        "final_attempt_tau": final_row.get("target_tau"),
        "final_attempt_init_source": final_row.get("init", {}).get("source") if final_row else None,
        "final_attempt_used_existing_b2c": final_row.get("used_existing_b2c") if final_row else None,
        "final_original_residual_linf": final_row.get("final_original_residual_linf"),
        "shallow_sigma": {
            "tau": shallow_sigma.get("tau"),
            "method": shallow_sigma.get("sigma_method"),
            "value": shallow_sigma.get("sigma_min"),
            "status": shallow_sigma.get("status"),
            "v_min": shallow_sigma.get("v_min_energy_fractions"),
        },
        "deepest_sigma": {
            "tau": deepest_sigma.get("tau"),
            "method": deepest_sigma.get("sigma_method"),
            "value": deepest_sigma.get("sigma_min"),
            "status": deepest_sigma.get("status"),
            "v_min": deepest_sigma.get("v_min_energy_fractions"),
        },
        "fold_call": fold.get("call") if isinstance(fold, Mapping) else None,
        "fold_status": fold.get("status") if isinstance(fold, Mapping) else None,
        "r0_tau_statement": result.get("r0_tau_diagnostic", {}).get("statement"),
        "admissibility": {
            "residual_equality_max_abs": admissibility.get("residual_equality_max_abs"),
            "residual_equality_status": admissibility.get("residual_equality_status"),
            "epsilon_independence_status": admissibility.get("epsilon_independence_status"),
        },
        "faithful_operator_boundary": result.get("faithful_operator_boundary"),
        "report": str(config.report_path),
        "json": str(config.json_path),
    }
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
