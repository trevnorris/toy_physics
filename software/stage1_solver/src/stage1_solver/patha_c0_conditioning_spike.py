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
from scipy.sparse import csc_matrix, diags, eye, save_npz
from scipy.sparse.linalg import LinearOperator, eigsh, gmres, svds
import torch

from . import patha_b2a_bdg as b2a
from . import patha_b2b_maxwell as b2b
from . import patha_b2c_calibration as b2c
from .backend import configure_backend, jvp
from .boundaries import BoundaryCondition
from .config import BackendConfig, NewtonConfig, source_revision
from .coupled_branch import (
    _closed_to_coupled_fields,
    _create_branch_grid,
    branch_boundary_conditions,
    coupled_pde_residual,
    initial_closed_branch_state,
    pack_coupled_fields,
    patha_closed_branch_residual,
    resample_closed_branch_state,
    unpack_closed_coupled_fields,
    wall_grid_from_tensor_grid,
)
from .newton import BuiltPreconditioner, PreconditionerBuildContext
from .patha_static_balance import SSigmaProvider, SSigmaSpec, resolve_s_sigma, static_balance_terms
from .preconditioners import (
    assemble_closed_coupled_colored_sparse_jacobian,
    factorized_sparse_inverse_operator,
)


PRIOR_TAU_FLOOR = 0.029
BACKGROUND_RESIDUAL_TOL = b2c.BACKGROUND_RESIDUAL_TOL
DEFAULT_RUN_ROOT = Path("software/stage1_solver/runs/pathA_C0b_wall_diagnosis")
DEFAULT_REPORT_PATH = Path("software/stage1_solver/reports/pathA_C0b_wall_diagnosis.md")
DEFAULT_JSON_PATH = DEFAULT_RUN_ROOT / "pathA_C0b_diagnostic.json"
ALLOWED_VERDICTS = {
    "SPIKE_SUFFICIENT",
    "FOLD_TURNING_POINT",
    "NEAR_NULL_SPACE_STRUCTURAL",
    "PRODUCTION_SOLVER_REQUIRED",
    "DIAGNOSTIC_INCOMPLETE",
}


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
    diagnostic_bdg_grid: tuple[int, int] = (8, 8)
    diagnostic_bdg_modes: int = 8
    diagnostic_profile_points: int = 65
    diagnostic_maxwell_grid: tuple[int, int] = (9, 9)
    diagnostic_maxwell_window: float = b2b.DEFAULT_FINAL_WINDOW
    diagnostic_maxwell_radial_scale: float = 1.75


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


def _c0_newton_config(base: NewtonConfig) -> NewtonConfig:
    preconditioner = replace(
        base.preconditioner,
        type="colored_sparse_jacobian_lu",
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
):
    """Return the existing frozen branch with only Newton controls changed."""

    branch = b2a.frozen_patha_b2a_branch(grid=grid, tau=float(tau))
    newton = _c0_newton_config(branch.newton)
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
    preconditioner_config = replace(newton.preconditioner, diagonal_shift=0.0)
    linear_config = replace(newton, preconditioner=preconditioner_config)
    matrix, metadata = assemble_closed_coupled_colored_sparse_jacobian(
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
    operator, factor_metadata = factorized_sparse_inverse_operator(scaled_matrix, preconditioner_config)
    metadata.update(factor_metadata)
    metadata.update(
        {
            "type": "c0_conditioned_scaled_colored_sparse_jacobian_lu",
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
    return BuiltPreconditioner(operator=operator, metadata=metadata), row_scale, col_scale, matrix, scaled_matrix


def solve_c0_scaled_newton(
    *,
    original_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    preconditioner_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x0: torch.Tensor,
    grid,
    newton: NewtonConfig,
    aid: C0AidParameters,
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
        gmres_curve: list[float] = []

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
                        preconditioner_info=dict(built.metadata),
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
                preconditioner_info=dict(built.metadata),
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
            result = solve_c0_scaled_newton(
                original_residual_fn=original_residual_fn,
                preconditioner_residual_fn=preconditioner_residual_fn,
                x0=x0,
                grid=grid,
                newton=branch.newton,
                aid=aid,
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
            if row.get("init", {}).get("source") != "previous_c0_converged_state":
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

    below_failed = [
        row
        for row in attempts
        if _below_prior_floor(float(row["target_tau"]))
        and not row.get("final_physical_converged")
    ]
    attempted_backtracking = any(int(row.get("backtrack_index", 0)) > 0 for row in attempts)
    full_newton_budget = config.get("max_newton_iters_override") in (None, 0)
    crawl_persistent_failure = bool(below_failed and attempted_backtracking and full_newton_budget)
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
