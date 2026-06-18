"""Path-A chunk 1c validation for the closed self-consistent solve."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from .backend import configure_backend
from .config import BackendConfig, BranchSmokeConfig, HarnessConfig, WallGridSpec, stable_hash
from .conservation_diagnostics import (
    ConservationDiagnosticsConfig,
    evaluate_conservation_state,
)
from .convergence import (
    _volume_restrict_2d,
    observed_order_from_three,
    richardson_estimate,
    validate_refinement_ladder,
)
from .coupled_branch import (
    ClosedCoupledFields,
    CoupledFields,
    _create_branch_grid,
    _max_rss_mb,
    branch_boundary_conditions,
    confinement_wall_to_matter_coefficient_torch,
    initial_closed_branch_state,
    pack_coupled_fields,
    patha_closed_branch_residual,
    patha_closed_wall_terms,
    resample_closed_branch_state,
    run_branch_continuation,
    unpack_closed_coupled_fields,
)
from .grid import TensorProductGrid, WallGrid
from .newton import solve_newton_jvp
from .patha_closed_newton import (
    default_closed_branch_config,
    default_closed_s_sigma_spec,
    target_token_scan,
)
from .patha_static_balance import (
    SSigmaProvider,
    SSigmaSpec,
    resolve_s_sigma,
    wall_center_gradient,
)
from .preconditioners import make_closed_coupled_colored_sparse_jacobian_lu_factory


def _stage_root() -> Path:
    return Path(__file__).resolve().parents[2]


def validation_scan_paths(report_path: Path | None = None) -> list[Path]:
    root = _stage_root()
    paths = [
        root / "src" / "stage1_solver" / "patha_closed_validation.py",
        root / "src" / "stage1_solver" / "patha_closed_newton.py",
        root / "src" / "stage1_solver" / "coupled_branch.py",
        root / "src" / "stage1_solver" / "patha_static_balance.py",
        root / "tests" / "test_patha_closed_validation.py",
    ]
    if report_path is not None:
        paths.append(report_path)
    return paths


def default_validation_branch_config() -> BranchSmokeConfig:
    branch = default_closed_branch_config()
    return replace(
        branch,
        newton=replace(
            branch.newton,
            preconditioner=replace(
                branch.newton.preconditioner,
                rebuild_policy="once_per_continuation",
            ),
        ),
    )


@dataclass(frozen=True)
class PathAClosedValidationConfig:
    name: str = "patha_chunk1c_self_consistent_validation"
    branch: BranchSmokeConfig = field(default_factory=default_validation_branch_config)
    s_sigma_spec: SSigmaSpec | None = None
    run_root: str = "software/stage1_solver/runs/patha_chunk1c_self_consistent_validation"
    report_path: str = "software/stage1_solver/reports/patha_chunk1c_self_consistent_validation.md"
    convergence_levels: tuple[tuple[int, int], ...] = ((4, 4), (8, 8), (16, 16))
    refinement_ratio: int = 2
    min_levels: int = 3
    expected_order: float = 2.0
    max_level_seconds: float = 590.0
    nontrivial_term_floor: float = 1.0e-8
    independent_source_abs_tol: float = 1.0e-12
    independent_flux_abs_tol: float = 1.0e-12
    independent_flux_h2_tolerance_factor: float = 8.0
    raw_difference_floor: float = 1.0e-10
    conservation_baseline_factor: float = 4.0
    conservation_floor: float = 1.0e-10
    conservation: ConservationDiagnosticsConfig = field(
        default_factory=ConservationDiagnosticsConfig
    )

    def resolved_s_sigma_spec(self) -> SSigmaSpec:
        return self.s_sigma_spec or default_closed_s_sigma_spec(self.branch)

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["s_sigma_spec"] = self.resolved_s_sigma_spec().to_dict()
        return payload

    def config_hash(self) -> str:
        payload = self.to_dict()
        payload.pop("run_root", None)
        payload.pop("report_path", None)
        return stable_hash(payload)


def _tensor_linf(values: torch.Tensor) -> float:
    return float(torch.max(torch.abs(values)).detach().cpu().item())


def _wall_l2(values: torch.Tensor, wall_grid: WallGrid) -> float:
    return float(
        torch.sqrt(torch.sum(values * values * wall_grid.cell_widths)).detach().cpu().item()
    )


def independent_source_recompute(
    fields: ClosedCoupledFields,
    grid: TensorProductGrid,
    branch: BranchSmokeConfig,
) -> torch.Tensor:
    """Parallel recompute of the exact reduced source kernel for edit-drift checks."""

    density = fields.psi_real**2 + fields.psi_imag**2
    k1 = confinement_wall_to_matter_coefficient_torch(grid, branch, radius=fields.r0)
    return torch.sum(grid.radial_shell_volumes[:, None] * (-k1 * density), dim=0)


def independent_flux_divergence_recompute(
    r0: torch.Tensor,
    wall_grid: WallGrid,
    provider: SSigmaProvider,
    *,
    lower_value: float,
) -> torch.Tensor:
    """Parallel rebuild of the conservative flux kernel for edit-drift checks."""

    if r0.ndim != 1 or r0.shape[0] != wall_grid.spec.nw:
        raise ValueError("R0 must live on wall centers")
    lower = torch.as_tensor(lower_value, dtype=r0.dtype, device=r0.device)
    face_values = torch.empty(wall_grid.spec.nw + 1, dtype=r0.dtype, device=r0.device)
    face_values[0] = lower
    face_values[1:-1] = 0.5 * (r0[:-1] + r0[1:])
    face_values[-1] = r0[-1]
    t_faces = provider.T_w(face_values, wall_grid.w_faces)
    fluxes = torch.zeros(wall_grid.spec.nw + 1, dtype=r0.dtype, device=r0.device)
    fluxes[0] = t_faces[0] * (
        (-46.0 / 15.0) * lower
        + (15.0 / 4.0) * r0[0]
        - (5.0 / 6.0) * r0[1]
        + (3.0 / 20.0) * r0[2]
    ) / wall_grid.dw
    fluxes[1:-1] = t_faces[1:-1] * (r0[1:] - r0[:-1]) / wall_grid.dw
    fluxes[-1] = torch.zeros((), dtype=r0.dtype, device=r0.device)
    return -(fluxes[1:] - fluxes[:-1]) / wall_grid.cell_widths


def nonconservative_flux_divergence_discretization(
    r0: torch.Tensor,
    wall_grid: WallGrid,
    provider: SSigmaProvider,
) -> torch.Tensor:
    """Evaluate ``-d_w(T_w(R0, w) d_w R0)`` with center gradients only."""

    if r0.ndim != 1 or r0.shape[0] != wall_grid.spec.nw:
        raise ValueError("R0 must live on wall centers")
    gradient = wall_center_gradient(r0, wall_grid)
    center_flux = provider.T_w(r0, wall_grid.w_centers) * gradient
    return -wall_center_gradient(center_flux, wall_grid)


def balance_validation_diagnostic(
    state: torch.Tensor,
    grid: TensorProductGrid,
    branch: BranchSmokeConfig,
    spec: SSigmaSpec,
    *,
    final_residual_linf: float,
    final_tolerance: float,
    nontrivial_floor: float,
    source_abs_tol: float,
    flux_abs_tol: float,
) -> dict[str, Any]:
    fields, _ = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    provider = resolve_s_sigma(spec)
    terms = patha_closed_wall_terms(fields, grid, branch, s_sigma=provider)
    wall_grid = WallGrid.create(
        WallGridSpec(w_min=grid.spec.w_min, w_max=grid.spec.w_max, nw=grid.spec.nw),
        dtype=grid.dtype,
        device=grid.device,
    )
    term_tensors = {
        "flux_divergence": terms.flux_divergence,
        "gradient_square": terms.gradient_square,
        "potential_gradient": terms.potential_gradient,
        "source": terms.source,
    }
    term_rows = [
        {
            "term": name,
            "linf": _tensor_linf(values),
            "l2": _wall_l2(values, wall_grid),
        }
        for name, values in term_tensors.items()
    ]
    dominant = max(row["linf"] for row in term_rows)
    for row in term_rows:
        row["relative_to_dominant"] = row["linf"] / max(dominant, 1.0e-300)
    residual_linf = _tensor_linf(terms.residual)
    residual_l2 = _wall_l2(terms.residual, wall_grid)
    residual_relative = residual_linf / max(dominant, 1.0e-300)

    recomputed_source = independent_source_recompute(fields, grid, branch)
    recomputed_flux = independent_flux_divergence_recompute(
        fields.r0,
        wall_grid,
        provider,
        lower_value=branch.r_mouth,
    )
    source_diff = recomputed_source - terms.source
    flux_diff = recomputed_flux - terms.flux_divergence
    solver_floor = max(final_residual_linf, final_tolerance)
    checks = {
        "balance_terms_nontrivial": all(row["linf"] > nontrivial_floor for row in term_rows),
    }
    identity_checks = {
        "independent_source_recompute_matches_not_a_physics_gate": (
            _tensor_linf(source_diff) <= source_abs_tol
        ),
        "independent_flux_recompute_matches_not_a_physics_gate": (
            _tensor_linf(flux_diff) <= flux_abs_tol
        ),
        "balance_residual_relative_below_solver_floor_not_a_physics_gate": (
            residual_relative <= solver_floor
        ),
    }
    return {
        "term_rows": term_rows,
        "dominant_linf": dominant,
        "nontrivial_floor_reference": nontrivial_floor,
        "residual_linf": residual_linf,
        "residual_l2": residual_l2,
        "residual_relative_to_dominant": residual_relative,
        "solver_floor_reference": solver_floor,
        "newton_final_residual_linf": final_residual_linf,
        "newton_final_tolerance": final_tolerance,
        "source_recompute_max_abs_diff": _tensor_linf(source_diff),
        "flux_recompute_max_abs_diff": _tensor_linf(flux_diff),
        "checks": checks,
        "identity_checks": identity_checks,
        "passed": all(checks.values()),
        "diagnostic_note": (
            "The residual-floor row is reported only as a solve diagnostic because "
            "the Newton residual already contains the wall residual. The source and "
            "conservative-flux recomputes are parallel reconstructions of the same "
            "kernels: they pin edit drift and wiring, but do not establish operator "
            "correctness as physics gates."
        ),
    }


def _closed_to_coupled_state(
    state: torch.Tensor,
    grid: TensorProductGrid,
) -> torch.Tensor:
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    assert mu is not None
    return pack_coupled_fields(
        CoupledFields(
            psi_real=fields.psi_real,
            psi_imag=fields.psi_imag,
            a0=fields.a0,
            ar=fields.ar,
            aw=fields.aw,
        ),
        mu,
    )


def _solve_closed_level(
    *,
    branch: BranchSmokeConfig,
    spec: SSigmaSpec,
    grid: TensorProductGrid,
    initial_state: torch.Tensor | None,
    grid_name: str,
) -> tuple[torch.Tensor, dict[str, Any]]:
    provider = resolve_s_sigma(spec)
    boundaries = branch_boundary_conditions(branch)
    state = (
        initial_state.detach()
        if initial_state is not None
        else initial_closed_branch_state(grid, branch, dtype=grid.dtype, device=grid.device)
    )
    started = time.perf_counter()
    stages: list[dict[str, Any]] = []
    converged = True
    message = "closed validation continuation completed"
    stop_conditions: list[dict[str, Any]] = []
    shared_preconditioner_factory = make_closed_coupled_colored_sparse_jacobian_lu_factory(grid)

    for eos_K in branch.continuation_K_values:
        fields_before, _ = unpack_closed_coupled_fields(
            state,
            grid,
            has_chemical_potential=True,
        )
        if float(torch.min(fields_before.r0).detach().cpu().item()) <= 0.0:
            converged = False
            message = f"stopped before eos_K={eos_K}: R0 became nonpositive"
            stop_conditions.append({"condition": "R0_nonpositive", "tripped": True})
            break
        residual_fn = lambda x, eos_K=eos_K: patha_closed_branch_residual(
            x,
            grid,
            branch,
            eos_K=eos_K,
            boundaries=boundaries,
            s_sigma=provider,
        )
        result = solve_newton_jvp(
            residual_fn,
            state,
            branch.newton,
            preconditioner_factory=(
                make_closed_coupled_colored_sparse_jacobian_lu_factory(grid)
                if branch.newton.preconditioner.rebuild_policy == "once_per_newton_solve"
                else shared_preconditioner_factory
            ),
        )
        state = result.x.detach()
        fields_after, _ = unpack_closed_coupled_fields(
            state,
            grid,
            has_chemical_potential=True,
        )
        r0_min = float(torch.min(fields_after.r0).detach().cpu().item())
        r0_max = float(torch.max(fields_after.r0).detach().cpu().item())
        stages.append(
            {
                "eos_K": eos_K,
                "converged": result.converged,
                "iterations": result.iterations,
                "initial_residual_norm": result.initial_residual_norm,
                "final_residual_norm": result.final_residual_norm,
                "tolerance": result.tolerance,
                "message": result.message,
                "gmres_iterations": [
                    row.gmres_iterations
                    for row in result.history
                    if row.gmres_iterations is not None
                ],
                "r0_min": r0_min,
                "r0_max": r0_max,
            }
        )
        if r0_min <= 0.0:
            converged = False
            message = f"stopped after eos_K={eos_K}: R0 became nonpositive"
            stop_conditions.append({"condition": "R0_nonpositive", "tripped": True})
            break
        if not result.converged:
            converged = False
            message = f"continuation stopped at eos_K={eos_K}: {result.message}"
            break

    final_stage_index = max(0, min(len(stages), len(branch.continuation_K_values)) - 1)
    final_eos_K = branch.continuation_K_values[final_stage_index]
    final_residual = patha_closed_branch_residual(
        state,
        grid,
        branch,
        eos_K=final_eos_K,
        boundaries=boundaries,
        s_sigma=provider,
    ).detach()
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    density = fields.psi_real**2 + fields.psi_imag**2
    final_residual_linf = _tensor_linf(final_residual)
    final_tolerance = stages[-1]["tolerance"] if stages else branch.newton.residual_atol
    elapsed = time.perf_counter() - started
    if not any(row["condition"] == "R0_nonpositive" for row in stop_conditions):
        stop_conditions.append({"condition": "R0_nonpositive", "tripped": False})
    summary = {
        "grid": grid_name,
        "nr": grid.spec.nr,
        "nw": grid.spec.nw,
        "dof": int(state.numel()),
        "wall_clock_seconds": elapsed,
        "peak_memory_mb": _max_rss_mb(),
        "converged": converged,
        "message": message,
        "final_eos_K": final_eos_K,
        "final_residual_linf": final_residual_linf,
        "final_tolerance": final_tolerance,
        "final_mass": float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
        "chemical_potential": float(mu.detach().cpu().item()) if mu is not None else None,
        "r0_min": float(torch.min(fields.r0).detach().cpu().item()),
        "r0_max": float(torch.max(fields.r0).detach().cpu().item()),
        "stages": stages,
        "stop_conditions": stop_conditions,
    }
    return state, summary


def _restrict_wall_1d(
    fine_values: torch.Tensor,
    fine_grid: TensorProductGrid,
    coarse_grid: TensorProductGrid,
) -> torch.Tensor:
    ratio = fine_grid.spec.nw // coarse_grid.spec.nw
    if fine_grid.spec.nw != coarse_grid.spec.nw * ratio or ratio < 1:
        raise ValueError("fine wall grid must be an integer refinement of coarse grid")
    blocks = fine_values.reshape(coarse_grid.spec.nw, ratio)
    widths = torch.full_like(fine_values, fine_grid.dw).reshape(coarse_grid.spec.nw, ratio)
    return torch.sum(blocks * widths, dim=1) / torch.sum(widths, dim=1)


def conservative_restrict_closed_fields(
    fine_state: torch.Tensor,
    fine_grid: TensorProductGrid,
    coarse_grid: TensorProductGrid,
) -> dict[str, torch.Tensor]:
    fields, _ = unpack_closed_coupled_fields(
        fine_state,
        fine_grid,
        has_chemical_potential=True,
    )
    return {
        "psi_real": _volume_restrict_2d(fields.psi_real, fine_grid, coarse_grid),
        "psi_imag": _volume_restrict_2d(fields.psi_imag, fine_grid, coarse_grid),
        "a0": _volume_restrict_2d(fields.a0, fine_grid, coarse_grid),
        "ar": _volume_restrict_2d(fields.ar, fine_grid, coarse_grid),
        "aw": _volume_restrict_2d(fields.aw, fine_grid, coarse_grid),
        "r0": _restrict_wall_1d(fields.r0, fine_grid, coarse_grid),
    }


def closed_raw_field_self_difference(
    coarse_state: torch.Tensor,
    coarse_grid: TensorProductGrid,
    fine_state: torch.Tensor,
    fine_grid: TensorProductGrid,
) -> dict[str, float]:
    coarse_fields, _ = unpack_closed_coupled_fields(
        coarse_state,
        coarse_grid,
        has_chemical_potential=True,
    )
    restricted = conservative_restrict_closed_fields(fine_state, fine_grid, coarse_grid)
    field_map = {
        "psi_real": coarse_fields.psi_real,
        "psi_imag": coarse_fields.psi_imag,
        "a0": coarse_fields.a0,
        "ar": coarse_fields.ar,
        "aw": coarse_fields.aw,
    }
    diff_density = torch.zeros_like(coarse_grid.cell_volumes)
    ref_density = torch.zeros_like(coarse_grid.cell_volumes)
    linf = 0.0
    for name, coarse_values in field_map.items():
        diff = coarse_values - restricted[name]
        diff_density = diff_density + diff**2
        ref_density = ref_density + coarse_values**2
        linf = max(linf, _tensor_linf(diff))
    field_l2 = float(
        torch.sqrt(torch.sum(diff_density * coarse_grid.cell_volumes)).detach().cpu().item()
    )
    field_reference = float(
        torch.sqrt(torch.sum(ref_density * coarse_grid.cell_volumes)).detach().cpu().item()
    )
    wall_widths = torch.full_like(coarse_fields.r0, coarse_grid.dw)
    r0_diff = coarse_fields.r0 - restricted["r0"]
    r0_l2 = float(torch.sqrt(torch.sum(r0_diff * r0_diff * wall_widths)).detach().cpu().item())
    r0_reference = float(
        torch.sqrt(torch.sum(coarse_fields.r0 * coarse_fields.r0 * wall_widths))
        .detach()
        .cpu()
        .item()
    )
    linf = max(linf, _tensor_linf(r0_diff))
    combined_l2 = math.sqrt(field_l2 * field_l2 + r0_l2 * r0_l2)
    combined_reference = math.sqrt(field_reference * field_reference + r0_reference * r0_reference)
    return {
        "raw_field_l2_change": combined_l2,
        "raw_field_relative_l2_change": combined_l2 / max(combined_reference, 1.0e-300),
        "raw_field_linf_change": linf,
        "matter_gauge_l2_change": field_l2,
        "matter_gauge_relative_l2_change": field_l2 / max(field_reference, 1.0e-300),
        "r0_l2_change": r0_l2,
        "r0_relative_l2_change": r0_l2 / max(r0_reference, 1.0e-300),
        "r0_linf_change": _tensor_linf(r0_diff),
    }


def closed_surrogate_observables(
    state: torch.Tensor,
    grid: TensorProductGrid,
    branch: BranchSmokeConfig,
    *,
    final_residual_linf: float,
    balance_residual_relative: float,
) -> dict[str, float]:
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    density = fields.psi_real**2 + fields.psi_imag**2
    raw_density = (
        fields.psi_real**2
        + fields.psi_imag**2
        + fields.a0**2
        + fields.ar**2
        + fields.aw**2
    )
    wall_widths = torch.full_like(fields.r0, grid.dw)
    return {
        "density_mass": float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
        "raw_matter_gauge_l2_norm": float(
            torch.sqrt(torch.sum(raw_density * grid.cell_volumes)).detach().cpu().item()
        ),
        "r0_l2_norm": float(
            torch.sqrt(torch.sum(fields.r0 * fields.r0 * wall_widths)).detach().cpu().item()
        ),
        "r0_min": float(torch.min(fields.r0).detach().cpu().item()),
        "r0_max": float(torch.max(fields.r0).detach().cpu().item()),
        "r0_range": float((torch.max(fields.r0) - torch.min(fields.r0)).detach().cpu().item()),
        "chemical_potential": float(mu.detach().cpu().item()) if mu is not None else float("nan"),
        "final_residual_linf": float(final_residual_linf),
        "balance_residual_relative_to_dominant": float(balance_residual_relative),
        "mass_request": float(branch.mass),
    }


def _self_convergence_rows(
    solved_levels: list[dict[str, Any]],
    *,
    ratio: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index in range(1, len(solved_levels)):
        coarse = solved_levels[index - 1]
        fine = solved_levels[index]
        diff = closed_raw_field_self_difference(
            coarse["state"],
            coarse["grid_object"],
            fine["state"],
            fine["grid_object"],
        )
        rows.append(
            {
                "coarse_grid": coarse["grid"],
                "fine_grid": fine["grid"],
                "coarse_dof": coarse["dof"],
                "fine_dof": fine["dof"],
                **diff,
                "observed_order_from_raw_l2_change": None,
                "observed_order_from_r0_l2_change": None,
            }
        )
    for index in range(1, len(rows)):
        rows[index]["observed_order_from_raw_l2_change"] = observed_order_from_three(
            rows[index - 1]["raw_field_l2_change"],
            0.0,
            rows[index]["raw_field_l2_change"],
            ratio,
        )
        rows[index]["observed_order_from_r0_l2_change"] = observed_order_from_three(
            rows[index - 1]["r0_l2_change"],
            0.0,
            rows[index]["r0_l2_change"],
            ratio,
        )
    return rows


def _observable_table(
    level_rows: list[dict[str, Any]],
    *,
    ratio: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if not level_rows:
        return [], []
    names = list(level_rows[0]["surrogate_values"].keys())
    summary_rows: list[dict[str, Any]] = []
    table_rows: list[dict[str, Any]] = []
    for name in names:
        values = [float(row["surrogate_values"][name]) for row in level_rows]
        changes = [None]
        for left, right in zip(values, values[1:]):
            changes.append(abs(right - left))
        orders: list[float | None] = [None, None]
        for index in range(2, len(values)):
            orders.append(observed_order_from_three(values[index - 2], values[index - 1], values[index], ratio))
        last_order = next((order for order in reversed(orders) if order is not None), None)
        continuum = (
            richardson_estimate(values[-2], values[-1], ratio, last_order)
            if len(values) >= 2
            else None
        )
        summary_rows.append(
            {
                "observable": name,
                "finest_grid": level_rows[-1]["grid"],
                "finest_value": values[-1],
                "last_observed_order": last_order,
                "richardson_estimate": continuum,
                "finest_error_estimate": (
                    abs(values[-1] - continuum) if continuum is not None else None
                ),
            }
        )
        for index, row in enumerate(level_rows):
            table_rows.append(
                {
                    "observable": name,
                    "grid": row["grid"],
                    "value": values[index],
                    "successive_change": changes[index],
                    "observed_order": orders[index],
                }
            )
    return summary_rows, table_rows


def _raw_difference_decreases_or_at_floor(
    self_rows: list[dict[str, Any]],
    *,
    floor: float,
) -> bool:
    if len(self_rows) < 2:
        return False
    raw_ok = (
        self_rows[-1]["raw_field_relative_l2_change"]
        <= self_rows[-2]["raw_field_relative_l2_change"]
        or self_rows[-1]["raw_field_relative_l2_change"] <= floor
    )
    r0_ok = (
        self_rows[-1]["r0_relative_l2_change"] <= self_rows[-2]["r0_relative_l2_change"]
        or self_rows[-1]["r0_relative_l2_change"] <= floor
    )
    return bool(raw_ok and r0_ok)


def _balance_stable(level_rows: list[dict[str, Any]]) -> bool:
    if len(level_rows) < 3:
        return False
    for row in level_rows:
        balance = row["surrogate_values"]["balance_residual_relative_to_dominant"]
        floor = max(row["final_residual_linf"], row["final_tolerance"])
        if balance > floor:
            return False
    return True


def _flux_discretization_row(
    *,
    state: torch.Tensor,
    grid: TensorProductGrid,
    branch: BranchSmokeConfig,
    spec: SSigmaSpec,
    final_residual_linf: float,
    final_tolerance: float,
    h2_tolerance_factor: float,
) -> dict[str, Any]:
    fields, _ = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    provider = resolve_s_sigma(spec)
    terms = patha_closed_wall_terms(fields, grid, branch, s_sigma=provider)
    wall_grid = WallGrid.create(
        WallGridSpec(w_min=grid.spec.w_min, w_max=grid.spec.w_max, nw=grid.spec.nw),
        dtype=grid.dtype,
        device=grid.device,
    )
    nonconservative = nonconservative_flux_divergence_discretization(
        fields.r0,
        wall_grid,
        provider,
    )
    interior = slice(1, -1)
    diff = nonconservative[interior] - terms.flux_divergence[interior]
    conservative_reference = _tensor_linf(terms.flux_divergence[interior])
    nonconservative_reference = _tensor_linf(nonconservative[interior])
    reference = max(conservative_reference, nonconservative_reference, 1.0e-300)
    solver_floor = max(final_residual_linf, final_tolerance)
    refinement_tolerance = max(
        solver_floor,
        h2_tolerance_factor * grid.dw * grid.dw * reference,
    )
    interior_diff = _tensor_linf(diff)
    return {
        "grid": f"nr_{grid.spec.nr}_nw_{grid.spec.nw}",
        "nw": grid.spec.nw,
        "spacing": grid.dw,
        "interior_cells_compared": max(grid.spec.nw - 2, 0),
        "interior_exclusion": "mouth cell and exit cell excluded",
        "interior_max_abs_diff": interior_diff,
        "conservative_flux_linf_interior": conservative_reference,
        "nonconservative_flux_linf_interior": nonconservative_reference,
        "h2_scaled_tolerance": refinement_tolerance,
        "solver_floor_reference": solver_floor,
        "observed_order": None,
        "decreases_or_at_floor": None,
        "finest_within_tolerance": None,
    }


def _flux_discretization_convergence_rows(
    solved_levels: list[dict[str, Any]],
    *,
    branch: BranchSmokeConfig,
    spec: SSigmaSpec,
    ratio: int,
    h2_tolerance_factor: float,
) -> list[dict[str, Any]]:
    rows = [
        _flux_discretization_row(
            state=level["state"],
            grid=level["grid_object"],
            branch=branch,
            spec=spec,
            final_residual_linf=level["final_residual_linf"],
            final_tolerance=level["final_tolerance"],
            h2_tolerance_factor=h2_tolerance_factor,
        )
        for level in solved_levels
    ]
    for index, row in enumerate(rows):
        if index == 0:
            row["decreases_or_at_floor"] = True
            continue
        previous = rows[index - 1]
        row["observed_order"] = observed_order_from_three(
            previous["interior_max_abs_diff"],
            0.0,
            row["interior_max_abs_diff"],
            ratio,
        )
        row["decreases_or_at_floor"] = (
            row["interior_max_abs_diff"] <= previous["interior_max_abs_diff"]
            or row["interior_max_abs_diff"] <= row["solver_floor_reference"]
        )
    if rows:
        rows[-1]["finest_within_tolerance"] = (
            rows[-1]["interior_max_abs_diff"] <= rows[-1]["h2_scaled_tolerance"]
        )
    return rows


def _flux_discretization_converges(rows: list[dict[str, Any]]) -> bool:
    if len(rows) < 2:
        return False
    return bool(
        all(row["decreases_or_at_floor"] for row in rows[1:])
        and rows[-1]["finest_within_tolerance"]
    )


def _strip_runtime(level: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in level.items()
        if key not in {"state", "grid_object"}
    }


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return list(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _run_convergence_ladder(
    config: PathAClosedValidationConfig,
    *,
    dtype: torch.dtype,
) -> dict[str, Any]:
    validate_refinement_ladder(config.convergence_levels, config.refinement_ratio)
    branch = config.branch
    spec = config.resolved_s_sigma_spec()
    previous_grid: TensorProductGrid | None = None
    previous_state: torch.Tensor | None = None
    solved_levels: list[dict[str, Any]] = []
    stop_reason = "completed configured ladder"

    for level_index, level in enumerate(config.convergence_levels):
        grid = _create_branch_grid(branch, level, dtype=dtype, device="cpu")
        initial = (
            resample_closed_branch_state(previous_state, previous_grid, grid, branch)
            if previous_state is not None and previous_grid is not None
            else None
        )
        state, summary = _solve_closed_level(
            branch=branch,
            spec=spec,
            grid=grid,
            initial_state=initial,
            grid_name=f"closed_l{level_index}_nr_{level[0]}_nw_{level[1]}",
        )
        balance = balance_validation_diagnostic(
            state,
            grid,
            branch,
            spec,
            final_residual_linf=summary["final_residual_linf"],
            final_tolerance=summary["final_tolerance"],
            nontrivial_floor=config.nontrivial_term_floor,
            source_abs_tol=config.independent_source_abs_tol,
            flux_abs_tol=config.independent_flux_abs_tol,
        )
        summary["balance_validation"] = balance
        summary["surrogate_values"] = closed_surrogate_observables(
            state,
            grid,
            branch,
            final_residual_linf=summary["final_residual_linf"],
            balance_residual_relative=balance["residual_relative_to_dominant"],
        )
        solved_levels.append(
            {
                **summary,
                "level": level_index,
                "spacing": max(grid.dr, grid.dw),
                "state": state.detach(),
                "grid_object": grid,
            }
        )
        if not summary["converged"]:
            stop_reason = summary["message"]
            break
        if summary["r0_min"] <= 0.0:
            stop_reason = f"stopped after {summary['grid']}: R0 became nonpositive"
            break
        previous_grid = grid
        previous_state = state.detach()
        if summary["wall_clock_seconds"] > config.max_level_seconds:
            stop_reason = (
                f"stopped after {summary['grid']} because the level exceeded "
                f"{config.max_level_seconds:.1f}s"
            )
            break

    level_rows = [_strip_runtime(level) for level in solved_levels]
    self_rows = _self_convergence_rows(solved_levels, ratio=config.refinement_ratio)
    flux_discretization_rows = _flux_discretization_convergence_rows(
        solved_levels,
        branch=branch,
        spec=spec,
        ratio=config.refinement_ratio,
        h2_tolerance_factor=config.independent_flux_h2_tolerance_factor,
    )
    observable_summary, observable_table = _observable_table(
        level_rows,
        ratio=config.refinement_ratio,
    )
    pass_checks = {
        "minimum_three_levels": len(level_rows) >= config.min_levels,
        "all_completed_levels_converged": all(row["converged"] for row in level_rows),
        "fixed_ratio_ladder": True,
        "raw_field_self_difference_decreases_or_at_floor": _raw_difference_decreases_or_at_floor(
            self_rows,
            floor=config.raw_difference_floor,
        ),
        "independent_flux_discretization_converges": _flux_discretization_converges(
            flux_discretization_rows
        ),
    }
    identity_checks = {
        "balance_residual_relative_stable_under_refinement_not_a_physics_gate": (
            _balance_stable(level_rows)
        ),
    }
    return {
        "level_rows": level_rows,
        "self_convergence_rows": self_rows,
        "flux_discretization_rows": flux_discretization_rows,
        "observable_summary": observable_summary,
        "observable_table": observable_table,
        "stop_reason": stop_reason,
        "pass_checks": pass_checks,
        "identity_checks": identity_checks,
        "passed": all(pass_checks.values()),
        "solved_levels": solved_levels,
    }


def _metric_max(rows: list[dict[str, Any]], key: str) -> float:
    return max((abs(float(row[key])) for row in rows), default=0.0)


def _conservation_metrics(diagnostics: dict[str, Any]) -> dict[str, float]:
    return {
        "local_l2_relative_max": _metric_max(
            diagnostics["local_residual_rows"],
            "interior_l2_relative",
        ),
        "local_linf_relative_max": _metric_max(
            diagnostics["local_residual_rows"],
            "interior_linf_relative",
        ),
        "independent_gauss_relative_max": _metric_max(
            diagnostics["gauss_closure_rows"],
            "relative_residual",
        ),
        "operator_maxwell_absolute_max": _metric_max(
            diagnostics["maxwell_residual_closure_rows"],
            "absolute_residual",
        ),
        "fv_budget_closure_absolute_max": _metric_max(
            diagnostics["budget_rows"],
            "closure_absolute",
        ),
    }


def _conservation_metric_activity() -> dict[str, dict[str, str]]:
    return {
        "local_l2_relative_max": {
            "activity": "STRUCTURALLY ZERO",
            "basis": "stationary current sectors have zero transport signal on this branch",
        },
        "local_linf_relative_max": {
            "activity": "STRUCTURALLY ZERO",
            "basis": "stationary current sectors have zero transport signal on this branch",
        },
        "independent_gauss_relative_max": {
            "activity": "LIVE",
            "basis": "center-gradient Gauss reconstruction independent of operator fluxes",
        },
        "operator_maxwell_absolute_max": {
            "activity": "SOLVER-FLOOR IDENTITY",
            "basis": "operator Gauss residual is bounded by the converged Newton residual",
        },
        "fv_budget_closure_absolute_max": {
            "activity": "STRUCTURALLY ZERO",
            "basis": "finite-volume budget identity closes algebraically for zero flux/sink",
        },
    }


def _run_conservation_validation(
    config: PathAClosedValidationConfig,
    finest: dict[str, Any],
) -> dict[str, Any]:
    branch = config.branch
    grid = finest["grid_object"]
    closed_state = finest["state"]
    closed_coupled_state = _closed_to_coupled_state(closed_state, grid)
    frozen_config = HarnessConfig(
        branch=branch,
        run_root=str(Path(config.run_root) / "frozen_geometry_baseline"),
        report_path=str(Path(config.run_root) / "frozen_geometry_baseline.md"),
    )
    frozen_state, frozen_summary = run_branch_continuation(
        frozen_config,
        grid,
        initial_state=closed_coupled_state,
        grid_name=f"frozen_baseline_nr_{grid.spec.nr}_nw_{grid.spec.nw}",
    )
    closed_diag = evaluate_conservation_state(
        mode="closed_r0_solved",
        level_index=0,
        state=closed_coupled_state,
        grid=grid,
        cfg=branch,
        study=config.conservation,
    )
    frozen_diag = evaluate_conservation_state(
        mode="frozen_geometry_baseline",
        level_index=0,
        state=frozen_state,
        grid=grid,
        cfg=branch,
        study=config.conservation,
    )
    closed_metrics = _conservation_metrics(closed_diag)
    frozen_metrics = _conservation_metrics(frozen_diag)
    thresholds = {
        key: max(
            config.conservation_floor,
            config.conservation_baseline_factor * frozen_metrics[key],
        )
        for key in closed_metrics
    }
    max_region_volume = max(
        (float(row["region_volume"]) for row in closed_diag["maxwell_residual_closure_rows"]),
        default=1.0,
    )
    thresholds["operator_maxwell_absolute_max"] = max(
        thresholds["operator_maxwell_absolute_max"],
        finest["final_residual_linf"] * max_region_volume,
        frozen_summary["final_residual_linf"] * max_region_volume,
    )
    budget_identity_closed = (
        closed_metrics["fv_budget_closure_absolute_max"]
        <= config.conservation.budget_closure_abs_tol
    )
    budget_identity_frozen = (
        frozen_metrics["fv_budget_closure_absolute_max"]
        <= config.conservation.budget_closure_abs_tol
    )
    metrics_within = all(closed_metrics[key] <= thresholds[key] for key in closed_metrics)
    null_labels = all(
        row["label"] in {"null diagnostic", "solver-floor diagnostic"}
        for row in closed_diag["local_residual_rows"]
    )
    pass_checks = {
        "closed_conservation_consistent_with_frozen_baseline": metrics_within,
        "stationary_current_sectors_labelled_null_or_floor": null_labels,
    }
    identity_checks = {
        "fv_budget_identity_roundoff_not_a_physics_gate": (
            budget_identity_closed and budget_identity_frozen
        ),
        "closed_operator_gauss_residual_at_solver_floor_not_a_physics_gate": (
            closed_metrics["operator_maxwell_absolute_max"]
            <= thresholds["operator_maxwell_absolute_max"]
        ),
    }
    return {
        "closed_metrics": closed_metrics,
        "frozen_metrics": frozen_metrics,
        "thresholds": thresholds,
        "metric_activity": _conservation_metric_activity(),
        "baseline_factor": config.conservation_baseline_factor,
        "baseline_factor_note": (
            "The 4x frozen-baseline envelope is a sanity bound for the closed solve. "
            "The live conservation signal is the independent center-gradient Gauss "
            "reconstruction; structurally-zero channels are disclosed separately."
        ),
        "r0_framing_note": (
            "Number/charge and independent Gauss closure are evaluated on the solved "
            "matter/gauge fields. The energy-density confinement term uses frozen "
            "geometry because conservation_diagnostics._energy_density calls "
            "confinement_potential_torch without radius; closed and frozen branches "
            "therefore compare that channel apples-to-apples."
        ),
        "null_floor_label_note": (
            "null_floor_label keys on transport-signal norm, not residual magnitude; "
            "that is correct for this isotropic stationary branch and should be "
            "revisited for a current-carrying branch."
        ),
        "closed_diagnostics": closed_diag,
        "frozen_diagnostics": frozen_diag,
        "frozen_solve": {
            "grid": frozen_summary["grid"],
            "dof": frozen_summary["dof"],
            "final_residual_linf": frozen_summary["final_residual_linf"],
            "converged": frozen_summary["converged"],
            "message": frozen_summary["message"],
            "wall_clock_seconds": frozen_summary["wall_clock_seconds"],
        },
        "pass_checks": pass_checks,
        "identity_checks": identity_checks,
        "passed": all(pass_checks.values()),
    }


def run_patha_closed_validation(
    config: PathAClosedValidationConfig | None = None,
) -> dict[str, Any]:
    if config is None:
        config = PathAClosedValidationConfig()
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(BackendConfig())
    convergence = _run_convergence_ladder(config, dtype=dtype)
    solved_levels = convergence["solved_levels"]
    conservation = (
        _run_conservation_validation(config, solved_levels[-1])
        if solved_levels
        else {
            "pass_checks": {"closed_conservation_consistent_with_frozen_baseline": False},
            "identity_checks": {},
            "passed": False,
        }
    )
    finest = solved_levels[-1] if solved_levels else None
    balance = finest["balance_validation"] if finest is not None else None
    report_path = Path(config.report_path)
    token_scan = target_token_scan(validation_scan_paths(report_path if report_path.exists() else None))
    balance_checks = balance["checks"] if balance is not None else {"balance_validation_ran": False}
    balance_identity_checks = balance["identity_checks"] if balance is not None else {}
    convergence_identity_checks = convergence.get("identity_checks", {})
    conservation_identity_checks = conservation.get("identity_checks", {})
    identity_checks = {
        **balance_identity_checks,
        **convergence_identity_checks,
        **conservation_identity_checks,
    }
    pass_checks = {
        **balance_checks,
        **convergence["pass_checks"],
        **conservation["pass_checks"],
        "target_token_scan_clean": bool(token_scan["passed"]),
    }
    stop_conditions = {
        "edit_drift_pin_mismatch_not_a_physics_gate": not (
            balance_identity_checks.get(
                "independent_source_recompute_matches_not_a_physics_gate",
                False,
            )
            and balance_identity_checks.get(
                "independent_flux_recompute_matches_not_a_physics_gate",
                False,
            )
        ),
        "independent_flux_discretization_not_converging": not convergence["pass_checks"].get(
            "independent_flux_discretization_converges",
            False,
        ),
        "balance_terms_degenerate": not balance_checks.get("balance_terms_nontrivial", False),
        "fewer_than_three_levels_completed": not convergence["pass_checks"]["minimum_three_levels"],
        "r0_nonpositive": any(row["r0_min"] <= 0.0 for row in convergence["level_rows"]),
        "closed_conservation_materially_worse_than_frozen": not conservation["pass_checks"].get(
            "closed_conservation_consistent_with_frozen_baseline",
            False,
        ),
    }
    public_convergence = {
        key: value for key, value in convergence.items() if key != "solved_levels"
    }
    results = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "s_sigma_digest": config.resolved_s_sigma_spec().digest(),
        "method": {
            "placeholder_family": config.resolved_s_sigma_spec().family,
            "refinement_ratio": config.refinement_ratio,
            "convergence_levels": list(config.convergence_levels),
            "preconditioner": asdict(config.branch.newton.preconditioner),
            "balance_framing": (
                "The Newton residual contains the wall residual, so residual-floor "
                "rows are reported as solve diagnostics only. Counted balance evidence "
                "is the non-degenerate term check plus an interior-only comparison "
                "between conservative face-flux divergence and a separate "
                "non-conservative center-gradient discretization."
            ),
            "source_framing": (
                "The return source is an exact reduced wall-to-matter kernel, not a "
                "differential operator with an alternate discretization. Its parallel "
                "recompute is an edit-drift pin; correctness rests on reciprocity with "
                "the forward wall-to-matter force, the chunk-1a dual-engine MMS, and "
                "the fidelity audit."
            ),
            "conservation_framing": (
                "Number/charge and the independent Gauss closure are evaluated on the "
                "solved matter/gauge fields. The energy-density confinement term uses "
                "frozen geometry because confinement_potential_torch is called without "
                "radius in conservation_diagnostics._energy_density; both branches use "
                "that same channel for an apples-to-apples comparison."
            ),
            "null_floor_label_note": (
                "null_floor_label keys on transport-signal norm rather than residual "
                "magnitude; this is correct for the isotropic stationary branch and "
                "needs revisiting for a current-carrying branch."
            ),
        },
        "balance_validation": balance,
        "convergence": public_convergence,
        "conservation": conservation,
        "target_token_scan": token_scan,
        "pass_checks": pass_checks,
        "identity_checks": identity_checks,
        "stop_conditions": stop_conditions,
        "passed": all(pass_checks.values()),
    }
    table_path = Path(config.run_root) / config.name / "closed_validation_table.json"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    results["machine_readable_table"] = str(table_path)
    table_path.write_text(
        json.dumps(results, indent=2, sort_keys=True, default=_json_default),
        encoding="utf-8",
    )
    return results


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.6e}"
    if isinstance(value, (tuple, list)):
        return "[" + ", ".join(str(item) for item in value) + "]"
    return str(value)


def _table(headers: list[str], rows: list[dict[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(header)) for header in headers) + " |")
    return "\n".join(lines)


def _gate_rows(checks: dict[str, bool]) -> list[dict[str, Any]]:
    return [{"gate": key, "status": "PASS" if value else "FAIL"} for key, value in checks.items()]


def _conservation_metric_rows(conservation: dict[str, Any]) -> list[dict[str, Any]]:
    closed = conservation.get("closed_metrics", {})
    frozen = conservation.get("frozen_metrics", {})
    thresholds = conservation.get("thresholds", {})
    activity = conservation.get("metric_activity", {})
    return [
        {
            "metric": key,
            "activity": activity.get(key, {}).get("activity"),
            "closed_r0_solved": closed.get(key),
            "frozen_geometry_baseline": frozen.get(key),
            "threshold": thresholds.get(key),
            "basis": activity.get(key, {}).get("basis"),
            "status": (
                "PASS"
                if key in closed and key in thresholds and closed[key] <= thresholds[key]
                else "FAIL"
            ),
        }
        for key in closed
    ]


def _gauss_closure_rows(conservation: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "surface": row["surface"],
            "reconstruction": row["reconstruction"],
            "relative_residual": row["relative_residual"],
            "absolute_residual": row["absolute_residual"],
            "surface_flux": row["surface_flux"],
            "enclosed_mu0_charge": row["enclosed_mu0_charge"],
        }
        for row in conservation.get("closed_diagnostics", {}).get("gauss_closure_rows", [])
    ]


def write_patha_closed_validation_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    config = results["config"]
    branch = config["branch"]
    balance = results["balance_validation"]
    convergence = results["convergence"]
    conservation = results["conservation"]
    lines: list[str] = []
    lines.append("# Path-A Chunk 1c Self-Consistent Closed Validation")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append(f"S_Sigma digest: `{results['s_sigma_digest']}`")
    lines.append("")
    lines.append(
        "Scope: target-blind engineering validation of the closed placeholder solve. "
        "No calibration, physical packet, coefficient extraction, or export is performed."
    )
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append(f"Closed placeholder label: {branch['placeholder_label']}")
    lines.append(
        f"Refinement ladder: `{config['convergence_levels']}` with ratio "
        f"`{config['refinement_ratio']}`."
    )
    lines.append(
        "Preconditioner rebuild policy for validation ladder: "
        f"`{branch['newton']['preconditioner']['rebuild_policy']}`. The true matrix-free "
        "JVP residual remains the solve operator."
    )
    lines.append("")
    lines.append("## Self-Consistent Balance Terms")
    lines.append("")
    lines.append(
        "The counted balance check here is term non-degeneracy. Residual-floor rows are "
        "solve diagnostics because the Newton residual already contains the wall residual. "
        "The source and conservative flux recomputes are parallel reconstructions of the "
        "same kernels: useful edit-drift pins, but not physics gates."
    )
    lines.append("")
    if balance is not None:
        lines.append(
            _table(
                ["term", "linf", "l2", "relative_to_dominant"],
                balance["term_rows"],
            )
        )
        lines.append("")
        lines.append(
            "Dominant term Linf "
            f"`{balance['dominant_linf']:.6e}`; balance residual Linf "
            f"`{balance['residual_linf']:.6e}`; residual relative to dominant "
            f"`{balance['residual_relative_to_dominant']:.6e}`; solver-floor reference "
            f"`{balance['solver_floor_reference']:.6e}`. The nontrivial-term floor "
            f"`{balance['nontrivial_floor_reference']:.6e}` is a degeneracy tripwire, "
            "not a calibrated non-degeneracy scale."
        )
        lines.append(
            "Edit-drift pins: source max_abs_diff "
            f"`{balance['source_recompute_max_abs_diff']:.6e}`; flux-divergence "
            f"max_abs_diff `{balance['flux_recompute_max_abs_diff']:.6e}`."
        )
        lines.append(
            "For the source, no alternate discretization is meaningful because it is an "
            "exact reduced source. Correctness is supported by reciprocity with the "
            "forward wall-to-matter force, the chunk-1a dual-engine MMS, and the "
            "fidelity audit; the row above only pins wiring/edit drift."
        )
    lines.append("")
    lines.append("Counted balance gates:")
    lines.append(_table(["gate", "status"], _gate_rows(balance["checks"] if balance else {})))
    lines.append("")
    lines.append("Balance identity checks, reported but not counted as physics gates:")
    lines.append(_table(["gate", "status"], _gate_rows(balance.get("identity_checks", {}) if balance else {})))
    lines.append("")
    lines.append("## Closed Grid-Convergence Ladder")
    lines.append("")
    lines.append(
        _table(
            [
                "level",
                "grid",
                "dof",
                "wall_clock_seconds",
                "converged",
                "final_residual_linf",
                "final_tolerance",
                "r0_min",
                "r0_max",
                "message",
            ],
            convergence["level_rows"],
        )
    )
    lines.append("")
    lines.append("Raw-field self-difference, including R0:")
    lines.append("")
    lines.append(
        _table(
            [
                "coarse_grid",
                "fine_grid",
                "raw_field_relative_l2_change",
                "r0_relative_l2_change",
                "raw_field_linf_change",
                "r0_linf_change",
                "observed_order_from_raw_l2_change",
                "observed_order_from_r0_l2_change",
            ],
            convergence["self_convergence_rows"],
        )
    )
    lines.append("")
    lines.append("Independent flux discretization check:")
    lines.append("")
    lines.append(
        "The non-conservative reconstruction uses center gradients only: "
        "`g = wall_center_gradient(R0)`, `q = T_w(R0, w_center) * g`, "
        "`-wall_center_gradient(q)`. It is compared with the conservative face-flux "
        "operator on interior wall cells only; the mouth one-sided stencil and "
        "zero-traction exit cell are excluded because those closures intentionally differ."
    )
    lines.append("")
    lines.append(
        _table(
            [
                "grid",
                "spacing",
                "interior_cells_compared",
                "interior_max_abs_diff",
                "observed_order",
                "h2_scaled_tolerance",
                "decreases_or_at_floor",
                "finest_within_tolerance",
            ],
            convergence["flux_discretization_rows"],
        )
    )
    lines.append("")
    lines.append("R0-aware surrogate observables:")
    lines.append("")
    lines.append(
        _table(
            [
                "observable",
                "finest_grid",
                "finest_value",
                "last_observed_order",
                "richardson_estimate",
                "finest_error_estimate",
            ],
            convergence["observable_summary"],
        )
    )
    lines.append("")
    lines.append(f"Laptop-limit stop reason: {convergence['stop_reason']}.")
    lines.append("")
    lines.append("Convergence gates:")
    lines.append(_table(["gate", "status"], _gate_rows(convergence["pass_checks"])))
    lines.append("")
    lines.append("Convergence identity checks, reported but not counted as physics gates:")
    lines.append(_table(["gate", "status"], _gate_rows(convergence.get("identity_checks", {}))))
    lines.append("")
    lines.append("## Closed Conservation Diagnostics")
    lines.append("")
    lines.append(
        "Number/charge and the independent Gauss closure are evaluated on the solved "
        "matter/gauge fields. The energy-density confinement term uses frozen geometry "
        "because `conservation_diagnostics._energy_density` calls "
        "`confinement_potential_torch` without `radius=`; the closed-vs-frozen "
        "comparison is apples-to-apples for that channel because both branches use the "
        "same frozen confinement term."
    )
    lines.append("")
    if "baseline_factor_note" in conservation:
        lines.append(conservation["baseline_factor_note"])
        lines.append(
            f"Closed/frozen metric threshold factor: `{conservation['baseline_factor']:.1f}`."
        )
        lines.append(conservation["null_floor_label_note"])
        lines.append("")
    if "frozen_solve" in conservation:
        frozen = conservation["frozen_solve"]
        lines.append(
            "Frozen baseline solve: "
            f"grid `{frozen['grid']}`, converged={frozen['converged']}, "
            f"final_residual_linf={frozen['final_residual_linf']:.6e}, "
            f"wall_clock={frozen['wall_clock_seconds']:.6e}s."
        )
        lines.append("")
    lines.append(
        _table(
            [
                "metric",
                "activity",
                "closed_r0_solved",
                "frozen_geometry_baseline",
                "threshold",
                "basis",
                "status",
            ],
            _conservation_metric_rows(conservation),
        )
    )
    lines.append("")
    lines.append("Independent center-gradient Gauss closure rows:")
    lines.append("")
    lines.append(
        _table(
            [
                "surface",
                "reconstruction",
                "relative_residual",
                "absolute_residual",
                "surface_flux",
                "enclosed_mu0_charge",
            ],
            _gauss_closure_rows(conservation),
        )
    )
    lines.append("")
    lines.append("Conservation gates:")
    lines.append(_table(["gate", "status"], _gate_rows(conservation.get("pass_checks", {}))))
    lines.append("")
    lines.append("Identity checks, reported but not counted as physics gates:")
    lines.append(_table(["gate", "status"], _gate_rows(conservation.get("identity_checks", {}))))
    lines.append("")
    lines.append("## Counted Gates")
    lines.append("")
    lines.append(_table(["gate", "status"], _gate_rows(results["pass_checks"])))
    lines.append("")
    lines.append("## Identity / Not Physics Gates")
    lines.append("")
    lines.append(_table(["gate", "status"], _gate_rows(results.get("identity_checks", {}))))
    lines.append("")
    lines.append("## Stop Conditions")
    lines.append("")
    lines.append(
        _table(
            ["condition", "status"],
            [
                {
                    "condition": key,
                    "status": "TRIPPED" if value else "not tripped",
                }
                for key, value in results["stop_conditions"].items()
            ],
        )
    )
    lines.append("")
    lines.append("## Reproduction")
    lines.append("")
    lines.append("```bash")
    lines.append(
        "timeout 600 env PYTHONPATH=software/stage1_solver/src "
        "python -m stage1_solver.patha_closed_validation"
    )
    lines.append(
        "timeout 600 env PYTHONPATH=software/stage1_solver/src "
        "pytest software/stage1_solver/tests -q"
    )
    lines.append("```")
    lines.append("")
    lines.append(f"Machine-readable table: `{results['machine_readable_table']}`.")
    lines.append(
        f"Target-token scan: {'PASS' if results['target_token_scan']['passed'] else 'FAIL'} "
        f"(scanned {results['target_token_scan']['path_count']} Path-A 1c files)."
    )
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def main() -> int:
    config = PathAClosedValidationConfig()
    results = run_patha_closed_validation(config)
    report_path = write_patha_closed_validation_report(results, config.report_path)
    token_scan = target_token_scan(validation_scan_paths(report_path))
    results["target_token_scan"] = token_scan
    results["pass_checks"]["target_token_scan_clean"] = bool(token_scan["passed"])
    results["passed"] = all(results["pass_checks"].values())
    table_path = Path(results["machine_readable_table"])
    table_path.write_text(
        json.dumps(results, indent=2, sort_keys=True, default=_json_default),
        encoding="utf-8",
    )
    write_patha_closed_validation_report(results, config.report_path)
    print(f"wrote Path-A chunk 1c validation report: {report_path}")
    if not results["passed"]:
        print("Path-A chunk 1c validation gate failed")
        return 1
    print("Path-A chunk 1c validation gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
