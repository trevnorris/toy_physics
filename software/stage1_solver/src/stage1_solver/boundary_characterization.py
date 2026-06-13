"""Step 5 stationary boundary-control characterization.

This module studies whether finite-domain boundary choices contaminate the
stationary coupled-branch interior.  It does not introduce a wave problem or an
extraction map; the genuine outgoing-wave reflection coefficient belongs to the
later tangent study.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .backend import configure_backend
from .config import BranchSmokeConfig, HarnessConfig, TensorGridSpec
from .convergence import OBSERVABLE_LABELS, step4_preconditioner_config
from .coupled_branch import (
    _create_branch_grid,
    _matter_number_current,
    branch_boundary_conditions,
    confinement_potential_torch,
    coupled_branch_residual,
    run_branch_continuation,
    tensor_center_gradient_r,
    tensor_center_gradient_w,
    target_blind_surrogate_observables,
    unpack_coupled_fields,
)
from .grid import TensorProductGrid


@dataclass(frozen=True)
class InteriorWindow:
    r_min: float = 0.0
    r_max: float = 1.2
    w_min: float = 0.2
    w_max: float = 1.0

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class BoundaryCharacterizationConfig:
    name: str = "step5_boundary_characterization"
    spacing: float = 0.05
    fixed_w_max: float = 1.6
    truncation_r_max_values: tuple[float, ...] = (2.8, 3.2, 3.6)
    impedance_r_max: float = 3.2
    impedance_alpha_scales: tuple[float, ...] = (0.5, 1.0, 2.0)
    interior: InteriorWindow = field(default_factory=InteriorWindow)
    step4_reference_relative_floor: float = 3.477501e-4
    threshold_floor_multiplier: float = 4.0
    sponge_width: float = 0.4
    sponge_matter_strength: float = 100.0
    sponge_gauge_strength: float = 100.0
    sponge_power: int = 2

    @property
    def relative_error_threshold(self) -> float:
        return self.threshold_floor_multiplier * self.step4_reference_relative_floor

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "spacing": self.spacing,
            "fixed_w_max": self.fixed_w_max,
            "truncation_r_max_values": list(self.truncation_r_max_values),
            "impedance_r_max": self.impedance_r_max,
            "impedance_alpha_scales": list(self.impedance_alpha_scales),
            "interior": self.interior.to_dict(),
            "step4_reference_relative_floor": self.step4_reference_relative_floor,
            "threshold_floor_multiplier": self.threshold_floor_multiplier,
            "relative_error_threshold": self.relative_error_threshold,
            "sponge_width": self.sponge_width,
            "sponge_matter_strength": self.sponge_matter_strength,
            "sponge_gauge_strength": self.sponge_gauge_strength,
            "sponge_power": self.sponge_power,
        }


@dataclass(frozen=True)
class BoundarySolve:
    sweep: str
    setting: str
    variable_value: float
    config: HarnessConfig
    grid: TensorProductGrid
    state: torch.Tensor
    summary: dict[str, Any]
    interior_surrogates: dict[str, float]
    full_domain_surrogates: dict[str, float]


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _grid_count(extent: float, spacing: float) -> int:
    count = int(round(extent / spacing))
    if abs(count * spacing - extent) > 1.0e-12:
        raise ValueError(f"extent {extent} is not an integer multiple of spacing {spacing}")
    return count


def _branch_with_step5_preconditioner(
    branch: BranchSmokeConfig,
    study: BoundaryCharacterizationConfig,
    *,
    sponge_enabled: bool,
) -> BranchSmokeConfig:
    return replace(
        branch,
        name="step5_boundary_characterization_engineering_smoke",
        sponge_enabled=sponge_enabled,
        sponge_width=study.sponge_width if sponge_enabled else 0.0,
        sponge_matter_strength=study.sponge_matter_strength if sponge_enabled else 0.0,
        sponge_gauge_strength=study.sponge_gauge_strength if sponge_enabled else 0.0,
        sponge_power=study.sponge_power,
        newton=replace(branch.newton, preconditioner=step4_preconditioner_config()),
    )


def _truncation_branch(
    base: BranchSmokeConfig,
    study: BoundaryCharacterizationConfig,
    r_max: float,
) -> BranchSmokeConfig:
    nr = _grid_count(r_max, study.spacing)
    nw = _grid_count(study.fixed_w_max - base.w_min, study.spacing)
    return replace(
        base,
        r_max=r_max,
        w_max=study.fixed_w_max,
        solve_grid=(nr, nw),
    )


def _impedance_branch(
    base: BranchSmokeConfig,
    study: BoundaryCharacterizationConfig,
    scale: float,
) -> BranchSmokeConfig:
    nr = _grid_count(study.impedance_r_max, study.spacing)
    nw = _grid_count(study.fixed_w_max - base.w_min, study.spacing)
    return replace(
        base,
        r_max=study.impedance_r_max,
        w_max=study.fixed_w_max,
        solve_grid=(nr, nw),
        matter_exit_impedance_alpha=base.matter_exit_impedance_alpha * scale,
        a0_radial_impedance_alpha=base.a0_radial_impedance_alpha * scale,
        a0_exit_impedance_alpha=base.a0_exit_impedance_alpha * scale,
    )


def truncation_sweep_branches(
    base: BranchSmokeConfig,
    study: BoundaryCharacterizationConfig,
    *,
    sponge_enabled: bool = False,
) -> list[BranchSmokeConfig]:
    preconditioned = _branch_with_step5_preconditioner(
        base,
        study,
        sponge_enabled=sponge_enabled,
    )
    return [_truncation_branch(preconditioned, study, r_max) for r_max in study.truncation_r_max_values]


def impedance_sweep_branches(
    base: BranchSmokeConfig,
    study: BoundaryCharacterizationConfig,
    *,
    sponge_enabled: bool = False,
) -> list[BranchSmokeConfig]:
    preconditioned = _branch_with_step5_preconditioner(
        base,
        study,
        sponge_enabled=sponge_enabled,
    )
    return [_impedance_branch(preconditioned, study, scale) for scale in study.impedance_alpha_scales]


def _interior_indices(grid: TensorProductGrid, interior: InteriorWindow) -> tuple[torch.Tensor, torch.Tensor]:
    r_mask = (grid.r_centers >= interior.r_min) & (grid.r_centers <= interior.r_max)
    w_mask = (grid.w_centers >= interior.w_min) & (grid.w_centers <= interior.w_max)
    if int(torch.count_nonzero(r_mask).detach().cpu().item()) == 0:
        raise ValueError("interior window contains no radial cells")
    if int(torch.count_nonzero(w_mask).detach().cpu().item()) == 0:
        raise ValueError("interior window contains no axial cells")
    return torch.nonzero(r_mask, as_tuple=True)[0], torch.nonzero(w_mask, as_tuple=True)[0]


def _interior_tensor(values: torch.Tensor, r_idx: torch.Tensor, w_idx: torch.Tensor) -> torch.Tensor:
    return values.index_select(0, r_idx).index_select(1, w_idx)


def _interior_fields(
    state: torch.Tensor,
    grid: TensorProductGrid,
    interior: InteriorWindow,
) -> tuple[dict[str, torch.Tensor], torch.Tensor]:
    fields, _ = unpack_coupled_fields(state, grid, has_chemical_potential=True)
    r_idx, w_idx = _interior_indices(grid, interior)
    volumes = _interior_tensor(grid.cell_volumes, r_idx, w_idx)
    return {
        "psi_real": _interior_tensor(fields.psi_real, r_idx, w_idx),
        "psi_imag": _interior_tensor(fields.psi_imag, r_idx, w_idx),
        "a0": _interior_tensor(fields.a0, r_idx, w_idx),
        "ar": _interior_tensor(fields.ar, r_idx, w_idx),
        "aw": _interior_tensor(fields.aw, r_idx, w_idx),
    }, volumes


def _assert_same_interior_coordinates(
    left_grid: TensorProductGrid,
    right_grid: TensorProductGrid,
    interior: InteriorWindow,
) -> None:
    left_r, left_w = _interior_indices(left_grid, interior)
    right_r, right_w = _interior_indices(right_grid, interior)
    left_r_values = left_grid.r_centers.index_select(0, left_r).detach().cpu().numpy()
    right_r_values = right_grid.r_centers.index_select(0, right_r).detach().cpu().numpy()
    left_w_values = left_grid.w_centers.index_select(0, left_w).detach().cpu().numpy()
    right_w_values = right_grid.w_centers.index_select(0, right_w).detach().cpu().numpy()
    if left_r_values.shape != right_r_values.shape or not np.allclose(
        left_r_values,
        right_r_values,
        rtol=0.0,
        atol=1.0e-13,
    ):
        raise ValueError("interior radial cells do not match across sweep settings")
    if left_w_values.shape != right_w_values.shape or not np.allclose(
        left_w_values,
        right_w_values,
        rtol=0.0,
        atol=1.0e-13,
    ):
        raise ValueError("interior axial cells do not match across sweep settings")


def interior_solution_difference(
    state: torch.Tensor,
    grid: TensorProductGrid,
    reference_state: torch.Tensor,
    reference_grid: TensorProductGrid,
    interior: InteriorWindow,
) -> dict[str, float]:
    """Measure-consistent raw-field difference on a fixed interior window."""

    _assert_same_interior_coordinates(grid, reference_grid, interior)
    fields, volumes = _interior_fields(state, grid, interior)
    ref_fields, ref_volumes = _interior_fields(reference_state, reference_grid, interior)
    if not torch.allclose(volumes, ref_volumes, rtol=0.0, atol=1.0e-13):
        raise ValueError("interior control volumes do not match across sweep settings")

    diff_density = torch.zeros_like(volumes)
    ref_density = torch.zeros_like(volumes)
    linf = 0.0
    per_field_relative: dict[str, float] = {}
    for name, values in fields.items():
        diff = values - ref_fields[name]
        diff_density = diff_density + diff**2
        ref_density = ref_density + ref_fields[name] ** 2
        linf = max(linf, float(torch.max(torch.abs(diff)).detach().cpu().item()))
        field_l2 = torch.sqrt(torch.sum(diff**2 * volumes))
        field_ref = torch.sqrt(torch.sum(ref_fields[name] ** 2 * volumes))
        per_field_relative[f"{name}_relative_l2_change"] = float(
            (field_l2 / torch.clamp(field_ref, min=1.0e-300)).detach().cpu().item()
        )

    l2 = float(torch.sqrt(torch.sum(diff_density * volumes)).detach().cpu().item())
    reference = float(torch.sqrt(torch.sum(ref_density * volumes)).detach().cpu().item())
    return {
        "interior_l2_change": l2,
        "interior_linf_change": linf,
        "interior_signal_l2_reference": reference,
        "interior_relative_l2_change": l2 / max(reference, 1.0e-300),
        **per_field_relative,
    }


def interior_surrogate_observables(
    state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    *,
    eos_K: float,
    final_residual_linf: float | None,
    interior: InteriorWindow,
) -> dict[str, float]:
    """Interior-restricted counterparts of the step-4 raw-field surrogates."""

    fields, mu = unpack_coupled_fields(state, grid, has_chemical_potential=True)
    r_idx, w_idx = _interior_indices(grid, interior)
    volumes = _interior_tensor(grid.cell_volumes, r_idx, w_idx)
    density = fields.psi_real**2 + fields.psi_imag**2
    jr_number, jw_number = _matter_number_current(fields, grid, cfg)
    grad_real_r = tensor_center_gradient_r(fields.psi_real, grid)
    grad_imag_r = tensor_center_gradient_r(fields.psi_imag, grid)
    grad_real_w = tensor_center_gradient_w(fields.psi_real, grid)
    grad_imag_w = tensor_center_gradient_w(fields.psi_imag, grid)

    density_i = _interior_tensor(density, r_idx, w_idx)
    psi_real_i = _interior_tensor(fields.psi_real, r_idx, w_idx)
    psi_imag_i = _interior_tensor(fields.psi_imag, r_idx, w_idx)
    a0_i = _interior_tensor(fields.a0, r_idx, w_idx)
    ar_i = _interior_tensor(fields.ar, r_idx, w_idx)
    aw_i = _interior_tensor(fields.aw, r_idx, w_idx)
    jr_i = _interior_tensor(jr_number, r_idx, w_idx)
    jw_i = _interior_tensor(jw_number, r_idx, w_idx)
    gradient_density = (
        _interior_tensor(grad_real_r, r_idx, w_idx) ** 2
        + _interior_tensor(grad_imag_r, r_idx, w_idx) ** 2
        + _interior_tensor(grad_real_w, r_idx, w_idx) ** 2
        + _interior_tensor(grad_imag_w, r_idx, w_idx) ** 2
    )
    potential = _interior_tensor(confinement_potential_torch(grid, cfg), r_idx, w_idx)
    raw_field_norm_density = psi_real_i**2 + psi_imag_i**2 + a0_i**2 + ar_i**2 + aw_i**2
    energy_like_density = (
        (cfg.hbar**2 / (2.0 * cfg.particle_mass)) * gradient_density
        + potential * density_i
        + 0.25 * eos_K * density_i**5
        + 0.5 * (a0_i**2 + ar_i**2 + aw_i**2)
    )
    result = {
        "density_mass": float(torch.sum(density_i * volumes).detach().cpu().item()),
        "peak_density": float(torch.max(density_i).detach().cpu().item()),
        "min_density": float(torch.min(density_i).detach().cpu().item()),
        "raw_field_l2_norm": float(
            torch.sqrt(torch.sum(raw_field_norm_density * volumes)).detach().cpu().item()
        ),
        "scalar_gauge_l2": float(torch.sqrt(torch.sum(a0_i**2 * volumes)).detach().cpu().item()),
        "spatial_gauge_l2": float(
            torch.sqrt(torch.sum((ar_i**2 + aw_i**2) * volumes)).detach().cpu().item()
        ),
        "spatial_current_l2": float(
            torch.sqrt(torch.sum((jr_i**2 + jw_i**2) * volumes)).detach().cpu().item()
        ),
        "field_energy_like_integral": float(torch.sum(energy_like_density * volumes).detach().cpu().item()),
        "chemical_potential": float(mu.detach().cpu().item()) if mu is not None else float("nan"),
    }
    if final_residual_linf is not None:
        result["final_residual_linf"] = float(final_residual_linf)
    return result


def surrogate_change_rows(
    values: dict[str, float],
    reference_values: dict[str, float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, reference in reference_values.items():
        value = values[name]
        change = abs(value - reference)
        reference_scale = max(abs(reference), 1.0e-300)
        rows.append(
            {
                "observable": name,
                "label": OBSERVABLE_LABELS.get(name, name),
                "value": value,
                "reference_value": reference,
                "absolute_change": change,
                "relative_change": change / reference_scale,
            }
        )
    return rows


def _solve_setting(
    *,
    parent_config: HarnessConfig,
    branch: BranchSmokeConfig,
    study: BoundaryCharacterizationConfig,
    sweep: str,
    setting: str,
    variable_value: float,
    dtype: torch.dtype,
) -> BoundarySolve:
    run_config = replace(
        parent_config,
        branch=branch,
        run_root=str(Path(parent_config.run_root) / sweep),
    )
    grid = _create_branch_grid(
        branch,
        branch.solve_grid,
        dtype=dtype,
        device=run_config.backend.device,
    )
    state, summary = run_branch_continuation(
        run_config,
        grid,
        grid_name=setting,
    )
    final_K = branch.continuation_K_values[-1]
    boundaries = branch_boundary_conditions(branch)
    final_residual = coupled_branch_residual(
        state,
        grid,
        branch,
        eos_K=final_K,
        boundaries=boundaries,
    )
    final_residual_linf = float(torch.max(torch.abs(final_residual)).detach().cpu().item())
    interior_values = interior_surrogate_observables(
        state,
        grid,
        branch,
        eos_K=final_K,
        final_residual_linf=final_residual_linf,
        interior=study.interior,
    )
    full_values = target_blind_surrogate_observables(
        state,
        grid,
        branch,
        eos_K=final_K,
        final_residual_linf=final_residual_linf,
    )
    return BoundarySolve(
        sweep=sweep,
        setting=setting,
        variable_value=variable_value,
        config=run_config,
        grid=grid,
        state=state.detach(),
        summary=summary,
        interior_surrogates=interior_values,
        full_domain_surrogates=full_values,
    )


def _solve_rows(solves: list[BoundarySolve]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for solve in solves:
        rows.append(
            {
                "sweep": solve.sweep,
                "setting": solve.setting,
                "variable_value": solve.variable_value,
                "nr": solve.grid.spec.nr,
                "nw": solve.grid.spec.nw,
                "r_max": solve.grid.spec.r_max,
                "w_max": solve.grid.spec.w_max,
                "dof": solve.summary["dof"],
                "converged": solve.summary["converged"],
                "message": solve.summary["message"],
                "wall_clock_seconds": solve.summary["wall_clock_seconds"],
                "final_residual_linf": solve.summary["final_residual_linf"],
                "manifest": solve.summary["manifest"],
                "preconditioner": solve.summary["preconditioner"],
                "boundaries": solve.summary["boundaries"],
            }
        )
    return rows


def _comparison_rows(
    solves: list[BoundarySolve],
    reference: BoundarySolve,
    study: BoundaryCharacterizationConfig,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for solve in solves:
        solution_metric = interior_solution_difference(
            solve.state,
            solve.grid,
            reference.state,
            reference.grid,
            study.interior,
        )
        observable_rows = surrogate_change_rows(
            solve.interior_surrogates,
            reference.interior_surrogates,
        )
        max_observable_relative = max(
            row["relative_change"]
            for row in observable_rows
            if np.isfinite(row["relative_change"])
        )
        rows.append(
            {
                "setting": solve.setting,
                "reference_setting": reference.setting,
                "variable_value": solve.variable_value,
                **solution_metric,
                "relative_to_step4_reference_floor": (
                    solution_metric["interior_relative_l2_change"]
                    / study.step4_reference_relative_floor
                ),
                "max_interior_surrogate_relative_change": max_observable_relative,
                "interior_surrogate_changes": observable_rows,
                "converged": solve.summary["converged"],
                "message": solve.summary["message"],
            }
        )
    return rows


def _max_nonreference_relative(rows: list[dict[str, Any]]) -> float:
    return max(
        (
            row["interior_relative_l2_change"]
            for row in rows
            if row["setting"] != row["reference_setting"]
        ),
        default=0.0,
    )


def _run_sweep_pair(
    *,
    config: HarnessConfig,
    study: BoundaryCharacterizationConfig,
    dtype: torch.dtype,
    label: str,
    sponge_enabled: bool,
) -> dict[str, Any]:
    truncation_branches = truncation_sweep_branches(
        config.branch,
        study,
        sponge_enabled=sponge_enabled,
    )
    impedance_branches = impedance_sweep_branches(
        config.branch,
        study,
        sponge_enabled=sponge_enabled,
    )

    truncation_solves = [
        _solve_setting(
            parent_config=config,
            branch=branch,
            study=study,
            sweep=f"{label}_truncation_placement",
            setting=f"{label}_rmax_{branch.r_max:.1f}_nr_{branch.solve_grid[0]}_nw_{branch.solve_grid[1]}",
            variable_value=branch.r_max,
            dtype=dtype,
        )
        for branch in truncation_branches
    ]
    impedance_solves = [
        _solve_setting(
            parent_config=config,
            branch=branch,
            study=study,
            sweep=f"{label}_impedance_coefficients",
            setting=f"{label}_alpha_scale_{scale:.1f}_nr_{branch.solve_grid[0]}_nw_{branch.solve_grid[1]}",
            variable_value=scale,
            dtype=dtype,
        )
        for branch, scale in zip(impedance_branches, study.impedance_alpha_scales)
    ]

    truncation_reference = max(truncation_solves, key=lambda solve: solve.variable_value)
    baseline_scale = min(
        impedance_solves,
        key=lambda solve: abs(solve.variable_value - 1.0),
    )
    truncation_comparisons = _comparison_rows(truncation_solves, truncation_reference, study)
    impedance_comparisons = _comparison_rows(impedance_solves, baseline_scale, study)
    max_truncation_relative = _max_nonreference_relative(truncation_comparisons)
    max_impedance_relative = _max_nonreference_relative(impedance_comparisons)
    max_boundary_relative = max(max_truncation_relative, max_impedance_relative)
    all_converged = all(
        solve.summary["converged"] for solve in [*truncation_solves, *impedance_solves]
    )
    passed = all_converged and max_boundary_relative <= study.relative_error_threshold
    metric = {
        "max_relative_l2": max_boundary_relative,
        "reference": "fixed-window interior raw-field L2 signal magnitude",
        "step4_reference_relative_floor": study.step4_reference_relative_floor,
        "relative_error_threshold": study.relative_error_threshold,
        "ratio_to_threshold": max_boundary_relative / study.relative_error_threshold,
        "ratio_to_step4_reference_floor": (
            max_boundary_relative / study.step4_reference_relative_floor
        ),
    }
    return {
        "label": label,
        "sponge_enabled": sponge_enabled,
        "truncation_solves": truncation_solves,
        "impedance_solves": impedance_solves,
        "truncation_reference": truncation_reference.setting,
        "impedance_reference": baseline_scale.setting,
        "truncation_comparisons": truncation_comparisons,
        "impedance_comparisons": impedance_comparisons,
        "max_truncation_relative": max_truncation_relative,
        "max_impedance_relative": max_impedance_relative,
        "metric": metric,
        "all_converged": all_converged,
        "passed": passed,
    }


def _sweep_pair_payload(pair: dict[str, Any]) -> dict[str, Any]:
    solves = [*pair["truncation_solves"], *pair["impedance_solves"]]
    return {
        "label": pair["label"],
        "sponge_enabled": pair["sponge_enabled"],
        "solve_rows": _solve_rows(solves),
        "truncation": {
            "settings": [solve.setting for solve in pair["truncation_solves"]],
            "reference": pair["truncation_reference"],
            "comparisons": pair["truncation_comparisons"],
            "max_nonreference_relative_l2": pair["max_truncation_relative"],
        },
        "impedance": {
            "settings": [solve.setting for solve in pair["impedance_solves"]],
            "reference": pair["impedance_reference"],
            "comparisons": pair["impedance_comparisons"],
            "max_nonreference_relative_l2": pair["max_impedance_relative"],
        },
        "boundary_error_metric": pair["metric"],
        "all_converged": pair["all_converged"],
        "passed": pair["passed"],
    }


def run_step5(
    config: HarnessConfig | None = None,
    study: BoundaryCharacterizationConfig | None = None,
) -> dict[str, Any]:
    if config is None:
        config = HarnessConfig(
            run_root="software/stage1_solver/runs/step5_boundary_characterization",
            report_path="software/stage1_solver/reports/step5_boundary_characterization.md",
        )
    if study is None:
        study = BoundaryCharacterizationConfig()
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(config.backend)

    robin_pair = _run_sweep_pair(
        config=config,
        study=study,
        dtype=dtype,
        label="robin_only",
        sponge_enabled=False,
    )
    sponge_pair = None
    if not robin_pair["passed"]:
        sponge_pair = _run_sweep_pair(
            config=config,
            study=study,
            dtype=dtype,
            label="sponge",
            sponge_enabled=True,
        )
    final_pair = sponge_pair if sponge_pair is not None else robin_pair
    robin_adequate = robin_pair["passed"]
    sponge_required = not robin_adequate
    sponge_effective = sponge_pair["passed"] if sponge_pair is not None else False
    passed = robin_adequate or sponge_effective
    solve_rows = _solve_rows(
        [
            *robin_pair["truncation_solves"],
            *robin_pair["impedance_solves"],
            *([] if sponge_pair is None else sponge_pair["truncation_solves"]),
            *([] if sponge_pair is None else sponge_pair["impedance_solves"]),
        ]
    )

    table_path = Path(config.run_root) / study.name / "boundary_characterization_table.json"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "study": study.to_dict(),
        "method": {
            "discipline": (
                "engineering-smoke placeholders only; stationary branch; raw-field "
                "target-blind interior diagnostics"
            ),
            "fixed_resolution": (
                f"uniform cell spacing {study.spacing}; truncation changes add or remove "
                "outer radial cells only"
            ),
            "placement_policy": (
                "The acceptance sweep uses radial placements at or beyond 2.8. "
                "Nearer stress placements down to 2.4 showed residual boundary-zone "
                "contamination and are not treated as adequate placements."
            ),
            "interior_definition": study.interior.to_dict(),
            "difference_norm": (
                "sqrt(sum_fields integral_interior (u - u_ref)^2 dV) with "
                "dV=4*pi*r^2 dr dw, divided by the reference interior raw-field L2"
            ),
            "reference_scale": (
                "interior raw-field L2 magnitude; threshold is a fixed multiple of "
                "the step-4 raw-field relative self-difference"
            ),
            "truncation_reference": final_pair["truncation_reference"],
            "impedance_reference": final_pair["impedance_reference"],
            "preconditioner": asdict(step4_preconditioner_config()),
            "sponge": {
                "used": sponge_pair is not None,
                "reason": (
                    "Robin-only exceeded the stationary interior-contamination criterion."
                    if sponge_pair is not None
                    else "Robin-only met the stationary interior-contamination criterion."
                ),
                "parameters": {
                    "width": study.sponge_width,
                    "matter_strength": study.sponge_matter_strength,
                    "gauge_strength": study.sponge_gauge_strength,
                    "power": study.sponge_power,
                }
                if sponge_pair is not None
                else None,
            },
            "wave_reflection_deferral": (
                "The outgoing-wave reflection coefficient is deferred to the step-8 "
                "linearized tangent because the stationary isotropic branch has no "
                "propagating waves to reflect."
            ),
        },
        "solve_rows": solve_rows,
        "robin_only": _sweep_pair_payload(robin_pair),
        "sponge": _sweep_pair_payload(sponge_pair) if sponge_pair is not None else None,
        "truncation": {
            "settings": [solve.setting for solve in final_pair["truncation_solves"]],
            "comparisons": final_pair["truncation_comparisons"],
            "max_nonreference_relative_l2": final_pair["max_truncation_relative"],
        },
        "impedance": {
            "settings": [solve.setting for solve in final_pair["impedance_solves"]],
            "comparisons": final_pair["impedance_comparisons"],
            "max_nonreference_relative_l2": final_pair["max_impedance_relative"],
        },
        "boundary_error_metric": final_pair["metric"],
        "verdict": {
            "robin_impedance_adequate": robin_adequate,
            "sponge_required": sponge_required,
            "sponge_effective": sponge_effective,
            "statement": (
                "Robin impedance exits are adequate for the stationary branch interior study; "
                "no sponge was added."
                if robin_adequate
                else (
                    "Robin-only boundary treatment exceeded the stationary criterion; "
                    "the smooth sponge layer is required and brought the measured "
                    "interior boundary error below threshold."
                    if sponge_effective
                    else (
                        "Robin-only boundary treatment exceeded the stationary criterion, "
                        "and the configured sponge did not bring the measured interior "
                        "boundary error below threshold."
                    )
                )
            ),
        },
        "pass_checks": {
            "minimum_three_truncation_settings": len(final_pair["truncation_solves"]) >= 3,
            "minimum_three_impedance_settings": len(final_pair["impedance_solves"]) >= 3,
            "all_final_solves_converged": final_pair["all_converged"],
            "uses_colored_sparse_jacobian_lu": (
                all(
                    solve.config.branch.newton.preconditioner.type == "colored_sparse_jacobian_lu"
                    for solve in [
                        *final_pair["truncation_solves"],
                        *final_pair["impedance_solves"],
                    ]
                )
            ),
            "boundary_error_below_threshold": (
                final_pair["metric"]["max_relative_l2"] <= study.relative_error_threshold
            ),
            "robin_adequate_or_sponge_effective": passed,
        },
        "passed": passed,
        "machine_readable_table": str(table_path),
    }
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


def _comparison_table_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "setting": row["setting"],
            "reference_setting": row["reference_setting"],
            "variable_value": row["variable_value"],
            "interior_relative_l2_change": row["interior_relative_l2_change"],
            "relative_to_step4_reference_floor": row["relative_to_step4_reference_floor"],
            "max_interior_surrogate_relative_change": row[
                "max_interior_surrogate_relative_change"
            ],
            "interior_linf_change": row["interior_linf_change"],
            "converged": row["converged"],
        }
        for row in rows
    ]


def write_step5_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    study = results["study"]
    method = results["method"]
    metric = results["boundary_error_metric"]
    verdict = results["verdict"]

    lines: list[str] = []
    lines.append("# Step 5 Boundary-Control Characterization")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append("")
    lines.append(
        "**Discipline:** engineering smoke, placeholder parameters, stationary branch, "
        "target-blind. This report uses raw-field interior diagnostics only."
    )
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append(f"Truncation sweep radial extents: `{study['truncation_r_max_values']}`.")
    lines.append(
        f"Impedance coefficient scales: `{study['impedance_alpha_scales']}` "
        f"on fixed `r_max={study['impedance_r_max']}`."
    )
    lines.append(f"Fixed resolution: {method['fixed_resolution']}.")
    lines.append(method["placement_policy"])
    lines.append(f"Interior window: `{method['interior_definition']}`.")
    lines.append(f"Difference norm: {method['difference_norm']}.")
    lines.append(f"Reference scale: {method['reference_scale']}.")
    lines.append(
        "Threshold: "
        f"{_fmt(study['relative_error_threshold'])} relative to the interior signal "
        f"({study['threshold_floor_multiplier']} x the step-4 reference floor "
        f"{_fmt(study['step4_reference_relative_floor'])})."
    )
    lines.append(
        "Pass criterion: the interior raw-field L2 relative change is the governing "
        "boundary-characterization pass criterion; per-surrogate-observable relative "
        "changes are diagnostic/informational only and are not thresholded."
    )
    lines.append("")
    lines.append("Preconditioner:")
    lines.append("```yaml")
    for key, value in method["preconditioner"].items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Solve Rows")
    lines.append("")
    lines.append(
        _table(
            [
                "sweep",
                "setting",
                "variable_value",
                "nr",
                "nw",
                "r_max",
                "w_max",
                "dof",
                "wall_clock_seconds",
                "final_residual_linf",
                "converged",
                "message",
            ],
            results["solve_rows"],
        )
    )
    lines.append("")
    if results["sponge"] is not None:
        lines.append("## Robin-Only Diagnostic")
        lines.append("")
        lines.append(
            "Robin-only exceeded the criterion, so the same sweeps were rerun with the "
            "smooth sponge layer enabled. The failed diagnostic is retained here."
        )
        lines.append("")
        lines.append("Robin-only truncation:")
        lines.append(
            _table(
                [
                    "setting",
                    "reference_setting",
                    "variable_value",
                    "interior_relative_l2_change",
                    "relative_to_step4_reference_floor",
                    "max_interior_surrogate_relative_change",
                    "interior_linf_change",
                    "converged",
                ],
                _comparison_table_rows(results["robin_only"]["truncation"]["comparisons"]),
            )
        )
        lines.append("")
        lines.append("Robin-only impedance:")
        lines.append(
            _table(
                [
                    "setting",
                    "reference_setting",
                    "variable_value",
                    "interior_relative_l2_change",
                    "relative_to_step4_reference_floor",
                    "max_interior_surrogate_relative_change",
                    "interior_linf_change",
                    "converged",
                ],
                _comparison_table_rows(results["robin_only"]["impedance"]["comparisons"]),
            )
        )
        lines.append("")
        robin_metric = results["robin_only"]["boundary_error_metric"]
        lines.append(
            "Robin-only achieved max interior relative L2 boundary error: "
            f"{_fmt(robin_metric['max_relative_l2'])}; threshold ratio "
            f"{_fmt(robin_metric['ratio_to_threshold'])}."
        )
        lines.append(
            "Physical note: the Robin-only impedance sensitivity is driven by the "
            "long-range A0 scalar-potential tail: the matter field psi is tightly "
            "confined and decays fast, but the A0 Poisson-type potential has a "
            "slowly-decaying tail that reaches the exit, so at fixed truncation the "
            "interior remains sensitive to how the Robin impedance terminates that "
            "tail, whereas moving the truncation outward only lowers the already-small "
            "tail amplitude."
        )
        lines.append("")
    lines.append("## Truncation-Placement Sweep")
    lines.append("")
    lines.append(
        _table(
            [
                "setting",
                "reference_setting",
                "variable_value",
                "interior_relative_l2_change",
                "relative_to_step4_reference_floor",
                "max_interior_surrogate_relative_change",
                "interior_linf_change",
                "converged",
            ],
            _comparison_table_rows(results["truncation"]["comparisons"]),
        )
    )
    lines.append("")
    lines.append(
        "Max non-reference truncation relative L2 change: "
        f"{_fmt(results['truncation']['max_nonreference_relative_l2'])}."
    )
    lines.append("")
    lines.append("## Impedance-Coefficient Sweep")
    lines.append("")
    lines.append(
        _table(
            [
                "setting",
                "reference_setting",
                "variable_value",
                "interior_relative_l2_change",
                "relative_to_step4_reference_floor",
                "max_interior_surrogate_relative_change",
                "interior_linf_change",
                "converged",
            ],
            _comparison_table_rows(results["impedance"]["comparisons"]),
        )
    )
    lines.append("")
    lines.append(
        "Max non-reference impedance relative L2 change: "
        f"{_fmt(results['impedance']['max_nonreference_relative_l2'])}."
    )
    lines.append("")
    lines.append("## Boundary-Error Metric")
    lines.append("")
    lines.append(
        "Achieved max interior relative L2 boundary error: "
        f"{_fmt(metric['max_relative_l2'])}; threshold ratio "
        f"{_fmt(metric['ratio_to_threshold'])}; step-4-floor ratio "
        f"{_fmt(metric['ratio_to_step4_reference_floor'])}."
    )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(verdict["statement"])
    lines.append(
        "Sponge/absorber used: "
        f"{_fmt(method['sponge']['used'])}; parameters: {_fmt(method['sponge']['parameters'])}."
    )
    lines.append("")
    lines.append("## Deferral Note")
    lines.append("")
    lines.append(method["wave_reflection_deferral"])
    lines.append("")
    lines.append("## Machine-Readable Output")
    lines.append("")
    lines.append(f"Boundary characterization table: `{results['machine_readable_table']}`")
    lines.append("")
    lines.append("Manifests:")
    for row in results["solve_rows"]:
        lines.append(f"- {row['setting']}: `{row['manifest']}`")
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path
