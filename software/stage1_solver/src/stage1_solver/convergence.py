"""Step 4 target-blind grid self-convergence study for the coupled branch."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .backend import configure_backend
from .config import HarnessConfig, PreconditionerConfig
from .coupled_branch import (
    _create_branch_grid,
    branch_boundary_conditions,
    coupled_branch_residual,
    resample_branch_state,
    run_branch_continuation,
    target_blind_surrogate_observables,
    unpack_coupled_fields,
)
from .grid import TensorProductGrid


OBSERVABLE_LABELS = {
    "density_mass": "density mass integral",
    "peak_density": "peak density",
    "min_density": "minimum density",
    "raw_field_l2_norm": "raw field L2 norm",
    "scalar_gauge_l2": "A0 L2 norm",
    "spatial_gauge_l2": "spatial gauge L2 norm",
    "spatial_current_l2": "spatial current L2 norm",
    "field_energy_like_integral": "field-energy-like integral",
    "chemical_potential": "chemical potential",
    "final_residual_linf": "final residual Linf",
}


def step4_preconditioner_config() -> PreconditionerConfig:
    return PreconditionerConfig(
        type="colored_sparse_jacobian_lu",
        side="left",
        rebuild_policy="once_per_continuation",
        stencil_radius=3,
        color_separation=7,
        factorization="splu",
        diagonal_shift=0.0,
        drop_tolerance=0.0,
        fill_factor=10.0,
        permutation="COLAMD",
    )


def with_step4_preconditioner(config: HarnessConfig) -> HarnessConfig:
    branch = config.branch
    return replace(
        config,
        branch=replace(
            branch,
            newton=replace(branch.newton, preconditioner=step4_preconditioner_config()),
        ),
    )


def validate_refinement_ladder(levels: tuple[tuple[int, int], ...], ratio: int) -> None:
    if len(levels) < 3:
        raise ValueError("Step 4 requires at least three refinement levels")
    if ratio <= 1:
        raise ValueError("refinement ratio must exceed one")
    for previous, current in zip(levels, levels[1:]):
        if current[0] != ratio * previous[0] or current[1] != ratio * previous[1]:
            raise ValueError(
                "Step 4 requires a fixed refinement ratio in both dimensions; "
                f"got {previous} -> {current} with ratio {ratio}"
            )


def _volume_restrict_2d(
    fine_values: torch.Tensor,
    fine_grid: TensorProductGrid,
    coarse_grid: TensorProductGrid,
) -> torch.Tensor:
    r_ratio = fine_grid.spec.nr // coarse_grid.spec.nr
    w_ratio = fine_grid.spec.nw // coarse_grid.spec.nw
    if (
        fine_grid.spec.nr != coarse_grid.spec.nr * r_ratio
        or fine_grid.spec.nw != coarse_grid.spec.nw * w_ratio
        or r_ratio < 1
        or w_ratio < 1
    ):
        raise ValueError("fine grid must be an integer refinement of coarse grid")
    blocks = fine_values.reshape(coarse_grid.spec.nr, r_ratio, coarse_grid.spec.nw, w_ratio)
    volume_blocks = fine_grid.cell_volumes.reshape(
        coarse_grid.spec.nr,
        r_ratio,
        coarse_grid.spec.nw,
        w_ratio,
    )
    restricted = torch.sum(blocks * volume_blocks, dim=(1, 3)) / torch.sum(
        volume_blocks,
        dim=(1, 3),
    )
    return restricted


def conservative_restrict_fields(
    fine_state: torch.Tensor,
    fine_grid: TensorProductGrid,
    coarse_grid: TensorProductGrid,
) -> dict[str, torch.Tensor]:
    fields, _ = unpack_coupled_fields(fine_state, fine_grid, has_chemical_potential=True)
    return {
        "psi_real": _volume_restrict_2d(fields.psi_real, fine_grid, coarse_grid),
        "psi_imag": _volume_restrict_2d(fields.psi_imag, fine_grid, coarse_grid),
        "a0": _volume_restrict_2d(fields.a0, fine_grid, coarse_grid),
        "ar": _volume_restrict_2d(fields.ar, fine_grid, coarse_grid),
        "aw": _volume_restrict_2d(fields.aw, fine_grid, coarse_grid),
    }


def raw_field_self_difference(
    coarse_state: torch.Tensor,
    coarse_grid: TensorProductGrid,
    fine_state: torch.Tensor,
    fine_grid: TensorProductGrid,
) -> dict[str, float]:
    coarse_fields, _ = unpack_coupled_fields(
        coarse_state,
        coarse_grid,
        has_chemical_potential=True,
    )
    restricted = conservative_restrict_fields(fine_state, fine_grid, coarse_grid)
    coarse_by_name = {
        "psi_real": coarse_fields.psi_real,
        "psi_imag": coarse_fields.psi_imag,
        "a0": coarse_fields.a0,
        "ar": coarse_fields.ar,
        "aw": coarse_fields.aw,
    }
    diff_density = torch.zeros_like(coarse_grid.cell_volumes)
    norm_density = torch.zeros_like(coarse_grid.cell_volumes)
    linf = 0.0
    for name, coarse_values in coarse_by_name.items():
        diff = coarse_values - restricted[name]
        diff_density = diff_density + diff**2
        norm_density = norm_density + coarse_values**2
        linf = max(linf, float(torch.max(torch.abs(diff)).detach().cpu().item()))
    l2 = float(torch.sqrt(torch.sum(diff_density * coarse_grid.cell_volumes)).detach().cpu().item())
    reference = float(
        torch.sqrt(torch.sum(norm_density * coarse_grid.cell_volumes)).detach().cpu().item()
    )
    return {
        "raw_field_l2_change": l2,
        "raw_field_linf_change": linf,
        "raw_field_relative_l2_change": l2 / max(reference, 1.0e-300),
    }


def observed_order_from_three(
    coarse_value: float,
    middle_value: float,
    fine_value: float,
    ratio: int,
) -> float | None:
    first = abs(coarse_value - middle_value)
    second = abs(middle_value - fine_value)
    if first <= 0.0 or second <= 0.0:
        return None
    order = math.log(first / second) / math.log(float(ratio))
    if not math.isfinite(order):
        return None
    return order


def richardson_estimate(
    previous_value: float,
    finest_value: float,
    ratio: int,
    order: float | None,
) -> float | None:
    if order is None or order <= 0.0:
        return None
    denominator = float(ratio) ** order - 1.0
    if denominator <= 0.0 or not math.isfinite(denominator):
        return None
    estimate = finest_value + (finest_value - previous_value) / denominator
    return estimate if math.isfinite(estimate) else None


def _diagnosis(observable: str, verdict: str) -> str:
    if verdict == "null diagnostic":
        return "This raw channel is identically zero on the completed placeholder branch; no order is measured."
    if observable == "density_mass":
        return "Newton mass constraint; differences read the solver floor, not discretization."
    if observable in {"final_residual_linf", "driven_residual_linf"}:
        return "Newton/GMRES stopping floor; useful as the numerical floor diagnostic."
    if observable in {"peak_density", "min_density"}:
        return "Pointwise extrema are sensitive to cell-center placement and coarse throat resolution."
    if observable in {"spatial_gauge_l2", "spatial_current_l2", "scalar_gauge_l2"}:
        return "Coupled gauge/current response can show reduced order from open Robin boundaries."
    if observable == "field_energy_like_integral":
        return "Gradient-weighted integral; boundary and coarse-grid throat terms can reduce order."
    if verdict == "drifts":
        return "Successive differences do not shrink on the completed ladder."
    return "Raw-field integral on the shared coupled branch."


def _verdict(observable: str, order: float | None, expected_order: float) -> str:
    if observable in {"density_mass", "final_residual_linf", "driven_residual_linf"}:
        return "solver-floor diagnostic"
    if order is None:
        return "drifts"
    if order >= 0.8 * expected_order:
        return "expected-order convergence"
    if order > 0.25:
        return "reduced-order convergence"
    return "drifts"


def _summarize_observables(
    level_rows: list[dict[str, Any]],
    *,
    ratio: int,
    expected_order: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if not level_rows:
        return [], []
    observable_names = list(level_rows[0]["surrogate_values"].keys())
    summary_rows: list[dict[str, Any]] = []
    table_rows: list[dict[str, Any]] = []
    dofs = [row["dof"] for row in level_rows]
    for observable in observable_names:
        values = [float(row["surrogate_values"][observable]) for row in level_rows]
        null_diagnostic = max(abs(value) for value in values) <= 1.0e-14
        changes = [None]
        for left, right in zip(values, values[1:]):
            changes.append(abs(right - left))
        orders: list[float | None] = [None, None]
        for index in range(2, len(values)):
            orders.append(
                observed_order_from_three(
                    values[index - 2],
                    values[index - 1],
                    values[index],
                    ratio,
                )
            )
        last_order = next((order for order in reversed(orders) if order is not None), None)
        if null_diagnostic:
            continuum = values[-1]
            errors = [0.0 for _ in values]
            verdict = "null diagnostic"
        else:
            continuum = (
                richardson_estimate(values[-2], values[-1], ratio, last_order)
                if len(values) >= 2
                else None
            )
            errors = [abs(value - continuum) if continuum is not None else None for value in values]
            verdict = _verdict(observable, last_order, expected_order)
        summary_rows.append(
            {
                "observable": observable,
                "label": OBSERVABLE_LABELS.get(observable, observable),
                "finest_grid": level_rows[-1]["grid"],
                "finest_dof": dofs[-1],
                "finest_value": values[-1],
                "last_observed_order": last_order,
                "richardson_estimate": continuum,
                "finest_error_estimate": errors[-1],
                "verdict": verdict,
                "diagnosis": _diagnosis(observable, verdict),
            }
        )
        for index, row in enumerate(level_rows):
            table_rows.append(
                {
                    "observable": observable,
                    "label": OBSERVABLE_LABELS.get(observable, observable),
                    "level": row["level"],
                    "grid": row["grid"],
                    "dof": row["dof"],
                    "value": values[index],
                    "successive_change": changes[index],
                    "observed_order": orders[index],
                    "richardson_estimate": continuum,
                    "error_estimate": errors[index],
                    "verdict": verdict,
                    "diagnosis": _diagnosis(observable, verdict),
                }
            )
    return summary_rows, table_rows


def _self_convergence_rows(
    solved_levels: list[dict[str, Any]],
    *,
    ratio: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index in range(1, len(solved_levels)):
        coarse = solved_levels[index - 1]
        fine = solved_levels[index]
        diff = raw_field_self_difference(
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
                "observed_order_from_l2_change": None,
                "observed_order_from_relative_l2_change": None,
            }
        )
    for index in range(1, len(rows)):
        rows[index]["observed_order_from_l2_change"] = observed_order_from_three(
            rows[index - 1]["raw_field_l2_change"],
            0.0,
            rows[index]["raw_field_l2_change"],
            ratio,
        )
        rows[index]["observed_order_from_relative_l2_change"] = observed_order_from_three(
            rows[index - 1]["raw_field_relative_l2_change"],
            0.0,
            rows[index]["raw_field_relative_l2_change"],
            ratio,
        )
    return rows


def _noise_floor(
    level_rows: list[dict[str, Any]],
    self_rows: list[dict[str, Any]],
    *,
    requested_mass: float,
    noise_floor_ratio: float,
) -> dict[str, Any]:
    residual_floor = max(
        (float(row["surrogate_values"]["final_residual_linf"]) for row in level_rows),
        default=None,
    )
    mass_floor = max(
        (abs(float(row["surrogate_values"]["density_mass"]) - requested_mass) for row in level_rows),
        default=None,
    )
    field_floor_reached = False
    if len(self_rows) >= 2:
        previous = self_rows[-2]["raw_field_relative_l2_change"]
        current = self_rows[-1]["raw_field_relative_l2_change"]
        field_floor_reached = current >= previous / noise_floor_ratio
    return {
        "solver_residual_floor_linf": residual_floor,
        "mass_constraint_floor": mass_floor,
        "field_self_difference_floor_reached": field_floor_reached,
        "last_raw_field_relative_l2_change": (
            self_rows[-1]["raw_field_relative_l2_change"] if self_rows else None
        ),
        "preliminary_numerical_floor": max(
            value for value in (residual_floor, mass_floor) if value is not None
        )
        if residual_floor is not None or mass_floor is not None
        else None,
    }


def _resolution_sizing(
    observable_summary: list[dict[str, Any]],
    *,
    thresholds: tuple[float, ...],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for summary in observable_summary:
        order = summary["last_observed_order"]
        error = summary["finest_error_estimate"]
        is_sizable = summary["verdict"] in {
            "expected-order convergence",
            "reduced-order convergence",
        }
        if not is_sizable or order is None or order <= 0.0 or error is None or error <= 0.0:
            sizing = {f"dof_for_error_{threshold:.0e}": None for threshold in thresholds}
        else:
            sizing = {}
            for threshold in thresholds:
                if error <= threshold:
                    required = summary["finest_dof"]
                else:
                    required = summary["finest_dof"] * (error / threshold) ** (2.0 / order)
                sizing[f"dof_for_error_{threshold:.0e}"] = int(math.ceil(required))
        rows.append(
            {
                "observable": summary["observable"],
                "label": summary["label"],
                "last_observed_order": order,
                "finest_error_estimate": error,
                "finest_dof": summary["finest_dof"],
                **sizing,
                "verdict": summary["verdict"],
            }
        )
    return rows


def _strip_runtime_objects(level: dict[str, Any]) -> dict[str, Any]:
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
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def run_step4(config: HarnessConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = HarnessConfig(
            run_root="software/stage1_solver/runs/step4_convergence_study",
            report_path="software/stage1_solver/reports/step4_convergence_study.md",
        )
    config = with_step4_preconditioner(config)
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(config.backend)
    study = config.convergence
    validate_refinement_ladder(study.levels, study.refinement_ratio)

    previous_grid: TensorProductGrid | None = None
    previous_state: torch.Tensor | None = None
    solved_levels: list[dict[str, Any]] = []
    stop_reason = "completed configured ladder"

    for level_index, level in enumerate(study.levels):
        grid = _create_branch_grid(config.branch, level, dtype=dtype, device=config.backend.device)
        initial = (
            resample_branch_state(previous_state, previous_grid, grid, config.branch)
            if previous_state is not None and previous_grid is not None
            else None
        )
        state, summary = run_branch_continuation(
            config,
            grid,
            initial_state=initial,
            grid_name=f"convergence_l{level_index}_nr_{level[0]}_nw_{level[1]}",
        )
        final_K = config.branch.continuation_K_values[-1]
        boundaries = branch_boundary_conditions(config.branch)
        final_residual = coupled_branch_residual(
            state,
            grid,
            config.branch,
            eos_K=final_K,
            boundaries=boundaries,
        )
        final_residual_linf = float(torch.max(torch.abs(final_residual)).detach().cpu().item())
        surrogate_values = target_blind_surrogate_observables(
            state,
            grid,
            config.branch,
            eos_K=final_K,
            final_residual_linf=final_residual_linf,
        )
        gmres_counts = [
            count
            for stage in summary["stages"]
            for count in stage.get("gmres_iterations", [])
            if count is not None
        ]
        level_record = {
            "level": level_index,
            "grid": summary["grid"],
            "nr": summary["nr"],
            "nw": summary["nw"],
            "spacing": max(grid.dr, grid.dw),
            "dof": summary["dof"],
            "wall_clock_seconds": summary["wall_clock_seconds"],
            "peak_memory_mb": summary["peak_memory_mb"],
            "newton_iterations": sum(stage["iterations"] for stage in summary["stages"]),
            "final_residual_linf": final_residual_linf,
            "gmres_iterations": gmres_counts,
            "gmres_max": max(gmres_counts) if gmres_counts else None,
            "gmres_mean": float(np.mean(gmres_counts)) if gmres_counts else None,
            "converged": summary["converged"],
            "message": summary["message"],
            "manifest": summary["manifest"],
            "surrogate_values": surrogate_values,
            "state": state.detach(),
            "grid_object": grid,
        }
        solved_levels.append(level_record)
        if not summary["converged"]:
            stop_reason = summary["message"]
            break
        previous_grid = grid
        previous_state = state
        if summary["wall_clock_seconds"] > study.max_level_seconds:
            stop_reason = (
                f"stopped after {summary['grid']} because the level exceeded "
                f"{study.max_level_seconds:.1f}s"
            )
            break

    level_rows = [_strip_runtime_objects(level) for level in solved_levels]
    self_rows = _self_convergence_rows(solved_levels, ratio=study.refinement_ratio)
    observable_summary, observable_table = _summarize_observables(
        level_rows,
        ratio=study.refinement_ratio,
        expected_order=study.expected_order,
    )
    noise_floor = _noise_floor(
        level_rows,
        self_rows,
        requested_mass=config.branch.mass,
        noise_floor_ratio=study.noise_floor_ratio,
    )
    sizing = _resolution_sizing(
        observable_summary,
        thresholds=study.illustrative_error_levels,
    )
    if level_rows:
        last = level_rows[-1]
        next_level = (
            last["nr"] * study.refinement_ratio,
            last["nw"] * study.refinement_ratio,
        )
        next_dof = int(5 * next_level[0] * next_level[1] + 1)
    else:
        next_level = None
        next_dof = None
    laptop_limit = {
        "max_completed_grid": level_rows[-1]["grid"] if level_rows else None,
        "max_completed_dof": level_rows[-1]["dof"] if level_rows else None,
        "stop_reason": stop_reason,
        "next_ratio_grid": next_level,
        "next_ratio_dof": next_dof,
        "note": (
            "Direct sparse LU is the CPU-laptop limiter; the next 2x level is the "
            "engineering sizing point for a scalable preconditioner."
        ),
    }
    pass_checks = {
        "minimum_three_levels": len(level_rows) >= study.min_levels,
        "all_completed_levels_converged": all(row["converged"] for row in level_rows),
        "fixed_ratio_ladder": True,
        "used_colored_sparse_jacobian_lu": (
            config.branch.newton.preconditioner.type == "colored_sparse_jacobian_lu"
        ),
    }
    results = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "method": {
            "refinement_ratio": study.refinement_ratio,
            "levels_configured": study.levels,
            "restriction": (
                "Nested finite-volume control volumes; fine fields are volume-averaged "
                "onto the next coarser grid before raw-field differences."
            ),
            "difference_norm": (
                "sqrt(sum_fields integral (u_h - R u_h/2)^2 dV) with dV=4*pi*r^2 dr dw; "
                "Linf is also reported."
            ),
            "placeholder_label": config.branch.placeholder_label,
            "preconditioner": asdict(config.branch.newton.preconditioner),
        },
        "level_rows": level_rows,
        "self_convergence_rows": self_rows,
        "observable_summary": observable_summary,
        "observable_table": observable_table,
        "noise_floor": noise_floor,
        "resolution_sizing": sizing,
        "laptop_limit": laptop_limit,
        "pass_checks": pass_checks,
        "passed": all(pass_checks.values()),
    }
    table_path = Path(config.run_root) / study.name / "convergence_table.json"
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


def write_step4_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    method = results["method"]
    config = results["config"]
    convergence = config["convergence"]
    lines: list[str] = []
    lines.append("# Step 4 Coupled Branch Grid-Convergence Study")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append("")
    lines.append(
        "**Discipline:** engineering smoke, placeholder parameters, not a physical packet, "
        "target-blind. Surrogates are raw-field functionals only; no extraction map is used."
    )
    lines.append("")
    lines.append("## Ladder And Norm")
    lines.append("")
    lines.append(f"Refinement ratio: `{method['refinement_ratio']}x` in both grid directions.")
    lines.append(f"Configured levels: `{convergence['levels']}`.")
    lines.append(f"Restriction: {method['restriction']}")
    lines.append(f"Difference norm: {method['difference_norm']}")
    lines.append("")
    lines.append("Preconditioner:")
    lines.append("```yaml")
    for key, value in method["preconditioner"].items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Level Performance")
    lines.append("")
    lines.append(
        _table(
            [
                "level",
                "grid",
                "dof",
                "spacing",
                "wall_clock_seconds",
                "peak_memory_mb",
                "newton_iterations",
                "final_residual_linf",
                "gmres_max",
                "gmres_mean",
                "converged",
                "message",
            ],
            results["level_rows"],
        )
    )
    lines.append("")
    lines.append("## Raw-Field Self-Convergence")
    lines.append("")
    lines.append(
        _table(
            [
                "coarse_grid",
                "fine_grid",
                "raw_field_l2_change",
                "raw_field_relative_l2_change",
                "raw_field_linf_change",
                "observed_order_from_l2_change",
                "observed_order_from_relative_l2_change",
            ],
            results["self_convergence_rows"],
        )
    )
    lines.append("")
    lines.append("## Surrogate Observable Verdicts")
    lines.append("")
    lines.append(
        _table(
            [
                "label",
                "finest_grid",
                "finest_value",
                "last_observed_order",
                "richardson_estimate",
                "finest_error_estimate",
                "verdict",
                "diagnosis",
            ],
            results["observable_summary"],
        )
    )
    lines.append("")
    lines.append("## Surrogate Per-Level Table")
    lines.append("")
    lines.append(
        _table(
            [
                "label",
                "grid",
                "value",
                "successive_change",
                "observed_order",
                "richardson_estimate",
                "error_estimate",
                "verdict",
            ],
            results["observable_table"],
        )
    )
    lines.append("")
    lines.append("## Numerical Floor")
    lines.append("")
    floor = results["noise_floor"]
    lines.append(
        "Solver floor read: "
        f"residual Linf <= {_fmt(floor['solver_residual_floor_linf'])}, "
        f"mass-constraint floor <= {_fmt(floor['mass_constraint_floor'])}."
    )
    lines.append(
        "Finest raw-field relative self-difference: "
        f"{_fmt(floor['last_raw_field_relative_l2_change'])}; "
        f"successive-difference floor reached: {_fmt(floor['field_self_difference_floor_reached'])}."
    )
    lines.append(
        f"Preliminary numerical floor for scalar diagnostics: "
        f"{_fmt(floor['preliminary_numerical_floor'])}."
    )
    lines.append("")
    lines.append("## Resolution Sizing")
    lines.append("")
    sizing_headers = [
        "label",
        "last_observed_order",
        "finest_error_estimate",
        "finest_dof",
    ] + [f"dof_for_error_{threshold:.0e}" for threshold in convergence["illustrative_error_levels"]] + [
        "verdict"
    ]
    lines.append(_table(sizing_headers, results["resolution_sizing"]))
    lines.append("")
    limit = results["laptop_limit"]
    lines.append(
        f"Laptop limit read: max completed `{limit['max_completed_grid']}` / "
        f"{limit['max_completed_dof']} DOF. Stop reason: {limit['stop_reason']}. "
        f"Next 2x grid would be {limit['next_ratio_grid']} / {limit['next_ratio_dof']} DOF."
    )
    lines.append(limit["note"])
    lines.append("")
    lines.append("## Machine-Readable Output")
    lines.append("")
    lines.append(f"Convergence table: `{results['machine_readable_table']}`")
    lines.append("")
    lines.append("Manifests:")
    for row in results["level_rows"]:
        lines.append(f"- {row['grid']}: `{row['manifest']}`")
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path
