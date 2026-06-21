"""Path-A C0f2 chunked timing rerun.

This module is a reporting/instrumentation wrapper around the C0 closed
Newton crawl.  It keeps the original ``patha_closed_branch_residual`` as the
only convergence and merit arbiter, uses the same C0 default Newton controls,
and adds timing/checkpoint/stagnation diagnostics for the C0f2 measurement run.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import time
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
import torch

from . import patha_c0_conditioning_spike as c0
from . import patha_c0f_globalization_probe as c0f
from .backend import configure_backend, jvp
from .config import BackendConfig, NewtonConfig, source_revision
from .coupled_branch import (
    initial_closed_branch_state,
    patha_closed_branch_residual,
    resample_closed_branch_state,
)


DEFAULT_C0F2_RUN_ROOT = Path("runs/pathA_C0f2_timing_rerun")
DEFAULT_C0F2_REPORT_PATH = Path("reports/pathA_C0f2_timing_rerun.md")
DEFAULT_C0F2_JSON_PATH = DEFAULT_C0F2_RUN_ROOT / "pathA_C0f2_timing_rerun.json"

DEFAULT_C0F2_DEPTH_SEQUENCE: tuple[float, ...] = (
    0.03,
    0.0295,
    0.02925,
    0.029125,
    0.0291125,
    0.0291,
    0.0290875,
    0.029075,
    0.0290625,
    0.02905,
    0.0290375,
    0.029025,
    0.0290125,
    0.029,
    0.0289375,
    0.028875,
    0.0288125,
    0.02875,
    0.0286875,
    0.028625,
    0.0285625,
    0.0285,
    0.0284375,
    0.028375,
    0.0283125,
    0.02825,
    0.0281875,
    0.028125,
    0.0280625,
    0.028,
    0.0279375,
    0.027875,
    0.0278125,
    0.02775,
)


@dataclass(frozen=True)
class C0f2Config:
    run_root: Path = DEFAULT_C0F2_RUN_ROOT
    report_path: Path = DEFAULT_C0F2_REPORT_PATH
    json_path: Path = DEFAULT_C0F2_JSON_PATH
    depth_sequence: tuple[float, ...] = DEFAULT_C0F2_DEPTH_SEQUENCE
    per_tau_attempt_timeout_seconds: float = 3600.0
    script_wall_budget_seconds: float = 570.0
    max_nominal_steps_per_run: int = 1
    b2c_linf_tolerance: float = c0.BACKGROUND_RESIDUAL_TOL
    required_unblock_tau: float = 0.028
    extrapolate_target_tau: float = 0.02
    stagnation_patience: int = 3
    stagnation_min_relative_decrease: float = 1.0e-4
    fold_growth_factor_threshold: float = 10.0
    merit_alpha_min_power: int = 20
    seed_from_c0f_checkpoint: bool = False
    c0f_checkpoint_path: Path = Path(
        "runs/pathA_C0f_globalization_probe/default_crawl/pathA_C0f_default_crawl_checkpoint.json"
    )


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _config_to_dict(config: C0f2Config) -> dict[str, Any]:
    data = asdict(config)
    for key in ("run_root", "report_path", "json_path", "c0f_checkpoint_path"):
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


def _checkpoint_path(config: C0f2Config) -> Path:
    return config.run_root / "pathA_C0f2_timing_checkpoint.json"


def _crawl_config(config: C0f2Config) -> c0.C0Config:
    return c0.C0Config(
        run_root=config.run_root,
        report_path=config.run_root / "c0_crawl_report.md",
        json_path=config.run_root / "c0_crawl.json",
        depth_sequence=config.depth_sequence,
        prefer_existing_b2c_background_predictor=False,
    )


def _state_artifact_path(config: C0f2Config, *, tau: float, attempt_index: int) -> Path:
    return (
        config.run_root
        / "states"
        / f"attempt_{attempt_index:03d}_tau_{c0._format_tau(float(tau))}.npz"
    )


def _load_checkpoint(config: C0f2Config) -> dict[str, Any] | None:
    path = _checkpoint_path(config)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    if config.seed_from_c0f_checkpoint and config.c0f_checkpoint_path.exists():
        old = json.loads(config.c0f_checkpoint_path.read_text(encoding="utf-8"))
        rows = []
        for index, row in enumerate(old.get("tau_attempts", [])):
            if not row.get("final_physical_converged"):
                continue
            seeded = dict(row)
            seeded["attempt_index"] = index
            seeded["seeded_from"] = str(config.c0f_checkpoint_path)
            seeded["timing_instrumented"] = False
            seeded.setdefault("cost_breakdown", _empty_cost_breakdown())
            rows.append(seeded)
        return _base_payload(
            config=config,
            crawl_config=_crawl_config(config),
            tau_attempts=rows,
            phase_status="seeded_from_prior_genuine_c0f_checkpoint",
            started=time.perf_counter(),
        )
    return None


def _empty_cost_breakdown() -> dict[str, Any]:
    return {
        "residual_eval_seconds": 0.0,
        "residual_eval_count": 0,
        "jacobian_assembly_seconds": 0.0,
        "jacobian_assembly_count": 0,
        "jacobian_assembly_jvp_count": 0,
        "linear_solve_seconds": 0.0,
        "linear_solve_count": 0,
        "linear_solve_jvp_seconds": 0.0,
        "linear_solve_jvp_count": 0,
        "line_search_trial_count": 0,
        "line_search_exhaustions": 0,
        "gmres_iterations": 0,
    }


def _sum_cost(left: Mapping[str, Any], right: Mapping[str, Any]) -> dict[str, Any]:
    keys = set(left) | set(right)
    total: dict[str, Any] = {}
    for key in keys:
        lval = left.get(key, 0.0)
        rval = right.get(key, 0.0)
        if isinstance(lval, (int, float)) and isinstance(rval, (int, float)):
            total[key] = lval + rval
    return total


def _timed_residual(
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x: torch.Tensor,
    profile: dict[str, Any],
) -> torch.Tensor:
    started = time.perf_counter()
    try:
        return residual_fn(x)
    finally:
        profile["residual_eval_seconds"] += time.perf_counter() - started
        profile["residual_eval_count"] += 1


def _relative_decrease(old_norm: float, new_norm: float) -> float:
    if old_norm <= 0.0 or not math.isfinite(old_norm) or not math.isfinite(new_norm):
        return math.nan
    return float((old_norm - new_norm) / old_norm)


def _solve_c0_scaled_newton_profiled(
    *,
    original_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    preconditioner_residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x0: torch.Tensor,
    grid,
    newton: NewtonConfig,
    aid: c0.C0AidParameters,
    attempt_deadline: float | None,
    stagnation_patience: int,
    stagnation_min_relative_decrease: float,
) -> tuple[c0.C0NewtonResult, dict[str, Any], str]:
    if newton.linear_solver != "gmres_jvp":
        raise ValueError(f"unsupported linear solver {newton.linear_solver!r}")
    if newton.line_search != "armijo":
        raise ValueError(f"unsupported line search {newton.line_search!r}")

    profile = _empty_cost_breakdown()
    x = x0.detach().clone()
    r = _timed_residual(original_residual_fn, x, profile).detach()
    initial_norm = c0._norm(r, newton.residual_norm)
    tolerance = max(newton.residual_atol, newton.residual_rtol * initial_norm)
    history = [c0.C0NewtonIteration(iteration=0, residual_norm=initial_norm)]
    if initial_norm <= tolerance:
        return (
            c0.C0NewtonResult(
                x=x,
                converged=True,
                iterations=0,
                initial_residual_norm=initial_norm,
                final_residual_norm=initial_norm,
                tolerance=tolerance,
                history=history,
                message="initial residual met tolerance",
            ),
            profile,
            "converged_initial",
        )

    stagnant_count = 0
    for iteration in range(1, newton.max_newton_iters + 1):
        if attempt_deadline is not None and time.perf_counter() >= attempt_deadline:
            final_norm = c0._norm(r, newton.residual_norm)
            return (
                c0.C0NewtonResult(
                    x=x,
                    converged=False,
                    iterations=iteration - 1,
                    initial_residual_norm=initial_norm,
                    final_residual_norm=final_norm,
                    tolerance=tolerance,
                    history=history,
                    message="attempt wall-clock budget reached after completed Newton iteration",
                ),
                profile,
                "attempt_budget_reached",
            )

        x_for_jvp = x.detach().clone()
        rhs = -r.detach().cpu().numpy().astype(np.float64)
        build_started = time.perf_counter()
        built, row_scale, col_scale, _matrix, _scaled_matrix = c0.build_scaled_preconditioner(
            original_residual_fn=original_residual_fn,
            preconditioner_residual_fn=preconditioner_residual_fn,
            x=x_for_jvp,
            rhs=rhs,
            grid=grid,
            newton=newton,
            aid=aid,
            iteration=iteration,
        )
        profile["jacobian_assembly_seconds"] += time.perf_counter() - build_started
        profile["jacobian_assembly_count"] += 1
        profile["jacobian_assembly_jvp_count"] += int(built.metadata.get("jvp_count", 0))

        dim = rhs.size
        rhs_scaled = row_scale * rhs
        gmres_curve: list[float] = []

        def matvec(vector_np: np.ndarray) -> np.ndarray:
            direction_np = col_scale * np.asarray(vector_np, dtype=np.float64)
            direction = torch.as_tensor(direction_np, dtype=x.dtype, device=x.device)
            jvp_started = time.perf_counter()
            try:
                jv = jvp(original_residual_fn, x_for_jvp, direction)
            finally:
                profile["linear_solve_jvp_seconds"] += time.perf_counter() - jvp_started
                profile["linear_solve_jvp_count"] += 1
            return row_scale * jv.detach().cpu().numpy().astype(np.float64)

        def callback(residual_norm: float) -> None:
            gmres_curve.append(float(residual_norm))

        linear_op = LinearOperator((dim, dim), matvec=matvec, dtype=np.float64)
        solve_started = time.perf_counter()
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
        profile["linear_solve_seconds"] += time.perf_counter() - solve_started
        profile["linear_solve_count"] += 1
        profile["gmres_iterations"] += len(gmres_curve)
        if not np.all(np.isfinite(step_scaled)):
            return (
                c0.C0NewtonResult(
                    x=x,
                    converged=False,
                    iterations=iteration - 1,
                    initial_residual_norm=initial_norm,
                    final_residual_norm=c0._norm(r, newton.residual_norm),
                    tolerance=tolerance,
                    history=history,
                    message=f"GMRES produced non-finite step, info={gmres_info}",
                ),
                profile,
                "gmres_nonfinite_step",
            )

        step_np = col_scale * np.asarray(step_scaled, dtype=np.float64)
        step = torch.as_tensor(step_np, dtype=x.dtype, device=x.device)
        step_norm = c0._norm(step, "l2")
        x_norm = max(c0._norm(x, "l2"), 1.0)
        if step_norm <= newton.step_atol + newton.step_rtol * x_norm:
            r_now = _timed_residual(original_residual_fn, x, profile).detach()
            final_norm = c0._norm(r_now, newton.residual_norm)
            return (
                c0.C0NewtonResult(
                    x=x,
                    converged=final_norm <= tolerance,
                    iterations=iteration - 1,
                    initial_residual_norm=initial_norm,
                    final_residual_norm=final_norm,
                    tolerance=tolerance,
                    history=history,
                    message="step tolerance reached",
                ),
                profile,
                "step_tolerance",
            )

        old_norm = c0._norm(r, newton.residual_norm)
        alpha = 1.0
        accepted = False
        best_x = x
        best_r = r
        best_norm = old_norm
        accepted_alpha = None
        line_search_exhausted = True
        for _ in range(newton.max_line_search_iters):
            candidate_x = x + alpha * step
            candidate_r = _timed_residual(original_residual_fn, candidate_x, profile).detach()
            profile["line_search_trial_count"] += 1
            candidate_norm = c0._norm(candidate_r, newton.residual_norm)
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
                line_search_exhausted = False
                break
            alpha *= newton.line_search_shrink

        if line_search_exhausted:
            profile["line_search_exhaustions"] += 1
        if not accepted:
            if newton.accept_best_line_search_decrease and best_norm < old_norm:
                x = best_x
                r = best_r
                accepted_alpha = accepted_alpha if accepted_alpha is not None else alpha
            else:
                history.append(
                    c0.C0NewtonIteration(
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
                return (
                    c0.C0NewtonResult(
                        x=x,
                        converged=False,
                        iterations=iteration - 1,
                        initial_residual_norm=initial_norm,
                        final_residual_norm=old_norm,
                        tolerance=tolerance,
                        history=history,
                        message="line search failed to reduce original residual",
                    ),
                    profile,
                    "line_search_exhausted",
                )

        new_norm = c0._norm(r, newton.residual_norm)
        rel_decrease = _relative_decrease(old_norm, new_norm)
        row_info = dict(built.metadata)
        row_info.update(
            {
                "jacobian_assembly_seconds_iteration": profile["jacobian_assembly_seconds"],
                "linear_solve_seconds_cumulative": profile["linear_solve_seconds"],
                "residual_eval_seconds_cumulative": profile["residual_eval_seconds"],
                "relative_residual_decrease": rel_decrease,
            }
        )
        history.append(
            c0.C0NewtonIteration(
                iteration=iteration,
                residual_norm=new_norm,
                step_norm=step_norm,
                line_search_alpha=accepted_alpha,
                gmres_info=int(gmres_info),
                gmres_iterations=len(gmres_curve),
                gmres_residual_curve=gmres_curve,
                preconditioner_info=row_info,
            )
        )
        if new_norm <= tolerance:
            return (
                c0.C0NewtonResult(
                    x=x,
                    converged=True,
                    iterations=iteration,
                    initial_residual_norm=initial_norm,
                    final_residual_norm=new_norm,
                    tolerance=tolerance,
                    history=history,
                    message="residual tolerance reached on original residual",
                ),
                profile,
                "converged",
            )

        if (
            math.isfinite(rel_decrease)
            and rel_decrease < float(stagnation_min_relative_decrease)
        ):
            stagnant_count += 1
        else:
            stagnant_count = 0
        if stagnant_count >= int(stagnation_patience):
            return (
                c0.C0NewtonResult(
                    x=x,
                    converged=False,
                    iterations=iteration,
                    initial_residual_norm=initial_norm,
                    final_residual_norm=new_norm,
                    tolerance=tolerance,
                    history=history,
                    message=(
                        "Newton stagnation detector: residual relative decrease "
                        f"< {float(stagnation_min_relative_decrease):.3e} for "
                        f"{int(stagnation_patience)} consecutive iterations"
                    ),
                ),
                profile,
                "newton_stagnation",
            )

    final_norm = c0._norm(r, newton.residual_norm)
    return (
        c0.C0NewtonResult(
            x=x,
            converged=False,
            iterations=newton.max_newton_iters,
            initial_residual_norm=initial_norm,
            final_residual_norm=final_norm,
            tolerance=tolerance,
            history=history,
            message="maximum Newton iterations reached",
        ),
        profile,
        "max_newton_iters",
    )


def _epsilon_attempt_row(
    *,
    tau: float,
    eos_K: float,
    aid_index: int,
    aid: c0.C0AidParameters,
    result: c0.C0NewtonResult,
    start_policy: str,
    metrics: Mapping[str, float],
    cost_breakdown: Mapping[str, Any],
    stop_reason: str,
) -> dict[str, Any]:
    row = c0._epsilon_attempt_row(
        tau=tau,
        eos_K=eos_K,
        aid_index=aid_index,
        aid=aid,
        result=result,
        start_policy=start_policy,
        metrics=metrics,
    )
    row["cost_breakdown"] = dict(cost_breakdown)
    row["stop_reason"] = stop_reason
    row["timing_instrumented"] = True
    return row


def _execute_tau_attempt_profiled(
    *,
    config: C0f2Config,
    crawl_config: c0.C0Config,
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
    branch, provider, grid, boundaries = c0._branch_context(
        tau=float(target_tau),
        config=crawl_config,
        dtype=dtype,
    )
    used_existing_b2c = False
    if previous_state is None:
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
        if crawl_config.use_wall_predictor:
            state, predictor = c0._apply_c0_wall_predictor(
                state=state,
                grid=grid,
                tau=float(target_tau),
                grid_level=crawl_config.grid,
            )
            initialization["wall_predictor"] = predictor

    clean_attempt_state = state.detach().clone()
    final_only = bool(
        crawl_config.eos_final_only_after_first
        and initialization.get("source") != "default_closed_branch_state"
    )
    continuation_values = (
        (float(branch.continuation_K_values[-1]),)
        if final_only
        else tuple(float(value) for value in branch.continuation_K_values)
    )
    epsilon_attempts: list[dict[str, Any]] = []
    stage_converged = True
    stage_message = "C0f2 tau attempt completed"
    stage_started = time.perf_counter()
    deadline = stage_started + float(config.per_tau_attempt_timeout_seconds)
    state_for_metrics = clean_attempt_state
    attempt_cost = _empty_cost_breakdown()
    stop_reason = "not_started"

    for eos_K in continuation_values:
        clean_eos_state = state.detach().clone()
        last_converged_epsilon_state: torch.Tensor | None = None
        zero_epsilon_result: c0.C0NewtonResult | None = None
        last_result: c0.C0NewtonResult | None = None
        for aid_index, aid in enumerate(crawl_config.epsilon_schedule):
            if last_converged_epsilon_state is None:
                x0 = clean_eos_state
                start_policy = "clean_tau_attempt_state"
            else:
                x0 = last_converged_epsilon_state
                start_policy = "previous_converged_epsilon_state"

            def original_residual_fn(x: torch.Tensor, eos_K: float = eos_K) -> torch.Tensor:
                return patha_closed_branch_residual(
                    x,
                    grid,
                    branch,
                    eos_K=float(eos_K),
                    boundaries=boundaries,
                    s_sigma=provider,
                )

            def preconditioner_residual_fn(
                x: torch.Tensor,
                eos_K: float = eos_K,
                aid: c0.C0AidParameters = aid,
            ) -> torch.Tensor:
                return c0.c0_preconditioner_residual(
                    x,
                    grid,
                    branch,
                    eos_K=float(eos_K),
                    boundaries=boundaries,
                    s_sigma=provider,
                    aid=aid,
                )

            result, cost, stop_reason = _solve_c0_scaled_newton_profiled(
                original_residual_fn=original_residual_fn,
                preconditioner_residual_fn=preconditioner_residual_fn,
                x0=x0,
                grid=grid,
                newton=branch.newton,
                aid=aid,
                attempt_deadline=deadline,
                stagnation_patience=config.stagnation_patience,
                stagnation_min_relative_decrease=config.stagnation_min_relative_decrease,
            )
            attempt_cost = _sum_cost(attempt_cost, cost)
            state_for_metrics = result.x.detach()
            metrics = c0._state_metrics(state_for_metrics, grid)
            epsilon_attempts.append(
                _epsilon_attempt_row(
                    tau=float(target_tau),
                    eos_K=float(eos_K),
                    aid_index=aid_index,
                    aid=aid,
                    result=result,
                    start_policy=start_policy,
                    metrics=metrics,
                    cost_breakdown=cost,
                    stop_reason=stop_reason,
                )
            )
            last_result = result
            if result.converged:
                last_converged_epsilon_state = result.x.detach()
            if aid.core_epsilon == 0.0 and aid.k1_radius_epsilon == 0.0:
                zero_epsilon_result = result
                if result.converged:
                    state = result.x.detach()
            if not result.converged and stop_reason in {
                "attempt_budget_reached",
                "line_search_exhausted",
                "newton_stagnation",
            }:
                break
        if zero_epsilon_result is None or not zero_epsilon_result.converged:
            stage_converged = False
            failed = zero_epsilon_result or last_result
            failed_message = "zero-epsilon solve did not run"
            if failed is not None:
                failed_message = failed.message
            stage_message = (
                f"tau={target_tau:.12e}, eos_K={float(eos_K):.12e}, "
                f"epsilon=0 failed or was not reached: {failed_message}"
            )
            break

    final_eos = float(continuation_values[-1])

    def final_residual_fn(x: torch.Tensor) -> torch.Tensor:
        return patha_closed_branch_residual(
            x,
            grid,
            branch,
            eos_K=final_eos,
            boundaries=boundaries,
            s_sigma=provider,
        )

    final_state = state if stage_converged else state_for_metrics
    final_started = time.perf_counter()
    final_residual = final_residual_fn(final_state).detach()
    attempt_cost["residual_eval_seconds"] += time.perf_counter() - final_started
    attempt_cost["residual_eval_count"] += 1
    final_linf = float(torch.max(torch.abs(final_residual)).detach().cpu().item())
    final_l2 = float(torch.linalg.vector_norm(final_residual).detach().cpu().item())
    metrics = c0._state_metrics(final_state, grid)
    final_physical_converged = bool(stage_converged and final_linf <= c0.BACKGROUND_RESIDUAL_TOL)
    state_path = _state_artifact_path(config, tau=float(target_tau), attempt_index=attempt_index)
    c0._save_state_artifact(state_path, final_state)
    elapsed = time.perf_counter() - stage_started
    row = {
        "attempt_index": int(attempt_index),
        "tau": float(target_tau),
        "target_tau": float(target_tau),
        "nominal_target_tau": float(nominal_target_tau),
        "delta_tau": None if delta_tau is None else float(delta_tau),
        "backtrack_index": int(backtrack_index),
        "start_from_tau": previous_tau,
        "target_relative_to_prior_floor": float(target_tau) / c0.PRIOR_TAU_FLOOR,
        "below_prior_floor": c0._below_prior_floor(float(target_tau)),
        "stage_converged": bool(stage_converged),
        "final_original_residual_linf": final_linf,
        "final_original_residual_l2": final_l2,
        "final_physical_converged": final_physical_converged,
        "b2c_background_tolerance": c0.BACKGROUND_RESIDUAL_TOL,
        "message": stage_message,
        "elapsed_seconds": elapsed,
        "wall_clock_seconds": elapsed,
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
            "core_epsilon_schedule": [aid.core_epsilon for aid in crawl_config.epsilon_schedule],
            "k1_radius_epsilon_schedule": [
                aid.k1_radius_epsilon for aid in crawl_config.epsilon_schedule
            ],
            "min_rho_during_final": metrics["min_rho"],
            "min_R0_during_final": metrics["min_R0"],
            "k1_clamp_active_in_path": any(
                aid.k1_radius_epsilon > 0.0 and metrics["min_R0"] < aid.k1_radius_epsilon
                for aid in crawl_config.epsilon_schedule
            ),
            "final_aids_inactive": True,
        },
        "cost_breakdown": attempt_cost,
        "timing_instrumented": True,
        "stop_reason": stop_reason,
        "stalled": bool(
            (not final_physical_converged)
            and stop_reason
            in {"attempt_budget_reached", "line_search_exhausted", "newton_stagnation"}
        ),
    }
    return row, final_state.detach(), grid


def _last_converged_state(
    rows: Sequence[Mapping[str, Any]],
    *,
    crawl_config: c0.C0Config,
    dtype: torch.dtype,
) -> tuple[torch.Tensor | None, Any | None, float | None]:
    converged = [row for row in rows if row.get("final_physical_converged")]
    if not converged:
        return None, None, None
    deepest = min(converged, key=lambda row: float(row["target_tau"]))
    tau = float(deepest["target_tau"])
    _branch, _provider, grid, _boundaries = c0._branch_context(
        tau=tau,
        config=crawl_config,
        dtype=dtype,
    )
    state = c0._load_state_artifact(c0f._resolve_existing_path(str(deepest["state_artifact"])), dtype=dtype)
    return state.detach().clone(), grid, tau


def _selected_next_nominal(
    rows: Sequence[Mapping[str, Any]],
    sequence: Sequence[float],
) -> float | None:
    converged = [row for row in rows if row.get("final_physical_converged")]
    deepest = None if not converged else min(float(row["target_tau"]) for row in converged)
    for tau in sequence:
        if deepest is None or float(tau) < deepest - 5.0e-13:
            return float(tau)
    return None


def _base_payload(
    *,
    config: C0f2Config,
    crawl_config: c0.C0Config,
    tau_attempts: Sequence[Mapping[str, Any]],
    phase_status: str,
    started: float,
) -> dict[str, Any]:
    return {
        "schema": "stage1_pathA_C0f2_timing_rerun/v1",
        "source_revision": source_revision(),
        "phase_status": phase_status,
        "config": _config_to_dict(config),
        "crawl_config": c0._config_to_dict(crawl_config),
        "tau_attempts": [dict(row) for row in tau_attempts],
        "elapsed_seconds": time.perf_counter() - started,
        "scope_guard": _scope_guard(crawl_config),
    }


def _write_checkpoint(
    *,
    config: C0f2Config,
    crawl_config: c0.C0Config,
    rows: Sequence[Mapping[str, Any]],
    phase_status: str,
    started: float,
) -> dict[str, Any]:
    payload = _base_payload(
        config=config,
        crawl_config=crawl_config,
        tau_attempts=rows,
        phase_status=phase_status,
        started=started,
    )
    path = _checkpoint_path(config)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n", encoding="utf-8")
    return payload


def _run_crawl_chunk(
    *,
    config: C0f2Config,
    crawl_config: c0.C0Config,
    dtype: torch.dtype,
    existing_rows: Sequence[Mapping[str, Any]],
    started: float,
) -> tuple[list[dict[str, Any]], str]:
    rows = [dict(row) for row in existing_rows]
    last_state, last_grid, deepest_tau = _last_converged_state(rows, crawl_config=crawl_config, dtype=dtype)
    attempt_index = len(rows)
    nominal_steps_completed = 0
    chunk_started = time.perf_counter()
    phase_status = "crawl_in_progress"

    while nominal_steps_completed < int(config.max_nominal_steps_per_run):
        nominal_tau = _selected_next_nominal(rows, config.depth_sequence)
        if nominal_tau is None:
            phase_status = "crawl_complete_depth_sequence_exhausted"
            break
        start_tau = deepest_tau
        full_delta = None if start_tau is None else float(nominal_tau) - float(start_tau)
        target_tau = float(nominal_tau)
        backtrack_index = 0
        accepted_this_nominal = False
        exhausted_this_nominal = False

        while True:
            delta_tau = None if start_tau is None else float(target_tau) - float(start_tau)
            row, final_state, final_grid = _execute_tau_attempt_profiled(
                config=config,
                crawl_config=crawl_config,
                dtype=dtype,
                target_tau=float(target_tau),
                nominal_target_tau=float(nominal_tau),
                delta_tau=delta_tau,
                backtrack_index=backtrack_index,
                attempt_index=attempt_index,
                previous_state=last_state,
                previous_grid=last_grid,
                previous_tau=deepest_tau,
            )
            rows.append(row)
            attempt_index += 1
            _write_checkpoint(
                config=config,
                crawl_config=crawl_config,
                rows=rows,
                phase_status="crawl_in_progress",
                started=started,
            )
            if row["final_physical_converged"]:
                last_state = final_state.detach().clone()
                last_grid = final_grid
                deepest_tau = float(row["target_tau"])
                accepted_this_nominal = True
                break
            if start_tau is None or full_delta is None:
                exhausted_this_nominal = True
                break
            if (
                backtrack_index >= crawl_config.max_tau_backtracks
                or abs(delta_tau or 0.0) <= crawl_config.min_tau_backtrack_delta
            ):
                exhausted_this_nominal = True
                break
            backtrack_index += 1
            target_tau = float(start_tau) + float(full_delta) / (2.0**backtrack_index)

        nominal_steps_completed += 1
        if exhausted_this_nominal and not accepted_this_nominal:
            phase_status = "crawl_stalled_backtrack_exhausted"
            break
        if deepest_tau is not None and deepest_tau <= min(config.depth_sequence) + 5.0e-13:
            phase_status = "crawl_complete_depth_sequence_exhausted"
            break
        if time.perf_counter() - chunk_started >= float(config.script_wall_budget_seconds):
            phase_status = "crawl_chunk_budget_reached"
            break
    return rows, phase_status


def _zero_attempt(row: Mapping[str, Any]) -> dict[str, Any] | None:
    return c0f._final_zero_epsilon_attempt(row)


def _representative_attempt(row: Mapping[str, Any]) -> dict[str, Any] | None:
    zero = _zero_attempt(row)
    if zero is not None:
        return zero
    attempts = row.get("epsilon_attempts", [])
    return dict(attempts[-1]) if attempts else None


def _line_search_backtracks(row: Mapping[str, Any]) -> int:
    attempt = _representative_attempt(row)
    alphas = c0f._accepted_alphas(attempt)
    branch = c0.c0_frozen_branch(tau=float(row["target_tau"]), grid=c0.C0Config().grid)
    return c0f._line_search_backtrack_count(alphas, shrink=float(branch.newton.line_search_shrink))


def _per_tau_timing(rows: Sequence[Mapping[str, Any]], *, tolerance: float) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        zero = _zero_attempt(row)
        representative = _representative_attempt(row)
        cost = row.get("cost_breakdown", {})
        linf = row.get("final_original_residual_linf")
        l2 = row.get("final_original_residual_l2")
        solver_converged = bool(zero.get("converged")) if zero is not None else bool(row.get("stage_converged"))
        table.append(
            {
                "attempt_index": int(row.get("attempt_index", index)),
                "tau": float(row["target_tau"]),
                "nominal_tau": float(row.get("nominal_target_tau", row["target_tau"])),
                "start_tau": row.get("start_from_tau"),
                "tau_backtrack_index": int(row.get("backtrack_index", 0)),
                "wall_clock_seconds": float(row.get("wall_clock_seconds", row.get("elapsed_seconds", math.nan))),
                "solver_converged": solver_converged,
                "b2c_linf_pass": bool(linf is not None and math.isfinite(float(linf)) and float(linf) <= tolerance),
                "accepted_default_success": bool(
                    solver_converged
                    and linf is not None
                    and math.isfinite(float(linf))
                    and float(linf) <= tolerance
                ),
                "linf": None if linf is None else float(linf),
                "l2": None if l2 is None else float(l2),
                "newton_iterations": int(representative.get("iterations", 0))
                if representative is not None
                else None,
                "line_search_backtracks": _line_search_backtracks(row),
                "jacobian_assembly_seconds": float(cost.get("jacobian_assembly_seconds", math.nan)),
                "linear_solve_seconds": float(cost.get("linear_solve_seconds", math.nan)),
                "residual_eval_seconds": float(cost.get("residual_eval_seconds", math.nan)),
                "linear_solve_jvp_seconds": float(cost.get("linear_solve_jvp_seconds", math.nan)),
                "jacobian_assembly_jvp_count": int(cost.get("jacobian_assembly_jvp_count", 0)),
                "linear_solve_jvp_count": int(cost.get("linear_solve_jvp_count", 0)),
                "gmres_iterations": int(cost.get("gmres_iterations", 0)),
                "stop_reason": row.get("stop_reason"),
                "message": row.get("message"),
                "state_artifact": row.get("state_artifact"),
                "init_source": row.get("init", {}).get("source"),
                "used_existing_b2c": bool(row.get("used_existing_b2c")),
                "timing_instrumented": bool(row.get("timing_instrumented")),
            }
        )
    return table


def _deepest_linf_pass_tau(table: Sequence[Mapping[str, Any]]) -> float | None:
    passed = [float(row["tau"]) for row in table if row.get("b2c_linf_pass")]
    return min(passed) if passed else None


def _deepest_accepted_tau(table: Sequence[Mapping[str, Any]]) -> float | None:
    passed = [float(row["tau"]) for row in table if row.get("accepted_default_success")]
    return min(passed) if passed else None


def _reached_tau(table: Sequence[Mapping[str, Any]], tau: float) -> bool:
    return any(row.get("accepted_default_success") and float(row["tau"]) <= float(tau) + 5.0e-13 for row in table)


def _first_stalled_row(rows: Sequence[Mapping[str, Any]], table: Sequence[Mapping[str, Any]], config: C0f2Config) -> dict[str, Any] | None:
    if _reached_tau(table, config.required_unblock_tau):
        return None
    table_by_index = {int(row["attempt_index"]): row for row in table}
    for index, row in enumerate(rows):
        summary = table_by_index.get(int(row.get("attempt_index", index)), {})
        if not summary.get("accepted_default_success"):
            out = dict(row)
            out["_attempt_index"] = int(row.get("attempt_index", index))
            out["_selection_reason"] = "first_c0f2_attempt_without_solver_converged_and_linf_pass"
            return out
    return None


def _extrapolate(table: Sequence[Mapping[str, Any]], *, target_tau: float) -> dict[str, Any]:
    measured = [
        row
        for row in table
        if row.get("accepted_default_success")
        and row.get("timing_instrumented")
        and math.isfinite(float(row.get("wall_clock_seconds", math.nan)))
    ]
    if len(measured) < 2:
        return {
            "status": "NOT_MEASURED",
            "reason": "fewer_than_two_measured_converged_steps",
            "assumed_target_tau": float(target_tau),
            "model": "recent_median_per_fine_tau_step",
        }
    ordered = sorted(measured, key=lambda row: float(row["tau"]), reverse=True)
    step_sizes = [
        abs(float(right["tau"]) - float(left["tau"]))
        for left, right in zip(ordered, ordered[1:])
        if abs(float(right["tau"]) - float(left["tau"])) > 0.0
    ]
    deepest = min(float(row["tau"]) for row in measured)
    recent = sorted(measured, key=lambda row: float(row["tau"]))[: min(5, len(measured))]
    recent_times = [float(row["wall_clock_seconds"]) for row in recent]
    median_time = float(np.median(recent_times))
    median_step = float(np.median(step_sizes)) if step_sizes else math.nan
    remaining_depth = max(0.0, deepest - float(target_tau))
    remaining_steps = math.ceil(remaining_depth / median_step) if median_step > 0.0 else math.nan
    measured_elapsed = float(sum(float(row["wall_clock_seconds"]) for row in measured))
    remaining_seconds = float(remaining_steps * median_time) if math.isfinite(remaining_steps) else math.nan
    return {
        "status": "MEASURED",
        "assumed_target_tau": float(target_tau),
        "model": "recent_median_per_fine_tau_step",
        "model_detail": (
            "Uses the median wall-clock over the deepest measured converged steps "
            "and the median accepted tau spacing; no claim of asymptotic validity."
        ),
        "deepest_measured_tau": deepest,
        "median_recent_step_seconds": median_time,
        "median_tau_step": median_step,
        "remaining_steps": int(remaining_steps) if math.isfinite(remaining_steps) else None,
        "measured_elapsed_seconds": measured_elapsed,
        "remaining_seconds": remaining_seconds,
        "estimated_total_seconds_from_tau_0p03": measured_elapsed + remaining_seconds,
    }


def _slow_vs_fold_call(
    *,
    table: Sequence[Mapping[str, Any]],
    fold: Mapping[str, Any],
    merit: Mapping[str, Any],
    config: C0f2Config,
) -> dict[str, Any]:
    reached_required = _reached_tau(table, config.required_unblock_tau)
    if reached_required:
        return {
            "call": "WALL_WAS_CONFIG",
            "reason": "accepted crawl reached required tau with Linf tolerance",
            "reached_required_tau": True,
            "fold_call": fold.get("call"),
        }
    if fold.get("call") == "FOLD_DETECTED":
        return {
            "call": "FOLD_DETECTED",
            "reason": "fold growth gate and smaller-tau failure gate both held",
            "fold_growth_factor": fold.get("growth_factor"),
            "fold_growth_threshold": fold.get("growth_threshold"),
        }
    if merit.get("status") == "MEASURED":
        if merit.get("any_alpha_reduces_true_l2"):
            return {
                "call": "SLOW_BUT_DESCENT_DIRECTION_EXISTS",
                "reason": "merit sweep found true residual decrease for at least one alpha",
                "linear_rel_resid": merit.get("linear_rel_resid"),
                "best_alpha": merit.get("best_alpha"),
                "best_actual_l2": merit.get("best_actual_l2"),
                "fold_call": fold.get("call"),
                "fold_growth_factor": fold.get("growth_factor"),
            }
        return {
            "call": "GLOBALIZATION_INSUFFICIENT_OR_FOLD_GEOMETRY",
            "reason": "merit sweep found no alpha reducing true residual",
            "linear_rel_resid": merit.get("linear_rel_resid"),
            "fold_call": fold.get("call"),
            "fold_growth_factor": fold.get("growth_factor"),
        }
    return {
        "call": "DIAGNOSTIC_INCOMPLETE",
        "reason": "required tau not reached and stall merit evidence not measured",
        "fold_call": fold.get("call"),
        "merit_status": merit.get("status"),
    }


def _scope_guard(crawl_config: c0.C0Config) -> dict[str, Any]:
    default = c0.C0Config()
    branch = c0.c0_frozen_branch(tau=default.depth_sequence[0], grid=default.grid)
    return {
        "genuine_continuation": True,
        "prefer_existing_b2c_background_predictor": bool(crawl_config.prefer_existing_b2c_background_predictor),
        "used_existing_b2c_permitted": False,
        "warm_start_rule": "previous_c0_converged_state",
        "single_arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
        "max_newton_iters": int(branch.newton.max_newton_iters),
        "max_tau_backtracks": int(default.max_tau_backtracks),
        "line_search": branch.newton.line_search,
        "max_line_search_iters": int(branch.newton.max_line_search_iters),
        "operators_touched_by_c0f2": False,
        "frozen_physics_touched_by_c0f2": False,
        "physical_export_guard_touched_by_c0f2": False,
        "solver_logic_changed_by_c0f2": False,
        "trust_region_dogleg_lm_implemented": False,
        "pseudo_arclength_implemented": False,
        "gauge_deflation_implemented": False,
        "depth_continuation": "tau_only",
    }


def _finalize_result(
    *,
    config: C0f2Config,
    crawl_config: c0.C0Config,
    rows: Sequence[Mapping[str, Any]],
    phase_status: str,
    started: float,
) -> dict[str, Any]:
    table = _per_tau_timing(rows, tolerance=config.b2c_linf_tolerance)
    crawl_result = {
        "tau_attempts": [dict(row) for row in rows],
    }
    fold = c0f._compute_fold_diagnostic(
        crawl_result,
        table,
        crawl_config=crawl_config,
        dtype=configure_backend(BackendConfig()),
        growth_threshold=config.fold_growth_factor_threshold,
    )
    stalled_row = _first_stalled_row(rows, table, config)
    if stalled_row is not None and phase_status.startswith("crawl_stalled"):
        merit = c0f._run_merit_sweep(
            stalled_row,
            crawl_config=crawl_config,
            run_root=config.run_root,
            dtype=configure_backend(BackendConfig()),
            alpha_min_power=config.merit_alpha_min_power,
        )
    else:
        merit = {"status": "SKIPPED", "reason": "no stalled row requiring merit sweep"}
    extrapolation = _extrapolate(table, target_tau=config.extrapolate_target_tau)
    if phase_status.startswith("crawl_stalled"):
        extrapolation = dict(extrapolation)
        extrapolation["blocked_by_measured_stall"] = True
        extrapolation["validity"] = "conditional_trend_only_not_a_completion_prediction"
    slow_vs_fold = _slow_vs_fold_call(
        table=table,
        fold=fold,
        merit=merit,
        config=config,
    )
    result = _base_payload(
        config=config,
        crawl_config=crawl_config,
        tau_attempts=rows,
        phase_status=phase_status,
        started=started,
    )
    result.update(
        {
            "per_tau_timing": table,
            "deepest_linf_pass_tau": _deepest_linf_pass_tau(table),
            "deepest_accepted_tau": _deepest_accepted_tau(table),
            "reached_tau_0p028": _reached_tau(table, config.required_unblock_tau),
            "fold_diagnostic": fold,
            "merit_sweep": merit,
            "slow_vs_fold_call": slow_vs_fold,
            "extrapolation": extrapolation,
            "report_path": str(config.report_path),
            "json_path": str(config.json_path),
            "checkpoint_path": str(_checkpoint_path(config)),
            "files_changed": [
                "src/stage1_solver/patha_c0f2_timing_rerun.py",
                "scripts/pathA_C0f2_timing_rerun.py",
                "tests/test_patha_c0f2_timing_rerun.py",
                "reports/pathA_C0f2_timing_rerun.md",
                "runs/pathA_C0f2_timing_rerun/pathA_C0f2_timing_rerun.json",
            ],
        }
    )
    config.json_path.parent.mkdir(parents=True, exist_ok=True)
    config.json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    _write_report(result, config.report_path)
    return result


def _write_report(result: Mapping[str, Any], path: Path) -> None:
    lines: list[str] = []
    lines.append("# Path-A C0f2 Timing Rerun")
    lines.append("")
    lines.append(f"Phase status: **{result.get('phase_status')}**")
    lines.append(f"Slow-vs-fold call: **{result.get('slow_vs_fold_call', {}).get('call')}**")
    lines.append("")
    lines.append("## Per-Tau Timing Table")
    lines.append("")
    rows = []
    for row in result.get("per_tau_timing", []):
        rows.append(
            {
                "tau": row.get("tau"),
                "wall_s": row.get("wall_clock_seconds"),
                "iters": row.get("newton_iterations"),
                "bt": row.get("line_search_backtracks"),
                "Linf": row.get("linf"),
                "L2": row.get("l2"),
                "J_asm_s": row.get("jacobian_assembly_seconds"),
                "lin_solve_s": row.get("linear_solve_seconds"),
                "resid_s": row.get("residual_eval_seconds"),
                "JVP_asm": row.get("jacobian_assembly_jvp_count"),
                "JVP_lin": row.get("linear_solve_jvp_count"),
                "status": row.get("stop_reason"),
            }
        )
    lines.append(
        _markdown_table(
            [
                "tau",
                "wall_s",
                "iters",
                "bt",
                "Linf",
                "L2",
                "J_asm_s",
                "lin_solve_s",
                "resid_s",
                "JVP_asm",
                "JVP_lin",
                "status",
            ],
            rows,
        )
    )
    lines.append("")
    lines.append("## Depth And Extrapolation")
    lines.append("")
    lines.append("```yaml")
    lines.append(f"deepest_linf_pass_tau: {result.get('deepest_linf_pass_tau')}")
    lines.append(f"deepest_accepted_tau: {result.get('deepest_accepted_tau')}")
    lines.append(f"reached_tau_0p028: {result.get('reached_tau_0p028')}")
    for key, value in result.get("extrapolation", {}).items():
        lines.append(f"extrapolation_{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Slow Vs Fold Evidence")
    lines.append("")
    lines.append("```yaml")
    for key, value in result.get("slow_vs_fold_call", {}).items():
        lines.append(f"{key}: {value}")
    lines.append("fold:")
    fold = result.get("fold_diagnostic", {})
    for key in ("status", "call", "reason", "growth_factor", "growth_threshold", "growth_condition"):
        lines.append(f"  {key}: {fold.get(key)}")
    merit = result.get("merit_sweep", {})
    lines.append("merit:")
    for key in (
        "status",
        "reason",
        "tau",
        "initial_l2",
        "initial_linf",
        "linear_rel_resid",
        "any_alpha_reduces_true_l2",
        "best_alpha",
        "best_actual_l2",
        "diagnostic_interpretation",
    ):
        lines.append(f"  {key}: {merit.get(key)}")
    lines.append("```")
    merit_rows = [
        {
            "alpha": row.get("alpha"),
            "actual_l2": row.get("actual_l2"),
            "pred_l2": row.get("predicted_l2"),
            "actual_linf": row.get("actual_linf"),
            "pred_linf": row.get("predicted_linf"),
            "reduces": row.get("reduces_true_l2"),
        }
        for row in merit.get("alpha_rows", [])
    ]
    if merit_rows:
        lines.append("")
        lines.append(
            _markdown_table(
                ["alpha", "actual_l2", "pred_l2", "actual_linf", "pred_linf", "reduces"],
                merit_rows,
            )
        )
    lines.append("")
    lines.append("## Genuineness And Scope")
    lines.append("")
    lines.append("```yaml")
    for key, value in result.get("scope_guard", {}).items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append(f"Machine JSON: `{result.get('json_path')}`")
    lines.append(f"Checkpoint JSON: `{result.get('checkpoint_path')}`")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_c0f2_timing_rerun(config: C0f2Config | None = None) -> dict[str, Any]:
    if config is None:
        config = C0f2Config()
    started = time.perf_counter()
    config.run_root.mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(BackendConfig())
    crawl_config = _crawl_config(config)
    checkpoint = _load_checkpoint(config)
    existing_rows = [] if checkpoint is None else [dict(row) for row in checkpoint.get("tau_attempts", [])]
    if config.max_nominal_steps_per_run <= 0 and checkpoint is not None:
        rows = existing_rows
        phase_status = str(checkpoint.get("phase_status", "crawl_in_progress"))
    else:
        rows, phase_status = _run_crawl_chunk(
            config=config,
            crawl_config=crawl_config,
            dtype=dtype,
            existing_rows=existing_rows,
            started=started,
        )
    _write_checkpoint(
        config=config,
        crawl_config=crawl_config,
        rows=rows,
        phase_status=phase_status,
        started=started,
    )
    return _finalize_result(
        config=config,
        crawl_config=crawl_config,
        rows=rows,
        phase_status=phase_status,
        started=started,
    )


def _stdout_summary(result: Mapping[str, Any]) -> str:
    lines: list[str] = []
    lines.append("Per-tau timing table:")
    for row in result.get("per_tau_timing", []):
        lines.append(
            "  "
            f"tau={float(row['tau']):.12e} "
            f"wall={_fmt(row.get('wall_clock_seconds'))}s "
            f"iters={row.get('newton_iterations')} "
            f"bt={row.get('line_search_backtracks')} "
            f"Linf={_fmt(row.get('linf'))} "
            f"L2={_fmt(row.get('l2'))} "
            f"Jasm={_fmt(row.get('jacobian_assembly_seconds'))}s "
            f"lin={_fmt(row.get('linear_solve_seconds'))}s "
            f"resid={_fmt(row.get('residual_eval_seconds'))}s "
            f"status={row.get('stop_reason')}"
        )
    lines.append(
        f"Deepest tau with Linf<=1e-6: {_fmt(result.get('deepest_linf_pass_tau'))}; "
        f"reached tau=0.028: {result.get('reached_tau_0p028')}"
    )
    extra = result.get("extrapolation", {})
    lines.append(
        "Extrapolated full-run estimate: "
        f"target_tau={extra.get('assumed_target_tau')} "
        f"model={extra.get('model')} "
        f"total_seconds={_fmt(extra.get('estimated_total_seconds_from_tau_0p03'))} "
        f"remaining_seconds={_fmt(extra.get('remaining_seconds'))}"
    )
    call = result.get("slow_vs_fold_call", {})
    fold = result.get("fold_diagnostic", {})
    merit = result.get("merit_sweep", {})
    lines.append(
        "Slow-vs-fold call: "
        f"{call.get('call')} ({call.get('reason')}); "
        f"fold_growth={_fmt(fold.get('growth_factor'))}; "
        f"merit_status={merit.get('status')} "
        f"linear_rel_resid={_fmt(merit.get('linear_rel_resid'))} "
        f"any_alpha_reduces={merit.get('any_alpha_reduces_true_l2')}"
    )
    guard = result.get("scope_guard", {})
    lines.append(
        "Genuine continuation confirmed: "
        f"prefer_existing={guard.get('prefer_existing_b2c_background_predictor')}, "
        "warm-start from prior converged state, no B2c cold-load across resumes; "
        "Single-Arbiter honored; operators/frozen/export/solver logic untouched; "
        "no trust-region/dogleg/LM, pseudo-arclength, or deflation."
    )
    lines.append("Files changed: " + ", ".join(str(item) for item in result.get("files_changed", [])))
    lines.append(f"Artifacts: report={result.get('report_path')}, json={result.get('json_path')}")
    return "\n".join(lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Path-A C0f2 chunked timing rerun.")
    parser.add_argument("--run-root", type=Path, default=DEFAULT_C0F2_RUN_ROOT)
    parser.add_argument("--report-path", type=Path, default=DEFAULT_C0F2_REPORT_PATH)
    parser.add_argument("--json-path", type=Path, default=DEFAULT_C0F2_JSON_PATH)
    parser.add_argument("--per-tau-attempt-timeout-seconds", type=float, default=C0f2Config().per_tau_attempt_timeout_seconds)
    parser.add_argument("--script-wall-budget-seconds", type=float, default=C0f2Config().script_wall_budget_seconds)
    parser.add_argument("--max-nominal-steps-per-run", type=int, default=C0f2Config().max_nominal_steps_per_run)
    parser.add_argument("--extrapolate-target-tau", type=float, default=C0f2Config().extrapolate_target_tau)
    parser.add_argument("--stagnation-patience", type=int, default=C0f2Config().stagnation_patience)
    parser.add_argument("--stagnation-min-relative-decrease", type=float, default=C0f2Config().stagnation_min_relative_decrease)
    parser.add_argument("--seed-from-c0f-checkpoint", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    config = C0f2Config(
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        per_tau_attempt_timeout_seconds=float(args.per_tau_attempt_timeout_seconds),
        script_wall_budget_seconds=float(args.script_wall_budget_seconds),
        max_nominal_steps_per_run=int(args.max_nominal_steps_per_run),
        extrapolate_target_tau=float(args.extrapolate_target_tau),
        stagnation_patience=int(args.stagnation_patience),
        stagnation_min_relative_decrease=float(args.stagnation_min_relative_decrease),
        seed_from_c0f_checkpoint=bool(args.seed_from_c0f_checkpoint),
    )
    result = run_c0f2_timing_rerun(config)
    print(_stdout_summary(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
