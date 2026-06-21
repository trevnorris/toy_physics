"""Path-A C0f globalization-first probe.

This module drives the existing C0 closed-Newton crawl with code defaults and
adds reporting-only diagnostics required by the C0f directive.  It does not
modify the Newton solver, the physical residual, the operators, frozen physics,
or the physical export guard.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import argparse
import json
import math
from pathlib import Path
import signal
import time
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.sparse import save_npz
import torch

from . import patha_c0_conditioning_spike as c0
from .backend import configure_backend
from .config import BackendConfig, source_revision
from .coupled_branch import patha_closed_branch_residual, unpack_closed_coupled_fields


DEFAULT_C0F_RUN_ROOT = Path("runs/pathA_C0f_globalization_probe")
DEFAULT_C0F_REPORT_PATH = Path("reports/pathA_C0f_globalization_probe.md")
DEFAULT_C0F_JSON_PATH = DEFAULT_C0F_RUN_ROOT / "pathA_C0f_globalization_probe.json"

C0F_VERDICTS = {
    "WALL_WAS_CONFIG",
    "FOLD_DETECTED",
    "GLOBALIZATION_INSUFFICIENT",
    "DIAGNOSTIC_INCOMPLETE",
}


@dataclass(frozen=True)
class C0fConfig:
    run_root: Path = DEFAULT_C0F_RUN_ROOT
    report_path: Path = DEFAULT_C0F_REPORT_PATH
    json_path: Path = DEFAULT_C0F_JSON_PATH
    c0b_json_path: Path = Path("runs/pathA_C0b_wall_diagnosis/pathA_C0b_diagnostic.json")
    fold_growth_factor_threshold: float = 10.0
    b2c_linf_tolerance: float = c0.BACKGROUND_RESIDUAL_TOL
    wall_was_config_target_tau: float = 0.028
    merit_alpha_min_power: int = 20
    mixed_control_epsilon: float = 0.1
    mixed_control_k_values: tuple[int, ...] = (1, 2, 4, 7)
    gauge_fraction_threshold: float = 2.0e-2
    transverse_fraction_threshold: float = 5.0e-2
    per_tau_attempt_timeout_seconds: float = 540.0


class _AttemptTimeout(RuntimeError):
    pass


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _config_to_dict(config: C0fConfig) -> dict[str, Any]:
    data = asdict(config)
    for key in ("run_root", "report_path", "json_path", "c0b_json_path"):
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


def _stage_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_existing_path(path: str | Path) -> Path:
    raw = Path(path)
    if raw.is_absolute():
        return raw
    candidates = [Path.cwd() / raw, _stage_root() / raw]
    if raw.parts[:2] == ("software", "stage1_solver"):
        candidates.append(Path(__file__).resolve().parents[4] / raw)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _final_zero_epsilon_attempt(row: Mapping[str, Any]) -> dict[str, Any] | None:
    zero = [
        dict(attempt)
        for attempt in row.get("epsilon_attempts", [])
        if abs(float(attempt.get("core_epsilon", attempt.get("epsilon", math.nan)))) <= 0.0
        and abs(float(attempt.get("k1_radius_epsilon", math.nan))) <= 0.0
    ]
    return zero[-1] if zero else None


def _accepted_alphas(attempt: Mapping[str, Any] | None) -> list[float]:
    if attempt is None:
        return []
    values: list[float] = []
    for alpha in attempt.get("line_search_alphas", []):
        if alpha is not None and math.isfinite(float(alpha)):
            values.append(float(alpha))
    if values:
        return values
    for row in attempt.get("newton_history", []):
        alpha = row.get("line_search_alpha")
        if alpha is not None and math.isfinite(float(alpha)):
            values.append(float(alpha))
    return values


def _line_search_backtrack_count(
    alphas: Sequence[float],
    *,
    shrink: float,
) -> int:
    if not alphas:
        return 0
    if not (0.0 < shrink < 1.0):
        return 0
    total = 0
    log_shrink = math.log(shrink)
    for alpha in alphas:
        if alpha <= 0.0:
            continue
        total += max(0, int(round(math.log(float(alpha)) / log_shrink)))
    return int(total)


def _original_residual_l2(
    *,
    row: Mapping[str, Any],
    crawl_config: c0.C0Config,
    dtype: torch.dtype,
) -> float:
    if not row.get("state_artifact"):
        return math.nan
    tau = float(row["target_tau"])
    branch, provider, grid, boundaries = c0._branch_context(tau=tau, config=crawl_config, dtype=dtype)
    state = c0._load_state_artifact(_resolve_existing_path(str(row["state_artifact"])), dtype=dtype)
    residual = patha_closed_branch_residual(
        state,
        grid,
        branch,
        eos_K=float(branch.continuation_K_values[-1]),
        boundaries=boundaries,
        s_sigma=provider,
    ).detach()
    return float(torch.linalg.vector_norm(residual).cpu().item())


def _per_tau_crawl_table(
    crawl_result: Mapping[str, Any],
    *,
    crawl_config: c0.C0Config,
    dtype: torch.dtype,
    b2c_linf_tolerance: float,
) -> list[dict[str, Any]]:
    branch = c0.c0_frozen_branch(
        tau=float(c0.C0Config().depth_sequence[0]),
        grid=crawl_config.grid,
        max_newton_iters=crawl_config.max_newton_iters_override,
    )
    shrink = float(branch.newton.line_search_shrink)
    rows: list[dict[str, Any]] = []
    for index, row in enumerate(crawl_result.get("tau_attempts", [])):
        zero = _final_zero_epsilon_attempt(row)
        alphas = _accepted_alphas(zero)
        linf_raw = row.get("final_original_residual_linf")
        linf = float(linf_raw) if linf_raw is not None else math.nan
        solver_converged = bool(zero.get("converged")) if zero is not None else bool(row.get("stage_converged"))
        try:
            l2 = _original_residual_l2(row=row, crawl_config=crawl_config, dtype=dtype)
        except Exception as exc:  # pragma: no cover - diagnostic reports failures.
            l2 = math.nan
            l2_error = f"{type(exc).__name__}: {exc}"
        else:
            l2_error = None
        rows.append(
            {
                "attempt_index": int(index),
                "tau": float(row["target_tau"]),
                "nominal_tau": float(row.get("nominal_target_tau", row["target_tau"])),
                "tau_backtrack_index": int(row.get("backtrack_index", 0)),
                "solver_converged": solver_converged,
                "b2c_linf_pass": bool(math.isfinite(linf) and linf <= b2c_linf_tolerance),
                "accepted_default_success": bool(
                    solver_converged and math.isfinite(linf) and linf <= b2c_linf_tolerance
                ),
                "linf": linf,
                "l2": l2,
                "l2_recompute_error": l2_error,
                "newton_iterations": int(zero.get("iterations", row.get("stage_iterations", 0)))
                if zero is not None
                else None,
                "line_search_backtracks": _line_search_backtrack_count(alphas, shrink=shrink),
                "smallest_alpha": float(min(alphas)) if alphas else None,
                "zero_epsilon_final": bool(
                    zero is not None
                    and float(zero.get("core_epsilon", zero.get("epsilon", math.nan))) == 0.0
                    and float(zero.get("k1_radius_epsilon", math.nan)) == 0.0
                ),
                "message": row.get("message"),
                "state_artifact": row.get("state_artifact"),
                "init_source": row.get("init", {}).get("source"),
            }
        )
    return rows


def _state_lanes(state: torch.Tensor, grid) -> dict[str, np.ndarray]:
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    assert mu is not None
    psi = torch.cat([fields.psi_real.reshape(-1), fields.psi_imag.reshape(-1)])
    a_lanes = torch.cat([fields.a0.reshape(-1), fields.ar.reshape(-1), fields.aw.reshape(-1)])
    return {
        "full": state.detach().cpu().numpy().astype(np.float64, copy=True),
        "psi": psi.detach().cpu().numpy().astype(np.float64, copy=True),
        "A": a_lanes.detach().cpu().numpy().astype(np.float64, copy=True),
        "R0": fields.r0.detach().cpu().numpy().astype(np.float64, copy=True),
        "mu": np.asarray([float(mu.detach().cpu().item())], dtype=np.float64),
    }


def _normalized_delta_over_tau(
    left: np.ndarray,
    right: np.ndarray,
    delta_tau: float,
) -> float:
    denom = max(float(np.linalg.norm(right)), np.finfo(np.float64).tiny)
    return float(np.linalg.norm((right - left) / delta_tau) / denom)


def _fold_growth_from_values(
    values: Sequence[float],
    *,
    threshold: float,
) -> dict[str, Any]:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if len(finite) < 2:
        return {
            "status": "NOT_MEASURED",
            "reason": "fewer_than_two_converged_intervals",
            "growth_condition": False,
            "growth_factor": math.nan,
            "monotone": False,
        }
    window = finite[-min(3, len(finite)) :]
    monotone = all(left <= right for left, right in zip(window, window[1:]))
    growth_factor = float(window[-1] / max(window[0], np.finfo(np.float64).tiny))
    return {
        "status": "MEASURED",
        "reason": None,
        "growth_condition": bool(monotone and growth_factor >= float(threshold)),
        "growth_factor": growth_factor,
        "monotone": bool(monotone),
        "window": window,
        "threshold": float(threshold),
    }


def _smaller_tau_failed_all_backtracks(
    crawl_rows: Sequence[Mapping[str, Any]],
    per_tau_rows: Sequence[Mapping[str, Any]],
    *,
    crawl_config: c0.C0Config,
) -> dict[str, Any]:
    successful = [row for row in per_tau_rows if row.get("accepted_default_success")]
    if not successful:
        return {
            "failed_all_backtracks": False,
            "reason": "no_converged_reference_tau",
            "candidate_nominal_tau": None,
        }
    deepest_success = min(float(row["tau"]) for row in successful)
    failures = [
        dict(row)
        for row in crawl_rows
        if float(row.get("target_tau", math.inf)) < deepest_success - 5.0e-13
        and not bool(row.get("final_physical_converged"))
    ]
    if not failures:
        return {
            "failed_all_backtracks": False,
            "reason": "no_failed_next_smaller_tau_after_deepest_converged",
            "candidate_nominal_tau": None,
        }
    by_nominal: dict[float, list[dict[str, Any]]] = {}
    for row in failures:
        by_nominal.setdefault(float(row.get("nominal_target_tau", row["target_tau"])), []).append(row)
    nominal = min(by_nominal)
    group = by_nominal[nominal]
    max_bt = max(int(row.get("backtrack_index", 0)) for row in group)
    min_delta = min(abs(float(row.get("delta_tau") or 0.0)) for row in group)
    exhausted = bool(
        max_bt >= int(crawl_config.max_tau_backtracks)
        or min_delta <= float(crawl_config.min_tau_backtrack_delta)
    )
    return {
        "failed_all_backtracks": exhausted,
        "reason": None if exhausted else "smaller_tau_failed_but_backtrack_exhaustion_not_observed",
        "candidate_nominal_tau": nominal,
        "attempted_target_taus": [float(row["target_tau"]) for row in group],
        "max_backtrack_index": int(max_bt),
        "max_tau_backtracks": int(crawl_config.max_tau_backtracks),
        "min_abs_delta_tau": float(min_delta),
        "min_tau_backtrack_delta": float(crawl_config.min_tau_backtrack_delta),
    }


def _compute_fold_diagnostic(
    crawl_result: Mapping[str, Any],
    per_tau_rows: Sequence[Mapping[str, Any]],
    *,
    crawl_config: c0.C0Config,
    dtype: torch.dtype,
    growth_threshold: float,
) -> dict[str, Any]:
    source_rows = list(crawl_result.get("tau_attempts", []))
    successful_by_index = {
        int(row["attempt_index"]): row
        for row in per_tau_rows
        if row.get("accepted_default_success")
    }
    ordered = [
        source_rows[index]
        for index in range(len(source_rows))
        if index in successful_by_index
    ]
    if len(ordered) < 2:
        return {
            "status": "NOT_MEASURED",
            "call": "DIAGNOSTIC_INCOMPLETE",
            "reason": "fewer_than_two_consecutive_converged_intervals",
            "intervals": [],
        }

    loaded: list[tuple[dict[str, Any], Any, dict[str, np.ndarray]]] = []
    for row in ordered:
        tau = float(row["target_tau"])
        _branch, _provider, grid, _boundaries = c0._branch_context(
            tau=tau,
            config=crawl_config,
            dtype=dtype,
        )
        state = c0._load_state_artifact(_resolve_existing_path(str(row["state_artifact"])), dtype=dtype)
        loaded.append((dict(row), grid, _state_lanes(state, grid)))

    intervals: list[dict[str, Any]] = []
    for (left_row, _left_grid, left), (right_row, _right_grid, right) in zip(loaded, loaded[1:]):
        left_tau = float(left_row["target_tau"])
        right_tau = float(right_row["target_tau"])
        delta_tau = right_tau - left_tau
        if abs(delta_tau) <= np.finfo(np.float64).tiny:
            continue
        lane_breakdown = {
            lane: _normalized_delta_over_tau(left[lane], right[lane], delta_tau)
            for lane in ("psi", "R0", "mu")
        }
        lane_breakdown["A"] = _normalized_delta_over_tau(left["A"], right["A"], delta_tau)
        intervals.append(
            {
                "from_tau": left_tau,
                "to_tau": right_tau,
                "delta_tau": float(delta_tau),
                "normalized_full": _normalized_delta_over_tau(
                    left["full"],
                    right["full"],
                    delta_tau,
                ),
                "lane_breakdown": lane_breakdown,
            }
        )

    growth = _fold_growth_from_values(
        [row["normalized_full"] for row in intervals],
        threshold=growth_threshold,
    )
    failure = _smaller_tau_failed_all_backtracks(
        source_rows,
        per_tau_rows,
        crawl_config=crawl_config,
    )
    if growth["status"] != "MEASURED":
        call = "DIAGNOSTIC_INCOMPLETE"
        reason = growth.get("reason")
    elif growth["growth_condition"] and failure["failed_all_backtracks"]:
        call = "FOLD_DETECTED"
        reason = None
    elif growth["growth_condition"]:
        call = "FOLD_RISK"
        reason = "growth_condition_without_all-backtracks_failure"
    else:
        call = "FOLD_RISK"
        reason = "growth_condition_false"
    return {
        "status": "MEASURED" if intervals else "NOT_MEASURED",
        "call": call,
        "reason": reason,
        "growth_factor": growth.get("growth_factor"),
        "growth_threshold": float(growth_threshold),
        "monotone_growth_window": growth.get("monotone"),
        "growth_condition": bool(growth.get("growth_condition")),
        "smaller_tau_failure_condition": failure,
        "intervals": intervals,
    }


def _default_config_report(crawl_config: c0.C0Config) -> dict[str, Any]:
    default_config = c0.C0Config()
    branch = c0.c0_frozen_branch(tau=default_config.depth_sequence[0], grid=default_config.grid)
    newton = branch.newton
    return {
        "source": "stage1_solver.patha_c0_conditioning_spike.C0Config plus _c0_newton_config",
        "actual_crawl_paths": {
            "run_root": str(crawl_config.run_root),
            "json_path": str(crawl_config.json_path),
        },
        "code_default_knobs": {
            "grid": list(default_config.grid),
            "depth_sequence": list(default_config.depth_sequence),
            "max_tau_backtracks": int(default_config.max_tau_backtracks),
            "min_tau_backtrack_delta": float(default_config.min_tau_backtrack_delta),
            "max_newton_iters_override": default_config.max_newton_iters_override,
            "max_newton_iters": int(newton.max_newton_iters),
            "residual_norm": newton.residual_norm,
            "residual_atol": float(newton.residual_atol),
            "residual_rtol": float(newton.residual_rtol),
            "line_search": newton.line_search,
            "line_search_c1": float(newton.line_search_c1),
            "line_search_shrink": float(newton.line_search_shrink),
            "max_line_search_iters": int(newton.max_line_search_iters),
            "accept_best_line_search_decrease": bool(newton.accept_best_line_search_decrease),
            "gmres_rtol": float(newton.gmres_rtol),
            "gmres_atol": float(newton.gmres_atol),
            "gmres_restart": int(newton.gmres_restart),
            "gmres_maxiter": int(newton.gmres_maxiter),
            "epsilon_schedule": [asdict(aid) for aid in default_config.epsilon_schedule],
            "eos_final_only_after_first": bool(default_config.eos_final_only_after_first),
            "continuation_K_values_at_tau_0p03": list(branch.continuation_K_values),
            "use_wall_predictor": bool(default_config.use_wall_predictor),
            "prefer_existing_b2c_background_predictor": bool(
                default_config.prefer_existing_b2c_background_predictor
            ),
        },
        "path_only_aid_confirmation": {
            "physical_residual_entry_point": (
                "stage1_solver.coupled_branch.patha_closed_branch_residual"
            ),
            "preconditioner_residual_entry_point": (
                "stage1_solver.patha_c0_conditioning_spike.c0_preconditioner_residual"
            ),
            "newton_residual_merit_line_search_convergence": "original physical residual",
            "jvp_residual": "original physical residual",
            "epsilon_aids_scope": "preconditioner residual path only",
            "jacobi_scaling_scope": "linear system row/column scaling only; recovered step checked by original residual",
            "wall_predictor_scope": "initial guess for tau continuation only",
            "eos_continuation_scope": "path continuation; final accepted solve uses final EOS K",
            "final_accepted_solve": "zero core_epsilon and zero k1_radius_epsilon",
        },
    }


def _crawl_result_payload(
    *,
    config: c0.C0Config,
    sequence_rows: Sequence[Mapping[str, Any]],
    started: float,
    phase_status: str,
) -> dict[str, Any]:
    deepest_rows = [row for row in sequence_rows if row.get("final_physical_converged")]
    deepest_converged_row = (
        min(deepest_rows, key=lambda row: float(row["target_tau"])) if deepest_rows else None
    )
    schema_errors = c0._validate_tau_attempts_schema(sequence_rows)
    return {
        "schema": "stage1_pathA_C0f_default_crawl/v1",
        "source_revision": source_revision(),
        "phase_status": phase_status,
        "verdict": "DIAGNOSTIC_INCOMPLETE",
        "verdict_support": {
            "reason": "c0f_diagnostics_phase_not_run",
        },
        "prior_tau_floor": c0.PRIOR_TAU_FLOOR,
        "b2c_background_tolerance": c0.BACKGROUND_RESIDUAL_TOL,
        "deepest_converged_tau": None
        if deepest_converged_row is None
        else float(deepest_converged_row["target_tau"]),
        "deepest_converged_R0_min": None
        if deepest_converged_row is None
        else deepest_converged_row["metrics"]["min_R0"],
        "elapsed_seconds": time.perf_counter() - started,
        "config": c0._config_to_dict(config),
        "tau_attempts": [dict(row) for row in sequence_rows],
        "depth_sequence": [dict(row) for row in sequence_rows],
        "tau_attempt_schema_errors": schema_errors,
        "aids_admissibility": {
            "status": "NOT_MEASURED",
            "reason": "c0f_uses_default_crawl_only_for_C0f-0",
        },
        "sigma_diagnostics": {
            "status": "NOT_MEASURED",
            "reason": "not_part_of_C0f-0_default_crawl",
        },
        "fold_diagnostic": {
            "status": "NOT_MEASURED",
            "reason": "c0f_numeric_fold_detector_runs_after_crawl",
        },
        "r0_tau_diagnostic": c0._diagnose_tau_vs_r0(sequence_rows),
        "faithful_operator_boundary": c0._faithful_operator_boundary(),
        "recommended_next_step": "Run C0f reporting diagnostics on the saved default-crawl states.",
    }


def _write_crawl_checkpoint(
    *,
    checkpoint_path: Path,
    config: c0.C0Config,
    sequence_rows: Sequence[Mapping[str, Any]],
    started: float,
    phase_status: str,
) -> None:
    payload = _crawl_result_payload(
        config=config,
        sequence_rows=sequence_rows,
        started=started,
        phase_status=phase_status,
    )
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )


def _resume_state_from_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    crawl_config: c0.C0Config,
    dtype: torch.dtype,
) -> tuple[torch.Tensor | None, Any | None, float | None, int]:
    converged = [row for row in rows if row.get("final_physical_converged")]
    if not converged:
        return None, None, None, 0
    deepest = min(converged, key=lambda row: float(row["target_tau"]))
    tau = float(deepest["target_tau"])
    _branch, _provider, grid, _boundaries = c0._branch_context(
        tau=tau,
        config=crawl_config,
        dtype=dtype,
    )
    state = c0._load_state_artifact(_resolve_existing_path(str(deepest["state_artifact"])), dtype=dtype)
    floor_failures = sum(
        1
        for row in rows
        if float(row.get("target_tau", math.inf)) <= c0.PRIOR_TAU_FLOOR + 5.0e-13
        and not row.get("final_physical_converged")
    )
    return state.detach().clone(), grid, tau, int(floor_failures)


def _timeout_tau_attempt_row(
    *,
    target_tau: float,
    nominal_target_tau: float,
    delta_tau: float | None,
    backtrack_index: int,
    attempt_index: int,
    previous_tau: float | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    return {
        "tau": float(target_tau),
        "target_tau": float(target_tau),
        "nominal_target_tau": float(nominal_target_tau),
        "delta_tau": None if delta_tau is None else float(delta_tau),
        "backtrack_index": int(backtrack_index),
        "attempt_index": int(attempt_index),
        "start_from_tau": previous_tau,
        "target_relative_to_prior_floor": float(target_tau) / c0.PRIOR_TAU_FLOOR,
        "below_prior_floor": c0._below_prior_floor(float(target_tau)),
        "stage_converged": False,
        "final_original_residual_linf": None,
        "final_physical_converged": False,
        "b2c_background_tolerance": c0.BACKGROUND_RESIDUAL_TOL,
        "message": (
            "NOT_MEASURED: default tau attempt exceeded internal "
            f"{float(timeout_seconds):.1f}s per-script budget"
        ),
        "elapsed_seconds": float(timeout_seconds),
        "init": {
            "source": "previous_c0_converged_state" if previous_tau is not None else "default_closed_branch_state",
            "previous_tau": previous_tau,
        },
        "initialization": {
            "source": "previous_c0_converged_state" if previous_tau is not None else "default_closed_branch_state",
            "previous_tau": previous_tau,
        },
        "used_existing_b2c": False,
        "continuation_K_values": [],
        "epsilon_attempts": [],
        "substeps": [],
        "metrics": {},
        "state_artifact": None,
        "linear_diagnostics": {
            "status": "NOT_MEASURED",
            "reason": "tau_attempt_timeout",
        },
        "observables": {
            "available": False,
            "reason": "tau_attempt_timeout",
        },
        "floor_activation": {
            "core_epsilon_schedule": [aid.core_epsilon for aid in c0.C0Config().epsilon_schedule],
            "k1_radius_epsilon_schedule": [
                aid.k1_radius_epsilon for aid in c0.C0Config().epsilon_schedule
            ],
            "final_aids_inactive": None,
        },
        "not_measured": True,
        "timeout_seconds": float(timeout_seconds),
    }


def _execute_tau_attempt_with_timeout(
    *,
    timeout_seconds: float,
    **kwargs: Any,
) -> tuple[dict[str, Any], torch.Tensor | None, Any | None, bool]:
    def _handle_timeout(_signum: int, _frame: Any) -> None:
        raise _AttemptTimeout()

    previous = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _handle_timeout)
    signal.setitimer(signal.ITIMER_REAL, max(1.0, float(timeout_seconds)))
    try:
        row, state, grid = c0._execute_tau_attempt(**kwargs)
        return row, state, grid, False
    except _AttemptTimeout:
        row = _timeout_tau_attempt_row(
            target_tau=float(kwargs["target_tau"]),
            nominal_target_tau=float(kwargs["nominal_target_tau"]),
            delta_tau=kwargs.get("delta_tau"),
            backtrack_index=int(kwargs["backtrack_index"]),
            attempt_index=int(kwargs["attempt_index"]),
            previous_tau=kwargs.get("previous_tau"),
            timeout_seconds=float(timeout_seconds),
        )
        return row, None, None, True
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, previous)


def _run_default_crawl_resumable(
    config: c0.C0Config,
    *,
    dtype: torch.dtype,
    timeout_seconds: float,
) -> dict[str, Any]:
    if config.json_path.exists():
        return json.loads(config.json_path.read_text(encoding="utf-8"))

    started = time.perf_counter()
    checkpoint_path = config.run_root / "pathA_C0f_default_crawl_checkpoint.json"
    if checkpoint_path.exists():
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        sequence_rows = [dict(row) for row in checkpoint.get("tau_attempts", [])]
    else:
        sequence_rows = []

    last_converged_state, last_converged_grid, deepest_converged_tau, floor_or_below_failures = (
        _resume_state_from_rows(sequence_rows, crawl_config=config, dtype=dtype)
    )
    attempt_index = len(sequence_rows)
    timed_out = False

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
            row, final_state, final_grid, timed_out = _execute_tau_attempt_with_timeout(
                timeout_seconds=timeout_seconds,
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
            _write_crawl_checkpoint(
                checkpoint_path=checkpoint_path,
                config=config,
                sequence_rows=sequence_rows,
                started=started,
                phase_status="crawl_in_progress",
            )
            if timed_out:
                floor_or_below_failures += int(float(row["target_tau"]) <= c0.PRIOR_TAU_FLOOR + 5.0e-13)
                break
            if row["final_physical_converged"]:
                assert final_state is not None
                last_converged_state = final_state.detach().clone()
                last_converged_grid = final_grid
                deepest_converged_tau = float(row["target_tau"])
                break

            if float(row["target_tau"]) <= c0.PRIOR_TAU_FLOOR + 5.0e-13:
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

        if timed_out:
            break
        if floor_or_below_failures >= config.max_depth_failures_after_floor:
            break

    result = _crawl_result_payload(
        config=config,
        sequence_rows=sequence_rows,
        started=started,
        phase_status="crawl_incomplete_timeout" if timed_out else "crawl_complete",
    )
    config.json_path.parent.mkdir(parents=True, exist_ok=True)
    config.json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    checkpoint_path.write_text(
        json.dumps(result, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    return result


def _accepted_through_target(
    per_tau_rows: Sequence[Mapping[str, Any]],
    *,
    target_tau: float,
) -> bool:
    return any(
        row.get("accepted_default_success") and float(row["tau"]) <= float(target_tau) + 5.0e-13
        for row in per_tau_rows
    )


def _deepest_linf_pass_tau(per_tau_rows: Sequence[Mapping[str, Any]]) -> float | None:
    passed = [float(row["tau"]) for row in per_tau_rows if row.get("b2c_linf_pass")]
    return min(passed) if passed else None


def _deepest_accepted_tau(per_tau_rows: Sequence[Mapping[str, Any]]) -> float | None:
    accepted = [float(row["tau"]) for row in per_tau_rows if row.get("accepted_default_success")]
    return min(accepted) if accepted else None


def _first_stalled_row(
    crawl_result: Mapping[str, Any],
    per_tau_rows: Sequence[Mapping[str, Any]],
    *,
    target_tau: float,
) -> dict[str, Any] | None:
    if _accepted_through_target(per_tau_rows, target_tau=target_tau):
        return None
    by_index = {int(row["attempt_index"]): row for row in per_tau_rows}
    for index, source in enumerate(crawl_result.get("tau_attempts", [])):
        summary = by_index.get(index, {})
        if not summary.get("accepted_default_success", False):
            row = dict(source)
            row["_attempt_index"] = index
            row["_selection_reason"] = "first_default_crawl_attempt_without_solver_converged_and_linf_pass"
            return row
    return None


def _assemble_unscaled_original_matrix(
    *,
    state: torch.Tensor,
    grid,
    branch,
    provider,
    boundaries,
    matrix_path: Path,
) -> tuple[np.ndarray, np.ndarray, Any]:
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

    residual = residual_fn(state).detach().cpu().numpy().astype(np.float64)
    aid = c0.C0AidParameters(core_epsilon=0.0, k1_radius_epsilon=0.0)
    _built, _row_scale, _col_scale, matrix, _scaled = c0.build_scaled_preconditioner(
        original_residual_fn=residual_fn,
        preconditioner_residual_fn=residual_fn,
        x=state.detach(),
        rhs=-residual,
        grid=grid,
        newton=branch.newton,
        aid=aid,
        iteration=0,
    )
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    save_npz(matrix_path, matrix)
    return residual, matrix.toarray().astype(np.float64, copy=False), residual_fn


def _run_merit_sweep(
    stalled_row: Mapping[str, Any] | None,
    *,
    crawl_config: c0.C0Config,
    run_root: Path,
    dtype: torch.dtype,
    alpha_min_power: int,
) -> dict[str, Any]:
    if stalled_row is None:
        return {
            "status": "SKIPPED",
            "reason": "default_crawl_reached_required_target_or_no_stalled_row",
        }
    if not stalled_row.get("state_artifact"):
        return {
            "status": "NOT_MEASURED",
            "reason": "first_stalled_or_incomplete_tau_has_no_saved_state",
            "tau": float(stalled_row.get("target_tau", math.nan)),
            "selection_reason": stalled_row.get("_selection_reason"),
        }
    tau = float(stalled_row["target_tau"])
    try:
        branch, provider, grid, boundaries = c0._branch_context(tau=tau, config=crawl_config, dtype=dtype)
        state = c0._load_state_artifact(
            _resolve_existing_path(str(stalled_row["state_artifact"])),
            dtype=dtype,
        )
        matrix_path = run_root / "matrices" / (
            f"c0f_merit_tau_{c0._format_tau(tau)}_attempt_{int(stalled_row.get('_attempt_index', 0)):03d}.npz"
        )
        residual, dense, residual_fn = _assemble_unscaled_original_matrix(
            state=state,
            grid=grid,
            branch=branch,
            provider=provider,
            boundaries=boundaries,
            matrix_path=matrix_path,
        )
        initial_l2 = float(np.linalg.norm(residual))
        initial_linf = float(np.max(np.abs(residual)))
        try:
            step = np.linalg.solve(dense, -residual)
            solve_method = "dense_np_linalg_solve"
        except np.linalg.LinAlgError:
            step, _resid, _rank, _singular = np.linalg.lstsq(dense, -residual, rcond=None)
            solve_method = "dense_np_linalg_lstsq_fallback"
        linear_residual = dense @ step + residual
        linear_rel_resid = float(
            np.linalg.norm(linear_residual) / max(initial_l2, np.finfo(np.float64).tiny)
        )
        j_step = dense @ step
        step_t = torch.as_tensor(step, dtype=state.dtype, device=state.device)
        rows = []
        for power in range(int(alpha_min_power) + 1):
            alpha = float(2.0 ** (-power))
            predicted = residual + alpha * j_step
            trial = state.detach() + alpha * step_t
            actual = residual_fn(trial).detach().cpu().numpy().astype(np.float64)
            actual_l2 = float(np.linalg.norm(actual))
            predicted_l2 = float(np.linalg.norm(predicted))
            actual_linf = float(np.max(np.abs(actual)))
            predicted_linf = float(np.max(np.abs(predicted)))
            rows.append(
                {
                    "alpha": alpha,
                    "power": int(power),
                    "actual_l2": actual_l2,
                    "actual_linf": actual_linf,
                    "predicted_l2": predicted_l2,
                    "predicted_linf": predicted_linf,
                    "actual_l2_ratio": actual_l2 / max(initial_l2, np.finfo(np.float64).tiny),
                    "predicted_l2_ratio": predicted_l2 / max(initial_l2, np.finfo(np.float64).tiny),
                    "reduces_true_l2": bool(actual_l2 < initial_l2),
                    "finite": bool(np.all(np.isfinite(actual))),
                }
            )
        reducing = [row for row in rows if row["reduces_true_l2"]]
        best = min(rows, key=lambda row: row["actual_l2"]) if rows else None
        return {
            "status": "MEASURED",
            "tau": tau,
            "selection_reason": stalled_row.get("_selection_reason"),
            "state_artifact": stalled_row.get("state_artifact"),
            "matrix_path": str(matrix_path),
            "solve_method": solve_method,
            "initial_l2": initial_l2,
            "initial_linf": initial_linf,
            "step_norm": float(np.linalg.norm(step)),
            "linear_rel_resid": linear_rel_resid,
            "alpha_rows": rows,
            "any_alpha_reduces_true_l2": bool(reducing),
            "smallest_reducing_alpha": float(min(row["alpha"] for row in reducing)) if reducing else None,
            "best_alpha": None if best is None else float(best["alpha"]),
            "best_actual_l2": None if best is None else float(best["actual_l2"]),
            "diagnostic_interpretation": (
                "LOCALIZED_GLOBALIZATION_EVIDENCE"
                if reducing
                else "NO_ALPHA_REDUCED_TRUE_RESIDUAL"
            ),
        }
    except Exception as exc:  # pragma: no cover - diagnostic reports failures.
        return {
            "status": "NOT_MEASURED",
            "reason": f"{type(exc).__name__}: {exc}",
            "tau": tau,
        }


def _classify_transverse_fraction(value: float | None, config: C0fConfig) -> str:
    if value is None or not math.isfinite(float(value)):
        return "NOT_MEASURED"
    if float(value) <= float(config.gauge_fraction_threshold):
        return "GAUGE_BY_1_MINUS_P"
    if float(value) >= float(config.transverse_fraction_threshold):
        return "TRANSVERSE_CANDIDATE_BY_1_MINUS_P"
    return "MIXED_BY_1_MINUS_P"


def _unit(values: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(values))
    if norm <= np.finfo(np.float64).tiny:
        return np.asarray(values, dtype=np.float64)
    return np.asarray(values, dtype=np.float64) / norm


def _orthogonalize_to_basis(values: np.ndarray, basis: np.ndarray) -> np.ndarray:
    vector = np.asarray(values, dtype=np.float64)
    if basis.shape[1] == 0:
        return vector
    return vector - basis @ (basis.T @ vector)


def _mixed_control_rows(
    *,
    a_subspace: Mapping[str, Any],
    coupled_subspace: Mapping[str, Any],
    state: torch.Tensor,
    grid,
    branch,
    config: C0fConfig,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    r_scaled = grid.r_centers[:, None] / float(grid.spec.r_max)
    w_span = float(grid.spec.w_max) - float(grid.spec.w_min)
    w_scaled = (grid.w_centers[None, :] - float(grid.spec.w_min)) / w_span
    max_k = max(1, min(int(grid.spec.nr), int(grid.spec.nw)) // 2)
    k_values = [max(1, min(int(k), max_k)) for k in config.mixed_control_k_values]
    seen: set[int] = set()
    for k in k_values:
        if k in seen:
            continue
        seen.add(k)
        lambda_field = torch.sin(k * math.pi * r_scaled) * torch.cos((k + 1) * math.pi * w_scaled)
        chi = torch.sin((k + 1) * math.pi * r_scaled) * torch.sin(k * math.pi * w_scaled)
        gradient = _unit(c0._c0d_scalar_gradient_vector(lambda_field, grid))
        coupled_gradient = _unit(c0._c0e_coupled_gauge_vector(lambda_field, state, grid, branch))
        transverse = _unit(
            _orthogonalize_to_basis(
                c0._c0e_stream_function_a_vector(chi, grid),
                np.asarray(a_subspace["basis"], dtype=np.float64),
            )
        )
        eps = float(config.mixed_control_epsilon)
        mixed_a = _unit(gradient + eps * transverse)
        mixed_coupled = _unit(coupled_gradient + eps * transverse)
        p_g = c0._c0d_projection_fraction(np.asarray(a_subspace["basis"], dtype=np.float64), mixed_a)
        p_cpl = c0._c0d_projection_fraction(
            np.asarray(coupled_subspace["basis"], dtype=np.float64),
            mixed_coupled,
        )
        curl = c0._c0e_curl_fraction_metrics(mixed_a, grid, branch)
        expected = float(eps * eps / (1.0 + eps * eps))
        rows.append(
            {
                "k": int(k),
                "epsilon": eps,
                "expected_orthogonal_fraction": expected,
                "one_minus_p_g": float(max(0.0, 1.0 - p_g)) if math.isfinite(p_g) else math.nan,
                "one_minus_p_cpl": float(max(0.0, 1.0 - p_cpl)) if math.isfinite(p_cpl) else math.nan,
                "p_g": float(p_g),
                "p_cpl": float(p_cpl),
                "raw_curl_unweighted": curl["unweighted"]["curl_fraction"],
                "raw_curl_cell_volume_weighted": curl["cell_volume_weighted"]["curl_fraction"],
                "boundary_curl_fraction_weighted": curl["cell_volume_weighted"][
                    "boundary_curl_fraction_of_norm"
                ],
            }
        )
    return rows


def _run_gauge_reconfirmation(config: C0fConfig, *, dtype: torch.dtype) -> dict[str, Any]:
    try:
        c0b_json_path = _resolve_existing_path(config.c0b_json_path)
        if not c0b_json_path.exists():
            return {
                "status": "NOT_MEASURED",
                "reason": f"missing_c0b_json:{c0b_json_path}",
            }
        c0b_result = json.loads(c0b_json_path.read_text(encoding="utf-8"))
        c0e_config = c0.C0eConfig(
            c0b_json_path=c0b_json_path,
            run_all_available_tau_sources=False,
        )
        rows = c0._c0e_available_artifact_rows(c0b_result, c0e_config)
        if not rows:
            return {
                "status": "NOT_MEASURED",
                "reason": "missing_saved_state_or_matrix_with_measured_linear_diagnostics",
                "c0b_source_json": str(c0b_json_path),
            }
        artifact = c0._c0e_analyze_artifact(
            row=rows[0],
            c0b_result=c0b_result,
            dtype=dtype,
            config=c0e_config,
        )
        tau = float(artifact["tau"])
        c0b_grid = tuple(int(value) for value in c0b_result.get("config", {}).get("grid", c0.C0Config().grid))
        branch, _provider, grid, _boundaries, _eos_K, _residual_fn = c0._c0c_residual_context(
            tau=tau,
            grid_shape=c0b_grid,
            dtype=dtype,
        )
        state = c0._load_state_artifact(Path(artifact["state_artifact"]), dtype=dtype)
        a_subspace = c0._c0d_build_gauge_subspace(
            grid,
            branch,
            c0.C0dConfig(
                gradient_rank_rtol=c0e_config.gradient_rank_rtol,
                harmonic_weighted_divergence_rtol=c0e_config.harmonic_weighted_divergence_rtol,
            ),
        )
        coupled_subspace = c0._c0e_build_coupled_subspace(state, grid, branch, c0e_config)
        mode_rows = []
        for mode in artifact.get("maxwell_modes", []):
            p_g = float(mode.get("a_only_p_g_fraction", math.nan))
            p_g_a = float(mode.get("a_only_p_g_a_normalized_fraction", math.nan))
            p_cpl = float(mode.get("coupled_capture_fraction", math.nan))
            one_minus_p_g = max(0.0, 1.0 - p_g) if math.isfinite(p_g) else math.nan
            one_minus_p_g_a = max(0.0, 1.0 - p_g_a) if math.isfinite(p_g_a) else math.nan
            one_minus_p_cpl = max(0.0, 1.0 - p_cpl) if math.isfinite(p_cpl) else math.nan
            classifier_value = one_minus_p_g_a if math.isfinite(one_minus_p_g_a) else one_minus_p_g
            curls = mode.get("curl_fraction", {})
            weighted = curls.get("cell_volume_weighted", {})
            mode_rows.append(
                {
                    "mode_index": int(mode["mode_index"]),
                    "sigma": float(mode["sigma"]),
                    "one_minus_p_g": one_minus_p_g,
                    "one_minus_p_g_a_normalized": one_minus_p_g_a,
                    "one_minus_p_cpl": one_minus_p_cpl,
                    "a_normalized_transverse_fraction": one_minus_p_g_a,
                    "classification": _classify_transverse_fraction(classifier_value, config),
                    "raw_curl_unweighted": curls.get("unweighted", {}).get("curl_fraction"),
                    "raw_curl_cell_volume_weighted": weighted.get("curl_fraction"),
                    "boundary_curl_fraction_weighted": weighted.get("boundary_curl_fraction_of_norm"),
                    "interior_curl_fraction_weighted": (
                        1.0 - float(weighted.get("boundary_curl_fraction_of_norm"))
                        if weighted.get("boundary_curl_fraction_of_norm") is not None
                        else math.nan
                    ),
                    "spatial_a_energy_fraction": mode.get("spatial_a_energy_fraction"),
                }
            )
        phase = artifact.get("phase_mode", {})
        phase_row = {
            "mode_index": int(phase.get("mode_index", 0)),
            "sigma": phase.get("sigma"),
            "one_minus_phase_capture": (
                max(0.0, 1.0 - float(phase.get("phase_capture_fraction", math.nan)))
                if phase.get("phase_capture_fraction") is not None
                else math.nan
            ),
            "one_minus_p_cpl": (
                max(0.0, 1.0 - float(phase.get("coupled_capture_fraction", math.nan)))
                if phase.get("coupled_capture_fraction") is not None
                else math.nan
            ),
            "one_minus_p_g": None,
            "classification": "PHASE_GAUGE_BY_PHASE_AND_COUPLED_CAPTURE",
            "raw_curl_reported_only": {
                "unweighted_norm": phase.get("z_f_rw_norm_unweighted"),
                "cell_volume_weighted_norm": phase.get("z_f_rw_norm_cell_volume_weighted"),
            },
        }
        mixed_controls = _mixed_control_rows(
            a_subspace=a_subspace,
            coupled_subspace=coupled_subspace,
            state=state,
            grid=grid,
            branch=branch,
            config=config,
        )
        return {
            "status": "MEASURED",
            "c0b_source_json": str(c0b_json_path),
            "tau": tau,
            "state_artifact": artifact.get("state_artifact"),
            "matrix_path": artifact.get("matrix_path"),
            "classifier": "dimensionless 1-P_G / 1-P_cpl; raw curl reported only",
            "thresholds": {
                "gauge_fraction_max": float(config.gauge_fraction_threshold),
                "transverse_candidate_fraction_min": float(config.transverse_fraction_threshold),
            },
            "phase_mode": phase_row,
            "maxwell_modes": mode_rows,
            "mixed_controls": mixed_controls,
            "old_c0e_labels_ignored": True,
            "raw_curl_used_as_classifier": False,
        }
    except Exception as exc:  # pragma: no cover - diagnostic reports failures.
        return {
            "status": "NOT_MEASURED",
            "reason": f"{type(exc).__name__}: {exc}",
        }


def _determine_c0f_verdict(
    *,
    per_tau_rows: Sequence[Mapping[str, Any]],
    fold: Mapping[str, Any],
    merit: Mapping[str, Any],
    config: C0fConfig,
) -> tuple[str, dict[str, Any], str]:
    accepted_target = _accepted_through_target(
        per_tau_rows,
        target_tau=config.wall_was_config_target_tau,
    )
    fold_call = str(fold.get("call"))
    fold_measured = fold.get("status") == "MEASURED"
    bounded = bool(
        fold_measured
        and not bool(fold.get("growth_condition"))
        and fold_call != "FOLD_DETECTED"
    )
    deepest_accepted = _deepest_accepted_tau(per_tau_rows)
    deepest_linf = _deepest_linf_pass_tau(per_tau_rows)
    if accepted_target and bounded:
        return (
            "WALL_WAS_CONFIG",
            {
                "accepted_default_crawl_through_tau": float(config.wall_was_config_target_tau),
                "deepest_accepted_tau": deepest_accepted,
                "deepest_linf_pass_tau": deepest_linf,
                "fold_call": fold_call,
                "fold_growth_factor": fold.get("growth_factor"),
                "bounded_by_numeric_rule": True,
            },
            "Continue the default crawl toward the physics-target depth; keep gauge deflation deferred as optional insurance.",
        )
    if fold_call == "FOLD_DETECTED":
        return (
            "FOLD_DETECTED",
            {
                "fold_growth_factor": fold.get("growth_factor"),
                "fold_growth_threshold": fold.get("growth_threshold"),
                "smaller_tau_failure_condition": fold.get("smaller_tau_failure_condition"),
            },
            "Gate C0g to pseudo-arclength continuation; do not implement it in C0f.",
        )
    if not accepted_target and merit.get("status") == "MEASURED":
        if not bool(merit.get("any_alpha_reduces_true_l2")) and bounded:
            return (
                "GLOBALIZATION_INSUFFICIENT",
                {
                    "accepted_default_crawl_through_target": False,
                    "deepest_accepted_tau": deepest_accepted,
                    "merit_any_alpha_reduces_true_l2": False,
                    "linear_rel_resid": merit.get("linear_rel_resid"),
                    "fold_call": fold_call,
                },
                "Localize the predicted-vs-actual gap before solver rebuild; C0g should target the reported blocker.",
            )
        return (
            "DIAGNOSTIC_INCOMPLETE",
            {
                "accepted_default_crawl_through_target": False,
                "deepest_accepted_tau": deepest_accepted,
                "merit_any_alpha_reduces_true_l2": merit.get("any_alpha_reduces_true_l2"),
                "merit_interpretation": merit.get("diagnostic_interpretation"),
                "fold_call": fold_call,
                "reason": "default crawl stalled but C0f substantive verdict gates do not fit",
            },
            "Do not declare the wall fixed; use the merit-sweep localization to define the gated C0g follow-up.",
        )
    return (
        "DIAGNOSTIC_INCOMPLETE",
        {
            "accepted_default_crawl_through_target": bool(accepted_target),
            "deepest_accepted_tau": deepest_accepted,
            "deepest_linf_pass_tau": deepest_linf,
            "fold_status": fold.get("status"),
            "fold_call": fold_call,
            "merit_status": merit.get("status"),
            "reason": "required fold or merit evidence not measured for a substantive verdict",
        },
        "Rerun the missing bounded diagnostic evidence before selecting C0g.",
    )


def _write_report(result: Mapping[str, Any], path: Path) -> None:
    lines: list[str] = []
    lines.append("# Path-A C0f Globalization Probe")
    lines.append("")
    lines.append(f"C0f-3 verdict: **{result.get('verdict')}**")
    lines.append("")
    lines.append("## Default Config")
    lines.append("")
    defaults = result.get("default_config", {})
    knobs = defaults.get("code_default_knobs", {})
    lines.append("```yaml")
    for key in (
        "grid",
        "depth_sequence",
        "max_newton_iters",
        "max_newton_iters_override",
        "max_tau_backtracks",
        "line_search",
        "max_line_search_iters",
        "line_search_shrink",
        "epsilon_schedule",
        "use_wall_predictor",
        "eos_final_only_after_first",
        "continuation_K_values_at_tau_0p03",
    ):
        lines.append(f"{key}: {knobs.get(key)}")
    lines.append("path_only_aid_confirmation:")
    for key, value in defaults.get("path_only_aid_confirmation", {}).items():
        lines.append(f"  {key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## C0f-0 Default Crawl")
    lines.append("")
    crawl_rows = [
        {
            "tau": row.get("tau"),
            "solver_converged": row.get("solver_converged"),
            "b2c_linf_pass": row.get("b2c_linf_pass"),
            "Linf": row.get("linf"),
            "L2": row.get("l2"),
            "iters": row.get("newton_iterations"),
            "backtracks": row.get("line_search_backtracks"),
            "alpha": row.get("smallest_alpha"),
            "tau_bt": row.get("tau_backtrack_index"),
        }
        for row in result.get("per_tau_crawl", [])
    ]
    lines.append(
        _markdown_table(
            [
                "tau",
                "solver_converged",
                "b2c_linf_pass",
                "Linf",
                "L2",
                "iters",
                "backtracks",
                "alpha",
                "tau_bt",
            ],
            crawl_rows,
        )
    )
    lines.append("")
    lines.append("## C0f-0a Numeric Fold Detector")
    lines.append("")
    fold = result.get("fold_diagnostic", {})
    lines.append("```yaml")
    for key in (
        "status",
        "call",
        "reason",
        "growth_factor",
        "growth_threshold",
        "monotone_growth_window",
        "growth_condition",
        "smaller_tau_failure_condition",
    ):
        lines.append(f"{key}: {fold.get(key)}")
    lines.append("```")
    fold_rows = []
    for row in fold.get("intervals", []):
        lanes = row.get("lane_breakdown", {})
        fold_rows.append(
            {
                "from": row.get("from_tau"),
                "to": row.get("to_tau"),
                "full": row.get("normalized_full"),
                "psi": lanes.get("psi"),
                "R0": lanes.get("R0"),
                "mu": lanes.get("mu"),
                "A": lanes.get("A"),
            }
        )
    lines.append(_markdown_table(["from", "to", "full", "psi", "R0", "mu", "A"], fold_rows))
    lines.append("")
    lines.append("## C0f-1 Merit Sweep")
    lines.append("")
    merit = result.get("merit_sweep", {})
    lines.append("```yaml")
    for key in (
        "status",
        "reason",
        "tau",
        "initial_l2",
        "initial_linf",
        "linear_rel_resid",
        "any_alpha_reduces_true_l2",
        "smallest_reducing_alpha",
        "diagnostic_interpretation",
    ):
        lines.append(f"{key}: {merit.get(key)}")
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
        lines.append(
            _markdown_table(
                ["alpha", "actual_l2", "pred_l2", "actual_linf", "pred_linf", "reduces"],
                merit_rows,
            )
        )
    lines.append("")
    lines.append("## C0f-2 1-P_G Gauge Re-Confirm")
    lines.append("")
    gauge = result.get("gauge_reconfirm", {})
    lines.append("```yaml")
    for key in (
        "status",
        "reason",
        "tau",
        "classifier",
        "thresholds",
        "old_c0e_labels_ignored",
        "raw_curl_used_as_classifier",
        "state_artifact",
        "matrix_path",
    ):
        lines.append(f"{key}: {gauge.get(key)}")
    lines.append("phase_mode: " + str(gauge.get("phase_mode")))
    lines.append("```")
    mode_rows = []
    for mode in gauge.get("maxwell_modes", []):
        mode_rows.append(
            {
                "mode": mode.get("mode_index"),
                "sigma": mode.get("sigma"),
                "1-P_G": mode.get("one_minus_p_g"),
                "1-P_G_A": mode.get("one_minus_p_g_a_normalized"),
                "1-P_cpl": mode.get("one_minus_p_cpl"),
                "boundary_curl": mode.get("boundary_curl_fraction_weighted"),
                "class": mode.get("classification"),
            }
        )
    lines.append(
        _markdown_table(
            ["mode", "sigma", "1-P_G", "1-P_G_A", "1-P_cpl", "boundary_curl", "class"],
            mode_rows,
        )
    )
    lines.append("")
    control_rows = [
        {
            "k": row.get("k"),
            "eps": row.get("epsilon"),
            "expected": row.get("expected_orthogonal_fraction"),
            "1-P_G": row.get("one_minus_p_g"),
            "1-P_cpl": row.get("one_minus_p_cpl"),
            "curl_w": row.get("raw_curl_cell_volume_weighted"),
            "boundary_curl": row.get("boundary_curl_fraction_weighted"),
        }
        for row in gauge.get("mixed_controls", [])
    ]
    lines.append("Mixed controls (exact gradient + epsilon transverse):")
    lines.append("")
    lines.append(
        _markdown_table(
            ["k", "eps", "expected", "1-P_G", "1-P_cpl", "curl_w", "boundary_curl"],
            control_rows,
        )
    )
    lines.append("")
    lines.append("## Verdict Support")
    lines.append("")
    lines.append("```yaml")
    for key, value in result.get("verdict_support", {}).items():
        lines.append(f"{key}: {value}")
    lines.append(f"recommended_next_step: {result.get('recommended_next_step')}")
    lines.append("```")
    lines.append("")
    lines.append("## Scope Guard")
    lines.append("")
    lines.append("```yaml")
    for key, value in result.get("scope_guard", {}).items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Verification")
    lines.append("")
    lines.append("```yaml")
    verification = result.get("verification", {})
    for key, value in verification.items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append(f"Machine artifact: `{result.get('json_path')}`")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _stdout_summary(result: Mapping[str, Any]) -> str:
    defaults = result.get("default_config", {}).get("code_default_knobs", {})
    lines: list[str] = []
    lines.append(
        "Default config: "
        f"max_newton_iters={defaults.get('max_newton_iters')}, "
        f"max_tau_backtracks={defaults.get('max_tau_backtracks')}, "
        f"depth_sequence={defaults.get('depth_sequence')}, "
        f"line_search={defaults.get('line_search')}({defaults.get('max_line_search_iters')}), "
        f"epsilon_schedule={defaults.get('epsilon_schedule')}; path-only aids, final zero-epsilon."
    )
    lines.append("Per-tau crawl:")
    for row in result.get("per_tau_crawl", []):
        lines.append(
            "  "
            f"tau={float(row['tau']):.12e} "
            f"solver_converged={row.get('solver_converged')} "
            f"b2c_linf_pass={row.get('b2c_linf_pass')} "
            f"Linf={_fmt(row.get('linf'))} "
            f"iters={row.get('newton_iterations')} "
            f"backtracks={row.get('line_search_backtracks')} "
            f"alpha={_fmt(row.get('smallest_alpha'))}"
        )
    lines.append(f"Deepest tau with Linf<=1e-6: {_fmt(result.get('deepest_linf_pass_tau'))}")
    fold = result.get("fold_diagnostic", {})
    lines.append(
        "Fold detector: "
        f"growth_factor={_fmt(fold.get('growth_factor'))}, "
        f"call={fold.get('call')}"
    )
    intervals = fold.get("intervals", [])
    if intervals:
        lanes = intervals[-1].get("lane_breakdown", {})
        lines.append(
            "  last lane breakdown: "
            f"psi={_fmt(lanes.get('psi'))}, R0={_fmt(lanes.get('R0'))}, "
            f"mu={_fmt(lanes.get('mu'))}, A={_fmt(lanes.get('A'))}"
        )
    merit = result.get("merit_sweep", {})
    if merit.get("status") == "MEASURED":
        lines.append(
            "Merit sweep: "
            f"linear_rel_resid={_fmt(merit.get('linear_rel_resid'))}, "
            f"any_alpha_reduces={merit.get('any_alpha_reduces_true_l2')}, "
            f"best_alpha={_fmt(merit.get('best_alpha'))}, "
            f"best_actual_l2={_fmt(merit.get('best_actual_l2'))}"
        )
    else:
        lines.append(f"Merit sweep: {merit.get('status')} ({merit.get('reason')})")
    gauge = result.get("gauge_reconfirm", {})
    if gauge.get("status") == "MEASURED":
        pieces = []
        for mode in gauge.get("maxwell_modes", []):
            pieces.append(
                f"mode {mode.get('mode_index')}:1-P_G_A="
                f"{_fmt(mode.get('one_minus_p_g_a_normalized'))}"
            )
        lines.append("Gauge re-confirm: " + ", ".join(pieces))
        phase = gauge.get("phase_mode", {})
        lines.append(
            "  phase/mode split: "
            f"phase 1-P_cpl={_fmt(phase.get('one_minus_p_cpl'))}; "
            "modes 1/3/4 vs mode 2 classified by 1-P_G, not raw curl."
        )
    else:
        lines.append(f"Gauge re-confirm: {gauge.get('status')} ({gauge.get('reason')})")
    lines.append(f"C0f-3 verdict: {result.get('verdict')}")
    lines.append(f"Recommended next step: {result.get('recommended_next_step')}")
    gates = result.get("verification", {}).get("chunk_gates", {})
    lines.append(
        "Guard confirmation: operators/frozen/export/SOLVER-LOGIC untouched; "
        "no new numerics, no deflation; chunk gates "
        f"{gates.get('status', 'NOT_RUN')}."
    )
    files = result.get("files_changed", [])
    if files:
        lines.append("Files changed: " + ", ".join(files))
    lines.append(f"Artifacts: report={result.get('report_path')}, json={result.get('json_path')}")
    return "\n".join(lines)


def run_c0f_globalization_probe(config: C0fConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = C0fConfig()
    started = time.perf_counter()
    config.run_root.mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(BackendConfig())
    crawl_config = c0.C0Config(
        run_root=config.run_root / "default_crawl",
        report_path=config.run_root / "default_crawl_report.md",
        json_path=config.run_root / "pathA_C0f_default_crawl.json",
    )
    crawl_result = _run_default_crawl_resumable(
        crawl_config,
        dtype=dtype,
        timeout_seconds=config.per_tau_attempt_timeout_seconds,
    )
    per_tau = _per_tau_crawl_table(
        crawl_result,
        crawl_config=crawl_config,
        dtype=dtype,
        b2c_linf_tolerance=config.b2c_linf_tolerance,
    )
    fold = _compute_fold_diagnostic(
        crawl_result,
        per_tau,
        crawl_config=crawl_config,
        dtype=dtype,
        growth_threshold=config.fold_growth_factor_threshold,
    )
    stalled_row = _first_stalled_row(
        crawl_result,
        per_tau,
        target_tau=config.wall_was_config_target_tau,
    )
    merit = _run_merit_sweep(
        stalled_row,
        crawl_config=crawl_config,
        run_root=config.run_root,
        dtype=dtype,
        alpha_min_power=config.merit_alpha_min_power,
    )
    gauge = _run_gauge_reconfirmation(config, dtype=dtype)
    verdict, support, next_step = _determine_c0f_verdict(
        per_tau_rows=per_tau,
        fold=fold,
        merit=merit,
        config=config,
    )
    result = {
        "schema": "stage1_pathA_C0f_globalization_probe/v1",
        "source_revision": source_revision(),
        "config": _config_to_dict(config),
        "default_config": _default_config_report(crawl_config),
        "crawl_json_path": str(crawl_config.json_path),
        "per_tau_crawl": per_tau,
        "deepest_linf_pass_tau": _deepest_linf_pass_tau(per_tau),
        "deepest_accepted_tau": _deepest_accepted_tau(per_tau),
        "fold_diagnostic": fold,
        "merit_sweep": merit,
        "gauge_reconfirm": gauge,
        "verdict": verdict,
        "verdict_support": support,
        "recommended_next_step": next_step,
        "scope_guard": {
            "existing_closed_newton_crawl_used": True,
            "solver_logic_changed_by_c0f": False,
            "operators_touched_by_c0f": False,
            "frozen_physics_touched_by_c0f": False,
            "physical_export_guard_touched_by_c0f": False,
            "trust_region_dogleg_lm_implemented": False,
            "pseudo_arclength_implemented": False,
            "gauge_deflation_implemented": False,
            "depth_continuation": "tau_only",
            "single_arbiter_residual": "stage1_solver.coupled_branch.patha_closed_branch_residual",
        },
        "verification": {
            "chunk_gates": {
                "status": "NOT_RUN",
                "reason": "external verification commands are recorded after script execution",
            }
        },
        "files_changed": [
            "src/stage1_solver/patha_c0f_globalization_probe.py",
            "scripts/pathA_C0f_globalization_probe.py",
            "tests/test_patha_c0f_globalization_probe.py",
            "reports/pathA_C0f_globalization_probe.md",
            "runs/pathA_C0f_globalization_probe/pathA_C0f_globalization_probe.json",
        ],
        "elapsed_seconds": time.perf_counter() - started,
        "report_path": str(config.report_path),
        "json_path": str(config.json_path),
    }
    config.json_path.parent.mkdir(parents=True, exist_ok=True)
    config.json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    _write_report(result, config.report_path)
    return result


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Path-A C0f globalization probe.")
    parser.add_argument("--run-root", type=Path, default=DEFAULT_C0F_RUN_ROOT)
    parser.add_argument("--report-path", type=Path, default=DEFAULT_C0F_REPORT_PATH)
    parser.add_argument("--json-path", type=Path, default=DEFAULT_C0F_JSON_PATH)
    parser.add_argument("--c0b-json-path", type=Path, default=C0fConfig().c0b_json_path)
    parser.add_argument(
        "--fold-growth-factor-threshold",
        type=float,
        default=C0fConfig().fold_growth_factor_threshold,
    )
    parser.add_argument(
        "--per-tau-attempt-timeout-seconds",
        type=float,
        default=C0fConfig().per_tau_attempt_timeout_seconds,
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    config = C0fConfig(
        run_root=args.run_root,
        report_path=args.report_path,
        json_path=args.json_path,
        c0b_json_path=args.c0b_json_path,
        fold_growth_factor_threshold=float(args.fold_growth_factor_threshold),
        per_tau_attempt_timeout_seconds=float(args.per_tau_attempt_timeout_seconds),
    )
    result = run_c0f_globalization_probe(config)
    print(_stdout_summary(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
