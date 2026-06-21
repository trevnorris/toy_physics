#!/usr/bin/env python3
"""Finalize C0g B-1/B-2 after bounded full-budget attempt timeouts."""

from __future__ import annotations

import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from stage1_solver import patha_c0_conditioning_spike as c0
from stage1_solver.backend import configure_backend
from stage1_solver.config import BackendConfig, source_revision


def _anchor_row(
    *,
    config: c0.C0gBuildB1B2Config,
    payload: dict,
    tau: float,
    attempt_index: int,
    dtype,
    start_from_tau: float | None,
) -> tuple[dict, object, object]:
    existing = c0._c0g_per_tau_row_for_sequence(payload, tau)
    if existing is None or not bool(existing.get("solver_converged")):
        raise RuntimeError(f"missing C0f2 anchor for tau={tau:.12e}")
    state = c0._c0g_state_vector_from_path(existing["state_artifact"], dtype=dtype)
    residual_linf = c0._c0g_saved_state_residual_linf(
        state=state,
        tau=tau,
        config=config,
        dtype=dtype,
    )
    branch, _provider, grid, _boundaries = c0._branch_context(
        tau=tau,
        config=c0.C0Config(grid=config.grid),
        dtype=dtype,
    )
    metrics = c0._state_metrics(state, grid)
    state_path = c0._resolve_output_path(
        c0._state_artifact_path(
            c0.C0Config(run_root=config.run_root / "gauge_fixed_crawl"),
            tau=tau,
            attempt_index=attempt_index,
        )
    )
    c0._save_state_artifact(state_path, state)
    row = {
        "tau": tau,
        "target_tau": tau,
        "nominal_target_tau": tau,
        "delta_tau": None if start_from_tau is None else tau - float(start_from_tau),
        "backtrack_index": 0,
        "start_from_tau": start_from_tau,
        "target_relative_to_prior_floor": tau / c0.PRIOR_TAU_FLOOR,
        "below_prior_floor": c0._below_prior_floor(tau),
        "stage_converged": residual_linf <= c0.BACKGROUND_RESIDUAL_TOL,
        "final_original_residual_linf": residual_linf,
        "final_physical_converged": residual_linf <= c0.BACKGROUND_RESIDUAL_TOL,
        "b2c_background_tolerance": c0.BACKGROUND_RESIDUAL_TOL,
        "message": "existing C0f2 genuine warm-start converged state accepted as measured anchor",
        "elapsed_seconds": 0.0,
        "init": {
            "source": "c0f2_genuine_warm_start_converged_state",
            "path": existing.get("state_artifact"),
            "original_init_source": existing.get("init_source"),
        },
        "initialization": {
            "source": "c0f2_genuine_warm_start_converged_state",
            "path": existing.get("state_artifact"),
        },
        "used_existing_b2c": False,
        "gauge_fix_solver_invoked": False,
        "max_newton_iters": 20,
        "max_newton_iters_override": config.max_newton_iters_override,
        "max_tau_backtracks": int(config.max_tau_backtracks),
        "min_tau_backtrack_delta": float(config.min_tau_backtrack_delta),
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
        "measurement_status": "MEASURED",
    }
    return row, state.detach().clone(), grid


def main() -> int:
    started = time.perf_counter()
    config = c0.C0gBuildB1B2Config()
    dtype = configure_backend(BackendConfig())
    payload = c0._c0g_load_c0f2_payload(c0.C0gConfig(c0f2_json_path=config.c0f2_json_path))
    run_root = c0._resolve_output_path(config.run_root / "gauge_fixed_crawl")
    run_root.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    row, state, grid = _anchor_row(
        config=config,
        payload=payload,
        tau=0.03,
        attempt_index=0,
        dtype=dtype,
        start_from_tau=None,
    )
    rows.append(row)
    last_tau = 0.03
    row, state, grid = _anchor_row(
        config=config,
        payload=payload,
        tau=0.0295,
        attempt_index=1,
        dtype=dtype,
        start_from_tau=last_tau,
    )
    rows.append(row)
    last_tau = 0.0295

    proof_existing = c0._c0g_per_tau_row_for_sequence(payload, config.proof_tau)
    if proof_existing is None or not bool(proof_existing.get("solver_converged")):
        raise RuntimeError("missing proof tau C0f2 state")
    branch, _provider, proof_grid, _boundaries = c0._branch_context(
        tau=float(config.proof_tau),
        config=c0.C0Config(grid=config.grid),
        dtype=dtype,
    )
    proof_base = c0._c0g_state_vector_from_path(proof_existing["state_artifact"], dtype=dtype)
    proof_seed, perturbation = c0._c0g_apply_path_only_proof_perturbation(
        proof_base,
        proof_grid,
        phase_angle=float(config.proof_global_phase_radians),
        r0_amplitude=float(config.proof_r0_perturbation),
    )
    proof_row, proof_state, proof_grid = c0._c0g_single_gauge_fixed_attempt(
        x0=proof_seed,
        tau=float(config.proof_tau),
        config=config,
        dtype=dtype,
        attempt_index=2,
        start_from_tau=last_tau,
        init={
            "source": "c0f2_converged_state_non_gauge_perturbed_for_gauge_fixed_proof",
            "path": proof_existing.get("state_artifact"),
            "perturbation": perturbation,
        },
        max_newton_iters=config.max_newton_iters_override,
    )
    proof_row["proof_solve_role"] = "genuine_gauge_fixed_newton_solve"
    proof_row["proof_seed_description"] = (
        "C0f2 proof-tau state plus global phase and non-gauge r0 perturbation; not a pure global-phase orbit seed"
    )
    proof_row["measurement_status"] = "MEASURED"
    rows.append(proof_row)
    last_tau = float(config.proof_tau)

    row, state, grid = _anchor_row(
        config=config,
        payload=payload,
        tau=0.029125,
        attempt_index=3,
        dtype=dtype,
        start_from_tau=last_tau,
    )
    rows.append(row)
    last_tau = 0.029125

    timeout_rows = [
        {
            "tau": 0.029124,
            "target_tau": 0.029124,
            "nominal_target_tau": 0.029124,
            "delta_tau": 0.029124 - last_tau,
            "backtrack_index": 0,
            "start_from_tau": last_tau,
            "target_relative_to_prior_floor": 0.029124 / c0.PRIOR_TAU_FLOOR,
            "below_prior_floor": c0._below_prior_floor(0.029124),
            "stage_converged": False,
            "final_original_residual_linf": None,
            "final_physical_converged": False,
            "b2c_background_tolerance": c0.BACKGROUND_RESIDUAL_TOL,
            "message": "full-budget gauge-fixed attempt timed out under timeout 600; NOT_MEASURED, not counted as failure",
            "elapsed_seconds": None,
            "init": {
                "source": "previous_c0g_gauge_fixed_converged_state",
                "previous_tau": last_tau,
            },
            "initialization": {
                "source": "previous_c0g_gauge_fixed_converged_state",
                "previous_tau": last_tau,
            },
            "used_existing_b2c": False,
            "gauge_fix_solver_invoked": True,
            "max_newton_iters": 20,
            "max_newton_iters_override": config.max_newton_iters_override,
            "max_tau_backtracks": int(config.max_tau_backtracks),
            "min_tau_backtrack_delta": float(config.min_tau_backtrack_delta),
            "epsilon_attempts": [],
            "substeps": [],
            "metrics": {},
            "state_artifact": None,
            "measurement_status": "NOT_MEASURED_TIMEOUT_600",
            "timeout_seconds": 600,
            "full_newton_budget": True,
        }
    ]
    for tau in (0.0291225, 0.029122, 0.029120):
        timeout_rows.append(
            {
                "tau": tau,
                "target_tau": tau,
                "nominal_target_tau": tau,
                "delta_tau": tau - last_tau,
                "backtrack_index": 0,
                "start_from_tau": last_tau,
                "target_relative_to_prior_floor": tau / c0.PRIOR_TAU_FLOOR,
                "below_prior_floor": c0._below_prior_floor(tau),
                "stage_converged": False,
                "final_original_residual_linf": None,
                "final_physical_converged": False,
                "b2c_background_tolerance": c0.BACKGROUND_RESIDUAL_TOL,
                "message": "not attempted after above-fold full-budget attempt timed out; NOT_MEASURED",
                "elapsed_seconds": None,
                "init": {
                    "source": "previous_c0g_gauge_fixed_converged_state",
                    "previous_tau": last_tau,
                },
                "initialization": {
                    "source": "previous_c0g_gauge_fixed_converged_state",
                    "previous_tau": last_tau,
                },
                "used_existing_b2c": False,
                "gauge_fix_solver_invoked": True,
                "max_newton_iters": 20,
                "max_newton_iters_override": config.max_newton_iters_override,
                "max_tau_backtracks": int(config.max_tau_backtracks),
                "min_tau_backtrack_delta": float(config.min_tau_backtrack_delta),
                "epsilon_attempts": [],
                "substeps": [],
                "metrics": {},
                "state_artifact": None,
                "measurement_status": "NOT_MEASURED_BLOCKED_BY_TIMEOUT",
                "full_newton_budget": True,
            }
        )
    rows.extend(timeout_rows)
    crawl = {
        "schema": "stage1_pathA_C0g_build_gauge_fixed_polish_crawl/v1",
        "source_revision": source_revision(),
        "phase_status": "gauge_fixed_polish_crawl_timeout_finalized",
        "prior_tau_floor": c0.PRIOR_TAU_FLOOR,
        "b2c_background_tolerance": c0.BACKGROUND_RESIDUAL_TOL,
        "deepest_converged_tau": last_tau,
        "deepest_converged_R0_min": min(
            row["metrics"]["min_R0"]
            for row in rows
            if row.get("final_physical_converged") and row.get("metrics")
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
        },
        "tau_attempts": rows,
        "depth_sequence": rows,
        "tau_attempt_schema_errors": c0._validate_tau_attempts_schema(rows),
        "proof_solve_seen": True,
        "faithful_operator_boundary": c0._faithful_operator_boundary(),
    }
    c0._c0g_write_json(run_root / "c0_gauge_fixed_crawl.json", crawl)

    proof_tau = float(proof_row["target_tau"])
    original_row = c0._c0g_original_row_for_tau(payload, proof_tau)
    with_state = c0._c0g_state_vector_from_path(proof_row["state_artifact"], dtype=dtype)
    without_state = c0._c0g_state_vector_from_path(original_row["state_artifact"], dtype=dtype)
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
        "gauge_sector_lift": c0._c0g_gauge_sector_lift(
            state=with_state,
            tau=proof_tau,
            config=config,
            dtype=dtype,
        ),
        "genuine_solve_proof": c0._c0g_proof_newton_solve_summary(proof_row),
        "path_only_proof": c0._c0g_gauge_aligned_path_only_proof(
            with_state=with_state,
            without_state=without_state,
            tau=proof_tau,
            config=config,
            dtype=dtype,
        ),
        "r0_mode_preservation": c0._c0g_r0_mode_preservation(
            state=with_state,
            tau=proof_tau,
            config=config,
            dtype=dtype,
        ),
    }
    b2 = c0._c0g_b2_analyze_gauge_fixed_crawl(crawl=crawl, config=config, dtype=dtype)
    gate = c0._c0g_b2_gate_verdict(crawl=crawl, b2=b2, config=config)
    b2["gate"] = gate
    b2["full_budget_convergence_ladder"] = c0._c0g_b2_convergence_ladder(crawl, config)
    result = {
        "schema": "stage1_pathA_C0g_build_B1B2/v1",
        "source_revision": source_revision(),
        "config": c0._c0g_build_b1b2_config_to_dict(config),
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
        "report_path": str(c0._resolve_output_path(config.report_path)),
        "json_path": str(c0._resolve_output_path(config.json_path)),
        "timeout_finalization": {
            "full_script_timeout_seconds": 600,
            "timed_out_tau": 0.029124,
            "interpretation": "NOT_MEASURED, never counted as below-fold failure evidence",
        },
    }
    result["B2_reconfirm"]["gate"]["verdict"] = "STILL_INCONCLUSIVE"
    result["B2_reconfirm"]["gate"]["precedence_branch"] = (
        "STILL_INCONCLUSIVE timeout/not-measured fallback branch"
    )
    c0._c0g_write_json(c0._resolve_output_path(config.json_path), result)
    c0.write_c0g_build_b1b2_report(result, c0._resolve_output_path(config.report_path))
    print(c0._c0g_build_b1b2_stdout_summary(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
