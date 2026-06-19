"""Path-A B2c calibration and soft-wall floor reporting.

This module is intentionally a thin integration layer over B1/B2a/B2b:
closed backgrounds, BdG moments, and Maxwell transfer are re-used from the
validated engines. The new B2c surface is the direct-coefficient assembly,
stable-side root logic, confirmation gate, error propagation, and reporting.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.optimize import brentq, least_squares, root
import sympy as sp
import torch

from . import patha_b2a_bdg as b2a
from . import patha_b2b_maxwell as b2b
from . import patha_extraction as pe
from .coupled_branch import (
    ClosedCoupledFields,
    pack_closed_coupled_fields,
    patha_closed_wall_terms,
    unpack_closed_coupled_fields,
)


DEFAULT_RUN_ROOT = Path("software/stage1_solver/runs/patha_b2c_calibration")
DEFAULT_REPORT_PATH = Path("software/stage1_solver/reports/patha_b2c_calibration_report.md")
DEFAULT_B2A_VALIDATED = (
    b2a.DEFAULT_RUN_ROOT / "bundles" / "patha_b2a_validated_bdg_bundle_tau_1.json"
)
DEFAULT_B2B_VALIDATED = (
    b2b.DEFAULT_RUN_ROOT / "bundles" / "patha_b2b_validated_maxwell_transfer_tau_1.json"
)
TARGET_P0 = 54.0 / 5.0
TOL_RNORM = 1.0e-6
BACKGROUND_RESIDUAL_TOL = 1.0e-6
DIRECT_KEYS = ("K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4")
ROOT_KEYS = ("D0", "P0", "R_norm", "R_pole", "P2", "P4")
HELD_OUT_KEYS = ("R_pole", "P2", "P4")
OBSERVABLE_ERROR_KEYS = ("D0", "R_norm", "R_pole", "P2", "P4")
ASSEMBLY_ABS_TOL = 1.0e-11
ASSEMBLY_REL_TOL = 1.0e-10
MODEL_COEFF_KEYS = ("B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4", "M")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _stable_hash(obj: Any) -> str:
    text = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )


def _format_tau(tau: float) -> str:
    return b2a._format_tau(float(tau))


def _parse_grid(text: str) -> tuple[int, int]:
    parts = [int(part) for part in text.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("grid must be NR,NW")
    return parts[0], parts[1]


def _grid_from_background(background: Mapping[str, Any], *, dtype: torch.dtype):
    grid_level = (int(background["grid"]["nr"]), int(background["grid"]["nw"]))
    tau = float(background["constants"]["tau"])
    branch = b2a.frozen_patha_b2a_branch(grid=grid_level, tau=tau)
    return b2a._create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")


def _closed_state_from_background(
    background: Mapping[str, Any],
    *,
    dtype: torch.dtype,
) -> tuple[torch.Tensor, Any]:
    grid = _grid_from_background(background, dtype=dtype)

    def tensor(key: str) -> torch.Tensor:
        return torch.as_tensor(background["fields"][key], dtype=dtype, device="cpu")

    fields = ClosedCoupledFields(
        psi_real=tensor("psi_R0"),
        psi_imag=tensor("psi_I0"),
        a0=tensor("A_00"),
        ar=tensor("A_r0"),
        aw=tensor("A_w0"),
        r0=tensor("R0_w"),
    )
    mu = torch.tensor(float(background["solver"]["chemical_potential"]), dtype=dtype, device="cpu")
    return pack_closed_coupled_fields(fields, mu), grid


def _warm_start_payload(
    warm_start_background: Path | None,
) -> tuple[torch.Tensor | None, Any | None, dict[str, Any] | None]:
    if warm_start_background is None:
        return None, None, None
    dtype = b2a.configure_backend(b2a.BackendConfig())
    background = _load_json(warm_start_background)
    state, grid = _closed_state_from_background(background, dtype=dtype)
    info = {
        "path": str(warm_start_background),
        "tau": float(background["constants"]["tau"]),
        "residual_linf": float(background["residuals"]["closed_stationary_linf"]),
        "content_hash": background["content_hash"],
    }
    return state, grid, info


def _final_eos_only_values(tau: float, grid_level: tuple[int, int]) -> tuple[float, ...]:
    branch = b2a.frozen_patha_b2a_branch(grid=grid_level, tau=float(tau))
    return (float(branch.continuation_K_values[-1]),)


def _apply_wall_radius_predictor(
    *,
    state: torch.Tensor,
    grid: Any,
    tau: float,
    grid_level: tuple[int, int],
) -> tuple[torch.Tensor, dict[str, Any]]:
    branch = b2a.frozen_patha_b2a_branch(grid=grid_level, tau=float(tau))
    provider = b2a.resolve_s_sigma(b2a.frozen_s_sigma_spec(float(tau)))
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    assert mu is not None
    dtype = state.dtype
    x0 = fields.r0.detach().cpu().numpy().astype(np.float64)

    def residual(r0_values: np.ndarray) -> np.ndarray:
        r0 = torch.as_tensor(r0_values, dtype=dtype, device="cpu")
        trial_fields = ClosedCoupledFields(
            psi_real=fields.psi_real,
            psi_imag=fields.psi_imag,
            a0=fields.a0,
            ar=fields.ar,
            aw=fields.aw,
            r0=r0,
        )
        return (
            patha_closed_wall_terms(trial_fields, grid, branch, s_sigma=provider)
            .residual.detach()
            .cpu()
            .numpy()
            .astype(np.float64)
        )

    initial = residual(x0)
    initial_linf = float(np.max(np.abs(initial)))
    solved = root(residual, x0, method="hybr", options={"xtol": 1.0e-10, "maxfev": 2000})
    method = "scipy.root(hybr)"
    candidate = np.asarray(solved.x, dtype=np.float64)
    candidate_residual = residual(candidate)
    candidate_linf = float(np.max(np.abs(candidate_residual)))
    success = bool(solved.success and candidate_linf < initial_linf and np.all(candidate > 0.0))
    message = str(solved.message)
    if not success:
        bounded = least_squares(
            residual,
            x0,
            bounds=(0.05, 2.0),
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=2000,
        )
        bounded_candidate = np.asarray(bounded.x, dtype=np.float64)
        bounded_linf = float(np.max(np.abs(residual(bounded_candidate))))
        if bounded.success and bounded_linf < candidate_linf and np.all(bounded_candidate > 0.0):
            method = "scipy.least_squares(bounds)"
            candidate = bounded_candidate
            candidate_linf = bounded_linf
            success = candidate_linf < initial_linf
            message = str(bounded.message)
    if not success:
        return state, {
            "applied": False,
            "method": method,
            "initial_wall_linf": initial_linf,
            "predicted_wall_linf": candidate_linf,
            "message": message,
        }
    predicted_fields = ClosedCoupledFields(
        psi_real=fields.psi_real,
        psi_imag=fields.psi_imag,
        a0=fields.a0,
        ar=fields.ar,
        aw=fields.aw,
        r0=torch.as_tensor(candidate, dtype=dtype, device="cpu"),
    )
    predicted_state = pack_closed_coupled_fields(predicted_fields, mu)
    return predicted_state, {
        "applied": True,
        "method": method,
        "initial_wall_linf": initial_linf,
        "predicted_wall_linf": candidate_linf,
        "r0_min": float(np.min(candidate)),
        "r0_max": float(np.max(candidate)),
        "message": message,
    }


def _max_abs_rel(left: Sequence[float], right: Sequence[float]) -> dict[str, float]:
    l_arr = np.asarray(left, dtype=np.float64)
    r_arr = np.asarray(right, dtype=np.float64)
    diff = np.abs(l_arr - r_arr)
    denom = np.maximum(np.abs(r_arr), 1.0e-300)
    return {
        "max_abs": float(np.max(diff)) if diff.size else 0.0,
        "max_rel": float(np.max(diff / denom)) if diff.size else 0.0,
    }


def _strict_abs_rel(diff: Mapping[str, float]) -> bool:
    return bool(float(diff["max_abs"]) <= ASSEMBLY_ABS_TOL and float(diff["max_rel"]) <= ASSEMBLY_REL_TOL)


def direct_coefficients_from_packets(
    *,
    bdg_packet: Mapping[str, Any],
    maxwell_packet: Mapping[str, Any],
) -> dict[str, float]:
    wall = bdg_packet["wall"]
    bdg = bdg_packet["bdg_moments"]
    coeff = {
        "K": float(wall["K"]),
        "M": float(wall["M"]),
        "B0": float(bdg["B0"]),
        "B2": float(bdg["B2"]),
        "B4": float(bdg["B4"]),
    }
    coeff.update({key: float(maxwell_packet["coefficients"][key]) for key in b2b.COEFF_KEYS})
    missing = sorted(set(DIRECT_KEYS) - set(coeff))
    if missing:
        raise ValueError(f"missing direct coefficients: {missing}")
    return {key: float(coeff[key]) for key in DIRECT_KEYS}


def assemble_primary(coefficients: Mapping[str, float]) -> dict[str, float]:
    coeff = {key: float(coefficients[key]) for key in DIRECT_KEYS}
    packet = pe.lane_extract({"direct_coefficients": coeff})
    residuals = pe.observable_residuals(
        P0=float(packet["P0"]),
        D0=float(packet["D0"]),
        B4=coeff["B4"],
        Z4=coeff["Z4"],
        M=coeff["M"],
        B2=coeff["B2"],
        Z2=coeff["Z2"],
    )
    return {
        **{key: float(packet[key]) for key in packet if isinstance(packet[key], (int, float, np.floating))},
        "R_norm": float(residuals["R_norm"]),
        "R_pole": float(residuals["R_pole"]),
    }


def assemble_secondary(coefficients: Mapping[str, float]) -> dict[str, float]:
    """Independent assembly by symbolic amplitude-series expansion in x=omega^2.

    This deliberately does not call B1 ``lane_extract`` or copy its P2/P4
    formulas. It builds the squared response amplitude from the independently
    derived ``u2/u4`` series, expands through x^2, and reads the coefficients.
    """

    vals = {key: sp.Float(float(coefficients[key]), 50) for key in DIRECT_KEYS}
    x = sp.symbols("x")
    d0 = vals["K"] - vals["B0"] - vals["Z0"]
    d2 = -(vals["M"] + vals["B2"] + vals["Z2"])
    d4 = -(vals["B4"] + vals["Z4"])
    u2 = -d2 / d0
    u4 = (d2**2 - d0 * d4) / d0**2
    numer = vals["N0"] + vals["N2"] * x + vals["N4"] * x**2
    amplitude = 1 + u2 * x + u4 * x**2
    series = sp.series((numer / d0) * amplitude**2, x, 0, 3).removeO().expand()
    p0 = series.coeff(x, 0)
    p2 = series.coeff(x, 1)
    p4 = series.coeff(x, 2)
    return {
        "D0": float(d0),
        "D2": float(d2),
        "D4": float(d4),
        "P0": float(p0),
        "P2": float(p2),
        "P4": float(p4),
        "R_norm": float(p0 - sp.Rational(54, 5)),
        "R_pole": float(d0 * (vals["B4"] + vals["Z4"]) - sp.Integer(3) * (vals["M"] + vals["B2"] + vals["Z2"]) ** 2),
    }


def compare_assembly(primary: Mapping[str, float], secondary: Mapping[str, float]) -> dict[str, Any]:
    per_key = {
        key: _max_abs_rel([float(primary[key])], [float(secondary[key])])
        for key in ROOT_KEYS
    }
    all_diff = _max_abs_rel(
        [float(primary[key]) for key in ROOT_KEYS],
        [float(secondary[key]) for key in ROOT_KEYS],
    )
    return {
        "passed": bool(_strict_abs_rel(all_diff) and all(_strict_abs_rel(row) for row in per_key.values())),
        "diffs": all_diff,
        "per_key": per_key,
        "criteria": f"abs<={ASSEMBLY_ABS_TOL:.1e} AND rel<={ASSEMBLY_REL_TOL:.1e} on {ROOT_KEYS}",
    }


def require_stable_side(observables: Mapping[str, float]) -> None:
    if not math.isfinite(float(observables["D0"])) or float(observables["D0"]) <= 0.0:
        raise ValueError(f"unstable Schur side rejected: D0={observables['D0']!r}")


def assemble_dual(coefficients: Mapping[str, float]) -> dict[str, Any]:
    primary = assemble_primary(coefficients)
    secondary = assemble_secondary(coefficients)
    dual = compare_assembly(primary, secondary)
    stable = bool(float(primary["D0"]) > 0.0)
    return {
        "primary": primary,
        "secondary": secondary,
        "dual_engine": dual,
        "stable_side": {"D0_positive": stable, "D0": float(primary["D0"])},
    }


def frozen_background_root(coefficients: Mapping[str, float]) -> float:
    kappa_hat = pe.hooke_kappa_hat()
    return (
        float(coefficients["B0"]) + float(coefficients["Z0"]) + float(coefficients["N0"]) / TARGET_P0
    ) / kappa_hat


def schur_critical_tau(coefficients: Mapping[str, float]) -> float:
    return (float(coefficients["B0"]) + float(coefficients["Z0"])) / pe.hooke_kappa_hat()


def naturalness(coefficients: Mapping[str, float], tau: float) -> dict[str, float]:
    k = float(coefficients["K"])
    denom_sum = float(coefficients["B0"]) + float(coefficients["Z0"])
    d0 = k - denom_sum
    cancel_fraction = abs(d0) / max(abs(k), abs(denom_sum), 1.0e-300)
    return {
        "tau": float(tau),
        "abs_ln_tau": abs(math.log(float(tau))),
        "K_over_B0_plus_Z0": k / denom_sum if denom_sum != 0.0 else math.inf,
        "D0": d0,
        "cancellation_fraction": cancel_fraction,
        "leading_cancelled_digits": max(0.0, -math.log10(max(cancel_fraction, 1.0e-300))),
    }


def coefficient_relative_errors(
    *,
    b2a_validated: Mapping[str, Any],
    b2b_validated: Mapping[str, Any],
) -> dict[str, float]:
    b_budget = b2a_validated["error_budget"]["B_moments"]
    b_rel: dict[str, float] = {}
    for key in ("B0", "B2", "B4"):
        pieces = [
            float(b_budget.get("modal_truncation_rel", {}).get(key, 0.0)),
            float(b_budget.get("spatial_confirmation_rel", {}).get(key, 0.0)),
            float(b_budget.get("spatial_ladder_rel", {}).get(key, 0.0)),
        ]
        b_rel[key] = math.sqrt(sum(piece * piece for piece in pieces))

    maxwell_rel = {
        key: float(value)
        for key, value in b2b_validated["error_budget"]["Maxwell_ZN_rel"].items()
    }
    rel = {
        "K": 0.0,
        "M": 0.0,
        **b_rel,
        **maxwell_rel,
    }
    return {key: float(rel.get(key, 0.0)) for key in DIRECT_KEYS}


def coefficient_absolute_errors(
    coefficients: Mapping[str, float],
    rel_errors: Mapping[str, float],
) -> dict[str, float]:
    return {key: abs(float(coefficients[key])) * float(rel_errors.get(key, 0.0)) for key in DIRECT_KEYS}


def propagate_observable_errors(
    coefficients: Mapping[str, float],
    rel_errors: Mapping[str, float],
) -> dict[str, Any]:
    base = assemble_primary(coefficients)
    variances = {key: 0.0 for key in OBSERVABLE_ERROR_KEYS}
    contributions: dict[str, dict[str, float]] = {}
    for coeff_key in DIRECT_KEYS:
        delta = abs(float(coefficients[coeff_key])) * float(rel_errors.get(coeff_key, 0.0))
        if delta == 0.0:
            continue
        plus = dict(coefficients)
        minus = dict(coefficients)
        plus[coeff_key] = float(plus[coeff_key]) + delta
        minus[coeff_key] = float(minus[coeff_key]) - delta
        plus_obs = assemble_primary(plus)
        minus_obs = assemble_primary(minus)
        contributions[coeff_key] = {}
        for obs_key in OBSERVABLE_ERROR_KEYS:
            spread = 0.5 * abs(float(plus_obs[obs_key]) - float(minus_obs[obs_key]))
            variances[obs_key] += spread * spread
            contributions[coeff_key][obs_key] = spread
    errors = {key: math.sqrt(value) for key, value in variances.items()}
    root_tau_variance = 0.0
    kappa_hat = pe.hooke_kappa_hat()
    for coeff_key, derivative in (
        ("B0", 1.0 / kappa_hat),
        ("Z0", 1.0 / kappa_hat),
        ("N0", 1.0 / (TARGET_P0 * kappa_hat)),
    ):
        delta = abs(float(coefficients[coeff_key])) * float(rel_errors.get(coeff_key, 0.0))
        root_tau_variance += (derivative * delta) ** 2
    return {
        "base_observables": {key: float(base[key]) for key in OBSERVABLE_ERROR_KEYS},
        "absolute_errors": errors,
        "coefficient_contributions": contributions,
        "tau_root_absolute_error_local": math.sqrt(root_tau_variance),
        "method": "central finite-difference coefficient perturbation using recorded B2a/B2b relative budgets",
    }


def _background_success(background: Mapping[str, Any]) -> bool:
    return bool(
        background["residuals"]["self_consistent"]
        and float(background["residuals"]["closed_stationary_linf"]) <= BACKGROUND_RESIDUAL_TOL
    )


def _evaluation_path(run_root: Path, tau: float, kind: str) -> Path:
    return run_root / "evaluations" / f"patha_b2c_{kind}_tau_{_format_tau(tau)}.json"


def make_evaluation_bundle(
    *,
    tau: float,
    kind: str,
    background: Mapping[str, Any],
    bdg_packet: Mapping[str, Any],
    maxwell_packet: Mapping[str, Any],
    paths: Mapping[str, str],
    source_note: str,
) -> dict[str, Any]:
    coeff = direct_coefficients_from_packets(bdg_packet=bdg_packet, maxwell_packet=maxwell_packet)
    assembled = assemble_dual(coeff)
    primary = assembled["primary"]
    bundle: dict[str, Any] = {
        "schema": "stage1_patha_b2c_tau_evaluation/v1",
        "scope": "one R_norm(tau) evaluation from direct B1 coefficients; target-blind for held-out observables",
        "tau": float(tau),
        "kind": kind,
        "source_note": source_note,
        "paths": dict(paths),
        "r1_trace": {
            "closed_background": {
                "content_hash": background["content_hash"],
                "converged": bool(background["residuals"]["self_consistent"]),
                "residual_linf": float(background["residuals"]["closed_stationary_linf"]),
                "solver_message": background["solver"]["message"],
                "solver_initialization": background["solver"].get("initialization", {}),
            },
            "bdg": {
                "content_hash": bdg_packet["content_hash"],
                "engine": bdg_packet.get("engine", "validated_bundle"),
                "modes_exported": int(
                    bdg_packet.get("diagnostics", {}).get(
                        "modes_exported",
                        len(bdg_packet.get("bdg_modes", [])),
                    )
                ),
            },
            "maxwell": {
                "content_hash": maxwell_packet["content_hash"],
                "engine": maxwell_packet.get("engine", "validated_bundle"),
                "grid": maxwell_packet.get("grid"),
            },
        },
        "direct_coefficients": coeff,
        "observables": {
            "D0": float(primary["D0"]),
            "P0": float(primary["P0"]),
            "R_norm": float(primary["R_norm"]),
            "R_pole": float(primary["R_pole"]),
            "P2": float(primary["P2"]),
            "P4": float(primary["P4"]),
        },
        "assembly": assembled,
        "stable_side": assembled["stable_side"],
        "target_blind": {
            "calibration_target": "R_norm=0 only",
            "held_out_reported_only": list(HELD_OUT_KEYS),
            "no_held_out_targets_used": True,
        },
    }
    bundle["content_hash"] = _stable_hash(bundle)
    return bundle


def seed_tau1_from_validated(
    *,
    run_root: Path = DEFAULT_RUN_ROOT,
    b2a_validated_path: Path = DEFAULT_B2A_VALIDATED,
    b2b_validated_path: Path = DEFAULT_B2B_VALIDATED,
) -> tuple[Path, dict[str, Any]]:
    bdg = _load_json(b2a_validated_path)
    maxwell = _load_json(b2b_validated_path)
    background = _load_json(Path(bdg["background"]["path"]))
    bundle = make_evaluation_bundle(
        tau=float(bdg["tau"]),
        kind="seed_validated",
        background=background,
        bdg_packet=bdg,
        maxwell_packet=maxwell,
        paths={
            "background": str(bdg["background"]["path"]),
            "bdg": str(b2a_validated_path),
            "maxwell": str(b2b_validated_path),
        },
        source_note="Reused validated real tau=1 B2a/B2b production bundles; no frozen-background shortcut.",
    )
    out = _evaluation_path(run_root, float(bundle["tau"]), "seed_validated")
    _write_json(out, bundle)
    return out, bundle


def probe_background(
    *,
    tau: float,
    run_root: Path = DEFAULT_RUN_ROOT,
    background_grid: tuple[int, int] = b2a.DEFAULT_BACKGROUND_GRID,
    warm_start_background: Path | None = None,
    warm_start_final_only: bool = False,
    warm_start_wall_predictor: bool = False,
    newton_jvp_epsilon: float | None = None,
) -> tuple[Path, dict[str, Any]]:
    started = time.perf_counter()
    initial_state, initial_grid, warm_info = _warm_start_payload(warm_start_background)
    wall_predictor = None
    if warm_start_wall_predictor and initial_state is not None and initial_grid is not None:
        initial_state, wall_predictor = _apply_wall_radius_predictor(
            state=initial_state,
            grid=initial_grid,
            tau=float(tau),
            grid_level=background_grid,
        )
    continuation_values = (
        _final_eos_only_values(float(tau), background_grid) if warm_start_final_only else None
    )
    background = b2a.make_background_bundle(
        tau=float(tau),
        grid_level=background_grid,
        initial_state=initial_state,
        initial_grid=initial_grid,
        continuation_K_values=continuation_values,
        initialization_note=None
        if warm_info is None
        else f"B2c tau-homotopy warm start from tau={warm_info['tau']:.12e}",
        newton_overrides={}
        if newton_jvp_epsilon is None
        else {"finite_difference_jvp_epsilon": float(newton_jvp_epsilon)},
    )
    elapsed = time.perf_counter() - started
    record = {
        "schema": "stage1_patha_b2c_background_probe/v1",
        "tau": float(tau),
        "background_grid": list(background_grid),
        "elapsed_seconds": float(elapsed),
        "converged": _background_success(background),
        "solver_converged_flag": bool(background["solver"]["converged"]),
        "solver_message": background["solver"]["message"],
        "residual_linf": float(background["residuals"]["closed_stationary_linf"]),
        "self_consistent": bool(background["residuals"]["self_consistent"]),
        "r0_min": float(background["convergence_evidence"]["r0_min"]),
        "r0_max": float(background["convergence_evidence"]["r0_max"]),
        "background_content_hash": background["content_hash"],
        "warm_start": warm_info,
        "warm_start_final_eos_only": bool(warm_start_final_only),
        "warm_start_wall_predictor": wall_predictor,
        "newton_jvp_epsilon": None if newton_jvp_epsilon is None else float(newton_jvp_epsilon),
    }
    record["content_hash"] = _stable_hash(record)
    bg_path = run_root / "backgrounds" / f"patha_b2c_background_tau_{_format_tau(tau)}.json"
    probe_path = run_root / "background_probes" / f"patha_b2c_background_probe_tau_{_format_tau(tau)}.json"
    _write_json(bg_path, background)
    _write_json(probe_path, record)
    return probe_path, record


def record_timeout_probe(
    *,
    tau: float,
    seconds: float,
    run_root: Path = DEFAULT_RUN_ROOT,
) -> tuple[Path, dict[str, Any]]:
    record = {
        "schema": "stage1_patha_b2c_background_probe/v1",
        "tau": float(tau),
        "elapsed_seconds": float(seconds),
        "converged": False,
        "solver_converged_flag": False,
        "solver_message": f"external timeout after {seconds:g} seconds",
        "residual_linf": None,
        "self_consistent": False,
        "timeout": True,
    }
    record["content_hash"] = _stable_hash(record)
    path = run_root / "background_probes" / f"patha_b2c_background_probe_tau_{_format_tau(tau)}.json"
    _write_json(path, record)
    return path, record


def evaluate_tau(
    *,
    tau: float,
    run_root: Path = DEFAULT_RUN_ROOT,
    kind: str = "resolved",
    background_grid: tuple[int, int] = b2a.DEFAULT_BACKGROUND_GRID,
    bdg_grid: tuple[int, int] = (10, 10),
    maxwell_grid: tuple[int, int] = b2b.DEFAULT_FINAL_GRID,
    profile_points: int = b2a.DEFAULT_PROFILE_POINTS,
    modes: int = b2a.DEFAULT_BDG_MODES,
    window: float = b2b.DEFAULT_FINAL_WINDOW,
    radial_scale: float = b2b.DEFAULT_FINAL_TRUNCATION,
    warm_start_background: Path | None = None,
    warm_start_final_only: bool = False,
    warm_start_wall_predictor: bool = False,
    background_bundle: Path | None = None,
    newton_jvp_epsilon: float | None = None,
) -> tuple[Path, dict[str, Any]]:
    tau = float(tau)
    background_path = run_root / "backgrounds" / f"patha_b2c_background_tau_{_format_tau(tau)}.json"
    bdg_path = run_root / "bdg" / f"patha_b2c_bdg_tau_{_format_tau(tau)}_nr_{bdg_grid[0]}_nw_{bdg_grid[1]}.json"
    maxwell_path = run_root / "maxwell" / (
        f"patha_b2c_maxwell_tau_{_format_tau(tau)}_nr_{maxwell_grid[0]}_nw_{maxwell_grid[1]}"
        f"_w_{window:.6g}_rs_{radial_scale:.6g}.json"
    )

    started = time.perf_counter()
    warm_info = None
    wall_predictor = None
    if background_bundle is not None:
        background = _load_json(background_bundle)
        if abs(float(background["constants"]["tau"]) - tau) > max(1.0e-15, 1.0e-12 * abs(tau)):
            raise ValueError(f"background bundle tau does not match requested tau={tau}")
        background_path = background_bundle
    else:
        initial_state, initial_grid, warm_info = _warm_start_payload(warm_start_background)
        if warm_start_wall_predictor and initial_state is not None and initial_grid is not None:
            initial_state, wall_predictor = _apply_wall_radius_predictor(
                state=initial_state,
                grid=initial_grid,
                tau=tau,
                grid_level=background_grid,
            )
        continuation_values = _final_eos_only_values(tau, background_grid) if warm_start_final_only else None
        background = b2a.make_background_bundle(
            tau=tau,
            grid_level=background_grid,
            initial_state=initial_state,
            initial_grid=initial_grid,
            continuation_K_values=continuation_values,
            initialization_note=None
            if warm_info is None
            else f"B2c tau-homotopy warm start from tau={warm_info['tau']:.12e}",
            newton_overrides={}
            if newton_jvp_epsilon is None
            else {"finite_difference_jvp_epsilon": float(newton_jvp_epsilon)},
        )
        _write_json(background_path, background)
    if not _background_success(background):
        failure = {
            "schema": "stage1_patha_b2c_tau_evaluation_failure/v1",
            "tau": tau,
            "kind": kind,
            "stage": "closed_background",
            "elapsed_seconds": time.perf_counter() - started,
            "background_path": str(background_path),
            "residual_linf": float(background["residuals"]["closed_stationary_linf"]),
            "self_consistent": bool(background["residuals"]["self_consistent"]),
            "solver_converged_flag": bool(background["solver"]["converged"]),
            "solver_message": background["solver"]["message"],
            "warm_start": warm_info,
            "warm_start_final_eos_only": bool(warm_start_final_only),
            "warm_start_wall_predictor": wall_predictor,
            "newton_jvp_epsilon": None if newton_jvp_epsilon is None else float(newton_jvp_epsilon),
        }
        failure["content_hash"] = _stable_hash(failure)
        fail_path = run_root / "failures" / f"patha_b2c_failure_tau_{_format_tau(tau)}.json"
        _write_json(fail_path, failure)
        raise RuntimeError(f"closed background did not converge at tau={tau}: {failure}")

    wall_input, wall_path = b2a.make_wall_input(
        background,
        profile_points_count=profile_points,
        out_dir=run_root / "wall_inputs",
    )
    bdg_packet = b2a.solve_bdg_python(
        background,
        wall_input,
        nr=int(bdg_grid[0]),
        nw=int(bdg_grid[1]),
        modes_to_export=int(modes),
    )
    _write_json(bdg_path, bdg_packet)
    maxwell_packet = b2b.transfer_for_grid(
        background,
        bdg_packet,
        nr=int(maxwell_grid[0]),
        nw=int(maxwell_grid[1]),
        window=float(window),
        radial_scale=float(radial_scale),
        mode_count=int(modes),
        engine="python_primary_second_order",
        discretization="primary_second_order",
    )
    _write_json(maxwell_path, maxwell_packet)
    bundle = make_evaluation_bundle(
        tau=tau,
        kind=kind,
        background=background,
        bdg_packet=bdg_packet,
        maxwell_packet=maxwell_packet,
        paths={
            "background": str(background_path),
            "wall_input": "" if wall_path is None else str(wall_path),
            "bdg": str(bdg_path),
            "maxwell": str(maxwell_path),
        },
        source_note="Fresh B2c full per-tau re-solve: closed background, B2a BdG, and B2b Maxwell transfer.",
    )
    if background_bundle is not None:
        bundle["source_note"] = (
            "B2c per-tau re-solve from a precomputed residual-gated closed background, followed by fresh "
            "B2a BdG and B2b Maxwell transfer."
        )
        bundle["content_hash"] = _stable_hash(bundle)
    if wall_predictor is not None:
        bundle["r1_trace"]["closed_background"]["warm_start_wall_predictor"] = wall_predictor
        bundle["content_hash"] = _stable_hash(bundle)
    out = _evaluation_path(run_root, tau, kind)
    _write_json(out, bundle)
    return out, bundle


def _load_evaluations(run_root: Path) -> list[dict[str, Any]]:
    paths = sorted((run_root / "evaluations").glob("patha_b2c_*_tau_*.json"))
    return [_load_json(path) for path in paths]


def _load_background_probes(run_root: Path) -> list[dict[str, Any]]:
    paths = sorted((run_root / "background_probes").glob("patha_b2c_background_probe_tau_*.json"))
    return [_load_json(path) for path in paths]


def _model_coefficients(samples: Sequence[Mapping[str, Any]], tau: float) -> dict[str, float]:
    if not samples:
        raise ValueError("need at least one sample")
    sorted_samples = sorted(samples, key=lambda row: float(row["tau"]))
    taus = np.asarray([float(row["tau"]) for row in sorted_samples], dtype=np.float64)
    if tau < taus[0] or tau > taus[-1]:
        raise ValueError(f"tau={tau} outside model sample range [{taus[0]}, {taus[-1]}]")
    coeffs = [row["direct_coefficients"] for row in sorted_samples]
    log_taus = np.log(taus)
    log_tau = math.log(float(tau))
    out: dict[str, float] = {"K": pe.hooke_kappa_hat() * float(tau)}
    for key in MODEL_COEFF_KEYS:
        values = np.asarray([float(coeff[key]) for coeff in coeffs], dtype=np.float64)
        out[key] = float(np.interp(log_tau, log_taus, values))
    return {key: float(out[key]) for key in DIRECT_KEYS}


def _model_rnorm(samples: Sequence[Mapping[str, Any]], tau: float, *, secondary: bool = False) -> float:
    coeff = _model_coefficients(samples, tau)
    obs = assemble_secondary(coeff) if secondary else assemble_primary(coeff)
    require_stable_side(obs)
    return float(obs["R_norm"])


def _stable_brackets(samples: Sequence[Mapping[str, Any]]) -> list[tuple[float, float]]:
    rows = sorted(samples, key=lambda row: float(row["tau"]))
    brackets: list[tuple[float, float]] = []
    previous: Mapping[str, Any] | None = None
    for row in rows:
        obs = row["observables"]
        if float(obs["D0"]) <= 0.0 or not math.isfinite(float(obs["R_norm"])):
            previous = None
            continue
        if previous is not None:
            prev_r = float(previous["observables"]["R_norm"])
            curr_r = float(obs["R_norm"])
            if prev_r == 0.0:
                brackets.append((float(previous["tau"]), float(previous["tau"])))
            elif prev_r * curr_r < 0.0:
                brackets.append((float(previous["tau"]), float(row["tau"])))
        previous = row
    return brackets


def root_find_primary(samples: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    brackets = _stable_brackets(samples)
    if not brackets:
        raise ValueError("no stable-side R_norm sign bracket in resolved samples")
    lo, hi = brackets[0]
    if lo == hi:
        root = lo
    else:
        root = float(brentq(lambda t: _model_rnorm(samples, t, secondary=False), lo, hi, xtol=1.0e-18, rtol=1.0e-15))
    coeff = _model_coefficients(samples, root)
    obs = assemble_primary(coeff)
    require_stable_side(obs)
    return {
        "method": "scipy.brentq on log-linear coefficient model",
        "tau": root,
        "bracket": [lo, hi],
        "direct_coefficients": coeff,
        "observables": {key: float(obs[key]) for key in ("D0", "P0", "R_norm", "R_pole", "P2", "P4")},
    }


def root_find_secondary(samples: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    brackets = _stable_brackets(samples)
    if not brackets:
        raise ValueError("no stable-side R_norm sign bracket in resolved samples")
    lo, hi = brackets[0]
    if lo == hi:
        root = lo
    else:
        flo = _model_rnorm(samples, lo, secondary=True)
        fhi = _model_rnorm(samples, hi, secondary=True)
        if flo * fhi > 0.0:
            raise ValueError("secondary root finder lost the sign bracket")
        for _ in range(160):
            mid = 0.5 * (lo + hi)
            fmid = _model_rnorm(samples, mid, secondary=True)
            if abs(fmid) <= 1.0e-12 or abs(hi - lo) <= max(1.0e-18, 1.0e-15 * abs(mid)):
                lo = hi = mid
                break
            if flo * fmid <= 0.0:
                hi = mid
                fhi = fmid
            else:
                lo = mid
                flo = fmid
        root = 0.5 * (lo + hi)
    coeff = _model_coefficients(samples, root)
    obs = assemble_secondary(coeff)
    require_stable_side(obs)
    return {
        "method": "independent bisection plus symbolic-series assembly",
        "tau": float(root),
        "bracket": brackets[0],
        "direct_coefficients": coeff,
        "observables": {key: float(obs[key]) for key in ("D0", "P0", "R_norm", "R_pole", "P2", "P4")},
    }


def confirmation_gate(
    *,
    candidate_tau: float,
    confirmation_bundle: Mapping[str, Any],
    located_coefficients: Mapping[str, float],
    coefficient_abs_errors: Mapping[str, float],
    tol_rnorm: float = TOL_RNORM,
) -> dict[str, Any]:
    obs = confirmation_bundle["observables"]
    stable = float(obs["D0"]) > 0.0
    rnorm_ok = abs(float(obs["R_norm"])) <= float(tol_rnorm)
    coeff_rows = {}
    coeff_ok = True
    for key in DIRECT_KEYS:
        diff = abs(float(confirmation_bundle["direct_coefficients"][key]) - float(located_coefficients[key]))
        allowed = max(float(coefficient_abs_errors.get(key, 0.0)), 1.0e-15)
        coeff_rows[key] = {"abs_diff": diff, "allowed": allowed, "passed": diff <= allowed}
        coeff_ok = coeff_ok and diff <= allowed
    tau_ok = abs(float(confirmation_bundle["tau"]) - float(candidate_tau)) <= max(
        1.0e-15, 1.0e-12 * abs(float(candidate_tau))
    )
    return {
        "passed": bool(stable and rnorm_ok and coeff_ok and tau_ok),
        "stable_side": stable,
        "rnorm_ok": rnorm_ok,
        "coefficients_within_recorded_bars": coeff_ok,
        "tau_matches_candidate": tau_ok,
        "candidate_tau": float(candidate_tau),
        "confirmation_tau": float(confirmation_bundle["tau"]),
        "tol_Rnorm": float(tol_rnorm),
        "R_norm": float(obs["R_norm"]),
        "D0": float(obs["D0"]),
        "coefficient_rows": coeff_rows,
    }


def drift_between(
    low_eval: Mapping[str, Any],
    high_eval: Mapping[str, Any],
) -> dict[str, Any]:
    low = low_eval["direct_coefficients"]
    high = high_eval["direct_coefficients"]
    per = {}
    for key in ("B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4"):
        per[key] = abs(float(low[key]) - float(high[key])) / max(abs(float(high[key])), 1.0e-300)
    return {
        "tau_low": float(low_eval["tau"]),
        "tau_high": float(high_eval["tau"]),
        "per_coefficient_relative_to_high_tau": per,
        "max_rel": max(per.values()) if per else 0.0,
    }


def _local_edge_estimate(coefficients: Mapping[str, float]) -> dict[str, float]:
    b_plus_z = float(coefficients["B0"]) + float(coefficients["Z0"])
    n0 = float(coefficients["N0"])
    kappa_hat = pe.hooke_kappa_hat()
    tau_star = (b_plus_z + n0 / TARGET_P0) / kappa_hat
    d0_star = n0 / TARGET_P0
    cancellation_fraction = d0_star / max(abs(kappa_hat * tau_star), abs(b_plus_z), 1.0e-300)
    return {
        "tau_crit_local_coefficients": b_plus_z / kappa_hat,
        "tau_star_local_coefficients": tau_star,
        "D0_at_tau_star_local_coefficients": d0_star,
        "abs_ln_tau_star_local_coefficients": abs(math.log(tau_star)),
        "K_over_B0_plus_Z0_at_tau_star_local_coefficients": (kappa_hat * tau_star) / b_plus_z,
        "cancellation_fraction_at_tau_star_local_coefficients": cancellation_fraction,
        "leading_cancelled_digits_at_tau_star_local_coefficients": max(
            0.0,
            -math.log10(max(cancellation_fraction, 1.0e-300)),
        ),
        "B0_plus_Z0": b_plus_z,
        "N0": n0,
        "kappa_hat": kappa_hat,
    }


def _coefficient_edge_terms(eval_bundle: Mapping[str, Any]) -> tuple[float, float, float]:
    coeff = eval_bundle["direct_coefficients"]
    return (
        float(eval_bundle["tau"]),
        float(coeff["B0"]) + float(coeff["Z0"]),
        float(coeff["N0"]),
    )


def _positive_power_law_value(prefactor: float, exponent: float, tau: float) -> float:
    if tau <= 0.0:
        return math.nan
    log_value = math.log(prefactor) + exponent * math.log(tau)
    if log_value > 700.0:
        return math.inf
    if log_value < -745.0:
        return 0.0
    return math.exp(log_value)


def _power_law_edge_estimate(
    low_eval: Mapping[str, Any],
    high_eval: Mapping[str, Any],
    *,
    tau_floor: float | None,
) -> dict[str, Any]:
    tau_low, s_low, n_low = _coefficient_edge_terms(low_eval)
    tau_high, s_high, n_high = _coefficient_edge_terms(high_eval)
    if tau_low <= 0.0 or tau_high <= 0.0 or tau_low == tau_high or min(s_low, s_high, n_low, n_high) <= 0.0:
        return {
            "source_taus": [tau_low, tau_high],
            "tau_estimate": None,
            "status": "degenerate_pair",
        }
    s_exponent = math.log(s_high / s_low) / math.log(tau_high / tau_low)
    n_exponent = math.log(n_high / n_low) / math.log(tau_high / tau_low)
    s_prefactor = s_low / (tau_low**s_exponent)
    n_prefactor = n_low / (tau_low**n_exponent)
    kappa_hat = pe.hooke_kappa_hat()

    def residual(tau: float) -> float:
        return (
            kappa_hat * tau
            - _positive_power_law_value(s_prefactor, s_exponent, tau)
            - _positive_power_law_value(n_prefactor, n_exponent, tau) / TARGET_P0
        )

    lo = max(min(tau_low, tau_high) * 1.0e-6, 1.0e-12)
    hi = max(max(tau_low, tau_high) * 10.0, 1.0)
    roots: list[float] = []
    prev_tau = lo
    prev_val = residual(prev_tau)
    for idx in range(1, 800):
        tau = math.exp(math.log(lo) + (math.log(hi) - math.log(lo)) * idx / 799.0)
        val = residual(tau)
        if math.isfinite(prev_val) and math.isfinite(val) and prev_val * val <= 0.0:
            left = prev_tau
            right = tau
            f_left = prev_val
            for _ in range(80):
                mid = 0.5 * (left + right)
                f_mid = residual(mid)
                if not math.isfinite(f_mid):
                    right = mid
                    continue
                if f_left * f_mid <= 0.0:
                    right = mid
                    val = f_mid
                else:
                    left = mid
                    f_left = f_mid
            roots.append(0.5 * (left + right))
            break
        prev_tau = tau
        prev_val = val

    tau_estimate = roots[0] if roots else None
    if tau_estimate is None:
        status = "no_positive_root_found"
    elif tau_floor is not None and tau_estimate >= tau_floor:
        status = "above_floor_degenerate_not_a_valid_bound"
    else:
        status = "below_floor_fit_dependent"
    return {
        "source_taus": [tau_low, tau_high],
        "tau_estimate": tau_estimate,
        "status": status,
        "exponents": {
            "B0_plus_Z0": s_exponent,
            "N0": n_exponent,
        },
        "method": "two-point power-law fit of B0+Z0 and N0; solve K=B0+Z0+N0/P0_target",
    }


def _edge_estimate_fit_spread(
    *,
    evaluations: Sequence[Mapping[str, Any]],
    resolved_low: Mapping[str, Any] | None,
    floor: Mapping[str, Any],
) -> dict[str, Any] | None:
    if resolved_low is None:
        return None
    tau_floor = None if floor.get("tau_floor") is None else float(floor["tau_floor"])
    local_edge = _local_edge_estimate(resolved_low["direct_coefficients"])
    estimates: list[dict[str, Any]] = [
        {
            "fit": "frozen_deepest_coefficients",
            "tau_estimate": local_edge["tau_star_local_coefficients"],
            "status": "below_floor_fit_dependent" if tau_floor is None or local_edge["tau_star_local_coefficients"] < tau_floor else "above_floor",
            "source_taus": [float(resolved_low["tau"])],
            "method": "hold deepest measured B0/Z0/N0 fixed; solve K=B0+Z0+N0/P0_target",
        }
    ]
    sorted_evals = sorted(evaluations, key=lambda row: float(row["tau"]))
    seed = max(sorted_evals, key=lambda row: float(row["tau"]))
    near_003 = min(sorted_evals, key=lambda row: abs(float(row["tau"]) - 0.03))
    estimates.append(
        {
            "fit": "wide_pair_power_law_tau_1_to_0p03",
            **_power_law_edge_estimate(near_003, seed, tau_floor=tau_floor),
        }
    )
    full_resolved = [
        ev
        for ev in sorted_evals
        if ev.get("kind") != "seed_validated" and ev["r1_trace"]["closed_background"]["converged"]
    ]
    if len(full_resolved) >= 2:
        estimates.append(
            {
                "fit": "close_pair_power_law_deepest_two",
                **_power_law_edge_estimate(full_resolved[0], full_resolved[1], tau_floor=tau_floor),
            }
        )
    finite = [
        float(row["tau_estimate"])
        for row in estimates
        if row.get("tau_estimate") is not None and math.isfinite(float(row["tau_estimate"]))
    ]
    span_orders = math.log10(max(finite) / min(finite)) if len(finite) >= 2 and min(finite) > 0.0 else None
    return {
        "purpose": "Record fit dependence of drift-aware edge estimates; this does not pin tau_star.",
        "rigorous_upper_bound": {
            "tau_star_lt": tau_floor,
            "source": floor.get("floor_statement"),
        },
        "target_P0": TARGET_P0,
        "estimates": estimates,
        "finite_estimate_span_orders_including_degenerate": span_orders,
        "tau_star_pinned_by_this_run": False,
    }


def _nearest_failed_below(
    *,
    tau_floor: float | None,
    probes: Sequence[Mapping[str, Any]],
) -> dict[str, Any] | None:
    if tau_floor is None:
        return None
    failed = [
        probe
        for probe in probes
        if not bool(probe.get("converged", False)) and float(probe["tau"]) < float(tau_floor)
    ]
    if not failed:
        return None
    probe = max(failed, key=lambda row: float(row["tau"]))
    return {
        "tau": float(probe["tau"]),
        "residual_linf": probe.get("residual_linf"),
        "message": probe.get("solver_message", ""),
        "background_grid": probe.get("background_grid"),
        "r0_min": probe.get("r0_min"),
        "r0_max": probe.get("r0_max"),
    }


def _maxwell_gridcheck_for_eval(eval_bundle: Mapping[str, Any]) -> dict[str, Any] | None:
    paths = eval_bundle.get("paths", {})
    production_path = paths.get("maxwell")
    if not production_path:
        return None
    production = _load_json(Path(production_path))
    run_root = Path(production_path).parents[1]
    tau_label = _format_tau(float(eval_bundle["tau"]))
    candidates = sorted((run_root / "maxwell").glob(f"patha_b2c_maxwell_tau_{tau_label}_*_gridcheck.json"))
    if not candidates:
        return None
    refined = _load_json(candidates[-1])
    rel = {}
    for key in b2b.COEFF_KEYS:
        prod = float(production["coefficients"][key])
        ref = float(refined["coefficients"][key])
        rel[key] = abs(ref - prod) / max(abs(ref), 1.0e-300)
    return {
        "production_path": production_path,
        "refined_path": str(candidates[-1]),
        "production_grid": production["grid"],
        "refined_grid": refined["grid"],
        "relative_change_vs_refined": rel,
        "max_rel": max(rel.values()) if rel else 0.0,
    }


def _spot_convergence_for_eval(eval_bundle: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if eval_bundle is None:
        return None
    bdg_path = eval_bundle.get("paths", {}).get("bdg")
    modal = None
    if bdg_path:
        bdg = _load_json(Path(bdg_path))
        meta = bdg["diagnostics"]["modal_truncation"]
        modal = {
            "bdg_path": bdg_path,
            "grid": bdg["grid"],
            "exported_mode_count": int(meta["exported_mode_count"]),
            "all_positive_mode_count": int(meta["all_positive_mode_count"]),
            "max_truncation_error_at_export": float(meta["max_truncation_error_at_export"]),
            "tolerance": float(meta["tolerance"]),
            "passed": float(meta["max_truncation_error_at_export"]) <= float(meta["tolerance"]),
        }
    maxwell = _maxwell_gridcheck_for_eval(eval_bundle)
    return {
        "tau": float(eval_bundle["tau"]),
        "modal": modal,
        "maxwell_grid_refinement": maxwell,
    }


def _edge_diagnostic(
    *,
    resolved_low: Mapping[str, Any] | None,
    floor: Mapping[str, Any],
    probes: Sequence[Mapping[str, Any]],
) -> dict[str, Any] | None:
    if resolved_low is None:
        return None
    obs = resolved_low["observables"]
    coeff = resolved_low["direct_coefficients"]
    local_edge = _local_edge_estimate(coeff)
    failed = _nearest_failed_below(tau_floor=floor.get("tau_floor"), probes=probes)
    d0 = float(obs["D0"])
    classification = "numerical"
    evidence = (
        "No D0 or Jacobian-smallest-mode measurement exists at the failing tau because that solve did not "
        "converge. The deepest residual-gated full bundle has D0=O(0.1) and K/(B0+Z0) far from one; the lower "
        "near-floor failures are smooth line-search/continuation stalls with residuals rising monotonically as tau drops. "
        "This is consistent with a numerical (continuation/conditioning) floor, not a diagnosed physical "
        "marginal-stability edge."
    )
    if d0 <= 1.0e-3:
        classification = "physical_or_marginal"
        evidence = "The deepest residual-gated full bundle is already close to D0=0."
    return {
        "deepest_full_evaluation": {
            "tau": float(resolved_low["tau"]),
            "kind": resolved_low.get("kind"),
            "R_norm": float(obs["R_norm"]),
            "D0": d0,
            "P0": float(obs["P0"]),
            "K_over_B0_plus_Z0": naturalness(coeff, float(resolved_low["tau"]))["K_over_B0_plus_Z0"],
            "background_residual_linf": float(
                resolved_low["r1_trace"]["closed_background"]["residual_linf"]
            ),
        },
        "nearest_failed_probe_below_floor": failed,
        "local_edge_estimate_from_deepest_coefficients": local_edge,
        "physical_vs_numerical": classification,
        "evidence": evidence,
        "no_bundle_extrapolated_below_tau": float(resolved_low["tau"]),
    }


def negative_control_injected_drift() -> dict[str, Any]:
    base = {
        "K": pe.hooke_kappa_hat() * 1.0,
        "M": 1.0,
        "B0": 4.0e-5,
        "B2": 1.0e-6,
        "B4": 1.0e-7,
        "Z0": 2.0e-8,
        "Z2": -1.0e-8,
        "Z4": 5.0e-9,
        "N0": 2.0e-5,
        "N2": -5.0e-9,
        "N4": 3.0e-9,
    }
    frozen = frozen_background_root(base)
    low_coeff = dict(base)
    high_coeff = dict(base)
    low_coeff["B0"] = base["B0"] * 1.40
    low_coeff["Z0"] = base["Z0"] * 1.20
    high_coeff["B0"] = base["B0"] * 1.35
    high_coeff["Z0"] = base["Z0"] * 1.15
    low_root = frozen_background_root(low_coeff)
    high_root = frozen_background_root(high_coeff)
    low_tau = low_root * (1.0 - 1.0e-5)
    high_tau = high_root * (1.0 + 2.0e-4)
    for tau, coeff in ((low_tau, low_coeff), (high_tau, high_coeff)):
        coeff["K"] = pe.hooke_kappa_hat() * tau
    samples = []
    for tau, coeff in ((low_tau, low_coeff), (high_tau, high_coeff)):
        obs = assemble_primary(coeff)
        samples.append({"tau": tau, "direct_coefficients": coeff, "observables": obs})
    resolved = root_find_primary(samples)["tau"]
    return {
        "failed_as_expected": abs(resolved - frozen) / max(abs(frozen), 1.0e-300) > 1.0e-2,
        "frozen_tau": frozen,
        "resolved_tau_under_injected_drift": resolved,
        "relative_shift": abs(resolved - frozen) / max(abs(frozen), 1.0e-300),
        "wrong_answer": "frozen-background root finder that ignores tau-dependent {B,Z,N} drift",
    }


def confirmation_negative_control() -> dict[str, Any]:
    coeff = {
        "K": 1.1e-3,
        "M": 1.0,
        "B0": 1.0e-3,
        "B2": 1.0e-6,
        "B4": 1.0e-7,
        "Z0": 1.0e-8,
        "Z2": -1.0e-8,
        "Z4": 5.0e-9,
        "N0": 2.0e-8,
        "N2": -5.0e-9,
        "N4": 3.0e-9,
    }
    obs = assemble_primary(coeff)
    bundle = {"tau": 1.0e-4, "direct_coefficients": coeff, "observables": obs}
    gate = confirmation_gate(
        candidate_tau=1.0e-4,
        confirmation_bundle=bundle,
        located_coefficients=coeff,
        coefficient_abs_errors={key: 1.0 for key in DIRECT_KEYS},
        tol_rnorm=1.0e-6,
    )
    return {
        "failed_as_expected": not gate["passed"],
        "gate_passed": gate["passed"],
        "R_norm": gate["R_norm"],
        "wrong_answer": "mis-located tau accepted without a fresh |R_norm| gate",
    }


def stable_side_negative_control() -> dict[str, Any]:
    coeff = {
        "K": 9.0e-4,
        "M": 1.0,
        "B0": 1.0e-3,
        "B2": 1.0e-6,
        "B4": 1.0e-7,
        "Z0": 1.0e-8,
        "Z2": -1.0e-8,
        "Z4": 5.0e-9,
        "N0": 2.0e-8,
        "N2": -5.0e-9,
        "N4": 3.0e-9,
    }
    obs = assemble_primary(coeff)
    rejected = False
    try:
        require_stable_side(obs)
    except ValueError:
        rejected = True
    return {
        "failed_as_expected": rejected,
        "D0": obs["D0"],
        "wrong_answer": "unstable-side D0<=0 root accepted",
    }


def _tau_floor_from_records(
    evaluations: Sequence[Mapping[str, Any]],
    probes: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    converged_rows: list[dict[str, Any]] = []
    failed_rows: list[dict[str, Any]] = []
    for ev in evaluations:
        source = "validated_seed" if ev.get("kind") == "seed_validated" else "full_evaluation"
        row = {
            "tau": float(ev["tau"]),
            "source": source,
            "converged": bool(ev["r1_trace"]["closed_background"]["converged"]),
            "residual_linf": float(ev["r1_trace"]["closed_background"]["residual_linf"]),
            "message": ev["r1_trace"]["closed_background"]["solver_message"],
        }
        (converged_rows if row["converged"] else failed_rows).append(row)
    for probe in probes:
        row = {
            "tau": float(probe["tau"]),
            "source": "background_probe",
            "converged": bool(probe.get("converged", False)),
            "residual_linf": probe.get("residual_linf"),
            "message": probe.get("solver_message", ""),
        }
        (converged_rows if row["converged"] else failed_rows).append(row)
    tau_floor_row = min(converged_rows, key=lambda row: row["tau"]) if converged_rows else None
    return {
        "tau_floor": None if tau_floor_row is None else float(tau_floor_row["tau"]),
        "tau_floor_row": tau_floor_row,
        "converged_rows": sorted(converged_rows, key=lambda row: row["tau"]),
        "failed_rows": sorted(failed_rows, key=lambda row: row["tau"]),
        "floor_statement": (
            "lowest converged tau among executed production-grid closed-background probes/evaluations; "
            "no background/BdG/Maxwell quantity is extrapolated below it"
        ),
    }


def calibrate_from_records(
    *,
    run_root: Path = DEFAULT_RUN_ROOT,
    report_path: Path = DEFAULT_REPORT_PATH,
    b2a_validated_path: Path = DEFAULT_B2A_VALIDATED,
    b2b_validated_path: Path = DEFAULT_B2B_VALIDATED,
) -> dict[str, Any]:
    evaluations = _load_evaluations(run_root)
    if not evaluations:
        raise RuntimeError("no B2c evaluations found; run seed-tau1 and at least one resolved evaluate")
    probes = _load_background_probes(run_root)
    b2a_validated = _load_json(b2a_validated_path)
    b2b_validated = _load_json(b2b_validated_path)
    rel_errors = coefficient_relative_errors(
        b2a_validated=b2a_validated,
        b2b_validated=b2b_validated,
    )

    seed = min(evaluations, key=lambda ev: abs(float(ev["tau"]) - 1.0))
    seed_coeff = seed["direct_coefficients"]
    tau_frozen = frozen_background_root(seed_coeff)
    tau_crit_frozen = schur_critical_tau(seed_coeff)
    floor = _tau_floor_from_records(evaluations, probes)
    full_resolved = [
        ev for ev in evaluations if ev.get("kind") != "seed_validated" and ev["r1_trace"]["closed_background"]["converged"]
    ]
    resolved_low = min(full_resolved, key=lambda ev: float(ev["tau"])) if full_resolved else None
    high = max(evaluations, key=lambda ev: float(ev["tau"]))
    drift = None if resolved_low is None else drift_between(resolved_low, seed)

    root_primary = root_secondary = root_agreement = confirmation = None
    status = "root_not_bracketed"
    finding = "No stable-side sign bracket was found in the resolved sample set."
    try:
        root_primary = root_find_primary(evaluations)
        root_secondary = root_find_secondary(evaluations)
        tau_diff = abs(float(root_primary["tau"]) - float(root_secondary["tau"]))
        root_agreement = {
            "passed": tau_diff <= max(1.0e-12 * abs(float(root_primary["tau"])), 1.0e-15),
            "tau_abs_diff": tau_diff,
            "primary_tau": float(root_primary["tau"]),
            "secondary_tau": float(root_secondary["tau"]),
        }
        status = "root_model_located_unconfirmed"
        finding = "A stable-side model root was found, but a fresh confirmation bundle is still required."
    except ValueError as exc:
        if resolved_low is not None:
            status = "root_bounded_by_resolved_floor"
            finding = (
                "No stable-side R_norm sign bracket exists in the real converged re-solve set. "
                f"The deepest full bundle is tau={float(resolved_low['tau']):.12e} with "
                f"R_norm={float(resolved_low['observables']['R_norm']):.12e}; any root, if present on the "
                "stable branch, is only bounded above by that residual-gated floor and is not confirmed."
            )
        else:
            finding = str(exc)

    floor_eval = resolved_low if resolved_low is not None else seed
    propagated = propagate_observable_errors(floor_eval["direct_coefficients"], rel_errors)
    natural = naturalness(floor_eval["direct_coefficients"], float(floor_eval["tau"]))
    negative_controls = {
        "frozen_vs_resolved_injected_drift": negative_control_injected_drift(),
        "confirmation_gate_mislocated_tau": confirmation_negative_control(),
        "stable_side_filter": stable_side_negative_control(),
    }
    no_cant_fail = {
        "passed": all(control["failed_as_expected"] for control in negative_controls.values()),
        "negative_controls": negative_controls,
    }
    b2b_dual_engine_max_rel = float(b2b_validated["error_budget"].get("dual_engine_max_rel", math.nan))
    edge = _edge_diagnostic(resolved_low=resolved_low, floor=floor, probes=probes)
    edge_fit_spread = _edge_estimate_fit_spread(
        evaluations=evaluations,
        resolved_low=resolved_low,
        floor=floor,
    )
    spot_convergence = _spot_convergence_for_eval(resolved_low)

    bundle: dict[str, Any] = {
        "schema": "stage1_patha_b2c_calibration/v1",
        "status": status,
        "finding": finding,
        "target": {"R_norm": 0.0, "P0_target": TARGET_P0, "tol_Rnorm": TOL_RNORM},
        "tau_floor": floor,
        "frozen_negative_control": {
            "tau_frozen": tau_frozen,
            "tau_crit_frozen": tau_crit_frozen,
            "source_tau": float(seed["tau"]),
            "source_evaluation_hash": seed["content_hash"],
            "below_tau_floor": None if floor["tau_floor"] is None else tau_frozen < float(floor["tau_floor"]),
        },
        "resolved_floor_evaluation": {
            "tau": float(floor_eval["tau"]),
            "content_hash": floor_eval["content_hash"],
            "observables": floor_eval["observables"],
            "direct_coefficients": floor_eval["direct_coefficients"],
            "naturalness": natural,
        },
        "root_primary": root_primary,
        "root_secondary": root_secondary,
        "root_agreement": root_agreement,
        "confirmation": confirmation,
        "edge_diagnostic": edge,
        "edge_estimate_fit_spread": edge_fit_spread,
        "spot_convergence": spot_convergence,
        "drift": drift,
        "error_bars": {
            "recorded_relative_errors": rel_errors,
            "propagated_at_resolved_floor": propagated,
            "recorded_b2b_dual_engine_max_rel": b2b_dual_engine_max_rel,
            "recorded_budget_sources": {
                "B2a": str(b2a_validated_path),
                "B2b": str(b2b_validated_path),
            },
        },
        "target_blind_firewall": {
            "calibrated_only_on_R_norm": True,
            "held_out_keys": list(HELD_OUT_KEYS),
            "held_out_tuning_used": False,
            "new_degrees_of_freedom": 0,
            "note": (
                "When no confirmed tau* exists above tau_floor, held-out tau* predictions are not emitted; "
                "the floor observables are diagnostic, not a tuned prediction."
            ),
        },
        "dual_engine": {
            "evaluation_assembly_all_passed": all(ev["assembly"]["dual_engine"]["passed"] for ev in evaluations),
            "evaluation_rows": [
                {
                    "tau": float(ev["tau"]),
                    "kind": ev["kind"],
                    "passed": bool(ev["assembly"]["dual_engine"]["passed"]),
                    "max_abs": float(ev["assembly"]["dual_engine"]["diffs"]["max_abs"]),
                    "max_rel": float(ev["assembly"]["dual_engine"]["diffs"]["max_rel"]),
                }
                for ev in sorted(evaluations, key=lambda row: float(row["tau"]))
            ],
            "root_agreement": root_agreement,
        },
        "negative_controls": no_cant_fail,
        "evaluations": [
            {
                "tau": float(ev["tau"]),
                "kind": ev["kind"],
                "content_hash": ev["content_hash"],
                "R_norm": float(ev["observables"]["R_norm"]),
                "D0": float(ev["observables"]["D0"]),
                "stable": bool(ev["stable_side"]["D0_positive"]),
                "background_residual_linf": float(ev["r1_trace"]["closed_background"]["residual_linf"]),
            }
            for ev in sorted(evaluations, key=lambda row: float(row["tau"]))
        ],
    }
    bundle["content_hash"] = _stable_hash(bundle)
    bundle_path = run_root / "bundles" / "patha_b2c_calibration_bundle.json"
    _write_json(bundle_path, bundle)
    write_report(report_path=report_path, calibration=bundle, bundle_path=bundle_path)
    return {"bundle_path": str(bundle_path), "report_path": str(report_path), "status": status}


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.12e}"
    return str(value)


def _fmt_approx(value: Any, digits: int = 2) -> str:
    if value is None:
        return "-"
    return f"{float(value):.{digits}e}"


def _markdown_table(headers: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(header)) for header in headers) + " |")
    return "\n".join(lines)


def write_report(
    *,
    report_path: Path,
    calibration: Mapping[str, Any],
    bundle_path: Path,
) -> None:
    floor = calibration["tau_floor"]
    floor_eval = calibration["resolved_floor_evaluation"]
    obs = floor_eval["observables"]
    natural = floor_eval["naturalness"]
    frozen = calibration["frozen_negative_control"]
    errors = calibration["error_bars"]["propagated_at_resolved_floor"]
    drift = calibration["drift"]
    edge = calibration.get("edge_diagnostic")
    edge_spread = calibration.get("edge_estimate_fit_spread")
    spot = calibration.get("spot_convergence")
    p0_rows = []
    n0_over_target_values: list[float] = []
    k_values: list[float] = []
    kappa_from_records = float(floor_eval["direct_coefficients"]["K"]) / float(floor_eval["tau"])
    for row in sorted(calibration["evaluations"], key=lambda item: float(item["tau"]), reverse=True):
        p0 = float(row["R_norm"]) + TARGET_P0
        orders = math.log10(TARGET_P0 / p0) if p0 > 0.0 else math.inf
        n0_over_target = p0 * float(row["D0"]) / TARGET_P0
        k_value = kappa_from_records * float(row["tau"])
        n0_over_target_values.append(n0_over_target)
        k_values.append(k_value)
        p0_rows.append(
            {
                "tau": f"{float(row['tau']):.6g}",
                "P0": f"{p0:.3e}",
                "orders_below_target": f"{orders:.2f}",
                "N0/10.8": f"{n0_over_target:.3e}",
                "K": f"{k_value:.3e}",
            }
        )
    lines: list[str] = []
    lines.append("# Path-A B2c Calibration")
    lines.append("")
    lines.append(f"Overall B2c status: `{calibration['status']}`")
    lines.append(f"Machine bundle: `{bundle_path}`")
    lines.append(f"Bundle content hash: `{calibration['content_hash']}`")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(
        "Measured, fit-independent verdict: MISS. At every converged tau, the directly measured `P0` is far below "
        f"the target `54/5={TARGET_P0:.1f}` on the stable side; this is measured, not extrapolated."
    )
    lines.append("")
    lines.append(_markdown_table(["tau", "P0", "orders_below_target", "N0/10.8", "K"], p0_rows))
    lines.append("")
    lines.append(
        f"The measured deficit is {min(float(row['orders_below_target']) for row in p0_rows):.1f}-"
        f"{max(float(row['orders_below_target']) for row in p0_rows):.1f} decimal orders, i.e. the "
        "~7-to-9 order deficit is the headline result. Algebraically `P0=N0/D0`, so any stable root must "
        "satisfy `D0(tau*)=N0(tau*)/10.8`. In the measured records `N0/10.8` is only "
        f"`{min(n0_over_target_values):.3e}` to `{max(n0_over_target_values):.3e}` while `K` is "
        f"`{min(k_values):.3e}` to `{max(k_values):.3e}`; reaching the target therefore requires "
        "`D0 << K`, a knife-edge cancellation, regardless of the precise fitted tau*."
    )
    lines.append(
        "The `tau=1` structural win, meaning no fragile cancellation, does not survive calibration. The GR anchor "
        "can only be approached by riding the `D0 -> 0` stability boundary, so this is not a natural calibration."
    )
    if edge_spread is not None:
        estimates = {row["fit"]: row for row in edge_spread["estimates"]}
        frozen_est = estimates.get("frozen_deepest_coefficients", {})
        wide_est = estimates.get("wide_pair_power_law_tau_1_to_0p03", {})
        close_est = estimates.get("close_pair_power_law_deepest_two", {})
        lines.append(
            f"The rigorous edge statement is only the upper bound `tau* < tau_floor ~= "
            f"{float(floor['tau_floor']):.3g}`. Drift-aware estimates are fit-dependent: frozen deepest "
            f"coefficients give `tau* ~ {_fmt_approx(frozen_est.get('tau_estimate'))}`, the wide "
            f"`tau=1` to `tau=0.03` power-law pair gives `tau* ~ {_fmt_approx(wide_est.get('tau_estimate'))}`, "
            f"and the deepest close-pair fit is `{close_est.get('status', 'unavailable')}` "
            f"(`tau* ~ {_fmt_approx(close_est.get('tau_estimate'))}`). Thus tau* is not pinned by this run."
        )
    lines.append(calibration["finding"])
    if edge is not None:
        lines.append(
            f"Numerical-floor evidence: {edge['evidence']}"
        )
    if calibration.get("root_agreement") is None:
        lines.append(
            "The real dual-engine root finder was attempted on the real converged sample set, but no "
            "stable-side sign bracket exists, so there is no `root_agreement` to report."
        )
    else:
        agree = calibration["root_agreement"]
        lines.append(
            f"Real root-finder agreement: primary tau `{agree['primary_tau']:.12e}`, "
            f"secondary tau `{agree['secondary_tau']:.12e}`, passed `{agree['passed']}`."
        )
    if floor["tau_floor_row"] is not None:
        floor_tau = float(floor["tau_floor_row"]["tau"])
        lines.append(
            f"Numerics status: the requested `1e-3`-`1e-2` band was not reached; the deepest converged "
            f"tau is `{floor_tau:.6g}`, about `{floor_tau / 1.0e-2:.1f}x` above the top of that band. "
            "Warm-start machinery is functional but converted no new converged point here: the deepest full "
            "evaluation converges cold, and the floor descent from `0.03` to the current floor is only about "
            f"`{(0.03 - floor_tau) / 0.03:.1%}`. This is a known limitation, not a material improvement; pinning "
            "tau* would require a small-tau preconditioner/linear-solver effort, deferred and not required for "
            "the MISS verdict."
        )
    lines.append("")
    lines.append("No background, BdG, or Maxwell coefficient bundle is extrapolated below the recorded floor.")
    lines.append("")
    lines.append("## Tau Floor")
    lines.append("")
    if floor["tau_floor_row"] is not None:
        row = floor["tau_floor_row"]
        lines.append(
            f"Lowest converged closed-background probe/evaluation: `tau_floor={float(row['tau']):.12e}` "
            f"with residual `{_fmt(row['residual_linf'])}` from `{row['source']}`."
        )
    lines.append("")
    lines.append("Converged floor/probe rows:")
    lines.append("")
    lines.append(_markdown_table(["tau", "source", "residual_linf", "message"], floor["converged_rows"]))
    lines.append("")
    lines.append("Failed or dropped tau rows:")
    lines.append("")
    lines.append(_markdown_table(["tau", "source", "residual_linf", "message"], floor["failed_rows"]))
    lines.append("")
    lines.append("## Real Resolved Evaluations")
    lines.append("")
    lines.append(
        _markdown_table(
            ["tau", "kind", "R_norm", "D0", "stable", "background_residual_linf"],
            calibration["evaluations"],
        )
    )
    lines.append("")
    lines.append("## Frozen Negative Control")
    lines.append("")
    lines.append(
        f"Frozen-background root from the validated `tau=1` coefficients: "
        f"`tau_frozen ~ {_fmt_approx(frozen['tau_frozen'])}`; the frozen Schur edge is also "
        f"`~{_fmt_approx(frozen['tau_crit_frozen'])}`."
    )
    lines.append(
        f"That prior is below the resolved floor: `{frozen['below_tau_floor']}`. "
        "It is recorded only as the R6 negative control and is not used by the final status logic."
    )
    analytic_kappa = pe.hooke_kappa_hat()
    record_kappa = float(floor_eval["direct_coefficients"]["K"]) / float(floor_eval["tau"])
    lines.append(
        f"Kappa provenance note: local analytic formulas use continuum `kappa_hat={analytic_kappa:.10f}`; "
        f"direct `K/tau` in the B2a-seeded records is `{record_kappa:.10f}` "
        f"(relative difference `{abs(record_kappa - analytic_kappa) / analytic_kappa:.2e}`, discretization scale)."
    )
    if drift is not None:
        lines.append(
            f"Measured resolved coefficient drift from `tau={drift['tau_high']:.12e}` to "
            f"`tau={drift['tau_low']:.12e}` has max relative movement "
            f"`{drift['max_rel']:.12e}` across `{{B,Z,N}}`."
        )
        lines.append("")
        lines.append(
            _markdown_table(
                ["coefficient", "relative_drift"],
                [
                    {"coefficient": key, "relative_drift": value}
                    for key, value in drift["per_coefficient_relative_to_high_tau"].items()
                ],
            )
        )
    lines.append("")
    lines.append("")
    lines.append("## Edge Diagnostic")
    lines.append("")
    if edge is None:
        lines.append("No residual-gated full bundle was available for an edge diagnostic.")
    else:
        deepest = edge["deepest_full_evaluation"]
        local_edge = edge["local_edge_estimate_from_deepest_coefficients"]
        failed = edge["nearest_failed_probe_below_floor"]
        lines.append(
            f"Deepest full bundle: `tau={deepest['tau']:.12e}`, "
            f"`R_norm={deepest['R_norm']:.12e}`, `D0={deepest['D0']:.12e}`, "
            f"`P0={deepest['P0']:.12e}`, `K/(B0+Z0)={deepest['K_over_B0_plus_Z0']:.12e}`, "
            f"closed residual `{deepest['background_residual_linf']:.12e}`."
        )
        if failed is not None:
            lines.append(
                f"Nearest failed probe below the floor: `tau={failed['tau']:.12e}`, "
                f"residual `{_fmt(failed['residual_linf'])}`, message `{failed['message']}`."
            )
        lines.append(
            f"Deepest frozen-coefficient local edge estimate, reported only as one fit in the spread block: "
            f"`tau* ~ {_fmt_approx(local_edge['tau_star_local_coefficients'])}`, "
            f"`D0(tau*) ~ {_fmt_approx(local_edge['D0_at_tau_star_local_coefficients'])}`, "
            f"`K/(B0+Z0) ~ {local_edge['K_over_B0_plus_Z0_at_tau_star_local_coefficients']:.4g}`, "
            f"cancellation digits `~{local_edge['leading_cancelled_digits_at_tau_star_local_coefficients']:.1f}`. "
            "This is fit-dependent and does not pin the edge."
        )
        if edge_spread is not None:
            lines.append("")
            lines.append("Fit-spread estimates recorded in the machine bundle:")
            lines.append("")
            lines.append(
                _markdown_table(
                    ["fit", "tau_estimate", "status"],
                    [
                        {
                            "fit": row["fit"],
                            "tau_estimate": None
                            if row.get("tau_estimate") is None
                            else _fmt_approx(row["tau_estimate"]),
                            "status": row.get("status"),
                        }
                        for row in edge_spread["estimates"]
                    ],
                )
            )
    lines.append("")
    lines.append("## Deepest Spot Convergence")
    lines.append("")
    if spot is None:
        lines.append("No deepest-τ spot convergence packet was available.")
    else:
        modal = spot.get("modal")
        maxwell = spot.get("maxwell_grid_refinement")
        if modal is not None:
            lines.append(
                f"BdG modal sweep at deepest τ: exported K=`{modal['exported_mode_count']}`, "
                f"all-positive modes `{modal['all_positive_mode_count']}`, max truncation "
                f"`{modal['max_truncation_error_at_export']:.12e}` against tolerance "
                f"`{modal['tolerance']:.12e}`, passed `{modal['passed']}`."
            )
        if maxwell is not None:
            lines.append(
                f"Maxwell grid refinement `{maxwell['production_grid']['nr']}x{maxwell['production_grid']['nw']}` "
                f"to `{maxwell['refined_grid']['nr']}x{maxwell['refined_grid']['nw']}` max relative Z/N movement "
                f"`{maxwell['max_rel']:.12e}`."
            )
            lines.append("")
            lines.append(
                _markdown_table(
                    ["coefficient", "relative_change"],
                    [
                        {"coefficient": key, "relative_change": value}
                        for key, value in maxwell["relative_change_vs_refined"].items()
                    ],
                )
            )
        if maxwell is None and modal is None:
            lines.append("No modal or Maxwell spot convergence data were found.")
    lines.append("")
    lines.append("## Naturalness At Bound")
    lines.append("")
    lines.append(
        f"At the deepest residual-gated upper bound `tau={float(floor_eval['tau']):.12e}`: "
        f"`|ln tau|={natural['abs_ln_tau']:.12e}`, "
        f"`K/(B0+Z0)={natural['K_over_B0_plus_Z0']:.12e}`, "
        f"`D0={natural['D0']:.12e}`, cancellation fraction "
        f"`{natural['cancellation_fraction']:.12e}` (`{natural['leading_cancelled_digits']:.3f}` digits)."
    )
    lines.append(
        "Because no real root is confirmed, this is bound/floor naturalness, not naturalness at a calibrated root."
    )
    lines.append("")
    lines.append("## Held-Out Quantities")
    lines.append("")
    lines.append(
        "No calibrated `tau*` is confirmed in this run, so `R_pole/P2/P4` are not emitted as calibrated held-out "
        "predictions. The deepest resolved full bundle is reported below as a diagnostic only."
    )
    lines.append("")
    lines.append(
        _markdown_table(
            ["tau", "R_norm", "D0", "R_pole", "P2", "P4"],
            [{"tau": floor_eval["tau"], **obs}],
        )
    )
    lines.append("")
    lines.append("## Error Bars")
    lines.append("")
    lines.append(
        "Recorded B2a/B2b budgets were propagated by coefficient perturbation at the deepest residual-gated "
        "bundle. `K` is exact for the frozen Hooke family; no error bar here converts the bound into a confirmed root, "
        "and no local tau* error bar is reported because tau* is not pinned by this run."
    )
    lines.append("")
    lines.append(
        _markdown_table(
            ["quantity", "value", "abs_error"],
            [
                {
                    "quantity": key,
                    "value": errors["base_observables"].get(key),
                    "abs_error": errors["absolute_errors"].get(key),
                }
                for key in OBSERVABLE_ERROR_KEYS
            ],
        )
    )
    lines.append(
        f"Recorded B2b cross-engine max relative spread (reported, not hidden): "
        f"`{calibration['error_bars']['recorded_b2b_dual_engine_max_rel']:.12e}`."
    )
    lines.append("")
    lines.append("## Dual Engine")
    lines.append("")
    lines.append(
        "B2c's new assembly was checked two ways: B1 `lane_extract` plus `observable_residuals`, and an independent "
        "SymPy expansion of the squared response-amplitude series through `x=omega^2` order two. Root finding, when bracketed, uses "
        "Brent on the primary model and independent bisection on the secondary model."
    )
    lines.append("")
    lines.append(
        _markdown_table(
            ["tau", "kind", "passed", "max_abs", "max_rel"],
            calibration["dual_engine"]["evaluation_rows"],
        )
    )
    lines.append("")
    lines.append("## Gates And Negative Controls")
    lines.append("")
    controls = calibration["negative_controls"]["negative_controls"]
    lines.append(
        _markdown_table(
            ["gate", "failed_as_expected", "wrong_answer"],
            [
                {
                    "gate": key,
                    "failed_as_expected": value["failed_as_expected"],
                    "wrong_answer": value["wrong_answer"],
                }
                for key, value in controls.items()
            ],
        )
    )
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("Commands used by this route, each intended to stay under `timeout 600`:")
    lines.append("")
    lines.append("```bash")
    lines.append("env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration seed-tau1")
    lines.append("timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.03 --background-grid 16,16")
    lines.append("timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.0295 --warm-start-background software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p03.json --warm-start-final-only --warm-start-wall-predictor")
    lines.append("timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.02925 --warm-start-background software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p0295.json --warm-start-final-only --warm-start-wall-predictor")
    lines.append("timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.029125 --warm-start-background software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p02925.json --warm-start-final-only --warm-start-wall-predictor")
    lines.append("timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration probe-background --tau 0.0290625 --warm-start-background software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p029125.json --warm-start-final-only --warm-start-wall-predictor")
    lines.append("timeout 600 env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration evaluate --tau 0.029125 --kind resolved_deepest --background-bundle software/stage1_solver/runs/patha_b2c_calibration/backgrounds/patha_b2c_background_tau_0p029125.json")
    lines.append("env PYTHONPATH=software/stage1_solver/src python -m stage1_solver.patha_b2c_calibration calibrate")
    lines.append("timeout 600 env PYTHONPATH=software/stage1_solver/src pytest -q software/stage1_solver/tests/test_patha_b2c_calibration.py")
    lines.append("```")
    lines.append("")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    sub = parser.add_subparsers(dest="command", required=True)

    seed = sub.add_parser("seed-tau1")
    seed.add_argument("--b2a-validated", type=Path, default=DEFAULT_B2A_VALIDATED)
    seed.add_argument("--b2b-validated", type=Path, default=DEFAULT_B2B_VALIDATED)

    probe = sub.add_parser("probe-background")
    probe.add_argument("--tau", type=float, required=True)
    probe.add_argument("--background-grid", type=_parse_grid, default=b2a.DEFAULT_BACKGROUND_GRID)
    probe.add_argument("--warm-start-background", type=Path)
    probe.add_argument("--warm-start-final-only", action="store_true")
    probe.add_argument("--warm-start-wall-predictor", action="store_true")
    probe.add_argument("--newton-jvp-epsilon", type=float)

    timeout = sub.add_parser("record-timeout")
    timeout.add_argument("--tau", type=float, required=True)
    timeout.add_argument("--seconds", type=float, default=600.0)

    evaluate = sub.add_parser("evaluate")
    evaluate.add_argument("--tau", type=float, required=True)
    evaluate.add_argument("--kind", default="resolved")
    evaluate.add_argument("--background-grid", type=_parse_grid, default=b2a.DEFAULT_BACKGROUND_GRID)
    evaluate.add_argument("--bdg-grid", type=_parse_grid, default=(10, 10))
    evaluate.add_argument("--maxwell-grid", type=_parse_grid, default=b2b.DEFAULT_FINAL_GRID)
    evaluate.add_argument("--profile-points", type=int, default=b2a.DEFAULT_PROFILE_POINTS)
    evaluate.add_argument("--modes", type=int, default=b2a.DEFAULT_BDG_MODES)
    evaluate.add_argument("--window", type=float, default=b2b.DEFAULT_FINAL_WINDOW)
    evaluate.add_argument("--radial-scale", type=float, default=b2b.DEFAULT_FINAL_TRUNCATION)
    evaluate.add_argument("--warm-start-background", type=Path)
    evaluate.add_argument("--warm-start-final-only", action="store_true")
    evaluate.add_argument("--warm-start-wall-predictor", action="store_true")
    evaluate.add_argument("--background-bundle", type=Path)
    evaluate.add_argument("--newton-jvp-epsilon", type=float)

    calibrate = sub.add_parser("calibrate")
    calibrate.add_argument("--report-path", type=Path, default=DEFAULT_REPORT_PATH)
    calibrate.add_argument("--b2a-validated", type=Path, default=DEFAULT_B2A_VALIDATED)
    calibrate.add_argument("--b2b-validated", type=Path, default=DEFAULT_B2B_VALIDATED)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.command == "seed-tau1":
        path, bundle = seed_tau1_from_validated(
            run_root=args.run_root,
            b2a_validated_path=args.b2a_validated,
            b2b_validated_path=args.b2b_validated,
        )
        print(json.dumps({"path": str(path), "tau": bundle["tau"], "R_norm": bundle["observables"]["R_norm"]}, sort_keys=True))
        return 0
    if args.command == "probe-background":
        path, record = probe_background(
            tau=args.tau,
            run_root=args.run_root,
            background_grid=args.background_grid,
            warm_start_background=args.warm_start_background,
            warm_start_final_only=args.warm_start_final_only,
            warm_start_wall_predictor=args.warm_start_wall_predictor,
            newton_jvp_epsilon=args.newton_jvp_epsilon,
        )
        print(json.dumps({"path": str(path), "converged": record["converged"], "residual_linf": record["residual_linf"]}, sort_keys=True))
        return 0
    if args.command == "record-timeout":
        path, record = record_timeout_probe(tau=args.tau, seconds=args.seconds, run_root=args.run_root)
        print(json.dumps({"path": str(path), "converged": record["converged"], "timeout": True}, sort_keys=True))
        return 0
    if args.command == "evaluate":
        path, bundle = evaluate_tau(
            tau=args.tau,
            run_root=args.run_root,
            kind=args.kind,
            background_grid=args.background_grid,
            bdg_grid=args.bdg_grid,
            maxwell_grid=args.maxwell_grid,
            profile_points=args.profile_points,
            modes=args.modes,
            window=args.window,
            radial_scale=args.radial_scale,
            warm_start_background=args.warm_start_background,
            warm_start_final_only=args.warm_start_final_only,
            warm_start_wall_predictor=args.warm_start_wall_predictor,
            background_bundle=args.background_bundle,
            newton_jvp_epsilon=args.newton_jvp_epsilon,
        )
        print(json.dumps({"path": str(path), "tau": bundle["tau"], "R_norm": bundle["observables"]["R_norm"]}, sort_keys=True))
        return 0
    if args.command == "calibrate":
        result = calibrate_from_records(
            run_root=args.run_root,
            report_path=args.report_path,
            b2a_validated_path=args.b2a_validated,
            b2b_validated_path=args.b2b_validated,
        )
        print(json.dumps(result, sort_keys=True))
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
