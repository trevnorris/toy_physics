"""M1c physical frozen-run driver.

This module keeps the physical freeze and the numerical run separate:
``freeze`` writes the content-only freeze sheet and hash; ``solve`` refuses to
run unless that frozen sheet already exists and hashes back to its directory.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import torch

from .backend import configure_backend, library_versions
from .config import BackendConfig, BranchSmokeConfig, HarnessConfig, NewtonConfig, source_revision
from .coupled_branch import run_branch_continuation
from .grid import TensorProductGrid
from .m1c_background_export import _make_bundle
from .config import TensorGridSpec


PRIMARY_GRID = (10, 8)
CONVERGENCE_LEVELS = ((6, 4), (8, 6), PRIMARY_GRID)
RESIDUAL_CEILING = 5.0e-8
M1C_PREP_N0_INTERP_RELATIVE_FLOOR = 0.015591236611409905
FROZEN_ROOT = Path("software/stage1_solver/frozen/m1c")
RUN_ROOT = Path("software/stage1_solver/runs/m1c_physical")
REPORT_PATH = "software/stage1_solver/reports/m1c_physical_run.md"
SCHEMA = "stage1_m1c_physical_wp1_background/v1"


def canonical_json(data: Any) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_root() -> Path:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        if out:
            return Path(out)
    except (OSError, subprocess.CalledProcessError):
        pass
    return Path(__file__).resolve().parents[4]


def _git_status_short(root: Path) -> list[str]:
    try:
        out = subprocess.check_output(
            ["git", "-C", str(root), "status", "--short"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return ["unavailable"]
    return [line for line in out.splitlines() if line.strip()]


def _source_hashes(root: Path) -> dict[str, str]:
    rel_paths = [
        "software/stage1_solver/src/stage1_solver/config.py",
        "software/stage1_solver/src/stage1_solver/coupled_branch.py",
        "software/stage1_solver/src/stage1_solver/m1c_background_export.py",
        "software/stage1_solver/src/stage1_solver/m1c_physical_run.py",
        "software/stage1_solver/mathematica/mt15_06_m1c_prep_crossengine.wls",
        "research/pde_audit/scripts/stage_v2_22b_solver_handoff_validator.py",
        "research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py",
        "research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py",
        "research/pde_audit/scripts/stage_v2_22c_end_to_end_smoke_pipeline.py",
    ]
    hashes: dict[str, str] = {}
    for rel in rel_paths:
        path = root / rel
        hashes[rel] = sha256_file(path) if path.exists() else "missing"
    return hashes


def frozen_branch_config(grid: tuple[int, int] = PRIMARY_GRID) -> BranchSmokeConfig:
    """Return the decision-07 physical branch config.

    Values not named as physics targets in decision-07 are still frozen here as
    solver/closure controls so they cannot drift after residual evaluation.
    """

    return BranchSmokeConfig(
        name="m1c_decision_07_physical_frozen_effective_closure",
        placeholder_label=(
            "PHYSICAL M1c decision-07 structured effective-closure branch; "
            "target-blind frozen values; not BranchSmokeConfig placeholders"
        ),
        r_max=2.0,
        w_min=0.0,
        w_max=37.0 / 20.0,
        solve_grid=grid,
        ladder_levels=CONVERGENCE_LEVELS,
        hbar=1.0,
        particle_mass=1.0,
        gauge_charge=0.35,
        mu0=1.0,
        xi=1.0,
        continuation_K_values=(0.05, 0.15, 0.3, 0.5),
        mass=1.0,
        localization_width=0.75,
        localization_floor=0.8,
        localization_amplitude=0.45,
        geometry_profile="cubic_smoothstep",
        r_mouth=1.0,
        r_exit=1.5,
        geometry_decay_length=1.0,
        radial_wall_strength=0.65,
        axial_trap_strength=0.12,
        matter_mouth_boundary="dirichlet",
        a0_mouth_boundary="dirichlet",
        matter_exit_impedance_alpha=0.25,
        a0_radial_impedance_alpha=0.5,
        a0_exit_impedance_alpha=0.4,
        sponge_enabled=False,
        sponge_width=0.0,
        sponge_matter_strength=0.0,
        sponge_gauge_strength=0.0,
        sponge_power=2,
        initial_mu=1.0,
        max_ladder_level_seconds=120.0,
        newton=NewtonConfig(
            residual_atol=5.0e-8,
            residual_rtol=5.0e-9,
            step_atol=1.0e-11,
            step_rtol=1.0e-10,
            max_newton_iters=16,
            gmres_rtol=5.0e-8,
            gmres_atol=1.0e-10,
            gmres_restart=96,
            gmres_maxiter=12,
            max_line_search_iters=18,
            accept_best_line_search_decrease=True,
            finite_difference_jvp_epsilon=1.0e-5,
        ),
    )


def harness_config(grid: tuple[int, int], run_root: Path) -> HarnessConfig:
    return HarnessConfig(
        backend=BackendConfig(),
        branch=frozen_branch_config(grid),
        run_root=str(run_root),
        report_path=REPORT_PATH,
    )


def build_freeze_sheet() -> dict[str, Any]:
    root = repo_root()
    branch = frozen_branch_config(PRIMARY_GRID)
    return {
        "schema": "m1c_gate_a_decision_07_freeze/v1",
        "parent_action_status": "effective_closure",
        "branch_identity": {
            "name": "m1c_decision_07_structured_physical_falsification_run",
            "decision_record": "software/stage1_solver/decisions/07_gate_a_freeze_sheet.md",
            "gate": "GATE-A GO",
            "claim_scope": "pre_registered_effective_closure_branch",
            "m1c_physical_frozen_run": True,
            "target_blind": True,
            "no_post_residual_refit": True,
            "physical_export_permitted_guard": "not_touched",
            "firewalled_simulation_tree": "not_touched",
        },
        "geometry": {
            "a": 1.0,
            "R_mouth": 1.0,
            "L": 37.0 / 20.0,
            "L_exact": "37/20",
            "R_exit": 3.0 / 2.0,
            "R_exit_exact": "3/2",
            "R0_family": "cubic_smoothstep",
            "R0_formula": "R0(w)=1+(1/2)*(3*x^2-2*x^3), x=w/L",
            "boundary_class": "mouth Dirichlet / open-impedance(Robin) exit / no hard cap",
            "epsilon_r": 1.0 / 20.0,
            "ell_radial_support": 1.0 / 20.0,
            "radial_support_layer_active": False,
        },
        "wall_constitutive": {
            "mu_eta": "1",
            "T_w": "1",
            "T_Omega": "1",
            "K_eta": "1",
            "positivity": {"mu_eta_gt_0": True, "T_w_gt_0": True, "K_eta_plus_6_T_Omega": 7.0},
        },
        "constants": {
            "G": 1.0,
            "c": 1.0,
            "c_s": 1.0,
            "mhat0": 1.0,
            "S_port": 1.0,
            "theta_tail": 1.0,
            "tail_sector_active": False,
            "gauge": "H=Z",
            "EOS": "h=(5*K/4)*rho^4",
        },
        "source_map": {
            "chi_Q_input_frozen": False,
            "relation_frozen": "N_Q=chi_Q^-1",
            "chi_Q_extraction": "chi_Q=P0_target/(mhat0^2*S_port*(N0/D0))",
            "N_Q_extraction": "N_Q=mhat0^2*S_port*(N0/D0)/P0_target",
        },
        "extraction_formulas": {
            "D0": "K-B0-Z0",
            "D2": "-(M+B2+Z2)",
            "D4": "-(B4+Z4)",
            "P0": "N0/D0",
            "P2": "(D0*N2-2*D2*N0)/D0^2",
            "P4": "(D0^2*N4-2*D0*(D2*N2+D4*N0)+3*D2^2*N0)/D0^3",
            "R_pole": "D0*(B4+Z4)-3*(M+B2+Z2)^2",
            "R_norm": "mhat0^2*S_port*N0/D0-54*G*c_s^5/(5*a^5*c^5)",
            "P0_target": "54*G*c_s^5/(5*a^5*c^5)",
        },
        "solver": {
            "torch_wp1_entry_point": "stage1_solver.coupled_branch.run_branch_continuation",
            "background_export_entry_point": "stage1_solver.m1c_physical_run",
            "mathematica_chain": "software/stage1_solver/mathematica/mt15_06_m1c_prep_crossengine.wls",
            "v2_chain": "research/pde_audit/scripts/stage_v2_22c_end_to_end_smoke_pipeline.py",
            "backend": asdict(BackendConfig()),
            "primary_grid": list(PRIMARY_GRID),
            "convergence_levels": [list(level) for level in CONVERGENCE_LEVELS],
            "residual_ceiling": RESIDUAL_CEILING,
            "branch_config": asdict(branch),
        },
        "validation_protocol": {
            "freeze_before_solve": True,
            "frozen_packet_tracked_path": "software/stage1_solver/frozen/m1c/<freeze_hash>/",
            "transient_output_path": "software/stage1_solver/runs/m1c_physical/<freeze_hash>/",
            "torch_residual_required_linf_lt": RESIDUAL_CEILING,
            "wp1_background_resolution_estimate": "three-level target-blind surrogate Richardson estimate",
            "derived_chain_checks": [
                "current_frechet_matches_step8c",
                "outgoing_flux_positive",
                "open_not_hard_cap",
                "pure_gauge_zero_physical_transfer",
                "basis_invariance",
                "v2_09_regression",
                "green_residuals_small",
                "bdg_residuals_small",
                "N0_positive",
            ],
            "forbidden_field_scan": "stage_v2_22b recursive target-field rejection with injected R_norm packet",
            "byte_reproducibility": "two reruns must reproduce freeze hash, packet hash, observable hash, and R_norm",
            "budget_scope": "numerical-only; excludes free_choice/posit uncertainty and model-form uncertainty",
        },
        "source_revision": {
            "git_head": source_revision(),
            "working_tree_status_at_freeze": _git_status_short(root),
            "source_file_sha256": _source_hashes(root),
            "library_versions": library_versions(),
        },
    }


def write_freeze() -> tuple[str, Path]:
    root = repo_root()
    sheet = build_freeze_sheet()
    text = canonical_json(sheet)
    freeze_hash = sha256_text(text)
    out_dir = root / FROZEN_ROOT / freeze_hash
    out_dir.mkdir(parents=True, exist_ok=True)
    freeze_path = out_dir / "freeze_sheet.json"
    freeze_path.write_text(text, encoding="utf-8")
    (out_dir / "freeze_hash.txt").write_text(freeze_hash + "\n", encoding="utf-8")
    return freeze_hash, freeze_path


def load_freeze(freeze_hash: str) -> tuple[dict[str, Any], Path]:
    root = repo_root()
    freeze_path = root / FROZEN_ROOT / freeze_hash / "freeze_sheet.json"
    if not freeze_path.exists():
        raise FileNotFoundError(f"freeze sheet not found: {freeze_path}")
    text = freeze_path.read_text(encoding="utf-8")
    actual = sha256_text(text)
    if actual != freeze_hash:
        raise ValueError(f"freeze hash mismatch: dir={freeze_hash} content={actual}")
    return json.loads(text), freeze_path


def _full_content_hash(payload: dict[str, Any]) -> str:
    return sha256_text(canonical_json({k: v for k, v in payload.items() if k != "content_hash"}))


def solve_background(freeze_hash: str, grid: tuple[int, int]) -> dict[str, Any]:
    root = repo_root()
    _, freeze_path = load_freeze(freeze_hash)
    run_root = root / RUN_ROOT / freeze_hash / f"wp1_{grid[0]}x{grid[1]}"
    config = harness_config(grid, run_root)
    dtype = configure_backend(config.backend)
    cfg = config.branch
    tensor_grid = TensorProductGrid.create(
        TensorGridSpec(
            r_max=cfg.r_max,
            nr=grid[0],
            w_min=cfg.w_min,
            w_max=cfg.w_max,
            nw=grid[1],
        ),
        dtype=dtype,
        device=config.backend.device,
    )
    state, summary = run_branch_continuation(
        config,
        tensor_grid,
        grid_name=f"m1c_physical_nr_{grid[0]}_nw_{grid[1]}",
    )
    bundle = _make_bundle(config=config, grid=tensor_grid, state=state, summary=summary)
    bundle["schema"] = SCHEMA
    bundle["scope"] = (
        "M1c physical frozen decision-07 effective-closure WP1 background; "
        "freeze was written and hashed before solve"
    )
    bundle["freeze_hash"] = freeze_hash
    bundle["freeze_sheet_sha256"] = freeze_hash
    bundle["freeze_sheet_path"] = str(freeze_path.relative_to(root))
    bundle["m1c_physical_frozen_run"] = True
    bundle["claim_scope"] = "pre_registered_effective_closure_branch"
    bundle["content_hash"] = _full_content_hash(bundle)
    residual = float(bundle["residuals"]["coupled_stationary_linf"])
    if not summary["converged"] or residual >= RESIDUAL_CEILING:
        raise RuntimeError(
            f"WP1 solve failed frozen residual floor at {grid}: "
            f"converged={summary['converged']} residual={residual:.6e}"
        )
    return bundle


def _estimate_triplet(level_rows: list[dict[str, Any]], observable: str) -> dict[str, Any]:
    vals = [float(row["surrogate_values"][observable]) for row in level_rows]
    hs = [float(row["h"]) for row in level_rows]
    d01 = abs(vals[0] - vals[1])
    d12 = abs(vals[1] - vals[2])
    order = None
    if d01 > 0.0 and d12 > 0.0 and hs[0] != hs[1]:
        order_value = math.log(d01 / d12) / math.log(hs[0] / hs[1])
        if math.isfinite(order_value):
            order = order_value
    if order is not None and order > 0.0 and hs[1] > hs[2]:
        denom = (hs[1] / hs[2]) ** order - 1.0
        finest_error = d12 / denom if denom > 0.0 else d12
    else:
        finest_error = d12
    rel = finest_error / max(abs(vals[2]), 1.0e-30)
    return {
        "observable": observable,
        "values": vals,
        "h_values": hs,
        "observed_order": order,
        "finest_error_estimate": finest_error,
        "finest_relative_error_estimate": rel,
    }


def write_solve_outputs(freeze_hash: str) -> dict[str, Any]:
    root = repo_root()
    out_dir = root / FROZEN_ROOT / freeze_hash
    load_freeze(freeze_hash)
    bundles = [solve_background(freeze_hash, level) for level in CONVERGENCE_LEVELS]
    paths: list[str] = []
    for bundle in bundles:
        nr = int(bundle["grid"]["nr"])
        nw = int(bundle["grid"]["nw"])
        path = out_dir / f"wp1_background_{nr}x{nw}.json"
        path.write_text(json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        paths.append(str(path.relative_to(root)))
    level_rows = []
    for bundle in bundles:
        grid = bundle["grid"]
        level_rows.append(
            {
                "grid": [grid["nr"], grid["nw"]],
                "h": max(float(grid["dr"]), float(grid["dw"])),
                "residual_linf": bundle["residuals"]["coupled_stationary_linf"],
                "content_hash": bundle["content_hash"],
                "surrogate_values": bundle["convergence_evidence"]["surrogate_values"],
            }
        )
    observable_names = list(level_rows[-1]["surrogate_values"].keys())
    estimates = [_estimate_triplet(level_rows, name) for name in observable_names]
    governing = max(
        (
            row
            for row in estimates
            if row["observable"] not in {"density_mass", "final_residual_linf"}
        ),
        key=lambda row: row["finest_relative_error_estimate"],
    )
    convergence = {
        "schema": "m1c_wp1_background_resolution_estimate/v1",
        "freeze_hash": freeze_hash,
        "levels": level_rows,
        "estimates": estimates,
        "governing_background_resolution_relative_floor": governing["finest_relative_error_estimate"],
        "governing_observable": governing["observable"],
        "scope": "target-blind WP1 raw-field surrogate background-resolution estimate; not a model-form uncertainty",
    }
    convergence_path = out_dir / "wp1_background_resolution.json"
    convergence_path.write_text(json.dumps(convergence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema": "m1c_physical_wp1_solve_manifest/v1",
        "freeze_hash": freeze_hash,
        "freeze_sheet_path": str((out_dir / "freeze_sheet.json").relative_to(root)),
        "background_paths": paths,
        "primary_background_path": paths[-1],
        "primary_background_content_hash": bundles[-1]["content_hash"],
        "primary_residual_linf": bundles[-1]["residuals"]["coupled_stationary_linf"],
        "wp1_background_resolution_path": str(convergence_path.relative_to(root)),
        "wp1_background_resolution_hash": _full_content_hash(convergence),
    }
    manifest_path = out_dir / "wp1_solve_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def summarize(freeze_hash: str) -> dict[str, Any]:
    root = repo_root()
    out_dir = root / FROZEN_ROOT / freeze_hash
    freeze, _ = load_freeze(freeze_hash)
    packet_path = out_dir / "m1c_v2_22b_physical_frozen_packet.json"
    diagnostics_path = out_dir / "m1c_physical_derived_diagnostics.json"
    observable_path = out_dir / "v2_22c_observable_packet.json"
    pipeline_path = out_dir / "v2_22c_pipeline_report.json"
    tolerance_path = out_dir / "v2_22c_tolerance_budget.json"
    forbidden_path = out_dir / "forbidden_R_norm_validation_report.json"
    convergence_path = out_dir / "wp1_background_resolution.json"
    packet = _load_json(packet_path)
    diagnostics = _load_json(diagnostics_path)
    observable = _load_json(observable_path)
    pipeline = _load_json(pipeline_path)
    tolerance = _load_json(tolerance_path)
    forbidden = _load_json(forbidden_path)
    convergence = _load_json(convergence_path)
    coeff = diagnostics["direct_coefficients"]
    constants = packet["constants"]
    d0 = coeff["K"] - coeff["B0"] - coeff["Z0"]
    p0 = coeff["N0"] / d0
    target = 54.0 * constants["G"] * constants["c_s"] ** 5 / (
        5.0 * constants["a"] ** 5 * constants["c"] ** 5
    )
    transfer = constants["mhat0"] ** 2 * constants["S_port"] * p0
    chi_q = target / transfer
    n_q = 1.0 / chi_q
    branch = observable["branches"][0]
    interp = diagnostics["interpolation_method_sensitivity"]
    interp_n0 = interp["relative_deltas_order2_vs_order1"]["N0"]["relative_abs_delta"]
    interp_n0_budget = max(float(interp_n0), M1C_PREP_N0_INTERP_RELATIVE_FLOOR)
    mesh_rows = diagnostics["spike2_transfer"]["mesh_convergence_rows"]
    mesh_n0_rel = abs(mesh_rows[-1]["N0"] - mesh_rows[-2]["N0"]) / max(abs(mesh_rows[-1]["N0"]), 1.0e-30)
    mesh_z0_rel = abs(mesh_rows[-1]["Z0"] - mesh_rows[-2]["Z0"]) / max(abs(mesh_rows[-1]["Z0"]), 1.0e-30)
    operator = diagnostics["derived_maxwell_transfer"]["operator_gauge_residual_metrics"]
    rel_components = {
        "background_resolution": convergence["governing_background_resolution_relative_floor"],
        "interp_method_N0_budget": interp_n0_budget,
        "transfer_mesh_N0_successive": mesh_n0_rel,
        "transfer_mesh_Z0_successive": mesh_z0_rel,
        "continuity": operator["max_continuity_residual"],
        "boundary_hardcap_contrast": operator["hard_cap_to_open_boundary_norm_ratio"],
    }
    combined_relative_rss = math.sqrt(sum(float(v) ** 2 for v in rel_components.values()))
    combined_relative_max = max(float(v) for v in rel_components.values())
    abs_rnorm_budget_rss = combined_relative_rss * abs(transfer)
    abs_rnorm_budget_max = combined_relative_max * abs(transfer)
    summary = {
        "schema": "m1c_physical_frozen_summary/v1",
        "freeze_hash": freeze_hash,
        "branch_identity": freeze["branch_identity"],
        "tracked_artifact_dir": str(out_dir.relative_to(root)),
        "hashes": {
            "freeze_sheet_sha256": freeze_hash,
            "packet_content_hash": packet["content_hash"],
            "packet_file_sha256": sha256_file(packet_path),
            "diagnostics_content_hash": diagnostics["content_hash"],
            "diagnostics_file_sha256": sha256_file(diagnostics_path),
            "observable_packet_hash": pipeline["observable_packet_hash"],
            "observable_file_sha256": sha256_file(observable_path),
            "tolerance_budget_hash": pipeline["hashes"]["tolerance_budget_hash"],
        },
        "derived_coefficients": coeff,
        "physical_observables": {
            "D0": d0,
            "P0": p0,
            "P0_transfer": transfer,
            "P0_target": target,
            "R_norm": branch["residuals"]["R_norm"],
            "R_pole": branch["residuals"]["R_pole"],
            "P2": branch["residuals"]["R_P2"],
            "P4": branch["residuals"]["R_P4"],
            "chi_Q": chi_q,
            "N_Q": n_q,
            "chi_Q_minus_1": chi_q - 1.0,
        },
        "gates": {
            "v2_mechanical_pipeline_pass": pipeline["mechanical_pipeline_pass"],
            "open_gate_pass": branch["pass_flags"]["open_gate_pass"],
            "stability_gate_pass": branch["pass_flags"]["stability_gate_pass"],
            "target_packet_pass": branch["pass_flags"]["target_packet_pass"],
            "forbidden_R_norm_rejected": (forbidden["validation_pass"] is False)
            and forbidden["error_count"] > 0,
            "forbidden_issues": forbidden["issues"],
        },
        "wp1": {
            "primary_residual_linf": diagnostics["torch_background"]["residual_linf"],
            "primary_grid": diagnostics["torch_background"]["grid"],
            "background_content_hash": diagnostics["torch_background"]["content_hash"],
            "background_resolution": convergence,
        },
        "section_J_numerical_error_budget": {
            "scope": (
                "NUMERICAL ONLY. Excludes free_choice/posit uncertainty and model-form uncertainty; "
                "conditional on the frozen target-blind decision-07 effective-closure posits."
            ),
            "components_relative": rel_components,
            "interp_method_detail": {
                "physical_order2_vs_order1_N0_relative": interp_n0,
                "m1c_prep_carried_N0_relative_floor": M1C_PREP_N0_INTERP_RELATIVE_FLOOR,
                "budget_uses_max": interp_n0_budget,
            },
            "solver_newton_linf_absolute": diagnostics["torch_background"]["residual_linf"],
            "green_residual_absolute": operator["max_green_residual"],
            "bdg_eigen_residual_absolute": operator["max_bdg_eigen_residual"],
            "bdg_amp_residual_absolute": operator["max_bdg_amp_residual"],
            "bdg_phase_residual_absolute": operator["max_bdg_phase_residual"],
            "combined_relative_rss": combined_relative_rss,
            "combined_relative_max": combined_relative_max,
            "absolute_R_norm_budget_rss": abs_rnorm_budget_rss,
            "absolute_R_norm_budget_max": abs_rnorm_budget_max,
            "resolution_scope": "modest CPU run; tighter background and transfer ladders are deferred to GPU/cluster.",
        },
        "v2_tolerance_budget": tolerance,
    }
    summary["content_hash"] = _full_content_hash(summary)
    summary_path = out_dir / "m1c_physical_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="M1c physical frozen-run helper")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("freeze")
    solve = sub.add_parser("solve")
    solve.add_argument("freeze_hash")
    summ = sub.add_parser("summarize")
    summ.add_argument("freeze_hash")
    args = parser.parse_args()

    if args.cmd == "freeze":
        freeze_hash, freeze_path = write_freeze()
        print(f"freeze_hash: {freeze_hash}")
        print(f"freeze_sheet: {freeze_path}")
        return 0
    if args.cmd == "solve":
        manifest = write_solve_outputs(args.freeze_hash)
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0
    if args.cmd == "summarize":
        summary = summarize(args.freeze_hash)
        print(json.dumps(summary["physical_observables"], indent=2, sort_keys=True))
        print(f"summary_hash: {summary['content_hash']}")
        return 0
    raise AssertionError(args.cmd)


if __name__ == "__main__":
    raise SystemExit(main())
