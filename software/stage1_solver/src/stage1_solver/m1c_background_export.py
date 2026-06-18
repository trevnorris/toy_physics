"""M1c-prep target-blind WP1 background exporter.

This exporter deliberately reuses the validated Step-3 coupled branch
continuation machinery.  It writes a Mathematica-facing background bundle for
M1c-prep only: placeholder parameters, no Gate A freeze, no physical claim.
"""

from __future__ import annotations

from dataclasses import asdict
import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .backend import configure_backend
from .config import HarnessConfig, TensorGridSpec, source_revision, stable_hash, stable_json
from .coupled_branch import (
    branch_boundary_conditions,
    coupled_branch_residual,
    geometry_radius_torch,
    localization_weight_torch,
    run_branch_continuation,
    unpack_coupled_fields,
)
from .grid import TensorProductGrid


DEFAULT_RUN_ROOT = "software/stage1_solver/runs/m1c_prep_background_export"
DEFAULT_REPORT_PATH = "software/stage1_solver/reports/m1c_prep_background_export.md"
DEFAULT_RESIDUAL_CEILING = 5.0e-8


def _numpy_list(value: torch.Tensor | np.ndarray) -> Any:
    if isinstance(value, torch.Tensor):
        value = value.detach().cpu().numpy()
    return np.asarray(value, dtype=np.float64).tolist()


def _stage_evidence(summary: dict[str, Any]) -> list[dict[str, Any]]:
    """Keep deterministic convergence evidence; exclude timing/RSS noise."""

    rows: list[dict[str, Any]] = []
    for stage in summary["stages"]:
        rows.append(
            {
                "eos_K": stage["eos_K"],
                "converged": stage["converged"],
                "iterations": stage["iterations"],
                "initial_residual_norm": stage["initial_residual_norm"],
                "final_residual_norm": stage["final_residual_norm"],
                "tolerance": stage["tolerance"],
                "message": stage["message"],
                "gmres_iterations": stage["gmres_iterations"],
                "residual_history": stage["residual_history"],
                "newton_history": [
                    {
                        key: value
                        for key, value in row.items()
                        if key
                        not in {
                            "wall_clock_seconds",
                            "elapsed_seconds",
                            "preconditioner_setup_seconds",
                        }
                    }
                    for row in stage.get("newton_history", [])
                ],
            }
        )
    return rows


def _component_residual_norms(
    residual: torch.Tensor,
    *,
    nr: int,
    nw: int,
) -> dict[str, float]:
    blocks = residual[:-1].reshape(5, nr, nw)
    names = ("psi_R", "psi_I", "A_0", "A_r", "A_w")
    result = {
        name: float(torch.max(torch.abs(blocks[idx])).detach().cpu().item())
        for idx, name in enumerate(names)
    }
    result["mass_constraint"] = float(torch.abs(residual[-1]).detach().cpu().item())
    return result


def _make_bundle(
    *,
    config: HarnessConfig,
    grid: TensorProductGrid,
    state: torch.Tensor,
    summary: dict[str, Any],
) -> dict[str, Any]:
    cfg = config.branch
    fields, chemical_potential = unpack_coupled_fields(
        state,
        grid,
        has_chemical_potential=True,
    )
    if chemical_potential is None:
        raise RuntimeError("coupled branch state did not include chemical potential")

    boundaries = branch_boundary_conditions(cfg)
    final_eos_K = cfg.continuation_K_values[len(summary["stages"]) - 1]
    residual = coupled_branch_residual(
        state,
        grid,
        cfg,
        eos_K=final_eos_K,
        boundaries=boundaries,
    )
    residual_linf = float(torch.max(torch.abs(residual)).detach().cpu().item())
    residual_l2 = float(torch.linalg.vector_norm(residual).detach().cpu().item())
    density = fields.psi_real**2 + fields.psi_imag**2
    r0 = geometry_radius_torch(grid.w_centers, cfg)
    z = localization_weight_torch(grid.w_centers, cfg)

    payload: dict[str, Any] = {
        "schema": "stage1_m1c_prep_wp1_background/v1",
        "scope": "M1c-prep target-blind machinery background; placeholder parameters; no Gate A; not physical",
        "source_revision": source_revision(),
        "solver": {
            "engine": "torch",
            "entry_point": "stage1_solver.coupled_branch.run_branch_continuation",
            "residual_entry_point": "stage1_solver.coupled_branch.coupled_branch_residual",
            "validated_step": "Step-3 coupled matter-Maxwell branch continuation",
            "grid_name": summary["grid"],
            "converged": bool(summary["converged"]),
            "message": summary["message"],
            "dof": int(summary["dof"]),
            "placeholder_label": cfg.placeholder_label,
            "final_eos_K": float(final_eos_K),
            "chemical_potential": float(chemical_potential.detach().cpu().item()),
        },
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
            "r_faces": _numpy_list(grid.r_faces),
            "r_centers": _numpy_list(grid.r_centers),
            "w_faces": _numpy_list(grid.w_faces),
            "w_centers": _numpy_list(grid.w_centers),
            "radial_face_areas": _numpy_list(grid.radial_face_areas),
            "radial_shell_volumes": _numpy_list(grid.radial_shell_volumes),
            "cell_volumes": _numpy_list(grid.cell_volumes),
        },
        "fields": {
            "psi_R0": _numpy_list(fields.psi_real),
            "psi_I0": _numpy_list(fields.psi_imag),
            "A_00": _numpy_list(fields.a0),
            "A_r0": _numpy_list(fields.ar),
            "A_w0": _numpy_list(fields.aw),
        },
        "derived": {
            "rho0": _numpy_list(density),
            "R0_w": _numpy_list(r0),
            "Z_w": _numpy_list(z),
            "mass": float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
            "A_00_linf": float(torch.max(torch.abs(fields.a0)).detach().cpu().item()),
            "A_r0_linf": float(torch.max(torch.abs(fields.ar)).detach().cpu().item()),
            "A_w0_linf": float(torch.max(torch.abs(fields.aw)).detach().cpu().item()),
        },
        "residuals": {
            "coupled_stationary_linf": residual_linf,
            "coupled_stationary_l2_unweighted": residual_l2,
            "component_linf": _component_residual_norms(
                residual,
                nr=grid.spec.nr,
                nw=grid.spec.nw,
            ),
            "newton_floor_reference": "~1e-8",
            "smoke_background_residual_reference": "~243 for the old analytic M1b smoke rho0",
            "self_consistent": bool(residual_linf < DEFAULT_RESIDUAL_CEILING),
        },
        "convergence_evidence": {
            "continuation_stages": _stage_evidence(summary),
            "final_mass": summary["final_mass"],
            "surrogate_values": summary["surrogate_values"],
            "boundaries": summary["boundaries"],
        },
        "config": {
            "backend": asdict(config.backend),
            "branch": asdict(cfg),
        },
        "byte_reproducibility": {
            "hash_algorithm": "sha256",
            "canonical_json": "sort_keys=True,separators=(',',':')",
            "excluded_live_fields": [
                "wall_clock_seconds",
                "peak_memory_mb",
                "manifest_path",
                "preconditioner_setup_seconds",
                "timestamps",
            ],
            "content_hash_covers": "all fields except this content_hash value",
        },
    }
    payload["content_hash"] = stable_hash(payload)
    return payload


def export_m1c_prep_background(
    *,
    run_root: str = DEFAULT_RUN_ROOT,
    grid_level: tuple[int, int] | None = None,
    residual_ceiling: float = DEFAULT_RESIDUAL_CEILING,
) -> tuple[Path, dict[str, Any]]:
    config = HarnessConfig(
        run_root=run_root,
        report_path=DEFAULT_REPORT_PATH,
    )
    dtype = configure_backend(config.backend)
    cfg = config.branch
    level = grid_level or cfg.solve_grid
    grid = TensorProductGrid.create(
        TensorGridSpec(
            r_max=cfg.r_max,
            nr=int(level[0]),
            w_min=cfg.w_min,
            w_max=cfg.w_max,
            nw=int(level[1]),
        ),
        dtype=dtype,
        device=config.backend.device,
    )
    state, summary = run_branch_continuation(
        config,
        grid,
        grid_name=f"m1c_prep_nr_{grid.spec.nr}_nw_{grid.spec.nw}",
    )
    bundle = _make_bundle(config=config, grid=grid, state=state, summary=summary)
    residual = float(bundle["residuals"]["coupled_stationary_linf"])
    if not summary["converged"] or residual >= residual_ceiling:
        raise RuntimeError(
            "M1c-prep background solve did not reach the self-consistent residual floor: "
            f"converged={summary['converged']} residual={residual:.6e} ceiling={residual_ceiling:.6e}"
        )

    out_dir = Path(run_root) / "background_bundles"
    out_dir.mkdir(parents=True, exist_ok=True)
    content_hash = bundle["content_hash"]
    out_path = out_dir / f"m1c_prep_wp1_background_{content_hash}.json"
    latest_path = out_dir / "m1c_prep_wp1_background_latest.json"
    text = json.dumps(bundle, indent=2, sort_keys=True)
    out_path.write_text(text + "\n", encoding="utf-8")
    latest_path.write_text(text + "\n", encoding="utf-8")

    return out_path, bundle


def main() -> int:
    parser = argparse.ArgumentParser(description="Export target-blind M1c-prep WP1 background JSON")
    parser.add_argument("--run-root", default=DEFAULT_RUN_ROOT)
    parser.add_argument("--grid", default=None, help="Override solve grid as NR,NW; default uses BranchSmokeConfig.solve_grid")
    parser.add_argument("--residual-ceiling", type=float, default=DEFAULT_RESIDUAL_CEILING)
    args = parser.parse_args()

    grid_level = None
    if args.grid:
        parts = [int(part) for part in args.grid.split(",")]
        if len(parts) != 2:
            raise SystemExit("--grid must be NR,NW")
        grid_level = (parts[0], parts[1])

    path, bundle = export_m1c_prep_background(
        run_root=args.run_root,
        grid_level=grid_level,
        residual_ceiling=args.residual_ceiling,
    )
    print(f"background_path: {path}")
    print(f"content_hash: {bundle['content_hash']}")
    print(f"grid: {bundle['grid']['nr']}x{bundle['grid']['nw']}")
    print(f"coupled_stationary_linf: {bundle['residuals']['coupled_stationary_linf']:.12e}")
    print(f"A_00_linf: {bundle['derived']['A_00_linf']:.12e}")
    print(f"A_r0_linf: {bundle['derived']['A_r0_linf']:.12e}")
    print(f"A_w0_linf: {bundle['derived']['A_w0_linf']:.12e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
