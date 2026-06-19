"""Path-A B2a background export, BdG cross-check, and validation.

This module stays inside the B2a scope: it exports the closed Path-A
background, derives matter-sector BdG modes/couplings, and validates the
independent Mathematica derivation. It deliberately does not assemble Maxwell
transfer terms or downstream response residuals.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.linalg import eig
import torch

from .backend import configure_backend
from .config import BackendConfig, PreconditionerConfig, TensorGridSpec, source_revision
from .coupled_branch import (
    ClosedCoupledFields,
    _create_branch_grid,
    branch_boundary_conditions,
    confinement_potential_torch,
    confinement_wall_to_matter_coefficient_torch,
    initial_closed_branch_state,
    localization_weight_torch,
    patha_closed_branch_residual,
    resample_closed_branch_state,
    unpack_closed_coupled_fields,
)
from .newton import solve_newton_jvp
from .patha_closed_newton import default_closed_branch_config
from .patha_extraction import bdg_moments, solve_l2_wall_eigenproblem
from .patha_static_balance import SSigmaSpec, resolve_s_sigma
from .preconditioners import make_closed_coupled_colored_sparse_jacobian_lu_factory


DEFAULT_RUN_ROOT = Path("software/stage1_solver/runs/patha_b2a_bdg_derivation")
DEFAULT_REPORT_PATH = Path("software/stage1_solver/reports/patha_b2a_bdg_derivation_report.md")
SMOKE_RESIDUAL_REFERENCE = 243.39250922131095
FROZEN_A = 1.0
FROZEN_L = 37.0 / 20.0
DEFAULT_BACKGROUND_GRID = (16, 16)
DEFAULT_BDG_MODES = 30
DEFAULT_PROFILE_POINTS = 129
MOMENT_KEYS = ("B0", "B2", "B4")
MODAL_SWEEP_COUNTS = (3, 5, 8, 15, 30)
MODAL_TRUNCATION_TOL = 1.0e-4
DUAL_ENGINE_TOLERANCES = {
    "varpi": {"abs": 1.0e-9, "rel": 1.0e-9},
    "coupling": {"abs": 1.0e-12, "rel": 1.0e-8},
    "moments": {"abs": 1.0e-16, "rel": 1.0e-9},
}
CONSUMER_MOMENT_TOLERANCES = {"abs": 1.0e-16, "rel": 1.0e-9}


def _full_stable_hash(obj: Any) -> str:
    text = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _numpy_list(value: torch.Tensor | np.ndarray | Sequence[float]) -> Any:
    if isinstance(value, torch.Tensor):
        value = value.detach().cpu().numpy()
    return np.asarray(value, dtype=np.float64).tolist()


def _format_tau(tau: float) -> str:
    text = f"{float(tau):.12g}".replace("-", "m").replace(".", "p")
    return text


def _parse_grid(text: str) -> tuple[int, int]:
    parts = [int(part) for part in text.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("grid must be NR,NW")
    return parts[0], parts[1]


def frozen_patha_b2a_branch(*, grid: tuple[int, int], tau: float):
    """Return the closed-Newton branch with frozen Path-A geometry.

    The wall family is controlled by ``tau`` through ``SSigmaSpec``. The
    remaining matter/trap constants are the existing closed-Newton solver
    constants and are recorded in the exported bundle.
    """

    base = default_closed_branch_config()
    preconditioner = replace(
        base.newton.preconditioner,
        type="colored_sparse_jacobian_lu",
        side="left",
        rebuild_policy="every_newton_step",
        stencil_radius=3,
        color_separation=7,
        factorization="splu",
        diagonal_shift=1.0e-10,
        drop_tolerance=0.0,
        fill_factor=10.0,
        permutation="COLAMD",
    )
    newton = replace(
        base.newton,
        residual_atol=2.0e-9,
        residual_rtol=2.0e-9,
        step_atol=1.0e-12,
        step_rtol=1.0e-11,
        max_newton_iters=18,
        gmres_rtol=2.0e-9,
        gmres_atol=1.0e-11,
        gmres_restart=192,
        gmres_maxiter=12,
        max_line_search_iters=18,
        accept_best_line_search_decrease=True,
        finite_difference_jvp_epsilon=1.0e-5,
        preconditioner=preconditioner,
    )
    return replace(
        base,
        name=f"patha_b2a_closed_frozen_hooke_tau_{_format_tau(tau)}",
        placeholder_label=(
            "Path-A B2a closed solve on frozen homogeneous_isotropic_hooke_v1; "
            "matter/trap constants inherited from closed Newton solver and recorded"
        ),
        r_max=1.6,
        w_min=0.0,
        w_max=FROZEN_L,
        solve_grid=grid,
        ladder_levels=(grid,),
        geometry_profile="cubic_smoothstep",
        r_mouth=FROZEN_A,
        r_exit=FROZEN_A,
        newton=newton,
    )


def frozen_s_sigma_spec(tau: float) -> SSigmaSpec:
    return SSigmaSpec.homogeneous_isotropic_hooke(
        tau=float(tau),
        a=FROZEN_A,
        w_min=0.0,
        w_max=FROZEN_L,
    )


def _closed_component_residual_norms(
    residual: torch.Tensor,
    *,
    nr: int,
    nw: int,
) -> dict[str, float]:
    n = nr * nw
    blocks = residual[: 5 * n].reshape(5, nr, nw)
    names = ("psi_R", "psi_I", "A_0", "A_r", "A_w")
    result = {
        name: float(torch.max(torch.abs(blocks[idx])).detach().cpu().item())
        for idx, name in enumerate(names)
    }
    wall = residual[5 * n : 5 * n + nw]
    result["wall_R0"] = float(torch.max(torch.abs(wall)).detach().cpu().item())
    result["mass_constraint"] = float(torch.abs(residual[-1]).detach().cpu().item())
    return result


def _field_linf(fields: ClosedCoupledFields) -> dict[str, float]:
    return {
        "psi_real": float(torch.max(torch.abs(fields.psi_real)).detach().cpu().item()),
        "psi_imag": float(torch.max(torch.abs(fields.psi_imag)).detach().cpu().item()),
        "a0": float(torch.max(torch.abs(fields.a0)).detach().cpu().item()),
        "ar": float(torch.max(torch.abs(fields.ar)).detach().cpu().item()),
        "aw": float(torch.max(torch.abs(fields.aw)).detach().cpu().item()),
        "r0": float(torch.max(torch.abs(fields.r0)).detach().cpu().item()),
    }


def solve_closed_background(
    *,
    tau: float,
    grid_level: tuple[int, int],
    initial_state: torch.Tensor | None = None,
    initial_grid: Any | None = None,
    continuation_K_values: Sequence[float] | None = None,
    initialization_note: str | None = None,
    newton_overrides: Mapping[str, Any] | None = None,
) -> tuple[Any, Any, Any, dict[str, Any]]:
    dtype = configure_backend(BackendConfig())
    branch = frozen_patha_b2a_branch(grid=grid_level, tau=tau)
    if newton_overrides:
        branch = replace(branch, newton=replace(branch.newton, **dict(newton_overrides)))
    spec = frozen_s_sigma_spec(tau)
    provider = resolve_s_sigma(spec)
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    boundaries = branch_boundary_conditions(branch)
    if initial_state is None:
        state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
        initialization_source = "default_closed_branch_state"
    else:
        state = initial_state.detach().clone().to(dtype=dtype, device="cpu")
        if initial_grid is not None and (
            int(initial_grid.spec.nr) != int(grid.spec.nr)
            or int(initial_grid.spec.nw) != int(grid.spec.nw)
            or abs(float(initial_grid.spec.r_max) - float(grid.spec.r_max)) > 0.0
            or abs(float(initial_grid.spec.w_max) - float(grid.spec.w_max)) > 0.0
        ):
            state = resample_closed_branch_state(state, initial_grid, grid, branch)
        initialization_source = "warm_start_state"
    shared_preconditioner_factory = make_closed_coupled_colored_sparse_jacobian_lu_factory(grid)
    continuation_values = tuple(
        float(value) for value in (
            branch.continuation_K_values if continuation_K_values is None else continuation_K_values
        )
    )
    if not continuation_values:
        raise ValueError("continuation_K_values must contain at least one EOS value")

    stages: list[dict[str, Any]] = []
    converged = True
    message = "closed B2a continuation completed"
    started = time.perf_counter()
    previous_grid = grid

    for eos_K in continuation_values:
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
            preconditioner_factory=shared_preconditioner_factory,
        )
        state = result.x.detach()
        fields_after, _ = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
        stages.append(
            {
                "eos_K": float(eos_K),
                "converged": bool(result.converged),
                "iterations": int(result.iterations),
                "initial_residual_norm": float(result.initial_residual_norm),
                "final_residual_norm": float(result.final_residual_norm),
                "tolerance": float(result.tolerance),
                "message": result.message,
                "gmres_iterations": [
                    row.gmres_iterations
                    for row in result.history
                    if row.gmres_iterations is not None
                ],
                "residual_history": [float(row.residual_norm) for row in result.history],
                "r0_min": float(torch.min(fields_after.r0).detach().cpu().item()),
                "r0_max": float(torch.max(fields_after.r0).detach().cpu().item()),
            }
        )
        if not result.converged:
            converged = False
            message = f"continuation stopped at eos_K={eos_K}: {result.message}"
            break
        previous_grid = grid

    final_eos_K = float(stages[-1]["eos_K"])
    final_residual_fn = lambda x: patha_closed_branch_residual(
        x,
        grid,
        branch,
        eos_K=final_eos_K,
        boundaries=boundaries,
        s_sigma=provider,
    )
    residual = final_residual_fn(state).detach()
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    density = fields.psi_real**2 + fields.psi_imag**2
    residual_linf = float(torch.max(torch.abs(residual)).detach().cpu().item())
    residual_l2 = float(torch.linalg.vector_norm(residual).detach().cpu().item())
    r0_min = float(torch.min(fields.r0).detach().cpu().item())
    if r0_min <= 0.0:
        converged = False
        message = "closed solve produced nonpositive R0"
    summary = {
        "branch": branch,
        "s_sigma_spec": spec,
        "source_revision": source_revision(),
        "converged": bool(converged),
        "message": message,
        "stages": stages,
        "final_eos_K": final_eos_K,
        "wall_clock_seconds": time.perf_counter() - started,
        "final_residual_linf": residual_linf,
        "final_residual_l2": residual_l2,
        "component_linf": _closed_component_residual_norms(
            residual,
            nr=grid.spec.nr,
            nw=grid.spec.nw,
        ),
        "final_mass": float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
        "chemical_potential": float(mu.detach().cpu().item()) if mu is not None else None,
        "r0_min": r0_min,
        "r0_max": float(torch.max(fields.r0).detach().cpu().item()),
        "field_linf": _field_linf(fields),
        "previous_grid": previous_grid.to_dict(),
        "initialization": {
            "source": initialization_source,
            "note": initialization_note,
            "continuation_K_values": list(continuation_values),
            "newton_overrides": {} if newton_overrides is None else dict(newton_overrides),
        },
    }
    return grid, state, fields, summary


def make_background_bundle(
    *,
    tau: float,
    grid_level: tuple[int, int],
    initial_state: torch.Tensor | None = None,
    initial_grid: Any | None = None,
    continuation_K_values: Sequence[float] | None = None,
    initialization_note: str | None = None,
    newton_overrides: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    grid, _state, fields, summary = solve_closed_background(
        tau=tau,
        grid_level=grid_level,
        initial_state=initial_state,
        initial_grid=initial_grid,
        continuation_K_values=continuation_K_values,
        initialization_note=initialization_note,
        newton_overrides=newton_overrides,
    )
    branch = summary["branch"]
    spec = summary["s_sigma_spec"]
    density = fields.psi_real**2 + fields.psi_imag**2
    z = localization_weight_torch(grid.w_centers, branch)
    potential = confinement_potential_torch(grid, branch, radius=fields.r0)
    drive = confinement_wall_to_matter_coefficient_torch(grid, branch, radius=fields.r0)

    payload: dict[str, Any] = {
        "schema": "stage1_patha_b2a_closed_background/v1",
        "scope": (
            "Path-A B2a machine handoff: closed self-consistent background on frozen "
            "homogeneous_isotropic_hooke_v1; target-blind; no Maxwell transfer or response residual"
        ),
        "source_revision": summary["source_revision"],
        "solver": {
            "engine": "torch",
            "entry_point": "stage1_solver.patha_b2a_bdg.solve_closed_background",
            "residual_entry_point": "stage1_solver.coupled_branch.patha_closed_branch_residual",
            "closed_newton_source": "stage1_solver.patha_closed_newton + stage1_solver.coupled_branch",
            "converged": summary["converged"],
            "message": summary["message"],
            "final_eos_K": summary["final_eos_K"],
            "chemical_potential": summary["chemical_potential"],
            "wall_clock_seconds": summary["wall_clock_seconds"],
            "initialization": summary["initialization"],
        },
        "geometry": {
            "a": FROZEN_A,
            "L": FROZEN_L,
            "w_min": 0.0,
            "w_max": FROZEN_L,
            "r_mouth": float(branch.r_mouth),
            "initial_r_exit": float(branch.r_exit),
            "radial_domain_r_max": float(branch.r_max),
        },
        "constants": {
            "eos_K": float(summary["final_eos_K"]),
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
        "s_sigma_spec": spec.to_dict(),
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
            "R0_w": _numpy_list(fields.r0),
        },
        "derived": {
            "rho0": _numpy_list(density),
            "R0_w": _numpy_list(fields.r0),
            "Z_w": _numpy_list(z),
            "V_conf": _numpy_list(potential),
            "wall_drive_k1": _numpy_list(drive),
            "mass": summary["final_mass"],
            "field_linf": summary["field_linf"],
        },
        "residuals": {
            "closed_stationary_linf": summary["final_residual_linf"],
            "closed_stationary_l2_unweighted": summary["final_residual_l2"],
            "component_linf": summary["component_linf"],
            "smoke_background_residual_reference": SMOKE_RESIDUAL_REFERENCE,
            "self_consistent": bool(
                summary["converged"]
                and summary["final_residual_linf"] < 1.0e-6
                and summary["final_residual_linf"] < SMOKE_RESIDUAL_REFERENCE / 1.0e6
            ),
        },
        "convergence_evidence": {
            "continuation_stages": summary["stages"],
            "final_mass": summary["final_mass"],
            "r0_min": summary["r0_min"],
            "r0_max": summary["r0_max"],
            "branch_config": asdict(branch),
        },
        "byte_reproducibility": {
            "hash_algorithm": "sha256",
            "canonical_json": "sort_keys=True,separators=(',',':')",
            "content_hash_covers": "all fields except this content_hash value",
        },
    }
    payload["content_hash"] = _full_stable_hash(payload)
    return payload


def export_background(
    *,
    tau: float,
    grid_level: tuple[int, int],
    run_root: Path,
) -> tuple[Path, dict[str, Any]]:
    bundle = make_background_bundle(tau=tau, grid_level=grid_level)
    out_dir = run_root / "backgrounds"
    out_dir.mkdir(parents=True, exist_ok=True)
    tau_label = _format_tau(tau)
    out_path = out_dir / f"patha_b2a_closed_background_tau_{tau_label}_{bundle['content_hash']}.json"
    latest_path = out_dir / f"patha_b2a_closed_background_tau_{tau_label}_latest.json"
    text = json.dumps(bundle, indent=2, sort_keys=True, default=_json_default)
    out_path.write_text(text + "\n", encoding="utf-8")
    latest_path.write_text(text + "\n", encoding="utf-8")
    return out_path, bundle


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _interp1_clamped(x_old: np.ndarray, y_old: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    return np.interp(x_new, x_old, y_old, left=y_old[0], right=y_old[-1])


def _interp2_clamped(
    r_old: np.ndarray,
    w_old: np.ndarray,
    values: np.ndarray,
    r_new: np.ndarray,
    w_new: np.ndarray,
) -> np.ndarray:
    rows = np.empty((r_old.size, w_new.size), dtype=np.result_type(values, np.float64))
    for i in range(r_old.size):
        rows[i, :] = np.interp(w_new, w_old, values[i, :], left=values[i, 0], right=values[i, -1])
    out = np.empty((r_new.size, w_new.size), dtype=np.result_type(values, np.float64))
    for j in range(w_new.size):
        out[:, j] = np.interp(r_new, r_old, rows[:, j], left=rows[0, j], right=rows[-1, j])
    return out


def _flatten_w_major(values: np.ndarray) -> np.ndarray:
    return np.asarray(values).T.reshape(-1)


def _unflatten_w_major(values: np.ndarray, *, nr: int, nw: int) -> np.ndarray:
    return np.asarray(values).reshape(nw, nr).T


def _trapz(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.trapz(y, x))


def _weighted_norm(x: np.ndarray, weight: np.ndarray, f: np.ndarray) -> float:
    return _trapz(x, weight * f * f)


def _make_bdg_grid(background: Mapping[str, Any], *, nr: int, nw: int) -> dict[str, np.ndarray | float | int]:
    r_max = float(background["grid"]["r_max"])
    L = float(background["geometry"]["L"])
    r_faces = np.linspace(0.0, r_max, nr + 1, dtype=np.float64)
    w_faces = np.linspace(0.0, L, nw + 1, dtype=np.float64)
    r_centers = 0.5 * (r_faces[:-1] + r_faces[1:])
    w_centers = 0.5 * (w_faces[:-1] + w_faces[1:])
    radial_face_areas = 4.0 * math.pi * r_faces**2
    radial_shell_volumes = (4.0 * math.pi / 3.0) * (r_faces[1:] ** 3 - r_faces[:-1] ** 3)
    return {
        "nr": nr,
        "nw": nw,
        "r_faces": r_faces,
        "w_faces": w_faces,
        "r_centers": r_centers,
        "w_centers": w_centers,
        "radial_face_areas": radial_face_areas,
        "radial_shell_volumes": radial_shell_volumes,
        "dr": float(r_faces[1] - r_faces[0]),
        "dw": float(w_faces[1] - w_faces[0]),
    }


def _cell_index(i: int, j: int, nr: int) -> int:
    return j * nr + i


def assemble_laplacian_matrix(grid: Mapping[str, Any]) -> np.ndarray:
    nr = int(grid["nr"])
    nw = int(grid["nw"])
    dr = float(grid["dr"])
    dw = float(grid["dw"])
    radial_volumes = np.asarray(grid["radial_shell_volumes"], dtype=np.float64)
    radial_areas = np.asarray(grid["radial_face_areas"], dtype=np.float64)
    n_cells = nr * nw
    lap = np.zeros((n_cells, n_cells), dtype=np.float64)
    for j in range(nw):
        for i in range(nr):
            row = _cell_index(i, j, nr)
            vol = radial_volumes[i]
            if i > 0:
                coeff = radial_areas[i] / (dr * vol)
                lap[row, _cell_index(i - 1, j, nr)] += coeff
                lap[row, row] -= coeff
            if i < nr - 1:
                coeff = radial_areas[i + 1] / (dr * vol)
                lap[row, _cell_index(i + 1, j, nr)] += coeff
                lap[row, row] -= coeff
            else:
                coeff = radial_areas[i + 1] / (dr * vol)
                lap[row, row] -= 2.0 * coeff

            if j == 0:
                lap[row, _cell_index(i, j + 1, nr)] += 1.0 / dw**2
                lap[row, row] -= 1.0 / dw**2
            elif j == nw - 1:
                lap[row, _cell_index(i, j - 1, nr)] += 1.0 / dw**2
                lap[row, row] -= 1.0 / dw**2
            else:
                lap[row, _cell_index(i, j - 1, nr)] += 1.0 / dw**2
                lap[row, _cell_index(i, j + 1, nr)] += 1.0 / dw**2
                lap[row, row] -= 2.0 / dw**2
    return lap


def _background_arrays_on_bdg_grid(background: Mapping[str, Any], grid: Mapping[str, Any]) -> dict[str, np.ndarray]:
    r_old = np.asarray(background["grid"]["r_centers"], dtype=np.float64)
    w_old = np.asarray(background["grid"]["w_centers"], dtype=np.float64)
    r_new = np.asarray(grid["r_centers"], dtype=np.float64)
    w_new = np.asarray(grid["w_centers"], dtype=np.float64)

    def field2(name: str) -> np.ndarray:
        return _interp2_clamped(
            r_old,
            w_old,
            np.asarray(background["fields"][name], dtype=np.float64),
            r_new,
            w_new,
        )

    r0 = _interp1_clamped(
        w_old,
        np.asarray(background["fields"]["R0_w"], dtype=np.float64),
        w_new,
    )
    return {
        "psi_real": field2("psi_R0"),
        "psi_imag": field2("psi_I0"),
        "a0": field2("A_00"),
        "ar": field2("A_r0"),
        "aw": field2("A_w0"),
        "r0": r0,
    }


def assemble_bdg_operator_python(
    background: Mapping[str, Any],
    *,
    nr: int,
    nw: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    grid = _make_bdg_grid(background, nr=nr, nw=nw)
    arrays = _background_arrays_on_bdg_grid(background, grid)
    constants = background["constants"]
    hbar = float(constants["hbar"])
    particle_mass = float(constants["particle_mass"])
    gauge_charge = float(constants["gauge_charge"])
    eos_K = float(constants["eos_K"])
    mu = float(background["solver"]["chemical_potential"])
    radial_wall_strength = float(constants["radial_wall_strength"])
    axial_trap_strength = float(constants["axial_trap_strength"])
    L = float(background["geometry"]["L"])
    ell_factor = 6.0

    r_centers = np.asarray(grid["r_centers"], dtype=np.float64)
    w_centers = np.asarray(grid["w_centers"], dtype=np.float64)
    r_2d = r_centers[:, None]
    w_2d = w_centers[None, :]
    r0_2d = arrays["r0"][None, :]
    psi = arrays["psi_real"] + 1j * arrays["psi_imag"]
    rho = np.abs(psi) ** 2
    a0 = arrays["a0"]

    v_conf = radial_wall_strength * (r_2d / r0_2d) ** 4 + 0.5 * axial_trap_strength * (
        (w_2d - 0.5 * L) / L
    ) ** 2
    drive = 4.0 * radial_wall_strength * r_2d**4 / r0_2d**5

    laplacian = assemble_laplacian_matrix(grid)
    bdg_radii = np.tile(r_centers, nw)
    angular = np.diag(ell_factor / bdg_radii**2)
    kinetic = -(hbar**2 / (2.0 * particle_mass)) * (laplacian - angular)

    rho_vec = _flatten_w_major(rho)
    psi_vec = _flatten_w_major(psi)
    a0_vec = _flatten_w_major(a0)
    v_conf_vec = _flatten_w_major(v_conf)
    enthalpy = (5.0 * eos_K / 4.0) * rho_vec**4
    dh_drho = 5.0 * eos_K * rho_vec**3
    rho_dh = rho_vec * dh_drho
    single_diag = v_conf_vec + gauge_charge * a0_vec - mu + enthalpy
    diag_block = kinetic.astype(np.complex128) + np.diag(single_diag + rho_dh)
    pair_diag = dh_drho * psi_vec**2
    pair_block = np.diag(pair_diag)
    bdg = np.block(
        [
            [diag_block, pair_block],
            [-np.conjugate(pair_block), -np.conjugate(diag_block)],
        ]
    )
    diagnostics = {
        "grid": grid,
        "arrays": arrays,
        "psi_vec": psi_vec,
        "rho_vec": rho_vec,
        "drive_vec": _flatten_w_major(drive),
        "v_conf_vec": v_conf_vec,
        "cell_volume_weights": np.tile(np.asarray(grid["radial_shell_volumes"], dtype=np.float64), nw)
        * float(grid["dw"]),
        "spatial_gauge_linf": float(max(np.max(np.abs(arrays["ar"])), np.max(np.abs(arrays["aw"])))),
    }
    return bdg, diagnostics


def make_wall_input(
    background: Mapping[str, Any],
    *,
    profile_points_count: int,
    out_dir: Path | None = None,
) -> tuple[dict[str, Any], Path | None]:
    L = float(background["geometry"]["L"])
    w_old = np.asarray(background["grid"]["w_centers"], dtype=np.float64)
    r0_old = np.asarray(background["fields"]["R0_w"], dtype=np.float64)
    profile_points = np.linspace(0.0, L, profile_points_count, dtype=np.float64)
    r0_profile = _interp1_clamped(w_old, r0_old, profile_points)
    spec = SSigmaSpec.from_dict(background["s_sigma_spec"])
    wall = solve_l2_wall_eigenproblem(profile_points, r0_profile, spec)
    weight = np.ones_like(profile_points)
    wall_payload: dict[str, Any] = {
        "schema": "stage1_patha_b2a_b1_wall_input/v1",
        "background_content_hash": background["content_hash"],
        "source": "stage1_solver.patha_extraction.solve_l2_wall_eigenproblem",
        "profile_points": profile_points.tolist(),
        "R0_profile": r0_profile.tolist(),
        "weight": weight.tolist(),
        "chi": wall.chi.tolist(),
        "K": float(wall.K),
        "M": float(wall.M),
        "normalization_chi_T_W_chi": float(wall.normalization),
        "trapz_weight_norm": _weighted_norm(profile_points, weight, wall.chi),
        "orientation_integral": float(wall.orientation_integral),
        "mode_index": int(wall.mode_index),
    }
    wall_payload["content_hash"] = _full_stable_hash(wall_payload)
    out_path = None
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        tau_label = _format_tau(float(background["constants"]["tau"]))
        out_path = out_dir / f"patha_b2a_wall_input_tau_{tau_label}.json"
        out_path.write_text(
            json.dumps(wall_payload, indent=2, sort_keys=True, default=_json_default) + "\n",
            encoding="utf-8",
        )
    return wall_payload, out_path


def _rotate_by_max_component(vec: np.ndarray) -> np.ndarray:
    idx = int(np.argmax(np.abs(vec)))
    if abs(vec[idx]) == 0.0:
        return vec
    return np.exp(-1j * np.angle(vec[idx])) * vec


def _stieltjes_gap(moments: Mapping[str, float]) -> float:
    return float(moments["B0"] * moments["B4"] - moments["B2"] ** 2)


def _moments_with_stieltjes(modes: Sequence[Mapping[str, Any]]) -> dict[str, float]:
    moments = {key: float(value) for key, value in bdg_moments(modes).items()}
    moments["B0_B4_minus_B2_squared"] = _stieltjes_gap(moments)
    return moments


def _moment_rel_errors(current: Mapping[str, float], reference: Mapping[str, float]) -> dict[str, float]:
    return {
        key: _relative_change(float(current[key]), float(reference[key]))
        for key in MOMENT_KEYS
    }


def _build_modal_truncation_metadata(
    modal_modes: Sequence[Mapping[str, Any]],
    *,
    modes_to_export: int,
    tolerance: float = MODAL_TRUNCATION_TOL,
) -> dict[str, Any]:
    if not modal_modes:
        raise RuntimeError("modal truncation sweep requires at least one positive mode")
    if modes_to_export > len(modal_modes):
        raise RuntimeError(
            f"cannot export {modes_to_export} modes from only {len(modal_modes)} positive modes"
        )

    all_count = len(modal_modes)
    sweep_counts = [count for count in MODAL_SWEEP_COUNTS if count < all_count]
    if all_count not in sweep_counts:
        sweep_counts.append(all_count)

    all_moments = _moments_with_stieltjes(modal_modes)
    rows: list[dict[str, Any]] = []
    previous_moments: dict[str, float] | None = None
    for count in sweep_counts:
        moments = _moments_with_stieltjes(modal_modes[:count])
        rel_previous = (
            None
            if previous_moments is None
            else _moment_rel_errors(moments, previous_moments)
        )
        rel_all = _moment_rel_errors(moments, all_moments)
        rows.append(
            {
                "M": int(count),
                "label": "all_positive" if count == all_count else str(count),
                "B0": moments["B0"],
                "B2": moments["B2"],
                "B4": moments["B4"],
                "B0_B4_minus_B2_squared": moments["B0_B4_minus_B2_squared"],
                "rel_change_from_previous": rel_previous,
                "max_rel_change_from_previous": (
                    None if rel_previous is None else max(rel_previous.values())
                ),
                "rel_error_vs_all_positive": rel_all,
                "max_rel_error_vs_all_positive": max(rel_all.values()),
            }
        )
        previous_moments = moments

    selected_moments = _moments_with_stieltjes(modal_modes[:modes_to_export])
    selected_error = _moment_rel_errors(selected_moments, all_moments)
    eligible = [
        row["M"]
        for row in rows
        if row["label"] != "all_positive"
        and row["max_rel_error_vs_all_positive"] <= tolerance
    ]
    smallest_count_meeting_tolerance = min(eligible) if eligible else all_count

    return {
        "sweep_counts": sweep_counts,
        "sweep": rows,
        "all_positive_mode_count": int(all_count),
        "all_positive_moments": all_moments,
        "exported_mode_count": int(modes_to_export),
        "exported_moments": selected_moments,
        "truncation_error_at_export": selected_error,
        "max_truncation_error_at_export": max(selected_error.values()),
        "tolerance": float(tolerance),
        "smallest_count_meeting_tolerance": int(smallest_count_meeting_tolerance),
        "selection_rule": (
            "exported K must have max_n |B_n(K)-B_n(all)|/|B_n(all)| <= tolerance; "
            "all-positive is the reference tail check"
        ),
    }


def solve_bdg_python(
    background: Mapping[str, Any],
    wall_input: Mapping[str, Any],
    *,
    nr: int,
    nw: int,
    modes_to_export: int = DEFAULT_BDG_MODES,
) -> dict[str, Any]:
    bdg, diag = assemble_bdg_operator_python(background, nr=nr, nw=nw)
    evals, evecs = eig(bdg)
    order = np.argsort(evals.real)
    pairs = [(evals[idx], evecs[:, idx]) for idx in order]
    positive = [
        (val, vec)
        for val, vec in pairs
        if val.real > 1.0e-8 and abs(val.imag) < 1.0e-7
    ]
    if len(positive) < modes_to_export:
        raise RuntimeError(f"found only {len(positive)} positive real BdG modes")

    grid = diag["grid"]
    nr = int(grid["nr"])
    nw = int(grid["nw"])
    n_cells = nr * nw
    cell_weights = np.asarray(diag["cell_volume_weights"], dtype=np.float64)
    psi_vec = np.asarray(diag["psi_vec"], dtype=np.complex128)
    drive_vec = np.asarray(diag["drive_vec"], dtype=np.float64)
    radial_volumes = np.asarray(grid["radial_shell_volumes"], dtype=np.float64)
    w_centers = np.asarray(grid["w_centers"], dtype=np.float64)
    points = np.asarray(wall_input["profile_points"], dtype=np.float64)
    chi = np.asarray(wall_input["chi"], dtype=np.float64)
    weight = np.asarray(wall_input["weight"], dtype=np.float64)

    modes: list[dict[str, Any]] = []
    modal_spectrum: list[dict[str, Any]] = []
    residuals: list[float] = []
    for mode_idx, (eigenvalue, raw_vec) in enumerate(positive, start=1):
        rotated = _rotate_by_max_component(raw_vec)
        u0 = rotated[:n_cells]
        v0 = rotated[n_cells:]
        symp_before = float(np.real(np.sum(cell_weights * (np.abs(u0) ** 2 - np.abs(v0) ** 2))))
        scale = 1.0 / math.sqrt(max(abs(symp_before), 1.0e-300))
        mode_vec = scale * rotated
        u = mode_vec[:n_cells]
        v = mode_vec[n_cells:]
        symp_after = float(np.real(np.sum(cell_weights * (np.abs(u) ** 2 - np.abs(v) ** 2))))
        residual_abs = float(np.linalg.norm(bdg @ mode_vec - eigenvalue * mode_vec))
        residual_rel = residual_abs / max(
            float(np.linalg.norm(bdg @ mode_vec)),
            abs(complex(eigenvalue)) * float(np.linalg.norm(mode_vec)),
            1.0,
        )

        density_response = np.real(np.conjugate(psi_vec) * u + psi_vec * v)
        drive_response = drive_vec * density_response
        drive_by_grid = _unflatten_w_major(drive_response, nr=nr, nw=nw)
        density_by_w = np.sum(radial_volumes[:, None] * drive_by_grid, axis=0)
        density_on_profile = _interp1_clamped(w_centers, density_by_w, points)
        signed_overlap = -_trapz(points, chi * density_on_profile)
        coupling = abs(signed_overlap)
        varpi = float(eigenvalue.real / float(background["constants"]["hbar"]))
        modal_spectrum.append(
            {
                "name": f"patha_b2a_bdg_mode_{mode_idx}",
                "varpi": varpi,
                "coupling": float(coupling),
                "raw_eigenvalue": float(eigenvalue.real),
                "eigenvalue_imag_abs": float(abs(eigenvalue.imag)),
                "relative_eigen_residual": float(residual_rel),
            }
        )
        if mode_idx > modes_to_export:
            continue

        residuals.append(residual_rel)
        raw_phi = density_on_profile / weight
        phi_norm = _weighted_norm(points, weight, raw_phi)
        if not math.isfinite(phi_norm) or phi_norm <= 0.0:
            raise RuntimeError(f"mode {mode_idx} has nonpositive axial profile norm")
        phi = raw_phi / math.sqrt(phi_norm)
        adapter_overlap = _trapz(points, weight * chi * phi)
        if adapter_overlap < 0.0:
            phi = -phi
            adapter_overlap = -adapter_overlap
        coupling = abs(signed_overlap)
        if adapter_overlap <= 1.0e-300:
            raise RuntimeError(f"mode {mode_idx} has near-zero wall/profile overlap")
        lambda_b = coupling / adapter_overlap
        modes.append(
            {
                "name": f"patha_b2a_bdg_mode_{mode_idx}",
                "lambda_B": float(lambda_b),
                "varpi": varpi,
                "profile": phi.astype(np.float64).tolist(),
                "profile_values": phi.astype(np.float64).tolist(),
                "overlap_I_eta_phi": float(adapter_overlap),
                "coupling": float(coupling),
                "wall_drive_overlap_signed": float(signed_overlap),
                "raw_eigenvalue": float(eigenvalue.real),
                "eigenvalue_imag_abs": float(abs(eigenvalue.imag)),
                "relative_eigen_residual": float(residual_rel),
                "symplectic_norm_before_abs_scaling": float(symp_before),
                "symplectic_norm_after_abs_scaling": float(symp_after),
                "profile_norm_weighted": float(_weighted_norm(points, weight, phi)),
            }
        )

    moments = _moments_with_stieltjes(modes)
    modal_truncation = _build_modal_truncation_metadata(
        modal_spectrum,
        modes_to_export=modes_to_export,
    )
    packet: dict[str, Any] = {
        "schema": "stage1_patha_b2a_bdg_bundle/v1",
        "engine": "python_independent_crosscheck",
        "background_content_hash": background["content_hash"],
        "wall_input_hash": wall_input["content_hash"],
        "tau": float(background["constants"]["tau"]),
        "geometry": {
            "a": float(background["geometry"]["a"]),
            "L": float(background["geometry"]["L"]),
        },
        "constants": {
            key: float(background["constants"][key])
            for key in [
                "eos_K",
                "hbar",
                "particle_mass",
                "gauge_charge",
                "tau",
                "radial_wall_strength",
                "axial_trap_strength",
            ]
        },
        "grid": {
            "nr": int(grid["nr"]),
            "nw": int(grid["nw"]),
            "dr": float(grid["dr"]),
            "dw": float(grid["dw"]),
            "r_centers": np.asarray(grid["r_centers"], dtype=np.float64).tolist(),
            "w_centers": np.asarray(grid["w_centers"], dtype=np.float64).tolist(),
            "profile_points": points.tolist(),
        },
        "wall": {
            "source": wall_input["source"],
            "K": float(wall_input["K"]),
            "M": float(wall_input["M"]),
            "chi": np.asarray(wall_input["chi"], dtype=np.float64).tolist(),
            "weight": weight.tolist(),
            "normalization_chi_T_W_chi": float(wall_input["normalization_chi_T_W_chi"]),
            "trapz_weight_norm": float(wall_input["trapz_weight_norm"]),
            "orientation_integral": float(wall_input["orientation_integral"]),
        },
        "bdg_modes": modes,
        "bdg_moments": moments,
        "diagnostics": {
            "operator_matrix_dimension": int(bdg.shape[0]),
            "positive_real_modes_found": int(len(positive)),
            "modes_exported": int(modes_to_export),
            "max_relative_eigen_residual": float(max(residuals)),
            "modal_truncation": modal_truncation,
            "spatial_gauge_linf": float(diag["spatial_gauge_linf"]),
            "operator_term_map": [
                {
                    "term": "kinetic",
                    "implementation": "-(hbar^2/(2m))*(laplacian-l(l+1)/r^2), l(l+1)=6",
                    "source": "notes/moving_throat_pde_program_compact.md:1406 and existing mt15_02 form",
                },
                {
                    "term": "EOS",
                    "implementation": "h=(5K/4)rho^4, dh/drho=5K rho^3",
                    "source": "notes/moving_throat_pde_program_compact.md:578",
                },
                {
                    "term": "wall drive",
                    "implementation": "delta V_conf=-4 V_radial r^4/R0^5 eta",
                    "source": "notes/moving_throat_pde_program_compact.md:1080",
                },
            ],
        },
    }
    packet["content_hash"] = _full_stable_hash(packet)
    return packet


def solve_python_command(
    *,
    background_path: Path,
    nr: int,
    nw: int,
    profile_points: int,
    modes: int,
    run_root: Path,
) -> tuple[Path, Path, dict[str, Any]]:
    background = _load_json(background_path)
    wall_input, wall_path = make_wall_input(
        background,
        profile_points_count=profile_points,
        out_dir=run_root / "wall_inputs",
    )
    packet = solve_bdg_python(background, wall_input, nr=nr, nw=nw, modes_to_export=modes)
    out_dir = run_root / "python"
    out_dir.mkdir(parents=True, exist_ok=True)
    tau_label = _format_tau(float(background["constants"]["tau"]))
    out_path = out_dir / f"patha_b2a_python_tau_{tau_label}_nr_{nr}_nw_{nw}.json"
    out_path.write_text(
        json.dumps(packet, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    if wall_path is None:
        raise RuntimeError("wall input path was not written")
    return out_path, wall_path, packet


def _max_diffs(left: Sequence[float], right: Sequence[float]) -> dict[str, float]:
    left_arr = np.asarray(left, dtype=np.float64)
    right_arr = np.asarray(right, dtype=np.float64)
    diff = np.abs(left_arr - right_arr)
    denom = np.maximum(np.abs(right_arr), 1.0e-300)
    return {
        "max_abs": float(np.max(diff)) if diff.size else 0.0,
        "max_rel": float(np.max(diff / denom)) if diff.size else 0.0,
    }


def _mode_values(packet: Mapping[str, Any], key: str) -> list[float]:
    return [float(mode[key]) for mode in packet["bdg_modes"]]


def _moment_values(packet: Mapping[str, Any]) -> list[float]:
    m = packet["bdg_moments"]
    return [float(m[key]) for key in MOMENT_KEYS]


def compare_engine_packets(python_packet: Mapping[str, Any], mma_packet: Mapping[str, Any]) -> dict[str, Any]:
    py_modes = python_packet["bdg_modes"]
    mma_modes = mma_packet["bdg_modes"]
    if len(py_modes) != len(mma_modes):
        raise ValueError("Python/MMA mode count mismatch")
    return {
        "varpi": _max_diffs(_mode_values(python_packet, "varpi"), _mode_values(mma_packet, "varpi")),
        "coupling": _max_diffs(
            _mode_values(python_packet, "coupling"),
            _mode_values(mma_packet, "coupling"),
        ),
        "moments": _max_diffs(_moment_values(python_packet), _moment_values(mma_packet)),
    }


def _strict_diff_pass(diff: Mapping[str, float], *, abs_tol: float, rel_tol: float) -> bool:
    return bool(float(diff["max_abs"]) <= abs_tol and float(diff["max_rel"]) <= rel_tol)


def _dual_engine_gate(dual_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    max_dual = {
        "varpi_abs": max(row["varpi"]["max_abs"] for row in dual_rows),
        "varpi_rel": max(row["varpi"]["max_rel"] for row in dual_rows),
        "coupling_abs": max(row["coupling"]["max_abs"] for row in dual_rows),
        "coupling_rel": max(row["coupling"]["max_rel"] for row in dual_rows),
        "moments_abs": max(row["moments"]["max_abs"] for row in dual_rows),
        "moments_rel": max(row["moments"]["max_rel"] for row in dual_rows),
    }
    per_quantity = {
        "varpi": {
            "max_abs": max_dual["varpi_abs"],
            "max_rel": max_dual["varpi_rel"],
            **DUAL_ENGINE_TOLERANCES["varpi"],
        },
        "coupling": {
            "max_abs": max_dual["coupling_abs"],
            "max_rel": max_dual["coupling_rel"],
            **DUAL_ENGINE_TOLERANCES["coupling"],
        },
        "moments": {
            "max_abs": max_dual["moments_abs"],
            "max_rel": max_dual["moments_rel"],
            **DUAL_ENGINE_TOLERANCES["moments"],
        },
    }
    passed_by_quantity = {
        key: _strict_diff_pass(
            values,
            abs_tol=float(values["abs"]),
            rel_tol=float(values["rel"]),
        )
        for key, values in per_quantity.items()
    }
    return {
        "passed": bool(all(passed_by_quantity.values())),
        "max_dual": max_dual,
        "per_quantity": per_quantity,
        "passed_by_quantity": passed_by_quantity,
        "criteria": (
            "AND gate: varpi abs/rel <= 1e-9/1e-9, "
            "coupling abs/rel <= 1e-12/1e-8, "
            "B moments abs/rel <= 1e-16/1e-9"
        ),
    }


def _packet_moments_from_modes(packet: Mapping[str, Any]) -> dict[str, float]:
    return {key: float(value) for key, value in bdg_moments(packet["bdg_modes"]).items()}


def _modal_truncation_gate(
    final_packet: Mapping[str, Any],
    confirmation_packet: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    final_meta = final_packet["diagnostics"]["modal_truncation"]
    checks = [
        {
            "grid": f"{final_packet['grid']['nr']}x{final_packet['grid']['nw']}",
            "max_truncation_error_at_export": float(final_meta["max_truncation_error_at_export"]),
            "tolerance": float(final_meta["tolerance"]),
            "passed": bool(
                float(final_meta["max_truncation_error_at_export"])
                <= float(final_meta["tolerance"])
            ),
        }
    ]
    confirmation_meta = None
    if confirmation_packet is not None:
        confirmation_meta = confirmation_packet["diagnostics"]["modal_truncation"]
        checks.append(
            {
                "grid": f"{confirmation_packet['grid']['nr']}x{confirmation_packet['grid']['nw']}",
                "max_truncation_error_at_export": float(
                    confirmation_meta["max_truncation_error_at_export"]
                ),
                "tolerance": float(confirmation_meta["tolerance"]),
                "passed": bool(
                    float(confirmation_meta["max_truncation_error_at_export"])
                    <= float(confirmation_meta["tolerance"])
                ),
            }
        )

    recomputed = _packet_moments_from_modes(final_packet)
    reported = {key: float(final_packet["bdg_moments"][key]) for key in MOMENT_KEYS}
    reproduction = _max_diffs(
        [recomputed[key] for key in MOMENT_KEYS],
        [reported[key] for key in MOMENT_KEYS],
    )
    reproduction_passed = _strict_diff_pass(
        reproduction,
        abs_tol=CONSUMER_MOMENT_TOLERANCES["abs"],
        rel_tol=CONSUMER_MOMENT_TOLERANCES["rel"],
    )
    return {
        "passed": bool(all(check["passed"] for check in checks) and reproduction_passed),
        "checks": checks,
        "final": final_meta,
        "confirmation": confirmation_meta,
        "bundle_mode_reproduction": {
            "scope": "shape check only: bdg_modes[] reproduces stored B_n through the B1 formula",
            "recomputed": recomputed,
            "reported": reported,
            "diffs": reproduction,
            "passed": reproduction_passed,
        },
        "criteria": (
            "max_n |B_n(K)-B_n(all-positive)|/|B_n(all-positive)| <= "
            f"{MODAL_TRUNCATION_TOL:.1e} at the export grid and finer confirmation grid"
        ),
    }


def _b1_cross_engine_consumer_check(
    *,
    mma_packet: Mapping[str, Any],
    python_packet: Mapping[str, Any],
) -> dict[str, Any]:
    consumer_from_mma_modes = _packet_moments_from_modes(mma_packet)
    python_converged = {key: float(python_packet["bdg_moments"][key]) for key in MOMENT_KEYS}
    diffs = _max_diffs(
        [consumer_from_mma_modes[key] for key in MOMENT_KEYS],
        [python_converged[key] for key in MOMENT_KEYS],
    )
    passed = _strict_diff_pass(
        diffs,
        abs_tol=CONSUMER_MOMENT_TOLERANCES["abs"],
        rel_tol=CONSUMER_MOMENT_TOLERANCES["rel"],
    )
    return {
        "consumer_from_mma_modes": consumer_from_mma_modes,
        "python_converged_moments": python_converged,
        "diffs": diffs,
        "passed": passed,
        "criteria": "B1 bdg_moments(MMA bdg_modes) vs Python converged B_n, abs AND rel",
    }


def _max_abs_rel_arrays(left: Sequence[float], right: Sequence[float]) -> dict[str, float]:
    return _max_diffs(left, right)


def _background_tau_movement(
    background1: Mapping[str, Any],
    background2: Mapping[str, Any],
) -> dict[str, Any]:
    psi1 = np.asarray(background1["fields"]["psi_R0"], dtype=np.float64) + 1j * np.asarray(
        background1["fields"]["psi_I0"], dtype=np.float64
    )
    psi2 = np.asarray(background2["fields"]["psi_R0"], dtype=np.float64) + 1j * np.asarray(
        background2["fields"]["psi_I0"], dtype=np.float64
    )
    rho1 = np.abs(psi1) ** 2
    rho2 = np.abs(psi2) ** 2
    r0_1 = np.asarray(background1["fields"]["R0_w"], dtype=np.float64)
    r0_2 = np.asarray(background2["fields"]["R0_w"], dtype=np.float64)
    a0_1 = np.asarray(background1["fields"]["A_00"], dtype=np.float64)
    a0_2 = np.asarray(background2["fields"]["A_00"], dtype=np.float64)
    return {
        "rho0": _max_abs_rel_arrays(rho1.ravel(), rho2.ravel()),
        "R0_w": _max_abs_rel_arrays(r0_1.ravel(), r0_2.ravel()),
        "A_00": _max_abs_rel_arrays(a0_1.ravel(), a0_2.ravel()),
        "summary": "matter density, wall radius profile, and scalar potential movement for tau=1 to tau=2",
    }


def _spatial_error_bars(
    packet: Mapping[str, Any],
    reference_packet: Mapping[str, Any],
) -> dict[str, Any]:
    moment_rel = {
        key: _relative_change(
            float(packet["bdg_moments"][key]),
            float(reference_packet["bdg_moments"][key]),
        )
        for key in MOMENT_KEYS
    }
    varpi_rel = [
        _relative_change(a, b)
        for a, b in zip(_mode_values(packet, "varpi"), _mode_values(reference_packet, "varpi"))
    ]
    varpi_primary3_rel = [
        _relative_change(a, b)
        for a, b in zip(_mode_values(packet, "varpi")[:3], _mode_values(reference_packet, "varpi")[:3])
    ]
    return {
        "packet_grid": f"{packet['grid']['nr']}x{packet['grid']['nw']}",
        "reference_grid": f"{reference_packet['grid']['nr']}x{reference_packet['grid']['nw']}",
        "varpi_max_rel": max(varpi_rel) if varpi_rel else 0.0,
        "varpi_exported_max_rel": max(varpi_rel) if varpi_rel else 0.0,
        "varpi_primary3_max_rel": max(varpi_primary3_rel) if varpi_primary3_rel else 0.0,
        "B_moments_rel": moment_rel,
        "B_moments_max_rel": max(moment_rel.values()),
        "scope": "relative difference between exported grid and the spatial reference grid",
    }


def _default_packet_path(run_root: Path, engine: str, tau: float, nr: int, nw: int) -> Path:
    tau_label = _format_tau(tau)
    if engine == "python":
        return run_root / "python" / f"patha_b2a_python_tau_{tau_label}_nr_{nr}_nw_{nw}.json"
    if engine == "mma":
        return run_root / "mathematica" / f"patha_b2a_mma_tau_{tau_label}_nr_{nr}_nw_{nw}.json"
    raise ValueError(engine)


def _default_background_latest(run_root: Path, tau: float) -> Path:
    return run_root / "backgrounds" / f"patha_b2a_closed_background_tau_{_format_tau(tau)}_latest.json"


def _grid_rows_from_defaults(run_root: Path, tau: float, grids: Sequence[tuple[int, int]]) -> list[dict[str, Any]]:
    rows = []
    for nr, nw in grids:
        rows.append(
            {
                "nr": nr,
                "nw": nw,
                "python": _load_json(_default_packet_path(run_root, "python", tau, nr, nw)),
                "mma": _load_json(_default_packet_path(run_root, "mma", tau, nr, nw)),
            }
        )
    return rows


def _relative_change(new: float, old: float) -> float:
    return abs(new - old) / max(abs(new), 1.0e-300)


def _grid_convergence(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    table: list[dict[str, Any]] = []
    previous: Mapping[str, Any] | None = None
    for row in rows:
        packet = row["mma"]
        moments = packet["bdg_moments"]
        values = {
            "grid": f"{row['nr']}x{row['nw']}",
            "varpi": _mode_values(packet, "varpi"),
            "B0": float(moments["B0"]),
            "B2": float(moments["B2"]),
            "B4": float(moments["B4"]),
        }
        if previous is not None:
            prev_packet = previous["mma"]
            prev_m = prev_packet["bdg_moments"]
            values["rel_change_varpi_max"] = max(
                _relative_change(a, b)
                for a, b in zip(_mode_values(packet, "varpi"), _mode_values(prev_packet, "varpi"))
            )
            primary_varpi_pairs = list(
                zip(_mode_values(packet, "varpi")[:3], _mode_values(prev_packet, "varpi")[:3])
            )
            values["rel_change_varpi_primary3_max"] = max(
                _relative_change(a, b)
                for a, b in primary_varpi_pairs
            )
            for key in MOMENT_KEYS:
                values[f"rel_change_{key}"] = _relative_change(
                    float(moments[key]),
                    float(prev_m[key]),
                )
            values["rel_change_B_max"] = max(
                float(values[f"rel_change_{key}"])
                for key in MOMENT_KEYS
            )
        else:
            values["rel_change_varpi_max"] = None
            values["rel_change_varpi_primary3_max"] = None
            for key in MOMENT_KEYS:
                values[f"rel_change_{key}"] = None
            values["rel_change_B_max"] = None
        table.append(values)
        previous = row
    final = table[-1]
    passed = bool(
        final["rel_change_varpi_primary3_max"] is not None
        and final["rel_change_varpi_primary3_max"] <= 0.08
        and final["rel_change_B_max"] <= 0.20
    )
    return {
        "table": table,
        "spatial_error_bars": {
            "grid": final["grid"],
            "previous_grid": table[-2]["grid"] if len(table) >= 2 else None,
            "varpi_max_rel": final["rel_change_varpi_max"],
            "varpi_exported_max_rel": final["rel_change_varpi_max"],
            "varpi_primary3_max_rel": final["rel_change_varpi_primary3_max"],
            "B_moments_rel": {key: final[f"rel_change_{key}"] for key in MOMENT_KEYS},
            "B_moments_max_rel": final["rel_change_B_max"],
            "scope": "finest-vs-previous relative change on the refinement ladder",
        },
        "passed": passed,
        "criteria": "finest-vs-previous max rel change: primary-three varpi <= 0.08 and B moments <= 0.20; full exported varpi max is recorded as an error bar",
        "order_note": "No Richardson order is claimed on this non-doubling 6->8->10 ladder.",
    }


def _tau_sensitivity(tau1_packet: Mapping[str, Any], tau2_packet: Mapping[str, Any]) -> dict[str, Any]:
    varpi_diff = _max_diffs(_mode_values(tau1_packet, "varpi"), _mode_values(tau2_packet, "varpi"))
    moments_diff = _max_diffs(_moment_values(tau1_packet), _moment_values(tau2_packet))
    max_rel = max(varpi_diff["max_rel"], moments_diff["max_rel"])
    return {
        "tau1": float(tau1_packet["tau"]),
        "tau2": float(tau2_packet["tau"]),
        "varpi": varpi_diff,
        "moments": moments_diff,
        "max_rel": max_rel,
        "passed": bool(max_rel >= 1.0e-5),
        "criteria": "max relative movement in {varpi,B0,B2,B4} >= 1e-5",
    }


def validate_and_report(
    *,
    run_root: Path,
    report_path: Path,
    tau1: float,
    tau2: float,
    grids: Sequence[tuple[int, int]],
    final_grid: tuple[int, int],
    modal_check_grid: tuple[int, int] | None,
) -> dict[str, Any]:
    background1 = _load_json(_default_background_latest(run_root, tau1))
    background2 = _load_json(_default_background_latest(run_root, tau2))
    rows = _grid_rows_from_defaults(run_root, tau1, grids)
    tau2_python = _load_json(_default_packet_path(run_root, "python", tau2, *final_grid))
    tau2_mma = _load_json(_default_packet_path(run_root, "mma", tau2, *final_grid))
    final_row = rows[-1]
    final_mma = final_row["mma"]
    modal_check_row = None
    if modal_check_grid is not None and modal_check_grid != final_grid:
        modal_check_row = {
            "nr": modal_check_grid[0],
            "nw": modal_check_grid[1],
            "python": _load_json(_default_packet_path(run_root, "python", tau1, *modal_check_grid)),
            "mma": _load_json(_default_packet_path(run_root, "mma", tau1, *modal_check_grid)),
        }

    dual_rows = []
    dual_input_rows = [
        *rows,
        {"nr": final_grid[0], "nw": final_grid[1], "python": tau2_python, "mma": tau2_mma},
    ]
    if modal_check_row is not None:
        dual_input_rows.append(modal_check_row)
    for row in dual_input_rows:
        dual_rows.append(
            {
                "tau": float(row["python"]["tau"]),
                "grid": f"{row['nr']}x{row['nw']}",
                **compare_engine_packets(row["python"], row["mma"]),
            }
        )

    dual_gate = _dual_engine_gate(dual_rows)

    all_packets = [row["mma"] for row in rows] + [row["python"] for row in rows] + [tau2_mma, tau2_python]
    if modal_check_row is not None:
        all_packets.extend([modal_check_row["mma"], modal_check_row["python"]])
    max_eigen_residual = max(
        max(float(mode["relative_eigen_residual"]) for mode in packet["bdg_modes"])
        for packet in all_packets
    )
    eigen_passed = bool(max_eigen_residual <= 1.0e-8)

    grid = _grid_convergence(rows)
    final_moments = final_mma["bdg_moments"]
    stieltjes = float(final_moments["B0_B4_minus_B2_squared"])
    structural_passed = bool(stieltjes >= -1.0e-14)
    tau_move = _tau_sensitivity(final_mma, tau2_mma)
    tau_move["background"] = _background_tau_movement(background1, background2)
    modal_gate = _modal_truncation_gate(
        final_mma,
        None if modal_check_row is None else modal_check_row["mma"],
    )
    b1_check = _b1_cross_engine_consumer_check(
        mma_packet=final_mma,
        python_packet=final_row["python"],
    )
    spatial_confirmation = (
        None
        if modal_check_row is None
        else _spatial_error_bars(final_mma, modal_check_row["mma"])
    )
    error_budget = {
        "B_moments": {
            "spatial_ladder_rel": grid["spatial_error_bars"]["B_moments_rel"],
            "spatial_ladder_max_rel": grid["spatial_error_bars"]["B_moments_max_rel"],
            "spatial_confirmation_rel": (
                None if spatial_confirmation is None else spatial_confirmation["B_moments_rel"]
            ),
            "spatial_confirmation_max_rel": (
                None if spatial_confirmation is None else spatial_confirmation["B_moments_max_rel"]
            ),
            "modal_truncation_rel": final_mma["diagnostics"]["modal_truncation"][
                "truncation_error_at_export"
            ],
            "modal_truncation_max_rel": final_mma["diagnostics"]["modal_truncation"][
                "max_truncation_error_at_export"
            ],
        },
        "varpi": {
            "spatial_ladder_exported_max_rel": grid["spatial_error_bars"]["varpi_exported_max_rel"],
            "spatial_ladder_primary3_max_rel": grid["spatial_error_bars"]["varpi_primary3_max_rel"],
            "spatial_confirmation_max_rel": (
                None if spatial_confirmation is None else spatial_confirmation["varpi_max_rel"]
            ),
            "spatial_confirmation_primary3_max_rel": (
                None if spatial_confirmation is None else spatial_confirmation["varpi_primary3_max_rel"]
            ),
        },
        "spatial_ladder": grid["spatial_error_bars"],
        "spatial_confirmation": spatial_confirmation,
        "modal_truncation": {
            "grid": f"{final_mma['grid']['nr']}x{final_mma['grid']['nw']}",
            "K": int(final_mma["diagnostics"]["modes_exported"]),
            "all_positive_mode_count": int(
                final_mma["diagnostics"]["modal_truncation"]["all_positive_mode_count"]
            ),
            "tolerance": float(final_mma["diagnostics"]["modal_truncation"]["tolerance"]),
            "rel": final_mma["diagnostics"]["modal_truncation"]["truncation_error_at_export"],
            "max_rel": final_mma["diagnostics"]["modal_truncation"][
                "max_truncation_error_at_export"
            ],
        },
    }
    self_consistency_passed = bool(
        background1["residuals"]["self_consistent"]
        and background2["residuals"]["self_consistent"]
        and abs(float(background1["geometry"]["a"]) - FROZEN_A) <= 1.0e-15
        and abs(float(background1["geometry"]["L"]) - FROZEN_L) <= 1.0e-15
        and abs(float(background2["geometry"]["a"]) - FROZEN_A) <= 1.0e-15
        and abs(float(background2["geometry"]["L"]) - FROZEN_L) <= 1.0e-15
    )

    gates = [
        {
            "gate": "self_consistency",
            "passed": self_consistency_passed,
            "catches": "old smoke profile, stale background bundle, or M1a L=2/R_exit geometry",
        },
        {
            "gate": "dual_engine_agreement",
            "passed": dual_gate["passed"],
            "catches": "engine-specific matrix assembly, interpolation, normalization, overlap bug, or tiny-B abs-only false pass",
        },
        {
            "gate": "eigen_residual",
            "passed": eigen_passed,
            "catches": "wrong eigenvector/eigenvalue pair exported from either engine",
        },
        {
            "gate": "grid_convergence",
            "passed": grid["passed"],
            "catches": "under-resolved BdG spectrum or moments being shipped as converged",
        },
        {
            "gate": "modal_truncation",
            "passed": modal_gate["passed"],
            "catches": "under-truncated B_n moment sum with a fat high-mode tail",
        },
        {
            "gate": "structural_sanity",
            "passed": structural_passed,
            "catches": "sign/index/exponent bug in B0/B2/B4 moment assembly",
        },
        {
            "gate": "tau_sensitivity",
            "passed": tau_move["passed"],
            "catches": "stale background reuse instead of re-solving at the requested tau",
        },
        {
            "gate": "b1_consumer_moments",
            "passed": b1_check["passed"],
            "catches": "MMA packet whose B1-consumed modes no longer reproduce the independently assembled Python converged B_n",
        },
    ]
    passed = all(gate["passed"] for gate in gates)

    final_bundle = {
        "schema": "stage1_patha_b2a_validated_bdg_bundle/v1",
        "source": "Mathematica mt15_02 Path-A B2a output after Python dual-engine validation",
        "background": {
            "path": str(_default_background_latest(run_root, tau1)),
            "content_hash": background1["content_hash"],
            "residual_linf": float(background1["residuals"]["closed_stationary_linf"]),
            "smoke_residual_reference": SMOKE_RESIDUAL_REFERENCE,
        },
        "tau": float(tau1),
        "geometry": final_mma["geometry"],
        "constants": final_mma["constants"],
        "grid": final_mma["grid"],
        "wall": final_mma["wall"],
        "bdg_modes": final_mma["bdg_modes"],
        "bdg_moments": final_mma["bdg_moments"],
        "diagnostics": final_mma["diagnostics"],
        "error_budget": error_budget,
        "validation": {
            "passed": passed,
            "gates": gates,
            "dual_engine_rows": dual_rows,
            "dual_engine_gate": dual_gate,
            "grid_convergence": grid,
            "modal_truncation": modal_gate,
            "tau_sensitivity": tau_move,
            "b1_cross_engine_consumer_check": b1_check,
        },
    }
    final_bundle["content_hash"] = _full_stable_hash(final_bundle)
    bundle_dir = run_root / "bundles"
    bundle_dir.mkdir(parents=True, exist_ok=True)
    bundle_path = bundle_dir / f"patha_b2a_validated_bdg_bundle_tau_{_format_tau(tau1)}.json"
    bundle_path.write_text(
        json.dumps(final_bundle, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )

    write_report(
        report_path=report_path,
        passed=passed,
        background1=background1,
        background2=background2,
        final_bundle=final_bundle,
        bundle_path=bundle_path,
        gates=gates,
        dual_rows=dual_rows,
        dual_gate=dual_gate,
        max_eigen_residual=max_eigen_residual,
        grid=grid,
        modal_gate=modal_gate,
        modal_check_packet=None if modal_check_row is None else modal_check_row["mma"],
        tau_move=tau_move,
        b1_check=b1_check,
        error_budget=error_budget,
    )

    result = {
        "passed": passed,
        "report_path": str(report_path),
        "bundle_path": str(bundle_path),
        "gates": gates,
        "final_bundle_hash": final_bundle["content_hash"],
    }
    return result


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.12e}"
    if isinstance(value, list):
        return "[" + ", ".join(_fmt(float(v)) for v in value) + "]"
    return str(value)


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
    passed: bool,
    background1: Mapping[str, Any],
    background2: Mapping[str, Any],
    final_bundle: Mapping[str, Any],
    bundle_path: Path,
    gates: Sequence[Mapping[str, Any]],
    dual_rows: Sequence[Mapping[str, Any]],
    dual_gate: Mapping[str, Any],
    max_eigen_residual: float,
    grid: Mapping[str, Any],
    modal_gate: Mapping[str, Any],
    modal_check_packet: Mapping[str, Any] | None,
    tau_move: Mapping[str, Any],
    b1_check: Mapping[str, Any],
    error_budget: Mapping[str, Any],
) -> None:
    modes = final_bundle["bdg_modes"]
    moments = final_bundle["bdg_moments"]
    max_dual = dual_gate["max_dual"]
    modal_meta = final_bundle["diagnostics"]["modal_truncation"]
    lines: list[str] = []
    lines.append("# Path-A B2a BdG Derivation")
    lines.append("")
    lines.append(f"Overall B2a gate: {'PASS' if passed else 'FAIL'}")
    lines.append(f"Validated bundle: `{bundle_path}`")
    lines.append(f"Bundle content hash: `{final_bundle['content_hash']}`")
    lines.append("")
    lines.append("## Load-Bearing Sources")
    lines.append("")
    lines.append("- Decision-12 requires Reading A, deriving spectra on the Path-A background, not inheriting the M1b packet: `software/stage1_solver/decisions/12_pathA_b2_derive_backreaction_bundle.md:18` and `:69`.")
    lines.append("- Decision-11 freezes the Hooke family and geometry `a=1`, `L=37/20`: `software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md:45` and `:74`.")
    lines.append("- Decision-11 clarifies numeric `varpi`/mode data are derived outputs, not frozen numbers: `software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md:166`.")
    lines.append("- Canonical EOS is `P=K rho^5`, `U=(K/4)rho^5`, `h=(5K/4)rho^4`: `notes/moving_throat_pde_program_compact.md:576`.")
    lines.append("- Canonical BdG and wall source are the compact BdG skeleton and `delta V_conf`: `notes/moving_throat_pde_program_compact.md:1406` and `:1424`.")
    lines.append("")
    lines.append("## Background")
    lines.append("")
    lines.append(
        f"`tau=1` closed residual Linf `{float(background1['residuals']['closed_stationary_linf']):.12e}` "
        f"versus old smoke residual `{SMOKE_RESIDUAL_REFERENCE:.12e}`."
    )
    lines.append(
        f"`tau=2` closed residual Linf `{float(background2['residuals']['closed_stationary_linf']):.12e}`."
    )
    lines.append(
        f"Frozen geometry used by the exported background: `a={background1['geometry']['a']}`, "
        f"`L={background1['geometry']['L']}`."
    )
    lines.append("")
    lines.append("## Derived Tau-1 Bundle")
    lines.append("")
    lines.append(
        _markdown_table(
            ["mode", "varpi", "coupling", "overlap", "lambda_B", "eigen_residual"],
            [
                {
                    "mode": mode["name"],
                    "varpi": float(mode["varpi"]),
                    "coupling": float(mode["coupling"]),
                    "overlap": float(mode["overlap_I_eta_phi"]),
                    "lambda_B": float(mode["lambda_B"]),
                    "eigen_residual": float(mode["relative_eigen_residual"]),
                }
                for mode in modes
            ],
        )
    )
    lines.append("")
    lines.append(
        f"`B0={float(moments['B0']):.12e}`, `B2={float(moments['B2']):.12e}`, "
        f"`B4={float(moments['B4']):.12e}`."
    )
    lines.append(
        f"Structural value `B0*B4-B2^2={float(moments['B0_B4_minus_B2_squared']):.12e}`."
    )
    lines.append("")
    lines.append("## Modal Truncation")
    lines.append("")
    lines.append(
        f"Exported mode count `K={int(modal_meta['exported_mode_count'])}` was checked against "
        f"`{int(modal_meta['all_positive_mode_count'])}` positive modes with tolerance "
        f"`{float(modal_meta['tolerance']):.1e}`."
    )
    lines.append(f"Modal gate criterion: {modal_gate['criteria']}.")
    lines.append(
        _markdown_table(
            [
                "grid",
                "M",
                "B0",
                "B2",
                "B4",
                "rel_prev_max",
                "rel_all_B0",
                "rel_all_B2",
                "rel_all_B4",
                "rel_all_max",
            ],
            [
                {
                    "grid": f"{final_bundle['grid']['nr']}x{final_bundle['grid']['nw']}",
                    "M": row["label"],
                    "B0": row["B0"],
                    "B2": row["B2"],
                    "B4": row["B4"],
                    "rel_prev_max": row["max_rel_change_from_previous"],
                    "rel_all_B0": row["rel_error_vs_all_positive"]["B0"],
                    "rel_all_B2": row["rel_error_vs_all_positive"]["B2"],
                    "rel_all_B4": row["rel_error_vs_all_positive"]["B4"],
                    "rel_all_max": row["max_rel_error_vs_all_positive"],
                }
                for row in modal_meta["sweep"]
            ],
        )
    )
    if modal_check_packet is not None:
        check_meta = modal_check_packet["diagnostics"]["modal_truncation"]
        lines.append("")
        lines.append(
            f"Finer-grid modal confirmation `{modal_check_packet['grid']['nr']}x{modal_check_packet['grid']['nw']}`:"
        )
        lines.append(
            _markdown_table(
                [
                    "grid",
                    "M",
                    "B0",
                    "B2",
                    "B4",
                    "rel_prev_max",
                    "rel_all_B0",
                    "rel_all_B2",
                    "rel_all_B4",
                    "rel_all_max",
                ],
                [
                    {
                        "grid": f"{modal_check_packet['grid']['nr']}x{modal_check_packet['grid']['nw']}",
                        "M": row["label"],
                        "B0": row["B0"],
                        "B2": row["B2"],
                        "B4": row["B4"],
                        "rel_prev_max": row["max_rel_change_from_previous"],
                        "rel_all_B0": row["rel_error_vs_all_positive"]["B0"],
                        "rel_all_B2": row["rel_error_vs_all_positive"]["B2"],
                        "rel_all_B4": row["rel_error_vs_all_positive"]["B4"],
                        "rel_all_max": row["max_rel_error_vs_all_positive"],
                    }
                    for row in check_meta["sweep"]
                ],
            )
        )
    lines.append("")
    lines.append("## Dual Engine")
    lines.append("")
    lines.append(
        "Scope: the Python and Mathematica paths independently assemble the BdG matrix, solve the eigensystem, "
        "and perform the overlap quadrature. The closed background and B1 wall mode `chi` are shared Python "
        "inputs by directive pathA_09; those common-mode inputs were validated earlier by the Path-A closed "
        "background checks in chunks 1b/1c and the B1 analytic `chi` oracle."
    )
    lines.append("")
    lines.append(
        _markdown_table(
            [
                "tau",
                "grid",
                "varpi_abs",
                "varpi_rel",
                "coupling_abs",
                "coupling_rel",
                "moments_abs",
                "moments_rel",
            ],
            [
                {
                    "tau": row["tau"],
                    "grid": row["grid"],
                    "varpi_abs": row["varpi"]["max_abs"],
                    "varpi_rel": row["varpi"]["max_rel"],
                    "coupling_abs": row["coupling"]["max_abs"],
                    "coupling_rel": row["coupling"]["max_rel"],
                    "moments_abs": row["moments"]["max_abs"],
                    "moments_rel": row["moments"]["max_rel"],
                }
                for row in dual_rows
            ],
        )
    )
    lines.append("")
    lines.append(
        "Max dual-engine diffs: "
        f"varpi abs/rel `{max_dual['varpi_abs']:.3e}`/`{max_dual['varpi_rel']:.3e}`, "
        f"c abs/rel `{max_dual['coupling_abs']:.3e}`/`{max_dual['coupling_rel']:.3e}`, "
        f"B abs/rel `{max_dual['moments_abs']:.3e}`/`{max_dual['moments_rel']:.3e}`."
    )
    lines.append(f"Dual-engine criterion: {dual_gate['criteria']}.")
    lines.append(f"Max eigen residual across exported packets: `{max_eigen_residual:.12e}`.")
    lines.append("")
    lines.append("## Grid Convergence")
    lines.append("")
    lines.append(f"Criterion: {grid['criteria']}.")
    lines.append(
        _markdown_table(
            [
                "grid",
                "varpi",
                "B0",
                "B2",
                "B4",
                "rel_change_varpi_max",
                "rel_change_varpi_primary3_max",
                "rel_change_B0",
                "rel_change_B2",
                "rel_change_B4",
                "rel_change_B_max",
            ],
            grid["table"],
        )
    )
    lines.append(grid["order_note"])
    lines.append("")
    lines.append("## Tau Sensitivity")
    lines.append("")
    lines.append(
        f"`tau=1` vs `tau=2` max relative movement `{tau_move['max_rel']:.12e}` "
        f"({tau_move['criteria']})."
    )
    lines.append(
        f"Varpi max abs/rel `{tau_move['varpi']['max_abs']:.12e}`/`{tau_move['varpi']['max_rel']:.12e}`; "
        f"B max abs/rel `{tau_move['moments']['max_abs']:.12e}`/`{tau_move['moments']['max_rel']:.12e}`."
    )
    background_move = tau_move["background"]
    lines.append(
        f"Background movement: rho0 max abs/rel `{background_move['rho0']['max_abs']:.12e}`/"
        f"`{background_move['rho0']['max_rel']:.12e}`, R0(w) max abs/rel "
        f"`{background_move['R0_w']['max_abs']:.12e}`/`{background_move['R0_w']['max_rel']:.12e}`, "
        f"A00 max abs/rel `{background_move['A_00']['max_abs']:.12e}`/"
        f"`{background_move['A_00']['max_rel']:.12e}`."
    )
    lines.append(
        "B2c design input: doubling tau moves the matter background at sub-percent scale in this frozen "
        "Hooke family; tau leverage on `R_norm` is expected to be dominated by exact `K=tau*kappahat` and "
        "the wall/Maxwell sectors, not by the BdG `B_n` drift."
    )
    lines.append("")
    lines.append("## B1 Consumer Check")
    lines.append("")
    lines.append(
        f"`patha_extraction.bdg_moments` applied to MMA-engine `bdg_modes[]` versus Python-engine converged "
        f"`B_n` max abs/rel "
        f"`{b1_check['diffs']['max_abs']:.3e}`/`{b1_check['diffs']['max_rel']:.3e}`."
    )
    lines.append("")
    lines.append("## Error Budget")
    lines.append("")
    lines.append(
        "Recorded B-moment spatial ladder rel errors: "
        + ", ".join(
            f"{key} `{float(value):.3e}`"
            for key, value in error_budget["B_moments"]["spatial_ladder_rel"].items()
        )
        + f"; max `{float(error_budget['B_moments']['spatial_ladder_max_rel']):.3e}`."
    )
    if error_budget["B_moments"]["spatial_confirmation_rel"] is not None:
        lines.append(
            "Recorded B-moment finer-grid confirmation rel errors: "
            + ", ".join(
                f"{key} `{float(value):.3e}`"
                for key, value in error_budget["B_moments"]["spatial_confirmation_rel"].items()
            )
            + f"; max `{float(error_budget['B_moments']['spatial_confirmation_max_rel']):.3e}`."
        )
    lines.append(
        "Recorded B-moment modal truncation rel errors at K: "
        + ", ".join(
            f"{key} `{float(value):.3e}`"
            for key, value in error_budget["B_moments"]["modal_truncation_rel"].items()
        )
        + f"; max `{float(error_budget['B_moments']['modal_truncation_max_rel']):.3e}`."
    )
    lines.append(
        f"Recorded varpi spatial ladder max rel: exported `{float(error_budget['varpi']['spatial_ladder_exported_max_rel']):.3e}`, "
        f"primary-three `{float(error_budget['varpi']['spatial_ladder_primary3_max_rel']):.3e}`"
        + (
            ""
            if error_budget["varpi"]["spatial_confirmation_max_rel"] is None
            else (
                f"; finer-grid confirmation exported `{float(error_budget['varpi']['spatial_confirmation_max_rel']):.3e}`, "
                f"primary-three `{float(error_budget['varpi']['spatial_confirmation_primary3_max_rel']):.3e}`"
            )
        )
        + "."
    )
    lines.append("")
    lines.append("## Gates And Wrong Answers Caught")
    lines.append("")
    lines.append(
        _markdown_table(
            ["gate", "status", "catches"],
            [
                {
                    "gate": gate["gate"],
                    "status": "PASS" if gate["passed"] else "FAIL",
                    "catches": gate["catches"],
                }
                for gate in gates
            ],
        )
    )
    lines.append("")
    lines.append("## Adaptations")
    lines.append("")
    lines.append("- The old inline smoke `rho0` and M1a `L=2.0`, `R_exit=1.65` were replaced by the exported closed Path-A background and frozen `a=1`, `L=1.85`.")
    lines.append("- The BdG operator form was kept to the canonical matter-sector terms: l=2 kinetic, confinement, quintic enthalpy plus `rho dh/drho`, `-mu`, `q A0`, and the wall-drive kernel.")
    lines.append("- The v2_22a overlap algebra is used with B1's flat-normalized wall mode. The Path-A bundle records both `lambda_B`, `overlap_I_eta_phi`, and `coupling` so B1 can either adapt profiles or call `bdg_moments` directly.")
    lines.append("- The radial volume convention follows the Path-A finite-volume grid (`4*pi*r^2 dr`) and the closed wall-return source, rather than the old smoke script's ad hoc radial volume normalization.")
    lines.append("- No Maxwell transfer, `R_norm`, `R_pole`, `P2`, `P4`, or root-find is assembled in this chunk.")
    lines.append("")
    lines.append("## Files Created Or Modified")
    lines.append("")
    lines.append("- `software/stage1_solver/src/stage1_solver/patha_b2a_bdg.py`")
    lines.append("- `software/stage1_solver/mathematica/mt15_02_bdg_wall_derivation.wls`")
    lines.append("- `software/stage1_solver/mathematica/mt15_02_patha_b2a_bdg_wall_derivation.wls`")
    lines.append("- `software/stage1_solver/tests/test_patha_b2a_bdg.py`")
    lines.append(f"- `{report_path}`")
    run_root = bundle_path.parents[1]
    generated_files = sorted(
        path
        for path in run_root.rglob("*")
        if path.is_file()
        and path.parent.name in {"backgrounds", "wall_inputs", "python", "mathematica", "bundles"}
    )
    for path in generated_files:
        lines.append(f"- `{path}`")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Path-A B2a BdG derivation helpers")
    sub = parser.add_subparsers(dest="command", required=True)

    export = sub.add_parser("export-background")
    export.add_argument("--tau", type=float, required=True)
    export.add_argument("--grid", type=_parse_grid, default=DEFAULT_BACKGROUND_GRID)
    export.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)

    solve = sub.add_parser("solve-python")
    solve.add_argument("--background", type=Path, required=True)
    solve.add_argument("--nr", type=int, required=True)
    solve.add_argument("--nw", type=int, required=True)
    solve.add_argument("--profile-points", type=int, default=DEFAULT_PROFILE_POINTS)
    solve.add_argument("--modes", type=int, default=DEFAULT_BDG_MODES)
    solve.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)

    validate = sub.add_parser("validate-report")
    validate.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    validate.add_argument("--report-path", type=Path, default=DEFAULT_REPORT_PATH)
    validate.add_argument("--tau1", type=float, default=1.0)
    validate.add_argument("--tau2", type=float, default=2.0)
    validate.add_argument("--grids", default="6,6;8,8;10,10")
    validate.add_argument("--final-grid", type=_parse_grid, default=(10, 10))
    validate.add_argument("--modal-check-grid", type=_parse_grid, default=(12, 12))

    args = parser.parse_args()
    if args.command == "export-background":
        path, bundle = export_background(tau=args.tau, grid_level=args.grid, run_root=args.run_root)
        print(f"background_path: {path}")
        print(f"content_hash: {bundle['content_hash']}")
        print(f"tau: {bundle['constants']['tau']}")
        print(f"grid: {bundle['grid']['nr']}x{bundle['grid']['nw']}")
        print(f"closed_stationary_linf: {bundle['residuals']['closed_stationary_linf']:.12e}")
        print(f"self_consistent: {bundle['residuals']['self_consistent']}")
        return 0 if bundle["residuals"]["self_consistent"] else 1

    if args.command == "solve-python":
        path, wall_path, packet = solve_python_command(
            background_path=args.background,
            nr=args.nr,
            nw=args.nw,
            profile_points=args.profile_points,
            modes=args.modes,
            run_root=args.run_root,
        )
        print(f"python_packet_path: {path}")
        print(f"wall_input_path: {wall_path}")
        print(f"content_hash: {packet['content_hash']}")
        print(f"tau: {packet['tau']}")
        print(f"grid: {args.nr}x{args.nw}")
        print(f"varpi: {[mode['varpi'] for mode in packet['bdg_modes']]}")
        print(f"B0_B2_B4: {[packet['bdg_moments'][key] for key in ('B0', 'B2', 'B4')]}")
        print(f"max_eigen_residual: {packet['diagnostics']['max_relative_eigen_residual']:.12e}")
        return 0

    if args.command == "validate-report":
        grids = [_parse_grid(piece) for piece in args.grids.split(";") if piece]
        result = validate_and_report(
            run_root=args.run_root,
            report_path=args.report_path,
            tau1=args.tau1,
            tau2=args.tau2,
            grids=grids,
            final_grid=args.final_grid,
            modal_check_grid=args.modal_check_grid,
        )
        print(f"report_path: {result['report_path']}")
        print(f"bundle_path: {result['bundle_path']}")
        print(f"final_bundle_hash: {result['final_bundle_hash']}")
        print(f"passed: {result['passed']}")
        for gate in result["gates"]:
            print(f"{gate['gate']}: {'PASS' if gate['passed'] else 'FAIL'}")
        return 0 if result["passed"] else 1

    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
