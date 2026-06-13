"""Step 8a grouped-real P2 tangent operator and static surrogate solve.

The tangent is deliberately wired as the Jacobian of the Step-3 coupled
residual, with only the compact-program l=2 angular pieces appended.  The
physical matter/gauge-to-wall source is not present here: compact lines
1377-1381 mark that coupled renormalization open, while lines 1441-1451 keep
``S_eta^(psi,A)`` schematic.  Step 8a therefore uses only a target-blind
external surrogate forcing.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .backend import configure_backend, jvp
from .boundaries import BoundaryCondition
from .config import (
    BranchSmokeConfig,
    HarnessConfig,
    P2TangentConfig,
    PreconditionerConfig,
    TensorGridSpec,
    WallGridSpec,
)
from .convergence import observed_order_from_three, richardson_estimate
from .coupled_branch import (
    CoupledFields,
    _create_branch_grid,
    branch_boundary_conditions,
    coupled_pde_residual,
    geometry_radius_torch,
    localization_weight_torch,
    pack_coupled_fields,
    resample_branch_state,
    run_branch_continuation,
    unpack_coupled_fields,
)
from .grid import TensorProductGrid, WallGrid
from .manifest import write_manifest
from .mms import weighted_l2_norm
from .newton import NewtonResult, solve_newton_jvp
from .operators import (
    spherical_ell_factor,
    wall_s_eta_operator,
)
from .preconditioners import make_p2_tangent_colored_sparse_jacobian_lu_factory


SOURCE_CITATIONS = {
    "linearized_skeleton": "notes/moving_throat_pde_program_compact.md lines 1383-1455",
    "wall_modal_split": "notes/moving_throat_pde_program_compact.md lines 1298-1308",
    "p2_degeneracy": "notes/moving_throat_pde_program_compact.md lines 1340-1364",
    "open_backreaction": "notes/moving_throat_pde_program_compact.md lines 1377-1381",
    "wall_to_matter": "notes/moving_throat_pde_program_compact.md lines 1075-1089",
    "bulk_angular_reduction": "Step-8a directive bulk angular reduction paragraph",
    "maxwell_scalarization": (
        "Step-8a engineering-smoke retained-lane scalarization; the full vector-harmonic "
        "Maxwell l=2 reduction is deferred."
    ),
}

P2_OBSERVABLE_LABELS = {
    "total_response_l2": "total raw response L2",
    "matter_interior_response_l2": "matter perturbation interior L2",
    "scalar_gauge_response_l2": "A0 perturbation L2",
    "spatial_gauge_response_l2": "spatial gauge perturbation L2",
    "wall_eta_l2": "wall eta L2",
    "wall_lower_endpoint_eta_abs": "lower-wall eta absolute value",
    "static_residual_linf": "static tangent residual Linf",
}


@dataclass(frozen=True)
class P2TangentFields:
    psi_real: torch.Tensor
    psi_imag: torch.Tensor
    a0: torch.Tensor
    ar: torch.Tensor
    aw: torch.Tensor
    eta: torch.Tensor


def _cell_count(grid: TensorProductGrid) -> int:
    return grid.spec.nr * grid.spec.nw


def pack_p2_tangent_fields(fields: P2TangentFields) -> torch.Tensor:
    return torch.cat(
        [
            fields.psi_real.reshape(-1),
            fields.psi_imag.reshape(-1),
            fields.a0.reshape(-1),
            fields.ar.reshape(-1),
            fields.aw.reshape(-1),
            fields.eta.reshape(-1),
        ]
    )


def unpack_p2_tangent_fields(
    state: torch.Tensor,
    grid: TensorProductGrid,
) -> P2TangentFields:
    n = _cell_count(grid)
    expected = 5 * n + grid.spec.nw
    if state.numel() != expected:
        raise ValueError(f"Expected P2 tangent state with {expected} entries, got {state.numel()}")
    shape = (grid.spec.nr, grid.spec.nw)
    return P2TangentFields(
        psi_real=state[0:n].reshape(shape),
        psi_imag=state[n : 2 * n].reshape(shape),
        a0=state[2 * n : 3 * n].reshape(shape),
        ar=state[3 * n : 4 * n].reshape(shape),
        aw=state[4 * n : 5 * n].reshape(shape),
        eta=state[5 * n : 5 * n + grid.spec.nw],
    )


def zero_p2_tangent_state(grid: TensorProductGrid) -> torch.Tensor:
    return torch.zeros(5 * _cell_count(grid) + grid.spec.nw, dtype=grid.dtype, device=grid.device)


def _wall_grid_from_tensor_grid(grid: TensorProductGrid) -> WallGrid:
    return WallGrid.create(
        WallGridSpec(w_min=grid.spec.w_min, w_max=grid.spec.w_max, nw=grid.spec.nw),
        dtype=grid.dtype,
        device=grid.device,
    )


def p2_wall_coefficients(
    grid: TensorProductGrid,
    p2_cfg: P2TangentConfig,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Engineering-smoke wall packet for the live ``K_eta+6 T_Omega`` term.

    Compact lines 1298-1308 and 1441-1451 pin the operator form.  The numeric
    coefficients remain Step-8a free-choice placeholders, labelled in the
    config and report.
    """

    wall_grid = _wall_grid_from_tensor_grid(grid)
    length = grid.spec.w_max - grid.spec.w_min
    x_faces = (wall_grid.w_faces - grid.spec.w_min) / length
    x_centers = (wall_grid.w_centers - grid.spec.w_min) / length
    t_w_faces = p2_cfg.placeholder_t_w_base + p2_cfg.placeholder_t_w_sine_amp * torch.sin(
        2.0 * math.pi * x_faces
    )
    t_omega_centers = p2_cfg.placeholder_t_omega_base + (
        p2_cfg.placeholder_t_omega_cosine_amp * torch.cos(2.0 * math.pi * x_centers)
    )
    bump = 256.0 * x_centers**4 * (1.0 - x_centers) ** 4
    k_eta_centers = p2_cfg.placeholder_k_eta_base + p2_cfg.placeholder_k_eta_bump_amp * bump
    return t_w_faces, t_omega_centers, k_eta_centers


def wall_to_matter_coefficient_torch(
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
) -> torch.Tensor:
    """Return the Step-3 ratio-confinement derivative used by ``delta V_conf``.

    Compact lines 1075-1089 pin the wall-to-matter coupling direction through
    ``delta V_conf = -(V_wall'/ell_c) eta`` for a signed-distance confinement
    ``V(Sigma=r-R)``.  The Step-3 engineering-smoke residual instead uses the
    ratio placeholder ``V_radial (r/R_0(w))^4`` plus an ``R_0``-independent
    axial term.  Linearizing that implemented radial term under the displaced
    wall radius ``R=R_0+eta`` gives ``dV/dR_0 = -4 V_radial r^4/R_0^5`` and
    therefore ``delta V_conf = -coeff * eta`` with the positive coefficient
    returned here.
    """

    r = grid.r_centers[:, None]
    radius = geometry_radius_torch(grid.w_centers, cfg)[None, :]
    return 4.0 * cfg.radial_wall_strength * r**4 / radius**5


def p2_mode_residual(
    mode_state: torch.Tensor,
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    p2_cfg: P2TangentConfig,
    *,
    eos_K: float,
    boundaries: Any,
) -> torch.Tensor:
    """Full Step-8a l=2 mode residual before taking its Jacobian.

    The first five lanes reuse ``coupled_pde_residual`` from Step 3 and then
    append the l=2 angular reductions.  The sixth lane is the static wall
    operator with ``K_eta + l(l+1)T_Omega``.  No matter/gauge-to-wall source is
    inserted because compact lines 1377-1381 leave that kernel open.
    """

    mode = unpack_p2_tangent_fields(mode_state, grid)
    background_fields, mu = unpack_coupled_fields(
        background_state.detach(),
        grid,
        has_chemical_potential=True,
    )
    if mu is None:
        raise ValueError("P2 tangent requires a background branch chemical potential")

    base_field_state = pack_coupled_fields(background_fields)
    candidate_fields = CoupledFields(
        psi_real=background_fields.psi_real + mode.psi_real,
        psi_imag=background_fields.psi_imag + mode.psi_imag,
        a0=background_fields.a0 + mode.a0,
        ar=background_fields.ar + mode.ar,
        aw=background_fields.aw + mode.aw,
    )
    candidate_field_state = pack_coupled_fields(candidate_fields)
    base_residual = coupled_pde_residual(
        base_field_state,
        grid,
        cfg,
        eos_K=eos_K,
        chemical_potential=mu,
        boundaries=boundaries,
    )
    candidate_residual = coupled_pde_residual(
        candidate_field_state,
        grid,
        cfg,
        eos_K=eos_K,
        chemical_potential=mu,
        boundaries=boundaries,
    )
    delta_residual = candidate_residual - base_residual

    ell = spherical_ell_factor(p2_cfg.spherical_l)
    inv_r2 = 1.0 / (grid.r_centers[:, None] ** 2)

    # Compact lines 1298-1308 plus Delta_S2 Y_lm=-l(l+1)Y_lm: the matter
    # kinetic block gains +(hbar^2/2m) l(l+1) delta_psi/r^2 in the residual.
    matter_angular = (cfg.hbar**2 / (2.0 * cfg.particle_mass)) * ell * inv_r2

    # Step-3 ratio-confinement wall-to-matter source.  The reverse
    # matter/gauge-to-wall source is explicitly not supplied in Step 8a.
    delta_v_conf = -wall_to_matter_coefficient_torch(grid, cfg) * mode.eta[None, :]

    weight_centers = localization_weight_torch(grid.w_centers, cfg)
    maxwell_angular = ell * weight_centers[None, :] * inv_r2

    matter_re = (
        delta_residual[0]
        + matter_angular * mode.psi_real
        + delta_v_conf * background_fields.psi_real
    )
    matter_im = (
        delta_residual[1]
        + matter_angular * mode.psi_imag
        + delta_v_conf * background_fields.psi_imag
    )
    maxwell_0 = delta_residual[2] + maxwell_angular * mode.a0
    maxwell_r = delta_residual[3] + maxwell_angular * mode.ar
    maxwell_w = delta_residual[4] + maxwell_angular * mode.aw
    pde_residual = torch.stack([matter_re, matter_im, maxwell_0, maxwell_r, maxwell_w], dim=0)

    wall_grid = _wall_grid_from_tensor_grid(grid)
    t_w_faces, t_omega_centers, k_eta_centers = p2_wall_coefficients(grid, p2_cfg)
    wall_residual = wall_s_eta_operator(
        mode.eta,
        wall_grid,
        t_w_faces=t_w_faces,
        t_omega_centers=t_omega_centers,
        k_eta_centers=k_eta_centers,
        spherical_l=p2_cfg.spherical_l,
        lower_bc=BoundaryCondition.neumann(0.0),
        upper_bc=BoundaryCondition.neumann(0.0),
    )
    return torch.cat([pde_residual.reshape(-1), wall_residual.reshape(-1)])


def apply_p2_tangent(
    direction: torch.Tensor,
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    p2_cfg: P2TangentConfig,
    *,
    eos_K: float,
    boundaries: Any,
) -> torch.Tensor:
    """Apply ``T_l=2`` as the JVP of ``p2_mode_residual`` at zero mode state."""

    zero = torch.zeros_like(direction)
    residual_fn = lambda mode: p2_mode_residual(
        mode,
        background_state,
        grid,
        cfg,
        p2_cfg,
        eos_K=eos_K,
        boundaries=boundaries,
    )
    return jvp(residual_fn, zero, direction)


def p2_static_linear_residual(
    state: torch.Tensor,
    forcing: torch.Tensor,
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    p2_cfg: P2TangentConfig,
    *,
    eos_K: float,
    boundaries: Any,
) -> torch.Tensor:
    return (
        apply_p2_tangent(
            state,
            background_state,
            grid,
            cfg,
            p2_cfg,
            eos_K=eos_K,
            boundaries=boundaries,
        )
        - forcing
    )


def p2_tangent_fd_check(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    p2_cfg: P2TangentConfig,
    *,
    eos_K: float,
    boundaries: Any,
    epsilon: float,
    seed: int,
) -> dict[str, float]:
    generator = torch.Generator(device=grid.device)
    generator.manual_seed(seed)
    direction = torch.randn(
        zero_p2_tangent_state(grid).shape,
        dtype=grid.dtype,
        device=grid.device,
        generator=generator,
    )
    direction = direction / torch.linalg.norm(direction)
    tangent = apply_p2_tangent(
        direction,
        background_state,
        grid,
        cfg,
        p2_cfg,
        eos_K=eos_K,
        boundaries=boundaries,
    ).detach()
    zero = zero_p2_tangent_state(grid)
    residual_fn = lambda mode: p2_mode_residual(
        mode,
        background_state,
        grid,
        cfg,
        p2_cfg,
        eos_K=eos_K,
        boundaries=boundaries,
    )
    fd = (
        residual_fn((zero + epsilon * direction).detach()).detach()
        - residual_fn((zero - epsilon * direction).detach()).detach()
    ) / (2.0 * epsilon)
    diff = tangent - fd
    absolute = float(torch.linalg.norm(diff).detach().cpu().item())
    tangent_norm = float(torch.linalg.norm(tangent).detach().cpu().item())
    fd_norm = float(torch.linalg.norm(fd).detach().cpu().item())
    return {
        "epsilon": float(epsilon),
        "absolute_residual": absolute,
        "relative_residual": absolute / max(1.0, tangent_norm),
        "tangent_norm": tangent_norm,
        "fd_norm": fd_norm,
    }


def p2_surrogate_forcing(
    grid: TensorProductGrid,
    p2_cfg: P2TangentConfig,
) -> torch.Tensor:
    """Target-blind smooth forcing that excites all retained Step-8a lanes."""

    r = grid.r_centers[:, None] / grid.spec.r_max
    x = (grid.w_centers[None, :] - grid.spec.w_min) / (grid.spec.w_max - grid.spec.w_min)
    radial = r**2 * torch.exp(-0.75 * r**2)
    axial = 1.0 + 0.15 * torch.cos(math.pi * x) + 0.10 * torch.sin(2.0 * math.pi * x)
    amp = p2_cfg.surrogate_force_amplitude
    fields = P2TangentFields(
        psi_real=amp * radial * axial,
        psi_imag=0.7 * amp * radial * (1.0 + 0.2 * torch.sin(math.pi * x)),
        a0=0.5 * amp * radial * (1.0 - 0.1 * torch.cos(2.0 * math.pi * x)),
        ar=0.4 * amp * radial * (0.8 + x),
        aw=-0.3 * amp * radial * (1.1 - 0.2 * x),
        eta=(0.6 * amp)
        * (
            1.0
            + 0.25
            * torch.sin(
                math.pi
                * (grid.w_centers - grid.spec.w_min)
                / (grid.spec.w_max - grid.spec.w_min)
            )
        ),
    )
    return pack_p2_tangent_fields(fields)


def p2_response_observables(
    state: torch.Tensor,
    grid: TensorProductGrid,
    *,
    static_residual_linf: float | None = None,
) -> dict[str, float]:
    fields = unpack_p2_tangent_fields(state, grid)
    volumes = grid.cell_volumes
    wall_grid = _wall_grid_from_tensor_grid(grid)
    matter = torch.sqrt(
        torch.sum((fields.psi_real**2 + fields.psi_imag**2) * volumes)
    )
    interior = (
        (grid.r_centers[:, None] < 0.75 * grid.spec.r_max)
        & (
            grid.w_centers[None, :]
            < grid.spec.w_min + 0.75 * (grid.spec.w_max - grid.spec.w_min)
        )
    )
    matter_density = fields.psi_real**2 + fields.psi_imag**2
    matter_interior = torch.sqrt(torch.sum(matter_density[interior] * volumes[interior]))
    scalar = torch.sqrt(torch.sum(fields.a0**2 * volumes))
    spatial = torch.sqrt(torch.sum((fields.ar**2 + fields.aw**2) * volumes))
    wall = torch.sqrt(torch.sum(fields.eta**2 * wall_grid.cell_widths))
    total = torch.sqrt(matter**2 + scalar**2 + spatial**2 + wall**2)
    result = {
        "total_response_l2": float(total.detach().cpu().item()),
        "matter_interior_response_l2": float(matter_interior.detach().cpu().item()),
        "scalar_gauge_response_l2": float(scalar.detach().cpu().item()),
        "spatial_gauge_response_l2": float(spatial.detach().cpu().item()),
        "wall_eta_l2": float(wall.detach().cpu().item()),
        "wall_lower_endpoint_eta_abs": float(torch.abs(fields.eta[0]).detach().cpu().item()),
    }
    if static_residual_linf is not None:
        result["static_residual_linf"] = float(static_residual_linf)
    return result


def _newton_history_rows(result: NewtonResult) -> list[dict[str, Any]]:
    rows = []
    for row in result.history:
        payload = asdict(row)
        payload.pop("preconditioner_setup_seconds", None)
        if payload.get("preconditioner_info"):
            payload["preconditioner_info"] = {
                key: value
                for key, value in payload["preconditioner_info"].items()
                if key not in {"factor_nnz_l", "factor_nnz_u"}
            }
        rows.append(payload)
    return rows


def solve_static_p2_tangent(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    grid_name: str,
) -> tuple[torch.Tensor, dict[str, Any]]:
    p2_cfg = config.p2_tangent
    branch_cfg = config.branch
    boundaries = branch_boundary_conditions(branch_cfg)
    final_K = branch_cfg.continuation_K_values[-1]
    forcing = p2_surrogate_forcing(grid, p2_cfg)
    residual_fn = lambda x: p2_static_linear_residual(
        x,
        forcing,
        background_state,
        grid,
        branch_cfg,
        p2_cfg,
        eos_K=final_K,
        boundaries=boundaries,
    )
    x0 = zero_p2_tangent_state(grid)
    result = solve_newton_jvp(
        residual_fn,
        x0,
        p2_cfg.newton,
        preconditioner_factory=make_p2_tangent_colored_sparse_jacobian_lu_factory(grid),
    )
    residual = residual_fn(result.x).detach()
    residual_linf = float(torch.max(torch.abs(residual)).detach().cpu().item())
    summary = {
        "grid": grid_name,
        "nr": grid.spec.nr,
        "nw": grid.spec.nw,
        "dof": int(result.x.numel()),
        "converged": result.converged,
        "iterations": result.iterations,
        "initial_residual_norm": result.initial_residual_norm,
        "final_residual_norm": result.final_residual_norm,
        "tolerance": result.tolerance,
        "message": result.message,
        "newton_history": _newton_history_rows(result),
        "preconditioner": asdict(p2_cfg.newton.preconditioner),
        "forcing_l2": float(weighted_l2_norm(forcing, torch.ones_like(forcing)).detach().cpu().item()),
        "static_residual_linf": residual_linf,
        "surrogate_values": p2_response_observables(
            result.x,
            grid,
            static_residual_linf=residual_linf,
        ),
    }
    manifest = write_manifest(
        run_root=config.run_root,
        benchmark_name=p2_cfg.name,
        grid_name=grid_name,
        config=config.to_dict(),
        mesh=grid.to_dict(),
        results=summary,
        config_hash=config.config_hash(),
        solver_controls=asdict(p2_cfg.newton),
    )
    summary["manifest"] = str(manifest)
    return result.x.detach(), summary


def _branch_preconditioner_config() -> PreconditionerConfig:
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


def with_step8a_preconditioners(config: HarnessConfig) -> HarnessConfig:
    branch = config.branch
    return replace(
        config,
        branch=replace(
            branch,
            newton=replace(branch.newton, preconditioner=_branch_preconditioner_config()),
        ),
    )


def _strip_runtime_level(level: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in level.items()
        if key not in {"grid_object", "background_state", "tangent_state"}
    }


def _background_summary(summary: dict[str, Any]) -> dict[str, Any]:
    return {
        "grid": summary["grid"],
        "nr": summary["nr"],
        "nw": summary["nw"],
        "dof": summary["dof"],
        "converged": summary["converged"],
        "message": summary["message"],
        "final_residual_linf": summary["final_residual_linf"],
        "final_mass": summary["final_mass"],
        "chemical_potential": summary["chemical_potential"],
        "placeholder_label": summary["placeholder_label"],
        "preconditioner": summary["preconditioner"],
        "boundaries": summary["boundaries"],
        "manifest": summary["manifest"],
    }


def _verdict(observable: str, order: float | None, expected_order: float) -> str:
    if observable == "static_residual_linf":
        return "solver-floor diagnostic"
    if order is None:
        return "drifts"
    if order >= 0.8 * expected_order:
        return "expected-order convergence"
    if order > 0.25:
        return "reduced-order convergence"
    return "drifts"


def _summarize_p2_observables(
    level_rows: list[dict[str, Any]],
    *,
    ratio: int,
    expected_order: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if len(level_rows) < 1:
        return [], []
    names = list(level_rows[0]["surrogate_values"].keys())
    summary_rows: list[dict[str, Any]] = []
    table_rows: list[dict[str, Any]] = []
    for observable in names:
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
                "label": P2_OBSERVABLE_LABELS.get(observable, observable),
                "finest_grid": level_rows[-1]["grid"],
                "finest_dof": level_rows[-1]["dof"],
                "finest_value": values[-1],
                "last_observed_order": last_order,
                "richardson_estimate": continuum,
                "finest_error_estimate": errors[-1],
                "verdict": verdict,
                "floor_or_null_label": verdict,
            }
        )
        for index, row in enumerate(level_rows):
            table_rows.append(
                {
                    "observable": observable,
                    "label": P2_OBSERVABLE_LABELS.get(observable, observable),
                    "level": row["level"],
                    "grid": row["grid"],
                    "dof": row["dof"],
                    "value": values[index],
                    "successive_change": changes[index],
                    "observed_order": orders[index],
                    "richardson_estimate": continuum,
                    "error_estimate": errors[index],
                    "verdict": verdict,
                }
            )
    return summary_rows, table_rows


def run_p2_static_ladder(config: HarnessConfig) -> list[dict[str, Any]]:
    dtype = configure_backend(config.backend)
    levels: list[dict[str, Any]] = []
    previous_grid: TensorProductGrid | None = None
    previous_background: torch.Tensor | None = None
    for level_index, level in enumerate(config.p2_tangent.convergence_levels):
        grid = _create_branch_grid(config.branch, level, dtype=dtype, device=config.backend.device)
        initial = (
            resample_branch_state(previous_background, previous_grid, grid, config.branch)
            if previous_background is not None and previous_grid is not None
            else None
        )
        background_state, background = run_branch_continuation(
            config,
            grid,
            initial_state=initial,
            grid_name=f"p2_background_l{level_index}_nr_{level[0]}_nw_{level[1]}",
        )
        tangent_state, solve = solve_static_p2_tangent(
            background_state,
            grid,
            config,
            grid_name=f"p2_static_l{level_index}_nr_{level[0]}_nw_{level[1]}",
        )
        levels.append(
            {
                "level": level_index,
                "grid": solve["grid"],
                "nr": grid.spec.nr,
                "nw": grid.spec.nw,
                "spacing": max(grid.dr, grid.dw),
                "dof": solve["dof"],
                "background": _background_summary(background),
                "solve": solve,
                "converged": background["converged"] and solve["converged"],
                "surrogate_values": solve["surrogate_values"],
                "grid_object": grid,
                "background_state": background_state.detach(),
                "tangent_state": tangent_state.detach(),
            }
        )
        previous_grid = grid
        previous_background = background_state.detach()
    return levels


def p2_tangent_dense_matrix(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
) -> np.ndarray:
    p2_cfg = config.p2_tangent
    branch_cfg = config.branch
    boundaries = branch_boundary_conditions(branch_cfg)
    final_K = branch_cfg.continuation_K_values[-1]
    size = zero_p2_tangent_state(grid).numel()
    columns: list[np.ndarray] = []
    for column in range(size):
        basis = torch.zeros(size, dtype=grid.dtype, device=grid.device)
        basis[column] = 1.0
        applied = apply_p2_tangent(
            basis,
            background_state,
            grid,
            branch_cfg,
            p2_cfg,
            eos_K=final_K,
            boundaries=boundaries,
        )
        columns.append(applied.detach().cpu().numpy().astype(np.float64, copy=False))
    return np.column_stack(columns)


def p2_wellposedness_check(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
) -> dict[str, float]:
    matrix = p2_tangent_dense_matrix(background_state, grid, config)
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    return {
        "smallest_singular_value": float(np.min(singular_values)),
        "largest_singular_value": float(np.max(singular_values)),
        "condition_number": float(np.max(singular_values) / np.min(singular_values)),
        "state_size": int(matrix.shape[0]),
    }


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _digest_sanitize(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _digest_sanitize(item)
            for key, item in value.items()
            if key not in {"manifest", "machine_readable_table"}
        }
    if isinstance(value, list):
        return [_digest_sanitize(item) for item in value]
    return value


def _diagnostics_digest(results: dict[str, Any]) -> str:
    payload = {
        "mms": _digest_sanitize(results["mms"]),
        "tangent_fd_check": _digest_sanitize(results["tangent_fd_check"]),
        "wellposedness": _digest_sanitize(results["wellposedness"]),
        "level_rows": _digest_sanitize(results["level_rows"]),
        "observable_summary": _digest_sanitize(results["observable_summary"]),
        "pass_checks": results["pass_checks"],
        "asserted_checks": results["asserted_checks"],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def run_step8a(config: HarnessConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = HarnessConfig(
            run_root="software/stage1_solver/runs/step8a_p2_tangent",
            report_path="software/stage1_solver/reports/step8a_p2_tangent.md",
        )
    config = with_step8a_preconditioners(config)
    Path(config.run_root).mkdir(parents=True, exist_ok=True)

    from .mms_benchmarks import run_p2_centrifugal_mms, run_p2_maxwell_angular_mms

    mms = {
        "p2_centrifugal": run_p2_centrifugal_mms(config),
        "p2_maxwell_angular": run_p2_maxwell_angular_mms(config),
    }
    runtime_levels = run_p2_static_ladder(config)
    level_rows = [_strip_runtime_level(level) for level in runtime_levels]
    observable_summary, observable_table = _summarize_p2_observables(
        level_rows,
        ratio=config.p2_tangent.refinement_ratio,
        expected_order=config.p2_tangent.expected_order,
    )
    fd_level = runtime_levels[0]
    boundaries = branch_boundary_conditions(config.branch)
    final_K = config.branch.continuation_K_values[-1]
    tangent_fd_check = p2_tangent_fd_check(
        fd_level["background_state"],
        fd_level["grid_object"],
        config.branch,
        config.p2_tangent,
        eos_K=final_K,
        boundaries=boundaries,
        epsilon=config.p2_tangent.newton.finite_difference_jvp_epsilon,
        seed=config.jacobian_check_seed,
    )
    wellposedness = p2_wellposedness_check(
        fd_level["background_state"],
        fd_level["grid_object"],
        config,
    )
    non_floor_observables = [
        row
        for row in observable_summary
        if row["verdict"] not in {"null diagnostic", "solver-floor diagnostic"}
    ]
    pass_checks = {
        "tangent_matches_fd_jacobian": (
            tangent_fd_check["relative_residual"] <= config.p2_tangent.tangent_fd_relative_tol
            or tangent_fd_check["absolute_residual"] <= config.p2_tangent.tangent_fd_absolute_tol
        ),
        "p2_centrifugal_mms_order": mms["p2_centrifugal"]["passed"],
        "p2_maxwell_angular_mms_order": mms["p2_maxwell_angular"]["passed"],
        "static_tangent_solves_converged": all(row["converged"] for row in level_rows),
        "p2_static_operator_wellposed": (
            wellposedness["smallest_singular_value"]
            >= config.p2_tangent.smallest_singular_min
        ),
        "surrogate_observable_orders": bool(non_floor_observables)
        and all(
            row["last_observed_order"] is not None
            and row["last_observed_order"] >= config.p2_tangent.min_observable_order
            for row in non_floor_observables
        ),
    }
    asserted_checks = {
        "target_blind_surrogate_forcing_only_not_a_physics_gate": True,
        "physical_export_permitted_is_false_not_a_physics_gate": True,
        "matter_gauge_to_wall_source_deferred_not_a_physics_gate": True,
        "wall_constitutive_placeholders_labelled_not_a_physics_gate": True,
    }
    asserted_check_notes = {
        "target_blind_surrogate_forcing_only_not_a_physics_gate": (
            "The Step-8a RHS is generated by p2_surrogate_forcing and reads no target, "
            "reference packet, or extraction map."
        ),
        "physical_export_permitted_is_false_not_a_physics_gate": (
            "Step 8a emits no physical packet and does not import the firewalled physical model."
        ),
        "matter_gauge_to_wall_source_deferred_not_a_physics_gate": (
            "Compact lines 1377-1381 leave S_eta^(psi,A) open; Step 8a does not invent it."
        ),
        "wall_constitutive_placeholders_labelled_not_a_physics_gate": (
            config.p2_tangent.placeholder_label
        ),
    }
    results: dict[str, Any] = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "source_citations": SOURCE_CITATIONS,
        "method": {
            "scope": (
                "Step 8a static l=2 grouped-real P2 tangent only; no absorber, "
                "no driven omega response, no physical response packet."
            ),
            "tangent_definition": (
                "T_l=2 is apply_p2_tangent = JVP[p2_mode_residual](0, delta_u)."
            ),
            "surrogate_forcing": (
                "target-blind smooth external forcing f_ext that excites matter, "
                "Maxwell, and wall lanes; not S_eta^(psi,A)."
            ),
            "sponge_note": (
                "The Step-8a sponge remains a pointwise residual term.  If Step 8b "
                "adds a differential absorber, that stencil must enter MMS forcing first."
            ),
            "placeholder_label": config.p2_tangent.placeholder_label,
            "wall_placeholder_coefficients": {
                "mu_eta": config.p2_tangent.placeholder_mu_eta,
                "T_w_base": config.p2_tangent.placeholder_t_w_base,
                "T_w_sine_amp": config.p2_tangent.placeholder_t_w_sine_amp,
                "T_Omega_base": config.p2_tangent.placeholder_t_omega_base,
                "T_Omega_cosine_amp": config.p2_tangent.placeholder_t_omega_cosine_amp,
                "K_eta_base": config.p2_tangent.placeholder_k_eta_base,
                "K_eta_bump_amp": config.p2_tangent.placeholder_k_eta_bump_amp,
            },
        },
        "mms": mms,
        "tangent_fd_check": tangent_fd_check,
        "wellposedness": wellposedness,
        "level_rows": level_rows,
        "observable_summary": observable_summary,
        "observable_table": observable_table,
        "pass_checks": pass_checks,
        "asserted_checks": asserted_checks,
        "asserted_check_notes": asserted_check_notes,
        "passed": all(pass_checks.values()),
    }
    results["diagnostics_digest"] = _diagnostics_digest(results)
    table_path = Path(config.run_root) / config.p2_tangent.name / "p2_tangent_table.json"
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


def write_step8a_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# Step 8a P2 Tangent Operator")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append(f"Diagnostics digest: `{results['diagnostics_digest']}`")
    lines.append("")
    lines.append(
        "**Scope:** engineering-smoke, target-blind, static omega=0 grouped-real P2 tangent only. "
        "The observables below are raw-field numerical surrogates, not physical WP3 observables; "
        "there is no response-packet export and no extraction map."
    )
    lines.append("")
    lines.append("## Transliteration Sources")
    lines.append("")
    for key, value in results["source_citations"].items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("## Pinned Vs Open")
    lines.append("")
    lines.append(
        "Pinned in this build: centrifugal matter/scalar angular reduction, live "
        "`[K_eta+6T_Omega]` wall term, and the Step-3 ratio-confinement "
        "wall-to-matter `delta V_conf` linearized with respect to the displaced "
        "radius `R=R0+eta`."
    )
    lines.append("")
    lines.append(
        "Engineering-smoke retained Maxwell angular block: exact for scalar lanes "
        "`A0` and `Aw`; componentwise scalarization for the radial vector lane `Ar`, "
        "omitting the `-2 A_r/r^2` vector-Laplacian self-term and the truncated "
        "transverse vector-harmonic lanes. The full vector-harmonic `l=2` Maxwell "
        "reduction is a standard derivable downstream task and is deferred."
    )
    lines.append("")
    lines.append(
        "Genuinely open and deferred: the physical matter/gauge-to-wall "
        "`S_eta^(psi,A)` source in compact lines 1377-1381; Step 8a uses a "
        "target-blind surrogate `f_ext`."
    )
    lines.append("")
    lines.append("## MMS For New Operator Pieces")
    lines.append("")
    for key in ("p2_centrifugal", "p2_maxwell_angular"):
        section = results["mms"][key]
        lines.append(f"### {section['name']}")
        lines.append("")
        lines.append(section["description"])
        lines.append(f"Source: {section['continuum_source']}")
        lines.append(f"Forcing: {section['forcing_derivation']}")
        lines.append("")
        lines.append(
            _table(
                ["grid", "spacing", "error", "observed_order", "reference_norm"],
                section["rows"],
            )
        )
        lines.append("")
        lines.append("Checks:")
        for check, passed in section["pass_checks"].items():
            lines.append(f"- {check}: {'PASS' if passed else 'FAIL'}")
        lines.append("")
    lines.append("## Tangent Jacobian Gate")
    lines.append("")
    fd = results["tangent_fd_check"]
    lines.append(
        "Gate A central-FD check: "
        f"relative={fd['relative_residual']:.6e}, "
        f"absolute={fd['absolute_residual']:.6e}, epsilon={fd['epsilon']:.6e}."
    )
    lines.append("")
    well = results["wellposedness"]
    lines.append(
        "Static l=2 well-posedness on the small grid: "
        f"sigma_min={well['smallest_singular_value']:.6e}, "
        f"condition={well['condition_number']:.6e}."
    )
    lines.append("")
    lines.append("## Static Tangent Solves")
    lines.append("")
    solve_rows = [
        {
            "level": row["level"],
            "grid": row["grid"],
            "dof": row["dof"],
            "background_residual": row["background"]["final_residual_linf"],
            "iterations": row["solve"]["iterations"],
            "final_residual_norm": row["solve"]["final_residual_norm"],
            "static_residual_linf": row["solve"]["static_residual_linf"],
            "converged": row["converged"],
        }
        for row in results["level_rows"]
    ]
    lines.append(
        _table(
            [
                "level",
                "grid",
                "dof",
                "background_residual",
                "iterations",
                "final_residual_norm",
                "static_residual_linf",
                "converged",
            ],
            solve_rows,
        )
    )
    lines.append("")
    lines.append("## Surrogate Observable Convergence")
    lines.append("")
    lines.append(
        "These are target-blind raw-field functionals of `delta_u`; they are not "
        "physical P2 response observables."
    )
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
                "floor_or_null_label",
            ],
            results["observable_summary"],
        )
    )
    lines.append("")
    lines.append("## Counted Checks")
    lines.append("")
    for key, value in results["pass_checks"].items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")
    lines.append("")
    lines.append("## Asserted Checks")
    lines.append("")
    lines.append("These are excluded from `passed` and carry `_not_a_physics_gate` suffixes.")
    for key, value in results["asserted_checks"].items():
        note = results["asserted_check_notes"].get(key, "")
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'} - {note}")
    lines.append("")
    lines.append("## Provenance And Reproduction")
    lines.append("")
    lines.append(f"Machine-readable table: `{results['machine_readable_table']}`.")
    lines.append(
        "Run with: `PYTHONPATH=software/stage1_solver/src timeout 600 "
        "python -m stage1_solver.p2_tangent_harness`."
    )
    lines.append(
        "Forward note for Step 8b: if the absorber becomes a differential stencil, "
        "it must be included in the MMS continuum forcing before certification."
    )
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path
