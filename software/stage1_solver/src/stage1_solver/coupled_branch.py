"""Coupled matter-Maxwell engineering smoke for Step 3.

This module composes the Step-1/Step-2 operators into a single stationary
isotropic residual.  Numeric parameter values here are explicitly
engineering-smoke placeholders; they are not a physical branch packet.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
import resource
import time
from typing import Any

import numpy as np
import sympy as sp
import torch
import torch.nn.functional as torch_functional

from .backend import configure_backend, tensor
from .boundaries import BoundaryCondition
from .config import (
    BranchSmokeConfig,
    CoupledBranchMMSConfig,
    HarnessConfig,
    PreconditionerConfig,
    TensorGridSpec,
)
from .grid import TensorProductGrid
from .manifest import write_manifest
from .mms import run_convergence_study, weighted_l2_norm
from .newton import finite_difference_jvp_check, solve_newton_jvp
from .operators import (
    axisymmetric_vector_divergence,
    localized_maxwell_operator,
    tensor_center_gradient_r,
    tensor_center_gradient_w,
    tensor_laplacian,
)
from .physics import quintic_enthalpy
from .preconditioners import make_coupled_colored_sparse_jacobian_lu_factory


@dataclass(frozen=True)
class CoupledBoundarySet:
    matter_radial_outer: BoundaryCondition
    matter_w_lower: BoundaryCondition
    matter_w_upper: BoundaryCondition
    a0_radial_outer: BoundaryCondition
    a0_w_lower: BoundaryCondition
    a0_w_upper: BoundaryCondition

    def to_dict(self) -> dict[str, dict[str, float | str]]:
        return {
            "matter_radial_outer": self.matter_radial_outer.to_dict(),
            "matter_w_lower": self.matter_w_lower.to_dict(),
            "matter_w_upper": self.matter_w_upper.to_dict(),
            "a0_radial_outer": self.a0_radial_outer.to_dict(),
            "a0_w_lower": self.a0_w_lower.to_dict(),
            "a0_w_upper": self.a0_w_upper.to_dict(),
        }


@dataclass(frozen=True)
class CoupledFields:
    psi_real: torch.Tensor
    psi_imag: torch.Tensor
    a0: torch.Tensor
    ar: torch.Tensor
    aw: torch.Tensor


def _torch_from_numpy(values: np.ndarray, *, dtype: torch.dtype, device: str) -> torch.Tensor:
    return tensor(np.asarray(values), dtype=dtype, device=device)


def _cell_count(grid: TensorProductGrid) -> int:
    return grid.spec.nr * grid.spec.nw


def pack_coupled_fields(fields: CoupledFields, chemical_potential: torch.Tensor | None = None) -> torch.Tensor:
    pieces = [
        fields.psi_real.reshape(-1),
        fields.psi_imag.reshape(-1),
        fields.a0.reshape(-1),
        fields.ar.reshape(-1),
        fields.aw.reshape(-1),
    ]
    if chemical_potential is not None:
        pieces.append(chemical_potential.reshape(1))
    return torch.cat(pieces)


def unpack_coupled_fields(
    state: torch.Tensor,
    grid: TensorProductGrid,
    *,
    has_chemical_potential: bool,
) -> tuple[CoupledFields, torch.Tensor | None]:
    n = _cell_count(grid)
    expected = 5 * n + (1 if has_chemical_potential else 0)
    if state.numel() != expected:
        raise ValueError(f"Expected state with {expected} entries, got {state.numel()}")
    shape = (grid.spec.nr, grid.spec.nw)
    fields = CoupledFields(
        psi_real=state[0:n].reshape(shape),
        psi_imag=state[n : 2 * n].reshape(shape),
        a0=state[2 * n : 3 * n].reshape(shape),
        ar=state[3 * n : 4 * n].reshape(shape),
        aw=state[4 * n : 5 * n].reshape(shape),
    )
    chemical_potential = state[5 * n] if has_chemical_potential else None
    return fields, chemical_potential


def localization_weight_torch(w: torch.Tensor, cfg: BranchSmokeConfig | CoupledBranchMMSConfig) -> torch.Tensor:
    midpoint = 0.5 * (cfg.w_min + cfg.w_max)
    return cfg.localization_floor + cfg.localization_amplitude * torch.exp(
        -((w - midpoint) / cfg.localization_width) ** 2
    )


def geometry_radius_torch(w: torch.Tensor, cfg: BranchSmokeConfig | CoupledBranchMMSConfig) -> torch.Tensor:
    s = w - cfg.w_min
    return cfg.r_exit + (cfg.r_mouth - cfg.r_exit) * torch.exp(-s / cfg.geometry_decay_length)


def confinement_potential_torch(
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig | CoupledBranchMMSConfig,
) -> torch.Tensor:
    r = grid.r_centers[:, None]
    w = grid.w_centers[None, :]
    radius = geometry_radius_torch(grid.w_centers, cfg)[None, :]
    length = cfg.w_max - cfg.w_min
    midpoint = 0.5 * (cfg.w_min + cfg.w_max)
    radial_wall = cfg.radial_wall_strength * (r / radius) ** 4
    axial = 0.5 * cfg.axial_trap_strength * ((w - midpoint) / length) ** 2
    return radial_wall + axial


def boundary_sponge_profile_torch(
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig | CoupledBranchMMSConfig,
) -> torch.Tensor:
    """Smooth stationary boundary-layer profile; disabled unless configured.

    This is a coefficient layer, not a field clamp: it adds damping-like mass
    terms near the outer radial and upper axial exits while leaving the
    interior untouched.
    """

    if not getattr(cfg, "sponge_enabled", False):
        return torch.zeros(
            (grid.spec.nr, grid.spec.nw),
            dtype=grid.dtype,
            device=grid.device,
        )
    width = float(getattr(cfg, "sponge_width", 0.0))
    if width <= 0.0:
        raise ValueError("sponge_width must be positive when sponge_enabled=True")
    power = int(getattr(cfg, "sponge_power", 2))
    if power < 1:
        raise ValueError("sponge_power must be at least 1")
    r = grid.r_centers[:, None]
    w = grid.w_centers[None, :]
    radial = torch.clamp((r - (cfg.r_max - width)) / width, min=0.0, max=1.0)
    axial = torch.clamp((w - (cfg.w_max - width)) / width, min=0.0, max=1.0)
    radial = radial**power
    axial = axial**power
    return 1.0 - (1.0 - radial) * (1.0 - axial)


def branch_boundary_conditions(cfg: BranchSmokeConfig) -> CoupledBoundarySet:
    return CoupledBoundarySet(
        matter_radial_outer=BoundaryCondition.dirichlet(0.0),
        matter_w_lower=BoundaryCondition.neumann(0.0),
        matter_w_upper=BoundaryCondition.robin(
            alpha=cfg.matter_exit_impedance_alpha,
            beta=1.0,
            gamma=0.0,
        ),
        a0_radial_outer=BoundaryCondition.robin(
            alpha=cfg.a0_radial_impedance_alpha,
            beta=1.0,
            gamma=0.0,
        ),
        a0_w_lower=BoundaryCondition.dirichlet(0.0),
        a0_w_upper=BoundaryCondition.robin(
            alpha=cfg.a0_exit_impedance_alpha,
            beta=1.0,
            gamma=0.0,
        ),
    )


def mms_boundary_conditions() -> CoupledBoundarySet:
    zero_flux = BoundaryCondition.neumann(0.0)
    return CoupledBoundarySet(
        matter_radial_outer=zero_flux,
        matter_w_lower=zero_flux,
        matter_w_upper=zero_flux,
        a0_radial_outer=zero_flux,
        a0_w_lower=zero_flux,
        a0_w_upper=zero_flux,
    )


def _matter_number_current(
    fields: CoupledFields,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig | CoupledBranchMMSConfig,
) -> tuple[torch.Tensor, torch.Tensor]:
    grad_real_r = tensor_center_gradient_r(fields.psi_real, grid)
    grad_imag_r = tensor_center_gradient_r(fields.psi_imag, grid)
    grad_real_w = tensor_center_gradient_w(fields.psi_real, grid)
    grad_imag_w = tensor_center_gradient_w(fields.psi_imag, grid)
    density = fields.psi_real**2 + fields.psi_imag**2
    phase_current_r = fields.psi_real * grad_imag_r - fields.psi_imag * grad_real_r
    phase_current_w = fields.psi_real * grad_imag_w - fields.psi_imag * grad_real_w
    jr = (cfg.hbar / cfg.particle_mass) * phase_current_r - (
        cfg.gauge_charge / cfg.particle_mass
    ) * fields.ar * density
    jw = (cfg.hbar / cfg.particle_mass) * phase_current_w - (
        cfg.gauge_charge / cfg.particle_mass
    ) * fields.aw * density
    return jr, jw


def coupled_pde_residual(
    field_state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig | CoupledBranchMMSConfig,
    *,
    eos_K: float,
    chemical_potential: torch.Tensor | float,
    boundaries: CoupledBoundarySet,
) -> torch.Tensor:
    fields, _ = unpack_coupled_fields(field_state, grid, has_chemical_potential=False)
    psi = fields.psi_real.to(torch.complex128) + 1j * fields.psi_imag.to(torch.complex128)
    lap = tensor_laplacian(
        psi,
        grid,
        boundaries.matter_radial_outer,
        boundaries.matter_w_lower,
        boundaries.matter_w_upper,
    )
    grad_r = tensor_center_gradient_r(psi, grid)
    grad_w = tensor_center_gradient_w(psi, grid)
    div_a = axisymmetric_vector_divergence(fields.ar, fields.aw, grid)
    a_dot_grad = fields.ar * grad_r + fields.aw * grad_w
    alpha = cfg.gauge_charge / cfg.hbar
    covariant_lap = (
        lap
        - 1j * alpha * (2.0 * a_dot_grad + div_a * psi)
        - (alpha**2) * (fields.ar**2 + fields.aw**2) * psi
    )
    density = fields.psi_real**2 + fields.psi_imag**2
    potential = confinement_potential_torch(grid, cfg)
    sponge = boundary_sponge_profile_torch(grid, cfg)
    mu = torch.as_tensor(chemical_potential, dtype=fields.psi_real.dtype, device=fields.psi_real.device)
    matter = (
        -(cfg.hbar**2 / (2.0 * cfg.particle_mass)) * covariant_lap
        + potential * psi
        + getattr(cfg, "sponge_matter_strength", 0.0) * sponge * psi
        + quintic_enthalpy(density, eos_K) * psi
        + cfg.gauge_charge * fields.a0 * psi
        - mu * psi
    )

    weight_centers = localization_weight_torch(grid.w_centers, cfg)
    weight_faces = localization_weight_torch(grid.w_faces, cfg)
    maxwell = localized_maxwell_operator(
        fields.a0,
        fields.ar,
        fields.aw,
        grid,
        xi=cfg.xi,
        weight_w_centers=weight_centers,
        weight_w_faces=weight_faces,
        a0_radial_outer_bc=boundaries.a0_radial_outer,
        a0_w_lower_bc=boundaries.a0_w_lower,
        a0_w_upper_bc=boundaries.a0_w_upper,
    )
    jr_number, jw_number = _matter_number_current(fields, grid, cfg)
    charge_density = cfg.gauge_charge * density
    charge_current_r = cfg.gauge_charge * jr_number
    charge_current_w = cfg.gauge_charge * jw_number
    maxwell_residual = torch.stack(
        [
            maxwell[0]
            - cfg.mu0 * charge_density
            + getattr(cfg, "sponge_gauge_strength", 0.0) * sponge * fields.a0,
            maxwell[1]
            - cfg.mu0 * charge_current_r
            + getattr(cfg, "sponge_gauge_strength", 0.0) * sponge * fields.ar,
            maxwell[2]
            - cfg.mu0 * charge_current_w
            + getattr(cfg, "sponge_gauge_strength", 0.0) * sponge * fields.aw,
        ],
        dim=0,
    )
    return torch.stack(
        [
            torch.real(matter),
            torch.imag(matter),
            maxwell_residual[0],
            maxwell_residual[1],
            maxwell_residual[2],
        ],
        dim=0,
    )


def coupled_branch_residual(
    state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    *,
    eos_K: float,
    boundaries: CoupledBoundarySet,
) -> torch.Tensor:
    fields, chemical_potential = unpack_coupled_fields(
        state,
        grid,
        has_chemical_potential=True,
    )
    assert chemical_potential is not None
    field_state = pack_coupled_fields(fields)
    pde = coupled_pde_residual(
        field_state,
        grid,
        cfg,
        eos_K=eos_K,
        chemical_potential=chemical_potential,
        boundaries=boundaries,
    )
    density = fields.psi_real**2 + fields.psi_imag**2
    mass = torch.sum(density * grid.cell_volumes) - cfg.mass
    return torch.cat([pde.reshape(-1), mass.reshape(1)])


def _renormalize_state_mass(
    state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
) -> torch.Tensor:
    fields, mu = unpack_coupled_fields(state, grid, has_chemical_potential=True)
    density = fields.psi_real**2 + fields.psi_imag**2
    mass = torch.sum(density * grid.cell_volumes)
    scale = torch.sqrt(torch.as_tensor(cfg.mass, dtype=state.dtype, device=state.device) / mass)
    fields = CoupledFields(
        psi_real=fields.psi_real * scale,
        psi_imag=fields.psi_imag * scale,
        a0=fields.a0,
        ar=fields.ar,
        aw=fields.aw,
    )
    assert mu is not None
    return pack_coupled_fields(fields, mu)


def initial_branch_state(
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    *,
    dtype: torch.dtype,
    device: str,
) -> torch.Tensor:
    r = grid.r_centers[:, None]
    w = grid.w_centers[None, :]
    radius = geometry_radius_torch(grid.w_centers, cfg)[None, :]
    length = cfg.w_max - cfg.w_min
    midpoint = 0.5 * (cfg.w_min + cfg.w_max)
    psi_real = torch.exp(-0.5 * (r / radius) ** 2) * (
        1.0 + 0.05 * torch.cos(torch.pi * (w - midpoint) / length)
    )
    psi_imag = torch.zeros_like(psi_real)
    zeros = torch.zeros_like(psi_real)
    fields = CoupledFields(psi_real=psi_real, psi_imag=psi_imag, a0=zeros, ar=zeros, aw=zeros)
    mu = torch.as_tensor(cfg.initial_mu, dtype=dtype, device=device)
    return _renormalize_state_mass(pack_coupled_fields(fields, mu), grid, cfg)


def resample_branch_state(
    state: torch.Tensor,
    old_grid: TensorProductGrid,
    new_grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
) -> torch.Tensor:
    fields, mu = unpack_coupled_fields(state, old_grid, has_chemical_potential=True)
    assert mu is not None
    stacked = torch.stack(
        [fields.psi_real, fields.psi_imag, fields.a0, fields.ar, fields.aw],
        dim=0,
    )
    resized = torch_functional.interpolate(
        stacked.unsqueeze(0),
        size=(new_grid.spec.nr, new_grid.spec.nw),
        mode="bilinear",
        align_corners=False,
    ).squeeze(0)
    new_fields = CoupledFields(
        psi_real=resized[0],
        psi_imag=resized[1],
        a0=resized[2],
        ar=resized[3],
        aw=resized[4],
    )
    return _renormalize_state_mass(pack_coupled_fields(new_fields, mu), new_grid, cfg)


def _max_rss_mb() -> float:
    # Linux reports ru_maxrss in KiB.
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0


def _newton_history_rows(result: Any) -> list[dict[str, Any]]:
    return [asdict(row) for row in result.history]


def _coupled_preconditioner_factory(config: HarnessConfig, grid: TensorProductGrid):
    preconditioner = config.branch.newton.preconditioner
    if preconditioner.type == "none":
        return None
    if preconditioner.type == "colored_sparse_jacobian_lu":
        return make_coupled_colored_sparse_jacobian_lu_factory(grid)
    raise ValueError(f"Unsupported coupled preconditioner {preconditioner.type!r}")


def target_blind_surrogate_observables(
    state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    *,
    eos_K: float,
    final_residual_linf: float | None = None,
) -> dict[str, float]:
    """Extraction-free raw-field diagnostics for engineering-smoke studies."""

    fields, mu = unpack_coupled_fields(state, grid, has_chemical_potential=True)
    density = fields.psi_real**2 + fields.psi_imag**2
    volumes = grid.cell_volumes
    jr_number, jw_number = _matter_number_current(fields, grid, cfg)
    grad_real_r = tensor_center_gradient_r(fields.psi_real, grid)
    grad_imag_r = tensor_center_gradient_r(fields.psi_imag, grid)
    grad_real_w = tensor_center_gradient_w(fields.psi_real, grid)
    grad_imag_w = tensor_center_gradient_w(fields.psi_imag, grid)
    potential = confinement_potential_torch(grid, cfg)
    gradient_density = grad_real_r**2 + grad_imag_r**2 + grad_real_w**2 + grad_imag_w**2
    raw_field_norm_density = (
        fields.psi_real**2
        + fields.psi_imag**2
        + fields.a0**2
        + fields.ar**2
        + fields.aw**2
    )
    energy_like_density = (
        (cfg.hbar**2 / (2.0 * cfg.particle_mass)) * gradient_density
        + potential * density
        + 0.25 * eos_K * density**5
        + 0.5 * (fields.a0**2 + fields.ar**2 + fields.aw**2)
    )
    result = {
        "density_mass": float(torch.sum(density * volumes).detach().cpu().item()),
        "peak_density": float(torch.max(density).detach().cpu().item()),
        "min_density": float(torch.min(density).detach().cpu().item()),
        "raw_field_l2_norm": float(
            torch.sqrt(torch.sum(raw_field_norm_density * volumes)).detach().cpu().item()
        ),
        "scalar_gauge_l2": float(
            torch.sqrt(torch.sum(fields.a0**2 * volumes)).detach().cpu().item()
        ),
        "spatial_gauge_l2": float(
            torch.sqrt(torch.sum((fields.ar**2 + fields.aw**2) * volumes)).detach().cpu().item()
        ),
        "spatial_current_l2": float(
            torch.sqrt(torch.sum((jr_number**2 + jw_number**2) * volumes)).detach().cpu().item()
        ),
        "field_energy_like_integral": float(
            torch.sum(energy_like_density * volumes).detach().cpu().item()
        ),
        "chemical_potential": float(mu.detach().cpu().item()) if mu is not None else float("nan"),
    }
    if final_residual_linf is not None:
        result["final_residual_linf"] = float(final_residual_linf)
    return result


def run_branch_continuation(
    config: HarnessConfig,
    grid: TensorProductGrid,
    *,
    initial_state: torch.Tensor | None = None,
    grid_name: str,
) -> tuple[torch.Tensor, dict[str, Any]]:
    cfg = config.branch
    boundaries = branch_boundary_conditions(cfg)
    state = (
        initial_state
        if initial_state is not None
        else initial_branch_state(grid, cfg, dtype=grid.dtype, device=grid.device)
    )
    started = time.perf_counter()
    stages: list[dict[str, Any]] = []
    converged = True
    message = "continuation completed"
    shared_preconditioner_factory = _coupled_preconditioner_factory(config, grid)

    for eos_K in cfg.continuation_K_values:
        residual_fn = lambda x, eos_K=eos_K: coupled_branch_residual(
            x,
            grid,
            cfg,
            eos_K=eos_K,
            boundaries=boundaries,
        )
        result = solve_newton_jvp(
            residual_fn,
            state,
            cfg.newton,
            preconditioner_factory=(
                _coupled_preconditioner_factory(config, grid)
                if cfg.newton.preconditioner.rebuild_policy == "once_per_newton_solve"
                else shared_preconditioner_factory
            ),
        )
        state = result.x.detach()
        gmres_counts = [
            row.gmres_iterations for row in result.history if row.gmres_iterations is not None
        ]
        stages.append(
            {
                "eos_K": eos_K,
                "converged": result.converged,
                "iterations": result.iterations,
                "initial_residual_norm": result.initial_residual_norm,
                "final_residual_norm": result.final_residual_norm,
                "tolerance": result.tolerance,
                "message": result.message,
                "preconditioner": asdict(cfg.newton.preconditioner),
                "gmres_iterations": gmres_counts,
                "residual_history": [row.residual_norm for row in result.history],
                "newton_history": _newton_history_rows(result),
            }
        )
        if not result.converged:
            converged = False
            message = f"continuation stopped at eos_K={eos_K}: {result.message}"
            break

    elapsed = time.perf_counter() - started
    fields, mu = unpack_coupled_fields(state, grid, has_chemical_potential=True)
    final_residual = coupled_branch_residual(
        state,
        grid,
        cfg,
        eos_K=cfg.continuation_K_values[min(len(stages), len(cfg.continuation_K_values)) - 1],
        boundaries=boundaries,
    )
    final_density = fields.psi_real**2 + fields.psi_imag**2
    summary = {
        "name": cfg.name,
        "grid": grid_name,
        "nr": grid.spec.nr,
        "nw": grid.spec.nw,
        "dof": int(5 * grid.spec.nr * grid.spec.nw + 1),
        "converged": converged,
        "message": message,
        "wall_clock_seconds": elapsed,
        "peak_memory_mb": _max_rss_mb(),
        "final_residual_linf": float(torch.max(torch.abs(final_residual)).detach().cpu().item()),
        "final_mass": float(torch.sum(final_density * grid.cell_volumes).detach().cpu().item()),
        "chemical_potential": float(mu.detach().cpu().item()) if mu is not None else None,
        "stages": stages,
        "placeholder_label": cfg.placeholder_label,
        "preconditioner": asdict(cfg.newton.preconditioner),
        "boundaries": boundaries.to_dict(),
    }
    summary["surrogate_values"] = target_blind_surrogate_observables(
        state,
        grid,
        cfg,
        eos_K=cfg.continuation_K_values[min(len(stages), len(cfg.continuation_K_values)) - 1],
        final_residual_linf=summary["final_residual_linf"],
    )
    manifest = write_manifest(
        run_root=config.run_root,
        benchmark_name=cfg.name,
        grid_name=grid_name,
        config=config.to_dict(),
        mesh=grid.to_dict(),
        results=summary,
        config_hash=config.config_hash(),
        solver_controls=asdict(cfg.newton),
    )
    summary["manifest"] = str(manifest)
    return state, summary


def _coupled_mms_expressions() -> dict[str, Any]:
    r, radius = sp.symbols("r R", positive=True, real=True)
    w, w_min, w_max = sp.symbols("w w_min w_max", real=True)
    hbar, mass, charge, mu0, xi, eos_k, chemical_potential = sp.symbols(
        "hbar mass q mu0 xi eos_K chemical_potential",
        positive=True,
        real=True,
    )
    loc_width, loc_floor, loc_amp = sp.symbols(
        "lambda_z Z_floor Z_amp",
        positive=True,
        real=True,
    )
    r_mouth, r_exit, geom_decay = sp.symbols(
        "R_mouth R_exit ell_R",
        positive=True,
        real=True,
    )
    radial_wall_strength, axial_trap_strength = sp.symbols(
        "V_radial V_axial",
        positive=True,
        real=True,
    )
    length = w_max - w_min
    x = (w - w_min) / length
    midpoint = sp.Rational(1, 2) * (w_min + w_max)
    rb = (1 - (r / radius) ** 2) ** 4
    wb = 256 * x**4 * (1 - x) ** 4
    u = 1 + sp.Rational(1, 8) * rb + sp.Rational(3, 100) * wb
    v = sp.Rational(1, 10) * rb * wb + sp.Rational(1, 20) * (r / radius) ** 2 * rb
    a0 = sp.Rational(1, 7) * rb + sp.Rational(1, 11) * wb
    ar = r * rb * (sp.Rational(1, 9) + sp.Rational(1, 13) * wb)
    aw = sp.Rational(1, 8) * rb * wb + sp.Rational(1, 17) * (r / radius) ** 2 * rb

    z = loc_floor + loc_amp * sp.exp(-((w - midpoint) / loc_width) ** 2)
    r0 = r_exit + (r_mouth - r_exit) * sp.exp(-(w - w_min) / geom_decay)
    potential = radial_wall_strength * (r / r0) ** 4 + sp.Rational(1, 2) * axial_trap_strength * (
        (w - midpoint) / length
    ) ** 2
    rho = u**2 + v**2
    enthalpy = sp.Rational(5, 4) * eos_k * rho**4

    def lap(scalar: sp.Expr) -> sp.Expr:
        return sp.diff(scalar, r, 2) + 2 * sp.diff(scalar, r) / r + sp.diff(scalar, w, 2)

    div_a = sp.diff(r**2 * ar, r) / r**2 + sp.diff(aw, w)
    a2 = ar**2 + aw**2
    alpha = charge / hbar
    b_re = 2 * (ar * sp.diff(u, r) + aw * sp.diff(u, w)) + div_a * u
    b_im = 2 * (ar * sp.diff(v, r) + aw * sp.diff(v, w)) + div_a * v
    d2_re = lap(u) + alpha * b_im - alpha**2 * a2 * u
    d2_im = lap(v) - alpha * b_re - alpha**2 * a2 * v
    scalar_multiplier = potential + enthalpy + charge * a0 - chemical_potential
    matter_re = -(hbar**2 / (2 * mass)) * d2_re + scalar_multiplier * u
    matter_im = -(hbar**2 / (2 * mass)) * d2_im + scalar_multiplier * v

    jr_number = (hbar / mass) * (u * sp.diff(v, r) - v * sp.diff(u, r)) - (
        charge / mass
    ) * ar * rho
    jw_number = (hbar / mass) * (u * sp.diff(v, w) - v * sp.diff(u, w)) - (
        charge / mass
    ) * aw * rho
    scalar_lap = z * (sp.diff(a0, r, 2) + 2 * sp.diff(a0, r) / r) + sp.diff(
        z * sp.diff(a0, w), w
    )
    f_rw = sp.diff(aw, r) - sp.diff(ar, w)
    maxwell_0 = -scalar_lap - mu0 * charge * rho
    maxwell_r = -sp.diff(z * f_rw, w) + (1 / xi) * sp.diff(z * div_a, r) - (
        mu0 * charge * jr_number
    )
    maxwell_w = sp.diff(r**2 * z * f_rw, r) / r**2 + (1 / xi) * sp.diff(
        z * div_a, w
    ) - (mu0 * charge * jw_number)

    field_matrix = sp.Matrix([u, v, a0, ar, aw])
    residual_matrix = sp.Matrix([matter_re, matter_im, maxwell_0, maxwell_r, maxwell_w])
    args = (
        r,
        w,
        radius,
        w_min,
        w_max,
        hbar,
        mass,
        charge,
        mu0,
        xi,
        eos_k,
        chemical_potential,
        loc_width,
        loc_floor,
        loc_amp,
        r_mouth,
        r_exit,
        geom_decay,
        radial_wall_strength,
        axial_trap_strength,
    )
    return {
        "fields": sp.lambdify(args, field_matrix, "numpy"),
        "residual": sp.lambdify(args, residual_matrix, "numpy"),
        "field_expr": field_matrix,
        "residual_expr": residual_matrix,
    }


_COUPLED_MMS = _coupled_mms_expressions()


def _coupled_mms_args(
    rr: np.ndarray,
    ww: np.ndarray,
    cfg: CoupledBranchMMSConfig,
) -> tuple[Any, ...]:
    return (
        rr,
        ww,
        cfg.r_max,
        cfg.w_min,
        cfg.w_max,
        cfg.hbar,
        cfg.particle_mass,
        cfg.gauge_charge,
        cfg.mu0,
        cfg.xi,
        cfg.eos_K,
        cfg.chemical_potential,
        cfg.localization_width,
        cfg.localization_floor,
        cfg.localization_amplitude,
        cfg.r_mouth,
        cfg.r_exit,
        cfg.geometry_decay_length,
        cfg.radial_wall_strength,
        cfg.axial_trap_strength,
    )


def run_coupled_branch_mms(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.mms.coupled_branch
    full_config = config.to_dict()
    boundaries = mms_boundary_conditions()

    def build_level(level: tuple[int, int]):
        nr, nw = level
        grid = TensorProductGrid.create(
            TensorGridSpec(
                r_max=cfg.r_max,
                nr=nr,
                w_min=cfg.w_min,
                w_max=cfg.w_max,
                nw=nw,
            ),
            dtype=dtype,
            device=config.backend.device,
        )
        spacing = max(grid.dr, grid.dw)
        return grid, f"nr_{nr}_nw_{nw}", spacing, grid.to_dict(), grid.cell_volumes

    def evaluate_level(grid: TensorProductGrid):
        rr, ww = np.meshgrid(
            grid.r_centers.detach().cpu().numpy(),
            grid.w_centers.detach().cpu().numpy(),
            indexing="ij",
        )
        args = _coupled_mms_args(rr, ww, cfg)
        fields_np = np.asarray(_COUPLED_MMS["fields"](*args), dtype=np.float64).squeeze()
        exact_np = np.asarray(_COUPLED_MMS["residual"](*args), dtype=np.float64).squeeze()
        fields = CoupledFields(
            psi_real=_torch_from_numpy(fields_np[0], dtype=dtype, device=config.backend.device),
            psi_imag=_torch_from_numpy(fields_np[1], dtype=dtype, device=config.backend.device),
            a0=_torch_from_numpy(fields_np[2], dtype=dtype, device=config.backend.device),
            ar=_torch_from_numpy(fields_np[3], dtype=dtype, device=config.backend.device),
            aw=_torch_from_numpy(fields_np[4], dtype=dtype, device=config.backend.device),
        )
        discrete = coupled_pde_residual(
            pack_coupled_fields(fields),
            grid,
            cfg,
            eos_K=cfg.eos_K,
            chemical_potential=cfg.chemical_potential,
            boundaries=boundaries,
        )
        exact = _torch_from_numpy(exact_np, dtype=dtype, device=config.backend.device)
        jr, jw = _matter_number_current(fields, grid, cfg)
        density = fields.psi_real**2 + fields.psi_imag**2
        diagnostics = {
            "boundaries": boundaries.to_dict(),
            "density_min": float(torch.min(density).detach().cpu().item()),
            "density_max": float(torch.max(density).detach().cpu().item()),
            "spatial_current_l2": float(
                weighted_l2_norm(torch.stack([jr, jw], dim=0), grid.cell_volumes)
                .detach()
                .cpu()
                .item()
            ),
            "spatial_gauge_l2": float(
                weighted_l2_norm(torch.stack([fields.ar, fields.aw], dim=0), grid.cell_volumes)
                .detach()
                .cpu()
                .item()
            ),
        }
        return discrete, exact, diagnostics

    result = run_convergence_study(
        name=cfg.name,
        description=(
            "Full coupled stationary matter plus localized-Maxwell residual on the (r,w) grid: "
            "A_i enters D_i and the Maxwell vector blocks are sourced by q_star*j_psi."
        ),
        continuum_source=(
            "compact lines 556-583, 638-648, 651-659, 674-689; "
            "prereg lines 44-56."
        ),
        manufactured_field=str(_COUPLED_MMS["field_expr"]),
        forcing_derivation=(
            "SymPy evaluated the continuum coupled operator directly: covariant "
            "axisymmetric D_iD_i for matter and H=Z localized Maxwell minus the "
            "matter charge/current source."
        ),
        levels=cfg.grid_levels,
        build_level=build_level,
        evaluate_level=evaluate_level,
        config=full_config,
        run_root=config.run_root,
        min_observed_order=cfg.min_observed_order,
        final_error_max=cfg.final_error_max,
        config_hash=config.config_hash(),
    )
    result_dict = asdict(result)
    final = result_dict["rows"][-1]
    result_dict["pass_checks"]["cross_sector_gauge_nonzero"] = final["spatial_gauge_l2"] > 1.0e-3
    result_dict["pass_checks"]["cross_sector_current_nonzero"] = final["spatial_current_l2"] > 1.0e-3
    result_dict["passed"] = all(result_dict["pass_checks"].values())
    return result_dict


def _create_branch_grid(
    cfg: BranchSmokeConfig,
    level: tuple[int, int],
    *,
    dtype: torch.dtype,
    device: str,
) -> TensorProductGrid:
    nr, nw = level
    return TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.r_max, nr=nr, w_min=cfg.w_min, w_max=cfg.w_max, nw=nw),
        dtype=dtype,
        device=device,
    )


def run_resolution_ladder(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.branch
    rows: list[dict[str, Any]] = []
    previous_grid: TensorProductGrid | None = None
    previous_state: torch.Tensor | None = None
    stop_reason = "completed configured ladder"

    for level in cfg.ladder_levels:
        grid = _create_branch_grid(cfg, level, dtype=dtype, device=config.backend.device)
        if previous_state is not None and previous_grid is not None:
            initial = resample_branch_state(previous_state, previous_grid, grid, cfg)
        else:
            initial = None
        state, summary = run_branch_continuation(
            config,
            grid,
            initial_state=initial,
            grid_name=f"ladder_nr_{level[0]}_nw_{level[1]}",
        )
        gmres_counts = [
            count
            for stage in summary["stages"]
            for count in stage.get("gmres_iterations", [])
            if count is not None
        ]
        rows.append(
            {
                "grid": summary["grid"],
                "nr": summary["nr"],
                "nw": summary["nw"],
                "dof": summary["dof"],
                "wall_clock_seconds": summary["wall_clock_seconds"],
                "peak_memory_mb": summary["peak_memory_mb"],
                "newton_iterations": sum(stage["iterations"] for stage in summary["stages"]),
                "final_residual_linf": summary["final_residual_linf"],
                "gmres_iterations": gmres_counts,
                "converged": summary["converged"],
                "message": summary["message"],
                "manifest": summary["manifest"],
            }
        )
        if not summary["converged"]:
            stop_reason = summary["message"]
            break
        previous_grid = grid
        previous_state = state
        if summary["wall_clock_seconds"] > cfg.max_ladder_level_seconds:
            stop_reason = (
                f"stopped after {summary['grid']} because the level exceeded "
                f"{cfg.max_ladder_level_seconds:.1f}s"
            )
            break

    return {"rows": rows, "stop_reason": stop_reason}


def run_step3(config: HarnessConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = HarnessConfig(
            run_root="software/stage1_solver/runs/step3_coupled_branch_smoke",
            report_path="software/stage1_solver/reports/step3_coupled_branch_smoke.md",
        )
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(config.backend)
    coupled_mms = run_coupled_branch_mms(config)
    grid = _create_branch_grid(
        config.branch,
        config.branch.solve_grid,
        dtype=dtype,
        device=config.backend.device,
    )
    final_state, solve = run_branch_continuation(
        config,
        grid,
        grid_name=f"solve_nr_{config.branch.solve_grid[0]}_nw_{config.branch.solve_grid[1]}",
    )
    boundaries = branch_boundary_conditions(config.branch)
    final_K = config.branch.continuation_K_values[-1]
    residual_fn = lambda x: coupled_branch_residual(
        x,
        grid,
        config.branch,
        eos_K=final_K,
        boundaries=boundaries,
    )
    jacobian_check = finite_difference_jvp_check(
        residual_fn,
        final_state,
        epsilon=config.branch.newton.finite_difference_jvp_epsilon,
        seed=config.jacobian_check_seed,
    )
    jacobian_passed = (
        jacobian_check["relative_residual"] <= config.jacobian_check_rel_tol
        or jacobian_check["absolute_residual"] <= config.jacobian_check_abs_tol
    )
    ladder = run_resolution_ladder(config)
    stop_flags = [
        {
            "item": "isotropic wall-to-branch static balance",
            "status": "STOP_AND_FLAG",
            "source_gap": (
                "compact lines 6783-6785 keep the wall PDE as an effective closure unless "
                "S_Sigma[R] is promoted; NONLINEAR_PROTOCOL_V2 lines 14-19 state that "
                "the full coupled physical residual equations are not frozen.  For WP1, "
                "eta=0, so this smoke prescribes R_0(w) and does not invent a wall force law."
            ),
        }
    ]
    results = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "coupled_mms": coupled_mms,
        "solve": solve,
        "jacobian_check": jacobian_check,
        "jacobian_passed": jacobian_passed,
        "ladder": ladder,
        "stop_flags": stop_flags,
        "passed": coupled_mms["passed"] and solve["converged"] and jacobian_passed,
    }
    return results


def _with_branch_preconditioner(
    branch: BranchSmokeConfig,
    preconditioner: PreconditionerConfig,
    *,
    ladder_levels: tuple[tuple[int, int], ...] | None = None,
) -> BranchSmokeConfig:
    return replace(
        branch,
        ladder_levels=ladder_levels if ladder_levels is not None else branch.ladder_levels,
        newton=replace(branch.newton, preconditioner=preconditioner),
    )


def _ladder_with_gmres_stats(ladder: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in ladder["rows"]:
        counts = row.get("gmres_iterations", [])
        rows.append(
            {
                **row,
                "gmres_max": max(counts) if counts else None,
                "gmres_mean": float(np.mean(counts)) if counts else None,
                "gmres_total": int(np.sum(counts)) if counts else 0,
            }
        )
    return rows


def _run_solution_comparison(
    unpreconditioned_config: HarnessConfig,
    preconditioned_config: HarnessConfig,
) -> dict[str, Any]:
    dtype = configure_backend(unpreconditioned_config.backend)
    grid = _create_branch_grid(
        unpreconditioned_config.branch,
        unpreconditioned_config.branch.solve_grid,
        dtype=dtype,
        device=unpreconditioned_config.backend.device,
    )
    grid_name = (
        f"comparison_nr_{unpreconditioned_config.branch.solve_grid[0]}"
        f"_nw_{unpreconditioned_config.branch.solve_grid[1]}"
    )
    unpreconditioned_state, unpreconditioned_summary = run_branch_continuation(
        unpreconditioned_config,
        grid,
        grid_name=f"{grid_name}_unpreconditioned",
    )
    preconditioned_state, preconditioned_summary = run_branch_continuation(
        preconditioned_config,
        grid,
        grid_name=f"{grid_name}_preconditioned",
    )
    difference = (preconditioned_state - unpreconditioned_state).detach()
    linf = float(torch.max(torch.abs(difference)).detach().cpu().item())
    l2 = float(torch.linalg.norm(difference).detach().cpu().item())
    tolerance = max(
        1.0e-6,
        10.0 * unpreconditioned_config.branch.newton.residual_atol,
        10.0 * preconditioned_config.branch.newton.residual_atol,
    )
    return {
        "grid": grid_name,
        "dof": int(unpreconditioned_state.numel()),
        "solution_linf_difference": linf,
        "solution_l2_difference": l2,
        "tolerance": tolerance,
        "passed": (
            unpreconditioned_summary["converged"]
            and preconditioned_summary["converged"]
            and linf <= tolerance
        ),
        "unpreconditioned": unpreconditioned_summary,
        "preconditioned": preconditioned_summary,
    }


def run_step3b_preconditioner(config: HarnessConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = HarnessConfig(
            run_root="software/stage1_solver/runs/step3b_preconditioner",
            report_path="software/stage1_solver/reports/step3b_preconditioner.md",
        )
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    shared_ladder = (
        (8, 6),
        (10, 8),
        (12, 10),
        (16, 12),
        (20, 14),
    )
    extended_ladder = shared_ladder + (
        (24, 16),
        (28, 20),
        (32, 24),
        (40, 28),
    )
    no_preconditioner = PreconditionerConfig(type="none")
    colored_lu_preconditioner = PreconditionerConfig(
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
    baseline_branch = _with_branch_preconditioner(
        config.branch,
        no_preconditioner,
        ladder_levels=shared_ladder,
    )
    preconditioned_branch = _with_branch_preconditioner(
        config.branch,
        colored_lu_preconditioner,
        ladder_levels=extended_ladder,
    )
    baseline_config = replace(
        config,
        run_root=str(Path(config.run_root) / "unpreconditioned"),
        branch=baseline_branch,
    )
    preconditioned_config = replace(
        config,
        run_root=str(Path(config.run_root) / "colored_sparse_jacobian_lu"),
        branch=preconditioned_branch,
    )

    coupled_mms = run_coupled_branch_mms(preconditioned_config)
    comparison = _run_solution_comparison(baseline_config, preconditioned_config)

    dtype = configure_backend(preconditioned_config.backend)
    grid = _create_branch_grid(
        preconditioned_config.branch,
        preconditioned_config.branch.solve_grid,
        dtype=dtype,
        device=preconditioned_config.backend.device,
    )
    boundaries = branch_boundary_conditions(preconditioned_config.branch)
    solve_state, solve = run_branch_continuation(
        preconditioned_config,
        grid,
        grid_name=(
            f"solve_nr_{preconditioned_config.branch.solve_grid[0]}"
            f"_nw_{preconditioned_config.branch.solve_grid[1]}"
        ),
    )
    final_K = preconditioned_config.branch.continuation_K_values[-1]
    residual_fn = lambda x: coupled_branch_residual(
        x,
        grid,
        preconditioned_config.branch,
        eos_K=final_K,
        boundaries=boundaries,
    )
    jacobian_check = finite_difference_jvp_check(
        residual_fn,
        solve_state,
        epsilon=preconditioned_config.branch.newton.finite_difference_jvp_epsilon,
        seed=preconditioned_config.jacobian_check_seed,
    )
    jacobian_passed = (
        jacobian_check["relative_residual"] <= preconditioned_config.jacobian_check_rel_tol
        or jacobian_check["absolute_residual"] <= preconditioned_config.jacobian_check_abs_tol
    )
    baseline_ladder = run_resolution_ladder(baseline_config)
    preconditioned_ladder = run_resolution_ladder(preconditioned_config)
    baseline_rows = _ladder_with_gmres_stats(baseline_ladder)
    preconditioned_rows = _ladder_with_gmres_stats(preconditioned_ladder)
    shared_dofs = {5 * nr * nw + 1 for nr, nw in shared_ladder}
    shared_preconditioned_rows = [
        row for row in preconditioned_rows if row["dof"] in shared_dofs
    ]
    baseline_growth = (
        baseline_rows[-1]["gmres_max"] / baseline_rows[0]["gmres_max"]
        if baseline_rows and baseline_rows[0]["gmres_max"]
        else None
    )
    preconditioned_growth = (
        shared_preconditioned_rows[-1]["gmres_max"]
        / shared_preconditioned_rows[0]["gmres_max"]
        if shared_preconditioned_rows and shared_preconditioned_rows[0]["gmres_max"]
        else None
    )
    max_preconditioned_dof = max(
        (row["dof"] for row in preconditioned_rows if row["converged"]),
        default=0,
    )
    scaling_passed = (
        preconditioned_growth is not None
        and baseline_growth is not None
        and preconditioned_growth < 0.5 * baseline_growth
        and max(
            row["gmres_max"] for row in shared_preconditioned_rows if row["gmres_max"] is not None
        )
        <= 0.5
        * max(row["gmres_max"] for row in baseline_rows if row["gmres_max"] is not None)
    )
    higher_resolution_passed = max_preconditioned_dof > max(
        row["dof"] for row in baseline_rows if row["converged"]
    )
    results = {
        "config": preconditioned_config.to_dict(),
        "config_hash": preconditioned_config.config_hash(),
        "baseline_config": baseline_config.to_dict(),
        "baseline_config_hash": baseline_config.config_hash(),
        "coupled_mms": coupled_mms,
        "solution_comparison": comparison,
        "solve": solve,
        "jacobian_check": jacobian_check,
        "jacobian_passed": jacobian_passed,
        "baseline_ladder": {
            **baseline_ladder,
            "rows": baseline_rows,
        },
        "preconditioned_ladder": {
            **preconditioned_ladder,
            "rows": preconditioned_rows,
        },
        "scaling": {
            "baseline_growth": baseline_growth,
            "preconditioned_growth_on_shared_ladder": preconditioned_growth,
            "max_preconditioned_dof": max_preconditioned_dof,
            "scaling_passed": scaling_passed,
            "higher_resolution_passed": higher_resolution_passed,
        },
        "passed": (
            coupled_mms["passed"]
            and comparison["passed"]
            and solve["converged"]
            and jacobian_passed
            and scaling_passed
            and higher_resolution_passed
        ),
    }
    return results


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.6e}"
    if isinstance(value, list):
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


def write_step3b_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    config = results["config"]
    branch = config["branch"]
    preconditioner = branch["newton"]["preconditioner"]
    mms = results["coupled_mms"]
    comparison = results["solution_comparison"]
    jac = results["jacobian_check"]
    solve = results["solve"]
    baseline_rows = results["baseline_ladder"]["rows"]
    preconditioned_rows = results["preconditioned_ladder"]["rows"]
    pre_by_dof = {row["dof"]: row for row in preconditioned_rows}
    shared_rows = []
    for base in baseline_rows:
        pre = pre_by_dof.get(base["dof"])
        if pre is None:
            continue
        shared_rows.append(
            {
                "grid": base["grid"],
                "dof": base["dof"],
                "baseline_gmres_max": base["gmres_max"],
                "preconditioned_gmres_max": pre["gmres_max"],
                "baseline_gmres_mean": base["gmres_mean"],
                "preconditioned_gmres_mean": pre["gmres_mean"],
                "baseline_seconds": base["wall_clock_seconds"],
                "preconditioned_seconds": pre["wall_clock_seconds"],
            }
        )
    scaling = results["scaling"]
    if scaling["scaling_passed"]:
        d2_read = (
            "The Krylov counts flatten when the linearized coupled operator is "
            "preconditioned by an assembled sparse Jacobian inverse.  This supports "
            "the design premise that the physics residual can remain in the torch-owned "
            "operator stack, but direct sparse LU is a CPU smoke method, not the "
            "production preconditioner.  The production path should replace this "
            "inverse with structured multigrid or a PETSc-class algebraic equivalent."
        )
    else:
        d2_read = (
            "The preconditioned counts did not flatten enough in this run.  That points "
            "toward the PETSc-class multigrid escape hatch rather than relying on this "
            "CPU preconditioner shape."
        )

    lines: list[str] = []
    lines.append("# Step 3b Coupled Branch JFNK Preconditioner")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Preconditioned config hash: `{results['config_hash']}`")
    lines.append(f"Baseline config hash: `{results['baseline_config_hash']}`")
    lines.append("")
    lines.append(
        "**Discipline:** engineering smoke, placeholder parameters, not a physical packet, "
        "target-blind. No field-to-coefficient export is performed."
    )
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append(
        "Interface: `solve_newton_jvp` accepts a left-preconditioner factory and passes the "
        "resulting SciPy `LinearOperator` as GMRES `M`.  The residual and JVP used by GMRES "
        "are unchanged."
    )
    lines.append(
        "Concrete preconditioner: assembled sparse Jacobian from the coupled residual JVP, "
        "colored by the radius-3 local stencil, factored with SciPy SuperLU, then used as an "
        "inverse preconditioner.  No new dependency was added."
    )
    lines.append("")
    lines.append("```yaml")
    for key, value in preconditioner.items():
        lines.append(f"{key}: {value}")
    lines.append("```")
    lines.append("")
    lines.append("## Correctness Preservation")
    lines.append("")
    lines.append(
        f"Fixed-grid solution comparison `{comparison['grid']}` / {comparison['dof']} DOF: "
        f"linf difference={comparison['solution_linf_difference']:.6e}, "
        f"l2 difference={comparison['solution_l2_difference']:.6e}, "
        f"tolerance={comparison['tolerance']:.6e}, "
        f"status={'PASS' if comparison['passed'] else 'FAIL'}."
    )
    lines.append(
        "Coupled residual JVP vs centered finite difference on the preconditioned solve: "
        f"relative={jac['relative_residual']:.6e}, absolute={jac['absolute_residual']:.6e}, "
        f"epsilon={jac['epsilon']:.6e}, "
        f"status={'PASS' if results['jacobian_passed'] else 'FAIL'}."
    )
    lines.append("")
    lines.append(
        _table(
            [
                "grid",
                "spacing",
                "error",
                "observed_order",
                "reference_norm",
                "spatial_gauge_l2",
                "spatial_current_l2",
            ],
            mms["rows"],
        )
    )
    lines.append("")
    lines.append("MMS checks:")
    for key, value in mms["pass_checks"].items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")
    lines.append("")
    lines.append("## Before After GMRES Curve")
    lines.append("")
    lines.append(
        _table(
            [
                "grid",
                "dof",
                "baseline_gmres_max",
                "preconditioned_gmres_max",
                "baseline_gmres_mean",
                "preconditioned_gmres_mean",
                "baseline_seconds",
                "preconditioned_seconds",
            ],
            shared_rows,
        )
    )
    lines.append("")
    lines.append(
        f"Baseline max-GMRES growth on the shared ladder: "
        f"{scaling['baseline_growth']:.6e}. Preconditioned growth on the same ladder: "
        f"{scaling['preconditioned_growth_on_shared_ladder']:.6e}."
    )
    lines.append("")
    lines.append("## Extended Preconditioned Ladder")
    lines.append("")
    lines.append(
        _table(
            [
                "grid",
                "dof",
                "wall_clock_seconds",
                "peak_memory_mb",
                "newton_iterations",
                "final_residual_linf",
                "gmres_iterations",
                "gmres_max",
                "gmres_mean",
                "converged",
                "message",
            ],
            preconditioned_rows,
        )
    )
    lines.append("")
    lines.append(
        f"New maximum converged DOF on this laptop run: {scaling['max_preconditioned_dof']}. "
        f"Ladder stop reason: {results['preconditioned_ladder']['stop_reason']}."
    )
    lines.append("")
    lines.append("## Main Preconditioned Solve")
    lines.append("")
    lines.append(
        f"`{solve['grid']}` final residual linf={solve['final_residual_linf']:.6e}, "
        f"wall-clock={solve['wall_clock_seconds']:.6e}s, "
        f"peak RSS={solve['peak_memory_mb']:.6e} MB, manifest=`{solve['manifest']}`."
    )
    lines.append("")
    lines.append("## Scaling Trigger Read")
    lines.append("")
    lines.append(d2_read)
    lines.append("")
    lines.append("## Manifests")
    lines.append("")
    for row in mms["rows"]:
        lines.append(f"- coupled MMS {row['grid']}: `{row['manifest']}`")
    lines.append(f"- comparison unpreconditioned: `{comparison['unpreconditioned']['manifest']}`")
    lines.append(f"- comparison preconditioned: `{comparison['preconditioned']['manifest']}`")
    lines.append(f"- main preconditioned solve: `{solve['manifest']}`")
    for row in baseline_rows:
        lines.append(f"- baseline ladder {row['grid']}: `{row['manifest']}`")
    for row in preconditioned_rows:
        lines.append(f"- preconditioned ladder {row['grid']}: `{row['manifest']}`")
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def write_step3_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    config = results["config"]
    branch = config["branch"]
    coupled = config["mms"]["coupled_branch"]
    solve = results["solve"]
    mms = results["coupled_mms"]
    jac = results["jacobian_check"]
    lines: list[str] = []
    lines.append("# Step 3 Coupled Branch Engineering Smoke")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append("")
    lines.append(
        "**Discipline:** engineering smoke, placeholder parameters, not a physical packet, target-blind. "
        "No field-to-coefficient export is performed."
    )
    lines.append("")
    lines.append("## Coupled Residual")
    lines.append("")
    lines.append(
        "Matter sector: stationary gauged quintic GNLS on `(r,w)`, "
        "`[-hbar^2/(2m) D_iD_i + V_conf + (5K/4)rho^4 + q_star A0 - mu] psi`; "
        "`D_i = d_i - i(q_star/hbar) A_i`, expanded with the axisymmetric full radial measure."
    )
    lines.append(
        "Maxwell sector: localized H=Z operator over `(A0, Ar, Aw)` with "
        "`d_M(Z F^{MN}) + xi^-1 d^N(Z d.A) = mu0 J_psi^N`; "
        "`J_psi^0=q_star rho`, `J_psi^i=q_star j^i`, and "
        "`j^i=(hbar/m) Im(psi* D_i psi)`.  No external matter source is added."
    )
    lines.append(
        "Sources: compact lines 556-583, 638-648, 651-659, 674-689; "
        "prereg lines 44-56; NONLINEAR_PROTOCOL_V2 lines 35-38."
    )
    lines.append("")
    lines.append("## Placeholder Parameters")
    lines.append("")
    lines.append(f"Label: {branch['placeholder_label']}")
    lines.append("")
    lines.append("```yaml")
    for key in (
        "hbar",
        "particle_mass",
        "gauge_charge",
        "mu0",
        "xi",
        "continuation_K_values",
        "mass",
        "localization_width",
        "localization_floor",
        "localization_amplitude",
        "r_mouth",
        "r_exit",
        "geometry_decay_length",
        "radial_wall_strength",
        "axial_trap_strength",
        "solve_grid",
        "ladder_levels",
    ):
        lines.append(f"{key}: {branch[key]}")
    lines.append("```")
    lines.append("")
    lines.append(
        "`R_0(w)=R_exit+(R_mouth-R_exit) exp(-(w-w_min)/geometry_decay_length)` is prescribed "
        "for this isotropic smoke; it satisfies the mouth value and finite open-exit Robin class."
    )
    lines.append("")
    lines.append("## Coupled MMS")
    lines.append("")
    lines.append(mms["description"])
    lines.append(f"Forcing: {mms['forcing_derivation']}")
    lines.append("MMS placeholders:")
    lines.append("```yaml")
    for key in (
        "hbar",
        "particle_mass",
        "gauge_charge",
        "mu0",
        "xi",
        "eos_K",
        "chemical_potential",
        "localization_width",
        "localization_floor",
        "localization_amplitude",
        "r_mouth",
        "r_exit",
        "geometry_decay_length",
        "radial_wall_strength",
        "axial_trap_strength",
    ):
        lines.append(f"{key}: {coupled[key]}")
    lines.append("```")
    lines.append("")
    lines.append(
        _table(
            [
                "grid",
                "spacing",
                "error",
                "observed_order",
                "reference_norm",
                "spatial_gauge_l2",
                "spatial_current_l2",
            ],
            mms["rows"],
        )
    )
    lines.append("")
    lines.append("MMS checks:")
    for key, value in mms["pass_checks"].items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")
    lines.append("")
    lines.append("## Newton And Continuation")
    lines.append("")
    stage_rows = []
    for stage in solve["stages"]:
        stage_rows.append(
            {
                "eos_K": stage["eos_K"],
                "converged": stage["converged"],
                "iterations": stage["iterations"],
                "initial_residual_norm": stage["initial_residual_norm"],
                "final_residual_norm": stage["final_residual_norm"],
                "gmres_iterations": stage["gmres_iterations"],
                "message": stage["message"],
            }
        )
    lines.append(_table(
        [
            "eos_K",
            "converged",
            "iterations",
            "initial_residual_norm",
            "final_residual_norm",
            "gmres_iterations",
            "message",
        ],
        stage_rows,
    ))
    lines.append("")
    lines.append(
        f"Main solve `{solve['grid']}`: final residual linf={solve['final_residual_linf']:.6e}, "
        f"wall-clock={solve['wall_clock_seconds']:.6e}s, peak RSS={solve['peak_memory_mb']:.6e} MB, "
        f"manifest=`{solve['manifest']}`."
    )
    lines.append("")
    lines.append("## Jacobian Check")
    lines.append("")
    lines.append(
        "Coupled residual JVP vs centered finite difference: "
        f"relative={jac['relative_residual']:.6e}, absolute={jac['absolute_residual']:.6e}, "
        f"epsilon={jac['epsilon']:.6e}, status={'PASS' if results['jacobian_passed'] else 'FAIL'}."
    )
    lines.append("")
    lines.append("## Resolution Ladder")
    lines.append("")
    lines.append(
        _table(
            [
                "grid",
                "dof",
                "wall_clock_seconds",
                "peak_memory_mb",
                "newton_iterations",
                "final_residual_linf",
                "gmres_iterations",
                "converged",
                "message",
            ],
            results["ladder"]["rows"],
        )
    )
    lines.append("")
    lines.append(f"Ladder stop reason: {results['ladder']['stop_reason']}")
    lines.append("")
    lines.append("## STOP And Flag")
    lines.append("")
    for item in results["stop_flags"]:
        lines.append(f"- {item['item']}: {item['status']}. {item['source_gap']}")
    lines.append("")
    lines.append("## Manifests")
    lines.append("")
    for row in mms["rows"]:
        lines.append(f"- coupled MMS {row['grid']}: `{row['manifest']}`")
    lines.append(f"- main continuation: `{solve['manifest']}`")
    for row in results["ladder"]["rows"]:
        lines.append(f"- ladder {row['grid']}: `{row['manifest']}`")
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path
