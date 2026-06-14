"""Step 8b driven P2 tangent operator with a complex absorbing exit layer.

The background-dependent piece remains the Step-8a JVP.  This module adds only
the compact §4.7 temporal terms and an additive CAP, then solves the resulting
linear complex system directly.  The physical matter/gauge-to-wall source
``S_eta^(psi,A)`` is still absent; all drives are target-blind surrogate
``f_ext`` data.
"""

from __future__ import annotations

from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import csc_matrix, coo_matrix
from scipy.sparse.linalg import splu
import sympy as sp
import torch

from .backend import configure_backend
from .boundaries import BoundaryCondition
from .config import BranchSmokeConfig, HarnessConfig, TensorGridSpec
from .coupled_branch import (
    _create_branch_grid,
    branch_boundary_conditions,
    localization_weight_torch,
    resample_branch_state,
    run_branch_continuation,
)
from .grid import TensorProductGrid, WallGrid
from .manifest import write_manifest
from .mms import run_convergence_study
from .mms_benchmarks import (
    _P2_CENTRIFUGAL,
    _P2_MAXWELL_ANGULAR,
    _WALL,
    _torch_from_numpy,
)
from .newton import PreconditionerBuildContext
from .operators import (
    localized_maxwell_operator_with_spherical_l,
    tensor_laplacian_with_spherical_l,
    wall_s_eta_operator,
)
from .p2_tangent import (
    SOURCE_CITATIONS,
    P2TangentFields,
    _background_summary,
    _json_default,
    _strip_runtime_level,
    _summarize_p2_observables,
    _table,
    apply_p2_tangent,
    p2_surrogate_forcing,
    p2_tangent_fd_check,
    p2_wall_coefficients,
    pack_p2_tangent_fields,
    unpack_p2_tangent_fields,
    with_step8a_preconditioners,
    zero_p2_tangent_state,
)
from .preconditioners import assemble_p2_tangent_colored_sparse_jacobian


DRIVEN_SOURCE_CITATIONS = {
    **SOURCE_CITATIONS,
    "time_derivative_terms": (
        "notes/moving_throat_pde_program_compact.md lines 1383-1455: "
        "matter i*hbar*d_t doublet, boxed Maxwell equation, boxed wall equation."
    ),
    "absorber_reference": (
        "docs/boundary_and_noise_methods.md sections 0, 3, 5: Step-8 CAP decision "
        "and a generic CAP in the spirit of Manolopoulos 2002; this implementation "
        "does not use the wavelength-sized pole profile and claims no analytic "
        "reflection-free guarantee."
    ),
    "maxwell_temporal_truncation": (
        "The driven Maxwell temporal block is an engineering-smoke truncation: only "
        "the diagonal -Z(w)*omega^2*A_N self-term is retained.  The omega-linear "
        "A0<->Ai temporal cross-couplings and gauge-term temporal pieces are dropped; "
        "the A0 coefficient is a placeholder, not the correct N=0 gauge-term reduction."
    ),
}

DRIVEN_OBSERVABLE_LABELS = {
    "total_response_l2": "total raw complex response L2",
    "matter_interior_response_l2": "matter perturbation interior L2",
    "scalar_gauge_response_l2": "A0 perturbation L2",
    "spatial_gauge_response_l2": "spatial gauge perturbation L2",
    "wall_eta_l2": "wall eta L2",
    "wall_lower_endpoint_eta_abs": "lower-wall eta absolute value",
    "driven_residual_linf": "driven tangent residual Linf",
}

_RATIONAL_CAP_C = 2.62206
MAXWELL_TEMPORAL_TRUNCATION_NOTE = (
    "Retains only the diagonal -Z(w)*omega^2*A_N self-term.  It drops the omega-linear "
    "A0<->Ai temporal cross-couplings and the gauge-term temporal pieces; the A0 lane's "
    "-Z*omega^2 coefficient is an engineering-smoke placeholder, not the correct N=0 "
    "gauge-term reduction.  This temporal truncation is separate from the Step-8a "
    "spatial Ar scalarization and separate from the open S_eta^(psi,A) source."
)
GENERIC_CAP_NOTE = (
    "Generic smooth additive -i*sigma(w) upper-exit CAP using a normalized polynomial/"
    "rational ramp inspired by the Manolopoulos rational form, but not the wavelength-"
    "sized pole profile, not the Manolopoulos transmission-free profile, and not "
    "assigned an analytic reflection bound; reflection is only the measured "
    "domain-extension contamination metric."
)


def step8b_default_config(
    *,
    run_root: str = "software/stage1_solver/runs/step8b_driven_absorber",
    report_path: str = "software/stage1_solver/reports/step8b_driven_absorber.md",
) -> HarnessConfig:
    """Bounded CPU config for the Step-8b acceptance harness."""

    base = HarnessConfig(run_root=run_root, report_path=report_path)
    branch = replace(
        base.branch,
        solve_grid=(4, 4),
        ladder_levels=((4, 4),),
        continuation_K_values=(0.05,),
        newton=replace(
            base.branch.newton,
            residual_atol=1.0e-8,
            residual_rtol=1.0e-8,
            gmres_rtol=1.0e-8,
            gmres_atol=1.0e-10,
            gmres_restart=96,
            gmres_maxiter=8,
            max_newton_iters=10,
        ),
    )
    p2_tangent = replace(
        base.p2_tangent,
        convergence_levels=((4, 4), (8, 8), (16, 16)),
        fd_check_grid=(4, 4),
        wellposed_grid=(4, 4),
        refinement_ratio=2,
        min_observable_order=1.45,
        newton=replace(
            base.p2_tangent.newton,
            residual_atol=1.0e-9,
            residual_rtol=1.0e-9,
            gmres_rtol=1.0e-9,
            gmres_atol=1.0e-11,
            gmres_restart=128,
            gmres_maxiter=8,
        ),
    )
    p2_driven = replace(
        base.p2_driven,
        convergence_levels=((4, 4), (8, 8), (16, 16)),
        reflection_grid=(6, 6),
        response_table_grid=(6, 6),
        reflection_fit_window=(0.25, 0.75),
        min_observable_order=1.45,
    )
    mms = replace(
        base.mms,
        p2_driven=replace(
            base.mms.p2_driven,
            tensor_grid_levels=((24, 20), (48, 40), (96, 80)),
            wall_grid_levels=(32, 64, 128),
            tensor_final_error_max=1.0e-2,
            wall_final_error_max=2.0e-3,
        ),
    )
    return replace(base, branch=branch, p2_tangent=p2_tangent, p2_driven=p2_driven, mms=mms)


def _smooth_polynomial_rational_shape_numpy(x: np.ndarray) -> np.ndarray:
    c = _RATIONAL_CAP_C
    a = 1.0 - 16.0 / c**3
    b = (1.0 - 17.0 / c**3) / c**2
    clipped = np.clip(x, 0.0, 1.0)
    raw = a * clipped - b * clipped**3 + 4.0 / (c - clipped) ** 2 - 4.0 / (c + clipped) ** 2
    edge = a - b + 4.0 / (c - 1.0) ** 2 - 4.0 / (c + 1.0) ** 2
    return raw / edge


def _smooth_polynomial_rational_shape_torch(x: torch.Tensor) -> torch.Tensor:
    c = _RATIONAL_CAP_C
    a = 1.0 - 16.0 / c**3
    b = (1.0 - 17.0 / c**3) / c**2
    clipped = torch.clamp(x, min=0.0, max=1.0)
    raw = a * clipped - b * clipped**3 + 4.0 / (c - clipped) ** 2 - 4.0 / (c + clipped) ** 2
    edge = a - b + 4.0 / (c - 1.0) ** 2 - 4.0 / (c + 1.0) ** 2
    return raw / edge


def p2_cap_profile_torch(
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    enabled: bool | None = None,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Return the generic smooth additive CAP strengths for the upper ``w`` exit."""

    driven = config.p2_driven
    use_cap = driven.cap_enabled if enabled is None else enabled
    if not use_cap:
        cell = torch.zeros((grid.spec.nr, grid.spec.nw), dtype=grid.dtype, device=grid.device)
        wall = torch.zeros(grid.spec.nw, dtype=grid.dtype, device=grid.device)
        return cell, wall
    width = float(driven.cap_width)
    if width <= 0.0 or width >= (grid.spec.w_max - grid.spec.w_min):
        raise ValueError("CAP width must be positive and smaller than the axial domain")
    x = (grid.w_centers - (grid.spec.w_max - width)) / width
    wall_sigma = driven.cap_strength * _smooth_polynomial_rational_shape_torch(x)
    cell_sigma = wall_sigma[None, :].expand(grid.spec.nr, grid.spec.nw)
    return cell_sigma, wall_sigma


def p2_cap_profile_numpy(
    w_centers: np.ndarray,
    *,
    w_max: float,
    width: float,
    strength: float,
    enabled: bool,
) -> np.ndarray:
    if not enabled:
        return np.zeros_like(w_centers, dtype=np.float64)
    x = (w_centers - (w_max - width)) / width
    return strength * _smooth_polynomial_rational_shape_numpy(x)


def p2_driven_surrogate_forcing(
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
) -> torch.Tensor:
    """Target-blind complex surrogate drive extending the Step-8a forcing."""

    base = p2_surrogate_forcing(grid, config.p2_tangent).to(torch.complex128)
    fields = unpack_p2_tangent_fields(base, grid)
    length = grid.spec.w_max - grid.spec.w_min
    x = (grid.w_centers - grid.spec.w_min) / length
    source_center = grid.spec.w_min + 0.28 * length
    source_width = 0.16 * length
    source = torch.exp(-((grid.w_centers - source_center) / source_width) ** 2).to(torch.complex128)
    phase = torch.exp(1j * torch.as_tensor(omega, dtype=torch.float64) * (0.15 + x)).to(
        torch.complex128
    )
    eta = fields.eta + 0.8 * config.p2_tangent.surrogate_force_amplitude * source * phase
    r_profile = (grid.r_centers[:, None] / grid.spec.r_max) ** 2
    axial = (source * phase)[None, :]
    driven_fields = P2TangentFields(
        psi_real=fields.psi_real + 0.25 * config.p2_tangent.surrogate_force_amplitude * r_profile * axial,
        psi_imag=fields.psi_imag - 0.20j * config.p2_tangent.surrogate_force_amplitude * r_profile * axial,
        a0=fields.a0,
        ar=fields.ar,
        aw=fields.aw,
        eta=eta,
    )
    return pack_p2_tangent_fields(driven_fields)


def _complexified_static_tangent(
    direction: torch.Tensor,
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    p2_cfg: Any,
    *,
    eos_K: float,
    boundaries: Any,
) -> torch.Tensor:
    real_part = apply_p2_tangent(
        torch.real(direction).to(grid.dtype),
        background_state,
        grid,
        cfg,
        p2_cfg,
        eos_K=eos_K,
        boundaries=boundaries,
    ).to(torch.complex128)
    imag_values = torch.imag(direction)
    if float(torch.linalg.norm(imag_values).detach().cpu().item()) == 0.0:
        return real_part
    imag_part = apply_p2_tangent(
        imag_values.to(grid.dtype),
        background_state,
        grid,
        cfg,
        p2_cfg,
        eos_K=eos_K,
        boundaries=boundaries,
    ).to(torch.complex128)
    return real_part + 1j * imag_part


def p2_driven_frequency_terms(
    state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
    cap_enabled: bool | None = None,
) -> torch.Tensor:
    """Apply retained driven-frequency terms and the additive CAP.

    Maxwell caveat: the temporal block below is only the diagonal
    ``-Z(w)*omega^2*A_N`` self-term.  It omits the omega-linear A0<->Ai temporal
    cross-couplings and the gauge-term temporal pieces; the A0 lane coefficient
    is a placeholder rather than the correct N=0 gauge-term reduction.
    """

    fields = unpack_p2_tangent_fields(state.to(torch.complex128), grid)
    branch = config.branch
    p2_cfg = config.p2_tangent
    weight = localization_weight_torch(grid.w_centers, branch).to(torch.complex128)
    cell_sigma, wall_sigma = p2_cap_profile_torch(grid, config, enabled=cap_enabled)
    cell_sigma_c = cell_sigma.to(torch.complex128)
    wall_sigma_c = wall_sigma.to(torch.complex128)
    omega_value = torch.as_tensor(float(omega), dtype=torch.float64, device=state.device)

    matter_re = 1j * branch.hbar * omega_value * fields.psi_imag
    matter_im = -1j * branch.hbar * omega_value * fields.psi_real
    maxwell_factor = -weight[None, :] * omega_value**2
    wall_factor = -p2_cfg.placeholder_mu_eta * omega_value**2

    terms = P2TangentFields(
        psi_real=matter_re - 1j * cell_sigma_c * fields.psi_real,
        psi_imag=matter_im - 1j * cell_sigma_c * fields.psi_imag,
        a0=maxwell_factor * fields.a0 - 1j * cell_sigma_c * fields.a0,
        ar=maxwell_factor * fields.ar - 1j * cell_sigma_c * fields.ar,
        aw=maxwell_factor * fields.aw - 1j * cell_sigma_c * fields.aw,
        eta=wall_factor * fields.eta - 1j * wall_sigma_c * fields.eta,
    )
    return pack_p2_tangent_fields(terms)


def apply_p2_driven_tangent(
    direction: torch.Tensor,
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
    cap_enabled: bool | None = None,
) -> torch.Tensor:
    boundaries = branch_boundary_conditions(config.branch)
    static = _complexified_static_tangent(
        direction.to(torch.complex128),
        background_state,
        grid,
        config.branch,
        config.p2_tangent,
        eos_K=config.branch.continuation_K_values[-1],
        boundaries=boundaries,
    )
    return static + p2_driven_frequency_terms(
        direction.to(torch.complex128),
        grid,
        config,
        omega=omega,
        cap_enabled=cap_enabled,
    )


def p2_driven_linear_residual(
    state: torch.Tensor,
    forcing: torch.Tensor,
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
    cap_enabled: bool | None = None,
) -> torch.Tensor:
    return apply_p2_driven_tangent(
        state,
        background_state,
        grid,
        config,
        omega=omega,
        cap_enabled=cap_enabled,
    ) - forcing.to(torch.complex128)


def p2_driven_fd_check(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
    epsilon: float,
    seed: int,
) -> dict[str, float]:
    generator = torch.Generator(device=grid.device)
    generator.manual_seed(seed)
    shape = zero_p2_tangent_state(grid).shape
    real = torch.randn(shape, dtype=grid.dtype, device=grid.device, generator=generator)
    imag = torch.randn(shape, dtype=grid.dtype, device=grid.device, generator=generator)
    direction = (real + 1j * imag).to(torch.complex128)
    direction = direction / torch.linalg.norm(direction)
    zero = torch.zeros(shape, dtype=torch.complex128, device=grid.device)
    forcing = torch.zeros_like(zero)
    residual_fn = lambda mode: p2_driven_linear_residual(
        mode,
        forcing,
        background_state,
        grid,
        config,
        omega=omega,
        cap_enabled=True,
    )
    applied = apply_p2_driven_tangent(
        direction,
        background_state,
        grid,
        config,
        omega=omega,
        cap_enabled=True,
    )
    fd = (residual_fn(zero + epsilon * direction) - residual_fn(zero - epsilon * direction)) / (
        2.0 * epsilon
    )
    diff = applied - fd
    absolute = float(torch.linalg.norm(diff).detach().cpu().item())
    applied_norm = float(torch.linalg.norm(applied).detach().cpu().item())
    return {
        "epsilon": float(epsilon),
        "absolute_residual": absolute,
        "relative_residual": absolute / max(1.0, applied_norm),
        "applied_norm": applied_norm,
        "fd_norm": float(torch.linalg.norm(fd).detach().cpu().item()),
    }


def _assemble_static_sparse_matrix(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
) -> tuple[csc_matrix, dict[str, Any]]:
    boundaries = branch_boundary_conditions(config.branch)
    final_K = config.branch.continuation_K_values[-1]
    zero = zero_p2_tangent_state(grid)
    residual_fn = lambda x: apply_p2_tangent(
        x,
        background_state,
        grid,
        config.branch,
        config.p2_tangent,
        eos_K=final_K,
        boundaries=boundaries,
    )
    matrix, metadata = assemble_p2_tangent_colored_sparse_jacobian(
        PreconditionerBuildContext(
            residual_fn=residual_fn,
            x=zero,
            rhs=np.zeros(zero.numel(), dtype=np.float64),
            iteration=1,
            config=config.p2_tangent.newton,
        ),
        grid,
    )
    return matrix, metadata


def _frequency_sparse_delta(
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
    cap_enabled: bool,
) -> csc_matrix:
    n = grid.spec.nr * grid.spec.nw
    size = 5 * n + grid.spec.nw
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []

    def add_diag(offset: int, values: np.ndarray) -> None:
        idx = offset + np.arange(values.size, dtype=np.int64)
        rows.append(idx)
        cols.append(idx)
        data.append(values.astype(np.complex128, copy=False))

    def add_block(row_offset: int, col_offset: int, values: np.ndarray) -> None:
        idx = np.arange(values.size, dtype=np.int64)
        rows.append(row_offset + idx)
        cols.append(col_offset + idx)
        data.append(values.astype(np.complex128, copy=False))

    weight = localization_weight_torch(grid.w_centers, config.branch).detach().cpu().numpy()
    z_flat = np.tile(weight[None, :], (grid.spec.nr, 1)).reshape(-1)
    cell_sigma, wall_sigma = p2_cap_profile_torch(grid, config, enabled=cap_enabled)
    sigma_flat = cell_sigma.detach().cpu().numpy().reshape(-1)
    wall_sigma_np = wall_sigma.detach().cpu().numpy()
    omega2 = float(omega) ** 2
    matter_value = np.full(n, config.branch.hbar * float(omega), dtype=np.float64)
    add_block(0, n, 1j * matter_value)
    add_block(n, 0, -1j * matter_value)
    for lane in range(5):
        add_diag(lane * n, -1j * sigma_flat)
    for lane in (2, 3, 4):
        add_diag(lane * n, -z_flat * omega2)
    wall_offset = 5 * n
    add_diag(
        wall_offset,
        -config.p2_tangent.placeholder_mu_eta * omega2 * np.ones(grid.spec.nw)
        - 1j * wall_sigma_np,
    )
    return coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(size, size),
        dtype=np.complex128,
    ).tocsc()


def p2_driven_sparse_matrix(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
    cap_enabled: bool = True,
) -> tuple[csc_matrix, dict[str, Any]]:
    static, metadata = _assemble_static_sparse_matrix(background_state, grid, config)
    matrix = static.astype(np.complex128) + _frequency_sparse_delta(
        grid,
        config,
        omega=omega,
        cap_enabled=cap_enabled,
    )
    metadata = dict(metadata)
    metadata.update(
        {
            "complex_direct_solve": True,
            "omega": float(omega),
            "cap_enabled": bool(cap_enabled),
            "cap_width": config.p2_driven.cap_width if cap_enabled else 0.0,
            "cap_strength": config.p2_driven.cap_strength if cap_enabled else 0.0,
            "cap_profile": config.p2_driven.cap_profile if cap_enabled else "off",
            "cap_profile_note": GENERIC_CAP_NOTE if cap_enabled else "CAP disabled",
            "maxwell_temporal_truncation": MAXWELL_TEMPORAL_TRUNCATION_NOTE,
            "matrix_nnz_complex": int(matrix.nnz),
        }
    )
    return matrix, metadata


def p2_driven_response_observables(
    state: torch.Tensor,
    grid: TensorProductGrid,
    *,
    residual_linf: float | None = None,
) -> dict[str, float]:
    fields = unpack_p2_tangent_fields(state.to(torch.complex128), grid)
    volumes = grid.cell_volumes.to(torch.float64)
    wall_grid = WallGrid.create(
        config_wall_grid_spec(grid),
        dtype=grid.dtype,
        device=grid.device,
    )

    def norm2(values: torch.Tensor, weights: torch.Tensor) -> torch.Tensor:
        return torch.sqrt(torch.sum(torch.real(torch.conj(values) * values) * weights))

    matter_density = torch.real(torch.conj(fields.psi_real) * fields.psi_real) + torch.real(
        torch.conj(fields.psi_imag) * fields.psi_imag
    )
    interior = (
        (grid.r_centers[:, None] < 0.75 * grid.spec.r_max)
        & (
            grid.w_centers[None, :]
            < grid.spec.w_min + 0.75 * (grid.spec.w_max - grid.spec.w_min)
        )
    )
    matter = torch.sqrt(torch.sum(matter_density * volumes))
    matter_interior = torch.sqrt(torch.sum(matter_density[interior] * volumes[interior]))
    scalar = norm2(fields.a0, volumes)
    spatial = torch.sqrt(
        torch.sum(
            (
                torch.real(torch.conj(fields.ar) * fields.ar)
                + torch.real(torch.conj(fields.aw) * fields.aw)
            )
            * volumes
        )
    )
    wall = norm2(fields.eta, wall_grid.cell_widths)
    total = torch.sqrt(matter**2 + scalar**2 + spatial**2 + wall**2)
    result = {
        "total_response_l2": float(total.detach().cpu().item()),
        "matter_interior_response_l2": float(matter_interior.detach().cpu().item()),
        "scalar_gauge_response_l2": float(scalar.detach().cpu().item()),
        "spatial_gauge_response_l2": float(spatial.detach().cpu().item()),
        "wall_eta_l2": float(wall.detach().cpu().item()),
        "wall_lower_endpoint_eta_abs": float(torch.abs(fields.eta[0]).detach().cpu().item()),
    }
    if residual_linf is not None:
        result["driven_residual_linf"] = float(residual_linf)
    return result


def config_wall_grid_spec(grid: TensorProductGrid):
    from .config import WallGridSpec

    return WallGridSpec(w_min=grid.spec.w_min, w_max=grid.spec.w_max, nw=grid.spec.nw)


def solve_driven_p2_tangent(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
    grid_name: str,
    cap_enabled: bool = True,
) -> tuple[torch.Tensor, dict[str, Any]]:
    matrix, matrix_metadata = p2_driven_sparse_matrix(
        background_state,
        grid,
        config,
        omega=omega,
        cap_enabled=cap_enabled,
    )
    forcing = p2_driven_surrogate_forcing(grid, config, omega=omega)
    forcing_np = forcing.detach().cpu().numpy().astype(np.complex128, copy=False)
    factor = splu(matrix, permc_spec=config.p2_tangent.newton.preconditioner.permutation)
    solution_np = factor.solve(forcing_np)
    residual_np = matrix @ solution_np - forcing_np
    state = torch.as_tensor(solution_np, dtype=torch.complex128, device=grid.device)
    residual_linf = float(np.max(np.abs(residual_np)))
    residual_l2 = float(np.linalg.norm(residual_np))
    summary = {
        "grid": grid_name,
        "nr": grid.spec.nr,
        "nw": grid.spec.nw,
        "dof": int(state.numel()),
        "omega": float(omega),
        "cap_enabled": bool(cap_enabled),
        "cap_profile": config.p2_driven.cap_profile if cap_enabled else "off",
        "cap_width": config.p2_driven.cap_width if cap_enabled else 0.0,
        "cap_strength": config.p2_driven.cap_strength if cap_enabled else 0.0,
        "converged": residual_linf <= config.p2_tangent.newton.residual_atol * 100.0,
        "iterations": 1,
        "initial_residual_norm": float(np.linalg.norm(forcing_np)),
        "final_residual_norm": residual_l2,
        "driven_residual_linf": residual_linf,
        "tolerance": config.p2_tangent.newton.residual_atol,
        "message": "complex sparse direct solve",
        "newton_history": [],
        "preconditioner": asdict(config.p2_tangent.newton.preconditioner),
        "matrix": {
            **matrix_metadata,
            "factor_nnz_l": int(factor.L.nnz),
            "factor_nnz_u": int(factor.U.nnz),
        },
        "forcing_l2": float(np.linalg.norm(forcing_np)),
        "surrogate_values": p2_driven_response_observables(
            state,
            grid,
            residual_linf=residual_linf,
        ),
    }
    manifest = write_manifest(
        run_root=config.run_root,
        benchmark_name=config.p2_driven.name,
        grid_name=grid_name,
        config=config.to_dict(),
        mesh=grid.to_dict(),
        results=summary,
        config_hash=config.config_hash(),
        solver_controls={
            "linear_solver": "scipy.sparse.linalg.splu",
            "omega": float(omega),
            "cap_enabled": bool(cap_enabled),
            "cap_profile": summary["cap_profile"],
            "cap_profile_note": GENERIC_CAP_NOTE if cap_enabled else "CAP disabled",
            "maxwell_temporal_truncation": MAXWELL_TEMPORAL_TRUNCATION_NOTE,
        },
    )
    summary["manifest"] = str(manifest)
    return state.detach(), summary


def run_p2_driven_ladder(config: HarnessConfig, *, omega: float) -> list[dict[str, Any]]:
    dtype = configure_backend(config.backend)
    levels: list[dict[str, Any]] = []
    previous_grid: TensorProductGrid | None = None
    previous_background: torch.Tensor | None = None
    for level_index, level in enumerate(config.p2_driven.convergence_levels):
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
            grid_name=f"p2_driven_background_l{level_index}_nr_{level[0]}_nw_{level[1]}",
        )
        driven_state, solve = solve_driven_p2_tangent(
            background_state,
            grid,
            config,
            omega=omega,
            grid_name=f"p2_driven_l{level_index}_nr_{level[0]}_nw_{level[1]}",
            cap_enabled=True,
        )
        levels.append(
            {
                "level": level_index,
                "grid": solve["grid"],
                "nr": grid.spec.nr,
                "nw": grid.spec.nw,
                "spacing": max(grid.dr, grid.dw),
                "dof": solve["dof"],
                "omega": float(omega),
                "background": _background_summary(background),
                "solve": solve,
                "converged": background["converged"] and solve["converged"],
                "surrogate_values": solve["surrogate_values"],
                "grid_object": grid,
                "background_state": background_state.detach(),
                "driven_state": driven_state.detach(),
            }
        )
        previous_grid = grid
        previous_background = background_state.detach()
    return levels


def p2_driven_wellposedness_check(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
) -> dict[str, float]:
    matrix, _metadata = p2_driven_sparse_matrix(
        background_state,
        grid,
        config,
        omega=omega,
        cap_enabled=True,
    )
    dense = matrix.toarray()
    singular_values = np.linalg.svd(dense, compute_uv=False)
    return {
        "smallest_singular_value": float(np.min(singular_values)),
        "largest_singular_value": float(np.max(singular_values)),
        "condition_number": float(np.max(singular_values) / np.min(singular_values)),
        "state_size": int(dense.shape[0]),
        "omega": float(omega),
    }


def p2_omega_zero_static_limit_check(
    background_state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
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
    boundaries = branch_boundary_conditions(config.branch)
    static = apply_p2_tangent(
        direction,
        background_state,
        grid,
        config.branch,
        config.p2_tangent,
        eos_K=config.branch.continuation_K_values[-1],
        boundaries=boundaries,
    ).to(torch.complex128)
    driven = apply_p2_driven_tangent(
        direction.to(torch.complex128),
        background_state,
        grid,
        config,
        omega=0.0,
        cap_enabled=False,
    )
    diff = driven - static
    absolute = float(torch.linalg.norm(diff).detach().cpu().item())
    static_norm = float(torch.linalg.norm(static).detach().cpu().item())
    return {
        "absolute_residual": absolute,
        "relative_residual": absolute / max(1.0, static_norm),
        "static_norm": static_norm,
    }


def p2_matter_free_particle_dispersion_check(config: HarnessConfig) -> dict[str, float | str | bool]:
    """Check that the corrected matter time block gives hbar*omega=+E(k).

    For a free kinetic carrier with energy ``E(k)>0``, the positive branch
    lane vector ``(lane0,lane1)=(1,i)`` should satisfy
    ``E I + hbar*omega [[0,i],[-i,0]] = 0`` at ``hbar*omega=E``.
    The old negated sign gives ``2E`` on the same branch.
    """

    dtype = configure_backend(config.backend)
    grid = TensorProductGrid.create(
        TensorGridSpec(
            r_max=config.branch.r_max,
            nr=4,
            w_min=config.branch.w_min,
            w_max=config.branch.w_max,
            nw=2,
        ),
        dtype=dtype,
        device=config.backend.device,
    )
    n = grid.spec.nr * grid.spec.nw
    state = torch.zeros(zero_p2_tangent_state(grid).shape, dtype=torch.complex128, device=grid.device)
    state[0:n] = 1.0
    state[n : 2 * n] = 1.0j
    wave_number = math.pi / max(config.branch.r_max, 1.0e-300)
    kinetic_energy = config.branch.hbar**2 * wave_number**2 / (2.0 * config.branch.particle_mass)
    omega = kinetic_energy / config.branch.hbar
    terms = p2_driven_frequency_terms(state, grid, config, omega=omega, cap_enabled=False)
    static_kinetic = torch.zeros_like(state)
    static_kinetic[0 : 2 * n] = kinetic_energy * state[0 : 2 * n]
    residual = static_kinetic[0 : 2 * n] + terms[0 : 2 * n]
    residual_norm = float(torch.linalg.norm(residual).detach().cpu().item())
    reference_norm = float(torch.linalg.norm(static_kinetic[0 : 2 * n]).detach().cpu().item())
    old_sign_time = torch.zeros_like(state)
    old_sign_time[0:n] = -1j * config.branch.hbar * omega * state[n : 2 * n]
    old_sign_time[n : 2 * n] = 1j * config.branch.hbar * omega * state[0:n]
    old_residual = static_kinetic[0 : 2 * n] + old_sign_time[0 : 2 * n]
    old_relative = float(
        torch.linalg.norm(old_residual).detach().cpu().item() / max(reference_norm, 1.0e-300)
    )
    relative = residual_norm / max(reference_norm, 1.0e-300)
    return {
        "branch_vector": "lane0=1,lane1=i",
        "wave_number": float(wave_number),
        "kinetic_energy": float(kinetic_energy),
        "omega": float(omega),
        "hbar_omega": float(config.branch.hbar * omega),
        "relative_residual": relative,
        "old_negated_sign_relative_residual": old_relative,
        "passed": relative <= 1.0e-14 and old_relative > 1.0,
    }


def driven_observable_convergence_gate(
    observable_summary: list[dict[str, Any]],
    config: HarnessConfig,
) -> dict[str, Any]:
    """Gate driven response observables, preferring the clean all-observable claim."""

    gated = [
        row
        for row in observable_summary
        if row["verdict"] not in {"null diagnostic", "solver-floor diagnostic"}
    ]
    threshold = float(config.p2_driven.min_observable_order)
    all_failing = [
        row
        for row in gated
        if row["last_observed_order"] is None or row["last_observed_order"] < threshold
    ]
    all_gate = {
        "mode": "all_non_floor_observables",
        "threshold": threshold,
        "gated_observable_names": [row["observable"] for row in gated],
        "failing_observable_names": [row["observable"] for row in all_failing],
        "gated_count": len(gated),
        "failed_count": len(all_failing),
        "passed": bool(gated) and not all_failing,
    }
    scoped_names = tuple(config.p2_driven.scoped_convergence_observables)
    rows_by_name = {row["observable"]: row for row in gated}
    scoped_rows = [rows_by_name[name] for name in scoped_names if name in rows_by_name]
    missing_scoped = [name for name in scoped_names if name not in rows_by_name]
    scoped_threshold = float(config.p2_driven.scoped_min_observable_order)
    scoped_failing = [
        row
        for row in scoped_rows
        if row["last_observed_order"] is None or row["last_observed_order"] < scoped_threshold
    ]
    scoped_gate = {
        "mode": "resolution_limited_scoped_observables",
        "threshold": scoped_threshold,
        "gated_observable_names": [row["observable"] for row in scoped_rows],
        "missing_observable_names": missing_scoped,
        "failing_observable_names": [row["observable"] for row in scoped_failing] + missing_scoped,
        "excluded_observable_names": [
            row["observable"] for row in gated if row["observable"] not in scoped_names
        ],
        "gated_count": len(scoped_rows),
        "failed_count": len(scoped_failing) + len(missing_scoped),
        "rationale": config.p2_driven.scoped_convergence_rationale,
        "passed": bool(scoped_rows) and not scoped_failing and not missing_scoped,
    }
    selected = all_gate if all_gate["passed"] else scoped_gate
    return {
        **selected,
        "clean_all_observable_gate": all_gate,
        "scoped_gate": scoped_gate,
        "resolution_limited": not all_gate["passed"],
        "rationale": (
            "The clean gate first evaluates all non-null, non-solver-floor driven "
            "response observables at min_observable_order.  If that fails on the "
            "omega=6 ladder, the configured scoped fallback may pass only the "
            "predeclared observables in scoped_convergence_observables at "
            "scoped_min_observable_order while reporting every excluded lane."
        ),
    }


def driven_resolution_diagnostics(
    runtime_levels: list[dict[str, Any]],
    config: HarnessConfig,
    *,
    omega: float,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for row in runtime_levels:
        grid = row["grid_object"]
        k = _wall_effective_wave_number(grid, config, omega=omega)
        wavelength = (2.0 * math.pi / k) if k > 0.0 else float("inf")
        points_per_wavelength = wavelength / grid.dw if math.isfinite(wavelength) else float("inf")
        rows.append(
            {
                "grid": row["grid"],
                "nw": int(grid.spec.nw),
                "dw": float(grid.dw),
                "wall_effective_k": float(k),
                "wavelength": float(wavelength),
                "points_per_wavelength": float(points_per_wavelength),
            }
        )
    return rows


def _wall_effective_wave_number(
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
) -> float:
    t_w_faces, t_omega, k_eta = p2_wall_coefficients(grid, config.p2_tangent)
    t_w_center = 0.5 * (t_w_faces[:-1] + t_w_faces[1:])
    w0, w1 = config.p2_driven.reflection_fit_window
    mask = (grid.w_centers >= w0) & (grid.w_centers <= w1)
    if int(torch.count_nonzero(mask).detach().cpu().item()) < 3:
        mask = torch.ones_like(grid.w_centers, dtype=torch.bool)
    ell = config.p2_tangent.spherical_l * (config.p2_tangent.spherical_l + 1)
    restoring = torch.mean(k_eta[mask] + ell * t_omega[mask])
    tension = torch.mean(t_w_center[mask])
    k2 = (config.p2_tangent.placeholder_mu_eta * omega**2 - float(restoring)) / float(tension)
    return math.sqrt(max(k2, 0.0))


def wall_lane_reflection_coefficient(
    state: torch.Tensor,
    grid: TensorProductGrid,
    config: HarnessConfig,
    *,
    omega: float,
) -> dict[str, float]:
    fields = unpack_p2_tangent_fields(state.to(torch.complex128), grid)
    k = _wall_effective_wave_number(grid, config, omega=omega)
    w0, w1 = config.p2_driven.reflection_fit_window
    mask = (grid.w_centers >= w0) & (grid.w_centers <= w1)
    if int(torch.count_nonzero(mask).detach().cpu().item()) < 3 or k <= 0.0:
        return {
            "reflection_coefficient": float("inf"),
            "outgoing_amplitude": 0.0,
            "incoming_amplitude": float("inf"),
            "effective_k": float(k),
            "fit_cells": int(torch.count_nonzero(mask).detach().cpu().item()),
            "fit_residual_l2": float("inf"),
        }
    w = grid.w_centers[mask].detach().cpu().numpy()
    eta = fields.eta[mask].detach().cpu().numpy().astype(np.complex128, copy=False)
    design = np.column_stack([np.exp(1j * k * w), np.exp(-1j * k * w)])
    coeffs, residuals, _rank, _svals = np.linalg.lstsq(design, eta, rcond=None)
    outgoing = float(abs(coeffs[0]))
    incoming = float(abs(coeffs[1]))
    fitted = design @ coeffs
    fit_residual = float(np.linalg.norm(eta - fitted) / max(np.linalg.norm(eta), 1.0e-300))
    return {
        "reflection_coefficient": incoming / max(outgoing, 1.0e-300),
        "outgoing_amplitude": outgoing,
        "incoming_amplitude": incoming,
        "effective_k": float(k),
        "fit_cells": int(mask.detach().cpu().numpy().sum()),
        "fit_residual_l2": fit_residual,
    }


def _wall_probe_matrix_and_source(
    config: HarnessConfig,
    *,
    w_max: float,
    nw: int,
    cap_enabled: bool,
) -> tuple[WallGrid, csc_matrix, np.ndarray, dict[str, float]]:
    dtype = configure_backend(config.backend)
    grid = WallGrid.create(
        config_wall_grid_spec_from_values(config.branch.w_min, w_max, nw),
        dtype=dtype,
        device=config.backend.device,
    )
    p2_cfg = config.p2_tangent
    driven = config.p2_driven
    omega = driven.propagating_omega
    t_w_value = p2_cfg.placeholder_t_w_base
    t_omega_value = p2_cfg.placeholder_t_omega_base
    k_eta_value = p2_cfg.placeholder_k_eta_base
    t_w = torch.full((nw + 1,), t_w_value, dtype=dtype, device=config.backend.device)
    t_omega = torch.full((nw,), t_omega_value, dtype=dtype, device=config.backend.device)
    k_eta = torch.full((nw,), k_eta_value, dtype=dtype, device=config.backend.device)
    lower = BoundaryCondition.neumann(0.0)
    upper = BoundaryCondition.neumann(0.0)
    columns = []
    for idx in range(nw):
        basis = torch.zeros(nw, dtype=dtype, device=config.backend.device)
        basis[idx] = 1.0
        applied = wall_s_eta_operator(
            basis,
            grid,
            t_w_faces=t_w,
            t_omega_centers=t_omega,
            k_eta_centers=k_eta,
            spherical_l=p2_cfg.spherical_l,
            lower_bc=lower,
            upper_bc=upper,
        )
        columns.append(applied.detach().cpu().numpy())
    static = np.column_stack(columns).astype(np.complex128)
    sigma = p2_cap_profile_numpy(
        grid.w_centers.detach().cpu().numpy(),
        w_max=config.branch.w_max,
        width=driven.cap_width,
        strength=driven.cap_strength,
        enabled=cap_enabled,
    )
    diagonal = -p2_cfg.placeholder_mu_eta * omega**2 - 1j * sigma
    matrix = csc_matrix(static + np.diag(diagonal))
    w = grid.w_centers.detach().cpu().numpy()
    length = config.branch.w_max - config.branch.w_min
    source_center = config.branch.w_min + 0.28 * length
    source_width = 0.08 * length
    k2 = (
        p2_cfg.placeholder_mu_eta * omega**2
        - (k_eta_value + p2_cfg.spherical_l * (p2_cfg.spherical_l + 1) * t_omega_value)
    ) / t_w_value
    k = math.sqrt(max(k2, 0.0))
    source = np.exp(-((w - source_center) / source_width) ** 2) * np.exp(1j * k * w)
    diagnostics = {
        "effective_k": float(k),
        "source_center": float(source_center),
        "source_width": float(source_width),
        "cap_width": float(driven.cap_width if cap_enabled else 0.0),
        "cap_strength": float(driven.cap_strength if cap_enabled else 0.0),
    }
    return grid, matrix, source.astype(np.complex128), diagnostics


def _solve_wall_probe(
    config: HarnessConfig,
    *,
    w_max: float,
    nw: int,
    cap_enabled: bool,
) -> tuple[WallGrid, np.ndarray, dict[str, float]]:
    grid, matrix, source, diagnostics = _wall_probe_matrix_and_source(
        config,
        w_max=w_max,
        nw=nw,
        cap_enabled=cap_enabled,
    )
    factor = splu(matrix, permc_spec=config.p2_tangent.newton.preconditioner.permutation)
    solution = factor.solve(source)
    residual = matrix @ solution - source
    diagnostics.update(
        {
            "residual_linf": float(np.max(np.abs(residual))),
            "residual_l2": float(np.linalg.norm(residual)),
            "factor_nnz_l": int(factor.L.nnz),
            "factor_nnz_u": int(factor.U.nnz),
            "nw": int(nw),
            "w_max": float(w_max),
        }
    )
    return grid, solution, diagnostics


def _wall_probe_contamination_metric(
    base_grid: WallGrid,
    base_solution: np.ndarray,
    reference_grid: WallGrid,
    reference_solution: np.ndarray,
    config: HarnessConfig,
) -> dict[str, float]:
    interior_end = min(
        config.p2_driven.reflection_fit_window[1],
        config.branch.w_max - config.p2_driven.cap_width - 0.25 * config.p2_driven.cap_width,
    )
    base_w = base_grid.w_centers.detach().cpu().numpy()
    ref_w = reference_grid.w_centers.detach().cpu().numpy()
    mask = (base_w >= config.p2_driven.reflection_fit_window[0]) & (base_w <= interior_end)
    if int(mask.sum()) < 3:
        mask = base_w <= interior_end
    ref_values = reference_solution[: base_solution.size][mask]
    diff = base_solution[mask] - ref_values
    weights = base_grid.cell_widths.detach().cpu().numpy()[mask]
    numerator = math.sqrt(float(np.sum(np.abs(diff) ** 2 * weights)))
    denominator = math.sqrt(float(np.sum(np.abs(ref_values) ** 2 * weights)))
    if not np.allclose(base_w[mask], ref_w[: base_solution.size][mask], rtol=0.0, atol=1.0e-13):
        raise ValueError("wall probe base/reference coordinates do not align")
    return {
        "reflection_coefficient": numerator / max(denominator, 1.0e-300),
        "interior_l2_change": numerator,
        "interior_signal_l2_reference": denominator,
        "fit_cells": int(mask.sum()),
        "window_w_min": float(base_w[mask][0]),
        "window_w_max": float(base_w[mask][-1]),
    }


def run_reflection_study(config: HarnessConfig) -> dict[str, Any]:
    base_nw = max(64, int(config.p2_driven.reflection_grid[1]) * 8)
    base_w_max = config.branch.w_max
    base_dw = (base_w_max - config.branch.w_min) / base_nw
    extension_cells = max(8, int(round(config.p2_driven.cap_width / base_dw)))
    reference_nw = base_nw + extension_cells
    reference_w_max = base_w_max + extension_cells * base_dw
    omega = config.p2_driven.propagating_omega

    absorbed_grid, absorbed_state, absorbed_solve = _solve_wall_probe(
        config,
        w_max=base_w_max,
        nw=base_nw,
        cap_enabled=True,
    )
    absorbed_ref_grid, absorbed_ref_state, absorbed_ref_solve = _solve_wall_probe(
        config,
        w_max=reference_w_max,
        nw=reference_nw,
        cap_enabled=True,
    )
    control_grid, control_state, control_solve = _solve_wall_probe(
        config,
        w_max=base_w_max,
        nw=base_nw,
        cap_enabled=False,
    )
    control_ref_grid, control_ref_state, control_ref_solve = _solve_wall_probe(
        config,
        w_max=reference_w_max,
        nw=reference_nw,
        cap_enabled=False,
    )
    absorbed = _wall_probe_contamination_metric(
        absorbed_grid,
        absorbed_state,
        absorbed_ref_grid,
        absorbed_ref_state,
        config,
    )
    control = _wall_probe_contamination_metric(
        control_grid,
        control_state,
        control_ref_grid,
        control_ref_state,
        config,
    )
    absorbed.update(absorbed_solve)
    control.update(control_solve)
    floor = min(
        config.p2_driven.reflection_max_relative,
        config.p2_driven.reflection_floor_multiplier
        * config.p2_driven.step4_reference_relative_floor,
    )
    contrast = control["reflection_coefficient"] / max(
        absorbed["reflection_coefficient"],
        1.0e-300,
    )
    return {
        "method": (
            "Domain-extension contamination metric on a standalone reduced 1D "
            "wall_s_eta_operator proxy for the retained wall lane: solve the proxy "
            "wall frequency operator on a base domain and on an aligned extended "
            "domain, then measure the fixed-window weighted L2 difference divided "
            "by the extended-domain signal.  This does not characterize matter, "
            "Maxwell, or full six-lane coupled-operator reflection."
        ),
        "non_circularity": (
            "The coordinate window, source, frequency, grid spacing, wall coefficients, "
            "and norm are fixed before absorber on/off is compared; the reflecting "
            "control uses the identical metric with sigma=0."
        ),
        "cap_profile_note": GENERIC_CAP_NOTE,
        "floor_note": (
            "The reflection floor multiplier is 8.0 times the Step-4 raw-field floor; "
            "it is a conditioning-motivated free choice, not a derived bound."
        ),
        "omega": float(omega),
        "grid": f"wall_nw_{base_nw}_reference_nw_{reference_nw}",
        "absorbed": absorbed,
        "absorbed_reference": absorbed_ref_solve,
        "reflecting_control": control,
        "reflecting_control_reference": control_ref_solve,
        "target_blind_floor": floor,
        "control_contrast": float(contrast),
        "pass_checks": {
            "absorbed_below_floor": absorbed["reflection_coefficient"] <= floor,
            "control_materially_higher": contrast >= config.p2_driven.reflection_control_min_contrast,
            "control_fails_floor": control["reflection_coefficient"] > floor,
            "propagating_wall_k_real": absorbed["effective_k"] > 0.0,
        },
    }


def _driven_first_time_lane_expressions() -> dict[str, Any]:
    lane0, lane1, time, hbar, omega = sp.symbols("lane0 lane1 t hbar omega", real=True)
    top = (lane0 + sp.I * lane1) * sp.exp(-sp.I * omega * time)
    bottom = (lane0 - sp.I * lane1) * sp.exp(sp.I * omega * time)
    doublet_time = sp.Matrix(
        [
            sp.I * hbar * sp.diff(top, time),
            sp.I * hbar * sp.diff(bottom, time),
        ]
    ).subs(time, 0)
    lane_time = sp.Matrix(
        [
            (doublet_time[0] + doublet_time[1]) / 2,
            (doublet_time[0] - doublet_time[1]) / (2 * sp.I),
        ]
    )
    return {
        "lane0": sp.lambdify((lane0, lane1, hbar, omega), lane_time[0], "numpy"),
        "lane1": sp.lambdify((lane0, lane1, hbar, omega), lane_time[1], "numpy"),
        "expr": sp.simplify(lane_time),
    }


def _driven_second_time_factor_expression() -> dict[str, Any]:
    time, omega = sp.symbols("t omega", real=True)
    factor = sp.diff(sp.exp(-sp.I * omega * time), time, 2).subs(time, 0)
    return {
        "factor": sp.lambdify((omega,), factor, "numpy"),
        "expr": sp.simplify(factor),
    }


_DRIVEN_FIRST_TIME = _driven_first_time_lane_expressions()
_DRIVEN_SECOND_TIME = _driven_second_time_factor_expression()


def _run_tensor_matter_frequency_cap_mms(
    config: HarnessConfig,
    *,
    matter_time_sign: int = 1,
) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.mms.p2_driven
    zero_flux = BoundaryCondition.neumann(0.0)
    cap_config = replace(
        config,
        p2_driven=replace(
            config.p2_driven,
            cap_width=cfg.cap_width,
            cap_strength=cfg.cap_strength,
            cap_enabled=True,
        ),
    )
    full_config = cap_config.to_dict()

    def build_level(level: tuple[int, int]):
        nr, nw = level
        grid = TensorProductGrid.create(
            TensorGridSpec(r_max=cfg.r_max, nr=nr, w_min=cfg.w_min, w_max=cfg.w_max, nw=nw),
            dtype=dtype,
            device=config.backend.device,
        )
        return grid, f"nr_{nr}_nw_{nw}", max(grid.dr, grid.dw), grid.to_dict(), grid.cell_volumes

    def evaluate_level(grid: TensorProductGrid):
        rr, ww = np.meshgrid(
            grid.r_centers.detach().cpu().numpy(),
            grid.w_centers.detach().cpu().numpy(),
            indexing="ij",
        )
        x_np = _P2_CENTRIFUGAL["field"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        y_np = 0.37 * x_np
        lap_x_np = _P2_CENTRIFUGAL["operator"](
            rr,
            ww,
            cfg.r_max,
            cfg.w_min,
            cfg.w_max,
            cfg.spherical_l,
        )
        lap_y_np = 0.37 * lap_x_np
        bg_np = 0.85 + 0.11 * x_np + 0.07 * (1.0 + np.cos(np.pi * ww / cfg.w_max))
        anomalous_np = 5.0 * cfg.anomalous_eos_K * bg_np**8
        sigma = p2_cap_profile_numpy(
            grid.w_centers.detach().cpu().numpy(),
            w_max=cfg.w_max,
            width=cfg.cap_width,
            strength=cfg.cap_strength,
            enabled=True,
        )[None, :]
        scale = -(cfg.hbar**2 / (2.0 * cfg.particle_mass))
        time_x = _DRIVEN_FIRST_TIME["lane0"](x_np, y_np, cfg.hbar, cfg.omega)
        time_y = _DRIVEN_FIRST_TIME["lane1"](x_np, y_np, cfg.hbar, cfg.omega)
        exact_x = (
            scale * lap_x_np
            + anomalous_np * x_np
            + time_x
            - 1j * sigma * x_np
        )
        exact_y = (
            scale * lap_y_np
            - anomalous_np * y_np
            + time_y
            - 1j * sigma * y_np
        )
        x = _torch_from_numpy(x_np, dtype=dtype, device=config.backend.device)
        y = _torch_from_numpy(y_np, dtype=dtype, device=config.backend.device)
        discrete_x = scale * tensor_laplacian_with_spherical_l(
            x,
            grid,
            cfg.spherical_l,
            zero_flux,
            zero_flux,
            zero_flux,
        ).to(torch.complex128)
        discrete_y = scale * tensor_laplacian_with_spherical_l(
            y,
            grid,
            cfg.spherical_l,
            zero_flux,
            zero_flux,
            zero_flux,
        ).to(torch.complex128)
        anomalous = _torch_from_numpy(anomalous_np, dtype=dtype, device=config.backend.device).to(
            torch.complex128
        )
        state = torch.zeros(
            zero_p2_tangent_state(grid).shape,
            dtype=torch.complex128,
            device=config.backend.device,
        )
        n = grid.spec.nr * grid.spec.nw
        state[0:n] = x.to(torch.complex128).reshape(-1)
        state[n : 2 * n] = y.to(torch.complex128).reshape(-1)
        if matter_time_sign == 1:
            driven_terms = unpack_p2_tangent_fields(
                p2_driven_frequency_terms(
                    state,
                    grid,
                    cap_config,
                    omega=cfg.omega,
                    cap_enabled=True,
                ),
                grid,
            )
            driven_x = driven_terms.psi_real
            driven_y = driven_terms.psi_imag
        elif matter_time_sign == -1:
            sigma_t = _torch_from_numpy(sigma, dtype=dtype, device=config.backend.device).to(
                torch.complex128
            )
            driven_x = -1j * cfg.hbar * cfg.omega * y.to(torch.complex128) - 1j * sigma_t * x
            driven_y = 1j * cfg.hbar * cfg.omega * x.to(torch.complex128) - 1j * sigma_t * y
        else:
            raise ValueError("matter_time_sign must be +1 or -1")
        discrete = torch.stack(
            [
                discrete_x + anomalous * x.to(torch.complex128) + driven_x,
                discrete_y - anomalous * y.to(torch.complex128) + driven_y,
            ],
            dim=0,
        )
        exact = _torch_from_numpy(
            np.stack([exact_x, exact_y], axis=0),
            dtype=torch.complex128,
            device=config.backend.device,
        )
        return discrete, exact, {
            "omega": cfg.omega,
            "cap_width": cfg.cap_width,
            "time_sign": matter_time_sign,
            "anomalous_term": "g=h'(rho_bg)*psi_bg^2=5*K*psi_bg^8, real nonzero",
            "anomalous_min": float(np.min(anomalous_np)),
            "anomalous_max": float(np.max(anomalous_np)),
        }

    return asdict(
        run_convergence_study(
            name=f"{cfg.name}_matter",
            description=(
                "Matter BdG time term in real/imag lanes plus a generic additive CAP; "
                "the differential carrier is the P2 scalar operator with a nonzero "
                "anomalous h'(rho_bg)*psi_bg^2*delta_psi* block."
            ),
            continuum_source="compact lines 1406-1418 and boundary methods sections 3, 5.",
            manufactured_field=str(_P2_CENTRIFUGAL["field_expr"]),
            forcing_derivation=(
                "SymPy applied i*hbar*d_t to the manufactured doublet "
                "(delta_psi,delta_psi*)=(x+i*y, x-i*y) with exp(-i*omega*t), "
                "then transformed back to real/imag lanes.  The exact carrier also "
                "includes an independently evaluated anomalous g*delta_psi* block and "
                "the additive -i*sigma CAP coefficient."
            ),
            levels=cfg.tensor_grid_levels,
            build_level=build_level,
            evaluate_level=evaluate_level,
            config=full_config,
            run_root=config.run_root,
            min_observed_order=cfg.min_observed_order,
            final_error_max=cfg.tensor_final_error_max,
            config_hash=config.config_hash(),
        )
    )


def _run_tensor_maxwell_frequency_cap_mms(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.mms.p2_driven
    cap_config = replace(
        config,
        p2_driven=replace(
            config.p2_driven,
            cap_width=cfg.cap_width,
            cap_strength=cfg.cap_strength,
            cap_enabled=True,
        ),
    )
    full_config = cap_config.to_dict()
    zero_flux = BoundaryCondition.neumann(0.0)

    def build_level(level: tuple[int, int]):
        nr, nw = level
        grid = TensorProductGrid.create(
            TensorGridSpec(r_max=cfg.r_max, nr=nr, w_min=cfg.w_min, w_max=cfg.w_max, nw=nw),
            dtype=dtype,
            device=config.backend.device,
        )
        return grid, f"nr_{nr}_nw_{nw}", max(grid.dr, grid.dw), grid.to_dict(), grid.cell_volumes

    def evaluate_level(grid: TensorProductGrid):
        r_np = grid.r_centers.detach().cpu().numpy()
        w_np = grid.w_centers.detach().cpu().numpy()
        rr, ww = np.meshgrid(r_np, w_np, indexing="ij")
        a0_np = _P2_MAXWELL_ANGULAR["a0"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        ar_np = _P2_MAXWELL_ANGULAR["ar"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        aw_np = _P2_MAXWELL_ANGULAR["aw"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        static_np = np.asarray(
            _P2_MAXWELL_ANGULAR["operator"](
                rr,
                ww,
                cfg.r_max,
                cfg.w_min,
                cfg.w_max,
                cfg.spherical_l,
                cfg.localization_width,
                cfg.xi,
            ),
            dtype=np.float64,
        ).squeeze()
        weight_np = _P2_MAXWELL_ANGULAR["weight"](w_np, cfg.localization_width)
        temporal_weight_np = (
            localization_weight_torch(grid.w_centers, cap_config.branch).detach().cpu().numpy()
        )
        second_time_factor = complex(_DRIVEN_SECOND_TIME["factor"](cfg.omega))
        sigma = p2_cap_profile_numpy(
            w_np,
            w_max=cfg.w_max,
            width=cfg.cap_width,
            strength=cfg.cap_strength,
            enabled=True,
        )
        lanes_np = np.stack([a0_np, ar_np, aw_np], axis=0)
        exact_np = static_np + (temporal_weight_np[None, None, :] * second_time_factor) * lanes_np
        exact_np = exact_np - 1j * sigma[None, None, :] * lanes_np
        weight_centers = _torch_from_numpy(weight_np, dtype=dtype, device=config.backend.device)
        weight_faces = _torch_from_numpy(
            _P2_MAXWELL_ANGULAR["weight"](
                grid.w_faces.detach().cpu().numpy(),
                cfg.localization_width,
            ),
            dtype=dtype,
            device=config.backend.device,
        )
        discrete_static = localized_maxwell_operator_with_spherical_l(
            _torch_from_numpy(a0_np, dtype=dtype, device=config.backend.device),
            _torch_from_numpy(ar_np, dtype=dtype, device=config.backend.device),
            _torch_from_numpy(aw_np, dtype=dtype, device=config.backend.device),
            grid,
            spherical_l=cfg.spherical_l,
            xi=cfg.xi,
            weight_w_centers=weight_centers,
            weight_w_faces=weight_faces,
            a0_radial_outer_bc=zero_flux,
            a0_w_lower_bc=zero_flux,
            a0_w_upper_bc=zero_flux,
        ).to(torch.complex128)
        state = torch.zeros(
            zero_p2_tangent_state(grid).shape,
            dtype=torch.complex128,
            device=config.backend.device,
        )
        n = grid.spec.nr * grid.spec.nw
        state[2 * n : 3 * n] = _torch_from_numpy(a0_np, dtype=dtype, device=config.backend.device).to(
            torch.complex128
        ).reshape(-1)
        state[3 * n : 4 * n] = _torch_from_numpy(ar_np, dtype=dtype, device=config.backend.device).to(
            torch.complex128
        ).reshape(-1)
        state[4 * n : 5 * n] = _torch_from_numpy(aw_np, dtype=dtype, device=config.backend.device).to(
            torch.complex128
        ).reshape(-1)
        driven_terms = unpack_p2_tangent_fields(
            p2_driven_frequency_terms(
                state,
                grid,
                cap_config,
                omega=cfg.omega,
                cap_enabled=True,
            ),
            grid,
        )
        driven = torch.stack([driven_terms.a0, driven_terms.ar, driven_terms.aw], dim=0)
        exact = _torch_from_numpy(exact_np, dtype=torch.complex128, device=config.backend.device)
        return discrete_static + driven, exact, {
            "omega": cfg.omega,
            "cap_width": cfg.cap_width,
            "temporal_factor": str(_DRIVEN_SECOND_TIME["expr"]),
            "temporal_truncation": MAXWELL_TEMPORAL_TRUNCATION_NOTE,
            "spatial_scalarization_caveat": "Step-8a Ar is componentwise engineering smoke.",
        }

    return asdict(
        run_convergence_study(
            name=f"{cfg.name}_maxwell",
            description=(
                "Temporal Maxwell engineering-smoke diagonal -Z(w)*omega^2 retained-lane "
                "self-term plus additive CAP.  This temporal truncation is additional "
                "to the Step-8a spatial Ar scalarization caveat."
            ),
            continuum_source="compact lines 1423-1434 and boundary methods sections 3, 5.",
            manufactured_field=(
                f"A0={_P2_MAXWELL_ANGULAR['a0_expr']}; "
                f"Ar={_P2_MAXWELL_ANGULAR['ar_expr']}; "
                f"Aw={_P2_MAXWELL_ANGULAR['aw_expr']}"
            ),
            forcing_derivation=(
                "SymPy evaluated d_t^2 exp(-i*omega*t)=-omega^2 and applied only "
                "the retained diagonal Z(w)*d_t^2 A_N self-term plus -i*sigma*A_N.  "
                "The omega-linear A0<->Ai terms and gauge temporal pieces are "
                "explicitly not included."
            ),
            levels=cfg.tensor_grid_levels,
            build_level=build_level,
            evaluate_level=evaluate_level,
            config=full_config,
            run_root=config.run_root,
            min_observed_order=cfg.min_observed_order,
            final_error_max=cfg.tensor_final_error_max,
            config_hash=config.config_hash(),
        )
    )


def _run_wall_frequency_cap_mms(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.mms.p2_driven
    cap_config = replace(
        config,
        p2_driven=replace(
            config.p2_driven,
            cap_width=cfg.cap_width,
            cap_strength=cfg.cap_strength,
            cap_enabled=True,
        ),
    )
    full_config = cap_config.to_dict()
    lower = BoundaryCondition.neumann(
        float(_WALL["lower_outward_derivative"](cfg.w_min, cfg.w_max))
    )
    upper = BoundaryCondition.neumann(
        float(_WALL["upper_outward_derivative"](cfg.w_min, cfg.w_max))
    )

    def build_level(nw: int):
        grid = WallGrid.create(
            config_wall_grid_spec_from_values(cfg.w_min, cfg.w_max, nw),
            dtype=dtype,
            device=config.backend.device,
        )
        return grid, f"nw_{nw}", grid.dw, grid.to_dict(), grid.cell_widths

    def evaluate_level(grid: WallGrid):
        w_centers = grid.w_centers.detach().cpu().numpy()
        w_faces = grid.w_faces.detach().cpu().numpy()
        eta_np = _WALL["eta"](w_centers, cfg.w_min, cfg.w_max)
        static_np = _WALL["operator"](
            w_centers,
            cfg.w_min,
            cfg.w_max,
            cfg.spherical_l,
            1.1,
            0.2,
            0.8,
            0.1,
            0.9,
            0.15,
        )
        sigma = p2_cap_profile_numpy(
            w_centers,
            w_max=cfg.w_max,
            width=cfg.cap_width,
            strength=cfg.cap_strength,
            enabled=True,
        )
        second_time_factor = complex(_DRIVEN_SECOND_TIME["factor"](cfg.omega))
        exact_np = static_np + cfg.mu_eta * second_time_factor * eta_np - 1j * sigma * eta_np
        eta = _torch_from_numpy(eta_np, dtype=dtype, device=config.backend.device)
        t_w = _torch_from_numpy(
            _WALL["t_w"](w_faces, cfg.w_min, cfg.w_max, 1.1, 0.2),
            dtype=dtype,
            device=config.backend.device,
        )
        t_omega = _torch_from_numpy(
            _WALL["t_omega"](w_centers, cfg.w_min, cfg.w_max, 0.8, 0.1),
            dtype=dtype,
            device=config.backend.device,
        )
        k_eta = _torch_from_numpy(
            _WALL["k_eta"](w_centers, cfg.w_min, cfg.w_max, 0.9, 0.15),
            dtype=dtype,
            device=config.backend.device,
        )
        discrete = wall_s_eta_operator(
            eta,
            grid,
            t_w_faces=t_w,
            t_omega_centers=t_omega,
            k_eta_centers=k_eta,
            spherical_l=cfg.spherical_l,
            lower_bc=lower,
            upper_bc=upper,
        ).to(torch.complex128)
        sigma_t = _torch_from_numpy(sigma, dtype=dtype, device=config.backend.device).to(
            torch.complex128
        )
        discrete = discrete - cfg.mu_eta * cfg.omega**2 * eta.to(torch.complex128)
        discrete = discrete - 1j * sigma_t * eta.to(torch.complex128)
        exact = _torch_from_numpy(exact_np, dtype=torch.complex128, device=config.backend.device)
        return discrete, exact, {
            "omega": cfg.omega,
            "cap_width": cfg.cap_width,
            "temporal_factor": str(_DRIVEN_SECOND_TIME["expr"]),
        }

    return asdict(
        run_convergence_study(
            name=f"{cfg.name}_wall",
            description="Wall mu_eta*d_t^2 term and additive CAP on the modal wall operator.",
            continuum_source="compact lines 1441-1451 and boundary methods sections 3, 5.",
            manufactured_field=str(_WALL["eta_expr"]),
            forcing_derivation=(
                "SymPy evaluated d_t^2 exp(-i*omega*t)=-omega^2, giving the "
                "mu_eta*d_t^2 eta term as -mu_eta*omega^2*eta, plus -i*sigma*eta."
            ),
            levels=cfg.wall_grid_levels,
            build_level=build_level,
            evaluate_level=evaluate_level,
            config=full_config,
            run_root=config.run_root,
            min_observed_order=cfg.min_observed_order,
            final_error_max=cfg.wall_final_error_max,
            config_hash=config.config_hash(),
        )
    )


def config_wall_grid_spec_from_values(w_min: float, w_max: float, nw: int):
    from .config import WallGridSpec

    return WallGridSpec(w_min=w_min, w_max=w_max, nw=nw)


def run_p2_driven_mms(config: HarnessConfig) -> dict[str, Any]:
    sections = {
        "matter_frequency_cap": _run_tensor_matter_frequency_cap_mms(config),
        "maxwell_frequency_cap": _run_tensor_maxwell_frequency_cap_mms(config),
        "wall_frequency_cap": _run_wall_frequency_cap_mms(config),
    }
    return {
        "sections": sections,
        "passed": all(section["passed"] for section in sections.values()),
    }


def run_response_vs_omega(
    config: HarnessConfig,
    *,
    background_state: torch.Tensor | None = None,
    grid: TensorProductGrid | None = None,
) -> list[dict[str, Any]]:
    dtype = configure_backend(config.backend)
    if grid is None:
        level = config.p2_driven.response_table_grid
        grid = _create_branch_grid(config.branch, level, dtype=dtype, device=config.backend.device)
    if background_state is None:
        background_state, _background = run_branch_continuation(
            config,
            grid,
            grid_name=f"p2_response_background_nr_{grid.spec.nr}_nw_{grid.spec.nw}",
        )
    rows = []
    for omega in config.p2_driven.drive_frequencies:
        state, solve = solve_driven_p2_tangent(
            background_state,
            grid,
            config,
            omega=omega,
            grid_name=f"p2_response_omega_{omega:.6g}_nr_{grid.spec.nr}_nw_{grid.spec.nw}",
            cap_enabled=True,
        )
        propagating_k = _wall_effective_wave_number(grid, config, omega=omega)
        rows.append(
            {
                "omega": float(omega),
                "wall_effective_k": propagating_k,
                "propagating_wall_lane": propagating_k > 0.0,
                **solve["surrogate_values"],
            }
        )
    return rows


def _driven_digest(results: dict[str, Any]) -> str:
    def sanitize(value: Any) -> Any:
        if isinstance(value, dict):
            return {
                key: sanitize(item)
                for key, item in value.items()
                if key not in {"manifest", "machine_readable_table"}
            }
        if isinstance(value, list):
            return [sanitize(item) for item in value]
        return value

    payload = {
        "mms": sanitize(results["mms"]),
        "static_tangent_fd_check": sanitize(results["static_tangent_fd_check"]),
        "driven_fd_check": sanitize(results["driven_fd_check"]),
        "omega_zero_static_limit": sanitize(results["omega_zero_static_limit"]),
        "matter_dispersion": sanitize(results["matter_dispersion"]),
        "wellposedness": sanitize(results["wellposedness"]),
        "level_rows": sanitize(results["level_rows"]),
        "observable_summary": sanitize(results["observable_summary"]),
        "convergence_gate": sanitize(results["convergence_gate"]),
        "resolution_diagnostics": sanitize(results["resolution_diagnostics"]),
        "reflection": sanitize(results["reflection"]),
        "response_vs_omega": sanitize(results["response_vs_omega"]),
        "pass_checks": results["pass_checks"],
        "asserted_checks": results["asserted_checks"],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def run_step8b(config: HarnessConfig | None = None) -> dict[str, Any]:
    if config is None:
        config = step8b_default_config()
    config = with_step8a_preconditioners(config)
    Path(config.run_root).mkdir(parents=True, exist_ok=True)

    mms = run_p2_driven_mms(config)
    runtime_levels = run_p2_driven_ladder(config, omega=config.p2_driven.primary_omega)
    level_rows = [_strip_runtime_level(level) for level in runtime_levels]
    for row in level_rows:
        row.pop("driven_state", None)
    observable_summary, observable_table = _summarize_p2_observables(
        level_rows,
        ratio=config.p2_tangent.refinement_ratio,
        expected_order=config.p2_tangent.expected_order,
    )
    fd_level = runtime_levels[0]
    boundaries = branch_boundary_conditions(config.branch)
    static_tangent_fd_check = p2_tangent_fd_check(
        fd_level["background_state"],
        fd_level["grid_object"],
        config.branch,
        config.p2_tangent,
        eos_K=config.branch.continuation_K_values[-1],
        boundaries=boundaries,
        epsilon=config.p2_tangent.newton.finite_difference_jvp_epsilon,
        seed=config.jacobian_check_seed,
    )
    driven_fd_check = p2_driven_fd_check(
        fd_level["background_state"],
        fd_level["grid_object"],
        config,
        omega=config.p2_driven.primary_omega,
        epsilon=config.p2_tangent.newton.finite_difference_jvp_epsilon,
        seed=config.jacobian_check_seed + 1,
    )
    omega_zero_static_limit = p2_omega_zero_static_limit_check(
        fd_level["background_state"],
        fd_level["grid_object"],
        config,
        seed=config.jacobian_check_seed + 2,
    )
    matter_dispersion = p2_matter_free_particle_dispersion_check(config)
    wellposedness = p2_driven_wellposedness_check(
        fd_level["background_state"],
        fd_level["grid_object"],
        config,
        omega=config.p2_driven.primary_omega,
    )
    reflection = run_reflection_study(config)
    response_vs_omega = run_response_vs_omega(
        config,
        background_state=fd_level["background_state"],
        grid=fd_level["grid_object"],
    )
    convergence_gate = driven_observable_convergence_gate(observable_summary, config)
    resolution_diagnostics = driven_resolution_diagnostics(
        runtime_levels,
        config,
        omega=config.p2_driven.primary_omega,
    )
    pass_checks = {
        "static_background_jacobian_gate_a_reverified": (
            static_tangent_fd_check["relative_residual"]
            <= config.p2_tangent.tangent_fd_relative_tol
            or static_tangent_fd_check["absolute_residual"]
            <= config.p2_tangent.tangent_fd_absolute_tol
        ),
        "driven_operator_internal_fd_check": (
            driven_fd_check["relative_residual"] <= config.p2_tangent.tangent_fd_relative_tol
            or driven_fd_check["absolute_residual"] <= config.p2_tangent.tangent_fd_absolute_tol
        ),
        "omega_zero_recovers_static_operator": (
            omega_zero_static_limit["relative_residual"] <= config.p2_tangent.tangent_fd_relative_tol
            or omega_zero_static_limit["absolute_residual"]
            <= config.p2_tangent.tangent_fd_absolute_tol
        ),
        "matter_time_sign_dispersion_sanity": bool(matter_dispersion["passed"]),
        "new_frequency_and_cap_terms_mms_order": mms["passed"],
        "driven_tangent_solves_converged": all(row["converged"] for row in level_rows),
        "complex_driven_operator_wellposed": (
            wellposedness["smallest_singular_value"]
            >= config.p2_driven.smallest_singular_min
        ),
        "scoped_surrogate_observable_convergence": convergence_gate["passed"],
        "reflection_metric_below_target_blind_floor": reflection["pass_checks"][
            "absorbed_below_floor"
        ],
        "reflecting_control_has_teeth": (
            reflection["pass_checks"]["control_materially_higher"]
            and reflection["pass_checks"]["control_fails_floor"]
        ),
        "propagating_drive_present": any(row["propagating_wall_lane"] for row in response_vs_omega),
    }
    asserted_checks = {
        "target_blind_surrogate_forcing_only_not_a_physics_gate": True,
        "physical_export_permitted_is_false_not_a_physics_gate": True,
        "matter_gauge_to_wall_source_deferred_not_a_physics_gate": True,
        "maxwell_radial_lane_scalarization_labelled_not_a_physics_gate": True,
        "maxwell_temporal_truncation_labelled_not_a_physics_gate": True,
        "raw_response_vs_omega_not_fitted_not_a_physics_gate": True,
    }
    asserted_check_notes = {
        "target_blind_surrogate_forcing_only_not_a_physics_gate": (
            "The driven RHS extends p2_surrogate_forcing with smooth numerical source data."
        ),
        "physical_export_permitted_is_false_not_a_physics_gate": (
            "Step 8b writes only run manifests/tables and a report; it does not import the firewalled model."
        ),
        "matter_gauge_to_wall_source_deferred_not_a_physics_gate": (
            "The compact S_eta^(psi,A) source remains open and is not invented here."
        ),
        "maxwell_radial_lane_scalarization_labelled_not_a_physics_gate": (
            "The Step-8a spatial Maxwell Ar lane remains componentwise and is not the full "
            "vector-harmonic l=2 reduction."
        ),
        "maxwell_temporal_truncation_labelled_not_a_physics_gate": (
            MAXWELL_TEMPORAL_TRUNCATION_NOTE
        ),
        "raw_response_vs_omega_not_fitted_not_a_physics_gate": (
            "The omega table reports raw magnitudes only and performs no coefficient extraction."
        ),
    }
    results: dict[str, Any] = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "source_citations": DRIVEN_SOURCE_CITATIONS,
        "method": {
            "representation": (
                "complex128 packed P2 fields with the same six lanes as Step 8a; "
                "the real Step-8a Jacobian is complexified without changing its entries."
            ),
            "matter_time_term": (
                "The +/-hbar*omega BdG doublet is represented in real/imag lanes as "
                "deltaR_re=+i*hbar*omega*psi_imag and deltaR_im=-i*hbar*omega*psi_real."
            ),
            "maxwell_time_term": MAXWELL_TEMPORAL_TRUNCATION_NOTE,
            "wall_time_term": "-mu_eta*omega^2*eta.",
            "absorber": GENERIC_CAP_NOTE,
            "drive_frequencies": list(config.p2_driven.drive_frequencies),
        },
        "mms": mms,
        "static_tangent_fd_check": static_tangent_fd_check,
        "driven_fd_check": driven_fd_check,
        "omega_zero_static_limit": omega_zero_static_limit,
        "matter_dispersion": matter_dispersion,
        "wellposedness": wellposedness,
        "level_rows": level_rows,
        "observable_summary": observable_summary,
        "observable_table": observable_table,
        "convergence_gate": convergence_gate,
        "resolution_diagnostics": resolution_diagnostics,
        "reflection": reflection,
        "response_vs_omega": response_vs_omega,
        "pass_checks": pass_checks,
        "asserted_checks": asserted_checks,
        "asserted_check_notes": asserted_check_notes,
        "passed": all(pass_checks.values()),
    }
    results["diagnostics_digest"] = _driven_digest(results)
    table_path = Path(config.run_root) / config.p2_driven.name / "p2_driven_table.json"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    results["machine_readable_table"] = str(table_path)
    table_path.write_text(
        json.dumps(results, indent=2, sort_keys=True, default=_json_default),
        encoding="utf-8",
    )
    return results


def write_step8b_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# Step 8b Driven P2 Tangent And Absorber")
    lines.append("")
    if results["passed"] and results["convergence_gate"].get("resolution_limited"):
        status = "PASS_WITH_RESOLUTION_LIMIT"
    else:
        status = "PASS" if results["passed"] else "FAIL"
    lines.append(f"Overall engineering gate: {status}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append(f"Diagnostics digest: `{results['diagnostics_digest']}`")
    lines.append("")
    lines.append(
        "**Scope:** engineering-smoke, target-blind driven tangent. Raw magnitudes are "
        "diagnostics only; no physical packet, no extraction map, and no fitted response coefficients."
    )
    lines.append("")
    lines.append("## Sources")
    lines.append("")
    for key, value in results["source_citations"].items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("## Driven Operator")
    lines.append("")
    for key in (
        "representation",
        "matter_time_term",
        "maxwell_time_term",
        "wall_time_term",
        "absorber",
    ):
        lines.append(f"- {key}: {results['method'][key]}")
    lines.append("")
    lines.append(
        "The omega->0, CAP-off sanity check gives "
        f"relative={results['omega_zero_static_limit']['relative_residual']:.6e}. "
        "This is a weak static-path consistency check that is zero by construction "
        "when omega=0 and CAP is off; it is not an omega/CAP physics certification."
    )
    lines.append(
        "Matter sign dispersion sanity check: "
        f"hbar*omega={results['matter_dispersion']['hbar_omega']:.6e}, "
        f"E(k)={results['matter_dispersion']['kinetic_energy']:.6e}, "
        f"relative={results['matter_dispersion']['relative_residual']:.6e}; "
        "the old negated sign would give "
        f"relative={results['matter_dispersion']['old_negated_sign_relative_residual']:.6e}."
    )
    lines.append("")
    lines.append("## Frequencies")
    lines.append("")
    lines.append(
        "Drive frequencies are target-blind numerical probes selected from near-static, "
        "intermediate, and wall-propagating regimes."
    )
    lines.append(
        _table(
            ["omega", "wall_effective_k", "propagating_wall_lane", "total_response_l2", "wall_eta_l2"],
            results["response_vs_omega"],
        )
    )
    lines.append("")
    lines.append("## MMS")
    lines.append("")
    for section in results["mms"]["sections"].values():
        lines.append(f"### {section['name']}")
        lines.append("")
        lines.append(section["description"])
        lines.append(f"Source: {section['continuum_source']}")
        lines.append(f"Forcing: {section['forcing_derivation']}")
        lines.append(
            _table(
                ["grid", "spacing", "error", "observed_order", "reference_norm"],
                section["rows"],
            )
        )
        lines.append("")
    lines.append("## Gate A And Well-Posedness")
    lines.append("")
    sfd = results["static_tangent_fd_check"]
    dfd = results["driven_fd_check"]
    well = results["wellposedness"]
    lines.append(
        "Static background JVP gate-A re-check: "
        f"relative={sfd['relative_residual']:.6e}, absolute={sfd['absolute_residual']:.6e}. "
        "This is internal consistency only; omega/CAP physics is certified by MMS and fidelity review."
    )
    lines.append(
        "Driven assembled operator FD check: "
        f"relative={dfd['relative_residual']:.6e}, absolute={dfd['absolute_residual']:.6e}."
    )
    lines.append(
        "Complex small-grid singular values: "
        f"sigma_min={well['smallest_singular_value']:.6e}, "
        f"condition={well['condition_number']:.6e}. "
        f"The configured wellposedness floor {results['config']['p2_driven']['smallest_singular_min']:.1e} "
        "is a loose numerical floor, not a sharp spectral margin."
    )
    lines.append("")
    lines.append("## Driven Solve Convergence")
    lines.append("")
    clean_gate = results["convergence_gate"]["clean_all_observable_gate"]
    scoped_gate = results["convergence_gate"]["scoped_gate"]
    if results["convergence_gate"].get("resolution_limited"):
        lines.append(
            "Resolution-limited fallback is active: the clean all-non-floor-observable "
            f"gate at bar {clean_gate['threshold']:.3f} FAILS.  The scoped convergence "
            f"claim uses the predeclared bar {scoped_gate['threshold']:.3f}; this is a "
            "reduced engineering claim, not a clean all-observable convergence pass."
        )
        lines.append(
            "Clean-gate failing observables: "
            + ", ".join(clean_gate["failing_observable_names"])
            + "."
        )
        lines.append(
            "Scoped gated observables: "
            + ", ".join(scoped_gate["gated_observable_names"])
            + "."
        )
        lines.append(
            "Scoped exclusions reported as under-resolved/drifting: "
            + ", ".join(scoped_gate["excluded_observable_names"])
            + "."
        )
        lines.append("Scoped rationale: " + scoped_gate["rationale"])
    else:
        lines.append(
            "The surrogate-observable convergence gate counts every non-null, "
            "non-solver-floor driven response observable at the configured bar "
            f"{clean_gate['threshold']:.3f}.  No post-run subset is selected."
        )
    lines.append("")
    lines.append(
        "Wall-lane wavelength resolution for omega=6 shows why the individual "
        "lane diagnostics are resolution-limited on this CPU-bounded ladder."
    )
    lines.append(
        _table(
            ["grid", "nw", "dw", "wall_effective_k", "wavelength", "points_per_wavelength"],
            results["resolution_diagnostics"],
        )
    )
    lines.append("")
    solve_rows = [
        {
            "level": row["level"],
            "grid": row["grid"],
            "omega": row["omega"],
            "dof": row["dof"],
            "background_residual": row["background"]["final_residual_linf"],
            "driven_residual_linf": row["solve"]["driven_residual_linf"],
            "converged": row["converged"],
        }
        for row in results["level_rows"]
    ]
    lines.append(
        _table(
            [
                "level",
                "grid",
                "omega",
                "dof",
                "background_residual",
                "driven_residual_linf",
                "converged",
            ],
            solve_rows,
        )
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
    lines.append("## Reflection")
    lines.append("")
    refl = results["reflection"]
    lines.append(refl["method"])
    lines.append(refl["non_circularity"])
    lines.append(refl["cap_profile_note"])
    lines.append(refl["floor_note"])
    reflection_rows = [
        {"case": "absorbed", **refl["absorbed"]},
        {"case": "reflecting_control", **refl["reflecting_control"]},
    ]
    lines.append(
        _table(
            [
                "case",
                "reflection_coefficient",
                "interior_l2_change",
                "interior_signal_l2_reference",
                "effective_k",
                "fit_cells",
            ],
            reflection_rows,
        )
    )
    lines.append(
        f"Target-blind floor: {refl['target_blind_floor']:.6e}; "
        f"control contrast: {refl['control_contrast']:.6e}."
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
    lines.append("## Provenance")
    lines.append("")
    lines.append(f"Machine-readable table: `{results['machine_readable_table']}`.")
    lines.append(
        "Run with: `PYTHONPATH=software/stage1_solver/src timeout 600 "
        "python -m stage1_solver.step8b_harness`."
    )
    lines.append(
        "Forward note for Step 8c: the OPEN matter/gauge-to-wall source remains deferred, "
        "and physical extraction/export gates remain untouched."
    )
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path
