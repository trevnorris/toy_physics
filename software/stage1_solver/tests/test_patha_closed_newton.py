from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import torch

from stage1_solver.backend import configure_backend
from stage1_solver.boundaries import BoundaryCondition
from stage1_solver.config import (
    BackendConfig,
    BranchSmokeConfig,
    PreconditionerConfig,
    TensorGridSpec,
    WallGridSpec,
)
from stage1_solver.coupled_branch import (
    branch_boundary_conditions,
    confinement_wall_to_matter_coefficient_torch,
    initial_closed_branch_state,
    patha_closed_branch_residual,
    patha_radial_reduced_return_source,
    unpack_closed_coupled_fields,
)
from stage1_solver.grid import TensorProductGrid, WallGrid
from stage1_solver.newton import finite_difference_jvp_check
from stage1_solver.p2_tangent import wall_to_matter_coefficient_torch
from stage1_solver.patha_closed_newton import (
    PathAClosedNewtonConfig,
    default_closed_branch_config,
    placeholder_provider_derivative_check,
    target_token_scan,
)
from stage1_solver.patha_static_balance import SSigmaSpec, static_balance_terms
from stage1_solver.preconditioners import (
    assemble_closed_coupled_colored_sparse_jacobian,
)
from stage1_solver.newton import PreconditionerBuildContext


def _small_grid(dtype: torch.dtype) -> TensorProductGrid:
    return TensorProductGrid.create(
        TensorGridSpec(r_max=1.4, nr=4, w_min=0.0, w_max=1.0, nw=4),
        dtype=dtype,
        device="cpu",
    )


def _small_branch() -> BranchSmokeConfig:
    return replace(
        default_closed_branch_config(),
        r_max=1.4,
        w_max=1.0,
        solve_grid=(4, 4),
        continuation_K_values=(0.08,),
        newton=replace(
            default_closed_branch_config().newton,
            preconditioner=PreconditionerConfig(
                type="colored_sparse_jacobian_lu",
                side="left",
                rebuild_policy="every_newton_step",
                stencil_radius=1,
                color_separation=3,
                factorization="splu",
            ),
        ),
    )


def test_static_balance_supports_natural_zero_traction_exit() -> None:
    dtype = configure_backend(BackendConfig())
    grid = WallGrid.create(
        WallGridSpec(w_min=0.0, w_max=1.0, nw=5),
        dtype=dtype,
        device="cpu",
    )
    values = torch.linspace(0.9, 1.0, grid.spec.nw, dtype=dtype)
    source = torch.zeros_like(values)
    terms = static_balance_terms(
        values,
        grid,
        s_sigma=SSigmaSpec.smooth_positive_placeholder(w_min=0.0, w_max=1.0),
        source=source,
        lower_bc=BoundaryCondition.dirichlet(1.0),
        upper_bc=BoundaryCondition.neumann(0.0),
    )
    assert torch.isfinite(terms.residual).all()
    assert terms.face_fluxes[-1].item() == 0.0


def test_return_source_reduces_full_density_with_shared_kernel() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _small_grid(dtype)
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    fields, _ = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)

    assert wall_to_matter_coefficient_torch is confinement_wall_to_matter_coefficient_torch
    density = fields.psi_real**2 + fields.psi_imag**2
    k1 = confinement_wall_to_matter_coefficient_torch(grid, branch, radius=fields.r0)
    manual = torch.sum(grid.radial_shell_volumes[:, None] * (-k1 * density), dim=0)
    reduced = patha_radial_reduced_return_source(fields, grid, branch)

    assert torch.allclose(reduced, manual)
    assert torch.max(reduced).item() < 0.0


def test_closed_residual_jvp_matches_finite_difference() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _small_grid(dtype)
    spec = SSigmaSpec.smooth_positive_placeholder(w_min=branch.w_min, w_max=branch.w_max)
    boundaries = branch_boundary_conditions(branch)
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    residual_fn = lambda x: patha_closed_branch_residual(
        x,
        grid,
        branch,
        eos_K=branch.continuation_K_values[0],
        boundaries=boundaries,
        s_sigma=spec,
    )
    check = finite_difference_jvp_check(
        residual_fn,
        state,
        epsilon=branch.newton.finite_difference_jvp_epsilon,
        seed=20260617,
    )
    assert check["relative_residual"] <= 1.0e-6 or check["absolute_residual"] <= 1.0e-6


def test_closed_colored_preconditioner_uses_wall_plus_mass_layout() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _small_branch()
    grid = _small_grid(dtype)
    spec = SSigmaSpec.smooth_positive_placeholder(w_min=branch.w_min, w_max=branch.w_max)
    boundaries = branch_boundary_conditions(branch)
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    residual_fn = lambda x: patha_closed_branch_residual(
        x,
        grid,
        branch,
        eos_K=branch.continuation_K_values[0],
        boundaries=boundaries,
        s_sigma=spec,
    )
    rhs = -residual_fn(state).detach().cpu().numpy()
    matrix, metadata = assemble_closed_coupled_colored_sparse_jacobian(
        PreconditionerBuildContext(
            residual_fn=residual_fn,
            x=state,
            rhs=rhs,
            iteration=1,
            config=branch.newton,
        ),
        grid,
    )
    expected = 5 * grid.spec.nr * grid.spec.nw + grid.spec.nw + 1
    assert matrix.shape == (expected, expected)
    assert metadata["layout"] == "5*cells+nw+1"
    assert metadata["state_size"] == expected


def test_placeholder_provider_derivative_check_and_scan(tmp_path: Path) -> None:
    dtype = configure_backend(BackendConfig())
    config = PathAClosedNewtonConfig(branch=_small_branch())
    check = placeholder_provider_derivative_check(
        config.resolved_s_sigma_spec(),
        dtype=dtype,
        device="cpu",
    )
    assert check["max_relative"] <= config.derivative_relative_tol
    assert target_token_scan([tmp_path / "missing.py"])["passed"] is True
