from __future__ import annotations

import math

import torch

from stage1_solver.backend import configure_backend
from stage1_solver.boundaries import BoundaryCondition
from stage1_solver.config import (
    BackendConfig,
    CubicGPEConfig,
    NewtonConfig,
    RadialGridSpec,
    TensorGridSpec,
)
from stage1_solver.grid import RadialGrid, TensorProductGrid
from stage1_solver.newton import finite_difference_jvp_check
from stage1_solver.operators import integrate, radial_fluxes, radial_laplacian, tensor_laplacian
from stage1_solver.physics import cubic_gpe_residual, gaussian_initial_state


def test_radial_constant_has_zero_laplacian_and_origin_flux() -> None:
    dtype = configure_backend(BackendConfig())
    grid = RadialGrid.create(RadialGridSpec(r_max=2.0, nr=32), dtype=dtype, device="cpu")
    values = torch.ones(grid.spec.nr, dtype=dtype)
    bc = BoundaryCondition.neumann(0.0)
    lap = radial_laplacian(values, grid, bc)
    fluxes = radial_fluxes(values, grid, bc)
    assert torch.max(torch.abs(lap)).item() < 1.0e-13
    assert abs(fluxes[0].item()) == 0.0


def test_discrete_mass_normalization() -> None:
    dtype = configure_backend(BackendConfig())
    grid = RadialGrid.create(RadialGridSpec(r_max=8.0, nr=32), dtype=dtype, device="cpu")
    cfg = CubicGPEConfig()
    state = gaussian_initial_state(grid, cfg)
    mass = integrate(state[:-1] * state[:-1], grid).item()
    assert math.isclose(mass, cfg.mass, rel_tol=0.0, abs_tol=1.0e-14)


def test_tensor_grid_constant_has_zero_laplacian() -> None:
    dtype = configure_backend(BackendConfig())
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=2.0, nr=16, w_min=-1.0, w_max=1.0, nw=8),
        dtype=dtype,
        device="cpu",
    )
    values = torch.ones((grid.spec.nr, grid.spec.nw), dtype=dtype)
    zero_flux = BoundaryCondition.neumann(0.0)
    lap = tensor_laplacian(values, grid, zero_flux, zero_flux, zero_flux)
    assert torch.max(torch.abs(lap)).item() < 1.0e-13


def test_jvp_probe_is_consistent_on_small_grid() -> None:
    dtype = configure_backend(BackendConfig())
    grid = RadialGrid.create(RadialGridSpec(r_max=8.0, nr=16), dtype=dtype, device="cpu")
    cfg = CubicGPEConfig()
    bc = BoundaryCondition.dirichlet(cfg.outer_boundary_value)
    state = gaussian_initial_state(grid, cfg)
    residual_fn = lambda x: cubic_gpe_residual(x, grid, cfg, bc)
    check = finite_difference_jvp_check(
        residual_fn,
        state,
        epsilon=NewtonConfig().finite_difference_jvp_epsilon,
        seed=123,
    )
    assert check["relative_residual"] < 1.0e-6
