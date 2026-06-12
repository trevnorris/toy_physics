"""Physics residuals built on the backend-neutral discrete operators."""

from __future__ import annotations

import math

import torch

from .boundaries import BoundaryCondition
from .config import CubicGPEConfig, LinearEigenConfig
from .grid import RadialGrid
from .operators import integrate, radial_laplacian


def harmonic_potential(r: torch.Tensor, omega: float) -> torch.Tensor:
    return 0.5 * (omega**2) * r * r


def linear_harmonic_operator(
    field: torch.Tensor,
    grid: RadialGrid,
    config: LinearEigenConfig,
    outer_bc: BoundaryCondition,
) -> torch.Tensor:
    return -0.5 * radial_laplacian(field, grid, outer_bc) + harmonic_potential(
        grid.r_centers, config.omega
    ) * field


def cubic_gpe_residual(
    state: torch.Tensor,
    grid: RadialGrid,
    config: CubicGPEConfig,
    outer_bc: BoundaryCondition,
) -> torch.Tensor:
    n = grid.spec.nr
    field = state[:n]
    mu = state[n]
    lap = radial_laplacian(field, grid, outer_bc)
    potential = harmonic_potential(grid.r_centers, config.omega)
    pde = -0.5 * lap + potential * field + config.coupling_g * field**3 - mu * field
    mass_residual = integrate(field * field, grid) - config.mass
    return torch.cat([pde, mass_residual.reshape(1)])


def gaussian_initial_state(grid: RadialGrid, config: CubicGPEConfig) -> torch.Tensor:
    """Harmonic-oscillator profile, normalized by the discrete radial measure."""

    amp = math.sqrt(config.mass) * (config.omega / math.pi) ** 0.75
    field = amp * torch.exp(-0.5 * config.omega * grid.r_centers * grid.r_centers)
    norm = torch.sqrt(integrate(field * field, grid) / config.mass)
    field = field / norm
    mu = torch.as_tensor(config.initial_mu, dtype=field.dtype, device=field.device)
    return torch.cat([field, mu.reshape(1)])
