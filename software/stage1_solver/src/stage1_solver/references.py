"""Independent references for the analytic and nonlinear benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_bvp
from scipy.linalg import eig
import torch

from .boundaries import BoundaryCondition
from .config import CubicGPEConfig, LinearEigenConfig
from .grid import RadialGrid
from .operators import radial_laplacian_matrix


def linear_harmonic_robin_bc(config: LinearEigenConfig) -> BoundaryCondition:
    return BoundaryCondition.robin(alpha=config.omega * config.r_max, beta=1.0, gamma=0.0)


def linear_harmonic_exact_field(r: np.ndarray, config: LinearEigenConfig) -> np.ndarray:
    return np.exp(-0.5 * config.omega * r * r)


def discrete_linear_harmonic_ground_state(
    grid: RadialGrid,
    config: LinearEigenConfig,
) -> tuple[float, np.ndarray, float]:
    """Solve the discrete linear eigenproblem built from the FV operator."""

    bc = linear_harmonic_robin_bc(config)
    lap = radial_laplacian_matrix(grid, bc)
    r = grid.numpy_centers()
    hamiltonian = -0.5 * lap + np.diag(0.5 * (config.omega**2) * r * r)
    eigenvalues, eigenvectors = eig(hamiltonian)
    real_mask = np.abs(eigenvalues.imag) < 1.0e-9
    if not np.any(real_mask):
        raise RuntimeError("No real eigenvalues found in linear benchmark")
    real_values = eigenvalues.real[real_mask]
    real_vectors = eigenvectors[:, real_mask].real
    index = int(np.argmin(real_values))
    mu = float(real_values[index])
    field = real_vectors[:, index]
    volumes = grid.numpy_volumes()
    norm = math.sqrt(float(np.sum(volumes * field * field)))
    field = field / norm
    exact = linear_harmonic_exact_field(r, config)
    exact = exact / math.sqrt(float(np.sum(volumes * exact * exact)))
    if float(np.sum(volumes * field * exact)) < 0.0:
        field = -field
    residual = hamiltonian @ field - mu * field
    residual_linf = float(np.max(np.abs(residual)))
    return mu, field, residual_linf


@dataclass(frozen=True)
class CubicReference:
    mu: float
    r: np.ndarray
    field: np.ndarray
    mass: float
    max_ode_residual: float
    solver_nodes: int
    solver_status: int
    solver_message: str

    def evaluate(self, r_eval: np.ndarray) -> np.ndarray:
        return np.interp(r_eval, self.r, self.field)


def solve_cubic_gpe_reference(config: CubicGPEConfig) -> CubicReference:
    """High-accuracy SciPy BVP reference independent of the torch Newton path."""

    eps = config.reference_eps
    r = np.linspace(eps, config.r_max, config.reference_nodes, dtype=np.float64)
    amp = math.sqrt(config.mass) * (config.omega / math.pi) ** 0.75
    field_guess = amp * np.exp(-0.5 * config.omega * r * r)
    mass_guess = cumulative_trapezoid(4.0 * math.pi * r * r * field_guess * field_guess, r, initial=0.0)
    if mass_guess[-1] > 0.0:
        field_guess *= math.sqrt(config.mass / mass_guess[-1])
        mass_guess = cumulative_trapezoid(
            4.0 * math.pi * r * r * field_guess * field_guess, r, initial=0.0
        )
    derivative_guess = -config.omega * r * field_guess
    y_guess = np.vstack([field_guess, derivative_guess, mass_guess])

    def ode(r_vals: np.ndarray, y: np.ndarray, p: np.ndarray) -> np.ndarray:
        mu = p[0]
        radius = np.maximum(r_vals, eps)
        potential = 0.5 * (config.omega**2) * r_vals * r_vals
        dy0 = y[1]
        dy1 = 2.0 * (potential + config.coupling_g * y[0] * y[0] - mu) * y[0]
        dy1 = dy1 - 2.0 * y[1] / radius
        dy2 = 4.0 * math.pi * r_vals * r_vals * y[0] * y[0]
        return np.vstack([dy0, dy1, dy2])

    def bc(ya: np.ndarray, yb: np.ndarray, p: np.ndarray) -> np.ndarray:
        return np.array(
            [
                ya[1],
                ya[2],
                yb[0] - config.outer_boundary_value,
                yb[2] - config.mass,
            ],
            dtype=np.float64,
        )

    solution = solve_bvp(
        ode,
        bc,
        r,
        y_guess,
        p=np.array([config.initial_mu], dtype=np.float64),
        tol=config.reference_tol,
        bc_tol=config.reference_bc_tol,
        max_nodes=config.reference_max_nodes,
        verbose=0,
    )
    if not solution.success:
        raise RuntimeError(f"Cubic reference BVP failed: {solution.message}")

    dense_r = np.linspace(eps, config.r_max, max(4000, config.reference_nodes), dtype=np.float64)
    dense_y = solution.sol(dense_r)
    mass = float(np.trapz(4.0 * math.pi * dense_r * dense_r * dense_y[0] * dense_y[0], dense_r))
    ode_res = ode(dense_r, dense_y, solution.p)
    numerical_derivative = np.vstack(
        [
            np.gradient(dense_y[0], dense_r, edge_order=2),
            np.gradient(dense_y[1], dense_r, edge_order=2),
            np.gradient(dense_y[2], dense_r, edge_order=2),
        ]
    )
    max_ode_residual = float(np.max(np.abs(numerical_derivative - ode_res)))
    return CubicReference(
        mu=float(solution.p[0]),
        r=dense_r,
        field=dense_y[0],
        mass=mass,
        max_ode_residual=max_ode_residual,
        solver_nodes=int(solution.x.size),
        solver_status=int(solution.status),
        solver_message=str(solution.message),
    )


def reference_tensor(reference: CubicReference, grid: RadialGrid) -> torch.Tensor:
    values = reference.evaluate(grid.numpy_centers())
    return torch.as_tensor(values, dtype=grid.dtype, device=grid.r_centers.device)
