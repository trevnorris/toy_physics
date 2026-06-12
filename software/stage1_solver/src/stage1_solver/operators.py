"""Conservative radial finite-volume operators."""

from __future__ import annotations

import numpy as np
import torch

from .boundaries import BoundaryCondition, ghost_value
from .grid import RadialGrid, TensorProductGrid


def radial_fluxes(
    values: torch.Tensor,
    grid: RadialGrid,
    outer_bc: BoundaryCondition,
) -> torch.Tensor:
    """Return area-weighted radial derivative fluxes on all faces."""

    if values.ndim != 1 or values.shape[0] != grid.spec.nr:
        raise ValueError("radial_fluxes expects a one-dimensional cell field")
    fluxes = torch.zeros(grid.spec.nr + 1, dtype=values.dtype, device=values.device)
    dr = grid.dr
    interior_grad = (values[1:] - values[:-1]) / dr
    fluxes[1:-1] = grid.face_areas[1:-1] * interior_grad
    fluxes[0] = 0.0
    right_ghost = ghost_value(values[-1], dr, outer_bc)
    fluxes[-1] = grid.face_areas[-1] * (right_ghost - values[-1]) / dr
    return fluxes


def radial_laplacian(
    values: torch.Tensor,
    grid: RadialGrid,
    outer_bc: BoundaryCondition,
) -> torch.Tensor:
    """Conservative approximation to ``r^-2 d_r(r^2 d_r values)``.

    The denominator is the exact shell volume, so the operator is a finite
    volume cell-average divergence. The r=0 regularity condition is the zero
    area of the inner face, giving zero radial flux without a special stencil.
    """

    fluxes = radial_fluxes(values, grid, outer_bc)
    return (fluxes[1:] - fluxes[:-1]) / grid.cell_volumes


def radial_laplacian_matrix(grid: RadialGrid, outer_bc: BoundaryCondition) -> np.ndarray:
    """Dense matrix matching :func:`radial_laplacian` for homogeneous BCs."""

    n = grid.spec.nr
    dr = grid.dr
    face_areas = grid.face_areas.detach().cpu().numpy()
    volumes = grid.cell_volumes.detach().cpu().numpy()
    mat = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        if i > 0:
            coeff = face_areas[i] / dr
            mat[i, i - 1] += coeff / volumes[i]
            mat[i, i] -= coeff / volumes[i]
        if i < n - 1:
            coeff = face_areas[i + 1] / dr
            mat[i, i + 1] += coeff / volumes[i]
            mat[i, i] -= coeff / volumes[i]
        else:
            c_ghost = homogeneous_right_ghost_coefficient(dr, outer_bc)
            coeff = face_areas[-1] / dr
            mat[i, i] += coeff * (c_ghost - 1.0) / volumes[i]
    return mat


def tensor_laplacian(
    values: torch.Tensor,
    grid: TensorProductGrid,
    radial_outer_bc: BoundaryCondition,
    w_lower_bc: BoundaryCondition,
    w_upper_bc: BoundaryCondition,
) -> torch.Tensor:
    """Conservative ``(r,w)`` Laplacian on the tensor-product grid."""

    if values.shape != (grid.spec.nr, grid.spec.nw):
        raise ValueError("tensor_laplacian expects shape (nr, nw)")

    nr, nw = values.shape
    radial_fluxes_rw = torch.zeros((nr + 1, nw), dtype=values.dtype, device=values.device)
    radial_fluxes_rw[1:-1, :] = (
        grid.radial_face_areas[1:-1, None]
        * grid.dw
        * (values[1:, :] - values[:-1, :])
        / grid.dr
    )
    radial_fluxes_rw[0, :] = 0.0
    radial_ghost = ghost_value(values[-1, :], grid.dr, radial_outer_bc)
    radial_fluxes_rw[-1, :] = (
        grid.radial_face_areas[-1] * grid.dw * (radial_ghost - values[-1, :]) / grid.dr
    )

    w_fluxes = torch.zeros((nr, nw + 1), dtype=values.dtype, device=values.device)
    lower_ghost = ghost_value(values[:, 0], grid.dw, w_lower_bc)
    upper_ghost = ghost_value(values[:, -1], grid.dw, w_upper_bc)
    w_fluxes[:, 0] = grid.radial_shell_volumes * (values[:, 0] - lower_ghost) / grid.dw
    w_fluxes[:, 1:-1] = (
        grid.radial_shell_volumes[:, None] * (values[:, 1:] - values[:, :-1]) / grid.dw
    )
    w_fluxes[:, -1] = grid.radial_shell_volumes * (upper_ghost - values[:, -1]) / grid.dw

    divergence = (
        radial_fluxes_rw[1:, :]
        - radial_fluxes_rw[:-1, :]
        + w_fluxes[:, 1:]
        - w_fluxes[:, :-1]
    )
    return divergence / grid.cell_volumes


def homogeneous_right_ghost_coefficient(spacing: float, bc: BoundaryCondition) -> float:
    if bc.kind == "dirichlet":
        if bc.gamma != 0.0:
            raise ValueError("Matrix assembly requires homogeneous Dirichlet data")
        return -1.0
    if bc.kind == "robin":
        if bc.gamma != 0.0:
            raise ValueError("Matrix assembly requires homogeneous Robin data")
        denom = 0.5 * bc.alpha + bc.beta / spacing
        return -(0.5 * bc.alpha - bc.beta / spacing) / denom
    raise ValueError(f"Unsupported boundary kind {bc.kind!r}")


def integrate(values: torch.Tensor, grid: RadialGrid) -> torch.Tensor:
    return torch.sum(values * grid.cell_volumes)


def weighted_l2_error(values: torch.Tensor, reference: torch.Tensor, grid: RadialGrid) -> torch.Tensor:
    diff = values - reference
    return torch.sqrt(torch.sum(diff * diff * grid.cell_volumes))


def max_abs_radial_current(real_field: torch.Tensor, grid: RadialGrid) -> float:
    """Current diagnostic for real stationary benchmark fields."""

    zeros = torch.zeros_like(real_field)
    psi = real_field.to(torch.complex128) + 1j * zeros.to(torch.complex128)
    grad = torch.zeros_like(psi)
    grad[1:-1] = (psi[2:] - psi[:-2]) / (2.0 * grid.dr)
    grad[0] = (psi[1] - psi[0]) / grid.dr
    grad[-1] = (psi[-1] - psi[-2]) / grid.dr
    current = torch.imag(torch.conj(psi) * grad)
    return float(torch.max(torch.abs(current)).detach().cpu().item())
