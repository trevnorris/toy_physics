"""Conservative radial finite-volume operators."""

from __future__ import annotations

import numpy as np
import torch

from .boundaries import BoundaryCondition, ghost_value
from .grid import RadialGrid, TensorProductGrid, WallGrid


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


def tensor_weighted_laplacian(
    values: torch.Tensor,
    weight_w_centers: torch.Tensor,
    weight_w_faces: torch.Tensor,
    grid: TensorProductGrid,
    radial_outer_bc: BoundaryCondition,
    w_lower_bc: BoundaryCondition,
    w_upper_bc: BoundaryCondition,
) -> torch.Tensor:
    """Conservative ``div(weight(w) grad values)`` on the ``(r,w)`` grid."""

    if values.shape != (grid.spec.nr, grid.spec.nw):
        raise ValueError("tensor_weighted_laplacian expects shape (nr, nw)")
    if weight_w_centers.shape != (grid.spec.nw,):
        raise ValueError("weight_w_centers expects shape (nw,)")
    if weight_w_faces.shape != (grid.spec.nw + 1,):
        raise ValueError("weight_w_faces expects shape (nw + 1,)")

    nr, nw = values.shape
    radial_fluxes_rw = torch.zeros((nr + 1, nw), dtype=values.dtype, device=values.device)
    radial_fluxes_rw[1:-1, :] = (
        grid.radial_face_areas[1:-1, None]
        * grid.dw
        * weight_w_centers[None, :]
        * (values[1:, :] - values[:-1, :])
        / grid.dr
    )
    radial_fluxes_rw[0, :] = 0.0
    radial_ghost = ghost_value(values[-1, :], grid.dr, radial_outer_bc)
    radial_fluxes_rw[-1, :] = (
        grid.radial_face_areas[-1]
        * grid.dw
        * weight_w_centers
        * (radial_ghost - values[-1, :])
        / grid.dr
    )

    w_fluxes = torch.zeros((nr, nw + 1), dtype=values.dtype, device=values.device)
    lower_ghost = ghost_value(values[:, 0], grid.dw, w_lower_bc)
    upper_ghost = ghost_value(values[:, -1], grid.dw, w_upper_bc)
    w_fluxes[:, 0] = (
        grid.radial_shell_volumes
        * weight_w_faces[0]
        * (values[:, 0] - lower_ghost)
        / grid.dw
    )
    w_fluxes[:, 1:-1] = (
        grid.radial_shell_volumes[:, None]
        * weight_w_faces[1:-1][None, :]
        * (values[:, 1:] - values[:, :-1])
        / grid.dw
    )
    w_fluxes[:, -1] = (
        grid.radial_shell_volumes
        * weight_w_faces[-1]
        * (upper_ghost - values[:, -1])
        / grid.dw
    )

    divergence = (
        radial_fluxes_rw[1:, :]
        - radial_fluxes_rw[:-1, :]
        + w_fluxes[:, 1:]
        - w_fluxes[:, :-1]
    )
    return divergence / grid.cell_volumes


def radial_center_gradient(values: torch.Tensor, grid: RadialGrid) -> torch.Tensor:
    """Second-order gradient at radial cell centers."""

    if values.ndim != 1 or values.shape[0] != grid.spec.nr:
        raise ValueError("radial_center_gradient expects a one-dimensional cell field")
    if values.shape[0] < 3:
        raise ValueError("Need at least three cells for a second-order boundary gradient")
    grad = torch.empty_like(values)
    grad[1:-1] = (values[2:] - values[:-2]) / (2.0 * grid.dr)
    grad[0] = (-3.0 * values[0] + 4.0 * values[1] - values[2]) / (2.0 * grid.dr)
    grad[-1] = (3.0 * values[-1] - 4.0 * values[-2] + values[-3]) / (2.0 * grid.dr)
    return grad


def tensor_center_gradient_r(values: torch.Tensor, grid: TensorProductGrid) -> torch.Tensor:
    if values.shape != (grid.spec.nr, grid.spec.nw):
        raise ValueError("tensor_center_gradient_r expects shape (nr, nw)")
    if grid.spec.nr < 3:
        raise ValueError("Need at least three radial cells")
    grad = torch.empty_like(values)
    grad[1:-1, :] = (values[2:, :] - values[:-2, :]) / (2.0 * grid.dr)
    grad[0, :] = (-3.0 * values[0, :] + 4.0 * values[1, :] - values[2, :]) / (
        2.0 * grid.dr
    )
    grad[-1, :] = (3.0 * values[-1, :] - 4.0 * values[-2, :] + values[-3, :]) / (
        2.0 * grid.dr
    )
    return grad


def tensor_center_gradient_w(values: torch.Tensor, grid: TensorProductGrid) -> torch.Tensor:
    if values.shape != (grid.spec.nr, grid.spec.nw):
        raise ValueError("tensor_center_gradient_w expects shape (nr, nw)")
    if grid.spec.nw < 3:
        raise ValueError("Need at least three w cells")
    grad = torch.empty_like(values)
    grad[:, 1:-1] = (values[:, 2:] - values[:, :-2]) / (2.0 * grid.dw)
    grad[:, 0] = (-3.0 * values[:, 0] + 4.0 * values[:, 1] - values[:, 2]) / (
        2.0 * grid.dw
    )
    grad[:, -1] = (3.0 * values[:, -1] - 4.0 * values[:, -2] + values[:, -3]) / (
        2.0 * grid.dw
    )
    return grad


def radial_current(
    psi: torch.Tensor,
    grid: RadialGrid,
    *,
    hbar: float,
    particle_mass: float,
    gauge_charge: float = 0.0,
    gauge_potential: torch.Tensor | None = None,
) -> torch.Tensor:
    """Gauge-covariant radial current ``(hbar/m) Im(psi* D_r psi)``."""

    if psi.ndim != 1 or psi.shape[0] != grid.spec.nr:
        raise ValueError("radial_current expects a one-dimensional cell field")
    grad = radial_center_gradient(psi, grid)
    covariant = grad
    if gauge_potential is not None:
        if gauge_potential.shape != psi.shape:
            raise ValueError("gauge_potential shape must match psi")
        covariant = covariant - 1j * (gauge_charge / hbar) * gauge_potential * psi
    return (hbar / particle_mass) * torch.imag(torch.conj(psi) * covariant)


def _face_average_with_linear_extrapolation(values: torch.Tensor, axis: int) -> torch.Tensor:
    """Face values from cell centers, second-order at exterior faces."""

    if axis == 0:
        faces = torch.empty(
            (values.shape[0] + 1, values.shape[1]), dtype=values.dtype, device=values.device
        )
        faces[1:-1, :] = 0.5 * (values[:-1, :] + values[1:, :])
        faces[0, :] = 1.5 * values[0, :] - 0.5 * values[1, :]
        faces[-1, :] = 1.5 * values[-1, :] - 0.5 * values[-2, :]
        return faces
    if axis == 1:
        faces = torch.empty(
            (values.shape[0], values.shape[1] + 1), dtype=values.dtype, device=values.device
        )
        faces[:, 1:-1] = 0.5 * (values[:, :-1] + values[:, 1:])
        faces[:, 0] = 1.5 * values[:, 0] - 0.5 * values[:, 1]
        faces[:, -1] = 1.5 * values[:, -1] - 0.5 * values[:, -2]
        return faces
    raise ValueError("axis must be 0 or 1")


def axisymmetric_vector_divergence(
    radial_component: torch.Tensor,
    w_component: torch.Tensor,
    grid: TensorProductGrid,
) -> torch.Tensor:
    """Cell-average ``r^-2 d_r(r^2 A_r) + d_w A_w``."""

    expected = (grid.spec.nr, grid.spec.nw)
    if radial_component.shape != expected or w_component.shape != expected:
        raise ValueError("axisymmetric_vector_divergence expects two (nr, nw) fields")
    radial_faces = _face_average_with_linear_extrapolation(radial_component, axis=0)
    radial_faces[0, :] = 0.0
    radial_flux = grid.radial_face_areas[:, None] * grid.dw * radial_faces
    w_faces = _face_average_with_linear_extrapolation(w_component, axis=1)
    w_flux = grid.radial_shell_volumes[:, None] * w_faces
    divergence = (
        radial_flux[1:, :]
        - radial_flux[:-1, :]
        + w_flux[:, 1:]
        - w_flux[:, :-1]
    )
    return divergence / grid.cell_volumes


def radial_divergence_from_center_flux(
    radial_flux_density: torch.Tensor,
    grid: TensorProductGrid,
) -> torch.Tensor:
    """Conservative ``r^-2 d_r(r^2 flux)`` for a centered radial flux density."""

    if radial_flux_density.shape != (grid.spec.nr, grid.spec.nw):
        raise ValueError("radial_divergence_from_center_flux expects shape (nr, nw)")
    radial_faces = _face_average_with_linear_extrapolation(radial_flux_density, axis=0)
    radial_faces[0, :] = 0.0
    flux = grid.radial_face_areas[:, None] * grid.dw * radial_faces
    return (flux[1:, :] - flux[:-1, :]) / grid.cell_volumes


def localized_maxwell_operator(
    a0: torch.Tensor,
    ar: torch.Tensor,
    aw: torch.Tensor,
    grid: TensorProductGrid,
    *,
    xi: float,
    weight_w_centers: torch.Tensor,
    weight_w_faces: torch.Tensor,
    a0_radial_outer_bc: BoundaryCondition,
    a0_w_lower_bc: BoundaryCondition,
    a0_w_upper_bc: BoundaryCondition,
) -> torch.Tensor:
    """Stationary axisymmetric localized Maxwell operator with ``H=Z``.

    The retained components are ordered ``(A_0, A_r, A_w)``. The sign of the
    static scalar-potential equation follows the mostly-plus convention,
    so the ``A_0`` block is ``-div(Z grad A_0)``.
    """

    expected = (grid.spec.nr, grid.spec.nw)
    if a0.shape != expected or ar.shape != expected or aw.shape != expected:
        raise ValueError("localized_maxwell_operator expects three (nr, nw) fields")
    z = weight_w_centers[None, :]
    scalar_block = -tensor_weighted_laplacian(
        a0,
        weight_w_centers,
        weight_w_faces,
        grid,
        a0_radial_outer_bc,
        a0_w_lower_bc,
        a0_w_upper_bc,
    )
    divergence = axisymmetric_vector_divergence(ar, aw, grid)
    weighted_divergence = z * divergence
    f_rw = tensor_center_gradient_r(aw, grid) - tensor_center_gradient_w(ar, grid)
    weighted_f_rw = z * f_rw
    ar_block = -tensor_center_gradient_w(weighted_f_rw, grid) + (
        1.0 / xi
    ) * tensor_center_gradient_r(weighted_divergence, grid)
    aw_block = radial_divergence_from_center_flux(weighted_f_rw, grid) + (
        1.0 / xi
    ) * tensor_center_gradient_w(weighted_divergence, grid)
    return torch.stack([scalar_block, ar_block, aw_block], dim=0)


def wall_s_eta_operator(
    values: torch.Tensor,
    grid: WallGrid,
    *,
    t_w_faces: torch.Tensor,
    t_omega_centers: torch.Tensor,
    k_eta_centers: torch.Tensor,
    spherical_l: int,
    lower_bc: BoundaryCondition,
    upper_bc: BoundaryCondition,
) -> torch.Tensor:
    """Stationary modal wall operator from ``S_eta^(2)`` in flat ``dw`` convention."""

    if values.ndim != 1 or values.shape[0] != grid.spec.nw:
        raise ValueError("wall_s_eta_operator expects a one-dimensional wall field")
    if t_w_faces.shape != (grid.spec.nw + 1,):
        raise ValueError("t_w_faces expects shape (nw + 1,)")
    if t_omega_centers.shape != values.shape or k_eta_centers.shape != values.shape:
        raise ValueError("wall restoring coefficients must match values")

    fluxes = torch.zeros(grid.spec.nw + 1, dtype=values.dtype, device=values.device)
    lower_ghost = ghost_value(values[0], grid.dw, lower_bc)
    upper_ghost = ghost_value(values[-1], grid.dw, upper_bc)
    fluxes[0] = t_w_faces[0] * (values[0] - lower_ghost) / grid.dw
    fluxes[1:-1] = t_w_faces[1:-1] * (values[1:] - values[:-1]) / grid.dw
    fluxes[-1] = t_w_faces[-1] * (upper_ghost - values[-1]) / grid.dw
    ell_factor = float(spherical_l * (spherical_l + 1))
    restoring = (k_eta_centers + ell_factor * t_omega_centers) * values
    return -(fluxes[1:] - fluxes[:-1]) / grid.cell_widths + restoring


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
