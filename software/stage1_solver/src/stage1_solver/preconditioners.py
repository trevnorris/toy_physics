"""Preconditioners for matrix-free Newton-Krylov solves.

The Step-3b coupled-branch preconditioner assembles a sparse Jacobian from the
same JVP used by GMRES.  Graph coloring keeps the number of JVP probes bounded
by the local stencil width, and the resulting sparse factorization is used only
as a Krylov preconditioner.  The residual and matrix-free JVP remain the source
of truth for the Newton solve.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.sparse import csc_matrix, coo_matrix, eye
from scipy.sparse.linalg import LinearOperator, spilu, splu
import torch

from .backend import jvp
from .config import PreconditionerConfig
from .grid import TensorProductGrid
from .newton import BuiltPreconditioner, PreconditionerBuildContext


@dataclass(frozen=True)
class CoupledSparseJacobianLayout:
    field_count: int
    cell_count: int
    state_size: int
    nr: int
    nw: int


@dataclass(frozen=True)
class P2TangentSparseJacobianLayout:
    cell_field_count: int
    cell_count: int
    wall_count: int
    state_size: int
    nr: int
    nw: int


@dataclass(frozen=True)
class ClosedCoupledSparseJacobianLayout:
    cell_field_count: int
    cell_count: int
    wall_count: int
    state_size: int
    nr: int
    nw: int


def _validate_coloring_config(config: PreconditionerConfig) -> None:
    if config.type != "colored_sparse_jacobian_lu":
        raise ValueError(f"Unsupported preconditioner type {config.type!r}")
    _validate_sparse_factorization_config(config)
    if config.stencil_radius < 0:
        raise ValueError("stencil_radius must be non-negative")
    if config.color_separation <= 2 * config.stencil_radius:
        raise ValueError("color_separation must exceed 2 * stencil_radius")


def _validate_autodiff_sparse_config(config: PreconditionerConfig) -> None:
    if config.type != "autodiff_sparse_jacobian_lu":
        raise ValueError(f"Unsupported preconditioner type {config.type!r}")
    _validate_sparse_factorization_config(config)


def _validate_sparse_factorization_config(config: PreconditionerConfig) -> None:
    if config.side != "left":
        raise ValueError(f"Unsupported preconditioner side {config.side!r}")
    if config.rebuild_policy not in {
        "every_newton_step",
        "once_per_newton_solve",
        "once_per_continuation",
    }:
        raise ValueError(f"Unsupported rebuild policy {config.rebuild_policy!r}")
    if config.factorization not in {"splu", "spilu"}:
        raise ValueError(f"Unsupported factorization {config.factorization!r}")
    if config.drop_tolerance < 0.0:
        raise ValueError("drop_tolerance must be non-negative")
    if config.fill_factor <= 0.0:
        raise ValueError("fill_factor must be positive")


def _matrix_from_dense_jacobian(
    jacobian: np.ndarray,
    *,
    drop_tolerance: float,
) -> csc_matrix:
    if jacobian.ndim != 2 or jacobian.shape[0] != jacobian.shape[1]:
        raise ValueError("Jacobian must be a square two-dimensional array")
    if drop_tolerance > 0.0:
        row, col = np.nonzero(np.abs(jacobian) >= drop_tolerance)
    else:
        row, col = np.nonzero(jacobian != 0.0)
    data = jacobian[row, col].astype(np.float64, copy=False)
    matrix = coo_matrix(
        (data, (row.astype(np.int64), col.astype(np.int64))),
        shape=jacobian.shape,
        dtype=np.float64,
    ).tocsc()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    return matrix


def _color_groups(layout: CoupledSparseJacobianLayout, separation: int) -> list[list[int]]:
    colors: dict[int, list[int]] = {}
    plane_colors = separation * separation
    for field in range(layout.field_count):
        field_offset = field * layout.cell_count
        color_offset = field * plane_colors
        for i in range(layout.nr):
            for j in range(layout.nw):
                cell = i * layout.nw + j
                color = color_offset + (i % separation) * separation + (j % separation)
                colors.setdefault(color, []).append(field_offset + cell)
    colors[layout.field_count * plane_colors] = [layout.field_count * layout.cell_count]
    return [colors[key] for key in sorted(colors)]


def _p2_tangent_color_groups(
    layout: P2TangentSparseJacobianLayout,
    separation: int,
) -> list[list[int]]:
    colors: dict[int, list[int]] = {}
    plane_colors = separation * separation
    for field in range(layout.cell_field_count):
        field_offset = field * layout.cell_count
        color_offset = field * plane_colors
        for i in range(layout.nr):
            for j in range(layout.nw):
                cell = i * layout.nw + j
                color = color_offset + (i % separation) * separation + (j % separation)
                colors.setdefault(color, []).append(field_offset + cell)
    wall_color_offset = layout.cell_field_count * plane_colors
    wall_offset = layout.cell_field_count * layout.cell_count
    for j in range(layout.wall_count):
        color = wall_color_offset + (j % separation)
        colors.setdefault(color, []).append(wall_offset + j)
    return [colors[key] for key in sorted(colors)]


def _closed_coupled_color_groups(
    layout: ClosedCoupledSparseJacobianLayout,
    separation: int,
) -> list[list[int]]:
    colors: dict[int, list[int]] = {}
    plane_colors = separation * separation
    for field in range(layout.cell_field_count):
        field_offset = field * layout.cell_count
        color_offset = field * plane_colors
        for i in range(layout.nr):
            for j in range(layout.nw):
                cell = i * layout.nw + j
                color = color_offset + (i % separation) * separation + (j % separation)
                colors.setdefault(color, []).append(field_offset + cell)
    wall_color_offset = layout.cell_field_count * plane_colors
    wall_offset = layout.cell_field_count * layout.cell_count
    for j in range(layout.wall_count):
        color = wall_color_offset + (j % separation)
        colors.setdefault(color, []).append(wall_offset + j)
    colors[wall_color_offset + separation] = [layout.state_size - 1]
    return [colors[key] for key in sorted(colors)]


def _local_residual_rows(
    *,
    layout: CoupledSparseJacobianLayout,
    cell: int,
    radius: int,
) -> np.ndarray:
    i = cell // layout.nw
    j = cell % layout.nw
    i0 = max(0, i - radius)
    i1 = min(layout.nr, i + radius + 1)
    j0 = max(0, j - radius)
    j1 = min(layout.nw, j + radius + 1)
    cells = np.array(
        [ii * layout.nw + jj for ii in range(i0, i1) for jj in range(j0, j1)],
        dtype=np.int64,
    )
    rows = [
        row_field * layout.cell_count + cells
        for row_field in range(layout.field_count)
    ]
    return np.concatenate(rows)


def _p2_cell_residual_rows(
    *,
    layout: P2TangentSparseJacobianLayout,
    cell: int,
    radius: int,
) -> np.ndarray:
    i = cell // layout.nw
    j = cell % layout.nw
    i0 = max(0, i - radius)
    i1 = min(layout.nr, i + radius + 1)
    j0 = max(0, j - radius)
    j1 = min(layout.nw, j + radius + 1)
    cells = np.array(
        [ii * layout.nw + jj for ii in range(i0, i1) for jj in range(j0, j1)],
        dtype=np.int64,
    )
    rows = [
        row_field * layout.cell_count + cells
        for row_field in range(layout.cell_field_count)
    ]
    return np.concatenate(rows)


def _p2_wall_column_rows(
    *,
    layout: P2TangentSparseJacobianLayout,
    wall_index: int,
    radius: int,
) -> np.ndarray:
    j0 = max(0, wall_index - radius)
    j1 = min(layout.nw, wall_index + radius + 1)
    cells = np.array(
        [i * layout.nw + j for i in range(layout.nr) for j in range(j0, j1)],
        dtype=np.int64,
    )
    cell_rows = [
        row_field * layout.cell_count + cells
        for row_field in range(layout.cell_field_count)
    ]
    wall_rows = (
        layout.cell_field_count * layout.cell_count
        + np.arange(j0, j1, dtype=np.int64)
    )
    return np.concatenate([*cell_rows, wall_rows])


def _closed_cell_column_rows(
    *,
    layout: ClosedCoupledSparseJacobianLayout,
    cell: int,
    radius: int,
) -> np.ndarray:
    cell_rows = _p2_cell_residual_rows(
        layout=P2TangentSparseJacobianLayout(
            cell_field_count=layout.cell_field_count,
            cell_count=layout.cell_count,
            wall_count=layout.wall_count,
            state_size=layout.cell_field_count * layout.cell_count + layout.wall_count,
            nr=layout.nr,
            nw=layout.nw,
        ),
        cell=cell,
        radius=radius,
    )
    j = cell % layout.nw
    j0 = max(0, j - radius)
    j1 = min(layout.nw, j + radius + 1)
    wall_offset = layout.cell_field_count * layout.cell_count
    wall_rows = wall_offset + np.arange(j0, j1, dtype=np.int64)
    return np.concatenate([cell_rows, wall_rows])


def _closed_wall_column_rows(
    *,
    layout: ClosedCoupledSparseJacobianLayout,
    wall_index: int,
    radius: int,
) -> np.ndarray:
    j0 = max(0, wall_index - radius)
    j1 = min(layout.nw, wall_index + radius + 1)
    cells = np.array(
        [i * layout.nw + j for i in range(layout.nr) for j in range(j0, j1)],
        dtype=np.int64,
    )
    cell_rows = [
        row_field * layout.cell_count + cells
        for row_field in range(layout.cell_field_count)
    ]
    wall_offset = layout.cell_field_count * layout.cell_count
    wall_rows = wall_offset + np.arange(j0, j1, dtype=np.int64)
    return np.concatenate([*cell_rows, wall_rows])


def _constraint_row_entries(
    x: torch.Tensor,
    grid: TensorProductGrid,
    layout: CoupledSparseJacobianLayout,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    volumes = grid.cell_volumes.detach().cpu().numpy().reshape(-1).astype(np.float64)
    state = x.detach().cpu().numpy().astype(np.float64, copy=False)
    psi_real = state[: layout.cell_count]
    psi_imag = state[layout.cell_count : 2 * layout.cell_count]
    row = np.full(2 * layout.cell_count, layout.state_size - 1, dtype=np.int64)
    col = np.concatenate(
        [
            np.arange(layout.cell_count, dtype=np.int64),
            np.arange(layout.cell_count, 2 * layout.cell_count, dtype=np.int64),
        ]
    )
    data = np.concatenate([2.0 * psi_real * volumes, 2.0 * psi_imag * volumes])
    return row, col, data


def assemble_coupled_colored_sparse_jacobian(
    context: PreconditionerBuildContext,
    grid: TensorProductGrid,
) -> tuple[csc_matrix, dict[str, Any]]:
    config = context.config.preconditioner
    _validate_coloring_config(config)
    cell_count = grid.spec.nr * grid.spec.nw
    layout = CoupledSparseJacobianLayout(
        field_count=5,
        cell_count=cell_count,
        state_size=5 * cell_count + 1,
        nr=grid.spec.nr,
        nw=grid.spec.nw,
    )
    if context.x.numel() != layout.state_size:
        raise ValueError(
            f"Expected coupled state with {layout.state_size} entries, got {context.x.numel()}"
        )

    color_groups = _color_groups(layout, config.color_separation)
    row_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []
    jvp_count = 0

    for columns in color_groups:
        direction_np = np.zeros(layout.state_size, dtype=np.float64)
        direction_np[columns] = 1.0
        direction = torch.as_tensor(direction_np, dtype=context.x.dtype, device=context.x.device)
        probe = (
            jvp(context.residual_fn, context.x, direction)
            .detach()
            .cpu()
            .numpy()
            .astype(np.float64, copy=False)
        )
        jvp_count += 1
        for column in columns:
            if column == layout.state_size - 1:
                rows = np.arange(0, 2 * layout.cell_count, dtype=np.int64)
            else:
                cell = column % layout.cell_count
                rows = _local_residual_rows(
                    layout=layout,
                    cell=cell,
                    radius=config.stencil_radius,
                )
            row_chunks.append(rows)
            col_chunks.append(np.full(rows.shape, column, dtype=np.int64))
            data_chunks.append(probe[rows])

    constraint_rows, constraint_cols, constraint_data = _constraint_row_entries(
        context.x,
        grid,
        layout,
    )
    row_chunks.append(constraint_rows)
    col_chunks.append(constraint_cols)
    data_chunks.append(constraint_data)

    rows = np.concatenate(row_chunks)
    cols = np.concatenate(col_chunks)
    data = np.concatenate(data_chunks)
    if config.drop_tolerance > 0.0:
        keep = np.abs(data) >= config.drop_tolerance
        rows = rows[keep]
        cols = cols[keep]
        data = data[keep]
    matrix = coo_matrix(
        (data, (rows, cols)),
        shape=(layout.state_size, layout.state_size),
        dtype=np.float64,
    ).tocsc()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    if config.diagonal_shift != 0.0:
        matrix = matrix + config.diagonal_shift * eye(
            layout.state_size,
            format="csc",
            dtype=np.float64,
        )
    metadata = {
        "type": config.type,
        "side": config.side,
        "rebuild_policy": config.rebuild_policy,
        "stencil_radius": config.stencil_radius,
        "color_separation": config.color_separation,
        "active_color_count": len(color_groups),
        "jvp_count": jvp_count,
        "state_size": layout.state_size,
        "matrix_nnz": int(matrix.nnz),
        "matrix_density": float(matrix.nnz / (layout.state_size * layout.state_size)),
        "factorization": config.factorization,
        "diagonal_shift": config.diagonal_shift,
        "drop_tolerance": config.drop_tolerance,
        "fill_factor": config.fill_factor,
        "permutation": config.permutation,
    }
    return matrix, metadata


def assemble_p2_tangent_colored_sparse_jacobian(
    context: PreconditionerBuildContext,
    grid: TensorProductGrid,
) -> tuple[csc_matrix, dict[str, Any]]:
    config = context.config.preconditioner
    _validate_coloring_config(config)
    cell_count = grid.spec.nr * grid.spec.nw
    layout = P2TangentSparseJacobianLayout(
        cell_field_count=5,
        cell_count=cell_count,
        wall_count=grid.spec.nw,
        state_size=5 * cell_count + grid.spec.nw,
        nr=grid.spec.nr,
        nw=grid.spec.nw,
    )
    if context.x.numel() != layout.state_size:
        raise ValueError(
            f"Expected P2 tangent state with {layout.state_size} entries, got {context.x.numel()}"
        )

    color_groups = _p2_tangent_color_groups(layout, config.color_separation)
    row_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []
    jvp_count = 0
    wall_offset = layout.cell_field_count * layout.cell_count

    for columns in color_groups:
        direction_np = np.zeros(layout.state_size, dtype=np.float64)
        direction_np[columns] = 1.0
        direction = torch.as_tensor(direction_np, dtype=context.x.dtype, device=context.x.device)
        probe = (
            jvp(context.residual_fn, context.x, direction)
            .detach()
            .cpu()
            .numpy()
            .astype(np.float64, copy=False)
        )
        jvp_count += 1
        for column in columns:
            if column < wall_offset:
                cell = column % layout.cell_count
                rows = _p2_cell_residual_rows(
                    layout=layout,
                    cell=cell,
                    radius=config.stencil_radius,
                )
            else:
                rows = _p2_wall_column_rows(
                    layout=layout,
                    wall_index=column - wall_offset,
                    radius=config.stencil_radius,
                )
            row_chunks.append(rows)
            col_chunks.append(np.full(rows.shape, column, dtype=np.int64))
            data_chunks.append(probe[rows])

    rows = np.concatenate(row_chunks)
    cols = np.concatenate(col_chunks)
    data = np.concatenate(data_chunks)
    if config.drop_tolerance > 0.0:
        keep = np.abs(data) >= config.drop_tolerance
        rows = rows[keep]
        cols = cols[keep]
        data = data[keep]
    matrix = coo_matrix(
        (data, (rows, cols)),
        shape=(layout.state_size, layout.state_size),
        dtype=np.float64,
    ).tocsc()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    if config.diagonal_shift != 0.0:
        matrix = matrix + config.diagonal_shift * eye(
            layout.state_size,
            format="csc",
            dtype=np.float64,
        )
    metadata = {
        "type": config.type,
        "side": config.side,
        "rebuild_policy": config.rebuild_policy,
        "stencil_radius": config.stencil_radius,
        "color_separation": config.color_separation,
        "active_color_count": len(color_groups),
        "jvp_count": jvp_count,
        "state_size": layout.state_size,
        "cell_field_count": layout.cell_field_count,
        "wall_count": layout.wall_count,
        "matrix_nnz": int(matrix.nnz),
        "matrix_density": float(matrix.nnz / (layout.state_size * layout.state_size)),
        "factorization": config.factorization,
        "diagonal_shift": config.diagonal_shift,
        "drop_tolerance": config.drop_tolerance,
        "fill_factor": config.fill_factor,
        "permutation": config.permutation,
    }
    return matrix, metadata


def assemble_closed_coupled_colored_sparse_jacobian(
    context: PreconditionerBuildContext,
    grid: TensorProductGrid,
) -> tuple[csc_matrix, dict[str, Any]]:
    config = context.config.preconditioner
    _validate_coloring_config(config)
    cell_count = grid.spec.nr * grid.spec.nw
    layout = ClosedCoupledSparseJacobianLayout(
        cell_field_count=5,
        cell_count=cell_count,
        wall_count=grid.spec.nw,
        state_size=5 * cell_count + grid.spec.nw + 1,
        nr=grid.spec.nr,
        nw=grid.spec.nw,
    )
    if context.x.numel() != layout.state_size:
        raise ValueError(
            f"Expected closed coupled state with {layout.state_size} entries, "
            f"got {context.x.numel()}"
        )

    color_groups = _closed_coupled_color_groups(layout, config.color_separation)
    row_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []
    jvp_count = 0
    wall_offset = layout.cell_field_count * layout.cell_count
    mass_column = layout.state_size - 1

    for columns in color_groups:
        direction_np = np.zeros(layout.state_size, dtype=np.float64)
        direction_np[columns] = 1.0
        direction = torch.as_tensor(direction_np, dtype=context.x.dtype, device=context.x.device)
        probe = (
            jvp(context.residual_fn, context.x, direction)
            .detach()
            .cpu()
            .numpy()
            .astype(np.float64, copy=False)
        )
        jvp_count += 1
        for column in columns:
            if column < wall_offset:
                cell = column % layout.cell_count
                rows = _closed_cell_column_rows(
                    layout=layout,
                    cell=cell,
                    radius=config.stencil_radius,
                )
            elif column < mass_column:
                rows = _closed_wall_column_rows(
                    layout=layout,
                    wall_index=column - wall_offset,
                    radius=config.stencil_radius,
                )
            else:
                rows = np.arange(0, 2 * layout.cell_count, dtype=np.int64)
            row_chunks.append(rows)
            col_chunks.append(np.full(rows.shape, column, dtype=np.int64))
            data_chunks.append(probe[rows])

    constraint_rows, constraint_cols, constraint_data = _constraint_row_entries(
        context.x,
        grid,
        CoupledSparseJacobianLayout(
            field_count=layout.cell_field_count,
            cell_count=layout.cell_count,
            state_size=layout.state_size,
            nr=layout.nr,
            nw=layout.nw,
        ),
    )
    row_chunks.append(constraint_rows)
    col_chunks.append(constraint_cols)
    data_chunks.append(constraint_data)

    rows = np.concatenate(row_chunks)
    cols = np.concatenate(col_chunks)
    data = np.concatenate(data_chunks)
    if config.drop_tolerance > 0.0:
        keep = np.abs(data) >= config.drop_tolerance
        rows = rows[keep]
        cols = cols[keep]
        data = data[keep]
    matrix = coo_matrix(
        (data, (rows, cols)),
        shape=(layout.state_size, layout.state_size),
        dtype=np.float64,
    ).tocsc()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    if config.diagonal_shift != 0.0:
        matrix = matrix + config.diagonal_shift * eye(
            layout.state_size,
            format="csc",
            dtype=np.float64,
        )
    metadata = {
        "type": config.type,
        "side": config.side,
        "rebuild_policy": config.rebuild_policy,
        "stencil_radius": config.stencil_radius,
        "color_separation": config.color_separation,
        "active_color_count": len(color_groups),
        "jvp_count": jvp_count,
        "state_size": layout.state_size,
        "cell_field_count": layout.cell_field_count,
        "wall_count": layout.wall_count,
        "matrix_nnz": int(matrix.nnz),
        "matrix_density": float(matrix.nnz / (layout.state_size * layout.state_size)),
        "factorization": config.factorization,
        "diagonal_shift": config.diagonal_shift,
        "drop_tolerance": config.drop_tolerance,
        "fill_factor": config.fill_factor,
        "permutation": config.permutation,
        "layout": "5*cells+nw+1",
    }
    return matrix, metadata


def assemble_closed_coupled_autodiff_sparse_jacobian(
    context: PreconditionerBuildContext,
    grid: TensorProductGrid,
    *,
    autodiff_mode: str = "jacfwd",
) -> tuple[csc_matrix, dict[str, Any]]:
    """Assemble the closed residual Jacobian by exact torch autodiff.

    This is a sparse storage path for the unchanged residual.  It materializes
    the dense Jacobian once with torch's batched AD, then keeps only nonzero
    entries in CSC form.  The dense wall/mass and mu-column couplings are
    therefore included without graph-color aliasing.
    """

    config = context.config.preconditioner
    _validate_autodiff_sparse_config(config)
    cell_count = grid.spec.nr * grid.spec.nw
    layout = ClosedCoupledSparseJacobianLayout(
        cell_field_count=5,
        cell_count=cell_count,
        wall_count=grid.spec.nw,
        state_size=5 * cell_count + grid.spec.nw + 1,
        nr=grid.spec.nr,
        nw=grid.spec.nw,
    )
    if context.x.numel() != layout.state_size:
        raise ValueError(
            f"Expected closed coupled state with {layout.state_size} entries, "
            f"got {context.x.numel()}"
        )
    x = context.x.detach()
    if autodiff_mode == "jacfwd":
        jacobian_t = torch.func.jacfwd(context.residual_fn)(x)
    elif autodiff_mode == "jacrev":
        jacobian_t = torch.func.jacrev(context.residual_fn)(x)
    else:
        raise ValueError(f"Unsupported autodiff mode {autodiff_mode!r}")
    if tuple(jacobian_t.shape) != (layout.state_size, layout.state_size):
        raise ValueError(
            "Closed residual Jacobian has unexpected shape "
            f"{tuple(jacobian_t.shape)}, expected {(layout.state_size, layout.state_size)}"
        )
    jacobian = jacobian_t.detach().cpu().numpy().astype(np.float64, copy=False)
    matrix = _matrix_from_dense_jacobian(
        jacobian,
        drop_tolerance=config.drop_tolerance,
    )
    if config.diagonal_shift != 0.0:
        matrix = matrix + config.diagonal_shift * eye(
            layout.state_size,
            format="csc",
            dtype=np.float64,
        )
    metadata = {
        "type": config.type,
        "side": config.side,
        "rebuild_policy": config.rebuild_policy,
        "assembly_method": "torch_func_dense_autodiff_to_sparse_csc",
        "autodiff_mode": autodiff_mode,
        "stencil_radius": None,
        "color_separation": None,
        "active_color_count": 0,
        "jvp_count": 0,
        "state_size": layout.state_size,
        "cell_field_count": layout.cell_field_count,
        "wall_count": layout.wall_count,
        "matrix_nnz": int(matrix.nnz),
        "matrix_density": float(matrix.nnz / (layout.state_size * layout.state_size)),
        "factorization": config.factorization,
        "diagonal_shift": config.diagonal_shift,
        "drop_tolerance": config.drop_tolerance,
        "fill_factor": config.fill_factor,
        "permutation": config.permutation,
        "layout": "5*cells+nw+1",
        "dense_wall_mass_mu_couplings_included": True,
        "residual_equivalent_to": "context.residual_fn",
    }
    return matrix, metadata


def assemble_closed_coupled_sparse_jacobian(
    context: PreconditionerBuildContext,
    grid: TensorProductGrid,
) -> tuple[csc_matrix, dict[str, Any]]:
    if context.config.preconditioner.type == "colored_sparse_jacobian_lu":
        return assemble_closed_coupled_colored_sparse_jacobian(context, grid)
    if context.config.preconditioner.type == "autodiff_sparse_jacobian_lu":
        return assemble_closed_coupled_autodiff_sparse_jacobian(context, grid)
    raise ValueError(f"Unsupported preconditioner type {context.config.preconditioner.type!r}")


def factorized_sparse_inverse_operator(
    matrix: csc_matrix,
    config: PreconditionerConfig,
) -> tuple[LinearOperator, dict[str, Any]]:
    if config.factorization == "splu":
        factor = splu(matrix, permc_spec=config.permutation)
    elif config.factorization == "spilu":
        factor = spilu(
            matrix,
            drop_tol=config.drop_tolerance,
            fill_factor=config.fill_factor,
            permc_spec=config.permutation,
        )
    else:
        raise ValueError(f"Unsupported factorization {config.factorization!r}")

    def solve(vector: np.ndarray) -> np.ndarray:
        return factor.solve(np.asarray(vector, dtype=np.float64))

    operator = LinearOperator(matrix.shape, matvec=solve, dtype=np.float64)
    metadata = {
        "factor_nnz_l": int(factor.L.nnz),
        "factor_nnz_u": int(factor.U.nnz),
    }
    return operator, metadata


def make_coupled_colored_sparse_jacobian_lu_factory(
    grid: TensorProductGrid,
) -> Callable[[PreconditionerBuildContext], BuiltPreconditioner]:
    cached: BuiltPreconditioner | None = None

    def factory(context: PreconditionerBuildContext) -> BuiltPreconditioner:
        nonlocal cached
        config = context.config.preconditioner
        _validate_coloring_config(config)
        if config.rebuild_policy != "every_newton_step" and cached is not None:
            metadata = dict(cached.metadata)
            metadata["reused"] = True
            return BuiltPreconditioner(operator=cached.operator, metadata=metadata)
        matrix, metadata = assemble_coupled_colored_sparse_jacobian(context, grid)
        operator, factor_metadata = factorized_sparse_inverse_operator(matrix, config)
        metadata.update(factor_metadata)
        metadata["reused"] = False
        built = BuiltPreconditioner(operator=operator, metadata=metadata)
        if config.rebuild_policy != "every_newton_step":
            cached = built
        return built

    return factory


def make_p2_tangent_colored_sparse_jacobian_lu_factory(
    grid: TensorProductGrid,
) -> Callable[[PreconditionerBuildContext], BuiltPreconditioner]:
    cached: BuiltPreconditioner | None = None

    def factory(context: PreconditionerBuildContext) -> BuiltPreconditioner:
        nonlocal cached
        config = context.config.preconditioner
        _validate_coloring_config(config)
        if config.rebuild_policy != "every_newton_step" and cached is not None:
            metadata = dict(cached.metadata)
            metadata["reused"] = True
            return BuiltPreconditioner(operator=cached.operator, metadata=metadata)
        matrix, metadata = assemble_p2_tangent_colored_sparse_jacobian(context, grid)
        operator, factor_metadata = factorized_sparse_inverse_operator(matrix, config)
        metadata.update(factor_metadata)
        built = BuiltPreconditioner(operator=operator, metadata=metadata)
        cached = built
        return built

    return factory


def make_closed_coupled_colored_sparse_jacobian_lu_factory(
    grid: TensorProductGrid,
) -> Callable[[PreconditionerBuildContext], BuiltPreconditioner]:
    cached: BuiltPreconditioner | None = None

    def factory(context: PreconditionerBuildContext) -> BuiltPreconditioner:
        nonlocal cached
        config = context.config.preconditioner
        _validate_coloring_config(config)
        if config.rebuild_policy != "every_newton_step" and cached is not None:
            metadata = dict(cached.metadata)
            metadata["reused"] = True
            return BuiltPreconditioner(operator=cached.operator, metadata=metadata)
        matrix, metadata = assemble_closed_coupled_colored_sparse_jacobian(context, grid)
        operator, factor_metadata = factorized_sparse_inverse_operator(matrix, config)
        metadata.update(factor_metadata)
        metadata["reused"] = False
        built = BuiltPreconditioner(operator=operator, metadata=metadata)
        if config.rebuild_policy != "every_newton_step":
            cached = built
        return built

    return factory


def make_closed_coupled_autodiff_sparse_jacobian_lu_factory(
    grid: TensorProductGrid,
) -> Callable[[PreconditionerBuildContext], BuiltPreconditioner]:
    cached: BuiltPreconditioner | None = None

    def factory(context: PreconditionerBuildContext) -> BuiltPreconditioner:
        nonlocal cached
        config = context.config.preconditioner
        _validate_autodiff_sparse_config(config)
        if config.rebuild_policy != "every_newton_step" and cached is not None:
            metadata = dict(cached.metadata)
            metadata["reused"] = True
            return BuiltPreconditioner(operator=cached.operator, metadata=metadata)
        matrix, metadata = assemble_closed_coupled_autodiff_sparse_jacobian(context, grid)
        operator, factor_metadata = factorized_sparse_inverse_operator(matrix, config)
        metadata.update(factor_metadata)
        metadata["reused"] = False
        built = BuiltPreconditioner(operator=operator, metadata=metadata)
        if config.rebuild_policy != "every_newton_step":
            cached = built
        return built

    return factory


def make_closed_coupled_sparse_jacobian_lu_factory(
    grid: TensorProductGrid,
) -> Callable[[PreconditionerBuildContext], BuiltPreconditioner]:
    cached: BuiltPreconditioner | None = None

    def factory(context: PreconditionerBuildContext) -> BuiltPreconditioner:
        nonlocal cached
        config = context.config.preconditioner
        _validate_sparse_factorization_config(config)
        if config.type not in {
            "colored_sparse_jacobian_lu",
            "autodiff_sparse_jacobian_lu",
        }:
            raise ValueError(f"Unsupported preconditioner type {config.type!r}")
        if config.rebuild_policy != "every_newton_step" and cached is not None:
            metadata = dict(cached.metadata)
            metadata["reused"] = True
            return BuiltPreconditioner(operator=cached.operator, metadata=metadata)
        matrix, metadata = assemble_closed_coupled_sparse_jacobian(context, grid)
        operator, factor_metadata = factorized_sparse_inverse_operator(matrix, config)
        metadata.update(factor_metadata)
        metadata["reused"] = False
        built = BuiltPreconditioner(operator=operator, metadata=metadata)
        if config.rebuild_policy != "every_newton_step":
            cached = built
        return built

    return factory
