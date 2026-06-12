"""Logically rectangular grid primitives.

Step 1 exercises the radial line of the later ``(r, w)`` grid. The radial
measure and face-area bookkeeping are already the production-transferable part:
weighted divergence terms are represented by conservative face fluxes.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np
import torch

from .backend import tensor
from .config import RadialGridSpec, TensorGridSpec


@dataclass(frozen=True)
class RadialGrid:
    spec: RadialGridSpec
    dtype: torch.dtype
    device: str
    r_faces: torch.Tensor
    r_centers: torch.Tensor
    face_areas: torch.Tensor
    cell_volumes: torch.Tensor
    dr: float

    @staticmethod
    def create(spec: RadialGridSpec, *, dtype: torch.dtype, device: str) -> "RadialGrid":
        if spec.r_min != 0.0:
            raise ValueError("The r=0 regularity treatment requires r_min=0 for step 1")
        if spec.nr < 4:
            raise ValueError("Need at least 4 radial cells")
        r_faces_np = np.linspace(spec.r_min, spec.r_max, spec.nr + 1, dtype=np.float64)
        dr = float(r_faces_np[1] - r_faces_np[0])
        r_centers_np = 0.5 * (r_faces_np[:-1] + r_faces_np[1:])
        face_areas_np = 4.0 * math.pi * r_faces_np**2
        cell_volumes_np = (4.0 * math.pi / 3.0) * (r_faces_np[1:] ** 3 - r_faces_np[:-1] ** 3)
        return RadialGrid(
            spec=spec,
            dtype=dtype,
            device=device,
            r_faces=tensor(r_faces_np, dtype=dtype, device=device),
            r_centers=tensor(r_centers_np, dtype=dtype, device=device),
            face_areas=tensor(face_areas_np, dtype=dtype, device=device),
            cell_volumes=tensor(cell_volumes_np, dtype=dtype, device=device),
            dr=dr,
        )

    def to_dict(self) -> dict[str, float | int | str]:
        data = asdict(self.spec)
        data["dr"] = self.dr
        return data

    def numpy_centers(self) -> np.ndarray:
        return self.r_centers.detach().cpu().numpy().copy()

    def numpy_volumes(self) -> np.ndarray:
        return self.cell_volumes.detach().cpu().numpy().copy()


@dataclass(frozen=True)
class TensorProductGrid:
    spec: TensorGridSpec
    dtype: torch.dtype
    device: str
    r_faces: torch.Tensor
    r_centers: torch.Tensor
    w_faces: torch.Tensor
    w_centers: torch.Tensor
    radial_face_areas: torch.Tensor
    radial_shell_volumes: torch.Tensor
    cell_volumes: torch.Tensor
    dr: float
    dw: float

    @staticmethod
    def create(spec: TensorGridSpec, *, dtype: torch.dtype, device: str) -> "TensorProductGrid":
        if spec.r_min != 0.0:
            raise ValueError("The r=0 regularity treatment requires r_min=0 for step 1")
        if spec.nr < 4 or spec.nw < 2:
            raise ValueError("Need at least 4 radial cells and 2 w cells")
        r_faces_np = np.linspace(spec.r_min, spec.r_max, spec.nr + 1, dtype=np.float64)
        w_faces_np = np.linspace(spec.w_min, spec.w_max, spec.nw + 1, dtype=np.float64)
        dr = float(r_faces_np[1] - r_faces_np[0])
        dw = float(w_faces_np[1] - w_faces_np[0])
        r_centers_np = 0.5 * (r_faces_np[:-1] + r_faces_np[1:])
        w_centers_np = 0.5 * (w_faces_np[:-1] + w_faces_np[1:])
        radial_face_areas_np = 4.0 * math.pi * r_faces_np**2
        radial_shell_volumes_np = (4.0 * math.pi / 3.0) * (
            r_faces_np[1:] ** 3 - r_faces_np[:-1] ** 3
        )
        cell_volumes_np = radial_shell_volumes_np[:, None] * np.full((1, spec.nw), dw)
        return TensorProductGrid(
            spec=spec,
            dtype=dtype,
            device=device,
            r_faces=tensor(r_faces_np, dtype=dtype, device=device),
            r_centers=tensor(r_centers_np, dtype=dtype, device=device),
            w_faces=tensor(w_faces_np, dtype=dtype, device=device),
            w_centers=tensor(w_centers_np, dtype=dtype, device=device),
            radial_face_areas=tensor(radial_face_areas_np, dtype=dtype, device=device),
            radial_shell_volumes=tensor(radial_shell_volumes_np, dtype=dtype, device=device),
            cell_volumes=tensor(cell_volumes_np, dtype=dtype, device=device),
            dr=dr,
            dw=dw,
        )

    def to_dict(self) -> dict[str, float | int | str]:
        data = asdict(self.spec)
        data["dr"] = self.dr
        data["dw"] = self.dw
        return data
