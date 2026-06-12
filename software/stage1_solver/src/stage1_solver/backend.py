"""Torch backend utilities kept separate from the physics residuals."""

from __future__ import annotations

import importlib.metadata
import platform
import random
from typing import Callable

import numpy as np
import scipy
import torch

from .config import BackendConfig


def resolve_dtype(dtype_name: str) -> torch.dtype:
    if dtype_name != "float64":
        raise ValueError(f"Stage-1 validation requires float64, got {dtype_name!r}")
    return torch.float64


def configure_backend(config: BackendConfig) -> torch.dtype:
    dtype = resolve_dtype(config.dtype)
    random.seed(config.seed)
    np.random.seed(config.seed)
    torch.manual_seed(config.seed)
    torch.set_default_dtype(dtype)
    torch.set_num_threads(config.torch_num_threads)
    torch.use_deterministic_algorithms(config.deterministic_algorithms)
    if config.device != "cpu" and not torch.cuda.is_available():
        raise ValueError(f"Requested device {config.device!r}, but CUDA is unavailable")
    return dtype


def tensor(data, *, dtype: torch.dtype, device: str) -> torch.Tensor:
    return torch.as_tensor(data, dtype=dtype, device=torch.device(device))


def jvp(
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x: torch.Tensor,
    direction: torch.Tensor,
) -> torch.Tensor:
    _, jv = torch.func.jvp(residual_fn, (x,), (direction,))
    return jv


def library_versions() -> dict[str, str]:
    versions = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "torch": torch.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    for package in ("sympy", "h5py", "matplotlib"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions
