"""Explicit configuration objects for the Stage-1 validation harness."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
from typing import Any


def _json_default(value: Any) -> Any:
    if hasattr(value, "tolist"):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def stable_json(data: Any) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"), default=_json_default)


def stable_hash(data: Any) -> str:
    return hashlib.sha256(stable_json(data).encode("utf-8")).hexdigest()[:16]


@dataclass(frozen=True)
class BackendConfig:
    dtype: str = "float64"
    device: str = "cpu"
    seed: int = 20260612
    deterministic_algorithms: bool = True
    torch_num_threads: int = 1


@dataclass(frozen=True)
class NewtonConfig:
    residual_atol: float = 1.0e-11
    residual_rtol: float = 1.0e-10
    step_atol: float = 1.0e-13
    step_rtol: float = 1.0e-12
    max_newton_iters: int = 24
    residual_norm: str = "linf"
    linear_solver: str = "gmres_jvp"
    gmres_rtol: float = 1.0e-11
    gmres_atol: float = 1.0e-13
    gmres_restart: int = 256
    gmres_maxiter: int = 8
    line_search: str = "armijo"
    line_search_c1: float = 1.0e-4
    line_search_shrink: float = 0.5
    max_line_search_iters: int = 16
    accept_best_line_search_decrease: bool = True
    finite_difference_jvp_epsilon: float = 1.0e-6


@dataclass(frozen=True)
class RadialGridSpec:
    r_max: float
    nr: int
    r_min: float = 0.0
    measure: str = "4*pi*r^2"


@dataclass(frozen=True)
class TensorGridSpec:
    r_max: float
    nr: int
    w_min: float
    w_max: float
    nw: int
    r_min: float = 0.0
    measure: str = "4*pi*r^2 dr dw"


@dataclass(frozen=True)
class LinearEigenConfig:
    name: str = "linear_harmonic_robin"
    r_max: float = 8.0
    grid_levels: tuple[int, ...] = (128, 256, 512, 1024)
    omega: float = 1.0
    exact_mu: float = 1.5
    min_observed_order: float = 1.85
    final_eigenvalue_error_max: float = 2.0e-5
    final_field_l2_error_max: float = 1.0e-5


@dataclass(frozen=True)
class CubicGPEConfig:
    name: str = "cubic_gpe_harmonic_mass"
    r_max: float = 8.0
    grid_levels: tuple[int, ...] = (64, 128, 256)
    omega: float = 1.0
    coupling_g: float = 1.0
    mass: float = 1.0
    outer_boundary_value: float = 0.0
    initial_mu: float = 2.0
    reference_eps: float = 1.0e-5
    reference_nodes: int = 500
    reference_tol: float = 1.0e-8
    reference_bc_tol: float = 1.0e-10
    reference_max_nodes: int = 20000
    min_observed_order: float = 1.75
    final_mu_error_max: float = 3.0e-4
    final_field_l2_error_max: float = 2.0e-4
    max_mass_residual: float = 1.0e-10
    max_current_abs: float = 1.0e-12


@dataclass(frozen=True)
class HarnessConfig:
    backend: BackendConfig = field(default_factory=BackendConfig)
    newton: NewtonConfig = field(default_factory=NewtonConfig)
    linear: LinearEigenConfig = field(default_factory=LinearEigenConfig)
    cubic: CubicGPEConfig = field(default_factory=CubicGPEConfig)
    run_root: str = "software/stage1_solver/runs/step1_gpe_benchmark"
    report_path: str = "software/stage1_solver/reports/step1_gpe_benchmark_validation.md"
    jacobian_check_seed: int = 1729
    jacobian_check_rel_tol: float = 1.0e-6
    jacobian_check_abs_tol: float = 2.0e-7

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def config_hash(self) -> str:
        return stable_hash(self.to_dict())
