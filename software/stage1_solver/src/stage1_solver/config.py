"""Explicit configuration objects for the Stage-1 validation harness."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
import subprocess
from typing import Any


def _json_default(value: Any) -> Any:
    if hasattr(value, "tolist"):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def stable_json(data: Any) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"), default=_json_default)


def stable_hash(data: Any) -> str:
    return hashlib.sha256(stable_json(data).encode("utf-8")).hexdigest()[:16]


def source_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def config_hash_payload(config: dict[str, Any]) -> dict[str, Any]:
    payload = {
        key: value
        for key, value in config.items()
        if key not in {"run_root", "report_path"}
    }
    payload["source_revision"] = source_revision()
    return payload


def config_hash_from_dict(config: dict[str, Any]) -> str:
    return stable_hash(config_hash_payload(config))


@dataclass(frozen=True)
class BackendConfig:
    dtype: str = "float64"
    device: str = "cpu"
    seed: int = 20260612
    deterministic_algorithms: bool = True
    torch_num_threads: int = 1


@dataclass(frozen=True)
class PreconditionerConfig:
    type: str = "none"
    side: str = "left"
    rebuild_policy: str = "every_newton_step"
    stencil_radius: int = 3
    color_separation: int = 7
    factorization: str = "splu"
    diagonal_shift: float = 0.0
    drop_tolerance: float = 0.0
    fill_factor: float = 10.0
    permutation: str = "COLAMD"


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
    finite_difference_jvp_epsilon: float = 1.0e-4
    preconditioner: PreconditionerConfig = field(default_factory=PreconditionerConfig)


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
class WallGridSpec:
    w_min: float
    w_max: float
    nw: int
    measure: str = "flat dw"


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
class QuinticMatterMMSConfig:
    name: str = "mms_quintic_matter_radial"
    r_max: float = 2.0
    grid_levels: tuple[int, ...] = (32, 64, 128, 256)
    hbar: float = 0.9
    particle_mass: float = 1.2
    eos_K: float = 0.37
    trap_omega: float = 0.41
    chemical_potential: float = 0.73
    min_observed_order: float = 1.85
    final_error_max: float = 5.0e-4


@dataclass(frozen=True)
class TensorLaplacianMMSConfig:
    name: str = "mms_tensor_laplacian_2d"
    r_max: float = 2.0
    w_min: float = -0.75
    w_max: float = 0.85
    grid_levels: tuple[tuple[int, int], ...] = (
        (24, 20),
        (48, 40),
        (96, 80),
        (192, 160),
        (384, 320),
        (768, 640),
    )
    min_observed_order: float = 1.85
    final_error_max: float = 2.0e-4


@dataclass(frozen=True)
class CurrentMMSConfig:
    name: str = "mms_complex_current_radial"
    r_max: float = 2.0
    grid_levels: tuple[int, ...] = (32, 64, 128, 256, 512)
    hbar: float = 0.9
    particle_mass: float = 1.2
    gauge_charge: float = 0.0
    min_observed_order: float = 1.85
    final_error_max: float = 2.0e-5
    min_current_norm: float = 1.0e-1


@dataclass(frozen=True)
class MaxwellMMSConfig:
    name: str = "mms_localized_maxwell_h_equals_z"
    r_max: float = 2.0
    w_min: float = -0.75
    w_max: float = 0.85
    grid_levels: tuple[tuple[int, int], ...] = (
        (24, 20),
        (48, 40),
        (96, 80),
        (192, 160),
        (384, 320),
    )
    xi: float = 1.3
    mu0: float = 1.0
    localization_width: float = 0.9
    min_observed_order: float = 1.75
    final_error_max: float = 8.0e-3


@dataclass(frozen=True)
class CoupledBranchMMSConfig:
    name: str = "mms_coupled_branch_matter_maxwell"
    r_max: float = 2.0
    w_min: float = 0.0
    w_max: float = 1.6
    grid_levels: tuple[tuple[int, int], ...] = (
        (12, 10),
        (24, 20),
        (48, 40),
        (96, 80),
    )
    hbar: float = 1.0
    particle_mass: float = 1.0
    gauge_charge: float = 0.35
    mu0: float = 1.0
    xi: float = 1.0
    eos_K: float = 0.45
    chemical_potential: float = 0.9
    localization_width: float = 0.75
    localization_floor: float = 0.8
    localization_amplitude: float = 0.45
    r_mouth: float = 1.2
    r_exit: float = 0.9
    geometry_decay_length: float = 0.8
    radial_wall_strength: float = 0.65
    axial_trap_strength: float = 0.12
    min_observed_order: float = 1.75
    final_error_max: float = 1.5e-2


@dataclass(frozen=True)
class WallMMSConfig:
    name: str = "mms_wall_s_eta_2"
    w_min: float = -0.75
    w_max: float = 0.85
    grid_levels: tuple[int, ...] = (32, 64, 128, 256, 512)
    spherical_l: int = 2
    # MMS-only structural-certification placeholders; NOT the physical wall
    # packet, which is `free_choice` (prereg §E) and is frozen per-run at solve time.
    mms_only_placeholder_mu_eta: float = 1.0
    mms_only_placeholder_t_w_base: float = 1.1
    mms_only_placeholder_t_w_sine_amp: float = 0.2
    mms_only_placeholder_t_omega_base: float = 0.8
    mms_only_placeholder_t_omega_cosine_amp: float = 0.1
    mms_only_placeholder_k_eta_base: float = 0.9
    mms_only_placeholder_k_eta_bump_amp: float = 0.15
    min_observed_order: float = 1.85
    final_error_max: float = 2.0e-4


@dataclass(frozen=True)
class ManufacturedSolutionsConfig:
    matter: QuinticMatterMMSConfig = field(default_factory=QuinticMatterMMSConfig)
    tensor: TensorLaplacianMMSConfig = field(default_factory=TensorLaplacianMMSConfig)
    current: CurrentMMSConfig = field(default_factory=CurrentMMSConfig)
    maxwell: MaxwellMMSConfig = field(default_factory=MaxwellMMSConfig)
    coupled_branch: CoupledBranchMMSConfig = field(default_factory=CoupledBranchMMSConfig)
    wall: WallMMSConfig = field(default_factory=WallMMSConfig)


@dataclass(frozen=True)
class BranchSmokeConfig:
    name: str = "step3_coupled_branch_engineering_smoke"
    placeholder_label: str = (
        "engineering-smoke placeholders only; not a physical packet; target-blind"
    )
    r_max: float = 2.0
    w_min: float = 0.0
    w_max: float = 1.6
    solve_grid: tuple[int, int] = (10, 8)
    ladder_levels: tuple[tuple[int, int], ...] = (
        (8, 6),
        (10, 8),
        (12, 10),
        (16, 12),
        (20, 14),
    )
    hbar: float = 1.0
    particle_mass: float = 1.0
    gauge_charge: float = 0.35
    mu0: float = 1.0
    xi: float = 1.0
    continuation_K_values: tuple[float, ...] = (0.05, 0.15, 0.3, 0.5)
    mass: float = 1.0
    localization_width: float = 0.75
    localization_floor: float = 0.8
    localization_amplitude: float = 0.45
    r_mouth: float = 1.2
    r_exit: float = 0.9
    geometry_decay_length: float = 0.8
    radial_wall_strength: float = 0.65
    axial_trap_strength: float = 0.12
    matter_exit_impedance_alpha: float = 0.25
    a0_radial_impedance_alpha: float = 0.5
    a0_exit_impedance_alpha: float = 0.4
    initial_mu: float = 1.0
    max_ladder_level_seconds: float = 120.0
    newton: NewtonConfig = field(
        default_factory=lambda: NewtonConfig(
            residual_atol=5.0e-8,
            residual_rtol=5.0e-9,
            step_atol=1.0e-11,
            step_rtol=1.0e-10,
            max_newton_iters=16,
            gmres_rtol=5.0e-8,
            gmres_atol=1.0e-10,
            gmres_restart=96,
            gmres_maxiter=12,
            max_line_search_iters=18,
            accept_best_line_search_decrease=True,
            finite_difference_jvp_epsilon=1.0e-5,
        )
    )


@dataclass(frozen=True)
class BranchConvergenceConfig:
    name: str = "step4_coupled_branch_grid_convergence"
    refinement_ratio: int = 2
    levels: tuple[tuple[int, int], ...] = (
        (6, 4),
        (12, 8),
        (24, 16),
        (48, 32),
        (96, 64),
    )
    min_levels: int = 3
    expected_order: float = 2.0
    illustrative_error_levels: tuple[float, ...] = (1.0e-3, 1.0e-4)
    max_level_seconds: float = 180.0
    noise_floor_ratio: float = 1.25


@dataclass(frozen=True)
class HarnessConfig:
    backend: BackendConfig = field(default_factory=BackendConfig)
    newton: NewtonConfig = field(default_factory=NewtonConfig)
    linear: LinearEigenConfig = field(default_factory=LinearEigenConfig)
    cubic: CubicGPEConfig = field(default_factory=CubicGPEConfig)
    mms: ManufacturedSolutionsConfig = field(default_factory=ManufacturedSolutionsConfig)
    branch: BranchSmokeConfig = field(default_factory=BranchSmokeConfig)
    convergence: BranchConvergenceConfig = field(default_factory=BranchConvergenceConfig)
    run_root: str = "software/stage1_solver/runs/step2_manufactured_solutions"
    report_path: str = "software/stage1_solver/reports/step2_manufactured_solutions_validation.md"
    jacobian_check_seed: int = 1729
    jacobian_check_rel_tol: float = 1.0e-6
    jacobian_check_abs_tol: float = 2.0e-7

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def config_hash(self) -> str:
        return config_hash_from_dict(self.to_dict())
