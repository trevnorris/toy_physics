from __future__ import annotations

from dataclasses import replace
import json
import math

import torch

from stage1_solver.backend import configure_backend
from stage1_solver.boundaries import BoundaryCondition
from stage1_solver.config import (
    BackendConfig,
    CubicGPEConfig,
    HarnessConfig,
    NewtonConfig,
    RadialGridSpec,
    TensorGridSpec,
)
from stage1_solver.grid import RadialGrid, TensorProductGrid
from stage1_solver.manifest import write_manifest
from stage1_solver.newton import finite_difference_jvp_check
from stage1_solver.operators import (
    integrate,
    radial_current,
    radial_fluxes,
    radial_laplacian,
    tensor_laplacian,
)
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


def test_config_hash_excludes_cosmetic_output_paths() -> None:
    cfg = HarnessConfig()
    moved = replace(
        cfg,
        run_root="software/stage1_solver/runs/moved_elsewhere",
        report_path="software/stage1_solver/reports/moved_elsewhere.md",
    )
    changed_physics = replace(
        cfg,
        mms=replace(
            cfg.mms,
            matter=replace(cfg.mms.matter, eos_K=cfg.mms.matter.eos_K * 1.1),
        ),
    )
    assert cfg.config_hash() == moved.config_hash()
    assert cfg.config_hash() != changed_physics.config_hash()


def test_manifest_records_full_config_but_narrow_hash(tmp_path) -> None:
    cfg = HarnessConfig(run_root=str(tmp_path / "runs"), report_path=str(tmp_path / "report.md"))
    manifest = write_manifest(
        run_root=str(tmp_path / "runs"),
        benchmark_name="hash_hygiene",
        grid_name="smoke",
        config=cfg.to_dict(),
        mesh={"nr": 4},
        results={"ok": True},
    )
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["config_hash"] == cfg.config_hash()
    assert payload["config"]["run_root"] == str(tmp_path / "runs")
    assert payload["config"]["report_path"] == str(tmp_path / "report.md")


def test_complex_radial_current_is_nonzero() -> None:
    dtype = configure_backend(BackendConfig())
    grid = RadialGrid.create(RadialGridSpec(r_max=2.0, nr=128), dtype=dtype, device="cpu")
    wave_number = 0.7
    psi = torch.exp(1j * wave_number * grid.r_centers).to(torch.complex128)
    current = radial_current(psi, grid, hbar=1.0, particle_mass=2.0)
    expected = torch.full_like(grid.r_centers, 0.5 * wave_number)
    error = torch.sqrt(torch.sum((current - expected) ** 2 * grid.cell_volumes)).item()
    assert torch.max(torch.abs(current)).item() > 0.1
    assert error < 2.0e-4
