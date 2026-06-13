from __future__ import annotations

from dataclasses import asdict, replace
import json
import math

import torch

from stage1_solver.boundary_characterization import (
    BoundaryCharacterizationConfig,
    InteriorWindow,
    impedance_sweep_branches,
    interior_solution_difference,
    truncation_sweep_branches,
)
from stage1_solver.backend import configure_backend
from stage1_solver.backend import jvp as torch_jvp
from stage1_solver.boundaries import BoundaryCondition
from stage1_solver.config import (
    BackendConfig,
    BranchSmokeConfig,
    CubicGPEConfig,
    HarnessConfig,
    NewtonConfig,
    PreconditionerConfig,
    RadialGridSpec,
    TensorGridSpec,
)
from stage1_solver.convergence import (
    _volume_restrict_2d,
    observed_order_from_three,
    richardson_estimate,
    validate_refinement_ladder,
)
from stage1_solver.conservation_diagnostics import (
    ConservationRegion,
    _sponge_sink_densities,
    independent_gauss_face_fluxes,
    null_floor_label,
)
from stage1_solver.coupled_branch import (
    CoupledFields,
    boundary_sponge_profile_torch,
    branch_boundary_conditions,
    coupled_branch_residual,
    coupled_pde_residual,
    initial_branch_state,
    pack_coupled_fields,
)
from stage1_solver.grid import RadialGrid, TensorProductGrid
from stage1_solver.manifest import write_manifest
from stage1_solver.newton import PreconditionerBuildContext, finite_difference_jvp_check
from stage1_solver.operators import (
    integrate,
    radial_current,
    radial_fluxes,
    radial_laplacian,
    tensor_flux_divergence,
    tensor_laplacian,
    tensor_vector_face_fluxes,
)
from stage1_solver.preconditioners import assemble_coupled_colored_sparse_jacobian
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


def test_conservative_tensor_divergence_theorem_on_hand_field() -> None:
    dtype = configure_backend(BackendConfig())
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=2.0, nr=10, w_min=0.0, w_max=1.2, nw=6),
        dtype=dtype,
        device="cpu",
    )
    rr = grid.r_centers[:, None]
    ww = grid.w_centers[None, :]
    radial = 0.2 + 0.03 * rr + 0.07 * ww + 0.01 * rr * ww
    axial = -0.1 + 0.05 * rr**2 - 0.02 * ww
    radial_flux, w_flux = tensor_vector_face_fluxes(radial, axial, grid)
    divergence = tensor_flux_divergence(radial_flux, w_flux, grid)

    r0, r1 = 2, 8
    w0, w1 = 1, 5
    volume_integral = torch.sum(divergence[r0:r1, w0:w1] * grid.cell_volumes[r0:r1, w0:w1])
    surface_flux = (
        torch.sum(radial_flux[r1, w0:w1])
        - torch.sum(radial_flux[r0, w0:w1])
        + torch.sum(w_flux[r0:r1, w1])
        - torch.sum(w_flux[r0:r1, w0])
    )
    assert abs((volume_integral - surface_flux).item()) < 1.0e-13


def test_gauss_closure_recovers_known_quadratic_charge_field() -> None:
    dtype = configure_backend(BackendConfig())
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=2.0, nr=12, w_min=0.0, w_max=1.0, nw=6),
        dtype=dtype,
        device="cpu",
    )
    source = 0.7
    a0 = -(source / 6.0) * grid.r_centers[:, None] ** 2 + torch.zeros(
        (grid.spec.nr, grid.spec.nw),
        dtype=dtype,
    )
    cfg = BranchSmokeConfig(
        r_max=2.0,
        w_min=0.0,
        w_max=1.0,
        solve_grid=(12, 6),
        localization_floor=1.0,
        localization_amplitude=0.0,
    )
    gauss_radial_flux, gauss_w_flux = independent_gauss_face_fluxes(a0, grid, cfg)
    region = ConservationRegion(r_start=0, r_stop=7, w_start=1, w_stop=5, label="known")
    lhs = (
        torch.sum(gauss_radial_flux[region.r_stop, region.w_start : region.w_stop])
        - torch.sum(gauss_radial_flux[region.r_start, region.w_start : region.w_stop])
        + torch.sum(gauss_w_flux[region.r_start : region.r_stop, region.w_stop])
        - torch.sum(gauss_w_flux[region.r_start : region.r_stop, region.w_start])
    )
    r0 = grid.r_faces[region.r_start].item()
    r1 = grid.r_faces[region.r_stop].item()
    w0 = grid.w_faces[region.w_start].item()
    w1 = grid.w_faces[region.w_stop].item()
    analytic_volume = (4.0 * math.pi / 3.0) * (r1**3 - r0**3) * (w1 - w0)
    rhs = source * analytic_volume
    assert abs((lhs - rhs).item()) < 1.0e-13


def test_conservative_restriction_preserves_constant_on_nested_tensor_grid() -> None:
    dtype = configure_backend(BackendConfig())
    coarse = TensorProductGrid.create(
        TensorGridSpec(r_max=2.0, nr=6, w_min=0.0, w_max=1.0, nw=4),
        dtype=dtype,
        device="cpu",
    )
    fine = TensorProductGrid.create(
        TensorGridSpec(r_max=2.0, nr=12, w_min=0.0, w_max=1.0, nw=8),
        dtype=dtype,
        device="cpu",
    )
    restricted = _volume_restrict_2d(torch.ones((12, 8), dtype=dtype), fine, coarse)
    assert torch.max(torch.abs(restricted - 1.0)).item() < 1.0e-14


def test_refinement_ladder_requires_fixed_ratio() -> None:
    validate_refinement_ladder(((6, 4), (12, 8), (24, 16)), 2)
    try:
        validate_refinement_ladder(((6, 4), (12, 8), (30, 16)), 2)
    except ValueError as exc:
        assert "fixed refinement ratio" in str(exc)
    else:
        raise AssertionError("non-ratio ladder should fail validation")


def test_richardson_order_and_extrapolate_second_order_sequence() -> None:
    continuum = 3.0
    c = 0.5
    values = [continuum + c * h * h for h in (0.4, 0.2, 0.1)]
    order = observed_order_from_three(values[0], values[1], values[2], 2)
    assert order is not None
    assert abs(order - 2.0) < 1.0e-12
    estimate = richardson_estimate(values[1], values[2], 2, order)
    assert estimate is not None
    assert abs(estimate - continuum) < 1.0e-12


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


def test_coupled_branch_jvp_probe_is_consistent_on_tiny_grid() -> None:
    dtype = configure_backend(BackendConfig())
    cfg = BranchSmokeConfig(solve_grid=(6, 4), ladder_levels=((6, 4),))
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.r_max, nr=6, w_min=cfg.w_min, w_max=cfg.w_max, nw=4),
        dtype=dtype,
        device="cpu",
    )
    state = initial_branch_state(grid, cfg, dtype=dtype, device="cpu")
    boundaries = branch_boundary_conditions(cfg)
    residual_fn = lambda x: coupled_branch_residual(
        x,
        grid,
        cfg,
        eos_K=cfg.continuation_K_values[0],
        boundaries=boundaries,
    )
    check = finite_difference_jvp_check(
        residual_fn,
        state,
        epsilon=cfg.newton.finite_difference_jvp_epsilon,
        seed=321,
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


def test_coupled_branch_residual_shape_and_mass_constraint() -> None:
    dtype = configure_backend(BackendConfig())
    cfg = BranchSmokeConfig(solve_grid=(6, 4), ladder_levels=((6, 4),))
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.r_max, nr=6, w_min=cfg.w_min, w_max=cfg.w_max, nw=4),
        dtype=dtype,
        device="cpu",
    )
    state = initial_branch_state(grid, cfg, dtype=dtype, device="cpu")
    residual = coupled_branch_residual(
        state,
        grid,
        cfg,
        eos_K=cfg.continuation_K_values[0],
        boundaries=branch_boundary_conditions(cfg),
    )
    assert residual.shape == state.shape
    assert abs(residual[-1].item()) < 1.0e-13


def test_colored_sparse_jacobian_matches_coupled_jvp_on_tiny_grid() -> None:
    dtype = configure_backend(BackendConfig())
    preconditioner = PreconditionerConfig(
        type="colored_sparse_jacobian_lu",
        side="left",
        rebuild_policy="every_newton_step",
        stencil_radius=3,
        color_separation=7,
        factorization="splu",
    )
    base_cfg = BranchSmokeConfig(solve_grid=(4, 3), ladder_levels=((4, 3),))
    cfg = replace(base_cfg, newton=replace(base_cfg.newton, preconditioner=preconditioner))
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.r_max, nr=4, w_min=cfg.w_min, w_max=cfg.w_max, nw=3),
        dtype=dtype,
        device="cpu",
    )
    state = initial_branch_state(grid, cfg, dtype=dtype, device="cpu")
    boundaries = branch_boundary_conditions(cfg)
    residual_fn = lambda x: coupled_branch_residual(
        x,
        grid,
        cfg,
        eos_K=cfg.continuation_K_values[0],
        boundaries=boundaries,
    )
    matrix, metadata = assemble_coupled_colored_sparse_jacobian(
        PreconditionerBuildContext(
            residual_fn=residual_fn,
            x=state,
            rhs=torch.zeros_like(state).detach().cpu().numpy(),
            iteration=1,
            config=cfg.newton,
        ),
        grid,
    )
    generator = torch.Generator(device="cpu")
    generator.manual_seed(777)
    direction = torch.randn(state.shape, dtype=dtype, generator=generator)
    direct = torch_jvp(residual_fn, state, direction).detach().cpu().numpy()
    assembled = matrix @ direction.detach().cpu().numpy()
    relative = math.sqrt(float(((assembled - direct) ** 2).sum())) / max(
        1.0,
        math.sqrt(float((direct**2).sum())),
    )
    assert metadata["active_color_count"] == state.numel()
    assert relative < 1.0e-10


def test_colored_sparse_jacobian_matches_coupled_jvp_with_compression() -> None:
    dtype = configure_backend(BackendConfig())
    preconditioner = PreconditionerConfig(
        type="colored_sparse_jacobian_lu",
        side="left",
        rebuild_policy="every_newton_step",
        stencil_radius=3,
        color_separation=7,
        factorization="splu",
    )
    base_cfg = BranchSmokeConfig(solve_grid=(14, 14), ladder_levels=((14, 14),))
    cfg = replace(base_cfg, newton=replace(base_cfg.newton, preconditioner=preconditioner))
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.r_max, nr=14, w_min=cfg.w_min, w_max=cfg.w_max, nw=14),
        dtype=dtype,
        device="cpu",
    )
    state = initial_branch_state(grid, cfg, dtype=dtype, device="cpu")
    boundaries = branch_boundary_conditions(cfg)
    residual_fn = lambda x: coupled_branch_residual(
        x,
        grid,
        cfg,
        eos_K=cfg.continuation_K_values[0],
        boundaries=boundaries,
    )
    matrix, metadata = assemble_coupled_colored_sparse_jacobian(
        PreconditionerBuildContext(
            residual_fn=residual_fn,
            x=state,
            rhs=torch.zeros_like(state).detach().cpu().numpy(),
            iteration=1,
            config=cfg.newton,
        ),
        grid,
    )
    expected_color_bound = 5 * preconditioner.color_separation**2 + 1
    assert metadata["active_color_count"] == expected_color_bound
    assert metadata["active_color_count"] < state.numel()

    relative_errors = []
    for seed in (778, 779, 780):
        generator = torch.Generator(device="cpu")
        generator.manual_seed(seed)
        direction = torch.randn(state.shape, dtype=dtype, generator=generator)
        direct = torch_jvp(residual_fn, state, direction).detach().cpu().numpy()
        assembled = matrix @ direction.detach().cpu().numpy()
        relative = math.sqrt(float(((assembled - direct) ** 2).sum())) / max(
            1.0,
            math.sqrt(float((direct**2).sum())),
        )
        relative_errors.append(relative)
    assert max(relative_errors) < 1.0e-10


def test_coupled_residual_changes_when_spatial_gauge_is_enabled() -> None:
    dtype = configure_backend(BackendConfig())
    cfg = BranchSmokeConfig()
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.r_max, nr=6, w_min=cfg.w_min, w_max=cfg.w_max, nw=4),
        dtype=dtype,
        device="cpu",
    )
    rr = grid.r_centers[:, None]
    ww = grid.w_centers[None, :]
    psi_real = torch.ones((grid.spec.nr, grid.spec.nw), dtype=dtype)
    psi_imag = 0.05 * rr * (ww - cfg.w_min)
    zeros = torch.zeros_like(psi_real)
    spatial_gauge = 0.1 * torch.ones_like(psi_real)
    base = pack_coupled_fields(
        CoupledFields(psi_real=psi_real, psi_imag=psi_imag, a0=zeros, ar=zeros, aw=zeros)
    )
    gauged = pack_coupled_fields(
        CoupledFields(
            psi_real=psi_real,
            psi_imag=psi_imag,
            a0=zeros,
            ar=spatial_gauge,
            aw=0.5 * spatial_gauge,
        )
    )
    boundaries = branch_boundary_conditions(cfg)
    base_residual = coupled_pde_residual(
        base,
        grid,
        cfg,
        eos_K=cfg.continuation_K_values[0],
        chemical_potential=cfg.initial_mu,
        boundaries=boundaries,
    )
    gauged_residual = coupled_pde_residual(
        gauged,
        grid,
        cfg,
        eos_K=cfg.continuation_K_values[0],
        chemical_potential=cfg.initial_mu,
        boundaries=boundaries,
    )
    assert torch.linalg.norm(gauged_residual - base_residual).item() > 1.0e-2


def test_boundary_sponge_profile_and_residual_wiring() -> None:
    dtype = configure_backend(BackendConfig())
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=2.0, nr=8, w_min=0.0, w_max=1.0, nw=4),
        dtype=dtype,
        device="cpu",
    )
    base_cfg = BranchSmokeConfig(
        r_max=2.0,
        w_min=0.0,
        w_max=1.0,
        solve_grid=(8, 4),
        sponge_enabled=False,
        sponge_width=0.0,
        sponge_matter_strength=0.0,
        sponge_gauge_strength=0.0,
    )
    off_profile = boundary_sponge_profile_torch(grid, base_cfg)
    assert torch.count_nonzero(off_profile).item() == 0

    sponge_cfg = replace(
        base_cfg,
        sponge_enabled=True,
        sponge_width=0.5,
        sponge_matter_strength=3.0,
        sponge_gauge_strength=5.0,
        sponge_power=2,
    )
    profile = boundary_sponge_profile_torch(grid, sponge_cfg)
    assert torch.min(profile).item() >= 0.0
    assert torch.max(profile).item() <= 1.0

    interior = (grid.r_centers[:, None] <= 1.5) & (grid.w_centers[None, :] <= 0.5)
    assert torch.count_nonzero(profile[interior]).item() == 0
    assert torch.all(profile[1:, -1] >= profile[:-1, -1])
    assert torch.all(profile[-1, 1:] >= profile[-1, :-1])

    def expected_profile(r_value: float, w_value: float) -> float:
        radial = max(0.0, min(1.0, (r_value - 1.5) / 0.5)) ** 2
        axial = max(0.0, min(1.0, (w_value - 0.5) / 0.5)) ** 2
        return 1.0 - (1.0 - radial) * (1.0 - axial)

    assert math.isclose(
        profile[-2, 2].item(),
        expected_profile(grid.r_centers[-2].item(), grid.w_centers[2].item()),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    )
    assert math.isclose(
        profile[-1, -1].item(),
        expected_profile(grid.r_centers[-1].item(), grid.w_centers[-1].item()),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    )

    rr = grid.r_centers[:, None]
    ww = grid.w_centers[None, :]
    fields = CoupledFields(
        psi_real=1.0 + 0.1 * rr + 0.2 * ww,
        psi_imag=0.3 + 0.05 * rr + torch.zeros_like(ww),
        a0=0.2 + 0.1 * ww + torch.zeros_like(rr),
        ar=0.4 + 0.03 * rr + torch.zeros_like(ww),
        aw=-0.2 + 0.07 * ww + torch.zeros_like(rr),
    )
    state = pack_coupled_fields(fields)
    boundaries = branch_boundary_conditions(base_cfg)
    residual_off = coupled_pde_residual(
        state,
        grid,
        base_cfg,
        eos_K=base_cfg.continuation_K_values[0],
        chemical_potential=base_cfg.initial_mu,
        boundaries=boundaries,
    )
    residual_on = coupled_pde_residual(
        state,
        grid,
        sponge_cfg,
        eos_K=sponge_cfg.continuation_K_values[0],
        chemical_potential=sponge_cfg.initial_mu,
        boundaries=boundaries,
    )
    difference = residual_on - residual_off
    expected = torch.stack(
        [
            sponge_cfg.sponge_matter_strength * profile * fields.psi_real,
            sponge_cfg.sponge_matter_strength * profile * fields.psi_imag,
            sponge_cfg.sponge_gauge_strength * profile * fields.a0,
            sponge_cfg.sponge_gauge_strength * profile * fields.ar,
            sponge_cfg.sponge_gauge_strength * profile * fields.aw,
        ],
        dim=0,
    )
    assert torch.allclose(difference, expected, rtol=0.0, atol=1.0e-12)
    assert torch.count_nonzero(difference[:, interior]).item() == 0


def test_conservation_sponge_accounting_is_interior_zero_with_expected_outer_sink() -> None:
    dtype = configure_backend(BackendConfig())
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=2.0, nr=8, w_min=0.0, w_max=1.0, nw=4),
        dtype=dtype,
        device="cpu",
    )
    cfg = BranchSmokeConfig(
        r_max=2.0,
        w_min=0.0,
        w_max=1.0,
        solve_grid=(8, 4),
        gauge_charge=0.35,
        sponge_enabled=True,
        sponge_width=0.5,
        sponge_matter_strength=3.0,
        sponge_gauge_strength=5.0,
        sponge_power=2,
    )
    rr = grid.r_centers[:, None]
    ww = grid.w_centers[None, :]
    fields = CoupledFields(
        psi_real=1.0 + 0.1 * rr + torch.zeros_like(ww),
        psi_imag=0.2 + 0.05 * ww + torch.zeros_like(rr),
        a0=0.3 + 0.02 * rr + torch.zeros_like(ww),
        ar=0.1 + 0.03 * ww + torch.zeros_like(rr),
        aw=-0.2 + 0.01 * rr + torch.zeros_like(ww),
    )
    profile = boundary_sponge_profile_torch(grid, cfg)
    sinks = _sponge_sink_densities(fields, grid, cfg)
    interior = (grid.r_centers[:, None] <= 1.5) & (grid.w_centers[None, :] <= 0.5)
    assert torch.sum(sinks["number"][interior] * grid.cell_volumes[interior]).item() == 0.0
    expected_number = torch.sum(
        cfg.sponge_matter_strength
        * profile
        * (fields.psi_real**2 + fields.psi_imag**2)
        * grid.cell_volumes
    )
    expected_energy = torch.sum(
        (
            cfg.sponge_matter_strength * profile * (fields.psi_real**2 + fields.psi_imag**2)
            + cfg.sponge_gauge_strength
            * profile
            * (fields.a0**2 + fields.ar**2 + fields.aw**2)
        )
        * grid.cell_volumes
    )
    assert expected_number.item() > 0.0
    assert expected_energy.item() > expected_number.item()
    assert torch.allclose(
        torch.sum(sinks["number"] * grid.cell_volumes),
        expected_number,
        rtol=0.0,
        atol=1.0e-13,
    )
    assert torch.allclose(
        torch.sum(sinks["charge"] * grid.cell_volumes),
        cfg.gauge_charge * expected_number,
        rtol=0.0,
        atol=1.0e-13,
    )
    assert torch.allclose(
        torch.sum(sinks["energy"] * grid.cell_volumes),
        expected_energy,
        rtol=0.0,
        atol=1.0e-13,
    )


def test_conservation_null_diagnostic_guard_uses_step4_vocabulary() -> None:
    assert null_floor_label(0.0) == "null diagnostic"
    assert null_floor_label(5.0e-15, floor=1.0e-14) == "null diagnostic"
    assert null_floor_label(2.0e-14, floor=1.0e-14) == "solver-floor diagnostic"


def test_boundary_interior_solution_difference_uses_fixed_volume_window() -> None:
    dtype = configure_backend(BackendConfig())
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=1.5, nr=6, w_min=0.0, w_max=1.0, nw=4),
        dtype=dtype,
        device="cpu",
    )
    ones = torch.ones((grid.spec.nr, grid.spec.nw), dtype=dtype)
    zeros = torch.zeros_like(ones)
    reference = pack_coupled_fields(
        CoupledFields(psi_real=ones, psi_imag=zeros, a0=zeros, ar=zeros, aw=zeros),
        torch.tensor(1.0, dtype=dtype),
    )
    shifted = pack_coupled_fields(
        CoupledFields(psi_real=1.1 * ones, psi_imag=zeros, a0=zeros, ar=zeros, aw=zeros),
        torch.tensor(1.0, dtype=dtype),
    )
    metric = interior_solution_difference(
        shifted,
        grid,
        reference,
        grid,
        InteriorWindow(r_min=0.0, r_max=1.0, w_min=0.25, w_max=0.75),
    )
    assert abs(metric["interior_relative_l2_change"] - 0.1) < 1.0e-13
    assert metric["interior_signal_l2_reference"] > 0.0


def test_boundary_sweep_configs_change_only_the_tested_variable() -> None:
    base = BranchSmokeConfig(r_max=2.0, w_max=1.6)
    study = BoundaryCharacterizationConfig(
        spacing=0.1,
        fixed_w_max=1.6,
        truncation_r_max_values=(1.8, 2.0, 2.2),
        impedance_alpha_scales=(0.5, 1.0, 2.0),
    )
    truncation = [asdict(branch) for branch in truncation_sweep_branches(base, study)]
    truncation_common = []
    for row in truncation:
        row = dict(row)
        row.pop("r_max")
        row.pop("solve_grid")
        truncation_common.append(row)
    assert truncation_common[0] == truncation_common[1] == truncation_common[2]

    impedance = [asdict(branch) for branch in impedance_sweep_branches(base, study)]
    impedance_common = []
    for row in impedance:
        row = dict(row)
        row.pop("matter_exit_impedance_alpha")
        row.pop("a0_radial_impedance_alpha")
        row.pop("a0_exit_impedance_alpha")
        impedance_common.append(row)
    assert impedance_common[0] == impedance_common[1] == impedance_common[2]
