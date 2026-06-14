from __future__ import annotations

import ast
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import re

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
    P2CentrifugalMMSConfig,
    P2DrivenConfig,
    P2DrivenMMSConfig,
    P2MaxwellAngularMMSConfig,
    P2TangentConfig,
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
    confinement_potential_torch,
    coupled_branch_residual,
    coupled_pde_residual,
    geometry_radius_torch,
    initial_branch_state,
    _matter_number_current,
    pack_coupled_fields,
    run_branch_continuation,
)
from stage1_solver.error_budget import (
    CombinationSensitivity,
    ErrorFloorScales,
    _observable_set_matches_step4,
    combine_uncertainty,
    compose_observable_budget,
    recorded_prior_results,
    run_step7,
)
from stage1_solver.grid import RadialGrid, TensorProductGrid
from stage1_solver.manifest import write_manifest
from stage1_solver.newton import PreconditionerBuildContext, finite_difference_jvp_check
from stage1_solver.mms_benchmarks import run_p2_centrifugal_mms, run_p2_maxwell_angular_mms
from stage1_solver.operators import (
    integrate,
    radial_current,
    radial_fluxes,
    radial_laplacian,
    tensor_flux_divergence,
    tensor_laplacian,
    tensor_vector_face_fluxes,
)
from stage1_solver.p2_tangent import (
    _diagnostics_digest,
    P2TangentFields,
    apply_p2_tangent,
    p2_tangent_fd_check,
    p2_wellposedness_check,
    wall_to_matter_coefficient_torch,
    with_step8a_preconditioners,
    zero_p2_tangent_state,
)
from stage1_solver.p2_driven_absorber import (
    _driven_digest,
    _run_tensor_matter_frequency_cap_mms,
    apply_p2_driven_tangent,
    driven_observable_convergence_gate,
    p2_matter_free_particle_dispersion_check,
    p2_cap_profile_torch,
    p2_driven_frequency_terms,
    p2_omega_zero_static_limit_check,
    run_p2_driven_mms,
    run_reflection_study,
)
from stage1_solver.p2_conservation_response import (
    Step8CStudyConfig,
    _current_rows,
    _sampling_stability_passes,
    _status as _step8c_status,
    _step8c_digest,
    continuity_source_projection,
    linearized_phasor_number_current,
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


def test_error_budget_floor_enforcement_carries_solver_floor() -> None:
    floors = ErrorFloorScales(
        solver_abs=1.0e-8,
        discretization_relative_floor=1.0e-4,
        boundary_relative_floor=0.0,
        conservation_relative_floor=0.0,
    )
    row = compose_observable_budget(
        {
            "observable": "density_mass",
            "label": "density mass integral",
            "finest_grid": "synthetic",
            "finest_dof": 1,
            "finest_value": 1.0,
            "finest_error_estimate": 1.0e-12,
            "verdict": "solver-floor diagnostic",
        },
        floors,
    )
    assert row["u_disc"] == 1.0e-12
    assert row["u_total"] == floors.solver_abs
    assert row["floor_limited"] is True
    assert row["verdict"] == "solver-floor diagnostic"


def test_error_budget_monotonic_rss_composition_bounds_components() -> None:
    base = combine_uncertainty(
        solver_floor=1.0e-8,
        discretization_abs=3.0e-6,
        boundary_abs=0.0,
        conservation_abs=0.0,
    )
    with_boundary = combine_uncertainty(
        solver_floor=1.0e-8,
        discretization_abs=3.0e-6,
        boundary_abs=4.0e-6,
        conservation_abs=0.0,
    )
    assert with_boundary["rss_total"] > base["rss_total"]
    assert with_boundary["rss_total"] >= 4.0e-6
    assert with_boundary["max_alternative"] <= with_boundary["rss_total"]
    assert with_boundary["rss_total"] <= with_boundary["sum_bound"]


def test_error_budget_relative_conservation_floor_scope_mapping() -> None:
    floors = ErrorFloorScales(
        solver_abs=1.0e-8,
        discretization_relative_floor=1.0e-4,
        boundary_relative_floor=0.01,
        conservation_relative_floor=0.02,
    )
    coupled = compose_observable_budget(
        {
            "observable": "scalar_gauge_l2",
            "label": "A0 L2 norm",
            "finest_grid": "synthetic",
            "finest_dof": 1,
            "finest_value": 2.0,
            "finest_error_estimate": 0.0,
            "verdict": "expected-order convergence",
        },
        floors,
    )
    uncoupled = compose_observable_budget(
        {
            "observable": "chemical_potential",
            "label": "chemical potential",
            "finest_grid": "synthetic",
            "finest_dof": 1,
            "finest_value": 2.0,
            "finest_error_estimate": 0.0,
            "verdict": "expected-order convergence",
        },
        floors,
    )
    assert coupled["u_conservation"] == 0.04
    assert uncoupled["u_conservation"] == 0.0
    assert coupled["u_boundary"] == uncoupled["u_boundary"] == 0.02


def test_error_budget_null_sector_guard_has_no_manufactured_uncertainty() -> None:
    floors = ErrorFloorScales(
        solver_abs=1.0e-8,
        discretization_relative_floor=1.0e-4,
        boundary_relative_floor=0.01,
        conservation_relative_floor=0.02,
    )
    row = compose_observable_budget(
        {
            "observable": "spatial_gauge_l2",
            "label": "spatial gauge L2 norm",
            "finest_grid": "synthetic",
            "finest_dof": 1,
            "finest_value": 0.0,
            "finest_error_estimate": 0.0,
            "verdict": "null diagnostic",
        },
        floors,
    )
    assert row["u_total"] == 0.0
    assert row["u_solver"] == 0.0
    assert row["relative_uncertainty"] is None
    assert row["dominant_component"] == "null"
    assert row["verdict"] == "null diagnostic"


def test_error_budget_asserted_items_are_reported_but_not_counted(tmp_path) -> None:
    recorded = recorded_prior_results()
    recorded["step4"]["passed"] = False
    recorded["step5"]["passed"] = False
    recorded["step6"]["passed"] = False
    cfg = HarnessConfig(
        run_root=str(tmp_path / "run"),
        report_path=str(tmp_path / "report.md"),
    )

    results = run_step7(
        cfg,
        step4_results=recorded["step4"],
        step5_results=recorded["step5"],
        step6_results=recorded["step6"],
    )

    assert set(results["pass_checks"]) == {
        "non_null_uncertainties_floor_at_solver",
        "solver_floor_limited_rows_flagged",
        "null_sectors_remain_null",
        "conservation_floor_scoped",
        "boundary_floor_scoped",
        "observable_set_matches_step4",
    }
    assert all(results["pass_checks"].values())
    assert results["passed"] is True
    assert all(key.endswith("_not_a_physics_gate") for key in results["asserted_checks"])
    assert results["asserted_checks"]["prior_step4_passed_not_a_physics_gate"] is False
    assert results["asserted_checks"]["prior_step5_passed_not_a_physics_gate"] is False
    assert results["asserted_checks"]["prior_step6_passed_not_a_physics_gate"] is False


def test_error_budget_observable_set_check_is_computed() -> None:
    recorded = recorded_prior_results()
    assert _observable_set_matches_step4(recorded["step4"]) is True

    dropped = recorded_prior_results()["step4"]
    dropped["observable_summary"] = [
        row for row in dropped["observable_summary"] if row["observable"] != "density_mass"
    ]

    assert _observable_set_matches_step4(dropped) is False


def test_error_budget_digest_is_deterministic(tmp_path) -> None:
    cfg_a = HarnessConfig(
        run_root=str(tmp_path / "run_a"),
        report_path=str(tmp_path / "report_a.md"),
    )
    cfg_b = HarnessConfig(
        run_root=str(tmp_path / "run_b"),
        report_path=str(tmp_path / "report_b.md"),
    )
    first = run_step7(cfg_a, sensitivity=CombinationSensitivity(material_ratio_threshold=1.05))
    second = run_step7(cfg_b, sensitivity=CombinationSensitivity(material_ratio_threshold=1.05))
    assert first["diagnostics_digest"] == second["diagnostics_digest"]
    assert first["component_floors"] == second["component_floors"]
    assert first["observable_budget"] == second["observable_budget"]


def _tiny_p2_harness_config(tmp_path) -> HarnessConfig:
    base = HarnessConfig(
        run_root=str(tmp_path / "run"),
        report_path=str(tmp_path / "report.md"),
    )
    branch = BranchSmokeConfig(
        solve_grid=(4, 4),
        ladder_levels=((4, 4),),
        continuation_K_values=(0.05,),
        newton=replace(
            BranchSmokeConfig().newton,
            residual_atol=1.0e-8,
            residual_rtol=1.0e-8,
            gmres_rtol=1.0e-8,
            gmres_atol=1.0e-10,
            gmres_restart=96,
            gmres_maxiter=8,
            max_newton_iters=10,
        ),
    )
    p2 = replace(
        P2TangentConfig(),
        convergence_levels=((4, 4), (8, 8), (16, 16)),
        fd_check_grid=(4, 4),
        wellposed_grid=(4, 4),
        min_observable_order=0.0,
        newton=replace(
            P2TangentConfig().newton,
            residual_atol=1.0e-9,
            residual_rtol=1.0e-9,
            gmres_rtol=1.0e-9,
            gmres_atol=1.0e-11,
            gmres_restart=128,
            gmres_maxiter=8,
        ),
    )
    mms = replace(
        base.mms,
        p2_centrifugal=replace(
            P2CentrifugalMMSConfig(),
            grid_levels=((24, 20), (48, 40), (96, 80)),
            final_error_max=3.0e-3,
        ),
        p2_maxwell_angular=replace(
            P2MaxwellAngularMMSConfig(),
            grid_levels=((24, 20), (48, 40), (96, 80)),
            final_error_max=4.0e-3,
        ),
    )
    return replace(base, branch=branch, p2_tangent=p2, mms=mms)


def _tiny_step8b_harness_config(tmp_path) -> HarnessConfig:
    base = _tiny_p2_harness_config(tmp_path)
    p2_driven = replace(
        P2DrivenConfig(),
        drive_frequencies=(0.05, 1.5, 6.0),
        primary_omega=6.0,
        propagating_omega=6.0,
        convergence_levels=((4, 4), (8, 8), (16, 16)),
        reflection_grid=(8, 8),
        response_table_grid=(6, 6),
        reflection_fit_window=(0.25, 0.75),
        cap_width=0.8,
        cap_strength=200.0,
        min_observable_order=1.45,
    )
    mms = replace(
        base.mms,
        p2_driven=replace(
            P2DrivenMMSConfig(),
            tensor_grid_levels=((24, 20), (48, 40), (96, 80)),
            wall_grid_levels=(32, 64, 128),
            tensor_final_error_max=1.0e-2,
            wall_final_error_max=2.0e-3,
        ),
    )
    return replace(base, p2_driven=p2_driven, mms=mms)


def _tiny_p2_background() -> tuple[HarnessConfig, TensorProductGrid, torch.Tensor]:
    dtype = configure_backend(BackendConfig())
    branch = BranchSmokeConfig(
        solve_grid=(4, 4),
        ladder_levels=((4, 4),),
        continuation_K_values=(0.05,),
        newton=replace(
            BranchSmokeConfig().newton,
            residual_atol=1.0e-8,
            residual_rtol=1.0e-8,
            gmres_rtol=1.0e-8,
            gmres_atol=1.0e-10,
            gmres_restart=96,
            gmres_maxiter=8,
            max_newton_iters=10,
        ),
    )
    cfg = with_step8a_preconditioners(HarnessConfig(branch=branch))
    grid = TensorProductGrid.create(
        TensorGridSpec(
            r_max=cfg.branch.r_max,
            nr=4,
            w_min=cfg.branch.w_min,
            w_max=cfg.branch.w_max,
            nw=4,
        ),
        dtype=dtype,
        device="cpu",
    )
    state, summary = run_branch_continuation(cfg, grid, grid_name="test_p2_background")
    assert summary["converged"] is True
    return cfg, grid, state


def test_p2_tangent_matches_fd_of_mode_residual() -> None:
    cfg, grid, background = _tiny_p2_background()
    boundaries = branch_boundary_conditions(cfg.branch)
    check = p2_tangent_fd_check(
        background,
        grid,
        cfg.branch,
        cfg.p2_tangent,
        eos_K=cfg.branch.continuation_K_values[-1],
        boundaries=boundaries,
        epsilon=cfg.p2_tangent.newton.finite_difference_jvp_epsilon,
        seed=1234,
    )
    assert (
        check["relative_residual"] <= cfg.p2_tangent.tangent_fd_relative_tol
        or check["absolute_residual"] <= cfg.p2_tangent.tangent_fd_absolute_tol
    )


def test_wall_to_matter_coefficient_matches_confinement_radius_fd() -> None:
    dtype = configure_backend(BackendConfig())
    cfg = BranchSmokeConfig(
        r_max=1.7,
        w_min=0.0,
        w_max=1.3,
        radial_wall_strength=0.73,
        axial_trap_strength=0.19,
        r_mouth=1.25,
        r_exit=0.82,
        geometry_decay_length=0.67,
    )
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.r_max, nr=6, w_min=cfg.w_min, w_max=cfg.w_max, nw=5),
        dtype=dtype,
        device="cpu",
    )
    eps = 1.0e-6
    plus = replace(cfg, r_mouth=cfg.r_mouth + eps, r_exit=cfg.r_exit + eps)
    minus = replace(cfg, r_mouth=cfg.r_mouth - eps, r_exit=cfg.r_exit - eps)
    fd_dv_d_radius = (
        confinement_potential_torch(grid, plus) - confinement_potential_torch(grid, minus)
    ) / (2.0 * eps)

    coeff = wall_to_matter_coefficient_torch(grid, cfg)
    radius = geometry_radius_torch(grid.w_centers, cfg)[None, :]
    expected_shape = 4.0 * cfg.radial_wall_strength * grid.r_centers[:, None] ** 4 / radius**5
    torch.testing.assert_close(coeff, expected_shape, rtol=0.0, atol=1.0e-14)
    torch.testing.assert_close(fd_dv_d_radius, -coeff, rtol=0.0, atol=1.0e-6)


def test_p2_new_mms_operator_pieces_are_second_order(tmp_path) -> None:
    cfg = _tiny_p2_harness_config(tmp_path)
    centrifugal = run_p2_centrifugal_mms(cfg)
    maxwell = run_p2_maxwell_angular_mms(cfg)
    assert centrifugal["passed"] is True
    assert maxwell["passed"] is True
    assert centrifugal["rows"][-1]["observed_order"] > 1.85
    assert maxwell["rows"][-1]["observed_order"] > 1.75


def test_p2_driven_frequency_and_cap_mms_are_second_order(tmp_path) -> None:
    cfg = _tiny_step8b_harness_config(tmp_path)
    result = run_p2_driven_mms(cfg)
    assert result["passed"] is True
    assert set(result["sections"]) == {
        "matter_frequency_cap",
        "maxwell_frequency_cap",
        "wall_frequency_cap",
    }
    for section in result["sections"].values():
        assert section["rows"][-1]["observed_order"] > 1.75


def test_p2_driven_matter_mms_rejects_old_negated_time_sign(tmp_path) -> None:
    cfg = _tiny_step8b_harness_config(tmp_path)
    wrong = _run_tensor_matter_frequency_cap_mms(cfg, matter_time_sign=-1)
    correct = _run_tensor_matter_frequency_cap_mms(cfg, matter_time_sign=1)
    assert correct["passed"] is True
    assert wrong["passed"] is False
    assert wrong["rows"][-1]["error"] > 0.5


def test_p2_static_operator_is_wellposed_on_small_grid() -> None:
    cfg, grid, background = _tiny_p2_background()
    check = p2_wellposedness_check(background, grid, cfg)
    assert check["state_size"] == zero_p2_tangent_state(grid).numel()
    assert check["smallest_singular_value"] > cfg.p2_tangent.smallest_singular_min


def test_p2_driven_omega_zero_recovers_static_operator() -> None:
    cfg, grid, background = _tiny_p2_background()
    check = p2_omega_zero_static_limit_check(background, grid, cfg, seed=2027)
    assert check["relative_residual"] <= cfg.p2_tangent.tangent_fd_relative_tol
    assert check["absolute_residual"] <= cfg.p2_tangent.tangent_fd_absolute_tol


def test_p2_driven_matter_bdg_sign_structure_is_offdiagonal() -> None:
    dtype = configure_backend(BackendConfig())
    base = HarnessConfig()
    cfg = replace(base, p2_driven=replace(base.p2_driven, cap_enabled=False))
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.branch.r_max, nr=4, w_min=cfg.branch.w_min, w_max=cfg.branch.w_max, nw=4),
        dtype=dtype,
        device="cpu",
    )
    n = grid.spec.nr * grid.spec.nw
    state = torch.zeros(zero_p2_tangent_state(grid).shape, dtype=torch.complex128)
    state[:n] = 3.0 - 0.25j
    state[n : 2 * n] = 2.0 + 0.5j
    omega = 1.7
    terms = p2_driven_frequency_terms(state, grid, cfg, omega=omega, cap_enabled=False)
    assert torch.max(torch.abs(terms[:n] - 1j * cfg.branch.hbar * omega * state[n : 2 * n])).item() < 1.0e-14
    assert torch.max(torch.abs(terms[n : 2 * n] + 1j * cfg.branch.hbar * omega * state[:n])).item() < 1.0e-14


def test_p2_driven_matter_dispersion_check_rejects_old_sign() -> None:
    check = p2_matter_free_particle_dispersion_check(HarnessConfig())
    assert check["passed"] is True
    assert check["relative_residual"] < 1.0e-14
    assert check["old_negated_sign_relative_residual"] > 1.0


def test_p2_cap_profile_is_exit_only_and_complex_absorbing() -> None:
    dtype = configure_backend(BackendConfig())
    cfg = HarnessConfig()
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=2.0, nr=6, w_min=0.0, w_max=1.6, nw=8),
        dtype=dtype,
        device="cpu",
    )
    cell_sigma, wall_sigma = p2_cap_profile_torch(grid, cfg)
    assert torch.count_nonzero(cell_sigma[:, grid.w_centers < grid.spec.w_max - cfg.p2_driven.cap_width]).item() == 0
    assert torch.max(wall_sigma).item() > 0.0
    assert torch.all(wall_sigma[1:] >= wall_sigma[:-1])
    assert cfg.p2_driven.cap_profile == "smooth_polynomial_rational_cap"


def test_step8c_linearized_phasor_current_matches_real_branch_fd() -> None:
    dtype = configure_backend(BackendConfig())
    cfg = BranchSmokeConfig(r_max=1.4, w_min=0.0, w_max=1.2, gauge_charge=0.37)
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=cfg.r_max, nr=8, w_min=cfg.w_min, w_max=cfg.w_max, nw=6),
        dtype=dtype,
        device="cpu",
    )
    rr = grid.r_centers[:, None]
    ww = grid.w_centers[None, :]
    background = CoupledFields(
        psi_real=1.0 + 0.07 * rr + 0.03 * ww,
        psi_imag=0.11 + 0.02 * rr - 0.04 * ww,
        a0=torch.zeros((grid.spec.nr, grid.spec.nw), dtype=dtype),
        ar=0.04 * rr + 0.01 * ww,
        aw=-0.03 * rr + 0.02 * ww,
    )
    perturb = P2TangentFields(
        psi_real=0.05 * rr**2 - 0.02 * ww,
        psi_imag=-0.03 * rr + 0.04 * ww**2,
        a0=torch.zeros((grid.spec.nr, grid.spec.nw), dtype=dtype),
        ar=0.02 * rr * (1.0 + ww),
        aw=-0.015 * (1.0 + rr) * ww,
        eta=torch.zeros(grid.spec.nw, dtype=dtype),
    )
    eps = 1.0e-6

    def shifted(sign: float) -> CoupledFields:
        return CoupledFields(
            psi_real=background.psi_real + sign * eps * perturb.psi_real,
            psi_imag=background.psi_imag + sign * eps * perturb.psi_imag,
            a0=background.a0 + sign * eps * perturb.a0,
            ar=background.ar + sign * eps * perturb.ar,
            aw=background.aw + sign * eps * perturb.aw,
        )

    plus_r, plus_w = _matter_number_current(shifted(1.0), grid, cfg)
    minus_r, minus_w = _matter_number_current(shifted(-1.0), grid, cfg)
    fd_r = (plus_r - minus_r) / (2.0 * eps)
    fd_w = (plus_w - minus_w) / (2.0 * eps)
    lin_r, lin_w = linearized_phasor_number_current(background, perturb, grid, cfg)
    torch.testing.assert_close(torch.real(lin_r), fd_r, rtol=2.0e-8, atol=2.0e-8)
    torch.testing.assert_close(torch.real(lin_w), fd_w, rtol=2.0e-8, atol=2.0e-8)
    assert torch.max(torch.abs(torch.imag(lin_r))).item() == 0.0
    assert torch.max(torch.abs(torch.imag(lin_w))).item() == 0.0


def test_step8c_continuity_projection_uses_background_phase() -> None:
    dtype = configure_backend(BackendConfig())
    cfg = BranchSmokeConfig(hbar=2.0)
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=1.0, nr=4, w_min=0.0, w_max=1.0, nw=2),
        dtype=dtype,
        device="cpu",
    )
    ones = torch.ones((grid.spec.nr, grid.spec.nw), dtype=dtype)
    background = CoupledFields(
        psi_real=2.0 * ones,
        psi_imag=0.5 * ones,
        a0=torch.zeros_like(ones),
        ar=torch.zeros_like(ones),
        aw=torch.zeros_like(ones),
    )
    residual = P2TangentFields(
        psi_real=(1.0 + 2.0j) * ones.to(torch.complex128),
        psi_imag=(3.0 - 4.0j) * ones.to(torch.complex128),
        a0=torch.zeros_like(ones, dtype=torch.complex128),
        ar=torch.zeros_like(ones, dtype=torch.complex128),
        aw=torch.zeros_like(ones, dtype=torch.complex128),
        eta=torch.zeros(grid.spec.nw, dtype=torch.complex128),
    )
    projected = continuity_source_projection(background, residual, cfg)
    expected = (2.0 / cfg.hbar) * (
        background.psi_real.to(torch.complex128) * residual.psi_imag
        - background.psi_imag.to(torch.complex128) * residual.psi_real
    )
    torch.testing.assert_close(projected, expected, rtol=0.0, atol=0.0)


def test_step8c_current_ratio_reports_static_null_sentinel() -> None:
    dtype = configure_backend(BackendConfig())
    cfg = HarnessConfig()
    grid = TensorProductGrid.create(
        TensorGridSpec(r_max=1.0, nr=4, w_min=0.0, w_max=1.0, nw=3),
        dtype=dtype,
        device="cpu",
    )
    ones = torch.ones((grid.spec.nr, grid.spec.nw), dtype=dtype)
    background = CoupledFields(
        psi_real=ones,
        psi_imag=torch.zeros_like(ones),
        a0=torch.zeros_like(ones),
        ar=torch.zeros_like(ones),
        aw=torch.zeros_like(ones),
    )
    perturbation = P2TangentFields(
        psi_real=torch.zeros_like(ones, dtype=torch.complex128),
        psi_imag=torch.zeros_like(ones, dtype=torch.complex128),
        a0=torch.zeros_like(ones, dtype=torch.complex128),
        ar=torch.zeros_like(ones, dtype=torch.complex128),
        aw=torch.zeros_like(ones, dtype=torch.complex128),
        eta=torch.zeros(grid.spec.nw, dtype=torch.complex128),
    )
    jr = torch.ones_like(ones, dtype=torch.complex128)
    jw = torch.zeros_like(ones, dtype=torch.complex128)

    rows = _current_rows(
        level_index=0,
        grid=grid,
        config=cfg,
        background=background,
        perturbation=perturbation,
        chemical_potential=torch.as_tensor(2.0, dtype=dtype),
        jr=jr,
        jw=jw,
    )

    assert rows[0]["phasor_current_l2"] > 0.0
    assert rows[0]["static_branch_current_l2"] == 0.0
    assert {row["non_null_vs_static_ratio"] for row in rows} == {"static_exactly_null"}


def test_step8c_sampling_stability_uses_primary_relative_rule_only() -> None:
    study = Step8CStudyConfig(coefficient_sampling_relative_tol=0.35)
    unstable_primary = [
        {
            "level": 2,
            "functional": "f",
            "coefficient": "taylor_c0",
            "counts_for_stability": True,
            "relative_change": 0.36,
            "absolute_change": 0.0,
        }
    ]
    assert _sampling_stability_passes(unstable_primary, study) is False

    stable_primary_with_failed_stress = [
        {
            "level": 2,
            "functional": "f",
            "coefficient": "taylor_c0",
            "counts_for_stability": True,
            "relative_change": 0.01,
            "absolute_change": 1.0,
        },
        {
            "level": 2,
            "functional": "f",
            "coefficient": "taylor_c1",
            "counts_for_stability": False,
            "relative_change": 1.29,
            "absolute_change": 0.0,
        },
    ]
    assert _sampling_stability_passes(stable_primary_with_failed_stress, study) is True


def test_step8c_source_balance_checks_are_identity_not_physics_gates() -> None:
    root = Path(__file__).resolve().parents[1]
    source = (root / "src/stage1_solver/p2_conservation_response.py").read_text(encoding="utf-8")
    pass_block = source.split("pass_checks = {", 1)[1].split("asserted_checks = {", 1)[0]
    identity_block = source.split("identity_checks = {", 1)[1].split("source_balance_note =", 1)[0]

    assert "source_balance_closes_on_finest_grid_not_a_physics_gate" in identity_block
    assert "source_balance_residual_decreases_not_a_physics_gate" in identity_block
    assert "source_balance_closes_on_finest_grid" not in pass_block
    assert "source_balance_residual_decreases" not in pass_block
    assert "number_divergence cancels algebraically" in source


def test_step8c_status_reports_resolution_limited_pass() -> None:
    assert _step8c_status({"passed": True, "resolution_limited": True}) == "PASS_WITH_RESOLUTION_LIMIT"
    assert _step8c_status({"passed": True, "resolution_limited": False}) == "PASS"
    assert _step8c_status({"passed": False, "resolution_limited": True}) == "FAIL"


def test_p2_driven_convergence_gate_fails_synthetic_drifting_observable() -> None:
    cfg = replace(
        HarnessConfig(),
        p2_driven=replace(
            HarnessConfig().p2_driven,
            min_observable_order=1.45,
            scoped_min_observable_order=1.45,
            scoped_convergence_observables=("total_response_l2",),
        ),
    )
    summary = [
        {
            "observable": "total_response_l2",
            "verdict": "drifts",
            "last_observed_order": 0.1,
        },
        {
            "observable": "driven_residual_linf",
            "verdict": "solver-floor diagnostic",
            "last_observed_order": None,
        },
    ]
    gate = driven_observable_convergence_gate(summary, cfg)
    assert gate["passed"] is False
    assert gate["gated_observable_names"] == ["total_response_l2"]
    assert gate["failing_observable_names"] == ["total_response_l2"]


def test_p2_solve_path_is_target_blind() -> None:
    root = Path(__file__).resolve().parents[1]
    solve_path_files = [
        root / "src/stage1_solver/p2_tangent.py",
        root / "src/stage1_solver/p2_tangent_harness.py",
        root / "src/stage1_solver/p2_driven_absorber.py",
        root / "src/stage1_solver/step8b_harness.py",
        root / "src/stage1_solver/p2_conservation_response.py",
        root / "src/stage1_solver/step8c_harness.py",
    ]
    solve_path_sources = {
        path: path.read_text(encoding="utf-8")
        for path in solve_path_files
    }
    solve_path_text = "\n".join(
        solve_path_sources.values()
    )
    forbidden = [
        "R_norm",
        "R_pole",
        "chi_Q",
        "χ_Q",
        "N_Q",
        "10.8",
        "54.0",
        "54/5",
        "54 G",
        "physical_nonlinear_model",
    ]
    for token in forbidden:
        assert token not in solve_path_text
    for path, source in solve_path_sources.items():
        tree = ast.parse(source, filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                module_parts = {
                    part
                    for alias in node.names
                    for part in alias.name.split(".")
                }
            elif isinstance(node, ast.ImportFrom):
                module_parts = set((node.module or "").split(".")) if node.module else set()
                module_parts.update(
                    part
                    for alias in node.names
                    for part in alias.name.split(".")
                )
            else:
                continue
            assert module_parts.isdisjoint({"benchmarks", "references", "targets"}), path
        for pattern in (
            r"(?<![\w.])10\.8(?![\w.])",
            r"(?<![\w.])54\.0(?![\w.])",
            r"(?<![\w.])54\s*/\s*5(?![\w.])",
        ):
            assert re.search(pattern, source) is None, (path, pattern)


def test_p2_digest_is_deterministic(tmp_path) -> None:
    results = {
        "mms": {
            "p2_centrifugal": {
                "rows": [{"grid": "a", "error": 1.0e-3, "manifest": str(tmp_path / "a")}],
                "passed": True,
            },
            "p2_maxwell_angular": {
                "rows": [{"grid": "a", "error": 2.0e-3, "manifest": str(tmp_path / "b")}],
                "passed": True,
            },
        },
        "tangent_fd_check": {"relative_residual": 1.0e-11, "absolute_residual": 1.0e-12},
        "wellposedness": {"smallest_singular_value": 1.0, "condition_number": 3.0},
        "level_rows": [
            {
                "grid": "tiny",
                "solve": {"manifest": str(tmp_path / "solve_a"), "final_residual_norm": 1.0e-12},
                "background": {"manifest": str(tmp_path / "background_a")},
            }
        ],
        "observable_summary": [{"observable": "total_response_l2", "last_observed_order": 2.0}],
        "pass_checks": {"computed_gate": True},
        "asserted_checks": {"target_blind_not_a_physics_gate": True},
    }
    moved = json.loads(json.dumps(results))
    moved["mms"]["p2_centrifugal"]["rows"][0]["manifest"] = str(tmp_path / "elsewhere")
    moved["level_rows"][0]["solve"]["manifest"] = str(tmp_path / "solve_b")
    first = _diagnostics_digest(results)
    second = _diagnostics_digest(moved)
    assert first == second
    assert len(first) == 16


def test_p2_driven_digest_is_deterministic(tmp_path) -> None:
    results = {
        "mms": {"sections": {"matter_frequency_cap": {"passed": True, "rows": []}}},
        "static_tangent_fd_check": {"relative_residual": 1.0e-12},
        "driven_fd_check": {"relative_residual": 2.0e-12},
        "omega_zero_static_limit": {"relative_residual": 0.0},
        "matter_dispersion": {"relative_residual": 0.0},
        "wellposedness": {"smallest_singular_value": 1.0},
        "level_rows": [{"grid": "a", "solve": {"manifest": str(tmp_path / "a")}}],
        "observable_summary": [{"observable": "total_response_l2", "last_observed_order": 2.0}],
        "convergence_gate": {"passed": True, "gated_observable_names": ["total_response_l2"]},
        "resolution_diagnostics": [{"grid": "a", "points_per_wavelength": 12.0}],
        "reflection": {"absorbed": {"reflection_coefficient": 1.0e-3}},
        "response_vs_omega": [{"omega": 1.0, "total_response_l2": 2.0}],
        "pass_checks": {"computed_gate": True},
        "asserted_checks": {"target_blind_not_a_physics_gate": True},
    }
    moved = json.loads(json.dumps(results))
    moved["level_rows"][0]["solve"]["manifest"] = str(tmp_path / "b")
    first = _driven_digest(results)
    second = _driven_digest(moved)
    assert first == second
    assert len(first) == 16


def test_step8c_digest_is_deterministic(tmp_path) -> None:
    results = {
        "study": Step8CStudyConfig().to_dict(),
        "solve_rows": [
            {
                "grid": "a",
                "omega": 0.1,
                "driven_residual_linf": 1.0e-12,
                "manifest": str(tmp_path / "solve_a"),
                "background_manifest": str(tmp_path / "background_a"),
            }
        ],
        "current_rows": [{"sector": "number", "phasor_current_l2": 1.0}],
        "gauss_closure_rows": [{"surface": "nested", "relative_residual": 1.0e-3}],
        "source_balance_rows": [
            {"sector": "number", "interior_balance_l2_relative": 1.0e-4}
        ],
        "response_rows": [{"omega": 0.1, "scalar_gauge_cap_free_l2": 2.0}],
        "coefficient_summary": [{"functional": "scalar", "coefficient": "taylor_c0"}],
        "coefficient_budget": [{"functional": "scalar", "u_total": 1.0e-3}],
        "omega_stability_rows": [{"functional": "scalar", "relative_change": 1.0e-2}],
        "pass_checks": {"computed_gate": True},
        "identity_checks": {"identity_not_a_physics_gate": True},
        "identity_check_notes": {"identity_not_a_physics_gate": "not counted"},
        "asserted_checks": {"target_blind_not_a_physics_gate": True},
        "resolution_limited": True,
    }
    moved = json.loads(json.dumps(results))
    moved["solve_rows"][0]["manifest"] = str(tmp_path / "solve_b")
    moved["solve_rows"][0]["background_manifest"] = str(tmp_path / "background_b")
    first = _step8c_digest(results)
    second = _step8c_digest(moved)
    assert first == second
    assert len(first) == 16


def test_p2_reflection_metric_has_absorber_control_teeth(tmp_path) -> None:
    cfg = _tiny_step8b_harness_config(tmp_path)
    reflection = run_reflection_study(cfg)
    assert reflection["pass_checks"]["propagating_wall_k_real"] is True
    assert reflection["absorbed"]["reflection_coefficient"] < reflection["target_blind_floor"]
    assert reflection["reflecting_control"]["reflection_coefficient"] > reflection["target_blind_floor"]
    assert reflection["control_contrast"] >= cfg.p2_driven.reflection_control_min_contrast


def test_p2_centrifugal_terms_change_residual_vs_l0() -> None:
    cfg, grid, background = _tiny_p2_background()
    boundaries = branch_boundary_conditions(cfg.branch)
    generator = torch.Generator(device="cpu")
    generator.manual_seed(2026)
    direction = torch.randn(
        zero_p2_tangent_state(grid).shape,
        dtype=grid.dtype,
        generator=generator,
    )
    l2 = apply_p2_tangent(
        direction,
        background,
        grid,
        cfg.branch,
        cfg.p2_tangent,
        eos_K=cfg.branch.continuation_K_values[-1],
        boundaries=boundaries,
    )
    l0_cfg = replace(cfg.p2_tangent, spherical_l=0)
    l0 = apply_p2_tangent(
        direction,
        background,
        grid,
        cfg.branch,
        l0_cfg,
        eos_K=cfg.branch.continuation_K_values[-1],
        boundaries=boundaries,
    )
    assert torch.linalg.norm(l2 - l0).item() > 1.0e-3
