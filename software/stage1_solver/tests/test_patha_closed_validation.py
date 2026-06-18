from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import torch

from stage1_solver.backend import configure_backend
from stage1_solver.config import BackendConfig, TensorGridSpec, WallGridSpec
from stage1_solver.coupled_branch import (
    ClosedCoupledFields,
    initial_closed_branch_state,
    pack_closed_coupled_fields,
    patha_closed_wall_terms,
    unpack_closed_coupled_fields,
)
from stage1_solver.grid import TensorProductGrid, WallGrid
from stage1_solver.patha_closed_newton import default_closed_s_sigma_spec, target_token_scan
from stage1_solver.patha_closed_validation import (
    balance_validation_diagnostic,
    conservative_restrict_closed_fields,
    default_validation_branch_config,
    independent_flux_divergence_recompute,
    independent_source_recompute,
    nonconservative_flux_divergence_discretization,
    validation_scan_paths,
)
from stage1_solver.patha_static_balance import resolve_s_sigma


def _grid(dtype: torch.dtype, nr: int = 4, nw: int = 4) -> TensorProductGrid:
    return TensorProductGrid.create(
        TensorGridSpec(r_max=1.4, nr=nr, w_min=0.0, w_max=1.0, nw=nw),
        dtype=dtype,
        device="cpu",
    )


def _branch():
    base = default_validation_branch_config()
    return replace(
        base,
        r_max=1.4,
        w_max=1.0,
        solve_grid=(4, 4),
        continuation_K_values=(0.08,),
    )


def test_independent_balance_recomputes_match_operator_terms() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _branch()
    grid = _grid(dtype)
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    fields, _ = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    spec = default_closed_s_sigma_spec(branch)
    provider = resolve_s_sigma(spec)
    terms = patha_closed_wall_terms(fields, grid, branch, s_sigma=provider)
    wall_grid = WallGrid.create(
        WallGridSpec(w_min=grid.spec.w_min, w_max=grid.spec.w_max, nw=grid.spec.nw),
        dtype=dtype,
        device="cpu",
    )

    source = independent_source_recompute(fields, grid, branch)
    flux = independent_flux_divergence_recompute(
        fields.r0,
        wall_grid,
        provider,
        lower_value=branch.r_mouth,
    )
    nonconservative_flux = nonconservative_flux_divergence_discretization(
        fields.r0,
        wall_grid,
        provider,
    )

    assert torch.allclose(source, terms.source)
    assert torch.allclose(flux, terms.flux_divergence)
    assert nonconservative_flux.shape == terms.flux_divergence.shape
    assert torch.isfinite(nonconservative_flux[1:-1]).all()

    diagnostic = balance_validation_diagnostic(
        state,
        grid,
        branch,
        spec,
        final_residual_linf=1.0,
        final_tolerance=1.0,
        nontrivial_floor=0.0,
        source_abs_tol=1.0e-12,
        flux_abs_tol=1.0e-12,
    )
    assert set(diagnostic["checks"]) == {"balance_terms_nontrivial"}
    assert (
        diagnostic["identity_checks"][
            "independent_source_recompute_matches_not_a_physics_gate"
        ]
        is True
    )
    assert (
        diagnostic["identity_checks"][
            "independent_flux_recompute_matches_not_a_physics_gate"
        ]
        is True
    )


def test_conservative_restrict_closed_fields_includes_r0_channel() -> None:
    dtype = configure_backend(BackendConfig())
    branch = _branch()
    coarse_grid = _grid(dtype, nr=4, nw=4)
    fine_grid = _grid(dtype, nr=8, nw=8)
    fine_state = initial_closed_branch_state(fine_grid, branch, dtype=dtype, device="cpu")
    fine_fields, mu = unpack_closed_coupled_fields(fine_state, fine_grid, has_chemical_potential=True)
    assert mu is not None
    r0 = torch.arange(fine_grid.spec.nw, dtype=dtype)
    fine_state = pack_closed_coupled_fields(
        ClosedCoupledFields(
            psi_real=fine_fields.psi_real,
            psi_imag=fine_fields.psi_imag,
            a0=fine_fields.a0,
            ar=fine_fields.ar,
            aw=fine_fields.aw,
            r0=r0,
        ),
        mu,
    )

    restricted = conservative_restrict_closed_fields(fine_state, fine_grid, coarse_grid)

    assert torch.allclose(
        restricted["r0"],
        torch.tensor([0.5, 2.5, 4.5, 6.5], dtype=dtype),
    )
    assert set(restricted) == {"psi_real", "psi_imag", "a0", "ar", "aw", "r0"}


def test_validation_scan_paths_are_target_token_clean(tmp_path: Path) -> None:
    report = tmp_path / "patha_chunk1c_self_consistent_validation.md"
    report.write_text("target-blind placeholder validation\n", encoding="utf-8")
    scan = target_token_scan(validation_scan_paths(report))
    assert scan["passed"] is True
