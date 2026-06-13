"""Step 6 stationary conservation diagnostics for the coupled branch."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .backend import configure_backend
from .boundary_characterization import InteriorWindow, _interior_indices, _interior_tensor
from .config import BranchSmokeConfig, HarnessConfig
from .convergence import observed_order_from_three, step4_preconditioner_config, validate_refinement_ladder
from .coupled_branch import (
    CoupledFields,
    _create_branch_grid,
    _matter_number_current,
    boundary_sponge_profile_torch,
    branch_boundary_conditions,
    confinement_potential_torch,
    localization_weight_torch,
    resample_branch_state,
    run_branch_continuation,
    tensor_center_gradient_r,
    tensor_center_gradient_w,
    unpack_coupled_fields,
)
from .grid import TensorProductGrid
from .operators import (
    tensor_flux_divergence,
    tensor_vector_face_fluxes,
    tensor_weighted_gradient_fluxes,
)


NULL_DIAGNOSTIC = "null diagnostic"
SOLVER_FLOOR_DIAGNOSTIC = "solver-floor diagnostic"


@dataclass(frozen=True)
class ConservationRegion:
    r_start: int
    r_stop: int
    w_start: int
    w_stop: int
    label: str

    def to_dict(self, grid: TensorProductGrid) -> dict[str, Any]:
        return {
            "label": self.label,
            "r_start": self.r_start,
            "r_stop": self.r_stop,
            "w_start": self.w_start,
            "w_stop": self.w_stop,
            "r_face_min": float(grid.r_faces[self.r_start].detach().cpu().item()),
            "r_face_max": float(grid.r_faces[self.r_stop].detach().cpu().item()),
            "w_face_min": float(grid.w_faces[self.w_start].detach().cpu().item()),
            "w_face_max": float(grid.w_faces[self.w_stop].detach().cpu().item()),
        }


@dataclass(frozen=True)
class ConservationDiagnosticsConfig:
    name: str = "step6_conservation_diagnostics"
    levels: tuple[tuple[int, int], ...] = ((6, 4), (12, 8), (24, 16))
    refinement_ratio: int = 2
    interior: InteriorWindow = field(default_factory=InteriorWindow)
    sponge_width: float = 0.4
    sponge_matter_strength: float = 100.0
    sponge_gauge_strength: float = 100.0
    sponge_power: int = 2
    gauss_r_faces: tuple[float, ...] = (2.0 / 3.0, 1.0)
    gauss_w_min: float = 0.4
    gauss_w_max: float = 0.8
    null_current_floor: float = 1.0e-14
    gauss_solver_floor_relative: float = 1.0e-11
    budget_closure_abs_tol: float = 1.0e-10
    sponge_interior_abs_tol: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "levels": list(self.levels),
            "refinement_ratio": self.refinement_ratio,
            "interior": self.interior.to_dict(),
            "sponge_width": self.sponge_width,
            "sponge_matter_strength": self.sponge_matter_strength,
            "sponge_gauge_strength": self.sponge_gauge_strength,
            "sponge_power": self.sponge_power,
            "gauss_r_faces": list(self.gauss_r_faces),
            "gauss_w_min": self.gauss_w_min,
            "gauss_w_max": self.gauss_w_max,
            "null_current_floor": self.null_current_floor,
            "gauss_solver_floor_relative": self.gauss_solver_floor_relative,
            "budget_closure_abs_tol": self.budget_closure_abs_tol,
            "sponge_interior_abs_tol": self.sponge_interior_abs_tol,
        }


def null_floor_label(signal_norm: float, *, floor: float = 1.0e-14) -> str:
    """Return the Step-4 null/floor vocabulary for a conservation channel."""

    return NULL_DIAGNOSTIC if abs(signal_norm) <= floor else SOLVER_FLOOR_DIAGNOSTIC


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _branch_for_step6(
    base: BranchSmokeConfig,
    study: ConservationDiagnosticsConfig,
    *,
    sponge_enabled: bool,
) -> BranchSmokeConfig:
    label = "sponge_on" if sponge_enabled else "sponge_off"
    return replace(
        base,
        name=f"{study.name}_{label}",
        sponge_enabled=sponge_enabled,
        sponge_width=study.sponge_width if sponge_enabled else 0.0,
        sponge_matter_strength=study.sponge_matter_strength if sponge_enabled else 0.0,
        sponge_gauge_strength=study.sponge_gauge_strength if sponge_enabled else 0.0,
        sponge_power=study.sponge_power,
        newton=replace(base.newton, preconditioner=step4_preconditioner_config()),
    )


def _region_from_interior(grid: TensorProductGrid, interior: InteriorWindow) -> ConservationRegion:
    r_idx, w_idx = _interior_indices(grid, interior)
    return ConservationRegion(
        r_start=int(r_idx[0].detach().cpu().item()),
        r_stop=int(r_idx[-1].detach().cpu().item()) + 1,
        w_start=int(w_idx[0].detach().cpu().item()),
        w_stop=int(w_idx[-1].detach().cpu().item()) + 1,
        label="step5_interior_window",
    )


def _first_face_at_or_above(faces: torch.Tensor, value: float) -> int:
    indices = torch.nonzero(faces >= value - 1.0e-12, as_tuple=True)[0]
    if int(indices.numel()) == 0:
        raise ValueError("no grid face lies above the requested lower bound")
    return int(indices[0].detach().cpu().item())


def _last_face_at_or_below(faces: torch.Tensor, value: float) -> int:
    indices = torch.nonzero(faces <= value + 1.0e-12, as_tuple=True)[0]
    if int(indices.numel()) == 0:
        raise ValueError("no grid face lies below the requested upper bound")
    return int(indices[-1].detach().cpu().item())


def _nearest_face(faces: torch.Tensor, value: float) -> int:
    index = int(torch.argmin(torch.abs(faces - value)).detach().cpu().item())
    spacing = float(torch.min(torch.diff(faces)).detach().cpu().item())
    if abs(float(faces[index].detach().cpu().item()) - value) > 0.51 * spacing:
        raise ValueError("requested surface is not represented by a nearby grid face")
    return index


def closed_interior_gauss_regions(
    grid: TensorProductGrid,
    interior: InteriorWindow,
    *,
    r_faces: tuple[float, ...] = (2.0 / 3.0, 1.0),
    w_min: float = 0.4,
    w_max: float = 0.8,
) -> list[ConservationRegion]:
    """Choose two nested closed control volumes whose faces lie in the interior."""

    r_start = _first_face_at_or_above(grid.r_faces, interior.r_min)
    w_start = _nearest_face(grid.w_faces, w_min)
    w_stop = _nearest_face(grid.w_faces, w_max)
    if r_start != 0:
        raise ValueError("the Step-6 Gauss regions assume the regular r=0 inner face")
    if w_stop - w_start < 1:
        raise ValueError("need at least one w cell inside the Gauss window")
    if not (interior.w_min <= w_min < w_max <= interior.w_max):
        raise ValueError("requested Gauss w faces must lie inside the interior window")

    stops = [_nearest_face(grid.r_faces, value) for value in r_faces]
    if len(set(stops)) < 2:
        raise ValueError("need at least two distinct radial Gauss surfaces")
    regions = []
    for index, r_stop in enumerate(stops):
        r_value = float(grid.r_faces[r_stop].detach().cpu().item())
        if not (interior.r_min < r_value <= interior.r_max):
            raise ValueError("requested Gauss radial face must lie inside the interior window")
        regions.append(
            ConservationRegion(
                r_start=r_start,
                r_stop=r_stop,
                w_start=w_start,
                w_stop=w_stop,
                label=f"nested_surface_{index}",
            )
        )
    return regions


def _region_flux_integral(
    radial_flux: torch.Tensor,
    w_flux: torch.Tensor,
    region: ConservationRegion,
) -> torch.Tensor:
    return (
        torch.sum(radial_flux[region.r_stop, region.w_start : region.w_stop])
        - torch.sum(radial_flux[region.r_start, region.w_start : region.w_stop])
        + torch.sum(w_flux[region.r_start : region.r_stop, region.w_stop])
        - torch.sum(w_flux[region.r_start : region.r_stop, region.w_start])
    )


def _region_volume_integral(
    values: torch.Tensor,
    grid: TensorProductGrid,
    region: ConservationRegion,
) -> torch.Tensor:
    local = values[region.r_start : region.r_stop, region.w_start : region.w_stop]
    volumes = grid.cell_volumes[region.r_start : region.r_stop, region.w_start : region.w_stop]
    return torch.sum(local * volumes)


def _region_volume(grid: TensorProductGrid, region: ConservationRegion) -> torch.Tensor:
    return torch.sum(
        grid.cell_volumes[region.r_start : region.r_stop, region.w_start : region.w_stop]
    )


def _face_average(values: torch.Tensor, axis: int) -> torch.Tensor:
    if axis == 0:
        faces = torch.empty(
            (values.shape[0] + 1, values.shape[1]),
            dtype=values.dtype,
            device=values.device,
        )
        faces[1:-1, :] = 0.5 * (values[:-1, :] + values[1:, :])
        faces[0, :] = values[0, :]
        faces[-1, :] = values[-1, :]
        return faces
    if axis == 1:
        faces = torch.empty(
            (values.shape[0], values.shape[1] + 1),
            dtype=values.dtype,
            device=values.device,
        )
        faces[:, 1:-1] = 0.5 * (values[:, :-1] + values[:, 1:])
        faces[:, 0] = values[:, 0]
        faces[:, -1] = values[:, -1]
        return faces
    raise ValueError("axis must be 0 or 1")


def independent_gauss_face_fluxes(
    a0: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Surface flux from center-gradient reconstruction, independent of the solver flux."""

    grad_r_centers = tensor_center_gradient_r(a0, grid)
    grad_w_centers = tensor_center_gradient_w(a0, grid)
    grad_r_faces = _face_average(grad_r_centers, axis=0)
    grad_w_faces = _face_average(grad_w_centers, axis=1)
    z_centers = localization_weight_torch(grid.w_centers, cfg)
    z_faces = localization_weight_torch(grid.w_faces, cfg)
    radial_flux = -grid.radial_face_areas[:, None] * grid.dw * z_centers[None, :] * grad_r_faces
    w_flux = -grid.radial_shell_volumes[:, None] * z_faces[None, :] * grad_w_faces
    return radial_flux, w_flux


def _operator_gauss_face_fluxes(
    a0: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
) -> tuple[torch.Tensor, torch.Tensor]:
    boundaries = branch_boundary_conditions(cfg)
    z_centers = localization_weight_torch(grid.w_centers, cfg)
    z_faces = localization_weight_torch(grid.w_faces, cfg)
    grad_radial_flux, grad_w_flux = tensor_weighted_gradient_fluxes(
        a0,
        z_centers,
        z_faces,
        grid,
        boundaries.a0_radial_outer,
        boundaries.a0_w_lower,
        boundaries.a0_w_upper,
    )
    return -grad_radial_flux, -grad_w_flux


def _weighted_l2(values: torch.Tensor, volumes: torch.Tensor) -> torch.Tensor:
    return torch.sqrt(torch.sum(values * values * volumes))


def _energy_density(
    fields: CoupledFields,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    *,
    eos_K: float,
) -> torch.Tensor:
    density = fields.psi_real**2 + fields.psi_imag**2
    alpha = cfg.gauge_charge / cfg.hbar
    grad_real_r = tensor_center_gradient_r(fields.psi_real, grid)
    grad_imag_r = tensor_center_gradient_r(fields.psi_imag, grid)
    grad_real_w = tensor_center_gradient_w(fields.psi_real, grid)
    grad_imag_w = tensor_center_gradient_w(fields.psi_imag, grid)
    cov_real_r = grad_real_r + alpha * fields.ar * fields.psi_imag
    cov_imag_r = grad_imag_r - alpha * fields.ar * fields.psi_real
    cov_real_w = grad_real_w + alpha * fields.aw * fields.psi_imag
    cov_imag_w = grad_imag_w - alpha * fields.aw * fields.psi_real
    kinetic = (cfg.hbar**2 / (2.0 * cfg.particle_mass)) * (
        cov_real_r**2 + cov_imag_r**2 + cov_real_w**2 + cov_imag_w**2
    )
    potential = confinement_potential_torch(grid, cfg) * density
    internal = 0.25 * eos_K * density**5

    z = localization_weight_torch(grid.w_centers, cfg)[None, :]
    grad_a0_r = tensor_center_gradient_r(fields.a0, grid)
    grad_a0_w = tensor_center_gradient_w(fields.a0, grid)
    f_rw = tensor_center_gradient_r(fields.aw, grid) - tensor_center_gradient_w(fields.ar, grid)
    em = (z / (2.0 * cfg.mu0)) * (grad_a0_r**2 + grad_a0_w**2 + f_rw**2)
    return kinetic + potential + internal + em


def _sponge_sink_densities(
    fields: CoupledFields,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
) -> dict[str, torch.Tensor]:
    profile = boundary_sponge_profile_torch(grid, cfg)
    matter_real = cfg.sponge_matter_strength * profile * fields.psi_real
    matter_imag = cfg.sponge_matter_strength * profile * fields.psi_imag
    gauge_a0 = cfg.sponge_gauge_strength * profile * fields.a0
    gauge_ar = cfg.sponge_gauge_strength * profile * fields.ar
    gauge_aw = cfg.sponge_gauge_strength * profile * fields.aw
    number = fields.psi_real * matter_real + fields.psi_imag * matter_imag
    gauge_quadratic = fields.a0 * gauge_a0 + fields.ar * gauge_ar + fields.aw * gauge_aw
    return {
        "number": number,
        "charge": cfg.gauge_charge * number,
        "energy": number + gauge_quadratic,
    }


def _local_residual_rows(
    *,
    mode: str,
    level_index: int,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    fields: CoupledFields,
    density: torch.Tensor,
    energy_density: torch.Tensor,
    number_divergence: torch.Tensor,
    charge_divergence: torch.Tensor,
    energy_divergence: torch.Tensor,
    number_current_norm: float,
    charge_current_norm: float,
    energy_flux_norm: float,
    study: ConservationDiagnosticsConfig,
) -> list[dict[str, Any]]:
    r_idx, w_idx = _interior_indices(grid, study.interior)
    volumes = _interior_tensor(grid.cell_volumes, r_idx, w_idx)
    interior_volume = float(torch.sum(volumes).detach().cpu().item())
    totals = {
        "number": float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
        "charge": abs(cfg.gauge_charge)
        * float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
        "energy": float(torch.sum(energy_density * grid.cell_volumes).detach().cpu().item()),
    }
    channels = {
        "number": (number_divergence, number_current_norm),
        "charge": (charge_divergence, charge_current_norm),
        "energy": (energy_divergence, energy_flux_norm),
    }
    rows = []
    for sector, (divergence, signal_norm) in channels.items():
        interior_values = _interior_tensor(divergence, r_idx, w_idx)
        l2_abs = float(_weighted_l2(interior_values, volumes).detach().cpu().item())
        linf_abs = float(torch.max(torch.abs(interior_values)).detach().cpu().item())
        aggregate = max(abs(totals[sector]), 1.0e-300)
        linf_reference = max(aggregate / max(interior_volume, 1.0e-300), 1.0e-300)
        rows.append(
            {
                "mode": mode,
                "level": level_index,
                "grid": f"nr_{grid.spec.nr}_nw_{grid.spec.nw}",
                "sector": sector,
                "interior_l2_abs": l2_abs,
                "interior_linf_abs": linf_abs,
                "interior_l2_relative": l2_abs / aggregate,
                "interior_linf_relative": linf_abs / linf_reference,
                "target_blind_reference_total": totals[sector],
                "transport_signal_l2": signal_norm,
                "label": null_floor_label(signal_norm, floor=study.null_current_floor),
            }
        )
    return rows


def _gauss_rows(
    *,
    mode: str,
    level_index: int,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    fields: CoupledFields,
    density: torch.Tensor,
    study: ConservationDiagnosticsConfig,
    reconstruction: str,
) -> list[dict[str, Any]]:
    if reconstruction == "independent_center_gradient":
        gauss_radial_flux, gauss_w_flux = independent_gauss_face_fluxes(fields.a0, grid, cfg)
    elif reconstruction == "operator_flux":
        gauss_radial_flux, gauss_w_flux = _operator_gauss_face_fluxes(fields.a0, grid, cfg)
    else:
        raise ValueError(f"unknown Gauss reconstruction {reconstruction!r}")
    source = cfg.mu0 * cfg.gauge_charge * density
    rows = []
    for region in closed_interior_gauss_regions(
        grid,
        study.interior,
        r_faces=study.gauss_r_faces,
        w_min=study.gauss_w_min,
        w_max=study.gauss_w_max,
    ):
        lhs = _region_flux_integral(gauss_radial_flux, gauss_w_flux, region)
        rhs = _region_volume_integral(source, grid, region)
        residual = lhs - rhs
        rhs_abs = float(abs(rhs.detach().cpu().item()))
        region_volume = float(_region_volume(grid, region).detach().cpu().item())
        rows.append(
            {
                "mode": mode,
                "level": level_index,
                "grid": f"nr_{grid.spec.nr}_nw_{grid.spec.nw}",
                "reconstruction": reconstruction,
                "surface": region.label,
                **region.to_dict(grid),
                "region_volume": region_volume,
                "surface_flux": float(lhs.detach().cpu().item()),
                "enclosed_mu0_charge": float(rhs.detach().cpu().item()),
                "absolute_residual": float(residual.detach().cpu().item()),
                "relative_residual": float(
                    (abs(residual) / max(rhs_abs, 1.0e-300)).detach().cpu().item()
                ),
            }
        )
    return rows


def _budget_rows(
    *,
    mode: str,
    level_index: int,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    fields: CoupledFields,
    density: torch.Tensor,
    energy_density: torch.Tensor,
    number_flux: tuple[torch.Tensor, torch.Tensor],
    charge_flux: tuple[torch.Tensor, torch.Tensor],
    energy_flux: tuple[torch.Tensor, torch.Tensor],
    number_divergence: torch.Tensor,
    charge_divergence: torch.Tensor,
    energy_divergence: torch.Tensor,
    study: ConservationDiagnosticsConfig,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    budget_region = _region_from_interior(grid, study.interior)
    sinks = _sponge_sink_densities(fields, grid, cfg)
    totals = {
        "number": torch.sum(density * grid.cell_volumes),
        "charge": cfg.gauge_charge * torch.sum(density * grid.cell_volumes),
        "energy": torch.sum(energy_density * grid.cell_volumes),
    }
    fluxes = {
        "number": number_flux,
        "charge": charge_flux,
        "energy": energy_flux,
    }
    divergences = {
        "number": number_divergence,
        "charge": charge_divergence,
        "energy": energy_divergence,
    }
    rows = []
    sponge_rows = []
    for sector in ("number", "charge", "energy"):
        sink = sinks[sector]
        radial_flux, w_flux = fluxes[sector]
        net_outward_flux = _region_flux_integral(radial_flux, w_flux, budget_region)
        interior_sponge = _region_volume_integral(sink, grid, budget_region)
        local_balance_residual = _region_volume_integral(
            divergences[sector] - sink,
            grid,
            budget_region,
        )
        closure = net_outward_flux - interior_sponge - local_balance_residual
        reference_total = max(abs(float(totals[sector].detach().cpu().item())), 1.0e-300)
        total_sponge = torch.sum(sink * grid.cell_volumes)
        rows.append(
            {
                "mode": mode,
                "level": level_index,
                "grid": f"nr_{grid.spec.nr}_nw_{grid.spec.nw}",
                "sector": sector,
                **budget_region.to_dict(grid),
                "net_outward_flux": float(net_outward_flux.detach().cpu().item()),
                "interior_sponge_absorbed": float(interior_sponge.detach().cpu().item()),
                "interior_local_balance_residual": float(
                    local_balance_residual.detach().cpu().item()
                ),
                "closure_absolute": float(closure.detach().cpu().item()),
                "closure_relative": float(
                    (abs(closure) / reference_total).detach().cpu().item()
                ),
                "target_blind_reference_total": float(totals[sector].detach().cpu().item()),
            }
        )
        sponge_rows.append(
            {
                "mode": mode,
                "level": level_index,
                "grid": f"nr_{grid.spec.nr}_nw_{grid.spec.nw}",
                "sector": sector,
                "interior_absorbed": float(interior_sponge.detach().cpu().item()),
                "total_absorbed": float(total_sponge.detach().cpu().item()),
                "outer_layer_absorbed": float((total_sponge - interior_sponge).detach().cpu().item()),
                "interior_zero": (
                    abs(float(interior_sponge.detach().cpu().item()))
                    <= study.sponge_interior_abs_tol
                ),
            }
        )
    return rows, sponge_rows


def evaluate_conservation_state(
    *,
    mode: str,
    level_index: int,
    state: torch.Tensor,
    grid: TensorProductGrid,
    cfg: BranchSmokeConfig,
    study: ConservationDiagnosticsConfig,
) -> dict[str, Any]:
    fields, chemical_potential = unpack_coupled_fields(state, grid, has_chemical_potential=True)
    assert chemical_potential is not None
    density = fields.psi_real**2 + fields.psi_imag**2
    eos_K = cfg.continuation_K_values[-1]
    e_density = _energy_density(fields, grid, cfg, eos_K=eos_K)

    jr_number, jw_number = _matter_number_current(fields, grid, cfg)
    number_radial_flux, number_w_flux = tensor_vector_face_fluxes(jr_number, jw_number, grid)
    number_divergence = tensor_flux_divergence(number_radial_flux, number_w_flux, grid)
    charge_divergence = cfg.gauge_charge * number_divergence
    charge_flux = (cfg.gauge_charge * number_radial_flux, cfg.gauge_charge * number_w_flux)

    energy_radial_flux = chemical_potential * number_radial_flux
    energy_w_flux = chemical_potential * number_w_flux
    energy_divergence = tensor_flux_divergence(energy_radial_flux, energy_w_flux, grid)

    number_current_norm = float(
        _weighted_l2(torch.sqrt(jr_number**2 + jw_number**2), grid.cell_volumes)
        .detach()
        .cpu()
        .item()
    )
    charge_current_norm = abs(cfg.gauge_charge) * number_current_norm
    energy_flux_norm = float(
        _weighted_l2(
            torch.sqrt((chemical_potential * jr_number) ** 2 + (chemical_potential * jw_number) ** 2),
            grid.cell_volumes,
        )
        .detach()
        .cpu()
        .item()
    )

    local_rows = _local_residual_rows(
        mode=mode,
        level_index=level_index,
        grid=grid,
        cfg=cfg,
        fields=fields,
        density=density,
        energy_density=e_density,
        number_divergence=number_divergence,
        charge_divergence=charge_divergence,
        energy_divergence=energy_divergence,
        number_current_norm=number_current_norm,
        charge_current_norm=charge_current_norm,
        energy_flux_norm=energy_flux_norm,
        study=study,
    )
    gauss_rows = _gauss_rows(
        mode=mode,
        level_index=level_index,
        grid=grid,
        cfg=cfg,
        fields=fields,
        density=density,
        study=study,
        reconstruction="independent_center_gradient",
    )
    maxwell_residual_rows = _gauss_rows(
        mode=mode,
        level_index=level_index,
        grid=grid,
        cfg=cfg,
        fields=fields,
        density=density,
        study=study,
        reconstruction="operator_flux",
    )
    budget_rows, sponge_rows = _budget_rows(
        mode=mode,
        level_index=level_index,
        grid=grid,
        cfg=cfg,
        fields=fields,
        density=density,
        energy_density=e_density,
        number_flux=(number_radial_flux, number_w_flux),
        charge_flux=charge_flux,
        energy_flux=(energy_radial_flux, energy_w_flux),
        number_divergence=number_divergence,
        charge_divergence=charge_divergence,
        energy_divergence=energy_divergence,
        study=study,
    )
    return {
        "local_residual_rows": local_rows,
        "gauss_closure_rows": gauss_rows,
        "maxwell_residual_closure_rows": maxwell_residual_rows,
        "budget_rows": budget_rows,
        "sponge_rows": sponge_rows,
    }


def _add_group_orders(
    rows: list[dict[str, Any]],
    *,
    group_keys: tuple[str, ...],
    value_key: str,
    ratio: int,
    output_key: str = "observed_order",
    floor_value: float | None = None,
) -> None:
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(tuple(row[key] for key in group_keys), []).append(row)
    for group_rows in groups.values():
        group_rows.sort(key=lambda row: row["level"])
        group_is_floor = (
            floor_value is not None
            and max(abs(float(row[value_key])) for row in group_rows) <= floor_value
        )
        for index, row in enumerate(group_rows):
            row[f"{output_key}_note"] = SOLVER_FLOOR_DIAGNOSTIC if group_is_floor else None
            if index < 2 or group_is_floor:
                row[output_key] = None
            else:
                row[output_key] = observed_order_from_three(
                    group_rows[index - 2][value_key],
                    group_rows[index - 1][value_key],
                    row[value_key],
                    ratio,
                )


def _diagnostics_digest(results: dict[str, Any]) -> str:
    payload = {
        "local_residual_rows": results["local_residual_rows"],
        "gauss_closure_rows": results["gauss_closure_rows"],
        "maxwell_residual_closure_rows": results["maxwell_residual_closure_rows"],
        "budget_rows": results["budget_rows"],
        "sponge_rows": results["sponge_rows"],
        "pass_checks": results["pass_checks"],
        "identity_checks": results["identity_checks"],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _solve_ladder(
    *,
    config: HarnessConfig,
    study: ConservationDiagnosticsConfig,
    dtype: torch.dtype,
    sponge_enabled: bool,
) -> dict[str, Any]:
    mode = "sponge_on" if sponge_enabled else "sponge_off"
    branch = _branch_for_step6(config.branch, study, sponge_enabled=sponge_enabled)
    run_config = replace(config, branch=branch, run_root=str(Path(config.run_root) / mode))
    previous_grid: TensorProductGrid | None = None
    previous_state: torch.Tensor | None = None
    solve_rows = []
    local_rows: list[dict[str, Any]] = []
    gauss_rows: list[dict[str, Any]] = []
    budget_rows: list[dict[str, Any]] = []
    sponge_rows: list[dict[str, Any]] = []
    maxwell_rows: list[dict[str, Any]] = []

    for level_index, level in enumerate(study.levels):
        grid = _create_branch_grid(branch, level, dtype=dtype, device=run_config.backend.device)
        initial = (
            resample_branch_state(previous_state, previous_grid, grid, branch)
            if previous_state is not None and previous_grid is not None
            else None
        )
        state, summary = run_branch_continuation(
            run_config,
            grid,
            initial_state=initial,
            grid_name=f"{mode}_l{level_index}_nr_{level[0]}_nw_{level[1]}",
        )
        solve_rows.append(
            {
                "mode": mode,
                "level": level_index,
                "grid": f"nr_{level[0]}_nw_{level[1]}",
                "dof": summary["dof"],
                "final_residual_linf": summary["final_residual_linf"],
                "converged": summary["converged"],
                "message": summary["message"],
                "manifest": summary["manifest"],
            }
        )
        diagnostics = evaluate_conservation_state(
            mode=mode,
            level_index=level_index,
            state=state,
            grid=grid,
            cfg=branch,
            study=study,
        )
        local_rows.extend(diagnostics["local_residual_rows"])
        gauss_rows.extend(diagnostics["gauss_closure_rows"])
        maxwell_rows.extend(diagnostics["maxwell_residual_closure_rows"])
        budget_rows.extend(diagnostics["budget_rows"])
        sponge_rows.extend(diagnostics["sponge_rows"])
        if not summary["converged"]:
            break
        previous_grid = grid
        previous_state = state.detach()

    return {
        "solve_rows": solve_rows,
        "local_residual_rows": local_rows,
        "gauss_closure_rows": gauss_rows,
        "maxwell_residual_closure_rows": maxwell_rows,
        "budget_rows": budget_rows,
        "sponge_rows": sponge_rows,
    }


def _gauss_decreases(rows: list[dict[str, Any]]) -> bool:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        if row["mode"] != "sponge_on":
            continue
        grouped.setdefault((row["mode"], row["surface"]), []).append(row)
    if not grouped:
        return False
    for group_rows in grouped.values():
        group_rows.sort(key=lambda row: row["level"])
        values = [abs(row["relative_residual"]) for row in group_rows]
        if len(values) < 3 or not values[-1] < values[0]:
            return False
    return True


def _gauss_order_passes(rows: list[dict[str, Any]], *, minimum_order: float) -> bool:
    finest = [
        row
        for row in rows
        if row["mode"] == "sponge_on"
        and row["level"] == max(r["level"] for r in rows if r["mode"] == "sponge_on")
    ]
    return bool(finest) and all(
        row["observed_order"] is not None and row["observed_order"] >= minimum_order
        for row in finest
    )


def _maxwell_residual_below_solver_floor(
    rows: list[dict[str, Any]],
    solve_rows: list[dict[str, Any]],
) -> bool:
    floors = {
        (row["mode"], row["level"]): row["final_residual_linf"]
        for row in solve_rows
    }
    checked = False
    for row in rows:
        if row["mode"] != "sponge_on":
            continue
        floor = floors[(row["mode"], row["level"])] * row["region_volume"]
        residual = abs(row["absolute_residual"])
        if not (0.0 < residual <= floor):
            return False
        checked = True
    return checked


def run_step6(
    config: HarnessConfig | None = None,
    study: ConservationDiagnosticsConfig | None = None,
) -> dict[str, Any]:
    if config is None:
        config = HarnessConfig(
            run_root="software/stage1_solver/runs/step6_conservation_diagnostics",
            report_path="software/stage1_solver/reports/step6_conservation_diagnostics.md",
        )
    if study is None:
        study = ConservationDiagnosticsConfig()
    validate_refinement_ladder(study.levels, study.refinement_ratio)
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(config.backend)

    off = _solve_ladder(config=config, study=study, dtype=dtype, sponge_enabled=False)
    on = _solve_ladder(config=config, study=study, dtype=dtype, sponge_enabled=True)
    solve_rows = [*off["solve_rows"], *on["solve_rows"]]
    local_rows = [*off["local_residual_rows"], *on["local_residual_rows"]]
    gauss_rows = [*off["gauss_closure_rows"], *on["gauss_closure_rows"]]
    maxwell_rows = [*off["maxwell_residual_closure_rows"], *on["maxwell_residual_closure_rows"]]
    budget_rows = [*off["budget_rows"], *on["budget_rows"]]
    sponge_rows = [*off["sponge_rows"], *on["sponge_rows"]]

    _add_group_orders(
        local_rows,
        group_keys=("mode", "sector"),
        value_key="interior_l2_relative",
        ratio=study.refinement_ratio,
    )
    _add_group_orders(
        gauss_rows,
        group_keys=("mode", "reconstruction", "surface"),
        value_key="relative_residual",
        ratio=study.refinement_ratio,
    )
    _add_group_orders(
        maxwell_rows,
        group_keys=("mode", "reconstruction", "surface"),
        value_key="relative_residual",
        ratio=study.refinement_ratio,
        floor_value=study.gauss_solver_floor_relative,
    )
    _add_group_orders(
        budget_rows,
        group_keys=("mode", "sector"),
        value_key="closure_relative",
        ratio=study.refinement_ratio,
    )

    max_budget_closure = max((abs(row["closure_absolute"]) for row in budget_rows), default=float("inf"))
    identity_checks = {
        "fv_budget_identity_roundoff_not_a_physics_gate": (
            max_budget_closure <= study.budget_closure_abs_tol
        ),
        "sponge_support_excludes_interior_window_not_a_physics_gate": all(
            row["interior_zero"] for row in sponge_rows
        ),
    }
    pass_checks = {
        "minimum_three_levels": len(study.levels) >= 3,
        "all_solves_converged": all(row["converged"] for row in solve_rows),
        "two_nested_gauss_surfaces": (
            len({row["surface"] for row in gauss_rows if row["mode"] == "sponge_on"}) >= 2
        ),
        "independent_gauss_decreases_under_refinement": _gauss_decreases(gauss_rows),
        "independent_gauss_observed_order_at_least_one": _gauss_order_passes(
            gauss_rows,
            minimum_order=1.0,
        ),
        "operator_maxwell_residual_positive_below_solver_floor": (
            _maxwell_residual_below_solver_floor(maxwell_rows, solve_rows)
        ),
        "sponge_on_absorbs_outer_layer": all(
            row["outer_layer_absorbed"] > 0.0
            for row in sponge_rows
            if row["mode"] == "sponge_on"
        ),
        "null_current_sectors_labelled": all(
            row["label"] == NULL_DIAGNOSTIC
            for row in local_rows
            if row["mode"] == "sponge_on" and row["level"] == len(study.levels) - 1
        ),
    }
    results: dict[str, Any] = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "study": study.to_dict(),
        "method": {
            "stationary_framing": (
                "No time marching is performed. Conservation is diagnosed as "
                "bulk finite-volume divergence plus stationary budget closure."
            ),
            "levels": list(study.levels),
            "refinement_ratio": study.refinement_ratio,
            "interior_window": study.interior.to_dict(),
            "number_current": "stage1_solver.coupled_branch._matter_number_current",
            "divergence": (
                "area-weighted FV face fluxes from operators.tensor_vector_face_fluxes, "
                "operators.tensor_weighted_gradient_fluxes, and operators.tensor_flux_divergence"
            ),
            "gauss_reconstruction": (
                "independent center-gradient reconstruction of E=-grad(A0), averaged to "
                "closed-surface faces and multiplied by Z; differs from the solver's "
                "tensor_weighted_gradient_fluxes operator"
            ),
            "measure": "grid.cell_volumes, grid.radial_face_areas, grid.radial_shell_volumes",
            "sponge": "boundary_sponge_profile_torch and sponge_*_strength config fields",
            "eos": "U=(K/4)*rho^5; h=(5K/4)*rho^4 remains in physics.quintic_enthalpy",
            "virial_identity": "omitted; no independent term-by-term parent-action derivation included",
            "target_blind": (
                "No benchmark target values, extraction coefficients, or physical packet export are used."
            ),
        },
        "solve_rows": solve_rows,
        "local_residual_rows": local_rows,
        "gauss_closure_rows": gauss_rows,
        "maxwell_residual_closure_rows": maxwell_rows,
        "budget_rows": budget_rows,
        "sponge_rows": sponge_rows,
        "pass_checks": pass_checks,
        "identity_checks": identity_checks,
        "passed": all(pass_checks.values()),
    }
    results["diagnostics_digest"] = _diagnostics_digest(results)
    table_path = Path(config.run_root) / study.name / "conservation_diagnostics_table.json"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    results["machine_readable_table"] = str(table_path)
    table_path.write_text(
        json.dumps(results, indent=2, sort_keys=True, default=_json_default),
        encoding="utf-8",
    )
    return results


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.6e}"
    if isinstance(value, (tuple, list)):
        return "[" + ", ".join(str(item) for item in value) + "]"
    return str(value)


def _table(headers: list[str], rows: list[dict[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(header)) for header in headers) + " |")
    return "\n".join(lines)


def _finest_rows(rows: list[dict[str, Any]], *, mode: str) -> list[dict[str, Any]]:
    levels = [row["level"] for row in rows if row["mode"] == mode]
    if not levels:
        return []
    finest = max(levels)
    return [row for row in rows if row["mode"] == mode and row["level"] == finest]


def write_step6_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    study = results["study"]
    lines: list[str] = []
    lines.append("# Step 6 Conservation Diagnostics")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append(f"Diagnostics digest: `{results['diagnostics_digest']}`")
    lines.append("")
    lines.append(
        "**Stationary framing:** no time-marching drift loop is run. Conservation is "
        "reported as conservative finite-volume bulk divergence plus stationary budget "
        "closure, with explicit sponge accounting."
    )
    lines.append("")
    lines.append(
        "**Discipline:** engineering smoke, placeholder parameters, target-blind, no "
        "field-to-coefficient export, and no physical packet."
    )
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append(f"Levels: `{study['levels']}` with refinement ratio `{study['refinement_ratio']}`.")
    lines.append(f"Interior window: `{study['interior']}`.")
    lines.append(f"Number current: `{results['method']['number_current']}`.")
    lines.append(f"Divergence: {results['method']['divergence']}.")
    lines.append(f"Gauss reconstruction: {results['method']['gauss_reconstruction']}.")
    lines.append(f"Measures: {results['method']['measure']}.")
    lines.append(f"Sponge accounting: {results['method']['sponge']}.")
    lines.append(f"EOS: {results['method']['eos']}.")
    lines.append(f"Virial/Pohozaev identity: {results['method']['virial_identity']}.")
    lines.append("")
    lines.append("## Solve Rows")
    lines.append("")
    lines.append(
        _table(
            [
                "mode",
                "level",
                "grid",
                "dof",
                "final_residual_linf",
                "converged",
                "message",
            ],
            results["solve_rows"],
        )
    )
    lines.append("")
    lines.append("## Local Conservation Residuals")
    lines.append("")
    lines.append(
        "The mass, charge-current, and energy-flux sectors are null diagnostics on this "
        "isotropic stationary branch because the spatial current and spatial gauge "
        "field vanish by symmetry. The current-carrying conservation test is deferred "
        "to the WP3 tangent in step 8."
    )
    lines.append("")
    lines.append(
        _table(
            [
                "mode",
                "grid",
                "sector",
                "interior_l2_relative",
                "interior_linf_relative",
                "transport_signal_l2",
                "label",
                "observed_order",
            ],
            _finest_rows(results["local_residual_rows"], mode="sponge_on"),
        )
    )
    lines.append("")
    lines.append("## Gauss-Law Closure")
    lines.append("")
    lines.append(
        "This is an independent reconstruction check of the localized Gauss law: "
        "the closed-surface flux uses center-gradient `E=-grad(A0)` values averaged "
        "to faces, not the solver's weighted-gradient face flux. The comparison is "
        "`surface flux of Z F^{i0}` versus `mu0*q*int rho dV` on nested interior "
        "control volumes."
    )
    lines.append("")
    lines.append(
        _table(
            [
                "mode",
                "level",
                "grid",
                "surface",
                "surface_flux",
                "enclosed_mu0_charge",
                "absolute_residual",
                "relative_residual",
                "observed_order",
                "observed_order_note",
            ],
            results["gauss_closure_rows"],
        )
    )
    lines.append("")
    lines.append("Operator-flux Maxwell residual closure:")
    lines.append("")
    lines.append(
        "This second table is not independent: it restates the discrete Gauss law as "
        "the integrated Maxwell-a0 Newton residual over the same control volumes. "
        "It is retained only as a solver-floor check."
    )
    lines.append("")
    lines.append(
        _table(
            [
                "mode",
                "level",
                "grid",
                "surface",
                "absolute_residual",
                "relative_residual",
                "observed_order",
                "observed_order_note",
            ],
            [row for row in results["maxwell_residual_closure_rows"] if row["mode"] == "sponge_on"],
        )
    )
    lines.append("")
    lines.append("## Budgets")
    lines.append("")
    lines.append(
        "Budget closure here is a finite-volume consistency identity, not an "
        "independent physical conservation balance: `net outward flux = interior "
        "sponge absorbed + local balance residual` is exact to roundoff because the "
        "same conservative FV face fluxes generate both the divergence and the "
        "surface term. The interior sponge term is structurally zero because the "
        "configured sponge support lies outside the interior window; the real sponge "
        "measurement is the outer-layer absorbed amount in the Sponge Accounting table."
    )
    lines.append("")
    lines.append(
        _table(
            [
                "mode",
                "grid",
                "sector",
                "net_outward_flux",
                "interior_local_balance_residual",
                "closure_absolute",
                "closure_relative",
                "observed_order",
            ],
            [
                row
                for row in results["budget_rows"]
                if row["level"] == max(r["level"] for r in results["budget_rows"])
            ],
        )
    )
    lines.append("")
    lines.append("## Sponge Accounting")
    lines.append("")
    lines.append(
        _table(
            [
                "mode",
                "grid",
                "sector",
                "interior_absorbed",
                "outer_layer_absorbed",
                "total_absorbed",
                "interior_zero",
            ],
            [
                row
                for row in results["sponge_rows"]
                if row["level"] == max(r["level"] for r in results["sponge_rows"])
            ],
        )
    )
    lines.append("")
    lines.append(
        "The sponge-on/off split isolates the configured sponge terms. The interior "
        "absorbed amount is zero in every sector, which is the conservation-side "
        "restatement of Step 5's interior-zero sponge property."
    )
    lines.append("")
    lines.append("## Pass Checks")
    lines.append("")
    for key, value in results["pass_checks"].items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")
    lines.append("")
    lines.append("Identity checks, reported but not counted as physics gates:")
    for key, value in results["identity_checks"].items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")
    lines.append("")
    lines.append("## Limitations")
    lines.append("")
    lines.append(
        "Mass/charge spatial-current residuals and the energy-flux divergence are "
        "null/floor diagnostics here, not current-carrying conservation tests. The "
        "isotropic real-amplitude branch has no phase winding, no transverse vector "
        "potential, and therefore no spatial transport to conserve. The current-carrying "
        "test is intentionally deferred to step 8."
    )
    lines.append("")
    lines.append("## Reproduction")
    lines.append("")
    lines.append("```bash")
    lines.append(
        "timeout 600 env PYTHONPATH=software/stage1_solver/src "
        "python -m stage1_solver.conservation_harness"
    )
    lines.append(
        "timeout 600 env PYTHONPATH=software/stage1_solver/src "
        "pytest software/stage1_solver/tests/test_stage1_solver.py"
    )
    lines.append("```")
    lines.append("")
    lines.append(f"Machine-readable table: `{results['machine_readable_table']}`.")
    lines.append(
        "Target-blindness: no benchmark targets, extraction coefficients, or export "
        "packet paths are imported by the Step 6 diagnostics."
    )
    lines.append(
        "Export guard: the external physical export guard remains untouched; no "
        "physical packet is emitted."
    )
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path
