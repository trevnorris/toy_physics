"""Manufactured-solution benchmarks for Stage-1 production operators."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp
import torch

from .backend import configure_backend, tensor
from .boundaries import BoundaryCondition
from .config import HarnessConfig, RadialGridSpec, TensorGridSpec, WallGridSpec
from .grid import RadialGrid, TensorProductGrid, WallGrid
from .manifest import write_manifest
from .mms import run_convergence_study, weighted_l2_norm
from .newton import finite_difference_jvp_check
from .operators import (
    localized_maxwell_operator,
    radial_current,
    tensor_laplacian,
    wall_s_eta_operator,
)
from .physics import quintic_gnls_operator


def _torch_from_numpy(values: np.ndarray, *, dtype: torch.dtype, device: str) -> torch.Tensor:
    return tensor(np.asarray(values), dtype=dtype, device=device)


def _radial_symbols() -> tuple[sp.Symbol, sp.Symbol]:
    return sp.symbols("r R", positive=True, real=True)


def _matter_expressions() -> dict[str, Any]:
    r, radius = _radial_symbols()
    hbar, mass, eos_k, omega, chemical_potential = sp.symbols(
        "hbar mass eos_K omega chemical_potential", positive=True, real=True
    )
    compact_bump = (1 - (r / radius) ** 2) ** 4
    field = 1 + sp.Rational(11, 50) * compact_bump + sp.Rational(3, 100) * (
        r / radius
    ) ** 2 * compact_bump
    lap = sp.diff(field, r, 2) + 2 * sp.diff(field, r) / r
    forcing = (
        -(hbar**2 / (2 * mass)) * lap
        + sp.Rational(1, 2) * omega**2 * r**2 * field
        + (5 * eos_k / 4) * field**9
        - chemical_potential * field
    )
    return {
        "field": sp.lambdify((r, radius), field, "numpy"),
        "forcing": sp.lambdify(
            (r, radius, hbar, mass, eos_k, omega, chemical_potential),
            forcing,
            "numpy",
        ),
        "outer_derivative": sp.lambdify((radius,), sp.diff(field, r).subs(r, radius), "numpy"),
        "field_expr": field,
        "forcing_expr": sp.simplify(forcing),
    }


def _tensor_expressions() -> dict[str, Any]:
    r, radius = _radial_symbols()
    w, w_min, w_max = sp.symbols("w w_min w_max", real=True)
    length = w_max - w_min
    x = (w - w_min) / length
    radial_bump = (1 - (r / radius) ** 2) ** 4
    axial_bump = 256 * x**4 * (1 - x) ** 4
    radial = sp.Rational(1, 5) * radial_bump + sp.Rational(1, 25) * (
        r / radius
    ) ** 2 * radial_bump
    axial = sp.Rational(3, 20) * axial_bump
    field = radial + axial
    lap = sp.diff(radial, r, 2) + 2 * sp.diff(radial, r) / r + sp.diff(axial, w, 2)
    return {
        "field": sp.lambdify((r, w, radius, w_min, w_max), field, "numpy"),
        "forcing": sp.lambdify((r, w, radius, w_min, w_max), lap, "numpy"),
        "r_outer_derivative": sp.lambdify(
            (radius,), sp.diff(radial, r).subs(r, radius), "numpy"
        ),
        "w_lower_outward_derivative": sp.lambdify(
            (w_min, w_max), -sp.diff(axial, w).subs(w, w_min), "numpy"
        ),
        "w_upper_outward_derivative": sp.lambdify(
            (w_min, w_max), sp.diff(axial, w).subs(w, w_max), "numpy"
        ),
        "field_expr": field,
        "forcing_expr": sp.simplify(lap),
    }


def _current_expressions() -> dict[str, Any]:
    r, radius = _radial_symbols()
    hbar, mass = sp.symbols("hbar mass", positive=True, real=True)
    amplitude = 1 + sp.Rational(1, 20) * r**2 + sp.Rational(3, 25) * sp.cos(
        sp.pi * r**2 / radius**2
    )
    phase = sp.Rational(2, 5) * r + sp.Rational(17, 100) * sp.sin(
        sp.Rational(17, 10) * r
    )
    current = (hbar / mass) * amplitude**2 * sp.diff(phase, r)
    return {
        "amplitude": sp.lambdify((r, radius), amplitude, "numpy"),
        "phase": sp.lambdify((r,), phase, "numpy"),
        "current": sp.lambdify((r, radius, hbar, mass), current, "numpy"),
        "amplitude_expr": amplitude,
        "phase_expr": phase,
        "current_expr": sp.simplify(current),
    }


def _maxwell_expressions() -> dict[str, Any]:
    r, radius = _radial_symbols()
    w, w_min, w_max = sp.symbols("w w_min w_max", real=True)
    lam, xi = sp.symbols("lambda xi", positive=True, real=True)
    length = w_max - w_min
    x = (w - w_min) / length
    z = sp.exp(-(w / lam) ** 2)
    radial_bump = (1 - (r / radius) ** 2) ** 4
    axial_bump = 256 * x**4 * (1 - x) ** 4

    a0_r = sp.Rational(1, 5) * radial_bump
    a0_w = sp.Rational(3, 20) * axial_bump
    a0 = a0_r + a0_w

    ar = r * radial_bump * (sp.Rational(19, 100) + sp.Rational(7, 100) * axial_bump)
    aw = (
        radial_bump * (sp.Rational(11, 100) + sp.Rational(9, 100) * axial_bump)
        + sp.Rational(3, 50) * (r / radius) ** 2 * radial_bump * axial_bump
    )

    scalar_lap = z * (sp.diff(a0, r, 2) + 2 * sp.diff(a0, r) / r) + sp.diff(
        z * sp.diff(a0, w), w
    )
    divergence = sp.diff(r**2 * ar, r) / r**2 + sp.diff(aw, w)
    f_rw = sp.diff(aw, r) - sp.diff(ar, w)
    o0 = -scalar_lap
    o_r = -sp.diff(z * f_rw, w) + (1 / xi) * sp.diff(z * divergence, r)
    o_w = sp.diff(r**2 * z * f_rw, r) / r**2 + (1 / xi) * sp.diff(z * divergence, w)
    operator = sp.Matrix([o0, o_r, o_w])

    return {
        "a0": sp.lambdify((r, w, radius, w_min, w_max), a0, "numpy"),
        "ar": sp.lambdify((r, w, radius, w_min, w_max), ar, "numpy"),
        "aw": sp.lambdify((r, w, radius, w_min, w_max), aw, "numpy"),
        "operator": sp.lambdify((r, w, radius, w_min, w_max, lam, xi), operator, "numpy"),
        "weight": sp.lambdify((w, lam), z, "numpy"),
        "a0_r_outer_derivative": sp.lambdify(
            (radius,), sp.diff(a0_r, r).subs(r, radius), "numpy"
        ),
        "a0_w_lower_outward_derivative": sp.lambdify(
            (w_min, w_max), -sp.diff(a0_w, w).subs(w, w_min), "numpy"
        ),
        "a0_w_upper_outward_derivative": sp.lambdify(
            (w_min, w_max), sp.diff(a0_w, w).subs(w, w_max), "numpy"
        ),
        "a0_expr": a0,
        "ar_expr": ar,
        "aw_expr": aw,
        "operator_expr": operator,
    }


def _wall_expressions() -> dict[str, Any]:
    w, w_min, w_max = sp.symbols("w w_min w_max", real=True)
    ell = sp.symbols("ell", integer=True, nonnegative=True)
    t_w_base, t_w_amp = sp.symbols("T_w_base T_w_amp", positive=True, real=True)
    t_omega_base, t_omega_amp = sp.symbols(
        "T_Omega_base T_Omega_amp", positive=True, real=True
    )
    k_eta_base, k_eta_amp = sp.symbols("K_eta_base K_eta_amp", positive=True, real=True)
    length = w_max - w_min
    x = (w - w_min) / length
    bump = 256 * x**4 * (1 - x) ** 4
    eta = 1 + sp.Rational(1, 5) * bump
    t_w = t_w_base + t_w_amp * sp.sin(2 * sp.pi * x)
    t_omega = t_omega_base + t_omega_amp * sp.cos(2 * sp.pi * x)
    k_eta = k_eta_base + k_eta_amp * bump
    operator = -sp.diff(t_w * sp.diff(eta, w), w) + (
        k_eta + ell * (ell + 1) * t_omega
    ) * eta
    return {
        "eta": sp.lambdify((w, w_min, w_max), eta, "numpy"),
        "t_w": sp.lambdify((w, w_min, w_max, t_w_base, t_w_amp), t_w, "numpy"),
        "t_omega": sp.lambdify(
            (w, w_min, w_max, t_omega_base, t_omega_amp), t_omega, "numpy"
        ),
        "k_eta": sp.lambdify((w, w_min, w_max, k_eta_base, k_eta_amp), k_eta, "numpy"),
        "operator": sp.lambdify(
            (
                w,
                w_min,
                w_max,
                ell,
                t_w_base,
                t_w_amp,
                t_omega_base,
                t_omega_amp,
                k_eta_base,
                k_eta_amp,
            ),
            operator,
            "numpy",
        ),
        "lower_outward_derivative": sp.lambdify(
            (w_min, w_max), -sp.diff(eta, w).subs(w, w_min), "numpy"
        ),
        "upper_outward_derivative": sp.lambdify(
            (w_min, w_max), sp.diff(eta, w).subs(w, w_max), "numpy"
        ),
        "eta_expr": eta,
        "t_w_expr": t_w,
        "t_omega_expr": t_omega,
        "k_eta_expr": k_eta,
        "operator_expr": operator,
    }


_MATTER = _matter_expressions()
_TENSOR = _tensor_expressions()
_CURRENT = _current_expressions()
_MAXWELL = _maxwell_expressions()
_WALL = _wall_expressions()


def run_matter_mms(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.mms.matter
    full_config = config.to_dict()
    outer_derivative = float(_MATTER["outer_derivative"](cfg.r_max))
    outer_bc = BoundaryCondition.neumann(outer_derivative)

    def build_level(nr: int):
        grid = RadialGrid.create(
            RadialGridSpec(r_max=cfg.r_max, nr=nr),
            dtype=dtype,
            device=config.backend.device,
        )
        return grid, f"nr_{nr}", grid.dr, grid.to_dict(), grid.cell_volumes

    def evaluate_level(grid: RadialGrid):
        r = grid.numpy_centers()
        field_np = _MATTER["field"](r, cfg.r_max)
        forcing_np = _MATTER["forcing"](
            r,
            cfg.r_max,
            cfg.hbar,
            cfg.particle_mass,
            cfg.eos_K,
            cfg.trap_omega,
            cfg.chemical_potential,
        )
        field = _torch_from_numpy(field_np, dtype=dtype, device=config.backend.device)
        exact = _torch_from_numpy(forcing_np, dtype=dtype, device=config.backend.device)
        discrete = quintic_gnls_operator(field, grid, cfg, outer_bc)
        density = torch.real(torch.conj(field) * field)
        diagnostics = {
            "density_min": float(torch.min(density).detach().cpu().item()),
            "density_max": float(torch.max(density).detach().cpu().item()),
            "outer_boundary": outer_bc.to_dict(),
        }
        return discrete, exact, diagnostics

    result = run_convergence_study(
        name=cfg.name,
        description="Radial stationary gauged-GNLS operator with A=0 and h(rho)=5K*rho^4/4.",
        continuum_source="compact lines 556-583 and 638-648.",
        manufactured_field=str(_MATTER["field_expr"]),
        forcing_derivation="SymPy applied the continuum radial operator to the manufactured field.",
        levels=cfg.grid_levels,
        build_level=build_level,
        evaluate_level=evaluate_level,
        config=full_config,
        run_root=config.run_root,
        min_observed_order=cfg.min_observed_order,
        final_error_max=cfg.final_error_max,
        config_hash=config.config_hash(),
    )

    final_grid = RadialGrid.create(
        RadialGridSpec(r_max=cfg.r_max, nr=cfg.grid_levels[-1]),
        dtype=dtype,
        device=config.backend.device,
    )
    r = final_grid.numpy_centers()
    field = _torch_from_numpy(_MATTER["field"](r, cfg.r_max), dtype=dtype, device=config.backend.device)
    forcing = _torch_from_numpy(
        _MATTER["forcing"](
            r,
            cfg.r_max,
            cfg.hbar,
            cfg.particle_mass,
            cfg.eos_K,
            cfg.trap_omega,
            cfg.chemical_potential,
        ),
        dtype=dtype,
        device=config.backend.device,
    )
    residual_fn = lambda x: quintic_gnls_operator(x, final_grid, cfg, outer_bc) - forcing
    jacobian_check = finite_difference_jvp_check(
        residual_fn,
        field,
        epsilon=config.newton.finite_difference_jvp_epsilon,
        seed=config.jacobian_check_seed,
    )
    result_dict = asdict(result)
    result_dict["jacobian_check"] = jacobian_check
    result_dict["pass_checks"]["jacobian_relative"] = (
        jacobian_check["relative_residual"] <= config.jacobian_check_rel_tol
    )
    result_dict["pass_checks"]["jacobian_absolute"] = (
        jacobian_check["absolute_residual"] <= config.jacobian_check_abs_tol
    )
    result_dict["passed"] = all(result_dict["pass_checks"].values())
    return result_dict


def run_tensor_laplacian_mms(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.mms.tensor
    full_config = config.to_dict()
    r_outer = float(_TENSOR["r_outer_derivative"](cfg.r_max))
    w_lower = float(_TENSOR["w_lower_outward_derivative"](cfg.w_min, cfg.w_max))
    w_upper = float(_TENSOR["w_upper_outward_derivative"](cfg.w_min, cfg.w_max))
    radial_bc = BoundaryCondition.neumann(r_outer)
    w_lower_bc = BoundaryCondition.neumann(w_lower)
    w_upper_bc = BoundaryCondition.neumann(w_upper)

    def build_level(level: tuple[int, int]):
        nr, nw = level
        grid = TensorProductGrid.create(
            TensorGridSpec(r_max=cfg.r_max, nr=nr, w_min=cfg.w_min, w_max=cfg.w_max, nw=nw),
            dtype=dtype,
            device=config.backend.device,
        )
        spacing = max(grid.dr, grid.dw)
        return grid, f"nr_{nr}_nw_{nw}", spacing, grid.to_dict(), grid.cell_volumes

    def evaluate_level(grid: TensorProductGrid):
        rr, ww = np.meshgrid(grid.r_centers.detach().cpu().numpy(), grid.w_centers.detach().cpu().numpy(), indexing="ij")
        field_np = _TENSOR["field"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        forcing_np = _TENSOR["forcing"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        field = _torch_from_numpy(field_np, dtype=dtype, device=config.backend.device)
        exact = _torch_from_numpy(forcing_np, dtype=dtype, device=config.backend.device)
        discrete = tensor_laplacian(field, grid, radial_bc, w_lower_bc, w_upper_bc)
        diagnostics = {
            "radial_outer_bc": radial_bc.to_dict(),
            "w_lower_bc": w_lower_bc.to_dict(),
            "w_upper_bc": w_upper_bc.to_dict(),
        }
        return discrete, exact, diagnostics

    return asdict(
        run_convergence_study(
            name=cfg.name,
            description="Genuine 2D tensor-grid scalar Laplacian with nontrivial radial and w dependence and all face boundary operators applied.",
            continuum_source="Step-1 FV tensor operator carrying compact/prereg full radial measure convention.",
            manufactured_field=str(_TENSOR["field_expr"]),
            forcing_derivation="SymPy evaluated r^-2*d_r(r^2*d_r u)+d_w^2 u.",
            levels=cfg.grid_levels,
            build_level=build_level,
            evaluate_level=evaluate_level,
            config=full_config,
            run_root=config.run_root,
            min_observed_order=cfg.min_observed_order,
            final_error_max=cfg.final_error_max,
            config_hash=config.config_hash(),
        )
    )


def run_current_mms(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.mms.current
    full_config = config.to_dict()

    def build_level(nr: int):
        grid = RadialGrid.create(
            RadialGridSpec(r_max=cfg.r_max, nr=nr),
            dtype=dtype,
            device=config.backend.device,
        )
        return grid, f"nr_{nr}", grid.dr, grid.to_dict(), grid.cell_volumes

    def evaluate_level(grid: RadialGrid):
        r = grid.numpy_centers()
        amplitude = _CURRENT["amplitude"](r, cfg.r_max)
        phase = _CURRENT["phase"](r)
        psi_np = amplitude * np.exp(1j * phase)
        current_np = _CURRENT["current"](r, cfg.r_max, cfg.hbar, cfg.particle_mass)
        psi = _torch_from_numpy(psi_np, dtype=torch.complex128, device=config.backend.device)
        exact = _torch_from_numpy(current_np, dtype=dtype, device=config.backend.device)
        discrete = radial_current(
            psi,
            grid,
            hbar=cfg.hbar,
            particle_mass=cfg.particle_mass,
            gauge_charge=cfg.gauge_charge,
            gauge_potential=None,
        )
        diagnostics = {
            "gauge_field": "A_r=0",
            "current_l2": float(weighted_l2_norm(exact, grid.cell_volumes).detach().cpu().item()),
        }
        return discrete, exact, diagnostics

    result = run_convergence_study(
        name=cfg.name,
        description="Complex radial manufactured field with nonzero A=0 current.",
        continuum_source="compact lines 651-659.",
        manufactured_field=f"amplitude={_CURRENT['amplitude_expr']}; phase={_CURRENT['phase_expr']}",
        forcing_derivation="SymPy evaluated (hbar/m)*amplitude^2*d_r phase.",
        levels=cfg.grid_levels,
        build_level=build_level,
        evaluate_level=evaluate_level,
        config=full_config,
        run_root=config.run_root,
        min_observed_order=cfg.min_observed_order,
        final_error_max=cfg.final_error_max,
        config_hash=config.config_hash(),
    )
    result_dict = asdict(result)
    final_current_norm = result_dict["rows"][-1]["current_l2"]
    result_dict["pass_checks"]["nonzero_current"] = final_current_norm >= cfg.min_current_norm
    result_dict["passed"] = all(result_dict["pass_checks"].values())
    return result_dict


def run_maxwell_mms(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    cfg = config.mms.maxwell
    full_config = config.to_dict()
    a0_r_outer = float(_MAXWELL["a0_r_outer_derivative"](cfg.r_max))
    a0_w_lower = float(_MAXWELL["a0_w_lower_outward_derivative"](cfg.w_min, cfg.w_max))
    a0_w_upper = float(_MAXWELL["a0_w_upper_outward_derivative"](cfg.w_min, cfg.w_max))
    radial_bc = BoundaryCondition.neumann(a0_r_outer)
    w_lower_bc = BoundaryCondition.neumann(a0_w_lower)
    w_upper_bc = BoundaryCondition.neumann(a0_w_upper)

    def build_level(level: tuple[int, int]):
        nr, nw = level
        grid = TensorProductGrid.create(
            TensorGridSpec(r_max=cfg.r_max, nr=nr, w_min=cfg.w_min, w_max=cfg.w_max, nw=nw),
            dtype=dtype,
            device=config.backend.device,
        )
        spacing = max(grid.dr, grid.dw)
        return grid, f"nr_{nr}_nw_{nw}", spacing, grid.to_dict(), grid.cell_volumes

    def evaluate_level(grid: TensorProductGrid):
        r_np = grid.r_centers.detach().cpu().numpy()
        w_np = grid.w_centers.detach().cpu().numpy()
        rr, ww = np.meshgrid(r_np, w_np, indexing="ij")
        a0_np = _MAXWELL["a0"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        ar_np = _MAXWELL["ar"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        aw_np = _MAXWELL["aw"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max)
        exact_np = np.asarray(
            _MAXWELL["operator"](rr, ww, cfg.r_max, cfg.w_min, cfg.w_max, cfg.localization_width, cfg.xi),
            dtype=np.float64,
        )
        exact_np = np.squeeze(exact_np)
        weight_centers = _torch_from_numpy(
            _MAXWELL["weight"](w_np, cfg.localization_width),
            dtype=dtype,
            device=config.backend.device,
        )
        weight_faces = _torch_from_numpy(
            _MAXWELL["weight"](grid.w_faces.detach().cpu().numpy(), cfg.localization_width),
            dtype=dtype,
            device=config.backend.device,
        )
        discrete = localized_maxwell_operator(
            _torch_from_numpy(a0_np, dtype=dtype, device=config.backend.device),
            _torch_from_numpy(ar_np, dtype=dtype, device=config.backend.device),
            _torch_from_numpy(aw_np, dtype=dtype, device=config.backend.device),
            grid,
            xi=cfg.xi,
            weight_w_centers=weight_centers,
            weight_w_faces=weight_faces,
            a0_radial_outer_bc=radial_bc,
            a0_w_lower_bc=w_lower_bc,
            a0_w_upper_bc=w_upper_bc,
        )
        exact = _torch_from_numpy(exact_np, dtype=dtype, device=config.backend.device)
        diagnostics = {
            "gauge_fixing": "H=Z",
            "a0_radial_outer_bc": radial_bc.to_dict(),
            "a0_w_lower_bc": w_lower_bc.to_dict(),
            "a0_w_upper_bc": w_upper_bc.to_dict(),
        }
        return discrete, exact, diagnostics

    return asdict(
        run_convergence_study(
            name=cfg.name,
            description="Stationary axisymmetric localized Maxwell operator on (r,w), components A0/Ar/Aw, H=Z.",
            continuum_source="compact lines 590-630 and 674-689; prereg section D gauge and mixed-channel rows.",
            manufactured_field=(
                f"A0={_MAXWELL['a0_expr']}; Ar={_MAXWELL['ar_expr']}; Aw={_MAXWELL['aw_expr']}"
            ),
            forcing_derivation="SymPy evaluated the displayed H=Z Maxwell equation after stationary axisymmetric reduction.",
            levels=cfg.grid_levels,
            build_level=build_level,
            evaluate_level=evaluate_level,
            config=full_config,
            run_root=config.run_root,
            min_observed_order=cfg.min_observed_order,
            final_error_max=cfg.final_error_max,
            config_hash=config.config_hash(),
        )
    )


def run_wall_mms(config: HarnessConfig) -> dict[str, Any]:
    cfg = config.mms.wall
    full_config = config.to_dict()
    lower_derivative = float(_WALL["lower_outward_derivative"](cfg.w_min, cfg.w_max))
    upper_derivative = float(_WALL["upper_outward_derivative"](cfg.w_min, cfg.w_max))
    lower_bc = BoundaryCondition.neumann(lower_derivative)
    upper_bc = BoundaryCondition.neumann(upper_derivative)

    def t_w_values(w_values: np.ndarray) -> np.ndarray:
        return _WALL["t_w"](
            w_values,
            cfg.w_min,
            cfg.w_max,
            cfg.mms_only_placeholder_t_w_base,
            cfg.mms_only_placeholder_t_w_sine_amp,
        )

    def t_omega_values(w_values: np.ndarray) -> np.ndarray:
        return _WALL["t_omega"](
            w_values,
            cfg.w_min,
            cfg.w_max,
            cfg.mms_only_placeholder_t_omega_base,
            cfg.mms_only_placeholder_t_omega_cosine_amp,
        )

    def k_eta_values(w_values: np.ndarray) -> np.ndarray:
        return _WALL["k_eta"](
            w_values,
            cfg.w_min,
            cfg.w_max,
            cfg.mms_only_placeholder_k_eta_base,
            cfg.mms_only_placeholder_k_eta_bump_amp,
        )

    def build_level(nw: int):
        grid = WallGrid.create(
            WallGridSpec(w_min=cfg.w_min, w_max=cfg.w_max, nw=nw),
            dtype=dtype,
            device=config.backend.device,
        )
        return grid, f"nw_{nw}", grid.dw, grid.to_dict(), grid.cell_widths

    dtype = configure_backend(config.backend)

    def evaluate_level(grid: WallGrid):
        w_centers = grid.w_centers.detach().cpu().numpy()
        w_faces = grid.w_faces.detach().cpu().numpy()
        eta_np = _WALL["eta"](w_centers, cfg.w_min, cfg.w_max)
        exact_np = _WALL["operator"](
            w_centers,
            cfg.w_min,
            cfg.w_max,
            cfg.spherical_l,
            cfg.mms_only_placeholder_t_w_base,
            cfg.mms_only_placeholder_t_w_sine_amp,
            cfg.mms_only_placeholder_t_omega_base,
            cfg.mms_only_placeholder_t_omega_cosine_amp,
            cfg.mms_only_placeholder_k_eta_base,
            cfg.mms_only_placeholder_k_eta_bump_amp,
        )
        eta = _torch_from_numpy(eta_np, dtype=dtype, device=config.backend.device)
        exact = _torch_from_numpy(exact_np, dtype=dtype, device=config.backend.device)
        discrete = wall_s_eta_operator(
            eta,
            grid,
            t_w_faces=_torch_from_numpy(t_w_values(w_faces), dtype=dtype, device=config.backend.device),
            t_omega_centers=_torch_from_numpy(
                t_omega_values(w_centers), dtype=dtype, device=config.backend.device
            ),
            k_eta_centers=_torch_from_numpy(
                k_eta_values(w_centers), dtype=dtype, device=config.backend.device
            ),
            spherical_l=cfg.spherical_l,
            lower_bc=lower_bc,
            upper_bc=upper_bc,
        )
        diagnostics = {
            "spherical_l": cfg.spherical_l,
            "ell_factor": cfg.spherical_l * (cfg.spherical_l + 1),
            "lower_bc": lower_bc.to_dict(),
            "upper_bc": upper_bc.to_dict(),
        }
        return discrete, exact, diagnostics

    result = asdict(
        run_convergence_study(
            name=cfg.name,
            description="Stationary modal wall operator from S_eta^(2) in densitized flat-dw convention.",
            continuum_source=(
                "research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md "
                "lines 198-225."
            ),
            manufactured_field=str(_WALL["eta_expr"]),
            forcing_derivation=(
                "SymPy varied the stationary modal action to "
                "-d_w(T_w d_w eta)+[K_eta+ell(ell+1)T_Omega]eta."
            ),
            levels=cfg.grid_levels,
            build_level=build_level,
            evaluate_level=evaluate_level,
            config=full_config,
            run_root=config.run_root,
            min_observed_order=cfg.min_observed_order,
            final_error_max=cfg.final_error_max,
            config_hash=config.config_hash(),
        )
    )
    result["placeholder_coefficients"] = {
        "label": (
            "MMS-only structural-certification placeholders; NOT the physical wall packet, "
            "which is `free_choice` (prereg §E) and is frozen per-run at solve time."
        ),
        "mu_eta": cfg.mms_only_placeholder_mu_eta,
        "T_w(w)": str(_WALL["t_w_expr"]),
        "T_Omega(w)": str(_WALL["t_omega_expr"]),
        "K_eta(w)": str(_WALL["k_eta_expr"]),
    }
    result["methodology_note"] = (
        "Previous STOP resolved: MMS certifies discretization of the pinned form; "
        "these placeholders are not physical constitutive values."
    )
    return result


def run_all_mms_benchmarks(config: HarnessConfig) -> dict[str, Any]:
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    matter = run_matter_mms(config)
    tensor_result = run_tensor_laplacian_mms(config)
    current = run_current_mms(config)
    maxwell = run_maxwell_mms(config)
    wall = run_wall_mms(config)
    sections = {
        "matter": matter,
        "tensor": tensor_result,
        "current": current,
        "maxwell": maxwell,
        "wall": wall,
    }
    return {
        "sections": sections,
        "passed": all(section["passed"] for section in sections.values()),
    }
