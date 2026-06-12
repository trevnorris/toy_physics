"""Known-answer benchmark suite and validation gates."""

from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .backend import configure_backend, tensor
from .boundaries import BoundaryCondition
from .config import HarnessConfig, RadialGridSpec
from .grid import RadialGrid
from .manifest import write_manifest
from .newton import finite_difference_jvp_check, solve_newton_jvp
from .operators import integrate, max_abs_radial_current, radial_fluxes, weighted_l2_error
from .physics import cubic_gpe_residual, gaussian_initial_state
from .references import (
    discrete_linear_harmonic_ground_state,
    linear_harmonic_exact_field,
    linear_harmonic_robin_bc,
    reference_tensor,
    solve_cubic_gpe_reference,
)


def observed_orders(errors: list[float], spacings: list[float]) -> list[float | None]:
    orders: list[float | None] = [None]
    for i in range(1, len(errors)):
        if errors[i] <= 0.0 or errors[i - 1] <= 0.0:
            orders.append(None)
        else:
            orders.append(math.log(errors[i - 1] / errors[i]) / math.log(spacings[i - 1] / spacings[i]))
    return orders


def run_linear_benchmark(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    rows: list[dict[str, Any]] = []
    eigen_errors: list[float] = []
    field_errors: list[float] = []
    spacings: list[float] = []
    bc = linear_harmonic_robin_bc(config.linear)
    full_config = config.to_dict()

    for nr in config.linear.grid_levels:
        grid = RadialGrid.create(
            RadialGridSpec(r_max=config.linear.r_max, nr=nr),
            dtype=dtype,
            device=config.backend.device,
        )
        mu, field_np, discrete_residual = discrete_linear_harmonic_ground_state(grid, config.linear)
        field = tensor(field_np, dtype=dtype, device=config.backend.device)
        exact_np = linear_harmonic_exact_field(grid.numpy_centers(), config.linear)
        exact = tensor(exact_np, dtype=dtype, device=config.backend.device)
        exact = exact / torch.sqrt(integrate(exact * exact, grid))
        if float(integrate(field * exact, grid).detach().cpu().item()) < 0.0:
            field = -field
        eigen_error = abs(mu - config.linear.exact_mu)
        field_error = float(weighted_l2_error(field, exact, grid).detach().cpu().item())
        mass = float(integrate(field * field, grid).detach().cpu().item())
        fluxes = radial_fluxes(field, grid, bc)
        origin_flux = float(abs(fluxes[0].detach().cpu().item()))
        current_abs = max_abs_radial_current(field, grid)
        row = {
            "nr": nr,
            "dr": grid.dr,
            "computed_mu": mu,
            "exact_mu": config.linear.exact_mu,
            "eigenvalue_error": eigen_error,
            "field_l2_error": field_error,
            "discrete_eigen_residual_linf": discrete_residual,
            "mass_drift": abs(mass - 1.0),
            "origin_flux_abs": origin_flux,
            "current_max_abs": current_abs,
            "boundary": bc.to_dict(),
        }
        manifest = write_manifest(
            run_root=config.run_root,
            benchmark_name=config.linear.name,
            grid_name=f"nr_{nr}",
            config=full_config,
            mesh=grid.to_dict(),
            results=row,
        )
        row["manifest"] = str(manifest)
        rows.append(row)
        eigen_errors.append(eigen_error)
        field_errors.append(field_error)
        spacings.append(grid.dr)

    eigen_orders = observed_orders(eigen_errors, spacings)
    field_orders = observed_orders(field_errors, spacings)
    for row, eig_order, fld_order in zip(rows, eigen_orders, field_orders):
        row["eigenvalue_order"] = eig_order
        row["field_l2_order"] = fld_order

    final = rows[-1]
    pass_checks = {
        "eigen_order": min(o for o in eigen_orders[1:] if o is not None) >= config.linear.min_observed_order,
        "field_order": min(o for o in field_orders[1:] if o is not None) >= config.linear.min_observed_order,
        "final_eigenvalue_error": final["eigenvalue_error"] <= config.linear.final_eigenvalue_error_max,
        "final_field_l2_error": final["field_l2_error"] <= config.linear.final_field_l2_error_max,
        "mass_drift": final["mass_drift"] <= 5.0e-14,
        "current_abs": final["current_max_abs"] <= 1.0e-14,
    }
    return {
        "name": config.linear.name,
        "description": "Linear radial harmonic oscillator ground state with exact Robin boundary.",
        "reference": "closed form exp(-omega*r^2/2), mu=3*omega/2",
        "rows": rows,
        "pass_checks": pass_checks,
        "passed": all(pass_checks.values()),
    }


def run_cubic_benchmark(config: HarnessConfig) -> dict[str, Any]:
    dtype = configure_backend(config.backend)
    reference = solve_cubic_gpe_reference(config.cubic)
    rows: list[dict[str, Any]] = []
    mu_errors: list[float] = []
    field_errors: list[float] = []
    spacings: list[float] = []
    final_residual_fn = None
    final_state = None
    full_config = config.to_dict()
    outer_bc = BoundaryCondition.dirichlet(config.cubic.outer_boundary_value)

    for nr in config.cubic.grid_levels:
        grid = RadialGrid.create(
            RadialGridSpec(r_max=config.cubic.r_max, nr=nr),
            dtype=dtype,
            device=config.backend.device,
        )
        residual_fn = lambda state, grid=grid: cubic_gpe_residual(state, grid, config.cubic, outer_bc)
        x0 = gaussian_initial_state(grid, config.cubic)
        newton = solve_newton_jvp(residual_fn, x0, config.newton)
        if not newton.converged:
            raise RuntimeError(f"Newton failed for nr={nr}: {newton.message}")
        state = newton.x.detach()
        field = state[:nr]
        mu = float(state[nr].detach().cpu().item())
        ref_field = reference_tensor(reference, grid)
        if float(integrate(field * ref_field, grid).detach().cpu().item()) < 0.0:
            field = -field
        mu_error = abs(mu - reference.mu)
        field_error = float(weighted_l2_error(field, ref_field, grid).detach().cpu().item())
        mass_residual = float(abs((integrate(field * field, grid) - config.cubic.mass).detach().cpu().item()))
        current_abs = max_abs_radial_current(field, grid)
        fluxes = radial_fluxes(field, grid, outer_bc)
        origin_flux = float(abs(fluxes[0].detach().cpu().item()))
        row = {
            "nr": nr,
            "dr": grid.dr,
            "computed_mu": mu,
            "reference_mu": reference.mu,
            "mu_error": mu_error,
            "field_l2_error": field_error,
            "newton_converged": newton.converged,
            "newton_iterations": newton.iterations,
            "newton_initial_residual_linf": newton.initial_residual_norm,
            "newton_final_residual_linf": newton.final_residual_norm,
            "newton_tolerance": newton.tolerance,
            "mass_drift": mass_residual,
            "origin_flux_abs": origin_flux,
            "current_max_abs": current_abs,
            "boundary": outer_bc.to_dict(),
        }
        manifest = write_manifest(
            run_root=config.run_root,
            benchmark_name=config.cubic.name,
            grid_name=f"nr_{nr}",
            config=full_config,
            mesh=grid.to_dict(),
            results=row,
        )
        row["manifest"] = str(manifest)
        rows.append(row)
        mu_errors.append(mu_error)
        field_errors.append(field_error)
        spacings.append(grid.dr)
        final_residual_fn = residual_fn
        final_state = state

    if final_residual_fn is None or final_state is None:
        raise RuntimeError("No cubic benchmark rows were produced")
    jacobian_check = finite_difference_jvp_check(
        final_residual_fn,
        final_state,
        epsilon=config.newton.finite_difference_jvp_epsilon,
        seed=config.jacobian_check_seed,
    )
    mu_orders = observed_orders(mu_errors, spacings)
    field_orders = observed_orders(field_errors, spacings)
    for row, mu_order, fld_order in zip(rows, mu_orders, field_orders):
        row["mu_order"] = mu_order
        row["field_l2_order"] = fld_order

    final = rows[-1]
    pass_checks = {
        "mu_order": min(o for o in mu_orders[1:] if o is not None) >= config.cubic.min_observed_order,
        "field_order": min(o for o in field_orders[1:] if o is not None) >= config.cubic.min_observed_order,
        "final_mu_error": final["mu_error"] <= config.cubic.final_mu_error_max,
        "final_field_l2_error": final["field_l2_error"] <= config.cubic.final_field_l2_error_max,
        "mass_drift": final["mass_drift"] <= config.cubic.max_mass_residual,
        "current_abs": final["current_max_abs"] <= config.cubic.max_current_abs,
        "jacobian_relative": jacobian_check["relative_residual"] <= config.jacobian_check_rel_tol,
        "jacobian_absolute": jacobian_check["absolute_residual"] <= config.jacobian_check_abs_tol,
    }
    return {
        "name": config.cubic.name,
        "description": "Mass-constrained stationary cubic GPE in a radial harmonic trap.",
        "reference": "independent SciPy solve_bvp radial BVP with unknown chemical potential",
        "reference_details": {
            "mu": reference.mu,
            "mass": reference.mass,
            "solver_nodes": reference.solver_nodes,
            "solver_status": reference.solver_status,
            "solver_message": reference.solver_message,
            "max_ode_residual_diagnostic": reference.max_ode_residual,
        },
        "rows": rows,
        "jacobian_check": jacobian_check,
        "pass_checks": pass_checks,
        "passed": all(pass_checks.values()),
    }


def run_all_benchmarks(config: HarnessConfig) -> dict[str, Any]:
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    linear = run_linear_benchmark(config)
    cubic = run_cubic_benchmark(config)
    return {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "linear": linear,
        "cubic": cubic,
        "passed": linear["passed"] and cubic["passed"],
    }
