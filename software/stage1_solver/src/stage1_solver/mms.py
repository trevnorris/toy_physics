"""Reusable manufactured-solution convergence harness."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable, Iterable

import torch

from .manifest import write_manifest


@dataclass(frozen=True)
class MMSLevelResult:
    grid: str
    spacing: float
    error: float
    reference_norm: float
    observed_order: float | None
    mesh: dict[str, Any]
    diagnostics: dict[str, Any]
    manifest: str


@dataclass(frozen=True)
class MMSResult:
    name: str
    description: str
    continuum_source: str
    manufactured_field: str
    forcing_derivation: str
    rows: list[dict[str, Any]]
    pass_checks: dict[str, bool]
    passed: bool


def observed_orders(errors: list[float], spacings: list[float]) -> list[float | None]:
    orders: list[float | None] = [None]
    for i in range(1, len(errors)):
        if errors[i] <= 0.0 or errors[i - 1] <= 0.0:
            orders.append(None)
        else:
            orders.append(math.log(errors[i - 1] / errors[i]) / math.log(spacings[i - 1] / spacings[i]))
    return orders


def _weighted_square(values: torch.Tensor, weights: torch.Tensor) -> torch.Tensor:
    magnitude = torch.real(torch.conj(values) * values)
    while weights.ndim < magnitude.ndim:
        weights = weights.unsqueeze(0)
    return magnitude * weights


def weighted_l2_norm(values: torch.Tensor, weights: torch.Tensor) -> torch.Tensor:
    return torch.sqrt(torch.sum(_weighted_square(values, weights)))


def run_convergence_study(
    *,
    name: str,
    description: str,
    continuum_source: str,
    manufactured_field: str,
    forcing_derivation: str,
    levels: Iterable[Any],
    build_level: Callable[[Any], tuple[Any, str, float, dict[str, Any], torch.Tensor]],
    evaluate_level: Callable[[Any], tuple[torch.Tensor, torch.Tensor, dict[str, Any]]],
    config: dict[str, Any],
    run_root: str,
    min_observed_order: float,
    final_error_max: float,
    config_hash: str,
) -> MMSResult:
    rows: list[dict[str, Any]] = []
    errors: list[float] = []
    spacings: list[float] = []
    pending: list[tuple[dict[str, Any], dict[str, Any]]] = []

    for level in levels:
        grid, grid_name, spacing, mesh, weights = build_level(level)
        discrete, exact, diagnostics = evaluate_level(grid)
        error_tensor = discrete - exact
        error = float(weighted_l2_norm(error_tensor, weights).detach().cpu().item())
        reference_norm = float(weighted_l2_norm(exact, weights).detach().cpu().item())
        row = {
            "grid": grid_name,
            "spacing": spacing,
            "error": error,
            "reference_norm": reference_norm,
            **diagnostics,
        }
        rows.append(row)
        pending.append((row, mesh))
        errors.append(error)
        spacings.append(spacing)

    orders = observed_orders(errors, spacings)
    for (row, mesh), order in zip(pending, orders):
        row["observed_order"] = order
        manifest = write_manifest(
            run_root=run_root,
            benchmark_name=name,
            grid_name=row["grid"],
            config=config,
            mesh=mesh,
            results=row,
            config_hash=config_hash,
        )
        row["manifest"] = str(manifest)

    finite_orders = [order for order in orders[1:] if order is not None]
    final = rows[-1]
    pass_checks = {
        "observed_order": bool(finite_orders) and min(finite_orders) >= min_observed_order,
        "final_error": final["error"] <= final_error_max,
    }
    return MMSResult(
        name=name,
        description=description,
        continuum_source=continuum_source,
        manufactured_field=manufactured_field,
        forcing_derivation=forcing_derivation,
        rows=rows,
        pass_checks=pass_checks,
        passed=all(pass_checks.values()),
    )
