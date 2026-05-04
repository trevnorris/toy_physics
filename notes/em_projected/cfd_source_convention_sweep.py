#!/usr/bin/env python3
"""Sweep source/current conventions for a runtime-monitor snapshot.

This is a diagnostic companion to cfd_runtime_monitor_postprocess.py. It does
not change the verdict. It asks whether the continuity/Poisson failures could
be explained by a simple sign, scale, uniform-offset, or current-vs-velocity
schema convention mismatch.
"""
from __future__ import annotations

import argparse
import math
import pathlib
from dataclasses import dataclass
from typing import Any

import numpy as np

from cfd_runtime_monitor_postprocess import (
    compute_s_rho,
    detect_source_schema,
    divergence3,
    gradient3,
    integrate_w,
    laplacian3,
    load_snapshot,
    rms,
    solve_phi_from_velocity,
    uniform_spacing,
)


@dataclass(frozen=True)
class Fields3D:
    rho: np.ndarray
    jx: np.ndarray
    jy: np.ndarray
    jz: np.ndarray
    source: np.ndarray
    dt_rho: np.ndarray
    phi: np.ndarray | None
    rho0: float
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray


def project_fields(data: dict[str, np.ndarray]) -> Fields3D:
    schema = detect_source_schema(data)
    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    z = np.asarray(data["z"], dtype=float)

    if schema == "full_4d":
        rho4 = np.asarray(data["rho"], dtype=float)
        jx4 = np.asarray(data["jx"], dtype=float)
        jy4 = np.asarray(data["jy"], dtype=float)
        jz4 = np.asarray(data["jz"], dtype=float)
        jw = np.asarray(data["jw"], dtype=float)
        w = np.asarray(data["w"], dtype=float)
        W = np.asarray(data["W"], dtype=float)
        dWdw = np.asarray(data.get("dWdw", np.gradient(W, w, edge_order=2)), dtype=float)
        rho = integrate_w(rho4, W, w)
        jx = integrate_w(jx4, W, w)
        jy = integrate_w(jy4, W, w)
        jz = integrate_w(jz4, W, w)
        source = compute_s_rho(jw, W, dWdw, w)
    else:
        rho = np.asarray(data["rho_brane"], dtype=float)
        jx = np.asarray(data["Jx_brane"], dtype=float)
        jy = np.asarray(data["Jy_brane"], dtype=float)
        jz = np.asarray(data["Jz_brane"], dtype=float)
        source = np.asarray(data["S_rho"], dtype=float)

    dt_rho = np.asarray(data.get("dt_rho", np.zeros_like(rho)), dtype=float)
    phi = np.asarray(data["phi3"], dtype=float) if "phi3" in data else None
    rho0 = float(data.get("rho0", float(np.mean(rho))))
    return Fields3D(rho, jx, jy, jz, source, dt_rho, phi, rho0, x, y, z)


def fit_scale(target: np.ndarray, source: np.ndarray) -> tuple[float, float]:
    denom = float(np.sum(source * source))
    if denom == 0.0:
        return (math.nan, math.nan)
    scale = float(np.sum(target * source) / denom)
    rel = rms(target - scale * source) / max(rms(source), 1e-30)
    return scale, rel


def fit_scale_offset(target: np.ndarray, source: np.ndarray) -> tuple[float, float, float]:
    n = float(source.size)
    sum_s = float(np.sum(source))
    sum_s2 = float(np.sum(source * source))
    sum_t = float(np.sum(target))
    sum_ts = float(np.sum(target * source))
    mat = np.array([[sum_s2, sum_s], [sum_s, n]], dtype=float)
    rhs = np.array([sum_ts, sum_t], dtype=float)
    try:
        scale, offset = np.linalg.solve(mat, rhs)
    except np.linalg.LinAlgError:
        return (math.nan, math.nan, math.nan)
    rel = rms(target - scale * source - offset) / max(rms(source), 1e-30)
    return float(scale), float(offset), rel


def summarize_interpretation(
    label: str,
    fields: Fields3D,
    *,
    stored_as_velocity: bool,
    periodic_xyz: bool,
) -> dict[str, Any]:
    dx = uniform_spacing(fields.x)
    dy = uniform_spacing(fields.y)
    dz = uniform_spacing(fields.z)
    spacings = (dx, dy, dz)

    if stored_as_velocity:
        vx = fields.jx
        vy = fields.jy
        vz = fields.jz
        jx = fields.rho * vx
        jy = fields.rho * vy
        jz = fields.rho * vz
    else:
        jx = fields.jx
        jy = fields.jy
        jz = fields.jz
        rho_floor = np.maximum(np.abs(fields.rho), 1e-12)
        vx = jx / rho_floor
        vy = jy / rho_floor
        vz = jz / rho_floor

    div_j = divergence3(jx, jy, jz, spacings, periodic_xyz)
    phi = fields.phi
    if phi is None:
        if not periodic_xyz:
            raise ValueError("phi3 is required for nonperiodic convention sweeps")
        phi = solve_phi_from_velocity(vx, vy, vz, spacings)

    grad_rho = gradient3(fields.rho, spacings, periodic_xyz)
    grad_phi = gradient3(phi, spacings, periodic_xyz)
    grad_rho_dot_v = grad_rho[0] * vx + grad_rho[1] * vy + grad_rho[2] * vz
    lap_phi = laplacian3(phi, spacings, periodic_xyz)

    source = fields.source
    cont_target = fields.dt_rho + div_j
    exact_target = fields.rho * lap_phi + fields.dt_rho + grad_rho_dot_v
    lin_target = fields.rho0 * lap_phi
    exact_target_phi_flip = -fields.rho * lap_phi + fields.dt_rho - grad_rho_dot_v
    lin_target_phi_flip = -fields.rho0 * lap_phi

    cont_scale, cont_rel_fit = fit_scale(cont_target, source)
    pois_scale, pois_rel_fit = fit_scale(exact_target, source)
    lin_scale, lin_rel_fit = fit_scale(lin_target, source)
    cont_scale_off, cont_offset, cont_rel_off = fit_scale_offset(cont_target, source)
    pois_scale_off, pois_offset, pois_rel_off = fit_scale_offset(exact_target, source)
    flip_scale, flip_rel_fit = fit_scale(exact_target_phi_flip, source)
    lin_flip_scale, lin_flip_rel_fit = fit_scale(lin_target_phi_flip, source)

    return {
        "label": label,
        "continuity_rel": rms(cont_target - source) / max(rms(source), 1e-30),
        "poisson_exact_rel": rms(exact_target - source) / max(rms(source), 1e-30),
        "poisson_lin_rel": rms(lin_target - source) / max(rms(source), 1e-30),
        "continuity_best_source_scale": cont_scale,
        "continuity_best_source_scale_rel": cont_rel_fit,
        "poisson_best_source_scale": pois_scale,
        "poisson_best_source_scale_rel": pois_rel_fit,
        "poisson_lin_best_source_scale": lin_scale,
        "poisson_lin_best_source_scale_rel": lin_rel_fit,
        "continuity_best_scale_offset": (cont_scale_off, cont_offset, cont_rel_off),
        "poisson_best_scale_offset": (pois_scale_off, pois_offset, pois_rel_off),
        "poisson_phi_sign_flip_best_scale": flip_scale,
        "poisson_phi_sign_flip_best_scale_rel": flip_rel_fit,
        "poisson_lin_phi_sign_flip_best_scale": lin_flip_scale,
        "poisson_lin_phi_sign_flip_best_scale_rel": lin_flip_rel_fit,
    }


def print_block(result: dict[str, Any]) -> None:
    print(f"\n[{result['label']}]")
    print("  continuity rel                    =", result["continuity_rel"])
    print("  exact Poisson rel                 =", result["poisson_exact_rel"])
    print("  linear Poisson rel                =", result["poisson_lin_rel"])
    print("  continuity best S scale, rel      =", (result["continuity_best_source_scale"], result["continuity_best_source_scale_rel"]))
    print("  exact Poisson best S scale, rel   =", (result["poisson_best_source_scale"], result["poisson_best_source_scale_rel"]))
    print("  linear Poisson best S scale, rel  =", (result["poisson_lin_best_source_scale"], result["poisson_lin_best_source_scale_rel"]))
    print("  continuity best S+offset          =", result["continuity_best_scale_offset"])
    print("  exact Poisson best S+offset       =", result["poisson_best_scale_offset"])
    print("  exact Poisson phi-sign flip       =", (result["poisson_phi_sign_flip_best_scale"], result["poisson_phi_sign_flip_best_scale_rel"]))
    print("  linear Poisson phi-sign flip      =", (result["poisson_lin_phi_sign_flip_best_scale"], result["poisson_lin_phi_sign_flip_best_scale_rel"]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", help="runtime-monitor snapshot .npz")
    parser.add_argument("--nonperiodic", action="store_true", help="use nonperiodic x/y/z derivatives; requires phi3")
    args = parser.parse_args()

    data = load_snapshot(pathlib.Path(args.snapshot))
    fields = project_fields(data)
    schema = detect_source_schema(data)
    source_rms = rms(fields.source)
    print("CFD SOURCE CONVENTION SWEEP")
    print("  schema       =", schema)
    print("  source rms   =", source_rms)
    print("  rho mean     =", float(np.mean(fields.rho)))
    print("  has phi3     =", fields.phi is not None)
    print("  periodic xyz =", not args.nonperiodic)

    current_result = summarize_interpretation(
        "stored Jx/Jy/Jz are brane currents",
        fields,
        stored_as_velocity=False,
        periodic_xyz=not args.nonperiodic,
    )
    velocity_result = summarize_interpretation(
        "stored Jx/Jy/Jz are velocities",
        fields,
        stored_as_velocity=True,
        periodic_xyz=not args.nonperiodic,
    )
    print_block(current_result)
    print_block(velocity_result)

    print("\nReading guide:")
    print("  If a best S scale is near -1, suspect a source sign convention.")
    print("  If it is near a constant not 1 and the rel drops sharply, suspect source normalization.")
    print("  If S+offset drops sharply, suspect the uniform compensation/background source term.")
    print("  If the velocity interpretation drops sharply, mx/my/mz were likely velocities, not currents.")
    print("  If none drops near the 0.05 threshold, this is a real solver/projection residual.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
