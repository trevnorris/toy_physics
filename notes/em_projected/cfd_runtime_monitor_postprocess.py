#!/usr/bin/env python3
"""Postprocess projected CFD snapshots into gravity-analog monitor diagnostics."""
from __future__ import annotations

import argparse
import json
import math
import pathlib
from typing import Any

import numpy as np


def require_keys(data: dict[str, np.ndarray], keys: list[str]) -> None:
    missing = [key for key in keys if key not in data]
    if missing:
        raise KeyError(f"missing required snapshot keys: {missing}")


def uniform_spacing(axis: np.ndarray) -> float:
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError("coordinate axes must be one-dimensional with at least two entries")
    diffs = np.diff(axis)
    if not np.allclose(diffs, diffs[0], rtol=1e-10, atol=1e-12):
        raise ValueError("runtime monitor currently expects uniform x/y/z/w axes")
    return float(diffs[0])


def first_derivative(field: np.ndarray, spacing: float, axis: int, periodic: bool) -> np.ndarray:
    if periodic:
        return (np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis)) / (2.0 * spacing)
    return np.gradient(field, spacing, axis=axis, edge_order=2)


def second_derivative(field: np.ndarray, spacing: float, axis: int, periodic: bool) -> np.ndarray:
    if periodic:
        return (np.roll(field, -1, axis=axis) - 2.0 * field + np.roll(field, 1, axis=axis)) / (spacing * spacing)
    return np.gradient(np.gradient(field, spacing, axis=axis, edge_order=2), spacing, axis=axis, edge_order=2)


def gradient3(field: np.ndarray, spacings: tuple[float, float, float], periodic: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dx, dy, dz = spacings
    return (
        first_derivative(field, dx, 0, periodic),
        first_derivative(field, dy, 1, periodic),
        first_derivative(field, dz, 2, periodic),
    )


def divergence3(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, spacings: tuple[float, float, float], periodic: bool) -> np.ndarray:
    dx, dy, dz = spacings
    return (
        first_derivative(vx, dx, 0, periodic)
        + first_derivative(vy, dy, 1, periodic)
        + first_derivative(vz, dz, 2, periodic)
    )


def laplacian3(field: np.ndarray, spacings: tuple[float, float, float], periodic: bool) -> np.ndarray:
    grad = gradient3(field, spacings, periodic)
    return divergence3(grad[0], grad[1], grad[2], spacings, periodic)


def integrate_w(field: np.ndarray, W: np.ndarray, w: np.ndarray) -> np.ndarray:
    return np.trapz(W.reshape((1,) * (field.ndim - 1) + (W.size,)) * field, x=w, axis=-1)


def compute_s_rho(jw: np.ndarray, W: np.ndarray, dWdw: np.ndarray, w: np.ndarray) -> np.ndarray:
    boundary_term = W[-1] * jw[..., -1] - W[0] * jw[..., 0]
    bulk_term = np.trapz(dWdw.reshape((1,) * (jw.ndim - 1) + (dWdw.size,)) * jw, x=w, axis=-1)
    return -boundary_term + bulk_term


def solve_phi_from_velocity(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, spacings: tuple[float, float, float]) -> np.ndarray:
    nx, ny, nz = vx.shape
    dx, dy, dz = spacings
    div_v = divergence3(vx, vy, vz, spacings, periodic=True)
    kx = 2.0 * math.pi * np.fft.fftfreq(nx, d=dx)
    ky = 2.0 * math.pi * np.fft.fftfreq(ny, d=dy)
    kz = 2.0 * math.pi * np.fft.fftfreq(nz, d=dz)
    k2 = (
        kx[:, None, None] ** 2
        + ky[None, :, None] ** 2
        + kz[None, None, :] ** 2
    )
    div_hat = np.fft.fftn(div_v)
    phi_hat = np.zeros_like(div_hat, dtype=np.complex128)
    nonzero = k2 > 0.0
    phi_hat[nonzero] = -div_hat[nonzero] / k2[nonzero]
    return np.fft.ifftn(phi_hat).real


def auto_center(source: np.ndarray, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float, float]:
    weights = np.abs(source)
    total = float(np.sum(weights))
    if total == 0.0:
        return (float(x.mean()), float(y.mean()), float(z.mean()))
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    return (
        float(np.sum(weights * X) / total),
        float(np.sum(weights * Y) / total),
        float(np.sum(weights * Z) / total),
    )


def shell_average(field: np.ndarray, r: np.ndarray, bins: int, r_min: float = 0.0, r_max: float | None = None) -> dict[str, list[float]]:
    if r_max is None:
        r_max = float(np.max(r))
    edges = np.linspace(r_min, r_max, bins + 1)
    shell_means: list[float] = []
    shell_radii: list[float] = []
    shell_counts: list[float] = []
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (r >= left) & (r < right)
        values = field[mask]
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            shell_means.append(float("nan"))
            shell_radii.append(0.5 * (left + right))
            shell_counts.append(0.0)
            continue
        shell_means.append(float(np.mean(finite)))
        shell_radii.append(0.5 * (left + right))
        shell_counts.append(float(finite.size))
    return {"r": shell_radii, "mean": shell_means, "count": shell_counts}


def tail_slice(values: list[float], counts: list[float], tail_fraction: float) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    count_array = np.asarray(counts, dtype=float)
    valid = np.isfinite(array) & (count_array > 0.0)
    valid_idx = np.flatnonzero(valid)
    if valid_idx.size == 0:
        return np.asarray([], dtype=float)
    tail_len = max(3, int(math.ceil(tail_fraction * valid_idx.size)))
    chosen = valid_idx[-tail_len:]
    return array[chosen]


def max_abs(array: np.ndarray) -> float:
    return float(np.max(np.abs(array)))


def rms(array: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(array))))


def detect_source_schema(data: dict[str, np.ndarray]) -> str:
    full4d = {"rho", "jx", "jy", "jz", "jw", "W", "x", "y", "z", "w"}
    projected3d = {"rho_brane", "Jx_brane", "Jy_brane", "Jz_brane", "S_rho", "x", "y", "z"}
    if full4d.issubset(data):
        return "full_4d"
    if projected3d.issubset(data):
        return "projected_3d"
    raise KeyError(
        "snapshot does not match a supported schema; expected either full 4D keys "
        f"{sorted(full4d)} or projected 3D keys {sorted(projected3d)}"
    )


def analyze_snapshot(
    data: dict[str, np.ndarray],
    *,
    c_probe: float = 1.0,
    bins: int = 24,
    tail_fraction: float = 0.35,
    periodic_xyz: bool = True,
    center: tuple[float, float, float] | None = None,
) -> dict[str, Any]:
    source_schema = detect_source_schema(data)
    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    z = np.asarray(data["z"], dtype=float)
    dx = uniform_spacing(x)
    dy = uniform_spacing(y)
    dz = uniform_spacing(z)
    spacings = (dx, dy, dz)

    if source_schema == "full_4d":
        rho = np.asarray(data["rho"], dtype=float)
        jx = np.asarray(data["jx"], dtype=float)
        jy = np.asarray(data["jy"], dtype=float)
        jz = np.asarray(data["jz"], dtype=float)
        jw = np.asarray(data["jw"], dtype=float)
        W = np.asarray(data["W"], dtype=float)
        w = np.asarray(data["w"], dtype=float)
        dWdw = np.asarray(data.get("dWdw", np.gradient(W, w, edge_order=2)), dtype=float)

        if rho.shape != jx.shape or rho.shape != jy.shape or rho.shape != jz.shape or rho.shape != jw.shape:
            raise ValueError("rho and all current components must have the same 4D shape")
        if rho.ndim != 4:
            raise ValueError("rho and current components must have shape (nx, ny, nz, nw)")
        if rho.shape[-1] != W.size or rho.shape[-1] != w.size:
            raise ValueError("w-axis length must match rho[..., nw] and weight arrays")

        _ = uniform_spacing(w)
        rho_brane = integrate_w(rho, W, w)
        Jx_brane = integrate_w(jx, W, w)
        Jy_brane = integrate_w(jy, W, w)
        Jz_brane = integrate_w(jz, W, w)
        S_rho = compute_s_rho(jw, W, dWdw, w)
    else:
        rho_brane = np.asarray(data["rho_brane"], dtype=float)
        Jx_brane = np.asarray(data["Jx_brane"], dtype=float)
        Jy_brane = np.asarray(data["Jy_brane"], dtype=float)
        Jz_brane = np.asarray(data["Jz_brane"], dtype=float)
        S_rho = np.asarray(data["S_rho"], dtype=float)
        if rho_brane.shape != Jx_brane.shape or rho_brane.shape != Jy_brane.shape or rho_brane.shape != Jz_brane.shape or rho_brane.shape != S_rho.shape:
            raise ValueError("projected 3D rho/J/S arrays must all have the same shape")
        if rho_brane.ndim != 3:
            raise ValueError("projected 3D arrays must have shape (nx, ny, nz)")

    rho_floor = np.maximum(np.abs(rho_brane), 1e-12)
    vx = Jx_brane / rho_floor
    vy = Jy_brane / rho_floor
    vz = Jz_brane / rho_floor

    divJ = divergence3(Jx_brane, Jy_brane, Jz_brane, spacings, periodic_xyz)
    dt_rho = np.asarray(data.get("dt_rho", np.zeros_like(rho_brane)), dtype=float)
    if dt_rho.shape != rho_brane.shape:
        raise ValueError("dt_rho must have shape (nx, ny, nz)")
    R_cont = dt_rho + divJ - S_rho

    phi3_input = data.get("phi3")
    if phi3_input is None:
        if not periodic_xyz:
            raise ValueError("phi3 must be provided when periodic_xyz is false")
        phi3 = solve_phi_from_velocity(vx, vy, vz, spacings)
    else:
        phi3 = np.asarray(phi3_input, dtype=float)
        if phi3.shape != rho_brane.shape:
            raise ValueError("phi3 must have shape (nx, ny, nz)")

    lap_phi = laplacian3(phi3, spacings, periodic_xyz)
    grad_rho = gradient3(rho_brane, spacings, periodic_xyz)
    grad_phi = gradient3(phi3, spacings, periodic_xyz)
    grad_rho_dot_v = grad_rho[0] * vx + grad_rho[1] * vy + grad_rho[2] * vz
    R_pois_exact = rho_brane * lap_phi - S_rho + dt_rho + grad_rho_dot_v
    rho0 = float(data.get("rho0", float(np.mean(rho_brane))))
    R_pois_lin = rho0 * lap_phi - S_rho

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    if center is None:
        center = auto_center(S_rho, x, y, z)
    cx, cy, cz = center
    rx = X - cx
    ry = Y - cy
    rz = Z - cz
    r = np.sqrt(rx * rx + ry * ry + rz * rz)
    r_safe = np.where(r > 0.0, r, np.nan)
    dphi_dr = (grad_phi[0] * rx + grad_phi[1] * ry + grad_phi[2] * rz) / r_safe
    Q_r = 4.0 * math.pi * r * r * dphi_dr
    mu_eff2 = np.full_like(phi3, np.nan, dtype=float)
    denom_mask = np.abs(phi3) > 1e-12
    mu_eff2[denom_mask] = (lap_phi[denom_mask] - S_rho[denom_mask] / rho0) / phi3[denom_mask]

    shell_Q = shell_average(Q_r, r, bins=bins, r_min=0.0)
    shell_mu = shell_average(mu_eff2, r, bins=bins, r_min=0.0)
    tail_Q = tail_slice(shell_Q["mean"], shell_Q["count"], tail_fraction)
    tail_mu = tail_slice(shell_mu["mean"], shell_mu["count"], tail_fraction)
    q_plateau_cv = float(np.std(tail_Q) / max(abs(float(np.mean(tail_Q))), 1e-12)) if tail_Q.size else float("nan")
    mu_tail_median = float(np.median(tail_mu)) if tail_mu.size else float("nan")

    summary: dict[str, Any] = {
        "schema": {
            "source_schema": source_schema,
            "supported_required_keys": {
                "full_4d": ["rho", "jx", "jy", "jz", "jw", "W", "x", "y", "z", "w"],
                "projected_3d": ["rho_brane", "Jx_brane", "Jy_brane", "Jz_brane", "S_rho", "x", "y", "z"],
            },
            "optional_keys": ["dWdw", "dt_rho", "phi3", "rho0", "N_probe", "Phi_eff"],
        },
        "center": {"x": cx, "y": cy, "z": cz},
        "rho0": rho0,
        "max_abs_S_rho": max_abs(S_rho),
        "rms_S_rho": rms(S_rho),
        "max_abs_R_cont": max_abs(R_cont),
        "rms_R_cont": rms(R_cont),
        "max_abs_R_pois_exact": max_abs(R_pois_exact),
        "rms_R_pois_exact": rms(R_pois_exact),
        "max_abs_R_pois_lin": max_abs(R_pois_lin),
        "rms_R_pois_lin": rms(R_pois_lin),
        "Q_r_tail_mean": float(np.mean(tail_Q)) if tail_Q.size else float("nan"),
        "Q_r_tail_cv": q_plateau_cv,
        "mu_eff2_tail_median": mu_tail_median,
        "shell_Q_r": shell_Q,
        "shell_mu_eff2": shell_mu,
    }

    if "N_probe" in data and "Phi_eff" in data:
        N_probe = np.asarray(data["N_probe"], dtype=float)
        Phi_eff = np.asarray(data["Phi_eff"], dtype=float)
        if N_probe.shape != rho_brane.shape or Phi_eff.shape != rho_brane.shape:
            raise ValueError("N_probe and Phi_eff must have shape (nx, ny, nz)")
        alpha_fit = np.full_like(Phi_eff, np.nan, dtype=float)
        mask = np.abs(Phi_eff) > 1e-12
        alpha_fit[mask] = -c_probe * c_probe * (N_probe[mask] - 1.0) / Phi_eff[mask]
        shell_alpha = shell_average(alpha_fit, r, bins=bins, r_min=0.0)
        tail_alpha = tail_slice(shell_alpha["mean"], shell_alpha["count"], tail_fraction)
        summary.update(
            {
                "alpha_fit_tail_mean": float(np.mean(tail_alpha)) if tail_alpha.size else float("nan"),
                "alpha_fit_tail_std": float(np.std(tail_alpha)) if tail_alpha.size else float("nan"),
                "shell_alpha_fit": shell_alpha,
            }
        )

    return summary


def load_snapshot(path: pathlib.Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as raw:
        return {key: raw[key] for key in raw.files}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="snapshot .npz file")
    parser.add_argument("--output-json", help="write summary JSON to this path")
    parser.add_argument("--bins", type=int, default=24)
    parser.add_argument("--tail-fraction", type=float, default=0.35)
    parser.add_argument("--c-probe", type=float, default=1.0)
    parser.add_argument("--nonperiodic", action="store_true", help="use nonperiodic x/y/z derivatives; requires phi3 in the snapshot")
    parser.add_argument("--center", nargs=3, type=float, metavar=("X0", "Y0", "Z0"))
    args = parser.parse_args()

    snapshot = load_snapshot(pathlib.Path(args.input))
    summary = analyze_snapshot(
        snapshot,
        c_probe=args.c_probe,
        bins=args.bins,
        tail_fraction=args.tail_fraction,
        periodic_xyz=not args.nonperiodic,
        center=tuple(args.center) if args.center else None,
    )

    print("CFD RUNTIME MONITOR SUMMARY")
    print("  max |R_cont|      =", summary["max_abs_R_cont"])
    print("  max |R_Pois_exact|=", summary["max_abs_R_pois_exact"])
    print("  max |R_Pois_lin|  =", summary["max_abs_R_pois_lin"])
    print("  Q_r tail mean     =", summary["Q_r_tail_mean"])
    print("  Q_r tail cv       =", summary["Q_r_tail_cv"])
    print("  mu_eff^2 tail med =", summary["mu_eff2_tail_median"])
    if "alpha_fit_tail_mean" in summary:
        print("  alpha_fit tail    =", summary["alpha_fit_tail_mean"])
        print("  alpha_fit tail std=", summary["alpha_fit_tail_std"])

    if args.output_json:
        pathlib.Path(args.output_json).write_text(json.dumps(summary, indent=2, sort_keys=True))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
