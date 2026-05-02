#!/usr/bin/env python3
"""Self-test the CFD runtime monitor postprocessor on synthetic exact fields."""
from __future__ import annotations

import hashlib
import json
import math

import numpy as np

from cfd_runtime_monitor_postprocess import (
    analyze_snapshot,
    divergence3,
    gradient3,
    laplacian3,
)


def assert_small(label: str, value: float, tol: float) -> None:
    if not np.isfinite(value) or abs(value) > tol:
        raise AssertionError(f"{label} failed: {value} > {tol}")


def assert_close(label: str, value: float, target: float, tol: float) -> None:
    if not np.isfinite(value) or abs(value - target) > tol:
        raise AssertionError(f"{label} failed: {value} vs {target} (tol {tol})")


def make_periodic_consistency_snapshot() -> dict[str, np.ndarray]:
    nx = ny = nz = 24
    nw = 5
    x = np.linspace(0.0, 1.0, nx, endpoint=False)
    y = np.linspace(0.0, 1.0, ny, endpoint=False)
    z = np.linspace(0.0, 1.0, nz, endpoint=False)
    w = np.linspace(0.0, 1.0, nw)
    W = np.ones_like(w)
    dWdw = np.zeros_like(w)
    rho0 = 2.5
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    phi = 0.07 * np.cos(2.0 * math.pi * X) * np.sin(2.0 * math.pi * Y) * np.cos(4.0 * math.pi * Z)
    spacings = (float(x[1] - x[0]), float(y[1] - y[0]), float(z[1] - z[0]))
    grad_phi = gradient3(phi, spacings, periodic=True)
    lap_phi = laplacian3(phi, spacings, periodic=True)
    Jx = rho0 * grad_phi[0]
    Jy = rho0 * grad_phi[1]
    Jz = rho0 * grad_phi[2]
    S = rho0 * lap_phi
    jw = (1.0 - w.reshape(1, 1, 1, nw)) * S[..., None]
    snapshot = {
        "rho": rho0 * np.ones((nx, ny, nz, nw), dtype=float),
        "jx": np.repeat(Jx[..., None], nw, axis=3),
        "jy": np.repeat(Jy[..., None], nw, axis=3),
        "jz": np.repeat(Jz[..., None], nw, axis=3),
        "jw": jw,
        "W": W,
        "dWdw": dWdw,
        "x": x,
        "y": y,
        "z": z,
        "w": w,
        "dt_rho": np.zeros((nx, ny, nz), dtype=float),
        "rho0": rho0,
        "phi3": phi,
        "N_probe": 1.0 - 2.0 * phi,
        "Phi_eff": phi,
    }
    return snapshot


def make_radial_snapshot(mu: float | None) -> dict[str, np.ndarray]:
    nx = ny = nz = 64
    nw = 4
    x = np.linspace(-1.6, 1.6, nx)
    y = np.linspace(-1.6, 1.6, ny)
    z = np.linspace(-1.6, 1.6, nz)
    w = np.linspace(0.0, 1.0, nw)
    W = np.ones_like(w)
    dWdw = np.zeros_like(w)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    eps = 0.08
    phi = -1.0 / np.sqrt(r * r + eps * eps)
    if mu is not None:
        phi = phi * np.exp(-mu * np.sqrt(r * r + eps * eps))
    rho0 = 1.0
    zeros4 = np.zeros((nx, ny, nz, nw), dtype=float)
    snapshot = {
        "rho": rho0 * np.ones((nx, ny, nz, nw), dtype=float),
        "jx": zeros4.copy(),
        "jy": zeros4.copy(),
        "jz": zeros4.copy(),
        "jw": zeros4.copy(),
        "W": W,
        "dWdw": dWdw,
        "x": x,
        "y": y,
        "z": z,
        "w": w,
        "dt_rho": np.zeros((nx, ny, nz), dtype=float),
        "rho0": rho0,
        "phi3": phi,
        "N_probe": 1.0 - 2.0 * phi,
        "Phi_eff": phi,
    }
    return snapshot


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_cfd_runtime_postprocessor",
        "pre_target_freeze": True,
        "target_blind": False,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "synthetic self-test for the CFD-side runtime monitor postprocessor",
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    periodic = analyze_snapshot(make_periodic_consistency_snapshot(), c_probe=1.0, bins=18, tail_fraction=0.35, periodic_xyz=True)
    assert_small("periodic max|R_cont|", periodic["max_abs_R_cont"], 1e-11)
    assert_small("periodic max|R_Pois_exact|", periodic["max_abs_R_pois_exact"], 1e-11)
    assert_small("periodic max|R_Pois_lin|", periodic["max_abs_R_pois_lin"], 1e-11)
    assert_close("periodic alpha_fit tail", periodic["alpha_fit_tail_mean"], 2.0, 5e-3)

    newton = analyze_snapshot(
        make_radial_snapshot(mu=None),
        c_probe=1.0,
        bins=26,
        tail_fraction=0.3,
        periodic_xyz=False,
        center=(0.0, 0.0, 0.0),
    )
    yukawa = analyze_snapshot(
        make_radial_snapshot(mu=1.4),
        c_probe=1.0,
        bins=26,
        tail_fraction=0.3,
        periodic_xyz=False,
        center=(0.0, 0.0, 0.0),
    )

    if not (newton["Q_r_tail_cv"] < yukawa["Q_r_tail_cv"]):
        raise AssertionError("Newton tail should plateau better than Yukawa")
    if not (abs(newton["mu_eff2_tail_median"]) < 0.2):
        raise AssertionError("Newton tail should stay near massless in the exterior")
    if not (yukawa["mu_eff2_tail_median"] > 1.0):
        raise AssertionError("Yukawa tail should show a positive exterior mass scale")
    assert_close("radial alpha_fit tail", newton["alpha_fit_tail_mean"], 2.0, 3e-2)

    print("STEP 34 CFD RUNTIME POSTPROCESSOR SELF-TEST")
    print("Validated the CFD-side runtime monitor postprocessor on synthetic exact fields.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Synthetic periodic consistency snapshot:")
    print("  max |R_cont| =", periodic["max_abs_R_cont"])
    print("  max |R_Pois_exact| =", periodic["max_abs_R_pois_exact"])
    print("  max |R_Pois_lin| =", periodic["max_abs_R_pois_lin"])
    print("  alpha_fit tail =", periodic["alpha_fit_tail_mean"])
    print("Exterior falsifier calibration:")
    print("  Newton Q_r tail cv =", newton["Q_r_tail_cv"])
    print("  Newton mu_eff^2 tail median =", newton["mu_eff2_tail_median"])
    print("  Yukawa Q_r tail cv =", yukawa["Q_r_tail_cv"])
    print("  Yukawa mu_eff^2 tail median =", yukawa["mu_eff2_tail_median"])
    print("  Newton alpha_fit tail =", newton["alpha_fit_tail_mean"])
    print("Snapshot schema for the real simulation:")
    print("  required keys = rho, jx, jy, jz, jw, W, x, y, z, w")
    print("  optional keys = dWdw, dt_rho, phi3, rho0, N_probe, Phi_eff")
    print("  monitor CLI = python cfd_runtime_monitor_postprocess.py snapshot.npz --output-json summary.json")
    print("Interpretation:")
    print("  The postprocessor now computes S_rho, brane continuity residuals, exact/linearized Poisson residuals, exterior Q_r plateaus, mu_eff^2 tails, and alpha_fit tails directly from a CFD snapshot.")
    print("  That is the shortest path from the hardcoded Branch-B patch to an actual simulation-side falsification run.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
